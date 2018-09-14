/*
 *  autoclean_maps.cc  -  automatically clean problematic regions from maps
 *
 *  Copyright (c) 2008-2018 Pittsburgh Supercomputing Center,
 *                          Carnegie Mellon University
 *
 *  This file is part of AlignTK.
 *
 *  AlignTK is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AlignTK is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AlignTK.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Acknowledgements:
 *     Development of this code was supported in part by
 *       NIH NCRR grant 5P41RR006009 and
 *       NIH NIGMS grant P41GM103712
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>
#include <stdarg.h>
#include <sys/time.h>
#include <string>
#include <vector>
#include <algorithm>
#include "imio.h"
#include "invert.h"

using std::string;
using std::max;
using std::min;
using std::vector;

#define MAX_LINE_LENGTH		255

struct Pair
{
  char *imageName;
  int imageMinX, imageMaxX, imageMinY, imageMaxY; 
  char *refName;
  int refMinX, refMaxX, refMinY, refMaxY;
  char *pairName;
  bool modified;
};

struct MapElementAndDistortion
{
  int x, y;
  float distortion;
};

bool lesserDistortion (const MapElementAndDistortion &m0,
		       const MapElementAndDistortion &m1)
{ return(m0.distortion < m1.distortion); }

#define VALID_BIT	0x01
#define REJECT_BIT	0x02	// red
#define ACCEPT_BIT	0x04	// blue
#define THRESHOLD_BIT	0x08	// orange or green
#define CLUSTER_BIT	0x10	// green

unsigned char *image = 0;
int imageWidth, imageHeight;
int imageMinX, imageMaxX;
int imageMinY, imageMaxY;

int mapWidth, mapHeight;
int mapOffsetX, mapOffsetY;
int mapLevel;
int mapFactor;
MapElement *map = 0;
InverseMap *inverseMap = 0;

unsigned char *mapMask = 0;
int *cluster = 0;
vector<int> seeds;       // cluster seeds

float avgOrthogonal;
float avgDiagonal;

char extension[8];
char pairsFile[PATH_MAX];
char imagesFile[PATH_MAX];
char inputName[PATH_MAX];
char mapsName[PATH_MAX];
char corrName[PATH_MAX];
char outputName[PATH_MAX];
char logName[PATH_MAX];
int reductionFactor = 1;
Pair* pairs;
int nPairs;
bool readOnly = false;
float oThreshold = 0.1;
float dThreshold = 0.1;
FILE *logFile = NULL;

void ProcessPair (int index);
void SelectSeed ();
void UpdateMask ();
void Log (const char *fmt, ...);
void Error (const char *fmt, ...);
char *GetTimestamp (char *timestamp, size_t size);

int
main (int argc, char **argv)
{
  FILE *f;
  int w, h;
  int m;
  int v;
  int i;
  int n;
  char tc;
  int error;
  int pos;
  char fn[PATH_MAX];
  int len;
  int z;
  bool maximizeWindow = false;
  int windowFactor = 1;
  int windowWidth, windowHeight;
  char imgn[PATH_MAX], refn[PATH_MAX], pairn[PATH_MAX];
  int imgMinX, imgMaxX, imgMinY, imgMaxY;
  int refMinX, refMaxX, refMinY, refMaxY;

  error = 0;
  strcpy(extension, "tif");
  pairsFile[0] = '\0';
  imagesFile[0] = '\0';
  inputName[0] = '\0';
  mapsName[0] = '\0';
  corrName[0] = '\0';
  outputName[0] = '\0';
  logName[0] = '\0';
  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-big") == 0)
      maximizeWindow = true;
    else if (strcmp(argv[i], "-input") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-input error\n");
	    break;
	  }
	strcpy(inputName, argv[i]);
      }
    else if (strcmp(argv[i], "-maps") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-maps error\n");
	    break;
	  }
	strcpy(mapsName, argv[i]);
      }
    else if (strcmp(argv[i], "-corr") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-corr error\n");
	    break;
	  }
	strcpy(corrName, argv[i]);
      }
    else if (strcmp(argv[i], "-output") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-output error\n");
	    break;
	  }
	strcpy(outputName, argv[i]);
      }
    else if (strcmp(argv[i], "-log") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-log error\n");
	    break;
	  }
	strcpy(logName, argv[i]);
      }	
    else if (strcmp(argv[i], "-pairs") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-pairs error\n");
	    break;
	  }
	strcpy(pairsFile, argv[i]);
      }
    else if (strcmp(argv[i], "-images") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-images error\n");
	    break;
	  }
	strcpy(imagesFile, argv[i]);
      }
    else if (strcmp(argv[i], "-reduction") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%d", &reductionFactor) != 1)
	  {
	    error = 1;
	    fprintf(stderr, "-reduction error\n");
	    break;
	  }
      }
    else if (strcmp(argv[i], "-orthogonal_threshold") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%f", &oThreshold) != 1)
	  {
	    error = 1;
	    fprintf(stderr, "-orthogonal_threshold error\n");
	    break;
	  }
      }
    else if (strcmp(argv[i], "-diagonal_threshold") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%f", &dThreshold) != 1)
	  {
	    error = 1;
	    fprintf(stderr, "-diagonal_threshold error\n");
	    break;
	  }
      }
    else if (strcmp(argv[i], "-readonly") == 0)
      readOnly = true;
    else if (strcmp(argv[i], "-pgm") == 0)
      strcpy(extension, "pgm");
    else
      error = 1;

  if (argc == 1 || error)
    {
      fprintf(stderr, "Usage: autoclean_maps -input images_file_prefix\n");
      fprintf(stderr, "                  -maps maps_file_prefix\n");
      fprintf(stderr, "                  -output maps_file_prefix\n");
      fprintf(stderr, "                  [-pairs pairfile]\n");
      fprintf(stderr, "                  [-images imagesfile]\n");
      fprintf(stderr, "                  [-reduction image_reduction_factor]\n");
      fprintf(stderr, "                  [-orthogonal_threshold value]\n");
      fprintf(stderr, "                  [-diagonal_threshold value]\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (inputName[0] == '\0' || mapsName[0] == '\0' ||
      pairsFile[0] == '\0' && imagesFile[0] == '\0')
    {
      fprintf(stderr, "-input, -maps, and -pairs (or -images) parameters must be specified.\n");
      exit(1);
    }
  if (!readOnly && (corrName[0] == '\0' || outputName[0] == '\0'))
    {
      fprintf(stderr, "-corr and -output parameters must be specified if not readonly.\n");
      exit(1);
    }

  pairs = 0;
  if (pairsFile[0] != '\0')
    {
      char line[MAX_LINE_LENGTH+1];

      f = fopen(pairsFile, "r");
      if (f == NULL)
	{
	  fprintf(stderr, "Could not open pairs file %s\n", pairsFile);
	  exit(1);
	}
      
      
      while (fgets(line, MAX_LINE_LENGTH, f) != NULL)
	{
	  if (line[0] == '\0' || line[0] == '#')
	    continue;
	  if (sscanf(line, "%s %d %d %d %d %s %d %d %d %d %s",
		     imgn, &imgMinX, &imgMaxX, &imgMinY, &imgMaxY,
		     refn, &refMinX, &refMaxX, &refMinY, &refMaxY,
		     pairn) != 11)
	    {
	      if (sscanf(line, "%s %s %s", imgn, refn, pairn) != 3)
		{
		  fprintf(stderr, "Invalid line in pairs file %s:\n%s\n", pairsFile, line);
		  exit(1);
		}
	      imgMinX = -1;
	      imgMaxX = -1;
	      imgMinY = -1;
	      imgMaxY = -1;
	      refMinX = -1;
	      refMaxX = -1;
	      refMinY = -1;
	      refMaxY = -1;
	    }

	  if ((nPairs & 1023) == 0)
	    pairs = (Pair *) realloc(pairs, (nPairs + 1024) * sizeof(Pair));
	  pairs[nPairs].imageName = (char *) malloc(strlen(imgn)+1);
	  strcpy(pairs[nPairs].imageName, imgn);
	  pairs[nPairs].imageMinX = imgMinX;
	  pairs[nPairs].imageMaxX = imgMaxX;
	  pairs[nPairs].imageMinY = imgMinY;
	  pairs[nPairs].imageMaxY = imgMaxY;
	  pairs[nPairs].refName = (char *) malloc(strlen(refn)+1);
	  strcpy(pairs[nPairs].refName, refn);
	  pairs[nPairs].refMinX = refMinX;
	  pairs[nPairs].refMaxX = refMaxX;
	  pairs[nPairs].refMinY = refMinY;
	  pairs[nPairs].refMaxY = refMaxY;
	  pairs[nPairs].pairName = (char *) malloc(strlen(pairn)+1);
	  strcpy(pairs[nPairs].pairName, pairn);
	  ++nPairs;
	}
      fclose(f);
      printf("%d pairs listed in pairs file.\n", nPairs);
    }
  else
    {
      f = fopen(imagesFile, "r");
      if (f == NULL)
	{
	  fprintf(stderr, "Could not open images file %s\n", imagesFile);
	  exit(1);
	}
      char line[MAX_LINE_LENGTH+1];
      while (fgets(line, MAX_LINE_LENGTH, f) != NULL)
	{
	  n = sscanf(line, "%s%d%d%d%d",
		     imgn, &imgMinX, &imgMaxX, &imgMinY, &imgMaxY);
	  if (n != 1 && n != 5)
	    {
	      fprintf(stderr, "Invalid line in images file %s\n", imagesFile);
	      exit(1);
	    }
	  if (n == 1)
	    {
	      imgMinX = -1;
	      imgMaxX = -1;
	      imgMinY = -1;
	      imgMaxY = -1;
	    }
	  if ((nPairs & 1023) == 0)
	    pairs = (Pair *) realloc(pairs, (nPairs + 1024) * sizeof(Pair));
	  pairs[nPairs].imageName = (char *) malloc(strlen(imgn)+1);
	  strcpy(pairs[nPairs].imageName, imgn);
	  pairs[nPairs].imageMinX = imgMinX;
	  pairs[nPairs].imageMaxX = imgMaxX;
	  pairs[nPairs].imageMinY = imgMinY;
	  pairs[nPairs].imageMaxY = imgMaxY;
	  pairs[nPairs].refName = (char *) malloc(1);
	  pairs[nPairs].refName[0] = '\0';
	  pairs[nPairs].refMinX = -1;
	  pairs[nPairs].refMaxX = -1;
	  pairs[nPairs].refMinY = -1;
	  pairs[nPairs].refMaxY = -1;
	  pairs[nPairs].pairName = (char *) malloc(1);
	  pairs[nPairs].pairName[0] = '\0';
	  ++nPairs;
	}
      fclose(f);
      printf("%d pairs listed in pairs file.\n", nPairs);
    }

  // go through all pairs one-by-one
  for (i = 0; i < nPairs; ++i)
    ProcessPair(i);

  return(0);
}

void
ProcessPair (int index)
{
  int n;
  int i;
  int x, y;
  float *warp;
  FILE *f;
  int refSection;
  char fn[PATH_MAX];
  char errorMsg[PATH_MAX + 256];
  char imgn[PATH_MAX], refn[PATH_MAX];
  double sumOrthogonal, sumDiagonal;
  double distX, distY;
  double distA, distB;

  if (image != 0)
    {
      free(image);
      image = 0;
    }

  /* read in the image */
  if (pairs[index].imageMinX >= 0)
    imageMinX = pairs[index].imageMinX / reductionFactor;
  else
    imageMinX = -1;
  if (pairs[index].imageMaxX >= 0)
    imageMaxX = pairs[index].imageMaxX / reductionFactor;
  else
    imageMaxX = -1;
  if (pairs[index].imageMinY >= 0)
    imageMinY = pairs[index].imageMinY / reductionFactor;
  else
    imageMinY = -1;
  if (pairs[index].imageMaxY >= 0)
    imageMaxY = pairs[index].imageMaxY / reductionFactor;
  else
    imageMaxY = -1;

  /* read in the map */
  sprintf(fn, "%s%s.map", mapsName,
	  (pairs[index].pairName[0] != '\0' ?
	   pairs[index].pairName : pairs[index].imageName));
  Log("Reading map %s\n", fn);
  if (!ReadMap(fn, &map, &mapLevel, &mapWidth, &mapHeight,
	       &mapOffsetX, &mapOffsetY, imgn, refn, errorMsg))
    {
      fprintf(stderr, "Error reading map %s:\n  %s\n", fn, errorMsg);
      exit(1);
    }

  /* construct the inverse map */
  mapFactor = (1 << mapLevel);
  inverseMap = InvertMap(map, mapWidth, mapHeight);

  /* compute the average ratios */
  sumOrthogonal = 0.0;
  sumDiagonal = 0.0;
  n = 0;
  for (y = 0; y < mapHeight-1; ++y)
    for (x = 0; x < mapWidth-1; ++x)
      {
	if (map[y * mapWidth + x].c == 0.0 || map[y * mapWidth + x + 1].c == 0.0 ||
	    map[(y + 1) * mapWidth + x].c == 0.0 || map[(y + 1) * mapWidth + x + 1].c == 0.0)
	  continue;
	++n;
	distX = 0.5 * hypot(map[y * mapWidth + x].x - map[y * mapWidth + x + 1].x,
		      map[y * mapWidth + x].y - map[y * mapWidth + x + 1].y) +
	  0.5 * hypot(map[(y + 1) * mapWidth + x].x - map[(y + 1) * mapWidth + x + 1].x,
		      map[(y + 1) * mapWidth + x].y - map[(y + 1) * mapWidth + x + 1].y);
	distY = 0.5 * hypot(map[y * mapWidth + x].x - map[(y + 1) * mapWidth + x].x,
			    map[y * mapWidth + x].y - map[(y + 1) * mapWidth + x].y) +
	  0.5 * hypot(map[y * mapWidth + x + 1].x - map[(y + 1) * mapWidth + x + 1].x,
		      map[y * mapWidth + x + 1].y - map[(y + 1) * mapWidth + x + 1].y);
	sumOrthogonal += distY / distX;
	distA = hypot(map[y * mapWidth + x].x - map[(y + 1) * mapWidth + x + 1].x,
		      map[y * mapWidth + x].y - map[(y + 1) * mapWidth + x + 1].y);
	distB = hypot(map[y * mapWidth + x + 1].x - map[(y + 1) * mapWidth + x].x,
		      map[y * mapWidth + x + 1].y - map[(y + 1) * mapWidth + x].y);
	sumDiagonal += distB/distA;
      }
  avgOrthogonal = sumOrthogonal / n;
  avgDiagonal = sumDiagonal / n;
  printf("avg orthogonal ratio = %f\n", avgOrthogonal);
  printf("avg diagonal ratio = %f\n", avgDiagonal);

  /* set up the map mask */
  printf("mapMask = %p\n", image);
  if (mapMask != 0)
    free(mapMask);
  mapMask = (unsigned char *) malloc(mapHeight * mapWidth * sizeof(unsigned char));
  memset(mapMask, 0, mapHeight * mapWidth * sizeof(unsigned char));
  printf("cluster = %p\n", image);
  if (cluster != 0)
    free(cluster);
  cluster = (int *) malloc(mapHeight * mapWidth * sizeof(int));
  for (int y = 0; y < mapHeight; ++y)
    for (int x = 0; x < mapWidth; ++x)
      cluster[y * mapWidth + x] = -1;

  seeds.clear();

  sprintf(fn, "%s%s.corr", corrName,
	  (pairs[index].pairName[0] != '\0' ?
	   pairs[index].pairName : pairs[index].imageName));
  f = fopen(fn, "r");
  if (f != NULL)
    {
      double oThreshold, dThreshold;
      if (fscanf(f, "%lf %lf", &oThreshold, &dThreshold) != 2)
	{
	  fprintf(stderr, "Invalid first line in %s\n", fn);
	  exit(1);
	}

      while (fscanf(f, "%d %d", &x, &y) == 2 && x >= 0 && y >= 0)
	seeds.push_back(y*mapWidth+x);
      while (fscanf(f, "%d %d", &x, &y) == 2 && x >= 0 && y >= 0)
	mapMask[y*mapWidth+x] |= REJECT_BIT;
      while (fscanf(f, "%d %d", &x, &y) == 2 && x >= 0 && y >= 0)
	mapMask[y*mapWidth+x] |= ACCEPT_BIT;
      fclose(f);
    }
  else
    {
      SelectSeed();
    }
  UpdateMask();

  Log("Set mapFactor to %d  (%d %d)\n", mapFactor, imageWidth, mapWidth);
  Log("image dim = (%d %d) imageMin = (%d %d)  mapOffset = (%d %d)\n",
      imageWidth, imageHeight, imageMinX, imageMinY,
      mapOffsetX, mapOffsetY);

  sprintf(fn, "%s%s.corr", corrName,
	  (pairs[index].pairName[0] != '\0' ?
	   pairs[index].pairName : pairs[index].imageName));
  if ((f = fopen(fn, "w")) == NULL)
    Error("Could not open file %s for writing.\n", fn);
  fprintf(f, "%f %f\n",
	  oThreshold, dThreshold);
  for (int i = 0; i < seeds.size(); ++i)
    fprintf(f, "%d %d\n", seeds[i] % mapWidth, seeds[i] / mapWidth);
  fprintf(f, "-1 -1\n");
  for (int y = 0; y < mapHeight-1; ++y)
    for (int x = 0; x < mapWidth-1; ++x)
      if (mapMask[y*mapWidth+x] & REJECT_BIT)
	fprintf(f, "%d %d\n", x, y);
  fprintf(f, "-1 -1\n");
  for (int y = 0; y < mapHeight-1; ++y)
    for (int x = 0; x < mapWidth-1; ++x)
      if (mapMask[y*mapWidth+x] & ACCEPT_BIT)
	fprintf(f, "%d %d\n", x, y);
  fprintf(f, "-1 -1\n");
  fclose(f);

  MapElement *outMap = (MapElement*) malloc(mapHeight * mapWidth *
					    sizeof(MapElement));
  for (int y = 0; y < mapHeight; ++y)
    for (int x = 0; x < mapWidth; ++x)
      {
	outMap[y*mapWidth + x].x = map[y*mapWidth + x].x;
	outMap[y*mapWidth + x].y = map[y*mapWidth + x].y;
	if (x > 0 && y > 0 && mapMask[(y-1)*mapWidth + x-1] & CLUSTER_BIT ||
	    x > 0 && mapMask[y*mapWidth + x-1] & CLUSTER_BIT ||
	    y > 0 && mapMask[(y-1)*mapWidth + x] & CLUSTER_BIT ||
	    mapMask[y*mapWidth + x] & CLUSTER_BIT)
	  outMap[y*mapWidth + x].c = map[y*mapWidth + x].c;
	else
	  outMap[y*mapWidth + x].c = 0.0;
      }
  sprintf(fn, "%s%s.map", outputName,
	  (pairs[index].pairName[0] != '\0' ?
	   pairs[index].pairName : pairs[index].imageName));
  if (!WriteMap(fn, outMap, mapLevel, mapWidth, mapHeight,
		mapOffsetX, mapOffsetY,
		pairs[index].imageName,
		(pairs[index].refName[0] != '\0' ?
		 pairs[index].refName : pairs[index].imageName),
		UncompressedMap, errorMsg))
    Error("Could not write map %s\n", fn);
  free(outMap);
}

void
SelectSeed ()
{
  // find the top 10% of map elements in terms of low distortion score;
  vector<MapElementAndDistortion> m;
  MapElementAndDistortion me;
  int x, y;
  double distX, distY;
  double distA, distB;
  float diagonal, orthogonal;
  for (y = 0; y < mapHeight-1; ++y)
    for (x = 0; x < mapWidth-1; ++x)
      {
	if (map[y * mapWidth + x].c == 0.0 || map[y * mapWidth + x + 1].c == 0.0 ||
	    map[(y + 1) * mapWidth + x].c == 0.0 || map[(y + 1) * mapWidth + x + 1].c == 0.0)
	  continue;
	distX = 0.5 * hypot(map[y * mapWidth + x].x - map[y * mapWidth + x + 1].x,
			    map[y * mapWidth + x].y - map[y * mapWidth + x + 1].y) +
	  0.5 * hypot(map[(y + 1) * mapWidth + x].x - map[(y + 1) * mapWidth + x + 1].x,
		      map[(y + 1) * mapWidth + x].y - map[(y + 1) * mapWidth + x + 1].y);
	distY = 0.5 * hypot(map[y * mapWidth + x].x - map[(y + 1) * mapWidth + x].x,
			    map[y * mapWidth + x].y - map[(y + 1) * mapWidth + x].y) +
	  0.5 * hypot(map[y * mapWidth + x + 1].x - map[(y + 1) * mapWidth + x + 1].x,
		      map[y * mapWidth + x + 1].y - map[(y + 1) * mapWidth + x + 1].y);
	orthogonal = distY / distX;
	distA = hypot(map[y * mapWidth + x].x - map[(y + 1) * mapWidth + x + 1].x,
		      map[y * mapWidth + x].y - map[(y + 1) * mapWidth + x + 1].y);
	distB = hypot(map[y * mapWidth + x + 1].x - map[(y + 1) * mapWidth + x].x,
		      map[y * mapWidth + x + 1].y - map[(y + 1) * mapWidth + x].y);
	diagonal = distB / distA;
	me.x = x;
	me.y = y;
	me.distortion = fabs(orthogonal - avgOrthogonal) / avgOrthogonal +
	  fabs(diagonal - avgDiagonal) / avgDiagonal;
	m.push_back(me);
      }
  sort(m.begin(), m.end(), lesserDistortion);

  //   find the average location of those
  double sumX, sumY;
  double avgX, avgY;
  int n = (m.size() + 9) / 10;
  if (n == 0)
    return;
  for (int i = 0; i < n; ++i)
    {
      sumX += m[i].x;
      sumY += m[i].y;
    }
  avgX = sumX / n;
  avgY = sumY / n;

  //   pick the element that is located closest to the average
  double closestDist = 1.0e30;
  int closest = -1;
  double dist;
  for (int i = 0; i < n; ++i)
    {
      dist = hypot(m[i].x - avgX, m[i].y - avgY);
      if (dist < closestDist)
	{
	  closest = i;
	  closestDist = dist;
	}
    }
  if (closest >= 0)
    seeds.push_back(m[closest].y * mapWidth + m[closest].x);
}

void
UpdateMask ()
{
  int x, y;
  double distX, distY;
  double distA, distB;
  bool overThreshold;
  float diagonal, orthogonal;

  for (y = 0; y < mapHeight-1; ++y)
    for (x = 0; x < mapWidth-1; ++x)
      {
	if (map[y*mapWidth+x].c == 0.0 || map[y*mapWidth+x+1].c == 0.0 ||
	    map[(y+1)*mapWidth+x].c == 0.0 || map[(y+1)*mapWidth+x+1].c == 0.0)
	  {
	    mapMask[y * mapWidth + x] = 0;
	    continue;
	  }
	mapMask[y * mapWidth + x] |= VALID_BIT;

	if (mapMask[y*mapWidth + x] & REJECT_BIT)
	  {
	    mapMask[y*mapWidth + x] &= ~THRESHOLD_BIT;
	    continue;
	  }
	if (mapMask[y*mapWidth + x] & ACCEPT_BIT)
	  {
	    mapMask[y*mapWidth + x] |= THRESHOLD_BIT;
	    continue;
	  }

	overThreshold = false;

	distX = 0.5 * hypot(map[y * mapWidth + x].x - map[y * mapWidth + x + 1].x,
			    map[y * mapWidth + x].y - map[y * mapWidth + x + 1].y) +
	  0.5 * hypot(map[(y + 1) * mapWidth + x].x - map[(y + 1) * mapWidth + x + 1].x,
		      map[(y + 1) * mapWidth + x].y - map[(y + 1) * mapWidth + x + 1].y);
	distY = 0.5 * hypot(map[y * mapWidth + x].x - map[(y + 1) * mapWidth + x].x,
			    map[y * mapWidth + x].y - map[(y + 1) * mapWidth + x].y) +
	  0.5 * hypot(map[y * mapWidth + x + 1].x - map[(y + 1) * mapWidth + x + 1].x,
		      map[y * mapWidth + x + 1].y - map[(y + 1) * mapWidth + x + 1].y);
	orthogonal = distY / distX;
	if (fabs(orthogonal - avgOrthogonal) / avgOrthogonal >= oThreshold)
	  overThreshold = true;
	distA = hypot(map[y * mapWidth + x].x - map[(y + 1) * mapWidth + x + 1].x,
		      map[y * mapWidth + x].y - map[(y + 1) * mapWidth + x + 1].y);
	distB = hypot(map[y * mapWidth + x + 1].x - map[(y + 1) * mapWidth + x].x,
		      map[y * mapWidth + x + 1].y - map[(y + 1) * mapWidth + x].y);
	diagonal = distB / distA;
	if (fabs(diagonal - avgDiagonal) / avgDiagonal >= dThreshold)
	  overThreshold = true;

	/* update mask */
	if (overThreshold)
	  mapMask[y * mapWidth + x] &= ~THRESHOLD_BIT;
	else
	  mapMask[y * mapWidth + x] |= THRESHOLD_BIT;
      }

  // clear all the cluster bits
  for (y = 0; y < mapHeight-1; ++y)
    for (x = 0; x < mapWidth-1; ++x)
      {
	mapMask[y * mapWidth + x] &= ~CLUSTER_BIT;
	cluster[y * mapWidth + x] = -1;
      }

  unsigned int *stack = (unsigned int *) malloc(mapHeight * mapWidth * sizeof(int));
  printf("updateMask: seeds.size = %zd\n", seeds.size());
  for (int i = 0; i < seeds.size(); ++i)
    {
      y = seeds[i] / mapWidth;
      x = seeds[i] - y * mapWidth;
      if (cluster[y * mapWidth + x] >= 0)
	continue;
      cluster[y * mapWidth + x] = i;
      int sp = 0;
      stack[sp++] = y*mapWidth+x;
      int count = 0;
      while (sp > 0)
	{
	  ++count;
	  --sp;
	  y = stack[sp] / mapWidth;
	  x = stack[sp] - y * mapWidth;
	  mapMask[y * mapWidth + x] |= CLUSTER_BIT;
	  
	  /* check up */
	  if (y > 0 && (mapMask[(y-1)*mapWidth+x] & THRESHOLD_BIT) &&
	      cluster[(y-1)*mapWidth+x] < 0)
	    {
	      cluster[(y-1)*mapWidth+x] = i;
	      stack[sp++] = (y-1)*mapWidth+x;
	    }
	  /* check right */
	  if (x < mapWidth-2 && (mapMask[y*mapWidth+x+1] & THRESHOLD_BIT) &&
	      cluster[y*mapWidth+x+1] < 0)
	    {
	      cluster[y*mapWidth+x+1] = i;
	      stack[sp++] = y*mapWidth+x+1;
	    }
	  /* check down */
	  if (y < mapHeight-2 && (mapMask[(y+1)*mapWidth+x] & THRESHOLD_BIT) &&
	      cluster[(y+1)*mapWidth+x] < 0)
	    {
	      cluster[(y+1)*mapWidth+x] = i;
	      stack[sp++] = (y+1)*mapWidth+x;
	    }
	  /* check left */
	  if (x > 0 && (mapMask[y*mapWidth+x-1] & THRESHOLD_BIT) &&
	      cluster[y*mapWidth+x-1] < 0)
	    {
	      cluster[y*mapWidth+x-1] = i;
	      stack[sp++] = y*mapWidth+x-1;
	    }
	}
      printf("flood fill marked %d elements\n", count);
    }
  free(stack);
}

void Log (const char *fmt, ...)
{
  va_list args;
  char timestamp[32];

  if (logName[0] == '\0')
    return;

  if (logFile == NULL)
    {
      logFile = fopen(logName, "w");
      if (logFile == NULL)
	Error("Could not open log file %s\n", logName);
    }

  va_start(args, fmt);
  fprintf(logFile, "%s: ", GetTimestamp(timestamp, 32));
  vfprintf(logFile, fmt, args);
  va_end(args);
  fflush(logFile);
}

void
Error (const char *fmt, ...)
{
  va_list args;
  char timestamp[32];

  if (logFile != NULL)
    {
      va_start(args, fmt);
      fprintf(logFile, "%s: ERROR: ", GetTimestamp(timestamp, 32));
      vfprintf(logFile, fmt, args);
      va_end(args);
      fflush(logFile);
    }
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  fflush(stderr);
  abort();
}

char *
GetTimestamp (char *timestamp, size_t size)
{
  struct timeval current_time;
  struct tm lt;
  size_t len;

  if (size < 24 ||
      gettimeofday(&current_time, NULL) != 0 ||
      strftime(timestamp, size, "%y%m%d %H:%M:%S",
	       localtime_r(&(current_time.tv_sec), &lt)) == 0)
    return(NULL);
  len = strlen(timestamp);
  if (size - len > 8)
    sprintf(&timestamp[len], ".%.6d", (int) current_time.tv_usec);
  return(timestamp);
}
