/*
 *  compose_maps.c  -  mathematically compose two maps into a single map
 *
 *  Copyright (c) 2010-2019 National Resource for Biomedical
 *                          Supercomputing,
 *                          Pittsburgh Supercomputing Center,
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
 *  HISTORY
 *    2010  Written by Greg Hood (ghood@psc.edu)
 *    2019  Added -extrapolate option
 */

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <errno.h>
#include <float.h>

#include "imio.h"
#include "invert.h"

int Extrapolate (float *prx, float *pry, float *prc,
		 int ix, int iy, float arrx, float arry,
		 MapElement* map, int mw, int mh, float threshold,
		 float confidenceThreshold);
void RecoverFromFold (int level);
void Error (char *fmt, ...);

int
main (int argc, char **argv)
{
  int invert;
  float thresholdC;
  float extrapolationDistance;

  char map1Name[PATH_MAX];
  MapElement *map1;
  int mLevel;
  int mw, mh;
  int mxMin, myMin;
  char imgName[PATH_MAX], refName[PATH_MAX];
  
  char map2Name[PATH_MAX];
  MapElement *map2;
  int mLevel2;
  int mw2, mh2;
  int mxMin2, myMin2;
  char imgName2[PATH_MAX], refName2[PATH_MAX];
  
  InverseMap *invMap2;

  char outputMapName[PATH_MAX];
  MapElement *omap;

  int x, y;
  float x1, y1, c1;
  float xv, yv;
  float xp, yp;
  int ix, iy;
  float rx, ry, rc;
  float rx00, rx01, rx10, rx11;
  float ry00, ry01, ry10, ry11;
  float rc00, rc01, rc10, rc11;
  float rrx, rry;
  float extrapolationThreshold;

  int dx, dy;
  double weight;
  double totalWeight;
  double trx, try;

  int i;
  int error;
  char msg[PATH_MAX+256];

  error = 0;
  invert = 0;
  thresholdC = 0.0;
  extrapolationDistance = 0.0;
  map1Name[0] = '\0';
  map2Name[0] = '\0';
  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-map1") == 0)
      {
	if (++i == argc)
          {
            error = 1;
            break;
          }
	if (map1Name[0] != '\0')
	  {
	    fprintf(stderr, "Only one -map1 argument allowed.\n");
	    error = 1;
	    break;
	  }
        strcpy(map1Name, argv[i]);
      }
    else if (strcmp(argv[i], "-inverse_map2") == 0)
      {
        if (++i == argc)
          {
            error = 1;
            break;
          }
	if (map2Name[0] != '\0')
	  {
	    fprintf(stderr, "Only one of -map2 and -inverse_map2 allowed.\n");
	    error = 1;
	    break;
	  }
        strcpy(map2Name, argv[i]);
	invert = 1;
      }
    else if (strcmp(argv[i], "-map2") == 0)
      {
        if (++i == argc)
          {
            error = 1;
            break;
          }
	if (map2Name[0] != '\0')
	  {
	    fprintf(stderr, "Only one of -map2 and -inverse_map2 allowed.\n");
	    error = 1;
	    break;
	  }
	invert = 0;
        strcpy(map2Name, argv[i]);
      }
    else if (strcmp(argv[i], "-output") == 0)
      {
        if (++i == argc)
          {
            error = 1;
            break;
          }
        strcpy(outputMapName, argv[i]);
      }
    else if (strcmp(argv[i], "-threshold_c") == 0)
      {
        if (++i == argc || sscanf(argv[i], "%f", &thresholdC) != 1)
          {
            error = 1;
            break;
          }
      }
    else if (strcmp(argv[i], "-extrapolate") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%f", &extrapolationDistance) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else
      {
	error = 1;
	break;
      }

  if (error)
    {
      if (i < argc)
        fprintf(stderr, "Invalid option: %s\n", argv[i]);
      else
        fprintf(stderr, "Incomplete option: %s\n", argv[i-1]);
      fprintf(stderr, "\n");

      fprintf(stderr, "Usage: compose_maps -map1 <map_name>\n");
      fprintf(stderr, "            [-map2 <map_name>]\n");
      fprintf(stderr, "            [-inverse_map2 <map_name>]\n");
      fprintf(stderr, "            -output <output_map_name>]\n");
      fprintf(stderr, "            [-threshold_c <float>]\n");
      fprintf(stderr, "            [-extrapolate <npixels>]\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (map1Name[0] == '\0')
    Error("-map1 parameter must be specified.\n");
  if (map2Name[0] == '\0')
    Error("-map2 or -inverse_map2 parameter must be specified.\n");
  if (outputMapName[0] == '\0')
    Error("-output parameter must be specified.\n");

  /* check that not both -inverse_map2 and -extrapolate were given */
  if (invert && extrapolationDistance != 0.0)
    Error("-extrapolate cannot be used with -inverse_map2 option\n");

  if (!ReadMap(map1Name, &map1, &mLevel, &mw, &mh,
	       &mxMin, &myMin, imgName, refName, msg))
    Error("Could not read map %s:\n%s\n", map1Name, msg);

  if (!ReadMap(map2Name, &map2, &mLevel2, &mw2, &mh2,
	       &mxMin2, &myMin2, imgName2, refName2, msg))
    Error("Could not read map %s:\n%s\n", map2Name, msg);

  omap = (MapElement*) malloc(mw*mh*sizeof(MapElement));

  if (invert)
    {
      invMap2 = InvertMap(map2, mw2, mh2);

      for (y = 0; y < mh; ++y)
	for (x = 0; x < mw; ++x)
	  {
	    x1 = map1[y*mw+x].x * (1 << mLevel) / (1 << mLevel2);
	    y1 = map1[y*mw+x].y * (1 << mLevel) / (1 << mLevel2);
	    c1 = map1[y*mw+x].c;

	    // next line is temporary to correct for bug in warp3
	    if (c1 != 0.0 && c1 < 0.001)
	      c1 = 0.001;

	    if (c1 == 0.0)
	      {
		omap[y*mw+x].x = 0.0;
		omap[y*mw+x].y = 0.0;
		omap[y*mw+x].c = 0.0;
		continue;
	      }
	    xp = 0.0;
	    yp = 0.0;
	    if (Invert(invMap2, &xp, &yp, x1, y1))
	      {
		//		printf("Invert of (%f %f) produced (%f %f)\n",
		//		       x1, y1, xp, yp);
		omap[y*mw+x].x = xp * (1 << mLevel2) / (1 << mLevel);
		omap[y*mw+x].y = yp * (1 << mLevel2) / (1 << mLevel);
		omap[y*mw+x].c = c1;
	      }
	    else
	      {
		//		printf("Invert of (%f %f) failed\n", x1, y1);
		omap[y*mw+x].x = 0.0;
		omap[y*mw+x].y = 0.0;
		omap[y*mw+x].c = 0.0;
	      }
	  }
    }
  else
    {
      extrapolationThreshold = extrapolationDistance / (1 << mLevel2);
      for (y = 0; y < mh; ++y)
	for (x = 0; x < mw; ++x)
	  {
	    x1 = map1[y*mw+x].x * (1 << mLevel);
	    y1 = map1[y*mw+x].y * (1 << mLevel);
	    c1 = map1[y*mw+x].c;

	    // next line is temporary to correct for bug in warp3
	    if (c1 != 0.0 && c1 < 0.001)
	      c1 = 0.001;

	    if (c1 == 0.0)
	      {
		omap[y*mw+x].x = 0.0;
		omap[y*mw+x].y = 0.0;
		omap[y*mw+x].c = 0.0;
		continue;
	      }
	    if (c1 < 0.0)
	      Error("c1 is negative! %f\n", c1);
	    xv = x1 / (1 << mLevel2);
	    yv = y1 / (1 << mLevel2);
	    ix = ((int) floor(xv)) - mxMin2;
	    iy = ((int) floor(yv)) - myMin2;
	    rrx = xv - (mxMin2 + ix);
	    rry = yv - (myMin2 + iy);

	    if (ix < -1 || ix == -1 && rrx < 0.999 ||
		ix == mw2-1 && rrx > 0.001 || ix >= mw2 ||
		iy < -1 || iy == -1 && rry < 0.999 ||
		iy == mh2-1 && rry > 0.001 || iy >= mh2)
	      {
		if (extrapolationThreshold == 0.0 ||
		    !Extrapolate(&omap[y*mw+x].x,
				 &omap[y*mw+x].y,
				 &omap[y*mw+x].c,
				 ix, iy, rrx, rry,
				 map2, mw2, mh2,
				 extrapolationThreshold,
				 thresholdC))
		  {
		    omap[y*mw+x].x = 0.0;
		    omap[y*mw+x].y = 0.0;
		    omap[y*mw+x].c = 0.0;
		  }
		omap[y*mw+x].x *= (float) (1 << mLevel2) / (float) (1 << mLevel);
		omap[y*mw+x].y *= (float) (1 << mLevel2) / (float) (1 << mLevel);
		continue;
	      }

	    while (ix < 0)
	      {
		++ix;
		rrx -= 1.0;
	      }
	    while (iy < 0)
	      {
		++iy;
		rry -= 1.0;
	      }
	    while (ix >= mw2-1)
	      {
		--ix;
		rrx += 1.0;
	      }
	    while (iy >= mh2-1)
	      {
		--iy;
		rry += 1.0;
	      }

	    rx00 = map2[iy*mw2+ix].x;
	    ry00 = map2[iy*mw2+ix].y;
	    rc00 = map2[iy*mw2+ix].c;
	    rx01 = map2[(iy+1)*mw2+ix].x;
	    ry01 = map2[(iy+1)*mw2+ix].y;
	    rc01 = map2[(iy+1)*mw2+ix].c;
	    rx10 = map2[iy*mw2+ix+1].x;
	    ry10 = map2[iy*mw2+ix+1].y;
	    rc10 = map2[iy*mw2+ix+1].c;
	    rx11 = map2[(iy+1)*mw2+ix+1].x;
	    ry11 = map2[(iy+1)*mw2+ix+1].y;
	    rc11 = map2[(iy+1)*mw2+ix+1].c;

	    rc = rc00;
	    if (rc01 < rc)
	      rc = rc01;
	    if (rc10 < rc)
	      rc = rc10;
	    if (rc11 < rc)
	      rc = rc11;
	    if (rc < 0.0)
	      Error("rc is negative! %f\n", rc);
	    
	    if (rc == 0.0)
	      {
		if (extrapolationThreshold == 0.0 ||
		    !Extrapolate(&omap[y*mw+x].x,
				 &omap[y*mw+x].y,
				 &omap[y*mw+x].c,
				 ix, iy, rrx, rry,
				 map2, mw2, mh2,
				 extrapolationThreshold,
				 thresholdC))
		  {
		    omap[y*mw+x].x = 0.0;
		    omap[y*mw+x].y = 0.0;
		    omap[y*mw+x].c = 0.0;
		  }
		omap[y*mw+x].x *= (float) (1 << mLevel2) / (float) (1 << mLevel);
		omap[y*mw+x].y *= (float) (1 << mLevel2) / (float) (1 << mLevel);
		continue;
	      }

	    rx = rx00 * (rrx - 1.0) * (rry - 1.0)
	      - rx10 * rrx * (rry - 1.0) 
	      - rx01 * (rrx - 1.0) * rry
	      + rx11 * rrx * rry;
	    ry = ry00 * (rrx - 1.0) * (rry - 1.0)
	      - ry10 * rrx * (rry - 1.0) 
	      - ry01 * (rrx - 1.0) * rry
	      + ry11 * rrx * rry;
	    xp = rx * (1 << mLevel2);
	    yp = ry * (1 << mLevel2);

	    if (isnan(xp) || isnan(yp))
	      Error("nan in xp or yp at (%d,%d).\n", x, y);
	    if (isnan(rc) || isnan(c1))
	      Error("nan in input maps at (%d,%d).\n", x, y);
	    omap[y*mw+x].x = xp / (1 << mLevel);
	    omap[y*mw+x].y = yp / (1 << mLevel);
	    if (c1 < rc)
	      rc = c1;
	    omap[y*mw+x].c = rc;
	  }
    }

  if (thresholdC != 0.0)
    for (y = 0; y < mh; ++y)
      for (x = 0; x < mw; ++x)
	if (omap[y*mw+x].c >= thresholdC)
	  omap[y*mw+x].c = 1.0;
	else
	  omap[y*mw+x].c = 0.0;

  if (!WriteMap(outputMapName, omap, mLevel, mw, mh,
		mxMin, myMin,
		imgName, invert ? imgName2 : refName2,
		UncompressedMap, msg))
    Error("Could not write map %s:\n%s\n",
	  outputMapName, msg);

  free(map1);
  free(map2);
  if (invert)
    FreeInverseMap(invMap2);
  free(omap);

  return(0);
}

int
Extrapolate (float *prx, float *pry, float *prc,
	     int ix, int iy, float arrx, float arry,
	     MapElement* map, int mw, int mh, float threshold,
	     float confidenceThreshold)
{
  float rx = 0.0;
  float ry = 0.0;
  float rc = 0.0;
  float totalWeight = 0.0;
  int maxDelta = ((int) ceil(threshold)) + 2;
  float rrx, rry;
  int dx, dy;
  int ixv, iyv;
  float rx00, rx01, rx10, rx11;
  float ry00, ry01, ry10, ry11;
  float rc00, rc01, rc10, rc11;
  float weight;
  float rcMin;
  float d;
  float xn, yn;
  float xf, yf;
  float closestDistance;

  closestDistance = FLT_MAX;
  for (dy = -maxDelta; dy <= maxDelta; ++dy)
    {
      iyv = iy + dy;
      if (iyv < 0 || iyv >= mh-1)
	continue;
      rry = arry - dy;

      for (dx = -maxDelta; dx <= maxDelta; ++dx)
	{
	  ixv = ix + dx;
	  if (ixv < 0 || ixv >= mw-1 ||
	      map[iyv*mw+ixv].c <= confidenceThreshold ||
	      map[(iyv+1)*mw+ixv].c <= confidenceThreshold ||
	      map[iyv*mw+ixv+1].c <= confidenceThreshold ||
	      map[(iyv+1)*mw+ixv+1].c <= confidenceThreshold)
	    continue;

	  rrx = arrx - dx;

	  rx00 = map[iyv*mw+ixv].x;
	  ry00 = map[iyv*mw+ixv].y;
	  rc00 = map[iyv*mw+ixv].c;
	  rx01 = map[(iyv+1)*mw+ixv].x;
	  ry01 = map[(iyv+1)*mw+ixv].y;
	  rc01 = map[(iyv+1)*mw+ixv].c;
	  rx10 = map[iyv*mw+ixv+1].x;
	  ry10 = map[iyv*mw+ixv+1].y;
	  rc10 = map[iyv*mw+ixv+1].c;
	  rx11 = map[(iyv+1)*mw+ixv+1].x;
	  ry11 = map[(iyv+1)*mw+ixv+1].y;
	  rc11 = map[(iyv+1)*mw+ixv+1].c;

	  if (dx == 0 && dy == 0)
	    weight = M_SQRT2;
	  else
	    weight = 1.0 / hypotf((float) dx, (float) dy);

	  // compute d, the weight to be given to the far-field vs
	  //    near-field extrapolation;   d is roughly the distance
	  //    from the map element to the point to be extrapolated,
	  //    but bounded between 0.0 and 1.0
	  if (rrx >= 1.0)
	    if (rry >= 1.0)
	      if (rrx >= 2.0 || rry >= 2.0)
		d = 1.0;
	      else
		d = hypotf(rrx - 1.0, rry - 1.0);
	    else if (rry < 0.0)
	      if (rrx >= 2.0 || rry <= -1.0)
		d = 1.0;
	      else
		d = hypotf(rrx - 1.0, rry);
	    else
	      d = rrx - 1.0;
	  else if (rrx < 0.0)
	    if (rry >= 1.0)
	      if (rrx <= -1.0 || rry >= 2.0)
		d = 1.0;
	      else
		d = hypotf(rrx, rry - 1.0);
	    else if (rry < 0.0)
	      if (rrx <= -1.0 || rry <= -1.0)
		d = 1.0;
	      else
		d = hypotf(rrx, rry);
	    else
	      d = -rrx;
	  else
	    if (rry >= 1.0)
	      d = rry - 1.0;
	    else if (rry < 0.0)
	      d = -rry;
	    else
	      d = 0.0;
	  if (d > 1.0)
	    d = 1.0;
	  // near-field extrapolation
	  if (d < 1.0)
	    {
	      xn = rx00 * (rrx - 1.0) * (rry - 1.0)
		- rx10 * rrx * (rry - 1.0) 
		- rx01 * (rrx - 1.0) * rry
		+ rx11 * rrx * rry;
	      yn = ry00 * (rrx - 1.0) * (rry - 1.0)
		- ry10 * rrx * (rry - 1.0) 
		- ry01 * (rrx - 1.0) * rry
		+ ry11 * rrx * rry;
	    }
	  else
	    xn = yn = 0.0;
	  // far-field extrapolation
	  if (d > 0.0)
	    {
	      xf = 0.5 * ((rx10 - rx00) + (rx11 - rx01)) * (rrx - 0.5) +
		0.5 * ((rx01 - rx00) + (rx11 - rx10)) * (rry - 0.5) +
		0.25 * (rx00 + rx10 + rx01 + rx11);
	      yf = 0.5 * ((ry01 - ry00) + (ry11 - ry10)) * (rry - 0.5) +
		0.5 * ((ry10 - ry00) + (ry11 - ry01)) * (rrx - 0.5) +
		0.25 * (ry00 + ry10 + ry01 + ry11);
	    }
	  else
	    xf = yf = 0.0;
	  // actual extrapolation is a weighted average of the near- and
	  //    far-field extrapolations
	  rx += weight * ((1.0 - d) * xn + d * xf);
	  ry += weight * ((1.0 - d) * yn + d * yf);
	  rcMin = rc00;
	  if (rc01 < rcMin)
	    rcMin = rc01;
	  if (rc10 < rcMin)
	    rcMin = rc10;
	  if (rc11 < rcMin)
	    rcMin = rc11;
	  rc += weight * rcMin;

	  totalWeight += weight;

	  // compute d, the distance from the point to this map element
	  if (rrx >= 1.0)
	    if (rry >= 1.0)
	      d = hypotf(rrx - 1.0, rry - 1.0);
	    else if (rry < 0.0)
	      d = hypotf(rrx - 1.0, rry);
	    else
	      d = rrx - 1.0;
	  else if (rrx < 0.0)
	    if (rry >= 1.0)
	      d = hypotf(rrx, rry - 1.0);
	    else if (rry < 0.0)
	      d = hypotf(rrx, rry);
	    else
	      d = -rrx;
	  else
	    if (rry >= 1.0)
	      d = rry - 1.0;
	    else if (rry < 0.0)
	      d = -rry;
	    else
	      d = 0.0;

	  if (d < closestDistance)
	    closestDistance = d;
	}
    }
  if (totalWeight == 0.0 ||
      closestDistance > threshold)
    return(0);
  *prx = rx / totalWeight;
  *pry = ry / totalWeight;
  *prc = rc / totalWeight;
  return(1);
}

void Error (char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  fflush(stderr);
  abort();
}
