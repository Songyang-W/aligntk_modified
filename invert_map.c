/*
 *  invert_map.c  -  invert an image map
 *
 *  Copyright (c) 2018      Pittsburgh Supercomputing Center,
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
 *    2018  Written by Greg Hood (ghood@psc.edu)
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

#include "imio.h"
#include "invert.h"

void Error (char *fmt, ...);

int
main (int argc, char **argv)
{
  char mapName[PATH_MAX];
  MapElement *map;
  int mLevel;
  int mw, mh;
  int mxMin, myMin;
  char imgName[PATH_MAX], refName[PATH_MAX];
  float scale;

  InverseMap *invMap;

  char outputMapName[PATH_MAX];
  int oLevel;
  int omw, omh;
  int omxMin, omyMin;
  int omxMax, omyMax;
  MapElement *omap;
  float oScale;

  int x, y;
  float xv, yv, cv;
  float xp, yp;

  float xMin, xMax, yMin, yMax;

  int i;
  int error;
  char msg[PATH_MAX+256];
  unsigned char *valid;

  error = 0;
  oLevel = -1;
  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-input") == 0)
      {
	if (++i == argc)
          {
            error = 1;
            break;
          }
        strcpy(mapName, argv[i]);
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
    else if (strcmp(argv[i], "-level") == 0)
      {
        if (++i == argc || sscanf(argv[i], "%d", &oLevel) != 1)
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

      fprintf(stderr, "Usage: invert_map -input <map_name>\n");
      fprintf(stderr, "            -output <output_map_name>\n");
      fprintf(stderr, "            [-level <output_map_level>]\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (mapName[0] == '\0')
    Error("-input parameter must be specified.\n");
  if (outputMapName[0] == '\0')
    Error("-output parameter must be specified.\n");

  if (!ReadMap(mapName, &map, &mLevel, &mw, &mh,
	       &mxMin, &myMin, imgName, refName, msg))
    Error("Could not read map %s:\n%s\n", mapName, msg);

  scale = (double) (1 << mLevel);
  if (oLevel < 0)
    oLevel = mLevel;
  oScale = (double) (1 << oLevel);

  /* find the bounding box */
  xMin = 1000000000.0;
  xMax = -1000000000.0;
  yMin = 1000000000.0;
  yMax = -1000000000.0;
  for (y = 0; y < mh; ++y)
    for (x = 0; x < mw; ++x)
      {
	cv = map[y*mw + x].c;
	if (cv == 0.0)
	  continue;
	xv = map[y*mw + x].x;
	yv = map[y*mw + x].y;
	if (xv < xMin)
	  xMin = xv;
	if (xv > xMax)
	  xMax = xv;
	if (yv < yMin)
	  yMin = yv;
	if (yv > yMax)
	  yMax = yv;
      }
  if (xMin >= xMax || yMin >= yMax)
    Error("No valid map blocks found.\n");

  omxMin = (int) floorf(xMin * scale / oScale);
  omxMax = (int) ceilf(xMax * scale / oScale);
  omyMin = (int) floorf(yMin * scale / oScale);
  omyMax = (int) ceilf(yMax * scale / oScale);
  omw = omxMax - omxMin + 1;
  omh = omyMax - omyMin + 1;
  omap = (MapElement*) malloc(omw*omh*sizeof(MapElement));

  invMap = InvertMap(map, mw, mh);

  for (y = 0; y < omh; ++y)
    for (x = 0; x < omw; ++x)
      {
	xv = (x + omxMin) * oScale / scale;
	yv = (y + omyMin) * oScale / scale;
	xp = 0.0;
	yp = 0.0;
	if (Invert(invMap, &xp, &yp, xv, yv))
	  {
	    //		printf("Invert of (%f %f) produced (%f %f)\n",
	    //		       x1, y1, xp, yp);
	    omap[y*omw+x].x = xp * scale / oScale;
	    omap[y*omw+x].y = yp * scale / oScale;
	    omap[y*omw+x].c = 1.0;
	  }
	else
	  {
	    //		printf("Invert of (%f %f) failed\n", x1, y1);
	    omap[y*omw+x].x = 0.0;
	    omap[y*omw+x].y = 0.0;
	    omap[y*omw+x].c = 0.0;
	  }
    }

  // check that all map elements are valid
  valid = (unsigned char *) malloc(omh*omw);
  memset(valid, 0, omh*omw);
  for (y = 0; y < omh-1; ++y)
    for (x = 0; x < omw-1; ++x)
      if (omap[y*omw+x].c != 0.0 &&
	  omap[y*omw+x+1].c != 0.0 &&
	  omap[(y+1)*omw+x].c != 0.0 &&
	  omap[(y+1)*omw+x+1].c != 0.0)
	valid[y*omw+x] = 1;
  for (y = 0; y < omh; ++y)
    for (x = 0; x < omw; ++x)
      if (x > 0 && valid[y*omw+x-1] ||
	  y > 0 && valid[(y-1)*omw+x] ||
	  x > 0 && y > 0 && valid[(y-1)*omw+x-1] ||
	  valid[y*omw+x])
	omap[y*omw+x].c = 1.0;
      else
	omap[y*omw+x].c = 0.0;
  

  if (!WriteMap(outputMapName, omap, oLevel, omw, omh,
		omxMin, omyMin,
		refName, imgName,
		UncompressedMap, msg))
    Error("Could not write map %s:\n%s\n",
	  outputMapName, msg);

  free(map);
  FreeInverseMap(invMap);
  free(omap);
  free(valid);

  return(0);
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
