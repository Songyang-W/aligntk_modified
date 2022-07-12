/*
 * best_translation - find best translation approximating the given nonlinear map
 *
 *  This file is part of the High Throughput Electron Microscopy
 *  Lab's distribution of AlignTK.
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
 *    -
 *
 *
 *  HISTORY
 *    2022 Written by Jasper Phelps (jasper.phelps@epfl.ch)
 */

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <errno.h>

#include "imio.h"

#define LINE_LENGTH	255

int
main (int argc, char **argv)
{
  FILE *f;
  int mw, mh;
  int ix, iy;
  int x, y;
  int i;
  char msg[PATH_MAX + 1024];
  int mLevel;
  int mxMin, myMin;
  char imgName[PATH_MAX], refName[PATH_MAX];
  int ix0, iy0;
  double v;
  MapElement *map;
  int mw1, mh1;
  int n;
  float xv, yv;
  int mFactor;
  MapElement rmap[4];
  double stx, sty;
  double tx, ty;
  double xMax, yMax;
  int rLevel;
  int rFactor;
  int error;
  char inputName[PATH_MAX];
  char outputName[PATH_MAX];
  
  error = 0;
  inputName[0] = '\0';
  outputName[0] = '\0';
  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-input") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-input error\n");
	    break;
	  }
	strcpy(inputName, argv[i]);
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
    else
      {
	fprintf(stderr, "Unknown option: %s\n", argv[i]);
	error = 1;
      }

  if (error)
    {
      if (i >= argc)
	fprintf(stderr, "Incomplete option: %s\n\n", argv[i-1]);

      fprintf(stderr, "Usage: best_rigid -input in.map -output out.map\n");
      exit(1);
    }
  if (inputName[0] == '\0' || outputName[0] == '\0')
    {
      fprintf(stderr, "Both -input and -output map files must be specified.\n");
      exit(1);
    }

  if (!ReadMap(inputName, &map, &mLevel, &mw, &mh,
	       &mxMin, &myMin,
	       imgName, refName,
	       msg))
    {
      fprintf(stderr, "Could not read map %s:\n  error: %s\n",
	      inputName, msg);
      exit(1);
    }
  mFactor = 1 << mLevel;

  /* go through all map points and determine translation */
  mw1 = mw - 1;
  mh1 = mh - 1;
  n = 0;
  stx = 0.0;
  sty = 0.0;
  n = 0;
  for (y = 0; y < mh; ++y)
    for (x = 0; x < mw; ++x)
      {
	if (map[y*mw+x].c == 0.0)
	  continue;
	xv = mFactor * (x + mxMin);
	yv = mFactor * (y + myMin);
	stx += mFactor * map[y*mw+x].x - xv;
	sty += mFactor * map[y*mw+x].y - yv;
	++n;
      }
  tx = stx / n;
  ty = sty / n;

  xMax = mFactor * (mw + mxMin);
  yMax = mFactor * (mh + myMin);
  rLevel = 0;
  rFactor = 1;
  while (xMax > rFactor || yMax > rFactor)
    rFactor = 1 << (++rLevel);
  xMax = rFactor;
  yMax = rFactor;

  rmap[0].x = (tx +    0) / rFactor;
  rmap[0].y = (ty +    0) / rFactor;
  rmap[0].c = 1.0;
  rmap[1].x = (tx + xMax) / rFactor;
  rmap[1].y = (ty +    0) / rFactor;
  rmap[1].c = 1.0;
  rmap[2].x = (tx +    0) / rFactor;
  rmap[2].y = (ty + yMax) / rFactor;
  rmap[2].c = 1.0;
  rmap[3].x = (tx + xMax) / rFactor;
  rmap[3].y = (ty + yMax) / rFactor;
  rmap[3].c = 1.0;

  printf("translation is  tx = %f  ty = %f\n", tx, ty);
  
  /* write new map out */
  if (!WriteMap(outputName, rmap, rLevel, 2, 2, 0, 0,
		imgName, refName,
		UncompressedMap, msg))
    {
      fprintf(stderr, "Could not write translation map %s:\n%s\n",
	      outputName, msg);
      exit(1);
    }

  /* deallocate all map data structures */
  free(map);
}
