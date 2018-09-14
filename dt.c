/*
 * dt.c -- defines computeDistance, a function to compute the
 *         Euclidean distance transform of a binary image; algorithm
 *         adapted from A. Meijster, J.B.T.M. Roerdink and
 *         W.H. Hesselink, A general algorithm for computing distance
 *         transforms in linear time. In: Mathematical Morphology and
 *         its Apllications to Image and Signal Processing, Kluwer
 *         Acad. Publ., 2000, pp. 331-340.
 *
 *  Copyright (c) 2009-2010 National Resource for Biomedical
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
 *       NIH NCRR grant 5P41RR006009
 *
 *  HISTORY
 *    2009     Written by Greg Hood (ghood@psc.edu)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dt.h"


void
computeDistance (int type,
		 int nx, int ny,
		 unsigned char *mask,
		 float *dist)
{
  long long *g;
  long long mw;
  long long q;
  long long *s, *t;
  long long u;
  long long w;
  long long z;
  int x, y;
  long long lnx;

  // make sure there at least one marked pixel in mask; otherwise,
  //   all distances are infinity
  lnx = nx;
  mw = (nx + 7) >> 3;
  for (y = 0; y < ny; ++y)
    for (x = 0; x < mw; ++x)
      if (mask[y*mw+x] != 0)
	goto nonempty;
  for (y = 0; y < ny; ++y)
    for (x = 0; x < nx; ++x)
      dist[y*lnx+x] = 1.0e30;
  return;

 nonempty:
  g = (long long *) malloc(lnx * ny * sizeof(long long));
  s = (long long *) malloc(nx * sizeof(long long));
  t = (long long *) malloc(nx * sizeof(long long));
  for (x = 0; x < nx; ++x)
    {
      if ((mask[x >> 3] & (0x80 >> (x & 7))) != 0)
	g[x] = 0;
      else
	g[x] = nx + ny;
      for (y = 1; y < ny; ++y)
	if ((mask[y * mw + (x >> 3)] & (0x80 >> (x & 7))) != 0)
	  g[y*lnx + x] = 0;
	else
	  g[y*lnx + x] = g[(y - 1) * lnx + x] + 1;
      for (y = ny - 2; y >= 0; --y)
	if (g[(y + 1) * lnx + x] < g[y*lnx + x])
	  g[y*lnx + x] = g[(y + 1) * lnx + x] + 1;
    }

  for (y = 0; y < ny; ++y)
    {
      q = 0;
      s[0] = 0;
      t[0] = 0;
      for (u = 1; u < nx; ++u)
	{
	  switch (type)
	    {
	    case EUCLIDEAN_DISTANCE:
	    case EUCLIDEAN_DISTANCE_SQUARED:
	      while (q >= 0 &&
		     (t[q] - s[q]) * (t[q] - s[q]) +
		     g[y*lnx + s[q]] * g[y*lnx + s[q]] >
		     (t[q] - u) * (t[q] - u) +
		     g[y*lnx + u] * g[y*lnx + u])
		--q;
	      break;
	    case MANHATTAN_DISTANCE:
	      while (q >= 0 &&
		     llabs(t[q] - s[q]) + g[y*lnx + s[q]] >
		     llabs(t[q] - u) + g[y*lnx + u])
		--q;
	      break;
	    case CHESSBOARD_DISTANCE:
	      while (q >= 0)
		{
		  w = llabs(t[q] - s[q]);
		  if (g[y*lnx + s[q]] > w)
		    w = g[y*lnx + s[q]];
		  z = llabs(t[q] - u);
		  if (g[y*lnx + u] > z)
		    z = g[y*lnx + u];
		  if (w <= z)
		    break;
		  --q;
		}
	      break;
	    }
	  if (q < 0)
	    {
	      q = 0;
	      s[0] = u;
	    }
	  else
	    {
	      switch (type)
		{
		case EUCLIDEAN_DISTANCE:
		case EUCLIDEAN_DISTANCE_SQUARED:
		  w = (u * u - s[q] * s[q] +
		       g[y*lnx+u] * g[y*lnx+u] -
		       g[y*lnx+s[q]] * g[y*lnx+s[q]]) /
		    (2 * (u - s[q])) + 1;
		  break;
		case MANHATTAN_DISTANCE:
		  if (g[y*lnx + u] >= g[y*lnx + s[q]] + u - s[q])
		    w = 1000000000;
		  else if (g[y*lnx + s[q]] > g[y*lnx + u] + u - s[q])
		    w = -1000000000;
		  else
		    w = (g[y*lnx + u] - g[y*lnx + s[q]] + u + s[q]) / 2 + 1;
		  break;
		case CHESSBOARD_DISTANCE:
		  if (g[y*lnx + s[q]] <= g[y*lnx + u])
		    {
		      w = s[q] + g[y*lnx + u];
		      z = (s[q] + u) / 2;
		      if (z > w)
			w = z;
		    }
		  else
		    {
		      w = u - g[y*lnx + s[q]];
		      z = (s[q] + u) / 2;
		      if (z < w)
			w = z;
		    }
		  ++w;
		  break;
		}
	      if (w < lnx)
		{
		  ++q;
		  s[q] = u;
		  t[q] = w;
		}
	    }
	}
      for (u = nx-1; u >= 0; --u)
	{
	  switch (type)
	    {
	    case EUCLIDEAN_DISTANCE:
	      dist[y * lnx + u] = sqrt((double) ((u - s[q]) * (u - s[q]) +
						 g[y*lnx + s[q]] * g[y*lnx + s[q]]));
	      break;
	    case EUCLIDEAN_DISTANCE_SQUARED:
	      dist[y * lnx + u] = (u - s[q]) * (u - s[q]) +
		g[y*lnx + s[q]] * g[y*lnx + s[q]];
	      break;
	    case MANHATTAN_DISTANCE:
	      dist[y * lnx + u] = llabs(u - s[q]) + g[y*lnx + s[q]];
	      break;
	    case CHESSBOARD_DISTANCE:
	      w = llabs(u - s[q]);
	      if (g[y*lnx + s[q]] > w)
		w = g[y*lnx + s[q]];
	      dist[y * lnx + u] = w;
	      break;
	    }
	  if (u == t[q])
	    --q;
	}
    }
  free(g);
  free(s);
  free(t);
}
