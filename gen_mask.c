/*
 * gen_mask.c  - generate mask from an image
 *
 *  Copyright (c) 2008-2022 Pittsburgh Supercomputing Center,
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
 *    2008-2022  Written by Greg Hood (ghood@psc.edu)
 */

#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>

#include "imio.h"

#define THRESHOLD_METHOD	0
#define RANGE_METHOD		1
#define BOUNDARY_FILL_METHOD	2

/* FORWARD DECLARATIONS */
void FloodFill (unsigned char *image, unsigned char *mask,
		unsigned char background,
		uint32_t x, uint32_t y, uint64_t w, uint64_t h);
uint64_t ClusterFloodFill (uint64_t w, uint64_t h,
			   int64_t *cluster,
			   unsigned char *image, unsigned char *mask,
			   int64_t id, uint32_t ix, uint32_t iy);
void Error (char *fmt, ...);


int
main (int argc, char **argv)
{
  int error;
  int method;
  int threshold;
  int lowerThreshold, upperThreshold;
  int erode;
  float clusterThreshold;  // a percentage
  char inputName[PATH_MAX];
  char outputName[PATH_MAX];
  char msg[PATH_MAX+256];
  int iw, ih;
  uint64_t w, h;
  int i;
  uint32_t x, y;
  uint64_t bpl;
  uint64_t count[256];
  unsigned char *newMask;
  unsigned char *tmpMask;
  unsigned char background = 0;
  unsigned char *image = 0;
  unsigned char *mask = 0;
  unsigned char *bitMask = 0;
  uint64_t nClusters = 0;
  int64_t *cluster = 0;
  uint64_t *clusterCount = 0;
  unsigned char *useCluster = 0;
  uint64_t totalCount;

  error = 0;
  method = THRESHOLD_METHOD;
  threshold = 1;
  lowerThreshold = 0;
  upperThreshold = 255;
  erode = 0;
  clusterThreshold = 0.0;
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
    else if (strcmp(argv[i], "-threshold") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &threshold) != 1)
	  {
	    fprintf(stderr, "-threshold error\n");
	    error = 1;
	    break;
	  }
	method = THRESHOLD_METHOD;
      }
    else if (strcmp(argv[i], "-range") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d-%d",
				  &lowerThreshold,
				  &upperThreshold) != 2)
	  {
	    fprintf(stderr, "-range error\n");
	    error = 1;
	    break;
	  }
	method = RANGE_METHOD;
      }
    else if (strcmp(argv[i], "-boundary-fill") == 0)
      method = BOUNDARY_FILL_METHOD;
    else if (strcmp(argv[i], "-erode") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &erode) != 1)
	  {
	    fprintf(stderr, "-erode error\n");
	    error = 1;
	    break;
	  }
      }
    else
      {
	fprintf(stderr, "Unknown option: %s\n", argv[i]);
	error = 1;
	break;
      }

  if (argc == 1 || error)
    {
      if (argc > 1 && i >= argc)
	fprintf(stderr, "Incomplete option: %s\n\n", argv[i-1]);

      fprintf(stderr, "Usage: gen_mask -input image.tif -output mask.pbm\n");
      fprintf(stderr, "               [-threshold int_value ]\n");
      fprintf(stderr, "               [-range lower_int_value-upper_int_value ]\n");
      fprintf(stderr, "               [-boundary-fill]\n");
      fprintf(stderr, "               [-erode int_iterations]\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (inputName[0] == '\0' || outputName[0] == '\0')
    {
      fprintf(stderr, "Both -input and -output parameters must be specified.\n");
      exit(1);
    }


  if (!ReadImage(inputName, &image, &iw, &ih, -1, -1, -1, -1, msg))
    Error("%s\n", msg);
  w = iw;
  h = ih;
  mask = (unsigned char *) malloc(h*w);

  if (method == THRESHOLD_METHOD)
    {
      for (y = 0; y < ih; ++y)
	for (x = 0; x < iw; ++x)
	  mask[y*w + x] = image[y*w + x] >= threshold;
    }
  else if (method == RANGE_METHOD)
    {
      for (y = 0; y < ih; ++y)
	for (x = 0; x < iw; ++x)
	  mask[y*w + x] =
	    image[y*w + x] >= lowerThreshold &&
	    image[y*w + x] <= upperThreshold;
    }
  else if (method == BOUNDARY_FILL_METHOD)
    {
      for (x = 0; x < iw; ++x)
	{
	  y = 0;
	  ++count[image[y*w + x]];
	  y = ih-1;
	  ++count[image[y*w + x]];
	}
      for (y = 1; y < ih-1; ++y)
	{
	  x = 0;
	  ++count[image[y*w + x]];
	  x = iw-1;
	  ++count[image[y*w + x]];
	}

      /* find the most common value and consider that the background */
      background = 0;
      for (i = 1; i < 256; ++i)
	if (count[i] > count[background])
	  background = i;

      /* flood fill in from all perimeter pixels with that value */
      memset(mask, 1, h*w);
      for (x = 0; x < iw; ++x)
	{
	  y = 0;
	  if (image[y*w + x] == background && mask[y*w + x])
	    FloodFill(image, mask, background, x, y, w, h);
	  y = ih-1;
	  if (image[y*w + x] == background && mask[y*w + x])
	    FloodFill(image, mask, background, x, y, w, h);
	}
      for (y = 1; y < ih-1; ++y)
	{
	  x = 0;
	  if (image[y*w + x] == background && mask[y*w + x])
	    FloodFill(image, mask, background, x, y, w, h);
	  x = iw-1;
	  if (image[y*w + x] == background && mask[y*w + x])
	    FloodFill(image, mask, background, x, y, w, h);
	}
    }

  if (erode > 0)
    newMask = (unsigned char *) malloc(h * w);
  for (i = 0; i < erode; ++i)
    {
      for (y = 0; y < ih; ++y)
        for (x = 0; x < iw; ++x)
          newMask[y*w+x] = mask[y*w+x] &&
            (x == 0 || mask[y*w+x-1]) &&
            (x >= iw-1 || mask[y*w+x+1]) &&
            (y == 0 || mask[(y-1)*w+x]) &&
            (y >= ih-1 || mask[(y+1)*w+x]);
      tmpMask = mask;
      mask = newMask;
      newMask = tmpMask;
    }
  if (erode > 0)
    free(newMask);

  /* use only pixels from clusters that are bigger than the
     given percentage threshold */
  if (clusterThreshold > 0.0)
    {
      // clear all the cluster bits
      cluster = (int64_t *) malloc(h * w * sizeof(int64_t));
      for (y = 0; y < h; ++y)
	for (x = 0; x < w; ++x)
	  cluster[y * w + x] = -1;

      clusterCount = (uint64_t *) malloc(h * w * sizeof(uint64_t));
      useCluster = (unsigned char *) malloc(h * w * sizeof(unsigned char));
      nClusters = 0;
      totalCount = 0;
      for (y = 0; y < ih; ++y)
	for (x = 0; x < iw; ++x)
	  {
	    if (!mask[y * w + x] || cluster[y * w + x] >= 0)
	      continue;
	    clusterCount[nClusters] = ClusterFloodFill(w, h, cluster, image, mask, nClusters, x, y);
	    printf("Cluster %llu has %llu elements\n", nClusters, count[nClusters]);
	    totalCount += clusterCount[nClusters];
	    ++nClusters;
	  }

      for (i = 0; i < nClusters; ++i)
	useCluster[i] = (100.0 * count[nClusters]) / totalCount >= clusterThreshold;

      for (y = 0; y < ih; ++y)
	for (x = 0; x < iw; ++x)
	  if (mask[y*w+x] &&
	      (cluster[y*w+x] < 0 || !useCluster[cluster[y*w+x]]))
	    mask[y*w+x] = 0;

      free(useCluster);
      free(clusterCount);
      free(cluster);
    }

  /* write out the mask */
  bpl = (w+7) >> 3;
  bitMask = (unsigned char *) malloc(h * bpl);
  memset(bitMask, 0, h * bpl);
  for (y = 0; y < ih; ++y)
    for (x = 0; x < iw; ++x)
      if (mask[y*w+x])
	bitMask[y * bpl + (x >> 3)] |= 0x80 >> (x & 7);
  free(mask);
  if (!WriteBitmap(outputName, bitMask,
		   iw, ih, UncompressedBitmap, msg))
    Error("%s\n", msg);
  free(bitMask);
  exit(0);
}

void
FloodFill (unsigned char *image, unsigned char *mask,
	   unsigned char background,
	   uint32_t x, uint32_t y, uint64_t w, uint64_t h)
{
  int32_t y1;
  uint32_t left, right;
  uint64_t sp;
  uint64_t *stack;

  sp = 0;
  stack = (uint64_t *) malloc(h * w * sizeof(uint64_t));
  stack[sp++] = y*w+x;
  while (sp > 0)
    {
      --sp;
      y = stack[sp] / w;
      x = stack[sp] - y * w;
      y1 = y;
      while (y1 >= 0 &&
             image[y1*w + x] == background &&
             mask[y1*w + x])
        --y1;
      ++y1;
      left = 0;
      right = 0;
      while (y1 < h &&
             image[y1*w + x] == background &&
             mask[y1*w + x])
        {
          mask[y1*w + x] = 0;
          if (!left &&
              x > 0 &&
	      image[y1*w + (x-1)] == background &&
	      mask[y1*w + (x-1)])
            {
              stack[sp++] = y1 * w + (x-1);
              left = 1;
            }
          else if (left &&
                   x > 0 &&
                   (image[y1*w + (x-1)] != background ||
                    !mask[y1*w + (x-1)]))
	    left = 0;

          if (!right &&
              x < w-1 &&
              image[y1*w + (x+1)] == background &&
              mask[y1*w + (x+1)])
            {
              stack[sp++] = y1 * w + (x+1);
              right = 1;
            }
          else if (right &&
                   x < w-1 &&
                   (image[y1*w + (x+1)] != background ||
                    !mask[y1*w + (x+1)]))
            right = 0;
          ++y1;
        }
    }
  free(stack);
}

uint64_t
ClusterFloodFill (uint64_t w, uint64_t h,
		  int64_t *cluster, unsigned char *image, unsigned char *mask,
		  int64_t id, uint32_t ix, uint32_t iy)
{
  uint64_t *stack = (uint64_t *) malloc(h * w * sizeof(uint64_t));
  uint64_t sp = 0;
  uint64_t count = 0;
  uint32_t x, y;
  cluster[iy * w + ix] = id;
  stack[sp++] = iy*w+ix;
  while (sp > 0)
    {
      ++count;
      --sp;
      y = stack[sp] / w;
      x = stack[sp] - y * w;
      
      /* check up */
      if (y > 0 && mask[(y-1)*w+x] && cluster[(y-1)*w+x] < 0)
        {
          cluster[(y-1)*w+x] = id;
          stack[sp++] = (y-1)*w+x;
        }
      /* check right */
      if (x < w-2 && mask[y*w+x+1] && cluster[y*w+x+1] < 0)
        {
          cluster[y*w+x+1] = id;
          stack[sp++] = y*w+x+1;
        }
      /* check down */
      if (y < h-2 && mask[(y+1)*w+x] && cluster[(y+1)*w+x] < 0)
        {
          cluster[(y+1)*w+x] = id;
          stack[sp++] = (y+1)*w+x;
        }
      /* check left */
      if (x > 0 && mask[y*w+x-1] && cluster[y*w+x-1] < 0)
        {
          cluster[y*w+x-1] = id;
          stack[sp++] = y*w+x-1;
        }
    }
  free(stack);
  return(count);
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
