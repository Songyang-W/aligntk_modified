/*
 *  apply_map.c  -  applies maps to images
 *
 *  Copyright (c) 2009-2013 Pittsburgh Supercomputing Center,
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
 *    2009  Written by Greg Hood (ghood@psc.edu)
 *    2011  Fixed out-of-range array indices when rendering multi-gigabyte
 *             images (ghood@psc.edu)
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>
#include <stdarg.h>
#include <errno.h>
#include <sched.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "imio.h"
#include "invert.h"
#include "dt.h"
#include <iostream>
#include <fstream>

#define LINE_LENGTH   255
#define MAX_LABEL_LENGTH  255
#define QUOTE(str)    #str
#define EXPAND_AND_QUOTE(str) QUOTE(str)

/* TYPES */
typedef struct Image
{
        int next; /* index of next image in hash bucket */
        char *name; /* name of this image */
        int width, height; /* width and height in pixels */
        unsigned char *image; /* the image bytes (NULL if not loaded) */
        unsigned char *mask; /* the image mask, if present */
        unsigned char *dist; /* the distance array -- each element holds the distance
                                from the corresponding image pixel to the closest
                                masked pixel (in units of 1/64 pixels) */
        float minX, maxX; /* bounds of the area this image covers in the */
        float minY, maxY; /*   final image */
        int needed; /* true if this image is needed for the current set of
                       tiles */
        int mapBytes;   /* size of maps in bytes */
        int mLevel; /* map level */
        int mw, mh; /* map width and height */
        int mxMin, myMin; /* map offset in x and y */
        MapElement *map; /* map of this image into the final image */
        InverseMap *invMap; /* inverse map that translates points in the final image
                               into points in this image */
        /* intensity map is optional */
        int imapLevel; /* intensity map level */
        int imapw, imaph; /* intensity map width and height */
        MapElement *imap; /* intensity map */
        MapElement *targetMap;/* map of where pixels from this image ended up
                                 in target image */
        time_t mtime;   /* the modification time for this image or its map */
} Image;

/* GLOBAL VARIABLES */
int resume = 0;
char imageListName[PATH_MAX];
char imageName[PATH_MAX];
char imagesName[PATH_MAX];
char extension[PATH_MAX];
char masksName[PATH_MAX];
char fontFileName[PATH_MAX];
char mapsName[PATH_MAX];
char imapsName[PATH_MAX];
char outputName[PATH_MAX];
char sourceMapName[PATH_MAX];
char targetMapsName[PATH_MAX];
int overlay = 0;
int blend = 1;
int margin = -1;
int tileWidth = -1;
int tileHeight = -1;
int memoryLimit = 1024;   /* 1 GB */
int tree = 0;
float blackValue = 0.0;
float whiteValue = 255.0;
float mapScale = 1.0;
float imapScale = 1.0;
float maskScale = 1.0;
int compress = 0;
int reductionFactor = 1;
int sourceMapLevel = 6;
int targetMapsLevel = 6;
int update = 0;
float rotation = 0.0;
float rotationX = 0.0;
float rotationY = 0.0;

int nImages = 0;
Image *images = 0;
int *imageHashTable = 0;
unsigned char *canvas = 0;
unsigned short *weight = 0;
size_t canvasSize = 0;
size_t canvasWidth, canvasHeight;
size_t oldCanvasWidth = 0;
int canvasMinX, canvasMinY;
size_t weightWidth, weightHeight;
int weightMinX, weightMinY;
size_t tw, th;
int oMinX, oMinY, oMaxX, oMaxY;
size_t oWidth, oHeight;
float range;
unsigned char *out;
size_t outWidth, outHeight;
int outMinY;
int nProcessed = 0;
size_t imageMem = 0;
int sourceMapSize;
MapElement *sourceMap = 0;
size_t sourceMapWidth, sourceMapHeight;
int sourceMapFactor;
int sourceMapMask;
int targetMapsFactor;
char label[MAX_LABEL_LENGTH+1];
char labelName[MAX_LABEL_LENGTH+1];
unsigned char *font = 0;
int fontWidth, fontHeight;
float cosRot, sinRot;

#define DIR_HASH_SIZE 8192
char *dirHash[DIR_HASH_SIZE];

/* FORWARD DECLARATIONS */
void PaintImage (int i, int minX, int maxX, int minY, int maxY);
void WriteTiles (int col, int startRow, int endRow, char *iName);
void Error (char *fmt, ...);
unsigned int Hash (char *s);
int CreateDirectories (char *fn);
void PrintUsage ();

int
main (int argc, char **argv, char **envp)
{
        int i, j;
        int n;
        int error;
        char fn[PATH_MAX];
        MapElement *map;
        int x, y;
        float minX, minY, maxX, maxY;
        int iMinX, iMinY, iMaxX, iMaxY;
        float xv, yv;
        float spacing;
        float ax, bx, ay, by;
        int mLevel;
        int mw, mh;
        int mxMin, myMin;
        char imName0[PATH_MAX];
        char imName1[PATH_MAX];
        char msg[PATH_MAX+256];
        char msg2[PATH_MAX+256];
        FILE *f;
        int imagesSize;
        char line[LINE_LENGTH+1];
        int dir;
        int nItems;
        int width, height;
        int rows, cols;
        int ixv, iyv;
        float rx, ry;
        float rrx, rry;
        float rx00, rx01, rx10, rx11, ry00, ry01, ry10, ry11;
        float rv;
        cpu_set_t cpumask;
        int nOutputImages;
        int startImage, endImage;
        int imapLevel;
        int imapw, imaph;
        int imapXMin, imapYMin;
        char imapName0[PATH_MAX], imapName1[PATH_MAX];
        MapElement *imap;
        int regionWidth, regionHeight, regionOffsetX, regionOffsetY;
        int labelWidth, labelHeight, labelOffsetX, labelOffsetY;
        size_t *increase, *decrease;
        size_t memoryRequired;
        size_t maxMemoryRequired;
        int startX, endX;
        int startY, endY;
        int tx, ty;
        int tCol;
        int ix, iy;
        int dx, dy;
        int cx;
        int nx, ny;
        int oi;
        int hs;
        int hi;
        int startRow, endRow;
        int iv;
        int sum;
        int offset;
        int a;
        int valid;
        struct stat sb;
        int updateOutput;
        int labelMinX, labelMaxX, labelMinY, labelMaxY;
        int charMinX, charMaxX;
        int ci;
        float sx, ex, sy, ey;
        int isx, iex, isy, iey;
        float v;
        float wx, wy, ws;
        int len;
        int ilv;
        float rxp, ryp;

        error = 0;
        imageListName[0] = '\0';
        imageName[0] = '\0';
        imagesName[0] = '\0';
        masksName[0] = '\0';
        strcpy(fontFileName, EXPAND_AND_QUOTE(FONT_FILE));
        mapsName[0] = '\0';
        imapsName[0] = '\0';
        outputName[0] = '\0';
        sourceMapName[0] = '\0';
        targetMapsName[0] = '\0';
        labelName[0] = '\0';
        strcpy(extension, "tif");
        regionWidth = regionHeight = regionOffsetX = regionOffsetY = -1;
        labelWidth = labelHeight = labelOffsetX = labelOffsetY = -1;
        for (i = 1; i < argc; ++i)
                if (strcmp(argv[i], "-image_list") == 0)
                {
                        if (++i == argc)
                        {
                                error = 1;
                                break;
                        }
                        strcpy(imageListName, argv[i]);
                }
                else if (strcmp(argv[i], "-image") == 0)
                {
                        if (++i == argc)
                        {
                                error = 1;
                                break;
                        }
                        strcpy(imageName, argv[i]);
                }
                else if (strcmp(argv[i], "-images") == 0)
                {
                        if (++i == argc)
                        {
                                error = 1;
                                break;
                        }
                        strcpy(imagesName, argv[i]);
                }
                else if (strcmp(argv[i], "-pgm") == 0)
                        strcpy(extension, "pgm");
                else if (strcmp(argv[i], "-masks") == 0)
                {
                        if (++i == argc)
                        {
                                error = 1;
                                break;
                        }
                        strcpy(masksName, argv[i]);
                }
                else if (strcmp(argv[i], "-mask_scale") == 0)
                {
                        if (++i == argc || sscanf(argv[i], "%f", &maskScale) != 1)
                        {
                                error = 1;
                                break;
                        }
                }
                else if (strcmp(argv[i], "-maps") == 0)
                {
                        if (++i == argc)
                        {
                                error = 1;
                                break;
                        }
                        strcpy(mapsName, argv[i]);
                }
                else if (strcmp(argv[i], "-map_scale") == 0)
                {
                        if (++i == argc || sscanf(argv[i], "%f", &mapScale) != 1)
                        {
                                error = 1;
                                break;
                        }
                }
                else if (strcmp(argv[i], "-imaps") == 0)
                {
                        if (++i == argc)
                        {
                                error = 1;
                                break;
                        }
                        strcpy(imapsName, argv[i]);
                }
                else if (strcmp(argv[i], "-imap_scale") == 0)
                {
                        if (++i == argc || sscanf(argv[i], "%f", &imapScale) != 1)
                        {
                                error = 1;
                                break;
                        }
                }
                else if (strcmp(argv[i], "-output") == 0)
                {
                        if (++i == argc)
                        {
                                error = 1;
                                break;
                        }
                        strcpy(outputName, argv[i]);
                }
                else if (strcmp(argv[i], "-overlay") == 0)
                        overlay = 1;
                else if (strcmp(argv[i], "-blend") == 0)
                        blend = 1;
                else if (strcmp(argv[i], "-mosaic") == 0)
                        blend = 0;
                else if (strcmp(argv[i], "-margin") == 0)
                {
                        if (++i == argc || sscanf(argv[i], "%d", &margin) != 1)
                        {
                                error = 1;
                                break;
                        }
                }
                else if (strcmp(argv[i], "-tile") == 0)
                {
                        if (++i == argc || sscanf(argv[i], "%dx%d",
                                                  &tileWidth, &tileHeight) != 2)
                        {
                                error = 1;
                                break;
                        }
                }
                else if (strcmp(argv[i], "-memory") == 0)
                {
                        if (++i == argc || sscanf(argv[i], "%d", &memoryLimit) != 1)
                        {
                                error = 1;
                                break;
                        }
                }
                else if (strcmp(argv[i], "-tree") == 0)
                        tree = 1;
                else if (strcmp(argv[i], "-black") == 0)
                {
                        if (++i == argc || sscanf(argv[i], "%f", &blackValue) != 1)
                        {
                                error = 1;
                                break;
                        }
                }
                else if (strcmp(argv[i], "-white") == 0)
                {
                        if (++i == argc || sscanf(argv[i], "%f", &whiteValue) != 1)
                        {
                                error = 1;
                                break;
                        }
                }
                else if (strcmp(argv[i], "-compress") == 0)
                        compress = 1;
                else if (strcmp(argv[i], "-resume") == 0)
                        resume = 1;
                else if (strcmp(argv[i], "-region") == 0)
                {
                        if (++i == argc || sscanf(argv[i], "%dx%d%d%d",
                                                  &regionWidth, &regionHeight,
                                                  &regionOffsetX, &regionOffsetY) != 4)
                        {
                                error = 1;
                                break;
                        }
                }
                else if (strcmp(argv[i], "-reduction") == 0)
                {
                        if (++i == argc || sscanf(argv[i], "%d", &reductionFactor) != 1)
                        {
                                error = 1;
                                break;
                        }
                }
                else if (strcmp(argv[i], "-source_map") == 0)
                {
                        if (++i == argc)
                        {
                                error = 1;
                                break;
                        }
                        strcpy(sourceMapName, argv[i]);
                }
                else if (strcmp(argv[i], "-source_map_level") == 0)
                {
                        if (++i == argc || sscanf(argv[i], "%d", &sourceMapLevel) != 1)
                        {
                                error = 1;
                                break;
                        }
                }
                else if (strcmp(argv[i], "-target_maps") == 0)
                {
                        if (++i == argc)
                        {
                                error = 1;
                                break;
                        }
                        strcpy(targetMapsName, argv[i]);
                }
                else if (strcmp(argv[i], "-target_maps_level") == 0)
                {
                        if (++i == argc || sscanf(argv[i], "%d", &targetMapsLevel) != 1)
                        {
                                error = 1;
                                break;
                        }
                }
                else if (strcmp(argv[i], "-update") == 0)
                        update = 1;
                else if (strcmp(argv[i], "-label") == 0)
                {
                        if (++i == argc)
                        {
                                error = 1;
                                break;
                        }
                        strcpy(labelName, argv[i]);
                }
                else if (strcmp(argv[i], "-label_location") == 0)
                {
                        if (++i == argc || sscanf(argv[i], "%dx%d%d%d",
                                                  &labelWidth, &labelHeight,
                                                  &labelOffsetX, &labelOffsetY) != 4)
                        {
                                error = 1;
                                break;
                        }
                }
                else if (strcmp(argv[i], "-rotation") == 0)
                {
                        if (++i == argc || sscanf(argv[i], "%f", &rotation) != 1)
                        {
                                error = 1;
                                break;
                        }
                        cosRot = cos(rotation * M_PI / 180.0);
                        sinRot = sin(rotation * M_PI / 180.0);
                }
                else if (strcmp(argv[i], "-rotation_center") == 0)
                {
                        if (++i == argc || sscanf(argv[i], "%f,%f",
                                                  &rotationX, &rotationY) != 2)
                        {
                                error = 1;
                                break;
                        }
                }
                else if (strcmp(argv[i], "-font") == 0)
                {
                        if (++i == argc)
                        {
                                error = 1;
                                break;
                        }
                        strcpy(fontFileName, argv[i]);
                }
                else error = 1;

        if (error)
        {
                fprintf(stderr, "Usage: apply_map -image_list list_file -images image_prefix\n");
                fprintf(stderr, "              -maps map_prefix -output file_prefix\n");
                fprintf(stderr, "              [-map_scale scaling_factor]\n");
                fprintf(stderr, "              [-masks mask_prefix]\n");
                fprintf(stderr, "              [-masks_scale scaling_factor]\n");
                fprintf(stderr, "              [-imaps imaps_prefix]\n");
                fprintf(stderr, "              [-imap_scale scaling_factor]\n");
                fprintf(stderr, "              [-blend]\n");
                fprintf(stderr, "              [-mosaic]\n");
                fprintf(stderr, "              [-margin margin_thickness_in_pixels]\n");
                fprintf(stderr, "              [-tile tilewidthxtileheight]\n");
                fprintf(stderr, "              [-memory memory_limit_in_MB]\n");
                fprintf(stderr, "              [-tree]\n");
                fprintf(stderr, "              [-black black_value]\n");
                fprintf(stderr, "              [-white white_value]\n");
                fprintf(stderr, "              [-compress]\n");
                fprintf(stderr, "              [-region WxH+X+Y]\n");
                fprintf(stderr, "              [-reduction reduction_factor]\n");
                fprintf(stderr, "              [-resume]\n");
                fprintf(stderr, "              [-rotation CCW_rotation_in_degrees]\n");
                fprintf(stderr, "              [-rotation_center x,y]\n");
                fprintf(stderr, "              [-source_map map_name]\n");
                fprintf(stderr, "              [-source_map_level level]\n");
                fprintf(stderr, "              [-target_maps map_name_prefix]\n");
                fprintf(stderr, "              [-target_maps_level level]\n");
                fprintf(stderr, "              [-update]\n");
                fprintf(stderr, "              [-label WxH+Y+Y]\n");
                exit(1);
        }

        /* check that at least minimal parameters were supplied */
        if ((imageListName[0] == '\0' && imageName[0] == '\0') ||
            imagesName[0] == '\0' ||
            mapsName[0] == '\0' ||
            outputName[0] == '\0')
        {
                fprintf(stderr, "-image_list (or -image), -images, -maps, and -output parameters must be specified.\n");
                exit(1);
        }
        range = whiteValue - blackValue;
        if (range == 0.0)
                Error("White value cannot be same as black value\n");

        /* load the font (if necessary) */
        if (labelWidth > 0)
        {
                if (!ReadImage(fontFileName, &font, &fontWidth, &fontHeight, -1, -1, -1, -1, msg))
                {
                        // if we can't find the default font file, try looking in the
                        //   current directory for font.pgm; this is a hack
                        //   to allow align to run before being formally installed
                        if (strcmp(fontFileName, EXPAND_AND_QUOTE(FONT_FILE)) != 0 ||
                            !ReadImage("font.pgm", &font, &fontWidth, &fontHeight,
                                       -1, -1, -1, -1, msg2))
                                Error("Could not load font: %s\n", msg);
                }
                fontHeight = fontHeight / 95;
        }

        if (imageName[0] != '\0')
        {
                images = (Image *) malloc(sizeof(Image));
                images[0].name = (char *) malloc(strlen(imageName) + 1);
                strcpy(images[0].name, imageName);
                images[0].width = -1;
                images[0].height = -1;
                images[0].image = NULL;
                images[0].mask = NULL;
                images[0].dist = NULL;
                images[0].map = NULL;
                images[0].invMap = NULL;
                images[0].imap = NULL;
                images[0].targetMap = NULL;
                nImages = 1;
                imagesSize = 1;
        }
        else
        {
                /* read the image list */
                nImages = 0;
                imagesSize = 0;
                f = fopen(imageListName, "r");
                if (f == NULL)
                        Error("Could not open file %s for reading\n", imageListName);
                while (fgets(line, LINE_LENGTH, f) != NULL)
                {
                        if (line[0] == '\0' || line[0] == '#')
                                continue;
                        width = -1;
                        height = -1;
                        nItems = sscanf(line, "%s%d%d", imageName, &width, &height);
                        if (nItems != 1 && nItems != 3)
                                Error("Malformed line in %s:\n%s\n", imageListName, line);

                        if (nImages >= imagesSize)
                        {
                                imagesSize = (imagesSize > 0) ? imagesSize * 2 : 64;
                                images = (Image*) realloc(images, imagesSize * sizeof(Image));
                        }
                        images[nImages].name = (char*) malloc(strlen(imageName) + 1);
                        strcpy(images[nImages].name, imageName);
                        images[nImages].width = width;
                        images[nImages].height = height;
                        images[nImages].image = NULL;
                        images[nImages].mask = NULL;
                        images[nImages].dist = NULL;
                        images[nImages].map = NULL;
                        images[nImages].invMap = NULL;
                        images[nImages].imap = NULL;
                        images[nImages].targetMap = NULL;
                        ++nImages;
                }
                fclose(f);
                images = (Image *) realloc(images, nImages * sizeof(Image));
                printf("nImages = %d\n", nImages);
        }
        for (i = 0; i < nImages; ++i)
        {
                /* check that the image exists and get its modification time */
                /* FIX:  TEMPORARY HACK */
                sprintf(fn, "%s%s.%s", imagesName, images[i].name, extension);
                if (stat(fn, &sb) != 0)
                {
                        sprintf(fn, "%s%s.tif", imagesName, images[i].name);
                        if (stat(fn, &sb) != 0)
                                Error("Could not stat file %s\n", fn);
                }
                if (S_ISDIR(sb.st_mode))
                        Error("Image %s is a directory.\n", fn);
                images[i].mtime = sb.st_mtime;

                /* read the image size */
                if (images[i].width < 0 || images[i].height < 0)
                {
                        if (!ReadImageSize(fn, &width, &height, msg))
                                Error("Could not determine image size of %s:\n%s\n", fn, msg);
                        images[i].width = width;
                        images[i].height = height;
                }
        }

        printf("Previewing maps: ");
        fflush(stdout);
        oMinX = 1000000000;
        oMaxX = -1000000000;
        oMinY = 1000000000;
        oMaxY = -1000000000;
        map = NULL;
        imap = NULL;
        sourceMapFactor = 1 << sourceMapLevel;
        sourceMapMask = sourceMapFactor - 1;
        targetMapsFactor = 1 << targetMapsLevel;
        for (i = 0; i < nImages; ++i)
        {
                sprintf(fn, "%s%s.map", mapsName, images[i].name);
                if (stat(fn, &sb) != 0)
                        Error("Could not stat map %s\n", fn);
                if (sb.st_mtime > images[i].mtime)
                        images[i].mtime = sb.st_mtime;
                if (!ReadMap(fn, &map, &mLevel,
                             &mw, &mh, &mxMin, &myMin,
                             imName0, imName1,
                             msg))
                        Error("Could not read map %s:\n  error: %s\n",
                              fn, msg);

                spacing = (1 << mLevel) * mapScale;
                minX = 1000000000;
                maxX = -1000000000;
                minY = 1000000000;
                maxY = -1000000000;
                for (y = 0; y < mh-1; ++y)
                        for (x = 0; x < mw-1; ++x)
                        {
                                if (map[y*mw+x].c == 0.0 ||
                                    map[y*mw+x+1].c == 0.0 ||
                                    map[(y+1)*mw+x].c == 0.0 ||
                                    map[(y+1)*mw+x+1].c == 0.0)
                                        continue;
                                for (dy = 0; dy < 2; ++dy)
                                        for (dx = 0; dx < 2; ++dx)
                                        {
                                                rx = map[(y+dy)*mw+x+dx].x * spacing;
                                                ry = map[(y+dy)*mw+x+dx].y * spacing;
                                                if (rotation != 0.0)
                                                {
                                                        rxp = rx - rotationX;
                                                        ryp = ry - rotationY;
                                                        rx = cosRot * rxp + sinRot * ryp + rotationX;
                                                        ry = -sinRot * rxp + cosRot * ryp + rotationY;
                                                }
                                                if (rx < minX)
                                                        minX = rx;
                                                if (rx > maxX)
                                                        maxX = rx;
                                                if (ry < minY)
                                                        minY = ry;
                                                if (ry > maxY)
                                                        maxY = ry;
                                        }
                        }
                free(map);
                images[i].mapBytes = mw * mh * (sizeof(MapElement) + sizeof(InverseMapElement) +
                                                sizeof(unsigned char)) +
                                     5 * 5 * sizeof(long long) +
                                     sizeof(InverseMap) +
                                     256;

                printf("IMAGE BOUNDS[%d] = %f %f %f %f\n", i, minX, maxX, minY, maxY);
                images[i].minX = minX;
                images[i].maxX = maxX;
                images[i].minY = minY;
                images[i].maxY = maxY;

                if (minX < oMinX)
                        oMinX = (int) floor(minX);
                if (maxX > oMaxX)
                        oMaxX = (int) ceil(maxX);
                if (minY < oMinY)
                        oMinY = (int) floor(minY);
                if (maxY > oMaxY)
                        oMaxY = (int) ceil(maxY);

                if (imapsName[0] != '\0')
                {
                        sprintf(fn, "%s%s.map", imapsName, images[i].name);
                        if (!ReadMap(fn, &imap, &imapLevel,
                                     &imapw, &imaph, &imapXMin, &imapYMin,
                                     imapName0, imapName1,
                                     msg))
                                Error("Could not read map %s:\n  error: %s\n",
                                      fn, msg);
                        images[i].mapBytes += imapw * imaph * sizeof(MapElement);
                        free(imap);
                }

                if (targetMapsName[0] != '\0')
                        images[i].mapBytes +=
                                ((images[i].width + targetMapsFactor - 1) / targetMapsFactor) *
                                ((images[i].height + targetMapsFactor - 1) / targetMapsFactor) *
                                sizeof(MapElement);

                if ((nProcessed % 50) == 0 && nProcessed != 0)
                        printf(" %d\n    ", nProcessed);
                printf(".");
                fflush(stdout);
                ++nProcessed;
        }

        printf("\nAll images previewed.\n");

        std::ofstream file_obj;
        file_obj.open("test_saved_image_metadata.dat",std::ios::app);
        file_obj.write((char*)&nImages, sizeof(nImages));
        file_obj.write((char*)images, nImages*sizeof(Image));

        printf("done.\n");
        fflush(stdout);
        return(0);
}

void
PaintImage (int i, int minX, int maxX, int minY, int maxY)
{
        int j;
        int iw, ih;
        int x, y;
        InverseMap *invMap;
        float xv, yv;
        int ixv, iyv;
        float rx, ry;
        float rrx, rry;
        float rv, dv;
        size_t maskBytes;
        size_t mbpl;
        unsigned char *mask;
        unsigned char *dist;
        float *distance;
        float dst;
        int idst;
        float cx, cy;
        float offset;
        float r00, r01, r10, r11;
        float d00, d01, d10, d11;
        int nx, ny;
        int pMinX, pMaxX, pMinY, pMaxY;
        unsigned char *image;
        float w;
        float d, d2;
        int row, col;
        int ix, iy;
        float mFactor;
        int offsetX, offsetY;
        int v;
        size_t newCanvasSize;
        int mxMin, myMin;
        char imName0[PATH_MAX];
        char imName1[PATH_MAX];
        char msg[PATH_MAX+256];
        char fn[PATH_MAX];
        struct stat statBuf;
        int complete;
        int iixv, iiyv;
        float rb00, rb01, rb10, rb11;
        float rw00, rw01, rw10, rw11;
        float rb, rw;
        MapElement *imap;
        float imapFactor;
        int imapw, imaph;
        int imapXMin, imapYMin;
        char imapName0[PATH_MAX], imapName1[PATH_MAX];
        float xvi, yvi;
        int sourceMapMask;
        int smi;
        int targetMapSize;
        MapElement *targetMap;
        int targetMapWidth, targetMapHeight;
        int targetMapMask;
        int tmi;
        int testX, testY;
        int mw, mh;
        MapElement *map;
        float spacing;
        float rxp, ryp;
        unsigned char *imask;
        size_t imbpl;
        int warned;

        /* read in map if necessary */
        if (images[i].map == NULL)
        {
                sprintf(fn, "%s%s.map", mapsName, images[i].name);
                if (!ReadMap(fn, &(images[i].map), &(images[i].mLevel),
                             &(images[i].mw), &(images[i].mh),
                             &(images[i].mxMin), &(images[i].myMin),
                             imName0, imName1,
                             msg))
                        Error("Could not read map %s:\n  error: %s\n",
                              fn, msg);
                if (rotation != 0.0)
                {
                        map = images[i].map;
                        spacing = (1 << images[i].mLevel) * mapScale;
                        mw = images[i].mw;
                        mh = images[i].mh;
                        for (y = 0; y < mh; ++y)
                                for (x = 0; x < mw; ++x)
                                {
                                        rx = map[y*mw+x].x * spacing;
                                        ry = map[y*mw+x].y * spacing;
                                        rxp = rx - rotationX;
                                        ryp = ry - rotationY;
                                        rx = cosRot * rxp + sinRot * ryp + rotationX;
                                        ry = -sinRot * rxp + cosRot * ryp + rotationY;
                                        map[y*mw+x].x = rx / spacing;
                                        map[y*mw+x].y = ry / spacing;
                                }
                }

                imageMem += images[i].mapBytes;
        }
        if (images[i].invMap == NULL)
                images[i].invMap = InvertMap(images[i].map,
                                             images[i].mw,
                                             images[i].mh);

        /* read in image if necessary */
        if (images[i].image == NULL)
        {
                sprintf(fn, "%s%s", imagesName, images[i].name);
                if (!ReadImage(fn, &(images[i].image),
                               &iw, &ih,
                               -1, -1, -1, -1,
                               msg))
                        Error("Could not read image %s:\n  error: %s\n",
                              fn, msg);
                if (iw != images[i].width ||
                    ih != images[i].height)
                        Error("Dimensions of image %s do not match those in images list.\n",
                              fn);
                imageMem += images[i].width * images[i].height;

                if (targetMapsName[0] != '\0')
                {
                        targetMapWidth = (iw + targetMapsFactor - 1) / targetMapsFactor;
                        targetMapHeight = (ih + targetMapsFactor - 1) / targetMapsFactor;
                        targetMapSize = targetMapWidth * targetMapHeight;
                        //	  printf("malloc %zu bytes for targetMap\n",
                        //		 targetMapSize * sizeof(MapElement));
                        images[i].targetMap =
                                (MapElement *) malloc(targetMapSize * sizeof(MapElement));
                        if (images[i].targetMap == 0)
                                Error("malloc of targetMap failed; errno = %d\n", errno);
                        memset(images[i].targetMap, 0,
                               targetMapSize * sizeof(MapElement));
                }
                else
                        images[i].targetMap = NULL;
        }
        targetMapWidth = (images[i].width + targetMapsFactor - 1) / targetMapsFactor;
        targetMapHeight = (images[i].height + targetMapsFactor - 1) / targetMapsFactor;
        targetMapSize = targetMapWidth * targetMapHeight;
        targetMap = images[i].targetMap;
        targetMapMask = targetMapsFactor - 1;

        /* read in mask if necessary */
        if (images[i].mask == NULL)
                if (masksName[0] != '\0')
                {
                        sprintf(fn, "%s%s", masksName, images[i].name);
                        if (!ReadBitmap(fn, &(images[i].mask),
                                        &iw, &ih,
                                        -1, -1, -1, -1,
                                        msg))
                                Error("Could not read mask %s:\n  error: %s\n",
                                      fn, msg);
                        if (maskScale == 1.0 &&
                            (iw != images[i].width ||
                             ih != images[i].height) ||
                            maskScale < 1.0 &&
                                        ((int) floor(maskScale * iw + 0.5) != images[i].width ||
                                         (int) floor(maskScale * ih + 0.5) != images[i].height) ||
                                        maskScale > 1.0 &&
                            ((int) floor(images[i].width / maskScale) != iw ||
                             (int) floor(images[i].height / maskScale) != ih))
                                Error("Dimensions of mask %s do not match those in images list (%f %d %d %d %d)\n",
                                      fn, maskScale, iw, ih, images[i].width, images[i].height);
                        mbpl = (images[i].width + 7) / 8;
                        maskBytes = mbpl * images[i].height;
                        if (maskScale != 1.0)
                        {
                                mask = (unsigned char *) malloc(maskBytes);
                                memset(mask, 0, maskBytes);
                                imask = images[i].mask;
                                imbpl = (iw + 7) / 8;
                                for (y = 0; y < images[i].height; ++y)
                                {
                                        iy = (int) floor(y / maskScale + 0.001);
                                        for (x = 0; x < images[i].width; ++x)
                                        {
                                                ix = (int) floor(x / maskScale + 0.001);
                                                if (imask[iy*imbpl + (ix >> 3)] & (0x80 >> (ix & 7)))
                                                        mask[y*mbpl + (x >> 3)] |= 0x80 >> (x & 7);
                                        }
                                }
                                free(images[i].mask);
                                images[i].mask = mask;
                        }

                        mask = images[i].mask;
                        for (j = 0; j < maskBytes; ++j)
                                mask[j] ^= 0xff;
                        imageMem += maskBytes;
                }
                else
                {
                        maskBytes = ((images[i].width + 7) / 8) * images[i].height;
                        //	printf("malloc %zu bytes for mask\n", maskBytes);
                        mask = images[i].mask = (unsigned char *) malloc(maskBytes);
                        if (mask == 0)
                                Error("malloc of mask failed; errno = %d\n", errno);
                        memset(images[i].mask, 0, maskBytes);
                        imageMem += maskBytes;
                }

        /* read in intensity map if necessary */
        if (images[i].imap == NULL)
                if (imapsName[0] != '\0')
                {
                        sprintf(fn, "%s%s.map", imapsName, images[i].name);
                        if (!ReadMap(fn, &(images[i].imap), &(images[i].imapLevel),
                                     &(images[i].imapw), &(images[i].imaph),
                                     &imapXMin, &imapYMin,
                                     imapName0, imapName1,
                                     msg))
                                Error("Could not read map %s:\n  error: %s\n",
                                      fn, msg);
                        if (imapXMin != 0 || imapYMin != 0)
                                Error("Can not handle partial intensity map: %s\n", fn);
                        imap = images[i].imap;
                        imageMem += images[i].imapw * images[i].imaph * sizeof(MapElement);
                }

        /* compute distance table if necessary */
        if (images[i].dist == NULL)
        {
                iw = images[i].width;
                ih = images[i].height;
                //      printf("malloc %zu bytes for dist\n", iw * ih * sizeof(unsigned char));
                images[i].dist = (unsigned char*) malloc(iw * ih * sizeof(unsigned char));
                if (images[i].dist == 0)
                        Error("malloc of images[i].dist failed; errno = %d\n", errno);
                distance = (float*) malloc(iw * ih * sizeof(float));
                if (distance == 0)
                        Error("malloc of distance failed; errno = %d\n", errno);
                computeDistance(EUCLIDEAN_DISTANCE,
                                iw, ih, images[i].mask,
                                distance);

                /* if distance from edge is less, use that */
                dist = images[i].dist;
                for (y = 0; y < ih; ++y)
                        for (x = 0; x < iw; ++x)
                        {
                                dst = distance[y*iw+x];
                                d = x + 1;
                                if (d < dst)
                                        dst = d;
                                d = iw - x;
                                if (d < dst)
                                        dst = d;
                                d = y + 1;
                                if (d < dst)
                                        dst = d;
                                d = ih - y;
                                if (d < dst)
                                        dst = d;
                                idst = (int) floor(64.0 * dst);
                                if (idst > 255)
                                        idst = 255;
                                else if (idst < 0)
                                        idst = 0;
                                dist[y*iw + x] = idst;
                        }
                free(distance);
                imageMem += iw * ih;
        }

        /* paint the image on the canvas */
        pMinX = minX;
        pMaxX = maxX;
        pMinY = minY;
        pMaxY = maxY;
        if (images[i].minX > pMinX)
                pMinX = images[i].minX;
        if (images[i].maxX < pMaxX)
                pMaxX = images[i].maxX;
        if (images[i].minY > pMinY)
                pMinY = images[i].minY;
        if (images[i].maxY < pMaxY)
                pMaxY = images[i].maxY;

        mFactor = (1 << images[i].mLevel) * mapScale;
        invMap = images[i].invMap;
        image = images[i].image;
        mask = images[i].mask;
        imap = images[i].imap;
        if (imap != NULL)
        {
                imapFactor = (1 << images[i].imapLevel) * imapScale;
                imapw = images[i].imapw;
                imaph = images[i].imaph;
        }
        dist = images[i].dist;
        iw = images[i].width;
        ih = images[i].height;
        mxMin = images[i].mxMin;
        myMin = images[i].myMin;
        //  printf("mxMin = %d myMin = %d\n", mxMin, myMin);
        cx = (iw - 1) / 2.0;
        cy = (ih - 1) / 2.0;
        mbpl = (iw + 7) / 8;
        offset = -1000000.0;
        printf("Rendering %s from x=%d to %d y=%d to %d\n",
               images[i].name, pMinX, pMaxX, pMinY, pMaxY);
        testX = (pMinX + pMaxX) / 2;
        testY = (pMinY + pMaxY) / 2;
        warned = 0;
        for (y = pMinY; y <= pMaxY; ++y)
                for (x = pMinX; x <= pMaxX; ++x)
                {
                        //	if (y == testY && x == testX)
                        //	  printf("TEST STARTED\n");
                        if (!Invert(invMap, &xv, &yv, (x + 0.5) / mFactor, (y + 0.5) / mFactor))
                        {
                                //	    if (y == testY && x == testX)
                                //	      printf("TEST INVERT FAILED %f %f %f %f\n", (x + 0.5)/mFactor,
                                //		     (y+0.5)/mFactor, xv, yv);
                                continue;
                        }
                        xv = (xv + mxMin) * mFactor - 0.5;
                        yv = (yv + myMin) * mFactor - 0.5;
                        //	if (y == testY && x == testX)
                        //	  printf("TEST INVERT OK %f %f %f %f\n", (x + 0.5)/mFactor,
                        //		 (y+0.5)/mFactor, xv, yv);
                        ixv = (int) floor(xv);
                        iyv = (int) floor(yv);
                        if (ixv < -1 || ixv >= iw ||
                            iyv < -1 || iyv >= ih)
                                continue;
                        rrx = xv - ixv;
                        rry = yv - iyv;
                        if (ixv >= 0 && iyv >= 0 &&
                            !(mask[iyv * mbpl + (ixv >> 3)] & (0x80 >> (ixv & 7))))
                        {
                                r00 = image[iyv * iw + ixv];
                                d00 = dist[iyv * iw + ixv];
                                if (targetMap != NULL && ((ixv | iyv) & targetMapMask) == 0)
                                {
                                        tmi = (iyv >> targetMapsLevel) * targetMapWidth +
                                              (ixv >> targetMapsLevel);
                                        if (tmi >= targetMapSize)
                                                Error("tmi out-of-bounds: %d %d\n%d %d\n",
                                                      tmi, targetMapSize,
                                                      targetMapWidth, targetMapHeight);
                                        targetMap[tmi].x = (float) x;
                                        targetMap[tmi].y = (float) y;
                                        targetMap[tmi].c = 1.0;
                                }
                        }
                        else
                        {
                                r00 = 0.0;
                                d00 = 0.0;
                        }
                        if (ixv >= 0 && iyv+1 < ih &&
                            !(mask[(iyv + 1) * mbpl + (ixv >> 3)] & (0x80 >> (ixv & 7))))
                        {
                                r01 = image[(iyv + 1) * iw + ixv];
                                d01 = dist[(iyv + 1) * iw + ixv];
                                if (targetMap != NULL && ((ixv | (iyv+1)) & targetMapMask) == 0)
                                {
                                        tmi = ((iyv + 1) >> targetMapsLevel) * targetMapWidth +
                                              (ixv >> targetMapsLevel);
                                        if (tmi >= targetMapSize)
                                                Error("tmi out-of-bounds: %d %d\n%d %d\n",
                                                      tmi, targetMapSize,
                                                      targetMapWidth, targetMapHeight);
                                        targetMap[tmi].x = (float) x;
                                        targetMap[tmi].y = (float) y;
                                        targetMap[tmi].c = 1.0;
                                }
                        }
                        else
                        {
                                r01 = 0.0;
                                d01 = 0.0;
                        }
                        if (ixv+1 < iw && iyv >= 0 &&
                            !(mask[iyv * mbpl + ((ixv+1) >> 3)] & (0x80 >> ((ixv+1) & 7))))
                        {
                                r10 = image[iyv * iw + ixv + 1];
                                d10 = dist[iyv * iw + ixv + 1];
                                if (targetMap != NULL && (((ixv + 1) | iyv) & targetMapMask) == 0)
                                {
                                        tmi = (iyv >> targetMapsLevel) * targetMapWidth +
                                              ((ixv + 1) >> targetMapsLevel);
                                        if (tmi >= targetMapSize)
                                                Error("tmi out-of-bounds: %d %d\n%d %d\n",
                                                      tmi, targetMapSize,
                                                      targetMapWidth, targetMapHeight);
                                        targetMap[tmi].x = (float) x;
                                        targetMap[tmi].y = (float) y;
                                        targetMap[tmi].c = 1.0;
                                }
                        }
                        else
                        {
                                r10 = 0.0;
                                d10 = 0.0;
                        }
                        if (ixv+1 < iw && iyv+1 < ih &&
                            !(mask[(iyv+1) * mbpl + ((ixv+1) >> 3)] & (0x80 >> ((ixv+1) & 7))))
                        {
                                r11 = image[(iyv + 1) * iw + ixv + 1];
                                d11 = dist[(iyv + 1) * iw + ixv + 1];
                                if (targetMap != NULL && (((ixv + 1) | (iyv+1)) & targetMapMask) == 0)
                                {
                                        tmi = ((iyv + 1) >> targetMapsLevel) * targetMapWidth +
                                              ((ixv + 1) >> targetMapsLevel);
                                        if (tmi >= targetMapSize)
                                                Error("tmi out-of-bounds: %d %d\n%d %d\n",
                                                      tmi, targetMapSize,
                                                      targetMapWidth, targetMapHeight);
                                        targetMap[tmi].x = (float) x;
                                        targetMap[tmi].y = (float) y;
                                        targetMap[tmi].c = 1.0;
                                }
                        }
                        else
                        {
                                r11 = 0.0;
                                d11 = 0.0;
                        }
                        rv = r00 * (rrx - 1.0) * (rry - 1.0)
                             - r10 * rrx * (rry - 1.0)
                             - r01 * (rrx - 1.0) * rry
                             + r11 * rrx * rry;
                        dv = d00 * (rrx - 1.0) * (rry - 1.0)
                             - d10 * rrx * (rry - 1.0)
                             - d01 * (rrx - 1.0) * rry
                             + d11 * rrx * rry;
                        if (dv <= 0.0)
                                continue;

                        if (imap != NULL)
                        {
                                /* lookup xv,yv in intensity map, if present */
                                xvi = (xv + 0.5) / imapFactor;
                                yvi = (yv + 0.5) / imapFactor;
                                iixv = (int) floor(xvi);
                                iiyv = (int) floor(yvi);
                                rrx = xvi - iixv;
                                rry = yvi - iiyv;
                                while (iixv < 0)
                                {
                                        ++iixv;
                                        rrx -= 1.0;
                                }
                                while (iixv >= imapw-1)
                                {
                                        --iixv;
                                        rrx += 1.0;
                                }
                                while (iiyv < 0)
                                {
                                        ++iiyv;
                                        rry -= 1.0;
                                }
                                while (iiyv >= imaph-1)
                                {
                                        --iiyv;
                                        rry += 1.0;
                                }
                                rb00 = imap[iiyv*imapw+iixv].x;
                                rw00 = imap[iiyv*imapw+iixv].y;
                                rb01 = imap[(iiyv+1)*imapw+iixv].x;
                                rw01 = imap[(iiyv+1)*imapw+iixv].y;
                                rb10 = imap[iiyv*imapw+iixv+1].x;
                                rw10 = imap[iiyv*imapw+iixv+1].y;
                                rb11 = imap[(iiyv+1)*imapw+iixv+1].x;
                                rw11 = imap[(iiyv+1)*imapw+iixv+1].y;
                                rb = rb00 * (rrx - 1.0) * (rry - 1.0)
                                     - rb10 * rrx * (rry - 1.0)
                                     - rb01 * (rrx - 1.0) * rry
                                     + rb11 * rrx * rry;
                                rw = rw00 * (rrx - 1.0) * (rry - 1.0)
                                     - rw10 * rrx * (rry - 1.0)
                                     - rw01 * (rrx - 1.0) * rry
                                     + rw11 * rrx * rry;
                                if (rw <= rb)
                                {
                                        if (!warned)
                                        {
                                                printf("Warning: rw level (%f) is less than rb level (%f) for image %s at (%d %d)\n",
                                                       rw, rb, images[i].name, ixv, iyv);
                                                warned = 1;
                                        }
                                        rw = 255.0;
                                        rb = 0.0;
                                }
                                rv = (rv - rb) / (rw - rb);
                        }

                        if (dv < 64.0)
                                w = 65535 - (int) floor(4.0 * dv);
                        else
                                w = (int) floor(hypotf(xv - cx, yv - cy));
                        if (w < weight[(y - minY) * weightWidth + x - weightMinX])
                        {
                                weight[(y - minY) * weightWidth + x - weightMinX] = w;
                                if (imap != NULL)
                                        v = (int) floor(255.0 * rv + 0.5);
                                else
                                        v = (int) floor(255.0 * (rv - blackValue) / range + 0.5);
                                if (v <= 0)
                                {
                                        if (w == 65535)
                                                v = 0;
                                        else
                                                v = 1;
                                }
                                else if (v > 255)
                                        v = 255;
                                canvas[(y - minY) * canvasWidth + x - canvasMinX] = v;

                                if (sourceMap != NULL &&
                                    (((x - oMinX) | (y - oMinY)) & sourceMapMask) == 0)
                                {
                                        smi = ((y - oMinY) >> sourceMapLevel) * sourceMapWidth +
                                              ((x - oMinX) >> sourceMapLevel);
                                        if (smi >= sourceMapSize)
                                                Error("smi out-of-bounds\n");
                                        sourceMap[smi].x = (float) ixv;
                                        sourceMap[smi].y = (float) iyv;
                                        sourceMap[smi].c = (float) i;
                                }
                        }
                }

        /* free up if no longer required */
        if (images[i].maxX <= maxX || maxX >= oMaxX)
        {
                if (images[i].targetMap != NULL)
                {
                        sprintf(fn, "%s%s.map", targetMapsName, images[i].name);
                        if (!CreateDirectories(fn))
                                Error("Could not create directories for target map file %s\n",
                                      fn);
                        if (!WriteMap(fn, images[i].targetMap, targetMapsLevel,
                                      targetMapWidth, targetMapHeight,
                                      0, 0,
                                      images[i].name, outputName,
                                      UncompressedMap, msg))
                                Error("Could not write target map %s:\n%s\n", fn, msg);
                        //	  printf("freeing %zu bytes from targetMap\n",
                        //		 targetMapHeight * targetMapWidth * sizeof(MapElement));
                        free(images[i].targetMap);
                        images[i].targetMap = NULL;
                }
                if (images[i].invMap != NULL)
                {
                        FreeInverseMap(images[i].invMap);
                        images[i].invMap = NULL;
                }
                if (images[i].map != NULL)
                {
                        free(images[i].map);
                        images[i].map = NULL;
                        imageMem -= images[i].mapBytes; // this accounts for the inverse
                                                        // map and target map as well
                }
                if (images[i].imap != NULL)
                {
                        free(images[i].imap);
                        images[i].imap = NULL;
                        imageMem -= images[i].imapw * images[i].imaph * sizeof(MapElement);
                }
                if (images[i].image != NULL)
                {
                        //	  printf("freeing %zu bytes from image\n",
                        //		 images[i].height * images[i].width * sizeof(unsigned char));
                        free(images[i].image);
                        images[i].image = NULL;
                        imageMem -= images[i].width * images[i].height;
                }
                if (images[i].mask != NULL)
                {
                        //	  printf("freeing %zu bytes from mask\n",
                        //		 (size_t) ((images[i].width + 7) / 8) * images[i].height);
                        free(images[i].mask);
                        images[i].mask = NULL;
                        imageMem -= ((images[i].width + 7) / 8) * images[i].height;
                }
                if (images[i].dist != NULL)
                {
                        //	  printf("freeing %zu bytes from dist\n",
                        //		 (size_t) (images[i].width  * images[i].height));
                        free(images[i].dist);
                        images[i].dist = NULL;
                        imageMem -= images[i].width * images[i].height;
                }
        }
}

void
WriteTiles (int col, int startRow, int endRow, char *iName)
{
        int row;
        char fn[PATH_MAX];
        char msg[PATH_MAX+256];

        for (row = startRow; row <= endRow; ++row)
        {
                if (overlay)
                {
                        if (tileWidth < 0 && tileHeight < 0)
                                sprintf(fn, "%s.tif", outputName);
                        else
                                sprintf(fn, "%sc%.2d%sr%.2d.tif", outputName, col+1,
                                        tree ? "/" : "", row+1);
                }
                else
                {
                        if (tileWidth < 0 && tileHeight < 0)
                                sprintf(fn, "%s%s.tif",
                                        outputName, iName);
                        else
                                sprintf(fn, "%s%s/c%.2d%sr%.2d.tif",
                                        outputName, iName,
                                        col+1, tree ? "/" : "", row+1);
                }
                if (!CreateDirectories(fn))
                        Error("Could not create directories for output file %s\n", fn);
                if (!WriteImage(fn, &out[(row - startRow) * th * tw], (int) tw, (int) th,
                                compress ? HDiffDeflateImage : UncompressedImage,
                                msg))
                        Error("Could not write output file %s:\n  error: %s\n",
                              fn, msg);
                if ((nProcessed % 50) == 0 && nProcessed != 0)
                        printf(" %d\n    ", nProcessed);
                printf(".");
                fflush(stdout);
                ++nProcessed;
        }
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

unsigned int
Hash (char *s)
{
        unsigned int v = 0;
        char *p = s;
        while (*p != '\0')
                v = 37 * v + *p++;
        return(v);
}

int
CreateDirectories (char *fn)
{
        int pos;
        char dn[PATH_MAX];
        char *slash;
        int len;
        struct stat statBuf;
        int hv;
        int i;

        for (pos = 0;
             (slash = strchr(&fn[pos], '/')) != NULL;
             pos = len + 1)
        {
                len = slash - fn;
                if (len == 0)
                        continue;
                strncpy(dn, fn, len);
                dn[len] = '\0';

                /* check if in hash table */
                hv = 0;
                for (i = 0; i < len; ++i)
                        hv = 239*hv + dn[i];
                hv &= DIR_HASH_SIZE-1;
                i = hv;
                while (dirHash[i] != NULL)
                {
                        if (strcmp(dirHash[i], dn) == 0)
                                break;
                        i = (i + 1) & (DIR_HASH_SIZE-1);
                        if (i == hv)
                                Error("Directory hash table is full!\n");
                }
                if (dirHash[i] != NULL)
                        continue;
                dirHash[i] = (char *) malloc(strlen(dn)+1);
                if (dirHash[i] == 0)
                        Error("malloc of dirHash[i] failed; errno = %d\n", errno);
                strcpy(dirHash[i], dn);

                if (stat(dn, &statBuf) == 0)
                {
                        if (S_ISDIR(statBuf.st_mode) ||
                            S_ISLNK(statBuf.st_mode))
                                continue;
                        Error("Output path component %s is not a directory\n", dn);
                }
                if (errno != ENOENT)
                        Error("Could not stat directory %s\n", dn);

                if (mkdir(dn, 0775) != 0)
                        Error("Could not create directory %s\n", dn);
        }
        return(1);
}

void
PrintUsage ()
{
        struct rusage usage;

        printf("imageMem = %zu\n", imageMem);
        getrusage(RUSAGE_SELF, &usage);
        printf("getrusage_maxrss = %ld\n", usage.ru_maxrss);
}
