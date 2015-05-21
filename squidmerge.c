//
// Program to combine seperate tile images that are from a common squid.
//
// -------------------------- LICENSE -----------------------------------
//
// This file is part of the LibSQUID software libraray.
//
// LibSQUID is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// LibSQUID is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with LibSQUID.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright 2014 James Wren and Los Alamos National Laboratory
//

#define _GNU_SOURCE 

#include <errno.h>
#include <fcntl.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#include <libsquid.h>
#include <libsquid_wcs.h>

#define SQUIDMERGE_VERSION 0.6.0
#define SQUIDMERGE_MAJOR 0
#define SQUIDMERGE_MINOR 6
#define SQUIDMERGE_PATCH 0
#define SQUIDMERGE_RELEASE ""

void printUsage(char* progName) {
   printf("Calling sequence: %s [OPTION]... infile0 ... infileN outfile\n",progName);
   printf("Mandatory arguments to long options are mandatory for short options too.\n");
   printf("Options:\n");
   printf("\t-c, --compressed\n");
   printf("\t-i, --initPixels value\n");
   printf("\t-m, --mode type\n");
   printf("\t\tMode types:\n");
   printf("\t\t\tmagicNull value tolerance\n");
   printf("\t\t\tthreshold value\n");
   printf("\t-v, --verbose\n");
}

static int isCompressed;
static int isVerbose;

int main(int argc, char *argv[]) {
   typedef enum {
      PixelValidationMode_MagicNull,
      PixelValidationMode_Threshold
   } EnumPixelValidationMode;
   EnumPixelValidationMode pixelValidationMode;
   double pixelInitValue;
   double pixelMagicNullValue;
   double pixelMagicNullTolerance;
   double pixelThresholdValue;
   int optRetval, optionIndex; // for getopt call
   const char* shortOptions = "ci:m:v";
   static struct option longOptions[] = {
      {"compressed", no_argument, &isCompressed, 1},
      {"initPixels", required_argument, 0, 'i'},
      {"mode",       required_argument, 0, 'm'},
      {"verbose",    no_argument, &isVerbose, 1},
      {0,0,0,0}
   };
   int i;
   fitsfile *bptr, *fptr; // base and input pointers
   int* contribArr;  //used to track contributing pixels
   double *bimg, *fimg; // base and input images
   int status=0; // status of cfitsio function calls
   int bitpix; // image parameters
   long bnaxes[2], fnaxes[2]; // image dimensions
   int bnaxis, fnaxis; // number of image axes
   long bnpix, fnpix; // number of pixels in img
   long bmid, fmid; // map id
   long pix,ohit; // pixel index, hits
   double contribution; // contribution in percent
   double coverage; // coverage in percent
   char tmpname[FNAME_MAX],outname[FNAME_MAX];
   char compressedStr[FNAME_MAX];
   fitsfile *ofptr; // output file pointer

   //set default operation
   pixelValidationMode = PixelValidationMode_Threshold;
   isCompressed        = 0;
   isVerbose           = 0;
   pixelInitValue      = 0.0;
   pixelThresholdValue = 10.0;

   //parse command line options
   while (1) {
      optionIndex = 0;
      optRetval = getopt_long(argc, argv, shortOptions, longOptions, &optionIndex);
      if (-1 == optRetval) {
         break;
      }

      switch (optRetval) {
         case 'c':
            //compressed
            isCompressed = 1;
            break;
         case 'i':
            //initPixels
            pixelInitValue = atof(optarg);
            break;
         case 'm':
            //mode
            if (0 == strcmp("magicNull", optarg)) {
               if (argc < (optind+2)) {
                  printf("Insufficient arguments for mode %s\n", optarg);
                  return -1;
               }
               pixelValidationMode     = PixelValidationMode_MagicNull;
               pixelMagicNullValue     = atof(argv[optind]);
               pixelMagicNullTolerance = atof(argv[optind+1]);
               optind += 2;
            }
            else if (0 == strcmp("threshold", optarg)) {
               if (argc < (optind+1)) {
                  printf("Insufficient arguments for mode %s\n", optarg);
                  return -1;
               }
               pixelValidationMode = PixelValidationMode_Threshold;
               pixelThresholdValue = atof(argv[optind]);
               optind += 1;
            }
            else {
               printf("Unhandled argument for option mode:  %s\n", optarg);
               return -1;
            }
            break;
         case 'v':
            //verbose
            isVerbose = 1;
            break;
         case '?':
            //error already reported by getopt_long
            break;
         case 0:
            if (0 == strcmp("compressed", longOptions[optionIndex].name)) {
               isCompressed = 1;
            }
            else if (0 == strcmp("verbose", longOptions[optionIndex].name)) {
               isVerbose = 1;
            }
            else {
               printf("Unhandled long option:\n");
               printf("\tOption %s\n", longOptions[optionIndex].name);
               printf("\tArgument %s\n", optarg);
            }
            break;
         default:
            printf("Unhandled option:\t%i\t%c\n", optRetval, (char)optRetval);
            printUsage(argv[0]);
            return -1;
      }
   }

   //show operation modes
   if (0 != isVerbose) {
      printf("pixelInitValue:\t\t%e\n", pixelInitValue);
      printf("pixelValidationMode:\t");
      switch (pixelValidationMode) {
         case PixelValidationMode_MagicNull:
            printf("MagicNull\n");
            printf("\tpixelMagicNullValue:\t\t%e\n", pixelMagicNullValue);
            printf("\tpixelMagicNullTolerance:\t%e\n", pixelMagicNullTolerance);
            break;
         case PixelValidationMode_Threshold:
            printf("Threshold\n");
            printf("\tpixelThresholdValue:\t%e\n", pixelThresholdValue);
            break;
         default:
            printf("\nUnhandled pixelValidationMode:\t%i\n", pixelValidationMode);
            printUsage(argv[0]);
            return -1;
            break;
      }
   }

   //check remaining arguments
   if (argc < (optind+2)) {
      printf("Insufficient arguments for input and output fnames\n");
      printUsage(argv[0]);
      return -1;
   }

   //------------------------------
   //Load original image and header
   //------------------------------

   //prepare base image fname
   if (0 == isCompressed) {
      snprintf(tmpname, FNAME_MAX, "%s", argv[optind]); // for uncompressed (see also inside of for loop)
   }
   else {
      snprintf(tmpname, FNAME_MAX, "%s[1]", argv[optind]); // for compressed (see also inside of for loop)
   }

   //open base image
   if (0 != isVerbose) {
      printf("opening base image\n");
      printf("%s\n",tmpname);
   }
   if (fits_open_file(&bptr, tmpname, READONLY, &status)) {
      fprintf(stderr,"fits_open_file failed in %s\n",argv[0]);
      fits_report_error(stderr, status);
      return -1;
   }

   //read base image dimensions
   if (fits_get_img_param(bptr, 2, &bitpix, &bnaxis, bnaxes, &status)) {
      fprintf(stderr,"fits_get_img_param failed in %s\n",argv[0]);
      fits_report_error(stderr, status);
      return -1;
   }
   bnpix = bnaxes[0]*bnaxes[1];
   if (0 != isVerbose) {
      printf("base image naxes are %ld %ld\n", bnaxes[0], bnaxes[1]);
      printf("base image npix is %ld\n", bnpix);
   }

   //read base image MAPID
   if (fits_read_key(bptr, TLONG, "MAPID", &bmid, NULL, &status)) {
      printf("base image keyword MAPID not in header.\n");
      fits_report_error(stderr, status);
      return -1;
   }
   if (0 != isVerbose) {
      printf("base image MAPID is %i\n", bmid);
   }

   //allocate contribution array
   contribArr=malloc(bnpix*sizeof(int));
   if (NULL == contribArr) {
      fprintf(stderr,"malloc contribArr failed in %s, %s\n",argv[0],strerror(errno));
      return -1;
   }
   //initialize each element in contribution array
   for (pix = 0; pix < bnpix; pix++) {
      contribArr[pix]=0;
   }

   //allocate base image array
   bimg=malloc(bnpix*sizeof(double));
   if (NULL == bimg) {
      fprintf(stderr,"malloc bimg failed in %s, %s\n",argv[0],strerror(errno));
      return -1;
   }
   //initialize each pixel in base image
   for (pix = 0; pix < bnpix; pix++) {
      bimg[pix]=pixelInitValue;
   }

   //-----------------------------------------------------
   // loop over image arguments and add them to base image
   //-----------------------------------------------------
   for (i=optind; i<=(argc-2); i++) {
      // Load current image and header

      //prepare current image fname
      if (0 == isCompressed) {
         snprintf(tmpname, FNAME_MAX, "%s", argv[i]); // for uncompressed
      }
      else {
         snprintf(tmpname, FNAME_MAX, "%s[1]", argv[i]); // for compressed
      }

      //open current image
      if (0 != isVerbose) {
         printf("opening current image\n");
      }
      printf("%s\n",argv[i]);
      if (fits_open_file(&fptr, tmpname, READONLY, &status)) {
         fits_report_error(stderr, status);
         return -1;
      }

      //read current image dimensions
      if (fits_get_img_param(fptr, 2, &bitpix, &fnaxis, fnaxes, &status)) {
         fits_report_error(stderr, status);
         return -1;
      }
      if (0 != isVerbose) {
         printf("current image naxes are %ld %ld\n",fnaxes[0],fnaxes[1]);
      }

      //allocate current image array
      fnpix = fnaxes[0]*fnaxes[1];
      if (0 != isVerbose) {
         printf("current image npix is %ld\n",fnpix);
      }
      fimg=malloc(fnpix*sizeof(double));
      if (fimg == NULL) {
         fprintf(stderr,"malloc failed in %s, %s\n",argv[0],strerror(errno));
         return -1;
      }

      //read current image pixels
      if (fits_read_img(fptr, TDOUBLE, 1, fnpix, NULL, fimg, NULL, &status)) {
         fits_report_error(stderr, status);
         return -1;
      }

      //read current image MAPID
      if (fits_read_key(fptr, TLONG, "MAPID", &fmid, NULL, &status)) {
         printf("current image keyword MAPID not in header.\n");
         fits_report_error(stderr, status);
         return -1;
      }

      // check if we are in the same map id and images are same size
      if ((bmid == fmid)&&(bnpix == fnpix)) {
         //keep count of current frame's pixels as used to calculate contribution
         ohit = 0;
         // overwrite base image with valid input image pixels
         for (pix = 0; pix < bnpix; pix++) {
            // only evaluate pixels needing contribution
            if (0 == contribArr[pix]) {
               //test pixel for contribution according to pixel validation mode
               if (PixelValidationMode_Threshold == pixelValidationMode) {
                  //simple thresholding
                  if (pixelThresholdValue > bimg[pix]) {
                     if (pixelThresholdValue < fimg[pix]) {
                        bimg[pix]=fimg[pix];
                        ++contribArr[pix];
                        ++ohit;
                     }
                  }
               }
               else if (PixelValidationMode_MagicNull == pixelValidationMode){
                  //regard magicNull value
                  if (!(pixelMagicNullTolerance > fabs(pixelMagicNullValue - fimg[pix]))) {
                     bimg[pix]=fimg[pix];
                     ++contribArr[pix];
                     ++ohit;
                  }
               }
               else {
                  printf("Unhandled pixel validation mode:  %i\n", pixelValidationMode);
                  return -1;
               }
            }
         }

         contribution=100.0*(ohit/((double)bnpix));
         if (0 != isVerbose) {
            printf("contribution=%.1f\n",contribution);
         }
      } else {
         // skip adding this image
         printf("skipping!\n");
      }

      if (fits_close_file(fptr, &status)) {
         //printf("error freeing original image\n");
         //fits_report_error(stderr, status);
         //if (status != 412){return -1;};
         fits_clear_errmsg();
         status=0;
      }
      free(fimg);
   }

   //count pixels with contribution
   ohit=0;
   for (pix = 0; pix < bnpix; pix++) {
      if (0 != contribArr[pix]) {
         ++ohit;
      }
   }

   //calculate coverage of contributing pixels
   coverage=100.0*(ohit/((double)bnpix));
   if (0 != isVerbose) {
      printf("coverage=%.1f\n",coverage);
   }

   //-------------------
   // Write output image
   //-------------------
   strncpy(outname, argv[argc-1], FNAME_MAX);
   printf("output file is %s\n",outname);
   strncpy(compressedStr,"",FNAME_MAX);

   unlink(outname); //remove file if it already exists

   //prepare output fname
   snprintf(tmpname, FNAME_MAX, "%s%s", outname, compressedStr);
   if (0 != isVerbose) {
      printf("writing %s\n",tmpname);
   }

   //open output image
   if (fits_create_file(&ofptr, tmpname, &status)) {
      fits_report_error(stderr, status);
      return -1;
   }

   //copy base image hdu
   if (fits_copy_hdu(bptr, ofptr, 0, &status)) {
      fits_report_error(stderr, status);
      return -1;
   }

   //update COVERAGE keyword
   if (fits_update_key(ofptr, TDOUBLE, "COVERAGE", &coverage,
            "Percent Populated Pixels", &status)) {
      fits_report_error(stderr, status);
      return -1;
   }

   //write pixel data
   if (fits_write_img(ofptr, TDOUBLE, 1, bnpix, bimg, &status)) {
      fits_report_error(stderr, status);
      if (status != 412){return -1;};
      status=0;
   }

   //done writing
   if (fits_close_file(ofptr, &status)) {
      fits_report_error(stderr, status);
      if (status != 412){return -1;};
      status=0;
   }

   //--------------------------------
   // Free memory from original image
   //--------------------------------
   if (fits_close_file(bptr, &status)) {
      if (0 != isVerbose) {
         printf("error freeing original image\n");
         fits_report_error(stderr, status);
         if (status != 412){return -1;};
      }
      fits_clear_errmsg();
      status=0;
   }
   free(bimg);
   free(contribArr);

   return 0;
}
