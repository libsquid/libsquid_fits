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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#include <libsquid.h>
#include <libsquid_wcs.h>

#define SQUIDMERGE_VERSION 0.5.0
#define SQUIDMERGE_MAJOR 0
#define SQUIDMERGE_MINOR 5
#define SQUIDMERGE_PATCH 0
#define SQUIDMERGE_RELEASE ""

int main(int argc, char *argv[]) {
  int i;
  fitsfile *bptr, *fptr; // base and input pointers
  double *bimg, *fimg; // base and input images
  int status=0; // status of cfitsio function calls
  int bitpix; // image parameters
  long bnaxes[2], fnaxes[2]; // image dimensions
  int bnaxis, fnaxis; // number of image axes
  long bnpix, fnpix; // number of pixels in img
  long bmid, fmid; // map id
  long pix,omiss; // pixel index, misses
  double coverage; // coverage in percent
  char tmpname[FNAME_MAX],outname[FNAME_MAX];
  char compstr[FNAME_MAX];
  fitsfile *ofptr; // output file pointer

  if (argc < 3) {
    printf("Calling sequence: %s infile1 ... infileN outfile\n",argv[0]);
    return(0);
  }
  //strncpy(outname,"test.fit",FNAME_MAX);
  strncpy(outname,argv[argc-1],FNAME_MAX);
  printf("output file is %s\n",outname);
  strncpy(compstr,"",FNAME_MAX);

  // Load original image and header
  //snprintf(tmpname,FNAME_MAX,"%s[1]",argv[1]); // for compressed (see also inside of for loop)
  snprintf(tmpname,FNAME_MAX,"%s",argv[1]); // for uncompressed (see also inside of for loop)
  //printf("opening %s\n",tmpname);
  if (fits_open_file(&bptr, tmpname, READONLY, &status)) {
    fprintf(stderr,"fits_open_file failed in %s\n",argv[0]);
    fits_report_error(stderr, status);
    exit(-1);
  }
  if (fits_get_img_param(bptr, 2, &bitpix, &bnaxis, bnaxes, &status)) {
    fprintf(stderr,"fits_get_img_param failed in %s\n",argv[0]);
    fits_report_error(stderr, status);
    exit(-1);
  }
  //printf("parent naxes are %ld %ld\n",bnaxes[0],bnaxes[1]);
  bnpix = bnaxes[0]*bnaxes[1];
  //printf("parent npix is %ld\n",bnpix);
  
  bimg=malloc(bnpix*sizeof(double)); // allocate image array
  if (bimg == NULL) {
    fprintf(stderr,"malloc failed in %s, %s\n",argv[0],strerror(errno));
    exit(-1);
  }
  if (fits_read_img(bptr, TDOUBLE, 1, bnpix, NULL, bimg, NULL, &status)) {
    printf("error in fits_read_img\n");
    fits_report_error(stderr, status);
    exit(-1);
  }
  if (fits_read_key(bptr, TLONG, "MAPID", &bmid, NULL, &status)) {
    printf("keyword MAPID not in header.\n");
    fits_report_error(stderr, status);
    exit(-1);
  }
  if (fits_read_key(bptr, TDOUBLE, "COVERAGE", &coverage, NULL, &status)) {
    printf("keyword COVERAGE not in header.\n");
    fits_report_error(stderr, status);
    exit(-1);
  }

  //***************************************************************
  // loop over image arguments and add them to base image
  //***************************************************************
  for (i=2; i < argc-1; i++) {

    printf("%s\n",argv[i]);

    // Load original image and header
    //snprintf(tmpname,FNAME_MAX,"%s[1]",argv[i]); // for compressed
    snprintf(tmpname,FNAME_MAX,"%s",argv[i]); // for uncompressed
    if (fits_open_file(&fptr, tmpname, READONLY, &status)) {
      fits_report_error(stderr, status);
      exit(-1);
    }
    if (fits_get_img_param(fptr, 2, &bitpix, &fnaxis, fnaxes, &status)) {
      fits_report_error(stderr, status);
      exit(-1);
    }
    //printf("naxes are %ld %ld\n",fnaxes[0],fnaxes[1]);
    fnpix = fnaxes[0]*fnaxes[1];
    //printf("npix is %ld\n",fnpix);

    fimg=malloc(fnpix*sizeof(double)); // allocate image array
    if (fimg == NULL) {
      fprintf(stderr,"malloc failed in %s, %s\n",argv[0],strerror(errno));
      exit(-1);
    }
    if (fits_read_img(fptr, TDOUBLE, 1, fnpix, NULL, fimg, NULL, &status)) {
      fits_report_error(stderr, status);
      exit(-1);
    }
    if (fits_read_key(fptr, TLONG, "MAPID", &fmid, NULL, &status)) {
      printf("keyword MAPID not in header.\n");
      fits_report_error(stderr, status);
      exit(-1);
    }

    // check if we are in the same map id and images are same size
    if ((bmid == fmid)&&(bnpix == fnpix)) {
      // add input image pixels to base image
      omiss=0;
      for (pix = 0; pix < bnpix; pix++) {
	if (bimg[pix] < 10) {
	  if (fimg[pix] > 10) {
	    bimg[pix]=fimg[pix];
	  } else {
	    bimg[pix]=0;
	    omiss++;
	  }
	}
      }
    } else {
      // skip adding this image
      printf("skipping!\n");
    }

    if (fits_close_file(fptr, &status)) {
      //printf("error freeing original image\n");
      //fits_report_error(stderr, status);
      //if (status != 412) exit(-1);
      fits_clear_errmsg();
      status=0;
    }
    free(fimg);
  }

  coverage=100.0*(bnpix-omiss)/bnpix;
  //printf("coverage=%.1f\n",coverage);

  //**********************************************************
  // Write output image
  //**********************************************************

  unlink(outname); // remove file if it already exists
  snprintf(tmpname,FNAME_MAX,"%s%s",outname,compstr);
  //printf("writing %s\n",tmpname);
  if (fits_create_file(&ofptr, tmpname, &status)) { // if you want to avoid compression
    fits_report_error(stderr, status);
    exit(-1);
  }
  if (fits_copy_hdu(bptr,ofptr,0,&status)) { // copy original image
    fits_report_error(stderr, status);
    exit(-1);
  }
  // Update COVERAGE keyword
  if (fits_update_key(ofptr, TDOUBLE, "COVERAGE", &coverage,
		      "Percent Populated Pixels", &status)) {
    fits_report_error(stderr, status);
    return(-1);
  }
  if (fits_write_img(ofptr, TDOUBLE, 1, bnpix, bimg, &status)) {
      fits_report_error(stderr, status);
      if (status != 412) exit(-1);
      status=0;
    }
  if (fits_close_file(ofptr, &status)) {
    fits_report_error(stderr, status);
    if (status != 412) exit(-1);
    status=0;
  }

  //***********************************************
  // Free memory from original image
  //***********************************************
  if (fits_close_file(bptr, &status)) {
    //printf("error freeing original image\n");
    //fits_report_error(stderr, status);
    //if (status != 412) exit(-1);
    fits_clear_errmsg();
    status=0;
  }
  free(bimg);

  return(0);
}
