//
// Program to break up a fits image with valid wcs into SQUID tiles.
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

#define SQUIDFITS_VERSION "0.6.1"
#define SQUIDFITS_MAJOR 0
#define SQUIDFITS_MINOR 6
#define SQUIDFITS_PATCH 1 
#define SQUIDFITS_RELEASE ""

int main(int argc, char *argv[]) {
  fitsfile *fptr, *ofptr; // pointers to fits images
  char infile[FNAME_MAX]; // input file name
  char *bname,*bname0,*bnamehead; // input file basename
  char bnamehead2[FNAME_MAX], outname[FNAME_MAX], outnamecomp[FNAME_MAX]; // output file name
  extern char *optarg;
  extern int optind;
  int optc,opterror; // for getopt call
  int do_compress=0; // compress output image
  int itype; // interpolation type, 0=NEAREST, 1=BILINEAR (default), 2=CSPLINE, 3=CCONVOL
  int projection; // projection type, 0=TSC, 1=CSC, 2=QSC, 3=HSC
  int pixtype=0; // pixel type, 1=INT16, 2=UINT16, 3=INT32, 4=UINT32, 5=FLOAT, 6=DOUBLE
  int pixouttype; // cfitsio image type write code
  int flipx; // flip x axis
  long i; // loop counter
  squid_type hx,hy,tside;
  double hxd,hyd;
  int naxis; // number of image axes
  int status=0; // status of cfitsio function calls
  long naxes[2], npix, onaxes[2]; //number of pix per axis and total number of pix
  int bitpix; // image parameters
  double cdelt; // deg/pix
  double xc,yc,rac,decc,prac,pdecc,pra0,pdec0,sdist; // ra and dec in degrees
  struct wcsprm *wcs, *hpxwcs; // wcs coordinate struct
  long opix; // output pix index, x+y*nrows;
  double interpval; // interpolated pixel value
  int interpout; // interpolation return value
  double *img; // fits image
  int16_t *outimg_int16; // output image, short
  uint16_t *outimg_uint16; // output image, unsigned short
  int32_t *outimg_int32; // output image, int
  uint32_t *outimg_uint32; // output image, unsigned int
  float *outimg_float; // output image float
  double *outimg_double; // output image double
  long outpix; // number of pixels in output image
  double ra,dec; // image coords in deg
  double ix,iy; // input image coords
  char *header; // full fits header string
  int nkeyrec, nreject, nwcs; // # header keys, # keys rejected by wcspih, # wcs found
  double pixcrd[2],imgcrd[2],phi[1],theta[1],wcor[2]; // args to wcslib functions
  int wstat[1];
  squid_type hpix; // squid tile nside, and pixel number
  int hres1; // squid res of full thumbnail and thumbnail pix
  int nhid; // number of chars in healpix id (log10(Npix))
  char hfmt[FNAME_MAX],hfmt2[FNAME_MAX]; // strings for contructing variable format
  long hidcount;
  squid_type hidarr[ARR_MAX], hidix;
  double parentx, parenty; // x,y from parent image of helpix tile center
  double parentr; // sph dist of healpix center to parent center
  long omiss; // output image population count
  double coverage; // percent output image populated
  char compstr[FNAME_MAX]; // compression string for filename
  char dateobs[100];
  //struct wcserr fixinfo[NWCSFIX];
  double tside0, fluxratio; // double version of tside and flux ratio param
  int rval; // cfitsio return val
  long bzval; // bzero value for header

  // Mask variables
  int            use_mask;
  char           fname_mask[FNAME_MAX];
  fitsfile*      fits_mask;
  int            naxis_mask;
  long           naxes_mask[2];
  int            bitpix_mask;
  unsigned char* buf_mask;
  unsigned char  pixel_mask;
  long           coords_mask[2];

  // Init mask
  use_mask      = 0;
  fname_mask[0] = 0;
  fits_mask     = NULL;
  naxis_mask    = 0;
  naxes_mask[0] = 0;
  naxes_mask[1] = 0;
  bitpix_mask   = 0;
  buf_mask      = NULL;
  pixel_mask    = 0;

  // Parse command line opts
  opterror=0;
  hres1=-1;
  tside=0;
  itype=BILINEAR; // default of BILINEAR interpolation
  flipx=0; // don't flip x axis by default
  projection=HSC; // default HSC projection
  while ((optc = getopt(argc, argv, "cfi:m:p:r:t:b:")) != -1)
    switch (optc) {
    case 'c':
      do_compress=1;
      printf("Compressing fits output.\n");
      break;
    case 'f':
      flipx=1;
      break;
    case 'i':
      // interpolation type, 0=NEAREST, 1=BILINEAR (default), 2=CSPLINE, 3=CCONVOL
      itype=atoi(optarg);
      if (itype == NEAREST) {
        printf("Interpolation type set to NEAREST.\n");
      } else if (itype == BILINEAR) {
        printf("Interpolation type set to BILINEAR.\n");
      } else if (itype == CSPLINE) {
        printf("Interpolation type set to CUBIC SPLINE.\n");
      } else if (itype == CCONVOL) {
        printf("Interpolation type set to CUBIC CONVOLUTION.\n");
      } else {
        printf("Interpolation type unknown, setting to BILINEAR.\n");
        itype=BILINEAR;
      }
      break;
    case 'm':
      // use pixel mask
      use_mask = 1;
      strncpy(fname_mask, optarg, FNAME_MAX);
      break;
    case 'p':
      // projection type, 0=TSC, 1=CSC, 2=QSC, 3=HSC (default)
      projection=atoi(optarg);
      if (projection == TSC) {
	printf("Projection is TSC.\n");
      } else if (projection == CSC) {
	printf("Projection is CSC.\n");
      } else if (projection == QSC) {
	printf("Projection is QSC.\n");
      } else if (projection == HSC) {
	printf("Projection is HSC.\n");
      } else {
	printf("Projection type unknown, setting to HSC.\n");
	projection=HSC;
      }
      break;
    case 'b':
      // image pixel type, 0=INT16, 1=UINT16, 2=INT32, 3=UINT32, 4=FLOAT, 5=DOUBLE
      // This is being done because currently cfisio (v3.36) isn't doing compression
      // when doing implicit type conversions.
      pixtype=atoi(optarg);
      if (pixtype == 1) {
         printf("Pixel type is INT16.\n");
         pixouttype=SHORT_IMG;
      } else if (pixtype == 2) {
         printf("Pixel type is UINT16.\n");
         pixouttype=USHORT_IMG;
      } else if (pixtype == 3) {
         printf("Pixel type is INT32.\n");
	 pixouttype=LONG_IMG;
      } else if (pixtype == 4) {
         printf("Pixel type is UINT32.\n");
	 pixouttype=ULONG_IMG;
      } else if (pixtype == 5) {
         printf("Pixel type is FLOAT.\n");
	 pixouttype=FLOAT_IMG;
      } else if (pixtype == 6) {
         printf("Pixel type is DOUBLE.\n");
	 pixouttype=DOUBLE_IMG;
      } else {
         printf("Pixel type is DOUBLE.\n");
         pixtype=6;
	 pixouttype=DOUBLE_IMG;
      }
    case 'r':
      hres1=atoi(optarg);
      printf("Setting tile resolution to %d\n",hres1);
      break;
    case 't':
      tside=(squid_type)atol(optarg);
      printf("Setting tile nside to %ld\n",(long)tside);
      break;
    case '?':
      if ((optopt == 'i')||
          (optopt == 'm')||
          (optopt == 'r')||
          (optopt == 't')) {
        printf("Option -%c requires an argument.\n", optopt);
      }
      printf("Bad option '-%c' detected.\n", optopt); 
      opterror=1;
      break;
    }
  if (pixtype == 0) {
    printf("Pixel type is DOUBLE.\n");
    pixtype=6;
    pixouttype=DOUBLE_IMG;
  }  
  if ((opterror)||((optind+1) > argc)) {
    printf("\n");
    printf("squidfits version %s%s\n",SQUIDFITS_VERSION,SQUIDFITS_RELEASE);
    printf("\n");
    printf("Example usage...\n");
    printf("%s [args] infile.fit\n",argv[0]);
    printf("\n");
    printf("\t-c\t\tcompress fits output images.\n");
    printf("\t-f\t\tflip x axis in output images.\n");
    printf("\t-i ival\t\tinterpolation type:\n");
    printf("\t\t\t0 = Nearest neighbor\n");
    printf("\t\t\t1 = Bilinear (default)\n");
    printf("\t\t\t2 = Cubic spline\n");
    printf("\t\t\t3 = Cubic convolution\n");
    printf("\t-m fname_mask\tfilename of mask.\n");
    printf("\t-p\t\tprojection type:\n");
    printf("\t\t\t0 = TSC\n");
    printf("\t\t\t1 = CSC\n");
    printf("\t\t\t2 = QSC\n");
    printf("\t\t\t3 = HSC (default)\n");
    printf("\t-b\t\toutput image data type\n");
    printf("\t\t\t1 = int16\n");
    printf("\t\t\t2 = uint16\n");
    printf("\t\t\t3 = int32\n");
    printf("\t\t\t4 = uint32\n");
    printf("\t\t\t5 = float\n");
    printf("\t\t\t6 = double (default)\n");
    printf("\t-r rval\t\ttile resolution parameter, default 3.\n");
    printf("\t-t tval\t\t# pix on side of tile, calculated if not set.\n");
    printf("\n");
    exit(-1);
  }
  strncpy(infile,argv[optind],FNAME_MAX);
  if (hres1 < 0) {
    printf("Setting tile resolution to default value of 3.\n");
    hres1=3; // whatever, should always set this on command line
  }

  // Load mask
  if (1 == use_mask) {
    // Open mask
    if (fits_open_file(&fits_mask, fname_mask, READONLY, &status)) {
      fprintf(stderr,"mask:  error in fits_open_file in %s\n",argv[0]);
      fits_report_error(stderr, status);
      exit(-1);
    }

    // Read image parameters
    if (fits_get_img_param(fits_mask, 2, &bitpix_mask, &naxis_mask, naxes_mask, &status)) {
      fprintf(stderr,"mask:  error in fits_get_img_param in %s\n",argv[0]);
      fits_report_error(stderr, status);
      exit(-1);
    }

    // Verify bitpix type
    if (BYTE_IMG != bitpix_mask) {
      fprintf(stderr,"mask:  error:  unhandeld bitpix:  %i\n", bitpix_mask);
      exit(-1);
    }

    // Read data
    buf_mask = malloc(naxes_mask[0] * naxes_mask[1] * sizeof(unsigned char));
    if (buf_mask == NULL) {
      fprintf(stderr,"malloc failed in %s, %s\n",argv[0],strerror(errno));
      exit(-1);
    }
    if (fits_read_img(fits_mask, TBYTE, 1, naxes_mask[0] * naxes_mask[1], NULL, buf_mask, NULL, &status)) {
      fprintf(stderr,"mask:  error in fits_read_img in %s\n",argv[0]);
      fits_report_error(stderr, status);
      exit(-1);
    }

    // Close mask
    fits_close_file(fits_mask, &status);
    fits_mask = NULL;
    status    = 0;
  }

  // Load original image and header
  if (fits_open_file(&fptr, infile, READONLY, &status)) {
    fits_report_error(stderr, status);
    exit(-1);
  }
  if (fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status)) {
    fits_report_error(stderr, status);
    exit(-1);
  }
  npix = naxes[0]*naxes[1];
  printf("parent npix is %ld\n",npix);
  img=malloc(npix*sizeof(double)); // allocate image array
  if (img == NULL) {
    fprintf(stderr,"malloc failed in %s, %s\n",argv[0],strerror(errno));
    exit(-1);
  }
  if (fits_read_img(fptr, TDOUBLE, 1, npix, NULL, img, NULL, &status)) {
    fits_report_error(stderr, status);
    exit(-1);
  }

  // Verify mask dimensions
  if (1 == use_mask) {
    if (naxes[0] != naxes_mask[0]) {
      fprintf(stderr,"mask to input size mismatch axis 0 in %s\n",argv[0]);
      exit(-1);
    }
    if (naxes[1] != naxes_mask[1]) {
      fprintf(stderr,"mask to input size mismatch axis 1 in %s\n",argv[0]);
      exit(-1);
    }
  }

  // Remove SCAMP PV* headers because they screw up wcslib
  for (i=0; i<1000; i++) {
    if (fits_delete_key(fptr,"PV?_*", &status)) {
      fits_clear_errmsg();
      status=0;
      break;
    }
  }
  // Preserve DATE-OBS, lost in wcspih
  if (fits_read_card(fptr, "DATE-OBS", dateobs, &status)) {
    fits_report_error(stderr, status);
    exit(-1);
  }
  // Now convert header into string for wcs parsing  
  if (fits_hdr2str(fptr, 1, NULL, 0, &header, &nkeyrec, &status)) {
    fits_report_error(stderr, status);
    exit(-1);
  }
  // Interpret the WCS keywords, deleting them afterwards
  if ((status = wcspih(header, nkeyrec, WCSHDR_all, -3, &nreject, &nwcs, &wcs))) {
    fprintf(stderr, "wcspih ERROR %d: %s.\n", status, wcshdr_errmsg[status]);
    fprintf(stderr,"%d rejected keywords\n",nreject);
    exit(-1);
  }

  // Get image center coordinates
  xc=round((double)naxes[0]/2.0);
  yc=round((double)naxes[1]/2.0);
  pixcrd[0]=xc;
  pixcrd[1]=yc;
  printf("xc=%.5f yc=%.5f\n",xc,yc);
  if ((status=wcsp2s(wcs,1,2,pixcrd,imgcrd,phi,theta,wcor,wstat)) > 0) {
    fprintf(stderr, "wcsp2s returned status=%d\n", status);
    exit(-1);
  }
  prac=wcor[0]; // parent center ra
  pdecc=wcor[1]; // parent center dec
  printf("parent rac=%f decc=%f\n",prac,pdecc);
  /* Just for testing reverse operation */
  if ((status=wcss2p(wcs,1,2,wcor,phi,theta,imgcrd,pixcrd,wstat)) > 0) {
    fprintf(stderr,"wcss2p returned status=%d\n", status);
    exit(-1);
  }
  pixcrd[0]=0.0;
  pixcrd[1]=0.0;
  if ((status=wcsp2s(wcs,1,2,pixcrd,imgcrd,phi,theta,wcor,wstat)) > 0) {
    fprintf(stderr, "wcsp2s returned status=%d\n", status);
    exit(-1);
  }
  pra0=wcor[0]; // parent corner ra
  pdec0=wcor[1]; // parent corner dec
  sphdist(prac*DD2R,pdecc*DD2R,pra0*DD2R,pdec0*DD2R,&sdist);
  cdelt=sdist/(sqrt(xc*xc+yc*yc)*DD2R);
  printf("cdelt=%.5f\n",cdelt);

  // Try getting squid ids
  wcs_getsquids(projection,wcs,cdelt,naxes,hres1,hidarr,&hidcount);
  printf("squid count=%ld\n",hidcount);

  // Set up filename base strings
  bname=basename(infile);
  bname0=bname; // save pointer before strtok call
  bnamehead=strtok(bname0,".");
  strncpy(bnamehead2,bnamehead,FNAME_MAX-50);
  if (projection == TSC) {
    strncpy(hfmt,"%s_tsc%0",FNAME_MAX);
  } else if (projection == CSC) {
    strncpy(hfmt,"%s_csc%0",FNAME_MAX);
  } else if (projection == QSC) {
    strncpy(hfmt,"%s_qsc%0",FNAME_MAX);
  } else if (projection == HSC) {
    strncpy(hfmt,"%s_hsc%0",FNAME_MAX);
  } else {
    // Should never get here, this should have been checked previously.
    fprintf(stderr,"unknown projection in %s\n",argv[0]);
    exit(-1);
  }
  nhid = ceil(log10(pow(2,2*hres1+4))); // number of possible digits in squid

  //-----------------------------------------------------------------------
  // Loop over squids and create a thumbnail image for each one
  //-----------------------------------------------------------------------
  squid_tside(projection,hres1,cdelt,&tside0); // calculate tside based on input cdelt
  if (tside == 0) { // wasn't set on command line, so use tside0
    tside=(squid_type)round(tside0);
  }
  fluxratio=pow((tside0/tside),2); // flux ratio based on area ratio
  printf("tside0=%.5f tside=%ld fluxratio=%.5f\n",tside0,(long)tside,fluxratio);
  outpix=(long)pow(tside,2); // number of pixels in output image
  printf("\n");
  for (hidix=0; hidix<hidcount; hidix++) {

    // find thumbnail squid
    hpix=hidarr[hidix]; // squid index
    printf("squid=%ld\n",(long)hpix);
    squid2sph(projection,hpix,&rac,&decc);
    printf("squid ra=%f dec=%f\n",rac/DD2R,decc/DD2R);
    if (wcs_rd2pix(wcs,rac/DD2R,decc/DD2R,&parentx,&parenty) < 0) {
      fprintf(stderr,"wcs_rd2pix failed in %s\n",argv[0]);
      exit(-1);
    }
    sphdist(rac,decc,prac*DD2R,pdecc*DD2R,&parentr);

    // create temporary image
    if (pixtype == 1) {
      outimg_int16=calloc((size_t)outpix,sizeof(int16_t)); // don't forget to free this!
      if (outimg_int16 == NULL) {
        fprintf(stderr,"calloc failed in %s, %s\n",argv[0],strerror(errno));
        exit(-1);
      }
    } else if (pixtype == 2) {
      outimg_uint16=calloc((size_t)outpix,sizeof(uint16_t)); // don't forget to free this!
      if (outimg_uint16 == NULL) {
        fprintf(stderr,"calloc failed in %s, %s\n",argv[0],strerror(errno));
        exit(-1);
      }
    } else if (pixtype == 3) {
      outimg_int32=calloc((size_t)outpix,sizeof(int32_t)); // don't forget to free this!
      if (outimg_int32 == NULL) {
        fprintf(stderr,"calloc failed in %s, %s\n",argv[0],strerror(errno));
        exit(-1);
      }
    } else if (pixtype == 4) {
      outimg_uint32=calloc((size_t)outpix,sizeof(uint32_t)); // don't forget to free this!
      if (outimg_uint32 == NULL) {
        fprintf(stderr,"calloc failed in %s, %s\n",argv[0],strerror(errno));
        exit(-1);
      }
    } else if (pixtype == 5) {
      outimg_float=calloc((size_t)outpix,sizeof(float)); // don't forget to free this!
      if (outimg_float == NULL) {
        fprintf(stderr,"calloc failed in %s, %s\n",argv[0],strerror(errno));
        exit(-1);
      }
    } else if (pixtype == 6) {
      outimg_double=calloc((size_t)outpix,sizeof(double)); // don't forget to free this!
      if (outimg_double == NULL) {
        fprintf(stderr,"calloc failed in %s, %s\n",argv[0],strerror(errno));
        exit(-1);
      }
    }
    // get wcs struct for squid tile
    if (tile_getwcs(projection,hpix,tside,&hpxwcs) == -1) {
      fprintf(stderr,"tile_getwcs failed in %s\n",argv[0]);
      exit(-1);
    }

    //
    // Image population loop
    //
    omiss=0;
    for (hy=0; hy<tside; hy++) {
      for (hx=0; hx<tside; hx++) {
        opix=hx+hy*tside;
        hxd=(double)hx+0.5; // squid pix start at 0.5,0.5
        hyd=(double)hy+0.5; // squid pix start at 0.5,0.5

	// Get ra,dec for squid tile x,y
        if (tile_xy2sph(projection,hpix,hxd,hyd,tside,&ra,&dec) < 0) {
          //fprintf(stderr, "tile_xy2sph a failed\n");
          omiss++;
          continue;
        }

	// Get image pixel coords from wcs info
        if (wcs_rd2pix(wcs,ra/DD2R,dec/DD2R,&ix,&iy) < 0) {
          //fprintf(stderr, "wcs_rd2pix failed\n");
          omiss++;
          continue;
        }
	// Verify the image pixel coords are actually within the image
	if ((ix < 1)||(ix > naxes[0])||(iy < 1)||(iy > naxes[1])) {
	  omiss++;
	  continue;
	}

        // Verify input pixel may contribute
        if (1 == use_mask) {
          // Drop fractional component and shift (1, 1) to (0, 0)
          coords_mask[0] = ((long)floor(ix)) - 1;
          coords_mask[1] = ((long)floor(iy)) - 1;
          // Is x coord in bounds
          if ((0 <= coords_mask[0]) && (naxes_mask[0] > coords_mask[0])) {
            // Is y coord in bounds
            if ((0 <= coords_mask[1]) && (naxes_mask[1] > coords_mask[1])) {
              // Get pixel from mask
              pixel_mask = buf_mask[coords_mask[0] + (coords_mask[1] * naxes_mask[0])];
              // If masked, skip
              if (0 != pixel_mask) {
                omiss++;
                continue;
              }
            }
          }
        }

        // Find interpolated pixel value.
        // Note: wcs pix start at 1,1.
        if (flipx) {
          interpout=interp_img(itype,naxes[0]-ix+1.0,iy-1.0,img,naxes[0],naxes[1],&interpval);
        } else {
          interpout=interp_img(itype,ix-1.0,iy-1.0,img,naxes[0],naxes[1],&interpval);
        }
        if (interpout < 0) {
          omiss++;
        } else {
	  if (pixtype == 1) {
            outimg_int16[opix]=interpval*fluxratio;
          } else if (pixtype == 2) {
            outimg_uint16[opix]=interpval*fluxratio;
          } else if (pixtype == 3) {
            outimg_int32[opix]=interpval*fluxratio;
          } else if (pixtype == 4) {
            outimg_uint32[opix]=interpval*fluxratio;
          } else if (pixtype == 5) {
            outimg_float[opix]=interpval*fluxratio;
          } else {
            outimg_double[opix]=interpval*fluxratio;
	  }
        }
      }
    }
    coverage=100.0*(outpix-omiss)/outpix;
    printf("coverage=%.2f\n",coverage);

    //------------------------------------------------------------------------------
    // Write output image
    //------------------------------------------------------------------------------
    // create output image
    if (do_compress) {
      /*
      if (coverage > 50.0) {
        snprintf(compstr,FNAME_MAX,"[compress HS %d,%d; s 2.0]",
        (long)(tside/2),(long)(tside/2)); // h compression lossy
      } else {
        // Do non-lossy here because of problems with current cfitsio where it
        // adds a GZIPPED section for all-zero regions.
        // This currently screws up pyfits and ds9.
        // This also could potentially go away with later versions of cfitsio
        //snprintf(compstr,FNAME_MAX,"[compress H %d,%d]",tside,tside); // h compression non-lossy
        snprintf(compstr,FNAME_MAX,"[compress R %d,%d]",tside,tside); // rice compression
      }
      */
      snprintf(compstr,FNAME_MAX,"[compress R %d,%d]",tside,tside); // rice compression
    } else {
      compstr[0]='\0';
    }
    snprintf(hfmt2,FNAME_MAX,"%s%dd.fit",hfmt,nhid);
    snprintf(outname,FNAME_MAX,hfmt2,bnamehead2,hpix);
    snprintf(outnamecomp,FNAME_MAX,"%s%s",outname,compstr);
    printf("outname=%s\n",outnamecomp);
    unlink(outname); // remove file if it already exists
    if (fits_create_file(&ofptr, outnamecomp, &status)) {
      fits_report_error(stderr, status);
      exit(-1);
    }
    onaxes[0]=tside;
    onaxes[1]=tside;
    if (fits_create_img(ofptr, pixouttype, 2, onaxes, &status)) {
      fits_report_error(stderr, status);
      exit(-1);
    }

    // Copy headers from original image and add wcs info
    if (tile_addwcs(projection,hpix,hpxwcs,header,ofptr) == -1) {
      printf("tile_addwcs failed!\n");
    }
    // Update COVERAGE keyword
    if (fits_update_key(ofptr, TDOUBLE, "COVERAGE", &coverage,
                        "Percent Populated Pixels", &status)) {
      fits_report_error(stderr, status);
      return(-1);
    }
    // Add parent filename
    if (fits_update_key(ofptr, TSTRING, "PARENT", bname,"Parent Filename", &status)) {
      fits_report_error(stderr, status);
      return(-1);
    }
    // Update parent x and y coords
    if (fits_update_key(ofptr, TDOUBLE, "PARENTX", &parentx,"Parent X Coord for SQUID Center", &status)) {
      fits_report_error(stderr, status);
      return(-1);
    }
    if (fits_update_key(ofptr, TDOUBLE, "PARENTY", &parenty,"Parent Y Coord for SQUID Center", &status)) {
      fits_report_error(stderr, status);
      return(-1);
    }
    if (fits_update_key(ofptr, TDOUBLE, "PARENTR", &parentr,"Sphdist Parent Center to SQUID Center", &status)) {
      fits_report_error(stderr, status);
      return(-1);
    }
    // Update filename header
    if (fits_update_key(ofptr, TSTRING, "FILENAME", outname,"", &status)) {
      fits_report_error(stderr, status);
      return(-1);
    }
    // Update compstring header
    if (fits_update_key(ofptr, TSTRING, "COMPSTR", compstr,"compression output", &status)) {
      fits_report_error(stderr, status);
      return(-1);
    }
    // Restore DATE-OBS, lost in wcspih
    if (fits_update_card(ofptr, "DATE-OBS", dateobs, &status)) {
      fits_report_error(stderr, status);
      return(-1);
    }
    // Update SQDVER keyword
    if (fits_update_key(ofptr, TSTRING, "SQDVER", SQUIDFITS_VERSION,
                        "SQUIDFITS version string", &status)) {
      fits_report_error(stderr, status);
      return(-1);
    }
    // Update LSQDVER keyword
    if (fits_update_key(ofptr, TSTRING, "LSQDVER", LIBSQUID_VERSION,
                        "LibSQUID version string", &status)) {
      fits_report_error(stderr, status);
      return(-1);
    }
    // Update LSQDWVER keyword
    if (fits_update_key(ofptr, TSTRING, "LSQDWVER", LIBSQUIDWCS_VERSION,
                        "LibSQUID_WCS version string", &status)) {
      fits_report_error(stderr, status);
      return(-1);
    }

    // write output image
    if (pixtype == 1) { // write short image
      bzval=0;
      if (fits_update_key(ofptr, TLONG, "BZERO", &bzval, "default", &status)) {
        fits_report_error(stderr, status);
        return(-1);
      }
      rval=fits_write_img(ofptr, TSHORT, 1, outpix, outimg_int16, &status);
      free(outimg_int16);
    } else if (pixtype == 2) { // write unsigned short image
      bzval=32768;
      if (fits_update_key(ofptr, TLONG, "BZERO", &bzval, "offset for unsigned short", &status)) {
        fits_report_error(stderr, status);
        return(-1);
      }
      rval=fits_write_img(ofptr, TUSHORT, 1, outpix, outimg_uint16, &status);
      free(outimg_uint16);
    } else if (pixtype == 3) { // write int image
      bzval=0;
      if (fits_update_key(ofptr, TLONG, "BZERO", &bzval,"default", &status)) {
        fits_report_error(stderr, status);
        return(-1);
      }
      rval=fits_write_img(ofptr, TLONG, 1, outpix, outimg_int32, &status);
      free(outimg_int32);
    } else if (pixtype == 4) { // write unsigned int image
      bzval=2147483648;
      if (fits_update_key(ofptr, TLONG, "BZERO", &bzval, "offset for unsigned int", &status)) {
        fits_report_error(stderr, status);
        return(-1);
      }
      rval=fits_write_img(ofptr, TUINT, 1, outpix, outimg_uint32, &status);
      free(outimg_uint32);
    } else if (pixtype == 5) { // write float image
      bzval=0;
      if (fits_update_key(ofptr, TLONG, "BZERO", &bzval,"default", &status)) {
        fits_report_error(stderr, status);
        return(-1);
      }
      rval=fits_write_img(ofptr, TFLOAT, 1, outpix, outimg_float, &status);
      free(outimg_float);
    } else { // write double image
      bzval=0;
      if (fits_update_key(ofptr, TLONG, "BZERO", &bzval,"default", &status)) {
        fits_report_error(stderr, status);
        return(-1);
      }
      rval=fits_write_img(ofptr, TDOUBLE, 1, outpix, outimg_double, &status);
      free(outimg_double);
    }
    if (rval) {
      fits_report_error(stderr, status);
      if (status != 412) exit(-1);
      status=0;
    }
    if (fits_close_file(ofptr, &status)) {
      fits_report_error(stderr, status);
      if (status != 412) exit(-1);
      status=0;
    }

    // free memory
    //wcsfree(hpxwcs);
    printf("\n");

  } // ***** end of squidfits id for loop *********
  
  // free stuff from original image
  if (fits_close_file(fptr, &status)) {
    //printf("error freeing original image\n");
    //fits_report_error(stderr, status);
    //if (status != 412) exit(-1);
    fits_clear_errmsg();
    status=0;
  }

  if (NULL != buf_mask) {
    free(buf_mask);
    buf_mask = NULL;
  }
  free(img);
  free(header);
  wcsfree(wcs);

  return(0);
}

