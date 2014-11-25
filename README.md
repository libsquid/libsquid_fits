libsquid_fits
=============

Provides utility program "squidfits" for re-projecting FITS images with WCS information onto SQUID tiles.

=================================================

This repository provides code for a utility called "squidfits".
The squidfits program will take a FITS image with a valid WCS header
and break it up into squid tile images (also FITS) which are re-projected
onto one of the four supported quad-cube projections.  Optionally, the
resulting squid tile images are compressed (non-lossy) as well.

An additional utility, "squidmerge", will merge squid tile images with
the same SQUID index.

The build system uses cmake.  The libsquid, libsquid_wcs, cfitsio, and
wcslib libraries are required for compilation.  See the INSTALL file for
more information.

For more information on the LibSQUID system, see http://libsquid.github.io,
the wiki, and the pdf documentation in the libsquid_doc repository.

This software is released under the LGPL license.  If you use it, please
let me know via email to jwren@lanl.gov.  Enjoy!
