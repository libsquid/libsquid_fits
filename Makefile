# -------------------------- LICENSE -----------------------------------
#
# This file is part of the LibSQUID software libraray.
#
#  LibSQUID is free software: you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License as published
#  by the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  LibSQUID is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with LibSQUID.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright 2014 James Wren and Los Alamos National Laboratory
#

TARGET_BINS = squidfits squidmerge

GCC     = gcc
CFLAGS  = -g -Wall -fPIC -I../libsquid -I../libsquid_wcs \
	$(shell pkg-config --cflags cfitsio) \
	$(shell pkg-config --cflags wcslib)
LDFLAGS = -L../libsquid -L../libsquid_wcs -lsquid -lsquid_wcs -lm \
	$(shell pkg-config --libs cfitsio) \
	$(shell pkg-config --libs wcslib)
LDFLAGS_STATIC = -L../libsquid -L../libsquid_wcs -lm \
	-Wl,-Bstatic -lsquid -lsquid_wcs \
	-Wl,-Bdynamic \
	$(shell pkg-config --libs cfitsio) \
	$(shell pkg-config --libs wcslib)

all: $(TARGET_BINS)

%.o: %.c
	$(GCC) -c $(CFLAGS) -o $@ $<

$(TARGET_BINS): % : %.o
	$(GCC) $(CFLAGS) $< $(LDFLAGS_STATIC) -o $@

clean:
	rm -f *.o $(TARGET_BINS)

