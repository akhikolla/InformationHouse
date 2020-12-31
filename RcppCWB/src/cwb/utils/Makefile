##   -*-Makefile-*-
## 
##  IMS Open Corpus Workbench (CWB)
##  Copyright (C) 1993-2006 by IMS, University of Stuttgart
##  Copyright (C) 2007-     by the respective contributers (see file AUTHORS)
## 
##  This program is free software; you can redistribute it and/or modify it
##  under the terms of the GNU General Public License as published by the
##  Free Software Foundation; either version 2, or (at your option) any later
##  version.
## 
##  This program is distributed in the hope that it will be useful, but
##  WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
##  Public License for more details (in the file "COPYING", or available via
##  WWW at http://www.gnu.org/copyleft/gpl.html).

## Read configuration settings and standard definitions
TOP = $(R_PACKAGE_SOURCE)
include $(TOP)/config.mk

# the following targets are available
#
#  all		compile tools
#  clean        delete object files and binaries
#  realclean    also deleted automatically generated parsers and dependencies
#  depend	update dependencies
#  install	install tool binaries
#  uninstall	uninstall tools from chosen location (currently not supported)
#  relase       install to binary release dir
#  size		print size of source code (line counts)
#

.PHONY: all clean realclean depend install uninstall size

## ----------------------------------------------------------------------
## CWB command-line utilities  headers / sources / binaries 

SRCS =  _cwb_makeall.c \
	_cwb_huffcode.c _cwb_compress_rdx.c \

OBJS =  _cwb_makeall.o \
	_cwb_huffcode.o _cwb_compress_rdx.o \


all: libcwb.a

libcwb.a: $(OBJS)
	@$(ECHO) "--------------------------------- CREATING ARCHIVE"
	$(RM) $@
	$(AR) cq $@ $^
#	$(AR) -s $@

_cwb_makeall.o: _cwb_makeall.c
	${CC} -c ${CFLAGS} -o _cwb_makeall.o _cwb_makeall.c

_cwb_huffcode.o: _cwb_huffcode.c
	${CC} -c ${CFLAGS} -o _cwb_huffcode.o _cwb_huffcode.c

_cwb_compress_rdx.o: _cwb_compress_rdx.c
	${CC} -c ${CFLAGS} -o _cwb_compress_rdx.o _cwb_compress_rdx.c

clean:
	$(RM) *.o *~

realclean:	clean
	-$(RM) depend.mk


