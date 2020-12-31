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


##  This is the main Makefile for the CWB. 


#
#  Welcome to the IMS Open Corpus Workbench.  See "INSTALL" for detailed
#  installation instructions.  For a quick start, select appropriate PLATFORM
#  and SITE configurations in "config.mk", then build with "make all".
#  
#  The following make targets are available:
#
#  all          compile the Corpus Library, CQP, command-line utilities and man pages
#  install      install CL, CQP, utilities and man pages
# [uninstall    uninstall everything from specified locations] -- currently not available
#  release      create a release file (binaries and documentation)
#  clean        clean up in all subdirectories 
#  realclean    also delete automatically generated parsers, dependencies, etc.
#
#  size         check size of sourcecode (in lines)
#  depend       update dependencies (only necessary after major changes)
#  tags         create an emacs ETAGS file
#
#  cl           to compile the Corpus Library only
#  cqp          to compile CQP only
#  utils        to compile the utilities only
#  man          to update the manpages from their .pod sources
#  doxygen      to update the HTML code documentation (requires doxygen to be installed)
#  

TOP = $(R_PACKAGE_SOURCE)
include $(TOP)/config.mk

# targets for external libraries; currently not used (Glib and PCRE assumed to be on the system already)
SUBDIRS = cl cqp # subdirectories that have their own makefiles
SRCDIRS = cl cqp utils CQi	# subdirectories containing C source code


.PHONY: clean realclean depend all test install uninstall release mingw-libgnurx-2.5.1 cl cqp utils man instutils tags size


cl:
	@$(ECHO) "--------------------------------- BUILDING CORPUS LIBRARY (CL)"
	cd cl; $(MAKE)

cqp:
	@$(ECHO) "--------------------------------- BUILDING CQP"
	cd cqp; $(MAKE)

utils:
	@$(ECHO) "--------------------------------- BUILDING UTILITIES"
	cd utils; $(MAKE)

man:
	@$(ECHO) "--------------------------------- BUILDING MANPAGES"
	cd man; $(MAKE)


size:
	for i in $(SUBDIRS) ;\
	do \
	$(ECHO) "--------------------------------- $$i"; \
	cd $$i; \
	$(MAKE) size; \
	done;

depend:
	for i in $(SUBDIRS) ; do cd $$i; $(MAKE) depend.mk; cd ..; done;

clean:
	for i in $(SUBDIRS); do if [ -f "$$i/Makefile" ]; then cd $$i; $(MAKE) clean; cd ..; fi; done;
	-$(RM) *~ config/*/*~ CQi/*~ CQi/*.o technical/*~ TAGS
	-$(RM) -rf build

