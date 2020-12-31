# RcppCWB 0.3.0

* To avoid warning when running R CMD check, the http://pcre.org is used rather than https://pcre.org in the DESCRIPTION and the README file.
* To overcome a somewhat dirty solution for multiple symbol definitions, adding the 
'fcommon' flag to the CFLAGS in the configure script has been removed. The C code 
has been modified such that multiple symbol definitions are omitted.
* The macOS image used for test on Travis CI is now 'xcode9.4'
* On Solaris, the configure script would define the flag "-Wl,--allow-multiple-definition" to be passed 
to the linker flags. The rework of the CWB includes and the inclusion of the header file 'env.h' makes it
possible to drop this flag. It was defined at a confusing place anyway.
* Using the compiler desired by the user (in Makeconf, Makevars file) is now there for all OSes.
* If pkg-config is not present on macOS, a warning is issued; the user gets the advice to use the brew 
package manager to install pkg-config.
* There is an explicit check in the configure script whether the dependencies ncurses, pcre and glib-2.0
are present. If not, a telling error with installation instructions is displayed.
* When unloading the package, the dynamic library RcppCWB.so is unloaded.
* When loading the package, CQP is initialized by default (call `cqp_initialize()`)


# RcppCWB 0.2.9

* Starting with GCC 10, the compiler defaults to -fno-common, resulting in error messages during the linker stage, see [the change log of the GCC compiler](https://gcc.gnu.org/gcc-10/changes.html). To address this issue, the -fcommon option is now used by default when compiling the CWB C files on Linux 64bit systems. The CWB code includes header files multiple times, causing multiple definitions.
* On Linux systems, the hard-coded definition as the preferred C compiler in the CWB configuration sripts will be replaced by what the CC variable defines (in ~/.R/Makevars or the Makeconf file, the result returned by R CMD config CC).
* Remaining bashisms have been removed from the cleanup file. The shebang line of the
cleanup and the configure file is now #!/bin/sh, to avoid any reliance on bash.

# RcppCWB 0.2.8

* There have been (minor) modifiations of the C code of the CWB so that compilation succeeds on Solaris.
* Using the '-C' flag in the CWB Makefiles has been replaced by 'cd cl' / 'cd cqp' to avoid dependence on GNU make. GNU make is still required, because of 'include' statements in the Makefiles.
* Removed an action on 'depend.mk' from 'cleanup' script to avoid error messages that depend.mk is not present when Makefiles are first loaded.
* Dummy depend.mk files will satisfy include statement in Makefiles when running 'make clean' (depend.mk files are created only when running depend.mk)
* For creating index of static archives (libcl, libcqb, libcwb), a call to 'ranlib' has been replaced by an equivalent 'ar -s' in the Makefiles, but commented out.
* In the platform-specific config files of the CWB, the '-march'-option has been taken out, to safeguard portability.
* To meet the requirements of the upcoming changes in the CRAN check process to use staged installs, the procedure to reset the paths in the test data within the package has been replaced throughout by using a temporary registry directory. The `get_tmp_registry()` will return the whereabouts of this directory.


# RcppCWB 0.2.7

* If glib-2.0 is not present on macOS, binaries of the static library and 
header files are downloaded from a GitHub repo. This prepares to get RcppCWB
pass macOS checking on CRAN machines.
* A slight modification of the C code will now prevent previous crashes resulting
from a faulty CQP syntax. The solution will not yet be effective for Windows
systems until we have recompiled the libcqp static library that is downloaded
during the installation process.
* A new C++-level function 'check_corpus' checks whether a given corpus is
available and is used by the `check_corpus()`-function. Problems with 
the previous implementation that relied on files in the registry directory to
ensure the presence of a corpus hopefully do not occur.
* Calling the 'find_readline.perl' utility script is omitted on macOS, so 
previous warning messages when running the makefile do not show up any more.

# RcppCWB 0.2.6

* Function `cl_charset_name()` is exposed, it will return the charset of a 
corpus. Faster than parsing the registry file again and again.
* A new `cl_delete_corpus()`-function can remove loaded corpora from memory.


# RcppCWB 0.2.5

* In Makevars.win, libiconv is explicitly linked, to make RcppCWB compatible with new
release of Rtools.
* regex in check_s_attribute() for parsing registry file improved so that it does not
produce an error if '# [attribute]' follows after declaration of s_attribute


# RcppCWB 0.2.4
* for linux and macOS, CWB 3.4.14 included, so that UTF-8 support is realized
* bug removed in check_cqp_query that would prevent special characters from working
in CQP queries
* check_strucs, check_cpos and check_id are checking for NAs now to avoid crashes
* cwb command line tools cwb-makeall, cwb-huffcode and cwb-compress-rdx exposed
  as cwb_makeall, cwb_huffcode and cwb_compress_rdx

# RcppCWB 0.2.3

* when loading the package, a check is performed to make sure that paths in the 
registry files point to the data files of the sample data (issues may occur when
installing binaries)
* auxiliary functions to check whether input to Rcpp-wrappers/C functions is valid
are now exported and documented
* more consistent validity checks of input to functions for structural attributes


# RcppCWB 0.2.2

* Compiling RcppCWB on unix-like systems (macOS, Linux) will work now without
the presence of glib (on Windows, the dependency persists).
* The presence of the bison parser is not required any more. The package includes 
the C source generated by the bison parser along with the original input files.
* Functionality to generate CWB-indexed corpora and to generate and manipulate
the registry file describing a corpus has been moved to a new package 'cwbtools'
(see https://www.github.com/PolMine/cwbtools) in order to maintain a clearly
defined scope of RcppCWB to expose functionality of the C code of the CWB.
* Minor intervention in function 'valid_subcorpus_name' to omit a -Wtautological-pointer-compare warning leading to a WARNING when checking package
for R 3.5.0 with option --as-cran

# RcppCWB 0.2.1

* In previous versions the drive of the working directory and of the 
registry/data directory had to be identical on Windows; this limitation 
does not persist;
* Some utility functions could be removed that were necessary to check the
identity of the drives of the working directory and the data.


# RcppCWB 0.2.0

* In addition to low-level functionality of the corpus library (CL), functions
of the Corpus Query Processor (CQP) are exposed, building  on C wrappers in the
rcqp package;
* The authors of the rcqp package (Bernard Desgraupes and Sylvain Loiseau) are
mentioned as package authors and as authors of functions using CQP, as the code
used to expose CQP functionality is a modified version of rcqp code;
* Extended package description explaining the rationale for developing the 
RcppCWB package;
* Documentation of functions has been rearranged, many examples have been 
included;
* Renaming of exposed functions of corpus library from cwb_... to cl_...;
* sanity checks in R wrappers for Rcpp functions.


# RcppCWB 0.1.7

* CWB source code included in package to be GPL compliant
* template to adjust HOME and INFO in registry file used (tools/setpaths.R)
* using VignetteBuilder has been removed
* definition of Rprintf in cwb/cl/macros.c

# RcppCWB 0.1.6

* now using configure/configure.win script in combination with setpaths.R


# RcppCWB 0.1.1

* vignette included that explains cross-compiling CWB for Windows
* check in struc2str to ensure that structure has attributes


# RcppCWB 0.1.0

* Windows compatibility (potentially still limited)
