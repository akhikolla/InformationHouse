/*
 *  IMS Open Corpus Workbench (CWB)
 *  Copyright (C) 1993-2006 by IMS, University of Stuttgart
 *  Copyright (C) 2007-     by the respective contributers (see file AUTHORS)
 *
 *  This program is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation; either version 2, or (at your option) any later
 *  version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 *  Public License for more details (in the file "COPYING", or available via
 *  WWW at http://www.gnu.org/copyleft/gpl.html).
 */

#include <stdio.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/time.h>

#include "../cl/globals.h"
#include "../cl/macros.h"

#include "cqp.h"
#include "options.h"
#include "output.h"
#include "symtab.h"
#include "tree.h"
#include "eval.h"
#include "corpmanag.h"
#include "regex2dfa.h"
#include "ranges.h"
#include "macro.h"
#include "targets.h"

#include "parser.tab.h"


/** File handle used by the CQP-query-language parser. */
extern FILE *yyin;
/** Activates the CQP-query-language parser. */
extern int yyparse (void);
/** restarts the CQP-query-language parser. */
extern void yyrestart(FILE *input_file);


/** Array of file handles. Allows nested execution of "included" text files full of commands. */
FILE *cqp_files[MAXCQPFILES];
/** index into cqp_files. @see cqp_files */
int cqp_file_p;

/**
 * Boolean: true iff cqp_parse_file() - the main query syntax parsing function -
 * is currently reading from the cqprc file handler.
 *
 * @see cqprc
 * @see cqp_parse_file
 */
int reading_cqprc = 0;

/**
 * Global error status for CQP (will be returned to the caller when CQP exits).
 *
 * For the moment, it's simply boolean (0 = all OK, anything else = error).
 * Later, we might want to define some error macros (maybe have the CDA_
 * macros as a subset for errors bubbling up from CL, plus *also* extra
 * macros for other errors).
 *
 * TODO actually get this variable set in various error conditions.
 * TODO work out whether cqpserver should exit with this error.
 */
int cqp_error_status = 0;

/* ======================================== Query Buffer Interface */

char QueryBuffer[QUERY_BUFFER_SIZE];        /**< buffer for queries */
int QueryBufferP = 0;                       /**< index into this buffer, for appending */
int QueryBufferOverflow = 0;                /**< flag which signals buffer overflows */





/**
 * This is the signal handler function that is used for the interrupt signal (CTRL+C).
 *
 * @param signum  The signal number (ignored).
 */
static void
sigINT_signal_handler(int signum)
{
  if (!signal_handler_is_installed)
    exit(1); /* make sure we abort if Ctrl-C is pressed a second time (even on platforms where signal handlers don't need to be reinstalled) */

  if (EvaluationIsRunning) {
    Rprintf("** Aborting evaluation ... (press Ctrl-C again to exit CQP)\n");
    EvaluationIsRunning = 0;
  }
  signal_handler_is_installed = 0;
}

/**
 * Installs the interrupt signal handler function with the OS.
 *
 * This function installs a Ctrl-C interrupt handler (clears
 * EvaluationIsRunning flag). The function installed is
 * sigINT_signal_handler.
 *
 * @see sigINT_signal_handler
 */
void
install_signal_handler(void)
{
  signal_handler_is_installed = 1; /* pretend it's installed even if it wasn't done properly, so CQP won't keep trying */
  if (signal(SIGINT, sigINT_signal_handler) == SIG_ERR) {
    cqpmessage(Warning, "Can't install interrupt handler.\n");
    signal(SIGINT, SIG_IGN);
  }
}



/**
 * Initialises the CQP program (or cqpserver or cqpcl).
 *
 * This function:
 * - initialises the global variables;
 * - initialises the built-in random number generator;
 * - initialises the macro database;
 * - parses the program options;
 * - reads the initialisation file;
 * - reads the macro initialisation file;
 * - and loads the default corpus, if any.
 *
 * @param argc  The argc from main()
 * @param argv  The argv from main()
 * @return      Always 1.
 */
int
initialize_cqp(int argc, char **argv)
{
  char *home = NULL;
  char init_file_fullname[CL_MAX_FILENAME_LENGTH];
#ifdef __MINGW__
  char *homedrive = NULL;
  char *homepath = NULL;
#endif

  /* file handle for initialisation files, if any */
  FILE *cqprc;

  extern int yydebug;

  /* initialize global variables */

  exit_cqp = 0;
  cqp_file_p = 0;

  corpuslist = NULL;

  eep = -1;

  /* intialise built-in random number generator */
  cl_randomize();

  /* initialise macro database */
  init_macros();

  /* parse program options */
  parse_options(argc, argv);

  /* let's always run stdout unbuffered */
  /*  if (batchmode || rangeoutput || insecure || !isatty(fileno(stdout))) */
  if (setvbuf(stdout, NULL, _IONBF, 0) != 0)
    perror("unbuffer stdout");

  yydebug = parser_debug;

  /* before we start looking for files, let's get the home directory, if we can,
   * so we don't have to detect it in more than one place. */
#ifndef __MINGW__
  home = (char *)getenv("HOME");
#else
  /* under Windows it is %HOMEDRIVE%%HOMEPATH% */
  if ((homepath = (char *)getenv("HOMEPATH")) != NULL && (homedrive = (char *)getenv("HOMEDRIVE")) != NULL )  {
    home = (char *)cl_malloc(256);
    sprintf(home, "%s%s", homedrive, homepath);
  }
#endif
  /* note that either way above, home is NULL if the needed env var(s) were not found. */


  /* read initialization file if specified via -I, or if we are in interactive mode */
  if (cqp_init_file ||
      (!child_process && (!batchmode || batchfd == NULL) && which_app != cqpserver)
      ) {

    /*
     * Read init file specified with -I <file>
     *   if no init file was specified, and we're not in batchmode, child mode, or cqpserver,
     *   looks for ~/.cqprc
     * Same with macro init file (-M <file> or ~/.cqpmacros), but ONLY if macros are enabled.
     */

    /*
     * allow interactive commands during processing of initialization file ???
     * (I don't think this is the case!!)
     */

    init_file_fullname[0] = '\0';

    /* read init file specified with -I , otherwise look for $HOME/.cqprc */
    if (cqp_init_file)
      sprintf(init_file_fullname, "%s", cqp_init_file);
    else if (home)
      sprintf(init_file_fullname, "%s%c%s", home, SUBDIR_SEPARATOR, CQPRC_NAME);

    if (init_file_fullname[0] != '\0') {
      if ((cqprc = fopen(init_file_fullname, "r")) != NULL) {

        reading_cqprc = 1;        /* not good for very much, really */
        if (!cqp_parse_file(cqprc, 1)) {
          Rprintf("Parse errors while reading %s, exiting.\n",
                  init_file_fullname);
          exit(1);
        }
        reading_cqprc = 0;

        /* fclose(cqprc);  was already closed by cqp_parse_file!! */
      }
      else if (cqp_init_file) {
        Rprintf("Can't read initialization file %s\n",
                init_file_fullname);
        exit(1);
      }
    }
  }

  if (!enable_macros && macro_init_file)
    cqpmessage(Warning, "Macros not enabled. Ignoring macro init file %s.", macro_init_file);

  if (enable_macros &&
      (macro_init_file ||
       (!child_process && (!batchmode || (batchfd == NULL)) && !(which_app == cqpserver))
       )
      ) {

    init_file_fullname[0] = '\0';

    /* read macro init file specified with -M , otherwise look for ~/.cqpmacros */
    if (macro_init_file)
      sprintf(init_file_fullname, "%s", macro_init_file);
    else if (home)
      sprintf(init_file_fullname, "%s%c%s", home, SUBDIR_SEPARATOR, CQPMACRORC_NAME);

    if (init_file_fullname[0] != '\0') {
      if ((cqprc = fopen(init_file_fullname, "r")) != NULL) {

        reading_cqprc = 1;        /* not good for very much, really */
        if (!cqp_parse_file(cqprc, 1)) {
          Rprintf("Parse errors while reading %s, exiting.\n",
                  init_file_fullname);
          exit(1);
        }
        reading_cqprc = 0;

        /* fclose(cqprc);  was already closed by cqp_parse_file!! */
      }
      else if (macro_init_file) {
        Rprintf("Can't read macro initialization file %s\n",
                init_file_fullname);
        exit(1);
      }
    }
  } /* ends if (!child_process || (batchmode ... ) ... ) */

  check_available_corpora(UNDEF);

  /* load the default corpus. */
  if ((default_corpus) && !set_current_corpus_name(default_corpus, 0)) {
    Rprintf("Can't set current corpus to default corpus %s, exiting.\n",
            default_corpus);
    exit(1);
  }

#ifndef __MINGW__
  if (signal(SIGPIPE, SIG_IGN) == SIG_IGN) {
    /* Rprintf("Couldn't install SIG_IGN for SIGPIPE signal\n"); */
    /* -- be silent about not being able to ignore the SIGPIPE signal, which often happens in slave mode */
    /* note that SIGPIPE does not seem to exist in signal.h under MinGW */
    signal(SIGPIPE, SIG_DFL);
  }
#endif

#ifdef __MINGW__
  /* due to how the home path was calculated, home contains a malloc'ed string */
  cl_free(home);
#endif

  return 1;
}


/**
 * Parses a stream for CQP query syntax.
 *
 * Note that cqp_parse_file() fclose()s fd unless it is STDOUT.
 *
 * @param fd                    File handle of the file to parse.
 * @param exit_on_parse_errors  Boolean: should CQP exit on parse errors?
 * @return                      Boolean: true = all ok, false = a problem.
 */
int
cqp_parse_file(FILE *fd, int exit_on_parse_errors)
{
  int ok, quiet;
  int cqp_status;

  ok = 1;
  quiet = silent || (fd != stdin);

  if (cqp_file_p < MAXCQPFILES) {

    cqp_files[cqp_file_p] = yyin;
    cqp_file_p++;

    yyin = fd;
    yyrestart(yyin);

    /* main read-from-parser loop */
    while (ok && !feof(fd) && !exit_cqp) {
      if (child_process && ferror(fd)) {
        /* in child mode, abort on read errors (to avoid hang-up when parent has died etc.) */
        Rprintf("READ ERROR -- aborting CQP session\n");
        break;
      }

      if (!quiet) {
        if (current_corpus != NULL)
          if (STREQ(current_corpus->name, current_corpus->mother_name))
            Rprintf("%s> ", current_corpus->name);
          else
            Rprintf("%s:%s[%d]> ",
                   current_corpus->mother_name,
                   current_corpus->name,
                   current_corpus->size);
        else
          Rprintf("[no corpus]> ");
      }

      cqp_status = yyparse();
      if ((cqp_status != 0) && exit_on_parse_errors)
        ok = 0;

      if (child_process && (cqp_status != 0) && !reading_cqprc) {
        Rprintf("PARSE ERROR\n"); /*  */
      }
      if (child_process && !reading_cqprc) {
#if 0
        /* empty lines after commands in child mode have been disabled as of version 2.2.b94 */
        Rprintf("\n");                /* print empty line as separator in child mode */
#endif
        fflush(stdout);
        fflush(stderr);
      }

    } /* endwhile */

    if (fd != stdin)
      fclose(fd);

    cqp_file_p--;
    yyin = cqp_files[cqp_file_p];

    if (save_on_exit)
      save_unsaved_subcorpora();

    return ok;
  } /* endif (cqp_file_p < MAXCQPFILES) */
  else {
    Rprintf("CQP: too many nested files (%d)\n", cqp_file_p);
    return 0;
  }
}

/**
 * Parses a string for CQP query syntax.
 *
 * @param s  The string to parse.
 * @return   Boolean: true = all ok, false = a problem.
 */
int
cqp_parse_string(char *s)
{
  int ok, len, abort;
  int cqp_status;

  ok = 1;
  abort = 0;
  len = strlen(s);

  cqp_input_string_position = 0;
  cqp_input_string = s;

  while (ok && (cqp_input_string_position < len) && !exit_cqp) {

    if (abort) {        /* trying to parse a second command -> abort with error */
      cqpmessage(Error, "Multiple commands on a single line not allowed in CQPserver mode.");
      ok = 0;
      break;
    }

    cqp_status = yyparse();
    if (cqp_status != 0)
      ok = 0;

    if (which_app == cqpserver)
      abort = 1;        /* only one command per line in CQPserver (security reasons) */

  } /* endwhile */

  cqp_input_string_position = 0;
  cqp_input_string = NULL;

  return ok;
}






/* ============================================================ */

/** Pointer for the interrupt callback function. */
InterruptCheckProc interruptCallbackHook = NULL;

/**
 * Sets the interrupt callback function.
 *
 * @param f  Pointer to the function to set as interrupt callback.
 * @return   Always 1.
 */
int
setInterruptCallback(InterruptCheckProc f)
{
  interruptCallbackHook = f;
  return 1;
}

/**
 * Calls the interrupt callback function, if set.
 */
void
CheckForInterrupts(void)
{
  if (interruptCallbackHook)
    interruptCallbackHook();
}

/* ============================================================ */
