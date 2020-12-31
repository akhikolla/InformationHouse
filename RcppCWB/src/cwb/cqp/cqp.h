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

#ifndef _cqp_h_
#define _cqp_h_

#define CQPRC_NAME ".cqprc"
#define CQPMACRORC_NAME ".cqpmacros"
/** The number of file handles CQP can store in its file-array (ie max number of nested files) @see cqp_parse_file */
#define MAXCQPFILES 20

/** Size of the CQP query buffer. */
#define QUERY_BUFFER_SIZE 2048

#include <stdio.h>

/** DEPRACATED means of storing a Boolean value  */
typedef int Boolean;

/** DEPRACATED macros for Boolean true and false */
#define True 1
/** DEPRACATED macros for Boolean true and false */
#define False 0
/* TODO In CWB 4 we should change both the above, and the more recent "int/1/0" convention, to C11-style bool/true/false */

/**
 * The "corpus yielding command type" type.
 *
 * Each possible value of the enumeration represents a particular "type"
 * of command that may potentially yield a (sub)corpus.
 */
typedef enum _cyctype {
  NoExpression,
  Query,                  /**< A query (yielding a query-result subcorpus) */
  Activation,             /**< A corpus-activation command. */
  SetOperation,
  Assignment
} CYCtype;

/** Global variable indicating type (CYC) of last expression */
extern CYCtype LastExpression;

extern int reading_cqprc;

extern int cqp_error_status;

/* ======================================== Query Buffer Interface */

/* ========== see parser.l:extendQueryBuffer() for details */
/* ========== initialization done in parse_actions.c:prepare_parse() */

extern char QueryBuffer[QUERY_BUFFER_SIZE];
extern int QueryBufferP;
extern int QueryBufferOverflow;

/* ======================================== Other global variables */

extern char *searchstr;         /* needs to be global, unfortunately */
extern int exit_cqp;                   /**< 1 iff exit-command was issued while parsing */

extern char *cqp_input_string;
extern int cqp_input_string_position;

extern int initialize_cqp(int argc, char **argv);

extern int cqp_parse_file(FILE *fd, int exit_on_parse_errors);

extern int cqp_parse_string(char *s);

/* ====================================================================== */

/**
 * Interrupt callback functions are of this type.
 */
typedef void (*InterruptCheckProc)(void);

/**
 * Boolean indicating that an interruptible process is currently running.
 *
 * The process in question is one that may be expected to be non-instantaneous.
 * This variable is turned off by the Ctrl+C interrupt handler.
 *
 * @see sigINT_signal_handler
 */
extern int EvaluationIsRunning;

extern int setInterruptCallback(InterruptCheckProc f);

extern void CheckForInterrupts(void);

extern int signal_handler_is_installed;

extern void install_signal_handler(void);


/* ====================================================================== */

#endif
