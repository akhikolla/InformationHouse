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

void Rprintf(const char *, ...);

#include <sys/stat.h>
#include <fcntl.h>

#include <glib.h>

#include <signal.h> /* added by Andreas Blaette  */
#include <sys/socket.h> /* added by Andreas Blaette */

#include "globals.h"
#include "fileutils.h"


/**
 * Gets the size of the specified file; returns EOF for error.
 *
 * @param filename  The file to size up.
 * @return          Size of file in bytes (or EOF if call to stat() failed)
 */
off_t
file_length(char *filename)
{
  struct stat stat_buf;
  if(stat(filename, &stat_buf) == EOF)
    return(EOF);
  else
    return(stat_buf.st_size);
}

/**
 * Gets the size of the specified file; returns EOF for error.
 *
 * As file_length, but the file is specified by file handle, not name.
 *
 * @param fd  The file to size up.
 * @return    Size of file in bytes.
 */
off_t
fd_file_length(FILE *fd)
{
  struct stat stat_buf;
  if (fstat(fileno(fd), &stat_buf) == EOF) return(EOF);
  else return(stat_buf.st_size);
}

/**
 * Gets the size of the specified file; returns EOF for error.
 *
 * As file_length, but the file is specified by number, not name.
 *
 * @see file_length
 * @param fileno  The file to size up.
 * @return        Size of file in bytes.
 */
off_t
fi_file_length(int fileno)
{
  struct stat stat_buf;
  if(fstat(fileno, &stat_buf) == EOF)
    return(EOF);
  else
    return(stat_buf.st_size);
}

/**
 * Gets the size of the specified file; returns EOF for error.
 *
 * Duplicates functionality of file_length, but return is long
 * instead of off_t.
 *
 * @see file_length
 * @param fname  The file to size up.
 * @return       Size of file in bytes.
 */
long
fprobe(char *fname)
{
  struct stat stat_buf;
  
  if(stat(fname, &stat_buf) == EOF) {
    return (long) EOF;
  }
  else
    /* stat_buf->st_mode holds the permission */
    /* we return the file size */
    return stat_buf.st_size;
}


/**
 * Checks whether the specified path indicates a directory.
 *
 * @param path  Path to check.
 * @return      Boolean. (Also false if there's an error.)
 */
int
is_directory(char *path)
{
  struct stat sBuf;

  if (stat(path, &sBuf) < 0) {
    return 0;
  }
  else {
    return S_ISDIR(sBuf.st_mode) ? 1 : 0;
  }
}

/**
 * Checks whether the specified path indicates a regular file.
 *
 * @param path  Path to check.
 * @return      Boolean. (Also false if there's an error.)
 */
int
is_file(char *path)
{
  struct stat sBuf;

  if (stat(path, &sBuf) < 0) {
    return 0;
  }
  else {
    return S_ISREG(sBuf.st_mode) ? 1 : 0;
  }
}

/**
 * Checks whether the specified path indicates a link.
 *
 * Note this function always returns false in Windows, because Windows
 * doesn't have Unix-style links. (.lnk files don't count.)
 *
 * @param path  Path to check.
 * @return      Boolean. (Also false if there's an error.)
 */
int
is_link(char *path)
{
#ifndef __MINGW__
  struct stat sBuf;
  
  if (stat(path, &sBuf) < 0) {
    return 0;
  }
  else {
    return S_ISLNK(sBuf.st_mode) ? 1 : 0;
  }
#else
  return 0;
#endif
}

/**
 * Implementation of automagic I/O streams.
 *
 * In order to make streams completely transparent to the caller
 * and return them as standard FILE* objects, the CL keeps a list
 * of all open streams with the necessary metadata information.
 */
CLStream open_streams;

/**
 * SIGPIPE handler and global status variable
 */
int cl_broken_pipe = 0;

static void
cl_handle_sigpipe(int signum) {
#ifndef __MINGW__
  cl_broken_pipe = 1;
  /* Rprintf("Handle broken pipe signal\n"); */

  if (signal(SIGPIPE, cl_handle_sigpipe) == SIG_ERR)
    perror("CL: Can't reinstall SIGPIPE handler (ignored)"); /* Is this still necessary on modern platforms? */
#endif
}

/** check whether stream type involves a pipe */
#define STREAM_IS_PIPE(type) (type == CL_STREAM_PIPE || type == CL_STREAM_GZIP || type == CL_STREAM_BZIP2)

/**
 * Open stream of specified (or guessed) type for reading or writing
 *
 * I/O streams opened with this function must always be closed with cl_close_stream()!
 *
 * @param filename  Filename or shell command
 * @param mode      Open for reading (CL_STREAM_READ) or writing (CL_STREAM_WRITE)
 * @param type      Type of stream (see above), or guess automagically from <filename> (CL_STREAM_MAGIC)
 *
 * @return          Standard C stream, or NULL on error
 */
FILE *
cl_open_stream(const char *filename, int mode, int type) {
  char *point, *mode_spec;
  int l = strlen(filename);
  FILE *handle;
  CLStream stream;
  char command[2 * CL_MAX_FILENAME_LENGTH]; /* may be longer than CL_MAX_FILENAME_LENGTH */
  char expanded_filename[2 * CL_MAX_FILENAME_LENGTH];

  if (l > CL_MAX_FILENAME_LENGTH) {
    Rprintf("CL: filename '%s' too long (limit: %d bytes)\n", filename, CL_MAX_FILENAME_LENGTH);
    cl_errno = CDA_EBUFFER;
    return NULL;
  }

  /* validate read/write/append mode */
  switch (mode) {
  case CL_STREAM_READ:
    mode_spec = "r";
    break;
  case CL_STREAM_WRITE:
    mode_spec = "w";
    break;
  case CL_STREAM_APPEND:
    mode_spec = "a";
    break;
  default:
    Rprintf("CL: invalid I/O stream mode = %d\n", mode);
    cl_errno = CDA_EARGS;
    return NULL;
  }

  /* apply magic */
  if (type == CL_STREAM_MAGIC || type == CL_STREAM_MAGIC_NOPIPE) {
    /* expand ~/ or $HOME/ */
    if ((strncmp(filename, "~/", 2) == 0) || (strncasecmp(filename, "$home/", 6) == 0)) {
      char *home = getenv("HOME");
      if (home && home[0] != '\0') {
        filename = (filename[0] == '~') ? filename + 2 : filename + 6;
        snprintf(expanded_filename, 2 * CL_MAX_FILENAME_LENGTH, "%s/%s", home, filename);
        filename = expanded_filename;
        l = strlen(filename); /* don't forget to update string length */
      }
    }
    /* guess type of stream */
    type = CL_STREAM_FILE; /* default */
    /* "-" = STDIN or STDOUT */
    if (strcmp(filename, "-") == 0) {
      type = CL_STREAM_STDIO;
    }
    else {
      point = (char *) filename + strspn(filename, " \t");
      /* " | ..." = read or write pipe to shell command */
      if (*point == '|') {
        if (type == CL_STREAM_MAGIC_NOPIPE) {
          cl_errno = CDA_EACCESS;
          return NULL;
        }
        type = CL_STREAM_PIPE;
        point++;
        filename = point + strspn(point, " \t");
      }
      /* *.gz = gzip-compressed file */
      else if (l > 3 && strcasecmp(filename + l - 3, ".gz") == 0) {
        type = CL_STREAM_GZIP;
      }
      /* *.bz2 = bzip2-compressed file */
      else if (l > 4 && strcasecmp(filename + l - 4, ".bz2") == 0) {
        type = CL_STREAM_BZIP2;
      }
    }
  }

  /* file access errors are delayed when reading/writing through pipe to gzip or bzip2,
   * so check first that we can access the file in the appropriate mode */
  if (type == CL_STREAM_GZIP || type == CL_STREAM_BZIP2) {
    handle = fopen(filename, mode_spec);
    if (handle == NULL) {
      cl_errno = CDA_EPOSIX;
      return NULL;
    }
    fclose(handle);
    handle = NULL;
  }

  /* open file or pipe */
  switch (type) {
  case CL_STREAM_FILE:
    handle = fopen(filename, mode_spec);
    break;
  case CL_STREAM_GZIP:
    point = g_shell_quote(filename);
    if (mode == CL_STREAM_APPEND) {
      sprintf(command, "gzip >> %s", point);
      mode_spec = "w";
    }
    else if (mode == CL_STREAM_WRITE)
      sprintf(command, "gzip > %s", point);
    else
      sprintf(command, "gzip -cd %s", point);
    handle = popen(command, mode_spec);
    g_free(point);
    break;
  case CL_STREAM_BZIP2:
    point = g_shell_quote(filename);
    if (mode == CL_STREAM_APPEND) {
      sprintf(command, "bzip2 >> %s", point);
      mode_spec = "w";
    }
    else if (mode == CL_STREAM_WRITE)
      sprintf(command, "bzip2 > %s", point);
    else
      sprintf(command, "bzip2 -cd %s", point);
    handle = popen(command, mode_spec);
    g_free(point);
    break;
  case CL_STREAM_PIPE:
    if (mode == CL_STREAM_APPEND) mode_spec = "w";
    handle = popen(filename, mode_spec);
    break;
  case CL_STREAM_STDIO:
    handle = (mode == CL_STREAM_READ) ? stdin : stdout;
    break;
  default:
    Rprintf("CL: invalid I/O stream type = %d\n", type);
    cl_errno = CDA_EARGS;
    return NULL;
  }
  if (handle == NULL) {
    cl_errno = CDA_EPOSIX;
    return NULL;
  }

  /* add to list of managed streams */
  stream = (CLStream) cl_malloc(sizeof(struct _CLStream));
  stream->handle = handle;
  stream->mode = mode;
  stream->type = type;
  stream->next = open_streams;
  open_streams = stream;

  /* install SIGPIPE handler if opening a pipe stream */
#ifndef __MINGW__
  if (STREAM_IS_PIPE(type)) {
    if (signal(SIGPIPE, cl_handle_sigpipe) == SIG_ERR)
      perror("CL: can't install SIGPIPE handler (ignored)");
  }
#endif

  cl_broken_pipe = 0;
  cl_errno = CDA_OK;
  return handle;
}

/**
 * Close I/O stream
 *
 * This function can only be used for FILE* objects opened with cl_open_stream()!
 *
 * @param stream    An I/O stream that has been opened with cl_open_stream()
 *
 * @return          0 on success, otherwise the error code returned by fclose() or pclose()
 */
int
cl_close_stream(FILE *handle) {
  CLStream stream = open_streams, point;
  int res = 0, was_pipe = 0;

  while (stream) {
    if (stream->handle == handle)
      break;
    stream = stream->next;
  }
  if (!stream) {
    Rprintf("CL: attempt to close non-managed I/O stream with cl_close_stream() [ignored]\n");
    return CDA_EATTTYPE;
  }

  /* close stream appropriately, depending on type */
  switch (stream->type) {
  case CL_STREAM_STDIO:
    break;
  case CL_STREAM_FILE:
    res = fclose(stream->handle);
    break;
  case CL_STREAM_GZIP:
  case CL_STREAM_BZIP2:
  case CL_STREAM_PIPE:
    res = pclose(stream->handle);
    was_pipe = 1;
    break;
  default:
    Rprintf("CL: internal error, managed I/O stream has invalid type = %d\n", stream->type);
    exit(1);
  }

  /* remove stream from list */
  if (stream == open_streams) {
    open_streams = stream->next;
  }
  else {
    for (point = open_streams; point->next != stream; point = point->next)
      /* pass */;
    point->next = stream->next;
  }
  cl_free(stream);

  /* if stream was a pipe, check whether we can uninstall the SIGPIPE handler */
#ifndef __MINGW__
  if (was_pipe) {
    int any_pipe = 0;
    for (point = open_streams; point; point = point->next)
      if (STREAM_IS_PIPE(point->type))
        any_pipe = 1;
    if (!any_pipe) {
      /* last pipe stream closed, uninstall SIGPIPE handler */
      if (signal(SIGPIPE, SIG_IGN) == SIG_ERR)
        perror("CL: can't uninstall SIGPIPE handler (ignored)");
    }
  }
#endif

  cl_broken_pipe = 0;
  cl_errno = (res) ? CDA_EPOSIX : CDA_OK;
  return res;
}
