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


#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <limits.h>

#include <sys/types.h>
#include <time.h>
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>

#include "../cl/globals.h"
#include "../cl/macros.h"
#include "../cl/storage.h"      /* NwriteInt() & NwriteInts() */
#include "../cl/lexhash.h"
/* byte order conversion functions taken from Corpus Library */
#include "../cl/endian.h"
#include "../cl/attributes.h"   /* for DEFAULT_ATT_NAME */


/* ---------------------------------------------------------------------- */

/** User privileges of new files (octal format) */
#define UMASK              0644

/* TODO the following belongs in the CL. */
/** Default string used as value of P-attributes when a value is missing ie if a tab-delimited field is empty */
#define UNDEF_VALUE "__UNDEF__"
/** Default string containing the characters that can function as field separators */
#define FIELDSEPS  "\t\n"

/** max number of s-attributes; also max number of p-attributes (-> could change this to implementation as a linked list) */
#define MAXRANGES 1024

/** nr of buckets of lexhashes used for checking duplicate errors (undeclared element and attribute names in XML tags) */
#define REP_CHECK_LEXHASH_SIZE 1000

/** Input buffer size. If we have XML tags with attributes, input lines can become pretty long
 * (but there's basically just a single buffer)
 */
#define MAX_INPUT_LINE_LENGTH  65536

/** Normal extension for CWB input text files. (must have exactly 4 characters; .gz/.bz2 may be added to this if the file is compressed.) */
#define DEFAULT_INFILE_EXTENSION ".vrt"

/* implicit knowledge about CL component files naming conventions */
#define STRUC_RNG  "%s" SUBDIR_SEP_STRING "%s.rng"            /**< CL naming convention for S-attribute RNG files */
#define STRUC_AVX  "%s" SUBDIR_SEP_STRING "%s.avx"            /**< CL naming convention for S-attribute AVX (attribute-value index) files */
#define STRUC_AVS  "%s" SUBDIR_SEP_STRING "%s.avs"            /**< CL naming convention for S-attribute AVS (attribute values) files */
#define POS_CORPUS "%s" SUBDIR_SEP_STRING "%s.corpus"         /**< CL naming convention for P-attribute Corpus files */
#define POS_LEX    "%s" SUBDIR_SEP_STRING "%s.lexicon"        /**< CL naming convention for P-attribute Lexicon files */
#define POS_LEXIDX "%s" SUBDIR_SEP_STRING "%s.lexicon.idx"    /**< CL naming convention for P-attribute Lexicon-index files */



/* ---------------------------------------------------------------------- */

/* global variables representing configuration */

char *field_separators = FIELDSEPS;     /**< string containing the characters that can function as field separators */
char *undef_value = UNDEF_VALUE;        /**< string used as value of P-attributes when a value is missing
                                             ie if a tab-delimited field is empty */
int debug = 0;                          /**< debug mode on or off? */
int silent = 0;                         /**< hide messages */
int verbose = 0;                        /**< show progress (this is _not_ the opposite of silent!) */
int xml_aware = 0;                      /**< substitute XML entities in p-attributes & ignore <? and <! lines */
int skip_empty_lines = 0;               /**< skip empty lines when encoding? */
unsigned line = 0;                      /**< corpus position currently being encoded (ie cpos of _next_ token) */
/* unsigned so it doesn't wrap after first 2^31 tokens and we can abort encoding when corpus size is exceeded */
int strip_blanks = 0;                   /**< strip leading and trailing blanks from input and token annotations */
cl_string_list input_files = NULL;      /**< list of input file(s) (-f option(s)) */
int nr_input_files = 0;                 /**< number of input files (length of list after option processing) */
int current_input_file = 0;             /**< index of input file currently being processed */
char *current_input_file_name = NULL;   /**< filename of current input file, for error messages */
FILE *input_fd = NULL;                  /**< file handle for current input file (or pipe) (text mode!) */
unsigned long input_line = 0;           /**< input line number (reset for each new file) for error messages */
char *registry_file = NULL;             /**< if set, auto-generate registry file named {registry_file}, listing declared attributes */
char *directory = NULL;                 /**< corpus data directory (no longer defaults to current directory) */
char *corpus_character_set = "latin1";  /**< character set label that is inserted into the registry file */
CorpusCharset encoding_charset;         /**< a charset object to be generated from corpus_character_set */
int clean_strings = 0;                  /**< clean up input strings by replacing invalid bytes with '?' */

/* ---------------------------------------------------------------------- */

/* cwb-encode encodes S-attributes and P-attributes, so there is an object-type and global array representing each. */

/**
 * Range object: represents an S-attribute being encoded, and holds some
 * information about the currently-being-processed instance of that S-attribute.
 *
 * TODO should probably be called an SAttr or SAttEncoder or something.
 */
typedef struct _Range {
  char *dir;                    /**< directory where this s-attribute is stored */
  char *name;                   /**< name of the s-attribute (range) */

  int in_registry;              /**< with "-R {reg_file}", this is set to 1 when the attribute is written to the registry
                                     (avoid duplicates) */

  int store_values;             /**< flag indicating whether to store values (does _not_ automatically apply to children, see below) */
  int feature_set;              /**< stored values are feature sets => validate and normalise format */
  int null_attribute;           /**< a NULL attribute ignores all corresponding XML tags, without checking structure or annotations */
  int automatic;                /**< automatic attributes are the 'children' used for recursion and element attributes below  */

  FILE *fd;                     /**< fd of rng component */
  FILE *avx;                    /**< fd of avx component (the attribute value index) */
  FILE *avs;                    /**< fd of avs component (the attribute values) */
  int offset;                   /**< string offset for next string (in avs component) */

  cl_lexhash lh;                /**< lexicon hash for attribute values */

  int has_children;             /**< whether attribute values of XML elements are stored in s-attribute 'children' */
  cl_lexhash el_attributes;     /**< maps XML element attribute names to the appropriate s-attribute 'children' (Range *) */
  cl_string_list el_atts_list;  /**< list of declared element attribute names, required by range_close() function */
  cl_lexhash el_undeclared_attributes; /**< remembers undeclared element attributes, so warnings will be issued only once */

  int max_recursion;            /**< maximum auto-recursion level; 0 = no recursion (maximal regions), -1 = assume flat structure */
  int recursion_level;          /**< keeps track of level of embedding when auto-recursion is activated */
  int element_drop_count;       /**< count how many recursive subelements were dropped because of the max_recursion limit */
  struct _Range **recursion_children;   /**< (usually very short) list of s-attribute 'children' for auto-recursion;
                                             use as array; recursion_children[0] points to self! */

  int is_open;                  /**< boolean: whether there is an open structure at the moment */
  int start_pos;                /**< if this->is_open, remember start position of current range */
  char *annot;                  /**< and annotation (if there is one) */

  int num;                      /**< number of current (if this->is_open) or next structure */

} Range;

/** A global array for keeping track of S-attributes being encoded. */
Range ranges[MAXRANGES];
/** @see ranges */
int range_ptr = 0;

/**
 * WAttr object: represents a P-attribute being encoded.
 *
 * TODO should probably be called a PAttr
 */
typedef struct {
  char *name;                   /**< TODO */
  cl_lexhash lh;                /**< String hash object containing the lexicon for the encoded P attrbute */
  int position;                 /**< Byte index of the lexicon file in progress; contains total number of bytes
                                     written so far (== the beginning of the -next- string that is written) */
  int feature_set;              /**< Boolean: is this a feature set attribute? => validate and normalise format */
  FILE *lex_fd;                 /**< file handle of lexicon component */
  FILE *lexidx_fd;              /**< file handle of lexicon index component */
  FILE *corpus_fd;              /**< file handle of corpus component */
} WAttr;

/** A global array for keeping track of P-attributes being encoded. */
WAttr wattrs[MAXRANGES];
/** @see wattrs */
int wattr_ptr = 0;

/* ---------------------------------------------------------------------- */

/**
 * lookup hash for undeclared s-attributes and s-attributes declared with -S that
 * have annotations (which will be ignored), so warnings are issued only once
 */
cl_lexhash undeclared_sattrs = NULL; 

/* ---------------------------------------------------------------------- */

/** name of the currently running program */
char *progname = NULL;


/* ======================================== helper functions */


/**
 * A replacement for the strtok() function which doesn't skip empty fields.
 *
 * @param s      The string to split.
 * @param delim  Delimiters to use in splitting.
 * @return       The next token from the string.
 */
char *
encode_strtok(register char *s, register const char *delim)
{
  register char *spanp;
  register int c, sc;
  char *tok;
  static char *last;
  

  if (s == NULL && (s = last) == NULL)
    return NULL;
  
  c = *s++;

  if (c == 0) {         /* no non-delimiter characters */
    last = NULL;
    return NULL;
  }
  tok = s - 1;
  
  while (1) {
    spanp = (char *)delim;
    do {
      if ((sc = *spanp++) == c) {
        if (c == 0)
          s = NULL;
        else
          s[-1] = 0;
        last = s;
        return (tok);
      }
    } while (sc != 0);
    c = *s++;
  }
  /* NOTREACHED */
  return NULL;
}





/* ======================================== print time */

#include <sys/types.h>
#include <sys/time.h>

/**
 * Prints a message plus the current time to the specified file/stream.
 *
 * @param stream  Stream to print to.
 * @param msg     Message to incorporate into the string that is printed.
 */
void
encode_print_time(FILE *stream, char *msg)
{
  time_t now;

  time(&now);

  fprintf(stream, "%s: %s\n", msg, ctime(&now));
}

/* ======================================== print error message and exit */

/**
 * Prints a usage message and exits the program.
 */
void 
encode_usage(void)
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:  %s -f <file> [options] -d <dir> [attribute declarations]\n", progname);
  fprintf(stderr, "        ... | %s [options] -d <dir> [attribute declarations]\n\n", progname);
  fprintf(stderr, "Reads verticalised text from stdin (or an input file with -f option) and \n");
  fprintf(stderr, "converts it to the CWB binary format. Each TAB-separated column is encoded as a\n");
  fprintf(stderr, "separate p-attribute. The first p-attribute is named \"word\" (unless changed\n");
  fprintf(stderr, "with -p), additional columns must be declared with -P flags. S-attributes can be\n");
  fprintf(stderr, "declared with -S (without annotations) or -V (with annotations) flags. In\n");
  fprintf(stderr, "the input data, they must appear as opening and closing XML tags on separate\n");
  fprintf(stderr, "lines. For each encoded attribute, one or more data files are created in the\n");
  fprintf(stderr, "current directory (or any directory specified with -d). After encoding, use\n");
  fprintf(stderr, "cwb-makeall to create the required index files and frequency lists, then\n");
  fprintf(stderr, "compress them with cwb-huffcode and cwb-compress-rdx (or preferably use the\n");
  fprintf(stderr, "cwb-make program from the CWB/Perl interface).\n\n");
  fprintf(stderr, "NB: If you re-encode an existing corpus, be sure to delete all old data files,\n");
  fprintf(stderr, "in particular the index and any compressed data files, before running\n");
  fprintf(stderr, "cwb-encode!\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Attribute declarations:\n");
  fprintf(stderr, "  -p <att>  change name of default p-attribute from \"word\" to <att>\n");
  fprintf(stderr, "  -p -      no default p-attribute (all must be declared with -P)\n");
  fprintf(stderr, "  -P <att>  declare additional p-attribute <att>\n");
  fprintf(stderr, "     * append / to mark as feature set => values will be validated and\n");
  fprintf(stderr, "       normalised\n");
  fprintf(stderr, "  -S <att>  declare s-attribute <att> without annotations\n");
  fprintf(stderr, "  -V <att>  declare s-attribute <att> with annotations\n");
  fprintf(stderr, "     * append :<n> for automatic renaming of nested regions, :0 to drop nested\n");
  fprintf(stderr, "       regions (highly recommended, otherwise every start tag will begin a new\n");
  fprintf(stderr, "       flat region)\n");
  fprintf(stderr, "     * attribute-value pairs in XML start tags can be auto-split into separate\n");
  fprintf(stderr, "       s-attributes; the relevant attribute names are appended with + signs\n");
  fprintf(stderr, "       (e.g., -S s:0+id+len stores XML tags like <s id=\"abc\" len=42> in\n");
  fprintf(stderr, "       attributes s, s_id and s_len)\n");
  fprintf(stderr, "     * use -V to store original attribute-value pairs as single string in\n");
  fprintf(stderr, "       addition to auto-splitting into individual s-attributes (e.g. \n");
  fprintf(stderr, "       -V s:0+id+len)\n");
  fprintf(stderr, "     * annotations and values of XML tag attributes can be feature sets; append\n");
  fprintf(stderr, "       / to relevant attribute name for format validation and normalisation\n");
  fprintf(stderr, "       (e.g. -S np:2+agr/+head)\n");
  fprintf(stderr, "  -0 <att>  declare null s-attribute <att> (discards tags)\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -d <dir>  directory for data files created by cwb-encode\n");
  fprintf(stderr, "     * this option always has to be specified (use -d . for current directory)\n");
  fprintf(stderr, "  -f <file> read input from <file> [default is stdin; -f may be used repeatedly]\n");
  fprintf(stderr, "     * gzipped files named *.gz will be decompressed automatically\n");
  fprintf(stderr, "     * alias -t <file> is provided for backward compatibility\n");
  fprintf(stderr, "  -F <dir>  read all files named *" DEFAULT_INFILE_EXTENSION
                               " or *" DEFAULT_INFILE_EXTENSION ".gz in directory <dir>\n");
  fprintf(stderr, "     * files will be added to the corpus in alphabetical order (ASCII)\n");
  fprintf(stderr, "     * it is not possible to scan subdirectories recursively\n");
/* uncomment the following lines when (if...) the -C and -r flags are implemented */
/* a different character for the "C" option would be needed as we are now using it for "clean" */
/*    fprintf(stderr, "  -C <id>   (re-)encode corpus <id> (using data path from registry)\n"); */
/*    fprintf(stderr, "  -r <dir>  set registry directory (for -C flag)\n"); */
  fprintf(stderr, "  -R <rf>   create registry entry (named <rf>) listing all encoded attributes\n");
  fprintf(stderr, "  -B        strip leading/trailing blanks from (input lines & token annotations)\n");
  fprintf(stderr, "  -x        XML-aware (replace XML entities and ignore <!.. and <?..)\n");
  fprintf(stderr, "  -s        skip empty lines in input data (recommended)\n");
  fprintf(stderr, "  -U <str>  insert <str> for missing columns [default: \"%s\"]\n", undef_value);
  fprintf(stderr, "  -b <n>    number of buckets in lexicon hash tables\n");
  fprintf(stderr, "  -c <charset> specify corpus character set (instead of the default latin1)\n");
  fprintf(stderr, "     * valid charsets: ascii ; latin1 .. latin9 ; arabic, greek, hebrew, cyrillic ; utf8\n");
  fprintf(stderr, "     * iso-8859-1 .. iso-8859-15 are also accepted, but converted to canonical names above\n");
  fprintf(stderr, "  -C        clean strings, replacing invalid bytes with '?' (not in UTF-8 mode)\n");
  fprintf(stderr, "  -v        verbose mode (show progress messages while encoding)\n");
  fprintf(stderr, "  -q        quiet mode (suppresses most warnings)\n");
  fprintf(stderr, "  -D        debug mode (quiet, sorry, quite the opposite :-)\n");
  fprintf(stderr, "  -h        this help page\n\n");
  fprintf(stderr, "Part of the IMS Open Corpus Workbench v" VERSION "\n\n");
  exit(2);
}

/**
 * Prints the input line number (and input filename, if applicable) on STDERR,
 * for error messages and warnings.
 */
void
encode_print_input_lineno(void)
{
  if (nr_input_files > 0 && current_input_file_name != NULL)
    fprintf(stderr, "file %s, line #%ld", current_input_file_name, input_line);
  else
    fprintf(stderr, "input line #%ld", input_line);
}

/**
 * Prints an error message to STDERR, automatically adding a
 * message on the location of the error in the corpus.
 *
 * Then exits the program.
 *
 * @param format  Format-specifying string of the error message.
 * @param ...     Additional arguments, printf-style.
 */
void
encode_error(char *format, ...)
{
  va_list ap;
  va_start(ap, format);

  if (format != NULL) {
    vfprintf(stderr, format, ap);
    fprintf(stderr, "\n");
  }
  else {
    fprintf(stderr, "Internal error. Aborted.\n");
  }
  if ((input_line > 0) || (current_input_file > 0)) {
    /* show location only if we've already been reading input */
    fprintf(stderr, "[location of error: ");
    encode_print_input_lineno();
    fprintf(stderr, "]\n");
  }
  exit(1);
}

/* =================================================== processing directories of input files */

/**
 * Get a list of files in a given directory.
 *
 * This function only lists files with .vrt or .vrt.gz extensions,
 * and only files identified  by POSIX stat() as "regular".
 *
 * (Note that .vrt is dependent on DEFAULT_INFILE_EXTENSION.)
 *
 * @see        DEFAULT_INFILE_EXTENSION
 * @param dir  Path of directory to look in.
 * @return     List of paths to files (*including* the directory name).
 *             Returned as a cl_string_list object.
 */
cl_string_list
encode_scan_directory(char *dir)
{
  DIR *dirp;
  struct dirent *dp;
  struct stat statbuf;
  int n_files = 0;
  int len_dir = strlen(dir);
  cl_string_list input_files = cl_new_string_list();
  

  dirp = opendir(dir);
  if (dirp == NULL) {
    perror("Can't access directory");
    encode_error("Failed to scan directory specified with -F %s -- aborted.\n", dir);
  }
  
  errno = 0;
  for (dp = readdir(dirp); dp != NULL; dp = readdir(dirp)) {
    char *name = dp->d_name;
    if (name != NULL) {
      int len_name = strlen(name);
      if ( (len_name >= 5 && (0 == strcasecmp(name + len_name - 4, DEFAULT_INFILE_EXTENSION)))
           || (len_name >= 8 && (0 == strcasecmp(name + len_name - 7, DEFAULT_INFILE_EXTENSION ".gz")))
           || (len_name >= 9 && (0 == strcasecmp(name + len_name - 8, DEFAULT_INFILE_EXTENSION ".bz2"))) )
      {
        char *full_name = (char *) cl_malloc(len_dir + len_name + 2);
        sprintf(full_name, "%s%c%s", dir, SUBDIR_SEPARATOR, name);
        if (stat(full_name, &statbuf) != 0) {
          perror("Can't stat file:");
          encode_error("Failed to access input file %s -- aborted.\n", full_name);
        }
        if (S_ISREG(statbuf.st_mode)) {
          cl_string_list_append(input_files, full_name);
          n_files++;
        }
        else {
          cl_free(full_name);
        }
      }
    }
  }

  if (errno != 0) {
    perror("Error reading directory");
    encode_error("Failed to scan directory specified with -F %s -- aborted.\n", dir);
  }
  if (n_files == 0) {
    fprintf(stderr, "Warning: No input files found in directory -F %s !!\n", dir);
  }
  closedir(dirp);
  
  cl_string_list_qsort(input_files);
  return(input_files);
}

/* =================================================== handling s-attributes and p-attributes */

/**
 * Gets the index (in the ranges array) of a specified S-attribute.
 *
 * @see         ranges
 * @param name  The S-attribute to search for.
 * @return      Index (as integer). -1 if the S-attribute is not found.
 */
int 
range_find(char *name)
{
  int i;

  for (i = 0; i < range_ptr; i++)
    if (strcmp(ranges[i].name, name) == 0)
      return i;
  return -1;
}


/**
 * Prints registry lines for a given s-attribute, and its children,
 * if any, to the specified file handle.
 *
 * @param rng            The s-attribute in question.
 * @param fd             Stream for the registry file to write the line to.
 * @param print_comment  Boolean: if true, a comment on the original XML tags is printed.
 */
void
range_print_registry_line(Range *rng, FILE *fd, int print_comment)
{
  Range *child;
  int i, n_atts;

  if (rng->in_registry)
    return;
  else 
    rng->in_registry = 1;               /* make a note that we've already handled the range */

  if (! rng->null_attribute) {

    if (print_comment) {
      /* print comment showing corresponding XML tags */
      fprintf(fd, "# <%s", rng->name);
      if (rng->has_children) {  /* if there are element attributes, show them in the order of declaration */
        n_atts = cl_string_list_size(rng->el_atts_list);
        for (i = 0; i < n_atts; i++) {
          fprintf(fd, " %s=\"..\"", cl_string_list_get(rng->el_atts_list, i));
        }
      }
      fprintf(fd, "> ... </%s>\n", rng->name);
      /* print comment showing hierarchical structure (if not flat) */
      if (rng->max_recursion == 0) {
        fprintf(fd, "# (no recursive embedding allowed)\n");
      }
      else if (rng->max_recursion > 0) {
        n_atts = rng->max_recursion;
        fprintf(fd, "# (%d levels of embedding: <%s>", n_atts, rng->name);
        for (i = 1; i <= n_atts; i++)
          fprintf(fd, ", <%s>", rng->recursion_children[i]->name);
        fprintf(fd, ").\n");
      }
    }

    /* print registry line for this s-attribute */
    if (rng->store_values)
      fprintf(fd, "STRUCTURE %-20s # [annotations]\n", rng->name);
    else
      fprintf(fd, "STRUCTURE %s\n", rng->name);
    
    /* print recursion children, then element attribute children */
    if (rng->max_recursion > 0) {
      n_atts = rng->max_recursion;
      for (i = 1; i <= n_atts; i++)
        range_print_registry_line(rng->recursion_children[i], fd, 0);
    }
    if (rng->has_children) {    /* element attribute children will print their recursion children as well */
      n_atts = cl_string_list_size(rng->el_atts_list);
      for (i = 0; i < n_atts; i++) {
        cl_lexhash_entry entry = cl_lexhash_find(rng->el_attributes,
                                                 cl_string_list_get(rng->el_atts_list, i));
        child = (Range *) entry->data.pointer;
        range_print_registry_line(child, fd, 0);
      }
    }
    
    if (print_comment)          /* print blank line after each att. declaration block headed by comment */
      fprintf(fd, "\n");
  }
}


/**
 * Creates a Range object to store a specified s-attribute
 * (and, if appropriate, does the same for children-attributes).
 *
 * The new Range object is placed in a global variable, but a pointer
 * is also returned. So you can ignore the return value or not, as
 * you prefer.
 *
 * This is the function where the command-line formalism for defining
 * s-attributes is defined.
 *
 * @see                   ranges
 *
 * @param name            The string from the user specifying the name of
 *                        this attribute, recursion and any "attributes"
 *                        of this XML element - e.g. "text:0+id"
 * @param directory       The directory where the CWB data files will go.
 * @param store_values    boolean: indicates whether this s-attribute was
 *                        specified with -V (true) or -S (false) when the
 *                        program was invoked.
 * @param null_attribute  boolean: this is a null attribute, i.e. an XML
 *                        element to be ignored.
 * @return                Pointer to the new Range object (which is a member
 *                        of the global ranges array).
 */
Range *
range_declare(char *name, char *directory, int store_values, int null_attribute)
{
  char buf[CL_MAX_LINE_LENGTH];
  Range *rng;
  char *p, *rec, *ea_start, *ea;
  cl_lexhash_entry entry;
  int i, is_feature_set;
  char *flag_SV = (store_values) ? "-V" : "-S";

  if (debug)
    fprintf(stderr, "ATT: %s %s\n", flag_SV, name);

  if (range_ptr >= MAXRANGES) {
    encode_error("Too many s-attributes declared (last was <%s>).", name);
  }

  if (directory == NULL)
    encode_error("Error: you must specify a directory for CWB data files with the -d option");

  rng = &ranges[range_ptr];     /* fill next entry in ranges[] */
  range_ptr++;                  /* must increment range pointer now, in case we have children */

  cl_strcpy(buf, name);
  /* check if recursion and/or element attributes are declared */
  if ((rec = strchr(buf, ':')) != NULL) {    /* recursion declaration ":<n>"  */
    *(rec++) = '\0';
    if (strchr(buf, '+'))       /* make sure recursion is declared _before_ element attributes */
      encode_error("Usage error: recursion depth must be declared before element attributes in %s %s !", flag_SV, name);
  }
  p = (rec != NULL) ? rec : buf; /* start looking for element attribute declarations from here */
  if ((ea_start = strchr(p, '+')) != NULL) { /* element att. declaration "+<ea>" */
    *(ea_start++) = '\0';
  }
  if (buf[strlen(buf)-1] == '/') {
    is_feature_set = 1;
    buf[strlen(buf)-1] = '\0';
    if (! store_values)
      encode_error("Usage error: feature set marker '/' is meaningless with -S flag in %s %s !", flag_SV, name);
    if (ea_start != NULL)
      encode_error("Usage error: values of s-attribute %s cannot be feature sets if element attributes are declared (%s %s).",
                   buf, flag_SV, name);
  }
  else {
    is_feature_set = 0;
  }
  /* now buf points to <name> rec points to <n> and ea_start to <ea> of the first element att.;
     all strings are NUL-terminated (ea_start has the form "<ea1>+<ea2>+...+<ea_n>" */

  rng->name = cl_strdup(buf);   /* name of the s-attribute */
  rng->dir = cl_strdup(directory);
  rng->in_registry = 0;
  rng->store_values = store_values;
  rng->feature_set = is_feature_set;
  rng->max_recursion = (rec) ? atoi(rec) : -1; /* set recursion depth: -1 = flat structure */
  rng->recursion_level = 0;
  rng->automatic = 0;

  if (null_attribute) {
    rng->null_attribute = 1;
    if (rec != NULL || ea_start != NULL)
      fprintf(stderr, "Warning: recursion and element attribute specificiers are ignored for null attributes (-0 %s).'n", name);
    return rng;                 /* stop initialisation here; other functions shouldn't do anything with this att */
  }
  else {
    rng->null_attribute = 0;
  }

  if (ea_start != NULL)
    ea_start = cl_strdup(ea_start); /* now buf can be re-used for pathnames below */

  /* open data files for this s-attribute (children will be added later) */
  /* create .rng component */
  sprintf(buf, STRUC_RNG, directory, rng->name);
  if ((rng->fd = fopen(buf, "wb")) == NULL) {
    perror(buf);
    encode_error("Can't write .rng file for s-attribute <%s>.", name);
  }
  if (rng->store_values) {
    /* create .avx and .avs components and initialise lexicon hash */
    sprintf(buf, STRUC_AVS, rng->dir, rng->name);
    if ((rng->avs = fopen(buf, "wb")) == NULL) {
      perror(buf);
      encode_error("Can't write .avs file for s-attribute <%s>.", name);
    }

    sprintf(buf, STRUC_AVX, rng->dir, rng->name);
    if ((rng->avx = fopen(buf, "wb")) == NULL) {
      perror(buf);
      encode_error("Can't write .avx file for s-attribute <%s>.", name);
    }

    rng->lh = cl_new_lexhash(10000);    /* typically, will only have moderate number of entries -> save memory */
  }
  else {
    rng->avs = NULL;
    rng->avx = NULL;
    rng->lh = NULL;
  }
  rng->offset = 0;
  rng->is_open = 0;
  rng->start_pos = 0;
  rng->annot = NULL;
  rng->num = 0;

  /* now that the range is initialised, declare its 'children' if necessary */
  /* recursion children */
  if (rng->max_recursion >= 0) {
    rng->recursion_children = (Range **) cl_calloc(rng->max_recursion + 1, sizeof(Range *));
    rng->recursion_children[0] = rng; /* zeroeth recursion level is stored in the att. itself */
    for (i = 1; i <= rng->max_recursion; i++) {
      /* recursion children have 'flat' structure, because recursion is handled explicitly */
      sprintf(buf, "%s%d%s", rng->name, i, is_feature_set ? "/" : "");
      rng->recursion_children[i] = range_declare(buf, rng->dir, rng->store_values, /*null*/ 0);
      rng->recursion_children[i]->automatic = 1; /* mark as automatically handled attribute */
    }
    rng->recursion_level = 0;
    rng->element_drop_count = 0;
  }
  /* element attributes children can handle recursion on their own */
  if (ea_start != NULL) {
    Range *att_ptr;

    rng->has_children = 1;
    rng->el_attributes = cl_new_lexhash(REP_CHECK_LEXHASH_SIZE);
    rng->el_atts_list = cl_new_string_list();
    rng->el_undeclared_attributes = cl_new_lexhash(REP_CHECK_LEXHASH_SIZE);
    ea = ea_start;
    while (ea != NULL) {
      if ((p = strchr(ea, '+')) != NULL)
        *p = '\0';              /* ea now points to NUL-terminated "<ea_i>" */

      if (rng->max_recursion >= 0) 
        sprintf(buf, "%s_%s:%d", rng->name, ea, rng->max_recursion);
      else
        sprintf(buf, "%s_%s", rng->name, ea);
      /* potential feature set marker (/) is passed on to the respective child attribute and handled there */

      if (ea[strlen(ea)-1] == '/')
        ea[strlen(ea)-1] = '\0';  /* remove feature set marker from element attribute name (used for lookup in encoding) */

      if (cl_lexhash_id(rng->el_attributes, ea) >= 0) 
        encode_error("Element attribute <%s %s=...> declared twice!", rng->name, ea);
      entry = cl_lexhash_add(rng->el_attributes, ea);
      att_ptr = range_declare(buf, rng->dir, 1, /*null*/ 0); /* element att. children always store value, of course */
      att_ptr->automatic = 1;   /* mark as automatically handled attribute */
      entry->data.pointer = att_ptr;
      cl_string_list_append(rng->el_atts_list, cl_strdup(ea)); /* make copy of name (for code cleanness) */

      if (p != NULL)
        ea = p + 1 ;
      else
        ea = NULL;              /* end of element att declarations */
    }
    cl_free(ea_start);          /* don't forget to free copy of element att declaration */
  }
  else {
    rng->has_children = 0;
  }

  return rng;
}

/**
 * Closes a currently open instance of an S-attribute.
 *
 * @param rng      Pointer to the S-attribute to close.
 * @param end_pos  The corpus position at which this instance closes.
 */
void
range_close(Range *rng, int end_pos)
{
  cl_lexhash_entry entry;
  int close_this_range = 0;     /* whether we actually have to close this range (may be skipped or delegated in recursion mode) */
  int i, n_children, l;

  if (rng->null_attribute)      /* do nothing for NULL attributes */
    return;

  if (rng->max_recursion >= 0) {                  /* recursive XML structure */
    rng->recursion_level--;     /* decrement level of nesting */
    if (rng->recursion_level < 0) {
      /* extra close tag (ignored) */
      if (!silent) {
        fprintf(stderr, "Surplus </%s> tag ignored (", rng->name);
        encode_print_input_lineno();
        fprintf(stderr, ").\n");
      }
      rng->recursion_level = 0;
    }
    else if (rng->recursion_level > rng->max_recursion) {
      /* deeply nested ranges are ignored */
    }
    else if (rng->recursion_level > 0) { 
      /* delegated to appropriate recursion child */
      range_close(rng->recursion_children[rng->recursion_level], end_pos);
    }
    else {                      /* rng->recursion_level == 0, i.e. actually applies to rng */
      close_this_range = 1;
    }
  }
  else {                        /* flat structure (traditional mode) */
    if (rng->is_open) {
      close_this_range = 1;     /* ok */
    }
    else {
      /* extra close tag (ignored) */
      if (!silent) {
        fprintf(stderr, "Close tag </%s> without matching open tag ignored (", rng->name);
        encode_print_input_lineno();
        fprintf(stderr, ").\n");
      }
    }
  }

  /* now close the range and write data to disk if we really have to */
  if (close_this_range) {
    if (end_pos >= rng->start_pos) {
      NwriteInt(rng->start_pos, rng->fd); /* write (start, end) to .rng component */
      NwriteInt(end_pos, rng->fd);
      if (rng->store_values) {
        if (rng->annot == NULL) /* shouldn't happen, just to be on the safe side ... */
          rng->annot = cl_strdup("");
        /* check annotation length & truncate if necessary */
        l = strlen(rng->annot);
        if (l >= CL_MAX_LINE_LENGTH) {
          if (!silent) {
            fprintf(stderr, "Value of <%s> region exceeds maximum string length (%d > %d chars), truncated (", 
                    rng->name, l, CL_MAX_LINE_LENGTH-1);
            encode_print_input_lineno();
            fprintf(stderr, ").\n");
          }
          rng->annot[CL_MAX_LINE_LENGTH-2] = '$'; /* truncation marker, as e.g. in Emacs */
          rng->annot[CL_MAX_LINE_LENGTH-1] = '\0';
        }
        /* check if annot is already in hash */
        if ((entry = cl_lexhash_find(rng->lh, rng->annot)) != NULL) {
          /* annotation is already in hash (and hence, stored in .avs component */
          int offset = entry->data.integer;
          /* write (range_num, offset) to .avx component */
          NwriteInt(rng->num, rng->avx); 
          NwriteInt(offset, rng->avx);
        }
        else {
          /* write annotation string to .avs component (at offset rng->offset) */
          fprintf(rng->avs, "%s%c", rng->annot, '\0'); 
          /* write (range_num, current offset) to .avx component */
          NwriteInt(rng->num, rng->avx);  /* this was intended for 'sparse' annotations, which I don't like (so they're no longer there) */
          NwriteInt(rng->offset, rng->avx);
          /* insert annotation string into lexicon hash (with offset_ptr as data ptr) */
          entry = cl_lexhash_add(rng->lh, rng->annot);
          entry->data.integer = rng->offset;
          /* update offset (string length + null byte) */
          rng->offset += strlen(rng->annot) + 1;
          /* check for integer overflow */
          if (rng->offset < 0)
            encode_error("Too many annotation values for <%s> regions (lexicon size > %d bytes)", rng->name, INT_MAX);
        }
        rng->num++;
        cl_free(rng->annot);
      }
      rng->is_open = 0;
    }  /* endif end_pos >= start_pos */
    else {                      
      rng->is_open = 0;      /* silently ignore empty region */
      cl_free(rng->annot);
    }
  }

  /* if range has element attribute children, send corresponding range_close() event to all children in the list
     (recursion and nesting errors will be handled by the children themselves) */
  if (rng->has_children) {
    n_children = cl_string_list_size(rng->el_atts_list);
    for (i = 0; i < n_children; i++) {
      entry = cl_lexhash_find(rng->el_attributes, 
                              cl_string_list_get(rng->el_atts_list, i));
      if (entry == NULL) {
        encode_error("Internal error in <%s>: rng->el_attributes inconsistent with rng->el_atts_list!", rng->name);
      }
      range_close((Range *) entry->data.pointer, end_pos);
    }
  }

  return;
}


/**
 * Opens an instance of the given S-attribute.
 *
 * If rng has element attribute children, range_open() will mess around
 * with the string annotation (otherwise not).
 *
 * @param rng        The S-attribute to open.
 * @param start_pos  The corpus position at which this instance begins.
 * @param annot      The annotation string (the XML element's att-val pairs).
 */
void
range_open(Range *rng, int start_pos, char *annot)
{
  cl_lexhash_entry entry;
  int open_this_range = 0;      /* whether we actually have to open this range (may be skipped or delegated in recursion mode) */
  int i, mark, point, n_children;
  char *el_att_name, *el_att_value;
  char quote_char;              /* quote char used for element attribute value ('"' or '\'') */

  if (rng->null_attribute)      /* do nothing for NULL attributes */
    return;

  if (rng->max_recursion >= 0) {                  /* recursive XML structure */
    if (rng->recursion_level > rng->max_recursion) {
      /* deeply nested ranges are ignored; count how many we've lost */
      rng->element_drop_count++;
    }
    else if (rng->recursion_level > 0) {
      /* delegate to appropriate recursion child (with same annotation) */
      range_open(rng->recursion_children[rng->recursion_level], start_pos,
                 (rng->store_values) ? annot : NULL);
      /* recursion children don't parse the annotation string, so annot will remain untouched;
         since recursion children always have the same -S or -V behaviour as the parent, we only
         pass the annotation string for -V attributes in order to avoid spurious warnings */
    }
    else {                      /* rng->recursion_level == 0, i.e. actually applies to rng */
      open_this_range = 1;
    }
    rng->recursion_level++;     /* increment level of nesting */
  }
  else {                                          /* flat structure (traditional mode) */
    if (rng->is_open) {
      /* if we assume flat structure, implicitly close a range that is already open */
      range_close(rng, start_pos - 1);
    }
    open_this_range = 1;        /* with flat structure, a start tag always opens a new range */
  }

  if (open_this_range) {
    rng->is_open = 1;
    rng->start_pos = line;

    if (annot == NULL)          /* shouldn't happen, but just to be on the safe side ... */
      annot = "";

    if (rng->store_values) {
      rng->annot = cl_strdup(annot); /* remember annotation for range_close(); must strdup because it's pointer into linebuf[] */
      /* don't warn about empty annotations, because that's explicitly allowed! */
      if (strip_blanks) {               /* annotation string may have trailing blanks */
        i = strlen(rng->annot) - 1;
        while ((i >= 0) && ((rng->annot[i] == ' ') || (rng->annot[i] == '\t')))
          rng->annot[i--] = '\0';
      }
      if (rng->feature_set) {
        char *token = cl_make_set(rng->annot, /*split*/ 0);
        if (token == NULL) {
          if (! silent) {
            fprintf(stderr, "Warning: '%s' is not a valid feature set for s-attribute %s, replaced by empty set | (", 
                            rng->annot, rng->name);
            encode_print_input_lineno();
            fprintf(stderr, ")\n");
          }
          token = cl_strdup("|"); /* rng->annot will be free()d later, so it must be an allocated string */
        }
        cl_free(rng->annot);
        rng->annot = token;
      }
    }
    else {
      /* warn about non-empty annotation string in -S attribute (unless annotation string is parsed), but only once */
      if ((!rng->has_children) && (*annot != '\0')) { 
        if (!cl_lexhash_freq(undeclared_sattrs, rng->name)) {
          if (!silent) {
            fprintf(stderr, "Annotations of s-attribute <%s> not stored (", rng->name);
            encode_print_input_lineno();
            fprintf(stderr, ", warning issued only once).\n");
          }
          cl_lexhash_add(undeclared_sattrs, rng->name); /* we can re-use the lookup hash for undeclared s-attributes :o) */
        }
      }
    }
  }
  
  /* if rng has element attribute children, try to parse the annotation string into
     XML attribute="value" or attribute=id pairs (destructively modifying the original)
     NB: there must not be any leading whitespace in annot
     NB: we don't bother about recursion here; the child attributes will take care of that themselves */
  if (rng->has_children) {
    /* we have to make sure that regions are opened for all declared element attributes, and that
       no element attribute occurs more than once in order to ensure proper nesting; */
    n_children = cl_string_list_size(rng->el_atts_list); /* use the integer data field of the el_attributes hash */
    for (i = 0; i < n_children; i++) {
      entry = cl_lexhash_find(rng->el_attributes, 
                              cl_string_list_get(rng->el_atts_list, i));
      entry->data.integer = 0;  /* initialise to 0, i.e. "not handled" */
    }

    mark = 0;                   /* mark and point are offsets into annot[] */
    while (annot[mark] != '\0') {
      point = mark;

      /* identify XML element attribute name (slightly relaxed attribute naming conventions) */
      while (cl_xml_is_name_char(annot[point]))
        point++;
      while ((annot[point] == ' ') || (annot[point] == '\t')) {
        annot[point] = '\0';    /* skip optional whitespace before '=' separator, and remove it from el.att. name */
        point++;
      }

      /* now annot[point] should be the separator '=' char */
      if (annot[point] != '=') {
        if (!silent) {
          fprintf(stderr, "Attributes of open tag <%s ...> ignored because of syntax error (``='' not found) (", rng->name);
          encode_print_input_lineno();
          fprintf(stderr, ").\n");
        }
        break;                  /* stop processing attributes */
      }
      annot[point] = '\0';      /* terminate el. attribute name in el_att_name = (annot+mark) */
      el_att_name = annot + mark;
      mark = point + 1;
      while ((annot[mark] == ' ') || (annot[mark] == '\t'))
        mark++; /* skip optional whitespace after '=' separator */

      /* now get the attribute value (either "value" or 'value' or id) */
      quote_char = annot[mark];
      if ((quote_char == '"') || (quote_char == '\'')) {    /* attribute="value" or attribute='value' format */
        mark++;                 /* assume it's well-formed XML and just look for next occurrence of quote_char */
        point = mark;
        while ((annot[point] != quote_char) && (annot[point] != '\0'))
          point++;
        if (annot[point] == '\0') { /* syntax error: missing end quote */
          if (!silent) {
            fprintf(stderr, "Attributes of open tag <%s ...> ignored because of syntax error (value missing end quote) (", rng->name);
            encode_print_input_lineno();
            fprintf(stderr, ").\n");
          }
          break;                /* stop processing attributes */
        }
        el_att_value = annot + mark;
        annot[point] = '\0';    /* terminate attribute value, and advance mark */
        mark = point + 1;
      }
      else {                    /* attribute=id format (accepts same id's as el.att. name) */
        point = mark;
        while (cl_xml_is_name_char(annot[point]))
          point++;
        el_att_value = annot + mark;
        if (annot[point] == '\0') { /* end of annot[] reached, don't advance mark beyond NUL byte */
          mark = point;
        }
        else {                  /* terminate attribute value, and advance mark */
          annot[point] = '\0';
          mark = point + 1;
        }
        if (strlen(el_att_value) == 0) { /* syntax error: attribute=id with empty value (not allowed) */
          if (!silent) {
            fprintf(stderr, "Attributes of open tag <%s ...> ignored because of syntax error (attribute=id with empty value (not allowed)) (", rng->name);
            encode_print_input_lineno();
            fprintf(stderr, ").\n");
          }
          break;                /* stop processing attributes */
        }
      }

      /* syntax check: el_att_name must be non-empty (values "" and '' are allowed) */
      if (strlen(el_att_name) == 0) {
        if (!silent) {
          fprintf(stderr, "Attributes of open tag <%s ...> ignored because of syntax error (empty attribute name)) (", rng->name);
          encode_print_input_lineno();
          fprintf(stderr, ").\n");
        }
        break;          /* stop processing attributes */
      }

      /* now delegate the attribute/value pair to the appropriate child attribute */
      entry = cl_lexhash_find(rng->el_attributes, el_att_name);
      if (entry == NULL) {      /* undeclared element attribute (ignored) */
        if (!cl_lexhash_freq(rng->el_undeclared_attributes, el_att_name)) {
          if (!silent) {
            fprintf(stderr, "Undeclared element attribute <%s %s=...> ignored (", rng->name, el_att_name);
            encode_print_input_lineno();
            fprintf(stderr, ", warning issued only once).\n");
          }
          cl_lexhash_add(rng->el_undeclared_attributes, el_att_name);
        }
      }
      else {                    /* declared element attribute -> decode XML entities in value and delegate to child */
        if (entry->data.integer) {
          /* attribute already handled, i.e. it must have occurred twice in start tag -> issue warning */
          if (!silent) {
            fprintf(stderr, "Duplicate attribute value <%s %s=... %s=...> ignored (", rng->name, el_att_name, el_att_name);
            encode_print_input_lineno();
            fprintf(stderr, ").\n");
          }
        }
        else {
          entry->data.integer = 1; /* mark el. att. as handled */
          cl_xml_entity_decode(el_att_value);
          range_open((Range *) entry->data.pointer, start_pos, el_att_value);
        }
      }

      while ((annot[mark] == ' ') || (annot[mark] == '\t'))
        mark++;                 /* skip whitespace before next attribute="value" pair */
    }

    /* phew. that was a bit of work; 
       and we still have to make sure that missing element attributes are encoded as empty strings  */
    for (i = 0; i < n_children; i++) {
      entry = cl_lexhash_find(rng->el_attributes, cl_string_list_get(rng->el_atts_list, i));
      if (entry->data.integer == 0) {
        range_open((Range *) entry->data.pointer, start_pos, "");
      }
    }

  } /* end if range has attribute children */

  return;
}

/**
 * Finds a p-attribute (in the global wattrs array).
 *
 * Returns the index (in wattrs) of the P-attribute with the given name.
 *
 * @see         wattrs
 * @param name  The P-attribute to search for.
 * @return      Index (as integer), or -1 if not found.
 */
int 
wattr_find(char *name)
{
  int i;

  for (i = 0; i < wattr_ptr; i++)
    if (strcmp(wattrs[i].name, name) == 0)
      return i;
  return -1;
}


/**
 * Sets up a new p-attribute, including opening corpus, lex and index file handles.
 *
 * Note: corpus_fd is a binary file, lex_fd is a text file(*), and lexidx_fd is
 * a binary file.
 *
 * (*) But lexicon items are delimited by '\0' not by '\n'. Therefore '\n' is never written,
 * so the text/binary distinction doesn't matter much.
 *
 * @param name        Identifier string of the p-attribute
 * @param directory   Directory in which CWB data files are to be created.
 * @param nr_buckets  Number of buckets in the lexhash of the new p-attribute (value passed to cl_new_lexhash() )
 * @return            Always 1.
 */
int 
wattr_declare(char *name, char *directory, int nr_buckets)
{
  char corname[CL_MAX_LINE_LENGTH];
  char lexname[CL_MAX_LINE_LENGTH];
  char idxname[CL_MAX_LINE_LENGTH];

  if (name == NULL)
    name = DEFAULT_ATT_NAME;

  /* TODO why is this a parameter rather than a global ... ? */
  if (directory == NULL)
    encode_error("Error: you must specify a directory for CWB data files with the -d option");
  /* This should be checked in option parsing, NOT here. (Specifically, in the -S and -P options.) */

  /* copy the name supplied as an argument (first removing feature set flag / if needful */
  wattrs[wattr_ptr].name = cl_strdup(name);
  if (name[strlen(name)-1] == '/') {
    wattrs[wattr_ptr].name[strlen(name)-1] = '\0';
    wattrs[wattr_ptr].feature_set = 1;
  }
  else {
    wattrs[wattr_ptr].feature_set = 0;
  }

  wattrs[wattr_ptr].lh = cl_new_lexhash(nr_buckets);
  
  wattrs[wattr_ptr].position = 0;

  /* We now create paths for each of the three files that this encoder generates.
   * The paths aren't stored in the Wattr - only the file handles from opening them. */
  sprintf(corname, POS_CORPUS, directory, wattrs[wattr_ptr].name);
  sprintf(lexname, POS_LEX,    directory, wattrs[wattr_ptr].name);
  sprintf(idxname, POS_LEXIDX, directory, wattrs[wattr_ptr].name);

  if ((wattrs[wattr_ptr].corpus_fd = fopen(corname, "wb")) == NULL) {
    perror(corname);
    encode_error("Can't write .corpus file for %s attribute.", name);
  }

  if ((wattrs[wattr_ptr].lex_fd = fopen(lexname, "w")) == NULL) {
    perror(lexname);
    encode_error("Can't write .lexicon file for %s attribute.", name);
  }

  if ((wattrs[wattr_ptr].lexidx_fd = fopen(idxname, "wb")) == NULL) {
    perror(idxname);
    encode_error("Can't write .lexicon.idx file for %s attribute.", name);
  }

  wattr_ptr++;

  return 1;
}

/**
 * Closes all three file handles for each of the wattr objects
 * in cwb-encode's global array.
 */
void
wattr_close_all(void)
{
  int i;

  for (i = 0; i < wattr_ptr; i++) {
    if (EOF == fclose(wattrs[i].lex_fd)) {
      perror("fclose() failed");
      encode_error("Error writing .lexicon file for %s attribute", wattrs[i].name);
    }
    if (EOF == fclose(wattrs[i].lexidx_fd)) {
      perror("fclose() failed");
      encode_error("Error writing .lexicon.idx file for %s attribute", wattrs[i].name);
    }
    if (EOF == fclose(wattrs[i].corpus_fd)) {
      perror("fclose() failed");
      encode_error("Error writing .corpus file for %s attribute", wattrs[i].name);
    }
  }
}




/**
 * Parses program options and sets global variables.
 *
 * @param argc  argc - passed from main()
 * @param argv  argv - passed from main()
 *
 */
void
encode_parse_options(int argc, char **argv)
{
  int c;
  extern char *optarg;
  extern int optind;
  struct stat dir_status;

  char *prefix = DEFAULT_ATT_NAME;

  int number_of_buckets = 0;    /* -> use CL default unless changed with -b <n> */
  int first_attr_declared = 0;  /* whether we have already declared the default 'word' attribute (useful for "-p -") */

  cl_string_list dir_files;   /* list of input files found in directory (-F option) */
  int i, l;

  while((c = getopt(argc, argv, "p:P:S:V:0:f:t:F:d:R:U:Bsb:c:CxvqhD")) != EOF)
    switch(c) {

      /* -B: strip leading and trailing blanks from tokens and annotations */
    case 'B':
      strip_blanks++;
      break;

      /* -v: show progress messages */
    case 'v':
      verbose++;
      break;

      /* -q: suppress warnings (quiet mode) */
    case 'q':
      silent++;
      break;

      /* -c: specifies a character set */
    case 'c':
      corpus_character_set = cl_charset_name_canonical(optarg);
      if (corpus_character_set == NULL)
        encode_error("Invalid character set specified with the -c flag!");
      break;

      /* -C: clean up strings (remove invalid bytes) */
    case 'C':
      clean_strings++;
      break;

      /* -x: translate XML entities and ignore declarations & comments */
    case 'x':
      xml_aware++;
      break;

      /* -p <att>: change name of first p-attribute ("-p -": skip first attribute) */
    case 'p':
      if (first_attr_declared)
        encode_error("Usage error: -p option used after -P <att>, or used twice.");
      prefix = optarg;
      if (strcmp(prefix, "-") != 0) {
        wattr_declare(prefix, directory, number_of_buckets);
      }
      first_attr_declared = 1;  /* even if we haven't _really_ declared it because it's "-" */
      break;

      /* -d <dir>: create files in this directory */
    case 'd':
      directory = optarg;
      /* Check if directory exists */
      if (stat(directory, &dir_status) != 0 || !(dir_status.st_mode & S_IFDIR)) {
        encode_error("Error: data directory '%s' does not exist.\nPlease create this directory first.",
                     directory);
      }
      break;

      /* -R <rf>: create registry file named <rf> */
    case 'R':
      if (registry_file != NULL)
        encode_error("Usage error: -R option used twice.");
      else {
        int size;
        int registry_is_ok = 1;
        int registry_is_canonical = 1;
        registry_file = optarg;

        /* Check for path ending in slash and for non-lowercase in last part of the filename;
         * allow EITHER possible value of SUBDIR_SEPARATOR */
        size = strlen(registry_file) - 1;
        if ((size < 0) || (registry_file[size] == '/') || (registry_file[size] == '\\'))
          encode_error("Usage error: invalid filename '%s' for registry entry", registry_file);

        while (size >= 0 && registry_file[size] != '/' && registry_file[size] != '\\') {
          char c = registry_file[size];
          if ((c >= 'A' && c <= 'Z') || c == '.' || c == '~')
            registry_is_ok = 0; /* uppercase characters, '.' and '~' are definitely not allowed */
          
          if (!( c == '_' || c == '-' || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') ))
            registry_is_canonical = 0; /* new canonical form allows only ASCII a-z, 0-9, _, - */

          size--;
        }
        
        if (!registry_is_ok)
          encode_error("Usage error: invalid filename '%s' for registry entry.\n"
              "Filename must not contain uppercase letters, '.' or '~'.", registry_file + size + 1);
        if (!registry_is_canonical)
          fprintf(stderr, "Warning: filename '%s' of registry entry not in canonical format.\n"
              "(Allowed characters: a-z, 0-9, -, _)\n", registry_file + size + 1);

        if (size >= 0) {
          /* the registry filename includes a directory part, so check that it exists and is indeed a directory */
          char sep = registry_file[size];
          registry_file[size] = 0; /* now registry_file holds the directory part as a NUL-terminated string */
          if (stat(registry_file, &dir_status) != 0 || !(dir_status.st_mode & S_IFDIR))
            encode_error("Error: registry directory '%s' does not exist.\nPlease create this directory first.", registry_file);
          registry_file[size] = sep;
        }
        
      }
      break;

      /* -f, -t: verticalised text input file */
    case 't':
    case 'f':
      cl_string_list_append(input_files, optarg);
      break;

      /* -F: read all files named *.vrt or *.vrt.gz in directory */
    case 'F':
      dir_files = encode_scan_directory(optarg);
      l = cl_string_list_size(dir_files);
      for (i = 0; i < l; i++) {
        cl_string_list_append(input_files, cl_string_list_get(dir_files, i));
      }
      cl_delete_string_list(dir_files); /* allocated strings have been moved into input_files, so don't free() them */
      break;

      /* -s: skip empty lines */
    case 's': 
      skip_empty_lines++;
      break;

      /* -b <n>: number of buckets */
    case 'b':
      number_of_buckets = atoi(optarg);
      break;

      /* -D: debug mode */
    case 'D':
      debug++;
      break;
      
      /* -S: declare s-attribute without annotations */
    case 'S':
      if (range_ptr < MAXRANGES) {
        if (range_find(optarg) == -1)
          range_declare(optarg, directory, /*annot*/ 0, /*null*/ 0);
        else 
          encode_error("Usage error: s-attribute <%s> declared twice!", optarg);
      }
      else 
        encode_error("Too many s-attributes (max. %d).", MAXRANGES);
      break;

      /* -V: declare s-attribute with annotations */
    case 'V':
      if (range_ptr < MAXRANGES) {
        if (range_find(optarg) == -1) {
          range_declare(optarg, directory, /*annot*/ 1, /*null*/ 0);
        }
        else 
          encode_error("Usage error: s-attribute <%s> declared twice!", optarg);
      }
      else
        encode_error("Too many s-attributes (max. %d).", MAXRANGES);
      break;

      /* -0: declare NULL s-attribute */
    case '0':
      if (range_ptr < MAXRANGES) {
        if (range_find(optarg) == -1) {
          range_declare(optarg, directory, /*annot*/ 0, /*null*/ 1);
        }
        else 
          encode_error("Usage error: s-attribute <%s> declared twice!", optarg);
      }
      else
        encode_error("Too many s-attributes (max. %d).", MAXRANGES);
      break;

      /* -P: declare additional p-attribute */
    case 'P':

      if (!first_attr_declared) { /* no word attribute declared yet */
        wattr_declare(prefix, directory, number_of_buckets);
        first_attr_declared = 1;
      }

      if (wattr_ptr < MAXRANGES) {
        if (wattr_find(optarg) == -1)
          wattr_declare(optarg, directory, number_of_buckets);
        else 
          encode_error("Usage error: %s attribute declared twice!", optarg);
      }
      else
        encode_error("Too many p-attributes (max. %d).", MAXRANGES);
      break;

      /* -U: default value for missing columns */
    case 'U':
      undef_value = optarg;
      break;

      /* default or -h: help page */
    case 'h':
    default:
      encode_usage();
      break;
    }
  
  /* if no attributes have been declared, declare the standard attribute */
  if (!first_attr_declared)     /* no word attribute declared yet */
    wattr_declare(prefix, directory, number_of_buckets);


  /* now, check the default and obligatory values */
  if (optind < argc) {
    
    fprintf(stderr, "%s:\n  Warning: additional arguments in command ignored:",
            progname);

    while (optind < argc) {
      fprintf(stderr, " %s", argv[optind]);
      optind++;
    }
    fprintf(stderr, "\n  (perhaps you forgot -P, -p, -S, or -V before an attribute name?)\n");
    
  }

}




/**
 * Processes a token data line.
 *
 * That is, it processes a line that is *not* an XML line.
 *
 * Note that this is destructive - the argument character
 * string will be changed *in situ* via an strtok-like mechanim.
 *
 * @param str  A string containing the line to process.
 */
void
encode_add_wattr_line(char *str)
{
  /* fc = field counter (current column number, zero indexed)
   * id = container for lexicon ID int.
   * length = temp holder for a strlen return. */
  int fc, id, length;
  /* field = the current column (string).
   * token = token we will store (same as field, except in case of feature sets).
   * Both are pointers to suitable chunks of the parameter string. */
  char *field, *token;

  cl_lexhash_entry entry;

  /* the following tokenization code messes around with the containts of the str parameter,
   * which (in the usage in this program) means changing linebuf[] in main()! */
  for (field = encode_strtok(str, field_separators), fc = 0;
       fc < wattr_ptr; 
       field = encode_strtok(NULL, field_separators), fc++) {
    /* LOOP across each column in the line... */
    
    if ((field != NULL) && strip_blanks) { /* need to strip both leading & trailing blanks from field values */
      length = strlen(field);
      while ((length > 0) && (field[length-1] == ' ')) {
        length--;
        field[length] = '\0';
      }
      while (*field == ' ') 
        field++;
    }
    if ((field != NULL) && (field[0] == '\0'))
      field = NULL;  /* field == NULL -> missing field; field == "" -> empty field; both inserted as __UNDEF__ */

    if ((field != NULL) && xml_aware) 
      cl_xml_entity_decode(field);

    if (field == NULL)          /* mustn't do this before cl_xml_entity_decode(), because undef_value is a constant */
      field = undef_value;

    if (wattrs[fc].feature_set) {
      token = cl_make_set(field, /*split*/ 0);
      if (token == NULL) {
        if (! silent) {
          fprintf(stderr, "Warning: '%s' is not a valid feature set for -P %s/, replaced by empty set | (", 
                          field, wattrs[fc].name);
          encode_print_input_lineno();
          fprintf(stderr, ")\n");
        }
        token = cl_strdup("|");
        /* so we always have to cl_free() token for feature set attributes,
         * because either cl_make_set or cl_strdup was used */
      }
    }
    else {
      token = field;
    }

    /* check annotation length & truncate if necessary (assumes it's ok to modify token[] destructively) */
    length = strlen(token);
    if (length >= CL_MAX_LINE_LENGTH) {
      if (!silent) {
        fprintf(stderr, "Value of p-attribute '%s' exceeds maximum string length (%d > %d chars), truncated (", 
                wattrs[fc].name, length, CL_MAX_LINE_LENGTH-1);
        encode_print_input_lineno();
        fprintf(stderr, ").\n");
      }
      token[CL_MAX_LINE_LENGTH-2] = '$'; /* truncation marker, as e.g. in Emacs */
      token[CL_MAX_LINE_LENGTH-1] = '\0';
    }

    id = cl_lexhash_id(wattrs[fc].lh, token);
    if (id < 0) {
      /* new entry -> write LEXIDX & LEXICON files */
      NwriteInt(wattrs[fc].position, wattrs[fc].lexidx_fd);
      wattrs[fc].position += strlen(token) + 1;
      if (wattrs[fc].position < 0)
        encode_error("Maximum size of .lexicon file exceeded for %s attribute (> %d bytes)", wattrs[fc].name, INT_MAX);
      if (EOF == fputs(token, wattrs[fc].lex_fd)) {
        perror("fputs() write error");
        encode_error("Error writing .lexicon file for %s attribute.", wattrs[fc].name);
      }
      if (EOF == putc('\0', wattrs[fc].lex_fd)) {
        perror("putc() write error");
        encode_error("Error writing .lexicon file for %s attribute.", wattrs[fc].name);
      }
      entry = cl_lexhash_add(wattrs[fc].lh, token);
      id = entry->id;
    }

    if (wattrs[fc].feature_set)
      cl_free(token); /* string has been allocated by cl_make_set(). See above.  */

    NwriteInt(id, wattrs[fc].corpus_fd);
  } /* end for loop (for each column in input data...) */
}
 
/**
 * Reads one input line into the specified buffer
 * (either from stdin, or from one or more input files).
 *
 * The input files are not passed to the function,
 * but are taken from the program global variables.
 *
 * This function returns False when the last input file
 * has been completely read, and automatically closes files.
 *
 * If the line that is read is not valid according to the
 * character set specified for the corpus, then an error
 * will be printed and the program shut down.
 *
 * @param buffer   Where to load the line to. Assumed to be
 *                 MAX_INPUT_LINE_LENGTH long.
 * @param bufsize  Not currently used, but should be
 *                 MAX_INPUT_LINE_LENGTH in case of future use!
 *
 * @return         boolean: true for all OK, false for a problem.
 */
int
encode_get_input_line(char *buffer, int bufsize)
{
  int ok;

  if (nr_input_files == 0) {
    /* read one line of text from stdin */
    ok = (NULL != fgets(buffer, MAX_INPUT_LINE_LENGTH, stdin));
  }
  else {
    /* input_fd is set to NULL in global initialisation. */
    if (! input_fd) {
      if (current_input_file >= nr_input_files)
        return 0;
      
      current_input_file_name = cl_string_list_get(input_files, current_input_file);

      input_fd = cl_open_stream(current_input_file_name, CL_STREAM_READ, CL_STREAM_MAGIC);
      if (input_fd == NULL) {
        cl_error(current_input_file_name);
        encode_error("Can't open input file %s!", current_input_file_name);
      }
      
      input_line = 0;
    } /* endif no input file is open */

    /* read one line of text from current input file */
    ok = (NULL != fgets(buffer, MAX_INPUT_LINE_LENGTH, input_fd));

    if (ok) {
      /* on first line of file, skip UTF8 byte-order-mark if present */
      if (input_line == 0 && encoding_charset == utf8)
        if (buffer[0] == (char)0xEF && buffer[1] == (char)0xBB && buffer[2] == (char)0xBF)
          cl_strcpy(buffer, (buffer+3));
    }
    else {
      /* assume we're at end of file -> close current input file, and try reading from next one */
      ok = (0 == cl_close_stream(input_fd));
      if (! ok) {
        fprintf(stderr, "ERROR reading from file %s (ignored).\n", current_input_file_name);
        cl_error(current_input_file_name);
      }

      /* use recursive tail call to open the next input file and read from it */
      input_fd = NULL;
      current_input_file++;
      return encode_get_input_line(buffer, bufsize);
    }
  } /* end of block to follow if we're not reading from stdin. */

  /* check encoding and standardise Unicode character composition */
  if (!cl_string_validate_encoding(buffer, encoding_charset, clean_strings))
    encode_error("Encoding error: an invalid byte or byte sequence for charset \"%s\" was encountered.\n",
        corpus_character_set);
  /* normalize UTF8 to precomposed form, but don't bother with the redundant function call otherwise */
  if (encoding_charset == utf8)
    cl_string_canonical(buffer, utf8, REQUIRE_NFC, MAX_INPUT_LINE_LENGTH);
  /* finally, get rid of C0 controls iff the user asked us to clean up strings */
  if (clean_strings)
    cl_string_zap_controls(buffer, encoding_charset, '?', 0, 0);
  /* note we DIDN'T zap tab and newline, because this string has yet to be column-split */

  return ok;
}

/**
 * Writes a registry file for the corpus that has been encoded.
 * Part of cwb-encode; not a library function.
 *
 * @param registry_file  String containing the path of the file to write.
 */
void
encode_generate_registry_file(char *registry_file)
{
  FILE *registry_fd;
  char *registry_id;          /* use last part of registry filename (i.e. string following last '/' character) */
  char *corpus_name = NULL;   /* name of the corpus == uppercase version of registry_id */
  char *info_file = NULL;     /* name of INFO file: <dir>/.info or <dir>/corpus-info.txt under win; see cl/globals.h */
  char *path = NULL;
  int i;

  if (debug)
    fprintf(stderr, "Writing registry file %s ...\n", registry_file);

  if ((registry_fd = fopen(registry_file, "w")) == NULL) {
    perror(registry_file);
    encode_error("Can't create registry entry in file %s!", registry_file);
  }

  i = strlen(registry_file) - 1;
  while ((i > 0) && (registry_file[i-1] != SUBDIR_SEPARATOR))
    i--;
  registry_id = registry_file + i;

  if (! cl_id_validate(registry_id))
      encode_error("%s is not a valid corpus ID! Can't create registry entry.", registry_id);
  /* enforce the "lowercase characters only" rule */
  cl_id_tolower(registry_id);

  i = strlen(directory) - 1;
  while ((i > 0) && (directory[i] == SUBDIR_SEPARATOR))
    directory[i--] = '\0';    /* remove trailing '/' from home directory */

  /* copy registry_id and convert it to uppercase */
  corpus_name = cl_strdup(registry_id);
  cl_id_toupper(corpus_name);

  info_file = (char *) cl_malloc(strlen(directory) + 1 + strlen(INFOFILE_DEFAULT_NAME) + 4); /* extra bytes as safety margin */
  sprintf(info_file, "%s%c%s", directory, SUBDIR_SEPARATOR, INFOFILE_DEFAULT_NAME);

  /* write header part for registry file */
  fprintf(registry_fd, "##\n## registry entry for corpus %s\n##\n\n", corpus_name);
  fprintf(registry_fd, "# long descriptive name for the corpus\n");
  fprintf(registry_fd, "NAME \"\"\n");
  fprintf(registry_fd, "# corpus ID (must be lowercase in registry!)\n");
  fprintf(registry_fd, "ID   %s\n", registry_id);
  fprintf(registry_fd, "# path to binary data files\n");
  path = cl_path_registry_quote(directory);
  fprintf(registry_fd, "HOME %s\n", path);
  cl_free(path);
  fprintf(registry_fd, "# optional info file (displayed by \"info;\" command in CQP)\n");
  path = cl_path_registry_quote(info_file);
  fprintf(registry_fd, "INFO %s\n\n", path);
  cl_free(path);
  fprintf(registry_fd, "# corpus properties provide additional information about the corpus:\n");
  /* lines marked with ##:: are NOT commented out, this is part of the normal registry format! */
  fprintf(registry_fd, "##:: charset  = \"%s\" # character encoding of corpus data\n", corpus_character_set);
  fprintf(registry_fd, "##:: language = \"??\"     # insert ISO code for language (de, en, fr, ...)\n");
  fprintf(registry_fd, "\n\n");

  /* insert p-attributes into registry file */
  fprintf(registry_fd, "##\n## p-attributes (token annotations)\n##\n\n");
  for (i = 0; i < wattr_ptr; i++) {
    fprintf(registry_fd, "ATTRIBUTE %s\n", wattrs[i].name);
  }
  fprintf(registry_fd, "\n\n");

  /* insert s-attributes into registry file */
  fprintf(registry_fd, "##\n## s-attributes (structural markup)\n##\n\n");
  for (i = 0; i < range_ptr; i++) {
    range_print_registry_line(&ranges[i], registry_fd, 1);
  }
  fprintf(registry_fd, "\n");
  fprintf(registry_fd, "# Yours sincerely, the Encode tool.\n");

  fclose(registry_fd);

  cl_free(corpus_name);
  cl_free(info_file);
}



/* *************** *\
 *      MAIN()     *
\* *************** */

/**
 * Main function for cwb-encode.
 *
 * As well as the entry point to the program, this contains
 * the main loop for each line of the corpus to be encoded.
 *
 * The string of each line is sent to one of a number of
 * different functions, depending on what is found in that string!
 *
 * @param argc   Number of command-line arguments.
 * @param argv   Command-line arguments.
 */
int 
main(int argc, char **argv)
{
  int i, j, k, rng, handled;

  char linebuf[MAX_INPUT_LINE_LENGTH];
  char *buf;                    /* 'virtual' buffer; may be advanced to skip leading blanks */
  char separator;
  
  int input_length;             /* length of input line */

  /* initialise global variables */
  progname = "cwb-encode";
  input_files = cl_new_string_list();

  /* parse command-line options */
  encode_parse_options(argc, argv);
  nr_input_files = cl_string_list_size(input_files);

  /* initialisation debug messages */
  if (debug) {
    cl_set_debug_level(1);
    if (nr_input_files > 0) {
      fprintf(stderr, "List of input files:\n");
      for (i = 0; i < nr_input_files; i++)
        fprintf(stderr, " - %s\n", cl_string_list_get(input_files, i));
    }
    else {
      fprintf(stderr, "Reading from standard input.\n");
    }
    encode_print_time(stderr, "Start");
  }

  /* initialise loop variables ... */
  encoding_charset = cl_charset_from_name(corpus_character_set);
  line = 0;
  input_line = 0;

  /* lookup hash for (undeclared) structural attributes (inserted as tokens into corpus) */
  undeclared_sattrs = cl_new_lexhash(REP_CHECK_LEXHASH_SIZE);

  /* MAIN LOOP: read one line of input and process it */
  while ( encode_get_input_line(linebuf, MAX_INPUT_LINE_LENGTH) ) {
    if (verbose && (line % 15000 == 0)) {
      printf("%" COMMA_SEP_THOUSANDS_CONVSPEC "9dk tokens processed\r", line >> 10);
      fflush(stdout);
    }

    input_line++;
    input_length = strlen(linebuf);
    if (input_length >= (MAX_INPUT_LINE_LENGTH - 1)) { /* buffer filled -> line may have been longer */
      encode_error("Input line too long (max: %d characters/bytes).", MAX_INPUT_LINE_LENGTH - 2);
    }

    /* remove trailing line break (LF or CR-LF) */
    cl_string_chomp(linebuf);

    buf = linebuf;
    if (strip_blanks) {
      while (*buf == ' ')
        /* strip leading blanks (trailing blanks will be erased during further processing) */
        buf++;
    }

    /* This bit runs UNLESS either (a) skip_empty_lines (-s) is active and this an empty line;
     * or (b) xml_aware (-x) is active and this line is an XML comment or declaration, i.e. <? or <!
     * To put it another way "if (this is a line that should be encoded)" ...  */
    if ( (! (skip_empty_lines && (buf[0] == '\0')) ) &&                            /* skip empty lines with -s  */
         (! (xml_aware && (buf[0] == '<') && ((buf[1] == '?') || (buf[1] == '!'))) ) /* skip XML declarations/comments with -x  */
      ) {
      /* skip empty lines with -s option (for an empty line, first character will usually be newline) */
      handled = 0;

      if (buf[0] == '<') {
        /* XML tag (may be declared or undeclared s-attribute, start or end tag) */
        k = (buf[1] == '/' ? 2 : 1);
        
        /* identify XML element name (according to slightly relaxed attribute naming conventions!) */
        i = k;                  
        while (cl_xml_is_name_char(buf[i]))
          i++;
        /* first non-valid XML element name character must be whitespace or '>' or '/' (for empty XML element) */
        if (! ((buf[i] == ' ') || (buf[i] == '\t') || (buf[i] == '>') || (buf[i] == '/')) ) 
          i = k;                /* no valid element name found */

        if (i > k) {
          /* looks like a valid XML tag */
          separator = buf[i];   /* terminate string containing element name, but remember original char */
          buf[i] = '\0';        /* so that we can reconstruct the line if we have to insert it literally after all */
          
          if ((rng = range_find(&buf[k])) >= 0) {
            /* good, it's a declared s-attribute and can be handled */
            handled = 1;

            if (ranges[rng].automatic) {
              if (!cl_lexhash_freq(undeclared_sattrs, &buf[k])) {
                fprintf(stderr, "explicit XML tag <%s%s> for implicit s-attribute ignored (", 
                        (k == 1) ? "" : "/", &buf[k]);
                encode_print_input_lineno();
                fprintf(stderr, ", warning issued only once).\n");
                cl_lexhash_add(undeclared_sattrs, &buf[k]); /* can reuse lexhash for undeclared attributes here */
              }
            }
            else {
              if (k == 1) {     /* XML start tag or empty tag */
                i++;            /* identify annotation string, i.e. tag attributes (if there are any) */
                while ((buf[i] == ' ') || (buf[i] == '\t')) /* skip whitespace between element name and first attribute */
                  i++;
                if (separator == '>') {
                  /* tag without annotations: check that there is no extraneous material on the line */
                  if (buf[i] != '\0') {
                    fprintf(stderr, "Warning: extra material after XML tag ignored (");
                    encode_print_input_lineno();
                    fprintf(stderr, ").\n");
                    buf[i] = '\0';
                  }
                }
                else {
                  j = i + strlen(buf+i); /* find '>' character marking end of tag (must be last character on line) */
                  while ((j > i) && (buf[j] == ' ' || buf[j] == '\t' || buf[j] == '\0'))
                    j--; /* set j to last non-blank character on line, which should be '>' */
                  if (buf[j] != '>') {
                    fprintf(stderr, "Malformed XML tag: missing > terminator at end of line (");
                    encode_print_input_lineno();
                    fprintf(stderr, ", annotations will be ignored).\n");
                    buf[i] = '\0'; /* so the annotation string passed to range_open() below is empty */
                  }
                  else {
                    if (buf[j-1] == '/') {
                      j--; /* empty tag: remove "/" from annotation string and handle as an open tag */
                      /* Note that this implicitly closes the previous instance of the empty tag:
                       *  - this means that we can work with empty elements by looking just at the "open-point" of each range;
                       *  - it also means that empty tags with metadata at the start of each text will automatically extend over the full text.
                       * However, the approach sketched here only works with "flat" s-attributes declared without recursion (even without :0). */
                    }
                    buf[j] = '\0';
                  }
                }
                /* start tag: open range */
                range_open(&ranges[rng], line, buf+i);
              }
              else {            /* XML end tag */
                if (separator != '>') {
                  fprintf(stderr, "Warning: no annotations allowed on XML close tag </%s ...> (", &buf[k]);
                  encode_print_input_lineno();
                  fprintf(stderr, ", ignored).\n");
                }
                range_close(&ranges[rng], line - 1); /* end tag belongs to previous line! */
              }
            }
          }
          else {
            /* no appropriate s-attribute declared -> insert tag literally */
            if (!silent) {
              if (!cl_lexhash_freq(undeclared_sattrs, &buf[k])) {
                fprintf(stderr,
                        "s-attribute <%s> not declared, inserted literally (", &buf[k]);
                encode_print_input_lineno();
                fprintf(stderr, ", warning issued only once).\n");
                cl_lexhash_add(undeclared_sattrs, &buf[k]);
              }
            }
            buf[i] = separator; /* restore original line, which will be interpreted as token line */
          }
        }
        /* malformed XML tag (no element name found) */
        else if (!silent) {
          fprintf(stderr, "Malformed tag %s, inserted literally (", buf);
          encode_print_input_lineno();
          fprintf(stderr, ").\n");
        }
      } /* endif line begins with < */
      
      /* if we haven't handled the line so far, it must be data for the positional attributes */
      if (!handled) {
        encode_add_wattr_line(buf);
        line++;                 /* line is now the corpus position of the next token that will be encoded */
        if (line >= CL_MAX_CORPUS_SIZE) {
          /* largest admissible corpus size should be 2^31 - 1 tokens, with maximal cpos = 2^31 - 2 */
          fprintf(stderr, "WARNING: Maximal corpus size has been exceeded.\n");
          fprintf(stderr, "         Input truncated to the first %d tokens (", CL_MAX_CORPUS_SIZE);
          encode_print_input_lineno();
          fprintf(stderr, ").\n");
          break;
        }
      }
    } /* endif (this is a line that should be encoded) */
  } /* endwhile (main loop for each line) */

  if (verbose) {
    printf("%50s\r", "");       /* clear progress line */
    printf("Total size: %" COMMA_SEP_THOUSANDS_CONVSPEC "d tokens (%.1fM)\n", line, ((float) line) / 1048576);
  }

  /* close open regions at end of input; then close file handles for s-attributes */
  for (i = 0; i < range_ptr; i++) {
    Range *rng = &ranges[i];

    if (! rng->null_attribute) { /* don't attempt to close NULL attribute */
    
      /* This is fairly tricky: When multiple end tags are missing for an attribute declared with recursion (even ":0"),
         we have to call range_close() repeatedly to ensure that the open region at the top level is really closed
         (which happens when rng->recursion_level == 1). At the same time, range_close() will also close the corresponding
         ranges of any implicitly defined attributes (used to resolve recursive embedding and element attributes).
         Therefore, the following code calls range_close() repeatedly until the current range is actually closed.
         It also relies on the ordering of the ranges[] array, where top level attributes always precede their children,
         so they should be closed automatically before cleanup reaches them. If the ordering were different, children might
         be closed directly at first, and the following attempt to close them automatically from within the range_close()
         function would produce highly confusing error messages. To be on the safe side (for some definition of safe :-),
         we _never_ close ranges for implicit attributes, and issue a warning if they're still open when cleanup reaches them. 
      */
      if (rng->automatic) {     /* implicitly generated s-attributes should have been closed automatically */
        if (!silent && rng->is_open) {
          fprintf(stderr, "Warning: implicit s-attribute <%s> open at end of input (should not have happened).\n",
                  rng->name);
        }
      }
      else {
        if (rng->is_open) {
          if (rng->recursion_level > 1) 
            fprintf(stderr, "Warning: %d missing </%s> tags inserted at end of input.\n", 
                    rng->recursion_level, rng->name);
          else
            fprintf(stderr, "Warning: missing </%s> tag inserted at end of input.\n", 
                    rng->name);
          
          /* close open region; this will automatically close children from recursion and element attributes;
             if multiple end tags are missing, we have to call range_close() repeatedly until we reach the top level */
          while (rng->is_open) { /* should _not_ create an infinite loop, I hope */
            range_close(rng, line - 1);
          }
        }

        if (!silent && (rng->max_recursion >= 0) && (rng->element_drop_count > 0)) {
          fprintf(stderr, "%7d <%s> regions dropped because of deep nesting.\n",
                  rng->element_drop_count, rng->name);
        }
      }

      /* close file handles for s-attributes */
      if (EOF == fclose(rng->fd)) {
        perror("fclose() failed");
        encode_error("Error writing .rng file for s-attribute <%s>", rng->name);
      }
      if (rng->store_values) {
        if (EOF == fclose(rng->avs)) {
          perror("fclose() failed");
          encode_error("Error writing .avs file for s-attribute <%s>", rng->name);
        }
        if (EOF == fclose(rng->avx)) {
          perror("fclose() failed");
          encode_error("Error writing .avx file for s-attribute <%s>", rng->name);
        }
      }

    }
  } /* endfor: closing each open region and s-attribute filehandle for each Range */

  /* close file handles for positional attributes */
  wattr_close_all();

  /* if registry_file != NULL, write appropriate registry entry to file named <registry_file> */
  if (registry_file != NULL) {
    encode_generate_registry_file(registry_file);
  }

  if (debug)
    encode_print_time(stderr, "Done");

  exit(0);
}

