/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     ID = 258,
     QID = 259,
     NQRID = 260,
     LABEL = 261,
     STRING = 262,
     FLAG = 263,
     TAGSTART = 264,
     TAGEND = 265,
     VARIABLE = 266,
     IPAddress = 267,
     IPSubnet = 268,
     INTEGER = 269,
     DOUBLEFLOAT = 270,
     FIELD = 271,
     FIELDLABEL = 272,
     ANCHORTAG = 273,
     ANCHORENDTAG = 274,
     SEARCH_STRATEGY = 275,
     TAB_SYM = 276,
     CAT_SYM = 277,
     DEFINE_SYM = 278,
     DIFF_SYM = 279,
     DISCARD_SYM = 280,
     EXPAND_SYM = 281,
     EXIT_SYM = 282,
     FLAT_SYM = 283,
     INTER_SYM = 284,
     JOIN_SYM = 285,
     SUBSET_SYM = 286,
     LEFT_SYM = 287,
     RIGHT_SYM = 288,
     SAVE_SYM = 289,
     SCATTER_SYM = 290,
     SHOW_SYM = 291,
     CD_SYM = 292,
     TO_SYM = 293,
     WITHIN_SYM = 294,
     SET_SYM = 295,
     EXEC_SYM = 296,
     CUT_SYM = 297,
     OCCURS_SYM = 298,
     INFO_SYM = 299,
     GROUP_SYM = 300,
     WHERE_SYM = 301,
     ESCAPE_SYM = 302,
     MEET_SYM = 303,
     UNION_SYM = 304,
     MU_SYM = 305,
     SORT_SYM = 306,
     COUNT_SYM = 307,
     ASC_SYM = 308,
     DESC_SYM = 309,
     REVERSE_SYM = 310,
     BY_SYM = 311,
     FOREACH_SYM = 312,
     ON_SYM = 313,
     YES_SYM = 314,
     OFF_SYM = 315,
     NO_SYM = 316,
     SLEEP_SYM = 317,
     REDUCE_SYM = 318,
     MAXIMAL_SYM = 319,
     WITH_SYM = 320,
     WITHOUT_SYM = 321,
     DELETE_SYM = 322,
     SIZE_SYM = 323,
     DUMP_SYM = 324,
     UNDUMP_SYM = 325,
     TABULATE_SYM = 326,
     NOT_SYM = 327,
     CONTAINS_SYM = 328,
     MATCHES_SYM = 329,
     GCDEL = 330,
     APPEND = 331,
     LET = 332,
     GET = 333,
     NEQ = 334,
     IMPLIES = 335,
     RE_PAREN = 336,
     EOL_SYM = 337,
     ELLIPSIS = 338,
     MATCHALL = 339,
     LCSTART = 340,
     LCEND = 341,
     LCMATCHALL = 342,
     EXTENSION = 343,
     PLUSEQ = 344,
     MINUSEQ = 345,
     UNLOCK_SYM = 346,
     USER_SYM = 347,
     HOST_SYM = 348,
     UNDEFINED_MACRO = 349,
     MACRO_SYM = 350,
     RANDOMIZE_SYM = 351,
     FROM_SYM = 352,
     INCLUSIVE_SYM = 353,
     EXCLUSIVE_SYM = 354,
     NULL_SYM = 355
   };
#endif
/* Tokens.  */
#define ID 258
#define QID 259
#define NQRID 260
#define LABEL 261
#define STRING 262
#define FLAG 263
#define TAGSTART 264
#define TAGEND 265
#define VARIABLE 266
#define IPAddress 267
#define IPSubnet 268
#define INTEGER 269
#define DOUBLEFLOAT 270
#define FIELD 271
#define FIELDLABEL 272
#define ANCHORTAG 273
#define ANCHORENDTAG 274
#define SEARCH_STRATEGY 275
#define TAB_SYM 276
#define CAT_SYM 277
#define DEFINE_SYM 278
#define DIFF_SYM 279
#define DISCARD_SYM 280
#define EXPAND_SYM 281
#define EXIT_SYM 282
#define FLAT_SYM 283
#define INTER_SYM 284
#define JOIN_SYM 285
#define SUBSET_SYM 286
#define LEFT_SYM 287
#define RIGHT_SYM 288
#define SAVE_SYM 289
#define SCATTER_SYM 290
#define SHOW_SYM 291
#define CD_SYM 292
#define TO_SYM 293
#define WITHIN_SYM 294
#define SET_SYM 295
#define EXEC_SYM 296
#define CUT_SYM 297
#define OCCURS_SYM 298
#define INFO_SYM 299
#define GROUP_SYM 300
#define WHERE_SYM 301
#define ESCAPE_SYM 302
#define MEET_SYM 303
#define UNION_SYM 304
#define MU_SYM 305
#define SORT_SYM 306
#define COUNT_SYM 307
#define ASC_SYM 308
#define DESC_SYM 309
#define REVERSE_SYM 310
#define BY_SYM 311
#define FOREACH_SYM 312
#define ON_SYM 313
#define YES_SYM 314
#define OFF_SYM 315
#define NO_SYM 316
#define SLEEP_SYM 317
#define REDUCE_SYM 318
#define MAXIMAL_SYM 319
#define WITH_SYM 320
#define WITHOUT_SYM 321
#define DELETE_SYM 322
#define SIZE_SYM 323
#define DUMP_SYM 324
#define UNDUMP_SYM 325
#define TABULATE_SYM 326
#define NOT_SYM 327
#define CONTAINS_SYM 328
#define MATCHES_SYM 329
#define GCDEL 330
#define APPEND 331
#define LET 332
#define GET 333
#define NEQ 334
#define IMPLIES 335
#define RE_PAREN 336
#define EOL_SYM 337
#define ELLIPSIS 338
#define MATCHALL 339
#define LCSTART 340
#define LCEND 341
#define LCMATCHALL 342
#define EXTENSION 343
#define PLUSEQ 344
#define MINUSEQ 345
#define UNLOCK_SYM 346
#define USER_SYM 347
#define HOST_SYM 348
#define UNDEFINED_MACRO 349
#define MACRO_SYM 350
#define RANDOMIZE_SYM 351
#define FROM_SYM 352
#define INCLUSIVE_SYM 353
#define EXCLUSIVE_SYM 354
#define NULL_SYM 355




/* Copy the first part of user declarations.  */
#line 1 "parser.y"

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


#include <sys/time.h>
#ifndef __MINGW__
#include <sys/resource.h>
#endif
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>

#include "../cl/globals.h"
#include "../cl/special-chars.h"
#include "../cl/attributes.h"
#include "../cl/macros.h"

#include "cqp.h"
#include "options.h"
#include "ranges.h"
#include "symtab.h"
#include "treemacros.h"
#include "tree.h"
#include "eval.h"
#include "corpmanag.h"
#include "regex2dfa.h"
#include "builtins.h"
#include "groups.h"
#include "targets.h"
#include "attlist.h"
#include "concordance.h"
#include "output.h"
#include "print-modes.h"
#include "variables.h"

#include "parse_actions.h"

/* CQPserver user authentication */
#include "../CQi/auth.h"

/* macro expansion */
#include "macro.h"
 
/* ============================================================ YACC IF */

extern int yychar;

extern int yylex(void);

void yyerror (char *s)
{
  /* Commented out to avoid a crash */
  /* cqpmessage(Error, "CQP Syntax Error: %s\n\t%s <--", s, QueryBuffer); */
  cqpmessage(Error, "CQP Syntax Error: %s", s); /* slight modification */
  generate_code = 0;
}

void warn_query_lock_violation(void) {
  if (which_app != cqpserver)
    Rprintf("WARNING: query lock violation attempted\n");
  query_lock_violation++;       /* this is for the CQPserver */
}

/* ============================================================ */

/* note: SYCHRONIZE is a windows API identifier, and it doesn't seem at all
necessary here - it is just defined, then tested.
So: commented out. AH 2/4/2010
#define SYNCHRONIZE
*/

void
synchronize(void)
{
/*#if defined(SYNCHRONIZE)*/
  int macro_status;

  /* delete macro buffers & disable macro expansion while sync'ing */
  delete_macro_buffers(1); /* print stack trace on STDERR */
  macro_status = enable_macros;
  enable_macros = 0;

  if (cqp_input_string != NULL) {
    Rprintf("Synchronizing to end of line ... \n");
    while (!(yychar <= 0))
      yychar = yylex();
  }
  else {
    Rprintf("Synchronizing until next ';'...\n");
    while (!(yychar <= 0 || yychar == ';'))
      yychar = yylex();
  }

  enable_macros = macro_status; /* reset enable_macros to previous value */
/*#endif*/
}

#define YYERROR_VERBOSE



/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 118 "parser.y"
{
  Evaltree           evalt;
  Constrainttree     boolt;
  enum b_ops         boolo;
  int                ival;
  double             fval;
  int                index;
  char              *strval;
  CorpusList        *cl;

  struct {
    int a, b;
  } intpair;

  Context            context;
  ActualParamList   *apl;

  enum ctxtdir       direction;

  struct Redir       redir;

  struct InputRedir  in_redir;

  struct {
    int ok;
    int ival;
    char *cval;
  }                  varval;

  struct {
    FieldType field;
    int inclusive;
  }                  base;

  struct {
    char *variableValue;
    char operator;
  }                  varsetting;

  struct {
    int mindist;
    int maxdist;
  }                  Distance;

  struct {
    FieldType anchor;
    int offset;
  }                  Anchor;

  struct {
    FieldType anchor1;
    int offset1;
    FieldType anchor2;
    int offset2;
  }                  AnchorPair;

  struct {
    char *name;
    int flags;
  }                  AttributeSpecification;

  RangeSetOp         rngsetop;

  SortClause         sortclause;

  FieldType          field;

  SearchStrategy     search_strategy;

  TabulationItem     tabulation_item;
}
/* Line 193 of yacc.c.  */
#line 485 "parser.tab.c"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 498 "parser.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  4
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   508

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  123
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  142
/* YYNRULES -- Number of rules.  */
#define YYNRULES  340
/* YYNRULES -- Number of states.  */
#define YYNSTATES  511

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   355

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   105,     2,     2,     2,   112,   104,     2,
     115,   113,   102,   103,   111,   110,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   121,   106,
     109,   107,   108,   114,   118,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   119,     2,   120,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   116,   101,   117,   122,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     4,     7,     9,    10,    11,    15,    16,
      21,    24,    25,    29,    30,    33,    35,    37,    41,    45,
      48,    50,    52,    53,    56,    58,    60,    62,    64,    66,
      68,    70,    72,    74,    76,    78,    80,    82,    84,    86,
      88,    90,    92,    94,    96,    98,   100,   103,   105,   109,
     117,   121,   125,   129,   131,   132,   135,   138,   140,   141,
     144,   146,   149,   152,   155,   158,   160,   163,   166,   171,
     175,   177,   179,   181,   183,   185,   186,   193,   196,   199,
     201,   204,   208,   213,   218,   223,   226,   229,   232,   235,
     239,   242,   244,   249,   254,   255,   267,   271,   272,   274,
     276,   277,   279,   281,   283,   285,   287,   289,   291,   294,
     297,   300,   302,   313,   321,   323,   325,   328,   331,   333,
     334,   335,   341,   342,   352,   354,   358,   361,   363,   367,
     370,   371,   375,   380,   386,   388,   389,   396,   399,   404,
     405,   407,   409,   410,   412,   413,   419,   425,   429,   434,
     436,   437,   442,   447,   451,   453,   457,   460,   464,   467,
     471,   479,   485,   487,   488,   489,   492,   496,   498,   502,
     504,   506,   508,   514,   518,   519,   524,   526,   527,   528,
     529,   530,   537,   541,   543,   546,   548,   551,   554,   557,
     560,   562,   566,   568,   570,   572,   576,   581,   587,   589,
     591,   594,   600,   602,   604,   606,   608,   609,   613,   615,
     616,   618,   619,   621,   623,   626,   628,   632,   634,   638,
     640,   642,   643,   646,   647,   648,   655,   656,   658,   659,
     663,   664,   667,   668,   670,   671,   673,   674,   676,   678,
     679,   684,   685,   687,   689,   690,   693,   695,   697,   698,
     700,   702,   704,   708,   712,   716,   720,   723,   725,   729,
     734,   736,   739,   742,   744,   745,   747,   750,   752,   755,
     757,   759,   761,   764,   769,   771,   773,   775,   777,   779,
     781,   783,   785,   787,   792,   794,   798,   800,   802,   804,
     806,   808,   815,   818,   820,   821,   827,   831,   834,   838,
     839,   843,   849,   854,   856,   858,   860,   861,   862,   868,
     871,   874,   877,   881,   882,   885,   886,   894,   902,   907,
     909,   910,   913,   917,   924,   927,   929,   931,   934,   935,
     937,   938,   940,   941,   943,   944,   946,   947,   949,   950,
     952
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int16 yyrhs[] =
{
     124,     0,    -1,    -1,   125,   126,    -1,   106,    -1,    -1,
      -1,   127,   131,   106,    -1,    -1,    91,    14,   128,   106,
      -1,   136,   106,    -1,    -1,   129,   135,   106,    -1,    -1,
     130,    27,    -1,     1,    -1,   132,    -1,     3,   107,   132,
      -1,     3,   107,   146,    -1,   133,   224,    -1,   228,    -1,
     147,    -1,    -1,   134,   189,    -1,   143,    -1,   137,    -1,
     138,    -1,   151,    -1,   181,    -1,   179,    -1,   156,    -1,
     154,    -1,   153,    -1,   161,    -1,   162,    -1,   163,    -1,
     173,    -1,   183,    -1,   184,    -1,   185,    -1,   186,    -1,
     166,    -1,   249,    -1,   253,    -1,   255,    -1,   257,    -1,
      47,   258,    -1,    82,    -1,    22,   227,   139,    -1,    22,
     227,   261,    14,   262,    14,   139,    -1,    22,   147,   139,
      -1,    22,     7,   139,    -1,    34,   227,   139,    -1,   140,
      -1,    -1,   108,     7,    -1,    76,     7,    -1,   142,    -1,
      -1,   109,     7,    -1,    36,    -1,    36,     3,    -1,    36,
     144,    -1,    36,    37,    -1,   144,   145,    -1,   145,    -1,
     103,     3,    -1,   110,     3,    -1,    97,   228,    38,     3,
      -1,   148,   228,   228,    -1,   149,    -1,    30,    -1,    49,
      -1,    29,    -1,    24,    -1,    -1,    31,   227,    46,    17,
     150,   211,    -1,    25,   152,    -1,   152,   228,    -1,   228,
      -1,    36,    11,    -1,    23,    11,   155,    -1,    23,    11,
      89,    11,    -1,    23,    11,    90,    11,    -1,    23,    11,
     107,    11,    -1,   109,     7,    -1,   107,     7,    -1,    89,
       7,    -1,    90,     7,    -1,    40,     3,   160,    -1,    40,
       3,    -1,    40,    -1,    40,   228,    16,    16,    -1,    40,
     228,    16,   100,    -1,    -1,    40,   228,    16,    20,   157,
     211,    39,   225,   220,     3,   158,    -1,    97,    16,   159,
      -1,    -1,    98,    -1,    99,    -1,    -1,    14,    -1,     7,
      -1,     3,    -1,    58,    -1,    59,    -1,    60,    -1,    61,
      -1,    14,     3,    -1,    41,     7,    -1,    44,   228,    -1,
      44,    -1,    45,   228,   264,     3,   164,   264,     3,   219,
     165,   139,    -1,    45,   228,   264,     3,   219,   165,   139,
      -1,    56,    -1,    57,    -1,    45,    56,    -1,    45,    57,
      -1,    26,    -1,    -1,    -1,    71,   228,   167,   169,   139,
      -1,    -1,    71,   228,   168,   261,    14,   262,    14,   169,
     139,    -1,   170,    -1,   169,   111,   170,    -1,   171,   172,
      -1,   264,    -1,   264,   263,   264,    -1,     3,   213,    -1,
      -1,    51,   227,   174,    -1,    51,   227,    96,   221,    -1,
      52,   227,   175,   219,   139,    -1,   175,    -1,    -1,    56,
       3,   213,   176,   177,   178,    -1,   260,   264,    -1,   260,
     264,   263,   264,    -1,    -1,    53,    -1,    54,    -1,    -1,
      55,    -1,    -1,    63,   227,    38,   222,   180,    -1,    63,
     227,    38,    64,    74,    -1,    42,   227,   222,    -1,    42,
     227,   222,   222,    -1,   112,    -1,    -1,    67,   227,    65,
      16,    -1,    67,   227,    66,    16,    -1,    67,   227,   182,
      -1,    14,    -1,    14,    83,    14,    -1,    62,    14,    -1,
      68,   228,   259,    -1,    68,    11,    -1,    69,   228,   139,
      -1,    69,   228,   261,    14,   262,    14,   139,    -1,    70,
       3,   188,   187,   141,    -1,    53,    -1,    -1,    -1,    65,
      16,    -1,    65,    16,    16,    -1,   190,    -1,   190,    51,
     174,    -1,   191,    -1,   193,    -1,   245,    -1,   192,   195,
     215,   219,   194,    -1,    88,     3,   113,    -1,    -1,    50,
     241,   194,   219,    -1,   105,    -1,    -1,    -1,    -1,    -1,
     196,   199,   197,   214,   198,   218,    -1,   199,   101,   200,
      -1,   200,    -1,   200,   201,    -1,   201,    -1,   202,   203,
      -1,   202,   102,    -1,   202,   103,    -1,   202,   114,    -1,
     202,    -1,   115,   199,   113,    -1,   207,    -1,   205,    -1,
     204,    -1,   116,   222,   117,    -1,   116,   111,   222,   117,
      -1,   116,   222,   111,   223,   117,    -1,    18,    -1,    19,
      -1,     9,   108,    -1,     9,   206,     7,   213,   108,    -1,
      10,    -1,   232,    -1,    79,    -1,   107,    -1,    -1,   208,
     209,   210,    -1,   118,    -1,    -1,     6,    -1,    -1,   211,
      -1,   212,    -1,     7,   213,    -1,    11,    -1,   119,   230,
     120,    -1,    84,    -1,    85,   230,    86,    -1,    87,    -1,
       8,    -1,    -1,    75,   230,    -1,    -1,    -1,   215,   121,
       3,   216,   217,   195,    -1,    -1,   105,    -1,    -1,    39,
     225,   226,    -1,    -1,    42,   222,    -1,    -1,    14,    -1,
      -1,    14,    -1,    -1,    14,    -1,   222,    -1,    -1,    26,
     225,    38,   226,    -1,    -1,    32,    -1,    33,    -1,    -1,
     220,     3,    -1,    14,    -1,   228,    -1,    -1,   229,    -1,
       5,    -1,     3,    -1,   230,    80,   230,    -1,   230,   101,
     230,    -1,   230,   104,   230,    -1,   115,   230,   113,    -1,
     105,   230,    -1,   231,    -1,   234,   236,   235,    -1,   234,
     232,     7,   213,    -1,   234,    -1,   233,    73,    -1,   233,
      74,    -1,    72,    -1,    -1,   240,    -1,   122,   240,    -1,
       3,    -1,   122,     3,    -1,    16,    -1,   237,    -1,   234,
      -1,     7,   213,    -1,    81,    11,   113,   213,    -1,    11,
      -1,    14,    -1,    15,    -1,   109,    -1,   108,    -1,   107,
      -1,    79,    -1,    77,    -1,    78,    -1,     3,   115,   238,
     113,    -1,   239,    -1,   238,   111,   239,    -1,   235,    -1,
       4,    -1,   242,    -1,   244,    -1,   210,    -1,   115,    48,
     241,   241,   243,   113,    -1,    14,    14,    -1,     3,    -1,
      -1,   115,    49,   241,   241,   113,    -1,    21,   246,   218,
      -1,   210,   247,    -1,   247,   248,   210,    -1,    -1,   116,
     222,   117,    -1,   116,   222,   111,   223,   117,    -1,   116,
     111,   222,   117,    -1,   102,    -1,   103,    -1,   114,    -1,
      -1,    -1,    92,     3,     7,   250,   251,    -1,    93,    12,
      -1,    93,    13,    -1,    93,   102,    -1,   115,   252,   113,
      -1,    -1,   252,     3,    -1,    -1,   254,    95,     3,   115,
      14,   113,   256,    -1,   254,    95,     3,   115,     7,   113,
     256,    -1,   254,    95,   109,     7,    -1,    23,    -1,    -1,
      36,    95,    -1,    36,    95,     3,    -1,    36,    95,     3,
     115,    14,   113,    -1,   256,     7,    -1,     7,    -1,    96,
      -1,    96,    14,    -1,    -1,    16,    -1,    -1,    58,    -1,
      -1,    97,    -1,    -1,    38,    -1,    -1,    83,    -1,    -1,
      16,    -1,    16,   119,    14,   120,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   364,   364,   364,   369,   370,   373,   373,   376,   375,
     386,   387,   387,   389,   389,   391,   400,   401,   403,   407,
     410,   411,   412,   412,   416,   417,   418,   419,   420,   421,
     422,   423,   424,   425,   426,   427,   428,   429,   430,   431,
     432,   433,   434,   435,   436,   437,   438,   443,   446,   449,
     452,   458,   462,   466,   467,   473,   477,   484,   485,   490,
     495,   498,   512,   514,   519,   521,   525,   526,   529,   532,
     533,   536,   537,   538,   539,   545,   542,   561,   564,   567,
     572,   577,   586,   592,   598,   609,   610,   611,   612,   615,
     641,   649,   652,   653,   654,   654,   671,   676,   682,   683,
     684,   687,   691,   695,   700,   701,   702,   703,   704,   711,
     714,   715,   718,   724,   732,   733,   734,   735,   740,   741,
     749,   748,   753,   752,   758,   760,   764,   776,   778,   783,
     786,   795,   807,   818,   830,   831,   834,   852,   853,   857,
     861,   862,   863,   866,   867,   872,   876,   880,   884,   890,
     891,   894,   900,   906,   914,   915,   920,   924,   925,   929,
     931,   935,   939,   940,   944,   945,   950,   961,   962,   972,
     973,   974,   977,   983,   984,   987,   991,   992,   995,   996,
     998,   995,  1002,  1004,  1007,  1009,  1012,  1020,  1025,  1030,
    1035,  1038,  1039,  1044,  1049,  1056,  1061,  1066,  1079,  1080,
    1083,  1084,  1086,  1089,  1090,  1091,  1092,  1095,  1101,  1102,
    1105,  1106,  1109,  1110,  1113,  1114,  1126,  1127,  1135,  1136,
    1144,  1175,  1178,  1179,  1184,  1183,  1189,  1192,  1193,  1196,
    1204,  1212,  1213,  1216,  1217,  1220,  1221,  1225,  1233,  1234,
    1237,  1243,  1250,  1251,  1252,  1255,  1256,  1259,  1260,  1264,
    1284,  1285,  1288,  1289,  1290,  1291,  1292,  1293,  1296,  1297,
    1306,  1309,  1310,  1313,  1314,  1317,  1318,  1319,  1320,  1321,
    1322,  1325,  1326,  1327,  1331,  1348,  1356,  1366,  1367,  1368,
    1369,  1370,  1371,  1374,  1377,  1379,  1396,  1407,  1410,  1411,
    1412,  1420,  1427,  1433,  1434,  1442,  1448,  1454,  1458,  1462,
    1465,  1466,  1477,  1478,  1479,  1480,  1481,  1485,  1485,  1487,
    1488,  1489,  1492,  1493,  1497,  1499,  1503,  1512,  1522,  1528,
    1529,  1533,  1536,  1540,  1549,  1559,  1563,  1564,  1567,  1572,
    1573,  1575,  1576,  1578,  1579,  1581,  1582,  1584,  1585,  1589,
    1590
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "ID", "QID", "NQRID", "LABEL", "STRING",
  "FLAG", "TAGSTART", "TAGEND", "VARIABLE", "IPAddress", "IPSubnet",
  "INTEGER", "DOUBLEFLOAT", "FIELD", "FIELDLABEL", "ANCHORTAG",
  "ANCHORENDTAG", "SEARCH_STRATEGY", "TAB_SYM", "CAT_SYM", "DEFINE_SYM",
  "DIFF_SYM", "DISCARD_SYM", "EXPAND_SYM", "EXIT_SYM", "FLAT_SYM",
  "INTER_SYM", "JOIN_SYM", "SUBSET_SYM", "LEFT_SYM", "RIGHT_SYM",
  "SAVE_SYM", "SCATTER_SYM", "SHOW_SYM", "CD_SYM", "TO_SYM", "WITHIN_SYM",
  "SET_SYM", "EXEC_SYM", "CUT_SYM", "OCCURS_SYM", "INFO_SYM", "GROUP_SYM",
  "WHERE_SYM", "ESCAPE_SYM", "MEET_SYM", "UNION_SYM", "MU_SYM", "SORT_SYM",
  "COUNT_SYM", "ASC_SYM", "DESC_SYM", "REVERSE_SYM", "BY_SYM",
  "FOREACH_SYM", "ON_SYM", "YES_SYM", "OFF_SYM", "NO_SYM", "SLEEP_SYM",
  "REDUCE_SYM", "MAXIMAL_SYM", "WITH_SYM", "WITHOUT_SYM", "DELETE_SYM",
  "SIZE_SYM", "DUMP_SYM", "UNDUMP_SYM", "TABULATE_SYM", "NOT_SYM",
  "CONTAINS_SYM", "MATCHES_SYM", "GCDEL", "APPEND", "LET", "GET", "NEQ",
  "IMPLIES", "RE_PAREN", "EOL_SYM", "ELLIPSIS", "MATCHALL", "LCSTART",
  "LCEND", "LCMATCHALL", "EXTENSION", "PLUSEQ", "MINUSEQ", "UNLOCK_SYM",
  "USER_SYM", "HOST_SYM", "UNDEFINED_MACRO", "MACRO_SYM", "RANDOMIZE_SYM",
  "FROM_SYM", "INCLUSIVE_SYM", "EXCLUSIVE_SYM", "NULL_SYM", "'|'", "'*'",
  "'+'", "'&'", "'!'", "';'", "'='", "'>'", "'<'", "'-'", "','", "'%'",
  "')'", "'?'", "'('", "'{'", "'}'", "'@'", "'['", "']'", "':'", "'~'",
  "$accept", "line", "@1", "command", "@2", "@3", "@4", "@5",
  "CorpusCommand", "UnnamedCorpusCommand", "CYCommand", "@6",
  "InteractiveCommand", "EOLCmd", "Cat", "Saving", "OptionalRedir",
  "Redir", "OptionalInputRedir", "InputRedir", "Showing",
  "AttributeSelections", "AttributeSelection", "TranslateExpr",
  "CorpusSetExpr", "SetOp", "SubsetExpr", "@7", "Discard", "DiscArgs",
  "VarPrintCmd", "VarDefCmd", "VariableValueSpec", "OptionSetCmd", "@8",
  "OptBase", "InclusiveExclusive", "VarValue", "ExecCmd", "InfoCmd",
  "GroupCmd", "GroupBy", "OptExpansion", "TabulateCmd", "@9", "@10",
  "TabulationItems", "TabulationItem", "TabulationRange",
  "OptAttributeSpec", "SortCmd", "OptionalSortClause", "SortClause",
  "SortBoundaries", "SortDirection", "OptReverse", "Reduction",
  "OptPercent", "Delete", "LineRange", "SleepCmd", "SizeCmd", "DumpCmd",
  "UndumpCmd", "OptAscending", "OptWithTargetKeyword", "Query", "AQuery",
  "StandardQuery", "EmbeddedModifier", "MUQuery", "OptKeep",
  "SearchPattern", "@11", "@12", "@13", "RegWordfExpr", "RegWordfTerm",
  "RegWordfFactor", "RegWordfPower", "Repeat", "AnchorPoint", "XMLTag",
  "RegexpOp", "NamedWfPattern", "OptTargetSign", "OptRefId",
  "WordformPattern", "ExtConstraint", "LookaheadConstraint",
  "OptionalFlag", "GlobalConstraint", "AlignmentConstraints", "@14",
  "OptNot", "SearchSpace", "CutStatement", "OptNumber", "OptInteger",
  "PosInt", "OptMaxNumber", "ReStructure", "OptDirection", "Description",
  "OptionalCID", "CID", "ID_OR_NQRID", "BoolExpr", "RelExpr", "MvalOp",
  "OptionalNot", "RelLHS", "RelRHS", "RelOp", "FunctionCall",
  "FunctionArgList", "SingleArg", "LabelReference", "MUStatement",
  "MeetStatement", "MeetContext", "UnionStatement", "TABQuery",
  "TabPatterns", "TabOtherPatterns", "OptDistance", "AuthorizeCmd", "@15",
  "OptionalGrants", "Grants", "Macro", "OptDEFINE_SYM", "ShowMacro",
  "MultiString", "RandomizeCmd", "OtherCommand", "OptionalFIELD", "OptON",
  "OptFROM", "OptTO", "OptELLIPSIS", "Anchor", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   124,    42,    43,    38,    33,    59,    61,    62,    60,
      45,    44,    37,    41,    63,    40,   123,   125,    64,    91,
      93,    58,   126
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint16 yyr1[] =
{
       0,   123,   125,   124,   124,   124,   127,   126,   128,   126,
     126,   129,   126,   130,   126,   126,   131,   131,   131,   132,
     133,   133,   134,   133,   135,   135,   135,   135,   135,   135,
     135,   135,   135,   135,   135,   135,   135,   135,   135,   135,
     135,   135,   135,   135,   135,   135,   135,   136,   137,   137,
     137,   137,   138,   139,   139,   140,   140,   141,   141,   142,
     143,   143,   143,   143,   144,   144,   145,   145,   146,   147,
     147,   148,   148,   148,   148,   150,   149,   151,   152,   152,
     153,   154,   154,   154,   154,   155,   155,   155,   155,   156,
     156,   156,   156,   156,   157,   156,   158,   158,   159,   159,
     159,   160,   160,   160,   160,   160,   160,   160,   160,   161,
     162,   162,   163,   163,   164,   164,   164,   164,   165,   165,
     167,   166,   168,   166,   169,   169,   170,   171,   171,   172,
     172,   173,   173,   173,   174,   174,   175,   176,   176,   176,
     177,   177,   177,   178,   178,   179,   179,   179,   179,   180,
     180,   181,   181,   181,   182,   182,   183,   184,   184,   185,
     185,   186,   187,   187,   188,   188,   188,   189,   189,   190,
     190,   190,   191,   192,   192,   193,   194,   194,   196,   197,
     198,   195,   199,   199,   200,   200,   201,   201,   201,   201,
     201,   202,   202,   202,   202,   203,   203,   203,   204,   204,
     205,   205,   205,   206,   206,   206,   206,   207,   208,   208,
     209,   209,   210,   210,   211,   211,   211,   211,   212,   212,
     213,   213,   214,   214,   216,   215,   215,   217,   217,   218,
     218,   219,   219,   220,   220,   221,   221,   222,   223,   223,
     224,   224,   225,   225,   225,   226,   226,   227,   227,   228,
     229,   229,   230,   230,   230,   230,   230,   230,   231,   231,
     231,   232,   232,   233,   233,   234,   234,   234,   234,   234,
     234,   235,   235,   235,   235,   235,   235,   236,   236,   236,
     236,   236,   236,   237,   238,   238,   239,   240,   241,   241,
     241,   242,   243,   243,   243,   244,   245,   246,   247,   247,
     248,   248,   248,   248,   248,   248,   248,   250,   249,   249,
     249,   249,   251,   251,   252,   252,   253,   253,   253,   254,
     254,   255,   255,   255,   256,   256,   257,   257,   258,   259,
     259,   260,   260,   261,   261,   262,   262,   263,   263,   264,
     264
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     1,     0,     0,     3,     0,     4,
       2,     0,     3,     0,     2,     1,     1,     3,     3,     2,
       1,     1,     0,     2,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     2,     1,     3,     7,
       3,     3,     3,     1,     0,     2,     2,     1,     0,     2,
       1,     2,     2,     2,     2,     1,     2,     2,     4,     3,
       1,     1,     1,     1,     1,     0,     6,     2,     2,     1,
       2,     3,     4,     4,     4,     2,     2,     2,     2,     3,
       2,     1,     4,     4,     0,    11,     3,     0,     1,     1,
       0,     1,     1,     1,     1,     1,     1,     1,     2,     2,
       2,     1,    10,     7,     1,     1,     2,     2,     1,     0,
       0,     5,     0,     9,     1,     3,     2,     1,     3,     2,
       0,     3,     4,     5,     1,     0,     6,     2,     4,     0,
       1,     1,     0,     1,     0,     5,     5,     3,     4,     1,
       0,     4,     4,     3,     1,     3,     2,     3,     2,     3,
       7,     5,     1,     0,     0,     2,     3,     1,     3,     1,
       1,     1,     5,     3,     0,     4,     1,     0,     0,     0,
       0,     6,     3,     1,     2,     1,     2,     2,     2,     2,
       1,     3,     1,     1,     1,     3,     4,     5,     1,     1,
       2,     5,     1,     1,     1,     1,     0,     3,     1,     0,
       1,     0,     1,     1,     2,     1,     3,     1,     3,     1,
       1,     0,     2,     0,     0,     6,     0,     1,     0,     3,
       0,     2,     0,     1,     0,     1,     0,     1,     1,     0,
       4,     0,     1,     1,     0,     2,     1,     1,     0,     1,
       1,     1,     3,     3,     3,     3,     2,     1,     3,     4,
       1,     2,     2,     1,     0,     1,     2,     1,     2,     1,
       1,     1,     2,     4,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     4,     1,     3,     1,     1,     1,     1,
       1,     6,     2,     1,     0,     5,     3,     2,     3,     0,
       3,     5,     4,     1,     1,     1,     0,     0,     5,     2,
       2,     2,     3,     0,     2,     0,     7,     7,     4,     1,
       0,     2,     3,     6,     2,     1,     1,     2,     0,     1,
       0,     1,     0,     1,     0,     1,     0,     1,     0,     1,
       4
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint16 yydefact[] =
{
       2,     4,     0,     0,     1,    15,    47,     0,     3,    22,
     320,     0,     0,     8,   251,   250,    74,    73,    71,   248,
      72,     0,    16,   241,   174,    21,     0,    70,    20,   249,
     248,   319,     0,   248,    60,    91,     0,   248,   111,     0,
     328,   248,   248,     0,   248,   248,     0,     0,     0,     0,
       0,     0,   326,     0,    25,    26,    24,    27,    32,    31,
      30,    33,    34,    35,    41,    36,    29,    28,    37,    38,
      39,    40,    42,    43,     0,    44,    45,    14,    10,     0,
      22,   251,     0,   247,     7,   244,    19,     0,     0,     0,
      23,   167,   169,   178,   170,   171,     0,    54,    54,    54,
       0,    77,    79,    54,    61,    80,    63,   321,     0,     0,
      62,    65,    90,     0,   109,     0,   110,     0,    46,   135,
       0,   156,     0,     0,   158,   330,    54,   164,   122,     0,
     309,   310,   311,   327,    12,     0,     9,     0,    17,    18,
       0,   242,   243,     0,   221,   215,   217,     0,   219,     0,
     299,   212,   213,   230,     0,   290,   177,   288,   289,     0,
     135,   226,   209,    69,     0,     0,    51,    53,    50,   333,
      48,     0,     0,     0,     0,     0,    81,    78,    52,   322,
      66,    67,    64,   103,   102,   101,   104,   105,   106,   107,
      89,     0,   237,   147,   339,     0,     0,   236,   131,   134,
     232,     0,   154,     0,     0,   153,   329,   157,   159,     0,
       0,   163,     0,   334,   307,     0,     0,     0,    75,   234,
     220,   214,   267,   287,   269,     0,     0,     0,     0,   257,
     260,   270,   265,     0,   306,   244,   296,     0,     0,   176,
     232,   173,   168,   232,   264,   202,   198,   199,   209,   208,
     179,   183,   185,   190,   194,   193,   192,   211,    56,    55,
     336,    87,    82,    88,    83,    86,    84,    85,     0,   108,
      92,    94,    93,   148,     0,   232,   221,   235,   132,     0,
      54,     0,   150,     0,   151,   152,   336,   165,   162,    58,
      54,   124,   130,   127,     0,   313,     0,   318,     0,     0,
     246,     0,   240,     0,   256,     0,   268,   266,     0,   218,
       0,     0,   263,   281,   282,   280,   279,   278,   277,     0,
       0,     0,   216,   303,   304,   305,     0,     0,   234,     0,
       0,   175,     0,   177,   204,   205,   200,     0,   203,     0,
     209,   223,   184,   187,   188,   189,     0,   186,   210,     0,
     335,     0,     0,     0,     0,     0,   114,   115,     0,   119,
     139,   231,   133,   146,   149,   145,   155,     0,   166,     0,
     161,    57,     0,   121,   221,   126,   337,     0,   336,   315,
     308,     0,     0,    68,    76,   245,   221,   274,   275,   276,
       0,   271,   286,     0,   284,   255,   252,   253,   254,   221,
     261,   262,   258,     0,     0,   298,   229,   294,     0,   224,
     172,   221,   191,   182,     0,   180,     0,     0,   207,    54,
     323,     0,   340,   116,   117,     0,   118,    54,   331,   142,
       0,    54,    59,   125,   129,   128,     0,     0,     0,     0,
     272,     0,     0,   283,   259,     0,   239,   300,   293,     0,
       0,   295,   228,     0,   222,   230,     0,   239,   195,    49,
     244,   232,   113,   140,   141,   144,   137,   160,     0,   314,
     312,   325,   317,   316,   221,   285,   302,   238,     0,   292,
     291,   227,   178,   201,   181,   196,     0,   234,   119,   143,
     136,     0,    54,   324,   273,   301,   225,   197,   233,     0,
      54,   138,   123,    97,   112,     0,    95,   100,    98,    99,
      96
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     2,     3,     8,     9,    79,    10,    11,    21,    22,
      23,    24,    53,    12,    54,    55,   166,   167,   370,   371,
      56,   110,   111,   139,    25,    26,    27,   299,    57,   101,
      58,    59,   176,    60,   353,   506,   510,   190,    61,    62,
      63,   358,   427,    64,   212,   213,   290,   291,   292,   375,
      65,   198,   199,   429,   465,   490,    66,   365,    67,   205,
      68,    69,    70,    71,   289,   211,    90,    91,    92,    93,
      94,   240,   161,   162,   341,   455,   250,   251,   252,   253,
     347,   254,   255,   337,   256,   257,   349,   155,   151,   152,
     221,   415,   243,   452,   482,   236,   280,   301,   278,   477,
     478,    86,   143,   302,    82,    83,    29,   228,   229,   319,
     320,   230,   392,   321,   231,   393,   394,   232,   156,   157,
     450,   158,    95,   153,   234,   327,    72,   295,   380,   437,
      73,    74,    75,   472,    76,   118,   207,   430,   171,   351,
     377,   293
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -292
static const yytype_int16 yypact[] =
{
      21,  -292,    44,   321,  -292,  -292,  -292,    61,  -292,   206,
     401,    71,    -5,  -292,    15,  -292,  -292,  -292,  -292,   192,
    -292,    47,  -292,   130,   113,  -292,   192,  -292,  -292,  -292,
     197,   204,   192,   192,    84,   297,   158,   192,   192,   192,
    -292,   192,   192,   218,   192,   192,   188,   192,   217,   192,
     237,    59,   236,   176,  -292,  -292,  -292,  -292,  -292,  -292,
    -292,  -292,  -292,  -292,  -292,  -292,  -292,  -292,  -292,  -292,
    -292,  -292,  -292,  -292,   190,  -292,  -292,  -292,  -292,   183,
     159,  -292,   246,  -292,  -292,   105,  -292,    39,    30,   293,
    -292,   252,  -292,  -292,  -292,  -292,   192,   -22,   -22,    49,
      96,   192,  -292,   -22,  -292,  -292,  -292,   304,   313,   314,
     103,  -292,   209,   302,  -292,   305,  -292,   307,  -292,    -3,
     264,  -292,   287,    38,  -292,   318,    49,   270,   320,   330,
    -292,  -292,  -292,  -292,  -292,    20,  -292,   192,  -292,  -292,
     324,  -292,  -292,   300,   339,  -292,  -292,    13,  -292,    13,
    -292,  -292,  -292,   310,   132,  -292,   248,  -292,  -292,   241,
     264,  -292,    60,  -292,   349,   351,  -292,  -292,  -292,  -292,
    -292,   345,   143,   234,   242,   353,  -292,  -292,  -292,   249,
    -292,  -292,  -292,  -292,  -292,   364,  -292,  -292,  -292,  -292,
    -292,    31,  -292,   305,   250,   371,   372,   363,  -292,  -292,
     336,    41,   296,   365,   366,  -292,  -292,  -292,  -292,   373,
     369,   327,   307,   289,  -292,   278,   388,   358,  -292,   383,
    -292,  -292,   283,  -292,  -292,    13,    13,   260,   138,  -292,
     201,  -292,  -292,   118,   145,   105,  -292,    30,    30,  -292,
     336,  -292,  -292,    -8,    69,  -292,  -292,  -292,    60,  -292,
     298,    55,  -292,   185,  -292,  -292,  -292,   394,  -292,  -292,
     377,  -292,  -292,  -292,  -292,  -292,  -292,  -292,   387,  -292,
    -292,  -292,  -292,  -292,   390,    40,   339,  -292,  -292,   305,
     -22,   333,   299,   396,  -292,  -292,   377,   402,  -292,   311,
       1,  -292,   416,    16,   407,   312,   203,  -292,   419,    24,
     422,   425,  -292,    11,  -292,   180,  -292,  -292,    13,  -292,
      13,    13,  -292,  -292,  -292,  -292,  -292,  -292,  -292,   423,
     232,    11,  -292,  -292,  -292,  -292,    -1,    39,   383,    30,
      30,  -292,   426,   248,  -292,  -292,  -292,   424,  -292,   -13,
      60,   357,  -292,  -292,  -292,  -292,    25,  -292,  -292,    39,
    -292,   420,   325,    24,   329,   220,  -292,  -292,   307,   418,
      22,  -292,  -292,  -292,  -292,  -292,  -292,   433,  -292,   443,
    -292,  -292,   307,  -292,   339,  -292,  -292,   307,   377,  -292,
    -292,   338,   341,  -292,  -292,  -292,   339,  -292,  -292,  -292,
     444,  -292,  -292,   -53,  -292,  -292,    51,   352,  -292,   339,
    -292,  -292,  -292,   305,   141,  -292,  -292,    88,   344,  -292,
    -292,   339,  -292,    55,    13,  -292,   305,   154,  -292,   -22,
    -292,   421,  -292,  -292,  -292,   430,  -292,   -22,  -292,   258,
     307,   -22,  -292,  -292,  -292,  -292,   445,    17,   451,   451,
    -292,   348,    11,  -292,  -292,   350,   305,  -292,  -292,   448,
     360,  -292,   361,   367,   153,   310,   359,   305,  -292,  -292,
     105,   336,  -292,  -292,  -292,   410,    16,  -292,   307,  -292,
    -292,  -292,   467,   467,   339,  -292,  -292,  -292,   362,  -292,
    -292,  -292,  -292,  -292,  -292,  -292,   368,   463,   418,  -292,
    -292,   307,     1,  -292,  -292,  -292,  -292,  -292,  -292,   475,
     -22,  -292,  -292,   384,  -292,   464,  -292,   215,  -292,  -292,
    -292
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -292,  -292,  -292,  -292,  -292,  -292,  -292,  -292,  -292,   403,
    -292,  -292,  -292,  -292,  -292,  -292,   -98,  -292,  -292,  -292,
    -292,  -292,   374,  -292,   452,  -292,  -292,  -292,  -292,  -292,
    -292,  -292,  -292,  -292,  -292,  -292,  -292,  -292,  -292,  -292,
    -292,  -292,    -2,  -292,  -292,  -292,    19,   116,  -292,  -292,
    -292,   331,   370,  -292,  -292,  -292,  -292,  -292,  -292,  -292,
    -292,  -292,  -292,  -292,  -292,  -292,  -292,  -292,  -292,  -292,
    -292,   156,    14,  -292,  -292,  -292,   244,   155,  -241,  -292,
    -292,  -292,  -292,  -292,  -292,  -292,  -292,   -83,  -259,  -292,
    -267,  -292,  -292,  -292,  -292,    43,  -232,    12,  -292,  -112,
      45,  -292,  -229,   172,   253,    10,  -292,  -142,  -292,   257,
    -292,  -291,   182,  -292,  -292,  -292,    62,   279,  -170,  -292,
    -292,  -292,  -292,  -292,  -292,  -292,  -292,  -292,  -292,  -292,
    -292,  -292,  -292,    66,  -292,  -292,  -292,  -292,   -93,  -262,
      42,  -115
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -339
static const yytype_int16 yytable[] =
{
     168,   170,   195,   193,   150,   178,   328,   233,   331,   360,
     342,   333,   391,   192,   222,   223,   222,   223,   386,    28,
     469,    -5,   387,   215,   367,   388,   389,   224,   208,   224,
     391,   144,  -338,   209,   279,   145,    96,   144,  -332,   192,
     384,   145,   102,   359,     4,   113,   144,   270,   116,   117,
     145,   271,   202,   196,   164,   192,   125,   126,   442,   128,
     443,  -209,  -209,  -334,   244,   245,  -209,   329,   330,   244,
     245,   130,   131,   246,   247,    13,  -206,   164,   246,   247,
     428,   273,   279,   304,   305,   355,   165,   104,   340,   282,
      28,   448,   390,   197,   421,   105,   356,   357,    77,   376,
     412,    78,   449,   203,   204,   281,   163,   434,   146,   165,
     403,   177,   372,   332,   146,   147,   436,   148,   225,   440,
     294,   106,    80,   146,   147,   164,   148,     1,   226,   216,
     470,   272,   444,   227,    87,   227,   416,   141,   142,  -209,
    -209,   312,  -209,   149,   453,   154,   169,   217,   334,   149,
     261,   391,   310,    84,   262,   311,    85,   165,   149,   407,
     408,   132,    81,    88,    15,   114,   396,   361,   397,   398,
     248,  -297,   342,   249,  -209,   248,   335,   336,   249,   107,
     237,   238,   362,    16,  -297,   172,   173,   108,    17,    18,
      19,    81,   373,    15,   109,    81,  -297,    15,   308,   124,
      81,    89,    15,   174,    97,   175,   108,   494,    20,    14,
     381,    15,   183,   109,   404,   100,   184,   382,   308,   310,
     127,    16,   311,   185,   309,  -251,    17,    18,    19,   488,
      16,   487,   121,   308,   417,    17,    18,    19,   322,   310,
     129,   263,   311,   425,   405,   264,    20,   323,   324,   265,
     133,  -297,   446,   266,   310,    20,   137,   311,   447,   325,
     308,   326,   435,   306,   223,   457,   418,   186,   187,   188,
     189,   458,   454,   312,  -264,  -264,   423,   424,   313,   314,
     315,   310,   134,    99,   311,   135,   103,   343,   344,   136,
     115,   445,   140,   395,   119,   120,   159,   122,   123,   345,
     112,   346,    15,   160,   456,   400,   401,   179,   316,   317,
     318,   463,   464,   508,   509,   466,   180,   181,   191,   192,
     196,   459,     5,   194,    -6,   201,    -6,    -6,    -6,   462,
      -6,    -6,    -6,   467,   206,   210,  -120,   214,   219,    -6,
      -6,   218,    -6,   -11,   -11,    -6,   -11,   220,   -13,   235,
      -6,    -6,    -6,   239,   241,   -11,   258,   -11,   259,   260,
     267,   -11,   -11,   -11,   268,   -11,   -11,   269,   -11,   274,
      -6,    -6,   -11,   -11,   275,   276,   501,   277,   279,   283,
     288,   284,   285,   -11,   -11,   287,   169,   286,   -11,   -11,
     -11,   -11,   -11,   296,   502,   297,   298,   300,   303,   340,
     348,   352,   504,     6,   354,    -6,    -6,   363,    -6,    -6,
     366,   364,     7,   -11,   -11,   350,   -11,   -11,   368,   374,
     369,   378,   383,    30,    31,  -233,    32,   379,   385,   409,
     399,   411,   414,   461,   419,    33,    -6,    34,   420,    -6,
      -6,    35,    36,    37,   426,    38,    39,   431,    40,   422,
     432,   438,    41,    42,   439,   441,   311,   451,   471,   468,
     460,   474,   479,    43,    44,   489,   481,   476,    45,    46,
      47,    48,    49,   480,   493,   483,   485,   498,   503,   495,
     507,   505,    98,   138,   182,   497,   500,   492,   433,   410,
     200,   242,   339,    50,    51,   413,   496,    52,   484,   499,
     406,   338,   486,   402,   475,   473,   307,     0,   491
};

static const yytype_int16 yycheck[] =
{
      98,    99,   117,   115,    87,   103,   235,   149,   240,   276,
     251,   243,   303,    14,     3,     4,     3,     4,     7,     9,
       3,     0,    11,     3,   286,    14,    15,    16,   126,    16,
     321,     7,    16,   126,    42,    11,    26,     7,    16,    14,
     299,    11,    32,   275,     0,    35,     7,    16,    38,    39,
      11,    20,    14,    56,    76,    14,    46,    47,   111,    49,
     113,     6,     7,    14,     9,    10,    11,   237,   238,     9,
      10,    12,    13,    18,    19,    14,     7,    76,    18,    19,
      58,   193,    42,   225,   226,    45,   108,     3,   101,   201,
      80,     3,    81,    96,   353,    11,    56,    57,    27,    83,
     113,   106,    14,    65,    66,    64,    96,   374,    84,   108,
     111,   101,   111,   121,    84,    85,   378,    87,   105,   386,
     213,    37,   107,    84,    85,    76,    87,   106,   115,   109,
     113,   100,   399,   122,    21,   122,   111,    32,    33,    84,
      85,    72,    87,   119,   411,   115,    97,   137,    79,   119,
       7,   442,   101,   106,    11,   104,    26,   108,   119,   329,
     330,   102,     3,    50,     5,     7,   308,   279,   310,   311,
     115,    26,   413,   118,   119,   115,   107,   108,   118,    95,
      48,    49,   280,    24,    39,    89,    90,   103,    29,    30,
      31,     3,   290,     5,   110,     3,    51,     5,    80,    11,
       3,    88,     5,   107,     7,   109,   103,   474,    49,     3,
       7,     5,     3,   110,   326,    11,     7,    14,    80,   101,
       3,    24,   104,    14,    86,    16,    29,    30,    31,   461,
      24,   460,    14,    80,   346,    29,    30,    31,   120,   101,
       3,     7,   104,   358,   327,    11,    49,   102,   103,     7,
      14,   106,   111,    11,   101,    49,    97,   104,   117,   114,
      80,   116,   377,     3,     4,   111,   349,    58,    59,    60,
      61,   117,   414,    72,    73,    74,    56,    57,    77,    78,
      79,   101,   106,    30,   104,    95,    33,   102,   103,   106,
      37,   403,    46,   113,    41,    42,     3,    44,    45,   114,
       3,   116,     5,    51,   416,    73,    74,     3,   107,   108,
     109,    53,    54,    98,    99,   430,     3,     3,    16,    14,
      56,   419,     1,    16,     3,    38,     5,     6,     7,   427,
       9,    10,    11,   431,    16,    65,    16,     7,    38,    18,
      19,    17,    21,    22,    23,    24,    25,     8,    27,    39,
      29,    30,    31,   105,   113,    34,     7,    36,     7,    14,
       7,    40,    41,    42,   115,    44,    45,     3,    47,   119,
      49,    50,    51,    52,     3,     3,   491,    14,    42,    83,
      53,    16,    16,    62,    63,    16,    97,    14,    67,    68,
      69,    70,    71,   115,   492,     7,    38,    14,   115,   101,
       6,    14,   500,    82,    14,    84,    85,    74,    87,    88,
      14,   112,    91,    92,    93,    38,    95,    96,    16,     3,
     109,    14,     3,    22,    23,     3,    25,   115,     3,     3,
       7,     7,    75,     3,    14,    34,   115,    36,   113,   118,
     119,    40,    41,    42,    26,    44,    45,    14,    47,   120,
       7,   113,    51,    52,   113,    11,   104,   113,     7,    14,
      39,   113,    14,    62,    63,    55,   105,   117,    67,    68,
      69,    70,    71,   113,     7,   108,   117,    14,     3,   117,
      16,    97,    30,    80,   110,   117,   488,   468,   372,   333,
     120,   160,   248,    92,    93,   340,   482,    96,   455,   487,
     328,   244,   457,   321,   442,   439,   227,    -1,   466
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint16 yystos[] =
{
       0,   106,   124,   125,     0,     1,    82,    91,   126,   127,
     129,   130,   136,    14,     3,     5,    24,    29,    30,    31,
      49,   131,   132,   133,   134,   147,   148,   149,   228,   229,
      22,    23,    25,    34,    36,    40,    41,    42,    44,    45,
      47,    51,    52,    62,    63,    67,    68,    69,    70,    71,
      92,    93,    96,   135,   137,   138,   143,   151,   153,   154,
     156,   161,   162,   163,   166,   173,   179,   181,   183,   184,
     185,   186,   249,   253,   254,   255,   257,    27,   106,   128,
     107,     3,   227,   228,   106,    26,   224,    21,    50,    88,
     189,   190,   191,   192,   193,   245,   228,     7,   147,   227,
      11,   152,   228,   227,     3,    11,    37,    95,   103,   110,
     144,   145,     3,   228,     7,   227,   228,   228,   258,   227,
     227,    14,   227,   227,    11,   228,   228,     3,   228,     3,
      12,    13,   102,    14,   106,    95,   106,    97,   132,   146,
      46,    32,    33,   225,     7,    11,    84,    85,    87,   119,
     210,   211,   212,   246,   115,   210,   241,   242,   244,     3,
      51,   195,   196,   228,    76,   108,   139,   140,   139,    97,
     139,   261,    89,    90,   107,   109,   155,   228,   139,     3,
       3,     3,   145,     3,     7,    14,    58,    59,    60,    61,
     160,    16,    14,   222,    16,   264,    56,    96,   174,   175,
     175,    38,    14,    65,    66,   182,    16,   259,   139,   261,
      65,   188,   167,   168,     7,     3,   109,   228,    17,    38,
       8,   213,     3,     4,    16,   105,   115,   122,   230,   231,
     234,   237,   240,   230,   247,    39,   218,    48,    49,   105,
     194,   113,   174,   215,     9,    10,    18,    19,   115,   118,
     199,   200,   201,   202,   204,   205,   207,   208,     7,     7,
      14,     7,    11,     7,    11,     7,    11,     7,   115,     3,
      16,    20,   100,   222,   119,     3,     3,    14,   221,    42,
     219,    64,   222,    83,    16,    16,    14,    16,    53,   187,
     169,   170,   171,   264,   261,   250,   115,     7,    38,   150,
      14,   220,   226,   115,   230,   230,     3,   240,    80,    86,
     101,   104,    72,    77,    78,    79,   107,   108,   109,   232,
     233,   236,   120,   102,   103,   114,   116,   248,   225,   241,
     241,   219,   121,   219,    79,   107,   108,   206,   232,   199,
     101,   197,   201,   102,   103,   114,   116,   203,     6,   209,
      38,   262,    14,   157,    14,    45,    56,    57,   164,   219,
     213,   222,   139,    74,   112,   180,    14,   262,    16,   109,
     141,   142,   111,   139,     3,   172,    83,   263,    14,   115,
     251,     7,    14,     3,   211,     3,     7,    11,    14,    15,
      81,   234,   235,   238,   239,   113,   230,   230,   230,     7,
      73,    74,   235,   111,   222,   210,   226,   241,   241,     3,
     194,     7,   113,   200,    75,   214,   111,   222,   210,    14,
     113,   211,   120,    56,    57,   264,    26,   165,    58,   176,
     260,    14,     7,   170,   213,   264,   262,   252,   113,   113,
     213,    11,   111,   113,   213,   222,   111,   117,     3,    14,
     243,   113,   216,   213,   230,   198,   222,   111,   117,   139,
      39,     3,   139,    53,    54,   177,   264,   139,    14,     3,
     113,     7,   256,   256,   113,   239,   117,   222,   223,    14,
     113,   105,   217,   108,   218,   117,   223,   225,   219,    55,
     178,   263,   169,     7,   213,   117,   195,   117,    14,   220,
     165,   264,   139,     3,   139,    97,   158,    16,    98,    99,
     159
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     Rprintf ("%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF Rprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF ("%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF ("\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF ("token %s (", yytname[yytype]);
  else
    YYFPRINTF ("nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF ("Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (" %d", *bottom);
  YYFPRINTF ("\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF ("Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      Rprintf ("   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      Rprintf ("\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF (("Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF (("Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF (("Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF (("Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF (("Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 364 "parser.y"
    { prepare_parse(); ;}
    break;

  case 3:
#line 365 "parser.y"
    { if (generate_code)
                                           addHistoryLine();
                                         resetQueryBuffer();
                                         YYACCEPT; ;}
    break;

  case 4:
#line 369 "parser.y"
    { YYACCEPT; ;}
    break;

  case 5:
#line 370 "parser.y"
    { YYACCEPT; ;}
    break;

  case 6:
#line 373 "parser.y"
    { prepare_input(); ;}
    break;

  case 7:
#line 374 "parser.y"
    { after_CorpusCommand((yyvsp[(2) - (3)].cl)); ;}
    break;

  case 8:
#line 376 "parser.y"
    { 
                      if ((yyvsp[(2) - (2)].ival) == query_lock) {
                        query_lock = 0;
                      }
                      else {
                        Rprintf("ALERT! Query lock violation.\n");
                        Rprintf("\n"); /* so CQP.pm won't block -- should no longer be needed after switching to .EOL. mechanism */
                        exit(1);
                      }
                    ;}
    break;

  case 11:
#line 387 "parser.y"
    {if (query_lock) {warn_query_lock_violation(); YYABORT;} ;}
    break;

  case 12:
#line 388 "parser.y"
    { ;}
    break;

  case 13:
#line 389 "parser.y"
    {if (query_lock) {warn_query_lock_violation(); YYABORT;} ;}
    break;

  case 14:
#line 390 "parser.y"
    { exit_cqp++; ;}
    break;

  case 15:
#line 391 "parser.y"
    { if (yychar == YYEMPTY) yychar = yylex(); /* so synchronize works if lookahead yychar is empty */
                						   synchronize();
                						   /* in case of syntax errors, don't save history file */
                                           resetQueryBuffer();
                                           YYABORT; /* Oli did this:   yyerrok; */
                                           /* but we don't want to continue processing a line, when part of it failed */
                                         ;}
    break;

  case 16:
#line 400 "parser.y"
    { (yyval.cl) = (yyvsp[(1) - (1)].cl); ;}
    break;

  case 17:
#line 402 "parser.y"
    { (yyval.cl) = in_CorpusCommand((yyvsp[(1) - (3)].strval), (yyvsp[(3) - (3)].cl)); ;}
    break;

  case 18:
#line 403 "parser.y"
    { (yyval.cl) = in_CorpusCommand((yyvsp[(1) - (3)].strval), after_CorpusSetExpr((yyvsp[(3) - (3)].cl))); ;}
    break;

  case 19:
#line 407 "parser.y"
    { (yyval.cl) = in_UnnamedCorpusCommand((yyvsp[(1) - (2)].cl)); ;}
    break;

  case 20:
#line 410 "parser.y"
    { if (query_lock) {warn_query_lock_violation(); YYABORT;} (yyval.cl) = ActivateCorpus((yyvsp[(1) - (1)].cl)); ;}
    break;

  case 21:
#line 411 "parser.y"
    { (yyval.cl) = after_CorpusSetExpr((yyvsp[(1) - (1)].cl)); ;}
    break;

  case 22:
#line 412 "parser.y"
    { prepare_Query(); ;}
    break;

  case 23:
#line 413 "parser.y"
    { (yyval.cl) = after_Query((yyvsp[(2) - (2)].cl)); ;}
    break;

  case 47:
#line 443 "parser.y"
    { Rprintf("-::-EOL-::-\n"); fflush(stdout); ;}
    break;

  case 48:
#line 447 "parser.y"
    { do_cat((yyvsp[(2) - (3)].cl), &((yyvsp[(3) - (3)].redir)), 0, -1); cl_free((yyvsp[(3) - (3)].redir).name); ;}
    break;

  case 49:
#line 450 "parser.y"
    { do_cat((yyvsp[(2) - (7)].cl), &((yyvsp[(7) - (7)].redir)), (yyvsp[(4) - (7)].ival), (yyvsp[(6) - (7)].ival)); cl_free((yyvsp[(7) - (7)].redir).name); ;}
    break;

  case 50:
#line 453 "parser.y"
    { if (generate_code) 
                              do_cat((yyvsp[(2) - (3)].cl), &((yyvsp[(3) - (3)].redir)), 0, -1);
                            cl_free((yyvsp[(3) - (3)].redir).name);
                            drop_temp_corpora();
                          ;}
    break;

  case 51:
#line 459 "parser.y"
    { do_echo((yyvsp[(2) - (3)].strval), &((yyvsp[(3) - (3)].redir))); cl_free((yyvsp[(3) - (3)].redir).name); cl_free((yyvsp[(2) - (3)].strval)); ;}
    break;

  case 52:
#line 463 "parser.y"
    { do_save((yyvsp[(2) - (3)].cl), &((yyvsp[(3) - (3)].redir))); cl_free((yyvsp[(3) - (3)].redir).name); ;}
    break;

  case 54:
#line 467 "parser.y"
    { (yyval.redir).name = NULL; /* will open STDOUT */
                                          (yyval.redir).mode = "w";
                                          (yyval.redir).stream = NULL;
                                        ;}
    break;

  case 55:
#line 473 "parser.y"
    { (yyval.redir).name = (yyvsp[(2) - (2)].strval);
                                          (yyval.redir).mode = "w";
                                          (yyval.redir).stream = NULL;
                                        ;}
    break;

  case 56:
#line 477 "parser.y"
    { (yyval.redir).name = (yyvsp[(2) - (2)].strval);
                                          (yyval.redir).mode = "a";
                                          (yyval.redir).stream = NULL;
                                        ;}
    break;

  case 58:
#line 485 "parser.y"
    { (yyval.in_redir).name = NULL; /* will open STDIN */
                                          (yyval.in_redir).stream = NULL;
                                        ;}
    break;

  case 59:
#line 490 "parser.y"
    { (yyval.in_redir).name = (yyvsp[(2) - (2)].strval);
                                          (yyval.in_redir).stream = NULL;
                                        ;}
    break;

  case 60:
#line 495 "parser.y"
    { 
                                          show_corpora_files(UNDEF);
                                        ;}
    break;

  case 61:
#line 498 "parser.y"
    { 
                                          if (strncasecmp((yyvsp[(2) - (2)].strval), "var", 3) == 0) {
                                            do_PrintAllVariables();
                                          }
                                          else if ((strncasecmp((yyvsp[(2) - (2)].strval), "sys", 3) == 0) || (strncasecmp((yyvsp[(2) - (2)].strval), "corp", 4) == 0)) {
                                            show_corpora_files(SYSTEM);
                                          }
                                          else if ((strncasecmp((yyvsp[(2) - (2)].strval), "sub", 3) == 0) || (strcasecmp((yyvsp[(2) - (2)].strval), "named") == 0) || (strcasecmp((yyvsp[(2) - (2)].strval), "queries") == 0)) {
                                            show_corpora_files(SUB);    
                                          }
                                          else {
                                            cqpmessage(Error, "show what?");
                                          }
                                        ;}
    break;

  case 63:
#line 515 "parser.y"
    { PrintContextDescriptor(&CD); ;}
    break;

  case 66:
#line 525 "parser.y"
    { do_attribute_show((yyvsp[(2) - (2)].strval), 1); ;}
    break;

  case 67:
#line 526 "parser.y"
    { do_attribute_show((yyvsp[(2) - (2)].strval), 0); ;}
    break;

  case 68:
#line 529 "parser.y"
    { if (query_lock) {warn_query_lock_violation(); YYABORT;} (yyval.cl) = do_translate((yyvsp[(2) - (4)].cl), (yyvsp[(4) - (4)].strval)); ;}
    break;

  case 69:
#line 532 "parser.y"
    { if (query_lock) {warn_query_lock_violation(); YYABORT;} (yyval.cl) = do_setop((yyvsp[(1) - (3)].rngsetop), (yyvsp[(2) - (3)].cl), (yyvsp[(3) - (3)].cl)); ;}
    break;

  case 71:
#line 536 "parser.y"
    { (yyval.rngsetop) = RUnion; ;}
    break;

  case 72:
#line 537 "parser.y"
    { (yyval.rngsetop) = RUnion; ;}
    break;

  case 73:
#line 538 "parser.y"
    { (yyval.rngsetop) = RIntersection; ;}
    break;

  case 74:
#line 539 "parser.y"
    { (yyval.rngsetop) = RDiff; ;}
    break;

  case 75:
#line 545 "parser.y"
    { 
                                          do_start_timer();
                                          prepare_do_subset((yyvsp[(2) - (4)].cl), (yyvsp[(4) - (4)].field));  
                                          next_environment();   /* create environment for pattern compilation (will be freed in prepare_input() before next command) */
                                        ;}
    break;

  case 76:
#line 550 "parser.y"
    { 
                                          if (generate_code) {
                                            (yyval.cl) = do_subset((yyvsp[(4) - (6)].field), (yyvsp[(6) - (6)].boolt));
                                            do_timing("Subset computed");
                                          }
                                          else 
                                            (yyval.cl) = NULL;
                                        ;}
    break;

  case 77:
#line 561 "parser.y"
    {  ;}
    break;

  case 78:
#line 564 "parser.y"
    { if ((yyvsp[(2) - (2)].cl))
                                            dropcorpus((yyvsp[(2) - (2)].cl));
                                        ;}
    break;

  case 79:
#line 567 "parser.y"
    { if ((yyvsp[(1) - (1)].cl))
                                            dropcorpus((yyvsp[(1) - (1)].cl));
                                        ;}
    break;

  case 80:
#line 572 "parser.y"
    { do_PrintVariableValue((yyvsp[(2) - (2)].strval)); 
                                    free((yyvsp[(2) - (2)].strval));
                                  ;}
    break;

  case 81:
#line 579 "parser.y"
    { do_SetVariableValue((yyvsp[(2) - (3)].strval), 
                                                        (yyvsp[(3) - (3)].varsetting).operator, 
                                                        (yyvsp[(3) - (3)].varsetting).variableValue); 
                                    free((yyvsp[(2) - (3)].strval));
                                    free((yyvsp[(3) - (3)].varsetting).variableValue);
                                  ;}
    break;

  case 82:
#line 586 "parser.y"
    {
                                    do_AddSubVariables((yyvsp[(2) - (4)].strval), /*add*/1, (yyvsp[(4) - (4)].strval));
                                    free((yyvsp[(2) - (4)].strval));
                                    free((yyvsp[(4) - (4)].strval));
                                  ;}
    break;

  case 83:
#line 592 "parser.y"
    {
                                    do_AddSubVariables((yyvsp[(2) - (4)].strval), /*sub*/0, (yyvsp[(4) - (4)].strval));
                                    free((yyvsp[(2) - (4)].strval));
                                    free((yyvsp[(4) - (4)].strval));
                                  ;}
    break;

  case 84:
#line 598 "parser.y"
    {
                                    char *temp = cl_strdup("");
                                    do_SetVariableValue((yyvsp[(2) - (4)].strval), '=', temp);         /* cheap trick, this is :o) */
                                    free(temp);
                                    do_AddSubVariables((yyvsp[(2) - (4)].strval), /*add*/1, (yyvsp[(4) - (4)].strval));
                                    free((yyvsp[(2) - (4)].strval));
                                    free((yyvsp[(4) - (4)].strval));
                                  ;}
    break;

  case 85:
#line 609 "parser.y"
    { (yyval.varsetting).variableValue = (yyvsp[(2) - (2)].strval); (yyval.varsetting).operator = '<'; ;}
    break;

  case 86:
#line 610 "parser.y"
    { (yyval.varsetting).variableValue = (yyvsp[(2) - (2)].strval); (yyval.varsetting).operator = '='; ;}
    break;

  case 87:
#line 611 "parser.y"
    { (yyval.varsetting).variableValue = (yyvsp[(2) - (2)].strval); (yyval.varsetting).operator = '+'; ;}
    break;

  case 88:
#line 612 "parser.y"
    { (yyval.varsetting).variableValue = (yyvsp[(2) - (2)].strval); (yyval.varsetting).operator = '-'; ;}
    break;

  case 89:
#line 615 "parser.y"
    { char *msg;

                                          if ((yyvsp[(3) - (3)].varval).cval != NULL && (yyvsp[(3) - (3)].varval).ival >= 0) {
                                            msg = set_context_option_value((yyvsp[(2) - (3)].strval), (yyvsp[(3) - (3)].varval).cval, (yyvsp[(3) - (3)].varval).ival);
                                          }
                                          else if ((yyvsp[(3) - (3)].varval).cval != NULL) {
                                            /* get rid of quotes at start and end of value */
                                            /* -- removed because quotes should be stripped by lexer ({string} rule in parser.l) */
                                            /*
                                            if (($3.cval[0] == '"') && ($3.cval[strlen($3.cval)-1] == '"')
                                                || ($3.cval[0] == '\'') && ($3.cval[strlen($3.cval)-1] == '\'') ) {
                                              

                                              $3.cval[strlen($3.cval)-1] = '\0';
                                              $3.cval = $3.cval + 1;
                                            }
                                            */
                                            msg = set_string_option_value((yyvsp[(2) - (3)].strval), (yyvsp[(3) - (3)].varval).cval);
                                          }
                                          else
                                            msg = set_integer_option_value((yyvsp[(2) - (3)].strval), (yyvsp[(3) - (3)].varval).ival);

                                          if (msg != NULL)
                                            cqpmessage(Warning,
                                                       "Option set error:\n%s", msg);
                                        ;}
    break;

  case 90:
#line 641 "parser.y"
    { int opt;

                                          if ((opt = find_option((yyvsp[(2) - (2)].strval))) >= 0)
                                            print_option_value(opt);
                                          else
                                            cqpmessage(Warning,
                                                     "Unknown option: ``%s''\n", (yyvsp[(2) - (2)].strval));
                                        ;}
    break;

  case 91:
#line 649 "parser.y"
    {
                                          print_option_values();
                                        ;}
    break;

  case 92:
#line 652 "parser.y"
    { do_set_target((yyvsp[(2) - (4)].cl), (yyvsp[(3) - (4)].field), (yyvsp[(4) - (4)].field)); ;}
    break;

  case 93:
#line 653 "parser.y"
    { do_set_target((yyvsp[(2) - (4)].cl), (yyvsp[(3) - (4)].field), NoField); ;}
    break;

  case 94:
#line 654 "parser.y"
    {
                                          if (generate_code) {
                                            old_query_corpus = query_corpus;
                                            query_corpus = (yyvsp[(2) - (4)].cl);  /* set query_corpus for compiling the ExtConstraint pattern */
                                            next_environment(); /* create environment for pattern compilation (will be freed in prepare_input() before next command) */
                                            do_start_timer();
                                          }
                                        ;}
    break;

  case 95:
#line 664 "parser.y"
    {
                                          do_set_complex_target((yyvsp[(2) - (11)].cl), (yyvsp[(3) - (11)].field), (yyvsp[(4) - (11)].search_strategy), (yyvsp[(6) - (11)].boolt), (yyvsp[(8) - (11)].direction), (yyvsp[(9) - (11)].ival), (yyvsp[(10) - (11)].strval), (yyvsp[(11) - (11)].base).field, (yyvsp[(11) - (11)].base).inclusive);
                                          if (generate_code) 
                                            do_timing("``set target ...'' completed");
                                        ;}
    break;

  case 96:
#line 671 "parser.y"
    {
                                          /* from (match|keyword|target) [inclusive|exclusive] */
                                          (yyval.base).field = (yyvsp[(2) - (3)].field);
                                          (yyval.base).inclusive = (yyvsp[(3) - (3)].ival);
                                        ;}
    break;

  case 97:
#line 676 "parser.y"
    { 
                                          (yyval.base).field = MatchField;
                                          (yyval.base).inclusive = 0;
                                        ;}
    break;

  case 98:
#line 682 "parser.y"
    { (yyval.ival) = 1; ;}
    break;

  case 99:
#line 683 "parser.y"
    { (yyval.ival) = 0; ;}
    break;

  case 100:
#line 684 "parser.y"
    { (yyval.ival) = 0; /* default is exclusive */ ;}
    break;

  case 101:
#line 687 "parser.y"
    { (yyval.varval).ival = (yyvsp[(1) - (1)].ival);
                                          (yyval.varval).cval = NULL;
                                          (yyval.varval).ok = 1;
                                        ;}
    break;

  case 102:
#line 691 "parser.y"
    { (yyval.varval).ival = -1;
                                          (yyval.varval).cval = (yyvsp[(1) - (1)].strval);
                                          (yyval.varval).ok = 1;
                                        ;}
    break;

  case 103:
#line 695 "parser.y"
    {
                                          (yyval.varval).ival = -1;
                                          (yyval.varval).cval = (yyvsp[(1) - (1)].strval);
                                          (yyval.varval).ok = 1;
                                        ;}
    break;

  case 104:
#line 700 "parser.y"
    { (yyval.varval).ival = 1; (yyval.varval).cval = NULL; (yyval.varval).ok = 1; ;}
    break;

  case 105:
#line 701 "parser.y"
    { (yyval.varval).ival = 1; (yyval.varval).cval = NULL; (yyval.varval).ok = 1; ;}
    break;

  case 106:
#line 702 "parser.y"
    { (yyval.varval).ival = 0; (yyval.varval).cval = NULL; (yyval.varval).ok = 1; ;}
    break;

  case 107:
#line 703 "parser.y"
    { (yyval.varval).ival = 0; (yyval.varval).cval = NULL; (yyval.varval).ok = 1; ;}
    break;

  case 108:
#line 704 "parser.y"
    {
                                          (yyval.varval).ival = (yyvsp[(1) - (2)].ival);
                                          (yyval.varval).cval = (yyvsp[(2) - (2)].strval);
                                          (yyval.varval).ok = 1;
                                        ;}
    break;

  case 109:
#line 711 "parser.y"
    { do_exec((yyvsp[(2) - (2)].strval)); ;}
    break;

  case 110:
#line 714 "parser.y"
    { do_info((yyvsp[(2) - (2)].cl)); ;}
    break;

  case 111:
#line 715 "parser.y"
    { do_info(current_corpus); ;}
    break;

  case 112:
#line 720 "parser.y"
    { 
                                  do_group((yyvsp[(2) - (10)].cl), (yyvsp[(3) - (10)].Anchor).anchor, (yyvsp[(3) - (10)].Anchor).offset, (yyvsp[(4) - (10)].strval), (yyvsp[(6) - (10)].Anchor).anchor, (yyvsp[(6) - (10)].Anchor).offset, (yyvsp[(7) - (10)].strval), (yyvsp[(8) - (10)].ival), (yyvsp[(9) - (10)].ival), (yyvsp[(5) - (10)].ival), &((yyvsp[(10) - (10)].redir))); 
                                  cl_free((yyvsp[(10) - (10)].redir).name);
                                ;}
    break;

  case 113:
#line 726 "parser.y"
    { 
                                  do_group2((yyvsp[(2) - (7)].cl), (yyvsp[(3) - (7)].Anchor).anchor, (yyvsp[(3) - (7)].Anchor).offset, (yyvsp[(4) - (7)].strval), (yyvsp[(5) - (7)].ival), (yyvsp[(6) - (7)].ival), &((yyvsp[(7) - (7)].redir)));
                                  cl_free((yyvsp[(7) - (7)].redir).name);
                                ;}
    break;

  case 114:
#line 732 "parser.y"
    { (yyval.ival) = 0; ;}
    break;

  case 115:
#line 733 "parser.y"
    { (yyval.ival) = 0; ;}
    break;

  case 116:
#line 734 "parser.y"
    { (yyval.ival) = 1; ;}
    break;

  case 117:
#line 735 "parser.y"
    { (yyval.ival) = 1; ;}
    break;

  case 118:
#line 740 "parser.y"
    { (yyval.ival) = 1; ;}
    break;

  case 119:
#line 741 "parser.y"
    { (yyval.ival) = 0; ;}
    break;

  case 120:
#line 749 "parser.y"
    { free_tabulation_list(); ;}
    break;

  case 121:
#line 751 "parser.y"
    { print_tabulation((yyvsp[(2) - (5)].cl), 0, INT_MAX, &((yyvsp[(5) - (5)].redir))); cl_free((yyvsp[(5) - (5)].redir).name); ;}
    break;

  case 122:
#line 753 "parser.y"
    { free_tabulation_list(); ;}
    break;

  case 123:
#line 755 "parser.y"
    { print_tabulation((yyvsp[(2) - (9)].cl), (yyvsp[(5) - (9)].ival), (yyvsp[(7) - (9)].ival), &((yyvsp[(9) - (9)].redir))); cl_free((yyvsp[(9) - (9)].redir).name); ;}
    break;

  case 124:
#line 759 "parser.y"
    { append_tabulation_item((yyvsp[(1) - (1)].tabulation_item)); ;}
    break;

  case 125:
#line 761 "parser.y"
    { append_tabulation_item((yyvsp[(3) - (3)].tabulation_item)); ;}
    break;

  case 126:
#line 765 "parser.y"
    {
                     (yyval.tabulation_item) = new_tabulation_item();
                     (yyval.tabulation_item)->attribute_name = (yyvsp[(2) - (2)].AttributeSpecification).name;
                     (yyval.tabulation_item)->flags   = (yyvsp[(2) - (2)].AttributeSpecification).flags;
                     (yyval.tabulation_item)->anchor1 = (yyvsp[(1) - (2)].AnchorPair).anchor1;
                     (yyval.tabulation_item)->offset1 = (yyvsp[(1) - (2)].AnchorPair).offset1;
                     (yyval.tabulation_item)->anchor2 = (yyvsp[(1) - (2)].AnchorPair).anchor2;
                     (yyval.tabulation_item)->offset2 = (yyvsp[(1) - (2)].AnchorPair).offset2;
                   ;}
    break;

  case 127:
#line 777 "parser.y"
    { (yyval.AnchorPair).anchor1 = (yyval.AnchorPair).anchor2 = (yyvsp[(1) - (1)].Anchor).anchor; (yyval.AnchorPair).offset1 = (yyval.AnchorPair).offset2 = (yyvsp[(1) - (1)].Anchor).offset; ;}
    break;

  case 128:
#line 779 "parser.y"
    { (yyval.AnchorPair).anchor1 = (yyvsp[(1) - (3)].Anchor).anchor; (yyval.AnchorPair).offset1 = (yyvsp[(1) - (3)].Anchor).offset;
                       (yyval.AnchorPair).anchor2 = (yyvsp[(3) - (3)].Anchor).anchor; (yyval.AnchorPair).offset2 = (yyvsp[(3) - (3)].Anchor).offset; ;}
    break;

  case 129:
#line 784 "parser.y"
    { (yyval.AttributeSpecification).name = (yyvsp[(1) - (2)].strval); (yyval.AttributeSpecification).flags = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 130:
#line 786 "parser.y"
    { (yyval.AttributeSpecification).name = NULL; (yyval.AttributeSpecification).flags = 0; ;}
    break;

  case 131:
#line 796 "parser.y"
    { 
                  int ok;
                  if ((yyvsp[(2) - (3)].cl) && generate_code) {
                    do_start_timer();
                    ok = SortSubcorpus((yyvsp[(2) - (3)].cl), (yyvsp[(3) - (3)].sortclause), 0, NULL);
                    FreeSortClause((yyvsp[(3) - (3)].sortclause));
                    do_timing("Query result sorted");
                    if (autoshow && ok && ((yyvsp[(2) - (3)].cl)->size > 0))
                      catalog_corpus((yyvsp[(2) - (3)].cl), NULL, 0, -1, GlobalPrintMode);
                  }
                ;}
    break;

  case 132:
#line 808 "parser.y"
    {
                  int ok;
                  if ((yyvsp[(2) - (4)].cl) && generate_code) {
                    do_start_timer();
                    ok = SortSubcorpusRandomize((yyvsp[(2) - (4)].cl), (yyvsp[(4) - (4)].ival));
                    do_timing("Query result randomized");
                    if (autoshow && ok && ((yyvsp[(2) - (4)].cl)->size > 0))
                      catalog_corpus((yyvsp[(2) - (4)].cl), NULL, 0, -1, GlobalPrintMode);
                  }
                ;}
    break;

  case 133:
#line 819 "parser.y"
    { 
                  int ok __attribute__((unused));
                  if ((yyvsp[(2) - (5)].cl) && generate_code) {
                    ok = SortSubcorpus((yyvsp[(2) - (5)].cl), (yyvsp[(3) - (5)].sortclause), ((yyvsp[(4) - (5)].ival) >= 1) ? (yyvsp[(4) - (5)].ival) : 1, &((yyvsp[(5) - (5)].redir)));
                    FreeSortClause((yyvsp[(3) - (5)].sortclause));
                    cl_free((yyvsp[(5) - (5)].redir).name);
                  }
                ;}
    break;

  case 134:
#line 830 "parser.y"
    { (yyval.sortclause) = (yyvsp[(1) - (1)].sortclause); ;}
    break;

  case 135:
#line 831 "parser.y"
    { (yyval.sortclause) = NULL; ;}
    break;

  case 136:
#line 835 "parser.y"
    {
                  if (generate_code) {
                    (yyval.sortclause) = cl_malloc(sizeof(SortClauseBuffer));
                    (yyval.sortclause)->attribute_name  = (yyvsp[(2) - (6)].strval);
                    (yyval.sortclause)->flags           = (yyvsp[(3) - (6)].ival);
                    (yyval.sortclause)->anchor1         = (yyvsp[(4) - (6)].AnchorPair).anchor1;
                    (yyval.sortclause)->offset1         = (yyvsp[(4) - (6)].AnchorPair).offset1;
                    (yyval.sortclause)->anchor2         = (yyvsp[(4) - (6)].AnchorPair).anchor2;
                    (yyval.sortclause)->offset2         = (yyvsp[(4) - (6)].AnchorPair).offset2;
                    (yyval.sortclause)->sort_ascending  = (yyvsp[(5) - (6)].ival);
                    (yyval.sortclause)->sort_reverse    = (yyvsp[(6) - (6)].ival);
                  }
                  else
                    (yyval.sortclause) = NULL;
                ;}
    break;

  case 137:
#line 852 "parser.y"
    { (yyval.AnchorPair).anchor1 = (yyval.AnchorPair).anchor2 = (yyvsp[(2) - (2)].Anchor).anchor; (yyval.AnchorPair).offset1 = (yyval.AnchorPair).offset2 = (yyvsp[(2) - (2)].Anchor).offset; ;}
    break;

  case 138:
#line 854 "parser.y"
    { (yyval.AnchorPair).anchor1 = (yyvsp[(2) - (4)].Anchor).anchor; (yyval.AnchorPair).offset1 = (yyvsp[(2) - (4)].Anchor).offset;
                              (yyval.AnchorPair).anchor2 = (yyvsp[(4) - (4)].Anchor).anchor; (yyval.AnchorPair).offset2 = (yyvsp[(4) - (4)].Anchor).offset; ;}
    break;

  case 139:
#line 857 "parser.y"
    { (yyval.AnchorPair).anchor1 = MatchField;    (yyval.AnchorPair).offset1 = 0;
                              (yyval.AnchorPair).anchor2 = MatchEndField; (yyval.AnchorPair).offset2 = 0; ;}
    break;

  case 140:
#line 861 "parser.y"
    { (yyval.ival) = 1; ;}
    break;

  case 141:
#line 862 "parser.y"
    { (yyval.ival) = 0; ;}
    break;

  case 142:
#line 863 "parser.y"
    { (yyval.ival) = 1; ;}
    break;

  case 143:
#line 866 "parser.y"
    { (yyval.ival) = 1; ;}
    break;

  case 144:
#line 867 "parser.y"
    { (yyval.ival) = 0; ;}
    break;

  case 145:
#line 873 "parser.y"
    {
                                          do_reduce((yyvsp[(2) - (5)].cl), (yyvsp[(4) - (5)].ival), (yyvsp[(5) - (5)].ival));
                                        ;}
    break;

  case 146:
#line 877 "parser.y"
    {
                                          RangeSetop((yyvsp[(2) - (5)].cl), RMaximalMatches, NULL, NULL);
                                        ;}
    break;

  case 147:
#line 881 "parser.y"
    {
                                          do_cut((yyvsp[(2) - (3)].cl), 0, (yyvsp[(3) - (3)].ival)-1);
                                        ;}
    break;

  case 148:
#line 885 "parser.y"
    {
                                          do_cut((yyvsp[(2) - (4)].cl), (yyvsp[(3) - (4)].ival), (yyvsp[(4) - (4)].ival));
                                        ;}
    break;

  case 149:
#line 890 "parser.y"
    { (yyval.ival) = 1; ;}
    break;

  case 150:
#line 891 "parser.y"
    { (yyval.ival) = 0; ;}
    break;

  case 151:
#line 895 "parser.y"
    {
                                          if ((yyvsp[(2) - (4)].cl) && generate_code) {
                                            do_delete_lines((yyvsp[(2) - (4)].cl), (yyvsp[(4) - (4)].field), SELECTED_LINES); 
                                          }
                                        ;}
    break;

  case 152:
#line 901 "parser.y"
    { 
                                          if ((yyvsp[(2) - (4)].cl) && generate_code) {
                                            do_delete_lines((yyvsp[(2) - (4)].cl), (yyvsp[(4) - (4)].field), UNSELECTED_LINES); 
                                          }
                                        ;}
    break;

  case 153:
#line 907 "parser.y"
    { 
                                          if ((yyvsp[(2) - (3)].cl) && generate_code) {
                                             do_delete_lines_num((yyvsp[(2) - (3)].cl), (yyvsp[(3) - (3)].Distance).mindist, (yyvsp[(3) - (3)].Distance).maxdist);
                                           }
                                         ;}
    break;

  case 154:
#line 914 "parser.y"
    { (yyval.Distance).mindist = (yyvsp[(1) - (1)].ival); (yyval.Distance).maxdist = (yyvsp[(1) - (1)].ival); ;}
    break;

  case 155:
#line 915 "parser.y"
    { (yyval.Distance).mindist = (yyvsp[(1) - (3)].ival); (yyval.Distance).maxdist = (yyvsp[(3) - (3)].ival); ;}
    break;

  case 156:
#line 920 "parser.y"
    { do_sleep((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 157:
#line 924 "parser.y"
    { do_size((yyvsp[(2) - (3)].cl), (yyvsp[(3) - (3)].field)); ;}
    break;

  case 158:
#line 925 "parser.y"
    { do_printVariableSize((yyvsp[(2) - (2)].strval)); free((yyvsp[(2) - (2)].strval)); ;}
    break;

  case 159:
#line 930 "parser.y"
    { do_dump((yyvsp[(2) - (3)].cl), 0, INT_MAX, &((yyvsp[(3) - (3)].redir))); cl_free((yyvsp[(3) - (3)].redir).name); ;}
    break;

  case 160:
#line 932 "parser.y"
    { do_dump((yyvsp[(2) - (7)].cl), (yyvsp[(4) - (7)].ival), (yyvsp[(6) - (7)].ival), &((yyvsp[(7) - (7)].redir))); cl_free((yyvsp[(7) - (7)].redir).name); ;}
    break;

  case 161:
#line 936 "parser.y"
    { do_undump((yyvsp[(2) - (5)].strval), (yyvsp[(3) - (5)].ival), !(yyvsp[(4) - (5)].ival), &((yyvsp[(5) - (5)].in_redir))); cl_free((yyvsp[(5) - (5)].in_redir).name); ;}
    break;

  case 162:
#line 939 "parser.y"
    { (yyval.ival) = 1; ;}
    break;

  case 163:
#line 940 "parser.y"
    { (yyval.ival) = 0; ;}
    break;

  case 164:
#line 944 "parser.y"
    { (yyval.ival) = 0; ;}
    break;

  case 165:
#line 946 "parser.y"
    { 
            if ((yyvsp[(2) - (2)].field) == TargetField) { (yyval.ival) = 1; }
            else { yyerror("Invalid extension anchor in undump command"); YYABORT; } 
          ;}
    break;

  case 166:
#line 951 "parser.y"
    { 
            if ( (((yyvsp[(2) - (3)].field) == TargetField) && ((yyvsp[(3) - (3)].field) == KeywordField))
                 || (((yyvsp[(3) - (3)].field) == TargetField) && ((yyvsp[(2) - (3)].field) == KeywordField))
               ) { (yyval.ival) = 2; }
            else { yyerror("Invalid extension anchor in undump command"); YYABORT; } 
          ;}
    break;

  case 168:
#line 963 "parser.y"
    {
                  if ((yyvsp[(1) - (3)].cl) && (yyvsp[(3) - (3)].sortclause) && (yyvsp[(1) - (3)].cl)->size > 0) {
                    SortSubcorpus((yyvsp[(1) - (3)].cl), (yyvsp[(3) - (3)].sortclause), 0, NULL);
                    FreeSortClause((yyvsp[(3) - (3)].sortclause));
                  }
                  (yyval.cl) = (yyvsp[(1) - (3)].cl);
                ;}
    break;

  case 172:
#line 980 "parser.y"
    { (yyval.cl) = do_StandardQuery((yyvsp[(4) - (5)].ival), (yyvsp[(5) - (5)].ival), (yyvsp[(1) - (5)].strval)); ;}
    break;

  case 173:
#line 983 "parser.y"
    { (yyval.strval) = (yyvsp[(2) - (3)].strval); ;}
    break;

  case 174:
#line 984 "parser.y"
    { (yyval.strval) = NULL; ;}
    break;

  case 175:
#line 988 "parser.y"
    { (yyval.cl) = do_MUQuery((yyvsp[(2) - (4)].evalt), (yyvsp[(3) - (4)].ival), (yyvsp[(4) - (4)].ival)); ;}
    break;

  case 176:
#line 991 "parser.y"
    { (yyval.ival) = 1; ;}
    break;

  case 177:
#line 992 "parser.y"
    { (yyval.ival) = 0; ;}
    break;

  case 178:
#line 995 "parser.y"
    { if (generate_code) { CurEnv->match_label = labellookup(CurEnv->labels, "match", LAB_DEFINED, 1); } ;}
    break;

  case 179:
#line 996 "parser.y"
    { within_gc = 1; 
                                          if (generate_code) { CurEnv->matchend_label = labellookup(CurEnv->labels, "matchend", LAB_DEFINED, 1); } ;}
    break;

  case 180:
#line 998 "parser.y"
    { within_gc = 0; ;}
    break;

  case 181:
#line 999 "parser.y"
    { do_SearchPattern((yyvsp[(2) - (6)].evalt), (yyvsp[(4) - (6)].boolt)); ;}
    break;

  case 182:
#line 1003 "parser.y"
    { (yyval.evalt) = reg_disj((yyvsp[(1) - (3)].evalt), (yyvsp[(3) - (3)].evalt)); ;}
    break;

  case 183:
#line 1004 "parser.y"
    { (yyval.evalt) = (yyvsp[(1) - (1)].evalt); ;}
    break;

  case 184:
#line 1008 "parser.y"
    { (yyval.evalt) = reg_seq((yyvsp[(1) - (2)].evalt), (yyvsp[(2) - (2)].evalt)); ;}
    break;

  case 185:
#line 1009 "parser.y"
    { (yyval.evalt) = (yyvsp[(1) - (1)].evalt); ;}
    break;

  case 186:
#line 1012 "parser.y"
    { if (generate_code) {
                                           (yyval.evalt) = (yyvsp[(2) - (2)].evalt);
                                           (yyval.evalt)->node.left = (yyvsp[(1) - (2)].evalt);
                                           (yyval.evalt)->node.right = NULL;
                                         }
                                         else
                                           (yyval.evalt) = NULL;
                                       ;}
    break;

  case 187:
#line 1020 "parser.y"
    { if (generate_code)
                                           NEW_EVALNODE((yyval.evalt), re_repeat, (yyvsp[(1) - (2)].evalt), NULL, 0, repeat_inf);
                                         else
                                           (yyval.evalt) = NULL;
                                       ;}
    break;

  case 188:
#line 1025 "parser.y"
    { if (generate_code)
                                           NEW_EVALNODE((yyval.evalt), re_repeat, (yyvsp[(1) - (2)].evalt), NULL, 1, repeat_inf);
                                         else
                                           (yyval.evalt) = NULL;
                                       ;}
    break;

  case 189:
#line 1030 "parser.y"
    { if (generate_code)
                                           NEW_EVALNODE((yyval.evalt), re_repeat, (yyvsp[(1) - (2)].evalt), NULL, 0, 1);
                                         else
                                           (yyval.evalt) = NULL;
                                       ;}
    break;

  case 190:
#line 1035 "parser.y"
    { (yyval.evalt) = (yyvsp[(1) - (1)].evalt); ;}
    break;

  case 191:
#line 1038 "parser.y"
    { (yyval.evalt) = (yyvsp[(2) - (3)].evalt); ;}
    break;

  case 192:
#line 1039 "parser.y"
    { if (generate_code)
                                            NEW_EVALLEAF((yyval.evalt), (yyvsp[(1) - (1)].index));
                                          else
                                            (yyval.evalt) = NULL;
                                        ;}
    break;

  case 193:
#line 1044 "parser.y"
    { if (generate_code)
                                            NEW_EVALLEAF((yyval.evalt), (yyvsp[(1) - (1)].index));
                                          else
                                            (yyval.evalt) = NULL;
                                        ;}
    break;

  case 194:
#line 1049 "parser.y"
    { if (generate_code)
                                            NEW_EVALLEAF((yyval.evalt), (yyvsp[(1) - (1)].index));
                                          else
                                            (yyval.evalt) = NULL;
                                        ;}
    break;

  case 195:
#line 1056 "parser.y"
    { if (generate_code)
                                            NEW_EVALNODE((yyval.evalt), re_repeat, NULL, NULL, (yyvsp[(2) - (3)].ival), (yyvsp[(2) - (3)].ival));
                                          else
                                            (yyval.evalt) = NULL;
                                        ;}
    break;

  case 196:
#line 1061 "parser.y"
    { if (generate_code) /* new syntax for consistency with TAB queries */
                                            NEW_EVALNODE((yyval.evalt), re_repeat, NULL, NULL, 0, (yyvsp[(3) - (4)].ival));
                                          else
                                            (yyval.evalt) = NULL;
                                        ;}
    break;

  case 197:
#line 1067 "parser.y"
    { if ((yyvsp[(4) - (5)].ival) != repeat_inf && (yyvsp[(4) - (5)].ival) < (yyvsp[(2) - (5)].ival)) {
                                        	yyerror("invalid repetition range (maximum < minimum)");
                                        	YYERROR;
                                          }
                                          if (generate_code)
                                            NEW_EVALNODE((yyval.evalt), re_repeat, NULL, NULL, (yyvsp[(2) - (5)].ival), (yyvsp[(4) - (5)].ival));
                                          else
                                            (yyval.evalt) = NULL;
                                        ;}
    break;

  case 198:
#line 1079 "parser.y"
    { (yyval.index) = do_AnchorPoint((yyvsp[(1) - (1)].field), 0); ;}
    break;

  case 199:
#line 1080 "parser.y"
    { (yyval.index) = do_AnchorPoint((yyvsp[(1) - (1)].field), 1); ;}
    break;

  case 200:
#line 1083 "parser.y"
    { (yyval.index) = do_XMLTag((yyvsp[(1) - (2)].strval), 0, 0, NULL, 0); ;}
    break;

  case 201:
#line 1085 "parser.y"
    { (yyval.index) = do_XMLTag((yyvsp[(1) - (5)].strval), 0, (yyvsp[(2) - (5)].ival), (yyvsp[(3) - (5)].strval), (yyvsp[(4) - (5)].ival)); ;}
    break;

  case 202:
#line 1086 "parser.y"
    { (yyval.index) = do_XMLTag((yyvsp[(1) - (1)].strval), 1, 0, NULL, 0); ;}
    break;

  case 203:
#line 1089 "parser.y"
    { (yyval.ival) = (yyvsp[(1) - (1)].ival); ;}
    break;

  case 204:
#line 1090 "parser.y"
    { (yyval.ival) = OP_EQUAL | OP_NOT; ;}
    break;

  case 205:
#line 1091 "parser.y"
    { (yyval.ival) = OP_EQUAL; ;}
    break;

  case 206:
#line 1092 "parser.y"
    { (yyval.ival) = OP_EQUAL; ;}
    break;

  case 207:
#line 1097 "parser.y"
    { (yyval.index) = do_NamedWfPattern((yyvsp[(1) - (3)].ival), (yyvsp[(2) - (3)].strval), (yyvsp[(3) - (3)].index)); ;}
    break;

  case 208:
#line 1101 "parser.y"
    { (yyval.ival) = 1; ;}
    break;

  case 209:
#line 1102 "parser.y"
    { (yyval.ival) = 0; ;}
    break;

  case 210:
#line 1105 "parser.y"
    { (yyval.strval) = (yyvsp[(1) - (1)].strval); ;}
    break;

  case 211:
#line 1106 "parser.y"
    { (yyval.strval) = NULL; ;}
    break;

  case 212:
#line 1109 "parser.y"
    { (yyval.index) = do_WordformPattern((yyvsp[(1) - (1)].boolt), 0); ;}
    break;

  case 213:
#line 1110 "parser.y"
    { (yyval.index) = do_WordformPattern((yyvsp[(1) - (1)].boolt), 1); ;}
    break;

  case 214:
#line 1113 "parser.y"
    { (yyval.boolt) = do_StringConstraint((yyvsp[(1) - (2)].strval), (yyvsp[(2) - (2)].ival)); ;}
    break;

  case 215:
#line 1114 "parser.y"
    { (yyval.boolt) = NULL;
                                        if (!FindVariable((yyvsp[(1) - (1)].strval))) {
                                          cqpmessage(Error, 
                                                     "%s: no such variable", 
                                                     (yyvsp[(1) - (1)].strval));
                                          generate_code = 0;
                                        }
                                        else {
                                          (yyval.boolt) = do_SimpleVariableReference((yyvsp[(1) - (1)].strval));
                                        }
                                        free((yyvsp[(1) - (1)].strval));
                                      ;}
    break;

  case 216:
#line 1126 "parser.y"
    { (yyval.boolt) = (yyvsp[(2) - (3)].boolt); ;}
    break;

  case 217:
#line 1127 "parser.y"
    { if (generate_code) {
                                          NEW_BNODE((yyval.boolt));
                                          (yyval.boolt)->constnode.type = cnode;
                                          (yyval.boolt)->constnode.val  = 1;
                                        }
                                      ;}
    break;

  case 218:
#line 1135 "parser.y"
    { (yyval.boolt) = (yyvsp[(2) - (3)].boolt); ;}
    break;

  case 219:
#line 1136 "parser.y"
    { if (generate_code) {
                                          NEW_BNODE((yyval.boolt));
                                          (yyval.boolt)->constnode.type = cnode;
                                          (yyval.boolt)->constnode.val  = 1;
                                        }
                                      ;}
    break;

  case 220:
#line 1144 "parser.y"
    { int flags, i;
                                          flags = 0;

                                          for (i = 0; (yyvsp[(1) - (1)].strval)[i] != '\0'; i++) {
                                            switch ((yyvsp[(1) - (1)].strval)[i]) {
                                            case 'c':
                                              flags |= IGNORE_CASE;
                                              break;
                                            case 'd':
                                              flags |= IGNORE_DIAC;
                                              break;
                                            case 'l':
                                              flags = IGNORE_REGEX; /* literal */
                                              break;
                                            default:
                                              cqpmessage(Warning, "Unknown flag %s%c (ignored)", "%", (yyvsp[(1) - (1)].strval)[i]);
                                              break;
                                            }
                                          }

                                          /* %l supersedes all others */
                                          if (flags & IGNORE_REGEX) {
                                            if (flags != IGNORE_REGEX) {
                                              cqpmessage(Warning, "%s and %s flags cannot be combined with %s (ignored)",
                                                         "%c", "%d", "%l");
                                            }
                                            flags = IGNORE_REGEX;
                                          }

                                          (yyval.ival) = flags;
                                        ;}
    break;

  case 221:
#line 1175 "parser.y"
    { (yyval.ival) = 0; ;}
    break;

  case 222:
#line 1178 "parser.y"
    { (yyval.boolt) = (yyvsp[(2) - (2)].boolt); ;}
    break;

  case 223:
#line 1179 "parser.y"
    { (yyval.boolt) = NULL; ;}
    break;

  case 224:
#line 1184 "parser.y"
    { prepare_AlignmentConstraints((yyvsp[(3) - (3)].strval)); ;}
    break;

  case 225:
#line 1186 "parser.y"
    { if (generate_code)
                                           CurEnv->negated = (yyvsp[(5) - (6)].ival);
                                       ;}
    break;

  case 226:
#line 1189 "parser.y"
    { ;}
    break;

  case 227:
#line 1192 "parser.y"
    { (yyval.ival) = 1; ;}
    break;

  case 228:
#line 1193 "parser.y"
    { (yyval.ival) = 0; ;}
    break;

  case 229:
#line 1197 "parser.y"
    { if (generate_code) {
                    CurEnv->search_context.direction = (yyvsp[(2) - (3)].direction);
                    CurEnv->search_context.type = (yyvsp[(3) - (3)].context).type;
                    CurEnv->search_context.size = (yyvsp[(3) - (3)].context).size;
                    CurEnv->search_context.attrib = (yyvsp[(3) - (3)].context).attrib;
                  }
                ;}
    break;

  case 230:
#line 1204 "parser.y"
    { if (generate_code) {
                                            CurEnv->search_context.type  = word;
                                            CurEnv->search_context.size  = hard_boundary;
                                            CurEnv->search_context.attrib = NULL;
                                          }
                                        ;}
    break;

  case 231:
#line 1212 "parser.y"
    { (yyval.ival) = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 232:
#line 1213 "parser.y"
    { (yyval.ival) = 0; ;}
    break;

  case 233:
#line 1216 "parser.y"
    { (yyval.ival) = (yyvsp[(1) - (1)].ival); ;}
    break;

  case 234:
#line 1217 "parser.y"
    { (yyval.ival) = 1; ;}
    break;

  case 235:
#line 1220 "parser.y"
    { (yyval.ival) = (yyvsp[(1) - (1)].ival); ;}
    break;

  case 236:
#line 1221 "parser.y"
    { (yyval.ival) = 0; ;}
    break;

  case 237:
#line 1225 "parser.y"
    { if ((yyvsp[(1) - (1)].ival) < 0) {
											yyerror("expected a non-negative integer value");
											YYERROR;
										  }
										  (yyval.ival) = (yyvsp[(1) - (1)].ival);
										;}
    break;

  case 238:
#line 1233 "parser.y"
    { (yyval.ival) = (yyvsp[(1) - (1)].ival); ;}
    break;

  case 239:
#line 1234 "parser.y"
    { (yyval.ival) = repeat_inf; ;}
    break;

  case 240:
#line 1238 "parser.y"
    { expansion.direction = (yyvsp[(2) - (4)].direction);
                                          expansion.type = (yyvsp[(4) - (4)].context).type;
                                          expansion.size = (yyvsp[(4) - (4)].context).size;
                                          expansion.attrib = (yyvsp[(4) - (4)].context).attrib;
                                        ;}
    break;

  case 241:
#line 1243 "parser.y"
    { expansion.direction = leftright;
                                          expansion.type = word;
                                          expansion.size = 0;
                                          expansion.attrib = NULL;
                                        ;}
    break;

  case 242:
#line 1250 "parser.y"
    { (yyval.direction) = left; ;}
    break;

  case 243:
#line 1251 "parser.y"
    { (yyval.direction) = right; ;}
    break;

  case 244:
#line 1252 "parser.y"
    { (yyval.direction) = leftright; ;}
    break;

  case 245:
#line 1255 "parser.y"
    { do_Description(&((yyval.context)), (yyvsp[(1) - (2)].ival), (yyvsp[(2) - (2)].strval)); ;}
    break;

  case 246:
#line 1256 "parser.y"
    { do_Description(&((yyval.context)), (yyvsp[(1) - (1)].ival), NULL); ;}
    break;

  case 247:
#line 1259 "parser.y"
    { (yyval.cl) = (yyvsp[(1) - (1)].cl); ;}
    break;

  case 248:
#line 1260 "parser.y"
    { (yyval.cl) = findcorpus("Last", UNDEF, 0); ;}
    break;

  case 249:
#line 1264 "parser.y"
    { CorpusList *cl;

                                          cqpmessage(Message, "CID: %s", (yyvsp[(1) - (1)].strval));

                                          if ((cl = findcorpus((yyvsp[(1) - (1)].strval), UNDEF, 1)) == NULL) {
                                            cqpmessage(Error,
                                                       "Corpus ``%s'' is undefined", (yyvsp[(1) - (1)].strval));
                                            generate_code = 0;
                                            (yyval.cl) = NULL;
                                          }
                                          else if (!access_corpus(cl)) {
                                            cqpmessage(Warning,
                                                       "Corpus ``%s'' can't be accessed", (yyvsp[(1) - (1)].strval));
                                            (yyval.cl) = NULL;
                                          }
                                          else
                                            (yyval.cl) = cl;
                                        ;}
    break;

  case 250:
#line 1284 "parser.y"
    { (yyval.strval) = (yyvsp[(1) - (1)].strval); ;}
    break;

  case 251:
#line 1285 "parser.y"
    { (yyval.strval) = (yyvsp[(1) - (1)].strval); ;}
    break;

  case 252:
#line 1288 "parser.y"
    { (yyval.boolt) = bool_implies((yyvsp[(1) - (3)].boolt), (yyvsp[(3) - (3)].boolt)); ;}
    break;

  case 253:
#line 1289 "parser.y"
    { (yyval.boolt) = bool_or((yyvsp[(1) - (3)].boolt), (yyvsp[(3) - (3)].boolt)); ;}
    break;

  case 254:
#line 1290 "parser.y"
    { (yyval.boolt) = bool_and((yyvsp[(1) - (3)].boolt), (yyvsp[(3) - (3)].boolt)); ;}
    break;

  case 255:
#line 1291 "parser.y"
    { (yyval.boolt) = (yyvsp[(2) - (3)].boolt); ;}
    break;

  case 256:
#line 1292 "parser.y"
    { (yyval.boolt) = bool_not((yyvsp[(2) - (2)].boolt)); ;}
    break;

  case 257:
#line 1293 "parser.y"
    { (yyval.boolt) = (yyvsp[(1) - (1)].boolt); ;}
    break;

  case 258:
#line 1296 "parser.y"
    { (yyval.boolt) = do_RelExpr((yyvsp[(1) - (3)].boolt), (yyvsp[(2) - (3)].boolo), (yyvsp[(3) - (3)].boolt)); ;}
    break;

  case 259:
#line 1298 "parser.y"
    {
                 if ((yyvsp[(2) - (4)].ival) & OP_NOT) {
                   (yyval.boolt) = do_RelExpr((yyvsp[(1) - (4)].boolt), cmp_neq, do_mval_string((yyvsp[(3) - (4)].strval), (yyvsp[(2) - (4)].ival), (yyvsp[(4) - (4)].ival)));
                 }
                 else {
                   (yyval.boolt) = do_RelExpr((yyvsp[(1) - (4)].boolt), cmp_eq,  do_mval_string((yyvsp[(3) - (4)].strval), (yyvsp[(2) - (4)].ival), (yyvsp[(4) - (4)].ival)));
                 }
              ;}
    break;

  case 260:
#line 1306 "parser.y"
    { (yyval.boolt) = do_RelExExpr((yyvsp[(1) - (1)].boolt)); ;}
    break;

  case 261:
#line 1309 "parser.y"
    {(yyval.ival) = OP_CONTAINS | (yyvsp[(1) - (2)].ival);;}
    break;

  case 262:
#line 1310 "parser.y"
    {(yyval.ival) = OP_MATCHES  | (yyvsp[(1) - (2)].ival);;}
    break;

  case 263:
#line 1313 "parser.y"
    {(yyval.ival) = OP_NOT;;}
    break;

  case 264:
#line 1314 "parser.y"
    {(yyval.ival) = 0;;}
    break;

  case 265:
#line 1317 "parser.y"
    { (yyval.boolt) = do_LabelReference((yyvsp[(1) - (1)].strval), 0); ;}
    break;

  case 266:
#line 1318 "parser.y"
    { (yyval.boolt) = do_LabelReference((yyvsp[(2) - (2)].strval), 1); ;}
    break;

  case 267:
#line 1319 "parser.y"
    { (yyval.boolt) = do_IDReference((yyvsp[(1) - (1)].strval), 0); ;}
    break;

  case 268:
#line 1320 "parser.y"
    { (yyval.boolt) = do_IDReference((yyvsp[(2) - (2)].strval), 1); ;}
    break;

  case 269:
#line 1321 "parser.y"
    { (yyval.boolt) = do_IDReference(cl_strdup(field_type_to_name((yyvsp[(1) - (1)].field))), 0); ;}
    break;

  case 270:
#line 1322 "parser.y"
    { (yyval.boolt) = (yyvsp[(1) - (1)].boolt); ;}
    break;

  case 271:
#line 1325 "parser.y"
    { (yyval.boolt) = (yyvsp[(1) - (1)].boolt); ;}
    break;

  case 272:
#line 1326 "parser.y"
    { (yyval.boolt) = do_flagged_string((yyvsp[(1) - (2)].strval), (yyvsp[(2) - (2)].ival)); ;}
    break;

  case 273:
#line 1328 "parser.y"
    { 
                                          (yyval.boolt) = do_flagged_re_variable((yyvsp[(2) - (4)].strval), (yyvsp[(4) - (4)].ival)); 
                                        ;}
    break;

  case 274:
#line 1331 "parser.y"
    { if (generate_code) {
                                            if (!FindVariable((yyvsp[(1) - (1)].strval))) {
                                              cqpmessage(Error, 
                                                         "%s: no such variable", 
                                                         (yyvsp[(1) - (1)].strval));
                                              generate_code = 0;
                                              (yyval.boolt) = NULL;
                                            }
                                            else {
                                              NEW_BNODE((yyval.boolt));
                                              (yyval.boolt)->type = var_ref;
                                              (yyval.boolt)->varref.varName = (yyvsp[(1) - (1)].strval);
                                            }
                                          }
                                          else
                                            (yyval.boolt) = NULL;
                                        ;}
    break;

  case 275:
#line 1348 "parser.y"
    { if (generate_code) {
                                            NEW_BNODE((yyval.boolt));
                                            (yyval.boolt)->type = int_leaf;
                                            (yyval.boolt)->leaf.ctype.iconst = (yyvsp[(1) - (1)].ival);
                                          }
                                          else
                                            (yyval.boolt) = NULL;
                                        ;}
    break;

  case 276:
#line 1356 "parser.y"
    { if (generate_code) {
                                            NEW_BNODE((yyval.boolt));
                                            (yyval.boolt)->type = float_leaf;
                                            (yyval.boolt)->leaf.ctype.fconst = (yyvsp[(1) - (1)].fval);
                                          }
                                          else
                                            (yyval.boolt) = NULL;
                                        ;}
    break;

  case 277:
#line 1366 "parser.y"
    { (yyval.boolo) = cmp_lt; ;}
    break;

  case 278:
#line 1367 "parser.y"
    { (yyval.boolo) = cmp_gt; ;}
    break;

  case 279:
#line 1368 "parser.y"
    { (yyval.boolo) = cmp_eq; ;}
    break;

  case 280:
#line 1369 "parser.y"
    { (yyval.boolo) = cmp_neq; ;}
    break;

  case 281:
#line 1370 "parser.y"
    { (yyval.boolo) = cmp_let; ;}
    break;

  case 282:
#line 1371 "parser.y"
    { (yyval.boolo) = cmp_get; ;}
    break;

  case 283:
#line 1374 "parser.y"
    { (yyval.boolt) = FunctionCall((yyvsp[(1) - (4)].strval), (yyvsp[(3) - (4)].apl)); ;}
    break;

  case 284:
#line 1377 "parser.y"
    { (yyval.apl) = (yyvsp[(1) - (1)].apl);
                                            ;}
    break;

  case 285:
#line 1380 "parser.y"
    { ActualParamList *last;

                                              if (generate_code) {
                                                assert((yyvsp[(1) - (3)].apl) != NULL);

                                                last = (yyvsp[(1) - (3)].apl);
                                                while (last->next != NULL)
                                                  last = last->next;
                                                last->next = (yyvsp[(3) - (3)].apl);
                                                (yyval.apl) = (yyvsp[(1) - (3)].apl);
                                              }
                                              else
                                                (yyval.apl) = NULL;
                                            ;}
    break;

  case 286:
#line 1396 "parser.y"
    { if (generate_code) {
                                             New((yyval.apl), ActualParamList);

                                             (yyval.apl)->param = (yyvsp[(1) - (1)].boolt);
                                             (yyval.apl)->next = NULL;
                                           }
                                           else
                                             (yyval.apl) = NULL;
                                         ;}
    break;

  case 287:
#line 1407 "parser.y"
    { (yyval.strval) = (yyvsp[(1) - (1)].strval); ;}
    break;

  case 288:
#line 1410 "parser.y"
    { (yyval.evalt) = (yyvsp[(1) - (1)].evalt); ;}
    break;

  case 289:
#line 1411 "parser.y"
    { (yyval.evalt) = (yyvsp[(1) - (1)].evalt); ;}
    break;

  case 290:
#line 1412 "parser.y"
    { if (generate_code) {
                                            NEW_EVALLEAF((yyval.evalt), (yyvsp[(1) - (1)].index));
                                          }
                                          else
                                            (yyval.evalt) = NULL;
                                        ;}
    break;

  case 291:
#line 1424 "parser.y"
    { (yyval.evalt) = do_MeetStatement((yyvsp[(3) - (6)].evalt), (yyvsp[(4) - (6)].evalt), &((yyvsp[(5) - (6)].context))); ;}
    break;

  case 292:
#line 1428 "parser.y"
    { (yyval.context).type = word;
                                          (yyval.context).size = (yyvsp[(1) - (2)].ival);
                                          (yyval.context).size2 = (yyvsp[(2) - (2)].ival);
                                          (yyval.context).attrib = NULL;
                                        ;}
    break;

  case 293:
#line 1433 "parser.y"
    { do_StructuralContext(&((yyval.context)), (yyvsp[(1) - (1)].strval)); ;}
    break;

  case 294:
#line 1434 "parser.y"
    { (yyval.context).type = word;
                                          (yyval.context).size = 1;
                                          (yyval.context).size2 = 1;
                                          (yyval.context).attrib = NULL;
                                        ;}
    break;

  case 295:
#line 1445 "parser.y"
    { (yyval.evalt) = do_UnionStatement((yyvsp[(3) - (5)].evalt), (yyvsp[(4) - (5)].evalt)); ;}
    break;

  case 296:
#line 1450 "parser.y"
    { (yyval.cl) = do_TABQuery((yyvsp[(2) - (3)].evalt)); ;}
    break;

  case 297:
#line 1455 "parser.y"
    { (yyval.evalt) = make_first_tabular_pattern((yyvsp[(1) - (2)].index), (yyvsp[(2) - (2)].evalt)); ;}
    break;

  case 298:
#line 1460 "parser.y"
    { (yyval.evalt) = add_tabular_pattern((yyvsp[(1) - (3)].evalt), &((yyvsp[(2) - (3)].context)), (yyvsp[(3) - (3)].index)); ;}
    break;

  case 299:
#line 1462 "parser.y"
    { (yyval.evalt) = NULL; ;}
    break;

  case 300:
#line 1465 "parser.y"
    { do_OptDistance(&((yyval.context)), (yyvsp[(2) - (3)].ival) + 1, (yyvsp[(2) - (3)].ival) + 1); ;}
    break;

  case 301:
#line 1467 "parser.y"
    { if ((yyvsp[(4) - (5)].ival) == repeat_inf)
                                        	do_OptDistance(&((yyval.context)), (yyvsp[(2) - (5)].ival) + 1, repeat_inf);
                                          else {
                                            if ((yyvsp[(4) - (5)].ival) < (yyvsp[(2) - (5)].ival)) {
                                        	  yyerror("invalid distance range (maximum < minimum)");
                                        	  YYERROR;
                                            }
                                            do_OptDistance(&((yyval.context)), (yyvsp[(2) - (5)].ival) + 1, (yyvsp[(4) - (5)].ival) + 1);
                                          } 
                                        ;}
    break;

  case 302:
#line 1477 "parser.y"
    { do_OptDistance(&((yyval.context)), 1, (yyvsp[(3) - (4)].ival) + 1); ;}
    break;

  case 303:
#line 1478 "parser.y"
    { do_OptDistance(&((yyval.context)), 1, repeat_inf); ;}
    break;

  case 304:
#line 1479 "parser.y"
    { do_OptDistance(&((yyval.context)), 2, repeat_inf); ;}
    break;

  case 305:
#line 1480 "parser.y"
    { do_OptDistance(&((yyval.context)), 1, 2); ;}
    break;

  case 306:
#line 1481 "parser.y"
    { do_OptDistance(&((yyval.context)), 1, 1); ;}
    break;

  case 307:
#line 1485 "parser.y"
    { add_user_to_list((yyvsp[(2) - (3)].strval), (yyvsp[(3) - (3)].strval)); ;}
    break;

  case 309:
#line 1487 "parser.y"
    { add_host_to_list((yyvsp[(2) - (2)].strval)); ;}
    break;

  case 310:
#line 1488 "parser.y"
    { add_hosts_in_subnet_to_list((yyvsp[(2) - (2)].strval)); ;}
    break;

  case 311:
#line 1489 "parser.y"
    { add_host_to_list(NULL); ;}
    break;

  case 314:
#line 1498 "parser.y"
    { add_grant_to_last_user((yyvsp[(2) - (2)].strval)); ;}
    break;

  case 316:
#line 1504 "parser.y"
    {
                                                if (enable_macros) 
                                                  define_macro((yyvsp[(3) - (7)].strval), (yyvsp[(5) - (7)].ival), NULL, (yyvsp[(7) - (7)].strval));  /* <macro.c> */
                                                else 
                                                  cqpmessage(Error, "CQP macros not enabled.");
                                                free((yyvsp[(7) - (7)].strval));  /* don't forget to free the allocated strings */
                                                free((yyvsp[(3) - (7)].strval));
                                                ;}
    break;

  case 317:
#line 1513 "parser.y"
    {
                                                if (enable_macros) 
                                                  define_macro((yyvsp[(3) - (7)].strval), 0, (yyvsp[(5) - (7)].strval), (yyvsp[(7) - (7)].strval));  /* <macro.c> */
                                                else 
                                                  cqpmessage(Error, "CQP macros not enabled.");
                                                free((yyvsp[(7) - (7)].strval));  /* don't forget to free the allocated strings */
                                                free((yyvsp[(5) - (7)].strval));
                                                free((yyvsp[(3) - (7)].strval));
                                                ;}
    break;

  case 318:
#line 1522 "parser.y"
    {
                                                load_macro_file((yyvsp[(4) - (4)].strval));
                                                free((yyvsp[(4) - (4)].strval));  /* don't forget to free the allocated string */
                                                ;}
    break;

  case 321:
#line 1533 "parser.y"
    {
                                          list_macros(NULL);
                                        ;}
    break;

  case 322:
#line 1536 "parser.y"
    {
                                          list_macros((yyvsp[(3) - (3)].strval));
                                          free((yyvsp[(3) - (3)].strval));
                                        ;}
    break;

  case 323:
#line 1541 "parser.y"
    {
                                          print_macro_definition((yyvsp[(3) - (6)].strval), (yyvsp[(5) - (6)].ival));
                                          free((yyvsp[(3) - (6)].strval));
                                        ;}
    break;

  case 324:
#line 1549 "parser.y"
    {
                                          int l1 = strlen((yyvsp[(1) - (2)].strval)), l2 = strlen((yyvsp[(2) - (2)].strval));
                                          char *s = (char *) cl_malloc(l1 + l2 + 2);
                                          strcpy(s, (yyvsp[(1) - (2)].strval)); s[l1] = ' ';
                                          strcpy(s+l1+1, (yyvsp[(2) - (2)].strval));
                                          s[l1+l2+1] = '\0';
                                          free((yyvsp[(1) - (2)].strval));
                                          free((yyvsp[(2) - (2)].strval));
                                          (yyval.strval) = s;
                                        ;}
    break;

  case 325:
#line 1559 "parser.y"
    { (yyval.strval) = (yyvsp[(1) - (1)].strval); ;}
    break;

  case 326:
#line 1563 "parser.y"
    { cl_randomize(); ;}
    break;

  case 327:
#line 1564 "parser.y"
    { cl_set_seed((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 329:
#line 1572 "parser.y"
    { (yyval.field) = (yyvsp[(1) - (1)].field); ;}
    break;

  case 330:
#line 1573 "parser.y"
    { (yyval.field) = NoField; ;}
    break;

  case 339:
#line 1589 "parser.y"
    { (yyval.Anchor).anchor = (yyvsp[(1) - (1)].field); (yyval.Anchor).offset = 0; ;}
    break;

  case 340:
#line 1590 "parser.y"
    { (yyval.Anchor).anchor = (yyvsp[(1) - (4)].field); (yyval.Anchor).offset = (yyvsp[(3) - (4)].ival); ;}
    break;


/* Line 1267 of yacc.c.  */
#line 4208 "parser.tab.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


#line 1593 "parser.y"



