/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

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
/* Line 1529 of yacc.c.  */
#line 321 "parser.tab.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

