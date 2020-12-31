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
     NAME_SYM = 258,
     ID_SYM = 259,
     INFO_SYM = 260,
     HOME_SYM = 261,
     ATTRIBUTE_SYM = 262,
     DIR_SYM = 263,
     CORPUS_SYM = 264,
     REVCORP_SYM = 265,
     REVCIDX_SYM = 266,
     FREQS_SYM = 267,
     LEXICON_SYM = 268,
     LEXIDX_SYM = 269,
     LEXSRT_SYM = 270,
     STRUCTURE_SYM = 271,
     ALIGNED_SYM = 272,
     DYNAMIC_SYM = 273,
     DOTS_SYM = 274,
     IGNORE_SYM = 275,
     ADMIN_SYM = 276,
     ACCESS_SYM = 277,
     USER_SYM = 278,
     GROUP_SYM = 279,
     ASSERT_SYM = 280,
     HOST_SYM = 281,
     PROPERTY_SYM = 282,
     IDENTIFIER = 283,
     STRING = 284,
     NUMBER = 285
   };
#endif
/* Tokens.  */
#define NAME_SYM 258
#define ID_SYM 259
#define INFO_SYM 260
#define HOME_SYM 261
#define ATTRIBUTE_SYM 262
#define DIR_SYM 263
#define CORPUS_SYM 264
#define REVCORP_SYM 265
#define REVCIDX_SYM 266
#define FREQS_SYM 267
#define LEXICON_SYM 268
#define LEXIDX_SYM 269
#define LEXSRT_SYM 270
#define STRUCTURE_SYM 271
#define ALIGNED_SYM 272
#define DYNAMIC_SYM 273
#define DOTS_SYM 274
#define IGNORE_SYM 275
#define ADMIN_SYM 276
#define ACCESS_SYM 277
#define USER_SYM 278
#define GROUP_SYM 279
#define ASSERT_SYM 280
#define HOST_SYM 281
#define PROPERTY_SYM 282
#define IDENTIFIER 283
#define STRING 284
#define NUMBER 285




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 69 "registry.y"
{
  char    *strval;
  int      ival;
  void    *args;
  void    *attr;

  IDList   idlist;

  struct {
    int status;
    char *path;
  } storage;
}
/* Line 1529 of yacc.c.  */
#line 123 "registry.tab.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE creglval;

