/* A Bison parser, made by GNU Bison 2.7.12-4996.  */

/* Bison interface for Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2013 Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

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

#ifndef YY_YY_CALC_BISON_H_INCLUDED
# define YY_YY_CALC_BISON_H_INCLUDED
/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     TKNUMl = 258,
     TKNUMf = 259,
     TKNUMd = 260,
     TKVAR = 261,
     TKNVAR = 262,
     TKIMAGE = 263,
     TKCOMMAND = 264,
     TKFUNC_d_d = 265,
     TKFUNC_dd_d = 266,
     TKFUNC_im_d = 267,
     NEG = 268
   };
#endif
/* Tokens.  */
#define TKNUMl 258
#define TKNUMf 259
#define TKNUMd 260
#define TKVAR 261
#define TKNVAR 262
#define TKIMAGE 263
#define TKCOMMAND 264
#define TKFUNC_d_d 265
#define TKFUNC_dd_d 266
#define TKFUNC_im_d 267
#define NEG 268



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{
/* Line 2053 of yacc.c  */
#line 21 "calc_bison.y"

  long     val_l;  /* long */  
  float    val_f;  /* float */
  double   val_d;  /* For returning numbers.     */
  char  *string;   /* For returning strings (variables, images)  */
  double (*fnctptr)();    /* pointer to function -> double */


/* Line 2053 of yacc.c  */
#line 92 "calc_bison.h"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE yylval;

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

#endif /* !YY_YY_CALC_BISON_H_INCLUDED  */
