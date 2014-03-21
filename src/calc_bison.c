/* A Bison parser, made by GNU Bison 2.7.12-4996.  */

/* Bison implementation for Yacc-like parsers in C
   
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
#define YYBISON_VERSION "2.7.12-4996"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
/* Line 371 of yacc.c  */
#line 1 "calc_bison.y"

#include <math.h>  /* For math functions, cos(), sin(), etc. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
  //#include "calc.h"  /* Contains definition of `symrec'        */

#include "Cfits.h"
#include "00CORE/00CORE.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_arith/COREMOD_arith.h"

extern DATA data;


 char calctmpimname[200];

  
/* Line 371 of yacc.c  */
#line 88 "calc_bison.c"

# ifndef YY_NULL
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULL nullptr
#  else
#   define YY_NULL 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "y.tab.h".  */
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
/* Line 387 of yacc.c  */
#line 21 "calc_bison.y"

  long     val_l;  /* long */  
  float    val_f;  /* float */
  double   val_d;  /* For returning numbers.     */
  char  *string;   /* For returning strings (variables, images)  */
  double (*fnctptr)();    /* pointer to function -> double */


/* Line 387 of yacc.c  */
#line 166 "calc_bison.c"
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

/* Copy the second part of user declarations.  */

/* Line 390 of yacc.c  */
#line 194 "calc_bison.c"

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
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef __attribute__
/* This feature is available in gcc versions 2.5 and later.  */
# if (! defined __GNUC__ || __GNUC__ < 2 \
      || (__GNUC__ == 2 && __GNUC_MINOR__ < 5))
#  define __attribute__(Spec) /* empty */
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif


/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(N) (N)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
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
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
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
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
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
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (YYID (0))
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   366

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  24
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  6
/* YYNRULES -- Number of rules.  */
#define YYNRULES  75
/* YYNRULES -- Number of states.  */
#define YYNSTATES  154

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   268

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      20,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      21,    22,    16,    15,    23,    14,     2,    17,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    13,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,    19,     2,     2,     2,     2,     2,
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
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    18
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     4,     7,     9,    12,    15,    18,    21,
      23,    27,    31,    35,    38,    42,    46,    48,    50,    54,
      58,    62,    66,    70,    74,    78,    82,    86,    90,    94,
      98,   102,   106,   110,   114,   118,   121,   125,   129,   133,
     137,   141,   147,   153,   159,   165,   169,   173,   175,   177,
     179,   183,   187,   191,   195,   199,   203,   207,   211,   215,
     219,   223,   227,   231,   235,   239,   243,   247,   251,   255,
     259,   263,   267,   271,   275,   279
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      25,     0,    -1,    -1,    25,    26,    -1,    20,    -1,    28,
      20,    -1,    27,    20,    -1,    29,    20,    -1,     1,    20,
      -1,     3,    -1,    27,    15,    27,    -1,    27,    14,    27,
      -1,    27,    16,    27,    -1,    14,    27,    -1,    27,    19,
      27,    -1,    21,    27,    22,    -1,     5,    -1,     6,    -1,
       6,    13,    27,    -1,     6,    13,    28,    -1,     7,    13,
      27,    -1,     7,    13,    28,    -1,    27,    15,    28,    -1,
      28,    15,    27,    -1,    28,    15,    28,    -1,    27,    14,
      28,    -1,    28,    14,    27,    -1,    28,    14,    28,    -1,
      27,    16,    28,    -1,    28,    16,    27,    -1,    28,    16,
      28,    -1,    27,    17,    27,    -1,    27,    17,    28,    -1,
      28,    17,    27,    -1,    28,    17,    28,    -1,    14,    28,
      -1,    27,    19,    28,    -1,    28,    19,    27,    -1,    28,
      19,    28,    -1,    10,    28,    22,    -1,    10,    27,    22,
      -1,    11,    28,    23,    28,    22,    -1,    11,    27,    23,
      28,    22,    -1,    11,    28,    23,    27,    22,    -1,    11,
      27,    23,    27,    22,    -1,    12,    29,    22,    -1,    21,
      28,    22,    -1,     7,    -1,     8,    -1,     9,    -1,     8,
      13,    29,    -1,     7,    13,    29,    -1,    29,    15,    29,
      -1,    29,    15,    28,    -1,    29,    15,    27,    -1,    28,
      15,    29,    -1,    27,    15,    29,    -1,    29,    14,    29,
      -1,    29,    14,    28,    -1,    29,    14,    27,    -1,    28,
      14,    29,    -1,    27,    14,    29,    -1,    29,    16,    29,
      -1,    29,    16,    28,    -1,    29,    16,    27,    -1,    28,
      16,    29,    -1,    27,    16,    29,    -1,    29,    17,    29,
      -1,    29,    17,    28,    -1,    29,    17,    27,    -1,    28,
      17,    29,    -1,    27,    17,    29,    -1,    29,    19,    27,
      -1,    29,    19,    28,    -1,    10,    29,    22,    -1,    21,
      29,    22,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint8 yyrline[] =
{
       0,    54,    54,    55,    59,    60,    65,    70,    74,    77,
      78,    79,    80,    81,    82,    83,    86,    87,    88,    89,
      90,    91,    92,    93,    94,    95,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   107,   108,   109,
     110,   111,   112,   113,   114,   115,   116,   120,   121,   122,
     123,   124,   125,   126,   127,   128,   129,   130,   131,   132,
     133,   134,   135,   136,   137,   138,   139,   140,   141,   142,
     143,   144,   145,   146,   147,   148
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "TKNUMl", "TKNUMf", "TKNUMd", "TKVAR",
  "TKNVAR", "TKIMAGE", "TKCOMMAND", "TKFUNC_d_d", "TKFUNC_dd_d",
  "TKFUNC_im_d", "'='", "'-'", "'+'", "'*'", "'/'", "NEG", "'^'", "'\\n'",
  "'('", "')'", "','", "$accept", "input", "line", "expl", "expd", "exps", YY_NULL
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,    61,    45,    43,    42,    47,   268,    94,
      10,    40,    41,    44
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    24,    25,    25,    26,    26,    26,    26,    26,    27,
      27,    27,    27,    27,    27,    27,    28,    28,    28,    28,
      28,    28,    28,    28,    28,    28,    28,    28,    28,    28,
      28,    28,    28,    28,    28,    28,    28,    28,    28,    28,
      28,    28,    28,    28,    28,    28,    28,    29,    29,    29,
      29,    29,    29,    29,    29,    29,    29,    29,    29,    29,
      29,    29,    29,    29,    29,    29,    29,    29,    29,    29,
      29,    29,    29,    29,    29,    29
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     1,     2,     2,     2,     2,     1,
       3,     3,     3,     2,     3,     3,     1,     1,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     2,     3,     3,     3,     3,
       3,     5,     5,     5,     5,     3,     3,     1,     1,     1,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3
};

/* YYDEFACT[STATE-NAME] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,     0,     9,    16,    17,    47,    48,    49,
       0,     0,     0,     0,     4,     0,     3,     0,     0,     0,
       8,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    13,    35,     0,     0,     0,
       0,     0,     0,     0,     0,     6,     0,     0,     0,     0,
       0,     5,     0,     0,     0,     0,     0,     7,    18,    19,
      20,    21,    51,    50,    40,    39,    74,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    45,    15,    46,    75,    11,    25,    61,    10,
      22,    56,    12,    28,    66,    31,    32,    71,    14,    36,
      26,    27,    60,    23,    24,    55,    29,    30,    65,    33,
      34,    70,    37,    38,    59,    58,    57,    54,    53,    52,
      64,    63,    62,    69,    68,    67,    72,    73,    20,    21,
      11,    25,    10,    22,    12,    28,    31,    32,     0,     0,
      26,    27,    23,    24,    29,    30,    33,    34,     0,     0,
      44,    42,    43,    41
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     1,    16,    32,    33,    19
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -20
static const yytype_int16 yypact[] =
{
     -20,   183,   -20,   -19,   -20,   -20,    -7,    -5,     0,   -20,
     202,    12,   202,    12,   -20,   202,   -20,   302,   309,   316,
     -20,    12,   202,   202,    42,    70,    96,     3,    12,    12,
     -12,    15,   323,   329,   129,     6,    16,   203,   212,   221,
     202,   202,   202,   202,    12,   -20,   202,   202,   202,   202,
      12,   -20,   202,   202,   202,   202,    12,   -20,   335,   341,
     323,   329,   347,   347,   -20,   -20,   -20,    12,   230,   239,
     248,   257,    12,    12,    12,    12,    12,    12,    12,    12,
      12,    12,   -20,   -20,   -20,   -20,    46,    51,    81,    46,
      51,    81,     6,    16,    17,     6,    16,    17,     6,    16,
      46,    51,    81,    46,    51,    81,     6,    16,    17,     6,
      16,    17,     6,    16,    46,    51,    81,    46,    51,    81,
       6,    16,    17,     6,    16,    17,     6,    16,   335,   341,
      87,   100,    87,   100,     6,    16,     6,    16,   266,   275,
      87,   100,    87,   100,     6,    16,     6,    16,   284,   293,
     -20,   -20,   -20,   -20
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -20,   -20,   -20,    -1,    80,   127
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
      17,    20,    72,    73,    74,    75,    21,    44,    22,    24,
      30,    76,    35,    23,    37,     4,    67,     5,     6,    27,
      58,    60,    28,    11,    12,    44,    13,    68,    70,    77,
      78,    79,    80,    29,    50,    50,    56,     0,    81,    86,
      89,    92,    95,    98,     0,   100,   103,   106,   109,   112,
       0,   114,   117,   120,   123,   126,    40,    41,    42,    43,
       0,    44,    42,    43,    64,    44,   128,    48,    49,     0,
      50,   130,   132,   134,   136,   138,   140,   142,   144,   146,
     148,    18,     0,     0,    46,    47,    48,    49,     0,    50,
      25,    31,    65,    36,     0,    38,     0,    54,    55,     0,
      56,    59,    61,    74,    75,     0,    44,     0,    69,    71,
      52,    53,    54,    55,     0,    56,    79,    80,    66,    50,
      87,    90,    93,    96,    99,     0,   101,   104,   107,   110,
     113,     0,   115,   118,   121,   124,   127,    26,     0,    34,
       0,     0,    39,    52,    53,    54,    55,   129,    56,    62,
      63,    82,   131,   133,   135,   137,   139,   141,   143,   145,
     147,   149,     0,     0,     0,     0,     0,    88,    91,    94,
      97,     0,     0,   102,   105,   108,   111,     0,     0,   116,
     119,   122,   125,     2,     3,     0,     4,     0,     5,     6,
       7,     8,     9,    10,    11,    12,     0,    13,     0,     0,
       0,     0,     0,    14,    15,     4,     0,     5,     6,     7,
       8,     9,    10,    11,    12,     0,    13,    40,    41,    42,
      43,     0,    44,    15,     0,    83,    46,    47,    48,    49,
       0,    50,     0,     0,    84,    52,    53,    54,    55,     0,
      56,     0,     0,    85,    72,    73,    74,    75,     0,    44,
       0,     0,    64,    77,    78,    79,    80,     0,    50,     0,
       0,    65,    72,    73,    74,    75,     0,    44,     0,     0,
      83,    77,    78,    79,    80,     0,    50,     0,     0,    84,
      72,    73,    74,    75,     0,    44,     0,     0,   150,    77,
      78,    79,    80,     0,    50,     0,     0,   151,    72,    73,
      74,    75,     0,    44,     0,     0,   152,    77,    78,    79,
      80,     0,    50,     0,     0,   153,    40,    41,    42,    43,
       0,    44,    45,    46,    47,    48,    49,     0,    50,    51,
      52,    53,    54,    55,     0,    56,    57,    40,    41,    42,
      43,     0,    44,    46,    47,    48,    49,     0,    50,    72,
      73,    74,    75,     0,    44,    77,    78,    79,    80,     0,
      50,    52,    53,    54,    55,     0,    56
};

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-20)))

#define yytable_value_is_error(Yytable_value) \
  YYID (0)

static const yytype_int8 yycheck[] =
{
       1,    20,    14,    15,    16,    17,    13,    19,    13,    10,
      11,    23,    13,    13,    15,     3,    13,     5,     6,     7,
      21,    22,    10,    11,    12,    19,    14,    28,    29,    14,
      15,    16,    17,    21,    19,    19,    19,    -1,    23,    40,
      41,    42,    43,    44,    -1,    46,    47,    48,    49,    50,
      -1,    52,    53,    54,    55,    56,    14,    15,    16,    17,
      -1,    19,    16,    17,    22,    19,    67,    16,    17,    -1,
      19,    72,    73,    74,    75,    76,    77,    78,    79,    80,
      81,     1,    -1,    -1,    14,    15,    16,    17,    -1,    19,
      10,    11,    22,    13,    -1,    15,    -1,    16,    17,    -1,
      19,    21,    22,    16,    17,    -1,    19,    -1,    28,    29,
      14,    15,    16,    17,    -1,    19,    16,    17,    22,    19,
      40,    41,    42,    43,    44,    -1,    46,    47,    48,    49,
      50,    -1,    52,    53,    54,    55,    56,    10,    -1,    12,
      -1,    -1,    15,    14,    15,    16,    17,    67,    19,    22,
      23,    22,    72,    73,    74,    75,    76,    77,    78,    79,
      80,    81,    -1,    -1,    -1,    -1,    -1,    40,    41,    42,
      43,    -1,    -1,    46,    47,    48,    49,    -1,    -1,    52,
      53,    54,    55,     0,     1,    -1,     3,    -1,     5,     6,
       7,     8,     9,    10,    11,    12,    -1,    14,    -1,    -1,
      -1,    -1,    -1,    20,    21,     3,    -1,     5,     6,     7,
       8,     9,    10,    11,    12,    -1,    14,    14,    15,    16,
      17,    -1,    19,    21,    -1,    22,    14,    15,    16,    17,
      -1,    19,    -1,    -1,    22,    14,    15,    16,    17,    -1,
      19,    -1,    -1,    22,    14,    15,    16,    17,    -1,    19,
      -1,    -1,    22,    14,    15,    16,    17,    -1,    19,    -1,
      -1,    22,    14,    15,    16,    17,    -1,    19,    -1,    -1,
      22,    14,    15,    16,    17,    -1,    19,    -1,    -1,    22,
      14,    15,    16,    17,    -1,    19,    -1,    -1,    22,    14,
      15,    16,    17,    -1,    19,    -1,    -1,    22,    14,    15,
      16,    17,    -1,    19,    -1,    -1,    22,    14,    15,    16,
      17,    -1,    19,    -1,    -1,    22,    14,    15,    16,    17,
      -1,    19,    20,    14,    15,    16,    17,    -1,    19,    20,
      14,    15,    16,    17,    -1,    19,    20,    14,    15,    16,
      17,    -1,    19,    14,    15,    16,    17,    -1,    19,    14,
      15,    16,    17,    -1,    19,    14,    15,    16,    17,    -1,
      19,    14,    15,    16,    17,    -1,    19
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    25,     0,     1,     3,     5,     6,     7,     8,     9,
      10,    11,    12,    14,    20,    21,    26,    27,    28,    29,
      20,    13,    13,    13,    27,    28,    29,     7,    10,    21,
      27,    28,    27,    28,    29,    27,    28,    27,    28,    29,
      14,    15,    16,    17,    19,    20,    14,    15,    16,    17,
      19,    20,    14,    15,    16,    17,    19,    20,    27,    28,
      27,    28,    29,    29,    22,    22,    22,    13,    27,    28,
      27,    28,    14,    15,    16,    17,    23,    14,    15,    16,
      17,    23,    22,    22,    22,    22,    27,    28,    29,    27,
      28,    29,    27,    28,    29,    27,    28,    29,    27,    28,
      27,    28,    29,    27,    28,    29,    27,    28,    29,    27,
      28,    29,    27,    28,    27,    28,    29,    27,    28,    29,
      27,    28,    29,    27,    28,    29,    27,    28,    27,    28,
      27,    28,    27,    28,    27,    28,    27,    28,    27,    28,
      27,    28,    27,    28,    27,    28,    27,    28,    27,    28,
      22,    22,    22,    22
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
   Once GCC version 2 has supplanted version 1, this can go.  However,
   YYFAIL appears to be in use.  Nevertheless, it is formally deprecated
   in Bison 2.4.2's NEWS entry, where a plan to phase it out is
   discussed.  */

#define YYFAIL		goto yyerrlab
#if defined YYFAIL
  /* This is here to suppress warnings from the GCC cpp's
     -Wunused-macros.  Normally we don't worry about that warning, but
     some users do, and we want to make it easy for users to remove
     YYFAIL uses, which will produce warnings from Bison 2.5.  */
#endif

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))

/* Error token number */
#define YYTERROR	1
#define YYERRCODE	256


/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
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
#  define YYFPRINTF fprintf
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
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
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
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  YYUSE (yytype);
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
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
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
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
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

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULL, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULL;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - Assume YYFAIL is not used.  It's too flawed to consider.  See
       <http://lists.gnu.org/archive/html/bison-patches/2009-12/msg00024.html>
       for details.  YYERROR is fine as it does not invoke this
       function.
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULL, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
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

  YYUSE (yytype);
}




/* The lookahead symbol.  */
int yychar;


#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval YY_INITIAL_VALUE(yyval_default);

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
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
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
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
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
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

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
        case 5:
/* Line 1787 of yacc.c  */
#line 60 "calc_bison.y"
    { 
printf("\t double: %.10g\n", (yyvsp[(1) - (2)].val_d)); 
data.cmdargtoken[data.cmdNBarg].type = 1; 
data.cmdargtoken[data.cmdNBarg].val.numf = (yyvsp[(1) - (2)].val_d);
}
    break;

  case 6:
/* Line 1787 of yacc.c  */
#line 65 "calc_bison.y"
    { 
printf("\t long:   %ld\n", (yyvsp[(1) - (2)].val_l)); 
data.cmdargtoken[data.cmdNBarg].type = 2; 
data.cmdargtoken[data.cmdNBarg].val.numl = (yyvsp[(1) - (2)].val_l);
}
    break;

  case 7:
/* Line 1787 of yacc.c  */
#line 70 "calc_bison.y"
    { if(data.Debug>0) {printf("\t string: %s\n", (yyvsp[(1) - (2)].string));}
    //data.cmdargtoken[data.cmdNBarg].type = 3;
sprintf(data.cmdargtoken[data.cmdNBarg].val.string, "%s", (yyvsp[(1) - (2)].string));
}
    break;

  case 8:
/* Line 1787 of yacc.c  */
#line 74 "calc_bison.y"
    { yyerrok;                  }
    break;

  case 9:
/* Line 1787 of yacc.c  */
#line 77 "calc_bison.y"
    { (yyval.val_l) = (yyvsp[(1) - (1)].val_l);        if(data.Debug>0){printf("this is a long\n");}}
    break;

  case 10:
/* Line 1787 of yacc.c  */
#line 78 "calc_bison.y"
    { (yyval.val_l) = (yyvsp[(1) - (3)].val_l) + (yyvsp[(3) - (3)].val_l);   if(data.Debug>0){printf("long + long\n");}}
    break;

  case 11:
/* Line 1787 of yacc.c  */
#line 79 "calc_bison.y"
    { (yyval.val_l) = (yyvsp[(1) - (3)].val_l) - (yyvsp[(3) - (3)].val_l);   if(data.Debug>0){printf("long - long\n");}}
    break;

  case 12:
/* Line 1787 of yacc.c  */
#line 80 "calc_bison.y"
    { (yyval.val_l) = (yyvsp[(1) - (3)].val_l) * (yyvsp[(3) - (3)].val_l);   if(data.Debug>0){printf("long * long\n");}}
    break;

  case 13:
/* Line 1787 of yacc.c  */
#line 81 "calc_bison.y"
    { (yyval.val_l) = -(yyvsp[(2) - (2)].val_l);       if(data.Debug>0){printf("-long\n");}}
    break;

  case 14:
/* Line 1787 of yacc.c  */
#line 82 "calc_bison.y"
    { (yyval.val_l) = (long) pow ((yyvsp[(1) - (3)].val_l),(yyvsp[(3) - (3)].val_l));  if(data.Debug>0){printf("long ^ long\n");}}
    break;

  case 15:
/* Line 1787 of yacc.c  */
#line 83 "calc_bison.y"
    { (yyval.val_l) = (yyvsp[(2) - (3)].val_l); }
    break;

  case 16:
/* Line 1787 of yacc.c  */
#line 86 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(1) - (1)].val_d);        if(data.Debug>0){printf("this is a double\n");}}
    break;

  case 17:
/* Line 1787 of yacc.c  */
#line 87 "calc_bison.y"
    { (yyval.val_d) = data.variable[variable_ID((yyvsp[(1) - (1)].string))].value;   }
    break;

  case 18:
/* Line 1787 of yacc.c  */
#line 88 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(3) - (3)].val_l); create_variable_ID((yyvsp[(1) - (3)].string), (yyvsp[(3) - (3)].val_l));  }
    break;

  case 19:
/* Line 1787 of yacc.c  */
#line 89 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(3) - (3)].val_d); create_variable_ID((yyvsp[(1) - (3)].string), (yyvsp[(3) - (3)].val_d));  }
    break;

  case 20:
/* Line 1787 of yacc.c  */
#line 90 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(3) - (3)].val_l); create_variable_ID((yyvsp[(1) - (3)].string), (yyvsp[(3) - (3)].val_l));  if(data.Debug>0){printf("creating long\n");}}
    break;

  case 21:
/* Line 1787 of yacc.c  */
#line 91 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(3) - (3)].val_d); create_variable_ID((yyvsp[(1) - (3)].string), (yyvsp[(3) - (3)].val_d));  if(data.Debug>0){printf("creating double\n");}}
    break;

  case 22:
/* Line 1787 of yacc.c  */
#line 92 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(1) - (3)].val_l) + (yyvsp[(3) - (3)].val_d);   if(data.Debug>0){printf("long + double\n");}}
    break;

  case 23:
/* Line 1787 of yacc.c  */
#line 93 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(1) - (3)].val_d) + (yyvsp[(3) - (3)].val_l);   if(data.Debug>0){printf("double + long\n");}}
    break;

  case 24:
/* Line 1787 of yacc.c  */
#line 94 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(1) - (3)].val_d) + (yyvsp[(3) - (3)].val_d);   if(data.Debug>0){printf("double + double\n");}}
    break;

  case 25:
/* Line 1787 of yacc.c  */
#line 95 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(1) - (3)].val_l) - (yyvsp[(3) - (3)].val_d);   if(data.Debug>0){printf("long - double\n");}}
    break;

  case 26:
/* Line 1787 of yacc.c  */
#line 96 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(1) - (3)].val_d) - (yyvsp[(3) - (3)].val_l);   if(data.Debug>0){printf("double - long\n");}}
    break;

  case 27:
/* Line 1787 of yacc.c  */
#line 97 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(1) - (3)].val_d) - (yyvsp[(3) - (3)].val_d);   if(data.Debug>0){printf("double - double\n");}}
    break;

  case 28:
/* Line 1787 of yacc.c  */
#line 98 "calc_bison.y"
    { (yyval.val_d) = (double) (yyvsp[(1) - (3)].val_l) * (yyvsp[(3) - (3)].val_d);   if(data.Debug>0){printf("long * double\n");}}
    break;

  case 29:
/* Line 1787 of yacc.c  */
#line 99 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(1) - (3)].val_d) * (yyvsp[(3) - (3)].val_l);   if(data.Debug>0){printf("double * long\n");}}
    break;

  case 30:
/* Line 1787 of yacc.c  */
#line 100 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(1) - (3)].val_d) * (yyvsp[(3) - (3)].val_d);   if(data.Debug>0){printf("double * double\n");}}
    break;

  case 31:
/* Line 1787 of yacc.c  */
#line 101 "calc_bison.y"
    { (yyval.val_d) = (double) (yyvsp[(1) - (3)].val_l) / (yyvsp[(3) - (3)].val_l);   if(data.Debug>0){printf("long / long\n");}}
    break;

  case 32:
/* Line 1787 of yacc.c  */
#line 102 "calc_bison.y"
    { (yyval.val_d) = (double) (yyvsp[(1) - (3)].val_l) / (yyvsp[(3) - (3)].val_d);   if(data.Debug>0){printf("long / double\n");}}
    break;

  case 33:
/* Line 1787 of yacc.c  */
#line 103 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(1) - (3)].val_d) / (yyvsp[(3) - (3)].val_l);   if(data.Debug>0){printf("double / long\n");}}
    break;

  case 34:
/* Line 1787 of yacc.c  */
#line 104 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(1) - (3)].val_d) / (yyvsp[(3) - (3)].val_d);   if(data.Debug>0){printf("double / double\n");}}
    break;

  case 35:
/* Line 1787 of yacc.c  */
#line 105 "calc_bison.y"
    { (yyval.val_d) = -(yyvsp[(2) - (2)].val_d);       if(data.Debug>0){printf("-double\n");}}
    break;

  case 36:
/* Line 1787 of yacc.c  */
#line 106 "calc_bison.y"
    { (yyval.val_d) = pow ((double) (yyvsp[(1) - (3)].val_l),(yyvsp[(3) - (3)].val_d));  if(data.Debug>0){printf("long ^ double\n");}}
    break;

  case 37:
/* Line 1787 of yacc.c  */
#line 107 "calc_bison.y"
    { (yyval.val_d) = pow ((yyvsp[(1) - (3)].val_d),(double) (yyvsp[(3) - (3)].val_l));  if(data.Debug>0){printf("double ^ long\n");}}
    break;

  case 38:
/* Line 1787 of yacc.c  */
#line 108 "calc_bison.y"
    { (yyval.val_d) = pow ((yyvsp[(1) - (3)].val_d),(yyvsp[(3) - (3)].val_d));  if(data.Debug>0){printf("double ^ double\n");}}
    break;

  case 39:
/* Line 1787 of yacc.c  */
#line 109 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(1) - (3)].fnctptr)((yyvsp[(2) - (3)].val_d));  if(data.Debug>0){printf("double=func(double)\n");}}
    break;

  case 40:
/* Line 1787 of yacc.c  */
#line 110 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(1) - (3)].fnctptr)((double) (yyvsp[(2) - (3)].val_l));  if(data.Debug>0){printf("double=func(double)\n");}}
    break;

  case 41:
/* Line 1787 of yacc.c  */
#line 111 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(1) - (5)].fnctptr)((yyvsp[(2) - (5)].val_d),(yyvsp[(4) - (5)].val_d));  if(data.Debug>0){printf("double=func(double,double)\n");}}
    break;

  case 42:
/* Line 1787 of yacc.c  */
#line 112 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(1) - (5)].fnctptr)((double) (yyvsp[(2) - (5)].val_l),(yyvsp[(4) - (5)].val_d));  if(data.Debug>0){printf("double=func(long->double,double)\n");}}
    break;

  case 43:
/* Line 1787 of yacc.c  */
#line 113 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(1) - (5)].fnctptr)((yyvsp[(2) - (5)].val_d),(double) (yyvsp[(4) - (5)].val_l));  if(data.Debug>0){printf("double=func(double,long->double)\n");}}
    break;

  case 44:
/* Line 1787 of yacc.c  */
#line 114 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(1) - (5)].fnctptr)((double) (yyvsp[(2) - (5)].val_l),(double) (yyvsp[(4) - (5)].val_l));  if(data.Debug>0){printf("double=func(long->double,long->double)\n");}}
    break;

  case 45:
/* Line 1787 of yacc.c  */
#line 115 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(1) - (3)].fnctptr)((yyvsp[(2) - (3)].string));  if(data.Debug>0){printf("double=func(image)\n");}}
    break;

  case 46:
/* Line 1787 of yacc.c  */
#line 116 "calc_bison.y"
    { (yyval.val_d) = (yyvsp[(2) - (3)].val_d);                         }
    break;

  case 47:
/* Line 1787 of yacc.c  */
#line 120 "calc_bison.y"
    {(yyval.string) = strdup((yyvsp[(1) - (1)].string));        data.cmdargtoken[data.cmdNBarg].type = 3; if(data.Debug>0){printf("this is a string (new variable/image)\n");}}
    break;

  case 48:
/* Line 1787 of yacc.c  */
#line 121 "calc_bison.y"
    {(yyval.string) = strdup((yyvsp[(1) - (1)].string));        data.cmdargtoken[data.cmdNBarg].type = 4; if(data.Debug>0){printf("this is a string (existing image)\n");}}
    break;

  case 49:
/* Line 1787 of yacc.c  */
#line 122 "calc_bison.y"
    {(yyval.string) = strdup((yyvsp[(1) - (1)].string));        data.cmdargtoken[data.cmdNBarg].type = 5; if(data.Debug>0){printf("this is a string (command)\n");}}
    break;

  case 50:
/* Line 1787 of yacc.c  */
#line 123 "calc_bison.y"
    {(yyval.string) = strdup((yyvsp[(1) - (3)].string));        delete_image_ID((yyvsp[(1) - (3)].string)); chname_image_ID((yyvsp[(3) - (3)].string),(yyvsp[(1) - (3)].string)); if(data.Debug>0){printf("changing name\n");}}
    break;

  case 51:
/* Line 1787 of yacc.c  */
#line 124 "calc_bison.y"
    {(yyval.string) = strdup((yyvsp[(1) - (3)].string));        chname_image_ID((yyvsp[(3) - (3)].string),(yyvsp[(1) - (3)].string)); if(data.Debug>0){printf("changing name\n");}}
    break;

  case 52:
/* Line 1787 of yacc.c  */
#line 125 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_add((yyvsp[(1) - (3)].string), (yyvsp[(3) - (3)].string), calctmpimname); (yyval.string) = strdup(calctmpimname); if(data.Debug>0){printf("image + image\n");}}
    break;

  case 53:
/* Line 1787 of yacc.c  */
#line 126 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_cstadd((yyvsp[(1) - (3)].string),(double) (yyvsp[(3) - (3)].val_d), calctmpimname); (yyval.string)=strdup(calctmpimname);  if(data.Debug>0){printf("image + double\n");}}
    break;

  case 54:
/* Line 1787 of yacc.c  */
#line 127 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_cstadd((yyvsp[(1) - (3)].string),(double) (yyvsp[(3) - (3)].val_l), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("image + long\n");}}
    break;

  case 55:
/* Line 1787 of yacc.c  */
#line 128 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_cstadd((yyvsp[(3) - (3)].string),(double) (yyvsp[(1) - (3)].val_d), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("double + image\n");}}
    break;

  case 56:
/* Line 1787 of yacc.c  */
#line 129 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_cstadd((yyvsp[(3) - (3)].string),(double) (yyvsp[(1) - (3)].val_l), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("long + image\n");}}
    break;

  case 57:
/* Line 1787 of yacc.c  */
#line 130 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_sub((yyvsp[(1) - (3)].string), (yyvsp[(3) - (3)].string), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("image + image\n");}}
    break;

  case 58:
/* Line 1787 of yacc.c  */
#line 131 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_cstadd((yyvsp[(1) - (3)].string),(double) -(yyvsp[(3) - (3)].val_d), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("image - double\n");}}
    break;

  case 59:
/* Line 1787 of yacc.c  */
#line 132 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_cstadd((yyvsp[(1) - (3)].string),(double) -(yyvsp[(3) - (3)].val_l), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("image - long\n");}}
    break;

  case 60:
/* Line 1787 of yacc.c  */
#line 133 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_cstadd((yyvsp[(3) - (3)].string),(double) -(yyvsp[(1) - (3)].val_d), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("double - image\n");}}
    break;

  case 61:
/* Line 1787 of yacc.c  */
#line 134 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_cstadd((yyvsp[(3) - (3)].string),(double) -(yyvsp[(1) - (3)].val_l), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("long - image\n");}}
    break;

  case 62:
/* Line 1787 of yacc.c  */
#line 135 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_mult((yyvsp[(1) - (3)].string), (yyvsp[(3) - (3)].string), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("image * image\n");}}
    break;

  case 63:
/* Line 1787 of yacc.c  */
#line 136 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_cstmult((yyvsp[(1) - (3)].string),(double) (yyvsp[(3) - (3)].val_d), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("image * double\n");}}
    break;

  case 64:
/* Line 1787 of yacc.c  */
#line 137 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_cstmult((yyvsp[(1) - (3)].string),(double) (yyvsp[(3) - (3)].val_l), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("image * long\n");}}
    break;

  case 65:
/* Line 1787 of yacc.c  */
#line 138 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_cstmult((yyvsp[(3) - (3)].string),(double) (yyvsp[(1) - (3)].val_d), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("double * image\n");}}
    break;

  case 66:
/* Line 1787 of yacc.c  */
#line 139 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_cstmult((yyvsp[(3) - (3)].string),(double) (yyvsp[(1) - (3)].val_l), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("long * image\n");}}
    break;

  case 67:
/* Line 1787 of yacc.c  */
#line 140 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_div((yyvsp[(1) - (3)].string), (yyvsp[(3) - (3)].string), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("image / image\n");}}
    break;

  case 68:
/* Line 1787 of yacc.c  */
#line 141 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_cstdiv((yyvsp[(1) - (3)].string),(double) (yyvsp[(3) - (3)].val_d), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("image - double\n");}}
    break;

  case 69:
/* Line 1787 of yacc.c  */
#line 142 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_cstdiv((yyvsp[(1) - (3)].string),(double) (yyvsp[(3) - (3)].val_l), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("image - long\n");}}
    break;

  case 70:
/* Line 1787 of yacc.c  */
#line 143 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_cstdiv((yyvsp[(3) - (3)].string),(double) (yyvsp[(1) - (3)].val_d), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("double - image\n");}}
    break;

  case 71:
/* Line 1787 of yacc.c  */
#line 144 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_cstdiv((yyvsp[(3) - (3)].string),(double) (yyvsp[(1) - (3)].val_l), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("long - image\n");}}
    break;

  case 72:
/* Line 1787 of yacc.c  */
#line 145 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_cstpow((yyvsp[(1) - (3)].string),(double) (yyvsp[(3) - (3)].val_l), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("image^long\n");}}
    break;

  case 73:
/* Line 1787 of yacc.c  */
#line 146 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_cstpow((yyvsp[(1) - (3)].string),(double) (yyvsp[(3) - (3)].val_d), calctmpimname); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("image^double\n");}}
    break;

  case 74:
/* Line 1787 of yacc.c  */
#line 147 "calc_bison.y"
    {sprintf(calctmpimname,"_tmpcalc%ld",data.calctmp_imindex); data.calctmp_imindex++; arith_image_function_d_d((yyvsp[(2) - (3)].string), calctmpimname, (yyvsp[(1) - (3)].fnctptr)); (yyval.string)=strdup(calctmpimname); if(data.Debug>0){printf("double_func(double)\n");}}
    break;

  case 75:
/* Line 1787 of yacc.c  */
#line 148 "calc_bison.y"
    { (yyval.string) = strdup((yyvsp[(2) - (3)].string));                         }
    break;


/* Line 1787 of yacc.c  */
#line 1984 "calc_bison.c"
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
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
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
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

  /* Else will try to reuse lookahead token after shifting the error
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
      if (!yypact_value_is_default (yyn))
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

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


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

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
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


/* Line 2050 of yacc.c  */
#line 153 "calc_bison.y"



#include <stdio.h>



yyerror (s)  /* Called by yyparse on error */
     char *s;
{
  printf ("PARSING ERROR ON COMMAND LINE ARG %ld: %s\n", data.cmdNBarg, s);
  data.parseerror = 1;
}
