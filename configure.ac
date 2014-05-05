#
# Process this file with autoconf to produce a configure script.
#
AC_INIT([Cfits], [3.2.00], [oliv.guyon@gmail.com])
AC_PREREQ(2.61)
#AC_CONFIG_HEADER([config.h])

# sanity check a Cfits file to see if we're really in the build directory
AC_CONFIG_SRCDIR([src/Cfits.c])

AC_CONFIG_MACRO_DIR([m4])
# Init Automake
AM_INIT_AUTOMAKE([foreign -Wall])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
# Init Libtool
LT_INIT

# find the c-complier
#AC_PROG_CC
AM_PROG_CC_C_O

# check for openMP
AC_OPENMP


# Support for Bison and Flex
AC_PROG_YACC
AM_PROG_LEX

# Support for pthreads
ACX_PTHREAD

#LIBS="$PTHREAD_LIBS $LIBS"
#CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
#CC="$PTHREAD_CC"

#
# AC_PREFIX_DEFAULT(~/Cfits) # howto coerce default "prefix" from /usr/local/ to ~/Cfits
#
# Coerce @prefix@ default from /usr/local/ to fullpath of dir configure is run in
AC_PREFIX_DEFAULT(`readlink -f .`)


CPPFLAGS="$CPPFLAGS -I/usr/local/include/"

#LDFLAGS="$LDFLAGS -L/usr/lib/"

AC_ARG_ENABLE(fftwmt,
 [AS_HELP_STRING(  --enable-fftwmt         Enable FFTW multi-thread )],
 [ fftwmt=yes
 ])




#
# Look for libraries 
#


AC_CHECK_LIB([m],[main],
		[],
		[ echo "Cannot find required library:math";
		  cat doc/configure.howto;
		  AC_MSG_FAILURE(["FATAL...Cannot find required library:math"]);
		  exit 0;
		])

AC_CHECK_LIB([rt], [clock],
		[],
		[ echo "Cannot find required library:rt";
		  cat doc/configure.howto;
		  AC_MSG_FAILURE(["FATAL...Cannot find required library:rt"]);
		  exit 0;
		])


AC_CHECK_LIB([gslcblas],[main],
		[],
		[ echo "Cannot find required library:gslcblas";
		  cat doc/configure.howto;
		  AC_MSG_FAILURE(["FATAL...Cannot find required library:gslcblas"]);
		  exit 0;
		])


AC_CHECK_LIB([gsl],[main],
		[],
		[ echo "Cannot find required library:gsl";
		  cat doc/configure.howto;
		  AC_MSG_FAILURE(["FATAL...Cannot find required library:gsl"]);
		  exit 0;
		])


AC_CHECK_LIB([cfitsio],[main],
		[],
		[ echo "Cannot find required static library:fitsio";
		  echo "did you install libcfitsio3-dev package ?";
		  cat doc/configure.howto;
		  AC_MSG_FAILURE(["FATAL...Cannot find required library:fitsio"]);
		  exit 0;
		])
		
		
AC_SEARCH_LIBS( [readline], [readline], 
                [AC_DEFINE([HAVE_READLINE_H], 1, "Define to 1 if you have <readline.h>") ], 
                [ echo "Cannot find required static library:readline";
		  echo "did you install libreadline<#>-dev package ?";
                  cat doc/configure.howto;
                  AC_MSG_FAILURE(["FATAL...Cannot find required library:readlline"]); 
                  exit 0;
                ])

AC_CHECK_LIB( [ncurses], [initscr], 
		[],
		[ echo "Cannot find required library:ncurses";
 		  echo "did you install ncurses package ?";
                  cat doc/configure.howto;
                  AC_MSG_FAILURE(["FATAL...Cannot find required library:ncurses"]); 
                ]
)

AC_CHECK_LIB( [fftw3f], [main], 
                [], 
                [ echo "Cannot find required library:fftw3f";
 		  echo "did you install fftw3f-dev package ?";
                  cat doc/configure.howto;
                  AC_MSG_FAILURE(["FATAL...Cannot find required library:fftw3f"]); 
                ]
)

AC_CHECK_LIB( [fftw3], [main], 
                [], 
                [ echo "Cannot find required library:fftw3";
		  echo "did you install fftw3-dev package ?";
                  cat doc/configure.howto;
                  AC_MSG_FAILURE(["FATAL...Cannot find required library:fftw3"]); 
                ]
)

AC_CHECK_LIB([gomp], [omp_get_num_threads])

#
# fftw multi-thread
#	
if   test "x$fftwmt" = "xyes" ; then
	AC_CHECK_LIB( [fftw3f_threads], [main],
        	        [], 
        	        [ echo "Cannot find required library:fftw3f_threads";
        	          cat doc/configure.howto;
        	          AC_MSG_FAILURE(["FATAL...Cannot find required library:fftw3f_threads"]); 
        	        ]
	)	
	AC_CHECK_LIB( [fftw3_threads], [main],
        	        [], 
        	        [ echo "Cannot find required library:fftw3_threads";
        	          cat doc/configure.howto;
        	          AC_MSG_FAILURE(["FATAL...Cannot find required library:fftw3f_threads"]); 
        	        ]
	)	
	
	AC_CHECK_LIB([pthread], [pthread_create],
		[],
		[ echo "Cannot find required library:pthread";
		  cat doc/configure.howto;
		  AC_MSG_FAILURE(["FATAL...Cannot find required library:pthread"]);
		  exit 0;
		])
	
	AC_DEFINE(FFTWMT)
fi


#AC_CHECK_LIB([usb], [usb_close], ,[AC_MSG_ERROR(
# [libusb library not found.  Check with vendor or see http://libusb.sourceforge.net])],
# [$USBLDFLAGS])


#AC_SEARCH_LIBS( [bob],      
#                [tom], 
#                [AC_DEFINE([HAVE_BOB_H], 1,
#                           "Define to 1 if you have <readline.h>")] , 
#                [AC_MSG_ERROR([$str]); exit 0] )

#AC_CHECK_HEADERS( [bob.h], 
#                  [], 
#                  [ echo "Cannot find required header:bob.h";
#                    cat doc/configure.howto;
#                    AC_MSG_FAILURE(["Cannot find required header:bob.h"]);
#                  ]
#)  

AC_CHECK_HEADERS( [fftw3.h], 
                  [], 
                  [ echo "Cannot find required header:fftw3.h";
                    cat doc/configure.howto;
                    AC_MSG_FAILURE(["Cannot find required header:fftw3.h"]);
                  ]
)  

AC_CHECK_HEADERS( [fitsio.h], 
                  [], 
                  [ echo "Cannot find required header:fitsio.h";
                    cat doc/configure.howto;
                    AC_MSG_FAILURE(["Cannot find required header:fitsio.h"]);
                  ]
)  




#AC_MSG_NOTICE([Looking for special Cfits libraries])
#AC_INCLUDES_DEFAULT
#AC_CHECK_LIB(gsl,gsl_ran_gaussian_ziggurat, , echo "FATAL..."; exit 0; )
#AC_PROG_LEX
#AC_PROG_YACC
#AC_PROG_INSTALL

#AC_MSG_CHECKING (feature-description)
#AC_MSG_NOTICE([eg:checking if stack overflow is detectable])
#AC_MSG_FAILURE (error-description, [exit-status])

# Create our own a symbol for execname until we can find the right one
EXECNAME=Cfits
AC_SUBST(EXECNAME)


AC_CONFIG_FILES(Makefile src/Makefile src/00CORE/Makefile src/COREMOD_arith/Makefile src/COREMOD_iofits/Makefile src/COREMOD_memory/Makefile src/COREMOD_tools/Makefile src/info/Makefile src/fft/Makefile src/statistic/Makefile src/image_gen/Makefile src/WFpropagate/Makefile src/ZernikePolyn/Makefile src/coronagraphs/Makefile src/OpticsMaterials/Makefile src/PIAACMCsimul/Makefile src/image_format/Makefile src/SCExAO_DM/Makefile src/AOsystSim/Makefile src/AOloopControl/Makefile)

AC_OUTPUT

