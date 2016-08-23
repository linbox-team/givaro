# Check for GMP
# Copyright(c)'1994-2009,2003,2013 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.
#
# Modified by Pascal Giorgi, 2003-12-03
# Modified by BB, 2013-5-22
# Hacked by Charles Bouillaguet, 2016-08-23

dnl Test for the GNU Multiprecision library and define GMP_CFLAGS and GMP_LIBS

AC_DEFUN([GIV_CHECK_GMP], [

	min_gmp_release=$1

	########## ./configure parameter
	AC_ARG_WITH(gmp, 
		[AC_HELP_STRING(
			[--with-gmp=<path>], 
			[Path to the GMP library. If unspecified, that means the library is reachable 
			 with the standard search path (/usr or /usr/local).]
		) ],
		[ if test "$withval" != ""; then
			GMP_HOME_PATH="$withval"
		  fi ],
		[GMP_HOME_PATH=${DEFAULT_CHECKING_PATH}]
	)

	######### Check for existence
	BACKUP_CFLAGS=${CFLAGS}
	BACKUP_CXXFLAGS=${CXXFLAGS}
	BACKUP_LIBS=${LIBS}

	AC_MSG_CHECKING(for GMP)

	for GMP_HOME in ${GMP_HOME_PATH}
	do
		if test -r "$GMP_HOME/include/gmp.h"; then

			if test "x$GMP_HOME" != "x/usr"; then
				GMP_CFLAGS="-I${GMP_HOME}/include"
				GMP_LIBS="-L${GMP_HOME}/lib"
			fi
			GMP_LIBS="$GMP_LIBS -lgmp"
			

			######### try to compile
			CXXFLAGS="${BACKUP_CFLAGS} ${GMP_CFLAGS}"
			LIBS="${BACKUP_LIBS} ${GMP_LIBS}"
			AC_LINK_IFELSE(
				[ 
					AC_LANG_PROGRAM(
						[[#include <gmp.h>]], 
						[[ mpz_t a; mpz_init (a); ]] )
				],
				[ 
					gmp_found="yes"
				 	break # exit the loop
				],
				[ gmp_found="no (gmp.h found but linking failed)" ]
			)
		else
	    		gmp_found="no (gmp.h not found)"
		fi
	done

	AC_MSG_RESULT(${gmp_found})

	if test "x$gmp_found" != "xyes" ; then
		echo '-------------------------------'
		AC_MSG_ERROR(ERROR: GMP not found/usable!)
	fi

	##### OK, we have found a working GMP. Check if it has c++ bindings, and is recent enough
	
	AC_MSG_CHECKING(for GMP version and cxx support)
	GMP_LIBS="$GMP_LIBS -lgmpxx"
	AC_LINK_IFELSE( 
		[ AC_LANG_PROGRAM(
			[[ #include <gmpxx.h> ]],
			[[ mpz_class a(2), b(3), c(5); ]]
		) ],
		,
		[ AC_MSG_ERROR(your GMP does not have c++ support. Compile GMP with --enable-cxx)]
	)

	AC_TRY_RUN(
		[ 
			#include <gmp.h>
			int main () {
				return (__GNU_MP_RELEASE < $min_gmp_release);
			}
		],
		,
		[AC_MSG_ERROR(your GMP is too old. GMP release >= $min_gmp_release needed)]
	)
	
	AC_MSG_RESULT(OK)
	AC_SUBST(GMP_CFLAGS)
	AC_SUBST(GMP_LIBS)
	AC_DEFINE(HAVE_GMP, 1 ,[Define if GMP is installed and OK])

	CFLAGS=${BACKUP_CFLAGS}
	CXXFLAGS=${BACKUP_CXXFLAGS}
	LIBS=${BACKUP_LIBS}
])
