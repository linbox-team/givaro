# Copyright(c)'2011 by The Givaro group
# Written by BB <bboyer@imag.fr>
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

dnl enable basic debug mode.
AC_DEFUN([AC_DEBUG],
[AC_MSG_CHECKING([whether to enable debugging options in the library])
  AC_ARG_ENABLE(debug,
[AC_HELP_STRING([--enable-debug=yes|no], [enable debugging options in library])],
      USE_DEBUG=$enableval,
      USE_DEBUG=no)
  AC_MSG_RESULT([$USE_DEBUG])
  AM_CONDITIONAL(DEBUG, [test x$USE_DEBUG = xyes])
  DBG=$USE_DEBUG
  AC_SUBST(DBG)dnl
]
)

AC_DEFUN([AC_PROFILE],
[AC_MSG_CHECKING([whether to enable profiling everything in the library])
  AC_ARG_ENABLE(profile,
[AC_HELP_STRING([--enable-profile=yes|no], [enable profiling options in library])],
      USE_PROFILE=$enableval,
      USE_PROFILE=no)
  AC_MSG_RESULT([$USE_PROFILE])
  AM_CONDITIONAL(PROFILE, [test $USE_PROFILE = yes])
  PROF=$USE_PROFILE
  AC_SUBST(PROF)dnl
]
)

dnl Enable warnings from compiler.
AC_DEFUN([AC_WARNINGS],
[AC_MSG_CHECKING([whether to enable warnings when compiling the library])
  AC_ARG_ENABLE(warnings,
[AC_HELP_STRING([--enable-warnings=yes|full|no],  [enable warnings when compiling the library.
If nothing or yes is given, more aggressive compiler warnings are passed to the compiler.
If full is given, we become paranoïd about warnings and treat them as errors.])],
      USE_WARNINGS=$enableval,
      USE_WARNINGS=no)
  AC_MSG_RESULT([$USE_WARNINGS])
  dnl  AM_CONDITIONAL(WARNINGS, [test $USE_WARNINGS = yes])
  WARN=$USE_WARNINGS
  AC_SUBST(WARN)dnl
]dnl
)dnl

CCNAM=""

AC_DEFUN([AC_COMPILER_NAME], [
		AC_MSG_CHECKING(for family name of compiler)

		dnl CHECKING for various compilers
		dnl ICC ?
		AC_TRY_RUN( [
           #ifdef __INTEL_COMPILER
   int main() { return 0 ; }
   #else
   pas intel
		   #endif],
		[ AC_MSG_RESULT(icc)
   CCNAM=icc
   AC_SUBST(CCNAM)
		])

dnl PATHSCALE > 4 ?
		AS_IF([ test -z "${CCNAM}"], [
			AC_TRY_RUN( [
				#ifdef __PATHSCALE__
				   int main() { return !(__PATHCC__ >= 4) ; }
			   #else
				   pas ekopath non plus.
				#endif], [
		AC_MSG_RESULT(eko)
		CCNAM=eko
		AC_SUBST(CCNAM) ])
		])

dnl CLANG > 3.1 ?
		AS_IF([ test -z "${CCNAM}"], [
			AC_TRY_RUN( [
				#ifdef __clang__
				   int main() { return !(__clang_major__  >=3 && __clang_minor__ >=1) ; }
			   #else
				   pas clang non plus.
				#endif], [
		AC_MSG_RESULT(clang31)
		CCNAM=clang31
		AC_SUBST(CCNAM) ])
		])

dnl CLANG > 3 ?
		AS_IF([ test -z "${CCNAM}"], [
			AC_TRY_RUN( [
				#ifdef __clang__
				   int main() { return !(__clang_major__  >=3) ; }
			   #else
				   pas clang non plus.
				#endif], [
		AC_MSG_RESULT(clang31)
		CCNAM=clang
		AC_SUBST(CCNAM) ])
		])


dnl GCC >= 4.8 ?
		AS_IF([ test -z "${CCNAM}"], [
			AC_TRY_RUN( [
				#ifdef __GNUC__
				   int main() { return !(__GNUC__ >= 5 || (__GNUC__ == 4  && __GNUC_MINOR__ > 7 )) ; }
				#else
				   pas gcc non plus ???
				#endif], [
		CCNOM=gcc
		AS_IF([ test -n "${CC}" ], [CCNOM="`$CC --version 2>&1|  awk 'NR<2{print $1}'`"])
		CCNAM=gcc48
		AC_SUBST(CCNAM)
		AC_MSG_RESULT($CCNOM)
		])
		])

dnl GCC > 4.2 ?
		AS_IF([ test -z "${CCNAM}"], [
			AC_TRY_RUN( [
				#ifdef __GNUC__
				   int main() { return !(__GNUC__ >= 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 2)) ; }
   #else
				   pas gcc non plus ???
				#endif], [
		CCNOM=gcc
		AS_IF([ test -n "${CC}" ], [CCNOM="`$CC --version 2>&1|  awk 'NR<2{print $1}'`"])
   CCNAM=gcc
   AC_SUBST(CCNAM)
		AC_MSG_RESULT($CCNOM)
   ])
		])

		dnl  autre ?

		AS_IF([ test -z "${CCNAM}"],
				[ AC_MSG_RESULT(unknown)
   CCNAM=unknown
   AC_SUBST(CCNAM)
				echo
				echo " *** unknow compiler. please file a bug "
				echo
				])
])

dnl compile library or make it mostly inlined headers ?

AC_DEFUN([AC_INLINE],
[AC_MSG_CHECKING([whether to inline or not most of the code ?])
  AC_ARG_ENABLE(inline,
[AC_HELP_STRING([--enable-inline],  [enable inlining most of the code])],
      USE_INLINE=$enableval
      AC_DEFINE(INLINE_ALL,1,[Define if you want most code inlined]) ,
      USE_INLINE=no)
  AC_MSG_RESULT([$USE_INLINE])
  AM_CONDITIONAL(GIVARO_INLINE_ALL, [test $USE_INLINE = yes])
  AC_SUBST(GIVARO_INLINE_ALL)
  echo $GIVARO_INLINE_ALL
  dnl  DBG=$USE_DEBUG
  dnl  AC_SUBST(DBG)dnl
]
)



