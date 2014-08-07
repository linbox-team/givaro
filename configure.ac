# Copyright(c)'1994-2011 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.
#/

AC_PREREQ([2.61])

AC_INIT([Givaro],[3.8.1],[Jean-Guillaume.Dumas@imag.fr],[givaro],
		[http://ljk.imag.fr/CASYS/LOGICIELS/givaro])

AC_CONFIG_MACRO_DIR([macros])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([1.8 gnu no-dependencies -Wall -Wno-portability])
AC_CONFIG_HEADERS([config.h])
AX_PREFIX_CONFIG_H(givaro-config.h, __GIVARO)
AC_PATH_PROG(RM, rm, $FALSE)
RM="$RM -f"


AM_MAINTAINER_MODE
AM_DEP_TRACK
AM_OUTPUT_DEPENDENCY_COMMANDS

AM_ACLOCAL_INCLUDE(macros)

# work around to fix the backward compatibility issue of automake 1.10 with 1.9 (pb with MKDIR_P)
AC_SUBST([MKDIR_P])

AC_LANG([C++])

echo "-----------------------------------------------"

AC_DEBUG
AC_PROFILE
AC_WARNINGS
AC_INLINE # DON'T TOUCH

echo "-----------------------------------------------"

# CFLAGS=${CFLAGS:-$DEFAULT_CFLAGS}
# CXXFLAGS=${CXXFLAGS:-$DEFAULT_CXXFLAGS}

######################################################
# Try and pass different flags according to compiler #
######################################################


# disable default -g -O2 CXXFLAGS
: ${CXXFLAGS=""}

#set CXX
AC_PROG_CXX

AC_COMPILER_NAME

if test "x$CCNAM" != "xgcc44" ; then
AX_CXX_COMPILE_STDCXX_11([ext],[optional])
fi

AC_SUBST([DEFAULT_CFLAGS])
AC_SUBST([DEBUG_CFLAGS])
AC_SUBST([TESTS_CFLAGS])

TESTS_CFLAGS="-O0"
DEBUG_CFLAGS="-g"
DEFAULT_CFLAGS=""
WARN_CFLAGS="-Wall"

#TODO use -fast for icc, -ipa for eko...
if test "x$DBG" = "xyes" ; then
	DEFAULT_CFLAGS="-O0 ${DEFAULT_CFLAGS} " #those are CXXFLAGS
	DEBUG_CFLAGS="${DEBUG_CFLAGS} -DDEBUG -DGIVARO_DEBUG"
else
	DEFAULT_CFLAGS="-O2 ${DEFAULT_CFLAGS} "
	DEBUG_CFLAGS="${DEBUG_CFLAGS} -DNDEBUG -UGIVARO_DEBUG -UDEBUG"
fi

if test "x$PROF" = "xyes" ; then
	DEFAULT_CFLAGS="${DEFAULT_CFLAGS} -pg"
fi

if test "x$WARN" = "xyes" -o "x$WARN" = "xfull" ; then
	case x${CCNAM} in
		xicc)
			WARN_CFLAGS="${WARN_CFLAGS} -Wcheck"
		WARN_CFLAGS="${WARN_CFLAGS} -Wall -Wno-unused-parameter -Wuninitialized -Wconversion -Wcast-qual -pedantic -Wshadow -Wpointer-arith  -Wwrite-strings -Wno-long-long"
		WARN_CFLAGS="${WARN_CFLAGS} -Wextra -ansi"
			;;
		xeko)
			WARN_CFLAGS="${WARN_CFLAGS} -Wno-unused-parameter"
			;;
		xgcc|xgcc44)
			WARN_CFLAGS="${WARN_CFLAGS} -Wextra -Wno-unused-parameter"
			if test "x${WARN}" = "xfull" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -Wuninitialized -Wconversion -Wcast-qual -pedantic -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wno-long-long  -Wno-variadic-macros -Wno-vla"
			fi
			if test "x${HAVE_CXX11}" = "x0" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -ansi"
			fi
		;;
		xgcc48)
			WARN_CFLAGS="${WARN_CFLAGS} -Wextra -Wno-unused-parameter"
			if test "x${WARN}" = "xfull" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -Wuninitialized -Wconversion -Wcast-qual -pedantic -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wno-long-long  -Wno-variadic-macros -Wno-vla"
				# WARN_CFLAGS="${WARN_CFLAGS} -fsanitize=address"
			fi
			if test "x${HAVE_CXX11}" = "x0" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -ansi"
			fi

			;;

		xclang)
			WARN_CFLAGS="${WARN_CFLAGS} -Wextra -Wno-unused-parameter"
			if test "x${WARN}" = "xfull" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -Wuninitialized -Wconversion -Wcast-qual -pedantic -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wno-long-long -Wno-vla-extension -Wno-variadic-macros"
				WARN_CFLAGS="${WARN_CFLAGS} -D__STRICT_ANSI__"
			fi
			;;
		xclang31)
			WARN_CFLAGS="${WARN_CFLAGS} -Wextra -Wno-unused-parameter"
			if test "x${WARN}" = "xfull" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -Wuninitialized -Wconversion -Wcast-qual -pedantic -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wno-long-long  -g -Wno-vla-extension -Wno-variadic-macros"
				WARN_CFLAGS="${WARN_CFLAGS} -D__STRICT_ANSI__"
				# WARN_CFLAGS="${WARN_CFLAGS} -fsanitize=address"
			fi
			;;

		*)
			echo
			echo "*******************************************************"
			echo "unsupported compiler ($CCNAM). Please file a bug."
			echo "*******************************************************"
			echo
			WARN_CFLAGS="${WARN_CFLAGS}"
	esac
fi


DEFAULT_CFLAGS="${DEFAULT_CFLAGS} ${WARN_CFLAGS} ${DEBUG_CFLAGS}"
TESTS_CFLAGS="${TESTS_CFLAGS} ${WARN_CFLAGS} ${DEBUG_CFLAGS}"


AC_HEADER_STDC
AC_PROG_LIBTOOL

# newer libtool...
LT_PREREQ([2.2])
LT_INIT


echo "-----------------------------------------------"
echo "        START  GIVARO CONFIG                   "
echo "-----------------------------------------------"

echo "-----------------------------------------------"

echo "-----------------------------------------------"

# Machine characteristics

AC_CHECK_SIZEOF(char, 8)
AC_CHECK_SIZEOF(short, 16)
AC_CHECK_SIZEOF(int, 32)
AC_CHECK_SIZEOF(long, 32)
AC_CHECK_SIZEOF(long long, 64)
AC_CHECK_SIZEOF(__int64, 64)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stddef.h stdlib.h string.h sys/time.h stdint.h pthread.h])


# check endianness of the architecture
AC_C_BIGENDIAN(
  [AC_DEFINE(HAVE_BIG_ENDIAN, 1, [Define that architecture uses big endian storage])],
  [AC_DEFINE(HAVE_LITTLE_ENDIAN, 1, [Define that architecture uses little endian storage])],
  [])

echo "Creating symbolic link for compilation"
if test ! -L src/kernel/system/givaro  ; then ln -s . src/kernel/system/givaro  ; fi
if test ! -L src/kernel/bstruct/givaro ; then ln -s . src/kernel/bstruct/givaro ; fi
if test ! -L src/kernel/integer/givaro ; then ln -s . src/kernel/integer/givaro ; fi
if test ! -L src/kernel/memory/givaro  ; then ln -s . src/kernel/memory/givaro  ; fi
if test ! -L src/kernel/rational/givaro; then ln -s . src/kernel/rational/givaro; fi
if test ! -L src/kernel/zpz/givaro     ; then ln -s . src/kernel/zpz/givaro     ; fi
if test ! -L src/kernel/gmp++/gmp++    ; then ln -s . src/kernel/gmp++/gmp++    ; fi
if test ! -L src/library/matrix/givaro ; then ln -s . src/library/matrix/givaro ; fi
if test ! -L src/library/poly1/givaro  ; then ln -s . src/library/poly1/givaro  ; fi
if test ! -L src/library/tools/givaro  ; then ln -s . src/library/tools/givaro  ; fi
if test ! -L src/library/vector/givaro ; then ln -s . src/library/vector/givaro ; fi



AC_DEFINE_UNQUOTED(INT8,  $GIVARO_INT8, Canonical 8-bit data type)
AC_DEFINE_UNQUOTED(INT16, $GIVARO_INT16, Canonical 16-bit data type)
AC_DEFINE_UNQUOTED(INT32, $GIVARO_INT32, Canonical 32-bit data type)
AC_DEFINE_UNQUOTED(INT64, $GIVARO_INT64, Canonical 64-bit data type)

echo "-----------------------------------------------"
# Feature checks

AC_LANG_CPLUSPLUS


echo "-----------------------------------------------"

GIV_CHECK_GMP(,,[
echo '*******************************************************************************'
echo ' ERROR: GMP not found!'
echo
echo ' GMP library compiled with --enable-cxx is required for this library to compile.'
echo ' Please make sure GMP is installed and specify its location with the option'
echo ' --with-gmp=<prefix> when running configure.'
echo '*******************************************************************************'
exit 1
])

GIV_DOC

CXXFLAGS="${GMP_CFLAGS} ${CXXFLAGS}"
LIBS="${GMP_LIBS}"

AC_SUBST(CXXFLAGS)

echo "-----------------------------------------------"
echo "          END GIVARO CONFIG                    "
echo "-----------------------------------------------"

AC_CONFIG_FILES([
Makefile
givaro-config
givaro-makefile
src/Makefile
src/kernel/Makefile
src/kernel/system/Makefile
src/kernel/memory/Makefile
src/kernel/gmp++/Makefile
src/kernel/integer/Makefile
src/kernel/bstruct/Makefile
src/kernel/rational/Makefile
src/kernel/zpz/Makefile
src/library/Makefile
src/library/vector/Makefile
src/library/tools/Makefile
src/library/poly1/Makefile
src/library/matrix/Makefile
macros/Makefile
examples/Makefile
examples/FiniteField/Makefile
examples/Integer/Makefile
examples/Rational/Makefile
examples/Polynomial/Makefile
tests/Makefile
docs/Makefile
])
AC_OUTPUT

echo "-----------------------------------------------"
