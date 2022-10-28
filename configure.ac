# Copyright(c)'1994-2011 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.
#/

AC_PREREQ([2.68])

AC_INIT([Givaro],[4.2.1],[http://github.com/linbox-team/givaro],[givaro],
		[https://casys.gricad-pages.univ-grenoble-alpes.fr/givaro])

AC_CONFIG_MACRO_DIR([macros])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADERS([config.h])

AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE([1.10 -Wall -Wno-portability foreign])
AX_PREFIX_CONFIG_H(givaro-config.h, __GIVARO)

# Set CXXFLAGS to an empty string if not defined, and save it. We must do it
# before calling AC_PROG_CXX that sets it to '-O2 -g' if not defined
: ${CXXFLAGS=""}
ORIGINAL_CXXFLAGS="${CXXFLAGS}"

# We set the language to C++
AC_LANG([C++])
AC_PROG_CXX

#
AM_PROG_AR

# Libtool
AC_PROG_LIBTOOL
LT_PREREQ([2.4.2])
LT_INIT

# Look for headers
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stddef.h stdlib.h string.h sys/time.h stdint.h pthread.h])

###########################################
AS_BOX([       START GIVARO CONFIG       ])
###########################################

AC_COMPILER_NAME

# We need a C++11 compiler now - AB 2014-12-12
# clang-3.8 does not support __float128 without explicitly passing it -std=c++11
AS_IF([test "x${CCNAM}" = "xclang38"],
      [AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])],
      [AX_CXX_COMPILE_STDCXX_11([ext],[mandatory])]
      )
REQUIRED_FLAGS="${CXX11FLAGS}"
# Add the c++11 flags for the configure compilations as clang needs it to work properly with a glibc++ compiled with gcc
CXXFLAGS="${CXX11FLAGS} ${CXXFLAGS}"

AS_ECHO([---------------------------------------])

AC_INLINE # DON'T TOUCH

# Set OPTIM_FLAGS, DEBUG_FLAGS depending on compiler and command line arguments
SET_FLAGS

# Append -march=native or -mcpu=native (if recognized by the compiler) to
# OPTIM_FLAGS if not present in CXXFLAGS and not cross-compiling and
# --without-archnative is not set
ARCH_FLAGS

# Append -mfpmath=sse to OPTIM_FLAGS on i386 and i686 architecture with SSE
FPMATH_FLAGS

AS_ECHO([---------------------------------------])
# Machine characteristics

# Size of some types
AC_CHECK_SIZEOF(char)
AC_CHECK_SIZEOF(short)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(long long)
AC_CHECK_SIZEOF(__int64_t)

# Looking for int128
AC_CHECK_TYPE([__int128_t], [AC_TRY_COMPILE([#include <type_traits>], [std::make_unsigned<__int128_t>::type y;],[AC_DEFINE(HAVE_INT128, 1, [Define that compiler allows int128_t types])])])


# check endianness of the architecture
AC_C_BIGENDIAN(
  [AC_DEFINE(HAVE_BIG_ENDIAN, 1, [Define that architecture uses big endian storage])],
  [AC_DEFINE(HAVE_LITTLE_ENDIAN, 1, [Define that architecture uses little endian storage])],
  [])

AC_MSG_NOTICE([Creating symbolic links for compilation])
for link in \
	src/kernel/system/givaro   \
	src/kernel/bstruct/givaro  \
	src/kernel/integer/givaro  \
	src/kernel/memory/givaro   \
	src/kernel/rational/givaro \
	src/kernel/recint/recint   \
	src/kernel/field/givaro    \
	src/kernel/ring/givaro     \
	src/kernel/gmp++/gmp++     \
	src/library/matrix/givaro  \
	src/library/poly1/givaro   \
	src/library/tools/givaro   \
	src/library/vector/givaro  \
	; do

	dir="$(dirname "$link")"
	AS_MKDIR_P([$dir])
	case "$srcdir$" in
	    /*) ln -sf          "$srcdir/$dir" "$link" ;;
	     *) ln -sf "../../../$srcdir/$dir" "$link" ;;
	esac
done

AS_ECHO([---------------------------------------])

# Look for GMP and add flags if necessary
GIV_CHECK_GMP(40000)
REQUIRED_FLAGS+=" ${GMP_CFLAGS}"

GIV_DOC

AS_ECHO([---------------------------------------])

CXXFLAGS="${ORIGINAL_CXXFLAGS}"
GIVARO_CXXFLAGS="${REQUIRED_FLAGS} ${OPTIM_FLAGS} ${DEBUG_FLAGS}"
AC_SUBST(GIVARO_CXXFLAGS)
AS_ECHO(["GIVARO_CXXFLAGS = ${GIVARO_CXXFLAGS}"])
AC_SUBST(GMP_LIBS)
AC_SUBST(REQUIRED_FLAGS)

###########################################
AS_BOX([        END GIVARO CONFIG        ])
###########################################

AC_CONFIG_FILES([
Makefile
givaro-makefile
src/Makefile
src/kernel/Makefile
src/kernel/system/Makefile
src/kernel/memory/Makefile
src/kernel/gmp++/Makefile
src/kernel/integer/Makefile
src/kernel/bstruct/Makefile
src/kernel/rational/Makefile
src/kernel/recint/Makefile
src/kernel/field/Makefile
src/kernel/ring/Makefile
src/library/Makefile
src/library/vector/Makefile
src/library/tools/Makefile
src/library/poly1/Makefile
src/library/matrix/Makefile
macros/Makefile
benchmarks/Makefile
examples/Makefile
examples/FiniteField/Makefile
examples/Integer/Makefile
examples/Rational/Makefile
examples/RecInt/Makefile
examples/Polynomial/Makefile
examples/Matrix/Makefile
tests/Makefile
docs/Makefile
givaro.pc
$DOXYGEN_CONFIG_FILES
])
AC_CONFIG_FILES([givaro-config],[chmod +x givaro-config])
AC_OUTPUT

AS_ECHO([---------------------------------------])
