# Copyright(c)'1994-2011 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.
#/

AC_PREREQ([2.68])

AC_INIT([Givaro],[4.0.4],[http://github.com/linbox-team/givaro],[givaro],
		[http://givaro.forge.imag.fr/])

AC_CONFIG_MACRO_DIR([macros])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADERS([config.h])

AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE([1.10 -Wall -Wno-portability foreign])
AX_PREFIX_CONFIG_H(givaro-config.h, __GIVARO)

AM_PROG_AR

AM_ACLOCAL_INCLUDE(macros)

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

# We need a C++11 compiler now - AB 2014-12-12
# clang-3.8 does not support __float128 without explicitly passing it -std=c++11
if test "x${CCNAM}" = "xclang38" ; then
   AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])
else
   AX_CXX_COMPILE_STDCXX_11([ext],[mandatory])
fi
REQUIRED_FLAGS="$CXX11FLAGS"
# Add the c++11 flags for the configure compilations as clang needs it to work properly with a glibc++ compiled with gcc
BACKUP_CXXFLAGS="${CXXFLAGS}"
CXXFLAGS="${CXXFLAGS} ${REQUIRED_FLAGS}"
TESTS_CFLAGS="-O0"
DEBUG_CFLAGS="-g"
DEFAULT_CFLAGS=""
WARN_CFLAGS="-Wall"

#TODO use -fast for icc, -ipa for eko...
if test "x$DBG" = "xyes" ; then
	DEFAULT_CFLAGS="-O0 ${DEFAULT_CFLAGS} " #those are CXXFLAGS
	DEBUG_CFLAGS="${DEBUG_CFLAGS} -UNDEBUG -DDEBUG -DGIVARO_DEBUG"
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
		xclang38)
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

AC_SUBST([DEFAULT_CFLAGS])
AC_SUBST([DEBUG_CFLAGS])
AC_SUBST([TESTS_CFLAGS])

AC_HEADER_STDC
AC_PROG_LIBTOOL

# newer libtool...
LT_PREREQ([2.4.2])
LT_INIT


echo "-----------------------------------------------"
echo "        START  GIVARO CONFIG                   "
echo "-----------------------------------------------"

# Checks which SIMD optimization flags to use and set SIMD_CFLAGS accordingly
INSTR_SET
AC_SUBST(SIMD_FLAGS)

# Veltkamp Split in modular-extended.inl requires to use the SSE floating point multiplier and not the fp87
# This is by default on x86_64 but not on i386.
arch=`echo $target | cut -d"-" -f1`
AS_IF([ test \( "x$arch" = "xi386" -o "x$arch" = "xi686" \) -a "x$HAVE_SSE" = "xyes" ],
   [
        AS_ECHO("Compiling for an $arch architecture: appending -mfpmath=sse flags")
        REQUIRED_FLAGS="${REQUIRED_FLAGS} -mfpmath=sse"
   ])
	
echo "-----------------------------------------------"

dnl With GCC's default ABI version, a __m128 or __m256 are the same types and therefore we cannot
dnl have overloads for both types without linking error.
AS_IF([test "x$CCNAM" = "xgcc48"],[REQUIRED_FLAGS="${REQUIRED_FLAGS} -fabi-version=6"],[])

AC_SUBST(REQUIRED_FLAGS)

# Machine characteristics

AC_CHECK_SIZEOF(char, 8)
AC_CHECK_SIZEOF(short, 16)
AC_CHECK_SIZEOF(int, 32)
AC_CHECK_SIZEOF(long, 32)
AC_CHECK_SIZEOF(long long, 64)
AC_CHECK_SIZEOF(__int64_t, 64)

AC_LANG([C++])

AC_CHECK_TYPE([__int128_t], [AC_TRY_COMPILE([#include <type_traits>], [std::make_unsigned<__int128_t>::type y;],[AC_DEFINE(HAVE_INT128, 1, [Define that compiler allows int128_t types])])])		   



# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stddef.h stdlib.h string.h sys/time.h stdint.h pthread.h])


# check endianness of the architecture
AC_C_BIGENDIAN(
  [AC_DEFINE(HAVE_BIG_ENDIAN, 1, [Define that architecture uses big endian storage])],
  [AC_DEFINE(HAVE_LITTLE_ENDIAN, 1, [Define that architecture uses little endian storage])],
  [])

AC_MSG_NOTICE([Creating symbolic link for compilation])
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
	AS_IF([ test ! -L "$link" ],
	  [ln -s "../../../$srcdir/$dir" "$link"])
done

AC_DEFINE_UNQUOTED(INT8,  $GIVARO_INT8, Canonical 8-bit data type)
AC_DEFINE_UNQUOTED(INT16, $GIVARO_INT16, Canonical 16-bit data type)
AC_DEFINE_UNQUOTED(INT32, $GIVARO_INT32, Canonical 32-bit data type)
AC_DEFINE_UNQUOTED(INT64, $GIVARO_INT64, Canonical 64-bit data type)

echo "-----------------------------------------------"
# Feature checks


echo "-----------------------------------------------"

GIV_CHECK_GMP(40000)

GIV_DOC
CXXFLAGS=${BACKUP_CXXFLAGS}
CXXFLAGS="${CXXFLAGS} ${GMP_CFLAGS}"
CXXFLAGS="${CXXFLAGS} ${REQUIRED_FLAGS}"
CXXFLAGS="${CXXFLAGS} ${SIMD_CFLAGS}"
LIBS="${GMP_LIBS}"
echo "CXXFLAGS=${CXXFLAGS}"

AC_SUBST(CXXFLAGS)

echo "-----------------------------------------------"
echo "          END GIVARO CONFIG                    "
echo "-----------------------------------------------"

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
macros/CodeChunk/Makefile
benchmarks/Makefile
examples/Makefile
examples/FiniteField/Makefile
examples/Integer/Makefile
examples/Rational/Makefile
examples/RecInt/Makefile
examples/Polynomial/Makefile
tests/Makefile
docs/Makefile
givaro.pc
])
AC_CONFIG_FILES([givaro-config],[chmod +x givaro-config])
AC_OUTPUT

echo "-----------------------------------------------"
