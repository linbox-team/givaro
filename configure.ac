# Copyright(c)'1994-2011 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.
#/

AC_PREREQ([2.68])

AC_INIT([Givaro],[4.1.1],[http://github.com/linbox-team/givaro],[givaro],
		[https://casys.gricad-pages.univ-grenoble-alpes.fr/givaro])

AC_CONFIG_MACRO_DIR([macros])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADERS([config.h])

AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE([1.10 -Wall -Wno-portability foreign])
AX_PREFIX_CONFIG_H(givaro-config.h, __GIVARO)

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

AC_DEBUG
AC_PROFILE
AC_WARNINGS
AC_INLINE # DON'T TOUCH

AS_ECHO([---------------------------------------])

# Set CXXFLAGS to an empty string if not defined, and save it
: ${CXXFLAGS=""}
ORIGINAL_CXXFLAGS="${CXXFLAGS}"

# Look for compiler name to try and pass different flags according to compiler
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

# --enable-debug ?
# __GIVARO_DEBUG will be (un)set in givaro-config.h
#TODO use -fast for icc, -ipa for eko...
AS_IF([test "x$USE_DEBUG" = "xyes"],
      [OPTIM_FLAGS="-O0"
        DEBUG_FLAGS="-Wall -g -UNDEBUG -DDEBUG"],
      [OPTIM_FLAGS="-O2"
        DEBUG_FLAGS="-Wall -DNDEBUG -UDEBUG"]
      )

# --enable-profile ?
AS_IF([test "x$PROF" = "xyes"], [ DEBUG_FLAGS+=" -pg" ])

# --enable-warnings ?
AS_IF([test "x$WARN" = "xyes" -o "x$WARN" = "xfull"],
      [AS_CASE([$CCNAM],
               [eko], [],
               [gcc*|icc*|clang*], [ DEBUG_FLAGS+=" -Wextra" ],
               [AS_BOX([Unsupported compiler ($CCNAM). Please file a bug.],[*])]
               )
      ])

AS_IF([test "x$WARN" = "xfull"],
      [AS_CASE([$CCNAM],
               [eko],
                [],
               [gcc*|icc*|clang*],
                [
                  DEBUG_FLAGS+=" -Wuninitialized -Wconversion -Wcast-qual "
                  DEBUG_FLAGS+=" -pedantic -Wshadow -Wpointer-arith "
                  DEBUG_FLAGS+=" -Wwrite-strings -Wno-long-long"
                  AS_CASE([$CCNAM],
                    [icc],
                        [ DEBUG_FLAGS+=" -Wcheck -ansi" ],
                    [gcc*],
                        [
                         DEBUG_FLAGS+=" -Wno-vla"
                         DEBUG_FLAGS+=" -Wcast-align -Wno-variadic-macros"
                        ],
                    [clang*],
                        [
                         DEBUG_FLAGS+=" -Wno-vla-extension -D__STRICT_ANSI__"
                         DEBUG_FLAGS+=" -Wcast-align -Wno-variadic-macros"
                        ])
                ],
               [AS_BOX([Unsupported compiler ($CCNAM). Please file a bug.],[*])]
               )
      ])


# Add -march=native to OPTIM_FLAGS if not present in CXXFLAGS and
# target==host and --no-marchnative is not set
AC_ARG_WITH(archnative, [AC_HELP_STRING([--without-arch-native],
    [do not use -march=native (default is to use it if -march is not present in
     CXXFLAGS)])])

AS_CASE([$CXXFLAGS],
        [*-march=*], [], # do nothing if already set in CXXFLAGS
        [AS_IF([test "x${with_archnative}" == "xno"],
            [], # do nothing if option is set to no
            [AS_IF([test "${host}" != "${build}" -o "${host}" != "${target}"],
                [AC_MSG_NOTICE("For efficiency you may want to add a '-march=...' flag in CXXFLAGS")],
                [AC_MSG_NOTICE("Adding '-march=native' to flags")
                 OPTIM_FLAGS+=" -march=native"])])])

BACKUP_CXXFLAGS="${CXXFLAGS}"
CXXFLAGS="${OPTIM_FLAGS} ${CXXFLAGS}"
# On i386 and i686, add the -mfpmath=sse flags if SSE is enable
AS_CASE([$target],
        [*i386*|*i686*],
            [AC_RUN_IFELSE([AC_LANG_PROGRAM([[]], [[#ifdef __SSE__
                                                      return 0;
                                                    #else
                                                      return 1;
                                                    #endif
                                                    ]])],
                [AC_MSG_NOTICE("Adding '-mfpmath=sse' to flags")
                 OPTIM_FLAGS+=" -mfpmath=sse"],
                [], # either the flag is not recognized by the compiler or SSE
                    # is not avail => do nothing
                [AC_MSG_NOTICE("If available you may want to add '-mfpmath=sse'
                 to flags")])] # cross-compilation case
            []) # not on i386 nor i686 => do nothing
CXXFLAGS="${BACKUP_CXXFLAGS}"

# With GCC 4.x, the default ABI version is 2. With this version, __m128 and
# __m256 are the same types and therefore we cannot have overloads for both
# types without linking error. It is fixed in ABI version 4.
# FIXME: Why do we set ABI version to 6 ? Why is it not done for GCC 4.9.x ?
AS_IF([test "x$CCNAM" = "xgcc48"],[REQUIRED_FLAGS+=" -fabi-version=6"])

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

AS_ECHO([---------------------------------------])

# Look for GMP and add flags if necessary
GIV_CHECK_GMP(40000)
REQUIRED_FLAGS+=" ${GMP_CFLAGS}"

GIV_DOC

AS_ECHO([---------------------------------------])

CXXFLAGS="${ORIGINAL_CXXFLAGS}"
GIVARO_CXXFLAGS="${REQUIRED_FLAGS} ${OPTIM_FLAGS} ${DEBUG_FLAGS}"
AC_SUBST(GIVARO_CXXFLAGS)
echo "GIVARO_CXXFLAGS = ${GIVARO_CXXFLAGS}"
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
tests/Makefile
docs/Makefile
givaro.pc
])
AC_CONFIG_FILES([givaro-config],[chmod +x givaro-config])
AC_OUTPUT

AS_ECHO([---------------------------------------])
