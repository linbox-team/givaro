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

AM_PROG_AR

AM_ACLOCAL_INCLUDE(macros)

AC_LANG([C++])

echo "-----------------------------------------------"

AC_DEBUG
AC_PROFILE
AC_WARNINGS
AC_INLINE # DON'T TOUCH

echo "-----------------------------------------------"

######################################################
# Try and pass different flags according to compiler #
######################################################


# Set CXXFLAGS to an empty string if not defined
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

#TODO use -fast for icc, -ipa for eko...
if test "x$DBG" = "xyes" ; then
    OPTIM_FLAGS="-O0"
    DEBUG_FLAGS="-Wall -g -UNDEBUG -DGIVARO_DEBUG -DDEBUG"
else
    OMPTIM_FLAGS="-O2"
    DEBUG_FLAGS="-Wall -DNDEBUG -UGIVARO_DEBUG -UDEBUG"
fi

if test "x$PROF" = "xyes" ; then
    DEBUG_FLAGS+=" -pg"
fi

if test "x$WARN" = "xyes" -o "x$WARN" = "xfull" ; then
    case x${CCNAM} in
        xeko)
            ;;
        xgcc|xgcc44|xgcc48|xicc|xclang|xclang38)
            DEBUG_FLAGS+=" -Wextra"
            if test "x${WARN}" = "xfull" ; then
                DEBUG_FLAGS+=" -Wuninitialized -Wconversion -Wcast-qual "
                DEBUG_FLAGS+=" -pedantic -Wshadow -Wpointer-arith "
                DEBUG_FLAGS+=" -Wwrite-strings -Wno-long-long"
                if test "x${CCNAM}" != "xicc" ; then
                    DEBUG_FLAGS+=" -Wcast-align -Wno-variadic-macros"
                fi
                case x${CCNAM} in
                    xicc)
                        DEBUG_FLAGS+=" -Wcheck -ansi"
                        ;;
                    xgcc|xgcc44|xgcc48)
                        DEBUG_FLAGS+=" -Wno-vla"
                        ;;
                    xclang|xclang38)
                        DEBUG_FLAGS+=" -Wno-vla-extension -D__STRICT_ANSI__"
                        ;;
                esac
            fi
            ;;

        *)
            echo
            echo "*******************************************************"
            echo "unsupported compiler ($CCNAM). Please file a bug."
            echo "*******************************************************"
            echo
    esac
fi

AC_SUBST([DEFAULT_CFLAGS])


#
AC_HEADER_STDC
AC_PROG_LIBTOOL

# newer libtool...
LT_PREREQ([2.4.2])
LT_INIT


echo "-----------------------------------------------"
echo "        START  GIVARO CONFIG                   "
echo "-----------------------------------------------"

# Add -march=native to OPTFLAGS if not present in CXXFLAGS and
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

OLDCXXFLAGS="$CXXFLAGS"
CXXFLAGS+=" $OPTIM_FLAGS $DEBUG_FLAGS"
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
CXXFLAGS="$OLDCXXFLAGS"
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
CXXFLAGS="${BACKUP_CXXFLAGS}"

GIV_CHECK_GMP(40000)
REQUIRED_FLAGS="${REQUIRED_FLAGS} ${GMP_CFLAGS}"

GIV_DOC

CXXFLAGS="${CXXFLAGS} ${REQUIRED_FLAGS} ${OPTIM_FLAGS} ${DEBUG_FLAGS}"
LIBS="${GMP_LIBS}"
echo "CXXFLAGS = ${CXXFLAGS}"

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
