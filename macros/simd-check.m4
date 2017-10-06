dnl Check for SIMD
dnl Created by BB, 2014-03-25
dnl modified by CP, 2016-07-11
dnl copied from fflas-ffpack by CP
dnl ========LICENCE========
dnl Copyright(c)'1994-2016 by The Givaro group
dnl This file is part of Givaro.
dnl Givaro is governed by the CeCILL-B license under French law
dnl and abiding by the rules of distribution of free software.
dnl see the COPYRIGHT file for more details.
dnl ========LICENCE========
dnl

dnl GIV_CHECK_SIMD
dnl
dnl turn on SSE4.1 AVX, AVX2 extensions if available


AC_DEFUN([SIMD_CHECK],
[
        arch=`echo $target | cut -d"-" -f1`
# if we are on a x86 (32 or 64 bits) with gcc>=4.8 then run the AX_CHECK_X86_FEATURES macro
        if [ test "x$arch" = "xx86_64" -o "x$arch" = "xi686" ]; then
                    archx86="yes"
        else
                    archx86="no"
        fi

        if [ test "x$CCNAM" != "xgcc48" -o "x$archx86" = "xno" ]; then
           CUSTOM_SIMD="yes"
        else
           CUSTOM_SIMD="no"
        fi

        AC_ARG_ENABLE(sse2,[AC_HELP_STRING([--disable-sse2], [ "Disable SSE2 instruction set (when available")])],[],[])
        AC_ARG_ENABLE(ssse3,[AC_HELP_STRING([--disable-ssse3], [ "Disable SSSE3 instruction set (when available")])],[],[])
        AC_ARG_ENABLE(sse4.1,[AC_HELP_STRING([--disable-sse4.1], [ "Disable SSE4.1 instruction set (when available")])],[],[])
        AC_ARG_ENABLE(sse4.2,[AC_HELP_STRING([--disable-sse4.2], [ "Disable SSE4.2 instruction set (when available")])],[],[])
        AC_ARG_ENABLE(avx,[AC_HELP_STRING([--disable-avx], [ "Disable AVX instruction set (when available")])],[],[])
        AC_ARG_ENABLE(avx2,[AC_HELP_STRING([--disable-avx2], [ "Disable AVX2 instruction set (when available")])],[],[])
        AC_ARG_ENABLE(fma,[AC_HELP_STRING([--disable-fma], [ "Disable FMA instruction set (when available")])],[],[])
        m4_foreach_w([simd_feature], [sse2 ssse3 sse4.1 sse4.2 avx avx2 fma],
        [
        AS_IF([ test "x$enable_[]simd_feature" != "xno" ],[
                dnl Autodetection of AVX instruction set enabled
                if [[ "x$CUSTOM_SIMD" = "xno" ]]; then
                                dnl gcc on x86_64 or i686: using AX_GCC_X86_SUPPORT
                                AX_GCC_X86_CPU_SUPPORTS(simd_feature,
                                 [SIMD_CFLAGS="$SIMD_CFLAGS -m[]simd_feature"], [])
                else
                        dnl Custom Check
                        dnl Intel compilers usually do not require option to enable avx
                        dnl Thus, we test with no option on
                        AC_MSG_CHECKING([for simd_feature])
                        CODE=`cat macros/CodeChunk/simd_feature.C`
                        BACKUP_CXXFLAGS=${CXXFLAGS}
                        dnl feature may be enabled by default, hence also testing with no flags (e.g. sse2 with gcc on x86_64) 
                        for switch_flags in "" "-m[]simd_feature"; do
                            CXXFLAGS="${BACKUP_CXXFLAGS} -O0 ${switch_flags}"
                            AC_TRY_RUN([ ${CODE} ],
                            [
                                feature_found="yes"
                                SIMD_CFLAGS="${SIMD_CFLAGS} ${switch_flags}"
                                break                       
                            ],
                            [ feature_found="no"
                            ],
                            [
                                echo "cross compiling...disabling"
                                feature_found="no"
                                break
                            ])
                         done
                         AS_IF([ test "x$feature_found" = "xno"],[AC_MSG_RESULT(no)],[AC_MSG_RESULT(yes)])
                         CXXFLAGS=${BACKUP_CXXFLAGS}
                fi
        ],[
                dnl Autodetection of AVX instruction set disabled
                AS_ECHO("simd_feature disabled")
        ])
        ])
])
