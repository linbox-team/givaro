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
	AS_IF([ test "x$arch" = "xx86_64" -o "x$arch" = "xi686" ],
		    [archx86="yes"],
		    [archx86="no"]
	     )
	     
        AS_ECHO($CCNAM)
        AS_ECHO($archx86)
	AS_IF([ test  "x$CCNAM" != "xgcc48" -o "x$archx86" = "xno" ],
		[CUSTOM_SIMD="yes"],
		[CUSTOM_SIMD="no"]
	)

	AS_ECHO($CUSTOM_SIMD)

	AC_ARG_ENABLE(avx,[AC_HELP_STRING([--disable-avx], [ Disable AVX 1 instruction set (when available)])],[],[])
	AS_ECHO(enable_avx=$enable_avx)
	if [[ "x$enable_avx" != "xno" ]]; then
	   	AS_ECHO([hehe])
		dnl Autodetection of AVX instruction set enabled
		AC_MSG_CHECKING([for AVX])
		if [[ "x$CUSTOM_SIMD" = "xno" ]]; then
				dnl gcc on x86_64 or i686: using AX_GCC_X86_SUPPORT
				AX_GCC_X86_CPU_SUPPORTS(avx,
				 [SIMD_CFLAGS="$SIMD_CFLAGS -mavx"
				 AC_MSG_RESULT(yes)
				 ],
				 [AC_MSG_RESULT(no)])
		else
			dnl Custom Check for AVX
		
		   	dnl Intel compilers usually do not require option to enable avx
		   	dnl Thus, we test with no option on
			CODE_AVX=`cat macros/CodeChunk/avx.C`
                        BACKUP_CXXFLAGS=${CXXFLAGS}
			for switch_avxflags in "" "-mavx"; do
		       	    CXXFLAGS="${BACKUP_CXXFLAGS} -O0 ${switch_avxflags}"
		       	    AC_TRY_RUN([ ${CODE_AVX} ],
		       	    [
				avx_found="yes"
				AC_MSG_RESULT(yes)
		        	SIMD_CFLAGS="${SIMD_CFLAGS} ${switch_avxflags}"
				break
		       	    ],
		       	    [ avx_found="no"
		            ],
		       	    [
				echo "cross compiling...disabling"
		        	avx_found="no"
		        	break
		            ])
		   	 done
			 AS_IF([ test "x$avx_found" = "xno"],[AC_MSG_RESULT(no)],[])
			 CXXFLAGS=${BACKUP_CXXFLAGS}
		fi
	else
		dnl Autodetection of AVX instruction set disabled
		AS_ECHO("AVX disabled")
	fi
])


AC_DEFUN([GIV_CHECK_SIMD],
[
	AC_ARG_ENABLE(simd,[AC_HELP_STRING([--disable-simd], [ Disable vectorized instructions: SSE4.1, AVX, AVX2])])
	AS_IF([ test  "x$enable_simd" != "xno" ],
	[
		AS_ECHO("SIMD enabled")
		arch=`echo $target | cut -d"-" -f1`
		# if we are on a x86 (32 or 64 bits) with gcc>=4.8 then run the AX_CHECK_X86_FEATURES macro
		AS_IF([test "x$arch" = "xx86_64" -o "x$arch" = "xi686"],
			    [archx86="yes"],
			    [archx86="no"]
		     )
		AS_IF([ test  "x$CCNAM" != "xgcc48" -o "x$archx86" = "xno" ],
		[
		   CUSTOM_SIMD="yes"
		   echo "Compiling with $CCNAM for a $arch target: running custom checks for SSE4.1 and AVX1,2"
		   AC_MSG_CHECKING(for SSE 4.1)
		   BACKUP_CXXFLAGS=${CXXFLAGS}
		   SSEFLAGS="-msse4.1"
		   CXXFLAGS="${BACKUP_CXXFLAGS} ${SSEFLAGS}"
		   CODE_SSE=`cat macros/CodeChunk/sse.C`
		   AC_TRY_RUN([ ${CODE_SSE} ],
			      [ sse_found="yes" ],
			       [ sse_found="no" ],
			       [ 
			       echo "cross compiling...disabling"
				 sse_found="no"
			       ])
	           AS_IF([ test "x$sse_found" = "xyes" ],
		   [
			AC_DEFINE(HAVE_SSE_INSTRUCTIONS,1,[Define if SSE instruction set is available])
			AC_DEFINE(HAVE_SSE2_INSTRUCTIONS,1,[Define if SSE2 instruction set is available])
			AC_DEFINE(HAVE_SSE3_INSTRUCTIONS,1,[Define if SSE3 instruction set is available])
			AC_DEFINE(HAVE_SSE4_1_INSTRUCTIONS,1,[Define if SSE4.1 instruction set is available])
			AC_MSG_RESULT(yes)
                   ],
		   [
			SSEFLAGS=""
			AC_MSG_RESULT(no)
		   ])
		   CXXFLAGS=${BACKUP_CXXFLAGS}
		   
		   dnl Check for AVX
		   AC_MSG_CHECKING(for AVX)
		   CODE_AVX=`cat macros/CodeChunk/avx.C`
		   dnl Intel compilers usually do not require option to enable avx
		   dnl Thus, we test with no option on
		   for switch_avxflags in "" "-mavx"; do
		       CXXFLAGS="${BACKUP_CXXFLAGS} -O0 ${switch_avxflags}"
		       AC_TRY_RUN([ ${CODE_AVX} ],
		       [
				avx_found="yes"
		        	AVXFLAGS=${switch_avxflags}
				break
		       ],
		       [ avx_found="no" ],
		       [
		        echo "cross compiling...disabling"
		        avx_found="no"
		        break
		       ])
		   done
			
		   dnl Is AVX found?
		   AS_IF([ test "x$avx_found" = "xyes" ],
		   [
			AC_MSG_RESULT(yes)
			AC_DEFINE(HAVE_AVX_INSTRUCTIONS,1,[Define if AVX is available])
			
	                dnl Check for AVX2
			AC_MSG_CHECKING(for AVX2)
			for switch_avx2flags in "" "-mfma -mavx2"; do
			    CXXFLAGS="${BACKUP_CXXFLAGS} -O0 ${switch_avx2flags}"
			    AC_TRY_RUN(
			    [
			        #define __try_avx2
				${CODE_AVX}
			    ],
			    [
			        avx2_found="yes"
			        AVX2FLAGS="${switch_avx2flags}"
			        break
		            ],
			    [ avx2_found="no" ],
			    [
			        echo "cross compiling...disabling"
			        avx2_found = "no"
			        break
			    ])
			done
				
	                dnl Is AVX2 found?
			AS_IF([ test "x$avx2_found" = "xyes" ],
			[
				AC_MSG_RESULT(yes)
				AC_DEFINE(HAVE_AVX2_INSTRUCTIONS,1,[Define if AVX2 is available])
				AVXFLAGS=${AVX2FLAGS}
			],
			[ AC_MSG_RESULT(no) ]
			)
		    ],
		    [
			dnl No AVX
		    	AC_MSG_RESULT(no)
		    ])
		
		    CXXFLAGS=${BACKUP_CXXFLAGS}
		],
		[ ])
	],[ AS_ECHO("SIMD disabled")
	    CUSTOM_SIMD="yes" ])
])
