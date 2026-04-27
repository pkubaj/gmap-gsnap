# ===========================================================================
#          http://www.gnu.org/software/autoconf-archive/ax_ext.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_EXT
#
# DESCRIPTION
#
#   Find supported SIMD extensions by requesting cpuid. When an SIMD
#   extension is found, the -m"simdextensionname" is added to SIMD_CFLAGS if
#   compiler supports it. For example, if "sse2" is available, then "-msse2"
#   is added to SIMD_CFLAGS.
#
#   This macro calls:
#
#     AC_SUBST(SIMD_CFLAGS)
#
#   And previously defined:
#
#     HAVE_MMX / HAVE_SSE / HAVE_SSE2 / HAVE_SSE3 / HAVE_SSSE3 / HAVE_SSE4_1 / HAVE_SSE4_2
#
#   But now defines:
#
#     HAVE_MM_POPCNT / HAVE_POPCNT / HAVE_LZCNT / HAVE_TZCNT / HAVE_PEXT
#
#
# LICENSE
#
#   Copyright (c) 2007 Christophe Tournayre <turn3r@users.sourceforge.net>
#   Copyright (c) 2013 Michael Petch <mpetch@capp-sysware.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#   Assumes AX_COMPILER_VENDOR has been called to define ax_cv_c_compiler_vendor

#serial 13

AC_DEFUN([AX_EXT],
[
  AC_REQUIRE([AC_CANONICAL_HOST])

  CFLAGS_ORIG=$CFLAGS

  case $host_cpu in
    arm*|aarch*)
        ax_make_arm=yes
        SIMD_ARM_CFLAGS=
    ;;

    i[[3456]]86*|x86_64*|amd64*)

      # AX_CHECK_COMPILE_FLAG tests CFLAGS EXTRA-FLAGS FLAG

      if test x"$ax_cv_cpu_has_sse2_ext" = xyes; then
        #AC_MSG_CHECKING(for sse2 compiler/linker support)
        CFLAGS=
        TEST_CFLAGS=-msse2
	AX_CHECK_COMPILE_FLAG([$TEST_CFLAGS], [ax_cv_compile_sse2_ext=yes], [ax_cv_ext_compile_problem=yes])
	if test x"$ax_cv_compile_sse2_ext" != xyes; then
	  AC_MSG_WARN([Your CPU supports SSE2 instructions but not your compiler.  Can you try another compiler or update yours?])
	else
          CFLAGS=$TEST_CFLAGS
          AC_MSG_CHECKING(for emmintrin.h header file)
          AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <emmintrin.h>])],
                         [ax_cv_link_emmintrin_h=yes],
		         [ax_cv_ext_linker_problem=yes])
          if test x"$ax_cv_link_emmintrin_h" != xyes; then
            AC_MSG_WARN([Your compiler supports SSE2 instructions but not your linker.])
          else
            AC_MSG_RESULT([yes])
            ax_make_sse2=yes
            SIMD_SSE2_CFLAGS=$CFLAGS
#  	    AC_DEFINE(HAVE_SSE2,1,[Define to 1 if you support SSE2 (Streaming SIMD Extensions 2) instructions]) -- Not used
          fi            
	fi
      fi


      if test x"$ax_cv_cpu_has_ssse3_ext" = xyes; then
        #AC_MSG_CHECKING(for ssse3 compiler/linker support)
        CFLAGS=
        if test x"$ax_cv_c_compiler_vendor" = xintel; then
          TEST_CFLAGS="-mssse3"
	else
          TEST_CFLAGS="$SIMD_SSE2_CFLAGS -mssse3"
        fi
        AX_CHECK_COMPILE_FLAG([$TEST_CFLAGS], [ax_cv_compile_ssse3_ext=yes], [ax_cv_ext_compile_problem=yes])
        if test x"$ax_cv_compile_ssse3_ext" != xyes; then
          AC_MSG_WARN([Your CPU supports SSSE3 instructions but not your compiler.  Can you try another compiler or update yours?])
        else
          CFLAGS=$TEST_CFLAGS
          AC_MSG_CHECKING(for tmmintrin.h header file)
          AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <tmmintrin.h>])],
                         [ax_cv_link_tmmintrin_h=yes],
		         [ax_cv_ext_linker_problem=yes])
          if test x"$ax_cv_link_tmmintrin_h" != xyes; then
            AC_MSG_WARN([Your compiler supports SSSE3 instructions but not your linker.])
          else
            AC_MSG_RESULT([yes])
            ax_make_ssse3=yes
	    SIMD_SSSE3_CFLAGS=$CFLAGS
            if test x"$ax_cv_c_compiler_vendor" != xintel; then
              SIMD_SSE2_CFLAGS="$SIMD_SSE2_CFLAGS -mno-ssse3"
            fi
#           AC_DEFINE(HAVE_SSSE3,1,[Define to 1 if you support SSSE3 (Supplemental Streaming SIMD Extensions 3) instructions]) -- Defines run-type
          fi            
        fi
      fi


      if test x"$ax_cv_cpu_has_sse41_ext" = xyes; then
        #AC_MSG_CHECKING(for sse41 compiler/linker support)
        CFLAGS=
        if test x"$ax_cv_c_compiler_vendor" = xintel; then
          TEST_CFLAGS="-msse4.1"
	else
          TEST_CFLAGS="$SIMD_SSSE3_CFLAGS -msse4.1"
        fi
	AX_CHECK_COMPILE_FLAG([$TEST_CFLAGS], [ax_cv_compile_sse41_ext=yes], [ax_cv_ext_compile_problem=yes])
	if test x"$ax_cv_compile_sse41_ext" != xyes; then
	  AC_MSG_WARN([Your CPU supports SSE4.1 instructions but not your compiler.  Can you try another compiler or update yours?])
	else
          CFLAGS=$TEST_CFLAGS
          AC_MSG_CHECKING(for smmintrin.h header file)
          AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <smmintrin.h>])],
                         [ax_cv_link_smmintrin_h=yes],
		         [ax_cv_ext_linker_problem=yes])
          if test x"$ax_cv_link_smmintrin_h" != xyes; then
            AC_MSG_WARN([Your compiler supports SSE4.1 instructions but not your linker.])
          else
            AC_MSG_RESULT([yes])
            ax_make_sse41=yes
            SIMD_SSE4_1_CFLAGS=$CFLAGS
            if test x"$ax_cv_c_compiler_vendor" != xintel; then
  	      SIMD_SSE2_CFLAGS="$SIMD_SSE2_CFLAGS -mno-sse4.1"
	      SIMD_SSSE3_CFLAGS="$SIMD_SSSE3_CFLAGS -mno-sse4.1"
            fi
#  	    AC_DEFINE(HAVE_SSE4_1,1,[Define to 1 if you support SSE4.1 (Streaming SIMD Extensions 4.1) instructions]) -- Not used
          fi            
	fi
      fi


      if test x"$ax_cv_cpu_has_sse42_ext" = xyes; then
        #AC_MSG_CHECKING(for sse42 compiler/linker support)
        CFLAGS=
        if test x"$ax_cv_c_compiler_vendor" = xintel; then
          TEST_CFLAGS="-march=corei7"
        else
          TEST_CFLAGS="$SIMD_SSE4_1_CFLAGS -msse4.2"
        fi
        AX_CHECK_COMPILE_FLAG([$TEST_CFLAGS], [ax_cv_compile_sse42_ext=yes], [ax_cv_ext_compile_problem=yes])
        if test x"$ax_cv_compile_sse42_ext" != xyes; then
          AC_MSG_WARN([Your CPU supports SSE4.2 instructions but not your compiler.  Can you try another compiler or update yours?])
        else
          CFLAGS=$TEST_CFLAGS
          AC_MSG_CHECKING(for nmmintrin.h header file)
          AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <nmmintrin.h>])],
                         [ax_cv_link_nmmintrin_h=yes],
		         [ax_cv_ext_linker_problem=yes])
          if test x"$ax_cv_link_nmmintrin_h" != xyes; then
            AC_MSG_WARN([Your compiler supports SSE4.2 instructions but not your linker.])
          else
            AC_MSG_RESULT([yes])
            ax_make_sse42=yes
            SIMD_SSE4_2_CFLAGS=$CFLAGS
            if test x"$ax_cv_c_compiler_vendor" != xintel; then
	      SIMD_SSE2_CFLAGS="$SIMD_SSE2_CFLAGS -mno-sse4.2"
	      SIMD_SSSE3_CFLAGS="$SIMD_SSSE3_CFLAGS -mno-sse4.2"
	      SIMD_SSE4_1_CFLAGS="$SIMD_SSE4_1_CFLAGS -mno-sse4.2"
            fi
          fi
        fi
      fi


      # Extra options for SSE42
      if test x"$ax_make_sse42" = xyes; then
        if test x"$ax_cv_cpu_has_popcnt_ext" != xyes; then
          AC_MSG_RESULT([no])
        else
          AC_MSG_CHECKING(for mm_popcnt run support)
          AC_RUN_IFELSE(
            [AC_LANG_PROGRAM([[#include <nmmintrin.h>]],
		             [[return (_mm_popcnt_u32(0xffffffffu) == 32) ? 0 : 9;]])],
            [AC_MSG_RESULT(yes)
   	     AC_DEFINE(HAVE_MM_POPCNT,1,[Define to 1 if you support Intel intrinsic _mm_popcnt_u32/64 instructions])],
            [AC_MSG_RESULT(no)])
        fi


        if test x"$ax_cv_cpu_has_bmi1_ext" = xyes; then
           CFLAGS=
	  if test x"$ax_cv_c_compiler_vendor" = xintel; then
	    TEST_CFLAGS="$SIMD_SSE4_2_CFLAGS"
	  else
	    TEST_CFLAGS="$SIMD_SSE4_2_CFLAGS -mabm"
	  fi  
	  AX_CHECK_COMPILE_FLAG([$TEST_CFLAGS], [ax_cv_compile_abm_ext=yes], [ax_cv_ext_compile_problem=yes])
	  if test x"$ax_cv_compile_abm_ext" != xyes; then
	    AC_MSG_WARN([Your CPU supports BMI1 instructions but not your compiler.  Can you try another compiler or update yours?])
	  else
	    CFLAGS=$TEST_CFLAGS
            AC_MSG_CHECKING(for abm linker support)
	    AC_LINK_IFELSE([AC_LANG_PROGRAM([])],
			   [ax_cv_link_abm=yes],
			   [ax_cv_ext_linker_problem=yes])
	    if test x"$ax_cv_link_abm" != xyes; then
	      AC_MSG_WARN([Your compiler supports -mabm but not your linker])
	    else
              AC_MSG_RESULT([yes])
	      SIMD_SSE4_2_CFLAGS=$CFLAGS

	      # Test for functionality with the -mabm compiler flag
	      AC_MSG_CHECKING(for _popcnt32 run support)
	      AC_RUN_IFELSE(
		[AC_LANG_PROGRAM([[#include <immintrin.h>]],
				 [[return (_popcnt32(0xffffffffu) == 32) ? 0 : 9;]])],
		[AC_MSG_RESULT(yes)
		 AC_DEFINE(HAVE_POPCNT,1,[Define to 1 if you support Intel intrinsic _popcnt instruction])],
		[AC_MSG_RESULT(no)])

	      AC_MSG_CHECKING(for _lzcnt run support)
	      AC_RUN_IFELSE(
		[AC_LANG_PROGRAM([[#include <immintrin.h>]],
				 [[return (_lzcnt_u32(0x0fffffffu) == 4) ? 0 : 9;]])],
		[AC_MSG_RESULT(yes)
		 AC_DEFINE(HAVE_LZCNT,1,[Define to 1 if you support Intel intrinsic _lzcnt instruction])],
		[AC_MSG_RESULT(no)])
	    fi
	  fi
        fi


	if test x"$ax_cv_cpu_has_bmi1_ext" = xyes; then
          #AC_MSG_CHECKING(for bmi compiler/linker support)
	  CFLAGS=
	  if test x"$ax_cv_c_compiler_vendor" = xintel; then
	    TEST_CFLAGS="$SIMD_SSE4_2_CFLAGS"
	  else
	    TEST_CFLAGS="$SIMD_SSE4_2_CFLAGS -mbmi"
	  fi  
	  AX_CHECK_COMPILE_FLAG([$TEST_CFLAGS], [ax_cv_compile_bmi_ext=yes], [ax_cv_ext_compile_problem=yes])
	  if test x"$ax_cv_compile_bmi_ext" != xyes; then
	    AC_MSG_WARN([Your CPU supports BMI1 instructions but not your compiler.  Can you try another compiler or update yours?])
	  else
	    CFLAGS=$TEST_CFLAGS
            AC_MSG_CHECKING(for bmi linker support)
	    AC_LINK_IFELSE([AC_LANG_PROGRAM([])],
			   [ax_cv_link_bmi=yes],
			   [ax_cv_ext_linker_problem=yes])
	    if test x"$ax_cv_link_bmi" != xyes; then
	      AC_MSG_WARN([Your compiler supports -mbmi but not your linker])
	    else
              AC_MSG_RESULT([yes])
	      SIMD_SSE4_2_CFLAGS=$CFLAGS

	      # Test for functionality with the -mbmi compiler flag
	      AC_MSG_CHECKING(for _tzcnt run support)
	      AC_RUN_IFELSE(
		[AC_LANG_PROGRAM([[#include <immintrin.h>]],
				 [[return (_tzcnt_u32(0xfffffff0u) == 4) ? 0 : 9;]])],
		[AC_MSG_RESULT(yes)
		 AC_DEFINE(HAVE_TZCNT,1,[Define to 1 if you support Intel intrinsic _tzcnt instruction])],
		[AC_MSG_RESULT(no)])
	    fi
	  fi
	fi
      fi


      
      if test x"$ax_cv_cpu_has_avx2_ext" = xyes; then
        #AC_MSG_CHECKING(for avx2 compiler/linker support)
        CFLAGS=
	if test x"$ax_cv_c_compiler_vendor" = xintel; then
	  TEST_CFLAGS="-march=core-avx2"
	else
          TEST_CFLAGS="$SIMD_SSE4_2_CFLAGS -mavx2"
        fi
        AX_CHECK_COMPILE_FLAG([$TEST_CFLAGS], [ax_cv_compile_avx2_ext=yes], [ax_cv_ext_compile_problem=yes])
        if test x"$ax_cv_compile_avx2_ext" != xyes; then
          AC_MSG_WARN([Your CPU supports AVX2 instructions but not your compiler.  Can you try another compiler or update yours?])
        else
          CFLAGS=$TEST_CFLAGS
          AC_MSG_CHECKING(for immintrin.h header file)
          AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <immintrin.h>])],
                         [ax_cv_link_immintrin_h=yes],
		         [ax_cv_ext_linker_problem=yes])
          if test x"$ax_cv_link_immintrin_h" != xyes; then
            AC_MSG_WARN([Your compiler supports AVX2 instructions but not your linker.])
          else
            AC_MSG_RESULT([yes])
            ax_make_avx2=yes
            SIMD_AVX2_CFLAGS=$CFLAGS
            if test x"$ax_cv_c_compiler_vendor" != xintel; then
  	      SIMD_SSE2_CFLAGS="$SIMD_SSE2_CFLAGS -mno-avx2"
	      SIMD_SSSE3_CFLAGS="$SIMD_SSSE3_CFLAGS -mno-avx2"
	      SIMD_SSE4_1_CFLAGS="$SIMD_SSE4_1_CFLAGS -mno-avx2"
	      SIMD_SSE4_2_CFLAGS="$SIMD_SSE4_2_CFLAGS -mno-avx2"
            fi
#  	    AC_DEFINE(HAVE_AVX2,1,[Define to 1 if you support AVX2 (Advanced Vector Extensions 2) instructions]) -- Defines run-type
          fi
        fi
      fi


      # Extra option for AVX2
      if test x"$ax_make_avx2" = xyes; then
        if test x"$ax_cv_cpu_has_bmi2_ext" = xyes; then
	  CFLAGS=      
	  if test x"$ax_cv_c_compiler_vendor" = xintel; then
	    TEST_CFLAGS="$SIMD_AVX2_CFLAGS"
	  else
	    TEST_CFLAGS="$SIMD_AVX2_CFLAGS -mbmi2"
	  fi
	  AX_CHECK_COMPILE_FLAG([$TEST_CFLAGS], [ax_cv_compile_bmi2_ext=yes], [ax_cv_ext_compile_problem=yes])
	  if test x"$ax_cv_compile_bmi2_ext" != xyes; then
	    AC_MSG_WARN([Your CPU supports BMI2 instructions but not your compiler.  Can you try another compiler or update yours?])
	  else
	    CFLAGS=$TEST_CFLAGS
            AC_MSG_CHECKING(for bmi2 linker support)
	    AC_LINK_IFELSE([AC_LANG_PROGRAM([])],
			   [ax_cv_link_bmi2=yes],
			   [ax_cv_ext_linker_problem=yes])
	    if test x"$ax_cv_link_bmi2" != xyes; then
	      AC_MSG_WARN([Your compiler supports -mbmi2 but not your linker.  Can you try another linker or update yours?])
	    else
              AC_MSG_RESULT([yes])
	      SIMD_AVX2_CFLAGS=$CFLAGS
	      if test x"$ax_cv_c_compiler_vendor" != xintel; then
		SIMD_SSE2_CFLAGS="$SIMD_SSE2_CFLAGS -mno-bmi2"
		SIMD_SSSE3_CFLAGS="$SIMD_SSSE3_CFLAGS -mno-bmi2"
		SIMD_SSE4_1_CFLAGS="$SIMD_SSE4_1_CFLAGS -mno-bmi2"
		SIMD_SSE4_2_CFLAGS="$SIMD_SSE4_2_CFLAGS -mno-bmi2"
	      fi

	      AC_MSG_CHECKING(for _pext run support)
	      AC_RUN_IFELSE(
		[AC_LANG_PROGRAM([[#include <immintrin.h>]],
				 [[return (_pext_u32(0x01234567u,0x0f33aa55) == 0x00001b0b) ? 0 : 9;]])],
		[AC_MSG_RESULT(yes)
		 AC_DEFINE(HAVE_PEXT,1,[Define to 1 if you support Intel intrinsic _pext instruction])],
		[AC_MSG_RESULT(no)])
	    fi
	  fi
	fi
      fi


      if test x"$ax_cv_cpu_has_avx512_ext" = xyes; then
        #AC_MSG_CHECKING(for avx512f and avx512cd support)
        CFLAGS=
        if test x"$ax_cv_c_compiler_vendor" = xintel; then
          TEST_CFLAGS="-xCOMMON-AVX512"
        else
          TEST_CFLAGS="$SIMD_AVX2_CFLAGS -mavx512f -mavx512cd"
        fi
        AX_CHECK_COMPILE_FLAG([$TEST_CFLAGS], [ax_cv_compile_avx512_ext=yes], [ax_cv_ext_compile_problem=yes])
        if test x"$ax_cv_compile_avx512_ext" != xyes; then
          AC_MSG_WARN([Your CPU supports AVX512 instructions but not your compiler.  Can you try another compiler or update yours?])
        else
          CFLAGS=$TEST_CFLAGS
          AC_MSG_CHECKING(for nmmintrin.h header file)
          AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <nmmintrin.h>])],
                         [ax_cv_link_nmmintrin_h=yes],
		         [ax_cv_ext_linker_problem=yes])
          if test x"$ax_cv_link_nmmintrin_h" != xyes; then
            AC_MSG_WARN([Your compiler supports AVX512 instructions but not your linker.  Can you try another linker or update yours?])
          else
            AC_MSG_RESULT([yes])
            ax_make_avx512=yes
            SIMD_AVX512_CFLAGS=$CFLAGS
            if test x"$ax_cv_c_compiler_vendor" != xintel; then
	      SIMD_SSE2_CFLAGS="$SIMD_SSE2_CFLAGS -mno-avx512f -mno-avx512cd"
	      SIMD_SSSE3_CFLAGS="$SIMD_SSSE3_CFLAGS -mno-avx512f -mno-avx512cd"
	      SIMD_SSE4_1_CFLAGS="$SIMD_SSE4_1_CFLAGS -mno-avx512f -mno-avx512cd"
	      SIMD_SSE4_2_CFLAGS="$SIMD_SSE4_2_CFLAGS -mno-avx512f -mno-avx512cd"
	      SIMD_AVX2_CFLAGS="$SIMD_AVX2_CFLAGS -mno-avx512f -mno-avx512cd"
            fi
#  	    AC_DEFINE(HAVE_AVX512,1,[Define to 1 if you support AVX512 (Advanced Vector Extensions 512) instructions]) -- Defines run-type
#  	    AC_DEFINE(HAVE_AVX512BW,1,[Define to 1 if you support AVX512BW (Advanced Vector Extensions 512) BW instructions]) -- Defines run-type
          fi
        fi
      fi

      # Extra options for AVX512
      if test x"$ax_make_avx512" = xyes; then
        CFLAGS=
        if test x"$ax_cv_c_compiler_vendor" = xintel; then
          TEST_CFLAGS="$SIMD_AVX512_CFLAGS"
        else
          TEST_CFLAGS="$SIMD_AVX512_CFLAGS -mavx512vl -mavx512bw"
        fi
        AX_CHECK_COMPILE_FLAG([$TEST_CFLAGS], [ax_cv_compile_avx512bw_ext=yes], [ax_cv_ext_compile_problem=yes])
        if test x"$ax_cv_compile_avx512bw_ext" != xyes; then
          AC_MSG_WARN([Your CPU supports AVX512BW instructions but not your compiler.  Can you try another compiler or update yours?])
        else
          CFLAGS=$TEST_CFLAGS
          AC_MSG_CHECKING(for avx512vl and avx512bw linker support)
          AC_LINK_IFELSE([AC_LANG_PROGRAM([])],
                         [ax_cv_link_avx512bw=yes],
		         [ax_cv_ext_linker_problem=yes])
          if test x"$ax_cv_link_avx512bw" != xyes; then
            AC_MSG_WARN([Your compiler supports AVX512BW instructions but not your linker.])
          else
            AC_MSG_RESULT([yes])
            SIMD_AVX512_CFLAGS=$CFLAGS
            if test x"$ax_cv_c_compiler_vendor" != xintel; then
	      SIMD_SSE2_CFLAGS="$SIMD_SSE2_CFLAGS -mno-avx512vl -mno-avx512bw"
	      SIMD_SSSE3_CFLAGS="$SIMD_SSSE3_CFLAGS -mno-avx512vl -mno-avx512bw"
	      SIMD_SSE4_1_CFLAGS="$SIMD_SSE4_1_CFLAGS -mno-avx512vl -mno-avx512bw"
	      SIMD_SSE4_2_CFLAGS="$SIMD_SSE4_2_CFLAGS -mno-avx512vl -mno-avx512bw"
	      SIMD_AVX2_CFLAGS="$SIMD_AVX2_CFLAGS -mno-avx512vl -mno-avx512bw"
            fi
#  	    AC_DEFINE(HAVE_AVX512,1,[Define to 1 if you support AVX512 (Advanced Vector Extensions 512) instructions]) -- Defines run-type
#  	    AC_DEFINE(HAVE_AVX512BW,1,[Define to 1 if you support AVX512BW (Advanced Vector Extensions 512) BW instructions]) -- Defines run-type
          fi
        fi
      fi

    ;;
  esac

  CFLAGS=$CFLAGS_ORIG

  AC_SUBST(SIMD_SSE2_CFLAGS)
  AC_SUBST(SIMD_SSSE3_CFLAGS)
  AC_SUBST(SIMD_SSE4_1_CFLAGS)
  AC_SUBST(SIMD_SSE4_2_CFLAGS)
  AC_SUBST(SIMD_AVX2_CFLAGS)
  AC_SUBST(SIMD_AVX512_CFLAGS)

])

