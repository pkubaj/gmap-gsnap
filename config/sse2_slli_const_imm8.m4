
AC_DEFUN([ACX_SSE2_SLLI_CONST_IMM8], [
  AC_REQUIRE([AC_CANONICAL_HOST])
  AC_LANG_SAVE
  AC_LANG([C])

  AC_MSG_CHECKING(if compiler requires an immediate in sse2 _mm_slli_epi32 command)

  case $host_cpu in
    arm*|aarch*)
      ax_cv_sse2_slli_const_imm8=yes
    ;;

    i[[3456]]86*|x86_64*|amd64*)
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <stdio.h>
#include <stdlib.h>
#include <emmintrin.h>]],
                   [[int nshift = rand() % 32;
__m128i shifted;
shifted = _mm_slli_epi32(_mm_set1_epi32(1),nshift);
]])],
  [ax_cv_sse2_slli_const_imm8=no],
  [ax_cv_sse2_slli_const_imm8=yes])
    ;;
  esac

  AC_MSG_RESULT($ax_cv_sse2_slli_const_imm8)
  if test "$ax_cv_sse2_slli_const_imm8" = yes; then
    AC_DEFINE(SSE2_SLLI_CONST_IMM8,1,[Define to 1 if your compiler requires an immediate in sse2 _mm_slli_epi32 command.])
  fi

AC_LANG_RESTORE
])dnl

