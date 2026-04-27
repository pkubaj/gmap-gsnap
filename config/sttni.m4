
AC_DEFUN([ACX_STTNI], [
  AC_REQUIRE([AC_CANONICAL_HOST])
  AC_LANG_SAVE
  AC_LANG([C])

  AC_MSG_CHECKING(if compiler has STTNI commands)

  case $host_cpu in
    arm*|aarch*)
      ax_cv_sttni=no
    ;;

    i[[3456]]86*|x86_64*|amd64*)
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <stdio.h>
#include <stdlib.h>
#include <nmmintrin.h>]],
                   [[int nshift = rand() % 32;
__m128i a, b, result;
result = _mm_cmpestrm(b, 8, a, 8, _SIDD_UWORD_OPS | _SIDD_CMP_EQUAL_ANY | _SIDD_BIT_MASK);
]])],
  [ax_cv_sttni=no],
  [ax_cv_sttni=yes])
    ;;
  esac

  AC_MSG_RESULT($ax_cv_sttni)
  if test "$ax_cv_sttni" = yes; then
    AC_DEFINE(HAVE_STTNI,1,[Define to 1 if your compiler has STTNI commands.])
  fi

AC_LANG_RESTORE
])dnl

