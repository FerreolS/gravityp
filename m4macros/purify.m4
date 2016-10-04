# PURIFY
#------------------------
# Checks for the INTROOT area
AC_DEFUN([CHECK_PURIFY],
[

    AC_MSG_CHECKING([for PURIFY])

    AC_ARG_ENABLE(purify,
      AC_HELP_STRING([--disable-purify],
        [disalbes the check for the PURIFY installation]),
        enable_purify=$enableval,
        enable_purify=yes)

    if test "x$enable_purify" = xyes ; then
      AC_MSG_RESULT([enabled])
      AC_CHECK_PROG([PURIFY_CMD], [purify], [purify],[NONE])
      if test "$PURIFY_CMD" = "NONE" ; then
        AC_MSG_WARN([Disabling purify support])
        enable_purify=no
      fi
    else
      AC_MSG_RESULT([disabled])
    fi

    AM_CONDITIONAL([PURIFY], [test "x$enable_purify" = "xyes"])
])
