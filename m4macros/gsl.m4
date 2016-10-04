# ESO_CHECK_GSL
#------------------
# Checks for the gsl library and header files.
AC_DEFUN([ESO_CHECK_GSL],
[

    AC_ARG_WITH(gsl,
                AC_HELP_STRING([--with-gsl],
                               [location where gsl is installed]),
                [
                    gsl_dir=$withval
                ])

    export PKG_CONFIG_PATH="$gsl_dir/lib/pkgconfig:$gsl_dir/lib64/pkgconfig:$GSLDIR/lib/pkgconfig/:$GSLDIR/lib64/pkgconfig/:$PKG_CONFIG_PATH"

    AC_SUBST(PKG_CONFIG_PATH)
    PKG_CHECK_MODULES([GSL], [gsl], 
                      AC_DEFINE([HAVE_GSL], [1], [1 if GSL present]),
                      AC_MSG_ERROR([No GSL available]))

])
