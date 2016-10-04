
# CPL_CHECK_CEXT
#---------------
# Checks for the C extension library and header files.
AC_DEFUN([CPL_CHECK_CEXT],
[

    AC_MSG_CHECKING([for libcext])

    cpl_cext_check_header="cxmacros.h"
    cpl_cext_check_lib="libcext.la"

    cpl_cext_includes=""
    cpl_cext_libraries=""

    AC_ARG_WITH(cext,
                AC_HELP_STRING([--with-cext],
                               [location where Cext is installed]),
                [
                    cpl_with_cext_includes=$withval/include
                    cpl_with_cext_libs=$withval/lib
                ])

    AC_ARG_WITH(cext-includes,
                AC_HELP_STRING([--with-cext-includes],
                               [location of the Cext header files]),
                cpl_with_cext_includes=$withval)

    AC_ARG_WITH(cext-libs,
                AC_HELP_STRING([--with-cext-libs],
                               [location of the Cext library]),
                cpl_with_cext_libs=$withval)

    AC_ARG_ENABLE(cext-test,
                  AC_HELP_STRING([--disable-cext-test],
                                 [disables checks for the Cext library and headers]),
                  cpl_enable_cext_test=$enableval,
                  cpl_enable_cext_test=yes)


    if test "x$cpl_enable_cext_test" = xyes; then

        # Check for the Cext includes

        if test -z "$cpl_with_cext_includes"; then
            cpl_cext_incdirs="/usr/local/include /usr/include"

            test -n "$CPLDIR" && cpl_cext_incdirs="$CPLDIR/include \
                                                   $cpl_cext_incdirs"
        else
            cpl_cext_incdirs="$cpl_with_cext_includes"
        fi

        ESO_FIND_FILE($cpl_cext_check_header, $cpl_cext_incdirs,
                      cpl_cext_includes)


        # Check for the Cext library

        if test -z "$cpl_with_cext_libs"; then
            cpl_cext_libdirs="/usr/local/lib /usr/lib"

            test -n "$CPLDIR" && cpl_cext_libdirs="$CPLDIR/lib \
                                                   $cpl_cext_libdirs"
        else
            cpl_cext_libdirs="$cpl_with_cext_libs"
        fi

        ESO_FIND_FILE($cpl_cext_check_lib, $cpl_cext_libdirs,
                      cpl_cext_libraries)


        if test x"$cpl_cext_includes" = xno || \
            test x"$cpl_cext_libraries" = xno; then
            cpl_cext_notfound=""

            if test x"$cpl_cext_includes" = xno; then
                if test x"$cpl_cext_libraries" = xno; then
                    cpl_cext_notfound="(headers and libraries)"
                else            
                    cpl_cext_notfound="(headers)"
                fi
            else
                cpl_cext_notfound="(libraries)"
            fi

            AC_MSG_ERROR([Cext $cpl_cext_notfound was not found on your system. Please check!])
        else
            AC_MSG_RESULT([libraries $cpl_cext_libraries, headers $cpl_cext_includes])
        fi


        # Set up the symbols

        CX_INCLUDES="-I$cpl_cext_includes"
        CX_LDFLAGS="-L$cpl_cext_libraries"
        LIBCEXT="-lcext"
    else
        AC_MSG_RESULT([disabled])
        AC_MSG_WARN([Cext checks have been disabled! This package may not build!])
        CX_INCLUDES=""
        CX_LDFLAGS=""
        LIBCEXT=""
    fi

    AC_SUBST(CX_INCLUDES)
    AC_SUBST(CX_LDFLAGS)
    AC_SUBST(LIBCEXT)

])
