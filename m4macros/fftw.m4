#
# FFTW_CREATE_SYMBOLS(build=[])
#-----------------------------
# Sets the Makefile symbols for the FFTW library.
AC_DEFUN([FFTW_CREATE_SYMBOLS],
[
   LIBFFTW='-lsfftw'
   AC_SUBST(LIBFFTW)
])


# FFTW_CHECK_LIBS
#---------------
# Checks for the FFTW libraries and header files.
AC_DEFUN([FFTW_CHECK_LIBS],
[

    AC_MSG_CHECKING([for FFTW])

    fftw_check_fftw_header="sfftw.h"
    fftw_check_fftw_lib="libsfftw.a"

    fftw_includes=""
    fftw_libraries=""

    AC_ARG_WITH(fftw,
                AC_HELP_STRING([--with-fftw],
                               [location where FFTW is installed]),
                [
                    fftw_with_fftw_includes=$withval/include
                    fftw_with_fftw_libs=$withval/lib
                ])

    AC_ARG_WITH(fftw-includes,
                AC_HELP_STRING([--with-fftw-includes],
                               [location of the FFTW header files]),
                fftw_with_fftw_includes=$withval)

    AC_ARG_WITH(fftw-libs,
                AC_HELP_STRING([--with-fftw-libs],
                               [location of the FFTW library]),
                fftw_with_fftw_libs=$withval)

    AC_ARG_ENABLE(fftw-test,
                  AC_HELP_STRING([--disable-fftw-test],
                                 [disables checks for the FFTW library and headers]),
                  fftw_enable_fftw_test=$enableval,
                  fftw_enable_fftw_test=yes)


    if test "x$fftw_enable_fftw_test" = xyes; then

        # Check for the FFTW includes

        if test -z "$fftw_with_fftw_includes"; then
            fftw_incdirs="/opt/fftw/include \
                          /usr/local/include \
                          /usr/local/fftw/include \
                          /usr/local/include/fftw \
                          /usr/include/fftw \
                          /usr/include"

            test -n "$FFTWDIR" && fftw_incdirs="$FFTWDIR/include $fftw_incdirs"
        else
            fftw_incdirs="$fftw_with_fftw_includes"
        fi

        ESO_FIND_FILE($fftw_check_fftw_header, $fftw_incdirs, fftw_includes)


        # Check for the FFTW libraries

        if test -z "$fftw_with_fftw_libs"; then
            fftw_libdirs="/opt/fftw/lib \
                         /usr/local/lib \
                         /usr/local/fftw/lib \
                         /usr/lib"

            test -n "$FFTWDIR" && fftw_libdirs="$FFTWDIR/lib $fftw_libdirs"
        else
            fftw_libdirs="$fftw_with_fftw_libs"
        fi

        ESO_FIND_FILE($fftw_check_fftw_lib, $fftw_libdirs, fftw_libraries)


        if test x"$fftw_includes" = xno || test x"$fftw_libraries" = xno; then
            fftw_notfound=""

            if test x"$fftw_includes" = xno; then
                if test x"$fftw_libraries" = xno; then
                    fftw_notfound="(headers and libraries)"
                else            
                    fftw_notfound="(headers)"
                fi
            else
                fftw_notfound="(libraries)"
            fi

            AC_MSG_ERROR([FFTW $fftw_notfound was not found on your system. Please check!])
        else
            AC_MSG_RESULT([libraries $fftw_libraries, headers $fftw_includes])
        fi

        # Set up the symbols

	AC_DEFINE_UNQUOTED(HAVE_FFTW, 1, [Define to 1 if you have FFTW v. 2.1.5])
        FFTW_INCLUDES="-I$fftw_includes"
        FFTW_LDFLAGS="-L$fftw_libraries"

        FFTW_CREATE_SYMBOLS
    else
        AC_MSG_RESULT([disabled])
        AC_MSG_WARN([FFTW checks have been disabled! This package may not build!])
        FFTW_INCLUDES=""
        FFTW_LDFLAGS=""
        LIBFFTW=""
    fi

    AC_SUBST(FFTW_INCLUDES)
    AC_SUBST(FFTW_LDFLAGS)
    AC_SUBST(LIBFFTW)

])

