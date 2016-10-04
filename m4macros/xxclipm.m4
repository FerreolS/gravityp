#
# XXCLIPM_CREATE_SYMBOLS(libdir=[],libname)
#------------------------------------------
# Sets the Makefile symbols for the <XX>CLIPM libraries. If libdir
# is provided, then the symbols are setup for building CPL; if libdir is
# omitted (default), then the symbols are set for using the libraries
# for external package development.
AC_DEFUN([XXCLIPM_CREATE_SYMBOLS],
[
    if test -z "$1"; then
    	LIBXXCLIPM="$LIBXXCLIPM -l$2"
    else
    	LIBXXCLIPM="$LIBXXCLIPM $1/$2"
    fi
    AC_SUBST(LIBXXCLIPM)
])

# _XXCLIPM_FIND_FILE_IN_PATH_VARS(path_variable_names,
#                                 search_file_name,
#                                 return_file_location,
#                                 return_path_variable_name)
#-----------------------------------------------------------
# Searches a file in multiple path variables
AC_DEFUN([_XXCLIPM_FIND_FILE_IN_PATH_VARS],
[
    for temp_pathvar in $1; do
        ## remember double square brackets in m4 macros
        temp_path=`echo ${!temp_pathvar} | sed "s|[[:;]]| |g"`
        if test -n "$temp_path"; then
            ESO_FIND_FILE($2, \
                          $temp_path, \
                          $3)
            if test "$3" != "no"; then
                $4=$temp_pathvar
                break
            fi
        fi
    done
])

# XXCLIPM_CHECK_LIB(libname)
#----------------------------
# Checks for the <xx>clipm libraries and header files in:
# - standard locations, and
# - $LD_LIBRARY_PATH, $SHLIB_PATH (HPUX), $LIBPATH (AIX), and
#   $DYLD_LIBRARY_PATH (Mac)(replacing "/lib" with "/include" for the headers).
# Sets the symbols:
# - XXCLIPM_INCLUDES
# - XXCLIPM_LDFLAGS
# - LIBXXCLIPM
# Can be called more than once, in this case the above symbols will contain
# lists of the required elements.
AC_DEFUN([XXCLIPM_CHECK_LIB],
[

    libname=$1
    LIBNAME=`echo $1 | \
             sed 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`

    AC_MSG_CHECKING([for $LIBNAME])

    xxclipm_check_xxclipm_header="$libname.h"
    xxclipm_check_xxclipm_lib="lib$libname.la"

    xxclipm_includes=""
    xxclipm_libraries=""
    xxclipm_with_xxclipm_includes=""
    xxclipm_with_xxclipm_libs=""
    
    # additional library search path variable names
    sharedlib_path_vars="LD_LIBRARY_PATH SHLIB_PATH LIBPATH DYLD_LIBRARY_PATH"
    found_path_var=""

    # the AC_ARG_WITH macro does not accept a variable as first argument,
    # therefore do a workaround:
    # - copy the value of with_$libname_suffix into with_xxxxclipm_suffix
    # - call AC_ARG_WITH with xxxxclipm-suffix
    # This is allowed according to chapter 12.1 of the autoconf documentation
    # ("Site Configuration" - "Working With External Software") which says:
    # "The option's argument is available to the shell commands action-if-given
    #  in the shell variable withval, which is actually just the value of the
    #  shell variable with_package, with any `-' characters changed into `_'.
    #  You may use that variable instead, if you wish."

    # _xxclipm_get_option(prefix, suffix)
    _xxclipm_get_option() {
	    echo "unset ${1}xxxxclipm${2};"
	    echo "varname=${1}\${libname}${2};"
	    echo "test -n \"\${!varname}\" && ${1}xxxxclipm${2}=\${!varname};"
    }
    
    eval `_xxclipm_get_option "with_" ""`
    AC_ARG_WITH(xxxxclipm,
                AC_HELP_STRING([--with-<xx>clipm],
                               [location where <XX>CLIPM is installed]),
                [
                    xxclipm_with_xxclipm_includes=$withval/include
                    xxclipm_with_xxclipm_libs=$withval/lib
                ])

    eval `_xxclipm_get_option "with_" "_includes"`
    AC_ARG_WITH(xxxxclipm-includes,
                AC_HELP_STRING([--with-<xx>clipm-includes],
                               [location of the <XX>CLIPM header files]),
                xxclipm_with_xxclipm_includes=$withval)

    eval `_xxclipm_get_option "with_" "_libs"`
    AC_ARG_WITH(xxxxclipm-libs,
                AC_HELP_STRING([--with-<xx>clipm-libs],
                               [location of the <XX>CLIPM library]),
                xxclipm_with_xxclipm_libs=$withval)

    eval `_xxclipm_get_option "enable_" "_test"`
    AC_ARG_ENABLE(xxxxclipm-test,
                  AC_HELP_STRING([--disable-<xx>clipm-test],
                                 [disables checks for the <XX>CLIPM \
                                  library and headers]),
                  xxclipm_enable_xxclipm_test=$enableval,
                  xxclipm_enable_xxclipm_test=yes)

    if test "x$xxclipm_enable_xxclipm_test" = xyes; then

        # Check for the <XX>CLIPM includes
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if test -z "$xxclipm_with_xxclipm_includes"; then
            xxclipm_with_xxclipm_includes="/opt/$libname/include \
                         /usr/local/include \
                         /usr/local/$libname/include \
                         /usr/local/include/$libname \
                         /usr/include/$libname \
                         /usr/include"
        fi

        ESO_FIND_FILE($xxclipm_check_xxclipm_header, \
                      $xxclipm_with_xxclipm_includes, \
                      xxclipm_includes)

        # look in the directories of sharedlib_path_vars,
        # call in function to harmlessly modify the paths
        _XXCLIPM_LIBPATH_INCLUDES() {
            for pathvar in $sharedlib_path_vars; do
                eval $pathvar=`echo ${!pathvar} | sed 's|/lib|/include|g'`
            done
            
            _XXCLIPM_FIND_FILE_IN_PATH_VARS($sharedlib_path_vars, \
                                            $xxclipm_check_xxclipm_header, \
                                            xxclipm_includes, \
                                            found_path_var)
            echo "xxclipm_includes=$xxclipm_includes;"
            echo "found_path_var=$found_path_var;"
        }
        if test "$xxclipm_includes" = "no"; then
            eval `_XXCLIPM_LIBPATH_INCLUDES`
        fi

        # Check for the <XX>CLIPM libraries
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if test -z "$xxclipm_with_xxclipm_libs"; then
            xxclipm_with_xxclipm_libs="/opt/$libname/lib \
                         /usr/local/lib \
                         /usr/local/$libname/lib \
                         /usr/lib"
        fi

        ESO_FIND_FILE($xxclipm_check_xxclipm_lib, \
                      $xxclipm_with_xxclipm_libs, \
                      xxclipm_libraries)

        # look in the directories of sharedlib_path_vars
        if test x"$xxclipm_libraries" = xno; then
	        _XXCLIPM_FIND_FILE_IN_PATH_VARS($sharedlib_path_vars, \
	                                        $xxclipm_check_xxclipm_lib, \
	                                        xxclipm_libraries, \
	                                        found_path_var)
        fi

        # Summarise
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if test x"$xxclipm_includes" = xno || test x"$xxclipm_libraries" = xno; then
            xxclipm_notfound=""

            if test x"$xxclipm_includes" = xno; then
                if test x"$xxclipm_libraries" = xno; then
                    xxclipm_notfound="(headers and libraries)"
                else            
                    xxclipm_notfound="(headers)"
                fi
            else
                xxclipm_notfound="(libraries)"
            fi

            AC_MSG_RESULT([failed])
            AC_MSG_ERROR([$LIBNAME $xxclipm_notfound was not found on your system. Please check!])
        else
            msg="libraries: $xxclipm_libraries, headers: $xxclipm_includes"
            if test -n "$found_path_var"; then
                msg="$msg (found in \$$found_path_var)"
            fi
            AC_MSG_RESULT([$msg])
        fi

        # Set up the symbols

        if test -z "`echo "$XXCLIPM_INCLUDES " | \
                     grep "$xxclipm_includes "`"; then
            XXCLIPM_INCLUDES="$XXCLIPM_INCLUDES -I$xxclipm_includes"
        fi
        if test -z "`echo "$XXCLIPM_LDFLAGS " | \
                     grep "$xxclipm_libraries "`"; then
            XXCLIPM_LDFLAGS="$XXCLIPM_LDFLAGS -L$xxclipm_libraries"
        fi
        
        XXCLIPM_CREATE_SYMBOLS([],[$libname])
    else
        AC_MSG_RESULT([disabled])
        AC_MSG_WARN([$LIBNAME checks have been disabled! This package may not build!])
    fi

    AC_SUBST(XXCLIPM_INCLUDES)
    AC_SUBST(XXCLIPM_LDFLAGS)
    AC_SUBST(LIBXXCLIPM)

])
