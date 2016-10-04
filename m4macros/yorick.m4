# ===========================================================================
#     http://www.gnu.org/software/autoconf-archive/ax_normalize_path.html
# ===========================================================================#
# SYNOPSIS
#
#   AX_NORMALIZE_PATH(VARNAME, [REFERENCE_STRING])
#
# DESCRIPTION
#
#   Perform some cleanups on the value of $VARNAME (interpreted as a path):
AU_ALIAS([ADL_NORMALIZE_PATH], [AX_NORMALIZE_PATH])
AC_DEFUN([AX_NORMALIZE_PATH],
[case ":[$]$1:" in
# change empty paths to '.'
  ::) $1='.' ;;
# strip trailing slashes
  :*[[\\/]]:) $1=`echo "[$]$1" | sed 's,[[\\/]]*[$],,'` ;;
  :*:) ;;
esac
# squeze repeated slashes
case ifelse($2,,"[$]$1",$2) in
# if the path contains any backslashes, turn slashes into backslashes
 *\\*) $1=`echo "[$]$1" | sed 's,\(.\)[[\\/]][[\\/]]*,\1\\\\,g'` ;;
# if the path contains slashes, also turn backslashes into slashes
 *) $1=`echo "[$]$1" | sed 's,\(.\)[[\\/]][[\\/]]*,\1/,g'` ;;
esac])


# LIB_CHECK_YORICK
#------------------
# Checks for the yorick executable.
AC_DEFUN([CHECK_YORICK],
[

    AC_ARG_WITH(yorick,
                AC_HELP_STRING([--with-yorick],
                               [location where yorick binary is installed]),
                [
                    yorick_dir=$withval
                ])

    AC_PATH_PROG(YORICK_BIN, yorick, [] ,[$yorick_dir:$PATH])
    if test "x$YORICK_BIN" = "x"; then
        AC_MSG_WARN([yorick executable not found])
    fi

    AX_NORMALIZE_PATH(YORICK_BIN)

    AC_SUBST(YORICK_BIN)

    AC_DEFINE_UNQUOTED(YORICK_BIN, "$YORICK_BIN",
                       [Absolute path to the Yorick executable])


])

# CHECK_YORICK_YETI
#------------------
# Checks for the yorick yeti plugin
AC_DEFUN([CHECK_YORICK_PLUGIN],
[
    AC_REQUIRE([CHECK_YORICK])

    AC_MSG_CHECKING([$1 plugin for yorick $YORICK_BIN])

    AC_LANG_PUSH(C)
    
    AC_RUN_IFELSE([AC_LANG_PROGRAM(
        [[
        #include <unistd.h>
        #include <sys/types.h>
        #include <stdio.h>
        #include <sys/wait.h>
        #include <stdlib.h>
        ]],
        [
        int status;
        pid_t pid;
        char *const program = YORICK_BIN;
        char *param[[]] = {YORICK_BIN, "-batch", $1, (char*)0};
        switch(pid=fork())
        {
            case 0: /* Child process (Yorick) */
                if (execv(program,param) == -1 ) /* execv itself failed */
                    exit(1);
                return 0;
            case -1: /* Fork failed */
                exit(2);
            default:
                wait(&status);
                if(status) /* Error in child process */
                    exit(3);
        }
        return 0;
        ])],
        [AC_MSG_RESULT([available])],
        [AC_MSG_WARN([$1 plugin for yorick not installed])])

    AC_LANG_POP(C)

])

# CHECK_YORICK_YETI
#------------------
# Checks for the yorick yeti plugin
AC_DEFUN([CHECK_YORICK_YETI],
[
    CHECK_YORICK_PLUGIN(["yeti.i"])
])

#
# CHECK_YORICK_OPTIMPACK
#------------------
# Checks for the OptimPack plugin
AC_DEFUN([CHECK_YORICK_OPTIMPACK],
[
    CHECK_YORICK_PLUGIN(["OptimPack1.i"])
])
