# GRAVI_SET_VERSION_INFO(VERSION, [CURRENT], [REVISION], [AGE])
#----------------------------------------------------------------
# Setup various version information, especially the libtool versioning
AC_DEFUN([GRAVI_SET_VERSION_INFO],
[
    gravi_version=`echo "$1" | sed -e 's/[[a-z,A-Z]].*$//'`

    gravi_major_version=`echo "$gravi_version" | \
        sed 's/\([[0-9]]*\).\(.*\)/\1/'`
    gravi_minor_version=`echo "$gravi_version" | \
        sed 's/\([[0-9]]*\).\([[0-9]]*\)\(.*\)/\2/'`
    gravi_micro_version=`echo "$gravi_version" | \
        sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`

    if test -z "$gravi_major_version"; then gravi_major_version=0
    fi

    if test -z "$gravi_minor_version"; then gravi_minor_version=0
    fi

    if test -z "$gravi_micro_version"; then gravi_micro_version=0
    fi

    GRAVI_VERSION="$gravi_version"
    GRAVI_MAJOR_VERSION=$gravi_major_version
    GRAVI_MINOR_VERSION=$gravi_minor_version
    GRAVI_MICRO_VERSION=$gravi_micro_version

    if test -z "$4"; then GRAVI_INTERFACE_AGE=0
    else GRAVI_INTERFACE_AGE="$4"
    fi

    GRAVI_BINARY_AGE=`expr 100 '*' $GRAVI_MINOR_VERSION + $GRAVI_MICRO_VERSION`
    GRAVI_BINARY_VERSION=`expr 10000 '*' $GRAVI_MAJOR_VERSION + \
                          $GRAVI_BINARY_AGE`

    AC_SUBST(GRAVI_VERSION)
    AC_SUBST(GRAVI_MAJOR_VERSION)
    AC_SUBST(GRAVI_MINOR_VERSION)
    AC_SUBST(GRAVI_MICRO_VERSION)
    AC_SUBST(GRAVI_INTERFACE_AGE)
    AC_SUBST(GRAVI_BINARY_VERSION)
    AC_SUBST(GRAVI_BINARY_AGE)

    AC_DEFINE_UNQUOTED(GRAVI_MAJOR_VERSION, $GRAVI_MAJOR_VERSION,
                       [GRAVI major version number])
    AC_DEFINE_UNQUOTED(GRAVI_MINOR_VERSION, $GRAVI_MINOR_VERSION,
                       [GRAVI minor version number])
    AC_DEFINE_UNQUOTED(GRAVI_MICRO_VERSION, $GRAVI_MICRO_VERSION,
                       [GRAVI micro version number])
    AC_DEFINE_UNQUOTED(GRAVI_INTERFACE_AGE, $GRAVI_INTERFACE_AGE,
                       [GRAVI interface age])
    AC_DEFINE_UNQUOTED(GRAVI_BINARY_VERSION, $GRAVI_BINARY_VERSION,
                       [GRAVI binary version number])
    AC_DEFINE_UNQUOTED(GRAVI_BINARY_AGE, $GRAVI_BINARY_AGE,
                       [GRAVI binary age])

    ESO_SET_LIBRARY_VERSION([$2], [$3], [$4])
])


# GRAVI_SET_PATHS
#------------------
# Define auxiliary directories of the installed directory tree.
AC_DEFUN([GRAVI_SET_PATHS],
[

    if test -z "$plugindir"; then
        plugindir='${libdir}/esopipes-plugins/${PACKAGE}-${VERSION}'
    fi

    if test -z "$privatelibdir"; then
        privatelibdir='${libdir}/${PACKAGE}-${VERSION}'
    fi

    if test -z "$apidocdir"; then
        apidocdir='${datadir}/doc/esopipes/${PACKAGE}-${VERSION}/html'
    fi

    if test -z "$pipedocsdir"; then
        pipedocsdir='${datadir}/doc/esopipes/${PACKAGE}-${VERSION}'
    fi

    if test -z "$configdir"; then
       configdir='${datadir}/${PACKAGE}/config'
    fi

    if test -z "$wkfextradir"; then
        wkfextradir='${datadir}/esopipes/${PACKAGE}-${VERSION}/reflex'
    fi

    if test -z "$wkfcopydir"; then
        wkfcopydir='${datadir}/reflex/workflows/${PACKAGE}-${VERSION}'
    fi

    AC_SUBST(plugindir)
    AC_SUBST(privatelibdir)
    AC_SUBST(apidocdir)
    AC_SUBST(pipedocsdir)
    AC_SUBST(configdir)
    AC_SUBST(wkfextradir)
    AC_SUBST(wkfcopydir)


    # Define a preprocesor symbol for the plugin search paths

    AC_DEFINE_UNQUOTED(GRAVI_PLUGIN_DIR, "${PACKAGE}/plugins",
                       [Plugin directory tree prefix])

    eval plugin_dir="$plugindir"
    plugin_path=`eval echo $plugin_dir | \
                sed -e "s/\/${PACKAGE}-${VERSION}.*$//"`

    AC_DEFINE_UNQUOTED(GRAVI_PLUGIN_PATH, "$plugin_path",
                       [Absolute path to the plugin directory tree])

])


# GRAVI_CREATE_SYMBOLS
#-----------------------
# Define include and library related makefile symbols
AC_DEFUN([GRAVI_CREATE_SYMBOLS],
[

    # Symbols for package include file and library search paths

    GRAVI_INCLUDES='-I$(top_srcdir)/gravi -I$(top_srcdir)/irplib'
    GRAVI_LDFLAGS='-L$(top_builddir)/gravi'

    # Library aliases

    LIBGRAVI='$(top_builddir)/gravi/libgravi.la'
    LIBIRPLIB='$(top_builddir)/irplib/libirplib.la'

    # Substitute the defined symbols

    AC_SUBST(GRAVI_INCLUDES)
    AC_SUBST(GRAVI_LDFLAGS)

    AC_SUBST(LIBGRAVI)
    AC_SUBST(LIBIRPLIB)

    # Check for CPL and user defined libraries
    AC_REQUIRE([CPL_CHECK_LIBS])
    AC_REQUIRE([ESO_CHECK_EXTRA_LIBS])

    all_includes='$(GRAVI_INCLUDES) $(CPL_INCLUDES) $(EXTRA_INCLUDES)'
    all_ldflags='$(GRAVI_LDFLAGS) $(CPL_LDFLAGS) $(EXTRA_LDFLAGS)'

    AC_SUBST(all_includes)
    AC_SUBST(all_ldflags)
])


# MIRA_SET_PATHS
#------------------
AC_DEFUN([MIRA_SET_PATHS],
[

    if test -z "$miradir"; then
        miradir='${libdir}/${PACKAGE}-${VERSION}/mira'
    fi

    AC_SUBST(miradir)

    # Define preprocesor symbols for the yorick and mira programs.

    eval mira_dir="$miradir"

    eval mira='${prefix}/lib/${PACKAGE}-${VERSION}/mira/mira-script.i'

    AC_DEFINE_UNQUOTED(YORICK_MIRA, "$mira",
                       [Absolute path to the Mira script])

])
