gravi_image_reconstruct recipe for the Gravity ESO pipeline
-----------------------------------------------------------

Instructions to add Yorick, Yeti and Optimpack to this recipe.

The sources of the packages are obtained from:
  http://yorick.sourceforge.net/
  http://www-obs.univ-lyon1.fr/labo/perso/eric.thiebaut/yeti.html
  http://www-obs.univ-lyon1.fr/labo/perso/eric.thiebaut/optimpack.html

The source code of Yorick must be unpackaged into the directory
  yorick
The source code of Yeti and Optimpack, which are plug-ins for Yorick,
should be placed inside the sub-directories yeti and optimpack of the
directory yorick

AUTOMAKE/AUTOCONF

To comply with ESO pipeline guideline's, all the Makefile distributed with
Yorick, Yeti and Optimpack have been removed and replaced by Makefile.am
files that we have created.

Other 3 files that we have modified in the ./yorick directory
are the following:
  - ./yorick/configure
      A line has been added at the beginning, to set up the NO_XLIB
      system flag, to build Yorick without X11 libraries.
  - ./yorick/ysite.sh
      This file has been replaced by ./yorick/ysite.sh.in
      which adds some AUTOCONF global variables to define the specific
      location where Yorick will be installed.
  - ./yorick/Makepkg
      Has been modified to substitute some AUTOMAKE standard targets
      by other targets that will avoid conflicts. The targets that
      were replaced are the following:
        install ---> install-data-local
        install-exe ---> install-exec-local
        uninstall ---> uninstall-local
        clean ---> clean-local
        distclean ---> distclean-local

One additional file has been copied to the ./yorick directory, named
"config.h" which comes from the source code of Yeti, and it should be placed
there (as described in the Yeti installation procedure) to be used in the
installation of Yeti.

UPGRADE

    If a new version of Yorick is going to be used, the yorick
    directory should be replaced with the new source code. Makefile
    files should be erased. The following shell command can be issued:
      for i in $(find -name Makefile -print); do rm $i; done

    The Makefile.am files from the ./yorick tree in the older
    gravi_image_reconstruct recipe should be copied into the new one.
    Notice that most of those files were created from the
    corresponding source Makefile, by adding some more specific
    building rules and, in some cases, incorporating the Makepkg file,
    rather than including it at run-time. (in those cases the "include
    Makepkg" command appears commented out, and the customized Makepkg
    follows). If the source Makefile changes, a comparison with the
    Makefile.am should be made to assess the necessary changes that
    must be made to Makefile.am.

    The 3 files described in the section automake/autoconf should also
    be included from the older recipe and compared to the corresponding
    files in the new Yorick source code.

Yeti and Optimpack: New versions should replace the yorick/yeti
    and yorick/optimpack trees. Makefile files should be
    removed and replaced by Makefile.am files, as explained in
    the upgrading of Yorick.

Jaime Villate <villate@fe.up.pt>. December 2011.
