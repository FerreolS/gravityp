gravi_image_reconstruct recipe for the Gravity ESO pipeline
-----------------------------------------------------------

- The pipeline template used was the file iiinstrument.tar.gz provided
  by ESO on November 2011 (the file name does not have a version
  number, but it differs from previous templates by its size in bytes:
  88206)

- All occurrences of iiinstrument were replaced by gravi and all
  occurrences of rrrecipe were replaced by gravi_image_reconstruct
  (respecting the combinations of upper/lower case in the original).

MODIFIED FILES

configure.ac
  * updated version name and e-mail address in the AC_INIT macro
  * execute the macro YORICK_SET_PATHS
  * the following files were added to the AC_CONFIG_FILES macro:
                recipes/Makefile
                yorick/Makefile
                yorick/yorick/Makefile
                yorick/optimpack/Makefile
                yorick/yeti/Makefile
                yorick/play/Makefile
                yorick/play/any/Makefile
                yorick/play/hacks/Makefile
                yorick/play/unix/Makefile
                yorick/play/x11/Makefile
                yorick/play/win/Makefile
                yorick/gist/Makefile
                yorick/matrix/Makefile
                yorick/fft/Makefile
                yorick/regexp/Makefile
                yorick/drat/Makefile
                yorick/hex/Makefile
                yorick/mpy/Makefile
                yorick/doc/Makefile
                yorick/ysite.sh
                mira/Makefile

acinclude.m4
  Added a YORICK_SET_PATHS macro at the end

Makefile.am
  * Added
      yorick mira
    to the variable COMPILE_FIRST

      dist_plugin_DATA = mira-script.i mira.i fft_utils.i fits.i fmin.i img.i \
                linop.i oifits.i options.i plot.i rgl.i utils.i

recipes/gravi_image_reconstruct.c
  Added the C functions specific to the recipe.

REMOVED FILES
  README

NEW FILES
  README_image_reconstruction.txt
  README_yorick_update.txt
  README_mira_update.txt

Jaime Villate <villate@fe.up.pt>. December 2011.
