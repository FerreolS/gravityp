gravi_image_reconstruct recipe for the Gravity ESO pipeline
-----------------------------------------------------------

Instructions to add Mira to this recipe.

Mira is distributed as a collection of text-file scripts written in the
high-level language of Yorick. They can be obtained from
  http://www-obs.univ-lyon1.fr/labo/perso/eric.thiebaut

Those scripts, identified by a name with an ".i" extension, have been
placed inside the ./mira directory.

That directory includes also a simple script "mira-script.i" that was
provided directly by Mira's author for this recipe's specific need of
running Mira in batch mode (as far as we know that scripts has not
been publicly released with the Mira package).

UPGRADE

In order to upgrade to a newer version, all the "*.i" scripts in the
mira/ directory should be replaced with the new ones and the
mira/Makefile.am file from the older gravity pipeline
should be copied into the new pipeline tree.

Jaime Villate <villate@fe.up.pt>. December 2011.
