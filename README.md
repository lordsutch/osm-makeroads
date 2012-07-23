osm-makeroads
=============

R code to "average" GPS traces from OpenStreetMap to derive usable GPX
tracks for roadways.

Dependencies:

* R (core)
* R packages: princurve, rgdal, maptools (all on CRAN)
* gpsbabel

To get started using Debian unstable (or presumably Ubuntu w/universe):

apt-get install r-cran-maptools gpsbabel libgdal-dev

R -e "install.packages('princurve', 'rgdal')"

makeroads.R is lightly commented.  You'll get a bunch of GPX files
named 'roadway-*.gpx', one per identified track, along with Rplots.pdf
containing a visual representation of each estimated track along with
the tracks that were used to construct it.
