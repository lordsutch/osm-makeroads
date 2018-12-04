osm-makeroads
=============

R code to "average" GPS traces from OpenStreetMap to derive usable GPX
tracks for roadways.

Dependencies:

* R (core)
* R packages: princurve, rgdal, maptools, geosphere, rgeos (all on CRAN)
* gpsbabel
* Also: libgdal and gfortran are required by some of the dependencies above.

To get started using Debian unstable (or presumably Ubuntu w/universe):

`apt-get install libgeos-dev r-cran-maptools gpsbabel libgdal-dev gfortran libproj-dev`

`R -e "install.packages(c('princurve', 'geosphere', 'rgeos', 'rgdal'), dependencies=TRUE)"`

makeroads.R and process.R are lightly commented.  Place them in the
same directory (and cd there), edit process.R to your liking, then use:

`R -f process.R`

You'll get a bunch of GPX files named 'roadway-*.gpx', one per
identified track, along with Rplots.pdf containing a visual
representation of each estimated track along with the tracks that were
used to construct it.

Recent Changes
==============

* This version uses rgeos to ridiculously improve the time taken to
  match the tracks.  (Data matching now takes an order of minutes
  versus hours in previous versions.)

To-Do
=====

* Use local principal curves (R LPCM package).
* Take advantage of parallelism (investigate R multicore?).
* Memory usage improvements.

Thanks
======

* A few problems with the instructions were pointed out by malenki:
  http://www.openstreetmap.org/user/malenki/diary/17976
* Instructions improved thanks to Valdas Kondrotas.
