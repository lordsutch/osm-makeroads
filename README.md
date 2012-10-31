osm-makeroads
=============

R code to "average" GPS traces from OpenStreetMap to derive usable GPX
tracks for roadways.

Dependencies:

* R (core)
* R packages: princurve, rgdal, maptools, geosphere (all on CRAN)
* gpsbabel
* Also: libgdal and gfortran are required by some of the dependencies above.

To get started using Debian unstable (or presumably Ubuntu w/universe):

apt-get install r-cran-maptools gpsbabel libgdal-dev gfortran libgdal-dev

R -e "install.packages('princurve', 'rgdal', 'geosphere', dependencies=TRUE)"

makeroads.R and process.R are lightly commented.  Place them in the
same directory (and cd there), edit process.R to your liking, then use:

R -f process.R

You'll get a bunch of GPX files named 'roadway-*.gpx', one per
identified track, along with Rplots.pdf containing a visual
representation of each estimated track along with the tracks that were
used to construct it.

To-Do
=====

* Optimize search algorithm to avoid comparing intensively comparing
  data from across the map (e.g. first pass to determine bbox of each
  track before we bother to intensively match).
* Take advantage of parallelism (investigate R multicore?).
* Memory usage improvements.

Thanks
======

* A few problems with the instructions were pointed out by malenki:
  http://www.openstreetmap.org/user/malenki/diary/17976
