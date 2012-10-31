## This is R code - http://www.r-project.org/

## Attempt at using principal curves to produce sensible routes for
## roadways using GPS traces from OpenStreetMap.
##
## Algorithm based on Chris Brunsdon, 2007, "Path Estimation from GPS
## Tracks," Proceedings of the 9th International Conference on
## GeoComputation, National Centre for Geocomputation, National
## University of Ireland, Maynooth, Eire.
##
## http://www.geocomputation.org/2007/1B-Algorithms_and_Architecture1/1B2.pdf

## (C) Copyright 2012 Chris Lawrence <lordsutch@gmail.com>
## You may freely use, redistribute, and modify this software under the
## terms of the GNU General Public License, version 2 or later.

## Bring in the generic code.  Some of this code should be moved there
## eventually.
source('makeroads.R')

savefile <- "US41-Hampton-Barnesville.Rdata"

## Set the bounding box here.
left <- -84.311
top <- 33.380
right <- -84.137
bottom <- 33.023

## Get GPS points in the area described by this bounding box
tracks <- getOSMtracks(left, top, right, bottom)

## To load GPX file contents instead...
## tracks <- getOSMtracksFiles("filename1.gpx", "filename2.gpx", etc...)

## What angle to split tracks at (degrees) - splits merged tracks
splitangle <- 90

## How close the angles need to be for ways to be "similar"
similarangle <- 90

## Max distance similar tracks can be apart (m) somewhere
## Realistically the same road usually gets within 1 meter at least once...
maxtrackdist <- 10

## Max distance a track can ever be away from the bundle before we split it
## in meters
maxseparation <- 40

## How much an opposite-angle track can be off before we reject it
## should be much smaller than above, since we need to account for medians
## in North America, 20m (60 ft) is reasonable; if your country typically has
## narrow medians, needs to be smaller.
maxoppositeseparation <- 20

## Interpolate distance (m) - ensure tracks have a point this often
## This copes with the "outlier points drag the track" issue with thinned GPX
## traces.
interpolate <- TRUE
interpolatedist <- 20

## How close to fit
delta <- 1/1800    ## Ensure curves are reconstructed within 2 deg second
maxit <- 50

debug <- TRUE

## If points are more than 250m apart, split the line
## Deals with overly-thinned traces/waypoint uploads
splitdistance <- 250

## Debugging
##tracks <- getOSMtracksFiles(list.files(pattern='file.*[.]gpx'))
##tracks <- getOSMtracksFiles(list.files(pattern='I69.*[.]gpx'))
length(tracks)

## Distance split
newtracks <- splitTracks(tracks)
length(newtracks)

## Split tracks by direction threshold
newtracks2 <- splitTracks2(newtracks)
length(newtracks2)

## Sort tracks by length
newtracks3 <- sortTracks(newtracks2)

## Simplify tracks - improves consolidateTracks performance immensely
stracks <- simplifyTracks(newtracks3)

## Find related tracks
ret <- consolidateTracks(stracks)

tracklist <- ret$tracklist
newtracks4 <- ret$tracks
trackforpoints <- ret$trackforpoints
tpoints <- ret$points

## Run a second pass
## p2tracks <- sortTracks(newtracks4)
## ret <- consolidateTracks(p2tracks)
## tracklist <- ret$tracklist
## newtracks4 <- ret$tracks

## Interpolate extra points to improve fit algorithm performance
newtracks5 <- interpolateTracks(newtracks4)

## Save the data for any future runs(?)
save(tracks, newtracks, newtracks2, newtracks3, newtracks4, newtracks5,
     tracklist, points, trackforpoints, file=savefile)

##load(savefile)

## Make the best fit tracks and save them as GPX files.
pdf(onefile=T, height=8, width=6)
for(tcount in seq_along(tracklist)) {
  group <- tracklist[[tcount]]
  mergedtracks <- NULL
  plot(newtracks5[[1]], xlim=c(left, right), ylim=c(bottom, top), asp=0.75,
       type='l', col='gray50')
  color <- 3
  for(t in group) {
    mergedtracks <- rbind(mergedtracks, newtracks5[[t]])
    lines(newtracks5[[t]], pch='.', col=color, lwd=0.2)
    if(nrow(newtracks5[[t]]) > 2) {
      text(newtracks5[[t]][1,], labels=t, col=color, cex=.5)
    }
    color <- color+1
  }
  thesepoints <- do.call(rbind, tpoints[trackforpoints == tcount])
  points(thesepoints, pch='x', col=2, cex=.3)
  mergedtracks <- rbind(mergedtracks, thesepoints)

  if(nrow(mergedtracks) < 50) {
    show(paste('Skipping small group', tcount))
    next
  }
  
  ## Fit the curve and render it onto the plot
  f <- fitpcurve(mergedtracks)
  lines(f$s[f$tag,], asp=0.8, lty=3, col=2)

  fname <- paste('roadway-', tcount, sep='')
  csvname <- tempfile(fileext='.csv')
  gpxname <- paste(fname, 'gpx', sep='.')

  points <- f$s[f$tag,]
  ## Export the curve as a series of points in a CSV file, dropping
  ## excess precision; OSM stores 6 decimal places
  write.csv(round(points, 6), file=csvname)

  ## Then run GPSBabel to convert the CSV to GPX
  ## Could probably do this with an R package but this is simple enough for now
  gpsbabel.out(csvname, gpxname)

  if(!debug) unlink(csvname)
}
dev.off()
