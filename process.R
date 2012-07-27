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

source('makeroads.R')

## Set the bounding box here.
left <- -81.91
right <- -81.86
top <- 34.50
bottom <- 34.61

## What angle to split tracks at (degrees) - splits merged tracks
splitangle <- 60

## How close the angles need to be for ways to be "similar"
similarangle <- 30

## Max distance similar tracks can be apart (m) somewhere
## Realistically the same road usually gets within 1 meter at least once...
maxtrackdist <- 10

## Max distance a track can ever be away from the bundle to be distinct (m)
maxseparation <- 100

## How much an opposite-angle track can be off before we reject it
## should be much smaller than above, since we need to account for medians
## in North America, 15m (45 ft) is reasonable; if your country typically has
## narrow medians, needs to be smaller.
maxoppositeseparation <- 15

## Interpolate distance (m) - ensure tracks have a point this often
## This copes with the "outlier points drag the track" issue with thinned GPX
## traces.  Downside: dist2Line gets a lot slower.
interpolate <- TRUE
interpolatedist <- 100

## How close to fit
delta <- 1/1800    ## Ensure curves are reconstructed within 2 deg second
maxit <- 50

debug <- TRUE

## How much distance to split tracks by (km)
## This deals with points w/o track info that tend to wrap around
## at the edge of the downloaded area.

## We'll say 25% of the diagonal distance is reasonable(?),
## with a minimum threshold of 1 km (1000m)
splitdistance <- max(distHaversine(c(top, left), c(bottom, right))*.25, 1000)

## Get GPS points in the area described by this bounding box
tracks <- getOSMtracks(left, bottom, right, top)

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

## Find related tracks
tracklist <- consolidateTracks(newtracks3)

## Interpolate extra points to improve fit algorithm performance
newtracks4 <- interpolateTracks(newtracks3)

## Save the data for any future runs(?)
save(tracks, newtracks, newtracks2, newtracks3, newtracks4, tracklist,
     file='I26.Rdata')

## Make the best fit tracks and save them as GPX files.
tcount <- 0
pdf(onefile=T, height=8, width=6)
for(group in tracklist) {
  tcount <- tcount + 1
  mergedtracks <- NULL
  plot(newtracks4[[1]], xlim=c(left, right), ylim=c(bottom, top), asp=0.75,
       type='l', col='gray50')
  color <- 3
  for(t in group) {
    mergedtracks <- rbind(mergedtracks, newtracks4[[t]])
    lines(newtracks4[[t]], pch='.', col=color, lwd=0.2)
    text(newtracks4[[t]][1,], labels=t, col=color)
    color <- color+1
  }

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
  ## Export the curve as a series of points in a CSV file
  write.csv(f$s[f$tag,], file=csvname)

  ## Then run GPSBabel to convert the CSV to GPX
  ## Could probably do this with an R package but this is simple enough for now
  gpsbabel.out(csvname, gpxname)

  if(!debug) unlink(csvname)
}
dev.off()
