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

## Filename to process
loadfile <- "logatec_region.gpx"

fileprefix <- "logatec"

savefile <- paste(fileprefix, "Rdata", sep='.')

## What angle to split tracks at (degrees) - splits merged tracks
splitangle <- 135

## How close the angles need to be for ways to be "similar" enough to compare
similarangle <- 60

## Max distance similar tracks can be apart (m) somewhere
## Realistically the same road usually gets within 1 meter at least once...
maxtrackdist <- 10

## Max distance a track can ever be away from the bundle before we split it
## in meters
maxseparation <- 20

## How much an opposite-angle track can be off before we reject it
## should be much smaller than above, since we need to account for medians
## in North America, 20m (60 ft) is reasonable; if your country typically has
## narrow medians, needs to be smaller.
maxoppositeseparation <- 15

## Interpolate distance (m) - ensure tracks have a point this often
## This copes with the "outlier points drag the track" issue with thinned GPX
## traces.
interpolate <- TRUE
interpolatedist <- 250

## If points are more than 30m apart, split the line
## Deals with overly-thinned traces/waypoint uploads
splitdistance <- 30

# Minimum track length before we turn it into a point
mintracklen <- 10

## How close to fit
delta <- 1/50000
maxit <- 50

debug <- TRUE

## Below here, there be dragons.... actual code follows

## Bring in the generic code.  Some of this code should be moved there
## eventually.
source('makeroads.R')

ret <- getOSMtracksFiles(loadfile)
tracks <- ret$tracks
points <- ret$points

bounds <- bbox(tracks)
minlat <- bounds['y','min']
maxlat <- bounds['y','max']
minlon <- bounds['x','min']
maxlon <- bounds['x','max']

UTMzone <- trunc((180+mean(minlon,maxlon))/6)+1
projection <- paste("+proj=utm +datum=WGS84 +zone=", UTMzone, sep="")
if(maxlon < 0)
  projection <- paste(projection, '+south')

# Project tracks to local UTM zone
ptracks <- spTransform(tracks, CRS(projection))
ppoints <- spTransform(points, CRS(projection))

## Put our original points in the point field
ppoints <- rbind(ppoints, SpatialPoints(coordinatesSL(ptracks),
                                        proj4string=CRS(projection)))
## Simplify the tracks using Douglas-Peucker to reduce overnoding
ptracks <- gSimplify(ptracks, 3, topologyPreserve=TRUE)

ret <- splitTracksAngle(ptracks, splitangle)
stracks <- ret$tracks
if(!is.null(ret$points)) {
  spoints <- rbind(ret$points, ppoints)
} else {
  spoints <- ppoints
}

## ftracks <- flattenSpatialLines(stracks)

ftracks <- sortTracks(flattenSpatialLines(stracks))

## Find related tracks
##system.time(ret <- consolidateTracks(ftracks[-11], spoints))
system.time(ret <- consolidateTracks(ftracks, spoints))

tracklist <- ret$tracklist
xtracks <- ret$tracks
trackforpoints <- ret$trackforpoints
xpoints <- ret$points

## Save the data for any future runs(?)
save(tracks, points, xtracks, xpoints, tracklist, trackforpoints,
     file=savefile)

##load(savefile)

## Make the best fit tracks and save them as GPX files.
pbounds <- bbox(xtracks)
pdf(onefile=T, height=8, width=6)
for(tcount in seq_along(tracklist)) {
  group <- tracklist[[tcount]]
  plot(xtracks, col='gray50', xlim=pbounds['x',], ylim=pbounds['y',])
  color <- 3
  for(t in group) {
    lines(xtracks[t], pch='.', col=color, lwd=0.2)
    text(coordinatesSL(xtracks[t])[1,], labels=t, col=color, cex=.5)
    color <- color+1
  }
  thesepoints <- xpoints[trackforpoints == tcount]
  ##points(thesepoints, pch='x', col=2, cex=.1)

  if(length(thesepoints) < 50) {
    show(paste('Skipping small group', tcount))
    next
  }

  ## Convert to coordinates
  ptracks <- coordinates(thesepoints)
  
  ## Fit the curve and render it onto the plot
  f <- fitpcurve(ptracks, projection)

  points <- f$s[f$tag,]

  lines(points, lty=3, col=2)

  points <- project(points, projection, inv=TRUE)
  
  ## Export the curve as a series of points in a CSV file, dropping
  ## excess precision; OSM stores 6 decimal places
  points <- data.frame(lon=points[,1], lat=points[,2])
  
  fname <- paste(fileprefix, '-roadway-', tcount, sep='')
  csvname <- tempfile(fileext='.csv')
  gpxname <- paste(fname, 'gpx', sep='.')

  write.csv(round(points, 6), file=csvname)

  ## Then run GPSBabel to convert the CSV to GPX
  ## Could probably do this with an R package but this is simple enough for now
  gpsbabel.out(csvname, gpxname)
  if(!debug) unlink(csvname)

  ## Also write point array
  fname <- paste(fileprefix, '-rawpoints-', tcount, sep='')
  csvname <- tempfile(fileext='.csv')
  gpxname <- paste(fname, 'gpx', sep='.')

  write.csv(coordinates(spTransform(thesepoints, CRS("+proj=latlong"))),
            file=csvname)

  gpsbabel.pointsout(csvname, gpxname)
  if(!debug) unlink(csvname)
}
dev.off()
