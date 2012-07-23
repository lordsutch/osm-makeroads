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

library(princurve)
library(rgdal)
library(maptools)

## Set the bounding box here.

left <- -81.01
right <- -80.97
top <- 34.40
bottom <- 34.20

## What angle to split tracks at (degrees) - splits merged tracks
splitangle <- 120

## How close the angles need to be for ways to be "similar"
similarangle <- 22.5

## Max distance similar tracks can be apart (km) somewhere
maxtrackdist <- 0.025

## How close to fit
delta <- 1/1800    ## Ensure curves are reconstructed within 2 deg second
maxit <- 50

debug <- TRUE

## Return distance in km
geodetic.distance <- function(point1, point2) {
  R <- 6371
  p1rad <- point1 * pi/180
  p2rad <- point2 * pi/180
  d <- sin(p1rad[2])*sin(p2rad[2])+cos(p1rad[2])*cos(p2rad[2])*cos(abs(p1rad[1]-p2rad[1]))
  d <- acos(d)
  as.double(R*d)
}

## How much distance to split tracks by (km)
## This deals with points w/o track info.
## We'll say 25% of the diagonal distance is reasonable(?),
## with a minimum threshold of 1 km
splitdistance <- max(geodetic.distance(c(top, left), c(bottom, right))*.25, 1)

## Turn SpatialLines stuff into coordinate tracks

## There's probably an easier way to do this... maybe just parse the GPX
## file directly?

fixupLines <- function(SpLinesOb) {
  res <- sapply(SpLinesOb@lines, function(z) { lapply(z@Lines, coordinates) })

  res <- lapply(res, function(z) {
    if(is.matrix(z))
      z
    else if(length(z) > 1) # Points without track
      matrix(c(z, recursive=TRUE), ncol=2, byrow=T)
    else # Complete track
      matrix(c(z, recursive=TRUE), ncol=2, byrow=F)
  })

  res <- lapply(res, function(z) {
    z <- as.data.frame(z)
    colnames(z) <- c('x', 'y')
    z
  })
}

parseGPXfile <- function(filename) {
  converted <- try(readOGR(filename, 'tracks'), TRUE)
  
  if(inherits(converted, 'try-error')) return(NULL)

  data <- as.SpatialLines.SLDF(converted)
  fixupLines(data)
}

## Get GPX data for an area from OSM API
fetchGPXpage <- function(left, bottom, right, top, page=0) {
  coords <- paste(left, bottom, right, top, sep=",")
  gpx.url <- paste('http://api.openstreetmap.org/api/0.6/trackpoints?bbox=',
                   coords, sep='')
  pageurl <- paste(gpx.url, page, sep='&page=')

  filename1 <- tempfile(fileext='.gpx')
  print(filename1)
  download.file(pageurl, filename1)

  converted <- parseGPXfile(filename1)
  
  if(!debug) unlink(filename1)
  converted
}

getOSMtracks <- function(left, bottom, right, top) {
  page <- 0
  tracks <- c()
  repeat {
    converted <- fetchGPXpage(left, bottom, right, top, page)

    if(!is.null(converted)) {
      tracks <- c(tracks, converted)
      page <- page+1
    } else {
      break
    }
  }

  tracks
}

## Debugging code
getOSMtracksFiles <- function(filenames) {
  tracks <- c()
  for(filename in filenames) {
    converted <- parseGPXfile(filename)

    if(!is.null(converted)) {
      tracks <- c(tracks, converted)
    } else {
      break
    }
  }

  tracks
}

## Get GPS points in the area described by this bounding box
tracks <- getOSMtracks(left, bottom, right, top)

## Debugging
##tracks <- getOSMtracksFiles(list.files(pattern='file.*[.]gpx'))
##tracks <- getOSMtracksFiles(list.files(pattern='I77.*[.]gpx'))
length(tracks)

## Split tracks by distance threshold
splitTracks <- function(tracks) {
  newtracks <- list()

  for(track in tracks) {
    track <- track[!duplicated(track),]
    dmat <- spDists(as.matrix(track), longlat=TRUE)

    thistrack <- NULL
    sofar <- 0
    for(r in 1:(nrow(track)-1)) {
      ##show(r)
      r1 <- track[r,]
      r2 <- track[r+1,]
      ##show(r1)
      ##show(r2)
      dist <- dmat[r,r+1]
      ##show(dist)
      thistrack <- rbind(thistrack, r1)
      if(nrow(thistrack) > (nrow(track)/5) && dist > splitdistance) {
        show('Making new track')
        show(r1)
        show(r2)
        show(dist)
        newtracks <- c(newtracks, list(thistrack))
        thistrack <- NULL
      }
    }
    thistrack <- rbind(thistrack, r2)
    newtracks <- c(newtracks, list(thistrack))
  }

  newtracks
}
newtracks <- splitTracks(tracks)
length(newtracks)

## Split tracks by direction threshold
splitTracks2 <- function(tracks) {
  newtracks <- list()

  for(track in tracks) {
    bearings <- trackAzimuth(as.matrix(track))
    
    thistrack <- NULL
    tracklen <- 0
    for(r in 1:(nrow(track)-1)) {
      b1 <- bearings[r]
      b2 <- bearings[r+1]
      ##b3 <- bearings[r+2]      
      avbearing <- gzAzimuth(as.matrix(track[1,]),
                             as.matrix(track[r,]))

      thistrack <- rbind(thistrack, track[r,])
      difference <- min((b1-avbearing) %% 360, (avbearing-b1) %% 360)
      seglen <- geodetic.distance(track[r,], track[r+1,])
      ##show(difference)

      if((nrow(thistrack) > (nrow(track)/5) || tracklen > splitdistance)
         && abs(difference) >= splitangle) {
        ## Assume a big turn between segments is a break in the track
        if(seglen > 0.01) {
          show('** Making new track')
          show(difference)
          newtracks <- c(newtracks, list(thistrack))
          thistrack <- NULL
          tracklen <- 0
        } else {
          show('** Ignoring small deviation')
        }
      }
      tracklen <- tracklen + seglen
    }
    thistrack <- rbind(thistrack, track[r+1,])
    newtracks <- c(newtracks, list(thistrack))
  }

  newtracks
}
newtracks2 <- splitTracks2(newtracks)
length(newtracks2)

## Sort the tracks by length
sortTracks <- function(tracks) {
  lengths <- NULL
  for(track in tracks) {
    len <- 0
    for(r in 1:(nrow(track)-1)) {
      len <- len + geodetic.distance(as.matrix(track[r,]),
                                     as.matrix(track[r+1,]))
    }
    lengths <- c(lengths, len)
    ##show(lengths)
  }
  x <- sort(lengths, index.return=T, decreasing=T)
  ##show(x$ix)
  tracks[x$ix]
}
newtracks3 <- sortTracks(newtracks2)

## Identify related tracks. Tracks are related if their bearings are similar
## and they are located close enough together.
consolidateTracks <- function(tracks) {
  bearings <- sapply(tracks, function(x) {
    gzAzimuth(as.matrix(x[1,]), as.matrix(x[nrow(x),]))
  })

  tracklist <- list()
  tracklist[[1]] <- c(1)
  for(t in 2:length(tracks)) {
    found <- FALSE
    show(t)
    for(v in 1:length(tracklist)) {
      closeenough <- FALSE
      for(st in 1:length(tracklist[[v]])) {
        comparetrack <- tracklist[[v]][st]
        bdiff <- abs(bearings[t] - bearings[comparetrack])
        if(bdiff > 180) bdiff <- 360-bdiff
        if(bdiff < similarangle) closeenough <- TRUE
      }
      if(closeenough) {
        closeenough <- FALSE
        show(tracklist[[v]])
        for(ctrack in tracklist[[v]]) {
          tinfo <- as.matrix(tracks[[ctrack]])
          show(ctrack)
          for(r in 1:nrow(tracks[[t]])) {
            dist <- spDistsN1(tinfo, as.matrix(tracks[[t]][r,]), longlat=TRUE)
            show(min(dist))
            if(min(dist) < maxtrackdist) {
              closeenough <- TRUE
              break
            }
          }
          if(closeenough) break
        }
        if(closeenough) {
          tracklist[[v]] <- c(tracklist[[v]], t)
          found <- TRUE
          break
        }
      }
    }
    if(!found) tracklist[[length(tracklist)+1]] <- c(t)
  }
  tracklist
}
tracklist <- consolidateTracks(newtracks3)

## Not working yet
weighted.loess <- function(x, y, f=2/3, iter=3, delta=0.01, weights=NULL) {
  show(y)
  loess(y ~ x, span=f, iterations=iter+1, degree=1,
        surface='direct', family='symmetric', weights=weights)$y
}

sdistances <- function(curve, track) {
  dvec <- vector()
  for(p in 1:nrow(track)) {
    dvec <- c(dvec, geodetic.distance(curve$s[p,], track[p,]))
  }
  ##dvec
  sdists <- dvec/sd(dvec)
  sdists
}

fitpcurve <- function(track) {
  olen <- nrow(track)
  bandwidth <- round(min(0.1, 25/olen), 3)
  show(paste('Fitting curve with', olen, 'points; bandwidth',
             bandwidth))
  curve <- principal.curve(as.matrix(track), trace=T, f=bandwidth, maxit=maxit,
                           delta=delta, iter=2, smoother='lowess')

  ## Screen out outliers and reestimate curve
  ## Really should use weights...
  d <- sdistances(curve, track)
  track <- track[-(d >= 4),]
  if(nrow(track) > 0 && nrow(track) < olen) {
    show(paste('Refitting curve with', nrow(track), 'points'))
    
    curve <- principal.curve(as.matrix(track), ## start=curve$s[-(d>=4),],
                             trace=T, f=bandwidth,
                             maxit=maxit,
                             delta=delta, iter=2, smoother='lowess')
  } else {
    show('Empty refit?')
  }
  curve
}

gpsbabel.out <- function(infile, outfile) {
  system2('gpsbabel', c('-t', '-i', 'unicsv', '-f', infile,
                        ## Remove nearby points that don't contribute much
                        ## ~4m is the width of a lane
                        '-x', 'simplify,error=0.003k',
                        ## Ensure we have a point every 250m regardless
                        '-x', 'interpolate,distance=0.25k',
                        '-o', 'gpx', '-F', outfile))
}

tcount <- 0
pdf(onefile=T, height=8, width=6)
for(group in tracklist) {
  tcount <- tcount + 1
  mergedtracks <- NULL
  plot(newtracks3[[1]], xlim=c(left, right), ylim=c(bottom, top), asp=0.75,
       type='l', col='gray50')
  color <- 3
  for(t in group) {
    mergedtracks <- rbind(mergedtracks, newtracks3[[t]])
    lines(newtracks3[[t]], pch='.', col=color, lwd=0.2)
    text(newtracks3[[t]][1,], labels=t, col=color)
    color <- color+1
  }

  if(nrow(mergedtracks) < 20) {
    show(paste('Skipping small group', tcount))
    next
  }
  
  ## Do the other stuff here...
  f <- fitpcurve(mergedtracks)
  lines(f$s[f$tag,], asp=0.8, lty=3, col=2)

  fname <- paste('roadway-', tcount, sep='')
  csvname <- paste(fname, 'csv', sep='.')
  gpxname <- paste(fname, 'gpx', sep='.')
  ## Export the curve as a CSV file
  write.csv(f$s[f$tag,], file=csvname)

  ## Then run GPSBabel to convert the CSV to GPX
  gpsbabel.out(csvname, gpxname)
}
dev.off()
