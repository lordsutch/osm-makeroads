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
library(geosphere)

## Set the bounding box here.

left <- -85.60
right <- -85.55
top <- 40.22
bottom <- 40.07

## What angle to split tracks at (degrees) - splits merged tracks
splitangle <- 60

## How close the angles need to be for ways to be "similar"
similarangle <- 30

## Max distance similar tracks can be apart (m) somewhere
## Realistically the same road usually gets within 1 meter at least once...
maxtrackdist <- 10

## Max distance a track can ever be away from the bundle to be distinct (m)
maxseparation <- 100

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

## Turn SpatialLines stuff into coordinate tracks

## There's probably an easier way to do this... maybe just parse the GPX
## file directly?

fixupLines <- function(SpLinesOb) {
  tracks <- list()
  
  for(linesOb in SpLinesOb@lines) {
    for(lineOb in linesOb@Lines) {
      coords <- coordinates(lineOb)
      coords <- data.frame(lon=coords[,1], lat=coords[,2])
      tracks <- c(tracks, list(coords))
    }
  }
  tracks
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
tracks <- getOSMtracksFiles(list.files(pattern='I69.*[.]gpx'))
length(tracks)

## Split tracks by distance threshold
splitTracks <- function(tracks) {
  newtracks <- list()

  for(track in tracks) {
    track <- track[!duplicated(track),]

    # Convert to same scale as geosphere stuff. Grr.
    dmat <- spDists(as.matrix(track), longlat=TRUE)*1000

    thistrack <- NULL
    sofar <- 0
    for(r in 1:(nrow(track)-1)) {
      ##show(r)
      r1 <- track[r,]
      r2 <- track[r+1,]
      ##show(r1)
      ##show(r2)
      ##show(dist)
      seglen <- dmat[r,r+1]
      thistrack <- rbind(thistrack, r1)
      if(seglen > splitdistance) {
        show('Making new track')
        show(r1)
        show(r2)
        show(seglen)
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
    ## Skip degenerate tracks
    if(nrow(track) < 2)
      next
      
    bearings <- trackAzimuth(as.matrix(track))
    
    thistrack <- track[1,]
    tracklen <- 1
    for(r in 2:(nrow(track)-2)) {
      b1 <- bearings[r]
      b2 <- bearings[r+1]
      #b3 <- bearings[r+2]      
      avbearing <- gzAzimuth(as.matrix(track[max(1,r-20),]),
                             as.matrix(track[r,]))

      thistrack <- rbind(thistrack, track[r,])
      difference <- min(min((b1-avbearing) %% 360, (avbearing-b1) %% 360),
                        min((b2-avbearing) %% 360, (avbearing-b2) %% 360))
      seglen <- distHaversine(track[r,], track[r+1,])
      
      if(nrow(thistrack) > 5 && abs(difference) >= splitangle) {
        ## Assume a big turn between segments is a break in the track
        ## If points within 10m, ignore (probably stationary GPS error)
        if(seglen > 10) {
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
    thistrack <- rbind(thistrack, track[r+1,], track[r+2,])
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
      seglen <- distHaversine(as.matrix(track[r,]), as.matrix(track[r+1,]))
      len <- len + seglen
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
    for(angle in c(0, 180)) {
      ## Try 180-degrees off as a last resort
      found <- FALSE
      for(v in 1:length(tracklist)) {
        closeenough <- FALSE
        for(st in 1:length(tracklist[[v]])) {
          comparetrack <- tracklist[[v]][st]
          bdiff <- abs(bearings[t] - bearings[comparetrack])
          bdiff <- (bdiff+angle) %% 360
          if(bdiff > 180) bdiff <- 360-bdiff
          show(bdiff)
          if(bdiff <= similarangle) closeenough <- TRUE
        }
        if(closeenough) {
          show(paste('Testing', t, 'against group', v))
          closeenough <- FALSE
          show(tracklist[[v]])
          tinfo <- NULL
          for(ctrack in tracklist[[v]]) {
            dists <- 0
            tinfo <- as.matrix(tracks[[ctrack]])
            d <- dist2Line(as.matrix(tracks[[t]]), tinfo)
            dist <- d[,1]
            if(min(dist) < maxtrackdist) {
              show(paste(t, 'is within',min(dist),'of',v))
              if(max(dist) > maxseparation) {
                show(paste(t, 'too far away', max(dist),'from track', ctrack))
                closeenough <- FALSE
                break
              } else {
                closeenough <- TRUE
              }
            }
          }
          if(closeenough) {
            tracklist[[v]] <- c(tracklist[[v]], t)
            found <- TRUE
            break
          }
        }
        if(found) break;
      }
      if(found) break;
    }
    if(!found) tracklist[[length(tracklist)+1]] <- c(t)
  }
  tracklist
}
tracklist <- consolidateTracks(newtracks3)

interpolateTracks <- function(tracks) {
  if(!interpolate) {
    tracks
  } else {
    newtracks <- list()
    tnum <- 0
    for(track in tracks) {
      newtrack <- NULL
      tnum <- tnum+1
      len <- 0
      for(r in 1:(nrow(track)-1)) {
        seglen <- distHaversine(as.matrix(track[r,]), as.matrix(track[r+1,]))
        if(seglen > 1000) {
          show(tnum)
          show(seglen)
        }
        len <- len + seglen
        newtrack <- rbind(newtrack, track[r,])
        if(seglen > interpolatedist) {
          points <- round(seglen/interpolatedist)+1
          show(paste(tnum, r, 'interpolating', points))
          newpoints <- gcIntermediate(as.matrix(track[r,]),
                                      as.matrix(track[r+1,]), points)
          newpoints <- data.frame(lon=newpoints[,1], lat=newpoints[,2])
          newtrack <- rbind(newtrack, newpoints)
        }
      }
      newtrack <- rbind(newtrack, track[nrow(track),])
      newtracks <- c(newtracks, list(newtrack))
    }
  }
  newtracks
}
newtracks4 <- interpolateTracks(newtracks3)

save(tracks, newtracks, newtracks2, newtracks3, newtracks4, tracklist,
     file='I69-SR67.Rdata')

##load('testing.Rdata')

showtracks <- function(tracks, basetrack, others) {
  plot(tracks[[basetrack]], xlim=c(left, right), ylim=c(bottom, top),
       asp=0.75, type='l', col='gray50')
  color <- 3
  for(t in others) {
    lines(tracks[[t]], pch='.', col=color, lwd=0.2)
    text(tracks[[t]][1,], labels=t, col=color)
    color <- color+1
  }
}

## Not working yet
weighted.loess <- function(x, y, f=2/3, iter=3, delta=0.01, weights=NULL) {
  show(y)
  loess(y ~ x, span=f, iterations=iter+1, degree=1,
        surface='direct', family='symmetric', weights=weights)$y
}

sdistances <- function(curve, track) {
  dvec <- NULL
  for(r in 1:nrow(track)) {
    dvec <- c(dvec, distHaversine(track[r,], curve$s[r,]))
  }
  sdists <- dvec/sd(dvec)
  sdists
}

fitpcurve <- function(track) {
  olen <- nrow(track)
  bandwidth <- round(min(0.1, 40/olen), 3)
  show(paste('Fitting curve with', olen, 'points; bandwidth',
             bandwidth))
  curve <- principal.curve(as.matrix(track), trace=T, f=bandwidth, maxit=maxit,
                           delta=delta, iter=0, smoother='lowess')

  ## Screen out outliers and reestimate curve
  ## Really should use weights...
  d <- sdistances(curve, track)
  track <- track[-(d >= 4),]
  if(nrow(track) > 0 && nrow(track) < olen) {
    show(paste('Refitting curve with', nrow(track), 'points'))
    
    curve <- principal.curve(as.matrix(track), ## start=curve$s[-(d>=4),],
                             trace=T, f=bandwidth,
                             maxit=maxit,
                             delta=delta, iter=0, smoother='lowess')
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
