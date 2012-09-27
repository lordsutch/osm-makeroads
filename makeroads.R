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
  download.file(pageurl, filename1, cacheOK=FALSE)

  converted <- parseGPXfile(filename1)
  
  if(!debug) unlink(filename1)
  converted
}

getOSMtracks <- function(left, bottom, right, top) {
  page <- 0
  tracks <- c()
  if(right < left) {
    temp <- right
    right <- left
    left <- temp
  }
  if(top < bottom) {
    temp <- top
    top <- bottom
    bottom <- temp
  }
  
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

## Split tracks by distance threshold
splitTracks <- function(tracks) {
  newtracks <- list()

  for(track in tracks) {
    track <- track[!duplicated(track),]

    # Convert to same scale as geosphere stuff. Grr.
    dmat <- spDists(as.matrix(track), longlat=TRUE)*1000

    thistrack <- NULL
    sofar <- 0
    if(nrow(track) == 1) {
      newtracks <- c(newtracks, list(track))
      next
    }

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
        sofar <- 0
      } else {
        sofar <- sofar + seglen
      }
    }
    thistrack <- rbind(thistrack, r2)
    newtracks <- c(newtracks, list(thistrack))
  }

  newtracks
}

## Split tracks by direction threshold
splitTracks2 <- function(tracks) {
  newtracks <- list()

  for(track in tracks) {
    ignorecount <- 0
    
    if(nrow(track) <= 2) {
      newtracks <- c(newtracks, list(track))
      next
    }
    
    ##bearings <- trackAzimuth(as.matrix(track))
    
    thistrack <- track[1,]
    tracklen <- 0
    for(r in 2:(nrow(track)-1)) {
      thistrack <- rbind(thistrack, track[r,])

      spos = max(nrow(thistrack)-2, 1)
      
      avbearing <- gzAzimuth(as.matrix(thistrack[spos,]),
                             as.matrix(track[r,]))
      avbearing2 <- gzAzimuth(as.matrix(thistrack[1,]),
                              as.matrix(track[r,]))
      
      difference <- min((avbearing2-avbearing) %% 360,
                        (avbearing-avbearing2) %% 360)

      seglen <- distHaversine(track[r,], track[r+1,])
      
      if(nrow(thistrack) > 1 && abs(difference) >= splitangle) {
        ## Assume a big turn between segments is a break in the track
        ## If points within 25m, ignore (probably stationary GPS error)
        show(difference)
        if(seglen > 25 || ignorecount > 5 || tracklen > 1000) {
          show('** Making new track')
          newtracks <- c(newtracks, list(thistrack))
          thistrack <- track[r,]
          tracklen <- 0
          ignorecount <- 0
        } else {
          ##show('** Ignoring small deviation')
          ignorecount <- ignorecount+1
        }
      }
      tracklen <- tracklen + seglen
    }
    thistrack <- rbind(thistrack, track[r+1,])
    newtracks <- c(newtracks, list(thistrack))
  }

  newtracks
}

## Sort the tracks by length
sortTracks <- function(tracks) {
  tcount <- length(tracks)
  lengths <- numeric(tcount)
  for(i in seq(tcount)) {
    track <- tracks[[i]]
    len <- 0
    if(nrow(track) >= 2) {
      for(r in 1:(nrow(track)-1)) {
        seglen <- distHaversine(as.matrix(track[r,]), as.matrix(track[r+1,]))
        len <- len + seglen
      }
    }
    lengths[i] <- len ## *log(nrow(track))
  }
  x <- sort(lengths, index.return=T, decreasing=T)
  ##show(x$ix)
  tracks[x$ix]
}

simplifyTracks <- function(tracks) {
  newtracks <- NULL
  for(track in tracks) {
    if(nrow(track) >= 3) {
      infile <- tempfile(fileext='.csv')
      outfile <- tempfile(fileext='.csv')
      
      write.csv(round(track, 6), file=infile)

      system2('gpsbabel', c('-t', '-i', 'unicsv', '-f', infile,
                            ## Remove nearby points that don't contribute much
                            ## ~4m is the width of a lane
                            '-x', 'simplify,error=0.001k',
                            '-o', 'unicsv', '-F', outfile))

      track <- read.csv(outfile)
      track <- data.frame(lon=track$Longitude, lat=track$Latitude)
      newtracks <- c(newtracks, list(track))
      unlink(infile)
      unlink(outfile)
    } else {
      newtracks <- c(newtracks, list(track))
    }
  }
  newtracks
}

calcbearings <- function(x) {
  if(nrow(x) >= 2)
    gzAzimuth(as.matrix(x[1,]), as.matrix(x[nrow(x),]))
  else
    0
}

## Identify related tracks. Tracks are related if their bearings are similar
## and they are located close enough together.
consolidateTracks <- function(tracks) {
  bearings <- sapply(tracks, calcbearings)

  tracklist <- list()
  tracklist[[1]] <- c(1)
  newtracks <- tracks
  loosepoints <- NULL
  t <- 2
  while(t <= length(newtracks)) {
    tryagain <- FALSE
    for(angle in c(0, 180)) {
      ## Try 180-degrees off as a last resort
      found <- FALSE
      if(angle == 180)
        separation <- maxoppositeseparation
      else
        separation <- maxseparation
        
      for(v in 1:length(tracklist)) {
        closeenough <- FALSE
        for(st in 1:length(tracklist[[v]])) {
          comparetrack <- tracklist[[v]][st]
          if(nrow(newtracks[[comparetrack]]) < 2) {
            closeenough <- TRUE
            break
          }
          
          bdiff <- abs(bearings[t] - bearings[comparetrack])
          bdiff <- (bdiff+angle) %% 360
          if(bdiff > 180) bdiff <- 360-bdiff

          if(bdiff <= similarangle) {
            closeenough <- TRUE
            break
          }
        }
        if(closeenough) {
          show(paste('Testing', t, 'against group', v))
          closeenough <- FALSE
          show(tracklist[[v]])
          for(ctrack in tracklist[[v]]) {
            if(nrow(newtracks[[ctrack]]) < 2) {
              if(nrow(newtracks[[t]]) < 2) {
                dist <- distHaversine(newtracks[[t]], newtracks[[ctrack]])
              } else {
                d <- dist2Line(newtracks[[ctrack]], newtracks[[t]])
                dist <- d[,1]
              }
            } else {
              d <- dist2Line(newtracks[[t]], newtracks[[ctrack]])
              dist <- d[,1]
            }
            show(dist)
            fdist <- dist
            if(length(dist) >= 4)
              dist <- dist[2:(length(dist)-1)] ## Ignore endpoints
            
            if(min(dist) < maxtrackdist) {
              show(paste(t, 'is within',min(dist),'of', ctrack))
              if(max(dist) > separation) {
                show(paste(t, 'too far away', max(dist),'from track', ctrack))
                runs <- rle(fdist <= separation)
                ## Consider splitting, but only if the overlap is long enough
                if(any(runs$lengths[runs$values] > 2)) {
                  closetracks <- NULL
                  fartracks <- NULL
                  points <- NULL
                  j <- 1
                  show(runs)
                  track <- newtracks[[t]]
                  for(i in 1:length(runs$lengths)) {
                    ##show(i)
                    ##show(runs$lengths[i])
                    trackseg <- track[j:(j+runs$lengths[i]-1),]
                    j <- j + nrow(trackseg)
                    ##show(nrow(trackseg))
                    if(nrow(trackseg) > 1) {
                      if(runs$values[i])
                        closetracks <- c(closetracks, list(trackseg))
                      else
                        fartracks <- c(fartracks, list(trackseg))
                    } else {
                      points <- c(points, list(trackseg))
                    }
                  }
                  ##print(length(closetracks))
                  ##print(length(fartracks))
                  
                  closeenough <- FALSE
                  xtracks <- c(newtracks[1:(t-1)], closetracks, fartracks)
                  xbearings <- bearings[1:(t-1)]
                  if(length(closetracks))
                    xbearings <- c(xbearings, sapply(closetracks, calcbearings))
                  if(length(fartracks))
                    xbearings <- c(xbearings, sapply(fartracks, calcbearings))
                  if(t < length(newtracks)) {
                    xtracks <- c(xtracks, newtracks[(t+1):length(newtracks)])
                    xbearings <- c(xbearings, bearings[(t+1):length(bearings)])
                  }
                  ## Points at the end
                  newtracks <- c(xtracks, points)
                  bearings <- c(xbearings, rep(0, length(points)))
                  ## loop again without incrementing t; why we're not using for
                  tryagain <- TRUE
                }
                ## Just fall through to the next group
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
    if(!found) {
      if(nrow(newtracks[[t]]) >= 2)
        tracklist[[length(tracklist)+1]] <- c(t)
      else
        loosepoints <- c(loosepoints, t)
    }
    if(!tryagain)
      t <- t + 1
  }
  ret <- list(tracklist=tracklist, tracks=newtracks, loosepoints=loosepoints)
  ret
  ##tracklist
}

interpolateTracks <- function(tracks) {
  if(!interpolate) {
    tracks
  } else {
    newtracks <- list()
    tnum <- 0
    for(track in tracks) {
      if(nrow(track) < 2) {
        newtracks <- c(newtracks, list(track))
        next
      }
      
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
          ##show(paste(tnum, r, 'interpolating', points))
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
  bandwidth <- round(min(0.1, 15/olen), 3)
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
