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

library(rgeos)
library(rgdal)
library(maptools)
library(geosphere)
##library(LPCM)
library(princurve)

parseGPXfile <- function(filename, tracksfrompoints=FALSE) {
  timesep <- 60  ## Split if points are more than 'timesep' seconds apart
  distsep <- 100 ## Only make tracks if points less than this distance apart

  converted <- try(readOGR(filename, 'track_points'), TRUE)

  if(inherits(converted, 'try-error')) return(NULL)

  proj <- proj4string(converted)
  ##as.SpatialLines.SLDF(converted)
  trackrows <- !is.na(converted$time)
  points <- coordinates(converted)[!trackrows,,drop=FALSE]
  lines <- NULL

  ## Try to assemble anonymized points into lines
  linenum <- 1
  if(tracksfrompoints && !is.null(points) && nrow(points) >= 2) {
    start <- 1
    linelist <- NULL
    npoints <- NULL
    
    for(j in seq_len(nrow(points)-1)) {
      dist <- distHaversine(points[j,], points[j+1,])
      if(dist > distsep) {
        seg <- points[start:j,,drop=FALSE]
        if(nrow(seg) > 1) {
          linelist <- c(linelist, Line(seg))
        } else {
          npoints <- rbind(npoints, seg)
        }
        start <- j+1
      }
    }
    
    if(!is.null(linelist)) {
      lines <- c(lines, Lines(linelist, linenum))
      linenum <- linenum+1
    }
    points <- npoints
  }
  
  alltracks <- converted[trackrows,]
  for(tnum in unique(alltracks$track_fid)) {
    thistrack <- alltracks[alltracks$track_fid == tnum,]
    for(tseg in unique(thistrack$track_seg_id)) {
      thisseg <- thistrack[thistrack$track_seg_id == tseg,]
      times <- strptime(as.character(thisseg$time), format="%Y/%m/%d %H:%M:%S",
                        tz='GMT')
      ltimes <- c(tail(times, -1), NA)
      tdiff <- (ltimes-times)
      l <- nrow(thisseg)

      splits <- which(tdiff > timesep)
      start <- 1
      linelist <- NULL
      for(seg in splits) {
        segpoints <- coordinates(thisseg)[start:seg,,drop=FALSE]
        ##str(segpoints)
        if(!is.null(segpoints)) {
          if(nrow(segpoints) > 1) {
            linelist <- c(linelist, Line(segpoints))
          } else {
            points <- rbind(points, segpoints)
          }
        }
        start <- seg+1
      }
      segpoints <- coordinates(thisseg)[start:nrow(thisseg),,drop=FALSE]
      if(!is.null(segpoints)) {
        if(nrow(segpoints) > 1) {
          linelist <- c(linelist, Line(segpoints))
        } else {
          points <- rbind(points, segpoints)
        }
      }
      
      if(!is.null(linelist)) {
        lines <- c(lines, Lines(linelist, linenum))
        linenum <- linenum+1
      }
    }
  }

  tracks <- NULL
  spoints <- NULL

  str(lines)
  str(points)
  
  if(!is.null(lines))
    tracks <- SpatialLines(lines, proj4string=CRS(proj))
  if(!is.null(points))
    spoints <- SpatialPoints(points, proj4string=CRS(proj))
  list(tracks=tracks, points=spoints)
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

  converted <- parseGPXfile(filename1, tracksfrompoints=closetracks)

  if(!debug) unlink(filename1)
  converted
}

getOSMtracks <- function(left, bottom, right, top) {
  page <- 0
  tracks <- NULL
  points <- NULL
  
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
      if(!is.null(converted$tracks)) {
        if(!is.null(tracks))
          tracks <- rbind(tracks, converted$tracks, makeUniqueIDs=TRUE)
        else
          tracks <- converted$tracks
      }
      if(!is.null(converted$points)) {
        if(!is.null(points))
          points <- rbind(points, converted$points)
        else
          points <- converted$points
      }
      page <- page+1
    } else {
      break
    }
  }
  list(tracks=tracks, points=points)
}

## Open GPX files directly (bypassing API download step)
getOSMtracksFiles <- function(...) {
  filenames <- list(...)
  tracks <- NULL
  points <- NULL
  for(filename in filenames) {
    converted <- parseGPXfile(filename, tracksfrompoints=closetracks)
    if(!is.null(converted)) {
      if(is.null(tracks)) {
        tracks <- converted$tracks
      } else {
        tracks <- rbind(tracks, converted$tracks, makeUniqueIDs=TRUE)
      }
      if(is.null(points)) {
        points <- converted$points
      } else {
        points <- rbind(points, converted$points)
      }
    } else {
      break
    }
  }
  list(tracks=tracks, points=points)
}

## Split tracks by direction threshold
splitTracksAngle <- function(tracks, splitangle, splitdistance) {
  minanglesplit <- 200 # Minimum track distance before we split a track

  ## Unproject back to latitude and longitude for this step
  uptracks <- spTransform(tracks, CRS("+proj=latlong"))
  i <- 0
  points <- NULL
  lineslist <- NULL
  for(t in seq_along(uptracks@lines)) {
    linelist <- NULL
    for(l in seq_along(uptracks@lines[[t]]@Lines)) {
      linesegs <- uptracks@lines[[t]]@Lines[[l]]
      ##str(linesegs)
      cds <- coordinates(linesegs)
      if(nrow(cds) < 3) {
        ## Use *projected* coordinates when reconstructing list
        linelist <- c(linelist, tracks@lines[[t]]@Lines[[l]])
      } else {
        pcds <- coordinates(tracks@lines[[t]]@Lines[[l]])
        d <- LineLength(uptracks@lines[[t]]@Lines[[l]], longlat=TRUE,
                        sum=FALSE)*1000
        az <- trackAzimuth(cds, type='abdali')
        if(length(az) > 3) {
          laz <- c(tail(az, -3), NA, NA, NA)
          ldiff <- pmin((az-laz) %% 360, (laz-az) %% 360)
          ldiff[is.na(ldiff)] <- 0
          str(ldiff)
        } else {
          ldiff <- rep(0, length(az))
        }
        repeat {
          ##str(d)
          firstangle <- az[1]
          
          adiff <- c(0,pmin((az-firstangle) %% 360, (firstangle-az) %% 360))
          adiff <- pmax(adiff, c(0,ldiff))

          distfromstart <- c(0, cumsum(d))
          ##str(distfromstart)
          x <- which(adiff > splitangle & distfromstart > minanglesplit)
          str(x)
          if(length(x) > 0 && x[[1]] == 1) {
            x <- x[-1]
          }
          if(length(x) < 1) {
            if(nrow(pcds) > 1)
              linelist <- c(linelist, Line(pcds))
            else
              points <- rbind(points, pcds)
            break
          }
          pos <- x[1]
          bit <- pcds[1:pos,,drop=FALSE]
          if(nrow(bit) > 1) {
            linelist <- c(linelist, Line(bit))
          } else {
            points <- rbind(points, bit)
          }
          cds <- cds[-(1:pos-1),,drop=FALSE]
          pcds <- pcds[-(1:pos-1),,drop=FALSE]
          az <- az[-(1:pos-1),drop=FALSE]
          ldiff <- ldiff[-(1:pos-1), drop=FALSE]
          d <- d[-(1:pos-1),drop=FALSE]
        }
      }
    }
    if(!is.null(linelist)) {
      newlines <- Lines(linelist, as.character(i))
      lineslist <- c(lineslist, newlines)
      i <- i+1
    }
  }

  if(!is.null(points)) {
    spoints <- SpatialPoints(points, proj4string=CRS(proj4string(tracks)))
  } else {
    spoints <- NULL
  }

  list(tracks=SpatialLines(lineslist, proj4string=CRS(proj4string(tracks))),
       points=spoints)
}

## Sort the tracks by length
sortTracks <- function(tracks, latlon=FALSE) {
  lengths <- SpatialLinesLengths(tracks, latlon)
  x <- sort(lengths, index.return=T, decreasing=T)
  ##show(x$ix)
  tracks[x$ix,]
}

## get the coordinates out of a spatiallines as a flat object
coordinatesSL <- function(x) {
  c0 <- coordinates(x)
  c1 <- lapply(c0, function(x) do.call('rbind', x))
  do.call('rbind', c1)
}

calcbearings <- function(x) {
  c2 <- coordinatesSL(x)
  gzAzimuth(c2[1,,drop=FALSE], c2[nrow(c2),,drop=FALSE],
            type='abdali')
}

findClosestTrackToPoint <- function(point, tracks, tracklist) {
  dists <- gDistance(tracks, point, byid=TRUE)
  
  mindist <- min(dists)
  closestTrack <- which.min(dists)
  show(closestTrack)
  
  if(mindist > maxpointdist)
    return(0)

  for(t in seq_along(tracklist)) {
    if(closestTrack %in% tracklist[[t]])
      return(t)
  }

  ## Failure mode
  show('Why are we here?')
  return(0)
}

flattenSpatialLines <- function(sl) {
  # flatten lines: Ensure each Lines object only has one Line in it.
  lineslist <- NULL
  label <- 1
  for(i in seq_along(sl)) {
    for(j in seq_along(sl[i]@lines)) {
      x <- sl[i]@lines[[j]]@Lines
      for(k in seq_along(x)) {
        l <- Line(coordinates(x[[k]]))
        lineslist <- c(lineslist, Lines(list(l), as.character(label)))
        label <- label+1
      }
    }
  }
  SpatialLines(lineslist, proj4string=CRS(proj4string(sl)))
}

## Identify related tracks. Tracks are related if their bearings are similar
## and they are located close enough together.
consolidateTracks <- function(tracks, points) {
  tracklist <- list()
  tracklist[[1]] <- c(1)

  proj <- proj4string(tracks)

  bounds <- bbox(gEnvelope(tracks))
  ylim <- c(bounds['y','min'], bounds['y','max'])
  xlim <- c(bounds['x','min'], bounds['x','max'])
  
  ##newtracks <- flattenSpatialLines(tracks)

  ## Sort by length
  ##newtracks <- sortTracks(newtracks)
  newtracks <- tracks
  
  bearings <- integer(length(newtracks))
  for(i in seq_along(newtracks)) {
    bearings[[i]] <- calcbearings(newtracks[i])
  }

  t <- 2
  while(t <= length(newtracks)) {
    tryagain <- FALSE
    found <- FALSE
    closeenough <- FALSE

    track <- newtracks[t]
    
    alist <- c(0, 180)
    for(angle in alist) {
      if(tryagain)
        break

      ## Try 180-degrees off as a last resort
      found <- FALSE
      if(angle == 180)
        separation <- maxoppositeseparation
      else
        separation <- maxseparation

      for(v in seq_along(tracklist)) {
        if(tryagain)
          break
        closeenough <- FALSE

        show(paste('Testing', t, 'of', length(newtracks),'against group', v))

        ztracks <- newtracks[tracklist[[v]]]
        ch <- gConvexHull(ztracks)
        buffzone <- gBuffer(ch, width=maxtrackdist)
        if(!gCrosses(track, buffzone)) next

        plot(buffzone, xlim=xlim, ylim=ylim)
        lines(ztracks, col='gray50')
        lines(track)

        possible <- NULL
        for(ctrack in tracklist[[v]]) {
          bdiff <- abs(bearings[t] - bearings[ctrack])
          if(is.na(bdiff) || (angle == 0 && bdiff <= similarangle)) {
            possible <- cbind(possible, ctrack)
            closeenough <- TRUE
            break
          } else if (angle == 180 &&
                     abs(180-bdiff) >= abs(180-similarangle)) {
            possible <- cbind(possible, ctrack)
            closeenough <- TRUE
            break
          }
        }

        if(!closeenough) next

        show('Angles close')
        closeenough <- FALSE

        for(ctrack in tracklist[[v]]) {
          buffzone <- gBuffer(newtracks[ctrack], width=maxtrackdist)
          fatzone <- gBuffer(newtracks[ctrack], width=separation)

          if(!gCrosses(track, buffzone))
            next

          show(paste("Considering track", ctrack))
          show(paste(t, 'is close enough to', v))
          
          ispartoutside <- gCrosses(track, fatzone)
          if(ispartoutside) {
            insidepart <- gIntersection(track, fatzone)
            outsidepart <- gDifference(track, fatzone)

            ##str(insidepart)
            ##str(outsidepart)
            
            lines(newtracks[ctrack], col='green')
            if(!is.null(outsidepart))
              lines(outsidepart, col='red')
            if(!is.null(insidepart))
              lines(insidepart, col='blue')
            
            show(paste(t, 'too far away from group', v))

            insidepart <- flattenSpatialLines(insidepart)
            outsidepart <- flattenSpatialLines(outsidepart)

            ll <- SpatialLinesLengths(insidepart)
            insidetracks <- (ll > 100)
            
            if(sum(SpatialLinesLengths(outsidepart)) < mintracklen) {
              show('Accepting small piece outside.')
              ## Forgive a small part outside the envelope
              closeenough <- TRUE
              break
            } else if(any(insidetracks)) {
              keepinside <- flattenSpatialLines(insidepart[insidetracks])

              if(length(insidepart) < length(insidepart[insidetracks]))
                outsidepart <- rbind(outsidepart,
                                     insidepart[-insidetracks],
                                     makeUniqueIDs=TRUE)
              ##outsidepart <- gLineMerge(outsidepart)
              outsidepart <- flattenSpatialLines(outsidepart)

              show(SpatialLinesLengths(keepinside))
              show(SpatialLinesLengths(outsidepart))
              
              oldtracks <- newtracks
              if(t > 1)
                newtracks <- rbind(oldtracks[1:(t-1)], keepinside,
                                   makeUniqueIDs=TRUE)
              else
                newtracks <- keepinside
              show(length(newtracks))
              
              if(t < length(oldtracks))
                newpart <- rbind(outsidepart,
                                 oldtracks[(t+1):length(oldtracks)],
                                 makeUniqueIDs=TRUE)
              else
                newpart <- outsidepart

              newpart <- sortTracks(newpart)
              newtracks <- rbind(newtracks, newpart, makeUniqueIDs=TRUE)
              
              bearings <- integer(length(newtracks))
              for(i in seq_along(newtracks)) {
                bearings[[i]] <- calcbearings(newtracks[i])
              }

              ntlist <- t:(t+length(keepinside)-1)
              show(ntlist)
              tracklist[[v]] <- c(tracklist[[v]], ntlist)
              t <- t+length(ntlist)-1
              show(t)
              tryagain <- FALSE
              found <- TRUE
              closeenough <- FALSE
              break
            } else {
              show("Too-short overlap; moving on to next group.")
              closeenough <- FALSE
              tryagain <- FALSE
              found <- FALSE
              break
            }
          } else {
            closeenough <- TRUE
            break
          }
        }
        if(closeenough) {
          tracklist[[v]] <- c(tracklist[[v]], t)
          found <- TRUE
          break
        }
        if(found) break
      }
      if(found) break
    }
    if(!tryagain && !found) {
      center <- gCentroid(track)
      if(!is.null(center)) {
        proj4string(center) <- CRS(proj)
        if(gCoveredBy(track, gBuffer(center, width=mintracklen))) {
          ## Degenerate track
          show(paste('Converting short track', t, 'to points.'))
          tpoints <- coordinatesSL(track)
          points <- rbind(points, SpatialPoints(tpoints, proj4string=CRS(proj)))
          newtracks <- newtracks[-t]
          bearings <- bearings[-t]
          tryagain <- TRUE
        } else {
          tracklist[[length(tracklist)+1]] <- c(t)
        }
      }
    }
    if(!tryagain)
      t <- t+1
    show(t)
  }

  points <- rbind(points, SpatialPoints(coordinatesSL(newtracks),
                                        proj4string=CRS(proj)))
  
  trackforpoints <- integer(0)
  if(length(points)) {
    show('Finding closest tracks for points')
    trackforpoints <- integer(length(points))
    for(i in seq_along(points)) {
      trackforpoints[[i]] <- findClosestTrackToPoint(points[i], newtracks, tracklist)
    }

    loosepoints <- points[trackforpoints == 0]
  } else {
    loosepoints <- NULL
  }

  ret <- list(tracklist=tracklist, trackforpoints=trackforpoints,
              tracks=newtracks, points=points, loosepoints=loosepoints)
  ret
  ##tracklist
}

showtracks <- function(tracks, basetrack, others) {
  pbounds <- bbox(tracks)
  plot(tracks[basetrack], col='gray50', xlim=pbounds['x',], ylim=pbounds['y',])
  color <- 3
  for(t in others) {
    lines(tracks[t], pch='.', col=color, lwd=0.2)
    text(tracks[t], labels=t, col=color)
    color <- color+1
  }
}

sdistances <- function(curve, track, projection) {
  ptrack <- project(as.matrix(track), projection, inv=TRUE)
  ctrack <- project(as.matrix(curve), projection, inv=TRUE)

  ##str(ptrack)
  ##str(ctrack)
  
  dvec <- distHaversine(ptrack, ctrack)
  str(dvec)
  str(sd(dvec))
  
  sdists <- dvec/sd(dvec)
  sdists
}

setupForFit <- function(points, tracklist, trackforpoints, tcount) {
  thesepoints <- points[trackforpoints == tcount]
  coordinates(thesepoints)
}

fitpcurve <- function(track, projection) {
  olen <- nrow(track)
  bandwidth <- round(min(0.1, 30/olen), 3)
  show(paste('Fitting curve with', olen, 'points; bandwidth',
             bandwidth))
  curve <- principal.curve(as.matrix(track), trace=T, f=bandwidth, maxit=maxit,
                           thresh=1/3600,
                           delta=delta, iter=1, smoother='lowess')

  ## Screen out outliers and reestimate curve
  ## Really should use weights...
  d <- sdistances(curve, track, projection)
  show(d)
  track <- track[-(d >= 4),]
  if(nrow(track) > 0 && nrow(track) < olen) {
    show(paste('Refitting curve with', nrow(track), 'points'))

    curve <- principal.curve(as.matrix(track), ## start=curve$s[-(d>=4),],
                             trace=T, f=bandwidth,
                             thresh=1/3600,
                             maxit=maxit,
                             delta=delta, iter=1, smoother='lowess')
  } else {
    show('Empty refit?')
  }
  curve[curve$tag,]
}

## Alternative using local principal curves algorithm. At present,
## more fiddly than the loess fit, so disabled.

## Could weight using hdop if available...
fitpcurve.lpc <- function(track, projection, bandwidth=0.14) {
  track <- as.matrix(track)
  olen <- nrow(track)
  show(paste('Fitting curve with', olen, 'points; bandwidth', bandwidth))

  x0 <- 0 #track[sample(olen, 10),]
  #str(x0)
  spoints <- max(floor(olen/50), 5)
  
  curve <- lpc(track, h=bandwidth, scaled=TRUE, pen=3, depth=3, x0=x0,
               control=lpc.control(mult=spoints))
  ##plot(curve)
  ucurve <- lpc.spline(curve, project=TRUE)
  ucurve <- unscale(ucurve)
  str(ucurve)
  
  d <- sdistances(ucurve$closest.coords, track, projection)
  show(which(d >= 3))
  weights <- pmin(pmax((4-d), 0.2), 1)
  str(d)
  str(weights)

  if(FALSE) { # any(weights) < 1
    show(paste('Refitting weighted curve'))
    
    curve <- lpc(track, h=bandwidth, scaled=TRUE, weights=weights,
                 depth=3, pen=3, x0=x0,
                 control=lpc.control(mult=spoints))
    ##plot(curve)
    ucurve <- lpc.spline(curve, project=TRUE)
    ucurve <- unscale(ucurve)
  }
  ucurve
}

gpsbabel.out <- function(infile, outfile) {
  system2('gpsbabel', c('-t', '-i', 'unicsv', '-f', infile,
                        ## Remove nearby points that don't contribute much
                        ## ~4m is the width of a lane
                        '-x', 'simplify,error=0.001k',
                        ## Ensure we have a point every 250m regardless
                        '-x', 'interpolate,distance=0.25k',
                        '-o', 'gpx', '-F', outfile))
}

gpsbabel.pointsout <- function(infile, outfile) {
  system2('gpsbabel', c('-w', '-i', 'unicsv', '-f', infile,
                        '-o', 'gpx', '-F', outfile))
}
