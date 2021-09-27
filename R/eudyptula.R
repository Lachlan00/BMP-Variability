# Penguin Cleaning
# modules
library(rgdal)
library(plyr)
library(dplyr)
library(stringr)
library(rlist)
library(ggplot2)
library(grid)
library(gridExtra)
library(sp)
library(geosphere)
library(data.table)
library(tools)
library(lubridate)
library(readxl)
library(crawl)
library(pbapply)
library(RANN)
#library(RGeostats)
library(terra)
library(ncdf4)
library(raster)

# load local 
source('R/visual.r')
source('R/utilities.r')

# __DEV_NOTES__
#
# Ocean Point function:
# ---------------------
# Tunes out the 'over()' method for determining if points are within a 
# shapefile polygon returns 'NA' if not on a polygon but if it does returns
# the variable data for that polygon. This is why we were getting the odd
# 76.65 value returned. Consequently code will break if other shape files used.
# Is kind of useful for cleaning land errors however so a slight adaptation may
# be useful. 
# UPDATE - I think this is fixed, I can't remember if I did or not (oh dear)...
# In any case this function could be improved greatly by using the much faster
# "point.in.polygon" function in "sp". Also would be much faster to use run 
# length encoding rather than looping over tracks with a sliding window.

################## 
# load kml track #
##################
load_kml <- function(dir){
  tryCatch(
    {
    track <- as.data.frame(readOGR(dir, verbose=FALSE))
    track <- plyr::rename(track, c('Name' = 'dtUTC', 
                             'coords.x1' = 'lon', 
                             'coords.x2' = 'lat'))
    track$Description <- NULL
    track$coords.x3 <- NULL
    track$id <- file_path_sans_ext(basename(dir))
    return(track)
    }, error = function(e){
      return(NULL)
    }
  )
}

####################################
# load all kml tracks in directory #
####################################
load_tracks <- function(dir, prog=TRUE){
  file_ls <- list.files(path=dir, pattern=".kml$", recursive=TRUE)
  fail_ls <- list()
  print(paste('Loading', length(file_ls), 'kml files..'))
  if (prog){
    pb <- txtProgressBar(min = 0, max = length(file_ls)-1, style = 3)
  }
  i <- 0
  for (file in file_ls){
    track <- load_kml(paste0(dir,file))
    if (is.null(track) == FALSE){
      if (i == 0){
        df <- track
      } else {
        df <- rbind(df, track)
      }
    } else {
      fail_ls <- list.append(fail_ls, file_path_sans_ext(basename(file)))
    }
    if (prog){
      setTxtProgressBar(pb, i)
    }
    i <- i + 1
  }
  if (length(fail_ls) > 0){
    cat('\nFailed to load files:')
    print(sapply(fail_ls, paste, collapse=":"))
  }
  message('Cleaning up data frame..')
  df$id <- as.character(sapply(as.character(df$id), function(x) strsplit(x,'_')[[1]][1]))
  df <- df[,c('id','dtUTC','lon','lat')]
  df$dtUTC <- as.POSIXct(df$dtUTC, format="%d/%m/%Y %H:%M:%S", tz='UTC')
  df$id <- as.factor(df$id)
  return(df)
}

####################
# ocean points all #
####################
ocean_points_trim <- function(df, shp="assets/shp/NSW-coastline/edited/polygon/NSW-coastline_WGS84", 
                              trim2land=FALSE, view_map=FALSE, return_ocean=FALSE){
  # by default trims to first and last point over ocean
  print('Loading coastal shapefile data..')
  coast <- readOGR(dsn=path.expand(str_sub(shp,1,-(nchar(basename(shp)))-1)), 
                   layer=basename(shp)) # load shapefile
  proj <- proj4string(coast) # grab projection information
  shp <- list(coast, proj) # pass as list
  df <- ocean_points(df, shp, trim2land=trim2land, view_map=view_map) # execute function
  if (!return_ocean){
    df$ocean <- NULL
  }
  return(df)
}

################
# ocean points #
################
ocean_points <- function(df, shp, maxwindow=600, view_map=FALSE, trim2land=FALSE){
  # ocean = 1
  # land = 0
  # Note - can't be bothered as it works fine but this function would probably be much 
  # more effcient if I just used `rle` (run length encoding) instead of the for loop.
  print('Calculating which points are over ocean..')
  ## determine if lat/lon pairs are on land or ocean
  ## ocean = 1, land = 0
  dat <- data.frame(lon = df$lon,
                    lat = df$lat)
  coordinates(dat) <- ~ lon + lat
  # use shp projection
  proj4string(dat) <- shp[[2]]
  # check if points on ocean
  ocean <- over(dat, shp[[1]])
  ocean <- ocean['area']
  ocean[is.na(ocean)] <- 1
  # add to dataset
  df$ocean <- as.factor(ocean$area)
  # filter strange error points
  df <- df[!(df$ocean==76.65),]
  df$ocean <- factor(df$ocean)
  ## check results on map
  if (view_map){
    print('Plotting all points..')
    ocean.land.plot(df, shp[[1]], 'All tracks')
  }
  
  # Find start and end of foraging trip for each track
  print('Finding start and end points of foraging trips..')
  track_ls <- split(df, df$id)
  df_ls <- list()
  pb <- txtProgressBar(min = 0, max = length(track_ls)-1, style = 3)
  prog_count <- 0
  for (track in track_ls){
    # reset track index
    row.names(track) <- NULL
    # check if any land points
    if (nrow(track[track$ocean == 0,]) > 0){
      # establish an appropriate window size
      window <- nrow(track[track$ocean == 1,]) %/% 6
      if (window > maxwindow){
        widnow <- maxwindow
      }
      # continue if there are land points
      ocean_count <- 0
      for (j in 1:length(track$ocean)){
        # reset ocean count if land point
        if (track$ocean[j] == 0){
          ocean_count <- 0
        } else if (ocean_count == window){
          land_lower <- j - (window + 1)
          break
        } else {
          ocean_count <- ocean_count + 1
        }
      }
      # if window not found set land_lower
      if (ocean_count < window){
        land_lower <- 0
      }
      # get last point by reverseing data.frame
      track_inv <- track[nrow(track):1, ]
      # reset counter
      ocean_count <- 0
      # start counting from end
      for (j in 1:length(track_inv$ocean)){
        # reset ocean count if land point
        if (track_inv$ocean[j] == 0){
          ocean_count <- 0
        } else if (ocean_count == window){
          land_upper <- length(track_inv$ocean) - j + window
          break
        } else {
          ocean_count <- ocean_count + 1
        }
      }
      # if window not found set land_upper
      if (ocean_count < window){
        land_upper <- nrow(track)-1
      }
      # trim data.frame
      if (trim2land){
        # check if there even was land lower section
        if (land_lower == 0){
          land_lower <- 1
        }
        # check if there even was land upper section
        if (land_upper == nrow(track)-1){
          land_upper <- land_upper - 1
        }
        # trim to land
        track <- track[(land_lower):(land_upper+2),]
      } else {
        # trim to ocean
        track <- track[(land_lower+1):(land_upper+1),]
      }
      # Visual check
      if (view_map == TRUE){
        ocean.land.plot(track, shp[[1]], track$id[1])
      }
      # reset index again
      row.names(track) <- NULL
    }
    # append
    df_ls <- list.append(df_ls, track)
    # update progressbar
    setTxtProgressBar(pb, prog_count)
    prog_count <- prog_count + 1
  } 
  df <- do.call("rbind", df_ls)
  return(df)
}

##################################
# make predictive path for track #
##################################
crawl_predict <- function(track, ts, view_map=FALSE){
  # track <- track[[1]]
  row.names(track) <- 1:nrow(track)
  # set initial parmaters
  initial = list(a=c(coordinates(track)[1,1],0,
                     coordinates(track)[1,2],0),
                 P=diag(c(10000^2,54000^2,10000^2,5400^2)))
  # fix value in model
  fixPar = c(log(5),NA,NA)
  ### Fit the model ###
  set.seed(420)
  fit1 <- suppressMessages(crwMLE(mov.model=~1, 
                                  err.model = list(x=~1),
                                  data=track, 
                                  Time.name="dtUTC",
                                  initial.state=initial,
                                  fixPar=fixPar,
                                  control=list(maxit=30, trace=0,REPORT=1),
                                  initialSANN=list(maxit=200, trace=0, REPORT=1)))
  # define min and max times in data
  predTime <- seq(ceiling_date(min(track$dtUTC), "min"), 
                  floor_date(max(track$dtUTC), "min"),
                  by=ts)
  # pedict locations
  predObj <- crwPredict(object.crwFit=fit1, 
                        predTime, 
                        return.type='flat')
  # view map
  if (view_map == TRUE){
    crwPredictPlot(predObj, "map")
    title(main=paste("Track", track$id[2], "predicted path"))
  }
  # make df of predicted points
  predObj2 <- split(predObj, predObj$locType) #split predObj
  df <- predObj2[2] # grab predicted ('p') points
  df <- ldply(df, data.frame) # coerce to df
  df <- df[c('id', 'dtUTC', 'mu.x', 'mu.y', 'speed')]
  df <- plyr::rename(df, c('mu.x' = 'lon', 'mu.y' = 'lat'))
  # prepare UTM coordinates matrix
  utmcoor <- SpatialPoints(cbind(df$lon,df$lat), 
                           proj4string=CRS("+proj=utm +zone=56H"))
  longlatcoor <- spTransform(utmcoor,CRS("+proj=longlat"))
  df$lon <- longlatcoor$coords.x1
  df$lat <- longlatcoor$coords.x2
  return(df)
}

#########################
# make predictive paths #
#########################
crawl_predict_all <- function(df, ts, view_map=FALSE){
  ts <- paste(ts,'min')
  coordinates(df) = ~lon+lat
  proj4string(df) <- CRS("+proj=longlat")
  df <- spTransform(df, CRS(paste("+proj=utm +zone=56H +datum=WGS84",
                                  "+ellps=WGS84 +units=m +no_defs")))
  track_ls <- split(df, df$id)
  print(paste('Building predicted paths for',length(track_ls),'tracks..'))
  track_ls <- pblapply(track_ls, function(x) crawl_predict(x, ts=ts, view_map=view_map))
  # merge dataframes
  df <- do.call("rbind", track_ls)
  row.names(df) <- 1:nrow(df)
  
  return(df)
}

#####################
# Split multi-trips #
#####################
split_trips <- function(df, landthreshold=80, maxwindow=600, trim2land=FALSE, view_map=FALSE, 
                        shp="assets/shp/NSW-coastline/edited/polygon/NSW-coastline_WGS84", 
                        return_ocean=FALSE){
  df_ls <- split(df, df$id)
  # scan for large periods on land
  print('Scanning tracks for multiple trips..')
  i <- 1
  multis <- list()
  for (track in df_ls){
    count <- as.integer(table(track$ocean)[1])
    if (count > landthreshold){
      multis <- list.append(multis, i)
    }
    i <- i + 1
  }
  if (length(multis) > 0){
    multi_ids <- lapply(multis, function(x) as.character(df_ls[[x]]$id[1]))
    print(paste(length(multis),'track(s) found with multiple trips:',
                paste(sapply(multi_ids, paste, collapse=":"), collapse=", ")))
  }
  # scan through and identify window
  # only works for multis with 2 trips
  print('Splitting track(s) into distinct trips..')
  multi_df_ls <- df_ls[unlist(multis)]
  i <- 1
  for (track in multi_df_ls){
    row.names(track) <- NULL
    # establish an appropriate window size
    window <- nrow(track[track$ocean == 0,]) %/% 6
    if (window > maxwindow){
      widnow <- maxwindow
    }
    # count land section
    land_count <- 0
    for (j in 1:length(track$ocean)){
      # reset ocean count if ocean point
      if (track$ocean[j] == 1){
        land_count <- 0
      } else if (land_count == window){
        mid_point <- j
        break
      } else {
        land_count <- land_count + 1
      }
    }
    track1 <- track[1:mid_point,]
    track2 <- track[mid_point:nrow(track),]
    row.names(track2) <- NULL
    
    # find and cut track1 tail
    track1_inv <- track1[nrow(track1):1,]
    ocean_count <- 0
    # start counting from end
    for (j in 1:length(track1_inv$ocean)){
      # reset ocean count if land point
      if (track1_inv$ocean[j] == 0){
        ocean_count <- 0
      } else if (ocean_count == window %/% 6){
        land_upper <- j - ocean_count 
        break
      } else {
        ocean_count <- ocean_count + 1
      }
    }
    if (trim2land){
      # trim to land
      track1 <- track1[0:(nrow(track1)-(land_upper-2)),]
    } else {
      # trim to ocean
      track1 <- track1[0:(nrow(track1)-(land_upper-1)),]
    }
    # set split id
    track1$id <- sapply(track1$id, function(x) paste0(x,'a'))
    
    # find and cut track2 head
    ocean_count <- 0
    for (j in 1:length(track2$ocean)){
      # reset ocean count if land point
      if (track2$ocean[j] == 0){
        ocean_count <- 0
      } else if (ocean_count == window %/% 6){
        land_lower <- j - ocean_count
        break
      } else {
        ocean_count <- ocean_count + 1
      }
    }
    if (trim2land){
      track2 <- track2[(land_lower-1):nrow(track2),]
    } else {
      track2 <- track2[land_lower:nrow(track2),]
    }
    # set split id
    track2$id <- sapply(track2$id, function(x) paste0(x,'b'))
    row.names(track2) <- NULL
    
    # check results
    if (view_map == TRUE){
      print('Visual check of results..')
      print('Loading coastal shapefile data..')
      coast <- readOGR(dsn=path.expand(str_sub(shp,1,-(nchar(basename(shp)))-1)), 
                       layer=basename(shp)) # load shapefile
      ocean.land.plot(track1,coast,track1$id[1])
      ocean.land.plot(track2,coast,track2$id[1])
    }
    
    # insert tracks back into dataframe
    print(paste('Merging new tracks',track1$id[1],'and',track2$id[1],'back into primary dataframe..'))
    track <- rbind(track1, track2)
    track$id <- factor(track$id)
    df_ls[[multis[[i]]]] <- track
    i <- i + 1 
  }
  df <- do.call("rbind", df_ls)
  row.names(df) <- NULL
  df$id <- as.factor(df$id)
  if (!return_ocean){
    df$ocean <- NULL
  }
  return(df)
}

#########################
# Clean high velocities #
#########################
# depracated velocity filter
legacy.speed_clean <- function(df, maxSpeed, hist=FALSE){
  # calculate distance between points
  df_ls <- split(df, df$id)
  print('Filtering points with erroneous velocities..')
  pb <- txtProgressBar(min = 0, max = length(df_ls)-1, style = 3)
  i <- 1
  for (df in df_ls){
    # creating lagging columns
    df$latLag <- lag(df$lat)
    df$lonLag <- lag(df$lon)
    df$dtUTCLag <- lag(df$dtUTC)
    # harversine distance
    df$distance <- mapply(function(lon1,lat1,lon2,lat2) distHaversine(c(lon1,lat1),c(lon2,lat2)),
                          df$lon,df$lat,df$lonLag,df$latLag)
    # timedelta
    df$td <- mapply(function(t1, t2) as.numeric(difftime(t2,t1,units='secs')), df$dtUTCLag, df$dtUTC)
    # calculate velocity 
    df$speed <- mapply(function(d, td) d/td, df$distance, df$td)
    # filter bad points
    df <- df[!(df$speed >= maxSpeed & is.na(df$speed) == FALSE),]
    df <- df[c('id','dtUTC','lat','lon','speed')]
    df_ls[[i]] <- df
    setTxtProgressBar(pb, i)
    i <- i + 1
  }
  df <- do.call("rbind", df_ls)
  row.names(df) <- 1:nrow(df)
  if (hist){
    hist(df$speed, col='grey', breaks=40, main="Historgram of penguin velocities (m/s)")
  }
  df <- df[c('id','dtUTC','lat','lon')]
  return(df)
}

#########################
# Clean high velocities #
#########################
speed.clean <- function(df, maxSpeed=20){
  df_ls <- split(df, df$id)
  # process each track
  i <- 1
  track_count <- 0
  for (track in df_ls){
    count_before <- nrow(track)
    # calculate distance, time and speed
    track <- track.distance(track)
    track <- track.time(track)
    track <- track.speed(track)
    # now begin filter process
    filter_incomplete <- TRUE
    n <- 0
    while (filter_incomplete){
      # filter
      track <- track[track$speed < maxSpeed & !is.na(track$speed),] # adding an NA???
      # recaclulate details
      track <- track.distance(track)
      track <- track.time(track)
      track <- track.speed(track)
      # permuate
      n <- n + 1
      # check if filtering complete
      if (nrow(track[track$speed > maxSpeed & !is.na(track$speed),]) == 0){
        filter_incomplete <- FALSE
      }
    }
    # points removed
    count_remove <- count_before-nrow(track)
    if (count_remove > 0){
      track_count <- track_count + 1
      track_portion <- round((count_remove/count_before)*100,0)
      if (track_portion < 25){
        warning = ''
      } else if (track_portion >= 25 && track_portion < 50){
        warning = ' *'
      } else if (track_portion >= 50 && track_portion < 75){
        warning = ' **'
      } else if (track_portion >= 75){
        warning = ' ***'
      }
      print(paste0(track$id[2],': Removed ',count_remove,' points of ',count_before,
                   ' in ',n,' iterations. (',track_portion,'%)',warning))
    }
    # place back in data frame 
    df_ls[[i]] <- track
    i <- i + 1
  }
  # report
  print(paste('Removed points with speeds above',maxSpeed,'m/s from',track_count,'tracks.'))
  # merge dataframe
  print('Merging data frames..')
  df <- do.call('rbind', df_ls)
  # remove dist and reset index
  df$dist <- NULL
  df$time <- NULL
  df$speed <- NULL
  row.names(df) <- NULL
  return(df)
}

######################################
# LEGACY: Read in and save dive data #
######################################
# dir='/Volumes/FieldData/master-data/penguins/tracks/2018/'
# fn <- 'penguin_dives_2018.rds'
# function
legacy.get.penguindive.data <- function(dir, fn, return=FALSE){
  print('Obtaining dive data, this is a slow process.. Plese wait..')
  file_ls <- list.files(path=paste0(dir, 'csv/'), pattern=".csv$")
  # make progress bar 
  pb <- txtProgressBar(min = 0, max = length(file_ls)-1, style = 3)
  setTxtProgressBar(pb, 0)
  df_ls <- list()
  j <- 1
  for (file in file_ls){
    dfp <- read.csv(paste0(paste0(dir, 'csv/'), file))
    dfp <- dfp[ ,c('TagID','Date','Time','Pressure')]
    # drop rows with no pressure
    dfp <- dfp[!is.na(dfp['Pressure']),]
    # convert factors to characters
    for (i in 1:3){
      dfp[,i] <- as.character(dfp[,i])
    }
    # make datetimes
    dfp$dtUTC <- paste(dfp$Date, dfp$Time)
    dfp$dtUTC <- as.POSIXct(strptime(dfp$dtUTC, "%d/%m/%Y %H:%M:%OS", tz="UTC"))
    # drop date and time columns
    dfp[,2:3] <- NULL
    # rename
    colnames(dfp)[1:2] <- c('id', 'pressure')
    # remove "_S[1-9]"
    dfp$id <- sapply(dfp$id, function(x) substr(x, 1, 6))
    # get depth
    dfd <- read.csv(paste0(paste0(dir, 'csv-depthconv/'), file))
    # drop rows with no depth
    dfd <- dfd[!is.na(dfd['Depth']),]
    # merge
    dfp$depth <- dfd$Depth
    # add to data frame list (also rearrange)
    df_ls[[j]] <- dfp[,c(1,3,2,4)]
    setTxtProgressBar(pb, j)
    j <- j + 1
  }
  # merge data.frames
  df <- do.call("rbind", df_ls)
  # reset index
  row.names(df) <- NULL
  # save as RDS
  saveRDS(df, paste0('data/', fn))
  # return if needed
  if (return == TRUE){
    return(df)
  }
}

####################################################
# Load and save penguin KML tracks from data drive #
####################################################
# dir='/Volumes/FieldData/master-data/penguins/tracks/'
load.save.tracks.source <- function(dir, label_L=TRUE, output_fn='penguin_tracks_lachlan_RAW.rds'){
  df <- load_tracks(dir)
  if (label_L){
    df$id <- paste0('L', df$id)
  }
  saveRDS(df, paste0('data/gps/', output_fn))
  print('Data saved..')
}

###########################################################
# Read in and standarise Gemma's legacy GPS tracking data #
###########################################################
process.legacy_tracks <- function(input_dir, output_dir, meta_out){
  # get all the tracks files paths
  path_ls <- list()
  name_ls <- list()
  # read the filenames
  file_ls <- list.files(path=input_dir, pattern=".csv$", recursive=FALSE)
  name_ls <- file_ls
  file_ls <- lapply(file_ls, function(x) paste0(input_dir, x))
  path_ls <- file_ls
  # remove ".csv" in names
  name_ls <- gsub(".csv", "", name_ls)
  # read in all of the data
  print(paste('Reading in', length(name_ls), 'legacy tracks..'))
  track_ls <- list()
  fail_ls <- list()
  i <- 1
  j <- 1
  for (path in path_ls){
    # check header for corrupt files
    header <- readLines(file(path, 'r'), n=1)
    corrupt <- grepl("MONTAGUE", header) | grepl("TIMEVALUE", header) 
    # if file not corrupt try to read
    if (!corrupt){
      print(paste('Loading track', i, 'of', length(path_ls)))
      track_ls[[i]] <- read.csv(path)
    } else {
      fail_ls[[j]] <- name_ls[i]
      j <- j + 1
    }
    i <- i + 1
  }
  print(paste(length(fail_ls),'tracks failed to load..'))
  # remove empty failed loads
  fail_bool <- unlist(lapply(track_ls, is.null))
  # drop failed tracks and names
  track_ls <- track_ls[!fail_bool]
  name_ls <- name_ls[!fail_bool]
  # now there is one messed up track with 2012 and 2015 dates. The 2012 dates are
  # scrambled so we will edit this one manually. "20151002-0517M"
  df <- track_ls[name_ls == "20151002-0517M"][[1]]
  # delete the 2012 rows
  df <- df[substr(df$Date,1,4) != '2012',]
  row.names(df) <- NULL
  track_ls[[which(name_ls == "20151002-0517M")]] <- df
  # process names column into trackid values
  id_ls <- name_ls
  id_ls <- unlist(lapply(id_ls, function(x) str_split(x, '-')[[1]][2]))
  # remove unwanted strings
  id_ls <- unlist(lapply(id_ls, function(x) sub("\\(accel\\)", "", x)))
  id_ls <- unlist(lapply(id_ls, function(x) sub("\\(1\\)", "", x)))
  # fix isolated case
  id_ls <- unlist(lapply(id_ls, function(x)sub("095314 LP 0526", "0526", x)))
  # create list of datetimes of first fixes (for id ordering)
  date_ls <- unlist(lapply(track_ls, function(x) as.character(x$Date[1])))
  time_ls <- unlist(lapply(track_ls, function(x) as.character(x$Time[1])))
  # remove white space
  time_ls <- unlist(lapply(time_ls, function(x) sub(" ", "", x)))
  # Date strings are all formatted differently so need to work out how
  # all are formatted. Just has to be long manual process of deduction unfortunately
  dates_altformatA <- which(nchar(date_ls) == 7 & as.numeric(substr(date_ls, 1, 1)) == 9) # %-m/%d/%y
  dates_altformatB <- which(nchar(date_ls) == 7)[!which(nchar(date_ls) == 7) %in% dates_altformatA][1:4] # %d/%-m/%y
  dates_altformatC <- which(nchar(date_ls) == 7)[!which(nchar(date_ls) == 7) %in% dates_altformatA][5:6] # %m/%-d/%y
  dates_altformatD <- which(nchar(date_ls) == 8 & as.numeric(substr(date_ls, 1, 2)) > 12) # %d/%m/%y
  dates_altformatE <- which(nchar(date_ls) == 8 & as.numeric(substr(date_ls, 4, 5)) > 12) # %m/%d/%y
  dates_altformatF <- which(nchar(date_ls) == 8)[!which(nchar(date_ls) == 8) %in% 
                                                   c(dates_altformatD, dates_altformatE)] # %m/%d/%y
  dates_altformatG <- which(nchar(date_ls) == 9) # %-d/%m/%Y
  dates_altformatH <- which(nchar(date_ls) == 10 & as.numeric(substr(date_ls, 7, 10)) == 2012) # %d/%m/%Y
  dates_altformatI <- which(nchar(date_ls) == 10 & as.numeric(substr(date_ls, 1, 4)) == 2012) # %Y/%d/%m # no longer exists
  dates_altformatJ <- which(nchar(date_ls) == 10)[!which(nchar(date_ls) == 10) %in% 
                                                    c(dates_altformatH, dates_altformatI)] # %Y/%m/%d
  # check we have the right number of dates
  date_format_no <- length(unique(c(dates_altformatA, dates_altformatB, dates_altformatC, 
                             dates_altformatD, dates_altformatE, dates_altformatF, dates_altformatG,
                             dates_altformatH, dates_altformatI, dates_altformatJ)))
  if (date_format_no != length(track_ls)){
    print('ERROR: Something has gone wrong with the date formatting. Check the source code.')
    print(paste(length(track_ls), 'tracks but we have determined', date_format_no, 'date formats...'))
    return()
  }
  # save format strings for later
  formatA = '%m/%d/%y'
  formatB = '%d/%m/%y'
  formatC = '%m/%d/%y'
  formatD = '%d/%m/%y'
  formatE = '%m/%d/%y'
  formatF = '%m/%d/%y'
  formatG = '%d/%m/%Y'
  formatH = '%d/%m/%Y'
  formatI = '%Y/%d/%m'
  formatJ = '%Y/%m/%d'
  # reformat these dates so all is the same
  date_ls[dates_altformatA] <- format(strptime(date_ls[dates_altformatA], formatA), '%Y-%m-%d')
  date_ls[dates_altformatB] <- format(strptime(date_ls[dates_altformatB], formatB), '%Y-%m-%d')
  date_ls[dates_altformatC] <- format(strptime(date_ls[dates_altformatC], formatC), '%Y-%m-%d')
  date_ls[dates_altformatD] <- format(strptime(date_ls[dates_altformatD], formatD), '%Y-%m-%d')
  date_ls[dates_altformatE] <- format(strptime(date_ls[dates_altformatE], formatE), '%Y-%m-%d')
  date_ls[dates_altformatF] <- format(strptime(date_ls[dates_altformatF], formatF), '%Y-%m-%d')
  date_ls[dates_altformatG] <- format(strptime(date_ls[dates_altformatG], formatG), '%Y-%m-%d')
  date_ls[dates_altformatH] <- format(strptime(date_ls[dates_altformatH], formatH), '%Y-%m-%d')
  date_ls[dates_altformatI] <- format(strptime(date_ls[dates_altformatI], formatI), '%Y-%m-%d')
  date_ls[dates_altformatJ] <- format(strptime(date_ls[dates_altformatJ], formatJ), '%Y-%m-%d')
  # finally check dates by cross referencing
  possibly_wrong <- date_ls[(substr(name_ls, 5, 6) != substr(date_ls, 6, 7)) & substr(date_ls, 9, 10) < 13]
  print('MANUAL CHECK!..')
  print('Here are all the dates that were possibly incorrectly formatted:')
  if (length(possibly_wrong) == 0){
    print('All look good!')
  } else {
    print(possibly_wrong)
  }
  # make datetime (time is local)
  dt_ls <- list()
  for (i in 1:length(id_ls)){
    dt <- as.POSIXct(paste(date_ls[i], time_ls[i]), tz="Australia/Sydney")
    dt_ls[[i]] <- dt
  }
  dt_ls <- unlist(dt_ls)
  # so because the burrow ids are the same from year to year we should whack on
  # the year to the start of the id codes so we don't say repeat tracks were 
  # between all the years
  id_ls_yr <- paste0(substr(date_ls, 1, 4), id_ls)
  # now we are going to create factors for burrow pairs to mat deployment codes
  # e.g. 54M and 54F
  pair_factors <- unlist(lapply(id_ls_yr, function(x) substr(x, 1, nchar(x)-1)))
  # report tracks and individuals
  print(paste('There are',length(id_ls),'tracks from',length(levels(as.factor(id_ls_yr))),
              'individuals in',length(levels(as.factor(pair_factors))),'pairs'))
  # get look up table for ids for multiple tracks
  id_counts <- summary(as.factor(id_ls_yr))
  # get the ordered index by timestamps
  order_idx <- order(dt_ls)
  # generate id names based on number of pairs
  new_name_set <- 1:length(levels(as.factor(pair_factors)))
  new_name_set <- str_pad(new_name_set, 3, pad = "0")
  # make the new ordered lists (by first timestamp)
  id_ls_yr <- id_ls_yr[order_idx]
  id_ls <- id_ls[order_idx]
  pair_factors <- pair_factors[order_idx]
  # do our best to determine sex from the id codes
  # first get the last 2 non space characters from the nestbox code
  sex_ls <- unlist(lapply(id_ls_yr, function(x) substr(gsub(" ", "", x), nchar(x)-2, nchar(x))))
  # then scan for capital M and F. If can't find then sex in unknown "U"
  for (i in 1:length(sex_ls)){
    if (grepl('M', sex_ls[i])){
      sex_ls[i] <- 'M'
    } else if (grepl('F', sex_ls[i])){
      sex_ls[i] <- 'F'
    } else {
      sex_ls[i] <- 'U'
    }
  }
  # that sorts out most but we still need to check the unknown remaining sexes
  # sort out Glen Gladys etc.
  males <- c('Glen', 'Axel', 'Ian')
  females <- c('Gladys', 'Rose')
  sex_ls[grepl(paste0(males, collapse = "|"), id_ls)] <- 'M'
  sex_ls[grepl(paste0(females, collapse = "|"), id_ls)] <- 'F'
  # nearly there, now we need to cross reference the metadata. We will cross reference
  # with deployment date and the first GPS fix. I have done this manually because it's
  # a bit messey
  males <- c('0512', 'E4')
  females <- c('0709', '0526')
  sex_ls[grepl(paste0(males, collapse = "|"), id_ls) & sex_ls == 'U'] <- 'M'
  sex_ls[grepl(paste0(females, collapse = "|"), id_ls) & sex_ls == 'U'] <- 'F'
  # now only three unnown sexes and I can live with that
  # now we iterate through the dataset and apply the new names
  new_names <- list()
  name_already_given <- list()
  pair_already_given <- list()
  j <- 0
  for (i in 1:length(id_ls_yr)){
    # if name not already given (not repeat track)
    if (!(pair_factors[i] %in% pair_already_given)){
      j <- j + 1
      new_names[i] <- paste0('G', new_name_set[j], tolower(sex_ls[i]), '01')
      pair_already_given[j] <- pair_factors[i]
      name_already_given[j] <- id_ls_yr[i]
      # else if repeat track
    } else {
      # find index where pair name is given
      k <- match(pair_factors[i], pair_already_given)
      # now check how many times the id has occured
      name_table <- table(id_ls_yr[1:i-1])
      occurance_no <- as.integer(name_table[names(name_table) == id_ls_yr[i]])
      # if first trip set
      if (length(occurance_no) == 0){
        occurance_no <- 0
      }
      # now make the new name
      new_names[i] <- paste0('G', new_name_set[k], tolower(sex_ls[i]), str_pad(occurance_no+1, 2, pad = "0"))
    }
  }
  # unlist the new name codes
  new_names <- unlist(new_names)
  # Now before we go saving these new names we first need to reformat all the datetimes
  # using the formats we determined before
  format_ls <- list()
  format_ls[dates_altformatA] <- formatA
  format_ls[dates_altformatB] <- formatB
  format_ls[dates_altformatC] <- formatC
  format_ls[dates_altformatD] <- formatD
  format_ls[dates_altformatE] <- formatE
  format_ls[dates_altformatF] <- formatF
  format_ls[dates_altformatG] <- formatG
  format_ls[dates_altformatH] <- formatH
  format_ls[dates_altformatI] <- formatI
  format_ls[dates_altformatJ] <- formatJ
  format_ls <- unlist(format_ls)
  # loop through the tracks and reformat
  # note that the process get's slower so could vectorise? nah, who cares.
  # it's fast enough
  i <- 1
  for (df in track_ls){
    print(paste('Reformatting date and time for track',i,'of',length(track_ls)))
    df <- df[complete.cases(df$Latitude), ] # fix bug where NAs throw out process
    df$Time <- unlist(lapply(df$Time, function(x) sub(" ", "", x)))
    df$dt <- as.POSIXct(paste(format(strptime(df$Date, format_ls[i]), '%Y-%m-%d'), df$Time), tz="Australia/Sydney")
    df$dtUTC <- df$dt
    attr(df$dtUTC, "tzone") <- "UTC"
    df[,c('Date', 'Time', 'Altitude', 'Speed', 'Course', 'Type', 'Distance', 'Essential', 'dt')] <- NULL
    names(df)[names(df) == 'Latitude'] <- 'lat'
    names(df)[names(df) == 'Longitude'] <- 'lon'
    track_ls[[i]] <- df
    i <- i + 1
  }
  # now we go through and add the id column and reorder the columns
  # first order the tracks properly
  track_ls <- track_ls[order_idx]
  # now iterate through and process the tracks again
  print('Standardising tracks..')
  i <- 1
  for (df in track_ls){
    df$id <- new_names[i]
    df <- df[,c('id', 'dtUTC', 'lon', 'lat')]
    track_ls[[i]] <- df
    i <- i + 1
  }
  # bind together
  print("Merging the tracks...")
  df <- do.call('rbind', track_ls)
  # last step is to filter our tracks that never left the island
  bad_tracks <- filter.notrack(df, minthreshold_km=5, map.large=FALSE, plot_warnings=FALSE)
  # map those that are being dumped
  tracks.map(df, bad_tracks, title='Tracks to be removed')
  # And create a meta data file so tracks can be linked back to their original files
  print('Making meta data..')
  meta_new_names <- unlist(lapply(track_ls, function(x) x$id[1]))
  meta_orig_id <- id_ls
  meta_files <- unlist(file_ls)[order_idx]
  # make dataframe
  df_meta <- data.frame(id=as.character(meta_new_names),
                        orig_id=id_ls,
                        year=as.factor(substr(id_ls_yr,1,4)),
                        track_file=meta_files)
  # add in raw GPS start and end times
  df_meta$gps_start_raw <- NA
  df_meta$gps_stop_raw <- NA
  for (i in 1:nrow(df_meta)){
    df_meta$gps_start_raw[i] <- min(df$dtUTC[df$id == df_meta$id[i]],na.rm=T)
    df_meta$gps_stop_raw[i] <-max(df$dtUTC[df$id == df_meta$id[i]],na.rm=T)
  }
  df_meta$gps_start_raw <- as.POSIXct(df_meta$gps_start_raw, origin='1970-01-01', tz='UTC')
  df_meta$gps_stop_raw <- as.POSIXct(df_meta$gps_stop_raw, origin='1970-01-01', tz='UTC')
  # save
  write.csv(df_meta, meta_out, row.names = FALSE)
  # finally save data
  # filter bad tracks
  df <- df[!df$id %in% bad_tracks,]
  # Now finally we can save the dataset as a RDS file
  print('Saving to file..')
  print(output_dir)
  saveRDS(df, output_dir)
  print('Done!')
}

##########################
# Get sex factor from id #
##########################
# needs to be updated!!
sex.factor <- function(df, id_sex_idx=5){
  df$sex <- as.factor(str_sub(df$id, id_sex_idx, id_sex_idx))
  return(df)
}

###################
# Read OCL tracks #
###################
read_oclTracks <- function(input_fn){
  # read
  df <- read.csv(input_fn)
  # standarise
  df$id <- df$track_id
  df$lon <- df$Longitude
  df$lat <- df$Latitude
  df$site <- df$site
  df$dtAEST <- as.POSIXct(df$DateTime, tz="Australia/Sydney", format='%FT%T')
  df$dtUTC <- df$dtAEST 
  attr(df$dtUTC, "tzone") <- "UTC"
  # drop
  df <- df[,c('id', 'dtUTC', 'dtAEST', 'lon', 'lat', 'site')]
  return(df)
}

######################################
# Determine if track left the island #
######################################
# Calculates the range of travel to filter out dead tracks
# Please note that this function is very specific to Gemma Carrolls tracks and should
# not be applied to other datasets
filter.notrack <- function(df, minthreshold_km=5, maxthreshold_km=90, plot_warnings=TRUE, 
                           map.large=TRUE, filter_repeat=5, maxSpeed=20){
  df_ls <- split(df, df$id)
  # now go through and check diagonal radius of travel
  warningBig <- list()
  warningSmall <- list()
  smallDist <- list()
  largeDist <- list()
  i <- 1
  small_i <- 1
  large_i <- 1
  # iterate
  for (track in df_ls){
    minlon <- min(track$lon)
    maxlon <- max(track$lon)
    minlat <- min(track$lat)
    maxlat <- max(track$lat)
    # use harversine to calc diagonal distance (returns meters so divde by 1000)
    dist <- distHaversine(c(minlon,minlat),c(maxlon,maxlat))/1000
    if (dist < minthreshold_km){
      warningSmall[small_i] <- i
      smallDist[small_i] <- dist
      small_i <- small_i + 1
    } else if (dist > maxthreshold_km){
      warningBig[large_i] <- i
      largeDist[large_i] <- dist
      large_i <- large_i + 1
    }
    i <- i + 1
  }
  # get suspect tracks
  tracksSmall <- df_ls[unlist(warningSmall)]
  tracksBig <- df_ls[unlist(warningBig)]
  # small_plots
  if (length(tracksSmall) > 0){
    if (plot_warnings){
      small_plot <-tracks.wrap.lineplots.dist(tracksSmall, title='Low distance', dist_ls=smallDist)
    }
    i <- 1
    print('Low distance')
    for (track in tracksSmall){
      print(paste0(track$id[1],': ',as.character(smallDist[i])))
      i <- i + 1
    }
  }
  # large plots
  if (length(tracksBig) > 0){
    if (plot_warnings){
      large_plot <- tracks.wrap.lineplots.dist(tracksBig, title='Large distance', dist_ls=largeDist)
    }
    # filter velocity
    tracksBig_filter <- tracksBig
    for (i in 1:filter_repeat){
      tracksBig_filter <- legacy.speed_clean(do.call('rbind', tracksBig_filter), maxSpeed=maxSpeed) # in m/s
      tracksBig_filter <- split(tracksBig_filter, tracksBig_filter$id)
    }
    # recalculate filter distance
    filterDist <- list()
    i <- 1
    for (track in tracksBig_filter){
      minlon <- min(track$lon)
      maxlon <- max(track$lon)
      minlat <- min(track$lat)
      maxlat <- max(track$lat)
      # use harversine to calc diagonal distance (returns meters so divde by 1000)
      dist <- distHaversine(c(minlon,minlat),c(maxlon,maxlat))/1000
      filterDist[i] <- dist
      i <- i + 1
    }
    # now for any that have greater 2* the max filter remove points beyond that distance
    persistent_error_idx <- which(filterDist > 2*maxthreshold_km)
    persistent_error <- tracksBig_filter[persistent_error_idx]
    if (length(persistent_error) > 0){
      i <- 1
      for (track in persistent_error){
        track$distOrig <- unlist(mapply(distHaversine, 
                                        mapply(c, rep(track$lon[1], length(track$lon)), 
                                               rep(track$lat[1], length(track$lat)), SIMPLIFY = FALSE), 
                                        mapply(c, track$lon, track$lat, SIMPLIFY = FALSE)))
        track$distOrig <- track$distOrig/1000
        track <- track[track$distOrig < 2*maxthreshold_km,]
        track$distOrig <- NULL
        tracksBig_filter[[persistent_error_idx[i]]] <- track
        # new filterDist
        minlon <- min(track$lon)
        maxlon <- max(track$lon)
        minlat <- min(track$lat)
        maxlat <- max(track$lat)
        # use harversine to calc diagonal distance (returns meters so divde by 1000)
        dist <- distHaversine(c(minlon,minlat),c(maxlon,maxlat))/1000
        filterDist[persistent_error_idx[i]] <- dist
      }
    }
    # make filter plot
    if (plot_warnings){
      large_plot_filter <- tracks.wrap.lineplots.dist(tracksBig_filter, 
                                                      title='Large distance error filtered', dist_ls=filterDist)
      # make map
      if (map.large){
        large_filter_map <- tracks.map(c(tracksBig_filter, tracksSmall), c())
      } else {
        large_filter_map <- tracks.map(c(tracksSmall), c())
      }
    }
    i <- 1
    print('Large distance')
    for (track in tracksBig){
      print(paste0(track$id[1],': ',as.character(largeDist[i])))
      i <- i + 1
    }
  }
  if (plot_warnings){
    do.call("grid.arrange", c(list(small_plot, large_plot, large_plot_filter, large_filter_map), ncol=2))
  }
  # Finally get the bad tracks but checking which points are all outside Montague
  mean_lon <- unlist(lapply(tracksSmall, function(x) mean(x$lon)))
  bad_tracks <- names(mean_lon[mean_lon > 150.2])
  return(bad_tracks)
}

###############
# Range Clean #
###############
range.clean <- function(df, lon_orig=150.2269, lat_orig=-36.25201, km=150){
  region = boxmaker(lon_orig=lon_orig, lat_orig=lat_orig, km=km)
  # remove all points outside of the region
  # points in region
  point_bool <- ((df$lon > region[1] & df$lon < region[2]) & (df$lat > region[3] & df$lat < region[4]))
  # filtered
  df <- df[point_bool,]
  return(df)
}

#######################
# duplicate GPS error #
#######################
# This function removes dupicate point GPS errors that 
# apprear after data gaps wheb using igotu tags
duplicate.error.clean <- function(df, gap_threshold=300){
  df_ls <- split(df, df$id)
  # process each track
  i <- 1
  track_count <- 0
  for (track in df_ls){
    count_before <- nrow(track)
    # calcultate distance travelled (distance from last point)
    track <- track.distance(track)
    # get time between points in seconds
    track <- track.time(track)
    # get speed
    track <- track.speed(track)
    # sub in first point
    track$dist[1] <- 1
    track$time[1] <- 1
    # remove points where dist is 0 and time is greater than gap_threshold
    track <- track[!((track$dist == 0) & (track$time > gap_threshold)),]
    # points removed
    count_remove <- count_before-nrow(track)
    if (count_remove > 0){
      track_count <- track_count + 1
      track_portion <- round((count_remove/count_before)*100,0)
      if (track_portion < 25){
        warning = ''
      } else if (track_portion >= 25 && track_portion < 50){
        warning = ' *'
      } else if (track_portion >= 50 && track_portion < 75){
        warning = ' **'
      } else if (track_portion >= 75){
        warning = ' ***'
      }
      print(paste0(track$id[1],': Removed ',count_remove,' points of ',count_before,
                   '. (',track_portion,'%)',warning))
    }
    # place back in data frame 
    df_ls[[i]] <- track
    i <- i + 1
  }
  # report
  print(paste('Removed duplication errors from',track_count,'tracks.'))
  # merge dataframe
  print('Merging data frames..')
  df <- do.call('rbind', df_ls)
  # remove dist and reset index
  df$dist <- NULL
  df$time <- NULL
  row.names(df) <- NULL
  
  return(df)
}

##################################
# Calculate distance for a track #
##################################
track.distance <- function(track){
  # offset coords (call shift from package so as not to mask functions)
  track$lon_lag <- data.table::shift(track$lon)
  track$lat_lag <- data.table::shift(track$lat)
  # calculate distance (haversine)
  track$dist <- mapply(distHaversine,
                       mapply(c, track$lon, track$lat, SIMPLIFY = FALSE),
                       mapply(c, track$lon_lag, track$lat_lag, SIMPLIFY = FALSE))
  # remove offsets
  track[,c('lon_lag', 'lat_lag')] <- NULL
  # return the track
  return(track)
}

################################################
# Calculate displacement for points on a track #
################################################
# intiial coords set to island lighthouse
track.displacement <- function(track, initCoords=c(150.2269, -36.25201)){
  # add initial coords to dataframe
  track$initLon <- initCoords[1]
  track$initLat <- initCoords[2]
  # calculate distance (haversine)
  track$disp <- mapply(distHaversine,
                       mapply(c, track$lon, track$lat, SIMPLIFY = FALSE),
                       mapply(c, track$initLon, track$initLat, SIMPLIFY = FALSE))
  # remove offsets
  track[,c('initLon', 'initLat')] <- NULL
  # return the track
  return(track)
}

#############################################
# Calculate time between points for a track #
#############################################
track.time <- function(track){
  # offset time
  track$dtUTC_lag <- data.table::shift(track$dtUTC)
  # calculate distance (haversine)
  track$time <- mapply(difftime, track$dtUTC, track$dtUTC_lag, units='secs')
  # remove offsets
  track$dtUTC_lag <- NULL
  # return the track
  return(track)
}

################################################
# Calculate time between points for all tracks #
################################################
tracks.time <- function(tracks){
  track.ls <- split(tracks, as.character(tracks$id))
  for (i in 1:length(track.ls)) {
    # offset time
    track.ls[[i]]$dtUTC_lag <- data.table::shift(track.ls[[i]]$dtUTC)
    # calculate distance (haversine)
    track.ls[[i]]$time <- mapply(difftime, track.ls[[i]]$dtUTC, 
                                 track.ls[[i]]$dtUTC_lag, units='secs')
    # remove offsets
    track.ls[[i]]$dtUTC_lag <- NULL
  }
  tracks <- do.call(rbind, track.ls)
  # return the track
  return(tracks)
}

#################################################
# Calculate velcoity between points for a track #
#################################################
# requires distance and time columns
track.speed <- function(df){
  # calc velocity
  df$speed <- df$dist / df$time
  return(df)
}
##################################################
# Track speed without previous values calculated #
##################################################
# wrapper function for speed when distance and time columns not present
track.speed.longcalc <- function(df){
  track <- df
  # calc velocity
  track <- track.distance(track)
  track <- track.time(track)
  track <- track.speed(track)
  df$speed <- track$speed
  return(df)
}

#################################
# Add ocean column for track(s) #
#################################
track.ocean <- function(df, 
                        shp="assets/shp/NSW-coastline/edited/polygon/NSW-coastline_WGS84", 
                        view_map=FALSE){
  # load shapefile data
  print('Loading coastal shapefile data..')
  coast <- readOGR(dsn=path.expand(str_sub(shp,1,-(nchar(basename(shp)))-1)), 
                   layer=basename(shp)) # load shapefile
  proj <- proj4string(coast) # grab projection information
  shp <- list(coast, proj) # pass as list
  # calculate distance (haversine)
  print('Calculating which points are over ocean..')
  ## determine if lat/lon pairs are on land or ocean
  ## ocean = 1, land = 0
  dat <- data.frame(lon = df$lon,
                    lat = df$lat)
  coordinates(dat) <- ~ lon + lat
  # use shp projection
  proj4string(dat) <- shp[[2]]
  # check if points on ocean
  ocean <- over(dat, shp[[1]])
  ocean <- ocean['area']
  ocean[is.na(ocean)] <- 1
  # add to dataset
  df$ocean <- as.factor(ocean$area)
  # filter strange error points
  df <- df[!(df$ocean==76.65),]
  df$ocean <- factor(df$ocean, levels=c(0,1))
  # make plot
  if (view_map){
    ocean.land.plot(df, shp[[1]], 'Ocean Points')
  }
  # return the track
  return(df)
}

########################
# Make local timestamp #
########################
track.dtAEST <- function(df, replace_dtUTC=FALSE, place_after_dtUTC=TRUE){
  dtAEST <- df$dtUTC
  attr(dtAEST, "tzone") <- "Australia/Sydney"
  if (place_after_dtUTC){
    insert_spot <- which(names(df)=="dtUTC")
    df <- data.frame(df[1:insert_spot],dtAEST=dtAEST,df[(insert_spot+1):ncol(df)])
  } else {
    df$dtAEST <- dtAEST
  }
  if (replace_dtUTC){
    df$dtUTC <- NULL
  }
  return(df)
}

###################################
# convert datetime vector to AEST #
###################################
# this changes timezone without changing timestamp string!
# i.e. integer value will change
convert2AEST <- function(dt){
  attr(dt, "tzone") <- "Australia/Sydney"
  return(dt)
}

##############################
# convert AEST vector to UTC #
##############################
# this changes timezone without changing timestamp string!
# i.e. integer value will change
convert2UTC <- function(dt){
  attr(dt, "tzone") <- "UTC"
  return(dt)
}

####################
# Shapefile filter #
####################
# removes all points directly around Branaguba as
# specified by custom shapefile
shapefile.filter <- function(df, shp="assets/shp/baranguba_radius/baranguba_radius"){
  print('Loading shapefile data..')
  zone <- readOGR(dsn=path.expand(str_sub(shp,1,-(nchar(basename(shp)))-1)), 
                   layer=basename(shp)) # load shapefile
  proj <- proj4string(zone) # grab projection information
  ## determine if lat/lon pairs are in the shape file
  ## out = 1, in = 0
  dat <- data.frame(lon = df$lon,
                    lat = df$lat)
  coordinates(dat) <- ~ lon + lat
  # use shp projection
  proj4string(dat) <- proj
  # check if points on ocean
  print('Removing points in exclusion zone..')
  in_zone <- over(dat, zone)
  in_zone <- as.vector(in_zone[['area']])
  in_zone <- !is.na(in_zone)
  # remove values in the exclusion zone
  df <- df[!in_zone,]
  return(df)
}

###########################
# Shapefile tracks delete #
###########################
# Deletes tracks that leave a spaefile region
inside.survey.zone <- function(tracks, threshold=0, UTM=FALSE, 
                               plot.map=FALSE, plot.title=NULL,
                               shp="./assets/shp/survey-coast/survey_coast_straight"){
  df <- tracks
  print('Loading shapefile data..')
  zone <- readOGR(dsn=path.expand(str_sub(shp,1,-(nchar(basename(shp)))-1)), 
                  layer=basename(shp)) # load shapefile
  proj <- proj4string(zone) # grab projection information
  # HACK SOLUTION
  zone@data$id <- 'inside'
  warning('HEY! This function is hacked to work only with ',
  '"./assets/shp/survey-coast/survey_coast_straight". Be careful if ',
  'using this with other data..')
  ## determine if lat/lon pairs are in the shape file
  ## in shape = FALSE, outshape=TRUE
  dat <- data.frame(lon = df$lon,
                    lat = df$lat)
  if (UTM){
    dat <- na.sub.FUN.df(dat, c('lon','lat'), FUN=UTM2lonlat.df)
  }
  
  # check the filtering
  na.idx <- which(!complete.cases(dat[,c('lon','lat')]))
  dat.idx <- which(complete.cases(dat[,c('lon','lat')]))
  # filter NAs
  dat.filt <- dat[dat.idx,]
  # apply the function
  dat.coord <- dat.filt
  coordinates(dat.coord) <- ~ lon + lat
  proj4string(dat.coord) <- proj
  inside <- !is.na(over(dat.coord, zone))
  dat.filt$inZone <- inside
  # sub back into first df
  df$inZone <- NA
  df[dat.idx,'inZone'] <- dat.filt$inZone

  # Split the dataframe
  df.ls <- split(df, df$id)
  # For each track see if it leaves survey area
  df.ls <- lapply(df.ls, function(df) sum(!df$inZone, na.rm=T) > threshold)
  inside <- names(df.ls)[!unlist(df.ls)]
  dropped <- names(df.ls)[!unique(tracks$id) %in% inside]
  message(length(dropped),' tracks dropped..')
  #print(dropped)
  
  # make before map if requested
  if (plot.map){
    tracks.before <- tracks
    N <- length(unique(tracks.before$id))
    custom.limits <- list(xmin=min(tracks.before$lon, na.rm=T),
                          xmax=max(tracks.before$lon, na.rm=T),
                          ymin=min(tracks.before$lat, na.rm=T),
                          ymax=max(tracks.before$lat, na.rm=T))
    m1 <- tracks.map(tracks.before, plot_line=F, legend=F, 
                     title=paste0('Before filtering (N = ',N,')'), zoom=10,
                     custom.limits=custom.limits) +
      geom_polygon(data=zone, aes(x=long, y=lat, group=group), 
                   fill='transparent', col='grey', size=1.8)
  }
  
  # filter
  tracks <- tracks[tracks$id %in% inside,]
  
  # make after map if requested
  if (plot.map){
    N <- length(unique(tracks$id))
    m2 <- tracks.map(tracks, plot_line=F, legend=F, 
                     title=paste0('After filtering (N = ',N,')'), zoom=10,
                     custom.limits=custom.limits) +
      geom_polygon(data=zone, aes(x=long, y=lat, group=group), 
                   fill='transparent', col='grey', size=1.8)
    m <- ggarrange(m1, m2)
    if (!is.null(plot.title)){
      m <- annotate_figure(m, top=plot.title)
    }
    print(m)
      
  }
  
  return(tracks)
}

##########################
# Remove all land points #
##########################
remove.land <- function(df, calc_ocean=FALSE){
  # if needed calculate ocean column
  if (calc_ocean){
    df <- track.ocean(df)
  }
  # return
  df <- df[df$ocean == 1,]
  if (calc_ocean){
    df$ocean <- NULL
  }
  return(df)
}

##########################
# Filter multi day trips #
##########################
tracks.filter.multiday <- function(df, max_days=1){
  i <- 1
  keep_ids <- list()
  before_count <- length(unique(df$id))
  # check if trip is longer than `max_days`
  for (track in split(df, df$id)){
    td <- as.numeric(difftime(track$dtUTC[nrow(track)], track$dtUTC[1], units='days'))
    if (td < max_days){
      keep_ids[[i]] <- track$id[1]
      i <- i + 1
    }
  }
  # filter
  keep_ids <- unlist(keep_ids)
  df <- df[df$id %in% keep_ids,]
  df <- rows.and.levels(df)
  after_count <- length(unique(df$id))
  message(before_count - after_count,' tracks dropped')
  return(df)
}

###################
# Temporal filter #
###################
# time as ISO date string e.g. '2019-10-23'
track.date.filter <- function(df, t_start='1900-01-01', t_end='2100-01-01'){
  t_start <- as.POSIXct(t_start, tz='UTC')
  t_end <- as.POSIXct(t_end, tz='UTC')
  # only keep tracks in this period
  i <- 1
  keeps <- list()
  for (track in split(df, df$id)){
    if (track$dtUTC[1] > t_start & track$dtUTC[nrow(track)] < t_end){
      keeps[[i]] <- track$id[1]
      i <- i + 1
    }
  }
  keeps <- unlist(keeps)
  df <- df[df$id %in% keeps,]
  df$id <- droplevels(df$id)
  return(df)
}

###################
# Tracks selector #
###################
# quick wrapper to select a track easily
track.selector <- function(df, track, reset_index=TRUE, droplevels=TRUE){
  df <- df[df$id %in% track,]
  if (reset_index){row.names(df) <- NULL}
  if (droplevels){df$id <- droplevels(df$id)}
  return(df)
}

######################
# Load Axy dive data #
######################
axy.data.source.2.RDS <- function(input_dir='/Volumes/LP_MstrData/master-data/penguins/tracks/',
                                  tracks_rds='./data/GPS/penguin_tracks_all_CLEAN.rds',
                                  output_fn_dives='./data/dive/Axy/penguins_depth_temp_axy.rds',
                                  output_dir_accel='./data/dive/Axy/accel/',
                                  atmospheric_pressure_fn='./data/BOM/atmospheric_pressure/atom_pressure.csv',
                                  trim2track=TRUE){
  # (1) Load penguin track data and find tracks with multiple trips so 
  # dive data can be split
  tracks_df <- readRDS(tracks_rds)
  # get track id character length
  tracks_id_len <- unlist(lapply(as.character(tracks_df$id), nchar))
  # find those greater than 7
  multi_ids <- unique(tracks_df$id[which(tracks_id_len > 7)])
  # remove Gs (Gemma used CEFAS not Axy tags)
  multi_ids <- factor(multi_ids[substr(multi_ids,1,1) != 'G'])
  # make into dataframe
  multi_ids_timerange <- data.frame(id = multi_ids)
  # get id without a, b string
  multi_ids_timerange$id_short <- substr(multi_ids_timerange$id, 1, 7)
  # start
  multi_ids_timerange$start <- as.POSIXct(unlist(lapply(as.character(multi_ids_timerange$id), 
                                                        function(x) min(tracks_df$dtUTC[tracks_df$id == x]))),
                                          origin='1970-01-01')
  # end
  multi_ids_timerange$end <- as.POSIXct(unlist(lapply(as.character(multi_ids_timerange$id),
                                                      function(x) max(tracks_df$dtUTC[tracks_df$id == x]))),
                                        origin='1970-01-01')
  # now get list of ids that need to be processed in main dive data base
  ids_2_process <- unique(multi_ids_timerange$id_short)
  
  # read in atmospheric pressure data for depth calibration
  atom_pressure <- read.csv(atmospheric_pressure_fn)
  atom_pressure$date <- as.character(atom_pressure$date)
  
  # (2) Load the AXY dive and acceleromtry data
  # get all csv files in directory
  file_ls <- list.files(input_dir, pattern=".csv$", recursive=TRUE)
  # Pressure
  file_ls_pressure <- file_ls[unlist(lapply(file_ls,
                                            function(x) grepl('/csv/', x, fixed=TRUE)))]
  # split the filenames to group multi file data
  file_ls_pressure <- split(file_ls_pressure, substr(basename(file_ls_pressure),1,6))
  # read in all the dataframes for pressure and depth
  df_dive_ls <- list()
  message(paste0('\nReading Axy dive data from "',input_dir,'".'))
  # setup progress bar
  pb <- txtProgressBar(min = 0, max = length(file_ls_pressure)-1, style = 3)
  for (i in 1:length(file_ls_pressure)){
    setTxtProgressBar(pb, i)
    # get filenames
    filename_pressure <- paste0(input_dir,file_ls_pressure[i][[1]])
    # check if multi file and then read in
    if (length(filename_pressure) == 1){
      # read in pressure
      dfp <- read.csv(filename_pressure)
    } else if (length(filename_pressure) > 1){
      dfp <- do.call('rbind', lapply(filename_pressure, read.csv))
    } else {
      message(paste('Error: Zero length files on iteration',i))
      next
    }
    # make datetimes
    dfp$dtUTC <- paste(dfp$Date, dfp$Time)
    dfp$dtUTC <- as.POSIXct(strptime(dfp$dtUTC, "%d/%m/%Y %H:%M:%OS", tz="UTC"))
    # fix id column
    dfp$TagID <- paste0('L',substr(dfp$TagID, 1, 6))
    # order by timestamps (incase multiple files scrambled "_S#")
    dfp <- dfp[order(dfp$dtUTC),]
    # trim the data to the GPS track period
    if (trim2track){
      # use substring due to a and b multi trips in data frame
      dfp <- dfp[dfp$dtUTC > min(tracks_df[substr(tracks_df$id,1,7) == dfp$TagID[1],'dtUTC']) &
                  dfp$dtUTC < max(tracks_df[substr(tracks_df$id,1,7) == dfp$TagID[1],'dtUTC']),]
    }
    if (nrow(dfp)  == 0){
      message(paste('Error: Bad filter for',dfp$TagID[1]))
    }
    # make separate acceleromtry data file as won't be dropping data
    df_accel <- dfp
    # remove uneeded data
    dfp <- dfp[ ,c('TagID','dtUTC','Pressure','Temp....C.')]
    df_accel <- df_accel[ ,c('TagID','dtUTC','X','Y','Z')]
    # drop rows with no pressure/depth
    dfp <- dfp[!is.na(dfp['Pressure']),]
    # rename columns
    colnames(dfp) <- c('id','dtUTC','pressure','temp')
    colnames(df_accel) <- c('id','dtUTC','X','Y','Z')
    # calulcate depth from pressure (skip if no results i.e. L041m01 where sensors did not work)
    if (nrow(dfp) > 0){
      # get local timestamps
      local_t <- format(track.dtAEST(dfp)$dtAEST, '%F')
      atmosphere <- unname(sapply(local_t, 
                                   function(x) atom_pressure[atom_pressure$date == x, 'pressure_mean']))
      dfp$depth <- pressure2depth(dfp$pressure, atmosphere)
      # zero negative depths
      dfp$depth[dfp$depth < 0] <- 0
    } else {
      # hack to add an empty column to an empty dataframe so later rbind is consistent
      depth <- data.frame(depth=1)
      depth <- depth[depth$depth == 2,]
      dfp <- cbind(dfp, depth)
    }
    # flip Y and X axis because I am an idiot and put the tags on upside down for 3 years..
    # At least I was consistent
    df_accel$X <- -df_accel$X
    df_accel$Y <- -df_accel$Y
    # reorder
    dfp <- dfp[,c('id','dtUTC','depth','pressure','temp')]
    
    # if a multi trip process the id names
    if (dfp$id[1] %in% ids_2_process){
      id_check_df <- multi_ids_timerange[as.character(multi_ids_timerange$id_short) == dfp$id[1],]
      id_check_df <- data.frame(id=rep(id_check_df$id,2), dt=c(id_check_df$start, id_check_df$end))
      # find the correct track id
      dfp$id <- unlist(lapply(dfp$dtUTC,
                              function(x) id_check_df$id[which.min(abs(x-id_check_df$dt))]))
      df_accel$id <- unlist(lapply(df_accel$dtUTC,
                                   function(x) id_check_df$id[which.min(abs(x-id_check_df$dt))]))
      df_accel$id <- droplevels(df_accel$id)
      # split accelerometry data and save
      df_accel <- split(df_accel, df_accel$id)
      # now trim dfp and df_accel to split track ids
      if (trim2track){
        dfp$id <- droplevels(dfp$id)
        dfp <- split(dfp, dfp$id)
        # split the multi trips
        for (j in 1:length(dfp)){
          df <- dfp[[j]]
          dfp[[j]] <- df[df$dtUTC > min(tracks_df[tracks_df$id == as.character(df$id[1]),'dtUTC']) &
                         df$dtUTC < max(tracks_df[tracks_df$id == as.character(df$id[1]),'dtUTC']),]
          df <- df_accel[[j]]
          df_accel[[j]] <- df[df$dtUTC > min(tracks_df[tracks_df$id == as.character(df$id[1]),'dtUTC']) &
                              df$dtUTC < max(tracks_df[tracks_df$id == as.character(df$id[1]),'dtUTC']),]
        }
        # rbind back dfp
        dfp <- do.call('rbind', dfp)
      }
      for (df in df_accel){
        row.names(df) <- NULL
        saveRDS(droplevels(df), paste0(output_dir_accel,df$id[1],'.rds'))
      }
    } else {
      row.names(df_accel) <- NULL
      saveRDS(droplevels(df_accel), paste0(output_dir_accel,df_accel$id[1],'.rds'))
    }
    # add to list and save accel to file
    df_dive_ls[[i]] <- dfp
  }
  # close progress bar
  close(pb)
  # combine
  dives_df <- do.call("rbind", df_dive_ls)
  row.names(dives_df) <- NULL
  # save dive data
  saveRDS(droplevels(dives_df), output_fn_dives)
  # Final hack step add on at end (to avoid re-running everything again!)
  # Subtract minimum value from dives to correct penguin data that is too deep (never goes to surface).
  # Usually just takes 14-17 cm off the dive data
  dives_df <- readRDS(output_fn_dives)
  dives_df <- split(dives_df, dives_df$id)
  for (i in 1:length(dives_df)){
    dives_df[[i]]$depth <- dives_df[[i]]$depth - min(dives_df[[i]]$depth)
  }
  dives_df <- do.call(rbind, dives_df)
  saveRDS(droplevels(dives_df), output_fn_dives)
}

###########################################
# match CEFAS dive data to GPS track data #
###########################################
# creates a linking file of GPS tracks to CEFAS filenames
match.cefas.data.2.gps <- function(input_dir='/Volumes/LP_MstrData/master-data/gemma/accelerometry/',
                                   output_fn='./data/misc/tracks_CEFAS_link.csv',
                                   verbose=FALSE){
  message("Matching penguin GPS tracks to correct CEFAS files..")
  # get all csv files in directory
  file_ls <- list.files(input_dir, pattern=".CSV$", recursive=TRUE)
  # get file basenames and extract data
  files_df <- data.frame(filename=basename(file_ls), filepath=file_ls)
  # read in gemma deployment meta data
  gemma_meta <- read.csv('./data/gps/penguin_tracks_gemma_RAW_meta.csv')
  gemma_meta_pre2105 <- gemma_meta[gemma_meta$year < 2015,]
  gemma_meta <- gemma_meta[gemma_meta$year >= 2015,]
  gemma_meta$track_file <- as.character(gemma_meta$track_file) 
  gemma_meta$filenames <- unlist(lapply(gemma_meta$track_file, function(x) unlist(strsplit(x, '/'))[[8]]))
  gemma_meta$filenames <- file_path_sans_ext(gemma_meta$filenames)
  # work out new id codes from meta data conversion
  files_df$orig_id <- unlist(lapply(strsplit(file_path_sans_ext(files_df$filename), '_'), function(x) x[3]))
  files_df$date <- unlist(lapply(strsplit(file_path_sans_ext(files_df$filename), '_'), function(x) x[2]))
  files_df$date <- as.POSIXct(strptime(files_df$date, "%d-%m-%y", tz="UTC"))
  # work out codes which match gemmas meta data
  files_df$meta_code <- paste0(as.character(year(files_df$date)),
                               as.character(str_pad(month(files_df$date),2,pad='0')),
                               as.character(str_pad(day(files_df$date),2,pad='0')),'-',
                               files_df$orig_id)
  gemma_meta$meta_code <- unlist(lapply(strsplit(gemma_meta$filenames, '-'), function(x) x[1]))
  gemma_meta$meta_code <- paste0(gemma_meta$meta_code,'-',gemma_meta$orig_id)

  # get accel start and stop times
  message('Extracting CEFAS metadata..')
  files_df$start <- NA
  files_df$stop <- NA
  pb <- txtProgressBar(min = 0, max = nrow(files_df), style = 3)
  for (i in 1:nrow(files_df)){
    setTxtProgressBar(pb, i)
    dt_range <- get.CEFAS.date.source(paste0(input_dir,files_df$filepath[i]))
    files_df$start[i] <- dt_range[[1]]
    files_df$stop[i] <- dt_range[[2]]
  }
  close(pb)
  files_df$accel_start <- as.POSIXct(files_df$start, origin='1970-01-01', tz='UTC')
  files_df$accel_stop <- as.POSIXct(files_df$stop, origin='1970-01-01', tz='UTC')
  files_df[,c('start','stop')] <- NULL
  files_df$Accelerometer <- sapply(files_df$filename, function(x) strsplit(as.character(x),'_')[[1]][1])
  
  # get start and end times for all gps tracks
  gemma_meta$start <- NA
  gemma_meta$stop <- NA
  gps_track_data <- readRDS('./data/gps/penguin_tracks_all_CLEAN.rds')
  bad_ids <- list()
  j <- 1
  message('Loading GPS tracks metadata..')
  for (i in 1:nrow(gemma_meta)){
    # check if track exists and if not store to erase from dataframes 
    # (for bad tracks that were filtered in legacy clean)
    if (!as.character(gemma_meta$id[i]) %in% as.character(gps_track_data$id)){
      bad_ids[[j]] <- as.character(gemma_meta$id[i])
      j <- j + 1
    } else {
    gemma_meta$start[i] <- min(gps_track_data$dtUTC[as.character(gps_track_data$id) == 
                                                      as.character(gemma_meta$id[i])])
    gemma_meta$stop[i] <- max(gps_track_data$dtUTC[as.character(gps_track_data$id) == 
                                                     as.character(gemma_meta$id[i])])
    }
  }
  gemma_meta$gps_start <- as.POSIXct(gemma_meta$start, origin='1970-01-01', tz='UTC')
  gemma_meta$gps_stop <- as.POSIXct(gemma_meta$stop, origin='1970-01-01', tz='UTC')
  gemma_meta[,c('start','stop')] <- NULL
  bad_ids <- unlist(bad_ids)
  
  ###### DO THE MATCHING ######
  message('Matching files..')
  gemma_meta$CEFAS_fn <- NA
  gemma_meta$error <- NA
  # Reverse approach
  for (i in 1:nrow(gemma_meta)){
    # first skip bad ids
    if (gemma_meta$id[i] %in% bad_ids){
      if (verbose){
        message('')
        message(i)
        message(paste('Skipping previously filtered track',gemma_meta$id[i]))
      }
      gemma_meta$error[i] <- 'bad_track'
      next
    }
    # filter first by id
    df <- files_df[files_df$orig_id == gemma_meta$orig_id[i],]
    if (nrow(df) < 1){
      if (verbose){
        message('')
        message(i)
        message(paste('Could not find CEFAS data for track', gemma_meta$id[i]))
        message('Error: No matches')
      }
      gemma_meta$error[i] <- 'no_match'
      next
    }
    # filter by year
    df <- df[as.character(year(df$date)) == as.character(gemma_meta$year[i]),]
    if (nrow(df) < 1){
      if (verbose){
        message('')
        message(i)
        message(paste('Could not find CEFAS data for track', gemma_meta$id[i]))
        message('Error: No matches')
      }
      gemma_meta$error[i] <- 'no_match'
      next
    }
    # check if track time range falls into CEFAS time range
    df$time_check <- gemma_meta$gps_start[i] > df$accel_start & 
                     gemma_meta$gps_stop[i] < df$accel_stop
    df <- df[df$time_check,]
    if (nrow(df) < 1){
      if (verbose){
        message('')
        message(i)
        message(paste('Could not find CEFAS data for track', gemma_meta$id[i]))
        message('Error: No matches')
      }
      gemma_meta$error[i] <- 'no_match'
      next
    }
    # if only 1 match we're good!
    if (nrow(df) == 1){
      gemma_meta$CEFAS_fn[i] <- as.character(df$filename[1])
    } else {
      # multi matches :(
      if (verbose){
        message('')
        message(i)
        message(paste('Could not find CEFAS data for track', gemma_meta$id[i]))
        message('Error: Too many matches')
      }
      gemma_meta$error[i] <- paste0('multi_match',nrow(df))
      if (nrow(df) != length(unique(df$Accelerometer))){
        if (verbose){
          message('CEFAS tag has multiple records of simeltaneous deployments (impossible?)')
        }
        gemma_meta$error[i] <- 'CEFAS_duplication'
      }
      if (verbose){
        message('Track')
        print(gemma_meta[i,])
        message('CEFAS Matches')
        print(df)
        message('Deployment Matches')
        print(dply_df)
      }
    }
  }
  # report status of the matches
  message(paste0('CEFAS acceleromtry files found for ',sum(!is.na(gemma_meta$CEFAS_fn)),' of ',
                 nrow(gemma_meta), ' (',round((sum(!is.na(gemma_meta$CEFAS_fn))/nrow(gemma_meta))*100,0),
                 '%) GPS tracks.'))
  if (verbose){
    message(paste0('Failed to find acceleromtry files for ',sum(is.na(gemma_meta$CEFAS_fn)),' GPS tracks'))
    message(paste0('Of these ',sum(is.na(gemma_meta$CEFAS_fn)),' tracks..'))
    message(paste0('\t',sum(gemma_meta$error == 'bad_track',na.rm=T),' were previously filtered tracks'))
    message(paste0('\t',sum(gemma_meta$error == 'no_match',na.rm=T),' had no matches found'))
    message('GPS tracks with no matching CEFAS files found:')
    print(gemma_meta[gemma_meta$error == 'no_match' & !is.na(gemma_meta$error),c('id', 'meta_code')])
  }
  write.csv(gemma_meta[,c('id','meta_code','CEFAS_fn')], output_fn,
            row.names=FALSE)
}

##################################
# Get CEFAS Start and Stop times #
##################################
get.CEFAS.date.source <- function(file_path, UTC=TRUE){
  # load data
  CEFAS <- readLines(file_path)
  ## grab primary data block (only one)
  idx1 <- grep("Data Block", CEFAS)
  idx2 <- grep("Fast Log Data Follows", CEFAS)
  data_block <- CEFAS[(idx1 + 7):(idx2 - 1)] # grab data block
  data_block <- data_block[data_block != ""] # remove empty rows
  # turn data block into dt list
  dt_ls <- unlist(lapply(data_block[c(1,length(data_block))], function(x) strsplit(x, ',')[[1]][1]))
  # get the start and end times
  t_start <- as.POSIXct(dt_ls[1], "%d/%m/%y %H:%M:%OS", tz='Australia/Sydney')
  t_stop <- as.POSIXct(dt_ls[2], "%d/%m/%y %H:%M:%OS", tz='Australia/Sydney')
  if (UTC){
    attr(t_start, "tzone") <- "UTC"
    attr(t_stop, "tzone") <- "UTC"
  }
  return(list(t_start, t_stop))
}

######################
# Extract CEFAS data #
######################
# Note: unlike axy reader this function does not split accel files for multi
# tracks (e.g. a and b etc..) because there are no multi trip files with accel data
CEFAS.data.source.2.RDS <- function(input_dir='/Volumes/LP_MstrData/master-data/gemma/accelerometry/',
                                    link_file='./data/misc/tracks_CEFAS_link.csv',
                                    output_dir='./data/dive/CEFAS/',
                                    tracks_rds='./data/GPS/penguin_tracks_all_CLEAN.rds',
                                    DST_fn='./data/misc/dst_dates.csv',
                                    trim2track=TRUE,
                                    bad_tracks=c('G126f06','G126f07')){
  message('IF RERUNNING THIS PROCESS REMEMBER TO UPDATE CODE TO USE MORE ROBUST DIVE CALC METHOD')
  message('Temporary soltutions at end should be temporary!')
  # Load penguin track data and find tracks with multiple trips so 
  # dive data can be split and trips trimmed
  tracks_df <- readRDS(tracks_rds)
  # remove all before 2015 (no accel data)
  tracks_df <- tracks_df[year(tracks_df$dtUTC) >= 2015,]
  # get track id character length
  tracks_id_len <- unlist(lapply(as.character(tracks_df$id), nchar))
  # find those greater than 7
  multi_ids <- unique(tracks_df$id[which(tracks_id_len > 7)])
  # remove Ls (use Axy not CEFAS tags)
  multi_ids <- factor(multi_ids[substr(multi_ids,1,1) != 'L'])
  # make into dataframe
  multi_ids_timerange <- data.frame(id = multi_ids)
  # get id without a, b string
  multi_ids_timerange$id_short <- substr(multi_ids_timerange$id, 1, 7)
  # start
  multi_ids_timerange$start <- as.POSIXct(unlist(lapply(as.character(multi_ids_timerange$id), 
                                                        function(x) min(tracks_df$dtUTC[tracks_df$id == x]))),
                                          origin='1970-01-01')
  # end
  multi_ids_timerange$end <- as.POSIXct(unlist(lapply(as.character(multi_ids_timerange$id),
                                                      function(x) max(tracks_df$dtUTC[tracks_df$id == x]))),
                                        origin='1970-01-01')
  # now get list of ids that need to be processed in main dive data base
  ids_2_process <- unique(multi_ids_timerange$id_short)
  
  # load DST data
  # NOTE: these are loaded in as UTC but they are actually Australia/Sydney
  DST <- read.csv(DST_fn)
  DST$DST_end <- as.POSIXct(DST$DST_end, tz='UTC')
  DST$DST_start <- as.POSIXct(DST$DST_start, tz='UTC')
  
  # load in the link file
  link_df <- read.csv(link_file)
  # read in all files in CEFAS directory
  file_ls <- list.files(input_dir, pattern=".CSV$", recursive=TRUE)
  # get file basenames and extract data
  files_df <- data.frame(filename=basename(file_ls), filepath=file_ls)
  # add the file path to the linking file
  link_df$filepath <- unlist(lapply(link_df$CEFAS_fn, function(fn) 
    ifelse(is.na(fn), NA, as.character(files_df$filepath[as.character(files_df$filename) == fn]))))
  # how many acceleromtry files don't have assoicated tracks?
  no_tracks <- files_df[!files_df$filename %in% link_df$CEFAS_fn,]
  if (nrow(no_tracks) > 0){
    warning(paste0(nrow(no_tracks),' CEFAS files not read as no associated GPS track was found.'))
  }
  # only consider complete cases for link_df
  link_df <- link_df[complete.cases(link_df),]
  row.names(link_df) <-  NULL
  # confirm no multi trips
  if (sum(ids_2_process %in% link_df$id) >= 1){
    message('Error: A multitrip exists but this function is not able to process it')
    message('Function needs to be recoded.')
    return('Error')
  }
  # for each track read in the CEFAS data and save the results
  message('Reading CEFAS data files..')
  pb <- txtProgressBar(min = 0, max = nrow(link_df), style = 3)
  for (i in 1:nrow(link_df)){
    df <- CEFAS.reader(paste0(input_dir,link_df$filepath[i]), DST)
    df$id <- link_df$id[i]
    if (df$id[1] %in% bad_tracks){
      setTxtProgressBar(pb, i)
      next()
    }
    df <- df[,c('id','dive_id','dtUTC','depth','temp','X','Y','Z')]
    # trim to GPS track
    if (trim2track){
      df <- df[df$dtUTC > min(tracks_df[as.character(tracks_df$id) == df$id[1],'dtUTC']) &
               df$dtUTC < max(tracks_df[as.character(tracks_df$id) == df$id[1],'dtUTC']),]
    }
    #saveRDS(df, paste0(output_dir,link_df$id[i],'.rds'))
    # Final hack step add on at end
    # Subtract minimum value from dives to correct penguin data that is too deep (never goes to surface).
    # Usually just takes 14-17 cm off the dive data
    #df$depth <- df$depth - min(df$depth)
    # UPDATE - Have stopped this because it mucks up dive calculations which
    # are predone in CEFAS
    saveRDS(df, paste0(output_dir,link_df$id[i],'.rds'))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  # some more final hackS so this does not need to be run again
  # Fix the weird data point in G103m01
  fix_df <- readRDS(paste0(output_dir,'G103m01','.rds'))
  fix_df[fix_df$dtUTC == as.POSIXct("2015-09-30 01:24:36.967", tz='UTC'),'depth'] <- 0.44
  # G103m01 was stuck in fast log mode so dives were not properly calculated. 
  # apply axy dive detection to the data
  # extract signal as a vector as this is much faster
  depth <- fix_df$depth
  dive_id_vec <- fix_df$dive_id
  dive_state <- FALSE
  dive_idx <- 0
  prev_depth <- 0
  # find crosses with for loop
  for (i in 1:length(depth)){
    # check for cross over
    if (!dive_state & depth[i] > 2 & prev_depth < 2){
      dive_state <- TRUE
      dive_idx <- dive_idx + 1
    }
    if (dive_state & depth[i] < 1.5 & prev_depth > 1.5){
      dive_state <- FALSE
    }
    dive_id_vec[i] <- ifelse(dive_state, dive_idx, 0)
    prev_depth <- depth[i]
  }
  fix_df$dive_id <- dive_id_vec
  # save RDS again
  saveRDS(fix_df, paste0(output_dir,'G103m01','.rds'))
  # G127m08 incorrectly moved to fast log mode
  fix_df <- readRDS(paste0(output_dir,'G127m08','.rds'))
  # dive ids 280 and 284 are corrupt
  # this is a temporary solution
  fix_df[fix_df$dive_id %in% c(280, 284),'dive_id'] <- 0
  saveRDS(fix_df, paste0(output_dir,'G127m08','.rds'))
}

################
# CEFAS Reader #
################
# Read a CEFAS file
# Note: Due to period od study DST corrections only considers DST start, not end. 
CEFAS.reader <- function(filepath, DST=NA, DST_fn='./data/misc/dst_dates.csv'){
  # if DST not prvided load in from source
  if (is.na(DST)){
    DST <- read.csv(DST_fn)
  }
  # read the file
  CEFAS <- readLines(filepath)
  # extract primary data block
  idx1 <- grep("Data Block", CEFAS)
  idx2 <- grep("Fast Log Data Follows", CEFAS)
  data_block <- CEFAS[(idx1 + 7):(idx2 - 1)]
  data_block <- data_block[data_block != ""] # remove empty rows
  # make into a dataframe
  data_block <- as.data.frame(do.call('rbind',strsplit(data_block, ',')))
  # set columns
  colnames(data_block) <- c('dtAEST','depth','temp','X','Y','Z')
  # make correct data types
  data_block[,2:6] <- apply(data_block[,2:6], 2, function(x) as.numeric(as.character(x)))
  # set dive id
  data_block$dive_id <- 0
  # extract fast log blocks
  ls_idx <- grep("Fast Log", CEFAS) # make list of idexes
  ls_idx <- ls_idx[2:(length(ls_idx))] # remove first line which is not a fast log
  fastlog_ls <- list() # make empty list to store blocks
  # place each fast log block into a list based on index
  for (i in 1:(length(ls_idx)-1)){
    fastlog_ls[i] <- as.data.frame(CEFAS[((ls_idx[i]):(ls_idx[(i+1)]-1))])
  }
  # add last fast log block
  fastlog_ls[length(ls_idx)] <- as.data.frame(
    CEFAS[ls_idx[length(ls_idx)]:length(CEFAS)])
  # if "Wet Table" in log then delete
  if ('Table of Wet Times follows' %in% fastlog_ls[[length(ls_idx)]]){
    wet_table_idx <- grep('Table of Wet Times follows', fastlog_ls[[length(ls_idx)]], fixed=TRUE)
    fastlog_ls[[length(ls_idx)]] <- fastlog_ls[[length(ls_idx)]][1:(wet_table_idx-1)]
  }
  # if data points avaliable = 0 drop fast log
  # else process int data frame
  for (i in 1:length(fastlog_ls)){
    if ('Data points available = 0' %in% fastlog_ls[[i]]){
      fastlog_ls[[i]] <- NA
    } else {
      # drop empty lines
      fastlog_ls[[i]] <- fastlog_ls[[i]][fastlog_ls[[i]] != ""]
      # cut to data
      fastlog_ls[[i]] <- fastlog_ls[[i]][7:length(fastlog_ls[[i]])]
      # turn into data frame
      fastlog_ls[[i]] <- as.data.frame(do.call('rbind',strsplit(as.character(fastlog_ls[[i]]), ',')))
      # set columns
      colnames(fastlog_ls[[i]]) <- c('dtAEST','depth','temp','X','Y','Z')
      # make correct data types
      fastlog_ls[[i]][,2:6] <- apply(fastlog_ls[[i]][,2:6], 2, function(x) as.numeric(as.character(x)))
      # set dive id
      fastlog_ls[[i]]$dive_id <- i
    }
  }
  # remove NAs
  fastlog_ls <- fastlog_ls[!is.na(fastlog_ls)]
  # bind with data block
  df <- do.call('rbind', fastlog_ls)
  df <- rbind(data_block, df)
  
  # load in as UTC so we can compensate for DST
  # Note: this is a fix for CEFAS tags not doing DST compensation
  df$FAKEUTC <- as.POSIXct(df$dtAEST, "%d/%m/%y %H:%M:%OS", tz='UTC')
  # DST compensation
  # get most recent clock set time
  clock_set <- CEFAS[grep('Clock Set', CEFAS)]
  clock_set <- max(as.POSIXct(sapply(clock_set, function(s) strsplit(s, ',')[[1]][2]),
                              "%d/%m/%y %H:%M:%OS", tz='UTC'))
  # get last timestamp in records
  last_time <- max(df$FAKEUTC, na.rm=TRUE)
  # get correct DST cross record
  dst_cross <- DST[DST$year == year(last_time), "DST_start"]
  # do correction if needed
  if (clock_set < dst_cross & last_time > dst_cross){
    # add hour to correct
    df[df$FAKEUTC >= dst_cross,'FAKEUTC'] <- df[df$FAKEUTC >= dst_cross,'FAKEUTC'] + hours(1)
  }
  # convert to characters
  df$dtAEST <- as.character(format(df$FAKEUTC, "%Y-%m-%d %H:%M:%OS3"))
  # read in the real time
  df$dtAEST <- as.POSIXct(df$dtAEST, "%Y-%m-%d %H:%M:%OS", tz='Australia/Sydney')
  # drop fake time
  df$FAKEUTC <- NULL
  # create dtUTC and drop local timestamp
  df$dtUTC <- df$dtAEST
  attr(df$dtUTC, "tzone") <- "UTC"
  df$dtAEST <- NULL
  # reorganise
  df <- df[,c('dive_id','dtUTC','depth','temp','X','Y','Z')]
  # order by dtUTC
  df <- df[order(df$dtUTC),]
  # remove duplicates timestamps (due to fast and slow log)
  df <- df[!duplicated(df$dtUTC, fromLast=TRUE),]
  # reset row names
  row.names(df) <- NULL
  # DONE!
  return(df)
}

################################
# add dive variables to tracks #
################################
# for each track point find nearest depth
# add depth, temp or both
# tracks with no dive data are dropped
tracks.attach.diveVars <- function(tracks, vars=c('depth','temp'), 
                                   method=c('segment_stats','segment_mean','join'), percentile=85,
                                   zero_neg_depths=TRUE, filter_nodive_tracks=TRUE,
                                   axy_dive_fn='./data/dive/Axy/penguins_depth_temp_axy.rds',
                                   CEFAS_dive_dir='./data/dive/CEFAS/'){
  # set method argument
  method <- match.arg(method)
  message('Attaching dive data to GPS tracks using method: "',method,'"')
  # split tracks into Axy and CEFAS
  L_tracks <- tracks[substr(tracks$id,1,1) == 'L',]
  G_tracks <- tracks[substr(tracks$id,1,1) == 'G',]
  # if using segment stats extend variable names
  if (method == 'segment_stats'){
    # extend vars to stats
    vars_stats <- c(paste0(vars,'_mean'),
                    paste0(vars,'_min'),
                    paste0(vars,'_max'),
                    paste0(vars,'_percen',percentile))
  }
  # AXY
  if (nrow(L_tracks) > 0){
    message('Matching Axy dive data..')
    # get Axy dive times
    AXY_dive <- readRDS(axy_dive_fn)
    # zero negative depths
    if (zero_neg_depths){
      AXY_dive$depth[AXY_dive$depth < 0 & !is.na(AXY_dive$depth)] <- 0
    }
    if (method == 'join'){
      # left join data by matching id and dtUTC
      L_tracks[,vars] <- left_join(L_tracks, AXY_dive, by=c('id','dtUTC'))[,vars]
    } else if (method == 'segment_mean'){
      # calulate mean dive period between track segments (back looking)
      L_tracks_ls <- split(L_tracks, droplevels(L_tracks$id))
      # process each
      pb <- txtProgressBar(min=0, max=length(L_tracks_ls), style = 3)
      for (i in 1:length(L_tracks_ls)){
        # calculate lag
        L_tracks_ls[[i]]$dtUTC_lag <- data.table::shift(L_tracks_ls[[i]]$dtUTC)
        # get dive data frame for track
        track_dive <- AXY_dive[AXY_dive$id == L_tracks_ls[[i]]$id[1],]
        # calulate mean variables for each segment
        for (j in 1:nrow(L_tracks_ls[[i]])){
          row <- L_tracks_ls[[i]][j,]
          var_cols <- track_dive[track_dive$dtUTC > row['dtUTC_lag'] &
                                 track_dive$dtUTC < row['dtUTC'],][,vars]
          L_tracks_ls[[i]][j,vars] <- lapply(var_cols, mean, na.rm=TRUE)
        }
        # first value get 2nd value
        L_tracks_ls[[i]][1,vars] <- L_tracks_ls[[i]][2,][vars]
        setTxtProgressBar(pb, i)
      }
      close(pb)
      # merge data
      L_tracks <- do.call('rbind', L_tracks_ls)
      L_tracks$dtUTC_lag <- NULL
    } else if (method == 'segment_stats'){
      # calulate mean dive period between track segments (back looking)
      L_tracks_ls <- split(L_tracks, droplevels(L_tracks$id))
      # process each
      pb <- txtProgressBar(min=0, max=length(L_tracks_ls), style = 3)
      for (i in 1:length(L_tracks_ls)){
        # calculate lag
        L_tracks_ls[[i]]$dtUTC_lag <- data.table::shift(L_tracks_ls[[i]]$dtUTC)
        # get dive data frame for track
        track_dive <- AXY_dive[AXY_dive$id == L_tracks_ls[[i]]$id[1],]
        # calulate mean variables for each segment
        for (j in 1:nrow(L_tracks_ls[[i]])){
          row <- L_tracks_ls[[i]][j,]
          var_cols <- suppressWarnings(track_dive[track_dive$dtUTC > row['dtUTC_lag'] &
                                   track_dive$dtUTC < row['dtUTC'],][,vars])
          # make sure not empty (numeric(0) or all NAs)
          if (length(var_cols) == 0 | sum(is.na(var_cols)) == length(var_cols)){
            L_tracks_ls[[i]][j,vars_stats] <- c(NA,NA,NA,NA)
          } else {
            # attach to dataframe
            L_tracks_ls[[i]][j,vars_stats] <- as.numeric(
              c(lapply(var_cols, mean, na.rm=TRUE),
                lapply(var_cols, min, na.rm=TRUE),
                lapply(var_cols, max, na.rm=TRUE),
                lapply(var_cols, quantile, probs=percentile/100, na.rm=TRUE)))
          }
        }
        # first value get 2nd value
        L_tracks_ls[[i]][1,vars_stats] <- L_tracks_ls[[i]][2,][vars_stats]
        setTxtProgressBar(pb, i)
      }
      close(pb)
      # merge data
      L_tracks <- do.call('rbind', L_tracks_ls)
      L_tracks$dtUTC_lag <- NULL
    }
  }
  
  # CEFAS
  if (nrow(G_tracks) > 0){
    # split tracks
    G_tracks_ls <- split(G_tracks, droplevels(G_tracks$id))
    # loop and read in accel data
    message('Matching CEFAS dive data to each accel file..')
    pb <- txtProgressBar(min=0, max=length(G_tracks_ls), style = 3)
    for (i in 1:length(G_tracks_ls)){
      # read file and if error drop track
      accel <- tryCatch(
        readRDS(paste0(CEFAS_dive_dir,G_tracks_ls[[i]]$id[1],'.rds')),
        error=function(e) e)
      if(inherits(accel, "error")){
        G_tracks_ls[[i]] <- NA
        setTxtProgressBar(pb, i)
        next
      }
      # zero negative depths
      if (zero_neg_depths){
        accel$depth[accel$depth < 0 & !is.na(accel$depth)] <- 0
      }
      if (method == 'join'){
        matches <- left_join(G_tracks_ls[[i]], accel, by=c('id','dtUTC'))
        # remove dups (due to zero block and fast logs duplicating data)
        G_tracks_ls[[i]][,vars] <- matches[!duplicated(matches$dtUTC),vars]
      } else if (method == 'segment_mean'){
        # calculate lag
        G_tracks_ls[[i]]$dtUTC_lag <- data.table::shift(G_tracks_ls[[i]]$dtUTC)
        for (j in 1:nrow(G_tracks_ls[[i]])){
          row <- G_tracks_ls[[i]][j,]
          var_cols <- accel[accel$dtUTC > row['dtUTC_lag'] &
                            accel$dtUTC < row['dtUTC'],][,vars]
          G_tracks_ls[[i]][j,vars] <- lapply(var_cols, mean, na.rm=TRUE)
        }
        # first NULL row gets second
        G_tracks_ls[[i]][1,vars] <- G_tracks_ls[[i]][2,][vars]
      } else if (method == 'segment_stats'){
        # calculate lag
        G_tracks_ls[[i]]$dtUTC_lag <- data.table::shift(G_tracks_ls[[i]]$dtUTC)
        for (j in 1:nrow(G_tracks_ls[[i]])){
          row <- G_tracks_ls[[i]][j,]
          var_cols <- accel[accel$dtUTC > row['dtUTC_lag'] &
                              accel$dtUTC < row['dtUTC'],][,vars]
          # make sure not empty (numeric(0) or all NAs)
          if (length(var_cols) == 0 | sum(is.na(var_cols)) == length(var_cols)){
            G_tracks_ls[[i]][j,vars_stats] <- c(NA,NA,NA,NA)
          } else {
            G_tracks_ls[[i]][j,vars_stats] <- as.numeric(
              c(lapply(var_cols, mean, na.rm=TRUE),
                lapply(var_cols, min, na.rm=TRUE),
                lapply(var_cols, max, na.rm=TRUE),
                lapply(var_cols, quantile, probs=percentile/100, na.rm=TRUE)))
          }
        }
        # first NULL row gets second
        G_tracks_ls[[i]][1,vars_stats] <- G_tracks_ls[[i]][2,][vars_stats]
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)
    # drop NAs and merge
    G_tracks_ls <- G_tracks_ls[!is.na(G_tracks_ls)]
    G_tracks <- do.call('rbind', G_tracks_ls)
    G_tracks$dtUTC_lag <- NULL
  }
  # merge all
  if (nrow(L_tracks) > 0 & nrow(G_tracks) > 0){
    tracks <- rbind(G_tracks, L_tracks)
  } else if (nrow(L_tracks) > 0){
    tracks <- L_tracks
  } else if (nrow(G_tracks) > 0){
    tracks <- G_tracks
  } else {
    message('Error: no tracks returned..')
    return(NULL)
  }
  row.names(tracks) <- NULL
  # replace NaN and Inf with NA
  tracks <- do.call(data.frame,lapply(tracks, function(x) replace(x, is.infinite(x),NA)))
  tracks <- do.call(data.frame,lapply(tracks, function(x) replace(x, is.nan(x),NA)))
  # remove tracks that are all NA
  if (filter_nodive_tracks){
    # get first col with a varname to check
    test_col <- colnames(tracks)[grep(vars[1], colnames(tracks))][1]
    tracks <- split(tracks, droplevels(tracks$id))
    for (i in 1:length(tracks)){
      if (sum(is.na(tracks[[i]][,test_col])) == nrow(tracks[[i]])){
        warning(tracks[[i]]$id,' has no dive data. Track dropped.')
        tracks[[i]] <- NA
      }
    }
    tracks <- tracks[!is.na(tracks)]
    tracks <- do.call('rbind', tracks)
  }
  # any depths below zero get zero
  if ('depth' %in% vars & zero_neg_depths){
    if (method == 'segment_stats'){
      # replace less than zero
      tracks$depth_mean[tracks$depth_mean < 0 & !is.na(tracks$depth_mean)] <- 0
      tracks$depth_min[tracks$depth_min < 0 & !is.na(tracks$depth_min)] <- 0
      tracks$depth_max[tracks$depth_max < 0 & !is.na(tracks$depth_max)] <- 0
      tracks[tracks[,paste0('depth_percen',percentile)] < 0 &
             !is.na(tracks[,paste0('depth_percen',percentile)]),
             paste0('depth_percen',percentile)] <- 0
    } else {
      tracks$depth[tracks$depth < 0 & !is.na(tracks$depth)] <- 0
    }
  }
  # reorganise stats when using segment_stats method
  if (method == 'segment_stats'){
    stat_cols <- (1:(length(vars)*4))+(length(colnames(tracks)) - length(vars)*4)
    new_order <- c(1:(length(colnames(tracks)) - length(vars)*4),
                   stat_cols[c(seq(1,length(stat_cols),2),seq(2,length(stat_cols),2))])
    tracks  <- tracks[,new_order]
  }
  row.names(tracks) <- NULL
  tracks$id <- droplevels(tracks$id)
  return(tracks)
}

######################################
# Calculate vertical speed from dive #
######################################
# returns a vector
dive.vertical.speed <- function(tracks, depth_col='depth'){
  # split tracks
  df_ls <- split(tracks, droplevels(tracks$id))
  # get progressbar length 
  pb_len <- sum(sapply(df_ls, nrow))
  pb_i = 0
  message('Caluating vertical speed..')
  pb <- txtProgressBar(min=0, max=pb_len, style = 3)
  # calc speed for each
  for (i in 1:length(df_ls)){
    # make lag columns
    df_ls[[i]]$dtUTC_lag <- data.table::shift(df_ls[[i]]$dtUTC)
    df_ls[[i]]$depth_lag <- data.table::shift(df_ls[[i]][,depth_col])
    df_ls[[i]]$vertcal_speed <- NA
    # calculate the vertical speed
    for (j in 1:nrow(df_ls[[i]])){
      row <- df_ls[[i]][j,]
      df_ls[[i]][j,'vertical_speed'] <- as.numeric(abs(row[depth_col] - row$depth_lag)) / 
        as.numeric(difftime(row$dtUTC, row$dtUTC_lag, unit='secs'))
      pb_i <- pb_i + 1
      setTxtProgressBar(pb, pb_i)
    }
  }
  close(pb)
  # merge
  vertical_speed <- unlist(lapply(df_ls, function(x) x$vertical_speed))
  return(as.numeric(vertical_speed))
}

##############################
# Zero negative depth values #
##############################
# Zero values that indicate a penguin is flying
zero.flying.penguins <- function(accel, depth_col='depth'){
  accel[accel[,depth_col] < 0 & !is.na(accel[,depth_col]),depth_col] <- 0
  return(accel)
}

#########################
# Pressure to depth (m) #
#########################
# Calculate pressure from depth
# input can be a vector or a single value
# Notes:
# mean seawater density from every survey is 1026 (except for 2015_S1 which is 1025)
# Equation: p = p0 + rho*h*g
# Where p is the pressure at a particular depth, p0 is the pressure of the atmosphere,  
# rho is the density of the fluid, g is the acceleration due to gravity, and h is the depth.
# Equation rearranged: h = (p - p0) / g*rho
# input pressure is millibars which needs to be convereted to Newtons per metre squared (*100)
pressure2depth <- function(p, p0, rho=1026, g=9.807){
  return((p*100 - p0*100) / (g*rho))
}

#####################################
# filter gaps in regularised tracks #
#####################################
# filter gaps in regularised gps tracks
# requires a colum with speed
tracks.reg.filter.gaps <- function(regdata, method=c('replace', 'remove'), round_factor=2,
                                   coordNames=c('lon','lat'), calc_speed=TRUE){
  # calc speed
  if (calc_speed){
    regdata <- track.speed.longcalc(regdata)
  }
  # set arguments
  method <- match.arg(method)
  # round speed
  regdata$r_speed <- sapply(regdata$speed, function(x) round(x, round_factor))
  # find the gaps
  len <- rle(regdata$r_speed)$lengths
  v <- NULL
  for(i in 1:length(len)){
    v <- c(v,rep(len[i],len[i]))
  }
  # either replace gaps with NAs or..
  if (method == 'replace')
    regdata[which(v>1), coordNames] <- NA
  # ..remove gaps
  if (method == 'remove')
    regdata[which(v>1), coordNames] <- NULL
  # clean up
  message(sum(len > 1),' gaps identified, ',sum(v > 1),' rows ',method,'d')
  regdata$r_speed <- NULL
  return(regdata)
}

############################
# Attach "DivesPerSegment" #
############################
# DEPRECATED
# calculates the number of dives completed per segement
# Axy dive start and dive end depth (m) is set to macth CEFAS fast log settings
legacy.tracks_divesPerSegement <- function(tracks, dive_start=2, dive_end=1.5,
                                   axy_dive_fn='./data/dive/Axy/penguins_depth_temp_axy.rds',
                                   CEFAS_dive_dir='./data/dive/CEFAS/'){
  # separate Axy and CEFAS tracks
  message('Extracting dive count per segment from Axy dive data..')
  # split tracks into Axy and CEFAS
  L_tracks <- tracks[substr(tracks$id,1,1) == 'L',]
  G_tracks <- tracks[substr(tracks$id,1,1) == 'G',]
  
  if (nrow(L_tracks) > 0){
    # Process L_Tracks
    AXY_dive <- readRDS(axy_dive_fn)
    L_tracks$id <- droplevels(L_tracks$id)
    L_tracks <- split(L_tracks, L_tracks$id)
    pb <- txtProgressBar(min=0, max=length(L_tracks), style = 3)
    setTxtProgressBar(pb, 1)
    for (i in 1:length(L_tracks)){
      track <- L_tracks[[i]]
      # get the dive portion
      dive_data <- AXY_dive[as.character(AXY_dive$id) == as.character(track$id[1]),]
      row.names(dive_data) <- NULL
      dive_data$dive_switch <- NA
      dive_data$dive_state <- NA
      # find dive start and dive end points
      # initial states
      prev_depth <- dive_data$depth[1]
      if (prev_depth > dive_start){
        dive_data$dive_switch[1] <- 'start'
        dive_data$dive_state[1] <- TRUE
      } else {
        dive_data$dive_state[1] <- FALSE
      }
      state_dive <- dive_data$dive_state[1]
      # iterations
      for (j in 2:nrow(dive_data)){
        # if not in dive state check for dive start
        if (!state_dive){
          if (dive_data$depth[j] > dive_start & prev_depth < dive_start){
            state_dive <- TRUE
            dive_data$dive_switch[j] <- 'start'
          }
        # if in dive state check for dive end  
        } else {
          if (dive_data$depth[j] < dive_end & prev_depth > dive_end){
            state_dive <- FALSE
            dive_data$dive_switch[j] <- 'end'
          }
        }
        # set data for next loop
        dive_data$dive_state[j] <- state_dive
        prev_depth <- dive_data$depth[j]
      }
      # calculate dives per segment
      # calculate lag
      track$dtUTC_lag <- data.table::shift(track$dtUTC)
      # calulate mean variables for each segment
      track$dive_count <- NA
      for (j in 2:nrow(track)){
        track$dive_count[j] <- sum(suppressWarnings(dive_data[dive_data$dtUTC > track[j,'dtUTC_lag'] &
                                      dive_data$dtUTC < track[j,'dtUTC'],]$dive_switch == 'start'),na.rm=T)
      }
      L_tracks[[i]]$dive_count <- track$dive_count
      setTxtProgressBar(pb, i)
    }
    close(pb)
    # merge all back together
    L_tracks <- do.call(rbind, L_tracks)
    L_tracks$id <- factor(L_tracks$id)
    row.names(L_tracks) <- NULL
  }
  
  # Now work on G_tracks
  message('Extracting dives count per segment from CEFAS dive data..')
  if (nrow(G_tracks) > 0){
    G_tracks$id <- droplevels(G_tracks$id)
    G_tracks <- split(G_tracks, G_tracks$id)
    pb <- txtProgressBar(min=0, max=length(G_tracks), style = 3)
    setTxtProgressBar(pb, 1)
    # process each track
    for (i in 1:length(G_tracks)){
      track <- G_tracks[[i]]
      # load the dive data
      dive_data <- readRDS(paste0(CEFAS_dive_dir,track$id[1],'.rds'))
      dive_data <- dive_data[dive_data$dive_id != 0,]
      # Because CEFAS splits dives just need to count the number of unique dive ids per segment
      track$dtUTC_lag <- data.table::shift(track$dtUTC)
      # calulate mean variables for each segment
      track$dive_count <- NA
      for (j in 2:nrow(track)){
        track$dive_count[j] <- length(unique(suppressWarnings(dive_data[dive_data$dtUTC > track[j,'dtUTC_lag'] &
                                                                dive_data$dtUTC < track[j,'dtUTC'],]$dive_id)))
      }
      G_tracks[[i]]$dive_count <- track$dive_count
      setTxtProgressBar(pb, i)
    }
    close(pb)
    # merge all back together
    G_tracks <- do.call(rbind, G_tracks)
    G_tracks$id <- factor(G_tracks$id)
    row.names(G_tracks) <- NULL
  }
  
  # merge and return final results
  if (nrow(L_tracks) > 0 & nrow(G_tracks) > 0){
    tracks <- rbind(G_tracks, L_tracks)
  } else if (nrow(L_tracks) > 0){
    tracks <- L_tracks
  } else if (nrow(G_tracks) > 0){
    tracks <- G_tracks
  } else {
    message('Error: no tracks returned..')
    return(NULL)
  }
  row.names(tracks) <- NULL
  tracks$id <- droplevels(tracks$id)
  tracks$dtUTC_lag <- NULL
  return(tracks)
}

#####################
# Attach dive stats #
#####################
# Calculate diving statistics per track segement
# - dive count per unit time
# - dive duration (sum, mean and max)
# - dive depth (mean and max)
#   **To be added**
#     - Temperature (min, mean, max, 85th percentile)
# If return_dive_df set to TRUE also returns a dataframe of dive summaries
# Diving dataframe returns all dives with
# - duration
# - depth (mean and max)
#   **To be added**
#     - Temperature (min, mean, max, variance, temp at max depth, temp at min depth)
# Axy dive start and dive end depth (m) is set to macth CEFAS fast log settings
# (2 m start and 1.5 m end)
tracks_diveStats <- function(tracks, dive_start=2, dive_end=1.5,
                             axy_dive_fn='./data/dive/Axy/penguins_depth_temp_axy.rds',
                             CEFAS_dive_dir='./data/dive/CEFAS/',
                             return_dive_df=FALSE,
                             quality_control=TRUE,
                             quality_ls=list(min_duration=2)){
  # separate Axy and CEFAS tracks
  message('Extracting dive statistics from Axy dive data..')
  # split tracks into Axy and CEFAS
  L_tracks <- tracks[substr(tracks$id,1,1) == 'L',]
  G_tracks <- tracks[substr(tracks$id,1,1) == 'G',]
  
  #------------------#
  # Process L tracks #
  #------------------#
  if (nrow(L_tracks) > 0){
    # Process L_Tracks
    AXY_dive <- readRDS(axy_dive_fn)
    L_tracks$id <- droplevels(L_tracks$id)
    L_tracks <- split(L_tracks, L_tracks$id)
    pb <- txtProgressBar(min=0, max=length(L_tracks), style = 3)
    setTxtProgressBar(pb, 1)
    L_dive_df_ls <- list()
    for (i in 1:length(L_tracks)){
      track <- L_tracks[[i]]
      # get the dive portion
      dive_data <- AXY_dive[as.character(AXY_dive$id) == as.character(track$id[1]),]
      row.names(dive_data) <- NULL
      dive_data$dive_switch <- NA
      dive_data$dive_state <- NA
      # find dive start and dive end points
      # initial states
      prev_depth <- dive_data$depth[1]
      if (prev_depth > dive_start){
        dive_data$dive_switch[1] <- 'start'
        dive_data$dive_state[1] <- TRUE
      } else {
        dive_data$dive_state[1] <- FALSE
      }
      state_dive <- dive_data$dive_state[1]
      # iterations
      for (j in 2:nrow(dive_data)){
        # if not in dive state check for dive start
        if (!state_dive){
          if (dive_data$depth[j] > dive_start & prev_depth < dive_start){
            state_dive <- TRUE
            dive_data$dive_switch[j] <- 'start'
          }
          # if in dive state check for dive end  
        } else {
          if (dive_data$depth[j] < dive_end & prev_depth > dive_end){
            state_dive <- FALSE
            dive_data$dive_switch[j] <- 'end'
          }
        }
        # set data for next loop
        dive_data$dive_state[j] <- state_dive
        prev_depth <- dive_data$depth[j]
      }
      
      # Calculate dive statistics
      # get indexes for dive calculations
      dive_idx <- which(!is.na(dive_data$dive_switch))
      # calculate lag
      track$dtUTC_lag <- data.table::shift(track$dtUTC)
      # Get dataframe of dives
      # make sure start end end count is equal
      switch_summary <- summary(factor(dive_data$dive_switch))
      dive_df <- data.frame(start=dive_data$dtUTC[dive_data$dive_switch == 'start' &
                                                              !is.na(dive_data$dive_switch)])
      if (switch_summary['start'] == switch_summary['end']){
        dive_df$end <- dive_data$dtUTC[dive_data$dive_switch == 'end' &
                                                   !is.na(dive_data$dive_switch)]
      } else {
        dive_df$end <- c(dive_data$dtUTC[dive_data$dive_switch == 'end' &
                                                   !is.na(dive_data$dive_switch)], NA)
      }
      # calculate duration
      dive_df$dive_duration <- as.numeric(difftime(dive_df$end,
                                                   dive_df$start,
                                                    units='secs'))
      # QUALITY CHECK
      # drop tracks less than 2 seconds in length
      if (quality_control)
        dive_df <- dive_df[dive_df$dive_duration > quality_ls$min_duration,]
      
      # calculate dive depth mean and max
      for (j in 1:nrow(dive_df)){
        sub_dive <- dive_data[dive_data$dtUTC > dive_df[j,'start'] &
                              dive_data$dtUTC < dive_df[j,'end'],]
        dive_df$depth_mean[j] <- suppressWarnings(mean(sub_dive$depth, na.rm=T))
        dive_df$depth_max[j] <- suppressWarnings(max(sub_dive$depth, na.rm=T))
      }
      dive_df$depth_mean[is.nan(dive_df$depth_mean)] <- 1.5
      dive_df$depth_max[is.infinite(dive_df$depth_max)] <- 1.5

      # Calculate averages for the track segments
      for (j in 2:nrow(track)){
        dive_sub <- suppressWarnings(dive_df[dive_df$start >= track[j,'dtUTC_lag'] &
                                               dive_df$end <= track[j,'dtUTC'],])
        # if no dives move on
        if (nrow(dive_sub) == 0){
          track[j,'diveDepth_mean'] <- NA
          track[j, 'diveDepth_max'] <- NA
          track[j, 'diveDuration_mean'] <- NA
          track[j, 'diveDuration_max'] <- NA
          track[j, 'diveDuration_sum'] <- NA
          track[j, 'diveCount'] <- 0
        } else {
          track[j,'diveDepth_mean'] <- mean(dive_sub$depth_mean)
          track[j, 'diveDepth_max'] <- max(dive_sub$depth_max)
          track[j, 'diveDuration_mean'] <- mean(dive_sub$dive_duration)
          track[j, 'diveDuration_max'] <- max(dive_sub$dive_duration)
          track[j, 'diveDuration_sum'] <- sum(dive_sub$dive_duration)
          track[j, 'diveCount'] <- nrow(dive_sub)
        }
      }
      L_tracks[[i]] <- track
      setTxtProgressBar(pb, i)
      if (return_dive_df){
        dive_df$id <- track$id[1]
        L_dive_df_ls[[i]] <- dive_df
      }
    }
    close(pb)
    # merge all back together
    L_tracks <- do.call(rbind, L_tracks)
    L_tracks <- rows.and.levels(L_tracks)
    L_dive_df <- do.call(rbind, L_dive_df_ls)
  }
  
  #------------------#
  # Process G_tracks #
  #------------------#
  message('Extracting dives statistics from CEFAS data..')
  G_dive_df_ls <- list()
  if (nrow(G_tracks) > 0){
    G_tracks$id <- droplevels(G_tracks$id)
    G_tracks <- split(G_tracks, G_tracks$id)
    pb <- txtProgressBar(min=0, max=length(G_tracks), style = 3)
    setTxtProgressBar(pb, 1)
    # process each track
    no_Gdata <- list()
    k <- 1
    for (i in 1:length(G_tracks)){
      track <- G_tracks[[i]]
      # load the dive data (if it exists)
      if (file.exists(paste0(CEFAS_dive_dir,track$id[1],'.rds'))){
        dive_data <- readRDS(paste0(CEFAS_dive_dir,track$id[1],'.rds'))
        dive_data <- dive_data[dive_data$dive_id != 0,]
        # Because CEFAS splits dives just need to count the number of unique dive ids per segment
        track$dtUTC_lag <- data.table::shift(track$dtUTC)
  
        # build dive data frame
        dive_df <- data.frame(start=aggregate(dive_data$dtUTC, 
                                              by=list(dive_id=dive_data$dive_id), 
                                              FUN=min)$x)
        dive_df$end  <- aggregate(dive_data$dtUTC, 
                                  by=list(dive_id=dive_data$dive_id), 
                                  FUN=max)$x
        # convert back to UTC
        attr(dive_df$start, 'tzone') <- 'UTC'
        attr(dive_df$end, 'tzone') <- 'UTC'
        # calculate durations
        dive_df$dive_duration <- as.numeric(difftime(dive_df$end,
                                                     dive_df$start,
                                                     units='secs'))
        
        # QUALITY CHECK
        # dives shorter than 2 seconds can fuck off. 
        if (quality_control)
          dive_df <- dive_df[dive_df$dive_duration > quality_ls$min_duration,]
        
        # calculate dive depth mean and max
        for (j in 1:nrow(dive_df)){
          sub_dive <- dive_data[dive_data$dtUTC > dive_df[j,'start'] &
                                  dive_data$dtUTC < dive_df[j,'end'],]
          dive_df$depth_mean[j] <- suppressWarnings(mean(sub_dive$depth, na.rm=T))
          dive_df$depth_max[j] <- suppressWarnings(max(sub_dive$depth, na.rm=T))
        }
        dive_df$depth_mean[is.nan(dive_df$depth_mean)] <- 1.5
        dive_df$depth_max[is.infinite(dive_df$depth_max)] <- 1.5
        
        # Calculate averages for the track segments
        for (j in 2:nrow(track)){
          dive_sub <- suppressWarnings(dive_df[dive_df$start >= track[j,'dtUTC_lag'] &
                                                 dive_df$end <= track[j,'dtUTC'],])
          # if no dives move on
          if (nrow(dive_sub) == 0){
            track[j,'diveDepth_mean'] <- NA
            track[j, 'diveDepth_max'] <- NA
            track[j, 'diveDuration_mean'] <- NA
            track[j, 'diveDuration_max'] <- NA
            track[j, 'diveDuration_sum'] <- NA
            track[j, 'diveCount'] <- 0
          } else {
            track[j,'diveDepth_mean'] <- mean(dive_sub$depth_mean)
            track[j, 'diveDepth_max'] <- max(dive_sub$depth_max)
            track[j, 'diveDuration_mean'] <- mean(dive_sub$dive_duration)
            track[j, 'diveDuration_max'] <- max(dive_sub$dive_duration)
            track[j, 'diveDuration_sum'] <- sum(dive_sub$dive_duration)
            track[j, 'diveCount'] <- nrow(dive_sub)
          }
        }
        G_tracks[[i]] <- track
        setTxtProgressBar(pb, i)
        if (return_dive_df){
          dive_df$id <- track$id[1]
          G_dive_df_ls[[i]] <- dive_df
        }
      } else {
        G_tracks[[i]] <- NA
        no_Gdata[k] <- as.character(track$id[1])
        k <- k + 1
        if (return_dive_df){
          G_dive_df_ls[[i]] <- NA
        }
      }
    }
    close(pb)
    # Clean failed 
    G_tracks <- G_tracks[!is.na(G_tracks)]
    if (return_dive_df){
      G_dive_df_ls <- G_dive_df_ls[!is.na(G_dive_df_ls)]
    }
    if (length(no_Gdata > 0)){
      warning('No dive data found for ',no_Gdata,'. Tracks dropped')
    }
    # merge all back together
    G_tracks <- do.call(rbind, G_tracks)
    G_tracks <- rows.and.levels(G_tracks)
    G_dive_df <- do.call(rbind, G_dive_df_ls)  
  }
  # merge and return final results
  if (nrow(L_tracks) > 0 & nrow(G_tracks) > 0){
    tracks <- rbind(G_tracks, L_tracks)
    dive_df <- rbind(G_dive_df, L_dive_df)
  } else if (nrow(L_tracks) > 0){
    tracks <- L_tracks
    dive_df <- L_dive_df
  } else if (nrow(G_tracks) > 0){
    tracks <- G_tracks
    dive_df <- G_dive_df
  } else {
    message('Error: no tracks returned..')
    return(NULL)
  }
  tracks <- rows.and.levels(tracks)
  tracks$dtUTC_lag <- NULL
  if (return_dive_df){
    # reorder dive_df
    dive_df <- dive_df[,c(6,1:5)]
    # keep only complete dive info
    dive_df <- dive_df[complete.cases(dive_df),]
    return(list(tracks, dive_df))
  } else {
    return(tracks)
  }
}

#######################################
# Filter tracks based on survey times #
#######################################
# window size is in days
# Update: window can be a integer or a named list (for each year)
# Update: balanceTime overwrites window list and will make the window size equal
# based on survey lengths
tracks.survey.time.filter <- function(tracks, window=5, report_only=FALSE, balanceTime=F,
                                      survey_time_fn='./data/transects/survey_times.csv'){
  # Get the mean datetime
  tracks_mean <- tracks[,c('id','dtUTC')]
  tracks_mean <- aggregate(tracks_mean, by=list(tracks_mean$id), FUN=mean)
  tracks_mean$id <- tracks_mean$Group.1
  tracks_mean$Group.1 <- NULL
  
  # Attach closest survey data
  tracks_mean <- tracks.nearest.survey(tracks_mean, survey_time_fn=survey_time_fn)
  tracks_mean <- rows.and.levels(tracks_mean, idcol='survey_id')
  # load survey data 
  survey_times <- load.survey.times()
  # Only keep relevant times
  survey_times <- survey_times[survey_times$id %in% 
                                 unique(as.character(tracks_mean$survey_id)),]
  survey_times <- rows.and.levels(survey_times)
  
  # if balance time true chnage the window
  if (balanceTime){
    survey_times_2 <- survey_times
    survey_times_2$type  <- c(rep('start', nrow(survey_times)/2), rep('end', nrow(survey_times)/2))
    survey_times_2 <- split(survey_times_2, survey_times_2$type)
    survey_times_join <- survey_times_2$start[,1:3]
    names(survey_times_join) <- c('id','year','start')
    survey_times_join$end <- survey_times_2$end$dtUTC
    survey_times_join$diff <- as.integer(ceiling(difftime(survey_times_join$end, 
                                                          survey_times_join$start,
                                                          units = 'day')))
    # Make new window list
    window <- max(survey_times_join$diff)  - survey_times_join$diff
    message('\nWindow vector for balanced time intervals')
    print(window)
  }
  
  # Integer version
  if (is.numeric(window)){
    # check to see if each track is within N days of survey
    tracks_mean$in_window <- NA
    for (i in 1:nrow(tracks_mean)){
      # get survey to check
      st <- survey_times[survey_times$id == tracks_mean$survey_id[i],]
      # check time diffs (neg is before, pos is after)
      st$diff <- as.numeric(difftime(st$dtUTC, tracks_mean$dtUTC[i], units='days'))
      # TRUE if between dates, or N days before (1) or N days after (2)
      tracks_mean$in_window[i] <- (st$diff[1] <= 0 & st$diff[2] >= 0) |
                                  (abs(st$diff[1]) <= window) | (abs(st$diff[2]) <= window)
    }
  # Named list version
  } else if (is.list(window)){
    # split and see if each track is within N days of survey
    tracks_mean$in_window <- NA
    tracks_mean.ls <- split(tracks_mean, factor(tracks_mean$survey_id))
    for (j in 1:length(tracks_mean.ls)){
      for (i in 1:nrow(tracks_mean.ls[[j]])){
        # get survey to check
        st <- survey_times[survey_times$id == tracks_mean.ls[[j]]$survey_id[i],]
        # check time diffs (neg is before, pos is after)
        st$diff <- as.numeric(difftime(st$dtUTC, tracks_mean.ls[[j]]$dtUTC[i], units='days'))
        # TRUE if between dates, or N days before (1) or N days after (2)
        tracks_mean.ls[[j]]$in_window[i] <- (st$diff[1] <= 0 & st$diff[2] >= 0) |
          (abs(st$diff[1]) <= window[substr(tracks_mean.ls[[j]]$survey_id[i],1,4)]) | 
          (abs(st$diff[2]) <= window[substr(tracks_mean.ls[[j]]$survey_id[i],1,4)])
      }
    }
    tracks_mean <- do.call(rbind, tracks_mean.ls)
    window='variable'
  }
  
  
  # report 
  pre_df <- aggregate(tracks_mean[,'dtUTC'], 
                      by=list(year(tracks_mean$dtUTC)), FUN=length)
  colnames(pre_df) <- c('Year', 'Tracks')
  report_df <- aggregate(tracks_mean[,'in_window'], 
                         by=list(year(tracks_mean$dtUTC)), FUN=sum)
  colnames(report_df) <- c('Year', 'Tracks')
  message('Window of ',window,' days will filter ',sum(pre_df$Tracks) - sum(report_df$Tracks), ' tracks ',
          '- Before: ',sum(pre_df$Tracks),', After: ',sum(report_df$Tracks))
  message('Tracks beofre:')
  print(pre_df)
  message('Tracks after:')
  print(report_df)
  
  # return results
  if (!report_only){
    message('Filter complete.')
    tracks_result <- tracks[tracks$id %in% tracks_mean$id[tracks_mean$in_window],]
    tracks_result$id <- droplevels(tracks_result$id)
    return(tracks_result)
  }
}

#########################################
# Link penguin tracks to nearest survey #
#########################################
tracks.nearest.survey <- function(tracks, survey_time_fn='./data/transects/survey_times.csv',
                                  dtUTCcol='dtUTC'){
  # load survey times
  survey_times <- read.csv(survey_time_fn)
  survey_times$start_UTC <- as.POSIXct(survey_times$start_UTC, tz='UTC')
  survey_times$end_UTC <- as.POSIXct(survey_times$end_UTC, tz='UTC')
  # stack columns
  survey_times_2 <- survey_times
  survey_times_2$start_UTC <- survey_times_2$end_UTC
  survey_times <- rbind(survey_times, survey_times_2)
  survey_times$end_UTC <- NULL
  colnames(survey_times)[3] <- 'dtUTC'
  # find which survey each track should pull data from
  track_survey_links <- data.frame(track_id=unique(tracks$id))
  track_survey_links$mean_dtUTC <- as.POSIXct(sapply(track_survey_links$track_id, 
                                                     function(id) mean(tracks[,dtUTCcol][tracks$id == id])),
                                              origin='1970-01-01', tz='UTC')
  # find survey from time
  track_survey_links$survey_id <- sapply(track_survey_links$mean_dtUTC,
                                         function(mean_dt) survey_times$id[which.min(abs(
                                           difftime(survey_times$dtUTC, mean_dt, units='hour')))])
  # attach to tracks
  tracks$survey_id <- sapply(tracks$id, function(id) 
    track_survey_links$survey_id[track_survey_links$track_id == id])
  return(tracks)
}

#####################
# Load survey times #
#####################
load.survey.times <- function(survey_time_fn='./data/transects/survey_times.csv', process=T){
  # load survey times
  survey_times <- read.csv(survey_time_fn)
  survey_times$start_UTC <- as.POSIXct(survey_times$start_UTC, tz='UTC')
  survey_times$end_UTC <- as.POSIXct(survey_times$end_UTC, tz='UTC')
  if (process){
    # stack columns
    survey_times_2 <- survey_times
    survey_times_2$start_UTC <- survey_times_2$end_UTC
    survey_times <- rbind(survey_times, survey_times_2)
    survey_times$end_UTC <- NULL
    colnames(survey_times)[3] <- 'dtUTC'
  }
  return(survey_times)
}

####################
# reset_id factors #
####################
# I do this alot so just a shorcut function
# Resets the index and drops extra id factor levels
rows.and.levels <- function(df, idcol='id'){
  df[, idcol] <- droplevels(df[, idcol])
  row.names(df) <- NULL
  return(df)
}

##############################################
# Remove tracks that overlap certain surveys #
##############################################
tracks.filter.surveys <- function(tracks, survey_rm, keep=FALSE){
  # find closest surveys (if haven't already)
  if (!'survey_id' %in% colnames(tracks))
    tracks <- tracks.nearest.survey(tracks)
  # keep instead of filter?
  if (keep){
    tracks <- tracks[tracks$survey_id %in% survey_rm,]
  } else {
    tracks <- tracks[!tracks$survey_id %in% survey_rm,]
  }
  return(tracks)
}

#################################
# Make relative datetime stamps #
#################################
# Note: AEST is calculated as GMT+10 and ignores dayloght savings
relatvie.timestamps <- function(dt, makeAEST=TRUE){
  # make local
  tzone <- "UTC"
  if (makeAEST){
    warning('AEST is calculated for Australia/Brisbane to ignore daylight savings.')
    tzone <- "Australia/Brisbane"
    attr(dt, "tzone") <- tzone
  }
  # make relative datetimes
  hours <- hour(dt)
  mins <- minute(dt)
  secs <- second(dt)
  dt <- ISOdatetime(2000, 1, 1, hours, mins, secs, tz=tzone)
  return(dt)
}

####################################
# Cumulative sum for track columns #
####################################
tracks.cumsum <- function(tracks, col){
  track_ls <- split(tracks, tracks$id)
  for (i in 1:length(track_ls)){
    track_ls[[i]][,paste0(col,'_cumsum')] <- cumsum(track_ls[[i]][,col])
  }
  tracks <- do.call(rbind, track_ls)
  tracks <- rows.and.levels(tracks)
  return(tracks)
}

################################
# Transform Sv mean to linear  #
################################
Sv_mean.linear <- function(Sv){
  return(10**(Sv/10))
}

################################
# Transform Sv mean to linear  #
################################
Sv_mean.log <- function(Sv){
  return(10*log10(Sv))
}

##############################
# Cumulative sum for Sv mean #
##############################
# needs a special function because more positive is greater density therefore 
# it doesn't sum right. Also needs to be transformed to linear from the 
# log domain
tracks.cumsum.Sv <- function(tracks, Sv_col='Sv_linear'){
  track_ls <- split(tracks, tracks$id)
  for (i in 1:length(track_ls)){
    Sv <- track_ls[[i]]$Sv_linear
    track_ls[[i]][,paste0(Sv_col,'_cumsum')] <- cumsum(Sv)
  }
  tracks <- do.call(rbind, track_ls)
  tracks <- rows.and.levels(tracks)
  return(tracks)
}

######################################
# Add track index to a set of tracks #
######################################
# Gives an index for each track
tracks.add.indexes <- function(tracks, label='rank'){
  tracks <- rows.and.levels(tracks)
  track_ls <- split(tracks, tracks$id)
  for (i in 1:length(track_ls)){
    track_ls[[i]][,label] <- 1:nrow(track_ls[[i]])
  }
  tracks <- do.call(rbind, track_ls)
  tracks <- rows.and.levels(tracks)
  return(tracks)
}

######################################
# Compute penguin summary statistics #
######################################
peng.summary.stats <- function(tracks, add.last.point=T, 
                               last.point=c(150.2269, -36.25201)){
  # dive tracks by their id
  track.ls <- split(tracks, tracks$id)
  # Create a dataframe to capture results
  sumStats <- data.frame(id=unique(tracks$id))
  sumStats$survey_id <- unlist(lapply(track.ls, function(df) df$survey_id[1]))
  sumStats$year <- as.factor(unlist(lapply(sumStats$survey_id, function(s) substr(s,1,4))))
  
  # __Trip Duration__ #
  # unit is hours
  message('Calculating trip duration..')
  sumStats$duration <- unlist(pblapply(track.ls, function(df) 
    round(sum(track.time(df)$time, na.rm=T), 2)/60/60))
  
  # __Maximum Displacement__ #
  # unit is km
  message('Calculating maximum displacement..')
  sumStats$displacement_max <- unlist(pblapply(track.ls, function(df) 
    round(max(track.displacement(df)$disp, na.rm=T), 0)/1000))
  
  # add return journey so distance not skewed
  message('Adding return journey point (',paste(last.point),') for incomplete trips..')
  if (add.last.point){
    for (i in 1:length(track.ls)){
      row <- track.ls[[i]][nrow(track.ls[[i]]),]
      row$lon <- last.point[1]
      row$lat <- last.point[2]
      track.ls[[i]] <- rbind(track.ls[[i]], row)
    }
  }
  
  # __Distance Travelled__ #
  # unit is km
  message('Calculating total distance travelled..')
  sumStats$distance_total <- unlist(pblapply(track.ls, function(df) 
    round(sum(track.distance(df)$dist, na.rm=T), 0)/1000))
  
  # Return the results
  return(sumStats)
}

##########################
# Load kriging databases #
##########################
load.krig.dbs <- function(var=c('temperature','salinity','Sv_mean'), 
                          surveys='all', dir='./kriging/output/', 
                          return_items=F, survey_col=F){
  # set arguments
  var <- match.arg(var)
  # set surveys
  if ('all' %in% surveys){
    surveys <- c('2015_S1','2016_S1','2016_S2','2017_S1','2017_S2',
                 '2018_S1','2018_S2','2019_S1','2019_S2')
  }
  # load databases
  message('Loading ',var,' databases')
  if (var == 'temperature'){
    dbs <- pblapply(surveys, function(s) readRDS(paste0(dir,'CTD/models/3Dkrig_CTDdb_temperature_',
                                                           s,'.rds')))
  } else if (var == 'salinity'){
    dbs <- pblapply(surveys, function(s) readRDS(paste0(dir,'CTD/models/3Dkrig_CTDdb_salinity_',
                                                           s,'.rds')))
  } else if (var == 'Sv_mean'){
    dbs <- pblapply(surveys, function(s) readRDS(paste0(dir,'agg/models/3Dkrig_aggdb_Sv_mean_',
                                                           s,'.rds')))
  }
  # add survey ids
  if (survey_col){
    for (i in 1:length(dbs)){
      dbs[[i]]@items$survey_id <- surveys[i]
      row.names(dbs[[i]]@items) <- NULL
    }
  }
  # convert to dataframe
  if (return_items){
    dbs <- lapply(dbs, function(db) db@items)
  }
  # name by survey
  names(dbs) <- surveys
  return(dbs)
}

##################################################################################
# Load data from kriging models for penguin tracks with the mean taken for depth #
##################################################################################
# if survey not provided survey is looked for in survey_id column
# If passing weights pass a dataframe with depth, freq and survey
db.extract.covar.depth <- function(tracks, depth, convert2UTM=TRUE,
                                   var=c('Sv_mean','temperature','salinity'), 
                                   method=c('mean','slice','weight'), 
                                   krig.dir='./kriging/output/',
                                   weight.df=NULL,
                                   survey=NULL){
  # Set variables
  var <- match.arg(var)
  method <- match.arg(method)
  
  # Get a list of surveys to read in
  if (!is.null(survey)){
    tracks$survey_id <- factor(survey)
    warning('Survey id provided. Survey_id columns will be overwritten if present.')
    
  }
  survey.ls <- as.character(unique(tracks$survey_id))
  
  # Load the databases
  dbs <- load.krig.dbs(var, surveys=survey.ls)
  
  # Create UTM coordinates (filter for 2015 before function)
  if (convert2UTM)
    tracks <- na.sub.FUN.df(tracks, c('lon','lat'), lonlat2UTM.df)
  
  # Prepare databases based on method
  dbs <- lapply(dbs, function(db) db@items[db@items$Polygon,])
  if (method == 'slice'){
    dbs <- lapply(dbs, function(db) db[db$x3 == depth,])
    # Take the mean over x3
    dbs <- lapply(dbs, function(db) aggregate(db, by=list(db$x1, db$x2), FUN=mean))
  } else if (method == 'mean') {
    dbs <- lapply(dbs, function(db) db[db$x3 <= depth,])
    # Take the mean over x3
    dbs <- lapply(dbs, function(db) aggregate(db, by=list(db$x1, db$x2), FUN=mean))
  } else if (method == 'weight'){
    dbs <- lapply(dbs, function(db) db[db$x3 <= max(weight.df$depth),])
    # Weight by dive count
    for (i in 1:length(dbs)){
      dbs[[i]] <- aggregate(dbs[[i]], by=list(dbs[[i]]$x1, dbs[[i]]$x2),
              FUN = function(x) weighted.mean(x, weight.df$freq[weight.df$survey_id == names(dbs)[i]]))
    }
  }
  
  # Remove redundant variables
  for (i in 1:length(dbs)){
    dbs[[i]] <- subset(dbs[[i]], select=-c(rank,x1,x2,x3,Polygon))
    names(dbs[[i]])[1:2] <- c('x1','x2')
    dbs[[i]]$survey_id <- names(dbs)[i]
  }
  # Merge into single dataframe
  df <- do.call(rbind, dbs)
  row.names(df) <- NULL
  df$survey_id <- as.factor(df$survey_id)
  
  # Remove track NAs for kd-tree
  # get na.idx
  na.idx <- which(is.na(tracks$lon))
  dat.idx <- which(!is.na(tracks$lon))
  # filter NAs
  track.filt <- tracks[dat.idx,]
  # Match index using kd-trees
  df.idx <- nn2(df[,c('x1','x2','survey_id')], 
                track.filt[,c('lon', 'lat','survey_id')], k=1)$nn.idx
  # sub back in NAs and get covariate data
  tracks[dat.idx, paste0('krig_',var)] <- df[df.idx, paste0('Kriging.',var,'.estim')]

  return(tracks)
}

##########################################
# Report number of tracks by some factor #
##########################################
tracksCount.report <- function(tracks, split.var){
  split.var <- as.factor(split.var)
  # split by var
  track.ls <- split(tracks, split.var)
  df.count <- data.frame(var=names(track.ls),
                         count=unlist(lapply(track.ls,
                                      function(df) length(unique(df$id)))))
  row.names(df.count) <- NULL
  message('Track count:')
  print(df.count)
}

##############################################
# Make tracks numebr equal along some factor #
##############################################
equalise.tracks <- function(tracks, split.var, seed=420){
  split.var <- as.factor(split.var)
  # split tracks by the variable (random string so no risk of overwrite)
  tracks$split.var.hdfkl <- split.var
  track.ls <- split(tracks, split.var)
  # count vars
  df.count.before <- data.frame(var=names(track.ls),
                                count=unlist(lapply(track.ls,
                                             function(df) length(unique(df$id)))))
  row.names(df.count.before) <- NULL
  # find the minimum value
  min.count <- min(df.count.before$count, na.rm=T)
  
  # randomly sample from each so we have the minimum
  id.ls <- lapply(track.ls, function(df) unique(df$id))
  set.seed(seed)
  id.ls <- lapply(id.ls, function(ids) ids[sample(1:length(ids), min.count)])
  # filter tracks
  tracks <- tracks[tracks$id %in% unlist(id.ls),]
  
  # split tracks by the variable again
  track.ls <- split(tracks, tracks$split.var.hdfkl)
  # count vars
  df.count.after <- data.frame(var=names(track.ls),
                               count=unlist(lapply(track.ls,
                                            function(df) length(unique(df$id)))))
  row.names(df.count.after) <- NULL
  tracks$split.var.hdfkl <- NULL
  # report
  message('Track count before:')
  print(df.count.before)
  message('Tracks count after:')
  print(df.count.after)
  tracks <- rows.and.levels(tracks)
  return(tracks)
}

#######################################################
# Function to calculate centre of gravity and inertia #
#######################################################
# Author: Gemma Carroll, 2020
# Modified by Lachlan Phillips, 2020
# Notes: 
# Here z would be your sv values
# Result should be a data frame with the center of mass, and the (4) endpoints 
# of the 2 inertia axes
# x and y are lon/lat
calc.CGI <- function(x, y, z=NA, w=1){
  # Center of gravity coordinates
  xg <- sum(x * z * w)/sum(z * w)
  yg <- sum(y * z * w)/sum(z * w)
  
  # Inertia
  dx <- x - xg
  dy <- y - yg
  d <- sqrt(dx^2 + dy^2)
  inert <- sum(z * w * (d^2))/sum(z * w)
  I <- inert
  
  # Weighted PCA
  M11 <- sum(dx^2 * z * w)
  M22 <- sum(dy^2 * z * w)
  M21 <- sum(dx * dy * z * w)
  M12 <- M21
  M <- matrix(c(M11, M12, M21, M22), byrow = T, ncol = 2)
  x1 <- eigen(M)$vectors[1, 1]
  y1 <- eigen(M)$vectors[2, 1]
  x2 <- eigen(M)$vectors[1, 2]
  y2 <- eigen(M)$vectors[2, 2]
  r1 <- eigen(M)$values[1]/(eigen(M)$values[1] + eigen(M)$values[2])
  
  # Principal axis coordinates
  e1 <- (y1/x1)^2
  sx1 <- x1/abs(x1)
  sy1 <- y1/abs(y1)
  sx2 <- x2/abs(x2)
  sy2 <- y2/abs(y2)
  xa <- xg + sx1 * sqrt((r1 * inert)/(1 + e1))
  ya <- yg + sy1 * sqrt((r1 * inert)/(1 + (1/e1)))
  xb <- 2 * xg - xa
  yb <- 2 * yg - ya
  xc <- xg + sx2 * sqrt(((1 - r1) * inert)/(1 + (1/e1)))
  yc <- yg + sy2 * sqrt(((1 - r1) * inert)/(1 + e1))
  xd <- 2 * xg - xc
  yd <- 2 * yg - yc
  
  resdf <- data.frame(name = c("centre_gravity","axis1","axis1","axis2","axis2"),
                    x = c(xg, xa, xb, xc, xd), 
                    y = c(yg, ya, yb, yc, yd))
  return(resdf)
}

########################
# Convert Depth Raster #
########################
# Load in depth raster and convert to CV
depthDB.CV <- function(filename='./data/bathymetry/raster/depthKrigAcoustics.tif',
                       method=c('relativeMean','totalMean','totalMedian'),
                       plotRaster=T,
                       depthFilter=NULL){
  # match arguments
  method <- match.arg(method)
  # load source raster
  bathy <- rast(filename)
  
  #---------------------#
  #   Vertcial Relief   #
  #---------------------#
  # Remove land
  message('Filtering land values..')
  bathy[bathy <= 0] <- NA
  # Filter depth
  if (!is.null(depthFilter)){
    message('Filtering depth > ',depthFilter,'m')
    bathy[bathy < -depthFilter] <- NA
  }
  # Caculate vertical relief and resize by factor of 10
  message('Calculating vertical relief and downsampling raster resolution..')
  if (method == 'relativeMean'){
    CV.raster <- terra::aggregate(bathy, fact=2, fun=function(x) (sd(x)/mean(x)))
  } else if (method == 'totalMean'){
    meanRast <- mean(raster[], na.rm=T)
    CV.raster <- terra::aggregate(bathy, fact=2, fun=function(x) (sd(x)/meanRast))
  } else if (method == 'totalMedian'){
    medianRast <- median(raster[], na.rm=T)
    CV.raster <- terra::aggregate(bathy, fact=2, fun=function(x) (sd(x)/medianRast))
  }
  # Plot if requetsed
  if (plotRaster){
    message('\nShowing plot.')
    plot(CV.raster, col=cmocean('haline')(100))
  }
  return(CV.raster)
}

####################
# NetCDF to raster #
####################
netCDF.2.raster <- function(ncfname='/Volumes/LP_MstrData/master-data/ocean/ROMS/highres/avg_2016.nc',
                            tif.template='./data/bathymetry/raster/depthKrigAcoustics.tif',
                            poly.fn='./kriging/assets/kriging_field/survey-coast-straight_poly.rds',
                            output.fn=NULL){
  # Load the NetCDF file
  message('Opening NetCDF file..')
  ncin <- nc_open(ncfname)
  nc_bathy <- as.vector(ncvar_get(ncin, 'h'))
  nc_lons <- as.vector(ncvar_get(ncin, 'lon_rho'))
  nc_lats <- as.vector(ncvar_get(ncin, 'lat_rho'))
  # make data frame
  nc_df <- data.frame(lon=nc_lons, lat=nc_lats, depth=nc_bathy)
  # convert to UTM
  nc_df <- lonlat2UTM.df(nc_df)
  # load shapefile 
  poly <- readRDS(poly.fn)
  poly <- lonlat2UTM.poly(poly)
  # Capture extent
  poly.extent <- c(min(poly@sets[[1]]@x), min(poly@sets[[1]]@y),
                   max(poly@sets[[1]]@x), max(poly@sets[[1]]@y))
  # Filter to range of template
  nc_df <- nc_df[nc_df$lon >= poly.extent[1] & nc_df$lon <= poly.extent[3] &
                 nc_df$lat >= poly.extent[2] & nc_df$lat <= poly.extent[4],]
  # load template 
  template <- terra::rast(tif.template)
  # set up an 'empty' raster, here via an extent object derived the data
  e <- extent(as.matrix(nc_df[,1:2]))
  r <- raster(e, ncol=ncol(template), nrow=nrow(template))
  # you need to provide a function 'fun' for when there are multiple points per cell
  x <- rasterize(as.matrix(nc_df[,1:2]), r, as.matrix(nc_df[,3]), fun=mean)
  # Load into regular grid
  x <- terra::rast(x)
  crs(x) <- "+proj=utm +zone=56H +ellps=WGS84"
  bathy <- terra::resample(x, template)

}

#####################################
# Load in acoustic aggregation data #
#####################################
process.aggregations <- function(agg.dir="/Volumes/LP_MstrData/master-data/survey/acoustics/acoustics_processed",
                                 out.fn="./data/acoustics/agg.rds"){
  # Becuase data is a bit of a mess right now this is a manual processing job
  # Do my best to retrieve all aggregations and associdated GPS
  # 2015_S1
  files <- data.frame(survey = '2015_S1',
                      files=c("2015_S1_schools.csv", "2015_S1.gps.csv"),
                      type=c('agg', 'gps'))
  # 2016_S1
  files <- rbind(files, 
                 data.frame(survey = '2016_S1',
                            files=c("2016_S1_schools.csv", "2016_S1_2_schools.csv",
                                    "2016_S1.gps.csv", "2016_S1_2.gps.gps.csv"),
                            type=c('agg', 'agg', 'gps', 'gps')))
  # 2016_S2
  files <- rbind(files, 
                 data.frame(survey = '2016_S2',
                            files=c("2016_S2_schools.csv", "2016_S2.gps.csv"),
                            type=c('agg', 'gps')))
  # 2017_S1
  files <- rbind(files, 
                 data.frame(survey = '2017_S1',
                            files=c("2017_S1_schools.csv", "2017_S1_2_schools.csv",
                                    "2017_S1.gps.csv", "2017_S1_2_.gps.gps.csv"),
                            type=c('agg', 'agg', 'gps', 'gps')))
  # 2017_S2
  files <- rbind(files, 
                 data.frame(survey = '2017_S2',
                            files=c("2017_S2_schools.csv", "2017_S2.gps.csv"),
                            type=c('agg', 'gps')))
  # 2018_S1
  files <- rbind(files, 
                 data.frame(survey = '2018_S1',
                            files=c("2018_S1_schools.csv"),
                            type=c('agg')))
  # 2018_S2
  files <- rbind(files, 
                 data.frame(survey = '2018_S2',
                            files=c("2018_S2_schools.csv", "2018_S2_.gps.csv"),
                            type=c('agg', 'gps')))
  # 2019_S1
  files <- rbind(files,
                 data.frame(survey = '2019_S1',
                            files=c("2019_S1_schools.csv", "2019_S1.gps.csv"),
                            type=c('agg', 'gps')))
  # 2019_S2
  list.files(paste0(agg.dir,"/2019_S2"), recursive = T)
  files <- rbind(files,
                 data.frame(survey = '2019_S2',
                            files=c("2019_S2_schools.csv", "2019_S2.gps.csv"),
                            type=c('agg', 'gps')))
  
  # Read in aggregation data
  agg <- split(files[files$type == 'agg',], files[files$type == 'agg','survey'])
  for (i in 1:length(agg)){
    schools <- list()
    for (j in 1:nrow(agg[[i]])){
      schools[[j]] <- read.csv(paste0(agg.dir,'/',names(agg)[i],'/',agg[[i]]$files[j]))
    }
    agg[[i]] <- do.call(rbind, schools)
  }
  
  # Match up data best we can
  # First remove Nb column for 2016_S2 (only survey that has it)
  agg[['2016_S2']]$Nb <- NULL
  
  # For other suverys that have less data add columns that 
  # they are missing from other years
  missing.cols <- names(agg[['2019_S2']])[(ncol(agg[['2016_S2']])+1):ncol(agg[['2019_S2']])]
  surveys.2.process <- names(agg)[lapply(agg, ncol) < 104]
  for (i in 1:length(surveys.2.process)){
    for (col.name in missing.cols){
      agg[[surveys.2.process[i]]][col.name] <- NA
    }
  }
  
  # Add survey names
  for (i in 1:length(agg)){
    agg[[i]]$survey <- names(agg)[i]
  }
  
  # Join all data
  agg <- do.call(rbind, agg)
  row.names(agg) <- NULL
  
  # Save
  saveRDS(agg, out.fn)
}


