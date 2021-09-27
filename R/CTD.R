# CTD processing
# functions for processing CTD data
library(dplyr)
library(measurements)
library(geosphere)
library(ggplot2)
library(lubridate)
library(gatepoints)
library(cmocean)
library(oce)

###################
# CTD cast reader #
###################
# read in CTD cast csv files
cast.reader <- function(input_dir='data/surveys/CTD/casts/',
                        meta_dir='data/surveys/CTD/meta/cast_meta.csv',
                        save_meta=FALSE, drop_no_location=TRUE,
                        add_manual_locations=TRUE, location_from_time=TRUE,
                        check_save=TRUE, drop_unloaded_meta=TRUE,
                        drop_non_transect=TRUE, drop_duplicates=TRUE){
  # read in the file names
  file_ls <- list.files(path=input_dir, pattern='*.csv')
  print(paste(length(file_ls),'cast files found..'))
  # if producing new meta load in meta section first
  if (save_meta){
    # check if file should be overwritten
    if (file.exists(meta_dir) & check_save){
      proceed <- askYesNo(paste0('"',meta_dir,'" already exists. Do you wish to ',
                                 'overwrite?'))
      if (!proceed){
        break
      }
    }
    # firstly go through and read in cast meta data
    meta <- lapply(file_ls, function(fn) readLines(paste0(input_dir, fn), n=27))
    # clean the strings
    meta <- lapply(meta, function(x) substr(x, 3, nchar(x)))
    meta_values <- lapply(meta, function(ls) 
      unname(sapply(ls, function(vect) strsplit(vect, ',')[[1]][2])))
    # make empty vector
    empty <- rep(NA, length(meta))
    # make data frame to hold data
    df_meta <- data.frame(device=empty,
                          file_name=empty,
                          cast_time_UTC=empty,
                          cast_time_local=empty,
                          sample_type=empty,
                          cast_data=empty,
                          location_source=empty,
                          default_latitude=empty,
                          default_altitude=empty,
                          start_latitude=empty,
                          start_longitude=empty,
                          start_altitude=empty,
                          start_gps_horizontal_error=empty,
                          start_gps_vertical_error=empty,
                          start_gps_satellite_no=empty,
                          end_latitude=empty,
                          end_longitude=empty,
                          end_altitude=empty,
                          end_gps_horizontal_error=empty,
                          end_gps_vertical_error=empty,
                          end_gps_satellite_no=empty,
                          cast_duration=empty,
                          samples_per_second=empty,
                          electronics_calibration_date=empty,
                          conductivity_calibration_date=empty,
                          temperature_calibration_date=empty,
                          pressure_calibration_date=empty)
    # populate data frame with meta data
    df_meta[,] <- data.frame(do.call(rbind, meta_values))
    # format data
    df_meta[,8:23] <- apply(df_meta[,8:23], 2, function(x) as.numeric(x))
    df_meta$cast_time_UTC <- as.POSIXct(df_meta$cast_time_UTC, tz='UTC')
    df_meta$cast_time_local <- as.POSIXct(df_meta$cast_time_local, tz='Australia/Sydney')
    # datetimes should be ordered due to filename structure but order 
    # anyway just to be sure
    df_meta <- df_meta[order(df_meta$cast_time_UTC),]
    # add manual locations
    if (add_manual_locations){
      df_meta <- add.manual.locations(df_meta)
    }
    # add transect and station data
    df_meta <- add.transect.station(df_meta)
    # work out remaining cast positons if there are any unknown
    if (location_from_time & 
        nrow(df_meta[df_meta$location_source == 'Manual' & 
             df_meta$sample_type != 'Invalid' &
             is.na(df_meta$start_latitude),]) > 0){
      df_meta <- location.from.time(df_meta)
    }
    # add the survey id information
    df_meta <- survey.splitter(df_meta, threshold_days=14)
    # save metadata
    write.csv(df_meta, meta_dir, row.names=FALSE)
    print('Meta data saved to file.')
  # if not making meta, then load from source
  } else {
    # read from source
    df_meta <- read.csv(meta_dir)
    # format
    df_meta[,8:23] <- apply(df_meta[,8:23], 2, function(x) as.numeric(x))
    df_meta$cast_time_UTC <- as.POSIXct(df_meta$cast_time_UTC, tz='UTC')
    df_meta$cast_time_local <- as.POSIXct(df_meta$cast_time_local, tz='Australia/Sydney')
  }
  
  # drop invalid casts
  warning(paste(nrow(df_meta[df_meta$sample_type == 'Invalid',]),
                'invalid casts not read.'))
  df_read <- df_meta[!df_meta$sample_type == 'Invalid',]
  # drop 'manual' casts with no GPS data
  if (drop_no_location){
    warning(paste(nrow(df_read[is.na(df_read$start_latitude),]),
                  'casts without GPS data not read.'))
    df_read <- df_read[!is.na(df_read$start_latitude),]
  }
  # drop non transect
  if (drop_non_transect){
    warning(paste(nrow(df_read[is.na(df_read$transect_id),]),
                  'casts not on transect not read.'))
    df_read <- df_read[!is.na(df_read$transect_id),]
  }
  # remove duplicates (subsequent in time)
  if (drop_duplicates){
    before_cast_no <- nrow(df_read)
    df_read <- suppressWarnings(CTD.remove.dup.casts(df_read))
    warning(paste(before_cast_no - nrow(df_read),
                  'casts that were same place and close in time not read.'))
  }
  # drop unread casts from meta
  if (drop_unloaded_meta){
    df_meta <- df_read
  }
  # reset index
  row.names(df_meta) <- NULL
  # Now read in CTD casts
  # make empty list to populate
  cast_ls <- as.list(rep(NA, nrow(df_read)))
  # iterate through and read data
  i <- 1
  print(paste('Reading',nrow(df_read),'cast files..'))
  pb <- txtProgressBar(min = 0, max = nrow(df_read), style = 3)
  setTxtProgressBar(pb, i)
  for (fn in df_read$file_name){
    lines <- readLines(paste0(input_dir, fn, '.csv'))
    lines <- lines[29:length(lines)]
    lines <- strsplit(lines, ',')
    # to data frame
    df <- data.frame(do.call(rbind, lines[2:length(lines)]))
    colnames(df) <- c('pressure', 'depth', 'temperature', 'conductivity',
                      'specific_conductance', 'salinity', 'sound_velocity', 'density')
    # units = (decibar, meter, celcius, microSiemens/cm, microSiemens/cm, PSU, m/s, kg/m^3)
    # add id
    df$id <- fn
    df <- df[,c(9,1:8)]
    # add transect and survey ids
    df$transect_id <- df_read[df_read$file_name == df$id[1], 'transect_id'][1]
    df$survey_id <- df_read[df_read$file_name == df$id[1], 'survey_id'][1]
    # add to list
    cast_ls[[i]] <- df
    i <- i +1
    setTxtProgressBar(pb, i)
  }
  close(pb)
  # return cast_ls
  print('Merging data frames..')
  df <- do.call('rbind', cast_ls)
  # format data
  df$id <- as.factor(df$id)
  df[,2:9] = apply(df[,2:9], 2, function(x) as.numeric(x))
  # Add names
  output <- list(df, df_meta)
  names(output) <- c('data', 'meta')
  return(output)
}

################################
# Report manual loaction casts #
################################
report.manual.casts <- function(){
  # read in data
  dat <- cast.reader()
  CTD_df <- dat[[1]]
  meta <- dat[[2]]
  # get manual non invalid casts
  meta_manual <- meta[meta$location_source == 'Manual' & meta$sample_type != 'Invalid',]
  row.names(meta_manual) <- NULL
  # group by id and get max depth
  depthmax <- CTD_df[CTD_df$id %in% meta_manual$file_name,] %>%
    group_by(id) %>%
    summarise(depthmax = max(depth))
  # add to meta manual
  meta_manual$depthmax <- depthmax$depthmax
  # report time stamps and file ids
  print(paste(nrow(meta_manual),'casts with manual location source'))
  return(meta_manual[,c('cast_time_UTC', 'file_name', 'cast_duration', 'depthmax')])
}

########################
# Transect and Station #
########################
# cutoff is in meters
add.transect.station <- function(df_meta, cutoff=500, 
                                 transect_grid_dir='data/transects/transect_master.csv'){
  # read in transect grid
  transect_grid <- read.csv(transect_grid_dir)
  # transform into data frame
  transect_df <- data.frame(transect=rep(1:9, each=5),
                            station=rep(1:5, times=9),
                            lon=stack(transect_grid[,seq(2,ncol(transect_grid),2)])$values,
                            lat=stack(transect_grid[,seq(1,ncol(transect_grid),2)])$values)
  # make ID names
  transect_df$id <- paste0('T',transect_df$transect,'S',transect_df$station)
  transect_df <- transect_df[,c(5,1:4)]
  # drop non existing stations
  transect_df <- transect_df[complete.cases(transect_df),]
  # add column to collect data
  df_meta$transect_id <- NA
  # setup progress
  print('Determining nearest transect station to casts..')
  pb <- txtProgressBar(min = 0, max = nrow(df_meta), style = 3)
  # for each cast find the nearest station
  for(i in 1:nrow(df_meta)){
    row <- df_meta[i,]
    # calcuate distance to each station
    # cast start coords
    station_dists_start <- mapply(distHaversine,
                                  mapply(c, rep(row$start_longitude, nrow(transect_df)),
                                         rep(row$start_latitude, nrow(transect_df)), 
                                         SIMPLIFY = FALSE),
                                  mapply(c, transect_df$lon, transect_df$lat, 
                                         SIMPLIFY = FALSE))
    # cast end coords
    station_dists_end <- mapply(distHaversine,
                                mapply(c, rep(row$end_longitude, nrow(transect_df)),
                                       rep(row$end_latitude, nrow(transect_df)), 
                                       SIMPLIFY = FALSE),
                                mapply(c, transect_df$lon, transect_df$lat, 
                                       SIMPLIFY = FALSE))
    # combine
    # get closest station
    station_dists <- c(station_dists_start, station_dists_end)
    nearest_dist <- min.sw(station_dists, na.rm=TRUE)
    # if within 1km of a transect then add id
    if (!is.na(nearest_dist) & nearest_dist <= cutoff){
      min_idx <- which.min(station_dists)
      # if idx greater than length of data frame (i.e. end point closer) then
      # subtract length of start dist
      if (min_idx > length(station_dists_start)){
        min_idx <- min_idx - length(station_dists_start)
      }
      df_meta$transect_id[i] <- transect_df[min_idx, 'id']
    }
    # progress
    setTxtProgressBar(pb, i)
  }
  close(pb)
  # format 
  df_meta$transect_id <- as.factor(df_meta$transect_id)
  # return metadata
  return(df_meta)
}

#############################
# Add missing location data #
#############################
add.manual.locations <- function(df_meta, manual_dir='data/CTD/processed/meta/manual_locations.csv'){
  print('Adding manual location records to meta data..')
  # read in manual loactions
  manual_locs <- read.csv(manual_dir)
  # populate decimal degree columns
  # start_lon_deg_dec_min to start_lon_dec_deg
  manual_locs[!is.na(manual_locs$start_lon_deg_dec_min),]$start_lon_dec_deg <- 
    as.numeric(
      conv_unit(manual_locs[!is.na(manual_locs$start_lon_deg_dec_min),]$start_lon_deg_dec_min,
              from='deg_dec_min', to='dec_deg'))
  # start_lat_deg_dec_min to start_lat_dec_deg
  manual_locs[!is.na(manual_locs$start_lat_deg_dec_min),]$start_lat_dec_deg <- 
    as.numeric(
      conv_unit(manual_locs[!is.na(manual_locs$start_lat_deg_dec_min),]$start_lat_deg_dec_min,
                from='deg_dec_min', to='dec_deg'))
  # end_lon_deg_dec_min to end_lon_dec_deg
  manual_locs[!is.na(manual_locs$end_lon_deg_dec_min),]$end_lon_dec_deg <- 
    as.numeric(
      conv_unit(manual_locs[!is.na(manual_locs$end_lon_deg_dec_min),]$end_lon_deg_dec_min,
                from='deg_dec_min', to='dec_deg'))
  # end_lat_deg_dec_min to end_lat_dec_deg
  manual_locs[!is.na(manual_locs$end_lat_deg_dec_min),]$end_lat_dec_deg <- 
    as.numeric(
      conv_unit(manual_locs[!is.na(manual_locs$end_lat_deg_dec_min),]$end_lat_deg_dec_min,
                from='deg_dec_min', to='dec_deg'))
  # round off the columns
  manual_locs$start_lon_dec_deg <- round(as.numeric(manual_locs$start_lon_dec_deg), 4)
  manual_locs$start_lat_dec_deg <- round(as.numeric(manual_locs$start_lat_dec_deg), 5)
  manual_locs$end_lon_dec_deg <- round(as.numeric(manual_locs$end_lon_dec_deg), 4)
  manual_locs$end_lat_dec_deg <- round(as.numeric(manual_locs$end_lat_dec_deg), 5)
  # now feed data back into the meta data frame
  # first check if dataframes are ordered correctly
  if (!(sum(factor(df_meta[df_meta$file_name %in% manual_locs$file_name,'file_name']) 
            == manual_locs[,'file_name'])
      == nrow(manual_locs))){
    print('Data frame IDs do not align, investigate source code and data bases..')
    break
  }
  # insert data in meta data frame
  df_meta[which(df_meta$file_name %in% manual_locs$file_name),
          c('start_latitude', 'start_longitude', 'end_latitude','end_longitude')] <-
    manual_locs[,c('start_lat_dec_deg', 'start_lon_dec_deg', 
                   'end_lat_dec_deg', 'end_lon_dec_deg')]
  # return
  return(df_meta)
}

###############################################
# Add unknown CTD locations from time records #
###############################################
location.from.time <- function(df_meta, check=TRUE){
  # expand transect and station
  df_meta <- expand.transect.ID(df_meta)
  # get the points we need to assess
  no_coords_idx <- which(df_meta$location_source == 'Manual' &
                           df_meta$sample_type != 'Invalid' &
                           is.na(df_meta$start_latitude))
  print(paste(length(no_coords_idx),'casts have no coordinates.'))
  # give factors
  df_meta$missing_coords <- 'N'
  df_meta$missing_coords[no_coords_idx] <- 'Y'
  # for passing center values
  center_df <- df_meta
  center_df$lat <- center_df$start_latitude
  center_df$lon <- center_df$start_longitude
  # load basemap
  basemap <- load.basemap(mapset='montague', zoom=11, 
                          suffix='_transects', center=center_df)
  # load transect grid
  transect_df <- load.transect.grid()
  # print the points before and after
  for(i in 1:length(no_coords_idx)){
    # subset the dataframe
    idx <- no_coords_idx[i]
    df <- df_meta[(idx-3):(idx+3), ]
    cat(paste('Displaying cast',i,'of',length(no_coords_idx),'|',df_meta[idx,'file_name'],'\n\n'))
    # print the data
    print(df[, c('cast_time_UTC','sample_type', 'start_longitude', 'start_latitude', 
                 'transect_id','missing_coords')])
    # if no transect IDs skip this section
    if (is.na(df_meta$transect_id[idx-1]) | is.na(df_meta$transect_id[idx-1])){
      print('Cannot determine a possible transect ID..')
      invisible(readline(prompt="Press [enter] to continue"))
      next
    }
    # determine transect ID
    # if previous and next are the same likely also the same
    if (df_meta$transect_id[idx-1] == df_meta$transect_id[idx+1]){
      missing_id <- df_meta$transect_id[idx-1]
    # if in ascending sequence then will sit between
    } else if ((df_meta$transect[idx-1] == df_meta$transect[idx-1]) & 
               ((df_meta$station[idx-1] - df_meta$station[idx+1]) %in% c(-2, 2))){
      missing_station <- max(df_meta$station[idx-1], df_meta$station[idx+1]) - 1
      missing_id <- paste0('T',df_meta$transect[idx-1],'S',missing_station)
    }
    print(paste('Missing ID determined to be', missing_id))
    # check if right
    if (check){
      add2df <- askYesNo(paste0('Add into dataframe?'))
    } else {
      add2df <- TRUE
    }
    # add in data
    if (add2df){
      df_meta$transect_id[idx] <- missing_id
      df_meta[idx, c('start_longitude', 'start_latitude')] <- 
        transect_df[transect_df$id == missing_id, c('lon', 'lat')]
      print('Data added.')
    }
    # make a plot
    df <- df_meta[(idx-3):(idx+3), ]
    print(plot.CTD.missing(df, basemap))
    invisible(readline(prompt="Showing map. Press [enter] to continue"))
  }
  # fix formatting
  df_meta[,c('transect','station','missing_coords')] <- NULL
  return(df_meta)
}

###############################################
# Transect ID to transect and station columns #
###############################################
expand.transect.ID <- function(df_meta){
  df_meta$transect <- as.numeric(substr(df_meta$transect_id, 2, 2))
  df_meta$station <- as.numeric(substr(df_meta$transect_id, 4, 4))
  return(df_meta)
}

######################
# Load transect grid #
######################
load.transect.grid <- function(transect_grid_dir='data/transects/transect_master.csv'){
  # read in transect grid
  transect_grid <- read.csv(transect_grid_dir)
  # transform into data frame
  transect_df <- data.frame(transect=rep(1:9, each=5),
                            station=rep(1:5, times=9),
                            lon=stack(transect_grid[,seq(2,ncol(transect_grid),2)])$values,
                            lat=stack(transect_grid[,seq(1,ncol(transect_grid),2)])$values)
  # make ID names
  transect_df$id <- paste0('T',transect_df$transect,'S',transect_df$station)
  transect_df <- transect_df[,c(5,1:4)]
  # drop non existing stations
  transect_df <- transect_df[complete.cases(transect_df),]
  return(transect_df)
}

############################
# Generate survey metadata #
############################
# breaks CTD cast data into survey work periods
# reports how many complete surveys were done in each work period
survey.splitter <- function(df_meta, threshold_days=14, make_survey_time_file=FALSE, 
                            survey_time_fn='./data/transects/survey_times.csv'){
  print('Splitting casts into distinct surveys..')
  # calc diff time for casts
  df_meta <- cast.difftime(df_meta)
  # set threshold in seconds
  threshold_secs <- threshold_days*24*60*60
  # create survey ids
  # get year of casts
  df_meta$year <- year(df_meta$cast_time_UTC)
  # get index of splits
  split_idx <- which(df_meta$difftime > threshold_secs)
  # 2019 is a special case so will divide this into separate surveys using
  # a split on September 29
  split_idx <- append(split_idx,
                      which(df_meta$cast_time_UTC > ISOdatetime(2019,9,29,0,0,0))[1])
  # append start and end indexes to splits
  split_idx <- append(split_idx, 1, after=0)
  split_idx <- append(split_idx, nrow(df_meta), after=length(split_idx))
  # put in survey splits as factors
  df_meta$survey_year_no <- NA
  survey_no <- 1
  for (i in 2:length(split_idx)){
    df_meta[split_idx[i-1]:split_idx[i], 'survey_year_no'] <- survey_no
    survey_no <- survey_no + 1
  }
  # now rename based on how many factors in the year
  for (df in split(df_meta, df_meta$year)){
    df_index <- as.numeric(rownames(df))
    df_meta[df_index[1]:df_index[length(df_index)], 'survey_year_no'] <-
      as.numeric(factor(df_meta[df_index[1]:df_index[length(df_index)], 'survey_year_no']))
  }
  # append as character to id
  df_meta$survey_id <- as.factor(paste0(df_meta$year,'_S',df_meta$survey_year_no))
  
  # check threshold
  plot(df_meta$difftime, ylab='Time between casts (secs)', main='Cast time splitting threshold')
  abline(h=14*24*60*60, col='red')
  # check survey gorupings
  plot(df_meta$cast_time_UTC, col=as.factor(df_meta$survey_id), ylab='Datetime UTC',
       main='Survey_ID splitting')
  legend("bottomright", legend=levels(df_meta$survey_id), pch=16,
         col=unique(df_meta$survey_id), cex = 0.75)
  # plot survey times
  print(ggplot(df_meta[!is.na(df_meta$transect_id),]) +
    geom_point(aes(x=cast_time_UTC, y=transect_id, color=as.factor(substr(transect_id,1,2)))) +
    facet_wrap(~survey_id, scales='free') +
    labs(color = "Transect"))
  if (make_survey_time_file){
    survey_df <- data.frame(id=unique(df_meta$survey_id))
    survey_df$year <- as.numeric(substr(survey_df$id,1,4))
    survey_df$start_UTC <- as.POSIXct(sapply(split(df_meta, df_meta$survey_id), 
                                             function(df) min(df$cast_time_UTC)),
                                      origin='1970-01-01', tz='UTC')
    survey_df$end_UTC <- as.POSIXct(sapply(split(df_meta, df_meta$survey_id),
                                           function(df) max(df$cast_time_UTC)),
                                    origin='1970-01-01', tz='UTC')
    write.csv(survey_df, survey_time_fn, row.names=FALSE)
  }
  # clean uneeded columns
  df_meta[,c('year','difftime','survey_year_no')] <- NULL
  return(df_meta)
}

####################################
# Calculate time between CTD casts #
####################################
cast.difftime <- function(df_meta){
  # offset time
  df_meta$cast_time_UTC_lag <- data.table::shift(df_meta$cast_time_UTC)
  # calculate time difference
  df_meta$difftime <- mapply(difftime, df_meta$cast_time_UTC, 
                             df_meta$cast_time_UTC_lag, units='secs')
  # remove offsets
  df_meta$cast_time_UTC_lag <- NULL
  # return the track
  return(df_meta)
}

##########################################
# min and max wrappers to avoid warnings #
##########################################
# min
min.sw <- function(x, na.rm=FALSE){
  return(suppressWarnings(min(x, na.rm=na.rm)))
}
# max
max.sw <- function(x, na.rm=FALSE){
  return(suppressWarnings(max(x, na.rm=na.rm)))
}

#################################################
# Calculate repeats of values subset by another #
#################################################
calc.repeats <- function(df, child_col, parent_col){
  # reset row names
  row.names(df) <- NULL
  # make list to hold values
  repeat_col <- as.list(rep(NA, length(split(df, df[,parent_col]))))
  index_ls <- repeat_col
  # break into smaller data frames and iterate
  i <- 1
  for (subdf in split(df, df[,parent_col])){
    # work out repeats
    rle_result <- rle(as.vector(subdf[,child_col]))
    rep_ls <- as.list(rep(NA, length(rle_result$lengths)))
    for (j in 1:length(rle_result$lengths)){
      rep_ls[[j]] <- rep(j, rle_result$lengths[j])
    }
    # add to list
    repeat_col[[i]] <- unlist(rep_ls)
    index_ls[[i]] <- as.numeric(rownames(subdf))
    i <- i + 1
  }
  # order
  index_ls <- unlist(index_ls)
  output <- rep(NA, length(index_ls))
  output[index_ls] <- unlist(repeat_col)
  # return
  return(as.factor(output))
}

################################################
# Add in survey grid lat/lon to CTD data frame #
################################################
CTD.add.grid.coords <- function(CTD_df){
  # load grid data
  grid_df <- load.transect.grid()
  # add lats
  CTD_df$lat <- sapply(CTD_df$transect_id, function(t_id) 
    grid_df[grid_df$id==t_id, 'lat'])
  # add lons
  CTD_df$lon <- sapply(CTD_df$transect_id, function(t_id) 
    grid_df[grid_df$id==t_id, 'lon'])
  # return data frame
  return(CTD_df)
}

############################
# Cast duplication removal #
############################
# Cleans cast meta data of casts on transects that were done at
# the same time (last is better to keep as first few will be a dud.. probably)
# Use subsequent index numbers to know if drops were successive
# UPDATE - use time instead of indexes (index fails with invalid)
CTD.remove.dup.casts <- function(CTD_meta, method='time', t_threshold=35,
                                 verbose=FALSE){
  # only consider casts with transect IDs
  CTD_multidrops <- CTD_meta[!is.na(CTD_meta$transect),]
  # create survey/tranect/date column (date is local to avoid date changes)
  CTD_multidrops$drop_id <- paste0(CTD_multidrops$survey_id,'_',CTD_multidrops$transect_id,
                             '_',strftime(CTD_multidrops$cast_time_local, format='%Y%m%d'))
  # find drops ids with multiple casts
  multidrops <- names(table(CTD_multidrops$drop_id)[table(CTD_multidrops$drop_id) > 1])
  CTD_multidrops <- CTD_multidrops[CTD_multidrops$drop_id %in% multidrops,]
  # for each multidrop find drops that are within in N minutes of each opther
  # and if they are save to a drop list
  j <- 1
  droplist <- list()
  for (i in 1:length(multidrops)){
  # UPDATE - use time instead of indexes (index fails with invalid casts)
    # get the drop to work on
    thedrop <- multidrops[i]
    # make subdf of meta data
    CTD_sub <- CTD_multidrops[CTD_multidrops$drop_id == thedrop,]
    # group into drops based on index numbers
    # method index
    if (method == 'index'){
      CTD_sub$index_diff <- mapply(function(i1, i2) i2 - i1, 
                                   data.table::shift(as.numeric(row.names(CTD_sub))),
                                   as.numeric(row.names(CTD_sub)))
      # replace NAs
      CTD_sub$index_diff[is.na(CTD_sub$index_diff)] <- 1
      # cycle through the drops and create suffix based on index 
      # (i.e. are they subsequent drops?)
      drop_suffix <- 1
      for (idx in 1:nrow(CTD_sub)){
        # if not successive increase drop suffix
        if (!CTD_sub$index_diff[idx] <= 1){
          drop_suffix <- drop_suffix + 1
        }
        # add suffix
        CTD_sub$drop_id[idx] <- paste0(CTD_sub$drop_id[idx],'-',drop_suffix)
      }
    # method time
    } else if (method == 'time'){
      CTD_sub$index_diff <- mapply(function(t1, t2) difftime(t2, t1, units='mins'), 
                                   data.table::shift(CTD_sub$cast_time_local),
                                   CTD_sub$cast_time_local)
      # replace NAs
      CTD_sub$index_diff[is.na(CTD_sub$index_diff)] <- 0
      # cycle through the drops and create suffix based on difftime 
      # (i.e. are within time threshold)
      drop_suffix <- 1
      for (idx in 1:nrow(CTD_sub)){
        # if not successive increase drop suffix
        if (CTD_sub$index_diff[idx] > t_threshold){
          drop_suffix <- drop_suffix + 1
        }
        # add suffix
        CTD_sub$drop_id[idx] <- paste0(CTD_sub$drop_id[idx],'-',drop_suffix)
      }
    } else {
      cat(paste0('Error: Unknow method "',method,'"'))
      return(NA)
    }
    if (verbose){
      print(CTD_sub[,c('file_name', 'cast_time_local', 'index_diff', 'drop_id')])
    }  
    # only keep drop ids now with multiple values
    multidrops_subsequent <- names(table(CTD_sub$drop_id)[table(CTD_sub$drop_id) > 1])
    # now grab the file names if there are casts to drop
    if (length(multidrops_subsequent) > 0){
      drops <- as.character(CTD_sub$file_name[CTD_sub$drop_id %in% multidrops_subsequent])
      # if multiple remove the last filename (one we keep)
      if (length(drops) > 1){
        drops <- drops[1:(length(drops)-1)]
      }
      # add to the drop list
      droplist[[j]] <- drops
      j = j + 1
    }
  }
  # now remove these from the meta data and return
  droplist <- unlist(droplist)
  CTD_meta_out <- CTD_meta[!CTD_meta$file_name %in% droplist,]
  # return processed meta data
  return(CTD_meta_out)
}

########################
# Survey Cast Selector #
########################
# Manual selection of casts in survey
# saves data in file to be pulled later
CTD.survey.cast.selector <- function(CTD_meta){
  # check if survey set exists
  if (file.exists('./data/transects/survey_set.csv')){
    if (!askYesNo(paste0('"./data/transects/survey_set.csv" ',
                    'already exists. Would you like to overwrite?'))){
      return(NA)
    }
  }
  # split casts by survey
  CTD_meta_survey_ls <- split(CTD_meta, CTD_meta$survey_id)
  # make interactive plot
  # setup the color levels
  col_levels <- length(levels(CTD_meta$transect_id))
  color_set <- cmocean('phase')(col_levels)
  # cycle through sub dataframes
  survey_set <- list()
  i <- 1
  for (survey_df in CTD_meta_survey_ls){
    plt_title <- survey_df$survey_id[1]
    # make the dataframe
    df <- data.frame(x=survey_df$cast_time_UTC, y=as.numeric(survey_df$transect_id))
    df_color <- color_set[df$y]
    # make plot
    x11(width=12, height=7.5) # open graphics window
    # make plot
    plot(df, pch=16, col=df_color, main=plt_title, ylim = c(1,col_levels),
         panel.first = grid(), yaxt="n", ylab=NA, xlab=NA)
    # text y ticks
    axis(2, at=1:col_levels, labels=levels(CTD_meta$transect_id), las=1)
    # interactive plot to select points
    selectedPoints <- as.numeric(fhs(df, mark=TRUE))
    # make a second plot to save
    df_color[-selectedPoints] = '#a9a9a9'
    png(paste0('./data/transects/survey_set_plots/',plt_title,'.png'))
    plot(df, pch=16, col=df_color, main=plt_title, ylim = c(1,col_levels),
         panel.first = grid(), yaxt="n", ylab=NA, xlab=NA)
    axis(2, at=1:col_levels, labels=levels(CTD_meta$transect_id), las=1)
    dev.off()
    # make a survey set out of selected casts
    survey_df <- survey_df[selectedPoints,]
    survey_set[[i]] <- survey_df
    # permutate
    i <- i + 1
    # check if any duplicates and provide warning
    if (max(table(survey_df$transect_id)) >  1){
      warning(paste('Duplicate cast stations detected in survey', plt_title))
    }
  }
  # close remaining windows
  graphics.off()
  # export the set
  survey_set <- do.call('rbind', survey_set)
  survey_set <- survey_set[,c('survey_id', 'transect_id', 'file_name', 'cast_time_UTC')]
  # save as csv
  write.csv(survey_set, './data/transects/survey_set.csv', row.names = FALSE)
  print('Saved survey set to file: ./data/transects/survey_set.csv')
}

###########################
# Filter survey CTD casts #
###########################
# go through and filter the CTD drops based on the survey filter file
# i.e. the selected survey points
cast.survey.filter <- function(CTD_data, survey_set='./data/transects/survey_set.csv'){
  # read in survey set data
  survey_set <- read.csv(survey_set)
  # filter datsets for casts in survey set
  # filter dataframe
  CTD_data[[1]] <- CTD_data[[1]][CTD_data[[1]]$id %in% survey_set$file_name,]
  # filter meta data
  CTD_data[[2]] <- CTD_data[[2]][CTD_data[[2]]$file_name %in% survey_set$file_name,]
  # return results
  return(CTD_data)
}

#############################################
# Get all CTD cast indexes in kriging feild #
#############################################
krig.get.cast.index <- function(t.df, tol=3,
                                transect_csv='./data/transects/transect_master.csv'){
  # Get all CTD cast indexes in kriging feild
  # load transects
  transect_df <- read.csv(transect_csv)
  # get lat and lons column names
  transect_colnames <- colnames(transect_df)
  latcols <- transect_colnames[which(substr(transect_colnames,12,16) == 'lat')]
  loncols <- transect_colnames[which(substr(transect_colnames,12,16) == 'lon')]
  # melt dataset
  transect_coords <- melt(transect_df[,latcols])
  transect_coords$lons <- melt(transect_df[,loncols])$value
  transect_coords$lats <- transect_coords$value
  transect_coords <- transect_coords[,c('lons','lats')]
  transect_coords <- transect_coords[complete.cases(transect_coords),]
  transect_coords <- paste0(round(transect_coords$lons,tol),'_',round(transect_coords$lats,tol))
  # filter cast datasets (replace with high value)
  cast_idx <- which(unlist(lapply(paste0(round(t.df$lon,tol),'_',round(t.df$lat,tol)),
                                  function(x) x %in% transect_coords)))
  return(cast_idx)
}
###################
# CTD time filter #
###################
CTD.time.filter <- function(CTD_dat, start, end, tz='UTC'){
  start <- as.POSIXct(start, tz=tz)
  end <- as.POSIXct(end, tz=tz)
  # filter meta
  CTD_dat[[2]] <- CTD_dat[[2]][CTD_dat[[2]]$cast_time_UTC > start &
                               CTD_dat[[2]]$cast_time_UTC < end,]
  # filter data
  CTD_dat[[1]] <- CTD_dat[[1]][as.character(CTD_dat[[1]]$id) %in% 
                               as.character(CTD_dat[[2]]$file_name),]
  return(CTD_dat)
}

############################
# Contruct OCE CTD objects #
############################
# downcast argument speicies keeping only downcast component
ctd.2.oce <- function(CTD, trim=T, trim.method='downcast',
                      transect.sections=F){
  oce <- list()
  # loop all casts
  message('Converting casts to oce format..')
  pb <- txtProgressBar(min = 0, max = nrow(CTD$meta), style = 3)
  for (i in 1:nrow(CTD$meta)){
    # get meta
    meta <- CTD$meta[i,]
    # get data subset
    data <- CTD$data[CTD$data$id == meta$file_name,]
    
    # Create oce object
    oce[[i]] <- as.ctd(salinity = data$salinity,
                       temperature = data$temperature,
                       pressure =  data$pressure,
                       conductivity = data$conductivity,
                       station = meta$transect_id,
                       startTime = meta$cast_time_UTC,
                       longitude = meta$start_longitude,
                       latitude = meta$start_latitude,
                       deploymentType = 'profile',
                       sampleInterval = meta$samples_per_second)
    
    # Add additional information
    oce[[i]]@data[["castAway_depth"]] <- data$depth
    oce[[i]]@metadata[["survey"]] <- meta$survey_id
    oce[[i]]@metadata[["id"]] <- meta$file_name
    
    #Progress
    setTxtProgressBar(pb, i)
  }
  
  # trim upcast and only keep downcast
  if (trim){
    oce <- lapply(oce, function(d) ctdTrim(d, method=trim.method))
  }
  
  # return data
  return(oce)
}

##################################
# OCE - Create transect sections #
##################################
oce.transect.sections <- function(oce, tz='Australia/Sydney', minN=3){
  # Extract transect and survey data
  survey.meta <- data.frame(id=sapply(oce, function(ctd) ctd@metadata$id),
                            survey=sapply(oce, function(ctd) ctd@metadata$survey),
                            station=sapply(oce, function(ctd) ctd@metadata$station),
                            date=sapply(oce, function(ctd) ctd@metadata$startTime),
                            idx=1:length(oce))
  survey.meta$date <- as.POSIXct(survey.meta$date, origin='1970-01-01', tz='UTC')
  attr(survey.meta$date, "tzone") <- tz
  survey.meta$date <- as.Date(survey.meta$date)
  # Group and make into sections
  sections <- split(survey.meta, paste0(survey.meta$survey,'_',
                                        substr(survey.meta$station,1,2),'_',
                                        survey.meta$date))
  # Drop secitons with less than N stations
  sections <- sections[sapply(sections, function(d) nrow(d) >= minN)]
  
  # Create the sections
  for (i in 1:length(sections)){
    id <- paste(sections[[i]]$survey[1],
                substr(sections[[i]]$station[1],1,2))
    sections[[i]] <- as.section(oce[sections[[i]]$idx])
    sections[[i]]@metadata$sectionId <- id
  }
  
  return(sections)
}

##############################
# CTD plot of missing values #
##############################
# used in location.from.time in CTD_functions
plot.CTD.missing <- function(df, basemap){
  transect_map <- basemap +
    geom_path(data=df, aes(x=start_longitude, y=start_latitude)) +
    geom_point(data=df, aes(x=start_longitude, y=start_latitude, color=missing_coords))
  return(transect_map)
}

#############################
# Plot subsequent CTD drops #
#############################
CTD.plot.station.repeats <- function(CTD_df, CTD_meta, split.year=F, 
                                     colour=c('datetime','repeats')){
  # match arg
  colour <- match.arg(colour)
  # create survey/station id
  if (split.year){
    CTD_meta$survey_station_id <- as.factor(paste0(CTD_meta$survey_id,'_',CTD_meta$transect_id))
  } else {
    CTD_meta$survey_station_id <- as.factor(CTD_meta$transect_id)
  }
  
  # drop all levels that do not have multiple repeats
  cast_multis <- table(CTD_meta$survey_station_id)
  cast_multis <- names(cast_multis[cast_multis > 1])
  CTD_meta <- CTD_meta[CTD_meta$survey_station_id %in% cast_multis,]
  # filter the main data frame
  CTD_df <- CTD_df[CTD_df$id %in% CTD_meta$file_name,]
  if (split.year){
    CTD_df$survey_station_id <- as.factor(paste0(CTD_df$survey_id,'_',CTD_df$transect_id))
  } else {
    CTD_df$survey_station_id <- as.factor(CTD_df$transect_id)
  }
  # work out repeat num
  CTD_df$repeat_no <- calc.repeats(CTD_df, 'id', 'survey_station_id')
  # if colour datetime replace repeat no with datetime
  if (colour == 'datetime'){
    CTD_df$repeat_no <- unlist(lapply(CTD_df$id, function(id)
      CTD_meta$cast_time_UTC[CTD_meta$file_name == id]))
  }
  # make the plot
  CTD_repeat_plot <- ggplot(data=CTD_df, aes(x=temperature, y=depth, color=repeat_no)) +
    geom_point(size=0.5) + scale_y_continuous(trans="reverse") +
    facet_wrap(~survey_station_id)
  return(CTD_repeat_plot)
}

#############################
# Plot subsequent CTD drops #
#############################
CTD.plot.all.casts <- function(CTD_df, survey_filter=NA, data='temperature',
                               facet_surveys=FALSE, colour_surveys=FALSE){
  # filter by survey
  if (!is.na(survey_filter)){
    CTD_df <- CTD_df[CTD_df$survey_id %in% survey_filter,]
  }
  # make the plot
  if (!colour_surveys){
    CTD_cast_plot <- ggplot(data=CTD_df, aes(x=CTD_df[,data], y=depth, color=id)) +
      geom_path(size=0.5) + scale_y_continuous(trans="reverse") +
      theme(legend.position='none')
  } else {
    CTD_cast_plot <- ggplot(data=CTD_df, aes(x=CTD_df[,data], y=depth, color=survey_id)) +
      geom_path(size=0.5) + scale_y_continuous(trans="reverse")
    #theme(legend.position='none')
  }
  # facet surveys?
  if (facet_surveys){
    CTD_cast_plot <- CTD_cast_plot + facet_wrap(~survey_id)
  }
  return(CTD_cast_plot)
}

########################
# Plot single CTD cast #
########################
CTD.plot.cast <- function(CTD, castID='random', survey_id=NULL, transect_id=NULL){
  data <- CTD$data
  meta <- CTD$meta
  # Pre filter data
  if (!is.null(survey_id)){
    data <- data[data$survey_id == survey_id,]
    meta <- meta[meta$survey_id == survey_id,]
  }
  if (!is.null(transect_id)){
    data <- data[data$transect_id == transect_id,]
    meta <- meta[meta$transect_id == transect_id,]
  }
  # If ID is random return random ID
  if (castID == 'random'){
    castID <- as.character(sample(unique(data$id), 1))
  }
  
  # Subset data
  data <- data[data$id == castID,]
  
  # Rescale salinity
  data$salt_rescale <- scales::rescale(data$salinity,
                                       to=c(min(data$temperature), 
                                            max(data$temperature)))
  # make title
  title <- paste0(meta$file_name[meta$file_name == castID],' \n',
                  meta$transect_id[meta$file_name == castID], 
                  strftime(meta$cast_time_UTC[meta$file_name == castID], ' | %Y-%m-%d'))
  
  print(castID)
  
  # Set colours
  colors <- c("T" = "darkred", "S" = "darkgreen")
  
  # make the plot
  p <- ggplot(data=data, aes(x=depth)) +
    geom_path(aes(y=temperature,  color='T'), size=1, alpha=.8) +
    geom_point(aes(y=temperature), size=1, alpha=.8, col='darkred') +
    geom_path(aes(y=salt_rescale, color='S'), size=1, alpha=.8) +
    geom_point(aes(y=salt_rescale), size=1, alpha=.8, col='darkgreen') +
    ylab('Temperature (\u00B0C)') + xlab('Depth (m)') +
    scale_color_manual(values = colors) + labs(color=NULL) +
    scale_y_continuous(sec.axis = sec_axis(~scales::rescale(.,
                                                            to=c(min(data$salinity), max(data$salinity))),
                                           name="Salinity (PSU)")) +
    coord_flip() +
    scale_x_reverse() +
    ggtitle(title)
  
  return(p)
}


########################
# Plot CTD survey path #
########################
plot.CTD.survey.path <- function(CTD_meta, split_dates=TRUE, color_time=FALSE){
  # load base map
  basemap <- load.basemap('montague', 11, suffix='transect', center=c(150.2,-36.325))
  # make map
  path_plot <- basemap +
    scale_y_continuous(limits = c(min(CTD_meta$start_latitude, na.rm=TRUE)-0.04, 
                                  max(CTD_meta$start_latitude, na.rm=TRUE)+0.04), expand = c(0, 0)) +
    scale_x_continuous(limits = c(min(CTD_meta$start_longitude, na.rm=TRUE)-0.09, 
                                  max(CTD_meta$start_longitude, na.rm=TRUE)+0.09), expand = c(0, 0))
  if (color_time){
    path_plot <- path_plot + 
      geom_path(data=CTD_meta, mapping=aes(x=start_longitude, y=start_latitude, color=cast_time_UTC),
                arrow=arrow(angle=30,length=unit(0.1,"inches"),type="open")) +
      geom_point(data=CTD_meta, mapping=aes(x=start_longitude, y=start_latitude, color=cast_time_UTC),
                 arrow=arrow(angle=30,length=unit(0.1,"inches"),type="open"))
  } else {
    path_plot <- path_plot + 
      geom_path(data=CTD_meta, mapping=aes(x=start_longitude, y=start_latitude),
                color='orange', arrow=arrow(angle=30,length=unit(0.1,"inches"),type="open")) +
      geom_point(data=CTD_meta, mapping=aes(x=start_longitude, y=start_latitude),
                 color='orange', arrow=arrow(angle=30,length=unit(0.1,"inches"),type="open"))
  }
  if (split_dates){
    path_plot <- path_plot + facet_wrap(~paste(date(cast_time_local),'(AEST)'))
  }
  return(path_plot)
}
