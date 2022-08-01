source('R/CTD.R')
source('R/eudyptula.R')

#------------------------------#
# Normalise data between range #
#------------------------------#
normalise <- function(v, a=0, b=1){
  norm <- a + ((v - min(v, na.rm=T)*(b - a))/(max(v) - min(v)))
  return(norm)
}

#----------------------#
# Nomralise data frame #
#----------------------#
# For each column in a data frame, normalise if numeric
normalise.df <- function(df, a=0, b=1){
  for (var in colnames(df)){
    if (is.numeric(df[,var])){
      df[,var] <- normalise(df[,var], a=a, b=b)
    }
  }
  return(df)
}

#-----------------------------------#
# Load agggretaion data and process #
#-----------------------------------#
load.agg <- function(fn='./data/surveys/acoustics/agg.rds', raw=F, verbose=F){
  agg <- readRDS(fn)
  if (!raw){
    agg <- agg.clean(agg, verbose=verbose)
  }
  return(agg)
}

#------------------------#
# Clean aggregation data #
#------------------------#
agg.clean <- function(agg, verbose=F, 
                      transect_lines_fn='./data/surveys/meta/transect_lines.csv'){
  # Clean out junk values
  agg <- agg[agg$NASC > 0,] # only NASC > than 0
  agg <- agg[agg$Lon_M <= 360,] # no error lons
  agg <- agg[agg$Lat_M < 0, ] # no northern hemisphere
  agg <- agg[agg$Corrected_area > 0,] # no negative 
  agg <- agg[agg$Corrected_thickness  > 0,] # no negative 
  agg <- agg[agg$Corrected_length > 0,] # no negative 
  agg <- agg[agg$Sv_mean < 0,] # no positive
  agg <- agg[agg$Sv_mean > -900,] # no missing
  
  # add datetimes
  agg$dtUTC <- as.POSIXct(paste0(agg$Date_M,' ',agg$Time_M), format="%Y%m%d %H:%M:%OS", tz='UTC')
  
  # Add cross shore distance
  # Find nearest transect for each point
  transect_lines <- read.csv(transect_lines_fn)
  agg$transect <- sapply(agg$Lat_M, function(lat) transect_lines$transect[which.min(abs(transect_lines$lat1 - lat))])
  # Calculate crossshore distance using haversin formula
  points <- transect_lines[agg$transect, c('lon1', 'lat1')]
  points$lon2 <- agg$Lon_M
  agg$crossShoreDist <- mapply(function(lon1, lon2, lat) distHaversine(c(lon1, lat), c(lon2, lat)),
                               points$lon1, points$lon2, points$lat1)
  # convert to km
  agg$crossShoreDist <- agg$crossShoreDist/1000
  
  # Outliers
  agg <- outlier.rm(agg, c("NASC", "Sv_mean", "Corrected_area"), verbose=verbose)
  
  return(agg)
}

#----------------------------------------#
# Orutlier removal for aggregatioon data #
#----------------------------------------#
outlier.rm <- function(df, variables, threshold=.99, verbose=T){
  rm.idx.all <- c()
  for (var in variables){
    if (verbose){
      plot(df[,var], main=var)
      abline(h=quantile(df[,var], threshold), col='red')
    }
    rm.idx <- which(df[,var] > quantile(df[,var], threshold))
    rm.idx.all <- c(rm.idx, rm.idx.all)
    if (verbose){
      message(var,' outlier threshold set to < ',quantile(df[,var], threshold))
    }
  }
  rm.idx.all <- rm.idx.all[!duplicated(rm.idx.all)]
  if (verbose){
    message(length(rm.idx.all),' rows removed for in total')
  }
  df <- df[-rm.idx.all,]
  return(df)
}

#------------------------------------------#
# Make CTD/agg data frame of each CTD cast #
#------------------------------------------#
# CTD_area_agg filters agg data so only agg data in CTD area is considered
# max.depth reduces all data to only consider swarms in top N m (use fule foor penguins)
CTD.agg <- function(CTD, agg, filter=TRUE, CTD_area_agg=FALSE, max.depth=NULL,
                    station_fn='./data/surveys/meta/transect_df.csv'){
  # Load transect station data
  stations <- read.csv(station_fn) 
  # make list of first and last station IDs
  end_points <- split(stations, stations$transect)
  end_points <- lapply(end_points, function(df) df[c(1, nrow(df)),])
  end_points <- do.call(rbind, end_points)
  end_points <- end_points$id
  
  # filter values
  if (!is.null(max.depth)){
    agg <- agg[agg$Depth_mean < max.depth,]
  }
  
  cast.df <- data.frame(survey=rep(NA, nrow(CTD$meta)))
  message('Extracting aggregation data for casts..')
  pb <- txtProgressBar(min = 0, max = nrow(CTD$meta), style = 3)
  for (i in 1:nrow(CTD$meta)){
    cast <- CTD$meta[i,]
    cast.data <- CTD$data[CTD$data$id == cast$file_name,]
    # To search efficiently, first filter survey, then time (2 hours) and then GPS range
    # Survey
    cast.agg <- agg[agg$survey == cast$survey_id,]
    if (nrow(cast.agg) == 0) next
    cast.agg <- cast.agg[cast.agg$dt > (cast$cast_time_UTC - hours(2)) &
                           cast.agg$dt < (cast$cast_time_UTC + hours(2)),]
    # Calc distance to cast
    if (nrow(cast.agg) == 0) next
    cast.agg$distance <- mapply(function(lon1, lat1, lon2, lat2) 
      distHaversine(c(lon1,lat1), c(lon2,lat2)),
      cast.agg$Lon_M, cast.agg$Lat_M, 
      cast$start_longitude, cast$start_latitude)
    # Filter distance 
    cast.agg <- cast.agg[cast.agg$distance <= 1875,]
    setTxtProgressBar(pb, i)
    if (nrow(cast.agg) == 0) next
    # Populate data frame
    # Cast meta
    cast.df[i, c('survey', 'cast_id', 'station_id', 'cast_time_UTC')] <- 
      cast[, c('survey_id', 'file_name', 'transect_id', 'cast_time_UTC')]
    cast.df[i, 'max_depth'] <- ceiling(max(cast.data$depth))
    # filter depth to CTD area if needed
    if (CTD_area_agg){
      cast.agg <- cast.agg[cast.agg$Depth_mean < ceiling(cast.df[i, 'max_depth']),]
    }
    if (nrow(cast.agg) == 0) next
    cast.df[i, 'swarm_count'] <- nrow(cast.agg)
    cast.df[i, 'swarm_mean_depth'] <- mean(cast.agg$Depth_mean, na.rm=T)
    # Cast data
    cast.df[i, c('temp_min', 'temp_mean', 'temp_max', 'salt_min', 'salt_mean', 'salt_max')] <-
      c(min(cast.data$temperature), mean(cast.data$temperature), max(cast.data$temperature),
        min(cast.data$salinity), mean(cast.data$salinity), max(cast.data$salinity))
    # Aggregation data
    cast.df[i, c('NASC_min', 'NASC_mean', 'NASC_max','Sv_mean_mean', 'Sv_mean_sum', 
                 'corrected_area_mean', 'corrected_area_sum')] <-
      c(min(cast.agg$NASC), mean(cast.agg$NASC), max(cast.agg$NASC), 
        Sv_mean.log(mean(Sv_mean.linear(cast.agg$Sv_mean))),
        Sv_mean.log(sum(Sv_mean.linear(cast.agg$Sv_mean))),
        mean(cast.agg$Corrected_area), sum(cast.agg$Corrected_area))
    # Calculate CTD area
    xdist <- ifelse(cast$transect_id %in% end_points, 1.875, 1.875*2)
    cast.df[i, 'CTD_area'] <- ceiling(max(cast.data$depth))*(xdist*1000)
  }
  close(pb)
  
  if (filter){
    cast.df <- cast.df[complete.cases(cast.df),]
  }
  
  # Calculate per unit area
  # change to km2
  cast.df$CTD_area <- cast.df$CTD_area/1000
  # swarms per unit area
  cast.df$area_swarm_count <- cast.df$swarm_count / cast.df$CTD_area
  # corrected area per unit area
  cast.df$area_corrected_area_sum <- (cast.df$corrected_area_sum/1000) / cast.df$CTD_area
  # Sv mean per unit area
  cast.df$area_Sv_mean_sum <- Sv_mean.log(Sv_mean.linear(cast.df$Sv_mean_sum)/cast.df$CTD_area)
  
  # Add cross shoroe distance
  cast.df$crossShoreDistance <- (as.numeric(substr(cast.df$station_id, 4, 4)) - 1) * 3.75
  
  return(cast.df)
}



