#  ROMS processing tools
library(stringr)
library(readr)
library(ncdf4)
library(lubridate)
library(terra)
library(spacetime)
library(ggplot2)
library(metR)
library(cmocean)
library(RANN)
source('./R/utilities.R')

##################
# Load ROMS data #
##################
# Load a data period/range
# Later add a way to extract depth
load.ROMS <- function(var='temp',
                      period=1994:2016,
                      downsample.factor=NULL,
                      ROMS.version=c('highres','surface'),
                      ROMS.dir='/Volumes/LP_MstrData/master-data/ocean/ROMS',
                      UTM=F,
                      season.rm=F,
                      output=c('dataframe', 'STFDF')){
  # Match args
  output <- match.arg(output)
  ROMS.version <- match.arg(ROMS.version)
  # Suppress Terra progress
  terraOptions(progress=0)
  # Set the direcotry path
  if (ROMS.version == 'highres'){
    ROMS.dir <- paste0(ROMS.dir,'/highres')
    s_rho <- 30
    ocean_count <- -1
  } else if (ROMS.version == 'surface'){
    ROMS.dir <- paste0(ROMS.dir,'/surface_subset/merged')
    message('Warning: surface subset is a single file and will take',
            ' a while to load into memory.\nPlease wait..')
    s_rho <- 1
    ocean_count <- 1
  }
  # Load the file list
  file.ls <- list.files(ROMS.dir, pattern='*.nc', full.name = T)
  if (ROMS.version == 'highres') {
    fn <- parse_number(as.vector(sapply(file.ls, function(x) strsplit(x, '/')[[1]][8])))
    # only keep files in period
    file.ls <- file.ls[which(fn %in% period)]
  }

  # First load in positional data
  message('Opening NetCDF metadata..')
  ncin <- nc_open(file.ls[1])
  nc.lons <- ncvar_get(ncin, 'lon_rho')
  nc.lats <- ncvar_get(ncin, 'lat_rho')
  dims <- c(nrow(nc.lons), ncol(nc.lons))
  # downsample if needed
  if (!is.null(downsample.factor)){
    message('Downsampling by a factor of ',downsample.factor)
    nc.lons <- terra::aggregate(
      rast(nc.lons), fact=downsample.factor, na.rm=T)
    nc.lats <- terra::aggregate(
      rast(nc.lats), fact=downsample.factor, na.rm=T)
    dims <- c(nrow(nc.lons), ncol(nc.lons))
    nc.lons <- terra::as.array(nc.lons)
    nc.lats <- terra::as.array(nc.lats)
  }
  nc.lons <- as.vector(nc.lons)
  nc.lats <- as.vector(nc.lats)
  
  # Lists to catch data
  nc.var=list()
  nc.time=list()
  nc.file=list()
  
  # Get time early if using surface subset
  if (ROMS.version == 'surface'){
    nc.time <- as.vector(ncvar_get(ncin, 'ocean_time'))
    nc.time <- as.POSIXct(nc.time, tz='UTC', origin='1990-01-01 00:00:00')
    # calculate time indexes in requested period
    time.index <- which(year(nc.time) %in% period)
    nc.time <- nc.time[time.index]
    file.ls <- rep(file.ls[1], length(time.index))
    j <- time.index[1]
  }
  
  # Loop each file and extract data
  message('Loading NETCDF data..')
  pb <- txtProgressBar(min = 0, max = length(file.ls), style = 3)
  if (ROMS.version == 'highres'){
    j <- 1
  }
  # Loop the file list
  # For surface subset loop the file for each time stamp and
  # use permutaed indexes to extract correct time
  for (i in 1:length(file.ls)){ # Testing on a few files for now
    if (ROMS.version == 'highres'){
      ncin <- nc_open(file.ls[i])
    }
    setTxtProgressBar(pb, i)
    # Get data (xi, eta, s, time)
    nc.var[[i]] <- ncvar_get(ncin, var, start=c(1,1,s_rho,j), count=c(-1,-1,1,ocean_count))
   
    if (ROMS.version == 'highres'){
      # Get time if using highres
      nc.time[[i]] <- as.vector(ncvar_get(ncin, 'ocean_time'))
      # Convert time to datetime (seconds since 1990-01-01 00:00:00)
      nc.time[[i]]  <- as.POSIXct(nc.time[[i]], tz='UTC', origin='1990-01-01 00:00:00')
      time.len <- length(nc.time[[i]])
      # Find timestamps not in current year (start or end) and drop
      current.year <- year(nc.time[[i]][183])
      drop.start <- year(nc.time[[i]][1]) != current.year
      drop.end <- year(nc.time[[i]][time.len]) != current.year
      if (drop.start){
        nc.time[[i]] <- nc.time[[i]][-1]
        nc.var[[i]] <- nc.var[[i]][,,-1]
      }
      if (drop.end){
        nc.time[[i]] <- nc.time[[i]][-time.len]
        nc.var[[i]] <- nc.var[[i]][,,-time.len]
      }
    }
    # Filter to polygon points
    # nc.var[[i]] <- apply(nc.var[[i]], 3, function(x) x[poly.idx])
    # day.count <- ncol(nc.var[[i]])
    # Downsample the matrix spatially
    if (!is.null(downsample.factor)){
      # Use terra raster processing to downsample
      nc.var[[i]] <- terra::as.array(
        terra::aggregate(rast(nc.var[[i]]), fact=downsample.factor, na.rm=T))
    }
    
    # Convert var to vector
    nc.var[[i]] <- as.vector(nc.var[[i]])
    
    # Permutate time if surface subset
    if (ROMS.version == 'surface'){
      j <- j + 1
    }

  }
  close(pb)
  
  # Unravel data
  day_count = sum(sapply(nc.time, length))
  coord_length = length(nc.lats)
  nc.var <- do.call(c, nc.var)

  # Convert to UTM
  # Convert nc.lon and nc.lay here instead
  if (UTM){
    message('Converting to UTM coordinates..')
    conversion.df <- lonlat2UTM.df(data.frame(lon=nc.lons, lat=nc.lats))
    nc.lons <- conversion.df$lon
    nc.lats <- conversion.df$lat
  }
  
  # Convert time so it can be easily unloaded
  if (ROMS.version == 'surface'){
    nc.time <- list(nc.time)
  }

  # Create a dataframe of the data
  if (output == 'dataframe' | season.rm){
    message('Converting extracted data to data frame..')
    ROMS <- data.frame(dt=rep(do.call(c, nc.time), each=coord_length),
                       lon=rep(nc.lons, day_count), 
                       lat=rep(nc.lats, day_count), 
                       var=nc.var)
  }
  
  if (season.rm){
    ROMS <- remove.seasonal.cycle(ROMS, var, UTM, show.plot=T)
    nc.var <- ROMS$var
  }
  
  if (output == 'STFDF'){
    message('Converting extracted data to spacetime STFDF object..')
    points <- data.frame(lon=nc.lons, lat=nc.lats)
    coordinates(points) <- c('lon','lat')
    points <- SpatialPoints(points)
    ROMS <- STFDF(points, do.call(c, nc.time), data.frame(var = nc.var))
  }
  
  # Return Terra progress to default
  terraOptions(progress=3)
  
  # Close connection
  nc_close(ncin)
  
  return(ROMS)
}

########################################
# Remove seasonal cycle from ROMS data #
########################################
remove.seasonal.cycle <- function(ROMS, var=c('temp', 'salt'), UTM=F, show.plot=T,
                                  CARS.dir='/Volumes/LP_MstrData/master-data/ocean/CARS'){
  # Match args
  var <- match.arg(var)
  # Load the correct CARS model
  fn  <- paste0(CARS.dir,'/', ifelse(var=='temp',
                                     'temperature_cars2009a.nc',
                                     'salinity_cars2009a.nc'))
  var.long <- ifelse(var=='temp','temperature','salinity')
  message('Loading CARS ',var.long,' data..')
  ncin <- nc_open(fn)
  # Get annual singal data (lon, lat, depth)
  nc.an_cos <- ncvar_get(ncin, 'an_cos', start=c(1,1,1), count=c(-1,-1,1))
  nc.an_sin <- ncvar_get(ncin, 'an_sin', start=c(1,1,1), count=c(-1,-1,1))
  nc.sa_cos <- ncvar_get(ncin, 'sa_cos', start=c(1,1,1), count=c(-1,-1,1))
  nc.sa_sin <- ncvar_get(ncin, 'sa_sin', start=c(1,1,1), count=c(-1,-1,1))
  # Extend lon lat into 2d arrays
  nc.lons <- rep(ncvar_get(ncin, 'lon'), ncol(nc.an_cos))
  nc.lats <- rep(ncvar_get(ncin, 'lat'), each=nrow(nc.an_cos))
  
  # close connection to disk
  nc_close(ncin)
  
  # Make dataframe
  CARS <- data.frame(lon=nc.lons, 
                     lat=nc.lats, 
                     doy=1,
                     an_cos=as.vector(nc.an_cos),
                     an_sin=as.vector(nc.an_sin),
                     sa_cos=as.vector(nc.sa_cos),
                     sa_sin=as.vector(nc.sa_sin))
  # Prefilter range to our area (avoid UTM conversion failure)
  CARS <- CARS[CARS$lon > 145 & CARS$lon < 165 & CARS$lat > -45 & CARS$lat < -25,]
  # Filter missing values
  CARS <- CARS[complete.cases(CARS),]
  # Convert to UTM coordinates
  if (UTM){
    CARS <- lonlat2UTM.df(CARS)
  }
  # Create a lon*lat*time dataframe with processed annual cycles
  CARS <- replicate(366, CARS, simplify=F)
  for (i in 1:length(CARS)){
    CARS[[i]]$doy <- i
  }
  CARS <- do.call(rbind, CARS)
  # Calculate signal
  CARS$t <- (2*pi)*(CARS$doy/366)
  CARS$signal <- CARS$an_cos*cos(CARS$t) + CARS$an_sin*sin(CARS$t) + 
                 CARS$sa_cos*cos(2*CARS$t) + CARS$sa_sin*sin(2*CARS$t)  
  
  # __Process ROMS__
  message('Removing seasonal signals from ROMS data..')
  # Add calculated day of year to ROMS dataset
  ROMS$doy <- as.numeric(strftime(ROMS$dt, format = "%j"))
  # Use KD tree to find vaulue to subtract from each value
  idx <- as.vector(nn2(CARS[,c('lon','lat','doy')], 
                       ROMS[,c('lon','lat','doy')], 
                       k=1)$nn.idx)
  
  # if plot change grab some before values
  if (show.plot){
    # get a value that is not NaN
    i <- sample.int(nrow(ROMS), 1)
    x <- ROMS$var[i]
    while (is.nan(x)){
      i <- sample.int(nrow(ROMS), 1)
      x <- ROMS$var[i]
    }
    idx.plotting <- ROMS$lon == ROMS$lon[i] & ROMS$lat == ROMS$lat[i]
    x.before <- ROMS$var[idx.plotting]
  }
  
  # extract sigal and subtract from ROMS value
  ROMS$var <- ROMS$var - CARS[idx,'signal']
  
  # Remove added columns
  ROMS$doy <- NULL
  
  # plot change
  if (show.plot){
    x.after <- ROMS$var[idx.plotting]
    # Make plotting data frame
    x.before <- data.frame(val=x.before, signal='before', x=1:length(x.before))
    x.after <- data.frame(val=x.after, signal='after', x=1:length(x.after))
    plot.df <- rbind(x.before, x.after)
    
    # plot
    p <- ggplot(plot.df, aes(x=x, y=val, color=signal)) +
      geom_point(alpha=.1) +
      geom_smooth()
    print(p)
  }
  
  # Return ROMS
  return(ROMS)
}


#####################
# Analyse variogram #
#####################
plot.vari <- function(vari, UTM.km=T, 
                      style=c('heatmap','contour','persp'),
                      theta=85, phi=25, expand=0.4, cmap='thermal'){
  # Match args
  style <- match.arg(style)
  # Convert UTM space to km
  xlab <- 'Distance (m)'
  if (UTM.km){
    vari$dist <- vari$dist/1000
    vari$spacelag <- vari$spacelag/1000
    vari$avgDist <- vari$avgDist/1000
    xlab <- 'Distance (km)'
  }
  # Heatmap plot (standard method)
  if (style == 'heatmap'){
    p <- plot(vari, xlab=xlab, ylab='Time lag (days)')
  # Contour plot
  } else if (style == 'contour'){
    p <- ggplot(vari, aes(x=avgDist, y=timelag)) + 
      geom_tile(aes(fill=gamma), width=51) + 
      geom_contour(aes(z=gamma), colour="black", size=0.4, alpha=0.4) + 
      geom_text_contour(aes(z=gamma),  colour="black" ) +
      labs(x = "Space lag (km)", 
           y = "Time lag (days)", 
           fill = "Variability") + 
      theme(panel.background = element_blank()) + 
      scale_fill_gradientn(colours=cmocean('haline')(100)) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) 
  } else if (style == 'persp'){
      tab <- xtabs(gamma ~ avgDist + timelag, data = vari)
      # downsample time (later add if greater than 150 or something)
      tab <- tab[, seq(1, ncol(tab), by=5)]
      z.facet.center <- (tab[-1, -1] + tab[-1, -ncol(tab)] + tab[-nrow(tab), -1] + 
                           tab[-nrow(tab), -ncol(tab)])/4
      # Range of the facet center on a 100-scale (number of colors)
      z.facet.range <- cut(z.facet.center, 100)
      p <- persp(tab, theta=theta, phi=phi, expand=expand, xlab='Space lag (km)', 
                 ylab='Time lag (days)',zlab='Variability', 
                 col=cmocean(cmap)(100)[z.facet.range])
  }
  return(p)
}

################
# FFT Analysis #
################
nff = function(x = NULL, n = NULL, up = 10L, plot = TRUE, add = FALSE, main = NULL, ...){
  #The direct transformation
  #The first frequency is DC, the rest are duplicated
  dff = fft(x)
  #The time
  t = seq(from = 1, to = length(x))
  #Upsampled time
  nt = seq(from = 1, to = length(x)+1-1/up, by = 1/up)
  #New spectrum
  ndff = array(data = 0, dim = c(length(nt), 1L))
  ndff[1] = dff[1] #Always, it's the DC component
  if(n != 0){
    ndff[2:(n+1)] = dff[2:(n+1)] #The positive frequencies always come first
    #The negative ones are trickier
    ndff[length(ndff):(length(ndff) - n + 1)] = dff[length(x):(length(x) - n + 1)]
  }
  #The inverses
  indff = fft(ndff/73, inverse = TRUE)
  idff = fft(dff/73, inverse = TRUE)
  if(plot){
    if(!add){
      plot(x = t, y = x, pch = 16L, xlab = "Time", ylab = "Measurement",
           main = ifelse(is.null(main), paste(n, "harmonics"), main))
      lines(y = Mod(idff), x = t, col = adjustcolor(1L, alpha = 0.5))
    }
    lines(y = Mod(indff), x = nt, ...)
  }
  ret = data.frame(time = nt, y = Mod(indff))
  return(ret)
}
