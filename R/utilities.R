# useful little functions
library(sp)
library(rgdal)
library(NISTunits)
library(htmlwidgets)
#library(FRK)

#################
# Spatial Tools #
#################
#---------------#
############
# Boxmaker #
############
# Calculate points directly north, south, east and 
# west a certain distance from given coordinates
boxmaker <- function(lon_orig=150.2269, lat_orig=-36.25201, km=150){
  # convert decimal degrees to radians
  lon_orig <- NISTdegTOradian(lon_orig)
  lat_orig <- NISTdegTOradian(lat_orig)
  # reverse harversine formula
  c = km / 6371.
  a = sin(c/2.)**2.
  dlat = 2. * asin(sqrt(a))
  dlon = 2. * asin(sqrt(a/(cos(lat_orig)**2.)))
  # convert back to decimal degrees
  lon_orig <- NISTradianTOdeg(lon_orig)
  lat_orig <- NISTradianTOdeg(lat_orig)
  dlon <- NISTradianTOdeg(dlon)
  dlat <- NISTradianTOdeg(dlat)
  # find coordinates
  north = lat_orig + dlat
  south = lat_orig - dlat
  east = lon_orig + dlon
  west = lon_orig - dlon
  # correct over the 0-360 degree line
  if (west > 360){
    west = west - 360
  }
  if (east > 360){
    east = east - 360
  }
  # round to 6 decimal places
  region = c(west, east, south, north)
  region = sapply(region, function(x) round(x, 6))
  
  # export region
  return(region)
}

################
# Circle maker #
################
circlemaker <- function(center, radius, nPoints=100){
  # centers: the data frame of centers with ID
  # radius: radius measured in kilometers
  # length per longitude changes with lattitude, so need correction
  # returns a polygon
  radiusLon <- radius / 111 / cos(center[2]/57.3) 
  radiusLat <- radius / 111
  angle <- seq(0,2*pi, length.out=nPoints)
  lon <- center[1] + radiusLon * cos(angle)
  lat <- center[2] + radiusLat * sin(angle)
  poly <- Polygon(cbind(lon, lat))
  return(poly)
}

#############
# Bug Fixes #
#############
# library(htmlwidgets)
saveWidgetFix <- function (widget,file,...){
  ## A wrapper to saveWidget which compensates for arguable BUG in
  ## saveWidget which requires `file` to be in current working
  ## directory.
  wd<-getwd()
  on.exit(setwd(wd))
  outDir<-dirname(file)
  file<-basename(file)
  setwd(outDir);
  saveWidget(widget,file=file,...)
}

# library(RGeoStats)
# Fixes issue of misnames variable in function, reported
# http://rgeostats.free.fr/forum/viewtopic.php?f=6&t=485
# don't add .shp to filename
polygon.read.format <- function (file, mode = "shape") 
{
  mode = toupper(mode)
  if (mode == "SHAPE") {
    if (requireNamespace("shapefiles", quietly = TRUE)) {
      shapefile = read.shapefile(file)
      number = length(shapefile$shp$shp)
      poly = NA
      for (ipol in 1:number) {
        polyset = shapefile$shp$shp[[ipol]]
        poly = polygon.create(polyset$points[, 1], polyset$points[, 
                                                                  2], polygon = poly)
      }
    }
    else {
      cat("You must download the package 'shapefiles' first")
    }
  }
  else {
    cat("The only format available are:\n")
    cat("SHAPE   : ShapeFile Polygont format\n")
    stop0("Unknown format : ", mode)
  }
  poly
}

########################
# Projection functions #
########################
#_________________#
###################
# UTM conversions #
###################
# lat lon to UTM
lonlat2UTM <- function(x,y,zone='56H'){
  xy = data.frame(x,y)
  coordinates(xy) <- c("x", "y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")
  res <- spTransform(xy, CRS(paste0("+proj=utm +zone=",zone," +ellps=WGS84")))
  return(data.frame(res@coords))
}
# UTM to lat lon
UTM2lonlat <- function(x,y,zone='56H'){
  xy = data.frame(x,y)
  coordinates(xy) <- c("x", "y")
  proj4string(xy) <-  CRS(paste0("+proj=utm +zone=",zone," +ellps=WGS84"))
  res <- spTransform(xy,CRS("+proj=longlat +datum=WGS84"))
  return(data.frame(res@coords))
}

###########################################
# wrappers for UTM data frame conversions #
###########################################
# UTM coordinates are rounded to avoid floating point errors
# lat lon to UTM
lonlat2UTM.df <- function(df, zone='56H'){
  out_UTM <- lonlat2UTM(df$lon, df$lat, zone=zone)
  df$lon <- round(out_UTM$x,1)
  df$lat <- round(out_UTM$y,1)
  return(df)
}
# UTM to lat lon
UTM2lonlat.df <- function(df, zone='56H'){
  out_decdeg <- UTM2lonlat(df$lon, df$lat, zone=zone)
  df$lon <- out_decdeg$x
  df$lat <- out_decdeg$y
  return(df)
}

#################################################
# wrappers for UTM rgeostats polygon conversion #
#################################################
# lon lat to UTM
lonlat2UTM.poly <- function(poly, zone='56H'){
  pdist = lonlat2UTM(poly@sets[[1]]@x, poly@sets[[1]]@y)
  poly@sets[[1]]@x <- pdist$x
  poly@sets[[1]]@y <- pdist$y
  return(poly)
}
# UTM to lon lat
UTM2lonlat.poly <- function(poly, zone='56H'){
  pdist = UTM2lonlat(poly@sets[[1]]@x, poly@sets[[1]]@y)
  poly@sets[[1]]@x <- pdist$x
  poly@sets[[1]]@y <- pdist$y
  return(poly)
}

##############################################
# wrappers for UTM penguin track conversions #
##############################################
# UTM coordinates are rounded to avoid floating point errors
# lat lon to UTM
lonlat2UTM.track <- function(df, zone='56H'){
  # make sure we do not add NA values (conversion will fail)
  na.idx <- which(is.na(df$lon))
  df_nafiltered <- df[!is.na(df$lon),]
  # do conversions
  out_UTM <- lonlat2UTM(df_nafiltered$lon, df_nafiltered$lat, zone=zone)
  x <- out_UTM$x
  y <- out_UTM$y
  # sub in NAs
  for (ii in na.idx){
    x <- c(x[1:(ii-1)],NA,x[ii:length(x)])
    y <- c(y[1:(ii-1)],NA,y[ii:length(y)])
  }
  # put in data frame
  df$lon <- round(x,1)
  df$lat <- round(y,1)
  return(df)
}
# UTM to lat lon
UTM2lonlat.track <- function(df, zone='56H'){
  #out_decdeg <- UTM2lonlat(df$lon, df$lat, zone=zone)
  #df$lon <- out_decdeg$x
  #df$lat <- out_decdeg$y
  df - 10 # throw error
  return(df)
}

####################
# Kres UTM km to m #
####################
# convert 3D kriging model output in UTM km to meters
kres_UTMkm2UTMm <- function(kres){
  kres@items$x1 <- kres@items$x1*1000
  kres@items$x2 <- kres@items$x2*1000
  kres@x0[1:2] <- kres@x0[1:2]*1000
  kres@dx[1:2] <- kres@dx[1:2]*1000
  return(kres)
}

#############
# DEBUGGING #
#############
#________________________________________________________#
##########################################################
# Search for duplicate coordinates in RGeostats database #
##########################################################
db_find_dup_coords <- function(db, colnames=c('lat','lon','rdepth')){
  coords <- as.data.frame(db@items[2:4])
  dups <- which(duplicated(coords))
  return(dups)
}

####################################################################
# Apply function that does not support NAs to a dataframe with NAs #
####################################################################
# will later make a version for vectors when needed
na.sub.FUN.df <- function(df, NA.cols, FUN, ...){
  # get na.idx
  na.idx <- which(!complete.cases(df[,NA.cols]))
  dat.idx <- which(complete.cases(df[,NA.cols]))
  # filter NAs
  df.filt <- df[dat.idx,]
  # apply the function
  df.filt <- FUN(df.filt, ...)
  # sub back into first df
  df[dat.idx,] <- df.filt
  return(df)
}

#####################
# Statistical Tools #
#####################
#------------------------#
##########
# NScore #
##########
nscore <- function(x) {                        
  # by Ashton Shortridge, 2008
  # https://www.youtube.com/watch?v=jfZKqO4TBtw
  # https://github.com/GeostatsGuy/geostatsr/blob/master/kriging_demo.R
  # Takes a vector of values x and calculates their normal scores. Returns 
  # a list with the scores and an ordered table of original values and
  # scores, which is useful as a back-transform table. See backtr().
  nscore <- qqnorm(x, plot.it = FALSE)$x  # normal score 
  trn.table <- data.frame(x=sort(x),nscore=sort(nscore))
  return (list(nscore=nscore, trn.table=trn.table))
}

backtr.nscore <- function(scores, nscore, tails='none', draw=TRUE) {
  # https://msu.edu/~ashton/temp/nscore.R
  # Given a vector of normal scores and a normal score object 
  # (from nscore), the function returns a vector of back-transformed 
  # values. One major issue is how to extrapolate to the tails. Options 
  # other than none may result in dramatically incorrect tail estimates!
  # tails options:
  # 'none' : No extrapolation; more extreme score values will revert 
  # to the original min and max values. 
  # 'equal' : Calculate magnitude in std deviations of the scores about 
  # initial data mean. Extrapolation is linear to these deviations. 
  # will be based upon deviations from the mean of the original 
  # hard data - possibly quite dangerous!
  # 'separate' :  This calculates a separate sd for values 
  # above and below the mean.
  
  if(tails=='separate') { 
    mean.x <- mean(nscore$trn.table$x)
    small.x <- nscore$trn.table$x < mean.x
    large.x <- nscore$trn.table$x > mean.x
    small.sd <- sqrt(sum((nscore$trn.table$x[small.x]-mean.x)^2)/
                       (length(nscore$trn.table$x[small.x])-1))
    large.sd <- sqrt(sum((nscore$trn.table$x[large.x]-mean.x)^2)/
                       (length(nscore$trn.table$x[large.x])-1))
    min.x <- mean(nscore$trn.table$x) + (min(scores) * small.sd)
    max.x <- mean(nscore$trn.table$x) + (max(scores) * large.sd)
    # check to see if these values are LESS extreme than the
    # initial data - if so, use the initial data.
    #print(paste('lg.sd is:',large.sd,'max.x is:',max.x,'max nsc.x is:',max(nscore$trn.table$x)))
    if(min.x > min(nscore$trn.table$x)) {min.x <- min(nscore$trn.table$x)}
    if(max.x < max(nscore$trn.table$x)) {max.x <- max(nscore$trn.table$x)}
  }
  if(tails=='equal') { # assumes symmetric distribution around the mean
    mean.x <- mean(nscore$trn.table$x)
    sd.x <- sd(nscore$trn.table$x)
    min.x <- mean(nscore$trn.table$x) + (min(scores) * sd.x)
    max.x <- mean(nscore$trn.table$x) + (max(scores) * sd.x)
    # check to see if these values are LESS extreme than the
    # initial data - if so, use the initial data.
    if(min.x > min(nscore$trn.table$x)) {min.x <- min(nscore$trn.table$x)}
    if(max.x < max(nscore$trn.table$x)) {max.x <- max(nscore$trn.table$x)}
  }
  if(tails=='none') {   # No extrapolation
    min.x <- min(nscore$trn.table$x)
    max.x <- max(nscore$trn.table$x)
  }
  min.sc <- min(scores)
  max.sc <- max(scores)
  x <- c(min.x, nscore$trn.table$x, max.x)
  nsc <- c(min.sc, nscore$trn.table$nscore, max.sc)
  
  if(draw) {plot(nsc,x, main='Transform Function')}
  back.xf <- approxfun(nsc,x) # Develop the back transform function
  val <- back.xf(scores)
  
  return(val)
}

###############################################
# Return a set of polygons from a custom grid #
###############################################
grid.polygons <- function(grid.df, ids=NULL){
  # iterate through h and v lines
  polys <- list()
  i <- 1
  for (h in min(grid.df$hline):(max(grid.df$hline)-1)){
    for (v in min(grid.df$vline):(max(grid.df$vline)-1)){
      polys[[i]] = list(c(h, v), c(h, v+1), c(h+1, v+1), c(h+1, v))
      i <- i + 1
    }
  }
  # convert to coordinates
  polys <- lapply(polys, function(poly) lapply(poly, function(p) 
    c(grid.df$lon[grid.df$hline == p[1] & grid.df$vline ==  p[2]],
      grid.df$lat[grid.df$hline == p[1] & grid.df$vline ==  p[2]])))
  # convert sublevel to dataframes
  polys <- lapply(polys, function(poly) 
    data.frame(lon=sapply(poly, function(p) p[1]),
               lat=sapply(poly, function(p) p[2])))
  # Give ids to polygons
  for (i in 1:length(polys)){
    if (!is.null(ids)){
      polys[[i]]$id <- ids[i]
    } else {
      polys[[i]]$id <- i
    }
  }
  
  # Convert to polygons (using great FRK wrapper function)
  polys <- do.call(rbind, polys)
  polys <- df_to_SpatialPolygons(polys,"id",c("lon","lat"), CRS())

  return(polys)
}

# Make normal
make_normal <- function(x, low=0, high=1){
  normalised = (x-min(x))/(max(x)-min(x))
  normalised = normalised*(high - low) + low
  return(normalised)
}

