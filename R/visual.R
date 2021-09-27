# visualisation functions
# list.of.packages <- c("ggplot2","ggmap", "viridis", 'base64enc',
#                       'plotly', 'cmocean','gganimate', 'lubridate')
# 
# # -> assuming you are using viridis with ggplot, this is now redundant as you can use 
# # scale_color_viridis_c for continuous or scale_color_viridis_d for discrete, obviously 
# # replace color with fill if raster or grid is used
# source("scripts/check_dep.R")
# check_dep(list.of.packages)

library(ggplot2)
library(ggmap)
library(ggnewscale)
library(viridis)
library(base64enc)
library(plotly)
library(cmocean)
library(gganimate)
library(lubridate)
library(gifski)
library(scales)
library(ggpubr)
library(grid)
library(patchwork)
library(lubridate)
source('R/utilities.r')

#####################
# Load/make basemap #
#####################
load.basemap <- function(mapset, zoom, suffix='', center=NA, force_generate=FALSE){
  # can pass a dataframe to center and will get middle point from lat/lon
  # make file path
  if (suffix != ''){
    suffix <- paste0('_',suffix)
  }
  mapfile <- paste0('assets/maps/',mapset,'/',mapset,'_zoom',zoom,suffix,'.rds')
  if(file.exists(mapfile) & !force_generate){
    message('Mapfile found..')
    basemap <- readRDS(mapfile)
  } else {
    # load key
    if (force_generate){
      message('Force generating map, loading from google maps API..')
    } else {
      message('No mapfile found, loading from google maps API..')
    }
    source("./keys/keys.R")
    # establish center
    if (!is.data.frame(center)){
      if (is.na(center)){
        # later I should make this be pulled from a metadata file in the mapset folder
        # probably only need to do this if I'm going to make more maps
        center <- c(150.2269, -36.25201)
      }
    }
    if (is.data.frame(center)){
      center = c(min(center$lon,na.rm=T)+(max(center$lon,na.rm=T)-min(center$lon,na.rm=T))/2,
                 min(center$lat,na.rm=T)+(max(center$lat,na.rm=T)-min(center$lat,na.rm=T))/2)
    }
    # make map
    register_google(key = rawToChar(base64decode(gmap_key)))
    basemap <- ggmap(get_googlemap(center = center,
                                   zoom = zoom, scale = 2,
                                   maptype ='satellite',
                                   color = 'color'))
    saveRDS(basemap,file=mapfile)
  }
  return(basemap)
}

#######################
# Line plot for track #
#######################
track.lineplot <- function(df, title=TRUE){
  trackplot <- ggplot(df, aes(x=lon, y=lat)) +
    geom_path() 
  if (title){
    trackplot <- trackplot + geom_title(df$id[1])
  }
  return(trackplot)
}

##################################
# Quick line plots of track list #
##################################
# input can be a list or dataframe
tracks.wrap.lineplots <- function(tracks, title='none'){
  # if need change list to dataframe
  if (is.vector(tracks)){
    tracks <- do.call('rbind', tracks)
  }
  # make factor
  tracks$id <- as.factor(tracks$id)
  # make the plot
  trackplot <- ggplot(tracks, aes(x=lon, y=lat, color=id)) +
    geom_path() + facet_wrap(tracks$id, scales="free") + 
    theme(legend.position="none")
  # if title add
  if (title != 'none'){
    make_title <- TRUE
  } else {
    make_title <- FALSE
  }
  if (make_title){
    trackplot <- trackplot + ggtitle(title)
  }
  return(trackplot)
}

###################################################################################
# Quick line plots of track list with diagonal distance shown in facet wrap title #
###################################################################################
# input can be a list or dataframe
tracks.wrap.lineplots.dist <- function(tracks, title='none', dist_ls){
  # if need change list to dataframe
  if (is.vector(tracks)){
    tracks <- do.call('rbind', tracks)
  }
  # make factor
  tracks$id <- as.factor(tracks$id)
  # load in dist data to id
  tracks$idDist <- as.integer(tracks$id)
  tracks$idDist <- as.factor(paste0(tracks$id, ' | ', round(unlist(dist_ls[tracks$idDist]), 2), ' km'))
  # make the plot
  trackplot <- ggplot(tracks, aes(x=lon, y=lat, color=id)) +
    geom_path() + facet_wrap(tracks$idDist, scales="free") + 
    theme(legend.position="none")
  # if title add
  if (title != 'none'){
    make_title <- TRUE
  } else {
    make_title <- FALSE
  }
  if (make_title){
    trackplot <- trackplot + ggtitle(title)
  }
  return(trackplot)
}

#########################################
# Plot specified tracks over a map file #
#########################################
# tracks_ids can be a vector or single string. If empty returns all in dataframe
# colfactor is an optional string for colouring
tracks.map <- function(df, track_ids=c(), plot_points=TRUE, plot_line=TRUE,
                       mapset="montague", suffix='',
                       zoom=9, legend=TRUE, title=NA, colfactor=NA, colfactor_alpha=FALSE, 
                       padlat=0.01, padlon=0.07, facet_survey=F, facet_year=F, 
                       custom.limits=NULL){
  # if need change list to dataframe
  if (is.vector(df)){
    df <- do.call('rbind', df)
  }
  # filter tracks
  if (length(track_ids) > 0){
    df <- df[df$id %in% track_ids,]
  }
  #  Add year if needed
  if (facet_year){
    if (is.null(df$year)){
      df$year <- year(df$dtUTC)
    }
  }
  # load or make map
  basemap <- load.basemap(mapset, zoom, suffix)
  # make map
  trackmap <- basemap +
    scale_y_continuous(limits = c(min(df$lat, na.rm=TRUE)-padlat, max(df$lat, na.rm=TRUE)+padlat), expand = c(0, 0)) +
    scale_x_continuous(limits = c(min(df$lon, na.rm=TRUE)-padlon, max(df$lon, na.rm=TRUE)+padlon), expand = c(0, 0))
  # add lines
  if (plot_line){
    if ((is.na(colfactor))){
      trackmap <- trackmap + geom_path(data=df, aes(x=lon,y=lat,color=id), size=0.5)
    } else {
      trackmap <- trackmap + geom_path(data=df, aes(x=lon,y=lat), size=0.5, color='grey')
    }
  }
  # add points
  if (plot_points){
    if ((is.na(colfactor))){
      trackmap <- trackmap + geom_point(data=df, aes(x=lon,y=lat,color=id), size=0.7)
    } else {
      if (!colfactor_alpha){
        trackmap <- trackmap + geom_point(data=df, aes(x=lon,y=lat,color=df[,colfactor]), size=0.7) +
          labs(color=colfactor)
      } else {
        alpha_val <- ((df[,colfactor]-min(df[,colfactor])) / (max(df[,colfactor])-min(df[,colfactor])))/4
        
        trackmap <- trackmap + geom_point(data=df, aes(x=lon,y=lat,color=df[,colfactor],alpha=alpha_val),
                                          size=0.7) +
          labs(color=colfactor) + scale_alpha(guide = 'none')
      }
    }
  }
  # use better colour scale if colfactor is numeric
  if ((!is.na(colfactor))){
    if (is.numeric(df[,colfactor])){
      trackmap <- trackmap + scale_colour_gradientn(colours = cmocean('thermal')(10))
    }
  }
  # Custom lonlat limits
  if (!is.null(custom.limits)){
    trackmap <- trackmap +
      scale_x_continuous(limits = c(custom.limits$xmin-padlon, custom.limits$xmax+padlon), expand = c(0, 0)) +
      scale_y_continuous(limits = c(custom.limits$ymin-padlat, custom.limits$ymax+padlat), expand = c(0, 0))
  }
  # legend remove
  if (!legend){
    trackmap <- trackmap + theme(legend.position = "none")
  }
  # title
  if (!is.na(title)){
    trackmap <- trackmap + ggtitle(title)
  }
  # facet survey
  if (facet_survey){
    trackmap <- trackmap + facet_wrap(~survey_id)
  }
  # facet year
  if (facet_year){
    trackmap <- trackmap + facet_wrap(~year)
  }
  return(trackmap)
}

##################################
# Quick line plots of track list #
##################################
# input can be a list or dataframe
map.points <- function(points, size=0.7, mapfile="./assets/maps/map_montague_zoom8.rds", zoom=8){
  # load or make map
  if(file.exists(mapfile)){
    print('Mapfile found..')
    basemap <- readRDS(mapfile)
  } else {
    # load key
    print('No mapfile found, loading from google maps API..')
    source("./keys/keys.R")
    # make map
    center = c(min(df$lon)+(max(points$lon)-min(points$lon))/2,min(points$lat)+(max(points$lat)-min(points$lat))/2)
    register_google(key = rawToChar(base64decode(gmap_key)))
    basemap <- ggmap(get_googlemap(center = center,
                                   zoom = zoom, scale = 2,
                                   maptype ='satellite',
                                   color = 'color'))
    saveRDS(basemap,file=mapfile)
  }
  # make map
  pointmap <- basemap +
    scale_y_continuous(limits = c(min(points$lat, na.rm=TRUE)-0.01, max(points$lat, na.rm=TRUE)+0.01), expand = c(0, 0)) +
    scale_x_continuous(limits = c(min(points$lon, na.rm=TRUE)-0.07, max(points$lon, na.rm=TRUE)+0.07), expand = c(0, 0)) +
    geom_point(data=points, aes(x=lon,y=lat, color='orange'), size=size) +
    theme(legend.title = element_blank())
  # return map
  return(pointmap)
}

#############################
# Ocean plot cheack wrapper #
#############################
ocean.land.plot.quick <- function(df, shp="assets/shp/NSW-coastline/edited/polygon/NSW-coastline_WGS84"){
  coast <- readOGR(dsn=path.expand(str_sub(shp,1,-(nchar(basename(shp)))-1)), 
                   layer=basename(shp)) # load shapefile
  df <- track.ocean(df)
  # make plot
  ocean.land.plot(df, coast, title='quick')
}

################################
# plot points as ocean or land #
################################
ocean.land.plot <- function(df, coast, title){
  # point colours
  color <- c('#E69F00','#56B4E9')
  labels <- c('land','ocean')
  # create bounding box 20km with lighthouse as center
  xmin <- 150.003861
  xmax <- 150.449939
  ymin <- -36.431874
  ymax <- -36.072146
  # make plot - requires 'mapproj'
  plot1 <- suppressMessages(
    ocean.land.plot.subfunc(xmin,xmax,ymin,ymax,coast,df,color,paste(title,'- Zoom1'),labels))
  # create bounding box 2km with lighthouse as center
  xmin <- 150.204596
  xmax <- 150.249204
  ymin <- -36.269996
  ymax <- -36.234024
  # make plot - requires 'mapproj'
  plot2 <- suppressMessages(
    ocean.land.plot.subfunc(xmin,xmax,ymin,ymax,coast,df,color,paste(title,'- Zoom2'),labels))
  # create bounding box 500m with lighthouse as center
  xmin <- 150.221324
  xmax <- 150.232476
  ymin <- -36.256507
  ymax <- -36.247513
  # make plot - requires 'mapproj'
  plot3 <- suppressMessages(
    ocean.land.plot.subfunc(xmin,xmax,ymin,ymax,coast,df,color,paste(title,'- Zoom3'),labels))
  # make plots
  plot_all <- grid.arrange(arrangeGrob(plot1,plot2,plot3,
                                       layout_matrix = cbind(c(1,NA),c(1,3),c(2,3),c(2,NA))))
}

######################
# ocean map function #
######################
ocean.land.plot.subfunc <- function(xmin, xmax, ymin, ymax, coast, df, color, title, labels){
  ggplot() + 
    geom_polygon(data=coast, aes(x=long, y=lat, group=group)) + 
    coord_fixed(1) +
    coord_map(xlim=c(xmin, xmax), ylim=c(ymin, ymax)) +
    geom_point(data=df, aes(x=lon, y=lat, group=ocean, color=ocean), 
               size=2, shape=18, na.rm=TRUE) +
    scale_color_manual(values=color, labels=labels, drop=F) +
    labs(x="Longitude",y="Latitude") + 
    theme(legend.title = element_blank()) +
    ggtitle(title)
}

#################
# check polygon #
#################
plot_polygon <- function(poly){
  # load or make map
  if(file.exists("./assets/maps/map_montague_zoom9.rds")){
    print('Mapfile found..')
    mapbase <- readRDS("./assets/maps/map_montague_zoom9.rds")
  } else {
    # load key
    print('No mapfile found, loading from google maps API..')
    source("./keys/keys.R")
    # make map
    register_google(key= rawToChar(base64decode(gmap_key)))
    mapbase <- ggmap(get_googlemap(center = c(lon = 150.2269, lat = -36.29201),
                                   zoom = 10, scale = 2,
                                   maptype ='satellite',
                                   color = 'color'))
    saveRDS(hmm_map,file="./assets/maps/map_montague_zoom9.rds")
  }
  # make map
  print('Plotting map..')
  df_poly <- data.frame(lon=poly@sets[[1]]@x, lat=poly@sets[[1]]@y)
  # make dataframe for geom_point
  gmap <- mapbase +
    geom_polygon(data = df_poly, aes(x = lon, y = lat), colour="white",fill=NA)
  return(gmap)
}

############
# Plot mpm #
############
plot_mpm <- function(df){
  # load or make map
  if(file.exists("./assets/maps/map_montague_zoom9.rds")){
    print('Mapfile found..')
    mpm_map <- readRDS("./assets/maps/map_montague_zoom9.rds")
  } else {
    # load key
    print('No mapfile found, loading from google maps API..')
    source("./keys/keys.R")
    # make map
    register_google(key= rawToChar(base64decode(gmap_key)))
    hmm_map <- ggmap(get_googlemap(center = c(lon = 150.2269, lat = -36.29201),
                                   zoom = 10, scale = 2,
                                   maptype ='satellite',
                                   color = 'color'))
    saveRDS(hmm_map,file="./assets/maps/map_montague_zoom9.rds")
  }
  # make map
  print('Plotting map..')
  # make dataframe for geom_point
  # point map
  gmap <- mpm_map + 
    geom_point(data=df, aes(x=lon,y=lat,color=mpm,alpha=1-mpm), size=0.6) +
    scale_colour_gradient2(high='blue', mid='yellow', low='white') +
    scale_alpha(guide = 'none') + 
    scale_y_continuous(limits = c(min(df$lat, na.rm=TRUE)-0.01, max(df$lat, na.rm=TRUE)+0.01), expand = c(0, 0)) +
    scale_x_continuous(limits = c(min(df$lon, na.rm=TRUE)-0.07, max(df$lon, na.rm=TRUE)+0.07), expand = c(0, 0))
  
  return(gmap)
}

########################
# check polygon points #
########################
check_polygon <- function(df, poly, mask){
  # load or make map
  if(file.exists("./assets/maps/map_montague_zoom9.rds")){
    print('Mapfile found..')
    mapbase <- readRDS("./assets/maps/map_montague_zoom9.rds")
  } else {
    # load key
    print('No mapfile found, loading from google maps API..')
    source("./keys/keys.R")
    # make map
    register_google(key= rawToChar(base64decode(gmap_key)))
    mapbase <- ggmap(get_googlemap(center = c(lon = 150.2269, lat = -36.29201),
                                   zoom = 10, scale = 2,
                                   maptype ='satellite',
                                   color = 'color'))
    saveRDS(hmm_map,file="./assets/maps/map_montague_zoom9.rds")
  }
  # make map
  print('Plotting map..')
  df$mask <- mask
  df_poly <- data.frame(lon=poly@sets[[1]]@x, lat=poly@sets[[1]]@y)
  # make dataframe for geom_point
  gmap <- mapbase + geom_point(data=df, aes(x=lon,y=lat,color=mask), size=0.3, alpha=0.7) +
    scale_y_continuous(limits = c(min(df$lat, na.rm=TRUE)-0.01, max(df$lat, na.rm=TRUE)+0.01), expand = c(0, 0)) +
    scale_x_continuous(limits = c(min(df$lon, na.rm=TRUE)-0.07, max(df$lon, na.rm=TRUE)+0.07), expand = c(0, 0)) +
    geom_polygon(data = df_poly, aes(x = lon, y = lat), colour="white",fill=NA)
  return(gmap)
}

###########################
# check dive/lat/lon data #
###########################
check.dive.latlon <- function(df, range){
  # keep only 25 tracks
  levels <- levels(df$id)[range]
  df <- df[df$id == levels,]
  # red if NA
  df$col <- 'depth'
  df[is.na(df$depth), 'col'] <- 'NA'
  ggplot(df, aes(x=lon, y=lat, color=col)) +
    geom_point() +
    facet_wrap(df$id)
}

####################################
# Pilot time plot (separate hours) #
####################################
plot_time <- function(df, sep_hour=FALSE){
  # load or make map
  if(file.exists("./assets/maps/hmm_mapfile.rds")){
    print('Mapfile found..')
    mpm_map <- readRDS("./assets/maps/hmm_mapfile.rds")
  } else {
    # load key
    print('No mapfile found, loading from google maps API..')
    source("./keys/keys.R")
    # make map
    register_google(key= rawToChar(base64decode(gmap_key)))
    hmm_map <- ggmap(get_googlemap(center = c(lon = 150.2269, lat = -36.29201),
                                   zoom = 10, scale = 2,
                                   maptype ='satellite',
                                   color = 'color'))
    saveRDS(hmm_map,file="./assets/maps/hmm_mapfile.rds")
  }
  # make map
  print('Plotting map..')
  # make dataframe for geom_point
  # point map
  gmap <- mpm_map + 
    scale_y_continuous(limits = c(min(df$lat, na.rm=TRUE)-0.01, max(df$lat, na.rm=TRUE)+0.01), expand = c(0, 0)) +
    scale_x_continuous(limits = c(min(df$lon, na.rm=TRUE)-0.07, max(df$lon, na.rm=TRUE)+0.07), expand = c(0, 0))
  if (sep_hour){
    gmap <- gmap + geom_point(data=df, aes(x=lon,y=lat,color='#ff6600'), size=0.8) +
      facet_wrap(~hour)
  } else {
    gmap <- gmap + geom_point(data=df, aes(x=lon,y=lat,color=hour), size=0.8) +
      scale_color_viridis()
  }
  return(gmap)
}

#############################
# plot speed graph of track #
#############################
track.speed.plot <- function(track, filter_speed=NA){
  track <- track.distance(track)
  track <- track.time(track)
  track <- track.speed(track)
  if (!is.na(filter)){
    track <- track[track$speed < filter_speed,]
    track <- track.distance(track)
    track <- track.time(track)
    track <- track.speed(track)
  }
  # make plot
  ggplot(data=track, aes(x=dtUTC, y=speed)) + 
    geom_line(na.rm=TRUE, color='grey') +
    geom_point(na.rm=TRUE, size=0.4) +
    labs(x="Datetime",y="ms-1") + 
    ggtitle(track$id[2])
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

##########################
# Animate penguin tracks #
##########################
# col_factor options: 'id', 'year', 'year_facet'
# style options: 'long_wake'
anim.peng.tracks <- function(df, track_ids=c(), nframes=150, col_factor='id', 
                             styles=c('default'), output_fn="./output/animations/peng_tracks.gif", 
                             save_gif=FALSE, custom_dim=NA, mapset='montague', zoom=10,
                             map_force_generate=FALSE, lat_range_custom=NA){
  # cutom_dim as vector c(height, width) in pixels
  # make relative local datetimes
  df <- track.dtAEST(df)
  # year as factor
  df$year <- as.factor(year(df$dtAEST))
  # filter tracks
  if (length(track_ids) > 0){
    df <- df[df$id %in% track_ids,]
  }
  # read in map data
  mapfile <- load.basemap(mapset=mapset, zoom=zoom, force_generate=map_force_generate)
  # work out if day one or two of the trip
  i <- 1
  df_ls <- split(df, df$id)
  for (df in df_ls){
    min_dt <- df$dtAEST[1]
    min_dt <- ISOdatetime(year(min_dt), month(min_dt), day(min_dt), 0, 0, 0)
    df$tripDay <- as.numeric(difftime(df$dtAEST,  min_dt, units='days') > 1) + 1
    df_ls[[i]] <- df
    i <- i + 1
  }
  df <- do.call('rbind', df_ls)
  row.names(df) <- NULL
  # make relative datetimes
  df$hour <- hour(df$dtAEST)
  df$min <- minute(df$dtAEST)
  df$datetime_relative <- ISOdatetime(2000, 01, df$tripDay, as.integer(df$hour), as.integer(df$min), 0)
  
  # make the map
  if (is.na(lat_range_custom)){
    peng.map <- mapfile + 
      scale_y_continuous(limits = c(min(df$lat, na.rm=TRUE)-0.01, max(df$lat, na.rm=TRUE)+0.01), expand = c(0, 0)) +
      scale_x_continuous(limits = c(min(df$lon, na.rm=TRUE)-0.07, max(df$lon, na.rm=TRUE)+0.07), expand = c(0, 0)) +
      labs(title = 'Relative Time: {strftime(frame_time, "%H:%M")}')
  } else {
    peng.map <- mapfile + 
      scale_y_continuous(limits = c(lat_range_custom[1], lat_range_custom[2]), expand = c(0, 0)) +
      scale_x_continuous(limits = c(min(df$lon, na.rm=TRUE)-0.07, max(df$lon, na.rm=TRUE)+0.07), expand = c(0, 0)) +
      labs(title = 'Relative Time: {strftime(frame_time, "%H:%M")}')
  }
  
  ## styles ##
  # id
  if (col_factor == 'id'){
    peng.map <- peng.map + geom_point(data=df, aes(x=lon,y=lat,color=id), size=0.7)
    # remove legend if more than 10 ids
    if (length(levels(droplevels(df$id))) > 10){
      peng.map <- peng.map + theme(legend.position = "none")
    }
  }
  # year
  if (col_factor == 'year'){
    peng.map <- peng.map + geom_point(data=df, aes(x=lon,y=lat,color=year), size=0.7)
  }
  # year_facet
  if (col_factor == 'year_facet'){
    peng.map <- peng.map + geom_point(data=df, aes(x=lon,y=lat,color=year), size=0.7) +
      facet_wrap(~year) + theme(legend.position = "none")
  }
  
  # make animation
  if ('long_wake' %in% styles){
    peng.anim <- peng.map +
      transition_time(df$datetime_relative) +
      shadow_wake(wake_length = 0.7, alpha=TRUE)
  } else {
    peng.anim <- peng.map +
      transition_time(df$datetime_relative) +
      shadow_wake(wake_length = 0.2, alpha=TRUE)
  }
  
  # output
  if (save_gif){
    if (is.na(custom_dim)){
      anim_save(output_fn, peng.anim, nframes=nframes)
    } else {
      anim_save(output_fn, peng.anim, nframes=nframes, height=custom_dim[1], width=custom_dim[2])
    }
  } else {
    # show
    return(peng.anim)
  }
}

###################
# Quick dive plot #
###################
dive.plot <- function(tracks, ids='all', covar=NA, covar_lab=NA, depth_col='depth',
                      remove_surface=FALSE, surface_thresh=0.2, downsample=NA,
                      facet_scale=c('free','fixed')){
  # downsample sets each lot to have a max number of rows
  # set argument matches
  facet_scale <- match.arg(facet_scale)
  # set all ids
  if (ids == 'all'){
    ids <- as.character(unique(tracks$id))
  }
  # fitler
  df <- tracks[tracks$id %in% ids,]
  df$id <- factor(df$id)
  # if reduce row numbers (for large files)
  if (!is.na(downsample)){
    df_ls <- split(df, df$id)
    for (i in 1:length(df_ls)){
      N <- as.integer(nrow(df_ls[[i]])/downsample)
      df_ls[[i]] <- df_ls[[i]][(seq(1,to=nrow(df_ls[[i]]),by=N)),]
    }
    df <- do.call('rbind', df_ls)
  }
  # check depth column exists
  if (!depth_col %in% colnames(df)){
    # get columns with depth
    depth_options <- colnames(df)[grep('depth', colnames(df))]
    warning('"',depth_col,'" not column not found. Using "',depth_options[length(depth_options)],'"')
    depth_col <- depth_options[length(depth_options)]
  }
  # change covar label if not set
  if (is.na(covar_lab)){
    covar_lab <- covar
  }
  # remove surface
  if (remove_surface){
    df <- df[df[,depth_col] > surface_thresh,]
  }
  if (!is.na(covar)){
    dive_plot <- ggplot(data=df, aes(x=dtUTC, y=df[,depth_col], color=df[,covar])) +
      geom_path(size=0.5) + geom_point(size=1.5) + 
      scale_y_continuous(trans="reverse") +
      scale_colour_gradientn(colours = cmocean('thermal')(10)) +
      xlab('Timestamp') +
      ylab('Depth (m)') +
      labs(color=covar_lab)
  } else {
    dive_plot <- ggplot(data=df, aes(x=dtUTC, y=df[,depth_col])) +
      geom_path(size=0.5) + geom_point(size=1.5) + 
      scale_y_continuous(trans="reverse") +
      xlab('Timestamp') +
      ylab('Depth (m)')
  }
  dive_plot <- dive_plot + facet_wrap(~id, scales=facet_scale)
  return(dive_plot)
}

#############################
# compare dive stat metrics #
#############################
dive.plot.comp.stats <- function(tracks, id, vars=c('depth','temp'), size=0.5){
  # filter
  df <- tracks[tracks$id == id,]
  # find percentile column
  percen_num <- as.numeric(gsub("[^0-9.]","",colnames(df)[grep('percen', colnames(df))][1]))
  # make depth plot
  dive_plot <- ggplot(data=df) +
    geom_path(mapping=aes(x=dtUTC, y=depth_mean, color="Mean"),size=size)  +
    geom_path(mapping=aes(x=dtUTC, y=depth_min, color="Min"),size=size)  +
    geom_path(mapping=aes(x=dtUTC, y=depth_max, color="Max"),size=size)  +
    geom_path(mapping=aes(x=dtUTC, y=df[,paste0('depth_percen',percen_num)],
                          colour=paste0(ordinal(percen_num),' Percentile')),size=size)  +
    scale_y_continuous(trans="reverse") +
    xlab('Timestamp') +
    ylab('Depth (m)') +
    ggtitle(paste0('Depth Metrics - ',df$id[1]))
  # make temp plot
  temp_plot <- ggplot(data=df) +
    geom_path(mapping=aes(x=dtUTC, y=temp_mean, color="Mean"),size=size)  +
    geom_path(mapping=aes(x=dtUTC, y=temp_min, color="Min"),size=size)  +
    geom_path(mapping=aes(x=dtUTC, y=temp_max, color="Max"),size=size)  +
    geom_path(mapping=aes(x=dtUTC, y=df[,paste0('temp_percen',percen_num)],
                          colour=paste0(ordinal(percen_num),' Percentile')),size=size)  +
    xlab('Timestamp') +
    ylab('Temperature (°C)') +
    ggtitle(paste0('Temperature Metrics - ',df$id[1]))
  # arrange
  p <- ggarrange(dive_plot, temp_plot, nrow=2)
  return(p)
}

#############################################
# CEFAS - Interactive dive plot with Plotly #
#############################################
# makes a high resolution interatcive plot for a CEFAS acceleromtry 
# used to investigate G126f06 deep dives
dive.plot.CEFAS.plotly <- function(ids,  snapshot=10000, downsample=NA, 
                                   zero_flying_penguins=TRUE, AEST=TRUE,
                                   output_dir='./output/plotly_html_dump/', 
                                   accel_dir='./data/dive/CEFAS/'){
  if (!is.na(snapshot)){
    warning('Snapshot set to: ',prettyNum(snapshot, big.mark=","),' points. Dive logs may be trimmed.')
  }
  for (id in ids){
    # load file
    accel <- readRDS(paste0(accel_dir,id,'.rds'))
    accel$id <- droplevels(accel$id)
    # zero negative depths
    if (zero_flying_penguins){
      accel <- zero.flying.penguins(accel)
    }
    # make snapshot
    if (!is.na(snapshot) & nrow(accel) > snapshot*2.1){
      accel <- accel[(floor(nrow(accel)/2)):(floor(nrow(accel)/2)+snapshot-1),]
    }
    # downsample dataframe
    if (!is.na(downsample) & downsample < nrow(accel)){
      N <- as.integer(nrow(accel)/downsample)
      accel <- accel[(seq(1,to=nrow(accel),by=N)),]
      row.names(accel) <- NULL
    }
    # convert to local timezone
    if (AEST){
      accel$dtUTC <- convert2AEST(accel$dtUTC)
      timestamp_title = 'Timestamp (AEST)'
    } else {
      timestamp_title = 'Timestamp (UTC)'
    }
    # make traces
    trace_time <- accel$dtUTC
    trace_depth <- accel$depth
    
    # make dataframe
    data <- data.frame(trace_time, trace_depth)
    
    # make plot
    fig_depth <- plot_ly(data, x=~trace_time, y=~trace_depth, 
                         name='Depth',type='scatter',mode='lines+markers') %>%
      layout(yaxis = list(autorange = "reversed", title = 'Depth (m)'),
             xaxis = list(range = c(min(trace_time),max(trace_time)),title = timestamp_title),
             title = paste0('CEFAS divelog for ',accel$id[1]))
    # save results
    output_fn <- paste0(output_dir,'CEFAS_divelog_',accel$id[1],'.html')
    saveWidgetFix(fig_depth, output_fn, selfcontained=TRUE)
    message('Plot saved to: ',output_fn)
  }
}

#############################################
# AXY - Interactive dive plot with Plotly #
#############################################
# makes a high resolution interatcive plot for a AXY acceleromtry
dive.plot.AXY.plotly <- function(ids,  snapshot=10000, downsample=NA, 
                                 zero_flying_penguins=TRUE, AEST=TRUE,
                                 output_dir='./output/plotly_html_dump/',
                                 AXY_dive_fn='./data/dive/AXY/penguins_depth_temp_axy.rds'){
  # load data
  dive_df <- readRDS(AXY_dive_fn)
  # convert to local time
  if (AEST){
    dive_df$dtUTC <- convert2AEST(dive_df$dtUTC)
    timestamp_title = 'Timestamp (AEST)'
  } else {
    timestamp_title = 'Timestamp (UTC)'
  }
  if (!is.na(snapshot)){
    warning('Snapshot set to: ',prettyNum(snapshot, big.mark=","),' points. Dive logs may be trimmed.')
  }
  for (id in ids){
    # filter to track
    df <- dive_df[dive_df$id == id,]
    # zero negative depths
    if (zero_flying_penguins){
      df <- zero.flying.penguins(df)
    }
    # make snapshot
    if (!is.na(snapshot) & nrow(df) > snapshot*2.1){
      df <- df[(floor(nrow(df)/2)):(floor(nrow(df)/2)+snapshot-1),]
    }
    # downsample dataframe
    if (!is.na(downsample) & downsample < nrow(df)){
      N <- as.integer(nrow(df)/downsample)
      df <- df[(seq(1,to=nrow(df),by=N)),]
      row.names(df) <- NULL
    }
    # make traces
    trace_time <- df$dtUTC
    trace_depth <- df$depth
    
    # make dataframe
    data <- data.frame(trace_time, trace_depth)
    # make plot
    fig_depth <- plot_ly(data, x=~trace_time, y=~trace_depth, 
                         name='Depth',type='scatter',mode='lines+markers') %>%
      layout(yaxis = list(autorange = "reversed", title = 'Depth (m)'),
             xaxis = list(range = c(min(trace_time),max(trace_time)),title = timestamp_title),
             title = paste0('Axy divelog for ',df$id[1]))
    # save results
    output_fn <- paste0(output_dir,'Axy_divelog_',df$id[1],'.html')
    saveWidgetFix(fig_depth, output_fn, selfcontained=TRUE)
    message('Plot saved to: ',output_fn)
  }
}

###############################################
# compare logged and krigged temperature data #
###############################################
plot.temp.log.vs.krig <- function(tracks, ids=NA, logger_temp_name='temp_percen85'){
  # if no ids selected sample 9 at random
  if (is.na(ids)){
    tracks <- tracks[tracks$id %in% sample(unique(tracks$id), 9),]
  } else {
    tracks <- tracks[tracks$id %in% ids,]
  }
  p <- ggplot(tracks, aes(x=dtUTC)) +
    geom_line(aes(y=tracks[,logger_temp_name]), color='#0373fc') +
    geom_line(aes(y=krig_temperature), color='#404040') +
    facet_wrap(~id, scales='free')
  return(p)
}

##############################
# Plot 2-state moveHMM model #
##############################
plot.HMM.2state <- function(m.HMM, type=c('density','split_points','combined_points'), zoom=9,
                                suffix='', lon_lim=NA, lat_lim=NA, grid_year=NA){
  # set arguments
  type <- match.arg(type)
  # load or make map
  basemap <- load.basemap('montague', zoom, suffix=suffix)
  # make map
  print('Plotting map..')
  # make dataframe for geom_point
  df <- m.HMM$data
  df$state <- as.character(viterbi(m.HMM))
  df$state[df$state == "1"] <- "Foraging"
  df$state[df$state == "2"] <- "Travel"
  df$state <- as.factor(df$state)
  df$grid_year <- grid_year
  # make the map
  if (type == 'combined_points'){
    p <- basemap + geom_point(data=df, aes(x=x, y=y, color=state), size=0.3, alpha=0.7) +
      scale_y_continuous(limits = c(min(df$y, na.rm=TRUE)-0.01, max(df$y, na.rm=TRUE)+0.01), expand = c(0, 0)) +
      scale_x_continuous(limits = c(min(df$x, na.rm=TRUE)-0.07, max(df$x, na.rm=TRUE)+0.07), expand = c(0, 0))
  }
  if (type == 'split_points'){
    p <- basemap + geom_point(data=df, aes(x=x, y=y, color=state), size=0.3, alpha=0.7) +
      facet_wrap(~ state) +
      scale_y_continuous(limits = c(min(df$y, na.rm=TRUE)-0.01, max(df$y, na.rm=TRUE)+0.01), expand = c(0, 0)) +
      scale_x_continuous(limits = c(min(df$x, na.rm=TRUE)-0.07, max(df$x, na.rm=TRUE)+0.07), expand = c(0, 0))
  }
  if (type == 'density'){
    p <- basemap + stat_density_2d(data=df, aes(x=x, y=y, fill=stat(nlevel), alpha=..level..), 
                                      geom="polygon", size = 0.01, bins=12) + 
      scale_fill_gradientn(colours = cmocean('thermal')(12)) +
      scale_y_continuous(limits = c(min(df$y, na.rm=TRUE)-0.01, max(df$y, na.rm=TRUE)+0.01), expand = c(0, 0)) +
      scale_x_continuous(limits = c(min(df$x, na.rm=TRUE)-0.07, max(df$x, na.rm=TRUE)+0.07), expand = c(0, 0)) +
      scale_alpha(guide = 'none')
    if (is.na(grid_year)){
      p <- p + facet_wrap(~state)
    } else {
      p <- p + facet_grid(cols=vars(state), rows=vars(grid_year))
    }
  }
  # custom boundaries
  if (!is.na(lon_lim))
    p <- p + scale_x_continuous(limits = c(lon_lim[1], lon_lim[2]), expand = c(0, 0))
  if (!is.na(lat_lim))
    p <- p + scale_y_continuous(limits = c(lat_lim[1], lat_lim[2]), expand = c(0, 0))
  return(p)
}

#######################################
# Plot HMM ofraging density by survey #
#######################################
plot.HMM.2state.density.survey <- function(m.HMM, survey_vector=NA, zoom=9, suffix='', 
                                           lon_lim=NA, lat_lim=NA){
  # load or make map
  basemap <- load.basemap('montague', zoom, suffix=suffix)
  # make map
  print('Plotting map..')
  # make dataframe for geom_point
  df <- m.HMM$data
  df$state <- as.character(viterbi(m.HMM))
  df$state[df$state == "1"] <- "Foraging"
  df$state[df$state == "2"] <- "Travel"
  df$state <- as.factor(df$state)
  df$survey <- survey_vector
  # get only foraging
  df <- df[df$state == 'Foraging',]
  # make the map
  p <- basemap + stat_density_2d(data=df, aes(x=x, y=y, fill=stat(nlevel), alpha=..level..), 
                                 geom="polygon", size = 0.01, bins=12) + 
    scale_fill_gradientn(colours = cmocean('thermal')(12)) +
    scale_y_continuous(limits = c(min(df$y, na.rm=TRUE)-0.01, max(df$y, na.rm=TRUE)+0.01), expand = c(0, 0)) +
    scale_x_continuous(limits = c(min(df$x, na.rm=TRUE)-0.07, max(df$x, na.rm=TRUE)+0.07), expand = c(0, 0)) +
    scale_alpha(guide = 'none')
  if (!is.na(survey_vector)){
    p <- p + facet_wrap(~survey)
  }
  # custom boundaries
  if (!is.na(lon_lim))
    p <- p + scale_x_continuous(limits = c(lon_lim[1], lon_lim[2]), expand = c(0, 0))
  if (!is.na(lat_lim))
    p <- p + scale_y_continuous(limits = c(lat_lim[1], lat_lim[2]), expand = c(0, 0))
  return(p)
}

#######################
# Investigate rafting #
#######################
# plot a moveHMM model and check where all the island rafting occurs
plot.check.rafting.HMM <- function(m.HMM, radius, center=c(150.2269, -36.25201)){
  # load basemap
  basemap <- load.basemap('montague', zoom=12)
  # process model data
  df <- m.HMM$data
  df$state <- as.character(viterbi(m.HMM))
  df$state[df$state == "1"] <- "ARS"
  df$state[df$state == "2"] <- "Travel"
  df$state <- as.factor(df$state)
  # plot points
  p <- basemap + geom_point(data=df, aes(x=x, y=y, color=state), size=0.3, alpha=0.7) +
    facet_wrap(~state) +
    geom_polygon(data=fortify(circlemaker(center, radius)), aes(long, lat), color='red', fill=NA)
  return(p)
}

##########################################
# Investigate distributions of variables #
##########################################
# can submit one variable or multiple
plot.distributions <- function(tracks, vars, titles=NULL){
  if (length(vars) == 1){
    p <- ggplot(tracks, aes(x=tracks[,vars], col=factor(year(dtUTC)))) +
      geom_density() +
      labs(col='Year') +
      ggtitle(titles)  +
      xlab(vars)
  } else {
    p_ls <- list()
    for (i in 1:length(vars)){
      dat <- tracks[,c('dtUTC',vars[i])]
      colnames(dat)[2] <- 'var'
      p_ls[[i]]  <- ggplot(dat, aes(x=var, col=factor(year(dtUTC)))) +
        geom_density() +
        labs(col='Year') +
        ggtitle(titles[i]) +
        xlab(vars[i])
    }
    p <- ggarrange(plotlist=p_ls, common.legend=T)
  }
  return(p)
}

#############################################
# Plot survey acoustic energy distributions #
#############################################
plot.acoustic.distribution <- function(agg.data_all){
  p1 <- ggplot(agg.data_all, aes(x=Sv_mean, fill=survey_id)) +
    geom_density(alpha=.5, linetype='dashed', size=.2)
  p2 <- ggplot(agg.data_all, aes(x=NASC, fill=survey_id)) +
    geom_density(alpha=.5, linetype='dashed', size=.2) +
    xlim(0, quantile(agg.data_all$NASC, .75))
  p3 <- (p1 | p2)
  return(p3)
}

###############################
# Plot summary stat box plots #
###############################
sumStats.boxplot <- function(sumStats){
  p1 <- ggplot(sumStats, aes(x=year, y=distance_total)) +
    geom_boxplot() +
    ggtitle('Mean Distance Travelled (km)')
  p2 <- ggplot(sumStats, aes(x=year, y=displacement_max)) +
    geom_boxplot() +
    ggtitle('Mean Maximum Displacement (km)')
  p3 <- ggplot(sumStats, aes(x=year, y=duration)) +
    geom_boxplot() +
    ggtitle('Mean Trip Duration (hours)')
  return(ggarrange(plotlist = list(p1,p2,p3)))
}

####################################################
# Plot step and angle distributions for HMM models #
####################################################
plot.step.angle.dist <- function(tracks){
  if(!exists('prepData.wrapper', inherits=T)){
    message('Loading HMM functions...')
    source('./scripts/HMM/HMM_functions.R')
  }
  data <- prepData.wrapper(tracks)
  step <- ggplot(data, aes(x=step, col=factor(year(dtUTC)))) +
    geom_density() +
    labs(col='Year') +
    xlab('Step length') +
    xlim(c(0,1))
  angle <- ggplot(data, aes(x=angle, col=factor(year(dtUTC)))) +
    geom_density() +
    labs(col='Year') +
    xlab('Relative turning angles')
  ggarrange(step, angle, common.legend=T)
}





