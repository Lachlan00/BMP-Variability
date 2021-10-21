# Prey field analysis
library(lubridate)
library(geosphere)
library(cmocean)
library(ggplot2)
library(oce)
library(ggcorrplot)
library(gridExtra)
source('R/CTD.R')
source('R/eudyptula.R')

# __ GOALS __
# Calculate T/S and acoustic metrics
# Make CTD correlation posts

# __ Settings __
# filter aggregations below max cast depth?
filter.agg.depth = TRUE

#------------------#
#  Pre Processing  #
#------------------#
# CTD data
CTD <- cast.reader()

# Load transect station data
stations <- read.csv('data/surveys/meta/transect_df.csv') 
# make list of firt and last station IDs
end_points <- split(stations, stations$transect)
end_points <- lapply(end_points, function(df) df[c(1, nrow(df)),])
end_points <- do.call(rbind, end_points)
end_points <- end_points$id

# Aggregation data
agg <- readRDS('data/surveys/acoustics/agg.rds')
# Create datetimes
agg$dt <- as.POSIXct(paste0(agg$Date_M,' ',agg$Time_M), format="%Y%m%d %H:%M:%OS", tz='UTC')
# Clean bad values
agg <- agg[agg$NASC > 0,]
agg <- agg[agg$Lon_M <= 360,]
agg <- agg[agg$Corrected_area > 0,]
agg <- agg[agg$Lon_M <= 360,]
agg <- agg[agg$Corrected_thickness  > -9999,]
agg <- agg[agg$Corrected_length  > -9999,]
agg <- agg[agg$Sv_mean < 0,]
agg <- agg[agg$Sv_mean > -999,]

# Extract agg data for each CTD cast
n <- rep(NA, nrow(CTD$meta))
cast.df <- data.frame(survey=n,
                      cast_id=n,
                      station_id=n,
                      cast_time_UTC=n,
                      max_depth=n, swarm_count=n,
                      temp_min=n, temp_mean=n, temp_max=n,
                      salt_min=n, salt_mean=n, salt_max=n,
                      NASC_min=n, NASC_mean=n, NASC_max=n,
                      Sv_mean=n, Sv_mean_sum=n, 
                      corrected_area_mean=n, corrected_area_sum=n,
                      CTD_area=n)
message('Extracting aggregation data for casts..')
pb <- txtProgressBar(min = 0, max = length(cast.df), style = 3)
for (i in 1:nrow(CTD$meta)){
  cast <- CTD$meta[i,]
  cast.data <- CTD$data[CTD$data$id == cast$file_name,]
  # To search efficiently, first filter survey, then time (2 hours) and then GPS range
  # Survey
  cast.agg <- agg[agg$survey == cast$survey_id,]
  if (nrow(cast.agg) == 0) next
  cast.agg <- cast.agg[cast.agg$dt > (cast$cast_time_UTC - hours(2)) &
                         cast.agg$dt < (cast$cast_time_UTC + hours(2)),]
  if (filter.agg.depth){
    cast.agg <- cast.agg[cast.agg$Depth_mean <= ceiling(max(cast.data$depth)),]
  }
  # Calc distance to cast
  if (nrow(cast.agg) == 0) next
  cast.agg$distance <- mapply(function(lon1, lat1, lon2, lat2) 
    distHaversine(c(lon1,lat1), c(lon2,lat2)),
    cast.agg$Lon_M, cast.agg$Lat_M, 
    cast$start_longitude, cast$start_latitude)
  # Filter distance 
  cast.agg <- cast.agg[cast.agg$distance <= 1875,]
  if (nrow(cast.agg) == 0) next
  # Populate data frame
  # Cast meta
  cast.df[i, c('survey', 'cast_id', 'station_id', 'cast_time_UTC')] <- 
    cast[, c('survey_id', 'file_name', 'transect_id', 'cast_time_UTC')]
  cast.df[i, 'max_depth'] <- ceiling(max(cast.data$depth))
  cast.df[i, 'swarm_count'] <- nrow(cast.agg)
  # Cast data
  cast.df[i, c('temp_min', 'temp_mean', 'temp_max', 'salt_min', 'salt_mean', 'salt_max')] <-
    c(min(cast.data$temperature), mean(cast.data$temperature), max(cast.data$temperature),
      min(cast.data$salinity), mean(cast.data$salinity), max(cast.data$salinity))
  # Aggregation data
  cast.df[i, c('NASC_min', 'NASC_mean', 'NASC_max','Sv_mean', 'Sv_mean_sum', 
               'corrected_area_mean', 'corrected_area_sum')] <-
    c(min(cast.agg$NASC), mean(cast.agg$NASC), max(cast.agg$NASC), 
      Sv_mean.log(mean(Sv_mean.linear(cast.agg$Sv_mean))),
      Sv_mean.log(sum(Sv_mean.linear(cast.agg$Sv_mean))),
      mean(cast.agg$Corrected_area), sum(cast.agg$Corrected_area))
  # Calculate CTD area
  xdist <- ifelse(cast$transect_id %in% end_points, 1.875, 1.875*2)
  cast.df[i, 'CTD_area'] <- ceiling(max(cast.data$depth))*(xdist*1000)
  
  # Progress
  setTxtProgressBar(pb, i)
}
close(pb)

# Filter casrs with no data
cast.df <- cast.df[complete.cases(cast.df),]

# Calculate acoustic energy per unit area
# chnage to km2
cast.df$CTD_area <- cast.df$CTD_area/1000
# swarms per unit area
cast.df$area_swarms <- cast.df$swarm_count / cast.df$CTD_area
# corrected area per unit area
cast.df$area_corrected_area_sum <- cast.df$corrected_area_sum / cast.df$CTD_area
cast.df$area_corrected_area_mean <- cast.df$corrected_area_mean / cast.df$CTD_area
# Sv mean per unit area
cast.df$area_Sv_mean <- Sv_mean.log(Sv_mean.linear(cast.df$Sv_mean)/cast.df$CTD_area)
cast.df$area_Sv_mean_sum <- Sv_mean.log(Sv_mean.linear(cast.df$Sv_mean_sum)/cast.df$CTD_area)


#------------#
#  Plotting  #
#------------#
cast.df$transect <- as.numeric(substr(cast.df$station_id, 2, 2))
cast.df$station <- as.numeric(substr(cast.df$station_id, 4, 4))
############
# Transect #
############
vars = c('area_swarms', 'area_corrected_area_sum', 'area_corrected_area_mean',
         'area_Sv_mean', 'area_Sv_mean_sum')
for (var in vars){
  cast.df['var'] = cast.df[var]
  p <- ggplot(cast.df, aes(x=temp_mean, y=var, col=as.factor(transect))) +
    geom_point(size=3) +
    facet_wrap(~survey) +
    theme_bw() +
    labs(x="Temperature (\u00B0C)",
         y=var) +
    ggtitle(var)
  print(p)
}

#############
#  Boxplot  #
#############
vars = c('area_swarms', 'area_corrected_area_sum', 'area_corrected_area_mean',
         'area_Sv_mean', 'area_Sv_mean_sum')
for (var in vars){
  cast.df['var'] = cast.df[var]
  p <- ggplot(cast.df, aes(x=as.factor(survey), y=var)) +
    geom_boxplot() +
    theme_bw() +
    labs(x="Survey",
         y=var) +
    ggtitle(var)
  ggsave(paste0('output/temp/boxplots_',var,'.png'), p, width=10, height=3)
}

######
# TS #
######
# for (var in vars){
#   cast.df['var'] = cast.df[var]
#   p <- ggplot(cast.df, aes(x=salt_mean, y=temp_mean, col=var)) +
#     geom_point(size=3) +
#     facet_wrap(~survey) +
#     theme_bw() +
#     scale_color_viridis() +
#     labs(x="Salinity (PSU)",
#          y="Temperature (\u00B0C)") +
#     ggtitle(var)
#   ggsave(paste0('output/temp/agg_',var,'.png'), p)
# }

####################
# area Sv mean sum #
####################
surveys <- unique(cast.df$survey)
for (i in 1:length(surveys)){
  p <- ggplot(cast.df[cast.df$survey == surveys[i],],
              aes(x=salt_mean, y=temp_mean,
                  col=area_Sv_mean_sum, z=area_Sv_mean)) +
    geom_density2d(col='grey', alpha=.4, size=0.5) +
    geom_point(size=1.5) +
    facet_wrap(~survey) +
    theme_bw() +
    scale_color_viridis(limits=c(min(cast.df$area_Sv_mean),
                                 max(cast.df$area_Sv_mean))) +
    lims(x=c(min(cast.df$salt_mean), max(cast.df$salt_mean)),
         y=c(min(cast.df$temp_mean), max(cast.df$temp_mean))) +
    labs(x="Salinity (PSU)",
         y="Temperature (\u00B0C)",
         col='Sv mean m2')
  ggsave(paste0('output/TS/',surveys[i],'_Sv_mean.png'), p,
         width=5, height=4)
}

################
# Correlations #
################
cor.df <- cast.df[,c('survey', 'temp_mean', 'salt_mean', 'Sv_mean')]
cor.ls <- split(cor.df, cor.df$survey)
cor.ls <- lapply(cor.ls, function(df) cor(df[,c('temp_mean', 'salt_mean', 'Sv_mean')]))
p.ls <- list()
for (i in 1:length(cor.ls)){
  p.ls[[i]] <- ggcorrplot(cor.ls[[i]],
             colors=c('#0b3487','#ffffff','#87130b'),
             lab=T, type='upper', title=names(cor.ls[i])) +
    theme(legend.position = 'none')
}
p <- do.call("grid.arrange", c(p.ls, ncol=floor(sqrt(length(p.ls)))))
ggsave('output/correlations/corr_all.png', p, height=12, width=12, bg='white')

# plot each
names(p.ls) <- unique(cor.df$survey)
for (i in 1:length(p.ls)){
  ggsave(paste0('output/correlations/corr_',names(p.ls[i]),'.png'), p.ls[[i]], bg='white',
         width=4, height=4)
}

