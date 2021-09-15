# Acoustics / CTD GAMS
library(ggplot2)
library(cmocean)
library(geosphere)
library(viridisLite)
library(mgcv)

#------------------#
#  Pre Processing  #
#------------------#
# Read in data
agg <- read.csv('./data/surveys/acoustics/agg_CTD_interporlations.csv')
# Clean data
keep <- c('Sv_mean', 'NASC', 'Sv_max', "Sv_min", 'Corrected_length',
          'Corrected_perimeter', 'Corrected_area', 'Corrected_thickness', 'survey',
          'Lon_M', 'Lat_M', 'Date_M', 'Time_M', 'Depth_mean','CTD_temp', 'CTD_salt')
agg <- agg[, keep]
# clean
agg <- agg[agg$Sv_mean  > -998,]
agg <- agg[agg$Lon_M <= 360,]
agg <- agg[agg$Corrected_thickness  > -9999,]
agg <- agg[agg$Corrected_length  > -9999,]
agg <- agg[agg$Sv_mean  < 0,]

# Make survey factor
agg$survey <- as.factor(agg$survey)

#--------------------------#
# Add cross-shore distance #
#--------------------------#
# Find nearest transect for each point
transect_lines <- read.csv('./data/surveys/meta/transect_lines.csv')
agg$transect <- sapply(agg$Lat_M, function(lat) transect_lines$transect[which.min(abs(transect_lines$lat1 - lat))])
# Calculate crossshore distance using haversin formula
points <- transect_lines[agg$transect, c('lon1', 'lat1')]
points$lon2 <- agg$Lon_M
agg$crossShoreDist <- mapply(function(lon1, lon2, lat) distHaversine(c(lon1, lat), c(lon2, lat)),
       points$lon1, points$lon2, points$lat1)
# convert to km
agg$crossShoreDist <- agg$crossShoreDist/1000

#------------#
#    GAMS    #
#------------#
# Thanks Martin for looking at this!
# GAMS to go here
head(agg)
m1=gam(Sv_mean~s(Depth_mean,k=5,bs='cr',by=survey)+survey,data=agg)
summary(m1)
plot(m1)

#-------------#
#    Plots    #
#-------------#
hist(agg$Sv_mean)
# Note: looking at  this plot it's clear there's a lot of duplicate points (by T and S)
# due to the resolution of the the interporlation
ggplot(agg, aes(x=CTD_salt, y=CTD_temp,colour=Sv_mean)) +
  geom_point(alpha=.2) +
  facet_wrap(~survey) +
  scale_colour_viridis_c()
  #scale_fill_cmocean(name='thermal')




