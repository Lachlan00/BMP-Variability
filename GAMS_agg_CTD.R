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
# data path - specific to Martin's laptop
dat_dir='/home/martin/Documents/lac_dat/surveys/'
#agg <- read.csv('./data/surveys/acoustics/agg_CTD_interporlations.csv')
agg <- read.csv(paste0(dat_dir,'acoustics/agg_CTD_interporlations.csv'))

transect_lines <- read.csv(paste0(dat_dir,'meta/transect_lines.csv'))

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
#transect_lines <- read.csv('./data/surveys/meta/transect_lines.csv')
#
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
# Martin: the process for this is as follows:
# 0. use AIC for model selection
# 1.fit a null model - make sure we are doing better than a mean
# 2. check the influence of temperature and salinity
# 3.  check model performance
# 4. with the 'best' model from step. 2 see if position is important.

#null:
m0=gam(Sv_mean~1+survey,data=agg)
#temperature:
m1=gam(Sv_mean~s(CTD_temp,k=5,bs='cr'),data=agg)
#salinty:
m2=gam(Sv_mean~s(CTD_salt,k=5,bs='cr'),data=agg)
#temperature and salinity:
m3=gam(Sv_mean~s(CTD_temp,k=5,bs='cr')+s(CTD_salt,k=5,bs='cr'),data=agg)
#check AIC:
AIC_v=sapply(list(m0,m1,m2,m3),function(x) AIC(x))
which.min(AIC_v) #m3 is the 'best' one:

#model diganostics for m3
hist(residuals(m3,tyep='deviance')) #not great
qq.gam(m3) #not great either
plot(m3$model$Sv_mean,fitted(m3)) #pretty ordinary

#lets try adding in survey:
m4=gam(Sv_mean~s(CTD_temp,k=5,bs='cr',by=survey)+
         s(CTD_salt,k=5,bs='cr',by=survey)+survey,data=agg)
AIC_too_v=sapply(list(m3,m4),function(x) AIC(x))

#model diagnostics for m4
hist(residuals(m4,type='deviance')) #not great - this suggests we are missing an explanatory variable
qq.gam(m4) #not great either
plot(m4$model$Sv_mean,fitted(m4)) #little better
plot(m4) #salinity doesn't look like its fitting very well
summary(m4)
AIC(m4)

m5=gam(Sv_mean~s(CTD_temp,k=5,bs='cr',by=survey)+
         s(CTD_salt,bs='cr',by=survey)+survey,data=agg)

#not too happy with the salinity modelling

#


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




