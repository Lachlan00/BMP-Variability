# Make plots of data 
library(lubridate)
library(geosphere)
library(cmocean)
library(ggplot2)
library(oce)
library(ggcorrplot)
source('R/CTD.R')
source('R/eudyptula.R')

# CTD data
CTD <- cast.reader()
# Check CTD profile
CTD.plot.cast(CTD, castID='CC1449004_20160917_041449')

# Aggregation data
agg <- readRDS('data/surveys/acoustics/agg.rds')
# Create datetimes
agg$dt <- as.POSIXct(paste0(agg$Date_M,' ',agg$Time_M), format="%Y%m%d %H:%M:%OS", tz='UTC')
# Clean bad values
agg <- agg[agg$NASC > 0,]
agg <- agg[agg$Lon_M <= 360,]
agg <- agg[agg$Corrected_area > 0,]


# Extract agg data for each CTD cast
n <- rep(NA, nrow(CTD$meta))
cast.df <- data.frame(survey=n,
                      cast_id=n,
                      station=n,
                      cast_time_UTC=n, 
                      temp_min=n, temp_mean=n, temp_max=n,
                      salt_min=n, salt_mean=n, salt_max=n,
                      NASC_min=n, NASC_mean=n, NASC_max=n,
                      Sv_mean=n)
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
  cast.df[i, c('survey', 'cast_id', 'station', 'cast_time_UTC')] <- 
    cast[, c('survey_id', 'file_name', 'transect_id', 'cast_time_UTC')]
  # Cast data
  cast.df[i, c('temp_min', 'temp_mean', 'temp_max', 'salt_min', 'salt_mean', 'salt_max')] <-
    c(min(cast.data$temperature), mean(cast.data$temperature), max(cast.data$temperature),
      min(cast.data$salinity), mean(cast.data$salinity), max(cast.data$salinity))
  # Aggregation data
  cast.df[i, c('NASC_min', 'NASC_mean', 'NASC_max','Sv_mean', 'Sv_mean_sum', 'correct_area_sum')] <-
    c(min(cast.agg$NASC), mean(cast.agg$NASC), max(cast.agg$NASC), 
      Sv_mean.log(mean(Sv_mean.linear(cast.agg$Sv_mean))),
      Sv_mean.log(sum(Sv_mean.linear(cast.agg$Sv_mean))),
      sum(cast.agg$Corrected_area))
  # Progress
  setTxtProgressBar(pb, i)
}
close(pb)

# Filter casrs with no data
cast.df <- cast.df[complete.cases(cast.df),]

# Remove outliers
plot(cast.df$NASC_mean)
abline(h=7e6, col='red')
cast.df <- cast.df[cast.df$NASC_mean < 7e6,]

# Plot results
#--------------------------#
#       Mean Sv mean       # 
#--------------------------#
p <- ggplot(cast.df, aes(x=salt_mean, y=temp_mean, col=Sv_mean)) +
  geom_point(size=3, alpha=0.8) +
  scale_color_viridis() +
  facet_wrap(~survey) +
  theme_bw() +
  labs(x="Salinity (PSU)", 
       y="Temperature (\u00B0C)", 
       col=expression(S[v]~mean)) +
  ggtitle('T/S and mean swarm Sv_mean')

ggsave('output/TS/TS_mean_swarm_Sv.png', p)

# Just 2015/2018
p <- ggplot(cast.df[cast.df$survey %in% c('2015_S1', '2018_S1'),],
            aes(x=salt_mean, y=temp_mean, col=Sv_mean)) +
  geom_point(size=3, alpha=0.8) +
  scale_color_viridis() +
  facet_wrap(~survey) +
  theme_bw() +
  labs(x="Salinity (PSU)", 
       y="Temperature (\u00B0C)", 
       col=expression(S[v]~mean)) +
  ggtitle('T/S and mean swarm Sv_mean')

ggsave('output/TS/Extreme_TS_mean_swarm_Sv.png', p)

# Just 2015/2018 - contour version
p <- ggplot(cast.df[cast.df$survey %in% c('2015_S1', '2018_S1'),],
            aes(x=salt_mean, y=temp_mean, col=Sv_mean)) +
  geom_point(size=3, alpha=0.8) +
  scale_color_viridis() +
  facet_wrap(~survey) +
  theme_bw() +
  labs(x="Salinity (PSU)", 
       y="Temperature (\u00B0C)", 
       col=expression(S[v]~mean)) +
  ggtitle('T/S and mean swarm Sv_mean')

ggsave('output/TS/Extreme_TS_mean_swarm_Sv.png', p)

#--------------------------#
#    Cumulative Sv mean    #
#--------------------------#
p <- ggplot(cast.df, aes(x=salt_mean, y=temp_mean, col=Sv_mean_sum)) +
  geom_point(size=3, alpha=0.8) +
  scale_color_viridis() +
  facet_wrap(~survey) +
  theme_bw() +
  labs(x="Salinity (PSU)",
       y="Temperature (\u00B0C)", 
       col=expression(S[v]~mean)) +
  ggtitle('T/S and cumulative swarm Sv_mean')

ggsave('output/TS/TS_cumulative_swarm_Sv.png', p)

# Just 2015/2018
p <- ggplot(cast.df[cast.df$survey %in% c('2015_S1', '2018_S1'),],
            aes(x=salt_mean, y=temp_mean, col=Sv_mean_sum)) +
  geom_point(size=3, alpha=0.8) +
  scale_color_viridis() +
  facet_wrap(~survey) +
  theme_bw() +
  labs(x="Salinity (PSU)",
       y="Temperature (\u00B0C)", 
       col=expression(S[v]~mean)) +
  ggtitle('T/S and cumulative swarm Sv_mean')

ggsave('output/TS/Extreme_TS_cumulative_swarm_Sv.png', p)

#----------------------#
#    Corrected Area    #
#----------------------#
p <- ggplot(cast.df, aes(x=salt_mean, y=temp_mean, col=correct_area_sum)) +
  geom_point(size=3, alpha=0.8) +
  scale_color_viridis() +
  facet_wrap(~survey) +
  theme_bw() +
  labs(x="Salinity (PSU)",
       y="Temperature (\u00B0C)", 
       col=expression(S[v]~mean)) +
  ggtitle('T/S and cumulative swarm corrected_area')

ggsave('output/TS/TS_cumulative_swarm_corrected_area.png', p)

# Just 2015/2018
p <- ggplot(cast.df[cast.df$survey %in% c('2015_S1', '2018_S1'),],
            aes(x=salt_mean, y=temp_mean, col=correct_area_sum)) +
  geom_point(size=3, alpha=0.8) +
  scale_color_viridis() +
  facet_wrap(~survey) +
  theme_bw() +
  labs(x="Salinity (PSU)",
       y="Temperature (\u00B0C)", 
       col=expression(S[v]~mean)) +
  ggtitle('T/S and cumulative swarm corrected_area')

ggsave('output/TS/Extreme_TS_cumulative_swarm_corrected_area.png', p)

#--------------------#
#    CTD Boxplots    #
#--------------------#
# Temperature
ggplot(CTD$data, aes(x=survey_id, y=temperature)) +
  geom_boxplot() + 
  theme_bw()

ggplot(CTD$data[CTD$data$survey_id %in% c('2015_S1', '2018_S1'),], aes(x=survey_id, y=temperature)) +
  geom_boxplot() + 
  theme_bw()

# Salinty
ggplot(CTD$data, aes(x=survey_id, y=salinity)) +
  geom_boxplot() + 
  theme_bw()

ggplot(CTD$data[CTD$data$survey_id %in% c('2015_S1', '2018_S1'),], aes(x=survey_id, y=salinity)) +
  geom_boxplot() + 
  theme_bw()


#--------------------------#
#    Correlation Matrix    #
#--------------------------#
cordf_2015_S1 <- cast.df[cast.df$survey == '2015_S1',
                         c('temp_mean', 'salt_mean', 'Sv_mean')]
cormat_2015_S1 <- cor(cordf_2015_S1)
cordf_2018_S1 <- cast.df[cast.df$survey == '2018_S1',
                         c('temp_mean', 'salt_mean', 'Sv_mean')]
cormat_2018_S1 <- cor(cordf_2018_S1)

ggcorrplot(cormat_2015_S1, colors=c('#0b3487','#ffffff','#87130b'), lab=T, type='upper', title='2015_S1')
ggcorrplot(cormat_2018_S1, colors=c('#0b3487','#ffffff','#87130b'), lab=T, type='upper', title='2018_S1')


#------------------------#
#    Swarm Histograms    #
#------------------------#
ggplot(agg[agg$survey %in% c('2015_S1', '2018_S1'),],
       aes(x=Sv_mean, col=survey, fill=survey)) +
  geom_density(alpha=.6) +
  xlim(-70, -25) +
  theme_bw()

ggplot(agg[agg$survey %in% c('2015_S1', '2018_S1'),],
       aes(x=Corrected_area, col=survey, fill=survey)) +
  geom_density(alpha=.6) +
  xlim(-10, 50) +
  theme_bw()

