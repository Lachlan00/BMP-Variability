# Make plots of data 
library(lubridate)
library(geosphere)
library(cmocean)
library(ggplot2)
library(oce)
library(ggcorrplot)
library(TMB)
library(reshape2)
source('R/acoustics_functions.R')

# NOTE:
# I did not filter aggregation data to only include
# CTD sampled area. Not sure if I should.

#-------------#
#  Load data  #
#-------------#
# CTD datas
CTD <- cast.reader()
# Aggregation data
agg <- load.agg()
# Extract agg data for each CTD cast
cast.df <- CTD.agg(CTD, agg)
# Calc mean for annual comparisons
mean.casts <- aggregate(cast.df, by=list(cast.df$survey), mean)
mean.casts[,2:4] <- NULL
names(mean.casts)[1] <- "survey"

# Temperature
p.temp <- ggplot(cast.df, aes(x=survey, y=temp_mean)) +
  geom_boxplot(outlier.shape = NA, fill='lightgrey') + 
  stat_summary(fun=mean, geom="line", aes(group=1), color='darkred') +
  theme_bw() +
  xlab(NULL) +
  ylab("Temperature (\u00B0C)")

# Salinty
p.salt <- ggplot(cast.df, aes(x=survey, y=salt_mean)) +
  geom_boxplot(outlier.shape = NA, fill='lightgrey') + 
  stat_summary(fun=mean, geom="line", aes(group=1), color='darkred') +
  theme_bw() +
  lims(y=c(34.75, 35.75)) +
  xlab(NULL) +
  ylab("Salinity (PSU)")

# Acoustics
# Swarms
p.swarms <- ggplot(cast.df, aes(x=survey, y=area_swarm_count)) +
  geom_boxplot(outlier.shape = NA, fill='lightgrey') + 
  stat_summary(fun=mean, geom="line", aes(group=1), color='darkred') +
  theme_bw() +
  xlab(NULL) +
  ylab("Swarms / km2") +
  scale_y_continuous(limits = quantile(cast.df$area_swarm_count, c(0.1, 0.9)))

# Prey avliability
p.available <- ggplot(cast.df, aes(x=survey, y=area_corrected_area_sum*100)) +
  geom_boxplot(outlier.shape = NA, fill='lightgrey') + 
  stat_summary(fun=mean, geom="line", aes(group=1), color='darkred') +
  theme_bw() +
  xlab(NULL) +
  ylab("Prey availablility (%)") +
  scale_y_continuous(limits = quantile(cast.df$area_corrected_area_sum*100, c(0.1, 0.9)))

# mean aggregation density
p.density <- ggplot(cast.df, aes(x=survey, y=Sv_mean_mean)) +
  geom_boxplot(outlier.shape = NA, fill='lightgrey') +
  stat_summary(fun=mean, geom="line", aes(group=1), color='darkred') +
  theme_bw() +
  xlab(NULL) +
  ylab("Swarm density (Sv mean)") +
  scale_y_continuous(limits = quantile(cast.df$Sv_mean_mean, c(0.1, 0.9)))

# Total biomass per unit area
p.biomass <- ggplot(cast.df, aes(x=survey, y=area_Sv_mean_sum)) +
  geom_boxplot(outlier.shape = NA, fill='lightgrey') + 
  stat_summary(fun=mean, geom="line", aes(group=1), color='darkred') +
  theme_bw() +
  xlab(NULL) +
  ylab("Biomass (Sv mean / km2)") +
  scale_y_continuous(limits = quantile(cast.df$area_Sv_mean_sum, c(0.1, 0.9)))

# Total biomass per unit area
p.size <- ggplot(cast.df, aes(x=survey, y=corrected_area_mean)) +
  geom_boxplot(outlier.shape = NA, fill='lightgrey') + 
  stat_summary(fun=mean, geom="line", aes(group=1), color='darkred') +
  theme_bw() +
  xlab(NULL) +
  ylab("Swarm size (m2)") +
  scale_y_continuous(limits = quantile(cast.df$corrected_area_mean, c(0.1, 0.9)))

# Mean depth of swarms
p.depth <- ggplot(cast.df, aes(x=survey, y=swarm_mean_depth)) +
  geom_boxplot(outlier.shape = NA, fill='lightgrey') + 
  stat_summary(fun=mean, geom="line", aes(group=1), color='darkred') +
  theme_bw() +
  xlab(NULL) +
  ylab("Mean swarm depth (m)") +
  scale_y_reverse(limits = rev(quantile(cast.df$swarm_mean_depth, c(0.1, 0.9))))

ggarrange(p.temp, p.salt, ncol=1)
ggarrange(p.swarms, p.available, p.biomass, p.density, p.size, p.depth, ncol=2, nrow=3)

# Calculate correlations between these variables
cor.dat <- mean.casts[,c('temp_mean', 'salt_mean', 'area_swarm_count', 'corrected_area_mean',
                         'swarm_mean_depth','Sv_mean_mean',
                         'area_corrected_area_sum', 'area_Sv_mean_sum')]
names(cor.dat) <- c('Temperature', 'Salinity', 'Swarm count', 'Swarm size', 'Swarm depth',
                    'Swarm density', 'Prey availability', 'Biomass')
# Correlation plot
ggcorrplot(cor(cor.dat),
           colors=c('#0b3487','#ffffff','#87130b'),
           lab=T, type='upper') +
  theme(legend.position = 'none')

