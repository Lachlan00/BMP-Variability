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
# CTD sampled area (vertically). Not sure if I should.

#-------------#
#  Load data  #
#-------------#
# CTD data
CTD <- cast.reader()
# Aggregation data
agg <- load.agg()
# Extract agg data for each CTD cast
cast.df <- CTD.agg(CTD, agg)
cast.df$year <- str_sub(cast.df$survey, 1, 4)
# Calc mean for annual comparisons
mean.casts <- aggregate(cast.df, by=list(cast.df$survey), mean)
mean.casts[,2:4] <- NULL
names(mean.casts)[1] <- "survey"

# Load VCUR data
VCUR <- read.csv('./data/ANMN/custom/VCUR_surveys.csv')
names(VCUR)[3] <- "BMP"
VCUR$survey <- factor(VCUR$survey, levels=unique(cast.df$survey))
VCUR$year <- str_sub(VCUR$survey, 1, 4)

# VCUR
p.VCUR <- ggplot(VCUR, aes(x=survey, y=VCUR, fill=BMP)) +
  geom_boxplot() + 
  #stat_summary(fun=mean, geom="line", aes(group=1), color='darkred') +
  theme_bw() +
  xlab(NULL) +
  #ylab("Southward Geostrophic Velocity (m/s)") +
  ylab(bquote('Meridional Current Velocity ('*ms^-1*')')) +
  theme(legend.position=c(.23,.85),
        legend.background = element_rect(colour = 'black', fill='white', 
                                         linetype='solid')) +
  guides(colour = guide_legend(override.aes = list(alpha = .05))) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", size=.5) +
  labs(fill=NULL) +
  scale_x_discrete(drop = FALSE)

# Temperature
p.temp <- ggplot(cast.df, aes(x=survey, y=temp_mean, fill=year)) +
  #geom_violin(fill='lightgrey') +
  geom_boxplot() +  #, width=.1) + 
  #stat_summary(fun=mean, geom="line", aes(group=1), color='darkred') +
  theme_bw() +
  xlab(NULL) +
  ylab("Temperature (\u00B0C)") + 
  theme(legend.position = "none")

# Salinty
p.salt <- ggplot(cast.df, aes(x=survey, y=salt_mean, fill=year)) +
  #geom_violin(fill='lightgrey') +
  geom_boxplot() + #, width=.1) + 
  #stat_summary(fun=mean, geom="line", aes(group=1), color='darkred') +
  theme_bw() +
  lims(y=c(34.75, 35.75)) +
  xlab(NULL) +
  ylab("Salinity (PSU)") + 
  theme(legend.position = "none")

# Acoustics
# Swarms
p.swarms <- ggplot(cast.df, aes(x=survey, y=area_swarm_count)) +
  geom_boxplot(outlier.shape = NA, fill='lightgrey') + 
  #stat_summary(fun=mean, geom="line", aes(group=1), color='darkred') +
  theme_bw() +
  xlab(NULL) +
  #ylab("Swarms / km2") +
  ylab(bquote('Swarms / '*km^2*'')) +
  scale_y_continuous(limits = quantile(cast.df$area_swarm_count, c(0.1, 0.9)))

# Prey avliability
p.available <- ggplot(cast.df, aes(x=survey, y=area_corrected_area_sum*100)) +
  geom_boxplot(outlier.shape = NA, fill='lightgrey') + 
  #stat_summary(fun=mean, geom="line", aes(group=1), color='darkred') +
  theme_bw() +
  xlab(NULL) +
  ylab("Swarm coverage (%)") +
  scale_y_continuous(limits = quantile(cast.df$area_corrected_area_sum*100, c(0.1, 0.9)))

# mean aggregation density
p.density <- ggplot(cast.df, aes(x=survey, y=Sv_mean_mean)) +
  geom_boxplot(outlier.shape = NA, fill='lightgrey') +
  #stat_summary(fun=mean, geom="line", aes(group=1), color='darkred') +
  theme_bw() +
  xlab(NULL) +
  ylab("Mean swarm density (Sv mean)") +
  scale_y_continuous(limits = quantile(cast.df$Sv_mean_mean, c(0.1, 0.9)))

# Total biomass per unit area
p.biomass <- ggplot(cast.df, aes(x=survey, y=area_Sv_mean_sum)) +
  geom_boxplot(outlier.shape = NA, fill='lightgrey') + 
  #stat_summary(fun=mean, geom="line", aes(group=1), color='darkred') +
  theme_bw() +
  xlab(NULL) +
  ylab("Biomass (Sv mean / km2)") +
  ylab(bquote("Biomass (Sv mean / "*km^2*")")) +
  scale_y_continuous(limits = quantile(cast.df$area_Sv_mean_sum, c(0.1, 0.9)))

# Total biomass per unit area
p.size <- ggplot(cast.df, aes(x=survey, y=corrected_area_mean)) +
  geom_boxplot(outlier.shape = NA, fill='lightgrey') + 
  #stat_summary(fun=mean, geom="line", aes(group=1), color='darkred') +
  theme_bw() +
  xlab(NULL) +
  #ylab("Swarm size (m2)") +
  ylab(bquote("Mean swarm size ("*m^2*")")) +
  scale_y_continuous(limits = quantile(cast.df$corrected_area_mean, c(0.1, 0.9)))

# Mean depth of swarms
p.depth <- ggplot(cast.df, aes(x=survey, y=swarm_mean_depth)) +
  geom_boxplot(outlier.shape = NA, fill='lightgrey') + 
  #stat_summary(fun=mean, geom="line", aes(group=1), color='darkred') +
  theme_bw() +
  xlab(NULL) +
  ylab("Mean swarm depth (m)") +
  scale_y_reverse(limits = rev(quantile(cast.df$swarm_mean_depth, c(0.1, 0.9))))

p.env <- ggarrange(p.VCUR, p.temp, p.salt, ncol=1)
ggsave("./figures/boxplots/env_boxplot.png", p.env, width=6, height=10.7)

# report stats for reporting
# VCUR
message('\nVCUR')
vcur.d <- aggregate(VCUR ~ survey + BMP, VCUR, mean)
vcur.d$sd  <- aggregate(VCUR ~ survey + BMP, VCUR, sd)$VCUR
print(vcur.d)
# Temperature
message('\nTEMP')
temp.d <- aggregate(temp_mean ~ survey, cast.df, mean)
temp.d$sd <- aggregate(temp_mean ~ survey, cast.df, sd)$temp_mean
print(temp.d)
# Salinity
message('\nSALT')
salt.d <- aggregate(salt_mean ~ survey, cast.df, mean)
salt.d$sd <- aggregate(salt_mean ~ survey, cast.df, sd)$salt_mean
print(salt.d)



ggarrange(p.swarms, p.available, p.biomass, p.density, p.size, p.depth, ncol=2, nrow=3)



# Calculate correlations between these variables
cor.dat <- mean.casts[,c('temp_mean', 'salt_mean', 'area_swarm_count', 'corrected_area_mean',
                         'swarm_mean_depth','Sv_mean_mean',
                         'area_corrected_area_sum', 'area_Sv_mean_sum')]
names(cor.dat) <- c('Temperature', 'Salinity', 'Swarm count', 'Swarm size', 'Swarm depth',
                    'Swarm density', 'Swarm coverage', 'Biomass')
# Correlation plot
p.corr.inter <- ggcorrplot(cor(cor.dat),
           colors=c('#0b3487','#ffffff','#87130b'),
           lab=T, type='upper') +
  theme(legend.position = 'none')
ggsave("./figures/corr/inter_corr.png", p.corr.inter)

