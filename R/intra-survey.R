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

#-------------#
#  Load data  #
#-------------#
# CTD datas
CTD <- cast.reader()
# Aggregation data
agg <- load.agg()
# Extract agg data for each CTD cast
cast.df <- CTD.agg(CTD, agg)

# Prey availability
p.available <- ggplot(cast.df, aes(x=salt_mean, y=temp_mean, col=area_corrected_area_sum)) +
  geom_point() + 
  theme_bw() +
  xlab("Salinity") +
  ylab("Temperature") +
  facet_wrap(~survey) +
  scale_color_viridis() +
  ggtitle("Prey availability")

# Mean aggregation density
p.density <- ggplot(cast.df, aes(x=salt_mean, y=temp_mean, col=Sv_mean_mean)) +
  geom_point() + 
  theme_bw() +
  xlab("Salinity") +
  ylab("Temperature") +
  facet_wrap(~survey) +
  scale_color_viridis() +
  ggtitle("Swarm density")

# Total biomass per unit area
p.biomass <- ggplot(cast.df, aes(x=salt_mean, y=temp_mean, col=area_Sv_mean_sum)) +
  geom_point() + 
  theme_bw() +
  xlab("Salinity") +
  ylab("Temperature") +
  facet_wrap(~survey) +
  scale_color_viridis() +
  ggtitle("Biomass (Sv mean / km2)")

ggarrange(p.biomass, p.density, legend='none')

# Calculate correlations between variables
# cor.dat <- cast.df[,c('survey', 'temp_mean', 'salt_mean', 'area_swarm_count',
#                       'corrected_area_mean', 'Sv_mean_mean', 'area_corrected_area_sum', 'area_Sv_mean_sum')]
# names(cor.dat) <- c('Survey', 'Temperature', 'Salinity', 'Swarm Count',
#                     'Swarm Size', 'Swarm Density', 'Prey Availability', 'Biomass')
cor.dat <- cast.df[,c('survey',
                      'temp_mean', 'salt_mean', 'area_swarm_count', 'corrected_area_mean',
                         'swarm_mean_depth','Sv_mean_mean',
                         'area_corrected_area_sum', 'area_Sv_mean_sum')]
names(cor.dat) <- c('Survey',
                    'Temperature', 'Salinity', 'Swarm count', 'Swarm size', 'Swarm depth',
                    'Swarm density', 'Prey availability', 'Biomass')
# cor.dat <- cast.df[,c('survey', 'temp_mean', 'salt_mean', 'area_Sv_mean_sum')]
# names(cor.dat) <- c('Survey', 'Temperature', 'Salinity', 'Biomass')

p.ls <- list()
i <- 1
for (survey in unique(cor.dat$Survey)){
  p <- ggcorrplot(cor(cor.dat[cor.dat$Survey == survey, 2:ncol(cor.dat)]),
             colors=c('#0b3487','#ffffff','#87130b'),
             lab=T, type='upper', lab_size=3) +
    theme(legend.position = 'none',
          axis.text.x = element_text(size=9),
          axis.text.y = element_text(size=9)) +
    ggtitle(survey)
  p.ls[[i]] <- p
  #ggsave(paste0('./output/correlations/corr_',survey,'.png'), p)
  i <- i + 1
}
p.corr.intra <- ggarrange(plotlist=p.ls)
ggsave("./figures/corr/intra_corr.png", p.corr.intra, width=14, height=14)


