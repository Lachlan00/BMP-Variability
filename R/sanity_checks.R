# Sanity check
setwd("/Users/lachlanphillips/Development/PhD/paper-repos/BMP-Variability")
source('R/CTD.R')
source('R/acoustics_functions.R')
library(reshape2)

#'------------------------------#
#' __1. Load CTD and Agg data__ #
#'------------------------------#
# CTD data
CTD <- cast.reader()
# Aggregation data
agg <- load.agg()

# Extract agg data for each CTD cast
cast.df <- CTD.agg(CTD, agg)
data.all <- cast.df[,c('survey', 'temp_mean', 'salt_mean',
                   'area_corrected_area_sum', 'Sv_mean_mean', 'swarm_mean_depth')]
data.all <- melt(data.all, id.vars='survey')
data.all$group = 'all'

# Extract agg data for each CTD cast
cast.df <- CTD.agg(CTD, agg, max.depth=30)
data.30 <- cast.df[,c('survey', 'temp_mean', 'salt_mean',
                   'area_corrected_area_sum', 'Sv_mean_mean', 'swarm_mean_depth')]
data.30 <- melt(data.30, id.vars='survey')
data.30$group = '30m'
#'------------------#
#' __1. Box Plots__ #
#'------------------#
data <- rbind(data.all, data.30)

ggplot(data, aes(x=survey, y=value, fill=group)) +
  geom_boxplot() +
  facet_wrap(~variable, scales='free') +
  theme_bw()

