# Vertical correlations
library(ggcorrplot)
library(plyr)
source('R/acoustics_functions.R')

# CTD data
CTD <- cast.reader()
# Aggregation data
agg <- load.agg()
surveys <- unique(CTD$data$survey_id)

# round both by depth at 5m bins
CTD.depth.bin <- CTD$data[,c('survey_id','depth','temperature','salinity')]
CTD.depth.bin$depth <- round_any(CTD.depth.bin$depth, 5)
CTD.depth.bin <- CTD.depth.bin[order(CTD.depth.bin$survey_id, CTD.depth.bin$depth),]

agg <- agg[,c('survey','Depth_mean','Corrected_area','Sv_mean')]
agg$Depth_mean <- round_any(agg$Depth_mean, 5)
names(agg)[1:2] <- c('survey_id', 'depth')
agg <- agg[order(agg$survey_id, agg$depth),]

# eliminate values that are two deep so we have equal stuff to deal with
agg.split <- split(agg, agg$survey_id)
CTD.depth.bin.split <- split(CTD.depth.bin, CTD.depth.bin$survey_id)
for (i in 1:length(surveys)){
  agg.split[[i]]$depth[agg.split[[i]]$depth > max(CTD.depth.bin.split[[i]]$depth)] <- max(CTD.depth.bin.split[[i]]$depth)
  CTD.depth.bin.split[[i]]$depth[CTD.depth.bin.split[[i]]$depth > max(agg.split[[i]]$depth)] <- max(agg.split[[i]]$depth)
  agg.split[[i]]$depth[agg.split[[i]]$depth < min(CTD.depth.bin.split[[i]]$depth)] <- min(CTD.depth.bin.split[[i]]$depth)
  CTD.depth.bin.split[[i]]$depth[CTD.depth.bin.split[[i]]$depth < min(agg.split[[i]]$depth)] <- min(agg.split[[i]]$depth)
}
agg <- do.call(rbind, agg.split)
row.names(agg) <- NULL
CTD.depth.bin <- do.call(rbind, CTD.depth.bin.split)
row.names(CTD.depth.bin) <- NULL


# calc agg statistcis
CTD.depth.bin <- aggregate(. ~ survey_id + depth, 
                           data=CTD.depth.bin, mean)
agg_sum <- aggregate(. ~ survey_id + depth, data=agg, sum)
agg_mean <-  aggregate(. ~ survey_id + depth, data=agg, mean)

# Calculate correlations
p.ls1 <- list()
p.ls2 <- list()
for (i in 1:length(surveys)){
  # check depths match
  if (!all(agg_sum$depth[agg_sum$survey_id == surveys[i]] == CTD.depth.bin$depth[CTD.depth.bin$survey_id == surveys[i]])){
    print("Oh SNAP!")
    break
  }
  df.corr <- data.frame(depth=agg_sum$depth[agg_sum$survey_id == surveys[i]],
                        temperature = CTD.depth.bin$temperature[CTD.depth.bin$survey_id == surveys[i]],
                        salinity = CTD.depth.bin$salinity[CTD.depth.bin$survey_id == surveys[i]],
                        sum_desnity = agg_sum$Sv_mean[agg_sum$survey_id == surveys[i]],
                        mean_density = agg_mean$Sv_mean[agg_mean$survey_id == surveys[i]])
  p.ls1[[i]] <- ggcorrplot(cor(df.corr),
            colors=c('#0b3487','#ffffff','#87130b'),
            lab=T, type='upper', p.mat=cor_pmat(df.corr)) +
    theme(legend.position = 'none') +
    ggtitle(surveys[i])
  p.ls2[[i]] <- ggplot(df.corr, aes(x=temperature, y=mean_density, col=depth)) +
    geom_smooth(col='darkgrey') +
    geom_point(size=3) + ggtitle(surveys[i]) +
    scale_color_viridis(direction=-1)
   
  
}
ggarrange(plotlist=p.ls1)
ggarrange(plotlist=p.ls2)





# ggplot(agg, aes(x=Sv_mean)) +
#   geom_density() +
#   facet_wrap(~survey_id)

