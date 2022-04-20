# Make plots of data 
library(lubridate)
library(geosphere)
library(cmocean)
library(ggplot2)
library(oce)
library(ggcorrplot)
library(TMB)
library(reshape2)
library(ncdf4)
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
# attach lon/lat data
cast.df$lon <- sapply(cast.df$cast_id, function(id) CTD$meta$start_longitude[CTD$meta$file_name == id])
cast.df$lat <- sapply(cast.df$cast_id, function(id) CTD$meta$start_latitude[CTD$meta$file_name == id])

cast.df <- aggregate(cast.df, by=list(paste0(cast.df$survey,'_',cast.df$station_id)), mean)
cast.df$survey <- substr(cast.df$Group.1, 1, 7)

# Load in coastline shapefile
fn = './assets/shp/NSW-coastline/edited/polygon/NSW-coastline_WGS84'
landShp <- readOGR(dsn=path.expand(str_sub(fn,1,-(nchar(basename(fn)))-1)), 
                   layer=basename(fn))

# make plots 

# ------------- CTD ------------- #
# Raw agg
xlim <- c(150.04, 150.275)
ylim <- c(-36.47, -36.15)
maxsize <- 7
minsize <- 1.2

# Swarmcount
p.coverage <- ggplot(cast.df, aes(x=lon, y=lat, col=temp_mean)) +
  # coastline
  geom_polygon(data=landShp, aes(x=long, y=lat, group=group), 
               fill='grey', col='#636363', size=.6) +
  # Bubble
  geom_point(shape=1, size=make_normal(cast.df$area_corrected_area_sum, minsize, maxsize)) +
  
  coord_map(xlim=xlim,  ylim=ylim) +
  theme_bw() + grids(linetype = "dashed") +
  #scale_color_cmocean(name='thermal') +
  facet_wrap(~survey) +
  #theme(legend.position = 'none') +
  ggtitle("Swarm Coverage (%)") +
  labs(x="Longitude", y="Latitude") +
  scale_color_cmocean(name='thermal')

# Mean depth
p.depth <- ggplot(cast.df, aes(x=lon, y=lat, col=temp_mean)) +
  # coastline
  geom_polygon(data=landShp, aes(x=long, y=lat, group=group), 
               fill='grey', col='#636363', size=.6) +
  # Bubble
  geom_point(shape=1, size=make_normal(cast.df$swarm_mean_depth, minsize, maxsize)) +
  coord_map(xlim=xlim,  ylim=ylim) +
  theme_bw() + grids(linetype = "dashed") +
  #scale_color_cmocean(name='thermal') +
  facet_wrap(~survey) +
  #theme(legend.position = 'none') +
  ggtitle("Mean depth") +
  labs(x="Longitude", y="Latitude") +
  scale_color_cmocean(name='thermal')

# Swarmcount
p.density <- ggplot(cast.df, aes(x=lon, y=lat, col=temp_mean)) +
  # coastline
  geom_polygon(data=landShp, aes(x=long, y=lat, group=group), 
               fill='grey', col='#636363', size=.6) +
  # bubbles
  geom_point(shape=1, size=make_normal(cast.df$Sv_mean_mean, minsize, maxsize)) +
  coord_map(xlim=xlim,  ylim=ylim) +
  theme_bw() + grids(linetype = "dashed") +
  #scale_color_cmocean(name='thermal') +
  facet_wrap(~survey) +
  #theme(legend.position = 'none') +
  ggtitle("Mean swarm density") +
  labs(x="Longitude", y="Latitude") +
  scale_color_cmocean(name='thermal')

width=8
height=12
ggsave("./figures/bubbleplots/bubble_depth.png", p.depth, width=width, height=height)
ggsave("./figures/bubbleplots/bubble_desnity.png", p.density, width=width, height=height)
ggsave("./figures/bubbleplots/bubble_coverage.png", p.coverage, width=width, height=height)


# ------------- Correaltions ------------- #
# ---- ALL ---- #
cor.dat <- cast.df[,c('survey',
                      'temp_mean', 'salt_mean', 
                      'area_corrected_area_sum', 'Sv_mean_mean', 'swarm_mean_depth',
                      'lon', 'lat')]
names(cor.dat) <- c('Survey',
                    'Temperature', 'Salinity',
                    'Agg. Coverage', 'Mean Agg. Density', 'Mean Agg. Depth',
                    'Longitude', 'Latitude')

p.ls <- list()
p.store <- list()
i <- 1
for (survey in unique(cor.dat$Survey)){
  p <- ggcorrplot(cor(cor.dat[cor.dat$Survey == survey, 2:ncol(cor.dat)]),
                  colors=c('#0b3487','#ffffff','#87130b'),
                  lab=T, type='upper', lab_size=3) +
    theme(legend.position = 'none',
          axis.text.x = element_text(size=9),
          axis.text.y = element_text(size=9)) +
    ggtitle(paste0("(",chartr("123456789", "ABCDEFGHI", i),") ", survey))
  p.ls[[i]] <- p
  # atore p values
  p.store[[survey]] <- cor_pmat(cor.dat[cor.dat$Survey == survey, 2:ncol(cor.dat)])
  i <- i + 1
}
p.corr <- ggarrange(plotlist=p.ls)
ggsave("./figures/corr/intra_corr_simple.png", p.corr, width=14, height=14, bg='white')

# # ---- T/S lon/lat ---- #
# cor.dat <- cast.df[,c('survey',
#                       'temp_mean', 'salt_mean', 
#                       'lon', 'lat')]
# names(cor.dat) <- c('Survey',
#                     'Temperature', 'Salinity', 
#                     'Longitude', 'Latitude')
# # cor.dat <- cast.df[,c('survey', 'temp_mean', 'salt_mean', 'area_Sv_mean_sum')]
# # names(cor.dat) <- c('Survey', 'Temperature', 'Salinity', 'Biomass')
# p.ls <- list()
# i <- 1
# for (survey in unique(cor.dat$Survey)){
#   p <- ggcorrplot(cor(cor.dat[cor.dat$Survey == survey, 2:ncol(cor.dat)]),
#                   colors=c('#0b3487','#ffffff','#87130b'),
#                   lab=T, type='upper', lab_size=3) +
#     theme(legend.position = 'none',
#           axis.text.x = element_text(size=9),
#           axis.text.y = element_text(size=9)) +
#     ggtitle(survey)
#   p.ls[[i]] <- p
#   #ggsave(paste0('./output/correlations/corr_',survey,'.png'), p)
#   i <- i + 1
# }
# p.corr <- ggarrange(plotlist=p.ls)
# ggsave("./figures/corr/intra_corr_simple_TSlonlat.png", p.corr, width=14, height=14)
# 
# # ---- Agg TS ---- #
# cor.dat <- cast.df[,c('survey',
#                       'temp_mean', 'salt_mean', 
#                       'area_corrected_area_sum', 'Sv_mean_mean', 'swarm_mean_depth')]
# names(cor.dat) <- c('Survey',
#                     'Temperature', 'Salinity', 
#                     'Swarm Coverage', 'Mean Swarm Desnty', 'Mean Swarm Depth')
# # cor.dat <- cast.df[,c('survey', 'temp_mean', 'salt_mean', 'area_Sv_mean_sum')]
# # names(cor.dat) <- c('Survey', 'Temperature', 'Salinity', 'Biomass')
# p.ls <- list()
# i <- 1
# for (survey in unique(cor.dat$Survey)){
#   p <- ggcorrplot(cor(cor.dat[cor.dat$Survey == survey, 2:ncol(cor.dat)]),
#                   colors=c('#0b3487','#ffffff','#87130b'),
#                   lab=T, type='upper', lab_size=3) +
#     theme(legend.position = 'none',
#           axis.text.x = element_text(size=9),
#           axis.text.y = element_text(size=9)) +
#     ggtitle(survey)
#   p.ls[[i]] <- p
#   #ggsave(paste0('./output/correlations/corr_',survey,'.png'), p)
#   i <- i + 1
# }
# p.corr <- ggarrange(plotlist=p.ls)
# ggsave("./figures/corr/intra_corr_simple_TSagg.png", p.corr, width=14, height=14)
# 
# # ---- Agg T ---- #
# cor.dat <- cast.df[,c('survey',
#                       'temp_mean', 
#                       'area_corrected_area_sum', 'Sv_mean_mean', 'swarm_mean_depth')]
# names(cor.dat) <- c('Survey',
#                     'Temperature',
#                     'Swarm Coverage', 'Mean Swarm Desnty', 'Mean Swarm Depth')
# # cor.dat <- cast.df[,c('survey', 'temp_mean', 'salt_mean', 'area_Sv_mean_sum')]
# # names(cor.dat) <- c('Survey', 'Temperature', 'Salinity', 'Biomass')
# p.ls <- list()
# i <- 1
# for (survey in unique(cor.dat$Survey)){
#   p <- ggcorrplot(cor(cor.dat[cor.dat$Survey == survey, 2:ncol(cor.dat)]),
#                   colors=c('#0b3487','#ffffff','#87130b'),
#                   lab=T, type='upper', lab_size=3) +
#     theme(legend.position = 'none',
#           axis.text.x = element_text(size=9),
#           axis.text.y = element_text(size=9)) +
#     ggtitle(survey)
#   p.ls[[i]] <- p
#   #ggsave(paste0('./output/correlations/corr_',survey,'.png'), p)
#   i <- i + 1
# }
# p.corr <- ggarrange(plotlist=p.ls)
# ggsave("./figures/corr/intra_corr_simple_Tagg.png", p.corr, width=14, height=14)


# ------------- SST ------------- #
# Load SST data
surveys <- unique(cast.df$survey)
df.sst <- list()
for (i in 1:length(surveys)){
  dir <- paste0("./data/IMOS/SST_L3S/",surveys[i],"/")
  fn <- list.files(dir, pattern="*.nc")[1]
  nc_data <- nc_open(paste0(dir, fn))
  lon <- ncvar_get(nc_data, "lon")
  lat <- ncvar_get(nc_data, "lat", verbose = F)
  time <- ncvar_get(nc_data, "time")
  sst <- ncvar_get(nc_data, "sea_surface_temperature", 
                   start=c(1,1,1), count=c(-1,-1,10)) - 273.15 # convert from kelvin
  # flatten sst
  sst <- as.vector(apply(sst, 1:2, function(x) mean(x, na.rm=T)))

  
  
  # # Find the day with  best coverage (least NAs)
  # # Shape is [lon, lat,  t]
  # t.idx <- which.min(sapply(1:length(time), function(t) sum(is.na(sst[,,t]))))
  # # reduce SST
  # sst <- sst[,,t.idx]
  # time.frame <- as.POSIXct(time[t.idx], origin="1981-01-01",  tz="UTC")
  
  # Make into ggplot dataframe
  df.sst[[i]] <- data.frame(survey=surveys[i],
                       lon=rep(lon, length(lat)),
                       lat=rep(lat, each=length(lon)),
                       sst=as.vector(sst))

}
df.sst <- do.call(rbind, df.sst)

ggplot(df.sst, aes(x=lon, y=lat, fill=sst)) +
  geom_tile()  +
  scale_fill_cmocean(name='thermal')  +
  facet_wrap(~survey)

xlim <- c(150.04, 150.275)
ylim <- c(-36.47, -36.15)
fn <- 'bubble_sst_finescale'
# xlim <- c(150.04, 150.5)
# ylim <- c(-36.67, -35.85)
# fn <- 'bubble_sst_mesoscale'
maxsize <- 9
minsize <- 2
mintemp <- 14
maxtmep <- 20

# output raster so we can use in QGIS
df.sst.EAC <- df.sst[df.sst$survey == "2015_S1", c('lon','lat','sst')]
cast.df$survey_lab <- paste0("(",chartr("123456789", "ABCDEFGHI", as.numeric(as.factor(cast.df$survey))),") ",
                             cast.df$survey)
df.sst$survey_lab <- paste0("(",chartr("123456789", "ABCDEFGHI", as.numeric(as.factor(df.sst$survey))),") ",
                            df.sst$survey)
# make df.sst smaller to speed plotting up
df.sst <- df.sst[df.sst$lat > ylim[1]-0.1 & 
                 df.sst$lat < ylim[2]+0.1 &
                 df.sst$lon > xlim[1]-0.1 & 
                 df.sst$lon < xlim[2]+0.1,]

# Filter to only show 2015_S1, 2018_S1 and 2018_S2
survey.filter <- c("2015_S1", "2018_S1", "2018_S2")
cast.df.f <- cast.df[cast.df$survey %in% survey.filter,]
df.sst.f <- df.sst[df.sst$survey %in% survey.filter,]

# Swarm Coverage
p.coverage <- ggplot(cast.df.f, aes(x=lon, y=lat)) +
  geom_tile(mapping=aes(x=lon, y=lat, fill=sst),
            data=df.sst.f,
            inherit.aes = F)  +
  scale_fill_cmocean(name='thermal') +
  guides(fill = "none") +
  # coastline
  geom_polygon(data=landShp, aes(x=long, y=lat, group=group), 
               fill='grey', col='#636363', size=.6) +
  # Bubble
  geom_point(shape=1,
             mapping=aes(size=make_normal(area_corrected_area_sum, minsize, maxsize))) +
  scale_size(name = "Coverage (%)",
             breaks = fivenum(make_normal(cast.df.f$area_corrected_area_sum, minsize, maxsize)),
             labels = round(fivenum(cast.df.f$area_corrected_area_sum)*100, 2)) +
  # Map stuff
  coord_map(xlim=xlim,  ylim=ylim) +
  theme_bw() + grids(linetype = "dashed") +
  facet_wrap(~survey) +
  ggtitle("(A) Aggregation coverage") +
  labs(x=NULL, y=NULL)

# Swarm Density
p.density <- ggplot(cast.df.f, aes(x=lon, y=lat)) +
  geom_tile(mapping=aes(x=lon, y=lat, fill=sst),
            data=df.sst.f,
            inherit.aes = F)  +
  scale_fill_cmocean(name='thermal') +
  guides(fill = "none") +
  # coastline
  geom_polygon(data=landShp, aes(x=long, y=lat, group=group), 
               fill='grey', col='#636363', size=.6) +
  # Bubble
  geom_point(shape=1,
             mapping=aes(size=make_normal(Sv_mean_mean, minsize, maxsize))) +
  scale_size(name = bquote(""*S[v]*" mean"),
             breaks = fivenum(make_normal(cast.df.f$Sv_mean_mean, minsize, maxsize)),
             labels = round(fivenum(cast.df.f$Sv_mean_mean),2)) +
  # Map stuff
  coord_map(xlim=xlim,  ylim=ylim) +
  theme_bw() + grids(linetype = "dashed") +
  #scale_color_cmocean(name='thermal') +
  facet_wrap(~survey) +
  ggtitle(bquote("(B) Mean aggregation density")) +
  labs(x=NULL, y=NULL)

# Mean depth
p.depth <- ggplot(cast.df.f, aes(x=lon, y=lat)) +
  geom_tile(mapping=aes(x=lon, y=lat, fill=sst),
            data=df.sst.f,
            inherit.aes = F)  +
  scale_fill_cmocean(name='thermal') +
  guides(fill = "none") +
  # coastline
  geom_polygon(data=landShp, aes(x=long, y=lat, group=group), 
               fill='grey', col='#636363', size=.6) +
  # Bubble
  geom_point(shape=1,
             mapping=aes(size=make_normal(swarm_mean_depth, minsize, maxsize))) +
  scale_size(name = "Depth (m)",
             breaks = fivenum(make_normal(cast.df.f$swarm_mean_depth, minsize, maxsize)),
             labels = round(fivenum(cast.df.f$swarm_mean_depth),0)) +
  coord_map(xlim=xlim,  ylim=ylim) +
  theme_bw() + grids(linetype = "dashed") +
  #scale_color_cmocean(name='thermal') +
  facet_wrap(~survey) +
  ggtitle("(C) Mean aggregation depth") +
  labs(x=NULL, y=NULL)

p.out <- ggarrange(p.coverage, p.density, p.depth, ncol=1)
#legendpos <- theme(legend.position = "right")
# width = 6
# height = 9
ggsave(paste0("./figures/bubbleplots/bubbleOut.png"), p.out, width=8, height=12, bg='white')


# ggsave(paste0("./figures/bubbleplots/",fn,"_SwarmDepth.png"), p.depth + legendpos, width=width, height=height)
# ggsave(paste0("./figures/bubbleplots/",fn,"_SwarmDesnity.png"), p.density + legendpos, width=width, height=height)
# ggsave(paste0("./figures/bubbleplots/",fn,"_SwarmCoverage.png"), p.coverage + legendpos, width=width, height=height)


# plot legebd
library(cowplot)
legend.get <- cowplot::get_legend(p.depth + 
                                    theme(legend.key.height = unit(3, "cm"),
                                          legend.position = "right") +
                                    guides(size = "none", fill=guide_colourbar(order = 1)) +
                                    labs(fill='SST (°C)'))
grid.newpage()
grid.draw(legend.get)
ggsave(paste0("./figures/bubbleplots/colorbar.png"), legend.get, width=2)







# ------- INterannual correlations -----

# ----------- Agg-T ----------- #
cor.dat <- cast.df[,c('survey',
                      'temp_mean', 
                      'area_corrected_area_sum', 'Sv_mean_mean', 'swarm_mean_depth')]
names(cor.dat) <- c('Survey',
                    'Temperature',
                    'Swarm Coverage', 'Mean Swarm Desnty', 'Mean Swarm Depth')
# Correlation plot
cor.dat$Survey <- NULL
p.corr.inter <- ggcorrplot(cor(cor.dat),
                           colors=c('#0b3487','#ffffff','#87130b'),
                           lab=T, type='upper') +
  theme(legend.position = 'none')
ggsave("./figures/corr/inter_corr_simple_AggT.png", p.corr.inter)

# ---- T/S lon/lat ---- #
cor.dat <- cast.df[,c('survey',
                      'temp_mean', 'salt_mean', 
                      'lon', 'lat')]
names(cor.dat) <- c('Survey',
                    'Temperature', 'Salinity', 
                    'Longitude', 'Latitude')
cor.dat$Survey <- NULL
# Correlation plot
p.corr.inter <- ggcorrplot(cor(cor.dat),
                           colors=c('#0b3487','#ffffff','#87130b'),
                           lab=T, type='upper') +
  theme(legend.position = 'none')
ggsave("./figures/corr/inter_corr_simple_AggTS-lonlat.png", p.corr.inter)

# ---- Agg TS ---- #
cor.dat <- cast.df[,c('survey',
                      'temp_mean', 'salt_mean', 
                      'area_corrected_area_sum', 'Sv_mean_mean', 'swarm_mean_depth')]
names(cor.dat) <- c('Survey',
                    'Temperature', 'Salinity', 
                    'Swarm Coverage', 'Mean Swarm Desnty', 'Mean Swarm Depth')
cor.dat$Survey <- NULL
# Correlation plot
p.corr.inter <- ggcorrplot(cor(cor.dat),
                           colors=c('#0b3487','#ffffff','#87130b'),
                           lab=T, type='upper') +
  theme(legend.position = 'none')
ggsave("./figures/corr/inter_corr_simple_AggTS.png", p.corr.inter)






# AND USING MEANS
mean.casts <- aggregate(cast.df, by=list(cast.df$survey), mean)
mean.casts[,2:4] <- NULL
names(mean.casts)[1] <- "survey"

# ----------- Agg-T ----------- #
cor.dat <- mean.casts[,c('survey',
                      'temp_mean', 
                      'area_corrected_area_sum', 'Sv_mean_mean', 'swarm_mean_depth')]
names(cor.dat) <- c('Survey',
                    'Temperature',
                    'Swarm Coverage', 'Mean Swarm Desnty', 'Mean Swarm Depth')
cor.dat$Survey <- NULL
# Correlation plot
p.corr.inter <- ggcorrplot(cor(cor.dat),
                           colors=c('#0b3487','#ffffff','#87130b'),
                           lab=T, type='upper') +
  theme(legend.position = 'none')
ggsave("./figures/corr/intermean_corr_simple_AggT.png", p.corr.inter)

# ---- T/S lon/lat ---- #
cor.dat <- mean.casts[,c('survey',
                      'temp_mean', 'salt_mean', 
                      'lon', 'lat')]
names(cor.dat) <- c('Survey',
                    'Temperature', 'Salinity', 
                    'Longitude', 'Latitude')
cor.dat$Survey <- NULL
# Correlation plot
p.corr.inter <- ggcorrplot(cor(cor.dat),
                           colors=c('#0b3487','#ffffff','#87130b'),
                           lab=T, type='upper') +
  theme(legend.position = 'none')
ggsave("./figures/corr/intermean_corr_simple_AggTS-lonlat.png", p.corr.inter)

# ---- Agg TS ---- #
cor.dat <- mean.casts[,c('survey',
                      'temp_mean', 'salt_mean', 
                      'area_corrected_area_sum', 'Sv_mean_mean', 'swarm_mean_depth')]
names(cor.dat) <- c('Survey',
                    'Temperature', 'Salinity', 
                    'Swarm Coverage', 'Mean Swarm Desnty', 'Mean Swarm Depth')
cor.dat$Survey <- NULL
# Correlation plot
p.corr.inter <- ggcorrplot(cor(cor.dat),
                           colors=c('#0b3487','#ffffff','#87130b'),
                           lab=T, type='upper') +
  theme(legend.position = 'none')
ggsave("./figures/corr/intermean_corr_simple_AggTS.png", p.corr.inter)




# FINAL INTERCOR
mean.casts <- aggregate(cast.df, by=list(cast.df$survey), mean)
mean.casts[,2:4] <- NULL
names(mean.casts)[1] <- "survey"

cor.dat <- mean.casts[,c('survey',
                      'temp_mean', 'salt_mean', 
                      'area_corrected_area_sum', 'Sv_mean_mean', 'swarm_mean_depth')]
names(cor.dat) <- c('Survey',
                    'Temperature', 'Salinity',
                    'Agg. Coverage', 'Mean Agg. Density', 'Mean Agg. Depth')
cor.dat$Survey <- NULL
# Correlation plot
p.corr.inter <- ggcorrplot(cor(cor.dat),
                           colors=c('#0b3487','#ffffff','#87130b'),
                           lab=T, type='upper', p.mat=cor_pmat(cor.d05at)) +
  theme(legend.position = 'none')
ggsave("./figures/corr/intermean_corr.png", p.corr.inter, width=7, height=7)


