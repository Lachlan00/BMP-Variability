# make map panel
library(rgdal)
library(stringr)
library(ggplot2)
library(ggpubr)
library(maps)
library(sp)
library(maptools)
source('R/eudyptula.R')

# load df.sst.EAC from bubble plot script

#####################################################
# EAC panel 2 #
#####################################################
# --- Data Prep --- #
# Load in transect data
transect.lines <- read.csv('data/transects/transect_lines.csv')
#  stations
stations <- read.csv('data/transects/transect_df.csv')
moorings <- read.csv('data/surveys/moorings.csv')
# Reshape for ggplot
transect.lines1 <- transect.lines[,1:3]
transect.lines2 <- transect.lines[,c(1,4,5)]
names(transect.lines2)  <- names(transect.lines1)
transect.lines <- rbind(transect.lines1, transect.lines2)
# Load in coastline shapefile
fn = './assets/shp/NSW-coastline/edited/polygon/NSW-coastline_WGS84'
landShp <- readOGR(dsn=path.expand(str_sub(fn,1,-(nchar(basename(fn)))-1)), 
                   layer=basename(fn))

# --- Plot Map --- #
xlim <- c(150.04, 150.33)
ylim <- c(-36.47, -36.15)
summary(df.sst.EAC)
station_labels <- stations[stations$station == 1,]
station_labels$id <- paste0("T",station_labels$transect)
# plot
p.panel <- ggplot() +
  # SST
  geom_tile(mapping=aes(x=lon, y=lat, fill=sst),
            data=df.sst.EAC,
            inherit.aes = F)  +
  scale_fill_cmocean(name='thermal', limits=c(12, 28))  +
  # transects
  geom_line(data=transect.lines, mapping=aes(x=lon1, y=lat1, group=transect), 
            size=.5, linetype='dashed') +
  # stations
  geom_point(data=stations, mapping=aes(x=lon, y=lat),
             color='black', fill='grey', size=4, shape=21) +
  geom_point(data=stations, mapping=aes(x=lon, y=lat),
             color='black', size=1.5, shape=1) +
  # moorings
  geom_point(data=moorings, mapping=aes(x=lon, y=lat),
             color='black', fill="green", size=4, shape=23) +
  # coastline
  geom_polygon(data=landShp, aes(x=long, y=lat, group=group), 
               fill='grey', col='#636363', size=.4) +
  # station labels
  geom_text(data=station_labels, mapping=aes(x=lon, y=lat, label=id), size=3,
            nudge_x=-.005, nudge_y=.008, fontface='bold') +
  # Mooring labels
  geom_label(data=moorings, mapping=aes(x=lon, y=lat, label=id), size=2.7, 
             fill="green", alpha=0.4, nudge_x=0, nudge_y=-.011) +
  # settings
  coord_map(xlim=xlim,  ylim=ylim) +
  theme_bw() + grids(linetype = "dashed") +
  theme(legend.position = "none") +
  labs(x=NULL, y=NULL, alpha=NULL, size=NULL, fill=NULL) +
  theme(text = element_text(size=14)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank())

# Save files
ggsave(p.panel, filename='figures/Map/panel2.png',
       bg="transparent",dpi=300, width=2000/300, height=2000/300) 
  

