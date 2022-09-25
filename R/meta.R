# Survey dimensions
library(geosphere)
df <- read.csv('./data/transects/transect_df.csv')
width <- distHaversine(c(min(df$lon), mean(df$lat)),  c(max(df$lon), mean(df$lat)))/1000
height <- distHaversine(c(mean(df$lon), min(df$lat)),  c(mean(df$lon), max(df$lat)))/1000

width * height
