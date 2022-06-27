setwd("/Users/lachlanphillips/Development/PhD/paper-repos/BMP-Variability")
source('R/CTD.R')
source('R/acoustics_functions.R')
source('R/modified_functions.R')
library(nlme)
library(cmocean)
library(ggplot2)
library(ggcorrplot)
library(reshape2)
library(cluster)
library(factoextra)
library(nnet)
library(scatterplot3d)
library(dbscan) # clustering for anisotropic data

# ggbiplot
# library(devtools)
# install_github("vqv/ggbiplot")
library(ggbiplot)
library(psych) # pair.panels function

# For editing ggplot objects
remove_geom <- function(ggplot2_object, geom_type) {
  # Delete layers that match the requested type.
  layers <- lapply(ggplot2_object$layers, function(x) {
    if (class(x$geom)[1] == geom_type) {
      NULL
    } else {
      x
    }
  })
  # Delete the unwanted layers.
  layers <- layers[!sapply(layers, is.null)]
  ggplot2_object$layers <- layers
  ggplot2_object
}


#'------------------------------#
#' __1. Load CTD and Agg data__ #
#'------------------------------#
# CTD data
CTD <- cast.reader()
# Aggregation data
agg <- load.agg()
# Extract agg data for each CTD cast
cast.df <- CTD.agg(CTD, agg)

# Add lon and lat
cast.df$lon <- sapply(cast.df$cast_id, function(id) CTD$meta$start_longitude[CTD$meta$file_name == id])
cast.df$lat <- sapply(cast.df$cast_id, function(id) CTD$meta$start_latitude[CTD$meta$file_name == id])


#'----------------------------------#
#' __2. Make Correlation datasets__ #
#'----------------------------------#
cor.dat <- cast.df[,c('temp_mean', 'salt_mean',
                      'area_corrected_area_sum', 'Sv_mean_mean', 'swarm_mean_depth')]
names(cor.dat) <- c('Temperature', 'Salinity',
                    'Agg. Coverage', 'Mean Agg. Density', 'Mean Agg. Depth')

# Make density linear (Correct?)
#cor.dat[,'Mean Agg. Density'] <- Sv_mean.linear(cor.dat[,'Mean Agg. Density'])


#'-------------------------------------------#
#' __3. Make Interannual Correlation Plots__ #
#'-------------------------------------------#
ggcorrplot(cor(cor.dat),
           colors=c('#0b3487','#ffffff','#87130b'),
           lab=T, type='upper', p.mat=cor_pmat(cor.dat)) +
  theme(legend.position = 'none')

pairs.panels(cor.dat,
             gap = 0,
             bg = as.factor(cast.df$survey),
             pch=21, lm=F, ci=T,
             ellipses=F, hist.col='grey')

# compare with aggregated years
cor.dat.years <- aggregate(cor.dat, by=list(cast.df$survey), mean)
pairs.panels(cor.dat.years[,-1],
             gap = 0,
             bg = as.factor(cor.dat.years$Group.1),
             pch=21, lm=T, ci=T,
             ellipses=F, hist.col='grey')


#'---------------------#
#' __4. PCA analysis__ #
#'---------------------#
# Build PCA model
pc <- prcomp(cor.dat, center = TRUE, scale. = TRUE)
ggbiplot(pc, ellipse=TRUE, obs.scale = 1, var.scale = 1, var.axes=TRUE,
         groups=cast.df$survey)
print(pc)
summary(pc)

# Check scatter of PCA (pros removal of multicollinearity)
pairs.panels(pc$x,
             gap = 0,
             bg = as.factor(cor.dat.years$Group.1),
             pch=21, lm=T, ci=T,
             ellipses=F, hist.col='grey')

#' __Part 2 - PCA axis Boxplot__ #
pca.box <- cbind(pc$x[,1], cast.df$survey)
pca.box <- as.data.frame(pca.box)
pca.box$V1 <- as.numeric(pca.box$V1)
pca.box$V2 <- as.factor(pca.box$V2)
names(pca.box) <- c('PCA1', 'Survey')

ggplot(pca.box, aes(x=Survey, y=PCA1)) + 
  geom_boxplot(fill='grey') +
  theme_bw() +
  labs(x=NULL, y='PCA axis 1')

#' __Part 3 - PCA spatial__ #
# Add PCA axis 1 to cast.df
cast.df$PCA1 <- as.numeric(pc$x[,1])

# Load grid 
stations <- read.csv('./data/transects/ghostStations.csv')
grid <- read.csv('./data/transects/oceanGrid.csv')
# Create a set of shape files
grid.polys <- grid.polygons(grid, ids=stations$id)

# get cell numbers
cast.df$cell <- as.numeric(factor(cast.df$station_id, levels=stations$id))
# calc ffor each survey
cell.PCA <- aggregate(cbind(PCA1, lon, lat) ~ survey + cell, data=cast.df, FUN=mean)
# make coords
cell.PCA$col <- cell.PCA$cell %% 5
cell.PCA$col[cell.PCA$col == 0] <- 5
cell.PCA$row <- ceiling(cell.PCA$cell / 5)

ggplot(cell.PCA, aes(x=col, y=row, fill=PCA1)) +
  geom_tile() +
  facet_wrap(~survey) +
  scale_fill_cmocean('thermal', direction=-1)

#'---------------------------------------#
#' __5. Multinominal regression on PCA__ #
#'---------------------------------------#
#'# Divide data into training and testing
ind <- sample(2, nrow(cor.dat),
              replace = TRUE,
              prob = c(0.8, 0.2))
cor.dat.train <- cor.dat[ind==1,]
cor.dat.test <- cor.dat[ind==2,]

# Build PCA model
pc.mn <- prcomp(cor.dat.train, center = TRUE, scale. = TRUE)

# Check preidctions on testing set
trn <- predict(pc.mn, cor.dat.train)
trn <- data.frame(trn, survey=cast.df$survey[ind==1])
tst <- predict(pc.mn, cor.dat.test)
tst <- data.frame(tst, survey=cast.df$survey[ind==2])

# Multinominal regression
mn.reg <- multinom(survey ~ PC1 + PC2 + PC3 + PC4 + PC5, data=trn)
summary(mn.reg)

# Confusion Matrices
# Training
p.trn <- predict(mn.reg, trn)
tab <- table(p.trn, trn$survey)
tab
print(paste0('Model accuracy: ',round((sum(diag(tab))/sum(tab))*100,2),'%'))
# Testing
p.tst <- predict(mn.reg, tst)
tab <- table(p.tst, tst$survey)
tab
print(paste0('Model accuracy: ',round((sum(diag(tab))/sum(tab))*100,2),'%'))


#'----------------------------------#
#' __6. K Means Clustering on PCA__ #
#'----------------------------------#
#' Note!! Earlier plots using fviz_cluster were actually doing a PCA on a PCA!
# Apply K means cluster to PCA
pc.transform <- as.data.frame(pc$x)
# detrmine optimal clusters
fviz_nbclust(pc.transform, kmeans, method = 'silhouette')
#fviz_nbclust(pc.transform, kmeans, method = 'wss')
# gap_stat <- clusGap(pc.transform, FUN = kmeans, nstart = 25, K.max = 15, B = 50)
# fviz_gap_stat(gap_stat)
# ccluster data
kmeans_dat = kmeans(pc.transform, centers = 2, nstart = 50)
# Plot clusters
p <- fviz_cluster(kmeans_dat, data = pc.transform, ellipse=T, ellipse.alpha=0)
p + geom_point(mapping=aes(color = factor(cast.df$survey)))

# Extract data and make a better plot
plot.data <- fviz_cluster.mod(kmeans_dat, data = cor.dat, 
                  ellipse=T, ellipse.alpha=0, shape=1, 
                  show.clust.cent=F)
plot.data$survey <- cast.df$survey

ggplot(plot.data) +
  geom_point(aes(x=x, y=y, color=survey, shape=cluster), size=2) + # 
  stat_ellipse(aes(x=x, y=y, group=cluster), type = "norm", 
               alpha=1, linetype='dashed', color='darkgrey') +
  labs()



p <- remove_geom(p, "GeomText")
p + ggtitle(NULL) + labs(col='Cluster',  shape='Cluster', fill='Cluster') +
  theme_bw()



#'---------------------#
#' __7. K means only__ #
#'---------------------#
# detrmine optimal clusters 
fviz_nbclust(cor.dat, kmeans, method = 'silhouette')
# ccluster data
kmeans_dat = kmeans(cor.dat, centers = 2, nstart = 50)
# Plot clusters
p <- fviz_cluster(kmeans_dat, data = cor.dat, ellipse=T) # This is wrong... 
p + geom_point(mapping=aes(colour = cast.df$survey))


#'------------------------------------#
#' __8. Mixed effects model for T/S__ #
#'------------------------------------#
head(cast.df)

fit.TS <- lm(temp_mean ~ salt_mean,data = cast.df, method='ML')
fit.TS.sur <- lme(temp_mean ~ salt_mean, random = ~ salt_mean | survey, data = cast.df, method='ML')
summary(fit.TS)
AIC(fit.TS)
summary(fit.TS.sur)
plot(fit.TS.sur)
#'-----------------------------#
#' __8. T/S plot time scales__ #
#'-----------------------------#
ggplot(cast.df, aes(x=salt_mean, y=temp_mean, fill=survey, color=survey, group=survey)) +
  geom_point(shape=3, alpha=.8) +
  geom_smooth(method='lm', alpha=.3) +
  geom_smooth(method='lm', mapping=aes(x=salt_mean, y=temp_mean), inherit.aes = F, color='black') +
  theme_bw() +
  labs(x='Mean Salinity (PSU)', y='Mean Temperature (°C)', color='Survey', fill='Survey')


#'----------------------------------------------------#
#' __9. PCA analysis on acoustics and T/S separated__ #
#'----------------------------------------------------#
# Build PCA model
cor.dat.aggs <- cor.dat[,-(1:2)]
pc.agg <- prcomp(cor.dat.aggs, center = TRUE, scale. = TRUE)
ggbiplot(pc.agg, ellipse=TRUE, obs.scale = 1, var.scale = 1, var.axes=TRUE,
         groups=cast.df$survey)
print(pc.agg)
summary(pc.agg)

# Build PCA model
cor.dat.TS <- cor.dat[,1:2]
pc.TS <- prcomp(cor.dat.TS, center = TRUE, scale. = TRUE)
ggbiplot(pc.TS, ellipse=TRUE, obs.scale = 1, var.scale = 1, var.axes=TRUE,
         groups=cast.df$survey)
print(pc.TS)
summary(pc.TS)

# Make boxplots
survey = cast.df$survey
cluster <- survey
cluster[cluster %in% c('2015_S1','2016_S1', '2016_S2')] <- 'A'
cluster[cluster != 'A']  <- 'B'
ggplot(melt(cbind(cor.dat.aggs, survey, cluster)), aes(x=survey, y=value, fill=cluster)) +
  geom_boxplot() + facet_wrap(~variable, scales='free') +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x=NULL)

cluster <- survey
cluster[cluster %in% c('2015_S1','2016_S1', '2016_S2')] <- 'A'
cluster[cluster %in% c('2017_S1','2017_S2', '2018_S1', '2018_S2')] <- 'B'
cluster[cluster %in% c('2019_S1', '2019_S2')] <- 'C'
ggplot(melt(cbind(cor.dat.aggs, survey, cluster)), aes(x=survey, y=value, fill=cluster)) +
  geom_boxplot() + facet_wrap(~variable, scales='free') +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x=NULL)

# YEAR
# Make boxplots
survey = cast.df$survey
survey <- substr(survey, 1, 4)
cluster <- survey
cluster[cluster %in% c('2015','2016')] <- 'A'
cluster[cluster != 'A']  <- 'B'
ggplot(melt(cbind(cor.dat.aggs, survey, cluster)), aes(x=survey, y=value, fill=cluster)) +
  geom_boxplot() + facet_wrap(~variable, scales='free') +
  geom_point()  +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x=NULL)

cluster <- cast.df$survey
cluster[cluster %in% c('2015_S1','2016_S1', '2016_S2')] <- 'A'
cluster[cluster != 'A']  <- 'B'
melt.casts <- melt(cbind(cast.df[,-c(2,3,4)],cluster), c('survey', 'cluster'))

# Plot the lot
ggplot(melt.casts, aes(x=survey, y=value, fill=cluster)) +
  geom_boxplot() + facet_wrap(~variable, scales='free') +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x=NULL)

