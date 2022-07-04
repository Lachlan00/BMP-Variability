setwd("/Users/lachlanphillips/Development/PhD/pape r-repos/BMP-Variability")
source('R/CTD.R')
source('R/acoustics_functions.R')
source('R/modified_functions.R')
library(nlme)
library(cmocean)
library(ggplot2)
library(ggcorrplot)
library(ggthemes)
library(ggpubr)
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

#'----------------------------------------------------#
#' __5. PCA analysis on acoustics and T/S separated__ #
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

#'---------------------------------------#
#' __6. Multinominal regression on PCA__ #
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
#' __7. K Means Clustering on PCA__ #
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
# Extract PCA data and make a better plot (this functions performs PCA)
plot.data <- fviz_cluster.mod(kmeans_dat, data = cor.dat, 
                  ellipse=T)
plot.data$survey <- cast.df$survey
plot.data$cluster <- gsub('1','B',plot.data$cluster)
plot.data$cluster <- gsub('2','A',plot.data$cluster)

# Add clusters to cast.df
cast.df$cluster <- plot.data$cluster

# Get axis from biplot
df.v <- ggbiplot.mod(pc, ellipse=TRUE, obs.scale = 1, var.scale = 1, var.axes=TRUE,
          groups=cast.df$survey)
# Adjust arrow size bigger
df.v[,1:2] <- df.v[,1:2]*1.0
# Adjust angles
df.v[1,'angle'] <- df.v[1,'angle'] + 10 # Temp
df.v[2,'angle'] <- df.v[2,'angle'] + 6 # Salt
df.v[3,'angle'] <- df.v[3,'angle'] - 4 # Coverage
df.v[4,'angle'] <- df.v[4,'angle'] - 7 # Density
df.v[5,'angle'] <- df.v[5,'angle'] + 0 # Depth

p.pca <- ggplot(plot.data) +
  geom_point(aes(x=x, y=y, color=survey, shape=cluster), size=1.5) +
  #scale_color_cmocean(name='phase',  discrete=T) +
  stat_ellipse(aes(x=x, y=y, group=cluster), type = "norm", 
               alpha=.5, linetype='dashed', size=.6, segments=51) +
  labs(x='Principal Component 1', y='Principal Component 2', color='Survey', shape='Cluster') +
  theme_bw() +
  # axis
  geom_segment(data = df.v,
               aes(x = 0, y = 0, xend = xvar, yend = yvar),
               arrow = arrow(length = unit(0.6, 'picas')), 
               color = muted('red')) +
  geom_text(data = df.v, 
            aes(label = varname, x = xvar, y = yvar, 
                angle = angle, hjust = hjust), 
            color = 'darkred', size = 3)

ggsave('./figures/PCA/pca-cluster.png', p.pca, width=9, height=6)


#'-------------------------------------#
#' __8. K Means Clustering on PCA TS__ #
#'-------------------------------------#
#' Note!! Earlier plots using fviz_cluster were actually doing a PCA on a PCA!
# Apply K means cluster to PCA
pc.TS.transform <- as.data.frame(pc.TS$x)
# detrmine optimal clusters
fviz_nbclust(pc.TS.transform, kmeans, method = 'silhouette')
#fviz_nbclust(pc.transform, kmeans, method = 'wss')
# gap_stat <- clusGap(pc.transform, FUN = kmeans, nstart = 25, K.max = 15, B = 50)
# fviz_gap_stat(gap_stat)
# ccluster data
kmeans_dat.TS = kmeans(pc.TS.transform, centers = 2, nstart = 50)
# Extract PCA data and make a better plot (this functions performs PCA)
plot.data.TS <- fviz_cluster.mod(kmeans_dat.TS, data = cor.dat[,1:2], 
                              ellipse=T)
plot.data.TS$survey <- cast.df$survey
plot.data.TS$cluster <- gsub('1','A', plot.data.TS$cluster)
plot.data.TS$cluster <- gsub('2','B', plot.data.TS$cluster)

# Add clusters to cast.df
cast.df$cluster_TS <- plot.data.TS$cluster

# Get axis from biplot
df.v.TS <- ggbiplot.mod(pc.TS, ellipse=TRUE, obs.scale = 1, var.scale = 1, var.axes=TRUE,
                     groups=cast.df$survey)
# Adjust arrow size bigger
df.v.TS[,1:2] <- df.v.TS[,1:2]*1.0

p.pca.TS <- ggplot(plot.data.TS) +
  geom_point(aes(x=x, y=y, color=survey, shape=cluster), size=1.5) +
  #scale_color_cmocean(name='phase',  discrete=T) +
  stat_ellipse(aes(x=x, y=y, group=cluster), type = "norm", 
               alpha=.5, linetype='dashed', size=.6, segments=51) +
  labs(x='Principal Component 1', y='Principal Component 2', color='Survey', shape='Cluster') +
  theme_bw() +
  # axis
  geom_segment(data = df.v.TS,
               aes(x = 0, y = 0, xend = xvar, yend = yvar),
               arrow = arrow(length = unit(0.6, 'picas')), 
               color = muted('red')) +
  geom_text(data = df.v.TS, 
            aes(label = varname, x = xvar, y = yvar, 
                angle = angle, hjust = hjust), 
            color = 'darkred', size = 3)

ggsave('./figures/PCA/pca_TS-cluster.png', p.pca.TS, width=9, height=6)

#'-----------------------------#
#' __9. T/S plot__ #
#'-----------------------------#
p.TS <- ggplot(cast.df, aes(x=salt_mean, y=temp_mean, fill=survey, color=survey, group=survey)) +
  geom_point(shape=3, alpha=.7) +
  geom_smooth(method='lm', alpha=.3, se=F) +
  # geom_smooth(method='lm', mapping=aes(x=salt_mean, y=temp_mean),
  #             inherit.aes = F, color='black', se=F, linetype='dashed', alpha=.0) +
  geom_line(mapping=aes(x=salt_mean, y=temp_mean), 
            stat="smooth", method = "lm",
            linetype ="dashed",
            alpha = 0.9, size=1.1, 
            inherit.aes = F) +
  theme_bw() +
  labs(x='Mean Salinity (PSU)', y='Mean Temperature (°C)', color='Survey', fill='Survey')
ggsave('./figures/TS/TS-scatter.png', p.TS, width=8, height=5.5)


#'--------------------------------------------#
#' __10. Plot Temp against Agg for Clusters__ #
#'--------------------------------------------#
point_alpha=.43
# AGGs
p1 <- ggplot(cast.df, aes(x=temp_mean, y=swarm_mean_depth, color=cluster, group=cluster)) +
  geom_point(shape=3, alpha=point_alpha) +
  geom_smooth(method='lm', alpha=.3, se=F) +
  theme_bw() +
  scale_y_reverse() + 
  labs(x='Mean Temperature (°C)', y='Mean Agg. Depth (m)', color='Cluster') +
  ggtitle('(A) Mean Aggregation Depth')

p2 <- ggplot(cast.df, aes(x=temp_mean, y=area_corrected_area_sum, color=cluster, group=cluster)) +
  geom_point(shape=3, alpha=point_alpha) +
  geom_smooth(method='lm', alpha=.3, se=F) +
  theme_bw() +
  labs(x='Mean Temperature (°C)', y='Mean Agg. Coverage (%)', color='Cluster') +
  ggtitle('(B) Mean Aggregation Coverage')

p3 <- ggplot(cast.df, aes(x=temp_mean, y=Sv_mean_mean, color=cluster, group=cluster)) +
  geom_point(shape=3, alpha=point_alpha) +
  geom_smooth(method='lm', alpha=.3, se=F) +
  theme_bw() +
  labs(x='Mean Temperature (°C)', y='Mean Agg. Density (Sv mean)', color='Cluster') +
  ggtitle('(C) Mean Aggregation Density')

p4 <- ggplot(cast.df, aes(x=temp_mean, y=area_Sv_mean_sum, color=cluster, group=cluster)) +
  geom_point(shape=3, alpha=point_alpha) +
  geom_smooth(method='lm', alpha=.3, se=F) +
  theme_bw() +
  labs(x='Mean Temperature (°C)', y='Total Agg. Density (Sv mean)', color='Cluster') +
  ggtitle('(D) Total Biomass')

p.ls <- list(p1, p2, p3, p4)
p.ls <- lapply(p.ls, function(p) p + theme(legend.position = 'none'))
p.agg <- ggarrange(plotlist = p.ls)#, common.legend = T, legend = 'bottom')
#names(cast.df)
ggsave('./figures/TS/agg-cluster-scatter.png', p.agg, width=7, height=6)


#'-------------------------------------------#
#' __11. Plot Temp against Agg for Surveys__ #
#'-------------------------------------------#
point_alpha=.75
# AGGs
p1 <- ggplot(cast.df, aes(x=temp_mean, y=swarm_mean_depth)) +
  geom_point(shape=3, alpha=point_alpha, mapping=aes(color=survey)) +
  geom_smooth(method='lm', alpha=.3, se=F, color='black') +
  theme_bw() +
  scale_y_reverse() + 
  labs(x='Mean Temperature (°C)', y='Mean Agg. Depth (m)', color='Survey') +
  ggtitle('(A) Mean Aggregation Depth')

p2 <- ggplot(cast.df, aes(x=temp_mean, y=area_corrected_area_sum)) +
  geom_point(shape=3, alpha=point_alpha, mapping=aes(color=survey)) +
  geom_smooth(method='lm', alpha=.3, se=F, color='black') +
  theme_bw() +
  labs(x='Mean Temperature (°C)', y='Mean Agg. Coverage (%)', color='Survey') +
  ggtitle('(B) Mean Aggregation Coverage')

p3 <- ggplot(cast.df, aes(x=temp_mean, y=Sv_mean_mean)) +
  geom_point(shape=3, alpha=point_alpha, mapping=aes(color=survey)) +
  geom_smooth(method='lm', alpha=.3, se=F, color='black') +
  theme_bw() +
  labs(x='Mean Temperature (°C)', y='Mean Agg. Density (Sv mean)', color='Survey') +
  ggtitle('(C) Mean Aggregation Density')

p4 <- ggplot(cast.df, aes(x=temp_mean, y=area_Sv_mean_sum)) +
  geom_point(shape=3, alpha=point_alpha, mapping=aes(color=survey)) +
  geom_smooth(method='lm', alpha=.3, se=F, color='black') +
  theme_bw() +
  labs(x='Mean Temperature (°C)', y='Total Agg. Density (Sv mean)', color='Survey') +
  ggtitle('(D) Total Biomass')

p.ls <- list(p1, p2, p3, p4)
p.ls <- lapply(p.ls, function(p) p + theme(legend.position = 'none'))
p.agg.survey <- ggarrange(plotlist = p.ls)
#names(cast.df)
ggsave('./figures/TS/agg-survey-scatter.png', p.agg.survey, width=7, height=6)

#'---------------------------------------------------------#
#' __12. Plot Temp against Agg and all data for Clusters__ #
#'---------------------------------------------------------#
point_alpha=.43
se=T
# AGGs
p1 <- ggplot(cast.df, aes(x=temp_mean, y=swarm_mean_depth)) +
  geom_point(shape=3, alpha=point_alpha, mapping=aes(color=cluster, group=cluster)) +
  geom_smooth(method='glm', alpha=.3, se=se, mapping=aes( color=cluster, group=cluster)) +
  geom_smooth(method='glm', alpha=.3, se=se, color='black', linetype='dashed') +
  theme_bw() +
  scale_y_reverse() + 
  labs(x='Mean Temperature (°C)', y='Mean Agg. Depth (m)', color='Cluster') +
  ggtitle('(A) Mean Aggregation Depth')

p2 <- ggplot(cast.df, aes(x=temp_mean, y=area_corrected_area_sum)) +
  geom_point(shape=3, alpha=point_alpha, mapping=aes(color=cluster, group=cluster)) +
  geom_smooth(method='glm', alpha=.3, se=se, mapping=aes( color=cluster, group=cluster)) +
  geom_smooth(method='glm', alpha=.3, se=se, color='black', linetype='dashed') +
  theme_bw() +
  labs(x='Mean Temperature (°C)', y='Mean Agg. Coverage (%)', color='Cluster') +
  ggtitle('(B) Mean Aggregation Coverage')

p3 <- ggplot(cast.df, aes(x=temp_mean, y=Sv_mean_mean)) +
  geom_point(shape=3, alpha=point_alpha, mapping=aes(color=cluster, group=cluster)) +
  geom_smooth(method='glm', alpha=.3, se=se, mapping=aes( color=cluster, group=cluster)) +
  geom_smooth(method='glm', alpha=.3, se=se, color='black', linetype='dashed') +
  theme_bw() +
  labs(x='Mean Temperature (°C)', y='Mean Agg. Density (Sv mean)', color='Cluster') +
  ggtitle('(C) Mean Aggregation Density')

p4 <- ggplot(cast.df, aes(x=temp_mean, y=area_Sv_mean_sum)) +
  geom_point(shape=3, alpha=point_alpha, mapping=aes(color=cluster, group=cluster)) +
  geom_smooth(method='glm', alpha=.3, se=se, mapping=aes(color=cluster, group=cluster)) +
  geom_smooth(method='glm', alpha=.3, se=se, color='black', linetype='dashed') +
  theme_bw() +
  labs(x='Mean Temperature (°C)', y='Total Agg. Density (Sv mean)', color='Cluster') +
  ggtitle('(D) Total Biomass')

p.ls <- list(p1, p2, p3, p4)
p.ls <- lapply(p.ls, function(p) p + theme(legend.position = 'none'))
p.agg <- ggarrange(plotlist = p.ls)#, common.legend = T, legend = 'bottom')
#names(cast.df)
#ggsave('./figures/TS/agg-cluster-scatter.png', p.agg, width=7, height=6)


#'-----------------------------#
#' __13. Agg coverage vs PCA__ #
#'-----------------------------#
ggplot(cast.df, aes(x=temp_mean, y=Sv_mean_mean)) +
  geom_point(shape=3, alpha=0.5) +
  geom_smooth(method='lm', alpha=.5, se=F, color='black') +
  theme_bw() +
  labs(x='Temperature', y='Mean Agg. Density (Sv mean)', color='Survey') +
  facet_wrap(~survey, scales='free')

ggplot(cast.df, aes(x=temp_mean, y=swarm_mean_depth)) +
  geom_point(shape=3, alpha=0.5) +
  scale_y_reverse() + 
  geom_smooth(method='lm', alpha=.5, se=F, color='black') +
  theme_bw() +
  labs(x='Temperature', y='Mean Agg. Coverage (%)', color='Survey') +
  facet_wrap(~survey, scales='free')


#'----------------------------#
#' __14. Statistical Models__ #
#'----------------------------#
#'
#'### MOVE TOO RMARKDOWN
#'
#' saveRDS(cast.df, './data/processed/cast_df_stats.rds')
#' library(sjPlot) # for statistical plotting
#' 
#' #' __GLMS__ #
#' # Agg depth
#' hist(log1p(cast.df$swarm_mean_depth))
#' fit.agg.depth <- glm(log1p(swarm_mean_depth) ~ temp_mean, data=cast.df, family='gaussian')
#' summary(fit.agg.depth)
#' tab_model(fit.agg.depth)
#' 
#' head(cast.df)
#' fit.1 <- lme(temp_mean ~ salt_mean, random = ~ 1 | survey, data = cast.df, method='ML')
#' fit.2 <- lme(temp_mean ~ salt_mean, random = ~ salt_mean | survey, data = cast.df, method='ML')
#' summary(fit.1)
#' summary(fit.2)
#' plot(fit.1)
#' plot(fit.2)
#' 
#' # agg models
#' fit.3 <- lme(Sv_mean_mean ~ temp_mean, random = ~ temp_mean | survey, data = cor.dat.sur, method='ML')
#' summary(fit.3)
#' plot(fit.3)
#' 
#' fit.4 <- lme(swarm_mean_depth ~ temp_mean, random = ~ temp_mean | survey, data = cor.dat.sur, method='ML')
#' summary(fit.4)
#' plot(fit.4)
#' 
#' fit.3 <- lme(Sv_mean_mean ~ temp_mean, random = ~ temp_mean | survey, data = cor.dat.sur, method='ML')
#' summary(fit.3)
#' plot(fit.3)

#'----------------------------------#
#' __N. Extra Experimental Stufff__ #
#'----------------------------------#
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


# Offshore temperature gradient
ggplot(cast.df, aes(x=crossShoreDistance, y=temp_mean)) +
  geom_point(shape=3, alpha=0.5) +
  geom_smooth(method='lm', alpha=.5, se=F) +
  theme_bw() +
  labs(x='Distance offshore (km)', y='Temperature', color='Cluster') +
  facet_wrap(~survey, scales='free')

# Offshore biomass gradient
ggplot(cast.df, aes(x=crossShoreDistance, y=area_Sv_mean_sum)) +
  geom_point(shape=3, alpha=0.5) +
  geom_smooth(method='lm', alpha=.5, se=F) +
  theme_bw() +
  labs(x='Distance offshore (km)', y='Total biomass (Sv mean)') +
  facet_wrap(~survey)

# Check total Biomass and color by offshore
ggplot(cast.df, aes(x=temp_mean, y=area_Sv_mean_sum, color=crossShoreDistance)) +
  geom_point(shape=16, alpha=1, size=2.5) +
  scale_color_cmocean(name='haline', start=.2, end=.8, direction=-1) +
  geom_smooth(method='lm', alpha=.5, se=F) +
  theme_bw() +
  labs(x='Temperature', y='Total biomass (Sv mean)', color='Distance offshore (km)') +
  theme(legend.position = 'top')

# Same but wiht longitude
ggplot(cast.df, aes(x=temp_mean, y=area_Sv_mean_sum, color=lon)) +
  geom_point(shape=16, alpha=1, size=2.5) +
  scale_color_cmocean(name='haline', start=.2, end=.8, direction=-1) +
  geom_smooth(method='lm', alpha=.5, se=F) +
  theme_bw() +
  labs(x='Temperature', y='Total biomass (Sv mean)', color='Longitude') +
  theme(legend.position = 'top')

# Interesting... 
plot(cast.df$lon, cast.df$crossShoreDistance, col=factor(substr(cast.df$station_id,1,2)))
plot(cast.df$lon, cast.df$lat)




ggplot(cast.df, aes(x=temp_mean, y=swarm_mean_depth)) +
  geom_point(size=2, shape=3) +
  geom_smooth(method='lm') +
  facet_wrap(~survey, scales='free')
