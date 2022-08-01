setwd("/Users/lachlanphillips/Development/PhD/paper-repos/BMP-Variability")
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
cor.dat.orignames <- cor.dat
names(cor.dat) <- c('Temperature', 'Salinity',
                    'Agg. Coverage', 'Mean Agg. Density', 'Mean Agg. Depth')


#'---------------------#
#' __3. PCA analysis__ #
#'---------------------#
# Build PCA model
pc <- prcomp(cor.dat, center = TRUE, scale. = TRUE)
ggbiplot(pc, ellipse=TRUE, obs.scale = 1, var.scale = 1, var.axes=TRUE,
         groups=cast.df$survey)
print(pc)
summary(pc)

#' # Check scatter of PCA (pros removal of multicollinearity)
#' pairs.panels(pc$x,
#'              gap = 0,
#'              bg = as.factor(cor.dat.years$Group.1),
#'              pch=21, lm=T, ci=T,
#'              ellipses=F, hist.col='grey')
#' 
#' #' __Part 2 - PCA axis Boxplot__ #
#' pca.box <- cbind(pc$x[,1], cast.df$survey)
#' pca.box <- as.data.frame(pca.box)
#' pca.box$V1 <- as.numeric(pca.box$V1)
#' pca.box$V2 <- as.factor(pca.box$V2)
#' names(pca.box) <- c('PCA1', 'Survey')
#' 
#' ggplot(pca.box, aes(x=Survey, y=PCA1)) + 
#'   geom_boxplot(fill='grey') +
#'   theme_bw() +
#'   labs(x=NULL, y='PCA axis 1')
#' 
#' #' __Part 3 - PCA spatial__ #
#' # Add PCA axis 1 to cast.df
#' cast.df$PCA1 <- as.numeric(pc$x[,1])
#' 
#' # Load grid 
#' stations <- read.csv('./data/transects/ghostStations.csv')
#' grid <- read.csv('./data/transects/oceanGrid.csv')
#' # Create a set of shape files
#' grid.polys <- grid.polygons(grid, ids=stations$id)
#' 
#' # get cell numbers
#' cast.df$cell <- as.numeric(factor(cast.df$station_id, levels=stations$id))
#' # calc ffor each survey
#' cell.PCA <- aggregate(cbind(PCA1, lon, lat) ~ survey + cell, data=cast.df, FUN=mean)
#' # make coords
#' cell.PCA$col <- cell.PCA$cell %% 5
#' cell.PCA$col[cell.PCA$col == 0] <- 5
#' cell.PCA$row <- ceiling(cell.PCA$cell / 5)
#' 
#' ggplot(cell.PCA, aes(x=col, y=row, fill=PCA1)) +
#'   geom_tile() +
#'   facet_wrap(~survey) +
#'   scale_fill_cmocean('thermal', direction=-1)


#'----------------------------------#
#' __3. K Means Clustering on PCA__ #
#'----------------------------------#
#' Note!! Earlier plots using fviz_cluster were actually doing a PCA on a PCA!
# Apply K means cluster to PCA
pc.transform <- as.data.frame(pc$x)
# detrmine optimal clusters (2 is optimal)
fviz_nbclust(pc.transform, kmeans, method = 'silhouette')
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
  geom_point(aes(x=x, y=y, color=survey, shape=cluster), size=2, alpha=.9) +
  #scale_color_cmocean(name='phase',  discrete=T) +
  #scale_shape_manual(values=c(3,17)) + 
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
            color = 'darkred', size = 3.5, alpha=1) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14))

ggsave('./figures/PCA/pca-cluster.png', p.pca, width=9, height=6)

# Cluster table
table(cast.df$survey, cast.df$cluster)

#'---------------------------------------#
#' __4. Multinominal regression on PCA__ #
#'---------------------------------------#
# # Add clusters to cor.dat
# cor.dat$cluster <- cast.df$cluster
# Divide data into training and testing
ind <- sample(2, nrow(cor.dat),
              replace = TRUE,
              prob = c(0.8, 0.2))
cor.dat.train <- cor.dat[ind==1,]
cor.dat.test <- cor.dat[ind==2,]

# Build PCA model
pc.mn <- prcomp(cor.dat.train, center = TRUE, scale. = TRUE)

# Check preidctions on testing set
trn <- predict(pc.mn, cor.dat.train)
trn <- data.frame(trn, survey=cast.df$survey[ind==1], cluster=cast.df$cluster[ind==1])
tst <- predict(pc.mn, cor.dat.test)
tst <- data.frame(tst, survey=cast.df$survey[ind==2], cluster=cast.df$cluster[ind==2])

# Multinominal regression
mn.reg.survey <- multinom(survey ~ PC1 + PC2 + PC3 + PC4 + PC5, data=trn)
summary(mn.reg.survey)
mn.reg.cluster <- multinom(cluster ~ PC1 + PC2 + PC3 + PC4 + PC5, data=trn)
summary(mn.reg.cluster)

# Confusion Matrices
# SURVEY
# Training
message('\NResults: SURVEY')
p.trn <- predict(mn.reg.survey, trn)
tab <- table(p.trn, trn$survey)
tab
print(paste0('Model accuracy (TRAIN): ',round((sum(diag(tab))/sum(tab))*100,2),'%'))
# Testing
p.tst <- predict(mn.reg.survey, tst)
tab <- table(p.tst, tst$survey)
tab
print(paste0('Model accuracy (TEST): ',round((sum(diag(tab))/sum(tab))*100,2),'%'))

# CLUSTER
# Training
message('\NResults: CLUSTER')
p.trn <- predict(mn.reg.cluster, trn)
tab <- table(p.trn, trn$cluster)
tab
print(paste0('Model accuracy (TRAIN): ',round((sum(diag(tab))/sum(tab))*100,2),'%'))
# Testing
p.tst <- predict(mn.reg.cluster, tst)
tab <- table(p.tst, tst$cluster)
tab
print(paste0('Model accuracy (TEST): ',round((sum(diag(tab))/sum(tab))*100,2),'%'))


#'-----------------------------#
#' __5. T/S plot__ #
#'-----------------------------#
p.TS <- ggplot(cast.df, aes(x=salt_mean, y=temp_mean, fill=survey, color=survey, group=survey)) +
  geom_point(shape=3, alpha=.7, size=2) +
  geom_smooth(method='glm', alpha=.3, se=F) +
  # geom_smooth(method='lm', mapping=aes(x=salt_mean, y=temp_mean),
  #             inherit.aes = F, color='black', se=F, linetype='dashed', alpha=.0) +
  geom_line(mapping=aes(x=salt_mean, y=temp_mean), 
            stat="smooth", method = "lm",
            linetype ="dashed",
            alpha = 0.9, size=1.1, 
            inherit.aes = F) +
  theme_bw() +
  labs(x='Mean Salinity (PSU)', y='Mean Temperature (°C)', color='Survey', fill='Survey') +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14))

ggsave('./figures/TS/TS-scatter.png', p.TS, width=8, height=5.5)


#'-------------------------------------------#
#' __6. Plot Temp against Agg for Clusters__ #
#'-------------------------------------------#
point_alpha=.45
se=T
se.alpha=.3
# AGGs
p1 <- ggplot(cast.df, aes(x=temp_mean, y=swarm_mean_depth, color=cluster, group=cluster)) +
  geom_point(shape=3, alpha=point_alpha) +
  geom_smooth(method='glm', alpha=se.alpha, se=se) + #, method.args = list(family = Gamma(link = "log"))) +
  theme_bw() +
  scale_y_reverse() + 
  labs(x='Mean Temperature (°C)', y='Mean Agg. Depth (m)', color='Cluster') +
  ggtitle('(A) Mean Aggregation Depth')
  # stat_ellipse(aes(x=temp_mean, y=swarm_mean_depth, group=survey), type = "norm", 
  #              alpha=.2, linetype='dashed', size=.4, segments=51, inherit.aes = F)

p2 <- ggplot(cast.df, aes(x=temp_mean, y=area_corrected_area_sum, color=cluster, group=cluster)) +
  geom_point(shape=3, alpha=point_alpha) +
  geom_smooth(method='glm', alpha=se.alpha, se=se) +
  theme_bw() +
  labs(x='Mean Temperature (°C)', y='Mean Agg. Coverage (%)', color='Cluster') +
  ggtitle('(B) Mean Aggregation Coverage')

p3 <- ggplot(cast.df, aes(x=temp_mean, y=Sv_mean_mean, color=cluster, group=cluster)) +
  geom_point(shape=3, alpha=point_alpha) +
  geom_smooth(method='glm', alpha=se.alpha, se=se) +
  theme_bw() +
  labs(x='Mean Temperature (°C)', y=bquote('Mean Agg. Density ('~S[v]~' mean)'), color='Cluster') +
  ggtitle('(C) Mean Aggregation Density')
  # stat_ellipse(aes(x=temp_mean, y=Sv_mean_mean, group=survey), type = "norm", 
  #               alpha=.2, linetype='dashed', size=.4, segments=51, inherit.aes = F)

p4 <- ggplot(cast.df, aes(x=temp_mean, y=area_Sv_mean_sum, color=cluster, group=cluster)) +
  geom_point(shape=3, alpha=point_alpha) +
  geom_smooth(method='glm', alpha=se.alpha, se=se) +
  theme_bw() +
  labs(x='Mean temperature (°C)', y=bquote('Agg. Density ('~S[v]~' mean / '~km^2~')'), color='Cluster') +
  ggtitle('(D) Biomass')

p1 <- p1 + labs(x='')
p3 <- p3 + labs(x='')
p.ls <- list(p1, p2, p3)#, p4)
p.ls <- lapply(p.ls, function(p) p + theme(legend.position = 'none'))
p.agg <- ggarrange(plotlist = p.ls, nrow=1)#, common.legend = T, legend = 'bottom')
#names(cast.df)
ggsave('./figures/TS/agg-cluster-scatter.png', p.agg, width=12, height=4)


#'----------------------#
#' __Cluster boxplots__ #
#'----------------------#
dat <- melt(cbind(cor.dat.orignames, cast.df[,c('survey','cluster')]), id.vars = c('survey','cluster'))
ggplot(dat, aes(x=cluster,  y=value, fill=cluster)) +
  geom_boxplot() +
  facet_wrap(~variable, scales='free')

anova(area_corrected_area_sum ~ cluster, data=dat)
