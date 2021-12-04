source('R/ROMS-tools.R')

# Load ROMS data
temp <- load.ROMS('temp', period=1994:2016, downsample.factor=10, ROMS.version = 'surface')
salt <- load.ROMS('salt', period=1994:2016, downsample.factor=10, ROMS.version = 'surface')

# Filter region
region <- list(lonmin=150.0056,
               latmin=-36.4825,
               lonmax=150.3628,
               latmax=-36.128)
temp <- temp[temp$lon >= region$lonmin & temp$lon <= region$lonmax &
             temp$lat >= region$latmin & temp$lat <= region$latmax,]
salt <- salt[salt$lon >= region$lonmin & salt$lon <= region$lonmax &
             salt$lat >= region$latmin & salt$lat <= region$latmax,]

# Aggregate by dates
temp <- aggregate(temp, by=list(temp$dt), FUN = function(x) mean(x, na.rm=T))
temp <- temp[,c('dt', 'var')]
salt <- aggregate(salt, by=list(salt$dt), FUN = function(x) mean(x, na.rm=T))
salt <- salt[,c('dt', 'var')]

# combine data
df <- temp
names(df)[2] <- 'temp'
df$salt <- salt$var

# Alt data load
if (FALSE){
  df <- read.csv('data/ANMN/custom/surface-series.csv')
  df <- df[613:nrow(df),]
  df <- df[complete.cases(df),]
  names(df) <- c('dt', 'temp', 'salt')
}
    
# function to take mean over specified time frame
mean_window <- function(df, tf){
  df$group <- rep(seq(from=1, to=floor(nrow(df)/tf)), each=tf)[1:nrow(df)]
  df <- aggregate(df, by=list(df$group), FUN = function(x) mean(x, na.rm=T))
  df <- df[,2:(ncol(df)-1)]
  return(df)
}

# Calc correlations
corr <- data.frame(tf=1:365, cor=NA)
for (i in 1:365){
  dat <- mean_window(df, i)
  corr$cor[i] <- cor(dat$temp, dat$salt)
}

plot(corr$tf[1:179], corr$cor[1:179], type='l')

