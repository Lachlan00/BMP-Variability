# Survey plots
source('R/IMOS.R')
library(lubridate)

# Get survey periods
surveys <- read.csv('data/surveys/meta/survey_times.csv')
surveys$start_UTC <- as.POSIXct(surveys$start_UTC, tz='UTC')
surveys$end_UTC <- as.POSIXct(surveys$end_UTC, tz='UTC')
padding <- days(3)
region = "Syd-Hob"

# Make images
base.dir <- paste0('./output/IMOS/surveys/',region,'/')
for (i in 1:nrow(surveys)){
  # Make dir
  message(surveys$id[i])
  dir.create(paste0(base.dir, surveys$id[i]), showWarnings=F)
  
  # Download images
  for (var in c('SST','chla')){
    dir.create(paste0(base.dir, surveys$id[i], '/', var), showWarnings=F)
    IMOS.OceanCurrent.download(product=var, region=region,
                               daterange=c(surveys$start_UTC[i]-padding, surveys$end_UTC[i]+padding),
                               outdir=paste0(base.dir, surveys$id[i],'/',var))
    
    # Animate
    message('Rendering animation')
    # Make gif
    system(paste0("convert ",base.dir,surveys$id[i],'/',var,
                  "/*.gif ",base.dir,surveys$id[i],"_",var,".gif"))
    # Chnage frame rate
    f.rate <- ifelse(var == 'SST', '10', '60')
    system(paste0("convert -delay ",f.rate,"x100 ",
                  base.dir,surveys$id[i],"_",var,".gif ",
                  base.dir,surveys$id[i],"_",var,".gif"))
    # # Gif to video
    # paste0("ffmpeg -f gif -i ",
    #        base.dir,surveys$id[i],"_",var,".gif ",
    #        base.dir,surveys$id[i],"_",var,".mp4")
  }
}
