# Functions to interact with IMOS
library(rvest)
library(xml2)
library(lubridate)

# Datetime not yet supported
IMOS.OceanCurrent.download <- function(product=c('SST','chla'), region='SNSW',
                                       daterange=NULL, datetimes=NULL,
                                       outdir='./output/IMOS/OceanCurrent', prefix=NULL){
  # Match arguments
  product <- match.arg(product)
  # Checks
  if (!is.null(datetimes) & !is.null(daterange)){
    stop('Please pass either a daterange or datetimes. Not both.')
  }
  
  # Set query periods
  if (!is.null(daterange)){
    daterange <- as.POSIXct(daterange, format='%Y-%m-%d', tz='UTC')
  }
  
  # Set parameters
  base.url <- 'http://oceancurrent.imos.org.au/'
  
  # Set region directory
  if (product == 'SST'){
    region_dir <- region
  } else if (product == 'chla') {
    region_dir <- paste0(region,'_chl')
  }
  
  # Retrive avliable dates
  if (!is.null(daterange)){
    date.set <- as.list(unique(format.Date(daterange, '%Y')))
  } else {
    date.set <- as.list(unique(format.Date(datetimes, '%Y')))
  }
  
  for (i in 1:length(date.set)){
    year <- date.set[i]
    if (as.numeric(year) > 2018)
      year <- 'index.html'
    date.page <- read_html(paste0(base.url,region_dir,'/',year))
    date.page <- xml_child(date.page, 2)
    if (year == 'index.html' | product != 'SST')
      date.page <- xml_child(date.page, length(xml_children(date.page)))
    date.href <- unlist(lapply(xml_children(date.page), function(elem) html_attr(elem, 'href')))
    date.href <- date.href[grepl('^[0-9]{10}\\.html', date.href)]
    date.set[[i]] <- as.POSIXct(date.href, format='%Y%m%d%H', tz='UTC')
  }
  date.set <- do.call("c", date.set)
  
  # Reduce date.set to dates within daterange
  if (!is.null(daterange)){
    date.set <- date.set[date.set > daterange[1] & date.set < daterange[2]]
  } else {
    date.set <- date.set[date(date.set) %in% unique(date(datetimes))]
  }
  
  # Create the request strings based on product
  # Any image <= 2018 has a year sub directory
  date.urls <- paste0(format(date.set,'%Y%m%d%H'),'.gif')
  date.urls[year(date.set) <= 2018] <- paste0(format(date.set[year(date.set) <= 2018],'%Y'),'/',
                                              date.urls[year(date.set) <= 2018])
  if (!is.null(prefix))
    prefix <- paste0(prefix,'_')
  # Make a dataframe of requests urls and output names
  df.download <- data.frame(url=paste0(base.url,region_dir,'/',date.urls),
                            fn=paste0(outdir,'/',prefix,region,'_',product,'_',
                                      format(date.set,'%Y%m%d%H.gif')),
                            stringsAsFactors=F)
  
  # Download images
  message('\nDownlaoding ',nrow(df.download),' IMOS OceanCurrent ',product,
          ' images for ',region,'...')
  pb <- txtProgressBar(max=nrow(df.download), style=3)
  for (i in 1:nrow(df.download)){
    download.file(df.download[i,1], df.download[i,2], mode='wb', quiet=T)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  message('Images saved to "',outdir,'"')
}
