CTD <- cast.reader()
surveys <- c("2015_S1","2016_S1","2016_S2",
             "2017_S1","2017_S2","2018_S1",
             "2018_S2","2019_S1","2019_S2")
for (survey in surveys){
  df <- CTD$meta[CTD$meta$survey_id == survey,]
  message("\n",survey)
  message("==========")
  message("Range")
  print((max(date(df$cast_time_local)) - min(date(df$cast_time_local))) + 1)
  message("Sample days")
  print(length(unique(date(df$cast_time_local))))
}
