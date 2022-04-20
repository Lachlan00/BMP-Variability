# CTD profiles
source("R/CTD.R")
CTD <- cast.reader()
surveys <- unique(CTD$meta$survey_id)
CTD.plot.survey.mean(CTD, c("2015_S1"))
CTD.plot.survey.mean(CTD, surveys) + facet_wrap(~survey_id, scales="free")

# mean for all years
mean.df <- aggregate(CTD$data[,c('temperature', 'salinity')], by=list(CTD$data$id), mean)
# put survey back in
mean.df$survey_id <- sapply(mean.df$Group.1, function(id) CTD$meta$survey_id[CTD$meta$file_name == id])

# mean for surveys
mean.df <- aggregate(mean.df[,c('temperature', 'salinity')], by=list(mean.df$survey_id), mean)
names(mean.df)[1] <- 'survey_id'

summary(mean.df)
