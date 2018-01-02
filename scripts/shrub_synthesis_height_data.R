##ARE TALLER SHRUBS MORE CLIMATE SENSITIVE?###
##Code for Team Shrub meeting 18-03-2015######
##Compiled by Sandra (from Isla's code)#######
##############################################


# Libraries ---------------------------------------------------------------

library(plyr)
library(nlme)
library(lme4)
library(ggplot2)

# Data --------------------------------------------------------------------

# Load relevant data
load("workspace_old/shrubhub.RData")
load("workspace_old/climate_data.RData")

siteinfosum <- read.csv("scripts/users/imyerssmith/shrub_synthesis/site_info.csv", 
                        header=TRUE, stringsAsFactors=FALSE)
siteinfosum <- as.data.frame(unique(cbind(siteinfosum$sgnum, siteinfosum$site, siteinfosum$site_code, siteinfosum$lat, siteinfosum$lon)))
colnames(siteinfosum) <- c("sgnum", "site", "site_code", "lat", "lon")

# Merge data frames
data <- merge(growth.full, climate_data)
data[data == "NaN"] <- "NA"

#Unfactorize data
unfactorize <- function(df){
  for(i in which(sapply(df, class) == "factor")) df[[i]] = as.character(df[[i]])
  return(df)
}

data <- unfactorize(data)

#Create the height subset
Hdata <- subset(data, !is.na(height))

#Correct height units
Hdata$height <- sub("m", "", Hdata$height)
Hdata$height <- sub("<", "", Hdata$height)
Hdata$height <- sub(" ", "", Hdata$height)
Hdata$height <- as.numeric(Hdata$height)
Hdata$height[Hdata$height>6] <- Hdata$height[Hdata$height>6]/100
#Some measurements are in cm, other in m: this converts the cm values into m 

heightdata <- cbind(Hdata$site, Hdata$site_code, Hdata$site_full, Hdata$subsite_full, Hdata$pi, Hdata$pi_full, Hdata$latitude_raw, Hdata$longitude_raw, Hdata$elevation_raw, Hdata$study_id, Hdata$soil_moisture, Hdata$herbivory, Hdata$species, Hdata$sex, Hdata$width, Hdata$width2, Hdata$height, Hdata$stem_age, Hdata$stemw, Hdata$meanrw)
heightdata <- as.data.frame(heightdata)
colnames(heightdata) <- c("site","site_code","site_full","subsite_full","pi","pi_full","latitude_raw","longitude_raw","elevation_raw","study_id","soil_moisture","herbivory","species","sex","width","width2","height","stem_age","stemw","meanrw")

MaxYear <- cbind(as.character(Hdata$site), as.numeric(as.character(Hdata$year)))
colnames(MaxYear) <-c("site", "year")
MaxYear <- as.data.frame(MaxYear)
MaxYear <- unfactorize(MaxYear)
MaxYear <- aggregate(as.numeric(as.character(MaxYear$year)), by=list(as.character(MaxYear$site)), FUN=max, NA.rm=TRUE)
colnames(MaxYear) <-c("site", "year_measured")
heightdata <- unique(heightdata)
heightdata1 <- merge(heightdata, MaxYear)

write.csv(heightdata1, file = "scripts/users/imyerssmith/shrub_synthesis/heightdata.csv", row.names=FALSE)
