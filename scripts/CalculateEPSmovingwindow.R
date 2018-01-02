# Calculate treering statistics
# Written by: Isla Myers-Smith (e-mail: isla.myers-smith@ed.ac.uk)
# 22 February 2015

# Load data ---------------------------------------------------------------
# Libraries
library(plyr)

load("workspace/shrubhub.RData")
load("workspace/climate_data.RData")
type <- gsub("[[:digit:]]", "", cru.cld$site)
cru <- cbind.data.frame(cru.cld[1:4], cru.cld[4], type, cru.cld[5], cru.dtr[5], cru.frs[5], cru.pre[5], cru.pet[5], cru.tmn[5], cru.tmp[5], cru.tmx[5], cru.vap[5], cru.wet[5])
colnames(cru) <- c("site", "snum", "year", "month", "season", "type", "cld", "dtr", "frs", "pre", "pet", "tmn", "tmp", "tmx", "vap", "wet")

# Merge data frames
dataall <- merge(growth.full, climate_data)
dataall[dataall == "NaN"] <- "NA"

unfactorize <- function(df){
  for(i in which(sapply(df, class) == "factor")) df[[i]] = as.character(df[[i]])
  return(df)
}

data <- unfactorize(dataall)

# Load soil moisture data
soilm <- read.csv("scripts/users/imyerssmith/shrub_synthesis/soil_moisture.csv", 
                  header=TRUE, stringsAsFactors=FALSE)

# Compute average climate summaries since 1975
climdata <- subset(climate_data, year < 2009 & year >= 1975)
climdata1 <- cbind.data.frame(climdata$site, climdata$snum, climdata$year, climdata$fdf, climdata$amt, climdata$junetmp, climdata$junetmx, climdata$junetmn, climdata$julytmp, climdata$julytmx, climdata$julytmn, climdata$augtmp, climdata$augtmx, climdata$augtmn, climdata$jjtmp, climdata$jjtmx, climdata$jjtmn, climdata$astmp, climdata$astmx, climdata$astmn, climdata$wintmp, climdata$wintmx, climdata$wintmn, climdata$mamtmp, climdata$mamtmx, climdata$mamtmn, climdata$junepre, climdata$julypre, climdata$augpre, climdata$jjpre, climdata$aspre, climdata$winpre, climdata$mampre, climdata$wetd, climdata$sumcloud, climdata$gscru.length)
names(climdata1) <- c("site", "snum", "year", "fdf", "amt", "junetmp", "junetmx", "junetmn", "julytmp", "julytmx", "julytmn", "augtmp", "augtmx", "augtmn", "jjtmp", "jjtmx", "jjtmn", "astmp", "astmx", "astmn", "wintmp", "wintmx", "wintmn", "mamtmp", "mamtmx", "mamtmn", "junepre", "julypre", "augpre", "jjpre", "aspre", "winpre", "mampre", "wetd", "sumcloud", "gscru.length")
samt <- ddply(climdata1, "site", summarise, amt = mean(amt, na.rm=TRUE))
sfdf <- ddply(climdata1, "site", summarise, fdf = mean(fdf, na.rm=TRUE))
sjjtmp <- ddply(climdata1, "site", summarise, jjtmp = mean(jjtmp, na.rm=TRUE))
swinpre <- ddply(climdata1, "site", summarise, winpre = mean(winpre, na.rm=TRUE))
swetd <- ddply(climdata1, "site", summarise, wetd = mean(wetd, na.rm=TRUE))
sjjpre <- ddply(climdata1, "site", summarise, jjpre = mean(jjpre, na.rm=TRUE))
ssumcloud <- ddply(climdata1, "site", summarise, sumcloud = mean(sumcloud, na.rm=TRUE))
sgs <- ddply(climdata1, "site", summarise, gs = mean(gscru.length, na.rm=TRUE))
climsum <- as.data.frame(cbind(as.character(samt$site), unlist(samt$amt), unlist(sfdf$fdf), unlist(sjjtmp$jjtmp), 
                               unlist(swinpre$winpre), unlist(swetd$wetd), unlist(sgs$gs), unlist(sjjpre$jjpre), 
                               unlist(ssumcloud$sumcloud)))
colnames(climsum) <- c("site", "samt", "sfdf", "sjjtmp", "swinpre", "swetd", "sgs", "sjjpre", "ssumcloud")

# Load mixed model comparison analysis
modelsum <- read.csv("scripts/users/imyerssmith/shrub_synthesis/normalized/model_comp.csv", 
                     header=TRUE, stringsAsFactors=FALSE)
siteinfo <- read.csv("scripts/users/imyerssmith/shrub_synthesis/site_info.csv", 
                     header=TRUE, stringsAsFactors=FALSE)
siteinfosum <- as.data.frame(unique(cbind(siteinfo$sgnum, siteinfo$site, siteinfo$site_code, siteinfo$lat, siteinfo$lon)))
colnames(siteinfosum) <- c("sgnum", "site", "site_code", "lat", "lon")
modelsum <- merge(modelsum, climsum, by.x = "site", by.y = "site", all.x = TRUE)
modelsum <- merge(modelsum, soilm)
modelsum <- merge(modelsum, siteinfo)

unfactorize(modelsum)

# Tree ring stats ---------------------------------------------------------
# calculate treering stats with moving window

data <- dataall
a <- 1
d <- 1
length <- length(unique(data$sgnum))
eps_all <- as.data.frame(array(0, c(length*3, 6)))
colnames(eps_all) <- c("year", "site", "sgnum", "samp.depth", "rbar", "eps")
b <- c(1950, 1970, 1990)
c <- c(1970, 1990, 2010)

for (a in 1:length){
  i <- 1
  l <- 1
  data <- subset(dataall, year >= b[a] & year < c[a])
  length1 <- length(unique(data$sgnum))
  length2 <- sum(!duplicated(cbind(data$sgnum, data$year)))
  chronology <- as.data.frame(array(0, c(length2, 4)))
  colnames(chronology) <- c("cstd", "samp.depth", "year", "sgnum")
  eps_values <- as.data.frame(array(0, c(length1, 6)))
  colnames(eps_values) <- c("year", "site", "sgnum", "samp.depth", "rbar", "eps")
  sgnums <- c(sort(unique(data$sgnum)))
  
  for (i in 1:length1){
    data <- subset(dataall, year >= b[a] & year < c[a])
    data <- subset(data, sgnum == sgnums[i])
    
    # Reformat data to replace missing rings indicated with "NA" to the value 0.001
    j <- 1
    k <- 1
    z <- 0
    z2 <- 0
    z3 <- as.data.frame(array(0, c(length(data$rw), 4)))
    colnames(z3) <- c("year", "shrub_id", "rw", "age")
    samples <- c(sort(unique(data$sample)))
    for (j in 1:length(samples)){
      data2 <- subset(data, sample == samples[j])
      
      yr <- seq(min(data2$year), max(data2$year))
      age <- seq(0, (max(data2$year)-min(data2$year)))
      shrubs <- rep(unique(data2$shrub_id)[[1]], length(yr))
      yrsh <- cbind(yr, age, shrubs)
      colnames(yrsh) <- c("year", "age", "shrub_id")
      z = data.frame(data2$year, as.character(data2$shrub_id), data2$rw)
      colnames(z) <- c("year", "shrub_id", "rw")
      z2 <- merge(z, yrsh, all = TRUE)
      z2$rw[is.na(z2$rw)] <- 0.001
      z2$age <- as.numeric(as.character(z2$age))
      z3[k:(k+length(z2$rw)-1), ] <- z2
      z3[k:(k+length(z2$rw)-1), 2] <- unique(data2$shrub_id)
      k <- k + length(z2$rw)
      j <- j+1 
    }
    
    # Reshape the data to fit the "Tucson" treering data format
    xyr <- z3
    xyr$age <- NULL
    xyr2 <- reshape(xyr, timevar = "shrub_id", idvar = "year", direction = "wide")
    rownames(xyr2) <- as.numeric(xyr2$year)
    xyr3 <- xyr2[order(-as.numeric(xyr2$year)), ]
    xyr3$year <- NULL
    
    xag <- z3
    xag$year <- NULL
    xag2 <- reshape(xag, timevar = "shrub_id", idvar = "age", direction = "wide")
    rownames(xag2) <- as.numeric(xag2$age)
    xag3 <- xag2[order(as.numeric(xag2$age)), ]
    xag3$age <- NULL 
    
    # Change detrending method from "Mean" to "Spline" for "FigS6chronologies_mean" or "FigS7chronologies_spline" figures
    xyr4 <- detrend(xyr3, method = c("Mean"))
    chronyr <- chron(xyr4, prefix = "c", biweight = TRUE, prewhiten = FALSE)
    chronyr$year <- as.vector(rownames(chronyr))
    colnames(chronyr) <- c("cstd", "samp.depth", "year")
    chronyr2 <- chronyr[chronyr$cstd > 0.01, ]
    jjtmp <- unique(cbind.data.frame(data$year, data$jjtmp))
    names(jjtmp) <- c("year", "jjtmp")
    chronyr3 <- merge(chronyr2, jjtmp, by.x = "year", by.y = "year", all.x = TRUE, all.y = FALSE)
    chronyr3 <- na.omit(chronyr3)
    
    xag4 <- detrend(xag3, method = c("Mean"))
    chronag <- chron(xag4, prefix = "c", biweight = TRUE, prewhiten = FALSE)
    chronag$age <- as.vector(rownames(chronag))
    colnames(chronag) <- c("cstd", "samp.depth", "age")
    chronag2 <- chronag[chronag$cstd > 0.01, ]
    
    # Calculate treering stats
    gini.coef <- round(gini.coef(xyr4), 1)
    sens <- round(sens1(xyr4), 1)
    xyr5 <- xyr4
    xyr5[is.na(xyr5)] <- 0
    names(xyr5) <- gsub("rw.", "", gsub("_", "", names(xyr4)))
    eps <- EPS.value(xyr5, stc=c(0, 8, 0))
    site <- as.character(unique(data$site_code))
    genus <- if(length(as.character(unique(data$speciescode))) == 1) { as.character(unique(data$speciescode)) } else { if(as.character(unique(data$genus)) == "S") { "SALspp" } else { if(as.character(unique(data$genus)) == "V") { "VACspp" } else { as.character(unique(data$genus)) } } }
    
    # Save chronology
    chron <- as.data.frame(cbind(chronyr2, rep(unique(data$sgnum), length(chronyr2$year))))
    chronology[l:(l+length(chron$year)-1), 1:4] <- chron
    eps_values[i, 1] <- b[a]
    eps_values[i, 2] <- site
    eps_values[i, 3] <- unique(data$sgnum)
    eps_values[i, 4] <- eps$n.trees[1]
    eps_values[i, 5] <- eps$rbar.tot[1]
    eps_values[i, 6] <- eps$eps[1]  
    eps_values <- subset(eps_values, year != 0) 
    
    l <- l + length(chron$year)
    i <- i+1
  }
  
  eps_all[d:(d+length(eps_values$year)-1), 1:6] <- eps_values
  d <- d + length(eps_values$year)
  a <- a+1
}
