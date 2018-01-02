# Shrub growth synthesis figures script:
# This code creates all the figures and supplementary analyses for the shrub synthesis manuscript.
# Written by: Isla Myers-Smith (e-mail: isla.myers-smith@ed.ac.uk)
# 26 June 2014

# Load data ---------------------------------------------------------------
# Libraries
library(plyr)

# Load data
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
modelsum <- read.csv("scripts/users/imyerssmith/shrub_synthesis/model_compn.csv", 
                     header=TRUE, stringsAsFactors=FALSE)
siteinfo <- read.csv("scripts/users/imyerssmith/shrub_synthesis/site_info.csv", 
                  header=TRUE, stringsAsFactors=FALSE)
siteinfosum <- as.data.frame(unique(cbind(siteinfo$sgnum, siteinfo$site, siteinfo$site_code, siteinfo$lat, siteinfo$lon)))
colnames(siteinfosum) <- c("sgnum", "site", "site_code", "lat", "lon")
modelsum <- merge(modelsum, climsum, by.x = "site", by.y = "site", all.x = TRUE)
modelsum <- merge(modelsum, soilm)
modelsum <- merge(modelsum, siteinfo)

unfactorize(modelsum)

# Climate data check ------------------------------------------------------

# Load station data and site info
stn <- read.csv("scripts/users/imyerssmith/shrub_synthesis/station_data_Mar2013.csv", 
                header=TRUE, stringsAsFactors=FALSE)
sitenames <- as.data.frame(unique(cbind(as.character(siteinfo$site), as.character(siteinfo$site_code))))
colnames(sitenames) <- c("site", "site_code")

# Climate data check
climcheck <- merge(stn, sitenames, by.x = "site_code", by.y = "site_code")
climcheck2 <- merge(climcheck, cru,  by.x = c("site", "year", "month"), by.y = c("site", "year", "month"))

i <- 1
j <- 1
length <- length(unique(climcheck2$Site))*12
ccsum <- as.data.frame(array(0, c(length, 5)))
colnames(ccsum) <- c("site", "month", "R", "p", "n")

for (i in 1:12){ 
  k <- j+19
  cc <- subset(climcheck2, month==i)
  rcclm <- lapply(split(cc, as.factor(cc$site_code)), function(cc) {summary(lm(cc$tmp~cc$stntmp))$adj.r.squared} ) 
  ncclm <- lapply(split(cc, as.factor(cc$site_code)), function(cc) {length(unique(cc$year))} ) 
  cclm <- lapply(split(cc, as.factor(cc$site_code)), function(cc) {summary.lm(lm(cc$tmp~cc$stntmp))$coefficients} ) 
  pcclm <- lapply(cclm, "[[", 8)
  ccsum[j:k, 1] <- unlist(names(cclm))
  ccsum[j:k, 2] <- i
  ccsum[j:k, 3] <- sqrt(unlist(rcclm))
  ccsum[j:k, 4] <- unlist(pcclm)  
  ccsum[j:k, 5] <- unlist(ncclm)
  i <- i + 1
  j <- j+20
}

write.table(ccsum, file="scripts/users/imyerssmith/shrub_synthesis/climatecheck.csv", sep=", ", row.names=FALSE)

# Linear models all figure ---------------------------------------------------------
png(file="scripts/users/imyerssmith/shrub_synthesis/FigS8lmall.png", width=2000, height=4000)
par(mfrow=c(12, 4), mar=c(8, 8, 8, 8), oma=c(2, 2, 2, 2), mgp=c(5, 2, 0))

i <- 1
j <- 1
data <- dataall
length <- length(unique(data$sgnum))

for (i in 1:length){
  data <- dataall
  data <- subset(data, sgnum == i)
  length2 <- length(unique(data$shrub_num))
  site <- as.character(unique(data$site_code))
  genus <- if(length(as.character(unique(data$speciescode))) == 1) { as.character(unique(data$speciescode)) } else { if(as.character(unique(data$genus)) == "S") { "SALspp" } else { if(as.character(unique(data$genus)) == "V") { "VACspp" } else { as.character(unique(data$genus)) } } }
  j <- 1
  
  par(xpd=TRUE)
  plot(data$rw ~ data$jjtmp, col = "#00000010", pch = 19, xlab = expression(paste("June-July Temperture (", degree, "C)")), ylab = "Stem Growth (mm)", xlim = c(min(na.omit(data$jjtmp)), max(na.omit(data$jjtmp))), ylim = c(0, max(data$rw)), cex.lab=2.8, cex.axis=2.5, cex=1)
  title(main=paste(bquote(.(site)), ":", bquote(.(genus))), cex.main = 2.5, font.main = 1)
  
  for (j in 1:length2){
    d <- data
    d <- subset(d, sample == j)
    d <- na.omit(cbind.data.frame(d$rw, d$jjtmp))
    names(d) <- c("rw", "jjtmp")
    lm <- lm(d$rw ~ d$jjtmp)
    par(xpd=FALSE)
    if(summary(lm)$coefficients[[8]] == "NA" | summary(lm)$coefficients[[8]] == "NaN") {
      abline(lm, col = "#00000050", lty = 1, lwd = 2)      
    } else {
      if(summary(lm)$coefficients[[8]]>0.05) {
        abline(lm, col = "#00000050", lty = 1, lwd = 2)
      } else {
        abline(lm, col = "#DF0101", lty = 1, lwd = 2)      
      }
    }
    j <- j + 1
  }  
  i <- i + 1
}

dev.off() 

# Chronology figures -------------------------------------------------------
library(dplR)
library(detrendeR)

# Change file name to either "FigS5chronologies_age", "FigS6chronologies_mean" or "FigS7chronologies_spline" and relevant code
# to make different versions of the figures
png(file="scripts/users/imyerssmith/shrub_synthesis/FigS5chronologies_age.png", width=2000, height=4000)
par(mfrow=c(12, 4), mar=c(8, 8, 8, 8), oma=c(2, 2, 2, 2), mgp=c(5, 2, 0))

i <- 1
l <- 1
data <- dataall
length <- length(unique(data$sgnum))
length2 <- sum(!duplicated(cbind(data$sgnum, data$year)))
chronology <- as.data.frame(array(0, c(length2, 4)))
colnames(chronology) <- c("cstd", "samp.depth", "year", "sgnum")
eps_values <- as.data.frame(array(0, c(length, 5)))
colnames(eps_values) <- c("site", "sgnum", "samp.depth", "rbar", "eps")

for (i in 1:length){
  data <- dataall
  data <- subset(data, sgnum == i)
  # Reformat data to replace missing rings indicated with "NA" to the value 0.001
  j <- 1
  k <- 1
  z3 <- as.data.frame(array(0, c(length(data$rw), 4)))
  colnames(z3) <- c("year", "shrub_id", "rw", "age")
  for (j in 1:length(unique(data$shrub_id))){
    data2 <- subset(data, sample == j)
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
    j <- j+1
    k <- k + length(z2$rw)
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
  eps_values[i, 1] <- site
  eps_values[i, 2] <- unique(data$sgnum)
  eps_values[i, 3] <- eps$n.trees[1]
  eps_values[i, 5] <- eps$rbar.tot[1]
  eps_values[i, 4] <- eps$eps[1]
  
  # Use this code to plot the "FigS5chronologies_age" figure
    par(xpd=TRUE)
    x <- as.numeric(as.character(chronag2$age))+5
    y <- chronag2$cstd
    z <- chronag2$samp.depth
    plot(y ~ x, col="#00000000", pch = 20, cex = 1, xlab = "Age (years)", ylab = "Stem Growth Index", cex.lab=2.8, cex.axis=2.5, xlim = c(max(x), 0), ylim = c(0, 3))
    lines(y ~ x, col="#000000", lwd = 2, xlim = c(max(x), 0), ylim = c(0, 3))
    title(main=paste(bquote(.(site)), ":", bquote(.(genus))), cex.main = 2.5, font.main = 1)
    par(new=TRUE)
    plot(z ~ x, col="#00000000", pch = 20, cex = 1, xlab="", ylab="", xaxt='n', yaxt='n', xlim = c(max(x), 0), ylim = c(0, (max(z)*5)))
    polygon(c(min(x), x, max(x)), c(0, z, 0),  col = "#00000050", border = NA)
    axis(4, cex=2.5, cex.lab=2.8, cex.axis=2.5)
    text(-(max(x)/5), ((max(z)*5)/2), "Sample depth", srt = 270, cex=2.8, xpd = TRUE)
  
 # Use this code to make the "FigS6chronologies_mean" figures
#     par(xpd=TRUE)
#     plot(chronyr3$cstd ~ chronyr3$year, col="#00000000", pch = 20, cex = 1, xlab = "Time (years)", ylab = "Stem Growth Index", cex.lab=2.8, cex.axis=2.5, xlim = c(1950, 2009))
#     par(xpd=FALSE)
#     lines(chronyr3$cstd ~ chronyr3$year, col="#00000080", lwd = 2)
#     title(main=paste(bquote(.(site)), ":", bquote(.(genus)), "rbar=", bquote(.(eps$rbar.tot)), ", EPS=", bquote(.(eps$eps))), cex.main = 2.5, font.main = 1)
#     par(new=TRUE)
#     plot(chronyr3$samp.depth ~ chronyr3$year, col="#00000000", pch = 20, cex = 1, xlab="", ylab="", xaxt='n', yaxt='n', xlim = c(1950, 2009), ylim = c(0, (max(chronyr3$samp.depth)*5)))
#     polygon(c(min(chronyr3$year), chronyr3$year, max(chronyr3$year)), c(0, chronyr3$samp.depth, 0),  col = "#00000050", border = NA)
#     axis(4, cex=2.5, cex.lab=2.8, cex.axis=2.5)
#     text(2021, ((max(chronyr3$samp.depth)*5)/2), "Sample depth", srt = 270, cex=2.8, xpd = TRUE)
  
  # Use this code to make the "FigS7chronologies_spline" figures
#     par(xpd=TRUE)
#     plot(chronyr3$jjtmp ~ chronyr3$year, col="#00000000", pch = 20, cex = 1, axes=F, xlab="Time (years)", ylab="Stem Growth Index", cex.lab=2.8, cex.axis=2.5, xlim = c(1950, 2009), ylim = c(min(chronyr3$jjtmp), max(chronyr3$jjtmp)))
#     par(xpd=FALSE)  
#     lines(chronyr3$jjtmp ~ chronyr3$year, col="#CC003380", lwd = 2)
#     axis(4, cex.lab=2.8, cex.axis=2.5)
#     text(2021, mean(chronyr3$jjtmp), expression(paste("Jun. - Jul. Temp. (", degree, "C)")), srt = 270, cex=2.8, xpd = TRUE)
#     par(new=TRUE)
#     plot(chronyr3$cstd ~ chronyr3$year, col="#00000000", pch = 20, cex = 1, xlab="", ylab="", cex.lab=2.8, cex.axis=2.5, xlim = c(1950, 2009))
#     lines(chronyr3$cstd ~ chronyr3$year, col="#00000080", lwd = 2)
#     title(main=paste(bquote(.(site)), ":", bquote(.(genus)), "rbar=", bquote(.(eps$rbar.tot)), ", EPS=", bquote(.(eps$eps))), cex.main = 2.5, font.main = 1)
  
  i <- i+1
  l <- l + length(chron$year)
}

dev.off()

# Chronology correlations -------------------------------------------------
chronology2 <- merge(chronology, siteinfo)
dendro <- merge(chronology2, climate_data, by.x = c("site", "year"), by.y = c("site", "year"), all.x = TRUE)
i <- 1
j <- 1
k <- 1
l <- 34
length <- length(unique(dendro$sgnum))
coeff <- as.data.frame(array(0, c(34, 4)))
colnames(coeff) <- c("model", "estimate", "se", "aic")
coeff2 <- as.data.frame(array(0, c(length*34, 5)))
colnames(coeff2) <- c("site", "model", "estimate", "se", "aic")

dendron <- dendro

# Normalize
normalize <- function(x) { scale(x, scale = TRUE) }

dendron <- ddply(dendron, .(site), transform, nfdf = normalize(fdf))
dendron <- ddply(dendron, .(site), transform, namt = normalize(amt))
dendron <- ddply(dendron, .(site), transform, njunetmp = normalize(junetmp))
dendron <- ddply(dendron, .(site), transform, njunetmx = normalize(junetmx))
dendron <- ddply(dendron, .(site), transform, njunetmn = normalize(junetmn))
dendron <- ddply(dendron, .(site), transform, njulytmp = normalize(julytmp))
dendron <- ddply(dendron, .(site), transform, njulytmx = normalize(julytmx))
dendron <- ddply(dendron, .(site), transform, njulytmn = normalize(julytmn))
dendron <- ddply(dendron, .(site), transform, naugtmp = normalize(augtmp))
dendron <- ddply(dendron, .(site), transform, naugtmx = normalize(augtmx))
dendron <- ddply(dendron, .(site), transform, naugtmn = normalize(augtmn))
dendron <- ddply(dendron, .(site), transform, njjtmp = normalize(jjtmp))
dendron <- ddply(dendron, .(site), transform, njjtmx = normalize(jjtmx))
dendron <- ddply(dendron, .(site), transform, njjtmn = normalize(jjtmn))
dendron <- ddply(dendron, .(site), transform, nastmp = normalize(astmp))
dendron <- ddply(dendron, .(site), transform, nastmx = normalize(astmx))
dendron <- ddply(dendron, .(site), transform, nastmn = normalize(astmn))
dendron <- ddply(dendron, .(site), transform, nwintmp = normalize(wintmp))
dendron <- ddply(dendron, .(site), transform, nwintmx = normalize(wintmx))
dendron <- ddply(dendron, .(site), transform, nwintmn = normalize(wintmn))
dendron <- ddply(dendron, .(site), transform, nmamtmp = normalize(mamtmp))
dendron <- ddply(dendron, .(site), transform, nmamtmx = normalize(mamtmx))
dendron <- ddply(dendron, .(site), transform, nmamtmn = normalize(mamtmn))
dendron <- ddply(dendron, .(site), transform, njunepre = normalize(junepre))
dendron <- ddply(dendron, .(site), transform, njulypre = normalize(julypre))
dendron <- ddply(dendron, .(site), transform, naugpre = normalize(augpre))
dendron <- ddply(dendron, .(site), transform, njjpre = normalize(jjpre))
dendron <- ddply(dendron, .(site), transform, naspre = normalize(aspre))
dendron <- ddply(dendron, .(site), transform, nwinpre = normalize(winpre))
dendron <- ddply(dendron, .(site), transform, nmampre = normalize(mampre))
dendron <- ddply(dendron, .(site), transform, nsumcloud = normalize(sumcloud))
dendron <- ddply(dendron, .(site), transform, nwetd = normalize(wetd))
dendron <- ddply(dendron, .(site), transform, ngscru.length = normalize(gscru.length))
dendron <- ddply(dendron, .(site), transform, nspring.tdd = normalize(spring.tdd))

for (i in 1:length){ 
  data <- dendron
  data <- subset(data, year >= 1950)
  data <- subset(data, sgnum==i)

  m1n <- glm(cstd ~ njjtmp, data = data)  
  m2n <- glm(cstd ~ njjtmx, data = data)
  m3n <- glm(cstd ~ njjtmn, data = data)
  m4n <- glm(cstd ~ njunetmp, data = data)
  m5n <- glm(cstd ~ njunetmx, data = data)
  m6n <- glm(cstd ~ njunetmn, data = data)
  m7n <- glm(cstd ~ njulytmp, data = data)
  m8n <- glm(cstd ~ njulytmx, data = data)
  m9n <- glm(cstd ~ njulytmn, data = data)
  m10n <- glm(cstd ~ nastmp, data = data)
  m11n <- glm(cstd ~ nastmx, data = data)
  m12n <- glm(cstd ~ nastmn, data = data)
  m13n <- glm(cstd ~ nmamtmp, data = data)
  m14n <- glm(cstd ~ nmamtmx, data = data)
  m15n <- glm(cstd ~ nmamtmn, data = data)
  m16n <- glm(cstd ~ nwintmp, data = data)
  m17n <- glm(cstd ~ nwintmx, data = data)
  m18n <- glm(cstd ~ nwintmn, data = data)
  m19n <- glm(cstd ~ njunepre, data = data)
  m20n <- glm(cstd ~ njulypre, data = data)
  m21n <- glm(cstd ~ naugpre, data = data)
  m22n <- glm(cstd ~ njjpre, data = data)
  m23n <- glm(cstd ~ naspre, data = data)
  m24n <- glm(cstd ~ nwinpre, data = data)
  m25n <- glm(cstd ~ nmampre, data = data)
  m26n <- glm(cstd ~ njunetmp + nwinpre, data = data)
  m27n <- glm(cstd ~ njunetmp + njunepre, data = data)
  m28n <- glm(cstd ~ njulytmp + nwinpre, data = data)
  m29n <- glm(cstd ~ njulytmp + njulypre, data = data)
  m30n <- glm(cstd ~ naugtmp + nwinpre, data = data)
  m31n <- glm(cstd ~ naugtmp + naugpre, data = data)
  m32n <- glm(cstd ~ njjtmp + nwinpre, data = data)
  m33n <- glm(cstd ~ njjtmp + njjpre, data = data)
  mNULLn <- glm(cstd ~ 1, data = data)
  
modeln <- list(summary(m1n), summary(m2n), summary(m3n), summary(m4n), summary(m5n), summary(m6n), summary(m7n), summary(m8n), summary(m9n), summary(m10n), summary(m11n), summary(m12n), summary(m13n), summary(m14n), summary(m15n), summary(m16n), summary(m17n), summary(m18n), summary(m19n), summary(m20n), summary(m21n), summary(m22n), summary(m23n), summary(m24n), summary(m25n), summary(m26n), summary(m27n), summary(m28n), summary(m29n), summary(m30n), summary(m31n), summary(m32n), summary(m33n), summary(mNULLn))

for (j in 1:34) { 
  coef <- as.data.frame(modeln[[j]]$coefficients)
  #model name
  coeff[j, 1] <- if(is.na(gsub("scale\\(", "", gsub(", scale = FALSE)", "", rownames(coef)[2]))))   { 
    "NULL" 
  } else if(is.na(rownames(coef)[3]))   { 
    gsub("scale\\(", "", gsub(", scale = FALSE)", "", rownames(coef)[2]))
  } else   { 
    paste(gsub("scale\\(", "", gsub(", scale = FALSE)", "", rownames(coef)[2])), "-", gsub("scale\\(", "", gsub(", scale = FALSE)", "", rownames(coef)[3]))) 
  }
  #estimate
  coeff[j, 2] <- coef[2, 1]
  #se
  coeff[j, 3] <- coef[2, 2]
  #AIC
  coeff[j, 4] <- modeln[[j]]$aic
  j <- j + 1
}
coeff2[k:l, 1] <- unique(data$sgnum)
coeff2[k:l, 2:5] <- coeff
i <- i + 1
j <- 1
k <- k + 34
l <- l + 34
}

modelsglm <- as.data.frame(cbind(coeff2$site, coeff2$model, coeff2$aic, coeff2$estimate, coeff2$se))
colnames(modelsglm) <- c("site", "model", "aic", "estimate", "se")

# Print models
save(modelsglm, file="scripts/users/imyerssmith/shrub_synthesis/modelsglm.RData")
write.csv(modelsglm, file = "scripts/users/imyerssmith/shrub_synthesis/modelsglm.csv", row.names=FALSE)

# select models for comparison
models_sort <- arrange(modelsglm, site, aic)
models_sort$site <- as.numeric(as.character(models_sort$site))
models_sort$model <- gsub("data\\$", "", models_sort$model)
best_model <- arrange(ddply(models_sort, .(site), function(x)x[which.min(x$aic), ]), site)
NULL_model <- subset(models_sort, model == "NULL")
NULL_model <- NULL_model[!duplicated(NULL_model$site), ]
summer_model <- subset(models_sort, model == "njjtmp" | model == "njjtmx" | model == "njjtmn" | model == "njunetmp" | model == "njunetmx" | model == "njunetmn" | model == "njulytmp" | model == "njulytmx" | model == "njulytmn")
summer_model <- arrange(ddply(summer_model, .(site), function(x)x[which.min(x$aic), ]), site)

# calculate delta aic values
model_compglm <- merge(best_model, NULL_model, by.x = "site", by.y = "site")
model_compglm <- merge(model_compglm, summer_model, by.x = "site", by.y = "site")
model_compglm$se.y <- NULL
model_compglm$model.y <- NULL
colnames(model_compglm) <- c("site", "nmodel_best", "naic_best", "nest_best", "nse_best", "naic_null", "nest_null", "nmodel_summer", "naic_summer", "nest_summer", "nse_summer")
model_compglm$naic_best <- as.numeric(as.character(model_compglm$naic_best))
model_compglm$naic_null <- as.numeric(as.character(model_compglm$naic_null))
model_compglm$naic_summer <- as.numeric(as.character(model_compglm$naic_summer))
model_compglm <- cbind(model_compglm, (model_compglm$naic_null-model_compglm$naic_best), (model_compglm$naic_null-model_compglm$naic_summer))
colnames(model_compglm) <- c("site", "nmodel_best", "naic_best", "nest_best", "nse_best", "naic_null", "nest_null", "nmodel_summer", "naic_summer", "nest_summer", "nse_summer", "ndaic", "ndaic_summer")

# If model comparison is less than 2 set delta aic to 0
model_compglm$ndaic[model_compglm$ndaic < 2] <- 0
model_compglm$ndaic_summer[model_compglm$ndaic_summer < 2] <- 0
model_compglm <- merge(siteinfo, model_compglm, by.x = "sgnum", by.y = "site")

save(model_compglm, file="scripts/users/imyerssmith/shrub_synthesis/model_compglm.RData")
write.csv(model_compglm, file = "scripts/users/imyerssmith/shrub_synthesis/model_compglm.csv", row.names=FALSE)

modelsumglm <- merge(model_compglm, climsum, by.x = "site", by.y = "site", all = FALSE)
modelsumglm <- merge(modelsumglm, soilm)
modelsumglm <- merge(modelsumglm, siteinfo)

# Map figure --------------------------------------------------------------
library(maps)
library(mapdata)
library(mapproj)
library(rgdal)
library(maptools)
library(raster)
library(sp)
library(scales)

# png(file="scripts/users/imyerssmith/shrub_synthesis/crumap.png", width=1000, height=800)
# #par(mfrow=c(1, 2), mar=c(10, 6, 0, 0), oma=c(0, 1, 1, 1), mgp=c(3.5, 1, 0))
# 
# crutif <- "scripts/users/imyerssmith/shrub_synthesis/pre1950-2010.tif"
# crutifr <- raster(crutif)
# test <- projectRaster(crutifr, crs="+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m")
# e <- extent(-3282524, 3272226, -3276109, 3283991)
# test2 <- crop(test, e)
# spplot(test2, col.regions=topo.colors(250))
# sp.points(laea_coord_sp, pch = 3)
# 
# # wrld <- map(plot=FALSE, interior=FALSE, wrap=TRUE, ylim=c(45, 90), xlim=c(-180, 180))
# # wrld_sp <- map2SpatialLines(wrld)
# # proj4string(wrld_sp) <- CRS("+proj=longlat")
# # laea_wrld_sp <- spTransform(wrld_sp, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
# # plot(laea_wrld_sp, col="gray")
# # data <- modelsum
# # coord <- SpatialPoints(data.frame(data$lon, data$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
# # laea_coord_sp <- spTransform(coord, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
# # plot(laea_coord_sp, pch=19, col="#000000", cex=1, add=TRUE)
# 
# dev.off()
# 
# spplot(test)

#load veg data
vegmap <- readOGR("scripts/users/imyerssmith/shrub_synthesis/cp_biozone_la_shp", "cavm_all polygon")
ogrInfo("scripts/users/imyerssmith/shrub_synthesis/cp_biozone_la_shp", "cavm_all polygon")
vegmap2 <- spTransform(vegmap, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m")) 

tiff(file="scripts/users/imyerssmith/shrub_synthesis/Fig1map.tiff", width=2000, height=2300, compression = "none")
par(mfrow=c(1,1), mar=c(0, 0, 0, 0), oma=c(2, 2, 2, 2), mgp=c(5, 2, 0))

#plot map
wrld <- map(plot=FALSE, interior=FALSE, wrap=TRUE, ylim=c(45, 90), xlim=c(-180, 180))
wrld_sp <- map2SpatialLines(wrld)
proj4string(wrld_sp) <- CRS("+proj=longlat")
laea_wrld_sp <- spTransform(wrld_sp, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_wrld_sp, col="black")

#Colour for the veg zones
colour <- vegmap$ZONE
colour[colour==0] = "white"
colour[colour==1] = "grey"
colour[colour==2] = "#F7FE2E60"
colour[colour==3] = "#BFFF0060"
colour[colour==4] = "#80FF0060"
colour[colour==5] = "#31B40460"
colour[colour==6] = "#298A0860"

# plot veg zones
plot(vegmap2, add = TRUE, col=colour, border = FALSE) 

data <- modelsum
data1 <- subset(data, (ndaic_summer > 0) & !is.na(as.numeric(as.character(nest_summer))))
data2 <- subset(data, (ndaic_summer == 0) & ndaic > 0 & (is.na(as.numeric(as.character(data$nest_summer)))))
data3 <- subset(data, (ndaic_summer == 0) & ndaic > 0 & (!is.na(as.numeric(as.character(data$nest_summer)))))
data4 <- subset(data, ndaic == 0 | is.na(as.numeric(as.character(R2MM_best))))

# create color palette from dark blue over pale grey to dark red
#install.packages('ColorBrewer')
require(RColorBrewer)
# my.pal <- brewer.pal(9,"Blues")[-1]
# neutral.middle.color <- rgb(.9,.9,.9)
# my.pal <- c(my.pal[length(my.pal):1],neutral.middle.color,brewer.pal(9,"Reds")[-1])

my.pal <- c("#0001FF", "#0A00F5", "#1400EB", "#1E00E1", "#2800D7", "#3300CD", "#3D00C3", "#4700BA", "#5100B0", "#5B00A6", "#66009C", "#700092", "#7A0088", "#84007F", "#8E0075", "#99006B", "#A30061", "#AD0057", "#B7004D", "#C10044", "#CC003A", "#D60030", "#E00026", "#EA001C", "#F40012", "#FE0009")

scale.slopes <- function(slope.vals,col.palette,center=F){
  min.slope.val <- min(data1$nest_best,na.rm=T)
  max.slope.val <- max(data1$nest_best,na.rm=T)
  if (center){
    #do you want to make sure that 0 values get the middle color?
    min.slope.val <- -1*max(abs(c(min.slope.val,max.slope.val)))
  }
  scaled.slopes <- (slope.vals - min.slope.val)/(max.slope.val- min.slope.val)
  scaled.slopes <- scaled.slopes*(length(col.palette)-1)+1
}

slope.col.indices <- scale.slopes(slope.vals=data1$nest_best, col.palette=my.pal,center=T)

coord1 <- SpatialPoints(data.frame(data1$lon, data1$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord1_sp <- spTransform(coord1, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord1_sp, add = TRUE, pch=1, col="#FFFFFF", cex=(rescale(data1$ndaic, to=c(10,32))), lwd=10)

coord2 <- SpatialPoints(data.frame(data2$lon, data2$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord2_sp <- spTransform(coord2, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord2_sp, add = TRUE, pch=1, col="#FFFFFF", cex=(rescale(data2$ndaic, to=c(10,32))), lwd=10)

coord3 <- SpatialPoints(data.frame(data3$lon, data3$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord3_sp <- spTransform(coord3, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord3_sp, add = TRUE, pch=1, col="#FFFFFF", cex=10, lwd=10)

plot(laea_coord3_sp, add = TRUE, pch=1, col="#000000", cex=10, lwd=10)

plot(laea_coord2_sp, add = TRUE, pch=1, col="#000000", cex=(rescale(data2$ndaic, to=c(10,32))), lwd=10)

plot(laea_coord1_sp, add = TRUE, pch=1, col=my.pal[slope.col.indices], cex=(rescale(data1$ndaic, to=c(10,32))), lwd=10)

coord4 <- SpatialPoints(data.frame(data4$lon, data4$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord4_sp <- spTransform(coord4, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord4_sp, add = TRUE, pch=4, col="#000000", cex=5, lwd=10)

leg.txt <- c("No clim. sens.", "Neg. summer temp. sens.", "Summer temp. sens. model with slope near zero", "Pos. summer temp. sens.", "Other best climate model")
leg.blank <- c("", "", "", "", "")
legend("bottomleft", legend=leg.blank, x.intersp=1.5, y.intersp=2, col=c("#000000", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"), pch=c(4, 1, 1, 1, 1), pt.cex=c(5, 10, 10, 10, 10), bty="n", cex = 3, pt.lwd=10)
legend("bottomleft", legend=leg.txt, x.intersp=1.5, y.intersp=2, col="#FFFFFF00", pt.cex=10, bty="n", cex = 3)
legend("bottomleft", legend=leg.blank, x.intersp=1.5, y.intersp=2, col=c("#FFFFFF00", "#0001FF", "#84007F", "#FE0009", "#000000"), pch=c(1,1,1,1,1), pt.cex=10, bty="n", cex = 3, pt.lwd=10, lty=NULL)

leg.txt <- c(expression(paste("Low sensitivity (", Delta, "AIC",  " < 5)")), expression(paste("High sensitivity (", Delta, "AIC",  " = 40)")))
legend("bottomright", legend=leg.blank, x.intersp=3, y.intersp=3, pch=c(1, 1), pt.cex=c(10, 28), bty="n", col="#FFFFFF00", cex = 3)
legend("bottomright", legend=leg.txt, x.intersp=3, y.intersp=3, pch=c(1, 1), pt.cex=c(10, 28), bty="n", col="#000000", pt.lwd=10, cex = 3)

leg.txt <- c("Bioclimatic", "Zones:", "Zone A", "Zone B", "Zone C", "Zone D", "Zone E")
legend("topleft", legend=leg.txt, x.intersp=1.5, y.intersp=1.5, pch = 15, col = c("#00000000", "#00000000", "#F7FE2E60", "#BFFF0060", "#80FF0060", "#31B40460", "#298A0860"), pt.cex = 10, border = FALSE, bty="n", cex = 3)

dev.off()

# Map pannel figure -------------------------------------------------------
png(file="scripts/users/imyerssmith/shrub_synthesis/FigS12map_allindices.png", width=1000, height=800)
par(mfrow=c(3, 2), mar=c(0, 0, 0, 0), oma=c(2, 2, 2, 2), mgp=c(5, 2, 0))

#plot map - delta AIC
wrld <- map(plot=FALSE, interior=FALSE, wrap=TRUE, ylim=c(45, 90), xlim=c(-180, 180))
wrld_sp <- map2SpatialLines(wrld)
proj4string(wrld_sp) <- CRS("+proj=longlat")
laea_wrld_sp <- spTransform(wrld_sp, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_wrld_sp, col="grey")

data <- modelsum
data1 <- subset(data, (ndaic_summer > 0) & !is.na(as.numeric(as.character(nest_summer))))
data2 <- subset(data, (ndaic_summer == 0) & ndaic > 0 & (is.na(as.numeric(as.character(data$nest_summer)))))
data3 <- subset(data, (ndaic_summer == 0) & ndaic > 0 & (!is.na(as.numeric(as.character(data$nest_summer)))))
data4 <- subset(data, ndaic == 0 | is.na(as.numeric(as.character(R2MM_best))))

coord1 <- SpatialPoints(data.frame(data1$lon, data1$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord_sp1 <- spTransform(coord1, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord_sp1, add = TRUE, pch=1, col="#FFFFFF", cex=(rescale(data1$ndaic, to=c(3,10))), lwd=3)

coord2 <- SpatialPoints(data.frame(data2$lon, data2$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord_sp2 <- spTransform(coord2, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord_sp2, add = TRUE, pch=1, col="#FFFFFF", cex=(rescale(data2$ndaic, to=c(3,10))), lwd=3)

coord3 <- SpatialPoints(data.frame(data3$lon, data3$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord_sp3 <- spTransform(coord3, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord_sp3, add = TRUE, pch=1, col="#FFFFFF", cex=3, lwd=3)

plot(laea_coord_sp3, add = TRUE, pch=1, col="#000000", cex=3, lwd=3)

plot(laea_coord_sp2, add = TRUE, pch=1, col="#000000", cex=(rescale(data2$ndaic, to=c(3,10))), lwd=3)

plot(laea_coord_sp1, add = TRUE, pch=1, col=my.pal[slope.col.indices], cex=(rescale(data1$ndaic, to=c(3,10))), lwd=3)

coord4 <- SpatialPoints(data.frame(data4$lon, data4$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord_sp4 <- spTransform(coord4, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord_sp4, add = TRUE, pch=4, col="#000000", cex=2, lwd=3)

leg.txt <- expression(paste(Delta, "AIC"))
legend("topleft", legend=leg.txt, bty="n", cex = 2)

#plot map - R2
wrld <- map(plot=FALSE, interior=FALSE, wrap=TRUE, ylim=c(45, 90), xlim=c(-180, 180))
wrld_sp <- map2SpatialLines(wrld)
proj4string(wrld_sp) <- CRS("+proj=longlat")
laea_wrld_sp <- spTransform(wrld_sp, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_wrld_sp, col="gray")

coord1 <- SpatialPoints(data.frame(data1$lon, data1$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord_sp1 <- spTransform(coord1, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord_sp1, add = TRUE, pch=1, col="#FFFFFF", cex=(rescale(data1$R2MM_best, to=c(3,10))), lwd=3)

coord2 <- SpatialPoints(data.frame(data2$lon, data2$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord_sp2 <- spTransform(coord2, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord_sp2, add = TRUE, pch=1, col="#FFFFFF", cex=(rescale(data2$R2MM_best, to=c(3,10))), lwd=3)

coord3 <- SpatialPoints(data.frame(data3$lon, data3$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord_sp3 <- spTransform(coord3, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord_sp3, add = TRUE, pch=1, col="#FFFFFF", cex=3, lwd=3)

plot(laea_coord_sp3, add = TRUE, pch=1, col="#000000", cex=3, lwd=3)

plot(laea_coord_sp2, add = TRUE, pch=1, col="#000000", cex=(rescale(data2$R2MM_best, to=c(3,10))), lwd=3)

plot(laea_coord_sp1, add = TRUE, pch=1, col=my.pal[slope.col.indices], cex=(rescale(data1$R2MM_best, to=c(3,10))), lwd=3)

coord4 <- SpatialPoints(data.frame(data4$lon, data4$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord_sp4 <- spTransform(coord4, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord_sp4, add = TRUE, pch=4, col="#000000", cex=2, lwd=3)

leg.txt <- expression(paste("R"^"2"))
legend("topleft", legend=leg.txt, bty="n", cex = 2)

#plot map - slope
wrld <- map(plot=FALSE, interior=FALSE, wrap=TRUE, ylim=c(45, 90), xlim=c(-180, 180))
wrld_sp <- map2SpatialLines(wrld)
proj4string(wrld_sp) <- CRS("+proj=longlat")
laea_wrld_sp <- spTransform(wrld_sp, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_wrld_sp, col="gray")

coord1 <- SpatialPoints(data.frame(data1$lon, data1$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord_sp1 <- spTransform(coord1, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord_sp1, add = TRUE, pch=1, col="#FFFFFF", cex=(rescale(data1$nest_best, to=c(3,10))), lwd=3)

coord2 <- SpatialPoints(data.frame(data2$lon, data2$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord_sp2 <- spTransform(coord2, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord_sp2, add = TRUE, pch=1, col="#FFFFFF", cex=(rescale(data2$nest_best, to=c(3,10))), lwd=3)

coord3 <- SpatialPoints(data.frame(data3$lon, data3$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord_sp3 <- spTransform(coord3, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord_sp3, add = TRUE, pch=1, col="#FFFFFF", cex=3, lwd=3)

plot(laea_coord_sp3, add = TRUE, pch=1, col="#000000", cex=3, lwd=3)

plot(laea_coord_sp2, add = TRUE, pch=1, col="#000000", cex=(rescale(data2$nest_best, to=c(3,10))), lwd=3)

plot(laea_coord_sp1, add = TRUE, pch=1, col=my.pal[slope.col.indices], cex=(rescale(data1$nest_best, to=c(3,10))), lwd=3)

coord4 <- SpatialPoints(data.frame(data4$lon, data4$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord_sp4 <- spTransform(coord4, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord_sp4, add = TRUE, pch=4, col="#000000", cex=2, lwd=3)

leg.txt <- "Slope"
legend("topleft", legend=leg.txt, bty="n", cex = 2)

#plot map - prop
wrld <- map(plot=FALSE, interior=FALSE, wrap=TRUE, ylim=c(45, 90), xlim=c(-180, 180))
wrld_sp <- map2SpatialLines(wrld)
proj4string(wrld_sp) <- CRS("+proj=longlat")
laea_wrld_sp <- spTransform(wrld_sp, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_wrld_sp, col="gray")

coord1 <- SpatialPoints(data.frame(data1$lon, data1$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord_sp1 <- spTransform(coord1, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord_sp1, add = TRUE, pch=1, col="#FFFFFF", cex=(rescale(data1$nprop, to=c(3,10))), lwd=3)

coord2 <- SpatialPoints(data.frame(data2$lon, data2$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord_sp2 <- spTransform(coord2, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord_sp2, add = TRUE, pch=1, col="#FFFFFF", cex=(rescale(data2$nprop, to=c(3,10))), lwd=3)

coord3 <- SpatialPoints(data.frame(data3$lon, data3$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord_sp3 <- spTransform(coord3, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord_sp3, add = TRUE, pch=1, col="#FFFFFF", cex=3, lwd=3)

plot(laea_coord_sp3, add = TRUE, pch=1, col="#000000", cex=3, lwd=3)

plot(laea_coord_sp2, add = TRUE, pch=1, col="#000000", cex=(rescale(data2$nprop, to=c(3,10))), lwd=3)

plot(laea_coord_sp1, add = TRUE, pch=1, col=my.pal[slope.col.indices], cex=(rescale(data1$nprop, to=c(3,10))), lwd=3)

coord4 <- SpatialPoints(data.frame(data4$lon, data4$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord_sp4 <- spTransform(coord4, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord_sp4, add = TRUE, pch=4, col="#000000", cex=2, lwd=3)

leg.txt <- "Prop. sens."
legend("topleft", legend=leg.txt, bty="n", cex = 2)

# legend

plot.new()

leg.txt <- c("No clim. sens.", "Neg. summer temp. sens.", "Pos. summer temp. sens.", "Other best climate model")
legend("bottomleft", legend=leg.txt, x.intersp=1.5, y.intersp=2, col=c("#000000", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00", "#FFFFFF00"), pch=c(4,1,1,1,1), pt.cex=c(2,5,5,5,5), bty="n", cex = 2, pt.lwd=3)
legend("bottomleft", legend=leg.txt, x.intersp=1.5, y.intersp=2, col=c("#FFFFFF00", "#0001FF", "#84007F", "#FE0009", "#000000"), pch=c(0,1,1,1,1), pt.cex=c(2,5,5,5,5), bty="n", pt.lwd=3, cex = 2)

plot.new()

leg.txt <- c("Low sensitivity", "High sensitivity")
legend("bottomright", legend=leg.txt, x.intersp=2.5, y.intersp=3, pch=c(1, 1), pt.cex=c(3, 10), bty="n", cex = 2, col="#000000", pt.lwd=3)

dev.off()

# Models figure -----------------------------------------------------------
models <- read.csv("scripts/users/imyerssmith/shrub_synthesis/modelsn.csv", 
                     header=TRUE, stringsAsFactors=FALSE)

modelspos <- subset(models, estimate >= 0)
modelsneg <- subset(models, estimate < 0)
modelsnull <- subset(models, model == "NULL")

modelspos <- rbind.data.frame(modelspos, modelsnull)
modelsneg <- rbind.data.frame(modelsneg, modelsnull)

# positive models
length <- length(unique(models$sgnum)) # Number of sites
i <- 1
j <- 1
k <- 1
l <- 1
deltaaic <- as.data.frame(array(0, c(33*length, 3)))
colnames(deltaaic) <- c("sgnum", "model", "daic")

for (i in 1:length) { 
  m <- modelspos
  m <- subset(m, sgnum == i)
  j <- 1
  len <- length(unique(m$model))-1
  l <- k + len - 1
  daic <- as.data.frame(array(0, c(len, 2)))
  colnames(daic) <- c("model", "daic")
  for (j in 1:len) {
    daic[j, 1] <- m$model[j]
    aic <- m[which(m$model == "NULL"),3] - m[j, 3]
    if(aic < 2) { aic <- NA }
    daic[j, 2] <- aic
    j <- j+1
  }
  deltaaic[k:l, 1] <- unique(m$sgnum)
  deltaaic[k:l, 2:3] <- daic
  i <- i+1
  k <- k+len
  l <- l+len
}

deltaaicpos2 <- na.omit(subset(deltaaic, sgnum != 0))

deltaaicpos2$variables <- deltaaicpos2$model; deltaaicpos2$modnum <- deltaaicpos2$model; deltaaicpos2$modtype <- deltaaicpos2$model

deltaaicpos2$variables[deltaaicpos2$variables == "njjtmp" ] <- "Jun.-Jul. mean temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "njjtmp" ] <- 3.5; deltaaicpos2$modtype[deltaaicpos2$modtype == "njjtmp" ] <- "tmp"
deltaaicpos2$variables[deltaaicpos2$variables == "njunetmp" ] <- "Jun. mean temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "njunetmp" ] <- 3; deltaaicpos2$modtype[deltaaicpos2$modtype == "njunetmp" ] <- "tmp"
deltaaicpos2$variables[deltaaicpos2$variables == "njulytmp" ] <- "Jul. mean temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "njulytmp" ] <- 4; deltaaicpos2$modtype[deltaaicpos2$modtype == "njulytmp" ] <- "tmp"
deltaaicpos2$variables[deltaaicpos2$variables == "njjtmx" ] <- "Jun.-Jul. max. temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "njjtmx" ] <- 3.5; deltaaicpos2$modtype[deltaaicpos2$modtype == "njjtmx" ] <- "tmx"
deltaaicpos2$variables[deltaaicpos2$variables == "njunetmx" ] <- "Jun. max. temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "njunetmx" ] <- 3; deltaaicpos2$modtype[deltaaicpos2$modtype == "njunetmx" ] <- "tmx"
deltaaicpos2$variables[deltaaicpos2$variables == "njulytmx" ] <- "Jul. max. temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "njulytmx" ] <- 4; deltaaicpos2$modtype[deltaaicpos2$modtype == "njulytmx" ] <- "tmx"
deltaaicpos2$variables[deltaaicpos2$variables == "njjtmn" ] <- "Jun.-Jul. min. temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "njjtmn" ] <- 3.5; deltaaicpos2$modtype[deltaaicpos2$modtype == "njjtmn" ] <- "tmn"
deltaaicpos2$variables[deltaaicpos2$variables == "njunetmn" ] <- "Jun. min. temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "njunetmn" ] <- 3; deltaaicpos2$modtype[deltaaicpos2$modtype == "njunetmn" ] <- "tmn"
deltaaicpos2$variables[deltaaicpos2$variables == "njulytmn" ] <- "Jul. min. temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "njulytmn" ] <- 4; deltaaicpos2$modtype[deltaaicpos2$modtype == "njulytmn" ] <- "tmn"
deltaaicpos2$variables[deltaaicpos2$variables == "nastmp" ] <- "Aug.-Sep. mean temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "nastmp" ] <- 5; deltaaicpos2$modtype[deltaaicpos2$modtype == "nastmp" ] <- "tmp"
deltaaicpos2$variables[deltaaicpos2$variables == "naugtmp" ] <- "Aug. mean temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "naugtmp" ] <- NA; deltaaicpos2$modtype[deltaaicpos2$modtype == "naugtmp" ] <- "tmp"
deltaaicpos2$variables[deltaaicpos2$variables == "nastmx" ] <- "Aug.-Sep. max. temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "nastmx" ] <- 5; deltaaicpos2$modtype[deltaaicpos2$modtype == "nastmx" ] <- "tmx"
deltaaicpos2$variables[deltaaicpos2$variables == "naugtmx" ] <- "Aug. max. temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "naugtmx" ] <- NA; deltaaicpos2$modtype[deltaaicpos2$modtype == "naugtmx" ] <- "tmx"
deltaaicpos2$variables[deltaaicpos2$variables == "nastmn" ] <- "Aug.-Sep. min. temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "nastmn" ] <- 5; deltaaicpos2$modtype[deltaaicpos2$modtype == "nastmn" ] <- "tmn"
deltaaicpos2$variables[deltaaicpos2$variables == "naugtmn" ] <- "Aug. min. temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "naugtmn" ] <- NA; deltaaicpos2$modtype[deltaaicpos2$modtype == "naugtmn" ] <- "tmn"
deltaaicpos2$variables[deltaaicpos2$variables == "nwintmp" ] <- "Winter mean temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "nwintmp" ] <- 1; deltaaicpos2$modtype[deltaaicpos2$modtype == "nwintmp" ] <- "tmp"
deltaaicpos2$variables[deltaaicpos2$variables == "nwintmn" ] <- "Winter min. temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "nwintmn" ] <- 1; deltaaicpos2$modtype[deltaaicpos2$modtype == "nwintmn" ] <- "tmn"
deltaaicpos2$variables[deltaaicpos2$variables == "nwintmx" ] <- "Winter max. temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "nwintmx" ] <- 1; deltaaicpos2$modtype[deltaaicpos2$modtype == "nwintmx" ] <- "tmx"
deltaaicpos2$variables[deltaaicpos2$variables == "nmamtmp" ] <- "Spring mean temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "nmamtmp" ] <- 2; deltaaicpos2$modtype[deltaaicpos2$modtype == "nmamtmp" ] <- "tmp"
deltaaicpos2$variables[deltaaicpos2$variables == "nmamtmn" ] <- "Spring min. temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "nmamtmn" ] <- 2; deltaaicpos2$modtype[deltaaicpos2$modtype == "nmamtmn" ] <- "tmn"
deltaaicpos2$variables[deltaaicpos2$variables == "nmamtmx" ] <- "Spring max. temp."; deltaaicpos2$modnum[deltaaicpos2$modnum == "nmamtmx" ] <- 2; deltaaicpos2$modtype[deltaaicpos2$modtype == "nmamtmx" ] <- "tmx"
deltaaicpos2$variables[deltaaicpos2$variables == "njjpre" ] <- "Jun.-Jul. precip."; deltaaicpos2$modnum[deltaaicpos2$modnum == "njjpre" ] <- 3.5; deltaaicpos2$modtype[deltaaicpos2$modtype == "njjpre" ] <- "pre"
deltaaicpos2$variables[deltaaicpos2$variables == "njunepre" ] <- "Jun. precip."; deltaaicpos2$modnum[deltaaicpos2$modnum == "njunepre" ] <- 3; deltaaicpos2$modtype[deltaaicpos2$modtype == "njunepre" ] <- "pre"
deltaaicpos2$variables[deltaaicpos2$variables == "njulypre" ] <- "Jul. precip."; deltaaicpos2$modnum[deltaaicpos2$modnum == "njulypre" ] <- 4; deltaaicpos2$modtype[deltaaicpos2$modtype == "njulypre" ] <- "pre"
deltaaicpos2$variables[deltaaicpos2$variables == "naspre" ] <- "Aug.-Sep. precip."; deltaaicpos2$modnum[deltaaicpos2$modnum == "naspre" ] <- 5; deltaaicpos2$modtype[deltaaicpos2$modtype == "naspre" ] <- "pre"
deltaaicpos2$variables[deltaaicpos2$variables == "naugpre" ] <- "Aug. precip."; deltaaicpos2$modnum[deltaaicpos2$modnum == "naugpre" ] <- NA; deltaaicpos2$modtype[deltaaicpos2$modtype == "naugpre" ] <- "pre"
deltaaicpos2$variables[deltaaicpos2$variables == "nwinpre" ] <- "Winter precip."; deltaaicpos2$modnum[deltaaicpos2$modnum == "nwinpre" ] <- 1; deltaaicpos2$modtype[deltaaicpos2$modtype == "nwinpre" ] <- "pre"
deltaaicpos2$variables[deltaaicpos2$variables == "nmampre" ] <- "Spring precip."; deltaaicpos2$modnum[deltaaicpos2$modnum == "nmampre" ] <- 2; deltaaicpos2$modtype[deltaaicpos2$modtype == "nmampre" ] <- "pre"
deltaaicpos2$variables[deltaaicpos2$variables == "njjtmp - njjpre" ] <- "JJ temp. - precip."; deltaaicpos2$modnum[deltaaicpos2$modnum == "njjtmp - njjpre" ] <- NA; deltaaicpos2$modtype[deltaaicpos2$modtype == "njjtmp - njjpre" ] <- NA
deltaaicpos2$variables[deltaaicpos2$variables == "njjtmp - nwinpre" ] <- "JJ temp. - Win. precip."; deltaaicpos2$modnum[deltaaicpos2$modnum == "njjtmp - nwinpre" ] <- NA; deltaaicpos2$modtype[deltaaicpos2$modtype == "njjtmp - nwinpre" ] <- NA
deltaaicpos2$variables[deltaaicpos2$variables == "njunetmp - njunepre" ] <- "Jun. temp. - precip."; deltaaicpos2$modnum[deltaaicpos2$modnum == "njunetmp - njunepre" ] <- NA; deltaaicpos2$modtype[deltaaicpos2$modtype == "njunetmp - njunepre" ] <- NA
deltaaicpos2$variables[deltaaicpos2$variables == "njunetmp - nwinpre" ] <- "Jun. temp. - Win. precip."; deltaaicpos2$modnum[deltaaicpos2$modnum == "njunetmp - nwinpre" ] <- NA; deltaaicpos2$modtype[deltaaicpos2$modtype == "njunetmp - nwinpre" ] <- NA
deltaaicpos2$variables[deltaaicpos2$variables == "njulytmp - njulypre" ] <- "Jul. temp. - precip."; deltaaicpos2$modnum[deltaaicpos2$modnum == "njulytmp - njulypre" ] <- NA; deltaaicpos2$modtype[deltaaicpos2$modtype == "njulytmp - njulypre" ] <- NA
deltaaicpos2$variables[deltaaicpos2$variables == "njulytmp - nwinpre" ] <- "Jul. temp. - Win. precip."; deltaaicpos2$modnum[deltaaicpos2$modnum == "njulytmp - nwinpre" ] <- NA; deltaaicpos2$modtype[deltaaicpos2$modtype == "njulytmp - nwinpre" ] <- NA
deltaaicpos2$variables[deltaaicpos2$variables == "naugtmp - naugpre" ] <- "Aug. temp. - precip."; deltaaicpos2$modnum[deltaaicpos2$modnum == "naugtmp - naugpre" ] <- NA; deltaaicpos2$modtype[deltaaicpos2$modtype == "naugtmp - naugpre" ] <- NA
deltaaicpos2$variables[deltaaicpos2$variables == "naugtmp - nwinpre" ] <- "Aug. temp. - Win. precip."; deltaaicpos2$modnum[deltaaicpos2$modnum == "naugtmp - nwinpre" ] <- NA; deltaaicpos2$modtype[deltaaicpos2$modtype == "naugtmp - nwinpre" ] <- NA
# deltaaicpos2$variables[deltaaicpos2$variables == "nspring.tdd" ] <- "Summer TDD"; # deltaaicpos2$modnum[deltaaicpos2$modnum == "nspring.tdd" ] <- NA; # deltaaicpos2$modtype[deltaaicpos2$modtype == "nmampre" ] <- NA
# deltaaicpos2$variables[deltaaicpos2$variables == "ngscru.length" ] <- "Growing season length"; # deltaaicpos2$modnum[deltaaicpos2$modnum == "ngscru.length" ] <- NA; # deltaaicpos2$modtype[deltaaicpos2$modtype == "nmampre" ] <- NA

daicdatapos <- subset(deltaaicpos2, deltaaicpos2$variables!="NA")
daicdatapos <- subset(daicdatapos, daicdatapos$model != "nspring.tdd" & daicdatapos$model != "ngscru.length")
sumdaicpos <- aggregate(daicdatapos$daic, by=list(daicdatapos$variables), FUN=sum, NA.rm=TRUE)
names(sumdaicpos) <- c("variables", "sumdaicpos")
sumdaicpos <- sumdaicpos[order(-sumdaicpos$sumdaicpos), ]
daicdatapos <- merge(daicdatapos, sumdaicpos)
daicdatapos <- daicdatapos[order(daicdatapos$sumdaicpos), ]

# negative models
length <- length(unique(models$sgnum)) # Number of sites
i <- 1
j <- 1
k <- 1
l <- 1
deltaaic <- as.data.frame(array(0, c(33*length, 3)))
colnames(deltaaic) <- c("sgnum", "model", "daic")

for (i in 1:length) { 
  m <- modelsneg
  m <- subset(m, sgnum == i)
  j <- 1
  len <- length(unique(m$model))-1
  l <- k + len - 1
  daic <- as.data.frame(array(0, c(len, 2)))
  colnames(daic) <- c("model", "daic")
  for (j in 1:len) {
    daic[j, 1] <- m$model[j]
    aic <- m[which(m$model == "NULL"),3] - m[j, 3]
    if(aic < 2) { aic <- NA }
    daic[j, 2] <- aic
    j <- j+1
  }
  deltaaic[k:l, 1] <- unique(m$sgnum)
  deltaaic[k:l, 2:3] <- daic
  i <- i+1
  k <- k+len
  l <- l+len
}

deltaaicneg2 <- na.omit(deltaaic)

deltaaicneg2$variables <- deltaaicneg2$model; deltaaicneg2$modnum <- deltaaicneg2$model; deltaaicneg2$modtype <- deltaaicneg2$model

deltaaicneg2$variables[deltaaicneg2$variables == "njjtmp" ] <- "Jun.-Jul. mean temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "njjtmp" ] <- 3.5; deltaaicneg2$modtype[deltaaicneg2$modtype == "njjtmp" ] <- "tmp"
deltaaicneg2$variables[deltaaicneg2$variables == "njunetmp" ] <- "Jun. mean temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "njunetmp" ] <- 3; deltaaicneg2$modtype[deltaaicneg2$modtype == "njunetmp" ] <- "tmp"
deltaaicneg2$variables[deltaaicneg2$variables == "njulytmp" ] <- "Jul. mean temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "njulytmp" ] <- 4; deltaaicneg2$modtype[deltaaicneg2$modtype == "njulytmp" ] <- "tmp"
deltaaicneg2$variables[deltaaicneg2$variables == "njjtmx" ] <- "Jun.-Jul. max. temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "njjtmx" ] <- 3.5; deltaaicneg2$modtype[deltaaicneg2$modtype == "njjtmx" ] <- "tmx"
deltaaicneg2$variables[deltaaicneg2$variables == "njunetmx" ] <- "Jun. max. temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "njunetmx" ] <- 3; deltaaicneg2$modtype[deltaaicneg2$modtype == "njunetmx" ] <- "tmx"
deltaaicneg2$variables[deltaaicneg2$variables == "njulytmx" ] <- "Jul. max. temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "njulytmx" ] <- 4; deltaaicneg2$modtype[deltaaicneg2$modtype == "njulytmx" ] <- "tmx"
deltaaicneg2$variables[deltaaicneg2$variables == "njjtmn" ] <- "Jun.-Jul. min. temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "njjtmn" ] <- 3.5; deltaaicneg2$modtype[deltaaicneg2$modtype == "njjtmn" ] <- "tmn"
deltaaicneg2$variables[deltaaicneg2$variables == "njunetmn" ] <- "Jun. min. temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "njunetmn" ] <- 3; deltaaicneg2$modtype[deltaaicneg2$modtype == "njunetmn" ] <- "tmn"
deltaaicneg2$variables[deltaaicneg2$variables == "njulytmn" ] <- "Jul. min. temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "njulytmn" ] <- 4; deltaaicneg2$modtype[deltaaicneg2$modtype == "njulytmn" ] <- "tmn"
deltaaicneg2$variables[deltaaicneg2$variables == "nastmp" ] <- "Aug.-Sep. mean temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "nastmp" ] <- 5; deltaaicneg2$modtype[deltaaicneg2$modtype == "nastmp" ] <- "tmp"
deltaaicneg2$variables[deltaaicneg2$variables == "naugtmp" ] <- "Aug. mean temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "naugtmp" ] <- NA; deltaaicneg2$modtype[deltaaicneg2$modtype == "naugtmp" ] <- "tmp"
deltaaicneg2$variables[deltaaicneg2$variables == "nastmx" ] <- "Aug.-Sep. max. temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "nastmx" ] <- 5; deltaaicneg2$modtype[deltaaicneg2$modtype == "nastmx" ] <- "tmx"
deltaaicneg2$variables[deltaaicneg2$variables == "naugtmx" ] <- "Aug. max. temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "naugtmx" ] <- NA; deltaaicneg2$modtype[deltaaicneg2$modtype == "naugtmx" ] <- "tmx"
deltaaicneg2$variables[deltaaicneg2$variables == "nastmn" ] <- "Aug.-Sep. min. temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "nastmn" ] <- 5; deltaaicneg2$modtype[deltaaicneg2$modtype == "nastmn" ] <- "tmn"
deltaaicneg2$variables[deltaaicneg2$variables == "naugtmn" ] <- "Aug. min. temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "naugtmn" ] <- NA; deltaaicneg2$modtype[deltaaicneg2$modtype == "naugtmn" ] <- "tmn"
deltaaicneg2$variables[deltaaicneg2$variables == "nwintmp" ] <- "Winter mean temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "nwintmp" ] <- 1; deltaaicneg2$modtype[deltaaicneg2$modtype == "nwintmp" ] <- "tmp"
deltaaicneg2$variables[deltaaicneg2$variables == "nwintmn" ] <- "Winter min. temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "nwintmn" ] <- 1; deltaaicneg2$modtype[deltaaicneg2$modtype == "nwintmn" ] <- "tmn"
deltaaicneg2$variables[deltaaicneg2$variables == "nwintmx" ] <- "Winter max. temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "nwintmx" ] <- 1; deltaaicneg2$modtype[deltaaicneg2$modtype == "nwintmx" ] <- "tmx"
deltaaicneg2$variables[deltaaicneg2$variables == "nmamtmp" ] <- "Spring mean temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "nmamtmp" ] <- 2; deltaaicneg2$modtype[deltaaicneg2$modtype == "nmamtmp" ] <- "tmp"
deltaaicneg2$variables[deltaaicneg2$variables == "nmamtmn" ] <- "Spring min. temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "nmamtmn" ] <- 2; deltaaicneg2$modtype[deltaaicneg2$modtype == "nmamtmn" ] <- "tmn"
deltaaicneg2$variables[deltaaicneg2$variables == "nmamtmx" ] <- "Spring max. temp."; deltaaicneg2$modnum[deltaaicneg2$modnum == "nmamtmx" ] <- 2; deltaaicneg2$modtype[deltaaicneg2$modtype == "nmamtmx" ] <- "tmx"
deltaaicneg2$variables[deltaaicneg2$variables == "njjpre" ] <- "Jun.-Jul. precip."; deltaaicneg2$modnum[deltaaicneg2$modnum == "njjpre" ] <- 3.5; deltaaicneg2$modtype[deltaaicneg2$modtype == "njjpre" ] <- "pre"
deltaaicneg2$variables[deltaaicneg2$variables == "njunepre" ] <- "Jun. precip."; deltaaicneg2$modnum[deltaaicneg2$modnum == "njunepre" ] <- 3; deltaaicneg2$modtype[deltaaicneg2$modtype == "njunepre" ] <- "pre"
deltaaicneg2$variables[deltaaicneg2$variables == "njulypre" ] <- "Jul. precip."; deltaaicneg2$modnum[deltaaicneg2$modnum == "njulypre" ] <- 4; deltaaicneg2$modtype[deltaaicneg2$modtype == "njulypre" ] <- "pre"
deltaaicneg2$variables[deltaaicneg2$variables == "naspre" ] <- "Aug.-Sep. precip."; deltaaicneg2$modnum[deltaaicneg2$modnum == "naspre" ] <- 5; deltaaicneg2$modtype[deltaaicneg2$modtype == "naspre" ] <- "pre"
deltaaicneg2$variables[deltaaicneg2$variables == "naugpre" ] <- "Aug. precip."; deltaaicneg2$modnum[deltaaicneg2$modnum == "naugpre" ] <- NA; deltaaicneg2$modtype[deltaaicneg2$modtype == "naugpre" ] <- "pre"
deltaaicneg2$variables[deltaaicneg2$variables == "nwinpre" ] <- "Winter precip."; deltaaicneg2$modnum[deltaaicneg2$modnum == "nwinpre" ] <- 1; deltaaicneg2$modtype[deltaaicneg2$modtype == "nwinpre" ] <- "pre"
deltaaicneg2$variables[deltaaicneg2$variables == "nmampre" ] <- "Spring precip."; deltaaicneg2$modnum[deltaaicneg2$modnum == "nmampre" ] <- 2; deltaaicneg2$modtype[deltaaicneg2$modtype == "nmampre" ] <- "pre"
deltaaicneg2$variables[deltaaicneg2$variables == "njjtmp - njjpre" ] <- "JJ temp. - precip."; deltaaicneg2$modnum[deltaaicneg2$modnum == "njjtmp - njjpre" ] <- NA; deltaaicneg2$modtype[deltaaicneg2$modtype == "njjtmp - njjpre" ] <- NA
deltaaicneg2$variables[deltaaicneg2$variables == "njjtmp - nwinpre" ] <- "JJ temp. - Win. precip."; deltaaicneg2$modnum[deltaaicneg2$modnum == "njjtmp - nwinpre" ] <- NA; deltaaicneg2$modtype[deltaaicneg2$modtype == "njjtmp - nwinpre" ] <- NA
deltaaicneg2$variables[deltaaicneg2$variables == "njunetmp - njunepre" ] <- "Jun. temp. - precip."; deltaaicneg2$modnum[deltaaicneg2$modnum == "njunetmp - njunepre" ] <- NA; deltaaicneg2$modtype[deltaaicneg2$modtype == "njunetmp - njunepre" ] <- NA
deltaaicneg2$variables[deltaaicneg2$variables == "njunetmp - nwinpre" ] <- "Jun. temp. - Win. precip."; deltaaicneg2$modnum[deltaaicneg2$modnum == "njunetmp - nwinpre" ] <- NA; deltaaicneg2$modtype[deltaaicneg2$modtype == "njunetmp - nwinpre" ] <- NA
deltaaicneg2$variables[deltaaicneg2$variables == "njulytmp - njulypre" ] <- "Jul. temp. - precip."; deltaaicneg2$modnum[deltaaicneg2$modnum == "njulytmp - njulypre" ] <- NA; deltaaicneg2$modtype[deltaaicneg2$modtype == "njulytmp - njulypre" ] <- NA
deltaaicneg2$variables[deltaaicneg2$variables == "njulytmp - nwinpre" ] <- "Jul. temp. - Win. precip."; deltaaicneg2$modnum[deltaaicneg2$modnum == "njulytmp - nwinpre" ] <- NA; deltaaicneg2$modtype[deltaaicneg2$modtype == "njulytmp - nwinpre" ] <- NA
deltaaicneg2$variables[deltaaicneg2$variables == "naugtmp - naugpre" ] <- "Aug. temp. - precip."; deltaaicneg2$modnum[deltaaicneg2$modnum == "naugtmp - naugpre" ] <- NA; deltaaicneg2$modtype[deltaaicneg2$modtype == "naugtmp - naugpre" ] <- NA
deltaaicneg2$variables[deltaaicneg2$variables == "naugtmp - nwinpre" ] <- "Aug. temp. - Win. precip."; deltaaicneg2$modnum[deltaaicneg2$modnum == "naugtmp - nwinpre" ] <- NA; deltaaicneg2$modtype[deltaaicneg2$modtype == "naugtmp - nwinpre" ] <- NA
# deltaaicneg2$variables[deltaaicneg2$variables == "nspring.tdd" ] <- "Summer TDD"; # deltaaicneg2$modnum[deltaaicneg2$modnum == "nspring.tdd" ] <- NA; # deltaaicneg2$modtype[deltaaicneg2$modtype == "nmampre" ] <- NA
# deltaaicneg2$variables[deltaaicneg2$variables == "ngscru.length" ] <- "Growing season length"; # deltaaicneg2$modnum[deltaaicneg2$modnum == "ngscru.length" ] <- NA; # deltaaicneg2$modtype[deltaaicneg2$modtype == "nmampre" ] <- NA

daicdataneg <- subset(deltaaicneg2, deltaaicneg2$variables!="NA")
daicdataneg <- subset(daicdataneg, daicdataneg$model != "nspring.tdd" & daicdataneg$model != "ngscru.length")
sumdaicneg <- aggregate(daicdataneg$daic, by=list(daicdataneg$variables), FUN=sum, NA.rm=TRUE)
names(sumdaicneg) <- c("variables", "sumdaicneg")
sumdaicneg <- sumdaicneg[order(-sumdaicneg$sumdaicneg), ]
daicdataneg <- merge(daicdataneg, sumdaicneg)
daicdataneg <- daicdataneg[order(daicdataneg$sumdaicneg), ]

tiff(file="scripts/users/imyerssmith/shrub_synthesis/Fig2models.tiff", width=2400, height=1800, compression = "none")
par(mfrow=c(1, 2), mar=c(20, 12, 0, 0), oma=c(4, 4, 4, 4), mgp=c(8, 2.5, 0))

modelnumspos <- as.data.frame(table(daicdatapos$modnum, daicdatapos$modtype))
names(modelnumspos) <- c("modnum", "modtype", "count")
modelnumsneg <- as.data.frame(table(daicdataneg$modnum, daicdataneg$modtype))
names(modelnumsneg) <- c("modnum", "modtype", "count")
plot(c(((as.numeric(as.character(modelnumspos$count))/46)*100),(-(as.numeric(as.character(modelnumsneg$count))/46)*100)), c(as.numeric(as.character(modelnumspos$modnum)),as.numeric(as.character(modelnumsneg$modnum))), xlab="", ylab="Percent models climate sensitive (%)", cex.axis=4, cex.lab=5, bty="l", pch=19, col="#FFFFFF00", 
     cex=5, xlim=c(0, 7), ylim=c(-25, 45), xaxt="n", xaxs = "i", yaxs = "i", lwd = 3, tck = -0.015) 
text(1:6, -26, srt = 45, adj = 1.2, labels=c("Winter", "Spring", "June", "July", "Autumn", "Winter"), cex = 4, xpd = TRUE)
axis(2, tick = T, labels = F, lwd = 3, tck = -0.015)
abline(h=-25, lwd = 3)
abline(v=0, lwd = 3)

# mean temp. positive
tmp <- subset(modelnumspos, modelnumspos$modtype == "tmp")
tmp$modnum <- as.numeric(as.character(tmp$modnum))
tmp <- rbind(tmp, c("6", "tmp", tmp$count[tmp$modnum == 1]))
tmp <- rbind(tmp, c("7", "tmp", tmp$count[tmp$modnum == 1]))
tmp <- rbind(tmp, c("0", "tmp", tmp$count[tmp$modnum == 1]))
tmp <- tmp[order(tmp$modnum), ]
tmpdata <- ((as.numeric(as.character(tmp$count))/46)*100)
polygon(c(tmpdata, rep(0, length.out = 9))~c(as.numeric(as.character(tmp$modnum)), rev(as.numeric(as.character(tmp$modnum)))), col = "#DF010135", border = NA)

# mean temp. negative
tmp <- subset(modelnumsneg, modelnumsneg$modtype == "tmp")
tmp$modnum <- as.numeric(as.character(tmp$modnum))
tmp <- rbind(tmp, c("6", "tmp", tmp$count[tmp$modnum == 1]))
tmp <- rbind(tmp, c("7", "tmp", tmp$count[tmp$modnum == 1]))
tmp <- rbind(tmp, c("0", "tmp", tmp$count[tmp$modnum == 1]))
tmp <- tmp[order(tmp$modnum), ]
tmpdata <- (-(as.numeric(as.character(tmp$count))/46)*100)
polygon(c(tmpdata, rep(0, length.out = 10))~c(as.numeric(as.character(tmp$modnum)), rev(as.numeric(as.character(tmp$modnum)))), col = "#DF010135", border = NA)

# max. temp. positive
tmx <- subset(modelnumspos, modelnumspos$modtype == "tmx")
tmx$modnum <- as.numeric(as.character(tmx$modnum))
tmx <- rbind(tmx, c("6", "tmx", tmx$count[tmx$modnum == 1]))
tmx <- rbind(tmx, c("7", "tmx", tmx$count[tmx$modnum == 1]))
tmx <- rbind(tmx, c("0", "tmx", tmx$count[tmx$modnum == 1]))
tmx <- tmx[order(tmx$modnum), ]
tmxdata <- ((as.numeric(as.character(tmx$count))/46)*100)
polygon(c(tmxdata, rep(0, length.out = 9))~c(as.numeric(as.character(tmx$modnum)), rev(as.numeric(as.character(tmx$modnum)))), col = "#DF010145", border = NA)

# max. temp. negative
tmx <- subset(modelnumsneg, modelnumsneg$modtype == "tmx")
tmx$modnum <- as.numeric(as.character(tmx$modnum))
tmx <- rbind(tmx, c("6", "tmx", tmx$count[tmx$modnum == 1]))
tmx <- rbind(tmx, c("7", "tmx", tmx$count[tmx$modnum == 1]))
tmx <- rbind(tmx, c("0", "tmx", tmx$count[tmx$modnum == 1]))
tmx <- tmx[order(tmx$modnum), ]
tmxdata <- (-(as.numeric(as.character(tmx$count))/46)*100)
polygon(c(tmxdata, rep(0, length.out = 10))~c(as.numeric(as.character(tmx$modnum)), rev(as.numeric(as.character(tmx$modnum)))), col = "#DF010145", border = NA)

# min. temp. positive
tmn <- subset(modelnumspos, modelnumspos$modtype == "tmn")
tmn$modnum <- as.numeric(as.character(tmn$modnum))
tmn <- rbind(tmn, c("6", "tmn", tmn$count[tmn$modnum == 1]))
tmn <- rbind(tmn, c("7", "tmn", tmn$count[tmn$modnum == 1]))
tmn <- rbind(tmn, c("0", "tmn", tmn$count[tmn$modnum == 1]))
tmn <- tmn[order(tmn$modnum), ]
tmndata <- ((as.numeric(as.character(tmn$count))/46)*100)
polygon(c(tmndata, rep(0, length.out = 9))~c(as.numeric(as.character(tmn$modnum)), rev(as.numeric(as.character(tmn$modnum)))), col = "#DF010125", border = NA)

# min. temp. negative
tmn <- subset(modelnumsneg, modelnumsneg$modtype == "tmn")
tmn$modnum <- as.numeric(as.character(tmn$modnum))
tmn <- rbind(tmn, c("6", "tmn", tmn$count[tmn$modnum == 1]))
tmn <- rbind(tmn, c("7", "tmn", tmn$count[tmn$modnum == 1]))
tmn <- rbind(tmn, c("0", "tmn", tmn$count[tmn$modnum == 1]))
tmn <- tmn[order(tmn$modnum), ]
tmndata <- (-(as.numeric(as.character(tmn$count))/46)*100)
polygon(c(tmndata, rep(0, length.out = 10))~c(as.numeric(as.character(tmn$modnum)), rev(as.numeric(as.character(tmn$modnum)))), col = "#DF010125", border = NA)

# precip. positive
pre <- subset(modelnumspos, modelnumspos$modtype == "pre")
pre$modnum <- as.numeric(as.character(pre$modnum))
pre <- rbind(pre, c("6", "pre", pre$count[pre$modnum == 1]))
pre <- rbind(pre, c("7", "pre", pre$count[pre$modnum == 1]))
pre <- rbind(pre, c("0", "pre", pre$count[pre$modnum == 1]))
pre <- pre[order(pre$modnum), ]
predata <- ((as.numeric(as.character(pre$count))/46)*100)
polygon(c(predata, rep(0, length.out = 9))~c(as.numeric(as.character(pre$modnum)), rev(as.numeric(as.character(pre$modnum)))), col = "#FFFFFF", border = NA)
polygon(c(predata, rep(0, length.out = 9))~c(as.numeric(as.character(pre$modnum)), rev(as.numeric(as.character(pre$modnum)))), col = "#5882FA80", border = NA)
abline(h=0, lwd = 3)

# precip. negative
pre <- subset(modelnumsneg, modelnumsneg$modtype == "pre")
pre$modnum <- as.numeric(as.character(pre$modnum))
pre <- rbind(pre, c("6", "pre", pre$count[pre$modnum == 1]))
pre <- rbind(pre, c("7", "pre", pre$count[pre$modnum == 1]))
pre <- rbind(pre, c("0", "pre", pre$count[pre$modnum == 1]))
pre <- pre[order(pre$modnum), ]
predata <- (-(as.numeric(as.character(pre$count))/46)*100)
polygon(c(predata, rep(0, length.out = 10))~c(as.numeric(as.character(pre$modnum)), rev(as.numeric(as.character(pre$modnum)))), col = "#FFFFFF", border = NA)
polygon(c(predata, rep(0, length.out = 10))~c(as.numeric(as.character(pre$modnum)), rev(as.numeric(as.character(pre$modnum)))), col = "#5882FA80", border = NA)
abline(h=0, lwd = 1)

leg.txt <- c("Positive Slope Models")
legend("topleft", legend=leg.txt, bty="n", cex = 4, inset=-0.1)

leg.txt <- c("Negative Slope Models")
legend("bottomleft", legend=leg.txt, bty="n", cex = 4, inset=-0.1)

plot.new()
leg.txt <- c("Min. Temp.", "Mean Temp.", "Max. Temp.", "Total Precip.")
legend("left", legend=leg.txt, x.intersp=0.5, y.intersp=1.4, fill = c("#DF010125", "#DF010150", "#DF010175", "#5882FA80"), border = c("#DF010125", "#DF010150", "#DF010175", "#5882FA80"), bty="n", cex = 4)

dev.off()

# All models figure -------------------------------------------------------

# All models

length <- length(unique(models$sgnum)) # Number of sites
i <- 1
j <- 1
k <- 1
l <- 33
deltaaic <- as.data.frame(array(0, c(33*length, 3)))
colnames(deltaaic) <- c("sgnum", "model", "daic")

for (i in 1:length) { 
  m <- models
  m <- subset(m, sgnum == i)
  j <- 1
  for (j in 1:33) {
    daic[j, 1] <- m$model[j]
    aic <- m[34, 3] - m[j, 3]
    if(aic < 2) { aic <- NA }
    daic[j, 2] <- aic
    j <- j+1
  }
  deltaaic[k:l, 1] <- unique(m$sgnum)
  deltaaic[k:l, 2:3] <- daic
  i <- i+1
  k <- k+33
  l <- l+33
}

deltaaic2 <- na.omit(deltaaic)

deltaaic2$variables <- deltaaic2$model; deltaaic2$modnum <- deltaaic2$model; deltaaic2$modtype <- deltaaic2$model

deltaaic2$variables[deltaaic2$variables == "njjtmp" ] <- "Jun.-Jul. mean temp."; deltaaic2$modnum[deltaaic2$modnum == "njjtmp" ] <- 3.5; deltaaic2$modtype[deltaaic2$modtype == "njjtmp" ] <- "tmp"
deltaaic2$variables[deltaaic2$variables == "njunetmp" ] <- "Jun. mean temp."; deltaaic2$modnum[deltaaic2$modnum == "njunetmp" ] <- 3; deltaaic2$modtype[deltaaic2$modtype == "njunetmp" ] <- "tmp"
deltaaic2$variables[deltaaic2$variables == "njulytmp" ] <- "Jul. mean temp."; deltaaic2$modnum[deltaaic2$modnum == "njulytmp" ] <- 4; deltaaic2$modtype[deltaaic2$modtype == "njulytmp" ] <- "tmp"
deltaaic2$variables[deltaaic2$variables == "njjtmx" ] <- "Jun.-Jul. max. temp."; deltaaic2$modnum[deltaaic2$modnum == "njjtmx" ] <- 3.5; deltaaic2$modtype[deltaaic2$modtype == "njjtmx" ] <- "tmx"
deltaaic2$variables[deltaaic2$variables == "njunetmx" ] <- "Jun. max. temp."; deltaaic2$modnum[deltaaic2$modnum == "njunetmx" ] <- 3; deltaaic2$modtype[deltaaic2$modtype == "njunetmx" ] <- "tmx"
deltaaic2$variables[deltaaic2$variables == "njulytmx" ] <- "Jul. max. temp."; deltaaic2$modnum[deltaaic2$modnum == "njulytmx" ] <- 4; deltaaic2$modtype[deltaaic2$modtype == "njulytmx" ] <- "tmx"
deltaaic2$variables[deltaaic2$variables == "njjtmn" ] <- "Jun.-Jul. min. temp."; deltaaic2$modnum[deltaaic2$modnum == "njjtmn" ] <- 3.5; deltaaic2$modtype[deltaaic2$modtype == "njjtmn" ] <- "tmn"
deltaaic2$variables[deltaaic2$variables == "njunetmn" ] <- "Jun. min. temp."; deltaaic2$modnum[deltaaic2$modnum == "njunetmn" ] <- 3; deltaaic2$modtype[deltaaic2$modtype == "njunetmn" ] <- "tmn"
deltaaic2$variables[deltaaic2$variables == "njulytmn" ] <- "Jul. min. temp."; deltaaic2$modnum[deltaaic2$modnum == "njulytmn" ] <- 4; deltaaic2$modtype[deltaaic2$modtype == "njulytmn" ] <- "tmn"
deltaaic2$variables[deltaaic2$variables == "nastmp" ] <- "Aug.-Sep. mean temp."; deltaaic2$modnum[deltaaic2$modnum == "nastmp" ] <- 5; deltaaic2$modtype[deltaaic2$modtype == "nastmp" ] <- "tmp"
deltaaic2$variables[deltaaic2$variables == "naugtmp" ] <- "Aug. mean temp."; deltaaic2$modnum[deltaaic2$modnum == "naugtmp" ] <- NA; deltaaic2$modtype[deltaaic2$modtype == "naugtmp" ] <- "tmp"
deltaaic2$variables[deltaaic2$variables == "nastmx" ] <- "Aug.-Sep. max. temp."; deltaaic2$modnum[deltaaic2$modnum == "nastmx" ] <- 5; deltaaic2$modtype[deltaaic2$modtype == "nastmx" ] <- "tmx"
deltaaic2$variables[deltaaic2$variables == "naugtmx" ] <- "Aug. max. temp."; deltaaic2$modnum[deltaaic2$modnum == "naugtmx" ] <- NA; deltaaic2$modtype[deltaaic2$modtype == "naugtmx" ] <- "tmx"
deltaaic2$variables[deltaaic2$variables == "nastmn" ] <- "Aug.-Sep. min. temp."; deltaaic2$modnum[deltaaic2$modnum == "nastmn" ] <- 5; deltaaic2$modtype[deltaaic2$modtype == "nastmn" ] <- "tmn"
deltaaic2$variables[deltaaic2$variables == "naugtmn" ] <- "Aug. min. temp."; deltaaic2$modnum[deltaaic2$modnum == "naugtmn" ] <- NA; deltaaic2$modtype[deltaaic2$modtype == "naugtmn" ] <- "tmn"
deltaaic2$variables[deltaaic2$variables == "nwintmp" ] <- "Winter mean temp."; deltaaic2$modnum[deltaaic2$modnum == "nwintmp" ] <- 1; deltaaic2$modtype[deltaaic2$modtype == "nwintmp" ] <- "tmp"
deltaaic2$variables[deltaaic2$variables == "nwintmn" ] <- "Winter min. temp."; deltaaic2$modnum[deltaaic2$modnum == "nwintmn" ] <- 1; deltaaic2$modtype[deltaaic2$modtype == "nwintmn" ] <- "tmn"
deltaaic2$variables[deltaaic2$variables == "nwintmx" ] <- "Winter max. temp."; deltaaic2$modnum[deltaaic2$modnum == "nwintmx" ] <- 1; deltaaic2$modtype[deltaaic2$modtype == "nwintmx" ] <- "tmx"
deltaaic2$variables[deltaaic2$variables == "nmamtmp" ] <- "Spring mean temp."; deltaaic2$modnum[deltaaic2$modnum == "nmamtmp" ] <- 2; deltaaic2$modtype[deltaaic2$modtype == "nmamtmp" ] <- "tmp"
deltaaic2$variables[deltaaic2$variables == "nmamtmn" ] <- "Spring min. temp."; deltaaic2$modnum[deltaaic2$modnum == "nmamtmn" ] <- 2; deltaaic2$modtype[deltaaic2$modtype == "nmamtmn" ] <- "tmn"
deltaaic2$variables[deltaaic2$variables == "nmamtmx" ] <- "Spring max. temp."; deltaaic2$modnum[deltaaic2$modnum == "nmamtmx" ] <- 2; deltaaic2$modtype[deltaaic2$modtype == "nmamtmx" ] <- "tmx"
deltaaic2$variables[deltaaic2$variables == "njjpre" ] <- "Jun.-Jul. precip."; deltaaic2$modnum[deltaaic2$modnum == "njjpre" ] <- 3.5; deltaaic2$modtype[deltaaic2$modtype == "njjpre" ] <- "pre"
deltaaic2$variables[deltaaic2$variables == "njunepre" ] <- "Jun. precip."; deltaaic2$modnum[deltaaic2$modnum == "njunepre" ] <- 3; deltaaic2$modtype[deltaaic2$modtype == "njunepre" ] <- "pre"
deltaaic2$variables[deltaaic2$variables == "njulypre" ] <- "Jul. precip."; deltaaic2$modnum[deltaaic2$modnum == "njulypre" ] <- 4; deltaaic2$modtype[deltaaic2$modtype == "njulypre" ] <- "pre"
deltaaic2$variables[deltaaic2$variables == "naspre" ] <- "Aug.-Sep. precip."; deltaaic2$modnum[deltaaic2$modnum == "naspre" ] <- 5; deltaaic2$modtype[deltaaic2$modtype == "naspre" ] <- "pre"
deltaaic2$variables[deltaaic2$variables == "naugpre" ] <- "Aug. precip."; deltaaic2$modnum[deltaaic2$modnum == "naugpre" ] <- NA; deltaaic2$modtype[deltaaic2$modtype == "naugpre" ] <- "pre"
deltaaic2$variables[deltaaic2$variables == "nwinpre" ] <- "Winter precip."; deltaaic2$modnum[deltaaic2$modnum == "nwinpre" ] <- 1; deltaaic2$modtype[deltaaic2$modtype == "nwinpre" ] <- "pre"
deltaaic2$variables[deltaaic2$variables == "nmampre" ] <- "Spring precip."; deltaaic2$modnum[deltaaic2$modnum == "nmampre" ] <- 2; deltaaic2$modtype[deltaaic2$modtype == "nmampre" ] <- "pre"
deltaaic2$variables[deltaaic2$variables == "njjtmp - njjpre" ] <- "JJ temp. - precip."; deltaaic2$modnum[deltaaic2$modnum == "njjtmp - njjpre" ] <- NA; deltaaic2$modtype[deltaaic2$modtype == "njjtmp - njjpre" ] <- NA
deltaaic2$variables[deltaaic2$variables == "njjtmp - nwinpre" ] <- "JJ temp. - Win. precip."; deltaaic2$modnum[deltaaic2$modnum == "njjtmp - nwinpre" ] <- NA; deltaaic2$modtype[deltaaic2$modtype == "njjtmp - nwinpre" ] <- NA
deltaaic2$variables[deltaaic2$variables == "njunetmp - njunepre" ] <- "Jun. temp. - precip."; deltaaic2$modnum[deltaaic2$modnum == "njunetmp - njunepre" ] <- NA; deltaaic2$modtype[deltaaic2$modtype == "njunetmp - njunepre" ] <- NA
deltaaic2$variables[deltaaic2$variables == "njunetmp - nwinpre" ] <- "Jun. temp. - Win. precip."; deltaaic2$modnum[deltaaic2$modnum == "njunetmp - nwinpre" ] <- NA; deltaaic2$modtype[deltaaic2$modtype == "njunetmp - nwinpre" ] <- NA
deltaaic2$variables[deltaaic2$variables == "njulytmp - njulypre" ] <- "Jul. temp. - precip."; deltaaic2$modnum[deltaaic2$modnum == "njulytmp - njulypre" ] <- NA; deltaaic2$modtype[deltaaic2$modtype == "njulytmp - njulypre" ] <- NA
deltaaic2$variables[deltaaic2$variables == "njulytmp - nwinpre" ] <- "Jul. temp. - Win. precip."; deltaaic2$modnum[deltaaic2$modnum == "njulytmp - nwinpre" ] <- NA; deltaaic2$modtype[deltaaic2$modtype == "njulytmp - nwinpre" ] <- NA
deltaaic2$variables[deltaaic2$variables == "naugtmp - naugpre" ] <- "Aug. temp. - precip."; deltaaic2$modnum[deltaaic2$modnum == "naugtmp - naugpre" ] <- NA; deltaaic2$modtype[deltaaic2$modtype == "naugtmp - naugpre" ] <- NA
deltaaic2$variables[deltaaic2$variables == "naugtmp - nwinpre" ] <- "Aug. temp. - Win. precip."; deltaaic2$modnum[deltaaic2$modnum == "naugtmp - nwinpre" ] <- NA; deltaaic2$modtype[deltaaic2$modtype == "naugtmp - nwinpre" ] <- NA
# deltaaic2$variables[deltaaic2$variables == "nspring.tdd" ] <- "Summer TDD"; # deltaaic2$modnum[deltaaic2$modnum == "nspring.tdd" ] <- NA; # deltaaic2$modtype[deltaaic2$modtype == "nmampre" ] <- NA
# deltaaic2$variables[deltaaic2$variables == "ngscru.length" ] <- "Growing season length"; # deltaaic2$modnum[deltaaic2$modnum == "ngscru.length" ] <- NA; # deltaaic2$modtype[deltaaic2$modtype == "nmampre" ] <- NA

daicdata <- subset(deltaaic2, deltaaic2$variables!="NA")
daicdata <- subset(daicdata, daicdata$model != "nspring.tdd" & daicdata$model != "ngscru.length")
sumdaic <- aggregate(daicdata$daic, by=list(daicdata$variables), FUN=sum, NA.rm=TRUE)
names(sumdaic) <- c("variables", "sumdaic")
sumdaic <- sumdaic[order(-sumdaic$sumdaic), ]
daicdata <- merge(daicdata, sumdaic)
daicdata <- daicdata[order(daicdata$sumdaic), ]

png(file="scripts/users/imyerssmith/shrub_synthesis/FigS11allmodels.png", width=800, height=1200)
par(mfrow=c(2, 1), mar=c(10, 6, 0, 0), oma=c(0, 1, 1, 1), mgp=c(3.5, 1, 0))

# Normalized data
colours <- c("#DF010180", "#DF010180", "#5F04B480", "#5F04B480", "#DF010180", "#DF010180", "#5F04B480", "#DF010180", "#DF010180", "#DF010180", "#DF010180", "#5F04B480", "#5F04B480", "#5F04B480", "#DF010180", "#DF010180", "#DF010180", "#DF010180", "#5F04B480", "#5F04B480", "#DF010180", "#DF010180", "#DF010180", "#DF010180", "#DF010180", "#5882FA80", "#DF010180", "#5882FA80", "#5882FA80", "#5882FA80", "#5882FA80")

colours2 <- c("#DF010180", "#DF010180", "#5F04B480", "#DF010180", "#DF010180", "#5F04B480", "#DF010180", "#DF010180", "#5F04B480", "#5F04B480", "#5F04B480", "#DF010180", "#DF010180", "#DF010180", "#5F04B480", "#DF010180", "#5F04B480", "#DF010180", "#5F04B480", "#DF010180", "#DF010180", "#DF010180", "#DF010180", "#DF010180", "#DF010180", "#DF010180", "#5882FA80", "#5882FA80", "#5882FA80", "#5882FA80", "#5882FA80")

barplot(sort(((table(daicdata$variables)/33)*100), decreasing = TRUE), ylab="Percent models climate sensitive (%)", col=colours, 
        xlab="", xaxt="n", cex.axis=1.25, cex.lab=1.5, ylim=c(0, 80), border = NA)
abline(h=0, lty = 1, lwd = 2)
text(seq(1, (length(unique(daicdata$variables))*1.2), by=1.2), 0-(80*0.05), srt = 45, adj = 0.9, labels=names(sort(((table(daicdata$variables)/33)*100), decreasing = TRUE)), xpd = TRUE)
leg.txt <- "A. Percent of site-by-genus combinations that are climate sensitive"
legend("topleft", inset=c(-0.02, 0), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 1.5)

daicorder <- factor(daicdata$variables, levels=sumdaic$variables)

stripchart(daicdata$daic~daicorder, vertical = TRUE, ylab=expression(paste("Climate sensitivity (", Delta, " AIC)")), col=colours2, 
        xlab="", cex.axis=1.25, cex.lab=1.5, pch=19, ylim=c(0, max(na.omit(as.numeric(as.character(daicdata$daic))))+20), las=2, frame.plot = FALSE, xaxt = "n")
abline(h=0, lty = 1, lwd = 2)
text(1:length(sumdaic$variables), 0-(max(daicdata$daic)*0.05), srt = 45, adj = 1, labels=sumdaic$variables, xpd = TRUE)
leg.txt <- "B. Climate sensitivity for all site-by-genus combinations"
legend("topleft", inset=c(-0.02, 0), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 1.5)

dev.off()

# Site models figure -------------------------------------------------------
png(file="scripts/users/imyerssmith/shrub_synthesis/FigS10sitemodels.png", width=2000, height=500)
par(mfrow=c(1, 2), mar=c(10, 6, 0, 0), oma=c(0, 1, 1, 1), mgp=c(3.5, 1, 0))

daicorder <- deltaaic[order(deltaaic$sgnum), ]
daicorder <- merge(daicorder, siteinfo, by.x = "sgnum", by.y = "sgnum", all.x = TRUE)

data1 <- subset(daicorder, model == "njjtmp" | model == "njjtmx" | model == "njjtmn" | model == "njunetmp" | model == "njunetmx" | model == "njunetmn" | model == "njulytmp" | model == "njulytmx" | model == "njulytmn")
data2 <- subset(daicorder, model == "nwintmp" | model == "nwintmn" | model == "nwintmx" | model == "nmamtmp" | model == "nmamtmn" | model == "nmamtmx" | model == "nastmp" | model == "naugtmp" | model == "nastmx" | model == "naugtmx" | model == "nastmn" | model == "naugtmn" | model == "njuneprvtmp" | model == "njuneprvtmp" | model == "njulyprvtmp" | model == "naugprvtmp" | model == "njjprvtmp" | model == "nasprvtmp" | model == "nwinprvtmp")
data3 <- subset(daicorder, model == "njjpre" | model == "njunepre" | model == "njulypre" | model == "naspre" | model == "naugpre" | model == "nwinpre" | model == "nmampre" | model == "njjprvpre" | model == "njuneprvpre" | model == "njulyprvpre" | model == "naugprvpre" | model == "nasprvpre" | model == "nwinprvpre" | model == "nmamprvpre")
data4 <- subset(daicorder, model != "njjtmp" & model != "njjtmx" & model != "njjtmn" & model != "njunetmp" & model != "njunetmx" & model != "njunetmn" & model != "njulytmp" & model != "njulytmx" & model != "njulytmn" & model != "njjpre" & model != "njunepre" & model != "njulypre" & model != "naspre" & model != "naugpre" & model != "nwinpre" & model != "nmampre" & model != "njjprvpre" & model != "njuneprvpre" & model != "njulyprvpre" & model != "naugprvpre" & model != "nasprvpre" & model != "nwinprvpre" & model != "nmamprvpre" & model != "nwintmp" & model != "nwintmn" & model != "nwintmx" & model != "nmamtmp" & model != "nmamtmn" & model != "nmamtmx" & model != "nastmp" & model != "naugtmp" & model != "nastmx" & model != "naugtmx" & model != "nastmn" & model != "naugtmn" & model != "njuneprvtmp" & model != "njuneprvtmp" & model != "njulyprvtmp" & model != "naugprvtmp" & model != "njjprvtmp" & model != "nasprvtmp" & model != "nwinprvtmp")

stripchart(daicorder$sgnum~daicorder$daic, vertical = TRUE, ylab=expression(paste("Climate sensitivity (", Delta, " AIC)")), col="#00000000", 
           xlab="", cex.axis=1.25, cex.lab=1.5, pch=19, ylim=c(0, max(na.omit(as.numeric(as.character(daicorder$daic))))+20), xlim=c(0, max(length(unique(daicorder$sgnum)))), las=2, frame.plot = FALSE, xaxt = "n")
sitecode <- unique(paste(as.character(na.omit(daicorder$site_code)), as.character(na.omit(daicorder$genus_code)), sep=':'))
abline(h=0, lty = 1, lwd = 2)
abline(h=2, lty = 2, lwd = 2)
points(data3$sgnum, data3$daic, pch=19, col = "#5882FA80")
points(data2$sgnum, data2$daic, pch=19, col = "#6E6E6E80")
points(data1$sgnum, data1$daic, pch=19, col = "#DF010180")
points(data4$sgnum, data4$daic, pch=19, col = "#5F04B480")
text(1:length(unique(daicorder$sgnum)), 0-(max(na.omit(daicorder$daic))*0.02), srt = 45, adj = 1.2, labels=unique(sitecode), xpd = TRUE)

plot.new()
leg.txt <- c("Summer temp. models", "Other temp. models", "Precip. models", "Other models")
legend("left", legend=leg.txt, x.intersp=1.5, y.intersp=1.4, col=c("#DF010180", "#6E6E6E80", "#5882FA80", "#5F04B480"), pch=19, pt.cex=1, bty="n", cex = 1)

dev.off()

# Model comparison figure -------------------------------------------------
png(file="scripts/users/imyerssmith/shrub_synthesis/FigS9modelcomp.png", width=1600, height=800)
par(mfrow=c(2, 2), mar=c(10, 6, 0, 0), oma=c(0, 1, 1, 1), mgp=c(3.5, 1, 0))

modelsum_s <- read.csv("scripts/users/imyerssmith/shrub_synthesis/standardized/model_compn.csv", 
                       header=TRUE, stringsAsFactors=FALSE)

data1 <- modelsum[order(modelsum$sgnum), ]
data2 <- modelsum_s[order(modelsum_s$sgnum), ]
data3 <- modelsumglm[order(modelsumglm$sgnum), ]

models_s <- read.csv("scripts/users/imyerssmith/shrub_synthesis/standardized/modelsn.csv", 
                     header=TRUE, stringsAsFactors=FALSE)

allmodels <- cbind.data.frame(factor(unique(models$model)), factor(unique(modelsglm$model)), seq(1:length(unique(models$model))))
names(allmodels) <- c("nmodels", "models", "count")
allmodels <- subset(allmodels, allmodels$models!="NULL")

data1a <- cbind.data.frame(data1$sgnum, data1$ndaic, data1$nmodel_best)
names(data1a) <- c("sgnum", "ndaic", "best_model")
data1am <- merge(allmodels, data1a, by.x = "nmodels", by.y = "best_model", all = TRUE)
data1am <- data1am[order(data1am$count), ]
data1am <- subset(data1am, data1am$models!="NULL")

data2a <- cbind.data.frame(data2$sgnum, data2$ndaic, data2$nmodel_best)
names(data2a) <- c("sgnum", "ndaic", "best_model")
data2am <- merge(allmodels, data2a, by.x = "nmodels", by.y = "best_model", all = TRUE)
data2am <- data2am[order(data2am$count), ]
data2am <- subset(data2am, data2am$models!="NULL")

data3a <- cbind.data.frame(data3$sgnum, data3$ndaic, data3$nmodel_best)
names(data3a) <- c("sgnum", "ndaic", "best_model")
data3am <- merge(allmodels, data3a, by.x = "models", by.y = "best_model", all = TRUE)
data3am <- data3am[order(data3am$count), ]
data3am <- subset(data3am, data3am$models!="NULL")

data1a2 <- cbind(rep("norm", length(row.names(data1am))), data1am)
names(data1a2) <- c("analysis", "nmodels", "models", "count", "sgnum", "ndaic")
data2a2 <- cbind(rep("stand", length(row.names(data2am))), data2am)
names(data2a2) <- c("analysis", "nmodels", "models", "count", "sgnum", "ndaic")
data3a2 <- cbind(rep("chron", length(row.names(data3am))), data3am)
names(data3a2) <- c("analysis", "nmodels", "models", "count", "sgnum", "ndaic")
max <- rbind(data1a2, data2a2, data3a2)
max1 <- ddply(max, .(sgnum), transform, maxndaic = max(ndaic))
max1 <- aggregate(max1$maxndaic, by=list(max1$sgnum), FUN=max, NA.rm=FALSE)
names(max1) <- c("sgnum", "maxndaic")

max2 <- na.omit(max)
max2 <- ddply(max2, .(count), transform, maxndaic = max(ndaic))
max2 <- aggregate(max2$maxndaic, by=list(max2$count), FUN=max, NA.rm=FALSE)
names(max2) <- c("count", "maxndaic")

stripchart(data1$sgnum ~ data1$ndaic, vertical = TRUE, ylab=expression(paste("Climate sensitivity (", Delta, " AIC)")), col="#00000000", 
           xlab="", cex.axis=1.25, cex.lab=1.5, pch=19, ylim=c(0, max(na.omit(as.numeric(as.character(data1$ndaic))))+20), xlim=c(0, max(length(unique(data1$sgnum)))), las=2, frame.plot = FALSE, xaxt = "n")
abline(h=0, lty = 1, lwd = 2)
abline(h=2, lty = 2, lwd = 2)
points(max1$sgnum, max1$maxndaic, type="h", col = "#6E6E6E80", lty = 2, lwd = 2)
points(data1$sgnum, data1$ndaic, pch=19, col = "#8904B180", cex = 2)
points(data2$sgnum, data2$ndaic, pch=19, col = "#086A8780", cex = 2)
points(data3$sgnum, data3$ndaic, pch=19, col = "#FF800080", cex = 2)
text(1:length(unique(data1$sgnum)), 0-(max(na.omit(data1$ndaic))*0.05), srt = 45, adj = 1.1, labels=unique(na.omit(sitecode)), xpd = TRUE)

plot.new()
leg.txt <- c("Normalized data", "Standardized data", "Chronology analysis")
legend("left", legend=leg.txt, x.intersp=1.5, y.intersp=1.4, col=c("#8904B180", "#086A8780", "#FF800080"), pch=19, pt.cex=2, bty="n", cex = 1.5)

modelslable <- cbind.data.frame(unique(daicdata$model), unique(daicdata$variables))
colnames(modelslable) <- c("nmodels", "variables")
modelslable <- merge(allmodels, modelslable)
modelslable <- modelslable[order(modelslable$count), ]

stripchart(data1am$ndaic~data1am$nmodels, vertical = TRUE, ylab=expression(paste("Climate sensitivity (", Delta, " AIC)")), col="#00000000", 
           xlab="", cex.axis=1.25, cex.lab=1.5, pch=19, ylim=c(0, max(na.omit(as.numeric(as.character(data1$ndaic))))+20), las=2, frame.plot = FALSE, xaxt = "n")
abline(h=0, lty = 1, lwd = 2)
abline(h=2, lty = 2, lwd = 2)
points(max2$count, max2$maxndaic, type="h", col = "#6E6E6E80", lty = 2, lwd = 2)
points(data1am$count, data1am$ndaic, pch=19, col = "#8904B180", cex = 2)
points(data2am$count, data2am$ndaic, pch=19, col = "#086A8780", cex = 2)
par(new=TRUE)
stripchart(data1am$ndaic~data1am$models, vertical = TRUE, ylab="", col="#00000000", xaxt='n', yaxt='n', 
           xlab="", cex.axis=1.25, cex.lab=1.5, pch=19, ylim=c(0, max(na.omit(as.numeric(as.character(data1$ndaic)))+10)), las=2, frame.plot = FALSE, xaxt = "n")
points(data3am$count, data3am$ndaic, pch=19, col = "#FF800080", cex = 2)
text(1:length(levels(modelslable$variables)), 0-(max(na.omit(data1am$ndaic))*0.05), srt = 45, adj = 1, labels=modelslable$variables, xpd = TRUE)

plot.new()

dev.off()

# Climate sensitivity figure ----------------------------------------------------------
library(geepack)
library(Zelig)
library(doBy)
library(quantreg)
library(betareg)

tiff(file="scripts/users/imyerssmith/shrub_synthesis/Fig3climsens.tiff", width=3000, height=3000, compression = "none")
par(mfrow=c(3, 3), mar=c(20, 12, 0, 0), oma=c(4, 4, 4, 4), mgp=c(8, 4, 0))

# Wet Days
data <- modelsum
x <- as.numeric(as.character(data$swetd))
y <- as.numeric(as.character(data$ndaic))
z <- abs(as.numeric(as.character(data$nest_summer)))
se <- as.numeric(as.character(data$nse_summer))
w <- as.numeric(as.character(data$nprop))
v <- as.numeric(as.character(data$R2MM_best))
v[v<0] <- 0
vxdata <- cbind.data.frame(x,v)
vxdata <- subset(vxdata, v != "NA")
v1 <- vxdata$v
x1 <- vxdata$x

xl <- as.vector(scale(x, center=min(x), scale=(diff(range(x))+(diff(range(x))*0.01))))+0.001
yl <- as.vector(scale(y, center=min(y), scale=diff(range(y))+(diff(range(y))*0.01)))+0.001
zxl <- na.omit(cbind.data.frame(x, z, se))
zxlx <- as.vector(scale(zxl$x, center=min(zxl$x), scale=diff(range(zxl$x))+(diff(range(zxl$x))*0.01)))+0.001
zxlz <- as.vector(scale(zxl$z, center=min(zxl$z), scale=diff(range(zxl$z))+(diff(range(zxl$z))*0.01)))+0.001
zxsel <- as.vector(scale(zxl$se, center=min(zxl$se), scale=diff(range(zxl$se))+(diff(range(zxl$se))*0.01)))+0.001
wl <- as.vector(scale(w, center=min(w), scale=diff(range(w))+(diff(range(w))*0.01)))+0.001
vl <- as.vector(scale(v, center=min(v), scale=diff(range(v))+(diff(range(v))*0.01)))+0.001
xl1 <- as.vector(scale(x1, center=min(x1), scale=(diff(range(x1))+(diff(range(x1))*0.01))))+0.001
vl1 <- as.vector(scale(v1, center=min(v1), scale=diff(range(v1))+(diff(range(v1))*0.01)))+0.001

plot(yl~xl, xlab="Wet Day Frequency (days)", ylab="Climate sensitivity", cex.lab=5.6, bty="l", pch=19, col="#FFFFFF00", cex=6, xlim=c(0, max(xl)), ylim=c(0, max(yl)), axes=FALSE)
axis(1, pos=0, cex.axis=5, lwd = 3, tck = -0.02, padj = 0.5)
axis(2, pos=0, cex.axis=5, lwd = 3, tck = -0.02)
leg.txt <- "A. Wet Days"
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 5)

# daic
yx <- betareg(yl~xl)
curve(predict(yx, data.frame(xl=x), type="resp"), add=TRUE, lwd=5, col="#DF010180")

s <- summary(yx)
p <- round(s$coefficients$mean[8], digits=2)
pdaic <- if(p<0.01) {paste(expression(p[AIC]), " < ", 0.01)} else {paste(expression(p[AIC]), " = ", p)}

ls1<-rq(yl~xl,tau=0.9)
ls2<-rq(yl~xl,tau=0.1)
px <- c(sort(xl), rev(sort(xl)))
py <- c(sort(ls1[6]$fitted.values), sort(rep(0,length(ls1[6]$fitted.values))))
segments(min(xl), min(ls1[6]$fitted.values), max(xl), max(ls1[6]$fitted.values), col = "#DF010150", lty = 2, lwd = 2, cex = 10)
polygon(px, py, col = "#DF010105", border = "#FFFFFF00")

# R2W
vx <- betareg(vl1~xl1)
curve(predict(vx, data.frame(xl1=x), type="resp"), add=TRUE, lwd=5, col="#066A3780")

s <- summary(vx)
p <- round(s$coefficients$mean[8], digits=2)
pr2 <- if(p<0.01) {paste(expression(p[R2]), " < ", 0.01)} else {paste(expression(p[R2]), " = ", p)}

ls1<-rq(vl1~xl1,tau=0.9)
ls2<-rq(vl1~xl1,tau=0.1)
px <- c(sort(xl1), rev(sort(xl1)))
py <- c(sort(ls1[6]$fitted.values), sort(rep(0,length(ls1[6]$fitted.values))))
segments(min(xl1), min(ls1[6]$fitted.values), max(xl1), max(ls1[6]$fitted.values), col = "#066A3750", lty = 2, lwd = 2, cex = 10)
polygon(px, py, col = "#066A3705", border = "#FFFFFF00")

# slope
zx <- betareg(zxlz ~ zxlx, weights=(1/(zxsel+0.1)^2))
curve(predict(zx, data.frame(zxlx=x), type="resp"), add=TRUE, lwd=5, col="#5882FA80")

s <- summary(zx)
p <- round(s$coefficients$mean[8], digits=2)
ps <- if(p<0.01) {paste(expression(p[slope]), " < ", 0.01)} else {paste(expression(p[slope]), " = ", p)}

ls1<-rq(zxlz ~ zxlx,tau=0.9)
ls2<-rq(zxlz ~ zxlx,tau=0.1)
px <- c(sort(zxlx), rev(sort(zxlx)))
py <- c(sort(ls1[6]$fitted.values), sort(rep(0,length(ls1[6]$fitted.values))))
segments(min(zxlx), min(ls1[6]$fitted.values), max(zxlx), max(ls1[6]$fitted.values), col = "#5882FA50", lty = 2, lwd = 2, cex = 10)
polygon(px, py, col = "#5882FA05", border = "#FFFFFF00")

# prop
wx <- betareg(wl~xl)
curve(predict(wx, data.frame(xl=x), type="resp"), add=TRUE, lwd=5, col="#FF800080")

s <- summary(wx)
p <- round(s$coefficients$mean[8], digits=2)
pp <- if(p<0.01) {paste(expression(p[prop.]), " < ", 0.01)} else {paste(expression(p[prop.]), " = ", p)}

ls1<-rq(wl~xl,tau=0.9)
ls2<-rq(wl~xl,tau=0.1)
px <- c(sort(xl), rev(sort(xl)))
py <- c(sort(ls1[6]$fitted.values), sort(rep(0,length(ls1[6]$fitted.values))))
segments(min(xl), min(ls1[6]$fitted.values), max(xl), max(ls1[6]$fitted.values), col = "#FF800050", lty = 2, lwd = 2, cex = 10)
polygon(px, py, col = "#FF800005", border = "#FFFFFF00")

legend("topright", inset=c(0.05, -0.05), legend=c(pdaic,pr2,ps,pp), x.intersp=1.5, y.intersp=1.4, bty="n", cex = 5)

# Soil Moisture
wc <- cbind(data$soilm, data$ndaic, data$ndaic_summer, as.numeric(as.character(data$nest_best)), as.numeric(as.character(data$nse_best)), as.numeric(as.character(data$R2MM_best)), as.numeric(as.character(data$nest_summer)), as.numeric(as.character(data$nse_summer)), as.numeric(as.character(data$R2MM_summer)), data$nprop)
colnames(wc) <- c("soilm", "ndaic", "ndaic_summer", "nest_best", "nse_best", "R2MM_best", "nest_summer", "nse_summer", "R2MM_summer", "nprop")
wc <- as.data.frame(na.omit(wc))

x <- as.numeric(as.character(wc$soilm))
y <- as.numeric(as.character(wc$ndaic))
z <- abs(as.numeric(as.character(wc$nest_summer)))
se <- as.numeric(as.character(wc$nse_summer))
w <- as.numeric(as.character(wc$nprop))
v <- as.numeric(as.character(wc$R2MM_best))
v[v<0] <- 0
vxdata <- cbind.data.frame(x,v)
vxdata <- subset(vxdata, v != "NA")
v1 <- vxdata$v
x1 <- vxdata$x

xl <- as.vector(scale(x, center=min(x), scale=(diff(range(x))+(diff(range(x))*0.01))))+0.001
yl <- as.vector(scale(y, center=min(y), scale=diff(range(y))+(diff(range(y))*0.01)))+0.001
zxl <- na.omit(cbind.data.frame(x, z, se))
zxlx <- as.vector(scale(zxl$x, center=min(zxl$x), scale=diff(range(zxl$x))+(diff(range(zxl$x))*0.01)))+0.001
zxlz <- as.vector(scale(zxl$z, center=min(zxl$z), scale=diff(range(zxl$z))+(diff(range(zxl$z))*0.01)))+0.001
zxsel <- as.vector(scale(zxl$se, center=min(zxl$se), scale=diff(range(zxl$se))+(diff(range(zxl$se))*0.01)))+0.001
wl <- as.vector(scale(w, center=min(w), scale=diff(range(w))+(diff(range(w))*0.01)))+0.001
vl <- as.vector(scale(v, center=min(v), scale=diff(range(v))+(diff(range(v))*0.01)))+0.001
xl1 <- as.vector(scale(x1, center=min(x1), scale=(diff(range(x1))+(diff(range(x1))*0.01))))+0.001
vl1 <- as.vector(scale(v1, center=min(v1), scale=diff(range(v1))+(diff(range(v1))*0.01)))+0.001

plot(yl~xl, xlab="Soil Moisture (%)", ylab="", cex.lab=5.6, bty="l", pch=19, col="#FFFFFF00", 
     cex=6, xlim=c(0, max(xl)), ylim=c(0, max(yl)), yaxt='n', axes=FALSE)
axis(1, pos=0, cex.axis=5, lwd = 3, tck = -0.02, padj = 0.5)
axis(2, pos=0, cex.axis=5, lwd = 3, tck = -0.02)
leg.txt <- "B. Soil Moisture"
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex=5)

# daic
yx <- betareg(yl~xl)
curve(predict(yx, data.frame(xl=x), type="resp"), add=TRUE, lwd=5, col="#DF010180")

s <- summary(yx)
p <- round(s$coefficients$mean[8], digits=2)
pdaic <- if(p<0.01) {paste(expression(p[AIC]), " < ", 0.01)} else {paste(expression(p[AIC]), " = ", p)}

ls1<-rq(yl~xl,tau=0.9)
ls2<-rq(yl~xl,tau=0.1)
px <- c(sort(xl), rev(sort(xl)))
py <- c(sort(ls1[6]$fitted.values), sort(rep(0,length(ls1[6]$fitted.values))))
segments(min(xl), min(ls1[6]$fitted.values), max(xl), max(ls1[6]$fitted.values), col = "#DF010150", lty = 2, lwd = 2, cex = 5)
polygon(px, py, col = "#DF010105", border = "#FFFFFF00")

# R2W
vx <- betareg(vl1~xl1)
curve(predict(vx, data.frame(xl1=x), type="resp"), add=TRUE, lwd=5, col="#066A3780")

s <- summary(vx)
p <- round(s$coefficients$mean[8], digits=2)
pr2 <- if(p<0.01) {paste(expression(p[R2]), " < ", 0.01)} else {paste(expression(p[R2]), " = ", p)}

ls1<-rq(vl1~xl1,tau=0.9)
ls2<-rq(vl1~xl1,tau=0.1)
px <- c(sort(xl1), rev(sort(xl1)))
py <- c(sort(ls1[6]$fitted.values), sort(rep(0,length(ls1[6]$fitted.values))))
segments(min(xl1), min(ls1[6]$fitted.values), max(xl1), max(ls1[6]$fitted.values), col = "#066A3750", lty = 2, lwd = 2, cex = 5)
polygon(px, py, col = "#066A3705", border = "#FFFFFF00")

# slope
zx <- betareg(zxlz ~ zxlx, weights=(1/(zxsel+0.1)^2))
curve(predict(zx, data.frame(zxlx=x), type="resp"), add=TRUE, lwd=5, col="#5882FA80")

s <- summary(zx)
p <- round(s$coefficients$mean[8], digits=2)
ps <- if(p<0.01) {paste(expression(p[slope]), " < ", 0.01)} else {paste(expression(p[slope]), " = ", p)}

ls1<-rq(zxlz ~ zxlx,tau=0.9)
ls2<-rq(zxlz ~ zxlx,tau=0.1)
px <- c(sort(zxlx), rev(sort(zxlx)))
py <- c(sort(ls1[6]$fitted.values), sort(rep(0,length(ls1[6]$fitted.values))))
segments(min(zxlx), min(ls1[6]$fitted.values), max(zxlx), max(ls1[6]$fitted.values), col = "#5882FA50", lty = 2, lwd = 2, cex = 5)
polygon(px, py, col = "#5882FA05", border = "#FFFFFF00")

# prop
wx <- betareg(wl~xl)
curve(predict(wx, data.frame(xl=x), type="resp"), add=TRUE, lwd=5, col="#FF800080")

s <- summary(wx)
p <- round(s$coefficients$mean[8], digits=2)
pp <- if(p<0.01) {paste(expression(p[prop.]), " < ", 0.01)} else {paste(expression(p[prop.]), " = ", p)}

ls1<-rq(wl~xl,tau=0.9)
ls2<-rq(wl~xl,tau=0.1)
px <- c(sort(xl), rev(sort(xl)))
py <- c(sort(ls1[6]$fitted.values), sort(rep(0,length(ls1[6]$fitted.values))))
segments(min(xl), min(ls1[6]$fitted.values), max(xl), max(ls1[6]$fitted.values), col = "#FF800050", lty = 2, lwd = 2, cex = 5)
polygon(px, py, col = "#FF800005", border = "#FFFFFF00")

legend("topright", inset=c(0.05, -0.05), legend=c(pdaic,pr2,ps,pp), x.intersp=1.5, y.intersp=1.4, bty="n", cex=5)

plot.new()

# Range Edge
x <- data$disrelat
y <- as.numeric(as.character(data$ndaic))
z <- abs(as.numeric(as.character(data$nest_summer)))
se <- as.numeric(as.character(data$nse_summer))
w <- as.numeric(as.character(data$nprop))
v <- as.numeric(as.character(data$R2MM_best))
v[v<0] <- 0
vxdata <- cbind.data.frame(x,v)
vxdata <- subset(vxdata, v != "NA")
v1 <- vxdata$v
x1 <- vxdata$x

xl <- as.vector(scale(x, center=min(x), scale=(diff(range(x))+(diff(range(x))*0.01))))+0.001
yl <- as.vector(scale(y, center=min(y), scale=diff(range(y))+(diff(range(y))*0.01)))+0.001
zxl <- na.omit(cbind.data.frame(x, z, se))
zxlx <- as.vector(scale(zxl$x, center=min(zxl$x), scale=diff(range(zxl$x))+(diff(range(zxl$x))*0.01)))+0.001
zxlz <- as.vector(scale(zxl$z, center=min(zxl$z), scale=diff(range(zxl$z))+(diff(range(zxl$z))*0.01)))+0.001
zxsel <- as.vector(scale(zxl$se, center=min(zxl$se), scale=diff(range(zxl$se))+(diff(range(zxl$se))*0.01)))+0.001
wl <- as.vector(scale(w, center=min(w), scale=diff(range(w))+(diff(range(w))*0.01)))+0.001
vl <- as.vector(scale(v, center=min(v), scale=diff(range(v))+(diff(range(v))*0.01)))+0.001
xl1 <- as.vector(scale(x1, center=min(x1), scale=(diff(range(x1))+(diff(range(x1))*0.01))))+0.001
vl1 <- as.vector(scale(v1, center=min(v1), scale=diff(range(v1))+(diff(range(v1))*0.01)))+0.001

plot(yl~xl, xlab="Dist. to N/Ele. Range Edge (Rel. lat.)", ylab="Climate sensitivity", cex.lab=5.6, bty="l", pch=19, col="#FFFFFF00", cex=6, xlim=c(0, max(xl)), ylim=c(0, max(yl)), axes=FALSE)
axis(1, pos=0, cex.axis=5, lwd = 3, tck = -0.02, padj = 0.5)
axis(2, pos=0, cex.axis=5, lwd = 3, tck = -0.02)
leg.txt <- "C. Range edge"
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex=5)

# daic
yx <- betareg(yl~xl)
curve(predict(yx, data.frame(xl=x), type="resp"), add=TRUE, lwd=5, col="#DF010180")

s <- summary(yx)
p <- round(s$coefficients$mean[8], digits=2)
pdaic <- if(p<0.01) {paste(expression(p[AIC]), " < ", 0.01)} else {paste(expression(p[AIC]), " = ", p)}

ls1<-rq(yl~xl,tau=0.9)
ls2<-rq(yl~xl,tau=0.1)
px <- c(sort(xl), rev(sort(xl)))
py <- c(sort(ls1[6]$fitted.values, decreasing = TRUE), sort(rep(0,length(ls1[6]$fitted.values)), decreasing = TRUE))
segments(min(xl), max(ls1[6]$fitted.values), max(xl), min(ls1[6]$fitted.values), col = "#DF010150", lty = 2, lwd = 2, cex = 5)
polygon(px, py, col = "#DF010105", border = "#FFFFFF00")

# R2W
vx <- betareg(vl1~xl1)
curve(predict(vx, data.frame(xl1=x), type="resp"), add=TRUE, lwd=5, col="#066A3780")

s <- summary(vx)
p <- round(s$coefficients$mean[8], digits=2)
pr2 <- if(p<0.01) {paste(expression(p[R2]), " < ", 0.01)} else {paste(expression(p[R2]), " = ", p)}

ls1<-rq(vl1~xl1,tau=0.9)
ls2<-rq(vl1~xl1,tau=0.1)
px <- c(sort(xl1), rev(sort(xl1)))
py <- c(sort(ls1[6]$fitted.values, decreasing = TRUE), sort(rep(0,length(ls1[6]$fitted.values)), decreasing = TRUE))
segments(min(xl1), max(ls1[6]$fitted.values), max(xl1), min(ls1[6]$fitted.values), col = "#066A3750", lty = 2, lwd = 2, cex = 5)
polygon(px, py, col = "#066A3705", border = "#FFFFFF00")

# slope
zx <- betareg(zxlz ~ zxlx, weights=(1/(zxsel+0.1)^2))
curve(predict(zx, data.frame(zxlx=x), type="resp"), add=TRUE, lwd=5, col="#5882FA80")

s <- summary(zx)
p <- round(s$coefficients$mean[8], digits=2)
ps <- if(p<0.01) {paste(expression(p[slope]), " < ", 0.01)} else {paste(expression(p[slope]), " = ", p)}

ls1<-rq(zxlz ~ zxlx,tau=0.9)
ls2<-rq(zxlz ~ zxlx,tau=0.1)
px <- c(sort(zxlx), rev(sort(zxlx)))
py <- c(sort(ls1[6]$fitted.values, decreasing = TRUE), sort(rep(0,length(ls1[6]$fitted.values)), decreasing = TRUE))
segments(min(zxlx), max(ls1[6]$fitted.values), max(zxlx), min(ls1[6]$fitted.values), col = "#5882FA50", lty = 2, lwd = 2, cex = 5)
polygon(px, py, col = "#5882FA05", border = "#FFFFFF00")

# prop
wx <- betareg(wl~xl)
curve(predict(wx, data.frame(xl=x), type="resp"), add=TRUE, lwd=5, col="#FF800080")

s <- summary(wx)
p <- round(s$coefficients$mean[8], digits=2)
pp <- if(p<0.01) {paste(expression(p[prop.]), " < ", 0.01)} else {paste(expression(p[prop.]), " = ", p)}

ls1<-rq(wl~xl,tau=0.9)
ls2<-rq(wl~xl,tau=0.1)
px <- c(sort(xl), rev(sort(xl)))
py <- c(sort(ls1[6]$fitted.values, decreasing = TRUE), sort(rep(0,length(ls1[6]$fitted.values)), decreasing = TRUE))
segments(min(xl), max(ls1[6]$fitted.values), max(xl), min(ls1[6]$fitted.values), col = "#FF800050", lty = 2, lwd = 2, cex = 5)
polygon(px, py, col = "#FF800005", border = "#FFFFFF00")

legend("topright", inset=c(0.05, -0.05), legend=c(pdaic,pr2,ps,pp), x.intersp=1.5, y.intersp=1.4, bty="n", cex=5)

# Canopy
x <- data$canopy
y <- as.numeric(as.character(data$ndaic))
z <- abs(as.numeric(as.character(data$nest_summer)))
se <- as.numeric(as.character(data$nse_summer))
w <- as.numeric(as.character(data$nprop))
v <- as.numeric(as.character(data$R2MM_best))
v[v<0] <- 0
vxdata <- cbind.data.frame(x,v)
vxdata <- subset(vxdata, v != "NA")
v1 <- vxdata$v
x1 <- vxdata$x

xl <- as.vector(scale(x, center=min(x), scale=(diff(range(x))+(diff(range(x))*0.01))))+0.001
yl <- as.vector(scale(y, center=min(y), scale=diff(range(y))+(diff(range(y))*0.01)))+0.001
zxl <- na.omit(cbind.data.frame(x, z, se))
zxlx <- as.vector(scale(zxl$x, center=min(zxl$x), scale=diff(range(zxl$x))+(diff(range(zxl$x))*0.01)))+0.001
zxlz <- as.vector(scale(zxl$z, center=min(zxl$z), scale=diff(range(zxl$z))+(diff(range(zxl$z))*0.01)))+0.001
zxsel <- as.vector(scale(zxl$se, center=min(zxl$se), scale=diff(range(zxl$se))+(diff(range(zxl$se))*0.01)))+0.001
wl <- as.vector(scale(w, center=min(w), scale=diff(range(w))+(diff(range(w))*0.01)))+0.001
vl <- as.vector(scale(v, center=min(v), scale=diff(range(v))+(diff(range(v))*0.01)))+0.001
xl1 <- as.vector(scale(x1, center=min(x1), scale=(diff(range(x1))+(diff(range(x1))*0.01))))+0.001
vl1 <- as.vector(scale(v1, center=min(v1), scale=diff(range(v1))+(diff(range(v1))*0.01)))+0.001

plot(yl~xl, xlab="Canopy Height (cm)", ylab="", cex.lab=5.6, bty="l", pch=19, col="#FFFFFF00", 
     cex=6, xlim=c(0, max(xl)), ylim=c(0, max(yl)), yaxt='n', axes=FALSE)
axis(1, pos=0, cex.axis=5, lwd = 3, tck = -0.02, padj = 0.5)
axis(2, pos=0, cex.axis=5, lwd = 3, tck = -0.02)
leg.txt <- "D. Canopy height"
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex=5)

# daic
yx <- betareg(yl~xl)
curve(predict(yx, data.frame(xl=x), type="resp"), add=TRUE, lwd=5, col="#DF010180")

s <- summary(yx)
p <- round(s$coefficients$mean[8], digits=2)
pdaic <- if(p<0.01) {paste(expression(p[AIC]), " < ", 0.01)} else {paste(expression(p[AIC]), " = ", p)}

ls1<-rq(yl~xl,tau=0.9)
ls2<-rq(yl~xl,tau=0.1)
px <- c(sort(xl), rev(sort(xl)))
py <- c(sort(ls1[6]$fitted.values), sort(rep(0,length(ls1[6]$fitted.values))))
segments(min(xl), min(ls1[6]$fitted.values), max(xl), max(ls1[6]$fitted.values), col = "#DF010150", lty = 2, lwd = 2, cex = 5)
polygon(px, py, col = "#DF010105", border = "#FFFFFF00")

# R2W
vx <- betareg(vl1~xl1)
curve(predict(vx, data.frame(xl1=x), type="resp"), add=TRUE, lwd=5, col="#066A3780")

s <- summary(vx)
p <- round(s$coefficients$mean[8], digits=2)
pr2 <- if(p<0.01) {paste(expression(p[R2]), " < ", 0.01)} else {paste(expression(p[R2]), " = ", p)}

ls1<-rq(vl1~xl1,tau=0.9)
ls2<-rq(vl1~xl1,tau=0.1)
px <- c(sort(xl1), rev(sort(xl1)))
py <- c(sort(ls1[6]$fitted.values), sort(rep(0,length(ls1[6]$fitted.values))))
segments(min(xl1), min(ls1[6]$fitted.values), max(xl1), max(ls1[6]$fitted.values), col = "#066A3750", lty = 2, lwd = 2, cex = 5)
polygon(px, py, col = "#066A3705", border = "#FFFFFF00")

# slope
zx <- betareg(zxlz ~ zxlx, weights=(1/(zxsel+0.1)^2))
curve(predict(zx, data.frame(zxlx=x), type="resp"), add=TRUE, lwd=5, col="#5882FA80")

s <- summary(zx)
p <- round(s$coefficients$mean[8], digits=2)
ps <- if(p<0.01) {paste(expression(p[slope]), " < ", 0.01)} else {paste(expression(p[slope]), " = ", p)}

ls1<-rq(zxlz ~ zxlx,tau=0.9)
ls2<-rq(zxlz ~ zxlx,tau=0.1)
px <- c(sort(zxlx), rev(sort(zxlx)))
py <- c(sort(ls1[6]$fitted.values), sort(rep(0,length(ls1[6]$fitted.values))))
segments(min(zxlx), min(ls1[6]$fitted.values), max(zxlx), max(ls1[6]$fitted.values), col = "#5882FA50", lty = 2, lwd = 2, cex = 5)
polygon(px, py, col = "#5882FA05", border = "#FFFFFF00")

# prop
wx <- betareg(wl~xl)
curve(predict(wx, data.frame(xl=x), type="resp"), add=TRUE, lwd=5, col="#FF800080")

s <- summary(wx)
p <- round(s$coefficients$mean[8], digits=2)
pp <- if(p<0.01) {paste(expression(p[prop.]), " < ", 0.01)} else {paste(expression(p[prop.]), " = ", p)}

ls1<-rq(wl~xl,tau=0.9)
ls2<-rq(wl~xl,tau=0.1)
px <- c(sort(xl), rev(sort(xl)))
py <- c(sort(ls1[6]$fitted.values), sort(rep(0,length(ls1[6]$fitted.values))))
segments(min(xl), min(ls1[6]$fitted.values), max(xl), max(ls1[6]$fitted.values), col = "#FF800050", lty = 2, lwd = 2, cex = 5)
polygon(px, py, col = "#FF800005", border = "#FFFFFF00")

legend("topright", inset=c(0.05, -0.05), legend=c(pdaic,pr2,ps,pp), x.intersp=1.5, y.intersp=1.4, bty="n", cex=5)

plot.new()

leg.txt <- c(expression(paste(Delta, "AIC")), expression(paste("R"^"2")), "Slope", "Prop. sens.")
legend("left", legend=leg.txt, x.intersp=1.5, y.intersp=1.4, col=c("#DF010180", "#066A3780", "#5882FA80", "#FF800080"),  lty = "solid", lwd = 5, seg.len = 0.8, bty="n", cex = 5.6)

# Genera
par(fig=c(0,0.48,0,0.3), mar=c(20, 16, 0, 0), oma=c(4, 4, 4, 4), mgp=c(12, 4, 0), new = TRUE, bty="l")
x <- as.factor(as.character(data$genus_name))
y <- as.numeric(as.character(data$ndaic))
z <- abs(as.numeric(as.character(data$nest_summer)))
se <- as.numeric(as.character(data$nse_summer))
w <- as.numeric(as.character(data$nprop))
v <- as.numeric(as.character(data$R2MM_best))
v[v<0] <- 0
vxdata <- cbind.data.frame(x,v)
vxdata <- subset(vxdata, v != "NA")
v1 <- vxdata$v
x1 <- vxdata$x

yl <- as.vector(scale(y, center=min(y), scale=diff(range(y))+(diff(range(y))*0.01)))+0.001
zxl <- na.omit(cbind.data.frame(x, z, se))
zxlx <- zxl$x
zxlz <- as.vector(scale(zxl$z, center=min(zxl$z), scale=diff(range(zxl$z))+(diff(range(zxl$z))*0.01)))+0.001
zxsel <- as.vector(scale(zxl$se, center=min(zxl$se), scale=diff(range(zxl$se))+(diff(range(zxl$se))*0.01)))+0.001
wl <- as.vector(scale(w, center=min(w), scale=diff(range(w))+(diff(range(w))*0.01)))+0.001
vl <- as.vector(scale(v, center=min(v), scale=diff(range(v))+(diff(range(v))*0.01)))+0.001
vl1 <- as.vector(scale(v1, center=min(v1), scale=diff(range(v1))+(diff(range(v1))*0.01)))+0.001

myl <- aggregate(yl, by=list(x), FUN=mean, na.rm=TRUE)
mylmx <- aggregate(yl, by=list(x), FUN=max, na.rm=TRUE)
mylmn <- aggregate(yl, by=list(x), FUN=min, na.rm=TRUE)
mvl <- aggregate(vl1, by=list(x1), FUN=mean, na.rm=TRUE)
mvlmx <- aggregate(vl1, by=list(x1), FUN=max, na.rm=TRUE)
mvlmn <- aggregate(vl1, by=list(x1), FUN=min, na.rm=TRUE)
mzl <- aggregate(zxlz, by=list(zxlx), FUN=mean, na.rm=TRUE)
mzlmx <- aggregate(zxlz, by=list(zxlx), FUN=max, na.rm=TRUE)
mzlmn <- aggregate(zxlz, by=list(zxlx), FUN=min, na.rm=TRUE)
mwl <- aggregate(wl, by=list(x), FUN=mean, na.rm=TRUE)
mwlmx <- aggregate(wl, by=list(x), FUN=max, na.rm=TRUE)
mwlmn <- aggregate(wl, by=list(x), FUN=min, na.rm=TRUE)

stripchart(c(myl$x,0)~c(1,2,3,4,5,6,7,8,9), ylab="Climate sensitivity", vertical = TRUE, col="#FFFFFF00", xlab="", cex.axis=5, cex.lab=5.6, pch=19, cex=5, xaxt="n", ylim=c(0, 1.1), lwd = 3, bty="l", tck = -0.02)
abline(h=-0.042, lwd = 3)
abline(v=0.69, lwd = 3)
text(1.5:(length(myl$Group.1)+0.5), 0-(max(yl)/10), srt = 45, adj = 1, labels=myl$Group.1, xpd = TRUE, cex=5)
leg.txt <- "E. Genera"
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex=5)

points(myl$x~as.factor(myl$Group.1), pch=19, col="#DF010180", cex=6)
lines(rbind(as.factor(myl$Group.1)[1],as.factor(myl$Group.1)[1]), rbind(mylmn$x[1],mylmx$x[1]), col = "#DF010150", lwd = 3)
lines(rbind(as.factor(myl$Group.1)[2],as.factor(myl$Group.1)[2]), rbind(mylmn$x[2],mylmx$x[2]), col = "#DF010150", lwd = 3)
lines(rbind(as.factor(myl$Group.1)[3],as.factor(myl$Group.1)[3]), rbind(mylmn$x[3],mylmx$x[3]), col = "#DF010150", lwd = 3)
lines(rbind(as.factor(myl$Group.1)[4],as.factor(myl$Group.1)[4]), rbind(mylmn$x[4],mylmx$x[4]), col = "#DF010150", lwd = 3)
lines(rbind(as.factor(myl$Group.1)[5],as.factor(myl$Group.1)[5]), rbind(mylmn$x[5],mylmx$x[5]), col = "#DF010150", lwd = 3)
lines(rbind(as.factor(myl$Group.1)[6],as.factor(myl$Group.1)[6]), rbind(mylmn$x[6],mylmx$x[6]), col = "#DF010150", lwd = 3)
lines(rbind(as.factor(myl$Group.1)[7],as.factor(myl$Group.1)[7]), rbind(mylmn$x[7],mylmx$x[7]), col = "#DF010150", lwd = 3)
lines(rbind(as.factor(myl$Group.1)[8],as.factor(myl$Group.1)[8]), rbind(mylmn$x[8],mylmx$x[8]), col = "#DF010150", lwd = 3)
points(c(mylmn$Group.1,mylmx$Group.1),c(mylmn$x,mylmx$x), pch = "–", col = "#DF010180", cex=6)

points(mvl$x~c(1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1), pch=19, col = "#066A3780", cex=6)
lines(rbind(1.1,1.1), rbind(mvlmn$x[1],mvlmx$x[1]), col = "#066A3780", lwd = 3)
lines(rbind(2.1,2.1), rbind(mvlmn$x[2],mvlmx$x[2]), col = "#066A3780", lwd = 3)
lines(rbind(3.1,3.1), rbind(mvlmn$x[3],mvlmx$x[3]), col = "#066A3780", lwd = 3)
lines(rbind(4.1,4.1), rbind(mvlmn$x[4],mvlmx$x[4]), col = "#066A3780", lwd = 3)
lines(rbind(5.1,5.1), rbind(mvlmn$x[5],mvlmx$x[5]), col = "#066A3780", lwd = 3)
lines(rbind(6.1,6.1), rbind(mvlmn$x[6],mvlmx$x[6]), col = "#066A3780", lwd = 3)
lines(rbind(7.1,7.1), rbind(mvlmn$x[7],mvlmx$x[7]), col = "#066A3780", lwd = 3)
lines(rbind(8.1,8.1), rbind(mvlmn$x[8],mvlmx$x[8]), col = "#066A3780", lwd = 3)
points(c(c(1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1),c(1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1)),c(mvlmn$x,mvlmx$x), pch = "–", col = "#066A3780", cex=6)

points(mzl$x~c(1.25,2.25,3.25,4.25,5.25,6.25,7.25,8.25), pch=19, col = "#5882FA80", cex=6)
lines(rbind(1.25,1.25), rbind(mzlmn$x[1],mzlmx$x[1]), col = "#5882FA80", lwd = 3)
lines(rbind(2.25,2.25), rbind(mzlmn$x[2],mzlmx$x[2]), col = "#5882FA80", lwd = 3)
lines(rbind(3.25,3.25), rbind(mzlmn$x[3],mzlmx$x[3]), col = "#5882FA80", lwd = 3)
lines(rbind(4.25,4.25), rbind(mzlmn$x[4],mzlmx$x[4]), col = "#5882FA80", lwd = 3)
lines(rbind(5.25,5.25), rbind(mzlmn$x[5],mzlmx$x[5]), col = "#5882FA80", lwd = 3)
lines(rbind(6.25,6.25), rbind(mzlmn$x[6],mzlmx$x[6]), col = "#5882FA80", lwd = 3)
lines(rbind(7.25,7.25), rbind(mzlmn$x[7],mzlmx$x[7]), col = "#5882FA80", lwd = 3)
lines(rbind(8.25,8.25), rbind(mzlmn$x[8],mzlmx$x[8]), col = "#5882FA80", lwd = 3)
points(c(c(1.25,2.25,3.25,4.25,5.25,6.25,7.25,8.25),c(1.25,2.25,3.25,4.25,5.25,6.25,7.25,8.25)),c(mzlmn$x,mzlmx$x), pch = "–", col = "#5882FA80", cex=6)

points(mwl$x~c(1.4,2.4,3.4,4.4,5.4,6.4,7.4,8.4), pch=19, col = "#FF800080", cex=6)
lines(rbind(1.4,1.4), rbind(mwlmn$x[1],mwlmx$x[1]), col = "#FF800080", lwd = 3)
lines(rbind(2.4,2.4), rbind(mwlmn$x[2],mwlmx$x[2]), col = "#FF800080", lwd = 3)
lines(rbind(3.4,3.4), rbind(mwlmn$x[3],mwlmx$x[3]), col = "#FF800080", lwd = 3)
lines(rbind(4.4,4.4), rbind(mwlmn$x[4],mwlmx$x[4]), col = "#FF800080", lwd = 3)
lines(rbind(5.4,5.4), rbind(mwlmn$x[5],mwlmx$x[5]), col = "#FF800080", lwd = 3)
lines(rbind(6.4,6.4), rbind(mwlmn$x[6],mwlmx$x[6]), col = "#FF800080", lwd = 3)
lines(rbind(7.4,7.4), rbind(mwlmn$x[7],mwlmx$x[7]), col = "#FF800080", lwd = 3)
lines(rbind(8.4,8.4), rbind(mwlmn$x[8],mwlmx$x[8]), col = "#FF800080", lwd = 3)
points(c(c(1.4,2.4,3.4,4.4,5.4,6.4,7.4,8.4),c(1.4,2.4,3.4,4.4,5.4,6.4,7.4,8.4)),c(mwlmn$x,mwlmx$x), pch = "–", col = "#FF800080", cex=6)

# Data type
par(fig=c(0.45,0.65,0,0.3), mar=c(20, 12, 0, 0), oma=c(4, 4, 4, 4), mgp=c(12, 4, 0), new = TRUE, bty="l")
x <- gsub("increment", "Increment", gsub("ring", "Ring width", data$datatype))
y <- as.numeric(as.character(data$ndaic))
z <- abs(as.numeric(as.character(data$nest_summer)))
se <- as.numeric(as.character(data$nse_summer))
w <- as.numeric(as.character(data$nprop))
v <- as.numeric(as.character(data$R2MM_best))
v[v<0] <- 0
vxdata <- cbind.data.frame(x,v)
vxdata <- subset(vxdata, v != "NA")
v1 <- vxdata$v
x1 <- vxdata$x

yl <- as.vector(scale(y, center=min(y), scale=diff(range(y))+(diff(range(y))*0.01)))+0.001
zxl <- na.omit(cbind.data.frame(x, z, se))
zxlx <- zxl$x
zxlz <- as.vector(scale(zxl$z, center=min(zxl$z), scale=diff(range(zxl$z))+(diff(range(zxl$z))*0.01)))+0.001
zxsel <- as.vector(scale(zxl$se, center=min(zxl$se), scale=diff(range(zxl$se))+(diff(range(zxl$se))*0.01)))+0.001
wl <- as.vector(scale(w, center=min(w), scale=diff(range(w))+(diff(range(w))*0.01)))+0.001
vl <- as.vector(scale(v, center=min(v), scale=diff(range(v))+(diff(range(v))*0.01)))+0.001
xl1 <- as.vector(scale(x1, center=min(x1), scale=(diff(range(x1))+(diff(range(x1))*0.01))))+0.001
vl1 <- as.vector(scale(v1, center=min(v1), scale=diff(range(v1))+(diff(range(v1))*0.01)))+0.001

myl <- aggregate(yl, by=list(x), FUN=mean, na.rm=TRUE)
mylmx <- aggregate(yl, by=list(x), FUN=max, na.rm=TRUE)
mylmn <- aggregate(yl, by=list(x), FUN=min, na.rm=TRUE)
mvl <- aggregate(vl1, by=list(x1), FUN=mean, na.rm=TRUE)
mvlmx <- aggregate(vl1, by=list(x1), FUN=max, na.rm=TRUE)
mvlmn <- aggregate(vl1, by=list(x1), FUN=min, na.rm=TRUE)
mzl <- aggregate(zxlz, by=list(zxlx), FUN=mean, na.rm=TRUE)
mzlmx <- aggregate(zxlz, by=list(zxlx), FUN=max, na.rm=TRUE)
mzlmn <- aggregate(zxlz, by=list(zxlx), FUN=min, na.rm=TRUE)
mwl <- aggregate(wl, by=list(x), FUN=mean, na.rm=TRUE)
mwlmx <- aggregate(wl, by=list(x), FUN=max, na.rm=TRUE)
mwlmn <- aggregate(wl, by=list(x), FUN=min, na.rm=TRUE)

stripchart(myl$x~c(1,2), vertical = TRUE, col="#FFFFFF00", ylab="", xlab="", cex.axis=5, cex.lab=5.6, pch=19, cex=5, xaxt="n", yaxt="n", xlim=c(0.75, 3), ylim=c(0, 1.1), lwd = 3, bty="l", tck = -0.02)
abline(h=-0.042, lwd = 3)
abline(v=0.662, lwd = 3)
text(1.5:2.5, 0-(max(yl)/10), srt = 45, adj = 1, labels=myl$Group.1, xpd = TRUE, cex=5)
leg.txt <- "F. Growth measure"
legend("topleft", inset=c(-0.15, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex=5)

points(myl$x~as.factor(myl$Group.1), pch=19, col = "#DF010180", cex=6)
lines(rbind(as.factor(myl$Group.1)[1],as.factor(myl$Group.1)[1]), rbind(mylmn$x[1],mylmx$x[1]), col = "#DF010150", lwd = 3)
lines(rbind(as.factor(myl$Group.1)[2],as.factor(myl$Group.1)[2]), rbind(mylmn$x[2],mylmx$x[2]), col = "#DF010150", lwd = 3)
points(c(as.factor(mylmn$Group.1),as.factor(mylmx$Group.1)),c(mylmn$x,mylmx$x), pch = "–", col = "#DF010180", cex=6)

points(mvl$x~c(1.1,2.1), pch=19, col = "#066A3780", cex=6)
lines(rbind(1.1,1.1), rbind(mvlmn$x[1],mvlmx$x[1]), col = "#066A3780", lwd = 3)
lines(rbind(2.1,2.1), rbind(mvlmn$x[2],mvlmx$x[2]), col = "#066A3780", lwd = 3)
points(c(c(1.1,2.1),c(1.1,2.1)),c(mvlmn$x,mvlmx$x), pch = "–", col = "#066A3780", cex=6)

points(mzl$x~c(1.25,2.25), pch=19, col = "#5882FA80", cex=6)
lines(rbind(1.25,1.25), rbind(mzlmn$x[1],mzlmx$x[1]), col = "#5882FA80", lwd = 3)
lines(rbind(2.25,2.25), rbind(mzlmn$x[2],mzlmx$x[2]), col = "#5882FA80", lwd = 3)
points(c(c(1.25,2.25),c(1.25,2.25)),c(mzlmn$x,mzlmx$x), pch = "–", col = "#5882FA80", cex=6)

points(mwl$x~c(1.4,2.4), pch=19, col = "#FF800080", cex=6)
lines(rbind(1.4,1.4), rbind(mwlmn$x[1],mwlmx$x[1]), col = "#FF800080", lwd = 3)
lines(rbind(2.4,2.4), rbind(mwlmn$x[2],mwlmx$x[2]), col = "#FF800080", lwd = 3)
points(c(c(1.4,2.4),c(1.4,2.4)),c(mwlmn$x,mwlmx$x), pch = "–", col = "#FF800080", cex=6)

dev.off()

# All indices - Moisture -------------------------------------------------
data <- modelsum
data$mcolour <- data$moist
data$mcolour[data$mcolour == "W"] <- "#5882FA80"
data$mcolour[data$mcolour == "M"] <- "#00800080"
data$mcolour[data$mcolour == "D"] <- "#FF800080"
data$mcolour[is.na(data$mcolour)] <- "#6E6E6E80"

data1 <- subset(modelsum, modelsum$R2MM_best != "NA")
data1$mcolour <- data1$moist
data1$mcolour[data1$mcolour == "W"] <- "#5882FA80"
data1$mcolour[data1$mcolour == "M"] <- "#00800080"
data1$mcolour[data1$mcolour == "D"] <- "#FF800080"
data1$mcolour[is.na(data1$mcolour)] <- "#6E6E6E80"

png(file="scripts/users/imyerssmith/shrub_synthesis/FigS14allindices_moist.png", width=2000, height=2000)
par(mfrow=c(4, 4), mar=c(10, 10, 4, 4), oma=c(2, 2, 2, 2), mgp=c(5, 2, 0))

# Wet Days
x <- as.numeric(as.character(data$swetd))
y <- as.numeric(as.character(data$ndaic))
z <- abs(as.numeric(as.character(data$nest_summer)))
se <- as.numeric(as.character(data$nse_summer))
w <- as.numeric(as.character(data$nprop))
v <- as.numeric(as.character(data$R2MM_best))
v1 <- as.numeric(as.character(data1$R2MM_best))
x1 <- as.numeric(as.character(data1$swetd))

# daic
plot(y~x, xlab="Wet Day Frequency (days)", ylab=expression(paste(Delta, "AIC")), cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data$mcolour, cex=3, xlim=c(0, max(x)), ylim=c(0, max(y)+(max(y)/20)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x, y)
polygon(x[conv.hull], y[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "A. Wet Days"
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# R2W
plot(v1~x1, xlab="Wet Day Frequency (days)", ylab=expression(paste("R"^"2")), cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data1$mcolour, cex=3, xlim=c(0, max(x1)), ylim=c(0, max(v1)+(max(v1)/20)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x1, v1)
polygon(x1[conv.hull], v1[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "B."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# slope
xz <- na.omit(cbind.data.frame(x, z, se))
plot(xz$z~xz$x, xlab="Wet Day Frequency (days)", ylab="Abs. model slopes", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data$mcolour, cex=3, xlim=c(0, max(xz$x)), ylim=c(0, max(xz$z)+(max(xz$z)/10)))
abline(h=0, lty = 2, lwd = 2)

zy <- xz$z
zx <- xz$x
er <- xz$se
g <- lm(zy ~ zx, weights=(1/er^2))
curve <- data.frame(curve(predict(g, data.frame(zx=x), type="resp"), add=TRUE), lwd=2)

p <- round(summary(g)[5][[1]][8], digits=2)
r2 <- round(summary(g)[9][[1]], digits=2)
p2 <- if(p<0.01) {bquote(C.~p~"<"~0.01~", "~R[2]~"="~.(r2))} else {bquote(C.~p~"="~.(p)~", "~R[2]~"="~.(r2))}
leg.txt <- p2
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

i <- 1
CI <- as.data.frame(array(0, c(20, 4)))
colnames(CI) <- c("x", "est", "down", "up")
range <- seq(min(x), max(x), length.out=20)

for (i in 1:20){
  CI[i, 1] <- as.numeric(range[i])
  CI[i, 2] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[2])
  CI[i, 3] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[7])
  CI[i, 4] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[8])
  i <- i+1
}

polygon(c(CI$x, rev(CI$x)), c(CI$up, rev(CI$down)), col = "#84848425", border = NA)

# prop
plot(w~x, xlab="Wet Day Frequency (days)", ylab="Proportion", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data$mcolour, cex=3, xlim=c(0, max(x)), ylim=c(0, max(w)+(max(w)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x, w)
polygon(x[conv.hull], w[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "D."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# Soil moisture
wc <- cbind(data$soilm, data$ndaic, data$ndaic_summer, as.numeric(as.character(data$nest_best)), as.numeric(as.character(data$nse_best)), as.numeric(as.character(data$R2MM_best)), as.numeric(as.character(data$nest_summer)), as.numeric(as.character(data$nse_summer)), as.numeric(as.character(data$R2MM_summer)), as.numeric(as.character(data$nprop)), as.character(data$mcolour))
colnames(wc) <- c("soilm", "ndaic", "ndaic_summer", "nest_best", "nse_best", "R2MM_best", "nest_summer", "nse_summer", "R2MM_summer", "nprop", "mcolour")
wc <- as.data.frame(na.omit(wc))

wc1 <- cbind(data1$soilm, data1$ndaic, data1$ndaic_summer, as.numeric(as.character(data1$nest_best)), as.numeric(as.character(data1$nse_best)), as.numeric(as.character(data1$R2MM_best)), as.numeric(as.character(data1$nest_summer)), as.numeric(as.character(data1$nse_summer)), as.numeric(as.character(data1$R2MM_summer)), as.numeric(as.character(data1$nprop)), as.character(data1$mcolour))
colnames(wc1) <- c("soilm", "ndaic", "ndaic_summer", "nest_best", "nse_best", "R2MM_best", "nest_summer", "nse_summer", "R2MM_summer", "nprop", "mcolour")
wc1 <- as.data.frame(na.omit(wc1))

x <- as.numeric(as.character(wc$soilm))
y <- as.numeric(as.character(wc$ndaic))
z <- abs(as.numeric(as.character(wc$nest_summer)))
se <- as.numeric(as.character(wc$nse_summer))
w <- as.numeric(as.character(wc$nprop))
v <- as.numeric(as.character(wc$R2MM_best))
v1 <- as.numeric(as.character(wc1$R2MM_best))
x1 <- as.numeric(as.character(wc1$soilm))

# daic
plot(y~x, xlab="Soil Moisture (%)", ylab=expression(paste(Delta, "AIC")), cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=as.character(wc$mcolour), cex=3, xlim=c(0, max(x)), ylim=c(0, max(y)+(max(y)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x, y)
polygon(x[conv.hull], y[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "E. Soil Moisture"
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# R2W
plot(v1~x1, xlab="Soil Moisture (%)", ylab=expression(paste("R"^"2")), cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=as.character(wc1$mcolour), cex=3, xlim=c(0, max(x1)), ylim=c(0, max(v1)+(max(v1)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x1, v1)
polygon(x1[conv.hull], v1[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "F."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# slope
xz <- na.omit(cbind.data.frame(x, z, se))
plot(xz$z~xz$x, xlab="Soil Moisture (%)", ylab="Abs. model slopes", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=as.character(wc$mcolour), cex=3, xlim=c(0, max(xz$x)), ylim=c(0, max(xz$z)+(max(xz$z)/10)))
abline(h=0, lty = 2, lwd = 2)

zy <- xz$z
zx <- xz$x
er <- xz$se
g <- lm(zy ~ zx, weights=(1/er^2))
curve <- data.frame(curve(predict(g, data.frame(zx=x), type="resp"), add=TRUE), lwd=2)

p <- round(summary(g)[5][[1]][8], digits=2)
r2 <- round(summary(g)[9][[1]], digits=2)
p2 <- if(p<0.01) {bquote(G.~p~"<"~0.01~", "~R[2]~"="~.(r2))} else {bquote(G.~p~"="~.(p)~", "~R[2]~"="~.(r2))}
leg.txt <- p2
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

i <- 1
CI <- as.data.frame(array(0, c(20, 4)))
colnames(CI) <- c("x", "est", "down", "up")
range <- seq(min(x), max(x), length.out=20)

for (i in 1:20){
  CI[i, 1] <- as.numeric(range[i])
  CI[i, 2] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[2])
  CI[i, 3] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[7])
  CI[i, 4] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[8])
  i <- i+1
}

polygon(c(CI$x, rev(CI$x)), c(CI$up, rev(CI$down)), col = "#84848425", border = NA)

# prop
plot(w~x, xlab="Soil Moisture (%)", ylab="Proportion", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=as.character(wc$mcolour), cex=3, xlim=c(0, max(x)), ylim=c(0, max(w)+(max(w)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x, w)
polygon(x[conv.hull], w[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "H."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# Rain
x <- as.numeric(as.character(data$sjjpre))
y <- as.numeric(as.character(data$ndaic))
z <- abs(as.numeric(as.character(data$nest_summer)))
se <- as.numeric(as.character(data$nse_summer))
w <- as.numeric(as.character(data$nprop))
v <- as.numeric(as.character(data$R2MM_best))
v1 <- as.numeric(as.character(data1$R2MM_best))
x1 <- as.numeric(as.character(data1$sjjpre))

# daic
plot(y~x, xlab="Summer Rain (mm)", ylab=expression(paste(Delta, "AIC")), cex.axis=2.5, cex.lab=2.8, bty="l",  pch=19, col=data$mcolour, cex=3, xlim=c(0, max(x)), ylim=c(0, max(y)+(max(y)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x, y)
polygon(x[conv.hull], y[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "I. Summer Rain"
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# R2W
plot(v1~x1, xlab="Summer Rain (mm)", ylab=expression(paste("R"^"2")), cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data1$mcolour, cex=3, xlim=c(0, max(x1)), ylim=c(0, max(v1)+(max(v1)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x1, v1)
polygon(x1[conv.hull], v1[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "J."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# slope
xz <- na.omit(cbind.data.frame(x, z, se))
plot(xz$z~xz$x, xlab="Summer Rain (mm)", ylab="Abs. model slopes", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data$mcolour, cex=3, xlim=c(0, max(xz$x)), ylim=c(0, max(xz$z)+(max(xz$z)/10)))
abline(h=0, lty = 2, lwd = 2)

zy <- xz$z
zx <- xz$x
er <- xz$se
g <- lm(zy ~ zx, weights=(1/er^2))
curve <- data.frame(curve(predict(g, data.frame(zx=x), type="resp"), add=TRUE), lwd=2)

p <- round(summary(g)[5][[1]][8], digits=2)
r2 <- round(summary(g)[9][[1]], digits=2)
p2 <- if(p<0.01) {bquote(K.~p~"<"~0.01~", "~R[2]~"="~.(r2))} else {bquote(K.~p~"="~.(p)~", "~R[2]~"="~.(r2))}
leg.txt <- p2
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

i <- 1
CI <- as.data.frame(array(0, c(20, 4)))
colnames(CI) <- c("x", "est", "down", "up")
range <- seq(min(x), max(x), length.out=20)

for (i in 1:20){
  CI[i, 1] <- as.numeric(range[i])
  CI[i, 2] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[2])
  CI[i, 3] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[7])
  CI[i, 4] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[8])
  i <- i+1
}

polygon(c(CI$x, rev(CI$x)), c(CI$up, rev(CI$down)), col = "#84848425", border = NA)

# prop
plot(w~x, xlab="Summer Rain (mm)", ylab="Proportion", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data$mcolour, cex=3, xlim=c(0, max(x)), ylim=c(0, max(w)+(max(w)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x, w)
polygon(x[conv.hull], w[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "L."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# Snow
x <- as.numeric(as.character(data$swinpre))
y <- as.numeric(as.character(data$ndaic))
z <- abs(as.numeric(as.character(data$nest_summer)))
se <- as.numeric(as.character(data$nse_summer))
w <- as.numeric(as.character(data$nprop))
v <- as.numeric(as.character(data$R2MM_best))
v1 <- as.numeric(as.character(data1$R2MM_best))
x1 <- as.numeric(as.character(data1$swinpre))

# daic
plot(y~x, xlab="Winter Snow (mm)", ylab=expression(paste(Delta, "AIC")), cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data$mcolour, cex=3, xlim=c(0, max(x)), ylim=c(0, max(y)+(max(y)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x, y)
polygon(x[conv.hull], y[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "M. Winter Snow"
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# R2W
plot(v1~x1, xlab="Winter Snow (mm)", ylab=expression(paste("R"^"2")), cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data1$mcolour, cex=3, xlim=c(0, max(x1)), ylim=c(0, max(v1)+(max(v1)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x1, v1)
polygon(x1[conv.hull], v1[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "N."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# slope
xz <- na.omit(cbind.data.frame(x, z, se))
plot(xz$z~xz$x, xlab="Winter Snow (mm)", ylab="Abs. model slopes", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data$mcolour, cex=3, xlim=c(0, max(xz$x)), ylim=c(0, max(xz$z)+(max(xz$z)/10)))
abline(h=0, lty = 2, lwd = 2)

zy <- xz$z
zx <- xz$x
er <- xz$se
g <- lm(zy ~ zx, weights=(1/er^2))
curve <- data.frame(curve(predict(g, data.frame(zx=x), type="resp"), add=TRUE), lwd=2)

p <- round(summary(g)[5][[1]][8], digits=2)
r2 <- round(summary(g)[9][[1]], digits=2)
p2 <- if(p<0.01) {bquote(O.~p~"<"~0.01~", "~R[2]~"="~.(r2))} else {bquote(O.~p~"="~.(p)~", "~R[2]~"="~.(r2))}
leg.txt <- p2
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

i <- 1
CI <- as.data.frame(array(0, c(20, 4)))
colnames(CI) <- c("x", "est", "down", "up")
range <- seq(min(x), max(x), length.out=20)

for (i in 1:20){
  CI[i, 1] <- as.numeric(range[i])
  CI[i, 2] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[2])
  CI[i, 3] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[7])
  CI[i, 4] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[8])
  i <- i+1
}

polygon(c(CI$x, rev(CI$x)), c(CI$up, rev(CI$down)), col = "#84848425", border = NA)

# prop
plot(w~x, xlab="Winter Snow (mm)", ylab="Proportion", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data$mcolour, cex=3, xlim=c(0, max(x)), ylim=c(0, max(w)+(max(w)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x, w)
polygon(x[conv.hull], w[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "P."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

dev.off()

# All indices - Temperatures ---------------------------------------------
data <- modelsum
data$pcolour <- as.numeric(as.character(data$samt))
data$pcolour[data$pcolour > 0] <- "#DF010180"
data$pcolour[data$pcolour != "#DF010180"] <- "#08088A80"

data1 <- subset(modelsum, modelsum$R2MM_best != "NA")
data1$pcolour <- as.numeric(as.character(data1$samt))
data1$pcolour[data1$pcolour > 0] <- "#DF010180"
data1$pcolour[data1$pcolour != "#DF010180"] <- "#08088A80"

png(file="scripts/users/imyerssmith/shrub_synthesis/FigS15allindices_temp.png", width=2000, height=2000)
par(mfrow=c(4, 4), mar=c(10, 10, 4, 4), oma=c(2, 2, 2, 2), mgp=c(5, 2, 0))

# amt
x <- as.numeric(as.character(data$samt))
y <- as.numeric(as.character(data$ndaic))
z <- abs(as.numeric(as.character(data$nest_summer)))
se <- as.numeric(as.character(data$nse_summer))
w <- as.numeric(as.character(data$nprop))
v <- as.numeric(as.character(data$R2MM_best))
v1 <- as.numeric(as.character(data1$R2MM_best))
x1 <- as.numeric(as.character(data1$samt))

# daic
plot(y~x, xlab=expression(paste("Annual Mean Temperature (", degree, "C)")), ylab=expression(paste(Delta, "AIC")), cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data$pcolour, cex=3, xlim=c(min(x), max(x)), ylim=c(0, max(y)+(max(y)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x, y)
polygon(x[conv.hull], y[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "A. Annual Mean Temperature"
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# R2W
plot(v1~x1, xlab=expression(paste("Annual Mean Temperature (", degree, "C)")), ylab=expression(paste("R"^"2")), cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data1$pcolour, cex=3, xlim=c(min(x1), max(x1)), ylim=c(0, max(v1)+(max(v1)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x1, v1)
polygon(x1[conv.hull], v1[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "B."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# slope
xz <- na.omit(cbind.data.frame(x, z, se))
plot(xz$z~xz$x, xlab=expression(paste("Annual Mean Temperature (", degree, "C)")), ylab="Abs. model slopes", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data$pcolour, cex=3, xlim=c(min(xz$x), max(xz$x)), ylim=c(0, max(xz$z)+(max(xz$z)/10)))
abline(h=0, lty = 2, lwd = 2)

zy <- xz$z
zx <- xz$x
er <- xz$se
g <- lm(zy ~ zx, weights=(1/er^2))
curve <- data.frame(curve(predict(g, data.frame(zx=x), type="resp"), add=TRUE), lwd=2)

p <- round(summary(g)[5][[1]][8], digits=2)
r2 <- round(summary(g)[9][[1]], digits=2)
p2 <- if(p<0.01) {bquote(C.~p~"<"~0.01~", "~R[2]~"="~.(r2))} else {bquote(C.~p~"="~.(p)~", "~R[2]~"="~.(r2))}
leg.txt <- p2
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

i <- 1
CI <- as.data.frame(array(0, c(20, 4)))
colnames(CI) <- c("x", "est", "down", "up")
range <- seq(min(x), max(x), length.out=20)

for (i in 1:20){
  CI[i, 1] <- as.numeric(range[i])
  CI[i, 2] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[2])
  CI[i, 3] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[7])
  CI[i, 4] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[8])
  i <- i+1
}

polygon(c(CI$x, rev(CI$x)), c(CI$up, rev(CI$down)), col = "#84848425", border = NA)

# prop
plot(w~x, xlab=expression(paste("Annual Mean Temperature (", degree, "C)")), ylab="Proportion", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data$pcolour, cex=3, xlim=c(min(x), max(x)), ylim=c(0, max(w)+(max(w)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x, w)
polygon(x[conv.hull], w[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "D."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# Growing season
gsl <- cbind(data$sgs, data$ndaic, data$ndaic_summer, as.numeric(as.character(data$nest_best)), as.numeric(as.character(data$nse_best)), as.numeric(as.character(data$R2MM_best)), as.numeric(as.character(data$nest_summer)), as.numeric(as.character(data$R2MM_summer)), as.numeric(as.character(data$nse_summer)), as.numeric(as.character(data$nprop)), as.character(data$pcolour))
colnames(gsl) <- c("sgs", "ndaic", "ndaic_summer", "nest_best", "nse_best", "R2MM_best", "nest_summer", "R2MM_summer", "nse_summer", "nprop", "pcolour")
gsl <- as.data.frame(na.omit(gsl))

gsl1 <- cbind(data1$sgs, data1$ndaic, data1$ndaic_summer, as.numeric(as.character(data1$nest_best)), as.numeric(as.character(data1$nse_best)), as.numeric(as.character(data1$R2MM_best)), as.numeric(as.character(data1$nest_summer)), as.numeric(as.character(data1$R2MM_summer)), as.numeric(as.character(data1$nse_summer)), as.numeric(as.character(data1$nprop)), as.character(data1$pcolour))
colnames(gsl1) <- c("sgs", "ndaic", "ndaic_summer", "nest_best", "nse_best", "R2MM_best", "nest_summer", "R2MM_summer", "nse_summer", "nprop", "pcolour")
gsl1 <- as.data.frame(na.omit(gsl1))

x <- as.numeric(as.character(gsl$sgs))
y <- as.numeric(as.character(gsl$ndaic))
z <- abs(as.numeric(as.character(gsl$nest_summer)))
se <- as.numeric(as.character(gsl$nse_summer))
w <- as.numeric(as.character(gsl$nprop))
v <- as.numeric(as.character(gsl$R2MM_best))
v1 <- as.numeric(as.character(gsl1$R2MM_best))
x1 <- as.numeric(as.character(gsl1$sgs))

vxdata <- cbind.data.frame(x,v)
vxdata <- subset(vxdata, v != "NA")
v1 <- vxdata$v
x1 <- vxdata$x

# daic
plot(y~x, xlab="Growing Season Length (weeks)", ylab=expression(paste(Delta, "AIC")), cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=as.character(gsl$pcolour), cex=3, xlim=c(min(x), max(x)), ylim=c(0, max(y)+(max(y)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x, y)
polygon(x[conv.hull], y[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "E. Growing Season Length"
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# R2W
plot(v1~x1, xlab="Growing Season Length (weeks)", ylab=expression(paste("R"^"2")), cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=as.character(gsl1$pcolour), cex=3, xlim=c(min(x1), max(x1)), ylim=c(0, max(v1)+(max(v1)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x1, v1)
polygon(x1[conv.hull], v1[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "F."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# slope
xz <- na.omit(cbind.data.frame(x, z, se))
plot(xz$z~xz$x, xlab="Growing Season Length (weeks)", ylab="Abs. model slopes", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=as.character(gsl$pcolour), cex=3, xlim=c(min(xz$x), max(xz$x)), ylim=c(0, max(xz$z)+(max(xz$z)/10)))
abline(h=0, lty = 2, lwd = 2)

zy <- xz$z
zx <- xz$x
er <- xz$se
g <- lm(zy ~ zx, weights=(1/er^2))
curve <- data.frame(curve(predict(g, data.frame(zx=x), type="resp"), add=TRUE), lwd=2)

p <- round(summary(g)[5][[1]][8], digits=2)
r2 <- round(summary(g)[9][[1]], digits=2)
p2 <- if(p<0.01) {bquote(G.~p~"<"~0.01~", "~R[2]~"="~.(r2))} else {bquote(G.~p~"="~.(p)~", "~R[2]~"="~.(r2))}
leg.txt <- p2
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

i <- 1
CI <- as.data.frame(array(0, c(20, 4)))
colnames(CI) <- c("x", "est", "down", "up")
range <- seq(min(x), max(x), length.out=20)

for (i in 1:20){
  CI[i, 1] <- as.numeric(range[i])
  CI[i, 2] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[2])
  CI[i, 3] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[7])
  CI[i, 4] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[8])
  i <- i+1
}

polygon(c(CI$x, rev(CI$x)), c(CI$up, rev(CI$down)), col = "#84848425", border = NA)

# prop
plot(w~x, xlab="Growing Season Length (weeks)", ylab="Proportion", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=as.character(gsl$pcolour), cex=3, xlim=c(min(x), max(x)), ylim=c(0, max(w)+(max(w)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x, w)
polygon(x[conv.hull], w[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "H."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# Frost day frequency
x <- as.numeric(as.character(data$sfdf))
y <- as.numeric(as.character(data$ndaic))
z <- abs(as.numeric(as.character(data$nest_summer)))
se <- as.numeric(as.character(data$nse_summer))
w <- as.numeric(as.character(data$nprop))
v <- as.numeric(as.character(data$R2MM_best))
v1 <- as.numeric(as.character(data1$R2MM_best))
x1 <- as.numeric(as.character(data1$sfdf))

# daic
plot(y~x, xlab="Frost Day Frequency (days)", ylab=expression(paste(Delta, "AIC")), cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data$pcolour, cex=3, xlim=c(min(x), max(x)), ylim=c(0, max(y)+(max(y)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x, y)
polygon(x[conv.hull], y[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "I. Frost Day Frequency"
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# R2W
plot(v1~x1, xlab="Frost Day Frequency (days)", ylab=expression(paste("R"^"2")), cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data1$pcolour, cex=3, xlim=c(min(x1), max(x1)), ylim=c(0, max(v1)+(max(v1)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x1, v1)
polygon(x1[conv.hull], v1[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "J."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# slope
xz <- na.omit(cbind.data.frame(x, z, se))
plot(xz$z~xz$x, xlab="Frost Day Frequency (days)", ylab="Abs. model slopes", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data$pcolour, cex=3, xlim=c(min(xz$x), max(xz$x)), ylim=c(0, max(xz$z)+(max(xz$z)/10)))
abline(h=0, lty = 2, lwd = 2)

zy <- xz$z
zx <- xz$x
er <- xz$se
g <- lm(zy ~ zx, weights=(1/er^2))
curve <- data.frame(curve(predict(g, data.frame(zx=x), type="resp"), add=TRUE), lwd=2)

p <- round(summary(g)[5][[1]][8], digits=2)
r2 <- round(summary(g)[9][[1]], digits=2)
p2 <- if(p<0.01) {bquote(K.~p~"<"~0.01~", "~R[2]~"="~.(r2))} else {bquote(K.~p~"="~.(p)~", "~R[2]~"="~.(r2))}
leg.txt <- p2
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

i <- 1
CI <- as.data.frame(array(0, c(20, 4)))
colnames(CI) <- c("x", "est", "down", "up")
range <- seq(min(x), max(x), length.out=20)

for (i in 1:20){
  CI[i, 1] <- as.numeric(range[i])
  CI[i, 2] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[2])
  CI[i, 3] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[7])
  CI[i, 4] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[8])
  i <- i+1
}

polygon(c(CI$x, rev(CI$x)), c(CI$up, rev(CI$down)), col = "#84848425", border = NA)

# prop
plot(w~x, xlab="Frost Day Frequency (days)", ylab="Proportion", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data$pcolour, cex=3, xlim=c(min(x), max(x)), ylim=c(0, max(w)+(max(w)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x, w)
polygon(x[conv.hull], w[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "L."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# Clouds
x <- as.numeric(as.character(data$ssumcloud))
y <- as.numeric(as.character(data$ndaic))
z <- abs(as.numeric(as.character(data$nest_summer)))
se <- as.numeric(as.character(data$nse_summer))
w <- as.numeric(as.character(data$nprop))
v <- as.numeric(as.character(data$R2MM_best))
v1 <- as.numeric(as.character(data1$R2MM_best))
x1 <- as.numeric(as.character(data1$ssumcloud))

# daic
plot(y~x, xlab="Cloud Cover (%)", ylab=expression(paste(Delta, "AIC")), cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data$pcolour, cex=3, xlim=c(min(x), max(x)), ylim=c(0, max(y)+(max(y)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x, y)
polygon(x[conv.hull], y[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "M. Cloud Cover"
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# R2W
plot(v1~x1, xlab="Cloud Cover (%)", ylab=expression(paste("R"^"2")), cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data1$pcolour, cex=3, xlim=c(min(x1), max(x1)), ylim=c(0, max(v1)+(max(v1)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x1, v1)
polygon(x1[conv.hull], v1[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "N."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# slope
xz <- na.omit(cbind.data.frame(x, z, se))
plot(xz$z~xz$x, xlab="Cloud Cover (%)", ylab="Abs. model slopes", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data$pcolour, cex=3, xlim=c(min(xz$x), max(xz$x)), ylim=c(0, max(xz$z)+(max(xz$z)/10)))
abline(h=0, lty = 2, lwd = 2)

zy <- xz$z
zx <- xz$x
er <- xz$se
g <- lm(zy ~ zx, weights=(1/er^2))
curve <- data.frame(curve(predict(g, data.frame(zx=x), type="resp"), add=TRUE), lwd=2)

p <- round(summary(g)[5][[1]][8], digits=2)
r2 <- round(summary(g)[9][[1]], digits=2)
p2 <- if(p<0.01) {bquote(O.~p~"<"~0.01~", "~R[2]~"="~.(r2))} else {bquote(O.~p~"="~.(p)~", "~R[2]~"="~.(r2))}
leg.txt <- p2
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

i <- 1
CI <- as.data.frame(array(0, c(20, 4)))
colnames(CI) <- c("x", "est", "down", "up")
range <- seq(min(x), max(x), length.out=20)

for (i in 1:20){
  CI[i, 1] <- as.numeric(range[i])
  CI[i, 2] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[2])
  CI[i, 3] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[7])
  CI[i, 4] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[8])
  i <- i+1
}

polygon(c(CI$x, rev(CI$x)), c(CI$up, rev(CI$down)), col = "#84848425", border = NA)

# prop
plot(w~x, xlab="Cloud Cover (%)", ylab="Proportion", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col=data$pcolour, cex=3, xlim=c(min(x), max(x)), ylim=c(0, max(w)+(max(w)/10)))
abline(h=0, lty = 2, lwd = 2)
conv.hull <- chull(x, w)
polygon(x[conv.hull], w[conv.hull], border = "#00000050", lty = 2, lwd = 2)
leg.txt <- "P."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

dev.off()

# All indices - Range limits ------------------------------------------------
png(file="scripts/users/imyerssmith/shrub_synthesis/FigS16allindices_range.png", width=2000, height=2000)
par(mfrow=c(4, 4), mar=c(10, 10, 4, 4), oma=c(2, 2, 2, 2), mgp=c(5, 2, 0))

# Range Edge
data <- modelsum
data_t <- subset(data, dwarf == "tall")
data_d <- subset(data, dwarf == "dwarf")

data1 <- subset(modelsum, modelsum$R2MM_best != "NA")
data1_t <- subset(data1, dwarf == "tall")
data1_d <- subset(data1, dwarf == "dwarf")

x <- data$disrelat
y <- as.numeric(as.character(data$ndaic))
z <- abs(as.numeric(as.character(data$nest_summer)))
se <- as.numeric(as.character(data$nse_summer))
w <- as.numeric(as.character(data$nprop))
v <- as.numeric(as.character(data$R2MM_best))
v1a <- as.numeric(as.character(data1$R2MM_best))
x1a <- data1$disrelat

x1 <- data_t$disrelat
y1 <- as.numeric(as.character(data_t$ndaic))
x2 <- data_d$disrelat
y2 <- as.numeric(as.character(data_d$ndaic))
z1 <- abs(as.numeric(as.character(data_t$nest_summer)))
z2 <- abs(as.numeric(as.character(data_d$nest_summer)))
se1 <- as.numeric(as.character(data_t$nse_summer))
se2 <- as.numeric(as.character(data_d$nse_summer))
w1 <- as.numeric(as.character(data_t$nprop))
w2 <- as.numeric(as.character(data_d$nprop))
v1 <- as.numeric(as.character(data_t$R2MM_best))
v2 <- as.numeric(as.character(data_d$R2MM_best)) 
v1t <- as.numeric(as.character(data1_t$R2MM_best))
x1t <- data1_t$disrelat
v2d <- as.numeric(as.character(data1_d$R2MM_best))
x2d <- data1_d$disrelat

# daic
plot(y~x, xlab="Dist. Range Edge (Rel. lat.)", ylab=expression(paste("Climate Sensitivity (", Delta, " AIC)")), cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col="#FFFFFF00", cex=3, xlim=c(0, max(x)), ylim=c(0, max(y)+(max(y)/10)))
conv.hull <- chull(x, y)
polygon(x[conv.hull], y[conv.hull], border = "#00000050", lty = 2, lwd = 2)
points(y1~x1, pch=19, col="#088A8580", cex=3)
points(y2~x2, pch=19, col="#DF01A580", cex=3)
abline(h=0, lty = 2, lwd = 2)
leg.txt <- "A. Distance to N/Ele. range edge"
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# R2W
plot(v1a~x1a, xlab="Dist. Range Edge (Rel. lat.)", ylab=expression(paste("R"^"2")), cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col="#FFFFFF00", cex=3, xlim=c(0, max(x1a)), ylim=c(0, max(v1a)+(max(v1a)/10)))
conv.hull <- chull(x1a, v1a)
polygon(x1a[conv.hull], v1a[conv.hull], border = "#00000050", lty = 2, lwd = 2)
points(v1t~x1t, pch=19, col="#088A8580", cex=3)
points(v2d~x2d, pch=19, col="#DF01A580", cex=3)
abline(h=0, lty = 2, lwd = 2)
leg.txt <- "B."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# slope
xz <- na.omit(cbind.data.frame(x, z, se))
plot(xz$z~xz$x, xlab="Dist. Range Edge (Rel. lat.)", ylab="Abs. model slopes", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col="#FFFFFF00", cex=3, xlim=c(min(xz$x), max(xz$x)), ylim=c(0, max(xz$z)+(max(xz$z)/10)))
points(z1~x1, pch=19, col="#088A8580", cex=3)
points(z2~x2, pch=19, col="#DF01A580", cex=3)
abline(h=0, lty = 2, lwd = 2)

zy <- xz$z
zx <- xz$x
er <- xz$se
g <- lm(zy ~ zx, weights=(1/er^2))
curve <- data.frame(curve(predict(g, data.frame(zx=x), type="resp"), add=TRUE), lwd=2)

p <- round(summary(g)[5][[1]][8], digits=2)
r2 <- round(summary(g)[9][[1]], digits=2)
p2 <- if(p<0.01) {bquote(G.~p~"<"~0.01~", "~R[2]~"="~.(r2))} else {bquote(G.~p~"="~.(p)~", "~R[2]~"="~.(r2))}
leg.txt <- p2
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

i <- 1
CI <- as.data.frame(array(0, c(20, 4)))
colnames(CI) <- c("x", "est", "down", "up")
range <- seq(min(x), max(x), length.out=20)

for (i in 1:20){
  CI[i, 1] <- as.numeric(range[i])
  CI[i, 2] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[2])
  CI[i, 3] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[7])
  CI[i, 4] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[8])
  i <- i+1
}

polygon(c(CI$x, rev(CI$x)), c(CI$up, rev(CI$down)), col = "#84848425", border = NA)

# prop
plot(w~x, xlab="Dist. Range Edge (Rel. lat.)", ylab="Proportion", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col="#FFFFFF00", cex=3, xlim=c(0, max(x)), ylim=c(0, max(w)+(max(w)/10)))
conv.hull <- chull(x, w)
polygon(x[conv.hull], w[conv.hull], border = "#00000050", lty = 2, lwd = 2)
points(w1~x1, pch=19, col="#088A8580", cex=3)
points(w2~x2, pch=19, col="#DF01A580", cex=3)
abline(h=0, lty = 2, lwd = 2)
leg.txt <- "D."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# Canopy
x <- data$canopy
x1 <- data_t$canopy
x2 <- data_d$canopy
vxdata <- cbind.data.frame(x,v)
vxdata <- subset(vxdata, v != "NA")
v1a <- vxdata$v
x1a <- vxdata$x
v1t <- as.numeric(as.character(data1_t$R2MM_best))
x1t <- data1_t$canopy
v2d <- as.numeric(as.character(data1_d$R2MM_best))
x2d <- data1_d$canopy

# daic
plot(y~x, xlab="Canopy Height (cm)", ylab=expression(paste("Climate Sensitivity (", Delta, " AIC)")), cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col="#FFFFFF00", cex=3, xlim=c(0, max(x)), ylim=c(0, max(y)+(max(y)/10)))
conv.hull <- chull(x, y)
polygon(x[conv.hull], y[conv.hull], border = "#00000050", lty = 2, lwd = 2)
points(y1~x1, pch=19, col="#088A8580", cex=3)
points(y2~x2, pch=19, col="#DF01A580", cex=3)
abline(h=0, lty = 2, lwd = 2)
leg.txt <- "E. Canopy height"
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# R2W
plot(v1a~x1a, xlab="Canopy Height (cm)", ylab=expression(paste("R"^"2")), cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col="#FFFFFF00", cex=3, xlim=c(0, max(x1a)), ylim=c(0, max(v1a)+(max(v1a)/10)))
conv.hull <- chull(x1a, v1a)
polygon(x1a[conv.hull], v1a[conv.hull], border = "#00000050", lty = 2, lwd = 2)
points(v1t~x1t, pch=19, col="#088A8580", cex=3)
points(v2d~x2d, pch=19, col="#DF01A580", cex=3)
abline(h=0, lty = 2, lwd = 2)
leg.txt <- "F."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# slope
xz <- na.omit(cbind.data.frame(x, z, se))
plot(xz$z~xz$x, xlab="Canopy Height (cm)", ylab="Abs. model slopes", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col="#FFFFFF00", cex=3, xlim=c(min(xz$x), max(xz$x)), ylim=c(0, max(xz$z)+(max(xz$z)/10)))
points(z1~x1, pch=19, col="#088A8580", cex=3)
points(z2~x2, pch=19, col="#DF01A580", cex=3)
abline(h=0, lty = 2, lwd = 2)

zy <- xz$z
zx <- xz$x
er <- xz$se
g <- lm(zy ~ zx, weights=(1/er^2))
curve <- data.frame(curve(predict(g, data.frame(zx=x), type="resp"), add=TRUE), lwd=2)

p <- round(summary(g)[5][[1]][8], digits=2)
r2 <- round(summary(g)[9][[1]], digits=2)
p2 <- if(p<0.01) {bquote(G.~p~"<"~0.01~", "~R[2]~"="~.(r2))} else {bquote(G.~p~"="~.(p)~", "~R[2]~"="~.(r2))}
leg.txt <- p2
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

i <- 1
CI <- as.data.frame(array(0, c(20, 4)))
colnames(CI) <- c("x", "est", "down", "up")
range <- seq(min(x), max(x), length.out=20)

for (i in 1:20){
  CI[i, 1] <- as.numeric(range[i])
  CI[i, 2] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[2])
  CI[i, 3] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[7])
  CI[i, 4] <- as.numeric(esticon(g, cm = c(1, as.numeric(range[i])), conf.int=TRUE)[8])
  i <- i+1
}

polygon(c(CI$x, rev(CI$x)), c(CI$up, rev(CI$down)), col = "#84848425", border = NA)

# prop
plot(w~x, xlab="Canopy Height (cm)", ylab="Proportion", cex.axis=2.5, cex.lab=2.8, bty="l", pch=19, col="#FFFFFF00", cex=3, xlim=c(0, max(x)), ylim=c(0, max(w)+(max(w)/10)))
conv.hull <- chull(x, w)
polygon(x[conv.hull], w[conv.hull], border = "#00000050", lty = 2, lwd = 2)
points(w1~x1, pch=19, col="#088A8580", cex=3)
points(w2~x2, pch=19, col="#DF01A580", cex=3)
abline(h=0, lty = 2, lwd = 2)
leg.txt <- "H."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# Genera
x <- data$genus_num
y <- data$ndaic
vxdata <- cbind.data.frame(x,v)
vxdata <- subset(vxdata, v != "NA")
v1a <- vxdata$v
x1a <- vxdata$x

data_t <- subset(data, dwarf == "tall")
x1 <- as.numeric(as.character(data_t$genus_num))
y1 <- as.numeric(as.character(data_t$ndaic))

data_d <- subset(data, dwarf == "dwarf")
x2 <- as.numeric(as.character(data_d$genus_num))
y2 <- as.numeric(as.character(data_d$ndaic))
v1t <- as.numeric(as.character(data1_t$R2MM_best))
x1t <- as.numeric(as.character(data1_t$genus_num))
v2d <- as.numeric(as.character(data1_d$R2MM_best))
x2d <- as.numeric(as.character(data1_d$genus_num))

# daic
par(bty="l")
stripchart(as.numeric(y)~as.factor(x), ylab=expression(paste("Climate Sensitivity (", Delta, " AIC)")), vertical = TRUE, cex.axis=0.5, col="#FFFFFF00", xlab="", cex.axis=2.5, cex.lab=2.8, pch=19, cex=3, xaxt="n", ylim=c(0, max(y)+(max(y)/10)))
points(y1~x1, pch=19, col="#088A8580", cex=3)
points(y2~x2, pch=19, col="#DF01A580", cex=3)
text(1:length(unique(data$genus_name)), 0-(max(y)/10), srt = 45, adj = 1, labels=unique(as.factor(data$genus_name)), xpd = TRUE, cex=2.5)
abline(h=0, lty = 2, lwd = 2)
leg.txt <- "I. Genera"
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# R2W
par(bty="l")
stripchart(as.numeric(v1a)~as.factor(x1a), ylab=expression(paste("R"^"2")), vertical = TRUE, cex.axis=0.5, col="#FFFFFF00", xlab="", cex.axis=2.5, cex.lab=2.8, pch=19, cex=3, xaxt="n", ylim=c(0, max(v1a)+(max(v1a)/10)))
points(v1t~x1t, pch=19, col="#088A8580", cex=3)
points(v2d~x2d, pch=19, col="#DF01A580", cex=3)
text(1:length(unique(data$genus_name)), 0-(max(v1a)/10), srt = 45, adj = 1, labels=unique(as.factor(data$genus_name)), xpd = TRUE, cex=2.5)
abline(h=0, lty = 2, lwd = 2)
leg.txt <- "J."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# slope
xz <- na.omit(cbind.data.frame(x, z))
par(bty="l")
stripchart(as.numeric(z)~as.factor(x), ylab="Abs. model slopes", vertical = TRUE, cex.axis=0.5, col="#FFFFFF00", xlab="", cex.axis=2.5, cex.lab=2.8, pch=19, cex=3, xaxt="n", ylim=c(0, max(xz$z)+(max(xz$z)/10)))
points(z1~x1, pch=19, col="#088A8580", cex=3)
points(z2~x2, pch=19, col="#DF01A580", cex=3)
text(1:length(unique(data$genus_name)), 0-(max(xz$z)/10), srt = 45, adj = 1, labels=unique(as.factor(data$genus_name)), xpd = TRUE, cex=2.5)
abline(h=0, lty = 2, lwd = 2)
leg.txt <- "K."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# prop
par(bty="l")
stripchart(as.numeric(w)~as.factor(x), ylab="Proportion", vertical = TRUE, cex.axis=0.5, col="#FFFFFF00", xlab="", cex.axis=2.5, cex.lab=2.8, pch=19, cex=3, xaxt="n", ylim=c(0, max(w)+(max(w)/10)))
points(w1~x1, pch=19, col="#088A8580", cex=3)
points(w2~x2, pch=19, col="#DF01A580", cex=3)
text(1:length(unique(data$genus_name)), 0-(max(w)/10), srt = 45, adj = 1, labels=unique(as.factor(data$genus_name)), xpd = TRUE, cex=2.5)
abline(h=0, lty = 2, lwd = 2)
leg.txt <- "L."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# Data type
x <- gsub("increment", "Increment", gsub("ring", "Ring width", data$datatype))
y <- as.numeric(as.character(data$ndaic))
vxdata <- cbind.data.frame(x,v)
vxdata <- subset(vxdata, v != "NA")
v1a <- vxdata$v
x1a <- vxdata$x

data_i <- subset(data, datatype == "increment")
data_i_t <- subset(data_i, dwarf == "tall")
data_i_d <- subset(data_i, dwarf == "dwarf")
x1t <- rep(1, length(data_i_t$datatype))
y1t <- as.numeric(as.character(data_i_t$ndaic))
x1d <- rep(1, length(data_i_d$datatype))
y1d <- as.numeric(as.character(data_i_d$ndaic))
z1t <- abs(as.numeric(as.character(data_i_t$nest_summer)))
z1d <- abs(as.numeric(as.character(data_i_d$nest_summer)))
z1t <- abs(as.numeric(as.character(data_i_t$nest_summer)))
z1d <- abs(as.numeric(as.character(data_i_d$nest_summer)))
se1t <- as.numeric(as.character(data_i_t$nse_summer))
se1d <- as.numeric(as.character(data_i_d$nse_summer))
w1t <- as.numeric(as.character(data_i_t$nprop))
w1d <- as.numeric(as.character(data_i_d$nprop))
v1t <- as.numeric(as.character(data_i_t$R2MM_best))
v1d <- as.numeric(as.character(data_i_d$R2MM_best))
vxdatat <- cbind.data.frame(x1t,v1t)
vxdatat <- subset(vxdatat, v1t != "NA")
v1t1 <- vxdatat$v1t
x1t1 <- vxdatat$x1t
vxdatad <- cbind.data.frame(x1d,v1d)
vxdatad <- subset(vxdatad, v1d != "NA")
v1d1 <- vxdatad$v1d
x1d1 <- vxdatad$x1d

data_r <- subset(data, datatype == "ring")
data_r_t <- subset(data_r, dwarf == "tall")
data_r_d <- subset(data_r, dwarf == "dwarf")
x2t <- rep(2, length(data_r_t$datatype))
y2t <- as.numeric(as.character(data_r_t$ndaic))
x2d <- rep(2, length(data_r_d$datatype))
y2d <- as.numeric(as.character(data_r_d$ndaic))
z2t <- abs(as.numeric(as.character(data_r_t$nest_summer)))
z2d <- abs(as.numeric(as.character(data_r_d$nest_summer)))
z2t <- abs(as.numeric(as.character(data_r_t$nest_summer)))
z2d <- abs(as.numeric(as.character(data_r_d$nest_summer)))
se2t <- as.numeric(as.character(data_r_t$nse_summer))
se2d <- as.numeric(as.character(data_r_d$nse_summer))
w2t <- as.numeric(as.character(data_r_t$nprop))
w2d <- as.numeric(as.character(data_r_d$nprop))
v2t <- as.numeric(as.character(data_r_t$R2MM_best))
v2d <- as.numeric(as.character(data_r_d$R2MM_best))
vxdatat <- cbind.data.frame(x2t,v2t)
vxdatat <- subset(vxdatat, v2t != "NA")
v2t1 <- vxdatat$v2t
x2t1 <- vxdatat$x2t
vxdatad <- cbind.data.frame(x2d,v2d)
vxdatad <- subset(vxdatad, v2d != "NA")
v2d1 <- vxdatad$v2d
x2d1 <- vxdatad$x2d

# daic
par(bty="l")
stripchart(as.numeric(y)~as.factor(x), vertical = TRUE, cex.axis=0.5, col="#FFFFFF00", ylab=expression(paste("Climate Sensitivity (", Delta, " AIC)")), xlab="", cex.axis=2.5, cex.lab=2.8, pch=19, cex=3, xaxt="n", xlim=c(0, 3), ylim=c(0, max(y)+(max(y)/10)))
points(y1t~x1t, pch=19, col="#088A8580", cex=3)
points(y1d~x1d, pch=19, col="#DF01A580", cex=3)
points(y2t~x2t, pch=19, col="#088A8580", cex=3)
points(y2d~x2d, pch=19, col="#DF01A580", cex=3)
text(1:length(unique(x)), 0-(max(y)/10), srt = 45, adj = 1, labels=unique(as.factor(x)), xpd = TRUE, cex=2.5)
abline(h=0, lty = 2, lwd = 2)
leg.txt <- "M. Growth measure"
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# R2W
par(bty="l")
stripchart(as.numeric(v1a)~as.factor(x1a), vertical = TRUE, cex.axis=0.5, col="#FFFFFF00", ylab=expression(paste("R"^"2")), xlab="", cex.axis=2.5, cex.lab=2.8, pch=19, cex=3, xaxt="n", xlim=c(0, 3), ylim=c(0, max(v1a)+(max(v1a)/10)))
points(v1t1~x1t1, pch=19, col="#088A8580", cex=3)
points(v1d1~x1d1, pch=19, col="#DF01A580", cex=3)
points(v2t1~x2t1, pch=19, col="#088A8580", cex=3)
points(v2d1~x2d1, pch=19, col="#DF01A580", cex=3)
text(1:length(unique(x1a)), 0-(max(v1a)/10), srt = 45, adj = 1, labels=unique(as.factor(x)), xpd = TRUE, cex=2.5)
abline(h=0, lty = 2, lwd = 2)
leg.txt <- "N."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# slope
xz <- na.omit(cbind.data.frame(x, z))
par(bty="l")
stripchart(abs(as.numeric(xz$z))~as.factor(xz$x), vertical = TRUE, cex.axis=0.5, col="#FFFFFF00", ylab="Abs. model slopes", xlab="", cex.axis=2.5, cex.lab=2.8, pch=19, cex=3, xaxt="n", xlim=c(0, 3), ylim=c(0, max(abs(xz$z))+(max(abs(xz$z))/10)))
points(abs(z1t)~x1t, pch=19, col="#088A8580", cex=3)
points(abs(z1d)~x1d, pch=19, col="#DF01A580", cex=3)
points(abs(z2t)~x2t, pch=19, col="#088A8580", cex=3)
points(abs(z2d)~x2d, pch=19, col="#DF01A580", cex=3)
text(1:length(unique(x)), 0-(max(xz$z)/10), srt = 45, adj = 1, labels=unique(as.factor(x)), xpd = TRUE, cex=2.5)
abline(h=0, lty = 2, lwd = 2)
leg.txt <- "O."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

# prop
par(bty="l")
stripchart(as.numeric(w)~as.factor(x), vertical = TRUE, cex.axis=0.5, col="#FFFFFF00", ylab="Proportion", xlab="", cex.axis=2.5, cex.lab=2.8, pch=19, cex=3, xaxt="n", xlim=c(0, 3), ylim=c(0, max(w)+(max(w)/10)))
points(w1t~x1t, pch=19, col="#088A8580", cex=3)
points(w1d~x1d, pch=19, col="#DF01A580", cex=3)
points(w2t~x2t, pch=19, col="#088A8580", cex=3)
points(w2d~x2d, pch=19, col="#DF01A580", cex=3)
text(1:length(unique(x)), 0-(max(w)/10), srt = 45, adj = 1, labels=unique(as.factor(x)), xpd = TRUE, cex=2.5)
abline(h=0, lty = 2, lwd = 2)
leg.txt <- "P."
legend("topleft", inset=c(-0.05, -0.05), legend=leg.txt, x.intersp=1.5, y.intersp=1.4, bty="n", cex = 2.5)

dev.off()
