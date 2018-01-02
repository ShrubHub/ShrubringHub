# Shrub growth synthesis mixed model script:
# This code runs the mixed model analysis for the shrub synthesis manuscript.
# Written by: Isla Myers-Smith (e-mail: isla.myers-smith@ed.ac.uk)
# 3 November 2014

# Libraries
library(plyr)
library(nlme)
library(Matrix)
library(stats)
library(lmtest)
library(MuMIn)
library(ggplot2)
library(gridExtra)

# Load data
load("workspace/shrubhub.RData")
load("workspace/climate_data.RData")

siteinfosum <- read.csv("scripts/users/imyerssmith/shrub_synthesis/site_info.csv", 
                        header=TRUE, stringsAsFactors=FALSE)
siteinfosum <- as.data.frame(unique(cbind(siteinfosum$sgnum, siteinfosum$site, siteinfosum$site_code, siteinfosum$lat, siteinfosum$lon)))
colnames(siteinfosum) <- c("sgnum", "site", "site_code", "lat", "lon")

# Normalize climate data

# Standardize
# normalize <- function(x) { scale(x, scale = FALSE) }

# Normalize
normalize <- function(x) { scale(x, scale = TRUE) }

climate_datan <- climate_data

climate_datan <- ddply(climate_datan, .(site), transform, nfdf = normalize(fdf))
climate_datan <- ddply(climate_datan, .(site), transform, namt = normalize(amt))
climate_datan <- ddply(climate_datan, .(site), transform, njunetmp = normalize(junetmp))
climate_datan <- ddply(climate_datan, .(site), transform, njunetmx = normalize(junetmx))
climate_datan <- ddply(climate_datan, .(site), transform, njunetmn = normalize(junetmn))
climate_datan <- ddply(climate_datan, .(site), transform, njulytmp = normalize(julytmp))
climate_datan <- ddply(climate_datan, .(site), transform, njulytmx = normalize(julytmx))
climate_datan <- ddply(climate_datan, .(site), transform, njulytmn = normalize(julytmn))
climate_datan <- ddply(climate_datan, .(site), transform, naugtmp = normalize(augtmp))
climate_datan <- ddply(climate_datan, .(site), transform, naugtmx = normalize(augtmx))
climate_datan <- ddply(climate_datan, .(site), transform, naugtmn = normalize(augtmn))
climate_datan <- ddply(climate_datan, .(site), transform, njjtmp = normalize(jjtmp))
climate_datan <- ddply(climate_datan, .(site), transform, njjtmx = normalize(jjtmx))
climate_datan <- ddply(climate_datan, .(site), transform, njjtmn = normalize(jjtmn))
climate_datan <- ddply(climate_datan, .(site), transform, nastmp = normalize(astmp))
climate_datan <- ddply(climate_datan, .(site), transform, nastmx = normalize(astmx))
climate_datan <- ddply(climate_datan, .(site), transform, nastmn = normalize(astmn))
climate_datan <- ddply(climate_datan, .(site), transform, nwintmp = normalize(wintmp))
climate_datan <- ddply(climate_datan, .(site), transform, nwintmx = normalize(wintmx))
climate_datan <- ddply(climate_datan, .(site), transform, nwintmn = normalize(wintmn))
climate_datan <- ddply(climate_datan, .(site), transform, nmamtmp = normalize(mamtmp))
climate_datan <- ddply(climate_datan, .(site), transform, nmamtmx = normalize(mamtmx))
climate_datan <- ddply(climate_datan, .(site), transform, nmamtmn = normalize(mamtmn))
climate_datan <- ddply(climate_datan, .(site), transform, njunepre = normalize(junepre))
climate_datan <- ddply(climate_datan, .(site), transform, njulypre = normalize(julypre))
climate_datan <- ddply(climate_datan, .(site), transform, naugpre = normalize(augpre))
climate_datan <- ddply(climate_datan, .(site), transform, njjpre = normalize(jjpre))
climate_datan <- ddply(climate_datan, .(site), transform, naspre = normalize(aspre))
climate_datan <- ddply(climate_datan, .(site), transform, nwinpre = normalize(winpre))
climate_datan <- ddply(climate_datan, .(site), transform, nmampre = normalize(mampre))
climate_datan <- ddply(climate_datan, .(site), transform, nsumcloud = normalize(sumcloud))
climate_datan <- ddply(climate_datan, .(site), transform, nwetd = normalize(wetd))
climate_datan <- ddply(climate_datan, .(site), transform, ngscru.length = normalize(gscru.length))
climate_datan <- ddply(climate_datan, .(site), transform, nspring.tdd = normalize(spring.tdd))

# Merge data frames
data <- merge(growth.full, climate_data)
data <- merge(data, climate_datan)
data[data == "NaN"] <- "NA"

unfactorize <- function(df){
  for(i in which(sapply(df, class) == "factor")) df[[i]] = as.character(df[[i]])
  return(df)
}

data <- unfactorize(data)

x <- 1950
y <- 1
j <- 1
k <- 1
lengthx <- max(unique(data$year)) - (min(unique(data$year))+20)
lengthsgnum <- length(unique(data$sgnum))
model_comp_time_slope <- as.data.frame(array(0, c((lengthsgnum*34), 24)))
samplesize <- as.data.frame(array(0, c(lengthsgnum*lengthx, 3)))
colnames(samplesize) <- c("timewindow", "sgnum", "nshrub")

for (x in 1950:(1950+lengthx)) { 
  norm_data <- data
  t1 <- x
  t2 <- x + 20
  
  # Subset by time window
  norm_data <- subset(data, year > t1 & year < t2)
  
  # Normalize the data
  # Standardize
  # normalize <- function(x) { scale(x, scale = FALSE) }
  
  # Normalize
  normalize <- function(x) { scale(x, scale = TRUE) }
  
  norm_data <- ddply(norm_data, .(shrub_num), transform, nrw = scale(rw, scale = TRUE))
  
# Mixed model analysis slopes
ndata <- cbind.data.frame(norm_data$site, norm_data$snum, norm_data$sgnum, norm_data$shrub_num, norm_data$year, norm_data$xyear, norm_data$nrw, norm_data$nfdf, norm_data$namt, norm_data$njunetmp, norm_data$njunetmx, norm_data$njunetmn, norm_data$njulytmp, norm_data$njulytmx, norm_data$njulytmn, norm_data$naugtmp, norm_data$naugtmx, norm_data$naugtmn, norm_data$njjtmp, norm_data$njjtmx, norm_data$njjtmn, norm_data$nastmp, norm_data$nastmx, norm_data$nastmn, norm_data$nwintmp, norm_data$nwintmx, norm_data$nwintmn, norm_data$nmamtmp, norm_data$nmamtmx, norm_data$nmamtmn, norm_data$njunepre, norm_data$njulypre, norm_data$naugpre, norm_data$njjpre, norm_data$naspre, norm_data$nwinpre, norm_data$nmampre)
names(ndata) <- c("site", "snum", "sgnum", "shrub_num", "year", "xyear", "nrw", "nfdf", "namt", "njunetmp", "njunetmx", "njunetmn", "njulytmp", "njulytmx", "njulytmn", "naugtmp", "naugtmx", "naugtmn", "njjtmp", "njjtmx", "njjtmn", "nastmp", "nastmx", "nastmn", "nwintmp", "nwintmx", "nwintmn", "nmamtmp", "nmamtmx", "nmamtmn", "njunepre", "njulypre", "naugpre", "njjpre", "naspre", "nwinpre", "nmampre")
ndata <- na.omit(ndata)
gs <- cbind.data.frame(norm_data$sgnum, norm_data$shrub_num, norm_data$year, norm_data$ngscru.length, norm_data$nspring.tdd)
names(gs) <- c("sgnum", "shrub_num", "year", "ngscru.length", "nspring.tdd")
ndata <- merge(ndata, gs)
mdata <- ndata
i <- 1
j <- 1
k <- 1
l <- 34
length <- length(unique(mdata$sgnum))
coeff <- as.data.frame(array(0, c(34, 7)))
colnames(coeff) <- c("model", "estimate", "se", "estimate2", "se2", "aic", "R2MM")
coeff2 <- as.data.frame(array(0, c(length*34, 8)))
colnames(coeff2) <- c("sgnum", "model", "estimate", "se", "estimate2", "se2", "aic", "R2MM")

 for (x in 1950:(1950+lengthx)) { 
  norm_data <- data
  t1 <- x
  t2 <- x + 20
  
  # Subset by time window
  norm_data <- subset(data, year > t1 & year < t2)
  
  # Normalize the data
  # Standardize
  # normalize <- function(x) { scale(x, scale = FALSE) }
  
  # Normalize
  # normalize <- function(x) { scale(x, scale = TRUE) }

  norm_data <- ddply(norm_data, .(shrub_num), transform, nrw = scale(rw, scale = TRUE))

# Mixed model analysis slopes
ndata <- cbind.data.frame(norm_data$site, norm_data$snum, norm_data$sgnum, norm_data$shrub_num, norm_data$year, norm_data$xyear, norm_data$nrw, norm_data$nfdf, norm_data$namt, norm_data$njunetmp, norm_data$njunetmx, norm_data$njunetmn, norm_data$njulytmp, norm_data$njulytmx, norm_data$njulytmn, norm_data$naugtmp, norm_data$naugtmx, norm_data$naugtmn, norm_data$njjtmp, norm_data$njjtmx, norm_data$njjtmn, norm_data$nastmp, norm_data$nastmx, norm_data$nastmn, norm_data$nwintmp, norm_data$nwintmx, norm_data$nwintmn, norm_data$nmamtmp, norm_data$nmamtmx, norm_data$nmamtmn, norm_data$njunepre, norm_data$njulypre, norm_data$naugpre, norm_data$njjpre, norm_data$naspre, norm_data$nwinpre, norm_data$nmampre)
names(ndata) <- c("site", "snum", "sgnum", "shrub_num", "year", "xyear", "nrw", "nfdf", "namt", "njunetmp", "njunetmx", "njunetmn", "njulytmp", "njulytmx", "njulytmn", "naugtmp", "naugtmx", "naugtmn", "njjtmp", "njjtmx", "njjtmn", "nastmp", "nastmx", "nastmn", "nwintmp", "nwintmx", "nwintmn", "nmamtmp", "nmamtmx", "nmamtmn", "njunepre", "njulypre", "naugpre", "njjpre", "naspre", "nwinpre", "nmampre")
ndata <- na.omit(ndata)
gs <- cbind.data.frame(norm_data$sgnum, norm_data$shrub_num, norm_data$year, norm_data$ngscru.length, norm_data$nspring.tdd)
names(gs) <- c("sgnum", "shrub_num", "year", "ngscru.length", "nspring.tdd")
ndata <- merge(ndata, gs)
mdata <- ndata
i <- 1
j <- 1
k <- 1
l <- 34
length <- length(unique(mdata$sgnum))
coeff <- as.data.frame(array(0, c(34, 7)))
colnames(coeff) <- c("model", "estimate", "se", "estimate2", "se2", "aic", "R2MM")
coeff2 <- as.data.frame(array(0, c(length*34, 8)))
colnames(coeff2) <- c("sgnum", "model", "estimate", "se", "estimate2", "se2", "aic", "R2MM")

# To change between maximum likelyhood (model comparison) and restricted maximum likelyhood (slope estimation)
# use method = "ML" or method = "REML"

for (i in 1:length) {
  mdata <- ndata
  if(all(unique(mdata$sgnum) != i)) { next } else { 
    mdata <- subset(mdata, sgnum == i)
    if(length(unique(mdata$shrub_num)) < 10) { next } else { 
      if(max(mdata$year)-min(mdata$year) < 10) { next } else { 
        m1n <- lme(nrw ~ njjtmp, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))  
        m2n <- lme(nrw ~ njjtmx, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m3n <- lme(nrw ~ njjtmn, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m4n <- lme(nrw ~ njunetmp, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m5n <- lme(nrw ~ njunetmx, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m6n <- lme(nrw ~ njunetmn, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m7n <- lme(nrw ~ njulytmp, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m8n <- lme(nrw ~ njulytmx, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m9n <- lme(nrw ~ njulytmn, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m10n <- lme(nrw ~ nastmp, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m11n <- lme(nrw ~ nastmx, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m12n <- lme(nrw ~ nastmn, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m13n <- lme(nrw ~ nmamtmp, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m14n <- lme(nrw ~ nmamtmx, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m15n <- lme(nrw ~ nmamtmn, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m16n <- lme(nrw ~ nwintmp, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m17n <- lme(nrw ~ nwintmx, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m18n <- lme(nrw ~ nwintmn, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m19n <- lme(nrw ~ njunepre, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m20n <- lme(nrw ~ njulypre, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m21n <- lme(nrw ~ naugpre, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m22n <- lme(nrw ~ njjpre, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m23n <- lme(nrw ~ naspre, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m24n <- lme(nrw ~ nwinpre, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m25n <- lme(nrw ~ nmampre, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m26n <- lme(nrw ~ njunetmp + nwinpre, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m27n <- lme(nrw ~ njunetmp + njunepre, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m28n <- lme(nrw ~ njulytmp + nwinpre, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m29n <- lme(nrw ~ njulytmp + njulypre, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m30n <- lme(nrw ~ naugtmp + nwinpre, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m31n <- lme(nrw ~ naugtmp + naugpre, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m32n <- lme(nrw ~ njjtmp + nwinpre, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        m33n <- lme(nrw ~ njjtmp + njjpre, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        mNULLn <- lme(nrw ~ 1, random=~1|xyear, cor=corAR1(), data = mdata, method = "REML", control=list(maxIter=1000))
        modeln <- list(summary(m1n), summary(m2n), summary(m3n), summary(m4n), summary(m5n), summary(m6n), summary(m7n), summary(m8n), summary(m9n), summary(m10n), summary(m11n), summary(m12n), summary(m13n), summary(m14n), summary(m15n), summary(m16n), summary(m17n), summary(m18n), summary(m19n), summary(m20n), summary(m21n), summary(m22n), summary(m23n), summary(m24n), summary(m25n), summary(m26n), summary(m27n), summary(m28n), summary(m29n), summary(m30n), summary(m31n), summary(m32n), summary(m33n), summary(mNULLn))
        R2modeln <- list(r.squaredGLMM(m1n), r.squaredGLMM(m2n), r.squaredGLMM(m3n), r.squaredGLMM(m4n), r.squaredGLMM(m5n), r.squaredGLMM(m6n), r.squaredGLMM(m7n), r.squaredGLMM(m8n), r.squaredGLMM(m9n), r.squaredGLMM(m10n), r.squaredGLMM(m11n), r.squaredGLMM(m12n), r.squaredGLMM(m13n), r.squaredGLMM(m14n), r.squaredGLMM(m15n), r.squaredGLMM(m16n), r.squaredGLMM(m17n), r.squaredGLMM(m18n), r.squaredGLMM(m19n), r.squaredGLMM(m20n), r.squaredGLMM(m21n), r.squaredGLMM(m22n), r.squaredGLMM(m23n), r.squaredGLMM(m24n), r.squaredGLMM(m25n), r.squaredGLMM(m26n), r.squaredGLMM(m27n), r.squaredGLMM(m28n), r.squaredGLMM(m29n), r.squaredGLMM(m30n), r.squaredGLMM(m31n), r.squaredGLMM(m32n), r.squaredGLMM(m33n), r.squaredGLMM(mNULLn))
        
        for (j in 1:34) { 
          coef <- as.data.frame(modeln[[j]]$tTable)
          #model name
          coeff[j, 1] <- if(is.na(rownames(coef)[2])) { 
            "NULL" 
          } else if(is.na(rownames(coef)[3])) { 
            rownames(coef)[2] 
          } else { 
            paste(rownames(coef)[2], "-", rownames(coef)[3]) 
          }
          #estimate
          coeff[j, 2] <- coef[2, 1]
          #se
          coeff[j, 3] <- coef[2, 2]
          if(is.na(coef[3, 1])) {
            #estimate2
            coeff[j, 4] <- NA
            #se2
            coeff[j, 5] <- NA 
          } else {
            #estimate2
            coeff[j, 4] <- coef[2, 1]
            #se2
            coeff[j, 5] <- coef[2, 2]
          }
          #AIC
          coeff[j, 6] <- modeln[[j]]$AIC
          coeff[j, 7] <- R2modeln[[j]][2]
          j <- j + 1
        }
      }
    }
  }
coeff2[k:l, 1] <- if (length(mdata$sgnum) == 0) { "NA" } else { unique(mdata$sgnum) }
coeffNA <- as.data.frame(array(NA, c(34, 7)))
coeff2[k:l, 2:8] <- if (length(mdata$sgnum) == 0) { coeffNA } else { coeff }
i <- i + 1
j <- 1
k <- k + 34
l <- l + 34
}

models <- as.data.frame(cbind(coeff2$sgnum, coeff2$model, coeff2$aic, coeff2$R2MM, coeff2$estimate, coeff2$se, coeff2$estimate2, coeff2$se2))
colnames(models) <- c("sgnum", "model", "aic", "R2MM", "estimate", "se", "estimate2", "se2")

# select models for comparison
modelsall <- subset(models, sgnum != 0)
models_sort <- arrange(modelsall, sgnum, aic)
models_sort$sgnum <- as.numeric(as.character(models_sort$sgnum))
models_sort$model <- gsub("data\\$", "", models_sort$model)
best_model <- arrange(ddply(models_sort, .(sgnum), function(x)x[which.min(x$aic), ]), sgnum)
NULL_model <- subset(models_sort, model == "NULL")
NULL_model <- NULL_model[!duplicated(NULL_model$sgnum), ]
summer_model <- subset(models_sort, model == "njjtmp" | model == "njjtmx" | model == "njjtmn" | model == "njunetmp" | model == "njunetmx" | model == "njunetmn" | model == "njulytmp" | model == "njulytmx" | model == "njulytmn" | model == "NULL")
summer_model <- arrange(ddply(summer_model, .(sgnum), function(x)x[which.min(x$aic), ]), sgnum)

# calculate delta aic values
model_comp <- merge(best_model, NULL_model, by.x = "sgnum", by.y = "sgnum")
model_comp <- merge(model_comp, summer_model, by.x = "sgnum", by.y = "sgnum")
model_comp$estimate.y <- NULL
model_comp$se.y <- NULL
model_comp$model.y <- NULL
model_comp$estimate2.y <- NULL
model_comp$se2.y <- NULL
colnames(model_comp) <- c("sgnum", "nmodel_best", "naic_best", "R2MM_best", "nest_best", "nse_best", "nest2_best", "nse2_best", "naic_null", "R2MM_null", "nmodel_summer", "naic_summer", "R2MM_summer", "nest_summer", "nse_summer", "nest2_summer", "nse2_summer")
model_comp$naic_best <- as.numeric(as.character(model_comp$naic_best))
model_comp$naic_null <- as.numeric(as.character(model_comp$naic_null))
model_comp$naic_summer <- as.numeric(as.character(model_comp$naic_summer))
model_comp <- cbind(model_comp, (model_comp$naic_null-model_comp$naic_best), (model_comp$naic_null-model_comp$naic_summer))
colnames(model_comp) <- c("sgnum", "nmodel_best", "naic_best", "R2MM_best", "nest_best", "nse_best", "nest2_best", "nse2_best", "naic_null", "R2MM_null", "nmodel_summer", "naic_summer", "R2MM_summer", "nest_summer", "nse_summer", "nest2_summer", "nse2_summer", "ndaic", "ndaic_summer")

# If model comparison is less than 2 set delta aic to 0
model_comp$ndaic[model_comp$ndaic < 2] <- 0
model_comp$ndaic_summer[model_comp$ndaic_summer < 2] <- 0
model_comp <- merge(siteinfosum, model_comp, by.x = "sgnum", by.y = "sgnum")

model_comp_time_slope[y:(y+length(model_comp$sgnum)-1), 1] <- t1
model_comp_time_slope[y:(y+length(model_comp$sgnum)-1), 2:4] <- cbind(as.character(model_comp$sgnum), as.character(model_comp$site), as.character(model_comp$site_code))
model_comp_time_slope[y:(y+length(model_comp$sgnum)-1), 5:6] <- cbind(as.numeric(as.character(model_comp$lat)), as.numeric(as.character(model_comp$lon)))
model_comp_time_slope[y:(y+length(model_comp$sgnum)-1), 7] <- as.character(model_comp$nmodel_best)
model_comp_time_slope[y:(y+length(model_comp$sgnum)-1), 8:15] <- cbind(as.numeric(as.character(model_comp$naic_best)), as.numeric(as.character(model_comp$R2MM_best)), as.numeric(as.character(model_comp$nest_best)), as.numeric(as.character(model_comp$nse_best)), as.numeric(as.character(model_comp$nest2_best)), as.numeric(as.character(model_comp$nse2_best)), as.numeric(as.character(model_comp$naic_null)), as.numeric(as.character(model_comp$R2MM_null)))
model_comp_time_slope[y:(y+length(model_comp$sgnum)-1), 16] <- as.character(model_comp$nmodel_summer)
model_comp_time_slope[y:(y+length(model_comp$sgnum)-1), 17:24] <- cbind(as.numeric(as.character(model_comp$naic_summer)), as.numeric(as.character(model_comp$R2MM_summer)), as.numeric(as.character(model_comp$nest_summer)), as.numeric(as.character(model_comp$nse_summer)), as.numeric(as.character(model_comp$nest2_summer)), as.numeric(as.character(model_comp$nse2_summer)), as.numeric(as.character(model_comp$ndaic)), as.numeric(as.character(model_comp$ndaic_summer)), as.numeric(as.character(model_comp$nprop)), as.numeric(as.character(model_comp$shrubn)))
colnames(model_comp_time_slope) <- c("timewindow", names(model_comp))
x <- x + 1
y <- y + length(model_comp$sgnum)
}

save(model_comp_time_slope, file="scripts/users/imyerssmith/shrub_synthesis/model_comp_time_slope.RData")
write.csv(model_comp_time_slope, file = "scripts/users/imyerssmith/shrub_synthesis/model_comp_time_slope.csv", row.names=FALSE)

