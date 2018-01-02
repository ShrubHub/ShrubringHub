# Shrub growth synthesis figures script:
# This code creates all the figures and supplementary analyses for the shrub synthesis manuscript.
# Written by: Isla Myers-Smith (e-mail: isla.myers-smith@ed.ac.uk)
# 26 June 2014

# Libraries ---------------------------------------------------------------

library(maps)
library(mapdata)
library(mapproj)
library(rgdal)
library(maptools)
library(raster)
library(sp)
library(scales)

# Load data ---------------------------------------------------------------

modelsumdata <- read.csv("scripts/users/imyerssmith/shrub_synthesis/modelsumdata.csv", header=TRUE, stringsAsFactors=FALSE)

# Map figure --------------------------------------------------------------
  
  tiff(file="scripts/users/imyerssmith/shrub_synthesis/testmap.tiff", width=2000, height=2300, compression = "none")
  par(mfrow=c(1,1), mar=c(0, 0, 0, 0), oma=c(2, 2, 2, 2), mgp=c(5, 2, 0))
  
  #plot map
  wrld <- map(plot=FALSE, interior=FALSE, wrap=TRUE, ylim=c(45, 90), xlim=c(-180, 180))
  wrld_sp <- map2SpatialLines(wrld)
  proj4string(wrld_sp) <- CRS("+proj=longlat")
  laea_wrld_sp <- spTransform(wrld_sp, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
  plot(laea_wrld_sp, col="black")
  
  data <- modelsumdata
  data1 <- subset(data, (ndaic_summer > 0) & !is.na(as.numeric(as.character(nest_summer))))
  data2 <- subset(data, (ndaic_summer == 0) & ndaic > 0 & (is.na(as.numeric(as.character(data$nest_summer)))))
  data3 <- subset(data, (ndaic_summer == 0) & ndaic > 0 & (!is.na(as.numeric(as.character(data$nest_summer)))))
  data4 <- subset(data, ndaic == 0 | is.na(as.numeric(as.character(R2MM_best))))
  
coord1 <- SpatialPoints(data.frame(data1$lon, data1$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord1_sp <- spTransform(coord1, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))

coord2 <- SpatialPoints(data.frame(data2$lon, data2$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord2_sp <- spTransform(coord2, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))

coord3 <- SpatialPoints(data.frame(data3$lon, data3$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord3_sp <- spTransform(coord3, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
#plot(laea_coord3_sp, add = TRUE, pch=1, col="#FFFFFF", cex=10, lwd=10)

plot(laea_coord3_sp, add = TRUE, pch=1, col="#000000", cex=10, lwd=10)
plot(laea_coord2_sp, add = TRUE, pch=1, col="#000000", cex=10, lwd=10)
plot(laea_coord1_sp, add = TRUE, pch=1, col="#000000", cex=10, lwd=10)

coord4 <- SpatialPoints(data.frame(data4$lon, data4$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
laea_coord4_sp <- spTransform(coord4, CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m"))
plot(laea_coord4_sp, add = TRUE, pch=4, col="#000000", cex=5, lwd=10)

leg.txt <- c("No clim. sens.", "Climate sensitive")
leg.blank <- c("", "")
legend("bottomleft", legend=leg.blank, x.intersp=1.5, y.intersp=2, col=c("#000000", "#FFFFFF00"), pch=c(4, 1), pt.cex=c(5, 10), bty="n", cex = 3, pt.lwd=10)
legend("bottomleft", legend=leg.txt, x.intersp=1.5, y.intersp=2, col="#FFFFFF00", pt.cex=10, bty="n", cex = 3)
legend("bottomleft", legend=leg.blank, x.intersp=1.5, y.intersp=2, col=c("#FFFFFF00", "#000000"), pch=c(1), pt.cex=10, bty="n", cex = 3, pt.lwd=10, lty=NULL)

dev.off()
