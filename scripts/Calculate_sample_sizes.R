data1 <- cbind.data.frame(data$snum, data$sgnum, data$shrub_num, data$year, data$rw)
names(data1) <- c("snum", "sgnum", "shrub_num", "year", "rw")

data2 <- unique(cbind.data.frame(data$snum, data$sgnum, data$shrub_num, data$year))
names(data2) <- c("snum", "sgnum", "shrub_num", "year")

head(data2)

round((table(data1$sgnum, round(data1$year, digits = -1))/10), digits = 0)
