# Shrub growth synthesis time window figures script:
# This code creates the time window figure for the shrub synthesis manuscript.
# Written by: Isla Myers-Smith (e-mail: isla.myers-smith@ed.ac.uk)
# 10 April 2015

model_comp_time <- read.csv("scripts/users/imyerssmith/shrub_synthesis/model_comp_time.csv", 
                            header=TRUE, stringsAsFactors=FALSE)
model_comp_time_slope <- read.csv("scripts/users/imyerssmith/shrub_synthesis/model_comp_time_slope.csv", 
                                  header=TRUE, stringsAsFactors=FALSE)
samplesize <- read.csv("scripts/users/imyerssmith/shrub_synthesis/samplesize.csv", 
                                  header=TRUE, stringsAsFactors=FALSE)
# Time window analysis
model_comp_time_mean <- model_comp_time
model_comp_time_mean$decade <- signif(model_comp_time$timewindow, digits = 3)
model_comp_time_mean <- subset(model_comp_time_mean, decade != 0)
mctm <- aggregate(cbind(model_comp_time_mean$ndaic, model_comp_time_mean$R2MM_best, model_comp_time_mean$nprop), by=list(model_comp_time_mean$decade, model_comp_time_mean$sgnum), FUN=mean, na.rm=TRUE)
colnames(mctm) <- c("decade", "sgnum", "mean_ndaic", "mean_R2MM_best", "mean_nprop")

# Time window analysis - slopes
model_comp_time_slope_mean <- model_comp_time_slope
model_comp_time_slope_mean$decade <- signif(model_comp_time_slope$timewindow, digits = 3)
model_comp_time_slope_mean <- subset(model_comp_time_slope_mean, decade != 0)
model_comp_time_slope_mean[model_comp_time_slope_mean=="NaN"] <- "NA"
mctm_slope <- aggregate(cbind(model_comp_time_slope_mean$nest_best, model_comp_time_slope_mean$nest2_best, model_comp_time_slope_mean$nest_summer, model_comp_time_slope_mean$nest2_summer), by=list(model_comp_time_slope_mean$decade, model_comp_time_slope_mean$sgnum), FUN=mean, na.rm=TRUE)
colnames(mctm_slope) <- c("decade", "sgnum", "mean_nest_best", "mean_nest2_best", "mean_nest_summer", "mean_nest2_summer")
mctm_slope <- merge(mctm_slope, model_comp_time_slope_mean)

siteinfo <- read.csv("scripts/users/imyerssmith/shrub_synthesis/site_info.csv", 
                     header=TRUE, stringsAsFactors=FALSE)
mctm_slope <- merge(siteinfo, mctm_slope)
mctm_slope$sitesp <- paste(mctm_slope$site_code,mctm_slope$genus_code, sep=":")
mctm_slope$sitesp <- factor(mctm_slope$sitesp, c("1AA:AVC","1AA:BETspp","1AA:SALspp","3HE:SALspp","4KL:SALspp","5DH:SALspp","6NL:AVC","7BI:CTE","8BN:BGL","9RE:SAR","10AH:CTE","11TL:CTE","12HC:CTE","13CB:CTE","14AX:SAR","15AB:CTE","16AD:CTE","17AL:CTE","18LA:SAR","19WB:BGL","20GQ:BGL","21NA:AVC","21NA:BGL","22AF:AVC","23WG:SALspp","24ZA:SAR","25ZF:SAR","26NF:SALspp","27MI:Betpub","28ST:VACspp","29CH:PMU","30SB:CTE","31SC:SPO","32AS:JNA","36VR:SRI","37LB:AVC","37LB:SRI","38YR:SRI","39YA:PMU","40KY:BNA","40KY:SPU"))
mctm <- merge(mctm, mctm_slope)

# Sample size
samplesize_mean <- samplesize
samplesize_mean$decade <- signif(samplesize_mean$timewindow, digits = 3)
samplesize_mean <- subset(samplesize_mean, decade != 0)
samplesize_mean <- aggregate(samplesize_mean$nshrub, by=list(samplesize_mean$decade, samplesize_mean$sgnum), FUN=mean, na.rm=TRUE)
colnames(samplesize_mean) <- c("decade", "sgnum", "mean_nshrub")
mctm <- merge(mctm, samplesize_mean)
mctm_decade <- aggregate(cbind(mctm$mean_nshrub, mctm$mean_ndaic, mctm$mean_R2MM_best, mctm$mean_nest_best, mctm$mean_nprop), by=list(mctm$decade, mctm$sgnum, mctm$sitesp), FUN=mean, na.rm=TRUE)
colnames(mctm_decade) <- c("decade", "sgnum", "sitesp", "mean_nshrub", "mean_ndaic", "mean_R2MM_best", "mean_nest_best", "mean_nprop")

# Plot time window analysis
plot <- qplot(decade+10, mean_R2MM_best, data = mctm_decade, geom = c("point","line"), alpha = I(5/10), xlim = c(1960, 1990), colour = sitesp) + 
  theme_classic(base_size = 12, base_family = "Helvetica") + 
  theme(axis.title.x = element_text(vjust=-0.35), axis.title.y = element_text(vjust=1.75)) +
  labs(x = "Midpoint of the 20-year time window", y = expression(paste("Mixed ","Model ",R^2))) +
  guides(col = guide_legend(ncol = 3)) + 
  theme(legend.title=element_blank()) +
  scale_fill_brewer(palette="RdYlGn")

ssizeplot <- qplot(decade+10, mean_nshrub, data = mctm_decade, geom = c("point","line"), alpha = I(5/10), xlim = c(1960, 1990), colour = sitesp) + 
  theme_classic(base_size = 12, base_family = "Helvetica") + 
  theme(axis.title.x = element_text(vjust=-0.35), axis.title.y = element_text(vjust=2), axis.ticks = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(0.875,1,-0.375,1.1),"cm")) +
  labs(x = "", y = "Number of shrubs") +
  theme(legend.position="none") +
  scale_fill_brewer(palette="RdYlGn")

daicplot <- qplot(decade+10, mean_ndaic, data = mctm_decade, geom = c("point","line"), alpha = I(5/10), xlim = c(1960, 1990), colour = sitesp) + 
  theme_classic(base_size = 12, base_family = "Helvetica") + 
  theme(axis.title.x = element_text(vjust=-0.35), axis.title.y = element_text(vjust=2.5), axis.ticks = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(0.5,1,0,1.3),"cm")) +
  labs(x = "", y = expression(paste(Delta," AIC"))) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="RdYlGn")

R2plot <- qplot(decade+10, mean_R2MM_best, data = mctm_decade, geom = c("point","line"), alpha = I(5/10), xlim = c(1960, 1990), colour = sitesp) + 
  theme_classic(base_size = 12, base_family = "Helvetica") + 
  theme(axis.title.x = element_text(vjust=-0.35), axis.title.y = element_text(vjust=2.15), axis.ticks = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(0.125,1,0.375,1.05),"cm")) +
  labs(x = "", y = expression(paste("Mixed ","Model ",R^2))) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="RdYlGn")

slopeplot <- qplot(decade+10, mean_nest_best, data = mctm_decade, geom = c("point","line"), alpha = I(5/10), xlim = c(1960, 1990), colour = sitesp) + 
  theme_classic(base_size = 12, base_family = "Helvetica") + 
  theme(axis.title.x = element_text(vjust=-0.35), axis.title.y = element_text(vjust=2.15), axis.ticks = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(-0.25,1,0.75,1.05),"cm")) +
  labs(x = "", y = "Slope") +
  theme(legend.position="none") +
  scale_fill_brewer(palette="RdYlGn")

propplot <- qplot(decade+10, mean_nprop, data = mctm_decade, geom = c("point","line"), alpha = I(5/10), xlim = c(1960, 1990), colour = sitesp) + 
  theme_classic(base_size = 12, base_family = "Helvetica") + 
  theme(axis.title.x = element_text(vjust=-0.4), axis.title.y = element_text(vjust=1.7), plot.margin = unit(c(-0.625,1,0.875,0.95),"cm")) +
  labs(x = "Midpoint of the 20-year time window", y = "Prop. temp. sens.") +
  theme(legend.position="none") +
  scale_fill_brewer(palette="RdYlGn")

#extract legend
g_legend<-function(p){
  tmp <- ggplotGrob(p)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

legend <- g_legend(plot)

blank <- rectGrob(gp=gpar(col=NA))

timewindowfigure <- arrangeGrob(ssizeplot, blank, daicplot, blank, R2plot, blank, slopeplot, g_legend(plot), propplot, ncol=2, nrow=5)
print(timewindowfigure)

ggsave(file="scripts/users/imyerssmith/shrub_synthesis/FigS13timewindowfigure.png", timewindowfigure, width=8, height=10)

