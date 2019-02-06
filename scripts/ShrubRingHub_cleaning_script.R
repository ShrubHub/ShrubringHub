##############################################
### ShrubHub ring width time series        ###
### Cleaning and correction of dataset     ###
### Sandra Angers-Blondin                  ###
### 14-03-2018                             ###
##############################################


# Libraries -----------------------------------------------------------------------------------

library(dplyr)     # a data manipulation package
library(tidyr)     # to reshape dataframes
library(ggplot2)   # to make pretty plots


# Load data -----------------------------------------------------------------------------------

load("data/shrubring_growth_full.RData")  # original Shrub Hub growth data

data <- growth.full

# Cleaning: heights ---------------------------------------------------------------------------

data[data == "NaN"] <- "NA"

# I am creating a qualitative height column keeping imprecise measurements such as < 1.5
data = data %>% mutate(height.raw = height)

# Then overwriting height as a numeric-only column

data$height <- gsub("[^0-9.\\<>]", "", data$height) ## removes the "m" symbol and extra spaces in the df, but leaves the "greater than" and "smaller than" symbol so the precise measurements can be separated from the imprecise ones
data$height <- as.numeric(data$height) # keeps only the numeric values, those with <> symbols become NAs

# Correct height units - we want height in meters

data <- mutate(data, height = ifelse(height > 6, height/100, height))



# Cleaning: width -----------------------------------------------------------------------------

# Check width distribution - clearly some are in cm and some in m
ggplot(distinct(data, site_full, sgnum, genus, shrub_num, width, width2), aes(x = as.factor(as.character(sgnum)), y = width)) +
  geom_boxplot(aes(colour = genus), alpha = 0.2, outlier.shape=NA, width = 0.4) +
  geom_jitter(aes(colour = genus), shape = 16, alpha = 0.5, size = 1, width = 0.2) +
  labs(x = "", y = "Width (m)") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust=0.5, hjust = 1))

# Correct units to meters for sgnum 8, 10, and 33

data <- data %>% mutate(width = ifelse(sgnum %in% c(8, 10, 33), width/100, width),
                        width2 = ifelse(sgnum %in% c(8, 10, 33), width2/100, width2))


# Cleaning: Stem Width ------------------------------------------------------------------------

# Check actual measurements - probably some in mm, some in cm
ggplot(distinct(data, site_full, sgnum, genus, shrub_num, stem), aes(x = as.factor(as.character(sgnum)), y = stem)) +
  geom_boxplot(aes(colour = genus), alpha = 0.2, outlier.shape=NA, width = 0.4) +
  geom_jitter(aes(colour = genus), shape = 16, alpha = 0.5, size = 1, width = 0.2) +
  labs(x = "", y = "Width") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust=0.5, hjust = 1))

# Correct units from cm to mm for sgnum 8, 10, 11, and 24

data <- data %>% mutate(stem = ifelse(sgnum %in% c(8, 10, 11, 24), stem*10, stem))

# Check if that makes sense given shrub height

ggplot(distinct(data, sgnum, genus, shrub_num, height, stem), aes(x = height, y = stem)) +
  geom_point(aes(colour = genus), alpha = 0.2) +
  labs(x = "Height (m)", y = "Stem width (mm)") +
  theme_classic() +
  facet_wrap(~as.factor(as.character(sgnum)), scales = "fixed", drop = TRUE)


# Cleaning: raw ring width --------------------------------------------------------------------

# Check actual measurements - probably some in mm, some in um
ggplot(distinct(data, site_full, sgnum, genus, shrub_num, rw), aes(x = as.factor(as.character(sgnum)), y = rw)) +
  geom_boxplot(aes(colour = genus), alpha = 0.2, outlier.shape=NA, width = 0.4) +
  geom_jitter(aes(colour = genus), shape = 16, alpha = 0.5, size = 1, width = 0.2) +
  labs(x = "", y = "Width") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust=0.5, hjust = 1))

method <- data %>% group_by(sgnum) %>% filter(length(unique(sampling_method)) == 1) %>% summarise(method = unique(sampling_method)) ## this loses the 3 sites for which there were more than 1 entry but fine just to check
# the missing 3:
unique(data[data$sgnum == "35",]$sampling_method)  # what's PR? Ramet!
unique(data[data$sgnum == "13",]$sampling_method)  # ring width OK
unique(data[data$sgnum == "45",]$sampling_method)  # ring width OK

## The following sgnums are not shrub rings:
### 12,15, 16, 17, 18,20, 21, 22, 34, 35, 37      mm all make sense for these


## The problematic sgnums are:
### 1, 3, 11, 26 ### divide by 10

# Correct the errors by dividing by 10 (two identical columns, value and rw)

data <- data %>% mutate(rw = ifelse(sgnum %in% c(1, 3, 11, 26), rw/10, rw),
                        value = ifelse(sgnum %in% c(1, 3, 11, 26), value/10, value))

# Create an approximate stem width ------------------------------------------------------------

# Overwriting the stemw column which had some errors 
data <- data %>% group_by(sgnum, shrub_num) %>% 
  mutate(stemw = ifelse(sampling_method %in% c("LS", "R ", "RC", "LS,RC,SS", "LS, RC, SI"),
                             2*(sum(rw) + 5*mean(rw)),
                             NA)
  )



# Check: looks ok in mm
ggplot(distinct(data, site_full, sgnum, genus, shrub_num, stemw), aes(x = as.factor(as.character(sgnum)), y = stemw)) +
  geom_boxplot(aes(colour = genus), alpha = 0.2, outlier.shape=NA, width = 0.4) +
  geom_jitter(aes(colour = genus), shape = 16, alpha = 0.5, size = 1, width = 0.2) +
  labs(x = "", y = "Stem width") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust=0.5, hjust = 1))


# Visualisation -------------------------------------------------------------------------------

# Check that heights are correct
ggplot(distinct(data, site_full, sgnum, genus, shrub_num, height), aes(x = as.factor(as.character(sgnum)), y = height)) +
  geom_boxplot(aes(colour = genus), alpha = 0.2, outlier.shape=NA, width = 0.4) +
  geom_jitter(aes(colour = genus), shape = 16, alpha = 0.5, size = 1, width = 0.2) +
  labs(x = "", y = "Height (m)") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust=0.5, hjust = 1))


# Check that widths are correct
ggplot(distinct(data, site_full, sgnum, genus, shrub_num, width, width2), aes(x = as.factor(as.character(sgnum)), y = width)) +
  geom_boxplot(aes(colour = genus), alpha = 0.2, outlier.shape=NA, width = 0.4) +
  geom_jitter(aes(colour = genus), shape = 16, alpha = 0.5, size = 1, width = 0.2) +
  labs(x = "", y = "Width (m)") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust=0.5, hjust = 1))

# Height vs width

ggplot(distinct(data, site_full, sgnum, genus, shrub_num, height, width, width2), aes(x = width, y = height)) +
  geom_point(aes(colour = genus), alpha = 0.2) +
  facet_wrap(~sgnum, drop = TRUE) +
  labs(x = "Width (m)", y = "Height (m)") +
  theme_classic()



# Save the new object -------------------------------------------------------------------------

growth.full.clean <- data %>% ungroup()  # remove internal grouping structure to avoid problems

save(growth.full.clean, file = "data/shrubring_growth_full_clean.RData")
