######################################################################
# white paper, this code, and data available at                      #
# title: "JUMP SMOOTHING ALGORITHM"                                  #
# https://github.com/tsiappoutas/smoothing                           #
# author: michael tsiappoutas                                        #
# email: tsiappoutas@hotmail.com; tsiappoutas@gmail.com              #
# last update: Oct. 2020                                             #
######################################################################

library(rgdal)
library(spdep)
library(dplyr)
library(choroplethr)
library(choroplethrMaps)
# install choroplethrZip
library(devtools)
install_github('arilamstein/choroplethrZip@v1.5.0')
library(choroplethrZip)
library(ggplot2)

options(scipen=6, digits=6)
memory.limit(size=40000)
setwd("C:/Documents/R/smoothing/data")

# CREATE nbr.pairs.csv FROM
# source: ftp://ftp2.census.gov/geo/tiger/TIGER2013/ZCTA5/tl_2013_us_zcta510.zip
# download and extract tl_2013_us_zcta510.shp on your C: drive
# this files is not included on GitHub because it's too big,
# but you don't really need it for this code, because nbs.pairs.csv 
# is on GitHub. I left the code here in case you want to download a more
# recent version of the shape file and redo nbs.pairs.csv.

shp <- st_read("C:/.../smoothing/data/tl_2013_us_zcta510.shp")
# identify neighbours for each poly
nbs <- setNames(poly2nb(shp), shp$ZCTA5CE10)
# convert to a binary neighbour matrix
nbs.mat <- nb2mat(nbs, zero.policy=TRUE, style='B')
# assign zip codes as dimension names
dimnames(nbs.mat) <- list(shp$ZCTA5CE10, shp$ZCTA5CE10)
nbs.list <- sapply(row.names(nbs.mat), function(x) names(which(nbs.mat[x, ] == 1)))
nbs.pairs <- data.frame(zipcode=rep(names(nbs.list), sapply(nbs.list, length)), 
                        neighbour=unlist(nbs.list))
# we save this once and we will load it when needed (it takes time to process)
write.csv(nbs.pairs, "C:/.../smoothing/data/nbs.pairs.csv", row.names=FALSE)

# CREATE DATA TO BE SMOOTHED
# source: https://catalog.data.gov/dataset/2010-census-populations-by-zip-code
data <- read.csv('2010_Census_Populations_by_Zip_Code.csv', header=TRUE)
data <- data[, c('Zip.Code', 'Total.Population', 'Median.Age', 'Average.Household.Size')]
names(data) <- c('zip', 'population', 'age', 'hh_size')
data$population <- ifelse(data$population < 500, median(data$population), data$population)
data$age <- ifelse(data$age == 0, median(data$age), data$age)
data$hh_size <- ifelse(data$hh_size == 0, median(data$hh_size), data$hh_size)
data$exposure <- data$population / max(data$population)
data$age_relativity <- data$age / mean(data$age)
data$hh_size_relativity <- data$hh_size / mean(data$hh_size)


# IMPORT ZIP NEIGHBOR PAIRS
nbs.pairs <- read.csv('nbs.pairs.csv', header=TRUE)

# FIRST SMOOTHING ITERATION
# merge data to be smoothed
pairs <- merge(nbs.pairs, data, by.x='zipcode', by.y='zip')
pairs <- pairs[, c("zipcode", "neighbour", "exposure", "age_relativity")]
names(pairs) <- c('zip', 'nbr', 'zip.exp', 'zip.rel')
# so far we have the exposure and relativity of the main zip code.
# we need the exposure and relativity for the neighboring zip codes as well, so
# merge for exposure and zip for neighbor zip code
pairs <- merge(pairs, data, by.x='nbr', by.y='zip')
pairs <- pairs[order(pairs$zip),]
# order variables, keep, rename
names(pairs)[8:9] <- c('nbr.exp', 'nbr.rel')
pairs <- pairs[, c('zip', 'nbr', 'zip.exp', 'zip.rel', 'nbr.exp', 'nbr.rel')]

# JUMP
# if the difference (jump) between zip.rel with the nbr.rel is small
# (less than a threshold), then use the zip.rel in the nbr.rel column
# and call it nbr.rel.jump
pairs1 <- pairs
threshold <- 0.095
pairs1$jump <- ifelse(abs(pairs$zip.rel - pairs$nbr.rel) <= threshold, 1, 0)
prop.table(table(pairs1$jump))
pairs1$nbr.rel.jump <- ifelse(pairs1$jump == 1, pairs1$zip.rel, pairs1$nbr.rel)

# CALCULATE WEIGHTED AVERAGES
wtd.avgs1 <- pairs1 %>%
  na.omit() %>%
  group_by(zip) %>%
  summarise(
    zip.rel               = round(mean(zip.rel, na.rm=TRUE), 4),
    zip.exp               = round(mean(zip.exp, na.rm=TRUE), 4),
    wtd.nbr.rel           = round(weighted.mean(nbr.rel.jump, nbr.exp, na.rm=TRUE), 4),
    avg.nbr.exp           = round(mean(nbr.exp, na.rm=TRUE), 4))
wtd.avgs1$smooth.rel1 <- wtd.avgs1$zip.rel * (wtd.avgs1$zip.exp / (wtd.avgs1$zip.exp + wtd.avgs1$avg.nbr.exp)) + wtd.avgs1$wtd.nbr.rel * (wtd.avgs1$avg.nbr.exp / (wtd.avgs1$zip.exp + wtd.avgs1$avg.nbr.exp))

# GOVERNMENT ZIPS
# the dataset below comes from library(choroplethrZip)
data(df_pop_zip)
df_pop_zip$zip.num <- as.numeric(df_pop_zip$region)
map.all <- merge(wtd.avgs1, df_pop_zip, by.x='zip', by.y='zip.num', all.x=TRUE)  

# MAP
map.rel <- map.all[, c('region', 'zip.rel')]
names(map.rel) <- c('region', 'value')

# below we create the actual map. I'm using 6037 as my county_zoom
# because you need to use the exact number you find in this data(zip.regions)
# under county.fips.numeric. But check out the zip_choropleth documentation.
# You can do selected zip codes, states, nationwide, metropolitan ares, etc.

zip_choropleth(map.rel,
               county_zoom = '6037',
               title      = 'Unsmoothed Relativities',
               legend     = 'Relativities')

# Now replace 'value' for the smooth relativity, and plot again
map <- map.all[, c('region', 'smooth.rel1')]
names(map) <- c('region', 'value')

zip_choropleth(map,
               county_zoom = '6037',
               title      = 'Smooth Relativities Iteration 1, jump=0.4, 6% jumped',
               legend     = 'Relativities')


# SECOND SMOOTHING ITERATION
# we merge the smoothed zips from the previous iteration, along with the main zips 
# and neighbors, then smooth again
pairs2 <- merge(pairs, wtd.avgs1[, c('zip', 'smooth.rel1')], by.x='zip', by.y='zip')
names(pairs2)[7] <- 'smooth.zip1' 
# merge neighbor zip with relativites and exposures
pairs2 <- merge(pairs2, wtd.avgs1[, c('zip', 'smooth.rel1')], by.x='nbr', by.y='zip')
names(pairs2)[8] <- 'smooth.nbr1'
pairs2 <- pairs2[order(pairs2$zip),]
pairs2 <- pairs2[, c('zip', 'nbr', 'zip.exp', 'nbr.exp', 'zip.rel', 'nbr.rel', 'smooth.zip1', 'smooth.nbr1')]

# redefine jump threshold to achieve about 50% jump rate
threshold <- 0.062
pairs2$jump <- ifelse(abs(pairs2$smooth.zip1 - pairs2$smooth.nbr1) <= threshold, 1, 0)
prop.table(table(pairs2$jump))
pairs2$nbr.rel.jump <- ifelse(pairs2$jump == 1, pairs2$smooth.zip1, pairs2$smooth.nbr1)

# CALCULATE WEIGHTED AVERAGES
wtd.avgs2 <- pairs2 %>%
  na.omit() %>%
  group_by(zip) %>%
  summarise(
    zip.rel               = round(mean(zip.rel, na.rm=TRUE), 4),
    zip.exp               = round(mean(zip.exp, na.rm=TRUE), 4),
    wtd.nbr.rel           = round(weighted.mean(nbr.rel.jump, nbr.exp, na.rm=TRUE), 4),
    avg.nbr.exp           = round(mean(nbr.exp, na.rm=TRUE), 4))
wtd.avgs2$smooth.rel2 <- wtd.avgs2$zip.rel * (wtd.avgs2$zip.exp / (wtd.avgs2$zip.exp + wtd.avgs2$avg.nbr.exp)) + wtd.avgs2$wtd.nbr.rel * (wtd.avgs2$avg.nbr.exp / (wtd.avgs2$zip.exp + wtd.avgs2$avg.nbr.exp))

# GOVERNMENT ZIPS
map.all <- merge(wtd.avgs2, df_pop_zip, by.x='zip', by.y='zip.num', all.x=TRUE)  

# MAP
map.rel <- map.all[, c('region', 'smooth.rel2')]
names(map.rel) <- c('region', 'value')
zip_choropleth(map.rel,
               county_zoom = '6037',
               title      = 'Second Iteration',
               legend     = 'Relativities')


# THIRD SMOOTHING ITERATION
# smooth the second iteration smooth relativities
pairs3 <- merge(pairs, wtd.avgs2[, c('zip', 'smooth.rel2')], by.x='zip', by.y='zip')
names(pairs3)[7] <- 'smooth.zip2' 
# merge neighbor zip with relativites and exposures
pairs3 <- merge(pairs3, wtd.avgs2[, c('zip', 'smooth.rel2')], by.x='nbr', by.y='zip')
names(pairs3)[8] <- 'smooth.nbr2'
pairs3 <- pairs3[order(pairs3$zip),]
pairs3 <- pairs3[, c('zip', 'nbr', 'zip.exp', 'nbr.exp', 'zip.rel', 'nbr.rel', 'smooth.zip2', 'smooth.nbr2')]

# redefine jump threshold to achieve about 50% jump rate
threshold <- 0.056
pairs3$jump <- ifelse(abs(pairs3$smooth.zip2 - pairs3$smooth.nbr2) <= threshold, 1, 0)
prop.table(table(pairs3$jump))
pairs3$nbr.rel.jump <- ifelse(pairs3$jump == 1, pairs3$smooth.zip2, pairs3$smooth.nbr2)

# CALCULATE WEIGHTED AVERAGES
wtd.avgs3 <- pairs3 %>%
  na.omit() %>%
  group_by(zip) %>%
  summarise(
    zip.rel               = round(mean(zip.rel, na.rm=TRUE), 4),
    zip.exp               = round(mean(zip.exp, na.rm=TRUE), 4),
    wtd.nbr.rel           = round(weighted.mean(nbr.rel.jump, nbr.exp, na.rm=TRUE), 4),
    avg.nbr.exp           = round(mean(nbr.exp, na.rm=TRUE), 4))
wtd.avgs3$smooth.rel3 <- wtd.avgs3$zip.rel * (wtd.avgs3$zip.exp / (wtd.avgs3$zip.exp + wtd.avgs3$avg.nbr.exp)) + wtd.avgs3$wtd.nbr.rel * (wtd.avgs3$avg.nbr.exp / (wtd.avgs3$zip.exp + wtd.avgs3$avg.nbr.exp))

# GOVERNMENT ZIPS
map.all <- merge(wtd.avgs3, df_pop_zip, by.x='zip', by.y='zip.num', all.x=TRUE)  

# MAP
map.rel <- map.all[, c('region', 'smooth.rel3')]
names(map.rel) <- c('region', 'value')
zip_choropleth(map.rel,
               county_zoom = '6037',
               title      = 'Third Iteration',
               legend     = 'Relativities')


# JUMP COEFFICIENT PLOT
# 4th iteration
pairs4 <- merge(pairs, wtd.avgs3[, c('zip', 'smooth.rel3')], by.x='zip', by.y='zip')
names(pairs4)[7] <- 'smooth.zip3' 
# merge neighbor zip with relativites and exposures
pairs4 <- merge(pairs4, wtd.avgs3[, c('zip', 'smooth.rel3')], by.x='nbr', by.y='zip')
names(pairs4)[8] <- 'smooth.nbr3'
pairs4 <- pairs4[order(pairs4$zip),]
pairs4 <- pairs4[, c('zip', 'nbr', 'zip.exp', 'nbr.exp', 'zip.rel', 'nbr.rel', 'smooth.zip3', 'smooth.nbr3')]
# redefine jump threshold to achieve about 50% jump rate
threshold <- 0.057
pairs4$jump <- ifelse(abs(pairs4$smooth.zip3 - pairs4$smooth.nbr3) <= threshold, 1, 0)
prop.table(table(pairs4$jump))
pairs4$nbr.rel.jump <- ifelse(pairs4$jump == 1, pairs4$smooth.zip3, pairs4$smooth.nbr3)

# CALCULATE WEIGHTED AVERAGES
wtd.avgs4 <- pairs4 %>%
  na.omit() %>%
  group_by(zip) %>%
  summarise(
    zip.rel               = round(mean(zip.rel, na.rm=TRUE), 4),
    zip.exp               = round(mean(zip.exp, na.rm=TRUE), 4),
    wtd.nbr.rel           = round(weighted.mean(nbr.rel.jump, nbr.exp, na.rm=TRUE), 4),
    avg.nbr.exp           = round(mean(nbr.exp, na.rm=TRUE), 4))
wtd.avgs4$smooth.rel4 <- wtd.avgs4$zip.rel * (wtd.avgs4$zip.exp / (wtd.avgs4$zip.exp + wtd.avgs4$avg.nbr.exp)) + wtd.avgs4$wtd.nbr.rel * (wtd.avgs4$avg.nbr.exp / (wtd.avgs4$zip.exp + wtd.avgs4$avg.nbr.exp))


# 5th iteration
pairs5 <- merge(pairs, wtd.avgs4[, c('zip', 'smooth.rel4')], by.x='zip', by.y='zip')
names(pairs5)[7] <- 'smooth.zip4' 
# merge neighbor zip with relativites and exposures
pairs5 <- merge(pairs5, wtd.avgs4[, c('zip', 'smooth.rel4')], by.x='nbr', by.y='zip')
names(pairs5)[8] <- 'smooth.nbr4'
pairs5 <- pairs5[order(pairs5$zip),]
pairs5 <- pairs5[, c('zip', 'nbr', 'zip.exp', 'nbr.exp', 'zip.rel', 'nbr.rel', 'smooth.zip4', 'smooth.nbr4')]
# redefine jump threshold to achieve about 50% jump rate
threshold <- 0.057
pairs5$jump <- ifelse(abs(pairs5$smooth.zip4 - pairs5$smooth.nbr4) <= threshold, 1, 0)
prop.table(table(pairs5$jump))

####################################################################
######################### END OF CODE ##############################
