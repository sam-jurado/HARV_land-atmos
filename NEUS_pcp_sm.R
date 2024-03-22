install.packages('dplyr')
install.packages("oce")
install.packages("tidyr")
install.packages("lubridate")
install.packages("ncdf4")
install.packages("raster")
install.packages("rgdal")
install.packages("mapview")
install.packages("usmap")
install.packages("ggmap")
install.packages("mapdata")


library(rhdf5)
library(dplyr)
library(oce)
library(tidyr)
library(lubridate)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal)
library(sf)
library(usmap) #import the package
library(ggplot2) #use ggplot2 to add layer for visualization


######CLEANING NC FILES######

setwd("/Users/jurado/Harvard_Forest")
nc_data <- nc_open('soilm.mon.mean.nc')


lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat", verbose = F)
soilm <- ncvar_get(nc_data, "soilm")

dim(soilm) #check dimensions

fillvalue <- ncatt_get(nc_data, "soilm", "_FillValue") #check what the fill value is
soilm[soilm == fillvalue$value] <- NA #replace fill values with NA




# extract some specific ranges
soilm.slice <- soilm[,,1] 


LatIdx <- which( nc_data$var$lat$vals > 46126718 | nc_data$var$lat$vals < 4812618)
LonIdx <- which( nc_data$dim$x$vals > 630779 & nc_data$dim$x$vals < 830779)
soilm.slice <- ncvar_get(nc_data, "soilm")[ LonIdx, LatIdx,1]
nc_close(nc_data)



#soilm.slice <- soilm[,,1] 
r <- raster(t(soilm.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

r <- flip(r, direction='y')
plot(r)

#################
library(terra)

r <- rast('soilm.mon.mean.nc')
r
soilm <- r["soilm"]
soilm

e <- c(0, 11313351, 0, 8976020 ) |> ext()

soilm_crop <- crop(soilm, e)
soilm_crop

plot(soilm_crop)






####Setting up a square#####

# set the radius for the plots
radius <- 100000 # radius in meters (100 km)

#Northing and Easting of Petersham

northing <- 4712618
easting <- 730779



# define the plot edges based upon the plot radius. 

yPlus <- northing+radius
xPlus <- easting+radius
yMinus <- northing-radius
xMinus <- easting-radius
ID= "HF Study Site"

# calculate polygon coordinates for each plot centroid. 
square=cbind(xMinus,yPlus,  # NW corner
             xPlus, yPlus,  # NE corner
             xPlus,yMinus,  # SE corner
             xMinus,yMinus, # SW corner
             xMinus,yPlus)  # NW corner again - close ploygon

#Square Coordinates#




# create spatial polygons from coordinates
polys <- SpatialPolygons(mapply(function(poly, id) {
  xy <- matrix(poly, ncol=2, byrow=TRUE)
  Polygons(list(Polygon(xy)), ID=id)
}, 
split(square, row(square)), ID),
proj4string=CRS(as.character("+proj=utm +zone=11 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")))


#US Map
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)

usa <- map_data("usa") 
states <- map_data("state")
rect <- data.frame(
  x = c(-73.40645, -70.97355, -70.97355, -73.40645),
  y = c(43.42425, 43.42425, 41.6257,  41.6257)
)




north_east <- subset(states, region %in% c("connecticut", "maine", 
                                           "massachusetts", "new hampshire", 
                                           "new jersey", "new york", "pennsylvania",
                                           "rhode island","vermont"))

gg1 <- ggplot() + geom_polygon(data = north_east, aes(x=long, y = lat, group = group)) + 
  geom_polygon(data = north_east, aes(x = long, y = lat, group = group), color = "white") + 
  coord_fixed(1.3)

gg1+
  geom_polygon(data = rect, aes(x, y, group = 1),fill = NA, color = "red")+
  ggtitle("Northeastern United States Study Area")












####Using Harv Data#####

install.packages("anytime")   # Install anytime package
library(lubridate)
library(dplyr)
library("anytime")   

setwd("/Users/jurado/Harvard_Forest")
df_sm <- read.csv("hf206-03-HK-understory.csv")

setwd("/Users/jurado/")

df_prcp <- df_prcp<- read_csv("Downloads/dailyrain_19642023.csv")
#Monthly

df_sm$date <- format(as.POSIXct(df_sm$datetime), format = "%Y-%m-%d") 
dfMONTHLY <- df_sm %>% group_by(date) %>% summarise(swc1 = mean(swc1_ave, na.rm = TRUE),
                                                    swc2 = mean(swc2_ave, na.rm = TRUE),
                                                    swc3 = mean(swc3_ave, na.rm = TRUE),
                                                    swc4 = mean(swc4_ave, na.rm = TRUE))

#Filter


dfMONTHLY$swc_avg <- rowMeans(dfMONTHLY[ , c(2,3,4,5)], na.rm=TRUE)


dfMONTHLY$date <- anydate(dfMONTHLY$date) 
plot(dfMONTHLY$date,dfMONTHLY$swc_avg,type="l", xlab = "Time", ylab = "Soil Water Content")
abline(lm(dfMONTHLY$swc_avg ~ dfMONTHLY$date), lty = 2, lwd =2)
title(expression(paste(bold('Soil Water Content'))))
subtitle = "HEM 06/01/09 - 01/01/23, Soil Plots 1-4"
mtext(subtitle)






###### % of precip falling in heavy rain events#####
#####Heavy rain events are top 1% in time span ####
setwd("/Users/jurado/")
df_prcp <- df_prcp<- read_csv("Downloads/dailyrain_19642023.csv")
#Monthly


start_date <- as.Date("1964-01-01")  # Start date
end_date <- as.Date("2023-12-31")    # End date
list_of_dates <- seq(start_date, end_date, by = "day")
df_prcp$date <- list_of_dates



#filter for top 1% of each year, and check what amount of rain that constituted for the year#

#make a new column with NA instead of zero to only consider rain events

df_prcp["prec_mm"][df_prcp["prec_mm"] == 0] <- NA


#for loop to get and store the top 1% per year.

perc_hvy <- c()
num_hvy <- c()
sum_hvy <- c()



Q <- c()
for(x in 1964:2023){
  x <- toString(x)
  print(x)
  df_prcp <- df_prcp[df_prcp$date >= paste(x,"01-01",sep = "-") &   df_prcp$date <= paste(x,"12-31",sep = "-"), ]
  q <- quantile(df_prcp$prec_mm, prob=c(.90), type=1, na.rm =TRUE)
  Q <- append(Q,q)
  prcp_sum <- sum(df_prcp$prec_mm, na.rm = TRUE)
  df_prcp <- df_prcp %>% filter( prec_mm >= q)
  
  prcp_H_sum <- sum(df_prcp$prec_mm, na.rm = TRUE)
  prcp_num <- length(df_prcp$prec_mm)
  
  percent_heavy <- prcp_H_sum/prcp_sum*100
  
  perc_hvy <- append(perc_hvy,percent_heavy)
  num_hvy <- append(num_hvy ,length(df_prcp$prec_mm))
  sum_hvy <- append(sum_hvy,prcp_H_sum)
  
  print(percent_heavy)
  df_prcp <- read.csv("Downloads/dailyrain_19642023.csv")
  df_prcp["prec_mm"][df_prcp["prec_mm"] == 0] <- NA
  start_date <- as.Date("1964-01-01")  # Start date
  end_date <- as.Date("2023-12-31")    # End date
  list_of_dates <- seq(start_date, end_date, by = "day")
  df_prcp$date <- list_of_dates
  print(q)
}

df_perc_hvy <- data.frame(perc_hvy)
df_perc_hvy$Year <- seq(1964,2023,1)
df_perc_hvy$Num <- num_hvy
df_perc_hvy$Tot <- sum_hvy
df_perc_hvy$frac<- df_perc_hvy$Tot/df_perc_hvy$Num 

plot(df_perc_hvy$Year,df_perc_hvy$perc_hvy,type="b", xlab = "Year", ylab = "Rain contributed by Heavy Storms [%]",lwd =2)
abline(lm(df_perc_hvy$perc_hvy ~ df_perc_hvy$Year), lwd =2, lty =2)
title("Yearly Precip. Contributions of Heavy Rainfall")
subtitle = "Fisher Station 01/01/64 - 12/31/23"
mtext(subtitle)
summary(lm(df_perc_hvy$perc_hvy ~ df_perc_hvy$Year))

install.packages("Kendall")
library(Kendall)

Kendall(df_perc_hvy$Year,df_perc_hvy$perc_hvy)

#####SINCE 2000
df_perc_hvy_2000 <- df_perc_hvy %>% filter(df_perc_hvy$Year > 2000)
abline(lm(df_perc_hvy_2000$perc_hvy ~ df_perc_hvy_2000$Year), lwd =2, lty =2, col = "red")
summary(lm(df_perc_hvy_2000$perc_hvy ~ df_perc_hvy_2000$Year))
Kendall(df_perc_hvy_2000$Year,df_perc_hvy_2000$perc_hvy)





df_perc_hvy$Year <- as.integer(df_perc_hvy$Year)
df_perc_hvy  <- df_perc_hvy [-1,]

bp = barplot(df_perc_hvy$Num,names.arg = seq(1965,2021,1), 
        main="Number of Storms",
        xlab = "Year",
        ylab = "Number of Heavy Storms",
        ylim = c(0,15))
lines(predict(lm(df_perc_hvy$Num~df_perc_hvy$Year)),col='black')
#Create the linear regression




summary(lm(df_perc_hvy$Tot~df_perc_hvy$Year))


lm <- data.frame(seq(1965,2021,1))
colnames(lm)[1] ="x"
lm$y = (lm$x*2.4782)-4640.5789

barplot(df_perc_hvy$Tot,names.arg = seq(1965,2021,1),
        main="Heavy Storms Total Precip.",
        xlab = "Year",
        ylab = "Heavy Storms Total Precip. [mm]")







barplot(df_perc_hvy$frac,names.arg = seq(1965,2021,1),
        main="Heavy Storms Average Precip.",
        xlab = "Year",
        ylab = "Average Heavy Storm Precip. [mm]")
lines(predict(lm(df_perc_hvy$frac~df_perc_hvy$Year)),col='black')




















#######Heavy Rains Top 10% over last 30yrs####



df_prcp <- df_prcp<- read_csv("Downloads/dailyrain_19642023.csv")
df_prcp["prec_mm"][df_prcp["prec_mm"] == 0] <- NA
df_prcp$date <- as.Date(df_prcp$date, "%m/%d/%Y")
df_prcp$date <- update(df_prcp$date, year = df_prcp$year)

perc_hvy <- c()
sum_hvy <- c()
perc_90 <- c()



for(x in 1964:2023){
  x <- toString(x)
  print(x)
  df_prcp <- df_prcp[df_prcp$date >= as.Date(paste(x,"01-01", sep ="-")) & df_prcp$date <= as.Date(paste(x,"12-31",sep = "-")), ]
  q <- quantile(df_prcp$prec_mm, prob=c(.90), type=1, na.rm =TRUE)
  perc_90  <- append(perc_90 ,q)
  prcp_sum <- sum(df_prcp$prec_mm, na.rm = TRUE)
  df_prcp <- df_prcp %>% filter( prec_mm >= q)
  prcp_H_sum <- sum(df_prcp$prec_mm, na.rm = TRUE)
  sum_hvy <- append(sum_hvy,prcp_H_sum)
  percent_heavy <- prcp_H_sum/prcp_sum*100
  perc_hvy <- append(perc_hvy,percent_heavy)
  print(percent_heavy)
  df_prcp <- 
    df_prcp <- df_prcp<- read_csv("Downloads/dailyrain_19642023.csv")
  df_prcp["prec_mm"][df_prcp["prec_mm"] == 0] <- NA
  df_prcp$date <- as.Date(df_prcp$date, "%m/%d/%Y")
  df_prcp$date <- update(df_prcp$date, year = df_prcp$year)
  print(q)
}




plot(perc_90, type = "l")
x = seq(1,60,1)
lm(perc_90~x)
summary(lm(perc_90~x))

df_perc_hvy <- data.frame(perc_hvy)
df_perc_hvy$Year <- seq(1964,2023,1)
df_perc_hvy$Tot <- sum_hvy

plot(df_perc_hvy$Year,df_perc_hvy$perc_hvy,type="b", xlab = "Year", ylab = "Rain Contributed by Heavy Storms [%]", lwd =2)
abline(lm(df_perc_hvy$perc_hvy ~ df_perc_hvy$Year), lty = 2, lwd = 2)
title(expression(paste(bold("Yearly Precip. Contributions of Heavy Rainfall "))))
subtitle = "HEM 01/01/64 - 12/31/23"
mtext(subtitle)

summary(lm(df_perc_hvy$perc_hvy ~ df_perc_hvy$Year))












barplot(df_perc_hvy$Tot,names.arg = seq(1965,2022,1),
        xlab = "Year",
        ylab = "Heavy Storms Total Precip. [mm]",
        ylim = c(0,1000))
title("Heavy Storms Total Precip")
subtitle = "HEM 1964-2021"
mtext(subtitle)









##########################By YEAR###############################

df_prcp <- read.csv("hf300-05-daily-m.csv")
df_prcp$date <- as.POSIXct(df_prcp$date ,format="%Y-%m-%d", tz= "EST")
df_prcp <- df_prcp %>% separate(date, sep="-", into = c("year", "month", "day"))

dfYEARLY<- df_prcp %>% group_by(year) %>% summarise(temp.avg = mean(airt, na.rm = TRUE),
                                                    prec.sum = sum(prec, na.rm = TRUE )/10,
                                                    airmin.avg = mean(airtmin, na.rm = TRUE)
                                                       )

plot(dfYEARLY$year,dfYEARLY$prec.sum, type = "l")

plot(dfYEARLY$year,dfYEARLY$temp.avg, type = "l")

plot(dfYEARLY$year,dfYEARLY$airmin.avg, type = "l")


df_sm <- read.csv("hf206-03-HK-understory.csv")




#################################Soil Moisture by YEAR#########################

df_sm <- read.csv("hf206-03-HK-understory.csv")



#Monthly

df_sm$date <- format(as.POSIXct(df_sm$datetime), format = "%Y-%m-%d") 
df_sm <- df_sm %>% separate(date, sep="-", into = c("year", "month", "day"))

dfMONTHLY <- df_sm %>% group_by(year) %>% summarise(swc1 = mean(swc1_ave, na.rm = TRUE),
                                                    swc2 = mean(swc2_ave, na.rm = TRUE),
                                                    swc3 = mean(swc3_ave, na.rm = TRUE),
                                                    swc4 = mean(swc4_ave, na.rm = TRUE))

dfMONTHLY$swc_avg <- rowMeans(dfMONTHLY[ , c(2,3,4,5)], na.rm=TRUE)

#######################Perc heavy#######################

percent_heavy_rain <- df_perc_hvy %>% mutate(rain_binned = cut(perc_hvy, breaks= 5))
df_table <- data.frame(table(percent_heavy_rain$rain_binned))
plot(df_perc_hvy$Year,df_perc_hvy$perc_hvy, type="l")
abline(h = )





###############################COMBINED GRAPH##################################
















