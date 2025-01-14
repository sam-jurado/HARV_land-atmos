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
install.packages("EarthSystemDiagnostics")

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
library(stringr)
library(readr)


##############################LCL FUNCTION######################################

# Version 1.0 released by David Romps on September 12, 2017.
# Version 1.1 vectorized lcl.R, released on May 24, 2021.
# 
# When using this code, please cite:
# 
# @article{16lcl,
#   Title   = {Exact expression for the lifting condensation level},
#   Author  = {David M. Romps},
#   Journal = {Journal of the Atmospheric Sciences},
#   Year    = {2017},
#   Month   = dec,
#   Number  = {12},
#   Pages   = {3891--3900},
#   Volume  = {74}
# }
#
# This lcl function returns the height of the lifting condensation level
# (LCL) in meters.  The inputs are:
# - p in Pascals
# - T in Kelvins
# - Exactly one of rh, rhl, and rhs (dimensionless, from 0 to 1):
#    * The value of rh is interpreted to be the relative humidity with
#      respect to liquid water if T >= 273.15 K and with respect to ice if
#      T < 273.15 K. 
#    * The value of rhl is interpreted to be the relative humidity with
#      respect to liquid water
#    * The value of rhs is interpreted to be the relative humidity with
#      respect to ice
# - return_ldl is an optional logical flag.  If true, the lifting deposition
#   level (LDL) is returned instead of the LCL. 
# - return_min_lcl_ldl is an optional logical flag.  If true, the minimum of the
#   LCL and LDL is returned.

install.packages("LambertW")
library(LambertW)

lcl <- Vectorize(function(p,T,rh=NULL,rhl=NULL,rhs=NULL,return_ldl=FALSE,return_min_lcl_ldl=FALSE) {
  
  # Parameters
  Ttrip <- 273.16     # K
  ptrip <- 611.65     # Pa
  E0v   <- 2.3740e6   # J/kg
  E0s   <- 0.3337e6   # J/kg
  ggr   <- 9.81       # m/s^2
  rgasa <- 287.04     # J/kg/K 
  rgasv <- 461        # J/kg/K 
  cva   <- 719        # J/kg/K
  cvv   <- 1418       # J/kg/K 
  cvl   <- 4119       # J/kg/K 
  cvs   <- 1861       # J/kg/K 
  cpa   <- cva + rgasa
  cpv   <- cvv + rgasv
  
  # The saturation vapor pressure over liquid water
  pvstarl <- function(T) {
    return( ptrip * (T/Ttrip)^((cpv-cvl)/rgasv) *
              exp( (E0v - (cvv-cvl)*Ttrip) / rgasv * (1/Ttrip - 1/T) ) )
  }
  
  # The saturation vapor pressure over solid ice
  pvstars <- function(T) {
    return( ptrip * (T/Ttrip)^((cpv-cvs)/rgasv) *
              exp( (E0v + E0s - (cvv-cvs)*Ttrip) / rgasv * (1/Ttrip - 1/T) ) )
  }
  
  if (is.null(p)) { stop('Must specify p') }
  if (is.null(T)) { stop('Must specify T') }
  
  # Calculate pv from rh, rhl, or rhs
  rh_counter <- 0
  if (!is.null(rh )) { rh_counter <- rh_counter + 1 }
  if (!is.null(rhl)) { rh_counter <- rh_counter + 1 }
  if (!is.null(rhs)) { rh_counter <- rh_counter + 1 }
  if (rh_counter != 1) {
    stop('Exactly one of rh, rhl, and rhs must be specified')
  }
  if (!is.null(rh)) {
    # The variable rh is assumed to be 
    # with respect to liquid if T > Ttrip and 
    # with respect to solid if T < Ttrip
    if (T > Ttrip) {
      pv <- rh * pvstarl(T)
      rhl <- rh
      rhs <- pv / pvstars(T)
    } else {
      pv <- rh * pvstars(T)
      rhl <- pv / pvstarl(T)
      rhs <- rh
    }
  } else if (!is.null(rhl)) {
    pv <- rhl * pvstarl(T)
    rhs <- pv / pvstars(T)
    if (T > Ttrip) {
      rh <- rhl
    } else {
      rh <- rhs
    }
  } else if (!is.null(rhs)) {
    pv <- rhs * pvstars(T)
    rhl <- pv / pvstarl(T)
    if (T > Ttrip) {
      rh <- rhl
    } else {
      rh <- rhs
    }
  }
  #p = p
  #pv = pv
  #if (pv > p) {
  #return(NA)
  #}
  
  # Calculate lcl and ldl
  qv <- rgasa*pv / (rgasv*p + (rgasa-rgasv)*pv)
  rgasm <- (1-qv)*rgasa + qv*rgasv
  cpm <- (1-qv)*cpa + qv*cpv
  if (rh==0) {
    return(cpm*T/ggr)
  }
  al  <- -(cpv-cvl)/rgasv + cpm/rgasm
  bl  <- -(E0v-(cvv-cvl)*Ttrip)/(rgasv*T)
  cl  <- pv/pvstarl(T)*exp(-(E0v-(cvv-cvl)*Ttrip)/(rgasv*T))
  as  <- -(cpv-cvs)/rgasv + cpm/rgasm
  bs  <- -(E0v+E0s-(cvv-cvs)*Ttrip)/(rgasv*T)
  cs  <- pv/pvstars(T)*exp(-(E0v+E0s-(cvv-cvs)*Ttrip)/(rgasv*T))
  lcl <- cpm*T/ggr*( 1 - bl/(al*W(bl/al*cl^(1/al),-1)) )
  ldl <- cpm*T/ggr*( 1 - bs/(as*W(bs/as*cs^(1/as),-1)) )
  
  # Make zeros exact
  if (rhl==1) { lcl <- 0 }
  if (rhs==1) { ldl <- 0 }
  
  # Return either lcl or ldl
  if (return_ldl & return_min_lcl_ldl) {
    stop('return_ldl and return_min_lcl_ldl cannot both be true')
  } else if (return_ldl) {
    return(ldl)
  } else if (return_min_lcl_ldl) {
    return(min(lcl,ldl))
  } else {
    return(lcl)
  }
  
},vectorize.args=c('p','T','rh','rhl','rhs'))




#########################ANALYSIS START########################################

#Combine into one big csv
setwd("/Users/jurado/Harvard_Forest")
df_pbl1 <- read.csv("hpbl2017_cor.csv")
df_pbl2 <- read.csv("hpbl2018_cor.csv")
df_pbl3 <- read.csv("hpbl2019_cor.csv")
df_pbl4 <- read.csv("hpbl2020_cor.csv")
df_pbl5 <- read.csv("hpbl2021_cor.csv")
df_pbl6 <- read.csv("hpbl2022_cor.csv")
df_pbl7 <- read.csv("hpbl2023_cor.csv")

grand_pbl<- do.call("rbind", list(df_pbl1, df_pbl2, df_pbl3,df_pbl4,df_pbl5,df_pbl6,df_pbl7))
grand_pbl$datetime <- as.POSIXct(grand_pbl$datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
grand_pbl$datetime <- format(grand_pbl$datetime, tz="EST",usetz=TRUE)
grand_pbl$datetime <- as.POSIXct(grand_pbl$datetime, format = "%Y-%m-%d %H:%M:%S", tz = "EST")

plot(grand_pbl$datetime,grand_pbl$hpbl,type="l")


#################LCL long term PLOT##########################

#LCL input should have no Nans, RH should be 1 maximum

setwd("/Users/jurado/Harvard_Forest")
ems_hourlyustfiltnongf_2023 <- read_csv("ems_hourlyustfiltnongf_2023.csv")

#Make Timestamps
ems_hourlyustfiltnongf_2023$year <- substr(as.character(ems_hourlyustfiltnongf_2023 $TIMESTAMP_END), start = 1, stop = 4)
ems_hourlyustfiltnongf_2023$year <- as.numeric(ems_hourlyustfiltnongf_2023$year) 
ems_hourlyustfiltnongf_2023$month <- substr(as.character(ems_hourlyustfiltnongf_2023$TIMESTAMP_END), start = 5, stop = 6)
ems_hourlyustfiltnongf_2023$month <- as.numeric(ems_hourlyustfiltnongf_2023$month) 
ems_hourlyustfiltnongf_2023$day <- substr(as.character(ems_hourlyustfiltnongf_2023$TIMESTAMP_END), start = 7, stop = 8)
ems_hourlyustfiltnongf_2023$day <- as.numeric(ems_hourlyustfiltnongf_2023$day) 
ems_hourlyustfiltnongf_2023$hour <- substr(as.character(ems_hourlyustfiltnongf_2023$TIMESTAMP_END), start = 9, stop = 10)
ems_hourlyustfiltnongf_2023$hour <- as.numeric(ems_hourlyustfiltnongf_2023$hour) 

ems_hourlyustfiltnongf_2023 <- ems_hourlyustfiltnongf_2023 %>% filter(ems_hourlyustfiltnongf_2023$USTAR >0)


#growing season only
ems_hourlyustfiltnongf_2023  <- ems_hourlyustfiltnongf_2023  %>% filter(ems_hourlyustfiltnongf_2023 $month > 5)

ems_hourlyustfiltnongf_2023  <- ems_hourlyustfiltnongf_2023  %>% filter(ems_hourlyustfiltnongf_2023 $month <10)


ems_hourlyustfiltnongf_2023$RH <- replace(ems_hourlyustfiltnongf_2023$RH,ems_hourlyustfiltnongf_2023$RH>100,100)



LCL <- data.frame(P = ems_hourlyustfiltnongf_2023$PA *1000)
LCL$RH <- ems_hourlyustfiltnongf_2023$RH/100
LCL$TEMP <- ems_hourlyustfiltnongf_2023$Tair + 273.15
LCL$datetime <-ems_hourlyustfiltnongf_2023$datetime

#GET LCL VALUES
LCL <- na.omit(LCL)
LCL$LCL <- lcl(p = LCL$P,T = LCL$TEMP,rh = LCL$RH)


split <- LCL$datetime %>% str_split_fixed("-",3)
split <- data.frame(split)

#rm(EF$datetime)
LCL$year <- split$X1
LCL$month <- split$X2

LCL$year<- as.numeric(LCL$year)
LCL$month <- as.numeric(LCL$month)
LCL <- filter(LCL, month <10)
LCL <- filter(LCL, month >5)

LCL_YEARLY <- LCL %>% group_by(year) %>% summarise(lcl.mean = mean(LCL,na.rm = TRUE),
                                                   lcl.err = sd(LCL, na.rm=TRUE)/sqrt(length((LCL))))
LCL_YEARLY$year <- as.numeric(LCL_YEARLY$year)


plot(LCL_YEARLY$year,LCL_YEARLY$lcl.mean, type= "b", xlab = "Year",
     ylab = "LCL Height [m]")

abline(lm(LCL_YEARLY$lcl.mean~LCL_YEARLY$year))


summary(lm(LCL_YEARLY$lcl.mean~LCL_YEARLY$year))


##############################ABL LONG TERM PLOT################################


setwd("/Users/jurado/Harvard_Forest")


nc_data <- nc_open('hpbl.mon.mean.nc')

lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat", verbose = F)
hpbl <- ncvar_get(nc_data, "hpbl")
dim(hpbl)
fillvalue <- ncatt_get(nc_data, "hpbl", "_FillValue") #check what the fill value is
hpbl[hpbl == fillvalue$value] <- NA #replace fill values with NA

#Find coordinates closest to HF

#lat: 42.60564
#lon: -72.1105
#timeseries slice, linearly interpolated to every 30 min to be combined with lcl neon measurements

hpbl.timeseries <- hpbl[260,138,] 
plot(seq(1,538,1),hpbl.timeseries,type = "l")

df_pbl <- data.frame(hpbl = hpbl.timeseries)
df_pbl$datetime <- seq(as.POSIXct("1979-01-01 00:00:00"), as.POSIXct("2023-10-01 00:00:00"), by = "month")

split <- df_pbl$datetime %>% str_split_fixed("-",3)
split <- data.frame(split)

#rm(EF$datetime)
df_pbl$year <- split$X1
df_pbl$month <- split$X2


df_pbl_YEARLY <- df_pbl %>% group_by(year) %>% summarise(hpbl.mean = mean(hpbl,na.rm = TRUE))

#CHANGE OVER 30 YEARS
plot(df_pbl_YEARLY$year,df_pbl_YEARLY$hpbl.mean, type= "b")


#CHANGE OVER YEARS just growing season

df_pbl$year<- as.numeric(df_pbl$year)
df_pbl$month <- as.numeric(df_pbl$month)
df_pbl <- filter( df_pbl, month <10)
df_pbl <- filter( df_pbl, month >5)
df_pbl_YEARLY <- df_pbl %>% group_by(year) %>% summarise(hpbl.mean = mean(hpbl,na.rm = TRUE),
                                                         hpbl.err = sd(hpbl, na.rm=TRUE)/sqrt(length((hpbl))))

df_pbl_YEARLY$year <- as.numeric(df_pbl_YEARLY$year)

plot(df_pbl_YEARLY$year,df_pbl_YEARLY$hpbl.mean, type= "b", xlab = "Year",
     ylab = "Boundary Layer Height [m]")

abline(lm(df_pbl_YEARLY$hpbl.mean~df_pbl_YEARLY$year))


summary(lm(df_pbl_YEARLY$hpbl.mean~df_pbl_YEARLY$year))
# R = 2.823
#p-value = .0001729



##########################COMBINED ABL LCL PLOTS########################################
PBL <- data.frame(year= df_pbl_YEARLY$year)
PBL$hpbl <- df_pbl_YEARLY$hpbl.mean
PBL$hpbl.err <- df_pbl_YEARLY$hpbl.err
PBL <- merge(PBL,LCL_YEARLY,all.x=TRUE)


plot(PBL$year,PBL$hpbl, type="b", ylim = c(350,1000), xlab = "Year", ylab = "z [m]",lwd =2, col = "brown3")
lines(PBL$year,PBL$lcl.mean, type = "b", pch=2,lwd =2, col="cadetblue4")


abline(lm(df_pbl_YEARLY$hpbl.mean~df_pbl_YEARLY$year), col = "black", lty= 2)
abline(lm(LCL_YEARLY$lcl.mean~LCL_YEARLY$year), col = "black", lty= 2)


arrows(x0=PBL$year,y0 = PBL$hpbl - PBL$hpbl.err, x1=PBL$year, y1=PBL$hpbl + PBL$hpbl.err,
       code=3, angle=90, length=0.05, lwd = 1.5, col="darkgrey")
arrows(x0=PBL$year,y0 = PBL$lcl.mean - PBL$lcl.err, x1=PBL$year, y1=PBL$lcl.mean + PBL$lcl.err,
       code=3, angle=90, length=0.05, lwd = 1.5, col="darkgrey")


title("Average ABL and LCL Heights")
subtitle = "Fisher Station & NARR Dataset, June - September 1979 - 2023"
mtext(subtitle)
legend("left",legend = c("ABL","LCL"),col = c("brown3","cadetblue4"), pch = c(1,2),
       pt.lwd = c(2,2), bty = "n")



summary(lm(df_pbl_YEARLY$hpbl.mean~df_pbl_YEARLY$year))
summary(lm(LCL_YEARLY$lcl.mean~LCL_YEARLY$year))

install.packages("Kendall")
library(Kendall)

Kendall(df_pbl_YEARLY$year,df_pbl_YEARLY$hpbl.mean)
Kendall(LCL_YEARLY$year,LCL_YEARLY$lcl.mean)




##########FINAL LCL, ABL AND SOIL MOISTURE ANOMALY COMBO CODE
####Use The other data set for LE and RH in Ha1_ET_fill
###'Combine with soil moisture to make a diurnal average of LCL, ABL, and LE/H first,
###'If favorable do a crossover and rain analysis


setwd("/Users/jurado/Downloads")
Ha1_ETfill_9222 <- read_csv("Ha1_LEfill_9222.csv")

Ha1_ETfill_9222$DateTime


Ha1_ETfill_9222$year <- substr(as.character(Ha1_ETfill_9222$TIMESTAMP_END), start = 1, stop = 4)
Ha1_ETfill_9222$year <- as.numeric(Ha1_ETfill_9222$year) 
Ha1_ETfill_9222$month <- substr(as.character(Ha1_ETfill_9222$TIMESTAMP_END), start = 5, stop = 6)
Ha1_ETfill_9222$month <- as.numeric(Ha1_ETfill_9222$month) 
Ha1_ETfill_9222$day <- substr(as.character(Ha1_ETfill_9222$TIMESTAMP_END), start = 7, stop = 8)
Ha1_ETfill_9222$day <- as.numeric(Ha1_ETfill_9222$day) 
Ha1_ETfill_9222$hour <- substr(as.character(Ha1_ETfill_9222$TIMESTAMP_END), start = 9, stop = 10)
Ha1_ETfill_9222$hour <- as.numeric(Ha1_ETfill_9222$hour) 

Ha1_ETfill_9222$date <- paste(as.character(Ha1_ETfill_9222$year),as.character(Ha1_ETfill_9222$month), sep = "-0")
Ha1_ETfill_9222$date <- paste(as.character(Ha1_ETfill_9222$date),as.character(Ha1_ETfill_9222$day), sep = "-")


Ha1_ETfill_9222 <- Ha1_ETfill_9222 %>% filter(Ha1_ETfill_9222$month > 5)

Ha1_ETfill_9222 <- Ha1_ETfill_9222 %>% filter(Ha1_ETfill_9222$month <10)
Ha1_ETfill_9222 <- Ha1_ETfill_9222 %>% filter(Ha1_ETfill_9222$year >2016)

colnames(Ha1_ETfill_9222)[10] <- "datetime"


#########Using ungapfilled EMS Data#########


setwd("/Users/jurado/Harvard_Forest")
hfall_9223 <- read_csv("hfall_9223.csv")


# Convert days since 1990 to POSIXct
dates <- as.POSIXct("1990-01-01", origin = "1990-01-01", tz = "EST") + 
  hfall_9223$seq_from_1990.days. * 86400



#interpolating data into even 30 min intervals to be mergable with soil data
hfall_9223$datetime <- dates


#RH
timestamps <- hfall_9223$datetime
RH <- hfall_9223$RH.27.9m_.  # Example temperature data corresponding to the timestamps




# Interpolate data to even 30-minute intervals
interpolated_timestamps <- seq(min(timestamps), max(timestamps), by = "30 min")
interpolated_RH <- approx(timestamps, RH, xout = interpolated_timestamps, na.rm=TRUE)$y

#Temp
temp <- hfall_9223$Ta.27.9m_deg.C  # Example temperature data corresponding to the timestamps

# Interpolate data to even 30-minute intervals
interpolated_timestamps <- seq(min(timestamps), max(timestamps), by = "30 min")
interpolated_temp <- approx(timestamps, temp, xout = interpolated_timestamps,na.rm = TRUE)$y

#Pressure
P <- hfall_9223$Pamb_Pa  # Example temperature data corresponding to the timestamps

# Interpolate data to even 30-minute intervals
interpolated_P <- approx(timestamps, P, xout = interpolated_timestamps,na.rm = TRUE)$y

#Sensible Heat
H <- hfall_9223$Fheat_W.m2  # Example temperature data corresponding to the timestamps

# Interpolate data to even 30-minute intervals
interpolated_H <- approx(timestamps, H, xout = interpolated_timestamps,na.rm = TRUE)$y


#Latent Heat
LE <- hfall_9223$FH2O_e.3mol.m2.s  # Example temperature data corresponding to the timestamps

# Interpolate data to even 30-minute intervals
interpolated_LE <- approx(timestamps, LE, xout = interpolated_timestamps,na.rm = TRUE)$y

RNET <- hfall_9223$RNET_W.m2
# Interpolate data to even 30-minute intervals
interpolated_RNET <- approx(timestamps, RNET, xout = interpolated_timestamps,na.rm = TRUE)$y


hfall_9223_clean <- data.frame("datetime" =interpolated_timestamps, 
                               "RH"= interpolated_RH,
                               "temp"=interpolated_temp,
                               "P"=interpolated_P,
                               "H" =interpolated_H ,
                               "LE"= interpolated_LE,
                               "RNET" = interpolated_RNET)

data_raw = nrow(hfall_9223_clean)


#####clean only relevant ones

hfall_9223_clean <- hfall_9223_clean %>%
  filter(!is.na(RH) & !is.na(temp) & !is.na(P))



#



data_clean = nrow(hfall_9223_clean)

percent_removed = (1 - data_clean/data_raw) *100

#21% of data is NA

hfall_9223_clean$RH[hfall_9223_clean$RH > 100] <- 100



LCL <- merge(hfall_9223_clean,df_anom,by="datetime")
LCL <- merge(LCL,grand_pbl,by="datetime")

#Need hours

LCL$hour <- hour(LCL$datetime)


#still need LCL
LCL_clean <- LCL[,c("datetime","VSWCAnom.mean.med","RH","hpbl","hour","P","temp","H","LE","RNET")]



#LCL input should have no Nans, RH should be 1 maximum
LCL_clean$LCL <- lcl(p = LCL_clean$P*1000,T = LCL_clean$temp+273.15,rh = LCL_clean$RH/100)






LCLDIURNAL <- LCL_clean %>% group_by(hour) %>% summarise(LCL.avg = mean(LCL,na.rm=TRUE),
                                                         hpbl.avg = mean(hpbl,na.rm=TRUE),
                                                         swc.avg = mean(VSWCAnom.mean.med,na.rm=TRUE),
                                                         H.avg = mean(H,na.rm =TRUE),
                                                         LE.avg = mean(LE,na.rm =TRUE)
                                                         
)




quartiles <- quantile(LCL_clean$VSWCAnom.mean.med, probs=c(.25, .90), na.rm = TRUE)
mean(LCL_clean$VSWCAnom.mean.med, na.rm= TRUE)

LCL_DRY <- LCL_clean %>% filter(VSWCAnom.mean.med < -.001981223)


LCL_WET <- LCL_clean %>% filter(VSWCAnom.mean.med > -.001981223)



plot(LCL_clean$H,LCL_clean$LCL)
abline(lm(LCL_clean$LCL~LCL_clean$H))
Kendall(LCL_clean$H,LCL_clean$LCL)
summary(lm(LCL_clean$H~LCL_clean$LCL))




LCLDIURNAL_DRY <- LCL_DRY %>% group_by(hour) %>% summarise(LCL.avg = mean(LCL,na.rm=TRUE),
                                                           LCL.sd = sd(LCL,na.rm=TRUE),
                                                           hpbl.avg = mean(hpbl,na.rm=TRUE),
                                                           hpbl.sd = sd(hpbl,na.rm=TRUE),
                                                           swc.avg = mean(VSWCAnom.mean.med,na.rm=TRUE),
                                                           H.avg = mean(H,na.rm =TRUE),
                                                           LE.avg = mean(LE,na.rm =TRUE)
)

LCLDIURNAL_WET <- LCL_WET %>% group_by(hour) %>%  summarise(LCL.avg = mean(LCL,na.rm=TRUE),
                                                            LCL.sd = sd(LCL,na.rm=TRUE),
                                                            hpbl.avg = mean(hpbl,na.rm=TRUE),
                                                            hpbl.sd = sd(hpbl,na.rm=TRUE),
                                                            swc.avg = mean(VSWCAnom.mean.med,na.rm=TRUE),
                                                            H.avg = mean(H,na.rm =TRUE),
                                                            LE.avg = mean(LE,na.rm =TRUE)
)



LCLDIURNAL_DRY$LCL.err <- LCLDIURNAL_DRY$LCL.sd/sqrt(nrow(LCLDIURNAL_DRY))
LCLDIURNAL_WET$LCL.err <- LCLDIURNAL_WET$LCL.sd/sqrt(nrow(LCLDIURNAL_WET))

LCLDIURNAL_DRY$hpbl.err <- LCLDIURNAL_DRY$hpbl.sd/sqrt(nrow(LCLDIURNAL_DRY))
LCLDIURNAL_WET$hpbl.err <- LCLDIURNAL_WET$hpbl.sd/sqrt(nrow(LCLDIURNAL_WET))
#looks good! Try again with quantiles of soil moisture anomaly 



########SHADED VERSION######

#this needs error bars
plot(LCLDIURNAL_DRY$hour,LCLDIURNAL_DRY$LCL.avg,col="brown", ylim = c(0,1500),type="l",
     lwd =2, xlab = "Hour of Day", ylab = "Height [m]")
lines(LCLDIURNAL_WET$hour,LCLDIURNAL_WET$LCL.avg, col="darkgreen",lwd=2, type="l")
lines(LCLDIURNAL_DRY$hour,LCLDIURNAL_DRY$hpbl.avg, col="brown",type="l",lty = 2,lwd=2)
lines(LCLDIURNAL_WET$hour,LCLDIURNAL_WET$hpbl.avg, col="darkgreen",lty=2,lwd=2, type="l")
legend('topleft', legend = c("ABL DRY","ABL WET", "LCL DRY","LCL WET"), lty = c(2,2,1,1),
       col = c("brown","darkgreen","brown","darkgreen"), lwd=2, bty ="n")

x_polygon1 <- c(LCLDIURNAL_DRY$hour, rev(LCLDIURNAL_DRY$hour))
y_polygon1 <- c(LCLDIURNAL_DRY$LCL.avg + LCLDIURNAL_DRY$LCL.err, rev(LCLDIURNAL_DRY$LCL.avg - LCLDIURNAL_DRY$LCL.err))
polygon(x_polygon1, y_polygon1, col = alpha("brown", 0.3), border = NA)

x_polygon2 <- c(LCLDIURNAL_WET$hour, rev(LCLDIURNAL_WET$hour))
y_polygon2 <- c(LCLDIURNAL_WET$LCL.avg + LCLDIURNAL_WET$LCL.err, rev(LCLDIURNAL_WET$LCL.avg- LCLDIURNAL_WET$LCL.err))
polygon(x_polygon2, y_polygon2, col = alpha("darkgreen", 0.3), border = NA)

x_polygon3 <- c(LCLDIURNAL_DRY$hour, rev(LCLDIURNAL_DRY$hour))
y_polygon3 <- c(LCLDIURNAL_DRY$hpbl.avg + LCLDIURNAL_DRY$hpbl.err, rev(LCLDIURNAL_DRY$hpbl.avg - LCLDIURNAL_DRY$hpbl.err))
polygon(x_polygon3, y_polygon3, col = alpha("brown", 0.3), border = NA)


x_polygon4 <- c(LCLDIURNAL_WET$hour, rev(LCLDIURNAL_WET$hour))
y_polygon4 <- c(LCLDIURNAL_WET$hpbl.avg + LCLDIURNAL_WET$hpbl.err, rev(LCLDIURNAL_WET$hpbl.avg - LCLDIURNAL_WET$hpbl.err))
polygon(x_polygon4, y_polygon4, col = alpha("darkgreen", 0.3), border = NA)

title("Diurnal Average of LCL and ABL Heights")
subtitle="HARV,EMS, & NARR June-September 2017-2023"
mtext(subtitle)


# Perform independent two-sample t-test
result <- t.test(LCL_DRY$LCL, LCL_WET$LCL)

# Print the result
print(result)


# Perform independent two-sample t-test
result <- t.test(LCL_DRY$hpbl, LCL_WET$hpbl)

# Print the result
print(result)

#What is H and LE like diurnally?
#not statistically different.

#What is the temperature and RH relationship?


t.test(LCL_DRY$temp,LCL_WET$temp)

t.test(LCL_DRY$RH,LCL_WET$RH)


t.test(LCL_DRY$H,LCL_WET$H)

t.test(LCL_DRY$LE,LCL_WET$LE)




####################ATMOSPHERIC COUPLING REGIME GRAPH###########################

#filter by IQR first

quartiles_LCL <- quantile(LCL_clean$LCL, probs=c(.25, .75), na.rm = TRUE)
IQR_LCL <- (quartiles_LCL[2] - quartiles_LCL[1])*1.5
high <- IQR_LCL + quartiles_LCL[2]

LCL_clean$LCL[LCL_clean$LCL > high] <- NA


######Take the max LCL and max ABL and the soil moisture of eachand plot it per day###


LCL_clean$date <- date(LCL_clean$datetime)

LCL_clean_Daily <- LCL_clean %>% group_by(date) %>% summarise(LCL.avg = mean(LCL,na.rm=TRUE),
                                                              hpbl.avg= mean(hpbl,na.rm=TRUE),
                                                              VSWCAnom.mean.med = mean(VSWCAnom.mean.med, na.rm=TRUE),
                                                              H.avg = mean(H, na.rm=TRUE),
                                                              LE.avg_moisture = mean(LE, na.rm=TRUE),
                                                              temp= mean(temp, na.rm=TRUE),
                                                              RH = mean(RH, na.rm=TRUE),
                                                              RNET = max(RNET, na.rm=TRUE),)



LCL_binned <- LCL_clean_Daily%>% mutate(swc_binned_rnet = cut( VSWCAnom.mean.med, breaks= 30))

binned_daily <- LCL_binned  %>% group_by(swc_binned_rnet) %>% summarise(RNET = mean(RNET, na.rm=TRUE),
                                                                        swc_med = median(VSWCAnom.mean.med, na.rm=TRUE),
                                                                        freq = n())

plot(binned_daily$swc_med,binned_daily$RNET)





LCL_binned <- LCL_clean%>% mutate(swc_binned = cut( VSWCAnom.mean.med, breaks= 100))

binned <- LCL_binned  %>% group_by(swc_binned) %>% summarise(LCL.avg = mean(LCL,na.rm=TRUE),
                                                             
                                                             LCL.med = median(LCL,na.rm=TRUE),
                                                             LCL.sd = sd(LCL,na.rm=TRUE),
                                                             hpbl.avg = mean(hpbl,na.rm=TRUE),
                                                             hpbl.med = median(hpbl,na.rm =TRUE),
                                                             temp.avg = mean(temp,na.rm=TRUE),
                                                             LE.avg = mean(LE,na.rm=TRUE),
                                                             H.avg = mean(H,na.rm=TRUE),
                                                             RH.avg = mean(RH,na.rm=TRUE),
                                                             RNET.avg = max(RNET,na.rm=TRUE),
                                                             
                                                             hpbl.sd = sd(hpbl,na.rm=TRUE),
                                                             swc_avg= mean(VSWCAnom.mean.med, na.rm=TRUE),
                                                             swc_med = median(VSWCAnom.mean.med, na.rm=TRUE),
                                                             freq = n())


binned$LCL.std.err <- binned$LCL.sd/(sqrt(binned$freq))
binned$hpbl.std.err <- binned$hpbl.sd/(sqrt(binned$freq))


#filter by low 
binned <- binned %>% filter(binned$freq > 10)

#gg scatter plot



ggplot(binned) +
  geom_point(aes(x = swc_med, y = LCL.avg, size = freq),shape = 17, color = "cadetblue4") + 
  scale_size_continuous(range = c(1, 5)) +  # Adjust the range of point sizes as needed
  geom_smooth(method = "loess", fill = "cadetblue4", alpha = 0.2, col = "cadetblue4",aes(x = swc_med, y = LCL.avg,weight = freq), level = .9)+
  labs(title = "ABL and LCL Heights by Soil Moisture Anomaly",
       x = "X-axis", y = "Y-axis", size = "Count") + 
  geom_point(aes(x = swc_med, y = hpbl.avg, size = freq),shape = 16, color = "brown3") +
  geom_smooth(method = "loess", fill = "brown3", col = "brown3",alpha = 0.2, aes(x = swc_med, y = hpbl.avg, weight = freq), level = .95)+
  scale_size_continuous(range = c(1, 5)) +  # Adjust the range of point sizes as needed
  labs(title = "ABL and LCL Heights by Soil Moisture Anomaly",
       x = "Soil Moisture Anomaly", y = "24hr Maximum Height [z], Net Radiation [W/m2]", size = "N",
       subtitle = "HARV, EMS, & NARR Dataset, June - September 2017-2023")+
  geom_point(aes(x = swc_med, y = RNET.avg, size = freq),shape = 18, color = "black") +
  geom_smooth(method = "loess", fill = "black", col = "black",alpha = 0.2, aes(x = swc_med, y = RNET.avg, weight = freq), level = .95)+
  scale_size_continuous(range = c(1, 5)) +  # Adjust the range of point sizes as needed
  labs(title = "ABL and LCL Heights by Soil Moisture Anomaly",
       x = "Soil Moisture Anomaly", y = "Height [z], Binned Maximum Net Radiation [W/m2]", size = "N",
       subtitle = "HARV, EMS, & NARR Dataset, June - September 2017-2023")+
  coord_cartesian(ylim = c(0, 1000)) +
  theme(plot.title = element_text(face = "bold",hjust = 0.5),  # Center the title
        plot.subtitle = element_text(hjust = 0.5),  # Remove major grid lines
        panel.grid.minor = element_blank())   # Center the subtitle







