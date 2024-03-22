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

GetNetCDFAtCoords <- function(filename, req.coords, req.var, time.var = "time",
                              lon.var = "longitude", lat.var = "latitude"){

  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("package 'ncdf4' is needed for this function to work. Please install it.
         Linux users may have to install the 3rd party libraries libnetcdf-dev
         and libnetcdff-dev before installing ncdf4",
         call. = FALSE)
  }
  
  #browser()
  nc1 <- nc_open(filename, readunlim = FALSE)
  print(nc1)
  
  lats <- ncdf4::ncvar_get(nc1, varid = lat.var)
  lons <- ncdf4::ncvar_get(nc1, varid = lon.var)
  lons <- ifelse(lons > 180, lons - 360, lons)
  
  nc.coords <- expand.grid(lons = lons, lats = lats)
  
  ind <- sapply(1:nrow(req.coords), function(x) {
    which.min(geosphere::distHaversine(req.coords[x, ], nc.coords))
  })
  
  print(ind)
  
  req.nc.coords <- nc.coords[ind,]
  
  lon.inds.to.get <- sapply(req.nc.coords[,1], function(x) which(lons == x))
  lat.inds.to.get <- sapply(req.nc.coords[,2], function(x) which(lats == x))
  
  coord.inds.to.get <- cbind(lon.inds.to.get, lat.inds.to.get)
  
  system.time(
    dat.out <- coord.inds.to.get %>% 
      tibble::as_tibble() %>% 
      dplyr::mutate(row = 1:n()) %>% 
      dplyr::rowwise() %>% 
      dplyr::do({
        tmp <- ncdf4::ncvar_get(
          nc1,
          varid = req.var,
          start = c(.$lon.inds.to.get, .$lat.inds.to.get, 1),
          count = c(1, 1, -1)
        )
        data.frame(location = as.character(.$row), date = 1:length(tmp),
                   req.var = tmp)
      })
  )
  
  #dat.out <- tibble::as_tibble(t(dat.out))
  
  dates <- ncdf4::ncvar_get(nc1, varid = time.var, start = c(1), count = c(-1))
  
  ncdf4::nc_close(nc1)
  
  req.coords$location <- as.character(1:nrow(req.coords))
  
  #dat.out <- tidyr::gather(dat.out, location, UQ(req.var), -date)
  
  dat.out <- dplyr::left_join(req.coords, dat.out)
  
  # dat.out <- dat.out %>% 
  #   group_by(location) %>% 
  #   pivot_wider(names_from = location, values_from = req.var)
  
  dat.out$date <- rep(dates, nrow(req.coords))
  
  
  return(dat.out)
}






setwd("/Users/jurado/Harvard_Forest")


nc_data <- nc_open('hpbl.2022.nc')

hpbl <- ncvar_get(nc_data, "hpbl")
dim(hpbl)
fillvalue <- ncatt_get(nc_data, "hpbl", "_FillValue") #check what the fill value is
hpbl[hpbl == fillvalue$value] <- NA #replace fill values with NA

#Find coordinates closest to HF



#lat: 42.60564
#lon: -72.1105
#timeseries slice, linearly interpolated to every 30 min to be combined with lcl neon measurements


#maybe should have been [260,138]

#actually used
#lon -72.11015
#lat 42.60564

hpbl.timeseries <- hpbl[260,138,] 
plot(seq(1,2920,1),hpbl.timeseries,type = "l")

df_pbl <- data.frame(hpbl = hpbl.timeseries)
df_pbl$datetime <- seq(as.POSIXct("2022-01-01 00:00:00"), as.POSIXct("2022-12-31 21:00"), by = "3 hours")


#Interpolate into half hourly data#

newdt <- seq.POSIXt(as.POSIXct("2022-01-01 00:00:00"), as.POSIXct("2022-12-31 21:00"), by='30 min')
df_pbl <- merge(df_pbl, data.frame(datetime=newdt), by.x='datetime', all.y=TRUE)

library(zoo)
df_pbl <- df_pbl %>%
  mutate(hpbl = na.approx(hpbl))

#Give time series date, write as a csv do again with other years and splice together
plot(df_pbl$datetime,df_pbl$hpbl, type="l")

write.csv(df_pbl, "/Users/jurado/Harvard_Forest/hpbl2022_cor.csv", row.names=FALSE)






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


#use EMS VPD data for LCL since it already has RH. Use NEON (df_new from SL_T.R) soil moisture for consistency

library(readr)
setwd("/Users/jurado/")
Ha1_ETfill_9222 <- read_csv("Downloads/Ha1_ETfill_9222.csv", na = "-9999")



#making time stamps
Ha1_ETfill_9222$year <- substr(as.character(Ha1_ETfill_9222$TIMESTAMP_END), start = 1, stop = 4)
Ha1_ETfill_9222$year <- as.numeric(Ha1_ETfill_9222$year) 
Ha1_ETfill_9222$month <- substr(as.character(Ha1_ETfill_9222$TIMESTAMP_END), start = 5, stop = 6)
Ha1_ETfill_9222$month <- as.numeric(Ha1_ETfill_9222$month) 
Ha1_ETfill_9222$day <- substr(as.character(Ha1_ETfill_9222$TIMESTAMP_END), start = 7, stop = 8)
Ha1_ETfill_9222$day <- as.numeric(Ha1_ETfill_9222$day) 
Ha1_ETfill_9222$hour <- substr(as.character(Ha1_ETfill_9222$TIMESTAMP_END), start = 9, stop = 10)
Ha1_ETfill_9222$hour <- as.numeric(Ha1_ETfill_9222$hour) 

Ha1_ETfill_9222$date <- paste(as.character(Ha1_ETfill_9222$year),as.character(Ha1_ETfill_9222$month), sep = "-0")

Ha1_ETfill_9222 <- Ha1_ETfill_9222 %>% filter(Ha1_ETfill_9222$month > 5)

Ha1_ETfill_9222 <- Ha1_ETfill_9222 %>% filter(Ha1_ETfill_9222$month <10)






setwd("/Users/jurado/Harvard_Forest")

df_prcp <- read.csv("hf001-06-daily-m.csv", na = "-9999", header=TRUE) #EMS data







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

##############################CROSSOVER CALCULATION#############################


folder_list = c("NEON.D01.HARV.DP1.00098.001.2017-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2017-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2017-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2017-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2018-06.basic.20230127T120753Z.RELEASE-2023 2",
                "NEON.D01.HARV.DP1.00098.001.2018-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2018-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2018-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2019-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2019-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2019-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2019-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2020-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2020-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2020-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2020-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2021-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2021-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2021-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2021-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2022-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00098.001.2022-07.basic.20221204T200354Z.PROVISIONAL",
                "NEON.D01.HARV.DP1.00098.001.2022-08.basic.20221204T234919Z.PROVISIONAL",
                "NEON.D01.HARV.DP1.00098.001.2022-09.basic.20221002T184526Z.PROVISIONAL",
                "NEON.D01.HARV.DP1.00098.001.2023-06.basic.20230706T011608Z.PROVISIONAL",
                "NEON.D01.HARV.DP1.00098.001.2023-07.basic.20230809T214353Z.PROVISIONAL",
                "NEON.D01.HARV.DP1.00098.001.2023-08.basic.20230905T021136Z.PROVISIONAL",
                "NEON.D01.HARV.DP1.00098.001.2023-09.basic.20231002T184625Z.PROVISIONAL 7")

file_list = c("NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2017-06.basic.20221204T182005Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2017-07.basic.20221204T205553Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2017-08.basic.20221204T231105Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2017-09.basic.20221204T185602Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2018-06.basic.20221204T210933Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2018-07.basic.20221204T200733Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2018-08.basic.20221204T213954Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2018-09.basic.20221204T223820Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2019-06.basic.20221204T221150Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2019-07.basic.20221204T191754Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2019-08.basic.20221204T204302Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2019-09.basic.20221204T202526Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2020-06.basic.20221204T182311Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2020-07.basic.20221204T201556Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2020-08.basic.20221204T181906Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2020-09.basic.20221204T195405Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2021-06.basic.20221204T190324Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2021-07.basic.20221204T183618Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2021-08.basic.20221204T193409Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2021-09.basic.20221204T184259Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2022-06.basic.20221204T194842Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2022-07.basic.20221204T200245Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2022-08.basic.20221204T234757Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2022-09.basic.20221002T184427Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2023-06.basic.20230706T011608Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2023-07.basic.20230809T214353Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2023-08.basic.20230905T021136Z.csv",
              "NEON.D01.HARV.DP1.00098.001.000.060.030.RH_30min.2023-09.basic.20231002T184625Z.csv")


df_RH <- data.frame()


#DATA SCRAPING LOOP#
for (x in 1:length(folder_list)){
  wd = paste("/Users/jurado/Harvard_Forest/filesToStack00098/",folder_list[x],sep ="") 
  setwd(wd)
  file = file_list[x]
  RH <- read.csv(file)
  df_RH <- rbind(df_RH,RH)
  
}
df_RH$datetime<- as.POSIXct(df_RH$endDateTime, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC") 

df_RH$datetime<- as.POSIXct(df_RH$datetime, format = "%Y-%m-%dT%H:%M:%S", tz = "EST") 

#P and temp and RH


LCL <- merge(df_anom,df_RH, by = "datetime")

LCL <- data.frame(P.mean = LCL$P.mean, tempRHMean = LCL$tempRHMean, RHMean = LCL$RHMean, datetime = LCL$datetime )



LCL <- na.omit(LCL) #check this to ensure its just the 

#LCL input should have no Nans, RH should be 1 maximum
LCL$LCL <- lcl(p = LCL$P.mean*1000,T = LCL$tempRHMean+273.15,rh = LCL$RHMean/100)



LCL <- merge(LCL,grand_pbl, by = "datetime")
LCL <- merge(LCL,df_new, by = "datetime")
LCL <- merge(LCL,df, by = "datetime")

split <- LCL$datetime%>% str_split_fixed(" ",2)
split <- data.frame(split)
LCL$date <- split$X1
LCL$hour <- split$X2



LCL$hour  <- as.POSIXct(LCL$hour , format = "%H:%M:%S")
LCL <- LCL %>% filter(hour >= as.POSIXct('2023-12-22 5:30:00') & hour <= as.POSIXct('2023-12-22 17:30:00'))

LCL$LCL <- ifelse(LCL$LCL<0, 0,LCL$LCL )




#####THIS IS FOR PARAMETERS FOR CLASS
LCLDAILY_DRY <- LCL %>% filter(VSWC.mean < .21 )
LCLDAILY_BUFFER <- LCL %>% filter(VSWC.mean  > .095 & VSWC.mean  < .21 )
LCLDAILY_TRANS <- LCL %>% filter(VSWC.mean  > .21 )


#This is daytime only
LCLDAILY <- LCL %>% group_by(date) %>% summarise(LCL.avg = mean(LCL, na.rm = TRUE),
                                                 LCL.sd = sd(LCL,na.rm=TRUE),
                                                 LCL.max = max(LCL,na.rm=TRUE),
                                               ABL.avg = mean(hpbl, na.rm = TRUE),
                                               ABL.max = max(hpbl,na.rm=TRUE),
                                               ABL.sd = sd(hpbl,na.rm=TRUE),
                                               swc.avg = mean(VSWC.mean, na.rm = TRUE),
                                               temp.avg = mean(tempRHMean, na.rm = TRUE),
                                               RH.avg = mean(RHMean, na.rm = TRUE),
                                               LE.avg = mean(LE, na.rm = TRUE),
                                               H.avg = mean(H, na.rm = TRUE),
                                               Rnet.avg = mean(Rnet, na.rm = TRUE),
                                               P.mean = mean(P.mean.x, na.rm = TRUE))




plot(LCL$datetime,LCL$hpbl)
points(LCL$datetime,LCL$LCL, col = "red")




#ADD in rain data, H data, LE data

df_prcp <- data.frame(df_prcp$date,df_prcp$prec)
colnames(df_prcp)[1] ="date"
colnames(df_prcp)[2] ="prec"


LCLDAILY <- merge(LCLDAILY,df_prcp, by ="date")
LCLDAILY$CROSSOVER <- NA



#crossover must be reimagined
####CROSSOVER Check#####

LCLDAILY$CROSSOVER <- ifelse(LCLDAILY$LCL.max<LCLDAILY$ABL.max, TRUE,FALSE)
LCLDAILY$CROSSOVER_PREC <- ifelse(LCLDAILY$CROSSOVER & LCLDAILY$prec > 0, TRUE,FALSE)

LCLDAILY_DRY <- LCLDAILY %>% filter(swc.avg < .21 )
LCLDAILY_BUFFER <- LCLDAILY %>% filter(swc.avg > .095 & swc.avg < .21 )
LCLDAILY_TRANS <- LCLDAILY %>% filter(swc.avg > .21 )


DRY_TEMP <- mean(LCLDAILY_DRY$temp.avg, na.rm = TRUE)
DRY_RH <- mean(LCLDAILY_DRY$RH.avg, na.rm = TRUE)

BUFFER_TEMP <- mean(LCLDAILY_BUFFER$temp.avg, na.rm = TRUE)
BUFFER_RH <- mean(LCLDAILY_BUFFER$RH.avg, na.rm = TRUE)


TRANS_TEMP <- mean(LCLDAILY_TRANS$temp.avg, na.rm = TRUE)
TRANS_RH <- mean(LCLDAILY_TRANS$RH.avg, na.rm = TRUE)





CROSSOVER_DRY_perc <- sum(LCLDAILY_DRY$CROSSOVER_PREC, na.rm=TRUE)/143
CROSSOVER_BUFFER_perc <- sum(LCLDAILY_BUFFER$CROSSOVER_PREC, na.rm=TRUE)/543
CROSSOVER_TRANS_perc <- sum(LCLDAILY_TRANS$CROSSOVER_PREC, na.rm=TRUE)/61

CROSSOVER_DRY_cross <- sum(LCLDAILY_DRY$CROSSOVER, na.rm=TRUE)/143
CROSSOVER_BUFFER_cross <- sum(LCLDAILY_BUFFER$CROSSOVER, na.rm=TRUE)/543
CROSSOVER_TRANS_cross <- sum(LCLDAILY_TRANS$CROSSOVER, na.rm=TRUE)/61






IQR1 <- IQR(LCLDAILY_DRY_prec$prec)*3
IQR2 <- IQR(LCLDAILY_BUFFER_prec$prec)*3
IQR3 <- IQR(LCLDAILY_TRANS_prec$prec)*3
quartiles <- quantile(LCLDAILY_TRANS_prec$prec, probs=c(.25, .75), na.rm = TRUE)





################################BINNED APPROACH################################
LCLDAILY <- LCLDAILY %>% mutate(swc_binned = cut(swc.avg, breaks=3))

binned <- LCLDAILY %>% group_by(swc_binned) %>% summarise(CROSS_PREC_SUM = sum(CROSSOVER_PREC, na.rm=TRUE),
                                                          CROSS_SUM = sum(CROSSOVER, na.rm=TRUE),
                                                          n = n())


binned$ratio_precip <- binned$CROSS_PREC_SUM/binned$n
binned$ratio_cross <- binned$CROSS_SUM/binned$n


#Binned by Regime, filtered for extreme outliers 3IQR)
LCLDAILY_DRY_prec <- LCLDAILY %>% filter(swc.avg < .21 & CROSSOVER_PREC == TRUE & prec < 45.3) 
LCLDAILY_BUFFER_prec<- LCLDAILY %>% filter(swc.avg > .095 & swc.avg < .21 & CROSSOVER_PREC == TRUE & prec < 36.6 )
LCLDAILY_TRANS_prec <- LCLDAILY %>% filter(swc.avg > .21 & CROSSOVER_PREC == TRUE & prec < 101.7)



############################CROSSOVER BAR PLOT##################################



Dry=c(CROSSOVER_DRY_perc*100,CROSSOVER_DRY_cross*100)
Buffer=c(CROSSOVER_BUFFER_perc*100,CROSSOVER_BUFFER_cross*100)
Transition=c(CROSSOVER_TRANS_perc*100,CROSSOVER_TRANS_cross*100)


data <- data.frame(Dry,Buffer,Transition)



# plotting multiple bar plots
barplot(as.matrix(data),
        
        # setting y label only
        # because x-label will be our
        # barplots name
        ylab="Frequency of Days [%]",
        xlab = "Evapotransipiration Regime",
        ylim = c(0,100),
        
        # to plot the bars vertically
        beside=TRUE,
        
        names.arg=c("DRY","BUFFER","TRANSITION")
)



#mean rain

Dry=mean(LCLDAILY_DRY$prec)
Buffer=mean(LCLDAILY_BUFFER$prec)
Transition=mean(LCLDAILY_TRANS$prec)


data <- data.frame(Dry,Buffer,Transition)


par(new = TRUE)  
barplot(as.matrix(data),
        
        # setting y label only
        # because x-label will be our
        # barplots name
        ylab="",
        xlab = "",
        ylim = c(0,100),
        
        # to plot the bars vertically
        beside=TRUE,
        
        names.arg=c("","",""),
        col = "white"
)

par(new = TRUE)  
barplot(as.matrix(data),
        
        # setting y label only
        # because x-label will be our
        # barplots name
        ylab="",
        xlab = "",
        ylim = c(0,100),
        
        # to plot the bars vertically
        beside=TRUE,
        
        names.arg=c("","",""),
        density = c(5,5,5)
        
)




axis(side = 4, at = pretty(seq(0, 20, 5))) 

legend("topleft",
       legend = c("Crossover Events","Crossover & Precip."," Mean Precipitation"),
       pch = 15,
       col = c("lightgrey","black","white"),
       bty = "n",
       cex=.9)

title("Crossover Events by ET Regime")
subtitle = "HARV June-September 2017-2023"
mtext(subtitle)
#need to add error bars for rain



######DAILY PLOT######

#isolate hours for diurnal average 
#is ABL aligned?
LCL$hour <- substr(as.character(LCL$datetime), start = 12, stop = 13)

LCLDAILY_DRY <- LCL %>% filter(VSWC.mean < .095 )
LCLDAILY_DRY <- LCLDAILY_DRY %>% group_by(hour) %>% summarise(H.avg = mean(H,na.rm=TRUE),
                                                     LE.avg = mean(LE,na.rm=TRUE),
                                                     ABL.avg = mean(hpbl,na.rm=TRUE),
                                                     LCL.avg = mean(LCL,na.rm=TRUE),
                                                     G.avg = mean(G,na.rm=TRUE),
                                                     Rnet.avg = mean(Rnet,na.rm=TRUE))

LCLDAILY_BUFFER <- LCL %>% filter(VSWC.mean > .095 & VSWC.mean < .21 )
LCLDAILY_BUFFER <- LCLDAILY_BUFFER %>% group_by(hour) %>% summarise(H.avg = mean(H,na.rm=TRUE),
                                                              LE.avg = mean(LE,na.rm=TRUE),
                                                              ABL.avg = mean(hpbl,na.rm=TRUE),
                                                              LCL.avg = mean(LCL,na.rm=TRUE),
                                                              G.avg = mean(G,na.rm=TRUE),
                                                              Rnet.avg = mean(Rnet,na.rm=TRUE))

LCLDAILY_TRANS <- LCL %>% filter(VSWC.mean > .21 )
LCLDAILY_TRANS <- LCLDAILY_TRANS%>% group_by(hour) %>% summarise(H.avg = mean(H,na.rm=TRUE),
                                                              LE.avg = mean(LE,na.rm=TRUE),
                                                              ABL.avg = mean(hpbl,na.rm=TRUE),
                                                              LCL.avg = mean(LCL,na.rm=TRUE),
                                                              G.avg = mean(G,na.rm=TRUE),
                                                              Rnet.avg = mean(Rnet,na.rm=TRUE))



###ADD ERROR BARS###
plot(LCLDAILY_DRY$hour,LCLDAILY_DRY$ABL.avg, ylim =c(300,1500), type ="l", lwd = 2, col = "orange4",
     xlab = "Hour", ylab = "Height [m]")
lines(LCLDAILY_DRY$hour,LCLDAILY_DRY$LCL.avg, lwd =2, col = "orange4", lty =2)
lines(LCLDAILY_DRY$hour,LCLDAILY_TRANS$ABL.avg, lty = 1, lwd = 2, col = "forestgreen")
lines(LCLDAILY_DRY$hour,LCLDAILY_TRANS$LCL.avg, lty = 2, lwd = 2, col = "forestgreen")
lines(LCLDAILY_DRY$hour,LCLDAILY_BUFFER$ABL.avg, lty = 1, lwd = 2, col = "orange3")
lines(LCLDAILY_DRY$hour,LCLDAILY_BUFFER$LCL.avg, lty = 2, lwd = 2, col = "orange3")
legend("topleft", legend = c("ABL Dry","ABL Trans.", "ABL Wet","LCL Dry","LCL Trans.", "LCL Wet"),
       lwd = 2, lty = c(1,1,1,2,2,2), col = c("orange4","orange3","forestgreen","orange4","orange3","forestgreen"),
       bty = "n")
title("Mean Daytime ABL and LCL Heights by Soil Moisture Regime")
subtitle = "HARV & EMS June-September 2017-2023"
mtext(subtitle)






plot(LCLDAILY_DRY$hour,LCLDAILY_DRY$LE.avg, type ="l", ylim = c(-50,250), col = "orange4", lwd=2,
     xlab ="Hour", ylab = "Flux [W/m^2]")
lines(LCLDAILY_DRY$hour,LCLDAILY_DRY$H.avg, lty = 2,col = "orange4", lwd=2)
lines(LCLDAILY_DRY$hour,LCLDAILY_TRANS$LE.avg, lty = 1,col = "forestgreen", lwd=2)
lines(LCLDAILY_DRY$hour,LCLDAILY_TRANS$H.avg, lty = 2,col = "forestgreen", lwd=2)
lines(LCLDAILY_DRY$hour,LCLDAILY_BUFFER$LE.avg, lty = 1,col = "orange3", lwd=2)
lines(LCLDAILY_DRY$hour,LCLDAILY_BUFFER$H.avg, lty = 2,col = "orange3", lwd=2)
legend("topleft", legend = c("LE Dry","LE Trans.", "LE Wet","H Dry","H Trans.", "H Wet"),
       lwd = 2, lty = c(1,1,1,2,2,2), col = c("orange4","orange3","forestgreen","orange4","orange3","forestgreen"),
       bty = "n")
title("Mean Daytime LE and H by Soil Moisture Regime")
subtitle = "HARV June-September 2017-2023"
mtext(subtitle)


TRANS_LCL <- lcl( p = 97225, 17.27+273.15, .7519)
DRY_LCL <- lcl( p = 97225, 20.36+273.15, .7009)

#TRANS LCL = 562.1m
#DRY LCL = 716.1m
#diff = 154 m



#temperature 







# Create data 
gfg <- data.frame(data = c(69.12,20.11,70.09,20.37,75.20,17.27),
                  grp = c("BUFFER","BUFFER","DRY", "DRY","TRANS","TRANS"),
                  subgroup = c("RH","Temperature","RH","Temperature","RH","Temperature"))
# Modifying data 
gfg <- reshape(gfg,idvar = "subgroup", 
               timevar = "grp", 
               direction = "wide") 

row.names(gfg) <- gfg$subgroup 
gfg <- gfg[ , 2:ncol(gfg)] 
colnames(gfg) <- c("group 1", "group 2", 
                   "group 3","group 4") 
gfg <- as.matrix(gfg) 

# Create grouped barplot 
barplot(height = gfg,beside = TRUE)













 #################LCL long term PLOT##########################

#LCL input should have no Nans, RH should be 1 maximum
library(readr)
hf004 <- read_csv("hf004-01-final.csv")

split <- hf004$datetime %>% str_split_fixed("-",3)
split <- data.frame(split)

#rm(EF$datetime)
hf004$year <- split$X1
hf004$month <- split$X2
hf004$year<- as.numeric(hf004$year)
hf004$month <- as.numeric(hf004$month)
hf004 <- filter( hf004, month <10)
hf004 <- filter( hf004, month >5)
hf004$rh_27_9m <- replace(hf004$rh_27_9m,hf004$rh_27_9m>100,100)



LCL <- data.frame(P = hf004$p_amb)
LCL$RH <- hf004$rh_27_9m/100
LCL$TEMP <- hf004$ta_27_9m+273.15
LCL$datetime <- hf004$datetime

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


#FILTER TO GROWING SEASON ONLY
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



#####ABL through Mechanical shear only

library(readxl)
setwd("/Users/jurado/")
AMF_US_MMS_FLUXNET <- read_excel("Downloads/AMF_US-MMS_FLUXNET.xlsx")

#Make Timestamps
AMF_US_MMS_FLUXNET$year <- substr(as.character(AMF_US_MMS_FLUXNET$TIMESTAMP_END), start = 1, stop = 4)
AMF_US_MMS_FLUXNET$year <- as.numeric(AMF_US_MMS_FLUXNET$year) 
AMF_US_MMS_FLUXNET$month <- substr(as.character(AMF_US_MMS_FLUXNET$TIMESTAMP_END), start = 5, stop = 6)
AMF_US_MMS_FLUXNET$month <- as.numeric(AMF_US_MMS_FLUXNET$month) 
AMF_US_MMS_FLUXNET$day <- substr(as.character(AMF_US_MMS_FLUXNET$TIMESTAMP_END), start = 7, stop = 8)
AMF_US_MMS_FLUXNET$day <- as.numeric(AMF_US_MMS_FLUXNET$day) 
AMF_US_MMS_FLUXNET$hour <- substr(as.character(AMF_US_MMS_FLUXNET$TIMESTAMP_END), start = 9, stop = 10)
AMF_US_MMS_FLUXNET$hour <- as.numeric(AMF_US_MMS_FLUXNET$hour) 

AMF_US_MMS_FLUXNET <- AMF_US_MMS_FLUXNET %>% filter(AMF_US_MMS_FLUXNET$USTAR >0)

AMF_US_MMS_FLUXNET$hpbl <- .3*(AMF_US_MMS_FLUXNET$USTAR)/(2*(7.27*10**-5)*sin(42*(pi/180)))

#growing season only
AMF_US_MMS_FLUXNET <- AMF_US_MMS_FLUXNET %>% filter(AMF_US_MMS_FLUXNET$month > 5)

AMF_US_MMS_FLUXNET <- AMF_US_MMS_FLUXNET %>% filter(AMF_US_MMS_FLUXNET$month <10)




ABL_YEARLY <- AMF_US_MMS_FLUXNET %>% group_by(year) %>% summarise(mean.hpbl = mean(hpbl, na.rm =TRUE),
                                                                  sd.hpbl = sd(hpbl, na.rm =TRUE),
                                                                  n = n())
ABL_YEARLY$std_err <- ABL_YEARLY$sd.hpbl/sqrt(ABL_YEARLY$n)


plot(ABL_YEARLY$year,ABL_YEARLY$mean.hpbl, ylim =c(200,1100))





############COMBINED PLOTS###############
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

plot(x, y, ylim=c(-3, 3), xlab="x", ylab="y", pch=16, cex=2)
# Add error bars
arrows(x0=x, y0=y-y.sd, x1=x, y1=y+y.sd, code=3, angle=90, length=0.1)







#Averaging into bins, but this time daily 
EF <- LCLDAILY %>% mutate(swc_binned = cut(swc.avg, breaks=30))

binned <- EF %>% group_by(swc_binned) %>% summarise(LCL.avg = mean(LCL.avg, na.rm=TRUE),
                                                    ABL.avg = mean(ABL.avg, na.rm = TRUE),
                                                    swc.avg = mean(swc.avg, na.rm=TRUE),
                                                    H.avg = mean(H.avg, na.rm=TRUE),
                                                    LE.avg = mean(LE.avg,na.rm = TRUE),
                                                    temp.avg = mean(temp.avg, na.rm = TRUE),
                                                    RH.avg = mean(RH.avg, na.rm = TRUE))

df_table <- data.frame(table(EF$swc_binned))
colnames(df_table)[1] ="swc_binned"
binned <- merge(binned,df_table,by = "swc_binned")

binned <- binned %>% filter(Freq > 3)








#SL_T REGIME

plot(binned$swc.avg,binned$LCL.avg, pch = 2,ylab = "Height [m]", xlab = "Soil Water Content",
     type ="p",ylim =c(400,1200), col = "cadetblue4", lwd = 2)
points(binned$swc.avg,binned$ABL.avg, pch = 1, col ="brown3" , lwd = 2)

y = weightedLowess(binned$swc.avg, binned$LCL.avg, weights = binned$Freq,
                   delta=NULL, npts = 200, span = 0.06, iterations = 3)

lines(binned$swc.avg,y$fitted, col="cadetblue4",lty = 2,lwd =2)

y = weightedLowess(binned$swc.avg, binned$ABL.avg, weights = binned$Freq,
                   delta=NULL, npts = 200, span = 0.06, iterations = 3)

lines(binned$swc.avg,y$fitted, col="brown3",lty = 2,lwd =2)
title("ABL and LCL Heights")
subtitle = "HARV & NARR Dataset, June - September 2016-2022"
mtext(subtitle)

plot(binned$swc.avg,binned$LCL.avg-binned$ABL.avg,)

plot(LCLDAILY$swc.avg,LCLDAILY$LCL.avg - LCLDAILY$ABL.avg)


#what are the days like when the LCL is higher than the ABL?

DELTA_LCLDAILY<- LCLDAILY %>% filter(LCLDAILY$LCL.avg - LCLDAILY$ABL.avg  > 0)


##########FINAL LCL, ABL AND SOIL MOISTURE ?ANOMALY COMBO CODE
####Use The other data set for LE and RH in Ha1_ET_fill
###'Combine with soil moisture to make a diurnal average of LCL, ABL, and LE/H first,
###'If favorable do a crossover and rain analysis

library(readr)


setwd("/Users/jurado/Downloads/")
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








#not the best merge, taking the actual values rather than average.
LCL <- merge(Ha1_ETfill_9222,df_anom,by="datetime")
LCL <- merge(LCL,grand_pbl,by="datetime")

#still need LCL
LCL_clean <- LCL[,c("datetime","VSWCAnom.mean.med","LE_f","RH_PI_F_1_2_1","hpbl","hour","PA_PI_F","TA_PI_F_1_2_1","H")]


LCL_clean <- LCL_clean[!is.na(LCL_clean$RH_PI_F_1_2_1), ]
LCL_clean <- LCL_clean[!is.na(LCL_clean$TA_PI_F_1_2_1), ]
LCL_clean <- LCL_clean[!is.na(LCL_clean$PA_PI_F), ]



#LCL input should have no Nans, RH should be 1 maximum
LCL_clean$LCL <- lcl(p = LCL_clean$PA_PI_F*1000,T = LCL_clean$TA_PI_F_1_2_1+273.15,rh = LCL_clean$RH_PI_F_1_2_1/100)


LCLDIURNAL <- LCL_clean %>% group_by(hour) %>% summarise(LCL.avg = mean(LCL,na.rm=TRUE),
                                                         hpbl.avg = mean(hpbl,na.rm=TRUE),
                                                         swc.avg = mean(VSWCAnom.mean.med,na.rm=TRUE),
                                                         LE.avg = mean(LE_f, na.rm =TRUE),
                                                         H.avg = mean(H, na.rm = TRUE))



plot(LCLDIURNAL$hour,LCLDIURNAL$H.avg)
lines(LCLDIURNAL$hour,LCLDIURNAL$LE.avg)


quartiles <- quantile(LCL_clean$VSWCAnom.mean.med, probs=c(.10, .90), na.rm = TRUE)

LCL_DRY <- LCL_clean %>% filter(VSWCAnom.mean.med < -.045)


LCL_WET <- LCL_clean %>% filter(VSWCAnom.mean.med > 0.0338)

LCLDIURNAL_DRY <- LCL_DRY %>% group_by(hour) %>% summarise(LCL.avg = mean(LCL,na.rm=TRUE),
                                                           LCL.sd = sd(LCL,na.rm=TRUE),
                                                           hpbl.avg = mean(hpbl,na.rm=TRUE),
                                                           hpbl.sd = sd(hpbl,na.rm=TRUE),
                                                           swc.avg = mean(VSWCAnom.mean.med,na.rm=TRUE),
                                                           LE.avg = mean(LE_f, na.rm =TRUE),
                                                           H.avg = mean(H, na.rm = TRUE))

LCLDIURNAL_WET <- LCL_WET %>% group_by(hour) %>%  summarise(LCL.avg = mean(LCL,na.rm=TRUE),
                                                            LCL.sd = sd(LCL,na.rm=TRUE),
                                                            hpbl.avg = mean(hpbl,na.rm=TRUE),
                                                            hpbl.sd = sd(hpbl,na.rm=TRUE),
                                                            swc.avg = mean(VSWCAnom.mean.med,na.rm=TRUE),
                                                            LE.avg = mean(LE_f, na.rm =TRUE),
                                                            H.avg = mean(H, na.rm = TRUE))


LCLDIURNAL_DRY$LCL.err <- LCLDIURNAL_DRY$LCL.sd/sqrt(nrow(LCLDIURNAL_DRY))
LCLDIURNAL_WET$LCL.err <- LCLDIURNAL_WET$LCL.sd/sqrt(nrow(LCLDIURNAL_WET))

LCLDIURNAL_DRY$hpbl.err <- LCLDIURNAL_DRY$hpbl.sd/sqrt(nrow(LCLDIURNAL_DRY))
LCLDIURNAL_WET$hpbl.err <- LCLDIURNAL_WET$hpbl.sd/sqrt(nrow(LCLDIURNAL_WET))
#looks good! Try again with quantiles of soil moisture anomaly 



#this needs error bars
plot(LCLDIURNAL_DRY$hour,LCLDIURNAL_DRY$LCL.avg,col="brown", ylim = c(0,1500),type="l",
     lwd =2, xlab = "Hour of Day", ylab = "Height [m]")
lines(LCLDIURNAL_WET$hour,LCLDIURNAL_WET$LCL.avg, col="darkgreen",lwd=2, type="l")
lines(LCLDIURNAL_DRY$hour,LCLDIURNAL_DRY$hpbl.avg, col="brown",type="l",lty = 2,lwd=2)
lines(LCLDIURNAL_WET$hour,LCLDIURNAL_WET$hpbl.avg, col="darkgreen",lty=2,lwd=2, type="l")
legend('topleft', legend = c("ABL DRY","ABL WET", "LCL DRY","LCL WET"), lty = c(2,2,1,1),
       col = c("brown","darkgreen","brown","darkgreen"), lwd=2, bty ="n")

arrows(LCLDIURNAL_DRY$hour, y0 = LCLDIURNAL_DRY$LCL.avg -  LCLDIURNAL_DRY$LCL.err , LCLDIURNAL_DRY$hour, y1= LCLDIURNAL_DRY$LCL.avg   +  LCLDIURNAL_DRY$LCL.err ,
       code=3, angle=90, length=0.05, lwd = 1.5, col = alpha("brown",.65))
arrows(LCLDIURNAL_WET$hour, y0 = LCLDIURNAL_WET$LCL.avg -  LCLDIURNAL_WET$LCL.err , LCLDIURNAL_WET$hour, y1= LCLDIURNAL_WET$LCL.avg   +  LCLDIURNAL_WET$LCL.err ,
       code=3, angle=90, length=0.05, lwd = 1.5, col = alpha("darkgreen",.65))
arrows(LCLDIURNAL_DRY$hour, y0 = LCLDIURNAL_DRY$hpbl.avg -  LCLDIURNAL_DRY$hpbl.err , LCLDIURNAL_DRY$hour, y1= LCLDIURNAL_DRY$hpbl.avg   +  LCLDIURNAL_DRY$hpbl.err ,
       code=3, angle=90, length=0.05, lwd = 1.5, col = alpha("brown",.65))
arrows(LCLDIURNAL_WET$hour, y0 = LCLDIURNAL_WET$hpbl.avg  -  LCLDIURNAL_WET$hpbl.err , LCLDIURNAL_WET$hour, y1= LCLDIURNAL_WET$hpbl.avg   +  LCLDIURNAL_WET$hpbl.err ,
       code=3, angle=90, length=0.05, lwd = 1.5,col = alpha("darkgreen",.65))


title("Diurnal Average of LCL and ABL Heights")
subtitle="HARV,EMS, & NARR June-September 2017-2022"
mtext(subtitle)


plot(LCLDIURNAL_DRY$hour,LCLDIURNAL_DRY$LE.avg,col="blue", ylim = c(0,200),type="l",
     lwd =2, xlab = "Hour of Day", ylab = "Height [m]")
lines(LCLDIURNAL_WET$hour,LCLDIURNAL_WET$LE.avg, col="darkgreen",lwd=2, type="l")
lines(LCLDIURNAL_DRY$hour,LCLDIURNAL_DRY$H.avg, col="brown",type="l",lty = 2,lwd=2)
lines(LCLDIURNAL_WET$hour,LCLDIURNAL_WET$H.avg, col="darkgreen",lty=2,lwd=2, type="l")
legend('topleft', legend = c("ABL DRY","ABL WET", "LCL DRY","LCL WET"), lty = c(2,2,1,1),
       col = c("brown","darkgreen","brown","darkgreen"), lwd=2, bty ="n")
title("Diurnal Average of LCL and ABL Heights")
subtitle="HARV & EMS June-September 2017-2022"
mtext(subtitle)




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
subtitle="HARV,EMS, & NARR June-September 2017-2022"
mtext(subtitle)


















####SHADING TEST####
# Sample time series data
set.seed(123)
dates <- seq(as.Date("2022-01-01"), as.Date("2022-04-10"), by = "days")
values <- rnorm(length(dates), mean = 10, sd = 2)
ts_data <- data.frame(date = dates, value = values)

# Calculate mean and standard deviation
mean_value <- mean(ts_data$value)
sd_value <- sd(ts_data$value)

# Plot the time series line
plot(ts_data$date, ts_data$value, type = "l", col = "blue", lwd = 2, main = "Shaded Area around Average Line", xlab = "Date", ylab = "Value")

# Add the average line
abline(h = mean_value, col = "red", lwd = 2)

# Create x and y coordinates for the polygon
x_polygon <- c(ts_data$date, rev(ts_data$date))
y_polygon <- c(ts_data$value + sd_value, rev(ts_data$value - sd_value))

# Add the shaded area
polygon(x_polygon, y_polygon, col = alpha("red", 0.3), border = NA)






#Check correlation between LCL and soil moisture, check if ABL is primarily governed by USTAR

LCL_binned <- LCL_clean%>% mutate(swc_binned = cut(VSWCAnom.mean.med, breaks= 15))

binned <- LCL_binned  %>% group_by(swc_binned) %>% summarise(LCL.avg = mean(LCL,na.rm=TRUE),
                                                             LCL.sd = sd(LCL,na.rm=TRUE),
                                                             hpbl.avg = mean(hpbl,na.rm=TRUE),
                                                             hpbl.sd = sd(hpbl,na.rm=TRUE),
                                                             swc_avg= mean(VSWCAnom.mean.med, na.rm=TRUE),
                                                             freq = n())


binned$LCL.std.err <- binned$LCL.sd/(sqrt(binned$freq))

binned$hpbl.std.err <- binned$hpbl.sd/(sqrt(binned$freq))

#This plot needs error bars

plot(binned$swc_avg,binned$LCL.avg, type = "b", pch =2, col = "cadetblue", lwd=2 ,
     ylim = c(0,1100),
     xlab = "Soil Moisture Anomaly",
     ylab = "Height [m]")
lines(binned$swc_avg,binned$hpbl.avg, type ="b",lwd =2, col = "brown3")
arrows(binned$swc_avg,y0 = binned$hpbl.avg - binned$hpbl.std.err, binned$swc_avg, y1=binned$hpbl.avg  + binned$hpbl.std.err,
       code=3, angle=90, length=0.05, lwd = 1.5, col="darkgrey")
arrows(binned$swc_avg,y0 = binned$LCL.avg - binned$LCL.std.err, binned$swc_avg, y1=binned$LCL.avg  + binned$LCL.std.err,
       code=3, angle=90, length=0.05, lwd = 1.5, col="darkgrey")
legend("bottomleft",legend = c("ABL","LCL"),col = c("brown3","cadetblue4"), pch = c(1,2),
       pt.lwd = c(2,2), bty = "n")
title("ABL and LCL Heights by Soil Moisture Anomaly")
subtitle = "HARV & NARR Dataset, June - September 2016-2022"
mtext(subtitle)

####CROSSOVER EVENTS####
LCL_prec <- merge(LCL_clean,Ha1_ETfill_9222,by="datetime", all=TRUE)

LCL_prec_DRY <- LCL_prec %>% filter(LCL_prec$VSWCAnom.mean.med < 0)
LCL_prec_DRY <- LCL_prec_DRY%>% filter(LCL_prec_DRY$P > 0)
LCL_prec_WET <-LCL_prec %>% filter(LCL_prec$VSWCAnom.mean.med > 0)
LCL_prec_WET <- LCL_prec_WET%>% filter(LCL_prec_WET$P > 0)

#make a column of days, get mean of LCL and ABL, get sum of rain for that day
LCL_prec$date <- substr(as.character(LCL_prec$datetime), start = 1, stop = 10)

LCL_DAILY <- LCL_prec %>% group_by(date) %>% summarise(LCL.avg = mean(LCL, na.rm=TRUE),
                                                       LCL.sd = sd(LCL, na.rm=TRUE),
                                                       hpbl.avg = mean(hpbl, na.rm=TRUE),
                                                       hpbl.sd = sd(hpbl, na.rm=TRUE),
                                                       prec = sum(P, na.rm=TRUE),
                                                       swc.avg = mean(VSWCAnom.mean.med,na.rm=TRUE,
                                                       LE.avg = mean(LE.f, na.rm = TRUE),
                                                       H.avg = mean(H, na.rm=TRUE))
                                                       )



#do two graphs, on a histogram of rain patterns, and one bar chart comparing crossover

LCL_DAILY$CROSSOVER <- ifelse(LCL_DAILY$LCL.avg<LCL_DAILY$hpbl.avg, TRUE,FALSE)
LCL_DAILY$CROSSOVER_PREC <- ifelse(LCL_DAILY$CROSSOVER & LCL_DAILY$prec > 0, TRUE,FALSE)


LCL_DAILY_DRY <- LCL_DAILY %>% filter(swc.avg < 0)

LCL_DAILY_WET <- LCL_DAILY %>% filter(swc.avg > 0)

LCL_CROSS_DRY_percent <- sum(LCL_DAILY_DRY$CROSSOVER ==TRUE)/length(LCL_DAILY_DRY$CROSSOVER)
xLCL_CROSS_WET_percent <-  sum(LCL_DAILY_WET$CROSSOVER ==TRUE)/length(LCL_DAILY_WET$CROSSOVER)

LCL_CROSS_DRY_PREC_percent <- sum(LCL_DAILY_DRY$CROSSOVER_PREC ==TRUE)/length(LCL_DAILY_DRY$CROSSOVER_PREC)
LCL_CROSS_WET_PREC_percent <-  sum(LCL_DAILY_WET$CROSSOVER_PREC ==TRUE)/length(LCL_DAILY_WET$CROSSOVER_PREC)

#Comparing crossover percents 











#Histogram of rain patterns
LCL_DAILY_DRY_prec <- LCL_DAILY_DRY %>% filter(prec > 1)

LCL_DAILY_WET_prec <- LCL_DAILY_WET %>% filter(prec > 1)



w <- hist(LCL_DAILY_WET_prec$prec,breaks=20)                     # Store histogram info
w$density = w$counts/sum(w$counts)*100
d <- hist(LCL_DAILY_DRY_prec$prec,breaks = 10)                     # Store histogram info
d$density = d$counts/sum(d$counts)*100

plot(d,freq=FALSE, col =alpha("brown",.8), xlim = c(0,100),ylim = c(0,70), main="",
     xlab = "Precipitation [mm]")
par(new = TRUE) 
plot(w,freq=FALSE, xlim = c(0,100),ylim = c(0,70),col = alpha("darkgreen",.6),
     xlab = "",ylab ="",main="")
legend("topright", legend = c("Wet","Dry"), pch = 15, col = c("darkgreen","brown"),
       bty ="n")
title("Precipitation Density Histogram")
subtitle="HARV & EMS June-September 2017-2023"
mtext(subtitle)
par(new = TRUE) 

DRY_med <- mean(LCL_DAILY_DRY_prec$prec)
WET_med <- mean(LCL_DAILY_WET_prec$prec)

plot(5.4,50, xlim = c(0,100),ylim = c(0,70),pch=3,col ="brown",
     xlab = "", ylab ="",lwd=2)
par(new = TRUE) 
plot(12.2,50, xlim = c(0,100),ylim = c(0,70),pch=3,col = "darkgreen",
     xlab = "", ylab ="",lwd=2)
par(new = TRUE) 
plot(10.89,30, xlim = c(0,100),ylim = c(0,70),pch=4,col ="brown",
     xlab = "", ylab ="",lwd=2)
par(new = TRUE) 
plot(17.27,30, xlim = c(0,100),ylim = c(0,70),pch=4,col = "darkgreen",
     xlab = "", ylab ="",lwd=2)

text(9, 55, substitute(paste(bold('MEDIAN'))), col ="darkgrey")
text(14, 35, substitute(paste(bold('MEAN'))), col ="darkgrey")


