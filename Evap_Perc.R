install.packages('BiocManager', dependencies = TRUE)
install.packages("neonUtilities", dependencies = TRUE)
install.packages("readr", dependencies = TRUE)
BiocManager::install('rhdf5')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
install.packages('dplyr')
install.packages("oce")
install.packages("tidyr")
install.packages("lubridate")
install.packages("naniar")
install.packages("Hmisc")
install.packages("FREddyPro")


library(neonUtilities)
library(readr)
library(BiocManager)
library(rhdf5)
library(dplyr)
library(oce)
library(tidyr)
library(lubridate)
library(naniar)
library(ggplot2)
library(stringr)
library(limma)
library(zoo)
library(stringr)




#function that takes lev filepath as a string and return a dataframe of that file
lev20_reader <- function(filepath){
  watercolumn <- read_csv(filepath)
  watercolumn <- data.frame(watercolumn)
  watercolumn <- tail(watercolumn, -5)
  watercolumn  <- data.frame(str_split_fixed(watercolumn$AERONET.Version.3., ',', 80))
  names(watercolumn) <- watercolumn[1,]
  watercolumn  <- watercolumn[-1,]
  return(watercolumn)
}


file_list<- list.files("~/Harvard_Forest/", pattern=".lev20", all.files=FALSE,
          full.names=TRUE)

df_EVAP <- data.frame()
for (x in 1:length(file_list)){
  df <- lev20_reader(file_list[x])
  df_EVAP <- rbind(df_EVAP,df)
}                                         
df_EVAP$day <- as.numeric(substr(as.character(df_EVAP$`Date(dd:mm:yyyy)`), start = 1, stop = 2))
df_EVAP$month <- as.numeric(substr(as.character(df_EVAP$`Date(dd:mm:yyyy)`), start = 4, stop = 5))
df_EVAP$year <- as.numeric(substr(as.character(df_EVAP$`Date(dd:mm:yyyy)`), start = 7, stop = 10))

df_EVAP$date <- paste(as.character(df_EVAP$year),as.character(df_EVAP$month), sep = "-0")
df_EVAP$date <- paste(as.character(df_EVAP$date),as.character(df_EVAP$day), sep = "-")
df_EVAP$date <- as.POSIXct(df_EVAP$date, format = "%Y-%m-%d", tz = "UTC")


library(readr)
setwd("/Users/jurado/")
Ha1_ETfill_9222 <- read_csv("Downloads/Ha1_LEfill_9222.csv")



setwd("/Users/jurado/")
EMS <- read_csv("Downloads/ems_hourlygapfilled_9123.csv")



#translating into EST
EMS$DateTime<- as.POSIXct(EMS$DateTime, tz ="EST5EDT,M3.2.0/2:00:00,M11.1.0/2:00:00")


EMS$year <- substr(as.character(EMS$DateTime), start = 1, stop = 4)
EMS$year <- as.numeric(EMS$year) 

EMS$month <- substr(as.character(EMS$DateTime), start = 6, stop = 7)
EMS$month <- as.numeric(EMS$month) 

EMS$day <- substr(as.character(EMS$DateTime), start = 9, stop = 10)
EMS$day <- as.numeric(EMS$day) 

EMS$hour <- substr(as.character(EMS$DateTime), start = 12, stop = 13)
EMS$hour <- as.numeric(EMS$hour) 
EMS$hour<-replace(EMS$hour, is.na(EMS$hour), 0) 


EMS <- EMS %>% filter(EMS$month > 5)

EMS <- EMS %>% filter(EMS$month <10)
EMS <- EMS%>% filter(EMS$year >2016)

EMS$date <- paste(as.character(EMS$year),as.character(EMS$month), sep = "-0")
EMS$date <- paste(as.character(EMS$date),as.character(EMS$day), sep = "-")


#for conversion,
EMS$ET_mm <- EMS$LE_f*(1*10**-6)*(1/2.26)*1000*(1/10000)*3600*10



#making time stamps
#Ha1_ETfill_9222$year <- substr(as.character(Ha1_ETfill_9222$TIMESTAMP_END), start = 1, stop = 4)
#Ha1_ETfill_9222$year <- as.numeric(Ha1_ETfill_9222$year) 
#Ha1_ETfill_9222$month <- substr(as.character(Ha1_ETfill_9222$TIMESTAMP_END), start = 5, stop = 6)
#Ha1_ETfill_9222$month <- as.numeric(Ha1_ETfill_9222$month) 
#Ha1_ETfill_9222$day <- substr(as.character(Ha1_ETfill_9222$TIMESTAMP_END), start = 7, stop = 8)
#Ha1_ETfill_9222$day <- as.numeric(Ha1_ETfill_9222$day) 
#Ha1_ETfill_9222$hour <- substr(as.character(Ha1_ETfill_9222$TIMESTAMP_END), start = 9, stop = 10)
#Ha1_ETfill_9222$hour <- as.numeric(Ha1_ETfill_9222$hour) 

#Ha1_ETfill_9222$date <- paste(as.character(Ha1_ETfill_9222$year),as.character(Ha1_ETfill_9222$month), sep = "-0")
#Ha1_ETfill_9222$date <- paste(as.character(Ha1_ETfill_9222$date),as.character(Ha1_ETfill_9222$day), sep = "-")

#Ha1_ETfill_9222 <- Ha1_ETfill_9222 %>% filter(Ha1_ETfill_9222$month > 5)

#Ha1_ETfill_9222 <- Ha1_ETfill_9222 %>% filter(Ha1_ETfill_9222$month <10)
#Ha1_ETfill_9222 <- Ha1_ETfill_9222 %>% filter(Ha1_ETfill_9222$year >2016)

#Ha1_ETfill_9222$ET_mm <- Ha1_ETfill_9222$LE_f*(1*10**-6)*(1/2.26)*1000*(1/10000)*3600*10




#turn all zeros into NA

ET_PRECIP<- EMS %>% group_by(date) %>% summarise(ET_sum = sum(ET_mm, na.rm =TRUE),
                                                             VPD = mean(VPD_f,na.rm=TRUE))
ET_PRECIP$ET_sum  <- ifelse(ET_PRECIP$ET_sum == 0,NA,ET_PRECIP$ET_sum)

ET_PRECIP$date <- as.POSIXct(ET_PRECIP$date, format = "%Y-%m-%d", tz = "EST")

ET_PRECIP_COL <- merge(ET_PRECIP,df_EVAP, by ="date", all=TRUE)

ET_PRECIP_COL$`Precipitable_Water(cm)` <- as.numeric(ET_PRECIP_COL$`Precipitable_Water(cm)`)

mean((ET_PRECIP_COL$ET_sum)/(ET_PRECIP_COL$`Precipitable_Water(cm)`), na.rm=TRUE)

#on average, land gives 17% of water to column, the rest is advected.


#merge with soil moisture, must use df_anom_daily from Soil_Moisture.R


ET_PRECIP_COL_SM <- merge(ET_PRECIP_COL,df_anom_daily, by = "date", all=TRUE)
ET_PRECIP_COL_SM$perc_evap <- ET_PRECIP_COL_SM$ET_sum/ET_PRECIP_COL_SM$`Precipitable_Water(cm)`*10
ET_PRECIP_COL_SM <- data.frame(date = ET_PRECIP_COL_SM$date,
                               ET = ET_PRECIP_COL_SM$ET_sum,
                               Precip_col = ET_PRECIP_COL_SM$`Precipitable_Water(cm)`*10,
                               VSWC.mean = ET_PRECIP_COL_SM$VSWC.mean,
                               VPD = ET_PRECIP_COL_SM$VPD)
ET_PRECIP_COL_SM$perc_evap <- (ET_PRECIP_COL_SM$ET)/(ET_PRECIP_COL_SM$Precip_col)*100











rbPal <- colorRampPalette(c('blue','red'))
ET_PRECIP_COL_SM$Col <- rbPal(10)[as.numeric(cut(ET_PRECIP_COL_SM$VPD,breaks = 15))]

ET_PRECIP_COL_SM$pch <- NA
ET_PRECIP_COL_SM$pch <- ifelse(ET_PRECIP_COL_SM$VPD > .8501, 17,ET_PRECIP_COL_SM$pch )


mean(ET_PRECIP_COL_SM$VPD,na.rm=TRUE)
sd(ET_PRECIP_COL_SM$VPD,na.rm=TRUE)
ET_PRECIP_COL_SM$VPD %>% quantile(prob=c(.25,.5,.75,.90,.95), type=1,na.rm = TRUE)
#avg RH = 75.23
#sd = 12.8667
#25% = 66.2
#50% = 75.6
#75% = 84.5

#avg VPD = .5975
#sd = .34
#25% = .3449860
#50% = .592322
#75% = .8501


ET_PRECIP_COL_SM$VSWC.mean %>% quantile(prob=c(.25,.5,.75,.90,.95), type=1,na.rm = TRUE)


VPD_High <- ET_PRECIP_COL_SM %>% filter(VPD > .8501 )
VPD_High <- na.omit(VPD_High)


#RH_High <- RH_High %>% filter(VSWC.mean <.06 )

VPD_Avg <- ET_PRECIP_COL_SM %>% filter(VPD > .3449860)
VPD_Avg <- VPD_Avg%>% filter( VPD <.8501)
VPD_Avg <- na.omit(VPD_Avg)

ET_PRECIP_COL_SM$pch <- ifelse(ET_PRECIP_COL_SM$VPD < .8501 & ET_PRECIP_COL_SM$VPD > .3449860,19,ET_PRECIP_COL_SM$pch)
ET_PRECIP_COL_SM$pch <- ifelse(ET_PRECIP_COL_SM$VPD < .3449860,15 ,ET_PRECIP_COL_SM$pch)


VPD_Low <- ET_PRECIP_COL_SM %>% filter(VPD < .3449860)
VPD_Low <- na.omit(VPD_Low)




reg_high <- lm(VPD_High$perc_evap~VPD_High$VSWC.mean)
reg_avg <- lm(VPD_Avg$perc_evap~VPD_Avg$VSWC.mean)
reg_low <- lm(VPD_Low$perc_evap~VPD_Low$VSWC.mean)
mean <- lm(ET_PRECIP_COL_SM$perc_evap~ET_PRECIP_COL_SM$VSWC.mean)


# Load necessary library for add.color.bar (e.g., raster or custom function)
install.packages('phytools')
library(phytools) # `add.color.bar()` is often included in the `raster` package

# Define color ramp palette and assign colors based on VPD
rbPal <- colorRampPalette(c('blue', 'red'))
ET_PRECIP_COL_SM$Col <- rbPal(10)[as.numeric(cut(ET_PRECIP_COL_SM$VPD, breaks = 15))]

# Base plot
plot(ET_PRECIP_COL_SM$VSWC.mean, ET_PRECIP_COL_SM$perc_evap, 
     col = alpha(ET_PRECIP_COL_SM$Col, 0.8), 
     pch = ET_PRECIP_COL_SM$pch, 
     xlab = "Soil Moisture Anomaly", 
     ylab = "% Water from ET", 
     ylim = c(0, 70))

# Add regression lines
abline(reg_high, lwd = 2, col = "red")
abline(reg_avg, lwd = 2, col = "purple")
abline(reg_low, lwd = 2, col = "blue")

# Add text annotations
text(0.08, 47, substitute(paste('75% VPD ')), col = "red")
text(0.1, 30, substitute(paste('IQR VPD ')), col = "purple")
text(0.12, 18, substitute(paste('25% VPD ')), col = "blue")

# Add title and subtitle
title("Precipitable Water from Evapotranspiration")
subtitle <- "HARV, EMS June-September 2017-2023"
mtext(subtitle)

# Add color bar
vpd_min <- min(ET_PRECIP_COL_SM$VPD, na.rm = TRUE)
vpd_max <- max(ET_PRECIP_COL_SM$VPD, na.rm = TRUE)
vpd_range <- c(vpd_min, vpd_max)

add.color.bar(
  leg = 0.15,                      # Position of the legend (adjust as needed)
  cols = rbPal(10),               # Use the defined color palette
  title = "Vapor Pressure Deficit [kPa]",           # Title of the color bar
  lims = vpd_range,               # Limits of the VPD variable
  digits = 1,                     # Number of decimal places for labels
  lwd = 4,                        # Width of the color bar
  outline = TRUE                  # Outline around the color bar
)




mean(VPD_Low$perc_evap, na.rm=TRUE)
mean(VPD_High$perc_evap, na.rm=TRUE)

summary(mean)

library(Kendall)
Kendall(VPD_High$VSWC.mean, VPD_High$perc_evap)
summary(reg_high)
Kendall(VPD_Avg$VSWC.mean, VPD_Avg$perc_evap)
summary(reg_avg)
Kendall(VPD_Low$VSWC.mean, VPD_Low$perc_evap)
summary(reg_low)


mean(ET_PRECIP_COL_SM$perc_evap, na.rm=TRUE)


#'Add 3 lines, a line that shows increase with dry air conditions, average air conditions,
#'and with wet air conditions

#'We need to consider both temperature and how much water was already in the column




plot(VPD_Low$VSWC.mean,VPD_Low$perc_evap, ylim = c(0,.7))

plot(VPD_Avg$VSWC.mean,VPD_Avg$perc_evap, ylim = c(0,.7))

plot(VPD_High$VSWC.mean,VPD_High$perc_evap, ylim = c(0,.7))














                                         
                                         
                                         
                                         
                                         
                                         
                                         