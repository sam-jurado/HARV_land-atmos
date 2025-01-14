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
install.packages("ClimClass")



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
library(ClimClass)
library(geosphere)


library(readr)
hf300_03_monthly_m <- read_csv("~/Harvard_Forest/hf300_03_monthly_m.csv")
View(hf300_03_monthly_m)

series <- hf300_03_monthly_m

series <- series[c(1:6)]
#'Names in series columns must include: year, month, Tn and Tx 
#'(minimum and maximum temperatures, respectively) or, as an alternative,
#' Tm (mean temperatures), and P (mandatory)


clim_norm <- climate(series, first.yr = NULL, last.yr = NULL, max.perc.missing =25)


thornt_lst <- thornthwaite(series, latitude = 42.527, clim_norm = clim_norm, first.yr = 1964,
             last.yr = NULL, quant = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1),
             snow.init = 20, Tsnow = -1, TAW = 128, fr.sn.acc = 0.95,
             snow_melt_coeff = 1)

thornt_Evap <- thornt_lst$W_balance$Et0
thornt_Evap$avg <- rowMeans(thornt_Evap)
thornt_Evap$month <- seq(1,12,1)

thornt_P <- thornt_lst$W_balance$Precipitation
thornt_P$avg <- rowMeans(thornt_P)
thornt_P$month <- seq(1,12,1)

thornt_surp <- thornt_lst$W_balance$Surplus
thornt_surp$avg <- rowMeans(thornt_surp)
thornt_surp$month <- seq(1,12,1)
  
thornt_def <- thornt_lst$W_balance$Deficit
thornt_def$avg <- rowMeans(thornt_def)
thornt_def$month <- seq(1,12,1)

thornt_diff <- data.frame( P = thornt_P$avg, Evap = thornt_Evap$avg,
                           Diff = thornt_P$avg-thornt_Evap$avg) 




plot(thornt_Evap$month,thornt_Evap$avg, type = "b", col = "brown3", lwd=2, xlab = "Month",
     ylab = "Evap. and Precip [mm]", ylim = c(0,240))
lines(thornt_P$month,thornt_P$avg, type="b",pch = 2, col = "cadetblue3", lwd = 2)
segments(x0 = 9.24,
         x1 = 9.24,
         y0 = 80,
         y1 = 102,
         lwd = 2,
         col = "black",
         lty =2)

title("Thornthwaite Water Balance Diagram")
subtitle = "HARV 1973 - 2023"
mtext(subtitle)
text(7, 121, substitute(paste(bold(italic('Ut')))), col ="black")
text(8.9, 95, substitute(paste(bold(italic('Rc')))), col ="black")
text(11, 40, substitute(paste(bold(italic('Sr')))), col ="black")
text(2, 40, substitute(paste(bold(italic('Sr')))), col ="black")
text(10.5, 240, substitute(paste('AWHC = 128 mm')), col ="black")



#drought year 2022
thornt <- data.frame(month = seq(1,12,1),Evap = thornt_lst$W_balance$Et0$`2022`)
thornt$P <- thornt_lst$W_balance$Precipitation$`2022`
thornt$diff <-thornt$P-thornt$Evap 

plot(seq(1,12,1),thornt$Evap, type="b", ylim = c(0,240), col = "brown3", lwd =2,xlab = "Month",
     ylab = "Evap. and Precip [mm]")
lines(seq(1,12,1),thornt$P, type = "b", pch =2, col = "cadetblue3", lwd=2)
title("Thornthwaite Water Balance Diagram")
subtitle = "HARV 2022 Drought Year"
mtext(subtitle)
segments(x0 = 9,
         x1 = 9,
         y0 = 78,
         y1 = 180,
         lwd = 2,
         col = "black",
         lty = 2)
segments(x0 = 7.64,
         x1 = 7.64,
         y0 = 78,
         y1 = 124,
         lwd = 2,
         col = "black",
         lty = 2)
text(5.5, 83, substitute(paste(bold(italic('Ut')))), col ="black")
text(8.72, 103, substitute(paste(bold(italic('Rc')))), col ="black")
text(11, 60, substitute(paste(bold(italic('Sr')))), col ="black")
text(2, 40, substitute(paste(bold(italic('Sr')))), col ="black")
text(8,90, substitute(paste(bold(italic('Df')))), col ="black")
text(10.5, 240, substitute(paste('AWHC = 128 mm')), col ="black")


#Flood year 2023

thornt_evap <- thornt_lst$W_balance$Et0$`2023`
thornt_P <- thornt_lst$W_balance$Precipitation$`2023`

plot(seq(1,12,1),thornt_evap, type="b", ylim = c(0,240), col ="brown3",lwd =2,
     xlab = "Month", ylab = "Evap. and Precip. [mm]")
lines(seq(1,12,1),thornt_P, type ="b", pch =2, col = "cadetblue3",lwd=2)
title("Thornthwaite Water Balance Diagram")
subtitle = "HARV 2023 Flood Year"
mtext(subtitle)
text(5, 87, substitute(paste(bold(italic('Ut')))), col ="black")
text(9, 125, substitute(paste(bold(italic('Sr')))), col ="black")
text(2, 30, substitute(paste(bold(italic('Sr')))), col ="black")
text(11, 240, substitute(paste('AWHC = 128 mm')), col ="black")

######What are the flood years and the dry years in the past 60 years?#####

library(readr)
hf300_01_annual_m <- read_csv("~/Harvard_Forest/hf300-01-annual-m.csv")

#top flood: 2018,2011,2008,2023,2005,2021


#Average Flood year Thornwaithe

thornt_lst <- thornthwaite(series, latitude = 42.527, clim_norm = clim_norm, first.yr = 1973,
                           last.yr = NULL, quant = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1),
                           snow.init = 20, Tsnow = -1, TAW = 128, fr.sn.acc = 0.95,
                           snow_melt_coeff = 1)


thornt_Evap <- thornt_lst$W_balance$Et0
thornt_Evap  <- thornt_Evap[c(46,39,36,51,33,49)]
thornt_Evap$avg <- rowMeans(thornt_Evap)
thornt_Evap$month <- seq(1,12,1)

thornt_P <- thornt_lst$W_balance$Precipitation
thornt_P  <- thornt_P [c(46,39,36,51,33,49)]
thornt_P$avg <- rowMeans(thornt_P)
thornt_P$month <- seq(1,12,1)



plot(thornt_Evap$month,thornt_Evap$avg, type = "b", col = "brown3", lwd=2, xlab = "Month",
     ylab = "Evap. and Precip [mm]", ylim = c(0,240))
lines(thornt_P$month,thornt_P$avg, type="b",pch = 2, col = "cadetblue3", lwd = 2)
title("Thornthwaite Water Balance Diagram")
subtitle = "HARV Highest 10% Flood Years"
mtext(subtitle)

text(7, 155, substitute(paste(bold(italic('Sr')))), col ="black")
text(10.5, 240, substitute(paste('AWHC = 128 mm')), col ="black")



#Average Dry Years 
#Driest Years: 1965,1964,1966,1980,1970,1968


thornt_Evap  <- thornt_Evap[c(1,2,3,5,7,17)]
thornt_Evap$avg <- rowMeans(thornt_Evap)
thornt_Evap$month <- seq(1,12,1)

thornt_P <- thornt_lst$W_balance$Precipitation
thornt_P  <- thornt_P[c(1,2,3,5,7,17)]
thornt_P$avg <- rowMeans(thornt_P)
thornt_P$month <- seq(1,12,1)

thornt_diff <- data.frame(  month = seq(1,12,1),diff = thornt_P$avg - thornt_Evap$avg)

plot(thornt_Evap$month,thornt_Evap$avg, type = "b", col = "brown3", lwd=2, xlab = "Month",
     ylab = "Evap. and Precip [mm]", ylim = c(0,240))
lines(thornt_P$month,thornt_P$avg, type="b",pch = 2, col = "cadetblue3", lwd = 2)
title("Thornthwaite Water Balance Diagram")
subtitle = "HARV Lowest 10% Drought Years"
mtext(subtitle)


segments(x0 = 11.75,
         x1 = 11.75,
         y0 = 70,
         y1 = 5,
         lwd = 2,
         col = "black",
         lty = 2)
segments(x0 = 7.18,
         x1 = 7.18,
         y0 = 125,
         y1 = 50,
         lwd = 2,
         col = "black",
         lty = 2)


text(11, 45, substitute(paste(bold(italic('Rc')))), col ="black")
text(6.5, 100, substitute(paste(bold(italic('Ut')))), col ="black")
text(7.75, 80, substitute(paste(bold(italic('Df')))), col ="black")
text(2, 25, substitute(paste(bold(italic('Sr')))), col ="black")

text(10.5, 240, substitute(paste('AWHC = 128 mm')), col ="black")
































