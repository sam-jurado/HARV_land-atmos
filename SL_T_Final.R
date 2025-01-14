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



########Up to 2023 Data#########


library(readr)
setwd("/Users/jurado/")
Ha1_ETfill_9222 <- read_csv("Downloads/ems_dailyCandETflux.csv")
df_prcp<- read_csv("Downloads/dailyrain_19642023.csv")


#separate into month, year


Ha1_ETfill_9222$month <- substr(as.character(Ha1_ETfill_9222$date), start = 6, stop = 7)
Ha1_ETfill_9222$month <- as.numeric(Ha1_ETfill_9222$month) 


df_prcp$month <- substr(as.character(df_prcp$date), start = 1, stop = 1)
df_prcp$month <- as.numeric(df_prcp$month)


#filter by year and month

df_prcp <- df_prcp %>% filter(df_prcp$month > 5)
df_prcp <- df_prcp %>% filter(df_prcp$month <10)
precip_EMS <- df_prcp  %>% group_by(year) %>% summarise(prec_sum = sum(prec_mm, na.rm =TRUE))

Ha1_ETfill_9222 <- Ha1_ETfill_9222 %>% filter(Ha1_ETfill_9222$month > 5)
Ha1_ETfill_9222 <- Ha1_ETfill_9222 %>% filter(Ha1_ETfill_9222$month <10)
ET_PRECIP<- Ha1_ETfill_9222 %>% group_by(wyear) %>% summarise(ET_sum = sum(ETd_mm, na.rm =TRUE))

colnames(ET_PRECIP)[colnames(ET_PRECIP) == "wyear"] <- "year"

# combine by date

ET_PRECIP <- merge(precip_EMS,ET_PRECIP, by= "year")

ET_PRECIP <- ET_PRECIP [-1, ]

ET_PRECIP$diff <- ET_PRECIP$prec_sum-ET_PRECIP$ET_sum
ET_PRECIP$col <- ifelse(ET_PRECIP$diff>=0,"cadetblue4","red3" )



###PLOT FINAL####


barplot(-ET_PRECIP$ET_sum,ET_PRECIP$year, ylim = c(-600,1000), col = "red", ylab = "Water Flux [mm]")
par(new=TRUE)
barplot(ET_PRECIP$prec_sum,ET_PRECIP$year, ylim = c(-600,1000), col = "cadetblue3")
par(new=TRUE)
barplot(ET_PRECIP$diff,ET_PRECIP$year, ylim = c(-600,1000),xlab = "Year", 
        col = ET_PRECIP$col, names.arg= ET_PRECIP$year)
legend("topleft", legend = c("Precipitation","Evaporation","+ Difference", "- Difference"),
       pch = 15, col = c("cadetblue3","red","cadetblue4","red3"), bty = "n")
title("Evaporative Deficit")
subtitle = "EMS June-September of 1992 - 2023"
mtext(subtitle)


install.packages("Kendall")
library(Kendall)
Kendall(ET_PRECIP$year,-ET_PRECIP$ET_sum)

ET_PRECIP_SURPLUS <- ET_PRECIP %>% filter(ET_PRECIP$col=="cadetblue4")
Kendall(ET_PRECIP_SURPLUS$year,ET_PRECIP_SURPLUS$diff)
summary(lm(ET_PRECIP_SURPLUS$diff~ET_PRECIP_SURPLUS$year))



############Correlation between LCL and Evaporative surplus / deficit


ET_PRECIP <- ET_PRECIP  %>% filter(year < 2023)
PBL_temp <- PBL %>% filter(year >1991 & year < 2023)
PBL_temp$dist_diff <- PBL_temp$hpbl - PBL_temp$lcl.mean

mean(ET_PRECIP$diff)


plot(PBL_temp$lcl.mean,ET_PRECIP$diff)
abline(h=0)
abline(h=120)
abline(v=590)
cor.test(PBL_temp$lcl.mean,ET_PRECIP$diff)








