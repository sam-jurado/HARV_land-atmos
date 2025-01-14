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

setwd("/Users/jurado/Harvard_Forest")
hf300_03_monthly_m <- read_csv("hf300_03_monthly_m.csv")
View(hf300_03_monthly_m)

#'What is the average amount of rain received from months 3-5, and subtract from
#'potential evapotranspiration to get how much ground water is being utilized from surplus
#'Do the same for months 6-8 for each year. Plot precip - evap for 3-5 on the x axis,
#'and AWHC - (precip - evap) for 6-8 on the y-axis to see if the drier springs have
#'drier summers. Color dots that run into deficiet (AWHC - (precip - evap) < 0) 



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

df_prec<- thornt_lst$W_balance$Precipitation
df_evap<- thornt_lst$W_balance$Et0

prec_list <- c()
df_spring <- data.frame(year = seq(1964,2023,1))
for (x in 1:length(df_prec)){
  prec_sum <- sum(df_prec[5,x])
  prec_list <- append(prec_sum,prec_list)
}

prec_list <- rev(prec_list)

#precipitation in the month of May
df_spring$prec <- prec_list

prec_list <- c()
for (x in 1:length(df_prec)){
  prec_sum <- sum(df_prec[6:8,x])
  prec_list <- append(prec_sum,prec_list)
}
prec_list <- rev(prec_list)

#precipitation throughout the months of June, July, and August
df_spring$prec_summer <- prec_list

prec_list <- c()
for (x in 1:length(df_evap)){
  prec_sum <- sum(df_evap[5,x])
  prec_list <- append(prec_sum,prec_list)
}
prec_list <- rev(prec_list)

#evaporation in the month of May
df_spring$evap_spring <- prec_list


prec_list <- c()
for (x in 1:length(df_evap)){
  prec_sum <- sum(df_evap[6:8,x])
  prec_list <- append(prec_sum,prec_list)
}
prec_list <- rev(prec_list)
#Total evaporation in the months of June,July,and August
df_spring$evap_summer <- prec_list

#Net water gained in May
df_spring$diff_spring = df_spring$prec-df_spring$evap_spring

#Net water gained in the summer
df_spring$diff_summer = df_spring$prec_summer- df_spring$evap_summer

#'Spring starts at 128mm of AWHC due to winter surplus
#'This 128mm is added to the evaporative loss of may, when evaporation overcomes precip for first time in year
#'Thus, this is the start of summer AW
df_spring$summer_AW = 128 + df_spring$diff_spring

#If above 128, water percolates out, so set to 128, otherwise keep value
df_spring$summer_AW <- ifelse(df_spring$summer_AW >128,128,df_spring$summer_AW)

#Summer AW is then combined with the net water gained (or lost) from the summer to get end of summer AW
df_spring$summer_AW = df_spring$summer_AW+df_spring$diff_summer

#This is all of the springs in which rain outpaced evaporation
df_spring_surplus <- df_spring %>% filter(df_spring$summer_AW > 128)


#This is the final summer AW value, all values above 128 are capped at 128 since water percolates out
df_spring$summer_AW <- ifelse(df_spring$summer_AW >128,128,df_spring$summer_AW)




plot(df_spring$diff_spring,df_spring$summer_AW)
abline(h=0)
abline(v=12) #spring mean


quantile(df_spring$diff_spring)
mean(df_spring$diff_spring)
median(df_spring$diff_spring)

#Dry springs are defined as springs below the mean (1964-2023 AW value (precip - evap of May))
#Wet Springs are defined as springs above the mean  (1964-2023 AW value (precip - evap of May))
Dry_spring <- df_spring %>% filter(df_spring$diff_spring<12.095)
Wet_spring <- df_spring %>% filter(df_spring$diff_spring>12.095)


12/34
#35.29% of summers following anomalously dry springs are in water deficit

15/60
#25% of summers are in deficit 

12/15
#80% of all summers that were in water deficit followed an anomalously dry spring

####Are temperatures increasing?

##comparison between dry springs and wet springs
12/34
#35.29% of summers following anomalously dry springs are in water deficit

4/26
#15.38% of summers following anomalously wet springs are in water deficit

#Median wet spring summer AWHC
median(Wet_spring$summer_AW)
#Mean wet spring summer AWHC
mean(Wet_spring$summer_AW)
sd(Wet_spring$summer_AW)/sqrt(26)


#Median all summer AWHC
median(df_spring$summer_AW)
#Mean all summer AWHC
mean(df_spring$summer_AW)
sd(df_spring$summer_AW)/sqrt(60)


#Median dry spring summer AWHC
median(Dry_spring$summer_AW)
#Mean dry spring summer AWHC
mean(Dry_spring$summer_AW)
sd(Dry_spring$summer_AW)/sqrt(34)



error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

counts <- c(mean(Wet_spring$summer_AW),mean(df_spring$summer_AW),mean(Dry_spring$summer_AW))
names(counts) <- c("Wet","Mean","Dry")
barplot(counts, main="Summer Available Water Column by Spring Moisture Regime",
        xlab="Spring Available Water",
        ylab ="Available Water [mm]" ,
        ylim = c(0,80))
subtitle = "Fisher Station 1964 - 2023"
mtext(subtitle
      )


install.packages("ggpubr")
library(ggpubr)


df_spring$regime <- ifelse(df_spring$diff_spring < 12.095,"Dry","Wet")
df_spring_deficit <- df_spring %>% filter(df_spring$summer_AW < 0)
df_spring$summer_AW<- ifelse(df_spring$summer_AW < 0,0,df_spring$summer_AW)


# Create a simple bar plot


ggbarplot(
  df_spring, x = "regime", y = "summer_AW", 
  add = c("mean_se", "jitter"),
  fill = "darkorange",
  xlab = "Spring Moisture Regime",
  ylab = "Summer Available Water Column [mm]",
  ylim = c(-100,250))+
  ggtitle(label= "Available Water Column by Spring Moisture Regime", subtitle ="Fisher Station 1964 - 2023")+
  theme(plot.margin = margin(1,1,3,1, "cm"))+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle= element_text(hjust = .5),)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey")  # Add horizontal line at y = 0

data1 <- subset(df_spring, regime == "Dry" )
data2 <- subset(df_spring, regime == "Wet" )

result <- t.test(data1$summer_AW, data2$summer_AW)
print(result)
#p = .058

ggbarplot(
  df_spring_deficit, x = "regime", y = "summer_AW", 
  add = c("mean_se", "jitter"),
  fill = "brown",
  xlab = "Spring Moisture Regime",
  ylab = "Summer Water Column Deficit [mm]",
  ylim = c(-100,250))+
  ggtitle(label= "Available Water Column by Spring Moisture Regime", subtitle ="Fisher Station 1964 - 2023")+
  theme(plot.margin = margin(1,1,3,1, "cm"))+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle= element_text(hjust = .5),)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey")  # Add horizontal line at y = 0


mean <- df_spring_deficit %>% filter(regime == "Wet")
mean(mean$summer_AW)


df_spring_surplus <- df_spring_surplus[order(df_spring_surplus$summer_AW),]

df_spring_surplus$regime <- ifelse(df_spring_surplus$diff_spring < 12.095,"Dry","Wet")
df_spring_surplus$summer_AW <-df_spring_surplus$summer_AW 

ggbarplot(
  df_spring_surplus, x = "regime", y = "summer_AW", 
  add = c("mean_se", "jitter"),
  fill = "darkgreen",
  xlab = "Spring Moisture Regime",
  ylab = "Summer Water Column Surplus [mm]",
  ylim = c(-100,350))+
  ggtitle(label= "Available Water Column by Spring Moisture Regime", subtitle ="Fisher Station 1964 - 2023")+
  theme(plot.margin = margin(1,1,3,1, "cm"))+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle= element_text(hjust = .5),)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey")  # Add horizontal line at y = 0


sd <- df_spring_deficit %>% filter(regime == "Wet")
sd(sd$summer_AW, na.rm=TRUE)


surplus_wet <- df_spring_surplus %>% filter(regime == "Wet")
surplus_dry <- df_spring_surplus %>% filter(regime == "Dry")

deficit_wet <- df_spring_deficit %>% filter(regime == "Wet")
deficit_dry <- df_spring_deficit %>% filter(regime == "Dry")



surplus_ttest <- t.test(surplus_wet$summer_AW, surplus_dry$summer_AW)
print(surplus_ttest)

deficit_ttest <- t.test(deficit_wet$summer_AW, deficit_dry$summer_AW)
print(deficit_ttest)
####Deficit wet -58.8, surplus dry 240.8, 337.6









