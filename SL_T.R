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



#######DATA ACQUISITION#####

setwd("/Users/jurado/Harvard_Forest")


#For downloading many months
zipsByProduct(dpID="DP4.00200.001", site="HARV",
              startdate="2023-08", enddate="2023-08", package="basic", check.size=FALSE)

#For downloading one month
getPackage("DP4.00200.001", site_code="HARV", 
           year_month="2023-09", package="basic")

#Must manually unzip file and compile .gz file

folder = "NEON.D01.HARV.DP4.00200.001.2022-06.basic.20230127T120753Z.RELEASE-2023"
file = "NEON.D01.HARV.DP4.00200.001.nsae.2022-06.basic.20221216T041952Z.h5"

wd = paste("/Users/jurado/Harvard_Forest/filesToStack00200",folder,sep = "/") 

setwd(wd)
h5ls(file)


#Reading In Data

LE <- h5read(file , "/HARV/dp04/data/fluxH2o/nsae")
Lwin <-h5read(file , "/HARV/dp01/data/radiNet/000_060_30m/radiLwIn")
Lwout <-h5read(file , "/HARV/dp01/data/radiNet/000_060_30m/radiLwOut")
Swin <-h5read(file , "/HARV/dp01/data/radiNet/000_060_30m/radiSwIn")
Swout <-h5read(file , "/HARV/dp01/data/radiNet/000_060_30m/radiSwOut")

#evaportive fraction
EF <- data.frame(Lwin$timeEnd, Lwin$mean-Lwout$mean+Swin$mean-Swout$mean, LE$flux, LE$flux/(Lwin$mean-Lwout$mean+Swin$mean-Swout$mean))
colnames(EF)[1] ="datetime"
colnames(EF)[2] ="Rnet"
colnames(EF)[3] ="LE"
colnames(EF)[4] ="EF"

rm(LE,Lwin,Lwout,Swin,Swout)

setwd("/Users/jurado/Harvard_Forest")
df_sm <- read.csv("hf206-03-HK-understory.csv")

df_sm$datetime <- anydate(df_sm$datetime) 
df_sm <- df_sm[df_sm$datetime >= "2022-06-01" &   df_sm$datetime <= "2022-06-30", ]
df_sm$swc_avg <- rowMeans(df_sm[ , c(18,19,31,32)], na.rm=TRUE)

EF$soilm <- df_sm$swc_avg


####Despike#####

EF$Rnet <- replace(EF$Rnet, abs(EF$Rnet) > 1000, NA)
EF$LE <- replace(EF$LE, abs(EF$LE) > 500, NA)
EF$EF <- replace(EF$EF, abs(EF$EF) > 2, NA)
EF$Rnet <- replace(EF$Rnet, EF$Rnet < 0, NA)


EF$Rnet <- runmed(EF$Rnet,k=5)
EF$LE <- runmed(EF$LE,k=5)
EF$EF <- runmed(EF$EF,k=5)


#Do diurnal Rnet and LE 
  

#####Datetime stuff####

EF$datetime <- as.POSIXct(EF$datetime,format="%Y-%m-%dT%H:%M", tz= "UTC")
EF$datetime <- as.POSIXct(EF$datetime,format="%Y-%m-%dT%H:%M", tz= "EST")

library(stringr)
split <- EF$datetime %>% str_split_fixed(" ",2)
split <- data.frame(split)

#rm(EF$datetime)
EF$date <- split$X1
EF$time <- split$X2

  
dfDAILY <- EF %>% group_by(time) %>% summarise(EF.avg = mean(EF, na.rm = TRUE),
                                               Rnet.avg = mean(Rnet, na.rm = TRUE),
                                               LE.avg = mean(LE, na.rm = TRUE))


par(mar = c(5, 4, 4, 4) + 0.3) 
plot(seq(0, 23.5, .5),dfDAILY$Rnet.avg, type="l",lwd=1, xlab = "Time of Day [Hr]",
     ylab = expression(paste("Flux ", "[","W/"*m^2,"]"))
)
lines(seq(0, 23.5, .5),dfDAILY$LE.avg, type="l",col="green")

par(new = TRUE)                             # Add new plot
plot(seq(0, 23.5, .5), dfDAILY$EF.avg, type = "l", col = "red",              # Create second plot without axes
     axes = FALSE, xlab = "", ylab = "",ylim = c(0, 1))
legend(15,1,legend=c("Net Radiation","Latent Heat","Evaporative Fraction"), col=c("black","green","red"),
       lty=c(1,1,1), ncol=1,cex=.75)
axis(side = 4, at = pretty(seq(0, 1, .1)))     # Add second axis
mtext("Evaporative Fraction", side = 4, line = 3) 
title(expression(paste("Surface Energy Balance")))
subtitle = "HARV 06/01/22 - 06/30/22"
mtext(subtitle)



#Hours of 4-20 is essentially day time I will use that










###########Data Scraper###########


folder_list = c("NEON.D01.HARV.DP4.00200.001.2017-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2017-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2017-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2017-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2018-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2018-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2018-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2018-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2019-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2019-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2019-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2019-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2020-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2020-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2020-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2020-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2021-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2021-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2021-08.basic.20230127T120753Z.RELEASE-2023 2",
                "NEON.D01.HARV.DP4.00200.001.2021-09.basic.20230127T120753Z.RELEASE-2023 2",
                "NEON.D01.HARV.DP4.00200.001.2022-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2022-07.basic.20220810T224035Z.PROVISIONAL",
                "NEON.D01.HARV.DP4.00200.001.2022-08.basic.20230225T225559Z.PROVISIONAL",
                "NEON.D01.HARV.DP4.00200.001.2022-09.basic.20230225T225624Z.PROVISIONAL",
                "NEON.D01.HARV.DP4.00200.001.2023-06.basic.20230713T225900Z.PROVISIONAL",
                "NEON.D01.HARV.DP4.00200.001.2023-07.basic.20230810T213732Z.PROVISIONAL")

file_list = c("NEON.D01.HARV.DP4.00200.001.nsae.2017-06.basic.20221215T022303Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2017-07.basic.20221215T023843Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2017-08.basic.20221215T025413Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2017-09.basic.20221215T030921Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2018-06.basic.20221215T090448Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2018-07.basic.20221215T092218Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2018-08.basic.20221215T093710Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2018-09.basic.20221215T095218Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2019-06.basic.20221215T121219Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2019-07.basic.20221215T122934Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2019-08.basic.20221215T124633Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2019-09.basic.20221215T130324Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2020-06.basic.20221215T153201Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2020-07.basic.20221215T154738Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2020-08.basic.20221215T160708Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2020-09.basic.20221215T162259Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2021-06.basic.20221216T005046Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2021-07.basic.20221216T010713Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2021-08.basic.20221216T012510Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2021-09.basic.20221216T014228Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2022-06.basic.20221216T041952Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2022-07.basic.20220809T001307Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2022-08.basic.20230224T000842Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2022-09.basic.20230224T002819Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2023-06.basic.20230712T165914Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2023-07.basic.20230809T174858Z.h5")


df <- data.frame()

for (x in 1:length(folder_list)){
  wd = paste("/Users/jurado/Harvard_Forest/filesToStack00200",folder_list[x],sep = "/")
  setwd(wd)
  h5ls(file_list[x])
  
  #Reading In Data
  
  LE <- h5read(file_list[x] , "/HARV/dp04/data/fluxH2o/nsae")
  Lwin <-h5read(file_list[x] , "/HARV/dp01/data/radiNet/000_060_30m/radiLwIn")
  Lwout <-h5read(file_list[x] , "/HARV/dp01/data/radiNet/000_060_30m/radiLwOut")
  Swin <-h5read(file_list[x] , "/HARV/dp01/data/radiNet/000_060_30m/radiSwIn")
  Swout <-h5read(file_list[x] , "/HARV/dp01/data/radiNet/000_060_30m/radiSwOut")
  
  EF <- data.frame(Lwin$timeEnd, Lwin$mean-Lwout$mean+Swin$mean-Swout$mean, LE$flux, LE$flux/(Lwin$mean-Lwout$mean+Swin$mean-Swout$mean))
  colnames(EF)[1] ="datetime"
  colnames(EF)[2] ="Rnet"
  colnames(EF)[3] ="LE"
  colnames(EF)[4] ="EF"
  
  rm(LE,Lwin,Lwout,Swin,Swout)

  df <- rbind(df,EF)
  
  print(wd)
  
}

##############################NEON RELATIVE HUMIDITY###########################
zipsByProduct(dpID="DP1.00098.001", site="HARV",
              startdate="2017-06", enddate="2023-09", package="basic", check.size=FALSE)


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

#in celsius
df_RH <- df_RH %>% filter(RHFinalQF == 0 )
df_RH$VPSat <- (610.7*10**(7.5*df_RH$tempRHMean/(237.3+df_RH$tempRHMean)))/1000

df_RH$VPair <- (610.7*10**(7.5*df_RH$tempRHMean/(237.3+df_RH$tempRHMean)))/1000 * (df_RH$RHMean/100)

df_RH$VPD <- df_RH$VPSat - df_RH$VPair 

df_RH$datetime <- as.POSIXct(df_RH$endDateTime, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC") 
df_RH$datetime<- as.POSIXct(df_RH$datetime, format = "%Y-%m-%dT%H:%M:%S", tz = "EST") 
plot(df_RH$datetime,df_RH$VPD,type="l")






##############################NEON SOIL MOISTURE###############################


####################################################
setwd("/Users/jurado/Harvard_Forest")
getPackage("DP1.00094.001", site_code="HARV", 
           year_month="2023-09", package="basic")
####################################################




#Soil Moisture Levels#
#501 = .06 m
#502 = .15 m
#503 = .26 m
#504 = .36 m
#505 = .46 m
#506 = .55 m

#has 5 plots


folder_list = c("NEON.D01.HARV.DP1.00094.001.2017-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2017-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2017-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2017-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2018-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2018-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2018-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2018-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2019-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2019-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2019-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2019-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2020-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2020-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2020-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2020-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2021-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2021-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2021-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2021-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2022-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00094.001.2022-07.basic.20221205T003615Z.PROVISIONAL",
                "NEON.D01.HARV.DP1.00094.001.2022-08.basic.20221204T232840Z.PROVISIONAL",
                "NEON.D01.HARV.DP1.00094.001.2022-09.basic.20221206T021548Z.PROVISIONAL",
                "NEON.D01.HARV.DP1.00094.001.2023-06.basic.20230706T002534Z.PROVISIONAL",
                "NEON.D01.HARV.DP1.00094.001.2023-07.basic.20230809T214706Z.PROVISIONAL",
                "NEON.D01.HARV.DP1.00094.001.2023-08.basic.20230905T014608Z.PROVISIONAL",
                "NEON.D01.HARV.DP1.00094.001.2023-09.basic.20231002T181910Z.PROVISIONAL")

year_month = c("2017-06",
               "2017-07",
               "2017-08",
               "2017-09",
               "2018-06",
               "2018-07",
               "2018-08",
               "2018-09",
               "2019-06",
               "2019-07",
               "2019-08",
               "2019-09",
               "2020-06",
               "2020-07",
               "2020-08",
               "2020-09",
               "2021-06",
               "2021-07",
               "2021-08",
               "2021-09",
               "2022-06",
               "2022-07",
               "2022-08",
               "2022-09",
               "2023-06",
               "2023-07",
               "2023-08",
               "2023-09"
              )

serial = c("20221204T191813Z","20221204T215858Z","20221204T215317Z","20221204T210145Z",
           "20221204T201244Z","20221204T225550Z","20221204T191621Z","20221204T204237Z",
            "20221204T214601Z","20221204T220225Z","20221204T202737Z","20221204T195808Z",
           "20221204T195909Z","20221204T210801Z","20221204T230357Z","20221204T231534Z",
           "20221204T200801Z","20221204T202434Z","20221204T184108Z","20221204T211221Z",
           "20221204T234321Z", "20221205T002358Z", "20221204T231617Z", "20221206T020434Z",
           "20230706T002534Z","20230809T214706Z","20230905T014608Z","20231002T181910Z")



plot_num = c("1")
plot_depth = c("2")

df_sm <- data.frame()

for (x in 1:length(folder_list)){
  
  wd = paste("/Users/jurado/Harvard_Forest/filesToStack00094",folder_list[x],sep = "/")
  setwd(wd)
  yr_mn = year_month[x]
  ser = serial[x]
    
    
  for (x in 1:length(plot_num)){
    plt_num = plot_num[x]
    for (x in 1:length(plot_depth)){
      file = paste("NEON.D01.HARV.DP1.00094.001.00",plt_num,sep ="") #plot_num
      file = paste(file,".50", sep ="")
      file = paste(file,plot_depth[x], sep ="")
      file = paste(file,".030.SWS_30_minute.", sep ="")
      file = paste(file,yr_mn, sep ="") #year and month
      file = paste(file,".basic.", sep ="") 
      file = paste(file,ser, sep ="") #serial
      file = paste(file,".csv", sep ="") #serial
      print(file)
      if (file.exists(file)){
        sm <- read.csv(file)
        sm$DEPTH <- as.numeric(plot_depth[x])
        df_sm <- rbind(df_sm,sm)
      }

    }
  }
}  


df_sm <- df_sm %>% filter(VSWCFinalQF == 0 )


#Select one median value from across all plots spatially at each depth
df_new <- df_sm %>%  group_by(endDateTime, DEPTH) %>%reframe(med_swc = median(VSWCMean, na.rm=TRUE),
                                                             max_swc = median(VSWCMaximum, na.rm=TRUE))
df_new$index <- seq(1,nrow(df_new),1)
#number of distinct half hour intervals in dataset
num_obv <-  df_sm %>%  group_by(endDateTime) %>%reframe(med_swc = median(VSWCMean, na.rm=TRUE),
                                                               sd_swc = sd(VSWCMean, na.rm=TRUE))
num_obv <- nrow(num_obv)

soil_col <- data.frame( DEPTH.interp = rep(seq(2,4,1), times=num_obv))

df_new <- df_new %>%
  complete(endDateTime, DEPTH)


#linearly interpolate soil water column
df_new <- df_new %>%
  mutate(med_swc = na.approx(med_swc))



#Sum of water content within soil column, median of the column depths
df_new <- df_new %>% group_by(endDateTime) %>% reframe(VSWC.mean = median(med_swc,na.rm = TRUE),
                                                       VSWC.max = median(max_swc, na.rm=TRUE))

df_new$datetime <- as.POSIXct(df_new$endDateTime, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC") 
df_new$datetime<- as.POSIXct(df_new$datetime, format = "%Y-%m-%dT%H:%M:%S", tz = "EST") 
plot(df_new$datetime,df_new$VSWC.col, type = "l")
plot(df_new$datetime,df_new$VSWC.mean, type = "l")







######FLUX DATA########


folder_list = c("NEON.D01.HARV.DP4.00200.001.2017-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2017-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2017-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2017-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2018-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2018-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2018-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2018-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2019-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2019-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2019-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2019-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2020-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2020-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2020-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2020-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2021-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2021-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2021-08.basic.20230127T120753Z.RELEASE-2023 2",
                "NEON.D01.HARV.DP4.00200.001.2021-09.basic.20230127T120753Z.RELEASE-2023 2",
                "NEON.D01.HARV.DP4.00200.001.2022-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP4.00200.001.2022-07.basic.20220810T224035Z.PROVISIONAL",
                "NEON.D01.HARV.DP4.00200.001.2022-08.basic.20230225T225559Z.PROVISIONAL",
                "NEON.D01.HARV.DP4.00200.001.2022-09.basic.20230225T225624Z.PROVISIONAL",
                "NEON.D01.HARV.DP4.00200.001.2023-06.basic.20230713T225900Z.PROVISIONAL",
                "NEON.D01.HARV.DP4.00200.001.2023-07.basic.20230810T213732Z.PROVISIONAL",
                "NEON.D01.HARV.DP4.00200.001.2023-08.basic.20230907T195057Z.PROVISIONAL",
                "NEON.D01.HARV.DP4.00200.001.2023-09.basic.20231013T205014Z.PROVISIONAL")

file_list = c("NEON.D01.HARV.DP4.00200.001.nsae.2017-06.basic.20221215T022303Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2017-07.basic.20221215T023843Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2017-08.basic.20221215T025413Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2017-09.basic.20221215T030921Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2018-06.basic.20221215T090448Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2018-07.basic.20221215T092218Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2018-08.basic.20221215T093710Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2018-09.basic.20221215T095218Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2019-06.basic.20221215T121219Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2019-07.basic.20221215T122934Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2019-08.basic.20221215T124633Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2019-09.basic.20221215T130324Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2020-06.basic.20221215T153201Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2020-07.basic.20221215T154738Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2020-08.basic.20221215T160708Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2020-09.basic.20221215T162259Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2021-06.basic.20221216T005046Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2021-07.basic.20221216T010713Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2021-08.basic.20221216T012510Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2021-09.basic.20221216T014228Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2022-06.basic.20221216T041952Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2022-07.basic.20220809T001307Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2022-08.basic.20230224T000842Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2022-09.basic.20230224T002819Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2023-06.basic.20230712T165914Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2023-07.basic.20230809T174858Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2023-08.basic.20230906T203557Z.h5",
              "NEON.D01.HARV.DP4.00200.001.nsae.2023-09.basic.20231011T203128Z.h5")


df <- data.frame()


#DATA SCRAPING LOOP#
for (x in 1:length(folder_list)){
  wd = paste("/Users/jurado/Harvard_Forest/filesToStack00200",folder_list[x],sep = "/")
  setwd(wd)
  h5ls(file_list[x])
  
  #Reading In Data
  
  
  LE <- h5read(file_list[x] , "/HARV/dp04/data/fluxH2o/nsae")
  LE_flag <- h5read(file_list[x] , "/HARV/dp04/qfqm/fluxH2o/nsae")
  H  <- h5read(file_list[x] , "/HARV/dp04/data/fluxTemp/nsae")
  H_flag  <- h5read(file_list[x] , "/HARV/dp04/qfqm/fluxTemp/nsae")
  T1 <- h5read(file_list[x] , "/HARV/dp01/data/tempAirLvl/000_040_30m/temp")
  T2 <- h5read(file_list[x] , "/HARV/dp01/data/tempAirTop/000_060_30m/temp")
  Q1 <- h5read(file_list[x] , "/HARV/dp01/data/h2oStor/000_040_30m/rtioMoleDryH2o")
  Q2 <- h5read(file_list[x] , "/HARV/dp01/data/h2oStor/000_060_30m/rtioMoleDryH2o")
  C1 <- h5read(file_list[x] , "/HARV/dp01/data/co2Stor/000_040_30m/rtioMoleDryCo2")
  C2 <- h5read(file_list[x] , "/HARV/dp01/data/co2Stor/000_060_30m/rtioMoleDryCo2")
  CO2 <- h5read(file_list[x] , "/HARV/dp04/data/fluxTemp/nsae")
  P <- h5read(file_list[x] , "/HARV/dp01/data/presBaro/000_025_30m/presAtm")
  Lwin <-h5read(file_list[x] , "/HARV/dp01/data/radiNet/000_060_30m/radiLwIn")
  Lwin_flag <-h5read(file_list[x] , "/HARV/dp01/qfqm/radiNet/000_060_30m/radiLwIn")
  Lwout_flag <-h5read(file_list[x] , "/HARV/dp01/qfqm/radiNet/000_060_30m/radiLwOut")
  Lwout <-h5read(file_list[x] , "/HARV/dp01/data/radiNet/000_060_30m/radiLwOut")
  Swin <-h5read(file_list[x] , "/HARV/dp01/data/radiNet/000_060_30m/radiSwIn")
  Swin_flag <-h5read(file_list[x] , "/HARV/dp01/qfqm/radiNet/000_060_30m/radiSwIn")
  Swout <-h5read(file_list[x] , "/HARV/dp01/data/radiNet/000_060_30m/radiSwOut")
  Swout_flag <-h5read(file_list[x] , "/HARV/dp01/qfqm/radiNet/000_060_30m/radiSwOut")
  G1 <- h5read(file_list[x] , "/HARV/dp01/data/fluxHeatSoil/001_501_30m/fluxHeatSoil")
  G1_flag <- h5read(file_list[x] , "/HARV/dp01/qfqm/fluxHeatSoil/001_501_30m/fluxHeatSoil")
  G2 <- h5read(file_list[x] , "/HARV/dp01/data/fluxHeatSoil/003_501_30m/fluxHeatSoil")
  G2_flag <- h5read(file_list[x] , "/HARV/dp01/qfqm/fluxHeatSoil/003_501_30m/fluxHeatSoil")
  G3 <- h5read(file_list[x] , "/HARV/dp01/data/fluxHeatSoil/005_501_30m/fluxHeatSoil")
  G3_flag <- h5read(file_list[x] , "/HARV/dp01/qfqm/fluxHeatSoil/005_501_30m/fluxHeatSoil")
  
  #quality filter
  Lwin$qfqm <- Lwin_flag$qfFinl
  #Lwin$mean <- ifelse(Lwin$qfqm >0,NA,Lwin$mean)
  Lwout$qfqm <- Lwout_flag$qfFinl
  #Lwout$mean <- ifelse(Lwout$qfqm >0,NA,Lwout$mean)
  Swin$qfqm <- Swin_flag$qfFinl
  #Swin$mean <- ifelse(Swin$qfqm >0,NA,Swin$mean)
  Swout$qfqm <- Swout_flag$qfFinl
  #Swout$mean <- ifelse(Swout$qfqm >0,NA,Swout$mean)
  H$qfqm <- H_flag$qfFinl
  #H$flux <- ifelse(H$qfqm > 0,NA,H$flux)
  LE$qfqm <- LE_flag$qfFinl
  #LE$flux <- ifelse(LE$qfqm >0,NA,LE$flux)
  G1$qfqm <- G1_flag$qfFinl
  G2$qfqm <- G2_flag$qfFinl
  G3$qfqm <- G3_flag$qfFinl
  #G1$mean <- ifelse(G1$qfqm >0,NA,G1$mean)
  #G2$mean <- ifelse(G1$qfqm >0,NA,G2$mean)
  #G3$mean <- ifelse(G1$qfqm >0,NA,G3$mean)
  G <- data.frame(G1$mean,G2$mean,G3$mean)
  G$G_mean <- rowMeans(G[,c(1,2,3)], na.rm=TRUE)
  
  Lwin$timeEnd <- as.POSIXct(Lwin$timeEnd, format = "%Y-%m-%dT%H:%M:%S", tz = "EST")

  
  EF <- data.frame(Lwin$timeEnd,Lwin$mean,Lwout$mean,Swin$mean,Swout$mean,LE$flux,H$flux,G$G_mean)
  EF$Rnet <- Lwin$mean-Lwout$mean+Swin$mean-Swout$mean

  colnames(EF)[1] ="datetime"
  colnames(EF)[2] ="Lwin"
  colnames(EF)[3] ="Lwout"
  colnames(EF)[4] ="Swin"
  colnames(EF)[5] = "Swout"
  colnames(EF)[6] = "LE"
  colnames(EF)[7] = "H"
  colnames(EF)[8] ="G"
  colnames(EF)[9] ="Rnet"
  
  
  
  EF$LE.bud <- EF$Rnet-EF$H-EF$G
  EF$P.mean <- P$mean
  
  rm(LE,Lwin,Lwout,Swin,Swout)
  
  df <- rbind(df,EF)
  
  print(wd)
  
}

#not doing EST for some reason, did HST for SWIN correction
df$datetime<- as.POSIXct(df$datetime, format = "%Y-%m-%dT%H:%M:%S", tz = "HST") 
####COMBINE FLUX AND SOIL MOISTURE DATA####

df_sm_flux <- merge(df,df_new,by = "datetime")

####DERIVED VALUES

#Revisit assumptions staring here


df_sm_flux$EF <- df_sm_flux$LE/df_sm_flux$Rnet
df_sm_flux$EF.bud <- df_sm_flux$LE.bud/df_sm_flux$Rnet

#LE BUDGET INFILL#

df_sm_flux$LE.combined <- coalesce(df_sm_flux$LE,df_sm_flux$LE.bud)
df_sm_flux$EF.bud <- df_sm_flux$LE.combined/df_sm_flux$Rnet




#Despike 



df_sm_flux$Rnet <- replace(df_sm_flux$Rnet, abs(df_sm_flux$Rnet) > 1000, NA)
df_sm_flux$LE <- replace(df_sm_flux$LE, abs(df_sm_flux$LE) > 800, NA)
df_sm_flux$LE.combined <- replace(df_sm_flux$LE.combined, abs(df_sm_flux$LE.combined) > 800, NA)
df_sm_flux$budget <- (df_sm_flux$Rnet-df_sm_flux$G)/(df_sm_flux$H+df_sm_flux$LE.combined)
#df_sm_flux$LE <- replace(df_sm_flux$LE, df_sm_flux$LE < 0, NA)
#df_sm_flux$Rnet <- replace(df_sm_flux$Rnet, df_sm_flux$Rnet< 0, NA)
#df_sm_flux$H <- replace(df_sm_flux$H, df_sm_flux$H< 0, NA)
#df_sm_flux$H <- replace(df_sm_flux$H, abs(df_sm_flux$H)> 800, NA)


#BY HOUR
library(stringr)
split <- df_sm_flux$datetime %>% str_split_fixed(" ",2)  
split <- data.frame(split)
df_sm_flux$hour <- split$X2
df_sm_flux$date <- split$X1
split <- df_sm_flux$datetime %>% str_split_fixed("-",3)  
split <- data.frame(split)
df_sm_flux$year <- split$X1
df_sm_flux$month <- split$X2
df_sm_flux$month_year <- paste(df_sm_flux$year,df_sm_flux$month, sep="-")


#FILTER

#Filter for all rows with an EF > 1 and <0
df_sm_flux <- df_sm_flux %>% filter(EF.bud <= 1 & EF.bud >= 0)
df_sm_flux <- df_sm_flux %>% filter(budget <= 2 & budget >= -2)
#quantify percent of data lost





#Filter for all rows during the day 
df_sm_flux$hour <- as.POSIXct(df_sm_flux$hour, format = "%H:%M:%S")
df_sm_flux <- df_sm_flux %>% filter(hour >= as.POSIXct('2023-12-22 5:30:00') & hour <= as.POSIXct('2023-12-22 17:30:00'))


#Create Seasonal Averages#

dfYEARLY<- df_sm_flux %>% group_by(year) %>% summarise(EF.avg = mean(EF.bud, na.rm = TRUE),
                                                        EF.sd = sd(EF.bud, na.rm = TRUE),
                                                       swc.avg = mean(VSWC.mean, na.rm =TRUE),
                                                       swc.sd = sd(VSWC.mean, na.rm = TRUE))

dfmonthly<- df_sm_flux %>% group_by(month_year) %>% summarise(EF.avg = mean(EF.bud, na.rm = TRUE),
                                                       EF.sd = sd(EF.bud, na.rm = TRUE),
                                                       swc.avg = mean(VSWC.mean, na.rm =TRUE),
                                                       swc.sd = sd(VSWC.mean, na.rm = TRUE))




#merge with df_sm_flux on datetime


df_sm_flux <- merge(df_sm_flux,df_RH,by = "datetime", all.x = TRUE)



#Averaging into bins
EF <- df_sm_flux %>% mutate(swc_binned = cut(VSWC.mean, breaks=75))

binned <- EF %>% group_by(swc_binned) %>% summarise(EF.avg = mean(EF.bud, na.rm = TRUE),
                                                    swc.avg = mean(VSWC.mean, na.rm=TRUE),
                                                    VPD.avg = mean(VPD, na.rm=TRUE),
                                                    LE.avg = mean(LE.combined, na.rm=TRUE))

df_table <- data.frame(table(EF$swc_binned))
colnames(df_table)[1] ="swc_binned"
binned <- merge(binned,df_table,by = "swc_binned")

binned <- binned %>% filter(Freq > 10)




#SL_T REGIME
plot(binned$swc.avg,binned$EF.avg, pch = 0,ylab = "Evaporative Fraction", xlab = "Soil Water Content", type ="p")

y = weightedLowess(binned$swc.avg, binned$EF.avg, weights = binned$Freq,
               delta=NULL, npts = 200, span = 0.20, iterations = 8)

lines(binned$swc.avg,y$fitted, col="black",lty = 2,lwd =2)

par(new = TRUE)                             # Add new plot
plot(binned$swc.avg, binned$VPD.avg, type = "p",pch =3 , col = "red",              # Create second plot without axes
     axes = FALSE, xlab = "", ylab = "",ylim = c(0,1.2))




y = weightedLowess(binned$swc.avg, binned$VPD.avg, weights = binned$Freq,
                   delta=NULL, npts = 200, span = 0.05, iterations = 8)
lines(binned$swc.avg,y$fitted, col="red",lty = 2,lwd =2)
axis(side = 4, at = pretty(seq(0, 1.2, .1))) 
mtext("Vapor Pressure Deficit [KPa]", side = 4, line = 3) 
title("Evapotranspiration Regimes")
subtitle = "HARV June-September 2017-2023"
mtext(subtitle)
rect(.21 , -.05, .38, 1.25, col=rgb(0,1,0,alpha=0.075), lty = 3 )
rect(.095,-.05, .21, 1.25, col=rgb(1,.647,0,alpha=0.075), lty = 3 )
rect(-.01, -.05, .095, 1.25, col=rgb(1,0,0,alpha=0.075), lty = 3 )
legend(.10,.45,legend=c("Evaporative Fraction","Vapor Pressure Deficit"), col=c("black","red"),
       pch=c(0,3), ncol=1,cex=.75)

text(.045, 0, substitute(paste(bold('WILT'))), col ="red")
text(.157, .06, substitute(paste(bold('FIELD'))), col ="darkorange")
text(.157, 0, substitute(paste(bold('CAPACITY'))), col ="darkorange")
text(.3, 0, substitute(paste(bold('TRANSITIONAL'))), col ="darkgreen")
text(.1, 1.15, substitute(paste(bold('DRY REGIME'))), col ="brown")
text(.33, -.3, "Soil Depth: .15 - .36m", col ="black", xpd=NA,cex=.85)
#Regime SWC

#DRY: 0 - .095
#BUFFER: .095 - .21
#TRANSITIONAL: .21 - .367
#WET: .285 + 

#SL_T REGIME GGPLOT

library(ggplot2)

#create scatterplot with point size based on value of qsec
ggplot(data=binned, aes(x=swc.avg, y=EF.avg,weight=Freq)) +
  geom_point(aes(size=Freq))+
  geom_smooth(method = "loess")+
  ggtitle("Evapotranspiration Regimes")+
  labs(y= "Evaporative Fraction", x = "Soil Water Content",
       subtitle = "HARV June-September 2017-2022",
       caption = "Soil Depths .05 - .55m")+
  theme_bw()

############################PRECIP AND EVAPORATION##############################

setwd("/Users/jurado/Harvard_Forest")

zipsByProduct(dpID="DP1.00006.001", site="HARV",
              startdate="2017-06", enddate="2023-09", package="basic", check.size=FALSE)


folder_list = c("NEON.D01.HARV.DP1.00006.001.2017-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00006.001.2017-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00006.001.2017-08.basic.20230127T120753Z.RELEASE-2023 2",
                "NEON.D01.HARV.DP1.00006.001.2017-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00006.001.2018-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00006.001.2018-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00006.001.2018-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00006.001.2018-09.basic.20230127T120753Z.RELEASE-2023 2",
                "NEON.D01.HARV.DP1.00006.001.2019-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00006.001.2019-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00006.001.2019-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00006.001.2019-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00006.001.2020-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00006.001.2020-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00006.001.2020-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00006.001.2020-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00006.001.2021-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00006.001.2021-07.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00006.001.2021-08.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00006.001.2021-09.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00006.001.2022-06.basic.20230127T120753Z.RELEASE-2023",
                "NEON.D01.HARV.DP1.00006.001.2022-07.basic.20230320T183321Z.PROVISIONAL",
                "NEON.D01.HARV.DP1.00006.001.2022-08.basic.20230320T185748Z.PROVISIONAL",
                "NEON.D01.HARV.DP1.00006.001.2022-09.basic.20230320T184852Z.PROVISIONAL",
                "NEON.D01.HARV.DP1.00006.001.2023-06.basic.20230706T000512Z.PROVISIONAL",
                "NEON.D01.HARV.DP1.00006.001.2023-07.basic.20230809T233828Z.PROVISIONAL",
                "NEON.D01.HARV.DP1.00006.001.2023-08.basic.20230906T015816Z.PROVISIONAL",
                "NEON.D01.HARV.DP1.00006.001.2023-09.basic.20231002T191037Z.PROVISIONAL")

file_list = c("NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2017-06.basic.20221204T183906Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2017-07.basic.20221204T215843Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2017-08.basic.20221204T200331Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2017-09.basic.20221204T220300Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2018-06.basic.20221204T233324Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2018-07.basic.20221204T231654Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2018-08.basic.20221204T221134Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2018-09.basic.20221204T200052Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2019-06.basic.20221204T191356Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2019-07.basic.20221204T222335Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2019-08.basic.20221204T214639Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2019-09.basic.20221204T211112Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2020-06.basic.20221204T214449Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2020-07.basic.20221204T234726Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2020-08.basic.20221204T224318Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2020-09.basic.20221204T193058Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2021-06.basic.20221204T192957Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2021-07.basic.20221204T204239Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2021-08.basic.20221204T234429Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2021-09.basic.20221204T232645Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2022-06.basic.20221204T200519Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2022-07.basic.20230320T183321Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2022-08.basic.20230320T185748Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2022-09.basic.20230320T184852Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2023-06.basic.20230706T000512Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2023-07.basic.20230809T233828Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2023-08.basic.20230906T015816Z.csv",
              "NEON.D01.HARV.DP1.00006.001.900.000.030.PRIPRE_30min.2023-09.basic.20231002T191037Z.csv")

df_precip <- data.frame()


#DATA SCRAPING LOOP#
for (x in 1:length(folder_list)){
  wd = paste("/Users/jurado/Harvard_Forest/filesToStack00006/",folder_list[x],sep ="") 
  setwd(wd)
  file = file_list[x]
  precip <- read.csv(file)
  df_precip <- rbind(df_precip,precip)
  
}

#in Celsius


df_precip$datetime <- as.POSIXct(df_precip$endDateTime, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC") 
df_precip$datetime<- as.POSIXct(df_precip$datetime, format = "%Y-%m-%dT%H:%M:%S", tz = "EST") 
plot(df_precip$datetime,df_precip$priPrecipBulk,type="l")

df_sm_flux_precip <- merge(df_sm_flux,df_precip,by = "datetime")



####Fisher Tower Precip###
df_precip <- read_csv("~/Harvard_Forest/hf001-10-15min-m.csv")



split <- df_precip$datetime %>% str_split_fixed("-",3)
split <- data.frame(split)

df_precip$year <- split$X1
df_precip$month <- as.numeric(split$X2)
df_precip$day <- split$X3



df_precip <- df_precip %>% filter(6 < month < 10)
#rip into years, convert LE, convert r
plot(df_)



####Evaporation Deficit Plot####

#LATENT HEAT DATA
setwd("/Users/jurado/Harvard_Forest")
library(readr)
hf004 <- read_csv("hf004-01-final.csv")



df_evap <- data.frame(hf004$f_h20) #in millimolePerMeterSquaredPerSecond
df_evap$datetime <- hf004$datetime
df_evap$kg_h2o <- (df_evap$hf004.f_h20/1000)*(18.01528/1000)


split <- df_evap$datetime %>% str_split_fixed(" ",2)
split <- data.frame(split)

########MAKE INTO DAYTIME ONLY VALUES
#rm(EF$datetime)
df_evap$date <- split$X1
df_evap$time <- split$X2

split <- df_evap$time %>% str_split_fixed(":",3)
split <- data.frame(split)
df_evap$hour <- split$X1

split <- df_evap$date %>% str_split_fixed("-",3)
split <- data.frame(split)
df_evap$year <- split$X1
df_evap$month <- split$X2




df_evap$year<- as.numeric(df_evap$year)
df_evap$month <- as.numeric(df_evap$month)
df_evap$hour <- as.numeric(df_evap$hour)
df_evap <- filter(df_evap, month <10)
df_evap <- filter(df_evap, month >5)

#Daily filter determined from daily mean across all time periods
df_evap <- filter(df_evap, hour <20)
df_evap <- filter(df_evap, hour  >5)

df_evap$date <- paste(df_evap$year,df_evap$month , sep="-")

df_evap <- na.omit(df_evap)
#sum into months, remove months that have >90% missing data
df_evap_daily <- df_evap %>% group_by(hour) %>% summarise(
                                                          mean_h2o_kg = mean(kg_h2o,na.rm=TRUE),
                                                         )
plot(df_evap_daily$hour,df_evap_daily$mean_h2o_kg)



df_evap_MONTHLY <- df_evap %>% group_by(date) %>% summarise(sum_h2o_kg = sum(kg_h2o,na.rm=TRUE),
                                                            num = n(),
                                                            mean_h2o_kg = mean(kg_h2o,na.rm=TRUE),
                                                            median_h2o_kg= median(kg_h2o,na.rm=TRUE),
                                                            std_err_h2o_kg = sd(kg_h2o, na.rm=TRUE)/sqrt(length((kg_h2o))))
df_evap_MONTHLY$date <- paste(df_evap_MONTHLY$date,"01",sep= "-")
df_evap_MONTHLY$date <- as.Date(df_evap_MONTHLY$date, "%Y-%m-%d")


df_evap_MONTHLY$data_perc <- df_evap_MONTHLY$num/732



df_evap_MONTHLY <- filter(df_evap_MONTHLY, data_perc> .1)


plot(df_evap_MONTHLY$date ,df_evap_MONTHLY$mean_h2o_kg, type="b")

arrows(x0=seq(1,119,1),y0 = df_evap_MONTHLY$mean_h2o_kg - df_evap_MONTHLY$std_err_h2o_kg, x1=seq(1,119,1),
       y1=df_evap_MONTHLY$mean_h2o_kg + df_evap_MONTHLY$std_err_h2o_kg,
       code=3, angle=90, length=0.05, lwd = 1.5, col="darkgrey")

abline(lm(df_evap_MONTHLY$mean_h2o_kg~seq(1,119,1)))
summary(lm(df_evap_MONTHLY$mean_h2o_kg~seq(1,119,1)))




#PRECIP DATA
hf300_03 <- read_csv("hf300-03-monthly-m.csv")

split <- hf300_03$date %>% str_split_fixed("-",2)
split <- data.frame(split)

#Split by month, get rid of years not needed.
hf300_03$year <- split$X1
hf300_03$month <- split$X2
hf300_03$year<- as.numeric(hf300_03$year)
hf300_03$month <- as.numeric(hf300_03$month)

hf300_03 <- filter(hf300_03, month <10)
hf300_03 <- filter(hf300_03, month >5)
hf300_03 <- filter(hf300_03, year >1992)
 
hf300_03$date <- paste(hf300_03$date,"01",sep= "-")
hf300_03$date <- as.Date(hf300_03$date, "%Y-%m-%d")

barplot(hf300_03$prec,hf300_03$date,)


####COMBINED PLOT###
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(df_evap_MONTHLY$date ,df_evap_MONTHLY$mean_h2o_kg*1000*1000, type="b", ylab="Average Water Flux [kg/s km]",
     ylim = c(0,90),xlab = "Date")
par(new=TRUE)
barplot(hf300_03$prec, ylim =c(0,500), ylab="", axes ="FALSE", border=FALSE)
axis(side = 4, at = pretty(seq(0, 300, 10))) 
par(new=TRUE)
plot(df_evap_MONTHLY$date ,df_evap_MONTHLY$mean_h2o_kg*1000*1000, type="b", ylab="",xlab="", lwd=2,
     ylim = c(0,90), col = "cadetblue")

mtext("Precipitation [mm]", side = 4, line = 3) 
arrows(x0=df_evap_MONTHLY$date,y0 = df_evap_MONTHLY$mean_h2o_kg*1000*1000 - df_evap_MONTHLY$std_err_h2o_kg*1000*1000, x1=df_evap_MONTHLY$date,
       y1=df_evap_MONTHLY$mean_h2o_kg*1000*1000 + df_evap_MONTHLY$std_err_h2o_kg*1000*1000,
       code=3, angle=90, length=0.05, lwd = 1.5, col = "cadetblue")
title("Latent Heat and Precipitation")
subtitle ="EMS Tower June - September 1992 - 2022"
mtext(subtitle)






df_RH <- data.frame(VPSat <- (610.7*10**(7.5*hf004$ta_27_9m/(237.3+hf004$ta_27_9m)))/1000)

df_RH$VPair <- (610.7*10**(7.5*hf004$ta_27_9m/(237.3+hf004$ta_27_9m)))/1000 * (hf004$rh_27_9m/100)

df_RH$VPD <- df_RH$VPSat - df_RH$VPair

df_RH$datetime <- hf004$datetime




plot(df_RH$datetime,df_RH$VPD)
abline(lm(df_RH$VPD~df_RH$datetime),col="red")
summary(lm(df_RH$VPD~df_RH$datetime))



#Gap Filled Evapotranspiration deficit 


library(readr)
setwd("/Users/jurado/")
Ha1_ETfill_9222 <- read_csv("Downloads/Ha1_LEfill_9222.csv")

hf300_03 <- read_csv("Harvard_Forest/hf300-03-monthly-m.csv")



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

Ha1_ETfill_9222$ET_mm <- Ha1_ETfill_9222$LE_f*(1*10**-6)*(1/2.26)*1000*(1/10000)*3600*10


ET_PRECIP<- Ha1_ETfill_9222 %>% group_by(year) %>% summarise(ET_sum = sum(ET_mm, na.rm =TRUE))



ET_PRECIP$date <- as.character(ET_PRECIP$date)

hf300_03$date <- as.character(hf300_03$date)


hf300_03$year <- as.numeric(substr(as.character(hf300_03$date), start = 1, stop = 4))
hf300_03$month <- as.numeric(substr(as.character(hf300_03$date), start = 6, stop = 7))
hf300_03<- hf300_03 %>% filter(hf300_03$month > 5)
hf300_03 <- hf300_03 %>% filter(hf300_03$month <10)
hf300_03_PRECIP<- hf300_03 %>% group_by(year) %>% summarise(prec_sum = sum(prec, na.rm =TRUE))


setwd("/Users/jurado/Harvard_Forest")
df_prcp <- read.csv("hf001-06-daily-m.csv", na = "-9999", header=TRUE) #EMS data
df_prcp$year <- substr(as.character(df_prcp$date), start = 1, stop = 4)
df_prcp$year <- as.numeric(df_prcp$year)
df_prcp$month <- substr(as.character(df_prcp$date), start = 6, stop = 7)
df_prcp$month <- as.numeric(df_prcp$month)
df_prcp$date <- paste(as.character(df_prcp$year),as.character(df_prcp$month), sep = "-0")



df_prcp <- df_prcp %>% filter(df_prcp$month > 5)

df_prcp <- df_prcp %>% filter(df_prcp$month <10)

precip_EMS <- df_prcp  %>% group_by(year) %>% summarise(prec_sum = sum(prec, na.rm =TRUE))



ET_PRECIP_1 <- merge(hf300_03_PRECIP,ET_PRECIP, by = "year")
ET_PRECIP_1$diff <- ET_PRECIP_1$prec-ET_PRECIP_1$ET_sum
ET_PRECIP_1$col <- ifelse(ET_PRECIP_1$diff>=0,"cadetblue4","red3" )
ET_PRECIP_1$year <- as.numeric(substr(as.character(ET_PRECIP_1$date), start = 1, stop = 4))



ET_PRECIP_2 <- merge(precip_EMS,df_NEON_ET , by = "year")
ET_PRECIP_2$diff <- ET_PRECIP_2$prec_sum-ET_PRECIP_2$ET_sum
ET_PRECIP_2$col <- ifelse(ET_PRECIP_2$diff>=0,"cadetblue4","red2" )
ET_PRECIP_2$year <- as.numeric(substr(as.character(ET_PRECIP_2$date), start = 1, stop = 4))

#ET is in mm
#prec is in mm



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




 #########################Mean Variables CLASS Model#############################

#' Lat:
#' Lon:
#' DOY
#' 
#' DIURNAL STARTING POINT
#' Initial ABL height
#' Initial Potential Temp
#' Initial Specific Humidity 
#' 
#' FOR DAYTIME ONLY
#' Average velocity friction
#' Average wind
#' Average air Temp for surface temp
#' Average Roughness length for momentum -
#' Average Roughness Length for Scalars -
#' 
#' ALL TIME
#' Average soil temp layers 1 and 2 
#' Average soil moisture layers 1 and 2
#' Average surface albedo
#' Surface pressure
#' 
#' 
#' Field capacity soil moisture 3 days after heavy rain event
#' Wilting point ?
#' Saturated ?
#' 
#' LAI 3.7




























