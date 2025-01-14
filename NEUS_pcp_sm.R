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
library(ggplot2) #use ggplot2 to add layer for visualizatio
library(readr)
library("anytime")   



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

#this is for checking how rain event frequency has changed over time 
#df_prcp["prec_mm"][df_prcp["prec_mm"] > 20] <- NA

#Which months are we looking at
#df_prcp["prec_mm"][df_prcp["month"] <5 | df_prcp["month"] >9] <- NA

#for loop to get and store the top 1% per year.

perc_hvy <- c()
num_hvy <- c()
sum_hvy <- c()
prec_number <- c()
Q <- c()



for(x in 1964:2023){
  x <- toString(x)
  print(x)
  df_prcp <- df_prcp[df_prcp$date >= paste(x,"01-01",sep = "-") &   df_prcp$date <= paste(x,"12-31",sep = "-"), ]
  q <- quantile(df_prcp$prec_mm, prob=c(.90), type=1, na.rm =TRUE)
  Q <- append(Q,q)
  
  prec_number <- append(prec_number,length(na.omit(df_prcp$prec_mm)))
  
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
  
  #df_prcp["prec_mm"][df_prcp["prec_mm"] > 20] <- NA
  
  #df_prcp["prec_mm"][df_prcp["month"] >5 & df_prcp["month"] <9] <- NA
  
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


#####NUMBER OF RAIN EVENTS####

plot(prec_number)



#####SINCE 2000
df_perc_hvy_2000 <- df_perc_hvy %>% filter(df_perc_hvy$Year > 2000)
abline(lm(df_perc_hvy_2000$perc_hvy ~ df_perc_hvy_2000$Year), lwd =2, lty =2, col = "red")
summary(lm(df_perc_hvy_2000$perc_hvy ~ df_perc_hvy_2000$Year))
Kendall(df_perc_hvy_2000$Year,df_perc_hvy_2000$perc_hvy)


############Top 10% over entire time period########


df_prcp <- df_prcp<- read_csv("Downloads/dailyrain_19642023.csv")
df_prcp["prec_mm"][df_prcp["prec_mm"] == 0] <- NA
df_prcp$date <- as.Date(df_prcp$date, "%m/%d/%Y")
df_prcp$date <- update(df_prcp$date, year = df_prcp$year)


perc_hvy <- c()
sum_hvy <- c()

#df_prcp <- df_prcp %>% filter( month > 4 &  month <10)
q <- quantile(df_prcp$prec_mm, prob=c(.90), type=1, na.rm =TRUE)

x=1964

for(x in 1964:2023){
  x <- toString(x)
  print(x)
  df_prcp <- df_prcp[df_prcp$date >= as.Date(paste(x,"01-01", sep ="-")) & df_prcp$date <= as.Date(paste(x,"12-31",sep = "-")), ]
  #df_prcp <- df_prcp %>% filter( month > 4 &  month <10)
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

# Create a new dataframe with year and total precipitation
yearly_precip <- df_prcp %>%
  mutate(year = format(as.Date(date), "%Y")) %>%  # Extract year from the date
  group_by(year) %>%                              # Group by year
  summarise(total_precip = sum(prec_mm, na.rm = TRUE))  # Calculate total precipitation




df_perc_hvy <- data.frame(perc_hvy)
df_perc_hvy$Year <- seq(1964,2023,1)
df_perc_hvy$Tot <- sum_hvy


# Adjust plot margins to ensure visibility of the right-hand y-axis title
par(mar = c(5, 4, 4, 5) + 0.5)

# Base plot: Percent of precipitation contributed by heavy storms
plot(df_perc_hvy$Year, df_perc_hvy$perc_hvy, 
     type = "l", 
     xlab = "Year", 
     ylab = "Precip. Contributed by Heavy Storms [%]", 
     lwd = 2, 
     ylim = c(18, 100), 
     pch = 16)

# Add linear trend line for percent contribution
abline(lm(df_perc_hvy$perc_hvy ~ df_perc_hvy$Year), 
       lty = 1, lwd = 1, col = alpha("black", 0.6))

# Add title and subtitle
title(expression(paste(bold("Yearly Contributions of Heavy Precipitation"))))
subtitle = "EMS, 1964 - 2023"
mtext(subtitle)

# Add new plot (total contribution) on the same figure
par(new = TRUE)
plot(df_perc_hvy$Year, df_perc_hvy$Tot, 
     type = "b", 
     lwd = 2, 
     pch = 15, 
     ylim = c(0, 2000), 
     ylab = "", 
     xlab = "", 
     axes = FALSE, 
     col = alpha("cadetblue4", 0.8))

# Add right-hand y-axis
axis(side = 4)
mtext("Precip. Totals [mm]", side = 4, line = 3)

# Add linear trend line for total contribution
abline(lm(df_perc_hvy$Tot ~ df_perc_hvy$Year), 
       lty = 2, lwd = 1, col = alpha("cadetblue4", 0.8))

# Overlay yearly total precipitation data
lines(yearly_precip$year, yearly_precip$total_precip, 
      type = "b", 
      col = "red", 
      lwd = 2, 
      pch = 17)

# Add a legend for all series
legend(1965, 1100, 
       legend = c("Percent", "Heavy Precip.", "Total Precip"), 
       lty = c(1, 1, 1), 
       lwd = c(2, 2, 2), 
       col = c("black", alpha("cadetblue4", 0.8), "red"), 
       pch = c(NA, 15, 17), 
       bty = "n")

########################################


yearly_precip$year <- as.numeric(yearly_precip$year)


# Adjust plot margins to ensure visibility of the right-hand y-axis title
par(mar = c(5, 4, 4, 5))

# Base plot: Total precipitation contributed by heavy storms
plot(df_perc_hvy$Year, df_perc_hvy$Tot, 
     type = "b", 
     lwd = 2, 
     pch = 15, 
     ylim = c(0, 2000), 
     xlab = "Year", 
     ylab = "Precipitation [mm]",
     col = alpha("cadetblue4", 0.8))

# Add linear trend line for total contribution by heavy storms
abline(lm(df_perc_hvy$Tot ~ df_perc_hvy$Year), 
       lty = 2, lwd = 1, col = alpha("cadetblue4", 0.8))



# Add linear trend line for yearly total precipitation
abline(lm(yearly_precip$total_precip~yearly_precip$year), 
       lty = 2, lwd = 1, col = "black")

# Add yearly total precipitation data (from yearly_precip dataframe)
lines(yearly_precip$year, yearly_precip$total_precip, 
      type = "b", 
      col = "black", 
      lwd = 2, 
      pch = 17)



# Add a legend for both series
legend("topleft", 
       legend = c("Heavy Precip. Contribution", "Yearly Total Precipitation"), 
       lty = c(2, 1), 
       lwd = c(2, 2), 
       col = c(alpha("cadetblue4", 0.8), "black"), 
       pch = c(15, 17), 
       bty = "n")


# Add title and subtitle
title(expression(paste(bold("Total Precip. and Heavy Storm Contributions"))))
subtitle = "EMS, 1964 - 2023"
mtext(subtitle)


# Statistical summaries
summary(lm(df_perc_hvy$perc_hvy ~ df_perc_hvy$Year))
Kendall(df_perc_hvy$Year, df_perc_hvy$perc_hvy)

summary(lm(df_perc_hvy$Tot ~ df_perc_hvy$Year))
Kendall(df_perc_hvy$Year, df_perc_hvy$Tot)

summary(lm(yearly_precip$total_precip~yearly_precip$year))
Kendall(yearly_precip$year, yearly_precip$total_precip)




mean(df_perc_hvy$Tot[51:60])
mean(df_perc_hvy$Tot[1:10])



