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


# Initialize an empty list to store dataframes
df_list <- list()
plot_depth = c("1")
# Loop through plot numbers (1 through 5)
for (plot in 1:5) {
  
  # Initialize the dataframe for each plot (df1 through df5)
  df_name <- paste0("df", plot)  # Create dynamic dataframe names like df1, df2, etc.
  
  df_sm <- data.frame()  # Initialize the dataframe for this plot
  
  # Loop through folder_list, year_month, serial and plot_depth as before
  for (x in 1:length(folder_list)) {
    
    # Set working directory for each folder
    wd = paste("/Users/jurado/Harvard_Forest/filesToStack00094", folder_list[x], sep = "/")
    setwd(wd)
    
    yr_mn = year_month[x]
    ser = serial[x]
    
    for (plt_num in plot) {  # Just use the current plot number, no need for a nested loop
      for (depth in plot_depth) {
        
        # Construct file path
        file <- paste0("NEON.D01.HARV.DP1.00094.001.00", plt_num, ".50", depth, ".030.SWS_30_minute.", yr_mn, ".basic.", ser, ".csv")
        
        # Check if file exists and read it
        if (file.exists(file)) {
          sm <- read.csv(file)
          sm$DEPTH <- as.numeric(depth)
          
          # Append the data to the dataframe
          df_sm <- rbind(df_sm, sm)
        }
      }
    }
  }
  
  # Filter the dataframe
  df_sm <- df_sm %>% filter(VSWCFinalQF == 0)
  
  # Add the datetime column
  df_sm$datetime <- as.POSIXct(df_sm$endDateTime, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
  df_sm$datetime <- as.POSIXct(df_sm$datetime, format = "%Y-%m-%dT%H:%M:%S", tz = "EST")
  
  # Dynamically assign the dataframe to the list
  assign(df_name, df_sm)
  
  # Optionally, save the dataframe in the list if you want to keep them in a list
  df_list[[df_name]] <- df_sm
}


#top 2 layers all useable across all sensors, top 1 layer more usable

startdate = as.POSIXct("2017-07-01T19:30:00EST",format="%Y-%m-%dT%H:%M:%S", tz = "EST")
enddate = as.POSIXct("2023-07-26T19:00:00EST",format="%Y-%m-%dT%H:%M:%S", tz = "EST")



plot(df5$datetime,df5$VSWCMean,type ="l", ylim = c(0,.4), xlim = c(startdate,enddate),col="red",lwd=2)
lines(df2$datetime,df2$VSWCMean,col="orange",lwd=2)
lines(df3$datetime,df3$VSWCMean,col= "darkgreen",lwd=2)
lines(df4$datetime,df4$VSWCMean,col="blue",lwd=2)
lines(df1$datetime,df1$VSWCMean,col="purple",lwd=2)



#All time mean per sensor
df1_mean <- mean(df1$VSWCMean, na.rm =TRUE)
df2_mean <- mean(df2$VSWCMean, na.rm =TRUE)
df3_mean <- mean(df3$VSWCMean, na.rm =TRUE)
df4_mean <- mean(df4$VSWCMean, na.rm =TRUE)
df5_mean <- mean(df5$VSWCMean, na.rm =TRUE)

#Subtract mean from values to get anomaly
df1$VSWCAnom1 <- df1$VSWCMean - df1_mean
df2$VSWCAnom2  <-df2$VSWCMean - df2_mean
df3$VSWCAnom3 <- df3$VSWCMean - df3_mean
df4$VSWCAnom4  <- df4$VSWCMean - df4_mean
df5$VSWCAnom5 <- df5$VSWCMean - df5_mean

plot(df5$datetime,df5$VSWCAnom5,type ="l", ylim = c(-.5,.5),col="red",lwd=2 )
lines(df2$datetime,df2$VSWCAnom2,col="orange",lwd=2 )
lines(df3$datetime,df3$VSWCAnom3,col="darkgreen",lwd=2 )
lines(df4$datetime,df4$VSWCAnom4,col="blue",lwd=2 )
lines(df1$datetime,df1$VSWCAnom1,col="purple",lwd=2 )


#Get Mean of all anomalies
df_anom <- merge(x = df1, y = df2, by = "datetime", all=TRUE)
df_anom <- merge(x = df_anom, y = df3, by = "datetime", all=TRUE)
df_anom <- merge(x = df_anom, y = df4, by = "datetime", all=TRUE)
df_anom <- merge(x = df_anom, y = df5, by = "datetime", all=TRUE)

################################################
########df_anom is an important Dataframe#######
################################################

df_anom <- df_anom[c('datetime', 'VSWCAnom1', 'VSWCAnom2', 'VSWCAnom3', 'VSWCAnom4', 'VSWCAnom5')]
data <- df_anom[c('VSWCAnom1', 'VSWCAnom2', 'VSWCAnom3', 'VSWCAnom4', 'VSWCAnom5')]
df_anom$VSWCAnom.mean <- rowMeans(data,na.rm=TRUE)


plot(df_anom$datetime,df_anom$VSWCAnom.mean, type= "l", xlim = c(startdate,enddate),  ylim = c(-.1,.2))
df_anom$VSWCAnom.mean.med <- runmed(df_anom$VSWCAnom.mean, k=21)
plot(df_anom$datetime,df_anom$VSWCAnom.mean.med,type = "l", xlim = c(startdate,enddate),ylim = c(-.1,.2), lwd =2)


plot(df_anom$datetime,df_anom$VSWCAnom.mean.med, type = "l")




##RESPONSE OF SOIL MOISTURE TO RAIN (Day after)##
hf300_05_daily_m <- read_csv("~/Harvard_Forest/hf300-05-daily-m.csv")
hf300_05_daily_m$date <- as.POSIXct(hf300_05_daily_m$date, format = "%Y-%m-%d", tz = "UTC") 

df_anom$date <- substr(as.character(df_anom$datetime), start = 1, stop = 10)
df_anom$date <- as.POSIXct(df_anom$date, format = "%Y-%m-%d", tz = "UTC") 



#######################################################################
###############df_anom_daily is necessary for Evap_Perc.R##############
#######################################################################



df_anom_daily <- df_anom %>% group_by(date) %>% summarise(VSWC.mean = mean(VSWCAnom.mean.med, na.rm =TRUE))



#UTC is on purpose, both data sets are aligned temporally 

startdate = as.POSIXct("2017-06-01",format="%Y-%m-%d", tz = "UTC")
enddate = as.POSIXct("2017-09-30",format="%Y-%m-%d", tz = "UTC")

plot(df_anom_daily$date,df_anom_daily$VSWC.mean, xlim =c(startdate,enddate),type="l")

###Soil moisture anomaly in response to rain, shifted by one day
df_anom_precip <- merge(hf300_05_daily_m,df_anom_daily,by ="date")

df_anom_precip <- df_anom_precip  %>% mutate(VSWC.mean.shift = lead(VSWC.mean, 1, default = NA))

plot(df_anom_precip$prec,df_anom_precip$VSWC.mean)
plot(df_anom_precip$prec,df_anom_precip$VSWC.mean.shift)



#####Calculating Soil Moisture Response###


diff_list <- c()
for (x in 1:length(df_anom_daily$VSWC.mean)){
  diff <- df_anom_daily$VSWC.mean[x+1] - df_anom_daily$VSWC.mean[x-1]
  diff_list <- append(diff_list, diff)
}

#add a NA at the start
diff_list<- c(NA,diff_list)

df_anom_daily <- df_anom_daily  %>% mutate(VSWC.mean.shift = lead(VSWC.mean, 1, default = NA))


df_anom_daily$VSWC.diff <- diff_list
df_anom_precip <- merge(hf300_05_daily_m,df_anom_daily,by ="date")



######################################################
####Start here after generating soil moisture data####
######################################################

#what are the percentiles of rainfall by daily event?

hf300_05_daily_m <- read_csv("~/Harvard_Forest/hf300-05-daily-m.csv")

hf300_05_daily_m$year <- substr(as.character(hf300_05_daily_m$date), start = 1, stop = 4)


hf300_05_daily_m <- hf300_05_daily_m %>% filter(year >1972 )
hf300_05_daily_m <- hf300_05_daily_m %>% filter(prec > 1 )

rain_quant <- hf300_05_daily_m$prec %>% quantile(prob=c(.25,.5,.75,.90,.95), type=1)


#####REFINED RAIN EVENT SESNSING####
#using df anom and

setwd('/Users/jurado/Harvard_Forest')

hf001_10_15min_m <- read_csv("hf001-10-15min-m.csv")
df_anom_sm <- df_anom

# Create a new data frame with precipitation sums at 30-minute intervals
hf001_10_30min_m <- hf001_10_15min_m  %>%
  mutate(datetime = floor_date(datetime, "30 minutes")) %>%
  group_by(datetime) %>%
  summarize(precipitation_sum = sum(prec, na.rm = TRUE))

df_anom_sm_prec <- merge(df_anom_sm,hf001_10_30min_m, by ="datetime")

# Assuming your data frame is named df and has columns "datetime" and "precipitation"

# Initialize the new column with NA
df_anom_sm_prec$is_raining <- NA

# Loop through each row to set the value of is_raining
for (i in 1:nrow(df_anom_sm_prec)) {
  if (df_anom_sm_prec$precipitation[i] > 0) {
    # Check previous 6 rows
    if (i > 23 && all(df_anom_sm_prec$precipitation[(i-23):(i-1)] == 0)) {
      df_anom_sm_prec$is_raining[(i-6)] <- "START"
    }
    # Check next 48 rows
    if (i <= (nrow(df_anom_sm_prec) - 23) && all(df_anom_sm_prec$precipitation[(i+1):(i+23)] == 0)) {
      df_anom_sm_prec$is_raining[(i+23)] <- "STOP"
    }
  }
}



# Initialize new DataFrame
results <- data.frame(datetime = as.POSIXct(character()), precipitation = numeric(), 
                      initial_soil_moisture = numeric(), final_soil_moisture = numeric())

# Loop to calculate initial soil moisture, final soil moisture, precipitation sum and datetime
for (i in 1:nrow(df_anom_sm_prec)) {
  if (!is.na(df_anom_sm_prec$is_raining[i]) && df_anom_sm_prec$is_raining[i] == "START") {
    # Calculate initial soil moisture (mean of 6 rows including and after START)
    if (i + 5 <= nrow(df_anom_sm_prec)) {
      initial_soil_moisture <- mean(df_anom_sm_prec$VSWCAnom.mean.med[i:(i + 5)], na.rm = TRUE)
    } else {
      initial_soil_moisture <- NA
    }
    
    # Find the corresponding STOP and calculate values until then
    stop_index <- which(df_anom_sm_prec$is_raining == "STOP" & seq_along(df_anom_sm_prec$is_raining) > i)
    if (length(stop_index) > 0) {
      stop_index <- stop_index[1]
      # Calculate precipitation sum between START and STOP
      precip_sum <- sum(df_anom_sm_prec$precipitation_sum[i:(stop_index - 1)], na.rm = TRUE)
      
      # Calculate final soil moisture (mean of 48 rows up to STOP)
      if (stop_index + 23 <= nrow(df_anom_sm_prec)) {
        final_soil_moisture <- mean(df_anom_sm_prec$VSWCAnom.mean.med[(stop_index):(stop_index + 23)], na.rm = TRUE)
      } else {
        final_soil_moisture <- NA
      }
      
      # Add row to results
      new_row <- data.frame(
        datetime = df_anom_sm_prec$datetime[i],
        precipitation = precip_sum,
        initial_soil_moisture = initial_soil_moisture,
        final_soil_moisture = final_soil_moisture
      )
      results <- rbind(results, new_row)
    }
  }
}

results$sm_diff <- results$final_soil_moisture-results$initial_soil_moisture
results$prec_eff <- results$sm_diff/results$precipitation
results <- results %>% filter(precipitation >0)

mean(df_anom_sm_prec$VSWCAnom.mean.med)


quantile(results$precipitation, probs=c(.70,.75,.90))


results<- results[order(results$precipitation), ]

# Fit a LOESS model
loess_model <- loess(sm_diff~ precipitation, data = results, span = 0.3)  # span is a smoothing parameter

# Generate predictions with standard errors
pred <- predict(loess_model, newdata = results, se = TRUE)

# Add predictions and confidence intervals to the data frame
results$y_pred <- pred$fit


####Logistic fit###

SM.logistic = nls(sm_diff ~ (K * S0) / (S0 + (K - S0) * exp(-r * precipitation)),
                  start = list(S0 = .001, r = 0.09, K = 0.2),
                  data = results, trace = TRUE)

summary(SM.logistic)

results_test <- results %>% filter(initial_soil_moisture >0)




# Obtain the predicted values
predicted <- predict(SM.logistic, newdata = results)

# Create a data frame with the predictor and predicted values
plot_data <- data.frame(precipitation = results$precipitation,
                        predicted = predicted)

plot(plot_data$precipitation,plot_data$predicted)


# Provide initial parameter estimates
start_vals <- list(a = 10, b = 0.1)


plot(results$precipitation, results$sm_diff,
     col = ifelse(results$initial_soil_moisture > 0, "deepskyblue3", "lightsalmon4"),
     lwd = 1.5,
     pch = ifelse(results$initial_soil_moisture > 0, 16, 17),  # Different shapes for wet (16) and dry (17)
     xlab = "Precipitation [mm]",
     ylab = "Soil Moisture Anomaly Difference",
     cex.lab = 1.2  # Increase the size of the labels
)

abline(h = 0, lwd = 1.5, lty = 2)
abline(v = 19.75, lty = 3, col = "darkgrey", lwd = 1.5)
abline(v = 36.87, lty = 3, col = "darkgrey", lwd = 1.5)

par(new = TRUE)

lines(plot_data$precipitation, plot_data$predicted, lwd = 2)
abline(h = 0, lwd = 1.5, lty = 2)
abline(v = 19.75, lty = 3, col = "darkgrey", lwd = 1.5)
abline(v = 36.87, lty = 3, col = "darkgrey", lwd = 1.5)

legend(80, -0.03,
       legend = c("Wetter", "Drier"),
       col = c("deepskyblue3", "lightsalmon4"),
       pch = c(16, 17),  # Update shapes in legend to match plot
       ncol = 1, cex = 1, bty = "n", pt.lwd = 1.5)

title("Soil Moisture Response to Rain Events")
subtitle = "HARV June - September of 2017-2023"
mtext(subtitle)

text(15, 0.1, substitute(paste('75%')), col = "darkgrey")
text(32, 0.1, substitute(paste('90%')), col = "darkgrey")
text(60, 0.1, substitute(paste('Top 10% of Storms')), col = "darkgrey")



results_test <- results


results_test <- results_test %>% filter(precipitation < 9)

median(results_test$precipitation, na.rm =TRUE)
median(results_test$sm_diff, na.rm = TRUE)

cor.test(results_test$precipitation,results_test$sm_diff)

summary(lm(results_test$sm_diff~results_test$precipitation))



