# Measure growing season yield for June-Sept

# Met data - daily
met = read_csv("~/Documents/flux_data/met/hf300-05-daily-m.csv")

met$year = lubridate::year(met$date)
met$month = lubridate::month(met$date)
met$wyear = ifelse(met$month ==12, met$year+1, met$year)

met_ann = met %>%
  group_by(wyear) %>%
  summarize(Pmm = sum(prec, na.rm=T)) %>%
  filter(wyear >= 2008)

# watershed area
blo_area = 650000 #m2

# Read in streamflow data
streams = read_csv("~/Documents/field-data/hf070-04-15min-Dec24.csv") %>%
  mutate(year = lubridate::year(datetime),
         month = lubridate::month(datetime),
         wyear = ifelse(month>=11, year+1, year))

# Relationship between growing season daily discharge ~ precip
hydro_flow = streams %>%
  mutate(date = lubridate::date(datetime)) %>%
  filter(month >= 6 & month <= 9) %>%
  dplyr::group_by(wyear, date) %>%
  dplyr::summarize(alo_dis = sum(al.dis,na.rm=T)*(60*15)/blo_area) %>%
  left_join(met) 

# 95%tile
quantile(hydro_flow$prec, 0.95,na.rm=T)
hydro_flow$P_95 = ifelse(hydro_flow$prec>=24.49,1,0)

# Flow on rolling median 7-day window
window = 7 # days to calculate rolling median
hydro_roll = hydro_flow %>%
  arrange(date) %>%
  group_by(wyear) %>%
  mutate(alodis_2 = RcppRoll::roll_mean(alo_dis, window, 
                                       align = "left", na.rm=T, fill = NA),
    P_mm = RcppRoll::roll_median(prec, window, 
                                   align = "right", na.rm=T, fill = NA),
         P_sum = RcppRoll::roll_sum(prec, window, 
                                       align = "right", na.rm=T, fill = NA),
         P95 = RcppRoll::roll_sum(P_95, window, 
                                     align = "right", na.rm=T, fill = NA),
    dis_diff = alo_dis-alodis_2)

# estimate logistic model for runoff (y) & precip (x)
runoff_log <- nls(dis_diff ~ (K*N0)/(N0 + (K - N0)*
                                      exp(-1*r*P_sum)),
                  data = hydro_roll, 
                  start = list(N0 = -1, r = 0.1, K = 5))
summary(runoff_log)

# calculate logistic model for plot
x = seq(min(hydro_roll$P_sum,na.rm=T),
        max(hydro_roll$P_sum,na.rm=T),0.25)
y = (coef(runoff_log)[3]*coef(runoff_log)[1])/
  (coef(runoff_log)[1] + (coef(runoff_log)[3]-coef(runoff_log)[1])*
                        exp(-1*coef(runoff_log)[2]*x))

# Data frame of model estimated logistic curve
df = data.frame(x,y)

ggplot() +
  geom_point(data = filter(hydro_roll, P_sum>0),
             aes(P_sum, dis_diff, color=as.factor(P95))) +
  geom_line(aes(x=x, y=y),color="black")+
  labs(x = "Prior 7-day total precip (mm)",
       y = "Runoff anomaly difference (mm)",
       title = "Runoff with precip. event size") +
  scale_color_manual(values = c("lightgrey","skyblue1","skyblue3",
                                 "skyblue4","navyblue"),
                                name = "p95 events") +
  theme_linedraw() +
  theme(legend.position="bottom")

# Statistical t-test for precip response
tmp = filter(hydro_roll, !is.na(P_sum),!is.na(dis_diff))
t.test(tmp$dis_diff[tmp$P_sum>50],
       tmp$dis_diff[tmp$P_sum>0 & tmp$P_sum<50])

aov_flow = aov(dis_diff ~ as.factor(P95), 
               data = hydro_roll)
summary(aov_flow)
TukeyHSD(aov_flow)
