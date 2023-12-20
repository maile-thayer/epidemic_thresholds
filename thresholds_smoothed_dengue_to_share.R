rm(list=ls())


# install.packages("tidyverse") #install these packages if you don't have them already
# install.packages("ggplot2")


require(tidyverse) #load packages needed
require(ggplot2)
require(lubridate)

setwd("~/Dengue/PAHO/example_threshold_code_and_data") #<-- set to where you save the downloaded folder with the code and data


################################################################################
#################################Load & prep data###############################
################################################################################

# Need two columns/variables: 
#    1) Weekdate
#    2) Total weekly case counts

#load and format dengue case data; create year and week variables
df <- read.csv('Data/San_Juan_data_example.csv')[,-1] %>%
  mutate(weekdate=as.Date(weekdate))  %>% #<-- format date variable
  mutate(week = week(weekdate),
         year = year(weekdate)) %>% #<-- create variables for week and year
  filter(week != 53) #<-- exclude last partial week if applicable


################################################################################
#############################Calculate thresholds###############################
################################################################################

### calculate ahead one year, going 1 week at a time
### can take a little while to run
threshold_evol <- tibble()
for (yr in 2000:2009) { # <-- set to years when you want to create thresholds -- need at least a few years of burn in, so don't set it at the first year of data you have 
  for (wk in 1:52) {
    this_nb <- MASS:::glm.nb(total_cases ~ 1,
                             data=filter(df, week == wk, weekdate < as.Date(paste0(yr, '-01-01'))))
    threshold_evol <- bind_rows(threshold_evol,
                                       tibble(
                                         year = yr,
                                         week = wk,
                                         q50 = qnbinom(0.5, size=this_nb$theta, 
                                                       mu=exp(this_nb$coefficients[[1]])), # median
                                         q75 = qnbinom(0.75, size=this_nb$theta, 
                                                       mu=exp(this_nb$coefficients[[1]])) #threshold-- need to decide what level to use***
                                       ))
  }
}
# takes a little while to run
# if you receive an error, you likely need to set the starting year later (lin 38)
# if you receive a warning, you can ignore it

# Join datasets 
df_thresh <- mutate(df, year = year(weekdate)) %>% 
  left_join(threshold_evol)

df_thresh$cases_missing <- df_thresh$total_cases #keep current unobserved weeks as missing cases for plotting later
df_thresh$cases_nonmissing <- df_thresh$total_cases
df_thresh$cases_nonmissing[which(is.na(df_thresh$cases_nonmissing))] <- 0 #replace cases that haven't been reported yet with 0s (can't have missing values in some of the functions later)

################################################################################
###################SMOOTH THRESHOLDS AND COUNT OUTBREAK WEEKS###################
################################################################################

################################################################################
# in Puerto Rico, epidemics are defined as those that go above threshold for 
# AT LEAST 2 consecutive weeks. If a country has a similar definition, we can add 
# that to the code. I have excluded that here, and "outbreak weeks" are defined
# as ANY week when cases are above the threshold
################################################################################

year_desired <- 2009 #<-- set as current year, or year of interest
prev_year <- year_desired -1

# smoothing triangle
smooth.mo <- c(1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1)/36 #<-- you may need to play around with this triangular smoothing filter

#smooth and wrap (reported) data at end to extend smoothed thresholds through 1 year ahead
df_thresh$q75.sm <- head((stats::filter(c(df_thresh$q75,df_thresh$q75[(which(df_thresh$year==year_desired))[1:5]]), smooth.mo, circular=T)),-5)
df_thresh$q50.sm <- head((stats::filter(c(df_thresh$q50,df_thresh$q50[(which(df_thresh$year==year_desired))[1:5]]), smooth.mo, circular=T)),-5)
df_thresh$q75.sm.round <- floor(df_thresh$q75.sm) #round 
df_thresh$q50.sm.round <- floor(df_thresh$q50.sm)

df_thresh$abovethresh.sm <- (df_thresh$cases_nonmissing>df_thresh$q75.sm) #<-- identify ocurences when cases exceed threshold
df_thresh$abovethresh <- (df_thresh$cases_nonmissing>df_thresh$q75)


# calculate excess weeks and cases
df_thresh$ex.all <- (df_thresh$cases_nonmissing - df_thresh$q75.sm)*(df_thresh$abovethresh.sm) 
df_thresh$ex.all <- round((df_thresh$ex.all), 0) #<-- round to whole number
df_thresh$ex.all[which(df_thresh$ex.all<0)] <- 0 #<-- replace negatives and NAs as 0
df_thresh$ex.all[which(is.na(df_thresh$ex.all))] <- 0

excess.table <- aggregate(cbind(df_thresh$abovethresh, df_thresh$ex.all), list(year=df_thresh$year), sum)
names(excess.table) <- c("year","no.epidemic.weeks", "excess.all")
excess.table <- excess.table[which(excess.table$year>= 2000),]


################################################################################
##############################Create plots and save data########################
################################################################################



# plot full dataset with thresholds
plot(select(df_thresh, weekdate, total_cases), type='l', ylab="Weekly cases", xlab="Weekdate",lwd=3, bty="n")
lines(select(df_thresh, weekdate, q75.sm), lty=3, col="red",lwd=2)
legend("topright" ,c("Reported cases","Threshold"),col=c("black","red"),lty=c(1,3),lwd=c(3,2))

# plot full dataset with thresholds and identify weeks when cases exceed threshold
plot(select(df_thresh, weekdate, total_cases), type='l', ylab="Weekly cases", xlab="Weekdate",lwd=3, bty="n")
lines(select(df_thresh, weekdate, q75.sm), lty=3, col="red",lwd=2)
points(df_thresh$weekdate[df_thresh$abovethresh],df_thresh$total_cases[df_thresh$abovethresh],col="blue",pch=16)
legend("topright" ,c("Reported cases","Threshold","Weeks when cases exceed threshold"),col=c("black","red","blue"),lty=c(1,3,NA),lwd=c(3,2,NA),pch=c(NA,NA,16),cex=.8)

# plot number of weeks when cases exceed threshold
ggplot(excess.table, aes(y=no.epidemic.weeks, x=factor(year))) +
  geom_bar(position=position_dodge(width = .75), stat = "identity", fill="steelblue") +
  theme_classic() +
  ggtitle("Number of weeks when dengue cases exceed threshold") +
  xlab("Year") +
  ylab("Number of weeks") +
  theme(axis.text.x = element_text(angle = 90, vjust=.5, hjust=1))


# plot year of interest only

plot(df_thresh$weekdate, df_thresh$q75.sm, 
     type='n', ylab="", xlab="", 
     ylim=(c(0, 100)), bty="n", col="red", axes=F, cex=2, xlim=as.Date(c(paste0(year_desired,"-01-01"),paste0(year_desired,"-12-31"))))
lines(df_thresh$weekdate, 
      df_thresh$q75.sm, col="red",lwd=2,lty=2)
lines(df_thresh$weekdate, 
      df_thresh$q50.sm, col="grey60",lwd=2,lty=2)
lines(df_thresh$weekdate, df_thresh$total_cases,col="black", lwd=3)
datemo <- as.Date(c(paste(year_desired, 1:12, '01', sep='-')))
axis.Date(1, at=seq(min(datemo), max(datemo), by="months"), format="%m-%Y",line=-.3)
axis(2, line=0.7)
title( #"2022-present dengue thresholds (smoothed)",
  xlab="Weekdate", ylab="Cases")
legend("topleft",c("Threshold","Historical median","Reported cases"),col=c("red","grey60","black"),lty=c(2,2,1),lwd=c(2,2,3),cex=.8)


# Export full dataset
df_thresh_save <- df_thresh[,-c(7,8,11,12,13,15)]
names(df_thresh_save)[5:8] <- c("historical_median","threshold","historical_median_smoothed","threshold_smoothed")
write.csv(df_thresh_save,paste0("thresholds_all_", Sys.Date(), ".csv"))

# Export only desired year results as a table
df_thresh_current <- df_thresh_save[which(df_thresh_save$year==year_desired),]
write.csv(df_thresh_current,paste0("thresholds_",year_desired,"_only_", Sys.Date(), ".csv"))


