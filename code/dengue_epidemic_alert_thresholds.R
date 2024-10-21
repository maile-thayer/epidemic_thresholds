rm(list=ls())

require(tidyverse)
require(readxl)
require(lubridate)
require(mgcv)
library(ggplot2)
library(RColorBrewer)

colors <- brewer.pal(8, name = 'Dark2')

setwd("//cdc.gov/locker/OID_DB_EPI_PII/Datasets/Projects/Dengue Threshold")


###############################################################################################################
#load confirmed/probable case data
df <- read_xlsx('matrix_dengue confirmed and probable cases_1986_to_2024-06-11.xlsx')[,-1] %>%
  # select(-OBSNUM) %>%
  pivot_longer(!matches('weekDate'), names_to='gcode', values_to='conf_cases') %>%
  group_by(weekDate) %>%
  summarize(conf_cases = sum(conf_cases)) %>%
  ungroup() %>%
  mutate(weekdate = as.Date(weekDate),
         week = week(weekdate),
         year = year(weekdate)) %>%
  filter(week != 53)


##################################################################################################

### calculate post-zika thresholds
threshold_evol <- tibble()
# for (yr in 2001:2024) {
for (yr in 2001:2024) {
  for (wk in 1:52) {
    this_nb <- MASS:::glm.nb(conf_cases ~ 1,
                             data=filter(df, week == wk, weekdate < as.Date(paste0(yr, '-01-01'))))
    threshold_evol <- bind_rows(threshold_evol,
                                tibble(
                                  year = yr,
                                  week = wk,
                                  q50 = qnbinom(0.5, size=this_nb$theta, 
                                                mu=exp(this_nb$coefficients[[1]])), # median
                                  q75 = qnbinom(0.75, size=this_nb$theta, 
                                                mu=exp(this_nb$coefficients[[1]])), #threshold
                                  q60 = qnbinom(0.60, size=this_nb$theta, 
                                                mu=exp(this_nb$coefficients[[1]])), #threshold low
                                  q90 = qnbinom(0.90, size=this_nb$theta, 
                                                mu=exp(this_nb$coefficients[[1]])), #threshold high
                                ))
  }
}


# Join datasets and plot
df_thresh <- mutate(df, year = year(weekdate)) %>% 
  left_join(threshold_evol)

df_thresh$cases_nonmissing <- df_thresh$conf_cases
df_thresh$conf_cases[c(2003:2028)] <- NA #change 0s in the future (after data are pulled) to missing for plotting

####################




##########################################################
######## SMOOTH THRESHOLDS AND COUNT 2+ WEEKS#############

smooth.filt5 <- c(1, 2, 3, 2, 1)/9 #7
smooth.filt7 <- c(1, 2, 3, 4, 3, 2, 1)/16 #7
smooth.filt9 <- c(1, 2, 3, 4, 5, 4, 3, 2, 1)/25 #9
smooth.filt <- smooth.filt11 <- c(1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1)/36 #11 
smooth.filt13 <- c(1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2, 1)/49 #13


smooth.filt <- rep(1/11,11)
smooth.filt <- c(11:1)/(sum(11:1))

# differing definitions for thresholds - smoothed
df_thresh$q75.sm <- head((stats::filter(c(df_thresh$q75,df_thresh$q75[(which(df_thresh$year==2024))[1:5]]), smooth.filt, circular=T)),-5)
df_thresh$q50.sm <- head((stats::filter(c(df_thresh$q50,df_thresh$q50[(which(df_thresh$year==2024))[1:5]]), smooth.filt, circular=T)),-5)
df_thresh$q90.sm <- head((stats::filter(c(df_thresh$q90,df_thresh$q90[(which(df_thresh$year==2024))[1:5]]), smooth.filt, circular=T)),-5)
df_thresh$q60.sm <- head((stats::filter(c(df_thresh$q60,df_thresh$q60[(which(df_thresh$year==2024))[1:5]]), smooth.filt, circular=T)),-5)


# smoothed threshold with different smoothing filters
df_thresh$q75.sm5 <- head((stats::filter(c(df_thresh$q75,df_thresh$q75[(which(df_thresh$year==2024))[1:2]]), smooth.filt5, circular=T)),-2)
df_thresh$q75.sm7 <- head((stats::filter(c(df_thresh$q75,df_thresh$q75[(which(df_thresh$year==2024))[1:3]]), smooth.filt7, circular=T)),-3)
df_thresh$q75.sm9 <- head((stats::filter(c(df_thresh$q75,df_thresh$q75[(which(df_thresh$year==2024))[1:4]]), smooth.filt9, circular=T)),-4)
df_thresh$q75.sm11 <- head((stats::filter(c(df_thresh$q75,df_thresh$q75[(which(df_thresh$year==2024))[1:5]]), smooth.filt11, circular=T)),-5)
df_thresh$q75.sm13 <- head((stats::filter(c(df_thresh$q75,df_thresh$q75[(which(df_thresh$year==2024))[1:6]]), smooth.filt13, circular=T)),-6)


# df_thresh$abovethresh.sm <- (df_thresh$cases_nonmissing>df_thresh$q60.sm)
# df_thresh$abovethresh <- (df_thresh$cases_nonmissing>df_thresh$q60)
df_thresh$abovethresh.sm <- (df_thresh$cases_nonmissing>df_thresh$q75.sm)
df_thresh$abovethresh <- (df_thresh$cases_nonmissing>df_thresh$q75)
# df_thresh$abovethresh.sm <- (df_thresh$cases_nonmissing>df_thresh$q90.sm)
# df_thresh$abovethresh <- (df_thresh$cases_nonmissing>df_thresh$q90)


df_thresh$count <-sequence(rle(as.character(df_thresh$abovethresh))$lengths)

df_thresh$count_cum <- with(rle(df_thresh$abovethresh), rep(values * cumsum(values & lengths >= 2),lengths))

init_y <- 1986 #starting year
thresh_init_y <- 2001 #first year of thresholds
curr_y <- 2024 #last observed year in dataset (ok if incomplete)

df_thresh$count_y <- NA
df_thresh$count_cum_y <- NA
for (y in 1:length(table(df_thresh$year))){
  df_thresh$count_y[which(df_thresh$year==(y+(init_y - 1)))] <-sequence(rle(as.character(df_thresh$abovethresh[which(df_thresh$year==(y+(init_y - 1)))]))$lengths)
  df_thresh$count_cum_y[which(df_thresh$year==(y+(init_y - 1)))] <- with(rle(df_thresh$abovethresh[which(df_thresh$year==(y+(init_y - 1)))]), rep(values * cumsum(values & lengths >= 2),lengths))
}

outbrk_y <- list()
# for (y in (thresh_init_y ):((init_y - 1)+length(table(df_thresh$year)))){
for (y in (thresh_init_y ):(curr_y)){
  outbrk_wk_tbl <- (table(df_thresh$count_cum_y[which(df_thresh$year==y)]))
  if (dimnames(outbrk_wk_tbl)[[1]][1]=="0") {outbrk_wk_tbl <- outbrk_wk_tbl[-1]}
  outbrk_num <- length(outbrk_wk_tbl)
  
  num_wks <- NULL
  if (outbrk_num>1) {
    for(i in 1:(outbrk_num)){
      num_wks <- c(num_wks,max(df_thresh$count_y[which((df_thresh$year==y) & 
                                                         df_thresh$count_cum_y == i)]))
    }} 
  if (outbrk_num == 1) {
    # if (names(outbrk_wk_tbl)=="1"){
    num_wks <- c(num_wks,max(df_thresh$count_y[which((df_thresh$year==y) & 
                                                       df_thresh$count_cum_y == 1)]))
  }#}
  if (outbrk_num == 0){
    num_wks <- c(num_wks,0)
    
  } 
  
  outbrk_y[[y-(thresh_init_y - 1)]] <- num_wks
  
  print(c(y,num_wks))
}

outbrk_y

ds_y <- data.frame(year=(thresh_init_y ):curr_y)
ds_y$num_outbrk_wks <- NA
ds_y$num_outbrks <- NA
ds_y$same_outbrk_y_b4 <- 0
ds_y$same_outbrk_y_after <- 0
for (i in 1:length((thresh_init_y ):curr_y)){
  ds_y$num_outbrk_wks[i] <- sum(outbrk_y[[i]])
  ds_y$num_outbrks[i] <- length(outbrk_y[[i]])
  if (i<length((thresh_init_y ):curr_y)){
    if ((df_thresh$count_y[which((df_thresh$year==(i+(thresh_init_y ))) & 
                                 (df_thresh$week==(52)))] > 0) &
        (df_thresh$count_y[which((df_thresh$year==(i+(thresh_init_y ))) & 
                                 (df_thresh$week==(1)))] > 0) & 
        (df_thresh$count_cum_y[which((df_thresh$year==(i+(thresh_init_y -1))) & 
                                     (df_thresh$week==(52)))] > 0 ) & 
        (df_thresh$count_cum_y[which((df_thresh$year==(i+(thresh_init_y -1 ))) & 
                                     (df_thresh$week==(1)))] > 0 )) {
      ds_y$same_outbrk_y_b4[i] <- ds_y$same_outbrk_y_after[i+1] <- 1
    }}
}

ds_y$num_outbrks[which(ds_y$num_outbrk_wks==0)] <- NA

r <- rle(df_thresh$abovethresh)

r$lengths >= 2 & r$values

r$values <- r$lengths >= 2 & r$values 

df_thresh$outbrk_flag <- inverse.rle(r)

# calculate excess weeks and cases
df_thresh$ex.all <- (df_thresh$cases_nonmissing - df_thresh$q75.sm)*(df_thresh$outbrk_flag) 

df_thresh$ex.all <- round((df_thresh$ex.all), 0)
df_thresh$ex.all[which(df_thresh$ex.all<0)] <- 0
df_thresh$ex.all[which(is.na(df_thresh$ex.all))] <- 0

excess.table <-df_thresh %>%
  group_by(year) %>%
  summarise(no.consec.epidemic.weeks = sum(outbrk_flag),
            excess.all = sum(ex.all))

excess.table <- full_join(excess.table,ds_y)


excess.table$no.consec.epidemic.weeks[which(excess.table$no.consec.epidemic.weeks<2)] <- 0



# excess.table.60=excess.table

excess.table.75=excess.table

# excess.table.90=excess.table

excess.table <- data.frame(year=rep(excess.table.60$year,3),
                           
                           no.consec.epidemic.weeks=c(excess.table.60$no.consec.epidemic.weeks, 
                                                      excess.table.75$no.consec.epidemic.weeks, 
                                                      excess.table.90$no.consec.epidemic.weeks),
                           
                           num_outbrk_wks=c(excess.table.60$num_outbrk_wks, 
                                            excess.table.75$num_outbrk_wks, 
                                            excess.table.90$num_outbrk_wks),
                           
                           threshold_percentile=c(rep("60th %ile",length(excess.table.60$year)),
                                                  rep("75th %ile",length(excess.table.60$year)),
                                                  rep("90th %ile",length(excess.table.60$year))))


excess.table_plot <- excess.table %>% 
  filter(year %in% c(2001:2024))
#########################################################################################################################################
# Figures

yr = 2024
df_thresh$weekdate <- as.Date(df_thresh$weekdate)

### Figure S1
ggplot(df_thresh, aes(x = weekdate)) +
  geom_col( data=df_thresh,
            aes(y = conf_cases), fill = "dodgerblue4", 
            col = "white",
            show.legend = TRUE,
            na.rm = TRUE
            # position = position_dodge(1)
            # width = 1
  ) + 
  geom_line(data=df_thresh, aes(x=weekdate, y=q75.sm5),color=colors[1], size=.5) + 
  geom_line(data=df_thresh, aes(x=weekdate, y=q75.sm7),color=colors[2], size=.5) +
  geom_line(data=df_thresh, aes(x=weekdate, y=q75.sm9),color=colors[3], size=.5) +
  geom_line(data=df_thresh, aes(x=weekdate, y=q75.sm11),color=colors[4], size=.5) +
  geom_line(data=df_thresh, aes(x=weekdate, y=q75.sm13),color=colors[5], size=.5) +
  scale_x_date( date_breaks = "1 month",
                labels=scales::date_format("%b"),
                limits = as.Date(c(paste0(yr,"-01-01"),paste0(yr,"-12-31")))) +
  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(#title="Reported dengue cases & thresholds: 2024",
    x ="Onset week", y = "Cases") 


colors2 <- c("Reported cases"="grey30",
             "60th percentile"=colors[1],
             "75th percentile"=colors[2],
             "90th percentile"=colors[3]
)


#### FIGURE 2
df_thresh$weekdate <- as.Date(df_thresh$weekdate)

ggplot(df_thresh, aes(x = weekdate)) +
  geom_col( data=df_thresh,
            aes(y = conf_cases), fill = "grey", col = "grey30", na.rm = TRUE
  ) +  
  geom_line(data=df_thresh, aes(x=weekdate, y=q60.sm, color = "60th percentile"), size=.5) + 
  geom_line(data=df_thresh, aes(x=weekdate, y=q75.sm, color = "75th percentile"), size=.5) + 
  geom_line(data=df_thresh, aes(x=weekdate, y=q90.sm, color = "90th percentile"), size=.5) + 
  theme_classic() + 
  labs(
    x ="Onset week", y = "Cases", color = "Legend") +
  scale_color_manual(values=colors2)+ theme(legend.position="bottom")



#### FIGURE 3
ggplot(data=excess.table_plot, aes(y=num_outbrk_wks, x=factor(year), fill=factor(threshold_percentile))) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2")+
  ggtitle("Number of (2+) weeks when cases exceed threshold \n by threshold definition") +
  xlab("Year") +
  ylab("Number of weeks") +
  guides(fill=guide_legend(title="Threshold level")) +
  theme(axis.text.x = element_text(angle = 90, vjust=.5, hjust=1))




df_thresh_ex <- df_thresh %>%
  filter(year %in% c(2022,2024))
latest.wk <- as.Date("2024-06-30")
df_thresh_ex$year <- factor(df_thresh_ex$year,, levels = unique(df_thresh_ex$year))
df_thresh_ex$delays_xmin <- as.Date(NA)
df_thresh_ex$delays_xmax <- as.Date(NA)
df_thresh_ex$delays_ymin <- NA
df_thresh_ex$delays_ymax <- NA

df_thresh_ex$delays_xmin[76:79] <- as.Date("2024-06-30")-21
df_thresh_ex$delays_xmax[76:79] <- as.Date("2024-06-30")
df_thresh_ex$delays_ymin[76:79] <- -Inf
df_thresh_ex$delays_ymax[76:79] <- Inf

df_thresh_ex$weekDate[1:52] <- df_thresh_ex$weekDate[53:104]
df_thresh_ex$weekdate[1:52] <- df_thresh_ex$weekdate[53:104]

### FIGURE 4


ggplot(df_thresh_ex, aes(x = weekdate)) +
  geom_col( data=df_thresh_ex,
            aes(y = conf_cases), fill = "dodgerblue4", 
            col = "white",
            show.legend = TRUE,
            na.rm = TRUE
            # position = position_dodge(1)
            # width = 1
  ) + 
  geom_line(data=df_thresh_ex, aes(x=weekdate, y=q50.sm),col="grey60", linetype = "dashed", size=.8) +
  geom_line(data=df_thresh_ex, aes(x=weekdate, y=q75.sm),color="firebrick", linetype = "dashed", size=.8) +
  # 
  # geom_rect(data=data.frame (xmin=(latest.wk-3), xmax=latest.wk, ymin=-Inf, ymax=Inf), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey60", alpha=0.4, inherit.aes = FALSE) +
  geom_rect(data=df_thresh_ex, aes(xmin=delays_xmin, xmax=delays_xmax, ymin=delays_ymin, ymax=delays_ymax), fill="grey60", alpha=0.1, inherit.aes = FALSE) +
  facet_grid(year ~ .)+
  # 
  scale_x_date( date_breaks = "1 month",
                labels=scales::date_format("%b"),
                limits = as.Date(c(paste0(yr,"-01-01"),paste0(yr,"-12-31")))) +

  theme_classic() +
  theme(
    # axis.text.x = element_text(angle = 90),
        strip.text = element_text(face="bold", size=14),
        strip.background = element_rect(fill="lightblue", colour="black")) +
  labs(#title="Reported dengue cases & thresholds: 2024",
    x ="Onset week", y = "Cases")



