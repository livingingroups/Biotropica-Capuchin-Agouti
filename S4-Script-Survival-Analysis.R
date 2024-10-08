###Survival analysis script of waiting times for capuchin and agouti at Coiba/Jicarón
rm(list = ls()) #this cleans the environment work space
shell("cls") #this cleans the console

#Load packages
library(dplyr)
library(lubridate)
library(chron)
library(survival)

#Read in the data 
#d.summed object created on S1-Script-Data-Prep
d.summed <- read.csv(file = "Processed_Data.csv", header = T, row.names = 1)
#Exclude all mainland observations
d.summed <- subset(d.summed, sites == "Island")
#Exclude all arboreal camera trap observations
d.summed <- subset(d.summed, cam_loc == "ground")


#Format deployment start and end dates
d.summed$Date.time <- as.POSIXct(d.summed$DateTime, "%Y-%m-%d %H:%M:%S", tz = "EST")


d.summed$Date <- format(as.POSIXct(d.summed$Date.time), format = "%d-%m-%Y")
d.summed$Time <- format(as.POSIXct(d.summed$Date.time), format = "%H:%M:%S")

d.summed$datetime=chron(d.summed$Date, times=d.summed$Time, 
                  format=c(dates="d-m-Y", times="h:m:s"))

d.summed$datetime.num=as.numeric(d.summed$datetime)


d.summed$Start.Date.time <- as.POSIXct(d.summed$dep_start, "%Y-%m-%d %H:%M:%S", tz = "EST")
d.summed$StartDate <- format(as.POSIXct(d.summed$Start.Date.time), format = "%d-%m-%Y")
d.summed$StartTime <- format(as.POSIXct(d.summed$Start.Date.time), format = "%H:%M:%S")

d.summed$datetimestart=chron(d.summed$StartDate, times=d.summed$StartTime, 
                       format=c(dates="d-m-Y", times="h:m:s"))
d.summed$datetimestart.num=as.numeric(d.summed$datetimestart)

d.summed$End.Date.time <- as.POSIXct(d.summed$dep_end, "%Y-%m-%d %H:%M:%S", tz = "EST")
d.summed$EndDate <- format(as.POSIXct(d.summed$End.Date.time), format = "%d-%m-%Y")
d.summed$EndTime <- format(as.POSIXct(d.summed$End.Date.time), format = "%H:%M:%S")

d.summed$datetimeend =chron(d.summed$EndDate, times=d.summed$EndTime, 
                      format=c(dates="d-m-Y", times="h:m:s"))
d.summed$datetimeend.num=as.numeric(d.summed$datetimeend )


##Subset species#
capuchin <- subset(d.summed, scientificName == "Cebus capucinus")
agouti <- subset(d.summed, scientificName == "Dasyprocta punctata")




####Obtaining waiting times####
#Agouti waiting time (for each agouti observation we get the time interval to the previous capuchin):#####
agouti_previous_capuchin <- data.frame(time = rep(NA_real_, nrow(agouti)), status= rep(NA_real_, nrow(agouti)))
for (i in 1:nrow(agouti)) {
  # subsetting those caps from the same location
  capuchin.temp=capuchin[capuchin$locationName==agouti[i,"locationName"],] 
  if (nrow(capuchin.temp)>0) { #test if there are any capuchin obs at that particular location
    deltas = as.numeric(capuchin.temp$datetime - agouti[i,"datetime"])
    #All of the signed differences between the time of the focal agouti event and all capuchin events at that camera
    #min time between capture of agouti and capuchin
    if (all(deltas > 0)) { #Tells us the waiting time for the agouti event is censored (i.e., no capuchins before that particular agouti)
      agouti_previous_capuchin$status[i] = 0 #indicates censoring
      agouti_previous_capuchin$time[i] = as.numeric(agouti[i, "datetime"] - agouti[i, "datetimestart"]) #calculated censoring time (difference between focal obs and the camera deployment start)
    } #if True then waiting time is censored
    if (!all(deltas > 0)) { #There is a capuchin obs before the focal agouti event, i.e., not censored
      agouti_previous_capuchin$status[i] = 1 
      agouti_previous_capuchin$time[i] = min(abs(deltas[deltas < 0])) #uncensored waiting time is the minimum of the abs value of the negative deltas
    }
    
  } 
  if (!nrow(capuchin.temp)>0) { #for situations in which there were no capuchin observations ever at the camera for the focal agouti event (waiting time is thus censored)
    agouti_previous_capuchin$status[i] = 0
    agouti_previous_capuchin$time[i] = as.numeric(agouti[i, "datetime"] - agouti[i, "datetimestart"]) #calculated censoring time (difference between focal obs and the camera deployment start)
  } #if True then waiting time is censored
}
agouti_previous_capuchin <- cbind(agouti_previous_capuchin, agouti$tool_site, agouti$tooluse)
colnames(agouti_previous_capuchin) <- c("time", "status", "tool_site", "tool_use")

median(subset(agouti_previous_capuchin, status==1 & tool_site == 0)$time) #7.222025 days (173.3 hrs)
median(subset(agouti_previous_capuchin, status==1 & tool_site == 1)$time) #1.168032 days (28.0 hrs)


max(subset(agouti_previous_capuchin, status==1 & tool_site == 0)$time) #82.8344 days (1988.0 hrs)
max(subset(agouti_previous_capuchin, status==1 & tool_site == 1)$time) #81.85068 days (1964.4 hrs)

min(subset(agouti_previous_capuchin, status==1 & tool_site == 0)$time) #0.00162037 days (0.04 hrs -> 2.333333 minutes)
min(subset(agouti_previous_capuchin, status==1 & tool_site == 1)$time) #0.005578704 (0.13 hrs -> 8.033334 minutes)

quantile(subset(agouti_previous_capuchin, status == 1 & tool_site == 0)$time, probs = 0.25) #2.396814 days (57.52354 hrs)
quantile(subset(agouti_previous_capuchin, status == 1 & tool_site == 1)$time, probs = 0.25) #0.4216667 days (10.12 hrs)

quantile(subset(agouti_previous_capuchin, status == 1 & tool_site == 0)$time, probs = 0.75) #17.1788  days (412.291 hrs)
quantile(subset(agouti_previous_capuchin, status == 1 & tool_site == 1)$time, probs = 0.75) #2.845677 days (68.29625 hrs)

###Mann Whitey U test to assess significance
table(agouti_previous_capuchin$tool_site) 

agouti_previous_capuchin_toolsite <- subset(agouti_previous_capuchin, status==1 & tool_site == 1) 
agouti_previous_capuchin_non_toolsite <- subset(agouti_previous_capuchin, status==1 & tool_site == 0)

wilcox.test(agouti_previous_capuchin_toolsite$time,agouti_previous_capuchin_non_toolsite$time,
            conf.int = TRUE, conf.level = 0.84) # P < 2.2e-16 Difference -5.20 84% CI (-6.25 - -4.24)


#Capuchin waiting time (for each agouti observation we get the next capuchin):
#Agouti passes (focal event) --> get time interval to consecutive capuchin visit
agouti_next_capuchin <- data.frame(time = rep(NA_real_, nrow(agouti)), status= rep(NA_real_, nrow(agouti)))
for (i in 1:nrow(agouti)) {
  # subsetting those caps from the same location
  capuchin.temp1=capuchin[capuchin$locationName==agouti[i,"locationName"],] 
  if (nrow(capuchin.temp1)>0) { #test if there are any capuchin obs at that particular location
    deltas = as.numeric(capuchin.temp1$datetime - agouti[i,"datetime"])
    #All of the signed differences bw the time of the focal agouti event and all capuchin events at that camera
    #min time between capture of agouti and capuchin
    if (all(deltas < 0)) { #Tells us the waiting time for the capuchin event is censored (i.e., no capuchins after the focal agouti event)
      agouti_next_capuchin$status[i] = 0 #indicates censoring
      agouti_next_capuchin$time[i] = as.numeric(agouti[i, "datetimeend"] - agouti[i, "datetime"]) #calculated censoring time (difference between focal obs and the camera deployment end)
    } #if True then waiting time is censored
    if (!all(deltas < 0)) { #There is a capuchin observation after agouti passes (i.e. not censored)
      agouti_next_capuchin$status[i] = 1 
      agouti_next_capuchin$time[i] = min(abs(deltas[deltas > 0])) #uncensored waiting time is the minimum of the abs value of the positive deltas
    }
    
  } 
  if (!nrow(capuchin.temp1)>0) { #for situations in which there were no capuchin observations ever at the camera for the focal agouti event (waiting time is thus censored)
    agouti_next_capuchin$status[i] = 0
    agouti_next_capuchin$time[i] = as.numeric(agouti[i, "datetimeend"] - agouti[i, "datetime"]) #calculated censoring time (difference between focal obs and the camera deployment start)
  } #if True then waiting time is censored
}
agouti_next_capuchin <- cbind(agouti_next_capuchin, agouti$tool_site, agouti$tooluse)
colnames(agouti_next_capuchin) <-c("time", "status", "tool_site", "tool_use")



median(subset(agouti_next_capuchin, status==1 & tool_site == 0)$time) ##5.300822 days (127.2 hrs)
median(subset(agouti_next_capuchin, status==1 & tool_site == 1)$time) #0.9265394 days (22.2 hrs)

max(subset(agouti_next_capuchin, status==1 & tool_site == 0)$time) #returns Inf, 61.31779 days (1,471.6 hrs)
max(subset(agouti_next_capuchin, status==1 & tool_site == 1)$time) #23.13365 (555.2 hrs)

min(subset(agouti_next_capuchin, status==1 & tool_site == 0)$time) #0.002361111 days (0.06 hrs -> 3.4 minutes)
min(subset(agouti_next_capuchin, status==1 & tool_site == 1)$time) #0.001319444 days (0.03 hrs -> 1.9 minutes)

quantile(subset(agouti_next_capuchin, status == 1 & tool_site == 0)$time, probs = 0.25) #2.11265 days (50.7036 hrs)
quantile(subset(agouti_next_capuchin, status == 1 & tool_site == 1)$time, probs = 0.25) #0.3260532 days (7.825277 hrs)

quantile(subset(agouti_next_capuchin, status == 1 & tool_site == 0)$time, probs = 0.75) #13.10436  days (314.5046 hrs)
quantile(subset(agouti_next_capuchin, status == 1 & tool_site == 1)$time, probs = 0.75) #2.392292 days (57.41501 hrs)

###Mann Whitey U test to assess significance
table(agouti_next_capuchin$tool_site) 

agouti_next_capuchin_toolsite <- subset(agouti_next_capuchin, status==1 & tool_site == 1)
agouti_next_capuchin_non_toolsite <- subset(agouti_next_capuchin, status==1 & tool_site == 0)
##remove inf observations for the nontool site sample
agouti_next_capuchin_non_toolsite <- agouti_next_capuchin_non_toolsite[agouti_next_capuchin_non_toolsite$time != Inf, ]

wilcox.test(agouti_next_capuchin_toolsite$time,agouti_next_capuchin_non_toolsite$time,
            conf.int = TRUE, conf.level = 0.84) # P < 2.2e-16 Difference -3.92, 84% CI => -4.66 - -3.23



#Histograms of waiting times SUPPLEMENTAL FIGURE S9
par(mfrow = c(2,2),
    mar = c(5, 5, 4, 2))
hist(agouti_previous_capuchin_toolsite$time, xlab = "Time (days)", main = "Agouti waiting time (Tool use area)",
     breaks = seq(0,82, by = 1), ylim = c(0,60),cex.axis = 2, cex.main=2, cex.lab=2)
segments(x0 =  median(agouti_previous_capuchin_toolsite$time), y0 = 0, x1 =  median(agouti_previous_capuchin_toolsite$time), y1 = 55, col = "red", lwd = 2)
text(1,60,substitute(paste(bold("(a)"))),cex=3)

hist(agouti_previous_capuchin_non_toolsite$time, xlab = "Time (days)", main = "Agouti waiting time (Non-tool use area)",
     breaks = seq(0,83, by = 1), ylim = c(0,500),cex.axis = 2, cex.main=2, cex.lab=2)
segments(x0 =  median(agouti_previous_capuchin_non_toolsite$time), y0 = 0, x1 =  median(agouti_previous_capuchin_non_toolsite$time), y1 = 485, col = "red", lwd = 2)
text(1,490,substitute(paste(bold("(b)"))),cex=3)

hist(agouti_next_capuchin_toolsite$time, xlab = "Time (days)", main = "Capuchin waiting time (Tool use area)",
     breaks = seq(0,25, by = 1), ylim = c(0,70),cex.axis = 2, cex.main=2, cex.lab=2)
segments(x0 =  median(agouti_next_capuchin_toolsite$time), y0 = 0, x1 =  median(agouti_next_capuchin_toolsite$time), y1 = 60, col = "red", lwd = 2)
text(1,68,substitute(paste(bold("(c)"))),cex=3)

hist(agouti_next_capuchin_non_toolsite$time, xlab = "Time (days)", main = "Capuchin waiting time (Non-tool use area)",
     breaks = seq(0,62, by = 1), ylim = c(0,500),cex.axis = 2, cex.main=2, cex.lab=2)
segments(x0 =  median(agouti_next_capuchin_non_toolsite$time), y0 = 0, x1 =  median(agouti_next_capuchin_non_toolsite$time), y1 = 485, col = "red", lwd = 2)
text(1,490,substitute(paste(bold("(d)"))),cex=3)
graphics.off()


###Survival plots####

##models
summary(m1 <- survfit(Surv(time, status) ~ tool_site, data=agouti_previous_capuchin))
summary(m2 <- survfit(Surv(time, status) ~ tool_site, data=agouti_next_capuchin))




###Create object
model_output1 <- summary(m1)
lower_nontool <- model_output1$lower[model_output1$strata=="tool_site=0"] #lower for nontool site
upper_nontool <- model_output1$upper[model_output1$strata=="tool_site=0"] 
time_nontool <- model_output1$time[model_output1$strata=="tool_site=0"] #times for nontool use

time_tool <- model_output1$time[model_output1$strata=="tool_site=1"] #times for nontool use
lower_tool <- model_output1$lower[model_output1$strata=="tool_site=1"]
upper_tool <- model_output1$upper[model_output1$strata=="tool_site=1"]

##Figure 4
par(mfrow= c(1,2))
plot(m1,
     lwd=2, bty="n",
     col=c("black", "tomato3"),
     xlab="time t (days)", ylab="Proportion of waiting times > t", main = "Agouti waiting time") #Capuchin waiting time in relation to agouti
lines(x= time_nontool, y=lower_nontool, col = "gray15", lwd= 0.5)
lines(x= time_nontool, y=upper_nontool, col = "gray15", lwd= 0.5)
lines(x= time_tool, y=lower_tool, col = "tomato3", lwd= 0.5)
lines(x= time_tool, y=upper_tool, col = "tomato3", lwd= 0.5)
mtext(substitute(paste(bold("(a)"))), cex = 1.1, at=1)



model_output2 <- summary(m2)
lower_nontool2 <- model_output2$lower[model_output2$strata=="tool_site=0"] #lower for nontool site
upper_nontool2 <- model_output2$upper[model_output2$strata=="tool_site=0"] 
time_nontool2 <- model_output2$time[model_output2$strata=="tool_site=0"]

time_tool2 <- model_output2$time[model_output2$strata=="tool_site=1"] #times for nontool use
lower_tool2 <- model_output2$lower[model_output2$strata=="tool_site=1"]
upper_tool2 <- model_output2$upper[model_output2$strata=="tool_site=1"]


plot(m2,
     lwd=2, bty="n",
     col=c("black", "tomato3"),
     xlab="time t (days)", ylab="Proportion of waiting times > t", main = "Capuchin waiting time")
legend( "topright",
        c("Tool-use area","Non-tool use area"),
        text.col=c("tomato3","black"),
        bty="n" )

lines(x= time_nontool2, y=lower_nontool2, col = "gray60", lwd = 0.5)
lines(x= time_nontool2, y=upper_nontool2, col = "gray60", lwd = 0.5)
lines(x= time_tool2, y=lower_tool2, col = "tomato3", lwd = 0.5)
lines(x= time_tool2, y=upper_tool2, col = "tomato3", lwd = 0.5)
mtext(substitute(paste(bold("(b)"))), cex = 1.1, at=1)


graphics.off()


#####Waiting times tool use no tool use at tool-use area only####
###Focal event capuchin passes --> time to the next agouti
###We are interested in tool using capuchins, hence we need to obtain waiting times with respect to the capuchin
###An event with a capuchin using tools may be followed by other capuchin events before the next agouti!
capuchin_tool <- subset(capuchin,tool_site=="1") #subset capuchin observations at tool_sites
agouti_previous_capuchin_tool_use_area_only <- data.frame(time = rep(NA_real_, nrow(capuchin_tool)), status= rep(NA_real_, nrow(capuchin_tool)))
for (i in 1:nrow(capuchin_tool)) {
  # subsetting those caps from the same location
  agouti.temp1=agouti[agouti$locationName==capuchin_tool[i,"locationName"],] 
  if (nrow(agouti.temp1)>0) { #test if there are any agouti obs at that particular location
    deltas = as.numeric(agouti.temp1$datetime - capuchin_tool[i,"datetime"])
    #All of the signed differences between the time of the focal capuchin event and all agouti events at that camera
    #min time between capture of agouti and capuchin
    if (all(deltas < 0)) { #Tells us the waiting time for the focal capuchin event is censored (i.e., no agouti after that particular capuchin)
      agouti_previous_capuchin_tool_use_area_only$status[i] = 0 #indicates censoring
      agouti_previous_capuchin_tool_use_area_only$time[i] = as.numeric(capuchin_tool[i, "datetimeend"] - capuchin_tool[i, "datetime"]) #calculated censoring time (difference between focal obs and the camera deployment start)
    } #if True then waiting time is censored
    if (!all(deltas < 0)) { #AThere is an agouti obs after the focal capuchin event, i.e., not censored
      agouti_previous_capuchin_tool_use_area_only$status[i] = 1 
      agouti_previous_capuchin_tool_use_area_only$time[i] = min(deltas[deltas > 0]) #uncensored waiting time is the minimum of the abs value of the positive deltas
    }
    
  } 
  if (!nrow(agouti.temp1)>0) { #for situations in which there were no agouti observations ever at the camera  (waiting time is thus censored)
    agouti_previous_capuchin_tool_use_area_only$status[i] = 0
    agouti_previous_capuchin_tool_use_area_only$time[i] = as.numeric(capuchin_tool[i, "datetimeend"] - capuchin_tool[i, "datetime"]) #calculated censoring time (difference between focal obs and the camera deployment end)
  } #if True then waiting time is censored
}
agouti_previous_capuchin_tool_use_area_only <- cbind(agouti_previous_capuchin_tool_use_area_only, capuchin_tool$tool_site, capuchin_tool$tooluse)
colnames(agouti_previous_capuchin_tool_use_area_only) <- c("time", "status", "tool_site", "tool_use")


median(subset(agouti_previous_capuchin_tool_use_area_only, status==1 & tool_site == 1 & tool_use == TRUE)$time) #38.84034 days (932.2 hrs)
median(subset(agouti_previous_capuchin_tool_use_area_only, status==1 & tool_site == 1 & tool_use == FALSE)$time) #28.13593 days (675.3 hrs)
       
max(subset(agouti_previous_capuchin_tool_use_area_only, status==1 & tool_site == 1 & tool_use == TRUE)$time) #181.5738 days (4357.8 hrs)
max(subset(agouti_previous_capuchin_tool_use_area_only, status==1 & tool_site == 1 & tool_use == FALSE)$time) #222.5342 days (5340.8 hrs)

min(subset(agouti_previous_capuchin_tool_use_area_only, status==1 & tool_site == 1 & tool_use == TRUE)$time) #0.5211458 days (12.5 hrs)
min(subset(agouti_previous_capuchin_tool_use_area_only, status==1 & tool_site == 1 & tool_use == FALSE)$time) # 0.005578704 days (0.13 hrs -> 8.0 min) 

quantile(subset(agouti_previous_capuchin_tool_use_area_only, status == 1 & tool_site == 1 & tool_use == TRUE)$time, probs = 0.25) #9.285127 days (222.843 hrs)
quantile(subset(agouti_previous_capuchin_tool_use_area_only, status == 1 & tool_site == 1 & tool_use == TRUE)$time, probs = 0.75) #59.71377 days (1433.13 hrs)

quantile(subset(agouti_previous_capuchin_tool_use_area_only, status == 1 & tool_site == 1 & tool_use == FALSE)$time, probs = 0.25) #9.050237 days (217.2057 hrs)
quantile(subset(agouti_previous_capuchin_tool_use_area_only, status == 1 & tool_site == 1 & tool_use == FALSE)$time, probs = 0.75) #56.90134  days (1365.632 hrs)


###Mann Whitey U test to assess significance
table(agouti_previous_capuchin_tool_use_area_only$tool_use) #8371 events with no tool use and 1841 with tool use

agouti_previous_capuchin_tool_use_area_only_tooluse <- subset(agouti_previous_capuchin_tool_use_area_only, status==1 & tool_use == 1)
agouti_previous_capuchin_tool_use_area_only_non_tooluse <- subset(agouti_previous_capuchin_tool_use_area_only, status==1 & tool_use == 0)

wilcox.test(agouti_previous_capuchin_tool_use_area_only_tooluse$time,agouti_previous_capuchin_tool_use_area_only_non_tooluse$time,
            conf.int = TRUE, conf.level = 0.84) # P = 0.07 #84% CI (0.58 - 4.67)


###Figure 5####

##model

summary(m3 <- survfit(Surv(time, status) ~ tool_use, data=agouti_previous_capuchin_tool_use_area_only))


model_output3 <- summary(m3)

lower_nontool3 <- model_output3$lower[model_output3$strata=="tool_use=FALSE"] #lower for nontool site
upper_nontool3 <- model_output3$upper[model_output3$strata=="tool_use=FALSE"] 
time_nontool3 <- model_output3$time[model_output3$strata=="tool_use=FALSE"]

time_tool3 <- model_output3$time[model_output3$strata=="tool_use=TRUE"] #times for nontool use
lower_tool3 <- model_output3$lower[model_output3$strata=="tool_use=TRUE"]
upper_tool3 <- model_output3$upper[model_output3$strata=="tool_use=TRUE"]


plot(m3,
     lwd=2, bty="n",
     col=c("black", "tomato3"),
     xlab="time t (days)", ylab="Proportion of waiting times > t", main = "Agouti waiting time at tool-use area")
legend( "topright",
        c("Tool use","No tool use"),
        text.col=c("tomato3","black"),
        bty="n" )

lines(x= time_nontool3, y=lower_nontool3, col = "gray60", lwd = 0.5)
lines(x= time_nontool3, y=upper_nontool3, col = "gray60", lwd = 0.5)
lines(x= time_tool3, y=lower_tool3, col = "tomato3", lwd = 0.5)
lines(x= time_tool3, y=upper_tool3, col = "tomato3", lwd = 0.5)
graphics.off()
