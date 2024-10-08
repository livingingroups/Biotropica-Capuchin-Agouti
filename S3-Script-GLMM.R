### Script for 
###Modeling rates of agouti Capuchin terrestrial events at tool-use and non tool-use areas.
###in Coiba and Jicaron

rm(list = ls()) #this cleans the environment work space
shell("cls") #this cleans the console

##Load packages
library(MASS)
library(psych)
library(glmmADMB)
library(MCMCglmm)
library(stringr)
library(ggplot2)

##Read in data
#d.summed object created on S1-Script-Data-Prep
d.summed <- read.csv(file = "Processed_Data.csv", header = T, row.names = 1)

###Analysis###

#Exclude all mainland observations
d.summed <- subset(d.summed, sites == "Island")
#Exclude all arboreal camera trap observations
d.summed <- subset(d.summed, cam_loc == "ground")
#Remove unused levels of categorical variables.
d.summed <- droplevels(d.summed)
#Create a variable to count Capuchin and Agouti events. 
d.summed$events <- 1
d.summed$tool_site <- as.factor(d.summed$tool_site)

#Count events within Tool use sites, Camera Sites, Species Name, 
#Deployment Duration 
#The number of events is the variable to be analyzed, to estimate "rates".  
d.summed <- aggregate(events ~ tool_site + locationName + scientificName + 
                        dep_length_hours, 
                      data = d.summed, FUN = "sum")

#Calculate duration in months, so that rates per month can be easily calculated. 
d.summed$duration.mo <- d.summed$dep_length_hours/(24*30)

#Declaring the levels of the factor ToolUse, so that ToolUse is mapped by R to numeric values
#in the preferred order.
table(d.summed$tool_site, as.numeric(d.summed$tool_site)) 
d.summed$tool_site <- factor(d.summed$tool_site, 
                             levels = c("1", "0"))
table(d.summed$tool_site, as.numeric(d.summed$tool_site)) 

#Table of the number of deployments per Camera_Site.
table(d.summed$locationName)

#Calculate raw rates per deployment.
est.rates <- d.summed$events/d.summed$duration.mo
#Summaries of raw rates by tool use
describeBy(est.rates, group=d.summed$tool_site, digits=4)  

###Counting the number of  events in the data set.
sum(d.summed$events) ### 16382 events in the primary data set. 
table(d.summed$tool_site, d.summed$events)

###Getting the effort of camera deployment duration by toolsite
aggregate(duration.mo ~ tool_site, data = d.summed, FUN = "sum")
d.summed$locationName <- as.factor(d.summed$locationName)

### FIT NEGATIVE BINOMIAL MIXED MODEL. 

#Camera deployment duration is treated as an offset. 
m1 <- glmmadmb(events ~ -1 + tool_site + scientificName + tool_site:scientificName + 
                 offset(log(duration.mo)), 
               data = d.summed, random = ~1 | locationName, family = "nbinom")
summary(m1)

##

##We are calculating a constrast and the standard error of the contrast.
#Coefficient vector coef.vec will produce the difference between average rate on tool sites vs nontool sites.
coef.vec <- c(1/3, 1/3, -1, 1/3)
contrast <- coef.vec %*% (m1$b[1:4])
var.contrast <- t(coef.vec) %*% vcov(m1) %*% coef.vec
se.contrast <- sqrt(var.contrast)
##Contrast +/- SE(contrast) : 7.23 +/- 0.35. Since our value is higher than zero, the average rates between sites are not the same. 

#Estimated rates from negative binomial GLMM. Compare with output of "describeBy".
#The rate is number of events per month 
exp(m1$b[1:4])
str(m1)
#Estimated rates and confidence Intervals.
rates.ub <- m1$b[1:4] + 1.96*m1$stdbeta
rates.lb <- m1$b[1:4] - 1.96*m1$stdbeta
output<- cbind(exp(m1$b[1:4]), exp(rates.lb), exp(rates.ub))
colnames(output) <- c("Est.Rate", "95% Lower.bd", "95% Upper.bd")
round(output,3)

#Estimated rates and 84% confidence Intervals.
rates.ub84 <- m1$b[1:4] + 1.41*m1$stdbeta
rates.lb84 <- m1$b[1:4] - 1.41*m1$stdbeta
output84 <- cbind(exp(m1$b[1:4]), exp(rates.lb84), exp(rates.ub84))
colnames(output84) <- c("Est.Rate", "84% Lower.bd", "84% Upper.bd")
round(output84,3)

#Goodness of Fit.
#Obtaining "Pearson" residuals to check the goodness of fit of the model. 
res.pearson <- m1$residuals/m1$sd.est

#Calling the predictions for each value within the model.
predictions <- fitted(m1)

#Residual plot: log(predicted value) on the x-axis, Pearson residual on the y-axis. 
plot(log(predictions), res.pearson, pty="s", # pty = "s" for a square plotting region
     xlab="log(predicted value)", ylab="Pearson residual")

sum((res.pearson)^2) #147.7683 vs 181 obs
#The sum of the squared pearson residuals is called the Pearson Statistic for the model. 
#It is a single number that tells about the goodness of fit. 
#The statistic should be close to the total number of observations used to fit the model.

###GRAPHICAL DISPLAY OF MODEL ALONG WITH OBSERVATIONS.

#Predicted values.  
#Make a dataframe for generating predictions. 
#Using median values for duration
pred.frame <- data.frame(
  tool_site = c("1", "0"),
  scientificName = c("Cebus capucinus","Dasyprocta punctata"),
  duration.mo=median(d.summed$duration.mo))

pred.frame <- expand.grid(scientificName = c("Cebus capucinus","Dasyprocta punctata"), tool_site=c("1", "0"),
                          duration.mo = median(d.summed$duration.mo))
#Declare levels of Array in pred.frame. 
pred.frame$tool_site <- factor(pred.frame$tool_site, 
                           levels = c("1", "0"))
preds <- predict(m1, newdata= pred.frame, type="link", se.fit = TRUE)
preds.mean <- preds$fit
preds.lb <- preds$fit - 1.96*preds$se.fit
preds.ub <- preds$fit + 1.96*preds$se.fit
preds.lb84 <- preds$fit - 1.41*preds$se.fit
preds.ub84 <- preds$fit + 1.41*preds$se.fit

pred.frame <- cbind(pred.frame,preds.mean,preds.lb,preds.ub,preds.lb84,preds.ub84)
names(pred.frame) <- c("scientificName", "tool_site", "duration.mo", "mean", "lb", "ub", "lb84", "ub84")
pred.frame$ToolSite <- ifelse(pred.frame$tool_site == 0, "Non-tool use area", "Tool-use area")

#frame with raw data for plotting
df_raw_pred <- data.frame(d.summed$tool_site, d.summed$scientificName, 
                          d.summed$events)
names(df_raw_pred) <- c("tool_site", "scientificName", "events")
df_raw_pred$ToolSite <- ifelse(df_raw_pred$tool_site == 0, "Non-tool use area", "Tool-use area")
df_raw_pred$log <- log(df_raw_pred$events)

###Plot (FIGURE 3)
library(dplyr)
# Filter the raw data for Agouti
df_raw_agouti <- df_raw_pred %>%
  filter(scientificName == "Dasyprocta punctata")
df_raw_agouti$Species <- rep("Agouti", nrow(df_raw_agouti))
# Filter the raw data for Capuchin
df_raw_monkey <- df_raw_pred %>%
  filter(scientificName == "Cebus capucinus")
df_raw_monkey$Species <- rep("Capuchin", nrow(df_raw_monkey))



ggplot(pred.frame, aes(x = ToolSite, y = mean)) + 
  scale_x_discrete("") + 
  geom_hline(yintercept = 0, color = "red") + 
  theme_bw() + 
  theme(text = element_text(size = 10)) +
  ylab("Estimated Monthly Rates of Events (events/month)")  +
  facet_grid(~factor(scientificName, levels = c("Dasyprocta punctata","Cebus capucinus"),
                     labels = c("Agouti", "Capuchin"))) +
  geom_point(data = df_raw_agouti, aes(x = ToolSite, y = log), color = "grey", shape = 19, alpha = 0.4) +
  geom_point(data = df_raw_monkey, aes(x = ToolSite, y = log), color = "grey", shape = 19, alpha = 0.4) +
  geom_pointrange(aes(ymax = ub84, ymin = lb84), color = "black", lwd = 0.7)+
  scale_y_continuous(breaks = c(0,2.302585, 4.60517,6.214608), labels=c('0','10','100','500'))

