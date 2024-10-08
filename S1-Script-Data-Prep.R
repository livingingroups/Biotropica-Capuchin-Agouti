##Script for data prep 
###This script is for data prep
###For analysis replication please load object generated at end of this script and begin on S2-Script
library(dplyr)
library(lubridate)
library(psych)
library(stringr)
library(camtrapR)

###Read in data####

d <- read.csv("Supp.DATA_CoibaBioblitz.csv",
              header = T) #Survey 1 - Coiba Bioblitz

str(d)

d2 <- read.csv("Supp.DATA_capuchin_terrestriality.csv",
               header = T) ##Survey 2 - Capuchin terrestriality
str(d2)

d3 <-read.csv("Supp.DATA_tool_users.csv",
              header = T) ##Survey 3 - Capuchin tool use project
str(d3)


#Apply 2 minute filter to Survey 3
# from the R package camtrapR

assessTemporalIndependence <- function(intable,
                                       deltaTimeComparedTo,
                                       columnOfInterest,     # species/individual column
                                       cameraCol,
                                       camerasIndependent,
                                       stationCol,
                                       minDeltaTime,
                                       removeNonIndependentRecords = TRUE)
{
  # check if all Exif DateTimeOriginal tags were read correctly
  if(any(is.na(intable$DateTimeOriginal))){
    which.tmp <- which(is.na(intable$DateTimeOriginal))
    if(length(which.tmp) == nrow(intable)) stop("Could not read any Exif DateTimeOriginal tag at station: ", paste(unique(intable[which.tmp, stationCol])), " Consider checking for corrupted Exif metadata.")
    warning(paste("Could not read Exif DateTimeOriginal tag of", length(which.tmp),"image(s) at station", paste(unique(intable[which.tmp, stationCol]), collapse = ", "), ". Will omit them. Consider checking for corrupted Exif metadata. \n",
                  paste(file.path(intable[which.tmp, "Directory"],
                                  intable[which.tmp, "FileName"]), collapse = "\n")), call. = FALSE, immediate. = TRUE)
    intable <- intable[-which.tmp ,]
    rm(which.tmp)
  }
  
  # prepare to add time difference between observations columns
  intable <- data.frame(intable,
                        delta.time.secs  = NA,
                        delta.time.mins  = NA,
                        delta.time.hours = NA,
                        delta.time.days  = NA)
  
  # introduce column specifying independence of records
  if(minDeltaTime == 0) {
    intable$independent <- TRUE    # all independent if no temporal filtering
  } else {
    intable$independent <- NA
  }
  
  
  for(xy in 1:nrow(intable)){     # for every record
    
    # set independent = TRUE if it is the 1st/only  record of a species / individual
    
    if(camerasIndependent == TRUE){
      if(intable$DateTimeOriginal[xy]  == min(intable$DateTimeOriginal[which(intable[, columnOfInterest] == intable[xy, columnOfInterest] &
                                                                             intable[, stationCol]       == intable[xy, stationCol] &
                                                                             intable[, cameraCol]        == intable[xy, cameraCol]) ])){    # cameras at same station assessed independently
        intable$independent[xy]       <- TRUE
        intable$delta.time.secs[xy]   <- 0
      }
    } else {
      if(intable$DateTimeOriginal[xy]  == min(intable$DateTimeOriginal[which(intable[, columnOfInterest] == intable[xy, columnOfInterest] &
                                                                             intable[, stationCol]       == intable[xy, stationCol]) ])){
        intable$independent[xy]       <- TRUE
        intable$delta.time.secs[xy]   <- 0
      }
    }
    
    if(is.na(intable$delta.time.secs[xy])) {   # if not the 1st/only record, calculate time difference to previous records of same species at this station
      
      if(deltaTimeComparedTo == "lastIndependentRecord"){
        
        if(camerasIndependent == TRUE){
          which_time2 <- which(intable[, columnOfInterest]       == intable[xy, columnOfInterest] &    # same species/individual
                                 intable[, stationCol]              == intable[xy, stationCol] &          # at same station
                                 intable[, cameraCol]               == intable[xy, cameraCol] &           # at same camera
                                 intable$independent                == TRUE &                             # independent (first or only record of a species at a station)
                                 intable$DateTimeOriginal           <  intable$DateTimeOriginal[xy])      # earlier than record xy
        } else {
          which_time2 <- which(intable[, columnOfInterest]       == intable[xy, columnOfInterest] &
                                 intable[, stationCol]             == intable[xy, stationCol] &
                                 intable$independent               == TRUE &
                                 intable$DateTimeOriginal          <  intable$DateTimeOriginal[xy])
        }
      }  else {
        if(camerasIndependent  == TRUE){
          which_time2 <- which(intable[, columnOfInterest]       == intable[xy, columnOfInterest] &
                                 intable[, stationCol]             == intable[xy, stationCol] &
                                 intable[, cameraCol]              == intable[xy, cameraCol] &
                                 intable$DateTimeOriginal          <  intable$DateTimeOriginal[xy])
        } else {
          which_time2 <- which(intable[, columnOfInterest]       == intable[xy, columnOfInterest] &
                                 intable[, stationCol]             == intable[xy, stationCol] &
                                 intable$DateTimeOriginal          <  intable$DateTimeOriginal[xy])
        }
      }
      
      # time difference to last (independent) record
      diff_tmp <- min(na.omit(difftime(time1 = intable$DateTimeOriginal[xy],            # delta time to last independent record
                                       time2 = intable$DateTimeOriginal[which_time2],
                                       units = "secs")))
      
      # save delta time in seconds
      intable$delta.time.secs[xy] <-  diff_tmp
      if(intable$delta.time.secs[xy] >= (minDeltaTime * 60) | intable$delta.time.secs[xy] == 0){
        intable$independent[xy] <- TRUE
      } else {
        intable$independent[xy] <- FALSE
      }
      
    }   # end   if(intable$DateTimeOriginal[xy] == min(...)} else {...}
  }     # end for(xy in 1:nrow(intable))
  
  
  if(removeNonIndependentRecords){
    # keep only independent records
    outtable <- intable[intable$delta.time.secs >= (minDeltaTime * 60) |
                          intable$delta.time.secs == 0,]
  } else {
    outtable <- intable
  }
  
  return(outtable)
}

d3$DateTimeOriginal <- d3$DateTime
records_filter_2min <- camtrapR:::assessTemporalIndependence(intable = d3, 
                                                             deltaTimeComparedTo = "lastIndepentRecord",
                                                             columnOfInterest = "scientificName",
                                                             stationCol = "locationName",
                                                             minDeltaTime = 2,
                                                             camerasIndependent = FALSE)
d3 <- records_filter_2min[, c("scientificName", "locationName", "DateTime", "dep_start", "dep_end", "cam_loc", "sites", "suntime", "dep_length_hours", "tool_site", "tooluse")]

#Set up all cameras on Survey 1 and Survey 2 as "non-tool sites" and "No tool use"
d$tool_site <- rep(0, nrow(d))
d2$tool_site <- rep(0, nrow(d2))
d$tooluse <- rep(FALSE, nrow(d))
d2$tooluse <- rep(FALSE, nrow(d2))

##Select relevant columns for activity patterns analysis

d <- d[, c("scientificName", "locationName", "DateTime", "dep_start", "dep_end", "suntime", "cam_loc", "sites", "dep_length_hours", "tool_site", "tooluse")]
d2 <- d2[, c("scientificName", "locationName", "DateTime", "dep_start", "dep_end", "suntime", "cam_loc", "sites", "dep_length_hours", "tool_site", "tooluse")]
d3 <- d3[, c("scientificName", "locationName", "DateTime", "dep_start", "dep_end", "suntime", "cam_loc", "sites", "dep_length_hours", "tool_site", "tooluse")]


##Merge datasets
d.summed <- rbind(d,d2,d3)
str(d.summed)

#save object as csv for subsequent analyses
#write.csv(d.summed, file = "Processed_Data.csv")



##Camera operation plots####
#Survey 1 (Figure S4)
d <- read.csv("Supp.DATA_CoibaBioblitz.csv",
               header = T) ##Coiba Bioblitz project
str(d)

# Convert the date columns to Date objects
d$dep_start <- as.Date(d$dep_start)
d$dep_end <- as.Date(d$dep_end)
d$cam_id <- as.integer(as.factor(d$locationName))

#Sort cameras by Island
d$Island <- ifelse(grepl("^Coiba-1-", d$locationName), "Coiba", "Jicaron")

#png(file= "D:/MPI/Biotropica_Supp/Figures/16.10.2023_CameraOP_Survey1.png",
    #width = 6, height = 6, units = "in", res = 600, pointsize = 10)

# Create a color palette for Island
colors <- c("Coiba" = "tomato3", "Jicaron" = "blue")
# Define the start and end of the x-axis range
x_axis_start <- as.Date("2015-01-01")  # Adjust the start date as needed
x_axis_end <- as.Date("2019-11-01")    # Adjust the end date as needed
# Increase left margin to allow for longer y-axis labels
par(mar = c(5, 8, 4, 2) + 0.1)
# Plot setup
plot(d$dep_start, as.numeric(d$locationName), type = "n", xlim = range(d$dep_start-40, d$dep_end+120),
     ylim = c(1, length(unique(d$locationName))), xlab = "Month (Year 2015)", ylab = " ", yaxt = "n", main = "Survey 1")

# Loop through unique location names
unique_locations <- unique(d$locationName)
for (i in 1:length(unique_locations)) {
  subset_data <- d[d$locationName == unique_locations[i], ]
  island_colors <- colors[as.character(subset_data$Island)]
  segments(subset_data$dep_start, i, subset_data$dep_end, i, col = island_colors, lwd = 3)
}
legend(max(d$dep_end)+70,26, c("Coiba" , "Jicarón"), pch=15, col=colors, box.col=NA, cex=1 )
axis(2, at = 1:length(unique(d$cam_id)), labels = sort(unique(d$locationName)), las = 1, cex.axis = 0.5)
# Set the x-axis labels for each month
axis(1, at = seq(x_axis_start, x_axis_end, by = "1 month"), labels = format(seq(x_axis_start, x_axis_end, by = "1 month"), "%b"))
#graphics.off()

#Survey 2 (Figure S5)
d2 <- read.csv("Supp.DATA_capuchin_terrestriality.csv",
               header = T) ##Capuchin terrestriality project
str(d2)
# Convert the date columns to Date objects
d2$dep_start <- as.Date(d2$dep_start)
d2$dep_end <- as.Date(d2$dep_end)
d2$cam_id <- as.integer(as.factor(d2$locationName))

# Sort cameras by Site
d2$Site <- ifelse(grepl("^COI-", d2$locationName), "Coiba",
                  ifelse(grepl("^JIC-", d2$locationName), "Jicaron", "BCI"))
#png(file= "D:/MPI/Biotropica_Supp/Figures/16.10.2023_CameraOP_Survey2.png",
    #width = 6, height = 6, units = "in", res = 600, pointsize = 10)

# Create a color palette for Island
colors <- c("Coiba" = "tomato3", "Jicaron" = "blue", "BCI" = "grey20")

# Specify the order of unique locations
unique_locations <- unique(d2$locationName)
unique_locations <- unique_locations[order(match(unique_locations, unique_locations[duplicated(unique_locations)]), unique_locations)]

# Define the start and end of the x-axis range
x_axis_start <- as.Date("2019-01-01")  # Adjust the start date as needed
x_axis_end <- as.Date("2019-11-01")    # Adjust the end date as needed
# Increase left margin to allow for longer y-axis labels
par(mar = c(5, 8, 4, 2) + 0.1)
# Plot setup
plot(d2$dep_start, as.numeric(d2$locationName), type = "n", xlim = range(d2$dep_start - 60, d2$dep_end + 60),
     ylim = c(1, length(unique_locations)), xlab = "Month (Year 2019)", ylab = "", yaxt = "n", main = "Survey 2")

# Loop through unique location names
for (i in 1:length(unique_locations)) {
  subset_data <- d2[d2$locationName == unique_locations[i], ]
  Site_colors <- colors[as.character(subset_data$Site)]
  
  # Define line type based on cam_loc
  lty <- ifelse("ground" %in% subset_data$cam_loc, "solid", "dashed")
  
  segments(subset_data$dep_start, i, subset_data$dep_end, i, col = Site_colors, lwd = 3, lty = lty)
}

axis(2, at = 1:length(unique_locations), padj = 0, cex.axis = 0.7, labels = FALSE, tck = -0.01)
text(y = 1:length(unique_locations), par("usr")[1], labels = unique_locations, pos = 2, xpd = TRUE, cex = 0.5)
legend(max(d2$dep_end) + 15, length(unique_locations), c("Coiba", "Jicarón", "BCI"), pch = 15, col = colors, box.col = NA, cex = 1)

# Set the x-axis labels for each month
axis(1, at = seq(x_axis_start, x_axis_end, by = "1 month"), labels = format(seq(x_axis_start, x_axis_end, by = "1 month"), "%b"))
#graphics.off()



#Survey 3 (Figure S6)
d3 <-read.csv("Supp.DATA_tool_users.csv",
              header = T) ##Survey 3 - Capuchin tool use project
str(d3)

# Convert dep_start and dep_end to POSIXct
d3$dep_start <- as.POSIXct(d3$dep_start)
d3$dep_end <- as.POSIXct(d3$dep_end)
d3$cam_id <- as.integer(as.factor(d3$locationName))

# Convert the date columns to Date objects
d3$dep_start <- as.Date(d3$dep_start)
d3$dep_end <- as.Date(d3$dep_end)

#png(file= "D:/MPI/Biotropica_Supp/Figures/24.10.2023_CameraOP_Survey3.png",
   # width = 6, height = 6, units = "in", res = 600, pointsize = 10)
# Create a color palette for tool_site
colors <- c("1" = "tomato3", "0" = "blue")
# Specify the order of unique locations
unique_locations <- unique(d3$locationName)
unique_locations <- unique_locations[order(match(unique_locations, unique_locations[duplicated(unique_locations)]), unique_locations)]
# Increase left margin to allow for longer y-axis labels
par(mar = c(5, 8, 4, 2) + 0.1)
# Plot setup
plot(d3$dep_start, as.numeric(d3$locationName), type = "n", xlim = range(d3$dep_start, d3$dep_end),
     ylim = c(1, length(unique(d3$locationName))), xlab = "Year", ylab = "", yaxt = "n", main = "Survey 3")

# Loop through unique location names
for (i in 1:length(unique_locations)) {
  subset_data <- d3[d3$locationName == unique_locations[i], ]
  tool_site_colors <- colors[as.character(subset_data$tool_site)]
  segments(subset_data$dep_start, i, subset_data$dep_end, i, col = tool_site_colors, lwd = 3)
}

axis(2, at = 1:length(unique(d3$cam_id)), padj = 0, cex.axis = 0.7, labels = FALSE, tck = -0.01)
text(y = 1:length(unique(d3$cam_id)), par("usr")[1], labels = unique_locations, pos = 2, xpd = TRUE, cex = 0.5)

legend(max(d3$dep_start)-480,60, c("Tool-use area" , "Non-tool-use area"), pch=15, col=colors, box.col=NA, cex=1 )
#graphics.off()



