###Script for Activity patterns of white-faced capuchins and agoutis in Coiba, Jicarón, and Barro Colorado
###Contrasting arboreal capuchin activity with agouti activity
###Terrestrial capuchin activity with agouti activity (excluding tool sites)
rm(list = ls()) #this cleans the environment work space
shell("cls") #this cleans the console

library(overlap)
library(activity)

###Read in data# 
#d.summed object created on S1-Script-Data-Prep
d.summed <- read.csv(file = "Processed_Data.csv", header = T, row.names = 1)

###Overlap coefficients#####

#Mainland capuchin (arboreal activity) vs agouti
agouti_capuchinBCI <- overlapEst(subset(d.summed, scientificName=="Dasyprocta punctata" & sites == "Mainland")$suntime,
                                 subset(d.summed, scientificName=="Cebus capucinus" & sites == "Mainland" & cam_loc == "canopy")$suntime,
                                 type="Dhat4") #0.8893958 
agouti_capuchinBCI_boot <- bootstrap(subset(d.summed, scientificName=="Dasyprocta punctata" & sites == "Mainland")$suntime,
                                     subset(d.summed, scientificName=="Cebus capucinus" & sites == "Mainland" & cam_loc == "canopy")$suntime,
                                     type="Dhat4", 10000)
bootCI(agouti_capuchinBCI, agouti_capuchinBCI_boot) #0.8146501 0.9520009
bootCI(agouti_capuchinBCI, agouti_capuchinBCI_boot, conf = 0.84) #0.8373690 0.9367403
#Island capuchin (arboreal activity) vs agouti
agouti_capuchin_island <- overlapEst(subset(d.summed, scientificName=="Dasyprocta punctata" & sites == "Island")$suntime,
                                     subset(d.summed, scientificName=="Cebus capucinus" & sites == "Island" & cam_loc == "canopy")$suntime,
                                     type="Dhat4") #0.7280255 
agouti_capuchin_island_boot <- bootstrap(subset(d.summed, scientificName=="Dasyprocta punctata" & sites == "Island")$suntime,
                                         subset(d.summed, scientificName=="Cebus capucinus" & sites == "Island" & cam_loc == "canopy")$suntime,
                                         type="Dhat4", 10000)
bootCI(agouti_capuchin_island, agouti_capuchin_island_boot) #0.6043022 0.8372302
bootCI(agouti_capuchin_island, agouti_capuchin_island_boot, conf = 0.84) #0.6421716 0.8075390




#Island capuchin (terrestrial activity minus tool sites) vs agouti
agouti_capuchin_island_ground <- overlapEst(subset(d.summed, scientificName=="Dasyprocta punctata" & sites == "Island" & tool_site == "0")$suntime,
                                            subset(d.summed, scientificName=="Cebus capucinus" & sites == "Island" & cam_loc == "ground" & tool_site == "0")$suntime,
                                            type="Dhat4") #0.7848011 
agouti_capuchin_island_ground_boot <- bootstrap(subset(d.summed, scientificName=="Dasyprocta punctata" & sites == "Island" & tool_site == "0")$suntime,
                                                subset(d.summed, scientificName=="Cebus capucinus" & sites == "Island" & cam_loc == "ground" & tool_site == "0")$suntime,
                                                type="Dhat4", 10000)
bootCI(agouti_capuchin_island_ground, agouti_capuchin_island_ground_boot) #0.7621876 0.8076651
bootCI(agouti_capuchin_island_ground, agouti_capuchin_island_ground_boot, conf = 0.84) #0.7682701 0.8016023

###Overlap plots####
png(file= "FIGURE2.png",
width = 9, height = 18, units = "in", res = 500, pointsize = 10)
par(mfrow= c(3,1),
    mar = c(5, 5, 4, 2)) #Creates FIGURE 2 

overlapPlot(subset(d.summed, scientificName=="Dasyprocta punctata" & sites == "Mainland")$suntime,
            subset(d.summed, scientificName=="Cebus capucinus" & sites == "Mainland" & cam_loc == "canopy")$suntime,
            bty = "n" , xaxt='n', olapcol = "white", ylim = c(0,0.15), 
            adjust = 2, rug = T, linewidth = c(2,2), linecol= c('tomato3', 'gray20'),
            ylab= "Density of activity", xlab = "", main= "Agouti - capuchin (canopy/BCI)",
            cex.axis = 2, cex.main=2, cex.lab=2)
axis(1, at=c(0, 6, 12, 18, 24),
     labels=c("midnight", "sunrise", "noon", "sunset", "midnight"), cex.axis = 2)
text(1,0.14,substitute(paste(bold("(a)"))),cex=3)
legend('topright', c("Agouti", "Capuchin"), lty=c(1,2), 
       col= c('tomato3', 'gray20'), bty='n', lwd = c(2,2),
       cex = 2)

overlapPlot(subset(d.summed, scientificName=="Dasyprocta punctata" & sites == "Island")$suntime,
            subset(d.summed, scientificName=="Cebus capucinus" & sites == "Island" & cam_loc == "canopy")$suntime,
            bty = "n" , xaxt='n', olapcol = "white", ylim = c(0,0.15), 
            adjust = 2, rug = T, linewidth = c(2,2), linecol= c('tomato3', 'gray20'),
            ylab= "Density of activity", xlab = "", main= "Agouti - capuchin (canopy/Island)",
            cex.axis = 2, cex.main=2, cex.lab=2)
axis(1, at=c(0, 6, 12, 18, 24),
     labels=c("midnight", "sunrise", "noon", "sunset", "midnight"),  cex.axis = 2)
text(1,0.14,substitute(paste(bold("(b)"))),cex=3)
legend('topright', c("Agouti", "Capuchin"), lty=c(1,2), 
       col= c('tomato3', 'gray20'), bty='n', lwd = c(2,2),
       cex = 2)

overlapPlot(subset(d.summed, scientificName=="Dasyprocta punctata" & sites == "Island" & tool_site == "0")$suntime,
            subset(d.summed, scientificName=="Cebus capucinus" & sites == "Island" & cam_loc == "ground" & tool_site == "0")$suntime,
            bty = "n" , xaxt='n', olapcol = "white", ylim = c(0,0.15), 
            adjust = 2, rug = T, linewidth = c(2,2), linecol= c('tomato3', 'gray20'),
            ylab= "Density of activity", xlab = "", main= "Agouti - capuchin (ground/Island)",
            cex.axis = 2, cex.main=2, cex.lab=2)
axis(1, at=c(0, 6, 12, 18, 24),
     labels=c("midnight", "sunrise", "noon", "sunset", "midnight"),  cex.axis = 2)
text(1,0.14,substitute(paste(bold("(c)"))),cex=3)
legend('topright', c("Agouti", "Capuchin"), lty=c(1,2), 
       col= c('tomato3', 'gray20'), bty='n', lwd = c(2,2),
       cex = 2)
graphics.off()

##For Supplemental Information####
#Agouti and capuchin (ground) overlap at tool-use area only
agouti_capuchin_island_toolsites <- overlapEst(subset(d.summed, scientificName=="Dasyprocta punctata" & sites == "Island" & tool_site == 1)$suntime,
                                               subset(d.summed, scientificName=="Cebus capucinus" & sites == "Island" & cam_loc == "ground" & tool_site == 1)$suntime,
                                               type="Dhat4") #0.739627  
agouti_capuchin_island_toolsites_boot <- bootstrap(subset(d.summed, scientificName=="Dasyprocta punctata" & sites == "Island" & tool_site == 1)$suntime,
                                                   subset(d.summed, scientificName=="Cebus capucinus" & sites == "Island" & cam_loc == "ground" & tool_site == 1)$suntime,
                                                   type="Dhat4", 10000)
bootCI(agouti_capuchin_island_toolsites, agouti_capuchin_island_toolsites_boot) #0.6043022 0.8372302
bootCI(agouti_capuchin_island_toolsites, agouti_capuchin_island_toolsites_boot, conf = 0.84) #0.6421716 0.8075390

###FIGURE S10
par(mar = c(5, 5, 4, 2))

overlapPlot(subset(d.summed, scientificName=="Dasyprocta punctata" & sites == "Island" & tool_site == "1")$suntime,
            subset(d.summed, scientificName=="Cebus capucinus" & sites == "Island" & cam_loc == "ground" & tool_site == "1")$suntime,
            bty = "n" , xaxt='n', olapcol = "white", ylim = c(0,0.15), 
            adjust = 2, rug = T, linewidth = c(2,2), linecol= c('tomato3', 'gray20'),
            ylab= "Density of activity", xlab = "", main= "Agouti - capuchin (tool-sites)",
            cex.axis = 2, cex.main=2, cex.lab=2)
axis(1, at=c(0, 6, 12, 18, 24),
     labels=c("midnight", "sunrise", "noon", "sunset", "midnight"),  cex.axis = 2)
legend('topright', c("Agouti", "Capuchin"), lty=c(1,2), 
       col= c('tomato3', 'gray20'), bty='n', lwd = c(2,2),
       cex = 2)
graphics.off()
