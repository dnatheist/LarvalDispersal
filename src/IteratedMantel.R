## @knitr IM

library(Hmisc)
library(ade4)

#This file is to calculate the data as before but mainly to iterate the mantel test from 1 to 5000 metres at 100m increments (to save time).

# The calculations that can happen outside the iteration are:
#  
  #Age from Otolith Length: 74.308*[MeanOtolithLength]-4.44361
  #Hatch DoY :  [Day of Year Caught]-[Age From Otolith Length]
  #Incubation: 20.67-0.667*[WaterTemp(DegC) Mean]
  #Spawnin:[Hatch]-[Incubation]
  #larv$nestdist<-larv$Distance.to.Angle.Crossing..m.-(300*(larv$Day.of.Year-(larv$hatchdoy+7)))
#
#

larv$ageOL<-74.308*larv$Mean.Otolith.Length.is.in.Millimetres.for.comparison.with.Adults-4.44361
larv$hatchdoy<-larv$Day.of.Year-larv$ageOL
larv$incTime<-20.67-0.667*larv$WaterTemp.DegC..Mean
larv$spawn<-larv$hatchdoy-larv$incTime
#larv$nestdist<-larv$Distance.to.Angle.Crossing..m.-(0*(larv$Day.of.Year-(larv$hatchdoy+7)
                                                       
#Create a MCsnps set with row names as a column.
MCchecklist<-row.names(MCsnps)
MCchecklist<-as.data.frame(MCchecklist)# 93 records

#remove a few more anomolies
MCchecklist1 <- as.data.frame(MCchecklist[-c(1:7), ])
# Keep every record in larv that is also in MCchecklist (i.e., the intersection).

larv_intersection <- larv[larv$Label %in% MCchecklist$MCchecklist,]
#Thanks: https://heuristically.wordpress.com/2009/10/08/delete-rows-from-r-data-frame/

larv<-larv_intersection
larv_intersection<-NULL

itmant <- matrix(nrow=5000, ncol=3) #Is 5000 DF to store result but NA omited later. They result from the increment 100 in the for loop below.

#Iteration begins here:
for (nd in seq(1,5000, by=100)){#To be 0:5000 eventually for(i in seq(1, 10, by = 2)) 

larv$nestdist<-larv$Distance.to.Angle.Crossing..m.-(nd*(larv$Day.of.Year-(larv$hatchdoy+7)))  

###########
# Create GenDist from code in the Murray Cod SNPS table
MCdm<-MCsnps[-c(1:7),] #remove non-numeric variables
MCdm <- dist(MCdm) # Create a Murray Cod distance matrix
MCdm<-as.matrix(MCdm)
MCdm<-as.data.frame(MCdm)
#This is to be used for plotting
###########
#Create Geographic Distance Matrix using Nest Distance
geodist<-data.frame(larv$Label,larv$nestdist)
row.names(geodist)<-geodist[,1]
geodist$larv.Label<-NULL
geodist<-na.omit(geodist)
#geodist<-geodist[complete.cases(geodist),]

GeoDistMat<-dist(geodist)
GeoDistMathm <- as.matrix(GeoDistMat)


#make sure both matrices are in correct order - rows and cols
#First sort MCdm

MCdm<-as.data.frame(MCdm)
MCdm$sort<-row.names(MCdm)
MCdm <- MCdm[order(MCdm$sort),]#sort row order
MCdm$sort<-NULL
MCdm<-MCdm[,order(names(MCdm))]#sort column order
MCdm<-as.matrix(MCdm)

#Second sort GeoDist
GeoDistMathm<-as.data.frame(GeoDistMathm)
GeoDistMathm$sort<-row.names(GeoDistMathm)
GeoDistMathm <- GeoDistMathm[order(GeoDistMathm$sort),]#sort row order
GeoDistMathm$sort<-NULL
GeoDistMathm<-GeoDistMathm[,order(names(GeoDistMathm))]#sort column order
GeoDistMathm<-as.matrix(GeoDistMathm)

mant<-mantel.rtest(as.dist(GeoDistMathm), as.dist(MCdm), nrepet = 9999)
#print(nd)
#print(mant$obs)
#print(mant$pvalue)

 itmant[nd,] <- c(nd, mant$obs, mant$pvalue)
  }

itmant<-na.omit(itmant)
itmantdf<-as.data.frame(itmant)
plot(itmantdf$V1,itmantdf$V2)


