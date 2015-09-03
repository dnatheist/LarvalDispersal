library(Hmisc)
library(ade4)

#This file is to use the mantel test to look for correlation between genetic distance and geographic distance. In this case by iterating the geographic distance from -2000 to 5000 metres at 10m increments, to find the mean distance that provides the best correlation between the two distance matrices. Given Isolation by Distance, the optimal correlation provides an estimate of the mean dispersal velocity.

# Firstly remove larvae from larv that do not have genetic analysis done.
# Create a MCsnps set with row names as a column.
MCchecklist<-row.names(MCsnps)
MCchecklist<-as.data.frame(MCchecklist)# 93 records

#remove a few more anomolies. These first seven rows were test larvae that are not for inclusion.
MCchecklist1 <- as.data.frame(MCchecklist[-c(1:7), ])

# Keep every record in larv that is also in MCchecklist (i.e., the intersection).
larv_intersection <- larv[larv$Label %in% MCchecklist$MCchecklist,]
#Thanks: https://heuristically.wordpress.com/2009/10/08/delete-rows-from-r-data-frame/

larv<-larv_intersection
rm(larv_intersection)

# Create GenDist from code in the Murray Cod SNPS table
MCdm<-MCsnps[-c(1:7),] #remove non-numeric variables
MCdm <- dist(MCdm) # Create a Murray Cod distance matrix
MCdm<-as.matrix(MCdm)
MCdm<-as.data.frame(MCdm)
#This is to be used for plotting. It does not vary inside the iteration (Genetic distance is fixed) while actual geographic distance (dispersal velocity) is unknown.

#First sort MCdm - just to be sure that the two matrices are ordered identically before mantel test.
MCdm<-as.data.frame(MCdm)
MCdm$sort<-row.names(MCdm)
MCdm <- MCdm[order(MCdm$sort),]#sort row order
MCdm$sort<-NULL
MCdm<-MCdm[,order(names(MCdm))]#sort column order
MCdm<-as.matrix(MCdm)

itmant <- matrix(nrow=21, ncol=3) #Is a 7002 row DF to store result but NA generated in for loop below are omited later. 

#Iteration begins here:
ptm <- proc.time()
for (broodCareDuration in seq(2,12, by=0.5)){# Possible brood care duration from 2 to 12 days, 0.5 day increments (Hump says 7?)
        
        larv$tmpDispDur<-larv$ageOL-broodCareDuration #begining an if statement. not working yet. May not be do-able?
        #if bc+dispersal duration > age then do it
        #else stop
        
        larv$tmpNestDist<-larv$Distance.to.Angle.Crossing..m.-(1320*(larv$Day.of.Year-(larv$hatchDoY+broodCareDuration))) 
        
        #Create Geographic Distance Matrix using Nest Distance
        geodist<-data.frame(larv$Label,larv$tmpNestDist)
        row.names(geodist)<-geodist[,1]
        geodist$larv.Label<-NULL
        geodist<-na.omit(geodist)
        GeoDistMat<-dist(geodist)
        GeoDistMathm <- as.matrix(GeoDistMat)
        
        #make sure GeoDist matrix is in correct order - rows and cols
        GeoDistMathm<-as.data.frame(GeoDistMathm)
        GeoDistMathm$sort<-row.names(GeoDistMathm)
        GeoDistMathm <- GeoDistMathm[order(GeoDistMathm$sort),]#sort row order
        GeoDistMathm$sort<-NULL
        GeoDistMathm<-GeoDistMathm[,order(names(GeoDistMathm))]#sort column order
        GeoDistMathm<-as.matrix(GeoDistMathm)
        
        mant<-mantel.rtest(as.dist(GeoDistMathm), as.dist(MCdm), nrepet = 999)
        
        itmant[broodCareDuration,] <- c(broodCareDuration, mant$obs, mant$pvalue)
}

larv$tmpNestDist<-NULL #tidy up temporary variables
rm(broodCareDuration)
rm(geodist)
rm(GeoDistMat)
rm(GeoDistMathm)
rm(mant)

proc.time() - ptm

itmant<-na.omit(itmant) #remove NAs 
itmantdf<-as.data.frame(itmant) #make into a data frame
plot(itmantdf$V1,itmantdf$V2,main = "Mantel Correlation Between Genetic Distance and Geographic Distance", xlab = "Geographic Distance (metres)", ylab = "Correlation with Genetic Distance Matrix (Mantel Test)")
rm(itmant)

###########
require(plotrix)
op <- par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
plot(itmantdf$V1, itmantdf$V2, col="black", pch=21, bg = "grey", cex = 2,
     xlim=c(0,200), ylim=c(0,.35), ylab="", xlab="", axes=F)
axis(1)
axis(2) 
reg1 <- lm(itmantdf$V2~itmantdf$V1)
ablineclip(reg1, lwd=2,x1 = .9, x2 = 1.2) 
par(las=0)
mtext("Daily Drift Distance (m)", side=1, line=2.5, cex=1.5)
mtext("Correlation", side=2, line=3.7, cex=1.5)

###########

broodCareDuration<- subset(itmantdf, V2==max(V2) , select = V1)


