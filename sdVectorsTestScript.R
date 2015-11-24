library('ProjectTemplate')
load.project()

source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")# load code of A2R function
library(ggplot2)
library(ggdendro)
library(ape)
library(dendextend)

################################### Clades and Dendrograms ########################################
maccullochellaDistMat <- dist(allsnps) #Create a distance matrix for all Maccullochella larvae

peeliiDistMat<-MCsnps[-c(1:7),]
peeliiDistMat <- dist(peeliiDistMat)
peeliiCluster <- hclust(peeliiDistMat) # Heirarchical Cluster, Murray cod only.


#A2Rplot(MChc, k = 12, boxes = FALSE, col.up = "gray50")
peeliiDendro <- color_branches(peeliiCluster, k = 12)
peeliiDendro <- color_labels(peeliiDendro, k = 12) #Use h or k

#ExtractRaceClades
raceDendro <- color_labels(peeliiCluster, k = 3)#Use h or k
raceCladeNumber<-cutree(raceDendro,k=3) #This is like using:dendextend:::cutree.dendrogram(dend1,h=70) h or k can be specified
temp<-as.data.frame(raceCladeNumber)
larv<-merge(larv,temp,all.x=TRUE,by="row.names")
larv$Row.names<-NULL

#Extract Family Clades

familyDendro<-color_labels(peeliiCluster,k=12) #This is like using:dendextend:::cutree.dendrogram(dend1,h=70) h or k can be specified
familyCladeNumber<-cutree(familyDendro,k=12) #This is like using:dendextend:::cutree.dendrogram(dend1,h=70) h or k can be specified
temp<-as.data.frame(familyCladeNumber)
row.names(larv)<-larv$Label
larv<-merge(larv,temp,all.x=TRUE,by="row.names")
larv$Row.names<-NULL

###################################Extract Siblings and those Without Siblings##########################

peeliiDistMat<-as.matrix(peeliiDistMat)
withAnySibs<-subset(peeliiDistMat, apply(peeliiDistMat, 1, function(peeliiDistMat){any(peeliiDistMat > 0 & peeliiDistMat < 26.5)})) #The siblings have a dissimilarity index of less than 26.5. This can be seen in the dendrogram and in the scatterplot for IBD.

withSibsLogical<-withAnySibs<26.5 #a logical DF to record T and F for sibling pairs
withAnySibs<-withSibsLogical*withAnySibs #This turns to 0 all the genetic similarity above 26.5 (no-sibs)

onlyWithSibs<-withAnySibs[!sapply(withAnySibs, function(x) all(x == 0))] #now remove those with no sibs
sibsLogical<-withAnySibs>0

#Make a Label list of those with sibs (data frame) with names 
withAnySibs<-row.names(withAnySibs) #make a list of those with sibs
withAnySibsDF<-as.data.frame(withAnySibs)
row.names(withAnySibsDF)<-withAnySibsDF[,1] # Make sensible row names (larva labels)
withAnySibsDF<-as.data.frame(withAnySibsDF) #make back to a vector

#Now a Label list of those without any sibs (check this - number seem right)
y<-row.names(peeliiDistMat)
withOutAnySibs<-as.data.frame(setdiff(y,withAnySibsDF$withAnySibs) )

######################################################################################

row.names(larv)<-larv$Label

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
#This is one part of the Mantel test given Isolation by Distance. It does not vary inside the iteration (Genetic distance is fixed) while actual geographic distance (dispersal velocity) is the unknown and is recalculated within the iteration below.

#First sort MCdm - just to be sure that the two matrices are ordered identically before mantel test.
MCdm<-as.data.frame(MCdm)
MCdm$sort<-row.names(MCdm)
MCdm <- MCdm[order(MCdm$sort),]#sort row order
MCdm$sort<-NULL
MCdm<-MCdm[,order(names(MCdm))]#sort column order
MCdmForSibs<-as.data.frame(MCdm) #need this matrix for sibling analysis .r
MCdm<-as.matrix(MCdm)

itmant <- matrix(nrow=7001, ncol=3) #Is a 7002 row DF to store result but NA generated in for loop below are omited later. 

#Initialise two vectors for variance analysis a la Shuvo suggestion
geoDistVector<-vector()
geoDistPerDayVector<-vector()
MCdmVector<-vector()
ldvRangeVector<-vector()
ldvRangeMatrix<-matrix(nrow = 86, ncol=86)
daysDispersalVector<-vector()
#Days available for dispersal
daysDispersal<-larv$Day.of.Year-(larv$hatchedDoY+7)


#Iteration begins here:
ptm <- proc.time()
for (ldvRange in seq(-2000,5000, by=10)){# Possible dispersal from 2km upstream to 5 km downstream per day, 10m increments 
        
        larv$tmpNestDist<-larv$Distance.to.Angle.Crossing..m.-(ldvRange*daysDispersal)
        
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
        
        #make Vectors for linear model
        MCdmVector<-append(MCdmVector,(as.vector(MCdm)))
        geoDistVector<-append(geoDistVector,(as.vector(GeoDistMathm)))
        ldvRangeVector<-append(ldvRangeVector,replicate(nrow(geodist),ldvRange))
        daysDispersalVector<-append(daysDispersalVector,replicate(nrow(geodist),daysDispersal))
        #geoDistPerDay<-geoDistVector/daysDispersal
        #geoDistPerDayVector<-append(geoDistPerDayVector,geoDistVector)
        #ldvRangeMatrix<-(replicate(length(ldvRange), ldvRange))
        
        #OLDappend(ldvRangeVector, (as.vector(ldvRange)))
        #OLDldvRangeMatrix<-(replicate(length(ldvRangeVector), ldvRangeVector))
        
        mant<-mantel.rtest(as.dist(GeoDistMathm), as.dist(MCdm), nrepet = 99)
        ndr<-ldvRange+2001 #to keep in positive for row numbers in itmant data frame
        itmant[ndr,] <- c(ldvRange, mant$obs, mant$pvalue)
}

head(geoDistVector, 100)
head(ldvRangeVector, 100)
head(daysDispersalVector,100)
ldvVector<-geoDistVector/daysDispersalVector
head(ldvVector,100)

ldvMat<-matrix(ldvVector, nrow=86, ncol=86)
lower.tri(ldvMat, diag = TRUE)
ldvMat[lower.tri(ldvMat,diag=TRUE)] <- NA
head(ldvMat, 3)
mean(ldvMat, na.rm = TRUE)


mod<-lm(as.vector(ldvMat)~1)
summary(mod)


#geoDistPerDayVector<-geoDistVector/daysDispersal
#larv$tmpNestDist<-NULL #tidy up temporary variables
#rm(ndr)
#rm(geodist)
#rm(GeoDistMat)
#rm(GeoDistMathm)
#rm(mant)

proc.time() - ptm

itmant<-na.omit(itmant) #remove NAs 
itmantdf<-as.data.frame(itmant) #make into a data frame
plot(itmantdf$V1,itmantdf$V2,main = "Mantel Correlation Between Pairwise Genetic Distance and Geographic Distance", xlab = "Geographic Distance (metres)", ylab = "Correlation with Genetic Distance Matrix (Mantel Test)",cex.main = 2.5,cex.lab=2)
rm(itmant)

###########
require(plotrix)
op <- par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
plot(itmantdf$V1, itmantdf$V2, col="black", pch=21, bg = "grey", cex = 2,
     xlim=c(0,5000), ylim=c(0,.25), ylab="", xlab="", axes=F)
axis(1)
axis(2) 
reg1 <- lm(itmantdf$V2~itmantdf$V1)
ablineclip(reg1, lwd=2,x1 = .9, x2 = 1.2) 
par(las=0)
mtext("Daily Drift Distance (m)", side=1, line=2.5, cex=1.5)
mtext("Correlation", side=2, line=3.7, cex=1.5)
###########

dispersalVelocity<- subset(itmantdf, V2==max(V2) , select = V1)
dispersalVelocity<-as.numeric(dispersalVelocity)
