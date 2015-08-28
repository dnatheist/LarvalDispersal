#This script is to pre-calculate all the items needed to populate the database, or is required for subsequent chapters.

# This file, needs to be run after go.R and before each various chapter .Rmd files. Most vectors and dataframe additions will be commented out as they have been calculated and transferred to the database

# Some are variables eg: dispersalVelocity

# Some are formula eg: aging.

# Some are objects: data frame, vectors, matrices, dendrograms

# Many will be comments in an R file as they will only be run once and transfered. The file could be ‘spin’ rather than knit and much would become part of subsequent chapter’s methods sections. The naming conventions would need to be good and the calculations would need to be removed from the chapter .Rmd and R files.  This will eliminate duplication of code, and document when and where things are done. The data needs to be in the database so it can be used by Tableau if required.
#
#Some will be done up front. Some by chapter.It is done in this one file so as to ensure code is not duplicated in many places.

source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")# load code of A2R function
library(ggplot2)
library(ggdendro)
library(ape)
library(dendextend)
library(Hmisc)
library(ade4)

################################### Aging ############################################

#Aging has already been calculated in the database (qslAgedata). But for the record they are: 

#ageOL - tbc etc
#hatchDoY
#incubTime
#spawnDoY

larv<-merge(larv,qslAgeData,by.x="LarvalID", by.y="LarvaID")
row.names(larv)<-larv$Label

######################### Larval Dispersal Velocity ###################################
dispersalVelocity<-1200 # This important variable needs to be up-to-date with best estimate from the Iterated Mantle procedure. It is an output of Iterated Mantel.

####################### DNA, SNPs, Alleles and Base Pairs #############################
df<-CopyOfDMac14.1567snps
countUniqueSnps<-length(unique(df$X..1))-2 # number of unique snps in the Maccullochella set.
polymorphismPosition<-as.numeric(df$X..4) # the relative position of the polymorphism along the 69 base pair fragment of DNA.
fragmentLengthDNA<-median(as.numeric(nchar(as.character(df$X..2)))) #The most common fragment length with DaRT SNPS (69 in 2013).

#regex etc to pulls a,c,t,g for counts.
dfBPchange<-df$X..3
dfBPchange<-dfBPchange[-c(1:5)]
head(dfBPchange)
dfBPchange<-as.data.frame(dfBPchange)
dfBPchange<-gsub(".*:","",dfBPchange$dfBPchange) #remove all position numbers before the ':'
#dfBPchange<-gsub(" ","",dfBPchange$dfBPchange) #remove all spaces
dfBPchange<-as.data.frame(dfBPchange)
dfBPchange<-dfBPchange[!(dfBPchange$dfBPchange==""), ]
dfBPchange<-as.data.frame(dfBPchange)
table(dfBPchange)


maccullochellaDistMat <- dist(allsnps) #Create a distance matrix for all Maccullochella larvae

peeliiDistMat<-MCsnps[-c(1:7),]
peeliiDistMat <- dist(peeliiDistMat)
peeliiCluster <- hclust(peeliiDistMat) # Heirarchical Cluster, Murray cod only.

################################### Clades and Dendrograms ########################################
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


############################ Trout Cod Larvae Appendix ###############################

macquariensisDistMat <- dist(TCsnps) #Create Trout cod Distance Matrix
macquariensisCluster <- hclust(macquariensisDistMat) #create cluster

######################################################################################

################## Dispersal Distances and Nest Location Variables ###################
larv<-merge(larv,qslGeneticsForNestChapter,by.x="LarvalID", by.y="LarvaID", all.x=TRUE)

larv$nestDist<-larv$Distance.to.Angle.Crossing..m.-(dispersalVelocity*(larv$Day.of.Year-(larv$hatchDoY+larv$incubTime))) # To create a distance using the previously calculated mean drift velocity (m/d available since leaving brood care)


# for(id in 1:nrow(larv)){
#         larv$nestdist[larv$Label %in% withAnySibsDF$withAnySibs[id]] <- 0
# }

# The following is to Create sibling Corrected Nest Dist ie: this is a mean of all the siblings nest distances - it becomes the nestDist. Not the individually calculated one from age/dispersalVelocity.

##Calculate a family group
families<-as.data.frame(cutree(peeliiDendro, h = 26.5))
row.names(larv)<-larv$Label
larv<-merge(larv,families, by="row.names", all.x=TRUE)

##Create a data frame from which to develop nestID names (not yet combining siblings)
larv$roundedHatchDoY<-round(larv$hatchDoY)
larv$nestSeg<-round(larv$nestDist/100)
larv$siteLetter<-substr(larv$Label, 1, 1)
names(larv)[names(larv)=="cutree(peeliiDendro, h = 26.5)"] <- "family"

famNo<-unique(larv$family)
famNo<-na.omit(famNo)

##Initialise a data frame for use in for loop
final <- data.frame(Date=as.Date(character()),
                 File=character(), 
                 User=character(), 
                 stringsAsFactors=FALSE) 

##Loop through calculating means (corrected nestDist)
for (nn in famNo){
        newsub<-subset(larv,family==nn, select=c(LarvalID,Label,nestDist,hatchDoY,incubTime))
        newsub$sibCorrectedNestDist<-paste(mean(as.numeric(newsub$nestDist)))
        final<-rbind(final,newsub)
        }

#Now to calculate the furtherst potential settling distance if they moved at average for more of their remaining drift days.
larv$potentialSettlingDistance<-larv$Distance.to.Angle.Crossing..m.+(20-larv$ageOL)*1200

######################################################################################

#remove larvae that do not have genetic analysis done.
#Creat a MCsnps set with row names as a column.
MCchecklist<-row.names(MCsnps)
MCchecklist<-as.data.frame(MCchecklist)# 93 records

#remove a few more anomolies
MCchecklist1 <- as.data.frame(MCchecklist[-c(1:7), ])
# Keep every record in larv that is also in MCchecklist (i.e., the intersection).

larv_intersection <- larv[larv$Label %in% MCchecklist$MCchecklist,]
#Thanks: https://heuristically.wordpress.com/2009/10/08/delete-rows-from-r-data-frame/

larv<-larv_intersection
larv_intersection<-NULL

######################################################################################

# Calculate Mantel based on Sibling Corrected Distance

geodist<-data.frame(final$Label,final$sibCorrectedNestDist)
row.names(geodist)<-geodist$final.Label
geodist$final.Label<-NULL
geodist<-na.omit(geodist)
geodist1000<-geodist #save this estimate for haplogroups distance plot (after the Iterated Mantel has changed it)

GeoDistMat<-dist(geodist)
GeoDistMathm <- as.matrix(GeoDistMat)
heatmap(GeoDistMathm)
geoclust<-hclust(GeoDistMat)
plot(geoclust)

#make sure both matrices are in correct order - rows and cols
#First sort MCdm

peeliiDistMat<-as.data.frame(peeliiDistMat)
peeliiDistMat$sort<-row.names(peeliiDistMat)
peeliiDistMat <- peeliiDistMat[order(peeliiDistMat$sort),]#sort row order
peeliiDistMat$sort<-NULL
peeliiDistMat<-peeliiDistMat[,order(names(peeliiDistMat))]#sort column order
peeliiDistMat<-as.matrix(peeliiDistMat)

#Second sort GeoDist
GeoDistMathm<-as.data.frame(GeoDistMathm)
GeoDistMathm$sort<-row.names(GeoDistMathm)
GeoDistMathm <- GeoDistMathm[order(GeoDistMathm$sort),]#sort row order
GeoDistMathm$sort<-NULL
GeoDistMathm<-GeoDistMathm[,order(names(GeoDistMathm))]#sort column order
GeoDistMathm<-as.matrix(GeoDistMathm)

heatmap(GeoDistMathm)
geoclust<-hclust(GeoDistMat)
plot(geoclust)

#hist(larv$nestdist)
#larv1<-larv#save this estimate for haplogroups distance plot (after the Iterated Mantel has changed it)

#Now that these various matrices, class 'dist' objects are created we can proceed for plot.

###Plot and Correlate genetic and geographic distance matrices

#First a regression model is calculated and then the plot:

#Linear Regression Model
reg<-lm(peeliiDistMat[lower.tri(peeliiDistMat)]~GeoDistMathm[lower.tri(GeoDistMathm)])
summary(reg)
plot(GeoDistMathm[lower.tri(GeoDistMathm)],peeliiDistMat[lower.tri(peeliiDistMat)])
abline(reg)
# Correlations with significance levels
rcorr(GeoDistMathm[lower.tri(GeoDistMathm)],peeliiDistMat[lower.tri(peeliiDistMat)])#(x, type="pearson") # type can be pearson or spearman

###Mantel Test
mant<-mantel.rtest(as.dist(GeoDistMathm), as.dist(peeliiDistMat), nrepet = 9999)

##This results in a slightly better r value and Monty Carlo test so it seems appropriate that sibling corrected nest distances are used in the iterated mantel rather than the uncorrected ones. That said the difference, if only, is expected to only be slight.

########################################################################################
# Autocorrelation Test? How to interpret this?

acf(as.dist(GeoDistMathm))
acf(as.dist(peeliiDistMat))

###########################################################################################