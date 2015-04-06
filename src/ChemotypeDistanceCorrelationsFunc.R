## @knitr chemotypeDistanceMatrices
library(clusterSim)
#This creates distance matrices for chemotype and geographic distance and then use a mantel test to see correlation
#First using Otolith Core chemistry.

allVars<-ChemAnalCore

allVars<-allVars[c(2,84,13,14,15,120:151)] #
allVars<-allVars[complete.cases(allVars),] #remove any nulls
row.names(allVars)<-allVars[,1]
allVars$Label<-NULL
colnames(allVars)[1]<-"SiteName"

allVars$SiteName<-as.factor(allVars$SiteName)

df<-allVars
chemoTypeDistMat <- dist(df)
hClusters <- hclust(chemoTypeDistMat)
plot(hClusters,labels=(df$SiteName), hang = -1, main="Otolith Core Chemistry")

#Make distance matrices for geographic distance as well
allVars<-ChemAnalCore
allVars<-allVars[c(2,84,13,14,15,107,120:151)]
allVars<-allVars[complete.cases(allVars),] #remove any nulls
allVars<-allVars[c(1,6)] # Distance from Angle Crossing changes from 107 above to 6
row.names(allVars)<-allVars[,1]
allVars$Label<-NULL

geoDist<-allVars
geoDist<-na.omit(geoDist)
#geoDistColl1000<-geodist #save this estimate for haplogroups distance plot (after the Iterated Mantel has changed it)

geoDistMat<-dist(geoDist)


#make sure both matrices are in correct order - rows and cols
#Check all is in order
as.matrix(geoDistMat)[1:5, 1:5] # zero distances in the first 5
as.matrix(chemoTypeDistMat)[1:5, 1:5]

#Conduct Mantel Test on Matrices
mant<-mantel.rtest(as.dist(geoDistMat), as.dist(chemoTypeDistMat), nrepet = 9999)
mant
plot(geoDistMat,chemoTypeDistMat, main="Otolith Core Chemistry")

#Even if the data is scaled there is no correlation to speak of.
#########################################
#allVars$SiteName<-NULL
#allVars <- droplevels(allVars)#Not sure why get error without this line
#allVarScaled<-data.Normalization (allVars,type="n7",normalization="column")
# allVarScaled$Label<-row.names(allVarScaled)
# sChemAnalCore<-subset(ChemAnalCore, select=c(Label,SiteName))
# allVars<-merge(allVarScaled, sChemAnalCore,by = "Label")
# allVars$Label<-NULL
# df<-allVarScaled
# distxy <- dist(df)
# hClusters <- hclust(distxy)
# plot(hClusters,labels=(df$Site.Name), hang = -1)
# df<-allVarScaled
###################################

#But now with important variables only. These include Delta13C, Delta15N, CN ratio, B, K, V, Na, Rb and were settled on by using Random Forest.

otoParts<-list(ChemAnalCore,ChemAnalE1,ChemAnalE2,ChemAnal0,ChemAnalAll4,ChemAnalNonCore)
lapply(otoParts, otoPartChem)

# otoPartChem <- function(otoParts){
#   
#   ImportantVars<-otoParts
#   
#   ImportantVars<-ImportantVars[c(2,84,13,14,15,121,122,127,131,139)] #
#   ImportantVars<-ImportantVars[complete.cases(ImportantVars),] #remove any nulls
#   #row.names(ImportantVars)<-ImportantVars[,1]
#   ImportantVars$Label<-NULL
#   colnames(ImportantVars)[1]<-"SiteName"
#   
#   ImportantVars$SiteName<-as.factor(ImportantVars$SiteName)
#   
#   df<-ImportantVars
#   print(df)
#   chemoTypeDistMat <- dist(df)
#   hClusters <- hclust(chemoTypeDistMat)
#   plot(hClusters,labels=(df$SiteName), hang = -1, main=deparse(substitute(otoParts)))
#   
#   #Make distance matrices for geographic distance as well
#   ImportantVars<-otoParts
#   ImportantVars<-ImportantVars[c(2,84,13,14,15,107,121,122,127,131,139)]
#   ImportantVars<-ImportantVars[complete.cases(ImportantVars),] #remove any nulls
#   ImportantVars<-ImportantVars[c(1,6)] # Distance from Angle Crossing changes from 107 above to 6
#   #row.names(ImportantVars)<-ImportantVars[,1]
#   ImportantVars$Label<-NULL
#   
#   geoDist<-ImportantVars
#   geoDist<-na.omit(geoDist)
#   #geoDistColl1000<-geodist #save this estimate for haplogroups distance plot (after the Iterated Mantel has changed it)
#   
#   geoDistMat<-dist(geoDist)
#   
#   #make sure both matrices are in correct order - rows and cols
#   #Check all is in order
#   print(as.matrix(geoDistMat)[1:5, 1:5]) # zero distances in the first 5
#   print(as.matrix(chemoTypeDistMat)[1:5, 1:5])
#   
#   #Conduct Mantel Test on Matrices
#   mantIV<-mantel.rtest(as.dist(geoDistMat), as.dist(chemoTypeDistMat), nrepet = 9999)
#   print(mantIV)
#   plot(geoDistMat,chemoTypeDistMat, main="Otolith Core Chemistry")
# }

##############################

#But now with important variables only AFTER SCALING. These include C13, N15, CN ratio, B, K, V, Na, Rb.
ImportantVars<-ChemAnalCore

ImportantVars<-ImportantVars[c(2,84,13,14,15,121,122,127,131,139)] #
ImportantVars<-ImportantVars[complete.cases(ImportantVars),] #remove any nulls
row.names(ImportantVars)<-ImportantVars[,1]
ImportantVars$Label<-NULL
colnames(ImportantVars)[1]<-"SiteName"

ImportantVars$SiteName<-as.factor(ImportantVars$SiteName)

ImportantVars$SiteName<-NULL
ImportantVars <- droplevels(ImportantVars)#Not sure why get error without this line
ImpVarScaled<-data.Normalization (ImportantVars,type="n4",normalization="column")
ImpVarScaled$Label<-row.names(ImpVarScaled)
sChemAnalCore<-subset(ChemAnalCore, select=c(Label,SiteName))
ImportantVars<-merge(ImpVarScaled, sChemAnalCore,by = "Label")
ImportantVars$Label<-NULL
df<-ImpVarScaled
distxy <- dist(df)
hClusters <- hclust(distxy)
plot(hClusters,labels=(df$Site.Name), hang = -1, main="Otolith Core Chemistry")
df<-ImpVarScaled

#df<-ImportantVars
chemoTypeDistMat <- dist(df)
hClusters <- hclust(chemoTypeDistMat)
plot(hClusters,labels=(df$SiteName), hang = -1, main="Otolith Core Chemistry")

#Make distance matrices for geographic distance as well
ImportantVars<-ChemAnalCore
ImportantVars<-ImportantVars[c(2,84,13,14,15,107,121,122,127,131,139)]
ImportantVars<-ImportantVars[complete.cases(ImportantVars),] #remove any nulls
ImportantVars<-ImportantVars[c(1,6)] # Distance from Angle Crossing changes from 107 above to 6
row.names(ImportantVars)<-ImportantVars[,1]
ImportantVars$Label<-NULL

geoDist<-ImportantVars
geoDist<-na.omit(geoDist)
#geoDistColl1000<-geodist #save this estimate for haplogroups distance plot (after the Iterated Mantel has changed it)

geoDistMat<-dist(geoDist)

#make sure both matrices are in correct order - rows and cols
#Check all is in order
as.matrix(geoDistMat)[1:5, 1:5] # zero distances in the first 5
as.matrix(chemoTypeDistMat)[1:5, 1:5]

#Conduct Mantel Test on Matrices
mantIVScaled<-mantel.rtest(as.dist(geoDistMat), as.dist(chemoTypeDistMat), nrepet = 9999)
mantIVScaled
plot(geoDistMat,chemoTypeDistMat, main="Otolith Core Chemistry")
