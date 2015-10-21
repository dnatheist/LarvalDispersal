#I think this needs re vamping with dplyr. It is too dependant on column position. 

otoPartChem <- function(otoParts){
  
  ImportantVars<-otoParts
  
  ImportantVars<-ImportantVars[c(2,84,13,14,15,121,122,127,131,139)] #
  ImportantVars<-ImportantVars[complete.cases(ImportantVars),] #remove any nulls
  #row.names(ImportantVars)<-ImportantVars[,1]
  ImportantVars$Label<-NULL
  colnames(ImportantVars)[1]<-"SiteName"
  
  ImportantVars$SiteName<-as.factor(ImportantVars$SiteName)
  
  df<-ImportantVars
  #print(df)
  chemoTypeDistMat <- dist(df)
  hClusters <- hclust(chemoTypeDistMat)
  plot(hClusters,labels=(df$SiteName), hang = -1, main=deparse(substitute(otoParts)))
  
  #Make distance matrices for geographic distance as well
  ImportantVars<-otoParts
  ImportantVars<-ImportantVars[c(2,84,13,14,15,107,121,122,127,131,139)]
  ImportantVars<-ImportantVars[complete.cases(ImportantVars),] #remove any nulls
  ImportantVars<-ImportantVars[c(1,6)] # Distance from Angle Crossing changes from 107 above to 6
  #row.names(ImportantVars)<-ImportantVars[,1]
  ImportantVars$Label<-NULL
  
  geoDist<-ImportantVars
  geoDist<-na.omit(geoDist)
  #geoDistColl1000<-geodist #save this estimate for haplogroups distance plot (after the Iterated Mantel has changed it)
  
  geoDistMat<-dist(geoDist)
  
  #make sure both matrices are in correct order - rows and cols
  #Check all is in order
  print(as.matrix(geoDistMat)[1:5, 1:5]) # zero distances in the first 5
  print(as.matrix(chemoTypeDistMat)[1:5, 1:5])
  
  #Conduct Mantel Test on Matrices
  mantIV<-mantel.rtest(as.dist(geoDistMat), as.dist(chemoTypeDistMat), nrepet = 9999)
  print(mantIV)
  plot(geoDistMat,chemoTypeDistMat, main="Otolith Core Chemistry")
}