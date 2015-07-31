## @knitr sumMatricesToFindNests

require(clusterSim) # to allow data normalization

# This code is to take the three distances matrices; (genetic distance, spatial distance, temporal distance) and add them and then cluster to give a better approximation of the nest an individual is likely to belong to. 
#However, on relfection I dont think this is valid. Not on dist matrices. should just be on three columns of values clade, nestdist, hatchdoy

#First need to create temporal distance from the larv data set so need label and hatch day of the year.

#Create Hatch DoY and then Temporal Distance Matrix

temporalDistance<-data.frame(larv2$Label,(larv2$Day.of.Year-(larv2$ageOL)))
colnames(temporalDistance)[2] <- "hatchDoY"
colnames(temporalDistance)[1] <- "label"
row.names(temporalDistance)<-temporalDistance$label

#Need to have only those which have genetic (and spatial distance) equivalents
geneticDistance<-MCsnps[-c(1:7),] # so there are 86 entries.
df<-merge(temporalDistance, geneticDistance, by="row.names")
df<-df[,2:3]
row.names(df)<-df$label
df$label<-NULL
temporalDistance<-df

#Make the Actual Temporal Distance Matrix
temporalDistance<-dist(temporalDistance)

# Create Spatial Distance Matrix based on Best Nest Estimate as previously determined via Iterated Mantel.
#spatialDistance<-data.frame(larv$Label,(larv$Distance.to.Angle.Crossing..m.-(bestNestEst*(larv$Day.of.Year-((larv$Day.of.Year-(74.308*larv$Mean.Otolith.Length.is.in.Millimetres.for.comparison.with.Adults-4.44361)))+7))))

spatialDistance<-data.frame(larv2$Label,(larv2$Distance.to.Angle.Crossing..m.-(dispersalVelocity*(larv2$Day.of.Year-((larv2$Day.of.Year-(larv2$ageOL)))+7))))

geneticDistance<-MCsnps[-c(1:7),] # so there are 86 entries.
row.names(spatialDistance)<-spatialDistance[,1]
df<-merge(spatialDistance, geneticDistance, by="row.names")
df<-df[,2:3]
row.names(df)<-df$larv2.Label
df$larv2.Label<-NULL
spatialDistance<-df

#Make the Actual Spatial Distance Matrix
spatialDistance<-dist(spatialDistance)

#Make the Actual Genetic Distance Matrix
geneticDistance <- dist(geneticDistance)

plot(hclust(geneticDistance), main="Genetic Distance")
plot(hclust(spatialDistance), main="Spatial Distance")
plot(hclust(temporalDistance), main="Temporal Distance")

#Add the distance matrices together after scaling. Need to scale of Distance dominates with its big numbers.
summedMatrices<- (data.Normalization (as.matrix(spatialDistance),type="n4",normalization="column"))+(data.Normalization (as.matrix(geneticDistance),type="n4",normalization="column"))+(data.Normalization (as.matrix(temporalDistance),type="n4",normalization="column"))


plot(hclust(as.dist(summedMatrices)), main="Summed Matrices")

#Weight Genetic Distance two times more than Spatial and Temporal distances.
weightedSummedMatrices<- (data.Normalization (as.matrix(spatialDistance),type="n4",normalization="column"))+2*(data.Normalization (as.matrix(geneticDistance),type="n4",normalization="column"))+(data.Normalization (as.matrix(temporalDistance),type="n4",normalization="column"))

plot(hclust(as.dist(weightedSummedMatrices)), main="Clade Weighted Summed Matrices")

#pam medoid clustering

require(cluster)
wsm.pam <- pam(weightedSummedMatrices,12)
table(wsm.pam$clustering)
plot(wsm.pam)

