#Take the three matrices (genetic distance, spatial distance, temporal distance) add the matrices and then cluster to give a better approximation of individual nests. 

#First need to create temporal distance from the larv data set so need label and hatch day of the year.

#Create Hatch DoY and Temporal Distance Matrix

TemporalDist<-data.frame(larv$Label,(larv$Day.of.Year-(74.308*larv$Mean.Otolith.Length.is.in.Millimetres.for.comparison.with.Adults-4.44361)))
colnames(TemporalDist)[2] <- "hatchDoY"
colnames(TemporalDist)[1] <- "label"
row.names(TemporalDist)<-TemporalDist$label

#Need to have only those which have genetic (and spatial distance) equivalents
MCdm<-MCsnps[-c(1:7),] # so there are 86 entries.
df<-merge(TemporalDist, MCdm, by="row.names")
df<-df[,2:3]
row.names(df)<-df$label
df$label<-NULL
TemporalDist<-df

#TemporalDist<-TemporalDist[complete.cases(TemporalDist),]


#Make the Actual Distance Matrix
TemporalDist<-dist(TemporalDist)

TemporalCluster<-hclust(TemporalDist)
plot(TemporalCluster)




geoDist<-data.frame(larv$Label,(larv$Distance.to.Angle.Crossing..m.-(1000*(larv$Day.of.Year-((larv$Day.of.Year-(74.308*larv$Mean.Otolith.Length.is.in.Millimetres.for.comparison.with.Adults-4.44361)))+7))))
MCdm<-MCsnps[-c(1:7),] # so there are 86 entries.
row.names(geoDist)<-geoDist[,1]
df<-merge(geoDist, MCdm, by="row.names")
df<-df[,2:3]
row.names(df)<-df$larv.Label
df$larv.Label<-NULL
geoDist<-df
geoDist<-dist(geoDist)

MCdm <- dist(MCdm)


#Add the distance matrices after scaling.Need to scale of Distance dominates with its big numbers.
scaleTest<- (data.Normalization (as.matrix(geoDist),type="n4",normalization="column"))+(data.Normalization (as.matrix(MCdm),type="n4",normalization="column"))+(data.Normalization (as.matrix(TemporalDist),type="n4",normalization="column"))

scaleTestCluster<-hclust(as.dist(scaleTest))
plot(scaleTestCluster)

#Weight clade more heavily than Spatial and Temporal seperation.
scaleTest<- (data.Normalization (as.matrix(geoDist),type="n4",normalization="column"))+2*(data.Normalization (as.matrix(MCdm),type="n4",normalization="column"))+(data.Normalization (as.matrix(TemporalDist),type="n4",normalization="column"))

scaleTestCluster<-hclust(as.dist(scaleTest))
plot(scaleTestCluster)

#Is this circular reasoning?
