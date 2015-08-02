#This file is to generate names for nests. By making rules to do it based only on location and day of hatch I can then compare to siblings list. This should give and idea of error maybe.


#Establish Nest Distance based on Dispersal Velocity
dispersalVelocity<-1200
larv1<-merge(larv,qslGeneticsForNestChapter, by.x="LarvalID", by.y="LarvaID")
rownames(CladeNamesToMerge) <- CladeNamesToMerge[,1]
CladeNamesToMerge$Label<-rownames(CladeNamesToMerge)
larv2<-merge(larv1,CladeNamesToMerge, by="Label")
larv2<-merge(larv2,qslAgeData, by.x="LarvalID", by.y="LarvaID")
# To create a distance using the previously calculated best estimate of drift velocity (m/d available since leaving brood care)
larv2$nestDist<-larv2$Distance.to.Angle.Crossing..m.-(dispersalVelocity*(larv2$Day.of.Year-(larv2$hatchDoY+7))) #this should come from DB perhaps

plot(larv2$nestDist,larv2$clade, xlab="Nest Distance (m)", ylab = "Clade")

#Create a data frame from which to develop nestID names (not yet combining siblings)
larv2$roundedHatchDoY<-round(larv2$hatchDoY)
larv2$nestSeg<-round(larv2$nestDist/100)
larv2$siteLetter<-substr(larv2$SiteName, 1, 1)
ndDF<-subset(larv2, select=c(LarvalID,nestSeg,roundedHatchDoY,siteLetter))

#create the name (still no sibs accounted for)
ndDF$nestIDnoSibs<-paste(ndDF$nestSeg,ndDF$siteLetter,ndDF$roundedHatchDoY, sep="")

#Output a file
#write.csv(ndDF, file="ndDF.csv")

#####################
#Make some matrices to compare. First segment of river, then hatchDoY
row.names(ndDF)<-ndDF[,1] # Make sensible row names (larvalID)
ndDF$LarvalID<-NULL

seg<-as.matrix(dist(subset(ndDF, select=c(nestSeg))))
doy<-as.matrix(dist(subset(ndDF, select=c(roundedHatchDoY))))

seg[lower.tri(seg, diag=TRUE)] <- NA

#Create logical matrices
segLogical<-seg<5 #Less than 500m apart
doyLogical<-doy<3 #Hatch day less than days apart

#Combine them
sameNestEst<-segLogical*doyLogical
sameNestEst[lower.tri(sameNestEst, diag=TRUE)] <- NA

sameNestEst<-as.data.frame(sameNestEst)
#Line to extract subset of DF which include at least one sibling in the dataframe
sameNestEstSubset<-subset(sameNestEst, apply(sameNestEst, 1, function(sameNestEst){any(sameNestEst >0)}))


rowSums(sameNestEstSubset,na.rm=TRUE)

##############
# col<-colnames(sameNestEstSubset)[apply(sameNestEstSubset, 2, function(u) any(u==1))]
# row<-row.names(sameNestEstSubset)[apply(sameNestEstSubset, 1, function(u) any(u==1))]
# intersect(row,col)
# union(row,col)
##############
  
#newdata <- sameNestEstSubset[ which(sameNestEstSubset, any(u==1), )]
#apply(sameNestEstSubset, 2, function(u) any(u==1))
#vv<-apply(sameNestEstSubset, 2, function(u) any(u==1))
################

ai<-which(sameNestEst ==1)
k <- arrayInd(ai, dim(sameNestEstSubset))

#sibEstPairs<-data.frame(x=1,y=2)
r<-as.data.frame(rownames(sameNestEstSubset)[k[,1]])
c<-as.data.frame(colnames(sameNestEstSubset)[k[,2]])

sibEstPairs<-as.data.frame(c(r,c))

#tidy up var names
names(sibEstPairs)[names(sibEstPairs) == "rownames.sameNestEstSubset..k...1.."] <- 'row'
names(sibEstPairs)[names(sibEstPairs) == "colnames.sameNestEstSubset..k...2.."] <- 'col'

#sort
sibEstPairs<-sibEstPairs[order(row),] #not yet working
