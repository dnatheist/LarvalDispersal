#This file is to generate names for nests. By making rules to do it based only on location and day of hatch I can then compare to siblings list. This should give and idea of error maybe.

#First run I used first letter of catch site name between distance and hatch day. Probably better off using a letter representing the year. 11,12,13 el,tw,th (or en,te,tn) perhaps.

# I suppose this is actually a classification problem. Perhaps should treat as such, do confusion matrix etc.


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
ndDF<-subset(larv2, select=c(Label,nestSeg,roundedHatchDoY,siteLetter))

#create the name (still no sibs accounted for)
ndDF$nestIDnoSibs<-paste(ndDF$nestSeg,ndDF$siteLetter,ndDF$roundedHatchDoY, sep="")


#Calculate a family group
families<-as.data.frame(cutree(dend2, h = 26.5))
larvalClades<-merge(larvalClades,families, by="row.names")
names(larvalClades)[names(larvalClades)=="cutree(dend2, h = 26.5)"] <- "Family"
names(larvalClades)[names(larvalClades)=="Row.names"] <- "Label"


#calc nest distance
larv$nestDist<-((larv$Distance.to.Angle.Crossing..m.-(dispersalVelocity*(larv$Day.of.Year-((larv$Day.of.Year-(larv$ageOL)))+larv$incubTime))))

new<-merge(larv,larvalClades, by="Label")
newsub<-new


#Create a data frame from which to develop nestID names (not yet combining siblings)
new$roundedHatchDoY<-round(new$hatchDoY)
new$nestSeg<-round(new$nestDist/100)
new$siteLetter<-substr(new$Label, 1, 1)

famNo<-unique(new$Family)
final<- newsub[FALSE,]

for (nn in famNo){
        newsub<-subset(new,Family==nn, select=c(LarvalID,Label,nestSeg,YearOnly,siteLetter,roundedHatchDoY,Family))
        newsub$nestID<-paste(round(mean(newsub$nestSeg)),newsub$siteLetter,round(mean(newsub$roundedHatchDoY)), sep="")
        final<-rbind(final,newsub)
}
plot(final$nestSeg,as.factor(final$nestID))
nrow(as.data.frame(unique(final$nestID)))

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
segLogical<-seg<7 #Less than 500m apart
doyLogical<-doy<5 #Hatch day less than days apart

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
names(sibEstPairs)[names(sibEstPairs) == "rownames.sameNestEstSubset..k...1.."] <- 'rown'
names(sibEstPairs)[names(sibEstPairs) == "colnames.sameNestEstSubset..k...2.."] <- 'coln'

#Order by row (effectively LarvalID) so each Larval ID has all siblings listed.
sibEstPairs<-sibEstPairs[order(sibEstPairs$rown),]

#Create an actual sibs list so we can compare the two
###########
# Recreate Genetic Distance
# Create a Murray Cod distance matrix
MCdm<-MCsnps[-c(1:7),]
MCdm <- dist(MCdm)
MCdm<-as.matrix(MCdm)
MCdm<-as.data.frame(MCdm)
#This is to be used for plotting
MCdmForSibs<-MCdm #need this matrix for sibling analysis .r

#Then create actual sibs list
sibsSubset<-subset(MCdmForSibs, apply(MCdmForSibs, 1, function(MCdmForSibs){any(MCdmForSibs > 0 & MCdmForSibs < 26.5)}))
#The siblings have a dissimilarity index of less than 26.5. This can be seen in the dendrogram and in the scatterplot for IBD.

sibsLogical<-sibsSubset<26.5 #a logical DF to record T and F for sibling pairs
#sibsSubset<-sibsLogical*sibsSubset #This turns to 0 all the genetic similarity above 26.5 (no-sibs)

#sibsSubset<-sibsSubset[!sapply(sibsSubset, function(x) all(x == 0))] #now remove those with no sibs
#
#sibSubset<-sibSubset<0

sibsLogical[lower.tri(sibsLogical, diag=TRUE)] <- NA

#########
ai<-which(sibsLogical ==TRUE)
k <- arrayInd(ai, dim(sibsLogical))

#sibEstPairs<-data.frame(x=1,y=2)
r<-as.data.frame(rownames(sibsLogical)[k[,1]])
c<-as.data.frame(colnames(sibsLogical)[k[,2]])
rm(sibsLogical)
sibActualPairs<-as.data.frame(c(r,c))

#tidy up var names
names(sibActualPairs)[names(sibActualPairs) == "rownames.sibsLogical..k...1.."] <- 'rown'
names(sibActualPairs)[names(sibActualPairs) == "colnames.sibsLogical..k...2.."] <- 'coln'

#Order by row (effectively LarvalID) so each Larval ID has all siblings listed.
sibActualPairs<-sibActualPairs[order(sibActualPairs$rown),]
###########

intersect(sibActualPairs$coln, sibEstPairs$coln)

########################


