## @knitr siblingAnalysis

## This file is to examine the siblings collected that are known to be siblings from DaRT sequencing. The idea is to learn:
#       1. What duration is hatching in wild nests? 
#       2. Create a list of larvae with siblings (withSibs) and maybe (withoutSibs?) for checking
#       3. Create a logical matrix of siblinf relationships to re-use (sibsLogical)

#These data can then be used to inform Nest assignments in cases where clade and location are the same.

#Line to extract subset of DF which include at least one sibling in the dataframe
sibsSubset<-subset(MCdmForSibs, apply(MCdmForSibs, 1, function(MCdmForSibs){any(MCdmForSibs > 0 & MCdmForSibs < 26.5)}))
#The siblings have a dissimilarity index of less than 26.5. This can be seen in the dendrogram and in the scatterplot for IBD.

sibsLogical<-sibsSubset<26.5 #a logical DF to record T and F for sibling pairs
sibsSubset<-sibsLogical*sibsSubset #This turns to 0 all the genetic similarity above 26.5 (no-sibs)

sibsSubset<-sibsSubset[!sapply(sibsSubset, function(x) all(x == 0))] #now remove those with no sibs
sibsLogical<-sibsSubset>0

#Make a list (data frame) with names 
withSib<-row.names(sibsSubset) #make a list of those with sibs
withSibDF<-as.data.frame(withSib)
row.names(withSibDF)<-withSibDF[,1] # Make sensible row names (larva labels)
withSibDF<-as.data.frame(withSibDF) #make back to a df

##Now to see how many days between hatch of sibs.

withSib2<-subset(larv, select=c(Label,hatchdoy)) #get labels and hatch DoY
withSib2<-withSib2[withSib2$Label %in% withSib,] #reduce df to the siblings only

#row.names(withSib2)<-withSib2[,2] #make label row names for making distance matrix
withSib2$Label<-NULL #tidy up redundant

#Then to make a distance matrix.
withSibDM<-dist(withSib2)
withSibDM2<-as.matrix(withSibDM)

###sort rows and columns of both before we can do matrix algebra
withSibDM2<-as.data.frame(withSibDM2)
withSibDM2$sort<-row.names(withSibDM2)
withSibDM2 <- withSibDM2[order(withSibDM2$sort),]#sort row order
withSibDM2$sort<-NULL
withSibDM2<-withSibDM2[,order(names(withSibDM2))]#sort column order
withSibDM2<-as.matrix(withSibDM2)

sibsLogical<-as.data.frame(sibsLogical)
sibsLogical$sort<-row.names(sibsLogical)
sibsLogical <- sibsLogical[order(sibsLogical$sort),]#sort row order
sibsLogical$sort<-NULL
sibsLogical<-sibsLogical[,order(names(sibsLogical))]#sort column order
sibsLogical<-as.matrix(sibsLogical)
####

sibsHatchDiff<-sibsLogical*withSibDM2
sibsHatchDiff[sibsHatchDiff == 0] <- NA
sibsHatchDiffMeans<-as.data.frame(colMeans(sibsHatchDiff, na.rm=TRUE))

#Get rid of NAs
sibsHatchDiffMeans<-na.omit(sibsHatchDiffMeans)

#Output something
write.csv(format(sibsHatchDiffMeans), file="sibsHatchDiffMeans.csv")

#Describe and Plot Summary Statistics regading the difference in hatch day-of-year between siblings
describe(sibsHatchDiffMeans)

hist(sibsHatchDiffMeans,breaks=seq(0,5,l=8),freq=FALSE,col="orange",main="Histogram", xlab="x",ylab="f(x)",yaxs="i",xaxs="i")
#So the three objective have been achieved. 
#       1. What duration is hatching in wild nests? Description above showns mean etc.
#       2. Create a list of larvae with siblings (withSibsDF) done
#       3. Create a logical matrix of siblinf relationships to re-use (sibsLogical) done

#Clean Up
rm(sibsSubset)
rm(sibsHatchDiff)
rm(sibsHatchDiffMeans)
rm(sibsLogical)
rm(withSib)
rm(withSib2)
rm(withSibDF)
rm(withSibDM)
rm(withSibDM2)
rm(MCdmForSibs)

#END