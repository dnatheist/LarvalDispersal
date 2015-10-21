## This is to create files for AMOVA and then run the AMOVA

require(ade4)
##Need three data frames: distance, samples, structures
DArTdm<-DMac14.1567DistMatrix
rm(DMac14.1567DistMatrix) #Tidy up

row.names(DArTdm)<-DArTdm[,1]
DArTdm$X.<-NULL
#Now we have amova matrix $distances 144x144 

#Now to make $samples and then $structures - below all unfinished!!!
siteGroupings<-read.csv("/data/siteGroupings.csv")
source("addNewData.r")
allowedVars <- c("twoGroups","threeGroups")
addNewData("dataNew.csv", myData, allowedVars)
# samples - sites x 6

# merge with larv to get site names

read.csv("data/DMac14-1567snps.csv")

DArTsnps<-DMac14.1567snps
rm(DMac14.1567snps)# Tidy up

larv<-qslLarvaeAgePlus #Make dataframe from larval data
rownames(larv) <- larv[,1] #Make Labels row names too
rm(qslLarvaeAgePlus) #Tidy up
lc<-merge(larv,allOtoChemData, by="LarvalID") 1st merge
# this just from munge files. Needs lots work 26/3/15######################################

# Extract more meaningful names (Labels) in place of Larval ID numbers
DArTsnps<-t(DArTsnps) # Just transpose this way larvae become observations and snps become the variables.
vec<-DArTsnps[,4]# List of DArT 'larval identifiers' for use recoding names with recoderFunc
# which translates the Dart identifiers into something more meaningful (Labels) from the database.

recoderFunc <- function(data, oldvalue, newvalue) {
  
  # convert any factors to characters
  
  if (is.factor(data))     data     <- as.character(data)
  if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
  if (is.factor(newvalue)) newvalue <- as.character(newvalue)
  
  # create the return vector
  
  newvec <- data
  
  # put recoded values into the correct position in the return vector
  
  for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]
  
  newvec
  
}

vec<-recoderFunc(vec, lc$LarvalRecords.LarvaID, lc$Label)

DArTsnps[,4]<-vec # Load more meaningful names (Labels) in place of Larval ID numbers

rownames(DArTsnps) <- DArTsnps[,4] #Make rownames meaningful
DArTsnps<-as.data.frame(DArTsnps) #Make a dataframe
DArTsnps$V4 <- NULL #Remove V4 as now redundant.

#New DArTsnps data frame now has correct labels 

#Now to delete data not required for distance matrix creation.
DArTsnps<-DArTsnps[-c(1:18),] # 18 columns of DArT calculations and allele sequences
DArTsnps<-DArTsnps[,-c(1:4)] # 4 columns of DArT descriptors

allsnps<-DArTsnps #Rename variable as it is now all SNPs (MC,TC and others)

#Create a data frame with only the Murray Cod SNPs
MCsnps <- allsnps[-c(9, 10, 50, 85, 2, 61), ]

#Create a data frame with only Trout cod and hybrid SNPs.
TCsnps <- allsnps[c(9, 10, 50, 85, 2, 61), ]
#remove TC and extraneous rows.




#structures - zones A (above sand slug), B (below sand slug and above red rocks/KP), C(below red rocks/KP)
#copy sites and recode 1-144 with A, B or C based on site names. 
