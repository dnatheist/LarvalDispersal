
#Set the path to data
fp<-file.path("C:\\Users\\s428825\\Google Drive\\PhD Analyses\\LarvalDispersal\\data\\")

DArTsnps<-DMac14.1567snps
#rm(DMac14.1567snps)# Tidy up

DArTdm<-DMac14.1567DistMatrix
rm(DMac14.1567DistMatrix) #Tidy up

larv<-qslLarvaeAgePlus #Make dataframe from larval data
rownames(larv) <- larv[,1] #Make Labels row names too
rm(qslLarvaeAgePlus) #Tidy up

# Extract more meaningful names (Labels) in place of Larval ID numbers
DArTsnps<-t(DArTsnps) # Just transpose this way larvae become observations and snps become the variables.
vec<-DArTsnps[,4]# List of DArT 'larval identifiers' for use recoding names with recoderFunc
# which translates the DaRT identifiers into something more meaningful (Labels) from the database.

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

vec<-recoderFunc(vec, larv$LarvalRecords.LarvaID, larv$Label)

DArTsnps[,4]<-vec # Load more meaningful names (Labels) in place of Larval ID numbers
rm(vec)# Tidy up

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


