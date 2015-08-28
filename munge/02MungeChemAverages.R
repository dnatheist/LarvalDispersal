#Draft code for extracting elemental averages and sd for each otolith core and shoving it in a dataframe for export to database, for subsequent analysis with delta C delta N etc.Seems to be working fine.

##Create new empty dataframe for the use of, with element headings.
EachOtoData<-allOtoChemData
## Empty it but keep column names
EachOtoData<- EachOtoData[which(is.na(EachOtoData$text)), ]

for (l in (unique(allOtoChemData[,2]))){ # for each larval otolith
  df<-subset(allOtoChemData, LarvalID == l & OtoPart == "Core")
  EachOtoData[l,3]<-df[1,3] #Add SiteName to EachOtoData DF
  EachOtoData[l,1] <- l
  for (e in 5:68){ # for each element
    mean<-mean(df[,e])
    #sd<-sd(df[,e])
    #n<-length(unique(allOtoChemData[,66]))
    EachOtoData[l,e]<- mean 
    EachOtoData$OtoPart<-"Core"
    
  } 
}
EachOtoData$LarvalID<-NULL
EachOtoData$Time<-NULL
colnames(EachOtoData)[1]<-"LarvalID" 
EachOtoData<-EachOtoData[,-c(3:34)] #Remove ratios
EachOtoData<-EachOtoData[complete.cases(EachOtoData),]   #final[complete.cases(final),]

write.csv(EachOtoData,file= ".\\reports\\ElemForEachOtolithCore.csv") # NB One "." and back slashes.
#Now to merge Labels and Sites
colnames(larv)[2]<-"LarvalID"
ChemAnalCore<-merge(larv,EachOtoData,by = "LarvalID")
#Clean up
rm(df)

#####################################################Repeat for edge1

##Create new empty dataframe for the use of, with element headings.
EachOtoData<-allOtoChemData
## Empty it but keep column names
EachOtoData<- EachOtoData[which(is.na(EachOtoData$text)), ]

for (l in (unique(allOtoChemData[,2]))){ # for each larval otolith
  df<-subset(allOtoChemData, LarvalID == l & OtoPart == "E1")
  EachOtoData[l,3]<-df[1,3] #Add SiteName to EachOtoData DF
  EachOtoData[l,1] <- l
  for (e in 5:68){ # for each element
    mean<-mean(df[,e])
    #sd<-sd(df[,e])
    #n<-length(unique(allOtoChemData[,66]))
    EachOtoData[l,e]<- mean 
    EachOtoData$OtoPart<-"E1"
    
  } 
}
EachOtoData$LarvalID<-NULL
EachOtoData$Time<-NULL
colnames(EachOtoData)[1]<-"LarvalID" 
EachOtoData<-EachOtoData[,-c(3:34)] #Remove ratios
EachOtoData<-EachOtoData[complete.cases(EachOtoData),]   #final[complete.cases(final),]

write.csv(EachOtoData,file= ".\\reports\\ElemForEachOtolithE1.csv") # NB One "." and back slashes.

#Now to merge Labels and Sites
colnames(larv)[2]<-"LarvalID"
ChemAnalE1<-merge(larv,EachOtoData,by = "LarvalID")
#Clean up
rm(df)


#####################################################Repeat for edge2

##Create new empty dataframe for the use of, with element headings.
EachOtoData<-allOtoChemData
## Empty it but keep column names
EachOtoData<- EachOtoData[which(is.na(EachOtoData$text)), ]

for (l in (unique(allOtoChemData[,2]))){ # for each larval otolith
  df<-subset(allOtoChemData, LarvalID == l & OtoPart == "E2")
  EachOtoData[l,3]<-df[1,3] #Add SiteName to EachOtoData DF
  EachOtoData[l,1] <- l
  for (e in 5:68){ # for each element
    mean<-mean(df[,e])
    #sd<-sd(df[,e])
    #n<-length(unique(allOtoChemData[,66]))
    EachOtoData[l,e]<- mean 
    EachOtoData$OtoPart<-"E2"
    
  } 
}
EachOtoData$LarvalID<-NULL
EachOtoData$Time<-NULL
colnames(EachOtoData)[1]<-"LarvalID" 
EachOtoData<-EachOtoData[,-c(3:34)] #Remove ratios
EachOtoData<-EachOtoData[complete.cases(EachOtoData),]   #final[complete.cases(final),]

write.csv(EachOtoData,file= ".\\reports\\ElemForEachOtolithE2.csv") # NB One "." and back slashes.

#Now to merge Labels and Sites
colnames(larv)[2]<-"LarvalID"
ChemAnalE2<-merge(larv,EachOtoData,by = "LarvalID")
#Clean up
rm(df)


#####################################################Repeat for 0

##Create new empty dataframe for the use of, with element headings.
EachOtoData<-allOtoChemData
## Empty it but keep column names
EachOtoData<- EachOtoData[which(is.na(EachOtoData$text)), ]

for (l in (unique(allOtoChemData[,2]))){ # for each larval otolith
  df<-subset(allOtoChemData, LarvalID == l & OtoPart == "0")
  EachOtoData[l,3]<-df[1,3] #Add SiteName to EachOtoData DF
  EachOtoData[l,1] <- l
  for (e in 5:68){ # for each element
    mean<-mean(df[,e])
    #sd<-sd(df[,e])
    #n<-length(unique(allOtoChemData[,66]))
    EachOtoData[l,e]<- mean 
    EachOtoData$OtoPart<-"0"
    
  } 
}
EachOtoData$LarvalID<-NULL
EachOtoData$Time<-NULL
colnames(EachOtoData)[1]<-"LarvalID" 
EachOtoData<-EachOtoData[,-c(3:34)] #Remove ratios
EachOtoData<-EachOtoData[complete.cases(EachOtoData),]   #final[complete.cases(final),]

write.csv(EachOtoData,file= ".\\reports\\ElemForEachOtolith0.csv") # NB One "." and back slashes.

#Clean up
rm(df)
rm(e)
rm(l)
rm(mean)
#Now to merge Labels and Sites
colnames(larv)[2]<-"LarvalID"
ChemAnal0<-merge(larv,EachOtoData,by = "LarvalID")
##############################

ChemAnal0<-ChemAnal0[-26,] # an extreme B outlier that was missed at data clean

# Now make two more:
ChemAnalAll4<-rbind(ChemAnal0,ChemAnalCore,ChemAnalE1,ChemAnalE2)
ChemAnalNonCore<-rbind(ChemAnal0,ChemAnalE1,ChemAnalE2)

rm(EachOtoData) #clean up last temporary file added line 20/8/2015



