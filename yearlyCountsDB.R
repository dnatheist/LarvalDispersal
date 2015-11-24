#conditional counts examples for inline r

#table(qslAllLarvaInfo$Label,qslAllLarvaInfo$YearOnly)



nrow(subset(qslAllLarvaInfo, qslAllLarvaInfo$YearOnly==2011))
nrow(subset(qslAllLarvaInfo, qslAllLarvaInfo$YearOnly==2012))
nrow(subset(qslAllLarvaInfo, qslAllLarvaInfo$YearOnly==2013))
               
nrow(qslAllLarvaInfo)               
               
               
