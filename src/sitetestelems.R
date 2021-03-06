
require(dplyr)
#This identifies the important Chemistry factors that could be used to predict site the larva are from.
#First using Core as the values. So need to combine first
ImportantVars<-ChemAnalCore

ImportantVars<- select(ImportantVars, SiteName, Li:Pb)#ImportantVars<-ImportantVars[c(84,120:151)] #
ImportantVars$SiteName[ImportantVars$SiteName==""]  <- NA
ImportantVars<-ImportantVars[complete.cases(ImportantVars),] #remove any nulls
colnames(ImportantVars)[1]<-"SiteName"
ImportantVars$SiteName<-as.factor(ImportantVars$SiteName)

ImportantVars <- droplevels(ImportantVars)#Not sure why get error without this line.

library(randomForest)
forest <- randomForest(SiteName ~.,data=ImportantVars, importance=TRUE)

forest

library(ggplot2)
forest.importance = as.data.frame(importance(forest, scale=FALSE))
forest.importance = forest.importance[,1:(ncol(forest.importance)-2)]
forest.importance$mean = rowMeans(forest.importance)

#forest.importance

ggplot(forest.importance, aes(x=row.names(forest.importance), y=mean)) +
        ylab('mean relative feature importance') +
        xlab('feature') +
        geom_bar(stat='identity')
##################################################################################
#Now to do the same but using Edge1 as the values. So need to combine first
ImportantVars<-ChemAnalE1

ImportantVars<- select(ImportantVars, Genetics.NestID, Li:Pb)#ImportantVars<-ImportantVars[c(84,120:151)] #
ImportantVars<-ImportantVars[complete.cases(ImportantVars),] #remove any nulls
colnames(ImportantVars)[1]<-"Genetics.NestID"
ImportantVars$Genetics.NestID<-as.factor(ImportantVars$Genetics.NestID)

ImportantVars <- droplevels(ImportantVars)#Not sure why get error without this line.

library(randomForest)
forest <- randomForest(Genetics.NestID ~.,data=ImportantVars, importance=TRUE)

forest

library(ggplot2)
forest.importance = as.data.frame(importance(forest, scale=FALSE))
forest.importance = forest.importance[,1:(ncol(forest.importance)-2)]
forest.importance$mean = rowMeans(forest.importance)

#forest.importance

ggplot(forest.importance, aes(x=row.names(forest.importance), y=mean)) +
        ylab('mean relative feature importance') +
        xlab('feature') +
        geom_bar(stat='identity')

################################################################################################
#Now to do the same but using Edge2 as the factor rather than Genetics.NestID. So need to combine first
ImportantVars<-ChemAnalE2

ImportantVars<- select(ImportantVars, Genetics.NestID, Li:Pb)#ImportantVars<-ImportantVars[c(84,120:151)] #
ImportantVars<-ImportantVars[complete.cases(ImportantVars),] #remove any nulls
colnames(ImportantVars)[1]<-"Genetics.NestID"
ImportantVars$Genetics.NestID<-as.factor(ImportantVars$Genetics.NestID)

ImportantVars <- droplevels(ImportantVars)#Not sure why get error without this line.

library(randomForest)
forest <- randomForest(Genetics.NestID ~.,data=ImportantVars, importance=TRUE)

forest

library(ggplot2)
forest.importance = as.data.frame(importance(forest, scale=FALSE))
forest.importance = forest.importance[,1:(ncol(forest.importance)-2)]
forest.importance$mean = rowMeans(forest.importance)

#forest.importance

ggplot(forest.importance, aes(x=row.names(forest.importance), y=mean)) +
        ylab('mean relative feature importance') +
        xlab('feature') +
        geom_bar(stat='identity')
##################################################################################################
#Now to do the same but using Edge0 values. So need to combine first
ImportantVars<-ChemAnal0

ImportantVars<- select(ImportantVars, Genetics.NestID, Li:Pb)#ImportantVars<-ImportantVars[c(84,120:151)] #
ImportantVars<-ImportantVars[complete.cases(ImportantVars),] #remove any nulls
colnames(ImportantVars)[1]<-"Genetics.NestID"
ImportantVars$Genetics.NestID<-as.factor(ImportantVars$Genetics.NestID)

ImportantVars <- droplevels(ImportantVars)#Not sure why get error without this line.

library(randomForest)
forest <- randomForest(Genetics.NestID ~.,data=ImportantVars, importance=TRUE)

forest

library(ggplot2)
forest.importance = as.data.frame(importance(forest, scale=FALSE))
forest.importance = forest.importance[,1:(ncol(forest.importance)-2)]
forest.importance$mean = rowMeans(forest.importance)

#forest.importance

ggplot(forest.importance, aes(x=row.names(forest.importance), y=mean)) +
        ylab('mean relative feature importance') +
        xlab('feature') +
        geom_bar(stat='identity')

##################################################################################################
#Now to do the same but using All4 values. So need to combine first
ImportantVars<-ChemAnalAll4

ImportantVars<- select(ImportantVars, Genetics.NestID, Li:Pb)#ImportantVars<-ImportantVars[c(84,120:151)] #
ImportantVars<-ImportantVars[complete.cases(ImportantVars),] #remove any nulls
colnames(ImportantVars)[1]<-"Genetics.NestID"
ImportantVars$Genetics.NestID<-as.factor(ImportantVars$Genetics.NestID)

ImportantVars <- droplevels(ImportantVars)#Not sure why get error without this line.

library(randomForest)
forest <- randomForest(Genetics.NestID ~.,data=ImportantVars, importance=TRUE)

forest

library(ggplot2)
forest.importance = as.data.frame(importance(forest, scale=FALSE))
forest.importance = forest.importance[,1:(ncol(forest.importance)-2)]
forest.importance$mean = rowMeans(forest.importance)

#forest.importance

ggplot(forest.importance, aes(x=row.names(forest.importance), y=mean)) +
        ylab('mean relative feature importance') +
        xlab('feature') +
        geom_bar(stat='identity')

##################################################################################################
#Now to do the same but using Non-Core values. So need to combine first
ImportantVars<-ChemAnalNonCore

ImportantVars<- select(ImportantVars, Genetics.NestID, Li:Pb)#ImportantVars<-ImportantVars[c(84,120:151)] #
ImportantVars<-ImportantVars[complete.cases(ImportantVars),] #remove any nulls
colnames(ImportantVars)[1]<-"Genetics.NestID"
ImportantVars$Genetics.NestID<-as.factor(ImportantVars$Genetics.NestID)

ImportantVars <- droplevels(ImportantVars)#Not sure why get error without this line.

library(randomForest)
forest <- randomForest(Genetics.NestID ~.,data=ImportantVars, importance=TRUE)

forest

library(ggplot2)
forest.importance = as.data.frame(importance(forest, scale=FALSE))
forest.importance = forest.importance[,1:(ncol(forest.importance)-2)]
forest.importance$mean = rowMeans(forest.importance)

#forest.importance

ggplot(forest.importance, aes(x=row.names(forest.importance), y=mean)) +
        ylab('mean relative feature importance') +
        xlab('feature') +
        geom_bar(stat='identity')