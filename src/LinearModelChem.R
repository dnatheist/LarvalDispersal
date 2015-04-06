#Linear Model based on most important elements. Need to confirm that these are the most Important Variables before finilasation.
require(MASS)
require(mvoutlier)

ImportantVars<-ChemAnalCore
ImportantVars<-ImportantVars[c(107,13,14,15,121,122,127,131,139)]
ImportantVars<-ImportantVars[complete.cases(ImportantVars),] #remove any nulls
ImportantVars <- droplevels(ImportantVars)#Not sure why get error without this line

# Normalise the Data
ImpVarScaled<-data.Normalization (ImportantVars,type="n4",normalization="column")
#Check it looks OK
print(ImpVarScaled)[10:20,]

#Create a linear model
model<-lm(ImpVarScaled$Distance.to.Angle.Crossing..m.~ImpVarScaled$B+ImpVarScaled$Na+ImpVarScaled$V+ImpVarScaled$Delta13C+ImpVarScaled$CNRatio+ImpVarScaled$V+ImpVarScaled$K)

summary(model)

anova(model)

par(mfrow=c(2,2))
plot(model)


#####This is new bit below
# Detect Outliers in the  Data
library(mvoutlier)
outliers <- 
  aq.plot(ImportantVars[c("Delta13C","CNRatio","B","K","V","Rb")])
outliers # show list of outliers

# Test Multivariate Normality 
require(mvShapiroTest)
m<- as.matrix(ImportantVars[c("Delta13C","CNRatio","B","K","V","Rb")])
mvShapiro.Test(m)

# Graphical Assessment of Multivariate Normality
x <- m # n x p numeric matrix
center <- colMeans(x) # centroid
n <- nrow(x); p <- ncol(x); cov <- cov(x); 
d <- mahalanobis(x,center,cov) # distances 
qqplot(qchisq(ppoints(n),df=p),d,
       main="QQ Plot Assessing Multivariate Normality",
       ylab="Mahalanobis D2")
abline(a=0,b=1)

# Linear Discriminant Analysis with Jacknifed Prediction 
library(MASS)
fit <- lda(Distance.to.Angle.Crossing..m. ~ Delta13C + CNRatio+B+K+V+Rb, data=ImportantVars, 
           na.action="na.omit", CV=TRUE)
fit # show results
# Assess the accuracy of the prediction
# percent correct for each category of G
ct <- table(ImportantVars$Distance.to.Angle.Crossing..m., fit$class)
diag(prop.table(ct, 1))
# total percent correct
sum(diag(prop.table(ct)))

