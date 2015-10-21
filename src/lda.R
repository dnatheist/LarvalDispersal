## @knitr lda

# Linear Discriminant Analysis with Jacknifed Prediction 
library(MASS)
fit <- lda(SiteName ~ Delta13C + CNRatio+B+K+V+Rb, data=ImportantVars, na.action="na.omit", CV=TRUE)
#fit # show results
# Assess the accuracy of the prediction
# percent correct for each category of G
ct <- table(ImportantVars$SiteName, fit$class)
diag(prop.table(ct, 1))
# total percent correct
sum(diag(prop.table(ct)))
