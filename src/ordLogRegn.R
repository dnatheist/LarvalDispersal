## @knitr ordLogRegn

require(clusterSim)
ImportantVars<-ChemAnalCore
ImportantVars<-ImportantVars[c(84,13,14,15,121,122,127,131,139)]
ImportantVars<-ImportantVars[complete.cases(ImportantVars),] #remove any nulls
ImportantVars <- droplevels(ImportantVars)#Not sure why get error without this line

# Normalise the Data
ImpVarScaled<-data.Normalization (ImportantVars,type="n4",normalization="column")
#Check it looks OK
print(ImpVarScaled)[10:20,]

ggplot(ImpVarScaled, aes(x = as.factor(SiteName), y = V)) +
  geom_boxplot(size = .75) +
  geom_jitter(alpha = .5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

#Ordinal Logistic Regression
## fit ordered logit model and store results 'model'
model <- polr(as.factor(SiteName) ~ B + Na + V + K+ Delta13C + CNRatio, data = ImpVarScaled, Hess=TRUE)

## view a summary of the model
summary(model)

## store table
(ctable <- coef(summary(model)))
## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
## combined table
(ctable <- cbind(ctable, "p value" = p))

#or just do confidence intervals
ci<-confint(model)


