library(datasets)
data(airquality)
fit<-lm(Ozone~ Wind+Temp+Solar.R, data=airquality)

library(xtable)
xt<-xtable(summary(fit))
print(xt, type="html")

# Linear Discriminant Analysis with Jacknifed Prediction 
library(MASS)

fit <- lda(G ~ x1 + x2 + x3, data=mtcars, na.action="na.omit", CV=TRUE)
fit # show results
# Assess the accuracy of the prediction
# percent correct for each category of G
ct <- table(mtcars$G, fit$class)
diag(prop.table(ct, 1))
# total percent correct
sum(diag(prop.table(ct)))




# leave-one-out and 6-fold cross-validation prediction error for 
# the mammals data set.
data(mammals, package="MASS")
mammals.glm <- glm(log(brain)~log(body),data=mammals)
cv.err <- cv.glm(mammals,mammals.glm)
cv.err.6 <- cv.glm(mammals, mammals.glm, K=6)

# As this is a linear model we could calculate the leave-one-out 
# cross-validation estimate without any extra model-fitting.
muhat <- mammals.glm$fitted
mammals.diag <- glm.diag(mammals.glm)
cv.err <- mean((mammals.glm$y-muhat)^2/(1-mammals.diag$h)^2)

# leave-one-out and 11-fold cross-validation prediction error for 
# the nodal data set.  Since the response is a binary variable an
# appropriate cost function is
cost <- function(r, pi=0) mean(abs(r-pi)>0.5)

nodal.glm <- glm(r~stage+xray+acid,binomial,data=nodal)
cv.err <- cv.glm(nodal, nodal.glm, cost, K=nrow(nodal))$delta 
cv.11.err <- cv.glm(nodal, nodal.glm, cost, K=11)$delta 

classDF <- data.frame(response = mydata$Y, predicted = round(fitted(mysteps),0))
xtabs(~ predicted + response, data = classDF)
wine.fl <- "http://archive.ics.uci.edu/ml/machine-learning-databases/wine/wine.data"
wine <- read.csv(wine.fl,header = F)
wine.names=c("Alcohol", "Malic acid", "Ash", "Alcalinity of ash", "Magnesium",
             "Total phenols", "Flavanoids", "Nonflavanoid phenols", "Proanthocyanins",
             "Color intensity", "Hue", "OD280/OD315 of diluted wines", "Proline")
colnames(wine)[2:14]=wine.names
colnames(wine)[1]="Class"

wine.lda = lda(Cultivar ~ .,data=wine)

We'll see that most of the modeling functions in R share many things in common. For example, to predict values based on a model, we pass the model object to the predict function along with a data frame containing the observations for which we want predictions:
> pred = predict(wine.lda,wine)

To see what's available from the call to predict, we can look at the names of the pred object:
        > names(pred)
[1] "class"     "posterior" "x"

The predicted classification is stored in the class component of the object returned by predict Now that we've got the predicted classes, we can see how well the classification went by making a cross-tabulation of the real Cultivar with our prediction, using the table function:
> table(wine$Cultivar,pred$class)
predclass
1  2  3
1 59  0  0
2  0 71  0
3  0  0 48

Before we get too excited about these results, remember the caution about predicting values based on models that were built using those values. The error rate we see in the table (0) is probably an overestimate of how good the classification rule is. We can use v-fold cross validation on the data, but it's a bit more complicated with linear discriminant analysis, since we'll actually have to rerun the lda command several times. We could write this out "by hand", but it would be useful to have a function that could do this for us. Here's such a function:
        vlda = function(v,formula,data,cl){
                require(MASS)
                grps = cut(1:nrow(data),v,labels=FALSE)[sample(1:nrow(data))]
                pred = lapply(1:v,function(i,formula,data){
                        omit = which(grps == i)
                        z = lda(formula,data=data[-omit,])
                        predict(z,data[omit,])
                },formula,data)
                
                wh = unlist(lapply(pred,function(pp)pp$class))
                table(wh,cl[order(grps)])
        }

#This function accepts four arguments: v, the number of folds in the cross classification, formula which is the formula used in the linear discriminant analysis, data which is the data frame to use, and cl, the classification variable (wine$Cultivar in this case).
#By using the sample function, we make sure that the groups that are used for cross-validation aren't influenced by the ordering of the data - notice how the classification variable (cl) is indexed by order(grps) to make sure that the predicted and actual values line up properly.
Applying this function to the wine data will give us a better idea of the actual error rate of the classifier:
tt = vlda(5,Cultivar~.,wine,wine$Cultivar) 
tt
wh   1  2  3
1 59  1  0
2  0 69  1
3  0  1 47
