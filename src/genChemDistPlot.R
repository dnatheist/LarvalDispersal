## @knitr genChemPlot

library(ade4)

MCdmPrep<-MCsnps[-c(1:7),]

test<-merge(MCdmPrep,ImpVarScaled,by = "row.names")
genDist<-dist(test[,c(1:21076)])
chemDist<-dist(test[,c(21077:21086)])

require(plotrix)
op <- par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
plot(genDist,chemDist, col="black", pch=21, bg = "grey", cex = 2,
     ylab="", xlab="", axes=F)
axis(1)
axis(2) 
reg1 <- lm(chemDist~genDist)
ablineclip(reg1, lwd=3) 
par(las=0)
mtext("Genetic Distance", side=1, line=2.5, cex=1.5)
mtext("Chemotype Distance", side=2, line=3.7, cex=1.5)

cor(genDist,chemDist)
rcorr(genDist[lower.tri(genDist)],chemDist[lower.tri(chemDist)])#(x, type="pearson") # type can be pearson or spearman

mant<-mantel.rtest(genDist, chemDist, nrepet = 9999)
mant

#Check all is in order
as.matrix(genDist)[1:5, 1:5]
as.matrix(chemDist)[1:5, 1:5]
