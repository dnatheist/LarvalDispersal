source("src/dart.r")
source("src/read.dart.r")

#library(devtools)
#dev_mode(on=T)
#dev_mode(on=F)

#Install non- cran versions (the older version of adgenet has a serious flaw and is not updated yet.)
#install_github("thibautjombart/adegenet")
#install_github("green-striped-gecko/PopGenReport")


mc.dart <- read.dart("data/DMac14-1567snps2.csv", topskip = 5)
gl.dart <- dart2genlight(mc.dart, covfilename = "data/qslDartCovariates.csv")

library(StAMPP)
system.time(snpfst <-stamppFst(gl.dart,nboots=1, percent=95, nclusters=8 ))


#pcoa

system.time(pca1 <- glPca(gl.dart  , parallel=F, nf=3))
s.class(pca1$scores, pop(gl.dart), col=rainbow(nlevels(pop(gl.dart))))




#####
mc.gen <- genlight2genind(gl.dart)

#% !Rnw weave = Sweave
library(PopGenReport)
pgr <- popgenreport(mc.gen )
