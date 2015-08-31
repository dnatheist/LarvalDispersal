require(aplpack)
faces()
faces(face.type=1)

############ Chemistry Analysis Core
data<-ChemAnalCore[1:25,c(1,120:151)]
row.names(data)<-data$LarvalID
faceTest<-data[,2:33]
faces(faceTest)

plot(data[1:25,2:3],bty="n")
a<-faces(data[1:25,],plot=FALSE)
plot.faces(a,data[1:25,2],data[1:25,3],width=.03,height=1)

############ muck around
faces(rbind(1:3,5:3,3:5,5:7))
data<-EachOtoData[1:25,c(1,3:4)]
faceTest<-data[1:25,2:3]
faces(faceTest)

plot(data[1:25,2:3],bty="n")
a<-faces(data[1:25,],plot=FALSE)
plot.faces(a,data[1:25,2],data[1:25,3],width=.03,height=.3)

plot(FaceTest[1:25,3:4],bty="n")
a<-faces(faceTest)
plot.faces(a,faceTest[1:25,3:4],faceTest[1:25,3:4],width=35,height=30)

set.seed(17)
faces(matrix(sample(1:1000,128,),16,8),main="random faces")

a<-faces(rbind(1:3,5:3,3:5,5:7),plot.faces=FALSE)
plot(0:5,0:5,type="n")
plot(a,x.pos=1:4,y.pos=1:4,1.5,0.7)


tmp<-withSibsLogical
#First sort
tmp<-as.data.frame(tmp)
tmp$sort<-row.names(tmp)
tmp <- tmp[order(tmp$sort),]#sort row order
tmp$sort<-NULL
tmp<-tmp[,order(names(tmp))]#sort column order
tmp<-as.matrix(tmp)
tmp[upper.tri(tmp)] <- NA

inds<-which(withSibsLogical == TRUE, arr.ind=TRUE)
int<-as.data.frame(rownames(withSibsLogical)[inds[,1]],colnames(withSibsLogical)[inds[,2]])

