## @knitr siblingAnalysis

#Line to extract subset of DF with just the siblings
#The siblings have a dissimilarity index of less than 26.5. This can be seen in the dendrogram and in the scatterplot for IBD. 

sibSub<-subset(MCdm, apply(MCdm, 1, function(MCdm){any(MCdm > 0 & MCdm < 26.5)}))
write.csv(format(sub), file="sibs.csv")

#I want to replace the index number with the difference between the two hatch dates for each sibling pair and creat a new df. In this way I can detemine an average hatch duration in-the-wild.
#sibHatchDiff<-

#added 22 July 2015. This can be refined to extract sibs as a list with hatch DoY and  look at days taken for sibs to hatch.In turn this can help ensure non sibs dont get mooshed together in the same nest if they are likely hatched too far apart, even if close geographically.