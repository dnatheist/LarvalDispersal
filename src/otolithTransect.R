# Plot an example of one otolith transect to illustrate core edge, etc. 
temp<-subset(allOtoChemData, LarvalID==186, select=c(Time,Mn,Ba,Sr,OtoPart))

p <- ggplot(temp, aes(Time, Mn)) + ggtitle("Otolith Ablation Transect for Manganese") + labs(colour = "Otolith Part")+ xlab("Ablation Time (ms)")+ ylab("Mn Concentration(ug/mg)")
p + geom_point()
p + geom_point(aes(colour = factor(OtoPart)))

# Clean Up
rm(temp)
rm(p)