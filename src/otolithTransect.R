# Plot an example of one otolith transect to illustrate core edge, etc. 
temp<-subset(allOtoChemData, LarvalID==186, select=c(Time,Mn,Ba,Sr,OtoPart))

p <- ggplot(temp, aes(Time, Mn)) + ggtitle("Otolith Ablation Transect for Manganese") + labs(colour = "Otolith Part")+ xlab("Ablation Time (ms)")+ ylab("Mn Concentration(ug/mg)")
p + geom_point()
p + geom_point(aes(colour = factor(OtoPart)))

# Clean Up
rm(temp)
rm(p)

##############scratchpad below

 



mynei %>%  #DPLYR required
        select(Emissions, year) %>%
        group_by(year) %>%
        summarise (total=sum(Emissions))%>%
        with(plot(year, total))

      
        plotTransect<-function(test2)
{
        p <- ggplot(test2, aes(Time, Mn)) + ggtitle("Otolith Ablation Transect for Manganese") + labs(colour = "Otolith Part")+ xlab("Ablation Time (ms)")+ ylab("Mn Concentration(ug/mg)")
        p+ geom_point()
        p + geom_point(aes(colour = factor(OtoPart)))
}

elList<-names(allOtoChemData[37:68])
larvIDs<-unique(allOtoChemData$LarvalID)
sapply(elList, nrow(allOtoChemData$elList))
subset.files<-subset(allOtoChemData)
data <- lapply(files, load.file)
names(data) <- elements


sapply(clist,nrow)

#p <- ggplot(temp, aes(Time, Mn)) + ggtitle("Otolith Ablation Transect for Manganese") + labs(colour = "Otolith Part")+ xlab("Ablation Time (ms)")+ ylab("Mn Concentration(ug/mg)")
#p + geom_point()
#p + geom_point(aes(colour = factor(OtoPart)))
#
#
testplot <- function(meansdf)
{
        p <- ggplot(meansdf, aes(fill=meansdf$condition, y=meansdf$means, x = meansdf$condition))
        p + geom_bar(position="dodge", stat="identity")
}
