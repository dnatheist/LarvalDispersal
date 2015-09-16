# Function to plot selection of examples otolith transects to illustrate core edge, and some elements etc. 

#elementList<-list("Mn","Ba","Sr")

#LarvalIDList<- list(186,187,188,189,190,191)

#lapply(LarvalIDList,otoTransectPlot)
        
otoTransectPlot <- function(LarvalIDList){
        
        temp<-subset(allOtoChemData, LarvalID==LarvalIDList, select=c(Time,Mn,Ba,Sr,OtoPart))
        head(temp)
        p <- ggplot(temp, aes(Time, Mn)) + ggtitle("Otolith Ablation Transect for ")+ labs(colour = "Otolith Part")+ xlab("Ablation Time (ms)")+ ylab("Concentration(ug/mg)")
        p + geom_point()
        p + geom_point(aes(colour = factor(OtoPart)))
        
        #rm(temp)
        }