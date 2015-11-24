# some ggplot error bar examples. 
tg <- ToothGrowth
head(tg)

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
        library(plyr)
        
        # New version of length which can handle NA's: if na.rm==T, don't count them
        length2 <- function (x, na.rm=FALSE) {
                if (na.rm) sum(!is.na(x))
                else       length(x)
        }
        
        # This does the summary. For each group's data frame, return a vector with
        # N, mean, and sd
        datac <- ddply(data, groupvars, .drop=.drop,
                       .fun = function(xx, col) {
                               c(N    = length2(xx[[col]], na.rm=na.rm),
                                 mean = mean   (xx[[col]], na.rm=na.rm),
                                 sd   = sd     (xx[[col]], na.rm=na.rm)
                               )
                       },
                       measurevar
        )
        
        # Rename the "mean" column    
        datac <- rename(datac, c("mean" = measurevar))
        
        datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
        
        # Confidence interval multiplier for standard error
        # Calculate t-statistic for confidence interval: 
        # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
        ciMult <- qt(conf.interval/2 + .5, datac$N-1)
        datac$ci <- datac$se * ciMult
        
        return(datac)
}


tgc <- summarySE(tg, measurevar="len", groupvars=c("supp","dose"))
tgc

# Standard error of the mean
ggplot(tgc, aes(x=dose, y=len, colour=supp)) + 
        geom_errorbar(aes(ymin=len-se, ymax=len+se), width=.1) +
        geom_line() +
        geom_point()


# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right

ggplot(tgc, aes(x=dose, y=len, colour=supp)) + 
        geom_errorbar(aes(ymin=len-se, ymax=len+se), width=.1, position=pd) +
        geom_line(position=pd) +
        geom_point(position=pd)


# Use 95% confidence interval instead of SEM
ggplot(tgc, aes(x=dose, y=len, colour=supp)) + 
        geom_errorbar(aes(ymin=len-ci, ymax=len+ci), width=.1, position=pd) +
        geom_line(position=pd) +
        geom_point(position=pd)

# Black error bars - notice the mapping of 'group=supp' -- without it, the error
# bars won't be dodged!
ggplot(tgc, aes(x=dose, y=len, colour=supp, group=supp)) + 
        geom_errorbar(aes(ymin=len-ci, ymax=len+ci), colour="black", width=.1, position=pd) +
        geom_line(position=pd) +
        geom_point(position=pd, size=3)


ggplot(tgc, aes(x=dose, y=len, colour=supp, group=supp)) + 
        geom_errorbar(aes(ymin=len-se, ymax=len+se), colour="black", width=.1, position=pd) +
        geom_line(position=pd) +
        geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
        xlab("Dose (mg)") +
        ylab("Tooth length") +
        scale_colour_hue(name="Supplement type",    # Legend label, use darker colors
                         breaks=c("OJ", "VC"),
                         labels=c("Orange juice", "Ascorbic acid"),
                         l=40) +                    # Use darker colors, lightness=40
        ggtitle("The Effect of Vitamin C on\nTooth Growth in Guinea Pigs") +
        expand_limits(y=0) +                        # Expand y range
        scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
        theme_bw() +
        theme(legend.justification=c(1,0),
              legend.position=c(1,0))               # Position legend in bottom right
