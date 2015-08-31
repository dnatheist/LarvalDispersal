####################################################################
## line_split.R
## R version 3.2.0 (2015-04-16) -- "Full of Ingredients)
## M. Cracknell
## 28/08/2015
##############################################################

# global settings
# make packages avaliable
# spatial data packages
library(sp)

# FUNCTIONS
#########################################################################
## split lines into segements of equal length (except final segment)

## modified from http://simonchamaille.net/sample-line-at-regular-interval/
## sample every distance interval along a line and get previous (or equal) node
## param cc a two-column matrix of x and y coords for the line
## param dist distance separating sample points
## returns coordinates of points (x= column 1, y = column 2) along the line separated by the given distance measured along the line
## also returns original line segments that point intersects (int = column 3)

sampleEvery <- function(cc, dist){
	## get line lengths
	lengths <- LineLength(cc, longlat = FALSE, sum = FALSE)
	## check precision
	if(any(abs(lengths) < .Machine$double.eps)) {
		wl <- which(abs(lengths) < .Machine$double.eps)
		cc <- cc[-(wl), ]
		lengths <- lengths[-(wl)]
	}
	## cumulative sum of lengths between line nodes
	csl <- c(0, cumsum(lengths))
	maxl <- csl[length(csl)]
	## point at length
	pts <- seq(0,sum(lengths),by=dist) #I am not interested in a sequence per se, but multiple, individual 
	                                   #distances.
	## get intervals
	int <- findInterval(pts, csl, all.inside = TRUE)
	## points intersect which line
	where <- (pts - csl[int])/diff(csl)[int]
	## get coords
	xy = cc[int, ,drop = FALSE] + where * (cc[int + 1, , drop = FALSE] - cc[int, , drop = FALSE])
	if (nrow(xy) < 1)
		return(NULL)
	return(cbind(xy,int))
}

####################################

## example

## coordinates of river meander line as a matrix
inCoords <-matrix(c(149.05563673736700,	-35.45112927380670,
                    149.05548650407800,	-35.45700190055300,
                    149.05887326947400,	-35.45963843202770,
                    149.06141632262600,	-35.46089587788500,
                    149.06439833043500,	-35.46344348678850,
                    149.06368945289200,	-35.46644646446390,
                    149.06212831006600,	-35.46921683868200,
                    149.05987352268900,	-35.47361586901360,
                    149.06001685303500,	-35.47511725395630,
                    149.06298218980300,	-35.47834293008770,
                    149.06638358320400,	-35.48099798038850,
                    149.06680909539600,	-35.48411619363670
                   ),ncol=2,byrow=TRUE)

#inCoords<-as.data.frame(inCoords)
## check and (NB swap x-y for lat-long)
plot(inCoords,xlim=c(149.055,149.07),ylim=c(-35.44,-35.5))
lines(inCoords)

## calculate points at distance dis
dis <- 0.03 # This is distance down the river from the start point (first entry in matrix). Ideally I'd like this to be a vector of n Distances. Can't be 'dist' as it is a reserved word?

## But The Earth is not a Sphere - it is an "Oblate Spheroid - it is 134.397 Km further around the Equator than it is around the Poles. So while 1 degree latitude ~ 110.6km, 1 degree longitude at -35S is about 91.2 km. So the distance in km needs to be converted to a 'decimal degree distance before use here' so does that mean it = km/91.2 for longitude component, and km/110.6 for latitude component?

pointsOut <- sampleEvery(inCoords,dis)

## plot
points(pointsOut[,-3],col="red",pch=20,cex=0.7)

#############################################################
