library(png)

#Replace the directory and file information with your info
ima <- readPNG("C:\\Users\\s428825\\Google Drive\\PhD Analyses\\LarvalDispersal\\images\\MurrProfile.png")

#Set up the plot area
plot(1:2, type='n', main="Murrumbidgee River Profile with Larval Clades")

#Get the plot information so the image will fill the plot box, and draw it
lim <- par()
rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
#grid()
points(1.2,1.6)

