library(tidyverse)
library(lubridate)
library(broom)
library(rgdal)
library(mgcv)


# No need to specify absolute file paths in an Rstudio project (if you open it as a project).
data_root<-file.path("data-raw")

#Read in bay-Delta shape outline shape file that Mike Beakes created
Delta.aut <- readOGR(file.path(data_root,"Bay_Delta_Poly_Outline3_UTM10", "Bay_Delta_Poly_Outline3_UTM10.shp"))

###################################################
################# Set up Boundary and Knots ############
###################################################

Delta.xy.aut <- tidy(Delta.aut)
head(Delta.xy.aut)

deltacoords <- Delta.xy.aut %>% dplyr::select(long,lat,piece)
names(deltacoords) <- c("x", "y", "piece")
borderlist <- split(deltacoords, deltacoords$piece)
names(borderlist)

Delta.xy.aut <- Delta.xy.aut %>% rename(x = long, y = lat)
#Delta.xy.aut$piece

border.aut <- lapply(borderlist, "[", c(1,2))
nr <- seq(1,length(borderlist))

border.aut <- lapply(nr, function(n) as.list.data.frame(border.aut[[n]]))

#Use custom knots
knots <- read.csv("custom_knots.csv")

#Read dataset
temp_dataset<-read.csv("temperature_dataset.csv")
str(temp_dataset)

############### Fit test model  ######################
m_test <- bam(Temperature_bottom ~ s(x, y, k = 5, bs = "so", xt = list(bnd = border.aut,nmax=500)),
              data = temp_dataset, method = "fREML", discrete=T, nthreads=3, knots = knots)

knots_edit<-knots[-5,]

m_test <- bam(Temperature_bottom ~ s(x, y, k = 5, bs = "so", xt = list(bnd = border.aut,nmax=500)),
              data = temp_dataset, method = "fREML", discrete=T, nthreads=4, knots = knots_edit)

knots_edit<-knots_edit[-6,]
knots_edit<-knots_edit[-9,]
knots_edit<-knots_edit[-10,]
knots_edit<-knots_edit[-14,]
knots_edit<-knots_edit[-14,]
knots_edit<-knots_edit[-17,]



plot(Delta.aut, col="grey")
points(knots, pch=21, bg="red")
points(knots_edit, pch=21, bg="yellow")


#text(knots, labels=rownames(knots))