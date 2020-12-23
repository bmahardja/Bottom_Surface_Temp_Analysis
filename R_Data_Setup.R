
#Install some github packages to pull Sam Bashevkin's integrated dataset

#library(devtools)
#devtools::install_github("sbashevkin/spacetools")
#install.packages("remotes")
#remotes::install_github("sbashevkin/deltareportr")


library(tidyverse)
library(deltareportr)
library(lubridate)
library(hms)
library(rgdal)
library(sf)
library(broom)
library(mgcv)
library(units)
library(rgeos)

#Set up a path for datasets and shapefiles
data_root<-file.path("data-raw")

###################################################
################# Setting Boundary ############
###################################################

#Read in bay-Delta shape outline shape file that Mike Beakes created
Delta.aut <- readOGR(file.path(data_root,"Bay_Delta_Poly_Outline3_UTM10", "Bay_Delta_Poly_Outline3_UTM10.shp"))

Delta.xy.aut <- tidy(Delta.aut)
head(Delta.xy.aut)

deltacoords <- Delta.xy.aut %>% dplyr::select(long,lat,piece)
names(deltacoords) <- c("x", "y", "piece")
borderlist <- split(deltacoords, deltacoords$piece)
names(borderlist)

Delta.xy.aut <- Delta.xy.aut %>% rename(x = long, y = lat)

border.aut <- lapply(borderlist, "[", c(1,2))
nr <- seq(1,length(borderlist))

border.aut <- lapply(nr, function(n) as.list.data.frame(border.aut[[n]]))

###################################################
#Load integrated temperature dataset and region shape files-------------------------------
###################################################

#Same steps from Sam's Temperature QAQC R script

# Load Delta Shapefile from Brian
Delta<-st_read(file.path(data_root,"Delta Subregions"))%>%
  filter(!SubRegion%in%c("South Bay", "San Francisco Bay", "San Pablo Bay", "Upper Yolo Bypass", 
                         "Upper Napa River", "Lower Napa River", "Carquinez Strait"))%>% # Remove regions outside our domain of interest
  dplyr::select(SubRegion)

# Load data
Data <- DeltaDater(Start_year = 1900, 
                   WQ_sources = c("EMP", "STN", "FMWT", "EDSM", "DJFMP", "SKT", "20mm", "Suisun", "Baystudy", "USBR", "USGS"), 
                   Variables = "Water quality", 
                   Regions = NULL)%>%
  filter(!is.na(Temperature) & !is.na(Datetime) & !is.na(Latitude) & !is.na(Longitude) & !is.na(Date))%>% #Remove any rows with NAs in our key variables
  filter(Temperature !=0)%>% #Remove 0 temps
  mutate(Temperature_bottom=if_else(Temperature_bottom>30, NA_real_, Temperature_bottom))%>% #Remove bad bottom temps
  filter(hour(Datetime)>=5 & hour(Datetime)<=20)%>% # Only keep data betwen 5AM and 8PM
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Convert to sf object
  st_transform(crs=st_crs(Delta))%>% # Change to crs of Delta
  st_join(Delta, join=st_intersects)%>% # Add subregions
  filter(!is.na(SubRegion))%>% # Remove any data outside our subregions of interest
  mutate(Datetime = with_tz(Datetime, tz="America/Phoenix"), #Convert to a timezone without daylight savings time
         Date = with_tz(Date, tz="America/Phoenix"),
         Julian_day = yday(Date), # Create julian day variable
         Month_fac=factor(Month), # Create month factor variable
         Source_fac=factor(Source),
         Year_fac=factor(Year))%>% #BM: Changed from Sam's original code to make it non-ordered
  mutate(Date_num = as.numeric(Date), # Create numeric version of date for models
         Time = as_hms(Datetime))%>% # Create variable for time-of-day, not date. 
  mutate(Time_num=as.numeric(Time))%>% # Create numeric version of time for models (=seconds since midnight)
  mutate_at(vars(Date_num, Longitude, Latitude, Time_num, Year, Julian_day), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T))) # Create centered and standardized versions of covariates

# Pull station locations for major monitoring programs
# This will be used to set a boundary for this analysis focused on well-sampled regions.
WQ_stations<-Data%>%
  filter(Source%in%c("FMWT", "STN", "SKT", "20mm", "EMP", "Suisun"))%>%
  group_by(StationID, Source)%>%
  summarise(N=n())%>% # Calculate how many times each station was sampled
  filter(N>50 & !StationID%in%c("20mm 918", "STN 918"))%>% # Only keep stations sampled >50 times when deciding which regions to retain. 
  # "20mm 918", "STN 918" are far south of the rest of the well-sampled sites and are not sampled year round, so we're removing them to exclude that far southern region
  st_join(Delta) # Add subregions

# Remove any subregions that do not contain at least one of these >50 samples stations from the major monitoring programs
Delta <- Delta%>%
  filter(SubRegion%in%unique(WQ_stations$SubRegion) | SubRegion=="Georgiana Slough") # Retain Georgiana Slough because it's surrounded by well-sampled regions
# Visualize sampling regions of major surveys

# Now filter data to only include this final set of subregions, and any stations outside the convex hull formed by the >50 samples stations from the major monitoring programs
Data<-Data%>%
  filter(SubRegion%in%unique(Delta$SubRegion))%>%
  st_join(WQ_stations%>%
            st_union()%>%
            st_convex_hull()%>% # Draws a hexagram or pentagram or similar around the outer-most points
            st_as_sf()%>%
            mutate(IN=TRUE),
          join=st_intersects)%>%
  filter(IN)%>%
  dplyr::select(-IN)

#Filter just those that have bottom temperature measurements
Data_subset<- Data %>% filter(!is.na(Temperature_bottom))

#Convert data to UTM since lat and long don't seem to work well with soap-film
cord.dec <- SpatialPoints(cbind(Data_subset$Longitude, -Data_subset$Latitude), proj4string = CRS("+proj=longlat +datum=WGS84"))
cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=10 +datum=WGS84"))
cord.UTM.data.frame <-as.data.frame(cord.UTM)
cord.UTM.data.frame$coords.x2<-abs(cord.UTM.data.frame$coords.x2)

Data_subset$x<-cord.UTM.data.frame$coords.x1
Data_subset$y<-cord.UTM.data.frame$coords.x2

#Also subset data to 2011 and on---
#As they have the best spatiotemporal coverage
#Data from previous years only come from fixed stations sampled once a month with limited geographical distribution
#See below
Data_subset_pre_2011 <- Data_subset %>% filter(year(Date)<2011)
str(Data_subset_pre_2011)
plot(Delta.aut, col="grey")
points(Data_subset_pre_2011$x,Data_subset_pre_2011$y, pch=21, bg="purple")

Data_subset_post_2011 <- Data_subset %>% filter(year(Date)>=2011)


#Remove temperatures that seem unreasonable for bottom and surface (those below 5 C and above 30 C)
summary(Data_subset_post_2011$Temperature_bottom)

Data_subset_post_2011<- Data_subset_post_2011 %>% filter(Temperature_bottom>=5&Temperature_bottom<=30)
Data_subset_post_2011<- Data_subset_post_2011 %>% filter(Temperature>=5&Temperature<=30)
hist(Data_subset_post_2011$Temperature_bottom)
hist(Data_subset_post_2011$Temperature)

#Keep only data points that are inside the boundaries
Data_subset_post_2011_inside <- Data_subset_post_2011[with(Data_subset_post_2011, inSide(bnd = border.aut, x, y)), ]
#Keep data points outside the boundaries just to see
Data_subset_post_2011_outside <- Data_subset_post_2011[!with(Data_subset_post_2011, inSide(bnd = border.aut, x, y)), ]

#Plot to show points and boundaries
plot(Delta.aut, col="grey")
points(Data_subset_post_2011_inside$x,Data_subset_post_2011_inside$y, pch=21, bg="purple")
points(Data_subset_post_2011_outside$x,Data_subset_post_2011_outside$y, pch=21, bg="red")
#Only 51 points are outside the boundaries

###################################################
################# Set up knots ############
###################################################

river_km_points <- readOGR(file.path(data_root,"Bay_Delta_Poly_Outline3_UTM10", "River_km_UTM10.shp"))

knots_interior<-data.frame(river_km_points)
knots_interior<-knots_interior %>% dplyr::select(coords.x1,coords.x2) %>%rename(x=coords.x1,y=coords.x2)

#Reduce knot number by 1/4
Nth.delete<-function(dataframe, n)dataframe[-(seq(n,to=nrow(dataframe),by=n)),]
knots_interior_reduced<-Nth.delete(knots_interior,2)
knots_interior_reduced<-Nth.delete(knots_interior_reduced,2)
knots_interior_reduced<-Nth.delete(knots_interior_reduced,2)

#Check outline and compare to interior knots
plot(Delta.aut, col="grey")
#points(knots_interior, pch=21, bg="yellow")
points(knots_interior_reduced, pch=21, bg="red")

# Editing out knots too close to the border with sf ----------------------------------------------------------
Delta.aut_sf<-st_as_sf(Delta.aut)
knots_sf<-st_as_sf(knots_interior_reduced, coords=c("x","y"), crs=st_crs(Delta.aut_sf), remove=F)

ggplot()+
  geom_sf(data=Delta.aut_sf)+
  geom_sf(data=knots_sf, color="red")+
  theme_bw()

distances<-Delta.aut_sf %>%
  st_cast(to = 'LINESTRING') %>% #TUrn polygon into linestring
  st_distance(y = knots_sf) # Get distance of each point from that perimeter linestring

#remove knots within 10 m of boundary
knots_sf_edited<-knots_sf[-which(distances<units::set_units(10, "m")),]

ggplot()+
  geom_sf(data=Delta.aut_sf)+
  geom_sf(data=knots_sf, color="red")+
  geom_sf(data=knots_sf_edited, color="blue")+
  theme_bw()

knots_interior_reduced<-knots_sf_edited%>%
  st_drop_geometry()


###################################################
################# Export data and knots ############
###################################################
#Data points compared to knots

#Plot to show points and boundaries
plot(Delta.aut, col="grey")
points(Data_subset_post_2011_inside$x,Data_subset_post_2011_inside$y, pch=21, bg="blue")
points(knots_interior_reduced, pch=21, bg="red")

Data_subset_post_2011_inside$geometry<-NULL
write.csv(Data_subset_post_2011_inside, file = "temperature_dataset.csv",row.names = F)
write.csv(knots_interior_reduced, file = "custom_knots.csv",row.names = F)
