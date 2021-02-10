
#Install some github packages to pull Sam Bashevkin's integrated dataset

#library(devtools)
#devtools::install_github("sbashevkin/spacetools")
#install.packages("remotes")
#remotes::install_github("sbashevkin/deltareportr")
#devtools::install_github("sbashevkin/discretewq")

require(gstat)
require(sp)
require(spacetime)
library(tidyverse)
library(discretewq)
library(mgcv)
library(lubridate)
library(hms)
library(sf)
library(stars)
require(patchwork)
require(geofacet)
require(dtplyr)
require(scales)
require(rgdal)
library(broom)

#Set up a path for datasets and shapefiles
data_root<-file.path("data-raw")
results_root<-file.path("results")

###################################################
#Load integrated temperature dataset and region shape files-------------------------------
###################################################

# Similar steps taken from Sam's R script

# Load Delta Shapefile from Brian
Delta<-st_read(file.path(data_root,"Delta subregions","EDSM_Subregions_03302020.shp"))%>%
  filter(!SubRegion%in%c("South Bay", "San Francisco Bay", "San Pablo Bay", "Upper Yolo Bypass", 
                         "Upper Napa River", "Lower Napa River", "Carquinez Strait"))%>% # Remove regions outside our domain of interest
  dplyr::select(SubRegion)
# Visualize regions
ggplot(data=Delta,aes(label=SubRegion))+geom_sf()+geom_sf_text()

# Load data from 'discretewq' package
Data <- wq()%>%
  filter(!is.na(Temperature) & !is.na(Datetime) & !is.na(Latitude) & !is.na(Longitude) & !is.na(Date) & !is.na(Temperature_bottom))%>% #Remove any rows with NAs in our key variables
  filter(Temperature !=0)%>% #Remove 0 temps
  filter(hour(Datetime)>=5 & hour(Datetime)<=20)%>% # Only keep data between 5AM and 8PM
  mutate(Datetime = with_tz(Datetime, tz="America/Phoenix"), #Convert to a timezone without daylight savings time
         Date = with_tz(Date, tz="America/Phoenix"),
         Time=as_hms(Datetime), # Create variable for time-of-day, not date. 
         Noon_diff=abs(hms(hours=12)-Time))%>% # Calculate difference from noon for each data point for later filtering
  group_by(Station, Source, Date)%>%
  filter(Noon_diff==min(Noon_diff))%>% # Select only 1 data point per station and date, choose data closest to noon
  filter(Time==min(Time))%>% # When points are equidistant from noon, select earlier point
  ungroup()%>%
  distinct(Date, Station, Source, .keep_all = TRUE)%>% # Finally, remove the ~10 straggling datapoints from the same time and station
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Convert to sf object
  st_transform(crs=st_crs(Delta))%>% # Change to crs of Delta
  st_join(Delta, join=st_intersects)%>% # Add subregions
  filter(!is.na(SubRegion))%>% # Remove any data outside our subregions of interest
  mutate(Julian_day = yday(Date), # Create julian day variable
         Month_fac=factor(Month), # Create month factor variable
         Source_fac=factor(Source),
         Year_fac=factor(Year))%>% 
  mutate(Date_num = as.numeric(Date))%>%  # Create numeric version of date; keep just in case we need it
  mutate(Time_num=as.numeric(Time)) # Create numeric version of time (=seconds since midnight); keep just in case we need it


# Pull station locations for major monitoring programs
# This will be used to set a boundary for this analysis focused on well-sampled regions.
WQ_stations<-Data%>%
  st_drop_geometry()%>%
  filter(Year>=2011)%>% #Filter just year 2011 and after because of the issue noted below of limited spatial scope prior to 2011
  group_by(StationID, Source, Latitude, Longitude)%>%
  summarise(N=n(), .groups="drop")%>% # Calculate how many times each station was sampled
  filter(N>25 & !StationID%in%c("20mm 918", "STN 918"))%>% # Only keep stations sampled >25 times when deciding which regions to retain. 
  # "20mm 918", "STN 918" are far south of the rest of the well-sampled sites and are not sampled year round, so we're removing them to exclude that far southern region
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Convert to sf object
  st_transform(crs=st_crs(Delta))%>%
  st_join(Delta) # Add subregions

# Remove any subregions that do not contain at least one of these >25 samples stations from the major monitoring programs
Delta <- Delta%>%
  filter(SubRegion%in%unique(WQ_stations$SubRegion)) 

# Write out Delta shapefile for later use in making figures and plots
st_write(Delta,append = F, file.path(results_root,paste0( "SFE_regions", "Final.shp")))

# Save WQ station file for setting boundaries later on
saveRDS(WQ_stations, file.path(results_root,"WQ_Station_boundaries.Rds"))

# Visualize sampling regions of major surveys
ggplot(data=Delta,aes(label=SubRegion))+geom_sf()+geom_sf_text()

# Now filter data to only include this final set of subregions, and any stations outside the convex hull formed by the >25 samples stations from the major monitoring programs
Data<-Data%>%
  filter(SubRegion%in%unique(Delta$SubRegion))%>%
  st_join(WQ_stations%>%
            st_union()%>%
            st_convex_hull()%>% # Draws a hexagram or pentagram or similar around the outer-most points
            st_as_sf()%>%
            mutate(IN=TRUE),
          join=st_intersects)%>%
  filter(IN)%>%
  dplyr::select(-IN)%>%
  mutate_at(vars(Date_num, Longitude, Latitude, Time_num, Year, Julian_day), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T))) # Create centered and standardized versions of covariates

#######Sam's script ends here

###################################################
################# Setting Boundary ################
###################################################

#Read in bay-Delta shape outline shape file that Mike Beakes created
Delta.aut <- readOGR(file.path(data_root,"Bay_Delta_Poly_Outline3_UTM10", "Bay_Delta_Poly_Outline_NoSSC_UTM10.shp"))

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

######################################################################################################
################# Add a few more columns, inspect data and subset to years with good coverage ################
######################################################################################################

# Add Water Year variable
Data$WaterYear<-ifelse(month(Data$Date)>=10,year(Data$Date)+1,year(Data$Date))

# Make water year factor
Data$WaterYear<-as.factor(Data$WaterYear)

# Add UTM coordinates
UTM_coordinates<-sf::st_coordinates(Data)
Data$x<-UTM_coordinates[,c("X")]
Data$y<-UTM_coordinates[,c("Y")]

# Also subset data to 2011 and on---
# As they have the best spatiotemporal coverage
# Data from previous years only come from fixed stations sampled once a month with limited geographical distribution
# See below
Data_subset_pre_2011 <- Data %>% filter(year(Date)<2011)
str(Data_subset_pre_2011)
plot(Delta.aut, col="grey")
points(Data_subset_pre_2011$x,Data_subset_pre_2011$y, pch=21, bg="purple")

# Use data since 2011
Data_subset_post_2011 <- Data %>% filter(year(Date)>=2011)

# Remove temperatures that seem unreasonable for bottom and surface (those below 5 C and above 30 C)
summary(Data_subset_post_2011$Temperature_bottom)
summary(Data_subset_post_2011$Temperature)
# There were data points with bottom temp at 0, and another with surface temp at 33 C while bottom temp were listed as 23 (clearly a mistake)
Data_subset_post_2011<- Data_subset_post_2011 %>% filter(Temperature_bottom>=5&Temperature_bottom<=30)
Data_subset_post_2011<- Data_subset_post_2011 %>% filter(Temperature>=5&Temperature<=30)
hist(Data_subset_post_2011$Temperature_bottom)
hist(Data_subset_post_2011$Temperature)

# Keep only data points that are inside the boundaries
Data_subset_post_2011_inside <- Data_subset_post_2011[with(Data_subset_post_2011, inSide(bnd = border.aut, x, y)), ]
# Keep data points outside the boundaries just to see
Data_subset_post_2011_outside <- Data_subset_post_2011[!with(Data_subset_post_2011, inSide(bnd = border.aut, x, y)), ]

# Plot to show points and boundaries
plot(Delta.aut, col="grey")
points(Data_subset_post_2011_inside$x,Data_subset_post_2011_inside$y, pch=21, bg="purple")
points(Data_subset_post_2011_outside$x,Data_subset_post_2011_outside$y, pch=21, bg="red")
# Only 1 point is outside the boundary, EMP site that's somehow located on land

###################################################
################# Set up knots ############
###################################################

## Set up a different knots set up using 20x20 points
N <- 20
gx <- seq(min(Delta.xy.aut[,1]), max(Delta.xy.aut[,1]), len = N)
gy <- seq(min(Delta.xy.aut[,2]), max(Delta.xy.aut[,2]), len = N)
gp <- expand.grid(gx, gy)
names(gp) <- c("x","y")
knots_grid <- gp[with(gp, inSide(border.aut, x, y)), ]
row.names(knots_grid)<-1:nrow(knots_grid)

# Plot in map just to evaluate the spread of knots
plot(Delta.aut, col="grey")
points(knots_grid, pch=21, bg="orange")
text(knots_grid, labels=rownames(knots_grid))

# Editing out knots too close to the border with sf ----------------------------------------------------------
Delta.aut_sf<-st_as_sf(Delta.aut)

knots_grid_sf<-st_as_sf(knots_grid, coords=c("x","y"), crs=st_crs(Delta.aut_sf), remove=F)

ggplot()+
  geom_sf(data=Delta.aut_sf)+
  geom_sf(data=knots_grid_sf, color="red")+
  theme_bw()

distances<-Delta.aut_sf %>%
  st_cast(to = 'LINESTRING') %>% #TUrn polygon into linestring
  st_distance(y = knots_grid_sf) # Get distance of each point from that perimeter linestring

#remove knots within 500 m of boundary
knots_grid_sf_edited<-knots_grid_sf[-which(distances<units::set_units(500, "m")),]

ggplot()+
  geom_sf(data=Delta.aut_sf)+
  geom_sf(data=knots_grid_sf, color="red")+
  geom_sf(data=knots_grid_sf_edited, color="blue")+
  theme_bw()

knots_grid_reduced<-knots_grid_sf_edited%>%
  st_drop_geometry()


######################################################################################################
################# Calculate temperature difference and create temperature anomaly calculation ############
######################################################################################################
###  Dataset is essentially final at this point 

# Calculate difference of bottom temperature from surface temperature as response variable
Data_subset_post_2011_inside$Temperature_difference <- Data_subset_post_2011_inside$Temperature_bottom-Data_subset_post_2011_inside$Temperature
hist(Data_subset_post_2011_inside$Temperature_difference)

# Temperature Anomaly Model with just julian day
# Julian day and temperature model
temperature_anomaly_GAM<- gam(Temperature ~ s(Julian_day_s,bs="cc",k=5),data=Data_subset_post_2011_inside)
summary(temperature_anomaly_GAM)
gam.check(temperature_anomaly_GAM)
plot(temperature_anomaly_GAM)
# k=5 seems to make the most sense based on evaluating plot wiggliness and very small changes in R^2 as we increase K
# Model is not used any longer because we would like to adjust temperature based on space in addition to time (i.e., season)


## Construct temperature anomaly model with longitude and latitude on top of julian day
temperature_anomaly_GAM_spatial<- gam(Temperature ~ te(x,y,Julian_day_s, d=c(2,1) ,bs=c("tp","cc"),k=c(10,5)),data=Data_subset_post_2011_inside)
summary(temperature_anomaly_GAM_spatial)
# R-sq.(adj) =  0.907   Deviance explained = 90.7%

#gam.check(temperature_anomaly_GAM_spatial)

# K index is low at 0.52, and edf is pretty close to k', but the goal was to remove collinearity, not fit
# k'  edf k-index p-value    
# 39.0 35.8    0.52  <2e-16 ***

#Add temperature anomaly term to the dataset
Data_subset_post_2011_inside$Temperature_prediction_spatial <-predict(temperature_anomaly_GAM_spatial,Data_subset_post_2011_inside)
Data_subset_post_2011_inside$Temperature_anomaly_spatial <-Data_subset_post_2011_inside$Temperature - Data_subset_post_2011_inside$Temperature_prediction_spatial 
hist(Data_subset_post_2011_inside$Temperature_anomaly_spatial)

# Save model for later use in plotting results
saveRDS(temperature_anomaly_GAM_spatial, file.path(results_root,"Temperature_anomaly_spatial_GAM.Rds"))

###################################################
################# Export data and knots ############
###################################################
#Data points compared to knots

#Plot to show points and boundaries
plot(Delta.aut, col="grey")
points(Data_subset_post_2011_inside$x,Data_subset_post_2011_inside$y, pch=21, bg="blue")
points(knots_grid_reduced, pch=21, bg="white")
points(Data_subset_post_2011_outside$x,Data_subset_post_2011_outside$y, pch=21, bg="yellow")

# Save boundaries file
saveRDS(border.aut, file.path(results_root,"Soap_film_boundaries.Rds"))

# Export out knots
write.csv(knots_grid_reduced, file = "knots_grid.csv",row.names = F)

# Export out data
saveRDS(Data_subset_post_2011_inside, "temperature_dataset.Rds")

## If we want csv instead below
Data_subset_post_2011_inside$geometry<-NULL
write.csv(Data_subset_post_2011_inside, file = "temperature_dataset.csv",row.names = F)
