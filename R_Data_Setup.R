
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

###################################################
#Load integrated temperature dataset and region shape files-------------------------------
###################################################

#Same steps from Sam's Temperature QAQC R script

is.even <- function(x) as.integer(x) %% 2 == 0

# Load Delta Shapefile from Brian
Delta<-st_read(file.path(data_root,"Delta subregions","EDSM_Subregions_03302020.shp"))%>%
  filter(!SubRegion%in%c("South Bay", "San Francisco Bay", "San Pablo Bay", "Upper Yolo Bypass", 
                         "Upper Napa River", "Lower Napa River", "Carquinez Strait"))%>% # Remove regions outside our domain of interest
  dplyr::select(SubRegion)

# Load data
Data <- wq()%>%
  filter(!is.na(Temperature) & !is.na(Datetime) & !is.na(Latitude) & !is.na(Longitude) & !is.na(Date) & !is.na(Temperature_bottom))%>% #Remove any rows with NAs in our key variables
  filter(Temperature !=0)%>% #Remove 0 temps
  mutate(Temperature_bottom=if_else(Temperature_bottom>30, NA_real_, Temperature_bottom))%>% #Remove bad bottom temps
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
  mutate(Date_num = as.numeric(Date))%>%  # Create numeric version of date for models
  mutate(Time_num=as.numeric(Time)) # Create numeric version of time for models (=seconds since midnight)


# Pull station locations for major monitoring programs
# This will be used to set a boundary for this analysis focused on well-sampled regions.
WQ_stations<-Data%>%
  st_drop_geometry()%>%
  filter(Source%in%c("FMWT", "STN", "SKT", "20mm", "EMP", "Suisun"))%>%
  group_by(StationID, Source, Latitude, Longitude)%>%
  summarise(N=n(), .groups="drop")%>% # Calculate how many times each station was sampled
  filter(N>50 & !StationID%in%c("20mm 918", "STN 918"))%>% # Only keep stations sampled >50 times when deciding which regions to retain. 
  # "20mm 918", "STN 918" are far south of the rest of the well-sampled sites and are not sampled year round, so we're removing them to exclude that far southern region
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Convert to sf object
  st_transform(crs=st_crs(Delta))%>%
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
  dplyr::select(-IN)%>%
  mutate(Group=if_else(is.even(Year), 1, 2))%>%
  mutate_at(vars(Date_num, Longitude, Latitude, Time_num, Year, Julian_day), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T))) # Create centered and standardized versions of covariates

######## Begin editing Sam's data here ###########################
#Create new dataframe to keep the original data
Data_subset<- Data

#Add Water Year variable
Data_subset$WaterYear<-ifelse(month(Data_subset$Date)>=10,year(Data_subset$Date)+1,year(Data_subset$Date))

#Add UTM coordinates
UTM_coordinates<-sf::st_coordinates(Data_subset)
Data_subset$x<-UTM_coordinates[,c("X")]
Data_subset$y<-UTM_coordinates[,c("Y")]

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
summary(Data_subset_post_2011$Temperature)
#There were data points with bottom temp at 0, and another with surface temp at 33 C while bottom temp were listed as 23 (clearly a mistake)
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
#Only 1 point is outside the boundary, EMP site that's somehow located on land

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

########################################################
# cut out data points close to boundaries as and save as a different dataset

ggplot()+
  geom_sf(data=Delta.aut_sf)+
  geom_sf(data=Data_subset_post_2011_inside, color="red")+
  theme_bw()

distances<-Delta.aut_sf %>%
  st_cast(to = 'LINESTRING') %>% #TUrn polygon into linestring
  st_distance(y = Data_subset_post_2011_inside) # Get distance of each point from that perimeter linestring


#remove knots within 400 m of boundary
Data_subset_post_2011_inside_edited<-Data_subset_post_2011_inside[-which(distances<units::set_units(400, "m")),]

ggplot()+
  geom_sf(data=Delta.aut_sf)+
  geom_sf(data=Data_subset_post_2011_inside, color="red")+
  geom_sf(data=Data_subset_post_2011_inside_edited, color="blue")+
  theme_bw()


Data_subset_post_2011_inside_edited<-Data_subset_post_2011_inside_edited%>%
  st_drop_geometry()

###################################################
################# Export data and knots ############
###################################################
#Data points compared to knots

#Plot to show points and boundaries
plot(Delta.aut, col="grey")
points(Data_subset_post_2011_inside$x,Data_subset_post_2011_inside$y, pch=21, bg="blue")
points(knots_grid_reduced, pch=21, bg="white")
points(Data_subset_post_2011_outside$x,Data_subset_post_2011_outside$y, pch=21, bg="yellow")

#Export out data
Data_subset_post_2011_inside$geometry<-NULL
write.csv(Data_subset_post_2011_inside, file = "temperature_dataset.csv",row.names = F)
write.csv(knots_grid_reduced, file = "knots_grid.csv",row.names = F)
write.csv(Data_subset_post_2011_inside_edited, file = "temperature_dataset_edited.csv",row.names = F)









#####Extra Codes below



##########################################
#CODE TO COMPARE DATASET
##########################################
temp_dataset<-read.csv("temperature_dataset_old.csv")

library(compareDF)
library(htmlTable)

new<-as.data.frame(unique(Data_subset_post_2011_inside[,c("Source","Station","Date","Latitude","Longitude")]))
new$Date<-as.Date(new$Date)

old<-as.data.frame(unique(temp_dataset[,c("Source","Station","Date","Latitude","Longitude")]))
old$Date<-as.Date(old$Date)


plot(old$y~old$x, pch=21, bg="blue")
points(new$y~new$x, pch=21, bg="red")

ctable = compare_df(new, old,c("Latitude","Longitude"))

print(ctable$comparison_table_diff)
print(ctable$html_output)
ctable$change_count
ctable$change_summary
ctable$comparison_df


write.csv(ctable$comparison_df, file = "temperature_dataset_mismatch.csv",row.names = F)




tabletest<-anti_join(old, new)


plot(old$Latitude~old$Longitude, pch=21, bg="blue")
points(tabletest$Latitude~tabletest$Longitude, pch=21, bg="yellow")
