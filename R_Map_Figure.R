###################################################
#Script to create maps-------------------------------
###################################################

require(dplyr)
require(sf)
require(ggplot2)
require(maps)
require(ggspatial)
require(viridis)
require(deltamapr)
require(ggrepel)

library(raster)
library(rgdal)
library(rgeos)
library(maptools)

data_root<-file.path("data-raw")
results_root<-file.path("results")

###############################

#Read Delta subregions
Delta<-st_read(file.path(data_root,"Delta subregions","EDSM_Subregions_03302020.shp"))

#Load lat and long from final temperature dataset
temp_dataset <- readRDS("temperature_dataset.Rds")
#temp_dataset<- temp_dataset %>% st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)

#Load data for continuous water quality stations
continuous_stations<-read.csv(file.path(data_root,"Continuous depth station locations.csv"))
continuous_stations<- continuous_stations %>% st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)

#Calculate sample size
stations_sample_size<-temp_dataset %>% mutate(Sample=1) %>% group_by(Source,Station,Latitude,Longitude) %>% summarise(SampleSize=sum(Sample))

#Use code below to check that the numbers from stations_sample_size is correct
#latlong_unique<-unique(temp_dataset[c("Latitude","Longitude")])

#Load boundaries we set in the data set up r code
WQ_stations<-readRDS(file.path(results_root,"WQ_Station_boundaries.Rds"))

#Remove any subregions not used in analysis
Delta <- Delta%>%  st_transform(crs=4326)%>%
  filter(SubRegion%in%unique(WQ_stations$SubRegion)) 

#Delta centroids
Delta_centroids <- st_point_on_surface(Delta)
Delta_centroids_coords<-as.data.frame(st_coordinates(Delta_centroids))

Delta_centroids <- cbind(Delta_centroids,Delta_centroids_coords)

Delta_centroids$nudge_x<-0
Delta_centroids$nudge_y<-0

x_range <- abs(Reduce("-", range(Delta_centroids$X)))
y_range <- abs(Reduce("-", range(Delta_centroids$Y)))

unique(Delta_centroids$SubRegion)

Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Upper Sacramento River Ship Channel"] <- -1 * 0.30 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Upper Sacramento River Ship Channel"] <- 1 * 0.15 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Lower Sacramento River Ship Channel"] <- -1 * 0.30 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Lower Sacramento River Ship Channel"] <- 1 * 0.30 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Liberty Island"] <- -1 * 0.30 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Liberty Island"] <- 1 * 0.30 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Cache Slough and Lindsey Slough"] <- -1 * 0.32 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Cache Slough and Lindsey Slough"] <- 1 * 0.25 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Lower Cache Slough"] <- -1 * 0.40 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Lower Cache Slough"] <- 1 * 0.20 * y_range

Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Sacramento River near Rio Vista"] <- -1 * 0.50 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Sacramento River near Rio Vista"] <- 1 * 0.28 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Suisun Marsh"] <- -1 * 0.1 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Suisun Marsh"] <- 1 * 0.1 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Grizzly Bay"] <- -1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Grizzly Bay"] <- 1 * 0.1 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="West Suisun Bay"] <- -1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="West Suisun Bay"] <- 1 * 0.1 * y_range

Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Mid Suisun Bay"] <- -1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Mid Suisun Bay"] <- -1 * 0.1 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Honker Bay"] <- -1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Honker Bay"] <- -1 * 0.1 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Confluence"] <- -1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Confluence"] <- -1 * 0.1 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Lower Sacramento River"] <- -1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Lower Sacramento River"] <- -1 * 0.1 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Lower San Joaquin River"] <- -1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Lower San Joaquin River"] <- -1 * 0.1 * y_range

Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Franks Tract"] <- -1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Franks Tract"] <- -1 * 0.30 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Holland Cut"] <- -1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Holland Cut"] <- -1 * 0.30 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Old River"] <- -1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Old River"] <- -1 * 0.30 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Mildred Island"] <- -1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Mildred Island"] <- -1 * 0.30 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Victoria Canal"] <- -1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Victoria Canal"] <- -1 * 0.15 * y_range

Delta_centroids$nudge_x[Delta_centroids$SubRegion=="San Joaquin River at Twitchell Island"] <- 1 * 0.35 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="San Joaquin River at Twitchell Island"] <- 1 * 0.30 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="San Joaquin River at Prisoners Pt"] <- 1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="San Joaquin River at Prisoners Pt"] <- 1 * 0.30 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Upper Mokelumne River"] <- 1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Upper Mokelumne River"] <- 1 * 0.30 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Lower Mokelumne River"] <- 1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Lower Mokelumne River"] <- 1 * 0.30 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Disappointment Slough"] <- 1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Disappointment Slough"] <- 1 * 0.20 * y_range

Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Middle Sacramento River"] <- 1 * 0.25 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Middle Sacramento River"] <- 1 * 0.10 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Steamboat and Miner Slough"] <- 1 * 0.35 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Steamboat and Miner Slough"] <- 1 * 0.09 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Sacramento River near Ryde"] <- 1 * 0.30 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Sacramento River near Ryde"] <- 1 * 0.20 * y_range

Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Middle River"] <- 1 * 0.30 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Middle River"] <- -1 * 0.20 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="San Joaquin River near Stockton"] <- 1 * 0.30 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="San Joaquin River near Stockton"] <- -1 * 0.20 * y_range

#Read water boundaries shape file
Water<-st_read(file.path(data_root,"Shapefiles for Map Figure","Hydro_poly_UTM10Copy.shp"))




crsLONGLAT <- "+proj=longlat +datum=WGS84 +no_defs"

fig1<-ggplot() + theme_bw()+
  #geom_sf(data = Water, fill="slategray1", color="slategray2") +
  geom_sf(data = Water, fill="cadetblue1", color="cadetblue1") +
  geom_sf(data=stations_sample_size,shape=19, size=2,aes(color=SampleSize))+
  geom_sf(data=continuous_stations, fill="red", size=1.5, shape=24)+
  geom_sf(data = Delta,color="navy",fill=NA) + 
  geom_label_repel(data=continuous_stations, aes(x=Longitude,y=Latitude,label=StationName), nudge_x=-0.2,nudge_y=-0.2,segment.alpha=0.7,color="red4")+
  geom_text_repel(data=Delta_centroids, aes(x=X,y=Y,label=SubRegion), nudge_x = Delta_centroids$nudge_x, nudge_y = Delta_centroids$nudge_y, 
                  segment.linetype="dotted", segment.alpha=0.7, color="blue3", size=3)+
  coord_sf(xlim = c(-122.3, -121.16), ylim = c(37.7, 38.61),crs=crsLONGLAT)  +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_y = unit(1.0, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "tr", width_hint = 0.5)+
  #scale_fill_gradient(low="white", high="blue") +
  scale_color_viridis(discrete=FALSE,name="N")
  #theme(legend.position = "none") +
  theme(axis.text.x = element_text(size=16, color="black"),axis.text.y = element_text(size=16, color="black"),axis.title.x=element_blank(),axis.title.y=element_blank())



tiff(filename=file.path(results_root,"Figure01_Map.tiff"), units="in",type="cairo", bg="white", height=10, 
    width=10, res=300, compression="lzw")
fig1
dev.off()


#Load US states map
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))%>%
  st_transform(crs=st_crs(SubRegions))
california<-filter(states, ID=="california")


####################


base<-deltamapr::WW_Delta%>%
  st_transform(crs=st_crs(Delta))%>%
  st_crop(Delta)

plot(select(base, geometry),reset=F, col="slategray1", border="slategray2")

plot(select(Delta, geometry), add=T, lwd=2)
points(as_tibble(st_coordinates(stations_sample_size)), col="black", pch=16)

Letter_locs<-locator()




###############################

#Read Delta subregions
Delta_subregions <- readOGR(file.path(data_root,"Delta subregions","EDSM_Subregions_03302020.shp"))
Delta_subregions <- spTransform(Delta_subregions, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#Load boundaries we set in the data set up r code
WQ_stations<-readRDS(file.path(results_root,"WQ_Station_boundaries.Rds"))

#Remove any subregions not used in analysis
Delta_subregions <- Delta_subregions[which(Delta_subregions$SubRegion %in% unique(WQ_stations$SubRegion)),]

#Load Delta subregions with sf package to get centroid
Delta_subregions_st <- st_read(file.path(data_root,"Delta subregions","EDSM_Subregions_03302020.shp")) %>%
  st_transform(crs=4326)%>%
  filter(SubRegion%in%unique(WQ_stations$SubRegion)) 

#Get centroid
Delta_subregions_centroid <- data.frame(SubRegion = Delta_subregions_st$SubRegion,
                                        st_coordinates(st_point_on_surface(Delta_subregions_st)))

#Read water boundaries shape file
Water <- readOGR(file.path(data_root,"Shapefiles for Map Figure","Hydro_poly_UTM10Copy.shp"))
Water <- spTransform(Water, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))



plot(Delta_subregions)




############################################
##Figures for BDSC 2021 Presentation Only
############################################
#Create map for just continuous stations


fig_cont<-ggplot() + theme_bw()+
  #geom_sf(data = Water, fill="slategray1", color="slategray2") +
  geom_sf(data = Water, fill="cadetblue1", color="cadetblue1") +
  geom_sf(data = Delta, alpha=0.1) + 
  geom_label_repel(data=continuous_stations,size=7, aes(x=Longitude,y=Latitude,label=StationName),segment.alpha=0.7,color="black")+
  geom_sf(data=continuous_stations, fill="red", size=5, shape=24)+
  coord_sf(xlim = c(-122.2, -121.2), ylim = c(37.8, 38.61),crs=crsLONGLAT)  +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_y = unit(1.0, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "tr", width_hint = 0.5)+
  #scale_fill_gradient(low="white", high="blue") +
  #theme(legend.position = "none") +
  theme(axis.text.x = element_text(size=16, color="black"),axis.text.y = element_text(size=16, color="black"),axis.title.x=element_blank(),axis.title.y=element_blank())

fig_cont

tiff(filename=file.path(results_root,"Figure01_Map_Continuous.tiff"), units="in",type="cairo", bg="white", height=10, 
     width=10, res=300, compression="lzw")
fig_cont
dev.off()

#Just the 3 letter code
fig_cont_short_label<-ggplot() + theme_bw()+
  #geom_sf(data = Water, fill="slategray1", color="slategray2") +
  geom_sf(data = Water, fill="cadetblue1", color="cadetblue1") +
  geom_sf(data = Delta, alpha=0.1) + 
  geom_label_repel(data=continuous_stations,size=7, aes(x=Longitude,y=Latitude,label=Station),segment.alpha=0.7,color="black")+
  geom_sf(data=continuous_stations, fill="red", size=5, shape=24)+
  coord_sf(xlim = c(-122.2, -121.2), ylim = c(37.8, 38.61),crs=crsLONGLAT)  +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_y = unit(1.0, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "tr", width_hint = 0.5)+
  #scale_fill_gradient(low="white", high="blue") +
  #theme(legend.position = "none") +
  theme(axis.text.x = element_text(size=16, color="black"),axis.text.y = element_text(size=16, color="black"),axis.title.x=element_blank(),axis.title.y=element_blank())

fig_cont_short_label

tiff(filename=file.path(results_root,"Figure01_Map_Continuous_shortlabel.tiff"), units="in",type="cairo", bg="white", height=10, 
     width=10, res=300, compression="lzw")
fig_cont_short_label
dev.off()