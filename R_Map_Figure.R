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
##Main Map Figure
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

#Read water boundaries shape file
Water<-st_read(file.path(data_root,"Shapefiles for Map Figure","Hydro_poly_UTM10Copy.shp"))

#Add crs
crsLONGLAT <- "+proj=longlat +datum=WGS84 +no_defs"

#Delta centroids
Delta_centroids <- st_point_on_surface(Delta)
Delta_centroids_coords<-as.data.frame(st_coordinates(Delta_centroids))
Delta_centroids <- cbind(Delta_centroids,Delta_centroids_coords)

#Add nudge data to customize label and text_repel for subregions
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
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Lower Sacramento River"] <- -1 * 0.25 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Lower Sacramento River"] <- -1 * 0.1 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Lower San Joaquin River"] <- -1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Lower San Joaquin River"] <- -1 * 0.1 * y_range

Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Franks Tract"] <- -1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Franks Tract"] <- -1 * 0.13 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="San Joaquin River at Prisoners Pt"] <- -1 * 0.38 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="San Joaquin River at Prisoners Pt"] <- -1 * 0.38 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Holland Cut"] <- -1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Holland Cut"] <- -1 * 0.30 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Old River"] <- -1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Old River"] <- -1 * 0.30 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Victoria Canal"] <- 1 * 0.15 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Victoria Canal"] <- -1 * 0.15 * y_range

Delta_centroids$nudge_x[Delta_centroids$SubRegion=="San Joaquin River at Twitchell Island"] <- 1 * 0.50 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="San Joaquin River at Twitchell Island"] <- 1 * 0.3 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Upper Mokelumne River"] <- 1 * 0.25 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Upper Mokelumne River"] <- 1 * 0.05 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Lower Mokelumne River"] <- 1 * 0.33 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Lower Mokelumne River"] <- 1 * 0.14 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Disappointment Slough"] <- 1 * 0.25 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Disappointment Slough"] <- 1 * 0.20 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="San Joaquin River near Stockton"] <- 1 * 0.30 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="San Joaquin River near Stockton"] <- 1 * 0.15 * y_range

Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Middle Sacramento River"] <- 1 * 0.3 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Middle Sacramento River"] <- 1 * 0.1 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Steamboat and Miner Slough"] <- 1 * 0.5 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Steamboat and Miner Slough"] <- 1 * 0.1 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Sacramento River near Ryde"] <- 1 * 0.40 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Sacramento River near Ryde"] <- 1 * 0.20 * y_range

Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Mildred Island"] <- 1 * 0.30 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Mildred Island"] <- -1 * 0.15 * y_range
Delta_centroids$nudge_x[Delta_centroids$SubRegion=="Middle River"] <- 1 * 0.30 * x_range
Delta_centroids$nudge_y[Delta_centroids$SubRegion=="Middle River"] <- -1 * 0.20 * y_range

#Adjust label nudge distance for continuous stations
continuous_stations$nudge_x<-0
continuous_stations$nudge_y<-0

continuous_stations$nudge_x[continuous_stations$Station=="MRZ"] <- -1 * 0.2 * x_range
continuous_stations$nudge_y[continuous_stations$Station=="MRZ"] <- -1 * 0.05 * y_range
continuous_stations$nudge_x[continuous_stations$Station=="MAL"] <- -1 * 0.0 * x_range
continuous_stations$nudge_y[continuous_stations$Station=="MAL"] <- -1 * 0.18 * y_range
continuous_stations$nudge_x[continuous_stations$Station=="ANH"] <- -1 * 0.0 * x_range
continuous_stations$nudge_y[continuous_stations$Station=="ANH"] <- -1 * 0.15 * y_range
continuous_stations$nudge_x[continuous_stations$Station=="RRI"] <- 1 * 0.25 * x_range
continuous_stations$nudge_y[continuous_stations$Station=="RRI"] <- -1 * 0.04 * y_range

#Create the map
fig1<-ggplot() + theme_bw()+
  geom_sf(data = Water, fill="cadetblue1", color="cadetblue1") +
  geom_sf(data=stations_sample_size,shape=19, size=2,aes(color=SampleSize))+
  geom_sf(data=continuous_stations, fill="red", size=2.2, shape=24)+
  geom_sf(data = Delta,color="navy",fill=NA) + 
  geom_label_repel(data=continuous_stations, aes(x=Longitude,y=Latitude,label=StationName), nudge_x=continuous_stations$nudge_x, nudge_y=continuous_stations$nudge_y
                   ,segment.alpha=0.7,color="red4", size=3.7,segment.linetype="dashed")+
  geom_text_repel(data=Delta_centroids, aes(x=X,y=Y,label=SubRegion), nudge_x = Delta_centroids$nudge_x, nudge_y = Delta_centroids$nudge_y, 
                  segment.linetype="dotted", segment.alpha=0.7, color="blue4", size=3.3)+
  coord_sf(xlim = c(-122.3, -121.10), ylim = c(37.65, 38.61),crs=crsLONGLAT)  +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_y = unit(1.0, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "tr", width_hint = 0.5)+
  #annotate(geom = "point", x = -121.81, y = 37.7, colour = "black", fill="red", size = 2.2,shape=24) + 
  annotate(geom = "text", x = -122.22, y = 37.85, label="San Francisco Bay",size=3.7,hjust="left") + 
  geom_segment(data=tibble(x=-122.22, y=37.85, xend=-122.3, yend=37.85), aes(x=x, y=y, xend=xend, yend=yend), arrow=arrow(length = unit(0.01, "npc")), size=0.5)+
  theme(legend.position = c(0.1,0.9),legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(size=12, color="black"),axis.text.y = element_text(size=12, color="black"),
        axis.title.x=element_blank(),axis.title.y=element_blank()
  )+
  scale_color_viridis(discrete=FALSE,name="Sample size",direction=-1, breaks=c(5, 25, 45, 65, 85))+
  guides(colour=guide_colourbar(ticks.colour = "black"))


#Print out the map
tiff(filename=file.path(results_root,"Figure01_Map.tiff"), units="in",type="cairo", bg="white", height=10, 
    width=11, res=300, compression="lzw")
fig1
dev.off()



###############################
##Supplementary Information - Map for soap-film smoother and knots
###############################

#Border outline for soap-film smoother
#Read in bay-Delta shape outline shape file that Mike Beakes created
Delta.aut <- st_read(file.path(data_root,"Bay_Delta_Poly_Outline3_UTM10", "Bay_Delta_Poly_Outline_NoSSC_UTM10.shp"))

#Add knots
knots_grid <- read.csv("knots_grid.csv")
knots_grid_sf<-st_as_sf(knots_grid, coords=c("x","y"), crs=st_crs(Water), remove=F)

#Create the map
fig_sup<-ggplot() + theme_bw()+
  geom_sf(data = Water, fill="cadetblue1", color="cadetblue1") +
  geom_sf(data=knots_grid_sf,shape=19, size=2,color="red")+
  geom_sf(data=Delta.aut,fill=NA)+
  coord_sf(xlim = c(-122.2, -121.20), ylim = c(37.65, 38.61),crs=crsLONGLAT)  +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_y = unit(1.0, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "tr", width_hint = 0.5)+
  #annotate(geom = "point", x = -121.81, y = 37.7, colour = "black", fill="red", size = 2.2,shape=24) + 
  #annotate(geom = "text", x = -122.22, y = 37.85, label="San Francisco Bay",size=3.5,hjust="left") + 
  geom_segment(data=tibble(x=-122.22, y=37.85, xend=-122.3, yend=37.85), aes(x=x, y=y, xend=xend, yend=yend), arrow=arrow(length = unit(0.01, "npc")), size=0.5)+
  theme(legend.position = c(0.1,0.9),legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(size=12, color="black"),axis.text.y = element_text(size=12, color="black"),
        axis.title.x=element_blank(),axis.title.y=element_blank()
  )+
  guides(colour=guide_colourbar(ticks.colour = "black"))


#Print out the map
tiff(filename=file.path(results_root,"FigureSupplementary_BoundaryKnotMap.tiff"), units="in",type="cairo", bg="white", height=10, 
     width=8, res=300, compression="lzw")
fig_sup
dev.off()




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





############################################
############################################
############################################
#Unused code below
############################################


#Load US states map
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))%>%
  st_transform(crs=st_crs(SubRegions))
california<-filter(states, ID=="california")


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


