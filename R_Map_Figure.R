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


data_root<-file.path("data-raw")
results_root<-file.path("results")

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
Delta <- Delta%>%
  filter(SubRegion%in%unique(WQ_stations$SubRegion)) 

#Read water boundaries shape file
Water<-st_read(file.path(data_root,"Shapefiles for Map Figure","Hydro_poly_UTM10Copy.shp"))

#Load US states map
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))%>%
  st_transform(crs=st_crs(SubRegions))
california<-filter(states, ID=="california")



crsLONGLAT <- "+proj=longlat +datum=WGS84 +no_defs"

fig1<-ggplot() + theme_bw()+
  #geom_sf(data = Water, fill="slategray1", color="slategray2") +
  geom_sf(data = Water, fill="cadetblue1", color="cadetblue1") +
  geom_sf(data = Delta, alpha=0.1) + 
  geom_sf(data=stations_sample_size, color="black",shape=21, size=2,aes(fill=SampleSize))+
  #geom_sf(data=continuous_stations, fill="red", size=1.5, shape=24)+
  coord_sf(xlim = c(-122.2, -121.2), ylim = c(37.8, 38.61),crs=crsLONGLAT)  +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_y = unit(1.0, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "tr", width_hint = 0.5)+
  #scale_fill_gradient(low="white", high="blue") +
  scale_fill_viridis(discrete=FALSE,name="N")
  #theme(legend.position = "none") +
  theme(axis.text.x = element_text(size=16, color="black"),axis.text.y = element_text(size=16, color="black"),axis.title.x=element_blank(),axis.title.y=element_blank())

fig1


tiff(filename=file.path(results_root,"Figure01_Map.tiff"), units="in",type="cairo", bg="white", height=10, 
    width=10, res=300, compression="lzw")
fig1
dev.off()



base<-deltamapr::WW_Delta%>%
  st_transform(crs=st_crs(Delta))%>%
  st_crop(Delta)

plot(select(base, geometry),reset=F, col="slategray1", border="slategray2")
plot(select(Delta, geometry), add=T, lwd=2)
points(as_tibble(st_coordinates(stations_sample_size)), col="black", pch=16)

Letter_locs<-locator()






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