require(dplyr)
require(sf)
require(ggplot2)
require(maps)
require(ggspatial)

data_root<-file.path("data-raw")

Delta<-st_read(file.path(data_root,"Delta subregions","EDSM_Subregions_03302020.shp"))

#Load lat and long from final temperature dataset
temp_dataset<-read.csv("temperature_dataset.csv")
temp_dataset<- temp_dataset %>% st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)

#Load data for continuous water quality stations
continuous_stations<-read.csv(file.path(data_root,"Continuous depth station locations.csv"))

latlong_unique<-unique(temp_dataset[c("Latitude","Longitude")])

Water<-st_read(file.path(data_root,"Shapefiles for Map Figure","Hydro_poly_UTM10Copy.shp"))

crsLONGLAT <- "+proj=longlat +datum=WGS84 +no_defs"

fig1<-ggplot() + theme_bw()+
  geom_sf(data = Water, fill="slategray1", color="slategray2") +
  geom_sf(data = Delta, alpha=0.1) + 
  geom_sf(data=latlong_unique, aes(fill="black"))+
  coord_sf(xlim = c(-122.5, -120.9), ylim = c(37.4, 38.61),crs=crsLONGLAT)  +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_y = unit(1.0, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "tr", width_hint = 0.5)+
  #theme(legend.position = "none") +
  theme(axis.text.x = element_text(size=16, color="black"),axis.text.y = element_text(size=16, color="black"),axis.title.x=element_blank(),axis.title.y=element_blank())

fig1

