###################################################
#Script to conduct model selection-------------------------------
###################################################

library(tidyverse)
library(lubridate)
library(broom)
library(rgdal)
library(mgcv)
library(AICcmodavg)
library(sf)
library(ggpubr)
library(gtable)
library(grid)
library(gridExtra)
library(patchwork)
library(png)

data_root<-file.path("data-raw")
results_root<-file.path("results")


######################################################################################################
################# Load prepared datasets ############
######################################################################################################

#Border outline for soap-film smoother
border.aut <- readRDS("Soap_film_boundaries.Rds")

temp_dataset <- readRDS("temperature_dataset.Rds")

#Add knots
knots_grid <- read.csv("knots_grid.csv")

#Read water boundaries shape file
Water<-st_read(file.path(data_root,"Shapefiles for Map Figure","Hydro_poly_UTM10Copy.shp"))
crsLONGLAT <- "+proj=longlat +datum=WGS84 +no_defs"

######################################################################################################
################# Reload model ############
######################################################################################################

#Read saved model results
model_05_soapfilm_xy_jd_ta<-readRDS(file.path(results_root,"Model_05_soapfilm.Rds"))

temperature_anomaly_GAM_spatial<-readRDS(file.path(results_root,"Temperature_anomaly_spatial_GAM.Rds"))

#Remove spatial component from dataset
temp_dataset$geometry<-NULL

#Create prediction column (not sure if necessary)
temp_dataset$Temperature_difference_prediction<-as.vector(predict(model_05_soapfilm_xy_jd_ta, newdata=temp_dataset, type="response", se.fit=F, n.threads=3))

hist(temp_dataset$Temperature_anomaly_spatial)

#Temperature anomaly data is fairly well-distributed with decent chunk of data between -2 and 2
#Add categories to temperature anomaly to ensure that we don't plot predictions in which we have no data for
temp_dataset$Temperature_anomaly_category<-0
temp_dataset$Temperature_anomaly_category<-ifelse(temp_dataset$Temperature_anomaly_spatial<=(-1.5),-1.5,temp_dataset$Temperature_anomaly_category)
temp_dataset$Temperature_anomaly_category<-ifelse(temp_dataset$Temperature_anomaly_spatial>=1.5,1.5,temp_dataset$Temperature_anomaly_category)


#Load Delta Shapefile
Delta<-st_read(file.path(results_root,"SFE_regionsFinal.shp"))
#Load WQ stations file
WQ_stations<-readRDS(file.path(results_root,"WQ_Station_boundaries.Rds"))
#Load Delta waterways from the map figure
Delta_water<-st_read(file.path(data_root,"Shapefiles for Map Figure","Hydro_poly_UTM10Copy.shp"))
#Load boundary lines from Mike Beakes
Delta_outline<-st_read(file.path(data_root,"Bay_Delta_Poly_Outline3_UTM10", "Bay_Delta_Poly_Outline_NoSSC_UTM10.shp"))

######################################################################################################
################# Figure with just four seasons ############
######################################################################################################

#Function from Sam to create prediction dataset
WQ_pred_grid<-function(Full_data=temp_dataset,
                       n=150, 
                       Temperature_anomaly_range=c(-1.5,0,1.5),
                       Julian_days_range=yday(ymd(paste("2001", c(1,4,7,10), "15", sep="-"))) #Jan, Apr, Jul, and Oct 15 for a non-leap year
){
  
  # Create point locations on a grid for predictions
  Points<-st_make_grid(Delta, n=n)%>%
    st_as_sf(crs=st_crs(Delta)) %>%
    st_filter(Delta_outline%>%
                st_transform(crs=st_crs(Delta)))%>% #Filter by the outline map from Mike Beakes
    st_filter(Delta_water%>% #Filter again for islands using the map figure water layer
                st_transform(crs=st_crs(Delta)))%>%
    st_join(WQ_stations%>% # Applying the same approach we did to the full data: remove any points outside the convex hull formed by major survey stations sampled >x times
              st_union()%>%
              st_convex_hull()%>%
              st_as_sf()%>%
              mutate(IN=TRUE),
            join=st_intersects)%>%
    filter(IN)%>%
    dplyr::select(-IN)%>%
    st_centroid()%>% # The prior grid was actually a set of polygons, this picks the center point of each
    st_coordinates() %>%
    as_tibble()%>%
    mutate(Location=1:nrow(.))
  
  # Create dataset for each year and season showing which regions were sampled and have the right temperature anomaly range
  Data_effort <- Full_data%>%
    group_by(SubRegion, Season,Temperature_anomaly_category)%>%
    summarise(N=n())%>%
    ungroup()%>%
    left_join(Delta, by="SubRegion")%>%
    dplyr::select(-geometry)
  
  # Create full dataset for predictions
  newdata<-expand.grid(Temperature_anomaly_spatial= Temperature_anomaly_range,
                       Location=1:nrow(Points),
                       Julian_day=Julian_days_range)%>% # Create all combinations of predictor variables
    left_join(Points, by="Location")%>% #Add Lat/Longs to each location
    mutate(Julian_day_s = (Julian_day-mean(Full_data$Julian_day, na.rm=T))/sd(Full_data$Julian_day, na.rm=T), # Standardize each variable based on full dataset for model
           Temperature_anomaly_category=Temperature_anomaly_spatial,
           Season=case_when(Julian_day<=80 | Julian_day>=356 ~ "Winter", # Create a variable for season
                            Julian_day>80 & Julian_day<=172 ~ "Spring",
                            Julian_day>=173 & Julian_day<=264 ~ "Summer",
                            Julian_day>=265 & Julian_day<=355 ~ "Fall"))%>%
    st_as_sf(coords=c("X", "Y"),crs=st_crs(Delta), remove=FALSE)%>%
    st_join(Delta, join = st_intersects)%>%
    filter(!is.na(SubRegion))%>% # Make sure all points are within our desired subregions
    left_join(Data_effort, by=c("SubRegion", "Season","Temperature_anomaly_category"))%>% # Use the Data_effort key created above to remove points in subregions that were not sampled that region, season, and year.
    filter(!is.na(N))%>%
    rename(x=X,y=Y)
  return(newdata)
}

#Use the function to create prediction data frame
newdata <- WQ_pred_grid()

#Add prediction from model
#model_predictions<-predict(model_05_soapfilm_xy_jd_ta, newdata=newdata, type="response",discrete=F, se.fit=TRUE) # Create predictions

#Save model prediction results for later use
#saveRDS(model_predictions, file.path(results_root,"Model_Predictions_for_Figures.Rds"))

#Load model predictions for faster script run
model_predictions<-readRDS(file.path(results_root,"Model_Predictions_for_Figures.Rds"))

#Reconfigure the data set
newdata<-newdata%>%
  mutate(Prediction=model_predictions$fit)%>%
  mutate(SE=model_predictions$se.fit,
         L95=Prediction-SE*1.96,
         U95=Prediction+SE*1.96)

#There are a bunch of NAs from model predictions, presumably new data being outside the boundaries?
newdata$Prediction

#Remove NAs
newdata_edit<-newdata[!is.na(newdata$Prediction),]
#Change temp anomaly category to factor
newdata_edit$Temperature_anomaly_category<-as.factor(newdata_edit$Temperature_anomaly_category)
#Change 1.5 C temp to +1.5 C for better clarity
levels(newdata_edit$Temperature_anomaly_category)[levels(newdata_edit$Temperature_anomaly_category)=="1.5"] <- "+1.5"


###########################################################################
#Create dataset for surface temperature "prediction" to see what the surface temperature looks like

newdata_edit$Temperature_prediction<-predict(temperature_anomaly_GAM_spatial,newdata_edit)
newdata_edit$Temperature_prediction<-newdata_edit$Temperature_prediction+as.numeric(newdata_edit$Temperature_anomaly_category)

#Create function to create ggplot for surface temperature
temperature_plot<-function(Full_data=newdata_edit,
                           Season_set="Winter",Season_name="Winter (Jan 15)"
){
  newplot<-ggplot(data=(Full_data %>% filter(Season==Season_set)))+
    geom_sf(aes(colour=Temperature_prediction),pch=15)+
    scale_colour_gradient(low = "yellow",high = "red",name=expression(""~degree * C *""))+
    theme_dark()+
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(size = 20, angle = 00), 
          axis.title.y = element_text(size = 18, angle = 90),
          strip.text = element_text(size = 16),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          legend.key.size = unit(1.2, 'cm'),
          plot.title = element_text(size=18)
    )+
    guides(colour=guide_colourbar(ticks.colour = "black"))+
    facet_grid(~Temperature_anomaly_category)+
    labs(y=Season_name,colour=NULL)
  return(newplot)
}

#Create surface temp plots for each season
plot_summer_surface<-temperature_plot(Season_set = "Summer",Season_name="Summer (Jul 15)")
plot_winter_surface<-temperature_plot(Season_set = "Winter",Season_name="Winter (Jan 15)")
plot_fall_surface<-temperature_plot(Season_set = "Fall",Season_name="Fall (Oct 15)")
plot_spring_surface<-temperature_plot(Season_set = "Spring",Season_name="Spring (Apr 15)")

#Add titles
#plot_winter_surface<-plot_winter_surface+labs(title="Temperature anomaly category")+theme(plot.title = element_text(hjust = 0.5))

plot_surface_full<-ggarrange(plot_winter_surface, plot_spring_surface, plot_summer_surface, plot_fall_surface, ncol=1, nrow=4, align = "hv",legend="right")

#Check limits for temperature difference predictions
summary(newdata_edit$Prediction)

###########################################################################
#Create function to create ggplot for model results based on season

model_results_plot<-function(Full_data=newdata_edit,
                             Season_set="Winter"
){
  newplot<-ggplot(data=(Full_data %>% filter(Season==Season_set)))+
    geom_sf(aes(colour=Prediction),pch=15)+
    #scale_colour_gradient2(limits=c(-3.5,0.5),low = "blue",high = "red",midpoint = 0)+
    scale_colour_viridis_c(name=expression(""~degree * C *""))+
    theme_dark()+
    facet_grid(~Temperature_anomaly_category)+
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(size = 20, angle = 00), 
          axis.title.y = element_text(size = 18, angle = 90),
          strip.text = element_text(size = 16),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          legend.key.size = unit(1.2, 'cm'),
          plot.title = element_text(size=18)
    )+
    guides(colour=guide_colourbar(ticks.colour = "black"))+
    labs(y=NULL,colour=NULL)
  return(newplot)
}

plot_summer_results<-model_results_plot(Season_set = "Summer")
plot_winter_results<-model_results_plot(Season_set = "Winter")
plot_spring_results<-model_results_plot(Season_set = "Spring")
plot_fall_results<-model_results_plot(Season_set = "Fall")

#plot_winter_results<-plot_winter_results+labs(title="Temperature anomaly category")+theme(plot.title = element_text(hjust = 0.5))

plot_results_full<-ggarrange(plot_winter_results, plot_spring_results, plot_summer_results, plot_fall_results, ncol=1, nrow=4, align = "hv",legend="right")

#Print out figure
tiff(filename=file.path(results_root,"Model_full_results.tiff"), units="in",type="cairo", bg="white", height=15, 
    width=20, res=300, pointsize=10,compression="lzw")
ggarrange(plot_surface_full, plot_results_full, ncol=2, nrow=1)
dev.off()




######################################################################################################
################# Suitability of Delta Smelt ############
######################################################################################################

#Extract july dataset from previous code/figure
deltasmelt_data<-newdata_edit %>% filter(Season == "Summer",Temperature_anomaly_spatial==0) %>% 
  mutate(Suitability_DSM=case_when(Temperature_prediction<20 ~ "<20 at surface", # Create a variable for smelt suitability index
                                   Temperature_prediction>=20 & Temperature_prediction<=25 ~ "20-25 at surface",
                                   Temperature_prediction>25 & Temperature_prediction+Prediction>25~ ">25 at surface and bottom",
                                   Temperature_prediction>25 & Temperature_prediction+Prediction<=25~ ">25 at surface but <=25 at bottom"))

#Create  different categories based on Brown et al. 25 C cutoff
deltasmelt_data$Suitability_DSM<-factor(deltasmelt_data$Suitability_DSM, levels = c( "20-25 at surface", ">25 at surface but <=25 at bottom",">25 at surface and bottom","<20 at surface"))

#Add custom color set
color_smelt<- c('#ffe119','#dcbeff','#e6194B','#f58231')

#Photo of Delta Smelt
pic_deltasmelt <- readPNG(file.path(data_root, "DSM_edit.png"), native = TRUE)

#Create Delta Smelt temperature suitability map
plot_deltasmelt <-ggplot(data=deltasmelt_data)+
  geom_sf(aes(colour=Suitability_DSM),pch=15)+
  geom_sf(data = Water, color=alpha("black",0.3),fill=NA) +
  #geom_raster(aes(x=x,y=y,fill=Suitability_DSM))+
  #geom_sf(data = deltasmelt_fill, color="red", fill="black") +
  scale_color_manual(values=color_smelt,labels = c(expression("20-25"~degree * C *" at surface"), expression(">25"~degree * C *" at surface, but"<="25"~degree * C * " at bottom"),expression(">25"~degree * C *" at surface and bottom")))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme_void()+
  #labs(title="Delta Smelt Temperature Suitability")+
  coord_sf(xlim = c(-122.2, -121.37), ylim = c(37.8, 38.61),crs=crsLONGLAT)  +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 7, angle = 00), 
        axis.title.y = element_text(size = 7, angle = 90),
        strip.text = element_text(size = 7),
        legend.text = element_text(size=9),
        #legend.key.size = unit(0.5, "lines"),
        #legend.key.width = unit(1.2,"lines"),
        legend.position = c(0.2, 0.1),
        legend.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"), 
        legend.background = element_rect(fill = "white",color="black"),
        legend.title = element_blank(),
        plot.title = element_text(size=8, hjust=0.5)
  )

#Add Delta Smelt picture to plot
plot_deltasmelt <- plot_deltasmelt+
  inset_element(p = pic_deltasmelt,
                left = 0.1,
                bottom = 0.7,
                right = 0.5,
                top = 0.9)

#Print out figure
#tiff(filename=file.path(results_root,"Figure_DeltaSmelt_Suitability.tiff"), units="in",type="cairo", bg="white", height=6, 
#     width=8, res=300, pointsize=10,compression="lzw")
#plot_deltasmelt
#dev.off()


######################################################################################################
################# Suitability of juvenile Chinook Salmon ############
######################################################################################################

#See Myrick and Cech 2004

chinooksalmon_data <- WQ_pred_grid(Julian_days_range=yday(ymd(paste("2001", "05", "15", sep="-")))) 
chinooksalmon_predictions <- predict(model_05_soapfilm_xy_jd_ta, newdata=chinooksalmon_data, type="response",discrete=F, se.fit=TRUE)

chinooksalmon_data$Temperature_prediction=predict(temperature_anomaly_GAM_spatial,chinooksalmon_data)

chinooksalmon_data<-chinooksalmon_data %>%
  mutate(Prediction=chinooksalmon_predictions$fit)%>%
  mutate(SE=chinooksalmon_predictions$se.fit,
         L95=Prediction-SE*1.96,
         U95=Prediction+SE*1.96) %>% 
  filter(Temperature_anomaly_spatial==0) %>% 
  mutate(Suitability_CHN=case_when(Temperature_prediction>=13 & Temperature_prediction<=16 & Temperature_prediction+Prediction<=16 ~ "13-16 at surface and bottom", # Create a variable for smelt suitability index
                                   Temperature_prediction>=13 & Temperature_prediction<=16 & Temperature_prediction+Prediction>16 ~ "13-16 at surface and >16 bottom", 
                                   Temperature_prediction>16 & Temperature_prediction+Prediction<=16 ~ ">16 at surface but <=16 at bottom",
                                   Temperature_prediction>16 & Temperature_prediction+Prediction>16 ~ ">16 at surface and bottom"))

summary(as.factor(chinooksalmon_data$Suitability_CHN))

#Remove NAs
chinooksalmon_data<-chinooksalmon_data[!is.na(chinooksalmon_data$Prediction),]

chinooksalmon_data$Suitability_CHN<-factor(chinooksalmon_data$Suitability_CHN, levels = c( "13-16 at surface and bottom","13-16 at surface and >16 bottom", ">16 at surface but <=16 at bottom",">16 at surface and bottom"))


color_salmon <- c('#4363d8','#42d4f4','#f032e6','#e6194B')


#Photo of Chinook Salmon
pic_chinooksalmon <- readPNG(file.path(data_root, "CHN_edit.png"), native = TRUE)


plot_chinooksalmon <-ggplot(data=chinooksalmon_data)+
  geom_sf(aes(colour=Suitability_CHN),pch=15)+
  geom_sf(data = Water, color=alpha("black",0.3),fill=NA) +
  scale_color_manual(values=color_salmon,labels = c(expression(""<="16"~degree * C *" at surface and bottom"), expression(""<="16"~degree * C *" at surface, but">"16"~degree * C * " at bottom"), expression(">16"~degree * C *" at surface, but"<="16"~degree * C * " at bottom"),expression(">16"~degree * C *" at surface and bottom")))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme_void()+
  #labs(title="Chinook Salmon Temperature Suitability")+
  coord_sf(xlim = c(-122.2, -121.37), ylim = c(37.8, 38.61),crs=crsLONGLAT)  +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 7, angle = 00), 
        axis.title.y = element_text(size = 7, angle = 90),
        strip.text = element_text(size = 7),
        legend.text = element_text(size=9),
        #legend.key.size = unit(0.5, "lines"),
        #legend.key.width = unit(1.2,"lines"),
        legend.position = c(0.2, 0.1),
        legend.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"), 
        legend.background = element_rect(fill = "white",color="black"),
        legend.title = element_blank(),
        plot.title = element_text(size=8, hjust=0.5)
  )

#Add Chinook Salmon picture
plot_chinooksalmon<- plot_chinooksalmon+
  inset_element(p = pic_chinooksalmon,
                left = 0.1,
                bottom = 0.7,
                right = 0.5,
                top = 0.9)


#Print out figure
#tiff(filename=file.path(results_root,"Figure_ChinookSalmon_Suitability.tiff"), units="in",type="cairo", bg="white", height=6, 
#     width=8, res=300, pointsize=10,compression="lzw")
#plot_chinooksalmon
#dev.off()

#Print out smelt and salmon figure side by side
tiff(filename=file.path(results_root,"Figure_Temperature_Suitability.tiff"), units="in",type="cairo", bg="white", height=8, 
     width=14, res=300, pointsize=10,compression="lzw")
ggarrange(plot_deltasmelt, plot_chinooksalmon, ncol=2, nrow=1,labels = c("A - Delta Smelt (July 15th)", "B - Chinook Salmon (May 15th)"),hjust=-0.1)
dev.off()

######################################################################################################
################# Figure for every month (12x3) ############
######################################################################################################

#Use the function to create prediction data frame
newdata_12_months <- WQ_pred_grid(Julian_days_range=yday(ymd(paste("2001", c(1:12), "15", sep="-"))))

#Add prediction from model
model_predictions_12_months<-predict(model_05_soapfilm_xy_jd_ta, newdata=newdata_12_months, type="response",discrete=F, se.fit=TRUE) # Create predictions

#Reconfigure the data set
newdata_12_months<-newdata_12_months%>%
  mutate(Prediction=model_predictions_12_months$fit)%>%
  mutate(SE=model_predictions_12_months$se.fit,
         L95=Prediction-SE*1.96,
         U95=Prediction+SE*1.96)

#Remove NAs
newdata_12_months<-newdata_12_months[!is.na(newdata_12_months$Prediction),]

#Change temp anomaly category to factor
newdata_12_months$Temperature_anomaly_category<-as.factor(newdata_12_months$Temperature_anomaly_category)
#Change 1.5 C temp to +1.5 C for better clarity
levels(newdata_12_months$Temperature_anomaly_category)[levels(newdata_12_months$Temperature_anomaly_category)=="1.5"] <- "+1.5"


#Create dataset for surface temperature "prediction" to see what the surface temperature looks like
newdata_12_months$Temperature_prediction<-predict(temperature_anomaly_GAM_spatial,newdata_12_months)
newdata_12_months$Temperature_prediction<-newdata_12_months$Temperature_prediction+as.numeric(newdata_12_months$Temperature_anomaly_category)

#Add month factor
newdata_12_months$Month<-as.factor(month(as.Date(newdata_12_months$Julian_day-1, origin = "2001-01-01"), label = TRUE, abbr = FALSE))
unique(newdata_12_months$Month)

#Create function to create ggplot for surface temperature
temperature_plot<-function(Full_data=newdata_12_months,
                           Month_set="January"
){
  newplot<-ggplot(data=(Full_data %>% filter(Month==Month_set)))+
    geom_sf(aes(colour=Temperature_prediction),pch=15)+
    scale_colour_gradient(low = "yellow",high = "red")+
    theme_dark()+
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(size = 22, angle = 00), 
          axis.title.y = element_text(size = 22, angle = 90),
          strip.text = element_text(size = 20),
          legend.text = element_text(size=18),
          legend.key.size = unit(1.5, 'cm'),
          plot.title = element_text(size=22)
    )+
    guides(colour=guide_colourbar(ticks.colour = "black"))+
    facet_grid(~Temperature_anomaly_category)+
    labs(y=Month_set,colour=NULL)
  return(newplot)
}

#Create surface temp plots for each month
plot_surface_01<-temperature_plot(Month_set = "January")
plot_surface_02<-temperature_plot(Month_set = "February")
plot_surface_03<-temperature_plot(Month_set = "March")
plot_surface_04<-temperature_plot(Month_set = "April")
plot_surface_05<-temperature_plot(Month_set = "May")
plot_surface_06<-temperature_plot(Month_set = "June")
plot_surface_07<-temperature_plot(Month_set = "July")
plot_surface_08<-temperature_plot(Month_set = "August")
plot_surface_09<-temperature_plot(Month_set = "September")
plot_surface_10<-temperature_plot(Month_set = "October")
plot_surface_11<-temperature_plot(Month_set = "November")
plot_surface_12<-temperature_plot(Month_set = "December")

#Add title
plot_surface_01<-plot_surface_01+labs(title="Surface temperature")

#Create function to create ggplot for model results based on month
model_results_plot<-function(Full_data=newdata_12_months,
                             Month_set="January"
){
  newplot<-ggplot(data=(Full_data %>% filter(Month==Month_set)))+
    geom_sf(aes(colour=Prediction),pch=15)+
    scale_colour_gradient2(limits=c(-3.5,0.5),low = "blue",high = "red",midpoint = 0)+
    theme_dark()+
    facet_grid(~Temperature_anomaly_category)+
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(size = 22, angle = 00), 
          axis.title.y = element_text(size = 22, angle = 90),
          strip.text = element_text(size = 20),
          legend.text = element_text(size=18),
          legend.key.size = unit(1.5, 'cm'),
          plot.title = element_text(size=22)
    )+
    labs(y=Month_set,colour=NULL)
  return(newplot)
}

#Create prediction plots for each month
plot_results_01<-model_results_plot(Month_set = "January")
plot_results_02<-model_results_plot(Month_set = "February")
plot_results_03<-model_results_plot(Month_set = "March")
plot_results_04<-model_results_plot(Month_set = "April")
plot_results_05<-model_results_plot(Month_set = "May")
plot_results_06<-model_results_plot(Month_set = "June")
plot_results_07<-model_results_plot(Month_set = "July")
plot_results_08<-model_results_plot(Month_set = "August")
plot_results_09<-model_results_plot(Month_set = "September")
plot_results_10<-model_results_plot(Month_set = "October")
plot_results_11<-model_results_plot(Month_set = "November")
plot_results_12<-model_results_plot(Month_set = "December")

#Add title
plot_results_01<-plot_results_01+labs(title="Bottom temperature difference from surface")

#Print out Jan-April figure
png(filename=file.path(results_root,"Model_full_results_01-04.png"), units="in",type="cairo", bg="white", height=18, 
    width=20, res=300, pointsize=20)
ggarrange(ggarrange(plot_surface_01, plot_surface_02, plot_surface_03, plot_surface_04, ncol=1, nrow=4), ggarrange(plot_results_01, plot_results_02, plot_results_03, plot_results_04, ncol=1, nrow=4), ncol=2, nrow=1)
dev.off()

#Print out May-August figure
png(filename=file.path(results_root,"Model_full_results_05-08.png"), units="in",type="cairo", bg="white", height=18, 
    width=20, res=300, pointsize=20)
ggarrange(ggarrange(plot_surface_05, plot_surface_06, plot_surface_07, plot_surface_08, ncol=1, nrow=4), ggarrange(plot_results_05, plot_results_06, plot_results_07, plot_results_08, ncol=1, nrow=4), ncol=2, nrow=1)
dev.off()

#Print out Sep-December figure
png(filename=file.path(results_root,"Model_full_results_09-12.png"), units="in",type="cairo", bg="white", height=18, 
    width=20, res=300, pointsize=20)
ggarrange(ggarrange(plot_surface_09, plot_surface_10, plot_surface_11, plot_surface_12, ncol=1, nrow=4), ggarrange(plot_results_09, plot_results_10, plot_results_11, plot_results_12, ncol=1, nrow=4), ncol=2, nrow=1)
dev.off()
