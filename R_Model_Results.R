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
                       n=100, 
                       Temperature_anomaly_range=c(-1.5,0,1.5),
                       Julian_days_range=yday(ymd(paste("2001", c(1,4,7,10), "15", sep="-"))) #Jan, Apr, Jul, and Oct 15 for a non-leap year
){
  
  # Create point locations on a grid for predictions
  Points<-st_make_grid(Delta, n=100)%>%
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
model_predictions<-predict(model_05_soapfilm_xy_jd_ta, newdata=newdata, type="response",discrete=F, se.fit=TRUE) # Create predictions

#Save model prediction results for later use
saveRDS(model_predictions, file.path(results_root,"Model_Predictions_for_Figures.Rds"))

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


#Create figure
plot_model_results<-ggplot(data=newdata_edit)+
  geom_sf(aes(colour=Prediction),pch=15)+
  scale_colour_gradient2(low = "blue",high = "red",midpoint = 0,breaks=seq(min(newdata_edit$Prediction),max(newdata_edit$Prediction),(max(newdata_edit$Prediction)-min(newdata_edit$Prediction))/5))+
  theme_dark()+
  facet_grid(Season~Temperature_anomaly_category)+
  theme(plot.title=element_text(size=28), 
        axis.text.x=element_text(size=21, color="black"), 
        axis.text.y = element_text(size=20, color="black"), 
        axis.title.x = element_text(size = 22, angle = 00), 
        axis.title.y = element_text(size = 22, angle = 90),
        strip.text = element_text(size = 20))+
  labs(x="Temperature Anomaly", y="Season")


png(filename=file.path(results_root,"Model_prediction_map.png"), units="in",type="cairo", bg="white", height=18, 
    width=20, res=500, pointsize=20)
plot_model_results
dev.off()


#Create dataset for surface temperature "prediction" to see what the surface temperature looks like
newdata_edit$Temperature_prediction<-predict(temperature_anomaly_GAM_spatial,newdata_edit)
newdata_edit$Temperature_prediction<-newdata_edit$Temperature_prediction+as.numeric(newdata_edit$Temperature_anomaly_category)

#Create function to create ggplot for surface temperature
temperature_plot<-function(Full_data=newdata_edit,
                           Season_set="Winter"
){
  newplot<-ggplot(data=(Full_data %>% filter(Season==Season_set)))+
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
    labs(y=Season_set,colour=NULL)
  return(newplot)
}

#Create surface temp plots for each season
plot_summer_surface<-temperature_plot(Season_set = "Summer")
plot_winter_surface<-temperature_plot(Season_set = "Winter")
plot_fall_surface<-temperature_plot(Season_set = "Fall")
plot_spring_surface<-temperature_plot(Season_set = "Spring")

#Add title
plot_spring_surface<-plot_spring_surface+labs(title="Surface temperature")


#Print out figure
png(filename=file.path(results_root,"Model_temperature_map.png"), units="in",type="cairo", bg="white", height=18, 
    width=20, res=500, pointsize=20)
ggarrange(plot_spring_surface, plot_summer_surface, plot_fall_surface, plot_winter_surface, ncol=1, nrow=4)
dev.off()

plot_surface_full<-ggarrange(plot_spring_surface, plot_summer_surface, plot_fall_surface, plot_winter_surface, ncol=1, nrow=4)


#Check limits for temperature difference predictions
summary(newdata_edit$Prediction)

#Create function to create ggplot for model results based on season
model_results_plot<-function(Full_data=newdata_edit,
                             Season_set="Winter"
){
  newplot<-ggplot(data=(Full_data %>% filter(Season==Season_set)))+
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
    labs(y=Season_set,colour=NULL)
  return(newplot)
}

plot_summer_results<-model_results_plot(Season_set = "Summer")
plot_winter_results<-model_results_plot(Season_set = "Winter")
plot_spring_results<-model_results_plot(Season_set = "Spring")
plot_fall_results<-model_results_plot(Season_set = "Fall")

plot_spring_results<-plot_spring_results+labs(title="Bottom temperature difference from surface")

plot_results_full<-ggarrange(plot_spring_results, plot_summer_results, plot_fall_results, plot_winter_results, ncol=1, nrow=4)

#Print out figure
png(filename=file.path(results_root,"Model_full_results.png"), units="in",type="cairo", bg="white", height=18, 
    width=20, res=300, pointsize=20)
ggarrange(plot_surface_full, plot_results_full, ncol=2, nrow=1)
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