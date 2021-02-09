library(tidyverse)
library(lubridate)
library(broom)
library(rgdal)
library(mgcv)
library(AICcmodavg)
library(sf)

source("soap_checker/soap_check.R")
data_root<-file.path("data-raw")
results_root<-file.path("results")

#Read in bay-Delta shape outline shape file that Mike Beakes created
Delta.aut <- readOGR(file.path(data_root,"Bay_Delta_Poly_Outline3_UTM10", "Bay_Delta_Poly_Outline_NoSSC_UTM10.shp"))

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

#Load knots
knots_grid <- read.csv("knots_grid.csv")


######################################################################################################
################# Read dataset and create temperature anomaly calculation ############
######################################################################################################

temp_dataset<-read.csv("temperature_dataset.csv")
str(temp_dataset)

# Calculate difference of bottom temperature from surface temperature as response variable
temp_dataset$Temperature_difference <- temp_dataset$Temperature_bottom-temp_dataset$Temperature
hist(temp_dataset$Temperature_difference)

#Temperature Anomaly Model---------------------------------------------
##Julian day and temperature model
temperature_anomaly_GAM<- gam(Temperature ~ s(Julian_day_s,bs="cc",k=5),data=temp_dataset)
summary(temperature_anomaly_GAM)
gam.check(temperature_anomaly_GAM)
plot(temperature_anomaly_GAM)

#Add term to the dataset
temp_dataset$Temperature_prediction <-predict(temperature_anomaly_GAM,temp_dataset)
temp_dataset$Temperature_anomaly <-temp_dataset$Temperature - temp_dataset$Temperature_prediction 
hist(temp_dataset$Temperature_anomaly)

############Summary of data
summary(temp_dataset$Temperature_difference)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-5.7600 -0.3000 -0.1000 -0.1669  0.0000  6.5600 
sd(temp_dataset$Temperature_difference)
#0.4546607

#Summarize by region for reference
temp_dataset$SubRegion=as.factor(temp_dataset$SubRegion)
temp_dataset$SampleSize=1

data_by_subregion<-temp_dataset %>% group_by(SubRegion) %>% summarize(SampleSize=sum(SampleSize),MeanTempDiff=mean(Temperature_difference),Temperature_anomaly_mean=mean(Temperature_anomaly),Temperature_anomaly_max=max(Temperature_anomaly),Temperature_anomaly_min=min(Temperature_anomaly))

data_by_subregion_month<-temp_dataset %>% group_by(SubRegion,Month) %>% summarise(SampleSize=sum(SampleSize),MeanTempDiff=mean(Temperature_difference),Temperature_anomaly_mean=mean(Temperature_anomaly),Temperature_anomaly_max=max(Temperature_anomaly),Temperature_anomaly_min=min(Temperature_anomaly))

data_by_subregion_year<-temp_dataset %>% group_by(SubRegion,Year) %>% summarise(SampleSize=sum(SampleSize),MeanTempDiff=mean(Temperature_difference),Temperature_anomaly_mean=mean(Temperature_anomaly),Temperature_anomaly_max=max(Temperature_anomaly),Temperature_anomaly_min=min(Temperature_anomaly))

data_by_year<-temp_dataset %>% group_by(Year) %>% summarise(SampleSize=sum(SampleSize),MeanTempDiff=mean(Temperature_difference),Temperature_anomaly_mean=mean(Temperature_anomaly),Temperature_anomaly_max=max(Temperature_anomaly),Temperature_anomaly_min=min(Temperature_anomaly))

data_by_month<-temp_dataset %>% group_by(Month) %>% summarise(SampleSize=sum(SampleSize),MeanTempDiff=mean(Temperature_difference),Temperature_anomaly_mean=mean(Temperature_anomaly),Temperature_anomaly_max=max(Temperature_anomaly),Temperature_anomaly_min=min(Temperature_anomaly))

#By inspecting the summarized data by subregion, it's apparent that temperature anomaly varies by region.
#South Delta is always warmer than North Delta as an example


##Construct model with longitude and latitude on top of julian day
temperature_anomaly_GAM_spatial<- gam(Temperature ~ te(x,y,Julian_day_s, d=c(2,1) ,bs=c("tp","cc"),k=c(10,5)),data=temp_dataset)
summary(temperature_anomaly_GAM_spatial)
#R-sq.(adj) =  0.907   Deviance explained = 90.7%

#K index is too low at 0.5, and edf is pretty close to k', but the goal was to remove collinearity, not fit
#k'  edf k-index p-value    
#39.0 35.8     0.5  <2e-16 ***
gam.check(temperature_anomaly_GAM_spatial)
plot(temperature_anomaly_GAM_spatial)

#Add term to the dataset
temp_dataset$Temperature_prediction_spatial <-predict(temperature_anomaly_GAM_spatial,temp_dataset)
temp_dataset$Temperature_anomaly_spatial <-temp_dataset$Temperature - temp_dataset$Temperature_prediction_spatial 
hist(temp_dataset$Temperature_anomaly_spatial)

#Make water year factor
temp_dataset$WaterYear<-as.factor(temp_dataset$WaterYear)

######################################################################################################
################# Model run with soap film smooth ############
######################################################################################################

model_05_soapfilm_xy_jd_ta <- bam(Temperature_difference ~  te(x, y,Julian_day_s,Temperature_anomaly_spatial, d=c(2,1,1), bs=c("sf", "cc","tp"), k=c(35,5,7),xt = list(list(bnd = border.aut,nmax=1000),NULL,NULL))+
                                    te(x, y, Julian_day_s,Temperature_anomaly_spatial, d=c(2,1,1), bs=c("sw","cc","tp"), k=c(35,5,7),xt = list(list(bnd = border.aut,nmax=1000),NULL,NULL)),
                                  data = temp_dataset, method="fREML", nthreads=3, knots =knots_grid)



#Save model results
saveRDS(model_05_soapfilm_xy_jd_ta, file.path(results_root,"Best_Model.Rds"))

######################################################################################################
################# Reload model ############
######################################################################################################

#Read saved model results
model_05_soapfilm_xy_jd_ta<-readRDS(file.path(results_root,"Best_Model.Rds"))

#Create prediction column (not sure if necessary)
temp_dataset$Temperature_difference_prediction<-as.vector(predict(model_05_soapfilm_xy_jd_ta, newdata=temp_dataset, type="response", se.fit=F, n.threads=3))

test2<-predict(model_05_soapfilm_xy_jd_ta, newdata=temp_dataset, type="response",discrete=F, se.fit=T, n.threads=3)

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


#Function from Sam to create prediction dataset
WQ_pred<-function(Full_data=temp_dataset,
                  n=100, 
                  Temperature_anomaly_range=c(-1.5,0,1.5),
                  Julian_days_range=yday(ymd(paste("2001", c(1,4,7,10), "15", sep="-"))) #Jan, Apr, Jul, and Oct 15 for a non-leap year
){
  
  # Create point locations on a grid for predictions
  Points<-st_make_grid(Delta, n=n)%>%
    st_as_sf(crs=st_crs(Delta)) %>%
    st_filter(Delta_water%>%
                st_transform(crs=st_crs(Delta)))%>% #Filter through the map figure file
    st_filter(Delta_outline%>%
                st_transform(crs=st_crs(Delta)))%>% #Filter through the outline shape
    st_join(WQ_stations%>% # Applying the same approach we did to the full data: remove any points outside the convex hull formed by major survey stations sampled >50 times
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
  
  # Create dataset for each year and season showing which subregions were sampled
  Data_effort <- Full_data%>%
    group_by(SubRegion, Season,Temperature_anomaly_category)%>%
    summarise(N=n())%>%
    ungroup()%>%
    left_join(Delta_subregions, by="SubRegion")%>%
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
    st_as_sf(coords=c("X", "Y"), crs="NAD83", remove=FALSE)%>% # Turn into sf object
    #st_transform(crs=st_crs(Delta))%>% # transform to crs of Delta shapefile
    st_join(Delta, join = st_intersects)%>%
    filter(!is.na(SubRegion))%>% # Make sure all points are within our desired subregions
    left_join(Data_effort, by=c("SubRegion", "Season","Temperature_anomaly_category"))%>% # Use the Data_effort key created above to remove points in subregions that were not sampled that region, season, and year.
    filter(!is.na(N))
  return(newdata)
}

#Use the function to create prediction data frame

newdata_bottom <- WQ_pred(Full_data=temp_dataset, 
                          Julian_days_range = yday(ymd(paste("2001", 1:12, "15", sep="-"))))





#############CODE TESTING BELOW

Points<-st_make_grid(Delta, n=100)%>%
  st_as_sf(crs=st_crs(Delta)) %>%
  #st_join(spacetools::Delta %>% # Joining a map of delta waterways (from my spacetools package) to ensure all these points are over water.
  #         dplyr::select(Shape_Area)%>%
  #         st_transform(crs=st_crs(Delta)))%>%
  st_filter(Delta_outline%>%
              st_transform(crs=st_crs(Delta)))%>%
  st_filter(Delta_water%>%
              st_transform(crs=st_crs(Delta)))%>%
  st_join(WQ_stations%>% # Applying the same approach we did to the full data: remove any points outside the convex hull formed by major survey stations sampled >50 times
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

Points

Data_effort <- temp_dataset%>%
  #st_drop_geometry()%>%
  group_by(SubRegion, Season,Temperature_anomaly_category)%>%
  summarise(N=n())%>%
  ungroup()%>%
  left_join(Delta, by="SubRegion")%>%
  dplyr::select(-geometry)

newdata<-expand.grid(Temperature_anomaly_spatial= c(-1.5,0,1.5),
                     Location=1:nrow(Points),
                     Julian_day=yday(ymd(paste("2001", c(1,4,7,10), "15", sep="-"))))%>% # Create all combinations of predictor variables
  left_join(Points, by="Location")%>% #Add Lat/Longs to each location
  mutate(Julian_day_s = (Julian_day-mean(temp_dataset$Julian_day, na.rm=T))/sd(temp_dataset$Julian_day, na.rm=T), # Standardize each variable based on full dataset for model
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


#Add prediction from model
model_predictions<-predict(model_05_soapfilm_xy_jd_ta, newdata=newdata, type="response",discrete=F, se.fit=TRUE, n.threads=3) # Create predictions

#Save model prediction results
saveRDS(model_predictions, file.path(results_root,"Model_Predictions_for_Figures.Rds"))



newdata<-newdata%>%
  mutate(Prediction=model_predictions$fit)%>%
  mutate(SE=model_predictions$se.fit,
         L95=Prediction-SE*1.96,
         U95=Prediction+SE*1.96)


#Create figures for each season, time of day, and residual temperature
newdata$Prediction
str(newdata)

newdata_edit<-newdata
newdata_edit<-newdata_edit[!is.na(newdata_edit$Prediction),]
newdata_edit$Prediction
str(newdata_edit)
newdata_edit$Temperature_anomaly_category<-as.factor(newdata_edit$Temperature_anomaly_category)


png(filename=file.path(results_root,"Model_prediction_map.png"), units="in",type="cairo", bg="white", height=18, 
    width=20, res=500, pointsize=20)
ggplot(data=newdata_edit)+
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
dev.off()

min(newdata$Prediction)

newdata_edit$Temperature_prediction<-predict(temperature_anomaly_GAM_spatial,newdata_edit)
newdata_edit$Temperature_prediction<-newdata_edit$Temperature_prediction+as.numeric(newdata_edit$Temperature_anomaly_category)

png(filename=file.path(results_root,"Model_temperature_map.png"), units="in",type="cairo", bg="white", height=18, 
    width=20, res=500, pointsize=20)
ggplot(data=newdata_edit)+
  geom_sf(aes(colour=Temperature_prediction),pch=15)+
  scale_colour_gradient(low = "yellow",high = "red",breaks=seq(min(newdata_edit$Temperature_prediction),max(newdata_edit$Temperature_prediction),(max(newdata_edit$Temperature_prediction)-min(newdata_edit$Temperature_prediction))/5))+
  theme_dark()+
  facet_grid(Season~Temperature_anomaly_category)+
  theme(plot.title=element_text(size=28), 
        axis.text.x=element_text(size=21, color="black"), 
        axis.text.y = element_text(size=20, color="black"), 
        axis.title.x = element_text(size = 22, angle = 00), 
        axis.title.y = element_text(size = 22, angle = 90),
        strip.text = element_text(size = 20))+
  labs(x="Temperature Anomaly", y="Season")
dev.off()
