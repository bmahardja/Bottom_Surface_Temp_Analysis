library(tidyverse)
library(deltareportr)
library(mgcv)
library(lubridate)
library(hms)
library(sf)
library(stars)
library(AICcmodavg)
library(scales)
library(cvTools)
require(patchwork)
require(geofacet)
require(gamm4)
require(dtplyr)

#Initial step to load integrated temperature dataset and region shape files-------------------------------
#Same steps from Sam's Temperature QAQC R script

# Load Delta Shapefile from Brian
Delta<-st_read("Delta Subregions")%>%
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

# Old code to visualize the prior steps

ggplot()+
  geom_sf(data=Delta, aes(fill=SubRegion))+
  #geom_sf_label(data=st_centroid(Delta)%>%st_transform(crs=4326), aes(label=SubRegion))+
  geom_sf(data=WQ_stations%>%st_union()%>%st_convex_hull(), alpha=0.1, color="red", size=2)+
  geom_sf(data=WQ_stations)
#Give all datasets the same ending year
#max_date <- Data%>%group_by(Source)%>%summarise(Date=max(Date))%>%pull(Date)%>%min()
#Data <- filter(Data, Year<=year(max_date) & !(Source=="SKT" & Field_coords))

#Subset data to 2011 and on-------------------------------------------------------
#As they have the best spatiotemporal coverage

#Create new data frame as to not affect Sam's original dataset
#Filter just those that have bottom temperature measurements
Data_subset<- Data %>% filter(!is.na(Temperature_bottom))
summary(Data_subset$Temperature_bottom)
hist(Data_subset$Temperature_bottom)

#Remove temperatures that seem unreasonable for bottom and surface (those below 5 C and above 30 C)
Data_subset<- Data_subset %>% filter(Temperature_bottom>=5&Temperature_bottom<=30)
Data_subset<- Data_subset %>% filter(Temperature>=5&Temperature<=30)
hist(Data_subset$Temperature_bottom)

#Subset just those with Depth data
Data_subset<- Data_subset %>% filter(!is.na(Depth))

#Summarize effort by Year and month
effort_summary<- Data_subset %>% mutate(Count=1) %>% group_by(Year_fac,Month) %>% summarise(Count=sum(Count))
#Since 2011, the sample size increased considerably each month
#It looks like there are only 10-20 data points per month for all years prior to 2011, with limited spatial extent
#Only stations from Carquinez to lower Sacramento and San Joaquin Rivers
example_dataset<-Data %>% filter(year(Date)<2011,!is.na(Temperature_bottom))
ggplot()+
  geom_sf(data=Delta, aes(fill=SubRegion))+
  #geom_sf_label(data=st_centroid(Delta)%>%st_transform(crs=4326), aes(label=SubRegion))+
  geom_sf(data=example_dataset%>%st_union()%>%st_convex_hull(), alpha=0.1, color="red", size=2)+
  geom_sf(data=example_dataset)
remove(example_dataset)


#Filter just 2011-2019 for faster model runs and less spatial autocorrelation (bias)
Data_subset <- Data_subset %>% filter(Year>2010)

#Visualize stations
png(filename="Map_bottom_stations_subset.png", units="in",type="cairo", bg="white", height=12, 
    width=16, res=300, pointsize=20)
ggplot()+
  geom_sf(data=Delta, aes(fill=SubRegion))+
  #geom_sf_label(data=st_centroid(Delta)%>%st_transform(crs=4326), aes(label=SubRegion))+
  geom_sf(data=Data_subset%>%st_union()%>%st_convex_hull(), alpha=0.1, color="red", size=2)+
  geom_sf(data=Data_subset)
dev.off()

#Data evaluation--------------------------------- 

#Look at relationship between bottom and surface temperature using OLS linear regression
temp_ols <- lm(Temperature_bottom~Temperature, data=Data_subset)

summary(temp_ols)

#R^2 = 0.99
#But coefficient is slightly lower than expected = 0.98



#Initial Model Runs--------------------------------
#Models with just latitude and longitude, and surface temperature OR julian day

# Calculate difference of bottom temperature from surface temperature as response variable
Data_subset$Temperature_difference <- Data_subset$Temperature_bottom-Data_subset$Temperature
hist(Data_subset$Temperature_difference)
#Standardize temperature covariate
Data_subset$Temperature_s<-(Data_subset$Temperature-mean(Data_subset$Temperature))/sd(Data_subset$Temperature)

#Start models

model_bottom_00 <- bam(Temperature_difference ~ te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("tp", "cc"), k=c(40, 25, 7)) + s(Time_num_s, k=5, bs="cc"),
                      data = Data_subset, method="fREML", discrete=T, nthreads=4)
gam.check(model_bottom_00)
summary(model_bottom_00)
AICc(model_bottom_00)
#AICc=6376.749
#R^2=0.106
#BIC=6805.083

model_bottom_01 <- bam(Temperature_difference ~ te(Longitude_s, Latitude_s, Temperature_s, d=c(2,1), bs=c("tp"), k=c(40, 25, 7)) + s(Time_num_s, k=5, bs="cc"),
                       data = Data_subset, method="fREML", discrete=T, nthreads=4)
gam.check(model_bottom_01)
summary(model_bottom_01)
AICc(model_bottom_01)
#AICc=6121.823
#R^2=0.141
#BIC=6715.499

model_bottom_02 <- bam(Temperature_difference ~ te(Longitude_s, Latitude_s, Temperature_s, d=c(2,1), bs=c("tp"), k=c(40, 25, 7)) + s(Time_num_s, k=5, bs="cc") + s(Depth, bs="tp", k=10),
                       data = Data_subset, method="fREML", discrete=T, nthreads=4)
gam.check(model_bottom_02)
summary(model_bottom_02)
AICc(model_bottom_02)
#AICc=6117.366
#Depth didn't add a whole lot to deviance explained, depth would be somewhat static anyways relative to lat and long
#R^2=0.142
#BIC=6711.553


#Combine temperature and time in the tensor smooth
model_bottom_03 <- bam(Temperature_difference ~  te(Latitude_s,Longitude_s,Temperature_s, Time_num_s, d=c(2,1,1), bs=c("tp", "tp","cc"), k=c(40, 25, 7, 5)),
                       data = Data_subset, method="fREML", discrete=T, nthreads=4)

summary(model_bottom_03)
plot(model_bottom_03)
gam.check(model_bottom_03)
AICc(model_bottom_03)
#AICc=5648.825
#R^2=0.213
#BIC=7034.833
#Changed R2 substantially


model_bottom_04 <- bam(Temperature_difference ~  te(Latitude_s,Longitude_s,Temperature_s, Time_num_s, d=c(2,1,1), bs=c("tp", "tp","cc"), k=c(40, 25, 7, 5))+ Year_fac,
                       data = Data_subset, method="fREML", discrete=T, nthreads=4)

summary(model_bottom_04)
gam.check(model_bottom_04)
AICc(model_bottom_04)
plot(model_bottom_04)
#AICc=5586.525
#R^2=0.221
#BIC=7024.829
#Year intercepts change R2 and AICc slightly
#Coefficients from year indicate differences in air temperature between years (2011 lower than almost all years, 2017 being higher than others)
#Perhaps consider using both Julian Day and temperature (after adjusting temperature and Julian day)


#Temperature Anomaly Model---------------------------------------------
##Julian day and temperature model
temperature_anomaly_GAM<- gam(Temperature ~ s(Julian_day_s,bs="cc",k=5),data=Data_subset)
summary(temperature_anomaly_GAM)
gam.check(temperature_anomaly_GAM)
plot(temperature_anomaly_GAM)

##Per discussion with Sam, construct model with longitude and latitude on top julian day instead
temperature_anomaly_GAM_final<- gam(Temperature ~ te(Latitude_s,Longitude_s,Julian_day_s, d=c(2,1) ,bs=c("tp","cc"),k=c(20,5)),data=Data_subset)
summary(temperature_anomaly_GAM_final)
#K index is too low at 0.72, and edf is pretty close to k', but let's go with this for now
gam.check(temperature_anomaly_GAM_final)
plot(temperature_anomaly_GAM_final)

#Add term to the dataset
Data_subset$Temperature_prediction <-predict(temperature_anomaly_GAM_final,Data_subset)
Data_subset$Temperature_anomaly <-Data_subset$Temperature - Data_subset$Temperature_prediction 
hist(Data_subset$Temperature_anomaly)

#From Sam: Look at TYPE=="TERMS" FOR PREDICT.BAM and feed newdata (give a range of julian day) - for further evaluation of temperature_anomaly_GAM_final

#Second batch of model runs---------------------------------------------
#Run models with Julian day and "temperature anomaly"

#model_bottom_01_r <- bam(Temperature_difference ~  te(Latitude_s,Longitude_s, Temperature_anomaly, Julian_day_s,Time_num_s, d=c(2,1,1,1), bs=c("tp", "tp","cc","cc"), k=c(40, 25, 7, 7, 5)),
#                       data = Data_subset, method="fREML", discrete=T, nthreads=4)
#summary(model_bottom_01_r)
#AICc(model_bottom_01_r)
#No longer relevant

model_bottom_02_r <- bam(Temperature_difference ~  te(Latitude_s,Longitude_s, Temperature_anomaly, Julian_day_s, d=c(2,1,1), bs=c("tp", "tp","cc"), k=c(40, 25, 7, 7)) + s(Time_num_s, k=5, bs="cc") ,
                         data = Data_subset, method="fREML", discrete=T, nthreads=4)
summary(model_bottom_02_r)
AICc(model_bottom_02_r)

#AICc=5510.468
#R^2=0.226
#BIC = 6799.042
#gam.check(model_bottom_02_r)
#plot(model_bottom_02_r)


model_bottom_03_r <- bam(Temperature_difference ~  te(Latitude_s,Longitude_s, Temperature_anomaly, Julian_day_s, d=c(2,1,1), bs=c("tp", "tp","cc"), k=c(40, 25, 7, 7)) ,
                         data = Data_subset, method="fREML", discrete=T, nthreads=4)
summary(model_bottom_03_r)
AICc(model_bottom_03_r)
#AICc=5494.158
#R^2=0.228
BIC(model_bottom_03_r)
#BIC = 6813.749

#Collinearity between residual temperature and Lat, Long, and Time of Day
plot(Data_subset$Residual_Temperature~Data_subset$Latitude_s)
plot(gam(Residual_Temperature ~ s(Latitude_s),data = Data_subset))

plot(Data_subset$Residual_Temperature~Data_subset$Longitude_s)
plot(gam(Residual_Temperature ~ s(Longitude_s),data = Data_subset))

plot(Data_subset$Residual_Temperature~Data_subset$Time_num_s)
plot(gam(Residual_Temperature ~ s(Time_num_s),data = Data_subset))

#Cross validation test ---------------------------------------------------------------
#with K=5
k <- 5 #the number of folds

set.seed(2020)
folds <- cvFolds(NROW(Data_subset), K=k)

Data_subset$holdoutpred <- rep(0,nrow(Data_subset))

for(i in 1:k){
  train <- Data_subset[folds$subsets[folds$which != i], ] #Set the training set
  validation <- Data_subset[folds$subsets[folds$which == i], ] #Set the validation set
  
  newlm <- bam(Temperature_difference ~  te(Latitude_s,Longitude_s, Temperature_anomaly, Julian_day_s, d=c(2,1,1), bs=c("tp", "tp","cc"), k=c(40, 25, 7, 7)) ,
               data = train, method="fREML", discrete=T, nthreads=4) #Get your new linear model (just fit on the train data)
  newpred <- predict(newlm,newdata=validation) #Get the predicitons for the validation set (from the model just fit on the train data)
  
  Data_subset[folds$subsets[folds$which == i], ]$holdoutpred <- newpred #Put the hold out prediction in the data set for later use
}

Data_subset$holdoutpred_CorrectSign<-ifelse(Data_subset$Temperature_difference==0,NA,(ifelse(Data_subset$Temperature_difference>0&Data_subset$holdoutpred>0|Data_subset$Temperature_difference<0&Data_subset$holdoutpred<0,1,0)))

summary(Data_subset$holdoutpred_CorrectSign)
#Predicted the correct sign 75% of the time (0.75)

#-----------------------------------------

#Attempt to create prediction map

#Function from Sam to create prediction dataset
WQ_pred<-function(Full_data=Data_subset,
                  Delta_subregions=Delta,
                  Delta_water=spacetools::Delta,
                  Stations = WQ_stations,
                  n=75, 
                  Temps_Anomaly=c(-2,0,2),
                  Julian_days=yday(ymd(paste("2001", c(1,4,7,10), "15", sep="-"))), #Jan, Apr, Jul, and Oct 15 for a non-leap year
                  Time_num=c(9*60*60,12*60*60,15*60*60) # 12PM x 60 seconds x 60 minutes
){
  
  # Create point locations on a grid for predictions
  Points<-st_make_grid(Delta_subregions, n=n)%>%
    st_as_sf(crs=st_crs(Delta_subregions))%>%
    st_join(Delta_water%>% # Joining a map of delta waterways (from my spacetools package) to ensure all these points are over water.
              dplyr::select(Shape_Area)%>%
              st_transform(crs=st_crs(Delta_subregions)))%>%
    filter(!is.na(Shape_Area))%>%
    st_join(Stations%>% # Applying the same approach we did to the full data: remove any pounts outside the convex hull formed by major survey stations sampled >50 times
              st_union()%>%
              st_convex_hull()%>%
              st_as_sf()%>%
              mutate(IN=TRUE),
            join=st_intersects)%>%
    filter(IN)%>%
    dplyr::select(-IN)%>%
    st_centroid()%>% # The prior grid was actually a set of polygons, this picks the center point of each
    st_transform(crs=4326)%>%
    st_coordinates()%>%
    as_tibble()%>%
    mutate(Location=1:nrow(.))%>%
    dplyr::select(Longitude=X, Latitude=Y, Location)
  
  # Create dataset for each year and season showing which subregions were sampled
  Data_effort <- Full_data%>%
    st_drop_geometry()%>%
    group_by(SubRegion, Season)%>%
    summarise(N=n())%>%
    ungroup()%>%
    left_join(Delta_subregions, by="SubRegion")%>%
    dplyr::select(-geometry)
  
  
  # Create full dataset for predictions
  newdata<-expand.grid(Temperature_anomaly= Temps_Anomaly,
                       Location=1:nrow(Points),
                       Julian_day=Julian_days,
                       Time_num=Time_num)%>% # Create all combinations of predictor variables
    left_join(Points, by="Location")%>% #Add Lat/Longs to each location
    mutate(Latitude_s=(Latitude-mean(Full_data$Latitude, na.rm=T))/sd(Full_data$Latitude, na.rm=T), # Standardize each variable based on full dataset for model
           Longitude_s=(Longitude-mean(Full_data$Longitude, na.rm=T))/sd(Full_data$Longitude, na.rm=T),
           Julian_day_s = (Julian_day-mean(Full_data$Julian_day, na.rm=T))/sd(Full_data$Julian_day, na.rm=T),
           Time_num_s=(Time_num-mean(Full_data$Time_num, na.rm=T))/sd(Full_data$Time_num, na.rm=T),
           Season=case_when(Julian_day<=80 | Julian_day>=356 ~ "Winter", # Create a variable for season
                            Julian_day>80 & Julian_day<=172 ~ "Spring",
                            Julian_day>173 & Julian_day<=264 ~ "Summer",
                            Julian_day>265 & Julian_day<=355 ~ "Fall"))%>%
    st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)%>% # Turn into sf object
    st_transform(crs=st_crs(Delta_subregions))%>% # transform to crs of Delta shapefile
    st_join(Delta_subregions, join = st_intersects)%>%
    filter(!is.na(SubRegion))%>% # Make sure all points are within our desired subregions
    left_join(Data_effort, by=c("SubRegion", "Season"))%>% # Use the Data_effort key created above to remove points in subregions that were not sampled that region, season, and year.
    filter(!is.na(N))
  return(newdata)
}

#Use the function to create prediction data frame
newdata_bottom <- WQ_pred(Full_data=Data_subset, 
                        Julian_days = yday(ymd(paste("2001", 1:12, "15", sep="-"))))


#Add prediction from model
model_03r_predictions<-predict(model_bottom_03_r, newdata=newdata_bottom, type="response", se.fit=TRUE, discrete=T, n.threads=3) # Create predictions


newdata<-newdata_bottom%>%
  mutate(Prediction=model_03r_predictions$fit)%>%
  mutate(SE=model_03r_predictions$se.fit,
         L95=Prediction-SE*1.96,
         U95=Prediction+SE*1.96)


#Create figures for each season, time of day, and residual temperature
png(filename=file.path("~/GIT Hub/WQ-discrete/Bottom_Surface_Temp_Results","Model_prediction_map_08-26-2020.png"), units="in",type="cairo", bg="white", height=18, 
    width=20, res=400, pointsize=20)
ggplot(data=newdata)+
  geom_sf(aes(colour=Prediction),pch=15)+
  scale_colour_gradient2(low = "blue",high = "red",midpoint = 0,breaks=seq(min(newdata$Prediction),max(newdata$Prediction),(max(newdata$Prediction)-min(newdata$Prediction))/5))+
  theme_dark()+
  facet_grid(Season~Temperature_anomaly)+
  theme(plot.title=element_text(size=28), 
        axis.text.x=element_text(size=21, color="black"), 
        axis.text.y = element_text(size=20, color="black"), 
        axis.title.x = element_text(size = 22, angle = 00), 
        axis.title.y = element_text(size = 22, angle = 90),
        strip.text = element_text(size = 20))+
  labs(x="Temperature Deviance from Julian Day", y="Season")
dev.off()


# Model error by region ---------------------------------------------------

mygrid <- data.frame(
  name = c("Cache Slough and Lindsey Slough", "Lower Sacramento River Ship Channel", "Liberty Island", "Suisun Marsh", "Middle Sacramento River", "Lower Cache Slough", "Steamboat and Miner Slough", "Upper Mokelumne River", "Lower Mokelumne River", "Georgiana Slough", "Sacramento River near Ryde", "Sacramento River near Rio Vista", "Grizzly Bay", "West Suisun Bay", "Mid Suisun Bay", "Honker Bay", "Confluence", "Lower Sacramento River", "San Joaquin River at Twitchell Island", "San Joaquin River at Prisoners Pt", "Disappointment Slough", "Lower San Joaquin River", "Franks Tract", "Holland Cut", "San Joaquin River near Stockton", "Mildred Island", "Middle River", "Old River", "Upper San Joaquin River", "Grant Line Canal and Old River", "Victoria Canal"),
  row = c(2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 8, 8, 8),
  col = c(4, 6, 5, 2, 8, 6, 7, 9, 9, 8, 7, 6, 2, 1, 2, 3, 4, 5, 6, 8, 9, 5, 6, 7, 9, 8, 8, 7, 9, 8, 7),
  code = c(" 1", " 2", " 3", " 8", " 4", " 5", " 6", " 7", " 9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "30", "29", "31"),
  stringsAsFactors = FALSE
)

Data_resid<-Data_subset%>%
  mutate(Residuals = model_bottom_03_r$residuals,Prediction = model_bottom_03_r$linear.predictors)

Resid_sum<-Data_resid%>%
  lazy_dt()%>%
  group_by(Year, Month, SubRegion)%>%
  summarise(Resid=mean(Residuals), SD=sd(Residuals))%>%
  ungroup()%>%
  as_tibble()

p_resid<-ggplot(Resid_sum)+
  geom_tile(aes(x=Year, y=Month, fill=Resid))+
  scale_fill_gradient2(high = muted("red"),
                       low = muted("blue"),
                       breaks=seq(-3,5.5, by=0.5),
                       guide=guide_colorbar(barheight=40))+
  scale_x_continuous(breaks=unique(Resid_sum$Year), labels = if_else((unique(Resid_sum$Year)/2)%% 2 == 0, as.character(unique(Resid_sum$Year)), ""))+
  scale_y_continuous(breaks=unique(Resid_sum$Month), labels = if_else(unique(Resid_sum$Month)%% 2 == 0, as.character(unique(Resid_sum$Month)), ""))+
  facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen())+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank(), panel.background = element_rect(fill="black"))

ggsave(plot=p_resid, file.path("~/GIT Hub/WQ-discrete/Bottom_Surface_Temp_Results","Residuals 2020-08-26.png"), device=png(), width=20, height=12, units="in")

#See how often the model predicts the right sign
#Remove data points where there's zero difference between surface and bottom temperature
Data_resid_sign<-Data_resid %>% filter(Temperature_difference!=0)

Data_resid_sign$CorrectSign<-ifelse(Data_resid_sign$Temperature_difference>0&Data_resid_sign$Prediction>0|Data_resid_sign$Temperature_difference<0&Data_resid_sign$Prediction<0,1,0)

summary(Data_resid_sign$CorrectSign)
#Mean 0.75

Data_resid_sign_sum<-Data_resid_sign%>%
  lazy_dt()%>%
  group_by(Year, Month, SubRegion)%>%
  summarise(CorrectSign=mean(CorrectSign))%>%
  ungroup()%>%
  as_tibble()



resid_sign<-ggplot(Data_resid_sign_sum)+
  geom_tile(aes(x=Year, y=Month, fill=CorrectSign))+
  scale_fill_gradient(low = "red",
                      high = "blue", breaks=seq(0,1, by=0.2),
                      guide=guide_colorbar(barheight=40))+
  scale_x_continuous(breaks=unique(Data_resid_sign_sum$Year), labels = if_else((unique(Data_resid_sign_sum$Year)/2)%% 2 == 0, as.character(unique(Data_resid_sign_sum$Year)), ""))+
  scale_y_continuous(breaks=unique(Data_resid_sign_sum$Month), labels = if_else(unique(Data_resid_sign_sum$Month)%% 2 == 0, as.character(unique(Data_resid_sign_sum$Month)), ""))+
  facet_geo(~SubRegion, grid=mygrid, labeller=label_wrap_gen())+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank(), panel.background = element_rect(fill="black"))

ggsave(plot=resid_sign, file.path("~/GIT Hub/WQ-discrete/Bottom_Surface_Temp_Results","Correct Prediction 2020-08-26.png"), device=png(), width=20, height=12, units="in")

#---------------------------------------------------------------
#Test dataset from pre-2011
test_dataset_hist<-Data %>% filter(year(Date)<2011,!is.na(Temperature_bottom))
test_dataset_hist$Temperature_difference <- test_dataset_hist$Temperature_bottom-test_dataset_hist$Temperature

test_dataset_hist$Temperature_prediction <-predict(temperature_anomaly_GAM_final,test_dataset_hist)
test_dataset_hist$Temperature_anomaly <-test_dataset_hist$Temperature - test_dataset_hist$Temperature_prediction 

test_dataset_hist$Temperature_difference_prediction<-predict(model_bottom_03_r,test_dataset_hist)

plot(test_dataset_hist$Temperature_difference~test_dataset_hist$Temperature_difference_prediction)

summary(lm(Temperature_difference~Temperature_difference_prediction,data=test_dataset_hist))
#R2 of 0.042

test_dataset_hist$Prediction_Residual <- test_dataset_hist$Temperature_difference_prediction-test_dataset_hist$Temperature_difference
plot(test_dataset_hist$Prediction_Residual~test_dataset_hist$Temperature)

#See how often the model predicts the right sign
#Remove data points where there's zero difference between surface and bottom temperature
test_dataset_sign<-test_dataset_hist %>% filter(Temperature_difference!=0)

test_dataset_sign$CorrectSign<-ifelse(test_dataset_sign$Temperature_difference>0&test_dataset_sign$Temperature_difference_prediction>0|test_dataset_sign$Temperature_difference<0&test_dataset_sign$Temperature_difference_prediction<0,1,0)

summary(test_dataset_sign$CorrectSign)
#Mean of 0.6384

#Try using model_bottom_03_r (since it has lowest BIC)
test_dataset_hist$Temperature_difference_prediction<-predict(model_bottom_03_r,test_dataset_hist)

summary(lm(Temperature_difference~Temperature_difference_prediction,data=test_dataset_hist))
#R2 of 0.05

test_dataset_hist$Prediction_Residual <- test_dataset_hist$Temperature_difference_prediction-test_dataset_hist$Temperature_difference
plot(test_dataset_hist$Prediction_Residual~test_dataset_hist$Temperature)

#See how often the model predicts the right sign
#Remove data points where there's zero difference between surface and bottom temperature
test_dataset_sign<-test_dataset_hist %>% filter(Temperature_difference!=0)

test_dataset_sign$CorrectSign<-ifelse(test_dataset_sign$Temperature_difference>0&test_dataset_sign$Temperature_difference_prediction>0|test_dataset_sign$Temperature_difference<0&test_dataset_sign$Temperature_difference_prediction<0,1,0)

summary(test_dataset_sign$CorrectSign)
#Mean of 0.6329

#---------------------------------------------------------------
#Regression Tree exploration
library(rpart)

rtree_01<- rpart(Temperature_difference ~Latitude_s + Longitude_s  + Temperature_s + Time_num_s,data=Data_subset,method="anova")

plot(rtree_01)
text(rtree_01, digits = 3)
printcp(rtree_01)
plotcp(rtree_01) # visualize cross-validation results
summary(rtree_01) # detailed summary of splits
rsq.rpart(rtree_01)


rtree_02<- rpart(Temperature_difference ~Latitude_s + Longitude_s  + Julian_day_s + Residual_Temperature+ Time_num_s,data=Data_subset,method="anova")
plot(rtree_02)
text(rtree_02, digits = 3)
plotcp(rtree_02) # visualize cross-validation results
summary(rtree_02) # detailed summary of splits
rsq.rpart(rtree_02)


Data_subset$Temperature_difference_code<-ifelse(Data_subset$Temperature_difference<(-0.5),"Colder by 0.5 C",ifelse(Data_subset$Temperature_difference>0.5,"Warmer by 0.5 C","Within 0.5 C of surface"))
Data_subset$Temperature_difference_code<-as.factor(Data_subset$Temperature_difference_code)

rcart_01<- rpart(Temperature_difference_code ~Latitude_s + Longitude_s  + Julian_day_s + Residual_Temperature+ Time_num_s,data=Data_subset,method="class")
plot(rcart_01)
summary(rcart_01)
text(rcart_01, digits = 3)


rcart_02<- rpart(Temperature_difference_code ~Latitude_s + Longitude_s  + Julian_day_s ,data=Data_subset,method="class")
plot(rcart_02,uniform=TRUE)

rcart_03<- rpart(Temperature_difference_code ~Latitude_s + Longitude_s ,data=Data_subset,method="class")
plot(rcart_03,uniform=TRUE)

##################################################################
