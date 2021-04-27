

require(sp)
require(gstat)
require(spacetime)
require(dplyr)
require(tidyr)
require(stringr)
require(dtplyr)
require(mgcv)
require(ggplot2)
require(geofacet)
require(lubridate)
require(hms)
require(scales)
require(patchwork)
require(sf)
require(stars)
require(itsadug)
require(purrr)
require(readr)
require(slider)
require(colorspace)
require(ggstance)
source("Utility_functions.R")

require(discretewq)
options(scipen=999)

#Set up a path for datasets and shapefiles
data_root<-file.path("data-raw")
results_root<-file.path("results")

# Create overall dataset --------------------------


# Load Delta Shapefile from Brian
Delta<-st_read(file.path(data_root,"Delta subregions","EDSM_Subregions_03302020.shp"))%>%
  filter(!SubRegion%in%c("Upper Yolo Bypass"))%>% # Remove regions outside our domain of interest
  dplyr::select(SubRegion)
# Visualize regions
ggplot(data=Delta,aes(label=SubRegion))+geom_sf()+geom_sf_text()

# Load data from 'discretewq' package
Data <- wq()%>%
  filter(!is.na(Salinity) & !is.na(Datetime) & !is.na(Latitude) & !is.na(Longitude) & !is.na(Date) & !is.na(Conductivity))%>% #Remove any rows with NAs in our key variables
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

#add dayflow data

dayflow<-read_csv(file.path(data_root, "Dayflow","dayflow-results-1997-2020.csv"), col_types = cols_only(Date="c", TOT="d", PREC="d"))%>%
  bind_rows(read_csv(file.path(data_root,"Dayflow", "dayflow-results-1984-1996.csv"), col_types = cols_only(Date="c", TOT="d", PREC="d")))%>%
  bind_rows(read_csv(file.path(data_root,"Dayflow", "dayflow-results-1970-1983.csv"), col_types = cols_only(Date="c", TOT="d", PREC="d")))%>%
  bind_rows(read_csv(file.path(data_root,"Dayflow", "dayflow-results-1956-1969.csv"), col_types = cols_only(Date="c", TOT="d", PREC="d")))%>%
  mutate(Date=parse_date_time(Date, "%m/%d/%Y", tz = "America/Los_Angeles"),
         PREC_acrefeetperday=(PREC*3600*24)/(43560),
         PREC_feetperday=if_else(Date>=parse_date_time("10/01/1980", "%m/%d/%Y", tz = "America/Los_Angeles"), PREC_acrefeetperday/682230, PREC_acrefeetperday/738000),
         PREC_CFS=(PREC_feetperday*682230*43560)/(3600*24))%>%
  arrange(Date)%>%
  mutate(TOT_mean30=slide_index_dbl(.$TOT, .$Date, .before=days(30), .f=mean, .complete = T),
         PREC_mean30=slide_index_dbl(.$PREC, .$Date, .before=days(30), .f=mean, .complete = T),
         Month=month(Date))%>%
  group_by(Month)%>%
  mutate(across(c(TOT_mean30, PREC_mean30), list(sd= ~sd(.x, na.rm=T), mean=~mean(.x, na.rm=T))))%>%
  ungroup()%>%
  select(-Month)



Data_edit <-Data%>%
  mutate(WY=Year)%>%
  mutate(WY=if_else(Month%in%10:12, WY+1, WY))%>%
  left_join(dayflow, by="Date")%>%
  filter(!is.na(TOT))%>%
  arrange(Station, Date)%>%
  group_by(Station)%>%
  mutate(Lag=Date-lag(Date, order_by = Date))%>%
  ungroup()%>%
  mutate(Start=if_else(is.na(Lag) | Lag>3600*24*60, TRUE, FALSE),
         Series_ID=1:n(),
         Series_ID=if_else(Start, Series_ID, NA_integer_),
         Series_ID=as.integer(as.factor(Series_ID)))%>%
  fill(Series_ID, .direction="down")%>%
  mutate(across(c(TOT_mean30, PREC_mean30), list(s=~(.x-mean(.x))/sd(.x))))%>%
  mutate(WY_s=(WY-mean(unique(WY)))/sd(unique(WY)),
         TOT_mean30_s_month=(TOT_mean30-TOT_mean30_mean)/TOT_mean30_sd,
         PREC_mean30_s_month=(PREC_mean30-PREC_mean30_mean)/PREC_mean30_sd)

str(Data_edit)

# Now filter data to only include this final set of subregions, and any stations outside the convex hull formed by the >25 samples stations from the major monitoring programs
Data_edit<-Data_edit%>%
  mutate_at(vars(Date_num, Longitude, Latitude, Time_num, Year, Julian_day), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T))) # Create centered and standardized versions of covariates


test_gam <- bam(Salinity ~ te(Latitude_s, Longitude_s, TOT_mean30_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13)) + 
                             te(Latitude_s, Longitude_s, TOT_mean30_s, d=c(2,1), bs=c("tp", "cc"), k=c(25, 13), by=WY_s) + 
                             s(Time_num_s, k=5), family=scat, data = Data_edit, method="fREML", discrete=T, nthreads=3)
summary(test_gam)
gam.check(test_gam)
plot(test_gam)


CC_pred<-predict(test_gam, newdata=Data_edit, type="terms", se.fit=TRUE, discrete=T, n.threads=4)
Data_edit$slope<-CC_pred$fit[,"te(TOT_mean30_s,Latitude_s,Longitude_s):WY_s"]

utils::View(Data_edit)

#Check out package 'rtide'

#New data
newdata<-Data_edit%>%
  mutate(Slope=test_gam$fit[,"te(TOT_mean30_s,Latitude_s,Longitude_s):WY_s"],
         Slope_se=test_gam$se.fit[,"te(TOT_mean30_s,Latitude_s,Longitude_s):WY_s"],
         Intercept=test_gam$fit[,"te(TOT_mean30_s,Latitude_s,Longitude_s)"]+test_gam$fit[,"s(Time_num_s)"])

test_gam$fit[,"te(TOT_mean30_s,Latitude_s,Longitude_s):WY_s"]
gam.fit(test_gam)
