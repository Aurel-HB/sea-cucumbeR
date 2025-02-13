#load the ncdf4 package 
library(ncdf4)
library(here)
library(CFtime)
library(lattice)
library(RColorBrewer)
library(sf)
library(viridis)
library(ggplot2)

#open a netCDF file 
ncin <- nc_open(paste(here(),
          "/environment_exploration/Environment_Data/cmems_mod_glo_phy_my_0.083deg.nc",
          sep=""))
print(ncin)

# get longitude and latitude
lon <- ncvar_get(ncin,"longitude")
nlon <- dim(lon)
head(lon)

lat <- ncvar_get(ncin,"latitude")
nlat <- dim(lat)
head(lat)

print(c(nlon,nlat))

# get time
time <- ncvar_get(ncin,"time")
time

tunits <- ncatt_get(ncin,"time","units")
tunits

nt <- dim(time)
nt

# get variable
variable <- "bottomT"
#variable <- "thetao"
tmp_array <- ncvar_get(ncin,variable)
dlname <- ncatt_get(ncin,variable,"long_name")
dunits <- ncatt_get(ncin,variable,"units")
fillvalue <- ncatt_get(ncin,variable,"_FillValue")
dim(tmp_array)

# get global attributes
title <- ncatt_get(ncin,0,"title")
institution <- ncatt_get(ncin,0,"institution")
datasource <- ncatt_get(ncin,0,"source")
references <- ncatt_get(ncin,0,"references")
history <- ncatt_get(ncin,0,"history")
Conventions <- ncatt_get(ncin,0,"Conventions")

# Reshaping from raster to rectangular ####
# load some packages
#library(lattice)
#library(RColorBrewer)

#Convert the time variable
# decode time
#convert time to CFtime class
cf <- CFtime(tunits$value, calendar = "gregorian", time)
cf

timestamps <- CFtimestamp(cf) # get character-string times
timestamps
class(timestamps)

#parse the string into date components
time_cf <- parse_timestamps(cf, timestamps)
time_cf
class(time_cf)

#Replace netCDF fillvalues with R NAs
# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA
length(na.omit(as.vector(tmp_array[,,1])))

#Get a single time slice of the data, create an R data frame,
#and write a .csv file
# get a single slice or layer (January)
m <- 1
tmp_slice <- tmp_array[,,m]
dim(tmp_slice)

# quick map
image(lon,lat,tmp_slice, col=rev(brewer.pal(10,"RdBu")))

# levelplot of the slice
grid <- expand.grid(lon=lon, lat=lat)
cutpts <- seq(-6,10,2)
levelplot(tmp_slice ~ lon * lat, data=grid, at=cutpts, cuts=11, pretty=T, 
          col.regions=(rev(brewer.pal(10,"RdBu"))))

# create dataframe -- reshape data
# matrix (nlon*nlat rows by 2 cols) of lons and lats
lonlat <- as.matrix(expand.grid(lon,lat))
dim(lonlat)
# vector of `tmp` values
tmp_vec <- as.vector(tmp_slice)
length(tmp_vec)

# create dataframe and add names
tmp_df01 <- data.frame(cbind(lonlat,tmp_vec))
names(tmp_df01) <- c("lon","lat",paste("Sea_floor_tmp",as.character(m),
                                       sep="_"))
head(na.omit(tmp_df01), 10)

# set path and filename
csvpath <- paste(here(),
                 "/environment_exploration/Environment_Data/",
                 sep="")
csvname <- "cru_tmp_1.csv"
csvfile <- paste(csvpath, csvname, sep="")
write.table(na.omit(tmp_df01),csvfile, row.names=FALSE, sep=",")


#Convert the whole array to a data frame, and calculate MTWA, 
#MTCO and the annual mean ####
# reshape the array into vector
tmp_vec_long <- as.vector(tmp_array)
length(tmp_vec_long)

# reshape the vector into a matrix
tmp_mat <- matrix(tmp_vec_long, nrow=nlon*nlat, ncol=nt)
dim(tmp_mat)
head(na.omit(tmp_mat))

# create a dataframe
lonlat <- as.matrix(expand.grid(lon,lat))
tmp_df02 <- data.frame(cbind(lonlat,tmp_mat))
names(tmp_df02) <- c("lon","lat","tmpJan","tmpFeb","tmpMar","tmpApr","tmpMay",
                     "tmpJun","tmpJul","tmpAug","tmpSep","tmpOct","tmpNov",
                     "tmpDec")
# options(width=96)
head(na.omit(tmp_df02, 20))
# create a dataframe without missing values
tmp_df03 <- na.omit(tmp_df02)
head(tmp_df03)

# transform dataframe in sf
#library(sf)
sea_floor_tmp <- st_as_sf(tmp_df03, crs = "EPSG:4326", coords = c("lon","lat"))
saveRDS(sea_floor_tmp, file = paste(here(), "/SIG/SIG_Data/sea_floor_tmp.rds",
        sep=""))

#map the tmp
ggplot(sea_floor_tmp['tmpJan'])+
  geom_sf(aes(color=tmpJan))+
  scale_color_viridis()+
  theme(aspect.ratio = 1,
        legend.title = element_blank(),
        title = element_text(color = "black",face = "bold"),
        plot.title = element_text( size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 8,hjust = 0.5),
        panel.border = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"),
        panel.background = element_rect(fill = "#d0d1e6"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"))+
  coord_sf(xlim = c(-58,-54), ylim = c(43,48), expand = FALSE)+
  xlab("")+ylab("")+
  labs(title = "Banc de Saint Pierre et Zone d'étude du projet")

# get the annual mean and MTWA and MTCO
tmp_df02$mtwa <- apply(tmp_df02[3:14],1,max) # mtwa
tmp_df02$mtco <- apply(tmp_df02[3:14],1,min) # mtco
tmp_df02$mat <- apply(tmp_df02[3:14],1,mean) # annual (i.e. row) means
head(na.omit(tmp_df02))
