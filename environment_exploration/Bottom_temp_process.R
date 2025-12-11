#This code extract the bottom temperature from a NetCDF file and process the
# data to keep the period between 2021 and 2025 in the study area : 
# (45.849N, 45.3656S, -56.1402E, -56.26W)

# load the package
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
                      "/environment_exploration/Environment_Data/row/",
                      "cmems_mod_glo_phy_my_0.083deg_P1M-m_1765382048849.nc",
                      sep=""))
print(ncin)

# extract dimension of the file ####
# get longitude and latitude
lon <- ncvar_get(ncin,"longitude")
nlon <- dim(lon)
head(lon)

lat <- ncvar_get(ncin,"latitude")
nlat <- dim(lat)
head(lat)

# get time
time <- ncvar_get(ncin,"time")
time

tunits <- ncatt_get(ncin,"time","units")
tunits
nt <- dim(time)

# get depth
depth <- ncvar_get(ncin,"depth")
ndepth <- dim(depth)
head(ndepth)


# get variable sea floor temperature
variable <- "bottomT"
tmp_array <- ncvar_get(ncin,variable)
dlname <- ncatt_get(ncin,variable,"long_name")
dunits <- ncatt_get(ncin,variable,"units")
fillvalue <- ncatt_get(ncin,variable,"_FillValue")
dim(tmp_array)

# Reshaping from raster to rectangular ####

#Convert the time variable
#convert time to CFtime class
cf <- CFtime(tunits$value, calendar = "gregorian", time)

timestamps <- as_timestamp(cf) # get character-string times
timestamps
class(timestamps)

#parse the string into date components
time_cf <- parse_timestamps(cf, timestamps)
time_cf
class(time_cf)

