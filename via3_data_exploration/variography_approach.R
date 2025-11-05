# This code create a first approach to variography with the HOLOTVSPM2025 data

library(gstlearn)
library(ggplot2)
library(here)
library(dplyr)

## Load data ####
data_abun <- readRDS(
  paste(here(),"/via3_data_exploration/Data/processed/data_abun_2025.rds",
        sep=""))

data_position <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/data_position_2025.rds",
  sep=""))

## Prepare data ####
filecsv <- "Density_2025.csv"
# creating Db oject for gstlearn analysis
csv <- CSVformat_create(flagHeader=TRUE, naString = "MISS")
dat.abun <- Db_createFromCSV(paste(here(),"/comparing_biomass_method/Data/",
                                   filecsv,sep=""), csv=csv)
dat.abun$setLocators(c("X","Y"), ELoc_X())
dat.abun$setLocator("intensity", ELoc_Z(), cleanSameLocator=TRUE)

dat.position <- Db_createFromCSV(paste(
  here(),"/via3_data_exploration/Data/processed/data_position_2025_133.csv",
  sep=""), csv=csv)
dat.position$setLocators(c("X","Y"), ELoc_X())
dat.position$setLocator("station", ELoc_Z(), cleanSameLocator=TRUE)

## work on the density per line ####
#Variogram Cloud
varioParamOmni = VarioParam_createOmniDirection(nlag = 100)
grid.cloud = db_vcloud(dat.abun, varioParamOmni)
grid.cloud$display()

p = plot.init()
p = p + plot.raster(grid.cloud, "Cloud.intensity*")
p = p + plot.geometry(asp=0)
plot.end(p)

#Experimental Variogram
varioParamOmni = VarioParam_createOmniDirection(nlag=40,dlag=1852,toldis = 0.1) 

varioexp = Vario(varioParamOmni)
err = varioexp$compute(dat.abun)
varioexp # print the contents of the variogram
plot.init() + plot.varmod(varioexp) # plot the variogram

plot.init() + plot.varmod(varioexp,drawPsize=TRUE)

# Automatic Model Fitting
fitmod = Model()
err = fitmod$fit(varioexp)
plot.init() + plot.varmod(varioexp, fitmod)
fitmod

# Model Fitting with pre-defined basic structures
ECov_printAll() #all  basic covariance structure
types = ECov_fromKeys(c("NUGGET","GAUSSIAN","SPHERICAL","MATERN"))
err = fitmod$fit(varioexp, types=types)
plot.init() + plot.varmod(varioexp, fitmod)
fitmod

# Model Fitting with constraints
types = ECov_fromKeys(c("NUGGET","CUBIC","SPHERICAL"))
constraints = Constraints()
err = constraints$addItemFromParamId(EConsElem_RANGE(),icov=1,
                                     type=EConsType_UPPER(),value=20.)
err = constraints$addItemFromParamId(EConsElem_SILL(),icov=1,
                                     type=EConsType_LOWER(),value=0.03)
err = fitmod$fit(varioexp, types=types, constraints=constraints,
                 optvar=Option_VarioFit(TRUE))

plot.init() + plot.varmod(varioexp, fitmod)
fitmod

constraints = Constraints() # equamity constraints
err = constraints$addItemFromParamId(EConsElem_RANGE(),icov=1,
                                     type=EConsType_EQUAL(),value=1000.)
err = constraints$addItemFromParamId(EConsElem_SILL(),icov=1,
                                     type=EConsType_EQUAL(),value=0.4)
err = fitmod$fit(varioexp, types=types, constraints=constraints,
                 optvar=Option_VarioFit(TRUE))

plot.init() + plot.varmod(varioexp, fitmod)
fitmod

# Directional Variograms
varioParamMulti = VarioParam_createMultiple(ndir=4, nlag=40, dlag=1852)
vario.4dir = Vario(varioParamMulti)
err = vario.4dir$compute(dat.abun)

plot.init() + plot.varmod(vario.4dir)

model.4dir = Model() # fit model
err = model.4dir$fit(vario.4dir,types=types)

plot.init() + plot.varmod(vario.4dir, model.4dir)

# Experimental variogram map
grid.vmap = db_vmap(dat.abun)

p1 = plot.init() + plot.raster(grid.vmap, flagLegend=TRUE, 
                               legendName="density")
p2 = plot.init() + plot.raster(grid.vmap, name="VMAP.intensity.Nb",
                               flagLegend=TRUE, legendName="# pairs")
ggarrange(p1,p2,nrow=1,ncol=2)

modelVM = Model() #fit a model directly on the experimental variogram map
err = modelVM$fitFromVMap(grid.vmap, types=types)
modelVM
