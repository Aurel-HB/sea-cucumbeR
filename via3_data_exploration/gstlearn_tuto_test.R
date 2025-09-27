# this code explore the structure of the PP of the survey HOLOTVSPM2025
# the package gstlearn is used to conduct the basic approach

library(gstlearn)
library(here)

#Global variables
verbose = TRUE
graphics = TRUE
err = OptCst_define(ECst_NTCOL(),6)


##############################
# Work on one Point Process
##############################

############
#Reading data
###########

data_position <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/data_position_2025.rds",
  sep=""))

write.csv(data_position %>% select(X,Y,station) %>% filter(station==133),
          paste(
            here(),"/data_position_2025.csv",
            sep=""), row.names = FALSE)

filepath = loadData("data_position_2025", "data_position_2025.csv")
mydb = Db_createFromCSV(filepath,CSVformat())
err = mydb$setLocator("X",ELoc_X(),0)
err = mydb$setLocator("Y",ELoc_X(),1)
err = mydb$setLocator("station",ELoc_Z())
if (verbose){
  dbfmt = DbStringFormat()
  dbfmt$setFlags(flag_extend = TRUE)
  mydb$display(dbfmt)
}
#Accessing to the variable names
cat("List of all variable names =",mydb$getAllNames())

#Extracting the vector containing the station variable in order to perform a selection

tabstation = mydb$getColumn("station")
selstation = as.numeric(tabstation == 133)
mydb$addSelection(selstation,"sel")


if (verbose){
  mydb$display()
}

#Display my Data, we used the selection sel because error if not

if (graphics){
  plot.init() + plot(mydb,nameColor="sel") + plot.decoration(title="Data Set")}


############
#Variograms
###########

#We first define the geometry of the variogram calculations
myVarioParamOmni = VarioParam()
mydir = DirParam_create(nlag=700,dlag=1)
myVarioParamOmni$addDir(mydir)

#We use the variogram definition in order to calculate the variogram cloud.
dbcloud = db_vcloud(db=mydb, varioparam=myVarioParamOmni)

#Representing the variogram cloud.
if (graphics){
  plot.init() + plot(dbcloud,name="Cloud*") + 
    plot.decoration(title="Variogram Cloud")
}

#Calculating the experimental omni-directional variogram
myVarioOmni = Vario(myVarioParamOmni)
err = myVarioOmni$compute(mydb, ECalcVario_VARIOGRAM())
if (verbose){
  myVarioOmni$display()
}

#The variogram is represented graphically for a quick check
if (graphics){
  plot.init() + plot.varmod(myVarioOmni) + 
    plot.decoration(title="Omni-directional Variogram for Pb")}

#Calculate a variogram in several directions
myvarioParam = VarioParam()
mydirs = DirParam_createMultiple(ndir=4, nlag=10, dlag=1.)
myvarioParam$addMultiDirs(mydirs)

myvario = Vario(myvarioParam)
myvario$compute(mydb, ECalcVario_VARIOGRAM())

if (verbose){
  myvario$display()}

if (graphics){
  plot.init() + plot.varmod(myvario) + 
    plot.decoration(title="Multi-Directional Variogram of Pb")
}

#Calculating the Variogram Map
myvmap = db_vmap(db=mydb,calcul_type=ECalcVario_VARIOGRAM(),nxx=c(20,20))
if (verbose){
  myvmap$display()
}

if (graphics){
  plot.init() + plot(myvmap, name="*Var") + 
    plot.decoration(title="Variogram Map") 
}

