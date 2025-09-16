# this script is a tutorial of the package gstlearn
library(gstlearn)

verbose = TRUE
graphics = TRUE
err = OptCst_define(ECst_NTCOL(),6)

#reading data
#dlfile = "https://soft.minesparis.psl.eu/gstlearn/data/Pollution/Pollution.dat"
#filepath = "Pollution.dat"
#download.file(dlfile, filepath, quiet=TRUE)
filepath = loadData("Pollution", "Pollution.dat")
mydb = Db_createFromCSV(filepath,CSVformat())
err = mydb$setLocator("X",ELoc_X(),0)
err = mydb$setLocator("Y",ELoc_X(),1)
err = mydb$setLocator("Zn",ELoc_Z())
if (verbose)
{
  dbfmt = DbStringFormat()
  dbfmt$setFlags(flag_extend = TRUE)
  mydb$display(dbfmt)
}

#Accessing to the variable names
cat("List of all variable names =",mydb$getAllNames())

#Extracting the vector containing the Zn variable in order to perform a selection

tabZn = mydb$getColumn("Zn")
selZn = as.numeric(tabZn < 20)
mydb$addSelection(selZn,"sel")

mydb$setLocator('Pb',ELoc_Z())

if (verbose){
   mydb$display()
}
 
#Display my Data (with samples represented by color and size)

if (graphics){
  plot.init() + plot(mydb,nameColor="Pb") + plot.decoration(title="Data Set")}



