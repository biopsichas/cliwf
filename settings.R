##------------------------------------------------------------------------------
## Settings to set paths to directories, files and RCPs & RCMs & Periods
##------------------------------------------------------------------------------
## !!!MAKE YOUR MODIFICATIONS ONLY IN THIS SECTION!!!

## Cores, which could to be used 
cores <- 10 ## if NULL max - 1 will be used
## Path to the SWAT+ setup database (sqlite)
db_path <- "data/cs6_project.sqlite"
## Define path where climate data are saved
## Data location in on OPTAIN Cloud WPs&Tasks>WP3>Task_3_2>local data>v1a>data
cli_dir <- "data/climate_data"
# Define the path to your SWAT project txtinout folder (after hard calibration)
setup_dir <- 'data/clean_setup'
# Define the path to your management schedule input file (from Micha's SWATfarmR
## input script)
mgt <- 'data/farmR_input.csv'
## Directory to save results
tmp_path <- "tmp"
## Define path to calibration.cal file
cal_file <- "data/calibration.cal" ## if NULL, no calibration.cal will be included

## names of RCP scenarios (same as folder names)
rcp <- c("rcp26", "rcp45", "rcp85")
## names of RCM used (same as folder names)
rcm <- c("1", "2", "3", "4", "5", "6")
## Periods to be prepared and used in modelling
periods <- list(c("H", "1988-01-01", "2020-12-31"),
                c("N", "2033-01-01", "2065-12-31"),
                c("E", "2066-01-01", "2098-12-31"))

## SWAT excutable name
swat_exe <- 'SWATp_jan_sept.exe'

## Outflow reach
outflow_reach <- 50

## Crop selection
crop_sel <- c("corn", "wbar", "csil", "fesc", "wwht", "soyb", "canp", "grap")


## Grain units

#If you want to use grain units to normalize the basin wide sum of crop yields
#by crop-specific nutritional values, please specify grain units for relevant crops
# The grain units must be applicable to dry mass!!!

## Values are provided on: 
## OPTAIN Cloud>WPs&Tasks>WP4>Tools to share>OPTAIN_crops_drymass_grain_units_v2.xlsx

grain_units <- data.frame('wbar' = 1.163, 
                          'csil' = 1.071, 
                          'wwht' = 1.209, 
                          'fesc' = 0.718,
                          'corn' = 1.071,
                          'soyb' = 1, 
                          'canp' = 1.3, 
                          'grap' = 1)

## Thresholds for nutrient and sediment concentrations for output analysis

# thresholds for nitrogen concentration (mg N/l) and phosphorus concentration (mg P/l)
# the number of days beyond these thresholds will later be calculated
# default value are the respective median values of reported threshold values 
# for very small siliceous rivers in lowland across Europe 
# (https://www.sciencedirect.com/science/article/pii/S0048969719338380)
# please check if this is appropriate for your case study (e.g. type of river)
# feel free to use other threshold values!!

threshold_N=2.3
threshold_P=0.082 

# threshold for sediment concentration (mg N/l) 
# the number of days beyond this threshold will later be calculated
# default value is 50 mg/l (missing reference), 
# if you know a reference please let me know (michael.strauch@ufz.de)
# feel free to use another threshold value!!

threshold_Sed=50

