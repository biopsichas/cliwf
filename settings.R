##------------------------------------------------------------------------------
## Settings to set paths to directories, files and RCPs & RCMs & Periods
##------------------------------------------------------------------------------
## !!!MAKE YOUR MODIFICATIONS ONLY IN THIS SECTION!!!

## Path to the SWAT+ setup database (sqlite)
db_path <- "D:/_WORKFLOW/Temp/buildr_project/cs4_project/cs4_project.sqlite"
## Define path where climate data are saved
## Data location in on OPTAIN Cloud WPs&Tasks>WP3>Task_3_2>local data>v1a>data
cli_dir <- "C:/Users/Svajunas/Documents/climate_data"
# Define the path to your SWAT project txtinout folder (after hard calibration)
setup_dir <- 'D:/Zglowiaczka/clean_setup4'
# Define the path to your management schedule input file (from Micha's SWATfarmR
## input script)
mgt <- 'D:/_WORKFLOW/Temp/farmR_input/farmR_input.csv'
##Directory to save results
tmp_path <- "tmp"

## names of RCP scenarios (same as folder names)
rcp <- c("rcp26", "rcp45", "rcp85")
## names of RCM used (same as folder names)
rcm <- c("1", "2", "3", "4", "5", "6")
## Periods to be prepared and used in modelling
periods <- list(c("H", "1988-01-01", "2020-12-31"),
                c("N", "2033-01-01", "2065-12-31"),
                c("E", "2066-01-01", "2098-12-31"))

## SWAT excutable name
swat_exe <- 'Rev_61_0_64rel.exe'

## Outflow reach
outflow_reach <- 43

## Crop selection
crop_sel <- c("wwht","trit","canp","barl","csil","corn","sgbt","onio","mint",
              "lett","crrt","fesc","alfa") 

## Grain units

#If you want to use grain units to normalize the basin wide sum of crop yields
#by crop-specific nutritional values, please specify grain units for relevant crops
# The grain units must be applicable to dry mass!!!

grain_units <- data.frame('alfa' = 1.163, 'csil' = 1.071, 'wwht' = 1.209, 
                          'barl' = 1.163, 'sgbt' = 1, 'canp' = 1, 'corn' = 1.071,
                          'crrt' = 1, 'fesc' = 1, 'lett' = 1, 'mint' = 1, 'onio' = 1,
                          'trit' =1)

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

