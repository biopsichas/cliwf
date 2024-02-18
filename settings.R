##------------------------------------------------------------------------------
## Settings to set paths to directories, files and RCPs & RCMs & Periods
##------------------------------------------------------------------------------
## !!!MAKE YOUR MODIFICATIONS ONLY IN THIS SECTION!!!

## Path to the SWAT+ setup database (sqlite)
db_path <- "D:/_WORKFLOW/Temp/buildr_project/cs4_project/cs4_project.sqlite"
## Define path where climate data are saved
## Data location in on OPTAIN Cloud WPs&Tasks>WP3>Task_3_2>local data>v1a>data
cli_dir <- "C:/Users/Svajunas/Documents/climate_data"
# Define the path to your SWAT project txtinout folder
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