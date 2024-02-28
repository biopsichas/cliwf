##------------------------------------------------------------------------------
## 1) Load (or install&load) libraries
##------------------------------------------------------------------------------

# library(remotes)
# !!!Please install SWATfarmR included in the workflow_script folder 
# Rstudio/Tools/Install packages/Install from Package Archive File
# remotes::install_github('chrisschuerz/SWATrunR@dev_run_scenario')
# remotes::install_github("biopsichas/SWATprepR")

##Load libraries
library(SWATfarmR)
library(SWATprepR)
library(SWATrunR)
library(tidyverse)
library(stringr)
library(ggpubr)
library(purrr)
library(doParallel)
library(foreach)
library(DBI)
library(RSQLite)

source('settings.R')
source('lib/functions.R')

##------------------------------------------------------------------------------
## 2) Prepare fresh management files
##------------------------------------------------------------------------------
## Saving reference setup files
tmp_setup_path <- paste0(tmp_path, "/", "setup")
## Saving reference result files
tmp_result_path <- paste0(tmp_path, "/", "results")

## !!!If the directory exists, delete the results directory. (Please be careful!!!)
if (file.exists(tmp_path)) unlink(tmp_path, recursive = TRUE)
dir.create(tmp_setup_path, recursive = TRUE)

## This is needed for write.exe to work (if it is right after SWATbuildR).
## Just writing a dot in project_config table.
db <- dbConnect(RSQLite::SQLite(), db_path)
project_config <- dbReadTable(db, 'project_config')
project_config$input_files_dir <- "."
dbWriteTable(db, 'project_config', project_config, overwrite = TRUE)
dbDisconnect(db)

## Coping.sqline database to tmp directory
file.copy(db_path, tmp_setup_path)

## Write out text files out of the database
exe_copy_run("lib", tmp_setup_path, "write.exe")

## Deleting unnecessary files
unlink(paste0(tmp_setup_path, "/", 
              setdiff(list.files(tmp_setup_path), c("plant.ini", "hru-data.hru", 
                                              "landuse.lum", "management.sch"))))

## Copying all your current setup files to tmp directory
file.copy(c(setdiff(list.files(path = setup_dir, full.names = TRUE), 
                    list.files(path = setup_dir, 
                               pattern = ".*.txt|.*.zip|.*success.fin|.*co2.out|.*simulation.out|.*plant.ini|.*hru-data.hru|.*landuse.lum|.*management.sch|.*.bak|.*.mgts|.*.farm|.*area_calc.out|.*checker.out|.*sqlite|.*diagnostics.out|.*erosion.out|.*files_out.out|.*.swf|.*.pcp|.*.tmp|.*.hmd|.*.slr|.*.wnd|.*.cli", full.names = TRUE)), paste(setup_dir, "soil_plant.ini", sep = "/")), 
          tmp_setup_path)

##------------------------------------------------------------------------------
## 3) Updating landuse.lum file
##------------------------------------------------------------------------------

## Backing up landuse.lum file
if(!file.exists(paste0(tmp_setup_path, "/", "landuse.lum.bak"))) {
  file.copy(from = paste0(tmp_setup_path, "/", "landuse.lum"),
            to = paste0(tmp_setup_path, "/", "landuse.lum", ".bak"), overwrite = TRUE)
}

## Updating it
source('lib/read_and_modify_landuse_lum.R')

##------------------------------------------------------------------------------
## 4) Run climate data preparing part
##------------------------------------------------------------------------------

## Files to update
files_to_update <- c("aquifer.con", "chandeg.con", "hru.con", "reservoir.con",
                     "rout_unit.con", "time.sim")

## If the directory exists, delete the results directory. 
if (file.exists(tmp_result_path)) {
  unlink(tmp_result_path, recursive = TRUE)
}
## Write results directory
dir.create(tmp_result_path, recursive = TRUE)

## Prepare parallelization
cores <- detectCores() - 1
cl <- makeCluster(cores,  outfile="")
registerDoParallel(cl)

## Run parallelization for climate data preparation
invisible(foreach(rcp = rcp, .packages = c("SWATprepR", "purrr")) %:% 
            foreach(rcm = rcm) %dopar% {write_cli(rcp, rcm)})

##Clean clusters after
stopCluster(cl)


##------------------------------------------------------------------------------
## 5) Prepare for SWATfarmR
##------------------------------------------------------------------------------

##Copy files from setup directory required for farmR and SWATrunR
m_dir <- list.dirs(tmp_result_path, recursive = TRUE)[-1]

##Copy files from setup directory required for SWATfarmR
walk(m_dir, ~file.copy(paste0(tmp_setup_path, "/",
                              c("file.cio", "plant.ini", "hru-data.hru", "landuse.lum",
                                "management.sch", "rout_unit.def", "rout_unit.ele",
                                "soils.sol", "topography.hyd", "fertilizer.frt",
                                "plants.plt",
                                "tillage.til")), .x), overwrite = FALSE)

##------------------------------------------------------------------------------
## 6) Run parallelized SWATfarmR calculation (very long)
##------------------------------------------------------------------------------

##Prepare parallelization
cl <- makeCluster(cores,  outfile="")
registerDoParallel(cl)

##Run parallelization
txt_info <- foreach(d = m_dir, .packages = c("SWATfarmR", "tidyverse",
                                 "stringr")) %dopar% {write_mgt(d, periods)}
##Clean after
stopCluster(cl)

##------------------------------------------------------------------------------
## 7) Make sure the right outputs is printed
##------------------------------------------------------------------------------

## Remove the SWIFT directories, which might be created by SWAT runs
m_dir <- m_dir[!endsWith(m_dir, "SWIFT")]

## Input for the print file
print_prt <- read_lines(paste0(tmp_setup_path,'/print.prt'), lazy = FALSE)
print_prt <- gsub(" y ", " n ", print_prt)

print_prt[10] <- "basin_wb                     n             y             n             y  "
print_prt[14] <- "basin_pw                     n             y             n             y  "
print_prt[15] <- "basin_aqu                    n             y             n             y  "
print_prt[33] <- "hru_wb                       n             y             n             y  "
print_prt[34] <- "hru_nb                       n             n             n             y  "
print_prt[35] <- "hru_ls                       n             n             n             y  "
print_prt[36] <- "hru_pw                       n             n             n             y  "
print_prt[42] <- "channel_sd                   n             n             n             y  "
print_prt[44] <- "reservoir                    n             n             n             y  "

write_lines(print_prt, paste0(tmp_setup_path,'/print.prt'))

## Updating in the result directories
overwrite_file("print.prt")

## Write object.prt
obj_prt_path <- paste0(tmp_setup_path,'/object.prt')
unlink(obj_prt_path )
write.table(paste0("object.prt", ": written by the climate-workflow on ", Sys.time()), 
            obj_prt_path, append = FALSE, sep = "\t", 
            dec = ".", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(paste(sprintf(c(rep('%12s', 4), '%20s'), 
                          c("ID", "OBJ_TYP", "OBJ_TYP_NO", "HYD_TYP", "FILENAME")), 
                          collapse = ' '), obj_prt_path , append = TRUE, 
            sep = "\t", dec = ".", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(paste(sprintf(c(rep('%12s', 4), '%20s'), 
                          c("1", "sdc", as.character(outflow_reach), "tot", "cha_day.out")), 
                  collapse = ' '), obj_prt_path , append = TRUE, 
            sep = "\t", dec = ".", row.names = FALSE, col.names = FALSE, quote = FALSE)

## Updating in the result directories
overwrite_file("object.prt")

file_cio <- readLines(paste0(tmp_setup_path, "/", "file.cio"))
file_cio[2] <- "simulation        time.sim          print.prt         object.prt        object.cnt        null "
writeLines(file_cio, paste0(tmp_setup_path, "/", "file.cio"))

overwrite_file("file.cio")

##------------------------------------------------------------------------------
## 8) Run all climate scenarios in parallel and collect results
##------------------------------------------------------------------------------

## Prepare parallelization
cl <- makeCluster(cores,  outfile="")
registerDoParallel(cl)

## Write files, if missing
txt_info <- foreach (d = m_dir) %dopar% {
  file.copy(paste0(tmp_path, "/", "setup/", setdiff(list.files(path = tmp_setup_path), 
                                                  list.files(path = d))), d)}
## Run parallelization
wd <- getwd()
txt_info <- foreach (d = m_dir) %dopar% {
  # run SWAT for all cal files in parallel
  setwd(paste(wd, d, sep='/'))
  system(swat_exe)
  
  #Collect output files
  cal <- paste(wd, tmp_path, "sim", tail(unlist(strsplit(d,'/')), n = 1), sep = "/")
  dir.create(cal, recursive = TRUE)
  
  files.out.aa <- dir(getwd(), pattern = 'aa')
  files.out.mon <- dir(getwd(), pattern = 'mon')
  files.out.day <- dir(getwd(), pattern = 'day')
  files.help <- c('hru.con', 'hru_agr.txt') # please provide hru_agr.txt file 
  ## (cropland hru names) in txt folder
  file.copy(c(files.out.aa,files.out.mon,files.out.day,files.help), cal, overwrite = T)
  file.remove(c(files.out.aa,files.out.mon,files.out.day,files.help))
  file.remove(list.files(path = paste(wd, d, sep='/'), 
                         pattern = ".*.txt|.*.out|.*.swf|.*.mgts|.*.farm", full.names = TRUE))
}
stopCluster(cl)

##------------------------------------------------------------------------------
## 9) Load Micha's functions for the output extraction and prepare results
##------------------------------------------------------------------------------

source('lib/calc_Indis.R')

## Direstory with the simulation results
path <- paste(tmp_path, "sim", sep = "/")

##------------------------------------------------------------------------------
## 10)  Output analysis from Micha (warranty provided by Micha ;)
##------------------------------------------------------------------------------

### In the following functions to calculate indicators are applied
### Please adjust function parameters (e.g. channel name, see also header 
## information of calc_Indis.R)
### In case an ensemble of calibration files is provided (in folder cal_files), 
## set ensemble=T
### The resulting dataframe will provide you the ensemble mean as well as the 
## ensemble minimum (lower) 
### and maximum (upper) of the respective indicator
### If no cal file ensemble can be provided, set ensemble=F 
### (but then make sure you have a calibration.cal with fitted parameters in 
## the txt folder)

### collect average annual output of water quantity and quality at outlet channel 
## (aggregated comparison)
r_dir <- list.dirs(path, recursive = TRUE)[-1]
rch <- sprintf("cha%03d", outflow_reach)
cha_aa_all <- ind_cha_aa(r_dir, rch) #adjust channel

###  collect indicators related to the daily dynamics Water, N, P, Sed
cha_day_all <- ind_cha_dayII(r_dir, rch, 'all') #adjust channel

###  collect HRU-based indicators related to water quality (average annual losses)
hru_aa_nb_all <- ind_hru_aa_nb(r_dir) #adjust channel

###  collect HRU-based indicators related to  water quantity (average annual values)
hru_aa_wb_all <- ind_hru_aa_wb(r_dir)

###  collect HRU-based indicators related to  water quantity (average annual for 
## specified months)
## please specify start and end months of interest for the soil water analysis
sw_periods <- list(c(5:9), 5, 6, 7, 8, 9) #this is an example for printing sw for 
## the period May to September and also for each single month in that period
hru_mon_all <- ind_hru_mon_wb(r_dir, period = sw_periods, nrows = 54) #might take a while

### collect cropping information for all scenarios - grain units and cultivated 
## hectare average annual
## define 1) path, 2) crop selection, 3) type of output: a) yield, b) ha, 4) 
## specify grain units equivalent for 
#  all of the selected crops (if you just keep the parameter 'grain_units', 
## there is already a parameterisation for 
#   'wwht', 'akgs', 'wbar', 'wira', 'csil', 'wiry', 'sgbt','barl'
# the measure list (measr.list) can be adapted to the measures you want to compare

#adjust
crop_aa_gu <- ind_bsn_aa_crp(r_dir, crop_sel, out_type = "yield", grain_units)
crop_aa_ha <- ind_bsn_aa_crp(r_dir, crop_sel, out_type = "ha", grain_units)

### collect cropping information for all scenarios - Crop specific average annual 
## yield and ha
crop_aa_all <- ind_bsn_aa_crp_ha_Y(r_dir, crop_sel)

##------------------------------------------------------------------------------
## 11)  Aggregate outputs into one dataframe for plotting
##------------------------------------------------------------------------------

### combine data of interest in one dataframe for plotting
df_plot <- cha_aa_all 
for(l in list(cha_day_all, hru_aa_nb_all, hru_aa_wb_all, hru_mon_all, crop_aa_gu, 
              crop_aa_ha, crop_aa_all)) df_plot <- left_join(df_plot, l, by = "scen_name")

df_plot_long <- pivot_longer(df_plot, cols = -scen_name, names_to = "indi", values_to = "value")

## Get long format for plotting
df_plot_long  <- df_plot_long[!grepl("_H$", df_plot_long$scen_name),] %>% 
  mutate(scen_base = gsub("_[[:alpha:]]*$","",scen_name)) %>% 
  left_join(df_plot_long[grepl("_H$", df_plot_long$scen_name),] %>% 
              mutate(scen_base = gsub("_[[:alpha:]]*$","",scen_name)) %>% 
              rename(value_base = value) %>% 
              select(-scen_name), by = c("scen_base", "indi")) %>% 
  mutate(value = round(100*(1 - value/value_base), 3)) %>% 
  select(-ends_with("base")) %>%
  mutate(scen_name = toupper(scen_name)) %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) %>% 
  mutate(scen_name = gsub("_[^_]*_", "_", scen_name)) %>% 
  separate(scen_name, into = c("RCP", "Period"), sep = "_") %>% 
  mutate(Period = case_when(Period == "N" ~ "Near future",
                         Period == "E" ~ "End century")) %>%
  mutate(Period = factor(Period, levels = c("Near future", "End century")))

##------------------------------------------------------------------------------
## 12)  Plotting your results
##------------------------------------------------------------------------------

## Print all available indicators
unique(df_plot_long$indi)

## Plotting selected indicators
throw_box(df_plot_long, c("precip", "snofall", "snomlt", "surq_gen", 
                          "latq", "wateryld", "et", "ecanopy", "eplant", "esoil",
                          "surq_cont", "cn"))

throw_box(df_plot_long, c("Q_mean", "Nload", "Pload", "Sedload", "Q_max", 
                          "Q_p95", "Q_p90", "Q_p50", "Q_p10", "Q_p05", "Q_min", 
                          "Q_maxmin", "Q_p95p05", "Q_p90p10", "Q_low_days", 
                          "Q_high_days"))
