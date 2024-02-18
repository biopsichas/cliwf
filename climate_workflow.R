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

source('settings.R')
source('functions.R')

##------------------------------------------------------------------------------
## 2) Prepare fresh management files
##------------------------------------------------------------------------------
## Saving reference setup files
tmp_setup_path <- paste0(tmp_path, "/", "setup")
## Saving reference result files
tmp_result_path <- paste0(tmp_path, "/", "results")

## If the directory exists, delete the results directory. (Please be careful!!!)
if (file.exists(tmp_path)) unlink(tmp_path, recursive = TRUE)
dir.create(tmp_setup_path, recursive = TRUE)

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

## Input for the print file
print_prt <- read_lines(paste0(tmp_setup_path,'/print.prt'), lazy = FALSE)
print_prt[3] <- "3           0         0         0         0         1         "
print_prt[41] <- "channel                      d             n             y             y  "
print_prt[42] <- "channel_sd                   d             n             y             y  "
write_lines(print_prt, paste0(tmp_setup_path,'/print.prt'))

## Updating in the result directories
file.remove(paste0(m_dir, "/", "print.prt"))
file.copy(paste0(tmp_setup_path, "/print.prt"), paste0(m_dir, "/", "print.prt"), overwrite = T)

##------------------------------------------------------------------------------
## 8) Run scenarios with save results (very long)
##------------------------------------------------------------------------------

## Prepare parallelization
cl <- makeCluster(cores,  outfile="")
registerDoParallel(cl)

wd <- getwd()
txt_info <- foreach (d = m_dir) %dopar% {
  file.copy(paste0(tmp_path, "/", "setup/", setdiff(list.files(path = tmp_setup_path), 
                                                    list.files(path = d))), d)
  # run SWAT for all cal files in parallel
  setwd(paste(wd, d, sep='/'))
  system(swat_exe)
  
  #file handling
  cal <- paste(wd, tmp_path, "sim", tail(unlist(strsplit(d,'/')), n = 1), sep = "/")
  dir.create(cal, recursive = TRUE)
  
  files.out.aa <- dir(getwd(), pattern = 'aa')
  files.out.mon <- dir(getwd(), pattern = 'mon')
  files.out.day <- dir(getwd(), pattern = 'day')
  files.help <- c('hru.con', 'hru_agr.txt') # please provide hru_agr.txt file (cropland hru names) in txt folder
  file.copy(c(files.out.aa,files.out.mon,files.out.day,files.help), cal, overwrite = T)
}
stopCluster(cl)

##------------------------------------------------------------------------------
## 9) Load Micha's functions for the output extraction and prepare results
##------------------------------------------------------------------------------

source('calc_Indis.R')

foo1(c('dplyr' , 'readr' , 'tidyverse', 'data.table', 'remotes', 'devtools', 
       'xts', 'dygraphs', 'R.utils', 'foreach', 'doParallel', 'data.table', 
       'ggplot2', 'fmsb'))

path <- paste(tmp_path, "sim", sep = "/")

##------------------------------------------------------------------------------
## 10)  Output analysis from Micha (just copied, not tested!!!!)
##------------------------------------------------------------------------------

### In the following functions to calculate indicators are applied
### Please adjust function parameters (e.g. channel name, see also header information of calc_Indis.R)
### In case an ensemble of calibration files is provided (in folder cal_files), set ensemble=T
### The resulting dataframe will provide you the ensemble mean as well as the ensemble minimum (lower) 
### and maximum (upper) of the respective indicator
### If no cal file ensemble can be provided, set ensemble=F 
### (but then make sure you have a calibration.cal with fitted parameters in the txt folder)

### collect average annual output of water quantity and quality at outlet channel (aggregated comparison)
cha_aa_all <- ind_cha_aa(path, 'cha0926', ensemble=T) #adjust channel

###  collect indicators related to the daily dynamics Water, N, P, Sed
cha_day_all <- ind_cha_dayII(path, 'cha0926', 'all', ensemble=T) #adjust channel

###  collect HRU-based indicators related to water quality (average annual losses)
hru_aa_nb_all <- ind_hru_aa_nb(path, a='agr', ensemble=T) #adjust channel

###  collect HRU-based indicators related to  water quantity (average annual values)
hru_aa_wb_all <- ind_hru_aa_wb(path, a='agr', ensemble=T)

###  collect HRU-based indicators related to  water quantity (average annual for specified months)
## please specify start and end months of interest for the soil water analysis
sw_periods <- list(c(5:9), 5, 6, 7, 8, 9) #this is an example for printing sw for the period May to September and also for each single month in that period
hru_mon_all <- ind_hru_mon_wb(path, period = sw_periods, a='agr', ensemble=T) #might take a while

### collect cropping information for all scenarios - grain units and cultivated hectare average annual
# define 1) path, 2) crop selection, 3) type of output: a) yield, b) ha, 4) specify grain units equivalent for 
#  all of the selected crops (if you just keep the parameter 'grain_units', there is already a parameterisation for 
#   'wwht', 'akgs', 'wbar', 'wira', 'csil', 'wiry', 'sgbt','barl'
# the measure list (measr.list) can be adapted to the measures you want to compare

crop_sel <- c('wwht', 'akgs', 'wbar', 'wira', 'csil', 'wiry', 'sgbt','barl') #adjust
crop_aa_gu <- ind_bsn_aa_crp(path, crop_sel, out_type = "yield", grain_units, ensemble=T)
crop_aa_ha <- ind_bsn_aa_crp(path, crop_sel, out_type = "ha", grain_units, ensemble=T)

### collect cropping information for all scenarios - Crop specific average annual yield and ha
crop_aa_all <- ind_bsn_aa_crp_ha_Y(path, crop_sel, ensemble=T)

##### ----------------
# prepare data for plotting
##### ----------------

### combine data of interest in one dataframe for plotting
df_plot <- data.frame(id=1:nrow(cha_aa_all), cha_aa_all[,c(1:5)])
df_plot <- merge(df_plot, cha_day_all[,c(1:16)], by = "scen_name")
df_plot <- merge(df_plot, hru_aa_nb_all[,c(1:6)], by = "scen_name")
df_plot <- merge(df_plot, hru_aa_wb_all[,c(1:3)], by = "scen_name")
df_plot <- merge(df_plot, hru_mon_all[,c(1:(grep('lower', names(hru_mon_all))[1]-1))], by = "scen_name")
df_plot <- merge(df_plot, crop_aa_gu[,c(1:2)], by = "scen_name")
df_plot <- merge(df_plot, crop_aa_ha[,c(1:2)], by = "scen_name")
df_plot <- merge(df_plot, crop_aa_all[,c(1:(grep('lower', names(crop_aa_all))[1]-1))], by = "scen_name")
df_plot <- df_plot[order(df_plot$id), ]

## calculate the scenario effects (absolute effects), by subtracting the status-quo outputs from the scenario outputs
df_plot[-1, -1:-2] <- df_plot[-1,-1:-2]-df_plot[rep(1, nrow(df_plot)-1), -1:-2]

## calculate the percentage change 
df_plot[-1, -1:-2] <- round((df_plot[-1,-1:-2]/df_plot[rep(1, nrow(df_plot)-1), -1:-2])*100, 3)

## prepare a long format dataframe, needed for specific plots (leaving out status quo)
df_plot_long <- melt(setDT(df_plot[-1, -2]), id.vars = c("scen_name"), variable.name = "indi")

## define specific order for plotting
## ! adapt this vector if you have other indicators (e.g. different sw_periods and crops) !
##   or want to modify the order in the plot!
colnames(df_plot)
df_plot_order <- c("Q_mean", "Q_max", "Q_p95", "Q_p90", "Q_p50", "Q_p10", "Q_p05", "Q_min", 
                   "Q_maxmin", "Q_p95p05", "Q_p90p10", "Q_low_days", "Q_high_days",
                   "perc", "sw", "sw_5_6_7_8_9", "sw_5", "sw_6", "sw_7", "sw_8", "sw_9", 
                   "Nload", "Nconc_days", "N_loss", "N_loss_ratio",
                   "Pload", "Pconc_days", "P_loss", "P_loss_ratio",
                   "Sedload", "Sedconc_days", "Sed_loss", 
                   "grain_units_aa", "crops_ha_aa", 
                   "wwht_ha", "akgs_ha", "wbar_ha", "wira_ha", "csil_ha", "wiry_ha", "sgbt_ha", 
                   "barl_ha", 
                   "wwht_yld_t_ha", "akgs_yld_t_ha", "wbar_yld_t_ha", "wira_yld_t_ha", 
                   "csil_yld_t_ha", "wiry_yld_t_ha", "sgbt_yld_t_ha", "barl_yld_t_ha")
# adapt the factors accordingly
df_plot_long$indi <- factor(df_plot_long$indi, levels = rev(df_plot_order))


##### ----------------
# Plotting barplots, spider web charts, profile plots
##### ----------------

### First rough plot:
## Plot profile-chart with percentage change
ggplot(data=df_plot_long, 
       aes(x=indi, y=value, group=scen_name, color=scen_name)) +
  geom_line() + 
  geom_point(show.legend = FALSE) +
  coord_flip() + 
  xlab("Group") +
  ylab("Value")

### Plot coloured vertical barplot
#Breaks for background rectangles
# ! @all: please modify the following if you changed the number of indicators in the plot 
# (check 'df_plot_order')
rects <- data.frame(xstart = c(0.5, 18.5, 21.5, 25.5, 29.5), 
                    xend = c(18.5, 21.5, 25.5, 29.5, 50.5), 
                    col = c("yellow", "tan4", "gray90", "seagreen1", "cadetblue1"))

## version 1 - separate bar for each scenario
ggplot() +
  geom_bar(data = df_plot_long, 
           aes(x = indi, y = value, fill = scen_name), 
           stat = "identity", 
           position = "dodge") + 
  labs(x="indicators", y="relative effects compared to status quo [%]") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_brewer(palette = "Dark2") + 
  scale_color_brewer(palette = "Dark2") + 
  coord_flip() + 
  geom_rect(data = rects, 
            aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf),
            fill = rects$col,
            alpha = 0.2) +
  geom_bar(data = df_plot_long, 
           aes(x = indi, y = value, fill = scen_name), 
           stat = "identity", 
           position = "dodge")


## version 2 - stacked bar + separate point for "all" scenario
ggplot() +
  geom_bar(data = subset(df_plot_long, df_plot_long$scen_name!="all"), 
           aes(x = indi, y = value, fill = scen_name), 
           stat = "identity", 
           position = "stack") + 
  labs(x="indicators", y="relative effects compared to status quo [%]") + 
  scale_fill_brewer(palette = "Dark2") + 
  scale_color_brewer(palette = "Dark2") + 
  coord_flip() + 
  geom_rect(data = rects, 
            aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf),
            fill = rects$col,
            alpha = 0.2) +
  geom_bar(data = subset(df_plot_long, df_plot_long$scen_name!="all"), 
           aes(x = indi, y = value, fill = scen_name), 
           stat = "identity", 
           position = "stack") +
  geom_point(data = subset(df_plot_long, df_plot_long$scen_name=="all"), 
             aes(x = indi, y = value, group = 1, color = 'all_scen'),
             size = 3,
             shape= 1,
             stroke = 2) +
  geom_hline(yintercept = 0, col = "black") + 
  scale_color_manual(labels = c('all_scen'), 
                     values = c('black')) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.title = element_blank())


### plot radar chart (spider web)

# function for adding MAx min values to the first rows
set_max_min <- function(df){
  df_max <- max(df, na.rm = TRUE)
  df_min <- min(df, na.rm = TRUE)
  df <- rbind(max = df_max, min = df_min, df)
  return(df)
}

#Prepare data (e.g. add max + min values; add row-names)
df_plot_spider <- df_plot[, -2]
row.names(df_plot_spider) <- df_plot_spider[,1]
df_plot_spider <- df_plot_spider[, -1]
# set status-quo to 0
df_plot_spider["status_quo",] <- 0
df_plot_spider <- set_max_min(df_plot_spider)
# reorder columns as defined above
df_plot_spider <- df_plot_spider[df_plot_order]


## finetuning radar chart
# select indicators you like to plot
plot_sel <- df_plot_spider[,c(1:34)]
# select colors for the scenarios
col_sel <- c("blue", "green", "brown", "skyblue", "forestgreen", "yellow", "black")
# re-calculate max_min
plot_sel <- plot_sel[c(-1,-2),]
plot_sel <- set_max_min(plot_sel)
# plot
op <- par(mar = c(1, 2, 2, 1))
radarchart(df = plot_sel,
           maxmin = T,
           axistype = 1,
           seg = 10,
           caxislabels = c(round(plot_sel[2,1],1),
                           rep(NA,9),      # needs to be the number of segments-1 )
                           round(plot_sel[1,1],1)),
           vlcex = 0.6, vlabels = colnames(plot_sel),
           title = "measure effects [%]" ,
           pcol = col_sel, #pfcol = scales::alpha(col_sel, 0.2),
           plwd = 2, plty = 1,
           cglcol = "grey", cglty = 1, cglwd = 0.8, axislabcol = "red",
)
# Add an horizontal legend
legend(
  x = "right", legend = rownames(df_plot_spider[-c(1,2),]), horiz = F,
  bty = "n", pch = 20 , col = col_sel,
  text.col = "black", cex = 1, pt.cex = 1.5
)
#par(op_bk)


## Create a separarte Radar Chart for each "Indicator category" (Water, N + P + Sed, Yield))
# select colors for the scenarios
col_sel <- c("blue", "green", "brown", "skyblue", "forestgreen", "yellow", "black")
# define titles and columns of the indicators you want to consider)
title_sel <- c("Water", "N_P_Sed", "Crops")
title_pos <- list(c(1:21), c(22:32), c(33:50))
# Reduce plot margin and split the screen in 3 parts
op <- par(mar = c(1, 1, 1, 1))
par(mfrow = c(1,3))
# Create the radar chart
for(i in 1:3){
  plot_sel <- df_plot_spider[,title_pos[[i]]]
  # re-calculate max_min
  plot_sel <- plot_sel[c(-1,-2),]
  plot_sel <- set_max_min(plot_sel)
  # plot
  radarchart(df = plot_sel,
             maxmin = T,
             axistype = 1,
             seg = 10,
             caxislabels = c(round(plot_sel[2,1],1),
                             rep(NA,9),      # needs to be the number of segments-1 )
                             round(plot_sel[1,1],1)),
             vlcex = 0.8, vlabels = colnames(plot_sel),
             title = paste0("measure effects [%] - ", title_sel[i]),
             pcol = col_sel, #pfcol = scales::alpha(col_sel, 0.2),
             plwd = 2, plty = 1,
             cglcol = "grey", cglty = 1, cglwd = 0.8, axislabcol = "red",
  )
}
#par(op_bk)

##TODO

# - find better solution for the axis lables
# - maybe fill/color the webs





