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
cores <- detectCores() - 2
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
## 7) Run scenarios with SWATrunR and save results (very long)
##------------------------------------------------------------------------------

res <- run_scenario(project_path = setup_dir,
                    scenario_path = tmp_result_path,
                    output = results_to_save,
                    version = "plus",
                    years_skip = 3,
                    n_thread = cores,
                    keep_folder = TRUE)

## Saving results just in case
saveRDS(res, paste0(tmp_result_path, "/", "scenario_results.rds"))
print("Finished climate runs!!!")

##------------------------------------------------------------------------------
## 7) Converting SWATrunR results into dataframe
##------------------------------------------------------------------------------

##Combine results into one dataframe
df <- map_dfr(names(res$simulation), function(x){
  map_dfr(names(res$simulation[[x]]), function(y){
    data.frame(res$simulation[[x]][[y]], scenario = x, variable = y)
  })
})

## Add missing columns, etc.
mcols <- c("values", "RCP", "RCM", "period")
df[mcols] <- data.frame(df$run_1, str_split_fixed(df$scenario, "_", 3))
df$period <- factor(df$period, levels = c(periods[[1]][[1]], periods[[2]][[1]],
                                          periods[[3]][[1]]))
df <- df[c("date", "variable", mcols)]

##------------------------------------------------------------------------------
## 8) Preparing plots for results
##------------------------------------------------------------------------------

##Function for plotting
plot_fig <- function(df, var_name){
  ggplot(df[df$variable == var_name,], aes(x = date, y = values, colour = RCP))+
    geom_line()+
    # geom_smooth(method=lm, formula = y ~ x)+
    labs(x = "Date", y = paste0("Units of ", toupper(var_name))) +
    theme_bw()+
    facet_grid(RCM ~ period,  scales = "free_x")
}

##Prepare combined plot
l <- map(unique(df$variable), ~plot_fig(df, .x))
names(l) <- unique(df$variable)

## Save it
saveRDS(l, paste0(tmp_result_path, "/", "plot_results.rds"))
print("Happy end :) Finished plot preparation preparation!!!")
