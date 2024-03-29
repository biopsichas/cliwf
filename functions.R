#' Function to copy exe between destinations and run it
#'
#' @param path_from character path of directory from which file should be taken
#' @param path_to character path to which file should be put
#' @param file_name character names of file to copy and run

exe_copy_run <- function(path_from, path_to, file_name){
  # Copy into the destination directory
  file.copy(from = paste0(path_from, "/", file_name), 
            to = paste0(path_to, "/", file_name), overwrite = TRUE)
  
  ##Reset working directory to setup location
  wd_base <- getwd()
  if (str_sub(getwd(), -nchar(path_to), -1) != path_to) setwd(path_to)
  
  ##Write files
  system(file_name)
  
  ##Reset back working directory
  setwd(wd_base)
}

## Function to prepare climate input for scenario
write_cli <- function(rcp, rcm){
  walk2(rcp, rcm, function(x, y){
    cli_lst <- load_swat_weather(paste(cli_dir, x, y, sep = "/"))
    ##Each period
    walk(periods, function(p){
      ##Creating new directory
      n_dir <- paste0(tmp_result_path, "/", x, "_", "rcm", y, "_", p[1])
      dir.create(n_dir, recursive = TRUE)
      ##Coping files to update
      file.copy(paste0(tmp_setup_path,"/", files_to_update), n_dir, overwrite = TRUE)
      ##Running function to prepare climate input for scenario
      prepare_climate(cli_lst, n_dir, p[2], p[3])
    })
  })
}

##Function to run farmR fully
write_mgt <- function(pth, p = periods){
  i <- 1
  while(str_sub(pth,-1) != p[[i]][[1]] & i < length(p)){
    i <- i + 1
  }
  if(str_sub(pth,-1) == p[[i]][[1]]){
    start_d <- as.numeric(substr(p[[i]][[2]], 1, 4))
    end_d <- as.numeric(substr(p[[i]][[3]], 1, 4))
  } else {
    stop("Folder's name last letter should correspond to period name letter.
         Please rename them to fit this!!!")
  }
  frm <- SWATfarmR::farmr_project$new(project_name = 'frm', project_path = pth)
  api <- variable_decay(frm$.data$variables$pcp, -5,0.8)
  asgn <- select(frm$.data$meta$hru_var_connect, hru, pcp)
  frm$add_variable(api, "api", asgn)
  frm$read_management(mgt, discard_schedule = TRUE)
  frm$schedule_operations(start_year = start_d, end_year = end_d, replace = 'all')
  frm$write_operations(start_year = start_d, end_year = end_d)
}
