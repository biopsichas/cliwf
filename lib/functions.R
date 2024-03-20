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
  ## If you your version farmR >= 4.0.0
  if(startsWith(as.character(packageVersion("SWATfarmR")), "4.")){
    frm <- SWATfarmR::farmr_project$new(project_name = 'frm', project_path = pth, 
                                        project_type = 'environment') 
    .GlobalEnv$frm <- frm
  } else if(startsWith(as.character(packageVersion("SWATfarmR")), "3.2.0")){
    frm <- SWATfarmR::farmr_project$new(project_name = 'frm', project_path = pth)
  } else {
    stop("SWATfarmR version should be > 4.0.0 or special version 3.2.0 
         from OPTAIN Cloud>WPs&Tasks>WP4>Task4.4>Tools to share>workflow_scripts>SWATfarmR_3.2.0.zip")
  }
  api <- variable_decay(frm$.data$variables$pcp, -5,0.8)
  asgn <- select(frm$.data$meta$hru_var_connect, hru, pcp)
  frm$add_variable(api, "api", asgn)
  frm$read_management(mgt, discard_schedule = TRUE)
  frm$schedule_operations(start_year = start_d, end_year = end_d, replace = 'all')
  frm$write_operations(start_year = start_d, end_year = end_d)
}

##Function to overwrite files in setup
overwrite_file <- function(file_name){
  file.remove(paste0(m_dir, "/", file_name))
  file.copy(paste0(tmp_setup_path, "/", file_name), paste0(m_dir, "/", file_name), overwrite = T)
}

##Function to plot results
throw_box <- function(df, vars, drop_outliers = FALSE, font_size = NULL){
  df <- df[df$indi %in% vars,] %>% 
    mutate(indi = toupper(indi))
  df$indi <- factor(df$indi, levels = toupper(vars))
  
  if(drop_outliers){
    nb_rec <- length(df$indi)
    # Filter out outliers
    df <- df %>%
      left_join(group_by(., indi) %>%
                  dplyr::summarise(l_lim = lim(value)[1], u_lim = lim(value)[2]), by = "indi") %>%
      filter(value >= l_lim & value <= u_lim)
    
    message("Number of outliers removed: ", nb_rec - length(df$indi))
  }
  
  # Plot
  fig <- ggplot(df, aes(x = Period, y = value, fill = RCP))+
    geom_boxplot(outlier.size=1, outlier.color = "grey30")+
    labs(x = "Period", y = "Change compared to baseline period [%]") +
    theme_bw()+
    facet_wrap(~indi,  scales = "free_y")
  
  if(length(font_size)>0){
    fig <- fig + 
      theme(text = element_text(size=font_size))
  }

  return(fig)
}

# Calculate whisker limits
lim <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower_limit <- q1 - 1.5 * iqr
  upper_limit <- q3 + 1.5 * iqr
  return(c(lower_limit, upper_limit))
}

read_model_output <- function(output_path, hru_shp_path, cha_shp_path, res_shp_path){
  # get buildR GIS data 
  hru <- read_sf(hru_shp_path)
  hru$area <- st_area(hru) %>% drop_units
  cha <- read_sf(cha_shp_path)
  res <- read_sf(res_shp_path)
  
  # get SWAT+ results for status quo 
  wb_aa_sq <- read_tbl2('hru_wb_aa.txt', output_path, 3)
  ls_aa_sq <- read_tbl2('hru_ls_aa.txt', output_path, 3)
  pw_aa_sq <- read_tbl2('hru_pw_aa.txt', output_path, 3)
  cha_aa_sd_sq <- read_tbl2('channel_sd_aa.txt', output_path, 3)
  
  # join results with GIS data
  hru_wb_sq <- left_join(hru, wb_aa_sq, 'name')
  hru_ls_sq <- left_join(hru, ls_aa_sq, 'name')
  hru_pw_sq <- left_join(hru, pw_aa_sq, 'name')
  cha_sq <- left_join(cha, cha_aa_sd_sq, 'name')
  
  # Unit conversions etc.
  
  # for hru fluxes provided in mm it might be appropriate to convert to l/s per hru
  # to get rid of extremely high values for small hrus located along flow concentration pathways
  
  for(i in c("surq_gen_lshru", "surq_ls_lshru", "surq_cha_lshru","surq_runon_lshru", "latq_lshru", "perc_lshru", "latq_ls_lshru")){
    hru_wb_sq[[i]] <- hru_wb_sq[[str_replace(i, "_lshru", "")]] * hru_wb_sq[[c("area")]] / (365.25 * 24 * 3600)
  }
  hru_wb_sq$surq_runon_cmahru <- hru_wb_sq$surq_runon * hru_wb_sq$area / 1000 #m³/s per hru
  
  # for hru fluxes provided in tons/ha or kg/ha it might be appropriate to convert to tons or kg per hru
  # to get rid of extremely high values for small hrus located along flow concentration pathways
  
  for(i in c("sedyld_tonshru", "surqno3_kghru", "lat3no3_kghru","sedorgp_kghru", "surqsolp_kghru", "sedminp_kghru")){
    hru_ls_sq[[i]] <- hru_ls_sq[[str_replace(i, c("_kghru|_tonshru"), "")]] * hru_ls_sq[[c("area")]] / 10000
  }
  
  for(i in c("percn_kghru", "pplnt_kghru", "nplt_kghru")){
    hru_pw_sq[[i]] <- hru_pw_sq[[str_replace(i, "_kghru", "")]] * hru_pw_sq[[c("area")]] / 10000
  }
  
  # total amounts of N and P lost from hrus
  hru_ls_sq$ntot <- hru_ls_sq$surqno3_kghru + hru_ls_sq$lat3no3_kghru + hru_pw_sq$percn_kghru
  hru_ls_sq$ptot <- hru_ls_sq$sedminp_kghru + hru_ls_sq$surqsolp_kghru
  hru_ls_sq$ntotha <- hru_ls_sq$surqno3 + hru_ls_sq$lat3no3 + hru_pw_sq$percn
  hru_ls_sq$ptotha <- hru_ls_sq$sedminp + hru_ls_sq$surqsolp
  
  # total N and P in channel network
  cha_sq$n_out <- cha_sq$orgn_out+cha_sq$no2_out+cha_sq$no3_out+cha_sq$nh3_out
  cha_sq$p_out <- cha_sq$solp_out+cha_sq$sedp_out
  
  # concentrations in mg/l
  cha_sq$no3_out_conc <- cha_sq$no3_out*1000/(cha_sq$flo_out*3600*24*365)
  cha_sq$p_out_conc <- (cha_sq$solp_out+cha_sq$sedp_out)*1000/(cha_sq$flo_out*3600*24*365)
  cha_sq$sed_out_conc <- cha_sq$sed_out*1000/(cha_sq$flo_out*3600*24*365)
  
  return(list(hru_wb_sq = hru_wb_sq, hru_ls_sq = hru_ls_sq, 
              hru_pw_sq = hru_pw_sq, cha_sq = cha_sq))
}

read_tbl2 <- function(file, run_path, n_skip) {
  if(length(run_path) == 1){
    if(!file.exists(paste0(run_path, '/', file))){
      stop(paste0('File "', file, '" does not exist in "', run_path, '"'))
    }
  }
  if(length(run_path) > 1){
    i <- FALSE
    for(f in run_path){
      if(file.exists(paste0(run_path[i], '/', file))){
        warning(paste0('File "', file, '" does not exist in "', run_path[i], '"'))
        i <- TRUE
      }
    }
    if(i) stop("Please remove paths without files")
  }
  rf <- function(file, run_path, n_skip){
    file_path <- paste0(run_path, '/', file)
    
    col_names <- read_lines(file = file_path, skip = 1, n_max = 1, lazy = FALSE) %>%
      str_trim(.) %>%
      str_split(., '[:space:]+') %>%
      unlist()
    
    name_duplicate <- table(col_names) %>%
      .[. > 1]
    if(length(name_duplicate) > 0) {
      for (i in 1:length(name_duplicate)) {
        col_names[col_names == names(name_duplicate[i])] <-
          paste0(names(name_duplicate[i]), 1:name_duplicate[i])
      }
    }
    
    fread(file_path, skip = n_skip, header = FALSE) %>%
      set_names(., col_names) %>%
      tibble(.)
  }
  
  if(length(run_path) == 1){
    df <- rf(file, run_path, n_skip)
  } else if (length(run_path) > 1){
    ls <- map(run_path, ~tryCatch(
      expr = {rf(file, .x, n_skip) %>% select(where(is.numeric))
      },
      error = function(e){ 
        cat(paste("Error in", .x))
      }))
    ls_ch <- map(run_path, ~tryCatch(
      expr = {rf(file, .x, n_skip) %>% select(where(is.character))
      }))
    name <- NULL
    i <- 1
    while(is.null(name) & i < length(ls_ch)){
      if(length(ls_ch[[i]][['name']]) > 0) name <- ls_ch[[i]][['name']]
      i <- i + 1
    }
    df <- reduce(ls, `+`)/length(run_path)
    df <-  bind_cols(name = as.data.frame(name), df)
    # message(paste0('Output values of file ', file , ' represent the mean of ', paste0(run_path, collapse = ', '), ' folders.'))
  }
  return(df)
}



