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
  # frm <- SWATfarmR::farmr_project$new(project_name = 'frm', project_path = pth)
  frm <- SWATfarmR::farmr_project$new(project_name = 'frm', project_path = pth, 
                                      project_type = 'environment') ## If you your version farmR >= 4.0.0
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
                  summarize(l_lim = lim(value)[1], u_lim = lim(value)[2]), by = "indi") %>%
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

