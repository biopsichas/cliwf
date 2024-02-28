#############################################################################################
##
##~~~~~~~~~~~~~~~ Functions calculating performance indicators ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Indicator     | Description                                                   | Function(parameter)[output]                    | Files (to be defined in print.prt)
##~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Q_mean        | mean discharge [m?/s]                                         | ind_cha_aa(path, channel)[1]                   | channel_sd_aa.txt
## Nload         | total N load [kg/yr]                                          | ind_cha_aa(path, channel)[2]                   | channel_sd_aa.txt
## Pload         | total P load [kg/yr]                                          | ind_cha_aa(path, channel)[3]                   | channel_sd_aa.txt
## Sedload       | total sediment load [tons/yr]                                 | ind_cha_aa(path, channel)[4]                   | channel_sd_aa.txt
##---------------|---------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------
## Q_max         | maximum daily discharge [m?/s]                                | ind_cha_day(path, channel, 'Q_max')[1]         | channel_sd_day.txt
## Q_p95         | 95 percentile daily discharge [m?/s]                          | ind_cha_day(path, channel, 'Q_p95')[2]         | channel_sd_day.txt
## Q_p90         | 90 percentile daily discharge [m?/s]                          | ind_cha_day(path, channel, 'Q_p90')[3]         | channel_sd_day.txt
## Q_p50         | 50 percentile daily discharge [m?/s]                          | ind_cha_day(path, channel, 'Q_p50')[4]         | channel_sd_day.txt
## Q_p10         | 10 percentile daily discharge [m?/s]                          | ind_cha_day(path, channel, 'Q_p10')[5]         | channel_sd_day.txt
## Q_p05         | 5 percentile daily discharge [m?/s]                           | ind_cha_day(path, channel, 'Q_p05')[6]         | channel_sd_day.txt
## Q_min         | minimum daily discharge [m?/s]                                | ind_cha_day(path, channel, 'Q_min')[7]         | channel_sd_day.txt
## Q_maxmin      | max/min daily discharge ratio [m?/s]                          | ind_cha_day(path, channel, 'Q_maxmin')[8]      | channel_sd_day.txt
## Q_p95p05      | 95 percentile/5 percentile daily discharge ratio [m?/s]       | ind_cha_day(path, channel, 'Q_p95p05')[9]      | channel_sd_day.txt
## Q_p90p10      | 90 percentile/10 percentile daily discharge ratio [m?/s]      | ind_cha_day(path, channel, 'Q_p90p10')[10]     | channel_sd_day.txt
## Q_low_days    | frequency daily discharge is below low flow threshold         | ind_cha_day(path, channel, 'Q_low_days', threshold_lowQ)[11]       | channel_sd_day.txt
## Q_high_days   | frequency daily discharge is below high flow threshold        | ind_cha_day(path, channel, 'Q_high_days', threshold_highQ)[12]     | channel_sd_day.txt
## Nconc_days    | frequency total N concentrations is below threshold           | ind_cha_day(path, channel, 'Nconc_days', threshold_N)[13]          | channel_sd_day.txt
## Pconc_days    | frequency total P concentrations is below threshold           | ind_cha_day(path, channel, 'Pconc_days', threshold_P)[14]          | channel_sd_day.txt
## Sedconc_days  | frequency total sediment concentrations is below threshold    | ind_cha_day(path, channel, 'Sedconc_days', threshold_Sed)[15]      | channel_sd_day.txt
##---------------|---------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------
## N_loss        | average annual Nitrogen loss from land objects [kg N/ha,yr]   | ind_hru_aa_nb(path, a)[1] # a=basin by default, a='agr' for cropland only (hru_agr.txt must be provided!) | hru_nb_aa.txt, hru_ls_aa.txt, hru_pw_aa.txt
## P_loss        | average annual Phosphorus loss from land objects [kg P/ha,yr] | ind_hru_aa_nb(path, a)[3] # a=basin by default, a='agr' for cropland only (hru_agr.txt must be provided!) | hru_nb_aa.txt, hru_ls_aa.txt, hru_pw_aa.txt
## Sed_loss      | average annual Sediment loss from land objects [tons/ha,yr]   | ind_hru_aa_nb(path, a)[5] # a=basin by default, a='agr' for cropland only (hru_agr.txt must be provided!) | hru_nb_aa.txt, hru_ls_aa.txt, hru_pw_aa.txt
## N_loss_ratio  | average annual Nitrogen loss/input ratio []                   | ind_hru_aa_nb(path, a)[2] # a=basin by default, a='agr' for cropland only (hru_agr.txt must be provided!) | hru_nb_aa.txt, hru_ls_aa.txt, hru_pw_aa.txt
## P_loss_ratio  | average annual Phosphorus loss/input ratio []                 | ind_hru_aa_nb(path, a)[4] # a=basin by default, a='agr' for cropland only (hru_agr.txt must be provided!) | hru_nb_aa.txt, hru_ls_aa.txt, hru_pw_aa.txt
##---------------|---------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------
## sw            | average annual soil moisture [mm/yr]                          | ind_hru_aa_wb(path, a)[1] # a=basin by default, a='agr' for cropland only (hru_agr.txt must be provided!) | hru_wb_aa.txt
## perc          | average annual percolation [mm/yr]                            | ind_hru_aa_wb(path, a)[2] # a=basin by default, a='agr' for cropland only (hru_agr.txt must be provided!) | hru_wb_aa.txt
##---------------|---------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------
## sw            | average soil moisture in period of interest [mm/yr]           | ind_hru_mon_wb(path, period, a)[1] # e.g. period=c(5:9), a=basin by default, a='agr' for cropland only (hru_agr.txt must be provided!) | hru_wb_mon.txt
##---------------|---------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------
## grain_units_aa| average annual sum of grain units in whole basin [gu]         | ind_bsn_aa_crp(path, grain_units, out_type='yield', crops_sel)[1]        | basin_crop_yld_aa.txt
## crops_ha_aa   | average annual area of cropland in whole basin [gu]           | ind_bsn_aa_crp(path, grain_units, out_type='area', crops_sel)[1]         | basin_crop_yld_aa.txt
## crop_yld_t_ha | average annual crop-specific yields [tons drymass/ha]         | ind_bsn_aa_crp_ha_Y(path, crops_sel)                                     | basin_crop_yld_aa.txt
## crop_ha       | average annual crop-specific area in whole basin [ha]         | ind_bsn_aa_crp_ha_Y(path, crops_sel)                                     | basin_crop_yld_aa.txt

foo1 <- function(x){
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      #  Load package after installing
      require( i , character.only = TRUE )
    }
  }
}

foo2 <- function(x){
  #  require returns TRUE invisibly if it was able to load package
  if( ! require( x , character.only = TRUE ) ){
    #  If package was not able to be loaded then re-install
    remotes::install_github("chrisschuerz/SWATfarmR")
    remotes::install_git('https://git.ufz.de/schuerz/SWATmeasR')
    #  Load package after installing
    require( x , character.only = TRUE )
  }
}

foo3 <- function(x){
  #  require returns TRUE invisibly if it was able to load package
  if( ! require( x , character.only = TRUE ) ){
    #  If package was not able to be loaded then re-install
    devtools::install_github("hzambran/hydroGOF")
    #  Load package after installing
    require( x , character.only = TRUE )
  }
}

is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

set2wd <- function(x) setwd(x)

# Read table function
read_tbl <- function(file, run_path, n_skip) {
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

# Indicators based on annual average channel output
ind_cha_aa <- function(path, channel, ensemble=F){
  
  if(ensemble==F){
    df_out <- data.frame(scen_name=sapply(strsplit(path, split ="/"),tail,1), 
                         Q_mean=NA, 
                         Nload=NA, 
                         Pload=NA, 
                         Sedload=NA)
    for (i in 1:length(path)) {
      # Read file
      channel_sd_aa <- read_tbl('channel_sd_aa.txt', path[i], 3)
      
      # Specify the columns you want to keep
      columns_to_keep <- c("jday", "mon", "day", "yr", "unit", "gis_id", "name", 
                           "flo_out", "sed_out", "orgn_out", "sedp_out", "no3_out", 
                           "solp_out", "nh3_out", "no2_out")
      
      # Create a new data frame with only the selected columns
      df_selected <- channel_sd_aa[, columns_to_keep]
      
      df_selected <- df_selected %>%
        mutate(total_N = no3_out + orgn_out + nh3_out + no2_out) %>% 
        mutate(total_P = solp_out + sedp_out)
      
      Q_mean <- round(df_selected$flo_out[which(df_selected$name==channel)],3)
      Nload <- round(df_selected$total_N[which(df_selected$name==channel)],3)
      Pload <- round(df_selected$total_P[which(df_selected$name==channel)],3)
      Sedload <- round(df_selected$sed_out[which(df_selected$name==channel)],3)
      
      df_out[i,2:5] <- c(Q_mean, Nload, Pload, Sedload)
    }
  }
  
  if(ensemble==T){
    df_out <- data.frame(scen_name=sapply(strsplit(path, split ="/"),tail,1),
                         Q_mean=NA, 
                         Nload=NA, 
                         Pload=NA, 
                         Sedload=NA,
                         Q_mean_lower=NA, 
                         Nload_lower=NA, 
                         Pload_lower=NA, 
                         Sedload_lower=NA,
                         Q_mean_upper=NA, 
                         Nload_upper=NA, 
                         Pload_upper=NA, 
                         Sedload_upper=NA)
    for (i in 1:length(path)) {
      # Read file
      path_ext <- dir(path[i], full.names = T)
      df_out2 <- data.frame(Q_mean=NA, 
                            Nload=NA, 
                            Pload=NA, 
                            Sedload=NA)
      for (k in 1:length(path_ext)){
        channel_sd_aa <- read_tbl('channel_sd_aa.txt', path_ext[k], 3)
        
        # Specify the columns you want to keep
        columns_to_keep <- c("jday", "mon", "day", "yr", "unit", "gis_id", "name", 
                             "flo_out", "sed_out", "orgn_out", "sedp_out", "no3_out", 
                             "solp_out", "nh3_out", "no2_out")
        
        # Create a new data frame with only the selected columns
        df_selected <- channel_sd_aa[, columns_to_keep]
        
        df_selected <- df_selected %>%
          mutate(total_N = no3_out + orgn_out + nh3_out + no2_out) %>% 
          mutate(total_P = solp_out + sedp_out)
        
        Q_mean <- round(df_selected$flo_out[which(df_selected$name==channel)],3)
        Nload <- round(df_selected$total_N[which(df_selected$name==channel)],3)
        Pload <- round(df_selected$total_P[which(df_selected$name==channel)],3)
        Sedload <- round(df_selected$sed_out[which(df_selected$name==channel)],3)
        
        df_out2[k,1:4] <- c(Q_mean, Nload, Pload, Sedload)
      }
      df_out[i,2:5] <- colMeans(df_out2)  
      df_out[i,6:9] <- lapply(df_out2, FUN='min')
      df_out[i,10:13] <- lapply(df_out2, FUN='max')
    }
    
  }
  
  return(df_out)
}

# Indicators based on daily channel output
# Option I: channel_sd_day (all channels will be printed, => try to avoid)
ind_cha_dayI <- function(path,
                         channel,
                         ind, 
                         threshold_lowQ=0.063,
                         threshold_highQ=2.01,
                         threshold_N=2.3, 
                         threshold_P=0.082, 
                         threshold_Sed=10){

  if(ensemble==F){
    df_out <- data.frame(scen_name=sapply(strsplit(path, split ="/"),tail,1), 
                         Q_max=NA, 
                         Q_p95=NA, 
                         Q_p90=NA,
                         Q_p50=NA, 
                         Q_p10=NA, 
                         Q_p05=NA, 
                         Q_min=NA,  
                         Q_maxmin=NA, 
                         Q_p95p05=NA, 
                         Q_p90p10=NA, 
                         Q_low_days=NA,
                         Q_high_days=NA,
                         Nconc_days=NA, 
                         Pconc_days=NA,
                         Sedconc_days=NA)
    
    for (i in 1:length(path)) {
      # Read file
      channel_sd_day <- read_tbl('channel_sd_day.txt', path[i], 3)
      
      if('Q_maxmin' %in% ind | 'Q_max' %in% ind  | 'Q_min' %in% ind | 'all' %in% ind){
        # Specify the columns you want to keep
        columns_to_keep <- c("jday", "mon", "day", "yr", "unit", "gis_id", "name", 
                             "flo_out")
        
        # Create a new data frame with only the selected columns
        df_selected <- channel_sd_day[, columns_to_keep]
        
        # Handle zero flow (defining 0.00000001 cm?/s)
        df_selected$flo_out[which(df_selected$flo_out==0)] <- 1e-40
        
        # Group data by channel (replace "name" with the actual column name for channel)
        max_min_ratio <- df_selected %>%
          group_by(name) %>% 
          summarise(
            max_discharge = max(flo_out, na.rm = TRUE),
            min_discharge = min(flo_out, na.rm = TRUE),
            extreme_streamflow_ratio = max_discharge / min_discharge
          )
        
        df_out[i,2] <- round(max_min_ratio$max_discharge[which(max_min_ratio$name==channel)],3)
        df_out[i,8] <- round(max_min_ratio$min_discharge[which(max_min_ratio$name==channel)],3)
        df_out[i,9] <- round(max_min_ratio$extreme_streamflow_ratio[which(max_min_ratio$name==channel)],3)
      }
      if('Q_p50' %in% ind | 'all' %in% ind){
        # Specify the columns you want to keep
        columns_to_keep <- c("jday", "mon", "day", "yr", "unit", "gis_id", "name", 
                             "flo_out")
        
        # Create a new data frame with only the selected columns
        df_selected <- channel_sd_day[, columns_to_keep]
        
        # Handle zero flow (defining 0.00000001 cm?/s)
        df_selected$flo_out[which(df_selected$flo_out==0)] <- 1e-40
        
        # Group data by channel (replace "name" with the actual column name for channel)
        p50 <- df_selected %>%
          group_by(name) %>% 
          summarise(
            p50_discharge = quantile(flo_out, probs = 0.50, na.rm = TRUE)
          )
        
        df_out[i,5] <- round(as.numeric(p50$p50_discharge[which(p50$name==channel)]),3)
      }
      if('Q_p95p05' %in% ind | 'Q_p95' %in% ind | 'Q_p05' %in% ind | 'all' %in% ind){
        # Specify the columns you want to keep
        columns_to_keep <- c("jday", "mon", "day", "yr", "unit", "gis_id", "name", 
                             "flo_out")
        
        # Create a new data frame with only the selected columns
        df_selected <- channel_sd_day[, columns_to_keep]
        
        # Handle zero flow (defining 0.00000001 cm?/s)
        df_selected$flo_out[which(df_selected$flo_out==0)] <- 1e-40
        # Group data by channel (replace "name" with the actual column name for channel)
        Q_p95p05 <- df_selected %>%
          group_by(name) %>%
          summarise(
            p05_discharge = quantile(flo_out, probs = 0.05, na.rm = TRUE),
            p95_discharge = quantile(flo_out, probs = 0.95, na.rm = TRUE),
            extreme_streamflow_ratio = p95_discharge / p05_discharge
          )
        df_out[i,3] <- round(as.numeric(Q_p95p05$p95_discharge[which(Q_p95p05$name==channel)]),3)
        df_out[i,7] <- round(as.numeric(Q_p95p05$p05_discharge[which(Q_p95p05$name==channel)]),3)
        df_out[i,10] <- round(as.numeric(Q_p95p05$extreme_streamflow_ratio[which(Q_p95p05$name==channel)]),3)
      }
      if('Q_p90p10' %in% ind | 'Q_p90' %in% ind | 'Q_p10' %in% ind | 'all' %in% ind){
        # Specify the columns you want to keep
        columns_to_keep <- c("jday", "mon", "day", "yr", "unit", "gis_id", "name", 
                             "flo_out")
        
        # Create a new data frame with only the selected columns
        df_selected <- channel_sd_day[, columns_to_keep]
        
        # Handle zero flow (defining 0.00000001 cm?/s)
        df_selected$flo_out[which(df_selected$flo_out==0)] <- 1e-40
        # Group data by channel (replace "name" with the actual column name for channel)
        Q_p90p10 <- df_selected %>%
          group_by(name) %>%
          summarise(
            p10_discharge = quantile(flo_out, probs = 0.10, na.rm = TRUE),
            p90_discharge = quantile(flo_out, probs = 0.90, na.rm = TRUE),
            extreme_streamflow_ratio = p90_discharge / p10_discharge
          )
        df_out[i,4] <- round(as.numeric(Q_p90p10$p90_discharge[which(Q_p90p10$name==channel)]),3)
        df_out[i,6] <- round(as.numeric(Q_p90p10$p10_discharge[which(Q_p90p10$name==channel)]),3)
        df_out[i,11] <- round(as.numeric(Q_p90p10$extreme_streamflow_ratio[which(Q_p90p10$name==channel)]),3)
      }
      if('Q_low_days' %in% ind | 'Q_high_days' %in% ind | 'Nconc_days' %in% ind | 'Pconc_days' %in% ind | 'Sedconc_days' %in% ind | 'all' %in% ind){
        # Specify the columns you want to keep
        columns_to_keep <- c("jday", "mon", "day", "yr", "unit", "gis_id", "name", 
                             "flo_out", "sed_out", "orgn_out", "sedp_out", "no3_out", 
                             "solp_out", "nh3_out", "no2_out")
        
        # Create a new data frame with only the selected columns
        df_selected <- channel_sd_day[, columns_to_keep]
        
        # Calculate the sum of no3, orgn, nh3, and no2
        df_selected <- df_selected %>%
          mutate(total_N = no3_out + orgn_out + nh3_out + no2_out) %>% 
          mutate(N_conc_mgl = ifelse(flo_out == 0, 0, (total_N * 1000) / (flo_out * 86400))) %>% 
          mutate(total_P = solp_out + sedp_out) %>% 
          mutate(P_conc_mgl = ifelse(flo_out == 0, 0, (total_P * 1000) / (flo_out * 86400))) %>% 
          mutate(sed_conc_mgl = ifelse(flo_out == 0, 0, (sed_out * 1e6) / (flo_out * 86400)))
        
        # Calculate the frequency of exceeding the thresholds for each "unit" (channel)
        frequency_summary_mean <- df_selected %>%
          group_by(name) %>%
          summarize(
            freq_below_threshold_lowQ = mean(flo_out < threshold_lowQ, na.rm = TRUE),
            freq_below_threshold_highQ = mean(flo_out < threshold_highQ, na.rm = TRUE),
            freq_below_threshold_N = mean(N_conc_mgl < threshold_N, na.rm = TRUE),
            freq_below_threshold_P = mean(P_conc_mgl < threshold_P, na.rm = TRUE),
            freq_below_threshold_Sed = mean(sed_conc_mgl < threshold_Sed, na.rm = TRUE)
          )
        df_out[i,14] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_N[which(frequency_summary_mean$name==channel)]),3)
        df_out[i,15] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_P[which(frequency_summary_mean$name==channel)]),3)
        df_out[i,16] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_Sed[which(frequency_summary_mean$name==channel)]),3)
        df_out[i,12] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_lowQ[which(frequency_summary_mean$name==channel)]),3)
        df_out[i,13] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_highQ[which(frequency_summary_mean$name==channel)]),3)
      }
  }
}
  
  if(ensemble==T){
    df_out <- data.frame(scen_name=sapply(strsplit(path, split ="/"),tail,1), 
                         Q_max=NA, 
                         Q_p95=NA, 
                         Q_p90=NA,
                         Q_p50=NA, 
                         Q_p10=NA, 
                         Q_p05=NA, 
                         Q_min=NA,  
                         Q_maxmin=NA, 
                         Q_p95p05=NA, 
                         Q_p90p10=NA, 
                         Q_low_days=NA,
                         Q_high_days=NA,
                         Nconc_days=NA, 
                         Pconc_days=NA,
                         Sedconc_days=NA,
                         Q_max_lower=NA,
                         Q_p95_lower=NA,
                         Q_p90_lower=NA,
                         Q_p50_lower=NA,
                         Q_p10_lower=NA,
                         Q_p05_lower=NA,
                         Q_min_lower=NA,
                         Q_maxmin_lower=NA,
                         Q_p95p05_lower=NA,
                         Q_p90p10_lower=NA,
                         Q_low_days_lower=NA,
                         Q_high_days_lower=NA,
                         Nconc_days_lower=NA,
                         Pconc_days_lower=NA,
                         Sedconc_days_lower=NA,
                         Q_max_upper=NA,
                         Q_p95_upper=NA,
                         Q_p90_upper=NA,
                         Q_p50_upper=NA,
                         Q_p10_upper=NA,
                         Q_p05_upper=NA,
                         Q_min_upper=NA,
                         Q_maxmin_upper=NA,
                         Q_p95p05_upper=NA,
                         Q_p90p10_upper=NA,
                         Q_low_days_upper=NA,
                         Q_high_days_upper=NA,
                         Nconc_days_upper=NA,
                         Pconc_days_upper=NA,
                         Sedconc_days_upper=NA)
    
    for (i in 1:length(path)) {
      # Read file
      path_ext <- dir(path[i], full.names = T)
      df_out2 <- data.frame(Q_max=NA,
                            Q_p95=NA,
                            Q_p90=NA,
                            Q_p50=NA,
                            Q_p10=NA,
                            Q_p05=NA,
                            Q_min=NA,
                            Q_maxmin=NA,
                            Q_p95p05=NA,
                            Q_p90p10=NA,
                            Q_low_days=NA,
                            Q_high_days=NA,
                            Nconc_days=NA,
                            Pconc_days=NA,
                            Sedconc_days=NA)
      for (k in 1:length(path_ext)){
        channel_sd_day <- read.table('channel_sd_day.txt', path_ext[k], 3)
        
        if('Q_maxmin' %in% ind | 'Q_max' %in% ind  | 'Q_min' %in% ind | 'all' %in% ind){
          
          # convert m?/day to m?/s
          channel_sd_day$flo <- channel_sd_day$flo/86400
          
          # Handle zero flow (defining 0.00000001 cm?/s)
          channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
          
          # Group data by channel (replace "name" with the actual column name for channel)
          max_min_ratio <- channel_sd_day %>%
            group_by(name) %>% 
            summarise(
              max_discharge = max(flo, na.rm = TRUE),
              min_discharge = min(flo, na.rm = TRUE),
              extreme_streamflow_ratio = max_discharge / min_discharge
            )
          
          df_out2[k,1] <- round(max_min_ratio$max_discharge[which(max_min_ratio$name==channel)],3)
          df_out2[k,7] <- round(max_min_ratio$min_discharge[which(max_min_ratio$name==channel)],3)
          df_out2[k,8] <- round(max_min_ratio$extreme_streamflow_ratio[which(max_min_ratio$name==channel)],3)
        }
        if('Q_p50' %in% ind | 'all' %in% ind){
          
          # convert m?/day to m?/s
          if('all' %in% ind == F) channel_sd_day$flo <- channel_sd_day$flo/86400
          
          # Handle zero flow (defining 0.00000001 cm?/s)
          channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
          
          # Group data by channel (replace "name" with the actual column name for channel)
          p50 <- channel_sd_day %>%
            group_by(name) %>% 
            summarise(
              p50_discharge = quantile(flo, probs = 0.50, na.rm = TRUE)
            )
          
          df_out2[k,4] <- round(as.numeric(p50$p50_discharge[which(p50$name==channel)]),3)
        }
        if('Q_p95p05' %in% ind | 'Q_p95' %in% ind | 'Q_p05' %in% ind | 'all' %in% ind){
          
          # convert m?/day to m?/s
          if('all' %in% ind == F) channel_sd_day$flo <- channel_sd_day$flo/86400
          
          # Handle zero flow (defining 0.00000001 cm?/s)
          channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
          
          # Group data by channel (replace "name" with the actual column name for channel)
          Q_p95p05 <- channel_sd_day %>%
            group_by(name) %>%
            summarise(
              p05_discharge = quantile(flo, probs = 0.05, na.rm = TRUE),
              p95_discharge = quantile(flo, probs = 0.95, na.rm = TRUE),
              extreme_streamflow_ratio = p95_discharge / p05_discharge
            )
          df_out2[k,2] <- round(as.numeric(Q_p95p05$p95_discharge[which(Q_p95p05$name==channel)]),3)
          df_out2[k,6] <- round(as.numeric(Q_p95p05$p05_discharge[which(Q_p95p05$name==channel)]),3)
          df_out2[k,9] <- round(as.numeric(Q_p95p05$extreme_streamflow_ratio[which(Q_p95p05$name==channel)]),3)
        }
        if('Q_p90p10' %in% ind | 'Q_p90' %in% ind | 'Q_p10' %in% ind | 'all' %in% ind){
          
          # convert m?/day to m?/s
          if('all' %in% ind == F) channel_sd_day$flo <- channel_sd_day$flo/86400
          
          # Handle zero flow (defining 0.00000001 cm?/s)
          channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
          
          # Group data by channel (replace "name" with the actual column name for channel)
          Q_p90p10 <- channel_sd_day %>%
            group_by(name) %>%
            summarise(
              p10_discharge = quantile(flo, probs = 0.10, na.rm = TRUE),
              p90_discharge = quantile(flo, probs = 0.90, na.rm = TRUE),
              extreme_streamflow_ratio = p90_discharge / p10_discharge
            )
          df_out2[k,3] <- round(as.numeric(Q_p90p10$p90_discharge[which(Q_p90p10$name==channel)]),3)
          df_out2[k,5] <- round(as.numeric(Q_p90p10$p10_discharge[which(Q_p90p10$name==channel)]),3)
          df_out2[k,10] <- round(as.numeric(Q_p90p10$extreme_streamflow_ratio[which(Q_p90p10$name==channel)]),3)
        }
        if('Q_low_days' %in% ind | 'Q_high_days' %in% ind | 'Nconc_days' %in% ind | 'Pconc_days' %in% ind | 'Sedconc_days' %in% ind | 'all' %in% ind){
          
          # convert m?/day to m?/s
          if('all' %in% ind == F) channel_sd_day$flo <- channel_sd_day$flo/86400
          
          # Handle zero flow (defining 0.00000001 cm?/s)
          channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
          
          # Calculate the sum of no3, orgn, nh3, and no2
          channel_sd_day <- channel_sd_day %>%
            mutate(total_N = no3 + orgn + nh3 + no2) %>% 
            mutate(N_conc_mgl = ifelse(flo == 0, 0, (total_N * 1000) / (flo * 86400))) %>% 
            mutate(total_P = solp + sedp) %>% 
            mutate(P_conc_mgl = ifelse(flo == 0, 0, (total_P * 1000) / (flo * 86400))) %>% 
            mutate(sed_conc_mgl = ifelse(flo == 0, 0, (sed * 1e6) / (flo * 86400)))
          
          # Calculate the frequency of exceeding the thresholds for each "unit" (channel)
          frequency_summary_mean <- channel_sd_day %>%
            group_by(name) %>%
            summarize(
              freq_below_threshold_lowQ = mean(flo < threshold_lowQ, na.rm = TRUE),
              freq_below_threshold_highQ = mean(flo < threshold_highQ, na.rm = TRUE),
              freq_below_threshold_N = mean(N_conc_mgl < threshold_N, na.rm = TRUE),
              freq_below_threshold_P = mean(P_conc_mgl < threshold_P, na.rm = TRUE),
              freq_below_threshold_Sed = mean(sed_conc_mgl < threshold_Sed, na.rm = TRUE)
            )
          df_out2[k,13] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_N[which(frequency_summary_mean$name==channel)]),3)
          df_out2[k,14] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_P[which(frequency_summary_mean$name==channel)]),3)
          df_out2[k,15] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_Sed[which(frequency_summary_mean$name==channel)]),3)
          df_out2[k,11] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_lowQ[which(frequency_summary_mean$name==channel)]),3)
          df_out2[k,12] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_highQ[which(frequency_summary_mean$name==channel)]),3)
        }
      }
      df_out[i,2:16] <- colMeans(df_out2)  
      df_out[i,17:31] <- lapply(df_out2, FUN='min')
      df_out[i,32:46] <- lapply(df_out2, FUN='max')
    }
  }

return(df_out)
}

# Indicators based on daily channel output
# Option b: only specific channels in cha.out (to be defined in object.prt file)
# !!! Don't forget to deactivate daily channel_sd printing in print.prt !!!
ind_cha_dayII <- function(path, 
                          channel, 
                          ind, 
                          threshold_lowQ=0.063,
                          threshold_highQ=2.01,
                          threshold_N=2.3, 
                          threshold_P=0.082, 
                          threshold_Sed=10,
                          hd=T,
                          ensemble=F){

  if(ensemble==F){
    df_out <- data.frame(scen_name=sapply(strsplit(path, split ="/"),tail,1), 
                         Q_max=NA, 
                         Q_p95=NA, 
                         Q_p90=NA,
                         Q_p50=NA, 
                         Q_p10=NA, 
                         Q_p05=NA, 
                         Q_min=NA,  
                         Q_maxmin=NA, 
                         Q_p95p05=NA, 
                         Q_p90p10=NA, 
                         Q_low_days=NA,
                         Q_high_days=NA,
                         Nconc_days=NA, 
                         Pconc_days=NA,
                         Sedconc_days=NA)
    
    for (i in 1:length(path)) {
      # Read file
      channel_sd_day <- read.table(paste0(path[i],'/cha_day.out'), h=T)
      names(channel_sd_day)[c(5,6)] <- c('type','name')
      
      if('Q_maxmin' %in% ind | 'Q_max' %in% ind  | 'Q_min' %in% ind | 'all' %in% ind){
        
        # convert m?/day to m?/s
        channel_sd_day$flo <- channel_sd_day$flo/86400
        
        # Handle zero flow (defining 0.00000001 cm?/s)
        channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
        
        # Group data by channel (replace "name" with the actual column name for channel)
        max_min_ratio <- channel_sd_day %>%
          group_by(name) %>% 
          summarise(
            max_discharge = max(flo, na.rm = TRUE),
            min_discharge = min(flo, na.rm = TRUE),
            extreme_streamflow_ratio = max_discharge / min_discharge
          )
        
        df_out[i,2] <- round(max_min_ratio$max_discharge[which(max_min_ratio$name==channel)],3)
        df_out[i,8] <- round(max_min_ratio$min_discharge[which(max_min_ratio$name==channel)],3)
        df_out[i,9] <- round(max_min_ratio$extreme_streamflow_ratio[which(max_min_ratio$name==channel)],3)
      }
      if('Q_p50' %in% ind | 'all' %in% ind){
        
        # convert m?/day to m?/s
        if('all' %in% ind == F) channel_sd_day$flo <- channel_sd_day$flo/86400
        
        # Handle zero flow (defining 0.00000001 cm?/s)
        channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
        
        # Group data by channel (replace "name" with the actual column name for channel)
        p50 <- channel_sd_day %>%
          group_by(name) %>% 
          summarise(
            p50_discharge = quantile(flo, probs = 0.50, na.rm = TRUE)
          )
        
        df_out[i,5] <- round(as.numeric(p50$p50_discharge[which(p50$name==channel)]),3)
      }
      if('Q_p95p05' %in% ind | 'Q_p95' %in% ind | 'Q_p05' %in% ind | 'all' %in% ind){
        
        # convert m?/day to m?/s
        if('all' %in% ind == F) channel_sd_day$flo <- channel_sd_day$flo/86400
        
        # Handle zero flow (defining 0.00000001 cm?/s)
        channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
        
        # Group data by channel (replace "name" with the actual column name for channel)
        Q_p95p05 <- channel_sd_day %>%
          group_by(name) %>%
          summarise(
            p05_discharge = quantile(flo, probs = 0.05, na.rm = TRUE),
            p95_discharge = quantile(flo, probs = 0.95, na.rm = TRUE),
            extreme_streamflow_ratio = p95_discharge / p05_discharge
          )
        df_out[i,3] <- round(as.numeric(Q_p95p05$p95_discharge[which(Q_p95p05$name==channel)]),3)
        df_out[i,7] <- round(as.numeric(Q_p95p05$p05_discharge[which(Q_p95p05$name==channel)]),3)
        df_out[i,10] <- round(as.numeric(Q_p95p05$extreme_streamflow_ratio[which(Q_p95p05$name==channel)]),3)
      }
      if('Q_p90p10' %in% ind | 'Q_p90' %in% ind | 'Q_p10' %in% ind | 'all' %in% ind){
        
        # convert m?/day to m?/s
        if('all' %in% ind == F) channel_sd_day$flo <- channel_sd_day$flo/86400
        
        # Handle zero flow (defining 0.00000001 cm?/s)
        channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
        
        # Group data by channel (replace "name" with the actual column name for channel)
        Q_p90p10 <- channel_sd_day %>%
          group_by(name) %>%
          summarise(
            p10_discharge = quantile(flo, probs = 0.10, na.rm = TRUE),
            p90_discharge = quantile(flo, probs = 0.90, na.rm = TRUE),
            extreme_streamflow_ratio = p90_discharge / p10_discharge
          )
        df_out[i,4] <- round(as.numeric(Q_p90p10$p90_discharge[which(Q_p90p10$name==channel)]),3)
        df_out[i,6] <- round(as.numeric(Q_p90p10$p10_discharge[which(Q_p90p10$name==channel)]),3)
        df_out[i,11] <- round(as.numeric(Q_p90p10$extreme_streamflow_ratio[which(Q_p90p10$name==channel)]),3)
      }
      if('Q_low_days' %in% ind | 'Q_high_days' %in% ind | 'Nconc_days' %in% ind | 'Pconc_days' %in% ind | 'Sedconc_days' %in% ind | 'all' %in% ind){
        
        # convert m?/day to m?/s
        if('all' %in% ind == F) channel_sd_day$flo <- channel_sd_day$flo/86400
        
        # Handle zero flow (defining 0.00000001 cm?/s)
        channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
        
        # Calculate the sum of no3, orgn, nh3, and no2
        channel_sd_day <- channel_sd_day %>%
          mutate(total_N = no3 + orgn + nh3 + no2) %>% 
          mutate(N_conc_mgl = ifelse(flo == 0, 0, (total_N * 1000) / (flo * 86400))) %>% 
          mutate(total_P = solp + sedp) %>% 
          mutate(P_conc_mgl = ifelse(flo == 0, 0, (total_P * 1000) / (flo * 86400))) %>% 
          mutate(sed_conc_mgl = ifelse(flo == 0, 0, (sed * 1e6) / (flo * 86400)))
        
        # Calculate the frequency of exceeding the thresholds for each "unit" (channel)
        frequency_summary_mean <- channel_sd_day %>%
          group_by(name) %>%
          summarize(
            freq_below_threshold_lowQ = mean(flo < threshold_lowQ, na.rm = TRUE),
            freq_below_threshold_highQ = mean(flo < threshold_highQ, na.rm = TRUE),
            freq_below_threshold_N = mean(N_conc_mgl < threshold_N, na.rm = TRUE),
            freq_below_threshold_P = mean(P_conc_mgl < threshold_P, na.rm = TRUE),
            freq_below_threshold_Sed = mean(sed_conc_mgl < threshold_Sed, na.rm = TRUE)
          )
        df_out[i,14] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_N[which(frequency_summary_mean$name==channel)]),3)
        df_out[i,15] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_P[which(frequency_summary_mean$name==channel)]),3)
        df_out[i,16] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_Sed[which(frequency_summary_mean$name==channel)]),3)
        df_out[i,12] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_lowQ[which(frequency_summary_mean$name==channel)]),3)
        df_out[i,13] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_highQ[which(frequency_summary_mean$name==channel)]),3)
      }
    }
  }
  
  if(ensemble==T){
    df_out <- data.frame(scen_name=sapply(strsplit(path, split ="/"),tail,1), 
                         Q_max=NA, 
                         Q_p95=NA, 
                         Q_p90=NA,
                         Q_p50=NA, 
                         Q_p10=NA, 
                         Q_p05=NA, 
                         Q_min=NA,  
                         Q_maxmin=NA, 
                         Q_p95p05=NA, 
                         Q_p90p10=NA, 
                         Q_low_days=NA,
                         Q_high_days=NA,
                         Nconc_days=NA, 
                         Pconc_days=NA,
                         Sedconc_days=NA,
                         Q_max_lower=NA,
                         Q_p95_lower=NA,
                         Q_p90_lower=NA,
                         Q_p50_lower=NA,
                         Q_p10_lower=NA,
                         Q_p05_lower=NA,
                         Q_min_lower=NA,
                         Q_maxmin_lower=NA,
                         Q_p95p05_lower=NA,
                         Q_p90p10_lower=NA,
                         Q_low_days_lower=NA,
                         Q_high_days_lower=NA,
                         Nconc_days_lower=NA,
                         Pconc_days_lower=NA,
                         Sedconc_days_lower=NA,
                         Q_max_upper=NA,
                         Q_p95_upper=NA,
                         Q_p90_upper=NA,
                         Q_p50_upper=NA,
                         Q_p10_upper=NA,
                         Q_p05_upper=NA,
                         Q_min_upper=NA,
                         Q_maxmin_upper=NA,
                         Q_p95p05_upper=NA,
                         Q_p90p10_upper=NA,
                         Q_low_days_upper=NA,
                         Q_high_days_upper=NA,
                         Nconc_days_upper=NA,
                         Pconc_days_upper=NA,
                         Sedconc_days_upper=NA)
    
    for (i in 1:length(path)) {
      # Read file
      path_ext <- dir(path[i], full.names = T)
      df_out2 <- data.frame(Q_max=NA,
                            Q_p95=NA,
                            Q_p90=NA,
                            Q_p50=NA,
                            Q_p10=NA,
                            Q_p05=NA,
                            Q_min=NA,
                            Q_maxmin=NA,
                            Q_p95p05=NA,
                            Q_p90p10=NA,
                            Q_low_days=NA,
                            Q_high_days=NA,
                            Nconc_days=NA,
                            Pconc_days=NA,
                            Sedconc_days=NA)
      for (k in 1:length(path_ext)){
        channel_sd_day <- read.table(paste0(path_ext[k],'/cha_day.out'), h=T)
        names(channel_sd_day)[c(5,6)] <- c('type','name')
        
        if('Q_maxmin' %in% ind | 'Q_max' %in% ind  | 'Q_min' %in% ind | 'all' %in% ind){
          
          # convert m?/day to m?/s
          channel_sd_day$flo <- channel_sd_day$flo/86400
          
          # Handle zero flow (defining 0.00000001 cm?/s)
          channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
          
          # Group data by channel (replace "name" with the actual column name for channel)
          max_min_ratio <- channel_sd_day %>%
            group_by(name) %>% 
            summarise(
              max_discharge = max(flo, na.rm = TRUE),
              min_discharge = min(flo, na.rm = TRUE),
              extreme_streamflow_ratio = max_discharge / min_discharge
            )
          
          df_out2[k,1] <- round(max_min_ratio$max_discharge[which(max_min_ratio$name==channel)],3)
          df_out2[k,7] <- round(max_min_ratio$min_discharge[which(max_min_ratio$name==channel)],3)
          df_out2[k,8] <- round(max_min_ratio$extreme_streamflow_ratio[which(max_min_ratio$name==channel)],3)
        }
        if('Q_p50' %in% ind | 'all' %in% ind){
          
          # convert m?/day to m?/s
          if('all' %in% ind == F) channel_sd_day$flo <- channel_sd_day$flo/86400
          
          # Handle zero flow (defining 0.00000001 cm?/s)
          channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
          
          # Group data by channel (replace "name" with the actual column name for channel)
          p50 <- channel_sd_day %>%
            group_by(name) %>% 
            summarise(
              p50_discharge = quantile(flo, probs = 0.50, na.rm = TRUE)
            )
          
          df_out2[k,4] <- round(as.numeric(p50$p50_discharge[which(p50$name==channel)]),3)
        }
        if('Q_p95p05' %in% ind | 'Q_p95' %in% ind | 'Q_p05' %in% ind | 'all' %in% ind){
          
          # convert m?/day to m?/s
          if('all' %in% ind == F) channel_sd_day$flo <- channel_sd_day$flo/86400
          
          # Handle zero flow (defining 0.00000001 cm?/s)
          channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
          
          # Group data by channel (replace "name" with the actual column name for channel)
          Q_p95p05 <- channel_sd_day %>%
            group_by(name) %>%
            summarise(
              p05_discharge = quantile(flo, probs = 0.05, na.rm = TRUE),
              p95_discharge = quantile(flo, probs = 0.95, na.rm = TRUE),
              extreme_streamflow_ratio = p95_discharge / p05_discharge
            )
          df_out2[k,2] <- round(as.numeric(Q_p95p05$p95_discharge[which(Q_p95p05$name==channel)]),3)
          df_out2[k,6] <- round(as.numeric(Q_p95p05$p05_discharge[which(Q_p95p05$name==channel)]),3)
          df_out2[k,9] <- round(as.numeric(Q_p95p05$extreme_streamflow_ratio[which(Q_p95p05$name==channel)]),3)
        }
        if('Q_p90p10' %in% ind | 'Q_p90' %in% ind | 'Q_p10' %in% ind | 'all' %in% ind){
          
          # convert m?/day to m?/s
          if('all' %in% ind == F) channel_sd_day$flo <- channel_sd_day$flo/86400
          
          # Handle zero flow (defining 0.00000001 cm?/s)
          channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
          
          # Group data by channel (replace "name" with the actual column name for channel)
          Q_p90p10 <- channel_sd_day %>%
            group_by(name) %>%
            summarise(
              p10_discharge = quantile(flo, probs = 0.10, na.rm = TRUE),
              p90_discharge = quantile(flo, probs = 0.90, na.rm = TRUE),
              extreme_streamflow_ratio = p90_discharge / p10_discharge
            )
          df_out2[k,3] <- round(as.numeric(Q_p90p10$p90_discharge[which(Q_p90p10$name==channel)]),3)
          df_out2[k,5] <- round(as.numeric(Q_p90p10$p10_discharge[which(Q_p90p10$name==channel)]),3)
          df_out2[k,10] <- round(as.numeric(Q_p90p10$extreme_streamflow_ratio[which(Q_p90p10$name==channel)]),3)
        }
        if('Q_low_days' %in% ind | 'Q_high_days' %in% ind | 'Nconc_days' %in% ind | 'Pconc_days' %in% ind | 'Sedconc_days' %in% ind | 'all' %in% ind){
          
          # convert m?/day to m?/s
          if('all' %in% ind == F) channel_sd_day$flo <- channel_sd_day$flo/86400
          
          # Handle zero flow (defining 0.00000001 cm?/s)
          channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
          
          # Calculate the sum of no3, orgn, nh3, and no2
          channel_sd_day <- channel_sd_day %>%
            mutate(total_N = no3 + orgn + nh3 + no2) %>% 
            mutate(N_conc_mgl = ifelse(flo == 0, 0, (total_N * 1000) / (flo * 86400))) %>% 
            mutate(total_P = solp + sedp) %>% 
            mutate(P_conc_mgl = ifelse(flo == 0, 0, (total_P * 1000) / (flo * 86400))) %>% 
            mutate(sed_conc_mgl = ifelse(flo == 0, 0, (sed * 1e6) / (flo * 86400)))
          
          # Calculate the frequency of exceeding the thresholds for each "unit" (channel)
          frequency_summary_mean <- channel_sd_day %>%
            group_by(name) %>%
            summarize(
              freq_below_threshold_lowQ = mean(flo < threshold_lowQ, na.rm = TRUE),
              freq_below_threshold_highQ = mean(flo < threshold_highQ, na.rm = TRUE),
              freq_below_threshold_N = mean(N_conc_mgl < threshold_N, na.rm = TRUE),
              freq_below_threshold_P = mean(P_conc_mgl < threshold_P, na.rm = TRUE),
              freq_below_threshold_Sed = mean(sed_conc_mgl < threshold_Sed, na.rm = TRUE)
            )
          df_out2[k,13] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_N[which(frequency_summary_mean$name==channel)]),3)
          df_out2[k,14] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_P[which(frequency_summary_mean$name==channel)]),3)
          df_out2[k,15] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_Sed[which(frequency_summary_mean$name==channel)]),3)
          df_out2[k,11] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_lowQ[which(frequency_summary_mean$name==channel)]),3)
          df_out2[k,12] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_highQ[which(frequency_summary_mean$name==channel)]),3)
        }
        
        }
      df_out[i,2:16] <- colMeans(df_out2)  
      df_out[i,17:31] <- lapply(df_out2, FUN='min')
      df_out[i,32:46] <- lapply(df_out2, FUN='max')
    }
  }
  
  return(df_out)
}

# water balance related indicators based on annual average hru output
ind_hru_aa_wb <- function(path, a = 'basin', ensemble = F){
  
  if(ensemble==F){
    df_out <- data.frame(scen_name=sapply(strsplit(path, split ="/"),tail,1), 
                         sw=NA, 
                         perc=NA)
    
    for (i in 1:length(path)) {
      # Read file
      hru_wb <- read_tbl('hru_wb_aa.txt', path[i], 3)
      hru_area <- read_tbl ('hru.con', path[i], 2)
      
      # Specify the columns you want to keep
      # Keep all hru_ls columns
      
      columns_to_keep <- c("jday", "mon", "day", "yr", "unit", "gis_id", "name", 
                           "sw_ave", "perc")
      
      # Create a new data frame with only the selected columns
      if(a == 'basin'){
        df_selected_hru_wb <- hru_wb[, columns_to_keep]
        idx <- c(1:dim(hru_wb)[1])
      }else{
        # Read in vector for agricultural area
        hru_agr <- read.table(paste0(path[i],'/hru_agr.txt'), h=T)
        idx <- hru_agr$hru_id
        df_selected_hru_wb <- hru_wb[idx, columns_to_keep]
      }
      
      # Add the HRU area to each HRU
      hru_wb <- hru_wb[idx,] %>% left_join(hru_area[idx,] %>% select(id, area), by = c("unit" = "id"))
      
      # Calculate the weighted values
      df_selected_hru_wb  <- df_selected_hru_wb  %>%
        mutate(weighted_sw = sw_ave * hru_wb$area,
               weighted_perc = perc * hru_wb$area)
      
      # Calculate the total weighted sum
      total_weighted_sw <- sum(df_selected_hru_wb$weighted_sw, na.rm = TRUE)
      total_weighted_perc <- sum(df_selected_hru_wb$weighted_perc, na.rm = TRUE)
      
      # Calculate the total area across all HRUs
      total_area <- sum(hru_wb$area, na.rm = TRUE)
      
      # Calculate the area-weighted averages
      df_out[i,2] <- round(total_weighted_sw / total_area,3)
      df_out[i,3] <- round(total_weighted_perc / total_area,3)
    }
  }   
  
  if(ensemble==T){
    df_out <- data.frame(scen_name=sapply(strsplit(path, split ="/"),tail,1), 
                         sw=NA, 
                         perc=NA,
                         sw_lower=NA,
                         perc_lower=NA,
                         sw_upper=NA,
                         perc_upper=NA)
    
    for (i in 1:length(path)) {
      path_ext <- dir(path[i], full.names = T)
      df_out2 <- data.frame(sw=NA, 
                            perc=NA)
      for (k in 1:length(path_ext)){
        # Read file
        hru_wb <- read_tbl('hru_wb_aa.txt', path_ext[k], 3)
        hru_area <- read_tbl ('hru.con', path_ext[k], 2)
        
        # Specify the columns you want to keep
        # Keep all hru_ls columns
        
        columns_to_keep <- c("jday", "mon", "day", "yr", "unit", "gis_id", "name", 
                             "sw_ave", "perc")
        
        # Create a new data frame with only the selected columns
        if(a == 'basin'){
          df_selected_hru_wb <- hru_wb[, columns_to_keep]
          idx <- c(1:dim(hru_wb)[1])
        }else{
          # Read in vector for agricultural area
          hru_agr <- read.table(paste0(path_ext[k],'/hru_agr.txt'), h=T)
          idx <- hru_agr$hru_id
          df_selected_hru_wb <- hru_wb[idx, columns_to_keep]
        }
        
        # Add the HRU area to each HRU
        hru_wb <- hru_wb[idx,] %>% left_join(hru_area[idx,] %>% select(id, area), by = c("unit" = "id"))
        
        # Calculate the weighted values
        df_selected_hru_wb  <- df_selected_hru_wb  %>%
          mutate(weighted_sw = sw_ave * hru_wb$area,
                 weighted_perc = perc * hru_wb$area)
        
        # Calculate the total weighted sum
        total_weighted_sw <- sum(df_selected_hru_wb$weighted_sw, na.rm = TRUE)
        total_weighted_perc <- sum(df_selected_hru_wb$weighted_perc, na.rm = TRUE)
        
        # Calculate the total area across all HRUs
        total_area <- sum(hru_wb$area, na.rm = TRUE)
        
        # Calculate the area-weighted averages
        df_out2[k,1] <- round(total_weighted_sw / total_area,3)
        df_out2[k,2] <- round(total_weighted_perc / total_area,3)
      }
      df_out[i,2:3] <- colMeans(df_out2)  
      df_out[i,4:5] <- lapply(df_out2, FUN='min')
      df_out[i,6:7] <- lapply(df_out2, FUN='max')
    }
  }
  
  return(df_out)
}

# nutrient and sediment related indicators based on annual average hru output
ind_hru_aa_nb <- function(path, a = 'basin', ensemble = F){
  
  if(ensemble==F){
    df_out <- data.frame(scen_name=sapply(strsplit(path, split ="/"),tail,1), 
                         N_loss=NA, 
                         P_loss=NA,
                         Sed_loss=NA,
                         N_loss_ratio=NA, 
                         P_loss_ratio=NA)
    
    for (i in 1:length(path)) {
      # Read file
      hru_ls <- read_tbl('hru_ls_aa.txt', path[i], 3)
      hru_nb <- read_tbl('hru_nb_aa.txt', path[i], 3)
      hru_pw <- read_tbl('hru_pw_aa.txt', path[i], 3)
      hru_area <- read_tbl ('hru.con', path[i], 2)
      
      # Specify the columns you want to keep
      # Keep all hru_ls columns
      
      columns_to_keep <- c("jday", "mon", "day", "yr", "unit", "gis_id", "name", 
                           "fertn", "fixn", "no3atmo", "nh4atmo", "fertp", "denit")
      
      
      # Create a new data frame with only the selected columns
      if(a == 'basin'){
        df_selected_hru_nb <- hru_nb[, columns_to_keep]
        idx <- c(1:dim(hru_nb)[1])
      }else{
        # Read in vector for agricultural area
        hru_agr <- read.table(paste0(path[i],'/hru_agr.txt'), h=T)
        idx <- hru_agr$hru_id
        df_selected_hru_nb <- hru_nb[idx, columns_to_keep]
      }
      
      # Calculate the sum of N inputs
      df_selected_hru_nb <- df_selected_hru_nb %>%
        mutate(N_inputs = fertn + fixn + no3atmo + nh4atmo)
      
      # Calculate N losses 
      # Add N_losses as a new column to df_selected_hru_nb
      
      hru_ls <- hru_ls[idx,] %>%
        mutate(N_losses = sedorgn + surqno3 + lat3no3 + tileno3 + hru_pw$percn[idx])
      
      # The P input is only fertp from hru_nb
      # Calculate the sum of P losses
      hru_ls <- hru_ls %>%
        mutate(P_losses = sedorgp + surqsolp + sedminp + tilelabp + lchlabp)
      
      # Add the HRU area to each HRU
      hru_ls <- hru_ls %>% left_join(hru_area[idx,] %>% select(id, area), by = c("unit" = "id"))
      
      # Calculate the weighted values for N and P inputs
      df_selected_hru_nb  <- df_selected_hru_nb  %>%
        mutate(weighted_N_inputs = N_inputs * hru_ls$area,
               weighted_P_inputs = fertp * hru_ls$area)
      # Calculate the weighted values for N, P and sediment losses
      hru_ls <- hru_ls %>%
        mutate(weighted_N_losses = N_losses * area,
               weighted_P_losses = P_losses * area,
               weighted_Sed_losses = sedyld * area)
      
      # Calculate the total weighted sum for N and P
      total_weighted_N_inputs <- sum(df_selected_hru_nb$weighted_N_inputs, na.rm = TRUE)
      total_weighted_N_losses <- sum(hru_ls$weighted_N_losses, na.rm = TRUE)
      total_weighted_P_inputs <- sum(df_selected_hru_nb$weighted_P_inputs, na.rm = TRUE)
      total_weighted_P_losses <- sum(hru_ls$weighted_P_losses, na.rm = TRUE)
      total_weighted_Sed_losses <- sum(hru_ls$weighted_Sed_losses, na.rm = TRUE)
      
      # Calculate the total area across all HRUs
      total_area <- sum(hru_ls$area, na.rm = TRUE)
      
      # Calculate the area-weighted averages for N and P
      area_weighted_average_N_inputs <- total_weighted_N_inputs / total_area
      area_weighted_average_P_inputs <- total_weighted_P_inputs / total_area
      area_weighted_average_N_losses <- total_weighted_N_losses / total_area
      area_weighted_average_P_losses <- total_weighted_P_losses / total_area
      area_weighted_average_Sed_losses <- total_weighted_Sed_losses / total_area
      
      df_out[i,2] <- round(area_weighted_average_N_losses,3)
      df_out[i,3] <- round(area_weighted_average_P_losses,3)
      df_out[i,4] <- round(area_weighted_average_Sed_losses,3)
      df_out[i,5] <- round(area_weighted_average_N_losses/area_weighted_average_N_inputs,3)
      df_out[i,6] <- round(area_weighted_average_P_losses/area_weighted_average_P_inputs,3)
      
    }
  }
  
  if(ensemble==T){
    df_out <- data.frame(scen_name=sapply(strsplit(path, split ="/"),tail,1), 
                         N_loss=NA, 
                         P_loss=NA,
                         Sed_loss=NA,
                         N_loss_ratio=NA, 
                         P_loss_ratio=NA,
                         N_loss_lower=NA, 
                         P_loss_lower=NA,
                         Sed_loss_lower=NA,
                         N_loss_ratio_lower=NA, 
                         P_loss_ratio_lower=NA,
                         N_loss_upper=NA, 
                         P_loss_upper=NA,
                         Sed_loss_upper=NA,
                         N_loss_ratio_upper=NA, 
                         P_loss_ratio_upper=NA)
    
    for (i in 1:length(path)) {
      path_ext <- dir(path[i], full.names = T)
      df_out2 <- data.frame(N_loss=NA, 
                            P_loss=NA,
                            Sed_loss=NA,
                            N_loss_ratio=NA, 
                            P_loss_ratio=NA)
      for (k in 1:length(path_ext)){
        # Read file
        hru_ls <- read_tbl('hru_ls_aa.txt', path_ext[k], 3)
        hru_nb <- read_tbl('hru_nb_aa.txt', path_ext[k], 3)
        hru_pw <- read_tbl('hru_pw_aa.txt', path_ext[k], 3)
        hru_area <- read_tbl ('hru.con', path_ext[k], 2)
        
        # Specify the columns you want to keep
        # Keep all hru_ls columns
        
        columns_to_keep <- c("jday", "mon", "day", "yr", "unit", "gis_id", "name", 
                             "fertn", "fixn", "no3atmo", "nh4atmo", "fertp", "denit")
        
        
        # Create a new data frame with only the selected columns
        if(a == 'basin'){
          df_selected_hru_nb <- hru_nb[, columns_to_keep]
          idx <- c(1:dim(hru_nb)[1])
        }else{
          # Read in vector for agricultural area
          hru_agr <- read.table(paste0(path_ext[k],'/hru_agr.txt'), h=T)
          idx <- hru_agr$hru_id
          df_selected_hru_nb <- hru_nb[idx, columns_to_keep]
        }
        
        # Calculate the sum of N inputs
        df_selected_hru_nb <- df_selected_hru_nb %>%
          mutate(N_inputs = fertn + fixn + no3atmo + nh4atmo)
        
        # Calculate N losses 
        # Add N_losses as a new column to df_selected_hru_nb
        
        hru_ls <- hru_ls[idx,] %>%
          mutate(N_losses = sedorgn + surqno3 + lat3no3 + tileno3 + hru_pw$percn[idx])
        
        # The P input is only fertp from hru_nb
        # Calculate the sum of P losses
        hru_ls <- hru_ls %>%
          mutate(P_losses = sedorgp + surqsolp + sedminp + tilelabp + lchlabp)
        
        # Add the HRU area to each HRU
        hru_ls <- hru_ls %>% left_join(hru_area[idx,] %>% select(id, area), by = c("unit" = "id"))
        
        # Calculate the weighted values for N and P inputs
        df_selected_hru_nb  <- df_selected_hru_nb  %>%
          mutate(weighted_N_inputs = N_inputs * hru_ls$area,
                 weighted_P_inputs = fertp * hru_ls$area)
        # Calculate the weighted values for N, P and sediment losses
        hru_ls <- hru_ls %>%
          mutate(weighted_N_losses = N_losses * area,
                 weighted_P_losses = P_losses * area,
                 weighted_Sed_losses = sedyld * area)
        
        # Calculate the total weighted sum for N and P
        total_weighted_N_inputs <- sum(df_selected_hru_nb$weighted_N_inputs, na.rm = TRUE)
        total_weighted_N_losses <- sum(hru_ls$weighted_N_losses, na.rm = TRUE)
        total_weighted_P_inputs <- sum(df_selected_hru_nb$weighted_P_inputs, na.rm = TRUE)
        total_weighted_P_losses <- sum(hru_ls$weighted_P_losses, na.rm = TRUE)
        total_weighted_Sed_losses <- sum(hru_ls$weighted_Sed_losses, na.rm = TRUE)
        
        # Calculate the total area across all HRUs
        total_area <- sum(hru_ls$area, na.rm = TRUE)
        
        # Calculate the area-weighted averages for N and P
        area_weighted_average_N_inputs <- total_weighted_N_inputs / total_area
        area_weighted_average_P_inputs <- total_weighted_P_inputs / total_area
        area_weighted_average_N_losses <- total_weighted_N_losses / total_area
        area_weighted_average_P_losses <- total_weighted_P_losses / total_area
        area_weighted_average_Sed_losses <- total_weighted_Sed_losses / total_area
        
        df_out2[k,1] <- round(area_weighted_average_N_losses,3)
        df_out2[k,2] <- round(area_weighted_average_P_losses,3)
        df_out2[k,3] <- round(area_weighted_average_Sed_losses,3)
        df_out2[k,4] <- round(area_weighted_average_N_losses/area_weighted_average_N_inputs,3)
        df_out2[k,5] <- round(area_weighted_average_P_losses/area_weighted_average_P_inputs,3)
      }
      df_out[i,2:6] <- colMeans(df_out2)  
      df_out[i,7:11] <- lapply(df_out2, FUN='min')
      df_out[i,12:16] <- lapply(df_out2, FUN='max')   
    }
  }
  
  return(df_out)
}

# water balance related indicators based on monthly hru outputs
ind_hru_mon_wb <- function(path, period = c(5:9), a = 'basin', ensemble = F, nrows = length(measr.list)){
  
  if(ensemble==F){
    sw_colnames <- paste0('sw_', map_chr(period, ~paste(.x, collapse = "_")))
    df_out <- data.frame(matrix(data=NA, nrow=nrows, ncol=length(sw_colnames)+1))
    names(df_out) <- c('scen_name', sw_colnames)
    df_out$scen_name <- sapply(strsplit(path, split ="/"),tail,1)
                         
    for (i in 1:length(path)) {
      # Read file
      hru_wb <- read_tbl('hru_wb_mon.txt', path[i], 3)
      hru_area <- read_tbl ('hru.con', path[i], 2)
      
      # Specify the columns you want to keep for hru_nb
      # Keep all hru_ls columns
      
      columns_to_keep <- c("jday", "mon", "day", "yr", "unit", "gis_id", "name", 
                           "sw_ave")
      
      for(k in 1:length(period)){
        # Create a new data frame with only the selected rows and columns
        if(a == 'basin'){
          df_selected_hru_wb <- hru_wb[which(hru_wb$mon %in% period[[k]]), columns_to_keep]
          idx <- c(1:dim(hru_area)[1])
        }else{
          # Read in vector for agricultural area
          hru_agr <- read.table(paste0(path[i],'/hru_agr.txt'), h=T)
          idx <- hru_agr$hru_id
          df_selected_hru_wb <- hru_wb[which(hru_wb$unit %in% idx & hru_wb$mon %in% period[[k]]), columns_to_keep]
        }
        
        # Calculate the area weighted values
        df_selected_hru_wb  <- df_selected_hru_wb  %>%
          left_join(hru_area %>% select(id, area), by = c("unit" = "id")) %>% 
          mutate(weighted_sw = sw_ave * area)
        
        # Number of months and years
        n_mon <- length(period[[k]])
        n_years <- length(unique(df_selected_hru_wb$yr))
        
        # Calculate the total weighted sum
        total_weighted_sw <- sum(df_selected_hru_wb$weighted_sw, na.rm = TRUE)/n_mon/n_years
        
        # Calculate the total area across all HRUs
        total_area <- sum(hru_area$area[idx], na.rm = TRUE)
        
        # Calculate the area-weighted averages for N and P
        df_out[i,k+1] <- round(total_weighted_sw / total_area,3)
      }
    }
      
  }
  
  if(ensemble==T){
    sw_colnames <- paste0('sw_', map_chr(period, ~paste(.x, collapse = "_")))
    df_out <- data.frame(matrix(data=NA, nrow=length(measr.list), ncol=length(sw_colnames)*3+1))
    names(df_out) <- c('scen_name', sw_colnames, paste0(sw_colnames, '_lower'), paste0(sw_colnames, '_upper'))
    df_out$scen_name <- sapply(strsplit(path, split ="/"),tail,1)
    
    for (i in 1:length(path)) {
      path_ext <- dir(path[i], full.names = T)
      df_out2 <- data.frame(matrix(data=NA, nrow=length(measr.list), ncol=length(sw_colnames)))
      names(df_out2) <- c(sw_colnames)
      for (j in 1:length(path_ext)){
        # Read file
        hru_wb <- read_tbl('hru_wb_mon.txt', path_ext[j], 3)
        hru_area <- read_tbl ('hru.con', path_ext[j], 2)
        
        # Specify the columns you want to keep for hru_nb
        # Keep all hru_ls columns
        
        columns_to_keep <- c("jday", "mon", "day", "yr", "unit", "gis_id", "name", 
                             "sw_ave")
        
        for(k in 1:length(period)){
          # Create a new data frame with only the selected rows and columns
          if(a == 'basin'){
            df_selected_hru_wb <- hru_wb[which(hru_wb$mon %in% period[[k]]), columns_to_keep]
            idx <- c(1:dim(hru_area)[1])
          }else{
            # Read in vector for agricultural area
            hru_agr <- read.table(paste0(path_ext[j],'/hru_agr.txt'), h=T)
            idx <- hru_agr$hru_id
            df_selected_hru_wb <- hru_wb[which(hru_wb$unit %in% idx & hru_wb$mon %in% period[[k]]), columns_to_keep]
          }
          
          # Calculate the area weighted values
          df_selected_hru_wb  <- df_selected_hru_wb  %>%
            left_join(hru_area %>% select(id, area), by = c("unit" = "id")) %>% 
            mutate(weighted_sw = sw_ave * area)
          
          # Number of months and years
          n_mon <- length(period[[k]])
          n_years <- length(unique(df_selected_hru_wb$yr))
          
          # Calculate the total weighted sum
          total_weighted_sw <- sum(df_selected_hru_wb$weighted_sw, na.rm = TRUE)/n_mon/n_years
          
          # Calculate the total area across all HRUs
          total_area <- sum(hru_area$area[idx], na.rm = TRUE)
          
          # Calculate the area-weighted averages for N and P
          df_out2[j,k] <- round(total_weighted_sw / total_area,3)
        }
      }
      df_out[i,2:(length(df_out2)+1)] <- colMeans(df_out2)  
      df_out[i,(length(df_out2)+2):(2*length(df_out2)+1)] <- lapply(df_out2, FUN='min')
      df_out[i,(2*length(df_out2)+2):(3*length(df_out2)+1)] <- lapply(df_out2, FUN='max')
      }
  }
  
  return(df_out)
}

# crop yield related indicators (grain units) based on annual average hru output
ind_bsn_aa_crp <- function(path, crop_sel, out_type, grain_units, ensemble = F){
  
  if (ensemble == F){
    if (out_type == "yield") {
      df_out <- data.frame(scen_name=sapply(strsplit(path, split ="/"),tail,1), grain_units_aa=NA)
    } else {
      df_out <- data.frame(scen_name=sapply(strsplit(path, split ="/"),tail,1), crops_ha_aa=NA)
    }
    for (i in 1:length(path)) {
      # read output file
      crop_aa <- read_tbl('basin_crop_yld_aa.txt', path[i], 2)
      # Index for reading
      if (out_type == "yield") { 
        crop_sel <- names(grain_units)[match(crop_sel, names(grain_units))]
        idx <- pmatch(crop_sel, crop_aa$plant_name, duplicates.ok = T)
        # Convert yield to grain units and sum it up
        crop_yld_gu <- round(sum(crop_aa$`yld(t)`[idx] * grain_units, na.rm = T),3)
        # collect grain unit values in df_out
        df_out[i,2] <- crop_yld_gu
      } else {
        idx <- pmatch(crop_sel, crop_aa$plant_name, duplicates.ok = T)
        # sum up hectare values and collect in df_out
        df_out[i,2] <- round(sum(crop_aa$`harv_area(ha)`[idx], na.rm = T),2)
      }
  }
  }
  
  if (ensemble == T){
    if (out_type == "yield") {
      df_out <- data.frame(scen_name=sapply(strsplit(path, split ="/"),tail,1), 
                           grain_units_aa=NA, grain_units_aa_lower=NA, grain_units_aa_upper=NA)
    } else {
      df_out <- data.frame(scen_name=sapply(strsplit(path, split ="/"),tail,1), 
                           crops_ha_aa=NA, crops_ha_aa_lower=NA, crops_ha_aa_upper=NA)
    }
    
    for (i in 1:length(path)) {
      path_ext <- dir(path[i], full.names = T)
      if (out_type == "yield") df_out2 <- data.frame(grain_units_aa=NA) else{
        df_out2 <- data.frame(crops_ha_aa=NA)
      }
      
      for (k in 1:length(path_ext)){
        # read output file
        crop_aa <- read_tbl('basin_crop_yld_aa.txt', path_ext[k], 2)
        if (out_type == "yield"){
        # Index for reading
          crop_sel <- names(grain_units)[match(crop_sel, names(grain_units))]
          idx <- pmatch(crop_sel, crop_aa$plant_name, duplicates.ok = T)
          # Convert yield to grain units and sum it up
          crop_yld_gu <- round(sum(crop_aa$`yld(t)`[idx] * grain_units, na.rm = T),3)
          # collect grain unit values in df_out
          df_out2[k,1] <- crop_yld_gu
        } else {
          idx <- pmatch(crop_sel, crop_aa$plant_name, duplicates.ok = T)
          # sum up hectare values and collect in df_out
          df_out2[k,1] <- round(sum(crop_aa$`harv_area(ha)`[idx], na.rm = T),2)
        }
      }
      df_out[i,2] <- colMeans(df_out2)  
      df_out[i,3] <- lapply(df_out2, FUN='min')
      df_out[i,4] <- lapply(df_out2, FUN='max')   
    }
  }
  
  return(df_out)
}

# crop yield related indicators (area and yield) based on annual average hru output
ind_bsn_aa_crp_ha_Y <- function(path, crop_sel, ensemble=F){
  
  if (ensemble == F){
    # Read files & prepare data
    df_out <- data.frame(matrix(ncol = length(crop_sel)*2+1, nrow = length(path)))
    colnames(df_out) <- c("scen_name", paste(crop_sel, "_ha", sep=""), paste(crop_sel, "_yld_t_ha", sep=""))
    df_out$scen_name <- sapply(strsplit(path, split ="/"),tail,1)
    for (i in 1:length(path)) {
      # read output file
      crop_aa <- read_tbl('basin_crop_yld_aa.txt', path[i], 2)
      # Index for reading
      idx <- pmatch(crop_sel, crop_aa$plant_name, duplicates.ok = T)
      # collect ha and yield values in df_out
      df_out[i,2:ncol(df_out)] <- round(c(crop_aa$`harv_area(ha)`[idx], crop_aa$`yld(t/ha)`[idx]), 2)
    }
  }
  
  if (ensemble == T){
    # Read files & prepare data
    df_out <- data.frame(matrix(ncol = length(crop_sel)*4+1, nrow = length(path)))
    colnames(df_out) <- c("scen_name", 
                          paste(crop_sel, "_ha", sep=""), 
                          paste(crop_sel, "_yld_t_ha", sep=""),
                          paste(crop_sel, "_yld_t_ha_lower", sep=""),
                          paste(crop_sel, "_yld_t_ha_upper", sep=""))
    df_out$scen_name <- sapply(strsplit(path, split ="/"),tail,1)
    for (i in 1:length(path)) {
      path_ext <- dir(path[i], full.names = T)
      df_out2 <- data.frame(matrix(ncol = length(crop_sel)*2, nrow = length(path)))
      for (k in 1:length(path_ext)){
        # read output file
        crop_aa <- read_tbl('basin_crop_yld_aa.txt', path_ext[k], 2)
        # Index for reading
        idx <- pmatch(crop_sel, crop_aa$plant_name, duplicates.ok = T)
        # collect ha and yield values in df_out
        df_out2[k,1:ncol(df_out2)] <- round(c(crop_aa$`harv_area(ha)`[idx], crop_aa$`yld(t/ha)`[idx]), 2)
      }
      df_out[i,2:(length(df_out2)+1)] <- colMeans(df_out2)  
      df_out[i,(length(df_out2)+2):(1.5*length(df_out2)+1)] <- lapply(df_out2[(0.5*length(df_out2)+1):length(df_out2)], FUN='min')
      df_out[i,(1.5*length(df_out2)+2):(2*length(df_out2)+1)] <- lapply(df_out2[(0.5*length(df_out2)+1):length(df_out2)], FUN='max')
    }
  }
 
  return(df_out)
}

# Load in observations
get_obs <- function(obs_path, obs_type, sim_period) {
  
  start_date <- ISOdate(sim_period[1], 1, 1)
  end_date <- ISOdate(sim_period[2], 12, 31)
  
  if(obs_type == 'q'){
    gauge_file_paths <- dir(obs_path, pattern = obs_type, full.names = T)
    gauge_file_names <- dir(obs_path, pattern = obs_type, full.names = F)
    
    pos1 <- unlist(gregexpr('cha', gauge_file_names))+3
    cha <- substring(gauge_file_names, pos1, pos1+3)
    
    obs <- setNames(replicate(length(cha),data.frame()),cha)
    
    for(i in 1:length(obs)){
      obs[[i]] <- read.csv(gauge_file_paths[i]) %>%
        mutate(date = as.Date(date)) %>% 
        filter(date >= as.Date(start_date),
               date <= as.Date(end_date))
    }
  }else{
    gauge_file_paths <- dir(obs_path, pattern = obs_type, full.names = T)
    gauge_file_names <- dir(obs_path, pattern = obs_type, full.names = F)
    flow_gauge_file_paths <- dir(obs_path, pattern = 'q', full.names = T)
    flow_gauge_file_names <- dir(obs_path, pattern = 'q', full.names = F)
    
    pos1 <- unlist(gregexpr('cha', gauge_file_names))+3
    pos1_flow <- unlist(gregexpr('cha', flow_gauge_file_names))+3
    cha <- substring(gauge_file_names, pos1, pos1+3)
    cha_flow <- substring(flow_gauge_file_names, pos1_flow, pos1_flow+3)
    
    idx <- which(cha %in% cha_flow)
    idx_flow <- which(cha_flow %in% cha)
    cha <- cha[idx]
    
    gauge_file_paths <- gauge_file_paths[idx]
    flow_gauge_file_paths <- flow_gauge_file_paths[idx_flow]
    
    obs <- obs_flow <- setNames(replicate(length(cha),data.frame()),cha)
    
    for(i in 1:length(obs)){
      obs_flow[[i]] <- read.csv(flow_gauge_file_paths[i]) %>%
        filter(date >= as.Date(start_date),
               date <= as.Date(end_date))
      obs[[i]] <- read.csv(gauge_file_paths[i]) %>%
        filter(date %in% obs_flow[[i]]$date)
      obs_flow[[i]] <- obs_flow[[i]] %>% 
        filter(date %in% obs[[i]]$date)
      if(obs_type == 'tss'){
        obs[[i]] <- obs[[i]] %>% 
          mutate(date = as.Date(date),
                 value = value/1e6 * obs_flow[[i]]$value * 86400)
      }else{
        obs[[i]] <- obs[[i]] %>% 
          mutate(date = as.Date(date),
                 value = value/1e3 * obs_flow[[i]]$value * 86400)
      }
    }
  }
  return(obs)
}

# Calculate model performance metrics
calc_performance <- function(sim, obs){
  
  idx <- which(sim$date %in% obs$date)
  
  # Calculate flow duration curves (FDC) for observed data and the simulations.
  fdc_obs <- calc_fdc(obs$value)
  fdc_sim <- calc_fdc(select(sim[idx,], - date))
  
  # Calculate the ratio of RSME and standard deviation for different segments
  # of the FDC (same as in the publications of the Kiel working group).
  rsr_fdc <- calc_fdc_rsr(fdc_sim, fdc_obs, c(5, 20, 70, 95))
  
  # I was testing if not FDC but different flow segments (keeping the time 
  # component) would be a better criterion.
  p <- c(5, 20, 70, 95)
  p_lbl <- c('p_0_5', 'p_5_20', 'p_20_70', 'p_70_95', 'p_95_100')
  fdc_thrs <- c(max(fdc_obs$value),
                approx(fdc_obs$p, fdc_obs$value, p)$y,
                -0.1)
  # Separate the hydrograph into high medium and low flows.
  obs_sep <- map2(fdc_thrs[1:(length(fdc_thrs) - 1)],
                  fdc_thrs[2:length(fdc_thrs)],
                  ~ mutate(obs, value = ifelse(value <= .x & value > .y, value, NA))) %>% 
    map2(., p_lbl, ~ set_names(.x, c('date', .y))) %>% 
    reduce(., left_join, by = 'date')
  
  # Calculate further criteria, e.g. NSE, KGE, pbias, whatever...
  nse_q <- map_dbl(select(sim[idx,], - date), ~NSE(.x, obs$value))
  rsr_vh <- map_dbl(select(sim[idx,], - date), ~rsr(.x, obs_sep$p_0_5))
  rsr_h <- map_dbl(select(sim[idx,], - date), ~rsr(.x, obs_sep$p_5_20))
  rsr_m <- map_dbl(select(sim[idx,], - date), ~rsr(.x, obs_sep$p_20_70))
  rsr_l <- map_dbl(select(sim[idx,], - date), ~rsr(.x, obs_sep$p_70_95))
  rsr_vl <- map_dbl(select(sim[idx,], - date), ~rsr(.x, obs_sep$p_95_100))
  kge_q <- map_dbl(select(sim[idx,], - date), ~KGE(.x, obs$value))
  pb_q  <- map_dbl(select(sim[idx,], - date), ~pbias(.x, obs$value))
  rd_q <- map_dbl(select(sim[idx,], - date), ~rd(.x, obs$value))
  
  # Put all of your selected criteria together in one objectives table.
  obj_tbl <- bind_cols(run = names(nse_q),
                       nse = nse_q,
                       rsr_vh = -rsr_vh,
                       rsr_h = -rsr_h,
                       rsr_m = -rsr_m,
                       rsr_l = -rsr_l,
                       rsr_vl = -rsr_vl,
                       kge = kge_q, 
                       pbias = pb_q,
                       rsr_0_5 = - rsr_fdc$p_0_5, 
                       rsr_5_20 = - rsr_fdc$p_5_20, 
                       rsr_20_70 = - rsr_fdc$p_20_70,
                       rsr_70_95 = - rsr_fdc$p_70_95,
                       rsr_95_100 = -rsr_fdc$p_95_100,
                       rd = rd_q) 
  return(obj_tbl)
}

calc_fdc <- function(x) {
  if(is.vector(x)) {
    x <- tibble(value = x)
  }
  
  n <- nrow(x)
  
  x %>%
    apply(., 2, sort, decreasing = TRUE) %>%
    as_tibble(.) %>%
    mutate(p = 100 * 1:n / (n + 1), .before = 1)
}

calc_fdc_rsr <- function(fdc_sim, fdc_obs, quantile_splits, out_tbl = 'long') {
  if(all(quantile_splits <= 1)) {
    quantile_splits <- 100 * quantile_splits
  }
  quantile_splits <- sort(unique(c(0, 100, quantile_splits)))
  p_cuts <- cut(fdc_obs$p, quantile_splits)
  
  obs <- split(select(fdc_obs, -p), p_cuts)
  sim <- split(select(fdc_sim, -p), p_cuts)
  
  rsr_list <- map2(sim, obs, ~ rsr_df(.x, .y[[1]]))
  
  if(out_tbl == 'long') {
    n_col <- length(quantile_splits) - 1
    col_names <- paste0('p_', quantile_splits[1:n_col],
                        '_',  quantile_splits[2:(n_col + 1)])
    rsr <- bind_cols(rsr_list) %>%
      set_names(col_names) %>%
      mutate(., run = names(fdc_sim)[2:ncol(fdc_sim)], .before = 1)
  } else {
    rsr <- rsr_list %>%
      bind_rows(.) %>%
      mutate(p = unique(p_cuts), .before = 1)
  }
  return(rsr)
}

rsr_df <- function(df_sim, v_obs) {
  map_dbl(df_sim, ~ rsr(.x, v_obs))
}

plot_selected_sim <- function(sim, obs = NULL, run_ids = NULL, run_sel = NULL, plot_bands = TRUE, 
                              x_label = 'Date', y_label = "Discharge (m<sup>3</sup> s<sup>-1</sup>)") {
  
  if(is.null(run_ids) & is.null(run_sel)) {
    stop("At least one of 'run_ids' or 'run_sel' must be provided.")
  }
  if(!is.Date(sim[[1]])){
    stop("The first column of 'sim' must by of type 'Date'.")
  }
  
  dy_tbl <- select(sim, date)
  
  if(!is.null(obs)) {
    if(!is.Date(obs[[1]])){
      stop("The first column of 'obs' must by of type 'Date'.")
    }
    names(obs) <- c('date', 'obs')
    dy_tbl <- left_join(dy_tbl, obs, by = 'date')
  } 
  
  nchar_run <- nchar(names(sim)[2]) - 4
  
  if(!is.null(run_sel)) {
    run_sel <- paste0('run_', sprintf(paste0('%0', nchar_run, 'd'), run_sel))
    dy_tbl <- add_column(dy_tbl, sim_sel = sim[, run_sel])
  }
  
  if(!is.null(run_ids)) {
    run_ids <- paste0('run_', sprintf(paste0('%0', nchar_run, 'd'), run_ids))
    sim_ids <- sim[, run_ids]
    if(plot_bands) {
      sim_upr <- apply(sim_ids, 1, max)
      sim_lwr <- apply(sim_ids, 1, min)
      dy_tbl <- add_column(dy_tbl, upr   = sim_upr, lwr   = sim_lwr,
                           upr_l = sim_upr, lwr_l = sim_lwr)
    } else {
      dy_tbl <- bind_cols(dy_tbl, sim_ids)
    }
  }
  
  
  dy_xts <- xts(dy_tbl, order.by = dy_tbl$date)
  
  dy_plt <- dy_xts %>% 
    dygraph(., xlab = x_label, ylab = y_label) %>% 
    dyRangeSelector(height = 30)
  
  if(!is.null(obs)) {
    dy_plt <- dy_plt %>% 
      dySeries('obs', color = 'black', drawPoints = TRUE, 
               pointSize = 2, strokeWidth = 0.75)
  }
  
  if(plot_bands) {
    if(!is.null(run_ids) & !is.null(run_sel)) {
      dy_plt <- dy_plt %>% 
        dySeries(c("lwr", "sim_sel", "upr"), label = "sim_sel", 
                 color =  "#A50F15", strokeWidth = 1.2, 
                 drawPoints = TRUE, pointSize = 2) %>% 
        dySeries("lwr_l", label = "lower", color =  "#CB181D", 
                 strokePattern = 'dashed') %>% 
        dySeries("upr_l", label = "upper", color =  "#CB181D", 
                 strokePattern = 'dashed')
    } else if (is.null(run_ids)) {
      dy_plt <- dy_plt %>% 
        dySeries("sim_sel", label = "sim_sel", 
                 color =  "#A50F15", strokeWidth = 1.2, 
                 drawPoints = TRUE, pointSize = 2) 
    } else {
      dy_plt <- dy_plt %>% 
        dySeries(c("lwr", "upr_l", "upr"), label = "upper", 
                 color =  "#CB181D", strokePattern = 'dashed') %>% 
        dySeries("lwr_l", label = "lower", color =  "#CB181D", 
                 strokePattern = 'dashed')
    }
  } else {
    col_pal <- c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
                 "#A65628", "#F781BF", "#999999", "#1B9E77", "#D9D9D9")
    if(!is.null(run_ids)) {
      for(i in 1:length(run_ids)) {
        dy_plt <- dy_plt %>% 
          dySeries(run_ids[i], label = run_ids[i], 
                   color =  col_pal[i], strokeWidth = 1)
      }
      dy_plt <- dy_plt %>% 
        dyHighlight(highlightSeriesBackgroundAlpha = 0.6,
                    hideOnMouseOut = TRUE)
    }
    if (!is.null(run_sel)) {
      dy_plt <- dy_plt %>% 
        dySeries("sim_sel", label = "sim_sel", 
                 color =  "#A50F15", strokeWidth = 1.2, 
                 drawPoints = TRUE, pointSize = 2) 
    } else {
      dy_plt <- dy_plt %>% 
        dySeries(c("lwr", "upr_l", "upr"), label = "upper", 
                 color =  "#CB181D", strokePattern = 'dashed') %>% 
        dySeries("lwr_l", label = "lower", color =  "#CB181D", 
                 strokePattern = 'dashed')
    }
  }
  
  return(dy_plt)
}

## Load required packages
foo1(c('dplyr' , 'readr' , 'tidyverse', 'data.table', 'remotes', 'devtools', 
       'xts', 'dygraphs', 'R.utils', 'foreach', 'doParallel', 'data.table', 
       'ggplot2', 'fmsb'))
