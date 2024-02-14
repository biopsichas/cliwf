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

## HRUs numbers to get
nb_hru <- 1:10240

##Results to extract from SWATrunR scenario runs
results_to_save <- list(channel_flo = define_output('channel_sd_day', 'flo_out', c(43, 52)), ##Your reaches
                        aqu_flo = define_output('basin_aqu_mon', 'flo', 1),
                        aqu_dep_wt = define_output('basin_aqu_mon', 'dep_wt', 1),
                        aqu_stor  = define_output('basin_aqu_mon', 'stor', 1),
                        aqu_rchrg = define_output('basin_aqu_mon', 'rchrg', 1),
                        aqu_seep = define_output('basin_aqu_mon', 'seep', 1),
                        aqu_revap = define_output('basin_aqu_mon', 'revap', 1),
                        bpw_lai = define_output("basin_pw_mon", "lai", 1),
                        bpw_bioms = define_output("basin_pw_mon", "bioms", 1),
                        bpw_yield = define_output("basin_pw_mon", "yield", 1),
                        bpw_strsw = define_output("basin_pw_mon", "strsw", 1),
                        bpw_strsa = define_output("basin_pw_mon", "strsa", 1),
                        bpw_strstmp = define_output("basin_pw_mon", "strstmp", 1),
                        bpw_strsn = define_output("basin_pw_mon", "strsn", 1),
                        hpw_lai = define_output("hru_pw_aa", "lai", nb_hru),
                        hpw_bioms = define_output("hru_pw_aa", "bioms", nb_hru),
                        hpw_yield = define_output("hru_pw_aa", "yield", nb_hru),
                        hpw_strsw = define_output("hru_pw_aa", "strsw", nb_hru),
                        hpw_strsa = define_output("hru_pw_aa", "strsa", nb_hru),
                        hpw_strstmp = define_output("hru_pw_aa", "strstmp", nb_hru),
                        hpw_strsn = define_output("hru_pw_aa", "strsn", nb_hru),
                        bwb_precip = define_output("basin_wb_mon", "precip", 1),
                        bwb_snofall = define_output("basin_wb_mon", "snofall", 1),
                        bwb_snomlt = define_output("basin_wb_mon", "snomlt", 1),
                        bwb_wateryld = define_output("basin_wb_mon", "wateryld", 1),
                        bwb_perc = define_output("basin_wb_mon", "perc", 1),
                        bwb_et = define_output("basin_wb_mon", "et", 1),
                        bwb_ecanopy = define_output("basin_wb_mon", "ecanopy", 1),
                        bwb_eplant = define_output("basin_wb_mon", "eplant", 1),
                        bwb_esoil = define_output("basin_wb_mon", "esoil", 1),
                        bwb_cn = define_output("basin_wb_mon", "cn", 1),
                        bwb_sw_ave = define_output("basin_wb_mon", "sw_ave", 1),
                        bwb_sw_300 = define_output("basin_wb_mon", "sw_300", 1),
                        bwb_snopack = define_output("basin_wb_mon", "snopack", 1),
                        bwb_pet = define_output("basin_wb_mon", "pet", 1),
                        bwb_qtile = define_output("basin_wb_mon", "qtile", 1),
                        bwb_surq_cha = define_output("basin_wb_mon", "surq_cha", 1),
                        bwb_latq_cha = define_output("basin_wb_mon", "latq_cha", 1),
                        hwb_precip = define_output("hru_wb_aa", "precip", nb_hru),
                        hwb_snofall = define_output("hru_wb_aa", "snofall", nb_hru),
                        hwb_snomlt = define_output("hru_wb_aa", "snomlt", nb_hru),
                        hwb_wateryld = define_output("hru_wb_aa", "wateryld", nb_hru),
                        hwb_perc = define_output("hru_wb_aa", "perc",nb_hru),
                        hwb_et = define_output("hru_wb_aa", "et", nb_hru),
                        hwb_ecanopy = define_output("hru_wb_aa", "ecanopy", nb_hru),
                        hwb_eplant = define_output("hru_wb_aa", "eplant", nb_hru),
                        hwb_esoil = define_output("hru_wb_aa", "esoil", nb_hru),
                        hwb_cn = define_output("hru_wb_aa", "cn", nb_hru),
                        hwb_sw_ave = define_output("hru_wb_aa", "sw_ave", nb_hru),
                        hwb_sw_300 = define_output("hru_wb_aa", "sw_300", nb_hru),
                        hwb_snopack = define_output("hru_wb_aa", "snopack", nb_hru),
                        hwb_pet = define_output("hru_wb_aa", "pet", nb_hru),
                        hwb_qtile = define_output("hru_wb_aa", "qtile", nb_hru),
                        hwb_surq_cha = define_output("hru_wb_aa", "surq_cha", nb_hru),
                        hwb_latq_cha = define_output("hru_wb_aa", "latq_cha", nb_hru))