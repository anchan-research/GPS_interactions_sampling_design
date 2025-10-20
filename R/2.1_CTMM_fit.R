# fitting IOU

required_packages <- c("lubridate", "dplyr", 
                       "ctmm", "doMC", "doParallel")

for (i in seq_along(required_packages)){
  p <- installed.packages()[,1]
  if(required_packages[i] %in% p){
    library(required_packages[i], character.only = T)
  }else{
    install.packages(required_packages[i])
    library(required_packages[i], character.only = T)
  } 
}

rm(i,p);gc()

# ******************************************************************************
# ******************************************************************************

# leave 2 cores for other system processes
cores_to_use <- detectCores() - 2
registerDoMC(cores_to_use)

# utm_zone
utm_zone <- "+proj=utm +zone=12 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

# ******************************************************************************
# ******************************************************************************

df <- read.csv("/path/to/file.csv") %>%
  dplyr::select(animal_id, t_, Latitude, Longitude) %>%
  rename(ID = 1, Date = 2, location.lat = 3, location.long = 4)


# ******************************************************************************
# ******************************************************************************
# Prep ----- 

# fit telemetry and semi-variance variogram
Unique_ID_res <- unique(df$ID)

# Structuring the data and creating telemetry object
ls_tele_res <- foreach(i = seq_along(Unique_ID_res)) %dopar% {
  
  # subset the data
  tt <- df %>%
    filter(ID == Unique_ID_res[i]) %>%
    filter(!duplicated(Date)) %>%
    arrange(Date) %>%
    droplevels()
  
  # create a telemetry object with correct utm and save:
  tele_tmp <- as.telemetry(tt,
                           timezone = "US/Mountain",
                           projection = utm_zone)
  
  # semi-variogram
  SVF <- variogram(tele_tmp)
  
  # package the data for output
  list(telemetry = tele_tmp, SVF = SVF)
}

# fit model
ctmm_res <- foreach(i = seq_along(ls_tele_res)) %dopar% {
  
  # fitting IOU model without following the recommended framework
  # but only taking the initial value from the guess for velocity 
  # and fixing the position autocorrelation to Inf using ctmm.select
  guess <- ctmm.guess(ls_tele_res[[i]]$telemetry,
                      interactive = FALSE)
  
  # ctmm fitting to IOU
  ctmm_obj <- ctmm.fit(ls_tele_res[[i]]$telemetry,
                       ctmm(tau = c(Inf, summary(guess)$CI[3, 2])))
  
  # package the data for download
  list(ctmm_object = ctmm_obj)
}

# ******************************************************************************
# ******************************************************************************
# Simulation -----

# initiate the index
sim_num <- seq(1, 20, 1)

start_time <- Sys.time()
for(i in seq_along(Unique_ID_res)){# loop through each individuals
  # subset the data
  id <- Unique_ID_res[[i]]
  
  tmp <- df %>%
    dplyr::filter(ID == id)
  
  # get the time for intervals
  min_t <- min(tmp$Date)
  max_t <- max(tmp$Date)
  interp_vals <- seq(min_t, max_t, by=60)
  ind <- (interp_vals >= min_t) & (interp_vals <= max_t)
  new_t <- interp_vals[ind]
  
  sim <- foreach(j = seq_along(sim_num)) %dopar% { # simulation loop
    traj_tmp <- ctmm::simulate(object = ls_tele_res[[i]]$telemetry, 
                               CTMM = ctmm_res[[i]]$ctmm_object,
                               t = new_t,
                               nsim = 1)
    df <- traj_tmp[,] %>%
      data.frame() %>%
      mutate(sim_num = j,
             t = new_t)
    # package data
    list(simulation = df)
  } # end of parallel simulation
  
  sim_df <- NULL
  
  for(k in seq_along(sim)){ # post simulation processing
    tmp <- sim[[k]]$simulation
    sim_df <- rbind(sim_df, tmp)
  } 
  
  # add usu_id
  sim_df$ID <- id
  
  # save the data
  saveRDS(sim_df,
          paste0("./processed_data/simulation/sim_traj/",
                 fix_rate,"/", id, ".rds"))
  
}# end of individual loop
