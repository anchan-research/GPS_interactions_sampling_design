# 2.5 Infectious period 

# packages -------
library(dplyr)
library(lubridate)
library(ggplot2)
library(amt)

utm_zone <- "+proj=utm +zone=12 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

# read in the clean data from the ctmm fit processed
df_30 <- read.csv("./path/to/file.csv") %>%
  dplyr::select(animal_id, t_, X_transformed, Y_transformed) %>%
  rename(x_ = 3, y_ = 4)

ids <- unique(df_30$animal_id)

# resample using for 1hr, 2hr, 4hr, 12h, 24h, 168hr (1 week - old VHF data)
sampling <- c(30, 1, 2, 4, 12, 24, 168)
resampled_df <- list()

for(i in seq_along(sampling)){
  if(i == 1){
    resampled_df[[i]] <- df_30 %>%
      nest(trk = -animal_id) %>%
      mutate(trk = lapply(trk, function(x) {
        x %>%
          make_track(.x = UTME,
                     .y = UTMN, 
                     .t = t_,
                     all_cols = TRUE) %>%
          track_resample(rate = minutes(sampling[i]),
                         tolerance = minutes(5))
      })) %>%
      unnest(cols = trk) %>%
      mutate(burst_ = paste0(animal_id, "_", burst_)) %>%
      dplyr::select(animal_id, t_, y_, x_) %>%
      rename(ID = animal_id, t = t_, x = x_, y = y_) %>%
      distinct()
  }else{
    resampled_df[[i]] <- df_30 %>%
      nest(trk = -animal_id) %>%
      mutate(trk = lapply(trk, function(x) {
        x %>%
          make_track(.x = UTME,
                     .y = UTMN, 
                     .t = t_,
                     all_cols = TRUE) %>%
          track_resample(rate = hours(sampling[i]),
                         tolerance = minutes(15))
      })) %>%
      unnest(cols = trk) %>%
      mutate(burst_ = paste0(animal_id, "_", burst_)) %>%
      dplyr::select(animal_id, t_, y_, x_) %>%
      rename(ID = animal_id, t = t_, x = x_, y = y_) %>%
      distinct()
  }
}

# create the pairwise dfs ------


sampling_save <- c("30mins", "1hr", "2hr", "4hr", "12hr", "24hr", "168hr")

for(t in seq_along(sampling_save)){
  tt <- resampled_df[[t]] %>%
    mutate(site = stringr::str_split_i(ID, pattern = "_", i = 2))
  
  sites <- unique(tt$site)
  out <- data.frame()
  
  for(k in seq_along(sites)){
    tmp <- tt %>%
      filter(site == sites[[k]])
    
    df_list <- list()
    id_tmp <- unique(tmp$ID)
    for(r in seq_along(id_tmp)){
      df_list[[r]] <- tmp %>%
        filter(ID == id_tmp[r])
    }
    
    
    for(i in 1:(length(df_list) - 1)){
      start_ind <- i + 1
      for(j in start_ind:length(df_list)){
        tmp_out <- dyad_dataframe(
          data_1 = df_list[[i]],
          data_2 = df_list[[j]])
        
        out <- rbind(out, tmp_out)
      }
    }
  }
  
  out <- out %>%
    dplyr::select(-site.y) %>%
    rename(site = site.x)
  
  # save
  saveRDS(out,
          file = paste0("./infectious_period/",
                        sampling_save[[t]],"_pairwise_df.rds"))
  
}

# segment through different infectious period (non-overlapping) ------
# 2-day, 7-day, 14-day, 1-month, 6-month
library(slider)

files <- list.files(path = "./infectious_period/",
                    pattern = "df.rds",
                    full.names = FALSE)

data_list <- lapply(paste0("./infectious_period/",
                           files),
                    readRDS)

# blocking function
infectious_period_block <- function(data, period, window_size, out_type="column"){
  require(slider); require(dplyr)
  period_id <-  block(x = data$index, 
                      i = data$t, 
                      period = period, 
                      every = window_size, 
                      origin = min(data$t))
  out_df <- data.frame()
  for(i in seq_along(period_id)){
    tmp_df <- data.frame(index = period_id[[i]])
    tmp_df$group <- i
    
    out_df <- rbind(out_df, tmp_df)
  }
  
  result <- data %>%
    left_join(out_df, by = "index")
  
  if(out_type == "column"){
    group <- result$group
    return(group)
  }else{
    return(result)
  }
}

# create an empty list to store the result
group_list <- list()

# loop through the items of the data_list
# corresponding to 1week, 1hr, 24hr, 2hr, 30mins, 4hr sampling rate:
for(i in seq_along(data_list)){
  dframe <- data_list[[i]] %>%
    tidyr::drop_na(distance) %>%
    group_by(dyad_id) %>%
    mutate(index = seq(1, n(), 1)) %>%
    ungroup()
  
  group_df <- dframe %>%
    group_by(dyad_id) %>%
    nest_by() %>%
    ungroup() %>%
    mutate(day_2 = purrr::map(data, function(x) infectious_period_block(x,
                                                                        period = "day",
                                                                        window_size = 2,
                                                                        out_type = "data_frame")),
           day_7 = purrr::map(data, function(x) infectious_period_block(x,
                                                                        period = "day",
                                                                        window_size = 7)),
           day_14 = purrr::map(data, function(x) infectious_period_block(x,
                                                                         period = "day",
                                                                         window_size = 14)),
           day_30 = purrr::map(data, function(x) infectious_period_block(x,
                                                                         period = "day",
                                                                         window_size = 30)),
           day_183 = purrr::map(data, function(x) infectious_period_block(x,
                                                                          period = "day",
                                                                          window_size = 183))) %>%
    select(-data) %>%
    tidyr::unnest(c(day_2, day_7, day_14, day_30, day_183)) %>%
    data.table::data.table() %>%
    rename(day_2 = group)
  
  group_list[[i]] <- group_df
}

saveRDS(group_list, 
        file = "./infectious_period/non_overlapping_list.rds")

# define different threshold distance
group_list <- readRDS("./infectious_period/non_overlapping_list.rds")

# interaction definition -----

interaction_only_def_ip <- function(x, dist_scale){
  require(dplyr)
  
  # definition for distance scale
  scale_cutoff <- dist_scale
  
  # data frame
  dframe <- x %>%
    # define contact 
    mutate(
      fission_fusion = case_when(
        # fission events are 0s
        distance > scale_cutoff ~ 0, 
        # fusion events are 1s
        distance <= scale_cutoff ~ 1)
    ) 
  
  # summary for each infectious period
  day_2 <- dframe %>%
    group_by(dyad_id, day_2) %>%
    summarise(N = n(),
              num_contact = sum(fission_fusion),
              p = num_contact/N) %>%
    filter(num_contact != 0) %>%
    rename(week_id = 2) %>%
    mutate(infect_period = "2day") %>%
    suppressMessages()
  
  day_7 <- dframe %>%
    group_by(dyad_id, day_7) %>%
    summarise(N = n(),
              num_contact = sum(fission_fusion),
              p = num_contact/N) %>%
    filter(num_contact != 0) %>%
    rename(week_id = 2) %>%
    mutate(infect_period = "7day") %>%
    suppressMessages()
  
  day_14 <- dframe %>%
    group_by(dyad_id, day_14) %>%
    summarise(N = n(),
              num_contact = sum(fission_fusion),
              p = num_contact/N) %>%
    filter(num_contact != 0) %>%
    rename(week_id = 2) %>%
    mutate(infect_period = "14day") %>%
    suppressMessages()
  
  day_30 <- dframe %>%
    group_by(dyad_id, day_30) %>%
    summarise(N = n(),
              num_contact = sum(fission_fusion),
              p = num_contact/N) %>%
    filter(num_contact != 0) %>%
    rename(week_id = 2) %>%
    mutate(infect_period = "30day") %>%
    suppressMessages()
  
  day_183 <- dframe %>%
    group_by(dyad_id, day_183) %>%
    summarise(N = n(),
              num_contact = sum(fission_fusion),
              p = num_contact/N) %>%
    filter(num_contact != 0) %>%
    rename(week_id = 2) %>%
    mutate(infect_period = "183day") %>%
    suppressMessages()
  
  out <- rbind(day_2, day_7, day_14, day_30, day_183) %>%
    mutate(contact_dist = scale_cutoff)
  
  return(out)
}

interaction_only_def_ip_1week <- function(x, dist_scale){
  require(dplyr)
  
  # definition for distance scale
  scale_cutoff <- dist_scale
  
  # data frame
  dframe <- x %>%
    # define contact 
    mutate(
      fission_fusion = case_when(
        # fission events are 0s
        distance > scale_cutoff ~ 0, 
        # fusion events are 1s
        distance <= scale_cutoff ~ 1)
    ) 
  
  # summary for each infectious period
  day_14 <- dframe %>%
    group_by(dyad_id, day_14) %>%
    summarise(N = n(),
              num_contact = sum(fission_fusion),
              p = num_contact/N) %>%
    filter(num_contact != 0) %>%
    rename(week_id = 2) %>%
    mutate(infect_period = "14day") %>%
    suppressMessages()
  
  day_30 <- dframe %>%
    group_by(dyad_id, day_30) %>%
    summarise(N = n(),
              num_contact = sum(fission_fusion),
              p = num_contact/N) %>%
    filter(num_contact != 0) %>%
    rename(week_id = 2) %>%
    mutate(infect_period = "30day") %>%
    suppressMessages()
  
  day_183 <- dframe %>%
    group_by(dyad_id, day_183) %>%
    summarise(N = n(),
              num_contact = sum(fission_fusion),
              p = num_contact/N) %>%
    filter(num_contact != 0) %>%
    rename(week_id = 2) %>%
    mutate(infect_period = "183day") %>%
    suppressMessages()
  
  out <- rbind(day_14, day_30, day_183) %>%
    mutate(contact_dist = scale_cutoff)
  
  return(out)
}

# run through each fix rate
# do not run the 1 week fix rate for 2 and 7 day infectious period

dist_scale <- seq(20, 300, 20)
fix_rate <- c("12hr", "1week" ,"1hr", "24hr", "2hr", "30mins", "4hr")

contact_list <- list()

for(i in seq_along(group_list)){
  if(i == 2){
    contact_list[[i]] <- do.call(rbind, 
                                 lapply(dist_scale, function(x) interaction_only_def_ip_1week(x = group_list[[i]],
                                                                                          dist_scale = x))) %>%
      mutate(fix_rate = fix_rate[[i]])
  }else{
    contact_list[[i]] <- do.call(rbind, 
                                 lapply(dist_scale, function(x) interaction_only_def_ip(x = group_list[[i]],
                                                                                    dist_scale = x))) %>%
      mutate(fix_rate = fix_rate[[i]])
  }
}

# calculate variance
# for each dyad, each inf period, each fix rate, each interaction_dist
# assuming binomial distribution

dat <- plt_data %>%
  dplyr::group_by(infect_period, dyad_id, interaction_dist, fix_rate) %>%
  summarise(N = sum(N),
            contact_sum = sum(num_contact),
            p_contact = contact_sum/N,
            mean = (contact_sum * (1 - p_contact))/p_contact,
            variance = (contact_sum*(1 - p_contact))/p_contact^2,
            cv = sqrt(variance)/mean) %>%
  ungroup()

# check the data
hist(dat$variance, breaks = 100)
hist(dat$cv, breaks = 100)

# export -------
saveRDS(dat,
        file = "./infectious_period/prep_data_contact_only.rds")