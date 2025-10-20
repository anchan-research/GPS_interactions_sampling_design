# 2.4 Sample size 

required_packages <- c("dplyr", "lubridate", "ggplot2",
                       "amt")

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

# plotting theme
theme_plt <-  ggplot2::theme(panel.background = element_rect(fill = "white",
                                                             inherit.blank = TRUE),
                             panel.grid.major = element_line(linetype = "solid",
                                                             linewidth = 0.5,
                                                             colour = "#efe9e8",
                                                             inherit.blank = FALSE),
                             plot.background = element_rect(color = "black",
                                                            size = 1),
                             axis.text.x = element_text(size = 20, angle = 0),
                             axis.text.y = element_text(size = 20), 
                             axis.line = element_line(colour = "black"),
                             axis.title.x = element_text(size = 22),
                             axis.title.y = element_text(size = 22),
                             plot.margin = ggplot2::margin(t = 1,  # Top margin
                                                           r = 2,  # Right margin
                                                           b = 1,  # Bottom margin
                                                           l = 2,  # Left margin
                                                           unit = "cm"),
                             title = element_text(size = 20),
                             legend.text = element_text(size = 22),
                             legend.title = element_text(size = 24),
                             strip.text = element_text(size = 20))


## Evaluation -----

# get random draws to process
set.seed(1000)
sample_size <- seq(2, 88, 1)

for(n in seq_along(sample_size)) {
  id_list <- list()
  for(i in 1:100){
    id_list[[i]] <- sample(final_ids, # vector collared animals
                           size = sample_size[[n]],
                           replace = FALSE)
  }
  out <- do.call(rbind, id_list)
  saveRDS(out,
          file = paste0(".sample_size/ids_",
                        sample_size[[n]],
                        ".rds"))
}
# ***************************************************************
# ***************************************************************

# ***************************************************************
# extraction to be run on the workstation
# run the first 9 on here
# ***************************************************************
# ***************************************************************

sample_sizes <- seq(2, 88, 1)
sample_files <- paste0("./sample_size/ids_", sample_sizes, ".rds")

ids_list <- lapply(sample_files,
                   readRDS)

# fix_rates <- c(0.5, 1, 2, 4, 12, 24, 168)
dist_scale <- seq(20, 300, 20)

# leave 2 cores for other system processes
library(doMC)
library(doParallel)
cores_to_use <- detectCores() - 2
registerDoMC(cores_to_use)

summary_data <- foreach(u = seq_along(ids_list)) %dopar% {# sample size loop
  result_df <- data.frame()
  for(k in 1:dim(ids_list[[u]])[1]) { # sampling variation loop
    
    tmp_id <- ids_list[[u]][k,]
    
    df_list <- list()
    for(r in seq_along(tmp_id)){
      df_list[[r]] <- df %>%
        filter(ID == tmp_id[r])
    }
    
    out <- data.frame()
    for(i in 1:(length(df_list) - 1)){
      start_ind <- i + 1
      for(j in start_ind:length(df_list)){
        tmp_out <- dyad_dataframe(
          data_1 = df_list[[i]],
          data_2 = df_list[[j]])
        
        out <- rbind(out, tmp_out)
      }
    }
    
    if(nrow(out)==0){
      next
    }else{
      out <- out %>%
        tidyr::drop_na(distance) %>%
        group_by(dyad_id) %>%
        arrange(t) %>%
        mutate(index = seq(1, n(), 1)) %>%
        ungroup()
      
      out <- out %>%
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
      
      # get the different def of distance and calculate variance in togetherness 
      contact_out <- do.call(rbind, 
                             lapply(dist_scale,
                                    function(x) contact_only_def_ip(x = out,
                                                                    dist_scale = x))) 
      
      # summarize
      summary_df <- contact_out %>%
        dplyr::group_by(dyad_id, infect_period, contact_dist) %>%
        summarise(N = sum(N),
                  contact_sum = sum(num_contact),
                  p_contact = contact_sum/N,
                  mean = (contact_sum * (1 - p_contact))/p_contact,
                  variance = (contact_sum*(1 - p_contact))/p_contact^2,
                  cv = sqrt(variance)/mean) %>%
        ungroup() %>%
        suppressMessages()
      
      # bind the result
      
      result_df <- rbind(result_df, summary_df)
    }
  }
  
  result_df <- result_df %>%
    mutate(sample_size = sample_sizes[[u]])
  
  saveRDS(result_df,
          file = paste0("./sample_size_",
                        sample_sizes[[u]],
                        ".rds"))
  return(result_df)
}

files_summary <- list.files(path = "./sample_size/",
                            pattern = "sample_size_",
                            full.names = TRUE)


summary_data_rc <- lapply(files_summary, readRDS)

# function to extract metric at sample size level
summary_function <- function(x) {
  out <- x %>%
    group_by(sample_size, infect_period, interaction_dist) %>%
    summarise(mean_mean = mean(mean),
              sd_mean = sd(mean),
              cv_mean = sd_mean/mean_mean,
              mean_cv = mean(cv, na.rm = TRUE),
              sd_cv = sd(cv, na.rm = TRUE),
              mean_var = mean(variance)) %>%
    suppressMessages()
  
  return(out)
}

final_out_rc <- do.call(rbind, lapply(
  summary_data_rc, summary_function
)) %>%
  mutate(sample_size = as.factor(sample_size),
         contact_dist = as.factor(interaction_dist))

plt_data_rc <- final_out_rc %>%
  filter(contact_dist == 40) %>%
  mutate(infect_period = case_when(
    infect_period == "2day" ~ "2 days",
    infect_period == "7day" ~ "7 days",
    infect_period == "14day" ~ "14 days",
    infect_period == "30day" ~ "30 days",
    infect_period == "183day" ~ "183 days",
  ),
  infect_period = factor(infect_period,
                         levels = c("2 days", "7 days", "14 days",
                                    "30 days", "183 days")))

# mean of the cv 
ggplot(plt_data_rc) +
  geom_point(aes(x = sample_size,
                 y = mean_cv),
             size = 2,
             alpha = 0.5) +
  facet_wrap(~infect_period, 
             axis.labels = "margins") +
  scale_x_discrete(breaks = seq(0, 45, 10)) +
  theme_plt +
  labs(x = "Sample size", 
       y = "Mean cv of togetherness") +
  theme(legend.position = "bottom")


## Simulation ------
popsize <- 1000
groupsize <- 10

dyad_replication_calculator <- function(proposed_samplesizes,
                                        # target_samegrp_dyads,
                                        groupsize,
                                        popsize,
                                        n_reps){
  
  # figure out number of groups
  # for every animal, pick a group at random
  # check to be sure no group has more than groupsize animals, if yes, repick for somebody
  # table number of animals in each group
  # replace table values by tabled value choose 2 (this accounts for multiple dyads within the same group)
  # sum those numbers up over all groups
  
  n_groups <- popsize / groupsize
  dyads_within_grps <- rep(NA, n_reps)
  for(i in 1:n_reps){
    group_assignments <- sample(x = 1:n_groups, size = proposed_samplesizes, replace = T)
    counts_within_grps <- table(group_assignments)
    dyads_within_grps[i] <- sum(choose(counts_within_grps, 2))
  }
  
  return(dyads_within_grps)
}

dyad_replication_calculator <- function(proposed_samplesizes,
                                        # target_samegrp_dyads,
                                        groupsize,
                                        popsize,
                                        n_reps){
  
  # figure out number of groups
  # for every animal, pick a group at random
  # check to be sure no group has more than groupsize animals, if yes, repick for somebody
  # table number of animals in each group
  # replace table values by tabled value choose 2 (this accounts for multiple dyads within the same group)
  # sum those numbers up over all groups
  
  simmed_group_sizes <- rpois(n = 100, lambda = groupsize)
  group_number <- seq(1:100)
  group_labels <- rep(group_number, each = simmed_group_sizes)[1:popsize]
  #  cumsum(simmed_group_sizes)
  #  n_groups <- popsize / groupsize
  dyads_within_grps <- rep(NA, n_reps)
  for(i in 1:n_reps){
    group_assignments <- sample(x = group_labels, size = proposed_samplesizes, replace = F)
    counts_within_grps <- table(group_assignments)
    dyads_within_grps[i] <- sum(choose(counts_within_grps, 2))
  }
  
  return(dyads_within_grps)
}



dyad_replication_calculator(proposed_samplesizes = 5,
                            # target_samegrp_dyads = 10,
                            groupsize = 10,
                            popsize = 1000,
                            n_reps = 100)

proposed_ss <- seq(5:200)
n_reps_in <- 10
dyads_within_grps <- dyads_within_grps_gs2 <- dyads_within_grps_gs50 <- proposed_ss_mat <- matrix(NA, nrow = length(proposed_ss), ncol = n_reps_in)
for(i in 1:length(proposed_ss)){
  dyads_within_grps[i, ] <- dyad_replication_calculator(proposed_samplesizes = proposed_ss[i],
                                                        # target_samegrp_dyads = 10,
                                                        groupsize = 10,
                                                        popsize = 1000,
                                                        n_reps = n_reps_in)
  dyads_within_grps_gs2[i, ] <- dyad_replication_calculator(proposed_samplesizes = proposed_ss[i],
                                                            # target_samegrp_dyads = 10,
                                                            groupsize = 2,
                                                            popsize = 1000,
                                                            n_reps = n_reps_in)
  dyads_within_grps_gs50[i, ] <- dyad_replication_calculator(proposed_samplesizes = proposed_ss[i],
                                                             # target_samegrp_dyads = 10,
                                                             groupsize = 50,
                                                             popsize = 1000,
                                                             n_reps = n_reps_in)
  proposed_ss_mat[i, ] <- rep(proposed_ss[i], n_reps_in)
}

dyad_count_vec <- as.vector(dyads_within_grps)
dyad_count_vec_gs2 <- as.vector(dyads_within_grps_gs2)
dyad_count_vec_gs50 <- as.vector(dyads_within_grps_gs50)
ss_vec <- as.vector(proposed_ss_mat)


median_dyads <- apply(dyads_within_grps, 1, median)
median_dyads_gs2 <- apply(dyads_within_grps_gs2, 1, median)
median_dyads_gs50 <- apply(dyads_within_grps_gs50, 1, median)
# ss where median is closest to 10:
closest_to_10 <- which.min(abs(median_dyads - 10))
closest_to_10_gs2 <- which.min(abs(median_dyads_gs2 - 10))
closest_to_10_gs50 <- which.min(abs(median_dyads_gs50 - 10))
#proposed_ss[closest_to_10]

closest_to_20 <- which.min(abs(median_dyads - 20))
closest_to_20_gs2 <- which.min(abs(median_dyads_gs2 - 20))
closest_to_20_gs50 <- which.min(abs(median_dyads_gs50 - 20))
#proposed_ss[closest_to_20]

pdf("CombinatoricsGroupThing_v2.pdf", height = 6, width = 10)
par(mfrow = c(1, 3), mar = c(4, 6, 2, 2))
plot(dyad_count_vec_gs2 ~ ss_vec, xlab = "Animals sampled", ylim = c(1, 1300),
     log = "y", 
     ylab = "Number of tracked dyads\nwithin the same group", las = 1)
abline(h = 10, lty = 2, col = "grey80")
abline(v = proposed_ss[closest_to_10_gs2])
abline(h = 20, lty = 2, col = "grey80")
abline(v = proposed_ss[closest_to_20_gs2])
title("Groups of 2, Population of 1000")

plot(dyad_count_vec ~ ss_vec, xlab = "Animals sampled", ylim = c(1, 1300),
     log = "y",
     ylab = "Number of tracked dyads\nwithin the same group", las = 1)
abline(h = 10, lty = 2, col = "grey80")
abline(v = proposed_ss[closest_to_10])
abline(h = 20, lty = 2, col = "grey80")
abline(v = proposed_ss[closest_to_20])
title("Groups of 10, Population of 1000")

plot(dyad_count_vec_gs50 ~ ss_vec, xlab = "Animals sampled", ylim = c(1, 1300),
     log = "y",
     ylab = "Number of tracked dyads\nwithin the same group", las = 1)
abline(h = 10, lty = 2, col = "grey80")
abline(v = proposed_ss[closest_to_10_gs50])
abline(h = 20, lty = 2, col = "grey80")
abline(v = proposed_ss[closest_to_20_gs50])
title("Groups of 50, Population of 1000")
dev.off()