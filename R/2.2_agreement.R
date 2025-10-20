# Agreement

calculate_distance <- function(x1, y1, x2, y2) {
  distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  return(distance)
}

# creating dyadic dataframe
dyad_dataframe <- function(data_1, data_2, file_path, save){
  
  sim_index <- unique(data_1$sim_num)
  
  id_1 <- unique(data_1$ID)
  id_2 <- unique(data_2$ID)
  
  time_frame <- c(pmax(min(data_1$t),min(data_2$t)),
                  pmin(max(data_1$t), max(data_2$t))
  )
  
  data_1 <- data_1 %>%
    filter(t >= time_frame[[1]] & t <= time_frame[[2]])
  
  data_2 <- data_2 %>%
    filter(t >= time_frame[[1]] & t <= time_frame[[2]])
  # this will have 20 list in the there
  contact <- list()
  
  for(i in seq_along(sim_index)){
    # go through each simulation of first ind. with the rest of the simulation
    df_1 <- data_1 %>%
      filter(sim_num == sim_index[[i]]) %>%
      rename(x_1 = x, y_1 = y, ID_1 = ID) %>%
      dplyr::select(t, x_1, y_1, sim_num, ID_1)
    
    # add for each loop 
    tmp_list <- foreach(j = seq_along(sim_index)) %dopar% {
      
      df_2 <- data_2 %>%
        filter(sim_num == sim_index[[j]]) %>%
        rename(ID_2 = ID, x_2 = x, y_2 = y, 
               sim_num_2 = sim_num) %>%
        dplyr::select(t, x_2, y_2, sim_num_2, ID_2)
      
      results <- df_1 %>%
        left_join(df_2, by = "t") %>%
        mutate(
          distance = calculate_distance(
            x1 = x_1,
            y1 = y_1,
            x2 = x_2,
            y2 = y_2),
          dyad_id = paste0(sim_num, "_", sim_num_2)
        )
      return(results)
    }
    contact[[i]] <- do.call(rbind, tmp_list)
  }
  
  output <- do.call(rbind, contact) %>%
    dplyr::select(-x_1, -x_2, -y_1, -y_2) %>%
    mutate(date_col = date(t))
  
  # save the output
  if(save){
    saveRDS(output, 
            paste0(file_path,"contact_", id_1, "&", id_2, ".rds"))
  }
  
  return(output)
}

# agreement summary
agreement_metric <- function(x, fix_rate, dist_scale){
  
  scale_cutoff <- dist_scale
  
  rate <- fix_rate
  
  # data frame
  agreement <- x %>%
    # define interaction distance 
    mutate(
      fission_fusion = case_when(
        # fission interaction events are 0s
        distance > scale_cutoff ~ 0, 
        # fusion interaction events are 1s
        distance <= scale_cutoff ~ 1),
    ) %>%
    group_by(t) %>%
    # sum conditional or not 
    summarise(sim_in_contact = sum(fission_fusion),
              sim_not_contact = 400 - sim_in_contact) %>%
    mutate(agreement = pmax(sim_in_contact, sim_not_contact)) 
  
  
  # for each pair calculate proportion of agreement 
  out_df <- data.frame(prop_agreement = sum(agreement$agreement)/(nrow(agreement)*400),
                       agreement = sum(agreement$agreement),
                       togetherness = sum(agreement$sim_in_contact),
                       N_agreement = nrow(agreement) * 400,
                       dyad_pair = paste0(unique(x$ID_1),"_",unique(x$ID_2)),
                       ID_1 = unique(x$ID_1),
                       ID_2 = unique(x$ID_2),
                       contact_dist = dist_scale,
                       fix_rate = rate)
  return(out_df)
}





