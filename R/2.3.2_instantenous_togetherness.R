#2.2 Instantaneous Togetherness
instan_togetherness <- function(x, fix_rate, dist_scale){
  require(dplyr)
  
  scale_cutoff <- dist_scale
  id_1 <- unique(x$ID_1)
  id_2 <- unique(x$ID_2)
  
  # data frame
  togetherness <- x %>%
    # define contact 
    mutate(
      fission_fusion = case_when(
        # fission events are 0s
        distance > scale_cutoff ~ 0, 
        # fusion events are 1s
        distance <= scale_cutoff ~ 1),
    ) %>%
    group_by(t) %>%
    summarise(sim_in_contact = sum(fission_fusion))
  
  togetherness$dyad_pair <- paste0(id_1, "_", id_2)
  togetherness$interaction_dist <- scale_cutoff
  togetherness$fix_rate <- fix_rate
  
  return(togetherness)
}
