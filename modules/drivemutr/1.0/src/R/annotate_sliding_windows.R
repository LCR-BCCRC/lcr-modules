split_region_into_bins = function(ssm_df_sub,
                                  group_column,
                                  window_size    = 20,
                                  max_range      = 50,
                                  slide_by       = 1,
                                  min_mutations_count = 6){
  
  # Sort by position (so sliding makes sense)
  maf <- ssm_df_sub %>% arrange(Start_Position)
  
  # Prepare a group vector: one value per row in maf
  group_num <- integer(nrow(maf))
  
  # Track which rows have already been assigned to a group
  assigned  <- logical(nrow(maf))  # all FALSE initially
  
  current_group   <- 1
  start_pos_group <- maf$Start_Position[1]
  total_window    <- integer(0)    # holds row indices not yet finalized
  
  for (pos in seq(min(maf$Start_Position), max(maf$Start_Position), by = slide_by)) {
    window_start <- pos
    window_end   <- pos + window_size
    span         <- pos - start_pos_group
    
    # Find rows in [window_start, window_end] that have NOT yet been assigned
    in_window <- which(
      !assigned & 
        maf$Start_Position >= window_start & 
        maf$Start_Position <= window_end
    )
    
    # Accumulate them into our current window
    total_window <- union(total_window, in_window)
    
    # If we exceed the number-of-mutations or distance threshold,
    # finalize this current window as one group and move on
    if (length(total_window) >= min_mutations_count && span > max_range) {
      # Assign all those rows to current_group
      group_num[total_window] <- current_group
      
      # Mark those rows as assigned so they're not reused in the future
      assigned[total_window] <- TRUE
      
      # Advance group
      current_group   <- current_group + 1
      start_pos_group <- pos + 1
      
      # Reset the window
      total_window <- integer(0)
    }
  }
  
  # After the loop ends, if there are rows still in total_window,
  # assign them to the final group.
  if (length(total_window) > 0) {
    group_num[total_window] <- current_group
    assigned[total_window]  <- TRUE
  }
  
  maf[[group_column]]= paste0(maf[[group_column]][1], "_", group_num)
  return(maf)
}  





split_region_custom_col = function(ssm_df,
                                   window_size         = 20,
                                   max_range           = 50,
                                   slide_by            = 1,
                                   min_mutations_count = 6,
                                   min_group_range     = 40) {

  ssm_df$row_id_tmp = 1:nrow(ssm_df)
  
  for (col in colnames(ssm_df)[startsWith(colnames(ssm_df), "group")]){
    
    ssm_df = ssm_df %>%
      ungroup() %>%
      group_by(!!rlang::sym(col)) %>%
      mutate(length_group = max(Start_Position) - min(Start_Position)) %>%
      ungroup()
    
    if (any(ssm_df$length_group >= min_group_range)){
      col_custom = paste0(col, "_custom")
      
      ssm_df = ssm_df %>%
        mutate(!!col_custom := as.character(!!rlang::sym(col)))
      
      groups_to_split = unique(
        ssm_df %>%
          filter(length_group >= min_group_range) %>%
          pull(!!rlang::sym(col))
      )
      
      for (grp in groups_to_split){
        ssm_df_sub = ssm_df %>%
          filter(!!rlang::sym(col) == grp)
        
        ssm_df_sub = split_region_into_bins(
          ssm_df_sub = ssm_df_sub,
          group_column = col_custom,
          window_size         = window_size,
          max_range           = max_range,
          slide_by            = slide_by,
          min_mutations_count = min_mutations_count
        )
        
        ssm_df[match(ssm_df_sub$row_id_tmp, ssm_df$row_id_tmp), col_custom] = ssm_df_sub[[col_custom]]
      }
    }
    
    ssm_df = ssm_df %>% dplyr::select(-length_group)
  }
  
  ssm_df = ssm_df %>% dplyr::select(-row_id_tmp)
  
  return(ssm_df)
}