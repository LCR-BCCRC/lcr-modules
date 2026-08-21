hierachical_grouping = function(ssm_df, lambda = 1){
  maf = ssm_df %>% dplyr::ungroup() %>% arrange(Start_Position)
  positions = unique(maf$Start_Position)
  
  dist_mat = dist(positions, method = "euclidean")
  hc_a = hclust(dist_mat, method = "centroid")
  
  dist_matrix <- as.matrix(dist_mat)
  
  h_values <- 1:100
  avg_dist <- numeric(length(h_values))
  
  for (i in seq_along(h_values)) {
    h_val <- h_values[i]
    
    cut_lab <- cutree(hc_a, h = h_val)
    clusters <- unique(cut_lab)
    h_output = data.frame(position = positions, group = cut_lab)
    h_range = h_output %>% group_by(group) %>% summarise(size = (max(position)-min(position)+1)) %>% pull(size)
    # If all values are in one cluster → skip
    if (length(unique(cut_lab)) > 1) {
      mean_pw_dist <- numeric(length(clusters))
      names(mean_pw_dist) <- clusters
      
      # Loop through clusters
      for (cl in clusters) {
        idx <- which(cut_lab == cl)
        k = length(idx)
        
        
        if (length(idx) > 1) {
          sub_dists <- dist_matrix[idx, idx]
          mean_pw_dist[as.character(cl)] <-  (k/h_range[cl]) * (1 / ( 1 + mean(sub_dists[lower.tri(sub_dists)]))) 
        } else {
          mean_pw_dist[as.character(cl)] <- 0  # No pairwise distance possible
        }
      }
      
      
      avg_dist[i] = mean(mean_pw_dist) /  (length(clusters))^log2(lambda*median(diff(positions))) 
      
    } else {
      avg_dist[i] <- NA  # Not meaningful for 1 cluster
    }
  }
  
  best_h = which.max(avg_dist)
  
  plotting_df = df <- data.frame(h = h_values, avg = avg_dist)
  p <- ggplot(df, aes(x = h, y = avg)) +
    geom_line() +
    geom_point(shape = 16) +
    labs(
      x = "Height (h)",
      y = "Score",
      title = "Clustering Quality vs. Tree Cut Height"
    ) +
    theme_minimal()
  
  cut_labels = cutree(hc_a, h = best_h)
  
  h_output = data.frame(position = positions, group = cut_labels)
  
  maf2 = left_join(maf, h_output, by = c("Start_Position" = "position"))
  
  return(list(maf = maf2, plot = p))
}