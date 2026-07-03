remove_identical_grouped_cols <- function(ssm_df, col_prefix = "^group", verbose = FALSE) {
  
  # columns to check
  grp_cols <- grep(col_prefix, names(ssm_df), value = TRUE)
  
  if (length(grp_cols) < 2) {
    message("Fewer than 2 matching columns found.")
    return(ssm_df)
  }
  
  parse_lambda <- function(x) {
    s <- sub("^group_lambda_", "", x)
    s <- gsub("_", ".", s)
    as.numeric(s)
  }
  
  lambdas <- setNames(vapply(grp_cols, parse_lambda, numeric(1)), grp_cols)
  
  # find sets of identical columns
  visited <- setNames(rep(FALSE, length(grp_cols)), grp_cols)
  identical_sets <- list()
  
  for (col in grp_cols) {
    if (visited[col]) next
    
    same <- grp_cols[
      vapply(grp_cols, function(x) identical(ssm_df[[col]], ssm_df[[x]]), logical(1))
    ]
    
    visited[same] <- TRUE
    
    if (length(same) > 1) {
      identical_sets[[length(identical_sets) + 1]] <- same
    }
  }
  
  if (length(identical_sets) == 0) {
    if (verbose){
      message("No identical group columns found. In: ", paste(unique(ssm_df$Hugo_Symbol), collapse = ", "))
    }  
    return(ssm_df)
  }
  
  cols_to_drop <- character()
  
  for (same in identical_sets) {
    vals <- lambdas[same]
    
    # keep 1 if present; otherwise keep largest absolute lambda
    keep <- if (any(vals == 1, na.rm = TRUE)) {
      same[which(vals == 1)[1]]
    } else {
      same[which.max(abs(vals))]
    }
    
    drop <- setdiff(same, keep)
    
    if (verbose){
      message(
        "In ", paste(unique(ssm_df$Hugo_Symbol), collapse = ", "),":\n",
        "Identical columns found: ",
        paste(same, collapse = ", "),
        "\nKeeping: ", keep,
        "\nRemoving: ", paste(drop, collapse = ", "),
        "\n"
      )
    }
    cols_to_drop <- c(cols_to_drop, drop)
  }
  
  ssm_df[, setdiff(names(ssm_df), cols_to_drop), drop = FALSE]
  
}