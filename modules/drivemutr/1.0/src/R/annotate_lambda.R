# Adjust lambda for column name
lambda_suffix <- function(x) {
  sapply(x, function(v) {
    if (v == 0) return("lambda_0")
    
    val_part <- gsub("\\.", "_", as.character(abs(v))) # Remove . and use _ instead
    
    paste0("lambda_",val_part)
  })
}

annotate_lambda = function(df, def, lambda = 1, add_df = FALSE){
  result <- df
  original_names <- names(result)
  has_lambda <- "lambda" %in% names(formals(def)) # Check if lambda is an argument of function
  list_data <- list()
  tmp_df <- NULL
  if (!has_lambda){
    stop(def, " Doesn't have lambda argument")
  }
  
  for (lam in lambda){
    
    tmp_df <- def(df, lambda = lam)
    
    if(is.data.frame(tmp_df)){
      new_cols = setdiff(names(tmp_df), original_names)
      
      if (length(new_cols)==0){
        next
      } else{
        
        suffix = lambda_suffix(lam)
        renamed_new_cols <- ifelse(
          grepl("_custom$", new_cols),
          sub("_custom$", paste0("_", suffix, "_custom"), new_cols),
          paste0(new_cols, "_", suffix)
        )
        
        names(tmp_df)[match(new_cols, names(tmp_df))] <- renamed_new_cols
        
        result <- suppressMessages(dplyr::left_join(
          result,
          tmp_df
        ))
      }
    }
    
    else if(is.list(tmp_df)){
      if (add_df == FALSE){
        maf_new_col <- setdiff(names(tmp_df[["maf"]]), names(df))
        suffix = lambda_suffix(lam)
        renamed_maf_new_cols <- paste0(maf_new_col, "_", suffix)
        names(tmp_df[["maf"]])[match(maf_new_col, names(tmp_df[["maf"]]))] <- renamed_maf_new_cols
        
        result <- suppressMessages(dplyr::left_join(
          result,
          tmp_df[["maf"]]
        ))
        list_data[[suffix]] = tmp_df[["plot"]]
      } else {
        if (length(tmp_df) == 0L) {
          next
        }
        list_data <- c(list_data, tmp_df)
      }
    }
  }
  
  if (is.data.frame(tmp_df)){
    
    return(result)
    
  } else if (is.list(tmp_df)){
    if (add_df == FALSE){
      return(list(maf = result, plot = list_data))
    }
  }
  if (add_df) {
    return(tibble::tibble(
      lambda_name = names(list_data),
      data = list_data
    ))
  }
}