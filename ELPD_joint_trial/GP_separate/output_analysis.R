
library(rlist)
library(this.path)

mydir <- this.path::here()
setwd(mydir)

##### Functions for comparison

compute_summary <- function(results, metric_path) {
  extract_metric <- function(res) {
    parts <- strsplit(metric_path, "\\$")[[1]]
    for (part in parts) {
      if (is.null(res)) return(NULL)
      res <- res[[part]]
    }
    return(res)
  }
  
  # Extract all metrics
  metric_list <- lapply(results, extract_metric)
  metric_list <- Filter(Negate(is.null), metric_list)
  
  # Abort early if empty or non-numeric
  if (length(metric_list) == 0) {
    stop("No valid metric data found for path: ", metric_path)
  }
  
  # Coverage: flatten and compute mean/se
  if (grepl("coverage", metric_path)) {
    flat_vals <- unlist(metric_list)
    suppressWarnings(flat_vals <- as.numeric(flat_vals))
    flat_vals <- flat_vals[!is.na(flat_vals)]
    if (length(flat_vals) == 0) {
      return(data.frame(Mean = NA_real_, SE = NA_real_))
    }
    return(data.frame(
      Mean = mean(flat_vals),
      SE = sd(flat_vals) / sqrt(length(flat_vals))
    ))
  }
  
  # Scalar numeric
  if (all(sapply(metric_list, function(x) is.numeric(x) && length(x) == 1))) {
    values <- unlist(metric_list)
    values <- values[!is.na(values)]
    if (length(values) == 0) {
      return(data.frame(Mean = NA_real_, SE = NA_real_))
    }
    return(data.frame(
      Mean = mean(values),
      SE = sd(values) / sqrt(length(values))
    ))
  }
  
  # Vector numeric
  if (is.numeric(metric_list[[1]]) && is.vector(metric_list[[1]])) {
    values_mat <- try(do.call(cbind, metric_list), silent = TRUE)
    if (inherits(values_mat, "try-error")) {
      stop("Failed to bind metric vectors: inconsistent lengths.")
    }
    mean_value <- rowMeans(values_mat, na.rm = TRUE)
    se_value <- apply(values_mat, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
    return(data.frame(Mean = mean_value, SE = se_value))
  }
  
  # Matrix numeric
  if (is.matrix(metric_list[[1]])) {
    dims <- dim(metric_list[[1]])
    n_rep <- length(metric_list)
    values_array <- array(unlist(metric_list), dim = c(dims[1], dims[2], n_rep))
    mean_matrix <- apply(values_array, c(1, 2), mean, na.rm = TRUE)
    se_matrix <- apply(values_array, c(1, 2), function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
    return(list(Mean = mean_matrix, SE = se_matrix))
  }
  
  stop("Unsupported or inconsistent data structure for metric: ", metric_path)
}


compute_min_ess_summary <- function(results, metric_name) {
  # Get the minimum ESS value from each result for the given metric
  min_ess <- sapply(results, function(res) {
    ess_vals <- res[[metric_name]]
    min(as.numeric(ess_vals))
  })
  
  # Compute mean and SE
  mean_val <- mean(min_ess)
  se_val <- sd(min_ess) / sqrt(length(min_ess))
  
  return(data.frame(
    Metric = metric_name,
    Mean_min_ESS = mean_val,
    SE = se_val
  ))
}







