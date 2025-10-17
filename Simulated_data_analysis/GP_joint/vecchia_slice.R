max_min <- function(locations) {
  n <- nrow(locations)
  m <- round(sqrt(n))  # Adaptive neighborhood size
  
  # Compute nearest neighbors
  NNall <- get.knn(locations, k = m)$nn.index
  
  # Find the first point closest to the centroid
  centroid <- colMeans(locations)
  dist_centroid <- sqrt(rowSums((locations - matrix(centroid, nrow = n, ncol = ncol(locations), byrow = TRUE))^2))
  first_point <- which.min(dist_centroid)
  
  # Initialize order and tracking
  ordered_indices <- integer(n)
  ordered_indices[1] <- first_point
  selected <- rep(FALSE, n)
  selected[first_point] <- TRUE
  
  # Store minimum distances efficiently
  min_dists <- rep(Inf, n)
  for (i in seq_len(n)) {
    min_dists[i] <- sqrt(sum((locations[i, ] - locations[first_point, ])^2))
  }
  
  # Max-Min Ordering with Dynamic Neighborhood
  for (k in 2:n) {
    # Find the next point as the one farthest from selected points
    next_point <- which.max(min_dists)
    ordered_indices[k] <- next_point
    selected[next_point] <- TRUE
    min_dists[next_point] <- -Inf  # Prevent re-selection
    
    # Dynamically update neighborhood size
    nneigh <- round(min(m, n / k))
    neighbors <- NNall[next_point, 1:nneigh]
    
    # Update min distances for all unselected points
    for (i in which(!selected)) {
      new_dist <- sqrt(sum((locations[i, ] - locations[next_point, ])^2))
      min_dists[i] <- min(min_dists[i], new_dist)
    }
  }
  
  return(ordered_indices)
}


# Neighbor matrix function
neighbor_matrix <- function(locs_ord, m) {
  n <- nrow(locs_ord)
  nn_matrix <- matrix(NA, nrow = n, ncol = m + 1)
  
    # Sequential computation using for loop
    for (i in 1:n) {
      succ_ind <- i:n
      dists <- sqrt(colSums((t(locs_ord[succ_ind, , drop = FALSE]) - locs_ord[i, ])^2))
      nearest_indices_within_subset <- order(dists)
      sorted_succ_ind <- succ_ind[nearest_indices_within_subset]
      sorted_neighbors <- sorted_succ_ind[1:min(m + 1, length(sorted_succ_ind))]
      row <- rep(NA, m + 1)
      row[1:length(sorted_neighbors)] <- sort(sorted_neighbors)
      nn_matrix[i, ] <- row
    }
  
  row.names(nn_matrix) <- NULL
  return(nn_matrix)
}



# Distance computation function returning a list of n matrices
dist.nn <- function(locs.ord, neighbor_matrix) {
  
  n <- nrow(locs.ord)
  
    dist_nn <- list()
    for (i in 1:n) {
      neighbors <- neighbor_matrix[i, ]
      neighbors <- neighbors[!is.na(neighbors)]  # Remove NA values
      
      if (length(neighbors) > 0) {
        dist_nn[[i]] <- rdist(locs.ord[neighbors, , drop = FALSE], locs.ord[neighbors, , drop = FALSE])
      } else {
        dist_nn[[i]] <- NULL
      }
    }
  
  return(dist_nn)
}



U.sgv <- function(dist.nn, neighbor_matrix, phi, nu, m) {
  
  n <- length(dist.nn)
  U.sparse <- spam(0, nrow = n, ncol = n)
  
    # Sequential version
    for (i in 1:n) {
      neighbors <- neighbor_matrix[i, ]
      neighbors <- neighbors[!is.na(neighbors)]
      
        
        K_n <- Matern(dist.nn[[i]], range = phi, smoothness = nu)
        prec_K_n <- chol2inv(chol(K_n))
        
        U_nn <- prec_K_n[1,]/sqrt(prec_K_n[1,1])
        
        U.sparse[i, neighbors] <- U_nn
  
    }
  
  return(U.sparse)
}


