read_matrix_from_folder <- function(folder, file) {
  full_filename <- paste(folder, file, sep="/")
  as.matrix(read.csv(full_filename, sep="\t", header=FALSE))
}

l2_norm <- function(x) {
  norm(x, type="2")
}

scale_factors_X_B_same_norm <- function(X, B) {
  # Scale so that columns of X and columns of B have same norm
  X_norms <- apply(X, MARGIN=2, FUN=l2_norm)
  B_norms <- apply(B, MARGIN=2, FUN=l2_norm)
  
  # Calculate the overall norm of each factor, then adjust the matrices so that for each
  # factor k X[, k] and B[, k] have the same L2 norm
  combined_norms <- X_norms * B_norms
  X_adjustment <- diag(sqrt(combined_norms) / X_norms)
  B_adjustment <- diag(sqrt(combined_norms) / B_norms)
  X_scaled <- X %*% X_adjustment
  B_scaled <- B %*% B_adjustment
  return(list("X_scaled"=X_scaled,
              "B_scaled"=B_scaled))
}

delete_below_threshold <- function(mat, threshold) {
  # Delete the elements that would be below the threshold if the matrix was scaled to have
  # columns with norm 1
  norms <- apply(mat, MARGIN=2, FUN=l2_norm)
  ifelse(abs(mat %*% diag(1 / norms)) > threshold, mat, 0)
}

read_thresholded_factor_matrices <- function(folder, threshold=0.01) {
  # Read in matrices
  X_raw <- read_matrix_from_folder(folder, "X.txt")
  B_raw <- read_matrix_from_folder(folder, "B.txt")
  
  # Scale so that columns of X and columns of B have same norm
  scaled <- scale_factors_X_B_same_norm(X_raw, B_raw)
    
  # Remove values smaller than threshold
  X_thresh <- delete_below_threshold(scaled$X_scaled, threshold)
  B_thresh <- delete_below_threshold(scaled$B_scaled, threshold)
  
  # Identify which factors have at least one gene and at least one sample
  X_nonempty <- (colSums(X_thresh != 0) != 0)
  B_nonempty <- (colSums(B_thresh != 0) != 0)
  both_nonempty <- X_nonempty & B_nonempty
  
  # Restrict to these non-empty factors
  X = X_thresh[, both_nonempty]
  B = B_thresh[, both_nonempty]
  
  return(list("X"=X,
              "B"=B,
              "K"=sum(both_nonempty)))
}

calc_factor_contribution <- function(biclustering, factor_index) {
  factor_contribution <- outer(biclustering$X[, factor_index],
        biclustering$B[, factor_index])
  return(factor_contribution)  
}

calc_nz_factor_contribution <- function(biclustering, factor_index) {
  nz_samples <- (biclustering$X[, factor_index] != 0)
  nz_genes <- (biclustering$B[, factor_index] != 0)
  factor_contribution <- outer(biclustering$X[nz_samples, factor_index],
                               biclustering$B[nz_genes, factor_index])
  return(factor_contribution)
}