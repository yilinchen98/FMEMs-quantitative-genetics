joint_FPCA_obj <- function(warping_data, C = 1, m = 3){
  #'
  #'@param warping_data fdawarp object from time_warping of aligned data
  #'@param C positive scaling parameter
  #'@param m number of principal components extracted
  
  f0 <- warping_data$f0 # original curve
  fn <- warping_data$fn # aligned curves
  time <- warping_data$time # time grid which curve data evaluated on
  gam <- warping_data$warping_functions # warping functions
  Tgam <- SqrtMean(gam) 
  mu_psi <- Tgam$mu #Karcher mean psi function
  vec <- Tgam$vec # shooting vectors
  
  n_col <- ncol(fn)
  n_row <- nrow(fn)
  ## joint curves and shooting vector
  g <- rbind(fn, C*vec)
  ## run FPCA on g
  fpcaobj <- prcomp(x=t(g), retx = TRUE, center = TRUE, rank. = m)
  eigenvecs <- fpcaobj$rotation # principal components
  scores <- fpcaobj$x # component scores
  
  mu_fn <- rowMeans(fn)
  fn_hat <- matrix(0, n_row, n_col) # estimated aligned curves
  for (i in 1:n_col ){
    fn_hat[,i] <- mu_fn + scores[i,] %*% t(eigenvecs[1:length(time),])
  }
  
  vec_hat <- matrix(0, n_row, n_col) # estimated shooting vectors
  for (i in 1:n_col){
    vec_hat[,i] <- (scores[i,]/C) %*% t(eigenvecs[length(time)+1:length(time),])
  }
  
  gam_hat <- v_to_gam(vec_hat) # estimated warping functions
  
  return(list(f0 = f0, fn_hat = fn_hat, vec_hat = vec_hat, gam_hat = gam_hat, eigenvecs = eigenvecs, time = time))
}

compute_C <- function(C, warping_data, m) {
  
  require(pracma)
  
  out.pca <- joint_FPCA_obj(warping_data, C = C, m = 3)
  f0 <- out.pca$f0
  fn_hat <- out.pca$fn_hat
  gam_hat <- out.pca$gam_hat
  n_col <- ncol(fn_hat)
  time <- out.pca$time
  d <- rep(0, n_col)
  for (i in 1:n_col) {
    tmp <- warp_f_gamma(fn_hat[,i], time, invertGamma(gam_hat[,i]))
    d[i] <- sum(trapz(time, (tmp - f0[, i])^2))
  }
  return(sum(d^2)/n_col)
}

C <- stats::optimize(compute_C, c(0, 100), warping_data = aligned_mass_process, m = m)$minimum
