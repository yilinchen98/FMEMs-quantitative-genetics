compute_C <- function(warping_data, m = 3){
  #'
  #'@param warping_data fdawarp object from time_warping of aligned data
  #'@param m number of principal components 
  #'
  #'@return C positive scaling parameter which joins the curve data and warping functions
  
  require(fdasrvf)
  
  f0 <- warping_data$f0 # original curve
  fn <- warping_data$fn # aligned curves
  time <- warping_data$time # time grid which curve data evaluated on
  gam <- warping_data$warping_functions # warping functions
  Tgam <- SqrtMean(gam) 
  vec <- Tgam$vec # shooting vectors
  
  jointFPCAobj <- function(fn, vec, C=1, m=m){
    n_col <- ncol(fn)
    n_row <- nrow(fn)
    time_new <- seq(0,2, length.out = n_row + nrow(vec))
    ## joint curves and shooting vector
    g <- rbind(fn, C*vec)
    ## run FPCA on g
    fpcaobj <- prcomp(x=t(g), retx = TRUE, center = TRUE, rank. = m)
    eigenvecs <- fpcaobj$rotation # principal components
    scores <- fpcaobj$x # component scores
    
    mu_fn <- rowMeans(fn) # mean of aligned curves
    fn_hat <- matrix(0, n_row, n_col) # estimated aligned curves
    for (i in 1:n_col ){
      fn_hat[,i] <- mu_fn + scores[i,] %*% t(eigenvecs[1:length(time),])
    }
    
    vec_hat <- matrix(0, n_row, n_col) # estimated shooting vectors
    for (i in 1:n_col){
      vec_hat[,i] <- (scores[i,]/C) %*% t(eigenvecs[length(time)+1:length(time),])
    }
    
    gam_hat <- v_to_gam(vec_hat) # estimated warping functions
    
    return(list(f0 = f0, fn = fn, fn_hat = fn_hat, vec_hat = vec_hat, gam_hat = gam_hat, eigenvecs = eigenvecs, time = time))
  }
  
  compute_C <- function(C, warping_data, m) {
    
    require(pracma)
    
    out.pca <- jointFPCAobj(fn, vec, C = C, m = m)
    f0 <- out.pca$f0
    fn_hat <- out.pca$fn_hat
    gam_hat <- out.pca$gam_hat
    n_col <- ncol(fn_hat)
    time <- out.pca$time
    d <- rep(0, n_col)
    for (i in 1:n_col) {
      tmp <- warp_f_gamma(fn_hat[,i], time, gam_hat[,i])
      d[i] <- sum(trapz(time, (tmp - fn[, i])^2))
    }
    return(sum(d^2)/n_col)
  }
  
  C <- stats::optimize(compute_C, c(0, 1000), warping_data = warping_data, m = m)$minimum
  
  return(C)
}