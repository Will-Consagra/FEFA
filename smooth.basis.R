
# This file will house the various PDE based smoothing methods

smooth.heat <- function(repr_object, M=NULL, B=NULL, niter=10, dt=0.1, bw=NULL)
{
  # Smoothing via the heat equation with dirichlet boundary conditions 
  # Arguments 
  # REPR ... functional representation object
  # M ... kxk sparse matrix, pairwise inner product of basis functions 
  # B ... kxk sparse matrix, pairwsie inner product of derivative of basis functions 
  # niter ... number of iterations for Euler time solver 
  # dt ... time step 
  # bw ... bandwidth parameter, will override niter and use equivalence between heat equation and Guassian smoother to compute niter and dt 
  # Returns 
  # List containing ...
  # REPR ... functional representation object with post-smoothing coefficients 
  # notes 
  # need to implement stabality check of proposed forward Euler method 
  
  # verify input types 
  if (class(repr_object)!="repr") {
    stop("repr_object argument must be of class 'repr', not '", class(repr_object), ",")
  }
  
  if (!is.numeric(dt)) {
    stop("'dt' not numeric")
  }
  
  if (!is.numeric(niter)) {
    stop("'niter' not numeric")
  } else if (niter%%1!=0) {
    niter <- ceil(niter)
  }
  
  if (!is.null(bw)) {
    tf <- 0.5*((bw/(4*qnorm(3/4)))^2)
    niter <- ceil(tf/dt)
  }
  
  basis_object <- repr_object$basis 
  
  # construct inner product matrices if not supplied 
  
  if (is.null(M)) {
    M <- eval.inprod(basis_object, nderiv=0)
  } 
  
  if (is.null(B)) {
    B <- eval.inprod(basis_object, nderiv=1)
  }
  
  # get indices of boundary elements 
  dOmega <- unique(as.vector(basis_object$chull))
  
  # allocate matrix for solution trajectory 
  It <- Matrix(0, nrow=dim(M)[1], ncol=niter+1)
  
  # define initial solution and impose dirichlet boundary conditions
  It[, 1] <- dirichlet.boundary(repr_object$coefs, dOmega) 
  
  # integrate solution 
  for (i in 1:niter) {
    dIt <- solve(M, -B%*%It[, i])
    temp <- dIt*dt + It[, i]
    temp <- dirichlet.boundary(temp, dOmega)
    It[, i+1] <- temp
  }
  
  smoothed_trajectories <- list(repr_object=repr_object, smooth_params=list(dt=dt, niter=niter, bw=bw),
                                mass_mat=M, stiff_mat=B, coef_trajectories=It)
  
  oldClass(smoothed_trajectories) <- "heat_smoother"
  
  return(smoothed_trajectories)
  
}