
representation.basis <- function(basis_object, X, Y) 
{
  # Represent a set of discrete observations to a smooth, pre-defined basis system
  # Arguments 
  # BASIS_OBJECT ... object of class "basis", defining the smooth basis system to be fit to observations 
  # X ... nxd dimensional matrix of the n-observations points in d-dimensional (covariate) space
  # Y ... nx1 vector of function observations at covariates X 
  # Returns 
  # List containing...
  # REPR ... functional representation object 
  # DF ... estimated degrees of freedom, computed via trace of projection operator 
  # BETA_HAT ... kx1 vector of regression coefficients for each basis element 
  # SSE ... error sum of squares 
  # GCV ... estimate of generalized cross validation error, not yet derived for FE-type basis systems 
  # PROJECTION_MAT ... nxn 'projection' matrix 
  # 
  # note: need additional functionality for dealing with numerical instabilities that can appear in matrix inversion steps 
  
  # check arguments 
  
  if (!is.numeric(Y)){
    stop("Observation vector Y must be numeric!")
  }
  
  # coerce Y -> type 'matrix'
  if (is.vector(Y)){
    Y <- as.matrix(Y)
  }
  
  if (class(basis_object)!="basis") {
    stop("basis_object argument must be of class 'basis', not '", class(basis_object), ",")
  }
  
  if (!is.matrix(X)) {
    stop("Covariate matrix X must be n x d dimensional 'matrix' class, not '", class(X), "'")
  }
  
  # validate point/seed dimensions are consistent
  if (ncol(X)!=ncol(basis_object$seeds)){
    stop("Dimension of points: ", ncol(X), " != Dimension of seeds: ", ncol(basis_object$seeds))
  }
  
  # validate number of observations and covariate vectors are the same 
  if (nrow(X)!=nrow(Y)){
    stop("Number of rows in covariate matrix : ", nrow(X), " != Number of observations: ", nrow(Y))
  }
  
  nobs <- nrow(X); nbasis <- nrow(basis_object$seeds)
  
  if (nbasis > nobs){
    stop("Representation via a basis function system of greater number than observations not defined!")
  }
  
  # compute matrix of basis function values
  
  Phi <- eval.basis(X, basis_object, nderiv=0) 
  
  # set up normal equations and solve for betahat, !note: singularity issues here -> add some regularization
  
  # betahat <- solve(t(Phi)%*%Phi)%*%t(Phi)%*%Y
  
  # using qr factorization to solve normal equations
  
  QR <- qr(Phi)
  betahat <- qr.coef(QR, Y)

  
  # compute projection matrix: unstable to compute solve(t(Phi)%*%Phi) directly
  
  S <- Phi%*%solve(t(Phi)%*%Phi)%*%t(Phi)
  
  # compute degrees of freedom
  
  df <- sum(diag(S))
  
  # compute error sum of squares
  
  Yhat <- Phi%*%betahat 
  SSE  <- sum((Y - Yhat)^2)
  
  # compute GCV-like index; not designed yet 

  gcv <- (nobs/(nobs - df))*(SSE/(nobs - df))
  
  # construct repr (representation) object
  
  repr_basis <- repr(betahat, basis_object)

  representation <- list(repr_basis=repr_basis, df=df, beta_hat=betahat,
                  SSE=SSE, GCV=gcv, projection_mat=S)
  
  
  return(representation)
  
}