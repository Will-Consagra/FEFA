
representation.basis <- function(basis_object, X, Y, pen_mat=NULL, lambda=0) 
{
  # Represent a set of discrete observations to a smooth, pre-defined basis system
  # Arguments 
  # BASIS_OBJECT ... object of class "basis", defining the smooth basis system to be fit to observations 
  # X ... nxd dimensional matrix of the n-observations points in d-dimensional (covariate) space
  # Y ... nx1 vector of function observations at covariates X 
  # PEN_MAT ... k x k dimensional penalty matrix for regularization of least-squares coefficient learning
  # LAMBDA ... numeric, strength of regularization
  # Returns 
  # List containing...
  # REPR ... functional representation object 
  # DF ... estimated degrees of freedom, computed via trace of projection operator 
  # BETA_HAT ... kx1 vector of regression coefficients for each basis element 
  # SSE ... error sum of squares 
  # GCV ... estimate of generalized cross validation error, not yet derived for FE-type basis systems 
  # PROJECTION_MAT ... nxn 'projection' matrix 
  # R_SQ ... coefficient of determination = 1 - SSE/SS_tot 
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
  
  # validate pen_mat type or create sparse 
  if (!is.null(pen_mat)) {
    if (!is.matrix(pen_mat)) {
      stop("'pen_mat' must be of matrix type!")
    } else if (nrow(pen_mat)!=ncol(pen_mat)) {
      stop("'pen_mat' must be square matrix!")
    } else if (nrow(pen_mat)!=nbasis) {
      stop("'pen_mat' must be of dimension num_basis x num_basis!")
    }
  } else {
    # eventually want to allow for sparse matrices when nbasis gets large
    #pen_mat <- sparseMatrix(i = 1:nbasis, j = 1:nbasis, x = 1) 
    #pen_mat <- eye(nbasis) 
    pen_mat <- diag(x=1, nrow=nbasis, ncol=nbasis)
  }
  
  # validate type of 'lambda' 
  if (!is.numeric(lambda)) {
    stop("'lambda' must be of type numeric")
  }
  
  # compute matrix of basis function values
  
  Phi <- eval.basis(X, basis_object, nderiv=0) 
  
  # set up normal equations 
  
  if (lambda > 0) {
    pen_mat_sqrt <- sqrtm(pen_mat)$B # consider more stable solution to compute sqrt(pen_mat) using eigen-decomp
    basis_mat <- rbind(Phi, sqrt(lambda)*pen_mat_sqrt)
    Y_tilde <- rbind(Y, matrix(0, nrow=nrow(pen_mat_sqrt), ncol=1))
  } else {
    basis_mat <- Phi
    Y_tilde <- Y
  }
  
  # use qr factorization to solve normal equations
  
  qr <- qr(basis_mat)
  betahat <- qr.coef(qr, Y_tilde)
  
  # compute projection matrix: unstable to compute solve(t(Phi)%*%Phi) directly
  
  S <- Phi%*%solve(t(basis_mat)%*%basis_mat)%*%t(Phi)
  
  # compute degrees of freedom
  
  df <- sum(diag(S))
  
  # compute error sum of squares
  
  Yhat <- Phi%*%betahat 
  SSE  <- sum((Y - Yhat)^2)
  
  # compute GCV 

  gcv <- (nobs/(nobs - df))*(SSE/(nobs - df))
  
  # compute coefficient of determination 
  
  r_sq <- 1 - (SSE/sum((Y-mean(Y))^2))
  
  # construct repr (representation) object
  
  repr_basis <- repr(betahat, basis_object)

  representation <- list(repr_basis=repr_basis, df=df, beta_hat=betahat,
                  SSE=SSE, GCV=gcv, projection_mat=S, r_sq=r_sq)
  
  
  return(representation)
  
}