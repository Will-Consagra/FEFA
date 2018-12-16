
eval.basis <- function(X, basis_object, nderiv=0) 
  {
  # Computes the basis matrix PHI(X); basis_object evaluated at points in X.
  # Arguments 
  # X ... k x d matrix of the k d-dimensional points where the basis is to be evaluated
  # BASIS_OBJECT ... object of class "basis"
  # NDERIV ... specifies the order of the derivative operator to be appled to basis functions prior to
  #             evaluation, e.g. n=0 -> phi_i(X), n=1 -> d(phi_i(X))/dx. For now, nderiv in (0, 1)
  # Returns
  # Phi ... k x n matrix of the evaluations of n-basis functions over the k points in X
  
  # check arguments 
  if (!is.matrix(X)) {
    stop("Evaluation points X must be k x d dimensional 'matrix' class, not '", class(X), "'")
  }
  
  if (class(basis_object)!="basis") {
    stop("basis_object argument must be of class 'basis', not '", class(basis_object), ",")
  }
  
  # validate point/seed dimensions are consistent
  if (ncol(X)!=ncol(basis_object$seeds)){
    stop("Dimension of points: ", ncol(X), " != Dimension of seeds: ", ncol(basis_object$seeds))
  }
  
  # check 'type' of piecewise basis and call appropriate matrix constructor function 
  type <- basis_object$type 
  if (type=="piecewise_linear") {
    Phi <- get.basis.matrix.piecewise.linear(X, basis_object$elements, 
                                             basis_object$seeds, basis_object$R,
                                             nderiv) 
  } else {
    stop("eval.basis not yet implemented for basis of type '", type, "'")
  }
  
  return(Phi)
  
}