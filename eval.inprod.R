
eval.inprod <- function(basis_object, nderiv=0) 
{
  # Wrapper function for computing pairwise inner product matrix of nderiv^{th} basis_object
  # BASIS_OBJECT ... object of class "basis"
  # NDERIV ... specifies the order of the derivative operator to be appled to basis functions prior to
  #             evaluation, e.g. n=0 -> phi_i(X), n=1 -> d(phi_i(X))/dx. For now, nderiv in (0, 1)
  # Returns
  # M^{nderiv} ... n x n matrix of the pairwise inner products of between the nderiv^{th} of the basis functions 
  # Note, eventually we may like to consider the pairwise inner product computation between 2 distinct basis_objects
  # over the space, basis_object1 and basis_object2
  
  # check arguments 
  
  if (class(basis_object)!="basis") {
    stop("basis_object argument must be of class 'basis', not '", class(basis_object), ",")
  }
  
  # check 'type' of piecewise basis and call appropriate matrix constructor function 
  type <- basis_object$type 
  if (type=="piecewise_linear") {
    if (nderiv > 1) {
      stop("Piecewise linear basis only up to 1st order differentiable!")
    } 
    
    #M <- innerprod(basis_obj)
  } else {
    stop("eval.inprod not yet implemented for basis of type '", type, "'")
  }
  
  return(Phi)
  
}