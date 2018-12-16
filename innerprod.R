
innerprod <- function(basis_object, nderiv=0) 
{
  # Computes the L2 pairwise inner product matrix of D^nderiv basis_object1, this function is not user facing
  # BASIS_OBJECT ... object of class "basis"
  # NDERIV ... specifies the order of the derivative operator to be appled to basis functions prior to
  #             evaluation, e.g. n=0 -> phi_i(X), n=1 -> d(phi_i(X))/dx. For now, nderiv in (0, 1)
  # Returns
  # M^{nderiv} ... n x n matrix of the pairwise inner products of between the nderiv^{th} of the basis functions 
  
  if (nderiv==0) {
    
  } else if (nderiv==1) {
    
  } else {
    stop("Implementation for computation of inner product matrices for D^{i}phi(X) for i > 1 not yet implemented!")
  }
}