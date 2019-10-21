
get.basis.matrix.piecewise.linear <- function(X, elements, seeds, R, nderiv=0) 
  {
  # Computes the basis matrix PHI(X); for piecewise linear basis {phi_i} evaluated at points in X, this 
  # function should not be user facing, i.e. all objects have been previously verified for type conformity
  # Arguments 
  # X ... k x d matrix of the k d-dimensional points where the basis is to be evaluated
  # ELEMENTS ... T x d+1 matrix of the tessellation of the space 
  # SEEDS ... n x d matrix of the coordinates of seeds (piecewise basis function 'centers')
  # R ... nxT dimensional (sparse) matrix defining the local ordering of node i in element j
  # NDERIV ... specifies the order of the derivative operator to be appled to basis functions prior to
  #             evaluation, e.g. n=0 -> phi_i(X), n=1 -> d(phi_i(X))/dx. For now, nderiv in (0, 1)
  # Returns
  # Phi ... k x n matrix of the evaluations of n-basis functions over the k points in X
  
  # Dimension of domain
  d <- ncol(seeds)
  
  # Associate the points in X with the tessellation elements which contain them
  coords <- tcovering(seeds, elements, X)
  
  # Evaluate points, note: consider using shape functions field in basis_object for computations 
  # will give same results but would make program usage more consistent
  if (nderiv==0) {
    # Allocate space for matrix PHI 
    Phi <- matrix(0, nrow=nrow(X), ncol=nrow(seeds))
    
    for (ix in 1:nrow(X)) {
      tix <- coords[["idx"]][ix] # tessellation element row index
      element_vec <- elements[tix, ] # ordered seed indices defining tix'th element
      eps <- coords[["p"]][ix, ] # barycentric coordinates of point X[ix, ]
      for (jx in element_vec) {
        Phi[ix, jx] <- eps[R[jx, tix]] # R[jx, tix] defines the linear lagrange shape function representation of
                                      # of the jix'th basis function over the tix'th tessellation element
      }
    } 
  } else if (nderiv==1) {
    # Allocate space for tensor Phi 
    Phi <- array(0, dim=c(nrow(X), nrow(seeds), d))
    
    for (ix in 1:nrow(X)) {
      tix <- coords[["idx"]][ix]
      element_vec <- elements[tix, ]
      for(jx in element_vec) {
        local_order <- R[jx, tix]
        if (local_order <= 2) {
          Phi[ix, jx, local_order] <- 1
        } else {
          Phi[ix, jx, ] <- rep(-1, d)
        }
      }
    }  
  } else {
      stop("Linear piecewise basis only first order differentiable!")
  }
    
  return(Phi)
}