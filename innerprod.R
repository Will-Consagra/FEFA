
innerprod <- function(basis_object, nderiv=0) 
{
  # Computes the L2 pairwise inner product matrix of D^nderiv basis_object1, this function is not user facing
  # BASIS_OBJECT ... object of class "basis"
  # NDERIV ... specifies the order of the derivative operator to be appled to basis functions prior to
  #             evaluation, e.g. n=0 -> phi_i(X), n=1 -> d(phi_i(X))/dx. For now, nderiv in (0, 1)
  # Returns
  # M^{nderiv} ... n x n matrix of the pairwise inner products of between the nderiv^{th} of the basis functions 
  
  if (nderiv==0) {
    
    V <- basis_object$shape_functions$V; 
    det_J <- basis_object$shape_functions$det_J
    quad_points <- basis_object$quadvals$quad_points
    quad_weight <- basis_object$quadvals$quad_weight
    Grid <- basis_object$seeds
    elements <- basis_object$elements

    S <- length(V)
    Q <- nrow(quad_points)
    num_nodes <- nrow(Grid)
    d <- ncol(Grid)
    
    M <- matrix(rep(0, num_nodes*num_nodes), nrow=num_nodes, ncol=num_nodes)
    
    for(ix in 1:nrow(elements)) {
      # Get tessellation element
      element_ix <- elements[ix, ]
      X_l <- index.to.element(element_ix, Grid)
      
      # Jacobian information to compute integral in local coordinates
      Jac_det <- det_J(X_l)
      j_l <- abs(Jac_det)*quad_weight
      
      for (i in 1:length(element_ix)) {
        for (j in 1:length(element_ix)) {
          for (q in 1:Q) {
            M[element_ix[i],element_ix[j]] <- M[element_ix[i],element_ix[j]] +
              (V[[i]](quad_points[q, ])*V[[j]](quad_points[q ,])*j_l[q])
          }
        }
      }
    }
    return(Matrix(M, sparse=TRUE))
    
  } else if (nderiv==1) {
     
    dV <- basis_object$shape_functions$dV
    det_J <- basis_object$shape_functions$det_J
    quad_points <- basis_object$quadvals$quad_points
    quad_weight <- basis_object$quadvals$quad_weight
    Grid <- basis_object$seeds
    elements <- basis_object$elements
    
    S <- length(dV)
    Q <- nrow(quad_points)
    num_nodes <- nrow(Grid)
    d <- ncol(Grid)
    
    B <- matrix(rep(0, num_nodes*num_nodes), nrow=num_nodes, ncol=num_nodes)
    
    for (ix in 1:nrow(elements)) {
      
      element_ix <- elements[ix, ]
      X_l <- index.to.element(element_ix, Grid)
      Jac_det <- det_J(X_l)
      j_l <- abs(Jac_det)*quad_weight
      
      for (i in 1:length(element_ix)) {
        for (j in 1:length(element_ix)) {
          for (q in 1:Q) {
            B[element_ix[i], element_ix[j]] <- B[element_ix[i], element_ix[j]] + 
              (dV[[i]](quad_points[q, ]) %*% dV[[j]](quad_points[q, ]))*j_l[q]
          }
        }
      }
    }
    return(Matrix(B, sparse=TRUE))
    
  } else {
    stop("Implementation for computation of inner product matrices for D^{i}phi(X) for i > 1 not yet implemented!")
  }
}