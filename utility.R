
tcovering <- function(seeds, elements, X) 
  {
  # Associates each d-dimensional point in X with its covering simplex in elements. 
  # This function is not user facing. 
  # Arguments
  # SEEDS ... n x d matrix of the coordinates of seeds (piecewise basis function 'centers')
  # ELEMENTS ... T x d+1 matrix of the tessellation of the space 
  # X ... k x d matrix of the k d-dimensional points
  # Returns 
  # coords ... 2 element list <idx, eps>, where idx is indext in ELEMENTS where point is found and eps is point in 
  #             barycentric coordinates with respect to the idx row in ELEMENTS
  
  d <- ncol(X)
  if (d!=ncol(seeds)) {
    stop("Dimensionality of point cloud X and seed points seeds are incompatible")
  }
  
  if (d==1) {
    # these are just intervals and the search is trivial
    stop("Matrix construction for d=1 not yet implemented!!!!")
  }
  else {
    # use the function from geometry, this may not be optimal
    coords <- tsearchn(seeds, elements, X, fast=TRUE)
  }
  return(coords)
}

