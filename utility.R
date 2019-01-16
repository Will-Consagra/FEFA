
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

##Map index arrays to d+1 x d matrices of points corresponding to vertices in array
##Must handle the d = 1, length(index) == 2, separately due to quirkiness of r 
index.to.element <- function(element_ix, Grid) {
  # Map index arrays to d+1 x d matrices of points corresponding to vertices in array
  # Must handle the d = 1, length(index) == 2, separately due to quirkiness of r
  # ELEMENT_IX ... numeric vector defining element via row indices in grid
  # Grid ... nxd matrix of coordinates
  X_l <- t(sapply(element_ix, function(e) Grid[e, ]))
  if (length(element_ix) == 2) {
    X_l <- t(X_l)
  } 
  return(X_l)
}

