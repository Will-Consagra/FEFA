
create.polynomial.basis <- function(seeds, boundary=NULL,
                         norder=1, quadord=1) 
  
  {
  # This is the constructor function for a linear-piecewise functional basis
  # Arguments
  # SEEDS ... a pxd dimensionsional matrix of knot points to center each piecewise basis element, 
  #             p = number of seeds, d = dimension of domain
  # BOUNDARY ... (not yet implemented) function defining the boundary in R^d, for now we consider convex hull
  #             of seeds to define boundary 
  # NORDER ... order of piecewise polynomial basis define on each element, for now fixed at norder=1
  # QUADORD ... order of integration for numerical quadrature rule, for now fixed at quadord=1 
  
  # set options warning level
  op <- options(warn=-1)
  options(op)
  
  # validate parameters
  
  if (!is.matrix(seeds)) {
    stop("seeds must be a matrix object! user specified seeds is of class ", class(seeds))
  } else {
    # domain dimension
    d <- ncol(seeds) 
  }
  
  if (norder == 1) {
    type <- "piecewise_linear"
  }
  else {
    stop('norder specified as ', norder, ' ... piecewise polynomial basis of > 1 have not yet been implemented')
  }
  
  if (quadord == 1) {
    if (d == 1) {
      quad_weight <- c(1/6, 2/3, 1/6)
      quad_points <- rbind(c(0), c(1/2), c(1))
    } else if (d == 2) {
      quad_weight <- rep(1/6, 3)
      quad_points <- rbind(c(1/2, 1/2), c(0, 1/2), c(1/2, 0))
    } else if (d == 3) {
      quad_weight <- c(1/240, 1/240, 1/240, 1/240, 3/80, 3/80, 3/80, 3/80)
      quad_points <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(0, 0, 0),
                           c(1/3, 1/3, 1/3), c(0, 1/3, 1/3), c(1/3, 0, 1/3), c(1/3, 1/3, 0))
    } else {
      stop("Numerical quadrature for d > 3 dimensional domain not yet implemented")
    }
    quadvals <- list("quad_weight"=quad_weight,
                     "quad_points"=quad_points)
  }
  
  else {
    stop('quadord specified as ', quadord, ' ... quadrature rules of order > 1 have not yet been implemented')
  }
  
  # construct tessellation of surface
  # for now this is restricted to be a d-simplicial complex, where simplices are associated to the space 
  # by the delaunay triangulation over the point cloud defined in seeds
  
  simplicial_complex = tryCatch(delaunayn(seeds), # note, probabaly should consider an object oriented approach to the tessellation
                              warning = function(w) {print(w)}, 
                              error = function(e) {stop(e)})
  
  # define the Jacobian mapping to reference element for purposes of numerical integration
  shape_functions <- get.shape.functions(d, ord=norder)
  
  # create r matirx: mapper between basis function i, tessellation element j and local function definition k
  # in the form R_ij <- k, where k \in {0, 1, ..., S} and 0 indicates the constant 0 function
  # note: enumeration technique for local shape functions is considered only for the linear functions over 
  # a reference simplex, higher order polynomials and alternative basis functions need to be considered seperately
  
  R <- matrix(0, nrow=nrow(seeds), ncol=nrow(simplicial_complex))
  for (jx in 1:nrow(simplicial_complex)) {
    element_jx <- simplicial_complex[jx, ]
    for (kx in 1:length(element_jx)) {
      ix <- element_jx[kx]
      R[ix, jx] <- kx
    }
  }
  R <- Matrix(R, sparse = TRUE)
  
  # create basis object 
  basis_object <- basis(seeds, simplicial_complex, shape_functions, R, type, quadvals)
  
  return(basis_object)
  
}
