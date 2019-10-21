

basis <- function(seeds, elements, shape_functions, R, chull, type, quadvals)
  
{
  # Generator (initialization) function for basis class, this function should not be user facing and 
  # accessible via create.basis constructor function
  # Arguments
  # SEEDS ... nxd dimensional matrix defining the coordinates of the basis functions
  # ELEMENTS ... Tx(d+1) dimensional matrix defining the tessellation of input domain
  # SHAPE_FUNCTIONS ... list of functions for computing local integrals 
  # R ... nxT dimensional (sparse) matrix defining the local ordering of node i in element j
  # CHULL ... R x d matrix, the R d-simplices defining the convexhull(seeds), given as indices into rows of seeds
  # TYPE ... Class of piecewise functions to be used as basis, for now fixed at "piecewise_linear"
  # QUADVALS ... list of two: <quad_weights, quad_points> defining the numerical quadrature rule
  
  if (type!="piecewise_linear") {
    stop("Basis function type ", type, " not yet implemented")
  }
  
  obj.call <- match.call()
  
  basis_object <- list(call=obj.call, seeds=seeds, elements=elements, 
                       shape_functions=shape_functions, R=R, chull=chull,
                       type=type, quadvals=quadvals)
  
  oldClass(basis_object) <- "basis"
  return(basis_object)
  
}

#  --------------------------------------------------------------------------
#                  Methods for basis class
#  --------------------------------------------------------------------------


print.basis <- function(x, ...) {
  warning("Print method not yet implemented")
}

summary.basisfd <- function(object, ...) {
  warning("Summary method not yet implemented")
}

"==.basis" <- function(basis1, basis2) {
  # Test of equality between two basis systems
  warning("Testing equaltiy not yet implemented")
}

"*.basis" <- function (basis1, basis2) {
  # Point-wise multiplication will expand basis1 by basis2, i.e. construct tessellation over seeds1 union seeds2
}
  
