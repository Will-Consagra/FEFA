
repr <- function(coef, basis_object)
{
  # This function constructs a functional-representation object from a basis function object and associated set of coefficients. 
  # Arguments 
  # COEF ... kx1 vector of the coefficients corresponding to each of the k basis functions. 
  # BASIS_OBJECT ... object of class basis 
  # Returns 
  # REPR ... representation (repr) object 
  
  # validate types 
  
  if(!is.numeric(coef)) {
    stop("coef is not numeric")
  } else if (is.vector(coef)) {
    coef <- as.matrix(coef)
  }
  
  if (class(basis_object)!="basis") {
    stop("basis_object argument must be of class 'basis', not '", class(basis_object), ",")
  }
  
  repobj <- list(coefs=coef, basis=basis_object)
  oldClass(repobj) <- "repr"
  repobj
  
}

#  --------------------------------------------------------------------------
#                  Methods for repr class
#  --------------------------------------------------------------------------


print.repr <- function(x, ...) {
  warning("Print method not yet implemented")
}

summary.repr <- function(object, ...) {
  warning("Summary method not yet implemented")
}

"+.repr" <- function(e1, e2){
  warning("Plus operation not yet implemented")
}

"-.repr" <- function(e1, e2){
  warning("Minus operation not yet implemented")
}

"*.repr" <- function (e1, e2) {
  warning("Point-wise multiplication of rep objects not yet implemented")
}


  