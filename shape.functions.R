
get.shape.functions <- function(dim, ord=1) 
  {
  # Generates a list of functions defining the functional basis and associated Jacobian  
  # and coordinate mappings on reference simplex
  # DIM: int, dimension of reference simplex
  # ORD: int, order of polynomial basis
  # note: should be able to refactor to build required functions dynamically as opposed to 
  # hardcoding for each dimension
  if (dim == 1) {
    # Linear shape functions on reference 1-Simplex
    # psi is barycentric coordinate: psi = e
    V <- list(
      "v1"= function(psi) {
        return(psi[1])
      },
      "v2"= function(psi) {
        return(1-psi[1])
      }
    )
    # Gradient of shape functions: dv_s/dpsi (psi)
    dV <- list(
      "dv1"= function(psi) {
        return(1)
      },
      "dv2"= function(psi) {
        return(-1)
      }
    )
    # Jacobian Function
    # X_l: <- 1x2 matrix of vertice vectors in proper ordering 
    J <- function(X_l) {
      x_0 <- X_l[1]; x_1 <- X_l[2]
      return(x_1 - x_0)
    }
    # Jacobian inverse
    Jinv <- function(X_l) {
      x_0 <- X_l[1]; x_1 <- X_l[2]
      return(1/(x_1 - x_0))
    }
    # Jacobian determinant
    det_J <- function(X_l) {
      x_0 <- X_l[1]; x_1 <- X_l[2]
      return(x_1 - x_0)
    }
    # Coordinate Map
    # X_l: Global coordinates of veritces of tesselation element 
    # psi: Local coordinates psi in (0, 1) 
    mu <- function(psi, X_l) {
      x_0 <- X_l[1]; x_1 <- X_l[2];
      return(psi[1]*x_0 + (1 - psi[1])*x_1)
    } 
  } else if (dim == 2) {
    # Linear shape functions on reference 2-Simplex
    # psi is barycentric coordinate: psi = (e_0, e_1)
    V <- list(
      "v1"= function(psi) {
        return(psi[1])
      },
      "v2"= function(psi) {
        return(psi[2])
        },
      "v3"= function(psi) {
        return(1 - psi[1] - psi[2])
      }
    )
    # Gradient of shape functions: dv_s/dpsi (psi)
    dV <- list(
      "dv1"= function(psi) {
        return(c(1, 0))
      },
      "dv2" = function(psi) {
        return(c(0, 1))
      },
      "dv3"= function(psi) {
        return(c(-1, -1))
      }
    )
    # Jacobian Function
    # X_l: <- 3x2 matrix of vertice vectors in proper ordering 
    J <- function(X_l) {
      x_0 <- X_l[1, ]; x_1 <- X_l[2, ]; x_2 <- X_l[3, ]
      return(rbind( c(x_0[1] - x_2[1], x_1[1] - x_2[1]),
                    c(x_0[2] - x_2[2], x_1[2] - x_2[2])
      ))
    }
    # Jacobian inverse
    Jinv <- function(X_l) {
      x_0 <- X_l[1, ]; x_1 <- X_l[2, ]; x_2 <- X_l[3, ]
      detJ <- ((x_0[1] - x_2[1])*(x_1[2] - x_2[2])) - ((x_1[1] - x_2[1])*(x_0[2] - x_2[2]))
      return((1/detJ)*rbind( c(x_1[2] - x_2[2], x_2[1] - x_1[1]),
                             c(x_2[2] - x_0[2], x_0[1] - x_2[1]) 
      ))
    }
    # Jacobian determinant
    det_J <- function(X_l) {
      x_0 <- X_l[1, ]; x_1 <- X_l[2, ]; x_2 <- X_l[3, ]
      return( ((x_0[1] - x_2[1])*(x_1[2] - x_2[2])) - ((x_1[1] - x_2[1])*(x_0[2] - x_2[2])) )
    }
    # Coordinate Map
    # X_l: Global coordinates of vertices of tesselation element 
    # psi: Local coordinates (barycentric coordinates): (e_0, e_1) 
    mu <- function(psi, X_l) {
      x_l0 <- X_l[1, ]; x_l1 <- X_l[2, ]; x_l2 <- X_l[3, ]
      return(psi[1]*x_l0 + psi[2]*x_l1 + (1-psi[2]-psi[1])*x_l2)
    } 
  } else if (dim == 3) {
    # Linear shape functions on reference 3-Simplex
    # psi is barycentric coordinate: psi = (e_0, e_1, e_2)
    V <- list(
      "v1"= function(psi) {
        return(psi[1])
      },
      "v2"= function(psi) {
        return(psi[2])
      },
      "v3"= function(psi) {
        return(psi[3])
      },
      "v4"= function(psi) {
        return(1 - psi[1] - psi[2] - psi[3])
      }
    )
    # Gradient of shape functions: dv_s/dpsi (psi)
    dV <- list(
      "dv1"= function(psi) {
        return(c(1, 0, 0))
      },
      "dv2"= function(psi) {
        return(c(0, 1, 0))
      },
      "dv3"= function(psi) {
        return(c(0, 0, 1))
      },
      "dv4"= function(psi) {
        return(c(-1, -1, -1))
      }
    )
    # Jacobian Function
    # X_l: <- 4x3 matrix of vertice vectors in proper ordering 
    J <- function(X_l) {
      x_0 <- X_l[1, ]; x_1 <- X_l[2, ]; x_2 <- X_l[3, ]; x_3 <- X_l[4, ]
      return(rbind(c(x_0[1]-x_3[1], x_0[2]-x_3[2], x_0[3]-x_3[3]),
                   c(x_1[1]-x_3[1], x_1[2]-x_3[2], x_1[3]-x_3[3]),
                   c(x_2[1]-x_3[1], x_2[2]-x_3[2], x_2[3]-x_3[3])))
    }
    # Jacobian inverse
    Jinv <- function(X_l) {
      x_0 <- X_l[1, ]; x_1 <- X_l[2, ]; x_2 <- X_l[3, ]; x_3 <- X_l[4, ]
      J <- rbind(c(x_0[1]-x_3[1], x_0[2]-x_3[2], x_0[3]-x_3[3]),
                 c(x_1[1]-x_3[1], x_1[2]-x_3[2], x_1[3]-x_3[3]),
                 c(x_2[1]-x_3[1], x_2[2]-x_3[2], x_2[3]-x_3[3]))
      return(solve(J))
    }
    # Jacobian determinant
    det_J <- function(X_l) {
      x_0 <- X_l[1, ]; x_1 <- X_l[2, ]; x_2 <- X_l[3, ]; x_3 <- X_l[4, ]
      J <- rbind(c(x_0[1]-x_3[1], x_0[2]-x_3[2], x_0[3]-x_3[3]),
                 c(x_1[1]-x_3[1], x_1[2]-x_3[2], x_1[3]-x_3[3]),
                 c(x_2[1]-x_3[1], x_2[2]-x_3[2], x_2[3]-x_3[3]))
      return( J[1,1]*(J[2,2]*J[3,3] - J[2,3]*J[3,2]) - 
                J[1,2]*(J[2,1]*J[3,3] - J[2,3]*J[3,1]) + 
                J[1,3]*(J[2,1]*J[3,2] - J[2,2]*J[3,1]) )
    }
    # Coordinate Map
    # X_l: Global coordinates of veritces of tesselation element 
    # psi: Local coordinates (barycentric coordinates): (e_0, e_1, e_2) 
    mu <- function(psi, X_l) {
      x_l0 <- X_l[1, ]; x_l1 <- X_l[2, ]; x_l2 <- X_l[3, ]; x_l3 <- X_l[4, ]
      return(psi[1]*x_l0 + psi[2]*x_l1 + psi[3]*x_l2 + (1-psi[3]-psi[2]-psi[1])*x_l3)
    }
  } else {
    stop("Selected dimension not yet implemented! Please choose dim in c(1, 2, 3)")
  }
  return(list("V"=V,
              "dV"=dV,
              "J"=J,
              "Jinv"=Jinv,
              "det_J"=det_J,
              "mu"=mu))
} 
