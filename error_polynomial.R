library(combinat)
library(SimplicialCubature)

A <- function(p) {
  if (p == 1) {
    return(rbind(c(4)))
  } else if (p == 2) {
    return(rbind(c(4, 0), c(3,1), c(2,2), c(1,3), c(0,4)))
  } else if (p == 3) {
    return(rbind(
      c(4,0,0), 
      c(3,1,0), c(3,0,1),
      c(2,2,0), c(2,0,2), c(2,1,1),
      c(1,3,0),c(1,0,3),c(1,2,1),c(1,1,2),
      c(0,4,0), c(0,0,4), c(0,3,1), c(0,1,3), c(0,2,2)
    ))
  } else {
    alpha1 <- c(4, rep(0, p-1))
    alpha2 <- c(c(3, 1), rep(0, p-2))
    alpha3 <- c(c(2, 2), rep(0, p-2))
    alpha4 <- c(c(2,1,1), rep(0, p-3))
    alpha5 <- c(c(1,1,1,1), rep(0, p-4))
    return(cbind(
      unique(t(do.call(cbind, permn(alpha1)))),
      unique(t(do.call(cbind, permn(alpha2)))),
      unique(t(do.call(cbind, permn(alpha3)))),
      unique(t(do.call(cbind, permn(alpha4)))),
      unique(t(do.call(cbind, permn(alpha5))))
    ))
  }
}

G <- function(alpha) {
  indices <- unlist(c(lapply(1:length(alpha), function (x) {rep(x, alpha[x])})))
  return(unique(t(do.call(cbind, permn(indices)))))
}

coef_ij <- function(H_i, H_j, Ap) {
  coefs <- rep(0, nrow(Ap))
  for (i in 1:nrow(Ap)) {
    alpha_i <- Ap[i, ]
    G_alpha <- G(alpha_i)
    coef <- 0 
    for (j in 1:nrow(G_alpha)) {
      g_alpha <- G_alpha[j, ]
      coef <- coef + H_i[g_alpha[1],g_alpha[2]]*H_j[g_alpha[3],g_alpha[4]]
    }
    coefs[i] <- coef
  }
  return(coefs)
}

construct_error_poly <- function(Sigma, H) {
  #Sigma: Covariance matrix of coefficints defining random function space
  #H: list of Hessian matrices for each basis function (must be at least length 1)
  k <- nrow(Sigma)
  p <- nrow(H[[1]])
  A_p <- A(p)
  coefs <- rep(0, nrow(A_p))
  for (i in 1:k) {
    for (j in 1:k) {
      H_i <- H[[i]]
      H_j <- H[[j]]
      gamma_ij <- coef_ij(H_i, H_j, A_p) 
      coefs <- coefs + Sigma[i,j]*gamma_ij
    }
  }
  return(coefs)
}

integrate_error_poly <- function(Tau, x_tau, H_xtau, Sigma) {
  # Tau: simplex, p x p+1 
  # x_tau: centriod of simplex, 1xp 
  # H_xtau: length k list of pxp hessians evaluated at xtau 
  # Sigma: kxk covariance matrix of coefficients 
  p = ncol(x_tau)
  A_p = A(p)
  coefs_xval <- construct_error_poly(Sigma, H_xtau)
  error_poly <- definePoly(coefs_xval, A_p)
  Tau_centered <- sweep(Tau, 1, x_tau, "-")
  result <- integrateSimplexPolynomial(error_poly, Tau_centered, method="GM" )
  return(result$integral)
}


## ToDo: Test out this code when simulating functions from random tensor basis!!!!

p <- 3 
Ap <- A(p)
alpha_1 <- Ap[1, ]
alpha_2 <- Ap[7, ]
g_alpha_1 <- G(alpha_1)
g_alpha_2 <- G(alpha_2)
