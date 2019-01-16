
get.rectangular.boundary <- function(Grid) 
  {
  # Computes the boundary points of a rectangular lattice 
  # GRID ... M x d matrix of lattice points 
  # RETURNS
  # DOMEGA ... list of indices corresponding to boundary points of GRID
  dOmega <- c()
  d <- ncol(Grid)
  if (d == 1) {
    dOmega <- c(1, nrow(Grid))
  } else if (d == 2) { ## Note: d == 2 works only for regular rectangular grids
    min_d1 <- min(Grid[, 1]); max_d1 <- max(Grid[, 1])
    min_d2 <- min(Grid[, 2]); max_d2 <- max(Grid[, 2])
    dOmega <- c(which(Grid[, 1] == min(Grid[, 1])), which(Grid[, 1] == max(Grid[, 1])),
                which(Grid[, 2] == min(Grid[, 2])), which(Grid[, 2] == max(Grid[, 2])))
  } else if (d == 3) {
    hull <- convhulln(Grid)
    for (ix in 1:nrow(hull)) {
      face <- hull[ix, ]
      P <- t(sapply(face, function (node) {Grid[node, ]}))
      common_dim <- which(apply(P, 2, function (x) {length(unique(x))}) == 1)
      plane_ix <- which(Grid[, common_dim] == P[1, common_dim])
      plane <- Grid[plane_ix, -common_dim]
      BP <- cart2bary(P[, -common_dim], plane)
      boundary_ix <- plane_ix[apply(BP >= 0, 1, function(row) {all(row)})]
      dOmega <- c(dOmega, boundary_ix)
    }
  }
  return(unique(dOmega))
}

dirichlet.boundary <- function(U, dOmega, val=0) 
  {
  # Enforces dirichlet boundary conditions
  U[dOmega] <- val
  return(U)
}

get.convex.hull <- function(X) 
  {
  # Computes the convex hull of X in the form of a matrix of indices of X
  # note: This function is a thin wrapper for convhulln from 'geometry' package, which is itself a thin wrapper 
  # for the underlying Qhull command qconvex.
  # X ... M x d matrix, point cloud of which we wish to determine the convex hull 
  # RETURNS
  # CHULL ... R x d matrix, the R d-simplices defining the convexhull(X), given as indices into rows of X 
  chull <- convhulln(X, options="Tv")
  return(chull)
}