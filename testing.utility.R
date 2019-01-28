
test.1 <- function(X, Y, seeds, nderiv=0, sliced=FALSE) 
{
  # Computes projection of data onto linear-piecewise basis spanned by elements over seeds and generates visualizations and fit statistics 
  # X ... NxP matrix of N-points in P dimensional space
  # Y ... Nx1 matrix of funtion observations at each of the P-dimensional points in X 
  # SEEDS ... MxK matrix of M-points in P dimensional space, seeds or knot-points of linear piecewise basis 
  # NDERIV ... int \in (0, 1), project onto basis or derivative of basis 
  # SLICED ... (applicable for 3D) integer if seeds are equispaced, defines the extra dimensional slice for plotting data 
  
  n <- nrow(X); p <- ncol(X)
  basis_obj <- create.polynomial.basis(seeds)
  if (p==2) {
    trimesh(basis_obj$elements, basis_obj$seeds, main="Tesselation")
    title(main="Tesselation")
  }
  rep_result <- representation.basis(basis_obj, X, Y) 
  
  MSE <- rep_result$SSE/n
  # display fit statistics 
  
  print(paste("Fit statistics...", paste0("SSS: ", rep_result$SSE), 
         paste0("MSE: ", MSE), paste0("df: ", rep_result$df),
         paste0("GCV: ", rep_result$GCV), collapse="\n"))
  
  # visualization of pairwise dimensions 

  if (p==2) {
    interpolation <- interp(basis_obj$seeds[,1], basis_obj$seeds[,2], rep_result$beta_hat)
    filled.contour(interpolation$x, interpolation$y, interpolation$z,
                   color.palette=colorRampPalette(c("lightyellow", "orange", "red"), bias=2),
                   xlab=paste0("dimension ", 1), ylab=paste0("dimension ", 2),
                   main="FE Approximation",
                   key.title = title(main="Value", cex.main=0.8))
  } else if (p==3 && sliced!=FALSE) {
    
    xslice <- basis_obj$seeds[,1]==sliced
    yslice <- basis_obj$seeds[,2]==sliced
    zslice <- basis_obj$seeds[,3]==sliced
  
    interpolation_xy <- interp(basis_obj$seeds[zslice,1], basis_obj$seeds[zslice,2], rep_result$beta_hat[zslice])
    interpolation_xz <- interp(basis_obj$seeds[yslice,1], basis_obj$seeds[yslice,3], rep_result$beta_hat[yslice])
    interpolation_yz <- interp(basis_obj$seeds[xslice,2], basis_obj$seeds[xslice,3], rep_result$beta_hat[xslice])
    
    filled.contour(interpolation_xy$x, interpolation_xy$y, interpolation_xy$z,
                   color.palette=colorRampPalette(c("lightyellow", "orange", "red"), bias=2),
                   xlab="d1", ylab="d2",
                   main="FE Approximation",
                   key.title = title(main="Value", cex.main=0.8))
    
    filled.contour(interpolation_xz$x, interpolation_xz$y, interpolation_xz$z,
                   color.palette=colorRampPalette(c("lightyellow", "orange", "red"), bias=2),
                   xlab="d1", ylab="d3",
                   main="FE Approximation",
                   key.title = title(main="Value", cex.main=0.8))
    
    filled.contour(interpolation_yz$x, interpolation_yz$y, interpolation_yz$z,
                   color.palette=colorRampPalette(c("lightyellow", "orange", "red"), bias=2),
                   xlab="d2", ylab="d3",
                   main="FE Approximation",
                   key.title = title(main="Value", cex.main=0.8))
  }
    
}

test.2 <- function(seeds) 
{
  # Computes the inner-product matrices of basis elements and first derivatives and gives some visualtions and statistics 
  # SEEDS ... MxK matrix of M-points in P dimensional space, seeds or knot-points of linear piecewise basis 
  
  basis_obj <- create.polynomial.basis(seeds)
  
  tm <- system.time(M <- eval.inprod(basis_obj, nderiv=0))
  tb <- system.time(B <- eval.inprod(basis_obj, nderiv=1))
  print(paste("Inner-product matrix construction analysis... ", 
              paste0("Number of basis functions: ", nrow(M)), 
              paste0("Dimension of parameter space: ", ncol(seeds)), 
              paste0("Time to construct basis inner product matrix ", tm[["elapsed"]]),
              paste0("Basis inner product matrix symmetric? ", isSymmetric(M)),
              paste0("Time to construct derivative-of-basis inner product matrix ", tb[["elapsed"]]),
              paste0("Derivative-of-basis inner product matrix symmetric? ", isSymmetric(B)),
              collapse="\n"))
  
  par(mfrow=c(2,1))
  
  image(M, main="Basis IP Matrix")
  image(B, main="Derivative-of-Basis IP Matrix")
  
}

test.3 <- function(X, Y, seeds, bw, sliced=FALSE)
{
  # Computes the heat-equation based smoothing of a functional representation with dirichlet boundary conditions 
  # X ... NxP matrix of N-points in P dimensional space
  # Y ... Nx1 matrix of funtion observations at each of the P-dimensional points in X 
  # SEEDS ... MxK matrix of M-points in P dimensional space, seeds or knot-points of linear piecewise basis 
  # BW ... smoothing parameter, eqiuvalent to bandwidth of gaussian kernel 
  # SLICED ... (applicable for 3D) integer if seeds are equispaced, defines the extra dimensional slice for plotting data 
  
  n <- nrow(X); p <- ncol(X)
  basis_obj <- create.polynomial.basis(seeds)
  if (p==2) {
    trimesh(basis_obj$elements, basis_obj$seeds, main="Tesselation")
    title(main="Tesselation")
  }
  representation <- representation.basis(basis_obj, X, Y) 
  
  smoothed_trajectories <- smooth.heat(representation$repr_basis, bw=bw, dt=0.01)
  
  trajectories <- smoothed_trajectories$coef_trajectories
  niter <- smoothed_trajectories$smooth_params$niter
  dt <- smoothed_trajectories$smooth_params$dt
  
  t1 <- 0; t2 <- floor((niter+1)/4)*dt; t3 <- dt*niter
  
  if (p==2) {
    interpolation1 <- interp(basis_obj$seeds[,1], basis_obj$seeds[,2],  trajectories[, 1])
    interpolation2 <- interp(basis_obj$seeds[,1], basis_obj$seeds[,2],  trajectories[, floor((niter+1)/4)])
    interpolation3 <- interp(basis_obj$seeds[,1], basis_obj$seeds[,2],  trajectories[, niter+1])
    
    filled.contour(interpolation1$x, interpolation1$y, interpolation1$z,
                   color.palette=colorRampPalette(c("lightyellow", "orange", "red"), bias=2),
                   xlab="X", ylab="Y",
                   main=paste0("tf = ", t1),
                   key.title = title(main="Value", cex.main=0.8)
    )
    filled.contour(interpolation2$x, interpolation2$y, interpolation2$z,
                   color.palette=colorRampPalette(c("lightyellow", "orange", "red"), bias=2),
                   xlab="X", ylab="Y",
                   main=paste0("tf = ", t2),
                   key.title = title(main="Value", cex.main=0.8)
    )
    filled.contour(interpolation3$x, interpolation3$y, interpolation3$z,
                   color.palette=colorRampPalette(c("lightyellow", "orange", "red"), bias=2),
                   xlab="X", ylab="Y",
                   main=paste0("tf = ", t3),
                   key.title = title(main="Value", cex.main=0.8)
    )
  } else if (p==3 && sliced!=FALSE) {
  
    zslice <- basis_obj$seeds[,3]==sliced

    interpolation1 <- interp(basis_obj$seeds[zslice,1], basis_obj$seeds[zslice,2],  trajectories[zslice, 1])
    interpolation2 <- interp(basis_obj$seeds[zslice,1], basis_obj$seeds[zslice,2],  trajectories[zslice, floor((niter+1)/4)])
    interpolation3 <- interp(basis_obj$seeds[zslice,1], basis_obj$seeds[zslice,2],  trajectories[zslice, niter+1])
    
    filled.contour(interpolation1$x, interpolation1$y, interpolation1$z,
                   color.palette=colorRampPalette(c("lightyellow", "orange", "red"), bias=2),
                   xlab="X", ylab="Y",
                   main=paste0("tf = ", t1),
                   key.title = title(main="Value", cex.main=0.8)
    )
    filled.contour(interpolation2$x, interpolation2$y, interpolation2$z,
                   color.palette=colorRampPalette(c("lightyellow", "orange", "red"), bias=2),
                   xlab="X", ylab="Y",
                   main=paste0("tf = ", t2),
                   key.title = title(main="Value", cex.main=0.8)
    )
    filled.contour(interpolation3$x, interpolation3$y, interpolation3$z,
                   color.palette=colorRampPalette(c("lightyellow", "orange", "red"), bias=2),
                   xlab="X", ylab="Y",
                   main=paste0("tf = ", t3),
                   key.title = title(main="Value", cex.main=0.8)
    )
  }
  
}

get.data <- function(d=2, res=100) 
{
  # Dataset builder for testing FEFA functionality 
  # Arguments 
  # D ... parameter space dimension of the data 
  # RES ... number of equispaced points in each direction of parameter space 
  # Returns 
  # X ... nxd dimensional matrix of parameter vectors for n observations 
  # Y ... nx1 dimensional matrix of function evaluations 
  
  if (d==2) {
    # Construct sample data, consider adding noise 
    x1 <- x2 <- seq(-4*pi, 4*pi, len = res)
    X <- as.matrix(expand.grid(X1=x1, X2=x2))
    Y <- cos(X[, 1]^2+X[, 2]^2)*exp(-1*sqrt(X[, 1]^2+X[, 2]^2)/6)
  } else if (d==3) {
    x1 <- x2 <- x3 <- seq(-4*pi, 4*pi, len = res)
    X <- as.matrix(expand.grid(X1=x1, X2=x2, X3=x3))
    Y <- cos(X[, 1]^2+X[, 2]^2+X[, 3]^2)*exp(-1*sqrt(X[, 1]^2+X[, 2]^2+X[, 3]^2)/6)
  } else {
    stop(paste0("d = ", d, " not yet implemented!"))
  }
  return(list(X=X, Y=Y))
}

plant.seeds <- function(d=2, rand=FALSE, res=100, fac=2, covar=-0.5, n=300) 
{
  # Create seed set for placing basis functions 
  # Arguments 
  # D ... dimension of paramter space 
  # RAND ... if TRUE, seeds placed according to multivariate normal distribution centered at origin
  # FAC ... if RAND==FALSE, knot spacing length
  # COVAR ... if RAND==TRUE, covariance of multivariate normal from which seeds are sampled
  # N ... if RAND==TRUE, number of seeds to be sampled from multivariate normal distribution
  # Returns 
  # SEEDS ... nxd matrix of seed locations
  if (rand) {
    if (d==2) {
      mu <- c(0.0, 0.0)
      std <- 4*pi
      sigma <- std*rbind(c(1, covar), c(covar, 1))
      samples <- mvrnorm(n, mu=mu, Sigma=sigma)
      vertex_set <- rbind(c(4*pi, 4*pi), 
                          c(4*pi, -4*pi), 
                          c(-4*pi, 4*pi), 
                          c(-4*pi, -4*pi))
      seeds <- rbind(samples, vertex_set) 
    } else if (d==3) {
      mu <- c(0.0, 0.0, 0.0)
      std <- 4*pi
      sigma <- std*rbind(c(1, covar, covar), c(covar, 1, covar), c(covar, covar, 1))
      samples <- mvrnorm(n, mu=mu, Sigma=sigma)
      vertex_set <- 4*pi*rbind(c(1, 1, 1), 
                          c(1, 1, -1), 
                          c(1, -1, 1),
                          c(-1, 1, 1), 
                          c(1, -1, -1),
                          c(-1, 1, -1),
                          c(-1, -1, 1),
                          c(-1, -1, -1))
      seeds <- rbind(samples, vertex_set) 
    } else {
      stop(paste0("d = ", d, " not yet implemented!"))
    } 
  } else {
    if (d==2) {
      x1s <- seq(-4*pi, 4*pi, len = res/fac)
      x2s <- seq(-4*pi, 4*pi, len = res/fac)
      seeds <- as.matrix(expand.grid(X1=x1s, X2=x2s))
    } else if (d==3) {
      x1s <- seq(-4*pi, 4*pi, len = res/fac)
      x2s <- seq(-4*pi, 4*pi, len = res/fac)
      x3s <- seq(-4*pi, 4*pi, len = res/fac)
      seeds <- as.matrix(expand.grid(X1=x1s, X2=x2s, X3=x3s))
    } else {
      stop(paste0("d = ", d, " not yet implemented!"))
    }
  }
  return(seeds)
}

main <- function(CASE=-1, d=2, res=100, bw=8) 
{
  # Central controller for tests 
  
  data <- get.data(d=d, res=res)
  X <- data$X; Y <- data$Y
  
  if (d==2) {
    #Z <- matrix(Y, nrow=length(x1), 
    #            ncol=length(x2))
    
    #image.plot(Z, axes = FALSE, main = "Original Function",
    #           xlab = expression(cos(r^2) * e^{-r/6}))
    
    interpolation <- interp(X[,1], X[,2], Y)
    
    filled.contour(interpolation$x, interpolation$y, interpolation$z,
                   color.palette=colorRampPalette(c("lightyellow", "orange", "red"), bias=2),
                   xlab=paste0("dimension ", 1), ylab=paste0("dimension ", 2),
                   main="Linear interpolation of image",
                   key.title = title(main="Value", cex.main=0.8))
    }
  
  seeds_equi <- plant.seeds(d=d, rand=FALSE, res=res, fac=2)
  seeds_rand <- plant.seeds(d=d, rand=TRUE, covar=-0.5, n=300)
  
  seeds <- seeds_equi
  
  knot_points <- unique(seeds[,1]) #equivalent in all directions when seeds <- seeds_equi
  sliced <- min(knot_points[knot_points > 0])
  
  if (CASE==1) {
    test.1(X, Y, seeds, nderiv=0, sliced=sliced)
  } else if (CASE==2) {
    test.2(seeds)
  } else if (CASE==3) {
    test.3(X, Y, seeds, bw, sliced=sliced)
  }
  
}

setwd("/axiom/home/wconsagra/Workspace/git/FEFA")
source("init.R")

set.seed(111)

main(CASE=3, d=3, res=20, bw=8)



