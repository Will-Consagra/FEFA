setwd("/axiom/home/wconsagra/Workspace/git/FEFA")
source("init.R")

set.seed(111)

## Testing construction of Phi and DPhi 
res <- 100
x1 <- x2 <- seq(-4*pi, 4*pi, len = res)
X <- as.matrix(expand.grid(X1=x1, X2=x2))
Y <- cos(X[, 1]^2+X[, 2]^2)*exp(-1*sqrt(X[, 1]^2+X[, 2]^2)/6)
Z <- matrix(Y, nrow=length(x1), 
               ncol=length(x2))
#image(z = Z, col  = gray((0:32)/32))
image(Z, axes = FALSE, main = "Math can be beautiful ...",
      xlab = expression(cos(r^2) * e^{-r/6}))

fac <- 4
x1s <- seq(-4*pi, 4*pi, len = res/fac)
x2s <- seq(-4*pi, 4*pi, len = res/fac)

seeds <- as.matrix(expand.grid(X1=x1s, 
                               X2=x2s))
basis_obj <- create.polynomial.basis(seeds)
#trimesh(basis_obj$elements, basis_obj$seeds)

Phi <- eval.basis(X, basis_obj, nderiv=0) 
## Little projection code, must formalize!
betahat <- solve(t(Phi)%*%Phi)%*%t(Phi)%*%Y

Zhat <- matrix(betahat, nrow=length(x1s), 
               ncol=length(x2s))
#image(z = Zhat, col  = gray((0:32)/32))
image(Zhat, axes = FALSE, main = "Fac = 4 ...",
      xlab = "Zhat")

fac <- 2
x1s <- seq(-4*pi, 4*pi, len = res/fac)
x2s <- seq(-4*pi, 4*pi, len = res/fac)

seeds <- as.matrix(expand.grid(X1=x1s, 
                               X2=x2s))
basis_obj <- create.polynomial.basis(seeds)
#trimesh(basis_obj$elements, basis_obj$seeds)

Phi <- eval.basis(X, basis_obj, nderiv=0) 
## Little projection code, must formalize!
betahat <- solve(t(Phi)%*%Phi)%*%t(Phi)%*%Y

Zhat <- matrix(betahat, nrow=length(x1s), 
               ncol=length(x2s))
#image(z = Zhat, col  = gray((0:32)/32))
image(Zhat, axes = FALSE, main = "Fac = 2 ...",
      xlab = "Zhat")


fac <- 1
x1s <- seq(-4*pi, 4*pi, len = res/fac)
x2s <- seq(-4*pi, 4*pi, len = res/fac)

seeds <- as.matrix(expand.grid(X1=x1s, 
                               X2=x2s))
basis_obj <- create.polynomial.basis(seeds)
#trimesh(basis_obj$elements, basis_obj$seeds)

Phi <- eval.basis(X, basis_obj, nderiv=0) 
## Little projection code, must formalize!
betahat <- solve(t(Phi)%*%Phi)%*%t(Phi)%*%Y

Zhat <- matrix(betahat, nrow=length(x1s), 
               ncol=length(x2s))
#image(z = Zhat, col  = gray((0:32)/32))
image(Zhat, axes = FALSE, main = "Fac = 1 ...",
      xlab = "Zhat")
