
equi.seeder <- function(V, r=5)
  
{
  # This function generates an equi-spaced lattice seed set over a rectangular domain 
  # Arguments 
  # V ... d x 2 dimensional matrix, each row of the form (min(x_i), max(x_i)), for each coordinate x_i
  # R ... number of seeds per dimension, i.e. marginal basis size 
  # Returns 
  # SEEDS ... R^d x d dimensional matrix of seed locations 
  
  if (!is.matrix(V)) {
    stop("'V' must be of type 'matrix'")
  }

  if (!is.numeric(r)) {
    stop("'r' must of of type 'numeric'")
  }
  
  X <- apply(V, 1, function(row){return(seq(row[1], row[2], len=r))})
  seeds <- as.matrix(expand.grid(lapply(seq_len(ncol(X)), function(j) c(X[,j]))))
  return(seeds)
    
}

runif.seeder <- function(V, N, repl=1, crit="max_min", include_boundary=TRUE) 
  
{
  # Generate a set of seeds sampled according to a d-dimensional uniform distribution over the rectange {[a_i, b_i]}_{i=1}^d
  # Arguments 
  # V ... d x 2 dimensional matrix, each row of the form (min(x_i), max(x_i)), for each coordinate x_i; i.e. bounds of marginal uniform distribution
  # N ... number of seeds to sow 
  # REPL ... number of replications 
  # CRIT ... selection criteria for replicants; can be one of the following 
  #           "max_min" ... maximum minimall pairwise distance between seeds
  #           "var" ... maximum sample variance among replicant's pairwise distances 
  # INCLUDE_BOUNDARY ... if TRUE, add the vertices of the hyperrectangle defined by V to the seed set. Note, this will have the effect of
  #                       generating N - 2^d random seed points
  # Returns 
  # SEEDS ... N x d dimensional matrix of seed locations 
  # note: probably want to insitute a minimum allowable distance between seeds to avoid numerical issues downstream: consider repulsive spring 
  #       force attached to each point, as modeled into tSNE 
  
  # validate input 
  
  if (!is.matrix(V)) {
    stop("'V' must be of type 'matrix'")
  }
  
  d <- nrow(V)
  
  if (!is.numeric(N)) {
    stop("'N' must of of type 'numeric'")
  }
  
  if (include_boundary) {
    N <- N - 2^d
  }
  
  # build seeds sample 
  
  replicants <- array(0, c(repl, N, d))
  D <-matrix(0, nrow=N*(N-1)/2, ncol=repl)
  for (r in 1:repl) {
    seeds <- apply(V, 1, function(row){return(runif(N, min=row[1], max=row[2]))})
    replicants[r,,] <- seeds
    D[, r] <- c(dist(seeds, method="euclidean"))
  }

  # select replicant based on previously defined criteria

  if (crit=="max_min") {
    ix <- which.max(apply(D, 2, min))
  } else if (crit=="var") {
    ix <- which.max(apply(D, 2, var))
  } else {
    stop(paste0("Selection criteria '", crit, "' not yet implemented"))
  }
  
  seeds <- replicants[ix,,]
  
  if (include_boundary) {
    
    helper1 <- function(row) {
      result <- c()
      for (i in 1:length(row)) {
        result <- c(result, V[i, row[i]])
      }
      return(result)
    }
    
    cix <- permutations(n=2, r=d, v=c(1,2), repeats.allowed = TRUE)
    vertex_set <- t(apply(cix, 1, helper1))
    seeds <- rbind(seeds, 
                   vertex_set)
    }

  return(seeds)
  
}
