library(Rcpp)
library(foreach)
library(doParallel)

Rcpp::sourceCpp("./summary.split.variable.cpp")

fetch.bart.variable <- function(bart) {
  tree.chain <- bart$trees
  chain.len <- length(tree.chain)
  m <- length(tree.chain[[1]])
  
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  res <- foreach(iter=icount(chain.len), .combine = rbind) %:%
    foreach(k=icount(m), .combine = rbind) %do% {
      tree.k <- tree.chain[[iter]][[k]]
      data.frame(k=k, C_summary_tree(tree.chain[[1]][[1]]))
    }
  stopCluster(cl)
  return(res)
}

res <- fetch.bart.variable(mc.result)
