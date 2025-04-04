library(Rcpp)
library(foreach)
library(doParallel)
library(ggplot2)
library(ggridges)

Rcpp::sourceCpp("./summary.split.variable.cpp")

plot.bart.variable <- function(bart, sum.by.tree=F) {
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
  if (sum.by.tree) {
    plt.var <- ggplot(res, aes(x = depth, y = variable, fill = variable)) +
      geom_density_ridges(stat='binline', alpha = 0.5, draw_baseline=T, binwidth=1, scale=0.8) +
      labs(fill = "Variable") + ggtitle('Variables of whole tree')+
      theme_minimal()
  } else {
    plt.var <- ggplot(res, aes(x = depth, y = variable, fill = variable)) +
      geom_density_ridges(stat='binline', alpha = 0.5, draw_baseline=T, binwidth=1, scale=0.8) +
      facet_wrap(~as.factor(k)) +
      labs(fill = "Variable") + ggtitle('Variables of each tree')+
      theme_minimal()
  }
  
  print(plt.var)
}

#plot.bart.variable(mc.result, T)

