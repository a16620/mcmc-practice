library(Rcpp)
library(foreach)
library(doParallel)

library(ggplot2)
library(ggridges)
library(ggpubr)

cpp.func <- new.env()
Rcpp::sourceCpp("./summary.split.variable.cpp", env = cpp.func)

fetch.bart.variable <- function(bart) {
  tree.chain <- bart$trees
  chain.len <- length(tree.chain)
  m <- length(tree.chain[[1]])
  
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  res <- foreach(iter=icount(chain.len), .combine = rbind) %:%
    foreach(k=icount(m), .combine = rbind) %do% {
      tree.k <- tree.chain[[iter]][[k]]
      data.frame(k=k, cpp.func$C_summary_tree(tree.k))
    }
  stopCluster(cl)
  return(res)
}

plot.bart.variable <- function(bart, sum.by.tree=F) {
  res <- fetch.bart.variable(bart)
  if (sum.by.tree) {
    plt.var <- ggplot(res, aes(x = depth, y = variable, fill = variable)) +
      geom_density_ridges(stat='binline', alpha = 0.5, draw_baseline=T, binwidth=1, scale=0.8) +
      labs(fill = "Variable") + ggtitle('Variables of whole tree')+
      theme_minimal()
    
    df <- res %>% group_by(variable) %>% summarise(freq=n(), .groups = 'drop')
    plt.var2 <- ggplot(df) + geom_bar(aes(x=variable, y=freq, fill=variable), stat='identity') + facet_wrap(~k) +
      theme_minimal()
  } else {
    plt.var <- ggplot(res, aes(x = depth, y = variable, fill = variable)) +
      geom_density_ridges(stat='binline', alpha = 0.5, draw_baseline=T, binwidth=1, scale=0.8) +
      facet_wrap(~as.factor(k)) +
      labs(fill = "Variable") + ggtitle('Variables of each tree')+
      theme_minimal()
    
    df <- res %>% group_by(k, variable) %>% summarise(freq=n(), .groups = 'drop')
    plt.var2 <- ggplot(df) + geom_bar(aes(x=variable, y=freq, fill=variable), stat='identity') + facet_wrap(~k) +
      theme_minimal()
  }
  
  print(ggarrange(plt.var, plt.var2))
}

