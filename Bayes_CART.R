library(data.tree)
library(rpart)
library(dplyr)
library(ggplot2)
library(ggpubr)

##################### Create #####################

BCART <- function(formula., data, method=NULL) {
  combine.df <- model.frame(formula., data)
  y.is.factor <- is.factor(combine.df[,1])
  
  method.list <- c("category", "normal", "normal2", "poisson")
  if (is.null(method)) {
    if (y.is.factor)
      method <- "category"
    else
      method <- "normal"
  }
  
  if (!(method %in% method.list)) {
    cat("method must be (", paste(method.list, collapse = ', '), ")\n", sep='')
    return(NULL)
  }
  
  cat(paste0("Model: ", format(formula.)), " method=", method, sep='')
  cat("\nPredictors: ", paste(colnames(combine.df)[-1], collapse=', '), '\n')
  
  model <- list()
  model$formula <- formula.
  model$method <- method
  model$tree <- CART.set.obs(Node$new("Root", full.obs=combine.df), 1:nrow(combine.df))
  
  return(model)
}

##################### MCMC #####################

extract.marginal.params <- function(model, uniform.dirichlet=0) {
  method <- model$method
  data <- model$tree$full.obs
  
  if (method == "category") {
    Y <- data[, all.vars(model$formula)[1]]
    if (uniform.dirichlet != 0) {
      a <- rep(uniform.dirichlet, length(levels(Y)))
    } else {
      a <- as.numeric(table(Y))
    }
    
    return(list(
      a=a,
      a.s=sum(a),
      a.p=lgamma(sum(a))-sum(lgamma(a)),
      y.cont=F
    ))
  } else if (method == "normal") {
    fit.tree <- rpart(model$formula, data=data, method="anova")
    var.est <- var(fit.tree$y)
    var.min <- data.frame(leaf=fit.tree$where, y=fit.tree$y) %>% group_by(leaf) %>% summarise(var=var(y)) %>% pull(var) %>% min()
    return(list(
      mu=mean(fit.tree$y),
      a=var.min/var.est,
      lambda=var.min,
      nu=length(unique(fit.tree$where)),
      y.cont=T
    ))
  }
}

MCMC.param <- function(model, marginal.only=T) {
  mcmc.param <- list()
  
  mcmc.param$split.prob <- function(level) {
    0.95*level**(-0.5)
  }
  
  marginal.param <- extract.marginal.params(model)
  if (model$method == "category")
    mcmc.param$tree.marginal <- function(tree) CART.prob.likelihood.marginal.category(tree, marginal.param)
  else if (model$method == "normal")
    mcmc.param$tree.marginal <- function(tree) CART.prob.likelihood.marginal.regression(tree, marginal.param)
  else if (model$method == "normal2")
    mcmc.param$tree.marginal <- function(tree) CART.prob.likelihood.marginal.regression2(tree, marginal.param)
  else if (model$method == "poisson")
    mcmc.param$tree.marginal <- NULL
  
  mcmc.param$use.latent <- !(marginal.only || model$method %in% c("category", "normal", "normal2"))

  return(mcmc.param)
}

do.MCMC <- function(model0, param=NULL, iteration=1000, burn.in=0,
                    moves=c("grow", "prune", "change", "swap"), prob.moves=NULL,
                    verbose=F, plot.trace=T, additional.criteria=NULL, additional.criteria.fun=NULL, model.select.criteria=NULL) {
  stopifnot(length(moves) > 0)
  stopifnot(length(additional.criteria) == length(additional.criteria.fun))

  if (is.null(param)) {
    param <- MCMC.param(model0)
  }
  
  if (model0$method == "category") {
    crit.loss <- "miscl"
    crit.loss.fn <- CART.get.tree.train.miscl
  } else if (model0$method %in% c("normal", "normal2")) {
    crit.loss <- "sse"
    crit.loss.fn <- CART.get.tree.train.SSE
  } else {
    crit.loss <- 'DIC'
    crit.loss.fn <- NULL
  }
  
  tree.marginal.lik <- param$tree.marginal
  
  criteria.name <- append(c("log.post", "n.leaf", crit.loss), additional.criteria)
  criteria.funs <- setNames(append(list(
    function(tree) {tree.marginal.lik(tree)+CART.prob.prior(tree, param$split.prob)},
    function(tree) {tree$leafCount},
    crit.loss.fn
  ), additional.criteria.fun), criteria.name)
  
  if (is.null(model.select.criteria))
    model.select.criteria <- crit.loss
  
  mtree <- Clone(model0$tree)
  #burn in
  if (burn.in > 0) {
    if (verbose) {
      cat("burn-in...")
    }
    burn.in.l.post <- tree.marginal.lik(mtree)
    for (iter in 1:burn.in) {
      mc.move <- sample(moves, 1, prob = prob.moves)
      if (mc.move == "grow") {
        mc.new <- CART.move.grow(mtree, param$split.prob)
      } else if (mc.move == "prune") {
        mc.new <- CART.move.prune(mtree, param$split.prob)
      } else if (mc.move == "change") {
        mc.new <- CART.move.change(mtree, param$split.prob, only.value = (runif(1) < 0.5))
      } else if (mc.move == "swap") {
        mc.new <- CART.move.swap(mtree, param$split.prob)
      }
      
      if (is.null(mc.new) || !CART.check.tree.ok(mc.new$tree.new)) {
        next
      }
      
      burn.in.l.post.new <- tree.marginal.lik(mc.new$tree.new)
      l.a.ratio.u <- mc.new$l.prob.rev+burn.in.l.post.new
      l.a.ratio.l <- mc.new$l.prob+burn.in.l.post
      
      if (log(runif(1)) <= (l.a.ratio.u-l.a.ratio.l)) {
        mtree <- mc.new$tree.new
        burn.in.l.post <- burn.in.l.post.new
      }
    }
    if (verbose) {
      cat("end!\n")
    }
  }
  
  #with record
  criteria.matrix <- `colnames<-`(matrix(0, ncol=length(criteria.name), nrow=iteration), criteria.name)
  criteria.matrix[1,] <- sapply(criteria.funs, function(cfn) {cfn(mtree)})
  
  selected.model <- mtree
  selected.model.score <- criteria.matrix[1, model.select.criteria]
  
  print.intval <- as.integer(iteration*0.1)
  
  accept.count <- setNames(rep(0, length(moves)), moves)
  fail.count <- setNames(rep(0, length(moves)), moves)
  for (iter in 2:iteration) {
    mc.move <- sample(moves, 1, prob = prob.moves)
    if (mc.move == "grow") {
      mc.new <- CART.move.grow(mtree, param$split.prob)
    } else if (mc.move == "prune") {
      mc.new <- CART.move.prune(mtree, param$split.prob)
    } else if (mc.move == "change") {
      mc.new <- CART.move.change(mtree, param$split.prob, only.value = (runif(1) < 0.5))
    } else if (mc.move == "swap") {
      mc.new <- CART.move.swap(mtree, param$split.prob)
    }
    
    if (is.null(mc.new) || !CART.check.tree.ok(mc.new$tree.new)) {
      fail.count[mc.move] <- fail.count[mc.move]+1
      criteria.new.tree <- criteria.matrix[iter-1,]
    } else {
      criteria.new.tree <- sapply(criteria.funs, function(cfn) {cfn(mtree)})
      
      l.a.ratio.u <- mc.new$l.prob.rev+criteria.new.tree[1]
      l.a.ratio.l <- mc.new$l.prob+criteria.matrix[iter-1,1]
      
      if (log(runif(1)) <= (l.a.ratio.u-l.a.ratio.l)) {
        mtree <- mc.new$tree.new
        accept.count[mc.move] <- accept.count[mc.move] + 1
      } else {
        criteria.new.tree <- criteria.matrix[iter-1,]
      }
    }
    criteria.matrix[iter,] <- criteria.new.tree
    
    if (selected.model.score >= criteria.new.tree[model.select.criteria]) {
      selected.model <- mtree
      selected.model.score <- criteria.new.tree[model.select.criteria]
    }
    
    if (verbose && iter %% print.intval == 0) {
      cat('============= #', as.character(iter), '(', format(sum(accept.count)/iter, digits = 3),')\n', sep='')
      print(colMeans(criteria.matrix[(iter-print.intval+1):(iter),]), quot=F)
    }
  }
  
  criteria.df <- as.data.frame(criteria.matrix)
  if (plot.trace) {
    ggarrange(ggplot(criteria.df) + geom_line(aes(x=1:iteration, y=log.post))+xlab('iter')+theme_classic(),
              ggplot() + geom_line(aes(x=1:iteration, y=criteria.df[,crit.loss]))+xlab('iter')+ylab(crit.loss)+theme_classic(),
              ggplot(criteria.df) + geom_line(aes(x=1:iteration, y=n.leaf))+xlab('iter')+theme_classic(),
              ncol = 1)
  }
  model0$tree <- selected.model
  return(list(
    history=criteria.df,
    error.move=fail.count,
    accept.move=accept.count,
    accept.rate=sum(accept.count)/iteration,
    model=model0,
    score=selected.model.score
  ))
}

##################### Deploy #####################

deploy <- function(model) {
  tree <- model$tree
  if (is.null(tree$full.obs))
    return()
  CART.set.predict(tree)
  CART.clean.obs(tree)
}

Update.prediction <- function(model, new.data) {
  tree <- model$tree
  tree$full.obs <- model.frame(model$formula, new.data)
  CART.set.obs(tree, 1:nrow(new.data))
  CART.update.obs(tree)
  if (!CART.check.tree.ok(tree))
    return(F)
  
  CART.set.predict(tree)
  CART.clean.obs(tree) 
  
  return(T)
}

##################### Probability #####################

CART.prob.likelihood.marginal.regression <- function(tree, params, verb=F) {
  c.a <- params$a
  c.b <- tree$leafCount
  c.c <- ifelse(is.null(params$c), 0, params$c)
  c.mu <- params$mu
  c.nu <- params$nu
  c.lambda <- params$lambda
  
  l.prob <- c.c+0.5*c.b*log(c.a)-0.5*sum(sapply(tree$leaves, function(node) log(length(node$obs.idx))))-
    0.5*(nrow(tree$full.obs)+c.nu)*log(sum(sapply(tree$leaves, function(node) {
      node.obs <- CART.get.obs(node)
      
      n.i <- nrow(node.obs)
      y.i <- node.obs[,1]
      s.i <- CART.get.node.train.SSE(node)
      t.i <- (mean(y.i)-c.mu)**2*c.a*n.i/(c.a+n.i)
      return(t.i+s.i)
    }))+c.nu*c.lambda)
  return(l.prob)
}

#V2
CART.prob.likelihood.marginal.regression2 <- function(tree, params) {
  c.a <- params$a
  c.c <- ifelse(is.null(params$c), 0, params$c)
  c.mu <- params$mu
  c.nu <- params$nu
  c.nuMlamb <- params$lambda*c.nu
  c.aMnu <- 0.5*log(c.a)-lgamma(c.nu/2)
  
  l.prob <- c.c+0.5*c.nu*log(c.nuMlamb)+
    sum(sapply(tree$leaves, function(node) {
      n.i <- length(node$obs.idx)
      c.aMnu-0.5*log(n.i+c.a)+lgamma((n.i+c.nu)/2)
    }))-
    0.5*sum(sapply(tree$leaves, function(node) {
      node.obs <- CART.get.obs(node)
      n.i <- nrow(node.obs)
      y.i <- node.obs[,1]
      s.i <- CART.get.node.train.SSE(node)
      t.i <- (mean(y.i)-c.mu)**2*c.a*n.i/(c.a+n.i)
      return((n.i+c.nu)*log(t.i+s.i+c.nuMlamb))
    }))
  return(l.prob)
}

CART.prob.likelihood.marginal.category <- function(tree, params, verb=F) {
  c.a <- params$a
  c.as <- params$a.s
  c.ap <- params$a.p
  l.prob <- sum(sapply(tree$leaves, function(node) {
    node.obs <- CART.get.obs(node)
    n.i <- as.numeric(table(node.obs[,1]))
    if (verb) {
      print(n.i)
    }
    return(c.ap+sum(lgamma(n.i+c.a))-lgamma(sum(n.i)+c.as))
  }))
  return(l.prob)
}

CART.prob.prior <- function(tree, split.prob) {
  node.info <- tree$Get(function(node) c(node$level, isLeaf(node)), simplify = T)
  node.prob <- split.prob(node.info[1,])
  sum(log(ifelse(node.info[2,], 1-node.prob, node.prob)))
}

CART.prob.select.rule <- function(obs, len.rule.values, is.categorical) {
  if (len.rule.values <= 1 || is.na(len.rule.values)) {
    return(0)
  }
  if (is.categorical) {
    return(-log(2**len.rule.values-2)-log(ncol(obs)-1))
  } else {
    return(-log(len.rule.values)-log(ncol(obs)-1))
  }
}

CART.prob.rule <- function(node) {
  obs <- CART.get.obs(node)
  rule.colname <- node$split.rule$rule.col
  if (is.factor(obs[,rule.colname])) {
    len.rule.values <- length(unique(obs[,rule.colname]))
    return(CART.prob.select.rule(obs, len.rule.values, T))
  } else {
    return(CART.prob.select.rule(obs, nrow(obs)-1, F))
  }
}

##################### Predict #####################

CART.set.predict <- function(tree) {
  if (is.factor(tree$full.obs[,1])) {
    Do(tree$leaves, function(node) {
      tb.val <- table(CART.get.obs(node)[,1])
      node$pred.val <- names(which.max(tb.val))
    })
  } else {
    Do(tree$leaves, function(node) {
      node$pred.val <- mean(CART.get.obs(node)[,1])
    })
  }
}

CART.get.predict <- function(tree, data) {
  pred.val <- numeric(nrow(data))
  for (i in 1:nrow(data)) {
    row <- data[i,,drop=F]
    node <- tree
    while(isNotLeaf(node)) {
      if (node$split.rule$fun(row)) {
        node <- node$children[[1]]
      } else {
        node <- node$children[[2]]
      }
    }
    pred.val[i] <- node$pred.val
  }
  return(pred.val)
}

CART.get.tree.train.SSE <- function(tree) {
  sum(sapply(tree$leaves, CART.get.node.train.SSE))
}

CART.get.node.train.SSE <- function(node) {
  obs.idx <- node$obs.idx
  if (length(obs.idx) == 0) {
    return(0)
  } else {
    y <- CART.get.obs(node)[,1]
    if (length(y) > 1) {
      return(var(y)*(length(y)-1))
    }
    return(0)
  }
}

CART.get.tree.train.miscl <- function(tree) {
  n.total <- nrow(tree$full.obs)
  sum(sapply(tree$leaves, function(node) {
    CART.get.node.train.miscl(node)*length(node$obs.idx)/n.total
  }))
}

CART.get.node.train.miscl <- function(node) {
  obs.idx <- node$obs.idx
  if (length(obs.idx) == 0) {
    return(0)
  } else {
    y <- CART.get.obs(node)[,1]
    p <- as.numeric(table(y))
    p <- p/sum(p)
    return(1-max(p))
  }
}

##################### Utility #####################

CART.check.tree.ok <- function(tree) {
  all(sapply(tree$leaves, function(node) {length(node$obs.idx) > 0}))
}

CART.rule2name <- function(obs, split.col, split.value) {
  obs.filt.col <- obs[,split.col]
  if (is.factor(obs.filt.col)) {
    lvs <- levels(obs.filt.col)
    lvs.key <- lvs %in% split.value
    subl <- lvs[lvs.key]
    subr <- lvs[!lvs.key]
    
    return(c(
      paste0(split.col, "={", paste(subl,collapse=","), '}'),
      paste0(split.col, "={", paste(subr,collapse=","), '}')
    ))
  } else {
    return(c(
      paste0(split.col, "<=", as.character(split.value)),
      paste0(split.col, ">", as.character(split.value))
    ))
  }
}

CART.get.obs <- function(node) {
  return(node$root$full.obs[node$obs.idx, , drop=F])
}

CART.set.obs <- function(node, obs.idx) {
  node$obs.idx <- obs.idx
  return(node)
}

CART.clean.obs <- function(tree) {
  tree$RemoveAttribute('full.obs')
  nodes <- Traverse(tree)
  Do(nodes, function(node) {
    node$RemoveAttribute('obs.idx')
  })
}

as.probability <- function(unp) {
  unp/sum(unp)
}

##################### Grow #####################

CART.select.rule <- function(node, col.candidates=NULL) {
  obs <- CART.get.obs(node)
  col.candidates <- if (is.null(col.candidates)) colnames(obs)[-1] else col.candidates
  while (length(col.candidates) > 0) {
    rule.colname <- sample(col.candidates, 1)
    if (is.factor(obs[,rule.colname])) {
      rule.values <- unique(obs[,rule.colname])
      len.values <- length(rule.values)
      if (len.values < 2) {
        col.candidates <- col.candidates[col.candidates != rule.colname]
        next
      }
      rule.value <- rule.values[sample(1:len.values, sample(1:(len.values-1), 1))]
      
      return(
        list(
          split.col=rule.colname,
          split.value=sort(rule.value)
        )
      )
    } else {
      rule.values <- sort(obs[,rule.colname])
      len.rule.values <- length(rule.values)
      
      if (len.rule.values <= 2 || rule.values[1] == rule.values[len.rule.values]) {
        col.candidates <- col.candidates[col.candidates != rule.colname]
        next
      }
      
      val.idx <- sample(1:(length(rule.values)-1), 1)
      val.aux <- 
      rule.value <- runif(1, min=rule.values[val.idx], max=rule.values[val.idx+1])
      if (rule.value == rule.values[val.idx+1]) {
        rule.value <- (rule.values[val.idx]+rule.values[val.idx+1])/2
      }
      return(
        list(
          split.col=rule.colname,
          split.value=rule.value
        )
      )
    }
  }
  return(NULL)
}

CART.set.rule <- function(node, rule.info) {
  rule <- function(obs) {
    obs.filt.col <- obs[,rule.info$split.col]
    if (is.factor(obs.filt.col)) {
      return(obs.filt.col %in% rule.info$split.value)
    }
    return(obs.filt.col <= rule.info$split.value)
  }
  
  obs <- CART.get.obs(node)
  obs.key <- rule(obs)
  if (all(obs.key) || !any(obs.key))
    return(F)
  
  rule.info$fun <- rule
  node$split.rule <- rule.info
  left.obs.idx <- node$obs.idx[obs.key]
  right.obs.idx <- node$obs.idx[!obs.key]
  
  child.names <- CART.rule2name(obs, rule.info$split.col, rule.info$split.value)
  CART.set.obs(node$AddChild(child.names[1]), left.obs.idx)
  CART.set.obs(node$AddChild(child.names[2]), right.obs.idx)
  
  return(T)
}

##################### Change, Swap #####################

CART.update.obs <- function(tree) {
  tree$Do(function(node) {
    obs.key <- node$split.rule$fun(CART.get.obs(node))
    
    left.obs.idx <- node$obs.idx[obs.key]
    right.obs.idx <- node$obs.idx[!obs.key]
    
    CART.set.obs(node$children[[1]], left.obs.idx)
    CART.set.obs(node$children[[2]], right.obs.idx)
  }, traversal="level", filterFun = isNotLeaf)
  return(CART.check.tree.ok(tree))
}

CART.update.rule.swap <- function(node, rule.info) {
  node$split.rule <- rule.info
  
  child.names <- CART.rule2name(CART.get.obs(node), rule.info$split.col, rule.info$split.value)
  names(node$children) <- child.names
  node$children[[1]]$name <- child.names[1]
  node$children[[2]]$name <- child.names[2]
}

CART.update.rule <- function(node, rule.info) {
  rule <- function(obs) {
    obs.filt.col <- obs[,rule.info$split.col]
    if (is.factor(obs.filt.col)) {
      return(obs.filt.col %in% rule.info$split.value)
    }
    return(obs.filt.col <= rule.info$split.value)
  }
  
  obs <- CART.get.obs(node)
  obs.key <- rule(obs)
  if (all(obs.key) || !any(obs.key))
    return(F)
  
  rule.info$fun <- rule
  node$split.rule <- rule.info
  
  child.names <- CART.rule2name(obs, rule.info$split.col, rule.info$split.value)
  names(node$children) <- child.names
  node$children[[1]]$name <- child.names[1]
  node$children[[2]]$name <- child.names[2]
  return(T)
}

CART.compare.rule <- function(rule1, rule2) {
  if (rule1$split.col != rule2$split.col)
    return(F)
  if (length(rule1$split.value) != length(rule2$split.value))
    return(F)
  if (any(rule1$split.value != rule1$split.value))
    return(F)
  
  return(T)
}

##################### MCMC moves #####################

CART.move.grow <- function(tree, split.prob) {
  tree.new <- Clone(tree)
  
  terminal.nodes <- tree.new$leaves
  node.level <- sapply(terminal.nodes, function(node) node$level)
  node.prob <- as.probability(split.prob(node.level))
  
  while (length(terminal.nodes) > 0) {
    selected.node.idx <- sample(1:length(node.prob), 1, prob = node.prob)
    selected.node <- terminal.nodes[[selected.node.idx]]
    rule <- CART.select.rule(selected.node)
    if (is.null(rule)) {
      terminal.nodes <- terminal.nodes[-selected.node.idx]
      node.prob <- as.probability(node.prob[-selected.node.idx])
      next
    }
    if (CART.set.rule(selected.node, rule)) {
      break
    }
  }
  
  if (length(terminal.nodes) == 0) {
    return(NULL)
  }
  
  node.prob.rev <- as.probability(1/node.prob)
  return(list(
    tree.new=tree.new,
    l.prob=log(node.prob[selected.node.idx])+CART.prob.rule(selected.node),
    l.prob.rev=log(node.prob.rev[selected.node.idx])
  ))
}

CART.move.prune <- function(tree, split.prob) {
  if (tree$leafCount == 1) {
    return(NULL)
  }
  tree.new <- Clone(tree)
  
  parent.of <- Traverse(tree.new, filterFun = function(node) {
    isNotLeaf(node) && isLeaf(node$children[[1]]) && isLeaf(node$children[[2]])
  })
  node.level <- sapply(parent.of, function(node) node$level)
  node.prob.rev <- as.probability(split.prob(node.level))
  node.prob <- as.probability(1/node.prob.rev)
  
  selected.node.idx <- sample(1:length(node.prob), 1, prob = node.prob)
  selected.node <- parent.of[[selected.node.idx]]
  
  if (Prune(selected.node, function(cn) F) != 2)
    print("Warning: Pruned not 2")
  l.prob.sel.split <- CART.prob.rule(selected.node)
  selected.node$RemoveAttribute('split.rule')
  
  return(list(
    tree.new=tree.new,
    l.prob=log(node.prob[selected.node.idx]),
    l.prob.rev=log(node.prob.rev[selected.node.idx])+l.prob.sel.split
  ))
}

CART.move.change <- function(tree, split.prob, only.value=F) {
  if (tree$leafCount == 1) {
    return(NULL)
  }
  tree.new <- Clone(tree)
  split.node <- Traverse(tree.new, filterFun = isNotLeaf)
  node.level <- sapply(split.node, function(node) node$level)
  node.prob <- as.probability(1/split.prob(node.level))
  
  if (length(split.node) == 0)
    return(NULL)
  ok <- F
  for (max.try in 1:1000) {
    selected.node.idx <- sample(1:length(node.prob), 1)#, prob = node.prob)
    selected.node <- split.node[[selected.node.idx]]
    if (only.value)
      rule <- CART.select.rule(selected.node, selected.node$split.rule$rule.col)
    else
      rule <- CART.select.rule(selected.node)
    if (is.null(rule))
      next
    
    subtree <- Clone(selected.node)
    subtree$full.obs <- tree.new$full.obs
    
    if (!CART.update.rule(subtree, rule) || !CART.update.obs(subtree)) {
      next
    }
    ok <- T
    if (isNotRoot(selected.node)) {
      subtree$RemoveAttribute('full.obs')
      parent.of.no <- selected.node$parent
      parent.of.no$RemoveChild(selected.node$name)
      parent.of.no$AddChildNode(subtree)
    } else {
      tree.new <- subtree
    }
    
    break
  }
  if (!ok)
    return(NULL)
  return(list(
    tree.new=tree.new,
    l.prob=0,
    l.prob.rev=0
  ))
}

CART.move.swap <- function(tree, split.prob) {
  if (tree$leafCount < 3) {
    return(NULL)
  }
  tree.new <- Clone(tree)
  split.node <- Traverse(tree.new, filterFun = function(node) {
    isNotLeaf(node) && (isNotLeaf(node$children[[1]]) || isNotLeaf(node$children[[2]]))
  })
  node.level <- sapply(split.node, function(node) node$level)
  node.prob <- as.probability(1/split.prob(node.level))
  
  if (length(split.node) == 0)
    return(NULL)
  ok <- F
  for (max.try in 1:1000) {
    selected.node.idx <- sample(1:length(node.prob), 1)#, prob = node.prob)
    selected.node <- split.node[[selected.node.idx]]
    rule.parent <- selected.node$split.rule
    
    sel.is.root <- isRoot(selected.node)
    if (sel.is.root) {
      subtree <- tree.new
    } else {
      subtree <- Clone(selected.node)
      subtree$full.obs <- tree.new$full.obs
    }
    
    child.isLeaf <- sapply(subtree$children, isNotLeaf)
    if (all(child.isLeaf) && CART.compare.rule(subtree$children[[1]]$split.rule, subtree$children[[2]]$split.rule)) {
      rule.child <- subtree$children[[1]]$split.rule
      CART.update.rule.swap(subtree$children[[1]], rule.parent)
      CART.update.rule.swap(subtree$children[[2]], rule.parent)
    } else {
      child.node <- subtree$children[[sample(1:2, 1, prob = as.probability(child.isLeaf))]]
      rule.child <- child.node$split.rule
      CART.update.rule.swap(child.node, rule.parent)
    }
    
    if (!CART.update.rule(subtree, rule.child) || !CART.update.obs(subtree)) {
      if (sel.is.root)
        tree.new <- Clone(tree)
      next
    }
    ok <- T
    if (!sel.is.root) {
      subtree$RemoveAttribute('full.obs')
      parent.of.no <- selected.node$parent
      parent.of.no$RemoveChild(selected.node$name)
      parent.of.no$AddChildNode(subtree)
    }
    
    break
  }
  if (!ok)
    return(NULL)
  return(list(
    tree.new=tree.new,
    l.prob=0,
    l.prob.rev=0
  ))
}

