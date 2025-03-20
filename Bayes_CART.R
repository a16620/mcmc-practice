library(data.tree)
library(rpart)
library(dplyr)
library(ggplot2)
library(ggpubr)

##################### Create #####################

BCART <- function(formula, data, method=NULL, print.model=F) {
  combine.df <- model.frame(formula, data)
  y.is.factor <- is.factor(combine.df[,1])
  
  method.list <- c("category", "normal")
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
  
  if (print.model) {
    cat(paste0("Model: ", format(formula)), " method=", method, sep='')
    cat("\nPredictors: ", paste(colnames(combine.df)[-1], collapse=', '), '\n')
  }
  
  model <- list()
  model$formula <- formula
  model$method <- method
  model$tree <- CART.set.obs(Node$new("Root", full.obs=combine.df, pred.dim=1), 1:nrow(combine.df))
  
  return(model)
}

lapBCART <- function(formula, rate=NULL, data, method=NULL, print.model=F) {
  combine.df <- model.frame(formula, data)
  if (is.null(rate) || !(rate %in% colnames(data))) {
    rate <- '.rate.1'
    combine.df <- data.frame(combine.df[, 1, drop=F], .rate.1=rep(1, nrow(combine.df)), combine.df[, 2:ncol(combine.df), drop=F])
  } else {
    ind <- which(rate == colnames(data))[1]
    combine.df <- data.frame(combine.df[, 1, drop=F], combine.df[, ind, drop=F], combine.df[, -c(1,ind), drop=F])
  }
  
  method.list <- c("poisson", "nb", "zip")
  if (is.null(method))
    method <- "poisson"
  if (!(method %in% method.list)) {
    cat("method must be (", paste(method.list, collapse = ', '), ")\n", sep='')
    return(NULL)
  }
  
  if (print.model) {
    cat(paste0("Model: ", format(formula)), " rate=", rate, " method=", method, sep='')
    cat("\nPredictors: ", paste(colnames(combine.df)[-1], collapse=', '), '\n')
  }
  
  model <- list()
  model$formula <- formula
  model$formula <- rate
  model$method <- method
  model$tree <- CART.set.obs(Node$new("Root", full.obs=combine.df, pred.dim=2), 1:nrow(combine.df))
  
  return(model)
}

##################### MCMC #####################

extract.marginal.params <- function(model, uniform=0) {
  method <- model$method
  data <- model$tree$full.obs
  
  if (method == "category") {
    Y <- data[, all.vars(model$formula)[1]]
    if (uniform != 0) {
      a <- rep(uniform, length(levels(Y)))
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
  } else if (method == "poisson") {
    if (uniform != 0) {
      a <- uniform
      b <- uniform
    } else {
      a <- sum(data[,1]) #shape
      b <- sum(data[,2]) #rate
    }
    
    return(list(
      a=a,
      b=b,
      C=a*log(b)-lgamma(a)
    ))
  }
}

MCMC.param <- function(model, split.param=NULL, marginal.param=NULL, marginal.only=T) {
  mcmc.param <- list()
  
  if (is.null(split.param)) {
    mcmc.param$split.prob <- function(level) {
      0.95*level**(-0.5)
    }
  } else {
    if (length(split.param) == 1) {
      mcmc.param$split.prob <- function(level) {
        split.param[1]**level
      }
    } else {
      mcmc.param$split.prob <- function(level) {
        split.param[1]*level**(-split.param[2])
      }
    }
  }
  
  if (is.null(marginal.param))
    marginal.param <- extract.marginal.params(model)
  if (model$method == "category") {
    mcmc.param$tree.marginal <- function(tree) CART.prob.likelihood.marginal.category(tree, marginal.param)
    mcmc.param$memo.writer <- CART.memo.category(marginal.param)
    
    mcmc.param$crit.loss <- "miscl"
    mcmc.param$crit.fn <- CART.get.tree.train.miscl
  } else if (model$method == "normal") {
    mcmc.param$tree.marginal <- function(tree) CART.prob.likelihood.marginal.regression(tree, marginal.param)
    mcmc.param$memo.writer <- CART.memo.regression(marginal.param)
    
    mcmc.param$crit.loss <- "SSE"
    mcmc.param$crit.fn <- CART.get.tree.train.SSE
  } else if (model$method == "poisson") {
    mcmc.param$tree.marginal <- function(tree) CART.prob.likelihood.marginal.poisson(tree, marginal.param)
    mcmc.param$memo.writer <- CART.memo.poisson(marginal.param)
    
    mcmc.param$crit.loss <- "DIC"
    mcmc.param$crit.fn <- CART.get.tree.DIC.poisson
  } else if (model$method == "nb") {
    mcmc.param$tree.marginal <- function(tree) CART.prob.likelihood.marginal.regression(tree, marginal.param)
    mcmc.param$memo.writer <- CART.memo.regression(marginal.param)
  } else if (model$method == "zip") {
    mcmc.param$tree.marginal <- function(tree) CART.prob.likelihood.marginal.regression(tree, marginal.param)
    mcmc.param$memo.writer <- CART.memo.regression(marginal.param)
  }
  
  mcmc.param$use.latent <- !(marginal.only || model$method %in% c("category", "normal", "poisson"))

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
  
  crit.loss <- param$crit.loss
  crit.loss.fn <- param$crit.fn
  
  tree.marginal.lik <- param$tree.marginal
  
  criteria.name <- append(c("log.post", "n.leaf", crit.loss), additional.criteria)
  criteria.funs <- setNames(append(list(
    function(tree) {tree.marginal.lik(tree)+CART.prob.prior(tree, param$split.prob)},
    function(tree) {tree$leafCount},
    crit.loss.fn
  ), additional.criteria.fun), criteria.name)
  
  if (is.null(model.select.criteria))
    model.select.criteria <- crit.loss
  
  split.prob <- param$split.prob
  memo.writer <- param$memo.writer
  
  mtree <- Clone(model0$tree)
  stopifnot(CART.check.tree.ok(mtree))
  mtree$Do(memo.writer)
  
  #burn in
  if (burn.in > 0) {
    if (verbose) {
      cat("burn-in...")
    }
    burn.in.l.post <- tree.marginal.lik(mtree)+CART.prob.prior(mtree, param$split.prob)
    for (iter in 1:burn.in) {
      mc.move <- sample(moves, 1, prob = prob.moves)
      if (mc.move == "grow") {
        mc.new <- CART.move.grow(mtree, split.prob, memo.writer)
      } else if (mc.move == "prune") {
        mc.new <- CART.move.prune(mtree, split.prob, memo.writer)
      } else if (mc.move == "change") {
        mc.new <- CART.move.change(mtree, split.prob, memo.writer, only.value = (runif(1) < 0.5))
      } else if (mc.move == "swap") {
        mc.new <- CART.move.swap(mtree, split.prob, memo.writer)
      }
      
      if (!is.null(mc.new) && CART.check.tree.ok(mc.new$tree.new)) {
        burn.in.l.post.new <- tree.marginal.lik(mc.new$tree.new)+CART.prob.prior(mc.new$tree.new, param$split.prob)
        l.a.ratio.u <- mc.new$l.prob.rev+burn.in.l.post.new
        l.a.ratio.l <- mc.new$l.prob+burn.in.l.post
        
        if (log(runif(1)) <= (l.a.ratio.u-l.a.ratio.l)) {
          mtree <- mc.new$tree.new
          burn.in.l.post <- burn.in.l.post.new
        }
      }
    }
    if (verbose) {
      cat("end!\n")
    }
  }
  
  if (!CART.check.tree.ok(mtree)) {
    stop('burn.in')
  }
  
  #with record
  criteria.matrix <- `colnames<-`(matrix(0, ncol=length(criteria.name), nrow=iteration), criteria.name)
  criteria.matrix[1,] <- sapply(criteria.funs, function(cfn) {cfn(mtree)})
  
  selected.model <- mtree
  selected.model.score <- criteria.matrix[1, model.select.criteria]
  
  print.intval <- if (iteration < 20) iteration else as.integer(iteration*0.1)
  
  accept.count <- setNames(rep(0, length(moves)), moves)
  fail.count <- setNames(rep(0, length(moves)), moves)
  for (iter in 2:iteration) {
    mc.move <- sample(moves, 1, prob = prob.moves)
    if (mc.move == "grow") {
      mc.new <- CART.move.grow(mtree, split.prob, memo.writer)
    } else if (mc.move == "prune") {
      mc.new <- CART.move.prune(mtree, split.prob, memo.writer)
    } else if (mc.move == "change") {
      mc.new <- CART.move.change(mtree, split.prob, memo.writer, only.value = (runif(1) < 0.5))
    } else if (mc.move == "swap") {
      mc.new <- CART.move.swap(mtree, split.prob, memo.writer)
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
        
        if (selected.model.score >= criteria.new.tree[model.select.criteria]) {
          selected.model <- mtree
          selected.model.score <- criteria.new.tree[model.select.criteria]
        }
      } else {
        criteria.new.tree <- criteria.matrix[iter-1,]
      }
    }
    criteria.matrix[iter,] <- criteria.new.tree
    
    if (verbose && iter %% print.intval == 0) {
      cat('============= #', as.character(iter), '(', format(sum(accept.count)/iter, digits = 3),')\n', sep='')
      iter.summary <- colMeans(criteria.matrix[(iter-print.intval+1):(iter),])
      print(iter.summary, quot=F)
    }
  }
  
  criteria.df <- as.data.frame(criteria.matrix)
  if (plot.trace) {
    print(ggarrange(ggplot(criteria.df) + geom_line(aes(x=1:iteration, y=log.post))+xlab('iter')+theme_classic(),
              ggplot() + geom_line(aes(x=1:iteration, y=criteria.df[,crit.loss]))+xlab('iter')+ylab(crit.loss)+theme_classic(),
              ggplot(criteria.df) + geom_line(aes(x=1:iteration, y=n.leaf))+xlab('iter')+theme_classic(),
              ncol = 1))
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
  obs <- CART.get.obs.deps(node)
  rule.colname <- node$split.rule$rule.col
  if (is.factor(obs[,rule.colname])) {
    len.rule.values <- length(unique(obs[,rule.colname]))
    return(CART.prob.select.rule(obs, len.rule.values, T))
  } else {
    return(CART.prob.select.rule(obs, nrow(obs)-1, F))
  }
}

##################### Probability marginal #####################

CART.prob.likelihood.marginal.regression <- function(tree, params) {
  leaf.val.sum <- rowSums(vapply(tree$leaves, function(node) node$memo, numeric(2))) #sum of (log(n.i), t.i+s.i)
  l.prob <- ifelse(is.null(params$c), 0, params$c)+
    0.5*tree$leafCount*log(params$a)-0.5*leaf.val.sum[1]-
    0.5*(nrow(tree$full.obs)+params$nu)*log(leaf.val.sum[2]+params$nu*params$lambda)
  return(l.prob)
}

CART.prob.likelihood.marginal.category <- function(tree, params) {
  l.prob <- sum(vapply(tree$leaves, function(node) node$memo, numeric(1)))
  return(l.prob)
}

CART.prob.likelihood.marginal.poisson <- function(tree, params) {
  l.prob <- sum(vapply(tree$leaves, function(node) {node$memo[1]}, numeric(1)))
  return(l.prob)
}

##################### Criteria #####################

CART.get.tree.train.SSE <- function(tree) {
  sum(vapply(tree$leaves, CART.get.node.train.SSE, numeric(1)))
}

CART.get.node.train.SSE <- function(node) {
  obs.idx <- node$obs.idx
  if (length(obs.idx) == 0) {
    return(0)
  } else {
    y <- CART.get.obs.pred(node)[,1]
    if (length(y) > 1) {
      return(var(y)*(length(y)-1))
    }
    return(0)
  }
}

CART.get.tree.train.miscl <- function(tree) {
  n.total <- nrow(tree$full.obs)
  sum(vapply(tree$leaves, function(node) {
    CART.get.node.train.miscl(node)*length(node$obs.idx)/n.total
  }, numeric(1)))
}

CART.get.node.train.miscl <- function(node) {
  obs.idx <- node$obs.idx
  if (length(obs.idx) == 0) {
    return(0)
  } else {
    y <- CART.get.obs.pred(node)[,1]
    p <- as.numeric(table(y))
    p <- p/sum(p)
    return(1-max(p))
  }
}

CART.get.tree.DIC.poisson <- function(tree) {
  sum(vapply(tree$leaves, function(node) node$memo[2], numeric(1)))
}

##################### Memo #####################

CART.memo.regression <- function(params) {
  c.mu <- params$mu
  c.a <- params$a
  
  function(node) {
    node.obs <- CART.get.obs.pred(node)
    #log(n.i), s.i+t.i
    n.i <- length(node$obs.idx)
    node$memo <- c(log(n.i), (mean(node.obs[,1])-c.mu)**2*c.a*n.i/(c.a+n.i)+CART.get.node.train.SSE(node))
  }
}

CART.memo.category <- function(params) {
  c.a <- params$a
  c.as <- params$a.s
  c.ap <- params$a.p
  
  function(node) {
    node.obs <- CART.get.obs.pred(node)
    n.i <- as.numeric(table(node.obs[,1]))
    node$memo <- c.ap+sum(lgamma(n.i+c.a))-lgamma(sum(n.i)+c.as)
  }
}

CART.memo.poisson <- function(params) {
  c.a <- params$a
  c.b <- params$b
  c.C <- params$C #a*log(b)-lgamma(a)
  
  function(node) {
    node.obs <- CART.get.obs.pred(node)
    N.sum <- sum(node.obs[,1])
    V.sum <- sum(node.obs[,2])
    NV.sum <- sum(node.obs[,1]*log(node.obs[,2]))
    lg.sum <- sum(lgamma(node.obs[,1]+1))
    N.sum.a <- N.sum+c.a
    l.post <- c.C+NV.sum-lg.sum+
      lgamma(N.sum.a)-(N.sum.a)*log(V.sum+c.b)
    nv.ratio <- N.sum.a/(V.sum+c.b)
    dic <- 2*(nv.ratio*V.sum-log(nv.ratio)*N.sum-NV.sum+lg.sum+2*(log(N.sum.a)-digamma(N.sum.a))*N.sum)
    node$memo <- c(l.post, dic, nv.ratio)
  }
}

##################### Utility #####################

CART.check.tree.ok <- function(tree) {
  all(vapply(tree$leaves, function(node) {length(node$obs.idx) > 0}, logical(1)))
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
    num.to.str <- format(split.value, digits=4)
    return(c(
      paste0(split.col, "<=", num.to.str),
      paste0(split.col, ">", num.to.str)
    ))
  }
}

CART.get.obs <- function(node) {
  return(node$root$full.obs[node$obs.idx, , drop=F])
}

CART.get.obs.pred <- function(node) {
  root <- node$root
  return(root$full.obs[node$obs.idx, 1:(root$pred.dim), drop=F])
}

CART.get.obs.deps <- function(node) {
  root <- node$root
  return(root$full.obs[node$obs.idx, -(1:(root$pred.dim)), drop=F])
}

CART.set.obs <- function(node, obs.idx) {
  node$obs.idx <- obs.idx
  return(node)
}

CART.clean.obs <- function(tree) {
  tree$RemoveAttribute('full.obs')
  subtree$RemoveAttribute('pred.dim')
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
  obs.deps <- CART.get.obs.deps(node)
  col.candidates <- if (is.null(col.candidates)) colnames(obs.deps) else col.candidates
  while (length(col.candidates) > 0) {
    rule.colname <- sample(col.candidates, 1)
    if (is.factor(obs.deps[,rule.colname])) {
      rule.values <- unique(obs.deps[,rule.colname])
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
      rule.values <- obs.deps[,rule.colname]
      filt.values <- rule.values[rule.values < max(rule.values)]
      len.filt.values <- length(filt.values)
      
      
      if (len.filt.values <= 2) {
        col.candidates <- col.candidates[col.candidates != rule.colname]
        next
      }
      
      val.idx <- sample(1:len.filt.values, 1)
      rule.value <- runif(1, min=filt.values[val.idx], max=min(rule.values[rule.values > filt.values[val.idx]]))
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
  
  rule.info$fun <- rule
  node$split.rule <- rule.info
  
  child.names <- CART.rule2name(CART.get.obs.deps(node), rule.info$split.col, rule.info$split.value)
  node$AddChild(child.names[1])
  node$AddChild(child.names[2])
  
  return(CART.update.obs(node))
}

##################### Change, Swap #####################

CART.update.obs <- function(tree) {
  tree$Do(function(node) {
    obs.key <- node$split.rule$fun(CART.get.obs.deps(node))
    
    left.obs.idx <- node$obs.idx[obs.key]
    right.obs.idx <- node$obs.idx[!obs.key]
    
    CART.set.obs(node$children[[1]], left.obs.idx)
    CART.set.obs(node$children[[2]], right.obs.idx)
  }, traversal="level", filterFun = isNotLeaf)
  return(CART.check.tree.ok(tree))
}

CART.update.rule.swap <- function(node, rule.info) {
  node$split.rule <- rule.info
  
  child.names <- CART.rule2name(CART.get.obs.deps(node), rule.info$split.col, rule.info$split.value)
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
  
  rule.info$fun <- rule
  node$split.rule <- rule.info
  
  child.names <- CART.rule2name(CART.get.obs.deps(node), rule.info$split.col, rule.info$split.value)
  names(node$children) <- child.names
  node$children[[1]]$name <- child.names[1]
  node$children[[2]]$name <- child.names[2]
  
  return(CART.update.obs(node))
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

CART.move.grow <- function(tree, split.prob, memo.writer) {
  tree.new <- Clone(tree)
  
  terminal.nodes <- tree.new$leaves
  node.level <- sapply(terminal.nodes, function(node) node$level)
  node.prob <- as.probability(split.prob(node.level))
  
  while (length(terminal.nodes) > 0) {
    selected.node.idx <- sample(1:length(node.prob), 1, prob = node.prob)
    selected.node <- terminal.nodes[[selected.node.idx]]
    rule <- CART.select.rule(selected.node)
    if (!is.null(rule) && CART.set.rule(selected.node, rule)) {
      selected.node$RemoveAttribute('memo')
      Do(selected.node$children, memo.writer)
      break
    }
    terminal.nodes <- terminal.nodes[-selected.node.idx]
    node.prob <- as.probability(node.prob[-selected.node.idx])
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

CART.move.prune <- function(tree, split.prob, memo.writer) {
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
  memo.writer(selected.node)
  
  return(list(
    tree.new=tree.new,
    l.prob=log(node.prob[selected.node.idx]),
    l.prob.rev=log(node.prob.rev[selected.node.idx])+l.prob.sel.split
  ))
}

CART.move.change <- function(tree, split.prob, memo.writer, only.value=F, max.try=100) {
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
  for (try in 1:max.try) {
    selected.node.idx <- sample(1:length(node.prob), 1)#, prob = node.prob)
    selected.node <- split.node[[selected.node.idx]]
    if (only.value)
      rule <- CART.select.rule(selected.node, selected.node$split.rule$rule.col)
    else
      rule <- CART.select.rule(selected.node)
    if (is.null(rule))
      next
    
    sel.isn.root <- isNotRoot(selected.node)
    subtree <- Clone(selected.node)
    if (sel.isn.root) {
      subtree$full.obs <- tree.new$full.obs
      subtree$pred.dim <- tree.new$pred.dim
    }
    
    if (!CART.update.rule(subtree, rule)) {
      next
    }
    ok <- T
    if (sel.isn.root) {
      subtree$RemoveAttribute('full.obs')
      subtree$RemoveAttribute('pred.dim')
      parent.of.no <- selected.node$parent
      parent.of.no$RemoveChild(selected.node$name)
      parent.of.no$AddChildNode(subtree)
    } else {
      tree.new <- subtree
    }
    subtree$Do(memo.writer, filterFun = isLeaf)
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

CART.move.swap <- function(tree, split.prob, memo.writer, max.try=100) {
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
  for (try in 1:max.try) {
    selected.node.idx <- sample(1:length(node.prob), 1)#, prob = node.prob)
    selected.node <- split.node[[selected.node.idx]]
    rule.parent <- selected.node$split.rule
    
    sel.isn.root <- isNotRoot(selected.node)
    subtree <- Clone(selected.node)
    if (sel.isn.root) {
      subtree$full.obs <- tree.new$full.obs
      subtree$pred.dim <- tree.new$pred.dim
    }
    
    child.isNotLeaf <- sapply(subtree$children, isNotLeaf)
    if (all(child.isNotLeaf) && CART.compare.rule(subtree$children[[1]]$split.rule, subtree$children[[2]]$split.rule)) {
      rule.child <- subtree$children[[1]]$split.rule
      CART.update.rule.swap(subtree$children[[1]], rule.parent)
      CART.update.rule.swap(subtree$children[[2]], rule.parent)
    } else {
      child.node <- subtree$children[[sample(1:2, 1, prob = as.probability(child.isNotLeaf))]]
      rule.child <- child.node$split.rule
      CART.update.rule.swap(child.node, rule.parent)
    }
    
    if (!CART.update.rule(subtree, rule.child)) {
      next
    }
    ok <- T
    if (sel.isn.root) {
      subtree$RemoveAttribute('full.obs')
      subtree$RemoveAttribute('pred.dim')
      parent.of.no <- selected.node$parent
      parent.of.no$RemoveChild(selected.node$name)
      parent.of.no$AddChildNode(subtree)
    } else {
      tree.new <- subtree
    }
    subtree$Do(memo.writer, filterFun = isLeaf)
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

