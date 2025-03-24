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
    combine.df <- data.frame(combine.df[, 1, drop=F], rep(1, nrow(combine.df)), combine.df[, 2:ncol(combine.df), drop=F])
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
  
  formula.txt <- paste0(rate, (if (is.null(rate)) NULL else '*'), format(formula))
  if (print.model) {
    cat(paste0("Model: ", formula.txt, " method=", method, sep=''))
    cat("\nPredictors: ", paste(colnames(combine.df)[-1], collapse=', '), '\n')
  }
  
  model <- list()
  model$formula <- formula
  model$rate <- rate
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
  } else if (method %in% c("poisson", "nb")) {
    if (uniform != 0) {
      a <- uniform
      b <- uniform
    } else {
      a <- mean(data[,1]) #shape
      b <- mean(data[,2]) #rate
    }
    
    return(list(
      a=a,
      b=b
    ))
  } else if (method == "zip") {
    if (uniform != 0) {
      a <- uniform
      b <- uniform
    } else {
      a <- mean(data[,1]) #shape
      b <- mean(data[,2]) #rate
    }
    
    return(list(
      a1=1,
      b1=1,
      a2=a,
      b2=b
    ))
  }
}

MCMC.param <- function(model, split.param=NULL, marginal.param=NULL, marginal.only=T) {
  mcmc.param <- list()
  
  if (is.null(split.param)) {
    split.param <- c(0.95, 0.5)
  }
  
  if (length(split.param) == 1) {
    log.gamma <- log(split.param)
    if (log.gamma > 0)
      stop("split.param r: 0<r<1")
    mcmc.param$split.prob.prior <- function(level) {
      l.p <- level*log.gamma
      l.p.inv <- ifelse(-0.69 > l.p, -l.p, log(-l.p))
      if (any(is.na(c(l.p, l.p.inv)))) {
        print(cbind(l.p, l.p.inv, level))
      }
      cbind(l.p, l.p.inv)
    }
    mcmc.param$split.prob.select <- function(level) {
      split.param**level
    }
  } else {
    log.alpha <- log(split.param[1])
    m.beta <- -split.param[2]
    mcmc.param$split.prob.prior <- function(level) {
      l.p <- log.alpha+m.beta*log(level)
      l.p.inv <- ifelse(-0.69 > l.p, -l.p, log(-l.p))
      cbind(l.p, l.p.inv)
    }
    mcmc.param$split.prob.select <- function(level) {
      level**m.beta
    }
  }
  
  if (is.null(marginal.param))
    marginal.param <- extract.marginal.params(model)
  if (model$method == "category") {
    mcmc.param$tree.marginal <- function(tree) CART.prob.likelihood.marginal.category(tree, marginal.param)
    mcmc.param$after.move <- CART.memo.category(marginal.param)
    mcmc.param$crit.loss <- "miscl"
    mcmc.param$crit.fn <- CART.criteria.train.miscl
    
  } else if (model$method == "normal") {
    mcmc.param$tree.marginal <- function(tree) CART.prob.likelihood.marginal.regression(tree, marginal.param)
    mcmc.param$after.move <- CART.memo.regression(marginal.param)
    mcmc.param$crit.loss <- "SSE"
    mcmc.param$crit.fn <- CART.criteria.train.SSE
    
  } else if (model$method == "poisson") {
    mcmc.param$tree.marginal <- function(tree) CART.prob.likelihood.marginal.memo(tree, marginal.param)
    mcmc.param$after.move <- CART.memo.poisson(marginal.param)
    mcmc.param$crit.loss <- "DIC"
    mcmc.param$crit.fn <- CART.criteria.augment.DIC
    
  } else if (model$method == "nb") {
    mcmc.param$tree.marginal <- function(tree) CART.prob.likelihood.marginal.memo(tree, marginal.param)
    mcmc.param$after.move <- CART.memo.NB(marginal.param)
    mcmc.param$crit.loss <- "DIC"
    mcmc.param$crit.fn <- CART.criteria.augment.DIC
    
  } else if (model$method == "zip") {
    mcmc.param$tree.marginal <- function(tree) CART.prob.likelihood.marginal.memo(tree, marginal.param)
    mcmc.param$after.move <- CART.memo.ZIP(marginal.param)
    mcmc.param$crit.loss <- "DIC"
    mcmc.param$crit.fn <- CART.criteria.augment.DIC
    
  }
  
  mcmc.param$use.augment <- model$method %in% c("nb", "zip")
  if (mcmc.param$use.augment) {
    if (model$method == "nb") {
      mcmc.param$augment.sampler <- CART.augment.sampler.NB
      mcmc.param$augment.init <- CART.augment.initialzer.NB(marginal.param)
      
    } else if (model$method == "zip") {
      mcmc.param$augment.sampler <- CART.augment.sampler.ZIP
      mcmc.param$augment.init <- CART.augment.initialzer.ZIP(marginal.param)
      
    }
  }
  
  return(mcmc.param)
}

do.MCMC <- function(model0, param=NULL, iteration=1000, burn.in=0,
                    moves=c("grow", "prune", "change", "swap"), prob.moves=NULL, adaptive.move=0,
                    verbose=F, plot.trace=T, additional.criteria=NULL, additional.criteria.fun=NULL, model.select.criteria=NULL) {
  stopifnot(length(moves) > 0)
  stopifnot(length(additional.criteria) == length(additional.criteria.fun))

  if (is.null(param)) {
    param <- MCMC.param(model0)
  }
  
  crit.loss <- param$crit.loss
  crit.loss.fn <- param$crit.fn
  
  tree.marginal.lik <- param$tree.marginal
  tree.prior.prob <- param$split.prob.prior
  
  default.crit.name <- c("log.post", "n.leaf", crit.loss)
  default.crit.fn <- list(
    function(tree) {tree.marginal.lik(tree)+CART.prob.prior(tree, tree.prior.prob)},
    function(tree) {tree$leafCount},
    crit.loss.fn
  )
  use.augment <- param$use.augment
  if (use.augment) {
    default.crit.name <- append(default.crit.name, "aug.llik")
    default.crit.fn <- append(default.crit.fn, CART.criteria.augment.llik)
  }
  
  criteria.name <- append(default.crit.name, additional.criteria)
  criteria.funs <- setNames(append(default.crit.fn, additional.criteria.fun), criteria.name)
  
  if (is.null(model.select.criteria))
    model.select.criteria <- crit.loss
  
  split.prob <- param$split.prob.select
  after.move <- param$after.move
  
  mtree <- Clone(model0$tree)
  stopifnot(CART.check.tree.ok(mtree))
  
  if (use.augment) {
    augment.sampler <- param$augment.sampler
    param$augment.init(mtree)
  } else {
    augment.sampler <- NULL
  }
  Do(mtree$leaves, after.move)
  
  if (adaptive.move > 0) {
    if (is.null(prob.moves))
      prob.moves <- as.probability(rep(1, length(moves)))
    prob.moves0 <- prob.moves
  }
  
  #burn in
  if (burn.in > 0) {
    if (verbose) {
      cat("burn-in...")
    }
    move.count <- setNames(rep(1, length(moves)), moves)
    accept.count <- setNames(rep(0, length(moves)), moves)
    burn.in.l.post <- tree.marginal.lik(mtree)+CART.prob.prior(mtree, tree.prior.prob)
    for (iter in 1:burn.in) {
      mc.move <- sample(moves, 1, prob = prob.moves)
      move.count[mc.move] <- move.count[mc.move]+1
      if (mc.move == "grow") {
        mc.new <- CART.move.grow(mtree, split.prob, after.move, augment.sampler)
      } else if (mc.move == "prune") {
        mc.new <- CART.move.prune(mtree, split.prob, after.move, augment.sampler)
      } else if (mc.move == "change") {
        mc.new <- CART.move.change(mtree, split.prob, after.move, augment.sampler, only.value = (runif(1) < 0.5))
      } else if (mc.move == "swap") {
        mc.new <- CART.move.swap(mtree, split.prob, after.move, augment.sampler)
      }
      
      if (!is.null(mc.new) && CART.check.tree.ok(mc.new$tree.new)) {
        burn.in.l.post.new <- tree.marginal.lik(mc.new$tree.new)+CART.prob.prior(mc.new$tree.new, tree.prior.prob)
        l.a.ratio.u <- mc.new$l.prob.rev+burn.in.l.post.new
        l.a.ratio.l <- mc.new$l.prob+burn.in.l.post

        if (log(runif(1)) <= (l.a.ratio.u-l.a.ratio.l)) {
          accept.count[mc.move] <- accept.count[mc.move] + 1
          mtree <- mc.new$tree.new
          burn.in.l.post <- burn.in.l.post.new
        } else if (use.augment) {
          CART.set.augment.obs(mtree, CART.get.augment.obs(mc.new$tree.new))
          burn.in.l.post <- tree.marginal.lik(mtree)+CART.prob.prior(mtree, tree.prior.prob)
        }
      }
      
      if (adaptive.move > 0) {
        prob.moves <- as.probability((prob.moves0+adaptive.move*accept.count/move.count)/(1+adaptive.move))
        #(prob.moves0+accept.count)/(iter+1), (prob.moves0+accept.count/iter)/2
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

  move.count <- setNames(rep(1, length(moves)), moves)
  accept.count <- setNames(rep(0, length(moves)), moves)
  fail.count <- setNames(rep(0, length(moves)), moves)
  for (iter in 2:iteration) {
    mc.move <- sample(moves, 1, prob = prob.moves)
    move.count[mc.move] <- move.count[mc.move]+1
    if (mc.move == "grow") {
      mc.new <- CART.move.grow(mtree, split.prob, after.move, augment.sampler)
    } else if (mc.move == "prune") {
      mc.new <- CART.move.prune(mtree, split.prob, after.move, augment.sampler)
    } else if (mc.move == "change") {
      mc.new <- CART.move.change(mtree, split.prob, after.move, augment.sampler, only.value = (runif(1) < 0.5))
    } else if (mc.move == "swap") {
      mc.new <- CART.move.swap(mtree, split.prob, after.move, augment.sampler)
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
      } else if (use.augment) {
        CART.set.augment.obs(mtree, CART.get.augment.obs(mc.new$tree.new))
        criteria.new.tree <- sapply(criteria.funs, function(cfn) {cfn(mtree)})
      }else {
        criteria.new.tree <- criteria.matrix[iter-1,]
      }
      
      if (selected.model.score >= criteria.new.tree[model.select.criteria]) {
        selected.model <- mtree
        selected.model.score <- criteria.new.tree[model.select.criteria]
      }
    }
    
    if (adaptive.move > 0) {
      prob.moves <- as.probability((prob.moves0+adaptive.move*accept.count/move.count)/(1+adaptive.move))
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
    plot.list <- list(ggplot(criteria.df) + geom_line(aes(x=1:iteration, y=log.post))+xlab('iter')+theme_classic(),
                  ggplot() + geom_line(aes(x=1:iteration, y=criteria.df[,crit.loss]))+xlab('iter')+ylab(crit.loss)+theme_classic(),
                  ggplot(criteria.df) + geom_line(aes(x=1:iteration, y=n.leaf))+xlab('iter')+theme_classic())
    if (use.augment)
      plot.list <- append(plot.list, list(ggplot(criteria.df) + geom_line(aes(x=1:iteration, y=aug.llik))+xlab('iter')+theme_classic()))
    plot.title <- paste0(model0$method, ':', model0$rate, (if (is.null(model0$rate)) NULL else '*'), format(model0$formula))
    print(annotate_figure(ggarrange(plotlist = plot.list, ncol = 1), top=text_grob(plot.title, face='bold')))
  }
  model0$tree <- selected.model
  return(list(
    history=criteria.df,
    move.error=fail.count,
    move.try=move.count-1,
    move.accept=accept.count,
    
    accept.rate=sum(accept.count)/iteration,
    model=model0,
    score=selected.model.score
  ))
}

##################### Deploy #####################

deploy <- function(model) {
  tree <- model$tree$root
  if (!is.null(tree$deploy)) {
    warning("Already deployed.")
    return()
  }
  if (model$method == "normal") {
    tree$deploy <- list(pred.dim=2, attr.name=c("mean", "var"))
    Do(tree$leaves, function (node) {
      node.obs <- CART.get.obs.pred(node)[,1]
      node$deploy.pred <- c(mean(node.obs), var(node.obs))
    })
  } else if (model$method == "category") {
    tree$deploy <- list(pred.dim=2, attr.name=c("class", "prob"))
    Do(tree$leaves, function (node) {
      node.obs <- CART.get.obs.pred(node)[,1]
      pred <- which.max(table(node.obs))
      node$deploy.pred <- c(names(pred), unname(pred))
    })
  } else if (model$method == "poisson") {
    tree$deploy <- list(pred.dim=1, attr.name=c("lambda"))
    Do(tree$leaves, function (node) {
      node$deploy.pred <- node$memo[3]
    })
  } else if (model$method == "nb") {
    tree$deploy <- list(pred.dim=2, attr.name=c("lambda", "K"))
    Do(tree$leaves, function (node) {
      node$deploy.pred <- node$memo[5:6]
    })
  } else if (model$method == "zip") {
    tree$deploy <- list(pred.dim=2, attr.name=c("lambda", "mu"))
    Do(tree$leaves, function (node) {
      node$deploy.pred <- node$memo[6:7]
    })
  }
  
  tree$RemoveAttribute('full.obs')
  tree$RemoveAttribute('pred.dim')
  if (!is.null(tree$augments))
    tree$RemoveAttribute('augments')
  
  tree$Do(function(node) {
    if (!is.null(node$memo))
      node$RemoveAttribute('memo')
    node$RemoveAttribute('obs.idx')
  })
}

predict <- function(model, obs) {
  tree <- model$tree
  tree$tmp.forward.idx <- 1:nrow(obs)
  pred.mat <- `colnames<-`(matrix(0, ncol=tree$deploy$pred.dim, nrow=nrow(obs)), tree$deploy$attr.name)
  
  tree$Do(function(node) {
    obs.key <- node$split.rule$fun(obs[node$tmp.forward.idx,])
    
    left.obs.idx <- node$tmp.forward.idx[obs.key]
    right.obs.idx <- node$tmp.forward.idx[!obs.key]
    
    node$children[[1]]$tmp.forward.idx <- left.obs.idx
    node$children[[2]]$tmp.forward.idx <- right.obs.idx
    
    node$RemoveAttribute('tmp.forward.idx')
  }, traversal="level", filterFun=isNotLeaf)
  
  for (node in tree$leaves) {
    pred.mat[node$tmp.forward.idx,] <- node$deploy.pred
    node$RemoveAttribute('tmp.forward.idx')
  }
  return(pred.mat)
}

##################### Probability #####################

CART.prob.prior <- function(tree, split.prob) {
  node.info <- tree$Get(function(node) c(node$level, isLeaf(node)), simplify = T)
  node.prob <- split.prob(node.info[1,])
  sum(ifelse(node.info[2,], node.prob[,1], node.prob[,2]))
}

CART.prob.rule <- function(node) {
  rule <- node$split.rule
  rule.val <- CART.get.obs.deps(node)
  len.rule.values <- nrow(rule.val)
  if (len.rule.values <= 1)
    return(0)
  
  lp.select.col <- -log(ncol(rule.val))
  rule.col.val <- rule.val[,rule$split.col]
  if (is.factor(rule.col.val)) {
    lp.select.val <- -log(max(2**length(levels(rule.col.val))-2, 1)) #unique or levels
  } else {
    lp.select.val <- 0.6931472-log(len.rule.values)-log(len.rule.values-1)+rule$l.aux
  }
  return(lp.select.val+lp.select.col)
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

CART.prob.likelihood.marginal.memo <- function(tree, params) {
  l.prob <- sum(vapply(tree$leaves, function(node) {node$memo[1]}, numeric(1)))
  return(l.prob)
}

##################### Criteria #####################

CART.criteria.train.SSE <- function(tree) {
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

CART.criteria.train.miscl <- function(tree) {
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

CART.criteria.augment.DIC <- function(tree) {
  sum(vapply(tree$leaves, function(node) node$memo[2], numeric(1)))
}

CART.criteria.augment.llik <- function(tree) {
  sum(vapply(tree$leaves, function(node) node$memo[3], numeric(1)))
}

##################### Augment #####################

CART.augment.param.sampling <- function(tree, sampler) {
  Do(tree$leaves, sampler)
}

CART.augment.sampler.NB <- function(node) {
  node.obs <- CART.get.obs.pred(node)
  N.t <- node.obs[,1]
  v.t <- node.obs[,2]
  K.t <- node$memo[6]
  lambda.t <- node$memo[4]
  
  xi.t <- rgamma(nrow(node.obs), K.t*v.t+N.t, K.t*v.t+lambda.t*v.t)
  CART.set.augment.obs(node, xi.t)
}

CART.augment.initialzer.NB <- function(params) {
  function (tree) {
    CART.init.augment(tree, 1)
    Do(tree$leaves, function(node) {
      node.obs <- CART.get.obs.pred(node)
      N.t <- node.obs[,1]
      v.t <- node.obs[,2]
      N.sum <- sum(N.t)
      v.sum <- sum(v.t)
      n.t <- nrow(node.obs)
      
      lambda.est.t <- N.sum/v.sum
      V.est.t2 <- 1/(n.t-1)*sum(v.t*(N.t/v.t-lambda.est.t)**2)
      K.t <- max(lambda.est.t**2/(V.est.t2-lambda.est.t)/(n.t-1)*(v.sum-sum(v.t**2)/v.sum), 0)
      
      if (K.t < 0)
        K.t <- -1/K.t
      else if (K.t == 0)
        K.t <- 1e-5
      
      lambda.t <- rgamma(1, params$a, params$b)
      
      xi.t <- rgamma(nrow(node.obs), K.t*v.t+N.t, K.t*v.t+lambda.t*v.t)
      CART.set.augment.obs(node, xi.t)
    })
  }
}

CART.augment.sampler.ZIP <- function(node) {
  node.obs <- CART.get.obs.pred(node)
  N.t <- node.obs[,1]
  v.t <- node.obs[,2]

  mu.t <- node$memo[4]
  lambda.t <- node$memo[5]
  n.t <- nrow(node.obs)
  
  odd.t <- mu.t*exp(-lambda.t*v.t)
  delta.t <- ifelse(N.t == 0, rbinom(n.t, 1, odd.t/(1+odd.t)), 1)
  phi.t <- rexp(n.t, mu.t*v.t+1)
  CART.set.augment.obs(node, cbind(delta.t, phi.t))
}

CART.augment.initialzer.ZIP <- function(params) {
  function (tree) {
    CART.init.augment(tree, 2)
    Do(tree$leaves, function(node) {
      node.obs <- CART.get.obs.pred(node)
      N.t <- node.obs[,1]
      v.t <- node.obs[,2]
      n.t <- nrow(node.obs)
      
      mu.t <- rgamma(1, params$a1, params$b1)
      lambda.t <- rgamma(1, params$a2, params$b2)
      
      odd.t <- mu.t*exp(-lambda.t*v.t)
      delta.t <- ifelse(N.t == 0, rbinom(n.t, 1, odd.t/(1+odd.t)), 1)
      phi.t <- rexp(n.t, mu.t*v.t+1)
      CART.set.augment.obs(node, cbind(delta.t, phi.t))
    })
  }
}

##################### Memo #####################

#[2] M1 M2
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

#[1] l.prob
CART.memo.category <- function(params) {
  c.a <- params$a
  c.as <- sum(c.a)
  c.ap <- lgamma(c.as)-sum(lgamma(c.a))
  
  function(node) {
    node.obs <- CART.get.obs.pred(node)
    n.i <- as.numeric(table(node.obs[,1]))
    node$memo <- c.ap+sum(lgamma(n.i+c.a))-lgamma(sum(n.i)+c.as)
  }
}

#[3] l.prob DIC post.lambda
CART.memo.poisson <- function(params) {
  c.a <- params$a
  c.b <- params$b
  c.C <- c.a*log(c.b)-lgamma(c.a)
  
  function(node) {
    node.obs <- CART.get.obs.pred(node)
    N.sum <- sum(node.obs[,1])
    V.sum <- sum(node.obs[,2])
    NV.sum <- sum(node.obs[,1]*log(node.obs[,2]))
    lg.sum <- sum(lgamma(node.obs[,1]+1))
    N.sum.a <- N.sum+c.a
    l.prob <- c.C+NV.sum-lg.sum+
      lgamma(N.sum.a)-(N.sum.a)*log(V.sum+c.b)
    nv.ratio <- N.sum.a/(V.sum+c.b)
    dic <- 2*(nv.ratio*V.sum-log(nv.ratio)*N.sum-NV.sum+lg.sum+2*(log(N.sum.a)-digamma(N.sum.a))*N.sum)
    node$memo <- c(l.prob, dic, nv.ratio)
  }
}

#[6] l.prob DIC data.l.prob lambda.t | post.lambda K.t
CART.memo.NB <- function(params) {
  c.a <- params$a
  c.b <- params$b
  c.C <- c.a*log(c.b)-lgamma(c.a)
  
  function(node) {
    node.obs <- CART.get.obs.pred(node)
    N.t <- node.obs[,1]
    v.t <- node.obs[,2]
    sum.N <- sum(N.t)
    sum.V <- sum(v.t)
    n.t <- nrow(node.obs)
    
    if (sum.N > 0) {
      lambda.est.t <- sum.N/sum.V
      V.est.t2 <- sum(v.t*(N.t/v.t-lambda.est.t)**2)/(n.t-1)
      K.t <- lambda.est.t**2/(V.est.t2-lambda.est.t)/(n.t-1)*(sum.V-sum(v.t**2)/sum.V)
      if (K.t < 0)
        K.t <- -1/K.t
      else if (K.t == 0)
        K.t <- 1e-5
    } else {
      K.t <- 1
    }
    
    xi.t <- CART.get.augment.obs(node)
  
    #memo
    log.v <- log(v.t)
    sum.v.log.v <- sum(v.t*log.v)
    sum.lg.N <- sum(lgamma(N.t+1))
    sum.lg.kv <- sum(lgamma(K.t*v.t))
    sum.n.l.v <- sum(N.t*log.v)
    M.1 <- K.t*(log(K.t)*sum.V+sum.v.log.v)
    M.3 <- sum((K.t*v.t+N.t-1)*log(xi.t))
    sum.n.a <- sum.N+c.a
    sum.xv.b <- sum(xi.t*v.t)+c.b
    
    #NB2
    l.prob <- c.C+M.1+sum.n.l.v-(sum.lg.N+sum.lg.kv)+M.3-K.t*(sum.xv.b-c.b)+lgamma(sum.n.a)-sum.n.a*log(sum.xv.b)
    nxv.ratio <- sum.n.a/sum.xv.b
    dic <- 2*( -sum.n.l.v+(log(sum.xv.b)-log(sum.n.a))*sum.N+nxv.ratio*(sum.xv.b-c.b)+
                 sum.lg.kv+sum.lg.N-M.1-M.3+K.t*(sum.xv.b-c.b)+
                 1+2*(log(sum.n.a)-digamma(sum.n.a))*sum.N )
    
    lambda.t <- rgamma(1, sum.n.a, sum.xv.b)
    
    l.prob.data <- sum.N*log(lambda.t)+sum.n.l.v-(lambda.t+K.t)*(sum.xv.b-c.b)+
      M.1+M.3-(sum.lg.N+sum.lg.kv)
    
    node$memo <- c(l.prob, dic, l.prob.data, lambda.t, (sum.n.a/sum.xv.b), K.t)
  }
}

#[7] l.prob DIC data.l.prob mu.t, lambda.t | post.lambda post.mu
CART.memo.ZIP <- function(params) {
  c.a1 <- params$a1
  c.b1 <- params$b1
  c.a2 <- params$a2
  c.b2 <- params$b2
  c.C <- c.a1*log(c.b1)-lgamma(c.a1)+c.a2*log(c.b2)-lgamma(c.a2)
  
  function(node) {
    obs <- CART.get.obs.pred(node)
    N.t <- obs[,1]
    v.t <- obs[,2]
    
    aug.t <- CART.get.augment.obs(node)
    delta.t <- aug.t[,1]
    phi.t <- aug.t[,2]
    
    sum.phi.t <- sum(phi.t)
    sum.delta <- sum(delta.t)
    sum.delta.a <- sum.delta+c.a1
    sum.delta.N <- sum(delta.t*N.t)
    sum.v.phi.b <- sum(v.t*phi.t)+c.b1
    sum.delta.lg.N <- sum(delta.t*lgamma(N.t+1))
    
    mu.bar <- sum.delta.a/(sum(phi.t*v.t)+c.b1)
    lambda.bar <- (sum.delta.N+c.a2)/(sum.delta+c.b2)
    
    #ZIP2
    l.prob <- c.C-sum.phi.t-sum(delta.t*log(v.t))-sum.delta.lg.N+
      lgamma(sum.delta.a)-sum.delta.a*log(sum.v.phi.b)+
      lgamma(sum.delta.N+c.a2)-(sum.delta.N+c.a2)*log(sum.delta+c.b2)
    dic <- 2*(
      -sum(log((N.t==0)/(1+mu.bar*v.t)+mu.bar*N.t*lambda.bar**N.t*exp(-lambda.bar)/(1+mu.bar*v.t)/factorial(N.t)))+
        2*((log(sum.delta.a)-digamma(sum.delta.a))*sum.delta+
             (log(sum.delta.N+c.a2)-digamma(sum.delta.N+c.a2))*sum.delta.N
        ))
    
    mu.t <- rgamma(1, sum.delta.a, sum.v.phi.b)
    lambda.t <- rgamma(1, sum.delta.N+c.a2, sum.delta+c.b2)
    
    l.prob.data <- -sum.phi.t-(sum.v.phi.b-c.b1)*mu.t+sum.delta*log(mu.t)+sum(delta.t*log(v.t))+
      sum.delta.N*log(lambda.t)-sum.delta.lg.N
    
    node$memo <- c(l.prob, dic, l.prob.data, mu.t, lambda.t, lambda.bar, mu.bar)
  }
}

##################### Utility #####################

CART.check.tree.ok <- function(tree) {
  all(vapply(tree$leaves, function(node) {length(node$obs.idx) >= 2}, logical(1)))
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

CART.get.augment.obs <- function(node) {
  root <- node$root
  return(root$augments[node$obs.idx, ])
}

CART.set.augment.obs <- function(node, values) {
  root <- node$root
  root$augments[node$obs.idx, ] <- values
}

CART.init.augment <- function(tree, ndim) {
  root <- tree$root
  root$augments <- matrix(0, nrow=nrow(root$full.obs), ncol=ndim)
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

as.probability <- function(unp) {
  unp/sum(unp)
}

##################### Grow #####################

CART.select.rule <- function(node, col.candidates=NULL, safe.variance=T) {
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
      if (safe.variance)
        rule.values <- rule.values[(min(rule.values) < rule.values) & (rule.values < max(rule.values))]
      len.rule.values <- length(rule.values)
      
      if (len.rule.values < 2 || min(rule.values) == max(rule.values)) {
        col.candidates <- col.candidates[col.candidates != rule.colname]
        next
      }
      
      smp <- sample(rule.values, 2)
      l.bound <- min(smp)
      u.bound <- max(smp)
      
      aux.param <- mean(rule.values >= mean(smp))
      u.aux <- rbeta(1, aux.param, 1)
      
      l.aux <- if (l.bound == u.bound) 0 else (dbeta(u.aux, aux.param, 1, log=T)-log(u.bound-l.bound))
      rule.value <- l.bound+(u.bound-l.bound)*u.aux
      #rule.value <- runif(1, l.bound, u.bound)
      #l.aux <- dunif(rule.value, l.bound, u.bound, log=T)
      return(
        list(
          split.col=rule.colname,
          split.value=rule.value,
          l.aux=l.aux
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

CART.move.grow <- function(tree, split.prob, after.move, augment.sampler) {
  tree.new <- Clone(tree)
  if (!is.null(augment.sampler))
    Do(tree.new$leaves, augment.sampler)
  
  terminal.nodes <- tree.new$leaves
  node.level <- unname(vapply(terminal.nodes, function(node) node$level, numeric(1)))
  node.prob <- as.probability(split.prob(node.level))
  
  while (length(terminal.nodes) > 0) {
    selected.node.idx <- sample(1:length(node.prob), 1, prob = node.prob)
    selected.node <- terminal.nodes[[selected.node.idx]]
    if (length(selected.node$obs.idx) > 5) {
      rule <- CART.select.rule(selected.node)
      if (!is.null(rule) && CART.set.rule(selected.node, rule)) {
        selected.node$RemoveAttribute('memo')
        Do(selected.node$children, after.move)
        break
      }
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

CART.move.prune <- function(tree, split.prob, after.move, augment.sampler) {
  if (tree$leafCount == 1) {
    return(NULL)
  }
  tree.new <- Clone(tree)
  if (!is.null(augment.sampler))
    Do(tree.new$leaves, augment.sampler)
  
  parent.of <- Traverse(tree.new, filterFun = function(node) {
    isNotLeaf(node) && isLeaf(node$children[[1]]) && isLeaf(node$children[[2]])
  })
  node.level <- unname(sapply(parent.of, function(node) node$level))
  node.prob.rev <- as.probability(split.prob(node.level))
  node.prob <- as.probability(1/node.prob.rev)
  
  selected.node.idx <- sample(1:length(node.prob), 1, prob = node.prob)
  selected.node <- parent.of[[selected.node.idx]]
  
  if (Prune(selected.node, function(cn) F) != 2)
    warning("Pruned not 2")
  l.prob.sel.split <- CART.prob.rule(selected.node)
  selected.node$RemoveAttribute('split.rule')
  after.move(selected.node)
  
  return(list(
    tree.new=tree.new,
    l.prob=log(node.prob[selected.node.idx]),
    l.prob.rev=log(node.prob.rev[selected.node.idx])+l.prob.sel.split
  ))
}

CART.move.change <- function(tree, split.prob, after.move, augment.sampler, only.value=F, max.try=100) {
  if (tree$leafCount == 1) {
    return(NULL)
  }
  tree.new <- Clone(tree)
  if (!is.null(augment.sampler))
    Do(tree.new$leaves, augment.sampler)
  
  split.node <- Traverse(tree.new, filterFun = isNotLeaf)
  #node.level <- unname(vapply(split.node, function(node) node$level, numeric(1)))
  #node.prob <- as.probability(1/split.prob(node.level))
  
  if (length(split.node) == 0)
    return(NULL)
  ok <- F
  for (try in 1:max.try) {
    selected.node.idx <- sample(1:length(split.node), 1)#, prob = node.prob)
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
    subtree$Do(after.move, filterFun = isLeaf)
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

CART.move.swap <- function(tree, split.prob, after.move, augment.sampler, max.try=100) {
  if (tree$leafCount < 3) {
    return(NULL)
  }
  tree.new <- Clone(tree)
  if (!is.null(augment.sampler))
    Do(tree.new$leaves, augment.sampler)
  
  split.node <- Traverse(tree.new, filterFun = function(node) {
    isNotLeaf(node) && (isNotLeaf(node$children[[1]]) || isNotLeaf(node$children[[2]]))
  })
  #node.level <- vapply(split.node, function(node) node$level, numeric(1))
  #node.prob <- as.probability(1/split.prob(node.level))
  
  if (length(split.node) == 0)
    return(NULL)
  ok <- F
  for (try in 1:max.try) {
    selected.node.idx <- sample(1:length(split.node), 1)#, prob = node.prob)
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
    subtree$Do(after.move, filterFun = isLeaf)
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

