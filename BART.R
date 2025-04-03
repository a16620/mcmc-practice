library(data.tree)
library(pracma)
library(invgamma)
library(parallel)

BART <- function(x, y, m, k, tree.prior, sigma.param=c(3, .9), col.prior=NA, iteration, burn.in=0, compress.type='list') {
  #Y transform
  y.range <- max(y)-min(y)
  y.shift <- min(y)
  y.transform <- (y-y.shift)/y.range-0.5
  
  #Sigma prior
  sigma.est <- var(lm(y~., data=x)$residual)
  sigma.nu2 <- sigma.param[1]/2
  sigma.q <- sigma.param[2]
  
  nr.C <- sigma.q*gamma(sigma.nu2)
  nr.K <- sigma.nu2/sigma.est

  sigma.lambda <- 1
  for (i in 1:300) {
    lambda.grad <- (gammainc(nr.K*sigma.lambda, sigma.nu2)[2]-nr.C)/(sigma.lambda**(sigma.nu2-1)*nr.K**sigma.nu2*exp(-sigma.lambda*nr.K))
    if (abs(lambda.grad) < 1e-2)
      break
    sigma.lambda <- sigma.lambda + 0.1*lambda.grad
  }
  sigma.prior <- c(sigma.nu2, unname(sigma.lambda)*sigma.nu2)
  
  #Node prior
  node.prior <- c(0, (0.5/k)**2/m)
  
  y.len <- length(y)
  
  if (tree.prior[2] > 0)
    tree.prior[2] <- -tree.prior[2]
  
  if (compress.type == 'matrix') {
    col.dict <- setNames(1:ncol(x), colnames(x))
    compressor <- compress.DT2(col.dict)
  } else if (compress.type == 'list') {
    compressor <- compress.DT
  } else {
    stop('compressor type error')
  }
  
  #MCMC
  R0 <- rep(0, y.len)
  mcmc.sigma <- rinvgamma(1, sigma.prior[1], sigma.prior[2])
  mcmc.trees <- list()
  
  if (is.na(col.prior)) {
    MC.MOVE.GROW <- BART.move.grow
    MC.MOVE.PRUNE <- BART.move.prune
    MC.MOVE.CHANGE <- BART.move.change
  } else {
    MC.MOVE.GROW <- BART.move.grow.diric
    MC.MOVE.PRUNE <- BART.move.prune.diric
    MC.MOVE.CHANGE <- BART.move.change.diric
  }
  
  for (i in 1:m) {
    mcmc.trees <- append(mcmc.trees, BART.create.single(R0, x, node.prior, mcmc.sigma, col.prior))
  }
  
  y.est <- vapply(mcmc.trees, BART.get.pred, numeric(y.len))
  y.residual.total <- y.transform-rowSums(y.est)
  
  for (iter in 1:burn.in) {
    for (k in 1:m) {
      y.residual <- y.residual.total+y.est[,k]
      tree.k <- mcmc.trees[[k]]
      
      mc.move <- sample(1:3, 1)
      if (mc.move == 1) { #grow
        mc.proposal <- MC.MOVE.GROW(tree.k, x, y.residual, tree.prior, node.prior, mcmc.sigma)
      } else if (mc.move == 2) { #prune
        mc.proposal <- MC.MOVE.PRUNE(tree.k, x, y.residual, tree.prior, node.prior, mcmc.sigma)
      } else if (mc.move == 3) { #change
        mc.proposal <- MC.MOVE.CHANGE(tree.k, x, y.residual, node.prior, mcmc.sigma)
      }
      
      tree.k.next <- mc.proposal$proposal
      mcmc.trees[[k]] <- tree.k.next
      if (log(runif(1)) <= mc.proposal$l.ratio) {
        y.est.k.new <- BART.get.pred(tree.k.next)
        y.residual.total <- y.residual.total + y.est[,k] - y.est.k.new
        y.est[,k] <- y.est.k.new
      }
    }
    
    mcmc.sigma <- BART.sample.sigma(m, y.residual.total, sigma.prior)
  }
  
  tree.store <- list()
  tree.compressed <- lapply(mcmc.trees, compressor)
  
  sigma.trace <- numeric(iteration)
  accept.count <- numeric(m)
  for (iter in 1:iteration) {
    for (k in 1:m) {
      y.residual <- y.residual.total+y.est[,k]
      tree.k <- mcmc.trees[[k]]
      
      mc.move <- sample(1:3, 1)
      if (mc.move == 1) { #grow
        mc.proposal <- BART.move.grow(tree.k, x, y.residual, tree.prior, node.prior, mcmc.sigma)
      } else if (mc.move == 2) { #prune
        mc.proposal <- BART.move.prune(tree.k, x, y.residual, tree.prior, node.prior, mcmc.sigma)
      } else if (mc.move == 3) { #change
        mc.proposal <- BART.move.change(tree.k, x, y.residual, node.prior, mcmc.sigma)
      }
      
      tree.k.next <- mc.proposal$proposal
      mcmc.trees[[k]] <- tree.k.next

      if (log(runif(1)) <= mc.proposal$l.ratio) {
        accept.count[k] <- accept.count[k] + 1
        y.est.k.new <- BART.get.pred(tree.k.next)
        y.residual.total <- y.residual.total + y.est[,k] - y.est.k.new
        y.est[,k] <- y.est.k.new
        
        tree.compressed[[k]] <- compressor(tree.k.next)
      }
    }
    
    mcmc.sigma <- BART.sample.sigma(m, y.residual.total, sigma.prior)
    
    sigma.trace[iter] <- mcmc.sigma
    tree.store[[iter]] <- tree.compressed
  }
  
  list(
    trees=tree.store,
    sigma.trace=sigma.trace,
    y.shift=y.shift,
    y.range=y.range
    )
}

Predict.compress <- function(bart, x, raw.mat=F, parallel=F) {
  if (parallel) {
    n.core <- parallel::detectCores()
    cl <- parallel::makeCluster(n.core)
    clusterExport(cl, varlist = c("x"), envir = environment())
    
    batch.pred <- parallel::parSapply(cl, bart$trees, function(trees) {
      #한 세트의 트리 각각 예측값을 얻는다.
      rowSums(vapply(trees, function(tree) {
        apply(x, 1, function(row) {
          cursor <- tree
          while (is.null(cursor$node.param)) {
            rule <- cursor$split.rule
            if (row[names(rule)] <= rule)
              cursor <- cursor$left
            else
              cursor <- cursor$right
          }
          return(cursor$node.param)
        })
      }, numeric(nrow(x))))
    })
    
    parallel::stopCluster(cl)
  } else {
    batch.pred <- vapply(bart$trees, function(trees) {
      #한 세트의 트리 각각 예측값을 얻는다.
      rowSums(vapply(trees, function(tree) {
        apply(x, 1, function(row) {
          cursor <- tree
          while (is.null(cursor$node.param)) {
            rule <- cursor$split.rule
            if (row[names(rule)] <= rule)
              cursor <- cursor$left
            else
              cursor <- cursor$right
          }
          return(cursor$node.param)
        })
      }, numeric(nrow(x))))
    }, numeric(nrow(x)))
  }
  
  if (raw.mat) {
    return(batch.pred)
  }
  
  batch.summary <- t(apply(batch.pred, 1, function(obs.preds) {
    c(mean(obs.preds), quantile(obs.preds, probs = c(.05, .95), names=F))
  }))
  
  `colnames<-`((batch.summary+0.5)*bart$y.range+bart$y.shift, c('mean', '5%', '95%'))
}

Predict.compress2 <- function(bart, x, raw.mat=F, parallel=F) {
  if (parallel) {
    n.core <- parallel::detectCores()
    cl <- parallel::makeCluster(n.core)
    clusterExport(cl, varlist = c("x"), envir = environment())
    
    batch.pred <- parallel::parSapply(cl, bart$trees, function(trees) {
      #한 세트의 트리 각각 예측값을 얻는다.
      rowSums(vapply(trees, function(tree) {
        apply(x, 1, function(row) {
          cursor <- tree[1,]
          while (!is.na(cursor[1])) {
            if (row[cursor[1]] <= cursor[2])
              cursor <- tree[cursor[3],]
            else
              cursor <- tree[cursor[3]+1,]
          }
          return(cursor[2])
        })
      }, numeric(nrow(x))))
    })
    
    parallel::stopCluster(cl)
  } else {
    batch.pred <- vapply(bart$trees, function(trees) {
      #한 세트의 트리 각각 예측값을 얻는다.
      rowSums(vapply(trees, function(tree) {
        apply(x, 1, function(row) {
          cursor <- tree[1,]
          while (!is.na(cursor[1])) {
            if (row[cursor[1]] <= cursor[2])
              cursor <- tree[cursor[3],]
            else
              cursor <- tree[cursor[3]+1,]
          }
          return(cursor[2])
        })
      }, numeric(nrow(x))))
    }, numeric(nrow(x)))
  }
  
  if (raw.mat) {
    return(batch.pred)
  }
  
  batch.summary <- t(apply(batch.pred, 1, function(obs.preds) {
    c(mean(obs.preds), quantile(obs.preds, probs = c(.05, .95), names=F))
  }))
  
  `colnames<-`((batch.summary+0.5)*bart$y.range+bart$y.shift, c('mean', '5%', '95%'))
}

BART.create.single <- function(y, x, lik.prior, sigma, diric.prior) {
  tree <- Node$new("Root", obs.idx = 1:length(y))
  BART.sample.mu(tree, y, lik.prior, sigma)
  if (!is.na(diric.prior))
    tree$col.alpha <- setNames(rep(diric.prior, ncol(x)), colnames(x))
  return(tree)
}

BART.move.grow <- function(tree, x, y, tree.prior, lik.prior, sigma) {
  tree.proposed <- Clone(tree)
  leaves <- tree.proposed$leaves
  n.grow <- length(leaves)
  can.grow <- leaves[vapply(leaves, function(node) {length(node$obs.idx) >= 2}, logical(1))]
  
  len.can.grow <- length(can.grow)
  while (len.can.grow > 0) {
    grow.node.idx <- sample(1:len.can.grow, 1)
    grow.node <- can.grow[[grow.node.idx]]
    
    rule.candidate <- rule.select(grow.node, x)
    
    len.rule.candidate <- length(rule.candidate)
    if (len.rule.candidate > 0) {
      rule.selected <- sample.safe.name(rule.candidate, 1)
      grow.node$split.rule <- rule.selected
      break
    }
    can.grow <- can.grow[-grow.node.idx]
    len.can.grow <- len.can.grow-1
  }
  
  if (len.can.grow > 0) {
    rule.string <- rule.to.name(rule.selected)
    rule.child.idx <- rule.split.idx(grow.node$split.rule, x, grow.node$obs.idx)
    
    grow.node$AddChild("left", rule.name=rule.string$left, obs.idx=rule.child.idx$left)
    grow.node$AddChild("right", rule.name=rule.string$right, obs.idx=rule.child.idx$right)
    
    n.prune <- sum(tree.proposed$Get(function(node) {
      isNotLeaf(node) && isLeaf(node$children[[1]]) && isLeaf(node$children[[2]])
    }, simplify = "array"))
    
    # rev = 1/n.prune for = 1/n.grow * 1/n.col * 1/n.row
    l.proposal <- log(n.grow)+log(ncol(x))+log(length(grow.node$obs.idx))-log(n.prune)
    #prior=(a, b), split(d)/(1-split(d)) * (1-split(d+1))^2
    grow.depth <- grow.node$level
    l.prior <- log(tree.prior[1])+tree.prior[2]*log(grow.depth)-log1p(-tree.prior[1]*grow.depth**tree.prior[2])+2*log1p(-tree.prior[1]*(grow.depth+1)**tree.prior[2])
    
    l.lik <- llik.node(grow.node, y, lik.prior, sigma)-llik.node(grow.node$children, y, lik.prior, sigma)
    
    BART.sample.mu(grow.node$children, y, lik.prior, sigma)
    grow.node$RemoveAttribute('node.param')
    #grow.node$RemoveAttribute('obs.idx')
  } else {
    l.proposal <- -Inf
    l.prior <- 0
    l.lik <- 0
  }
  
  list(
    l.ratio=l.proposal+l.prior+l.lik,
    proposal=tree.proposed
  )
}

BART.move.grow.diric <- function(tree, x, y, tree.prior, lik.prior, sigma) {
  tree.proposed <- Clone(tree)
  leaves <- tree.proposed$leaves
  n.grow <- length(leaves)
  can.grow <- leaves[vapply(leaves, function(node) {length(node$obs.idx) >= 2}, logical(1))]
  
  len.can.grow <- length(can.grow)
  while (len.can.grow > 0) {
    grow.node.idx <- sample(1:len.can.grow, 1)
    grow.node <- can.grow[[grow.node.idx]]
    
    rule.candidate <- rule.select(grow.node, x)
    
    len.rule.candidate <- length(rule.candidate)
    if (len.rule.candidate > 0) {
      rule.selected.dir <- sample.with.diric(rule.candidate, tree.proposed$col.alpha[names(rule.candidate)], sum.a=sum(tree.proposed$col.alpha))
      rule.selected <- rule.selected.dir$val
      
      grow.node$split.rule <- rule.selected
      break
    }
    can.grow <- can.grow[-grow.node.idx]
    len.can.grow <- len.can.grow-1
  }
  
  if (len.can.grow > 0) {
    rule.string <- rule.to.name(rule.selected)
    rule.child.idx <- rule.split.idx(grow.node$split.rule, x, grow.node$obs.idx)
    
    grow.node$AddChild("left", rule.name=rule.string$left, obs.idx=rule.child.idx$left)
    grow.node$AddChild("right", rule.name=rule.string$right, obs.idx=rule.child.idx$right)
    
    n.prune <- sum(tree.proposed$Get(function(node) {
      isNotLeaf(node) && isLeaf(node$children[[1]]) && isLeaf(node$children[[2]])
    }, simplify = "array"))
    
    cn <- names(rule.selected)
    tree.proposed$col.alpha[cn] <- tree.proposed$col.alpha[cn] + 1
    
    # rev = 1/n.prune for = 1/n.grow * 1/n.col * 1/n.row
    l.proposal <- log(n.grow)-rule.selected.dir$prob+log(length(grow.node$obs.idx))-log(n.prune)
    #prior=(a, b), split(d)/(1-split(d)) * (1-split(d+1))^2
    grow.depth <- grow.node$level
    l.prior <- log(tree.prior[1])+tree.prior[2]*log(grow.depth)-log1p(-tree.prior[1]*grow.depth**tree.prior[2])+2*log1p(-tree.prior[1]*(grow.depth+1)**tree.prior[2])
    
    l.lik <- llik.node(grow.node, y, lik.prior, sigma)-llik.node(grow.node$children, y, lik.prior, sigma)
    
    BART.sample.mu(grow.node$children, y, lik.prior, sigma)
    grow.node$RemoveAttribute('node.param')
    #grow.node$RemoveAttribute('obs.idx')
  } else {
    l.proposal <- -Inf
    l.prior <- 0
    l.lik <- 0
  }

  list(
    l.ratio=l.proposal+l.prior+l.lik,
    proposal=tree.proposed
  )
}

BART.move.prune <- function(tree, x, y, tree.prior, lik.prior, sigma) {
  tree.proposed <- Clone(tree)
  can.prune <- Traverse(tree.proposed, filterFun = function(node) {
    isNotLeaf(node) && isLeaf(node$children[[1]]) && isLeaf(node$children[[2]])
  })
  n.prune <- length(can.prune)
  if (n.prune == 0) {
    return(list(
      l.ratio=-Inf,
      proposal=tree
    ))
  }
  
  prune.node <- sample(can.prune, 1)[[1]]
  
  #prune.node$obs.idx <- c(prune.node$children$left$obs.idx, prune.node$children$right$obs.idx)
  grow.depth <- prune.node$level
  
  l.proposal <- log(n.prune)-log(tree.proposed$leafCount-1)-log(ncol(x))-log(length(prune.node$obs.idx))
  l.prior <- -(log(tree.prior[1])+tree.prior[2]*log(grow.depth)-log1p(-tree.prior[1]*grow.depth**tree.prior[2])+2*log1p(-tree.prior[1]*(grow.depth+1)**tree.prior[2]))
  
  l.lik <- llik.node(prune.node$children, y, lik.prior, sigma)-llik.node(prune.node, y, lik.prior, sigma)
  
  prune.node$RemoveAttribute('split.rule')
  prune.node$RemoveChild('left')
  prune.node$RemoveChild('right')
  BART.sample.mu(prune.node, y, lik.prior, sigma)
  
  list(
    l.ratio=l.proposal+l.prior+l.lik,
    proposal=tree.proposed
  )
}

BART.move.prune.diric <- function(tree, x, y, tree.prior, lik.prior, sigma) {
  tree.proposed <- Clone(tree)
  can.prune <- Traverse(tree.proposed, filterFun = function(node) {
    isNotLeaf(node) && isLeaf(node$children[[1]]) && isLeaf(node$children[[2]])
  })
  n.prune <- length(can.prune)
  if (n.prune == 0) {
    return(list(
      l.ratio=-Inf,
      proposal=tree
    ))
  }
  
  prune.node <- sample(can.prune, 1)[[1]]

  #prune.node$obs.idx <- c(prune.node$children$left$obs.idx, prune.node$children$right$obs.idx)
  grow.depth <- prune.node$level
  
  cn <- names(prune.node$split.rule)
  c.idx <- which(colnames(x) == cn)[1]
  tree.proposed$col.alpha[cn] <- tree.proposed$col.alpha[cn] - 1
  
  l.proposal <- log(n.prune)-log(tree.proposed$leafCount-1)+prob.sample.diric(c.idx, tree.proposed$col.alpha)-log(length(prune.node$obs.idx))
  l.prior <- -(log(tree.prior[1])+tree.prior[2]*log(grow.depth)-log1p(-tree.prior[1]*grow.depth**tree.prior[2])+2*log1p(-tree.prior[1]*(grow.depth+1)**tree.prior[2]))
  
  l.lik <- llik.node(prune.node$children, y, lik.prior, sigma)-llik.node(prune.node, y, lik.prior, sigma)

  prune.node$RemoveAttribute('split.rule')
  prune.node$RemoveChild('left')
  prune.node$RemoveChild('right')
  BART.sample.mu(prune.node, y, lik.prior, sigma)
    
  list(
    l.ratio=l.proposal+l.prior+l.lik,
    proposal=tree.proposed
  )
}

BART.move.change <- function(tree, x, y, lik.prior, sigma) {
  tree.proposed <- Clone(tree)
  can.change <- Traverse(tree.proposed, filterFun = isNotLeaf)
  len.can.change <- length(can.change)
  
  while (len.can.change > 0) {
    change.node.idx <- sample(1:len.can.change, 1)
    change.node <- can.change[[change.node.idx]]
    rule.rollback <- change.node$split.rule
    change.node.llik.before <- llik.node(change.node$leaves, y, lik.prior, sigma)
    
    for (try in 1:5) {
      rule.candidate <- rule.select(change.node, x)
      rule.shuffle <- sample.safe.name(rule.candidate, length(rule.candidate))
      
      ok <- F
      for (i in 1:length(rule.shuffle)) {
        change.node$split.rule <- rule.shuffle[i]
        if (update.subtree.obs(change.node, x)) {
          BART.sample.mu(change.node$leaves, y, lik.prior, sigma)
          
          rule.string <- rule.to.name(rule.shuffle[i])
          change.node$left$rule.name <- rule.string$left
          change.node$right$rule.name <- rule.string$right
          ok <- T
          break
        }
      }
      if (ok)
        break
    }
    
    if (ok)
      break
    
    change.node$split.rule <- rule.rollback
    update.subtree.obs(change.node, x)
    
    can.change <- can.change[-change.node.idx]
    len.can.change <- len.can.change-1
  }
  if (len.can.change == 0) {
    l.lik <- -Inf
  } else
    l.lik <- llik.node(change.node$leaves, y, lik.prior, sigma)-change.node.llik.before
  
  list(
    l.ratio=l.lik,
    proposal=tree.proposed
  )
}

BART.move.change.diric <- function(tree, x, y, lik.prior, sigma) {
  tree.proposed <- Clone(tree)
  can.change <- Traverse(tree.proposed, filterFun = isNotLeaf)
  len.can.change <- length(can.change)
  
  while (len.can.change > 0) {
    change.node.idx <- sample(1:len.can.change, 1)
    change.node <- can.change[[change.node.idx]]
    rule.rollback <- change.node$split.rule
    change.node.llik.before <- llik.node(change.node$leaves, y, lik.prior, sigma)
    
    alpha.before <- tree.proposed$col.alpha - as.integer(colnames(x) == names(rule.rollback))
    
    for (try in 1:50) {
      rule.candidate <- rule.select(change.node, x)
      rule.shuffle.dir <- sample.with.diric(rule.candidate, alpha.before[names(rule.candidate)], length(rule.candidate), sum(alpha.before))
      rule.shuffle <- rule.shuffle.dir$val
      
      ok <- F
      for (i in 1:length(rule.shuffle)) {
        change.node$split.rule <- rule.shuffle[i]
        if (update.subtree.obs(change.node, x)) {
          BART.sample.mu(change.node$leaves, y, lik.prior, sigma)
          
          rule.string <- rule.to.name(rule.shuffle[i])
          change.node$left$rule.name <- rule.string$left
          change.node$right$rule.name <- rule.string$right
          ok <- T
          break
        }
      }
      if (ok)
        break
    }
    
    if (ok)
      break
    
    change.node$split.rule <- rule.rollback
    update.subtree.obs(change.node, x)
    
    can.change <- can.change[-change.node.idx]
    len.can.change <- len.can.change-1
  }
  if (len.can.change == 0) {
    l.lik <- -Inf
    l.proposal <- 0
  } else {
    l.lik <- llik.node(change.node$leaves, y, lik.prior, sigma)-change.node.llik.before
    c.idx.reverse <-  which(colnames(x) == names(rule.rollback))[1]
    l.proposal <- prob.sample.diric(c.idx.reverse, alpha.before)-rule.shuffle.dir$prob[i]
    
    cn <- names(change.node$split.rule)
    alpha.before[cn] <- alpha.before[cn] + 1
    tree.proposed$col.alpha <- alpha.before
  }
  
  list(
    l.ratio=l.proposal+l.lik,
    proposal=tree.proposed
  )
}

BART.sample.mu <- function(nodes, y, lik.prior, sigma) {
  if (!is.list(nodes))
    nodes <- list(nodes)
  Do(nodes, function(node) {
    y.t <- y[node$obs.idx]
    n.t <- length(y.t)
    
    post.var <- 1/(n.t/sigma+1/lik.prior[2])
    post.mean <- post.var*(mean(y.t)*n.t/sigma+lik.prior[1]/lik.prior[2])

    node$node.param <- rnorm(1, post.mean, sqrt(post.var))
  })
}

BART.sample.sigma <- function(M, y.residual, sigma.prior) {
  par1 <- sigma.prior[1]+M*length(y.residual)/2
  par2 <- sigma.prior[2]+M*as.numeric(crossprod(y.residual))
  
  rinvgamma(1, par1, scale=par2)
}

BART.get.pred <- function(tree) {
  R <- numeric(length(tree$root$obs.idx))
  Do(tree$leaves, function(node) {
    R[node$obs.idx] <<- node$node.param
  })
  return(R)
}

update.subtree.obs <- function(tree, x) {
  tree$Do(function(node) {
    rule.child.idx <- rule.split.idx(node$split.rule, x, node$obs.idx)
    node$left$obs.idx <- rule.child.idx$left
    node$right$obs.idx <- rule.child.idx$right
  }, filterFun = isNotLeaf, traversal = "level")
  all(vapply(tree$leaves, function(node) length(node$obs.idx)>0, logical(1)))
}

llik.node <- function(nodes, y, lik.prior, sigma) {
  if (!is.list(nodes)) #단일 노드
    nodes <- list(nodes)
  sum(vapply(nodes, function(node) {
    dnorm(mean(y[node$obs.idx]), lik.prior[1], sqrt(sigma/length(node$obs.idx)+lik.prior[2]), log=T)
  }, numeric(1)))
}

rule.select <- function(node, x) {
  #Select rule
  rule.candidate <- setNames(vapply(colnames(x), function(cn) {
    split.values <- unique(x[node$obs.idx,cn])
    split.values <- split.values[split.values < max(split.values)]
    
    if (length(split.values) == 0)
      return(NA)
    else if (length(split.values) == 1)
      rule.value <- split.values
    else {
      smp <- sample(split.values, 2, replace = T)
      l.bound <- smp[1]; u.bound <- smp[2]
      if (l.bound > u.bound) {
        ll <- u.bound
        u.bound <- l.bound
        l.bound <- ll
      }
      
      aux.param <- mean(split.values[l.bound <= split.values & split.values <= u.bound] >= mean(smp))
      u.aux <- rbeta(1, 1-aux.param, aux.param)
      rule.value <- l.bound+(u.bound-l.bound)*u.aux
    }
    
    return(rule.value)
  }, numeric(1)), colnames(x))
  rule.candidate[!is.na(rule.candidate)]
}

rule.to.name <- function(rule) {
  cn <- names(rule)
  val <- trimws(format(rule, nsmall=3), whitespace = '.0')
  list(left=paste0(cn,'<=',val),right=paste0(cn,'>',val))
}

rule.split.idx <- function(rule, x, idx) {
  is.left <- x[idx,names(rule)] <= rule
  list(left=idx[is.left], right=idx[!is.left])
}

sample.safe.name <- function(x, n) {
  x[sample(1:length(x), n)]
}

compress.DT <- function(tree) {
  as.list(tree, keepOnly = c('split.rule', 'rule.name', 'node.param'))
}

compress.DT2 <- function(col.dict) {
  compressor <- function(tree) {
    trv <- Traverse(tree, traversal = "level")
    Set(trv, nid=1:length(trv))
    len.trv <- length(trv)
    
    split.col  <- numeric(len.trv)
    split.val  <- numeric(len.trv)
    child.left <- integer(len.trv)
    
    for (nid in 1:len.trv) {
      node <- trv[[nid]]
      if (isNotLeaf(node)) {
        split.col[nid]  <- col.dict[names(node$split.rule)]
        split.val[nid]  <- node$split.rule
        child.left[nid] <- node$left$nid
      } else {
        split.col[nid]  <- NA
        split.val[nid]  <- node$node.param
        child.left[nid] <- NA
      }
      
      node$RemoveAttribute('nid')
    }
    return(cbind(split.col,split.val,child.left))
  }
}

sample.with.diric <- function(x, a, n=1, sum.a=NA) {
  y <- rgamma(length(a), a, 1)
  prob.from.diric <- y/sum(y)
  
  x.idx <- sample(1:length(x), n, prob = prob.from.diric)
  if (is.na(sum.a))
    sum.a <- sum(a)
  
  list(val=x[x.idx], prob=lgamma(a[x.idx]+1)-lgamma(sum.a+1)+lgamma(sum.a)-lgamma(a[x.idx]))
}

prob.sample.diric <- function(x.idx, a) {
  sum.a <- sum(a)
  lgamma(a[x.idx]+1)-lgamma(sum.a+1)+lgamma(sum.a)-lgamma(a[x.idx])
}
