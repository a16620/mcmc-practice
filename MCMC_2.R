library(stats)
library(Rfast)
library(dplyr)

library(mvtnorm)

library(reshape2)
library(progress)
library(ggplot2)
library(ggpubr)

############# Simple MH sampler #############

bayes.single.mh.sampler <- function (target_density, proposal, initial.value, limit=Inf) {
  c <- 2.4
  rej <- 0
  while (limit > 0) {
    next_param <- proposal$sampler(initial.value, c)
    alpha <- target_density(next_param)/target_density(initial.value) * proposal$density(next_param, initial.value, c)/proposal$density(initial.value, next_param, c)
    if (alpha >= runif(1)) {
      return(next_param)
    } else {
      rej <- rej + 1
      if (rej > 20) {
        rej <- 0
        c <- c * 1.414
      }
    }
    limit <- limit-1
  }
  return(NA)
}

############# Sampler components #############

# (density(x, mean, c), sampler(mean, c))
make_proposal <- function(density., sampler=NULL, positive=F) {
  if (is.null(sampler)) {
    sampler <- function(mean, c) {
      return(bayes.single.mh.sampler(
        target_density = function(x) density.(x, mean, c),
        proposal = make_proposal_normal(),
        initial.value = ifelse(positive, rchisq(1, 11.52), rnorm(1, sd=11.52))
      ))
    }
  }
  return(list(
    density=density.,
    sampler=sampler
  ))
}

make_proposal_normal <- function(sd=1) {
  return(
    make_proposal(function(x, mean, c) dnorm(x, mean=mean, sd=sd*c), sampler=function(mean, c) rnorm(1, mean, sd*c))
  )
}

make_proposal_mvnormal <- function(d, cov=NA) {
  if (any(is.na(cov)))
    cov <- diag(rep(1, d))
  return(
    make_proposal(function(x, mean, c) mvtnorm::dmvnorm(x=x, mean=mean, sigma = cov*c**2), sampler=function(mean, c) mvtnorm::rmvnorm(1, mean=mean, sigma=cov*c**2))
  )
}

# (density(x), sampler())
make_prior <- function(density., sampler=NULL, positive=F) {
  if (is.null(sampler)) {
    sampler <- function() {
      return(bayes.single.mh.sampler(
        target_density = density.,
        proposal = make_proposal_normal(),
        initial.value = ifelse(positive, rchisq(1, 11.52), rnorm(1))
          ))
    }
  }
  return(list(
    density=density.,
    sampler=sampler
  ))
}

make_prior_normal <- function(mean=0, sd=1) {
  return(
    make_prior(function(x) dnorm(x, mean, sd), function() rnorm(1, mean, sd))
  )
}

make_prior_mvnormal <- function(d, mean=0, cov=NA) {
  if (any(is.na(cov)))
    cov <- diag(rep(1, d))
  if (length(mean) != d)
    mean <- rep(mean[1], d)
  return(
    make_prior(function(x) mvtnorm::dmvnorm(x=x, mean=mean, sigma = cov), function() mvtnorm::rmvnorm(1, mean=mean, sigma=cov))
  )
}

make_prior_hcauchy <- function(mean=0, sd=1) {
  return(
    make_prior(function(x) {
      ifelse(x>=mean, 0.636/(sd*(1+((x-mean)/sd)**2)), 0)
    }, positive = T)
  )
}

############# GLM MH sampler #############

#Data X,Y, likelihood(y_vector, mu_vector), prior list for beta, proposal list for beta
bayes.glm.mh.sampler <- function(x, y, likelihood, link.fn, prior, proposal, initial.beta=NULL, iter. = 2000, C.mult=2.4) {
  beta.order <- colnames(x)
  beta.count <- length(beta.order)
  beta.chain <- matrix(0L, ncol = beta.count, nrow = iter., byrow = T)
  colnames(beta.chain) <- beta.order
  
  #Mapping priors
  if (is.null(prior))
    prior <- list()
  if (is.null(prior$.mcmc.default))
    prior$.mcmc.default <- make_prior_normal()
  
  beta.prior <- lapply(beta.order, function(name) {
    pri <- prior[[name]]
    if (is.null(pri))
      return(prior$.mcmc.default)
    else
      return(pri)
  })
  
  #Mapping proposals
  if (is.null(proposal))
    proposal <- list()
  if (is.null(proposal$.mcmc.default))
    proposal$.mcmc.default <- make_proposal_normal()
  
  beta.proposal <- lapply(beta.order, function(name) {
    if (is.null(proposal[[name]]))
      return(proposal$.mcmc.default)
    else
      return(proposal[[name]])
  })
  
  #Set initial values
  if (is.null(initial.beta)) {
    beta.chain[1,] <- unlist(lapply(beta.prior, function(pri) pri$sampler()))
  } else if (is.list(initial.beta)) {
    beta.chain[1,] <- unlist(lapply(beta.order, function(b.name) {
      b.point <- initial.beta[[b.name]]
      if (is.null(b.point)) {
        pri <- beta.prior[[b.name]]
        return(ifelse(is.null(pri), beta.prior$.mcmc.default, pri)$sampler())
      }
    }))
  } else {
    beta.chain[1,] <- initial.beta
  }
  
  accept.count <- rep(1, beta.count)
  names(accept.count) <- beta.order
  accept.count.batch <- rep(1, beta.count)
  
  accept.interval <- 200
  C.min <- 0.001
  C.max <- 1000
  
  C.adj.rate  <- 0.23+0.6*exp(-beta.count)
  C.adj.lower <- C.adj.rate-0.015
  C.adj.upper <- C.adj.rate+0.015
  
  sampler.C <- rep(2.4/sqrt(beta.count), beta.count)
  #alpha.chain <- matrix(F, ncol = beta.count, nrow = iter., byrow = T)
  
  #Sampling
  warm.up <- T
  iter <- 2
  while (iter <= iter.) {
    beta.chain[iter,] <- beta.chain[iter-1,]
    aux.U <-  runif(beta.count) #한번에 샘플링해서 속도 향상
    for (i in 1:beta.count) {
      mcmc.proposal <- beta.proposal[[i]]
      mcmc.prior <- beta.prior[[i]]
      
      next.beta <- mcmc.proposal$sampler(beta.chain[iter, i], sampler.C[i])
      if (is.na(next.beta)) {
        print(beta.chain[iter,i])
        return()
      }
      prev.beta <- beta.chain[iter, i]
      
      mu.prev <- link.fn(rowSums(sweep(x, 2, beta.chain[iter,], `*`)))
      beta.chain[iter, i] <- next.beta
      mu.next <- link.fn(rowSums(sweep(x, 2, beta.chain[iter,], `*`)))
      
      log.ratio.likelihood <- likelihood(y, mu.next)-likelihood(y, mu.prev)
      log.ratio.prior <- log(mcmc.prior$density(next.beta))-log(mcmc.prior$density(prev.beta))
      log.ratio.proposal <- log(mcmc.proposal$density(prev.beta, next.beta, sampler.C[i]))-
        log(mcmc.proposal$density(next.beta, prev.beta, sampler.C[i]))
      
      log.alpha <- log.ratio.likelihood+log.ratio.prior+log.ratio.proposal
      if (is.na(log.ratio.prior) | is.na(log.ratio.likelihood)) {
        log.alpha <- -Inf
      }
      if (exp(min(log.alpha, 0)) >= aux.U[i]) { #Accpet
        accept.count.batch[i] <- accept.count.batch[i]+1
      } else {
        beta.chain[iter, i] <- prev.beta
      }
    }

    #Adjust sampling C
    if (iter %% accept.interval == 0) {
      a.rate <- accept.count.batch/accept.interval
      if (T) {
        sampler.C <- sampler.C * ifelse(a.rate < C.adj.lower, 0.7+0.3*a.rate/C.adj.lower, 1)*
          ifelse(a.rate > C.adj.upper, 1+0.4*(a.rate-C.adj.upper)/(1-C.adj.upper), 1)
      } else if (warm.up) {
        C.adaptive <- C.mult * sampler.C
        C.underflow <- min(sampler.C, 1) * 0.00001
        warm.up <- F
      } else {
        sampler.C <- C.adaptive*(colVars(beta.chain[(iter-accept.interval/2+1):iter, ,drop=F], std=F)+C.underflow)
      }
      accept.count <- accept.count + accept.count.batch
      accept.count.batch <- rep(0, beta.count)
    }
    
    iter <- iter+1
  }
  
  accept.count <- accept.count + accept.count.batch
  return(list(
    accept.rate=accept.count/iter.,
    params=beta.chain
  ))
}

bayesGLM <- function(formula., data, family, prior, proposal=NULL, chain.cnt=4, iteration=10000, cheat.beta=F, ...) {
  response <- data[, all.vars(formula.)[1]]
  deps <- model.matrix(formula., data)
  print(paste0("Model: ", format(formula.)), quote=F)
  
  density_dict <- list(
    gaussian =function(x,mean) sum(dnorm(x,mean, log = T)),
    binomial =function(x,mean) sum(dbinom(x, 1, mean, log = T)),
    poisson  =function(x,mean) sum(dpois(x, mean, log = T))
  )
  
  density.fn <- density_dict[[family$family]]
  
  initial.beta <- NULL
  if (cheat.beta) {
    cheat.glm <- glm(formula., data=data, family = family)
    initial.beta <- cheat.glm$coefficient
    if (is.null(proposal)) {
      proposal <- lapply(as.list(diag(vcov(cheat.glm))), function(V) make_proposal_normal(sd=sqrt(V)))
    }
  }
  
  #M-H
  #Output=b1:b2...:iter:chain
  chains <- data.frame()
  accept.rates <- data.frame()
  
  pb <- progress_bar$new(format = paste0("MCMC: [:bar] Chain: :current/:total time:[:elapsed/:eta]"),
                         total=chain.cnt, clear = F, width=60, show_after = 0)
  pb$tick(len=0)
  
  for (i in 1:chain.cnt) {
    chain.i <- suppressWarnings(bayes.glm.mh.sampler(deps, response, density.fn, family$linkinv, prior, proposal, iter.=iteration, initial.beta = initial.beta, ...))
    chains <- rbind(chains, cbind(chain.i$params, data.frame(.iter=1:iteration, .chain=i)))
    accept.rates <- rbind(accept.rates, chain.i$accept.rate)
    pb$tick()
  }
  
  deps.names <- colnames(deps)
  colnames(accept.rates) <- c(deps.names)
  accept.rates$.chain <- 1:chain.cnt
  
  print("Complete.", quote=F)
  print("Accept rates:", quote=F)
  print(accept.rates)
  
  return(chains)
}

############# General MH sampler #############
# blocks {
#   {
#     param: vector(string)
#     initial.value: vector(double)
#     prior, proposal
#     isGibbs
#   }
# }


bayes.block.mh.sampler <- function(blocks, likelihood, burn.in = 1000, iter.out = 2000, thin. = 1) {
  param.names <- NULL
  param.count.block <- NULL
  seq.block <- seq_along(blocks)
  for (b.idx in seq.block) {
    block <- blocks[[b.idx]]
    param.names <- append(param.names, block$param)
    param.count.block <- append(param.count.block, length(block$param))
  }
  
  if (is.null(param.names)) {
    print("이름을 찾을수 없음. 블록의 'param' attribute에 벡터로 지정하십시오.", quote=F)
    return()
  }
  
  param.count <- length(param.names)
  param.chain <- matrix(0.0, ncol=param.count, nrow = iter.out)
  colnames(param.chain) <- param.names
  
  block.count <- length(seq.block)
  block.isGibbs <- NULL
  
  param.current <- setNames(rep(0, param.count), param.names)
  for (b.idx in seq.block) {
    block <- blocks[[b.idx]]
    block.isGibbs <- append(block.isGibbs, ifelse(is.null(block$isGibbs), F, block$isGibbs))
    if (is.null(block$initial.value)) {
      param.current[block$param] <- block$prior$sampler()
    } else {
      param.current[block$param] <- block$initial.value
    }
  }
  
  sampler.C <- 2.4/param.count.block#rep(2.4/sqrt(param.count), block.count)
  sampler.C[block.isGibbs] <- 1
  
  accept.count <- rep(0, block.count)
  accept.count.batch <- rep(1, block.count)
  
  C.adjust.interval <- 200
  C.adj.rate  <- 0.23+0.6*exp(-param.count)
  C.adj.lower <- C.adj.rate-0.015
  C.adj.upper <- C.adj.rate+0.015
  
  #Sampling
  iter.total <- burn.in+iter.out*thin.
  
  iter <- 2
  chain.counter <- 1
  while (iter <= iter.total) {
    aux.U <-  runif(block.count)
    for (b.idx in seq.block) {
      block <- blocks[[b.idx]]
      
      if (block.isGibbs[b.idx]) {
        next.param <- block$proposal$sampler(param.current)
        param.current[block$param] <- next.param
        
        accept.count.batch[b.idx] <- accept.count.batch[b.idx]+1
      } else {
        C.current <- sampler.C[b.idx]
        prev.param <- param.current
        prev.param.subset <- prev.param[block$param]
        next.param <- block$proposal$sampler(prev.param.subset, C.current)
        param.current[block$param] <- next.param
        
        log.ratio.likelihood <- likelihood(param.current)-likelihood(prev.param)
        log.ratio.prior <- log(block$prior$density(next.param))-log(block$prior$density(prev.param.subset))
        log.ratio.proposal <- log(block$proposal$density(prev.param.subset, next.param, C.current))-
          log(block$proposal$density(next.param, prev.param.subset, C.current))
        
        log.alpha <- log.ratio.likelihood+log.ratio.prior+log.ratio.proposal
        if (is.na(log.ratio.prior) | is.na(log.ratio.likelihood)) {
          log.alpha <- -Inf
        }
        if (min(exp(log.alpha),1) >= aux.U[b.idx]) {
          accept.count.batch[b.idx] <- accept.count.batch[b.idx]+1
        } else { #rollback
          param.current <- prev.param
        }
      }
    }
    
    if (iter > burn.in & iter %% thin. == 0) {
      param.chain[chain.counter,] <- param.current
      chain.counter <- chain.counter + 1
    }
    
    #Adjust C
    if (iter %% C.adjust.interval == 0) {
      a.rate <- accept.count.batch/C.adjust.interval
      sampler.C <- sampler.C * ifelse(a.rate < C.adj.lower, 0.65+0.3*a.rate/C.adj.lower, 1)*
        ifelse(a.rate > C.adj.upper, 1.05+0.4*(a.rate-C.adj.upper)/(1-C.adj.upper), 1)
      sampler.C[block.isGibbs] <- 1
      
      accept.count <- accept.count + accept.count.batch
      accept.count.batch <- rep(0, block.count)
    }
    
    iter <- iter+1
  }
  
  accept.count <- accept.count + accept.count.batch
  return(list(
    accept.rate=accept.count/iter.total,
    params=param.chain
  ))
}

bayesBlockedMH <- function(likelihood, blocks.info, iteration=10000, chain.cnt=4, ...) {
  chains <- data.frame()
  accept.rates <- data.frame()
  
  pb <- progress_bar$new(format = paste0("MCMC: [:bar] Chain: :current/:total time:[:elapsed/:eta]"),
                         total=chain.cnt, clear = F, width=60, show_after = 0)
  pb$tick(len=0)
  
  for (i in 1:chain.cnt) {
    chain.i <- suppressWarnings(bayes.block.mh.sampler(blocks.info, likelihood, iter.out=iteration, ...))
    chains <- rbind(chains, cbind(chain.i$params, data.frame(.iter=1:iteration, .chain=i)))
    accept.rates <- rbind(accept.rates, chain.i$accept.rate)
    pb$tick()
  }
  
  colnames(accept.rates) <- paste0("Block#", 1:length(blocks.info))
  accept.rates$.chain <- 1:chain.cnt
  
  print("Complete.", quote=F)
  print("Accept rates:", quote=F)
  print(accept.rates)
  
  return(chains)
}

make_block <- function(param.name, prior, proposal=NULL, isGibbs=F, initial.value=NULL) {
  if (is.null(proposal)) {
    dp <- length(param.name)
    if (dp == 1) {
      proposal <- make_proposal_normal()
    } else {
      proposal <- make_proposal_mvnormal(d = dp)
    }
  }
  return(list(
    param=param.name,
    prior=prior,
    proposal=proposal,
    isGibbs=isGibbs,
    initial.value=initial.value
  ))
}

############# Joint MH sampler #############

bayes.joint.mh.sampler <- function(param.name, likelihood, prior, proposal, initial.param=NULL, iter. = 2000, C.mult=1) {
  param.count <- length(param.name)
  param.chain <- matrix(0L, ncol = param.count, nrow = iter., byrow = T)
  colnames(param.chain) <- param.name
  
  if (is.null(proposal))
    proposal <- make_proposal_mvnormal(param.count)
  
  accept.count <- 1
  accept.interval <- 100
  C.min <- 0.001
  C.max <- 1000
  
  C.adj.rate  <- 0.23+0.6*exp(-param.count)
  C.adj.lower <- C.adj.rate-0.015
  C.adj.upper <- C.adj.rate+0.015
  
  sampler.C <- C.mult*2.4/sqrt(param.count)
  
  if (is.null(initial.param)) {
    current.param <- prior$sampler()
  } else {
    current.param <- initial.param
  }
  param.chain[1,] <- current.param
  log.likelihood.current <- likelihood(setNames(current.param, param.name))
  
  #Sampling
  iter <- 2
  while (iter <= iter.) {
    next.param <- setNames(proposal$sampler(current.param, sampler.C), param.name)
    log.likelihood.new <- likelihood(next.param)
    
    log.ratio.likelihood <- log.likelihood.new - log.likelihood.current
    log.ratio.prior <- log(prior$density(next.param)) - log(prior$density(current.param))
    log.ratio.proposal <- log(proposal$density(current.param, next.param, sampler.C))-
      log(proposal$density(next.param, current.param, sampler.C))
    
    log.alpha <- log.ratio.likelihood+log.ratio.prior+log.ratio.proposal
    if (is.na(log.ratio.prior) | is.na(log.ratio.likelihood)) {
      log.alpha <- -Inf
    }
    if (exp(min(log.alpha, 0)) >= runif(1)) {
      accept.count <- accept.count+1
      current.param <- next.param
      log.likelihood.current <- likelihood(current.param)
    }
    param.chain[iter,] <- current.param
    
    #Adjust sampling sd
    if (iter %% accept.interval == 0) {
      a.rate <- accept.count/iter
      if (a.rate < 0.23) {
        a.rate <- a.rate * 0.7
      } else if (a.rate > 0.44) {
        a.rate <- a.rate * 1.4
      }
    }
    
    iter <- iter+1
  }
  
  return(list(
    accept.rate=accept.count/iter.,
    params=param.chain
  ))
}

bayesJointMH <- function(param.name, likelihood, prior, proposal=NULL, chain.cnt=4, iteration=10000, ...) {
  chains <- data.frame()
  accept.rates <- data.frame()
  
  pb <- progress_bar$new(format = paste0("MCMC: [:bar] Chain: :current/:total time:[:elapsed/:eta]"),
                         total=chain.cnt, clear = F, width=60, show_after = 0)
  pb$tick(len=0)
  
  for (i in 1:chain.cnt) {
    chain.i <- suppressWarnings(bayes.joint.mh.sampler(param.name, likelihood, prior, proposal, iter.=iteration, ...))
    chains <- rbind(chains, cbind(chain.i$params, data.frame(.iter=1:iteration, .chain=i)))
    accept.rates <- rbind(accept.rates, chain.i$accept.rate)
    pb$tick()
  }
  
  colnames(accept.rates) <- "Accept rate"
  accept.rates$.chain <- 1:chain.cnt
  
  print("Complete.", quote=F)
  print("Accept rates:", quote=F)
  print(accept.rates)
  
  return(chains)
}

############# Hamiltonian MC #############

bayes.joint.hmc.sampler <- function(param.names, log.joint.posterior, gradient.posterior=NULL, moment.M=NULL, 
                                    delta.t=0.01, step.count=100, initial.value=NULL, burn.in = 100, iter.out = 1000, thin. = 1) {
  param.count <- length(param.names)
  param.chain <- matrix(0.0, ncol=param.count, nrow = iter.out)
  colnames(param.chain) <- param.names
  
  if (is.null(moment.M)) {
    moment.M <- diag(rep(1, param.count))
  } else {
    moment.M <- unname(moment.M)
  }
  
  if (is.null(gradient.posterior)) {
    grad.delta.h <- 0.001
    gradient.posterior <- function(param) {
      grads <- rep(0, param.count)
      for (i in 1:param.count) {
        p.rollback <- param[i]
        param[i] <- p.rollback+grad.delta.h
        delta.f <- log.joint.posterior(param)
        param[i] <- p.rollback-grad.delta.h
        delta.f <- delta.f - log.joint.posterior(param)
        param[i] <- p.rollback
        grads[i] <- delta.f/(2*grad.delta.h)
      }
      return(grads)
    }
  }

  if (is.null(initial.value)) {
    initial.value <- matrix(rep(0, param.count), nrow=1)
  }
  param.current <- matrix(initial.value, nrow=1)
  colnames(param.current) <- param.names
  
  M.inv <- solve(moment.M)
  
  accept.count <- 0
  accept.count.batch <- 1
  
  C.adjust.interval <- 200
  C.adj.rate  <- 0.23+0.6*exp(-param.count)
  C.adj.lower <- C.adj.rate-0.015
  C.adj.upper <- C.adj.rate+0.015
  
  #Sampling
  iter.total <- burn.in+iter.out*thin.
  
  iter <- 2
  chain.counter <- 1
  while (iter <= iter.total) {
    moment.current <- rmvnorm(1, sigma = moment.M)
    param.rollback <- param.current
    
    log.accept.F <- log.joint.posterior(param.current[1,])
    log.accept.m <- moment.current%*%M.inv%*%t(moment.current) #moment.current=1xdim(params)
    
    for (lf.i in 1:step.count) {
      moment.halfstep <- moment.current + 0.5*delta.t*gradient.posterior(param.current[1,])
      param.current <- param.current + delta.t*moment.halfstep%*%M.inv
      moment.current <- moment.halfstep + 0.5*delta.t*gradient.posterior(param.current[1,])
      if (any(is.na(param.current)))
        break
    }

    log.accept.F <- log.joint.posterior(param.current[1,])-log.accept.F
    log.accept.m <- 0.5*(log.accept.m-moment.current%*%M.inv%*%t(moment.current))
    log.accept <- log.accept.F+log.accept.m
    if (is.na(log.accept.F) | is.na(log.accept.m)) {
      log.accept <- -Inf
    }
    
    if (exp(min(log.accept, 0)) >= runif(1)) { #accept
      accept.count.batch <- accept.count.batch + 1
    } else {
      param.current <- param.rollback
    }
    
    if (iter > burn.in & iter %% thin. == 0) {
      param.chain[chain.counter,] <- param.current
      chain.counter <- chain.counter + 1
    }
    
    #Adjust C
    if (iter %% C.adjust.interval == 0) {
      #a.rate <- accept.count.batch/C.adjust.interval
      #sampler.C <- sampler.C * ifelse(a.rate < C.adj.lower, 0.65+0.3*a.rate/C.adj.lower, 1)*
      #  ifelse(a.rate > C.adj.upper, 1.05+0.4*(a.rate-C.adj.upper)/(1-C.adj.upper), 1)
      #sampler.C[block.isGibbs] <- 1
      
      accept.count <- accept.count + accept.count.batch
      accept.count.batch <- 0
    }
    
    iter <- iter+1
  }
  
  accept.count <- accept.count + accept.count.batch
  return(list(
    accept.rate=accept.count/iter.total,
    params=param.chain
  ))
}

bayesJointHMC <- function(param.name, log.posterior, chain.cnt=4, iteration=1000, ...) {
  chains <- data.frame()
  accept.rates <- data.frame()
  
  pb <- progress_bar$new(format = paste0("MCMC: [:bar] Chain: :current/:total time:[:elapsed/:eta]"),
                         total=chain.cnt, clear = F, width=60, show_after = 0)
  pb$tick(len=0)
  
  for (i in 1:chain.cnt) {
    chain.i <- suppressWarnings(bayes.joint.hmc.sampler(param.name, log.posterior, iter.out=iteration, ...))
    chains <- rbind(chains, cbind(chain.i$params, data.frame(.iter=1:iteration, .chain=i)))
    accept.rates <- rbind(accept.rates, chain.i$accept.rate)
    pb$tick()
  }
  
  colnames(accept.rates) <- "Accept rate"
  accept.rates$.chain <- 1:chain.cnt
  
  print("Complete.", quote=F)
  print("Accept rates:", quote=F)
  print(accept.rates)
  
  return(chains)
}

############# Summary MCMC #############

summary.bayes.mcmc <- function(mcmc.out, burn.in=0.5, plot.draw=T, return.df=F) {
  iter <- max(mcmc.out$.iter)
  mcmc.after <- mcmc.out %>% filter(.iter > as.integer(iter*burn.in))
  
  #Summary statistics
  func.list <- list(mean=mean, sd=sd, `ci_5%`=function(x) quantile(x, 0.05), `ci_95%`=function(x) quantile(x, 0.95))
  summary.result <- apply(mcmc.after %>% dplyr::select(-.iter, -.chain), 2, function(column) {
    t(lapply(func.list, function(sm.fn) sm.fn(column)))
  })

  summary.df <- data.frame()
  df.col.names <- names(summary.result)
  for (cn in df.col.names) {
    summary.df <- rbind(summary.df, summary.result[[cn]])
  }
  row.names(summary.df) <- df.col.names
  
  B <- mcmc.after %>% dplyr::select(-.iter, -.chain) %>% summarise_all(var)*10000
  W <- mcmc.after %>% group_by(.chain) %>% summarise_all(var) %>% dplyr::select(-.iter, -.chain) %>% colMeans()
  R.hat <- ((iter-1)+W/B)/iter
  summary.df$R <- t(R.hat %>% round(digits = 2))
  
  print(summary.df)
  
  #Plot parameters
  if (plot.draw) {
    df_long = melt(mcmc.after, id.vars=c('.iter', '.chain'))
    tp <- ggplot(df_long, aes(x = .iter)) +
      geom_line(aes(y=value, color=factor(.chain)))+labs(color="Chain")+
      facet_wrap(~variable, ncol=1, scales = "free")+
      theme_minimal()+theme(axis.title.x=element_blank(), axis.title.y=element_blank())
    
    dp <- ggplot(df_long)+
      geom_density(aes(x=value), fill="skyblue")+
      facet_wrap(~variable, ncol=1, scales = "free")+
      theme_minimal()+theme(axis.title.x=element_blank(), axis.title.y=element_blank())
    print(ggarrange(dp, tp, ncol=2))
  }
  
  if (return.df) {
    return(summary.df)
  }
}
