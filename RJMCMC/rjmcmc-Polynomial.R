library(mvtnorm)
library(invgamma)

library(ggplot2)
library(reshape2)
library(dplyr)

#True 모델 설정
true.beta <- rnorm(4, sd=9)
true.sigma2 <- rchisq(1, 1)

#데이터 생성
X_raw <- rnorm(100)
X_max_k <- length(true.beta)+4

X <- sapply(0:X_max_k, function(p) X_raw**p)
Y <- rnorm(100, X[,1:length(true.beta), drop=F]%*%true.beta, sd=sqrt(true.sigma2))

X.mean <- colMeans(X)
Y.mean <- mean(Y)

#베타에 대한 사후 확률: 가변 길이 지원
lik.beta <- function(beta, sigma2) {
  sum(dnorm(Y-X[,1:length(beta), drop=F]%*%beta, sd=sqrt(sigma2), log = T))-log(sigma2)
}

#K(차원 수)에 대한 분포: Uniform하게 가정. 높을수록 적은 확률을 주려했지만 대신 K에 최대값 제한을 두었다
prob.k <- function(k) {
  if (k > X_max_k || k <= 0)
    return(-Inf)
  return(-log(X_max_k)) #uniform
}

#하나의 Matrix에 파라미터를 담기 위한 조정 함수
padding.vector <- function(val, d, default=NA) {
  a <- rep(default, d)
  a[1:length(val)] <- val
  return(a)
}

#MCMC 시작
beta.chain <- matrix(NA, ncol=2+X_max_k, nrow=1000)
colnames(beta.chain) <- c("K", "sigma", paste0("beta",1:X_max_k-1))

beta <- rnorm(1) #k=1
sigma2 <- rchisq(1,1)
k <- length(beta)
for (iter in 1:1000) {
  next.k <- k+sample(c(-1,0,1), 1, prob = c(0.33,0.33,0.33)) #K를 유지할지 변경할 지 결정
  if ((next.k-k) != 0 && !is.infinite(prob.k(next.k))) {
    if (next.k>k) {
      beta.new <- rnorm(1) #보조변수
      beta.next <- c(beta, beta.new)
      beta.next[1] <- beta[1]-X.mean[next.k]*beta.new #현재 수렴된 상태에서 Y값이 크게 달라지면 수락 확률이 너무 작아지기 때문에 조정한다
      prob.mod <- -dnorm(beta.new, log=T) #보조 변수에 대한 density를 추가로 곱해준다. 여기선 분모가 더 큰 차원이지 분자로 들어가도록(-)곱
    } else {
      beta.next <- beta[1:(length(beta)-1)] #가장 높은 차수의 베타를 제거
      beta.next[1] <- beta[1] + X.mean[k]*beta[length(beta)] #마찬가지로 보정 작업. 보정 작업은 자코비안을 통해 확률에도 반영됨
      prob.mod <- dnorm(beta[length(beta)], log=T)
    }
    
    log.prob <- lik.beta(beta.next, sigma2)-lik.beta(beta, sigma2)+prob.mod+prob.k(next.k)-prob.k(k) #MH sampling
    if (log.prob >= log(runif(1))) {
      beta <- beta.next
      k <- next.k
    }
  }
  
  #Gibbs
  Xt <- X[,1:length(beta), drop=F]
  Vb <- solve(t(Xt)%*%Xt)
  estB <- Vb%*%t(Xt)%*%Y
  beta <- as.numeric(rmvnorm(1, estB, Vb*sigma2))
  
  sigma2 <- as.numeric(rinvgamma(1, 50, crossprod(Y-X[,1:length(beta), drop=F]%*%beta)/2))
  beta.chain[iter,] <- padding.vector(c(k, sigma2, beta), 2+X_max_k)
}

#체인을 시각화
pldf <- data.frame(iter=1:1000, beta.chain) %>% select(-K, -sigma) %>% melt(id="iter", na.rm = T)
ggplot(pldf) + geom_line(aes(x=iter, y=value, color=variable))

#모수끼리 정량적 비교
abs(beta-true.beta)
#MSE 비교
mean((Y-X[,1:length(true.beta), drop=F]%*%true.beta)**2)
mean((Y-X[,1:length(beta), drop=F]%*%beta)**2)
