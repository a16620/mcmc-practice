library(mvtnorm)
library(invgamma)
library(ggplot2)

#데이터 생성
n <- 1000
X <- matrix(c(rep(1,n),sort(rnorm(n, sd=10))), ncol=2)
X <- cbind(X, X[,2]**2, X[,2]**3)

#Truncated polynomial이 적용된 열을 반환한다
make_truncated <- function(X, beg, deg=3) {
  ((X-beg)*as.integer(X>=beg))**deg
}

#True 모델 설정
true.sigma2 <- rchisq(1, 1)
knot_true <- sort(sample(seq(X[100,2], X[900,2], length.out=100), 2))
X_true <- cbind(X, make_truncated(X[,2], knot_true[1]), make_truncated(X[,2], knot_true[2]))
true.beta <- c(rnorm(2, sd=1),rnorm(4, sd=0.01))

############# 주의점: Y범위가 커지면 경우에 따라 MCMC에서의 베타 분산이 커진다.
Y <- rnorm(n, X_true%*%true.beta, sd=sqrt(true.sigma2))
plot(X[,2], Y)

#하나의 Matrix에 다른 길이의 파라미터 벡터를 담기 위한 길이 조정 함수
padding.vector <- function(val, d, default=NA) {
  a <- rep(default, d)
  a[1:length(val)] <- val
  return(a)
}

#베타 모델 사후 확률
lik.beta <- function(X_mat, beta, sigma2) {
  sum(dnorm(Y-X_mat%*%beta, sd=sqrt(sigma2), log = T))-log(sigma2)
}

#최대 knot 수 제한
max_K <- 5

#RJ MCMC
beta.chain <- matrix(NA, ncol=6+max_K, nrow=1000)
colnames(beta.chain) <- c("K", "sigma", paste0("beta",0:3), paste0("beta(knot)", 1:max_K))
knot.chain <- matrix(NA, ncol=max_K, nrow=1000)

X_knot <- NULL #현재 추가된 knot의 위치를 담을 벡터
cand_X <- X[,2]#knot추가시 가능한 knot의 후보. (중복 방지를 위해 관리한다)
X_mcmc <- X    #사후 확률 계산에 이용될 matrix. Trucated된 열이 추가된 상태. 계산 편의을 위해 유지한다

#초기 값
beta <- rnorm(4) #k=0
sigma2 <- rchisq(1,1)
k <- 0 #knot이 없는 상태에서 시작
for (iter in 1:1000) {
  next.k <- k+sample(c(-1,0,1), 1, prob = c(0.33,0.33,0.33)) #knot의 개수를 결정
  if ((next.k-k) != 0 && next.k >= 0 && next.k <= max_K) { 
    if (next.k>k) { #knot 증가
      beta.new <- rnorm(1)
      beta.next <- c(beta, beta.new) #Cubic spline은 값 연속+1차+2차 미분이 연속이므로 따로 보정하지 않아도 됨
      
      next.X.idx <- sample(1:length(cand_X), 1) #knot 후보에서 하나 선택. (n-k)개 중에 하나.
      X_knot.next <- append(X_knot, X_mcmc[next.X.idx,2])
      add.col <- make_truncated(X_mcmc[,2], X_mcmc[next.X.idx,2])
      next.X <- cbind(X_mcmc, add.col)
      
      beta.next[1] <- beta[1]-mean(add.col)*beta.new
      prob.mod <- log(k+1)-log(n-k)-dnorm(beta.new, log=T)
      cand.next <- cand_X[-next.X.idx]
    } else {
      rem.knot.X.idx <- sample(1:k, 1) #현재 knot중에 하나를 제거한다. k개 중에 하나.
      rem.knot.X <- X_knot[rem.knot.X.idx]
      X_knot.next <- X_knot[-rem.knot.X.idx]
      
      rem.col.idx <- 3+rem.knot.X.idx
      rem.col <- X_mcmc[,rem.col.idx]
      next.X <- X_mcmc[,-rem.col.idx]
      
      cand.next <- c(cand.next,X_mcmc[which.min(rem.col == 0)-1,2])
      
      beta.next <- beta[-rem.col.idx]
      beta.next[1] <- beta[1] + mean(rem.col)*beta[rem.col.idx]
      prob.mod <- log(n-k+1)-log(k)+dnorm(beta[length(beta)], log=T)
    }
    
    log.prob <- lik.beta(next.X, beta.next, sigma2)-lik.beta(X_mcmc, beta, sigma2)+prob.mod
    if (log.prob >= log(runif(1))) { #MH sampling
      beta <- beta.next
      X_mcmc <- next.X
      k <- next.k
      cand_X <- cand.next
      X_knot <- X_knot.next
    }
  }
  
  #Gibbs
  Vb <- solve(t(X_mcmc)%*%X_mcmc)
  estB <- Vb%*%t(X_mcmc)%*%Y
  beta <- as.numeric(rmvnorm(1, estB, Vb*sigma2))
  
  sigma2 <- as.numeric(rinvgamma(1, 50, crossprod(Y-X_mcmc%*%beta)/2))
  beta.chain[iter,] <- padding.vector(c(k, sigma2, beta), 6+max_K)
  if (!is.null(X_knot))
    knot.chain[iter,] <- padding.vector(X_knot, max_K)
}

#knot의 개수나 위치가 다르면 모수끼리 비교할 수 없음
#따라서 예측값에 대해 MSE를 비교한다.
mean((Y-X_true%*%true.beta)**2)
mean((Y-X_mcmc%*%beta)**2) #마지막 iterarion이 수렴된 상태임을 가정한다.

#MCMC 체인 시각화
pldf <- data.frame(iter=1:1000, beta.chain) %>% select(-K, -sigma) %>% melt(id="iter", na.rm = T)
ggplot(pldf) + geom_line(aes(x=iter, y=value, color=variable))

#save(beta.chain, knot.chain, file = "./rj2.Rdata")

#수렴된 부분(같은 knot가 유지되는 상태)만 추출한다
conv_chain <- tail(beta.chain, 100)[,3:(6+1)] #자르는 인덱스는 직접 조절해야 함.
conv_knot <- X_knot

pred_X <- seq(min(X[,2]), max(X[,2]), length.out=100)
pred_X_df <- matrix(c(rep(1,100),pred_X, pred_X**2, pred_X**3), ncol=4)
for (xk in 1:length(conv_knot)) {
  pred_X_df <- cbind(pred_X_df, make_truncated(pred_X_df[,2], conv_knot[xk]))
}

#단순히 체인의 모든 행에 y 예측값을 구한다
CI_tmp <- apply(conv_chain, 1, function(beta) {pred_X_df%*%beta})
CI <- t(apply(CI_tmp, 1, function(y_preds) { #95% 분위수 신뢰구간+평균
  qunt <- quantile(y_preds, c(.025, .975))
  return(c(qunt[1], mean(y_preds), qunt[2]))
}))
#ggplot을 위한 Data frame
predplotdf <- data.frame(Xt=X[,2], Yt=Y, Xp=pred_X, Yp=CI[,2], Yp_l=CI[,1], Yp_u=CI[,3])
#실제 데이터와 예측 값을 그린다. 빨강=예측, 검정=관측.
ggplot(predplotdf)+geom_line(aes(x=Xt,y=Yt))+geom_line(aes(x=Xp,y=Yp), color='red')+
  geom_ribbon(aes(x=Xp, ymin = Yp_l, ymax=Yp_u), fill='red', alpha=0.3)+
  theme_classic()


