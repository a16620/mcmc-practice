# mcmc-practice
R code for implementing MCMC and other simulations.

### How to use
It is just script, not R package. You can load script by 'source' function.
But it will mess up your environment. I recommend this way.
```R
BCART <- new.env()
source("./Bayes_CART.R", local=BCART)

model0 <- BCART$lapBCART(Y~., "exposure", example.data, method = "zip")
```

### Citations
#### For 'Bayes_GLM.R'
```
Scott A. Baldwin, Michael J. Larson,
An introduction to using Bayesian linear regression with clinical data,
Behaviour Research and Therapy, Volume 98, 2017, Pages 58-75
```

#### For 'Bayes_CART.R'
Normal regression and classification
```
Chipman, H. A., George, E. I., & McCulloch, R. E. (1998).
Bayesian CART Model Search.
Journal of the American Statistical Association, 93(443), 935–948.
```
Poisson, NegBinomial, ZIP
```
Yaojun Zhang, Lanpeng Ji, Georgios Aivaliotis, Charles Taylor,
Bayesian CART models for insurance claims frequency,
Insurance: Mathematics and Economics, Volume 114, 2024, Pages 108-131.
```
