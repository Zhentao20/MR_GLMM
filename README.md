## Overview

This project introduces a class of matrix-response generalized linear mixed models (MR-GLMMs) designed for longitudinal brain network analysis, where the responses are time-indexed square matrices. After applying canonical link transformations, the conditional means are modeled as linear combinations of time-varying covariates with matrix-valued coefficients, along with a matrix-valued intercept that incorporates random effects.

To capture key structural properties, the intercept term is assumed to be low-rank, while the coefficient matrices are assumed to be sparse. Parameter estimation is formulated as a constrained likelihood maximization problem, and an efficient Monte Carlo Expectation–Maximization (MCEM) algorithm is developed to compute the estimators.

This document provides a quick guide to the mrglmm( ) function, specifies the required input data structure for the algorithm to run properly, and explains each component of the function’s output.


## Notations
The function mrglmm( ) is corresponded to Rcpp function EM_grad_linear.cpp or EM_grad_logistic.cpp targeting identity and logit link. The notations in the code are:

* $N$: The total number of samples  
* $d$: The number of nodes/regions of the matrix response  
* $T$: The number of measurements for each subject  
* $p$: The number of covariates  

## Input Data Structure
The following are the data structure of the input for mrglmm( ):

- \(Y\): Matrix responses of all subjects along time points \(t\). A list with length \(t\). Each is a array with dimension \(d \times d \times N\).
- \(X\): Covariates information of all subjects. A list length \(p\). Each is a vector with length \(N \times t\) representing individual covariate information of all subjects across \(t\). 
- \(A\): To reorganize the data into a subject-centered format, we construct a new list \(A\). Each entry of \(A\) corresponds to a single subject and stores that subject’s sequence of matrix responses across time.Taking the simulation setting for example:
```r
A <- list()
for (i in 1:N) {
  A[[i]] <- list(Y[[1]][,,i], Y[[2]][,,i], Y[[3]][,,i], 
                 Y[[4]][,,i], Y[[5]][,,i])
}
```
where setting \(T=5\), \(A[[i]]\) contains all the matrices for subject \(i\). Each element inside \(A[[i]]\) is a \(d \times d\) matrix, one per time point. This format makes it easier for the algorithm to iterate over subjects and access their longitudinal matrix responses.

- \(Xt\):We reorganize \(X\) into a subject–time–structured list  
\(\;Xt = \{ Xt_{i,t} : i=1,\ldots,N;\; t=1,\ldots,T \}\),  where each element is defined as \(\;Xt_{i,t} = ( X^{(1)}_{(i-1)T+t}, \; X^{(2)}_{(i-1)T+t}, \; \ldots, \; X^{(p)}_{(i-1)T+t} )\).  Thus:  \(Xt_{i,t}\) corresponds to subject \(i\) at time \(t\). Each \(Xt_{i,t}\) is a \(p\)-dimensional covariate vector at that subject–time combination. Collectively, \(Xt\) is a list of length \(N \times T\), storing covariates in subject–time order. Taking the simulation setting for example:
```r
Xt <- list()
for (i in 1:N) {
  for (t in 1:5) {
    Xt[[(i-1)*5+t]] <- list(
      X[[1]][(i-1)*5+t],
      X[[2]][(i-1)*5+t],
      X[[3]][(i-1)*5+t],
      X[[4]][(i-1)*5+t],
      X[[5]][(i-1)*5+t]
    )
  }
}
```
where setting \(T=5\), \(Xt\) is a list of length \(N \times T\). Each element \(Xt_{i,t}\) corresponds to subject \(i\) at time \(t\), and is stored as a vector of length \(p\) containing all covariates.

- \(t\): The number of measurements for each subject.
- \(M\): The number of replications, which specifies how many Monte Carlo samples are drawn at each iteration to estimate the expectation(E-step).
- \(b0\): A list of length \(p\), where each element is a copy of the \(d \times d\) zero matrix. This provides empty placeholders for the matrix-valued coefficients matrix.  
- \(bt\): Similar structure as \(b0\) but used as initial value of the estimation for fixed slopes.
- \(sgamma0\): A zero matrix with dimension \(d \times d\). This provides empty placeholders for the matrix-valued standard deviation of random intercept. 
- \(sgammat\): Similar structure as \(sgamma0\) but used as initial value of the estimation for sd of random intercept. 
- \(se0\): A list of length \(N\), where each element is the same \(d \times d\) matrix filled with initial values. This list represents the initial error (noise) term assigned to each subject.
- \(tol\): EM stopping criterion.
- \(maxit\): The maximum number of iterations in EM algorithm.
- \(maxit2\): The maximum number of iterations in EM algorithm(symmetrization phase).
- \(t\):
- \(t\):
- \(t\):
- \(t\):
- \(t\):
 


```R
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.

