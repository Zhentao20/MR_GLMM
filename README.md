## Overview

This project introduces a class of matrix-response generalized linear mixed models (MR-GLMMs) designed for longitudinal brain network analysis, where the responses are time-indexed square matrices. After applying canonical link transformations, the conditional means are modeled as linear combinations of time-varying covariates with matrix-valued coefficients, along with a matrix-valued intercept that incorporates random effects.

To capture key structural properties, the intercept term is assumed to be low-rank, while the coefficient matrices are assumed to be sparse. Parameter estimation is formulated as a constrained likelihood maximization problem, and an efficient Monte Carlo Expectation–Maximization (MCEM) algorithm is developed to compute the estimators.

This document provides a quick guide to the `mrglmm_id()` and `mrglmm_logit()` functions, specifies the required input data structure for the algorithm to run properly, and explains each component of the function’s output.

## Installation and Setup

Before using the functions in this repository, make sure you have the following R packages installed:

```r
install.packages("Rcpp")
install.packages("RcppArmadillo")
```
Then, clone this repository:
```bash
git clone https://github.com/Zhentao20/MR_GLMM.git
cd MR_GLMM
```

## Notations

The functions `mrglmm_id()` and `mrglmm_logit()` interface with the Rcpp implementations in `EM_grad_linear.cpp` and `EM_grad_logistic.cpp`, corresponding to the identity and logit link functions, respectively.The notations in the code are:

-   $N$: The total number of samples
-   $d$: The number of nodes/regions of the matrix response
-   $T$: The number of measurements for each subject
-   $p$: The number of covariates

## Input Data Structure

The following are the data structure of the functions' input:

-   $Y$: Matrix responses of all subjects across time points $T$. This is a list of length $T$, where each element is an array of dimension $d \times d \times N$.
-   $X$: Covariate information for all subjects. This is a list of length $p$, where each element is a vector of length $N \times t$ representing covariate values for all subjects across $t$ time points.
-   $A$: A subject-centered reorganization of the data. Each entry of $A$ corresponds to a single subject and stores that subject’s sequence of matrix responses across time. For example, in the simulation setting:

``` r
A <- list()
for (i in 1:N) {
  A[[i]] <- list(Y[[1]][,,i], Y[[2]][,,i], Y[[3]][,,i], 
                 Y[[4]][,,i], Y[[5]][,,i])
}
```

where setting $T=5$, $A[[i]]$ contains all the matrices for subject $i$. Each element inside $A[[i]]$ is a $d \times d$ matrix, one per time point. This format makes it easier for the algorithm to iterate over subjects and access their longitudinal matrix responses.

-   $Xt$: We reorganize X into a subject–time–structured list
```math
Xt = \{ Xt_{i,t} : i = 1, …, N ;  t = 1, …, T \},
```
where each element is defined as
```math
Xt_{i,t} = ( X^{(1)}_{(i-1)T+t},  X^{(2)}_{(i-1)T+t}, …, X^{(p)}_{(i-1)T+t} ).
```
Thus:
- $Xt_{i,t}$ corresponds to subject $i$ at time $t$.
- Each $Xt_{i,t}$ is a $p$-dimensional covariate vector at that subject–time combination.
- Collectively, $Xt$ is a list of length $N \times T$, storing covariates in subject–time order.

Taking the simulation setting for example:

``` r
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

where setting $T=5$, $Xt$ is a list of length $N \times T$. Each element $Xt_{i,t}$ corresponds to subject $i$ at time $t$, and is stored as a vector of length $p$ containing all covariates.

-   $t$: The number of measurements for each subject.\
-   $M$: The number of replications, specifying how many Monte Carlo samples are drawn at each iteration in the E-step to approximate the expectation.
-   $b0$: A list of length $p$, where each element is a $d \times d$ zero matrix. This provides empty placeholders for the matrix-valued coefficient matrices.
-   $bt$: Same structure as $b0$, but used as the initial values for the fixed-slope estimates.
-   $sgamma0$: A $d \times d$ zero matrix serving as an empty placeholder for the matrix-valued standard deviation of the random intercept.
-   $sgammat$: Same structure as $sgamma0$, but used as the initial values for the standard deviation of the random intercept.
-   $se0$: A list of length $N$, where each element is the same $d \times d$ matrix filled with initial values. This list represents the initial error (noise) term assigned to each subject.
-   $tol$: Convergence tolerance used as the stopping criterion for the EM algorithm.
-   $tol1$: Convergence tolerance used as the stopping criterion for the backtracking line search.
-   $maxit$: The maximum number of iterations in the EM algorithm.
-   $maxit1$: The maximum number of backtracking line search process.
-   $maxit2$: The maximum number of iterations in the EM algorithm during the symmetrization phase.
-   $lhs$: Initial value of the Q-function used in backtracking line search.
-   $Abar$: The $d \times d$ mean matrix, obtained by averaging all subjects’ matrix responses across time points.
-   $\lambda$: A diagonal identity matrix, used as the initial value of $\lambda$ in the symmetrization phase.
-   $U$: True factor matrix $U$ from the factorization of the fixed intercept (used in simulation experiments).
-   $V$: True factor matrix $V$ from the factorization of the fixed intercept (used in simulation experiments).
-   $Btrue$: True value of the tensor parameters $B$ (used in simulation experiments).
-   $N$: The total number of samples.
-   $d$: The number of nodes/regions in the matrix response.
-   $p$: The number of covariates.
-   $r$: The assumed rank of $\Theta$.
-   $s$: The assumed sparsity level of $B$. For each covariate coefficient estimate, elements smaller than the top $s \times 100\text{\%}$ largest values are truncated to zero.
-   $step1$: Baseline step size in the gradient descent algorithm for estimating $U$ and $V$.
-   $step2$: Baseline step size in the gradient descent algorithm for estimating $B$.
