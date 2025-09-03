seed <- as.integer(abs(rnorm(1) * 100000))

#args <- commandArgs(trailingOnly = TRUE)
# r <- as.numeric(args[1])
# s <- as.numeric(args[2])
# filename1 <- paste("simu1",".RData", sep ="")
# filename2 <- paste("simu2",".RData", sep ="")
# filename3 <- paste("simu3",".RData", sep ="")
# filename4 <- paste("simu4",".RData", sep ="")

library(Rcpp)
library(RcppArmadillo)
library(lme4)
library(rpql)
library(Matrix)

r <- 2
s <- 0.1
print(r)
print(s)
d <- 10
set.seed(seed)
N <- 400
t <- matrix(rep(5,N),nrow=N,ncol = 1)
btrue <- list()

for (i in 1:5){
  btrue [[i]] <- matrix(0, nrow = d, ncol = d)
  
  # Calculate the number of elements that should be 1
  num_elements <- length(btrue [[i]])
  num_to_change <- floor(0.1 * num_elements)
  
  # Randomly select positions to change to 2
  positions <- sample(1:num_elements, size = num_to_change)
  
  btrue [[i]] [positions] <- 2
  
}


#generate the error of noise term and of random intercept
s2t_true <- matrix(as.numeric(rep(0.5,d^2)),nrow=d,ncol=d)
s2r_true <- matrix(as.numeric(rep(2,d^2)),nrow=d,ncol=d)


vec1 <- matrix(rnorm(d*r),nrow=d, ncol=r)
# Create matrices U and V each with two unique columns.(U and V are the same for simplification)
U <- vec1
V <- vec1

# Construct the fixed intercept theta
theta_f = U %*% t(V)

fun0 <- function(x){matrix(rnorm(d^2, mean=x, sd=s2r_true),nrow=d,ncol=d)}
theta_R_mean <- rep(0,d^2)
theta_R <- array(mapply(fun0,theta_R_mean),dim=c(d,d,5*N))

library(MASS)
mean_vec <- c(0,0,0,0) # Means for 4 non-time-varying variables

# Covariance matrix
cov_matrix <- matrix(c(1,   0.05, 0.05,  0.05,
                       0.05, 1,   0.05,  0.05,
                       0.05, 0.05, 1,    0.05,
                       0.05, 0.05, 0.05,  1), byrow = TRUE, nrow = 4)
sample <- mvrnorm(N, mu = mean_vec, Sigma = cov_matrix)

# Replicate each column of 'sample' 5 times by row
X <- apply(sample, 2, function(col) rep(col, each = 5))

# X is now a (5N) x 4 matrix; columns correspond to x1, x2, x3, x4
x1 <- X[, 1]
x2 <- X[, 2]
x3 <- X[, 3]
x4 <- X[, 4]
# x5 is set to be a time-varying variable. Create x5 by repeating -2:2 a total of N times. 
x5 <- matrix(rep(-2:2, N), ncol = 1)

# Collect x-vectors in a list
X_list <- list(x1, x2, x3, x4, x5)

bt_list <- Map(function(xvec, bmat) {
  # For each element of xvec, multiply the matrix bmat by that scalar
  # This produces a list of length(xvec) matrices
  mat_list <- lapply(xvec, function(x) bmat * x)
  
  # Combine those matrices along the 3rd dimension into a d x d x length(xvec) array
  simplify2array(mat_list)
  
}, X_list, btrue)

# bt_list is a list of length 5 containing 3D arrays with dimension d * d * (time * sample_size) 
bt1 <- bt_list[[1]]
bt2 <- bt_list[[2]]
bt3 <- bt_list[[3]]
bt4 <- bt_list[[4]]
bt5 <- bt_list[[5]]

#set mean
theta_F <- array(theta_f,dim=c(d,d,5*N))
mean_all <- theta_F + theta_R + bt1 + bt2 + bt3 + bt4 + bt5 
Btrue <- list(btrue[[1]],btrue[[2]],btrue[[3]],btrue[[4]],btrue[[5]])

fun3 <- function(x){matrix(rnorm(d^2,mean=x,sd=s2t_true),nrow=d,ncol=d,byrow=F)}
mean_inter <- array(apply(mean_all,3,fun3),dim = c(d,d,5*N))

t_list <- lapply(1:5, function(i) {
  slice_i <- mean_inter[, , seq(i, 5*N, by = 5)]     # Extract slices i, i+5, i+10, ...
  array(apply(slice_i, 3, fun3), dim = c(d, d, N))   # Apply fun3 along 3rd dimension, reshape to (d x d x N)
})

t1 <- t_list[[1]]
t2 <- t_list[[2]]
t3 <- t_list[[3]]
t4 <- t_list[[4]]
t5 <- t_list[[5]]

library(abind)

#list with length 5: each is a array with dimension d*d*sample_size
Y <- list(t1, t2, t3, t4, t5)

# Apply 'apply(..., c(1,2), mean)' to each 3D array in Y
Abar_list <- lapply(Y, function(x) apply(x, c(1,2), mean))

Abar1 <- Abar_list[[1]]
Abar2 <- Abar_list[[2]]
Abar3 <- Abar_list[[3]]
Abar4 <- Abar_list[[4]]
Abar5 <- Abar_list[[5]]

Abar <- as.matrix(Abar1+Abar2+Abar3+Abar4+Abar5)

X <- list(x1,x2,x3,x4,x5)

A <- list()
for (i in 1:N){
  A[[i]] = list(Y[[1]][,,i],Y[[2]][,,i],Y[[3]][,,i],Y[[4]][,,i],Y[[5]][,,i])
}

Xt <- list()

#Xt is length with sample size, each set including 5 repeated measurements is corresponded to one subject/participant
for(i in 1:N){
  for(t in 1:5){
    Xt[[(i-1)*5+t]] = list(X[[1]][(i-1)*5+t],X[[2]][(i-1)*5+t],X[[3]][(i-1)*5+t],X[[4]][(i-1)*5+t],X[[5]][(i-1)*5+t])
  }
}


#set initial value for estimation
b_ini <- matrix(N * rbinom(d^2, size = 1, prob = 0.1), nrow = d)
bt <- replicate(5, b_ini, simplify = FALSE)

b_ini0 <- matrix(0, nrow = d, ncol = d)
b0 <- replicate(5, b_ini0, simplify = FALSE)

#for simplicity
sgammat <- s2r_true
sgamma0 <- matrix(rnorm(n=d^2,mean = 2, sd= 0.5), nrow=d, byrow = F)
set <- s2t_true
num_matrices <- N
matrix_rows <- d
matrix_cols <- d
se0 <- rep(list(s2t_true), N)

tol = 1/10000
tol1 = 1/100

#number of iterations
maxit = 500
maxit2 = 200
maxit1 = 30
lambda <- diag(x=1, r,r)  
lhs = 1500
t <- rep(5,N)

M=20

#step1: step size for tuning Theta
step1 <- 0.04
#step2: step size for tuning B
step2 <- 0.4
p=5

#eps Eps iter iter1
sourceCpp("EM_grad_linear.cpp")
op <- mrglmm_id(Y, X, A, Xt, t, M, b0, sgamma0, sgammat, se0,  tol,  maxit, lhs, Abar, lambda, bt, tol1, maxit1, maxit2, s, U, V, Btrue,1/8,r,N,p,d,step1,step2)
save(op,file=filename4)

theta2 <- op[[9]] %*% lambda %*% t(op[[9]])
vecA <- as.vector(cbind(Btrue[[1]],Btrue[[2]],Btrue[[3]],Btrue[[4]],Btrue[[5]]))
vecB <- as.vector(cbind(op[[10]][[1]],op[[10]][[2]],op[[10]][[3]],op[[10]][[4]],op[[10]][[5]]))

b_error_net <- sqrt(sum((vecA - vecB)^2)) 
theta_error_net <- sqrt(sum((theta2 - theta_f)^2))

compute_sens_spec <- function(original, estimated) {
  
  # Define the counts
  TP <- sum(original != 0 & estimated != 0)
  FP <- sum(original == 0 & estimated != 0)
  TN <- sum(original == 0 & estimated == 0)
  FN <- sum(original != 0 & estimated == 0)
  
  # Calculate sensitivity and specificity
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  
  return(list(Sensitivity = sensitivity, Specificity = specificity))
}

sens_net <- compute_sens_spec(vecA,vecB)[[1]] 
spec_net <- compute_sens_spec(vecA,vecB)[[2]] 


TP <- sum(vecA  != 0 & vecB != 0)
FP <- sum(vecA  == 0 & vecB != 0)
TN <- sum(vecA  == 0 & vecB == 0)
FN <- sum(vecA  != 0 & vecB == 0)

F1_net <- 2*TP/(2*TP +FP + FN)
