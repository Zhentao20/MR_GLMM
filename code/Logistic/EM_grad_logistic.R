seed <- as.integer(abs(rnorm(1) * 100000))

#args <- commandArgs(trailingOnly = TRUE)
#r <- as.numeric(args[1])
#s <- as.numeric(args[2])

#-------program body-----  
library(Rcpp)
library(RcppArmadillo)


d <- 10
set.seed(7)
N <- 20
t <- matrix(rep(5,N),nrow=N,ncol = 1)
btrue <- list()

for (i in 1:5){
  btrue [[i]] <- matrix(0, nrow = d, ncol = d)
  
  # Calculate the number of elements that should be 1
  num_elements <- length(btrue [[i]])
  num_to_change <- floor(0.1 * num_elements)
  
  # Randomly select positions to change to 1
  positions <- sample(1:num_elements, size = num_to_change)
  
  # Change those positions to 1
  btrue [[i]] [positions] <- 2
}


s2t_true <- matrix(as.numeric(rep(0.5,d^2)),nrow=d,ncol=d)
s2r_true <- matrix(as.numeric(rep(0.5,d^2)),nrow=d,ncol=d)


vec1 <- matrix(rnorm(d*2),nrow=d, ncol=2)
# Create matrices U and V each with two unique columns.(U and V are the same for simplification)
U <- vec1
V <- vec1

# Compute the product
theta_f = U %*% t(V)

#set.seed(1)
fun0 <- function(x){matrix(rnorm(d^2, mean=x, sd=s2r_true),nrow=d,ncol=d)}
theta_R_mean <- rep(0,d^2)
theta_R <- array(mapply(fun0,theta_R_mean),dim=c(d,d,5*N))


set.seed(seed)
library(MASS)
mean_vec <- c(0,0,0,0) # Means for sex and age

# Covariance matrix
cov_matrix <- matrix(c(1,   0.2, 0.2,  0.2,
                       0.2, 1,   0.2,  0.2,
                       0.2, 0.2, 1,    0.2,
                       0.2, 0.2, 0.2,  1), byrow = TRUE, nrow = 4)
sample <- mvrnorm(N, mu = mean_vec, Sigma = cov_matrix)

x1 <- matrix(rep(sample[, 1],each=5), nrow = 5*N, ncol = 1)
x2 <- matrix(rep(sample[, 2],each=5), nrow = 5*N, ncol = 1)
x3 <- matrix(rep(sample[, 3],each=5), nrow = 5*N, ncol = 1)
x4 <- matrix(rep(sample[, 4],each=5), nrow = 5*N, ncol = 1)
x5 <- matrix(rep(c(-2:2),N), nrow = 5*N, ncol = 1)

#sex <- matrix(rnorm(200, mean=0, sd=1))
#age <- matrix(rnorm(200, mean=0, sd=1))

fun1 <- function(x){x*btrue[[1]]}
bt1 <- array((mapply(fun1,x1)),dim = c(d,d,5*N))
fun2 <- function(x){x*btrue[[2]]}
bt2 <- array((mapply(fun2,x2)),dim = c(d,d,5*N))
fun3 <- function(x){x*btrue[[3]]}
bt3 <- array((mapply(fun3,x3)),dim = c(d,d,5*N))
fun4 <- function(x){x*btrue[[4]]}
bt4 <- array((mapply(fun4,x4)),dim = c(d,d,5*N))
fun5 <- function(x){x*btrue[[5]]}
bt5 <- array((mapply(fun5,x5)),dim = c(d,d,5*N))
#mean
theta_F <- array(theta_f,dim=c(d,d,5*N))

mean_all <- theta_F + theta_R + bt1 + bt2 + bt3 + bt4 + bt5 
Btrue <- list(btrue[[1]],btrue[[2]],btrue[[3]],btrue[[4]],btrue[[5]])
#model <- glm(value ~ b1 + b2 + b3 + b4 + b5, family = binomial(link = "logit"), data = df)
funexp <- function(x){matrix(exp(x)/(1+exp(x)),nrow=d,ncol=d,byrow=F)}
pi <- array(apply(mean_all,3,funexp),dim = c(d,d,5*N))

#fun3 <- function(x){matrix(rnorm(d^2,mean=x,sd=s2t_true),nrow=d,ncol=d,byrow=F)}
fun3 <- function(x){matrix(rbinom(d^2,size=1,prob=x),nrow=d,ncol=d,byrow=F)}

pilist <- vector("list", 5)

for (i in 1:5) {
  pilist[[i]] <- pi[,,seq(i, N*5, by = 5)]
}
#make response symmetric
t1 <- array(apply(pilist[[1]],3,fun3),dim = c(d,d,N))
t2 <- array(apply(pilist[[2]],3,fun3),dim = c(d,d,N))
t3 <- array(apply(pilist[[3]],3,fun3),dim = c(d,d,N))
t4 <- array(apply(pilist[[4]],3,fun3),dim = c(d,d,N))
t5 <- array(apply(pilist[[5]],3,fun3),dim = c(d,d,N))

#Y:1000 x 4 x 4
library(abind)
# Y <- abind(t1,t2,t3,t4)

#list with length 4: each is a array with dimension 30x30x200
Y <- list(t1,t2,t3,t4,t5)
Abar1 <- apply(Y[[1]],c(1,2),mean)
Abar2 <- apply(Y[[2]],c(1,2),mean)
Abar3 <- apply(Y[[3]],c(1,2),mean)
Abar4 <- apply(Y[[4]],c(1,2),mean)
Abar5 <- apply(Y[[5]],c(1,2),mean)

Abar <- as.matrix(Abar1+Abar2+Abar3+Abar4+Abar5)

#list with length 2: each is a vector(matrix form) with length 200
# X <- abind(Sex, Age)
X <- list(x1,x2,x3,x4,x5)

A <- list()
for (i in 1:N){
  A[[i]] = list(Y[[1]][,,i],Y[[2]][,,i],Y[[3]][,,i],Y[[4]][,,i],Y[[5]][,,i])
}

Xt <- list()

#Xt is length with 1000, each set of 5 is corresponded to one subject
for(i in 1:N){
  for(t in 1:5){
  Xt[[(i-1)*5+t]] = list(X[[1]][(i-1)*5+t],X[[2]][(i-1)*5+t],X[[3]][(i-1)*5+t],X[[4]][(i-1)*5+t],X[[5]][(i-1)*5+t])
  }
}

b_sext <- matrix(100*rbinom(n=d^2,size=1,prob = 0.1),nrow=d, byrow = F) #30x30 matrix
bt <- list(b_sext,b_sext,b_sext,b_sext,b_sext)

b_sex0 <- matrix(rep(2,d^2),nrow=d, byrow = F) #30x30 matrix
b0 <- list(b_sex0,b_sex0,b_sex0,b_sex0,b_sex0)
sgammat <- s2r_true
sgamma0 <- matrix(rnorm(n=d^2,mean = 2, sd= 0.5), nrow=d, byrow = F)

set <- s2t_true

num_matrices <- N
matrix_rows <- d
matrix_cols <- d

se0 <- rep(list(s2t_true), N)


tol = 1/10000
tol1 = 1/100

#maxit: number of outer loop1, maxit1: number of inner loop,  maxit2: number of outer loop2
maxit = 20
maxit2 = 5
maxit1 = 2

s = 0.1

r = 2

lambda <- diag(x=1, r,r)  

lhs = 1500

t <- rep(5,N)
M=10


sourceCpp("EM_grad_logistic.cpp")
step1 <- 0.1
step2 <- 0.1
p=5
d=10

#eps, Eps, iter, iter1
op <- mrglmm_logit(Y, X, A, Xt, t, M, b0, sgamma0, sgammat, tol, maxit, lhs, Abar, lambda, bt, tol1, maxit1, maxit2, s, U, V, Btrue,1/8,r,N,p,d,step1,step2)

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

sens_net <- compute_sens_spec(vecA,vecB)[[1]] #1
spec_net <- compute_sens_spec(vecA,vecB)[[2]] #0.944444


TP <- sum(vecA  != 0 & vecB != 0)
FP <- sum(vecA  == 0 & vecB != 0)
TN <- sum(vecA  == 0 & vecB == 0)
FN <- sum(vecA  != 0 & vecB == 0)

F1_net <- 2*TP/(2*TP +FP + FN)
ebic <- op[[16]]
lhs <- op[[11]][length(op[[11]])]
ebica <- op[[18]]













