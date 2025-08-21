setwd("~/Desktop/final_results/logi_2_0.1/net")
load("theta_error_net.RData")
mean(result.1)
sd(result.1)

load("b_error_net.RData")
mean(result.2)
sd(result.2)
median(result.2)

load("sensitivity.RData")
mean(result.3)
sd(result.3)

load("specificity.RData")
mean(result.4)
sd(result.4)

load("F1_net.RData")
mean(result.5)
sd(result.5)

setwd("~/Desktop/final_results/logi_2_0.1/naive")
load("scad.RData")
#theta
mean(result.1)
sd(result.1)
median(result.1)
#b
mean(result.2)
sd(result.2)
median(result.2)
#sens
mean(result.3)
sd(result.3)
median(result.3)
#spec
mean(result.4)
sd(result.4)
median(result.4)
#F1
mean(result.5)
sd(result.5)
median(result.5)


load("lasso.RData")
#theta
mean(result.6)
sd(result.6)
median(result.6)
#b
mean(result.7)
sd(result.7)
median(result.7)
#sens
mean(result.8)
sd(result.8)
median(result.8)
#spec
mean(result.9)
sd(result.9)
median(result.9)
#F1
mean(result.10)
sd(result.10)
median(result.10)

