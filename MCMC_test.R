# MCMC test
#****************************************************************************************************
#  Test MH-MCMC
#****************************************************************************************************

library(MASS)
library(msm)
library(dplyr)
library(ggplot2)
library(coda)

general_van_param_test <- function (sdata, G) {
  
  x <- sdata$h_cm
  y <- sdata$VWC
  
  # G=10000      # should set at least 10000
  burin=2000     # should set to 2000
  
  n <- rep(NA,G)
  theta_s <- rep(NA,G)
  theta_r <- rep(NA,G)
  alpha <- rep(NA,G)
  # sdn <- rep(NA,G)
  # sds <- rep(NA,G)
  
  n[1] <- 1.5
  theta_s[1] <- 0.56
  theta_r[1] <- 0.18
  alpha[1] <- 0.049
  logratio <- 0
  
  sdn = 0.01 #0.01
  sds = 0.05 #0.05
  sdr = 0.025 #0.025
  sda = 0.005 #0.005
  
  for(g in 2:G){
    #--update n upper = 3.5 (uper = 2 for ID:89P00009)
    nnew <-  rtnorm(1, mean=n[g-1], sd=n[g-1]/10, lower=1, upper=10)
    
    meanold <- theta_r[g-1] + (theta_s[g-1]-theta_r[g-1])/ (1+(alpha[g-1]*x)^n[g-1])^(1-1/n[g-1])
    meannew <- theta_r[g-1] + (theta_s[g-1]-theta_r[g-1])/ (1+(alpha[g-1]*x)^nnew)^(1-1/nnew)
      
    logratio <- sum(dnorm(y, mean = meannew, sd = n[g-1]/10, log = TRUE))
    logratio <- logratio-sum(dnorm(y, mean = meanold, sd = n[g-1]/10, log = TRUE))
    logratio <- logratio+dtnorm(n[g-1], mean=nnew, sd=n[g-1]/10, lower=1, upper=10, log = TRUE)
    logratio <- logratio-dtnorm(nnew, mean=n[g-1], sd=n[g-1]/10, lower=1, upper=10, log = TRUE)
    
    if(log(runif(1)) < logratio){
      n[g] <- nnew
    }else {n[g] <- n[g-1]}
    
    # sdn[g] <- ifelse(g < 10, 0.01, sd(n[2:g-1]))
    
    #--update theta_s
    # if (y[1]>0.6) {var_s_upper <- 0.7}
    # else {var_s_upper <- 0.6}
    theta_snew <-  rtnorm(1, mean=theta_s[g-1], sd=theta_s[g-1]/10, lower=0.3003, upper=0.8)
    
    meanold <- theta_r[g-1] + (theta_s[g-1]-theta_r[g-1])/ (1+(alpha[g-1]*x)^n[g-1])^(1-1/n[g-1])
    meannew <- theta_r[g-1] + (theta_snew-theta_r[g-1])/ (1+(alpha[g-1]*x)^n[g-1])^(1-1/n[g-1])
    
    logratio <- sum(dnorm(y, mean = meannew, sd = theta_s[g-1]/10, log = TRUE))
    logratio <- logratio-sum(dnorm(y, mean = meanold, sd = theta_s[g-1]/10, log = TRUE))
    logratio <- logratio+dtnorm(theta_s[g-1], mean=theta_snew, sd=theta_s[g-1]/10, lower=0.3003, upper=0.8, log = TRUE)
    logratio <- logratio-dtnorm(theta_snew, mean=theta_s[g-1], sd=theta_s[g-1]/10, lower=0.3003, upper=0.8, log = TRUE)
    if(log(runif(1))<logratio){
      theta_s[g] <- theta_snew
    }else {theta_s[g] <- theta_s[g-1]}
    
    # sds[g] <- ifelse(g < 10, 0.05, sd(theta_s[2:g-1]))
    
    #--update theta_r
    theta_rnew <- rtnorm(1, mean= theta_r[g-1], sd=theta_r[g-1]/10, lower=0.0003, upper=0.3)
    
    meanold <- theta_r[g-1] + (theta_s[g-1]-theta_r[g-1])/ (1+(alpha[g-1]*x)^n[g-1])^(1-1/n[g-1])
    meannew <- theta_rnew + (theta_s[g-1]-theta_rnew)/ (1+(alpha[g-1]*x)^n[g-1])^(1-1/n[g-1])
    
    logratio <- sum(dnorm(y, mean = meannew, sd = theta_r[g-1]/10, log = TRUE))
    logratio <- logratio-sum(dnorm(y, mean = meanold, sd = theta_r[g-1]/10, log = TRUE))
    logratio <- logratio+dtnorm(theta_r[g-1], mean=theta_rnew, sd=theta_r[g-1]/10, lower=0.0003, upper=0.3, log = TRUE)
    logratio <- logratio-dtnorm(theta_rnew, mean=theta_r[g-1], sd=theta_r[g-1]/10, lower=0.0003, upper=0.3, log = TRUE)
    if(log(runif(1))<logratio){
      theta_r[g] <- theta_rnew
    }else {theta_r[g] <- theta_r[g-1]}
    
    
    #--update alpha
    alphanew <- rtnorm(1, mean=alpha[g-1], sd=alpha[g-1]/10, lower=0.00011, upper=1)
    
    meanold <- theta_r[g-1] + (theta_s[g-1]-theta_r[g-1])/ (1+(alpha[g-1]*x)^n[g-1])^(1-1/n[g-1])
    meannew <- theta_r[g-1] + (theta_s[g-1]-theta_r[g-1])/ (1+(alphanew*x)^n[g-1])^(1-1/n[g-1])
    
    logratio <- sum(dnorm(y, mean = meannew, sd = alpha[g-1]/10, log = TRUE))
    logratio <- logratio-sum(dnorm(y, mean = meanold, sd = alpha[g-1]/10, log = TRUE))
    logratio <- logratio+dtnorm(alpha[g-1], mean=alphanew, sd=alpha[g-1]/10, lower=0.00011, upper=1, log = TRUE)
    logratio <- logratio-dtnorm(alphanew, mean=alpha[g-1], sd=alpha[g-1]/10, lower=0.00011, upper=1, log = TRUE)
    if(log(runif(1))<logratio){
      alpha[g] <- alphanew
    }else {alpha[g] <- alpha[g-1]}
    
    print(g)
  }
  
  # put parameters in data frame
  thin = seq(burin,g,1)
  param <- cbind(n[thin], theta_s[thin], theta_r[thin], alpha[thin])
  colnames(param) <- c("n", "theta_s", "theta_r", "alpha")
  
  return (param)
}

# load data and get VG model parameters
lu_2008_2013_data <- read.csv("data/Lu_2008_2013_Data.csv")
soil_all_param1_test1 <- general_van_param_test (sdata = lu_2008_2013_data %>% filter(soil == "I"), G = 10000)

# trace plot
G = 10000
tibble(x = (2000:G),
       y = soil_all_param1_test1[, 1]) %>% 
  ggplot(aes(x, y)) +
  geom_line() +
  scale_x_continuous(breaks = seq(2000,G,5000))+
  theme(axis.title = element_blank())

tibble(x = (2000:G),
       y = soil_all_param1_test1[, 2]) %>% 
  ggplot(aes(x, y)) +
  geom_line() +
  scale_x_continuous(breaks = seq(2000,G,5000))+
  theme(axis.title = element_blank())

tibble(x = (2000:G),
       y = soil_all_param1_test1[, 3]) %>% 
  ggplot(aes(x, y)) +
  geom_line() +
  scale_x_continuous(breaks = seq(2000,G,5000))+
  theme(axis.title = element_blank())

tibble(x = (2000:G),
       y = soil_all_param1_test1[, 4]) %>% 
  ggplot(aes(x, y)) +
  geom_line() +
  scale_x_continuous(breaks = seq(2000,G,5000))+
  theme(axis.title = element_blank())

soil_all_param1_test1 = as.data.frame(soil_all_param1_test1)


# plot mean vs standard diviation
mov_mn_sdv <- function(sdata, row_n) {
  MN <- rep(NA, nrow(sdata)-5)
  SDV <- rep(NA, nrow(sdata)-5)
  
  for (i in 5:nrow(sdata)) {
    subdata = sdata[1:i,row_n]
    mn = mean(subdata)
    sdv = sd(subdata)
    
    MN[i-5] = mn
    SDV[i-5] = sdv
    
    print(i)
  }
  mov_m_sd <- as.data.frame(cbind(MN, SDV))
  return (mov_m_sd)
}

mov_m_sd <- mov_mn_sdv(soil_all_param1_test1, 1)
mov_m_sd_s <- mov_mn_sdv(soil_all_param1_test1, 2)
mov_m_sd_r <- mov_mn_sdv(soil_all_param1_test1, 3)
mov_m_sd_a <- mov_mn_sdv(soil_all_param1_test1, 4)

mov_m_sd %>% 
  ggplot(aes(x = 1:nrow(mov_m_sd), y = MN)) +
  geom_line() +
  geom_line(aes(x = 1:nrow(mov_m_sd), y = SDV), col = "red") +
  labs(x = expression(Sampling~time),
       y = expression(n))

mov_m_sd_s %>% 
  ggplot(aes(x = 1:nrow(mov_m_sd), y = MN)) +
  geom_line() +
  geom_line(aes(x = 1:nrow(mov_m_sd), y = SDV), col = "red") +
  labs(x = expression(Sampling~time),
       y = expression(theta_s)) 

mov_m_sd_r %>% 
  ggplot(aes(x = 1:nrow(mov_m_sd), y = MN)) +
  geom_line() +
  geom_line(aes(x = 1:nrow(mov_m_sd), y = SDV), col = "red") +
  labs(x = expression(Sampling~time),
       y = expression(theta_r)) 

mov_m_sd_a %>% 
  ggplot(aes(x = 1:nrow(mov_m_sd), y = MN)) +
  geom_line() +
  geom_line(aes(x = 1:nrow(mov_m_sd), y = SDV), col = "red") +
  labs(x = expression(Sampling~time),
       y = expression(alpha)) 


# geweke diagnostic
post=as.matrix(soil_all_param1_test1)
postmcmc=as.mcmc(post)
raftery.diag(postmcmc,q=0.025)
geweke.diag(post)

  