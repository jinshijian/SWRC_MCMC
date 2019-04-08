
library(MASS)
library(msm)
 

#****************************************************************************************************
#  calculate NCSS data van parameters
#****************************************************************************************************
# i= 1
van_cal_param <- function (sdata) {
  # x mean swc under certain pressure head, e.g. x <- c(60,100,330,15000)
  # y mean pressure head, e.g. y <- c(0.449, 0.429, 0.423, 0.4146, 0.127)  
  # mcmc_van (x=c(60,100,330,15000), y=c(0.429, 0.423, 0.4146, 0.127)) 
  for (i in 1:2) {
    x <- c(60,100,330,15000)
    y <- c(sdata[i,4],sdata[i,5],sdata[i,6],sdata[i,7])
    varLabsampnum <- paste(sdata[i,1])
    
    # van_params <- mcmc_van (x, y)
    
    G=10000      # should set at least 10000
    burin=2000     # should set to 2000
    
    n <- rep(NA,G)
    theta_s <- rep(NA,G)
    theta_r <- rep(NA,G)
    alpha <- rep(NA,G)
    
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
      nnew <-  rtnorm(1, mean=n[g-1], sd=sdn, lower=1, upper=3.5)
      
      meanold <- theta_r[g-1] + (theta_s[g-1]-theta_r[g-1])/ (1+(alpha[g-1]*x)^n[g-1])^(1-1/n[g-1])
      meannew <- theta_r[g-1] + (theta_s[g-1]-theta_r[g-1])/ (1+(alpha[g-1]*x)^nnew)^(1-1/nnew)
      
      logratio <- sum(dnorm(y, mean = meannew, sd = sdn, log = TRUE))
      logratio <- logratio-sum(dnorm(y, mean = meanold, sd = sdn, log = TRUE))
      logratio <- logratio+dtnorm(n[g-1], mean=nnew, sd=sdn, lower=1, upper=3.5, log = TRUE)
      logratio <- logratio-dtnorm(nnew, mean=n[g-1], sd=sdn, lower=1, upper=3.5, log = TRUE)
      
      if(log(runif(1)) < logratio){
        n[g] <- nnew
      }else {n[g] <- n[g-1]}
      
      #--update theta_s
      if (y[1]>0.6) {var_s_upper <- 0.7}
      else {var_s_upper <- 0.6}
      
      theta_snew <-  rtnorm(1, mean=theta_s[g-1], sd=sds, lower=0.3003, upper=var_s_upper)
      
      meanold <- theta_r[g-1] + (theta_s[g-1]-theta_r[g-1])/ (1+(alpha[g-1]*x)^n[g-1])^(1-1/n[g-1])
      meannew <- theta_r[g-1] + (theta_snew-theta_r[g-1])/ (1+(alpha[g-1]*x)^n[g-1])^(1-1/n[g-1])
      
      logratio <- sum(dnorm(y, mean = meannew, sd = sds, log = TRUE))
      logratio <- logratio-sum(dnorm(y, mean = meanold, sd = sds, log = TRUE))
      logratio <- logratio+dtnorm(theta_s[g-1], mean=theta_snew, sd=sds, lower=0.3003, upper=var_s_upper, log = TRUE)
      logratio <- logratio-dtnorm(theta_snew, mean=theta_s[g-1], sd=sds, lower=0.3003, upper=var_s_upper, log = TRUE)
      if(log(runif(1))<logratio){
        theta_s[g] <- theta_snew
      }else {theta_s[g] <- theta_s[g-1]}
      
      #--update theta_r
      theta_rnew <- rtnorm(1, mean= theta_r[g-1], sd=sdr, lower=0.0003, upper=0.3)
      
      meanold <- theta_r[g-1] + (theta_s[g-1]-theta_r[g-1])/ (1+(alpha[g-1]*x)^n[g-1])^(1-1/n[g-1])
      meannew <- theta_rnew + (theta_s[g-1]-theta_rnew)/ (1+(alpha[g-1]*x)^n[g-1])^(1-1/n[g-1])
      
      logratio <- sum(dnorm(y, mean = meannew, sd = sdr, log = TRUE))
      logratio <- logratio-sum(dnorm(y, mean = meanold, sd = sdr, log = TRUE))
      logratio <- logratio+dtnorm(theta_r[g-1], mean=theta_rnew, sd=sdr, lower=0.0003, upper=0.3, log = TRUE)
      logratio <- logratio-dtnorm(theta_rnew, mean=theta_r[g-1], sd=sdr, lower=0.0003, upper=0.3, log = TRUE)
      if(log(runif(1))<logratio){
        theta_r[g] <- theta_rnew
      }else {theta_r[g] <- theta_r[g-1]}
      
      
      #--update alpha
      alphanew <- rtnorm(1, mean=alpha[g-1], sd=sda, lower=0.00011, upper=0.1)
      
      meanold <- theta_r[g-1] + (theta_s[g-1]-theta_r[g-1])/ (1+(alpha[g-1]*x)^n[g-1])^(1-1/n[g-1])
      meannew <- theta_r[g-1] + (theta_s[g-1]-theta_r[g-1])/ (1+(alphanew*x)^n[g-1])^(1-1/n[g-1])
      
      logratio <- sum(dnorm(y, mean = meannew, sd = sda, log = TRUE))
      logratio <- logratio-sum(dnorm(y, mean = meanold, sd = sda, log = TRUE))
      logratio <- logratio+dtnorm(alpha[g-1], mean=alphanew, sd=sda, lower=0.00011, upper=0.1, log = TRUE)
      logratio <- logratio-dtnorm(alphanew, mean=alpha[g-1], sd=sda, lower=0.00011, upper=0.1, log = TRUE)
      if(log(runif(1))<logratio){
        alpha[g] <- alphanew
      }else {alpha[g] <- alpha[g-1]}
      
      # print(g)
    }
    
    # put parameters in data frame
    thin = seq(burin,g,1)
    
    n_mcmc <- mean(n[thin]) # vs. 1.5293
    theta_s_mcmc <- mean(theta_s[thin]) # vs. 0.5653
    theta_r_mcmc <- mean(theta_r[thin]) # vs. 0.1676
    alpha_mcmc <- mean(alpha[thin]) # vs. 0.0475
    
    # sd(n[thin])
    # sd(theta_s[thin])
    # sd(theta_r[thin])
    # sd(alpha[thin])
    
    # quantile(n,prob=c(0.975,0.025))
    # quantile(theta_s,prob=c(0.975,0.025))
    # quantile(theta_r,prob=c(0.975,0.025))
    # quantile(alpha,prob=c(0.975,0.025))
    
    #****************************************************************************************************
    # plot MCMC performance figures
    
    par( mar=c(2, 0.5, 0.5, 0.5)
         , mai=c(0.4, 0.3, 0.1, 0.1)  # by inches, inner margin
         , omi = c(0.1, 0.1, 0.3, 0.1)  # by inches, outer margin 
         , mgp = c(0, 0.5, 0) # set distance of axis
         , tcl = 0.2
         , cex.axis = 1
         , mfrow=c(5,2))
    
    # plot 2.1 plot Measured points and fit line
    plot ( log(x), y
           , pch = 19
           , col = "black"
           , las = 0
           , cex = 1.5
           , xlab = "" 
           , ylab = "" 
           , lwd = 2
           , xlim = c(0, 10)
           , ylim = c(0, 0.9)  )
    
    h <- seq(0,20000,1)
    theta_mcmc <- theta_r_mcmc + (theta_s_mcmc-theta_r_mcmc)/ (1+(alpha_mcmc*h)^n_mcmc)^(1-1/n_mcmc)
    # theta_mcmc <- mean(theta_r) + (mean(theta_s)-mean(theta_r))/ (1+(0.0475*h)^mean(n))^(1-1/mean(n))
    points(log(h),theta_mcmc, type="l", col="red", lwd=2, lty=2)
    
    mtext(side = 1, text = expression("Measured data and fitted line"), line = 2.5, cex=1, outer = F)
    
    #2.2 plot hist of n parameter
    truehist(n[thin], col='gray'
             ,xlab='', ylab='', main=''
             ,las=0  )
    
    lines(density(n[thin]), col='blue', lwd=2)
    
    box()
    
    mtext(side = 1, text = expression("Histgram of n"), line = 2.5, cex=1, outer = F)
    
    
    #2.3 plot random work trend
    
    plot(n[thin], xlab='', ylab='', type="l",main="")
    # title ("trace of n", adj=0.5, line=-1)
    mtext(side = 1, text = expression("Trace of n"), line = 2.5, cex=1, outer = F)
    
    plot(theta_s[thin],xlab='', ylab='', type="l",main="")
    # title ("trace of theta_s", adj=0.5, line=-1)
    mtext(side = 1, text = expression("Trace of Theta_s"), line = 2.5, cex=1, outer = F)
    
    plot(theta_r[thin],xlab='', ylab='', type="l",main="")
    # title ("trace of theta_r", adj=0.5, line=-1)
    mtext(side = 1, text = expression("Trace of theta_r"), line = 2.5, cex=1, outer = F)
    
    plot(alpha[thin],xlab='', ylab='', type="l",main="")
    # title ("trace of alpha", adj=0.5, line=-1)
    mtext(side = 1, text = expression("Trace of aptha"), line = 2.5, cex=1, outer = F)
    
    
    # 2.4. Auto-Correlation test
    #****************************************************************************************************
    
    nts <-ts(n[thin])
    theta_sts <-ts(theta_s[thin])
    theta_rts <-ts(theta_r[thin])
    alphats <-ts(alpha[thin])
    
    n.ccf <- acf(nts, lag.max = 300, plot = TRUE, main = "",xlab = "",ylab = "")
    mtext(side = 1, text = expression("Lag of n"), line = 2.5, cex=1, outer = F)
    mtext(side = 2, text = expression("ACF"), line = 1.5, cex=1, outer = F)
    
    theta_s.ccf <- acf(theta_sts, lag.max = 300, plot = TRUE, main = "",xlab = "",ylab = "")
    mtext(side = 1, text = expression("Lag of theta_s"), line = 2.5, cex=1, outer = F)
    mtext(side = 2, text = expression("ACF"), line = 1.5, cex=1, outer = F)
    
    theta_r.ccf <- acf(theta_rts, lag.max = 300, plot = TRUE, main = "",xlab = "",ylab = "")
    mtext(side = 1, text = expression("Lag of theta_r"), line = 2.5, cex=1, outer = F)
    mtext(side = 2, text = expression("ACF"), line = 1.5, cex=1, outer = F)
    
    alpha.ccf <- acf(alphats, lag.max = 300, plot = TRUE, main = "",xlab = "",ylab = "")
    mtext(side = 1, text = expression("Lag of alpha"), line = 2.5, cex=1, outer = F)
    mtext(side = 2, text = expression("ACF"), line = 1.5, cex=1, outer = F)
    
    mtext(side = 3, text = paste0(i, "-SampleNum-", varLabsampnum ), line = 0, cex=1.05, outer = T)
    
    sdata[i, "n_mcmc"] <- n_mcmc
    sdata[i, "theta_s_mcmc"] <- theta_s_mcmc
    sdata[i, "theta_r_mcmc"] <- theta_r_mcmc
    sdata[i, "alpha_mcmc"] <- alpha_mcmc
    
    print(paste0('**********', i))
  }
  
  return (sdata)
  
}

#****************************************************************************************************
#  test the van_cal_param function
#****************************************************************************************************

NCSSInv <- read.csv('NCSS_MCMC2_reorder.csv')
head(NCSSInv)
sapply(NCSSInv, class)
length(levels(NCSSInv$labsampnum))
length(unique(NCSSInv$labsampnum))
nrow(NCSSInv)


pdf(file='MCMC_rosseta2_final.pdf', width = 8, height = 12)
results <- van_cal_param (NCSSInv) 
dev.off()

write.csv(results,"NCSS_MCMC2_final.csv", row.names = F)

