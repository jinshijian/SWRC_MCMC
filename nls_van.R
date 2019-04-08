
# load data
# NCSSInv <- read.csv('NCSS_MCMC2_reorder.csv')
# head(NCSSInv)
# sapply(NCSSInv, class)
# length(levels(NCSSInv$labsampnum))
# length(unique(NCSSInv$labsampnum))
# nrow(NCSSInv)

# test data
x <- c(60,100,330,15000) 
y <- c(0.429, 0.423, 0.4146, 0.127)
log_x <- log(x)
test_data <- data.frame(y, log_x)

ls_m <- nls( y ~ r + (s-r)/ (1+(alpha*log_x)^n)^(1-1/n)
             , nls.control(maxiter=5000) , data = test_data
             , start = list(r = 0.1, s = 0.5, n = 1.2, alpha = 0.015), trace = TRUE )

