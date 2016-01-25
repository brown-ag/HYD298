foo=read.csv("../HW1/sebou.csv")
#Fit a normal distribution
#normalize data with first and second moment'
xbar_norm=mean(foo$maxQ)
ssd_norm=sd(foo$maxQ)
z.norm<-(foo$maxQ-xbar_norm)/ssd_norm
#plot normalized data versus normal reference line
qqnorm(z.norm)
abline(0,1)
#Statistical test to confirm the visually observed misfit
shapiro.test(foo$maxQ)
#Mean:
xbar_norm
#Variance:
ssd_norm^2
#Quantiles 1% and 99%
quantile(foo$maxQ,probs=c(0.01,0.99))

#Lognormal
lmaxq=log(foo$maxQ) #log-transform data vector
#Maximum Likelihood Estimators
mlemuhat=sum(lmaxq)/length(lmaxq)
mlemuhat
mlevarhat=sum((lmaxq-mlemuhat)^2)/length(lmaxq)
mlevarhat

#Method of Moments
#Variance:
lvarhat=log(1+ssd_norm^2/xbar_norm^2)
lvarhat
#Mean:
lmuhat=log(xbar_norm)-0.5*lvarhat
lmuhat

#Quantiles 1% and 99
quantile(lmaxq,probs=c(0.01,0.99))


