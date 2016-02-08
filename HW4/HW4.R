foo=read.csv("../HW1/sebou.csv")
xx=sort(foo$maxQ,decreasing=FALSE)
LL3 <- function(X, m, s, t)(1/((X-t)*s*(2*pi)^0.5))*exp(((-(log(X-t)-m)^2)/(2*s^2)))

dlnorm3 = function (x, meanlog = 0, sdlog = 1, threshold = 0) {
 dlnorm(x = x - threshold, meanlog = meanlog, sdlog = sdlog)
}
#3 PARAMETER LOGNORMAL
first=1
last=length(xx)
#stedinger lower quantile bound for threshold
chat=(xx[first]*xx[last]-(median(xx)^2))/(xx[first]+xx[last]-2*median(xx))
#realspace
plot(xx,dlnorm3(xx,meanlog=log(mean(xx)),sdlog=log(sd(xx)),threshold=chat),main="Sebou R.  Lognormal 3-parameter with real and log-space mean/variance",ylim=c(0,1E-3),ylab="Probability density",xlab="Maximum yearly discharge, m^3/s")
#quantiles
qlnorm(p=c(0.01,0.99),meanlog=mean(log(xx)),sdlog=sd(log(xx)))

#logspace
s_logs=sqrt(log(1+(sd(xx)/mean(xx))^2))
m_logs=(log(mean(xx))-0.5*(s_logs^2))
lines(xx,dlnorm3(xx,meanlog=m_logs,sdlog=s_logs,threshold=chat))
legend("topright", c("Real","Log"), pch=c(1,-1),lty=c(0,1), lwd=c(1,1), title = "Legend")
#quantiles
qlnorm(p=c(0.01,0.99),meanlog=m_logs,sdlog=s_logs)

#2-parameter Gamma
alpha_hat=(mean(xx)^2)/(sd(xx)^2)
beta_hat=(mean(xx)/(sd(xx)^2))
plot(xx,dgamma(xx,alpha_hat,beta_hat),ylim=c(0,1E-3),ylab="Probability density",xlab="Maximum yearly discharge, m^3/s",main="2 and 3 Parameter Gamma Distribution",sub="Method of Moments estimators")
qgam2=qgamma(shape=alpha_hat,scale=beta_hat,p=c(0.01,0.99))
data.frame(shape=alpha_hat,scale=beta_hat)
qgamma(p=c(0.01,0.99),alpha_hat,beta_hat)

#3-parameter Gamma
dgamma3 = function (x, aa = 1, bb = 2, tt = 0) {
  (bb^aa)*(x-tt)^(aa-1)*exp(-(x-tt)*bb)/gamma(aa)
}
qgamma3=function(x,mu,sigma,gam) {
  #Uses Wilson-Hilferty transformation
  mu+sigma*((2/gam*(1+(gam*x/6)-((gam^2)/36))^3)-(2/gam))
}
tau_hat=mean(xx)-2*(sd(xx)/skewness(xx))
alpha_hat=4/(skewness(xx)^2)
beta_hat=(2/(sd(xx)*skewness(xx)))
lines(xx,dgamma3(xx,alpha_hat,beta_hat,tau_hat))
qgamma3(qnorm(p=c(0.01,0.99)),mean(xx),sd(xx),skewness(xx))
legend("topright", c("2 Param.","3 Param."), pch=c(1,-1),lty=c(0,1), lwd=c(1,1), title = "Legend")
data.frame(shape=alpha_hat,scale=beta_hat,location=tau_hat)
