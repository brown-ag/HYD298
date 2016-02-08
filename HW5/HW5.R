foo=read.csv("sebou.csv")

xx=sort(foo$maxQ,decreasing=FALSE)
plot(density(xx),xlab="Discharge, cfs",ylab="Probability density",main="Empirical and modeled (Gumbel)\nprobability densities of Sebou R.\nmaximum annual flows",lty=1,lwd=2)
#Gumbel
#MM
alpha=sqrt(6*var(xx)/pi^2)
xii=mean(xx)-0.5772*alpha
gum_pdf=(1/alpha)*exp(-((xx-xii)/alpha)-exp(-(xx-xii)/alpha))
lines(xx,gum_pdf,col="BLUE",lty=2,lwd=2)
gum_qua=xii-alpha*log(-log(c(0.01,0.99)))
#PARAMETERS
xii
alpha
#QUANTILES
gum_qua

#L-Moments
pwm_gum = function(x,r) {
  #calculates probability weightetd moments of order r
  # used for calculating betas and associated L-moments
  cs=0
  n=length(x)
  for(i in seq(1,n-r)) {
    cs=cs+(choose(n-i,r)*x[i])/choose(n-1,r)
  }
  return(cs/n)
}

lambda1=mean(xx)
beta1=pwm_gum(sort(xx,decreasing=TRUE),1)
lambda2=2*beta1-lambda1
alpha=lambda2/log(2)
xii=mean(xx)-0.5772*alpha 
gum_pdf=(1/alpha)*exp(-((xx-xii)/alpha)-exp(-(xx-xii)/alpha))
gum_qua=xii-alpha*log(-log(c(0.01,0.99)))
lines(xx,gum_pdf,col="RED",lty=2,lwd=2)
#PARAMETERS
xii
alpha
#QUANTILES
gum_qua

#MLE
library(fitdistrplus)
dgumbel <- function(x, ps, al) {  exp((ps - x)/al - exp((ps - x)/al))/al }
pgumbel <- function(q, ps, al) {  exp(-exp(-((q - ps)/al))) }
qgumbel <- function(p, ps, al) {  ps-al*log(-log(p)) }
gumbel.fit <- fitdist(xx, "gumbel", start=list(ps=mean(xx), al=sd(xx)), method="mle")
xii=gumbel.fit$estimate['ps'][[1]]
alpha=gumbel.fit$estimate['al'][[1]]
gum_pdf=(1/alpha)*exp(-((xx-xii)/alpha)-exp(-(xx-xii)/alpha))
lines(xx,gum_pdf,col="GREEN",lty=2,lwd=2)
legend("topright", c("Empirical","MM","LM","MLE"), col=c('BLACK',"BLUE","RED","GREEN"), lty=c(1,2,2,2), lwd=2, title = "Legend")
#PARAMETERS:
xii
alpha
#QUANTILES:
gum_qua
quantile(gumbel.fit,probs=c(0.01,0.99))
plot(gumbel.fit)

#Q2: GEV
dgev=function (x, xi = 1, mu = 0, sigma = 1) {   
#tmp <- (1 + (xi * (x - mu))/sigma)   
#(as.numeric(tmp > 0) * (tmp^(-1/xi - 1) * exp(-tmp^(-1/xi))))/sigma 
if(xi==0) {
    ttt=(exp(-(x-mu)/sigma))
} else {
    ttt=(1+((x-mu)/sigma)*xi)^(-1/xi)
}
  print(ttt)
  1/sigma*(ttt^(xi+1))*exp(-ttt)  
}
#gev p fn
pgev=function (q, xi = 1, mu = 0, sigma = 1) {   exp(-(1 + (xi * (q - mu))/sigma)^(-1/xi)) }
#gev quantile fn
qgev=function (p, xi = 1, mu = 0, sigma = 1) {  mu + (sigma/xi) * ((-logb(p))^(-xi) - 1) }

plot(density(xx),xlab="Discharge, cfs",ylab="Probability density",
     main="Empirical and modeled (GEV)\nprobability densities of Sebou R.\n
     maximum annual flows",lty=1,lwd=2,ylim=c(0,2e-3))

##### L-Moments
pwm_gev=function(x,r) {
  rr=seq(r+1,length(x))
  sum(choose(rr-1,r)*x[rr]/choose(length(x),r+1))/(r+1)
}
lambda1=pwm_gev(xx,0) #=mean(xx)
beta1=pwm_gev(xx,1)
beta2=pwm_gev(xx,2)
lambda2=2*beta1-lambda1
lambda3=6*beta2-6*beta1+lambda1
tau3=lambda3/lambda2
cc=2/(tau3+3)-log(2)/log(3)
#calculate kappa first
cc=lambda2/(3*beta2-lambda1)-(log(2)/log(3))
kappa=7.8590*cc+2.9554*(cc^2)
alpha=kappa*lambda2/(gamma(1+kappa)*(1-2^-kappa))
xii=lambda1+(alpha/kappa)*(gamma(1+kappa)-1)
gev_cdf=exp(-(1-(kappa*(xx-xii)/alpha))^(1/kappa))
gev_qua=xii+(alpha/kappa)*(1-(-log(c(0.01,0.99)))^kappa)
lines(xx,dgev(xx,xi=kappa,mu=xii,sigma=alpha),col="BLUE",lty=2,lwd=3)
#Theoretical Upper bound
xii-(alpha/kappa) #note that density function is undefined above this value

#PARAMETERS
xii
alpha
kappa

#QUANTILES
gev_qua

##### Maximum Likelihood
#gev density fn

gev.fit <- fitdist(xx, "gev", start=list(mu=0, sigma=1, xi=1), method="mle")
xii=gev.fit$estimate['mu'][[1]]
alpha=gev.fit$estimate['sigma'][[1]]
kappa=gev.fit$estimate['xi'][[1]]
lines(xx,dgev(xx,xi=kappa,mu=xii,sigma=alpha),col="GREEN",lty=2,lwd=2)
legend("topright", c("Empirical","LM","MLE"), col=c('BLACK',"BLUE","GREEN"), lty=c(1,2,2), lwd=c(1,2,2), title = "Legend")

#PARAMETERS:
xii
alpha
kappa

#QUANTILES:
quantile(gev.fit,probs=c(0.01,0.99))
plot(gev.fit)

#GMLE
priork=function(kappa) {
  p=6
  q=9
  gamma(p)*gamma(q)*((0.5+kappa)^(p-1))*((0.5-kappa)^(q-1))/gamma(p+q)
}
gmlf=function(p,xx) {
  xii=p[1]
  alpha=p[2]
  kappa=p[3]
  yy=(1-(kappa/alpha)*(xx-xii))
  -length(xx)*log(alpha)+sum((1/kappa-1)*log(yy)-yy^(1/kappa))+log(priork(kappa))
}
minll=nlm(f=gmlf,p=c(xii,alpha,kappa),xx)
lines(xx,dgev(xx,mu=minll$estimate[1],sigma=minll$estimate[2],xi=minll$estimate[3]))
