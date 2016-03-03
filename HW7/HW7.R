#Load data
sebou=read.csv("../HW1/sebou.csv")
sebou$maxQ=sort(sebou$maxQ)
n=length(sebou$maxQ)
betar=function(x,r) {
  nn=length(x)
  out=0
  for(i in 1:(nn-r)) {
    out=out+(choose(nn-i,r)*x[i]/choose(n-1,r))
  }
  return (1/nn)*out
}
lambda2=2*betar(sebou$maxQ,1)-betar(sebou$maxQ,0)
lambda3=6*betar(sebou$maxQ,2)-6*betar(sebou$maxQ,1)+betar(sebou$maxQ,0)
tau3hat=lambda3/lambda2#L-skewness
Z=tau3hat/sqrt((0.1866/n)+0.8/(n^2))
Z

#The calculated Z value for the L-moment $\tau_3$ test of normality greatly exceeds the critical value for $\alpha=0.05$ therefore reject the null hypothesis of a normally distributed maximum annual flows for Sebou R. dataset.

#7-17g was completed on prior assignment


pp=function(i,nn) { (i-0.44)/(nn+(0.12)) }
peye=pp(1:n,n)
alpha=sqrt(6*var(sebou$maxQ)/pi^2)
xii=mean(sebou$maxQ)-0.5772*alpha
gum_qua=xii-alpha*log(-log(peye))
plot(gum_qua,sebou$maxQ)
abline(0,1)
cor(gum_qua,sebou$maxQ)
#Critical value for r with alpha=0.05 and n=40 is 0.9594. The Gumbel is not a good fit for Sebou R. dataset using the Weibull plotting position.

Z = (tau3hat-0.170)/sqrt(0.2326/n + 0.70/(n^2))
Z
#The calculated Z value for the Chowdhury $\tau_3$ test greatly exceeds the critical Z value for $\alpha=0.05$ therefore reject the null hypothesis of Gumbel distributed maximum annual flows for Sebou R. dataset.

