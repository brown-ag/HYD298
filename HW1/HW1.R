install.packages("moments")


library(boot)
#7.12 in Loucks & van Beek

#A
#Eq. 7.45 states that the absolute value of the sample skewness 
# is bounded by the square root of the sample size

skew_you=function(x,n) { #Calculates sample estimate of skewness
  xbar=sum(x)/n
  n*sum((x-xbar)^3)/((n-1)*(n-2))
}
coef_skew=function(x,n) { #Calculates sample estimate of coefficient of skewness
  skew_you(x,n)/(sd(x)^3)
}

foo=100
foo=c(100,rep(0,999))

skewed_dat=c(100,0,0,0,0) #Highly (positively) skewed sample; only 1 nonzero value
skewed_dat_neg=c(-10,0,0,0,0) #Highly (negatively) skewed sample; only 1 nonzero value
less_skewed_dat=c(10,10,0,0,0) #Less (positively) skewed sample; 2 nonzero values
random_dat=rnorm(5,5,2) #Random noisy sample

plot(skewed_dat,type="n",ylim=c(-10,10))
lines(skewed_dat,col="blue")
lines(skewed_dat_neg,col="purple")
lines(less_skewed_dat,col="red")
lines(random_dat,col="green")

coef_skew(skewed_dat,n=5)
coef_skew(skewed_dat_neg,n=5)
coef_skew(less_skewed_dat,n=5)
coef_skew(random_dat,n=5)

#Compare output of co_skew for above samples with the theoretical 
# upper bound to conclude the sample skewness coefficient is bounded
#above and below by +/- sqrt(n), resepctively, as demonstrated by Kirby (1974).
sqrt(5)

#B
#Similar logic can be applied to bounds on the coefficient of variation
coef_var = function(x) {
  sd(x)/mean(x)
}
#Use the positively skewed sample from part A
coef_var(skewed_dat)
sqrt(length(skewed_dat)-1)


