#Define functions to calculate Weibull, Hazen and Blom positions and the variance of the ith largest
ppwei=function(i,n)  { i/(n+1) }
pphaz=function(i,n)  { (i-(0.5))/(n) }
ppblom=function(i,n) { (i-(3/8))/(n+(1/4)) }
varpi=function(i,n)  { (i*(n-i-1))/(((n+1)^2)*(n+2)) }

#Generate a random lognormal with n=100 observations logmean=3 logsd=2; sort it from largest->smallest
lll=rlnorm(n=100,meanlog=3,sdlog=3)
lll=sort(lll,decreasing=TRUE)
nn=length(lll)
wei=ppwei(seq(1,length(lll)),length(lll))
haz=pphaz(seq(1,length(lll)),length(lll))
blom=ppblom(seq(1,length(lll)),length(lll))
sddd=(varpi(seq(1,length(lll)),length(lll)))

qiwei=sort(qlnorm(ppwei(1:nn,nn),meanlog=mean(log(lll)),sdlog=sd(log(lll))),decreasing=TRUE)
qihaz=sort(qlnorm(pphaz(1:nn,nn),meanlog=mean(log(lll)),sdlog=sd(log(lll))),decreasing=TRUE)
qiblom=sort(qlnorm(ppblom(1:nn,nn),meanlog=mean(log(lll)),sdlog=sd(log(lll))),decreasing=TRUE)
par(mfrow=c(1,3))
plot(qiwei,lll,xlab="Expected quantile - Weibull plotting position")
abline(0,1)
plot(qihaz,lll,xlab="Expected quantile - Hazen plotting position")
abline(0,1)
plot(qiblom,lll,xlab="Expected quantile - Blom plotting position")
abline(0,1)

mqwei=max(qiwei)
mqhaz=max(qihaz)
mqblo=max(qiblom)
#Comparison of the expected value for the largest observation using the three plotting positions and the actual largest extreme value. In general, these plotting positions fail to predict the magnitude of the 100 year flood (largest of 100 values in dataset). The relative magnitude of the offset of the different plotting positions 
data.frame(Weibull=mqwei,Hazen=mqhaz,Blom=mqblo,Actual=max(lll))

#Where n=i get a negative pi value
#Where i=n-1 get a pi value of zero
#In general, plotting positions do not differ greatly from one another. However, they are crude estimates of the exceedence probabilities of extreme values (e.g. smallest or largest obs)
#The actual exceedance probability for the largest observation lies between
data.frame(min=0.29/100,max=1.38/(100+2))
#Note that these are close to the estimates from the three plotting positions

sebou=read.csv("../HW1/sebou.csv")
sebou$maxQ=sort(sebou$maxQ)
ks.test(x=sebou$maxQ,"pnorm")
n=length(sebou$maxQ)
iii=2:(n-1)
#Komolgorov Smirnov test
#Reject null hypothesis that Sebou River maximum annual flows are normally distributed

#7.17c - Normal
qfit=qnorm(ppblom(1:n,n),mean=mean(sebou$maxQ),sd=sd(sebou$maxQ))
plot(qfit,sebou$maxQ)
abline(0,1)

#7.17d - Lognormal
qlfit=qlnorm(ppwei(1:n,n),meanlog=mean(log(sebou$maxQ)),sdlog=sd(log(sebou$maxQ)))
plot(qlfit,sebou$maxQ)
abline(0,1)

#7.17g
par(mfrow=c(1,2))
pi1=1-((0.5)^(1/n))
pi2n1=(iii-0.3175)/(n+0.365) #Filliben 1975 plotting position
piin=(0.5)^(1/n)
pii=c(pi1,pi2n1,piin)
mii=qnorm(pii,mean=mean(sebou$maxQ),sd=sd(sebou$maxQ))
miiln=qlnorm(pii,meanlog=mean(log(sebou$maxQ)),sdlog=sd(log(sebou$maxQ))) #test
cor1=cor(sebou$maxQ,mii)
cor2=cor(sebou$maxQ,miiln)
plot(mii,sebou$maxQ,sub=paste("Normal r =",signif(cor1,digits=4)))
abline(0,1)
plot(miiln,sebou$maxQ,sub=paste("Lognormal r =",signif(cor2,digits=4)))
abline(0,1)

#After consulting the full table of percent points from Filliben (1975) for n=41, the value of r at the 5% level is 0.972. Niether of the correlation coefficients exceed this value. The lognormal is likely adequate for the mid-range of flows, but becomes inaccurate at very low and very high flows. The normal is generally not a good fit at all, especially due to the issues with symmetry and necessarily positive values for discharge.
