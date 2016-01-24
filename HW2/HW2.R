#7.3a
#Probability of 1 or more 100 year floods in 5 years
1-dbinom(0,5,1/100)
#1 minus the probability of not observing a 100year flood in 5 years

#7.3b
#Probability of 1 or more 100 year floods in 100 years
1-dbinom(0,100,1/100)

#7.3c
1-dbinom(0,1*1000,1/100)
#Because sites are independent and we are only considering 1 year, 
# we have 1 year * 1000 sites opportunities to observe a 100 year flood event

p=c(5,10) #prices
pp=c(0.4,0.6) #probability of prices

q=c(30,55,80,100,120) #discrete quantities q 

pqbar5=c(0,0.15,0.3,0.35,0.2) #conditional probabilities of levels of q given price
pqbar10=c(0.2,0.3,0.4,0.1,0)

#7.4a 
#these give joint probabilities of p and q; from Pxy=Px|y*Py
jointpq=rbind(data.frame(joint=pqbar5*pp[1],price=p[1],qty=q),data.frame(joint=pqbar10*pp[2],price=p[2],qty=q))
most_likely=which(jointpq$joint==max(jointpq$joint))
mlrow=jointpq[most_likely,]
mlrow$price*mlrow$qty #most likely revenue, in dollars

#7.4b mean and variance of future water sales
mean_qty=sum(jointpq$joint*jointpq$qty)
mean_qty
var_qty=sum(((jointpq$qty-mean_sales)^2)*jointpq$joint)
var_qty

#7.4c median value and IQR
combined=jointpq$joint[which(jointpq$price==5)]+jointpq$joint[which(jointpq$price==10)] 
#probabilities assoc with q, ignoring p
cfor=cumsum(combined)
median_qty=q[which(cfor >= 0.5)[1]]
median_qty #median
firstq_qty=q[which(cfor >= 0.25)[1]]
firstq_qty #first quartile
thirdq_qty=q[which(cfor >= 0.75)[1]]
thirdq_qty #third quartile

#7.4d
for(i in p) {
  who=which(jointpq$price==i)
  print(paste("Price: $",i,"/unit has expected revenue ",sum(jointpq$joint[who]*jointpq$price[who]*jointpq$qty[who])))
}
#Price level of $10 will maximize revenue from water sales
