n=100
mus=rep(1,n)
delta=rep(0,n-1)
for(i in 2:n){
  mus[i]=(1+sqrt(1+4*mus[i-1]^2))/2
}


for(i in 1:length(delta)){
  delta[i]=(mus[i]-1)/mus[i+1]
}

plot(mus)
plot(delta)
points(d2,col="red")
d2=delta
for(i in 1:length(delta)){
  d2[i]=i/(i+3)
}
