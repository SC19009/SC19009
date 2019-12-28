## ----eval=FALSE---------------------------------------------------------------
#  BCa<-function(x,alpha,B){
#    n<-length(x)
#    bar<-mean(x)
#  
#    bar.boot<-numeric(B)
#    for (i in 1:B) {
#      j<-sample(1:n,n,replace=TRUE)
#      bar.boot[i]<-mean(x[j])
#    }
#  
#    bar.jack<-(sum(x)-x)/(n-1)
#    bar0.jack<-mean(bar.jack)
#  
#    z<-qnorm((1-alpha)/2)
#    z0<-qnorm(mean(bar.boot<bar))
#    alpha0<-sum((bar0.jack-bar.jack)^3)/(6*sum(abs(bar0.jack-bar.jack)^3))
#  
#    alpha1<-pnorm(z0+(z0+z)/(1-alpha0*(z0+z)))
#    alpha2<-pnorm(z0+(z0-z)/(1-alpha0*(z0-z)))
#  
#    lb.BCa<-quantile(bar.boot,alpha1,names=FALSE)
#    ub.BCa<-quantile(bar.boot,alpha2,names=FALSE)
#  
#    out<-c(lb.BCa,ub.BCa)
#    names(out)<-c("lb","ub")
#    return(out)
#  }

## -----------------------------------------------------------------------------
library(SC19009)
x<-rbeta(200,4,4)
BCa(x,0.95,500)

## ----eval=FALSE---------------------------------------------------------------
#  NumericVector jack(NumericVector x){
#    int n=sizeof(x); double s=sum(x); double bar0=s/n;
#    double bar[n]; double y[n]; double sum1=0; double sum2=0;
#  
#    for (int i=0;i<=n-1;i++){
#      bar[i]=(s-x[i])/(n-1);
#      sum1=sum1+bar[i];
#    }
#  
#    double bar1=sum1/n;
#    double jack=n*bar0-(n-1)*bar1;
#  
#    for (int i=0;i<=n-1;i++){
#      double z=bar[i]-bar1;
#      y[i]=z*z;
#      sum2=sum2+y[i];
#    }
#  
#    double var=(n-1)*sum2/n;
#  
#    NumericVector out= NumericVector::create(jack,var);
#    return (out);
#  }

## -----------------------------------------------------------------------------
y<-rnorm(500,4,4)
jack(y)

## -----------------------------------------------------------------------------
x<-rnorm(100,mean=0,sd=1)
y<-rnorm(100,mean=1,sd=2)
plot(c(x,y))
hist(c(x,y))
barplot(c(x,y))

## -----------------------------------------------------------------------------
name<-c("Tom","Tony","Tomy")
age<-c(20,30,40)
table(name,age)

## -----------------------------------------------------------------------------
x<-rpois(100,3)
ftable(x)

## -----------------------------------------------------------------------------
# Define the function to generate samples with varing sigma^2
Inv.sample<-function(n,sigma.square) {
  u<-runif(n)
  y<-sqrt(-2*sigma.square*log(1-u))
  return(y)
}

# Generate samples
x<-matrix(nrow=2000,ncol=6)
sigma.square<-c(0.01,0.1,0.5,1,2,5)  ## Choose different sigma^2
for (i in 1:6) x[,i]<-Inv.sample(2000,sigma.square[i])

# Visualization of samples and theoretical density
for (i in 1:6) {
  hist(x[,i],prob=TRUE,breaks=100,xlab="x",main=paste("sigma.square=",sigma.square[i])) ## histogram
  y<-seq(0,10,0.001)
  lines(y,y*exp(-y^2/(2*sigma.square[i]))/sigma.square[i]) ## theoretical density curve
}

## -----------------------------------------------------------------------------
# Define the function to generate samples with varying p1
Trans.sample<-function(p1) {
  x1<-rnorm(1000,mean=0,sd=1)
  x2<-rnorm(1000,mean=3,sd=1)
  r<-sample(c(1,0),1000,replace=TRUE,prob=c(p1,1-p1))
  x<-r*x1+(1-r)*x2
  return(x)
}

# Sampling with the given probabilty : p1=0.75
x0<-Trans.sample(0.75)
hist(x0,breaks=100,prob=TRUE,main="p1=0.75")

# Sampling with different p1
p1<-seq(0.1,0.9,0.1)
x<-matrix(nrow=1000,ncol=9)
y<-seq(-3,6,0.001)
for (i in 1:9) {
  x[,i]<-Trans.sample(p1[i])
  hist(x[,i],breaks=100,prob=TRUE,xlab="x",main=paste("p1=",p1[i]))
  lines(y,p1[i]*dnorm(y,mean=0,sd=1)+(1-p1[i])*dnorm(y,mean=3,sd=1)) ## superimposing density
}

## -----------------------------------------------------------------------------
# Define the function to generate Wishart samples with different n & d
Wis.sample<-function(d,n,sigma) {
  T<-matrix(0,nrow=d,ncol=d)
  if (n>(d+1) & (d+1)>1) {
    ## Common case when d>1
    if (d>1) {
      ### Generate lower triangular entries
      for (i in 2:d) {
        for (j in 1:i-1) T[i,j]<-rnorm(1,mean=0,sd=1)
      }

      ### Generate entries on the diagonal
      for (i in 1:d) {
        x<-rnorm(n-i+1,mean=0,sd=1)
        T[i,i]<-sqrt(sum(x^2))
      }
  
      ### Calculate the final result
      A=T%*%t(T)
      L<-t(chol(sigma)) #### Choleski factorization
      B=L%*%A%*%t(L)
      return(B)
    }
    ## Special case when d=1
    else {
      x<-rnorm(n,mean=0,sd=sigma)
      B<-sum(x^2)
      return(B)
    }
  }
  else warning("Input Error")
}

## -----------------------------------------------------------------------------
# General samples
set.seed(123) ## Omit this line to generate different samples if necessary
sigma<-matrix(c(1.2,.8,.8,.8,1.3,.8,.8,.8,1.4),nrow=3,ncol=3)
X1<-Wis.sample(3,5,sigma)

# Special sample with d=1
set.seed(1234)
sigma0<-5
X2<-Wis.sample(1,3,sigma0)

# Output
print(X1)
print(X2)

## -----------------------------------------------------------------------------
set.seed(1)
x<-runif(10000,0,pi/3)
MC<-mean(pi*sin(x)/3) # Monte Carlo integration
I<-1-cos(pi/3) # Exact value
print(c(MC,I))

## -----------------------------------------------------------------------------
set.seed(2)
x<-runif(10000,0,1)
y<-1-x
f<-function(x){
  y<-exp(x)/(1+x^2)
  return(y)
}
MC0<-mean(f(x))
MC<-mean(f(x)+f(y))/2
print(c(MC0,MC))
print(c(var(f(x)),var((f(x)+f(y))*0.5),100-100*var((f(x)+f(y))*0.5)/var(f(x))))

## -----------------------------------------------------------------------------
f<-function(x){          # Integrated function f
  y<-exp(-x)/(1+x^2)
  return(y)
}

g<-function(x,j){          # Inverse function of Gj
  y<--log(1-(1-exp(-1))*(x+j-1)/5)
  return(y)
}

h<-function(x){          # Density function of samples
  y<-5*exp(-x)/(1-exp(-1))
  return(y)
}

set.seed(3)
x<-y<-matrix(nrow=2000,ncol=5)
for (j in 1:5) {
  x[,j]<-runif(2000,0,1)
  y[,j]<-g(x[,j],j) # Random numbers from population Gj
}
z<-f(y)/h(y)   # Integrated items
I<-colMeans(z) # Integration on each subinterval
MC<-sum(I) # theta.hat
sd<-sd(z)  # standard variance of integrated items
k<-100*(1-sd/0.0970314) # Percentage of reduction in standard variance
print(c(MC,sd,k))

## -----------------------------------------------------------------------------
# t-interval
n<-20 # Required sample size
a<-0.05 # alpha
m<-10000 # Times to repeat MC experiment
I1<-c()
for (i in 1:m) {
  x<-rchisq(n,2)
  CI<-mean(x)+c(qt(a/2,n-1)*sd(x)/sqrt(n),-qt(a/2,n-1)*sd(x)/sqrt(n)) # Confidence interval
  I1[i]<-2<CI[2]&&2>CI[1]
}
p1.hat<-mean(I1)

# Simulation based on example 6.4
I2<-c()
for (i in 1:m){
  x<-rnorm(n,0,2)
  UCL<-(n-1)*var(x)/qchisq(a,n-1)
  I2[i]<-4<UCL
}
p2.hat<-mean(I2)
cbind(p1.hat,p2.hat)

## -----------------------------------------------------------------------------
# Sampling
n<-50 # Sample size to calculate skewness
k<-300 # Number of skewnesses to calculate their quantile
m<-10000 # Times to repeat MC experiment

a<-c(0.025,0.05,0.95,0.975) # alpha of quantiles

moment<-function(x,n){  # Function to calculate central moment
  y<-mean((x-mean(x))^n)
  return(y)
}

quantile.sk<-function(n,k){ # Function to calculate the quantile of skewness
  x<-matrix(nrow=k,ncol=n)
  skewness<-numeric(k)
  
  for (i in 1:k) {
    x[i,]<-rnorm(n,0,1) # Samples to calculate skewness
    skewness[i]<-moment(x[i,],3)/(moment(x[i,],2)^1.5)
  }
  q.est<-quantile(skewness,a) # Quantiles of skewness
  return(q.est)
}

q.est<-matrix(nrow=m,ncol=4)
for (j in 1:m) q.est[j,]<-quantile.sk(n,k) # m times of MC experiments for m quantiles

# Estimation
sd.sample<-apply(q.est,2,sd) # Actual standard error of quantiles
sd.formula<-sqrt(a*(1-a)/(n*dnorm(qnorm(a,0,sqrt(6/n)),0,sqrt(6/n))^2)) # Standard error calculated with the formula
rbind(sd.sample,sd.formula)

q.hat<-colMeans(q.est) # MC estimates of quantiles
q.norm<-qnorm(a,0,sqrt(6/n)) # Large sample approximation
rbind(q.hat,q.norm)

## -----------------------------------------------------------------------------
# Sample skewness function
sk<-function(x){
  m3<-mean((x-mean(x))^3)
  m2<-mean((x-mean(x))^2)
  y<-m3/m2^1.5
  return(y)
}

# Generate sample from Beta and t
n<-200 # Sample size
m<-5000 # Replication of MC experiment
theta<-c(5,10,15,20,25) # Parameters

p1<-p2<-numeric(5)
for (i in theta){
  R1<-R2<-numeric(m)
  for (j in 1:m){
    x1<-rbeta(n,i,i)
    x2<-rt(n,i)
    b1<-sk(x1)
    b2<-sk(x2)
    z<-qnorm(0.975,0,sqrt(6/n))
    R1[j]<-mean(abs(b1)>z)
    R2[j]<-mean(abs(b2)>z)
  }
  p1[i/5]<-mean(R1)
  p2[i/5]<-mean(R2)
}

# Power curve
par(mfrow=c(1,2))
plot(theta,p1,xlab="a",ylab="power",main="Beta(a,a)",type="l")
plot(theta,p2,xlab="v",ylab="power",main="t(v)",type="l")

## -----------------------------------------------------------------------------
m<-10000 # Replications of MC experiment
n<-100 # Sample size in each MC experiment
u<-1 # Population means are equivalent

# Chisq(1)
error.C<-replicate(m,expr={
  x<-rchisq(n,1)
  t<-sqrt(n)*(mean(x)-u)/sd(x)
  reject<-(abs(t)>qt(0.975,n-1))
})
Chi.rate<-mean(error.C)

# U(0,2)
error.U<-replicate(m,expr={
  x<-runif(n,0,2)
  t<-sqrt(n)*(mean(x)-u)/sd(x)
  reject<-(abs(t)>qt(0.975,n-1))
})
U.rate<-mean(error.U)

# Exp(1)
error.E<-replicate(m,expr={
  x<-rexp(n,1)
  t<-sqrt(n)*(mean(x)-u)/sd(x)
  reject<-(abs(t)>qt(0.975,n-1))
})
Exp.rate<-mean(error.E)

cbind(Chi.rate,U.rate,Exp.rate)

## -----------------------------------------------------------------------------
library(bootstrap)
x<-scor
Corr=cor(x) # Corrlation.hat
print(Corr)
plot(x)

n<-nrow(x)
B<-1e4
Corr.star<-matrix(nrow=B,ncol=4)

for (i in 1:B) {
index<-sample(1:n,n,replace=TRUE)
x.star<-x[index,] # Bootstrap sample
Cor<-cor(x.star)
Corr.star[i,]<-c(Cor[1,2],Cor[3,4],Cor[3,5],Cor[4,5])
}

se.hat<-apply(Corr.star,2,sd)
knitr::kable(data.frame(value=se.hat,row.names=c("se12.hat","se34.hat","se35.hat","se45.hat")))

## -----------------------------------------------------------------------------
n<-100 # Sample size
M<-1e3 # Repeats for MC experiments
B<-1e3 # Times of Bootstrap experiments in each MC experiment

skew<-function(x){ # Function to calculate sample skewness
  m3<-mean((x-mean(x))^3)
  m2<-mean((x-mean(x))^2)
  y<-m3/m2^1.5
  return(y)
}

set.seed(1)
I.normal<-I.basic<-I.percentile<-matrix(nrow=M,ncol=2) # Coverage indicator

for (i in 1:M){
  x<-rnorm(n) # Normal population
  y<-rchisq(n,5) # Chisquare population
  skew.x<-skew(x)
  skew.y<-skew(y)
  skstar.x<-skstar.y<-numeric(B)
  
  for (j in 1:B){
    x.star<-sample(x,n,replace=TRUE)
    y.star<-sample(y,n,replace=TRUE)
    skstar.x[j]<-skew(x.star)
    skstar.y[j]<-skew(y.star)
  }
  
  x.n<-qnorm(0.975)*sd(skstar.x)
  y.n<-qnorm(0.975)*sd(skstar.y)
  
  x.b<-quantile(skstar.x,c(0.975,0.025))
  y.b<-quantile(skstar.y,c(0.975,0.025))
  
  I.normal[i,]<-c((skew.x-x.n)<=0 & 0<=(skew.x+x.n),
                  (skew.y-y.n)<=sqrt(8/5) & sqrt(8/5)<=(skew.y+y.n))
  
  I.basic[i,]<-c(2*skew.x-x.b[1]<=0 & 0<=2*skew.x-x.b[2],
                 2*skew.y-y.b[1]<=sqrt(8/5) & sqrt(8/5)<=2*skew.y-y.b[2])
  
  I.percentile[i,]<-c(x.b[2]<=0 & 0<=x.b[1],
                      y.b[2]<=sqrt(8/5) & sqrt(8/5)<=y.b[1])
}

prob.hat<-rbind(colMeans(I.normal),colMeans(I.basic),colMeans(I.percentile))
colnames(prob.hat)<-c("N(0,1)","Chisquare(5)")
out<-as.data.frame(prob.hat,row.names<-c("Normal","Basic","Percentile"))
knitr::kable(out)

## -----------------------------------------------------------------------------
library(bootstrap)
x<-scor
n<-nrow(x)
lambda.hat<-eigen(cov(x))$values
theta.hat<-lambda.hat[1]/sum(lambda.hat)

theta<-numeric(n)
for (i in 1:n) {
  Cov<-cov(x[-i,])
  lambda<-eigen(Cov)$values
  theta[i]<-lambda[1]/sum(lambda)
}

theta.jack<-mean(theta)
bias.jack<-(n-1)*(theta.jack-theta.hat)
var.jack<-var(theta)*(n-1)^2/n
sd.jack<-sqrt(var.jack)
out<-data.frame(bias.jack,sd.jack)

knitr::kable(out)

## -----------------------------------------------------------------------------
library(DAAG)
library(lattice)
attach(ironslag)
n<-length(magnetic)

e.qubic<-numeric(n)

for (i in 1:n){
  x<-chemical[-i]
  y<-magnetic[-i]
  L.qubic<-lm(y~x+I(x^2)+I(x^3))
  attach(L.qubic)
  e.qubic[i]<-magnetic[i]-(coefficients[1]+coefficients[2]*chemical[i]+coefficients[3]*chemical[i]^2+coefficients[4]*chemical[i]^3)
  detach(L.qubic)
}

d.qubic<-mean(e.qubic^2) # error.hat for qubic polynomial model

m.l<-lm(magnetic ~ chemical) # linear
m.q<-lm(magnetic ~ chemical + I(chemical^2)) # quadratic
m.e<-lm(log(magnetic) ~ chemical) # exponential
m.qubic<-lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3)) # qubic

error.hat<-data.frame(Linear=c(19.55644),Quadratic=c(17.85248),
                Exponential=c(18.44188),Qubic=d.qubic)

adj.R<-c(summary(m.l)$adj.r.squared, # Linear model
         summary(m.q)$adj.r.squared, # Quadratic model
         summary(m.e)$adj.r.squared, # Exponential model
         summary(m.qubic)$adj.r.squared) # Qubic polynomial model

out<-rbind(error.hat,adj.R)
row.names(out)<-c("Cross validation","Adjusted R square")
knitr::kable(out)

## -----------------------------------------------------------------------------
# Sample size
n<-20
m<-30
N<-n+m
n0<-max(c(n,m))

# MC experiment for permutation test
# Equal variance
set.seed(12345)
x1<-rnorm(n,0,1)
y1<-rnorm(m,0,1)
z1<-c(x1,y1)

# Unequal variance
x2<-rnorm(n,0,1)
y2<-rnorm(m,0,2)
z2<-c(x2,y2)

count5test<-function(x,y) {
  X<-x-mean(x)
  Y<-y-mean(y)
  outx<-sum(X>max(Y))+sum(X<min(Y))
  outy<-sum(Y>max(X))+sum(Y<min(X))
  return(max(c(outx,outy))>5)
}

t1<-x1  # The comparison is conducted between larger-sized pupulation
t2<-x2
if (n<m) {
  t1<-y1
  t2<-y2
}

B<-10000
I1<-I2<-numeric(B)
for (i in 1:B){
  k<-sample(1:N,n0)
  I1[i]<-count5test(t1,z1[k])
  I2[i]<-count5test(t2,z2[k])
}

p1<-mean(I1)
p2<-mean(I2)
out<-cbind(p1,p2)
colnames(out)<-c("Equal variace","Unequal variance")
knitr::kable(out)

## -----------------------------------------------------------------------------
# Necessary Function
dCor<-function(z,ix){
n<-nrow(x)

D<-function(x){
d<-as.matrix(dist(x))
m<-rowMeans(d)
M<-mean(d)
a<-sweep(d,1,m)
b<-sweep(a,2,m)
return(b+M)
}

A<-D(z[,1:2])
B<-D(z[ix,3:4])
dCov<-sqrt(mean(A*B))
dVarX<-sqrt(mean(A*A))
dVarY<-sqrt(mean(B*B))
dCor<-sqrt(dCov/sqrt(dVarX*dVarY))
return(dCor)
}

# Testing
library(Ball)
library(mnormt)
library(boot)

k<-15
n<-k*(1:10) # Sample size
m<-99 # Times of permutation tests
B<-100 # Times of MC experiments
n0<-length(n)
pv1.dCor<-pv2.dCor<-pv1.ball<-pv2.ball<-numeric(n0) # Empirical power

for (i in n){
  x<-rmnorm(i,mean=rep(0,2),varcov=diag(1,2,2))
  e<-rmnorm(i,mean=rep(0,2),varcov=diag(1,2,2))
  y1<-x/4+e
  y2<-x/4*e
  z1<-cbind(x,y1)
  z2<-cbind(x,y2)
  
  p1.dCor<-p2.dCor<-p1.ball<-p2.ball<-numeric(B)
  for(j in 1:B){
  out1.dCor<-boot(z1,statistic=dCor,R=m,sim="permutation")
  out2.dCor<-boot(z2,statistic=dCor,R=m,sim="permutation")
  I1<-c(out1.dCor$t0,out1.dCor$t)
  I2<-c(out2.dCor$t0,out2.dCor$t)
  
  p1.dCor[j]<-mean(I1>=out1.dCor$t0)
  p2.dCor[j]<-mean(I2>=out2.dCor$t0)
  p1.ball[j]<-bcov.test(x,y1,R=m,seed=i*j)$p.value
  p2.ball[j]<-bcov.test(x,y2,R=m,seed=i*j)$p.value
  }
  
pv1.dCor[i/k]<-mean(p1.dCor<0.05)
pv2.dCor[i/k]<-mean(p2.dCor<0.05)
pv1.ball[i/k]<-mean(p1.ball<0.05)
pv2.ball[i/k]<-mean(p2.ball<0.05)
}

par(mfrow=c(1,2))
plot(n,pv1.dCor,xlab="n",ylab="Power",main="Model 1",type="l",ylim=c(0,1))
lines(n,pv1.ball,type="l",col="blue")

plot(n,pv2.dCor,xlab="n",ylab="Power",main="Model 2",type="l",ylim=c(0,1))
lines(n,pv2.ball,type="l",col="blue")

## -----------------------------------------------------------------------------
# Random walk Metropolis sampler
rw<-function(sigma,N){
  x<-numeric(N)
  u<-runif(N)
  x[1]<-20
  k<-0
  
  for (i in 2:N){
    y<-rnorm(1,x[i-1],sigma)
    
    if (u[i]<=(exp(-1*abs(y))/exp(-1*abs(x[i-1])))) {
      x[i]<-y
      k<-k+1
    }
    
    else x[i]<-x[i-1]
  }
  return(list(x=x,rate=k/N))
}

# Sampling from different sigma
sigma<-c(0.05,0.5,1,5)
N<-5000
x<-list(NULL)
length(x)<-4
rate<-numeric(4)
set.seed(111)
for (i in 1:4){
  out<-rw(sigma[i],N)
  x[[i]]<-out$x
  rate[i]<-out$rate
}

acc<-data.frame(sigma=sigma,Acceptance_rate=rate)
knitr::kable(acc)

for (i in 1:4) plot(1:N,x[[i]],type="l",xlab="n",ylab="X",main=paste("σ=",sigma[i]))

## -----------------------------------------------------------------------------
x<-seq(2,3,0.1)
y1<-log(exp(x))
y2<-exp(log(x))

all(y1==y2)
all.equal(y1,y2)

## -----------------------------------------------------------------------------
B<-function(a,k){
  g<-function(a,k){
    
    f<-function(u,k){
      f<-(k/(k+u^2))^(0.5*(k+1))
      return(f)
    }
    
    ck<-sqrt(a^2*k/(k+1-a^2))
    c<-2*exp(lgamma((k+1)/2)-lgamma(k/2))/sqrt(pi*k)
    I<-integrate(f,lower=0,upper=ck,k=k)$value
    g<-c*I
    return(g)
  }
  
  return(g(a,k-1)-g(a,k))
}

k<-as.integer(c(4:25,100,500,1000))
ak<-numeric(length(k))
for (i in 1:length(k)) ak[i]<-uniroot(B,lower=0.01,upper=2,k=k[i])$root

## -----------------------------------------------------------------------------
A<-function(a,k){
  sk<-function(a,k){
    ck<-sqrt(a^2*k/(k+1-a^2))
    sk<-1-pt(ck,k)
  }
  A<-sk(a,k-1)-sk(a,k)
  return(A)
}

Ak<-numeric(length(k))
for (i in 1:length(k)) Ak[i]<-uniroot(A,c(0.01,sqrt(k[i]-0.01)),k=k[i])$root

out<-rbind(k,round(Ak,3),round(ak,3))
row.names(out)<-c("k","A(k)","a(k)")
out<-as.data.frame(out)
knitr::kable(out)

## -----------------------------------------------------------------------------
l<-function(theta=numeric(2),pt,qt){
  na<-28;nb<-24;no<-41;nab<-70
  p<-theta[1]
  q<-theta[2]
  rt<-1-pt-qt
  c1<-pt^2/(pt^2+2*pt*rt)
  c2<-qt^2/(qt^2+2*qt*rt)
  
  l<-2*c1*na*log(p)+2*c2*nb*log(q)+2*no*log(1-p-q)+na*log(2*p*(1-p-q))-c1*na*log(2*p*(1-p-q))+nb*log(2*q*(1-p-q))-c2*nb*log(2*q*(1-p-q))+nab*log(2*p*q)
return(-l)
}

I<-1
lt<-numeric()
theta<-c(0.4,0.4)

while(I){
  theta0<-theta
  MLE<-optim(par=c(0.1,0.1),fn=l,pt=theta0[1],qt=theta0[2])
  theta<-MLE$par
  lt<-c(lt,-1*MLE$value)
  if(max(abs(theta0-theta))<1e-7) I<-0
}

theta.hat<-data.frame(p=theta[1],q=theta[2])
plot(1:length(lt),lt,xlab="t",ylab="log L",type="l")
knitr::kable(theta.hat)

## -----------------------------------------------------------------------------
formulas<-list(
mpg~disp,
mpg~I(1/disp),
mpg~disp+wt,
mpg~I(1/disp)+wt
)

# for loops
n<-length(formulas)
fit3.loop<-list(NULL)
length(fit3.loop)<-n

for (i in 1:n) fit3.loop[[i]]<-lm(formulas[[i]],mtcars)

# lapply()
fit3.lapply<-lapply(formulas,function(x) lm(x,mtcars))

print(fit3.loop)
print(fit3.lapply)

## -----------------------------------------------------------------------------
bootprintaps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})

n<-length(bootprintaps)

# for loops
fit4.loop<-list(NULL)
length(fit4.loop)<-n
for (i in 1:n) fit4.loop[[i]]<-lm(mpg~disp,data=bootprintaps[[i]])

# lapply()
fit4.lapply<-lapply(bootprintaps,function(x) lm(mpg~disp,data=x))

# lapply() without anonymous function
fit<-function(x) {
  out<-lm(mpg~disp,x)
  return(out)
}

fit4.named<-lapply(bootprintaps,fit)

print(fit4.loop)
print(fit4.lapply)
print(fit4.named)

## -----------------------------------------------------------------------------
rsq<-function(mod) summary(mod)$r.squared

model<-list(fit3.loop,fit3.lapply,fit4.loop,fit4.lapply,fit4.named)
rsquare<-lapply(model,function(x) lapply(x,rsq))

print(rsquare)

## -----------------------------------------------------------------------------
trials <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)

p.value<-sapply(trials, function(x) x$p.value)

# Extra answer
p.extra<-sapply(trials,"[[","p.value")

print(p.value)
print(p.extra)

## -----------------------------------------------------------------------------
library(parallel)
nc<-detectCores()
cl<-makeCluster(2)

n<-1:10000
fun<-function(n) sample(1:n,n,replace=T)

# mcsapply
mcsapply<-function(x,f) parSapply(cl,x,f)

# Compare runtime
system.time(mcsapply(n,fun))
system.time(sapply(n,fun))

stopCluster(cl)

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)
library(knitr)

# Random walk sampler in C
cppFunction(
'List rwC(double sigma, int N){
  NumericVector x(N); double y; double u; int k=0;

  x[0]=20;

  for(int i=1;i<=(N-1);i++){
    y=rnorm(1,x[i-1],sigma)[0];
    u=runif(1)[0];

    if (u<=(exp(-1*abs(y))/exp(-1*abs(x[i-1])))) {
    x[i]=y;
    k++;
    }
    
    else x[i]=x[i-1];
  }
  
  List out=List::create(x,k);
  return (out);
}')

# Random walk sampler in R
rwR<-function(sigma,N){
  x<-numeric(N)

  x[1]<-20
  k<-0
  
  for (i in 2:N){
    y<-rnorm(1,x[i-1],sigma)
    u<-runif(1)
    
    if (u<=(exp(-1*abs(y))/exp(-1*abs(x[i-1])))){
    x[i]<-y
    k<-k+1
    }
    
    else x[i]<-x[i-1]
  }
  return (list(x=x,k))
}

# Sampling from different sigma
sigma<-c(0.05,0.5,1,5)
N<-5000
xC<-xR<-list(NULL)
length(xC)<-length(xR)<-4
rateC<-rateR<-numeric(4)

set.seed(123)
for (i in 1:4){
  outC<-rwC(sigma[i],N)
  outR<-rwR(sigma[i],N)
  xC[[i]]<-outC[[1]]
  xR[[i]]<-outR[[1]]
  rateC[i]<-outC[[2]]/N
  rateR[i]<-outR[[2]]/N
}

kable(data.frame(sigma=sigma,rate.Rcpp=rateC,rate.R=rateR))
for (i in 1:4) {
  par(mfrow=c(1,2))
  plot(1:N,xC[[i]],type="l",xlab="n",ylab="X",main=paste("Rcpp: σ=",sigma[i]))
  plot(1:N,xR[[i]],type="l",xlab="n",ylab="X",main=paste("R: σ=",sigma[i]))
}

# Comparison
## Q-Q plot
for (i in 1:4) {
  qqplot(xC[[i]],xR[[i]],xlab="Rcpp",ylab="R",main=paste("Q-Q plot:σ=",sigma[i]),type="l")
  abline(c(1,1),col="blue")
  }

# Microbenchmark: only one comparison when sigma=0.05
# Because different sigma rendered similar runtime
runtime<-microbenchmark(rw.R=rwR(sigma[1],N),rw.C=rwC(sigma[1],N))
summary(runtime)

