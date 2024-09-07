packageName <- c("dfphase1","DescTools","moments","tolerance","lmomco")
for(i in 1:length(packageName)) {
  if(!(packageName[i] %in% rownames(installed.packages()))) {
    install.packages(packageName[i])
  }
}
lapply(packageName, require, character.only = TRUE)
library(dfphase1)
library(DescTools)
library(moments)
library(tolerance)
library(lmomco)

# Read the raw data
dt = read.csv("uci-secom.csv", header = F)

# There are 591 variables in the dataset.
# In the CSV file, the first column represents time stamps, 
# and the subsequent columns represent the quality characteristics.


prepro = function(x){
  #filters out the missing values
  xfull = x[!is.na(x)]
  xtrim=Trim(xfull, trim = 0.01)
  xfinal = xtrim
  return(xfinal)
}

# The 2nd, 25th, 158th variables are selected and 
# then named as x1, x2, x3, respectively.

x = x1 = prepro(dt$V2);x=x[1:53]
# x= x2 = prepro(dt$V25);x=x[1:255]
# x= x3 = prepro(dt$V190);x=x[1:1500]
# x= x4 = prepro(dt$V226);x=x[1:1000]

variable="X1"
Nsim=n=length(x)
mean=round(mean(x),2)
stdev=round(sd(x),2)
skewness=round(skewness(x),2)
kurtosis=round(kurtosis(x),2)
moments <- c(mean,stdev,skewness,kurtosis);moments

Set.Confid=0.7  # 1-p
z=qnorm((1+Set.Confid)/2);z
alpha=(1-pnorm(3));alpha # alpha= upper 3 sigma;
z.alpha=-qnorm(alpha);z.alpha
par(mfrow=c(1,1))

###############################################################
# Nonparametric Bernstein FOS method
###############################################################
Ber = function(x, pn = Set.Confid, alpha = 0.0027){
  library(lmomco)
  Rx = sort(x)
  n = length(x)
  n1=n+1
  
  # # one of two-sided directly
  fun = function(r, n, alpha = alpha, pn = pn ) pbeta(1-alpha, n-2*r+1, 2*r)-(1-pn)
  kr <- uniroot(fun, c(0, n1/2), n = n, alpha = alpha, pn = pn)$root
  urI=kr/n1
  usI=1-urI
  limit2=function(u){dat2bernqua(u, x, poly.type="Bernstein", bound.type = "Carv")}
  LCL2 = limit2(urI)
  UCL2 = limit2(usI)
  
  # # two of one-sided indirectly 
  findkl = function(r, n, alpha = alpha, pn = pn ) pbeta(alpha/2, r, n-r+1)-(1+pn)/2
  findku = function(r, n, alpha = alpha, pn = pn ) pbeta(1-alpha/2, r, n-r+1)-(1-pn)/2
  kl <- uniroot(findkl, c(0, n1), n = n, alpha = alpha, pn = pn)$root
  ku <- uniroot(findku, c(0, n1), n = n, alpha = alpha, pn = pn)$root
  ur = kl/n1
  us = ku/n1
  lu=dat2bernqua(f=c(1/n1,1-1/n1), x, poly.type="Bernstein", bound.type = "Carv")
  pp0=function(n){pnorm(n/m[1]-1, mean =0, sd =sdd[1])}
  p0=1
  limit1 = function(u){
    if (0<u && u<=1/(p0*n1)) {
      lu[1]+(Rx[2]-Rx[1])*log(p0*n1*u)-(Rx[2]-Rx[1])*(log(p0*n1*u))^2
    } else if (1/(p0*n1)<u && u<1-1/(p0*n1)) {
      dat2bernqua(u, x, poly.type="Bernstein", bound.type = "Carv")
    } else {
      lu[2]-(Rx[n]-Rx[n-1])*log(p0*n1*(1-u))+(Rx[n]-Rx[n-1])*(log(p0*n1*(1-u)))^2
    }
  }
  LCL1 = limit1(ur)
  UCL1 = limit1(us)
  
  mm=1/(p0*(alpha/2))-1
  
  LCL=if (n>mm) {
    LCL2
  } else {
    LCL1
  }
  
  UCL=if (n>mm) {
    UCL2
  } else {
    UCL1
  }
  
  return(list(UCL=UCL,LCL=LCL))
}

UCL.Q=Ber(x)$UCL
LCL.Q=Ber(x)$LCL

# COMBINED SPC CHARTS

LCL.FOS=LCL.Q
UCL.FOS=UCL.Q

LD=min(x,LCL.FOS)
UD=max(x,UCL.FOS)
lss=abs(UD-LD)/8
uss=abs(UD-LD)/1.5
plot(x,ylim=c(LD-lss,UD+uss),xlab="sample sequence", ylab="values",
     main=paste0(variable,": mean=",moments[1],", stdev=",moments[2],", skew=",moments[3],", kurt=",moments[4], ", n=",Nsim," (Pn=",round(Set.Confid,3),")"),cex.main=1.5,cex = 0.8,pch=20,lty=1)
lines(x)


abline(h=UCL.FOS,col="blue",lty=2,lwd=3) 
abline(h=LCL.FOS,col="blue",lty=2,lwd=3)
abline(h=mean,col="black",lty=2,lwd=1)

legend(x = "topright", legend = c("Bernstein: [UCL,LCL]"),
       col = c("blue"),
       lty = c(2), cex = 1.2, box.lty = 1,lwd=c(3))


#   UCL LCL
text(Nsim/5,UCL.FOS+1.2,round(UCL.FOS,digits = 2),cex = 1, pos = 3, col ="blue")
text(Nsim/5,LCL.FOS-1.2,round(LCL.FOS,digits = 2),cex = 1, pos = 1, col ="blue")



CL=rbind(LCL.FOS,UCL.FOS,UCL.Goed,LCL.Goed)
result=rbind(variable,n,CL,Set.Confid);result


