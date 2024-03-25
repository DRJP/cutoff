## remotes::install_github("DRJP/cutoff/cutoff2@penalties")
## remotes::install_github("DRJP/cutoff/cutoff2@penalties", build_manual=TRUE, build_vignettes=TRUE)
library(cutoff2)
rm(list=ls())

set.seed(4)
par(mfrow=c(2,2))
mu = c(4, 6)
sd = c(1, 4)
w  = c(0.5, 0.5)

#######################
## Plot 1: the model ##
#######################
xRange = range(qnorm(p=c(0.001,0.001,0.999,0.999),mu,sd))
xRange[1] = floor(xRange[1])
xRange[2] = ceiling(xRange[2])
curve(w[1]*dnorm(x,mu[1],sd[1]) + w[2]*dnorm(x,mu[2],sd[2]), xRange[1], xRange[2], n=1111,
      lwd=2, ylab="Density" ,xlab="MFI")
curve(w[1]*dnorm(x,mu[1],sd[1]), xRange[1], xRange[2], n=1111, lty=2, lwd=2 ,add=TRUE, col="blue")
curve(w[2]*dnorm(x,mu[2],sd[2]), xRange[1], xRange[2], n=1111, lty=2, lwd=2, add=TRUE, col="red")
title("Biological nonsense model")
legend(legend=c("model", "pos","neg"),"topright", lwd=c(2,2,2), col=c("black","red","blue"),
       bty="n" ,lty=c(1,2,2))
print(w)

######################
## Plot 2: the data ##
######################
n = 100
n1 = rbinom(1,n,w[1]); n2 = n - n1
y = c(rnorm(n1,mu[1],sd[1]), rnorm(n2,mu[2],sd[2]))
miny = floor(min(y))
maxy = ceiling(max(y))
hist(y, freq=FALSE, breaks=seq(miny, maxy, by=0.5), xlab="MFI", main="Simulated data")
curve(w[1]*dnorm(x,mu[1],sd[1])+w[2]*dnorm(x,mu[2],sd[2]),miny,maxy,add=TRUE,col="black", lwd=3)
curve(w[1]*dnorm(x,mu[1],sd[1]),miny,maxy,add=TRUE,col="black", lwd=2, lty=2)
curve(w[2]*dnorm(x,mu[2],sd[2]),miny,maxy,add=TRUE,col="black", lwd=2, lty=2)

#######################
## Unconstrained Fit ##
#######################
# Estimate parameters of finite mixture model:
(fit1 <- em(y,"normal","normal"))
hist(y, freq=FALSE, breaks=seq(miny, maxy, by=0.5), xlab="MFI", main="Unconstrained fit")
# Add the EM estimated finite mixture model:
lines(fit1, col="tomato", lwd=2)
# Estimate a cutoff from the fitted mixture model
(cut_off <- cutoff(fit1, whose="Titterington", nb=1000))
polygon(c(cut_off[-1],rev(cut_off[-1])),c(0,0,.55,.55), col=rgb(1, 0.39, 0.28,.2),border=NA)
abline(v=cut_off[-1],lty=2,col="tomato")
abline(v=cut_off[1],col="tomato")

###################
## Penalised fit ##
###################
# Estimate parameters of finite mixture model:
(fit2 <- em(y,"normal","normal", penaltyScale=1E8))
# Replot data
hist(y, freq=FALSE, breaks=seq(miny, maxy, by=0.5), xlab="MFI", main="Penalised fit")
# Add the penalised-EM estimated finite mixture model:
lines(fit2, col="red", lwd=2)
# Estimate a cutoff from the fitted mixture model
(cut_off <- cutoff(fit2, whose="Titterington", nb=1000))
polygon(c(cut_off[-1],rev(cut_off[-1])),c(0,0,.55,.55), col=rgb(1,0,0,.2),border=NA)
abline(v=cut_off[-1],lty=2,col="red")
abline(v=cut_off[1],col="red")

plot(y, fit2$pPositive)


#############################
## Add a shift to the data ##
#############################
yy = y - min(y) + 0.1
minyy = floor(min(yy)); maxyy=ceiling(max(yy))
(fit3 <- em(yy,"log-normal","log-normal", penaltyScale=0))
(fit4 <- em(yy,"log-normal","log-normal", penaltyScale=1E4))
(fit5 <- em(yy,"log-normal","log-normal", penaltyScale=1E8))
(fit5b <- em(yy,"log-normal","log-normal", penaltyScale=1E20))
# Replot data
par(mfrow=c(2,2))
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", main="Penalised fit"); lines(fit3, col="blue", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", main="Penalised fit"); lines(fit4, col="green", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", main="Penalised fit"); lines(fit5, col="red", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", main="Penalised fit"); lines(fit5b, col="red", lwd=2)


(fit6 <- em(yy,"weibull","weibull", penaltyScale=0))
(fit7 <- em(yy,"weibull","weibull", penaltyScale=1E4))
(fit8 <- em(yy,"weibull","weibull", penaltyScale=1E8))
# Replot data
par(mfrow=c(1,3))
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", main="Unconstrained fit"); lines(fit6, col="blue", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", main="Penalty = 1E4"); lines(fit7, col="green", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", main="Penalty = 1E8"); lines(fit8, col="red", lwd=2)


yy = y - min(y) + 0.1
penScale = 1E6
par(mfrow=c(4,4))
(fit1  <- em(yy,"normal","normal", penScale, TRUE))
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit1$deviance), main="normal-normal");           lines(fit1, col="blue", lwd=2)
(fit2  <- em(yy,"normal","weibull", penScale, TRUE))
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit2$deviance), main="normal-weibull");          lines(fit2, col="blue", lwd=2)
(fit3  <- em(yy,"normal","gamma", penScale, TRUE))
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit3$deviance), main="normal-gamma");            lines(fit3, col="blue", lwd=2)
(fit4  <- em(yy,"normal","log-normal", penScale, TRUE))
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit4$deviance), main="normal - log-normal");     lines(fit4, col="blue", lwd=2)
(fit5  <- em(yy,"weibull","normal", penScale, TRUE))
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit5$deviance), main="weibull-normal");          lines(fit5, col="blue", lwd=2)
(fit6  <- em(yy,"weibull","weibull", penScale, TRUE))
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit6$deviance), main="weibull-weibull");         lines(fit6, col="blue", lwd=2)
(fit7  <- em(yy,"weibull","gamma", penScale, TRUE))
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit7$deviance), main="weibull-gamma");           lines(fit7, col="blue", lwd=2)
(fit8  <- em(yy,"weibull","log-normal", penaltyScale=1E9))
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit8$deviance), main="weibull - log-normal");    lines(fit8, col="blue", lwd=2)
(fit9  <- em(yy,"gamma","normal", penScale, TRUE))
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit9$deviance), main="gamma-normal");            lines(fit9, col="blue", lwd=2)
(fit10 <- em(yy,"gamma","weibull", penScale, TRUE))
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit10$deviance), main="gamma-weibull");           lines(fit10, col="blue", lwd=2)
(fit11 <- em(yy,"gamma","gamma", penScale, TRUE))
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit11$deviance), main="gamma-gamma");             lines(fit11, col="blue", lwd=2)
(fit12 <- em(yy,"gamma","log-normal", penScale, TRUE, thresh=exp(-25)))
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit12$deviance), main="gamma-log-normal");        lines(fit12, col="blue", lwd=2)
(fit13 <- em(yy,"log-normal","normal", penScale, TRUE))
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit13$deviance), main="log-normal - normal");     lines(fit13, col="blue", lwd=2)
(fit14 <- em(yy,"log-normal","weibull", penScale, TRUE))
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit14$deviance), main="log-normal - weibull");    lines(fit14, col="blue", lwd=2)
(fit15 <- em(yy,"log-normal","gamma", penScale, TRUE, thresh=exp(-25)))
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit15$deviance), main="log-normal-gamma");        lines(fit15, col="blue", lwd=2)
(fit16 <- em(yy,"log-normal","log-normal", penScale, TRUE))
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit16$deviance), main="log-normal - log-normal"); lines(fit16, col="blue", lwd=2)

penScale = 1E10 # 1E6
fit10
(fit10 <- em(yy,"gamma","weibull", penaltyScale=penScale, TRUE))

par(mfrow=c(4,4))
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit1$deviance), main="normal-normal");           lines(fit1, col="blue", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit2$deviance), main="normal-weibull");          lines(fit2, col="blue", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit3$deviance), main="normal-gamma");            lines(fit3, col="blue", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit4$deviance), main="normal - log-normal");     lines(fit4, col="blue", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit5$deviance), main="weibull-normal");          lines(fit5, col="blue", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit6$deviance), main="weibull-weibull");         lines(fit6, col="blue", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit7$deviance), main="weibull-gamma");           lines(fit7, col="blue", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit8$deviance), main="weibull - log-normal");    lines(fit8, col="blue", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit9$deviance), main="gamma-normal");            lines(fit9, col="blue", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit10$deviance), main="gamma-weibull");           lines(fit10, col="blue", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit11$deviance), main="gamma-gamma");             lines(fit11, col="blue", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit12$deviance), main="gamma-log-normal");        lines(fit12, col="blue", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit13$deviance), main="log-normal - normal");     lines(fit13, col="blue", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit14$deviance), main="log-normal - weibull");    lines(fit14, col="blue", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit15$deviance), main="log-normal-gamma");        lines(fit15, col="blue", lwd=2)
hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit16$deviance), main="log-normal - log-normal"); lines(fit16, col="blue", lwd=2)
