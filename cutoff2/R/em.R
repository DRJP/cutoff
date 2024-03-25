#' Minus log-likelihood of the finite mixture model.
#'
#' @keywords internal
# This function returns minus the log-likelihood of the set of parameters "mu1",
# "sigma1", "mu2" and "sigma2", given dataset "data", distributions "D1" and
# "D2", and parameter "lambda".
mLL <- function(mu1,sigma1,mu2,sigma2,lambdaVec,data,D1,D2,d1,d2,p1,p2,q1,q2,penaltyScale) {
# "mu1", "mu2", "sigma1" and "sigma2" are the log of the parameter values.
# "lambdaVec" and "data" are two numeric vectors of the same length.
# "D1" and "D2" are two probability density functions.
# "p1" and "p2" are two cummulative distribution functions.
  # "q1" and "q2" are two quantile functions.
  params <- backTrans(D1, mu1, sigma1, D2, mu2, sigma2) # Transform from algorithm scale (-Inf, Inf) to biological scale
  names(params) <- c("mu1","sigma1","mu2","sigma2") # Could probably now remove this line, as naming done in backTrans should be OK
  lambda = mean(lambdaVec)
  # with(as.list(exp(params)), {
  with(as.list(params), {
    ## Expected values of complete data log likelihoods
    ECDLL <- (log(lambda)   + d1(data,mu1,sigma1,log=TRUE)) * lambdaVec +
             (log(1-lambda) + d2(data,mu2,sigma2,log=TRUE)) * (1 - lambdaVec)
    ## LEL <- log(lambda * d1(data,mu1,sigma1) + (1-lambda) * d2(data,mu2,sigma2))
    ## Penalties
    penalty <- 0
    if (penaltyScale > 0) {
      mean1 <- getMean(D1, mu1, sigma1)
      mean2 <- getMean(D2, mu2, sigma2)
      if (mean1 < mean2) {
        penalty <- calcPenalty(mu1,sigma1,mu2,sigma2,lambda,data,d1,d2,p1,p2,q1,q2,penaltyScale)
      } else {
        penalty <- calcPenalty(mu2,sigma2,mu1,sigma1,1-lambda,data,d2,d1,p2,p1,q2,q1,penaltyScale)
      }
    }
    return( -sum(ECDLL) - penalty)
    # return( -sum(LEL) - penalty)
  })
}


#' Generates a penalty against biologically non-sensical fits
#' to be added to the log-likelihood of the finite mixture model.
#'
#' @keywords internal
# This function returns the penalty to be added to the log-likelihood, given parameters "MU1",
# "SIGMA1", "MU2" and "SIGMA2", given dataset "data", density functions "d1" and
# "d2", parameter "LAMBDA", distribution functions p1, p2 and quantile functions q1 and q2.
calcPenalty <- function(MU1,SIGMA1,MU2,SIGMA2,LAMBDA,data,d1,d2,p1,p2,q1,q2,penaltyScale) {
  ## calcPenalty assumes MU1 < MU2
  w1 = LAMBDA
  w2 = (1-LAMBDA)
  #### Penalise LHS
  X             = seq(floor(q1(1E-11,MU1,SIGMA1)), ceiling(q1(1-1E-11,MU1,SIGMA1)), l=1111)
  iNonsenseLeft = w1*d1(X,MU1,SIGMA1) < w2*d2(X,MU2,SIGMA2)
  iNonsenseLeft = cumprod(iNonsenseLeft)==1
  ii            = sum(iNonsenseLeft)
  areaNonsenseLeft   = 0
  if (ii > 0) {
    areaNonsenseLeft = (w2*(p2(X[ii],MU2,SIGMA2))) - (w1*(p1(X[ii],MU1,SIGMA1)))
  }
  #### Penalise RHS
  X              = seq(floor(q2(1E-11,MU2,SIGMA2)), ceiling(q2(1-1E-11,MU2,SIGMA2)), l=1111)
  iNonsenseRight = w1*d1(X,MU1,SIGMA1) > w2*d2(X,MU2,SIGMA2)
  iNonsenseRight = rev(cumprod(rev(iNonsenseRight))==1)
  ii             = sum(!iNonsenseRight)
  areaNonsenseRight   = 0
  if (ii < length(X)) {
    areaNonsenseRight = w1*(1-p1(X[ii],MU1,SIGMA1)) - w2*(1-p2(X[ii],MU2,SIGMA2))
  }
  #### Total Penalty
  penalty      = penaltyScale * log(1-areaNonsenseLeft-areaNonsenseRight)
  return(penalty)
}

#-------------

#' Starting values of the parameter of a finite mixture model.
#'
#' @keywords internal
# This function calculates the starting values of parameter "lambda".
startval <- function(data,d1,d2) {
  #  require(tree) # for "tree".
  thresh <- tree::tree(data~data)$frame$yval[1]
  sel    <- data<thresh
  data1  <- data[sel]
  data2  <- data[!sel]
  lambda <- length(data1)/length(data)
  param1 <- MASS::fitdistr(data1,d1)$est
  param2 <- MASS::fitdistr(data2,d2)$est
  out    <- c(param1,param2,lambda)
  names(out) <- c("mu1","sigma1","mu2","sigma2","lambda")
  return(out)
}

#-------------

#' Expectation-Maximization estimation of a finite mixture model
#'
#' \code{em} returns points estimations of the parameters of a finite mixture
#'		model using the Expectation-Maximization (E-M) algorithm.
#'
#' The finite mixture model considered in this function is a mixture of two
#'  probability distributions that are one of the following: normal, log-normal,
#'  gamma or Weibull. Each of these distributions is defined by two parameters:
#'  a location and a scale parameter:
#' \tabular{lcc}{
#'              \tab location \tab scale \cr
#'   normal     \tab mean     \tab sd    \cr
#'   log-normal \tab meanlog  \tab sdlog \cr
#'   gamma      \tab shape    \tab rate  \cr
#'   Weibull    \tab shape    \tab scale
#' }
#' These parameters, together with the mixture parameter, are estimated by the
#'  Expection-Maximization algorithm.
#'
#' @param data A vector of real numbers, the data to model with a finite mixture model.
#' @param D1 First probability distribution in the finite mixture model.
#' @param D2 Second probability distribution in the finite mixture model. See Details.
#' @param threshold A numerical scalar indicating the value below which the E-M algorithm should stop.
#' @param penaltyScale A positive scale parameter to penalise biologically nonsensical solutions. Defaults to 0.
#' @param forceOrder Logical. If TRUE then models 1 and 2 are switched every time the maximisation step results in mu[1] > mu[2]. A warning message is printed each time this occurs.
#' @return A list with class \code{em} containing the following components:
#' 	\item{lambda}{a numerical vector of length \code{length(data)} containing, for each datum, the probability to belong to distribution \code{D1}.}
#'  \item{param}{the location (mu) and scale (sigma) parameters of the probability distributions \code{D1} and \code{D2}.}
#' 	\item{D1}{character scalar containing the name of the first probability distribution used in the finite mixture model.}
#' 	\item{D2}{character scalar containing the name of the second probability distribution used in the finite mixture model.}
#'  \item{deviance}{scalar value giving the devinace of the fit. Set in cutoff2 as twice the -ve log-likelihood (this doubling was omitted in cutoff). }
#'  \item{data}{the numerical vector of data used as input.}
#' 	\item{data_name}{character scalar containing the name of the dataset used as input.}
#' 	\item{out}{an object of class \code{mle2} that contains the maximum-likelihood estimates of parameters \code{mu1}, \code{lambda1}},
#' 	\item{threshold}{the input \code{threshold} argument value.}
#' @references
#' 	Chuong B. Do and Serafim Batzoglou (2008) What is the expectation
#'    maximization algorithm? Nature Biotechnology 26(8): 897-899.\cr
#'  \cr
#'  Peter Schlattmann (2009) Medical Applications of Finite Mixture Models.
#'    Springer-Verlag, Berlin.
#' @seealso \code{\link{confint.em}} method for calculating the confidence
#'    intervals of the parameters and \code{\link{cutoff}} for deriving a
#'    cut-off value.
#' @examples
#' # Measles IgG concentration data:
#' length(measles)
#' range(measles)
#' # Plotting the data:
#' hist(measles,100,FALSE,xlab="concentration",ylab="density",ylim=c(0,.55), main=NULL,col="grey")
#' # The kernel density:
#' lines(density(measles),lwd=1.5,col="blue")
#' # Estimating the parameters of the finite mixture model:
#' (measles_out <- em(measles,"normal","normal"))
#' # The confidence interval of the parameter estimates:
#' confint(measles_out,threshold=1e-64,nb=100,level=.95)
#' # Adding the E-M estimated finite mixture model:
#' lines(measles_out,lwd=1.5,col="red")
#' # The legend:
#' legend("topleft",leg=c("non-parametric","E-M"),col=c("blue","red"), lty=1,lwd=1.5,bty="n")
#'
#' # Example 2: using penalisation to avoid curve 2 dominating curve 1 at low values of x
#' set.seed(4)
#' par(mfrow=c(2,2))
#' mu = c(4, 6)
#' sd = c(1, 4)
#' w  = c(0.5, 0.5)
#'
#' #######################
#' ## Plot 1: the model ##
#' #######################
#' xRange = range(qnorm(p=c(0.001,0.001,0.999,0.999),mu,sd))
#' xRange[1] = floor(xRange[1])
#' xRange[2] = ceiling(xRange[2])
#' curve(w[1]*dnorm(x,mu[1],sd[1]) + w[2]*dnorm(x,mu[2],sd[2]), xRange[1], xRange[2], n=1111,
#' lwd=2, ylab="Density" ,xlab="MFI")
#' curve(w[1]*dnorm(x,mu[1],sd[1]), xRange[1], xRange[2], n=1111, lty=2, lwd=2 ,add=TRUE, col="blue")
#' curve(w[2]*dnorm(x,mu[2],sd[2]), xRange[1], xRange[2], n=1111, lty=2, lwd=2, add=TRUE, col="red")
#' title("Biological nonsense model")
#' legend(legend=c("model", "pos","neg"),"topright", lwd=c(2,2,2), col=c("black","red","blue"),
#' bty="n" ,lty=c(1,2,2))
#' print(w)
#'
#' ######################
#' ## Plot 2: the data ##
#' ######################
#' n = 100
#' n1 = rbinom(1,n,w[1]); n2 = n - n1
#' y = c(rnorm(n1,mu[1],sd[1]), rnorm(n2,mu[2],sd[2]))
#' miny = floor(min(y))
#' maxy = ceiling(max(y))
#' hist(y, freq=FALSE, breaks=seq(miny, maxy, by=0.5), xlab="MFI", main="Simulated data")
#' curve(w[1]*dnorm(x,mu[1],sd[1])+w[2]*dnorm(x,mu[2],sd[2]),miny,maxy,add=TRUE,col="black", lwd=3)
#' curve(w[1]*dnorm(x,mu[1],sd[1]),miny,maxy,add=TRUE,col="black", lwd=2, lty=2)
#' curve(w[2]*dnorm(x,mu[2],sd[2]),miny,maxy,add=TRUE,col="black", lwd=2, lty=2)
#'
#' #######################
#' ## Unconstrained Fit ##
#' #######################
#' # Estimate parameters of finite mixture model:
#' (fit1 <- em(y,"normal","normal"))
#' hist(y, freq=FALSE, breaks=seq(miny, maxy, by=0.5), xlab="MFI", main="Unconstrained fit")
#' # Add the EM estimated finite mixture model:
#' lines(fit1, col="tomato", lwd=2)
#' # Estimate a cutoff from the fitted mixture model
#' (cut_off <- cutoff(fit1, whose="Titterington", nb=1000))
#' polygon(c(cut_off[-1],rev(cut_off[-1])),c(0,0,.55,.55), col=rgb(1, 0.39, 0.28,.2),border=NA)
#' abline(v=cut_off[-1],lty=2,col="tomato")
#' abline(v=cut_off[1],col="tomato")
#'
#' ###################
#' ## Penalised fit ##
#' ###################
#' # Estimate parameters of finite mixture model:
#' (fit2 <- em(y,"normal","normal", penaltyScale=1E4))
#' # Replot data
#' hist(y, freq=FALSE, breaks=seq(miny, maxy, by=0.5), xlab="MFI", main="Penalised fit")
#' # Add the penalised-EM estimated finite mixture model:
#' lines(fit2, col="red", lwd=2)
#' # Estimate a cutoff from the fitted mixture model
#' (cut_off <- cutoff(fit2, whose="Titterington", nb=1000))
#' polygon(c(cut_off[-1],rev(cut_off[-1])),c(0,0,.55,.55), col=rgb(1,0,0,.2),border=NA)
#' abline(v=cut_off[-1],lty=2,col="red")
#' abline(v=cut_off[1],col="red")
#'
#' ###################################
#' ## Exploring different penalties ##
#' ###################################
#' par(mfrow=c(2,3))
#' fit0 <- em(y,"normal","normal")
#' fit2 <- em(y,"normal","normal", penaltyScale=1E2)
#' fit4 <- em(y,"normal","normal", penaltyScale=1E4)
#' fit6 <- em(y,"normal","normal", penaltyScale=1E6)
#' fit8 <- em(y,"normal","normal", penaltyScale=1E8)
#' fit10 <- em(y,"normal","normal", penaltyScale=1E10)
#' plot(y, fit0$pPositive, main="penaltyScale=0", xlab="Serology data", ylab="p(positive)", ylim=0:1)
#' abline(h=0)
#' plot(y, fit2$pPositive, main="penaltyScale=1E2", xlab="Serology data", ylab="p(positive)", ylim=0:1)
#' abline(h=0)
#' plot(y, fit4$pPositive, main="penaltyScale=1E4", xlab="Serology data", ylab="p(positive)")
#' abline(h=0)
#' plot(y, fit6$pPositive, main="penaltyScale=1E6", xlab="Serology data", ylab="p(positive)")
#' abline(h=0)
#' plot(y, fit8$pPositive, main="penaltyScale=1E8", xlab="Serology data", ylab="p(positive)")
#' abline(h=0)
#' plot(y, fit10$pPositive, main="penaltyScale=1E10", xlab="Serology data", ylab="p(positive)")
#' abline(h=0)
#'
#' ############################################
#' ## Comparing different model combinations ##
#' ############################################
#' \dontrun{
#' yy = y - min(y) + 0.1
#' penScale = 1E6
#' par(mfrow=c(4,4))
#' (fit1  <- em(yy,"normal","normal", penScale, TRUE))
#' hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit1$deviance), main="normal-normal");           lines(fit1, col="blue", lwd=2)
#' (fit2  <- em(yy,"normal","weibull", penScale, TRUE))
#' hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit2$deviance), main="normal-weibull");          lines(fit2, col="blue", lwd=2)
#' (fit3  <- em(yy,"normal","gamma", penScale, TRUE))
#' hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit3$deviance), main="normal-gamma");            lines(fit3, col="blue", lwd=2)
#' (fit4  <- em(yy,"normal","log-normal", penScale, TRUE))
#' hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit4$deviance), main="normal - log-normal");     lines(fit4, col="blue", lwd=2)
#' (fit5  <- em(yy,"weibull","normal", penScale, TRUE))
#' hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit5$deviance), main="weibull-normal");          lines(fit5, col="blue", lwd=2)
#' (fit6  <- em(yy,"weibull","weibull", penScale, TRUE))
#' hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit6$deviance), main="weibull-weibull");         lines(fit6, col="blue", lwd=2)
#' (fit7  <- em(yy,"weibull","gamma", penScale, TRUE))
#' hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit7$deviance), main="weibull-gamma");           lines(fit7, col="blue", lwd=2)
#' (fit8  <- em(yy,"weibull","log-normal", penaltyScale=1E9))
#' hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit8$deviance), main="weibull - log-normal");    lines(fit8, col="blue", lwd=2)
#' (fit9  <- em(yy,"gamma","normal", penScale, TRUE))
#' hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit9$deviance), main="gamma-normal");            lines(fit9, col="blue", lwd=2)
#' (fit10 <- em(yy,"gamma","weibull", penScale, TRUE))
#' hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit10$deviance), main="gamma-weibull");           lines(fit10, col="blue", lwd=2)
#' (fit11 <- em(yy,"gamma","gamma", penScale, TRUE))
#' hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit11$deviance), main="gamma-gamma");             lines(fit11, col="blue", lwd=2)
#' (fit12 <- em(yy,"gamma","log-normal", penScale, TRUE, thresh=exp(-25)))
#' hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit12$deviance), main="gamma-log-normal");        lines(fit12, col="blue", lwd=2)
#' (fit13 <- em(yy,"log-normal","normal", penScale, TRUE))
#' hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit13$deviance), main="log-normal - normal");     lines(fit13, col="blue", lwd=2)
#' (fit14 <- em(yy,"log-normal","weibull", penScale, TRUE))
#' hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit14$deviance), main="log-normal - weibull");    lines(fit14, col="blue", lwd=2)
#' (fit15 <- em(yy,"log-normal","gamma", penScale, TRUE, thresh=exp(-25)))
#' hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit15$deviance), main="log-normal-gamma");        lines(fit15, col="blue", lwd=2)
#' (fit16 <- em(yy,"log-normal","log-normal", penScale, TRUE))
#' hist(yy, freq=FALSE, breaks=seq(minyy, maxyy, by=0.5), xlab="MFI", sub=paste("deviance", fit16$deviance), main="log-normal - log-normal"); lines(fit16, col="blue", lwd=2)
#' }
#'
#' @export
# This function uses the EM algorithm to calculates parameters "lambda"
# (E step), "mu1", "sigma1", "mu2" and "sigma2" (M step).
em <- function(data, D1, D2, penaltyScale=0, forceOrder = FALSE, threshold=1e-64) {
  data_name <- unlist(strsplit(deparse(match.call()),"="))[2]
  data_name <- sub(",.*$","",gsub(" ","",data_name))
  start <- as.list(startval(data,D1,D2))
  d1 <- dHash[[D1]]
  d2 <- dHash[[D2]]
  p1 <- pHash[[D1]]
  p2 <- pHash[[D2]]
  q1 <- qHash[[D1]]
  q2 <- qHash[[D2]]
  lambda0 <- 0 # the previous value of lambda (scalar).
  # counter = 0 # For debugging only
  with(start, {
    # with(start, (abs(lambda0-mean(lambda))>threshold)  )
    while(abs(lambda0-mean(lambda))>threshold) {
      # counter = counter + 1
      # print(counter)
      # browser()
      lambda  <- mean(lambda)
      lambda0 <- lambda
      # Expectation step:
      distr1 <- lambda*d1(data,mu1,sigma1)
      distr2 <- (1-lambda)*d2(data,mu2,sigma2)
      lambda <- distr1/(distr1+distr2) # lambda is a vector.
      # Minimization step (maximum-likelihood parameters estimations):
      mLL2 <- function(mu1,sigma1,mu2,sigma2)
        return(mLL(mu1,sigma1,mu2,sigma2,lambda,data,D1,D2,d1,d2,p1,p2,q1,q2,penaltyScale))
      start <- as.list(trans(D1, mu1, sigma1, D2, mu2, sigma2))
      out   <- bbmle::mle2(mLL2,start,"Nelder-Mead", control=list(maxit=10000))
      # The following lines assign the MLE values to the corresponding parameters:
      # print( out)
      # print( log(abs(lambda0-mean(lambda))) )
      coef <- out@coef
      coef <- backTrans(D1, coef["mu1"], coef["sigma1"], D2, coef["mu2"], coef["sigma2"])
      coef_n <- names(coef)
      # names(coef) <- NULL # No longer needed. Causes problems in print.em
      for(i in 1:4) assign(coef_n[i], (coef[i]))
      #
      mean1 <- getMean(D1, mu1, sigma1)
      mean2 <- getMean(D2, mu2, sigma2)
      if(mean1 < mean2) {
        pPositive = 1 - lambda
      } else {
        pPositive = lambda
        if (forceOrder == TRUE) {
          D1old      <- D1
          D1         <- D2
          D2         <- D1old
          d1         <- dHash[[D1]]
          d2         <- dHash[[D2]]
          p1         <- pHash[[D1]]
          p2         <- pHash[[D2]]
          q1         <- qHash[[D1]]
          q2         <- qHash[[D2]]
          mu1_old    <- mu1
          mu1        <- mu2
          mu2        <- mu1_old
          sigma1_old <-sigma1
          sigma1     <-sigma2
          sigma2     <- sigma1_old
          lambda     <- 1 - lambda
          lambda0    <- 0 ## Forces at least one more iteration
          warning("Note: forceOrder = TRUE. Models 1 & 2 have been switched.")
        }
      }
    }
    #
    out <- list(
      pPositive=pPositive,
      lambda=lambda,
      param=coef, # exp(out@coef), # backTrans(D1, mu1, sigma1, D2, mu2, sigma2)
      D1=D1,
      D2=D2,
      deviance=2*out@min, # Was previously deviance = out@min (i.e. -log-lik), which is a non-standard definition.
      data=data,
      data_name=data_name,
      out=out,
      threshold=threshold)
    class(out) <- "em"
    return(out)
  })
}

#-------------

#' Print method of S3-class "em".
#'
#' @export
#' @method print em
print.em <- function(object) {
  hash <- list(
    normal=c("mean","sd"),
    "log-normal"=c("meanlog","sdlog"),
    gamma=c("shape","rate"),
    weibull=c("shape","scale")
  )
  param <- as.list(object$param)
  digits <- unlist(options("digits"))
  cat(paste0("\n        Finite mixture model fitting to dataset \"", object$data_name,"\":\n\n"))
  cat(paste("distribution 1:",object$D1,"\n"))
  cat(paste0("  probability        = ", round(mean(object$lambda),digits),"\n"))
  cat(paste0("  location (",hash[[object$D1]][1],")   = ", round(param$mu1,digits),"\n"))
  cat(paste0("  scale (",hash[[object$D1]][2],")      = ", round(param$sigma1,digits),"\n"))
  cat(paste0("  mean               = ", round(getMean(object$D1, param$mu1, param$sigma1),digits),"\n"))
  cat(paste0("  standard deviation = ", round(getSD(object$D1, param$mu1, param$sigma1),digits),"\n"))
  cat(paste("distribution 2:",object$D2,"\n"))
  cat(paste0("  probability        = ", 1-round(mean(object$lambda),digits),"\n"))
  cat(paste0("  location (",hash[[object$D2]][1],")   = ", round(param$mu2,digits),"\n"))
  cat(paste0("  scale (",hash[[object$D2]][2],")      = ", round(param$sigma2,digits),"\n"))
  cat(paste0("  mean               = ", round(getMean(object$D2, param$mu2, param$sigma2),digits),"\n"))
  cat(paste0("  standard deviation = ", round(getSD(object$D2, param$mu2, param$sigma2),digits),"\n"))
  cat(paste("deviance of the fitted model:", round(object$deviance,digits),"\n\n"))
}

#-------------

#' Lines method of S3-class "em".
#'
#' @param object An object generated by the em function.
#'
#' @export
#' @method lines em
lines.em <- function(object,...) {
# ...: parameter passed to the "line" function.
  with(object,with(as.list(param), {
    lambda <- mean(lambda)
    curve(lambda*dHash[[D1]](x,mu1,sigma1),add=T,n=512,lty=2,...) # curve(lambda*dnorm(x,mu1,sigma1),add=T,n=512,lty=2,...)
    curve((1-lambda)*dHash[[D2]](x,mu2,sigma2),add=T,n=512,lty=2,...)   # curve((1-lambda)*dnorm(x,mu2,sigma2),add=T,n=512,lty=2,...)
    curve(lambda*dHash[[D1]](x,mu1,sigma1) + (1-lambda)*dHash[[D2]](x,mu2,sigma2),add=T,n=512,...) # curve(lambda*dnorm(x,mu1,sigma1)+ (1-lambda)*dnorm(x,mu2,sigma2),add=T,n=512,...)
  }))
}

#-------------

#' Confint method of S3-class \code{em}.
#'
#' This method of the S3 \code{em} class calculates the confidence intervals of
#' the parameters of the fitted finite mixture model.
#'
#' Confidence intervals of the parameters of probability distributions are
#' calculated by the \code{confint} method of the S4 \code{confint} class of the
#' \code{bbmle} package with default values. See the help of this method for
#' technical details. The confidence interval of the \code{lambda} parameter is
#' calculated thanks to the information-based method of Oakes (1999). The
#' possible non-independance between lambda and the parameters of the probability
#' distributions is accounted for by Monte Carlo simulations where each iteration
#' consists in (i) sampling values of theses parameters in a multinormal
#' distribution, and (ii) applying the method of Oakes (1999). Samplings in the
#' multinormal distribution is performed by the \code{rmultinormal} function of
#' the \code{mc2d} package.
#' @param object an object of class em, produced via the function em.
#' @param nb Number of Monte Carlo simulations.
#' @param level The confidence level required.
#' @inheritParams em
#' @return A dataframe containing point estimates (first column) and confidence
#'   intervals (second and third columns) of each of the five parameters of the
#'   finite mixture model (five row, one per parameter).
#' @references David Oakes (1999) Direct calculation of the information matrix
#'   via the EM algorithm. J R Statist Soc B, 61: 479-482.
#' @export
#' @method confint em
# This function returns the parameter values and their confidence
# intervals from an output of the "em" function.
confint.em <- function(object,threshold=1e-64,nb=10,level=.95) {
  a <- coef_ci(object,level)
  b <- lambda_ci(object,threshold,nb,level)
  return(rbind(a,b))
}
