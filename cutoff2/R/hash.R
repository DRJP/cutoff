#' @noRd
dHash <- list(normal       = stats::dnorm,
              `log-normal` = stats::dlnorm,
              gamma        = stats::dgamma,
              weibull      = stats::dweibull)

#' @noRd
pHash <- list(normal       = stats::pnorm,
              `log-normal` = stats::plnorm,
              gamma        = stats::pgamma,
              weibull      = stats::pweibull)

#' @noRd
qHash <- list(normal       = stats::qnorm,
              `log-normal` = stats::qlnorm,
              gamma        = stats::qgamma,
              weibull      = stats::qweibull)

#' @noRd
transOne <- function(D, mu, sigma) {
  if(D=="log-normal") {
    return( c(mu, log(sigma)) ) # first parameter of dlnorm does not require transforming to (-Inf, Inf)
  } else {
    return( log(c(mu, sigma)) )
  }
}

#' @noRd
trans <- function(D1, mu1, sigma1, D2, mu2, sigma2) {
  #  browser()
  out <- c(transOne(D1, mu1, sigma1),
           transOne(D2, mu2, sigma2))
  names(out) = c("mu1", "sigma1", "mu2", "sigma2")
  return(out )
}

#' @noRd
backTransOne <- function(D, mu, sigma) {
  if(D=="log-normal") {
    return( c(mu, exp(sigma)) ) # first parameter of dlnorm does not require transforming to (-Inf, Inf)
  } else {
    return( exp(c(mu, sigma)) )
  }
}

#' @noRd
backTrans <- function(D1, mu1, sigma1, D2, mu2, sigma2) {
  # browser()
  out <- c(backTransOne(D1, mu1, sigma1),
           backTransOne(D2, mu2, sigma2))
  names(out) = c("mu1", "sigma1", "mu2", "sigma2")
  return(out )
}

#' @noRd
getMean <- function(D, mu, sigma) {
  # browser()
  if (D == "normal")
    Mean = mu
  if (D == "log-normal")
    Mean = exp(mu + sigma^2/2)
  if (D == "weibull")
    Mean = sigma*gamma(1+1/mu) # mu=proxy4(shape), sigma=proxy4(scale)
  if (D == "gamma")
    Mean = mu/sigma # mu=proxy4(shape), sigma=proxy4(rate)
  return(Mean)
}
