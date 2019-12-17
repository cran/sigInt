#' Summarize a fitted crisis signaling model
#' 
#' The default method for summarizing a \code{sigfit} object.
#'
#' Forms a block regression results table from fitted crisis signaling model.
#' 
#' @param object a fitted model of class \code{sigfit}
#' @param vcov a substitute variance covariance matrix 
#' @param ... Additional arguments (not currently used)
#' @return An object of class \code{summary.sigfit}.  This object
#' contains the information needed to print the summary
#' @seealso \code{\link{print.summary.sigfit}}.
#' @include convergence.r
#' @export
summary.sigfit <- function(object, vcov,...)
{
  if(object$method=="pl" && is.null(object$vcov) && missing(vcov)){
    cat("\nPsuedo-likelihood method was used with pl.vcov=FALSE and no user vcov is supplied\nOnly returning point estimates\n")
    
  }
  if(missing(vcov)){
    vcov <- object$vcov
    if(is.null(vcov)){
      vcov <- diag(rep(NA,length(object$coefficients)))
    }
  }
  ## make regression table
  par <- object$coefficients
  se <- sqrt(diag(vcov))
  z.stat <- par/se
  p.value <- 2 * pnorm(abs(z.stat), lower.tail = F)
  
  table.out <- cbind(par, se, z.stat, p.value)
  colnames(table.out) <- c("Estimate", "Std. Error", "z value",
                                  "Pr(>|z|)")
  output <- list(coef.table = table.out,
                 call = object$call,
                 maxlik.method = object$maxlik.method,
                 maxlik.code = object$maxlik.code,
                 maxlik.message = object$maxlik.message,
                 logLik = object$logLik,
                 ngames = length(object$Phat$PRhat),
                 method = object$method)
  output$eq.constraint <- object$eq.constraint
 
  
  class(output) <- "summary.sigfit"
  return(output)
}

#' Print the summary table for a \code{sigfit} object.
#' 
#' Prints the summary regression table for a model fitted with \code{\link{sigint}}.
#'
#' Prints the standard regression results table from a fitted strategic model,
#' along with the log-likelihood and number of games used in estimation.
#' @param x a \code{summary.sigfit} object.
#' @param ... Additional arguments (not currently used)
#' @include convergence.r
#' @export
print.summary.sigfit <- function(x, ...)
{
  METHOD <- ifelse(x$method=="pl", "Pseudo-likelihood", "Nested-psuedo-likelihood")
  
  cat("\nCall:\n")
  print(x$call)
  cat("\nMethod:\n")
  cat(METHOD, "\n")
  cat("\nCoefficients:\n")
  printCoefmat(x$coef.table)



  cat("\nLog-likelihood:", x$logLik)
  cat("\nAIC:", AIC(x))
  cat("\nMax. Equilibrium Constraint Violation:", max(abs(x$eq.constraint)))
 
  cat("\nNumber of Games:", x$ngames, "\n\n")
  cc <- convergenceCriterion(x$maxlik.method)
  if (!(x$maxlik.code %in% cc)) {
    cat("\nWarning: Model fitting did not converge\nCode:",
        x$maxlik.code, "\nMessage:", x$maxlik.message, "\n")
  }
  invisible(x)
}


#' @export
coef.sigfit <- function(object,...){
  return(object$coefficients)
}

#' @export
vcov.sigfit <- function(object,...){
  return(object$vcov)
}

#' @export
logLik.sigfit <- function(object,...){
  out <- object$logLik
  attr(out, "df") <- length(object$coefficients)
  attr(out, "nobs") <- length(object$Phat$PRhat)
  class(out) <- "logLik"
  return(out)
}

#' @export
logLik.summary.sigfit <- function(object,...){
  out <- object$logLik
  attr(out, "df") <- nrow(object$coefficients)
  attr(out, "nobs") <- object$ngames
  class(out) <- "logLik"
  return(out)
}

#' @export
AIC.sigfit <- function(object,...){
  k <- length(object$coefficients)
  L <- object$logLik
  return(-2*k - 2*L)
}


#' @export
AIC.summary.sigfit <- function(object,...){
  k <- nrow(object$coef.table)
  L <- object$logLik
  return(-2*k - 2*L)
}
