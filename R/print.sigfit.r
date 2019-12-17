#' @include convergence.r
#' @export
#' 
print.sigfit <- function(x, ...){
  oldx <- x
  METHOD <- ifelse(x$method=="pl", "Pseudo-likelihood", "Nested-pseudo-likelihood")
  
  cat("\nEstimated parameters of the Crisis Siginaling Model\n\nCALL:\n\n")
  print(x$call)
  cat("\nMethod:\n")
  cat(METHOD, "\n")
  cat("\nCOEFFICIENTS:\n")
  par <- x$coefficients
  mf.new <- x$model
  regr <- list()
  for(i in 1:7){
    regr[[i]] <- model.matrix(x$formulas, data=mf.new, rhs=i)
  }
  names(regr) <- c("SA", "VA", "CB", "barWA", "barWB", "bara", "VB")
  ncols <- sum(sapply(regr, ncol))
  ngames <- nrow(regr[[1]])
  u.names <- rep(names(regr), sapply(regr, ncol))
  
  
  ### NOTE this function --- and everything else --- is only for the 
  ### basic Lewis and Schultz (2003) model.  Does not directly extend to
  ### the extensions that estimate the covariances.  Predicated on 7 sets 
  ### of parameters
  
  idx0 <- lapply(regr, ncol)
  idx0 <- sapply(idx0, function(x){if(is.null(x)){0}else{x}})
  idx1 <- cumsum(idx0)
  idx0 <- idx1-idx0+1
  idx <- rbind(idx0, idx1)
  idx[,apply(idx, 2, function(x){x[1]>x[2]})] <- 0
  idx[,apply(idx, 2, function(x){x[1]==x[2]})] <- rbind(0,idx[1,apply(idx, 2, function(x){x[1]==x[2]})] )
  
  indx <- list(idx[1,1]:idx[2,1],
               idx[1,2]:idx[2,2],
               idx[1,3]:idx[2,3],
               idx[1,4]:idx[2,4],
               idx[1,5]:idx[2,5],
               idx[1,6]:idx[2,6],
               idx[1,7]:idx[2,7])
  indx <- lapply(indx, function(x){if(0 %in% x){return(x[length(x)])}else{return(x)}})
  
  
  ## print a table of coefficients for each utility 
  for (i in 1:length(regr)) {
    par.i <- par[indx[[i]]]
    u.names <- names(regr)[[i]]
    varnames <- colnames(regr[[i]])
    outcome <- switch(u.names,
                      SA = "SQ",
                      VA = "CD",
                      CB = "CD",
                      barWA = "SF",
                      barWB = "SF",
                      bara = "BD",
                      VB = "BD")
    Player <- switch(u.names,
                     SA = "A",
                     VA = "A",
                     CB = "B",
                     barWA = "A",
                     barWB = "B",
                     bara = "A",
                     VB = "B")
    
    if(Player == "A"){
      cat("\n--------------\n")
      cat("OUTCOME: ", outcome, "\n", sep="")
    }
    cat("Player ", Player, "'s utility \n", sep="")
    names(par.i) <- varnames
    tab <- data.frame(as.matrix(par.i))
    names(tab) <- " "
    if(nrow(tab)>0){
      print(tab)
      cat("\n")
    }else{
      cat("Fixed to 0\n")
    }
  }
  
  cc <- convergenceCriterion(x$maxlik.method)
  if (!(x$maxlik.code %in% cc)) {
    cat("\nWarning: Model fitting did not converge\nCode:",
        x$maxlik.code, "\nMessage:", x$maxlik.message)
  }
}