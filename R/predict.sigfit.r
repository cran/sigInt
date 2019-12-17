#' Predicted probabilities and comparative statics for signaling games
#' 
#'
#' This method uses a fitted model of class \code{sigfit} to compute
#' predicted probabilities or comparative statics. 
#' Users can provide either new data or a new parameters to generate 
#' counterfactuals of interest.
#' 
#' 
#' @param object a fitted model of class \code{sigfit}.
#' @param newdata data frame of covariates used to produce the predicted probabilities.
#' If this is left empty, the entire original data set is used.
#' When \code{newdata} is specified, \code{new.theta} should be either missing (use the coefficients
#' from the \code{sigfit} object) or a one row data frame. See "Details" for more information.
#' As with other \code{\link[stats]{predict}} methods, variable names must match those used to fit the model.
#' @param new.theta a data frame of alternative parameters for comparative statics.
#' When missing, the coefficients from the \code{object} are used.
#' When specified, each row should be a complete parameter vector. 
#' If \code{new.theta} is specified, then \code{newdata} must be a data frame with only one row.
#' Unlike other \code{\link[stats]{predict}} methods, column names do not matter here.
#' Instead, the columns must be the same order as the coefficients in \code{object}.
#' See "Details" and "Examples" for more information.
#' @param type whether to provide probabilities over actions 
#' (default, returns \eqn{p_C}, \eqn{p_R}, and \eqn{p_F}) 
#' or outcomes (returns \eqn{SQ}, \eqn{CD}, \eqn{SF}, and  \eqn{BD}).
#' @param na.action how to deal with \code{NA}s in \code{newdata}.
#' @param control list of options describing the grid search method. See "Details" for more information
#' @param parallel logical. Should the comparative statics be computed in parallel, requires the
#' \code{\link[parallel]{parallel}}
#' package be installed. Parallelization is done using \code{\link[parallel]{parSapply}}.
#' @param ... Additional arguments (not currently used)
#' @return An object of class \code{sigProb} containing three elements:
#' \describe{
#'  \item{\code{predicted}}{data frame of predicted probabilities. The first column of this data frame is 
#'  called \code{Row}, which corresponds to the rows in either \code{model} or \code{par}.
#'  In the event of multiple equilibria, this column allows for mapping data and parameters to 
#'  all computed equilibria.
#'  }
#'  \item{\code{model}}{data frame of covariates used to produce the predicted probabilities.}
#'  \item{\code{par}}{data frame of parameters used to produce the predicted probabilities.}
#'  }
#'
#' @details This function is used to consider comparative statics in the crisis signaling game.
#' The model of interest is fit using \code{\link{sigint}}.
#' How this function behaves largely depends on  how \code{newdata} and \code{new.theta} are
#' specified.
#'
#' When both \code{newdata} and \code{new.theta} are missing, all equilibria for every
#' observation used to fit the model are computed.
#' These equilibria are then used to calculated either choice probabilities (default,
#' \code{type = "action"}) the distribution over outcomes (\code{type = "outcomes"}).
#'
#' When only \code{newdata} is specified, then all equilibria are computed using the
#' data frame in \code{newdata} and the coefficients from \code{object}. This produces
#' standard comparative statics with respect to observed covariates.
#'
#' When \code{newdata} is specified and \code{new.theta} is a one row data frame,
#' then all equilibria are computed using the data frame in \code{newdata} and the
#' coefficients from \code{new.theta}.
#'
#' When \code{newdata} is a one row data frame and \code{new.theta} is specified,
#' then all equilibria are computed using the data frame in \code{newdata} and the
#' coefficients from \code{new.theta}. This is a comparative static on changing a
#' structural parameter in the model.
#'
#' If \code{new.theta} has more than one row, then \code{newdata} must be specified as a
#' data frame with only one row.  Anything else returns an error.
#'
#' Equilibria are computed using a line search method.
#' The \code{control} argument allows for user control over this process.
#' Users can specify a list with the following named elements
#' \describe{
#'  \item{gridsize}{Integer. The number of points considered in the line search (default, 1e4).
#'  More points makes it more likely that all equilibria are discovered, but can slow down the
#'  search. }
#'  \item{comp}{Logical. Should an equilibrium be computed when discovered?
#'             When \code{comp = FALSE} (default), the mean of the grid points surrounding the
#'             equilibrium is used as an approximate solution.
#'             When \code{comp = TRUE} the function \code{\link[stats]{uniroot}}
#'             is called to find a more precise solution to the equilibrium constraint problem.}
#'  \item{tol}{Numeric. When \code{comp = TRUE}, this is the tolerance used by
#'            \code{\link[stats]{uniroot}}}.
#' }
#'
#' When dealing with a larger problem, such as computing all equilibria for every observation,
#' it can be helpful to parallelize process.  If the user has the
#' (suggested) \code{\link[parallel]{parallel}} package, then the option \code{parallel = TRUE}
#' will use the function \code{\link[parallel]{parSapply}} is used.
#' @seealso \code{\link{plot.sigProb}}, \code{\link{generate.eq}}
#' @import Formula
#' @import stats
#' @import utils
#' @import stringr
#' @export
#'
#' @examples
#' data(sanctionsData)
#'
#' f1 <- sq+cd+sf+bd ~ sqrt(senderecondep) + senderdemocracy + contig + ally -1|#SA
#'                     anticipatedsendercosts|#VA
#'                     sqrt(targetecondep) + anticipatedtargetcosts + contig + ally|#CB
#'                     sqrt(senderecondep) + senderdemocracy + lncaprat | #barWA
#'                     targetdemocracy + lncaprat| #barWB
#'                     senderdemocracy| #bara
#'                     -1#VB
#
#' ## Using Nested-Pseudo Likelihood  with default first stage    
#' \dontrun{             
#' fit1 <- sigint(f1, data=sanctionsData, npl.trace=TRUE)
#' p.out <- predict(fit2, parallel=TRUE) #fitted choice probabilites for all observations
#' }
#'
#' ## Outcome probabilities for first five using PL method
#' Phat <- list(PRhat=sanctionsData$PRhat, PFhat=sanctionsData$PFhat)
#' fit2 <- sigint(f1, data=sanctionsData, method="pl", phat=Phat)
#' p1 <- predict(fit2, newdata=sanctionsData[1:5,], type="outcome")
#'
#' ## comparative static on \bar{a}, compute more precise equilibria with uniroot
#' new.theta <- data.frame(t(replicate(25, coef(fit2))))
#' new.theta[,19] <- seq(-6, 0, length=25)
#' p2 <- predict(fit2, newdata=sanctionsData[1,], new.theta=new.theta, control=list(comp=TRUE))
#'
#' 

predict.sigfit <- function(object, newdata, new.theta, type=c("actions", "outcomes"),
                           na.action=na.pass, control=list(), parallel=FALSE,...){
  
  control.default <- list(gridsize=1e4, comp=F, tol=1e-8)
  control <- modifyList(control.default, control)
  
  if(!missing(new.theta)){  #check for new.theta
    if(nrow(unique(new.theta)) != 1){ #is something in theta changing?
      if(missing(newdata) || is.null(newdata)){  #if theta is changing, there needs to be 1 newdata
        warning("When taking a comparative static on a parameter, 1 row of new data is required. Setting all covariates to zero.")
        newdata <- matrix(0, nrow=1, ncol=ncol(object$model))
        colnames(newdata) <- colnames(object$model)
        newdata <- as.data.frame(newdata)
      }else{
        if(nrow(unique(newdata)) != 1){
          stop("Varying both newdata and new.theta is not supported.")
        }else{
          newdata <- unique(newdata)
        }
      }
    }else{ #if theta is actually unique
      new.theta <- unique(new.theta)
    }
  }
  
  type <- match.arg(type)
  
  if (missing(newdata) || is.null(newdata)) {
    ## use original data if 'newdata' not supplied
    mf.new <- object$model
  } else {
    ## drop LHS
    formulas <- Formula(delete.response(terms(formula(object$formulas))))
    
    mf.new <- model.frame(formulas, data = newdata, na.action = na.action,
                          xlev = object$xlevels)
    
    ## check that variables are of the right classes
    Terms <- attr(object$model, "terms")
    if (!is.null(cl <- attr(Terms, "dataClasses"))){
      .checkMFClasses(cl, mf.new)
    }
  }
  
  
  #make the independent variables
  regr <- list()
  for(i in 1:7){
    regr[[i]] <- model.matrix(object$formulas, data=mf.new, rhs=i)
  }
  names(regr) <- c("SA", "VA", "CB", "barWA", "barWB", "bara", "VB")
  ncols <- sum(sapply(regr, ncol))
  ngames <- nrow(regr[[1]])
  u.names <- rep(names(regr), sapply(regr, ncol))
  
  
  if(missing(new.theta) || is.null(new.theta)){
    par <- as.data.frame(matrix(object$coefficients, nrow=1))
    npar <- 1
  }else{
    par <- new.theta
    npar <- nrow(new.theta)
  }
  colnames(par) <- names(object$coefficients)
  if(ncol(par) != length(object$coefficients)){
    stop("Differing number of parameters in new.theta and the  original model")
  }
  
  
  Ulist <- apply(as.matrix(par), 1, vec2U.regr, regr=regr, fixed.par=object$fixed.par)
  Ulist <- unlist(Ulist, recursive = F)
  if(npar>1){#Situation where we vary theta, so we can clean this up
    Ulist <- sapply(unique(names(Ulist)), 
                    function(x) unname(unlist(Ulist[names(Ulist)==x])), 
                    simplify=FALSE)
    # names(Ulist)   <- c(names(regr), "sig") #wrong order
    U.names <- do.call(rbind, str_split(names(Ulist), "\\."))
    U.names <- U.names[,ncol(U.names)]
    names(Ulist) <- U.names
  }else{
    # names(Ulist)   <- c(names(regr), "sig") #wrong order
    U.names <- do.call(rbind, str_split(names(Ulist), "\\."))
    U.names <- U.names[,ncol(U.names)]
    names(Ulist) <- U.names
    Ulist$sig <- rep(1, length(Ulist$SA))
  }
  
  
  out <- list()
  grid <- seq(from=0-.Machine$double.eps, to=1+.Machine$double.eps, 
              length.out=control$gridsize) #Actually include 0 and 1?
  length.out <- max(sapply(Ulist, length))
  
  # for(i in 1:length.out){
  
  compstat <- function(i){
    Ui <- lapply(Ulist, function(x){return(x[i])})
    fgrid <- const.jo(grid,Ui)
    sols <- which(tail(fgrid,-1)*head(fgrid,-1) <= 0)
    sols <- matrix(grid[c(sols, sols+1)], nrow=length(sols))
    
    # compute equilibria
    if (!control$comp){
      sols <- rowMeans(sols)
    } else{
      solver <- function(x){uniroot(function(x){const.jo(x,Ui)}, x, tol=control$tol)$root}     	
      sols <- apply(sols,1,solver)    
    }
    
    #return parameters of interest
    if (!length(sols)){
      return(cbind(i,NaN,NaN,NaN))
    } else {
      
      sols.pc <- g.jo(cStar.jo(sols,Ui),Ui)
      sols.pc[sols.pc <= .Machine$double.eps] <- .Machine$double.eps
      sols.pc[sols.pc >= 1- .Machine$double.eps] <- 1- .Machine$double.eps
      sols.pf <- h.jo(cStar.jo(sols,Ui),Ui)/sols.pc
      sols.pf[sols.pf <= .Machine$double.eps] <- .Machine$double.eps
      sols.pf[sols.pf >= 1- .Machine$double.eps] <- 1- .Machine$double.eps
      index <- is.nan(sols.pc*sols*sols.pf)
      sols.onset <- c(sols.pc[index], (sols.pc*sols*sols.pf)[!index])
      
      return(cbind(i,
                   sols,
                   sols.pc,
                   sols.onset,
                   sols.pf))
      
    }
  }
  
  if(parallel){
    if(!requireNamespace("parallel", quietly = TRUE)){
      warning("parallel option specified, but parallel package not found. 
Please install parallel to use this option.  Switching to parallel=FALSE")
      parallel <- FALSE
    }
  }
  if(parallel){
    clust <- parallel::makeCluster(parallel::detectCores())
    parallel::clusterExport(cl=clust,
                            varlist = c("Ulist",
                                        "const.jo",
                                        "cStar.jo",
                                        "g.jo",
                                        "h.jo",
                                        "f.jo",
                                        "grid",
                                        "control"),
                            envir = environment())
    parallel::clusterEvalQ(cl=clust, expr=library(pbivnorm))
    map <- function(X, FUN){
      parallel::parSapply(cl=clust, 
                          X, FUN, simplify=F)
    }
    on.exit(parallel::stopCluster(cl=clust))
  }else{
    map <- function(X, FUN){
      sapply(X, FUN, simplify=FALSE)
    }
  }
  out <- map(X=1:length.out, FUN=compstat)
  out <- do.call(rbind.data.frame, out)
  
  
  if(type=="actions"){
    output <- with(out, data.frame(Row=i,
                                   pc = sols.pc,
                                   pr = sols,
                                   pf = sols.pf))
  }else{
    output <- with(out, data.frame(Row=i,
                                   SQ = 1-sols.pc,
                                   CD = sols.pc*(1-sols),
                                   SF = sols.pc*sols*sols.pf,
                                   BD = sols.pc*sols*(1-sols.pf)))
    
  }
  colnames(mf.new) <- stringr::str_replace(string=colnames(mf.new),
                                           pattern=":", replacement=".")
  colnames(mf.new) <- stringr::str_replace(string=colnames(mf.new),
                                           pattern="\\(", replacement="")
  colnames(mf.new) <- stringr::str_replace(string=colnames(mf.new),
                                           pattern="\\)", replacement="")
  colnames(mf.new) <- stringr::str_replace(string=colnames(mf.new),
                                           pattern=" ", replacement="")    
  
  mf.new$Row <- 1:nrow(mf.new)
  
  colnames(par) <- stringr::str_replace(string=names(object$coef),
                                        pattern=":", replacement=".")
  colnames(par) <- stringr::str_replace(string=colnames(par),
                                        pattern="\\(", replacement="")
  colnames(par) <- stringr::str_replace(string=colnames(par),
                                        pattern="\\)", replacement="")
  colnames(par) <- stringr::str_replace(string=colnames(par),
                                        pattern=" ", replacement="")    
  par$Row <- 1:npar    
  output <- list(predicted = output,
                 model = mf.new,
                 par = par)
  if(any(table(output$predicted$Row))==2){
    warning("Only two equilibria found under these settings, consider a larger grid")
  }
  class(output) <- c("sigProb")
  return(output)
}
