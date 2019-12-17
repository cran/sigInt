#' Equilibrium analysis of the empirical crisis signaling game
#' 
#'
#' This method uses a formula and fixed data/parameters to allow for analysis 
#' of the crisis signaling game under specific settings.
#' This function is very similar to \code{predict.sigfit}, but it is designed 
#' for analysis outside of conducting counterfactuals on a fitted model.
#' 
#' @param formulas a Formula object with no left-hand side and 
#'   seven separate (7) right-hand sides. See "Details" and examples below.
#' @param data a data frame containing the variables in the model
#'   Each row of the data frame describes an individual game
#'   \eqn{d = 1, 2, ..., D}. Each row \eqn{d}  should be a summary of all of the
#'   within-game observations for game \eqn{d}.  If the model is all constants, 
#'   then this argument should be left empty.
#' @param theta a data frame with one or more rows where each row is a parameter vector. 
#' @param type whether to provide probabilities over actions 
#' (default, returns \eqn{p_C}, \eqn{p_R}, and \eqn{p_F}) 
#' or outcomes (returns \eqn{SQ}, \eqn{CD}, \eqn{SF}, and  \eqn{BD}).
#' @param na.action how to deal with \code{NA}s in \code{data}.
#' @param control list of options describing the grid search method. See "Details" for more information
#' @param parallel logical. Should the comparative statics be computed in parallel, requires the
#' \code{\link[parallel]{parallel}}
#' package be installed. Parallelization is done using \code{\link[parallel]{parSapply}}.
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
#' @details 
#' This function is used to consider comparative statics in the crisis signaling game, where the model
#' of interest has pre-defined parameters.
#' As such, it requires, at minimum, a seven-part formula and parameters.
#' How this function behaves has to do with how \code{data} and \code{theta} are specified.
#' 
#' When the model is all constants (every part of the \code{formula} argument is either \code{0} or \code{1}),
#' then \code{data} is ignored.
#' In these cases, equilibria are computed for every parameter vector, which are supplied 
#' as rows in a data frame to \code{theta}.
#' 
#' When there is one or more covariate in the model, then a data frame must be supplied to \code{data}.
#' In these cases both \code{data} or \code{theta} must have at least one row.
#' However, only one of these arguments can have multiple rows.  In other words, only
#' \code{data} or \code{theta} may vary, but not both.
#' 
#' For additional implementation details see \code{\link{predict.sigfit}}.
#' 
#' 
#' 
#' @examples 
#' ## An example with one covariate
#' ftest1 <-  ~ 0 | #SA
#'              1 | #VA
#'              0 | #CB
#'              1 | #barWA
#'              x1 | #barWB
#'              1| #bara
#'              1 #VB
#'              
#' theta <- data.frame(VA = 1, barWA = -1.9, barWB = -2.9,
#'                     barWB1 = 0.1, bara = -1.2, VB = 1)
#' data <- data.frame(x1 = seq(from = -1,to = 2, length.out = 101))
#' test <- generate.eq(ftest1, data = data, theta = theta)
#' plot(test, prob = "pr")
#' 
#' ## An example with all constants
#' ftest2 <-  ~ 0 | #SA
#'              1 | #VA
#'              0 | #CB
#'              1 | #barWA
#'              1 | #barWB  
#'              1 | #bara
#'              1 #VB
#'              
#' theta <- data.frame(VA = 1, barWA = -1.9, 
#'                     barWB = seq(-2.9, -2.2, length.out = 15),
#'                     bara = -1.2,
#'                     VB = 1)
#' test <- generate.eq(ftest2, theta = theta)
#' plot(test, prob = "pr")
#'
#' @seealso \code{\link{plot.sigProb}}, \code{\link{predict.sigfit}}
#' @import Formula
#' 
#'
#' @export
#'
generate.eq <- function(formulas, data, theta, type=c("actions", "outcomes"),
                        na.action=na.omit, control=list(), parallel=FALSE){
  
  
  control.default <- list(gridsize=1e4, comp=F, tol=1e-8)
  control <- modifyList(control.default, control)
  
  cl <- match.call()
  type <- match.arg(type)
  
  
  ####Formulas####
  ##process main formula
  formulas <- as.Formula(formulas)
  if (length(formulas)[2] != 7){
    stop("'formulas' should have seven components on the right-hand side")
  }
  
  
  
  ## drop LHS
  formulas <- Formula(delete.response(terms(formula(formulas))))
  
  ####Data####
  ## make the model frame
  mf <- match(c("data", "na.action"), names(cl), 0L)
  mf <- cl[c(1L, mf)]
  mf$formula <- formulas
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  #make the independent variables
  regr <- list()
  for(i in 1:7){
    regr[[i]] <- model.matrix(formulas, data=mf, rhs=i)
  }
  names(regr) <- c("SA", "VA", "CB", "barWA", "barWB", "bara", "VB")
  ncols <- sum(sapply(regr, ncol))
  ngames <- nrow(regr[[1]])
  u.names <- rep(names(regr), sapply(regr, ncol))
  regr.names <-  unlist(lapply(regr,colnames))
  regr.names <- paste(u.names, ":", regr.names, sep="")
  if(all( unlist(lapply(regr,colnames)) == "(Intercept)")){
    regr <- lapply(regr, function(x){
      if(ncol(x)==1 & nrow(x)==0){x <- rbind(x,1)}
      if(ncol(x)==0 & nrow(x)>=0){x <- matrix(nrow=1, ncol=0)}
      if(ncol(x)>0  & nrow(x)>0){x <- unique(x)}
      return(x)
    })
    mf <- data.frame()
  }else{
    
    
    if(missing(data) || missing(theta)){
      stop("data and theta are both required")
    }
    if((nrow(unique(theta)) != 1) && (nrow(unique(data)) != 1)){
      stop("Varying both data and theta is not supported.")
    }else{
      data <- unique(data)
      theta <- unique(theta)
    }
  }
  par <- theta; row.names(par) <- NULL
  npar <- nrow(theta)
  
  
  
  Ulist <- apply(as.matrix(par), 1, vec2U.regr, regr=regr,  fixed.par=list())
  Ulist <- unlist(Ulist, recursive = F)
  if(npar>1){#Situation where we vary theta, so we can clean this up
    Ulist <- sapply(unique(names(Ulist)), 
                    function(x) unname(unlist(Ulist[names(Ulist)==x])), 
                    simplify=FALSE)
    # names(Ulist)   <- c(names(regr), "sig")
  }else{
    # names(Ulist)   <- c(names(regr), "sig")
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
  colnames(mf) <- stringr::str_replace(string=colnames(mf),
                                       pattern=":", replacement=".")
  colnames(mf) <- stringr::str_replace(string=colnames(mf),
                                       pattern="\\(", replacement="")
  colnames(mf) <- stringr::str_replace(string=colnames(mf),
                                       pattern="\\)", replacement="")
  colnames(mf) <- stringr::str_replace(string=colnames(mf),
                                       pattern=" ", replacement="")    
  
  if(nrow(mf)>0){mf$Row <- 1:nrow(mf)}
  
  colnames(par) <- stringr::str_replace(string=regr.names,
                                        pattern=":", replacement=".")
  colnames(par) <- stringr::str_replace(string=colnames(par),
                                        pattern="\\(", replacement="")
  colnames(par) <- stringr::str_replace(string=colnames(par),
                                        pattern="\\)", replacement="")
  colnames(par) <- stringr::str_replace(string=colnames(par),
                                        pattern=" ", replacement="")    
  par$Row <- 1:npar    
  output <- list(predicted = output,
                 model = mf,
                 par = par)
  if(any(table(output$predicted$Row))==2){
    warning("Only two equilibria found under these settings, consider a larger grid")
  }
  class(output) <- c("sigProb")
  return(output) 
}



# 
# 
# selectEq <- function(X){
#   M <- length(X)
#   # select equilibria
#   Pstar <- c(0.2964518, 0.4715766, 0.8740314) # these are the only eq
#   Pstar <- rowSums(matrix(rep(Pstar, each=M),nrow=M) * cbind(X<1/3,X<2/3 & X>1/3, X>2/3))
#   
#   # compute equilibra 
#   f <- function(p){const.jo(p,U)}
#   grf <- function(p){diag(1-eval_gr_fh(p,U))}
#   out <- multiroot(f,Pstar, jacfunc=grf, jactype="fullusr", ctol=1e-10,rtol=1e-10,atol=1e-10)
#   Pstar <- out$root	
#   
#   return(Pstar)
# }
# 
# # this generate the datums
# genData.jo <- function(nObs, Pstar, U){
#   # Pstar is a vector of length K, where K is the number of games
#   # nObs is the number of observations for each game
#   
#   M <- length(Pstar)
#   EQ <- eqProbs(Pstar,U)
#   
#   Probs <- cbind(1-EQ[,2], EQ[,2]*(1-EQ[,1]), 	
#                  EQ[,2]*EQ[,1]*EQ[,3], EQ[,2]*EQ[,1]*(1-EQ[,3]))
#   
#   Data <- rmultinomial(M,nObs,Probs)
#   
#   return(t(Data))
# }