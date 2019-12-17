#' Estimating the parameters of the canonical discrete crisis bargaining game.
#' 
#' This function fits the Lewis and Schultz (2003) model to data using either 
#' the pseudo-likelihood (PL) or nested-pseudo likelihood (NPL) method from 
#' Crisman-Cox and Gibilisco (2018). Throughout, we refer to the data as
#' containing \eqn{D} games, where each game is  observed one or more times.
#' 
#' @param formulas a \code{Formula} object four variables on the left-hand side and 
#'   seven (7) separate right-hand sides. See "Details" and examples below.
#' @param data a data frame containing the variables used to fit the model.
#'   Each row of the data frame describes an individual game
#'   \eqn{d = 1, 2, ..., D}. Each row \eqn{d}  should be a summary of all of the
#'   within-game observations for game \eqn{d}. See "Details" for more
#'   information.
#' @param subset an optional logical expression to specify a subset of 
#'   observations to be used in fitting the model.
#' @param na.action how do deal with missing data (\code{NA}s).  Defaults to the 
#'   \code{na.action} setting of \code{\link[base]{options}} (typically \code{na.omit}).
#' @param fixed.par a list with up to seven (7) named elements for normalizing payoffs to non-zero values.  
#'    Names must match a payoff name as listed in "Details."
#'    Each named element should contain a single number that is the fixed (not estimated) value of that payoff.
#'    For example, to fix each side's victory-without-fighting payoff to 1
#'     use \code{fixed.par=list(VA=1, VB=1)} and set their portions of the \code{formulas} to zero.
#'    To normalize a payoff to zero, you only need to specify it has a zero in the \code{formulas}.
#' @param method whether to use the nested-pseudo-likelihood (\code{"npl"}, default) or
#'   the pseudo-likelihood method for fitting the model. See "Details" for more 
#'   information.
#' @param npl.maxit maximum number of outer-loop iterations to be used when fitting the NPL.
#'   See "Details" for more information.
#' @param npl.tol Convergence criteria for the NPL. When the estimates change by
#'   less than this amount, convergence is considered successful.
#' @param npl.trace logical. Should the NPL's progress be printed to screen?
#' @param start.beta starting values for the model coefficients as a single 
#'   vector. If missing, random values are drawn from a normal distribution with mean
#'   zero and standard deviation 0.05.
#' @param maxlik.method method used by  \code{\link[maxLik]{maxLik}}  to fit the
#'   model. Default is Newton-Raphson (\code{"NR"}). See  \code{\link[maxLik]{maxLik}} 
#'   for additional details. At this time only \code{"NR"}, \code{"BFGS"}, and \code{"Nelder-Mead"} are available.
#' @param phat a list containing two vectors: \code{PRhat} and \code{PFhat}.
#'   These are the first-stage estimates that \eqn{B} resists a threat and that \eqn{A} 
#'   follows through on a threat, respectively.
#'   If missing, they will be estimated by a  \code{\link[randomForest]{randomForest}}  with default options.
#'   See "Details" for more information.
#' @param phat.formulas if \code{phat} are  missing, you can supply formulas to 
#'   estimate them. Should be a Formulas object containing no left-hand side and
#'   1-2 right hand sides.  If  one right-hand side is given, the same 
#'   covariates are used to estimate both \code{PRhat} and \code{PFhat}.  
#'   Otherwise, the first RHS is used to generate \code{PRhat}, while the second RHS 
#'   generates \code{PFhat}.
#'   If no formulas are provided and \code{phat} is 
#'   missing, all the covariates used in formulas argument and used here. See 
#'   "Details" for more information.
#' @param  pl.vcov number of bootstrap iterations to generate \code{phat.vcov}. 
#'   If less than \code{0} or \code{FALSE} (default), the pseudo-likelihood 
#'   covariance is not estimated.  Only used if \code{method = "pl"}.
#' @param phat.vcov a covariance matrix for the estimates \code{PRhat} and \code{PFhat}.  
#'   If missing and \code{pl.vcov = TRUE} and \code{phat} is missing, it will be 
#'   estimated by bootstrapping the random forest used to fit \code{phat}.
#' @param seed integer.
#'   Used to set the seed for the random forest and for drawing the the starting values.
#'   The PL can be sensitive to starting value, so this makes results reproducible.
#'   The NPL is less sensitive, but we always recommend checking the first order conditions.
#' @param maxlik.options a list of options to be passed to 
#'   \code{\link[maxLik]{maxLik}}  for fitting the model.
#' @details The model corresponds to an extensive-form, 
#'   discrete-crisis-bargaining game from Lewis and Schultz (2003): 
#'   \preformatted{ 
#'   .       A 
#'   .      / \ 
#'   .     /   \ 
#'   .    /     \
#'   .   S_A     B 
#'   .    0     / \
#'   .         /   \
#'   .        /     \
#'   .      V_A      A 
#'   .      C_B     / \ 
#'   .             /   \
#'   .            /     \
#'   .     W_A + e_A    a + e_a 
#'   .     W_B + e_B    V_B} 
#'   If \eqn{A} chooses not to challenge \eqn{B},
#'   then the game ends at the leftmost node (\eqn{SQ}) and payoffs are
#'   \eqn{S_A} and 0 to players \eqn{A} and \eqn{B}, respectively. If \eqn{A}
#'   challenges \eqn{B}, \eqn{B} can concede or resist.  If \eqn{B} concedes,
#'   the game ends at \eqn{CD} with payoffs \eqn{V_A} and \eqn{C_B}.  However,
#'   if \eqn{B} resists, \eqn{A} decides to stand firm, which ends the game at
#'   \eqn{SF} with payoffs \eqn{W_A + \epsilon_A} and \eqn{W_B + \epsilon_B}.
#'   Finally, if \eqn{A} decides to back down in the face of \eqn{B}'s
#'   resistance, then the game ends at the rightmost node \eqn{BD}, with payoffs
#'   \eqn{a + \epsilon_a} and \eqn{V_B}.
#'   
#'   The seven right-hand formulas that are specified in the formula argument 
#'   correspond to the regressors to be placed in \eqn{S_A, V_A, C_B, W_A, W_B, 
#'   a}, and \eqn{V_B}, respectively. The model is unidentified if any regressor
#'   (including a constant term) is included in all the formulas for each player
#'   (Lewis and Schultz 2003). Often the easiest way to meet this requirement is
#'   set one formula per player  to 0. When an identification problem is
#'   detected, an error is issued. For example, the syntax for the formula
#'   argument could be:
#'   
#'   \code{formulas = sq + cd + sf + bd ~ x1 + 0 | x2 | x2 | x1 + x2 | x1 | 1 | 0)}
#' 
#'   Where: \itemize{
#'       \item \code{sq + cd + sf + bd} are the tallies of how many
#'   times each outcome is observed for each observation.  When the game is only
#'   observed once, that observation will be a 1 and three 0s.  When the game is
#'   observed multiple times, these variables should count the number of times 
#'   each outcome is observed.  They need to be in the order of \eqn{SQ}, 
#'   \eqn{CD}, \eqn{SF}, \eqn{BD}.
#'    \item \eqn{S_A} is a function of the variable \code{x1} and no constant term.
#'    \item \eqn{V_A} is a function of the variable \code{x2} and a constant term.
#'    \item \eqn{C_B} is a function of the variable \code{x2} and a constant term.
#'    \item \eqn{W_A} is a function  of the variables \code{x1}, \code{x2} and a constant term.
#'    \item \eqn{W_B} is a function of the variable  \code{x1} and a constant term.
#'   \item \eqn{a} is a constant term.
#'   \item \eqn{V_B} is fixed to 0 (or a non-zero value set by \code{fixed.par}. }
#
#'
#'
#'   Each row of the data frame should be a summary of the covariates and outcomes associated with that particular game.
#'   When each game is observed only once, then this will resemble an ordinary dyad-time data frame.
#'   However, if there are multiple observations per game, then each row should be a summary of all the data associated
#'   with that game.
#'   For example, if there are \eqn{D} games in the data, where each is observed \eqn{T_d} times, then the data frame
#'   should have \eqn{D} rows.
#'   The four columns making up the dependent variable will denote the frequencies of each outcome for game \eqn{d},
#'   such that \code{sq}\eqn{_d} + \code{cd}\eqn{_d} + \code{sf}\eqn{_d} + \code{bd}\eqn{_d = T_d}.
#'   The covariates in row \eqn{d} should be summary statistics for the exogenous variables (e.g., mean, median, mode, first observation).
#'   
#'    
#'   The model is first fit using a pseudo-likelihood estimator.  This approach 
#'   requires first stage estimation of the probability that \eqn{B} resists and
#'   the probability that \eqn{A} fights conditional on \eqn{B} choosing to 
#'   resist. These first stage estimates should be flexible and we recommend
#'   that users fit a flexible semi-parametric or non-parametric model to
#'   produce them. If these estimates are produced by the analyst prior to using
#'   this function, then they can be provided by providing a list to the
#'   \code{phat} argument. This list should contain two named elements \itemize{
#'   \item \code{PRhat} is the probability that \eqn{B} resists.  This should be
#'   a vector of probabilities with one estimated probability for each
#'   observation. \item \code{PFhat} is the probability that \eqn{A} stands firm
#'   conditional on \eqn{B} resisting.  This should be a vector of probabilities
#'   with one estimated probability for each observation. }
#'   
#'   
#'   If the user leaves the \code{phat} argument empty, then these first-stage 
#'   estimates are produced internally using the 
#'   \code{\link[randomForest]{randomForest}} function.
#'   Users wanting to use the
#'   random forest, can supply a formula for it using the argument 
#'   \code{phat.formulas}.
#'   This argument can take a formula with nothing on the 
#'   left-hand side and 1-2 right-hand sides.
#'   If two right-hand sides are 
#'   provided then the first is used to generate \code{PRhat}, and the second is
#'   used for \code{PFhat}.
#'   If only one right-hand side is provided, it is used
#'   for both. Some examples: \itemize{
#'        \item \code{phat.formulas = ~ x1 + x2} 
#'   predict \code{PRhat} and \code{PFhat} using \code{x1} and \code{x2}.
#'        \item \code{phat.formulas = ~ x1 + x2 | x1 + x2} predict \code{PRhat} and 
#'   \code{PFhat} using \code{x1} and \code{x2}
#'        \item \code{phat.formulas = ~ x1 + x2 | x1 } predict \code{PRhat}  using \code{x1} and \code{x2}, but
#'   predict \code{PFhat} using only \code{x1}. }
#'   If both \code{phat} and \code{phat.formula} are missing, then a random forest is fit using all the 
#'   exogenous variables listed in the formulas argument
#'   
#'   If \code{method = "npl"}, then estimation continues.  
#'   For each iteration of the NPL, the estimates of \code{PRhat} and \code{PFhat} are updated
#'   by one best-response iteration using the current parameter estimates.
#'   The model is then refit using these updated choice probabilities.
#'   This process continues until the maximum absolute change in
#'   parameters and choice probabilities is less than \code{npl.tol} (default, \code{1e-7}), or
#'   the number of outer iterations exceeds \code{npl.maxit} (default, \code{25}).
#'   In the latter case, a warning is produced.
#'   
#'   
#'   If pseudo-likelihood (\code{method="pl"}) is used, then 
#'   \code{pl.vcov} is checked.
#'   There are four possibilities here: 
#'   \itemize{
#'     \item \code{pl.vcov = FALSE} (default), then no covariance matrix or
#'   standard errors are returned, only the point estimates. 
#'     \item \code{pl.vcov > 0} and \code{phat.vcov} is supplied, 
#'   then \code{phat.vcov} is used to estimate the PL's covariance matrix. 
#'     \item  \code{pl.vcov > 0}, \code{phat.vcov} is missing, and \code{phat} 
#'   is missing, then the random forest used to estimate \code{PRhat} and
#'   \code{PFhat} is bootstrapped (simple, nonparametric bootstrap)  \code{pl.vcov} times.
#'     \item \code{pl.vcov > 0}, \code{phat.vcov} is missing, and \code{phat} is not
#'   missing, then an error is returned. 
#'   }
#'   
#' @return An object of class \code{sigfit}, containing:\describe{ 
#'   \item{\code{coefficients}}{A vector of estimated model parameters.}  
#'   \item{\code{vcov}}{Estimated variance-covariance matrix. When \code{pl.vcov = FALSE}, this slot is omitted.} 
#'   \item{\code{utilities}}{Each actor's utilities at the estimated values.}  
#'   \item{\code{fixed.par}}{The fixed utilities if specified in the call.}  
#'   \item{\code{logLik}}{Final log-likelihood value of the model.}
#'   \item{\code{gradient}}{First derivative values at the estimated parameters.}
#'   \item{\code{Phat}}{List of two elements 
#'      \itemize{ 
#'        \item \code{PRhat} The first stage estimates of the probability that 
#'      \eqn{B} resists (\code{method = "pl"}) or the final estimates that
#'      \eqn{B} resists (if \code{method = "npl"}) 
#'         \item \code{PFhat} The first stage estimates of the probability that 
#'      \eqn{A} stands firms given that \eqn{A} challenged (\code{method = "pl"}) or 
#'      the final estimates that \eqn{A} stands firms given that \eqn{A} challenged
#'      (if \code{method = "npl"})
#'      } 
#'      Note that \code{PRhat} will only be an equilibrium if \code{method = "npl"} and the NPL convergences
#'   }
#'   \item{\code{user.phat}}{Logical. Did the user provide phat?} 
#'   \item{\code{start.beta}}{The vector of starting values used in the PL optimization.}
#'   \item{\code{call}}{The call used to produce the object.}
#'   \item{\code{model}}{The data frame used to fit the model.}
#'   \item{\code{method}}{The method (\code{"pl"} or \code{"npl"}) used to fit the model.}
#'   \item{\code{maxlik.method}}{The optimization used by \code{maxLik} to fit the model.}
#'   \item{\code{maxlik.code}}{The convergence code returned by \code{maxLik}.}
#'   \item{\code{maxlik.message}}{The convergence message returned by \code{maxLik}.}
#'   }
#' Additionally, when \code{method = "npl"}, the following are also included in the \code{sigfit} object.\describe{
#'   \item{\code{npl.iter}}{Number of best response iterations used in fitting the NPL.}
#'   \item{\code{npl.eval}}{Maximum difference between the parameters at the last two NPL iterations. If the NPL method converged, this should be less than \code{npl.tol} specified in the function call.}
#'   \item{\code{eq.constraint}}{Maximum equilibrium constraint violation.}
#'   }
#' 
#' 
#' 
#' @references 
#' Casey Crisman-Cox and Michael Gibilisco. 2019. "Estimating 
#' Signaling Games in International Relations: Problems and Solutions."
#' \emph{Political Science Research and Methods}. Online First.
#' 
#' Jeffrey B. Lewis and Kenneth A. Schultz.  2003.  "Revealing
#' Preferences: Empirical Estimation of a Crisis Bargaining Game with
#' Incomplete Information."  \emph{Political Analysis} 11:345--367.
#'   
#' @examples 
#' data("sanctionsData")
#' f1 <- sq+cd+sf+bd ~ sqrt(senderecondep) + senderdemocracy + contig + ally -1|#SA
#'                     anticipatedsendercosts|#VA
#'                     sqrt(targetecondep) + anticipatedtargetcosts + contig + ally|#CB
#'                     sqrt(senderecondep) + senderdemocracy + lncaprat | #barWA
#'                     targetdemocracy + lncaprat| #barWB
#'                     senderdemocracy| #bara
#'                     -1#VB
#'
#' ## Using Nested-Pseudo Likelihood  with default first stage     
#' \dontrun{            
#' fit1 <- sigint(f1, data=sanctionsData, npl.trace=TRUE)
#' summary(fit1)
#' }
#' 
#' 
#' ## Using Pseudo Likelihood with user supplied first stage
#' Phat <- list(PRhat=sanctionsData$PRnpl, PFhat=sanctionsData$PFnpl)
#' fit2 <- sigint(f1, data=sanctionsData, method="pl", phat=Phat)
#' summary(fit2)
#' 
#' ## Using Pseudo Likelihood with user made first stage and user covariance
#' ## SIGMA is a bootstrapped first-stage covariance matrix (not provided)
#' \dontrun{
#' fit3 <- sigint(f1, data=sanctionsData, method="pl", phat=Phat, phat.vcov=SIGMA, pl.vcov=TRUE)
#' summary(fit3)
#' }
#' 
#' ## Using Pseudo Likelihood with default first stage and 
#' ## bootstrapped standard errors for the first stage covariance
#' \dontrun{
#' fit4 <- sigint(f1, data=sanctionsData, method="pl", pl.vcov=25) 
#' summary(fit4)
#' }
#' 
#' @import pbivnorm
#' @import randomForest
#' @import Formula
#' @import stats
#' @import utils
#' @import MASS
#' @import maxLik
#' @export
sigint <- function(formulas, data, subset, na.action,
                   fixed.par=list(),
                   method=c("npl", "pl"),
                   npl.maxit=25,
                   npl.tol=1e-7,
                   npl.trace=FALSE,
                   start.beta,
                   maxlik.method="NR",
                   phat,
                   phat.formulas,
                   pl.vcov=FALSE,
                   phat.vcov,
                   seed=12345,
                   maxlik.options=list()){
  
  
  
  
  #Basic prelminaries.
  maxlik.options.default <- list(tol=1e-10,
                                 reltol=1e-10,
                                 gradtol=1e-10,
                                 iterlim=100)
  maxlik.options <- modifyList(maxlik.options.default, maxlik.options)
  cl <- match.call()
  method <- match.arg(method)
  
  ####Formulas####
  ##process main formula
  formulas <- as.Formula(formulas)
  if (length(formulas)[2] != 7){
    stop("'formulas' should have seven components on the right-hand side")
  }
  
  
  
  ####Data####
  ## make the model frame
  mf <- match(c("data", "subset", "na.action"), names(cl), 0L)
  mf <- cl[c(1L, mf)]
  mf$formula <- formulas
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  #make the dependent variable
  Y <- model.part(formulas, mf, lhs = 1, drop = TRUE)
  Y <- t(as.matrix(Y)) ##Because that's how Mike set it up
  
  
  #make the independent variables
  regr <- list()
  for(i in 1:7){
    regr[[i]] <- model.matrix(formulas, data=mf, rhs=i)
  }
  names(regr) <- c("SA", "VA", "CB", "barWA", "barWB", "bara", "VB")
  ncols <- sum(sapply(regr, ncol))
  ngames <- nrow(regr[[1]])
  u.names <- rep(names(regr), sapply(regr, ncol))
  
  
  #identification check
  v.names <- lapply(regr, colnames)
  A.names <- unique(unlist(v.names[c(1,2,4,6)]))
  check.A <- colSums(do.call(rbind,lapply(v.names[c(1,2,4,6)], function(x){A.names %in% x})))==4
  B.names <- unique(unlist(v.names[c(3,5,7)]))
  check.B <- colSums(do.call(rbind,lapply(v.names[c(3,5,7)], function(x){B.names %in% x})))==3
  if(any(check.B) | any(check.A)){
    actors <- c()
    error.list <- list()
    error.A <- A.names[which(check.A)]; if(length(error.A)>0){actors <- c(actors, "A's"); error.list$A <- error.A}
    error.B <- B.names[which(check.B)]; if(length(error.B)>0){actors <- c(actors, "B's"); error.list$B <- error.B}
    error.A <- stringr::str_c(error.A, collapse = ", ")
    error.B <- stringr::str_c(error.B, collapse = ", ")
    
    error.message <- stringr::str_c(c("This model is unidentified.",
                                      paste("The following variables appear in all of", actors, "utilities:", error.list)),
                                    collapse = "\n")
    stop(error.message)                           
    
  }
  
  
  if(length(fixed.par)>0){
    if(any(sapply(regr[names(fixed.par)], ncol) != 0)){
      stop(paste("The following payoff is fixed but their formula component is not 0",
        names(which(sapply(regr[names(list(VA=1, VB=1))], ncol) != 0)), "\n"))
    }
  }
  
  
  
  #######First stage######
  init.seed <- runif(1) #just to make sure .Random exists
  old.seed <- .Random.seed
  on.exit(assign(".Random.seed", old.seed, envir=globalenv()))
  set.seed(seed)  
  ##Formula for first stage estimates.  Needed for ls and twostep
  if(missing(phat)){ #CHECK FOR USER VALUES
    USER.PHAT <- FALSE
    index1 <- colSums(Y[2:4,]) >= 1
    index2 <- colSums(Y[3:4,]) >= 1
    
    if(missing(phat.formulas)){ #NO USER VALUES, HOW ABOUT A FORMULA?
      #If neither, use all the variables in regr
      warning("No formula or first stage estimates found, using all observables in first stage random forest")
      
      X <- do.call(cbind, regr)
      X1 <- X2 <-  unique(X, MARGIN=2)
    }else{ #FORMULA PROVIDED, USE IT
      phat.formulas <- as.Formula(phat.formulas)
      if (length(phat.formulas)[2] ==1){
        phat.formulas <- as.Formula(formula(phat.formulas, rhs=c(1,1)))
      }
      ## make the model frame
      mf.phat <- match(c("data", "subset", "na.action"), names(cl), 0L)
      mf.phat <- cl[c(1L, mf.phat)]
      mf.phat$formula <- phat.formulas
      mf.phat$drop.unused.levels <- TRUE
      mf.phat[[1]] <- as.name("model.frame")
      mf.phat <- eval(mf.phat, parent.frame())
      
      X1 <- model.matrix(phat.formulas, data=mf.phat, rhs=1)
      X2 <- model.matrix(phat.formulas, data=mf.phat, rhs=2)
    }
    X1 <- X1[,which(colnames(X1) != "(Intercept)"), drop=FALSE]
    X2 <- X2[,which(colnames(X2) != "(Intercept)"), drop=FALSE]
    
    phat.X <- list(X1 = cbind((colSums(Y[3:4,])/colSums(Y[2:4,]))[index1], X1[index1,]),
                   X2 = cbind(((Y[3,])/colSums(Y[3:4,]))[index2], X2[index2,]))
    phat.X <- lapply(phat.X, as.data.frame)
    names(phat.X$X1) <- c("V1", colnames(X1))
    names(phat.X$X2) <- c("V1", colnames(X2))
    if(length(unique(phat.X$X1[,1])) == 2){
      phat.X$X1[,1] <- as.factor(phat.X$X1[,1])
      type1 <- "prob"
    }else{
      type1 <- "response"
    }
    if(length(unique(phat.X$X1[,1])) == 2){
      phat.X$X2[,1] <- factor(phat.X$X2[,1]); type2 <- "prob"
    }else{
      type2 <- "response"
    }
    

    m1 <- randomForest(y=phat.X$X1[,1], x=phat.X$X1[,-1, FALSE])
    m2 <- randomForest(y=phat.X$X2[,1], x=phat.X$X2[,-1, FALSE])
    
    phat <- list(PRhat = predict(m1, newdata=X1, type=type1),
                 PFhat = predict(m2, newdata=X2, type=type2))
    if(type1 == "prob"){phat$PRhat <- phat$PRhat[,2]}
    if(type2 == "prob"){phat$PFhat <- phat$PFhat[,2]}
    
  }else{ #VALUES PROVIDED USE THEM
    USER.PHAT <- TRUE
    if(! "list" %in%  class(phat)){stop("phat must be a list")}
    if(length(phat) !=2){stop("phat must be a list of length 2")}
    if(!all(c("PRhat", "PFhat") %in% names(phat))){stop("Elements in phat must be names 'PRhat' and 'PFhat'")}
    if(method=="pl"  & (pl.vcov>0) & missing(phat.vcov)){stop("PL covariance matrix can only be computed if phat covariance matrix is supplied")}
  }
  
  ######## PL Estimates #########
  
  # start values
  set.seed(seed)
  if(missing(start.beta)){
    start.beta <- rnorm(ncols)*0.05
  }
  names(start.beta) <- unlist(lapply(regr,colnames))
  names(start.beta) <- paste(u.names, ":", names(start.beta), sep="")
  
  
  fqll <- function(x){-QLL.jo(x,phat$PRhat,phat$PFhat,Y,regr, fixed.par)}
  grqll <- function(x){-eval_gr_qll(x, phat$PRhat, phat$PFhat,Y,regr, fixed.par)}
  
  #testing
  # cat("LL:", fqll(start.beta), "\n")
  # cat("Grad:", grqll(start.beta), "\n")
  # cat("Par:", start.beta, "\n")
  
  out <- maxLik::maxLik(start=start.beta, fqll,
                        gr=grqll,
                        method=maxlik.method,
                        control=maxlik.options)
  
  Phat0 = phat
  eq.const.pl <- const.jo(Phat0$PRhat, vec2U.regr(x=out$estimate, regr=regr, fixed.par=fixed.par))  
  ####### NPL Procedure #####
  if(method=="npl"){
    tol = npl.tol
    maxit = npl.maxit+1
    eval = 1000
    Phat = phat
    out.NPL = out
    iter = 0
    
    while(eval > tol & iter < maxit){
      Uk <- vec2U.regr(out.NPL$estimate, regr, fixed.par)
      Pk.F <- eqProbs(Phat$PRhat, Uk, RemoveZeros = T)[,3]
      Pk.R <- pnorm((Phat$PFhat*Uk$barWB + (1-Phat$PFhat)*Uk$VB - Uk$CB)/Phat$PFhat)
      Phat.k_1 <- Phat
      Phat <- list(PRhat = Pk.R, PFhat = Pk.F)
      Phat$PRhat <-  pmin(pmax(Phat$PRhat,
                               0.0001), .9999)
      Phat$PFhat <-  pmin(pmax(Phat$PFhat,
                               0.0001), .9999)
      fqll <- function(x){- QLL.jo(x,Phat$PRhat,Phat$PFhat,Y,regr, fixed.par)}
      grqll <- function(x){- eval_gr_qll(x, Phat$PRhat, Phat$PFhat,Y,regr, fixed.par)}
      out.NPL.k <- try(maxLik::maxLik(start=out.NPL$est, logLik=fqll, grad=grqll,
                                      method=maxlik.method, control=maxlik.options))
      if(class(out.NPL.k[[1]])=="character" || out.NPL.k$code==100){
        out.NPL <- out.NPL.k
        break
      }
      out.NPL.k$convergence <- out.NPL.k$code
      eval <- max( abs(c(out.NPL.k$est, unlist(Phat.k_1)) -c(out.NPL$est, unlist(Phat))))
      if(npl.trace){
        cat("Iteration ", iter, ": Convergence Criteria: ", eval, "\n")
      }
      out.NPL <- out.NPL.k
      iter <- iter + 1
    }
    out.NPL$iter <- iter
    if(class(out.NPL[[1]])=="character"|| out.NPL.k$code==100){
      stop("NPL failed. Consider using NR method, different starting values, or using only pl method")
    }
    if(iter == maxit){
      warning("NPL iterations exceeded")
    }
    Uk <- vec2U.regr(out.NPL$estimate, regr, fixed.par)
    Pk.F <- eqProbs(Phat$PRhat, Uk, RemoveZeros = T)[,3]
    Pk.R <- pnorm((Phat$PFhat*Uk$barWB + (1-Phat$PFhat)*Uk$VB - Uk$CB)/Phat$PFhat)
    Phat.k_1 <- Phat
    Phat <- list(PRhat = Pk.R, PFhat = Pk.F)
    eq.const <- const.jo(Phat$PRhat, Uk)
    
    ###NPL VCOV###
    Dtheta <- eval_gr_qll.i(out.NPL$estimate,Phat$PRhat, Phat$PFhat, Y, regr, fixed.par )
    Dtheta.theta <- crossprod(Dtheta)
    
    Dp <- eval_gr_qll.ip(Phat$PRhat, Phat$PFhat, out.NPL$est, Y, regr, fixed.par)
    Dtheta.p <- crossprod(Dtheta,Dp)
    
    JpPsi <- dPsiDp(Phat$PRhat, Phat$PFhat, out.NPL$est, Y, regr, fixed.par)
    JtPsi <- dPsi.dTheta(out.NPL$estimate,Phat$PRhat, Phat$PFhat, Y, regr, fixed.par)
    
    tJpPsi.inv <-  tryCatch(solve(diag(nrow(JpPsi)) - t(JpPsi)),
                           error=function(x){
                             warning("Possible singular matrix detected, using Pseudoinverse");
                             return(MASS::ginv(diag(nrow(JpPsi)) - t(JpPsi)))
                           })
    JpPsi.inv <-  tryCatch(solve(diag(nrow(JpPsi)) - (JpPsi)),
                           error=function(x){
                             warning("Possible singular matrix detected, using Pseudoinverse");
                             return(MASS::ginv(diag(nrow(JpPsi)) - (JpPsi)))
                           })
    
    topSlice <- tryCatch(solve(Dtheta.theta  +
                                 Dtheta.p %*% tJpPsi.inv %*%   JtPsi),
                         error=function(x){
                           warning("Possible singular matrix detected, using Pseudoinverse");
                           return(MASS::ginv(Dtheta.theta  +
                                               Dtheta.p %*% tJpPsi.inv %*%   JtPsi))
                         })
    bottomSlice <- tryCatch(solve(Dtheta.theta  +
                                     t(JtPsi) %*% JpPsi.inv %*% t(Dtheta.p)),
                         error=function(x){
                           warning("Possible singular matrix detected, using Pseudoinverse");
                           return(MASS::ginv(Dtheta.theta  +
                                               t(JtPsi) %*% JpPsi.inv %*% t(Dtheta.p)))
                         })
    
    
    # topSlice <- solve(Dtheta.theta  +
    #                     Dtheta.p %*% solve(diag(nrow(JpPsi)) - t(JpPsi)) %*%   JtPsi)
    # bottomSlice <- solve(Dtheta.theta  +
    #                        t(JtPsi) %*% solve(diag(nrow(JpPsi)) - t(JpPsi)) %*% t(Dtheta.p))
    vcov.NPL <- topSlice %*% Dtheta.theta %*% bottomSlice
    rownames(vcov.NPL) <- colnames(vcov.NPL) <- names(start.beta)
    # cat("gradient", out.NPL$gradient, "\n")
    # cat("gradient2",  grqll(out.NPL$estimate), "\n")
    # cat("gradient", out.NPL.k$gradient, "\n")
    #### Build npl.out ####
    npl.out <- list(coefficients = out.NPL$estimate,
                    fixed.par=fixed.par,
                    utilities = vec2U.regr(out.NPL$estimate,regr,fixed.par),
                    vcov=vcov.NPL,
                    logLik = out.NPL$maximum,
                    gradient = out.NPL$gradient,
                    Phat = Phat,
                    eq.constraint = eq.const,
                    user.phat=USER.PHAT,
                    start.values=start.beta,                    
                    call=cl,
                    formulas = formulas,
                    model=mf,
                    method=method,
                    maxlik.method=maxlik.method,
                    maxlik.code = out.NPL$code,
                    maxlik.message = out.NPL$message,
                    npl.iter=iter,      
                    npl.eval = eval)
    class(npl.out) <- "sigfit"
    return(npl.out)
  }else{
    #### PL output ###
    pl.out <- list(coefficients = out$estimate,
                   fixed.par=fixed.par,
                   utilities = vec2U.regr(out$estimate,regr,fixed.par),
                   logLik = out$maximum,
                   gradient = out$gradient,
                   Phat = Phat0,
                   eq.constraint = eq.const.pl,
                   user.phat=USER.PHAT,
                   start.values=start.beta,
                   call=cl,
                   formulas = formulas,
                   model=mf,
                   method=method,
                   maxlik.method=maxlik.method,
                   maxlik.code=out$code,
                   maxlik.message=out$message)
    
    #### PL VCOV ####
    if(!missing(phat.vcov)){pl.vcov <- TRUE}
    if(pl.vcov>0){
      if(USER.PHAT | !missing(phat.vcov)){
        SIGMA <- phat.vcov
      }else{
        if(USER.PHAT & missing(phat.vcov)){
          warning("User supplied phat, but phat.vcov is missing. Standard errors not returned")
        }else{
          Phat.boot <- matrix(0, ncol=ncol(Y)*2, nrow=pl.vcov)
          for(i in 1:pl.vcov){
            Y.boot <- 0
            while(mean(Y.boot) == 0 | mean(Y.boot) == 1){
              use <- sample(1:ncol(Y), replace=T)
              Y.boot <- Y[,use]
              X1.boot <- X1[use,,drop=FALSE]
              X2.boot <- X2[use,,drop=FALSE]
            }
            index1.boot <- colSums(Y.boot[2:4,]) >= 1
            index2.boot <- colSums(Y.boot[3:4,]) >= 1
            phat.X.boot <- list(X1 = cbind((colSums(Y.boot[3:4,])/colSums(Y.boot[2:4,]))[index1.boot],
                                           X1.boot[index1.boot,]),
                                X2 = cbind(((Y.boot[3,])/colSums(Y.boot[3:4,]))[index2.boot],
                                           X2.boot[index2.boot,]))
            
            
            
            phat.X.boot <- lapply(phat.X.boot, as.data.frame)
            phat.X.boot <- lapply(phat.X.boot, function(x){colnames(x) <- c("y", colnames(X1)); return(x)})
            if(length(unique(phat.X.boot$X1[,1])) == 2){
              phat.X.boot$X1[,1] <- as.factor(phat.X.boot$X1[,1])
              type1 <- "prob"
            }else{
              type1 <- "response"
            }
            if(length(unique(phat.X.boot$X1[,1])) == 2){
              phat.X.boot$X2[,1] <- factor(phat.X.boot$X2[,1]); type2 <- "prob"
            }else{
              type2 <- "response"
            }
            
            m1 <- randomForest(y=phat.X.boot$X1[,1], x=phat.X.boot$X1[,-1,drop=FALSE])
            m2 <- randomForest(y=phat.X.boot$X2[,1], x=phat.X.boot$X2[,-1,drop=FALSE])
            
            phat.out <- list(PRhat = predict(m1, newdata=X1, type=type1),
                             PFhat = predict(m2, newdata=X2, type=type2))
            
            if(type1 == "prob"){phat.out$PRhat <- phat.out$PRhat[,2]}
            if(type2 == "prob"){phat.out$PFhat <- phat.out$PFhat[,2]}
            
            
            
            Phat.boot[i,] <- unlist(phat.out)
          }
          SIGMA <- var(Phat.boot, na.rm=TRUE)
        }
      }
      Dtheta <- eval_gr_qll.i(out$estimate, Phat0$PRhat, Phat0$PFhat, Y, regr, fixed.par)
      # Dtheta <- numDeriv::hessian(QLL.jo, out$estimate, PRhat=Phat0$PRhat, PFhat=Phat0$PFhat, Y=Y, regr=regr)
      Dtheta.theta <- crossprod(Dtheta)
      Dtheta.theta.inv <- tryCatch(solve(Dtheta.theta),
                                   error=function(x){
                                     warning("OPG is possibly singular, using Pseudoinverse");
                                     return(MASS::ginv(Dtheta.theta))
                                   })
      
      Dp <- eval_gr_qll.ip(Phat0$PRhat, Phat0$PFhat, out$est, Y, regr, fixed.par) #NAs are appearing here
      Dtheta.p <- crossprod(Dtheta,Dp) 
      
      
      vcov.PL <- Dtheta.theta.inv + Dtheta.theta.inv %*% Dtheta.p %*% SIGMA %*% t(Dtheta.p) %*% Dtheta.theta.inv
      rownames(vcov.PL) <- colnames(vcov.PL) <- names(start.beta)
      
      pl.out$vcov <- vcov.PL
    }
    
    class(pl.out) <- "sigfit"
    return(pl.out)
  }
}
