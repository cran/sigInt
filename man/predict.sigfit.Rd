% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.sigfit.r
\name{predict.sigfit}
\alias{predict.sigfit}
\title{Predicted probabilities and comparative statics for signaling games}
\usage{
\method{predict}{sigfit}(
  object,
  newdata,
  new.theta,
  type = c("actions", "outcomes"),
  na.action = na.pass,
  control = list(),
  parallel = FALSE,
  ...
)
}
\arguments{
\item{object}{a fitted model of class \code{sigfit}.}

\item{newdata}{data frame of covariates used to produce the predicted probabilities.
If this is left empty, the entire original data set is used.
When \code{newdata} is specified, \code{new.theta} should be either missing (use the coefficients
from the \code{sigfit} object) or a one row data frame. See "Details" for more information.
As with other \code{\link[stats]{predict}} methods, variable names must match those used to fit the model.}

\item{new.theta}{a data frame of alternative parameters for comparative statics.
When missing, the coefficients from the \code{object} are used.
When specified, each row should be a complete parameter vector. 
If \code{new.theta} is specified, then \code{newdata} must be a data frame with only one row.
Unlike other \code{\link[stats]{predict}} methods, column names do not matter here.
Instead, the columns must be the same order as the coefficients in \code{object}.
See "Details" and "Examples" for more information.}

\item{type}{whether to provide probabilities over actions 
(default, returns \eqn{p_C}, \eqn{p_R}, and \eqn{p_F}) 
or outcomes (returns \eqn{SQ}, \eqn{CD}, \eqn{SF}, and  \eqn{BD}).}

\item{na.action}{how to deal with \code{NA}s in \code{newdata}.}

\item{control}{list of options describing the grid search method. See "Details" for more information}

\item{parallel}{logical. Should the comparative statics be computed in parallel, requires the
\code{\link[parallel]{parallel}}
package be installed. Parallelization is done using \code{\link[parallel]{parSapply}}.}

\item{...}{Additional arguments (not currently used)}
}
\value{
An object of class \code{sigProb} containing three elements:
\describe{
 \item{\code{predicted}}{data frame of predicted probabilities. The first column of this data frame is 
 called \code{Row}, which corresponds to the rows in either \code{model} or \code{par}.
 In the event of multiple equilibria, this column allows for mapping data and parameters to 
 all computed equilibria.
 }
 \item{\code{model}}{data frame of covariates used to produce the predicted probabilities.}
 \item{\code{par}}{data frame of parameters used to produce the predicted probabilities.}
 }
}
\description{
This method uses a fitted model of class \code{sigfit} to compute
predicted probabilities or comparative statics. 
Users can provide either new data or a new parameters to generate 
counterfactuals of interest.
}
\details{
This function is used to consider comparative statics in the crisis signaling game.
The model of interest is fit using \code{\link{sigint}}.
How this function behaves largely depends on  how \code{newdata} and \code{new.theta} are
specified.

When both \code{newdata} and \code{new.theta} are missing, all equilibria for every
observation used to fit the model are computed.
These equilibria are then used to calculated either choice probabilities (default,
\code{type = "action"}) the distribution over outcomes (\code{type = "outcomes"}).

When only \code{newdata} is specified, then all equilibria are computed using the
data frame in \code{newdata} and the coefficients from \code{object}. This produces
standard comparative statics with respect to observed covariates.

When \code{newdata} is specified and \code{new.theta} is a one row data frame,
then all equilibria are computed using the data frame in \code{newdata} and the
coefficients from \code{new.theta}.

When \code{newdata} is a one row data frame and \code{new.theta} is specified,
then all equilibria are computed using the data frame in \code{newdata} and the
coefficients from \code{new.theta}. This is a comparative static on changing a
structural parameter in the model.

If \code{new.theta} has more than one row, then \code{newdata} must be specified as a
data frame with only one row.  Anything else returns an error.

Equilibria are computed using a line search method.
The \code{control} argument allows for user control over this process.
Users can specify a list with the following named elements
\describe{
 \item{gridsize}{Integer. The number of points considered in the line search (default, 1e4).
 More points makes it more likely that all equilibria are discovered, but can slow down the
 search. }
 \item{comp}{Logical. Should an equilibrium be computed when discovered?
            When \code{comp = FALSE} (default), the mean of the grid points surrounding the
            equilibrium is used as an approximate solution.
            When \code{comp = TRUE} the function \code{\link[stats]{uniroot}}
            is called to find a more precise solution to the equilibrium constraint problem.}
 \item{tol}{Numeric. When \code{comp = TRUE}, this is the tolerance used by
           \code{\link[stats]{uniroot}}}.
}

When dealing with a larger problem, such as computing all equilibria for every observation,
it can be helpful to parallelize process.  If the user has the
(suggested) \code{\link[parallel]{parallel}} package, then the option \code{parallel = TRUE}
will use the function \code{\link[parallel]{parSapply}} is used.
}
\examples{
data(sanctionsData)

f1 <- sq+cd+sf+bd ~ sqrt(senderecondep) + senderdemocracy + contig + ally -1|#SA
                    anticipatedsendercosts|#VA
                    sqrt(targetecondep) + anticipatedtargetcosts + contig + ally|#CB
                    sqrt(senderecondep) + senderdemocracy + lncaprat | #barWA
                    targetdemocracy + lncaprat| #barWB
                    senderdemocracy| #bara
                    -1#VB
## Using Nested-Pseudo Likelihood  with default first stage    
\dontrun{             
fit1 <- sigint(f1, data=sanctionsData, npl.trace=TRUE)
p.out <- predict(fit2, parallel=TRUE) #fitted choice probabilites for all observations
}

## Outcome probabilities for first five using PL method
Phat <- list(PRhat=sanctionsData$PRhat, PFhat=sanctionsData$PFhat)
fit2 <- sigint(f1, data=sanctionsData, method="pl", phat=Phat)
p1 <- predict(fit2, newdata=sanctionsData[1:5,], type="outcome")

## comparative static on \bar{a}, compute more precise equilibria with uniroot
new.theta <- data.frame(t(replicate(25, coef(fit2))))
new.theta[,19] <- seq(-6, 0, length=25)
p2 <- predict(fit2, newdata=sanctionsData[1,], new.theta=new.theta, control=list(comp=TRUE))


}
\seealso{
\code{\link{plot.sigProb}}, \code{\link{generate.eq}}
}
