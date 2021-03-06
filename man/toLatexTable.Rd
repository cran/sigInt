% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toLatex.sigfit.r
\name{toLatexTable}
\alias{toLatexTable}
\title{Export a \code{sigfit} object into paper-ready LaTeX table}
\usage{
toLatexTable(
  ...,
  se.list,
  stars = c("default", "all", "none"),
  caption = "",
  label,
  align,
  digits = 2,
  se.note = "Standard errors in parenthesis",
  order,
  covariate.labels,
  model.names,
  dep.varnames,
  k = 1,
  print.xtable.options = list()
)
}
\arguments{
\item{...}{one or more models fit using \code{\link{sigint}}.}

\item{se.list}{an optional list where each element contains the standard errors of the models.
If included this list must include one element per model, even if that model's standard errors are unchanged.
This argument should be used when standard errors have been adjusted outside the model (e.g., a bootstrap).
If left empty the standard errors from the model's covariance matrix are used if available.
Within each element of this list, the order of standard errors must be the same order as the coefficients in the model.
Standard errors are not matched on name.}

\item{stars}{how should significance stars be used? \code{stars = "default"} returns a single star when \eqn{p < 0.05},
\code{stars = "all"} returns flags for \eqn{p < 0.1}, \eqn{p < 0.05}, and \eqn{p < 0.01}, and \code{stars = "none"} 
returns no stars at all.}

\item{caption}{a string to be used as the table's caption.}

\item{label}{a string to be used as the table's LaTeX \code{label} argument.}

\item{align}{a string to indicate the alignment of each column in the table.
Passed directly to the LaTeX \code{tabular} options.}

\item{digits}{how many digits after the decimal point should be displayed? Default 2.}

\item{se.note}{a string containing a note for the bottom of the table. 
Default is the common "Standard errors in parenthesis.}

\item{order}{a string vector describing the order that  the covariates should be in the table.
The default is to use the order they're listed in the \code{sigfit} object. 
When this vector is shorter than the coefficient vector, variables are first included by \code{order},
remaining variables are included based on their order in the model}

\item{covariate.labels}{a string vector of "nice" names for the variables appropriate for a published work.
If empty, the "ugly" names from the fitted model are used. Note that if \code{order} is specified then the
covariate labels must match the order in \code{order}.}

\item{model.names}{an optional vector of model names to include as column titles in the table. 
Should be either a single title or one title per model (repetition is allowed).
If only one title is given, that is centered over the table.}

\item{dep.varnames}{an optional vector to describe the dependent variable in the models.
Can be either a single variable name or one per model (repetition is allowed).}

\item{k}{an integer to start the counter for model numbers in the table.}

\item{print.xtable.options}{a list of options for \code{\link[xtable]{print.xtable}}.}
}
\description{
This method converts one or more fitted models of class \code{sigfit} into a publication-ready LaTeX table.
This conversion is performed by reformatting the models into a format that is fed to 
\code{\link[xtable]{xtable}}, which generates the actual LaTeX code.
}
\details{
This function produces a ready-to-use LaTeX table for \code{sigfit} objects.
Each column is its own model, along with model-based information.
The generation of LaTeX code is done by \code{\link[xtable]{xtable}} and so additional printing options can
be passed via \code{print.xtable.options}.
For a full list of these options see \code{\link[xtable]{print.xtable}}.
}
\examples{
data("sanctionsData")
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
}

## Using Pseudo Likelihood with user made first stage
Phat <- list(PRhat=sanctionsData$PRhat, PFhat=sanctionsData$PFhat)
fit2 <- sigint(f1, data=sanctionsData, method="pl", phat=Phat)

## Using Pseudo Likelihood with default first stage and bootstrapped standard errors
\dontrun{
fit3 <- sigint(f1, data=sanctionsData, method="pl", pl.vcov=25) 
}

## Simple regression table
toLatexTable(fit2)

## More options: multiple models and user supplied standard errors
\dontrun{
toLatexTable(fit1, fit2, fit3,
        se.list=list(sqrt(diag(vcov(fit1))),
                     sqrt(diag(vcov(fit3))),
                     sqrt(diag(vcov(fit3)))),
        stars="all",
        caption = "Economic Sanctions",
        label = "tab:sanctions",
        model.names = c("NPL", "PL", "PL"))
}         
        
\dontrun{
## More options, from print.xtable including printing to a file
toLatexTable(fit1, fit2, fit3,
        caption = "Economic Sanctions",
        label = "tab:sanctions",
        model.names = c("NPL", "PL", "PL"),
        print.xtable.options=list(file="myTable.tex",
                                  booktabs=TRUE))
}

}
\seealso{
\code{\link[xtable]{xtable}}, \code{\link[xtable]{print.xtable}}
}
