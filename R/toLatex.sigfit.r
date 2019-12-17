#' Export a \code{sigfit} object into paper-ready LaTeX table 
#'
#' This method converts one or more fitted models of class \code{sigfit} into a publication-ready LaTeX table.
#' This conversion is performed by reformatting the models into a format that is fed to 
#' \code{\link[xtable]{xtable}}, which generates the actual LaTeX code.
#' 
#' @param ... one or more models fit using \code{\link{sigint}}.
#' @param se.list an optional list where each element contains the standard errors of the models.
#' If included this list must include one element per model, even if that model's standard errors are unchanged.
#' This argument should be used when standard errors have been adjusted outside the model (e.g., a bootstrap).
#' If left empty the standard errors from the model's covariance matrix are used if available.
#' Within each element of this list, the order of standard errors must be the same order as the coefficients in the model.
#' Standard errors are not matched on name.
#' @param stars how should significance stars be used? \code{stars = "default"} returns a single star when \eqn{p < 0.05},
#' \code{stars = "all"} returns flags for \eqn{p < 0.1}, \eqn{p < 0.05}, and \eqn{p < 0.01}, and \code{stars = "none"} 
#' returns no stars at all.
#' @param caption a string to be used as the table's caption.
#' @param label a string to be used as the table's LaTeX \code{label} argument.
#' @param align a string to indicate the alignment of each column in the table.  
#' Passed directly to the LaTeX \code{tabular} options.
#' @param digits how many digits after the decimal point should be displayed? Default 2.
#' @param se.note a string containing a note for the bottom of the table. 
#' Default is the common "Standard errors in parenthesis.
#' @param order a string vector describing the order that  the covariates should be in the table.
#' The default is to use the order they're listed in the \code{sigfit} object. 
#' When this vector is shorter than the coefficient vector, variables are first included by \code{order},
#' remaining variables are included based on their order in the model
#' @param covariate.labels a string vector of "nice" names for the variables appropriate for a published work.
#' If empty, the "ugly" names from the fitted model are used. Note that if \code{order} is specified then the
#' covariate labels must match the order in \code{order}.
#' @param model.names an optional vector of model names to include as column titles in the table. 
#' Should be either a single title or one title per model (repetition is allowed).
#' If only one title is given, that is centered over the table.  
#' @param dep.varnames an optional vector to describe the dependent variable in the models.
#' Can be either a single variable name or one per model (repetition is allowed).
#' @param k an integer to start the counter for model numbers in the table.
#' @param print.xtable.options a list of options for \code{\link[xtable]{print.xtable}}.
#' @details 
#' This function produces a ready-to-use LaTeX table for \code{sigfit} objects.
#' Each column is its own model, along with model-based information.
#' The generation of LaTeX code is done by \code{\link[xtable]{xtable}} and so additional printing options can
#' be passed via \code{print.xtable.options}.
#' For a full list of these options see \code{\link[xtable]{print.xtable}}.
#' 
#' @seealso \code{\link[xtable]{xtable}}, \code{\link[xtable]{print.xtable}}
#' @examples
#' data("sanctionsData")
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
#' }
#' 
#' ## Using Pseudo Likelihood with user made first stage
#' Phat <- list(PRhat=sanctionsData$PRhat, PFhat=sanctionsData$PFhat)
#' fit2 <- sigint(f1, data=sanctionsData, method="pl", phat=Phat)
#' 
#' ## Using Pseudo Likelihood with default first stage and bootstrapped standard errors
#' \dontrun{
#' fit3 <- sigint(f1, data=sanctionsData, method="pl", pl.vcov=25) 
#' }
#' 
#' ## Simple regression table
#' toLatexTable(fit2)
#'
#' ## More options: multiple models and user supplied standard errors
#' \dontrun{
#' toLatexTable(fit1, fit2, fit3,
#'         se.list=list(sqrt(diag(vcov(fit1))),
#'                      sqrt(diag(vcov(fit3))),
#'                      sqrt(diag(vcov(fit3)))),
#'         stars="all",
#'         caption = "Economic Sanctions",
#'         label = "tab:sanctions",
#'         model.names = c("NPL", "PL", "PL"))
#'}         
#'         
#' \dontrun{
#' ## More options, from print.xtable including printing to a file
#' toLatexTable(fit1, fit2, fit3,
#'         caption = "Economic Sanctions",
#'         label = "tab:sanctions",
#'         model.names = c("NPL", "PL", "PL"),
#'         print.xtable.options=list(file="myTable.tex",
#'                                   booktabs=TRUE))
#'}
#'
#' @export

toLatexTable <- function(..., 
                           se.list, 
                           stars=c("default", "all", "none"),
                           caption="",
                           label, 
                           align, 
                           digits=2,
                           se.note="Standard errors in parenthesis",
                           order,
                           covariate.labels,
                           model.names,
                           dep.varnames,         
                           k=1,
                           print.xtable.options=list()
){
  
  mod.list <- list(...)
  
  print.xtable.default <- list(booktabs=FALSE,
                               sanitize.text.function=function(x){x},
                               include.rownames=FALSE,
                               caption.placement="top",
                               table.placement="h!",
                               include.colnames =FALSE)
  if(!is.null(print.xtable.options$tabular.environment) && print.xtable.options$tabular.environment=="longtable"){
    print.xtable.default$hline.after = c(-1,0)
  }
  print.xtable.options <- modifyList(print.xtable.default, print.xtable.options)
  if(missing(label)){
    label=paste("tab:",k, sep="")
  }

    coef.list <- lapply(mod.list, function(x){x$coef})    
    coef.list <- lapply(coef.list, 
                        function(x){
                            names(x) <- gsub(x=names(x), pattern = "_", replacement = ".")
                            return(x)
                        }
                        )
    
    
  stars <- match.arg(stars)
  if(missing(se.list)){
    se.list <- lapply(mod.list, 
                      function(x){
                        if(is.null(x$vcov)){
                          se <- NA*x$coef
                        }else{
                          se <- sqrt(diag(x$vcov))
                        }
                        se[is.nan(se)] <- NA
                        names(se) <- names(x$coef)
                        return(se)
                      }
    )
    
  }else{
      if(length(se.list) != length(mod.list)){
          stop("se.list must have one element per model")
      }
  }
    
  
  if(missing(model.names)){model.names <- NULL}
  if(missing(dep.varnames)){dep.varnames <- NULL}
  length <- 2*length(unique(unlist(lapply(coef.list,names))))
  if(missing(order)){
    order <- unique(unlist(lapply(coef.list, names)))
  }else{
    order <- unique(c(order, unlist(lapply(coef.list, names))))
  }
  order <- rep(order, each=2)
  order[seq(2, length(order), by=2)] <- paste(order[seq(2, length(order), by=2)], 
                                              "_SE", 
                                              sep="")
  
  output<-data.frame()
  for(i in 1:length(mod.list)){
    table <- coeftest2(coef.list[[i]], 
                       se.list[[i]])
    
    
    a.table<-data.frame()
    for(j in 1:nrow(table)){
      temp<-table[j, 1:2]
      temp<-t(t(temp))
      temp<-num2str(temp, digits=digits)
      temp<-matrix(stringr::str_trim(temp), nrow=nrow(temp))
      temp[2,]<- ifelse(temp[2,] != "NA" ,
                        paste("(", temp[2,], ")", sep=""),
                        "--")
      
      if(stars=="default"){
        if(!is.na(table[j, 4]) && table[j, 4]<0.05){
          temp[1,]<- paste(temp[1,] ,"$^*$", sep="")
        }
      }else{
        if(stars=="all"){
          if(!is.na(table[j, 4]) &&  (table[j, 4]<0.1 & table[j, 4]>0.05)){
            temp[1,]<- paste(temp[1,] 
                             ,"$^\\dagger$", sep="")
          }
          if(!is.na(table[j, 4]) &&  (table[j, 4]<0.05 & table[j, 4]>0.01)){
            temp[1,]<- paste(temp[1,],
                             "$^*$", sep="")
          }
          if(!is.na(table[j, 4]) &&  table[j, 4]<0.01){
            temp[1,]<- paste(temp[1,],
                             "$^{**}$", sep="")
          }
        }
      }
      row.names(temp)<-NULL
      var.name<-row.names(table)[j]
      colnames(temp)<-paste("Model.", i, sep="")
      row.names(temp)<-c(var.name, paste(var.name, "_SE", sep=""))
      a.table<-rbind(a.table, temp)  
    }
    if(i==1){
      a.output<-a.table
      ` `<-row.names(a.output)
      row.names(a.output)<-NULL
      a.output<-cbind(` `, a.output)
      a.output$` ` <- factor(a.output$` `, 
                               levels=order)
    }else{
      ` `<-row.names(a.table)
      row.names(a.table)<-NULL
      a.table<-cbind(` `, a.table) 
      a.table$` ` <- factor(a.table$` `, 
                            levels=order)
      
      
      a.output<-merge(a.output, 
                      a.table,
                      by= " ",
                      all=TRUE,
                      sort=TRUE)        
      a.output <- a.output[order(a.output[,1]),]
    }
    
  }
  a.output$` ` <- as.character(a.output$` `)
  a.output$` `[stringr::str_detect(a.output$` `, "_SE")]<-""
  output <- a.output
  
    
  
  if(!missing(covariate.labels)){
    output[output[,1]!="",1] <- covariate.labels
  }
  if(missing(align)){
    align <- stringr::str_c(c("rr", rep("c", ncol(output)-1)), collapse="")
  }else{
    align <- stringr::str_c("r", align, collapse="")
  }
 length <- nrow(output)
  info <- mod.sum(mod.list, stars, se.note, 
                  model.names, dep.varnames, k, length,
                  print.xtable.options$booktabs)
  
  print.xtable.options$add.to.row$pos <- append(print.xtable.options$add.to.row$pos,
                                                list(info[[3]][1],info[[3]][2]))
  print.xtable.options$add.to.row$command <- c(print.xtable.options$add.to.row$command,
                                               c(info[[1]], info[[2]]))

  suppressWarnings(
    x <- xtable::xtable(output,
                caption=caption,
                label=label, 
                align=align, 
                digits=digits))
  print.xtable.options <- modifyList(print.xtable.options, list(x=x))
  do.call(print, print.xtable.options)
  
  if(stringr::str_detect(align, "S")){
    cat("\nSI unit detected, remember to add the following to your TeX preamble:\n")
    cat("\n\\usepackage{siunitx}\n\\sisetup{\n\tinput-symbols=(),\n\ttable-align-text-post = false,\n\tgroup-digits=false,\n} ")
  }
}
num2str <- function(x, digits=2){formatC(x, digits=digits, format='f')}

mod.sum<-function(mod.list, 
                  stars,
                  seNote,
                  model.names,
                  dep.varnames,
                  k,
                  length,
                  booktabs){  ##Function to add in LogLik and N to xtable
  model.numbers <- paste("\\multicolumn{1}{c}{Model ",
                         k:(k+length(mod.list)-1),
                         "}",
                         sep="")
  model.numbers <- stringr::str_c(model.numbers,  collapse=" & " )
  model.numbers <- paste("&", model.numbers, "\\\\ \n")
  
  L.line <- ifelse(booktabs, "\\midrule", "\\hline")
  N.line <- ifelse(booktabs, "\\bottomrule", "\\hline")
  out.L<- paste(L.line, "\n Log $L$")
  out.N<-c("$N$ Games")
  modNums <- " "
  nCol <- length(mod.list)
  if(!is.null(dep.varnames)){ 
    if(length(dep.varnames)==1){
      dep.varnames <- paste("& \\multicolumn{", nCol, "}{c}{", dep.varnames, "}\\\\ \n")
    } else {
      if(length(dep.varnames)==nCol){
        counter <- 1; dep.varnames.out <- "& "
        for(j in 1:(nCol-1)){
          if(dep.varnames[j] == dep.varnames[j+1]){
            counter <- counter+1
          }else{
            dep.varnames.out <- stringr::str_c(dep.varnames.out,
                                               paste("\\multicolumn{", counter,"}{c}{", dep.varnames[j], "}& ", sep=""),
                                               collapse="")
            counter <- 1
          }
        }
        dep.varnames.out <- stringr::str_c(dep.varnames.out,
                                           paste("\\multicolumn{", counter,"}{c}{", dep.varnames[j+1], 
                                                 "}\\\\ \n", sep=""),
                                           collapse="")
        dep.varnames <- dep.varnames.out
      }else{
        warning("Unclear how to match dep.varnames, dropping this row")
        dep.varnames <- NULL
      }
    }
  }
  if(!is.null(model.names)){ 
    if(length(model.names)==1){
      model.names <- paste("& \\multicolumn{", nCol, "}{c}{", model.names, "}\\\\ \n")
    } else {
      if(length(model.names)==nCol){
        counter <- 1; model.names.out <- "& "
        for(j in 1:(nCol-1)){
          if(model.names[j] == model.names[j+1]){
            counter <- counter+1
          }else{
            model.names.out <- stringr::str_c(model.names.out,
                                               paste("\\multicolumn{", counter,"}{c}{", model.names[j], "}& ", sep=""),
                                               collapse="")
            counter <- 1
          }
        }
        model.names.out <- stringr::str_c(model.names.out,
                                           paste("\\multicolumn{", counter,"}{c}{", model.names[j+1], 
                                                 "}\\\\ \n", sep=""),
                                           collapse="")
        model.names <- model.names.out
      }else{
        warning("Unclear how to match model.names, dropping this row")
        model.names <- NULL
      }
    }
  }
  
  
  for(i in 1:length(mod.list)){
    mod.i<-mod.list[[i]]
    out.L <- stringr::str_c(out.L,
                            paste(" & ",
                                  "\\multicolumn{1}{c}{",
                                  num2str(mod.i$logLik),
                                  "}"))
    out.N <- stringr::str_c(out.N,
                            paste(" & ",
                                  "\\multicolumn{1}{c}{", 
                                  nrow(mod.i$model),
                                  "}"))
  }
  
  
  out <- stringr::str_c(out.L, "\\\\\n", 
                        out.N,   "\\\\", N.line,"\n")
  
  footnote  <- ifelse(stars=="all",
                      paste("\\multicolumn{", 
                            1+length(mod.list),
                            "}{l}{\\footnotesize{\\emph{Notes:} $^{**}p<0.01$; $^{*}p<0.05$; $^\\dagger p<0.1$ }} \\\\ \n\\multicolumn{",
                            1+length(mod.list),
                            "}{l}{\\footnotesize{", 
                            seNote,
                            "}}%",
                            sep=""),
                      ifelse(stars=="default",
                             paste("\\multicolumn{", 
                                   1+length(mod.list),
                                   "}{l}{\\footnotesize{\\emph{Notes:} $^{*}p<0.05$}}\\\\ \n\\multicolumn{",
                                   1+length(mod.list),
                                   "}{l}{\\footnotesize{", 
                                   seNote,
                                   "}}%",
                                   sep=""),
                             paste("\\multicolumn{", 
                                   1+length(mod.list),
                                   "}{l}{\\footnotesize{\\emph{Notes:} ",
                                   seNote,
                                   "}}%",
                                   sep="")))
  
  
  header <- stringr::str_c(c(dep.varnames, model.names, model.numbers), collapse="")
  
  
  
  out <- list(header=header,
              out=stringr::str_c(out, footnote), 
              length=c(0,length))
  
  return(out)
}

coeftest2<-function(beta, sd){
  z <- beta/sd
  p.value <- pnorm(abs(z), lower.tail=FALSE)*2
  results <- cbind(beta, sd,z, p.value)
  colnames(results) <- c("Estimate",
                         "Std. Error",
                         "Z Value",
                         "Pr(|z|)")
  return(results)
}
