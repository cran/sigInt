##' Economic Sanctions Threats and Outcomes
##' 
##' Dataset on economic sanctions threats and outcomes from 1970-2000
##'
##' These data were compiled using the Threat and Imposition of Sanctions (TIES),
##' data project (Morgan, Bapat, and Kobayashi 2014),
##' with additional data from the Correlates of War (COW),
##' and Polity IV datasets.
##' See Crisman-Cox and Gibilisco (2018) for more information.
##' The unit of
##' observation is the dyad-decade, and the variables are: 
##' \describe{
##' \item{\code{gameID}}{A dyad-decade identifier composed of COW country codes and the decade observed.}
##' \item{\code{dyadID}}{A dyad identifier composed of COW country codes}
##' \item{\code{tenyear}}{The observed decade}
##' \item{\code{code1}}{Challenger's COW code}
##' \item{\code{code2}}{Target's COW code}
##' \item{\code{sq}}{The number of status quo observations in this dyad decade}
##' \item{\code{cd}}{The number of times that the game ends with Challenge-Concede (Outcome \eqn{CD}) }
##' \item{\code{sf}}{The number of times that the game ends with Challenge-Resist
##' -Stand Firm (Outcome \eqn{SF}) }
##' \item{\code{bd}}{The number of times that the game ends with Challenge-Resist-Back Down
##'  (Outcome \eqn{BD}) }
##' \item{\code{senderecondep}}{Challenger's economic dependence (dyadic trade / 
##' Challenger's GDP per capita)}
##' \item{\code{senderdemocracy}}{Challenger's Polity score}
##' \item{\code{contig}}{Contiguity between states}
##' \item{\code{ally}}{Are these state allied? (indicator)}
##' \item{\code{anticipatedsendercosts}}{The Challenger's 
##' anticipated costs for enacting sanctions}
##' \item{\code{anticipatedtargetcosts}}{The Target's 
##' anticipated costs for being sanctions}
##' \item{\code{targetecondep}}{Target's economic dependence (dyadic trade /
##'  Target's GDP per capita)}
##' \item{\code{lncaprat}}{Ratio of the Challenger's military capability to the Target's (logged)}
##' \item{\code{targetdemocracy}}{Target's Polity score}
##' \item{\code{PRhat}}{Estimated probability that the Target resists a challenge (fit using a random forest)}
##' \item{\code{PFhat}}{Estimated probability that the Challenger stands firm given that 
##' it challenged (fit using a random forest)}
##' \item{\code{PRnpl}}{Estimated probability that the Target resists a challenge (taken from the last stage of NPL iteration)}
##' \item{\code{PFnpl}}{Estimated probability that the Challenger stands firm given that 
##' it challenged  (taken from the last stage of NPL iteration)}
##' }
##' @name sanctionsData
##' @usage data(sanctionsData)
##' @docType data
##' @references 
##' Barbieri, Katherine, Omar M. G. Keshk, and Brian Pollins. 2009.
##'  "TRADING DATA: Evaluating our Assumptions and Coding Rules." 
##'  Conflict Management and Peace Science. 26(5): 471-491.
##' 
##' Casey Crisman-Cox and Michael Gibilisco. 2018. "Estimating 
##' Signaling Games in International Relations: Problems and Solutions."
##' Unpublished Manuscript.
##'  
##' Gibler, Douglas M. 2009. International military alliances, 1648-2008. CQ Press.  
##'  
##'  Marshall, Monty G., and Keith Jaggers. 2013. "Polity IV Project." 
##'  \url{http://www.systemicpeace.org/polity/polity4.htm}.
##'  
##' Morgan, T. Clifton, Navin Bapat, and Yoshi Kobayashi. 2014.
##' "The Threat and Imposition of Sanctions: Updating the TIES dataset." 
##' Conflict Management and Peace Science 31(5): 541-558.   
##' 
##' Singer, J. David, Stuart Bremer, and John Stuckey. 1972. 
##' "Capability Distribution, Uncertainty, and Major Power War, 1820-1965."
##'  in Bruce Russett (ed) Peace, War, and Numbers, Beverly Hills: Sage, 19-48.
##'  
##'  Stinnett, Douglas M., Jaroslav Tir, Philip Schafer, Paul F. Diehl, and Charles Gochman. 2002. 
##'  "The Correlates of War Project Direct Contiguity Data, Version 3." 
##'  Conflict Management and Peace Science 19(2):58-66.
##' @seealso \code{\link{sigint}}
##' @keywords data
NULL


