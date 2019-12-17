############################################
###File to calculate the 1st derivatives ###
###CMLE FUNCTIONS ARE NOW HOUSED IN      ###
###cmleRegressFunctions3.r               ###
############################################

# casey, I changed this to remove dividing by 0 problems --- MBG

dbivnorm <- function(x1, x2, rho = 0){
  denom <- 1 - rho^2
  pdf <- 1/(2*pi * sqrt(denom)) * exp( -1/(2*denom)  * (x1^2 + x2^2 - 2 *rho * x1 *x2))
  return(pdf)
}

delPbivnorm <- function(x1, x2, rho=0){
  ## The partial derivative is from Wickens (1992) for 
  ## $\pnorm_2(x, y, rho)$ we have
  ## $Dx = \pnorm(x) * \pnorm((y-x*rho)/sqrt(1-rho^2))$
  ## This function always just put the one that is w.r.t.
  ## as the first argument.
  return(dnorm(x1)*pnorm((x2-x1*rho)/sqrt(1-rho^2)))
}


eval_gr_qll <- function(x,PRhat,PFhat,Y,regr, fixed.par){
  # here x = (theta,p)
  ##EQ[,1]= pR
  ##EQ[,2]= pC
  ##EQ[,3]= pF
  M <- dim(Y)[2]
  
  
  param <-vec2U.regr(x,regr, fixed.par)
  
  
  c <- cStar.jo(PRhat,param)
  pR <- f.jo(PFhat, param)
  pR[pR<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) 
  pC <- g.jo(c, param)
  pC[pC<=(.Machine$double.eps)^(.25)] <- (.Machine$double.eps)^(.25) #edited 7 March
  pC[pC>=1-(.Machine$double.eps)^(.25)] <- 1-(.Machine$double.eps)^(.25) #edited 7 March
  pF <- h.jo(c, param)/pC
  v1 <- (c-param$barWA)/param$sig
  v2 <- (c-param$bara)/param$sig
  d1 <- (param$barWA - param$bara)/(param$sig*sqrt(2))
  d2 <- (param$barWA - c)/(param$sig)
  
  r <- 1/sqrt(2)
  #pRratio <- (pR-1)/pR
  VApr <- (param$VA/pR - param$VA * (pR-1)/(pR^2))
  P1 <- pnorm(v1)
  P1[P1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P1[P1>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  P2 <- pnorm(v2)
  P2[P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P2[P2>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  D1 <- dnorm(v1)
  # D1[D1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  D2 <- dnorm(v2)
  # D2[D2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P1P2 <-  P1*P2 #more stable than pbiv w/rho=0
  P1P2[P1P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  Del_d1v1 <- delPbivnorm(d1, -v1, r)
  Del_v1d1 <- delPbivnorm(-v1, d1, r)
  PBdv <- pbivnorm(d1, -v1, rho=r)
  PBdv[PBdv<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) # LOOK AT THIS
  PBdv[PBdv>= (1-sqrt(.Machine$double.eps))] <- 1-sqrt(.Machine$double.eps) # LOOK AT THIS 
  WBinner <- as.numeric((-param$CB + param$VB*(-PFhat + 1) + param$barWB*PFhat)/PFhat)
  
  invpWBinner <- (pnorm(WBinner, lower.tail=F) )
  invpWBinner[invpWBinner<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) # LOOK AT THIS
  invpWBinner[invpWBinner>= (1-sqrt(.Machine$double.eps))] <- 1-sqrt(.Machine$double.eps) # LOOK AT THIS

  pWBinner <- (pnorm(WBinner) )
  pWBinner[pWBinner<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) # LOOK AT THIS
  pWBinner[pWBinner>= (1-sqrt(.Machine$double.eps))] <- 1-sqrt(.Machine$double.eps) # LOOK AT THIS

  dCB.denom <- (PFhat*(pnorm(WBinner, lower.tail=F))) 
  dCB.denom[dCB.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dCB.denom2 <- (PFhat*pWBinner)
  dCB.denom2[dCB.denom2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  
  dWA4.denom <- (pC*(1 - PBdv/pC)*pWBinner)  # LOOK AT THIS	
  dWA4.denom[dWA4.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) # AND THIS
  dWA <- rbind( -D1/P1,
                P2*D1/pC,
                (sqrt(2)*Del_d1v1/2 +Del_v1d1)/PBdv,
                (pC*((-sqrt(2)*Del_d1v1/2 - Del_v1d1)/pC + P2*PBdv*D1/pC**2)*pWBinner + (1 - PBdv/pC)*P2*pWBinner*D1)/dWA4.denom  # DITTO THIS!!!
  )
  
  dWA  <- apply(regr$barWA, 2, function(x){t(as.numeric(x)*t(dWA))*as.numeric(Y)})
  
  dWB.denom <- invpWBinner
  dWB.denom[dWB.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dWB.denom2 <- pWBinner
  dWB.denom2[dWB.denom2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dWB <- rbind(0,
               (-dnorm(WBinner)/dWB.denom),
               ( dnorm(WBinner)/dWB.denom2),
               ( dnorm(WBinner)/dWB.denom2)
  ) 
  
  dWB <- apply(regr$barWB, 2, function(x){t(as.numeric(x)*t(dWB))*as.numeric(Y)})
  
  
  dbara.denom <- (pC*(-pF + 1)*pWBinner)
  dbara.denom[dbara.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dbara.stuff <- (-2*P1*P2 + 2)
  dbara.stuff[dbara.stuff<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dbara <- rbind( -D2/P2,
                  P1*D2/pC,
                  - sqrt(2)*Del_d1v1/(2*PBdv),
                  ((pC*(sqrt(2)*Del_d1v1/dbara.stuff + pF*P1*D2/pC)*pWBinner + (-pF + 1)*P1*pWBinner*D2)/(dbara.denom))
                  
  )
  dbara <- apply(regr$bara, 2, function(x){t(as.numeric(x)*t(dbara))*as.numeric(Y)})
  
  dVA.denom <- (pC*(-pF + 1)*pWBinner)
  dVA.denom[dVA.denom<=.Machine$double.eps] <- .Machine$double.eps
  dVA <- rbind( ((PRhat - 1)*P1*D2/PRhat + (PRhat - 1)*P2*D1/PRhat)/(P1*P2),
                (-(PRhat - 1)*P1*D2/PRhat - (PRhat - 1)*P2*D1/PRhat)/pC,
                - (PRhat - 1)*Del_v1d1/(PRhat*PBdv),
                ((pC*(-pF*((PRhat - 1)*P1*D2/PRhat + (PRhat - 1)*P2*D1/PRhat)/pC + (PRhat - 1)*Del_v1d1/(PRhat*pC))*pWBinner + (-pF + 1)*(-(PRhat - 1)*P1*D2/PRhat - (PRhat - 1)*P2*D1/PRhat)*pWBinner)/(dVA.denom))
                
  )
  dVA  <- apply(regr$VA, 2, function(x){t(as.numeric(x)*t(dVA))*as.numeric(Y)})
  
  dVB <- rbind(0,
               -(-PFhat + 1)*dnorm(WBinner)/(dCB.denom),
               (-PFhat + 1)*dnorm(WBinner)/(dCB.denom2),
               (-PFhat + 1)*dnorm(WBinner)/(dCB.denom2)
  )  
  dVB  <- apply(regr$VB, 2, function(x){t(as.numeric(x)*t(dVB))*as.numeric(Y)})
  
  # dSA.denom <- (pC*(1 - PBdv/pC)*pWBinner)
  # dSA.denom[dSA.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dSA.denom2 <- (PRhat*(-pC+PBdv))
  dSA.denom2[abs(dSA.denom2)<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  dSA <- rbind( ((P1*D2 + P2*D1)/PRhat)/(P1P2),
                ((-P1*D2- P2*D1)/PRhat)/pC ,
                -Del_v1d1 /(PRhat*PBdv),
                (P1*D2 + P2*D1 - Del_v1d1)/dSA.denom2
  )
  dSA  <- apply(regr$SA, 2, function(x){t(as.numeric(x)*t(dSA))*as.numeric(Y)})
  

  dCB <- rbind(0,
               dnorm(WBinner)/dCB.denom,
               -dnorm(WBinner)/dCB.denom2,
               - dnorm(WBinner)/dCB.denom2
  )
  # dCB[2:3,][dnorm(WBinner) < 1e-5] <- 1e-5
  dCB  <- apply(regr$CB, 2, function(x){t(as.numeric(x)*t(dCB))*as.numeric(Y)})
  
  
  out <-  -cbind(dSA, dVA, dCB, dWA, dWB, dbara, dVB)
  return(colSums(out)) 
}









eval_gr_qll.i <- function(x,PRhat,PFhat,Y,regr, fixed.par){
  # here x = (theta,p)
  ##EQ[,1]= pR
  ##EQ[,2]= pC
  ##EQ[,3]= pF
  M <- dim(Y)[2]
  
  
  param <-vec2U.regr(x,regr, fixed.par)
  
  
  c <- cStar.jo(PRhat,param)
  pR <- f.jo(PFhat, param)
  pR[pR<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) 
  pC <- g.jo(c, param)
  pC[pC<=(.Machine$double.eps)^(.25)] <- (.Machine$double.eps)^(.25) #edited 7 March
  pC[pC>=1-(.Machine$double.eps)^(.25)] <- 1-(.Machine$double.eps)^(.25) #edited 7 March
  pF <- h.jo(c, param)/pC
  v1 <- (c-param$barWA)/param$sig
  v2 <- (c-param$bara)/param$sig
  d1 <- (param$barWA - param$bara)/(param$sig*sqrt(2))
  d2 <- (param$barWA - c)/(param$sig)
  
  r <- 1/sqrt(2)
  #pRratio <- (pR-1)/pR
  VApr <- (param$VA/pR - param$VA * (pR-1)/(pR^2))
  P1 <- pnorm(v1)
  P1[P1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P1[P1>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  P2 <- pnorm(v2)
  P2[P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P2[P2>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  D1 <- dnorm(v1)
  # D1[D1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  D2 <- dnorm(v2)
  # D2[D2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P1P2 <-  P1*P2 #more stable than pbiv w/rho=0
  P1P2[P1P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  
  Del_d1v1 <- delPbivnorm(d1, -v1, r)
  Del_v1d1 <- delPbivnorm(-v1, d1, r)
  PBdv <- pbivnorm(d1, -v1, rho=r)
  PBdv[PBdv<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) # LOOK AT THIS
  PBdv[PBdv>= (1-sqrt(.Machine$double.eps))] <- 1-sqrt(.Machine$double.eps) # LOOK AT THIS 
  WBinner <- as.numeric((-param$CB + param$VB*(-PFhat + 1) + param$barWB*PFhat)/PFhat)
  
  invpWBinner <- (pnorm(WBinner, lower.tail=F) )
  invpWBinner[invpWBinner<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) # LOOK AT THIS
  invpWBinner[invpWBinner>= (1-sqrt(.Machine$double.eps))] <- 1-sqrt(.Machine$double.eps) # LOOK AT THIS 
  
  pWBinner <- (pnorm(WBinner) )
  pWBinner[pWBinner<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) # LOOK AT THIS
  pWBinner[pWBinner>= (1-sqrt(.Machine$double.eps))] <- 1-sqrt(.Machine$double.eps) # LOOK AT THIS 
  
  
  dWA4.denom <- (pC*(1 - PBdv/pC)*pWBinner)  # LOOK AT THIS	
  dWA4.denom[dWA4.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) # AND THIS
  dWA <- rbind( -D1/P1,
                P2*D1/pC,
                (sqrt(2)*Del_d1v1/2 +Del_v1d1)/PBdv,
                (pC*((-sqrt(2)*Del_d1v1/2 - Del_v1d1)/pC + P2*PBdv*D1/pC**2)*pWBinner + (1 - PBdv/pC)*P2*pWBinner*D1)/dWA4.denom  # DITTO THIS!!!
  )
  
  dWA  <- apply(regr$barWA, 2, function(x){colSums(t(as.numeric(x)*t(dWA))*as.numeric(Y))})
  
  dWB.denom <- invpWBinner
  dWB.denom[dWB.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dWB.denom2 <- pWBinner
  dWB.denom2[dWB.denom2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dWB <- rbind(0,
               (-dnorm(WBinner)/dWB.denom),
               ( dnorm(WBinner)/dWB.denom2),
               ( dnorm(WBinner)/dWB.denom2)
  ) 
  
  dWB <- apply(regr$barWB, 2, function(x){colSums(t(as.numeric(x)*t(dWB))*as.numeric(Y))})
  
  
  dbara.denom <- (pC*(-pF + 1)*pWBinner)
  dbara.denom[dbara.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dbara.stuff <- (-2*P1*P2 + 2)
  dbara.stuff[dbara.stuff<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dbara <- rbind( -D2/P2,
                  P1*D2/pC,
                  - sqrt(2)*Del_d1v1/(2*PBdv),
                  ((pC*(sqrt(2)*Del_d1v1/dbara.stuff + pF*P1*D2/pC)*pWBinner + (-pF + 1)*P1*pWBinner*D2)/(dbara.denom))
                  
  )
  dbara <- apply(regr$bara, 2, function(x){colSums(t(as.numeric(x)*t(dbara))*as.numeric(Y))})
  
  dVA.denom <- (pC*(-pF + 1)*pWBinner)
  dVA.denom[dVA.denom<=.Machine$double.eps] <- .Machine$double.eps
  dVA <- rbind( ((PRhat - 1)*P1*D2/PRhat + (PRhat - 1)*P2*D1/PRhat)/(P1*P2),
                (-(PRhat - 1)*P1*D2/PRhat - (PRhat - 1)*P2*D1/PRhat)/pC,
                - (PRhat - 1)*Del_v1d1/(PRhat*PBdv),
                ((pC*(-pF*((PRhat - 1)*P1*D2/PRhat + (PRhat - 1)*P2*D1/PRhat)/pC + (PRhat - 1)*Del_v1d1/(PRhat*pC))*pWBinner + (-pF + 1)*(-(PRhat - 1)*P1*D2/PRhat - (PRhat - 1)*P2*D1/PRhat)*pWBinner)/(dVA.denom))
                
  )
  dVA  <- apply(regr$VA, 2, function(x){colSums(t(as.numeric(x)*t(dVA))*as.numeric(Y))})
  
  dVB <- rbind(0,
               -(-PFhat + 1)*dnorm(WBinner)/(PFhat*invpWBinner),
               (-PFhat + 1)*dnorm(WBinner)/(PFhat*pWBinner),
               (-PFhat + 1)*dnorm(WBinner)/(PFhat*pWBinner)
  )  
  dVB  <- apply(regr$VB, 2, function(x){colSums(t(as.numeric(x)*t(dVB))*as.numeric(Y))})
  
  
  
  # dSA.denom <- (pC*(1 - PBdv/pC)*pWBinner)
  # dSA.denom[dSA.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dSA.denom2 <- (PRhat*(-pC+PBdv))
  dSA.denom2[abs(dSA.denom2)<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  dSA <- rbind( ((P1*D2 + P2*D1)/PRhat)/(P1P2),
                ((-P1*D2- P2*D1)/PRhat)/pC ,
                -Del_v1d1 /(PRhat*PBdv),
                (P1*D2 + P2*D1 - Del_v1d1)/dSA.denom2
  )
  dSA  <- apply(regr$SA, 2, function(x){colSums(t(as.numeric(x)*t(dSA))*as.numeric(Y))})
  
  dCB.denom <- (PFhat*(pnorm(WBinner, lower.tail=F))) 
  dCB.denom[dCB.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dCB.denom2 <- (PFhat*pWBinner)
  dCB.denom2[dCB.denom2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dCB <- rbind(0,
               dnorm(WBinner)/dCB.denom,
               -dnorm(WBinner)/dCB.denom2,
               - dnorm(WBinner)/dCB.denom2
  )
  # dCB[2:3,][dnorm(WBinner) < 1e-5] <- 1e-5
  dCB  <- apply(regr$CB, 2, function(x){colSums(t(as.numeric(x)*t(dCB))*as.numeric(Y))})
  
  
  out <-  -cbind(dSA, dVA, dCB, dWA, dWB, dbara, dVB)
  return(out)
}








eval_gr_qll.ip <- function(PRhat, PFhat, x, Y,regr, fixed.par){
  # here x = (theta,p)
  ##EQ[,1]= pR
  ##EQ[,2]= pC
  ##EQ[,3]= pF
  M <- dim(Y)[2]
  
  param <-vec2U.regr(x,regr, fixed.par)
  
  
  c <- cStar.jo(PRhat,param)
  pR <- f.jo(PFhat, param)
  pR[pR<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) 
  pC <- g.jo(c, param)
  pC[pC<=(.Machine$double.eps)^(.25)] <- (.Machine$double.eps)^(.25) #edited 7 March
  pC[pC>=1-(.Machine$double.eps)^(.25)] <- 1-(.Machine$double.eps)^(.25) #edited 7 March
  pF <- h.jo(c, param)/pC
  v1 <- (c-param$barWA)/param$sig
  v2 <- (c-param$bara)/param$sig
  d1 <- (param$barWA - param$bara)/(param$sig*sqrt(2))
  d2 <- (param$barWA - c)/(param$sig)
  
  r <- 1/sqrt(2)
  #pRratio <- (pR-1)/pR
  VApr <- (param$VA/pR - param$VA * (pR-1)/(pR^2))
  P1 <- pnorm(v1)
  P1[P1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P1[P1>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  P2 <- pnorm(v2)
  P2[P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P2[P2>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  D1 <- dnorm(v1)
  # D1[D1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  D2 <- dnorm(v2)
  # D2[D2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P1P2 <-  P1*P2 #more stable than pbiv w/rho=0
  P1P2[P1P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  Del_d1v1 <- delPbivnorm(d1, -v1, r)
  Del_v1d1 <- delPbivnorm(-v1, d1, r)
  PBdv <- pbivnorm(d1, -v1, rho=r)
  PBdv[PBdv<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) # LOOK AT THIS
  PBdv[PBdv>= (1-sqrt(.Machine$double.eps))] <- 1-sqrt(.Machine$double.eps) # LOOK AT THIS 
  WBinner <- as.numeric((-param$CB + param$VB*(-PFhat + 1) + param$barWB*PFhat)/PFhat)
  invpWBinner <- (pnorm(WBinner, lower.tail=F) )
  invpWBinner[invpWBinner<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) # LOOK AT THIS
  invpWBinner[invpWBinner>= (1-sqrt(.Machine$double.eps))] <- 1-sqrt(.Machine$double.eps) # LOOK AT THIS 
  
  pWBinner <- (pnorm(WBinner) )
  pWBinner[pWBinner<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) # LOOK AT THIS
  pWBinner[pWBinner>= (1-sqrt(.Machine$double.eps))] <- 1-sqrt(.Machine$double.eps) # LOOK AT THIS 
  
  dPRhat  <- rbind(
    (P1*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D2 + P2*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D1)/(P1P2) ,
    (-P1*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D2 - P2*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D1)/pC,
    (-param$VA/PRhat + (param$SA - param$VA*(-PRhat + 1))/PRhat**2)*Del_v1d1/PBdv,
    (pC*(-(-param$VA/PRhat + (param$SA - param$VA*(-PRhat + 1))/PRhat**2)*Del_v1d1/pC - (P1*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D2+ P2*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D1)*PBdv/pC**2)*pWBinner + (1 - PBdv/pC)*(-P1*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D2- P2*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D1)*pWBinner)/(pC*(1 - PBdv/pC)*pWBinner) 
  )
  dPRhat <- diag(colSums(dPRhat*Y)) #should we make this sparse for bigger datasets?
  
  dPFhat  <- rbind(
    0,
    -((-param$VB +param$barWB)/PFhat + (param$CB - param$VB*(-PFhat + 1) -param$barWB*PFhat)/PFhat**2)*dnorm(WBinner)/invpWBinner,
    ((-param$VB + param$barWB)/PFhat + (param$CB - param$VB*(-PFhat + 1) - param$barWB*PFhat)/PFhat**2)*dnorm(WBinner)/pWBinner,
    ((-param$VB + param$barWB)/PFhat + (param$CB - param$VB*(-PFhat + 1) - param$barWB*PFhat)/PFhat**2)*dnorm(WBinner)/pWBinner
  )
  dPFhat <- diag(colSums(dPFhat*Y)) #should we make this sparse for bigger datasets?
  
  
  return(-cbind(dPRhat, dPFhat))
}


dPsiDp <- function(PRhat, PFhat, x, Y,regr, fixed.par){
  # here x = (theta,p)
  ##EQ[,1]= pR
  ##EQ[,2]= pC
  ##EQ[,3]= pF
  M <- dim(Y)[2]
  
  param <-vec2U.regr(x,regr, fixed.par)
  
  
  c <- cStar.jo(PRhat,param)
  pR <- f.jo(PFhat, param)
  pR[pR<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) 
  pC <- g.jo(c, param)
  pC[pC<=(.Machine$double.eps)^(.25)] <- (.Machine$double.eps)^(.25) #edited 7 March
  pC[pC>=1-(.Machine$double.eps)^(.25)] <- 1-(.Machine$double.eps)^(.25) #edited 7 March
  pF <- h.jo(c, param)/pC
  v1 <- (c-param$barWA)/param$sig
  v2 <- (c-param$bara)/param$sig
  d1 <- (param$barWA - param$bara)/(param$sig*sqrt(2))
  d2 <- (param$barWA - c)/(param$sig)
  
  r <- 1/sqrt(2)
  #pRratio <- (pR-1)/pR
  VApr <- (param$VA/pR - param$VA * (pR-1)/(pR^2))
  P1 <- pnorm(v1)
  P1[P1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P1[P1>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  P2 <- pnorm(v2)
  P2[P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P2[P2>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  D1 <- dnorm(v1)
  # D1[D1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  D2 <- dnorm(v2)
  # D2[D2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P1P2 <-  P1*P2 #more stable than pbiv w/rho=0
  P1P2[P1P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  Del_d1v1 <- delPbivnorm(d1, -v1, r)
  Del_v1d1 <- delPbivnorm(-v1, d1, r)
  PBdv <- pbivnorm(d1, -v1, rho=r)
  PBdv[PBdv<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) # LOOK AT THIS
  PBdv[PBdv>= (1-sqrt(.Machine$double.eps))] <- 1-sqrt(.Machine$double.eps) # LOOK AT THIS 
  WBinner <- as.numeric((-param$CB + param$VB*(-PFhat + 1) + param$barWB*PFhat)/PFhat)
  
  
  dPRhat <- ((-param$VB + param$barWB)/PFhat + (param$CB - param$VB*(-PFhat + 1) - param$barWB*PFhat)/PFhat**2)*dnorm(WBinner)
  dPRhat <- diag(dPRhat) #should we make this sparse for bigger datasets?
  
  dPFhat  <- (-param$VA/PRhat + (param$SA - param$VA*(-PRhat + 1))/PRhat**2)*Del_v1d1/pC + (P1*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D2 + P2*(param$VA/PRhat + (-param$SA + param$VA*(-PRhat + 1))/PRhat**2)*D1)*PBdv/pC**2
  dPFhat <- diag(dPFhat) #should we make this sparse for bigger datasets?
  ZERO <- matrix(0, nrow=M, ncol=M)
  
  return( rbind(cbind(ZERO, dPRhat), cbind(dPFhat, ZERO)) ) #sparsity gains
}




dPsi.dTheta <- function(x, PRhat, PFhat, Y, regr, fixed.par){
  # here x = (theta,p)
  ##EQ[,1]= pR
  ##EQ[,2]= pC
  ##EQ[,3]= pF
  M <- dim(Y)[2]
  
  param <-vec2U.regr(x,regr, fixed.par)
  
  
  c <- cStar.jo(PRhat,param)
  pR <- f.jo(PFhat, param)
  pR[pR<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) 
  pC <- g.jo(c, param)
  pC[pC<=(.Machine$double.eps)^(.25)] <- (.Machine$double.eps)^(.25) #edited 7 March
  pC[pC>=1-(.Machine$double.eps)^(.25)] <- 1-(.Machine$double.eps)^(.25) #edited 7 March
  pF <- h.jo(c, param)/pC
  v1 <- (c-param$barWA)/param$sig
  v2 <- (c-param$bara)/param$sig
  d1 <- (param$barWA - param$bara)/(param$sig*sqrt(2))
  d2 <- (param$barWA - c)/(param$sig)
  
  r <- 1/sqrt(2)
  #pRratio <- (pR-1)/pR
  VApr <- (param$VA/pR - param$VA * (pR-1)/(pR^2))
  P1 <- pnorm(v1)
  P1[P1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P1[P1>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  P2 <- pnorm(v2)
  P2[P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P2[P2>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  D1 <- dnorm(v1)
  # D1[D1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  D2 <- dnorm(v2)
  # D2[D2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P1P2 <-  P1*P2 #more stable than pbiv w/rho=0
  P1P2[P1P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  Del_d1v1 <- delPbivnorm(d1, -v1, r)
  Del_v1d1 <- delPbivnorm(-v1, d1, r)
  PBdv <- pbivnorm(d1, -v1, rho=r)
  PBdv[PBdv<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) # LOOK AT THIS
  PBdv[PBdv>= (1-sqrt(.Machine$double.eps))] <- 1-sqrt(.Machine$double.eps) # LOOK AT THIS 
  WBinner <- as.numeric((-param$CB + param$VB*(-PFhat + 1) + param$barWB*PFhat)/PFhat)
  
  
  # -cbind(dSA, dVA, dCB, dWA, dWB, dbara, dVB)
  dPRdSA <- matrix(0, nrow=M, ncol=ncol(regr$SA))
  dPRdVA <- matrix(0, nrow=M, ncol=ncol(regr$VA))
  dPRda <- matrix(0, nrow=M, ncol=ncol(regr$bara))
  dPRdWA <- matrix(0, nrow=M, ncol=ncol(regr$barWA))
  
  dPRdCB <- -dnorm(WBinner)/PFhat
  dPRdCB  <- apply(regr$CB, 2, function(x){t(as.numeric(x)*t(dPRdCB))})
  dPRdWB <- dnorm(WBinner)
  dPRdWB  <- apply(regr$barWB, 2, function(x){t(as.numeric(x)*t(dPRdWB))})
  dPRdVB <- (-PFhat + 1)*dnorm(WBinner)/PFhat
  dPRdVB  <- apply(regr$VB, 2, function(x){t(as.numeric(x)*t(dPRdVB))})
  
  dPFdCB <- matrix(0, nrow=M, ncol=ncol(regr$CB))
  dPFdWB <- matrix(0, nrow=M, ncol=ncol(regr$barWB))
  dPFdVB <- matrix(0, nrow=M, ncol=ncol(regr$VB))
  
  dPFdSA <- (P1*D2/PRhat + P2*D1/PRhat)*PBdv/pC**2 - Del_v1d1/(PRhat*pC)
  dPFdSA  <- apply(regr$SA, 2, function(x){t(as.numeric(x)*t(dPFdSA))})
  dPFdVA <- (P1*(PRhat - 1)*D2/PRhat + P2*(PRhat - 1)*D1/PRhat)*PBdv/pC**2 - (PRhat - 1)*Del_v1d1/(PRhat*pC)
  dPFdVA  <- apply(regr$VA, 2, function(x){t(as.numeric(x)*t(dPFdVA))})
  dPFdWA <- -P2*PBdv*D1/pC**2 + (sqrt(2)*Del_d1v1/2 + Del_v1d1)/pC
  dPFdWA  <- apply(regr$barWA, 2, function(x){t(as.numeric(x)*t(dPFdWA))})
  dPFda <- -P1*PBdv*D2/pC**2 - sqrt(2)*Del_d1v1/(-2*P1P2 + 2)
  dPFda  <- apply(regr$bara, 2, function(x){t(as.numeric(x)*t(dPFda))})
  
  return( rbind(cbind(dPRdSA, dPRdVA, dPRdCB, dPRdWA, dPRdWB, dPRda, dPRdVB),
                cbind(dPFdSA, dPFdVA, dPFdCB, dPFdWA, dPFdWB, dPFda, dPFdVB))
  )
}

