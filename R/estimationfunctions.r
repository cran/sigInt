vec2U.norm <- function(x,X=matrix(0)){
  param <- list(
    barWA = x[1],
    barWB = x[2] + X*x[3],
    bara = x[4],
    VA = x[5],
    VB = x[6],
    SA = 0,
    CB = 0,
    sig=1
  )	
  return(param)
}


vec2U.regr <- function(x,regr, fixed.par){
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
  
  
  param <- list(
    barWA = regr[[4]] %*% x[indx[[4]]],
    barWB = regr[[5]] %*% x[indx[[5]]],
    bara = regr[[6]]  %*% x[indx[[6]]],
    VA = regr[[2]] %*% x[indx[[2]]], 
    VB =  regr[[7]] %*% x[indx[[7]]], 
    SA = regr[[1]] %*% x[indx[[1]]],
    CB = regr[[3]] %*% x[indx[[3]]],
    sig=1
  )	
  if(length(fixed.par)>0){
    fixed.par <- lapply(fixed.par, rep, max(sapply(regr, nrow)))
    param <- modifyList(param, fixed.par)
  }
  param <- lapply(param, as.numeric)
  return(param)
}




f.jo <- function(p, U){
  return(pnorm((p*U$barWB + (1-p)*U$VB - U$CB)/(U$sig*p)))
}

cStar.jo <- function(p, U){
  return((U$SA - (1-p)*U$VA)/p)
}

h.jo <- function(c, U){
  d1 <- (U$barWA - U$bara)/(U$sig*sqrt(2))
  d2 <- (U$barWA - c)/(U$sig)
  
  return(
    pbivnorm(d1, d2,rho=1/sqrt(2)) 
  )
}

g.jo <- function(c,U){
  v1 <- (c-U$barWA)/U$sig
  v2 <- (c-U$bara)/U$sig
  return(1 - pnorm(v1)*pnorm(v2))
}

# this is the equilibrium constraint
const.jo <- function(p, U){
  c <- cStar.jo(p,U)
  g <- g.jo(c,U)
  g[g<=.Machine$double.eps]<-.Machine$double.eps
  j <- h.jo(c,U)/g
  return(p - f.jo(j,U)) 
} 


eqProbs <- function(p, U,RemoveZeros=F){
  ck <- cStar.jo(p,U)
  pC <- g.jo(ck, U)
  if (RemoveZeros){
    pC[pC <= .Machine$double.eps] <- .Machine$double.eps
  }
  pF <- h.jo(ck, U)/pC
  return(cbind(p, pC, pF))
  
}


QLL.jo  <- function(x,PRhat,PFhat,Y,regr, fixed.par){
  # here x = (theta)	
  U <- vec2U.regr(x,regr, fixed.par)
  PR <- f.jo(PFhat, U)
  PR[PR<=.Machine$double.eps] <- .Machine$double.eps 
  PC <- g.jo(cStar.jo(PRhat,U),U)
  PC[PC<=.Machine$double.eps] <- .Machine$double.eps
  PF <- h.jo(cStar.jo(PRhat,U),U)/PC
  
  
  OUT <- cbind(1-PC, PC*(1-PR), 	
               PC*PR*PF,PC*PR*(1-PF))
  OUT[OUT <= sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  QLL <- sum(log(t(OUT))*Y)
  return(-QLL)
}


QLL.jo_i  <- function(x,PRhat,PFhat,Y,regr, fixed.par){
  # here x = (theta)	
  U <- vec2U.regr(x,regr, fixed.par)
  PR <- f.jo(PFhat, U)
  PR[PR<=.Machine$double.eps] <- .Machine$double.eps 
  PC <- g.jo(cStar.jo(PRhat,U),U)
  PC[PC<=.Machine$double.eps] <- .Machine$double.eps
  PF <- h.jo(cStar.jo(PRhat,U),U)/PC
  
  
  OUT <- cbind(1-PC, PC*(1-PR), 	
               PC*PR*PF,PC*PR*(1-PF))
  OUT[OUT <= sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  QLL <- colSums(log(t(OUT))*Y)
  return(-QLL)
}



QLL.jo_ip  <- function(x,Phat,Y,regr, fixed.par){
  PRhat <- Phat[1:ncol(Y)]
  PFhat <- Phat[(ncol(Y)+1):length(Phat)]
  # here x = (theta)	
  U <- vec2U.regr(x,regr,  fixed.par)
  PR <- f.jo(PFhat, U)
  PR[PR<=.Machine$double.eps] <- .Machine$double.eps 
  PC <- g.jo(cStar.jo(PRhat,U),U)
  PC[PC<=.Machine$double.eps] <- .Machine$double.eps
  PF <- h.jo(cStar.jo(PRhat,U),U)/PC
  
  
  OUT <- cbind(1-PC, PC*(1-PR), 	
               PC*PR*PF,PC*PR*(1-PF))
  OUT[OUT <= sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  QLL <- colSums(log(t(OUT))*Y)
  return(-QLL)
}





signalEQC <- function(Uold, Unew, length.out=100, gridsize=1e6, comp=F, tol=1e-8){
	# This function computes all equilibria as you vary payoff from Uold to Unew
		# Uold, Unew: lists of model parameters
		# length.out: how many steps should we take?
		# gridsize: size of the grid search
		# comp: should we compute eq after grid search?
		# tol: if yes, what should the tolerance be?
	
	out <- list()
	grid <- seq(from=0, to=1, length.out=gridsize) #Actually include 0 and 1?
	
	 
  for (i in 0:(length.out-1)){
    
    Ui <- as.list(mapply("+", lapply(Uold, "*", 1-i/(length.out-1)), lapply(Unew, "*", i/(length.out-1))))
    fgrid <- const.jo(grid,Ui)
    sols <- which(tail(fgrid,-1)*head(fgrid,-1) <= 0)
    sols <- matrix(grid[c(sols, sols+1)], nrow=length(sols))
    
    # compute equilibria
    if (!comp){
		sols <- rowMeans(sols)
    } else{
    	solver <- function(x){uniroot(function(x){const.jo(x,Ui)}, x, tol=tol)$root}     	
    	sols <- apply(sols,1,solver)    
    }
	
	#return parameters of interest
	if (!length(sols)){
      out[[i+1]]  <- cbind(i,NaN,NaN,NaN)
    } else {
        
       sols.pc <- g.jo(cStar.jo(sols,Ui),Ui)
       sols.pf <- h.jo(cStar.jo(sols,Ui),Ui)/sols.pc
       index <- is.nan(sols.pc*sols*sols.pf)
       sols.onset <- c(sols.pc[index], (sols.pc*sols*sols.pf)[!index])
      
          out[[i+1]] <- cbind(i,
                        sols,
                        sols.pc,
                        sols.onset,
                        sols.pf)
        }

	}
	return(do.call(rbind,out))
}


