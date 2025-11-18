#' Fit MSAR Model via Expectation-Maximization Algorithm
#'
#' Estimates parameters of a Markov-Switching Autoregressive (MSAR) model using
#' the EM algorithm with optional LASSO regularization. Supports various model
#' configurations including homogeneous/non-homogeneous transitions and emissions.
#'
#' @param data 3D array of time series data (time × samples × variables).
#' @param theta Initial parameter object from \code{init.theta.MSAR_revised_2}.
#' @param MaxIter Integer. Maximum EM iterations. Default: 100.
#' @param eps Numeric. Convergence threshold for log-likelihood change. Default: 1e-5.
#' @param verbose Logical. Print iteration progress. Default: TRUE.
#' @param covar.emis Array of emission covariates (optional).
#' @param covar.trans Array of transition covariates (optional).
#' @param method Character. Estimation method for non-homogeneous models.
#' @param constraints Logical. Apply constraints on parameters. Default: FALSE.
#' @param reduct Logical. Use reduced M-step (sets small coefficients to zero). Default: FALSE.
#' @param K Matrix. Constraint matrix if constraints=TRUE.
#' @param d.y Vector. Constraint vector if constraints=TRUE.
#' @param ARfix Logical. Fix AR parameters. Default: FALSE.
#' @param penalty Character or Logical. Penalty type: FALSE, "LASSO", "ridge", "SCAD". Default: FALSE.
#' @param sigma.diag Logical. Force diagonal covariance matrices. Default: FALSE.
#' @param sigma.equal Logical. Force equal covariance across regimes. Default: FALSE.
#' @param lambda1 Numeric. LASSO/SCAD penalty parameter 1. Default: 0.1.
#' @param lambda2 Numeric. Ridge/SCAD penalty parameter 2. Default: 0.1.
#' @param a Numeric. SCAD tuning parameter. Default: 3.7.
#' @param ... Additional arguments passed to M-step functions.
#'
#' @return List containing:
#'   \describe{
#'     \item{theta}{Fitted parameter object with components:
#'       \itemize{
#'         \item A: List of lag-1 coefficient matrices for each regime
#'         \item A0: List of intercept vectors
#'         \item sigma: List of covariance matrices
#'         \item prior: Initial regime probabilities
#'         \item transmat: Regime transition matrix
#'       }
#'     }
#'     \item{loglik}{Final log-likelihood value}
#'     \item{BIC}{Bayesian Information Criterion (if computed)}
#'     \item{Npar}{Number of parameters (if computed)}
#'     \item{ll_history}{Vector of log-likelihood values across iterations}
#'     \item{FB}{Forward-backward output from final E-step}
#'     \item{converged}{Convergence status and values}
#'   }
#'
#' @details
#' **EM Algorithm**:
#'
#' The function alternates between:
#'
#' 1. **E-step** (\code{Estep.MSAR}): Computes regime probabilities given current parameters
#'    using forward-backward algorithm
#'
#' 2. **M-step**: Updates parameters to maximize expected log-likelihood. Choice depends on:
#'    \itemize{
#'      \item **label='HH'** (Homogeneous): \code{Mstep.hh.MSAR} or specialized versions
#'      \item **penalty="LASSO"**: First iteration uses \code{Mstep.hh.lasso.MSAR_patched_2},
#'        subsequent iterations use \code{Mstep.hh.reduct.MSAR_patched_2}
#'      \item **penalty="ridge"**: Ridge regression penalty
#'      \item **penalty="SCAD"**: Smoothly Clipped Absolute Deviation penalty
#'    }
#'
#' **Convergence**: Stops when |loglik[t] - loglik[t-1]| < eps or MaxIter reached.
#'
#' **LASSO Strategy**: Iteration 1 estimates sparse networks via LASSO. Iterations 2+
#' refine only non-zero edges (reduct method) for efficiency.
#'
#' **Label Types**:
#' \itemize{
#'   \item HH: Homogeneous transition & emission
#'   \item HN: Homogeneous transition, Non-homogeneous emission
#'   \item NH: Non-homogeneous transition, Homogeneous emission
#'   \item NN: Non-homogeneous transition & emission
#' }
#'
#' @note
#' \itemize{
#'   \item Requires initial parameters from \code{init.theta.MSAR_revised_2}
#'   \item LASSO is the standard penalty for network recovery applications
#'   \item Convergence is not guaranteed; may stop at local optimum
#'   \item Log-likelihood should be non-decreasing (EM guarantee)
#' }
#'
#' @seealso
#' \code{\link{init.theta.MSAR_revised_2}} for initialization
#' \code{\link{init_and_fit.MSAR_Lasso_2}} for combined init+fit with retry
#' \code{\link{Mstep.hh.lasso.MSAR_patched_2}} for LASSO M-step
#' \code{\link{Mstep.hh.reduct.MSAR_patched_2}} for reduced M-step
#' \code{\link{EM_converged_patched_2}} for convergence checking
#'
#' @examples
#' \dontrun{
#' # Initialize parameters
#' data_array <- array(rnorm(1000 * 4), dim = c(1000, 1, 4))
#' theta_init <- init.theta.MSAR_revised_2(data_array, M = 2, order = 1)
#'
#' # Fit with LASSO
#' fit <- fit.MSAR_revised_2(
#'   data = data_array,
#'   theta = theta_init,
#'   penalty = "LASSO",
#'   MaxIter = 200,
#'   eps = 1e-5,
#'   verbose = TRUE
#' )
#'
#' # Check convergence
#' print(fit$converged)
#' plot(fit$ll_history)  # Should show increasing log-likelihood
#' }
#'
#' @export
source("R/estimation/as.thetaMSAR_revised_2.R")
source("R/estimation/Mstep.hh.lasso.MSAR_patched_2.R")
source("R/estimation/Mstep.hh.reduct.MSAR_patched_2.R")
source("R/estimation/EM_converged_patched_2.R")

fit.MSAR_revised_2 <-
function(
    data,theta,MaxIter=100,eps=1e-5,verbose=TRUE,
    covar.emis=NULL,covar.trans=NULL,method=NULL,constraints=FALSE,reduct=FALSE,K=NULL,d.y=NULL,ARfix=FALSE,penalty=FALSE,sigma.diag=FALSE,sigma.equal=FALSE,lambda1=.1,lambda2=.1,a=3.7,...
) { 
  #browser()
  cl <- match.call()
  now <- Sys.time()
  if (missing(theta)) {
    stop('can not fit a MSAR model without initial value for theta')
  }
  
  att <- attributes(data)   
  att.theta <- attributes(theta)
  
  data <- as.array(data)
  
  T <- dim(data)[1]
  if(is.null(T) || is.na(T)){T <- length(data)}
  N.samples <- dim(data)[2]
  if(is.null(N.samples) || is.na(N.samples)){N.samples <- 1}
  d <- att.theta$NbComp
  if(is.null(d)|is.na(d)){d <- 1}
  
  data <- array(data,c(T,N.samples,d))
  
  order <- att.theta$order
  label <- att.theta$label
  M <- att.theta$NbRegimes
  
  if ((missing(covar.trans) && substr(label,1,1)=="N") || (missing(covar.emis) && substr(label,2,2)=="N")) {
    stop("Can not fit a non homogeneous MSAR without covariable")
  }
  if (length(covar.trans)==1) {
    Lag = covar.trans+1
    covar.trans = array(data[(1):(T-Lag+1),,],c(T-Lag+1,N.samples,d))
    data =  array(data[Lag:T,,],c(T-Lag+1,N.samples,d))
  }
  
  if (!is.null(covar.emis)) {
    if (is.null(dim(covar.emis)) ) {ncov.emis=1}
    else if (is.na(dim(covar.emis)[3])) {ncov.emis=1}
    else if (!is.na(dim(covar.emis)[3])) {ncov.emis = dim(covar.emis)[3]}
    covar.emis = array(covar.emis,c(T,N.samples,ncov.emis))
  }
  # si dim(covar.emis)[3] = NA, mattre sous forme d'array
  else {ncov.emis = 0}
  if (!is.null(covar.trans)) {ncov.trans = dim(covar.trans)[3]}
  else {ncov.trans = 0}
  
  BIC = NULL
  Npar = NULL
  
  cnt <- 0
  FB <- NHMSAR:::Estep.MSAR(data,theta,covar.trans=covar.trans,covar.emis=covar.emis)
  loglik = FB$loglik
  previous_loglik <- FB$loglik-1000 
  ll_history = NULL
  converged = EM_converged_patched_2(0,2*eps,eps);
  par = NULL
  while (converged[1]==0 && cnt < MaxIter) {
    cnt <- cnt+1
    
    if (verbose) {print(c("iteration ",cnt,"  loglik = ",loglik),quote = FALSE)}
    
    # -----------------------------
    # ...... M step
    if (label=='HH') {
      if (constraints == FALSE & penalty==FALSE & reduct==FALSE ) {
        par = NHMSAR:::Mstep.hh.MSAR(data,theta,FB,sigma.diag=sigma.diag,sigma.equal=sigma.equal) 
      } else if (constraints){
        par = NHMSAR:::Mstep.hh.MSAR.with.constraints(data,theta,FB,K=K,d.y=d.y) 
        attributes(theta)$n_par = M + M*(M-1) + 2*M*d # A and sigma diagonal
      } 
      else if (penalty=="ridge")	{
        par = NHMSAR:::Mstep.hh.SCAD.cw.MSAR(data,theta,FB,penalty="ridge",lambda1=0,lambda2=lambda2,par=par)
        #par = NHMSAR:::Mstep.hh.ridge.MSAR(data,theta,FB,lambda=lambda2)
      }
      else if (penalty=="LASSO")	{
        
        # 
        # call PATCHED !!!
        #
        
        if (cnt>1) {
          par = Mstep.hh.reduct.MSAR_patched_2(data,theta,FB,sigma.diag=sigma.diag)
        } 
        else {par = Mstep.hh.lasso.MSAR_patched_2(data,theta,FB)}
      }
      else if (penalty=="SCAD") {
        par = NHMSAR:::Mstep.hh.SCAD.MSAR(data,theta,FB,penalty="SCAD",lambda1=lambda1,lambda2=lambda2,par=par)
      }    
      #   		else if (penalty=="SIS") {
      #   			par = NHMSAR:::Mstep.hh.SIS.MSAR#(data,theta,FB,penalty="SCAD",lambda1=lambda1,lambda2=lambda2,par=par)
      #   		}
      else if (reduct) {
        
        # 
        # call PATCHED !!!
        #
        
        par = Mstep.hh.reduct.MSAR_patched_2(data,theta,FB,sigma.diag=sigma.diag)
      }
      
      if (order>0) {
        theta=list(par$A,par$A0,par$sigma,par$prior,par$transmat)
      } else {
        theta=list(par$A0,par$sigma,par$prior,par$transmat)
      } 		}                 
    else if (label=='HN') { 
      par = NHMSAR:::Mstep.hn.MSAR(data,theta,FB,covar=covar.emis,verbose=verbose) 
      if (order>0) {
        theta=list(par$A,par$A0,par$sigma,par$prior,par$transmat,par$par_emis)
      }       
      else{
        theta=list(par$A0,par$sigma,par$prior,par$transmat,par$par_emis)
      }
    }
    else if (label=='NH') { 
      par = NHMSAR:::Mstep.nh.MSAR(data,theta,FB,covar=covar.trans,method=method,ARfix=ARfix,reduct=reduct,sigma.diag=sigma.diag,
                          sigma.equal=sigma.equal,penalty=penalty,lambda1=lambda1,lambda2=lambda2,par=par) 
      if (order>0) {
        theta=list(par$A,par$A0,par$sigma,par$prior,par$transmat,par$par.trans)
      }          
      else{theta=list(par$A0,par$sigma,par$prior,par$transmat,par$par.trans)}
      
      
    }
    else if (label=='NN') { 
      par = NHMSAR:::Mstep.nn.MSAR(data,theta,FB,covar.emis=covar.emis,covar.trans=covar.trans,method=method) 
      if (order>0) {
        theta=list(par$A,par$A0,par$sigma,par$prior,par$transmat,par$par.trans,par$par.emis)
      }
      else{theta=list(par$A0,par$sigma,par$prior,par$transmat,par$par.trans,par$par.emis)}
      
    }
    ll_history[cnt] = loglik
    converged = EM_converged_patched_2(loglik, previous_loglik, eps)
    previous_loglik = loglik
    attributes(theta) = att.theta
    theta = as.thetaMSAR_revised_2(theta,label=label,ncov.emis = ncov.emis,ncov.trans=ncov.trans)             
    
    
    # -----------------------------
    # ...... E step
    FB = NHMSAR:::Estep.MSAR(data,theta,covar.emis=covar.emis,covar.trans=covar.trans)
    loglik = FB$loglik
    
    
  }
  #browser()
  if ( M>1) { 
    tr = NULL
    for (m in 1:M) {if (d==1) {tr[m] = theta$sigma[m]}
      else {tr[m] = sum(diag(theta$sigma[[m]]))}}
    i.tr = order(tr)
    theta$A0 = theta$A0[i.tr,]
    sigma.tmp = NULL
    A.tmp = theta$A
    for (m in 1:M) {
      sigma.tmp[[m]] = theta$sigma[[i.tr[m]]]
      if (order>0) {
        if (d>1) {
          for (o in 1:order ) {
            A.tmp[[m]][[o]] = theta$A[[i.tr[m]]][[o]]
          }
        }
        else {A.tmp[m,] = theta$A[i.tr[m],]}}
    }
    theta$A = A.tmp
    theta$sigma = sigma.tmp
    theta$prior = theta$prior[i.tr,]
    temp = theta$transmat[i.tr,i.tr]
    attributes(temp)$dimnames[[2]] <- attributes(theta$transmat)$
      dimnames[[2]]
    attributes(temp)$dimnames[[2]] <- attributes(theta$transmat)$dimnames[[2]]
    theta$transmat = temp
    if (substr(label,2,2)=="N") {tmp = theta$par.emis 
    for (j in 1:M) {theta$par.emis[[j]] = tmp[[i.tr[j]]]}
    }
    if (substr(label,1,1)=="N") {theta$par.trans = theta$par.trans[i.tr,]} 
    theta = as.thetaMSAR_revised_2(theta,label=label,ncov.emis = ncov.emis,ncov.trans=ncov.trans)
    FB$probS = FB$probS[,,i.tr]
  }
  
  if (penalty!="SCAD" ) {lambda1=rep(0,M)}
  npar = M*d+M*(M-1)
  if (substr(label,1,1)=="N") {npar = npar+M*length(theta$par.trans[1,])}
  if (substr(label,2,2)=="N") {npar = npar+M*length(theta$par.emis[[1]])} # modified 2019/11/21
  for (m in 1:M) {
    npar = npar+sum(abs(theta$A[[m]][[1]])>0)
    if (penalty!="SCAD" | max(abs(lambda1))==0) {npar = npar+sum(abs(theta$sigma[[m]][upper.tri(theta$sigma[[m]],diag=TRUE)])>0)}
    else { npar = npar+sum(abs(par$sigma.inv[[m]][upper.tri(par$sigma.inv[[m]],diag=TRUE)])>1e-5)} # 1e-5 : arbitrary level... 
  }
  if (sigma.diag){npar = (M-1)+M*(M-1)+M*d+M*d}
  if (sigma.equal){npar = (M-1)+M*(M-1)+M*d+d*(d+1)/2}
  if (sigma.diag & sigma.equal){npar = (M-1)+M*(M-1)+M*d+d}
  attributes(theta)$n_par = npar
  
  BIC = -2*ll_history[cnt]+ attributes(theta)$n_par*log(length(c(data)))
  ll.pen = NULL
  if (penalty=="SCAD") {
    a=3.7
    if (length(lambda1)==1) {lambda1 = matrix(lambda1,1,M)}
    if (length(lambda2)==1) {lambda2 = matrix(lambda2,1,M)}
    pen = 0
    for (m in 1:M){
      w = matrix(0,d,d)
      if (lambda1[m]>0) {wi = solve(theta$sigma[[m]])
      abs.S = abs(theta$sigma[[m]])
      w = lambda1[m]*abs.S
      wA = which(abs.S>lambda1[m] & abs.S<=a*lambda1[m])
      w[wA] = -(abs.S[wA]^2-2*a*lambda1[m]*abs.S[wA]+lambda1[m]^2)/2/(a-1)
      wA = which(abs.S>a*lambda1[m])
      w[wA] = (a+1)^2*lambda1[m]^2/2
      w=matrix(w,d,d)
      }
      pen = pen+sum((w-diag(diag(w)))*abs(theta$sigma[[m]]))
      omega = matrix(0,d,d)
      if (lambda2[m]>0) {
        abs.A = abs(theta$A[[m]][[1]])
        omega = lambda2[m]*abs.A
        wA = which(abs.A>lambda2[m] & abs.A<=a*lambda2[m])
        omega[wA] = -(abs.A[wA]^2-2*a*lambda2[m]*abs.A[wA]+lambda2[m]^2)/2/(a-1)
        wA = which(abs.A>a*lambda2[m])
        omega[wA] = (a+1)^2*lambda2[m]^2/2
        omega=matrix(omega,d,d)
      }
      pen = pen+sum(omega*abs(theta$A[[m]][[1]]))
    }
    ll.pen = (FB$loglik-((T-1)*N.samples)*pen)
  }
  res = list(theta=theta,ll_history=ll_history,Iter=cnt,Npar=Npar,BIC=BIC,smoothedprob=FB$probS,ll.pen = ll.pen)
  class(res) <- "MSAR"
  res$call = cl
  res
}
