#' Convert Parameter List to thetaMSAR Object
#'
#' Converts a list of MSAR parameters into a properly formatted thetaMSAR object
#' with appropriate attributes, dimensions, and names. Ensures consistency of
#' parameter structure for NHMSAR package functions.
#'
#' @param x List of MSAR parameters (from init or M-step).
#' @param label Character. Model type: "HH", "HN", "NH", or "NN". Default: "HH".
#' @param regime_names Character vector of regime names (optional).
#' @param ncov.emis Integer. Number of emission covariates.
#' @param ncov.trans Integer. Number of transition covariates.
#'
#' @return A thetaMSAR object with properly structured components:
#'   \describe{
#'     \item{A}{List of lag-1 coefficient matrices (if order > 0), one per regime}
#'     \item{A0}{Matrix of intercepts (M × d) with row/column names}
#'     \item{sigma}{List of covariance matrices (d × d), one per regime}
#'     \item{prior}{Matrix of initial probabilities (M × 1)}
#'     \item{transmat}{Transition matrix (M × M) with row/column names}
#'     \item{par.trans}{Transition covariates (if label includes "N" in position 1)}
#'     \item{par.emis}{Emission covariates (if label includes "N" in position 2)}
#'   }
#'   Plus attributes: class, label, NbRegimes, NbComp, order, n_par
#'
#' @details
#' **Function Purpose**:
#'
#' This function ensures that parameter lists conform to the thetaMSAR class
#' structure expected by NHMSAR functions. It:
#'
#' 1. **Validates Input**: Checks that x is a valid parameter list
#'
#' 2. **Adds Dimensions**: Reshapes parameters to correct matrix/array dimensions
#'
#' 3. **Adds Names**: Assigns row/column names like "Regime1", "Regime2", etc.
#'
#' 4. **Sets Attributes**: Adds class, label, counts (M, d, order)
#'
#' 5. **Computes n_par**: Counts total number of parameters (for BIC calculation)
#'
#' **Parameter Counting**:
#' \itemize{
#'   \item A matrices: M × d × d × order
#'   \item A0 intercepts: M × d
#'   \item sigma matrices: M × d × (d+1) / 2 (symmetric)
#'   \item prior: M - 1 (sum to 1 constraint)
#'   \item transmat: M × (M - 1) (row sum to 1 constraint)
#' }
#'
#' @note
#' \itemize{
#'   \item Called internally after each M-step to maintain structure
#'   \item Required before calling E-step functions
#'   \item Stops with error if input is not a valid parameter list
#'   \item The "_revised_2" suffix indicates version 2 (current)
#' }
#'
#' @seealso
#' \code{\link{init.theta.MSAR_revised_2}} which produces initial thetaMSAR objects
#' \code{\link{fit.MSAR_revised_2}} which calls this after each M-step
#'
#' @examples
#' \dontrun{
#' # After M-step, convert parameters
#' theta_list <- list(
#'   A = list(matrix(0.5, 2, 2)),
#'   A0 = c(1, 2),
#'   sigma = list(diag(2)),
#'   prior = c(0.5, 0.5),
#'   transmat = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2)
#' )
#'
#' # This would typically have attributes from previous theta
#' # Here we add them manually for illustration
#' attributes(theta_list)$NbRegimes <- 2
#' attributes(theta_list)$NbComp <- 2
#' attributes(theta_list)$order <- 1
#' attributes(theta_list)$label <- "HH"
#'
#' theta_obj <- as.thetaMSAR_revised_2(theta_list, label = "HH")
#' }
#'
#' @export
as.thetaMSAR_revised_2 <-
function(x,label='HH',regime_names=NULL,ncov.emis=ncov.emis,ncov.trans=ncov.trans) {
  if (!NHMSAR:::is.thetaMSAR(x)) {
    stop('as.thetaMSAR: your input is not like a theta at all')
  }
  
  dimname <- dimnames(x)
  att <- attributes(x)
  att$dimnames <- NULL
  att$names <- NULL
  att$class <- NULL
  
  label = att$label
  M = att$NbRegimes #number of regimes
  order=att$order
  d=att$NbComp
  
  dimname=c("A0","sigma","prior","transmat")
  
  n_par=NULL
  
  if(order>0){dimname=c("A",dimname)}
  if(substr(label,1,1)=="N"){dimname=c(dimname,'par.trans'); n_par=n_par+length(c(x$par.trans))/M}
  if(substr(label,2,2)=="N"){dimname=c(dimname,'par.emis'); n_par=n_par+length(c(x$par.emis))/M}
  
  names(x)=dimname
  
  x$A0 = matrix(x$A0,M,d)
  rownames(x$A0)=c(paste("Regime",1:M,sep=""))
  colnames(x$A0)=c(paste("A0",1:d,sep=""))
  
  # x$prior = matrix(x$prior,M,d) # patched: Not needed
  x$prior=matrix(x$prior,M,1)
  rownames(x$prior)=c(paste("Regime",1:M,sep=""))
  colnames(x$prior)=""
  
  rownames(x$transmat)=c(paste("Regime",1:M,sep=""))
  colnames(x$transmat)=c(paste("Regime",1:M,sep=""))
  
  if(substr(label,2,2)=="N"){
    names(x$par.emis)=c(paste('Regime',1:M,sep=""))
    for(i in 1:M){
      if (!is.null(dim(x$par.emis[[i]]))) {rownames(x$par.emis[[i]])=rep("",d)}
      #if (!is.null(dim(x$par.emis[[i]]))) {rownames(x$par.emis[[i]])=c(paste("Regime",1:M,sep=""))}
      #   		colnames(x$par.emis[[i]])=c(paste('coef.emis',1:max(2,ncov.emis+1),sep=""))
      #   		if (!is.null(dim(x$par.emis[[i]])[2])) {colnames(x$par.emis[[i]])=c(paste('coef.emis',1:(length(c(x$par.emis[[i]]))),sep=""))}
      if (!is.null(dim(x$par.emis[[i]])[2])) {
        colnames(x$par.emis[[i]])=c(paste('coef.emis',1:dim(x$par.emis[[i]])[2],sep=""))
      }   	
    }
  }
  
  if(substr(label,1,1)=="N"){
    rownames(x$par.trans)=c(paste("Regime",1:M,sep=""))
    colnames(x$par.trans)=c(paste('coef.trans',1:max(2,ncov.trans+1),sep=""))	
  }
  
  if(d==1){
    x$sigma = matrix(x$sigma,M,d)
    rownames(x$sigma)=c(paste("Regime",1:M,sep=""))
    colnames(x$sigma)=""
  }else{
    for(i in 1:M){
      colnames(x$sigma[[i]])=rep("",d)
      rownames(x$sigma[[i]])=rep("",d)
    }
  }
  
  if(order>0 && d==1){
    rownames(x$A)=c(paste('Regime',1:M,sep=""))
    colnames(x$A)=c(paste('A',1:order,sep=""))
  }
  else if (order>0){
    names(x$A)=c(paste("Regime",1:M,sep=""))
    for(i in 1:M){
      names(x$A[[i]])=c(paste('A',1:order,sep=""))
      for(j in 1:order){
        colnames(x$A[[i]][[j]])=rep("",d)
        rownames(x$A[[i]][[j]])=rep("",d)
      }
    }
    names(x$sigma)=c(paste('Regime',1:M,sep=""))
  }
  
  if (is.null(att$order)) { att$order <- order }
  if (is.null(att$n_par)) { att$n_par <- order*d^2+d+d^2+1+M }
  
  #     if (is.null(att$label)) { att$label <- label } # A VOIR !!!
  
  #    x <- array(x,thetadim,dimnames=dimname)
  for (a in names(att)) {
    attr(x,a) <- att[[a]]
  }
  
  class(x) <- "MSAR"
  return(x)
}
