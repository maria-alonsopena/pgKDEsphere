#' pi.kappa
#'
#' Function \code{pi.kappa} computes a plug-in type smoothing parameter for the parametrically guided (hyper)spherical kernel density estimator, equipped with a von Mises-Fisher guide.
#'
#' @param datax Matrix containing the data in cartesian coordinates, where the number of rows is the number of observations and the number of columns is the dimension of the Euclidean space where the sphere is embebed.
#' @param mu0 Vector containing the mean direction of the von Mises-Fisher guide.
#' @param tau0 Numerical value containing the concentration of the von Mises-Fisher guide.
#' @param guide Logical; if TRUE, the estimator with a von Mises-Fisher as guide is computed. If FALSE, the classical kernel density estimator without guide is computed (equivalent to uniform guide).
#' @return An object with class "sphkde" whose underlying structure is a list containing the following components:
#' @examples
#' library(Directional)
#' library(movMF)
#' # Data generation
#' n<-200
#' mu<-matrix(c(0,0,1,0,0,-1),ncol=3,byrow=TRUE)
#' k<-c(7,2)
#' probs<-c(0.85,0.15)
#' datax<-rmovMF(n,k*mu,alpha=probs)
#' # Estimation of parameters of a vMF
#' param<-vmf.mle(datax)
#' mu0<-param$mu
#' tau0<-param$kappa
#' # Selection of the smoothing parameter
#' kappa <- pi.kappa(datax,mu0,tau0)
#' @details
#' See Alonso-Pena et al. (2023) for details.
#' @references
#' Alonso-Pena, M., Claeskens, G. and Gijbels, I. (2023) Nonparametric estimation of densities on the hypersphere using a parametric guide. Under review.
#' @export



pi.kappa<-function(datax,mu0,tau0,guide=TRUE){


  if(is.matrix(datax)==FALSE){stop("Object datax must be a matrix")}
  datax <- na.omit(datax)

  if (guide!= TRUE & guide!= FALSE){stop("guide must be either TRUE or FALSE")}

  if (guide == FALSE){g<-bw_dir_ami(datax);kappa_opt<-1/g^2}else{
  	fit_mix <- DirStats::bic_vmf_mix(data = datax)
  	mu.mix <- fit_mix$best_fit$mu_hat
  	tau.mix <- fit_mix$best_fit$kappa_hat
 	pi.mix <- fit_mix$best_fit$p_hat
  	d <- ncol(datax) - 1
  	if(d==1){R<-int_cir(function(x){Rgr(x,mu.mix,tau.mix,pi.mix,mu0,tau0)})}
  	if(d==2){R<-int_sph(function(x){Rgr(x,mu.mix,tau.mix,pi.mix,mu0,tau0)})}
  	if(d>2){R<-int_hypsph(function(x){Rgr(x,mu.mix,tau.mix,pi.mix,mu0,tau0)},q=d)}

  	kappa_opt<-hopt(datax, R, fit_mix = fit_mix)
  }
  return(kappa_opt)
}
