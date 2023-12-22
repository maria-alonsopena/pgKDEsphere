#' sphkde.pg
#'
#' Function \code{sphkde.pg} computes the kernel density estimator for (hyper)spherical data with a parametric guide, which corresponds to the von Mises-Fisher model.
#'
#' @param datax Matrix containing the data in cartesian coordinates, where the number of rows is the number of observations and the number of columns is the dimension of the Euclidean space where the sphere is embebed.
#' @param kappa Smoothing parameter. It refers to the concentration when employing a von Mises-Fisher kernel.
#' @param eval.points Matrix containing the evaluation points for the estimation of the density.
#' @param guide Logical; if TRUE, the estimator with a von Mises-Fisher as guide is computed. If FALSE, the classical kernel density estimator without guide is computed (equivalent to uniform guide).
#' @return An object with class "sphkde" whose underlying structure is a list containing the following components: \item{estim}{ The estimated values of the density.}
#' \item{kappa}{The smoothing parameter used.The n coordinates of the points where the regression is estimated.}
#' \item{eval.points}{The points where the estimated density was evaluated.}
#' \item{data}{Original dataset.}
#' @examples
#' library(movMF)
#' n<-200
#' mu<-matrix(c(0,0,1,0,0,-1),ncol=3,byrow=TRUE)
#' k<-c(7,2)
#' probs<-c(0.85,0.15)
#' datax<-rmovMF(n,k*mu,alpha=probs)
#' est<-sphkde.pg(datax,guide=TRUE)
#' sphkde.plot(est,type="sph")
#' @details
#' See Alonso-Pena et al. (2023) for details.
#' @references
#' Alonso-Pena, M., Claeskens, G. and Gijbels, I. (2023) Nonparametric estimation of densities on the hypersphere using a parametric guide. Under review.
#' @export



sphkde.pg<-function(datax, kappa = NULL, eval.points = NULL, guide = TRUE){


  if(is.matrix(datax)==FALSE){stop("Object datax must be a matrix")}
  datax <- na.omit(datax)

  if (guide!= TRUE & guide!= FALSE){stop("guide must be either TRUE or FALSE")}

  d <- ncol(datax) - 1

  if (is.null(eval.points)) {
    if(d==1){
      grid<-seq(0,2*pi,length=400)
	eval.points<-to_cir(grid)
    }
    if(d==2){
	N<-300
	ele_grid<-seq(0,pi,length=N)
	ori<-list()

	ori[[1]]<-0
	ori[[N]]<-0
	for (i in 2:(N-1)){
	  ori[[i]]<- seq(0,2*pi,length=floor(2*N*sin(ele_grid[i])))
	}
	lonxi<-as.numeric(lapply(ori,length))
	both<-cbind(unlist(ori),rep(ele_grid,times=lonxi))
	eval.points<-to_sph(both[,1], both[,2]) # conversion of coordinates
    }
    if(d>2){
	grid<-list()
	for (j in 1:(d-1)){grid[[j]]<-seq(0,pi,length=7)}
	grid[[d]]<-seq(0,2*pi,length=7)
	allgrid<-expand.grid(grid)
	eval.points<-cart_sph(allgrid)
    }
  }else{
     if(is.matrix(eval.points)==FALSE){stop("Object eval.points must be a matrix")}
  }


  param<-vmf.mle(datax)
  mu0<-param$mu
  tau0<-param$kappa


  if (is.null(kappa)) {
    kappa <- pi.kappa(datax,mu0,tau0)
  }else{
	if(inherits(kappa,"numeric")){stop("kappa must be numeric")}else{if(length(kappa)<=0)stop("kappa must have positive length")}
  }


  if(guide==FALSE){tau0<-0}


  cl<-(kappa^(-1/2))^(d-1)*(2*pi)^((d + 1)/2)*besselI(1/(kappa^(-1/2))^2,(d - 1)/2,expon.scaled=TRUE)
  if(d==1){R<-int_cir(function(x){pgKDE_sph_vMF(x,datax,kappa^(-1/2),cl,mu0,tau0)})}
  if(d==2){norm<-int_sph(function(x){pgKDE_sph_vMF(x,datax,kappa^(-1/2),cl,mu0,tau0)})}
  if(d>2){norm<-int_hypsph(function(x){pgKDE_sph_vMF(x,datax,kappa^(-1/2),cl,mu0,tau0)},q=d)}

  est<-pgKDE_sph_vMF(eval.points,datax,kappa^(-1/2),cl,mu0,tau0)/norm

  ret<-list(estim = est, kappa = kappa,  eval.points = eval.points, data = datax)
  class(ret)<-"sphkde"
  return(ret)
}
