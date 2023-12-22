#' sphkde.plot
#'
#' Function \code{sphkde.plot} provides a graphical representation of the parametrically guided kernel density estimator for spherical and circular data. For circular data, both linear and circular representations are available. For spherical data, an interactive 3D spherical representation is provided.
#'
#' @param object Object of the class \code{sphkde}.
#' @param type Character string giving the desired type of plot. For circular data, it can be "sph" for a circular representation or "line" for a linear representation. For spherical data the value "sph" is required.
#' @param axis Logical; if TRUE, the axis are represented in the spherical representation. If FALSE, axis are not represented. Only for spherical representations.
#' @param shrink Numeric parameter that controls the size of the plotted circle in the circular representations. Default is 1.3. Larger values shrink the circle, while smaller values enlarge the circle.
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



sphkde.plot<-function(object, type = "sph", axis = TRUE, shrink = 1.2){

	if(inherits(object,"list")){stop("object must be of class sphkde")}

	d <- ncol(object$data) - 1


	if(type!="sph" & type!= "line"){stop("type must be either sph or line")}

	if (axis!= TRUE & axis!= FALSE){stop("axis must be either TRUE or FALSE")}


	if(d==1){
		if(type=="sph"){ptype<-"circ"}else if(type=="line"){ptype<-"line"}
	}

	if(d==2){
		if(type=="linear"){stop("Type linear is only available for dimension d=1")}
		if(type=="sph"){ptype<-"sph"}
	}

	if(d>2){stop("Graphical representation is only available for dimensions d=1 and d=2")}

	ev<-object$eval.points
	est<-object$estim

	if (ptype=="sph"){
		plot_sphkde(ev,est,axis)
	}

	if (ptype=="circ"){
		dat<-to_rad(object$data)
		gridc<-to_rad(ev)
		rose.diag(dat, bins=16, col="darkgrey", cex=1.5, prop=1.5,shrink = shrink)
		lines.circular(gridc,est,lwd=2,join=TRUE)
	}

	if (ptype=="line"){
		gridc<-to_rad(ev)
		plot(gridc,est,type="l",lwd=2,ylim=c(0,0.3),ylab="Density",xaxt="n" ,xlab="")
		axis(1,c(0,pi/2,pi,3*pi/2,2*pi),c(expression(0),expression(pi/2),expression(pi),expression(3*pi/2),expression(2*pi)))
	}

}
