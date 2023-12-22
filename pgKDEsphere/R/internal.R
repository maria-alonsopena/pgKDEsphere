
plot_sphkde<-function(grid,fgrid,axis){

  hdr.colors1<- colorRampPalette(c("blue","green" ),alpha=TRUE)(10)
  hdr.colors2<- colorRampPalette(c("green","yellow"),alpha=TRUE)(10)
  hdr.colors3<- colorRampPalette(c("yellow","orange"),alpha=TRUE)(10)
  hdr.colors4<- colorRampPalette(c("orange","red"),alpha=TRUE)(10)
  hdr.colors5<- colorRampPalette(c("red","darkred"),alpha=TRUE)(10)
  mycol<-c(hdr.colors1,hdr.colors2,hdr.colors3,hdr.colors4,hdr.colors5)

  interv<-seq(0,1,length=51)
  vector_int<-findInterval(fgrid, interv)
  vector_int[vector_int==51]<-50

  colors<-mycol[vector_int]

  zoom<-0.8
  windowRect<-c(500,50,0,0)
  windowRect[3]=windowRect[1]+256*2
  windowRect[4]=windowRect[2]+256*2

  open3d(zoom = zoom,  windowRect=windowRect)

  points3d(grid[,1],grid[,2],grid[,3],col=colors,radius=1,size=2)#ollo co size

  if(axis==TRUE){
    R<-1;Xaxis1<-seq(-R-0.5,R+0.5,length=100);Xaxis2<-rep(0,100);Xaxis3<-rep(0,100)
    points3d(Xaxis1,Xaxis2,Xaxis3,col=1)
    text3d(R+0.65,0,0,texts="x",col=1,cex=.9,font=2)

    R<-1;Yaxis1<-rep(0,100);Yaxis2<-seq(-R-0.5,R+0.5,length=100);Yaxis3<-rep(0,100)
    points3d(Yaxis1,Yaxis2,Yaxis3,col=1)
    text3d(0,R+0.65,0,texts="y",col=1,cex=.9,font=2)

    R<-1;Zaxis1<-rep(0,100);Zaxis2<-rep(0,100);Zaxis3<-seq(-R-0.5,R+0.5,length=100)
    points3d(Zaxis1,Zaxis2,Zaxis3,col=1)
    text3d(0,0,R+0.65,texts="z",col=1,cex=.9,font=2)
  }
}


cart_sph<-function(grid){

	# grid has d columns, t rows
	d<-ncol(grid)
	cart<-matrix(nrow=nrow(grid),ncol=d+1)
	cart[,1]<-cos(grid[,1])
	for (j in 2:(d-1)){ cart[,j]<-cos(grid[,j])*rowProds(as.matrix(sin(grid[,1:(j-1)]))) }
	cart[,d]<-sin(grid[,d])*rowProds(as.matrix(cos(grid[,1:(d-1)])))
	cart[,d+1]<-cos(grid[,d])*rowProds(as.matrix(cos(grid[,1:(d-1)])))

	return(cart)
}


Rgr<-function(x,mu.mix,tau.mix,pi.mix,mu0,tau0){


  p<-ncol(mu.mix); d<-p-1


  ctauj<-c_vMF(p = p, kappa = tau.mix)
  ctau0<-c_vMF(p = p, kappa = tau0)

  xmuj<-x %*% t(mu.mix)
  etxmuj<-exp(t(tau.mix*t(xmuj)))
  dvMFj<-t(pi.mix*ctauj*t(etxmuj))

  xmu0<-x%*%mu0
  g0<-ctau0*exp(tau0*xmu0)

  sumita<-colSums(t(mu.mix)*mu0)

  term<-rowSums(dvMFj*(tau0*xmuj - t(tau.mix*t(xmuj)) + matrix(tau.mix^2/d,nrow=nrow(x),ncol=length(tau.mix),byrow=TRUE)
                       - matrix(2*tau.mix*tau0*sumita/d,nrow=nrow(x),ncol=length(tau.mix),byrow=TRUE)
                       + tau0/d - t(tau.mix^2*t(xmuj^2/d))
                       + 2*d^(-1)*tau0*t(tau.mix*t(xmuj*as.numeric(xmu0)))  ))

  return((term)^2)
}


hopt<-function (datax, R, fit_mix){
  stopifnot(is.matrix(datax))
  data <- na.omit(datax)
  n <- nrow(datax)
  d <- ncol(datax) - 1
  lambda <- DirStats::lambda_L(q = d)
  b <- DirStats::b_L(q = d)
  dd <- DirStats::d_L(q = d)
  mu <- fit_mix$best_fit$mu_hat
  kappa <- fit_mix$best_fit$kappa_hat
  p <- fit_mix$best_fit$p_hat
  h <- ((d* dd)/(4 * b^2 * lambda * R * n))^(1/(d + 4))

  return(1/h^2)
}


