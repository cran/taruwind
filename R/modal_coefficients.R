# Independent R preparation for the same neutral modal operator.
# Matrices enter native code explicitly row-major; terrain rows south to north.
modal_coefficients <- function(h, dx, scales, z0, penalty, reference=100) {
  stopifnot(is.matrix(h), all(is.finite(h)), dx>0, all(scales>0), z0>0)
  ny <- nrow(h); nx <- ncol(h); k <- length(scales); alpha <- 1/scales
  gx <- (h[,c(2:nx,1),drop=FALSE]-h[,c(nx,1:(nx-1)),drop=FALSE])/(2*dx)
  gy <- (h[c(2:ny,1),,drop=FALSE]-h[c(ny,1:(ny-1)),,drop=FALSE])/(2*dx)
  g2 <- gx^2+gy^2
  C <- 1/outer(alpha,alpha,'+')
  Cbar <- C+penalty*mean(g2/(1+g2))/2
  Q <- (1+mean(g2))*(outer(alpha,alpha)*C+penalty*outer(alpha,alpha))
  invR <- backsolve(chol(Cbar),diag(k))
  eig <- eigen(t(invR)%*%Q%*%invR,symmetric=TRUE)
  P <- invR%*%eig$vectors
  # All currently supported prototype profiles have alpha*z0<1. Retain an
  # explicit guard rather than silently using a divergent large-x series.
  E1 <- function(x) {
    stopifnot(x>0,x<=1)
    term <- -x; sum <- term
    for(i in 2:200){term <- term*(-x)/i; delta <- term/i; sum <- sum+delta;if(abs(delta)<1e-17)break}
    -0.5772156649015328606-log(x)-sum
  }
  beta <- exp(alpha*z0)*vapply(alpha*z0,E1,0.)/log1p(reference/z0)
  list(alpha=alpha,C=as.vector(t(C)),P=as.vector(t(P)),lambda=eig$values,beta=beta)
}

modal_solve <- function(h,dx,scales,z0,penalty,u,v,height=10,reference=100,
                        row=as.integer(nrow(h)%/%2),col=as.integer(ncol(h)%/%2),
                        rtol=1e-6,maxiter=600,test=NULL) {
  cc <- modal_coefficients(h,dx,scales,z0,penalty,reference)
  do.call(modal_native,c(list(h=as.vector(t(h)),ny=nrow(h),nx=ncol(h),dx=dx),cc,
    list(penalty=penalty,u=u,v=v,height=height,z0=z0,reference=reference,
         row=row,col=col,rtol=rtol,maxiter=maxiter,test=test)))
}
