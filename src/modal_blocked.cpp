// Experimental batch target evaluator; source snapshot: operational_gate_20260908/modal.cpp.
// Same solve and physics; evaluate multiple coordinates/heights from one solution.
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#define POCKETFFT_CACHE_SIZE 16
#include "pocketfft_hdronly.h"
#include "blocked_apply.h"
// [[Rcpp::plugins(cpp17)]]
namespace {
using Vec=std::vector<double>;
struct Solver {
 int k,ny,nx;std::size_t N,D;double dx,gamma;
 Vec C,alpha,P,lam,beta,gx,gy,hx,hy,normal;
 Vec cp,cap,fx,fy,total,weighted,q,qx,qy,modes,back,denom;
 std::vector<std::complex<double>> hat;
 Solver(const Vec& h,int ky,int iy,int ix,double spacing,double penalty,
        Vec c,Vec a,Vec p,Vec l,Vec b):k(ky),ny(iy),nx(ix),N(std::size_t(iy)*ix),D(ky*N),dx(spacing),gamma(penalty),
        C(c),alpha(a),P(p),lam(l),beta(b),gx(N),gy(N),hx(N),hy(N),normal(N),
        cp(D),cap(D),fx(D),fy(D),total(N),weighted(N),q(N),qx(N),qy(N),modes(D),back(D),
        denom(k*ny*(nx/2+1)),hat(k*ny*(nx/2+1)) {
   for(int y=0;y<ny;++y)for(int x=0;x<nx;++x){auto i=idx(y,x);
     gx[i]=(h[idx(y,(x+1)%nx)]-h[idx(y,(x+nx-1)%nx)])/(2*dx);
     gy[i]=(h[idx((y+1)%ny,x)]-h[idx((y+ny-1)%ny,x)])/(2*dx);
     normal[i]=std::sqrt(1+gx[i]*gx[i]+gy[i]*gy[i]);hx[i]=gx[i]/normal[i];hy[i]=gy[i]/normal[i];
   }
   const int nf=nx/2+1;const double pi=std::acos(-1.);
   for(int a=0;a<k;++a)for(int y=0;y<ny;++y)for(int x=0;x<nf;++x){
     double sy=std::sin(2*pi*y/ny)/dx,sx=std::sin(2*pi*x/nx)/dx;
     denom[(std::size_t(a)*ny+y)*nf+x]=lam[a]+sy*sy+sx*sx;
   }
 }
 std::size_t idx(int y,int x)const{return std::size_t(y)*nx+x;}
 void apply(const Vec& v,Vec& out){
   modal_apply(k,ny,nx,dx,gamma,C.data(),alpha.data(),gx.data(),gy.data(),hx.data(),hy.data(),normal.data(),v.data(),cp.data(),cap.data(),fx.data(),fy.data(),total.data(),weighted.data(),q.data(),qx.data(),qy.data(),out.data());
 }
 void precondition(const Vec& v,Vec& out){
   std::fill(modes.begin(),modes.end(),0.);
   for(std::size_t start=0;start<N;start+=1024){auto end=std::min(N,start+1024);
    for(int j=0;j<k;++j)for(int i=0;i<k;++i){double coefficient=P[j*k+i];auto off=std::size_t(i)*N;auto in=std::size_t(j)*N;
     for(std::size_t t=start;t<end;++t)modes[off+t]+=coefficient*v[in+t];}}
   pocketfft::shape_t shape{std::size_t(k),std::size_t(ny),std::size_t(nx)},axes{1,2};
   pocketfft::stride_t sr{std::ptrdiff_t(N*sizeof(double)),std::ptrdiff_t(nx*sizeof(double)),std::ptrdiff_t(sizeof(double))};
   pocketfft::stride_t sc{std::ptrdiff_t(ny*(nx/2+1)*sizeof(std::complex<double>)),std::ptrdiff_t((nx/2+1)*sizeof(std::complex<double>)),std::ptrdiff_t(sizeof(std::complex<double>))};
   pocketfft::r2c<double>(shape,sr,sc,axes,pocketfft::FORWARD,modes.data(),hat.data(),1.,1);
   for(std::size_t t=0;t<hat.size();++t)hat[t]/=denom[t];
   pocketfft::c2r<double>(shape,sc,sr,axes,pocketfft::BACKWARD,hat.data(),back.data(),1./N,1);
   std::fill(out.begin(),out.end(),0.);
   for(std::size_t start=0;start<N;start+=1024){auto end=std::min(N,start+1024);
    for(int j=0;j<k;++j)for(int i=0;i<k;++i){double coefficient=P[i*k+j];auto off=std::size_t(i)*N;auto in=std::size_t(j)*N;
     for(std::size_t t=start;t<end;++t)out[off+t]+=coefficient*back[in+t];}}
 }
 double dot(const Vec& a,const Vec& b){double s=0;for(std::size_t i=0;i<D;++i)s+=a[i]*b[i];return s;}
 Rcpp::List solve(double u,double v,double tolerance,int maxiter,double height,double z0,double reference,double row,double col,Rcpp::Nullable<Rcpp::NumericMatrix> targets=R_NilValue){
   Vec rhs(D),solution(D,0.),r(D),z(D),direction(D),ad(D);double b2=0;
   for(int a=0;a<k;++a)for(std::size_t t=0;t<N;++t){auto i=std::size_t(a)*N+t;rhs[i]=-beta[a]*(u*gx[t]+v*gy[t]);b2+=rhs[i]*rhs[i];}
   r=rhs;int iteration=0;double rr=b2;
   if(b2>0){precondition(r,z);direction=z;double rz=dot(r,z);
     for(;iteration<maxiter;){
       apply(direction,ad);double pap=dot(direction,ad);if(!(pap>0)&&rr>0)Rcpp::stop("Non-positive CG curvature");
       double step=rz/pap;rr=0;
       for(std::size_t i=0;i<D;++i){solution[i]+=step*direction[i];r[i]-=step*ad[i];rr+=r[i]*r[i];}
       ++iteration;if(iteration%10==0)Rcpp::checkUserInterrupt();if(rr<=tolerance*tolerance*b2)break;
       precondition(r,z);double next=dot(r,z),weight=next/rz;
       for(std::size_t i=0;i<D;++i)direction[i]=z[i]+weight*direction[i];
       rz=next;
     }
   }
   apply(solution,ad);double residual2=0;for(std::size_t i=0;i<D;++i)residual2+=(ad[i]-rhs[i])*(ad[i]-rhs[i]);
   double residual=b2>0?std::sqrt(residual2/b2):0;
   if(residual>tolerance*1.01)Rcpp::stop("Actual residual %.6g exceeds tolerance",residual);
   if(targets.isNotNull()){
     Rcpp::NumericMatrix pts(targets),values(pts.nrow(),3);
     if(pts.ncol()!=3)Rcpp::stop("Targets must contain row, col, height");
     for(int j=0;j<pts.nrow();++j){
       double rr=pts(j,0),cc=pts(j,1),hh=pts(j,2);
       if(!std::isfinite(rr)||!std::isfinite(cc)||!std::isfinite(hh)||rr<0||cc<0||rr>ny-1||cc>nx-1||hh<=0)Rcpp::stop("Invalid batch target");
       int ya=int(std::floor(rr)),xa=int(std::floor(cc)),yb=std::min(ya+1,ny-1),xb=std::min(xa+1,nx-1);
       double wy=rr-ya,wx=cc-xa;
       int ys[4]={ya,ya,yb,yb},xs[4]={xa,xb,xa,xb};
       double ws[4]={(1-wy)*(1-wx),(1-wy)*wx,wy*(1-wx),wy*wx};
       for(int corner=0;corner<4;++corner){
         int y=ys[corner],x=xs[corner];auto t=idx(y,x);
         double profile=std::log1p(hh/z0)/std::log1p(reference/z0),vx=profile*u,vy=profile*v,vz=0.;
         for(int a=0;a<k;++a){auto off=std::size_t(a)*N;double ap=alpha[a]*solution[off+t];
           double ddx=(solution[off+idx(y,(x+1)%nx)]-solution[off+idx(y,(x+nx-1)%nx)])/(2*dx);
           double ddy=(solution[off+idx((y+1)%ny,x)]-solution[off+idx((y+ny-1)%ny,x)])/(2*dx);
           double w=std::exp(-alpha[a]*hh);vx+=w*(ddx+gx[t]*ap);vy+=w*(ddy+gy[t]*ap);vz-=w*ap;
         }
         values(j,0)+=ws[corner]*vx;values(j,1)+=ws[corner]*vy;values(j,2)+=ws[corner]*vz;
       }
     }
     return Rcpp::List::create(Rcpp::_["iterations"]=iteration,Rcpp::_["relative_residual"]=residual,Rcpp::_["targets"]=values);
   }
   const double profile=std::log1p(height/z0)/std::log1p(reference/z0);Vec vx(N,profile*u),vy(N,profile*v),vz(N,0.);double boundary_sum=0,boundary_max=0;
   for(int y=0;y<ny;++y)for(int x=0;x<nx;++x){auto t=idx(y,x);double boundary=0;
     for(int a=0;a<k;++a){auto off=std::size_t(a)*N;double ap=alpha[a]*solution[off+t];
       double ddx=(solution[off+idx(y,(x+1)%nx)]-solution[off+idx(y,(x+nx-1)%nx)])/(2*dx);
       double ddy=(solution[off+idx((y+1)%ny,x)]-solution[off+idx((y+ny-1)%ny,x)])/(2*dx);
       double weight=std::exp(-alpha[a]*height);vx[t]+=weight*(ddx+gx[t]*ap);vy[t]+=weight*(ddy+gy[t]*ap);vz[t]-=weight*ap;
       boundary-=hx[t]*ddx+hy[t]*ddy+normal[t]*ap;
     }
     boundary_sum+=boundary*boundary;boundary_max=std::max(boundary_max,std::abs(boundary));
   }
   // Same first-order target sampling as scipy.ndimage.map_coordinates.
   int y0=int(std::floor(row)),x0=int(std::floor(col));
   int y1=std::min(y0+1,ny-1),x1=std::min(x0+1,nx-1);
   double wy=row-y0,wx=col-x0;
   auto sample=[&](const Vec& field){return (1-wy)*((1-wx)*field[idx(y0,x0)]+wx*field[idx(y0,x1)])+wy*((1-wx)*field[idx(y1,x0)]+wx*field[idx(y1,x1)]);};
   return Rcpp::List::create(Rcpp::_["iterations"]=iteration,Rcpp::_["relative_residual"]=residual,
     Rcpp::_["target"]=Rcpp::NumericVector::create(sample(vx),sample(vy),sample(vz)),
     Rcpp::_["boundary_rms"]=std::sqrt(boundary_sum/N),Rcpp::_["boundary_max"]=boundary_max);
 }
};
}
// [[Rcpp::export]]
Rcpp::List modal_native(Rcpp::NumericVector h,int ny,int nx,double dx,Rcpp::NumericVector alpha,
 Rcpp::NumericVector C,Rcpp::NumericVector P,Rcpp::NumericVector lambda,Rcpp::NumericVector beta,
 double penalty,double u,double v,double height,double z0,double reference,double row,double col,
 double rtol=1e-6,int maxiter=600,Rcpp::Nullable<Rcpp::NumericVector> test=R_NilValue,Rcpp::Nullable<Rcpp::NumericMatrix> targets=R_NilValue){
 int k=alpha.size();
 if(ny<3||nx<3||h.size()!=double(ny)*nx||k<1||C.size()!=k*k||P.size()!=k*k||lambda.size()!=k||beta.size()!=k||row<0||row>=ny||col<0||col>=nx||dx<=0||rtol<=0||maxiter<1)Rcpp::stop("Invalid shape or solver arguments");
 for(double x:h)if(!std::isfinite(x))Rcpp::stop("Non-finite terrain");
 for(double x:{dx,penalty,u,v,height,z0,reference,row,col,rtol})if(!std::isfinite(x))Rcpp::stop("Non-finite scalar argument");
 if(penalty<0||height<=0||z0<=0||reference<=0||row>ny-1||col>nx-1)Rcpp::stop("Invalid physical or sampling arguments");
 for(double x:alpha)if(!std::isfinite(x)||x<=0)Rcpp::stop("Invalid inverse scale");
 for(double x:lambda)if(!std::isfinite(x)||x<=0)Rcpp::stop("Invalid preconditioner eigenvalue");
 for(double x:C)if(!std::isfinite(x))Rcpp::stop("Non-finite modal metric");
 for(double x:P)if(!std::isfinite(x))Rcpp::stop("Non-finite preconditioner");
 for(double x:beta)if(!std::isfinite(x))Rcpp::stop("Non-finite forcing coefficient");
 auto vec=[](Rcpp::NumericVector x){return Vec(x.begin(),x.end());};
 Solver solver(vec(h),k,ny,nx,dx,penalty,vec(C),vec(alpha),vec(P),vec(lambda),vec(beta));
 if(test.isNotNull()){
   Rcpp::NumericVector value(test);if(value.size()!=double(k)*ny*nx)Rcpp::stop("Test vector size mismatch");
   Vec x=vec(value),a(x.size()),p(x.size());solver.apply(x,a);solver.precondition(x,p);
   return Rcpp::List::create(Rcpp::_["apply"]=Rcpp::wrap(a),Rcpp::_["precondition"]=Rcpp::wrap(p));
 }
 return solver.solve(u,v,rtol,maxiter,height,z0,reference,row,col,targets);
}
