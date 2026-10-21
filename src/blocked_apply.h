// Same finite-difference modal operator, fused reusable work arrays.
// No changed physics, tolerances, parallel workers, or reduced modes.
#include <algorithm>
#include <cstddef>
extern "C" void modal_apply(int n,int ny,int nx,double dx,double penalty,
 const double* C,const double* alpha,const double* gx,const double* gy,
 const double* hx,const double* hy,const double* normal,const double* p,
 double* cp,double* cap,double* fx,double* fy,double* total,double* weighted,
 double* q,double* qx,double* qy,double* out) {
 const std::size_t N=std::size_t(ny)*nx;const double inv=1/(2*dx);
 std::fill(cp,cp+n*N,0.);std::fill(cap,cap+n*N,0.);
 std::fill(total,total+N,0.);std::fill(weighted,weighted+N,0.);
 for(std::size_t start=0;start<N;start+=1024){auto end=std::min(N,start+1024);
 for(int j=0;j<n;++j) {
   const double* v=p+j*N;const double aj=alpha[j];
   for(std::size_t k=start;k<end;++k){total[k]+=v[k];weighted[k]+=aj*v[k];}
   for(int i=0;i<n;++i){
     double* c=cp+i*N;double* a=cap+i*N;const double cij=C[i*n+j],ca=cij*aj;
     for(std::size_t k=start;k<end;++k){c[k]+=cij*v[k];a[k]+=ca*v[k];}
   }
 }}
 for(int y=0;y<ny;++y)for(int x=0;x<nx;++x){
   const std::size_t k=std::size_t(y)*nx+x;
   const auto xp=std::size_t(y)*nx+(x+1==nx?0:x+1),xm=std::size_t(y)*nx+(x==0?nx-1:x-1);
   const auto yp=std::size_t(y+1==ny?0:y+1)*nx+x,ym=std::size_t(y==0?ny-1:y-1)*nx+x;
   q[k]=-hx[k]*(total[xp]-total[xm])*inv-hy[k]*(total[yp]-total[ym])*inv-normal[k]*weighted[k];
   qx[k]=q[k]*hx[k];qy[k]=q[k]*hy[k];
   for(int i=0;i<n;++i){const auto o=std::size_t(i)*N;
     fx[o+k]=(cp[o+xp]-cp[o+xm])*inv+gx[k]*cap[o+k];
     fy[o+k]=(cp[o+yp]-cp[o+ym])*inv+gy[k]*cap[o+k];}
 }
 for(int y=0;y<ny;++y)for(int x=0;x<nx;++x){
   const std::size_t k=std::size_t(y)*nx+x;
   const auto xp=std::size_t(y)*nx+(x+1==nx?0:x+1),xm=std::size_t(y)*nx+(x==0?nx-1:x-1);
   const auto yp=std::size_t(y+1==ny?0:y+1)*nx+x,ym=std::size_t(y==0?ny-1:y-1)*nx+x;
   const double boundary=((qx[xp]-qx[xm])+(qy[yp]-qy[ym]))*inv;
   for(int i=0;i<n;++i){const auto o=std::size_t(i)*N;const double a=alpha[i];
     out[o+k]=-((fx[o+xp]-fx[o+xm])+(fy[o+yp]-fy[o+ym]))*inv
       +a*(gx[k]*fx[o+k]+gy[k]*fy[o+k]+cap[o+k])
       +penalty*(boundary-a*q[k]*normal[k]);}
 }
}
