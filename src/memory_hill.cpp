// Bolund supplied-boundary adapter for the unchanged exact memory integral.
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
// [[Rcpp::plugins(cpp17)]]
// Experimental native version of the existing exact raster-ray memory.
// Ascending cell-centre coordinates, rows south to north; no new physics.
// [[Rcpp::export]]
Rcpp::List memory_rose_native(Rcpp::NumericMatrix roughness,double x0,double y0,
 double dx,double dy,double x,double y,double source_z0,double height=10,
 double reference=100,double limit=10000,int sectors=72,
 Rcpp::Nullable<Rcpp::NumericVector> requested_directions=R_NilValue,
 bool supplied_boundary=false,double outside_z0=.0003,double coast_x=327,double coast_z0=.015){
 int ny=roughness.nrow(),nx=roughness.ncol();
 Rcpp::NumericVector requested;
 if(requested_directions.isNotNull()){requested=Rcpp::NumericVector(requested_directions);sectors=requested.size();}
 if(supplied_boundary&&(!std::isfinite(outside_z0)||!std::isfinite(coast_z0)||!std::isfinite(coast_x)||outside_z0<=0||coast_z0<=0||outside_z0>=std::min(height,reference)||coast_z0>=std::min(height,reference)))Rcpp::stop("Invalid supplied boundary");
 for(double a:{x0,y0,dx,dy,x,y,source_z0,height,reference,limit})if(!std::isfinite(a))Rcpp::stop("Non-finite input");
 if(nx<2||ny<2||dx<=0||dy<=0||source_z0<=0||height<=source_z0||reference<=source_z0||limit<=0||sectors<1)Rcpp::stop("Invalid geometry/profile");
 for(double a:roughness)if(!std::isfinite(a)||a<=0||a>=std::min(height,reference))Rcpp::stop("Invalid roughness cell");
 double left=x0-dx/2,bottom=y0-dy/2,right=left+nx*dx,top=bottom+ny*dy;
 if(x<left||x>right||y<bottom||y>top)Rcpp::stop("Target outside map");
 auto eq=[&](double z){return std::log(height/z)/std::log(reference/z);};
 const double denominator=eq(source_z0),pi=std::acos(-1.);
 Rcpp::NumericVector factors(sectors),tails(sectors),directions(sectors);
 for(int d=0;d<sectors;++d){
   double logsum=0,tailmax=0;directions[d]=requested_directions.isNotNull()?requested[d]:360.*d/sectors;
   if(!std::isfinite(directions[d]))Rcpp::stop("Invalid direction");
   for(double offset:{-2.5,0.,2.5}){
     double angle=(directions[d]+offset)*pi/180,vx=std::sin(angle),vy=std::cos(angle);
     if(std::abs(vx)<1e-14)vx=0;
     if(std::abs(vy)<1e-14)vy=0;
     double end=limit;
     if(!supplied_boundary){
     if(vx>0)end=std::min(end,(right-x)/vx);else if(vx<0)end=std::min(end,(left-x)/vx);
     if(vy>0)end=std::min(end,(top-y)/vy);else if(vy<0)end=std::min(end,(bottom-y)/vy);
     }
     if(end<=0)Rcpp::stop("No upwind support");
     std::vector<double> edges{0.,end};edges.reserve(nx+ny+4);
     if(vx)for(int i=0;i<=nx;++i){double t=(left+i*dx-x)/vx;if(t>0&&t<end)edges.push_back(t);}
     if(vy)for(int i=0;i<=ny;++i){double t=(bottom+i*dy-y)/vy;if(t>0&&t<end)edges.push_back(t);}
     if(supplied_boundary&&vx){double t=(coast_x-x)/vx;if(t>0&&t<end)edges.push_back(t);}
     std::sort(edges.begin(),edges.end());
     // Corner crossings can differ by one rounding unit. Do not integrate a
     // spurious zero-length segment just beyond the clipped map boundary.
     const double edge_tol=64*std::numeric_limits<double>::epsilon()*std::max(1.,end);
     edges.erase(std::unique(edges.begin(),edges.end(),[edge_tol](double a,double b){
       return std::abs(a-b)<=edge_tol;
     }),edges.end());
     if(edges.size()<2)Rcpp::stop("Insufficient positive ray support");
     if(std::abs(edges.back()-end)<=edge_tol)edges.back()=end;
     double scale=0,mem=0;
     for(std::size_t i=0;i+1<edges.size();++i){
       double middle=(edges[i]+edges[i+1])/2;
       int col=int(std::floor((x+vx*middle-left)/dx)),row=int(std::floor((y+vy*middle-bottom)/dy));
       bool inside=col>=0&&col<nx&&row>=0&&row<ny;
       if(!inside&&!supplied_boundary)Rcpp::stop("Ray midpoint outside clipped support");
       double z=inside?roughness(row,col):outside_z0;
       if(supplied_boundary&&x+vx*middle>coast_x)z=coast_z0;
       if(i==0)scale=z*std::pow(height/(.75*z),1.25);
       // Last mapped segment continues to infinity, identical to Python.
       double survival=std::exp(-edges[i]/scale);
       double next=i+2==edges.size()?0.:std::exp(-edges[i+1]/scale);
       mem+=(survival-next)*eq(z);
     }
     logsum+=std::log(mem);tailmax=std::max(tailmax,std::exp(-end/scale));
   }
   factors[d]=std::exp(logsum/3)/denominator;tails[d]=tailmax;
 }
 return Rcpp::List::create(Rcpp::_["direction"]=directions,Rcpp::_["factor"]=factors,Rcpp::_["tail_weight"]=tails);
}

// [[Rcpp::export]]
Rcpp::NumericVector wind_apply_native(Rcpp::NumericMatrix u,Rcpp::NumericMatrix v,
 Rcpp::NumericVector weights,Rcpp::NumericVector matrix,Rcpp::NumericVector rose,int mode=2){
 int n=u.nrow(),nr=rose.size();
 if(u.ncol()!=4||v.nrow()!=n||v.ncol()!=4||weights.size()!=4||matrix.size()!=4||nr<1||mode<0||mode>3)Rcpp::stop("Invalid application shape or mode");
 for(double a:weights)if(!std::isfinite(a))Rcpp::stop("Non-finite weights");
 for(double a:matrix)if(!std::isfinite(a))Rcpp::stop("Non-finite transfer matrix");
 for(double a:rose)if(!std::isfinite(a)||a<0)Rcpp::stop("Invalid surface factors");
 Rcpp::NumericVector out(n);const double pi=std::acos(-1.);
 for(int i=0;i<n;++i){double a=0,b=0;
   for(int j=0;j<4;++j){a+=weights[j]*u(i,j);b+=weights[j]*v(i,j);}
   if(!std::isfinite(a)||!std::isfinite(b)){out[i]=NA_REAL;continue;}
   double value=(mode==0||mode==3)?std::hypot(a,b):std::hypot(matrix[0]*a+matrix[1]*b,matrix[2]*a+matrix[3]*b);
   if(mode==2||mode==3){double d=std::atan2(-a,-b)*180/pi;if(d<0)d+=360;
     double t=d*nr/360;int j=int(std::floor(t));double f=t-j;j%=nr;
     value*=rose[j]*(1-f)+rose[(j+1)%nr]*f;
   }
   out[i]=value;
 }
 return out;
}
