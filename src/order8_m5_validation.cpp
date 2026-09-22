// Standalone DSMC/weak-form validation used for Figure 2 (Mach-5 shock).
//
// The two active_* functions are the single-generator theorem evaluations of
// the collisional increments for the order-eight Burnett observables
// (j,l)=(2,4) and (3,2). The angular loop near the end independently samples
// the conventional hard-sphere weak-form collision bracket. The historical
// output column name "tree" means "single-generator theorem prediction".
// No manuscript or supplementary-document source is required by this file.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>
using std::array; using std::vector;
struct P { double x; array<double,3> v; };
static constexpr double PI=3.1415926535897932384626433832795;
// PRL/Nature dimensionless outgoing-direction normalization.
// The DSMC dynamics below retains the physical diameter-one cross section PI;
// only the reported collision-production estimators use ETA_HS.
static constexpr double ETA_HS=0.88622692545275801364908374167057; // sqrt(pi)/2
static constexpr double LBOX=1.0;
static constexpr int NC=50;
static constexpr double DX=LBOX/NC;
static constexpr double GAM=5.0/3.0;
static constexpr double KN1=0.05;
static constexpr double N1=1.0/(std::sqrt(2.0)*PI*KN1);
static constexpr double T1=1.0;
inline double npdf(double z){ return std::exp(-0.5*z*z)/std::sqrt(2.0*PI); }
inline double ncdf(double z){ return 0.5*std::erfc(-z/std::sqrt(2.0)); }
double flux_plus(double U,double T){ double s=std::sqrt(T),a=U/s; return U*ncdf(a)+s*npdf(a); }
double flux_cdf_v(double v,double U,double T){ double s=std::sqrt(T),y0=-U/s,y=(v-U)/s; double num=U*(ncdf(y)-ncdf(y0))+s*(npdf(y0)-npdf(y)); return num/flux_plus(U,T); }
template<class RNG> double sample_flux_positive(double U,double T,RNG& rng){ std::uniform_real_distribution<double> unif(0,1); double q=unif(rng),s=std::sqrt(T),lo=0,hi=std::max(1.0,U+10*s); while(flux_cdf_v(hi,U,T)<q)hi*=2; for(int it=0;it<60;++it){double m=.5*(lo+hi); if(flux_cdf_v(m,U,T)<q)lo=m;else hi=m;} return .5*(lo+hi); }
template<class RNG> array<double,3> sample_incoming(bool left,double U,double T,RNG& rng){std::normal_distribution<double> gauss(0,std::sqrt(T)); array<double,3>v; v[0]=left?sample_flux_positive(U,T,rng):-sample_flux_positive(-U,T,rng);v[1]=gauss(rng);v[2]=gauss(rng);return v;}
inline double dot3(const array<double,3>&a,const array<double,3>&b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline double norm3(const array<double,3>&a){return std::sqrt(dot3(a,a));}
template<class RNG> void isotropic(array<double,3>&n,RNG&rng){std::uniform_real_distribution<double>u(0,1);double z=2*u(rng)-1,ph=2*PI*u(rng),r=std::sqrt(std::max(0.0,1-z*z));n={r*std::cos(ph),r*std::sin(ph),z};}

double solid_H(int l,const array<double,3>&c){double r2=dot3(c,c),xx=c[0];if(l==0)return 1;double Rm=1,R=xx;if(l==1)return R;for(int k=1;k<l;++k){double Rp=((2.0*k+1)*xx*R-k*r2*Rm)/(k+1.0);Rm=R;R=Rp;}double df=1;for(int k=1;k<=l;++k)df*=2*k-1;return R/std::sqrt(df);}
double laguerre(int j,double a,double t){if(j==0)return 1;if(j==1)return 1+a-t;double Lm=1,L=1+a-t;for(int k=1;k<j;++k){double Lp=((2.0*k+1+a-t)*L-(k+a)*Lm)/(k+1.0);Lm=L;L=Lp;}return L;}
double burnett(int j,int l,const array<double,3>&c){double s=dot3(c,c),a=l+.5;double lc=.5*(std::lgamma(j+1.0)+std::lgamma(a+1)-std::lgamma(a+j+1));double cj=((j%2)?-1:1)*std::exp(lc);return cj*laguerre(j,a,s/2)*solid_H(l,c);}

inline double active_24(double Cx,double Cy,double Cz,double Gx,double Gy,double Gz){
  const double x0 = pow(Gx, 4);
  const double x1 = pow(Gx, 6);
  const double x2 = pow(Gy, 4);
  const double x3 = pow(Gy, 6);
  const double x4 = pow(Gz, 4);
  const double x5 = pow(Gz, 6);
  const double x6 = Cy*Gy;
  const double x7 = Cx*Gx;
  const double x8 = 13178880*x7;
  const double x9 = Cz*Gz;
  const double x10 = pow(Cx, 2);
  const double x11 = pow(Gx, 2);
  const double x12 = x0*x10;
  const double x13 = pow(Gy, 2);
  const double x14 = 3294720*x10;
  const double x15 = 1080*x10;
  const double x16 = pow(Gz, 2);
  const double x17 = pow(Cx, 4);
  const double x18 = x11*x17;
  const double x19 = 27184*x17;
  const double x20 = pow(Cx, 6);
  const double x21 = 74240*x20;
  const double x22 = pow(Cy, 2);
  const double x23 = 3294720*x11;
  const double x24 = 480*x1;
  const double x25 = x2*x22;
  const double x26 = x22*x4;
  const double x27 = pow(Cy, 4);
  const double x28 = 15744*x0;
  const double x29 = x13*x27;
  const double x30 = pow(Cy, 6);
  const double x31 = 42240*x11;
  const double x32 = pow(Cz, 2);
  const double x33 = x2*x32;
  const double x34 = x32*x4;
  const double x35 = pow(Cz, 4);
  const double x36 = x16*x35;
  const double x37 = pow(Cz, 6);
  const double x38 = 823680*x11;
  const double x39 = 270*x11;
  const double x40 = 555*x0;
  const double x41 = 120*x1;
  const double x42 = x13*x4;
  const double x43 = x16*x2;
  const double x44 = pow(Gy, 3);
  const double x45 = 12960*x7;
  const double x46 = pow(Gy, 5);
  const double x47 = pow(Gx, 3);
  const double x48 = 5760*Cx*pow(Gx, 5);
  const double x49 = pow(Cy, 3);
  const double x50 = 207360*x7;
  const double x51 = pow(Cy, 5);
  const double x52 = pow(Gz, 3);
  const double x53 = pow(Gz, 5);
  const double x54 = pow(Cz, 3);
  const double x55 = pow(Cz, 5);
  const double x56 = pow(Cx, 3);
  const double x57 = 92160*pow(Cx, 5)*Gx;
  const double x58 = Cz*x52;
  const double x59 = x58*x6;
  const double x60 = Cy*x44;
  const double x61 = x60*x9;
  const double x62 = Gz*x54;
  const double x63 = x6*x62;
  const double x64 = Gy*x49;
  const double x65 = x64*x9;
  const double x66 = x6*x9;
  const double x67 = Cx*x47;
  const double x68 = 35520*x67;
  const double x69 = 172800*x7;
  const double x70 = 142080*x67;
  const double x71 = Gx*x56;
  const double x72 = 142080*x71;
  const double x73 = 76800*x47*x56;
  const double x74 = 568320*x71;
  const double x75 = 9888*x10;
  const double x76 = x16*x22;
  const double x77 = 1148160*x10;
  const double x78 = x10*x11;
  const double x79 = 38400*x78;
  const double x80 = 84480*x10;
  const double x81 = x13*x32;
  const double x82 = 13320*x78;
  const double x83 = 7200*x12;
  const double x84 = 3240*x10;
  const double x85 = x13*x22;
  const double x86 = 80000*x17;
  const double x87 = x16*x32;
  const double x88 = 65408*x18;
  const double x89 = x13*x16;
  const double x90 = x22*x32;
  const double x91 = x22*x81;
  const double x92 = x32*x76;
  const double x93 = 126720*x11;
  const double x94 = 16200*x11;
  const double x95 = 3240*x11;
  const double x96 = 13320*x0;
  const double x97 = 4440*x0;
  const double x98 = 336960*x89;
  const double x99 = 78528*x11;
  const double x100 = 26688*x11;
  const double x101 = 810*x11;
  const double x102 = 103680*x7;
  const double x103 = x32*x6;
  const double x104 = 25920*x7;
  const double x105 = x16*x6;
  const double x106 = x22*x9;
  const double x107 = 414720*x7;
  const double x108 = x13*x9;
  const double x109 = 103680*x10;
  const double x110 = 414720*x10;
  const double x111 = 25920*x11;
  const double x112 = 103680*x11;
  const double x113 = 311040*x7;
  const double x114 = 209856*x78;
  const double x115 = 105216*x11;
  const double x116 = 19440*x11*x89;
  return (1.0/922521600.0)*sqrt(30030)*(1048320*Cx*Cy*Gx*Gy*x16 + 4193280*Cx*Cy*Gx*Gy*x32 + 1048320*Cx*Cy*Gx*x44 + 798720*Cx*Cy*Gy*x47 + 1048320*Cx*Cz*Gx*Gz*x13 + 4193280*Cx*Cz*Gx*Gz*x22 + 1048320*Cx*Cz*Gx*x52 + 798720*Cx*Cz*Gz*x47 + 4193280*Cx*Gx*Gy*x49 + 4193280*Cx*Gx*Gz*x54 + 4193280*Cy*Cz*Gy*Gz*x10 + 1048320*Cy*Cz*Gy*Gz*x11 + 3294720*Cy*Cz*Gy*Gz + 8640*Cy*Cz*Gy*x53 + 8640*Cy*Cz*Gz*x46 + 17280*Cy*Cz*x44*x52 + 3194880*Cy*Gx*Gy*x56 + 138240*Cy*Gy*Gz*x55 + 115200*Cy*Gy*x52*x54 + 69120*Cy*Gz*x44*x54 - Cy*x45*x46 + 3194880*Cz*Gx*Gz*x56 + 138240*Cz*Gy*Gz*x51 + 69120*Cz*Gy*x49*x52 + 115200*Cz*Gz*x44*x49 - Cz*x45*x53 + 120*pow(Gx, 8) + 45*pow(Gy, 8) + 276480*Gy*Gz*x49*x54 - Gy*x50*x51 + 45*pow(Gz, 8) - Gz*x50*x55 + 26112*x0*x10*x22 + 26112*x0*x10*x32 + 24960*x0*x13 + 24960*x0*x16 + 116096*x0*x17 + 99840*x0*x22 + 99840*x0*x32 - 17760*x0*x66 - 1110*x0*x89 - 31488*x0*x90 + 274560*x0 + 13440*x1*x10 - 12480*x1 + 3264*x10*x11*x13*x32 + 599040*x10*x11*x13 + 3264*x10*x11*x16*x22 + 599040*x10*x11*x16 + 199680*x10*x11*x22 + 199680*x10*x11*x32 + 6589440*x10*x11 + 32064*x10*x13*x16*x22 + 32064*x10*x13*x16*x32 + 262080*x10*x13*x16 + 38400*x10*x13*x22*x32 + 948480*x10*x13*x22 + 122880*x10*x13*x35 + 38400*x10*x16*x22*x32 + 122880*x10*x16*x27 + 948480*x10*x16*x32 + 41952*x10*x2*x32 + 131040*x10*x2 + 41952*x10*x22*x4 + 131040*x10*x4 - x100*x13*x35 - x100*x16*x27 - x101*x42 - x101*x43 - x102*x13*x62 - x102*x16*x64 - x102*x22*x58 - x102*x32*x60 - x103*x70 - x103*x74 - x104*x13*x58 - x104*x16*x60 - x105*x68 - x105*x72 - x106*x70 - x106*x74 - x107*x22*x62 - x107*x32*x64 - x108*x68 - x108*x72 - x109*x59 - x109*x61 + 65520*x11*x13*x16 + 786240*x11*x13*x22 + 262080*x11*x13*x32 + 262080*x11*x16*x22 + 786240*x11*x16*x32 + 17920*x11*x17*x22 + 17920*x11*x17*x32 + 32760*x11*x2 + 148480*x11*x20 + 1597440*x11*x22*x32 + 798720*x11*x27 + 798720*x11*x35 + 32760*x11*x4 - x110*x63 - x110*x65 - x111*x59 - x111*x61 - x112*x63 - x112*x65 - x113*x6*x87 - x113*x85*x9 - x114*x85 - x114*x87 - x115*x91 - x115*x92 - x116*x22 - x116*x32 - 748800*x12 - x13*x14 + 128064*x13*x16*x22*x32 + 29472*x13*x16*x27 + 29472*x13*x16*x35 + 205920*x13*x16 + 62080*x13*x17*x32 + 1131520*x13*x17 - x13*x21 + 28800*x13*x22*x35 + 6480*x13*x22*x4 + 2471040*x13*x22 + 97920*x13*x27*x32 + 55680*x13*x30 + 10800*x13*x32*x4 + 823680*x13*x32 + 49920*x13*x35 - 13440*x13*x37 - x13*x38 - x13*x41 + 180*x13*x5 - x13*x83 - x13*x88 - x14*x16 - x15*x3 - x15*x5 + 62080*x16*x17*x22 + 1131520*x16*x17 + 10800*x16*x2*x22 + 6480*x16*x2*x32 - x16*x21 + 97920*x16*x22*x35 + 823680*x16*x22 + 28800*x16*x27*x32 + 49920*x16*x27 + 180*x16*x3 - 13440*x16*x30 + 2471040*x16*x32 + 55680*x16*x37 - x16*x38 - x16*x41 - x16*x83 - x16*x88 - 284160*x17*x66 - 54368*x17*x89 - 2263040*x18 - x19*x2 - x19*x4 + 29472*x2*x22*x32 + 43536*x2*x27 - 2544*x2*x35 + 270*x2*x4 - x2*x40 - x2*x45*x9 - x2*x82 + 102960*x2 - x22*x23 - x22*x24 + 5040*x22*x3 + 29472*x22*x32*x4 - x22*x35*x93 + 720*x22*x5 - x22*x98 - x23*x32 - x24*x32 - x25*x75 - x25*x94 - 280800*x25 - x26*x95 - 56160*x26 - x27*x28 - x27*x32*x93 - 2544*x27*x4 - x27*x50*x9 - x27*x79 - x28*x35 - x29*x80 - x29*x99 - 848640*x29 + 720*x3*x32 - x3*x39 - 4680*x3 - x30*x31 - x31*x37 + 5040*x32*x5 - x32*x98 - x33*x95 - 56160*x33 - x34*x75 - x34*x94 - 280800*x34 + 43536*x35*x4 - x35*x50*x6 - x35*x79 - x36*x80 - x36*x99 - 848640*x36 - x39*x5 - x4*x40 - x4*x45*x6 - x4*x82 + 102960*x4 - x42*x84 - 14040*x42 - x43*x84 - 14040*x43 - x44*x49*x69 - x48*x6 - x48*x9 - 4680*x5 - x52*x54*x69 - x57*x6 - x57*x9 - x58*x68 - x58*x72 - 449280*x59 - x6*x73 - x6*x8 - x60*x68 - x60*x72 - 449280*x61 - x62*x70 - x62*x74 - 1797120*x63 - x64*x70 - x64*x74 - 1797120*x65 - 426240*x66*x78 - x73*x9 - x76*x77 - x76*x97 - x77*x81 - 26640*x78*x89 - 76800*x78*x90 - x8*x9 - x81*x97 - x85*x86 - x85*x96 - x86*x87 - x87*x96 - 798720*x91 - 798720*x92);
}
inline double active_32(double Cx,double Cy,double Cz,double Gx,double Gy,double Gz){
  const double x0 = pow(Gx, 2);
  const double x1 = pow(Gx, 4);
  const double x2 = pow(Gx, 6);
  const double x3 = pow(Gy, 2);
  const double x4 = pow(Gy, 4);
  const double x5 = pow(Gy, 6);
  const double x6 = pow(Gz, 2);
  const double x7 = pow(Gz, 4);
  const double x8 = pow(Gz, 6);
  const double x9 = Cy*Gy;
  const double x10 = Cz*Gz;
  const double x11 = pow(Cx, 2);
  const double x12 = 348480*x11;
  const double x13 = 284*x11;
  const double x14 = pow(Cx, 4);
  const double x15 = x0*x14;
  const double x16 = 2928*x14;
  const double x17 = pow(Cx, 6);
  const double x18 = 5440*x17;
  const double x19 = pow(Cy, 2);
  const double x20 = 19008*x1;
  const double x21 = x19*x3;
  const double x22 = x19*x7;
  const double x23 = pow(Cy, 4);
  const double x24 = 63360*x0;
  const double x25 = pow(Cy, 6);
  const double x26 = pow(Cz, 2);
  const double x27 = x26*x4;
  const double x28 = x26*x6;
  const double x29 = pow(Cz, 4);
  const double x30 = pow(Cz, 6);
  const double x31 = 5*x0;
  const double x32 = 1980*x1;
  const double x33 = x3*x6;
  const double x34 = Cy*pow(Gy, 5);
  const double x35 = Cx*Gx;
  const double x36 = 240*x35;
  const double x37 = pow(Gx, 3);
  const double x38 = Cx*x37;
  const double x39 = pow(Gx, 5);
  const double x40 = pow(Cy, 5)*Gy;
  const double x41 = 3840*x35;
  const double x42 = Cz*pow(Gz, 5);
  const double x43 = pow(Cz, 5)*Gz;
  const double x44 = pow(Cx, 3);
  const double x45 = Gx*x44;
  const double x46 = pow(Cx, 5);
  const double x47 = pow(Gz, 3);
  const double x48 = pow(Gy, 3);
  const double x49 = pow(Cz, 3);
  const double x50 = pow(Cy, 3);
  const double x51 = x48*x50;
  const double x52 = 3200*x35;
  const double x53 = x47*x49;
  const double x54 = Cz*x47;
  const double x55 = Cy*x48;
  const double x56 = 1920*x55;
  const double x57 = Gz*x49;
  const double x58 = 7680*x57;
  const double x59 = Gy*x50;
  const double x60 = 7680*x59;
  const double x61 = x0*x11;
  const double x62 = 126720*x61;
  const double x63 = 2784*x11;
  const double x64 = 1824*x11;
  const double x65 = x23*x3;
  const double x66 = 4800*x11;
  const double x67 = 960*x11;
  const double x68 = x29*x6;
  const double x69 = 3168*x61;
  const double x70 = 432*x61;
  const double x71 = x3*x7;
  const double x72 = 852*x11;
  const double x73 = x4*x6;
  const double x74 = 4800*x14;
  const double x75 = 8640*x14;
  const double x76 = x19*x6;
  const double x77 = x26*x3;
  const double x78 = 336*x15;
  const double x79 = x0*x19;
  const double x80 = 22176*x0;
  const double x81 = 1632*x33;
  const double x82 = 15*x0;
  const double x83 = x26*x35;
  const double x84 = 480*x35;
  const double x85 = x19*x35;
  const double x86 = 1920*x54;
  const double x87 = 1920*x35;
  const double x88 = x11*x9;
  const double x89 = x10*x11;
  const double x90 = 480*x0;
  const double x91 = 1920*x0;
  const double x92 = 5760*x35;
  const double x93 = 5760*x11;
  const double x94 = 288*x61;
  const double x95 = 4608*x11;
  return (1.0/3548160.0)*sqrt(77)*(380160*Cx*Cy*Gx*Gy + 3840*Cx*Cy*Gy*x26*x37 + 960*Cx*Cy*Gy*x37*x6 + 1200*Cx*Cy*Gy*x39 + 960*Cx*Cy*x37*x48 + 380160*Cx*Cz*Gx*Gz + 3840*Cx*Cz*Gz*x19*x37 + 960*Cx*Cz*Gz*x3*x37 + 1200*Cx*Cz*Gz*x39 + 960*Cx*Cz*x37*x47 + 3840*Cx*Gy*x37*x50 + 3840*Cx*Gz*x37*x49 + 11520*Cy*Cz*Gy*Gz*x0*x11 + 480*Cy*Cz*Gy*Gz*x1 + 7680*Cy*Cz*Gy*Gz*x14 + 63360*Cy*Cz*Gy*x47 + 63360*Cy*Cz*Gz*x48 + 15360*Cy*Gx*Gy*x26*x44 + 3840*Cy*Gx*Gy*x44*x6 + 19200*Cy*Gx*Gy*x46 + 3840*Cy*Gx*x44*x48 + 253440*Cy*Gy*Gz*x49 + 16000*Cy*Gy*x37*x44 + 15360*Cz*Gx*Gz*x19*x44 + 3840*Cz*Gx*Gz*x3*x44 + 19200*Cz*Gx*Gz*x46 + 3840*Cz*Gx*x44*x47 + 253440*Cz*Gy*Gz*x50 + 16000*Cz*Gz*x37*x44 + 10*pow(Gx, 8) + 15360*Gx*Gy*x44*x50 + 15360*Gx*Gz*x44*x49 - 5*pow(Gy, 8) - 5*pow(Gz, 8) + 11520*x0*x11*x19*x26 + 5472*x0*x11*x19*x3 + 5760*x0*x11*x23 + 5472*x0*x11*x26*x6 + 5760*x0*x11*x29 + 696960*x0*x11 + 13440*x0*x14*x19 + 13440*x0*x14*x26 + 10880*x0*x17 + 4896*x0*x19*x26*x3 + 4896*x0*x19*x26*x6 + 9600*x0*x19*x29 + 432*x0*x19*x3*x6 + 96*x0*x19*x4 + 336*x0*x19*x7 + 316800*x0*x19 + 9600*x0*x23*x26 + 1968*x0*x23*x3 + 2928*x0*x23*x6 + 3200*x0*x25 + 432*x0*x26*x3*x6 + 336*x0*x26*x4 + 96*x0*x26*x7 + 316800*x0*x26 + 2928*x0*x29*x3 + 1968*x0*x29*x6 + 23760*x0*x3 + 3200*x0*x30 + 23760*x0*x6 - 443520*x0 + 4416*x1*x11*x19 + 4416*x1*x11*x26 + 708*x1*x11*x3 + 708*x1*x11*x6 - 57024*x1*x11 + 8032*x1*x14 + 3648*x1*x19*x26 + 756*x1*x19*x3 + 516*x1*x19*x6 + 1824*x1*x23 + 516*x1*x26*x3 + 756*x1*x26*x6 + 1824*x1*x29 + 30*x1*x3*x6 + 15*x1*x4 + 15*x1*x7 + 47520*x1 - x10*x21*x92 - x10*x23*x41 - 960*x10*x34 - x10*x36*x4 - 63360*x10*x38 - 15360*x10*x40 - 253440*x10*x45 - 12800*x10*x51 - x10*x55*x90 - x10*x59*x91 - 760320*x10*x9 + 63360*x11*x19*x3 + 63360*x11*x19*x6 + 856*x11*x2 + 63360*x11*x26*x3 + 63360*x11*x26*x6 + 44352*x11*x3*x6 + 22176*x11*x4 + 22176*x11*x7 - x12*x3 - x12*x6 - x13*x5 - x13*x8 + 95040*x14*x3 - 5856*x14*x33 + 95040*x14*x6 - 190080*x15 - x16*x4 - x16*x7 - x18*x3 - x18*x6 + 232*x19*x2 - x19*x20 + 63360*x19*x26*x3 + 63360*x19*x26*x6 - 1632*x19*x27 - x19*x28*x93 + 25344*x19*x3*x6 - x19*x4*x63 + 28512*x19*x4 - 428*x19*x5 + 31680*x19*x6 - x19*x62 - 8640*x19*x68 - 804*x19*x73 + 52*x19*x8 + 232*x2*x26 + 25*x2*x3 + 25*x2*x6 - 1320*x2 - x20*x26 - x21*x26*x93 - 10944*x21*x28 - 960*x21*x29 - x21*x6*x95 - 324*x21*x7 - x21*x74 - x21*x80 - 348480*x21 - 1632*x22*x26 - x22*x64 - 3168*x22 - x23*x24 - 960*x23*x28 + 95040*x23*x3 - 4016*x23*x4 - x23*x6*x67 - 31680*x23*x6 + 1104*x23*x7 - x23*x81 - x24*x29 - 5440*x25*x3 + 2240*x25*x6 + 25344*x26*x3*x6 + 31680*x26*x3 + 52*x26*x5 - x26*x62 - x26*x63*x7 - 8640*x26*x65 + 28512*x26*x7 - 804*x26*x71 - 126720*x26*x79 - 428*x26*x8 - x27*x64 - 3168*x27 - x28*x3*x95 - 324*x28*x4 - x28*x74 - x28*x80 - x28*x9*x92 - 348480*x28 - x29*x3*x67 - 31680*x29*x3 + 1104*x29*x4 - x29*x41*x9 + 95040*x29*x6 - 4016*x29*x7 - x29*x81 + 2240*x3*x30 - x3*x32 - x3*x54*x84 - x3*x57*x87 - x3*x69 + 1980*x3*x7 - x3*x78 - 20*x3*x8 + 221760*x3 - 5440*x30*x6 - x31*x5 - x31*x8 - x32*x6 - 864*x33*x61 - 47520*x33 - x34*x36 - x36*x42 - x36*x7*x9 - 63360*x38*x9 + 1980*x4*x6 - 30*x4*x7 - x4*x70 - 23760*x4 - x40*x41 - x41*x43 - 960*x42*x9 - 15360*x43*x9 - 253440*x45*x9 - 20*x5*x6 + 660*x5 - x51*x52 - x52*x53 - 12800*x53*x9 - x54*x56 - x54*x60 - x54*x9*x90 - x55*x58 - x55*x6*x84 - x56*x83 - x56*x89 - 30720*x57*x59 - x57*x9*x91 - x58*x85 - x58*x88 - x59*x6*x87 - x6*x69 - x6*x78 - 22176*x6*x79 + 221760*x6 - x60*x83 - x60*x89 - x65*x66 - x66*x68 - x7*x70 - 23760*x7 - x71*x72 - x71*x82 - x72*x73 - x73*x82 - x75*x76 - x75*x77 - x76*x94 - x77*x80 - x77*x94 + 660*x8 - x85*x86 - x86*x88);
}

struct Block {
  double direct[5][2]{};
  double tree[5][2]{};
  long long count[5][2]{};
  double relx[5]{};
  long long relcount[5]{};
};

int main(int argc,char**argv){
 if(argc<5){std::cerr<<"usage: exe Mach seed outprefix dt [burntime=1.5] [samptime=1.5] [Ntarget=18000] [pairs_per_station=300] [nang=12] [nblocks=20]\n";return 2;}
 double Mach=std::stod(argv[1]); int seed=std::stoi(argv[2]); std::string prefix=argv[3]; double dt=std::stod(argv[4]);
 double burntime=argc>5?std::stod(argv[5]):1.5, samptime=argc>6?std::stod(argv[6]):1.5; int Ntarget=argc>7?std::stoi(argv[7]):18000;
 int pairsPer=argc>8?std::stoi(argv[8]):300, nang=argc>9?std::stoi(argv[9]):12, nblocks=argc>10?std::stoi(argv[10]):20;
 double U1=Mach*std::sqrt(GAM*T1); double rr=((GAM+1)*Mach*Mach)/((GAM-1)*Mach*Mach+2); double p21=(2*GAM*Mach*Mach-(GAM-1))/(GAM+1); double T2=T1*p21/rr; double U2=U1/rr; double n2=N1*rr; double midrho=.5*(N1+n2);
 int burn=(int)std::llround(burntime/dt),samp=(int)std::llround(samptime/dt),sample_every=5; int nsamples=(samp+sample_every-1)/sample_every;
 std::mt19937_64 rng(0x9e3779b97f4a7c15ULL+1000003ULL*seed+((uint64_t)std::llround(Mach*1000))*9176ULL);std::uniform_real_distribution<double>unif(0,1);std::normal_distribution<double>stdn(0,1);
 double mass_est=0;for(int c=0;c<NC;++c){double xc=(c+.5)*DX,s=.5*(1+std::tanh((xc-.55)/.06));mass_est+=(N1+(n2-N1)*s)*DX;} double weight=mass_est/Ntarget;
 vector<P>ps;ps.reserve((size_t)(Ntarget*1.5));for(int c=0;c<NC;++c){double xc=(c+.5)*DX,s=.5*(1+std::tanh((xc-.55)/.06)),nn=N1+(n2-N1)*s,U=U1+(U2-U1)*s,T=T1+(T2-T1)*s;std::poisson_distribution<int>pois(nn*DX/weight);int np=pois(rng);for(int a=0;a<np;++a){P p;p.x=c*DX+unif(rng)*DX;double rt=std::sqrt(T);p.v={U+rt*stdn(rng),rt*stdn(rng),rt*stdn(rng)};ps.push_back(p);}}
 vector<vector<int>>bins(NC); vector<int> offsets={-6,-3,0,3,6}; vector<Block> blocks(nblocks); int sampleIndex=0;
 auto binp=[&](){for(auto&b:bins)b.clear();for(int i=0;i<(int)ps.size();++i){int c=(int)(ps[i].x/DX);if(c>=0&&c<NC)bins[c].push_back(i);}};
 auto inject=[&](double h){double JL=flux_plus(U1,T1),JR=flux_plus(-U2,T2);std::poisson_distribution<int>pL(N1*JL*h/weight),pR(n2*JR*h/weight);int nl=pL(rng),nr=pR(rng);for(int k=0;k<nl;++k){auto v=sample_incoming(true,U1,T1,rng);double age=unif(rng)*h;P p{v[0]*age,v};if(p.x<LBOX)ps.push_back(p);}for(int k=0;k<nr;++k){auto v=sample_incoming(false,U2,T2,rng);double age=unif(rng)*h;P p{LBOX+v[0]*age,v};if(p.x>=0)ps.push_back(p);}};
 auto stream=[&](double h){for(auto&p:ps)p.x+=p.v[0]*h;ps.erase(std::remove_if(ps.begin(),ps.end(),[](const P&p){return p.x<0||p.x>=LBOX;}),ps.end());inject(h);};
 auto collide=[&](){binp();for(int c=0;c<NC;++c){auto&b=bins[c];if(b.size()<2)continue;double sx=0,sy=0,sz=0;for(int ii:b){sx+=ps[ii].v[0];sy+=ps[ii].v[1];sz+=ps[ii].v[2];}array<double,3>uc={sx/b.size(),sy/b.size(),sz/b.size()};double vmax=0;for(int ii:b){array<double,3>dv={ps[ii].v[0]-uc[0],ps[ii].v[1]-uc[1],ps[ii].v[2]-uc[2]};vmax=std::max(vmax,norm3(dv));}double gmaj=std::max(1e-8,2.2*vmax);double meanCand=.5*b.size()*(b.size()-1)*weight*PI*gmaj*dt/DX;int ncand=(int)std::floor(meanCand);if(unif(rng)<meanCand-ncand)++ncand;std::uniform_int_distribution<int>pick(0,(int)b.size()-1);for(int cc=0;cc<ncand;++cc){int ia=pick(rng),ib=pick(rng);if(ia==ib)ib=(ib+1)%b.size();P&p=ps[b[ia]],&q=ps[b[ib]];array<double,3>Z,G;for(int d=0;d<3;++d){Z[d]=p.v[d]+q.v[d];G[d]=p.v[d]-q.v[d];}double g=norm3(G);if(g<1e-14||unif(rng)>std::min(1.0,g/gmaj))continue;array<double,3>n;isotropic(n,rng);for(int d=0;d<3;++d){double gg=g*n[d];p.v[d]=.5*(Z[d]+gg);q.v[d]=.5*(Z[d]-gg);}}}};
 auto sample=[&](){binp();vector<double>dens(NC);for(int c=0;c<NC;++c)dens[c]=bins[c].size()*weight/DX;double xm=.55,best=1e9;bool found=false;for(int c=0;c<NC-1;++c){double y0=dens[c]-midrho,y1=dens[c+1]-midrho;if(y0*y1<=0&&dens[c+1]!=dens[c]){double x0=(c+.5)*DX,x1=(c+1.5)*DX,xx=x0+(midrho-dens[c])*(x1-x0)/(dens[c+1]-dens[c]),e=std::abs(xx-.55);if(e<best){best=e;xm=xx;found=true;}}}if(!found)return;
    int bi=std::min(nblocks-1, sampleIndex*nblocks/std::max(1,nsamples)); ++sampleIndex; Block &B=blocks[bi];
    for(int kk=0;kk<5;++kk){double target=xm+offsets[kk]*DX;int c=(int)std::llround(target/DX-.5);c=std::max(0,std::min(NC-1,c));auto&b=bins[c];if(b.size()<20)continue;double sx=0,sy=0,sz=0;for(int ii:b){sx+=ps[ii].v[0];sy+=ps[ii].v[1];sz+=ps[ii].v[2];}array<double,3>u={sx/b.size(),sy/b.size(),sz/b.size()};double ss=0;for(int ii:b)for(int d=0;d<3;++d){double z=ps[ii].v[d]-u[d];ss+=z*z;}double T=ss/(3.0*b.size());if(T<=0)continue;double rt=std::sqrt(T), nn=b.size()*weight/DX, phys=nn*nn*rt;
      B.relx[kk]+=((c+.5)*DX-xm); B.relcount[kk]++;
      std::uniform_int_distribution<int>pick(0,(int)b.size()-1);
      for(int pp=0;pp<pairsPer;++pp){int ia=pick(rng),ib=pick(rng);if(ia==ib)ib=(ib+1)%b.size();array<double,3>c1,c2,C,G;for(int d=0;d<3;++d){c1[d]=(ps[b[ia]].v[d]-u[d])/rt;c2[d]=(ps[b[ib]].v[d]-u[d])/rt;C[d]=.5*(c1[d]+c2[d]);G[d]=c1[d]-c2[d];}double g=norm3(G);if(g<1e-12)continue;
        double act24=active_24(C[0],C[1],C[2],G[0],G[1],G[2]); double act32=active_32(C[0],C[1],C[2],G[0],G[1],G[2]);
        double tr24=-ETA_HS*g*act24*phys, tr32=-ETA_HS*g*act32*phys;
        double pre24=burnett(2,4,c1)+burnett(2,4,c2), pre32=burnett(3,2,c1)+burnett(3,2,c2);
        double d24=0,d32=0; for(int aa=0;aa<nang;++aa){array<double,3>n;isotropic(n,rng);array<double,3>p1,p2;for(int d=0;d<3;++d){double gg=g*n[d];p1[d]=C[d]+.5*gg;p2[d]=C[d]-.5*gg;}d24 += burnett(2,4,p1)+burnett(2,4,p2)-pre24; d32 += burnett(3,2,p1)+burnett(3,2,p2)-pre32;}d24/=nang;d32/=nang;
        B.tree[kk][0]+=tr24; B.direct[kk][0]+=ETA_HS*g*d24*phys; B.count[kk][0]++;
        B.tree[kk][1]+=tr32; B.direct[kk][1]+=ETA_HS*g*d32*phys; B.count[kk][1]++;
      }
    }
 };
 int total=burn+samp; auto t0=std::chrono::steady_clock::now(); for(int it=0;it<total;++it){stream(.5*dt);collide();stream(.5*dt);if(it>=burn&&((it-burn)%sample_every==0))sample();} auto t1c=std::chrono::steady_clock::now();
 std::ofstream fo(prefix+"_blocks.csv");fo<<"seed,block,offset,relx,obs,tree,direct,count\n"<<std::setprecision(14);for(int b=0;b<nblocks;++b)for(int k=0;k<5;++k)for(int o=0;o<2;++o){long long cnt=blocks[b].count[k][o];if(!cnt)continue;double rx=blocks[b].relcount[k]?blocks[b].relx[k]/blocks[b].relcount[k]:0;fo<<seed<<","<<b<<","<<offsets[k]<<","<<rx<<","<<(o==0?"j2l4":"j3l2")<<","<<blocks[b].tree[k][o]/cnt<<","<<blocks[b].direct[k][o]/cnt<<","<<cnt<<"\n";}fo.close();
 std::ofstream fm(prefix+"_meta.txt");fm<<std::setprecision(12)<<"Mach="<<Mach<<" seed="<<seed<<" dt="<<dt<<" burn="<<burntime<<" sample="<<samptime<<" pairsPer="<<pairsPer<<" nang="<<nang<<" nblocks="<<nblocks<<" etaHS="<<ETA_HS<<" runtime="<<std::chrono::duration<double>(t1c-t0).count()<<" Nfinal="<<ps.size()<<"\n";fm.close();
 std::cerr<<"done seed="<<seed<<" runtime="<<std::chrono::duration<double>(t1c-t0).count()<<" sec\n";
}
