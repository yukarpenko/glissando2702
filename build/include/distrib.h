/** \file distrib.h
 * Part of GLISSANDO 2
 * 
 */


#ifndef _GL_DISTRIB
  #define _GL_DISTRIB  


#include <TH1D.h>
#include <TH3D.h>

#include "counter.h"


//! Distribution of sources in space. 
/*!
Class for storage and basic operations (translation, rotation) of a distribution of "sources" in space. 
A source is a point with real weight and some additional integer property, e.g., charge.   
 */

class distr {
 public:
//! number of sources (points)
    Int_t    n;
//! x space coordinate 
    Float_t *x;
//! y space coordinate  
    Float_t *y;
//! z space coordinate    
    Float_t *z;
//! some integer property, here called "charge" 
/*! Depending on the situation, c is the electric charge of the nucleon in the nucleus, number of collisions of the nucleon, etc. */ 
    Int_t   *c; 
//! weight  
/*! Depending on the situation, the weight may describe the amount of the deposited energy, entropy, etc., in the source */ 
    Float_t *w;  
//! center-of-mass x coordinate of the distribution
    Float_t xcm, 
//! center-of-mass y coordinate of the distribution
    ycm, 
//! center-of-mass z coordinate of the distribution 
    zcm, 
//! mean squared radius of the distribution    
    msr, 
//! mean squared transverse radius of the distribution
    msrt, 
//! sum of weights of the distribution
    sumw;

//! constructor
    distr(
          Int_t k  /*!< number of sources, k>0 */  
         ){ 
           n=k;
           x=NULL;y=NULL;z=NULL;x=new Float_t[n];y=new Float_t[n];z=new Float_t[n];c=new Int_t[n];w=new Float_t[n];
	   for(Int_t i=0;i<n;i++){x[i]=0;y[i]=0;z[i]=0;c[i]=0;w[i]=1;};};

//! default constructor, 208 sources
/*! 208 corresponds to the number on nucleons in the 208Pb nucleus */
    distr(void){
           n=208;
           x=NULL;y=NULL;z=NULL;x=new Float_t[n];y=new Float_t[n];z=new Float_t[n];c=new Int_t[n];w=new Float_t[n];
	   for(Int_t i=0;i<n;i++){x[i]=0;y[i]=0;z[i]=0;c[i]=0;w[i]=1;};};

//! copying constructor
    distr(const distr & w1){ 
           n = w1.n;
           x=NULL;y=NULL;z=NULL;x=new Float_t[n];y=new Float_t[n];z=new Float_t[n];c=new Int_t[n];w=new Float_t[n];
	   for(Int_t i=0;i<n;i++){x[i]=w1.x[i];y[i]=w1.y[i];z[i]=w1.z[i];c[i]=w1.c[i];w[i]=w1.w[i];};};    

//! destructor
//   ~distr(){
// delete x; delete y; delete z; delete c; delete w;
//           }; 


//! generate sum of weights
    Float_t sum_w(){ 
           sumw=0;for(Int_t i=0;i<n;i++){sumw+=w[i];};return sumw;};

//! generate the center-of-mass x coordinate, no weights
    Float_t cmx(){ 
           xcm=0;for(Int_t i=0;i<n;i++){xcm+=x[i];};xcm/=n;return xcm;};


//! generate the center-of-mass y coordinate, no weights
    Float_t cmy(){
           ycm=0;for(Int_t i=0;i<n;i++){ycm+=y[i];};ycm/=n;return ycm;};

//! generate the center-of-mass z coordinate, no weights
    Float_t cmz(){
           zcm=0;for(Int_t i=0;i<n;i++){zcm+=z[i];};zcm/=n;return zcm;};

//! generate the center-of-mass x coordinate with weights
    Float_t cmx_w(){
           sum_w();xcm=0;for(Int_t i=0;i<n;i++){xcm+=x[i]*w[i];};xcm/=sumw;return xcm;};

//! generate the center-of-mass y coordinate with weights
    Float_t cmy_w(){
           sum_w();ycm=0;for(Int_t i=0;i<n;i++){ycm+=y[i]*w[i];};ycm/=sumw;return ycm;};

//! generate the center-of-mass z coordinate with weights
    Float_t cmz_w(){
           sum_w();zcm=0;for(Int_t i=0;i<n;i++){zcm+=z[i]*w[i];};zcm/=sumw;return zcm;};

//! translate in the x direction
    void shift_x(
                 Float_t xt /*!< displacement in the x direction */ 
                ){
           for(Int_t i=0;i<n;i++){x[i]+=xt;};};

//! translate in the y direction
    void shift_y(
                 Float_t yt /*!< displacement in the y direction */ 
                ){
           for(Int_t i=0;i<n;i++){y[i]+=yt;};};

//! translate in the z direction
    void shift_z(
                 Float_t zt /*!< displacement in the z direction */ 
                ){
           for(Int_t i=0;i<n;i++){z[i]+=zt;};};

//! translate to the cm reference frame in the x direction, no weights
    void shift_cmx(){ 
           cmx();for(Int_t i=0;i<n;i++){x[i]-=xcm;};};

//! translate to the cm reference frame in the y direction, no weights
    void shift_cmy(){
           cmy();for(Int_t i=0;i<n;i++){y[i]-=ycm;};};

//! translate to the cm reference frame in the z direction, no weights
    void shift_cmz(){
           cmz();for(Int_t i=0;i<n;i++){z[i]-=zcm;};};

//! translate to the cm reference frame in the x direction with weights
    void shift_cmx_w(){
           cmx_w();for(Int_t i=0;i<n;i++){x[i]-=xcm;};};

//! translate to the cm reference frame in the y direction with weights
    void shift_cmy_w(){
           cmy_w();for(Int_t i=0;i<n;i++){y[i]-=ycm;};};

//! translate to the cm reference frame in the z direction with weights
    void shift_cmz_w(){
           cmz_w();for(Int_t i=0;i<n;i++){z[i]-=zcm;};};

//! mean squared x, no weights
    Float_t msx(){ 
           cmx();
           msr=0;for(Int_t i=0;i<n;i++){msr+=(x[i]-xcm)*(x[i]-xcm);};msr/=n;return msr;};

//! mean squared y, no weights
    Float_t msy(){ 
           cmy();
           msr=0;for(Int_t i=0;i<n;i++){msr+=(y[i]-ycm)*(y[i]-ycm);};msr/=n;return msr;};

//! mean x*y, no weights
    Float_t mxy(){ 
           cmx();cmy();
           msr=0;for(Int_t i=0;i<n;i++){msr+=(x[i]-xcm)*(y[i]-ycm);};msr/=n;return msr;};

//! mean squared radius, no weights
    Float_t msrad(){ 
           cmx();cmy();cmz();
           msr=0;for(Int_t i=0;i<n;i++){msr+=(x[i]-xcm)*(x[i]-xcm)+(y[i]-ycm)*(y[i]-ycm)+(z[i]-zcm)*(z[i]-zcm);};msr/=n;return msr;};

//! mean squared transverse radius, no weights
    Float_t msrad_t(){
           cmx();cmy();
           msrt=0;for(Int_t i=0;i<n;i++){msrt+=(x[i]-xcm)*(x[i]-xcm)+(y[i]-ycm)*(y[i]-ycm);};msrt/=n;return msrt;};

//! mean squared radius with weights
    Float_t msrad_w(){
           cmx_w();cmy_w();cmz_w();
           msr=0;for(Int_t i=0;i<n;i++){msr+=w[i]*((x[i]-xcm)*(x[i]-xcm)+(y[i]-ycm)*(y[i]-ycm)+(z[i]-zcm)*(z[i]-zcm));};msr/=sumw;return msr;};

//! mean squared transverse radius with weights
    Float_t msrad_t_w(){
           cmx_w();cmy_w();
           msrt=0;for(Int_t i=0;i<n;i++){msrt+=w[i]*((x[i]-xcm)*(x[i]-xcm)+(y[i]-ycm)*(y[i]-ycm));};msrt/=sumw;return msrt;};

//! size - the weighted average of the distance from the orgin (in the cm frame)
    Float_t size(){
           cmx_w();cmy_w();
           Float_t siz=0;for(Int_t i=0;i<n;i++){siz+=w[i]*sqrt((x[i]-xcm)*(x[i]-xcm)+(y[i]-ycm)*(y[i]-ycm));};siz/=sumw;return siz;};

//! rotation angle maximizing the m-th cosine Fourier moment with weight r^2
/*! for m=2 this is the angle for passing to the "participant" frame */
    Float_t phrot(
                Int_t m /*!< rank of the Fourier moment, m=0,2,3,4,5,... (m=1 does not make sense in the cm frame) */
               ){
           Float_t sc=0; Float_t ss=0;
           for(Int_t i=0;i<n;i++){if(c[i]*c[i] > -1){Float_t r2=x[i]*x[i]+y[i]*y[i];Float_t ph=atan2(x[i],y[i]);
                                sc+=r2*cos(m*ph)*w[i]; ss+=r2*sin(m*ph)*w[i];};};
                                return atan2(ss,sc)/m;};

//! rotation angle maximizing the m-th cosine Fourier moment with weight r^k
/*! for m=2 this is the angle for passing to the "participant" frame */
    Float_t phrot(
                Int_t m,   /*!< rank of the Fourier moment, m=0,2,3,4,5,... (m=1 does not make sense in the cm frame) */
                Float_t k  /*!< power of transverse radius in the weight */
               ){
           Float_t sc=0; Float_t ss=0;
           for(Int_t i=0;i<n;i++){if(c[i]*c[i] > -1){Float_t r2=pow(x[i]*x[i]+y[i]*y[i],k/2.);Float_t ph=atan2(x[i],y[i]);
                                sc+=r2*cos(m*ph)*w[i]; ss+=r2*sin(m*ph)*w[i];};};
                                return atan2(ss,sc)/m;};

//! rotate in the transverse plane by the specified angle 
    void rotate(
                Float_t ph /*!< rotation angle in the transverse plane */
               ){
           for(Int_t i=0;i<n;i++){Float_t xp=x[i]*cos(ph)-y[i]*sin(ph);Float_t yp=y[i]*cos(ph)+x[i]*sin(ph);x[i]=xp;y[i]=yp;};};


//! rotate in the ZX plane by the theta angle
    void rotate_polar(
                      Float_t costh /*!< cosine of the rotation angle (theta) in the ZX plane */
               ){
            Float_t sinth=sqrt(1.-costh*costh);
            for(Int_t i=0;i<n;i++){Float_t zp=z[i]*costh-x[i]*sinth;
                                   Float_t xp=z[i]*sinth+x[i]*costh;
                                   z[i]=zp;x[i]=xp;};};

//! full event output to an external file	
    void writerds(
                  ofstream& eveout /*!< external file for the full event data */
                 ){
      eveout << n << endl;
      for(int i=0;i<n;i++){eveout << x[i] << " " << y[i] << " " << c[i] << " " << w[i] << endl;}
    };
	   
//! m-th cosine Fourier moment in the azimuthal angle in the transverse plane, weight r^2, limited charge range
/*! Limited charge range may be used to get the core and mantle (corona) distributions. */
    Float_t eps(
              Int_t m,    /*!< rank of the Fourier moment, m=0,2,3,4,5,... (m=1 does not make sense in the cm frame) */
              Int_t imin, /*!< lowest charge for the sources included */
              Int_t imax  /*!< highest charge for the sources included */
             ){
                      Float_t sc=0; Float_t ss=0;
                      for(Int_t i=0;i<n;i++){if((c[i]>=imin) && (c[i]<=imax)){Float_t r2=x[i]*x[i]+y[i]*y[i];Float_t ph=atan2(x[i],y[i]);
                                sc+=r2*cos(m*ph)*w[i]; ss+=r2*w[i];};};
                                return sc/ss;};

//! m-th cosine Fourier moment in the azimuthal angle in the transverse plane, weight r^2, no limit on charge
/*! Most popular is the second moment, known as eccentricity. The third moment is relevant for the triangular flow. */   
    Float_t eps(
              Int_t m /*!< rank of the Fourier moment, m=0,2,3,4,5,... (m=1 does not make sense in the cm frame) */
             ){ 
                      Float_t sc=0; Float_t ss=0;
                      for(Int_t i=0;i<n;i++){Float_t r2=x[i]*x[i]+y[i]*y[i];Float_t ph=atan2(x[i],y[i]);
                                sc+=r2*cos(m*ph)*w[i]; ss+=r2*w[i];};
                                return sc/ss;};

//! m-th sine Fourier moment in the azimuthal angle in the transverse plane, weight r^2, limited charge 
/*! Limited charge range may be used to get the core and mantle (corona) distributions. */
    Float_t eps_s(
                Int_t m,    /*!< rank of the Fourier moment, m=0,2,3,4,5,... (m=1 does not make sense in the cm frame) */
                Int_t imin, /*!< lowest charge for the sources included */
                Int_t imax  /*!< highest charge for the sources included */
               ){
                      Float_t sc=0; Float_t ss=0;
                      for(Int_t i=0;i<n;i++){if((c[i]>=imin) && (c[i]<=imax)){Float_t r2=x[i]*x[i]+y[i]*y[i];Float_t ph=atan2(x[i],y[i]);
                                sc+=r2*sin(m*ph)*w[i]; ss+=r2*w[i];};};
                                return sc/ss;};

//! m-th sine Fourier moment in the azimuthal angle in the transverse plane, weight r^2, no limit on charge   
    Float_t eps_s(
                Int_t m /*!< rank of the Fourier moment, m=2,3,4,5,... (m=1 does not make sense in the cm frame) */
               ){
                      Float_t sc=0; Float_t ss=0;
                      for(Int_t i=0;i<n;i++){Float_t r2=x[i]*x[i]+y[i]*y[i];Float_t ph=atan2(x[i],y[i]);
                                sc+=r2*sin(m*ph)*w[i]; ss+=r2*w[i];};
                                return sc/ss;};


//! m-th cosine Fourier moment in the azimuthal angle in the transverse plane, weight r^k, limited charge range
/*! Limited charge range may be used to get the core and mantle (corona) distributions. */
    Float_t eps(
              Int_t m,    /*!< rank of the Fourier moment, m=0,2,3,4,5,... (m=1 does not make sense in the cm frame) */
              Float_t k,  /*!< power of transverse radius in the weight */
              Int_t imin, /*!< lowest charge for the sources included */
              Int_t imax  /*!< highest charge for the sources included */
             ){
                      Float_t sc=0; Float_t ss=0;
                      for(Int_t i=0;i<n;i++){if((c[i]>=imin) && (c[i]<=imax)){Float_t r2=pow(x[i]*x[i]+y[i]*y[i],k/2.);Float_t ph=atan2(x[i],y[i]);
                                sc+=r2*cos(m*ph)*w[i]; ss+=r2*w[i];};};
                                return sc/ss;};

//! m-th cosine Fourier moment in the azimuthal angle in the transverse plane, weight r^k, no limit on charge
/*! Most popular is the second moment, known as eccentricity. The third moment is relevant for the triangular flow. */   
    Float_t eps(
              Int_t m,  /*!< rank of the Fourier moment, m=0,2,3,4,5,... (m=1 does not make sense in the cm frame) */
              Float_t k /*!< power of transverse radius in the weight */
             ){ 
                      Float_t sc=0; Float_t ss=0;
                      for(Int_t i=0;i<n;i++){Float_t r2=pow(x[i]*x[i]+y[i]*y[i],k/2.);Float_t ph=atan2(x[i],y[i]);
                                sc+=r2*cos(m*ph)*w[i]; ss+=r2*w[i];};
                                return sc/ss;};

//! m-th sine Fourier moment in the azimuthal angle in the transverse plane, weight r^k, limited charge 
/*! Limited charge range may be used to get the core and mantle (corona) distributions. */
    Float_t eps_s(
                Int_t m,    /*!< rank of the Fourier moment, m=0,2,3,4,5,... (m=1 does not make sense in the cm frame) */
                Float_t k,  /*!< power of transverse radius in the weight */
                Int_t imin, /*!< lowest charge for the sources included */
                Int_t imax  /*!< highest charge for the sources included */
               ){
                      Float_t sc=0; Float_t ss=0;
                      for(Int_t i=0;i<n;i++){if((c[i]>=imin) && (c[i]<=imax)){Float_t r2=pow(x[i]*x[i]+y[i]*y[i],k/2.);Float_t ph=atan2(x[i],y[i]);
                                sc+=r2*sin(m*ph)*w[i]; ss+=r2*w[i];};};
                                return sc/ss;};

//! m-th sine Fourier moment in the azimuthal angle in the transverse plane, weight r^k, no limit on charge   
    Float_t eps_s(
                Int_t m,  /*!< rank of the Fourier moment, m=2,3,4,5,... (m=1 does not make sense in the cm frame) */
                Float_t k /*!< power of transverse radius in the weight */
               ){
                      Float_t sc=0; Float_t ss=0;
                      for(Int_t i=0;i<n;i++){Float_t r2=pow(x[i]*x[i]+y[i]*y[i],k/2.);Float_t ph=atan2(x[i],y[i]);
                                sc+=r2*sin(m*ph)*w[i]; ss+=r2*w[i];};
                                return sc/ss;};
  
  

//! fill the histogram of the transverse distribution in cartesian coordinates, limited charge  
/*! Limited charge range may be used to get the core and mantle (corona) distributions. */
    void fill_xy(
                TH2D *xyh, /*!< 2-dim cartesian histogram in the transverse plane */ 
                Float_t fac, /*!< normalization factor */
                Int_t imin,  /*!< lowest charge for the sources included */
                Int_t imax   /*!< highest charge for the sources included */
                ){
           for(Int_t i=0;i<n;i++){if((c[i]*c[i]>=imin*imin) && (c[i]*c[i]<=imax*imax)){xyh->Fill(x[i],y[i],w[i]*fac);};};};

//! fill the histogram of the transverse distribution in cartesian coordinates, no limit on charge  
    void fill_xy(
                TH2D *xyh,  /*!< 2-dim cartesian histogram in the transverse plane */ 
                Float_t fac   /*!< normalization factor */
                ){
           for(Int_t i=0;i<n;i++){xyh->Fill(x[i],y[i],w[i]*fac);};};

//! fill the histogram of the transverse distribution in polar coordinates, limited charge  
/*! Limited charge range may be used to get the core and mantle (corona) distributions. */
    void fill_polar(
                   TH2D *polh,  /*!< 2-dim polar histogram in the transverse plane */ 
                   Int_t m,       /*!< rank of the Fourier moment */
                   Float_t fac,   /*!< normalization factor */
                   Int_t imin,    /*!< lowest charge for the sources included */
                   Int_t imax     /*!< highest charge for the sources included */
                   ){
           for(Int_t i=0;i<n;i++){if((c[i]>=imin) && (c[i]<=imax)){Float_t r=sqrt(x[i]*x[i]+y[i]*y[i]); Float_t ph=atan2(x[i],y[i]);
                                if(r>0){polh->Fill(r,ph,fac*w[i]*cos(m*ph)/r);};};};};

//! fill the histogram of the transverse sine distribution in polar coordinates, limited charge  
/*! Limited charge range may be used to get the core and mantle (corona) distributions. */
    void fill_polar_s(
                   TH2D *polh,  /*!< 2-dim polar sine histogram in the transverse plane */ 
                   Int_t m,       /*!< rank of the sine Fourier moment */
                   Float_t fac,   /*!< normalization factor */
                   Int_t imin,    /*!< lowest charge for the sources included */
                   Int_t imax     /*!< highest charge for the sources included */
                   ){
           for(Int_t i=0;i<n;i++){if((c[i]*c[i]>=imin*imin) && (c[i]*c[i]<=imax*imax)){Float_t r=sqrt(x[i]*x[i]+y[i]*y[i]); Float_t ph=atan2(x[i],y[i]);
                                if(r>0){polh->Fill(r,ph,fac*w[i]*sin(m*ph)/r);};};};};


//! fill the histogram of the transverse distribution in polar coordinates, no limit on charge
    void fill_polar(
                   TH2D *polh,   /*!< 2-dim polar histogram in the transverse plane */ 
                   Int_t m,        /*!< rank of the Fourier moment */
                   Float_t fac     /*!< normalization factor */
                   ){
           for(Int_t i=0;i<n;i++){Float_t r=sqrt(x[i]*x[i]+y[i]*y[i]); Float_t ph=atan2(x[i],y[i]);
                                if(r>0){polh->Fill(r,ph,fac*w[i]*cos(m*ph)/r);};};};

//! fill the histogram of the transverse sine distribution in polar coordinates, no limit on charge
    void fill_polar_s(
                   TH2D *polh,   /*!< 2-dim polar sine histogram in the transverse plane */ 
                   Int_t m,        /*!< rank of the sine Fourier moment */
                   Float_t fac     /*!< normalization factor */
                   ){
           for(Int_t i=0;i<n;i++){Float_t r=sqrt(x[i]*x[i]+y[i]*y[i]); Float_t ph=atan2(x[i],y[i]);
                                if(r>0){polh->Fill(r,ph,fac*w[i]*cos(m*ph)/r);};};};


//! substitution overloading
    distr& operator=(const distr& w1); 
};

//! substitution overloading
distr& distr::operator=(const distr &w1) 
{	
  if(this != &w1){                           
           n = w1.n;
           x=NULL;y=NULL;z=NULL;x=new Float_t[n];y=new Float_t[n];z=new Float_t[n];c=new Int_t[n];w=new Float_t[n];
	   for(Int_t i=0;i<n;i++){x[i]=w1.x[i];y[i]=w1.y[i];z[i]=w1.z[i];c[i]=w1.c[i];w[i]=w1.w[i];};
	   return *this;};	        
};


//! nucleus class
/*! Class to store distributions of nucleons in nuclei. */
class nucleus : public distr {
  
  private:
    bool g;
    Float_t cth, sth, phi, r;

  public:
//! constructor  
  nucleus(
          Int_t k /*!< number of nucleons in the nucleus (the mass naumber of the nucleus), 
                     k>0, k=1 - proton, k=2 - deuteron, k>2 - other nucleus */
         ) : distr(k){};

//! copying constructor 
  nucleus(const nucleus& w1) : distr(w1){};

//! substitution overloading
  nucleus& operator = (const nucleus& w1)
         {distr::operator=(w1);return *this;};

//! destructor
//  ~nucleus(){}; 

//! set the distribution for the proton 
/*! The proton is just placed at the origin */
  void set_proton(){x[0]=0.; y[0]=0.; z[0]=0.;}; 

//! set the distribution for the deuteron 
/*! Use the Hulthen distribution. */
  void set_deuteron(){x[0]=0.; y[0]=0.; z[0]=0.;
       phi=2*PI*los(); cth=2*los()-1; sth=sqrt(1.-cth*cth); r=rlos_hult(); 
       x[1]=r*cos(phi)*sth; y[1]=r*sin(phi)*sth; z[1]=r*cth;}; 

//! set randomly the distribution of nucleons in the nucleus A, use the rlosA() function, no correlations
  void set_random_A(){for(Int_t i=0;i<n;i++){phi=2*PI*los(); cth=2*los()-1; sth=sqrt(1.-cth*cth); r=rlosA(); 
       x[i]=r*cos(phi)*sth; y[i]=r*sin(phi)*sth; z[i]=r*cth;};}; 

//! set randomly the distribution of nucleons in the nucleus A, use the rlosA_hos() function, no correlations
void set_random_A_hos(){for(Int_t i=0;i<n;i++){phi=2*PI*los(); cth=2*los()-1; sth=sqrt(1.-cth*cth); r=rlosA_hos(); 
     x[i]=r*cos(phi)*sth; y[i]=r*sin(phi)*sth; z[i]=r*cth;};}; 	   

//! set randomly the distribution of nucleons in the nucleus B, use the rlosB() function, no correlations 
  void set_random_B(){for(Int_t i=0;i<n;i++){phi=2*PI*los(); cth=2*los()-1; sth=sqrt(1.-cth*cth); r=rlosB(); 
       x[i]=r*cos(phi)*sth; y[i]=r*sin(phi)*sth; z[i]=r*cth;};}; 

//! set randomly the distribution of nucleons in the nucleus B, use the rlosB_hos() function, no correlations 
  void set_random_B_hos(){for(Int_t i=0;i<n;i++){phi=2*PI*los(); cth=2*los()-1; sth=sqrt(1.-cth*cth); r=rlosB_hos(); 
       x[i]=r*cos(phi)*sth; y[i]=r*sin(phi)*sth; z[i]=r*cth;};}; 

//! set randomly the distribution of nucleons in the nucleus A with deformation, no correlations
  void set_random_A_def(){for(Int_t i=0;i<n;i++){phi=2*PI*los();   
                                                 r=rlosA_def(&cth,BETA2A,BETA4A); 
                                                 sth=sqrt(1.-cth*cth);
       x[i]=r*cos(phi)*sth; y[i]=r*sin(phi)*sth; z[i]=r*cth;};}; 

//! set randomly the distribution of nucleons in the nucleus B with deformation, no correlations 
  void set_random_B_def(){for(Int_t i=0;i<n;i++){phi=2*PI*los();   
                                                 r=rlosB_def(&cth,BETA2B,BETA4B); 
                                                 sth=sqrt(1.-cth*cth);
       x[i]=r*cos(phi)*sth; y[i]=r*sin(phi)*sth; z[i]=r*cth;};}; 

  void set_random_A_def(
                       Float_t d /*!< expulsion distance, nucleons cannot be closer to each other than d */
                       ){Int_t j=0;while(j<n){phi=2*PI*los();
                                              r=rlosA_def(&cth,BETA2A,BETA4A);
                                              sth=sqrt(1.-cth*cth);
       x[j]=r*cos(phi)*sth; y[j]=r*sin(phi)*sth; z[j]=r*cth; if(good_down(j,d)){
 j++;
}; };};


  void set_random_B_def(
                       Float_t d /*!< expulsion distance, nucleons cannot be closer to each other than d */
                       ){Int_t j=0;while(j<n){phi=2*PI*los();
                                              r=rlosA_def(&cth,BETA2B,BETA4B);
                                              sth=sqrt(1.-cth*cth);
       x[j]=r*cos(phi)*sth; y[j]=r*sin(phi)*sth; z[j]=r*cth; if(good_down(j,d)){j++;}; };};



//! set randomly the distribution of nucleons in the nucleus A with expulsion, use the rlosA() function
/*! The nucleon positions are subsequntly generated acording to the spherically-symmetric 
    Woods-Saxon distribution. If the nucleon happens to be generated closer than the expulsion distance d 
    to any of the previously generated nucleons, it is generated anew, until it is "good". 
    Since this leads to some swelling, the original distribution must be a bit narrower to 
    cancel neutrilize this effect (see our original paper for a detailed discussion). 
    The expulsion simulates in a simple manner the nuclear repulsion and generates the hard-core two-body 
    correlations. */
  void set_random_A(
                   Float_t d /*!< expulsion distance, nucleons cannot be closer to each other than d */
                   ){Int_t j=0;while(j<n){phi=2*PI*los(); cth=2*los()-1; sth=sqrt(1.-cth*cth); r=rlosA(); 
       x[j]=r*cos(phi)*sth; y[j]=r*sin(phi)*sth; z[j]=r*cth; if(good_down(j,d)){j++;};};};

//! similar as the above function but for harmonic oscillator shell model, use use the rlosA_hos() function.
  void set_random_A_hos(
                   Float_t d /*!< expulsion distance, nucleons cannot be closer to each other than d */
                   ){Int_t j=0;while(j<n){phi=2*PI*los(); cth=2*los()-1; sth=sqrt(1.-cth*cth); r=rlosA_hos(); 
       x[j]=r*cos(phi)*sth; y[j]=r*sin(phi)*sth; z[j]=r*cth; if(good_down(j,d)){j++;};};};

//! set randomly the distribution of nucleons in the nucleus B with expulsion, use the rlosB() function
/*! Same as set_random_A for the case of the nucleus B, which in general is different from A */
  void set_random_B(
                   Float_t d /*!< expulsion distance, nucleons cannot be closer to each other than d */
                   ){Int_t j=0;while(j<n){phi=2*PI*los(); cth=2*los()-1; sth=sqrt(1.-cth*cth); r=rlosB(); 
       x[j]=r*cos(phi)*sth; y[j]=r*sin(phi)*sth; z[j]=r*cth; if(good_down(j,d)){j++;};};};

//! similar as above function but for harmonic oscillator shell model, use use the rlosB_hos() function.
  void set_random_B_hos(
                   Float_t d /*!< expulsion distance, nucleons cannot be closer to each other than d */
                   ){Int_t j=0;while(j<n){phi=2*PI*los(); cth=2*los()-1; sth=sqrt(1.-cth*cth); r=rlosB_hos(); 
       x[j]=r*cos(phi)*sth; y[j]=r*sin(phi)*sth; z[j]=r*cth; if(good_down(j,d)){j++;};};};

//! set the distribution of nucleons in the nucleus from the tables generated earlier
/*! The nucleon position (x,y,z) is taken from the tables read earlier. The pointers x, y, z are 
    originally positioned at 
    px, py, pz at a place corresponding to the beginning of a randomly selected nucleus. 
    The tables with nuclear distributions have to be prepared externally, 
    see, e.g., http://www.phys.psu.edu/~malvioli/eventgenerator/ 
    for distributions involving correlations.
    */
  void set_file(
                Float_t *px, /*!< x coordinate of distributions read from files */ 
                Float_t *py, /*!< y coordinate of distributions read from files */ 
                Float_t *pz, /*!< z coordinate of distributions read from files */ 
                Int_t sn,    /*!< number of entries in the file (should be the mass number time the number of stored nuclei) */
                Int_t nu     /*!< mass number of the nucleus */ 
               ){  
                   Int_t count=raa.Integer(sn/nu)*nu;
                   x = (px+count); y = (py+count); z = (pz+count); // set the pointers
   }; 

//! set the distribution of nucleons in the nucleus from the tables generated earlier, killing correlations
/*! The nucleon position (x,y,z) is taken from the tables read earlier. The pointers x, y, z are positioned at 
    px, py, pz at a place corresponding to a completely randomly selected nucleon. This procedure kills any 
    correlations and is (probably) equivalent to the mixing technique.  
    The tables with nuclear distributions have to be prepared externally.
    */
  void set_file_uncor(
                Float_t *px, /*!< x coordinate of distributions read from files */ 
                Float_t *py, /*!< y coordinate of distributions read from files */ 
                Float_t *pz, /*!< z coordinate of distributions read from files */ 
                Int_t sn,    /*!< number of entries in the file (should be the mass number time the number of stored nuclei) */
                Int_t nu     /*!< mass number of the nucleus */ 
               ){  
                   Int_t count=raa.Integer(sn-nu);
                   x = (px+count); y = (py+count); z = (pz+count); // set the pointers
   }; 

//! distance between two nucleons 
  Float_t dist2(
             Int_t j1, /*!< index of nucleon 1 */
             Int_t j2  /*!< index of nucleon 2 */
             ){return (x[j1]-x[j2])*(x[j1]-x[j2])+(y[j1]-y[j2])*(y[j1]-y[j2])+(z[j1]-z[j2])*(z[j1]-z[j2]);};


//! the pair of nucleons is "good" when the distance between the nucleons is larger than the expulsion distance d
  bool good_pair(
                Int_t j1,  /*!< index of nucleon 1 */
                Int_t j2,  /*!< index of nucleon 2 */
                Float_t d  /*!< expulsion distance - nucleons cannot be closer to each other than d*/
                ){if(dist2(j1,j2)>d*d){return true;} else {return false;};};
  
//! nucleon j is "good" when the distance to all nucleons of index i with i<j is larger than the expulsion distance d
  bool good_down(
                Int_t j,   /*!< index of the nucleon */
                Float_t d  /*!< expulsion distance */
                ){if(j==0){return true;} else {g=true; for(Int_t i=0;i<j;i++){g=g && good_pair(i,j,d);};return g;};};

//! configuration of nucleons in the nucleus is "good" when the distances between all nucleons are larger than the expulsion distance d
  bool good_all(
               Float_t d   /*!< expulsion distance */
               ){g=true; for(Int_t i=1;i<n;i++){g=g && good_down(i,d);};return g;};

}; //class

#endif
