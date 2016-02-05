/** \file collision.h
 * Part of GLISSANDO 2
 * 
 */


#ifndef _GL_COLLISION
  #define _GL_COLLISION  

#include "distrib.h"
#include <TMath.h>

//! collision class
/*! 
Class performing the collision of two nuclei. Here the "source" is the position of the wounded nucleon (a nucleon that 
collided at least once) or the location of the binary collision, taken as the center-of-mass position in the nucleon pair. 
The weight, called RDS (relative deposited strength) depends on the adopted model. The "charge" is 0 for the binary 
collision, and equal to i for the wounded nucleons, where i>0 is the number of collisions experienced by the nucleon. 
Depending on the precompiler option _nnwp_, the wounding and binary collisions occur with a step-function distribution (_nnwp_=0),
Gaussian distribution (_nnwp_=1) or gamma approximation (_nnwp_=2). The latter is much more realistic.
*/

class collision : public distr {

  public:

  Float_t rd2; /*!< square of the distance between the nucleons */
  Int_t wc;    /*!< counter of the sources */
  Int_t *wwA,  /*!< wwA[i] is the number of collisions of the i-th nucleon from nucleus A with the nucleons from nucleus B */ 
        *wwB;  /*!< wwB[i] is the number of collisions of the i-th nucleon from nucleus A with the nucleons from nucleus A */
  Int_t nwA,   /*!< number of wounded (collided at least once) nucleons in nucleus A */
        nwB,   /*!< number of wounded (collided at least once) nucleons in nucleus B */
        nwAB,  /*!< total number of wounded nucleons (nwA+nwB) generated in the event */
        nzw,   /*!< number of sources with non-zero weight */
        nbin,  /*!< number of binary collisions in the event */
        nhotspot; /*!< number of hot spots */
  Float_t rpa; /*!< weight of the source (RDS) */
  
  #if(_weight_)
  TH1D *w_distr,      /*!< histogram for the collision profile for the wounded nucleons */
       *w_distr_bin;  /*!< histogram for the collision profile for the binary collisions */
  #endif

//! constructor
  collision(
           Int_t nnA,  /*!< mass number of nucleus A */ 
           Int_t nnB   /*!< mass number of nucleus B */ 
           ){
            Int_t nn=nnA+nnB+nnA*nnB;
  #if(_weight_)
           w_distr = new TH1D("w_distr", "wounding profile p(b^{2})", 100, 0.000001, 15.);
           w_distr->SetXTitle("NN distance squared [fm^{2}]");
           w_distr_bin = new TH1D("w_distr_bin", "binary profile  p(b^{2})", 100, 0.000001, 15.);
           w_distr_bin->SetXTitle("NN distance squared [fm^{2}]");
  #endif
           x=NULL;y=NULL;x=new Float_t[nn];y=new Float_t[nn];c=new Int_t[nn];w=new Float_t[nn];
           wwA=NULL;wwB=NULL;wwA=new Int_t[nnA];wwB=new Int_t[nnB];};


//! default constructor
/*! The default constructor assumes 208Pb-208Pb collisions. */
  collision(){
  #if(_weight_)
           w_distr = new TH1D("w_distr", "wounding profile  p(b^{2})", 100, 0.000001, 15.);
           w_distr->SetXTitle("NN distance squared [fm^{2}]");
           w_distr_bin = new TH1D("w_distr_bin", "binary profile  p(b^{2})", 100, 0.000001, 15.);
           w_distr_bin->SetXTitle("NN distance squared [fm^{2}]");
  #endif
           x=NULL;y=NULL;x=new Float_t[43680];y=new Float_t[43680];c=new Int_t[43680];w=new Float_t[43680];
           wwA=NULL;wwB=NULL;wwA=new Int_t[208];wwB=new Int_t[208];};

//! copying constructor 
  collision(const collision& w1) : distr(w1){};

//! substitution overloading
  collision& operator = (const collision& w1)
         {distr::operator=(w1);return *this;};

//! destructor
// ~collision(){
// delete x; delete y; delete z; delete c; delete w;
//              }; 

//! collision between two nuclei
/*! gen_RDS performs the collision of two nuclei, generating the wounded nucleon and the binary collisions, as well as
their RDS (relative deposited strength, or weight). It is the core function of the code, 
implementing the specific model mechanism of the collision. */
  void gen_RDS(
            const nucleus &nA,    /*!< nucleus A, nA>0, nA=1 - proton, nA=2 - deuteron, nA>2 - other nuclei */
            const nucleus &nB,    /*!< nucleus B, nB>0, nB=1 - proton, nB=2 - deuteron, nB>2 - other nuclei */
            Float_t d2,           /*!< wounding distance squared */
            Float_t dbin2,        /*!< binary-collison distance squared */
            Float_t mb            /*!< ratio of the wounding to binary cross sections */
            ){

	Float_t a     = 1./OMEGA;
    Float_t z;
	
   for(Int_t i=0;i<nA.n;i++){wwA[i]=0;};
   for(Int_t j=0;j<nB.n;j++){wwB[j]=0;};
   nwA=0; nwB=0; nzw=0; rpa=0; nbin=0; nhotspot=0; 
   wc=-1;

   for(Int_t i=0;i<nA.n;i++){
	for(Int_t j=0;j<nB.n;j++){
		rd2=(nA.x[i]-nB.x[j])*(nA.x[i]-nB.x[j])+(nA.y[i]-nB.y[j])*(nA.y[i]-nB.y[j]);

			#if(_nnwp_==0)
        		   if(rd2 < d2)
			#elif(_nnwp_==1)
			   if(los() < GA*exp(-GA*rd2/(d2)))
			#elif(_nnwp_==2)
			   z = GAMA*rd2/(d2*OMEGA);
			   if(los() < GAMA*(1.0 - TMath::Gamma(a,z)))   
			#endif

                {wwA[i]++;wwB[j]++;
                 #if(_weight_) 
                 w_distr->Fill(rd2,1);
                 #endif
    };};};

Float_t aux;

   for(Int_t i=0;i<nA.n;i++){if(wwA[i]>0) 
	       	       {
                        nwA++; // count the wounded nucleons in A
			wc++; // count the sources
                        c[wc]=wwA[i];
 			x[wc]=nA.x[i]; // store x-coordinate of the source
			if(DW>0){x[wc]+=disp(DW);}; // displace randomly
 			y[wc]=nA.y[i]; // store y-coordinate of the source
			if(DW>0){y[wc]+=disp(DW);}; // displace randomly
			w[wc]=(1-ALPHA)/2.0*dist(MODEL,Uw,Vw); // wounded nucleon RDS
           		rpa+=w[wc];
	        	if(w[wc]>0){nzw++;}; // count wounded if RDS > 0
                        };};

   for(Int_t i=0;i<nB.n;i++){if(wwB[i]>0) 
	       	       {nwB++; // count the wounded nucleons in B
			wc++; // count the sources
                        c[wc]=-wwB[i]; // by convention, negative numbers for nucleus B
 			x[wc]=nB.x[i]; // store x-coordinate of the source
			if(DW>0){x[wc]+=disp(DW);}; // displace randomly
 			y[wc]=nB.y[i]; // store y-coordinate of the source
			if(DW>0){y[wc]+=disp(DW);}; // displace randomly
			w[wc]=(1-ALPHA)/2.0*dist(MODEL,Uw,Vw); // wounded nucleon RDS
           		rpa+=w[wc];
	        	if(w[wc]>0){nzw++;}; // count wounded in nzw if RDS > 0
                       };};

if(ALPHA>0 || DOBIN==1){
   for(Int_t i=0;i<nA.n;i++){
	for(Int_t j=0;j<nB.n;j++){
		rd2=(nA.x[i]-nB.x[j])*(nA.x[i]-nB.x[j])+(nA.y[i]-nB.y[j])*(nA.y[i]-nB.y[j]);

			#if(_nnwp_==0)
        		   if(rd2 < d2)
			#elif(_nnwp_==1)
			   if(los() < GA*exp(-GA*rd2/(d2)))
			#elif(_nnwp_==2)
			   z = GAMA*rd2/(d2*OMEGA);
			   if(los() < GAMA*(1.0 - TMath::Gamma(a,z)))   
			#endif

			       {
                        #if(_weight_)
                                w_distr_bin->Fill(rd2,1);
                        #endif
                                nbin++; // count the binary collisions
                                wc++;   // count the sources
                                c[wc]=0; // be convention, 0 for binary collisions
				x[wc]=(nA.x[i]+nB.x[j])/2;
				if(DBIN>0){// displace randomly
					x[wc]+=disp(DBIN);
				};
                        	y[wc]=(nA.y[i]+nB.y[j])/2;
				if(DBIN>0){// displace randomly
					y[wc]+=disp(DBIN);
				};
		  		if(los()>mb){
					w[wc]=0;
				} 
				else {
					w[wc]=ALPHA/mb*dist(MODEL,Ubin,Vbin); // binary RDS
					nhotspot++; 
					rpa+=w[wc]; 
					if(w[wc]>0){nzw++;};  // count binary in nzw if RDS > 0
				};
			   };
		};
	};
};

  n=wc+1; nwAB=nwA+nwB;
};

}; //class


//! collision_rap class
/*! Collision class with an extra feature of generating the rapidity distribution. It can be used to create fully 3-dimensional 
distribution of sources, overlaying the (spatial) rapidity distribution over the transverse distribution. */
class collision_rap : public collision {

public:

TH3D *rap_distr;      /*!< histogram for storing the rapidity distribution */
 
//! constructor  
  collision_rap(
           Int_t nnA,  /*!< mass number of nucleus A */ 
           Int_t nnB   /*!< mass number of nucleus B */ 
           ) : collision(nnA,nnB)
          {
       rap_distr = new TH3D("rap_distr", "3 dim. x-y-rapidity distribution of particles", NBIN, -BTOT, BTOT,NBIN,-BTOT,BTOT,NBIN,-RAPRANGE,RAPRANGE);
           rap_distr->SetXTitle("x [fm]");  rap_distr->SetYTitle("y [fm]");  rap_distr->SetZTitle("#eta");  
         };

//! default constructor
  collision_rap() : collision()
          {
       rap_distr = new TH3D("rap_distr", "3 dim. x-y-rapidity distribution of particles", NBIN, -BTOT, BTOT,NBIN,-BTOT,BTOT,NBIN,-RAPRANGE,RAPRANGE);          
           rap_distr->SetXTitle("x [fm]");  rap_distr->SetYTitle("y [fm]");  rap_distr->SetZTitle("#eta");  
          };

//! copying constructor 
  collision_rap(const collision_rap& w1) : collision(w1){};

//! substitution overloading
  collision_rap& operator = (const collision_rap& w1)
         {collision::operator=(w1);return *this;};

//! destructor
   ~collision_rap(){
 delete x; delete y; delete z; delete c; delete w;
                   };


//! generate the rapidity distribution in a specified gap for the y spatial coordinate
void gen_rap(
            Int_t part,  /*!< number of particles per unit RDS */
            Float_t maxy /*!< maximum absolute value of the transverse y coordinate for histogramming in the x-y-rapidity histogram */
            ){
Float_t eta; // rapidity
for(Int_t i=0;i<n;i++) // loop over all n sources
  {
    Int_t np=w[i]*part;
    for(Int_t j=0;j<np;j++)
       {
         if(c[i]>0)
           eta=los_rap_A();
         else if(c[i]<0)          
           eta=los_rap_B();
         else
           eta=los_rap_bin();
  
         if(y[i]*y[i]< maxy*maxy){rap_distr->Fill(x[i],y[i],eta,w[i]);};       
       };
  };
   
             };
//! get the rapidity distribution at shifted rapidity in the Bialas-Czyz-Bozek model
void shift_rap(
          Float_t rr /*!< amount of the shift in rapidity */
              ){
                 for(Int_t i=0;i<n;i++) // loop over all n sources
                    {if(c[i]>0)
                        w[i]*=2*fg(rr)*fpm(rr); // modify the weight by the profile function for wounded A
                     else if(c[i]<0)
                        w[i]*=2*fg(rr)*fpm(-rr);// modify the weight by the profile function for wounded B
                     else
                        w[i]*=fg(rr);           // modify the weight by the profile function for binary
                     };
                };

}; //class

#endif
