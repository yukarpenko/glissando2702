/** \file functions2.h
 * Part of GLISSANDO 2
 * 
 */


#ifndef _GL_FUNCTION1
  #define _GL_FUNCTION1


#include <math.h>
#include <time.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>

#include "counter.h"

using namespace std; 


//! structure for output of the full event - transverse coordinates, weight, number of the event
typedef struct {
                Float_t X, //!< x coordinate
                        Y, //!< y coordinate
                        W; //!< z coordinate
                UInt_t KK; //!< number of the event
                } SOURCE; 

static SOURCE tSource;      //!< structure used when FULL=1 


/**********************************
parameters and their default values
**********************************/
Int_t      EVENTS=50000,    //!< number of generated events 
           NBIN=40,         //!< number of bins for histogramming in x, y, and r 
           FBIN=72,         //!< number of bins for histogramming in the azimuthal angle
           NUMA=208,        //!< mass number of nucleus A 
           NUMB=208,        //!< mass number of nucleus B 
           WMIN=2,          //!< minimum number of wounded nucleons to record the event 
           MODEL=0,         //!< switch for the superimposed multiplicity distribution: 0 - uniform, 1 - Poisson, 2 - Gamma, 3 - Negative Binomial
           DOBIN=0,         //!< 1 - count binary collisions even in the pure wounded-nucleon model. 0 - do not
           W0=2,            //!< minimum allowed number of wounded nucleons in the acceptance window
           W1=1000,         //!< maximum allowed number of wounded nucleons in the acceptance window
           SHIFT=1,         //!< 1 - shift the coordinates of the fireball to c.m. in the fixed-axes case (preferred), 0 - do not
           RET=0,           //!< 0 - fix-last algorithm (preferred), 1 - return-to-beginning algorithm for the generation of the nuclear distribution
           FULL=0,          //!< 1 - generate the full event tree (large output file), 0 - do not
           FILES=0,         //!< 1 - read distribution from files, 0 - do not
           NNWP=0,          //!< 0 - hard-sphere NN wounding profile, 1 - Gaussian NN wounding profile, 2 - gamma NN wounding profile
           NUMRAP=10,       //!< number of particles per unit weight generated in the whole rapidity range
           ARANK=2,         //!< rank of the Fourier moment for the forward-backward analysis
           PP=-1,           //!< power of the transverse radius in the Fourier moments 
           RO=0;            //!< rank of the rotation axes (0 - rotation rank = rank of the Fourier moment)
           
UInt_t     ISEED,        //!< read seed for the ROOT random number generator, if 0 - random seed generated
           ISEED1;       //!< copy of ISEED

Float_t      BMIN=0.,              //!< minimum value of the impact parameter in the acceptance window
             BMAX=25.,             //!< maximum value of the impact parameter in the acceptance window
             RDS0=0.,              //!< minimum value of the relative deposited strength (RDS) in the acceptance window 
             RDS1=100000,          //!< maximum value of the relative deposited strength (RDS) in the acceptance window 
             RWSA=6.407,            //!< Woods-Saxon radius for nucleus A
             AWSA=0.459,            //!< Woods-Saxon width for nucleus A
	     BETA2A=0.0,           //!< Deformation parameter beta2 of deformed Woods-Saxon distribution for nucleus A 
             BETA4A=0.0,           //!< Deformation parameter beta4 of deformed Woods-Saxon distribution for nucleus A  
             ROTA_THETA=-1.0,      //!< Parameter controlling the rotation of nucleus A in XZ plane (polar angle THETA, -1 means random rotation)
             ROTA_PHI=-1.0,        //!< Parameter controlling the rotation of nucleus A in XY plane (azimuthal angle PHI, -1 means random rotation)
             RWSB=6.407,            //!< Woods-Saxon radius for nucleus B
             AWSB=0.459,            //!< Woods-Saxon width for nucleus B
	     BETA2B=0.0,           //!< Deformation parameter beta2 of deformed Woods-Saxon distribution for nucleus B 
             BETA4B=0.0,           //!< Deformation parameter beta4 of deformed Woods-Saxon distribution for nucleus B 
             ROTB_THETA=-1.0,      //!< Parameter controlling the rotation of nucleus B in XZ plane (polar angle THETA, -1 means random rotation)
             ROTB_PHI=-1.0,        //!< Parameter controlling the rotation of nucleus B in XY plane (azimuthal angle PHI, -1 means random rotation)

             BTOT = fmax(RWSA,RWSB)+AWSA+AWSB,    //!< maximum coordinate value for some histograms

             WFA=0.,               //!< the w parameter for the Fermi distribution for nucleus A
             WFB=0.,               //!< the w parameter for the Fermi distribution for nucleus B
             SNN=73.5,              //!< NN "wounding" cross section in milibarns
             SBIN=73.5,             //!< NN binary cross section in milibarns
             ALPHA=0.15,            //!< the mixing parameter: 0 - wounded, 1 - binary, 0.145 - mixed (PHOBOS)
             Uw=2.,                //!< Poisson or Gamma parameters for superimposed distribution, wounded nucleons
             Ubin=2.,              //!< Poisson or Gamma parameters for superimposed distribution, binary collisions
             Vw=4.,                //!< Negative binomial variance, wounded nucleons
             Vbin=4.,              //!< Negative binomial variance, binary collisions
             PI=4.*atan(1.),       //!< the number pi
             CD=0.9,               //!< closest allowed distance (expulsion distance) between nucleons in the nucleus in fm (simulation of repulsion)
             DW=0.,                //!< dispersion of the location of the source for wounded nucleons (in fm)
             DBIN=0.,              //!< dispersion of the location of the source for binary collisions (in fm)
             GA=0.92,              //!< Gaussian wounding profile parameter (height at the origin)
             RAPRANGE=5.,          //!< range in rapidity
             ETA0=1.,		   //!< 2*ETA0 is the width of the plateau in eta
             ETAM=3.36,            //!< parameter of the Bialas-Czyz-Bozek model
             SIGETA=1.3,           //!< parameter controlling the width of the rapidity distribution
             MAXYRAP=10.,          //!< maximum absolute value of the y coordinate in the x-y-rapidity histogram 
             FBRAP=2.5,            //!< forward rapidity for the forward-backward analysis (backward rapidity = - FBRAP) 
             RCHA=5.66,            //!< harmonic oscillator shell model density mean squared charge radii of  12C-nucleus 
	     RCHB=5.66,            //!< harmonic oscillator shell model density mean squared charge radii of  12C-nucleus 
             RCHP=0.7714,          //!< harmonic oscillator shell model density mean squared charge radii of  proton
	     OMEGA=0.4,            //!< relative variance of cross-section fluctuations
             GAMA=1.0;             //!< gamma wounding profile parameter (height at the origin)
//! process the input file with parameters
void readpar(
            TString inpfile //!< name of the input file
            ){ 

FILE *in; 
in = fopen(inpfile,"rt");

char t[60];

char s1[]={"EVENTS"},
     s2[]={"NBIN"},
     s3[]={"NUMA"},
     s4[]={"NUMB"},
     s5[]={"WMIN"},
     s6[]={"BMIN"},
     s7[]={"BMAX"},
     s8[]={"RWSA"},
     s9[]={"AWSA"},
     s10[]={"BETA2A"},  
     s11[]={"BETA4A"},
     s12[]={"ROTA_THETA"},
     s13[]={"ROTA_PHI"},
     s14[]={"RWSB"},
     s15[]={"AWSB"},
     s16[]={"BETA2B"},  
     s17[]={"BETA4B"},
     s18[]={"ROTB_THETA"},
     s19[]={"ROTB_PHI"},
     s20[]={"SNN"},
     s21[]={"SBIN"},
     s22[]={"ALPHA"},
     s23[]={"Uw"},
     s24[]={"Ubin"},
     s25[]={"Vw"},
     s26[]={"Vbin"},
     s27[]={"CD"},
     s28[]={"MODEL"},
     s29[]={"ISEED"},
     s30[]={"BTOT"},
     s31[]={"W0"},
     s32[]={"W1"},
     s33[]={"RDS0"},
     s34[]={"RDS1"},
     s35[]={"SHIFT"},
     s36[]={"RET"},
     s37[]={"DW"},
     s38[]={"DBIN"},
     s39[]={"WFA"},
     s40[]={"WFB"},
     s41[]={"FULL"},
     s42[]={"FBIN"},
     s43[]={"DOBIN"},
     s44[]={"GA"},
     s45[]={"FILES"},
     s46[]={"ARANK"},
     s47[]={"PP"},
     s48[]={"RO"},
     s49[]={"MAXYRAP"},
     s50[]={"FBRAP"},
     s51[]={"RCHA"},  
     s52[]={"RCHB"},  
     s53[]={"RCHP"},
     s54[]={"OMEGA"},
     s55[]={"GAMA"};	 

	#if(_nnwp_==0)
        NNWP=0;
	#elif(_nnwp_==1)
	    NNWP=1;
    #elif(_nnwp_==2)
	    NNWP=2;
	#endif


cout << "parameters reset from default in " << inpfile << " :" << endl;

//! scan the input file for the parameters reset from the default values
Int_t fBTOT=0; double v; 
while (!feof(in)){
  fscanf(in,"%s %lf",t,&v);
  if ((!feof(in)) && (t[0]!='#') && (t[0]!='I')) cout << t << "\t"<< v <<endl;

  if ((!feof(in)) && (t[0]!='#') && (t[0]=='I')) cout << t << "\t"<< UInt_t(v) <<endl;
  if (!strcmp(t,s1)) {EVENTS=Int_t(v);goto dalej;} // number of events
  if (!strcmp(t,s2)) {NBIN=Int_t(v);goto dalej;}   // number of bins in histograms
  if (!strcmp(t,s3)) {NUMA=Int_t(v);goto dalej;}   // mass number of nucleus A
  if (!strcmp(t,s4)) {NUMB=Int_t(v);goto dalej;}   // mass number of nucleus B
  if (!strcmp(t,s5)) {WMIN=Int_t(v);goto dalej;}   // minimum number of wounded nucleons to record event
  if (!strcmp(t,s6)) {BMIN=v;goto dalej;}        // minimum impact parameter
  if (!strcmp(t,s7)) {BMAX=v;goto dalej;}        // maximum impact parameter
  if (!strcmp(t,s8)) {RWSA=v;goto dalej;}        // Woods-Saxon radius on nucleus A
  if (!strcmp(t,s9)) {AWSA=v;goto dalej;}        // Wood-Saxon surface thickness parameter for nucleus A
  if (!strcmp(t,s10)) {BETA2A=v;goto dalej;}     // Wood-Saxon deformation parameter beta2 for nucleus A  
  if (!strcmp(t,s11)) {BETA4A=v;goto dalej;}     // Wood-Saxon deformation parameter beta4 for nucleus A  
  if (!strcmp(t,s12)) {ROTA_THETA=v;goto dalej;}   // Parameter controlling rotation of nucleus A in XZ plane (polar angle THETA)
  if (!strcmp(t,s13)) {ROTA_PHI=v;goto dalej;}     // Parameter controlling rotation of nucleus A in XY plane (azimuthal angle PHI) 
  if (!strcmp(t,s14)) {RWSB=v;goto dalej;}       // Woods-Saxon radius on nucleus B
  if (!strcmp(t,s15)) {AWSB=v;goto dalej;}       // Wood-Saxon surface thickness parameter for nucleus B
  if (!strcmp(t,s16)) {BETA2B=v;goto dalej;}     // Wood-Saxon deformation parameter beta2 for nucleus B  
  if (!strcmp(t,s17)) {BETA4B=v;goto dalej;}     // Wood-Saxon deformation parameter beta4 for nucleus B  
  if (!strcmp(t,s18)) {ROTB_THETA=v;goto dalej;}   // Parameter controlling rotation of nucleus B in XZ plane (polar angle THETA)
  if (!strcmp(t,s19)) {ROTB_PHI=v;goto dalej;}     // Parameter controlling rotation of nucleus B in XY plane (azimuthal angle PHI)
  if (!strcmp(t,s20)) {SNN=v;goto dalej;}        // wounding cross section
  if (!strcmp(t,s21)) {SBIN=v;goto dalej;}       // binary cross section
  if (!strcmp(t,s22)) {ALPHA=v;goto dalej;}      // the mixing parameter
  if (!strcmp(t,s23)) {Uw=v;goto dalej;}         // parameter for the distribution superimposed over the wounded nucleons
  if (!strcmp(t,s24)) {Ubin=v;goto dalej;}       // parameter for the distribution superimposed over the binary collisions
  if (!strcmp(t,s25)) {Vw=v;goto dalej;}         // parameter for the distribution superimposed over the wounded nucleons
  if (!strcmp(t,s26)) {Vbin=v;goto dalej;}       // parameter for the distribution superimposed over the binary collisions
  if (!strcmp(t,s27)) {CD=v;goto dalej;}         // the expulsion distance
  if (!strcmp(t,s28)) {MODEL=Int_t(v);goto dalej;}     // switch for the superimposed distribution
  if (!strcmp(t,s29)) {ISEED=UInt_t(v);goto dalej;}  // seed for the random-number generator
  if (!strcmp(t,s30)) {fBTOT=1;BTOT=v;goto dalej;}   // range for histograms
  if (!strcmp(t,s31)) {W0=Int_t(v);goto dalej;}     // minimum number of wounded nucleons
  if (!strcmp(t,s32)) {W1=Int_t(v);goto dalej;}     // maximum number of wounded nucleons
  if (!strcmp(t,s33)) {RDS0=v;goto dalej;}        // minimum RDS (relative deposited strength, see the paper)
  if (!strcmp(t,s34)) {RDS1=v;goto dalej;}        // maximum RDS
  if (!strcmp(t,s35)) {SHIFT=Int_t(v);goto dalej;}  // parameter controlling the shift of the fixed-axes distributions to c.m. frame 
  if (!strcmp(t,s36)) {RET=Int_t(v);goto dalej;}    // parameter controlling nuclear density generation
  if (!strcmp(t,s37)) {DW=v;goto dalej;}          // width of the distribution of displacement of the location of source for wounded nucleons
  if (!strcmp(t,s38)) {DBIN=v;goto dalej;}        // width of the distribution of displacement of the location of source for binary collisions
  if (!strcmp(t,s39)) {WFA=v;goto dalej;}         // the w parameter of the Fermi distribution for nucleus A
  if (!strcmp(t,s40)) {WFB=v;goto dalej;}         // the w parameter of the Fermi distribution for nucleus B
  if (!strcmp(t,s41)) {FULL=Int_t(v);goto dalej;}   // parameter controlling generation of the full event tree 
  if (!strcmp(t,s42)) {FBIN=Int_t(v);goto dalej;}   // number of histogram bins in the the phi angle
  if (!strcmp(t,s43)) {DOBIN=Int_t(v);goto dalej;}  // if 1, generate binary collision also for the pure wounded-nucleon model
  if (!strcmp(t,s44)) {GA=v;goto dalej;}          // value of the Gaussian NN wounding profile at the origin
  if (!strcmp(t,s45)) {FILES=Int_t(v);goto dalej;}  // not used
  if (!strcmp(t,s46)) {ARANK=Int_t(v);goto dalej;}  // rank of the Fourier moment for the forward-backward analysis
  if (!strcmp(t,s47)) {PP=Int_t(v);goto dalej;}  // power of the transverse radius in the Fourier moments 
  if (!strcmp(t,s48)) {RO=Int_t(v);goto dalej;}  // rank of the rotation axes (0 - rotation rank = rank of the Fourier moment)
  if (!strcmp(t,s49)) {MAXYRAP=v;goto dalej;}          // maximum absolute value of the y coordinate in the x-y-rapidity histogram  
  if (!strcmp(t,s50)) {FBRAP=v;goto dalej;}          // rapidity for the forward-backward analysis 
  if (!strcmp(t,s51)) {RCHA=v;goto dalej;}          // harmonic oscillator shell model density mean squared charge radii of  A-nucleus
  if (!strcmp(t,s52)) {RCHB=v;goto dalej;}          // harmonic oscillator shell model density mean squared charge radii of  B-nucleus
  if (!strcmp(t,s53)) {RCHP=v;goto dalej;}          // harmonic oscillator shell model density mean squared charge radii of  proton 
  if (!strcmp(t,s54)) {OMEGA=v;goto dalej;}          // relative variance of cross-section fluctuations
  if (!strcmp(t,s55)) {GAMA=v;goto dalej;}          // gamma wounding profile parameter (height at the origin)
 dalej:;
 };

if(RET==1){cout << "Return to beginning algorithm no longer supported, use RET=0 in the input" << endl; exit(0);};

#if(_files_)
FILES=1;
#else
FILES=0;
#endif

//! correct wrong input
// if ((FILES==1) && CD*CD > 0.000001){CD=0; cout << "correction: CD=" << CD << endl;};
 if ((MODEL!=0) && (MODEL!=1) && (MODEL!=2) && (MODEL!=3)) {MODEL=0; cout << "correction: MODEL=" << MODEL << endl;};
 if (BMIN<0) {BMIN=0; cout << "correction: BMIN=" << BMIN << endl;};
 if (BMAX<BMIN) {BMAX=BMIN; cout << "correction: BMAX=" << BMAX << endl;};
 if ((ALPHA<0)||(ALPHA>1)) {ALPHA=0; cout << "correction: ALPHA=" << ALPHA << endl;};
 if ((SHIFT!=0) && (SHIFT!=1)) {SHIFT=0; cout << "correction: SHIFT=" << SHIFT << endl;};
 if (W0<WMIN) {W0=WMIN; cout << "correction: W0=" << W0 << endl;};
 if (W1<W0) {W1=W0; cout << "correction: W1=" << W1 << endl;};
 if (RDS1<RDS0) {RDS1=RDS0+1; cout << "correction: RDS1=" << RDS1 << endl;};
 if (Uw<0) {Uw=1; cout << "correction: Uw=" << Uw << endl;};
 if (Ubin<0) {Ubin=1; cout << "correction: Ubin=" << Ubin << endl;};
 if (Vw<0) {Vw=1; cout << "correction: Vw=" << Vw << endl;};
 if (Vbin<0) {Vbin=1; cout << "correction: Vbin=" << Vbin << endl;};
 if (CD<0) {CD=-CD; cout << "correction: CD=" << CD << endl;};
 if ( (NUMA<3) && (NUMB <3) )
    {cout << "GLISSANDO 2 is designed only for p(d)+A and A+A collisions of nuclei heavier than deuteron " << endl; 
exit(0); };
 if (((ROTA_THETA<0.0) || (ROTA_THETA>180.0)) && (ROTA_THETA!=-1.0))
    {cout << "ROTA_THETA out of range. It should be in range 0.0 - 180.0 or equal -1.0 (random rotation) "<< endl; 
exit(0); };
 if (((ROTA_PHI<0.0) || (ROTA_PHI>360.0)) && (ROTA_PHI!=-1.0)) 
    {cout << "ROTA_PHI out of range. It should be in range 0.0 - 360.0 or equal -1.0 (random rotation) "<< endl;
exit(0); }; 
 if (((ROTB_THETA<0.0) || (ROTB_THETA>180.0)) && (ROTB_THETA!=-1.0))
    {cout << "ROTB_THETA out of range. It should be in range 0.0 - 180.0 or equal -1.0 (random rotation) "<< endl;
exit(0); };
 if (((ROTB_PHI<0.0) || (ROTB_PHI>360.0)) && (ROTB_PHI!=-1.0))
    {cout << "ROTB_PHI out of range. It should be in range 0.0 - 360.0 or equal -1.0 (random rotation) "<< endl;
exit(0); };
//! see the paper for the discussion of parametrizatiations of the nuclear distributions
if (NUMA>16){
if ((RWSA<=0) && ((CD-0.9)*(CD-0.9)<0.001) && (RET==0)) {
   RWSA=1.1*pow(NUMA,1./3.)-0.656*pow(NUMA,-1./3.); AWSA=0.459;
   cout << "Woods-Saxon parameters: RWSA=" << RWSA << "fm, AWSA=0.459fm (see the paper)" << endl;};
}
if (NUMB>16){  
if ((RWSB<=0) && ((CD-0.9)*(CD-0.9)<0.001) && (RET==0)) {
   RWSB=1.1*pow(NUMB,1./3.)-0.656*pow(NUMB,-1./3.); AWSB=0.459;
   cout << "Woods-Saxon parameters: RWSB=" << RWSB << "fm, AWSB=0.459fm (see the paper)" << endl;};
}
/*
if (NUMA>16){  
if ((RWSA<=0) && ((CD-0.8)*(CD-0.8)<0.001) && (RET==0)) {
   RWSA=1.103*pow(NUMA,1./3.)-0.55*pow(NUMA,-1./3.); AWSA=0.455;
   cout << "Woods-Saxon parameters: RWSA=" << RWSA << "fm, AWSA=0.455fm (see the paper)" << endl;};
}
if (NUMB>16){  
if ((RWSB<=0) && ((CD-0.8)*(CD-0.8)<0.001) && (RET==0)) {
   RWSB=1.103*pow(NUMB,1./3.)-0.55*pow(NUMB,-1./3.); AWSB=0.455;
   cout << "Woods-Saxon parameters: RWSB=" << RWSB << "fm, AWSB=0.455fm (see the paper)" << endl;};
}
if (NUMA>16){  
if ((RWSA<=0) && ((CD-0.4)*(CD-0.4)<0.001) && (RET==0)) {
   RWSA=1.113*pow(NUMA,1./3.)-0.277*pow(NUMA,-1./3.); AWSA=0.45;
   cout << "Woods-Saxon parameters: RWSA=" << RWSA << "fm, AWSA=0.45fm (see the paper)" << endl;};
}
if (NUMB>16){  
if ((RWSB<=0) && ((CD-0.4)*(CD-0.4)<0.001) && (RET==0)) {
   RWSB=1.113*pow(NUMB,1./3.)-0.277*pow(NUMB,-1./3.); AWSB=0.45; 
   cout << "Woods-Saxon parameters: RWSB=" << RWSB << "fm, AWSA=0.45fm (see the paper)" << endl;};
}
*/
if (NUMA>16){  
if ((RWSA<=0) && (CD*CD<0.001) && (RET==0)) {
   RWSA=1.114*pow(NUMA,1./3.)-0.246*pow(NUMA,-1./3.); AWSA=0.45;
   cout << "Woods-Saxon parameters: RWSA=" << RWSA << "fm, AWSA=0.45fm (see the paper)" << endl;};
}
if (NUMB>16){  
if ((RWSB<=0) && (CD*CD<0.001) && (RET==0)) {
   RWSB=1.114*pow(NUMB,1./3.)-0.246*pow(NUMB,-1./3.); AWSB=0.45; 
   cout << "Woods-Saxon parameters: RWSB=" << RWSB << "fm, AWSA=0.45fm (see the paper)" << endl;};
}

#if(_files_==0)
if (NUMA>16){  
if ((RWSA<=0 || RWSB<=0)) {cout << "Correct radii! (no formula available for the chosen CD)"<< endl << endl; exit(0); };
}
if (NUMB>16){  
if ((RWSA<=0 || RWSB<=0)) {cout << "Correct radii! (no formula available for the chosen CD)"<< endl << endl; exit(0); };
}
#endif

if ((MODEL==3) && ((Vw<=Uw) || (Vbin<=Ubin))) {cout << "Correct variance of Negative Binomial distribution! (V>U)"<< endl << endl; exit(0); };

if (NUMA!=NUMB && SHIFT==0) {SHIFT=1; cout << "Reset to SHIFT=1 for collisions of different nuclei" << endl;};

#if(_nnwp_==2)
if ((OMEGA<=0.0) || (OMEGA>=1.0)){cout << "Correct OMEGA! OMEGA IN(0,1)"<< endl << endl; exit(0); };
if ((GAMA<=0.0) || (GAMA>1.0)){cout << "Correct GAMA! GAMA IN(0,1]"<< endl << endl; exit(0); };
#endif
// set the range for the histograms
if (fBTOT!=1) BTOT = fmax(RWSA,RWSB)+AWSA+AWSB;
};

//! echo parameters to the output
void echopar(){
cout << endl;
if(ISEED1==0){cout << "random ";} else { cout << "fixed ";};
cout << "seed: " << raa.GetSeed()  << ", number of events: " << EVENTS << endl;
#if(_profile_)
   cout << "generates the nuclear density profile for nuclei A and B" << endl;
#endif
cout << NUMA <<"+" << NUMB;

if((NUMA>16) && (FILES!=1)){cout << ", RA="  << RWSA << "fm, aA=" << AWSA << "fm";};
if((NUMB>16) && (FILES!=1)){cout << ", RB="  << RWSB << "fm, aB=" << AWSB << "fm";}
if((NUMA>16) && (FILES!=1)){cout << ", dA=" << CD  << "fm";};
if((NUMB>16) && (FILES!=1)){cout << ", dB=" << CD  << "fm";};

cout << endl;

if((BETA2A != 0) || (BETA4A != 0)) cout << "nuclear deformation parameters for nucleus A: BETA2=" << BETA2A << ", BETA4=" << BETA4A << endl;
if((BETA2A != 0) || (BETA4A != 0)) cout << "rotation parameters for nucleus A: ROTA_THETA=" << ROTA_THETA << ", ROTA_PHI=" << ROTA_PHI << endl;
if((BETA4B != 0) || (BETA4B != 0)) cout << "nuclear deformation parameters for nucleus B: BETA2=" << BETA2B << ", BETA4=" << BETA4B << endl;
if((BETA2A != 0) || (BETA4A != 0)) cout << "rotation parameters for nucleus B: ROTB_THETA=" << ROTB_THETA << ", ROTB_PHI=" << ROTB_PHI << endl;

if((WFA*WFA>0.00000001) && (FILES!=1)) {cout << ", wA="<< WFA;}
if((WFB*WFB>0.00000001) && (FILES!=1)) {cout << ", wB="<< WFB;}
cout << endl; 
if(ALPHA==0){cout << "wounded nucleon model: sig_w=" << SNN << "mb" << endl;};
if(ALPHA==0 && DOBIN==1){cout << "   (binary collisions counted)" << endl;};
if(ALPHA==0 && DOBIN!=1){cout << "   (binary collisions not counted)" << endl;};
if(ALPHA>0 && SBIN>=SNN){cout << "mixed model: sig_w=" << SNN << "mb, sig_bin=" << SBIN << "mb, alpha=" << ALPHA << endl;};
if(ALPHA>0 && SBIN<SNN){cout << "mixed model with hotspots: sig_w=" << SNN << "mb, sig_bin=" << SBIN << "mb, alpha=" << ALPHA << endl;};
if(MODEL==1){cout << "overlaid Poisson distribution with parameters " << Uw << " (wounded) and " << Ubin << " (binary)" << endl;};
if(MODEL==2){cout << "overlaid Gamma distribution with parameters " << Uw << " (wounded) and " << Ubin << " (binary)" << endl;};
if(MODEL==3){cout << "overlaid Negative binomial distribution with parameters " << Uw << ", "<< Vw << " (wounded) and " << Ubin << ", " << Vbin << " (binary)" << endl;};

// harmonic oscillator shell model
if ((NUMA<17) && (NUMA>2)){cout << "harmonic oscillator shell model density"<<endl; 
                           cout <<" mean squared charge radius of nucleus A: " <<RCHA << " fm^2"<<endl;};
if ((NUMB<17) && (NUMB>2)){cout << "harmonic oscillator shell model density"<<endl;
                           cout << "mean squared charge radius of nucleus B: " <<RCHB << " fm^2"<<endl;};
if ((NUMA<17) && (NUMA>2) && (FILES!=1)){cout << " dA=" << CD  << "fm ";};
if ((NUMB<17) && (NUMB>2) && (FILES!=1)){cout << " dB=" << CD  << "fm ";};

cout << endl;

#if(_nnwp_==0)
   cout << "hard sphere NN collision profile" << endl;
#elif(_nnwp_==1)
   cout << "Gaussian NN collision profile, A=" << GA << endl;
#elif(_nnwp_==2)
    cout << "gamma NN collision profile, G="<< GAMA <<" omega="<< OMEGA << endl;
#endif

if(RO==0){cout << "rank of rotation corresponds to the rank of the given Fourier moment" << endl;};
if(RO>0){cout << "rank of rotation fixed to " << RO << " for all Fourier moments" << endl;};
if(PP==-1) cout << "power of transverse radius in eccentricities = rank (see the paper)" << endl;
else cout << "power of transverse radius in eccentricities =" << PP << endl;
cout << "window: b_min=" << BMIN << "fm, b_max=" << BMAX << "fm";
if(W1!=1000 || W0!=2){cout << ", Nw_min=" << W0 << ", Nw_max=" << W1;};
if(RDS1!=100000 || RDS0!=0){cout << ", RDS_min=" << RDS0 << ", RDS_max=" << RDS1;};
cout << endl;
if(SHIFT==0){cout << "(fireball not shifted to its center of mass)" << endl;};
// if(SHIFT==1){cout << "(fixed-axes coordinates shifted to the c.m. frame of the fireball)" << endl;};
if(CD>0.0){
	if (RET==1){
	   cout << "return-to-beginning algorithm (slow, recommended to use RET=0)" << endl;} 
//	else {cout << "fix-last algorithm" << endl;};
}; 
if(DW>0.0 || DBIN >0.0) {cout << "source dispersion parameter: wounded=" << DW << "fm, binary=" << DBIN << "fm";}; 
cout << endl; 
if(FULL){cout << "full event tree generated (THIS GENERATES A LOT OF DATA, set FULL=0 to switch off)" << endl;};
};


/*************************************
 declaration of counters and variables 
*************************************/

counter2 
estd,		  //!< counter for epsilon standard (fixed-axes), <r^2 cos(2 phi)>/<r^2>
epart1,           //!< counter for epsilon participant (variable-axes), <r^3 cos(phi)>/<r^3>
epart,            //!< counter for epsilon participant (variable-axes), <r^n cos(2 phi)>/<r^n>
estd3,            //!< counter for fixed-axes <r^n cos(3 phi)>/<r^n>
epart3,           //!< counter for variable-axes <r^n cos(3 phi)>/<r^n>
estd4,            //!< counter for fixed-axes <r^n cos(4 phi)>/<r^n>
epart4,           //!< counter for variable-axes <r^n cos(4 phi)>/<r^n>
estd5,            //!< counter for fixed-axes <r^n cos(5 phi)>/<r^n>
epart5,           //!< counter for variable-axes <r^n cos(5 phi)>/<r^n>
estd6,            //!< counter for fixed-axes <r^n cos(6 phi)>/<r^n>
epart6,           //!< counter for variable-axes <r^n cos(6 phi)>/<r^n>
nwounded,         //!< counter for number of wounded nucleons
nbinary,          //!< counter for number of binary collisions
nhot,             //!< counter for number of hot-spots
nweight;          //!< counter for relative deposited strength (RDS)

counter_2D 
angles;           //!< counter for forward-backward reaction-plane angle correlations

Int_t evall,      //!< number of all attempted event
       roo,       //!< Fourier rank for the rotation axis
       ppp,       //!< power of the weight in eccentricty definition
        kk;       //!< number of the current event

Float_t 
d,         //!< the wounding distance
dbin,      //!< the binary-collision distance 
b,         //!< impact parameter 
sitot,     //!< the total A+B cross section in the acceptance window
sirad,     //!< equivalent hard-sphere radius for the cross section
rwA,       //!< number of wounded nucleons in A
rwB,       //!< number of wounded nucleons in B
rwAB,      //!< number of all wounded nucleons
rbin,      //!< number of binary collisions
rhotspot,  //!< number of hot-spots
rpa,       //!< relative deposited strength (RDS)
sizeav,	   //!< size
es,        //!< epsilon standard (fixed-axes), <r^2 cos(2 phi)>/<r^2>
ess,       //!< fixed-axes sine moment, <r^2 sin(2 phi)>/<r^2>
ep1,       //!< epsilon participant (variable-axes), <r^3 cos(phi)>/<r^3>
ep1s,      //!< variable-axes sine moment, <r^3 sin(phi)>/<r^3> 
ep,        //!< epsilon participant (variable-axes), <r^n cos(2 phi)>/<r^n>
eps,       //!< variable-axes sine moment, <r^n sin(2 phi)>/<r^n> 
es3,       //!< fixed-axes <r^n cos(3 phi)>/<r^n>
es3s,      //!< fixed-axes sine moment, <r^n sin(3 phi)>/<r^n>
ep3,       //!< variable-axes <r^2 cos(3 phi)>/<r^2>
ep3s,      //!< variable-axes sine moment, <r^n sin(3 phi)>/<r^n>
es4,       //!< fixed-axes <r^n cos(4 phi)>/<r^n>
es4s,      //!< fixed-axes sine moment, <r^n sin(4 phi)>/<r^n>
ep4,       //!< variable-axes <r^n cos(4 phi)>/<r^n>
ep4s,      //!< variable-axes sine moment, <r^n sin(4 phi)>/<r^n>
es5,       //!< fixed-axes <r^n cos(5 phi)>/<r^n>
es5s,      //!< fixed-axes sine moment, <r^n sin(5 phi)>/<r^n>
ep5,       //!< variable-axes <r^n cos(5 phi)>/<r^n>
ep5s,      //!< variable-axes sine moment, <r^n sin(5 phi)>/<r^n>
es6,       //!< fixed-axes <r^n cos(6 phi)>/<r^n>
es6s,      //!< fixed-axes sine moment, <r^n sin(6 phi)>/<r^n>
ep6,       //!< variable-axes <r^n cos(6 phi)>/<r^n>
ep6s,      //!< variable-axes sine moment, <r^n sin(6 phi)>/<r^n>
phirot,    //!< rotation angle maximizing the second Fourier moment 
phirot_plus,    //!< rotation angle maximizing the second Fourier moment - increased rapidity
phirot_minus,    //!< rotation angle maximizing the second Fourier moment - decreased rapidity
phirot3,   //!< rotation angle maximizing the third Fourier moment 
phirot4,   //!< rotation angle maximizing the fourth Fourier moment 
phirot5,   //!< rotation angle maximizing the fifth Fourier moment 
phirot6,   //!< rotation angle maximizing the sixth Fourier moment 
xx,        //!< center-of-mass x coordinate
yy,        //!< center-of-mass y coordinate
xeps,      //!< average es 
xseps,     //!< standard deviation of es
xepp,      //!< average ep 
xsepp;     //!< standard deviation of ep


//! reset the counters used to store physical quantities in the event
void reset_counters(                  
                   ){
estd.reset(); epart1.reset(); epart.reset(); estd3.reset(); epart3.reset(); estd4.reset(); epart4.reset(); estd5.reset(); epart5.reset(); 
estd6.reset(); epart6.reset(); nwounded.reset(); nbinary.reset(); nhot.reset(); nweight.reset(); angles.reset();
// nwAwB.reset(); 
                     };


//! class storing the trees and histograms
/*!
Class for storage of ROOT structures (trees, histograms) used for later off-line analysis within ROOT or other codes 
*/
class tr_his_c {
public:

// trees
    TTree *param,      //!< parameters 
          *phys,       //!< A+B cross section and other physical results
          *full_event, //!< full info on the event (positions and RDS of the sources)
          *tree;       //!< basic physical results

// histograms
	TH2D *xyhist,        //!< cartesian fixed-axes distribution
	     *xyhist_nuclA,  //!< cartesian fixed-axes distribution of density in nucleus A
         *xyhist_nuclB,  //!< cartesian fixed-axes distribution of density in nucleus B
         *rcostheta_nuclA,  //!< (r,cos(theta)) distribution of density in nucleus A
         *rcostheta_nuclB,  //!< (r,cos(theta)) distribution of density in nucleus B
         *xyhist_mantle, //!< cartesian fixed-axes mantle distribution
         *xyhist_core,   //!< cartesian fixed-axes core distribution
         *xyhistr,       //!< cartesian variable-axes distribution

             *c0hist, //!< polar fixed-axes distribution of cos(phi)
             *c2hist, //!< polar fixed-axes distribution of cos(2 phi)
             *c3hist, //!< polar fixed-axes distribution of cos(3 phi)
             *c4hist, //!< polar fixed-axes distribution of cos(4 phi) 
             *c5hist, //!< polar fixed-axes distribution of cos(5 phi)
             *c6hist, //!< polar fixed-axes distribution of cos(6 phi)

	     *c0rhist, //!< polar variable-axes distribution of cos(phi)
             *c2rhist, //!< polar variable-axes distribution of cos(2 phi) 
             *c3rhist, //!< polar variable-axes distribution of cos(3 phi) 
             *c4rhist, //!< polar variable-axes distribution of cos(4 phi) 
             *c5rhist, //!< polar variable-axes distribution of cos(5 phi) 
             *c6rhist, //!< polar variable-axes distribution of cos(6 phi)

	     *s3hist,  //!< polar fixed-axes distribution of sin(3 phi) 
             *s3rhist; //!< polar variable-axes distribution of sin(3 phi)
 
// histograms for the dependence of certain quantities on the total number of wounded nucleons
	TH1D *nx,      //!< center-of-mass x coordinate of the source distribution vs. Nw
             *nx2,     //!< square of cm x coordinate, then its variance, vs. Nw 
             *ny,      //!< center-of-mass y coordinate of the source distribution vs. Nw
             *ny2,     //!< square of cm y coordinate, then its variance, vs. Nw
             *nsize,   //!< size vs. Nw
             *nsize2,  //!< square of size, then its variance/size^2, vs. Nw
             *neps,    //!< fixed-axes eccentricity vs. Nw
             *neps2,   //!< square of fixed-axes eccentricity, then its variance, vs. Nw
             *neps4,   //!< fixed-axes fourth moment vs. Nw
             *nepsp,   //!< variable-axes eccentricity vs. Nw
             *nepsp2,  //!< square of variable-axes eccentricity, then its variance, vs. Nw
             *nepsp22, //!< fourth power of variable-axes eccentricity vs. Nw
             *nepsp4,  //!< variable-axes fourth moment vs. Nw
             *nuni,    //!< frequency of Nw, i.e. histogram of unity vs. Nw
             *nepsb,   //!< fixed-axes eccentricity vs. b 
             *neps2b,  //!< square of fixed-axes eccentricity, then its variance, vs. b
             *nepspb,  //!< variable-axes eccentricity vs. b 
             *nepsp2b, //!< square of variable-axes eccentricity, then its variance, vs. b
             *nunib,   //!< frequency of b, i.e. histogram of unity vs. b
             *nwb,     //!< number of wounded nucleons in nucleus B vs. total number of wounded nucleons
             *nw2b;    //!< square of the number of wounded nucleons in nucleus B, then its variance, vs. total number of wounded nucleons
 
// histograms for fluctuations of number of wounded nucleons and RDS
	TH1D *nwei,   //!< RDS vs. number of wounded nucleons in nucleus A
             *nwei2,  //!< square of RDS, then its variance, vs. number of wounded nucleons in nucleus A
             *ntarg,  //!< number of wounded nucleons in nucles B vs. number of wounded nucleons in nucleus A
             *ntarg2, //!< square of the number of wounded nucleons in nucles B, then its variance, vs. number of wounded nucleons in nucleus A
             *nbinar, //!< number of binary collisions vs. number of wounded nucleons in nucleus A 
             *nbinar2,//!< square of the number of binary collisions, then its variance, vs. number of wounded nucleons in nucleus A  
             *nunp;   //!< frequency of the number of wounded nucleons in nucleus A


// histograms for nuclear profiles, correlations, and the weight distribution
    TH1D  *radA,      //!< one-body radial distribution in the nucleus A 
	  *radB,      //!< one-body radial distribution in the nucleus B 
          *rrelA,     //!< distance between the pair of nucleons in the nucleus A 
	  *rrelB,     //!< distance between the pair of nucleons in the nucleus B 
          *rrel_u,   //!< uncorrelated distance between the pair of nucleons in the nucleus (one nucleon from A, the other one from B) 
          *weih,     //!< the distribution overlaid on the wounded nucleons
          *weih_bin, //!< the distribution overlaid over binary collisions
          *wpro;     //!< the wounding profile

// histogram for radial distributions
TH1D *c0hp,  //!< fixed-axes radial profile f_0
     *c2hp,  //!< fixed-axes radial profile f_2
     *c3hp,  //!< fixed-axes radial profile f_3
     *s3hp,  //!< fixed-axes radial profile for the sine moment, g_3
     *c4hp,  //!< fixed-axes radial profile f_4
     *c5hp,  //!< fixed-axes radial profile f_5
     *c6hp,  //!< fixed-axes radial profile f_6
     *c0rhp, //!< variable-axes radial profile f_0
     *c2rhp, //!< variable-axes radial profile f_2 
     *s3rhp, //!< variable-axes radial profile f_3 
     *c3rhp, //!< variable-axes radial profile for the sine moment, g_3 
     *c4rhp, //!< variable-axes radial profile f_4
     *c5rhp, //!< variable-axes radial profile f_5 
     *c6rhp; //!< variable-axes radial profile f_6

//! initialize the histograms
    void init(){
        param = new TTree("param","param tree"); //!< tree storing parameters
        param->Branch("EVENTS",&EVENTS,"EVENTS/I");
        param->Branch("NBIN",&NBIN,"NBIN/I");
        param->Branch("FBIN",&FBIN,"FBIN/I");
        param->Branch("NUMA",&NUMA,"NUMA/I");
	param->Branch("BETA2A",&BETA2A,"BETA2A/F");   
        param->Branch("BETA4A",&BETA4A,"BETA4A/F");   
        param->Branch("ROTA_THETA",&ROTA_THETA,"ROTA_THETA/F");  
        param->Branch("ROTA_PHI",&ROTA_PHI,"ROTA_PHI/F");        
        param->Branch("NUMB",&NUMB,"NUMB/I");
	param->Branch("BETA2B",&BETA2B,"BETA2B/F");   
        param->Branch("BETA4B",&BETA4B,"BETA4B/F");   
        param->Branch("ROTB_THETA",&ROTB_THETA,"ROTB_THETA/F");  
        param->Branch("ROTB_PHI",&ROTB_PHI,"ROTB_PHI/F");        
        param->Branch("WMIN",&WMIN,"WMIN/I");
        param->Branch("MODEL",&MODEL,"MODEL/I");
        param->Branch("W0",&W0,"W0/I");
        param->Branch("W1",&W1,"W1/I");
        param->Branch("RDS0",&RDS0,"RDS0/F");
        param->Branch("RDS1",&RDS1,"RDS1/F");
        param->Branch("ISEED",&ISEED,"ISEED/i");
        param->Branch("BMIN",&BMIN,"BMIN/F");
        param->Branch("BMAX",&BMAX,"BMAX/F");
        param->Branch("BTOT",&BTOT,"BTOT/F");
        param->Branch("RWSA",&RWSA,"RWSA/F");
        param->Branch("AWSA",&AWSA,"AWSA/F");
        param->Branch("RWSB",&RWSB,"RWSB/F");
        param->Branch("AWSB",&AWSB,"AWSB/F");
        param->Branch("WFA",&WFA,"WFA/F");
        param->Branch("WFB",&WFB,"WFB/F");
        param->Branch("SNN",&SNN,"SNN/F");
        param->Branch("SBIN",&SBIN,"SBIN/F");
        param->Branch("ALPHA",&ALPHA,"ALPHA/F");
        param->Branch("DOBIN",&DOBIN,"DOBIN/I");
        param->Branch("Uw",&Uw,"Uw/F");
        param->Branch("Ubin",&Ubin,"Ubin/F");
        param->Branch("Vw",&Vw,"Vw/F");
        param->Branch("Vbin",&Vbin,"Vbin/F");
        param->Branch("CD",&CD,"CD/F");
        param->Branch("SHIFT",&SHIFT,"SHIFT/I");
        param->Branch("RET",&RET,"RET/I");      
        param->Branch("DW",&DW,"DW/F");
        param->Branch("DBIN",&DBIN,"DBIN/F");
        param->Branch("NNWP",&NNWP,"NNWP/I");
        param->Branch("GA",&GA,"GA/F");
        param->Branch("FILES",&FILES,"FILES/I");
        param->Branch("ARANK",&ARANK,"ARANK/I");
        param->Branch("PP",&PP,"PP/I");
        param->Branch("RO",&RO,"RO/I");
        param->Branch("MAXYRAP",&MAXYRAP,"MAXYRAP/F");
        param->Branch("FBRAP",&FBRAP,"FBRAP/F");
	param->Branch("RCHA",&RCHA,"RCHA/F");  
	param->Branch("RCHB",&RCHB,"RCHB/F");  
	param->Branch("RCHP",&RCHP,"RCHP/F");  
    param->Branch("OMEGA",&OMEGA,"OMEGA/F");  
	param->Branch("GAMA",&GAMA,"GAMA/F");  
    param->Branch("ver",&ver,"ver/F");

        phys = new TTree("phys","physical results"); //!< tree storing some physical quantities 
        phys->Branch("sitot",&sitot,"sitot/F");
//        phys->Branch("eps_fixed",&xeps,"seps/F");
        phys->Branch("eps_variable",&xepp,"sepp/F");
//        phys->Branch("sigma_eps_fixed",&xseps,"xseps/F");
        phys->Branch("sigma_eps_variable",&xsepp,"xsepp/F");

        full_event = new TTree("full_event","full event"); //!< tree storing full info on the events
        full_event -> Branch("full_source",&tSource,"X:Y:W:KK");
     
        tree =  new TTree("events","event tree");         //!< tree storing basic information on events
        tree->Branch("nwA",&rwA,"nwA");    // wounded nucleons in A
        tree->Branch("nwB",&rwB,"nwB");    // wounded nucleons in B
        tree->Branch("nwAB",&rwAB,"nwAB"); // all wounded nucleons
        tree->Branch("nbin",&rbin,"nbin"); // binary collisions
        tree->Branch("npa",&rpa,"npa");    // source RDS
        tree->Branch("b",&b,"b");          // impact parameter

// fixed-axes and variable-axes cos and sin moments (some commented out)
	tree->Branch("size",&sizeav,"size"); // <r>
//        tree->Branch("es",&es,"es");         // < r^2 cos(2 phi) >
        tree->Branch("ep1",&ep1,"ep1");      // < r^3 cos((phi-phi*)) >
        tree->Branch("ep",&ep,"ep");         // < r^n cos(2 (phi-phi*)) >
//        tree->Branch("es3",&es3,"es3");      // < r^n cos(3 phi) >
        tree->Branch("ep3",&ep3,"ep3");      // < r^n cos(3 (phi-phi*)) >
//        tree->Branch("es4",&es4,"es4");      // < r^n cos(4 phi) >
        tree->Branch("ep4",&ep4,"ep4");      // < r^n cos(4 (phi-phi*)) >
//        tree->Branch("es5",&es5,"es5");      // < r^n cos(5 phi) >
        tree->Branch("ep5",&ep5,"ep5");      // < r^n cos(5 (phi-phi*)) >
//        tree->Branch("es6",&es6,"es6");      // < r^n cos(6 phi) >
        tree->Branch("ep6",&ep6,"ep6");      // < r^n cos(6 (phi-phi*)) >

// fixed-axes and variable-axes sin moments
//        tree->Branch("ess",&ess,"ess");    // < r^2 sin(2 phi) >
        tree->Branch("eps",&eps,"eps");    // < r^n sin(2 (phi-phi*)) >
//        tree->Branch("es3s",&es3s,"es3s"); // < r^2 sin(3 phi) >
        tree->Branch("ep3s",&ep3s,"ep3s"); // < r^n sin(3 (phi-phi*)) >
//        tree->Branch("es4s",&es4s,"es4s"); // < r^2 sin(4 phi) >
        tree->Branch("ep4s",&ep4s,"ep4s"); // < r^n sin(4 (phi-phi*)) >
//        tree->Branch("es5s",&es5s,"es5s"); // < r^2 sin(5 phi) >
        tree->Branch("ep5s",&ep5s,"ep5s"); // < r^n sin(5 (phi-phi*)) >
//        tree->Branch("es6s",&es6s,"es6s"); // < r^2 sin(6 phi) >
        tree->Branch("ep6s",&ep6s,"ep6s"); // < r^n sin(6 (phi-phi*)) >

// other quantities
        tree->Branch("phir",&phirot,"phir");    // the rotation angle phi*
        tree->Branch("phi2_plus",&phirot_plus,"phi2_plus");    // the rotation angle phi*, increased rapidity
        tree->Branch("phi2_minus",&phirot_minus,"phi2_minus");    // the rotation angle phi*, decreased rapidity
        tree->Branch("phir3",&phirot3,"phir3"); // the rotation angle phi3*
        tree->Branch("phir4",&phirot4,"phir4"); // the rotation angle phi4*
        tree->Branch("phir5",&phirot5,"phir5"); // the rotation angle phi5*
        tree->Branch("phir6",&phirot6,"phir6"); // the rotation angle phi6*
        tree->Branch("xx",&xx,"x");             // x c.m. coordinate
        tree->Branch("yy",&yy,"y");             // y c.m. coordinate

        radA =  new TH1D("radA", "one-body distribution, nucleus A", 600, 0., 15.);  
		radB =  new TH1D("radB", "one-body distribution, nucleus B", 600, 0., 15.);  
        rrelA = new TH1D("rrelA", "relative distance - correlated, nucleus A", 2500, 0.000001, 25.);
		rrelB = new TH1D("rrelB", "relative distance - correlated, nucleus B", 2500, 0.000001, 25.);
        rrel_u = new TH1D("rrel_u", "relative distance - uncorrelated", 2500, 0.000001, 25.);

	weih = new TH1D("weih", "source weight distribution, wounded", 500, -0.1, 1.5);
	weih_bin = new TH1D("weih_bin", "source weight distribution, binary", 500, -0.1, 1.5);
	wpro = new TH1D("wpro", "wounding profile", 100, 0.000001, 6.);

// histograms for the 2D profiles
	xyhist   = new TH2D("xyhist",  "fixed-axes source density", NBIN, -BTOT, BTOT,NBIN,-BTOT,BTOT);
        xyhist -> SetXTitle("x"); xyhist -> SetYTitle("y");
		
//   x-y histograms of density for single (A, B) nucleus
        xyhist_nuclA   = new TH2D("xyhist_nuclA",  "fixed-axes source density NUCLA", NBIN, -BTOT, BTOT,NBIN,-BTOT,BTOT);
        xyhist_nuclA -> SetXTitle("x"); xyhist_nuclA -> SetYTitle("y");
        xyhist_nuclB   = new TH2D("xyhist_nuclB",  "fixed-axes source density NUCLB", NBIN, -BTOT, BTOT,NBIN,-BTOT,BTOT);
        xyhist_nuclB -> SetXTitle("x"); xyhist_nuclB -> SetYTitle("y");
        rcostheta_nuclA   = new TH2D("rcostheta_nuclA",  "(r,cos(theta)) source density NUCLA", NBIN, 0.0, 10.0, NBIN,-1.0,1.0); 
        rcostheta_nuclA -> SetXTitle("r"); rcostheta_nuclA -> SetYTitle("cos(theta)"); 
        rcostheta_nuclB   = new TH2D("rcostheta_nuclB",  "(r,cos(theta)) source density NUCLB", NBIN, 0.0, 10.0, NBIN,-1.0,1.0); 
        rcostheta_nuclB -> SetXTitle("r"); rcostheta_nuclB -> SetYTitle("cos(theta)"); 
		
	xyhist_mantle   = new TH2D("xyhist_mantle",  "variable-axes mantle density", NBIN, -BTOT, BTOT,NBIN,-BTOT,BTOT);
        xyhist_mantle -> SetXTitle("x"); xyhist_mantle -> SetYTitle("y");
	xyhist_core   = new TH2D("xyhist_core",  "variable-axes core density", NBIN, -BTOT, BTOT,NBIN,-BTOT,BTOT);
        xyhist_core -> SetXTitle("x"); xyhist_core -> SetYTitle("y");
	xyhistr  = new TH2D("xyhistr", "variable-axes source density", NBIN, -BTOT, BTOT,NBIN,-BTOT,BTOT);
        xyhistr -> SetXTitle("x"); xyhistr -> SetYTitle("y");       

// histograms for the 2D profiles
// fixed-axes cos
	c0hist  = new TH2D("c0hist", "f_{0}", NBIN, 0.001, BTOT,FBIN,-PI,PI);
        c0hist -> SetXTitle("r"); c0hist -> SetYTitle("#phi   ");
	c2hist  = new TH2D("c2hist", "f_{2}", NBIN, 0.001, BTOT,FBIN,-PI,PI);
        c2hist -> SetXTitle("r"); c2hist -> SetYTitle("#phi   ");
	c3hist  = new TH2D("c3hist", "f_{3}", NBIN, 0.001, BTOT,FBIN,-PI,PI);
        c3hist -> SetXTitle("r"); c3hist -> SetYTitle("#phi   ");
	c4hist  = new TH2D("c4hist", "f_{4}", NBIN, 0.001, BTOT,FBIN,-PI,PI);
        c4hist -> SetXTitle("r"); c4hist -> SetYTitle("#phi   ");
	c5hist  = new TH2D("c5hist", "f_{5}", NBIN, 0.001, BTOT,FBIN,-PI,PI);
        c5hist -> SetXTitle("r"); c5hist -> SetYTitle("#phi   ");
	c6hist  = new TH2D("c6hist", "f_{6}", NBIN, 0.001, BTOT,FBIN,-PI,PI);
        c6hist -> SetXTitle("r"); c6hist -> SetYTitle("#phi   ");

// variable-axes cos
	c0rhist  = new TH2D("c0rhist", "f^{*}_{0}", NBIN, 0.001, BTOT,FBIN,-PI,PI);
        c0rhist -> SetXTitle("r"); c0rhist -> SetYTitle("#phi   ");
	c2rhist  = new TH2D("c2rhist", "f^{*}_{2}", NBIN, 0.001, BTOT,FBIN,-PI,PI);
        c2rhist -> SetXTitle("r"); c2rhist -> SetYTitle("#phi   ");
	c3rhist  = new TH2D("c3rhist", "f^{*}_{3}", NBIN, 0.001, BTOT,FBIN,-PI,PI);
        c3rhist -> SetXTitle("r"); c3rhist -> SetYTitle("#phi   ");
	c4rhist  = new TH2D("c4rhist", "f^{*}_{4}", NBIN, 0.001, BTOT,FBIN,-PI,PI);
        c4rhist -> SetXTitle("r"); c4rhist -> SetYTitle("#phi   ");
	c5rhist  = new TH2D("c5rhist", "f^{*}_{5}", NBIN, 0.001, BTOT,FBIN,-PI,PI);
        c5rhist -> SetXTitle("r"); c5rhist -> SetYTitle("#phi   ");
	c6rhist  = new TH2D("c6rhist", "f^{*}_{6}", NBIN, 0.001, BTOT,FBIN,-PI,PI);
        c6rhist -> SetXTitle("r"); c6rhist -> SetYTitle("#phi   ");

// fixed-axes and variable-axes sin
	s3hist   = new TH2D("s3hist",      "g_{3}", NBIN, 0.001, BTOT,FBIN,-PI,PI);
        s3hist -> SetXTitle("r"); s3hist -> SetYTitle("#phi   ");
	s3rhist  = new TH2D("s3rhist", "g^{*}_{3}", NBIN, 0.001, BTOT,FBIN,-PI,PI);
        s3rhist -> SetXTitle("r"); s3rhist -> SetYTitle("#phi   ");

// histograms for the dependence of certain quantities on the number of wounded nucleons
	nx = new TH1D("nx", "<x> vs. N_{w}", NUMA+NUMB, 0.5, NUMA+NUMB+0.5);
	nx2= new TH1D("nx2", "var(x) vs. N_{w}", NUMA+NUMB, 0.5, NUMA+NUMB+0.5); 
	ny = new TH1D("ny", "<y> vs. N_{w}", NUMA+NUMB, 0.5, NUMA+NUMB+0.5);
	ny2= new TH1D("ny2", "var(y) vs. N_{w}", NUMA+NUMB, 0.5, NUMA+NUMB+0.5); 
	nsize = new TH1D("nsize", "<r> vs. N_{w}", NUMA+NUMB, 0.5, NUMA+NUMB+0.5);
	nsize2 = new TH1D("nsize2", "var(<r>)/<<r>>^2 vs. N_{w}", NUMA+NUMB, 0.5, NUMA+NUMB+0.5);
	neps = new TH1D("neps", "#epsilon_{2} vs. N_{w}", NUMA+NUMB, 0.5, NUMA+NUMB+0.5);
	neps4 = new TH1D("neps4", "#epsilon_{4} vs. N_{w}", NUMA+NUMB, 0.5, NUMA+NUMB+0.5);
	neps2= new TH1D("neps2", "var(#epsilon)/#epsilon^{2} vs. N_{w}",  NUMA+NUMB, 0.5, NUMA+NUMB+0.5); 
	nepsp = new TH1D("nepsp", "#epsilon^{*}_{2} vs. N_{w}",  NUMA+NUMB, 0.5, NUMA+NUMB+0.5);
	nepsp22 = new TH1D("nepsp22", "#epsilon^{*4}_{2} vs. N_{w}",  NUMA+NUMB, 0.5, NUMA+NUMB+0.5);
	nepsp4 = new TH1D("nepsp4", "#epsilon^{*}_{4} vs. N_{w}",  NUMA+NUMB, 0.5, NUMA+NUMB+0.5);
	nepsp2= new TH1D("nepsp2", "var(#epsilon*)/#epsilon*^{2} vs. N_{w}",  NUMA+NUMB, 0.5, NUMA+NUMB+0.5); 
	nuni = new TH1D("nuni", "event multiplicity vs. N_{w}",  NUMA+NUMB, 0.5, NUMA+NUMB+0.5);
	nepsb = new TH1D("nepsb", "#epsilon vs. b", 200, BMIN, BMAX);
	neps2b= new TH1D("neps2b", "var(#epsilon)/#epsilon^{2} vs. b", 200, BMIN, BMAX); 
	nepspb = new TH1D("nepspb", "#epsilon* vs. b", 200, BMIN, BMAX);
	nepsp2b= new TH1D("nepsp2b", "var(#epsilon*)/#epsilon*^{2} vs. b", 200, BMIN, BMAX); 
	nunib = new TH1D("nunib", "event multiplicity vs. b", 200, BMIN, BMAX);
	nwb = new TH1D("nwb", "N_{w} vs. b", 200, BMIN, BMAX);
	nw2b= new TH1D("nw2b", "var(N_{w}) vs. b", 200, BMIN, BMAX); 


// histograms for fluctuations of number of wounded nucleons and RDS
	nwei = new TH1D("nwei", "RDS/N_{w}^{A} vs. N_{w}^{A}",  NUMA, 0.5, NUMA+0.5);
//	nwei = new TH1D("nwei", "RDS/(N_{w}/2) vs. N_{w}/2",  NUMA, 0.5, NUMA+0.5);
	nwei2= new TH1D("nwei2", "var(RDS)/RDS vs. N_{w}^{A}",  NUMA, 0.5, NUMA+0.5); 
	ntarg = new TH1D("ntarg", "N_{w}^{B} vs. N_{w}^{A}",  NUMA, 0.5, NUMA+0.5);
	ntarg2= new TH1D("ntarg2", "var(N_{w}^{B})/N_{w}^{B} vs. N_{w}^{A}",  NUMA, 0.5, NUMA+0.5); 
	nbinar = new TH1D("nbinar", "N_{bin} vs. N_{w}^{A}",  NUMA, 0.5, NUMA+0.5);
	nbinar2= new TH1D("nbinar2", "var(N_{bin})/N_{bin} vs. N_{w}^{A}",  NUMA, 0.5, NUMA+0.5); 
	nunp = new TH1D("nunp", "event multiplicity vs. N_{w}^{A}",  NUMA, 0.5, NUMA+0.5);
        };

//! fill trees param and phys
    void fill(){param->Fill();phys->Fill();}; 

//! fill the main tree
    void fill_tr(){tree->Fill();};

//! projects the 2-dim histograms with polar distributions on 1-dim histograms
    void proj(){
	c0hp = c0hist->ProjectionX("c0hp");
	c2hp = c2hist->ProjectionX("c2hp");
	c3hp = c3hist->ProjectionX("c3hp");
	s3hp = s3hist->ProjectionX("s3hp");
	c4hp = c4hist->ProjectionX("c4hp");        
	c5hp = c5hist->ProjectionX("c5hp");        
	c6hp = c6hist->ProjectionX("c6hp");
	c0rhp = c0rhist->ProjectionX("c0rhp");
	c2rhp = c2rhist->ProjectionX("c2rhp");
	s3rhp = s3rhist->ProjectionX("s3rhp");
	c3rhp = c3rhist->ProjectionX("c3rhp");        
	c4rhp = c4rhist->ProjectionX("c4rhp");
	c5rhp = c5rhist->ProjectionX("c5rhp");        
	c6rhp = c6rhist->ProjectionX("c6rhp");
	};

//! calculate eccentricity, size, etc. and their fluctuations vs. number of wounded nucleons or b 
    void fill_res(){        
        nx -> Fill(rwAB,xx);
        nx2 -> Fill(rwAB,xx*xx);
        ny -> Fill(rwAB,yy);
        ny2 -> Fill(rwAB,yy*yy);

	nsize ->  Fill(rwAB,sizeav);
	nsize2 -> Fill(rwAB,sizeav*sizeav);

        neps -> Fill(rwAB,es);
        neps4 -> Fill(rwAB,es4);
        neps2 -> Fill(rwAB,es*es);
        nepsp -> Fill(rwAB,ep);
        nepsp4 -> Fill(rwAB,ep4);
        nepsp2 -> Fill(rwAB,ep*ep);
        nepsp22 -> Fill(rwAB,ep*ep*ep*ep);

        nuni -> Fill(rwAB,1);

        nepsb -> Fill(b,es);
        neps2b -> Fill(b,es*es);
        nepspb -> Fill(b,ep);
        nepsp2b -> Fill(b,ep*ep);

	nwb -> Fill(b,rwAB);
	nw2b -> Fill(b,rwAB*rwAB);
        nunib -> Fill(b,1);

// for multiplicity fluctuations
        nwei->Fill(rwA,rpa); 
//        nwei->Fill(rwAB/2.,2.*rpa/rwAB);
        nwei2->Fill(rwA,rpa*rpa);
	ntarg->Fill(rwA,rwB);
	ntarg2->Fill(rwA,rwB*rwB);
	nbinar->Fill(rwA,rbin);
	nbinar2->Fill(rwA,rbin*rbin);
	nunp->Fill(rwA,1);};

//! write out trees param, phys, full_event, and the main tree
    void write(){param->Write();phys->Write();tree->Write();
                  if(FULL){full_event->Write();};
                };   

//! write out the radial density distribution an the pair distance distribution in the nucleus
    void write_r(){radA->Write();radB->Write();rrelA->Write();rrelB->Write();rrel_u->Write();
                   xyhist_nuclA->Write();xyhist_nuclB->Write();rcostheta_nuclA->Write();rcostheta_nuclB->Write();
                  };
//! write out the overlaid distributions
    void write_w(){weih->Write();weih_bin->Write();};
 
//! write out the wounding profile
    void write_wpro(){wpro->Write();};

//! write out the histograms with the 2-dim distributions and the radial distributions of the Fourier components of the source profiles
    void write_d(){
	c0hp->Write(); 
	c2hp->Write(); 
	c3hp->Write(); 
	s3hp->Write();
	c4hp->Write(); 
	c5hp->Write(); 
	c6hp->Write(); 
	c0rhp->Write();
	c2rhp->Write(); 
	c3rhp->Write();
	s3rhp->Write(); 
	c4rhp->Write(); 
	c5rhp->Write();
	c6rhp->Write(); 

	xyhist->Write(); 
	xyhist_mantle->Write(); 
	xyhist_core->Write(); 
	xyhistr->Write();

	c0hist->Write(); 
	c0rhist->Write();};
 
//! generate histograms of eccentricities and their variance, etc., vs. the number of wounded nucleons or b
   void gen(){
	nx->Divide(nx,nuni);
        nx->Write();
        nx2->Divide(nx2,nuni);
        nx->Multiply(nx,nx);
        nx2->Add(nx,-1);

        ny->Divide(ny,nuni);
        ny->Write();
        ny2->Divide(ny2,nuni);
        ny->Multiply(ny,ny);
        ny2->Add(ny,-1);

	nsize -> Divide(nsize,nuni);
	nsize -> Write();
	nsize2-> Divide(nsize2,nuni);
	nsize -> Multiply(nsize,nsize);
	nsize2-> Add(nsize,-1);
	nsize2-> Divide(nsize2,nsize);
	nsize2-> Write();

        neps->Divide(neps,nuni);
//        neps->Write();
        neps4->Divide(neps4,nuni);
//        neps4->Write();

        neps2->Divide(neps2,nuni);
        neps->Multiply(neps,neps);
        neps2->Add(neps,-1);
        neps2->Divide(neps2,neps);

        nepsp->Divide(nepsp,nuni);
        nepsp->Write();

        nepsp22->Divide(nepsp22,nuni);
        nepsp22->Write();

        nepsp4->Divide(nepsp4,nuni);
        nepsp4->Write();

        nepsp2->Divide(nepsp2,nuni);
        nepsp->Multiply(nepsp,nepsp);
        nepsp2->Add(nepsp,-1);
        nepsp2->Divide(nepsp2,nepsp);
        	
        nepsb->Divide(nepsb,nunib);
//        nepsb->Write();
        neps2b->Divide(neps2b,nunib);
        nepsb->Multiply(nepsb,nepsb);
        neps2b->Add(nepsb,-1);
        neps2b->Divide(neps2b,nepsb);

        nepspb->Divide(nepspb,nunib);
        nepspb->Write();
        nepsp2b->Divide(nepsp2b,nunib);
        nepspb->Multiply(nepspb,nepspb);
        nepsp2b->Add(nepspb,-1);
        nepsp2b->Divide(nepsp2b,nepspb);

	nwb -> Divide(nwb,nunib);
	nwb -> Write();
	nw2b -> Divide(nw2b,nunib);
	nwb -> Multiply(nwb,nwb);
	nw2b -> Add(nwb,-1);

        nx2->Write();
        ny2->Write();
//        neps2->Write();
        nepsp2->Write();
        nuni->Write();
//        neps2b->Write();
        nepsp2b->Write();
	nw2b->Write();
        nunib->Write();

	// for multiplicity fluctuations
	nwei->Divide(nwei,nunp);
	ntarg->Divide(ntarg,nunp);
	nbinar->Divide(nbinar,nunp);
        nwei->Write();
	ntarg->Write();
	nbinar->Write();
	   
	nwei2->Divide(nwei2,nunp);
        nwei2->Divide(nwei2,nwei);
	nwei2->Add(nwei,-1);
  	nwei2->Write();
	   
	ntarg2->Divide(ntarg2,nunp);
        ntarg2->Divide(ntarg2,ntarg);
	ntarg2->Add(ntarg,-1);
  	ntarg2->Write();
	   
	nbinar2->Divide(nbinar2,nunp);
        nbinar2->Divide(nbinar2,nbinar);
	nbinar2->Add(nbinar,-1);
  	nbinar2->Write();
	   
	nunp->Write();};
};



/**************************
  Helper and print blocks
**************************/

//! print the version or brief help
void helper(
   Int_t argc, //!< number of command line parameters
   char* str   //!< string parameter (-v for version, -h for brief help)
           )
{if(argc > 1){
 if(strcmp(str,"-v")==0){cout << "GLISSANDO 2 ver: " << ver << endl; exit(0);};
 if(strcmp(str,"-h")==0){cout << endl << "Usage:" << endl;
	#if(_files_!=1)
             cout << "glissando2 <input_file> <output_Root_file>  (for random nuclear distributions)" << endl;
	#else
	     cout << "glissando2 <input_file> <output_Root_file> <nucl_dictribution_A> <nucl_distribution_B>   (for nuclear distributions read from files)" << endl;
	#endif              
                          cout << endl
                                 << "for ver. 2 see http://arxiv.org.abs/xx13.xxxx, " << endl
                                 << "for ver. 1 see Comp. Phys. Comm. 180(2009)69, http://arxiv.org.abs/0710.5731, " << endl
                                 << endl << endl; exit(0);};
};};

//! print the header in the output
void header(){
cout << endl << "*********************************************************" << endl;
cout << "GLISSANDO 2 ver. " << ver << endl << 
"ver. 2: http://arxiv.org.abs/xx13.xxxx" << endl <<
"ver. 1: Computer Physics Communications 180(2009)69, http://arxiv.org.abs/0710.5731" << endl <<
"and Phys. Rev. C81(2010)064909 for implementation of the NN correlations" << endl; 
cout << "(tested with ROOT ver. 5.28--5.34)" << endl;
cout         << "**********************************************************" << endl;
cout << "Simulation of nucleus-nucleus collisions in Glauber models" << endl; 
cout << "----------------------------------------------------------" << endl;
#if(_profile_) 
cout << "generation of the nuclear profile and NN correlations" << endl; 
cout << "----------------------------------------------------------" << endl;
#endif
#if(_weight_) 
cout << "generation of the superimposed distribution and NN collision profiles" << endl; 
cout << "---------------------------------------------------------" << endl;
#endif
#if(_rapidity_) 
cout << "generation of rapidity distributions" << endl; 
cout << "---------------------------------------------------------" << endl;
#endif
};


//! print epilog to the output
void epilog(
           ){
// output  basic results to the console
      cout << "Some quantities for the specified b, N_w, and RDS window " << endl 
           << "(+/- gives the e-by-e standard deviation):" << endl;
      if(BMIN*BMIN<0.0001) {cout << "A+B cross section = " << sitot << "mb";};
      if((NUMA>2) && (BMIN*BMIN<0.0001)) {cout << ", equiv. hard-sphere radius = " << sirad << "fm";};
      cout << endl << "efficiency (accepted/all) = " << 100.*EVENTS/evall << "%" << endl;
      cout << "N_w = " << nwounded.mean()  << "+/-" << sqrt(nwounded.var()) << endl;
//      cout << "RDS_F = " << nwAwB.mean_x()  << "+/-" << sqrt(nwAwB.var_x()) << endl;
//      cout << "RDS_B = " << nwAwB.mean_y()  << "+/-" << sqrt(nwAwB.var_y()) << endl;
//      cout << "cor(F,B) = " << nwAwB.corr()  << endl;
      if(ALPHA>0 || DOBIN==1){cout << "N_bin = " << nbinary.mean()  << "+/-" << sqrt(nbinary.var());};
      if(SNN>SBIN && nhot.mean()>0.0){cout << ", N_hotspot = " << nhot.mean()  << "+/-" << sqrt(nhot.var());};
      cout << endl; 
      cout << "relative deposited strength (RDS) = " << nweight.mean()  << "+/-" << sqrt(nweight.var()) << endl << endl; 

      cout << "participant eccentricities:" << endl;
      cout << "eps_1 = " << epart1.mean() << "+/-" << sqrt(epart1.var()) << endl;
      cout << "eps_2 = " <<  epart.mean() << "+/-" << sqrt(epart.var())  << endl;
      cout << "eps_3 = " << epart3.mean() << "+/-" << sqrt(epart3.var()) << endl;
      cout << "eps_4 = " << epart4.mean() << "+/-" << sqrt(epart4.var()) << endl;

#if(_rapidity_)
      cout << "f-b correlation of rank-" << ARANK << " rho: " << angles.corr() << " cov: " << angles.cov() << endl;
#endif
              };

//! start the time measurement
Int_t time_start(){
time_t rawtime; 
struct tm * timeinfo;
time ( &rawtime ); 
timeinfo = localtime ( &rawtime ); 
Int_t ti=Int_t(rawtime);
cout << "Start: " << asctime (timeinfo); // time stamp
cout << "--------------------------------------" << endl;
return ti;
};

//! stop the time measurement and print stamp
void time_stop(
              Int_t ts //!< time at start
              ){
time_t rawtime; 
struct tm * timeinfo;
cout << endl; time ( &rawtime ); timeinfo = localtime ( &rawtime ); 
Int_t ti=Int_t(rawtime) - ts; Int_t tig=ti/3600; ti=ti-tig*3600; Int_t tim=ti/60; ti=ti-tim*60;
cout << "Finish: " << asctime (timeinfo) << "(" << tig << "h:" << tim <<"m:" << ti << "s)" << endl;
cout << "**************************************" << endl << endl;
}


/*****************************
random generator functions
*****************************/


//! random number generator using the built-in ROOT generator, uniform on (0,1)
Float_t los() {return raa.Uniform(0,1);}; 

//! random number generator for the Woods-Saxon (or Fermi) distribution - nucleus A
Float_t rlosA(){Int_t x=0;while(x==0){Float_t rr=3.*RWSA*los();Float_t b=RWSA*RWSA*los();
	       if(b<rr*rr*(1.+WFA*rr*rr/RWSA/RWSA)/(exp((rr-RWSA)/AWSA)+1)){x=1;return rr;};};}; // Woods-Saxon/Fermi

//! random number generator for the Woods-Saxon (or Fermi) distribution - nucleus B
Float_t rlosB(){Int_t x=0;while(x==0){Float_t rr=3.*RWSB*los();Float_t b=RWSB*RWSB*los();
	       if(b<rr*rr*(1.+WFB*rr*rr/RWSB/RWSB)/(exp((rr-RWSB)/AWSB)+1)){x=1;return rr;};};}; //Woods-Saxon/Fermi

/*! random number generator for the Woods-Saxon (deformed with beta2,beta4 parameters of deformation and spherical 
  harmonics Y20, Y40 (or Fermi) distribution - nucleus A */
Float_t rlosA_def(Float_t* cth_pointerA, Float_t beta2, Float_t beta4)
               {Int_t x=0;
                while(x==0)
                {Float_t rr=3.*RWSA*los(); 
                 Float_t ctheta =2*los()-1;
                 Float_t y20=sqrt(5.0/(4.0*PI))*((3.0*ctheta*ctheta-1.0)/2.0); 
                 Float_t y40=sqrt(9.0/(4.0*PI))*((35.0*ctheta*ctheta*ctheta*ctheta-30.0*ctheta*ctheta+3.0)/8.0); 
                 Float_t WS_def=rr*rr*(1.0+WFA*rr*rr/RWSA/RWSA)/(1 + exp( (rr-RWSA)/AWSA-RWSA*(beta2*y20+beta4*y40)/AWSA)  );
                 Float_t b=RWSA*RWSA*los();
//                 cout<<" ctheta= "<<ctheta<<endl;
                 if(b<WS_def) 
{x=1;
// cout<<" rlosA_def, b= "<<b<<" rr= "<<rr<<" ctheta= "<<ctheta<<endl;
*cth_pointerA=ctheta;
return rr;};
                };
               }; // Woods-Saxon DEFORMED/Fermi


/*! random number generator for the Woods-Saxon (deformed with beta2,beta4 parameters of deformation and spherical
  harmonics Y20, Y40 (or Fermi) distribution - nucleus B */
Float_t rlosB_def(Float_t* cth_pointerB, Float_t beta2, Float_t beta4)
               {Int_t x=0;
                while(x==0)
                {Float_t rr=3.*RWSB*los(); 
                 Float_t ctheta=2*los()-1;
                 Float_t y20=sqrt(5.0/(4.0*PI))*((3.0*ctheta*ctheta-1.0)/2.0);
                 Float_t y40=sqrt(9.0/(4.0*PI))*((35.0*ctheta*ctheta*ctheta*ctheta-30.0*ctheta*ctheta+3.0)/8.0); 
                 Float_t WS_def=rr*rr*(1.0+WFB*rr*rr/RWSB/RWSB) / (1 + exp( (rr-RWSB)/AWSB-RWSB*(beta2*y20+beta4*y40)/AWSB)  );
                 Float_t b=RWSB*RWSB*los();
//                 cout<<" ctheta= "<<ctheta<<endl;
                 if(b<WS_def)
{x=1;
// cout<<" rlosB_def, b= "<<b<<" rr= "<<rr<<" ctheta= "<<ctheta<<endl;
*cth_pointerB=ctheta;
return rr;};
                };
               }; // Woods-Saxon DEFORMED/Fermi

//! random number generator for the harmonic oscillator shell model density - nucleus A
/*! The harmonic oscillator shell distribution used to generate the distance 
between nucleons in light (2<NUMA<17) nuclei (Nuclear Sizes. L. R. B. Elton, Oxford University Press, New York, 1961.)*/		   
Float_t rlosA_hos(){
Int_t x=0;
while(x==0){
 Float_t rr=10.*los();
 Float_t b=10.*los();
 Float_t ddd2=2.5-(4./(Float_t)NUMA);
 ddd2=1./ddd2;
 ddd2=ddd2*(RCHA-RCHP);
Float_t hos=1.0;
 hos=hos*(1.+(((Float_t)NUMA-4.)/6.)*rr*rr/ddd2);
 hos=hos*exp(-rr*rr/ddd2);
 if(b<rr*rr*hos){
 x=1;return rr;
};
};
}; 	   

//! random number generator for the harmonic oscillator shell model density - nucleus B	
/*! The harmonic oscillator shell distribution used to generate the distance 
between nucleons in light (2<NUMB<17) nuclei (Nuclear Sizes. L. R. B. Elton, Oxford University Press, New York, 1961.)*/		   
Float_t rlosB_hos(){
Int_t x=0;
while(x==0){
 Float_t rr=10.*los();
 Float_t b=10.*los();
 Float_t ddd2=2.5-(4./(Float_t)NUMB);
 ddd2=1./ddd2;
 ddd2=ddd2*(RCHB-RCHP);
 Float_t hos=1.0;
 hos=hos*(1.+(((Float_t)NUMB-4.)/6.)*rr*rr/ddd2);
 hos=hos*exp(-rr*rr/ddd2);
 if(b<rr*rr*hos){
 x=1; return rr;
};
};
}; 
		   
//! random number generator for the Hulthen distribution
/*! The Hulthen distribution used to generate the distance between nucleons in the deuteron */
Float_t rlos_hult(){double const ah=.228;double const bh=1.118;Int_t x=0;while(x==0){Float_t rr=15.*los();
Float_t b=.25*los(); // 0.25 is an upper bound of the probability distribution
 if(b< 2*ah*bh*(ah+bh)*(exp(-2.*ah*rr)+exp(-2.*bh*rr)-2.*exp(-(ah+bh)*rr))/(ah-bh)/(ah-bh)){x=1;return rr;};};}; //Hulthen

//! rapidity distribution with the plateau 
Float_t fg(
           Float_t rr  //!< spatial rapidity
          ){return exp(-(sqrt(rr*rr)-ETA0)*(sqrt(rr*rr)-ETA0)/(2*SIGETA*SIGETA)*max(0.f, sqrt(rr*rr)-ETA0) );};

//! function f +/- from Bozek, arXiv:1002.4999v2 [nucl-th]
Float_t fpm(
            Float_t rr //!< spatial rapidity
           ) { if(rr<-ETAM) return 0; else if(rr<ETAM) return (rr+ETAM)/(2*ETAM); else return 1.; }

//! random number generator for the rapidity distribution - wounded nucleons from nucleus A
/*! Formula for the wounded quark rapidity profile from Bozek, arXiv:1002.4999v2 [nucl-th]. */
Float_t los_rap_A(){
Int_t x=0;while(x==0) {Float_t rr=raa.Uniform(-RAPRANGE,RAPRANGE); Float_t b=raa.Uniform(0.,1.);
	       if( b < fg(rr)*fpm(rr) ) {x=1;return rr;};};
                   };

//! random number generator for the rapidity distribution - wounded nucleons from nucleus B
/*! Formula for the wounded quark rapidity profile from Bozek, arXiv:1002.4999v2 [nucl-th]. */
Float_t los_rap_B(){
Int_t x=0;while(x==0) {Float_t rr=raa.Uniform(-RAPRANGE,RAPRANGE); Float_t b=raa.Uniform(0.,1.);
	       if( b < fg(rr)*fpm(-rr) ) {x=1;return rr;};};
                   };

//! random number generator for the rapidity distribution - binary collisions
/*! Formula for the wounded quark rapidity profile from Bozek, arXiv:1002.4999v2 [nucl-th]. */
Float_t los_rap_bin(){
Int_t x=0;while(x==0) {Float_t rr=raa.Uniform(-RAPRANGE,RAPRANGE); Float_t b=raa.Uniform(0.,1.);
	       if( b < fg(rr) ) {x=1;return rr;};};
                   };

//! random number generator for the Gamma distribution
Float_t gamgen(
            Float_t a //!< the parameter a in f(x) = x^(a-1) exp(-x)/Gamma(a)
            ){
Int_t k = Int_t(a);
Float_t b = a - k;
Float_t x1 = 0.0;
Float_t u=0.0;
Float_t v=0.0;
Float_t x2=0.0;
Float_t s;

if(k==0) goto blok;
s = 1.0 ;
for (Int_t i=1;i<=k;i++)
s=s*los();
x1 = -log(s);

if (b>1e-6) goto blok;
return x1;
 
blok:x2=-log(los());

bloka:u=pow(los(),1.0/b);
v=u+pow(los(),1.0/(1.0-b));
if (v>1.0) goto bloka;
return (x1+x2*u/v);
}//gamgen


//negative binomial distribution
Int_t negbin(Float_t m,Float_t v) {
Float_t p=1.0-(m/v);
Float_t r=m*(1.0-p)/p;
Float_t x=gamgen(r)*p/(1.0-p);
//Int_t n=ipois(x);
Int_t n=raa.Poisson(x);
return n;
}//negbin


Float_t dist(
           Int_t m,  //!< case: 0 - uniform, 1 - Poisson, 2 - Gamma distribution, 3 - Negative binomial distribution
           Float_t u, //!< average of the distribution in cases 1 and 2 and 3
           Float_t v //!< variance of the distribution in case 3, v>u

          ){
 switch(m) {
 case 0: { return 1; break; }
 case 1: { return raa.Poisson(u)/u;  break; } 
 case 2: { return gamgen(u)/u; break; }
 case 3: { return negbin(u,v)/u; break;}
}
}

//! random shift of the source location
/*! The location of the source may by shifted randomly when DW>0 od DBIN>0, with the Gaussian distribution of the width w. */
Float_t disp(
          Float_t w //! average shift of the source location, Gaussian distribution 
          ){return raa.Gaus(0.,w);};

#endif

