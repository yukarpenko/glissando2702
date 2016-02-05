/** \file glissando2.cxx
 * The main file of GLISSANDO 2
 * 
 */


/*! \mainpage
                                                                                         
            GLISSANDO 2 - GLauber Initial State Simulation AND mOre...       
                          ver. 2.702, 11 October 2013                            
                                                                              
  Homepage: http://www.ujk.edu.pl/homepages/mryb/GLISSANDO/index.html      

  arXiv:1310.5475 [nucl-th] 
																			 
  Authors: 
           - Wojciech Broniowski (Wojciech.Broniowski@ifj.edu.pl)              
           - Maciej Rybczynski   (Maciej.Rybczynski@ujk.edu.pl)            
           - Grzegorz Stefanek   (Grzegorz.Stefanek@ujk.edu.pl)   
           - Piotr Bozek         (Piotr.Bozek@ifj.edu.pl)                    
                                                                             
  
  For the detailed description of ver. 1 program and further references         
  to the description of the model, please, refer to 

  Computer Physics Communications 180(2009)69, arXiv:0710.5731 [nucl-th]     
  http://arxiv.org.abs/0710.5731                           
                                                                                                                                                          
  
  GLISSANDO is a Glauber Monte-Carlo generator for early-stages of relativistic
  heavy-ion collisions, written in c++ and interfaced to ROOT. Several models
  are implemented: the wounded-nucleon model, the binary collisions model,
  the mixed model, and the model with hot-spots. The original geometric
  distribution of sources (i.e., wounded nucleons or binary collisions) in
  the transverse plane can be superimposed with a statistical distribution
  simulating the dispersion in the generated transverse energy in each individual
  collision. The program generates inter alia the variable-axes (participant)
  two-dimensional profiles of the density of sources in the transverse plane
  and their Fourier components. These profiles can be used in further analyses
  of physical phenomena, such as the the jet quenching, event-by-event hydrodynamics,
  or analyses of the elliptic flow and its fluctuations. Characteristics of the event
  (multiplicities, eccentricities, Fourier shape coefficients, etc.) are evaluated
  and stored in a ROOT file for further off-line studies. A number of scripts
  is provided for that purpose. The code can also be used for the proton-nucleus 
  and deuteron-nucleus collisions. 
                                                                            

  Version 2 of GLISSANDO offers much more functionality than version 1, 
  moreover, it is fully object-oriented, providing the user with the flexibility 
  of inspecting and, if needed, modifying the code in a simple manner. New 
  features involve:                                                          

- Parametrization of shape of light nuclei, useful in particular for simulations 
  for the NA61 experiment 
- Generation of distributions of deformed nuclei 
  according to the deformed Woods-Saxon density with default deformation 
  parameters taken from P.Moller,J.R.Nix,W.D.Mayers, and W.J.Swiatecki 
  [Nucl.Data Tables 59,185,1995].See also W.Broniowski,M.Rybczynski,G.Stefanek
  [Phys.Rev. C87, 044908, 2013;arXiv:1211.2537]                                                                              
- The possibility of feeding into the simulations the nuclear distributions 
  accounting for the two-body NN correlations (read from external files, see
  Alvioli, Drescher and Strikman, [Phys. Lett. B680, 225, 2009], the distributions 
  can be found at http://www.phys.psu.edu/~malvioli/eventgenerator/ )  
- The use of the Gaussian NN wounding profile (which is more realistic than 
  the commonly-used hard-core wounding profile, see the analysis by Bialas 
  and Bzdak [Acta Phys. Polon. B38, 159,2007]) 
- The generation of the core-corona distributions (see Bozek 
  [Acta Phys. Polon. B36, 3071,2005] and Werner [Phys. Rev. Lett. 98, 152301, 
  2007], see also Becattini and Manninen [Phys. Lett. B673, 19, 2009] and Bozek 
  [Phys. Rev. C79, 054901, 2009]) 
- The analysis of the triangular shape deformation parameter and profiles, 
  relevant for the triangular flow, see Alver and Roland, [Phys. Rev. C81, 054905, 
  2010] and Alver, Gombeaud, Luzum,and Ollitrault, [arXiv:1007.5469] 
- Generation of rapidity distributions in the wounded-nucleon picture according 
  to the model of Bialas and Czyz [Acta Phys.Polon.B36:905-918,2005], as implemented 
  by Bozek [arXiv:1002.4999]. This allows to obtain the fully 3-dimensional
  distribution of matter in the early Glauber phase of the collision.        
                                             
  The reference manual for ver. 2, generated by Doxygen, is supplied at the home
  page. The full write-up of ver. 2 is under preparation.                   
                                                                             
  The code can be freely used and redistributed. However, if you decide to  
  make modifications, the authors would appreciate notification for the record.             
  Any publication or display of results obtained using GLISSANDO must        
  include a reference to our papers

Computer Physics Communications 180(2009)69, arXiv:0710.5731 [nucl-th]  

and arXiv:1310.5475 [nucl-th]                             

*/

#include <math.h>
#include <time.h>
#include <string.h>
#include <iostream>
#include <fstream>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>

/* current version of the code */

#define _VER_ 2.702

Float_t ver=_VER_;  //!< current version of the code
TRandom3 raa;       //!< ROOT random number generator 

#include "counter.h"
#include "functions2.h"
#include "distrib.h"
#include "collision.h"

using namespace std;

//! the main function of GLISSANDO 2
/*!
The main function of GLISSANDO 2 contains the basic structure of the Glauber Monte Carlo simulation, i.e., 
declarations and definitions of basic objects of the nucleus and collision classes, the main loop over events, 
evaluation of basic quantities, etc. It is meant to be taylored by the user to meet his needs. 
For speed of the execution, some switches of the code are controlled by the preprocessor variables 
(_nnwp_, _files_, _profile_, _weight_, _rapidity_).
*/

Int_t main(
        Int_t argc, //!< number of command line parameters
//! used for passing input and output file names
/*! first argument: <input.dat>
    second argument: <output.root>
    third argument: <nucl_A.dat>  (only when _files_=1)
    fourth argument: <nucl_B.dat> (only when _files_=1) */
        char* argv[] // command line parameters (file names)
        ) {

/***********************
 introductory stuff
***********************/

//! (the units for all dimensionful quantities in GLISSANDO 2 are powers of fm)

Int_t ts; 
//! start time
ts=time_start(); 

//! print basic info 
helper(argc, argv[1]);

//! print header
header();

//! set the input file
// the input file can be passed as the first command-line argument
TString inpfile("input/input.dat"); 
if (argc>1) inpfile = argv[1];

//! process input parameters
readpar(inpfile);

//! set the ROOT output file
// the ROOT output file with the results can be passed as the second command-line argument
TFile *fout;
TString outfile("output/glissando.root"); 
if (argc>2) outfile = argv[2];
fout = new TFile(outfile,"RECREATE"); fout->cd();
cout << endl << "generates Root output file " << outfile;

#if(_evout_)
//! output file for hydro
ofstream eveout("output/events.dat");
cout << endl << "generates events output file output/events.dat (large output!)" << outfile;
#endif

//! seed the ROOT random-number generator
raa.SetSeed(ISEED); 
ISEED1 = ISEED;
ISEED = raa.GetSeed();

//! reset counters used for some basic physical quantities
reset_counters(); 

//! declare and initialize trees and histograms for storage of data
tr_his_c tr_his; // basic storage object for ROOT trees and histograms
tr_his.init();   // initialize

//! set the minimum wounding and binary-collision distances (hard-sphere profile) or the Gaussian wounding parameters (Gaussian profile)
d=sqrt(SNN/PI/10.);         // the minimum wounding distance, 10 converts to fm^2 from mb
Float_t d2=d*d;             // square of the minimum wounding distance
// dbin=sqrt(SBIN/PI/10.);  // the minimum binary collision distance 
dbin=d;                     // the minimum binary collision distance (note that in our model choice it is equal to the wounding distance)
Float_t dbin2=dbin*dbin;    // square of the minimum binary collision distance 
Float_t mbin=SBIN/SNN;      // ratio for the hot-spot model (if mbin<1, only a fraction of mbin binary collisions is included, but with the RDS enhanced by the factor 1/mbin)

Int_t kks=Int_t(fmax(1,fmin(1000,EVENTS/10))); // for output of the progress of the run made once a while

//! echo basic parameters to the console
echopar();

/****************************************************************************************************************/
//! #if(_files_) then initialize the nucleon distributions from external tables
/****************************************************************************************************************/

#if(_files_) // initialize tables for the nucleon distributions

FILE *nucleusA, *nucleusB;

// the input files with nuclear configurations can be passed as the third and fourth command-line arguments

if (argc<4){cout << "Provide in the command line the input files with nuclear distributions! " << endl << endl; exit(0);}

TString nAfile("nucl/pb208.dat");
if (argc>3) nAfile = argv[3];
nucleusA = fopen(nAfile,"rt");

TString nBfile("nucl/pb208.dat");
if (argc>4) nBfile = argv[4]; else nBfile = argv[3];
nucleusB = fopen(nBfile,"rt");

// count number of configurations stored in files
cout << "counting the length of nuclear distribution files " << endl;

Int_t cfc=0,bfc;
while ((bfc=fgetc(nucleusA))!=EOF) cfc+=(bfc==10)?1:0;fseek(nucleusA,0,SEEK_SET);
Int_t snuA=cfc;
if(snuA % NUMA != 0){cout << "wrong number of lines in " << nAfile << ", check! - exiting" << endl << endl; exit(0);}
nucleusA = fopen(nAfile,"rt");

cfc=0; 
while ((bfc=fgetc(nucleusB))!=EOF) cfc+=(bfc==10)?1:0;fseek(nucleusB,0,SEEK_SET);
Int_t snuB=cfc;
if((argc>4) && (snuB % NUMB != 0)){cout << "wrong number of lines in " << nBfile << ", check! - exiting" << endl << endl; exit(0);}
nucleusB = fopen(nBfile,"rt");

Float_t x, y, z;
Int_t countA=0, countB=0, charge;

Float_t* posAx=NULL; Float_t* posBx=NULL; Float_t* posAy=NULL; Float_t* posBy=NULL; Float_t* posAz=NULL; Float_t* posBz=NULL;
posAx=new Float_t[snuA]; posBx=new Float_t[snuB];
posAy=new Float_t[snuA]; posBy=new Float_t[snuB];
posAz=new Float_t[snuA]; posBz=new Float_t[snuB];

cout  << "reading " << snuA/NUMA << " configurations for nucleus A from " << nAfile << endl;
for(Int_t i=0;i<snuA;i++){
                fscanf(nucleusA,"%f %f %f %d",&x,&y,&z,&charge);
		posAx[i]=x;
	        posAy[i]=y;
	        posAz[i]=z;};

if (argc>4) {
cout  << "reading " << snuB/NUMB << " configurations for nucleus B from " << nBfile << endl << endl;
for(Int_t i=0;i<snuB;i++){
                fscanf(nucleusB,"%f %f %f %d",&x,&y,&z,&charge);
		posBx[i]=x;
	        posBy[i]=y;
	        posBz[i]=z;};
            };
#endif

/****************************************************************************************************************
  simulation
****************************************************************************************************************/

//! declare nuclei A nad B
nucleus nucA(NUMA), nucB(NUMB); // colliding nuclei

//! declare the collision
#if(_rapidity_)
collision_rap collAB(NUMA,NUMB); // collision declaration for rapidity analysis
#else
collision collAB(NUMA,NUMB); // collision declaration, no rapidity analysis
#endif

//! --- start the main loop over events
kk=0; // counter of events when a collision occurs
for(Int_t k=1;kk<EVENTS;k++){ 

evall++; // attempted A+B collision (all cases, those where an NN collision occurs or not)

//! generate the distributions of nucleons in nuclei A and B
#if(_files_)
nucA.set_file(posAx, posAy, posAz, snuA, NUMA); 
if (argc>4)  
nucB.set_file(posBx, posBy, posBz, snuB, NUMB);
else 
if(NUMB>16 && BETA2B!=0) nucB.set_random_B_def(CD); 	   
  else if(NUMB>16) nucB.set_random_B(CD); 
    else if ((NUMB<17)&&(NUMB>2)) nucB.set_random_B_hos(CD); 
      else if(NUMB==2) nucB.set_deuteron(); 
        else nucB.set_proton();
#else

if(NUMA>16 && BETA2A!=0) nucA.set_random_A_def(CD); 
  else if(NUMA>16) nucA.set_random_A(CD); 
    else if ((NUMA<17)&&(NUMA>2)) nucA.set_random_A_hos(CD); 
      else if(NUMA==2) nucA.set_deuteron(); 
        else nucA.set_proton();

if(NUMB>16 && BETA2B!=0) nucB.set_random_B_def(CD); 	   
  else if(NUMB>16) nucB.set_random_B(CD); 
    else if ((NUMB<17)&&(NUMB>2)) nucB.set_random_B_hos(CD); 
      else if(NUMB==2) nucB.set_deuteron(); 
        else nucB.set_proton();

#endif

//! shift nuclei to the center-of-mass frame
nucA.shift_cmx(); nucA.shift_cmy(); nucA.shift_cmz(); 
nucB.shift_cmx(); nucB.shift_cmy(); nucB.shift_cmz();


//! rotate deformed nuclei by theta (zx plane) and phi (xy plane) angles 
//! The proper order of rotations is important. Opposite one (phi,theta) has no sense.
 
if(ROTA_THETA==-1.0)
{
   Float_t ethnucl=0.0;
   ethnucl=2*los()-1;
   nucA.rotate_polar(ethnucl);
}
   else
   {Float_t ethnucl=0.0;
    ethnucl=cos(PI*(ROTA_THETA)/180.0);
    nucA.rotate_polar(ethnucl);
   }

if(ROTA_PHI==-1.0)
{
   Float_t phnucl=0.0;
   phnucl=2*PI*los();
   nucA.rotate(phnucl);
}
   else
   {Float_t phnucl=0.0;
    phnucl=2*PI*(ROTA_PHI)/360.0;
    nucA.rotate(phnucl);
   }

if(ROTB_THETA==-1.0)
{   
   Float_t ethnucl=0.0;                              

   ethnucl=2*los()-1;
   nucB.rotate_polar(ethnucl);
}
   else
   {Float_t ethnucl=0.0;
    ethnucl=cos(PI*(ROTB_THETA)/180.0);
    nucB.rotate_polar(ethnucl);
   }

if(ROTB_PHI==-1.0)
{
   Float_t phnucl=0.0;              

   phnucl=2*PI*los();
   nucB.rotate(phnucl);
}
   else
   {Float_t phnucl=0.0;
    phnucl=2*PI*(ROTB_PHI)/360.0;
    nucB.rotate(phnucl);
   }

/*! #(if_profile_) generate the histograms of one-body density in x-y and r-cos(theta) coordinate systems
for nucleus A, and for the relative NN distance */ 
#if(_profile_)
Float_t r_nucl;  /*!< the radial position of the nucleon in the nucleus */
for(Int_t j=0;j<nucA.n;j++){

//! calculation of the radial position of the nucleon in the nucleus A
 r_nucl=sqrt(nucA.x[j]*nucA.x[j]+nucA.y[j]*nucA.y[j]+nucA.z[j]*nucA.z[j]);  
 tr_his.radA->Fill(sqrt(nucA.x[j]*nucA.x[j]+nucA.y[j]*nucA.y[j]+nucA.z[j]*nucA.z[j])); 
  
//! calculation of one-body density in x-y coordinate system in nucleus A
 tr_his.xyhist_nuclA->Fill(nucA.x[j]*1.0,nucA.y[j]*1.0);         
//! calculation of one-body density in r-cos(theta) coordinate system in nucleus A
 tr_his.rcostheta_nuclA->Fill(r_nucl*1.0,nucA.z[j]*1.0/r_nucl); 

    
          for(Int_t i=0;i<j;i++){
                           tr_his.rrelA->Fill(sqrt((nucA.x[j]-nucA.x[i])*(nucA.x[j]-nucA.x[i])+
                                                  (nucA.y[j]-nucA.y[i])*(nucA.y[j]-nucA.y[i])+
                                                  (nucA.z[j]-nucA.z[i])*(nucA.z[j]-nucA.z[i])));  
                      if(NUMA==NUMB){              
                         tr_his.rrel_u->Fill(sqrt((nucB.x[j]-nucA.x[i])*(nucB.x[j]-nucA.x[i])+
                                                  (nucB.y[j]-nucA.y[i])*(nucB.y[j]-nucA.y[i])+
                                                  (nucB.z[j]-nucA.z[i])*(nucB.z[j]-nucA.z[i])));  
                                     };
                               };
                         };
for(Int_t j=0;j<nucB.n;j++){

//! calculation of the radial position of the nucleon in the nucleus B
 r_nucl=sqrt(nucB.x[j]*nucB.x[j]+nucB.y[j]*nucB.y[j]+nucB.z[j]*nucB.z[j]);  
 tr_his.radB->Fill(sqrt(nucB.x[j]*nucB.x[j]+nucB.y[j]*nucB.y[j]+nucB.z[j]*nucB.z[j]));  

//! calculation of one-body density in x-y coordinate system in nucleus B
 tr_his.xyhist_nuclB->Fill(nucB.x[j]*1.0,nucB.y[j]*1.0);                    
//! calculation of one-body density in r-cos(theta) coordinate system in nucleus B
 tr_his.rcostheta_nuclB->Fill(r_nucl*1.0,nucB.z[j]*1.0/r_nucl);            

    
          for(Int_t i=0;i<j;i++){
                           tr_his.rrelB->Fill(sqrt((nucB.x[j]-nucB.x[i])*(nucB.x[j]-nucB.x[i])+
                                                  (nucB.y[j]-nucB.y[i])*(nucB.y[j]-nucB.y[i])+
                                                  (nucB.z[j]-nucB.z[i])*(nucB.z[j]-nucB.z[i])));  
                              };
                         };						 
#endif

//! generate the impact parameter b with the distribution proportional to b^2 in the range (BMAX, BMIN)
b=sqrt((BMAX*BMAX-BMIN*BMIN)*los()+BMIN*BMIN);

//! shift the coordinates of the nucleons in nucleus A such that the center of mass is at the point (b*NUMB/(NUMA+NUMB),0)
nucA.shift_x(b*NUMB/(NUMA+NUMB));

//! shift the coordinates of the nucleons in nucleus B such that the center of mass is at the point (-b*NUMA/(NUMA+NUMB),0)
nucB.shift_x(-b*NUMA/(NUMA+NUMB));

//! collide the nuclei, create the sources (wounded nucleons, binary collisions) and RDS
collAB.gen_RDS(nucA,nucB,d2,dbin2,mbin);

// Include the event only when the number of wounded nucleons is at least WMIN (=2 by default)
// and, additionally, lies between W0 and W1 and the relative deposited strength (RDS) lies between RDS0 and RDS1.
// We also request that the number of sources with non-zero RDS >= WMIN (relevant
// for the Poisson superposition, where RDS can be zero)

if((collAB.nwAB >= WMIN) && (collAB.nwAB >= W0) && (collAB.nwAB <= W1) && (collAB.rpa >= RDS0) && (collAB.rpa <= RDS1) && (collAB.nzw >= WMIN)){
	
kk++; // count the event in the specified window
	
// for asymmetric collisions necessarily use SHIFT=1 in the input
if(SHIFT==1){collAB.shift_cmx_w(); collAB.shift_cmy_w();};

//! #if(_rapidity_) generate the rapidity distribution

#if(_rapidity_)
collAB.gen_rap(NUMRAP,MAXYRAP); 

Float_t phref=collAB.phrot(ARANK);
collAB.rotate(phref);

collision_rap collAB_plus=collAB; 

collAB_plus.shift_rap(FBRAP);
collAB_plus.shift_cmx_w(); collAB_plus.shift_cmy_w();
Float_t pht2_plus=collAB_plus.phrot(ARANK); 
phirot_plus=pht2_plus;

collision_rap collAB_minus=collAB;
collAB_minus.shift_rap(-FBRAP);

collAB_minus.shift_cmx_w(); collAB_minus.shift_cmy_w();
Float_t pht2_minus=collAB_minus.phrot(ARANK);  
phirot_minus=pht2_minus;

angles.add(sin(ARANK*phirot_plus),sin(ARANK*phirot_minus)); // correlation between sines of forward and backward angles

collAB.rotate(-phref); // rotate back
#endif

//! #if(_weight_) generate the histograms for the NN collision profiles
#if(_weight_)
for(Int_t j=0;j<collAB.n;j++){if(collAB.c[j]!=0){tr_his.weih->Fill(collAB.w[j],1);}
                            else {tr_his.weih_bin->Fill(collAB.w[j],1);};};
#endif

//! generate various 2-dim histograms with the distributions of sources
collAB.fill_xy(tr_his.xyhist, 1./(4*BTOT*BTOT/NBIN/NBIN)/EVENTS); // all sources

collAB.fill_polar(tr_his.c0hist, 0, 1./(2*PI)/(BTOT/NBIN)/EVENTS);
collAB.fill_polar(tr_his.c2hist, 2, 1./(2*PI)/(BTOT/NBIN)/EVENTS);
collAB.fill_polar(tr_his.c3hist, 3, 1./(2*PI)/(BTOT/NBIN)/EVENTS);
collAB.fill_polar(tr_his.c4hist, 4, 1./(2*PI)/(BTOT/NBIN)/EVENTS);
collAB.fill_polar(tr_his.c5hist, 5, 1./(2*PI)/(BTOT/NBIN)/EVENTS);
collAB.fill_polar(tr_his.c6hist, 6, 1./(2*PI)/(BTOT/NBIN)/EVENTS);

collAB.fill_polar_s(tr_his.s3hist, 3, 1./(2*PI)/(BTOT/NBIN)/EVENTS);

//! generate the fixed-axes Fourier moments (obsolete)
//! weight r^2
es=collAB.eps(2,2); es3=collAB.eps(3,2); es4=collAB.eps(4,2); es5=collAB.eps(5,2); es6=collAB.eps(6,2);
ess=collAB.eps_s(2,2); es3s=collAB.eps_s(3,2); es4s=collAB.eps_s(4,2); es5s=collAB.eps_s(5,2); es6s=collAB.eps_s(6,2); 


#if(_evout_)
collAB.writerds(eveout);
#endif

//! generate the variable-axes Fourier moments (up to 6-th moment)

// rotate to the frame maximizing the m=1 moment a la Derek Teaney

Float_t pht1=collAB.phrot(1,3); 
phirot=pht1; 
collAB.rotate(pht1);
ep1=collAB.eps(1,3); ep1s=collAB.eps_s(1,3);
collAB.rotate(-pht1);

// rotate to the frame maximizing the m=2 moment 
roo=2;
if(RO==0){roo=2;};
ppp=PP; if(PP==-1) ppp=2;
Float_t pht2=collAB.phrot(roo,ppp); 
phirot=pht2; 
collAB.rotate(pht2);
collAB.fill_xy(tr_his.xyhistr, 1./(4*BTOT*BTOT/NBIN/NBIN)/EVENTS);
collAB.fill_polar(tr_his.c0rhist, 0, 1./(2*PI)/(BTOT/NBIN)/EVENTS);
collAB.fill_polar(tr_his.c2rhist, 2, 1./(2*PI)/(BTOT/NBIN)/EVENTS);
collAB.fill_xy(tr_his.xyhist_mantle, 1./(4*BTOT*BTOT/NBIN/NBIN)/EVENTS, 1, 1);  // corona sources
collAB.fill_xy(tr_his.xyhist_core, 1./(4*BTOT*BTOT/NBIN/NBIN)/EVENTS, 2, 1000); // core sources
ep=collAB.eps(2,ppp); eps=collAB.eps_s(2,ppp);
collAB.rotate(-pht2);

// rotate to the frame maximizing the m=3 moment 
roo=2;
if(RO==0){roo=3;};
ppp=PP; if(PP==-1) ppp=3;
Float_t pht3=collAB.phrot(roo,ppp);  
phirot3=pht3;
collAB.rotate(pht3);
collAB.fill_polar(tr_his.c3rhist, 3, 1./(2*PI)/(BTOT/NBIN)/EVENTS);
collAB.fill_polar_s(tr_his.s3rhist, 3, 1./(2*PI)/(BTOT/NBIN)/EVENTS);
ep3=collAB.eps(3,ppp); ep3s=collAB.eps_s(3,ppp); 
collAB.rotate(-pht3);

// rotate to the frame maximizing the m=4 moment 
roo=2;
if(RO==0){roo=4;};
ppp=PP; if(PP==-1) ppp=4;
Float_t pht4=collAB.phrot(roo,ppp);  
phirot4=pht4;
collAB.rotate(pht4);
collAB.fill_polar(tr_his.c4rhist, 4, 1./(2*PI)/(BTOT/NBIN)/EVENTS); 
ep4=collAB.eps(4,ppp); ep4s=collAB.eps_s(4,ppp); 
collAB.rotate(-pht4);

// rotate to the frame maximizing the m=5 moment 
roo=2;
if(RO==0){roo=5;};
ppp=PP; if(PP==-1) ppp=5;
Float_t pht5=collAB.phrot(roo,ppp); 
phirot5=pht5;
collAB.rotate(pht5);
collAB.fill_polar(tr_his.c5rhist, 5, 1./(2*PI)/(BTOT/NBIN)/EVENTS);
ep5=collAB.eps(5,ppp); ep5s=collAB.eps_s(5,ppp); 
collAB.rotate(-pht5);

// rotate to the frame maximizing the m=6 moment 
roo=2;
if(RO==0){roo=6;};
ppp=PP; if(PP==-1) ppp=6;
Float_t pht6=collAB.phrot(roo,ppp);  
phirot6=pht6;
collAB.rotate(pht6);
collAB.fill_polar(tr_his.c6rhist, 6, 1./(2*PI)/(BTOT/NBIN)/EVENTS); 
ep6=collAB.eps(6,ppp); ep6s=collAB.eps_s(6,ppp); 
collAB.rotate(-pht6);

//! get some basic properties of the event
rwA=collAB.nwA; rwB=collAB.nwB; rwAB=collAB.nwAB; rbin=collAB.nbin; rhotspot=collAB.nhotspot; rpa=collAB.rpa;
sizeav=collAB.size(); xx=collAB.cmx_w();  yy=collAB.cmy_w(); 

nwounded.add(rwAB); nbinary.add(rbin); nhot.add(rhotspot); nweight.add(rpa);
estd.add(es); epart1.add(ep1); epart.add(ep); estd3.add(es3); epart3.add(ep3); estd4.add(es4); epart4.add(ep4);

//! fill the data in trees and histograms
	tr_his.fill_res();
        tr_his.fill_tr();

//! if(FULL) write the full event info to the file (added for comparability reasons, takes a lot of space)
if(FULL){for(Int_t j=0;j<collAB.n;j++){
			tSource.X=collAB.x[j]; 
			tSource.Y=collAB.y[j]; 
			tSource.W=collAB.w[j]; 
			tSource.KK=kk; 
			tr_his.full_event->Fill();};};

}; // if((collAB.nwAB >= WMIN)...

// write the current event count
if((kk % kks)==0){cout << "\revent: " << kk << "   ("<< Int_t(Float_t(kk)/Float_t(EVENTS)*100.) << "%)              "; cout.flush();};

}; // end of loop over events 

//! --- end of main loop over events

cout << "\revent: " << kk << "     ("<< 100 << "%)              "; cout.flush();
cout << endl << endl;

/****************************************************************************************************************/
//! output of results
/****************************************************************************************************************/

//! the total cross section and the equivalent hard-sphere radius 

// the total cross section and the equivalent hard-sphere radius
sitot=PI*BMAX*BMAX*EVENTS/evall*10; // obtained NN cross section, 10 converts from fm^2 to mb
sirad=sqrt(sitot/PI/10.)/2.; // equivalent hard-sphere NN interaction radius in fm    

//! project out the marginal distribution in the radial variable (generate the radial Fourier profiles)
tr_his.proj();
tr_his.write_d();

//! generate and write some histograms with physical quantities
tr_his.gen();

// compute some physical info
xeps=estd.mean(); xseps=sqrt(estd.var()); xepp=epart.mean(); xsepp=sqrt(epart.var());

// write out the histograms 
        tr_his.fill();
        tr_his.write();
    #if(_profile_)
        tr_his.write_r();
    #endif

//! #if(_weight_) normalize to the wounding and the binary cross sections and write
 
    #if(_weight_)
        tr_his.write_w();

// normalize to the wounding cross section         
        Float_t intn=collAB.w_distr->Integral();
        collAB.w_distr->Scale(1./intn/PI*SNN*10/15);
        collAB.w_distr->Write();

// normalize to the binary cross section         
        intn=collAB.w_distr_bin->Integral();
        collAB.w_distr_bin->Scale(1./intn/PI*SBIN*10/15);
        collAB.w_distr_bin->Write();
    #endif

//! #if(_rapidity_) write rapidity distribution

    #if(_rapidity_)        
        collAB.rap_distr->Write();
    #endif

//! closing ROOT file
fout->Close();
#if(_evout_)
eveout.close();
#endif

//! write exit info
epilog();

//! stop time and print stamp
time_stop(ts);

}

