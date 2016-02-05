/** \file retrieve.cxx
 * auxilliary file, part of GLISSANDO 2
 * 
 */

#include <math.h>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

#include <fstream>
#include <iostream>
#include <sstream>


using namespace std;

//! structure for output of the full event - transverse coordinates, weight, number of the event
typedef struct {
                Float_t X, //!< x coordinate
                        Y, //!< y coordinate
                        W; //!< z coordinate
                UInt_t KK; //!< number of the event
                } SOURCE;


/*! The code serves as an example of using the full spatial distribution of sources, which 
    optionally may be generated in the simulation. */
Int_t main(
   Int_t argc, //!< number of command line parameters
   char **argv //!< name of the GLISSANDO output ROOT file        
          ){

// Default file name
TString inpfile("../output/glissando.root");
if (argc>1) inpfile = argv[1];
TFile f(inpfile);


  SOURCE buf;

  Int_t numev=0;            // number of event (all events)
  Int_t curev=1;            // current event + 1 (in each file)
  Int_t sw;                 // =1 when the event is changed, 0 otherwise
  Int_t ilow, ihigh;        // range of entries in the event

  // kinematic variables for particle 1 and 2
  Float_t x, y, w, b, xav, xxav=0, xx2av=0;
  Int_t iter, kk, m;
  Float_t mav=0, m2av=0;


// for a chain of files
  TChain *itree = new TChain("full_event");
  TChain *etree = new TChain("events");
  itree->Add(inpfile);
  // itree->Add("glissando_2.root"); // adding more files with data
  etree->Add(inpfile);
  // etree->Add("glissando_2.root");

  itree->SetBranchAddress("full_source", &buf); // set the branch to the sources
  etree->SetBranchAddress("b",&b);              // set the branch to the impact parameter b

  cout << "got data" <<endl;
  cout << "will print sample info only for first 10 and last 10 events" << endl;

// go over all entries
  ihigh=0;
  // the trick is that when the buf.KK, describing the number of event, changes, we know we pass to the new event
  for (iter=0; iter < itree->GetEntries(); iter++)
  { ilow=ihigh;
    m=0; xav=0;
	itree->GetEntry(iter);
    sw=0;
    if (curev != buf.KK){curev = buf.KK;sw=1;numev++;ihigh=iter;};

// do at the end of each event
  if(iter == (itree->GetEntries())-1 ){numev++;ihigh=iter+1;}; // in the last entry
  if(sw || iter == (itree->GetEntries())-1 ){
    etree->GetEntry(numev-1); // get the impact parameter
    if(numev<=10 || numev>=(etree->GetEntries())-10){ // test output for the first and last 10 events
  cout << "event: " << numev << "  range for particle labels: " << ilow << " - " << ihigh-1 << "  b= " << b << endl; 
  cout.flush();
  };

// loop over sources
for(Int_t j1=ilow;j1<ihigh;j1++){
  itree->GetEntry(j1); // here we have access to the sources
 m++; xav=xav+buf.X;
  } // end "loop over sources"

// a sample use of data
  xav=xav/m;
  mav=mav+m;
  m2av=m2av+m*m;
  xxav=xxav+xav;
  xx2av=xx2av+xav*xav;

  } //end "do at the end of each event"


 } // end "go over all entries"

cout << endl;
cout <<"# of sources=" << Float_t(mav)/numev << "+/-" << sqrt((Float_t)m2av/(numev-1)-(Float_t)mav/numev*mav/(numev-1)) << endl;
cout <<"x c.m. coordinate=" << xxav/numev << "+/-" << sqrt(xx2av/(numev-1)-xxav/numev*xxav/(numev-1)) << endl;
} // end main


