/** \file profile2_deformation_U.C
 * Script generating r-cos(theta) distributions of the nucleons in the nuleus 
 * (part of GLISSANDO 2)
 * 
 */


#include "label.C"

//! produces plots of the fixed-axes r-cos(theta) distributions of the nucleons in the nuleus.
/*! The plots are obtained as 2D r-cos(theta) histograms of nucleons positions*/
void profile2_deformation_U(
char *p, //!< name of the GLISSANDO input ROOT for the deformed case
char *ps //!< name of the GLISSANDO input ROOT for the spherical case
             ){

gROOT->Reset();

TString empty("");

// Default file name
TString inpfile("glissando.root");
 if (p!=empty) inpfile = p;
 cout << "reads from: " << inpfile << endl;;  

TFile *f = new TFile(inpfile);  

TH2D *r_costheta_hist = (TH2D*)f->Get("rcostheta_nuclA");

r_costheta_hist->SetTitle("r-cos(#theta) profile, nucleus U (deformed, #beta2 = 0.28, #beta4 = 0.093)");
r_costheta_hist->SetStats(kFALSE);

Float_t label_size=0.04;
Float_t title_size=0.05;
Float_t title_offset=1.5;

char y_title[]="cos(#theta)  ";
char x_title[]="r [fm]  ";


r_costheta_hist->GetXaxis()->SetLabelSize(label_size);
r_costheta_hist->GetXaxis()->SetTitleSize(title_size);
r_costheta_hist->GetXaxis()->SetTitleOffset(title_offset);
r_costheta_hist->GetXaxis()->SetTitle(x_title);

r_costheta_hist->GetYaxis()->SetLabelSize(label_size);
r_costheta_hist->GetYaxis()->SetTitleSize(title_size);
r_costheta_hist->GetYaxis()->SetTitleOffset(title_offset);
r_costheta_hist->GetYaxis()->SetTitle(y_title);

// Default file2 name
TString inpfile2("glissando.root");
 if (ps!=empty) inpfile2 = ps;
 cout << "reads from: " << inpfile2 << endl;;

TFile *f2 = new TFile(inpfile2);

TH2D *r_costheta_hist_sym = (TH2D*)f2->Get("rcostheta_nuclA");

r_costheta_hist_sym->SetTitle("r-cos(#theta) profile, nucleus U (spherical)");
r_costheta_hist_sym->SetStats(kFALSE);

r_costheta_hist_sym->GetXaxis()->SetLabelSize(label_size);
r_costheta_hist_sym->GetXaxis()->SetTitleSize(title_size);
r_costheta_hist_sym->GetXaxis()->SetTitleOffset(title_offset);
r_costheta_hist_sym->GetXaxis()->SetTitle(x_title);

r_costheta_hist_sym->GetYaxis()->SetLabelSize(label_size);
r_costheta_hist_sym->GetYaxis()->SetTitleSize(title_size);
r_costheta_hist_sym->GetYaxis()->SetTitleOffset(title_offset);
r_costheta_hist_sym->GetYaxis()->SetTitle(y_title);

gStyle->SetOptStat(0);
gStyle->SetStatBorderSize(0);
 
  TCanvas *c1= new TCanvas("c1", "r-cos(theta) distribution",1200, 600);
  c1->SetFillColor(0);
  c1->Divide(2,1);
  c1->cd(1);

label(inpfile);
 
r_costheta_hist->Draw("SURF");

  c1->cd(2);

label(inpfile2);

r_costheta_hist_sym->Draw("SURF");

c1->Modified();
c1->Update();
c1->SaveAs("rcostheta_spherical_deformed_U.eps");

}
