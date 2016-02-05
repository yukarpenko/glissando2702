/** \file interpolation.cxx
 * auxilliary file, part of GLISSANDO 2
 * 
 */

#include <math.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>

#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>


using namespace std;

//! name of the histogram from the GLISSANDO output ROOT file
char hname[100];

//! 1D Lagrange interpolation, calculation of value of the Lagrange polynomial f(x) at given x
Double_t w1D(Int_t n, Double_t dx, Double_t x[4], Double_t y[4]){
Int_t i,j;
Double_t dy,m,l;

dy=0;
for (i=0;i<n; i++) {
 l=1; m=1;
  for (j=0;j<n;j++) 
    if (i!=j) {
     m=m*(x[i]-x[j]);
     l=l*(dx-x[j]);
     }
     dy+=y[i]*l/m;
 }
 return dy;
}

//! 2D Lagrange interpolation, calculation of value of Lagrange polynomial (first order)
Double_t w2D(Double_t x,Double_t y,Double_t x0,Double_t y0,Double_t u0,Double_t x1,Double_t y1,Double_t u1,Double_t x2,Double_t y2,Double_t u2){

Double_t m = (x1*y2-x2*y1)-(x0*y2-x2*y0)+(x0*y1-x1*y0);
Double_t w0=u0*((x1*y2-x2*y1+(y1-y2)*x+(x2-x1)*y)/m);
Double_t w1=u1*((x2*y0-x0*y2+(y2-y0)*x+(x0-x2)*y)/m);
Double_t w2=u2*((x0*y1-x1*y0+(y0-y1)*x+(x1-x0)*y)/m);

return w0+w1+w2;
}


//! Two nodes. Linear regression method. Calculation of f(x)=a+x*b
Double_t lin(Double_t dx,Double_t x1,Double_t y1,Double_t x2,Double_t y2){
Double_t al=2.0*(x1*y1+x2*y2)-((x1+x2)*(y1+y2));
Double_t am=2.0*(x1*x1+x2*x2)-((x1+x2)*(x1+x2));
Double_t a=al/am;
Double_t b=(y1+y2-a*(x1+x2))/2.0;
return a*dx+b;
}

//! 1D interpolation, determination of the node values from the given histogram.
Double_t inter1D(Double_t dx,TFile & f,char *histname) {

Int_t i,j,ibin;
Double_t x[500],y[500];
Double_t xw[4],yw[4];

for (i=0;i<500;i++){
x[i]=0;
y[i]=0;
}

for (i=0;i<4;i++){
xw[i]=0;
yw[i]=0;
}

TH1F *h = (TH1F*)f.Get(histname); 
TAxis *xaxis = h->GetXaxis();

Int_t nr = xaxis->GetNbins();

for (i=1;i<=nr;i++) {
Double_t s = xaxis->GetBinCenter(i);
x[i-1]=s;
Double_t z = h->GetBinContent(i);
y[i-1]=z;
}

if ((dx<0)||(dx>x[nr-2])){
return -1.;
}

if (dx<=x[1]){ 
Double_t x1=x[0];
Double_t y1=y[0];
Double_t x2=x[1];
Double_t y2=y[1];
return lin(dx,x1,y1,x2,y2);
}

for (i=1;i<=nr-1;i++)
if ((dx>x[i]) && (dx<=x[i+1])){
  ibin=i;
  break;
}

 j=0;
for (i=ibin-1;i<=ibin+2;i++) {
 xw[j]=x[i];
 yw[j]=y[i];
 j++;
 }
return w1D(4,dx,xw,yw);
}

//! 2D interpolation, determination of the node values from the given histogram.
Double_t inter2D(Double_t dx, Double_t dy, TFile & f,char *histname) {

Double_t x[500],y[500],z[500][500];

Int_t i,j,ibin,jbin;

for (i=0;i<500;i++){
x[i]=0;
y[i]=0;
}

for (i=0;i<500;i++){
for (j=0;j<500;j++){
z[i][j]=0;
}
}

TH1F *h = (TH1F*)f.Get(histname); 

TAxis *xaxis = h->GetXaxis();
TAxis *yaxis = h->GetYaxis();

Int_t nrx = xaxis->GetNbins();
Int_t nry = yaxis->GetNbins();

for (i=1;i<=nrx;i++) {
Double_t sx = xaxis->GetBinCenter(i);
x[i-1]=sx;
}
for (i=1;i<=nry;i++){
Double_t sy = yaxis->GetBinCenter(i);
y[i-1]=sy;
}
for (i=1;i<=nrx;i++){
for (j=1;j<=nry;j++){
Double_t cont = h->GetBinContent(i,j);
z[i-1][j-1]=cont;
}
}

if ((dx<x[0])||(dx>x[nrx-1])){
return -1.;
}
if ((dy<y[0])||(dy>y[nrx-1])){
return -1.;
}

for (i=0;i<=nrx-1;i++){
for (j=0;j<nry-1;j++){
if ((dx>=x[i]) && (dx<x[i+1])){
if ((dy>=y[j]) && (dy<y[j+1])){
  ibin=i;
  jbin=j;
  break;
}
}
}
}

Double_t x1=x[ibin];
Double_t y1=y[jbin];
Double_t x2=x[ibin+1];
Double_t y2=y[jbin+1];

Double_t fun=lin(dx,x1,y1,x2,y2);

if (fun>=dy){ 
//cout << " bottom " << endl;
return w2D(dx,dy,x1,y1,z[ibin][jbin],x2,y1,z[ibin+1][jbin],x2,y2,z[ibin+1][jbin+1]);
}
else{
//cout << " top " << endl;
return w2D(dx,dy,x1,y1,z[ibin][jbin],x2,y2,z[ibin+1][jbin+1],x1,y2,z[ibin][jbin+1]);
}

}

//! 1D interpolation 
void d1(TFile & f){
Double_t x;
system("clear");
cout << "1D Interpolation" <<endl;
cout << endl;
cout << "Give name of 1D histogram (name from the GLISSANDO output file, e.g., c0rhp, c2rhp, etc. " <<endl;
cin>> hname;
cout << "Give value of rho: ";
cin >> x;
Double_t w = inter1D(x,f,hname);
if (w<=-1.) {
 cout <<"Out of range"<<endl;
 } else {
cout << "Interpolated value at " << x << " equals " << w <<endl;
}
}

//! 2D interpolation
void d2(TFile & f){
Double_t x,y;
system("clear");
cout << "2D Interpolation" <<endl;
cout << endl;
cout << "Give name of 2D histogram (name from the GLISSANDO output file, e.g., xyhistr" <<endl;
cin>> hname;
cout << "Give value of rho_x: ";
 cin >> x;
cout << "and"<<endl;
cout << " rho_y: ";
cin >> y;
Double_t w = inter2D(x,y,f,hname);
if (w<=-1.) {
 cout <<"Out of range"<<endl;
 } else {
cout << "Interpolated value at " << x <<", "<< y << " equals " <<w <<endl;
}
}

//! Start menu
void start(){
cout << "INTERPOLATION"<< endl;
cout << "ver 1.0" << endl;
cout << endl;

cout << "1 - One Dimensional Interpolation" << endl;
cout << "2 - Two Dimensional Interpolation" << endl;
cout << endl;
cout << "Type 1 or 2 (any other - Exit)"<< endl;
}

//! Determination of value of mean weight from the npa branch stored in GLISSANDO Root file
Double_t mean_weight(TString & inpfile){

Float_t npa;

// TTree *ev = (TTree*)f.Get("events");

// for a chain of files
  TChain *ev = new TChain("events");
  ev->Add(inpfile);

ev->SetBranchAddress("npa",&npa);

Double_t max=0;
Int_t nentries = (Int_t)ev->GetEntries();
for (Int_t i=0; i<nentries; i++) {
ev->GetEntry(i);
if (npa>max) max=npa;
}

TH1F *hnpa = new TH1F("hnpa","weight distribution",(Int_t)max,0,max);

for (Int_t i=0; i<nentries; i++) {
ev->GetEntry(i);
hnpa->Fill(npa);
}
Double_t sw=hnpa->GetMean();
return sw;
}

//! Calculation the normalization integral from the given 1-dim histogram stored in the GLISSANDO Root file. Rectangle method.
Double_t integral1D(Double_t nbc,TFile & f, char *histname){

TH1F *h = (TH1F*)f.Get(histname);
TAxis *xaxis = h->GetXaxis();
Int_t nr = xaxis->GetNbins();
Double_t width=xaxis->GetBinWidth(1);
Double_t start=xaxis->GetBinCenter(1);
Double_t end=xaxis->GetBinCenter(nr-2);

Double_t max=width*nr;
Double_t w=max/nbc;
Double_t xcalk=start;
Double_t c=0;
while (xcalk<end){
 c=c+(xcalk*w*inter1D(xcalk,f,histname));
 xcalk+=w;
 }
 c=2.0*3.1415926*c;

return c;
}

//! Calculation of value of integral: f(rho_x,rho_y) drho_x drho_y from the given histogram stored in the GLISSANDO Root file. Rectangle method.
Double_t integral2D(Double_t nbc, TFile & f, char *histname){

TH1F *h = (TH1F*)f.Get(histname);
TAxis *xaxis = h->GetXaxis();
TAxis *yaxis = h->GetYaxis();
Int_t nrx = xaxis->GetNbins();
Int_t nry = yaxis->GetNbins();

Double_t widthx=xaxis->GetBinWidth(1);
Double_t widthy=yaxis->GetBinWidth(1);

Double_t startx=xaxis->GetBinCenter(2);
Double_t endx=xaxis->GetBinCenter(nrx-1);
Double_t starty=yaxis->GetBinCenter(2);
Double_t endy=yaxis->GetBinCenter(nry-1);

Double_t maxx=widthx*nrx;
Double_t maxy=widthy*nry;

Double_t wx=maxx/nbc;
Double_t wy=maxy/nbc;

Double_t xcalk=startx;
Double_t ycalk=starty;

Double_t c1=0;
while (ycalk<endy){
 xcalk=startx;
 while (xcalk<endx) {
 c1=c1+(wx*wy*inter2D(xcalk,ycalk,f,histname));
 xcalk+=wx;
 }
 ycalk+=wy;
 }

return c1;
}

/*! The code computes the values of the one-dimensional or two-dimansional profiles at a specified point via
    Lagrange interpolation. Useful, if exporting of these data for other applications is needed. */
Int_t main(
        Int_t argc, //!< number of command line parameters
        char **argv //!< name of the GLISSANDO output ROOT file
       ){

Double_t nbc=100;
Double_t nbc2D=50;

// Default file name
TString inpfile("../output/glissando.root");
if (argc>1) inpfile = argv[1];
TFile f(inpfile);

char t, more='y';

cout << "Reading from " << inpfile << endl;

cout << "Test of normalization - wait a moment ..."<<endl;
cout <<"Mean RDS= "<< mean_weight(inpfile)<<endl;
cout <<"Integral over 2 Pi rho f0(rho) drho= "<<integral1D(nbc,f,"c0hp")<<endl;
cout <<"Integral over (xyhist) f(rho_x,rho_y) drho_x drho_y= "<<integral2D(nbc2D,f,"xyhist")<<endl;

start();
cin>>t;
if (t=='1'){
while ((more=='y')||(more=='Y')) {
d1(f);
cout << "Once more?(y/n)"<<endl;
cin>>more;
}
} else 
if (t=='2'){
while ((more=='y')||(more=='Y')){
d2(f);
cout << "Once more?(y/n)"<<endl;
cin>>more;
}
}
else {
cout  <<endl;
cout <<"You put number different than \"1\" or \"2\" - exitting"<<endl;
}
f.Close("R");
return 0;
}
