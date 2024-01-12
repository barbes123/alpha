//------------------------------------------------------------------
//(c)dUmon v.31OCT2017 Legnaro, fit 1-5, isotope neutron
//(c)dUmon v.17MARS13.2 Orsay, fit 1,2,isotope neutron
//------------------------------------------------------------------
#include <iostream>
#include <iomanip>
#include "TString.h"
#include <string.h>    
#include <fstream>
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <ctime>
#include "TStyle.h"
#include "TH1F.h"
#include "TF1.h"
#include "TFormula.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMath.h"
#include <TObjArray.h>
#include <TObjString.h>
#include <string>
#include <fstream>
#include <sstream>  
#include <TPaveText.h>

// #include "TROOT.h"
// #include "TChain.h"
// #include "TFile.h"
// #include <iostream>
// #include <string>
// #include <fstream>
// #include <stdlib.h> 
// #include <TSystem.h>
// #include <vector>
// #include <utility>
// #include <ctime>
// #include "TTree.h"
// #include "TF1.h"
// #include "TMath.h"
//#include <../physics/root/root_v5.30.05/include/TString.h>

using namespace std ;
TObjArray *toks;
bool debug = false;
bool blRebinPhiBG0 = false;

const double log_2 = 0.69314718;
const int increment_shift = 0;

double Tbeam, Tdecay, Tbg;
double TbeamPrime;

const int plot_option = 1; //1 - graph;0 - hist

double fit_lifetime2(Double_t *x, Double_t *par) 
{   
    //par[0]=Phi1; par[1]=lambda1; par[2]=Pn1
    //par[3]=Phi2; par[4]=lambda2; par[5]=Pn2
    //par[6]=NeutBG
 
    Double_t first_b, second_b, first_dec, second_dec, third_b, third_dec, fouth_b, fouth_dec, fith_b, fith_dec;
    Double_t  fit_5iso;
    first_b =  par[2]*(par[0] - par[0]/(exp((x[0] - Tbg)*par[1])));

    second_b = par[5]*(-par[2]+1)*(exp((x[0] - Tbg)*par[4])*par[4]*par[0] + exp((x[0] - Tbg)*(par[1] + par[4]))*(par[1] - par[4])*(par[0] + par[3]) +    
              exp((x[0] - Tbg)*par[1])*(par[4]*par[3] - par[1]*(par[0] + par[3])))/(exp((x[0] - Tbg)*(par[1] + par[4]))*(par[1] - par[4]));
	      
    first_dec = par[2]*(exp((-1)*(par[5] + ((x[0] - Tbg) - par[5]))* par[1])* (-1 + exp(par[5]* par[1]))* par[0]);

    second_dec = par[5]*(-par[2]+1)*(exp((-1)*(par[5] + ((x[0] - Tbg) - par[5]))* (par[1] + par[4]))* (exp((par[5] + ((x[0] - Tbg) - par[5]))* par[4])*par[4]*par[0] - exp(((x[0] - Tbg) - par[5])* par[4] +  par[5]* (par[1] + par[4]))*par[4]* par[0] +
              exp((par[5] + ((x[0] - Tbg) - par[5]))* par[1])*(par[4]* par[3] - par[1] *(par[0] + par[3])) + exp(((x[0] - Tbg) - par[5])* par[1] +  par[5]* (par[1] + par[4]))*
              ((-1)*par[4]* par[3] + par[1]*(par[0] + par[3]))))/(par[1] - par[4]);
	 TbeamPrime = par[5]+Tbg;
    return  (x[0]<=Tbg)*(par[3]+x[0]*par[4])+((x[0]>Tbg)&&(x[0]<=TbeamPrime))*(first_b + par[3]+x[0]*par[4]) + (x[0]>TbeamPrime)*(first_dec + par[3]+x[0]*par[4]);  
//      if (n==2)return  (x[0]<=Tbg)*(par[6])+((x[0]>Tbg)&&(x[0]<=TbeamPrime))*(first_b + second_b + par[6]) + (x[0]>TbeamPrime)*(first_dec + second_dec + par[6]);
      
};



void fitneutron(TString szInFile="", TString szParInput = "paramn.txt", string szIso="xxx", int n_data = 1, int rebin_const = 4) 
{
//   string szIso = "xxx";
  int number_of_lines = 3000;    
  double OX[number_of_lines], OY[number_of_lines], ERRY[number_of_lines], ERRX[number_of_lines];
//   const int n_data = 1;
  //const int n_param = 6;
  struct TData
  {
    double T12, errT12, Phi, errPhi, Pn, errPn ,lambda;
   // TString name;
  };
  
  TData data[n_data];
  
  double BG0, errBG0; double errTbeam;
  double BG1, errBG1;
  double rmin, rmax;
  
//   cout<<"----- szInFile -----"<< szInFile<<endl; 
  toks = szInFile.Tokenize(".");
  TString runnbr = ((TObjString* )toks->At(0))->GetString();
  cout<<"----- runnbr ----- "<< runnbr <<endl; 
//   toks = runnbr.Tokenize(".");
//   runnbr = ((TObjString* )toks->At(0))->GetString();
//   cout<<"----- runnbr ----- "<< runnbr<<endl;  
   
  string szTemp, szTemp1;
  
  ifstream fParam (szParInput, ios::in);
  
  fParam >> szTemp; cout<<szTemp<<endl;
  fParam >> szIso;  cout<<szIso<<endl;
  
  fParam >> szTemp; cout<<szTemp<<endl;
  fParam >> rmin;   cout<<rmin<<endl;
  fParam >> rmax;   cout<<rmax<<endl;
    
  fParam >> szTemp;cout<<szTemp<<endl;
  fParam >> Tbeam; cout<<Tbeam<<endl;
  fParam >> szTemp;  
  fParam >> errTbeam; cout<<"errTbeam "<<errTbeam<<endl;
  
  fParam >> szTemp;
  fParam >> Tdecay; cout<<"Tdecay "<<Tdecay<<endl;
  
  fParam >> szTemp;  
  fParam >> Tbg; cout<<"Tbg "<<Tbg<<endl;
  
  fParam >> szTemp;  
  fParam >> BG0; cout<<"BG "<<BG0<<endl;  
  fParam >> szTemp;  
  fParam >> errBG0; cout<<"errBG0 "<<errBG0<<endl;
  
  fParam >> szTemp;  
  fParam >> BG1; cout<<"BG1 "<<BG1<<endl;  
  fParam >> szTemp;  
  fParam >> errBG1; cout<<"errBG1 "<<errBG1<<endl;
  if (BG1 == 0) BG1=1e-15;
  
  
  TString szOutFile; 
  string szIsotope=szIso;  
  szOutFile = szIsotope + "_neutron"+".root";
  
  TbeamPrime = Tbeam + Tbg;
  TF1 *halflife2= new TF1("fit_lifetime2",fit_lifetime2, 0,TbeamPrime+Tdecay, 6);
// TF1 *halflife2= new TF1("fit_lifetime2",fit_lifetime2, 0,4, 4);    

  int par_number = 0;  
  for (int i=0;i<n_data;i++)
  {
    cout<<Form("-------------%i-------------",i)<<endl;
    par_number = i;  
    fParam >> szTemp1;
    fParam >> szTemp1;
    //cout<<szTemp1<<endl;
    fParam >> data[i].Phi; cout<<"Phi"<<i+1<<" "<<data[i].Phi<<endl;  
    fParam >> szTemp;  
    //cout<<szTemp<<endl;
    fParam >> data[i].errPhi; cout<<"errPhi"<<i+1<<" "<<data[i].errPhi<<endl;
    halflife2->SetParName(par_number,Form("PHI_%i",i+1));
    halflife2->SetParameter(par_number, data[i].Phi);  
    halflife2->SetParLimits (par_number,data[i].Phi*(1-data[i].errPhi), data[i].Phi*(1+data[i].errPhi));
  
    par_number++;
    fParam >> szTemp;
    fParam >> data[i].T12; cout<<"T12"<<i+1<<" "<<data[i].T12<<endl;
    data[i].lambda = log_2/data[i].T12; cout<<"lambda"<<i+1<<" "<<data[i].lambda<<endl;    
    fParam >> szTemp1;  
    //cout<<szTemp1<<endl;    
    fParam >> data[i].errT12; cout<<"errT12"<<i+1<<" "<<data[i].T12<<endl;
    halflife2->SetParName(par_number,Form("T1/2_%i",i+1));
    halflife2->SetParameter(par_number, data[i].lambda);
    halflife2->SetParLimits (par_number, data[i].lambda*(1-data[i].errT12), data[i].lambda*(1+data[i].errT12));

     par_number++;
    fParam >> szTemp1;  
    fParam >> data[i].Pn; cout<<"Pn"<<i+1<<" "<<data[i].Pn<<endl;  
    fParam >> szTemp1;      
    fParam >> data[i].errPn; cout<<"errPn"<<i+1<<" "<<data[i].errPn<<endl;  
    halflife2->SetParName(par_number,Form("Pn1_%i",i+1));
    halflife2->SetParameter(par_number, data[i].Pn);
    if (data[i].errPn == 0) {halflife2->FixParameter(par_number, data[i].Pn);}
    else{      
      halflife2->SetParLimits (par_number, data[i].Pn*(1-data[i].errPn), data[i].Pn*(1+data[i].errPn));
    };
  };
  
if (blRebinPhiBG0) {
 BG0*=rebin_const/4;
 data[0].Phi*= rebin_const/4;
 cout<<"Rebining "<<rebin_const/4<<endl;
 cout<<"BG "<<BG0<<" Phi "<<data[0].Phi<<endl;  
};
   
 par_number++;  
 halflife2->SetParName(par_number,"BG0");
 halflife2->SetParameter(par_number, BG0);
 halflife2->SetParLimits (par_number,BG0*(1-errBG0), BG0*(1+errBG0));
 par_number++;  
 halflife2->SetParName(par_number,"BG1");
 halflife2->SetParameter(par_number, BG1);
 halflife2->SetParLimits (par_number,TMath::Abs(BG1)*(1-errBG1), TMath::Abs(BG1)*(1+errBG1)); 
 cout<<"------------BG0 done---------"<<endl;
 par_number++;  
 halflife2->SetParName(par_number,"Tbeam");
 halflife2->SetParameter(par_number, Tbeam);
 halflife2->SetParLimits(par_number,Tbeam*(1-errTbeam), Tbeam*(1+errTbeam)); 
 cout<<"------------BG1 done---------"<<endl;
 fParam.close(); 
  
//   cout<<"Number of Points "<<number_of_lines<<endl;
  
 
  
 TFile* f1=new TFile(szInFile);
 cout<<"-----------file read--------"<<endl;
 TH1F* ne_time;// = f1->Get("Time/hTime__n_nim__sec");
 ne_time = (TH1F*) f1->Get("Time/hTime__n_nim__sec");
 ne_time->Rebin(rebin_const);
 number_of_lines=3000;//ne_time->GetMaximumBin();
 for (int i=0; i<=number_of_lines; i++)
  {
      OX[i] = ne_time->GetBinCenter(i);
      OY[i] = ne_time->GetBinContent(i);
      ERRY[i] = sqrt(OY[i]);
      ERRX[i] = 0;
    }; 
  
//   f1->Close();
  
  
  TCanvas      *c1 = new TCanvas("HalfLive",szInFile,200,10,700,500);
  c1->SetGrid();
  c1->cd();
  TGraphErrors *gr = new TGraphErrors(number_of_lines,OX,OY,ERRX,ERRY);   
  gr->SetName(Form("gr%s",szIso.c_str()));
  szIsotope = szIsotope + "_neutron activity";
  
  std::stringstream sTitleName ;
  sTitleName << szIsotope.c_str()<<"_"<<runnbr;
  gr->SetTitle(sTitleName.str().c_str());
  TFile        *f2 = new TFile(szOutFile,"RECREATE");
  

  int max_value=ne_time->GetMaximum()+ne_time->GetMaximum()*0.05;
  int min_value=ne_time->GetMinimum();
if (plot_option == 1) {
   gr->Fit(halflife2,"m","",rmin,rmax);
   gr->GetXaxis()->SetRangeUser(0,TbeamPrime+Tdecay);
   gr->GetXaxis()->SetTitle("Time, [s]");   
   gr->GetXaxis()->SetTitleSize(0.03);
   gr->GetXaxis()->CenterTitle(true);
   
   gr->GetYaxis()->SetTitle("Count/160#mus");   
   gr->GetYaxis()->SetTitleSize(0.03);  
   gr->GetYaxis()->CenterTitle(true);
   gr->GetYaxis()->SetTitleOffset(1.5);
   
   gStyle->SetOptFit();
   
   gr->GetYaxis()->SetRangeUser(BG0-20,max_value+70);
//    c1->cd();
   gr->Draw("ap*.");//PAL
   
   TPaveText *t12text=new TPaveText (.05,.1,.95,.8);
   t12text->AddText("test");
   t12text->Draw("same");
   
   
   gr->Write();
}else {
   ne_time->GetXaxis()->SetRangeUser(0,TbeamPrime+Tdecay);
   ne_time->Fit(halflife2,"r");
   ne_time->Draw();
   ne_time->Write();
};
  
  f2->Close();  

 
  std::stringstream OutputFile;
  OutputFile << szIso.c_str()<<"_run"<<runnbr<< ".fit" ;
  std::ofstream fdataout ;
  fdataout .open(OutputFile.str().c_str());
//   ofstream fdataout (Form("%s_fit.txt",szIso.c_str()));
  
  int pp = 0;
  double dblTmp;
  for (int i=0;i<n_data;i++)
  {
    cout<<"----- Fit -----"<<endl;
    fdataout<<"----- Fit -----"<<endl;
    
    cout<<"Phi"<<i+1<<" "<<halflife2->GetParameter(pp)<<" +/- "<<halflife2->GetParError(pp) <<endl;
    fdataout<<"Phi"<<i+1<<" "<<halflife2->GetParameter(pp)<<" +/- "<<halflife2->GetParError(pp) <<endl;
    pp++;
    dblTmp = log_2/halflife2->GetParameter(pp);
    cout<<"T12"<<i+1<<" "<<dblTmp<<" +/- "<< (dblTmp*(halflife2->GetParError(pp)/halflife2->GetParameter(pp)))<<endl;
    fdataout<<"T12"<<i+1<<" "<<dblTmp<<" +/- "<< (dblTmp*(halflife2->GetParError(pp)/halflife2->GetParameter(pp)))<<endl;
    pp++;
    cout<<"Pn"<<i+1<<" "<<halflife2->GetParameter(pp)<<" +/- "<<halflife2->GetParError(pp) <<endl;
    fdataout<<"Pn"<<i+1<<" "<<halflife2->GetParameter(pp)<<" +/- "<<halflife2->GetParError(pp) <<endl;
    pp++;
   };
   double chindf= halflife2->GetNDF()/(halflife2->GetNDF()*1.0);
   cout<<"NDF="<<halflife2->GetNDF()<<"; Chisquare="<< halflife2->GetChisquare()<<"; GetProb()="<<halflife2->GetProb() <<endl;
   cout<<"Chisquare/NDF="<<halflife2->GetNDF()/halflife2->GetNDF()<<endl;
   fdataout<<"NDF="<<halflife2->GetNDF()<<"; Chisquare="<< halflife2->GetChisquare()<<"; GetProb()="<<halflife2->GetProb() <<endl;
   fdataout<<"Chisquare/NDF="<<chindf<<endl; 
   cout<<"-----     -----"<<endl;
   cout<<"Rebining "<<rebin_const<<endl;
   cout<<"-----     -----"<<endl;
   fdataout<<"-----     -----"<<endl;
   fdataout<<"Rebining "<<rebin_const<<endl;
   fdataout<<"-----     -----"<<endl;
   
   fdataout<<szInFile<<" "<<szParInput<<std::endl;
}
   

