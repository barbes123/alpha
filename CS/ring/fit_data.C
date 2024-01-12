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
#include <string>

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

double MyFit(Double_t *x, Double_t *par) 
{  
    return  par[2]*(exp((-1)*(par[5] + ((x[0]) - par[5]))* par[1])* (-1 + exp(par[5]* par[1]))* par[0]);;
};



void fit_data() 
{
  int number_of_lines = 50;    
  double OX[number_of_lines], OY12[number_of_lines], OY13[number_of_lines], OY23[number_of_lines], ERRY[number_of_lines], ERRX[number_of_lines];
  
// TF1 *halflife2= new TF1("fit_lifetime2",fit_lifetime2, 0,TbeamPrime+Tdecay, 6);
 TF1 *fit1= new TF1("MyFit",MyFit, 0, 10, 4);    

 std::stringstream ringfile;
 ringfile<<"ring_sim.dat";
 std::ifstream infile(ringfile.str().c_str());
 
 int i = 0;
 
 while (infile.good()) {
      std::string oneline;
      std::getline(infile, oneline);
      if (!infile.good()) continue;
      if (oneline[0] == '#') continue; // ignore lines stating with #
      if (oneline.empty())   continue; // ignore empty lines

      std::istringstream is(oneline);
 
      int coinc_id = 0; int time_corr = 0;
      is >> OX[i] >> OY12[i] >> OY13[i] >> OY23[i];
      ERRX[i] = 0;
      ERRY[i] = 0;
      std::cout << OX[i]<<" " <<" "<< OY12[i] <<" "<< OY13[i]<<" " <<OY23[i]<<" \n";
      i++;
 
  }
  infile.close();
 

  TCanvas      *c1 = new TCanvas("Data","Data",200,10,700,500);
  c1->SetGrid();
  c1->cd();
//   TGraphErrors *gr = new TGraphErrors(number_of_lines,OX,OY12,ERRX,ERRY);
  TGraph *gr = new TGraphErrors(number_of_lines,OX,OY12);   
  gr->SetName("Data");
   gr->Draw("ap*.");//PAL

}
   

