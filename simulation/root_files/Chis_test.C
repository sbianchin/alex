#ifndef __CINT__
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <TROOT.h>
#include "TSystem.h"
#include "TFile.h"
#include "TProfile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TPolyLine.h"
#include "TLine.h"
#include "TArrow.h"
#include "TStyle.h"  
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TBranch.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TEllipse.h"
#include "TMarker.h"
#include "ANAPATH.h"
#include "CommonParameters.h"
#include "ADC_Thresholds.h"
#include "TDC_Windows.h"
#include "Cuts_and_Windows.h"
//#include "MWPC_Thr.h"
#include "MWPC_Thr2.h"
#include "Pedestals.h"
#endif


using namespace std;
 
void Chis_test(){

   gStyle->SetOptStat(0);

   float R_TARGET = 29.0; 
   TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0);
   TH2F *h_NULL = new TH2F("NULL", "NULL", 100,-50,50,100,-50,50); 

   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,700);


   int bar_number;
   double ChiS;
   int ndf;
   double a_fit, b_fit;
   double sum = 0;
   vector<double> vec_xx;        vec_xx.clear();
   vector<double> vec_yy;        vec_yy.clear();
   vector<double> vec_ex;        vec_ex.clear();
   vector<double> vec_ey;        vec_ey.clear();
   vector<int> barn;             barn.clear();
   vector<double> Y_expected;    Y_expected.clear();

   cout << endl;

   while(bar_number<256){

      cout << "Enter Bar Number :   ";  
      cin >> bar_number;   
      //test.push_back(i);

      if(bar_number>255){
         cout << "Thank you!" << endl;
         break;
      }
      
      vec_xx.push_back(Xloc[bar_number]);
      vec_yy.push_back(Yloc[bar_number]);
      vec_ex.push_back(1.6);
      vec_ey.push_back(1.6);
      barn.push_back(bar_number);
   }

   TGraph *gr = new TGraphErrors(vec_xx.size(),&vec_xx[0],&vec_yy[0],&vec_ex[0],&vec_ey[0]);
   //TGraph *gr = new TGraph(vec_xx.size(),&vec_xx[0],&vec_yy[0]);
   TF1 *gr_fit = new TF1("gr", "pol1");
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->GetYaxis()->SetTitle("Y title");
   gr->GetXaxis()->SetLimits(-50.,50.);
   gr->Draw("AP");
   gr->Fit("pol1","Q");
  
   gr_fit = gr->GetFunction("pol1");
   gr_fit->SetLineWidth(2);
   gr_fit->SetLineColor(2);

   a_fit = gr_fit->GetParameter(1);
   b_fit = gr_fit->GetParameter(0);
   ChiS = gr_fit->GetChisquare();
   ndf = gr_fit->GetNDF();


   cout << endl;
   cout << endl;

   cout << "BAR NUMBERS : " << endl;
   cout << fixed;
   for(unsigned int j=0; j<barn.size(); j++){
      Y_expected.push_back((a_fit*Xloc[barn[j]])+b_fit);
      sum += pow(Yloc[barn[j]]-Y_expected[j],2);
      cout << setw(4) << barn[j] << "  ";
      cout << setw(6) << setprecision(2) << Xloc[barn[j]] << "  ";
      cout << setw(6) << setprecision(2) << Yloc[barn[j]] << "  ";
      cout << setw(6) << setprecision(2) << Y_expected[j] << "  ";
      cout << setw(6) << setprecision(2) << Yloc[barn[j]]-Y_expected[j] << "  ";
      cout << endl;
   }
   cout << "SUM = " << setw(6) << setprecision(2) << sum << "  ";
   cout << endl;
   cout << endl;


   cout << endl;
   cout << endl;

   cout << "ChiS = " << ChiS << endl;
   cout << "NdF = " << ndf << endl;
   cout << endl;
   cout << "a_fit = " << a_fit << endl;  
   cout << "b_fit = " << b_fit << endl;  
   cout << endl;
   cout << endl;

   char ChiS_string[30]; 
   char ndf_string[30];          
   sprintf(ChiS_string,"#chi^{2} = %3.2f", ChiS);
   sprintf(ndf_string,"NdF = %3i", ndf);

   TLatex *tex_ChiS;
   tex_ChiS = new TLatex(-45.,43.,ChiS_string);
   tex_ChiS->SetTextSize(0.05);
   tex_ChiS->SetLineWidth(2);

   TLatex *tex_ndf;
   tex_ndf = new TLatex(-45.,35.,ndf_string);
   tex_ndf->SetTextSize(0.05);
   tex_ndf->SetLineWidth(2);


   c1->cd();
   h_NULL->Draw();
   ell_Target->Draw("same");
   gr->Draw("sameP");
   tex_ChiS->Draw("same");
   tex_ndf->Draw("same");
}

