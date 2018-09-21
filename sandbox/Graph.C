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

#include "intersect.cxx"
#include "C2_Strip_transform.h"
#include "Channel_to_Strip.h"
#include "SFT_functions.h"

using namespace std;
 
void Graph() {
   //Draw a simple graph
   // To see the output of this macro, click begin_html <a href="gif/graph.gif">here</a>. end_html
   //Author: Rene Brun

   float R_TARGET = 29.0; 
   TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0);
   TH2F *h_NULL = new TH2F("NULL", "NULL", 100,-50,50,100,-50,50); 

   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,500,500);

   //c1->SetFillColor(42);
   //c1->SetGrid();

   const Int_t n = 18;
   Double_t x[n], y[n];
   vector<double> vec_xx;
   vector<double> vec_yy;
   //for (Int_t i=0;i<n;i++) {
     //x[i] = 5;
     //y[i] = i;

     //x[i] = i*0.1;
     //y[i] = 10*sin(x[i]+0.2);
     //printf(" i %i %f %f \n",i,x[i],y[i]);
   //}
   
   vec_xx.push_back(7.75);          vec_yy.push_back(-26.35);
   vec_xx.push_back(7.75);          vec_yy.push_back(-26.35);
   vec_xx.push_back(7.75);          vec_yy.push_back(-23.25);
   vec_xx.push_back(7.75);          vec_yy.push_back(-23.25);
   vec_xx.push_back(7.75);          vec_yy.push_back(-20.15);
   vec_xx.push_back(7.75);          vec_yy.push_back(-10.85);
   vec_xx.push_back(10.85);         vec_yy.push_back(-17.05);
   vec_xx.push_back(10.85);         vec_yy.push_back(-10.85);
   //vec_xx.push_back(8.179);         vec_yy.push_back(-46.384);
   vec_xx.push_back(8.179);         vec_yy.push_back(-46.384);
   vec_xx.push_back(8.179);         vec_yy.push_back(-46.384);
   //vec_xx.push_back(20.15);         vec_yy.push_back(-1.55);
   //vec_xx.push_back(20.15);         vec_yy.push_back(1.55);
   //vec_xx.push_back(23.25);         vec_yy.push_back(-4.65);
   //vec_xx.push_back(23.25);         vec_yy.push_back(4.65);
   //vec_xx.push_back(23.25);         vec_yy.push_back(7.75);
   //vec_xx.push_back(23.25);         vec_yy.push_back(10.85);
   //vec_xx.push_back(26.35);         vec_yy.push_back(4.65);

   /*
   x[0] = 4.65;           y[0] = -26.35;
   x[1] = 4.65;           y[1] = -26.35;
   x[2] = 7.75;           y[2] = -23.25;
   x[3] = 7.75;           y[3] = -23.25;
   x[4] = 7.75;           y[4] = -20.15;
   x[5] = 10.85;          y[5] = -20.15;
   x[6] = 10.85;          y[6] = -17.05;
   x[7] = 13.95;          y[7] = -10.85;
   x[8] = 17.05;          y[8] = -10.85;
   x[9] = 17.05;          y[9] = -7.75;
   x[10] = 17.05;         y[10] = -4.65;
   x[11] = 20.15;         y[11] = -1.55;
   x[12] = 20.15;         y[12] = 1.55;
   x[13] = 23.25;         y[13] = -4.65;
   x[14] = 23.25;         y[14] = 4.65;
   x[15] = 23.25;         y[15] = 7.75;
   x[16] = 23.25;         y[16] = 10.85;
   x[17] = 26.35;         y[17] = 4.65;
  */

   //TGraph *gr = new TGraph(n,x,y);
   TGraph *gr = new TGraph(vec_xx.size(),&vec_xx[0],&vec_yy[0]);
   gr->SetLineColor(2);
   gr->SetLineWidth(4);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   //gr->SetTitle("a simple graph");
   //gr->GetXaxis()->SetTitle("X title");
   gr->GetYaxis()->SetTitle("Y title");
   gr->GetXaxis()->SetLimits(-50.,50.);
   //gr->GetYaxis()->SetRangeUser(-50.,50.);
   gr->Draw("AP");
   gr->Fit("pol1","Q");
   cout << endl;

   c1->cd();
   h_NULL->Draw();
   ell_Target->Draw("same");
   gr->Draw("sameP");

  

   // TCanvas::Update() draws the frame, after which one can change it
   //c1->Update();
   //c1->GetFrame()->SetFillColor(21);
   //c1->GetFrame()->SetBorderSize(12);
   //c1->Modified();
}

