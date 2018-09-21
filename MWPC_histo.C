#ifndef __CINT__
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>
// #include <TROOT.h>
// #include "TSystem.h"
// #include "TFile.h"
// #include "TProfile.h"
// #include "TChain.h"
#include "TH1F.h"
// #include "TH1D.h"
// #include "TH2F.h"
// #include "TF1.h"
// #include "TGaxis.h"
// #include "TRandom.h"
// #include "TNtuple.h"
#include "TCanvas.h"
// #include "TPolyLine.h"
// #include "TLine.h"
// #include "TArrow.h"
// #include "TStyle.h"  
// #include "TGraphErrors.h"
// #include "TGraph.h"
// #include "TBranch.h"
// #include "TLegend.h"
// #include "TLatex.h"
// #include "TEllipse.h"
// #include "TMarker.h"
// #include "ANAPATH.h"
// #include "CommonParameters.h"
// #include "ADC_Thresholds.h"
// #include "TDC_Windows.h"
// #include "Cuts_and_Windows.h"
// #include "MWPC_Thr.h"
#endif

#include "intersect.cxx"
#include <ctime>

using namespace std;


void MWPC_histo(Int_t Run_Number = 3994, Int_t first_event = 0, Int_t last_event = 100){


  // Set title for input file.
  char input_file_title[100];
  sprintf(input_file_title,"Event_Angle_Run%d__Event%d_to_Event%d.csv",Run_Number, first_event, last_event);
  
  ifstream input_file(input_file_title);




  const int n_bins = 10;
  const double bin_edges[n_bins+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

  // Histogram pointers
  TH1F *C2X_R = new TH1F("C2X Right", "C2X Right", n_bins, bin_edges);
  TH1F *C2Y_R = new TH1F("C2Y Right", "C2Y Right", n_bins, bin_edges);
  TH1F *C3X_R = new TH1F("C3X Right", "C3X Right", n_bins, bin_edges);
  TH1F *C3Y_R = new TH1F("C3Y Right", "C3Y Right", n_bins, bin_edges);
  TH1F *C4X_R = new TH1F("C4X Right", "C4X Right", n_bins, bin_edges);
  TH1F *C4Y_R = new TH1F("C4Y Right", "C4Y Right", n_bins, bin_edges);

  TH1F *C2X_L = new TH1F("C2X Left", "C2X Left", n_bins, bin_edges);
  TH1F *C2Y_L = new TH1F("C2Y Left", "C2Y Left", n_bins, bin_edges);
  TH1F *C3X_L = new TH1F("C3X Left", "C3X Left", n_bins, bin_edges);
  TH1F *C3Y_L = new TH1F("C3Y Left", "C3Y Left", n_bins, bin_edges);
  TH1F *C4X_L = new TH1F("C4X Left", "C4X Left", n_bins, bin_edges);
  TH1F *C4Y_L = new TH1F("C4Y Left", "C4Y Left", n_bins, bin_edges);  




  ////////////////////////
  //        TEST        //
  ////////////////////////

  TLine *C2X_R_arr[3];

  for(int i = 0; i<3; i++){
    C2X_R_arr[i] = new TLine(0,0,i,i);
  }





  // Read file, and fill histograms;
  string current_column;

  int column_counter = 0;
  string MWPC_side;

  while(input_file.good()){
    getline(input_file, current_column, ',');

    if(column_counter == 12){
      MWPC_side = current_column;
    }

    if(column_counter == 13){
      if(MWPC_side == "R"){
        C2X_R->Fill(double(stoi(current_column)) + 0.5);
      }
      else{
        C2X_L->Fill(double(stoi(current_column)) + 0.5);
      }
    }

    if(column_counter == 14){
      if(MWPC_side == "R"){
        C2Y_R->Fill(double(stoi(current_column)) + 0.5);
      }
      else{
        C2Y_L->Fill(double(stoi(current_column)) + 0.5);
      }
    }

    if(column_counter == 15){
      if(MWPC_side == "R"){
        C3X_R->Fill(double(stoi(current_column)) + 0.5);
      }
      else{
        C3X_L->Fill(double(stoi(current_column)) + 0.5);
      }
    }  

    if(column_counter == 16){
      if(MWPC_side == "R"){
        C3Y_R->Fill(double(stoi(current_column)) + 0.5);
      }
      else{
        C3Y_L->Fill(double(stoi(current_column)) + 0.5);
      }
    } 

    if(column_counter == 17){
      if(MWPC_side == "R"){
        C4X_R->Fill(double(stoi(current_column)) + 0.5);
      }
      else{
        C4X_L->Fill(double(stoi(current_column)) + 0.5);
      }
    }    

    if(column_counter == 18){
      if(MWPC_side == "R"){
        C4Y_R->Fill(double(stoi(current_column)) + 0.5);
      }
      else{
        C4Y_L->Fill(double(stoi(current_column)) + 0.5);
      }

      column_counter = 0; // Reset counter
    }            

    column_counter++;
  }


  TCanvas *c1;
  c1 = new TCanvas("Histo Display","MWPC Left/Right Histograms",50,50,1050,700);
  c1->Divide(4,3);
  c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

  c1->cd(1);       
  C2X_L->Draw();

  c1->cd(2);
  C2Y_L->Draw();

  c1->cd(3);
  C2X_R->Draw();

  c1->cd(4);
  C2Y_R->Draw();

  c1->cd(5);
  C3X_L->Draw();

  c1->cd(6);
  C3Y_L->Draw();

  c1->cd(7);
  C3X_R->Draw();

  c1->cd(8);
  C3Y_R->Draw();

  c1->cd(9);
  C4X_L->Draw();

  c1->cd(10);
  C4Y_L->Draw();

  c1->cd(11);
  C4X_R->Draw();

  c1->cd(12);
  C4Y_R->Draw();
 


  input_file.close();

  return;
}
