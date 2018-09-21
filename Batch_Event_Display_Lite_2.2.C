#ifndef __CINT__
#include <stdio.h>
#include <math.h>
#include <iomanip>
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
#endif

#include "intersect.cxx"
#include "C2_Strip_transform.h"
#include "Channel_to_Strip.h"
//#include "SFT_functions.h"
#include "Batch_Variables.h"
#include "Batch_Event_Display_functions.h"

//using namespace std;
 

//vector<double> GetZ(int Run_Number, int evt, double phi, int gap_to_fit, bool to_print, double C2X_centroid, double length_in_target); 
//vector<double> GetDeltaZ(int Run_Number, int evt, double phi, int gap_to_fit, bool to_print, double C2X_centroid, double length_in_target); 
//vector<double> GetTrackLength(int Run_Number, int evt, double phi, int gap_to_fit, bool to_print, double C2X_centroid, double length_in_target); 

void Batch_Event_Display_Lite_2_2(Int_t Run_Number=5, Int_t ievt_min=2, Int_t ievt_max=2, Int_t Switch_Display=1, Int_t Switch_Zlist=0){ 
  
  gROOT->SetBatch(1);

  gStyle->Clear();
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(11);
 
  Int_t TDC_min_TARGET = TARGET_TDC_min[0];
  Int_t TDC_max_TARGET = TARGET_TDC_max[0];
  Int_t ADC_TARGET_Thr = HG_TARGET_ADC_Thr[0];
  Int_t TDC_min_SFT = SFT_TDC_min[0];
  Int_t TDC_max_SFT = SFT_TDC_max[0];
  
  sprintf(ADC_cut,"(ADC >= %d)",HG_SFT_ADC_Thr[0]);

  sprintf(file_mapping,"../Mapping");

  sprintf(par_finput,"%s/%s",file_mapping,source_mapping);

  sprintf(par_finput2,"%s/MWPC_map2.txt",file_mapping);
  ifstream fdat(par_finput,ios::in);
  for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
  fdat.close();

  ifstream fdat2(par_finput2,ios::in);
  for(Int_t ii=0; ii<512; ii++) fdat2 >> par_temp2[ii];
  fdat2.close();
  
  sprintf(h_target_ADC_title,"Run %d (Event %d) -- ADC #geq %d)",Run_Number, ievt_min, ADC_TARGET_Thr);
  sprintf(h_target_TDC_title,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt_min, ADC_TARGET_Thr,TDC_min_TARGET,TDC_max_TARGET);

  sprintf(h_target_ADC_title2,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt_min, ADC_TARGET_Thr,TDC_min_TARGET, TDC_max_TARGET);
  sprintf(h_target_TDC_title2,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt_min, ADC_TARGET_Thr, TDC_min_TARGET,TDC_max_TARGET);
  sprintf(h_target_ADC_title3,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt_min, ADC_TARGET_Thr, TDC_min_TARGET,TDC_max_TARGET);
  sprintf(h_target_ADC_title4,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt_min, ADC_TARGET_Thr, TDC_min_TARGET,TDC_max_TARGET);
  
  sprintf(h_ADC_title,"(ADC offset = %d) | (%d < TDC < %d)  --  Run %d (Event %d)",SFT_ADC_Thr_HG_Offset, TDC_min_SFT,TDC_max_SFT,Run_Number,ievt_min);


  sprintf(h_target_TDC_copy_Name,"Event %d (Run %d)", ievt_min, Run_Number);
  sprintf(h_target_TDC_copy_Title,"ADC_HG #geq %d  |  %d #leq TDC #leq %d", ADC_TARGET_Thr, TDC_min_TARGET, TDC_max_TARGET);
 
  sprintf(h_TDC_selected_Name,"Event %d (Run %d)", ievt_min, Run_Number);
  sprintf(h_TDC_selected_Title,"ADC_HG #geq %d  |  %d #leq TDC #leq %d", ADC_TARGET_Thr+100, TDC_min_TARGET, TDC_max_TARGET);

  sprintf(h_GoodLG_Name,"Event %d (Run %d)", ievt_min, Run_Number);
  sprintf(h_GoodLG_Title,"ADC_LG #geq %d  |  %d #leq TDC #leq %d", ADC_TARGET_Thr,TDC_min_TARGET,TDC_max_TARGET);

  sprintf(h_TDC_selected2_Name,"Event %d (Run %d)", ievt_min, Run_Number);
  sprintf(h_TDC_selected2_Title,"ADC_HG #geq %d  |  %d #leq TDC #leq %d  (WEIGHTED)", ADC_TARGET_Thr+100,TDC_min_TARGET,TDC_max_TARGET);

  sprintf(h_GoodLG_weighted_Name,"Event %d (Run %d)", ievt_min, Run_Number);
  sprintf(h_GoodLG_weighted_Title,"ADC_LG #geq %d  |  %d #leq TDC #leq %d  (WEIGHTED)", ADC_TARGET_Thr,TDC_min_TARGET,TDC_max_TARGET);

  sprintf(h_Gap_Fibers_Name,"Event %d (Run %d)", ievt_min, Run_Number);
  sprintf(h_Gap_Fibers_Title,"ADC LG offset =  %d  |  %d #leq TDC #leq %d", TARGET_ADC_Thr_LG_Offset,TDC_min_TARGET,TDC_max_TARGET);


  TH2F *h_Circle = new TH2F("Test", "TARGET and TOF1", 500, -50, 50, 500, -50, 50);

  TH2F *h_target_ADC = new TH2F("Histo Fit 0",h_target_ADC_title,3000,-50,50,3000,-50,50);
  TH2F *h_target_TDC = new TH2F("Histo Fit 1",h_target_TDC_title,3000,-50,50,3000,-50,50);
  TH2F *h_target_TDC_copy = new TH2F(h_target_TDC_copy_Name,h_target_TDC_copy_Title,3000,-50,50,3000,-50,50);
  TH2F *h_TDC_selected = new TH2F(h_TDC_selected_Name, h_TDC_selected_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_TDC_selected2 = new TH2F(h_TDC_selected2_Name, h_TDC_selected2_Title, 500, -50, 50, 500, -50, 50);

  TH2F *h_target_ADC2 = new TH2F("Histo Fit 2",h_target_ADC_title2,3000,-50,50,3000,-50,50);
  TH2F *h_target_TDC2 = new TH2F("Histo Fit 3",h_target_TDC_title2,3000,-50,50,3000,-50,50);

  TH2F *h_target_ADC3 = new TH2F("Histo Fit 4",h_target_ADC_title3,3000,-50,50,3000,-50,50);
  TH2F *h_target_ADC4 = new TH2F("Histo Fit 5",h_target_ADC_title4,3000,-50,50,3000,-50,50);
  TH2F *h_target_ADCA = new TH2F("Histo Fit 6",h_target_ADC_title3,3000,-50,50,3000,-50,50);

  TH2F *h_Target = new TH2F("Test Target", "Target", 3000, -50, 50, 3000, -50, 50);
  TH2F *h_Target_Center = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);

  TH2F *h_kaon = new TH2F("Kaon", "Kaon", 500, -50, 50, 500, -50, 50);
  TH2F *h_kaon_copy = new TH2F("Kaon Copy", "Kaon Copy", 500, -50, 50, 500, -50, 50);
    
  TH2F *h_max = new TH2F("Max", "Max", 500, -50, 50, 500, -50, 50);
  TH2F *h_max_copy = new TH2F("Max Copy", "Max Copy", 500, -50, 50, 500, -50, 50);

  TH2F *h_GoodLG = new TH2F(h_GoodLG_Name, h_GoodLG_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_GoodLG_weighted = new TH2F(h_GoodLG_weighted_Name, h_GoodLG_weighted_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_GoodLG_copy = new TH2F("Good LG Copy", "Good LG Copy", 500, -50, 50, 500, -50, 50);

  TH2F *h_TOF1 = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_closest = new TH2F("h_TOF1_closest", "h_TOF1_closest", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_rotate = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50); // ROTATE_CHANGE
  TH2F *h_TOF1_closest_rotate = new TH2F("h_TOF1_closest", "h_TOF1_closest", 500, -50, 50, 500, -50, 50); //ROTATE_CHANGE
  
  TH2F *h_int_TDC = new TH2F("h_int_TDC", "h_int_TDC", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_selected = new TH2F("h_int_TDC_selected", "h_int_TDC_selected", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_selected_weighted = new TH2F("h_int_TDC_selected_weighted", "h_int_TDC_selected_weighted", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_Gap_Fibers = new TH2F("h_int_TDC_Gap_Fibers", "h_int_TDC_Gap_Fibers", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_Gap_Fibers_SFT = new TH2F("h_int_TDC_Gap_Fibers_SFT","h_int_TDC_Gap_Fibers_SFT", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_Gap_Fibers_rotate = new TH2F("h_int_TDC_Gap_Fibers", "h_int_TDC_Gap_Fibers", 500, -50, 50, 500, -50, 50);  // ROTATE_CHANGE
  TH2F *h_int_TDC_Gap_Fibers_SFT_rotate = new TH2F("h_int_TDC_Gap_Fibers_SFT","h_int_TDC_Gap_Fibers_SFT", 500, -50, 50, 500, -50, 50);  // ROTATE_CHANGE
  TH2F *h_int_TDC_SFT = new TH2F("h_int_TDC_SFT", "h_int_TDC_SFT", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_TARGET = new TH2F("h_int_TDC_TARGET", "h_int_TDC_TARGET", 500, -50, 50, 500, -50, 50);

  TH2F *h_int_GoodLG = new TH2F("h_int_GoodLG", "h_int_GoodLG", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_GoodLG_weighted = new TH2F("h_int_GoodLG_weighted", "h_int_GoodLG_weighted", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_GoodLG_SFT = new TH2F("h_int_GoodLG_SFT", "h_int_GoodLG_SFT", 500, -50, 50, 500, -50, 50);
    
  TH2F *h_TDC_Gap_Fibers = new TH2F(h_Gap_Fibers_Name, h_Gap_Fibers_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_K_Stop_Bars = new TH2F(h_Gap_Fibers_Name, h_Gap_Fibers_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_TDC_Gap_Fibers_copy = new TH2F(h_Gap_Fibers_Name, h_Gap_Fibers_Title, 500, -50, 50, 500, -50, 50); //ROTATE_CHANGE
  TH2F *h_TDC_Gap_Fibers_kaon = new TH2F(h_Gap_Fibers_Name, h_Gap_Fibers_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_TARGET_LG_Blue = new TH2F(h_Gap_Fibers_Name, h_Gap_Fibers_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_Centroid = new TH2F("Centroid", "Centroid", 500, -50, 50, 500, -50, 50);
  
  TEllipse *ell = new TEllipse(0, 0, R_TOF1, 0);
  TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0);
  TEllipse *ell_L1 = new TEllipse(0, 0, R_SFT_L1, 0);
  
  TLine *HorizontalAxis = new TLine(0., 0., 50., 0.);

  //Int_t TDC_LE_TARGET[256];

  sprintf(path_input,"%s",path_merged);          
  sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);
  sprintf(footer,"Event_Display_MS.C  --  Run %d ; Event %d",Run_Number,ievt_min);
  sprintf(run_string,"Run %d ; Event %d",Run_Number,ievt_min);
  sprintf(event_string,"Run %d ; Event %d",Run_Number,ievt_min);
  

  TChain *fChain= new TChain("Tree");   
  fChain->Add(Name_finput);   
  fChain->SetMakeClass(1);              

  fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);    fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
  fChain->SetBranchAddress("ADC_Low_TARGET",adc_low_target);      fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
  fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);        fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
  fChain->SetBranchAddress("TDC_TE_TARGET",tdc_te_target);        fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

  fChain->SetBranchAddress("ADC_TOF1U",ADC_tof1U);
  fChain->SetBranchAddress("ADC_TOF1D",ADC_tof1D);
  fChain->SetBranchAddress("TDC_TOF1U",TDC_tof1U);
  fChain->SetBranchAddress("TDC_TOF1D",TDC_tof1D); 
  
  fChain->SetBranchAddress("ADC_TOF2AO",ADC_tof2AO);
  fChain->SetBranchAddress("ADC_TOF2AI",ADC_tof2AI);
  fChain->SetBranchAddress("ADC_TOF2BO",ADC_tof2BO);
  fChain->SetBranchAddress("ADC_TOF2BI",ADC_tof2BI);    
  fChain->SetBranchAddress("TDC_TOF2AO",TDC_tof2AO);
  fChain->SetBranchAddress("TDC_TOF2AI",TDC_tof2AI);
  fChain->SetBranchAddress("TDC_TOF2BO",TDC_tof2BO);
  fChain->SetBranchAddress("TDC_TOF2BI",TDC_tof2BI);    

  fChain->SetBranchAddress("MWPCADC",MwpcADC);

  fChain->SetBranchAddress("ADC_C2X_R",adc_c2x_r);
  fChain->SetBranchAddress("ADC_C2X_L",adc_c2x_l);
  fChain->SetBranchAddress("ADC_C2Y_R",adc_c2y_r);
  fChain->SetBranchAddress("ADC_C2Y_L",adc_c2y_l);
  fChain->SetBranchAddress("ADC_C3X_R",adc_c3x_r);
  fChain->SetBranchAddress("ADC_C3X_L",adc_c3x_l);
  fChain->SetBranchAddress("ADC_C3Y_R",adc_c3y_r);
  fChain->SetBranchAddress("ADC_C3Y_L",adc_c3y_l);
  fChain->SetBranchAddress("ADC_C4X_R",adc_c4x_r);
  fChain->SetBranchAddress("ADC_C4X_L",adc_c4x_l);
  fChain->SetBranchAddress("ADC_C4Y_R",adc_c4y_r);
  fChain->SetBranchAddress("ADC_C4Y_L",adc_c4y_l);  

  fChain->SetBranchAddress("TDC_Ck", tdc_Ck);
  fChain->SetBranchAddress("TDC_Cpi", tdc_Cpi);

  fChain->SetBranchAddress("EvFlag", Event_flag);
 
  
  int count = 0; 

  ZZ.clear();
  DZ.clear();
  Track_Length.clear();
  Z_selected.clear();
  
  Z_selected_sorted.clear();

  sprintf(output,"RUN_%d_DATA.dat",Run_Number);
  ofstream fout;
  fout.open(output);

  sprintf(output2,"Run_%d_Zlist.dat",Run_Number);
  ofstream fout2;
  fout2.open(output2);

  cout << endl;



  for(int ivt=ievt_min; ivt<ievt_max+1; ivt++){
    fChain->GetEntry(ivt);  

    if((ivt-ievt_min)%10000==0 && (ivt-ievt_min)!=0) cout << "***  " << (ivt-ievt_min) << " events done!" << endl;

    for(int i=0; i<256; i++)  HG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Thr_HG[i]) + TARGET_ADC_Thr_HG_Offset;
    for(int i=0; i<256; i++)  LG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Thr_LG[i]) + TARGET_ADC_Thr_LG_Offset;
    for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;
    
    for (int j_TARGET=0; j_TARGET<256; j_TARGET++){
      ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-HG_TARGET_ADC_Thr[j_TARGET];
      ADC_High_TARGET_ped[j_TARGET]=adc_high_target[j_TARGET]-TARGET_ADC_Thr_HG[j_TARGET];
      ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-LG_TARGET_ADC_Thr[j_TARGET];
      ADC_Low_TARGET_ped[j_TARGET]=adc_low_target[j_TARGET]-TARGET_ADC_Thr_LG[j_TARGET]; 
      //TDC_LE_TARGET[j_TARGET]=tdc_le_target[j_TARGET][0];
    } 
  
    for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
      ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-HG_SFT_ADC_Thr[j_SFT];
    }

    for(int i=0; i<12; i++){
      ADC_TOF1[i] = ADC_tof1U[i]-TOF1U_ADC_Thr[i];
      ADC_TOF1[i+12] = ADC_tof1D[i]-TOF1D_ADC_Thr[i];
    }
    
    for (Int_t j_TOF2=0; j_TOF2<12; j_TOF2++) {
      ADC_TOF2AO[j_TOF2] = ADC_tof2AO[j_TOF2]-TOF2AO_ADC_Thr[j_TOF2];
      ADC_TOF2BO[j_TOF2] = ADC_tof2BO[j_TOF2]-TOF2AI_ADC_Thr[j_TOF2];
      ADC_TOF2AI[j_TOF2] = ADC_tof2AI[j_TOF2]-TOF2BO_ADC_Thr[j_TOF2];
      ADC_TOF2BI[j_TOF2] = ADC_tof2BI[j_TOF2]-TOF2BI_ADC_Thr[j_TOF2];
    }

    //for (Int_t j_MWPC=0; j_MWPC<512; j_MWPC++) {
      //MWPCADC[j_MWPC] = MwpcADC[j_MWPC]-MWPC_ADC_Thr[j_MWPC];
    //}

    for(Int_t j_C2=0; j_C2<56; j_C2++){
      ADC_C2X_R[j_C2] = double(adc_c2x_r[j_C2])-(ADC_C2X_Thr_R[j_C2]+MWPC_Thr_Offset_C2X);
      ADC_C2X_L[j_C2] = double(adc_c2x_l[j_C2])-(ADC_C2X_Thr_L[j_C2]+MWPC_Thr_Offset_C2X);
      //ADC_C2X_R[j_C2] = adc_c2x_r[j_C2]-ADC_C2X_Thr[j_C2];
      //ADC_C2X_L[j_C2] = adc_c2x_l[j_C2]-ADC_C2X_Thr[j_C2];
    }

    for(Int_t j_C3=0; j_C3<64; j_C3++){
      ADC_C3X_R[j_C3] = double(adc_c3x_r[j_C3])-(ADC_C3X_Thr_R[j_C3]+MWPC_Thr_Offset_C3X);
      ADC_C3X_L[j_C3] = double(adc_c3x_l[j_C3])-(ADC_C3X_Thr_L[j_C3]+MWPC_Thr_Offset_C3X);
      //ADC_C3X_R[j_C3] = adc_c3x_r[j_C3]-ADC_C3X_Thr[j_C3];
      //ADC_C3X_L[j_C3] = adc_c3x_l[j_C3]-ADC_C3X_Thr[j_C3];
    }

 
    for(Int_t j_C4=0; j_C4<72; j_C4++){
      ADC_C4X_R[j_C4] = double(adc_c4x_r[j_C4])-(ADC_C4X_Thr_R[j_C4]+MWPC_Thr_Offset_C4X);
      ADC_C4X_L[j_C4] = double(adc_c4x_l[j_C4])-(ADC_C4X_Thr_L[j_C4]+MWPC_Thr_Offset_C4X);
      //ADC_C4X_R[j_C4] = adc_c4x_r[j_C4]-ADC_C4X_Thr[j_C4];
      //ADC_C4X_L[j_C4] = adc_c4x_l[j_C4]-ADC_C4X_Thr[j_C4];
    }

    for(Int_t j_CY=0; j_CY<16; j_CY++){
      //ADC_C2Y_R[j_CY] = adc_c2y_r[j_CY]-ADC_C2Y_Thr[j_CY];
      //ADC_C2Y_L[j_CY] = adc_c2y_l[j_CY]-ADC_C2Y_Thr[j_CY];
      //ADC_C3Y_R[j_CY] = adc_c3y_r[j_CY]-ADC_C3Y_Thr[j_CY];
      //ADC_C3Y_L[j_CY] = adc_c3y_l[j_CY]-ADC_C3Y_Thr[j_CY];
      //ADC_C4Y_R[j_CY] = adc_c4y_r[j_CY]-ADC_C4Y_Thr[j_CY];
      //ADC_C4Y_L[j_CY] = adc_c4y_l[j_CY]-ADC_C4Y_Thr[j_CY];
      ADC_C2Y_R[j_CY] = double(adc_c2y_r[j_CY])-(ADC_C2Y_Thr_R[j_CY]+MWPC_Thr_Offset_C2Y);
      ADC_C2Y_L[j_CY] = double(adc_c2y_l[j_CY])-(ADC_C2Y_Thr_L[j_CY]+MWPC_Thr_Offset_C2Y);
      ADC_C3Y_R[j_CY] = double(adc_c3y_r[j_CY])-(ADC_C3Y_Thr_R[j_CY]+MWPC_Thr_Offset_C3Y);
      ADC_C3Y_L[j_CY] = double(adc_c3y_l[j_CY])-(ADC_C3Y_Thr_L[j_CY]+MWPC_Thr_Offset_C3Y);
      ADC_C4Y_R[j_CY] = double(adc_c4y_r[j_CY])-(ADC_C4Y_Thr_R[j_CY]+MWPC_Thr_Offset_C4Y);
      ADC_C4Y_L[j_CY] = double(adc_c4y_l[j_CY])-(ADC_C4Y_Thr_L[j_CY]+MWPC_Thr_Offset_C4Y);
    }    

    for (Int_t j_TDCTOF=0; j_TDCTOF<12; j_TDCTOF++) {
      TDC_TOF1U[j_TDCTOF] = TDC_tof1U[j_TDCTOF];
      TDC_TOF1D[j_TDCTOF] = TDC_tof1D[j_TDCTOF];
      TDC_TOF2AO[j_TDCTOF] = TDC_tof2AO[j_TDCTOF];
      TDC_TOF2BO[j_TDCTOF] = TDC_tof2BO[j_TDCTOF];
      TDC_TOF2AI[j_TDCTOF] = TDC_tof2AI[j_TDCTOF];
      TDC_TOF2BI[j_TDCTOF] = TDC_tof2BI[j_TDCTOF];
    }  
  
  
    

    //********* GOOD TARGET EVENTS
    Good_TARGET_Event = false;
    count_TARGET_evts = 0; 
    for(int i=0; i<256; i++)
    {
       if((ADC_High_TARGET[i]>=0 && tdc_le_target[i][0]>=TARGET_TDC_min[i] && tdc_le_target[i][0]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]>=0 && tdc_le_target[i][1]>=TARGET_TDC_min[i] && tdc_le_target[i][1]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]>=0 && tdc_le_target[i][2]>=TARGET_TDC_min[i] && tdc_le_target[i][2]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]>=0 && tdc_le_target[i][3]>=TARGET_TDC_min[i] && tdc_le_target[i][3]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]<0  && ADC_Low_TARGET[i]>0))
      { 
       count_TARGET_evts++;
      }
    }

    if(count_TARGET_evts >= n_hit) Good_TARGET_Event = true;

    //if(enable_cout == 9) Good_TARGET_Event = true;
    
    //********* GOOD TOF EVENTS
    for(int i=0; i<12; i++){
      Good_TOF1_ADC[i] = false;   
      Good_TOF2_ADC[i] = false;
      Good_TOF1_TDC[i] = false;   
      Good_TOF2_TDC[i] = false;
      Good_TOF1[i] = false;       
      Good_TOF2[i] = false;
      Good_TOFs[i] = false;
    }
    
    Good_TOF_Event = false;

  
 

    ////Event_On_Blacklist = false;
    
    //if(blacklist.fail()){
      //cout << "Error: Could not read blacklist file." << endl;
    //}
    //else{
    //    while(getline(blacklist,Current_Event)){
    //  sscanf(Current_Event.c_str(), "%d", &current_event);  
    //    if(current_event == ivt){
    //      Event_On_Blacklist = true;
    //      break;
    //    }
    //  }
    //}
    
    //blacklist.close();


    for(int i =0; i<12; i++){
      ADC_TOF1U[i] = ADC_TOF1[i];
      ADC_TOF1D[i] = ADC_TOF1[i+12];
    }

    for(int i=0; i<12; i++)
    {

      if(ADC_TOF1U[i]>=0 || ADC_TOF1D[i]>=0)  Good_TOF1_ADC[i] = true;
      
      if((TDC_TOF1U[i]>=TOF1U_TDC_min[i] && TDC_TOF1U[i]<=TOF1U_TDC_max[i]) ||
         (TDC_TOF1D[i]>=TOF1D_TDC_min[i] && TDC_TOF1D[i]<=TOF1D_TDC_max[i]))  Good_TOF1_TDC[i] = true;

      if(Good_TOF1_TDC[i] || Good_TOF1_ADC[i]) Good_TOF1[i] = true;

      if((ADC_TOF2AO[i]>=0 && ADC_TOF2AI[i]>=0) || (ADC_TOF2BO[i]>=0 && ADC_TOF2BI[i]>=0))  Good_TOF2_ADC[i] = true;

      if(((TDC_TOF2AO[i]>=TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<=TOF2AO_TDC_max[i])  &&
          (TDC_TOF2AI[i]>=TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<=TOF2AI_TDC_max[i])) ||
         ((TDC_TOF2BO[i]>=TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<=TOF2BO_TDC_max[i])  &&
          (TDC_TOF2BI[i]>=TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<=TOF2BI_TDC_max[i])))  Good_TOF2_TDC[i] = true;

      if(Good_TOF2_ADC[i] || Good_TOF2_TDC[i])
       Good_TOF2[i] = true;
    }


    for(int k=0; k<12; k++)
    {
      if(k!=0 && k!=11)
      {
        if((Good_TOF2[k] && Good_TOF1[k-1]) || (Good_TOF2[k] && Good_TOF1[k]) || (Good_TOF2[k] && Good_TOF1[k+1]))
        {
          Good_TOFs[k] = true;
        }
      }
    }

    if((Good_TOF2[0] && Good_TOF1[11]) || (Good_TOF2[0] && Good_TOF1[0]) || (Good_TOF2[0] && Good_TOF1[1]))  Good_TOFs[0] = true;

    if((Good_TOF2[11] && Good_TOF1[10]) || (Good_TOF2[11] && Good_TOF1[11]) || (Good_TOF2[11] && Good_TOF1[0]))  Good_TOFs[11] = true;

    for(int kk=0; kk<12; kk++){   
      if(Good_TOFs[kk]) Good_TOF_Event = true;
    }     

    //for(int i=0; i<12; i++){
      //Good_tof1[i] = Good_TOF1[i];
      //Good_tof2[i] = Good_TOF2[i];
    //}

    //if(enable_cout == 9) Good_TOF_Event = true;


    //********* GOOD MWPC EVENTS  
    count_C2X = 0;    count_C2Y = 0;
    count_C3X = 0;    count_C3Y = 0;
    count_C4X = 0;    count_C4Y = 0;
    Good_MWPC_Event = false;

        //for(int i=128; i<254; i++){
    //  if(i!=142 && i!=143 && i!=158 && i!=159 && i!=174 && i!=175 && i!=190 && i!=191 && 
    //     i!=206 && i!=207 && i!=223 && i!=224 && i!=238 && i!=239){
    //  
    //    if(MWPCADC[i]>=0)  count_C2X++; 
    //  } 
    //} 
    
    for(int i=0; i<56; i++){
      if(ADC_C2X_R[i]>0. || ADC_C2X_L[i]>0.) count_C2X++;
    }

    //for(int ii=96; ii<128; ii++){
    //  
    //  if(MWPCADC[ii]>0)  count_C2Y++;
    //}
    
    for(int ii=0; ii<16; ii++){
      if(ADC_C2Y_R[ii]>0. || ADC_C2Y_L[ii]>0.) count_C2Y++;
    }


    //for(int j=0; j<96; j++){
    //
    //  if(MWPCADC[j]>0){
    //    count_C3X++;
    //  }
    //}
    
  
    //for(int j=480; j<511; j++){
    // 
    //  if(MWPCADC[j]>0){
    //    count_C3X++;
    //  }
    //}

    for(int j=0; j<64; j++){
      if(ADC_C3X_R[j]>0. || ADC_C3X_L[j]>0.){
        count_C3X++;
      }
    }

 
    //for(int j=448; j<479; j++){
    //
    //  if(MWPCADC[j]>0) count_C3Y++;
    //}
    
    for(int jj=0; jj<16; jj++){
      if(ADC_C3Y_R[jj]>0. || ADC_C3Y_L[jj]>0.) count_C3Y++;
    }


    //for(int k=288; k<448; k++){
    //  if((k>=288 && k<=295) || (k>=304 && k<=311) || (k>=320 && k<=447)){
    //   if(MWPCADC[k]>0){
    //      count_C4X++;
    //    }
    //  }
    //}

    for(int k=0; k<72; k++){
      if(ADC_C4X_R[k]>0. || ADC_C4X_L[k]>0.) count_C4X++;
    }


    //for(int kk=256; kk<288; kk++){
    //  if(MWPCADC[kk]>0)  count_C4Y++;
    //}

    for(int kk=0; kk<16; kk++){
      if(ADC_C4Y_R[kk]>0. || ADC_C4Y_L[kk]>0.) count_C4Y++;
    }


    if(count_C2X>0 && count_C2Y>0 && count_C3X>0 && count_C3Y>0 && count_C4X>0 && count_C4Y>0)  Good_MWPC_Event = true;
    //cout << "TEST: " << count_C2X << "  " << count_C2Y << "  " << count_C3X << "  " << count_C3Y << "  " << count_C4X << "  " << count_C4Y << endl;

    //if(enable_cout == 9) Good_MWPC_Event = true;

    if(Good_TARGET_Event && Good_TOF_Event && Good_MWPC_Event && !Event_On_Blacklist)  Good_Event = true;
    //if(Good_TARGET_Event && Good_TOF_Event && Good_MWPC_Event)  Good_Event = true;



    //if(Good_Event) cout << "Event: "<< ievt << "   --  GOOD EVENT!" << endl;

    if(Good_TARGET_Event && Good_TOF_Event && !Good_MWPC_Event) continue;
    if(Good_TARGET_Event && !Good_TOF_Event && Good_MWPC_Event) continue;
    if(!Good_TARGET_Event && Good_TOF_Event && Good_MWPC_Event) continue;
    if(Good_TARGET_Event && !Good_TOF_Event && !Good_MWPC_Event) continue;
    if(!Good_TARGET_Event && Good_TOF_Event && !Good_MWPC_Event) continue;
    if(!Good_TARGET_Event && !Good_TOF_Event && Good_MWPC_Event) continue;
    if(!Good_TARGET_Event && !Good_TOF_Event && !Good_MWPC_Event) continue;
    
 
    //if(Event_On_Blacklist){
      //cout << "Event: "<< ivt << " is on the blacklist." << endl;
      //cout << " " << endl;
      //cout << ">>>>>  Please, choose another event" << endl;
      //cout << " " << endl;
      //break;
      //continue;
    //}


    if(!Good_Event) return;

    for(int j=0 ; j<128 ; j++){ 
      if(ADC_High_SFT[j]<0)      ADC_High_SFT_corr[j]=0; 
      if(ADC_High_SFT[j]>=0)     ADC_High_SFT_corr[j]=ADC_High_SFT[j]; 
    }

    for(Int_t ii=0; ii<128; ii++){
      if(ADC_High_SFT_corr[ii]>ADC_High_corr_max){
        ADC_High_corr_max=ADC_High_SFT_corr[ii];
      }
    }

    for(Int_t ii=0; ii<128; ii++){
      for (Int_t qq=0; qq<6; qq++) {
        if (tdc_le_sft[ii][qq] > SFT_TDC_min[ii] && tdc_le_sft[ii][qq] < SFT_TDC_max[ii]) has_TDC_SFT_hit[ii]++;
      }
    }

    //for(Int_t imark1=0; imark1<15; imark1++){
      //if((ADC_High_SFT_corr[imark1]!=0) && (has_TDC_SFT_hit[imark1] > 0)){  
        //h_ADC_L1_DS->Fill(imark1,ADC_High_SFT_corr[imark1]);
        //h_ADC_L1_DS->SetFillColor(2);
      //}

      //if((ADC_High_SFT_corr[imark1+64]!=0) && (has_TDC_SFT_hit[imark1+64] > 0)){  
        //h_ADC_L1_US->Fill(imark1+64,ADC_High_SFT_corr[imark1+64]);
        //h_ADC_L1_US->SetFillColor(2); 
      //}
    //}

    //for(Int_t imark2=0; imark2<15; imark2++){
      //if((ADC_High_SFT_corr[imark2+15]!=0) && (has_TDC_SFT_hit[imark2+15] > 0)){  
        //h_ADC_L2_DS->Fill(imark2+15,ADC_High_SFT_corr[imark2+15]);
        //h_ADC_L2_DS->SetFillColor(4);
      //}

      //if((ADC_High_SFT_corr[imark2+79]!=0) && (has_TDC_SFT_hit[imark2+79] > 0)){  
        //h_ADC_L2_US->Fill(imark2+79,ADC_High_SFT_corr[imark2+79]);
        //h_ADC_L2_US->SetFillColor(4);
      //}
    //}

    //for(Int_t imark3=0; imark3<2; imark3++){
      //if((ADC_High_SFT_corr[imark3+30]!=0) && (has_TDC_SFT_hit[imark3+30] > 0)){  
        //h_ADC_L3_DS->Fill(imark3+30,ADC_High_SFT_corr[imark3+30]);
        //h_ADC_L3_DS->SetFillColor(3);
      //}

      //if((ADC_High_SFT_corr[imark3+94]!=0) && (has_TDC_SFT_hit[imark3+94] > 0)){  
        //h_ADC_L3_US->Fill(imark3+94,ADC_High_SFT_corr[imark3+94]);
        //h_ADC_L3_US->SetFillColor(3);
      //}
    //}


    //for(Int_t imark3=2; imark3<17; imark3++){
      //if((ADC_High_SFT_corr[imark3+47]!=0) && (has_TDC_SFT_hit[imark3+47] > 0)){  
        //h_ADC_L3_DS->Fill(imark3+47,ADC_High_SFT_corr[imark3+47]);
        //h_ADC_L3_DS->SetFillColor(3);
      //}

      //if((ADC_High_SFT_corr[imark3+111]!=0) && (has_TDC_SFT_hit[imark3+111] > 0)){  
        //h_ADC_L3_US->Fill(imark3+111,ADC_High_SFT_corr[imark3+111]);
        //h_ADC_L3_US->SetFillColor(3);
      //}
    //}


    //for(Int_t imark4=0; imark4<17; imark4++){
      //if((ADC_High_SFT_corr[imark4+32]!=0) && (has_TDC_SFT_hit[imark4+32] > 0)){  
        //h_ADC_L4_DS->Fill(imark4+32,ADC_High_SFT_corr[imark4+32]);
        //h_ADC_L4_DS->SetFillColor(1);
      //}
  
      //if((ADC_High_SFT_corr[imark4+96]!=0) && (has_TDC_SFT_hit[imark4+96] > 0)){  
        //h_ADC_L4_US->Fill(imark4+96,ADC_High_SFT_corr[imark4+96]);
        //h_ADC_L4_US->SetFillColor(1);
      //}
    //}

  //////// DETERMINE FIBER WITH HIGHEST AND SECOND HIGHEST LOW GAIN AMPLITUDE


  Int_t TDC_average = -1;

  int max_index_all[256];
  int max_ADC_all[256];
  int max_index_flag;
  
  for(int i = 0; i<256; i++){
    max_ADC_all[i] = -100000000;
    max_index_all[i] = -1;
  }
  
  
  for(int j = 0; j<256; j++){
    for(int i=0; i<256; i++){
    max_index_flag = 0;  
  
    for(int k = 0; k<256; k++){
      if(i == max_index_all[k]){
      max_index_flag = 1;
      break;
      }
    }  
      if (max_index_flag == 1) continue;   
    
      else {
        if(ADC_Low_TARGET[i]>max_ADC_all[j]) {
          max_index_all[j] = i;
          max_ADC_all[j] = ADC_Low_TARGET[i];
        }
      }
    } 
  } 

  //Calculate TDC Average
  int TDC_LG_max = -1;
  int TDC_LG_max2 = -1;
  //int index_max1=0;
  //int index_max2=0;

  
  for(int i = 0; i<256; i++){
  for(int j=0; j<6; j++){  
    if(TDC_LG_max == -1){ 
      if(tdc_le_target[max_index_all[i]][j] >= TDC_Thr_min && tdc_le_target[max_index_all[i]][j]  <= TDC_Thr_max)
        TDC_LG_max = tdc_le_target[max_index_all[i]][j];
        //index_max1 = max_index_all[i];
    }
    else if(TDC_LG_max2 == -1){
      if(tdc_le_target[max_index_all[i]][j] >= TDC_Thr_min && tdc_le_target[max_index_all[i]][j] <= TDC_Thr_max)
        TDC_LG_max2 = tdc_le_target[max_index_all[i]][j];
        //index_max2 = max_index_all[i];
      }
  }     
  } 
  
  if(abs(TDC_LG_max - TDC_LG_max2) <= T_limit) TDC_average = (TDC_LG_max + TDC_LG_max2)/2;
  else if(TDC_LG_max > TDC_LG_max2) TDC_average = TDC_LG_max;
  else if(TDC_LG_max < TDC_LG_max2) TDC_average = TDC_LG_max2;
 
  //printf("                    Fiber  ADC-Ped   TDC\n");
  //printf("LG ADC 1st Max.  :  %4d    %4d    %4d\n",index_max1, ADC_Low_TARGET_ped[index_max1]+, TDC_LG_max);
  //printf("LG ADC 2nd Max.  :  %4d    %4d    %4d\n",index_max2, ADC_Low_TARGET_ped[index_max2], TDC_LG_max2);
  //printf("LG ADC 1st Max.  :  %4d    %4d    %4d\n",index_max1, ADC_Low_TARGET_ped[index_max1], TDC_LG_max);
  //printf("LG ADC 2nd Max.  :  %4d    %4d    %4d\n",index_max2, ADC_Low_TARGET_ped[index_max2], TDC_LG_max2);
  //printf("TDC Average  :  %4d\n",TDC_average);

  //cout << " " << endl;

  int kaon_TDC_min = TDC_average + TDC_Avg_Offset_min;
  int kaon_TDC_max = TDC_average + TDC_Avg_Offset_max;
  int TDC_min_Kstop = TDC_average + TDC_Kstop_Avg_Offset_min;
  int TDC_max_Kstop = TDC_average + TDC_Kstop_Avg_Offset_max;
  

  //Int_t TDC_min_TARGET = TARGET_TDC_min[0];
  //Int_t TDC_max_TARGET = TARGET_TDC_max[0];
  TDC_min_TARGET = kaon_TDC_min;
  TDC_max_TARGET = kaon_TDC_max;



    for(int i=0; i<256; i++){
      has_TDC_hit[i] = false;
      has_TDC_hit_Kstop[i] = false;
    }

    for(Int_t i=0; i<256; i++){
      for (Int_t k=0; k<4; k++) {
        if ((tdc_le_target[i][k]>=TDC_min_Kstop) && (tdc_le_target[i][k]<=TDC_max_Kstop)) has_TDC_hit_Kstop[i] = true;
        //if ((tdc_le_target[i][k]>=TARGET_TDC_min[i]) && (tdc_le_target[i][k]<=TARGET_TDC_max[i])) has_TDC_hit[i] = true;
        if ((tdc_le_target[i][k]>=kaon_TDC_min) && (tdc_le_target[i][k]<=kaon_TDC_max)){
          has_TDC_hit[i] = true;
        }
      }
    }


    sprintf(ch_ADC_cut_TARGET,"(ADC offset = %d)",TARGET_ADC_Thr_HG_Offset);
    sprintf(ch_ADC_and_TDC_cut,"(ADC offset = %d | %d #leq TDC #leq %d)",TARGET_ADC_Thr_HG_Offset,TDC_min_TARGET,TDC_max_TARGET);
    sprintf(ch_ADC_and_TDC_cut_Kstop,"(ADC: HG #geq %d, LG #geq %d | %d #leq TDC K Stop #leq %d)",HG_KAON,LG_KAON,TDC_min_Kstop,TDC_max_Kstop);

    Angle_ADC_cut = 0;


    x_inc = 0;
    y_inc = 0;
    hit_count = 0;
    count = 0; 

    max_index = 0;
    max_index2 = 0;
    max_index3 = 0;
    max_index4 = 0;

    max_ADC = -100000000;
    max_ADC2 = -100000000;
    max_ADC3 = -100000000;
    max_ADC4 = -100000000;  

    for(Int_t q=0; q<256; q++){
      if(ADC_High_TARGET[q]>max_ADC) {
        max_index = q;
        max_ADC = ADC_High_TARGET[q];
      }
    }

    for(Int_t q=0; q<256; q++){
      if (q == max_index) continue;   
      else {
        if(ADC_High_TARGET[q]>max_ADC2) {
          max_index2 = q;
          max_ADC2 = ADC_High_TARGET[q];
        }
      }
    }

    for(Int_t q=0; q<256; q++){
      if ((q == max_index) || (q == max_index2)) continue;    
      else {
        if(ADC_High_TARGET[q]>max_ADC3) {
          max_index3 = q;
          max_ADC3 = ADC_High_TARGET[q];
        }
      }
    }

    for(Int_t q=0; q<256; q++){
      if ((q == max_index) || (q == max_index2) || (q == max_index3)) continue;   
      else {
        if(ADC_High_TARGET[q]>max_ADC4) {
          max_index4 = q;
          max_ADC4 = ADC_High_TARGET[q];
        }
      }
    }

    x_cent = Xloc[max_index];   y_cent = Yloc[max_index];
  

    for(Int_t j=0; j<256; j++){
      hyp[j] = -1;
      hyp[j] = sqrt(pow(x_cent - Xloc[j],2) + pow(y_cent - Yloc[j],2));
    }

    for(int i=0; i<12; i++){
      has_ADC_TOF1_hit[i] = false;
      has_TDC_TOF1_hit[i] = false;
      has_ADC_TOF2_hit[i] = false;
      has_TDC_TOF2_hit[i] = false;
      //has_both_ADC_TOF1_hit[i] = false;
    }


    ///Set TOF2 Lines

    for(int i = 0; i < 12; i++){
      if ((ADC_TOF2AO[i]>0 && ADC_TOF2AI[i]>0) || (ADC_TOF2BO[i]>0 && ADC_TOF2BI[i]>0)) {has_ADC_TOF2_hit[i]=true;}
      if ((((TDC_TOF2AO[i]>TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<TOF2AO_TDC_max[i]) && (TDC_TOF2AI[i]>TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<TOF2AI_TDC_max[i]))) 
      ||  (((TDC_TOF2BO[i]>TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<TOF2BO_TDC_max[i]) && (TDC_TOF2BI[i]>TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<TOF2BI_TDC_max[i])))) {has_TDC_TOF2_hit[i]=true;}
    }

  
    //for(int i=0; i<12; i++){
      //if(ADC_TOF1U[i]>0 && ADC_TOF1D[i]>0)
      //has_both_ADC_TOF1_hit[i] = true;  
    //}
  
  
    ///Set TOF1 Lines

    for(int i = 0; i < 12; i++){
      if (ADC_TOF1U[i]>0 || ADC_TOF1D[i]>0) {has_ADC_TOF1_hit[i] = true;}
      if ((TDC_TOF1U[i]>TOF1U_TDC_min[i] && TDC_TOF1U[i]<TOF1U_TDC_max[i]) || (TDC_TOF1D[i]>TOF1D_TDC_min[i] && TDC_TOF1D[i]<TOF1D_TDC_max[i])) {has_TDC_TOF1_hit[i] = true;}
    }

  
    for(int i=0; i<12; i++){
      //gap_hit[i] = 0;
      ADC_TOF1_hit[i] = 0;
      ADCTDC_TOF1_hit[i] = 0;
      ADC_TOF2_hit[i] = 0;
      ADCTDC_TOF2_hit[i] = 0;
    }
 
    //for(int j=0; j<12; j++){
    //  for(int k=0; k<8; k++){
    //    if (ADC_High_TARGET[channel[j][k]] > 0 && has_TDC_hit[channel[j][k]]) gap_hit[j] = 1;
    //  }
    //}


    for (int k=0; k<12; k++) {
      if (has_ADC_TOF1_hit[k]) {
        if (has_TDC_TOF1_hit[k]) {ADCTDC_TOF1_hit[k]++;}
        else {ADC_TOF1_hit[k]++;}
      }
      if (has_ADC_TOF2_hit[k]) {
        if (has_TDC_TOF2_hit[k]) {ADCTDC_TOF2_hit[k]++;}
        else {ADC_TOF2_hit[k]++;}
      }
    }
  
    selected_TOF2 = 0;
  
    // Determine which TOF2 is hit
    for(int i = 0; i<12; i++){
      if(has_TDC_TOF2_hit[i] && has_ADC_TOF2_hit[i])
        selected_TOF2 = i + 1;  
    }

    if(selected_TOF2 == 0){
      for(int i = 0; i<12; i++){
        if(has_TDC_TOF2_hit[i] || has_ADC_TOF2_hit[i])
          selected_TOF2 = i+1;
      }
    }

    for(int i=0; i<12; i++) gap_counter[i] = 0;

 
    //// GAP SCORING !
  
    scoring_type = 2;           // scoring_type = 1  --->  Oscar's Method
                                // scoring_type = 2  --->  Sebastien's Method  TOF1[i]  TOF2[i-1], TOF2[i], TOF2[i+1]
                                // scoring_type = 3  --->  Sebastien's Method  TOF2[i]  TOF1[i-1], TOF1[i], TOF1[i+1]
    //Sebastien's Method
    //cout << "BIG TEST  !" << endl;

    if(scoring_type==2){
      for(int i=0; i<12; i++){
        if(ADC_TOF1U[i]>=0){
          gap_counter[i]++;
        }
      }

      for(int i=0; i<12; i++){
        if(ADC_TOF1D[i]>=0){
          gap_counter[i]++;
        }
      }

      for(int i=0; i<12; i++){
        if(TDC_TOF1U[i]>=TOF1U_TDC_min[i] && TDC_TOF1U[i]<=TOF1U_TDC_max[i]){
          gap_counter[i]++; 
        }
      }

      for(int i=0; i<12; i++){
        if(TDC_TOF1D[i]>=TOF1D_TDC_min[i] && TDC_TOF1D[i]<=TOF1D_TDC_max[i]){
          gap_counter[i]++;
        }
      }


      for(int i=0; i<12; i++){
        if(i!=0 && i!=11 && gap_counter[i]>0){

          if(ADC_TOF2AO[i]>=0) gap_counter[i]++;
          if(ADC_TOF2AO[i-1]>=0) gap_counter[i]++;
          if(ADC_TOF2AO[i+1]>=0) gap_counter[i]++;
 
          if(ADC_TOF2AI[i]>=0) gap_counter[i]++;
          if(ADC_TOF2AI[i-1]>=0) gap_counter[i]++;
          if(ADC_TOF2AI[i+1]>=0) gap_counter[i]++;   

          if(ADC_TOF2BO[i]>=0) gap_counter[i]++;
          if(ADC_TOF2BO[i-1]>=0) gap_counter[i]++;
          if(ADC_TOF2BO[i+1]>=0) gap_counter[i]++;   

          if(ADC_TOF2BI[i]>=0) gap_counter[i]++;
          if(ADC_TOF2BI[i-1]>=0) gap_counter[i]++;
          if(ADC_TOF2BI[i+1]>=0) gap_counter[i]++;

          if(TDC_TOF2AO[i]>=TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<=TOF2AO_TDC_max[i]) gap_counter[i]++;   
          if(TDC_TOF2AO[i-1]>=TOF2AO_TDC_min[i-1] && TDC_TOF2AO[i-1]<=TOF2AO_TDC_max[i-1]) gap_counter[i]++;   
          if(TDC_TOF2AO[i+1]>=TOF2AO_TDC_min[i+1] && TDC_TOF2AO[i+1]<=TOF2AO_TDC_max[i+1]) gap_counter[i]++;   

          if(TDC_TOF2AI[i]>=TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<=TOF2AI_TDC_max[i]) gap_counter[i]++;   
          if(TDC_TOF2AI[i-1]>=TOF2AI_TDC_min[i-1] && TDC_TOF2AI[i-1]<=TOF2AI_TDC_max[i-1]) gap_counter[i]++;   
          if(TDC_TOF2AI[i+1]>=TOF2AI_TDC_min[i+1] && TDC_TOF2AI[i+1]<=TOF2AI_TDC_max[i+1]) gap_counter[i]++;   

          if(TDC_TOF2BO[i]>=TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<=TOF2BO_TDC_max[i]) gap_counter[i]++;   
          if(TDC_TOF2BO[i-1]>=TOF2BO_TDC_min[i-1] && TDC_TOF2BO[i-1]<=TOF2BO_TDC_max[i-1]) gap_counter[i]++;   
          if(TDC_TOF2BO[i+1]>=TOF2BO_TDC_min[i+1] && TDC_TOF2BO[i+1]<=TOF2BO_TDC_max[i+1]) gap_counter[i]++;   

          if(TDC_TOF2BI[i]>=TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<=TOF2BI_TDC_max[i]) gap_counter[i]++;   
          if(TDC_TOF2BI[i-1]>=TOF2BI_TDC_min[i-1] && TDC_TOF2BI[i-1]<=TOF2BI_TDC_max[i-1]) gap_counter[i]++;   
          if(TDC_TOF2BI[i+1]>=TOF2BI_TDC_min[i+1] && TDC_TOF2BI[i+1]<=TOF2BI_TDC_max[i+1]) gap_counter[i]++;   
        }
      }

      if(gap_counter[0]>0){
        if(ADC_TOF2AO[0]>=0) gap_counter[0]++;
        if(ADC_TOF2AO[1]>=0) gap_counter[0]++;
        if(ADC_TOF2AO[11]>=0) gap_counter[0]++;

        if(ADC_TOF2AI[0]>=0) gap_counter[0]++;
        if(ADC_TOF2AI[1]>=0) gap_counter[0]++;
        if(ADC_TOF2AI[11]>=0) gap_counter[0]++;   

        if(ADC_TOF2BO[0]>=0) gap_counter[0]++;
        if(ADC_TOF2BO[1]>=0) gap_counter[0]++;
        if(ADC_TOF2BO[11]>=0) gap_counter[0]++;   

        if(ADC_TOF2BI[0]>=0) gap_counter[0]++;
        if(ADC_TOF2BI[1]>=0) gap_counter[0]++;
        if(ADC_TOF2BI[11]>=0) gap_counter[0]++;

        if(TDC_TOF2AO[0]>=TOF2AO_TDC_min[0] && TDC_TOF2AO[0]<=TOF2AO_TDC_max[0]) gap_counter[0]++;   
        if(TDC_TOF2AO[1]>=TOF2AO_TDC_min[1] && TDC_TOF2AO[1]<=TOF2AO_TDC_max[1]) gap_counter[0]++;   
        if(TDC_TOF2AO[11]>=TOF2AO_TDC_min[11] && TDC_TOF2AO[11]<=TOF2AO_TDC_max[11]) gap_counter[0]++;   

        if(TDC_TOF2AI[0]>=TOF2AI_TDC_min[0] && TDC_TOF2AI[0]<=TOF2AI_TDC_max[0]) gap_counter[0]++;   
        if(TDC_TOF2AI[1]>=TOF2AI_TDC_min[1] && TDC_TOF2AI[1]<=TOF2AI_TDC_max[1]) gap_counter[0]++;   
        if(TDC_TOF2AI[11]>=TOF2AI_TDC_min[11] && TDC_TOF2AI[11]<=TOF2AI_TDC_max[11]) gap_counter[0]++;   

        if(TDC_TOF2BO[0]>=TOF2BO_TDC_min[0] && TDC_TOF2BO[0]<=TOF2BO_TDC_max[0]) gap_counter[0]++;   
        if(TDC_TOF2BO[1]>=TOF2BO_TDC_min[1] && TDC_TOF2BO[1]<=TOF2BO_TDC_max[1]) gap_counter[0]++;   
        if(TDC_TOF2BO[11]>=TOF2BO_TDC_min[11] && TDC_TOF2BO[11]<=TOF2BO_TDC_max[11]) gap_counter[0]++;   

        if(TDC_TOF2BI[0]>=TOF2BI_TDC_min[0] && TDC_TOF2BI[0]<=TOF2BI_TDC_max[0]) gap_counter[0]++;   
        if(TDC_TOF2BI[1]>=TOF2BI_TDC_min[1] && TDC_TOF2BI[1]<=TOF2BI_TDC_max[1]) gap_counter[0]++;   
        if(TDC_TOF2BI[11]>=TOF2BI_TDC_min[11] && TDC_TOF2BI[11]<=TOF2BI_TDC_max[11]) gap_counter[0]++;   
      }

      if(gap_counter[11]>0){
        if(ADC_TOF2AO[11]>=0) gap_counter[11]++;
        if(ADC_TOF2AO[10]>=0) gap_counter[11]++;
        if(ADC_TOF2AO[0]>=0) gap_counter[11]++;

        if(ADC_TOF2AI[11]>=0) gap_counter[11]++;
        if(ADC_TOF2AI[10]>=0) gap_counter[11]++;
        if(ADC_TOF2AI[0]>=0) gap_counter[11]++;   

        if(ADC_TOF2BO[11]>=0) gap_counter[11]++;
        if(ADC_TOF2BO[10]>=0) gap_counter[11]++;
        if(ADC_TOF2BO[0]>=0) gap_counter[11]++;   

        if(ADC_TOF2BI[11]>=0) gap_counter[11]++;
        if(ADC_TOF2BI[10]>=0) gap_counter[11]++;
        if(ADC_TOF2BI[0]>=0) gap_counter[11]++;

        if(TDC_TOF2AO[11]>=TOF2AO_TDC_min[11] && TDC_TOF2AO[11]<=TOF2AO_TDC_max[11]) gap_counter[11]++;   
        if(TDC_TOF2AO[10]>=TOF2AO_TDC_min[10] && TDC_TOF2AO[10]<=TOF2AO_TDC_max[10]) gap_counter[11]++;   
        if(TDC_TOF2AO[0]>=TOF2AO_TDC_min[0] && TDC_TOF2AO[0]<=TOF2AO_TDC_max[0]) gap_counter[11]++;   

        if(TDC_TOF2AI[11]>=TOF2AI_TDC_min[11] && TDC_TOF2AI[11]<=TOF2AI_TDC_max[11]) gap_counter[11]++;   
        if(TDC_TOF2AI[10]>=TOF2AI_TDC_min[10] && TDC_TOF2AI[10]<=TOF2AI_TDC_max[10]) gap_counter[11]++;   
        if(TDC_TOF2AI[0]>=TOF2AI_TDC_min[0] && TDC_TOF2AI[0]<=TOF2AI_TDC_max[0]) gap_counter[11]++;   

        if(TDC_TOF2BO[11]>=TOF2BO_TDC_min[11] && TDC_TOF2BO[11]<=TOF2BO_TDC_max[11]) gap_counter[11]++;   
        if(TDC_TOF2BO[10]>=TOF2BO_TDC_min[10] && TDC_TOF2BO[10]<=TOF2BO_TDC_max[10]) gap_counter[11]++;   
        if(TDC_TOF2BO[0]>=TOF2BO_TDC_min[0] && TDC_TOF2BO[0]<=TOF2BO_TDC_max[0]) gap_counter[11]++;   

        if(TDC_TOF2BI[11]>=TOF2BI_TDC_min[11] && TDC_TOF2BI[11]<=TOF2BI_TDC_max[11]) gap_counter[11]++;   
        if(TDC_TOF2BI[10]>=TOF2BI_TDC_min[10] && TDC_TOF2BI[10]<=TOF2BI_TDC_max[10]) gap_counter[11]++;   
        if(TDC_TOF2BI[0]>=TOF2BI_TDC_min[0] && TDC_TOF2BI[0]<=TOF2BI_TDC_max[0]) gap_counter[11]++;   
      } 
    }
  

    high_gap_hit = 0;
    gap_to_fit = 0;
    score_max = 0;

    tof1_ties.clear();
    //tdc_tof1u_min = -10000;
    //tdc_tof1d_min = -10000;

    for (int k=0; k<12; k++) {
      //if((gap_counter[k]>=high_gap_hit) && (has_TDC_TOF1_hit[k] || has_ADC_TOF1_hit[k])){
      if(gap_counter[k]>=high_gap_hit){
        if(gap_counter[k] == high_gap_hit) tof1_ties.push_back(k); 
        else{
          tof1_ties.clear();
          tof1_ties.push_back(k);
        }
        high_gap_hit = gap_counter[k];
        //gap_to_fit = k+1;
      }
    }
    
    if(tof1_ties.size() > 1){  // only when there are more than 2 elements in tof1_ties!
      for(unsigned int j=0; j<tof1_ties.size(); j++){
        //cout << "POE:  " << j << "   " << tof1_ties[j] <<endl; 
        for(int jj=0; jj<8; jj++){
          //cout << "INT_CH:  " << channel[tof1_ties[j]][jj] << endl;
          if(ADC_High_TARGET[channel[tof1_ties[j]][jj]] > 0){
            //cout << "ALB:  " << channel[tof1_ties[j]][jj] << "    " << ADC_High_TARGET[channel[tof1_ties[j]][jj]] << endl;
            gap_counter[tof1_ties[j]]++;
          }
          if(ADC_Low_TARGET[channel[tof1_ties[j]][jj]] > 0){
            //cout << "ALB2:  " << channel[tof1_ties[j]][jj] << "    " << ADC_Low_TARGET[channel[tof1_ties[j]][jj]] << endl;
            gap_counter[tof1_ties[j]]++;
          }
        }
        //cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
      }
    }

    for(int k=0; k<12; k++){
      if(gap_counter[k] > score_max){
        score_max = gap_counter[k];
        gap_to_fit = k+1;
      }
    } 


    //if(!tof1_ties.empty()){
    //  for(vector<int>::iterator it = tof1_ties.begin(); it != tof1_ties.end(); it++){
    //    if(TDC_TOF1U[*it] >= tdc_tof1u_min && TDC_TOF1D[*it] >= tdc_tof1d_min){
    //      gap_to_fit = *it + 1;
    //      tdc_tof1u_min = TDC_TOF1U[*it];
    //      tdc_tof1d_min = TDC_TOF1D[*it];
    //    }
    //  }
    //}


    ////////////////////// ADDITIONAL


    // Determine the k_stop bars  
    for(int i=0; i<256; i++) k_stop_bar[i] = false;
    good_k_stop_bars.clear();
  
    X_BAR = 0;
    Y_BAR = 0;
  
    
    total_energy = 0.0;
    X_weights = 0.0;
    Y_weights = 0.0;   
    
    for(int i = 0; i<256; i++){
      //if(ADC_High_TARGET_ped[i] > HG_KAON && ADC_Low_TARGET_ped[i] > LG_KAON && has_TDC_hit_Kstop[i]){
      if(ADC_High_TARGET[i] > HG_KAON && ADC_Low_TARGET[i] > LG_KAON && has_TDC_hit_Kstop[i]){
      good_k_stop_bars.push_back(i);
      k_stop_bar[i] = true;
      }
    }    

  
    // Determine if a hit target has any hit neighbours
    
    for(int i=0; i<256; i++) TARGET_High_has_neighbours[i] = false;
    
    //for(int i=0; i<256; i++){
    //  if(k_stop_bar[i])  cout << "FFFF:  " << i << endl;
    //}

    // Determine if a hit target has any hit neighbours
    // bool TARGET_High_has_neighbours[256] = {false};
    for(int i = 0; i<256; i++){
      for(int j=0; j<8; j++){
        if((TARGET_neighbours[i][j] != -1 && ADC_High_TARGET[TARGET_neighbours[i][j]] >= Angle_ADC_cut && has_TDC_hit[TARGET_neighbours[i][j]] && !k_stop_bar[TARGET_neighbours[i][j]]) ||
          (TARGET_neighbours[i][j] != -1 && ADC_High_TARGET[TARGET_neighbours[i][j]] <0 && ADC_Low_TARGET[TARGET_neighbours[i][j]]>=0 && !k_stop_bar[TARGET_neighbours[i][j]])){
           TARGET_High_has_neighbours[i] = true;
             break;
          }
      }
    }


    for(Int_t i=0; i<256; i++){
      if(TARGET_High_has_neighbours[i]){
        //cout << "TEST:  " << i << endl;
        if(ADC_High_TARGET[i]>=Angle_ADC_cut && has_TDC_hit[i]){
          h_target_TDC->Fill(Xloc[i],Yloc[i]);
          h_target_TDC_copy->Fill(Xloc[i],Yloc[i]);
        
          if(!k_stop_bar[i]){
            h_TDC_Gap_Fibers->Fill(Xloc[i],Yloc[i]);
          }
          count++;
          //cout << "DEBUGGG:  " << h_TDC_Gap_Fibers->GetEntries() << endl;

        //cout << "DDEBUG: " << h_TDC_Gap_Fibers->GetEntries() << endl;

          if(IsIn(i,channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],   // Additional weight for fibers close to the edge if hit
                    channel[gap_to_fit-1][2], channel[gap_to_fit-1][3], 
                    channel[gap_to_fit-1][4], channel[gap_to_fit-1][5], 
                    channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){
            if(!k_stop_bar[i]){         
              h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
              h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
              }
          
          
            h_TDC_selected2->Fill(Xloc[i], Yloc[i]);
            h_TDC_selected2->Fill(Xloc[i], Yloc[i]);
  
            h_TDC_selected->Fill(Xloc[i], Yloc[i]);
            h_TDC_selected->Fill(Xloc[i], Yloc[i]);
          }
        }

        if(ADC_High_TARGET[i]>=Angle_ADC_cut){
          h_target_ADC->Fill(Xloc[i],Yloc[i]);
          x_inc = x_inc + Xloc[i];
          y_inc = y_inc + Yloc[i];
          hit_count++;
        }
    
        if(ADC_High_TARGET[i]<0 && ADC_Low_TARGET[i]>=0){
          if(!k_stop_bar[i])  
            h_TDC_Gap_Fibers->Fill(Xloc[i],Yloc[i]);
            h_TARGET_LG_Blue->Fill(Xloc[i],Yloc[i]);
            count++;
          }
      }
    }
    ///////// DETERMINE FIBER WITH HIGHEST AND SECOND HIGHEST LOW GAIN AMPLITUDE
    TDC_average = -1;

    for(int i=0; i<256; i++){
      max_index_all[i] = 999;
      max_ADC_all[i] = 999;
    }

    max_index_flag = 999;
  
    for(int i = 0; i<256; i++){
      max_ADC_all[i] = -100000000;
      max_index_all[i] = -1;
    }
  
  
    for(int j = 0; j<256; j++){
      for(int i=0; i<256; i++){
      max_index_flag = 0;  
    
      for(int k = 0; k<256; k++){
        if(i == max_index_all[k]){
        max_index_flag = 1;
        break;
        }
      }  
        if (max_index_flag == 1) continue;   
    
        else {
          if(ADC_Low_TARGET[i]>max_ADC_all[j]) {
            max_index_all[j] = i;
            max_ADC_all[j] = ADC_Low_TARGET[i];
          }
        }
      } 
    } 

    //Calculate TDC Average
    TDC_LG_max = -1;
    TDC_LG_max2 = -1;
    //index_max1 = 0;
    //index_max2 = 0;
  
    for(int i = 0; i<256; i++){
      for(int j=0; j<6; j++){  
        if(TDC_LG_max == -1){ 
          if(tdc_le_target[max_index_all[i]][j] >= TDC_Thr_min && tdc_le_target[max_index_all[i]][j]  <= TDC_Thr_max)
            TDC_LG_max = tdc_le_target[max_index_all[i]][j];
            //index_max1 = max_index_all[i];
        }
        else if(TDC_LG_max2 == -1){
          if(tdc_le_target[max_index_all[i]][j] >= TDC_Thr_min && tdc_le_target[max_index_all[i]][j] <= TDC_Thr_max)
            TDC_LG_max2 = tdc_le_target[max_index_all[i]][j];
            //index_max2 = max_index_all[i];
        }
      }     
    } 
  
  
    if(abs(TDC_LG_max - TDC_LG_max2) <= T_limit) TDC_average = (TDC_LG_max + TDC_LG_max2)/2;
    else if(TDC_LG_max > TDC_LG_max2) TDC_average = TDC_LG_max;
    else if(TDC_LG_max < TDC_LG_max2) TDC_average = TDC_LG_max2;
 

    for(int k=0; k<256; k++){
      for(int l=0; l<6; l++){
        TDC_LE_TARGET_corrected[k][l] = 9999;
        if(tdc_le_target[k][l]>-1)  TDC_LE_TARGET_corrected[k][l] = tdc_le_target[k][l] - TDC_average;
      }
    }

    for(int k=0; k<256; k++){
      for(int l=0; l<6; l++){
        if(TDC_LE_TARGET_corrected[k][l] == 9999) sprintf(TDC_LE_TARGET_corr[k][l],"----");
        if(TDC_LE_TARGET_corrected[k][l] != 9999) sprintf(TDC_LE_TARGET_corr[k][l],"%4d",TDC_LE_TARGET_corrected[k][l]);
      }    
    }  


    //// Don't do anything if the event has less than n_hits hits in the TARGET

    //if (enable_cout == 9) count = n_hit+1;
    if(count<n_hit){
      gROOT->Reset();
      gROOT->Clear();
      //cout << " >>>>  Event "<< ivt << " has fewer than "<< n_hit << " high gain hits in the TARGET (outliers and K stop removed)!" << endl;
      //cout << " >>>>  Please, choose another event" << endl;
      //cout << " " << endl;
      continue;
      //return; 
    }
  
    /// Print Values

    // Arrays for finding mwpc clustering

    for(int i=0; i<16; i++){
      C2Y_L[i] = 0;
      C2Y_R[i] = 0;
      C3Y_L[i] = 0;
      C3Y_R[i] = 0;
      C4Y_L[i] = 0;
      C4Y_R[i] = 0;
    }

    for(int i=0; i<56; i++){
      C2X_L[i] = 0;
      C2X_R[i] = 0;
    }

    for(int i=0; i<64; i++){
      C3X_L[i] = 0;
      C3X_R[i] = 0;
    }

    for(int i=0; i<72; i++){
      C4X_L[i] = 0;
      C4X_R[i] = 0;
    }

    //C2X_centroid = 0.0;


    //C2 Counters

    for(int i=0; i<56; i++){
      if (Good_Event && ADC_C2X_L[i]>0){
        C2X_L[i] = ADC_C2X_L[i];
      }
      if (Good_Event && ADC_C2X_R[i]>0){
        C2X_R[i] = ADC_C2X_R[i];
      }
    }

    for(int i=0; i<64; i++){
      if (Good_Event && ADC_C3X_L[i]>0){
        C3X_L[i] = ADC_C3X_L[i];
      }
      if (Good_Event && ADC_C3X_R[i]>0){
        C3X_R[i] = ADC_C3X_R[i];
      }
    }   

    for(int i=0; i<72; i++){
      if (Good_Event && ADC_C4X_L[i]>0){
        C4X_L[i] = ADC_C4X_L[i];
      }
      if (Good_Event && ADC_C4X_R[i]>0){
        C4X_R[i] = ADC_C4X_R[i];
      }
    } 

    for(int i=0; i<16; i++){
      if (Good_Event && ADC_C2Y_L[i]>0) {
        C2Y_L[i] = ADC_C2Y_L[i];
        }
      if (Good_Event && ADC_C2Y_R[i]>0) {
        C2Y_R[i] = ADC_C2Y_R[i];
        }
      if (Good_Event && ADC_C3Y_L[i]>0) {
        C3Y_L[i] = ADC_C3Y_L[i];
        }
      if (Good_Event && ADC_C3Y_R[i]>0) {
        C3Y_R[i] = ADC_C3Y_R[i];
        }
      if (Good_Event && ADC_C4Y_L[i]>0) {
        C4Y_L[i] = ADC_C4Y_L[i];
        }
      if (Good_Event && ADC_C4Y_R[i]>0) {
        C4Y_R[i] = ADC_C4Y_R[i];
        }                
    }        
   

    // Find clustering of mwpcs
    ////////////////////////////
    first_cluster = true;
    cluster_spacing = 0;
    cluster_length_count = 0;

    C2X_clusters = 0;
    C2Y_clusters = 0;
    C3X_clusters = 0;
    C3Y_clusters = 0;
    C4X_clusters = 0;
    C4Y_clusters = 0;

    C2X_cluster_index.clear(); // Hold the starting index of each cluster.
    C2X_cluster_length.clear(); // Hold the length of each cluster.
    C2Y_cluster_index.clear();
    C2Y_cluster_length.clear();

    C3X_cluster_index.clear();
    C3X_cluster_length.clear();
    C3Y_cluster_index.clear();
    C3Y_cluster_length.clear();    

    C4X_cluster_index.clear();
    C4X_cluster_length.clear();  
    C4Y_cluster_index.clear();
    C4Y_cluster_length.clear();  

    if(selected_TOF2 > 6){  // LEFT
      // count C2X clusters
      first_cluster = true;
      cluster_spacing = 0;
      cluster_length_count = 0;

      for(int i = 0; i<56; i++){
        if(C2X_L[i] > 0 && first_cluster){
          cluster_spacing = MWPC_cluster_separation + 1;
          first_cluster = false;
        }

        if(C2X_L[i] > 0){
          if(cluster_spacing > MWPC_cluster_separation){
            C2X_clusters++;
            C2X_cluster_index.push_back(i);
          }
          cluster_length_count++;
          cluster_spacing = 0;

          if(i+2 < 56 && C2X_L[i+1] <= 0 && C2X_L[i+2] <= 0){
            C2X_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 2 == 56 && C2X_L[i+1] <= 0){
            C2X_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 1 == 56){
            C2X_cluster_length.push_back(cluster_length_count);
          }

        }
        else{
          cluster_spacing++;
          if(i != 0 && C2X_L[i-1] <= 0){
            cluster_length_count = 0;
          }
          else if ( i != 0){
            cluster_length_count++;
          }
        }
      }  


    
      // count C2Y clusters
      first_cluster = true;
      cluster_spacing = 0;
      cluster_length_count = 0;

      for(int i = 0; i<16; i++){
        if(C2Y_L[i] > 0 && first_cluster){
          cluster_spacing = MWPC_cluster_separation + 1;
          first_cluster = false;
        }

        if(C2Y_L[i] > 0){
          if(cluster_spacing > MWPC_cluster_separation){
            C2Y_clusters++;
            C2Y_cluster_index.push_back(i);
          }
          cluster_length_count++;
          cluster_spacing = 0;

          if(i+2 < 16 && C2Y_L[i+1] <= 0 && C2Y_L[i+2] <= 0){
            C2Y_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 2 == 16 && C2Y_L[i+1] <= 0){
            C2Y_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 1 == 16){
            C2Y_cluster_length.push_back(cluster_length_count);
          } 

        }
        else{
          cluster_spacing++;
          if(i != 0 && C2Y_L[i-1] <= 0){
            cluster_length_count = 0;
          }
          else if ( i != 0){
            cluster_length_count++;
          }
        }
      }      

      // count C3X clusters
      first_cluster = true;
      cluster_spacing = 0;
      cluster_length_count = 0;

      for(int i = 0; i<64; i++){
        if(C3X_L[i] > 0 && first_cluster){
          cluster_spacing = MWPC_cluster_separation + 1;
          first_cluster = false;
        }

        if(C3X_L[i] > 0){
          if(cluster_spacing > MWPC_cluster_separation){
            C3X_clusters++;
            C3X_cluster_index.push_back(i);
          }
          cluster_length_count++;
          cluster_spacing = 0;

          if(i+2 < 64 && C3X_L[i+1] <= 0 && C3X_L[i+2] <= 0){
            C3X_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 2 == 64 && C3X_L[i+1] <= 0){
            C3X_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 1 == 64){
            C3X_cluster_length.push_back(cluster_length_count);
          }

        }
        else{
          cluster_spacing++;
          if(i != 0 && C3X_L[i-1] <= 0){
            cluster_length_count = 0;
          }
          else if ( i != 0){
            cluster_length_count++;
          }
        }
      }  

      // count C3Y clusters
      first_cluster = true;
      cluster_spacing = 0;
      cluster_length_count = 0;

      for(int i = 0; i<16; i++){
        if(C3Y_L[i] > 0 && first_cluster){
          cluster_spacing = MWPC_cluster_separation + 1;
          first_cluster = false;
        }

        if(C3Y_L[i] > 0){
          if(cluster_spacing > MWPC_cluster_separation){
            C3Y_clusters++;
            C3Y_cluster_index.push_back(i);
          }
          cluster_length_count++;
          cluster_spacing = 0;

          if(i+2 < 16 && C3Y_L[i+1] <= 0 && C3Y_L[i+2] <= 0){
            C3Y_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 2 == 16 && C3Y_L[i+1] <= 0){
            C3Y_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 1 == 16){
            C3Y_cluster_length.push_back(cluster_length_count);
          }

        }
        else{
          cluster_spacing++;
          if(i != 0 && C3Y_L[i-1] <= 0){
            cluster_length_count = 0;
          }
          else if ( i != 0){
            cluster_length_count++;
          }
        }
      }  

      // count C4X clusters
      first_cluster = true;
      cluster_spacing = 0;
      cluster_length_count = 0;

      for(int i = 0; i<72; i++){
        if(C4X_L[i] > 0 && first_cluster){
          cluster_spacing = MWPC_cluster_separation + 1;
          first_cluster = false;
        }

        if(C4X_L[i] > 0){
          if(cluster_spacing > MWPC_cluster_separation){
            C4X_clusters++;
            C4X_cluster_index.push_back(i);
          }
          cluster_length_count++;
          cluster_spacing = 0;

          if(i+2 < 72 && C4X_L[i+1] <= 0 && C4X_L[i+2] <= 0){
            C4X_cluster_length.push_back(cluster_length_count);
          }
          else if(i == 70 && C4X_L[71] <= 0){
            C4X_cluster_length.push_back(cluster_length_count);
          }
          else if(i == 71){
            C4X_cluster_length.push_back(cluster_length_count);
          }

        }
        else{
          cluster_spacing++;
          if(i != 0 && C4X_L[i-1] <= 0){
            cluster_length_count = 0;
          }
          else if ( i != 0){
            cluster_length_count++;
          }
        }
      }  


      // count C4Y clusters
      first_cluster = true;
      cluster_spacing = 0;
      cluster_length_count = 0;

      for(int i = 0; i<16; i++){
        if(C4Y_L[i] > 0 && first_cluster){
          cluster_spacing = MWPC_cluster_separation + 1;
          first_cluster = false;
        }

        if(C4Y_L[i] > 0){
          if(cluster_spacing > MWPC_cluster_separation){
            C4Y_clusters++;
            C4Y_cluster_index.push_back(i);
          }
          cluster_length_count++;
          cluster_spacing = 0;

          if(i+2 < 16 && C4Y_L[i+1] <= 0 && C4Y_L[i+2] <= 0){
            C4Y_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 2 == 16 && C4Y_L[i+1] <= 0){
            C4Y_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 1 == 16){
            C4Y_cluster_length.push_back(cluster_length_count);
          }

        }
        else{
          cluster_spacing++;
          if(i != 0 && C4Y_L[i-1] <= 0){
            cluster_length_count = 0;
          }
          else if ( i != 0){
            cluster_length_count++;
          }
        }
      }   
    }


    else{  // RIGHT
      // count C2X clusters
      first_cluster = true;
      cluster_spacing = 0;
      cluster_length_count = 0;

      for(int i = 0; i<56; i++){
        if(C2X_R[i] > 0 && first_cluster){
          cluster_spacing = MWPC_cluster_separation + 1;
          first_cluster = false;
        }

        if(C2X_R[i] > 0){
          if(cluster_spacing > MWPC_cluster_separation){
            C2X_clusters++;
            C2X_cluster_index.push_back(i);
          }
          cluster_length_count++;
          cluster_spacing = 0;

          if(i+2 < 56 && C2X_R[i+1] <= 0 && C2X_R[i+2] <= 0){
            C2X_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 2 == 56 && C2X_R[i+1] <= 0){
            C2X_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 1 == 56){
            C2X_cluster_length.push_back(cluster_length_count);
          }

        }
        else{
          cluster_spacing++;
          if(i != 0 && C2X_R[i-1] <= 0){
            cluster_length_count = 0;
          }
          else if ( i != 0){
            cluster_length_count++;
          }
        }
      }  


    
      // count C2Y clusters
      first_cluster = true;
      cluster_spacing = 0;
      cluster_length_count = 0;

      for(int i = 0; i<16; i++){
        if(C2Y_R[i] > 0 && first_cluster){
          cluster_spacing = MWPC_cluster_separation + 1;
          first_cluster = false;
        }

        if(C2Y_R[i] > 0){
          if(cluster_spacing > MWPC_cluster_separation){
            C2Y_clusters++;
            C2Y_cluster_index.push_back(i);
          }
          cluster_length_count++;
          cluster_spacing = 0;

          if(i+2 < 16 && C2Y_R[i+1] <= 0 && C2Y_R[i+2] <= 0){
            C2Y_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 2 == 16 && C2Y_R[i+1] <= 0){
            C2Y_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 1 == 16){
            C2Y_cluster_length.push_back(cluster_length_count);
          }

        }
        else{
          cluster_spacing++;
          if(i != 0 && C2Y_R[i-1] <= 0){
            cluster_length_count = 0;
          }
          else if ( i != 0){
            cluster_length_count++;
          }
        }
      }      

      // count C3X clusters
      first_cluster = true;
      cluster_spacing = 0;
      cluster_length_count = 0;

      for(int i = 0; i<64; i++){
        if(C3X_R[i] > 0 && first_cluster){
          cluster_spacing = MWPC_cluster_separation + 1;
          first_cluster = false;
        }

        if(C3X_R[i] > 0){
          if(cluster_spacing > MWPC_cluster_separation){
            C3X_clusters++;
            C3X_cluster_index.push_back(i);
          }
          cluster_length_count++;
          cluster_spacing = 0;

          if(i+2 < 64 && C3X_R[i+1] <= 0 && C3X_R[i+2] <= 0){
            C3X_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 2 == 64 && C3X_R[i+1] <= 0){
            C3X_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 1 == 64){
            C3X_cluster_length.push_back(cluster_length_count);
          }

        }
        else{
          cluster_spacing++;
          if(i != 0 && C3X_R[i-1] <= 0){
            cluster_length_count = 0;
          }
          else if ( i != 0){
            cluster_length_count++;
          }
        }
      }  

      // count C3Y clusters
      first_cluster = true;
      cluster_spacing = 0;
      cluster_length_count = 0;

      for(int i = 0; i<16; i++){
        if(C3Y_R[i] > 0 && first_cluster){
          cluster_spacing = MWPC_cluster_separation + 1;
          first_cluster = false;
        }

        if(C3Y_R[i] > 0){
          if(cluster_spacing > MWPC_cluster_separation){
            C3Y_clusters++;
            C3Y_cluster_index.push_back(i);
          }
          cluster_length_count++;
          cluster_spacing = 0;

          if(i+2 < 16 && C3Y_R[i+1] <= 0 && C3Y_R[i+2] <= 0){
            C3Y_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 2 == 16 && C3Y_R[i+1] <= 0){
            C3Y_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 1 == 16){
            C3Y_cluster_length.push_back(cluster_length_count);
          }

        }
        else{
          cluster_spacing++;
          if(i != 0 && C3Y_R[i-1] <= 0){
            cluster_length_count = 0;
          }
          else if ( i != 0){
            cluster_length_count++;
          }
        }
      }  

      // count C4X clusters
      first_cluster = true;
      cluster_spacing = 0;
      cluster_length_count = 0;

      for(int i = 0; i<72; i++){
        if(C4X_R[i] > 0 && first_cluster){
          cluster_spacing = MWPC_cluster_separation + 1;
          first_cluster = false;
        }

        if(C4X_R[i] > 0){
          if(cluster_spacing > MWPC_cluster_separation){
            C4X_clusters++;
            C4X_cluster_index.push_back(i);
          }
          cluster_length_count++;
          cluster_spacing = 0;

          if(i+2 < 72 && C4X_R[i+1] <= 0 && C4X_R[i+2] <= 0){
            C4X_cluster_length.push_back(cluster_length_count);
          }
          else if(i == 70 && C4X_R[71] <= 0){
            C4X_cluster_length.push_back(cluster_length_count);
          }
          else if(i == 71){
            C4X_cluster_length.push_back(cluster_length_count);
          }

        }
        else{
          cluster_spacing++;
          if(i != 0 && C4X_R[i-1] <= 0){
            cluster_length_count = 0;
          }
          else if ( i != 0){
            cluster_length_count++;
          }
        }
      }  


      // count C4Y clusters
      first_cluster = true;
      cluster_spacing = 0;
      cluster_length_count = 0;

      for(int i = 0; i<16; i++){
        if(C4Y_R[i] > 0 && first_cluster){
          cluster_spacing = MWPC_cluster_separation + 1;
          first_cluster = false;
        }

        if(C4Y_R[i] > 0){
          if(cluster_spacing > MWPC_cluster_separation){
            C4Y_clusters++;
            C4Y_cluster_index.push_back(i);
          }
          cluster_length_count++;
          cluster_spacing = 0;

          if(i+2 < 16 && C4Y_R[i+1] <= 0 && C4Y_R[i+2] <= 0){
            C4Y_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 2 == 16 && C4Y_R[i+1] <= 0){
            C4Y_cluster_length.push_back(cluster_length_count);
          }
          else if(i + 1 == 16){
            C4Y_cluster_length.push_back(cluster_length_count);
          }

        }
        else{
          cluster_spacing++;
          if(i != 0 && C4Y_R[i-1] <= 0){
            cluster_length_count = 0;
          }
          else if ( i != 0){
            cluster_length_count++;
          }
        }
      }    
    }    


        C2_centroid_num = 0;
        C2_centroid_den = 0;

        if(C2X_clusters == 1){
          if(selected_TOF2>=7){
            for(int i=0; i<56; i++){
              if(C2X_L[i] > 0){
                C2_centroid_num += C2X_L[i]*C2_ZLoc[i];
                C2_centroid_den += C2X_L[i];
              }
            }
          }
          else{
            for(int i=0; i<56; i++){
              if(C2X_R[i] > 0){
                C2_centroid_num += C2X_R[i]*C2_ZLoc[i];
                C2_centroid_den += C2X_R[i];
              }
            }
          }

          //C2X_centroid = double(C2_centroid_num)/C2_centroid_den;
        }


    TDC_Ck_counter = 0;
    TDC_Cpi_counter = 0;

      for(int ic = 0; ic < 14; ic++){

       
        if(((tdc_Ck[ic][0] >= TDC_Ck_min) && (tdc_Ck[ic][0] <= TDC_Ck_max))||
           ((tdc_Ck[ic][1] >= TDC_Ck_min) && (tdc_Ck[ic][1] <= TDC_Ck_max))||
           ((tdc_Ck[ic][2] >= TDC_Ck_min) && (tdc_Ck[ic][2] <= TDC_Ck_max))||
           ((tdc_Ck[ic][3] >= TDC_Ck_min) && (tdc_Ck[ic][3] <= TDC_Ck_max))||
           ((tdc_Ck[ic][4] >= TDC_Ck_min) && (tdc_Ck[ic][4] <= TDC_Ck_max))||
           ((tdc_Ck[ic][5] >= TDC_Ck_min) && (tdc_Ck[ic][5] <= TDC_Ck_max)))
          TDC_Ck_counter++;
     
        if(((tdc_Cpi[ic][0] >= TDC_Cpi_min) && (tdc_Cpi[ic][0] <= TDC_Cpi_max))||
           ((tdc_Cpi[ic][1] >= TDC_Cpi_min) && (tdc_Cpi[ic][1] <= TDC_Cpi_max))||
           ((tdc_Cpi[ic][2] >= TDC_Cpi_min) && (tdc_Cpi[ic][2] <= TDC_Cpi_max))||
           ((tdc_Cpi[ic][3] >= TDC_Cpi_min) && (tdc_Cpi[ic][3] <= TDC_Cpi_max))||
           ((tdc_Cpi[ic][4] >= TDC_Cpi_min) && (tdc_Cpi[ic][4] <= TDC_Cpi_max))||
           ((tdc_Cpi[ic][5] >= TDC_Cpi_min) && (tdc_Cpi[ic][5] <= TDC_Cpi_max)))
          TDC_Cpi_counter++; 
      }
    

    par0_ADC = 0;  par0_TDC = 0.;
    par1_ADC = 0;  par1_TDC = 0.;
    //phi_TDC = 0;

  
       
    // Compute energy weighted centroids


    for(vector<int>::iterator it = good_k_stop_bars.begin(); it != good_k_stop_bars.end(); it++){
      X_weights += ADC_Low_TARGET_ped[*it]*Xloc[*it];
      Y_weights += ADC_Low_TARGET_ped[*it]*Yloc[*it];
      total_energy += ADC_Low_TARGET_ped[*it];
    }  
  
    if(total_energy != 0){
      X_BAR = X_weights/total_energy;
      Y_BAR = Y_weights/total_energy;    
    }
    else{
      X_BAR = 999.999;
      Y_BAR = 999.999;
    }
  
    //Determine closest bar to the centroid, and draw in red
    closest_to_centroid = 1000;
    //closest_to_centroid_index = -1;
  
    for(int i = 0; i<256; i++){
      if(distance(X_BAR,Y_BAR,Xloc[i],Yloc[i]) <= closest_to_centroid){
        //closest_to_centroid_index = i;
        closest_to_centroid = distance(X_BAR,Y_BAR,Xloc[i],Yloc[i]);
      }
    }
  
    // 3rd window will show the condition for a kstop event
    for(Int_t icol=0; icol<256; icol++){
      if (has_TDC_hit_Kstop[icol]) {
        if(ADC_Low_TARGET_ped[icol]>=LG_KAON && ADC_High_TARGET_ped[icol] >=HG_KAON){ 
          //marker_ADCL_TARGET[icol]->SetMarkerColor(12); marker_ADCL_TARGET[icol]->Draw();
          h_K_Stop_Bars->Fill(Xloc[icol],Yloc[icol]);
        }     
      } 
    }   
  
    char ch_centroid[200];    

    if(X_BAR != -10000){
      sprintf(ch_centroid,"#bar{x} = %3.2f, #bar{y} = %3.2f",X_BAR, Y_BAR);
    }
    else{
      sprintf(ch_centroid,"No Centroid");
    }
  
    // Fill K Stop bars
    for(int i = 0; i<256; i++){
      if(k_stop_bar[i])
        h_TDC_Gap_Fibers_kaon->Fill(Xloc[i],Yloc[i]);
    }       
 
    TF1 *fit_line_ADC4;
    TF1 *fit_line_ADCA;
  
  
        // TARGET ROTATION (90 deg.)
        
    // ROTATE_CHANGE: any line tagged with ROTATE_CHANGE is used when tof1 is 6 or 12.
        
    gap_to_fit_rotate = gap_to_fit;
         
    if((gap_to_fit==12 || gap_to_fit==6 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){    
      for(int i=0; i<256; i++){
        ADC_High_TARGET_temp[i] = 9999;
        ADC_Low_TARGET_temp[i] = 9999;
        //adc_low_target_temp[i] = 9999;
        has_TDC_hit_temp[i] = false;
        has_TDC_hit_Kstop_temp[i] = false;
        HG_TARGET_ADC_Thr_temp[i] = 9999;
        LG_TARGET_ADC_Thr_temp[i] = 9999;
        TARGET_High_has_neighbours_temp[i] = 9999;
        k_stop_bar_temp[i] = false;
      }
  
      if(gap_to_fit==11)
        gap_to_fit = 2;
      else if(gap_to_fit==12)
        gap_to_fit = 3;
      else if(gap_to_fit==1)
        gap_to_fit = 4;
      else if(gap_to_fit==5)
        gap_to_fit = 8;
      else if(gap_to_fit==6)
        gap_to_fit = 9;
      else if(gap_to_fit==7)
        gap_to_fit = 10;
    

      h_TDC_Gap_Fibers->Reset();
      h_target_TDC->Reset();
      h_target_TDC_copy->Reset();
      h_TDC_selected->Reset();
      h_TDC_selected2->Reset();
      h_TARGET_LG_Blue->Reset();
    
      for(int i = 0; i<256; i++){
        ADC_High_TARGET_temp[i] = ADC_High_TARGET[i];
        ADC_Low_TARGET_temp[i] = ADC_Low_TARGET[i];
        HG_TARGET_ADC_Thr_temp[i] = HG_TARGET_ADC_Thr[i];
        LG_TARGET_ADC_Thr_temp[i] = LG_TARGET_ADC_Thr[i];
        has_TDC_hit_temp[i] = has_TDC_hit[i];
        has_TDC_hit_Kstop_temp[i] = has_TDC_hit_Kstop[i];
        //adc_low_target_temp[i] = adc_low_target[i];
        TARGET_High_has_neighbours_temp[i] = TARGET_High_has_neighbours[i];
        k_stop_bar_temp[i] = k_stop_bar[i];
      }
    
      for(int i=0; i<256; i++){
      ADC_High_TARGET[i] = ADC_High_TARGET_temp[TARGET_Rotated_index[i]];
      ADC_Low_TARGET[i] = ADC_Low_TARGET_temp[TARGET_Rotated_index[i]];
      HG_TARGET_ADC_Thr[i] = HG_TARGET_ADC_Thr_temp[TARGET_Rotated_index[i]];
      LG_TARGET_ADC_Thr[i] = LG_TARGET_ADC_Thr_temp[TARGET_Rotated_index[i]];
      has_TDC_hit[i] = has_TDC_hit_temp[TARGET_Rotated_index[i]];
      has_TDC_hit_Kstop[i] = has_TDC_hit_Kstop_temp[TARGET_Rotated_index[i]];
      adc_low_target[i] = adc_low_target[TARGET_Rotated_index[i]];
      TARGET_High_has_neighbours[i] = TARGET_High_has_neighbours_temp[TARGET_Rotated_index[i]];
      k_stop_bar[i] = k_stop_bar_temp[TARGET_Rotated_index[i]];
      }

      max_index = 0;
      max_index2 = 0;
      max_index3 = 0;
      max_index4 = 0;

      max_ADC = -100000000;
      max_ADC2 = -100000000;
      max_ADC3 = -100000000;
      max_ADC4 = -100000000;
    
      for(Int_t q=0; q<256; q++){
        if(ADC_High_TARGET[q]>max_ADC) {
          max_index = q;
          max_ADC = ADC_Low_TARGET[q];
        }
      }

      for(Int_t q=0; q<256; q++){
        if (q == max_index) continue;   
        else {
          if(ADC_High_TARGET[q]>max_ADC2) {
            max_index2 = q;
            max_ADC2 = ADC_Low_TARGET[q];
          }
        }
      }

      for(Int_t q=0; q<256; q++){
        if ((q == max_index) || (q == max_index2)) continue;    
        else {
          if(ADC_High_TARGET[q]>max_ADC3) {
            max_index3 = q;
            max_ADC3 = ADC_Low_TARGET[q];
          }
        }
      }

      for(Int_t q=0; q<256; q++){
        if ((q == max_index) || (q == max_index2) || (q == max_index3)) continue;   
        else {
          if(ADC_High_TARGET[q]>max_ADC4) {
            max_index4 = q;
            max_ADC4 = ADC_Low_TARGET[q];
          }
        }
      }
  
      x_cent = Xloc[max_index];
      y_cent = Yloc[max_index];
  
  
      for(Int_t j=0; j<256; j++){
        hyp[j] = sqrt(pow(x_cent - Xloc[j],2) + pow(y_cent - Yloc[j],2));
      }
    
      count = 0;
      hit_count = 0;
    
      for(Int_t i=0; i<256; i++){
        if(TARGET_High_has_neighbours[i]){ 
          if(ADC_High_TARGET[i]>=Angle_ADC_cut && has_TDC_hit[i]){
            h_target_TDC->Fill(Xloc[i],Yloc[i]);
            h_target_TDC_copy->Fill(Xloc[i],Yloc[i]);
          
            if(!k_stop_bar[i])
              h_TDC_Gap_Fibers->Fill(Xloc[i],Yloc[i]);
            
            count++;
          
          
            if(IsIn(i,channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],   // Additional weight for fibers close to the edge if hit
                      channel[gap_to_fit-1][2], channel[gap_to_fit-1][3], 
                      channel[gap_to_fit-1][4], channel[gap_to_fit-1][5], 
                      channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){
            
              if(!k_stop_bar[i]){     
                h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
                h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
              }
        
              h_TDC_selected2->Fill(Xloc[i], Yloc[i]);
              h_TDC_selected2->Fill(Xloc[i], Yloc[i]);
  
              h_TDC_selected->Fill(Xloc[i], Yloc[i]);
              h_TDC_selected->Fill(Xloc[i], Yloc[i]);
            }
          } 

          if(ADC_High_TARGET[i]>=Angle_ADC_cut){
            h_target_ADC->Fill(Xloc[i],Yloc[i]);
            x_inc += Xloc[i];
            y_inc += Yloc[i];
            hit_count++;
          }
    
          if(ADC_High_TARGET[i]<0 && ADC_Low_TARGET[i]>=0){
            if(!k_stop_bar[i])  
              h_TDC_Gap_Fibers->Fill(Xloc[i],Yloc[i]);
              h_TARGET_LG_Blue->Fill(Xloc[i],Yloc[i]);
              count++;
          }
        }
      }
    }


    /// Data counters
    has_data_TDC2 = 0;
    has_data_ADC2 = 0;
    has_data_ADC3 = 0;
    has_data_ADC4 = 0;
    has_data_ADCA = 0;

    h_target_ADC4->Fill(Xloc[max_index],Yloc[max_index]);
    h_target_ADC4->Fill(Xloc[max_index2],Yloc[max_index2]);
    h_target_ADC4->Fill(Xloc[max_index3],Yloc[max_index3]);
    h_target_ADC4->Fill(Xloc[max_index4],Yloc[max_index4]);

    ///// Select edge fiber for track fitting


    Xloc_gap = 0;
    Yloc_gap = 0;

    if (gap_to_fit == 0) {
      for (int j=0; j<8; j++) {
        if ((ADC_High_TARGET[channel[0][j]]>0) && (has_TDC_hit[channel[0][j]])) {
          Xloc_gap = Xloc[channel[0][j]];
          Yloc_gap = Yloc[channel[0][j]];
        }
      }
      if (Xloc_gap == 0) Xloc_gap = Xloc[channel[0][3]];
      if (Yloc_gap == 0) Yloc_gap = Yloc[channel[0][3]];
    }
    if (gap_to_fit == 1) {
      for (int j=0; j<8; j++) {
        if ((ADC_High_TARGET[channel[1][j]]>0) && (has_TDC_hit[channel[1][j]])) {
          Xloc_gap = Xloc[channel[1][j]];
          Yloc_gap = Yloc[channel[1][j]];
        }
      }
      if (Xloc_gap == 0) Xloc_gap = Xloc[channel[1][3]];
      if (Yloc_gap == 0) Yloc_gap = Yloc[channel[1][3]];
    }
    if (gap_to_fit == 2) {
      for (int j=0; j<8; j++) {
        if ((ADC_High_TARGET[channel[2][j]]>0) && (has_TDC_hit[channel[2][j]])) {
          Xloc_gap = Xloc[channel[2][j]];
          Yloc_gap = Yloc[channel[2][j]];
        }
      }
      if (Xloc_gap == 0) Xloc_gap = Xloc[channel[2][3]];
      if (Yloc_gap == 0) Yloc_gap = Yloc[channel[2][3]];
    }
    if (gap_to_fit == 3) {
      for (int j=0; j<8; j++) {
        if ((ADC_High_TARGET[channel[3][j]]>0) && (has_TDC_hit[channel[3][j]])) {
          Xloc_gap = Xloc[channel[3][j]];
          Yloc_gap = Yloc[channel[3][j]];
        } 
      }
      if (Xloc_gap == 0) Xloc_gap = Xloc[channel[3][3]];
      if (Yloc_gap == 0) Yloc_gap = Yloc[channel[3][3]];
    }
    if (gap_to_fit == 4) {
      for (int j=0; j<8; j++) {
        if ((ADC_High_TARGET[channel[4][j]]>0) && (has_TDC_hit[channel[4][j]])) {
          Xloc_gap = Xloc[channel[4][j]];
          Yloc_gap = Yloc[channel[4][j]];
        }
      }
      if (Xloc_gap == 0) Xloc_gap = Xloc[channel[4][3]];
      if (Yloc_gap == 0) Yloc_gap = Yloc[channel[4][3]];
    }
    if (gap_to_fit == 5) {
      for (int j=0; j<8; j++) {
        if ((ADC_High_TARGET[channel[5][j]]>0) && (has_TDC_hit[channel[5][j]])) {
          Xloc_gap = Xloc[channel[5][j]];
          Yloc_gap = Yloc[channel[5][j]];
        }
      }
      if (Xloc_gap == 0) Xloc_gap = Xloc[channel[5][3]];
      if (Yloc_gap == 0) Yloc_gap = Yloc[channel[5][3]];
    }
    if (gap_to_fit == 6) {
      for (int j=0; j<8; j++) {
        if ((ADC_High_TARGET[channel[6][j]]>0) && (has_TDC_hit[channel[6][j]])) {
          Xloc_gap = Xloc[channel[6][j]];
          Yloc_gap = Yloc[channel[6][j]];
        }
      }
      if (Xloc_gap == 0) Xloc_gap = Xloc[channel[6][3]];
      if (Yloc_gap == 0) Yloc_gap = Yloc[channel[6][3]];
    }
    if (gap_to_fit == 7) {
      for (int j=0; j<8; j++) {
        if ((ADC_High_TARGET[channel[7][j]]>0) && (has_TDC_hit[channel[7][j]])) {
          Xloc_gap = Xloc[channel[7][j]];
          Yloc_gap = Yloc[channel[7][j]];
        }
      }
      if (Xloc_gap == 0) Xloc_gap = Xloc[channel[7][3]];
      if (Yloc_gap == 0) Yloc_gap = Yloc[channel[7][3]];
    }
    if (gap_to_fit == 8) {
      for (int j=0; j<8; j++) {
        if ((ADC_High_TARGET[channel[8][j]]>0) && (has_TDC_hit[channel[8][j]])) {
          Xloc_gap = Xloc[channel[8][j]];
          Yloc_gap = Yloc[channel[8][j]];
        }
      }
      if (Xloc_gap == 0) Xloc_gap = Xloc[channel[8][3]];
      if (Yloc_gap == 0) Yloc_gap = Yloc[channel[8][3]];
    }
    if (gap_to_fit == 9) {
      for (int j=0; j<8; j++) {
        if ((ADC_High_TARGET[channel[9][j]]>0) && (has_TDC_hit[channel[9][j]])) {
          Xloc_gap = Xloc[channel[9][j]];
          Yloc_gap = Yloc[channel[9][j]];
        }
      }
      if (Xloc_gap == 0) Xloc_gap = Xloc[channel[9][3]];
      if (Yloc_gap == 0) Yloc_gap = Yloc[channel[9][3]];
    }
    if (gap_to_fit == 10) {
      for (int j=0; j<8; j++) {
        if ((ADC_High_TARGET[channel[10][j]]>0) && (has_TDC_hit[channel[10][j]])) {
          Xloc_gap = Xloc[channel[10][j]];
          Yloc_gap = Yloc[channel[10][j]];
        }
      }
      if (Xloc_gap == 0) Xloc_gap = Xloc[channel[10][3]];
      if (Yloc_gap == 0) Yloc_gap = Yloc[channel[10][3]];
    }
    if (gap_to_fit == 11) {
      for (int j=0; j<8; j++) {
        if ((ADC_High_TARGET[channel[11][j]]>0) && (has_TDC_hit[channel[11][j]])) {
          Xloc_gap = Xloc[channel[11][j]];
          Yloc_gap = Yloc[channel[11][j]];
        }
      }
      if (Xloc_gap == 0) Xloc_gap = Xloc[channel[11][3]];
      if (Yloc_gap == 0) Yloc_gap = Yloc[channel[11][3]];
    }

    /// Determine which gap point is closest to edge fiber

    xdistance1 = pow((TOF_Xloc[(gap_to_fit*3)]-Xloc_gap),2);
    ydistance1 = pow((TOF_Yloc[(gap_to_fit*3)]-Yloc_gap),2);

    xdistance2 = pow((TOF_Xloc[(gap_to_fit*3)+1]-Xloc_gap),2);
    ydistance2 = pow((TOF_Yloc[(gap_to_fit*3)+1]-Yloc_gap),2);

    xdistance3 = pow((TOF_Xloc[(gap_to_fit*3)+2]-Xloc_gap),2);
    ydistance3 = pow((TOF_Yloc[(gap_to_fit*3)+2]-Yloc_gap),2);

    xhyp1 = double(sqrt(double(xdistance1) + double(ydistance1)));
    xhyp2 = double(sqrt(double(xdistance2) + double(ydistance2)));
    xhyp3 = double(sqrt(double(xdistance3) + double(ydistance3)));


    if ((xhyp1 <= xhyp2) && (xhyp1 <= xhyp3)) {
      closest_gap_point = 0;
    }
    else if (xhyp2 <= xhyp3) {
      closest_gap_point = 1;
    }
    else {
      closest_gap_point = 2;
    }

    h_target_ADC3->Fill(TOF_Xloc[(gap_to_fit*3) + closest_gap_point],TOF_Yloc[(gap_to_fit*3) + closest_gap_point]);

    ////Fill TARGET tracking histograms with data

    for(Int_t i=0; i<256; i++){

      if ((i == max_index) || (i == max_index2) || (i == max_index3) || (i == max_index4)) continue;

      if((has_TDC_hit[i]) && (ADC_High_TARGET[i]>=Angle_ADC_cut) && (hyp[i] > 12.5)) {
        h_target_ADC2->Fill(Xloc[i],Yloc[i]);
        has_data_ADC2++;
      }
      if((has_TDC_hit[i]) && (ADC_High_TARGET[i]>=Angle_ADC_cut) && (adc_low_target[i]-ADC_cut_TARGET2>0)) {
        h_target_TDC2->Fill(Xloc[i],Yloc[i]);
        has_data_TDC2++;
      }
      if((has_TDC_hit[i]) && (ADC_High_TARGET[i]>=Angle_ADC_cut) && (adc_low_target[i]-ADC_cut_TARGET2>0)) {
        h_target_ADC3->Fill(Xloc[i],Yloc[i]);
        has_data_ADC3++;
      }
      if((has_TDC_hit[i]) && (ADC_High_TARGET[i]>=Angle_ADC_cut) && (hyp[i] < 12.5)) {
        h_target_ADC4->Fill(Xloc[i],Yloc[i]);
        has_data_ADC4++;
      } 
    }

    /*
    if (display_fits ==1){
      h_target_ADC2->SetMarkerStyle(25);
      h_target_ADC2->SetMarkerColor(4);

      h_target_TDC2->SetMarkerStyle(25);
      h_target_TDC2->SetMarkerColor(4);

      //cout << "  " << endl;

      if (has_data_TDC2 > 1){
        phi_TDC=(180/PI)*atan(par1_TDC);
        sprintf(angle_string_TDC,"Phi = %3.4f deg.",phi_TDC);
        tex_angle_TDC=new TLatex(0,-9,angle_string_TDC);
        tex_angle_TDC->SetTextSize(0.05);
      }
      else {
        cout << "Histo Fit 3 Is Empty" << endl;
      }

      h_target_ADC4->SetMarkerStyle(25);
      h_target_ADC4->SetMarkerColor(4);
    }  /// TO CHECK !!
    */
    
    //// Used for primary track fitting

    h_target_ADC3->SetMarkerStyle(25);
    h_target_ADC3->SetMarkerColor(4);

    if (has_data_ADC3 > 1){
    

      //xcoord = 0;
      //unique_x = 0;


      if (has_data_ADCA > 2) {


        h_target_ADCA->Fill(TOF_Xloc[(gap_to_fit*3) + closest_gap_point],TOF_Yloc[(gap_to_fit*3) + closest_gap_point]);
      
        /// Refit histogram with filtered data + closest gap point
        h_target_ADCA->Fit("pol1", "QCM");
        fit_line_ADCA = h_target_ADCA->GetFunction("pol1");

      
        /// Retrieve slope & intercept for refit line
        par0_ADC=fit_line_ADCA->GetParameter(0);
        par1_ADC=fit_line_ADCA->GetParameter(1);



        /// Ring intercept coordinates
        //determinant = 4*(pow(par0_ADC,2))*(pow(par1_ADC,2)) - 4*(pow(par1_ADC,2) + 1)*(pow(par0_ADC,2)-1600);
        //x_circle_int1 = (-2*(par0_ADC)*(par1_ADC) + sqrt(determinant))/(2*(pow(par1_ADC,2) + 1));
        //y_circle_int1 = (par1_ADC)*(x_circle_int1) + par0_ADC;

        //x_circle_int2 = (-2*(par0_ADC)*(par1_ADC) - sqrt(determinant))/(2*(pow(par1_ADC,2) + 1));
        //y_circle_int2 = (par1_ADC)*(x_circle_int2) + par0_ADC;

        //if (unique_x == 1) {
          //x_circle_int1 = xcoord;
          //y_circle_int1 = sqrt(1600-(pow(xcoord,2)));

          //x_circle_int2 = xcoord;
          //y_circle_int2 = sqrt(1600-(pow(xcoord,2)))*-1;
        //}

        //SFTxdistance1 = pow((x_circle_int1-Xloc_gap),2);
        //SFTydistance1 = pow((y_circle_int1-Yloc_gap),2);

        //SFTxdistance2 = pow((x_circle_int2-Xloc_gap),2);
        //SFTydistance2 = pow((y_circle_int2-Yloc_gap),2);

        //SFTxhyp1 = double(sqrt(double(SFTxdistance1) + double(SFTydistance1)));
        //SFTxhyp2 = double(sqrt(double(SFTxdistance2) + double(SFTydistance2)));

        //SFT_x_intercept = 999.;
        //SFT_y_intercept = 999.;

        //if (SFTxhyp1 < SFTxhyp2) {
          //SFT_x_intercept = x_circle_int1;
          //SFT_y_intercept = y_circle_int1;SFTxdistance1
        //}
        //else {
          //SFT_x_intercept = x_circle_int2;
          //SFT_y_intercept = y_circle_int2;
        //}


        /// Determine angle phi between centre of TARGET and SFT ring intercept
        //SFT_phi = 999.;

        //if (SFT_x_intercept >0) {
          //SFT_phi = (atan2(SFT_x_intercept,SFT_y_intercept))*57.2957795;
        //}
        //else {
          //SFT_phi = 360 + (atan2(SFT_x_intercept,SFT_y_intercept))*57.2957795;
        //}

        if (has_data_ADC4 > 1) {

          h_target_ADC4->Fit("pol1", "Q0CM");
          fit_line_ADC4 = h_target_ADC4->GetFunction("pol1");
          fit_line_ADC4->SetLineWidth(2);
          fit_line_ADC4->SetLineColor(1);

          x_intercept = double(double(par0_ADC - par0_TDC)/double(par1_TDC - par1_ADC));
          y_intercept = double(par1_ADC)*double(x_intercept) + par0_ADC;

          for(int i=0; i<256; i++){
            x_distance4[i] = 999.;
            y_distance4[i] = 999.;
            distances4[i] = 999.;
          }

          for (int q=0; q<256; q++) {
            x_distance4[q] = pow((Xloc[q]-x_intercept),2);
            y_distance4[q] = pow((Yloc[q]-y_intercept),2);
  
            distances4[q] = double(sqrt(double(x_distance4[q]) + double(y_distance4[q])));
          }

          min_distance = 9999.;

          for (int q=0; q<256; q++) {
            x_distance4[q] = pow((Xloc[q]-x_intercept),2);
            y_distance4[q] = pow((Yloc[q]-y_intercept),2);
  
            distances4[q] = double(sqrt(double(x_distance4[q]) + double(y_distance4[q])));
          }

          for (int q=0; q<256; q++) {
            if (distances4[q] < min_distance) {
              min_distance = distances4[q];
            }
          }

        
        }
      }
    }


   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    TF1 *fit_line_GoodLG;           TF1 *fit_line_GoodLG_weighed;
    TF1 *fit_line_TDC_selected;     TF1 *fit_line_TDC_selected2;
    TF1 *fit_line_TDC_Gap_Fibers;



    for(int i=0; i<12; i++){
      for(int j=0; j<3; j++){
        for(int k=0; k<2; k++){
          Gap[i][j][k] = 0;
        }
      }
    }
  
    a_fit_TDC_selected = 0.;               b_fit_TDC_selected = 0.;
    a_fit_GoodLG = 0.;                     b_fit_GoodLG = 0.;
    a_fit_GoodLG_weighted = 0.;            b_fit_GoodLG_weighted = 0.;
    a_fit_TDC_selected_weighted = 0.;      b_fit_TDC_selected_weighted = 0.;
    a_fit_TDC_Gap_Fibers = 0.;             b_fit_TDC_Gap_Fibers = 0.;

    //cout << "DEBUG4 :  " << a_fit_TDC_Gap_Fibers << "   " << b_fit_TDC_Gap_Fibers << endl;

    for(int g=0; g<12; g++)
    {
      Gap[g][0][0] = TOF_Xloc[3*g];
      Gap[g][1][0] = TOF_Xloc[3*g+1];
      Gap[g][2][0] = TOF_Xloc[3*g+2];
    
      Gap[g][0][1] = TOF_Yloc[3*g];
      Gap[g][1][1] = TOF_Yloc[3*g+1];
      Gap[g][2][1] = TOF_Yloc[3*g+2];
    }

   
    for(int l=0; l<256; l++)
    {
      if(adc_low_target[l]>=LG_TARGET_ADC_Thr[l] && has_TDC_hit[l])
      {
        h_GoodLG->Fill(Xloc[l], Yloc[l]);
        h_GoodLG_copy->Fill(Xloc[l], Yloc[l]);
        h_GoodLG_weighted->Fill(Xloc[l], Yloc[l]);
      }
    }

    for(int i=0; i<12; i++)
    {
      h_Circle->Fill(Gap[i][0][0], Gap[i][0][1]);
      h_Circle->Fill(Gap[i][1][0], Gap[i][1][1]);
      h_Circle->Fill(Gap[i][2][0], Gap[i][2][1]);
    }

    h_Target_Center->Fill(0., 0.);

    for(int t=0; t<256; t++)
    {
      h_Target->Fill(Xloc[t],Yloc[t]);
    }

    h_max->Fill(Xloc[max_index], Yloc[max_index]);        h_max_copy->Fill(Xloc[max_index], Yloc[max_index]);
 
    h_kaon->Fill(Xloc[max_index], Yloc[max_index]);       h_kaon_copy->Fill(Xloc[max_index], Yloc[max_index]);
    h_kaon->Fill(Xloc[max_index2], Yloc[max_index2]);     h_kaon_copy->Fill(Xloc[max_index2], Yloc[max_index2]);
    h_kaon->Fill(Xloc[max_index3], Yloc[max_index3]);     h_kaon_copy->Fill(Xloc[max_index3], Yloc[max_index3]);
    h_kaon->Fill(Xloc[max_index4], Yloc[max_index4]);     h_kaon_copy->Fill(Xloc[max_index4], Yloc[max_index4]);

    // ROTATE_CHANGE
    if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){  
      for(int g=0; g<3; g++){
        h_TOF1_rotate->Fill(Gap[gap_to_fit_rotate-1][g][0], Gap[gap_to_fit_rotate-1][g][1]);
        h_TDC_selected->Fill(Gap[gap_to_fit-1][g][0], Gap[gap_to_fit-1][g][1]);  
      } 
    }
    else{
      for(int g=0; g<3; g++){
        h_TOF1->Fill(Gap[gap_to_fit-1][g][0], Gap[gap_to_fit-1][g][1]);
        h_TDC_selected->Fill(Gap[gap_to_fit-1][g][0], Gap[gap_to_fit-1][g][1]);  
      } 
    }

    h_Circle->SetMarkerStyle(5);
    h_Circle->SetMarkerSize(1.2);
    h_Circle->SetLineWidth(2);

    h_Target->SetMarkerStyle(25);         h_Target_Center->SetMarkerStyle(5);
    h_Target->SetMarkerColor(1);          h_Target_Center->SetMarkerColor(1);
    h_Target->SetMarkerSize(0.8);         h_Target_Center->SetMarkerSize(1);

    h_GoodLG->SetMarkerStyle(21);         h_GoodLG_copy->SetMarkerStyle(21);          h_GoodLG_weighted->SetMarkerStyle(21);
    h_GoodLG->SetMarkerColor(2);          h_GoodLG_copy->SetMarkerColor(2);           h_GoodLG_weighted->SetMarkerColor(2);
    h_GoodLG->SetMarkerSize(0.8);         h_GoodLG_copy->SetMarkerSize(0.8);          h_GoodLG_weighted->SetMarkerSize(0.8);

    h_target_TDC->SetMarkerStyle(25);     h_target_TDC_copy->SetMarkerStyle(25);      h_TDC_selected->SetMarkerStyle(25);   h_TDC_selected2->SetMarkerStyle(25);
    h_target_TDC->SetMarkerColor(2);      h_target_TDC_copy->SetMarkerColor(2);       h_TDC_selected->SetMarkerColor(2);;   h_TDC_selected2->SetMarkerColor(2);
    h_target_TDC->SetMarkerSize(0.8);     h_target_TDC_copy->SetMarkerSize(0.8);      h_TDC_selected->SetMarkerSize(0.8);   h_TDC_selected2->SetMarkerSize(0.8);

    h_kaon->SetMarkerStyle(21);           h_kaon_copy->SetMarkerStyle(21);            h_TDC_Gap_Fibers->SetMarkerStyle(25);
    h_kaon->SetMarkerColor(1);            h_kaon_copy->SetMarkerColor(1);             h_TDC_Gap_Fibers->SetMarkerColor(2);
    h_kaon->SetMarkerSize(0.8);           h_kaon_copy->SetMarkerSize(0.8);            h_TDC_Gap_Fibers->SetMarkerSize(0.8);
  
    h_TDC_Gap_Fibers_kaon->SetMarkerStyle(21);
    h_TDC_Gap_Fibers_kaon->SetMarkerColor(kBlue-6);
    h_TDC_Gap_Fibers_kaon->SetMarkerSize(0.8);
  
    h_max->SetMarkerStyle(21);            h_max_copy->SetMarkerStyle(21);             h_TARGET_LG_Blue->SetMarkerStyle(25);
    h_max->SetMarkerColor(kViolet+1);     h_max_copy->SetMarkerColor(kViolet+1);      h_TARGET_LG_Blue->SetMarkerColor(4);
    h_max->SetMarkerSize(1);              h_max_copy->SetMarkerSize(1);               h_TARGET_LG_Blue->SetMarkerSize(0.8);

    h_TOF1->SetMarkerStyle(20);           h_TOF1_closest->SetMarkerStyle(20);         h_TDC_Gap_Fibers_copy->SetMarkerStyle(25); // ROTATE_CHANGE   
    h_TOF1->SetMarkerColor(2);            h_TOF1_closest->SetMarkerColor(3);          h_TDC_Gap_Fibers_copy->SetMarkerColor(2); // ROTATE_CHANGE
    h_TOF1->SetMarkerSize(1.5);           h_TOF1_closest->SetMarkerSize(1.5);         h_TDC_Gap_Fibers_copy->SetMarkerSize(0.8); // ROTATE_CHANGE

    h_TOF1_rotate->SetMarkerStyle(20);    h_TOF1_closest_rotate->SetMarkerStyle(20);  // ROTATE_CHANGE
    h_TOF1_rotate->SetMarkerColor(2);     h_TOF1_closest_rotate->SetMarkerColor(3);   // ROTATE_CHANGE
    h_TOF1_rotate->SetMarkerSize(1.5);    h_TOF1_closest_rotate->SetMarkerSize(1.5);  // ROTATE_CHANGE
  
    h_Centroid->SetMarkerStyle(5);
    h_Centroid->SetMarkerColor(1);         
    h_Centroid->SetMarkerSize(4); 
  

    ell->SetFillStyle(0);                 ell_Target->SetFillStyle(0);                 ell_L1->SetFillStyle(0);
    ell->SetLineColor(6);                 ell_Target->SetLineColor(1);                 ell_L1->SetLineColor(4);
    ell->SetLineWidth(1);                 ell_Target->SetLineWidth(1);                 ell_L1->SetLineWidth(1);                        

    Gap1l->SetLineWidth(10);
    Gap1l->SetLineColor(15);
  
    Gap2l->SetLineWidth(10);
    Gap2l->SetLineColor(15);
  
    Gap3l->SetLineWidth(10);
    Gap3l->SetLineColor(15);
  
    Gap4l->SetLineWidth(10);
    Gap4l->SetLineColor(15);
  
    Gap5l->SetLineWidth(10);
    Gap5l->SetLineColor(15);
  
    Gap6l->SetLineWidth(10);
    Gap6l->SetLineColor(15);
  
    Gap7l->SetLineWidth(10);
    Gap7l->SetLineColor(15);
  
    Gap8l->SetLineWidth(10);
    Gap8l->SetLineColor(15);
  
    Gap9l->SetLineWidth(10);
    Gap9l->SetLineColor(15);

    Gap10l->SetLineWidth(10);
    Gap10l->SetLineColor(15);

    Gap11l->SetLineWidth(10);
    Gap11l->SetLineColor(15);

    Gap12l->SetLineWidth(10);
    Gap12l->SetLineColor(15);




    //c2->cd(6);
    //h_TDC_selected->Draw();
    h_TDC_selected->Fit("pol1", "QCM");
    fit_line_TDC_selected = h_TDC_selected->GetFunction("pol1");
    fit_line_TDC_selected->SetLineWidth(2);
    fit_line_TDC_selected->SetLineColor(1);
    a_fit_TDC_selected=fit_line_TDC_selected->GetParameter(1);      
    b_fit_TDC_selected=fit_line_TDC_selected->GetParameter(0);
    //HorizontalAxis->Draw("same");

    if(h_GoodLG->GetEntries()!=0){
      h_GoodLG->Fit("pol1", "QCM");
      fit_line_GoodLG = h_GoodLG->GetFunction("pol1");
      fit_line_GoodLG->SetLineWidth(2);
      fit_line_GoodLG->SetLineColor(1);
      a_fit_GoodLG=fit_line_GoodLG->GetParameter(1);      
      b_fit_GoodLG=fit_line_GoodLG->GetParameter(0);
    }
    //else{
    //cout << "Empty Histogram" << endl;
    //}
  


    // Add weight to center point of selected TOF1 (gap_to_fit_rotate)
    for(int i = 0; i<3; i++){
      h_TDC_Gap_Fibers->Fill(Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
    }
  
    // Fit 1: Targets hit and center point of selected TOF1
      h_TDC_Gap_Fibers->Fit("pol1", "QCM");
      fit_line_TDC_Gap_Fibers = h_TDC_Gap_Fibers->GetFunction("pol1");
      fit_line_TDC_Gap_Fibers->SetLineWidth(2);
      fit_line_TDC_Gap_Fibers->SetLineColor(2);
      a_fit_TDC_Gap_Fibers = h_TDC_Gap_Fibers->GetFunction("pol1")->GetParameter(1);      
      b_fit_TDC_Gap_Fibers = h_TDC_Gap_Fibers->GetFunction("pol1")->GetParameter(0);
    
      //h_TDC_Gap_Fibers->Reset();
      //cout << "NB :  " << h_TDC_Gap_Fibers->GetEntries() << endl;
      //cout << "DEBUG5 :  " << a_fit_TDC_Gap_Fibers << "   " << b_fit_TDC_Gap_Fibers << endl;


    // ReFill histogram with points closer than max_dist
    h_TDC_Gap_Fibers->Reset();
    h_TDC_Gap_Fibers_copy->Reset();

    for(int i = 0; i<256; i++){   
      if(distance_to_line(Xloc[i],Yloc[i],a_fit_TDC_Gap_Fibers,b_fit_TDC_Gap_Fibers) <= max_dist && TARGET_High_has_neighbours[i] && !k_stop_bar[i]){   
          if(ADC_High_TARGET[i]>=Angle_ADC_cut && has_TDC_hit[i]){
          h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
          h_TDC_Gap_Fibers_copy->Fill(Xloc[TARGET_Rotated_index[i]],Yloc[TARGET_Rotated_index[i]]); // ROTATE_CHANGE
      
      
          if(IsIn(i,channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],   
                    channel[gap_to_fit-1][2], channel[gap_to_fit-1][3], 
                    channel[gap_to_fit-1][4], channel[gap_to_fit-1][5], 
                    channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){
            h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);     
            h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
            h_TDC_Gap_Fibers_copy->Fill(Xloc[TARGET_Rotated_index[i]],Yloc[TARGET_Rotated_index[i]]); // ROTATE_CHANGE
            h_TDC_Gap_Fibers_copy->Fill(Xloc[TARGET_Rotated_index[i]],Yloc[TARGET_Rotated_index[i]]); // ROTATE_CHANGE      
          }
        }
        if(ADC_High_TARGET[i]<0 && ADC_Low_TARGET[i]>=0){
          h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
          h_TDC_Gap_Fibers_copy->Fill(Xloc[TARGET_Rotated_index[i]],Yloc[TARGET_Rotated_index[i]]);  // ROTATE_CHANGE
        } 
      } 
    }
  
    if(h_TDC_Gap_Fibers->GetEntries() != 0){
      //Fit 2: Targets hit within max_distance of the best fit line from Fit 1
      h_TDC_Gap_Fibers->Fit("pol1", "QCM");
      fit_line_TDC_Gap_Fibers = h_TDC_Gap_Fibers->GetFunction("pol1");
      fit_line_TDC_Gap_Fibers->SetLineWidth(2);
      fit_line_TDC_Gap_Fibers->SetLineColor(2);
      a_fit_TDC_Gap_Fibers=fit_line_TDC_Gap_Fibers->GetParameter(1);      
      b_fit_TDC_Gap_Fibers=fit_line_TDC_Gap_Fibers->GetParameter(0);
      }
    //else
      //continue;
      //cout << "Empty FIT 2" << endl;
  


    for(int i=0; i<2; i++){
      x_int_TDC[i] = 999.;                           y_int_TDC[i] = 999.;
      x_int_TDC_selected[i] = 999.;                  y_int_TDC_selected[i] = 999.;
      x_int_TDC_selected_weighted[i] = 999.;         y_int_TDC_selected_weighted[i] = 999.;
      x_int_GoodLG[i] = 999.;                        y_int_GoodLG[i] = 999.;
      x_int_GoodLG_weighted[i] = 999.;               y_int_GoodLG_weighted[i] = 999.;
      x_int_TDC_Gap_Fibers[i] = 999.;                y_int_TDC_Gap_Fibers[i] = 999.;
      x_int_TDC_SFT[i] = 999.;                       y_int_TDC_SFT[i] = 999.;
      x_int_GoodLG_SFT[i] = 999.;                    y_int_GoodLG_SFT[i] = 999.;
      x_int_TDC_Gap_Fibers_SFT[i] = 999.;            y_int_TDC_Gap_Fibers_SFT[i] = 999.;
      x_int_TARGET[i] = 999.;                        y_int_TARGET[i] = 999.;
    }

    //cout << "DEBUG2 :  " << x_int_TDC_Gap_Fibers[0] << endl;

    x_int_TDC_selected[0] = intersectx1(a_fit_TDC_selected, b_fit_TDC_selected, R_TOF1);
    x_int_TDC_selected[1] = intersectx2(a_fit_TDC_selected, b_fit_TDC_selected, R_TOF1);
    y_int_TDC_selected[0] = y1_int(x_int_TDC_selected[0], a_fit_TDC_selected, b_fit_TDC_selected);
    y_int_TDC_selected[1] = y2_int(x_int_TDC_selected[1], a_fit_TDC_selected, b_fit_TDC_selected);


    x_int_GoodLG[0] = intersectx1(a_fit_GoodLG, b_fit_GoodLG, R_TOF1);
    x_int_GoodLG[1] = intersectx2(a_fit_GoodLG, b_fit_GoodLG, R_TOF1);
    y_int_GoodLG[0] = y1_int(x_int_GoodLG[0], a_fit_GoodLG, b_fit_GoodLG);
    y_int_GoodLG[1] = y2_int(x_int_GoodLG[1], a_fit_GoodLG, b_fit_GoodLG);

    //cout << "DEBUG3  :  " << a_fit_TDC_Gap_Fibers << "   " << b_fit_TDC_Gap_Fibers << "   " << R_TOF1 << endl;

    x_int_TDC_Gap_Fibers[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);
    x_int_TDC_Gap_Fibers[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);   
    y_int_TDC_Gap_Fibers[0] = y1_int(x_int_TDC_Gap_Fibers[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers); 
    y_int_TDC_Gap_Fibers[1] = y2_int(x_int_TDC_Gap_Fibers[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);

    x_int_TDC_Gap_Fibers_SFT[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_SFT_L1);
    x_int_TDC_Gap_Fibers_SFT[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_SFT_L1);   
    y_int_TDC_Gap_Fibers_SFT[0] = y1_int(x_int_TDC_Gap_Fibers_SFT[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers); 
    y_int_TDC_Gap_Fibers_SFT[1] = y2_int(x_int_TDC_Gap_Fibers_SFT[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);


    /// Selection of the Good Intersect Coordinates
    //////////////////////////////////////////////////////////////////////
  
    x_TDC_sel_intersect1 = 0.;   y_TDC_sel_intersect1 = 0.;
    x_GoodLG_intersect1 = 0.;    y_GoodLG_intersect1 = 0.;
    x_TDC_Gap_Fibers = 0.;       y_TDC_Gap_Fibers = 0.;
    //x_TDC_Gap_Fibers_SFT = 0.;   y_TDC_Gap_Fibers_SFT = 0.;

    for(int i=0; i<2; i++){
      dist1_TDC_selected[i] = 0.;
      dist1_GoodLG[i] = 0.;
      dist1_TDC_Gap_Fibers[i] = 0.;
      dist1_TDC_Gap_Fibers_SFT[i] = 0.;
    }

    for(int i=0; i<2; i++){
      dist1_TDC_selected[i] = distance(x_int_TDC_selected[i], y_int_TDC_selected[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_GoodLG[i] = distance(x_int_GoodLG[i], y_int_GoodLG[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
    }    
  
    if(dist1_TDC_selected[0] < dist1_TDC_selected[1])
    {
      x_TDC_sel_intersect1 = x_int_TDC_selected[0];
      y_TDC_sel_intersect1 = y_int_TDC_selected[0];
    }
    else if(dist1_TDC_selected[1] < dist1_TDC_selected[0])
    {
      x_TDC_sel_intersect1 = x_int_TDC_selected[1];
      y_TDC_sel_intersect1 = y_int_TDC_selected[1];
    }
    //else cout << "ERROR !" << endl;


    if(dist1_GoodLG[0] < dist1_GoodLG[1])
    {
      x_GoodLG_intersect1 = x_int_GoodLG[0];
      y_GoodLG_intersect1 = y_int_GoodLG[0];
    }
    else if(dist1_GoodLG[1] < dist1_GoodLG[0])
    {
      x_GoodLG_intersect1 = x_int_GoodLG[1];
      y_GoodLG_intersect1 = y_int_GoodLG[1];
    }
    //else cout << "ERROR !" << endl;

  
    if(dist1_TDC_Gap_Fibers[0] < dist1_TDC_Gap_Fibers[1])
    {
      x_TDC_Gap_Fibers = x_int_TDC_Gap_Fibers[0];
      y_TDC_Gap_Fibers = y_int_TDC_Gap_Fibers[0];
    }
    else if(dist1_TDC_Gap_Fibers[1] < dist1_TDC_Gap_Fibers[0])
    {
      x_TDC_Gap_Fibers = x_int_TDC_Gap_Fibers[1];
      y_TDC_Gap_Fibers = y_int_TDC_Gap_Fibers[1];
    }
    //else cout << "ERROR !" << endl;

    //if(dist1_TDC_Gap_Fibers_SFT[0] < dist1_TDC_Gap_Fibers_SFT[1])
    //{
      //x_TDC_Gap_Fibers_SFT = x_int_TDC_Gap_Fibers_SFT[0];
      //y_TDC_Gap_Fibers_SFT = y_int_TDC_Gap_Fibers_SFT[0];
    //}
    //else if(dist1_TDC_Gap_Fibers_SFT[1] < dist1_TDC_Gap_Fibers_SFT[0])
    //{
      //x_TDC_Gap_Fibers_SFT = x_int_TDC_Gap_Fibers_SFT[1];
      //y_TDC_Gap_Fibers_SFT = y_int_TDC_Gap_Fibers_SFT[1];
    //}
    //else cout << "ERROR !" << endl;

    /// Selection of the Good TOF1 Section
    //////////////////////////////////////////////////////////////////////
    for(int i=0; i<3; i++){
      dist2_TDC_selected[i] = 0.;
      dist2_GoodLG[i] = 0.;
    }

    dist2_TDC_selected_min = 1000.;
    dist2_GoodLG_min = 1000.;
    selected_TDC_selected = 0;


    for(int ii=0; ii<3; ii++)
    {
      dist2_TDC_selected[ii] = distance(x_TDC_Gap_Fibers, y_TDC_Gap_Fibers, Gap[gap_to_fit-1][ii][0], Gap[gap_to_fit-1][ii][1]);
      dist2_GoodLG[ii] = distance(x_GoodLG_intersect1, y_GoodLG_intersect1, Gap[gap_to_fit-1][ii][0], Gap[gap_to_fit-1][ii][1]);

      if(dist2_TDC_selected[ii] <= dist2_TDC_selected_min)
      {
        dist2_TDC_selected_min = dist2_TDC_selected[ii];
        selected_TDC_selected = ii;
      }

      if(dist2_GoodLG[ii] <= dist2_GoodLG_min)
      {
        dist2_GoodLG_min = dist2_GoodLG[ii];
      }
    }


    h_TOF1_closest->Fill(Gap[gap_to_fit_rotate-1][selected_TDC_selected][0], Gap[gap_to_fit_rotate-1][selected_TDC_selected][1]);

  

    //////////////////////////////////////////////////////////////////////


    h_int_TDC->Fill(x_int_TDC[0],y_int_TDC[0]);
    h_int_TDC->Fill(x_int_TDC[1],y_int_TDC[1]);


    h_int_TDC_selected->Fill(x_TDC_sel_intersect1, y_TDC_sel_intersect1);

    h_int_GoodLG->Fill(x_GoodLG_intersect1, y_GoodLG_intersect1);

    /// Add Weight on TOF1
    //////////////////////////////////////////////////////////////////
    for(int i = 0; i<5; i++)
      h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected][0], Gap[gap_to_fit-1][selected_TDC_selected][1]);
   
    for(int i = 0; i<5; i++)
      h_GoodLG_weighted->Fill(Gap[gap_to_fit-1][selected_TDC_selected][0], Gap[gap_to_fit-1][selected_TDC_selected][1]);
  
    for(int i = 0; i<3; i++)
      h_TDC_Gap_Fibers->Fill(Gap[gap_to_fit-1][selected_TDC_selected][0], Gap[gap_to_fit-1][selected_TDC_selected][1]);
    //////////////////////////////////////////////////////////////////


    h_int_TDC->SetMarkerStyle(20);            h_int_GoodLG->SetMarkerStyle(20);      h_int_GoodLG_weighted->SetMarkerStyle(20);
    h_int_TDC->SetMarkerColor(4);             h_int_GoodLG->SetMarkerColor(4);       h_int_GoodLG_weighted->SetMarkerColor(4);
    h_int_TDC->SetMarkerSize(0.8);            h_int_GoodLG->SetMarkerSize(0.8);      h_int_GoodLG_weighted->SetMarkerSize(0.8);

    h_int_TDC_selected->SetMarkerStyle(20);   h_int_TDC_selected_weighted->SetMarkerStyle(20);
    h_int_TDC_selected->SetMarkerColor(4);    h_int_TDC_selected_weighted->SetMarkerColor(4);
    h_int_TDC_selected->SetMarkerSize(0.8);   h_int_TDC_selected_weighted->SetMarkerSize(0.8);

    h_int_TDC_Gap_Fibers->SetMarkerStyle(20);   h_int_TDC_Gap_Fibers_SFT->SetMarkerStyle(20);
    h_int_TDC_Gap_Fibers->SetMarkerColor(4);    h_int_TDC_Gap_Fibers_SFT->SetMarkerColor(4);
    h_int_TDC_Gap_Fibers->SetMarkerSize(0.8);   h_int_TDC_Gap_Fibers_SFT->SetMarkerSize(0.8);
  
    h_int_TDC_Gap_Fibers_rotate->SetMarkerStyle(20);   h_int_TDC_Gap_Fibers_SFT_rotate->SetMarkerStyle(20);
    h_int_TDC_Gap_Fibers_rotate->SetMarkerColor(4);    h_int_TDC_Gap_Fibers_SFT_rotate->SetMarkerColor(4);
    h_int_TDC_Gap_Fibers_rotate->SetMarkerSize(0.8);   h_int_TDC_Gap_Fibers_SFT_rotate->SetMarkerSize(0.8);

    h_int_TDC_SFT->SetMarkerStyle(20);        h_int_GoodLG_SFT->SetMarkerStyle(20);
    h_int_TDC_SFT->SetMarkerColor(1);         h_int_GoodLG_SFT->SetMarkerColor(1);
    h_int_TDC_SFT->SetMarkerSize(0.8);        h_int_GoodLG_SFT->SetMarkerSize(0.8);

    h_int_TDC_SFT->SetMarkerStyle(20);
    h_int_TDC_SFT->SetMarkerColor(1);
    h_int_TDC_SFT->SetMarkerSize(0.8);

    ParError = 999.99;

    h_TDC_selected2->Fit("pol1", "QCM");
    fit_line_TDC_selected2 = h_TDC_selected2->GetFunction("pol1");
    fit_line_TDC_selected2->SetLineWidth(2);
    fit_line_TDC_selected2->SetLineColor(2);
    a_fit_TDC_selected_weighted=fit_line_TDC_selected2->GetParameter(1);      
    b_fit_TDC_selected_weighted=fit_line_TDC_selected2->GetParameter(0);

    h_GoodLG_weighted->Fit("pol1", "QCM");
    fit_line_GoodLG_weighed = h_GoodLG_weighted->GetFunction("pol1");
    fit_line_GoodLG_weighed->SetLineWidth(2);
    fit_line_GoodLG_weighed->SetLineColor(1);
    a_fit_GoodLG_weighted=fit_line_GoodLG_weighed->GetParameter(1);      
    b_fit_GoodLG_weighted=fit_line_GoodLG_weighed->GetParameter(0);

    h_TDC_Gap_Fibers->Fit("pol1", "QCM");
    fit_line_TDC_Gap_Fibers = h_TDC_Gap_Fibers->GetFunction("pol1");
    fit_line_TDC_Gap_Fibers->SetLineWidth(2);
    fit_line_TDC_Gap_Fibers->SetLineColor(2);
    a_fit_TDC_Gap_Fibers=fit_line_TDC_Gap_Fibers->GetParameter(1);      
    b_fit_TDC_Gap_Fibers=fit_line_TDC_Gap_Fibers->GetParameter(0);
    ParError=fit_line_TDC_Gap_Fibers->GetParError(1);
    ChiS = fit_line_TDC_Gap_Fibers->GetChisquare();
    ndf = fit_line_TDC_Gap_Fibers->GetNDF();
    prob = fit_line_TDC_Gap_Fibers->GetProb();

    x_int_TDC_selected_weighted[0] = intersectx1(a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted, R_TOF1);
    x_int_TDC_selected_weighted[1] = intersectx2(a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted, R_TOF1); 
    x_int_TDC_SFT[0] = intersectx1(a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted, R_SFT_L1);
    x_int_TDC_SFT[1] = intersectx2(a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted, R_SFT_L1);

    y_int_TDC_selected_weighted[0] = y1_int(x_int_TDC_selected_weighted[0], a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted);
    y_int_TDC_selected_weighted[1] = y2_int(x_int_TDC_selected_weighted[1], a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted);
    y_int_TDC_SFT[0] = y1_int(x_int_TDC_SFT[0], a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted);
    y_int_TDC_SFT[1] = y2_int(x_int_TDC_SFT[1], a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted);

    x_int_GoodLG_weighted[0] = intersectx1(a_fit_GoodLG_weighted, b_fit_GoodLG_weighted, R_TOF1);
    x_int_GoodLG_weighted[1] = intersectx2(a_fit_GoodLG_weighted, b_fit_GoodLG_weighted, R_TOF1);
    x_int_GoodLG_SFT[0] = intersectx1(a_fit_GoodLG_weighted, b_fit_GoodLG_weighted, R_SFT_L1);
    x_int_GoodLG_SFT[1] = intersectx2(a_fit_GoodLG_weighted, b_fit_GoodLG_weighted, R_SFT_L1);

    y_int_GoodLG_weighted[0] = y1_int(x_int_GoodLG_weighted[0], a_fit_GoodLG_weighted, b_fit_GoodLG_weighted);
    y_int_GoodLG_weighted[1] = y2_int(x_int_GoodLG_weighted[1], a_fit_GoodLG_weighted, b_fit_GoodLG_weighted);
    y_int_GoodLG_SFT[0] = y1_int(x_int_GoodLG_SFT[0], a_fit_GoodLG_weighted, b_fit_GoodLG_weighted);
    y_int_GoodLG_SFT[1] = y2_int(x_int_GoodLG_SFT[1], a_fit_GoodLG_weighted, b_fit_GoodLG_weighted);


    x_int_TDC_Gap_Fibers[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);
    x_int_TDC_Gap_Fibers[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);   
    x_int_TDC_Gap_Fibers_SFT[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_SFT_L1);
    x_int_TDC_Gap_Fibers_SFT[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_SFT_L1);   

    y_int_TDC_Gap_Fibers[0] = y1_int(x_int_TDC_Gap_Fibers[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers); 
    y_int_TDC_Gap_Fibers[1] = y2_int(x_int_TDC_Gap_Fibers[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);
    y_int_TDC_Gap_Fibers_SFT[0] = y1_int(x_int_TDC_Gap_Fibers_SFT[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers); 
    y_int_TDC_Gap_Fibers_SFT[1] = y2_int(x_int_TDC_Gap_Fibers_SFT[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);

    x_int_TARGET[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TARGET);
    x_int_TARGET[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TARGET);
    y_int_TARGET[0] = y1_int(x_int_TARGET[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);
    y_int_TARGET[1] = y2_int(x_int_TARGET[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);


    x_TDC_selected_weighted_intersect1 = 0.;   y_TDC_selected_weighted_intersect1 = 0.;
    x_GoodLG_weighted_intersect1 = 0.;         y_GoodLG_weighted_intersect1 = 0.;
    x_TDC_SFT_intersect1 = 0.;                 y_TDC_SFT_intersect1 = 0.;
    x_GoodLG_SFT_intersect1 = 0.;              y_GoodLG_SFT_intersect1 = 0.;
    x_TDC_Gap_Fibers_intersect1 = 0.;          y_TDC_Gap_Fibers_intersect1 = 0.;
    x_TDC_Gap_Fibers_SFT_intersect1 = 0.;      y_TDC_Gap_Fibers_SFT_intersect1 = 0.;
    x_TARGET_intersect = 0;                    y_TARGET_intersect = 0;
 
    //cout << " DEBUG 1 :  " << x_TDC_Gap_Fibers_intersect1 << endl;

    for(int i=0; i<2; i++){
      dist1_TDC_selected_weighted[i] = 0.;
      dist1_TDC_SFT[i] = 0.;
      dist1_GoodLG_weighted[i] = 0.;
      dist1_GoodLG_SFT[i] = 0.;
      dist1_TARGET_intersect[i] = 0.;
    }

    for(int i=0; i<2; i++)
    {
      dist1_TDC_selected_weighted[i] = distance(x_int_TDC_selected_weighted[i], y_int_TDC_selected_weighted[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TDC_SFT[i] = distance(x_int_TDC_SFT[i], y_int_TDC_SFT[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_GoodLG_weighted[i] = distance(x_int_GoodLG_weighted[i], y_int_GoodLG_weighted[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_GoodLG_SFT[i] = distance(x_int_GoodLG_SFT[i], y_int_GoodLG_SFT[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TARGET_intersect[i] = distance(x_int_TARGET[i], y_int_TARGET[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
    }    
  
    if(dist1_TDC_selected_weighted[0] < dist1_TDC_selected_weighted[1])
    {
      x_TDC_selected_weighted_intersect1 = x_int_TDC_selected_weighted[0];
      y_TDC_selected_weighted_intersect1 = y_int_TDC_selected_weighted[0];
    }
    else if(dist1_TDC_selected_weighted[1] < dist1_TDC_selected_weighted[0])
    {
      x_TDC_selected_weighted_intersect1 = x_int_TDC_selected_weighted[1];
      y_TDC_selected_weighted_intersect1 = y_int_TDC_selected_weighted[1];
    }
    //else cout << "ERROR !" << endl;


    if(dist1_GoodLG_weighted[0] < dist1_GoodLG_weighted[1])
    {
      x_GoodLG_weighted_intersect1 = x_int_GoodLG_weighted[0];
      y_GoodLG_weighted_intersect1 = y_int_GoodLG_weighted[0];
    }
    else if(dist1_GoodLG_weighted[1] < dist1_GoodLG_weighted[0])
    {
      x_GoodLG_weighted_intersect1 = x_int_GoodLG_weighted[1];
      y_GoodLG_weighted_intersect1 = y_int_GoodLG_weighted[1];
    }
    //else cout << "ERROR !" << endl;


    if(dist1_TDC_SFT[0] < dist1_TDC_SFT[1])
    {
      x_TDC_SFT_intersect1 = x_int_TDC_SFT[0];
      y_TDC_SFT_intersect1 = y_int_TDC_SFT[0];
    }
    else if(dist1_TDC_SFT[1] < dist1_TDC_SFT[0])
    {
      x_TDC_SFT_intersect1 = x_int_TDC_SFT[1];
      y_TDC_SFT_intersect1 = y_int_TDC_SFT[1];
    }
    //else cout << "ERROR !" << endl;


    if(dist1_GoodLG_SFT[0] < dist1_GoodLG_SFT[1])
    {
      x_GoodLG_SFT_intersect1 = x_int_GoodLG_SFT[0];
      y_GoodLG_SFT_intersect1 = y_int_GoodLG_SFT[0];
    }
    else if(dist1_GoodLG_SFT[1] < dist1_GoodLG_SFT[0])
    {
      x_GoodLG_SFT_intersect1 = x_int_GoodLG_SFT[1];
      y_GoodLG_SFT_intersect1 = y_int_GoodLG_SFT[1];
    }
    //else cout << "ERROR !" << endl;

    if(dist1_TDC_Gap_Fibers[0] < dist1_TDC_Gap_Fibers[1])
    {
      x_TDC_Gap_Fibers_intersect1 = x_int_TDC_Gap_Fibers[0];
      y_TDC_Gap_Fibers_intersect1 = y_int_TDC_Gap_Fibers[0];
      //cout << " DEBUG 1 :  " << x_TDC_Gap_Fibers_intersect1 << endl;
      //cout << " DEBUG 1 :  " << x_int_TDC_Gap_Fibers[0] << endl;
    }
    else if(dist1_TDC_Gap_Fibers[1] < dist1_TDC_Gap_Fibers[0])
    {
      x_TDC_Gap_Fibers_intersect1 = x_int_TDC_Gap_Fibers[1];
      y_TDC_Gap_Fibers_intersect1 = y_int_TDC_Gap_Fibers[1];
    }
    //else cout << "ERROR !" << endl;

    if(dist1_TDC_Gap_Fibers_SFT[0] < dist1_TDC_Gap_Fibers_SFT[1])
    {
      x_TDC_Gap_Fibers_SFT_intersect1 = x_int_TDC_Gap_Fibers_SFT[0];
      y_TDC_Gap_Fibers_SFT_intersect1 = y_int_TDC_Gap_Fibers_SFT[0];
    }
    else if(dist1_TDC_Gap_Fibers_SFT[1] < dist1_TDC_Gap_Fibers_SFT[0])
    {
      x_TDC_Gap_Fibers_SFT_intersect1 = x_int_TDC_Gap_Fibers_SFT[1];
      y_TDC_Gap_Fibers_SFT_intersect1 = y_int_TDC_Gap_Fibers_SFT[1];
    }
    //else cout << "ERROR !" << endl;


    if(dist1_TARGET_intersect[0] < dist1_TARGET_intersect[1])
    {
      x_TARGET_intersect = x_int_TARGET[1];
      y_TARGET_intersect = y_int_TARGET[1];
    }
    else if(dist1_TARGET_intersect[1] < dist1_TARGET_intersect[0])
    {
      x_TARGET_intersect = x_int_TARGET[0];
      y_TARGET_intersect = y_int_TARGET[0];
    }
    //else cout << "ERROR !" << endl;

 


    h_int_TDC_selected_weighted->Fill(x_TDC_selected_weighted_intersect1, y_TDC_selected_weighted_intersect1);
    h_int_GoodLG_weighted->Fill(x_GoodLG_weighted_intersect1, y_GoodLG_weighted_intersect1);
    h_int_TDC_SFT->Fill(x_TDC_SFT_intersect1, y_TDC_SFT_intersect1);
    h_int_GoodLG_SFT->Fill(x_GoodLG_SFT_intersect1, y_GoodLG_SFT_intersect1);
    
    if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
      h_int_TDC_Gap_Fibers_rotate->Fill(-y_TDC_Gap_Fibers_intersect1, x_TDC_Gap_Fibers_intersect1);
      h_int_TDC_Gap_Fibers_SFT_rotate->Fill(-y_TDC_Gap_Fibers_SFT_intersect1, x_TDC_Gap_Fibers_SFT_intersect1);   
    }
    else{
      h_int_TDC_Gap_Fibers->Fill(x_TDC_Gap_Fibers_intersect1, y_TDC_Gap_Fibers_intersect1);
      h_int_TDC_Gap_Fibers_SFT->Fill(x_TDC_Gap_Fibers_SFT_intersect1, y_TDC_Gap_Fibers_SFT_intersect1);   
    }
    h_int_TDC_TARGET->Fill(x_TARGET_intersect, y_TARGET_intersect);

    TLine *CenterLine_TDC = new TLine(0., 0., x_TDC_Gap_Fibers_intersect1, y_TDC_Gap_Fibers_intersect1);
    TLine *CenterLine_GoodLG = new TLine(0., 0., x_GoodLG_weighted_intersect1, y_GoodLG_weighted_intersect1);


    CenterLine_TDC->SetLineWidth(2);    CenterLine_GoodLG->SetLineWidth(2);   
    CenterLine_TDC->SetLineColor(3);    CenterLine_GoodLG->SetLineColor(3);
    HorizontalAxis->SetLineWidth(2);    HorizontalAxis->SetLineColor(1);


    /// Distance from K-stop bar to best fit line
    ///////////////////////////////////////////////////////////////////////////////
    // X_BAR, Y_BAR;

    dist_to_k_stop = 0.;

    if(X_BAR != -10000){
      if((gap_to_fit_rotate==6 || gap_to_fit_rotate==12 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
        dist_to_k_stop = distance_to_line(Y_BAR,-X_BAR,a_fit_TDC_Gap_Fibers,b_fit_TDC_Gap_Fibers);
      }
      else{
        dist_to_k_stop = distance_to_line(X_BAR,Y_BAR,a_fit_TDC_Gap_Fibers,b_fit_TDC_Gap_Fibers);
      }
    }
    else{
      dist_to_k_stop = -10000;
    }

    /// Angle Calculation
    ///////////////////////////////////////////////////////////////////////////////

    a_final_TDC = 0.;           a_final_GoodLG = 0.;             a_final_Gap_Fibers = 0.;       a_final_guide = 0.;
    alpha_TDC = 0.;             alpha_GoodLG = 0.;               alpha_Gap_Fibers = 0.;         alpha_guide = 0.;
    tanalpha_TDC = 0.;          tanalpha_GoodLG = 0.;            tanalpha_Gap_Fibers = 0.;      tanalpha_guide = 0.;
    angle_final_TDC = 0.;       angle_final_GoodLG = 0.;         angle_final_Gap_Fibers = 0.;   angle_final_guide = 0.;

    a_final_TDC = y_TDC_selected_weighted_intersect1 / x_TDC_selected_weighted_intersect1;
    a_final_GoodLG = y_GoodLG_weighted_intersect1 / x_GoodLG_weighted_intersect1;
    a_final_Gap_Fibers = y_TDC_Gap_Fibers_intersect1 / x_TDC_Gap_Fibers_intersect1;
    a_final_guide = (y_TDC_Gap_Fibers_intersect1 - y_TARGET_intersect) / (x_TDC_Gap_Fibers_intersect1 - x_TARGET_intersect);

    tanalpha_TDC = a_final_TDC;       tanalpha_GoodLG = a_final_GoodLG;       tanalpha_Gap_Fibers = a_final_Gap_Fibers;
    alpha_TDC = atan(tanalpha_TDC);   alpha_GoodLG = atan(tanalpha_GoodLG);   alpha_Gap_Fibers = atan(tanalpha_Gap_Fibers);

    tanalpha_guide = a_final_guide;
    alpha_guide = atan(tanalpha_guide);

    //Determination of Delta Phi
    Delta_phi = 999.99;  Delta_phi_deg = 999.99;
    Delta_phi = sin(alpha_guide)*cos(alpha_guide)*(ParError/a_final_guide);
    Delta_phi_deg = (180/PI)*Delta_phi;
    //cout << "DT:  " << Delta_phi_deg << endl; 


    Axis_Vector_Length = 10;
    TArrow *x_guide;
    TArrow *y_guide;
    TArrow *x_guide_rotate;
    TArrow *y_guide_rotate;


    x_guide = new TArrow(x_TARGET_intersect, y_TARGET_intersect, x_TARGET_intersect + Axis_Vector_Length, y_TARGET_intersect, 0.005, "|>");
    y_guide = new TArrow(x_TARGET_intersect, y_TARGET_intersect, x_TARGET_intersect, y_TARGET_intersect + Axis_Vector_Length, 0.005, "|>");

    x_guide->SetLineWidth(2);     y_guide->SetLineWidth(2);
    x_guide->SetLineColor(4);     y_guide->SetLineColor(4);
  
    x_guide_rotate = new TArrow(-y_TARGET_intersect, x_TARGET_intersect, -y_TARGET_intersect + Axis_Vector_Length, x_TARGET_intersect, 0.005, "|>");
    y_guide_rotate = new TArrow(-y_TARGET_intersect, x_TARGET_intersect, -y_TARGET_intersect, x_TARGET_intersect + Axis_Vector_Length, 0.005, "|>");

    x_guide_rotate->SetLineWidth(2);     y_guide_rotate->SetLineWidth(2);
    x_guide_rotate->SetLineColor(4);     y_guide_rotate->SetLineColor(4);    


    /// ###  KEITO'S ANALYSIS  ###

    x_Kstop = 0.;     x_SFT1 = 0.;
    y_Kstop = 0.;     y_SFT1 = 0.;
    x2_Keito = 0.;    y2_Keito = 0.;
    x3_Keito = 0.;    y3_Keito = 0.;

    length = 20;
    angle_Keito = 0.;

    //a_final_Keito = 0.;


    TArrow *x_guide_Keito;
    TArrow *y_guide_Keito;


    TH2F *h_Keito = new TH2F("h_Keito", "h_Keito", 500, -50, 50, 500, -50, 50);
    TH2F *h_help_point_Keito = new TH2F("h_help_point_Keito", "h_help_point_Keito", 500, -50, 50, 500, -50, 50);
    TH2F *h_help_point_Keito2 = new TH2F("h_help_point_Keito2", "h_help_point_Keito2", 500, -50, 50, 500, -50, 50);
    //TF1 *Track_Keito;
    //TF1 *Track_Keito2;

    h_Keito->SetMarkerStyle(21);      h_help_point_Keito->SetMarkerStyle(6);        h_help_point_Keito2->SetMarkerStyle(6);
    h_Keito->SetMarkerSize(0.8);      h_help_point_Keito->SetMarkerSize(1000);      h_help_point_Keito2->SetMarkerSize(1000); 
    h_Keito->SetMarkerColor(6);       h_help_point_Keito->SetMarkerColor(6);        h_help_point_Keito2->SetMarkerColor(6); 


    x_guide_Keito = new TArrow(x_SFT1, y_SFT1, x_SFT1 + Axis_Vector_Length, y_SFT1, 0.005, "|>");
    y_guide_Keito = new TArrow(x_SFT1, y_SFT1, x_SFT1, y_SFT1 + Axis_Vector_Length, 0.005, "|>");
  
    x_guide_Keito->SetLineWidth(2);     y_guide_Keito->SetLineWidth(2);
    x_guide_Keito->SetLineColor(4);     y_guide_Keito->SetLineColor(4);


    x2_Keito = x_SFT1 + length * cos(angle_Keito);
    y2_Keito = y_SFT1 + length * sin(angle_Keito);

    x3_Keito = x_Kstop + length * cos(angle_Keito);
    y3_Keito = y_Kstop + length * sin(angle_Keito);


    //a_final_Keito = y_SFT1 / x_SFT1;


    h_Keito->Fill(x_Kstop, y_Kstop);
    h_Keito->Fill(x_SFT1, y_SFT1);

    h_help_point_Keito->Fill(x_SFT1, y_SFT1);
    h_help_point_Keito->Fill(x2_Keito, y2_Keito);

    h_help_point_Keito2->Fill(x_Kstop, y_Kstop);
    h_help_point_Keito2->Fill(x3_Keito, y3_Keito);


    //### Option 2
    if(x_TDC_selected_weighted_intersect1 < 0) angle_final_TDC = 180. + alpha_TDC * (180./PI);

    if(x_TDC_selected_weighted_intersect1 >= 0){
      if(y_TDC_selected_weighted_intersect1 >= 0) angle_final_TDC = alpha_TDC * (180./PI);
      else angle_final_TDC = alpha_TDC * (180./PI) + 360.;
    }


    if(x_GoodLG_weighted_intersect1 < 0) angle_final_GoodLG = 180. + alpha_GoodLG * (180./PI);

    if(x_GoodLG_weighted_intersect1 >= 0){
      if(y_GoodLG_weighted_intersect1 >= 0) angle_final_GoodLG = alpha_GoodLG * (180./PI);
      else angle_final_GoodLG = alpha_GoodLG * (180./PI) + 360.;
    }


    if(x_TDC_Gap_Fibers_intersect1 < 0)  angle_final_Gap_Fibers = 180. + alpha_Gap_Fibers * (180./PI);

    if(x_TDC_Gap_Fibers_intersect1 >= 0){
      if(y_TDC_Gap_Fibers_intersect1 >= 0)  angle_final_Gap_Fibers = alpha_Gap_Fibers * (180./PI);
      else angle_final_Gap_Fibers = alpha_Gap_Fibers * (180./PI) + 360.;
    }
  
    //cout << "DEBUG :  " << x_TDC_Gap_Fibers_intersect1 << "   " <<  x_TARGET_intersect << "   " << angle_final_guide << "   " << alpha_guide << endl;
    if((x_TDC_Gap_Fibers_intersect1 - x_TARGET_intersect) < 0) angle_final_guide = 180. + alpha_guide * (180./PI);

    if((x_TDC_Gap_Fibers_intersect1 - x_TARGET_intersect) >= 0){
      if((y_TDC_Gap_Fibers_intersect1 - y_TARGET_intersect) >= 0)  angle_final_guide = alpha_guide * (180./PI);
      else angle_final_guide = alpha_guide * (180./PI) + 360.;
    }
    
    //cout << "DEBUG :  " << ivt << "   " << angle_final_guide << endl;
    //cout << "DEBUG :  " << x_TDC_Gap_Fibers_intersect1 << "   "
    //                  << x_TARGET_intersect << "   " 
    //                  << y_TDC_Gap_Fibers_intersect1 << "   "
    //                  << y_TARGET_intersect << "   "
    //                  << alpha_guide << "   " 
    //                  << angle_final_guide << endl;

  
    if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1) angle_final_guide += 90.;
    if(angle_final_guide > 360.) angle_final_guide -= 360.;

    //cout << "DEBUG :  " << ivt << "   " << angle_final_guide << endl;

    sprintf(Angle_TDC_string,"#phi = %3.2f#circ", angle_final_TDC);
    sprintf(Angle_GoodLD_string,"#phi = %3.2f#circ", angle_final_GoodLG);
    sprintf(Angle_Gap_Fibers_string,"#phi = %3.2f#circ  (OLD)", angle_final_Gap_Fibers);
    //sprintf(Angle_guide_string,"#phi = %3.2f#circ", angle_final_guide);
    sprintf(Angle_guide_string,"#phi = %3.2f#circ#pm%3.2f", angle_final_guide, Delta_phi_deg);
    sprintf(Angle_Keito_string,"#phi = %3.2f#circ (Keito)", angle_Keito * (180/PI));


    TLatex *tex_Angle_TDC;
    TLatex *tex_Angle_GoodLG;
    TLatex *tex_Angle_Gap_Fibers;
    TLatex *tex_Angle_guide;
    TLatex *tex_Angle_Keito;

    tex_Angle_TDC = new TLatex(-45.,43.,Angle_TDC_string);
    tex_Angle_TDC->SetTextSize(0.05);
    tex_Angle_TDC->SetLineWidth(2);

    tex_Angle_GoodLG = new TLatex(-45.,43.,Angle_GoodLD_string);
    tex_Angle_GoodLG->SetTextSize(0.05);
    tex_Angle_GoodLG->SetLineWidth(2);

    tex_Angle_Gap_Fibers = new TLatex(-45.,43.,Angle_Gap_Fibers_string);
    tex_Angle_Gap_Fibers->SetTextSize(0.05);
    tex_Angle_Gap_Fibers->SetLineWidth(2);

    tex_Angle_guide = new TLatex(-45.,43.,Angle_guide_string);
    //tex_Angle_guide->SetTextSize(0.05);
    tex_Angle_guide->SetTextSize(0.045);
    tex_Angle_guide->SetLineWidth(2);

    tex_Angle_Keito = new TLatex(-45.,37.,Angle_Keito_string);
    tex_Angle_Keito->SetTextSize(0.05);
    tex_Angle_Keito->SetLineWidth(2);

    ///////////////////////////////////////////////////////////////////////////////


    TLine *best_fit_rotate = new TLine(50*a_fit_TDC_Gap_Fibers - b_fit_TDC_Gap_Fibers,-50,-50*a_fit_TDC_Gap_Fibers - b_fit_TDC_Gap_Fibers,50); //ROTATE_CHANGE
    best_fit_rotate->SetLineWidth(2); //ROTATE_CHANGE
    best_fit_rotate->SetLineColor(kRed); // ROTATE_CHANGE
    h_Centroid->Fill(X_BAR,Y_BAR);

    //TGaxis *A1 = new TGaxis(-50,50,50,50,"pol1",510,"-U");
    //TGaxis *A2 = new TGaxis(50,-50,50,50,"pol1",510,"+U");
  
  
  
    TLatex *x_sft;
    TLatex *y_sft;
    TLatex *x_tof1;
    TLatex *y_tof1;    
  
  
   
    if((gap_to_fit_rotate==6 || gap_to_fit_rotate==12 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
      sprintf(X_SFT_String,"X(40) = %4.2f", -y_TDC_Gap_Fibers_SFT_intersect1);
      sprintf(Y_SFT_String,"Y(40) = %4.2f", x_TDC_Gap_Fibers_SFT_intersect1);
      sprintf(X_TOF1_String,"X(47.1) = %4.2f", -y_TDC_Gap_Fibers_intersect1);
      sprintf(Y_TOF1_String,"Y(47.1) = %4.2f", x_TDC_Gap_Fibers_intersect1);    
    }
    else{
      sprintf(X_SFT_String,"X(40) = %4.2f", x_TDC_Gap_Fibers_SFT_intersect1);
      sprintf(Y_SFT_String,"Y(40) = %4.2f", y_TDC_Gap_Fibers_SFT_intersect1);
      sprintf(X_TOF1_String,"X(47.1) = %4.2f", x_TDC_Gap_Fibers_intersect1);
      sprintf(Y_TOF1_String,"Y(47.1) = %4.2f", y_TDC_Gap_Fibers_intersect1);
    }

    x_sft = new TLatex(10.,42.,X_SFT_String);  
    x_sft->SetTextSize(0.05);
    x_sft->SetLineWidth(2);

    y_sft = new TLatex(10.,37.,Y_SFT_String);
    y_sft->SetTextSize(0.05);
    y_sft->SetLineWidth(2);
  
  
    x_tof1 = new TLatex(-45.,-40.,X_TOF1_String);  
    x_tof1->SetTextSize(0.05);
    x_tof1->SetLineWidth(2);

    y_tof1 = new TLatex(-45.,-45.,Y_TOF1_String);
    y_tof1->SetTextSize(0.05);
    y_tof1->SetLineWidth(2);   


    // K Stop line / length in target

    TLine *k_stop_line;

    if(X_BAR != -10000){
      if(angle_final_guide >= 90 && angle_final_guide <= 270){ // second x point in the left
        k_stop_line = new TLine(X_BAR, Y_BAR, -50, tan(angle_final_guide*M_PI/180.0)*(-50-X_BAR) + Y_BAR);
      }
      else{ // second x point in the right
        k_stop_line = new TLine(X_BAR, Y_BAR, 50, tan(angle_final_guide*M_PI/180.0)*(50-X_BAR) + Y_BAR);
      }
    }
    else{
      k_stop_line = new TLine(0,0,0,0);
    }

    k_stop_line->SetLineColor(3);
    k_stop_line->SetLineWidth(2);

    //length_in_target = 0.;


    if(X_BAR != -10000){
      x_tof1_intersect_1 = 0.;
      y_tof1_intersect_1 = 0.;
      x_tof1_intersect_2 = 0.;
      y_tof1_intersect_2 = 0.;
      //x_tof1_intersect = 0.;
      //y_tof1_intersect = 0.;

      alpha = angle_final_guide - 90.0;

      if((gap_to_fit_rotate==6 || gap_to_fit_rotate==12 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
        x_tof1_intersect_1 = y_TDC_Gap_Fibers_intersect1 - dist_to_k_stop*cos(alpha*M_PI/180.0);
        x_tof1_intersect_2 = y_TDC_Gap_Fibers_intersect1 + dist_to_k_stop*cos(alpha*M_PI/180.0);
        y_tof1_intersect_1 = -x_TDC_Gap_Fibers_intersect1 - dist_to_k_stop*sin(alpha*M_PI/180.0);
        y_tof1_intersect_2 = -x_TDC_Gap_Fibers_intersect1 + dist_to_k_stop*sin(alpha*M_PI/180.0);

        if(distance(X_BAR,Y_BAR,x_tof1_intersect_1,y_tof1_intersect_1) < distance(X_BAR,Y_BAR,x_tof1_intersect_2,y_tof1_intersect_2)){
          //x_tof1_intersect = x_tof1_intersect_1;
          //y_tof1_intersect = y_tof1_intersect_1;
        }
        else{
          //x_tof1_intersect = x_tof1_intersect_2;
          //y_tof1_intersect = y_tof1_intersect_2;
        }
      }
      else{

        x_tof1_intersect_1 = x_TDC_Gap_Fibers_intersect1 - dist_to_k_stop*cos(alpha*M_PI/180.0);
        x_tof1_intersect_2 = x_TDC_Gap_Fibers_intersect1 + dist_to_k_stop*cos(alpha*M_PI/180.0);
        y_tof1_intersect_1 = y_TDC_Gap_Fibers_intersect1 - dist_to_k_stop*sin(alpha*M_PI/180.0);
        y_tof1_intersect_2 = y_TDC_Gap_Fibers_intersect1 + dist_to_k_stop*sin(alpha*M_PI/180.0);

        if(distance(X_BAR,Y_BAR,x_tof1_intersect_1,y_tof1_intersect_1) < distance(X_BAR,Y_BAR,x_tof1_intersect_2,y_tof1_intersect_2)){
          //x_tof1_intersect = x_tof1_intersect_1;
          //y_tof1_intersect = y_tof1_intersect_1;
        }
        else{
          //x_tof1_intersect = x_tof1_intersect_2;
          //y_tof1_intersect = y_tof1_intersect_2;
        }
      }

      //length_in_target = distance(X_BAR,Y_BAR,x_tof1_intersect,y_tof1_intersect); //CORRECT
    }
    else{
      //cout << "No K-Stop for track length." << endl;
    }

    Z_selected.clear();

    //ZZ.clear(); 
    //ZZ = GetZ(Run_Number, ivt, angle_final_guide, gap_to_fit_rotate, true, C2X_centroid, length_in_target);

    if(Switch_Zlist > 0){
      //DZ.clear();
      //DZ = GetDeltaZ(Run_Number, ivt, angle_final_guide, gap_to_fit_rotate, true, C2X_centroid, length_in_target);

      //Track_Length.clear();
      //Track_Length = GetTrackLength(Run_Number, ivt, angle_final_guide, gap_to_fit_rotate, true, C2X_centroid, length_in_target);
    }  
      
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //TOF1_Z_range = 25.; 
    //TOF1_Z_min = 999.;
    //TOF1_Z_max = 999.;

    Z_TOF1[gap_to_fit_rotate-1] = 
    Slope_correction_TOF1[gap_to_fit_rotate-1]*((0.025*150*0.5*(TDC_TOF1U[gap_to_fit_rotate-1] - TDC_TOF1D[gap_to_fit_rotate-1])) + Offset_TOF1[gap_to_fit_rotate-1]);

    //TOF1_Z_min = Z_TOF1[gap_to_fit_rotate-1] - TOF1_Z_range;
    //TOF1_Z_max = Z_TOF1[gap_to_fit_rotate-1] + TOF1_Z_range;


    //for(unsigned i=0; i<ZZ.size(); i++){
    //  if(ZZ[i]>=TOF1_Z_min && ZZ[i]<=TOF1_Z_max){
    //    Z_selected.push_back(ZZ[i]);
    //  }
    //}


    //SORT VECTOR ELEMENTS
    //Z_selected_sorted.empty();
  
    //sort(Z_selected.begin(), Z_selected.end());

    //for(unsigned i=0; i<Z_selected.size(); i++) cout << "DEBUG:  " << Z_selected[i] << endl;

    //cout << "TEST:  " << Run_Number << "  " << ivt << "  " << gap_to_fit_rotate << "  " << angle_final_guide << "  " << Delta_phi_deg << "  "
    //     << X_BAR << "  " << Y_BAR << "  " << Z_TOF1[gap_to_fit_rotate-1] << "  ";
      

      //Switch_Display = 1;

      if(Switch_Display>0){ 
        cout << fixed;
        cout << setw(4) << Run_Number << "  ";
        cout << setw(7) << ivt << "  ";
        cout << setw(2) << gap_to_fit_rotate << "  ";
        cout << setw(7) << setprecision(3) << angle_final_guide << "  ";
        cout << setw(6) << setprecision(3) << Delta_phi_deg << "  ";
        cout << setw(8) << setprecision(3) << ChiS/ndf << "  ";
        cout << setw(7) << setprecision(3) << X_BAR << "  ";
        cout << setw(7) << setprecision(3) << Y_BAR << "  ";
        cout << setw(9) << setprecision(3) << Z_TOF1[gap_to_fit_rotate-1] << "  ";

        //for(unsigned i=0; i<Z_selected.size(); i++) cout << setw(7) << Z_selected[i] << "  ";
      }


    fout << fixed;
    fout << setw(4) << Run_Number << "  ";
    fout << setw(7) << ivt << "  ";
    fout << setw(2) << gap_to_fit_rotate << "  ";
    fout << setw(7) << setprecision(3) << angle_final_guide << "  ";
    fout << setw(6) << setprecision(3) << Delta_phi_deg << "  ";
    fout << setw(8) << setprecision(3) << ChiS/ndf << "  ";
    fout << setw(7) << setprecision(3) << X_BAR << "  ";
    fout << setw(7) << setprecision(3) << Y_BAR << "  ";
    fout << setw(9) << setprecision(3) << Z_TOF1[gap_to_fit_rotate-1] << "  ";

    //for(unsigned i=0; i<Z_selected.size(); i++)  fout << setw(7) << Z_selected[i] << "  ";

    fout << endl;

    //if(Switch_Zlist > 0){
    //  fout2 << fixed;
    //  fout2 << setw(4) << Run_Number << "  ";
    //  fout2 << setw(7) << ivt << "  ";

      //for(unsigned i=0; i<ZZ.size(); i++){
        //fout2 << setw(7) << ZZ[i] << "  ";
        //fout2 << setw(5) << DZ[i] << "  ";
      //}

    //  fout2 << endl;
    //}

    if(Switch_Display>0) cout << endl;

    //h_ADC_L1_DS->Reset();
    //h_ADC_L2_DS->Reset();
    //h_ADC_L3_DS->Reset();
    //h_ADC_L4_DS->Reset();

    //h_ADC_L1_US->Reset();
    //h_ADC_L2_US->Reset();
    //h_ADC_L3_US->Reset();
    //h_ADC_L4_US->Reset();



    h_Circle->Reset();

    h_target_ADC->Reset();
    h_target_TDC->Reset();
    h_target_TDC_copy->Reset();
    h_TDC_selected->Reset();
    h_TDC_selected2->Reset();

    h_target_ADC2->Reset();
    h_target_TDC2->Reset();

    h_target_ADC3->Reset();
    h_target_ADC4->Reset();
    h_target_ADCA->Reset();

    h_Target->Reset();
    h_Target_Center->Reset();

    h_kaon->Reset();
    h_kaon_copy->Reset();
    
    h_max->Reset();
    h_max_copy->Reset();

    h_GoodLG->Reset();
    h_GoodLG_weighted->Reset();
    h_GoodLG_copy->Reset();

    h_TOF1->Reset();
    h_TOF1_closest->Reset();
    h_TOF1_rotate->Reset();
    h_TOF1_closest_rotate->Reset();

    h_int_TDC->Reset();
    h_int_TDC_selected->Reset();
    h_int_TDC_selected_weighted->Reset();
    h_int_TDC_Gap_Fibers->Reset();
    h_int_TDC_Gap_Fibers_SFT->Reset();
    h_int_TDC_Gap_Fibers_rotate->Reset();
    h_int_TDC_Gap_Fibers_SFT_rotate->Reset();
    h_int_TDC_SFT->Reset();
    h_int_TDC_TARGET->Reset();

    h_int_GoodLG->Reset();
    h_int_GoodLG_weighted->Reset();
    h_int_GoodLG_SFT->Reset();
    
    h_TDC_Gap_Fibers->Reset();
    h_K_Stop_Bars->Reset();
    h_TDC_Gap_Fibers_copy->Reset();
    h_TDC_Gap_Fibers_kaon->Reset();
    h_TARGET_LG_Blue->Reset();
    h_Centroid->Reset();


  } // ENDLOOP OVER EVENTS! (BATCH LOOP)
  
  cout << endl;
  cout << endl;

  fout.close();
  fout2.close();

  if(!gROOT->IsBatch()) return;

  delete h_Circle;
  delete h_target_ADC;
  delete h_target_TDC;
  delete h_target_TDC_copy;
  delete h_TDC_selected;
  delete h_TDC_selected2;
  
  delete h_target_ADC2;
  delete h_target_TDC2;

  delete h_target_ADC3;
  delete h_target_ADC4;
  delete h_target_ADCA;

  delete h_Target;
  delete h_Target_Center;

  delete h_kaon;
  delete h_kaon_copy;
    
  delete h_max;
  delete h_max_copy;

  delete h_GoodLG;
  delete h_GoodLG_weighted;
  delete h_GoodLG_copy;

  delete h_TOF1;
  delete h_TOF1_closest;
  delete h_TOF1_rotate; // ROTATE_CHANGE
  delete h_TOF1_closest_rotate; //ROTATE_CHANGE

  delete h_int_TDC;
  delete h_int_TDC_selected;
  delete h_int_TDC_selected_weighted;
  delete h_int_TDC_Gap_Fibers;
  delete h_int_TDC_Gap_Fibers_SFT;
  delete h_int_TDC_Gap_Fibers_rotate;  // ROTATE_CHANGE
  delete h_int_TDC_Gap_Fibers_SFT_rotate;  // ROTATE_CHANGE
  delete h_int_TDC_SFT;
  delete h_int_TDC_TARGET;

  delete h_int_GoodLG;
  delete h_int_GoodLG_weighted;
  delete h_int_GoodLG_SFT;
  
  delete h_TDC_Gap_Fibers;
  delete h_K_Stop_Bars;
  delete h_TDC_Gap_Fibers_copy; //ROTATE_CHANGE
  delete h_TDC_Gap_Fibers_kaon;
  delete h_TARGET_LG_Blue;
  delete h_Centroid;

  delete ell;
  delete ell_Target;
  delete ell_L1;

  delete HorizontalAxis;


} // End void


/*  
vector<double> GetZ(int Run_Number, int evt, double phi, int gap_to_fit, bool to_print, double C2X_centroid, double length_in_target){
  Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];        
  Int_t adc_low_sft[128];           Int_t ADC_Low_SFT[128];   
  Int_t tdc_le_sft[128][16];        Int_t TDC_LE_SFT[128];         
  Int_t tdc_te_sft[128][16];        Int_t TDC_TE_SFT[128];  

  Int_t HG_SFT_ADC_Thr[128] = {0};
  Int_t LG_SFT_ADC_Thr[128] = {0};

  
  for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;
  for(int i=0; i<128; i++)  LG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_LG[i]) + SFT_ADC_Thr_LG_Offset;

  Double_t ADC_High_SFT_corr[128]; 

  Int_t has_TDC_SFT_hit[128] = {0};


  //char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!

  //char file_mapping[200];
  //sprintf(file_mapping,"../Mapping");

  //char par_finput[200];
  //sprintf(par_finput,"%s/%s",file_mapping,source_mapping);

  //Int_t par_temp[2][128];
  //ifstream fdat(par_finput,ios::in);
  //for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
  //fdat.close();

  char path_input[200];                   
  sprintf(path_input,"%s",path_merged);          


  char Name_finput[200];
  sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);


  TChain *fChain= new TChain("Tree");   
  fChain->Add(Name_finput);   
  fChain->SetMakeClass(1);              

  fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
  fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
  fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
  fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);


  fChain->GetEntry(evt);  

  for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
    ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-HG_SFT_ADC_Thr[j_SFT];
    TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
    ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-LG_SFT_ADC_Thr[j_SFT];
  }

  for(int j=0 ; j<128 ; j++){ 
    if(ADC_High_SFT[j]<0)     
      ADC_High_SFT_corr[j]=0; 
    else
      ADC_High_SFT_corr[j]=ADC_High_SFT[j]; 
  }  

  for(Int_t ii=0; ii<128; ii++){
    for (Int_t qq=0; qq<6; qq++) {
      if (tdc_le_sft[ii][qq] > SFT_TDC_min[ii] && tdc_le_sft[ii][qq] < SFT_TDC_max[ii]) has_TDC_SFT_hit[ii]++;
    }
  }


  //double sft_z_selected = 0.;
  vector<double> ZZ;

  //sft_z_selected = SFT_print(ADC_High_SFT_corr, has_TDC_SFT_hit, SFT_channel_to_fiber, evt, phi, to_print, 0, false, C2X_centroid, length_in_target);
  ZZ = Z_list(ADC_High_SFT_corr, has_TDC_SFT_hit, SFT_channel_to_fiber, phi, true, 0, false, C2X_centroid, length_in_target);


  //return sft_z_selected;
  return ZZ;
}


vector<double> GetDeltaZ(int Run_Number, int evt, double phi, int gap_to_fit, bool to_print, double C2X_centroid, double length_in_target){
  Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];        
  Int_t adc_low_sft[128];           Int_t ADC_Low_SFT[128];   
  Int_t tdc_le_sft[128][16];        Int_t TDC_LE_SFT[128];         
  Int_t tdc_te_sft[128][16];        Int_t TDC_TE_SFT[128];  

  Int_t HG_SFT_ADC_Thr[128] = {0};
  Int_t LG_SFT_ADC_Thr[128] = {0};

  
  for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;
  for(int i=0; i<128; i++)  LG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_LG[i]) + SFT_ADC_Thr_LG_Offset;

  Double_t ADC_High_SFT_corr[128]; 

  Int_t has_TDC_SFT_hit[128] = {0};


  char path_input[200];                   
  sprintf(path_input,"%s",path_merged);          


  char Name_finput[200];
  sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);


  TChain *fChain= new TChain("Tree");   
  fChain->Add(Name_finput);   
  fChain->SetMakeClass(1);              

  fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
  fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
  fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
  fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);


  fChain->GetEntry(evt);  

  for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
    ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-HG_SFT_ADC_Thr[j_SFT];
    TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
    ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-LG_SFT_ADC_Thr[j_SFT];
  }

  for(int j=0 ; j<128 ; j++){ 
    if(ADC_High_SFT[j]<0)     
      ADC_High_SFT_corr[j]=0; 
    else
      ADC_High_SFT_corr[j]=ADC_High_SFT[j]; 
  }  

  for(Int_t ii=0; ii<128; ii++){
    for (Int_t qq=0; qq<6; qq++) {
      if (tdc_le_sft[ii][qq] > SFT_TDC_min[ii] && tdc_le_sft[ii][qq] < SFT_TDC_max[ii]) has_TDC_SFT_hit[ii]++;
    }
  }


  //double sft_z_selected = 0.;
  vector<double> DZ;

  //sft_z_selected = SFT_print(ADC_High_SFT_corr, has_TDC_SFT_hit, SFT_channel_to_fiber, evt, phi, to_print, 0, false, C2X_centroid, length_in_target);
  DZ = Delta_Z_list(ADC_High_SFT_corr, has_TDC_SFT_hit, SFT_channel_to_fiber, phi, true, 0, false, C2X_centroid, length_in_target);


  //return sft_z_selected;
  return DZ;
}

vector<double> GetTrackLength(int Run_Number, int evt, double phi, int gap_to_fit, bool to_print, double C2X_centroid, double length_in_target){
  Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];        
  Int_t adc_low_sft[128];           Int_t ADC_Low_SFT[128];   
  Int_t tdc_le_sft[128][16];        Int_t TDC_LE_SFT[128];         
  Int_t tdc_te_sft[128][16];        Int_t TDC_TE_SFT[128];  

  Int_t HG_SFT_ADC_Thr[128] = {0};
  Int_t LG_SFT_ADC_Thr[128] = {0};

  
  for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;
  for(int i=0; i<128; i++)  LG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_LG[i]) + SFT_ADC_Thr_LG_Offset;

  Double_t ADC_High_SFT_corr[128]; 

  Int_t has_TDC_SFT_hit[128] = {0};


  char path_input[200];                   
  sprintf(path_input,"%s",path_merged);          


  char Name_finput[200];
  sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);


  TChain *fChain= new TChain("Tree");   
  fChain->Add(Name_finput);   
  fChain->SetMakeClass(1);              

  fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
  fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
  fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
  fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);


  fChain->GetEntry(evt);  

  for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
    ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-HG_SFT_ADC_Thr[j_SFT];
    TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
    ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-LG_SFT_ADC_Thr[j_SFT];
  }

  for(int j=0 ; j<128 ; j++){ 
    if(ADC_High_SFT[j]<0)     
      ADC_High_SFT_corr[j]=0; 
    else
      ADC_High_SFT_corr[j]=ADC_High_SFT[j]; 
  }  

  for(Int_t ii=0; ii<128; ii++){
    for (Int_t qq=0; qq<6; qq++) {
      if (tdc_le_sft[ii][qq] > SFT_TDC_min[ii] && tdc_le_sft[ii][qq] < SFT_TDC_max[ii]) has_TDC_SFT_hit[ii]++;
    }
  }


  //double sft_z_selected = 0.;
  vector<double> Track_Length;

  //sft_z_selected = SFT_print(ADC_High_SFT_corr, has_TDC_SFT_hit, SFT_channel_to_fiber, evt, phi, to_print, 0, false, C2X_centroid, length_in_target);
  Track_Length = Track_Length_list(ADC_High_SFT_corr, has_TDC_SFT_hit, SFT_channel_to_fiber, phi, true, 0, false, C2X_centroid, length_in_target);


  //return sft_z_selected;
  return Track_Length;
}
*/
