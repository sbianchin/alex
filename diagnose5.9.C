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
#include <chrono>

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
#include "ADC_TARGET_Pedestals.h"
#include "TDC_Windows.h"
#include "Cuts_and_Windows.h"
#include "MWPC_Thr2.h"
#include "Pedestals.h"
#include "intersect.cxx"
#include "C2_Strip_transform.h"
#include "Channel_to_Strip.h"
#include "SFT_functions.h"
// #include "Plot_Event_Display.C"
#endif

using namespace std;

void grab_event(Int_t Run_Number=5, Int_t ievt=0){
  char path_input[200];
  sprintf(path_input,"%s",path_merged);

  char Name_finput[200];
  sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);
  cout << "File opened:  " << Name_finput << endl;
  int tdc_trigger[2][16];
  Int_t adc_high_target[256];       Int_t ADC_High_TARGET[256];
  Int_t adc_low_target[256];        Int_t ADC_Low_TARGET[256];
  Int_t tdc_le_target[256][16];     Int_t TDC_LE_TARGET[256];
  Int_t tdc_te_target[256][16];
  Int_t TDC_LE_TARGET_corrected[256][6]={0};
  char TDC_LE_TARGET_corr[256][6][20];

  Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];        Double_t ADC_High_SFT_corr[128];
  Int_t adc_low_sft[128];           Int_t ADC_Low_SFT[128];
  Int_t tdc_le_sft[128][16];        Int_t TDC_LE_SFT[128];
  Int_t tdc_te_sft[128][16];

  Int_t ADC_TOF1[24];

  Int_t ADC_tof1U[12];              Int_t ADC_TOF1U[12];
  Int_t ADC_tof1D[12];              Int_t ADC_TOF1D[12];

  Int_t TDC_tof1U[12];              Int_t TDC_TOF1U[12];
  Int_t TDC_tof1D[12];              Int_t TDC_TOF1D[12];

  Int_t ADC_tof2AO[12];             Int_t ADC_TOF2AO[12];
  Int_t ADC_tof2BO[12];             Int_t ADC_TOF2BO[12];
  Int_t ADC_tof2AI[12];             Int_t ADC_TOF2AI[12];
  Int_t ADC_tof2BI[12];             Int_t ADC_TOF2BI[12];

  Int_t TDC_tof2AO[12];             Int_t TDC_TOF2AO[12];
  Int_t TDC_tof2BO[12];             Int_t TDC_TOF2BO[12];
  Int_t TDC_tof2AI[12];             Int_t TDC_TOF2AI[12];
  Int_t TDC_tof2BI[12];             Int_t TDC_TOF2BI[12];

  Int_t MwpcADC[512];

  Int_t adc_c2x_r[56];              double ADC_C2X_R[56];
  Int_t adc_c2x_l[56];              double ADC_C2X_L[56];
  Int_t adc_c2y_r[16];              double ADC_C2Y_R[16];
  Int_t adc_c2y_l[16];              double ADC_C2Y_L[16];

  Int_t adc_c3x_r[64];              double ADC_C3X_R[64];
  Int_t adc_c3x_l[64];              double ADC_C3X_L[64];
  Int_t adc_c3y_r[16];              double ADC_C3Y_R[16];
  Int_t adc_c3y_l[16];              double ADC_C3Y_L[16];

  Int_t adc_c4x_r[72];              double ADC_C4X_R[72];
  Int_t adc_c4x_l[72];              double ADC_C4X_L[72];
  Int_t adc_c4y_r[16];              double ADC_C4Y_R[16];
  Int_t adc_c4y_l[16];              double ADC_C4Y_L[16];

  Int_t tdc_Ck[14][16];
  Int_t tdc_Cpi[14][16];
  Bool_t Event_flag[40] = {false};



  TChain *fChain= new TChain("Tree");
  fChain->Add(Name_finput);
  fChain->SetMakeClass(1);
  int tdc_vt48[256][16];
  fChain->SetBranchAddress("VT48_TDC",tdc_vt48);
  fChain->SetBranchAddress("TDC_Trig",tdc_trigger); // tdc_trigger for TDC_diff calculation

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
  for (int ivt = 4750; ivt < 100000; ++ivt) {
    auto t1 = std::chrono::high_resolution_clock::now();
    fChain->GetEntry(ivt);
    auto t2 = std::chrono::high_resolution_clock::now();
    cout << "Time taken "
                << std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count()
                << " microseconds\n";
  }

  }

// void diagnose5_9(Int_t Run_Number=5) {
//   int i = 4750;
//
//   while (i < 20000) {
//     // if (i % 400 == 0) {
//     //   t1 = std::chrono::high_resolution_clock::now();
//     // }
//     // if (i % 400 == 399) {
//     //   t2 = std::chrono::high_resolution_clock::now();
//     //   cout << "Time taken for 100 events "
//     //               << std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count()
//     //               << " microseconds\n";
//     // }
//     cout << i << endl;
//     // Event_Display_5_8(Run_Number, i, 0, 0, 0, 0);
//     grab_event(Run_Number, i);
//     // cout << dist_to_k_stop;
//     ++i;
//   }
//
// }
