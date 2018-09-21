#ifndef __CINT__
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>
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
#include "MWPC_Thr.h"
#endif

#include "intersect.cxx"
#include "SFT_functions_test.h"

using namespace std;


// All lengths in millimeters

//                                          degrees
void SFT_Test(int Run_Number, int evt, double phi){

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


  char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!

  char file_mapping[200];
  sprintf(file_mapping,"../Mapping");

  char par_finput[200];
  sprintf(par_finput,"%s/%s",file_mapping,source_mapping);

  Int_t par_temp[2][128];
  ifstream fdat(par_finput,ios::in);
  for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
  fdat.close();

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

  double z = 0;

  z = SFT_print(ADC_High_SFT_corr, has_TDC_SFT_hit, SFT_channel_to_fiber, evt, phi, true, 0, false,0 ,0);


  return;
}
