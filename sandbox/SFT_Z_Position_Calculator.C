#ifndef __CINT__
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <TROOT.h>
#include "TSystem.h"
#include "TFile.h"
#include "TProfile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TRandom.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TPolyLine.h"
#include "TLine.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TBranch.h"
#include "TLegend.h"
#include "TLatex.h"
//#include "TSpectrum.h"
#include "TMarker.h"
//#include "Event_Display_MS.h"
#include "ANAPATH.h"
#include "Thresholds.h"
#endif
  
void SFT_Z_Position_Calculator(Int_t Run_Number=5, Int_t ievt=0, float phi = 45){ 

gStyle->SetOptStat(0);

char* source_mapping="SFT_Mapping_Oct14.txt";  // Mapping file !!!

Int_t adc_high_target[256];   		Int_t ADC_High_TARGET[256];    	Double_t ADC_High_TARGET_corr[256];   	
Int_t adc_low_target[256]; 		Int_t ADC_Low_TARGET[256]; 	//Double_t ADC_Low_TARGET_corr[256];
Int_t tdc_le_target[256][16];      	Int_t TDC_LE_TARGET[256];    	Double_t TDC_LE_TARGET_corr[256];   	
Int_t tdc_te_target[256][16]; 		Int_t TDC_TE_TARGET[256]; 	//Double_t TDC_TE_TARGET_corr[256];

Int_t adc_high_sft[128];   		Int_t ADC_High_SFT[128];    	Double_t fiber[128];   	
Int_t adc_low_sft[128]; 		Int_t ADC_Low_SFT[128]; 	Double_t fiberTDC[128];
Int_t tdc_le_sft[128][16];      	Int_t TDC_LE_SFT[128];    	//Double_t TDC_LE_SFT_corr[128];   	
Int_t tdc_te_sft[128][16]; 		Int_t TDC_TE_SFT[128]; 		//Double_t TDC_TE_SFT_corr[128];	

Int_t ADC_tof1[24];	Int_t ADC_TOF1[24];	
Int_t ADC_tof2[56];	Int_t ADC_TOF2[56];

Int_t TDC_tof1U[12];  Int_t TDC_TOF1U[12];
Int_t TDC_tof1D[12];  Int_t TDC_TOF1D[12];

Int_t TDC_tof2AO[12];   Int_t TDC_TOF2AO[12];
Int_t TDC_tof2BO[12];   Int_t TDC_TOF2BO[12];
Int_t TDC_tof2AI[12];   Int_t TDC_TOF2AI[12];
Int_t TDC_tof2BI[12];   Int_t TDC_TOF2BI[12];

Int_t MwpcADC[512];	Int_t MWPCADC[512];

 /////////////////////////////////////////////////////////////////////////////////////////////////////


char path_input[200];                   char file_mapping[200];
sprintf(path_input,path_merged);          sprintf(file_mapping,"../Mapping");
//sprintf(path_input,"/triumfcs/trshare/trek/E36/Data/April_2015/root");

char Name_finput[200];
sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);

////

char par_finput[200];
//sprintf(par_finput,"/triumfcs/trshare/trek/E36/Macros/April_2015/Mapping/%s",source_mapping);
sprintf(par_finput,"%s/%s",file_mapping,source_mapping);

char par_finput2[200];
//sprintf(par_finput,"/triumfcs/trshare/trek/E36/Macros/April_2015/Mapping/%s",source_mapping);
sprintf(par_finput2,"%s/MWPC_map.txt",file_mapping);

char par_finput3[200];
//sprintf(par_finput,"/triumfcs/trshare/trek/E36/Macros/April_2015/Mapping/%s",source_mapping);
sprintf(par_finput3,"%s/ADC_TARGET_Thresholds.txt",file_mapping);

char par_finput4[200];
//sprintf(par_finput,"/triumfcs/trshare/trek/E36/Macros/April_2015/Mapping/%s",source_mapping);
sprintf(par_finput4,"%s/ADC_SFT_Thresholds.txt",file_mapping);

char par_finput5[200];
//sprintf(par_finput,"/triumfcs/trshare/trek/E36/Macros/April_2015/Mapping/%s",source_mapping);
sprintf(par_finput5,"%s/ADC_TOF1_Thresholds.txt",file_mapping);

char par_finput6[200];
//sprintf(par_finput,"/triumfcs/trshare/trek/E36/Macros/April_2015/Mapping/%s",source_mapping);
sprintf(par_finput6,"%s/ADC_TOF2_Thresholds.txt",file_mapping);

////

cout << "   " << endl;
cout << Name_finput << endl;

cout << "   " << endl;
cout << "SFT Mapping File:   " << par_finput << endl;
cout << "" << endl;
cout << "MWPC Mapping File:   " << par_finput2 << endl;
cout << "" << endl;

////

Int_t par_temp[2][128];
ifstream fdat(par_finput,ios::in);
for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
fdat.close();

char par_temp2[512][50];
//char par_temp3[512][50];
ifstream fdat2(par_finput2,ios::in);
for(Int_t ii=0; ii<512; ii++) fdat2 >> par_temp2[ii];
fdat2.close();

Int_t par_temp_TARGET[2][256];
ifstream fdat3(par_finput3,ios::in);
for(Int_t ii=0; ii<256; ii++) fdat3 >> par_temp_TARGET[0][ii] >> par_temp_TARGET[1][ii];
fdat3.close();

Int_t par_temp_SFT[2][128];
ifstream fdat4(par_finput4,ios::in);
for(Int_t ii=0; ii<128; ii++) fdat4 >> par_temp_SFT[0][ii] >> par_temp_SFT[1][ii];
fdat4.close();

Int_t par_temp_TOF1[2][24];
ifstream fdat5(par_finput5,ios::in);
for(Int_t ii=0; ii<24; ii++) fdat5 >> par_temp_TOF1[0][ii] >> par_temp_TOF1[1][ii];
fdat5.close();

Int_t par_temp_TOF2[2][56];
ifstream fdat6(par_finput6,ios::in);
for(Int_t ii=0; ii<56; ii++) fdat6 >> par_temp_TOF2[0][ii] >> par_temp_TOF2[1][ii];
fdat6.close();

TChain *fChain= new TChain("Tree");		
fChain->Add(Name_finput);		
fChain->SetMakeClass(1);							

fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);		fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
fChain->SetBranchAddress("ADC_Low_TARGET",adc_low_target);		fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);		fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
fChain->SetBranchAddress("TDC_TE_TARGET",tdc_te_target);		fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

fChain->SetBranchAddress("ADC_TOF1",ADC_TOF1);
fChain->SetBranchAddress("ADC_TOF2",ADC_TOF2);

fChain->SetBranchAddress("TDC_TOF1U",TDC_tof1U);
fChain->SetBranchAddress("TDC_TOF1D",TDC_tof1D);
fChain->SetBranchAddress("TDC_TOF2AO",TDC_tof2AO);
fChain->SetBranchAddress("TDC_TOF2AI",TDC_tof2AI);
fChain->SetBranchAddress("TDC_TOF2BO",TDC_tof2BO);
fChain->SetBranchAddress("TDC_TOF2BI",TDC_tof2BI);    

fChain->SetBranchAddress("MWPCADC",MwpcADC);

//Int_t nentries;
Int_t nentries = (Int_t)fChain->GetEntries();//		Int_t nentries_SFT = (Int_t)fChain->GetEntries();
//if(nentries_TARGET==nentries_SFT) nentries = nentries_TARGET;
cout << "  " << endl;
cout << "****  Number of events: " << nentries << "  **** " <<endl;
cout << "  " << endl;

for(Int_t i=ievt; i<ievt+1; i++){
	fChain->GetEntry(i);	//	fChain_SFT->GetEntry(i);

  int has_data_TOF = 0;
  int has_data_TARGET = 0;
  int has_data_low_TARGET = 0;

  int C2X_hit = 0;
  int C2Y_hit = 0;
  int C3X_hit = 0;
  int C3Y_hit = 0;
  int C4X_hit = 0;
  int C4Y_hit = 0;

  int C2X_hitA = 0;
  int C2Y_hitA = 0;
  int C3X_hitA = 0;
  int C3Y_hitA = 0;
  int C4X_hitA = 0;
  int C4Y_hitA = 0;

  for (int j_TARGET=0; j_TARGET<256; j_TARGET++){
    ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-ADC_cut_TARGET;
    ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-ADC_cut_TARGET;
    TDC_LE_TARGET[j_TARGET]=tdc_le_target[j_TARGET][0];
    TDC_TE_TARGET[j_TARGET]=tdc_te_target[j_TARGET][0];   
//      cout << j_TARGET << "   " << ADC_High_TARGET[j_TARGET] << "   " << TDC_LE_TARGET[j_TARGET] << endl;
  } 
  	
 	for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
		ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-ADC_cut_SFT;
		ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-ADC_cut_SFT;
		TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
		TDC_TE_SFT[j_SFT]=tdc_te_sft[j_SFT][0];
	
	//	cout << j_SFT << "   " << ADC_High_SFT[j_SFT] << "   " << fiberTDC[j_SFT] << endl;
	}

  for(int j=0 ; j<128 ; j++){	
  	if(ADC_High_SFT[j]<0)   	 fiber[j]=0; 
  	if(ADC_High_SFT[j]>=0)  	 fiber[j]=ADC_High_SFT[j]; 
  }

  for (Int_t j_TOF1=0; j_TOF1<24; j_TOF1++) {
//    cout << j_TOF1 << " -- " << ADC_TOF1[j_TOF1] << endl;
    ADC_TOF1[j_TOF1] = ADC_TOF1[j_TOF1]-par_temp_TOF1[1][j_TOF1];
  }

  for (Int_t j_TOF2=0; j_TOF2<56; j_TOF2++) {
    ADC_TOF2[j_TOF2] = ADC_TOF2[j_TOF2]-par_temp_TOF2[1][j_TOF2];
  }

  for (Int_t j_MWPC=0; j_MWPC<512; j_MWPC++) {
    MWPCADC[j_MWPC] = MwpcADC[j_MWPC]-MWPC_thr;
  }

  for (Int_t j_TDCTOF=0; j_TDCTOF<12; j_TDCTOF++) {
      TDC_TOF1U[j_TDCTOF] = TDC_tof1U[j_TDCTOF];
      TDC_TOF1D[j_TDCTOF] = TDC_tof1D[j_TDCTOF];
      TDC_TOF2AO[j_TDCTOF] = TDC_tof2AO[j_TDCTOF];
      TDC_TOF2BO[j_TDCTOF] = TDC_tof2BO[j_TDCTOF];
      TDC_TOF2AI[j_TDCTOF] = TDC_tof2AI[j_TDCTOF];
      TDC_TOF2BI[j_TDCTOF] = TDC_tof2BI[j_TDCTOF];
    }

    //if ((TDC_TOF1U[0]>TDC_TOF1_min && TDC_TOF1U[0]<TDC_TOF1_max) || (TDC_TOF1D[0]>TDC_TOF1_min && TDC_TOF1D[0]<TDC_TOF1_max))

  ///

  for (Int_t p=0; p<256; p++) {
    if (ADC_High_TARGET[p] > 0) has_data_TARGET++;
    if (ADC_Low_TARGET[p] > 0) has_data_low_TARGET++;
  }

  bool has_ADC_TOF1_hit[12] = {false};
  bool has_TDC_TOF1_hit[12] = {false};
  bool has_ADC_TOF2_hit[12] = {false};
  bool has_TDC_TOF2_hit[12] = {false};

    for (int kk=0; kk<12; kk++) {
      if (kk == 6) {
        if ((ADC_TOF2[6]>0) || (ADC_TOF2[30]>0) || (ADC_TOF2[55]>0) || (ADC_TOF2[42]>0)) {has_ADC_TOF2_hit[6]=true;}
        if (((TDC_TOF2AO[6]>TDC_TOF2_min && TDC_TOF2AO[6] < TDC_TOF2_max) || (TDC_TOF2AI[6]>TDC_TOF2_min && TDC_TOF2AI[6] < TDC_TOF2_max)) 
          || ((TDC_TOF2BO[6]>TDC_TOF2_min && TDC_TOF2BO[6] < TDC_TOF2_max) || (TDC_TOF2BI[6]>TDC_TOF2_min && TDC_TOF2BI[6] < TDC_TOF2_max))) {has_TDC_TOF2_hit[6]=true;}

        if (ADC_TOF1[6]>0 || ADC_TOF1[18]>0) {has_ADC_TOF1_hit[6] = true;}
          if ((TDC_TOF1U[6]>TDC_TOF1_min && TDC_TOF1U[6]<TDC_TOF1_max) || (TDC_TOF1D[6]>TDC_TOF1_min && TDC_TOF1D[6]<TDC_TOF1_max)) {has_TDC_TOF1_hit[6] = true;}
        }

      else {
        if ((ADC_TOF2[kk]>0) || (ADC_TOF2[kk+24]>0) || (ADC_TOF2[kk+12]>0) || (ADC_TOF2[kk+36]>0)) {has_ADC_TOF2_hit[kk]=true;}
        if (((TDC_TOF2AO[kk]>TDC_TOF2_min && TDC_TOF2AO[kk] < TDC_TOF2_max) || (TDC_TOF2AI[kk]>TDC_TOF2_min && TDC_TOF2AI[kk] < TDC_TOF2_max)) 
          || ((TDC_TOF2BO[kk]>TDC_TOF2_min && TDC_TOF2BO[kk] < TDC_TOF2_max) || (TDC_TOF2BI[kk]>TDC_TOF2_min && TDC_TOF2BI[kk] < TDC_TOF2_max))) {has_TDC_TOF2_hit[kk]=true;}

        if (ADC_TOF1[kk]>0 || ADC_TOF1[kk+12]>0) {has_ADC_TOF1_hit[kk] = true;}
          if ((TDC_TOF1U[kk]>TDC_TOF1_min && TDC_TOF1U[kk]<TDC_TOF1_max) || (TDC_TOF1D[kk]>TDC_TOF1_min && TDC_TOF1D[kk]<TDC_TOF1_max)) {has_TDC_TOF1_hit[kk] = true;}
      }
    }

    for (int kk=0; kk<12; kk++) {
      if (kk == 0) {
        if (has_ADC_TOF2_hit[0] && has_TDC_TOF2_hit[0]) {
          if ((has_ADC_TOF1_hit[0] && has_TDC_TOF1_hit[0]) || (has_ADC_TOF1_hit[11] && has_TDC_TOF1_hit[11]) || (has_ADC_TOF1_hit[1] && has_TDC_TOF1_hit[1])) {
            has_data_TOF++;
          }
        }
      }
      else if (kk == 11) {
        if (has_ADC_TOF2_hit[11] && has_TDC_TOF2_hit[11]) {
          if ((has_ADC_TOF1_hit[11] && has_TDC_TOF1_hit[11]) || (has_ADC_TOF1_hit[10] && has_TDC_TOF1_hit[10]) || (has_ADC_TOF1_hit[0] && has_TDC_TOF1_hit[0])) {
            has_data_TOF++;
          }
        }
      }
      else {
        if (has_ADC_TOF2_hit[kk] && has_TDC_TOF2_hit[kk]) {
          if ((has_ADC_TOF1_hit[kk] && has_TDC_TOF1_hit[kk]) || (has_ADC_TOF1_hit[kk-1] && has_TDC_TOF1_hit[kk-1]) || (has_ADC_TOF1_hit[kk+1] && has_TDC_TOF1_hit[kk+1])) {
            has_data_TOF++;
          }
        }
      }
    }

  ///
  for (int q=0; q<512; q++) {

    //C3 Counters
    if (q >= 0 && q <= 95) {
      if (MWPCADC[q]>0) {
        C3X_hit++;
        if (q==0) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C3X_hitA++;}
        else if (q==95) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C3X_hitA++;}
        else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C3X_hitA++;}
      }
    }

    //C2 Counters
    if (q >= 96 && q <= 127) {
      if (MWPCADC[q]>0) {
        C2Y_hit++;
        if (q==96) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C2Y_hitA++;}
        else if (q==127) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C2Y_hitA++;}
        else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C2Y_hitA++;}
      }
    }

    if (q >= 128 && q <= 141) {
      if (MWPCADC[q]>0) {
        C2X_hit++;
        if (q==128) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C2X_hitA++;}
        else if (q==141) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C2X_hitA++;}
        else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C2X_hitA++;}
      }
    }
    if (q >= 144 && q <= 157) {
      if (MWPCADC[q]>0) {
        C2X_hit++;
        if (q==144) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C2X_hitA++;}
        else if (q==157) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C2X_hitA++;}
        else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C2X_hitA++;}
      }
    }
    if (q >= 160 && q <= 173) {
      if (MWPCADC[q]>0) {
        C2X_hit++;
        if (q==160) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C2X_hitA++;}
        else if (q==173) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C2X_hitA++;}
        else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C2X_hitA++;}
      }
    }
    if (q >= 176 && q <= 189) {
      if (MWPCADC[q]>0) {
        C2X_hit++;
        if (q==176) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C2X_hitA++;}
        else if (q==189) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C2X_hitA++;}
        else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C2X_hitA++;}
      }
    }
    if (q >= 192 && q <= 205) {
      if (MWPCADC[q]>0) {
        C2X_hit++;
        if (q==192) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C2X_hitA++;}
        else if (q==205) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C2X_hitA++;}
        else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C2X_hitA++;}
      }
    }
    if (q >= 208 && q <= 221) {
      if (MWPCADC[q]>0) {
        C2X_hit++;
        if (q==208) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C2X_hitA++;}
        else if (q==221) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C2X_hitA++;}
        else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C2X_hitA++;}
      }
    }
    if (q >= 224 && q <= 237) {
      if (MWPCADC[q]>0) {
        C2X_hit++;
        if (q==224) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C2X_hitA++;}
        else if (q==237) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C2X_hitA++;}
        else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C2X_hitA++;}
      }
    }
    if (q >= 240 && q <= 253) {
      if (MWPCADC[q]>0) {
        C2X_hit++;
        if (q==240) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C2X_hitA++;}
        else if (q==253) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C2X_hitA++;}
        else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C2X_hitA++;}
      }
    }

    //C4 Counters
    if (q >= 256 && q <= 287) {
      if (MWPCADC[q]>0) {
        C4Y_hit++;
        if (q==256) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C4Y_hitA++;}
        else if (q==287) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C4Y_hitA++;}
        else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C4Y_hitA++;}
      }
    }
    if (q >= 288 && q <= 295) {
      if (MWPCADC[q]>0) {
        C4X_hit++;
        if (q==288) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C4X_hitA++;}
        else if (q==295) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C4X_hitA++;}
        else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C4X_hitA++;}
      }
    }
    if (q >= 304 && q <= 311) {
      if (MWPCADC[q]>0) {
        C4X_hit++;
        if (q==304) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C4X_hitA++;}
        else if (q==311) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C4X_hitA++;}
        else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C4X_hitA++;}
      }
    }
    if (q >= 320 && q <= 447) {
      if (MWPCADC[q]>0) {
        C4X_hit++;
        if (q==320) {if (MWPCADC[q+1]>0) C4X_hitA++;}
        else if (q==447) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C4X_hitA++;}
        else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C4X_hitA++;}
      }
    }

    //C3 Counters
    if (q >= 448 && q <= 479) {
      if (MWPCADC[q]>0) {
        C3Y_hit++;
        if (q==448) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C3Y_hitA++;}
        else if (q==463) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C3Y_hitA++;}
        else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C3Y_hitA++;}
      }
    }
    if (q >= 480 && q <= 511) {
      if (MWPCADC[q]>0) {
        C3X_hit++;
        if (q==480) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C3X_hitA++;}
        else if (q==495) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C3X_hitA++;}
        else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C3X_hitA++;}
      }
    }
  }

  bool good_MWPC = false;
  bool good_TOF = false;
  bool good_TARGET = false;

  if ((C2X_hit >=1) && (C2Y_hit >=1) && (C3X_hit >=1) && 
      (C3Y_hit >=1) && (C4X_hit >=1) && (C4Y_hit >=1)) good_MWPC = true;

  if (has_data_TOF > 0) good_TOF = true;

  if ((has_data_TARGET >= 4) && (has_data_low_TARGET >= 3)) good_TARGET = true;

  for(Int_t ii=0; ii<128; ii++){
    fiber[ii]=ADC_High_SFT[par_temp[1][ii]];
    if(tdc_le_sft[par_temp[1][ii]][0]>TDC_min_SFT && tdc_le_sft[par_temp[1][ii]][0]<TDC_max_SFT)
      fiberTDC[ii]=tdc_le_sft[par_temp[1][ii]][0];
    else if(tdc_le_sft[par_temp[1][ii]][1]>TDC_min_SFT && tdc_le_sft[par_temp[1][ii]][1]<TDC_max_SFT)
      fiberTDC[ii]=tdc_le_sft[par_temp[1][ii]][1];
    else if(tdc_le_sft[par_temp[1][ii]][2]>TDC_min_SFT && tdc_le_sft[par_temp[1][ii]][2]<TDC_max_SFT)
      fiberTDC[ii]=tdc_le_sft[par_temp[1][ii]][2];
    else
      fiberTDC[ii]=tdc_le_sft[par_temp[1][ii]][3];
  }

    Int_t L1[15] = {0};
    Int_t L2[15] = {0};
    Int_t L3[17] = {0};
    Int_t L4[17] = {0};

  	for(Int_t imark1=0; imark1<15; imark1++){
  		if((fiber[imark1]!=0) && (fiberTDC[imark1]>TDC_min_SFT) && (fiberTDC[imark1]<TDC_max_SFT)){  
      L1[imark1]++;
  		}

  		if((fiber[imark1+64]!=0) && (fiberTDC[imark1+64]>TDC_min_SFT) && (fiberTDC[imark1+64]<TDC_max_SFT)){  
      L1[imark1]++; 
  		}
  	}

  	for(Int_t imark2=0; imark2<15; imark2++){
  		if((fiber[imark2+15]!=0) && (fiberTDC[imark2+15]>TDC_min_SFT) && (fiberTDC[imark2+15]<TDC_max_SFT)){  
      L2[imark2]++;
  		}

  		if((fiber[imark2+79]!=0) && (fiberTDC[imark2+79]>TDC_min_SFT) && (fiberTDC[imark2+79]<TDC_max_SFT)){  
      L2[imark2]++;
  		}
  	}

  	for(Int_t imark3=0; imark3<17; imark3++){
  		if((fiber[imark3+30]!=0) && (fiberTDC[imark3+30]>TDC_min_SFT) && (fiberTDC[imark3+30]<TDC_max_SFT)){  
      L3[imark3]++;
  		}

  		if((fiber[imark3+94]!=0) && (fiberTDC[imark3+94]>TDC_min_SFT) && (fiberTDC[imark3+94]<TDC_max_SFT)){  
      L3[imark3]++;
  		}
  	}


  	for(Int_t imark4=0; imark4<17; imark4++){
  		if((fiber[imark4+47]!=0) && (fiberTDC[imark4+47]>TDC_min_SFT) && (fiberTDC[imark4+47]<TDC_max_SFT)){  
      L4[imark4]++;
  		}
  	
  		if((fiber[imark4+111]!=0) && (fiberTDC[imark4+111]>TDC_min_SFT) && (fiberTDC[imark4+111]<TDC_max_SFT)){  
      L4[imark4]++;
  		}
  	}

    /*
    cout << "L1" << endl;
    for (int q=0; q<15; q++) {
      cout << q+1 << " : " << L1[q] << endl;
    }

    cout << "" << endl;
    cout << "L2" << endl;
    for (int q=0; q<15; q++) {
      cout << q+1 << " : " << L2[q] << endl;
    }

    cout << "" << endl;
    cout << "L3" << endl;
    for (int q=0; q<17; q++) {
      cout << q+1 << " : " <<  L3[q] << endl;
    }

    cout << "" << endl;
    cout << "L4" << endl;
    for (int q=0; q<17; q++) {
      cout << q+1 << " : " <<  L4[q] << endl;
    }
    */

    if ((good_TOF) && (good_MWPC) && (good_TARGET)) {
      
      cout << "Good TARGET, Gap, and MWPC - Performing Calculation:" << endl;
      cout << "" << endl;

      for (int q=0; q<15; q++) {
        if (L1[q]>0) {
          for (int r=0; r<15; r++) {
            if (L2[r]>0) {
              for (int s=0; s<17; s++) {
                if (L3[s]>0) {
                  for (int t=0; t<17; t++) {
                    if (L4[t]>0) {

                      float thr = 2.;
                      float pi = 3.14195;

                      //Constants
                      float z1 = 141.;
                      float z2 = 135.;
                      float z3 = 157.2;
                      float z4 = 168.;

                      float R1 = 40.25;
                      float R2 = 41.1;
                      float R3 = 42.1;
                      float R4 = 42.95;

                      float alpha1 = 0.0623082017;
                      float alpha2 = 0.0610864722;
                      float alpha3 = 0.067544185;
                      float alpha4 = 0.0661479228;

                      float phi1 = 165.;
                      float phi2 = 15.;
                      float phi3 = 75.;
                      float phi4 = 285.;

                      float Zpos1[18] = {0};
                      float Zpos2[18] = {0};
                      float Zpos3[16] = {0};
                      float Zpos4[14] = {0};

                      for (int n=0; n<18; n++) {
                        Zpos1[n] = z1 + (phi - phi1)*(pi/180)*(R1*tan(alpha1)) + (n-1)*15.7 + (15-q+1)*1.05 + 0.5;
                        Zpos2[n] = z2 + (phi - phi2)*(pi/180)*(R2*tan(alpha2)) + (n-1)*15.7 + (15-r+1)*1.05 + 0.5;
                        if (n<16) Zpos3[n] = z3 + (360 - phi + phi3)*(pi/180)*(R3*tan(alpha3)) + (n-1)*17.8 + (s)*1.05 + 0.5;
                        if (n<14) Zpos4[n] = z4 + (360 - phi + phi4)*(pi/180)*(R4*tan(alpha4)) + (n-1)*17.8 + (t)*1.05 + 0.5;
                      }

                      //for (int n=0; n<18; n++) {
                      //  cout << "" << endl;
                      //  cout << "N = " << n << endl;
                      //  cout << "Z1 with n = " << n << ": " << Zpos1[n] << endl;
                      //  cout << "Z2 with n = " << n << ": " << Zpos2[n] << endl;
                      //  cout << "Z3 with n = " << n << ": " << Zpos3[n] << endl;
                      //  cout << "Z4 with n = " << n << ": " << Zpos4[n] << endl;
                      //}

                      for (int i=1; i<18; i++) {
                        for (int j=1; j<18; j++) {
                          if (fabs(Zpos1[i]-Zpos2[j]) <= thr) {
                            for (int k=1; k<16; k++) {
                              if (fabs(Zpos1[i]-Zpos3[k]) <= thr) {
                                for (int h=1; h<14; h++) {
                                  if (fabs(Zpos1[i]-Zpos4[h]) <= thr) {
                                    cout << "Z1 -- N = " << i << ": " << Zpos1[i] << " | Z2 -- N = " << j << ": " << Zpos2[j] << " | Z3 -- N = " << k << ": " 
                                    << Zpos3[k] << " | Z4 -- N = " << h << ": " << Zpos4[h] << endl;
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    else cout << "Bad gap, TARGET or MWPC - Calculation not perfomed" << endl;
  }

}