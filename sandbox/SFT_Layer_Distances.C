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
  
void SFT_Layer_Distances(Int_t Run_Number=5, Int_t flag=0){ 

gStyle->SetOptStat(0);

char source_mapping[]="SFT_Mapping_Oct14.txt";  // Mapping file !!!

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

Int_t ADC_tof1U[12];  Int_t ADC_TOF1U[12];
Int_t ADC_tof1D[12];  Int_t ADC_TOF1D[12];

Int_t ADC_tof2AO[12];   Int_t ADC_TOF2AO[12];
Int_t ADC_tof2BO[12];   Int_t ADC_TOF2BO[12];
Int_t ADC_tof2AI[12];   Int_t ADC_TOF2AI[12];
Int_t ADC_tof2BI[12];   Int_t ADC_TOF2BI[12];

Int_t TDC_tof1U[12];  Int_t TDC_TOF1U[12];
Int_t TDC_tof1D[12];  Int_t TDC_TOF1D[12];

Int_t TDC_tof2AO[12];   Int_t TDC_TOF2AO[12];
Int_t TDC_tof2BO[12];   Int_t TDC_TOF2BO[12];
Int_t TDC_tof2AI[12];   Int_t TDC_TOF2AI[12];
Int_t TDC_tof2BI[12];   Int_t TDC_TOF2BI[12];

Int_t MwpcADC[512];	Int_t MWPCADC[512];

Int_t L1_L2count[15]={0};
Int_t L3_L4count[17]={0};

Int_t good_events = 0;

Int_t events[31] = {1,2,6,8,12,14,18,28,32,38,39,41,42,46,51,52,58,60,65,67,68,69,73,76,80,81,86,89,90,95,100};

 /////////////////////////////////////////////////////////////////////////////////////////////////////

char path_input[200];                   char file_mapping[200];
sprintf(path_input,"%s",path_merged);          sprintf(file_mapping,"../Mapping");
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

//fChain->SetBranchAddress("ADC_TOF1",ADC_tof1);
//fChain->SetBranchAddress("ADC_TOF2",ADC_tof2);
fChain->SetBranchAddress("ADC_TOF1U", ADC_tof1U);
fChain->SetBranchAddress("ADC_TOF1D", ADC_tof1D);

fChain->SetBranchAddress("ADC_TOF2AO", ADC_tof2AO);
fChain->SetBranchAddress("ADC_TOF2AI", ADC_tof2AI);
fChain->SetBranchAddress("ADC_TOF2BO", ADC_tof2BO);
fChain->SetBranchAddress("ADC_TOF2BI", ADC_tof2BI);

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


cout << "   " << endl;
if (flag!=0) nentries = flag;
for(Int_t i=0; i<nentries; i++){
	fChain->GetEntry(i);	//	fChain_SFT->GetEntry(i);

	if(i%10000==1)	cout<<"**** "<<i<<" events done"<<endl;

  int has_data_TARGET = 0;
  int has_data_low_TARGET = 0;
  int has_data_SFT = 0;
  int has_data_TOF = 0;
  int has_data_MWPC = 0;

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
  	
    //	cout << j_TARGET << "   " << ADC_High_TARGET[j_TARGET] << "   " << TDC_LE_TARGET[j_TARGET] << endl;
  }	


 	for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
		ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-ADC_cut_SFT;
		ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-ADC_cut_SFT;
		TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
		TDC_TE_SFT[j_SFT]=tdc_te_sft[j_SFT][0];
	
	//	cout << j_SFT << "   " << ADC_High_SFT[j_SFT] << "   " << fiberTDC[j_SFT] << endl;
	}

	//for (Int_t j_TOF1=0; j_TOF1<24; j_TOF1++) {
	//	ADC_TOF1[j_TOF1] = ADC_tof1[j_TOF1]-ADC_TOF1_thr;
	//}

	//for (Int_t j_TOF2=0; j_TOF2<56; j_TOF2++) {
	//	ADC_TOF2[j_TOF2] = ADC_tof2[j_TOF2]-ADC_TOF2_thr;
	//}

  for (Int_t j_TOF1=0; j_TOF1<12; j_TOF1++) {
    ADC_TOF1U[j_TOF1] = ADC_tof1U[j_TOF1]-par_temp_TOF1[1][j_TOF1];
    ADC_TOF1D[j_TOF1] = ADC_tof1D[j_TOF1]-par_temp_TOF1[1][j_TOF1+12];
  }

  for (Int_t j_TOF2=0; j_TOF2<12; j_TOF2++) {
    ADC_TOF2AO[j_TOF2] = ADC_tof2AO[j_TOF2]-par_temp_TOF2[1][j_TOF2];
    ADC_TOF2AI[j_TOF2] = ADC_tof2AI[j_TOF2]-par_temp_TOF2[1][j_TOF2+24];
    ADC_TOF2BO[j_TOF2] = ADC_tof2BO[j_TOF2]-par_temp_TOF2[1][j_TOF2+12];
    ADC_TOF2BI[j_TOF2] = ADC_tof2BI[j_TOF2]-par_temp_TOF2[1][j_TOF2+36];
  }


  for(int i=0; i<12; i++){
    ADC_TOF1[i] = ADC_TOF1U[i];
    ADC_TOF1[i+12] = ADC_TOF1D[i];
    ADC_TOF2[i] = ADC_TOF2AO[i];
    ADC_TOF2[i+12] = ADC_TOF2BO[i];
    ADC_TOF2[i+24] = ADC_TOF2AI[i];
    ADC_TOF2[i+36] = ADC_TOF2BI[i];
  }


	for (Int_t j_MWPC=0; j_MWPC<512; j_MWPC++) {
		MWPCADC[j_MWPC] = MwpcADC[j_MWPC];
	}

  for (Int_t j_TDCTOF=0; j_TDCTOF<12; j_TDCTOF++) {
    TDC_TOF1U[j_TDCTOF] = TDC_tof1U[j_TDCTOF];
    TDC_TOF1D[j_TDCTOF] = TDC_tof1D[j_TDCTOF];
    TDC_TOF2AO[j_TDCTOF] = TDC_tof2AO[j_TDCTOF];
    TDC_TOF2BO[j_TDCTOF] = TDC_tof2BO[j_TDCTOF];
    TDC_TOF2AI[j_TDCTOF] = TDC_tof2AI[j_TDCTOF];
    TDC_TOF2BI[j_TDCTOF] = TDC_tof2BI[j_TDCTOF];
  }

  for(int i=0 ; i<256 ; i++){	
  	if(ADC_High_TARGET[i]<0)     ADC_High_TARGET_corr[i]=0; 
  	if(ADC_High_TARGET[i]>=0)    ADC_High_TARGET_corr[i]=ADC_High_TARGET[i]; 
  	if(TDC_LE_TARGET[i]<0)  	 TDC_LE_TARGET_corr[i]=0; 
  	if(TDC_LE_TARGET[i]>=0)  	 TDC_LE_TARGET_corr[i]=TDC_LE_TARGET[i]; 	
  }

  for(int j=0 ; j<128 ; j++){	
  	if(ADC_High_SFT[j]<0)   	 fiber[j]=0; 
  	if(ADC_High_SFT[j]>=0)  	 fiber[j]=ADC_High_SFT[j]; 
  }

  for (Int_t p=0; p<256; p++) {
    if (ADC_High_TARGET[p] > 0) has_data_TARGET++;
    if (ADC_Low_TARGET[p] > 0) has_data_low_TARGET++;
  }

  for (int kk=0; kk<12; kk++) {
    if (kk == 0) {
      if ((ADC_TOF2[0]>0) || (ADC_TOF2[24]>0) || (ADC_TOF2[12]>0) || (ADC_TOF2[36]>0)) {
        if((ADC_TOF1[0] > 0) || (ADC_TOF1[12] > 0) || (ADC_TOF1[1] > 0) || (ADC_TOF1[13] > 0) || (ADC_TOF1[11] > 0) || (ADC_TOF1[23] > 0)) has_data_TOF++;
      } 
    }

    else if (kk == 6) {
      if ((ADC_TOF2[6]>0) || (ADC_TOF2[30]>0) || (ADC_TOF2[55]>0) || (ADC_TOF2[42]>0)) {
        if((ADC_TOF1[6] > 0) || (ADC_TOF1[18] > 0) || (ADC_TOF1[7] > 0) || (ADC_TOF1[19] > 0) || (ADC_TOF1[5] > 0) || (ADC_TOF1[17] > 0)) has_data_TOF++;
      }
    }

    else if (kk == 11) {
      if ((ADC_TOF2[11]>0) || (ADC_TOF2[35]>0) || (ADC_TOF2[23]>0) || (ADC_TOF2[47]>0)) {
        if((ADC_TOF1[11] > 0) || (ADC_TOF1[23] > 0) || (ADC_TOF1[0] > 0) || (ADC_TOF1[12] > 0) || (ADC_TOF1[10] > 0) || (ADC_TOF1[22] > 0)) has_data_TOF++;
      }
    }
  
    else {
      if ((ADC_TOF2[kk]>0) || (ADC_TOF2[kk+24]>0) || (ADC_TOF2[kk+12]>0) || (ADC_TOF2[kk+36]>0))  {
        if((ADC_TOF1[kk] > 0) || (ADC_TOF1[kk+12] > 0) || (ADC_TOF1[kk+1] > 0) || (ADC_TOF1[kk+13] > 0) || (ADC_TOF1[kk-1] > 0) || (ADC_TOF1[kk+11] > 0)) has_data_TOF++;
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

  bool good_TOF = false;
  bool good_TARGET = false;
  bool good_MWPC = false;

  if (has_data_TARGET >= 4 && has_data_low_TARGET >= 3) good_TARGET = true;
  if ((C2X_hit >=1) && (C2Y_hit >=1) && (C3X_hit >=1) && 
      (C3Y_hit >=1) && (C4X_hit >=1) && (C4Y_hit >=1)) good_MWPC = true;
  if (has_data_TOF > 0) good_TOF = true;
  
  //if ((good_TARGET) && (good_MWPC) && (good_TOF)) {
    //cout << i << endl;
   // good_events++;

  for(Int_t ii=0; ii<128; ii++){
    fiber[ii]=ADC_High_SFT[par_temp[1][ii]];
    if(tdc_le_sft[par_temp[1][ii]][0]>TDC_min_SFT && tdc_le_sft[par_temp[1][ii]][0]<TDC_max_SFT)
      fiberTDC[ii]=tdc_le_sft[par_temp[1][ii]][0];
    else if(tdc_le_sft[par_temp[1][ii]][1]>TDC_min_SFT && tdc_le_sft[par_temp[1][ii]][1]<TDC_max_SFT)
      fiberTDC[ii]=tdc_le_sft[par_temp[1][ii]][1];
    else if(tdc_le_sft[par_temp[1][ii]][2]>TDC_min_SFT && tdc_le_sft[par_temp[1][ii]][2]<TDC_max_SFT)
      fiberTDC[ii]=tdc_le_sft[par_temp[1][ii]][2];
    else if(tdc_le_sft[par_temp[1][ii]][3]>TDC_min_SFT && tdc_le_sft[par_temp[1][ii]][3]<TDC_max_SFT)
      fiberTDC[ii]=tdc_le_sft[par_temp[1][ii]][3];
    else if(tdc_le_sft[par_temp[1][ii]][4]>TDC_min_SFT && tdc_le_sft[par_temp[1][ii]][4]<TDC_max_SFT)
      fiberTDC[ii]=tdc_le_sft[par_temp[1][ii]][4];
    else if(tdc_le_sft[par_temp[1][ii]][5]>TDC_min_SFT && tdc_le_sft[par_temp[1][ii]][5]<TDC_max_SFT)
      fiberTDC[ii]=tdc_le_sft[par_temp[1][ii]][5];
    else
      fiberTDC[ii]=tdc_le_sft[par_temp[1][ii]][6];

    //  cout << j << "  " << count_target << "  " << ii << "   " << TDC_LE_sft[par_temp[1][ii]][0] << "   " << TDC_LE_sft[par_temp[1][ii]][1] << "   " << TDC_LE_sft[par_temp[1][ii]][2] << " ||  " << fiberTDC[ii] << endl;;
  
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

    for(Int_t imark3=0; imark3<2; imark3++){
      if((fiber[imark3+30]!=0) && (fiberTDC[imark3+30]>TDC_min_SFT) && (fiberTDC[imark3+30]<TDC_max_SFT)){  
      L3[imark3]++;
      }

      if((fiber[imark3+94]!=0) && (fiberTDC[imark3+94]>TDC_min_SFT) && (fiberTDC[imark3+94]<TDC_max_SFT)){  
      L3[imark3]++;
      }
    }


    for(Int_t imark3=2; imark3<17; imark3++){
      if((fiber[imark3+47]!=0) && (fiberTDC[imark3+47]>TDC_min_SFT) && (fiberTDC[imark3+47]<TDC_max_SFT)){  
      L3[18-imark3]++;
      }

      if((fiber[imark3+111]!=0) && (fiberTDC[imark3+111]>TDC_min_SFT) && (fiberTDC[imark3+111]<TDC_max_SFT)){  
      L3[18-imark3]++;
      }
    }


    for(Int_t imark4=0; imark4<17; imark4++){
      if((fiber[imark4+32]!=0) && (fiberTDC[imark4+32]>TDC_min_SFT) && (fiberTDC[imark4+32]<TDC_max_SFT)){  
      L4[16-imark4]++;
      }
    
      if((fiber[imark4+96]!=0) && (fiberTDC[imark4+96]>TDC_min_SFT) && (fiberTDC[imark4+96]<TDC_max_SFT)){  
      L4[16-imark4]++;
      }
    }

    Int_t L1_L2 = 0;
    Int_t L3_L4 = 0;

    for (int q=0; q<15; q++) {
      if (L1[q]>0) {
        for (int r=0; r<15; r++) {
          if (L2[r]>0) {
            L1_L2 = q-r;
            if (L1_L2 < 0) {
              L1_L2 = L1_L2 + 15;
              L1_L2count[L1_L2]++;
            }
            else {
              L1_L2count[L1_L2]++;
            }
          }
        }
      }
    }

    for (int q=0; q<17; q++) {
      if (L3[q]>0) {
        for (int r=0; r<17; r++) {
          if (L4[r]>0) {
            L3_L4 = q-r;
            if (L3_L4 < 0) {
              L3_L4 = L3_L4 + 17;
              L3_L4count[L3_L4]++;
            }
            else {
              L3_L4count[L3_L4]++;
            }
          }
        }
      }
    }
  //}
}

TH1D *h_L1L2count;  char Title_h_L1L2count[100];  char Name_h_L1L2count[100];
TH1D *h_L3L4count; char Title_h_L3L4count[100]; char Name_h_L3L4count[100];

double eff = double(double(good_events)/double(nentries))*100;

sprintf(Title_h_L1L2count,"SFT Layer 1-2 Distances - Run %d  |  ADC > %d , %d < TDC < %d", Run_Number, ADC_cut_SFT, TDC_min_SFT, TDC_max_SFT); 
sprintf(Name_h_L1L2count,"Layer 1-2 Distances"); 
h_L1L2count = new TH1D(Name_h_L1L2count,Title_h_L1L2count,17,-1,16);

sprintf(Title_h_L3L4count,"SFT Layer 3-4 Distances - Run %d  |  ADC > %d , %d < TDC < %d", Run_Number, ADC_cut_SFT, TDC_min_SFT, TDC_max_SFT); 
sprintf(Name_h_L3L4count,"Layer 3-4 Distances"); 
h_L3L4count = new TH1D(Name_h_L3L4count,Title_h_L3L4count,19,-1,18);

for (int i=0; i<15; i++) {
  h_L1L2count->Fill(i,L1_L2count[i]);
}

for (int i=0; i<17; i++) {
  h_L3L4count->Fill(i,L3_L4count[i]);
}

char Name_Can[100];     char Title_Can[100];
sprintf(Name_Can,"Run %d -- SFT_Layer_Distances.C",Run_Number);
sprintf(Title_Can,"Run %d  -- SFT_Layer_Distances.C",Run_Number);

TCanvas *c1;
c1 = new TCanvas(Name_Can,Title_Can,1000,500); 
c1->Divide(2,1);
c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

c1->cd(1);
h_L1L2count->Draw();

c1->cd(2);
h_L3L4count->Draw();

}
