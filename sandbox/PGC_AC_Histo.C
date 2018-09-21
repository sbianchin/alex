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
#include "ANAPATH.h"
#include "Thresholds.h"
//#include "Event_Display_MS.h"
#endif
  
void PGC_AC_Histo(Int_t Run_Number=5, Int_t time_window=1, Int_t flag=0){ 

gStyle->Clear();
TH1::AddDirectory(kFALSE);
gStyle->SetOptStat(1111111);

char source_mapping[]="SFT_Mapping_Oct14.txt";  // Mapping file !!!

Int_t adc_high_target[256];   		Int_t ADC_High_TARGET[256];    	Double_t ADC_High_TARGET_corr[256];   	
Int_t adc_low_target[256]; 			Int_t ADC_Low_TARGET[256]; 	//Double_t ADC_Low_TARGET_corr[256];
Int_t tdc_le_target[256][16];      	Int_t TDC_LE_TARGET[256];    	Double_t TDC_LE_TARGET_corr[256];   	
Int_t tdc_te_target[256][16]; 		Int_t TDC_TE_TARGET[256]; 	//Double_t TDC_TE_TARGET_corr[256];

Int_t adc_high_sft[128];   		Int_t ADC_High_SFT[128];    	Double_t ADC_High_SFT_corr[128];   	
Int_t adc_low_sft[128]; 		Int_t ADC_Low_SFT[128]; 	//Double_t ADC_Low_SFT_corr[128];
Int_t tdc_le_sft[128][16];      	Int_t TDC_LE_SFT[128];    	//Double_t TDC_LE_SFT_corr[128];   	
Int_t tdc_te_sft[128][16]; 		Int_t TDC_TE_SFT[128]; 		//Double_t TDC_TE_SFT_corr[128];	

Int_t ADC_TOF1[24];	
Int_t ADC_TOF2[56];

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

Int_t ADC_PGC[96];
Int_t TDC_PGC[84][16];

Int_t AC_ADC[24];
Int_t ADC_acu[12];
Int_t ADC_acd[12];

Int_t TDC_ACU[12][16]; 
Int_t TDC_ACD[12][16];

Int_t fiber[128]={-1};
Int_t TOF_Hit_Counter[12] = {0};

Int_t AC_Hit_Counter[24] = {0};
Int_t AC_Gap_Hit_Counter[12] = {0};

Int_t good_gap_counter[12] = {0};

Int_t PGC_Hit_Counter[96] = {0};
Int_t PGC_Gap_Counter[12]= {0};

int PGC_thr = 875;
int AC_thr = 800;

int TDC_AC_min = 700;
int TDC_AC_max = 1500;

int TDC_PGC_min = 700;
int TDC_PGC_max = 1500;

////////////////////////////////////////////////////////////////////////////////////////////////////


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

//fChain->SetBranchAddress("ADC_TOF1",ADC_TOF1);
//fChain->SetBranchAddress("ADC_TOF2",ADC_TOF2);

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

//fChain->SetBranchAddress("AC_ADC",AC_ADC);
fChain->SetBranchAddress("ADC_ACU",ADC_acu);
fChain->SetBranchAddress("ADC_ACD",ADC_acd);
fChain->SetBranchAddress("TDC_ACU",TDC_ACU);
fChain->SetBranchAddress("TDC_ACD",TDC_ACD);

fChain->SetBranchAddress("ADC_PGC",ADC_PGC);
fChain->SetBranchAddress("TDC_PGC",TDC_PGC);

fChain->SetBranchAddress("MWPCADC",MwpcADC);

cout << "TDC TOF1 Thresholds" << endl;
cout << TDC_TOF1_min << endl;
cout << TDC_TOF1_max << endl;

cout << "TDC TOF2 Thresholds" << endl;
cout << TDC_TOF2_min << endl;
cout << TDC_TOF2_max << endl;

//Int_t nentries;
Int_t nentries = (Int_t)fChain->GetEntries();//		Int_t nentries_SFT = (Int_t)fChain->GetEntries();
//if(nentries_TARGET==nentries_SFT) nentries = nentries_TARGET;
cout << "  " << endl;
cout << "****  Number of events: " << nentries << "  **** " <<endl;
cout << "  " << endl;


cout << "   " << endl;

if (flag!=0) nentries = flag;
for(Int_t i=0; i<nentries; i++){
//for(Int_t i=0; i<100; i++){
	fChain->GetEntry(i);	//	fChain_SFT->GetEntry(i);

	if(i%10000==1)	cout<<"**** "<<i<<" events done"<<endl;

//	cout << " " << endl;
//	cout << "********* Event " << i << endl;

	int has_data_TOF[12] = {0};
	int has_PGC_hit[12] = {0};
  	
  	for (int j_TARGET=0; j_TARGET<256; j_TARGET++){
  		ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-par_temp_TARGET[1][j_TARGET];
		ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-par_temp_TARGET[1][j_TARGET];
		TDC_LE_TARGET[j_TARGET]=tdc_le_target[j_TARGET][0];
		TDC_TE_TARGET[j_TARGET]=tdc_te_target[j_TARGET][0];
  	
//    	cout << j_TARGET << "   " << ADC_High_TARGET[j_TARGET] << "   " << TDC_LE_TARGET[j_TARGET] << endl;
  	}	


 	for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
		ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-par_temp_SFT[1][j_SFT];
		ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-par_temp_SFT[1][j_SFT];
		TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
		TDC_TE_SFT[j_SFT]=tdc_te_sft[j_SFT][0];
	
	//	cout << j_SFT << "   " << ADC_High_SFT[j_SFT] << "   " << TDC_LE_SFT[j_SFT] << endl;
	}

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


	//for (Int_t j_TOF1=0; j_TOF1<24; j_TOF1++) {
//		cout << j_TOF1 << " -- " << ADC_TOF1[j_TOF1] << endl;
	//	ADC_TOF1[j_TOF1] = ADC_TOF1[j_TOF1]-par_temp_TOF1[1][j_TOF1];
	//}

	//for (Int_t j_TOF2=0; j_TOF2<56; j_TOF2++) {
	//	ADC_TOF2[j_TOF2] = ADC_TOF2[j_TOF2]-par_temp_TOF2[1][j_TOF2];
	//}

	for (Int_t j_MWPC=0; j_MWPC<512; j_MWPC++) {
		MWPCADC[j_MWPC] = MwpcADC[j_MWPC]-MWPC_thr;
	}

	for (Int_t j_PGC=0; j_PGC<96; j_PGC++) {
		ADC_PGC[j_PGC] = ADC_PGC[j_PGC] - PGC_thr;
	}

	for (Int_t j_AC=0; j_AC<12; j_AC++) {
		AC_ADC[j_AC] = ADC_acu[j_AC]-AC_thr;
		AC_ADC[j_AC+12] = ADC_acd[j_AC]-AC_thr;
	}


	//for (Int_t j_AC=0; j_AC<24; j_AC++) {
	//	AC_ADC[j_AC] = AC_ADC[j_AC] - AC_thr;
	//}

	for (Int_t j_TDCTOF=0; j_TDCTOF<12; j_TDCTOF++) {
	    TDC_TOF1U[j_TDCTOF] = TDC_tof1U[j_TDCTOF];
	    TDC_TOF1D[j_TDCTOF] = TDC_tof1D[j_TDCTOF];
	    TDC_TOF2AO[j_TDCTOF] = TDC_tof2AO[j_TDCTOF];
	    TDC_TOF2BO[j_TDCTOF] = TDC_tof2BO[j_TDCTOF];
	    TDC_TOF2AI[j_TDCTOF] = TDC_tof2AI[j_TDCTOF];
	    TDC_TOF2BI[j_TDCTOF] = TDC_tof2BI[j_TDCTOF];
  	}

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

	if ((C2X_hit >=1) && (C2Y_hit >=1) && (C3X_hit >=1) && 
			(C3Y_hit >=1) && (C4X_hit >=1) && (C4Y_hit >=1)) good_MWPC = true;

  	/// Determine good gap without time window

  	if (time_window == 0) {
	  	for (int kk=0; kk<12; kk++) {
			if (kk == 0) {
				if ((ADC_TOF2[0]>0) || (ADC_TOF2[24]>0) || (ADC_TOF2[12]>0) || (ADC_TOF2[36]>0)) {
					if((ADC_TOF1[0] > 0) || (ADC_TOF1[12] > 0) || (ADC_TOF1[1] > 0) || (ADC_TOF1[13] > 0) || (ADC_TOF1[11] > 0) || (ADC_TOF1[23] > 0)) {TOF_Hit_Counter[0]++; has_data_TOF[0]++;}
				} 
			}

			/// Replace ADC TOF2[18] with ADC TOF2[55]
			//else if (kk == 6) {
			//	if ((ADC_TOF2[6]>0) || (ADC_TOF2[30]>0) || (ADC_TOF2[55]>0) || (ADC_TOF2[42]>0)) {
			//		if((ADC_TOF1[6] > 0) || (ADC_TOF1[18] > 0) || (ADC_TOF1[7] > 0) || (ADC_TOF1[19] > 0) || (ADC_TOF1[5] > 0) || (ADC_TOF1[17] > 0)) {TOF_Hit_Counter[6]++; has_data_TOF[6]++;}
			//	}
			//}

			else if (kk == 11) {
				if ((ADC_TOF2[11]>0) || (ADC_TOF2[35]>0) || (ADC_TOF2[23]>0) || (ADC_TOF2[47]>0)) {
					if((ADC_TOF1[11] > 0) || (ADC_TOF1[23] > 0) || (ADC_TOF1[0] > 0) || (ADC_TOF1[12] > 0) || (ADC_TOF1[10] > 0) || (ADC_TOF1[22] > 0)) {TOF_Hit_Counter[11]++; has_data_TOF[11]++;}
				}
			}
		
			else {
				if ((ADC_TOF2[kk]>0) || (ADC_TOF2[kk+24]>0) || (ADC_TOF2[kk+12]>0) || (ADC_TOF2[kk+36]>0))  {
					if((ADC_TOF1[kk] > 0) || (ADC_TOF1[kk+12] > 0) || (ADC_TOF1[kk+1] > 0) || (ADC_TOF1[kk+13] > 0) || (ADC_TOF1[kk-1] > 0) || (ADC_TOF1[kk+11] > 0)) {TOF_Hit_Counter[kk]++; has_data_TOF[kk]++;}
				}
			}
		}
	}

	bool has_ADC_TOF1_hit[12] = {false};
	bool has_TDC_TOF1_hit[12] = {false};
	bool has_ADC_TOF2_hit[12] = {false};
	bool has_TDC_TOF2_hit[12] = {false};

	/// Determine good gap with time window

	if (time_window != 0) {
		for (int kk=0; kk<12; kk++) {

			/// Replace ADC TOF2[18] with ADC TOF2[55]
			//if (kk == 6) {
			//	if ((ADC_TOF2[6]>0) || (ADC_TOF2[30]>0) || (ADC_TOF2[55]>0) || (ADC_TOF2[42]>0)) {has_ADC_TOF2_hit[6]=true;}
			//	if (((TDC_TOF2AO[6]>TDC_TOF2_min && TDC_TOF2AO[6] < TDC_TOF2_max) || (TDC_TOF2AI[6]>TDC_TOF2_min && TDC_TOF2AI[6] < TDC_TOF2_max)) 
			//    || ((TDC_TOF2BO[6]>TDC_TOF2_min && TDC_TOF2BO[6] < TDC_TOF2_max) || (TDC_TOF2BI[6]>TDC_TOF2_min && TDC_TOF2BI[6] < TDC_TOF2_max))) {has_TDC_TOF2_hit[6]=true;}

			//	if (ADC_TOF1[6]>0 || ADC_TOF1[18]>0) {has_ADC_TOF1_hit[6] = true;}
			//  	if ((TDC_TOF1U[6]>TDC_TOF1_min && TDC_TOF1U[6]<TDC_TOF1_max) || (TDC_TOF1D[6]>TDC_TOF1_min && TDC_TOF1D[6]<TDC_TOF1_max)) {has_TDC_TOF1_hit[6] = true;}
		  	//}

			//else {
				if ((ADC_TOF2[kk]>0) || (ADC_TOF2[kk+24]>0) || (ADC_TOF2[kk+12]>0) || (ADC_TOF2[kk+36]>0)) {has_ADC_TOF2_hit[kk]=true;}
				if (((TDC_TOF2AO[kk]>TDC_TOF2_min && TDC_TOF2AO[kk] < TDC_TOF2_max) || (TDC_TOF2AI[kk]>TDC_TOF2_min && TDC_TOF2AI[kk] < TDC_TOF2_max)) 
			    || ((TDC_TOF2BO[kk]>TDC_TOF2_min && TDC_TOF2BO[kk] < TDC_TOF2_max) || (TDC_TOF2BI[kk]>TDC_TOF2_min && TDC_TOF2BI[kk] < TDC_TOF2_max))) {has_TDC_TOF2_hit[kk]=true;}

				if (ADC_TOF1[kk]>0 || ADC_TOF1[kk+12]>0) {has_ADC_TOF1_hit[kk] = true;}
			  	if ((TDC_TOF1U[kk]>TDC_TOF1_min && TDC_TOF1U[kk]<TDC_TOF1_max) || (TDC_TOF1D[kk]>TDC_TOF1_min && TDC_TOF1D[kk]<TDC_TOF1_max)) {has_TDC_TOF1_hit[kk] = true;}
			//}
		}

		for (int kk=0; kk<12; kk++) {
			if (kk == 0) {
				if (has_ADC_TOF2_hit[0] && has_TDC_TOF2_hit[0]) {
					if ((has_ADC_TOF1_hit[0] && has_TDC_TOF1_hit[0]) || (has_ADC_TOF1_hit[11] && has_TDC_TOF1_hit[11]) || (has_ADC_TOF1_hit[1] && has_TDC_TOF1_hit[1])) {
						TOF_Hit_Counter[0]++;
						has_data_TOF[0]++;
					}
				}
			}
			else if (kk == 11) {
				if (has_ADC_TOF2_hit[11] && has_TDC_TOF2_hit[11]) {
					if ((has_ADC_TOF1_hit[11] && has_TDC_TOF1_hit[11]) || (has_ADC_TOF1_hit[10] && has_TDC_TOF1_hit[10]) || (has_ADC_TOF1_hit[0] && has_TDC_TOF1_hit[0])) {
						TOF_Hit_Counter[11]++;
						has_data_TOF[11]++;
					}
				}
			}
			else {
				if (has_ADC_TOF2_hit[kk] && has_TDC_TOF2_hit[kk]) {
					if ((has_ADC_TOF1_hit[kk] && has_TDC_TOF1_hit[kk]) || (has_ADC_TOF1_hit[kk-1] && has_TDC_TOF1_hit[kk-1]) || (has_ADC_TOF1_hit[kk+1] && has_TDC_TOF1_hit[kk+1])) {
						TOF_Hit_Counter[kk]++;
						has_data_TOF[kk]++;
					}
				}
			}
		}
	}

	bool has_TOF_hit = false;

	/*
	cout << "ADC PGC -- Event " << i << endl;
	for (int j=0; j<96; j++) {
			cout << ADC_PGC[j] << endl;
	}

	cout << "" << endl;
	cout << "ADC AC -- Event " << i << endl;
	for (int j=0; j<24; j++) {
			cout << AC_ADC[j] << endl;
	}
	cout << "" << endl;

	cout << "TDC ACU -- Event " << i << endl;
	for (int j=0; j<12; j++) {
			cout << TDC_ACU[j] << endl;
	}
	cout << "" << endl;

	cout << "TDC ACD -- Event " << i << endl;
	for (int j=0; j<12; j++) {
			cout << TDC_ACD[j] << endl;
	}
	cout << "" << endl;

	*/

	for (int q=0; q<12; q++) {
		if (has_data_TOF[q] > 0) has_TOF_hit = true;
	}

	bool has_TDC_ACU_hit[12] = {false};
	bool has_TDC_ACD_hit[12] = {false};
	bool has_TDC_PGC_hit[84] = {false};

	/// Filter TDC hits

	for (int q=0; q<12; q++) {
		for (int j=0; j<4; j++) {
			if (TDC_ACU[q][j] > TDC_AC_min && TDC_ACU[q][j] < TDC_AC_max) has_TDC_ACU_hit[q] = true;
			if (TDC_ACD[q][j] > TDC_AC_min && TDC_ACD[q][j] < TDC_AC_max) has_TDC_ACD_hit[q] = true;
		}
	}

	for (int q=0; q<84; q++) {
		for (int j=0; j<4; j++) {
			if (TDC_PGC[q][j] > TDC_PGC_min && TDC_PGC[q][j] < TDC_PGC_max) has_TDC_PGC_hit[q] = true;
		}
	}
	
	/// Fill Raw AC & PGC Histograms

	for (int j=0; j<96; j++) {
		if ((j==7) || (j==15) || (j==23) || (j==31) || (j==39) || (j==47) || (j==55) || (j==63) || (j==71) || (j==79)
		|| (j==87) || (j==95)) continue;         
		if (ADC_PGC[j] > 0) PGC_Hit_Counter[j]++;
	}
	for (int j=0; j<24; j++) {
		if (AC_ADC[j] > 0) AC_Hit_Counter[j]++;
	}

	/// Fill AC Gap Histogram

	for (int j=0; j<12; j++) {
		if (has_data_TOF[j] > 0 && good_MWPC) {
			if (((AC_ADC[j] > 0) && (has_TDC_ACU_hit[j])) || ((AC_ADC[j+12] > AC_thr) && (has_TDC_ACD_hit[j]))) AC_Gap_Hit_Counter[j]++;
		}
	}

	/// Fill PGC Gap Histogram

	for (int j=0; j<96; j++) {
		if (j>=0 && j<7) {
			if ((ADC_PGC[j] > 0) && (has_TDC_PGC_hit[j+21])) has_PGC_hit[0]++;
		}

		if (j>=8 && j<15) {
			if ((ADC_PGC[j] > 0) && (has_TDC_PGC_hit[j+6])) has_PGC_hit[1]++;
		}

		if (j>=16 && j<23) {
			if ((ADC_PGC[j] > 0) && (has_TDC_PGC_hit[j-9])) has_PGC_hit[2]++;
		}

		if (j>=24 && j<31) {
			if ((ADC_PGC[j] > 0) && (has_TDC_PGC_hit[j-24])) has_PGC_hit[3]++;
		}

		if (j>=32 && j<39) {
			if ((ADC_PGC[j] > 0) && (has_TDC_PGC_hit[j+3])) has_PGC_hit[4]++;
		}

		if (j>=40 && j<47) {
			if ((ADC_PGC[j] > 0) && (has_TDC_PGC_hit[j-12])) has_PGC_hit[5]++;
		}
		
		if (j>=48 && j<55) {
			if ((ADC_PGC[j] > 0) && (has_TDC_PGC_hit[j+15])) has_PGC_hit[6]++;
		}
		
		if (j>=56 && j<63) {
			if ((ADC_PGC[j] > 0) && (has_TDC_PGC_hit[j])) has_PGC_hit[7]++;
		}
			
		if (j>=64 && j<71) {
			if ((ADC_PGC[j] > 0) && (has_TDC_PGC_hit[j-15])) has_PGC_hit[8]++;	
		}

		if (j>=72 && j<79) {
			if ((ADC_PGC[j] > 0) && (has_TDC_PGC_hit[j-30])) has_PGC_hit[9]++;
		}

		if (j>=80 && j<87) {
			if ((ADC_PGC[j] > 0) && (has_TDC_PGC_hit[j-3])) has_PGC_hit[10]++;
		}

		if (j>=88 && j<95) {
			if ((ADC_PGC[j] > 0) && (has_TDC_PGC_hit[j-18])) has_PGC_hit[11]++;
		}
	}
	
	for (int j=0; j<12; j++) {
		if (has_PGC_hit[j] > 0 && good_MWPC) {
			PGC_Gap_Counter[j]++;
		}
	}

}

TH1D *h_count;	char Title_h_count[100];	char Name_h_count[100];
sprintf(Title_h_count,"ADC PGC Hit Multiplicity -- Run %d", Run_Number); 
sprintf(Name_h_count,"ADC PGC Hit Multiplicity -- Run %d", Run_Number);
h_count = new TH1D(Name_h_count,Title_h_count,98,0,98);

for (int i=0; i<96; i++) {
	h_count->Fill(i+1,PGC_Hit_Counter[i]);
}

TH1D *h_count2;	char Title_h_count2[100];	char Name_h_count2[100];
sprintf(Title_h_count2,"ADC AC Hit Multiplicity -- Run %d", Run_Number); 
sprintf(Name_h_count2,"ADC AC Hit Multiplicity -- Run %d", Run_Number);
h_count2 = new TH1D(Name_h_count2,Title_h_count2,26,0,26);


for (int i=0; i<24; i++) {
	h_count2->Fill(i+1,AC_Hit_Counter[i]);
}

TH1D *h_count3;	char Title_h_count3[100];	char Name_h_count3[100];
sprintf(Title_h_count3,"ADC PGC Gap Hit Multiplicity -- Run %d", Run_Number); 
sprintf(Name_h_count3,"ADC PGC Gap Hit Multiplicity -- Run %d", Run_Number);
h_count3 = new TH1D(Name_h_count3,Title_h_count3,14,0,14);

for (int i=0; i<12; i++) {
	h_count3->Fill(i+1,PGC_Gap_Counter[i]);
}

TH1D *h_count4;	char Title_h_count4[100];	char Name_h_count4[100];
sprintf(Title_h_count4,"ADC AC Gap Hit Multiplicity -- Run %d", Run_Number); 
sprintf(Name_h_count4,"ADC AC Gap Hit Multiplicity -- Run %d", Run_Number);
h_count4 = new TH1D(Name_h_count4,Title_h_count4,14,0,14);

for (int i=0; i<12; i++) {
//	cout << "xxxxxx   " << AC_Gap_Hit_Counter[i] << endl;
	h_count4->Fill(i+1,AC_Gap_Hit_Counter[i]);
}

/////////////////

char Name_Can[100];			char Title_Can[100];
sprintf(Name_Can,"ADC AC & PGC Hit Multiplicity -- Run %d", Run_Number);
sprintf(Title_Can,"ADC AC & PGC Hit Multiplicity -- Run %d", Run_Number);


TCanvas *c1;
c1 = new TCanvas(Name_Can,Title_Can,1000,1000);
c1->Divide(2,2); 
c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

c1->cd(1);
h_count->Draw();

c1->cd(2);
h_count2->Draw();

c1->cd(3);
h_count3->Draw();

c1->cd(4);
h_count4->Draw(); 


return;
}
