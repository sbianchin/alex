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
//#include "ListGoodEvents.h"
#endif
  
void TARGET_Hit_Multiplicity(Int_t Run_Number=5, Int_t flag=0, Int_t ievt=0){ 

gStyle->SetOptStat(0);

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
Int_t fiber[128]={-1};

Int_t HitCounter0 = 0;
Int_t HitCounter1 = 0;
Int_t HitCounter2 = 0;
Int_t HitCounter3 = 0;
Int_t HitCounter4 = 0;
Int_t HitCounter5 = 0;

Int_t LG_Mult0[20] = {0};
Int_t LG_Mult1[20] = {0};
Int_t LG_Mult2[20] = {0};
Int_t LG_Mult3[20] = {0};
Int_t LG_Mult4[20] = {0};
Int_t LG_Mult5[20] = {0};

int good_event_counter = 0;

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

char footer[100];
sprintf(footer,"Event_Display_MS.C  --  Run %d ; Event %d",Run_Number,ievt);

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

//for(Int_t jj=0; jj<256; jj++){
//	cout << par_temp3[0][jj] << "   " << par_temp3[1][jj] << endl;
//}

////

/*

TChain *fChain_TARGET = new TChain("Tree_TARGET");						TChain *fChain_SFT = new TChain("Tree_SFT");
fChain_TARGET->Add(Name_finput);										fChain_SFT->Add(Name_finput);
fChain_TARGET->SetMakeClass(1);											fChain_SFT->SetMakeClass(1);

fChain_TARGET->SetBranchAddress("ADC_High_TARGET",adc_high_target);		fChain_SFT->SetBranchAddress("ADC_High_SFT",adc_high_sft);
fChain_TARGET->SetBranchAddress("ADC_Low_TARGET",adc_low_target);		fChain_SFT->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
fChain_TARGET->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);			fChain_SFT->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
fChain_TARGET->SetBranchAddress("TDC_TE_TARGET",tdc_te_target);			fChain_SFT->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

*/

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

//	cout << " " << endl;
//	cout << "********* Event " << i << endl;

	int has_data_TARGET = 0;
	int has_data_TOF = 0;

	int LG0_count = 0;
	int LG1_count = 0;
	int LG2_count = 0;
	int LG3_count = 0;
	int LG4_count = 0;
	int LG5_count = 0;

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
  	
//    	cout << j_TARGET << "   " << ADC_High_TARGET[j_TARGET] << "   " << TDC_LE_TARGET[j_TARGET] << endl;
  	}	


 	for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
		ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-par_temp_SFT[1][j_SFT];
		ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-par_temp_SFT[1][j_SFT];
		TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
		TDC_TE_SFT[j_SFT]=tdc_te_sft[j_SFT][0];
	
	//	cout << j_SFT << "   " << ADC_High_SFT[j_SFT] << "   " << TDC_LE_SFT[j_SFT] << endl;
	}

	//for (Int_t j_TOF1=0; j_TOF1<24; j_TOF1++) {
//		cout << j_TOF1 << " -- " << ADC_TOF1[j_TOF1] << endl;
	//	ADC_TOF1[j_TOF1] = ADC_TOF1[j_TOF1]-par_temp_TOF1[1][j_TOF1];
	//}

	//for (Int_t j_TOF2=0; j_TOF2<56; j_TOF2++) {
	//	ADC_TOF2[j_TOF2] = ADC_TOF2[j_TOF2]-par_temp_TOF2[1][j_TOF2];
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

	/// Determine good gap condition with time window

	bool has_ADC_TOF1_hit[12] = {false};
	bool has_TDC_TOF1_hit[12] = {false};
	bool has_ADC_TOF2_hit[12] = {false};
	bool has_TDC_TOF2_hit[12] = {false};

	for (int kk=0; kk<12; kk++) {

		/// Replace ADC TOF2[18] with ADC TOF2[55]
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

	/// Determine MWPC chambers hit

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

	if ((has_data_TOF > 0) && (good_MWPC)) good_event_counter++;

	for (Int_t p=0; p<256; p++) {
		if (ADC_High_TARGET[p] > 0) has_data_TARGET++;
	}

	for (Int_t p=0; p<128; p++) {
		fiber[p]=ADC_High_SFT[par_temp[1][p]];
	}

	if ((has_data_TOF > 0) && (good_MWPC)) {
		if (has_data_TARGET == 0) {
			HitCounter0++;
			for (int i=0; i<256; i++) {
				if (ADC_Low_TARGET[i] > 0) LG0_count++;
			}
			LG_Mult0[LG0_count]++;
		}

		if (has_data_TARGET == 1) {
			HitCounter1++;
			for (int i=0; i<256; i++) {
				if (ADC_Low_TARGET[i] > 0) LG1_count++;
			}
			LG_Mult1[LG1_count]++;
		}

		if (has_data_TARGET == 2) {
			HitCounter2++;
			for (int i=0; i<256; i++) {
				if (ADC_Low_TARGET[i] > 0) LG2_count++;
			}
			LG_Mult2[LG2_count]++;
		}

		if (has_data_TARGET == 3) {
			HitCounter3++;
			for (int i=0; i<256; i++) {
				if (ADC_Low_TARGET[i] > 0) LG3_count++;
			}
			LG_Mult3[LG3_count]++;
		}

		if (has_data_TARGET == 4) {
			HitCounter4++;
			for (int i=0; i<256; i++) {
				if (ADC_Low_TARGET[i] > 0) LG4_count++;
			}
			LG_Mult4[LG4_count]++;
		}

		if (has_data_TARGET == 5) {
			HitCounter5++;
			for (int i=0; i<256; i++) {
				if (ADC_Low_TARGET[i] > 0) LG5_count++;
			}
			LG_Mult5[LG5_count]++;
		}
	}
}

double eff = (double(double(HitCounter0)/double(nentries)))*100;
double eff2 = (double(double(HitCounter1)/double(nentries)))*100;
double eff3 = (double(double(HitCounter2)/double(nentries)))*100;
double eff4 = (double(double(HitCounter3)/double(nentries)))*100;
double eff5 = (double(double(HitCounter4)/double(nentries)))*100;
double eff6 = (double(double(HitCounter5)/double(nentries)))*100;

double geff = (double(double(HitCounter0)/double(good_event_counter)))*100;
double geff2 = (double(double(HitCounter1)/double(good_event_counter)))*100;
double geff3 = (double(double(HitCounter2)/double(good_event_counter)))*100;
double geff4 = (double(double(HitCounter3)/double(good_event_counter)))*100;
double geff5 = (double(double(HitCounter4)/double(good_event_counter)))*100;
double geff6 = (double(double(HitCounter5)/double(good_event_counter)))*100;

cout << "" << endl;
cout << "Number of events with good gap & MWPC: " << good_event_counter << endl;

cout << "" << endl;
cout << "0 Target Hits: " << HitCounter0 << " -- % of total events: " << eff << "%" << " -- % of good gap events: " << geff << "%" << endl;
cout << "1 Target Hits: " << HitCounter1 << " -- % of total events: " << eff2 << "%" << " -- % of good gap events: " << geff2 << "%" << endl;
cout << "2 Target Hits: " << HitCounter2 << " -- % of total events: " << eff3 << "%" << " -- % of good gap events: " << geff3 << "%" << endl;
cout << "3 Target Hits: " << HitCounter3 << " -- % of total events: " << eff4 << "%" << " -- % of good gap events: " << geff4 << "%" << endl;
cout << "4 Target Hits: " << HitCounter4 << " -- % of total events: " << eff5 << "%" << " -- % of good gap events: " << geff5 << "%" << endl;
cout << "5 Target Hits: " << HitCounter5 << " -- % of total events: " << eff6 << "%" << " -- % of good gap events: " << geff6 << "%" << endl;
cout << "" << endl;

char Title_h_count_TARGET[7][100];	char Name_h_count_TARGET[7][100];
TH1D *h_count_TARGET[7];

for (int i=0; i<7; i++) {
	sprintf(Title_h_count_TARGET[i],"TARGET Multiplicities - Run %d | ADC > %d", Run_Number, ADC_cut_TARGET);
	sprintf(Name_h_count_TARGET[i],"Histogram Counter %d", i+1); 
	h_count_TARGET[i] = new TH1D(Name_h_count_TARGET[i],Title_h_count_TARGET[i],15,0,15);
	h_count_TARGET[i]->GetXaxis()->SetLabelOffset(999);
	h_count_TARGET[i]->GetXaxis()->SetLabelSize(0);
}

h_count_TARGET[0]->Fill(1,nentries);
h_count_TARGET[1]->Fill(3,HitCounter0);
h_count_TARGET[2]->Fill(5,HitCounter1);
h_count_TARGET[3]->Fill(7,HitCounter2);
h_count_TARGET[4]->Fill(9,HitCounter3);
h_count_TARGET[5]->Fill(11,HitCounter4);
h_count_TARGET[6]->Fill(13,HitCounter5);

char Title_h_count_LG[6][100];	char Name_h_count_LG[6][100];
TH1D *h_count_LG[6];

for (int i=0; i<6; i++) {
	sprintf(Title_h_count_LG[i],"LG ADC Multiplicities: %d HG TARGET Hit(s) -- Run %d | ADC > %d", i, Run_Number, ADC_cut_TARGET);
	sprintf(Name_h_count_LG[i],"Histogram LG Counter %d", i+1); 
	h_count_LG[i] = new TH1D(Name_h_count_LG[i],Title_h_count_LG[i],20,-1,19);
}

for (int i=0; i<20; i++) {
	h_count_LG[0]->Fill(i,LG_Mult0[i]);
	h_count_LG[1]->Fill(i,LG_Mult1[i]);
	h_count_LG[2]->Fill(i,LG_Mult2[i]);
	h_count_LG[3]->Fill(i,LG_Mult3[i]);
	h_count_LG[4]->Fill(i,LG_Mult4[i]);
	h_count_LG[5]->Fill(i,LG_Mult5[i]);
}

char Name_Can2[100];			char Title_Can2[100];
sprintf(Name_Can2,"TARGET_Hit_Multiplicity -- Run %d", Run_Number);
sprintf(Title_Can2,"TARGET_Hit_Multiplicity -- Run %d", Run_Number);

TCanvas *c1;
c1 = new TCanvas(Name_Can2,Title_Can2,1000,600); 
c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

TLegend *l1;
l1 = new TLegend(0.75,0.75,0.99,0.99);
l1->AddEntry(h_count_TARGET[0], "Number of Events", "F");
l1->AddEntry(h_count_TARGET[1], "0 Target Hits", "F");
l1->AddEntry(h_count_TARGET[2], "1 Target Hit", "F");
l1->AddEntry(h_count_TARGET[3], "2 Target Hits", "F");
l1->AddEntry(h_count_TARGET[4], "3 Target Hits", "F");
l1->AddEntry(h_count_TARGET[5], "4 Target Hits", "F");
l1->AddEntry(h_count_TARGET[6], "5 Target Hits", "F");

h_count_TARGET[0]->SetLineColor(1);
h_count_TARGET[0]->SetFillStyle(3004);
h_count_TARGET[0]->SetFillColor(1);
h_count_TARGET[0]->Draw();

h_count_TARGET[1]->SetLineColor(2);
h_count_TARGET[1]->SetFillStyle(3004);
h_count_TARGET[1]->SetFillColor(2);
h_count_TARGET[1]->Draw("same");

h_count_TARGET[2]->SetLineColor(3);
h_count_TARGET[2]->SetFillStyle(3004);
h_count_TARGET[2]->SetFillColor(3);
h_count_TARGET[2]->Draw("same");

h_count_TARGET[3]->SetLineColor(4);
h_count_TARGET[3]->SetFillStyle(3004);
h_count_TARGET[3]->SetFillColor(4);
h_count_TARGET[3]->Draw("same");

h_count_TARGET[4]->SetLineColor(5);
h_count_TARGET[4]->SetFillStyle(3004);
h_count_TARGET[4]->SetFillColor(5);
h_count_TARGET[4]->Draw("same");

h_count_TARGET[5]->SetLineColor(6);
h_count_TARGET[5]->SetFillStyle(3004);
h_count_TARGET[5]->SetFillColor(6);
h_count_TARGET[5]->Draw("same");

h_count_TARGET[6]->SetLineColor(7);
h_count_TARGET[6]->SetFillStyle(3004);
h_count_TARGET[6]->SetFillColor(7);
h_count_TARGET[6]->Draw("same");

l1->Draw("same");

char Name_Can[100];			char Title_Can[100];
sprintf(Name_Can,"TARGET Low Gain Hit Multiplicities -- Run %d", Run_Number);
sprintf(Title_Can,"TARGET Low Gain Hit Multiplicities -- Run %d", Run_Number);

TCanvas *c2;
c2 = new TCanvas(Name_Can,Title_Can,1200,500);
c2->Divide(3,2); 
c2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

for (int i=0; i<6; i++) {
	c2->cd(i+1);
	h_count_LG[i]->Draw();
}

}