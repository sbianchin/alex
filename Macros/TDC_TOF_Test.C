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
  
void TDC_TOF_Test(Int_t run_number=5, Int_t normalize = 0, Int_t flag=0, Int_t QUIET=0) { 

gStyle->SetOptStat(1111111);

Int_t adc_high_target[256]; 	
Int_t adc_low_target[256]; 		
Int_t tdc_le_target[256][16];      	  	
Int_t tdc_te_target[256][16]; 		

Int_t ADC_High_sft[128];
Int_t ADC_Low_sft[128];
Int_t TDC_LE_sft[128][16];
Int_t TDC_TE_sft[128][16];

Int_t ADC_tof1[24];
Int_t ADC_tof2[56];

Int_t ADC_tof1U[12];  Int_t ADC_TOF1U[12];
Int_t ADC_tof1D[12];  Int_t ADC_TOF1D[12];

Int_t ADC_tof2AO[12];   Int_t ADC_TOF2AO[12];
Int_t ADC_tof2BO[12];   Int_t ADC_TOF2BO[12];
Int_t ADC_tof2AI[12];   Int_t ADC_TOF2AI[12];
Int_t ADC_tof2BI[12];   Int_t ADC_TOF2BI[12];

Int_t TDC_tof1U[12];
Int_t TDC_tof1D[12];

Int_t TDC_tof2AO[12];
Int_t TDC_tof2BO[12];
Int_t TDC_tof2AI[12];
Int_t TDC_tof2BI[12];

Int_t MwpcADC[512];

Int_t TOF1_Mapping1[12] = {1,3,5,7,9,11,13,15,17,19,21,23};
Int_t TOF1_Mapping2[12] = {2,4,6,8,10,12,14,16,18,20,22,24};

Int_t TOF2_Mapping1[12] = {1,5,9,13,17,21,25,29,33,37,41,45};
Int_t TOF2_Mapping2[12] = {2,6,10,14,18,22,26,30,34,38,42,46};
Int_t TOF2_Mapping3[12] = {3,7,11,15,19,23,27,31,35,39,43,47};
Int_t TOF2_Mapping4[12] = {4,8,12,16,20,24,28,32,36,40,44,48};	

char path_input[200];                   char file_mapping[200];
sprintf(path_input,"%s",path_merged);          
//sprintf(file_mapping,path_mapping);
//sprintf(path_input,"/triumfcs/trshare/trek/E36/Data/April_2015/root");

char Name_finput[200];
sprintf(Name_finput,"%s/Run%dMS.root",path_input, run_number);

char path_output[200];
sprintf(path_output,"%s",path_histo);

//char footer[100];
//sprintf(footer,"Event_Display_MS.C  --  Run %d ; Event %d",Run_Number,ievt);

cout << "   " << endl;
cout << Name_finput << endl;

TChain *fChain= new TChain("Tree");		
fChain->Add(Name_finput);		
fChain->SetMakeClass(1);							

fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);		
fChain->SetBranchAddress("ADC_Low_TARGET",adc_low_target);
fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);		
fChain->SetBranchAddress("TDC_TE_TARGET",tdc_te_target);

fChain->SetBranchAddress("ADC_High_SFT",ADC_High_sft);
fChain->SetBranchAddress("ADC_Low_SFT",ADC_Low_sft);
fChain->SetBranchAddress("TDC_LE_SFT",TDC_LE_sft);
fChain->SetBranchAddress("TDC_TE_SFT",TDC_TE_sft);

//fChain->SetBranchAddress("ADC_TOF1",ADC_tof1);
//fChain->SetBranchAddress("ADC_TOF2",ADC_tof2);
fChain->SetBranchAddress("ADC_TOF1U", ADC_tof1U);
fChain->SetBranchAddress("ADC_TOF1D", ADC_tof1D);

fChain->SetBranchAddress("ADC_TOF2AO", ADC_tof2AO);
fChain->SetBranchAddress("ADC_TOF2AI", ADC_tof2AI);
fChain->SetBranchAddress("ADC_TOF2BO", ADC_tof2BO);
fChain->SetBranchAddress("ADC_TOF2BI", ADC_tof2BI);


fChain->SetBranchAddress("MWPCADC",MwpcADC);

fChain->SetBranchAddress("TDC_TOF1U",TDC_tof1U);
fChain->SetBranchAddress("TDC_TOF1D",TDC_tof1D);
fChain->SetBranchAddress("TDC_TOF2AO",TDC_tof2AO);
fChain->SetBranchAddress("TDC_TOF2AI",TDC_tof2AI);
fChain->SetBranchAddress("TDC_TOF2BO",TDC_tof2BO);
fChain->SetBranchAddress("TDC_TOF2BI",TDC_tof2BI);			

TH1D *h_ADC_tof1[24];   char Title_ADC_tof1[24][100];	char Name_ADC_tof1[24][100];     
TH1D *h_ADC_tof2[56];	char Title_ADC_tof2[56][100];	char Name_ADC_tof2[56][100];

TH1D *h_TDC_tof1U[12];   char Title_TDC_tof1U[12][100];	char Name_TDC_tof1U[12][100];
TH1D *h_TDC_tof1D[12];   char Title_TDC_tof1D[12][100];	char Name_TDC_tof1D[12][100];
TH1D *h_TDC_tof2AO[12];   char Title_TDC_tof2AO[12][100];	char Name_TDC_tof2AO[12][100];
TH1D *h_TDC_tof2AI[12];   char Title_TDC_tof2AI[12][100];	char Name_TDC_tof2AI[12][100];
TH1D *h_TDC_tof2BO[12];   char Title_TDC_tof2BO[12][100];	char Name_TDC_tof2BO[12][100];
TH1D *h_TDC_tof2BI[12];   char Title_TDC_tof2BI[12][100];	char Name_TDC_tof2BI[12][100];

TH1D *h_TDC_tof2AO_peak[12];   char Title_TDC_tof2AO_peak[12][100];	char Name_TDC_tof2AO_peak[12][100];
TH1D *h_TDC_tof2AI_peak[12];   char Title_TDC_tof2AI_peak[12][100];	char Name_TDC_tof2AI_peak[12][100];
TH1D *h_TDC_tof2BO_peak[12];   char Title_TDC_tof2BO_peak[12][100];	char Name_TDC_tof2BO_peak[12][100];
TH1D *h_TDC_tof2BI_peak[12];   char Title_TDC_tof2BI_peak[12][100];	char Name_TDC_tof2BI_peak[12][100];

TH1D *h_MWPCADC[512];	char Title_MWPCADC[512][100];	char Name_MWPCADC[512][100];  

for (int i=0; i<24; i++) {
	sprintf(Title_ADC_tof1[i],"Raw ADC (Ch. %d)  --  TOF1",i);
	sprintf(Name_ADC_tof1[i],"ADC_TOF1 (Ch. %d)",i);
}

for (int i=0; i<56; i++) {
	sprintf(Title_ADC_tof2[i],"Raw ADC (Ch. %d)  --  TOF2",i);
	sprintf(Name_ADC_tof2[i],"ADC_TOF2 (Ch. %d)",i);
}

for (int i=0; i<12; i++) {
	sprintf(Title_TDC_tof1U[i],"Raw TDC (Ch. %d)  --  TOF1",i);
	sprintf(Name_TDC_tof1U[i],"TDC_TOF1 (Ch. %d)",i);

	sprintf(Title_TDC_tof1D[i],"Raw TDC (Ch. %d)  --  TOF1",i+12);
	sprintf(Name_TDC_tof1D[i],"TDC_TOF1 (Ch. %d)",i+12);

	sprintf(Title_TDC_tof2AO[i],"Raw TDC (Ch. %d)  --  TOF2",i);
	sprintf(Name_TDC_tof2AO[i],"TDC_TOF2 (Ch. %d)",i);

	sprintf(Title_TDC_tof2AI[i],"Raw TDC (Ch. %d)  --  TOF2",i+24);
	sprintf(Name_TDC_tof2AI[i],"TDC_TOF2 (Ch. %d)",i+24);

	sprintf(Title_TDC_tof2BO[i],"Raw TDC (Ch. %d)  --  TOF2",i+12);
	sprintf(Name_TDC_tof2BO[i],"TDC_TOF2 (Ch. %d)",i+12);

	sprintf(Title_TDC_tof2BI[i],"Raw TDC (Ch. %d)  --  TOF2",i+36);
	sprintf(Name_TDC_tof2BI[i],"TDC_TOF2 (Ch. %d)",i+36);

	sprintf(Title_TDC_tof2AO_peak[i],"Raw TDC (Ch. %d)  --  TOF2 | Peak Cut",i);
	sprintf(Name_TDC_tof2AO_peak[i],"TDC_TOF2 (Ch. %d)",i);

	sprintf(Title_TDC_tof2AI_peak[i],"Raw TDC (Ch. %d)  --  TOF2 | Peak Cut",i+24);
	sprintf(Name_TDC_tof2AI_peak[i],"TDC_TOF2 (Ch. %d)",i+24);

	sprintf(Title_TDC_tof2BO_peak[i],"Raw TDC (Ch. %d)  --  TOF2 | Peak Cut",i+12);
	sprintf(Name_TDC_tof2BO_peak[i],"TDC_TOF2 (Ch. %d)",i+12);

	sprintf(Title_TDC_tof2BI_peak[i],"Raw TDC (Ch. %d)  --  TOF2 | Peak Cut",i+36);
	sprintf(Name_TDC_tof2BI_peak[i],"TDC_TOF2 (Ch. %d)",i+36);
}

for (int  i=0; i<24; i++) {
h_ADC_tof1[i] = new TH1D(Name_ADC_tof1[i],Title_ADC_tof1[i],1000,0,4100);
}

for (int  i=0; i<56; i++) {
h_ADC_tof2[i] = new TH1D(Name_ADC_tof2[i],Title_ADC_tof2[i],1000,0,4100);
}


for (int  i=0; i<12; i++) {
h_TDC_tof1U[i] = new TH1D(Name_TDC_tof1U[i],Title_TDC_tof1U[i],1000,0,4000);
h_TDC_tof1D[i] = new TH1D(Name_TDC_tof1D[i],Title_TDC_tof1D[i],1000,0,4000);
h_TDC_tof2AO[i] = new TH1D(Name_TDC_tof2AO[i],Title_TDC_tof2AO[i],1000,0,4000);
h_TDC_tof2AI[i] = new TH1D(Name_TDC_tof2AI[i],Title_TDC_tof2AI[i],1000,0,4000);
h_TDC_tof2BO[i] = new TH1D(Name_TDC_tof2BO[i],Title_TDC_tof2BO[i],1000,0,4000);
h_TDC_tof2BI[i] = new TH1D(Name_TDC_tof2BI[i],Title_TDC_tof2BI[i],1000,0,4000);

h_TDC_tof2AO_peak[i] = new TH1D(Name_TDC_tof2AO_peak[i],Title_TDC_tof2AO_peak[i],1000,0,4000);
h_TDC_tof2AI_peak[i] = new TH1D(Name_TDC_tof2AI_peak[i],Title_TDC_tof2AI_peak[i],1000,0,4000);
h_TDC_tof2BO_peak[i] = new TH1D(Name_TDC_tof2BO_peak[i],Title_TDC_tof2BO_peak[i],1000,0,4000);
h_TDC_tof2BI_peak[i] = new TH1D(Name_TDC_tof2BI_peak[i],Title_TDC_tof2BI_peak[i],1000,0,4000);
}

//read all entries and fill the histograms
Int_t nentries = (Int_t)fChain->GetEntries();
cout <<  "Total Number of Entries :     " <<  nentries << endl;

cout << "  " << endl;

cout << "TDC TOF1 min: " << TDC_TOF1_min << endl;
cout << "TDC TOF1 max: " << TDC_TOF1_max << endl;
cout << "TDC TOF2 min: " << TDC_TOF2_min << endl;
cout << "TDC TOF2 max: " << TDC_TOF2_max << endl;

cout << "  " << endl;
cout << "***************  TOF1  ***************" << endl;
cout << "  " << endl;

int TDC_TOF1U_count[12] = {0};
int TDC_TOF1D_count[12] = {0};
int TDC_TOF1U_count_peak[12] = {0};
int TDC_TOF1D_count_peak[12] = {0};

if(flag!=0) nentries=flag;

for (Int_t i=0; i<nentries; i++) {
//fChain_TOF1->GetEntry(i);
fChain->GetEntry(i);

if(QUIET==0){
	if(nentries<=30000){
		if(i%1000==1) cout<<"**** "<<i<<" events done"<<endl;
	}
	if(nentries>30000){
		if(i%10000==1) cout<<"**** "<<i<<" events done"<<endl;
	}
}

	/// Fill ADC TOF1 histograms with points within time window 


	for (Int_t j=0; j<12; j++) {
		if (TDC_tof1U[j] > TDC_TOF1_min && TDC_tof1U[j] < TDC_TOF1_max) {
			h_ADC_tof1[j]->Fill(ADC_tof1U[j]);
			TDC_TOF1U_count_peak[j]++;
		}
	}

	for (Int_t j=0; j<12; j++) {
		if (TDC_tof1D[j] > TDC_TOF1_min && TDC_tof1D[j] < TDC_TOF1_max) {
			h_ADC_tof1[j+12]->Fill(ADC_tof1D[j]);
			TDC_TOF1D_count_peak[j]++;
		}
	}

	for (Int_t j=0; j<12; j++) {
		h_TDC_tof1U[j]->Fill(TDC_tof1U[j]);
		h_TDC_tof1D[j]->Fill(TDC_tof1D[j]);
	}

	for (Int_t j=0; j<12; j++) {
		if (TDC_tof1U[j] >= 0) TDC_TOF1U_count[j]++;
		if (TDC_tof1D[j] >= 0) TDC_TOF1D_count[j]++;
	}

}


if(QUIET==1) cout << "OK !" << endl;

cout << "  " << endl;
cout << "  " << endl;
cout << "***************  TOF2  ***************" << endl;
cout << "  " << endl;

int TDC_TOF2AO_count[12] = {0};
int TDC_TOF2BO_count[12] = {0};
int TDC_TOF2AI_count[12] = {0};
int TDC_TOF2BI_count[12] = {0};

int TDC_TOF2AO_count_peak[12] = {0};
int TDC_TOF2BO_count_peak[12] = {0};
int TDC_TOF2AI_count_peak[12] = {0};
int TDC_TOF2BI_count_peak[12] = {0};

for (Int_t i=0; i<nentries; i++) {
//fChain_TOF2->GetEntry(i);
fChain->GetEntry(i);

if(QUIET==0){
	if(nentries<=30000){
		if(i%1000==1) cout<<"**** "<<i<<" events done"<<endl;
	}
	if(nentries>30000){
		if(i%10000==1) cout<<"**** "<<i<<" events done"<<endl;
	}
}	

	/// Fill ADC TOF2 histograms with points within time window 

	for (Int_t j=0; j<12; j++) {
		if (TDC_tof2AO[j] > TDC_TOF2_min && TDC_tof2AO[j] < TDC_TOF2_max) {h_ADC_tof2[j]->Fill(ADC_tof2AO[j]);}
	}

	for (Int_t j=0; j<12; j++) {
		if (TDC_tof2BO[j] > TDC_TOF2_min && TDC_tof2BO[j] < TDC_TOF2_max) h_ADC_tof2[j+12]->Fill(ADC_tof2BO[j]);
	}

	for (Int_t j=0; j<12; j++) {
		if (TDC_tof2AI[j] > TDC_TOF2_min && TDC_tof2AI[j] < TDC_TOF2_max) h_ADC_tof2[j+24]->Fill(ADC_tof2AI[j]);
	}

	for (Int_t j=0; j<12; j++) {
		if (TDC_tof2BI[j] > TDC_TOF2_min && TDC_tof2BI[j] < TDC_TOF2_max) h_ADC_tof2[j+36]->Fill(ADC_tof2BI[j]);
	}

	//if (TDC_tof2BO[6] > TDC_TOF2_min && TDC_tof2BO[6] < TDC_TOF2_max) h_ADC_tof2[55]->Fill(ADC_tof2[55]);

	/*
	for (Int_t j=0; j<56; j++) {
		h_ADC_tof2[j]->Fill(ADC_tof2[j]);
	}
	*/

	for (Int_t j=0; j<12; j++) {
	h_TDC_tof2AO[j]->Fill(TDC_tof2AO[j]);
	h_TDC_tof2AI[j]->Fill(TDC_tof2AI[j]);
	h_TDC_tof2BO[j]->Fill(TDC_tof2BO[j]);
	h_TDC_tof2BI[j]->Fill(TDC_tof2BI[j]);
	}

	for (Int_t j=0; j<12; j++) {
		if (TDC_tof2AO[j] >= 0) TDC_TOF2AO_count[j]++;
		if (TDC_tof2BO[j] >= 0) TDC_TOF2BO_count[j]++;
		if (TDC_tof2AI[j] >= 0) TDC_TOF2AI_count[j]++;
		if (TDC_tof2BI[j] >= 0) TDC_TOF2BI_count[j]++;
	}

	//// Fill peak histograms

	for (Int_t j=0; j<12; j++) {
		if (TDC_tof2AO[j] > TDC_TOF2_min && TDC_tof2AO[j] < TDC_TOF2_max) {
			h_TDC_tof2AO_peak[j]->Fill(TDC_tof2AO[j]);
			TDC_TOF2AO_count_peak[j]++;
		}
	}

	for (Int_t j=0; j<12; j++) {
		if (TDC_tof2BO[j] > TDC_TOF2_min && TDC_tof2BO[j] < TDC_TOF2_max) {
			h_TDC_tof2BO_peak[j]->Fill(TDC_tof2BO[j]);
			TDC_TOF2BO_count_peak[j]++;
		}
	}

	for (Int_t j=0; j<12; j++) {
		if (TDC_tof2AI[j] > TDC_TOF2_min && TDC_tof2AI[j] < TDC_TOF2_max) {
			h_TDC_tof2AI_peak[j]->Fill(TDC_tof2AI[j]);
			TDC_TOF2AI_count_peak[j]++;
		} 
	}

	for (Int_t j=0; j<12; j++) {
		if (TDC_tof2BI[j] > TDC_TOF2_min && TDC_tof2BI[j] < TDC_TOF2_max) {
			h_TDC_tof2BI_peak[j]->Fill(TDC_tof2BI[j]);
			TDC_TOF2BI_count_peak[j]++;
		}
	}
}

char TOF1U_Str[12][100];
char TOF1D_Str[12][100];

char TOF1Upeak_Str[12][100];
char TOF1Dpeak_Str[12][100];

char AO_Str[12][100];
char BO_Str[12][100];
char AI_Str[12][100];
char BI_Str[12][100];

char AOpeak_Str[12][100];
char BOpeak_Str[12][100];
char AIpeak_Str[12][100];
char BIpeak_Str[12][100];

for (int q=0; q<12; q++) {
	sprintf(TOF1U_Str[q], "%i", TDC_TOF1U_count[q]);
	sprintf(TOF1Upeak_Str[q],"%i",TDC_TOF1U_count_peak[q]);

	sprintf(TOF1D_Str[q], "%i", TDC_TOF1D_count[q]);
	sprintf(TOF1Dpeak_Str[q],"%i",TDC_TOF1D_count_peak[q]);

	////

	sprintf(AO_Str[q], "%i", TDC_TOF2AO_count[q]);
	sprintf(AOpeak_Str[q],"%i",TDC_TOF2AO_count_peak[q]);

	sprintf(BO_Str[q],"%i",TDC_TOF2BO_count[q]);
	sprintf(BOpeak_Str[q],"%i",TDC_TOF2BO_count_peak[q]);

	sprintf(AI_Str[q],"%i",TDC_TOF2AI_count[q]);
	sprintf(AIpeak_Str[q],"%i",TDC_TOF2AI_count_peak[q]);

	sprintf(BI_Str[q],"%i",TDC_TOF2BI_count[q]);
	sprintf(BIpeak_Str[q],"%i",TDC_TOF2BI_count_peak[q]);
}

int TOF1UMax[12] = {0};
int TOF1DMax[12] = {0};

int AOMax[12] = {0};
int BOMax[12] = {0};
int AIMax[12] = {0};
int BIMax[12] = {0};

float TOF1UMaxScl[12] = {0};
float TOF1UMaxScl2[12] = {0};
float TOF1DMaxScl[12] = {0};
float TOF1DMaxScl2[12] = {0};

float AOMaxScl[12] = {0};
float BOMaxScl[12] = {0};
float AIMaxScl[12] = {0};
float BIMaxScl[12] = {0};

float AOMaxScl2[12] = {0};
float BOMaxScl2[12] = {0};
float AIMaxScl2[12] = {0};
float BIMaxScl2[12] = {0};

for (Int_t j=0; j<12; j++) {
	TOF1UMax[j] = h_TDC_tof1U[j]->GetMaximum();
	TOF1DMax[j] = h_TDC_tof1D[j]->GetMaximum();

	AOMax[j] = h_TDC_tof2AO[j]->GetMaximum();
	BOMax[j] = h_TDC_tof2BO[j]->GetMaximum();
	AIMax[j] = h_TDC_tof2AI[j]->GetMaximum();
	BIMax[j] = h_TDC_tof2BI[j]->GetMaximum();
}

float TOF1_height_max = 0;
float TOF2_height_max = 0;

if (normalize != 0) {
	for (Int_t j=0; j<12; j++) {
		if (TOF1UMax[j] > TOF1_height_max) TOF1_height_max = TOF1UMax[j];
	}
	for (Int_t j=0; j<12; j++) {
		if (TOF1DMax[j] > TOF1_height_max) TOF1_height_max = TOF1DMax[j];
	}

//// 

	for (Int_t j=0; j<12; j++) {
		if (AOMax[j] > TOF2_height_max) TOF2_height_max = AOMax[j];
	}
	for (Int_t j=0; j<12; j++) {
		if (BOMax[j] > TOF2_height_max) TOF2_height_max = BOMax[j];
	}
	for (Int_t j=0; j<12; j++) {
		if (AIMax[j] > TOF2_height_max) TOF2_height_max = AIMax[j];
	}
	for (Int_t j=0; j<12; j++) {
		if (BIMax[j] > TOF2_height_max) TOF2_height_max = BIMax[j];
	}

	TOF1_height_max = TOF1_height_max * 1.1;
	TOF2_height_max = TOF2_height_max * 1.1;

//// Set Uniform Height

	for (Int_t j=0; j<12; j++) {
		TOF1UMaxScl[j] = TOF1_height_max*0.4;
		TOF1DMaxScl[j] = TOF1_height_max*0.4;
		TOF1UMaxScl2[j] = TOF1_height_max*0.6;
		TOF1DMaxScl2[j] = TOF1_height_max*0.6;

		AOMaxScl[j] = TOF2_height_max*0.4;
		BOMaxScl[j] = TOF2_height_max*0.4;
		AIMaxScl[j] = TOF2_height_max*0.4;
		BIMaxScl[j] = TOF2_height_max*0.4;

		AOMaxScl2[j] = TOF2_height_max*0.6;
		BOMaxScl2[j] = TOF2_height_max*0.6;
		AIMaxScl2[j] = TOF2_height_max*0.6;
		BIMaxScl2[j] = TOF2_height_max*0.6;
	}

}

if (normalize == 0) {
	for (Int_t j=0; j<12; j++) {
		TOF1UMaxScl[j] = TOF1UMax[j]*0.4;
		TOF1DMaxScl[j] = TOF1DMax[j]*0.4;
		TOF1UMaxScl2[j] = TOF1UMax[j]*0.6;
		TOF1DMaxScl2[j] = TOF1DMax[j]*0.6;

		AOMaxScl[j] = AOMax[j]*0.4;
		BOMaxScl[j] = BOMax[j]*0.4;
		AIMaxScl[j] = AIMax[j]*0.4;
		BIMaxScl[j] = BIMax[j]*0.4;

		AOMaxScl2[j] = AOMax[j]*0.6;
		BOMaxScl2[j] = BOMax[j]*0.6;
		AIMaxScl2[j] = AIMax[j]*0.6;
		BIMaxScl2[j] = BIMax[j]*0.6;
	}
}

char Name_Can_ADC_tof1B[100];			char Title_Can_ADC_tof1B[100];

sprintf(Name_Can_ADC_tof1B,"TOF1 ADCs & TDCs -- Run %d  (Gap 1 - 12)  |  1U, 1D, 2U, 2D, etc...",run_number);
sprintf(Title_Can_ADC_tof1B,"TOF1 ADCs & TDCs -- Run %d  (Gap 1 - 12)  |  1U, 1D, 2U, 2D, etc...",run_number);

TCanvas *c7;
	c7 = new TCanvas(Name_Can_ADC_tof1B,Title_Can_ADC_tof1B,1200,500); 
	c7->Divide(8,6);
	c7->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=0; ican<12; ican++){
		c7->cd(TOF1_Mapping1[ican])->SetLogy();

		char TOF1Title1[250];
		sprintf(TOF1Title1, "Raw ADC (Ch. %d) - TOF1 | Up-%d : Run %d", ican, ican+1, run_number);
//		h_ADC_tof1[ican]->SetAxisRange(0, TOF1_Graph_xmax,"X");
//		if(TOF1_Graph_ymax > 0) {  
//			h_ADC_tof1[ican]->SetAxisRange(0, TOF1_Graph_ymax,"Y");
//		}
		h_ADC_tof1[ican]->SetTitle(TOF1Title1);
		h_ADC_tof1[ican]->Draw();
	}

	for(Int_t ican=0; ican<12; ican++){
		c7->cd(TOF1_Mapping2[ican])->SetLogy();

		char TOF1Title2[250];
		sprintf(TOF1Title2, "Raw ADC (Ch. %d) - TOF1 | Down-%d : Run %d", ican+12, ican+1, run_number);
//		h_ADC_tof1[ican]->SetAxisRange(0, TOF1_Graph_xmax,"X");
//		if(TOF1_Graph_ymax > 0) {  
//			h_ADC_tof1[ican]->SetAxisRange(0, TOF1_Graph_ymax,"Y");
//		}
		h_ADC_tof1[ican+12]->SetTitle(TOF1Title2);
		h_ADC_tof1[ican+12]->Draw();
	}

	for(Int_t ican=0; ican<12; ican++){
		c7->cd(TOF1_Mapping1[ican]+24);

		char TOF1Title3[250];
		sprintf(TOF1Title3, "Raw TDC (Ch. %d) - TOF1 | Up-%d : Run %d", ican, ican+1, run_number);
//		h_ADC_tof1[ican]->SetAxisRange(0, TOF1_Graph_xmax,"X");

		if(normalize != 0) {  
			h_TDC_tof1U[ican]->SetAxisRange(0, TOF1_height_max,"Y");
		}

		h_TDC_tof1U[ican]->SetTitle(TOF1Title3);
		h_TDC_tof1U[ican]->Draw();
		TLatex *tex17 = new TLatex(2500,TOF1UMaxScl[ican],TOF1U_Str[ican]);
		tex17->SetTextSize(0.125);
		tex17->Draw("same");
		TLatex *tex18 = new TLatex(2500,TOF1UMaxScl2[ican],TOF1Upeak_Str[ican]);
		tex18->SetTextSize(0.125);
		tex18->Draw("same");
	}

	for(Int_t ican=0; ican<12; ican++){
		c7->cd(TOF1_Mapping2[ican]+24);

		char TOF1Title4[250];
		sprintf(TOF1Title4, "Raw TDC (Ch. %d) - TOF1 | Down-%d : Run %d", ican+12, ican+1, run_number);
//		h_ADC_tof1[ican]->SetAxisRange(0, TOF1_Graph_xmax,"X");

		if(normalize != 0) {  
			h_TDC_tof1D[ican]->SetAxisRange(0, TOF1_height_max,"Y");
		}

		h_TDC_tof1D[ican]->SetTitle(TOF1Title4);
		h_TDC_tof1D[ican]->Draw();
		TLatex *tex19 = new TLatex(2500,TOF1DMaxScl[ican],TOF1D_Str[ican]);
		tex19->SetTextSize(0.125);
		tex19->Draw("same");
		TLatex *tex20 = new TLatex(2500,TOF1DMaxScl2[ican],TOF1Dpeak_Str[ican]);
		tex20->SetTextSize(0.125);
		tex20->Draw("same");
	}

///////////////////////////////////////

char Name_Can_ADC_tof2B[100];			char Title_Can_ADC_tof2B[100];
char Name_Can_ADC_tof2C[100];			char Title_Can_ADC_tof2C[100];

sprintf(Name_Can_ADC_tof2B,"TOF2 ADCs & TDCs -- Run %d  (Gap 1 - 6)  |  1AO, 1AI, 1BO, 1BI, etc...",run_number);
sprintf(Title_Can_ADC_tof2B,"TOF2 ADCs & TDCs -- Run %d  (Gap 1 - 6)  |  1AO, 1AI, 1BO, 1BI, etc...",run_number);
sprintf(Name_Can_ADC_tof2C,"TOF2 ADCs & TDCs -- Run %d  (Gap 7 - 12)  |  7AO, 7AI, 7BO, 7BI, etc...",run_number);
sprintf(Title_Can_ADC_tof2C,"TOF2 ADCs & TDCs -- Run %d  (Gap 7 - 12)  |  7AO, 7AI, 7BO, 7BI, etc...",run_number);

TCanvas *c8;
	c8 = new TCanvas(Name_Can_ADC_tof2B,Title_Can_ADC_tof2B,1200,500); 
	c8->Divide(8,6);
	c8->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=0; ican<6; ican++){
		c8->cd(TOF2_Mapping1[ican])->SetLogy();

		char TOF2Title1[250];
		sprintf(TOF2Title1, "Raw ADC (Ch. %d) - TOF2 | OutA-%d : Run %d", ican, ican+1, run_number);
//		h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_xmax,"X");
//		if(TOF2_Graph_ymax > 0) {  
//			h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_ymax,"Y");
//		}
		h_ADC_tof2[ican]->SetTitle(TOF2Title1);
		h_ADC_tof2[ican]->Draw();
	}

	for(Int_t ican=0; ican<6; ican++){
		c8->cd(TOF2_Mapping2[ican])->SetLogy();
		
		char TOF2Title2[250];
		sprintf(TOF2Title2, "Raw ADC (Ch. %d) - TOF2 | InA-%d : Run %d", ican+24, ican+1, run_number);
//		h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_xmax,"X");
//		if(TOF2_Graph_ymax > 0) {  
//			h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_ymax,"Y");
//		}
		h_ADC_tof2[ican+24]->SetTitle(TOF2Title2);
		h_ADC_tof2[ican+24]->Draw();
	}
	
	for(Int_t ican=0; ican<6; ican++){
		c8->cd(TOF2_Mapping3[ican])->SetLogy();

		char TOF2Title3[250];
		sprintf(TOF2Title3, "Raw ADC (Ch. %d) - TOF2 | OutB-%d : Run %d", ican+12, ican+1, run_number);
//		h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_xmax,"X");
//		if(TOF2_Graph_ymax > 0) {  
//			h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_ymax,"Y");
//		}
		if (ican+12 == 18) {
			sprintf(TOF2Title3, "Raw ADC (Ch. 55) - TOF2 | OutB-7 : Run %d", run_number);
			h_ADC_tof2[55]->SetTitle(TOF2Title3);
			h_ADC_tof2[55]->Draw();
		}
		else {
			h_ADC_tof2[ican+12]->SetTitle(TOF2Title3);
			h_ADC_tof2[ican+12]->Draw();
		}
	}

	for(Int_t ican=0; ican<6; ican++){
		c8->cd(TOF2_Mapping4[ican])->SetLogy();

		char TOF2Title4[250];
		sprintf(TOF2Title4, "Raw ADC (Ch. %d) - TOF2 | InB-%d : Run %d", ican+36, ican+1, run_number);
//		h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_xmax,"X");
//		if(TOF2_Graph_ymax > 0) {  
//			h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_ymax,"Y");
//		}
		h_ADC_tof2[ican+36]->SetTitle(TOF2Title4);
		h_ADC_tof2[ican+36]->Draw();
	}

	for(Int_t ican=0; ican<6; ican++){
		c8->cd(TOF2_Mapping1[ican+6]);

		char TOF2Title5[250];
		sprintf(TOF2Title5, "Raw TDC (Ch. %d) - TOF2 | OutA-%d : Run %d", ican, ican+1, run_number);

		if(normalize != 0) {  
			h_TDC_tof2AO[ican]->SetAxisRange(0, TOF2_height_max,"Y");
		}

		h_TDC_tof2AO[ican]->SetTitle(TOF2Title5);
		h_TDC_tof2AO[ican]->Draw();
		//h_TDC_tof2AO_peak[ican]->SetLineColor(3);
		//h_TDC_tof2AO_peak[ican]->Draw("same");
		TLatex *tex1 = new TLatex(2500,AOMaxScl[ican],AO_Str[ican]);
		tex1->SetTextSize(0.125);
		tex1->Draw("same");
		TLatex *tex2 = new TLatex(2500,AOMaxScl2[ican],AOpeak_Str[ican]);
		tex2->SetTextSize(0.125);
		tex2->Draw("same");
	}

	for(Int_t ican=0; ican<6; ican++){
		c8->cd(TOF2_Mapping2[ican+6]);
		
		char TOF2Title6[250];
		sprintf(TOF2Title6, "Raw TDC (Ch. %d) - TOF2 | InA-%d : Run %d", ican+24, ican+1, run_number);

		if(normalize != 0) {  
			h_TDC_tof2AI[ican]->SetAxisRange(0, TOF2_height_max,"Y");
		}

		h_TDC_tof2AI[ican]->SetTitle(TOF2Title6);
		h_TDC_tof2AI[ican]->Draw();
		//h_TDC_tof2AI_peak[ican]->SetLineColor(3);
		//h_TDC_tof2AI_peak[ican]->Draw("same");
		TLatex *tex3 = new TLatex(2500,AIMaxScl[ican],AI_Str[ican]);
		tex3->SetTextSize(0.125);
		tex3->Draw("same");
		TLatex *tex4 = new TLatex(2500,AIMaxScl2[ican],AIpeak_Str[ican]);
		tex4->SetTextSize(0.125);
		tex4->Draw("same");
	}
	
	for(Int_t ican=0; ican<6; ican++){
		c8->cd(TOF2_Mapping3[ican+6]);

		char TOF2Title7[250];
		sprintf(TOF2Title7, "Raw TDC (Ch. %d) - TOF2 | OutB-%d : Run %d", ican+12, ican+1, run_number);

		if(normalize != 0) {  
			h_TDC_tof2BO[ican]->SetAxisRange(0, TOF2_height_max,"Y");
		}

		h_TDC_tof2BO[ican]->SetTitle(TOF2Title7);
		h_TDC_tof2BO[ican]->Draw();
		//h_TDC_tof2BO_peak[ican]->SetLineColor(3);
		//h_TDC_tof2BO_peak[ican]->Draw("same");
		TLatex *tex5 = new TLatex(2500,BOMaxScl[ican],BO_Str[ican]);
		tex5->SetTextSize(0.125);
		tex5->Draw("same");
		TLatex *tex6 = new TLatex(2500,BOMaxScl2[ican],BOpeak_Str[ican]);
		tex6->SetTextSize(0.125);
		tex6->Draw("same");
	}

	for(Int_t ican=0; ican<6; ican++){
		c8->cd(TOF2_Mapping4[ican+6]);

		char TOF2Title8[250];
		sprintf(TOF2Title8, "Raw TDC (Ch. %d) - TOF2 | InB-%d : Run %d", ican+36, ican+1, run_number);

		if(normalize != 0) {  
			h_TDC_tof2BI[ican]->SetAxisRange(0, TOF2_height_max,"Y");
		}

		h_TDC_tof2BI[ican]->SetTitle(TOF2Title8);
		h_TDC_tof2BI[ican]->Draw();
		//h_TDC_tof2BI_peak[ican]->SetLineColor(3);
		//h_TDC_tof2BI_peak[ican]->Draw("same");
		TLatex *tex7 = new TLatex(2500,BIMaxScl[ican],BI_Str[ican]);
		tex7->SetTextSize(0.125);
		tex7->Draw("same");
		TLatex *tex8 = new TLatex(2500,BIMaxScl2[ican],BIpeak_Str[ican]);
		tex8->SetTextSize(0.125);
		tex8->Draw("same");
	}

	TCanvas *c12;
	c12 = new TCanvas(Name_Can_ADC_tof2C,Title_Can_ADC_tof2C,1200,500); 
	c12->Divide(8,6);
	c12->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=6; ican<12; ican++){
		c12->cd(TOF2_Mapping1[ican-6])->SetLogy();

		char TOF2Title1[250];
		sprintf(TOF2Title1, "Raw ADC (Ch. %d) - TOF2 | OutA-%d : Run %d", ican, ican+1, run_number);
//		h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_xmax,"X");
//		if(TOF2_Graph_ymax > 0) {  
//			h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_ymax,"Y");
//		}
		h_ADC_tof2[ican]->SetTitle(TOF2Title1);
		h_ADC_tof2[ican]->Draw();
	}

	for(Int_t ican=6; ican<12; ican++){
		c12->cd(TOF2_Mapping2[ican-6])->SetLogy();
		
		char TOF2Title2[250];
		sprintf(TOF2Title2, "Raw ADC (Ch. %d) - TOF2 | InA-%d : Run %d", ican+24, ican+1, run_number);
//		h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_xmax,"X");
//		if(TOF2_Graph_ymax > 0) {  
//			h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_ymax,"Y");
//		}
		h_ADC_tof2[ican+24]->SetTitle(TOF2Title2);
		h_ADC_tof2[ican+24]->Draw();
	}
	
	for(Int_t ican=6; ican<12; ican++){
		c12->cd(TOF2_Mapping3[ican-6])->SetLogy();

		char TOF2Title3[250];
		sprintf(TOF2Title3, "Raw ADC (Ch. %d) - TOF2 | OutB-%d : Run %d", ican+12, ican+1, run_number);
//		h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_xmax,"X");
//		if(TOF2_Graph_ymax > 0) {  
//			h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_ymax,"Y");
//		}
		//if (ican+12 == 18) {
		//	sprintf(TOF2Title3, "Raw ADC (Ch. 55) - TOF2 | OutB-7 : Run %d", run_number);
		//	h_ADC_tof2[55]->SetTitle(TOF2Title3);
		//	h_ADC_tof2[55]->Draw();
		//}
		//else {
		
			h_ADC_tof2[ican+12]->SetTitle(TOF2Title3);
			h_ADC_tof2[ican+12]->Draw();
		//}
	}

	for(Int_t ican=6; ican<12; ican++){
		c12->cd(TOF2_Mapping4[ican-6])->SetLogy();

		char TOF2Title4[250];
		sprintf(TOF2Title4, "Raw ADC (Ch. %d) - TOF2 | InB-%d : Run %d", ican+36, ican+1, run_number);
//		h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_xmax,"X");
//		if(TOF2_Graph_ymax > 0) {  
//			h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_ymax,"Y");
//		}
		h_ADC_tof2[ican+36]->SetTitle(TOF2Title4);
		h_ADC_tof2[ican+36]->Draw();
	}

	for(Int_t ican=0; ican<6; ican++){
		c12->cd(TOF2_Mapping1[ican+6]);

		char TOF2Title5[250];
		sprintf(TOF2Title5, "Raw TDC (Ch. %d) - TOF2 | OutA-%d : Run %d", ican, ican+7, run_number);

		if(normalize != 0) {  
			h_TDC_tof2AO[ican+6]->SetAxisRange(0, TOF2_height_max,"Y");
		}

		h_TDC_tof2AO[ican+6]->SetTitle(TOF2Title5);
		h_TDC_tof2AO[ican+6]->Draw();
		//h_TDC_tof2AO_peak[ican+6]->SetLineColor(3);
		//h_TDC_tof2AO_peak[ican+6]->Draw("same");
		TLatex *tex9 = new TLatex(2500,AOMaxScl[ican+6],AO_Str[ican+6]);
		tex9->SetTextSize(0.125);
		tex9->Draw("same");
		TLatex *tex10 = new TLatex(2500,AOMaxScl2[ican+6],AOpeak_Str[ican+6]);
		tex10->SetTextSize(0.125);
		tex10->Draw("same");
	}

	for(Int_t ican=0; ican<6; ican++){
		c12->cd(TOF2_Mapping2[ican+6]);
		
		char TOF2Title6[250];
		sprintf(TOF2Title6, "Raw TDC (Ch. %d) - TOF2 | InA-%d : Run %d", ican+24, ican+7, run_number);

		if(normalize != 0) {  
			h_TDC_tof2AI[ican+6]->SetAxisRange(0, TOF2_height_max,"Y");
		}

		h_TDC_tof2AI[ican+6]->SetTitle(TOF2Title6);
		h_TDC_tof2AI[ican+6]->Draw();
		//h_TDC_tof2AI_peak[ican+6]->SetLineColor(3);
		//h_TDC_tof2AI_peak[ican+6]->Draw("same");
		TLatex *tex11 = new TLatex(2500,AIMaxScl[ican+6],AI_Str[ican+6]);
		tex11->SetTextSize(0.125);
		tex11->Draw("same");
		TLatex *tex12 = new TLatex(2500,AIMaxScl2[ican+6],AIpeak_Str[ican+6]);
		tex12->SetTextSize(0.125);
		tex12->Draw("same");
	}
	
	for(Int_t ican=0; ican<6; ican++){
		c12->cd(TOF2_Mapping3[ican+6]);

		char TOF2Title7[250];
		sprintf(TOF2Title7, "Raw TDC (Ch. %d) - TOF2 | OutB-%d : Run %d", ican+12, ican+7, run_number);

		if(normalize != 0) {  
			h_TDC_tof2BO[ican+6]->SetAxisRange(0, TOF2_height_max,"Y");
		}

		h_TDC_tof2BO[ican+6]->SetTitle(TOF2Title7);
		h_TDC_tof2BO[ican+6]->Draw();
		//h_TDC_tof2BO_peak[ican+6]->SetLineColor(3);
		//h_TDC_tof2BO_peak[ican+6]->Draw("same");
		TLatex *tex13 = new TLatex(2500,BOMaxScl[ican+6],AO_Str[ican+6]);
		tex13->SetTextSize(0.125);
		tex13->Draw("same");
		TLatex *tex14 = new TLatex(2500,BOMaxScl2[ican+6],AOpeak_Str[ican+6]);
		tex14->SetTextSize(0.125);
		tex14->Draw("same");
	}

	for(Int_t ican=0; ican<6; ican++){
		c12->cd(TOF2_Mapping4[ican+6]);

		char TOF2Title8[250];
		sprintf(TOF2Title8, "Raw TDC (Ch. %d) - TOF2 | InB-%d : Run %d", ican+36, ican+7, run_number);

		if(normalize != 0) {  
			h_TDC_tof2BI[ican+6]->SetAxisRange(0, TOF2_height_max,"Y");
		}

		h_TDC_tof2BI[ican+6]->SetTitle(TOF2Title8);
		h_TDC_tof2BI[ican+6]->Draw();
		//h_TDC_tof2BI_peak[ican+6]->SetLineColor(3);
		//h_TDC_tof2BI_peak[ican+6]->Draw("same");
		TLatex *tex15 = new TLatex(2500,BIMaxScl[ican+6],BI_Str[ican+6]);
		tex15->SetTextSize(0.125);
		tex15->Draw("same");
		TLatex *tex16 = new TLatex(2500,BIMaxScl2[ican+6],BIpeak_Str[ican+6]);
		tex16->SetTextSize(0.125);
		tex16->Draw("same");
	}		

return;

}

