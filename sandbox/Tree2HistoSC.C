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
  
void Tree2HistoSC(Int_t run_number=5, Int_t time_window = 0, Int_t n_MWPC_hits=1, Int_t flag=0, Int_t QUIET=0) { 

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

Int_t TDC_tof1U[12];  Int_t TDC_TOF1U[12];
Int_t TDC_tof1D[12];  Int_t TDC_TOF1D[12];

Int_t TDC_tof2AO[12];   Int_t TDC_TOF2AO[12];
Int_t TDC_tof2BO[12];   Int_t TDC_TOF2BO[12];
Int_t TDC_tof2AI[12];   Int_t TDC_TOF2AI[12];
Int_t TDC_tof2BI[12];   Int_t TDC_TOF2BI[12];

Int_t MWPCADC[512];

int good_event_counter = 0;	

char path_input[200];                   char file_mapping[200];
sprintf(path_input,"%s",path_merged);          sprintf(file_mapping,"../Mapping");
//sprintf(file_mapping,path_mapping);
//sprintf(path_input,"/triumfcs/trshare/trek/E36/Data/April_2015/root");

char Name_finput[200];
sprintf(Name_finput,"%s/Run%dMS.root",path_input, run_number);

char path_output[200];
sprintf(path_output,"%s",path_histo);

TFile *foutput;
char Name_foutput[200];
sprintf(Name_foutput,"%s/Hist_Run%dMSC.root",path_output, run_number);
foutput = new TFile(Name_foutput,"RECREATE");

char par_finput[200];
//sprintf(par_finput,"/triumfcs/trshare/trek/E36/Macros/April_2015/Mapping/%s",source_mapping);
sprintf(par_finput,"%s/SFT_Mapping_Oct14.txt",file_mapping);

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
cout << "   " << endl;

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

//char footer[100];
//sprintf(footer,"Event_Display_MS.C  --  Run %d ; Event %d",Run_Number,ievt);

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

fChain->SetBranchAddress("MWPCADC",MWPCADC);

fChain->SetBranchAddress("TDC_TOF1U",TDC_tof1U);
fChain->SetBranchAddress("TDC_TOF1D",TDC_tof1D);
fChain->SetBranchAddress("TDC_TOF2AO",TDC_tof2AO);
fChain->SetBranchAddress("TDC_TOF2AI",TDC_tof2AI);
fChain->SetBranchAddress("TDC_TOF2BO",TDC_tof2BO);
fChain->SetBranchAddress("TDC_TOF2BI",TDC_tof2BI);			

TH1D *h_ADC_High_TARGET[256];	char Title_ADC_High_TARGET[256][100];	char Name_ADC_High_TARGET[256][100];
TH1D *h_ADC_Low_TARGET[256];    char Title_ADC_Low_TARGET[256][100];	char Name_ADC_Low_TARGET[256][100];
TH1D *h_TDC_LE_TARGET[256];     char Title_TDC_LE_TARGET[256][100];	char Name_TDC_LE_TARGET[256][100];     
TH1D *h_TDC_TE_TARGET[256];	char Title_TDC_TE_TARGET[256][100];	char Name_TDC_TE_TARGET[256][100]; 

TH1D *h_ADC_High_SFT[256];	char Title_ADC_High_SFT[256][100];	char Name_ADC_High_SFT[256][100];
TH1D *h_ADC_Low_SFT[256];    	char Title_ADC_Low_SFT[256][100];	char Name_ADC_Low_SFT[256][100];
TH1D *h_TDC_LE_SFT[256];     	char Title_TDC_LE_SFT[256][100];	char Name_TDC_LE_SFT[256][100];     
TH1D *h_TDC_TE_SFT[256];	char Title_TDC_TE_SFT[256][100];	char Name_TDC_TE_SFT[256][100]; 

TH1D *h_ADC_tof1[24];   char Title_ADC_tof1[24][100];	char Name_ADC_tof1[24][100];     
TH1D *h_ADC_tof2[56];	char Title_ADC_tof2[56][100];	char Name_ADC_tof2[56][100];

TH1D *h_TDC_tof1U[12];   char Title_TDC_tof1U[12][100];	char Name_TDC_tof1U[12][100];
TH1D *h_TDC_tof1D[12];   char Title_TDC_tof1D[12][100];	char Name_TDC_tof1D[12][100];
TH1D *h_TDC_tof2AO[12];   char Title_TDC_tof2AO[12][100];	char Name_TDC_tof2AO[12][100];
TH1D *h_TDC_tof2AI[12];   char Title_TDC_tof2AI[12][100];	char Name_TDC_tof2AI[12][100];
TH1D *h_TDC_tof2BO[12];   char Title_TDC_tof2BO[12][100];	char Name_TDC_tof2BO[12][100];
TH1D *h_TDC_tof2BI[12];   char Title_TDC_tof2BI[12][100];	char Name_TDC_tof2BI[12][100];

TH1D *h_MWPCADC[512];	char Title_MWPCADC[512][100];	char Name_MWPCADC[512][100];  

for(int i=0; i<256; i++){
sprintf(Title_ADC_High_TARGET[i],"Raw ADC High Gain (Ch. %d)  --  TARGET",i); 
sprintf(Title_ADC_Low_TARGET[i],"Raw ADC Low Gain (Ch. %d)  --  TARGET",i); 
sprintf(Title_TDC_LE_TARGET[i],"Raw TDC (LE) (Ch. %d)  --  TARGET",i); 
sprintf(Title_TDC_TE_TARGET[i],"Raw TDC (TE) (Ch. %d)  --  TARGET",i); 

sprintf(Name_ADC_High_TARGET[i],"ADC_High (Ch. %d) - TARGET",i);
sprintf(Name_ADC_Low_TARGET[i],"ADC_Low (Ch. %d) - TARGET",i);
sprintf(Name_TDC_LE_TARGET[i],"TDC_LE (Ch. %d) - TARGET",i);
sprintf(Name_TDC_TE_TARGET[i],"TDC_TE (Ch. %d) - TARGET",i);
}

for(int i=0; i<128; i++){
sprintf(Title_ADC_High_SFT[i],"Raw ADC High Gain (Ch. %d)  --  SFT",i); 
sprintf(Title_ADC_Low_SFT[i],"Raw ADC Low Gain (Ch. %d)  --  SFT",i); 
sprintf(Title_TDC_LE_SFT[i],"Raw TDC (LE) (Ch. %d)  --  SFT",i); 
sprintf(Title_TDC_TE_SFT[i],"Raw TDC (TE) (Ch. %d)  --  SFT",i); 

sprintf(Name_ADC_High_SFT[i],"ADC_High (Ch. %d) - SFT",i);
sprintf(Name_ADC_Low_SFT[i],"ADC_Low (Ch. %d) - SFT",i);
sprintf(Name_TDC_LE_SFT[i],"TDC_LE (Ch. %d) - SFT",i);
sprintf(Name_TDC_TE_SFT[i],"TDC_TE (Ch. %d) - SFT",i);
}

for (int i=0; i<24; i++) {
	sprintf(Title_ADC_tof1[i],"Raw ADC (Ch. %d)  --  TOF1",i);
	sprintf(Name_ADC_tof1[i],"ADC_TOF1 (Ch. %d)",i);
}

for (int i=0; i<56; i++) {
	sprintf(Title_ADC_tof2[i],"Raw ADC (Ch. %d)  --  TOF2",i);
	sprintf(Name_ADC_tof2[i],"ADC_TOF2 (Ch. %d)",i);
}

for (int i=0; i<512; i++) {
	sprintf(Title_MWPCADC[i],"Raw ADC (Ch. %d)  --  MWPC",i);
	sprintf(Name_MWPCADC[i],"ADC_MWPC (Ch. %d)",i);
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
}

//create histograms

for (int i=0; i<256; i++){
h_ADC_High_TARGET[i] = new TH1D(Name_ADC_High_TARGET[i],Title_ADC_High_TARGET[i],1000,0,4100);
h_ADC_Low_TARGET[i] = new TH1D(Name_ADC_Low_TARGET[i],Title_ADC_Low_TARGET[i],1000,0,4100);
h_TDC_LE_TARGET[i] = new TH1D(Name_TDC_LE_TARGET[i],Title_TDC_LE_TARGET[i],1000,0,4100);
h_TDC_TE_TARGET[i] = new TH1D(Name_TDC_TE_TARGET[i],Title_TDC_TE_TARGET[i],1000,0,4100);
}

for (int i=0; i<128; i++){
h_ADC_High_SFT[i] = new TH1D(Name_ADC_High_SFT[i],Title_ADC_High_SFT[i],1000,0,4100);
h_ADC_Low_SFT[i] = new TH1D(Name_ADC_Low_SFT[i],Title_ADC_Low_SFT[i],1000,0,4100);
h_TDC_LE_SFT[i] = new TH1D(Name_TDC_LE_SFT[i],Title_TDC_LE_SFT[i],1000,0,4100);
h_TDC_TE_SFT[i] = new TH1D(Name_TDC_TE_SFT[i],Title_TDC_TE_SFT[i],1000,0,4100);
}

for (int  i=0; i<24; i++) {
h_ADC_tof1[i] = new TH1D(Name_ADC_tof1[i],Title_ADC_tof1[i],1000,0,4100);
}

for (int  i=0; i<56; i++) {
h_ADC_tof2[i] = new TH1D(Name_ADC_tof2[i],Title_ADC_tof2[i],1000,0,4100);
}

for (int  i=0; i<512; i++) {
h_MWPCADC[i] = new TH1D(Name_MWPCADC[i],Title_MWPCADC[i],1000,0,4100);
}

for (int  i=0; i<12; i++) {
h_TDC_tof1U[i] = new TH1D(Name_TDC_tof1U[i],Title_TDC_tof1U[i],1000,0,4100);
h_TDC_tof1D[i] = new TH1D(Name_TDC_tof1D[i],Title_TDC_tof1D[i],1000,0,4100);
h_TDC_tof2AO[i] = new TH1D(Name_TDC_tof2AO[i],Title_TDC_tof2AO[i],1000,0,4100);
h_TDC_tof2AI[i] = new TH1D(Name_TDC_tof2AI[i],Title_TDC_tof2AI[i],1000,0,4100);
h_TDC_tof2BO[i] = new TH1D(Name_TDC_tof2BO[i],Title_TDC_tof2BO[i],1000,0,4100);
h_TDC_tof2BI[i] = new TH1D(Name_TDC_tof2BI[i],Title_TDC_tof2BI[i],1000,0,4100);
}

//read all entries and fill the histograms
Int_t nentries = (Int_t)fChain->GetEntries();
cout <<  "Total Number of Entries :     " <<  nentries << endl;

if(flag!=0) nentries=flag;
for (Int_t i=0; i<nentries; i++) {
 //for (Int_t i=0; i<100000; i++) {
	fChain->GetEntry(i);

	if(QUIET==0){
		if(nentries<=30000){
			if(i%1000==1) cout<<"**** "<<i<<" events done"<<endl;
		}
		if(nentries>30000){
			if(i%10000==1) cout<<"**** "<<i<<" events done"<<endl;
		}
	}

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
		ADC_tof1[i] = ADC_TOF1U[i];
		ADC_tof1[i+12] = ADC_TOF1D[i];
		ADC_tof2[i] = ADC_TOF2AO[i];
		ADC_tof2[i+12] = ADC_TOF2BO[i];
		ADC_tof2[i+24] = ADC_TOF2AI[i];
		ADC_tof2[i+36] = ADC_TOF2BI[i];
	}

	//for (Int_t j_TOF1=0; j_TOF1<24; j_TOF1++) {
//		cout << j_TOF1 << " -- " << ADC_TOF1[j_TOF1] << endl;
	//	ADC_tof1[j_TOF1] = ADC_tof1[j_TOF1]-par_temp_TOF1[1][j_TOF1];
	//}

	//for (Int_t j_TOF2=0; j_TOF2<56; j_TOF2++) {
	//	ADC_tof2[j_TOF2] = ADC_tof2[j_TOF2]-par_temp_TOF2[1][j_TOF2];
	//}

	for (Int_t j_MWPC=0; j_MWPC<512; j_MWPC++) {
		MWPCADC[j_MWPC] = MWPCADC[j_MWPC]-MWPC_thr;
	}

	for (Int_t j_TDCTOF=0; j_TDCTOF<12; j_TDCTOF++) {
	    TDC_TOF1U[j_TDCTOF] = TDC_tof1U[j_TDCTOF];
	    TDC_TOF1D[j_TDCTOF] = TDC_tof1D[j_TDCTOF];
	    TDC_TOF2AO[j_TDCTOF] = TDC_tof2AO[j_TDCTOF];
	    TDC_TOF2BO[j_TDCTOF] = TDC_tof2BO[j_TDCTOF];
	    TDC_TOF2AI[j_TDCTOF] = TDC_tof2AI[j_TDCTOF];
	    TDC_TOF2BI[j_TDCTOF] = TDC_tof2BI[j_TDCTOF];
  	}

  	/// Determine good gap without time window

	if (time_window==0) {
		for (int kk=0; kk<12; kk++) {
			if (kk == 0) {
				if ((ADC_tof2[0]>0) || (ADC_tof2[24]>0) || (ADC_tof2[12]>0) || (ADC_tof2[36]>0)) {
					if((ADC_tof1[0] > 0) || (ADC_tof1[12] > 0) || (ADC_tof1[1] > 0) || (ADC_tof1[13] > 0) || (ADC_tof1[11] > 0) || (ADC_tof1[23] > 0)) has_data_TOF++;
				} 
			}

			/// Replace ADC TOF2[18] with ADC TOF2[55]
			else if (kk == 6) {
				if ((ADC_tof2[6]>0) || (ADC_tof2[30]>0) || (ADC_tof2[18]>0) || (ADC_tof2[42]>0)) {
					if((ADC_tof1[6] > 0) || (ADC_tof1[18] > 0) || (ADC_tof1[7] > 0) || (ADC_tof1[19] > 0) || (ADC_tof1[5] > 0) || (ADC_tof1[17] > 0)) has_data_TOF++;
				}
			}

			else if (kk == 11) {
				if ((ADC_tof2[11]>0) || (ADC_tof2[35]>0) || (ADC_tof2[23]>0) || (ADC_tof2[47]>0)) {
					if((ADC_tof1[11] > 0) || (ADC_tof1[23] > 0) || (ADC_tof1[0] > 0) || (ADC_tof1[12] > 0) || (ADC_tof1[10] > 0) || (ADC_tof1[22] > 0)) has_data_TOF++;
				}
			}
		
			else {
				if ((ADC_tof2[kk]>0) || (ADC_tof2[kk+24]>0) || (ADC_tof2[kk+12]>0) || (ADC_tof2[kk+36]>0))  {
					if((ADC_tof1[kk] > 0) || (ADC_tof1[kk+12] > 0) || (ADC_tof1[kk+1] > 0) || (ADC_tof1[kk+13] > 0) || (ADC_tof1[kk-1] > 0) || (ADC_tof1[kk+11] > 0)) has_data_TOF++;
				}
			}

		}
	}

	/// Determine good gap with time window

	bool has_ADC_TOF1_hit[12] = {false};
	bool has_TDC_TOF1_hit[12] = {false};
	bool has_ADC_TOF2_hit[12] = {false};
	bool has_TDC_TOF2_hit[12] = {false};

	if (time_window!=0) {
		for (int kk=0; kk<12; kk++) {
			
			/// Replace ADC TOF2[18] with ADC TOF2[55]
			//if (kk == 6) {
			//	if ((ADC_tof2[6]>0) || (ADC_tof2[30]>0) || (ADC_tof2[55]>0) || (ADC_tof2[42]>0)) {has_ADC_TOF2_hit[6]=true;}
			//	if (((TDC_TOF2AO[6]>TDC_TOF2_min && TDC_TOF2AO[6] < TDC_TOF2_max) || (TDC_TOF2AI[6]>TDC_TOF2_min && TDC_TOF2AI[6] < TDC_TOF2_max)) 
			//    || ((TDC_TOF2BO[6]>TDC_TOF2_min && TDC_TOF2BO[6] < TDC_TOF2_max) || (TDC_TOF2BI[6]>TDC_TOF2_min && TDC_TOF2BI[6] < TDC_TOF2_max))) {has_TDC_TOF2_hit[6]=true;}

			//	if (ADC_tof1[6]>0 || ADC_tof1[18]>0) {has_ADC_TOF1_hit[6] = true;}
			//  	if ((TDC_TOF1U[6]>TDC_TOF1_min && TDC_TOF1U[6]<TDC_TOF1_max) || (TDC_TOF1D[6]>TDC_TOF1_min && TDC_TOF1D[6]<TDC_TOF1_max)) {has_TDC_TOF1_hit[6] = true;}
		  	//}

			//else {
				if ((ADC_tof2[kk]>0) || (ADC_tof2[kk+24]>0) || (ADC_tof2[kk+12]>0) || (ADC_tof2[kk+36]>0)) {has_ADC_TOF2_hit[kk]=true;}
				if (((TDC_TOF2AO[kk]>TDC_TOF2_min && TDC_TOF2AO[kk] < TDC_TOF2_max) || (TDC_TOF2AI[kk]>TDC_TOF2_min && TDC_TOF2AI[kk] < TDC_TOF2_max)) 
			    || ((TDC_TOF2BO[kk]>TDC_TOF2_min && TDC_TOF2BO[kk] < TDC_TOF2_max) || (TDC_TOF2BI[kk]>TDC_TOF2_min && TDC_TOF2BI[kk] < TDC_TOF2_max))) {has_TDC_TOF2_hit[kk]=true;}

				if (ADC_tof1[kk]>0 || ADC_tof1[kk+12]>0) {has_ADC_TOF1_hit[kk] = true;}
			  	if ((TDC_TOF1U[kk]>TDC_TOF1_min && TDC_TOF1U[kk]<TDC_TOF1_max) || (TDC_TOF1D[kk]>TDC_TOF1_min && TDC_TOF1D[kk]<TDC_TOF1_max)) {has_TDC_TOF1_hit[kk] = true;}
			//}
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


	bool good_TOF = false;
	bool good_MWPC = false;

	if ((C2X_hit >=n_MWPC_hits) && (C2Y_hit >=n_MWPC_hits) && (C3X_hit >=n_MWPC_hits) && 
			(C3Y_hit >=n_MWPC_hits) && (C4X_hit >=n_MWPC_hits) && (C4Y_hit >=n_MWPC_hits)) good_MWPC = true;
	if (has_data_TOF > 0) good_TOF = true;

	if ((good_TOF) && (good_MWPC)) {

		good_event_counter++;

		for (Int_t j=0; j<256; j++){
			h_ADC_High_TARGET[j]->Fill(adc_high_target[j]);
			h_ADC_Low_TARGET[j]->Fill(adc_low_target[j]);
			h_TDC_LE_TARGET[j]->Fill(tdc_le_target[j][0]);
			h_TDC_TE_TARGET[j]->Fill(tdc_te_target[j][0]);

		}

		for (Int_t j=0; j<128; j++){
			h_ADC_High_SFT[j]->Fill(ADC_High_sft[j]);
			h_ADC_Low_SFT[j]->Fill(ADC_Low_sft[j]);
			h_TDC_LE_SFT[j]->Fill(TDC_LE_sft[j][0]);
			h_TDC_TE_SFT[j]->Fill(TDC_TE_sft[j][0]);
		}

		for (Int_t j=0; j<24; j++) {
			h_ADC_tof1[j]->Fill(ADC_tof1[j]);
		}

		for (Int_t j=0; j<12; j++) {
			h_TDC_tof1U[j]->Fill(TDC_tof1U[j]);
			h_TDC_tof1D[j]->Fill(TDC_tof1D[j]);
		}

		for (Int_t j=0; j<56; j++) {
			h_ADC_tof2[j]->Fill(ADC_tof2[j]);
		}

		for (Int_t j=0; j<12; j++) {
			h_TDC_tof2AO[j]->Fill(TDC_tof2AO[j]);
			h_TDC_tof2AI[j]->Fill(TDC_tof2AI[j]);
			h_TDC_tof2BO[j]->Fill(TDC_tof2BO[j]);
			h_TDC_tof2BI[j]->Fill(TDC_tof2BI[j]);
		}	


		for (Int_t j=0; j<512; j++) {
			h_MWPCADC[j]->Fill(MWPCADC[j]);
		}	
	}
}

if(QUIET==1) cout << "OK !" << endl;

//TFile *foutput;
//char Name_foutput[100];
//sprintf(Name_foutput,"%s/Hist_Run%dM.root",path_output, run_number);
//foutput = new TFile(Name_foutput,"RECREATE");

TDirectory *dir_Raw_TARGET = foutput->mkdir("RAW (TARGET)");
TDirectory *sub_Raw_ADC_High_TARGET = dir_Raw_TARGET->mkdir("Raw_ADC_High_TARGET");
TDirectory *sub_Raw_ADC_Low_TARGET = dir_Raw_TARGET->mkdir("Raw_ADC_Low_TARGET");
TDirectory *sub_Raw_TDC_LE_TARGET = dir_Raw_TARGET->mkdir("Raw_TDC_LE_TARGET");
TDirectory *sub_Raw_TDC_TE_TARGET = dir_Raw_TARGET->mkdir("Raw_TDC_TE_TARGET");

TDirectory *dir_Raw_SFT = foutput->mkdir("RAW (SFT)");
TDirectory *sub_Raw_ADC_High_SFT = dir_Raw_SFT->mkdir("Raw_ADC_High_SFT");
TDirectory *sub_Raw_ADC_Low_SFT = dir_Raw_SFT->mkdir("Raw_ADC_Low_SFT");
TDirectory *sub_Raw_TDC_LE_SFT = dir_Raw_SFT->mkdir("Raw_TDC_LE_SFT");
TDirectory *sub_Raw_TDC_TE_SFT = dir_Raw_SFT->mkdir("Raw_TDC_TE_SFT");

TDirectory *dir_Raw_TOF1 = foutput->mkdir("RAW (TOF1)");
TDirectory *sub_Raw_ADC_tof1 = dir_Raw_TOF1->mkdir("Raw_ADC_tof1");
TDirectory *sub_Raw_TDC_tof1U = dir_Raw_TOF1->mkdir("Raw_TDC_tof1U");
TDirectory *sub_Raw_TDC_tof1D = dir_Raw_TOF1->mkdir("Raw_TDC_tof1D");

TDirectory *dir_Raw_TOF2 = foutput->mkdir("RAW (TOF2)");
TDirectory *sub_Raw_ADC_tof2 = dir_Raw_TOF2->mkdir("Raw_ADC_tof2");
TDirectory *sub_Raw_TDC_tof2AO = dir_Raw_TOF2->mkdir("Raw_TDC_tof2AO");
TDirectory *sub_Raw_TDC_tof2AI = dir_Raw_TOF2->mkdir("Raw_TDC_tof2AI");
TDirectory *sub_Raw_TDC_tof2BO = dir_Raw_TOF2->mkdir("Raw_TDC_tof2BO");
TDirectory *sub_Raw_TDC_tof2BI = dir_Raw_TOF2->mkdir("Raw_TDC_tof2BI");

TDirectory *dir_Raw_MWPC = foutput->mkdir("RAW (MWPC)");
TDirectory *sub_Raw_MWPCADC = dir_Raw_MWPC->mkdir("Raw_ADC_MWPC");

//foutput->cd();

for (Int_t ii=0; ii<256; ii++){
sub_Raw_ADC_High_TARGET->cd();
h_ADC_High_TARGET[ii]->SetDirectory(sub_Raw_ADC_High_TARGET);
h_ADC_High_TARGET[ii]->Write();

sub_Raw_ADC_Low_TARGET->cd();
h_ADC_Low_TARGET[ii]->SetDirectory(sub_Raw_ADC_Low_TARGET);
h_ADC_Low_TARGET[ii]->Write();

sub_Raw_TDC_LE_TARGET->cd();
h_TDC_LE_TARGET[ii]->SetDirectory(sub_Raw_TDC_LE_TARGET);
h_TDC_LE_TARGET[ii]->Write();

sub_Raw_TDC_TE_TARGET->cd();
h_TDC_TE_TARGET[ii]->SetDirectory(sub_Raw_TDC_TE_TARGET);
h_TDC_TE_TARGET[ii]->Write();
}


for (Int_t ii=0; ii<128; ii++){
sub_Raw_ADC_High_SFT->cd();
h_ADC_High_SFT[ii]->SetDirectory(sub_Raw_ADC_High_SFT);
h_ADC_High_SFT[ii]->Write();

sub_Raw_ADC_Low_SFT->cd();
h_ADC_Low_SFT[ii]->SetDirectory(sub_Raw_ADC_Low_SFT);
h_ADC_Low_SFT[ii]->Write();

sub_Raw_TDC_LE_SFT->cd();
h_TDC_LE_SFT[ii]->SetDirectory(sub_Raw_TDC_LE_SFT);
h_TDC_LE_SFT[ii]->Write();

sub_Raw_TDC_TE_SFT->cd();
h_TDC_TE_SFT[ii]->SetDirectory(sub_Raw_TDC_TE_SFT);
h_TDC_TE_SFT[ii]->Write();
}

for (Int_t ii=0; ii<24; ii++){
sub_Raw_ADC_tof1->cd();
h_ADC_tof1[ii]->SetDirectory(sub_Raw_ADC_tof1);
h_ADC_tof1[ii]->Write();
}

for (Int_t ii=0; ii<12; ii++){
sub_Raw_TDC_tof1U->cd();
h_TDC_tof1U[ii]->SetDirectory(sub_Raw_TDC_tof1U);
h_TDC_tof1U[ii]->Write();

sub_Raw_TDC_tof1D->cd();
h_TDC_tof1D[ii]->SetDirectory(sub_Raw_TDC_tof1D);
h_TDC_tof1D[ii]->Write();
}

for (Int_t ii=0; ii<56; ii++){
sub_Raw_ADC_tof2->cd();
h_ADC_tof2[ii]->SetDirectory(sub_Raw_ADC_tof1);
h_ADC_tof2[ii]->Write();
}

for (Int_t ii=0; ii<12; ii++){
sub_Raw_TDC_tof2AO->cd();
h_TDC_tof2AO[ii]->SetDirectory(sub_Raw_TDC_tof2AO);
h_TDC_tof2AO[ii]->Write();

sub_Raw_TDC_tof2AI->cd();
h_TDC_tof2AI[ii]->SetDirectory(sub_Raw_TDC_tof2AI);
h_TDC_tof2AI[ii]->Write();

sub_Raw_TDC_tof2BO->cd();
h_TDC_tof2BO[ii]->SetDirectory(sub_Raw_TDC_tof2BO);
h_TDC_tof2BO[ii]->Write();

sub_Raw_TDC_tof2BI->cd();
h_TDC_tof2BI[ii]->SetDirectory(sub_Raw_TDC_tof2BI);
h_TDC_tof2BI[ii]->Write();
}

for (Int_t ii=0; ii<512; ii++){
sub_Raw_MWPCADC->cd();
h_MWPCADC[ii]->SetDirectory(sub_Raw_MWPCADC);
h_MWPCADC[ii]->Write();
}

cout << "" << endl;
cout << "Number of good events: " << good_event_counter << endl;
cout << "" << endl;
cout << "Histograms have been saved!" << endl;
cout << "" << endl;

foutput->Close();

return;

}

