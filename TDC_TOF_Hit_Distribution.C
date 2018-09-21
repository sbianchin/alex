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
  
void TDC_TOF_Hit_Distribution(Int_t run_number=5, Int_t flag=0, Int_t QUIET=0) { 

gStyle->Clear();
TH1::AddDirectory(kFALSE);
gStyle->SetOptStat(0);

char source_mapping[] = "SFT_Mapping_Oct14.txt";

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

Int_t ADC_TOF1[24];
Int_t ADC_TOF2[56];

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

char path_input[200];                   		char file_mapping[200];
sprintf(path_input,"%s",path_merged);          	sprintf(file_mapping,"../Mapping");
//sprintf(path_input,"/triumfcs/trshare/trek/E36/Data/April_2015/root");

char Name_finput[200];
sprintf(Name_finput,"%s/Run%dMS.root",path_input, run_number);

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

TH1D *h_TDC_tof1A[24];   char Title_TDC_tof1A[24][100];	char Name_TDC_tof1A[24][100];
TH1D *h_TDC_tof1B[24];   char Title_TDC_tof1B[24][100];	char Name_TDC_tof1B[24][100];
TH1D *h_TDC_tof1C[24];   char Title_TDC_tof1C[24][100];	char Name_TDC_tof1C[24][100];
TH1D *h_TDC_tof1D[24];   char Title_TDC_tof1D[24][100];	char Name_TDC_tof1D[24][100];

TH1D *h_TDC_tof2A[48];   char Title_TDC_tof2A[48][100];	char Name_TDC_tof2A[48][100];
TH1D *h_TDC_tof2B[48];   char Title_TDC_tof2B[48][100];	char Name_TDC_tof2B[48][100];
TH1D *h_TDC_tof2C[48];   char Title_TDC_tof2C[48][100];	char Name_TDC_tof2C[48][100];
TH1D *h_TDC_tof2D[48];   char Title_TDC_tof2D[48][100];	char Name_TDC_tof2D[48][100];

for (int i=0; i<24; i++) {
	sprintf(Title_TDC_tof1A[i],"Raw TDC (Ch. %d)  --  TOF1 - Counter A",i);
	sprintf(Name_TDC_tof1A[i],"TDC_TOF1 (Ch. %d) - Counter A",i);

	sprintf(Title_TDC_tof1B[i],"Raw TDC (Ch. %d)  --  TOF1 - Counter B",i);
	sprintf(Name_TDC_tof1B[i],"TDC_TOF1 (Ch. %d) - Counter B",i);

	sprintf(Title_TDC_tof1C[i],"Raw TDC (Ch. %d)  --  TOF1 - Counter C",i);
	sprintf(Name_TDC_tof1C[i],"TDC_TOF1 (Ch. %d) - Counter C",i);

	sprintf(Title_TDC_tof1D[i],"Raw TDC (Ch. %d)  --  TOF1 - Counter D",i);
	sprintf(Name_TDC_tof1D[i],"TDC_TOF1 (Ch. %d) - Counter D",i);
}

for (int i=0; i<48; i++) {
	sprintf(Title_TDC_tof2A[i],"Raw TDC (Ch. %d)  --  TOF2 - Counter A",i);
	sprintf(Name_TDC_tof2A[i],"TDC_TOF2 (Ch. %d) - Counter A",i);

	sprintf(Title_TDC_tof2B[i],"Raw TDC (Ch. %d)  --  TOF2 - Counter B ",i);
	sprintf(Name_TDC_tof2B[i],"TDC_TOF2 (Ch. %d) - Counter B",i);

	sprintf(Title_TDC_tof2C[i],"Raw TDC (Ch. %d)  --  TOF2 - Counter C",i);
	sprintf(Name_TDC_tof2C[i],"TDC_TOF2 (Ch. %d) - Counter C",i);

	sprintf(Title_TDC_tof2D[i],"Raw TDC (Ch. %d)  --  TOF2 - Counter D",i);
	sprintf(Name_TDC_tof2D[i],"TDC_TOF2 (Ch. %d) - Counter D",i);
}


for (int  i=0; i<24; i++) {
h_TDC_tof1A[i] = new TH1D(Name_TDC_tof1A[i],Title_TDC_tof1A[i],9,0,9);
h_TDC_tof1B[i] = new TH1D(Name_TDC_tof1B[i],Title_TDC_tof1B[i],9,0,9);
h_TDC_tof1C[i] = new TH1D(Name_TDC_tof1C[i],Title_TDC_tof1C[i],9,0,9);
h_TDC_tof1D[i] = new TH1D(Name_TDC_tof1D[i],Title_TDC_tof1D[i],9,0,9);
}

for (int  i=0; i<48; i++) {
h_TDC_tof2A[i] = new TH1D(Name_TDC_tof2A[i],Title_TDC_tof2A[i],9,0,9);
h_TDC_tof2B[i] = new TH1D(Name_TDC_tof2B[i],Title_TDC_tof2B[i],9,0,9);
h_TDC_tof2C[i] = new TH1D(Name_TDC_tof2C[i],Title_TDC_tof2C[i],9,0,9);
h_TDC_tof2D[i] = new TH1D(Name_TDC_tof2D[i],Title_TDC_tof2D[i],9,0,9);
}

for (int  i=0; i<24; i++) {
h_TDC_tof1A[i]->GetXaxis()->SetLabelOffset(999);
h_TDC_tof1A[i]->GetXaxis()->SetLabelSize(0);
h_TDC_tof1B[i]->GetXaxis()->SetLabelOffset(999);
h_TDC_tof1B[i]->GetXaxis()->SetLabelSize(0);
h_TDC_tof1C[i]->GetXaxis()->SetLabelOffset(999);
h_TDC_tof1C[i]->GetXaxis()->SetLabelSize(0);
h_TDC_tof1D[i]->GetXaxis()->SetLabelOffset(999);
h_TDC_tof1D[i]->GetXaxis()->SetLabelSize(0);
}

for (int  i=0; i<48; i++) {
h_TDC_tof2A[i]->GetXaxis()->SetLabelOffset(999);
h_TDC_tof2A[i]->GetXaxis()->SetLabelSize(0);
h_TDC_tof2B[i]->GetXaxis()->SetLabelOffset(999);
h_TDC_tof2B[i]->GetXaxis()->SetLabelSize(0);
h_TDC_tof2C[i]->GetXaxis()->SetLabelOffset(999);
h_TDC_tof2C[i]->GetXaxis()->SetLabelSize(0);
h_TDC_tof2D[i]->GetXaxis()->SetLabelOffset(999);
h_TDC_tof2D[i]->GetXaxis()->SetLabelSize(0);
}

//read all entries and fill the histograms
Int_t nentries = (Int_t)fChain->GetEntries();
cout <<  "Total Number of Entries :     " <<  nentries << endl;


cout << "  " << endl;
cout << "***************  TOF1  ***************" << endl;
cout << "  " << endl;

int TOF1CounterA[24] = {0};
int TOF1CounterB[24] = {0};
int TOF1CounterC[24] = {0};
int TOF1CounterD[24] = {0};

int TOF2CounterA[48] = {0};
int TOF2CounterB[48] = {0};
int TOF2CounterC[48] = {0};
int TOF2CounterD[48] = {0};

if(flag!=0) nentries=flag;

for (Int_t i=0; i<nentries; i++) {
//fChain_TOF1->GetEntry(i);
fChain->GetEntry(i);

	//for (Int_t j_TOF1=0; j_TOF1<24; j_TOF1++) {
	//   ADC_TOF1[j_TOF1] = ADC_tof1[j_TOF1]-par_temp_TOF1[1][j_TOF1];
	//}

	for (Int_t j_TOF1=0; j_TOF1<12; j_TOF1++) {
		ADC_TOF1U[j_TOF1] = ADC_tof1U[j_TOF1]-par_temp_TOF1[1][j_TOF1];
		ADC_TOF1D[j_TOF1] = ADC_tof1D[j_TOF1]-par_temp_TOF1[1][j_TOF1+12];
	}

	for(int i=0; i<12; i++){
		ADC_TOF1[i] = ADC_TOF1U[i];
		ADC_TOF1[i+12] = ADC_TOF1D[i];
	}


	if(QUIET==0){
		if(nentries<=30000){
			if(i%1000==1) cout<<"**** "<<i<<" events done"<<endl;
		}
		if(nentries>30000){
			if(i%10000==1) cout<<"**** "<<i<<" events done"<<endl;
		}
	}

	/// Filter TOF1 hits by distribution within peak

	for (int q=0; q<12; q++) {
		if (ADC_TOF1[q] > 0) {
			if (TDC_tof1U[q] == -1) TOF1CounterA[q]++;
			else if (TDC_tof1U[q] > TDC_TOF1_min && TDC_tof1U[q] < TDC_TOF1_max) TOF1CounterB[q]++;
			else if (TDC_tof1U[q] > 4000) TOF1CounterD[q]++;
			else {TOF1CounterC[q]++;}
		}
	}

	for (int q=12; q<24; q++) {
		if (ADC_TOF1[q] > 0) {
			if (TDC_tof1D[q-12] == -1) TOF1CounterA[q]++;
			else if (TDC_tof1D[q-12] > TDC_TOF1_min && TDC_tof1D[q-12] < TDC_TOF1_max) TOF1CounterB[q]++;
			else if (TDC_tof1D[q-12] > 4000) TOF1CounterD[q]++;
			else {TOF1CounterC[q]++;}
		}
	}

}


if(QUIET==1) cout << "OK !" << endl;

cout << "  " << endl;
cout << "  " << endl;
cout << "***************  TOF2  ***************" << endl;
cout << "  " << endl;

for (Int_t i=0; i<nentries; i++) {
//fChain_TOF2->GetEntry(i);
	fChain->GetEntry(i);

	//for (Int_t j_TOF2=0; j_TOF2<56; j_TOF2++) {
   	//	ADC_TOF2[j_TOF2] = ADC_tof2[j_TOF2]-par_temp_TOF2[1][j_TOF2];
	//}
	
	for (Int_t j_TOF2=0; j_TOF2<12; j_TOF2++) {
		ADC_TOF2AO[j_TOF2] = ADC_tof2AO[j_TOF2]-par_temp_TOF2[1][j_TOF2];
		ADC_TOF2AI[j_TOF2] = ADC_tof2AI[j_TOF2]-par_temp_TOF2[1][j_TOF2+24];
		ADC_TOF2BO[j_TOF2] = ADC_tof2BO[j_TOF2]-par_temp_TOF2[1][j_TOF2+12];
		ADC_TOF2BI[j_TOF2] = ADC_tof2BI[j_TOF2]-par_temp_TOF2[1][j_TOF2+36];
	}

	for(int i=0; i<12; i++){
		ADC_TOF2[i] = ADC_TOF2AO[i];
		ADC_TOF2[i+12] = ADC_TOF2BO[i];
		ADC_TOF2[i+24] = ADC_TOF2AI[i];
		ADC_TOF2[i+36] = ADC_TOF2BI[i];
	}

	if(QUIET==0){
		if(nentries<=30000){
			if(i%1000==1) cout<<"**** "<<i<<" events done"<<endl;
		}
		if(nentries>30000){
			if(i%10000==1) cout<<"**** "<<i<<" events done"<<endl;
		}
	}

	/// Filter TOF2 hits by distribution within peak

	for (int q=0; q<12; q++) {
		if (ADC_TOF2[q] > 0) {
			if (TDC_tof2AO[q] == -1) TOF2CounterA[q]++;
			else if (TDC_tof2AO[q] > TDC_TOF2_min && TDC_tof2AO[q] < TDC_TOF2_max) TOF2CounterB[q]++;

			else if (TDC_tof2AO[q] > 4000) TOF2CounterD[q]++;
			else {
				TOF2CounterC[q]++;
			}
		}
	}

	for (int q=12; q<24; q++) {
	//	if (q == 18) {
	//		if (ADC_TOF2[55] > 0) {
	//		if (TDC_tof2BO[q-12] == -1) TOF2CounterA[q]++;
	//		else if (TDC_tof2BO[q-12] > TDC_TOF2_min && TDC_tof2BO[q-12] < TDC_TOF2_max) TOF2CounterB[q]++;
	//		else if (TDC_tof2BO[q-12] > 4000) TOF2CounterD[q]++;
	//		else {TOF2CounterC[q]++;}
	//		}
	//	}
	//	else {
			if (ADC_TOF2[q] > 0) {
				if (TDC_tof2BO[q-12] == -1) TOF2CounterA[q]++;
				else if (TDC_tof2BO[q-12] > TDC_TOF2_min && TDC_tof2BO[q-12] < TDC_TOF2_max) TOF2CounterB[q]++;
				else if (TDC_tof2BO[q-12] > 4000) TOF2CounterD[q]++;
				else {TOF2CounterC[q]++;}
			}
	//	}
	}

	for (int q=24; q<36; q++) {
		if (ADC_TOF2[q] > 0) {
			if (TDC_tof2AI[q-24] == -1) TOF2CounterA[q]++;
			else if (TDC_tof2AI[q-24] > TDC_TOF2_min && TDC_tof2AI[q-24] < TDC_TOF2_max) TOF2CounterB[q]++;
			else if (TDC_tof2AI[q-24] > 4000) TOF2CounterD[q]++;
			else {TOF2CounterC[q]++;}
		}
	}

	for (int q=36; q<48; q++) {
		if (ADC_TOF2[q] > 0) {
			if (TDC_tof2BI[q-36] == -1) TOF2CounterA[q]++;
			else if (TDC_tof2BI[q-36] > TDC_TOF2_min && TDC_tof2BI[q-36] < TDC_TOF2_max) TOF2CounterB[q]++;
			else if (TDC_tof2BI[q-36] > 4000) TOF2CounterD[q]++;
			else {TOF2CounterC[q]++;}
		}
	}		
}

for (int  i=0; i<24; i++) {
	h_TDC_tof1A[i]->Fill(1,TOF1CounterA[i]);
	h_TDC_tof1B[i]->Fill(3,TOF1CounterB[i]);
	h_TDC_tof1C[i]->Fill(5,TOF1CounterC[i]);
	h_TDC_tof1D[i]->Fill(7,TOF1CounterD[i]);
}

for (int  i=0; i<48; i++) {
	h_TDC_tof2A[i]->Fill(1,TOF2CounterA[i]);
	h_TDC_tof2B[i]->Fill(3,TOF2CounterB[i]);
	h_TDC_tof2C[i]->Fill(5,TOF2CounterC[i]);
	h_TDC_tof2D[i]->Fill(7,TOF2CounterD[i]);
}

/// Set Histogram Maxima

int TOF1AMax[24] = {0};
int TOF1BMax[24] = {0};
int TOF1CMax[24] = {0};
int TOF1DMax[24] = {0};

int TOF2AMax[48] = {0};
int TOF2BMax[48] = {0};
int TOF2CMax[48] = {0};
int TOF2DMax[48] = {0};

for (Int_t j=0; j<24; j++) {
	TOF1AMax[j] = h_TDC_tof1A[j]->GetMaximum();
	TOF1BMax[j] = h_TDC_tof1B[j]->GetMaximum();
	TOF1CMax[j] = h_TDC_tof1C[j]->GetMaximum();
	TOF1DMax[j] = h_TDC_tof1D[j]->GetMaximum();
}

for (Int_t j=0; j<48; j++) {
	TOF2AMax[j] = h_TDC_tof2A[j]->GetMaximum();
	TOF2BMax[j] = h_TDC_tof2B[j]->GetMaximum();
	TOF2CMax[j] = h_TDC_tof2C[j]->GetMaximum();
	TOF2DMax[j] = h_TDC_tof2D[j]->GetMaximum();
}

float TOF1_height_max[24] = {0};
float TOF2_height_max[48] = {0};

for (Int_t j=0; j<24; j++) {
	if (TOF1AMax[j] > TOF1_height_max[j]) TOF1_height_max[j] = TOF1AMax[j];
	if (TOF1BMax[j] > TOF1_height_max[j]) TOF1_height_max[j] = TOF1BMax[j];
	if (TOF1CMax[j] > TOF1_height_max[j]) TOF1_height_max[j] = TOF1CMax[j];
	if (TOF1DMax[j] > TOF1_height_max[j]) TOF1_height_max[j] = TOF1DMax[j];
}

for (Int_t j=0; j<48; j++) {
	if (TOF2AMax[j] > TOF2_height_max[j]) TOF2_height_max[j] = TOF2AMax[j];
	if (TOF2BMax[j] > TOF2_height_max[j]) TOF2_height_max[j] = TOF2BMax[j];
	if (TOF2CMax[j] > TOF2_height_max[j]) TOF2_height_max[j] = TOF2CMax[j];
	if (TOF2DMax[j] > TOF2_height_max[j]) TOF2_height_max[j] = TOF2DMax[j];
}

for (Int_t j=0; j<24; j++) {
	TOF1_height_max[j] = TOF1_height_max[j]*1.1;
}

for (Int_t j=0; j<48; j++) {
	TOF2_height_max[j] = TOF2_height_max[j]*1.1;
}

////Draw TOF1 Histograms

char Name_Can_ADC_tof1B[100];			char Title_Can_ADC_tof1B[100];

sprintf(Name_Can_ADC_tof1B,"TOF1 TDC Distributions -- Run %d  (Ch. 0 - 23) -- Mapping Sorted",run_number);
sprintf(Title_Can_ADC_tof1B,"TOF1 TDC Distributions -- Run %d  (Ch. 0 - 23) -- Mapping Sorted",run_number);

char TOF1_String[100];
sprintf(TOF1_String, "TDC Hits %d < TDC < %d", TDC_TOF1_min, TDC_TOF1_max);

TLegend *l1;
l1 = new TLegend(0.75,0.75,0.99,0.99);
l1->AddEntry(h_TDC_tof1A[0], "Number of TDC Hits = -1", "F");
l1->AddEntry(h_TDC_tof1B[0], TOF1_String, "F");
l1->AddEntry(h_TDC_tof1C[0], "TDC Hits Outside Peak but >= 0", "F");
l1->AddEntry(h_TDC_tof1D[0], "TDC Hits > 4000", "F");

char TOF2_String[100];
sprintf(TOF2_String, "TDC Hits %d < TDC < %d", TDC_TOF2_min, TDC_TOF2_max);

TLegend *l2;
l2 = new TLegend(0.75,0.75,0.99,0.99);
l2->AddEntry(h_TDC_tof2A[0], "Number of TDC Hits = -1", "F");
l2->AddEntry(h_TDC_tof2B[0], TOF2_String, "F");
l2->AddEntry(h_TDC_tof2C[0], "TDC Hits Outside Peak but >= 0", "F");
l2->AddEntry(h_TDC_tof2D[0], "TDC Hits > 4000", "F");

TCanvas *c7;
	c7 = new TCanvas(Name_Can_ADC_tof1B,Title_Can_ADC_tof1B,900,333); 
	c7->Divide(6,4);
	c7->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=0; ican<12; ican++){
		c7->cd(TOF1_Mapping1[ican]);

		char TOF1Title3[250];
		sprintf(TOF1Title3, "TDC (Ch. %d) - TOF1 | Up-%d : Run %d", ican, ican+1, run_number);

		h_TDC_tof1B[ican]->SetTitle(TOF1Title3);
		h_TDC_tof1B[ican]->SetAxisRange(0, TOF1_height_max[ican],"Y");
		h_TDC_tof1B[ican]->SetLineColor(3);
		h_TDC_tof1B[ican]->SetFillStyle(3004);
		h_TDC_tof1B[ican]->SetFillColor(3);
		h_TDC_tof1B[ican]->Draw();

		h_TDC_tof1A[ican]->SetTitle(TOF1Title3);
		h_TDC_tof1A[ican]->SetAxisRange(0, TOF1_height_max[ican],"Y");
		h_TDC_tof1A[ican]->SetLineColor(2);
		h_TDC_tof1A[ican]->SetFillStyle(3004);
		h_TDC_tof1A[ican]->SetFillColor(2);
		h_TDC_tof1A[ican]->Draw("same");

		h_TDC_tof1C[ican]->SetTitle(TOF1Title3);
		h_TDC_tof1C[ican]->SetAxisRange(0, TOF1_height_max[ican],"Y");
		h_TDC_tof1C[ican]->SetLineColor(4);
		h_TDC_tof1C[ican]->SetFillStyle(3004);
		h_TDC_tof1C[ican]->SetFillColor(4);
		h_TDC_tof1C[ican]->Draw("same");

		h_TDC_tof1D[ican]->SetTitle(TOF1Title3);
		h_TDC_tof1D[ican]->SetAxisRange(0, TOF1_height_max[ican],"Y");
		h_TDC_tof1D[ican]->SetLineColor(1);
		h_TDC_tof1D[ican]->SetFillStyle(3004);
		h_TDC_tof1D[ican]->SetFillColor(1);
		h_TDC_tof1D[ican]->Draw("same");

		l1->Draw("same");
	}

	for(Int_t ican=0; ican<12; ican++){
		c7->cd(TOF1_Mapping2[ican]);

		char TOF1Title4[250];
		sprintf(TOF1Title4, "TDC (Ch. %d) - TOF1 | Down-%d : Run %d", ican+12, ican+1, run_number);
		h_TDC_tof1B[ican+12]->SetTitle(TOF1Title4);
		h_TDC_tof1B[ican+12]->SetAxisRange(0, TOF1_height_max[ican+12],"Y");
		h_TDC_tof1B[ican+12]->SetLineColor(3);
		h_TDC_tof1B[ican+12]->SetFillStyle(3004);
		h_TDC_tof1B[ican+12]->SetFillColor(3);
		h_TDC_tof1B[ican+12]->Draw();

		h_TDC_tof1A[ican+12]->SetTitle(TOF1Title4);
		h_TDC_tof1A[ican+12]->SetAxisRange(0, TOF1_height_max[ican+12],"Y");
		h_TDC_tof1A[ican+12]->SetLineColor(2);
		h_TDC_tof1A[ican+12]->SetFillStyle(3004);
		h_TDC_tof1A[ican+12]->SetFillColor(2);
		h_TDC_tof1A[ican+12]->Draw("same");

		h_TDC_tof1C[ican+12]->SetTitle(TOF1Title4);
		h_TDC_tof1C[ican+12]->SetAxisRange(0, TOF1_height_max[ican+12],"Y");
		h_TDC_tof1C[ican+12]->SetLineColor(4);
		h_TDC_tof1C[ican+12]->SetFillStyle(3004);
		h_TDC_tof1C[ican+12]->SetFillColor(4);
		h_TDC_tof1C[ican+12]->Draw("same");

		h_TDC_tof1D[ican+12]->SetTitle(TOF1Title4);
		h_TDC_tof1D[ican+12]->SetAxisRange(0, TOF1_height_max[ican+12],"Y");
		h_TDC_tof1D[ican+12]->SetLineColor(1);
		h_TDC_tof1D[ican+12]->SetFillStyle(3004);
		h_TDC_tof1D[ican+12]->SetFillColor(1);
		h_TDC_tof1D[ican+12]->Draw("same");

		l1->Draw("same");

	}

////Draw TOF2 Histograms

char Name_Can_ADC_tof2B[100];			char Title_Can_ADC_tof2B[100];

sprintf(Name_Can_ADC_tof2B,"TOF2 TDC Distributions -- Run %d  (Ch. 0 - 47) -- Mapping Sorted",run_number);
sprintf(Title_Can_ADC_tof2B,"TOF2 TDC Distributions -- Run %d  (Ch. 0 - 47) -- Mapping Sorted",run_number);

TCanvas *c8;
	c8 = new TCanvas(Name_Can_ADC_tof2B,Title_Can_ADC_tof2B,1200,500); 
	c8->Divide(8,6);
	c8->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=0; ican<12; ican++){
		c8->cd(TOF2_Mapping1[ican]);

		char TOF2Title5[250];
		sprintf(TOF2Title5, "Raw TDC (Ch. %d) - TOF2 | OutA-%d : Run %d", ican, ican+1, run_number);
		h_TDC_tof2B[ican]->SetTitle(TOF2Title5);
		h_TDC_tof2B[ican]->SetAxisRange(0, TOF2_height_max[ican],"Y");
		h_TDC_tof2B[ican]->SetLineColor(3);
		h_TDC_tof2B[ican]->SetFillStyle(3004);
		h_TDC_tof2B[ican]->SetFillColor(3);
		h_TDC_tof2B[ican]->Draw();

		h_TDC_tof2A[ican]->SetTitle(TOF2Title5);
		h_TDC_tof2A[ican]->SetAxisRange(0, TOF2_height_max[ican],"Y");
		h_TDC_tof2A[ican]->SetLineColor(2);
		h_TDC_tof2A[ican]->SetFillStyle(3004);
		h_TDC_tof2A[ican]->SetFillColor(2);
		h_TDC_tof2A[ican]->Draw("same");

		h_TDC_tof2C[ican]->SetTitle(TOF2Title5);
		h_TDC_tof2C[ican]->SetAxisRange(0, TOF2_height_max[ican],"Y");
		h_TDC_tof2C[ican]->SetLineColor(4);
		h_TDC_tof2C[ican]->SetFillStyle(3004);
		h_TDC_tof2C[ican]->SetFillColor(4);
		h_TDC_tof2C[ican]->Draw("same");

		h_TDC_tof2D[ican]->SetTitle(TOF2Title5);
		h_TDC_tof2D[ican]->SetAxisRange(0, TOF2_height_max[ican],"Y");
		h_TDC_tof2D[ican]->SetLineColor(1);
		h_TDC_tof2D[ican]->SetFillStyle(3004);
		h_TDC_tof2D[ican]->SetFillColor(1);
		h_TDC_tof2D[ican]->Draw("same");

		l2->Draw("same");
	}

	for(Int_t ican=0; ican<12; ican++){
		c8->cd(TOF2_Mapping2[ican]);
		
		char TOF2Title6[250];
		sprintf(TOF2Title6, "Raw TDC (Ch. %d) - TOF2 | InA-%d : Run %d", ican+24, ican+1, run_number);
		h_TDC_tof2B[ican+24]->SetTitle(TOF2Title6);
		h_TDC_tof2B[ican+24]->SetAxisRange(0, TOF2_height_max[ican+24],"Y");
		h_TDC_tof2B[ican+24]->SetLineColor(3);
		h_TDC_tof2B[ican+24]->SetFillStyle(3004);
		h_TDC_tof2B[ican+24]->SetFillColor(3);
		h_TDC_tof2B[ican+24]->Draw();

		h_TDC_tof2A[ican+24]->SetTitle(TOF2Title6);
		h_TDC_tof2A[ican+24]->SetAxisRange(0, TOF2_height_max[ican+24],"Y");
		h_TDC_tof2A[ican+24]->SetLineColor(2);
		h_TDC_tof2A[ican+24]->SetFillStyle(3004);
		h_TDC_tof2A[ican+24]->SetFillColor(2);
		h_TDC_tof2A[ican+24]->Draw("same");

		h_TDC_tof2C[ican+24]->SetTitle(TOF2Title6);
		h_TDC_tof2C[ican+24]->SetAxisRange(0, TOF2_height_max[ican+24],"Y");
		h_TDC_tof2C[ican+24]->SetLineColor(4);
		h_TDC_tof2C[ican+24]->SetFillStyle(3004);
		h_TDC_tof2C[ican+24]->SetFillColor(4);
		h_TDC_tof2C[ican+24]->Draw("same");

		h_TDC_tof2D[ican+24]->SetTitle(TOF2Title6);
		h_TDC_tof2D[ican+24]->SetAxisRange(0, TOF2_height_max[ican+24],"Y");
		h_TDC_tof2D[ican+24]->SetLineColor(1);
		h_TDC_tof2D[ican+24]->SetFillStyle(3004);
		h_TDC_tof2D[ican+24]->SetFillColor(1);
		h_TDC_tof2D[ican+24]->Draw("same");

		l2->Draw("same");
	}
	
	for(Int_t ican=0; ican<12; ican++){
		c8->cd(TOF2_Mapping3[ican]);

		char TOF2Title7[250];
		sprintf(TOF2Title7, "Raw TDC (Ch. %d) - TOF2 | OutB-%d : Run %d", ican+12, ican+1, run_number);
		h_TDC_tof2B[ican+12]->SetTitle(TOF2Title7);
		h_TDC_tof2B[ican+12]->SetAxisRange(0, TOF2_height_max[ican+12],"Y");
		h_TDC_tof2B[ican+12]->SetLineColor(3);
		h_TDC_tof2B[ican+12]->SetFillStyle(3004);
		h_TDC_tof2B[ican+12]->SetFillColor(3);
		h_TDC_tof2B[ican+12]->Draw();

		h_TDC_tof2A[ican+12]->SetTitle(TOF2Title7);
		h_TDC_tof2A[ican+12]->SetAxisRange(0, TOF2_height_max[ican+12],"Y");
		h_TDC_tof2A[ican+12]->SetLineColor(2);
		h_TDC_tof2A[ican+12]->SetFillStyle(3004);
		h_TDC_tof2A[ican+12]->SetFillColor(2);
		h_TDC_tof2A[ican+12]->Draw("same");

		h_TDC_tof2C[ican+12]->SetTitle(TOF2Title7);
		h_TDC_tof2C[ican+12]->SetAxisRange(0, TOF2_height_max[ican+12],"Y");
		h_TDC_tof2C[ican+12]->SetLineColor(4);
		h_TDC_tof2C[ican+12]->SetFillStyle(3004);
		h_TDC_tof2C[ican+12]->SetFillColor(4);
		h_TDC_tof2C[ican+12]->Draw("same");

		h_TDC_tof2D[ican+12]->SetTitle(TOF2Title7);
		h_TDC_tof2D[ican+12]->SetAxisRange(0, TOF2_height_max[ican+12],"Y");
		h_TDC_tof2D[ican+12]->SetLineColor(1);
		h_TDC_tof2D[ican+12]->SetFillStyle(3004);
		h_TDC_tof2D[ican+12]->SetFillColor(1);
		h_TDC_tof2D[ican+12]->Draw("same");

		l2->Draw("same");
	}

	for(Int_t ican=0; ican<12; ican++){
		c8->cd(TOF2_Mapping4[ican]);

		char TOF2Title8[250];
		sprintf(TOF2Title8, "Raw TDC (Ch. %d) - TOF2 | InB-%d : Run %d", ican+36, ican+1, run_number);
		h_TDC_tof2B[ican+36]->SetTitle(TOF2Title8);
		h_TDC_tof2B[ican+36]->SetAxisRange(0, TOF2_height_max[ican+36],"Y");
		h_TDC_tof2B[ican+36]->SetLineColor(3);
		h_TDC_tof2B[ican+36]->SetFillStyle(3004);
		h_TDC_tof2B[ican+36]->SetFillColor(3);
		h_TDC_tof2B[ican+36]->Draw();

		h_TDC_tof2A[ican+36]->SetTitle(TOF2Title8);
		h_TDC_tof2A[ican+36]->SetAxisRange(0, TOF2_height_max[ican+36],"Y");
		h_TDC_tof2A[ican+36]->SetLineColor(2);
		h_TDC_tof2A[ican+36]->SetFillStyle(3004);
		h_TDC_tof2A[ican+36]->SetFillColor(2);
		h_TDC_tof2A[ican+36]->Draw("same");

		h_TDC_tof2C[ican+36]->SetTitle(TOF2Title8);
		h_TDC_tof2C[ican+36]->SetAxisRange(0, TOF2_height_max[ican+36],"Y");
		h_TDC_tof2C[ican+36]->SetLineColor(4);
		h_TDC_tof2C[ican+36]->SetFillStyle(3004);
		h_TDC_tof2C[ican+36]->SetFillColor(4);
		h_TDC_tof2C[ican+36]->Draw("same");

		h_TDC_tof2D[ican+36]->SetTitle(TOF2Title8);
		h_TDC_tof2D[ican+36]->SetAxisRange(0, TOF2_height_max[ican+36],"Y");
		h_TDC_tof2D[ican+36]->SetLineColor(1);
		h_TDC_tof2D[ican+36]->SetFillStyle(3004);
		h_TDC_tof2D[ican+36]->SetFillColor(1);
		h_TDC_tof2D[ican+36]->Draw("same");

		l2->Draw("same");
	}

}

