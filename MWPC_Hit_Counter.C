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
  
void MWPC_Hit_Counter(Int_t Run_Number=5, Int_t ievt=0){ 

gStyle->SetOptStat(0);

char source_mapping[]="SFT_Mapping_Oct14.txt";  // Mapping file !!!

Int_t adc_high_target[256];   		Int_t ADC_High_TARGET[256];    	Double_t ADC_High_TARGET_corr[256];   	
Int_t adc_low_target[256]; 		Int_t ADC_Low_TARGET[256]; 	//Double_t ADC_Low_TARGET_corr[256];
Int_t tdc_le_target[256][16];      	Int_t TDC_LE_TARGET[256];    	Double_t TDC_LE_TARGET_corr[256];   	
Int_t tdc_te_target[256][16]; 		Int_t TDC_TE_TARGET[256]; 	//Double_t TDC_TE_TARGET_corr[256];

Int_t adc_high_sft[128];   		Int_t ADC_High_SFT[128];    	Double_t ADC_High_SFT_corr[128];   	
Int_t adc_low_sft[128]; 		Int_t ADC_Low_SFT[128]; 	//Double_t ADC_Low_SFT_corr[128];
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

Int_t MwpcADC[512];	Int_t MWPCADC[512];
Int_t fiber[128]={-1};

Int_t has_data_TARGET;
Int_t has_data_SFT;
Int_t has_data_TOF;
Int_t TOF1_counter;
Int_t TOF2_counter;
//Int_t has_data_MWPC;
Int_t good_event_counter = 0;
Int_t good_event_counter2 = 0; 
Int_t good_event_counter3 = 0; 
Int_t good_event_counter4 = 0;
Int_t Empty_TARGET_Counter = 0;
Int_t Empty_SFT_Counter = 0;   

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

char str1[100];
sprintf(str1,"Run_Number_%d_MWPC_Events_Cut_%d.txt",Run_Number,MWPC_thr);

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

fChain->SetBranchAddress("MWPCADC",MwpcADC);

//Int_t nentries;
Int_t nentries = (Int_t)fChain->GetEntries();//		Int_t nentries_SFT = (Int_t)fChain->GetEntries();
//if(nentries_TARGET==nentries_SFT) nentries = nentries_TARGET;
cout << "  " << endl;
cout << "****  Number of events: " << nentries << "  **** " <<endl;
cout << "  " << endl;

string filename = str1;
cout << filename << endl;

ofstream file (filename.c_str());
if (file.is_open()) {
	for(Int_t i=0; i<nentries; i++){
		fChain->GetEntry(i);	//	fChain_SFT->GetEntry(i);

		if(i%10000==1)	cout<<"**** "<<i<<" events done"<<endl;

		has_data_TARGET = 0;
		has_data_SFT = 0;
		has_data_TOF = 0;
	//	has_data_MWPC = 0;
		TOF1_counter = 0;
		TOF2_counter = 0;
	  	
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


		for (Int_t j_MWPC=0; j_MWPC<512; j_MWPC++) {
			MWPCADC[j_MWPC] = MwpcADC[j_MWPC]-MWPC_thr;
		}

		for (int q=0; q<24; q++) {
			if (ADC_TOF1[q]>0) TOF1_counter++;
		}

		for (int q=0; q<56; q++) {
			if (ADC_TOF2[q]>0) TOF2_counter++;
		}

		///

		for (Int_t p=0; p<256; p++) {
			if (ADC_High_TARGET[p] > 0) has_data_TARGET++;
		}

		for (Int_t p=0; p<128; p++) {
			fiber[p]=ADC_High_SFT[par_temp[1][p]];
			if (fiber[p] > 0) has_data_SFT++;
		}

		for (int kk=0; kk<12; kk++) {
			if (kk == 0) {
				if ((ADC_TOF2[0]>0 && ADC_TOF2[24]>0) || (ADC_TOF2[12]>0 && ADC_TOF2[36]>0)) {
					if((ADC_TOF1[0] > 0) || (ADC_TOF1[12] > 0) || (ADC_TOF1[1] > 0) || (ADC_TOF1[13] > 0) || (ADC_TOF1[11] > 0) || (ADC_TOF1[23] > 0)) has_data_TOF++;
				} 
			}

			else if (kk == 6) {
				if ((ADC_TOF2[6]>0 && ADC_TOF2[30]>0) || (ADC_TOF2[55]>0 && ADC_TOF2[42]>0)) {
					if((ADC_TOF1[6] > 0) || (ADC_TOF1[18] > 0) || (ADC_TOF1[7] > 0) || (ADC_TOF1[19] > 0) || (ADC_TOF1[5] > 0) || (ADC_TOF1[17] > 0)) has_data_TOF++;
				}
			}

			else if (kk == 11) {
				if ((ADC_TOF2[11]>0 && ADC_TOF2[35]>0) || (ADC_TOF2[23]>0 && ADC_TOF2[47]>0)) {
					if((ADC_TOF1[11] > 0) || (ADC_TOF1[23] > 0) || (ADC_TOF1[0] > 0) || (ADC_TOF1[12] > 0) || (ADC_TOF1[10] > 0) || (ADC_TOF1[22] > 0)) has_data_TOF++;
				}
			}
		
			else {
				if ((ADC_TOF2[kk]>0 && ADC_TOF2[kk+24]>0) || (ADC_TOF2[kk+12]>0 && ADC_TOF2[kk+36]>0)) {
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
			if (q >= 256 && q <= 271) {
				if (MWPCADC[q]>0) {
					C4Y_hit++;
					if (q==256) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C4Y_hitA++;}
					else if (q==271) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C4Y_hitA++;}
					else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C4Y_hitA++;}
				}
			}
			if (q >= 272 && q <= 287) {
				if (MWPCADC[q]>0) {
					C4Y_hit++;
					if (q==272) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C4Y_hitA++;}
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
			if (q >= 320 && q <= 335) {
				if (MWPCADC[q]>0) {
					C4X_hit++;
					if (q==320) {if (MWPCADC[q+1]>0) C4X_hitA++;}
					else if (q==335) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C4X_hitA++;}
					else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C4X_hitA++;}
				}
			}
			if (q >= 336 && q <= 351) {
				if (MWPCADC[q]>0) {
					C4X_hit++;
					if (q==336) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C4X_hitA++;}
					else if (q==351) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C4X_hitA++;}
					else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C4X_hitA++;}
				}
			}
			if (q >= 352 && q <= 367) {
				if (MWPCADC[q]>0) {
					C4X_hit++;
					if (q==352) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C4X_hitA++;}
					else if (q==367) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C4X_hitA++;}
					else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C4X_hitA++;}
				}
			}
			if (q >= 368 && q <= 383) {
				if (MWPCADC[q]>0) {
					C4X_hit++;
					if (q==368) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C4X_hitA++;}
					else if (q==383) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C4X_hitA++;}
					else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C4X_hitA++;}
				}
			}
			if (q >= 384 && q <= 399) {
				if (MWPCADC[q]>0) {
					C4X_hit++;
					if (q==384) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C4X_hitA++;}
					else if (q==399) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C4X_hitA++;}
					else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C4X_hitA++;}
				}
			}
			if (q >= 400 && q <= 415) {
				if (MWPCADC[q]>0) {
					C4X_hit++;
					if (q==400) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C4X_hitA++;}
					else if (q==415) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C4X_hitA++;}
					else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C4X_hitA++;}
				}
			}
			if (q >= 416 && q <= 431) {
				if (MWPCADC[q]>0) {
					C4X_hit++;
					if (q==416) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C4X_hitA++;}
					else if (q==431) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C4X_hitA++;}
					else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C4X_hitA++;}
				}
			}
			if (q >= 432 && q <= 447) {
				if (MWPCADC[q]>0) {
					C4X_hit++;
					if (q==432) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C4X_hitA++;}
					else if (q==447) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C4X_hitA++;}
					else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C4X_hitA++;}
				}
			}

			//C3 Counters
			if (q >= 448 && q <= 463) {
				if (MWPCADC[q]>0) {
					C3Y_hit++;
					if (q==448) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C3Y_hitA++;}
					else if (q==463) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C3Y_hitA++;}
					else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C3Y_hitA++;}
				}
			}
			if (q >= 464 && q <= 479) {
				if (MWPCADC[q]>0) {
					C3Y_hit++;
					if (q==464) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C3Y_hitA++;}
					else if (q==479) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C3Y_hitA++;}
					else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C3Y_hitA++;}
				}
			}
			if (q >= 480 && q <= 495) {
				if (MWPCADC[q]>0) {
					C3X_hit++;
					if (q==480) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C3X_hitA++;}
					else if (q==495) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C3X_hitA++;}
					else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C3X_hitA++;}
				}
			}
			if (q >= 496 && q <= 511) {
				if (MWPCADC[q]>0) {
					C3X_hit++;
					if (q==480) {if (MWPCADC[q+1]>0 || MWPCADC[q+2]>0) C3X_hitA++;}
					else if (q==495) {if (MWPCADC[q-1]>0 || MWPCADC[q-2]>0) C3X_hitA++;}
					else {if (MWPCADC[q-1]>0 || MWPCADC[q+1]>0) C3X_hitA++;}
				}
			}
		}
		///

		if (TOF1_counter > 0 && TOF2_counter > 0) {
			good_event_counter++;
			if (has_data_TOF > 0) { 
				good_event_counter2++;
				if (((C2X_hitA >= 2) && (C2Y_hitA >= 2)) && ((C3X_hitA >= 2) && (C3Y_hitA >= 2)) && ((C4X_hitA >= 2) && (C4Y_hitA >= 12))) {
					good_event_counter3++;
					if (has_data_TARGET > 10) good_event_counter4++;
				}
			}
		}

		if (((C2X_hitA >= 2) && (C2Y_hitA >= 2)) && ((C3X_hitA >= 2) && (C3Y_hitA >= 2)) && ((C4X_hitA >= 2) && (C4Y_hitA >= 12)) 
			&& (has_data_TOF > 0) && (has_data_TARGET > 10)) { 
			file << i;
			file << "\n";
		}
	}

file.close();

}

double eff = (double(double(good_event_counter)/double(nentries)))*100;
double eff2 = (double(double(good_event_counter2)/double(nentries)))*100;
double eff3 = (double(double(good_event_counter3)/double(nentries)))*100;
double eff4 = (double(double(good_event_counter4)/double(nentries)))*100;

cout << "" << endl;
cout << "TOF1 & TOF2 Hit: " << good_event_counter << " -- Percentage of total entries: " << eff << "%" <<endl;
cout << "Good Gap: " << good_event_counter2 << " -- Percentage of total entries: " << eff2 << "%" <<endl;
cout << "Good MWPC: " << good_event_counter3 << " -- Percentage of total entries: " << eff3 << "%" <<endl;
cout << "Good TARGET: " << good_event_counter4 << " -- Percentage of total entries: " << eff4 << "%" <<endl;
cout << "" << endl;

TH1D *h_count[5];	char Title_h_count[5][100];	char Name_h_count[5][100];
for (int i=0; i<5; i++) {
	sprintf(Title_h_count[i],"MWPC Hit Counter -- Run %d", Run_Number); 
	sprintf(Name_h_count[i],"Histogram Counter %d", i+1); 
	h_count[i] = new TH1D(Name_h_count[i],Title_h_count[i],11,0,11);
	h_count[i]->GetXaxis()->SetLabelOffset(999);
	h_count[i]->GetXaxis()->SetLabelSize(0);
}

h_count[0]->Fill(1,nentries);
h_count[1]->Fill(3,good_event_counter);
h_count[2]->Fill(5,good_event_counter2);
h_count[3]->Fill(7,good_event_counter3);
h_count[4]->Fill(9,good_event_counter4);

char Name_Can2[100];			char Title_Can2[100];
sprintf(Name_Can2,"MWPC Hit Counter -- Run %d", Run_Number);
sprintf(Title_Can2,"MWPC Hit Counter -- Run %d", Run_Number);

TCanvas *c1;
c1 = new TCanvas(Name_Can2,Title_Can2,1000,500); 

TLegend *l1;
l1 = new TLegend(0.75,0.75,0.99,0.99);
l1->AddEntry(h_count[0], "Total Entries", "F");
l1->AddEntry(h_count[1], "TOF1 & TOF2 Hit", "F");
l1->AddEntry(h_count[2], "Good Gap", "F");
l1->AddEntry(h_count[3], "Good MWPC", "F");
l1->AddEntry(h_count[4], "Good TARGET", "F");

h_count[0]->SetLineColor(1);
h_count[0]->SetFillStyle(3004);
h_count[0]->SetFillColor(1);
h_count[0]->Draw();

h_count[1]->SetLineColor(2);
h_count[1]->SetFillStyle(3004);
h_count[1]->SetFillColor(2);
h_count[1]->Draw("same");

h_count[2]->SetLineColor(3);
h_count[2]->SetFillStyle(3004);
h_count[2]->SetFillColor(3);
h_count[2]->Draw("same");

h_count[3]->SetLineColor(4);
h_count[3]->SetFillStyle(3004);
h_count[3]->SetFillColor(4);
h_count[3]->Draw("same");

h_count[4]->SetLineColor(6);
h_count[4]->SetFillStyle(3004);
h_count[4]->SetFillColor(6);
h_count[4]->Draw("same");

l1->Draw("same");

}




