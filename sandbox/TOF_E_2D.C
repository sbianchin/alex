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
  
void TOF_E_2D(Int_t run_number=5, Int_t flag=0, Int_t QUIET=0) { 


	gStyle->Clear();
	TH1::AddDirectory(kFALSE);
	gStyle->SetOptStat(1111111);

	Int_t TDC_min_TOF1 = TDC_TOF1_min;
	Int_t TDC_max_TOF1 = TDC_TOF1_max;
	Int_t TDC_min_TOF2 = TDC_TOF2_min;
	Int_t TDC_max_TOF2 = TDC_TOF2_max;

	///Adjust chart ranges

	Int_t TOF1_xmin = 0;
	Int_t TOF1_xmax = 4000;
	Int_t TOF1_ymin = 0;
	Int_t TOF1_ymax = 4000;

	Int_t TOF2_xmin = 0;
	Int_t TOF2_xmax = 4000;
	Int_t TOF2_ymin = 0;
	Int_t TOF2_ymax = 4000;

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

	//// Canvas board coordinates for TOF histograms

	Int_t TOF1_Mapping1[12] = {1,3,5,7,9,11,13,15,17,19,21,23};
	Int_t TOF1_Mapping2[12] = {2,4,6,8,10,12,14,16,18,20,22,24};

	Int_t TOF2_Mapping1[12] = {1,5,9,13,17,21,25,29,33,37,41,45};
	Int_t TOF2_Mapping2[12] = {2,6,10,14,18,22,26,30,34,38,42,46};
	Int_t TOF2_Mapping3[12] = {3,7,11,15,19,23,27,31,35,39,43,47};
	Int_t TOF2_Mapping4[12] = {4,8,12,16,20,24,28,32,36,40,44,48};	

	char path_input[200];                    	 char file_mapping[200];
	sprintf(path_input, "%s", path_merged);      sprintf(file_mapping,"../Mapping");
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

	TH2D *h_ADC_TDC_tof1[24];   char Title_ADC_TDC_tof1[24][100];	char Name_ADC_TDC_tof1[24][100];
	TH2D *h_ADC_TDC_tof1_peak[24];   char Title_ADC_TDC_tof1_peak[24][100];	char Name_ADC_TDC_tof1_peak[24][100];

	TH2D *h_ADC_TDC_tof2[48];   char Title_ADC_TDC_tof2[48][100];	char Name_ADC_TDC_tof2[48][100];
	TH2D *h_ADC_TDC_tof2_peak[48];   char Title_ADC_TDC_tof2_peak[48][100];	char Name_ADC_TDC_tof2_peak[48][100];

	char Title_ADC_TDC_TOF1_UP[12][100];		char Name_ADC_TDC_TOF1_UP[12][100];
	char Title_ADC_TDC_TOF1_DOWN[12][100];		char Name_ADC_TDC_TOF1_DOWN[12][100];
	char Title_ADC_TDC_TOF2_AO[12][100];		char Name_ADC_TDC_TOF2_AO[12][100];
	char Title_ADC_TDC_TOF2_AI[12][100];		char Name_ADC_TDC_TOF2_AI[12][100];	
	char Title_ADC_TDC_TOF2_BO[12][100];		char Name_ADC_TDC_TOF2_BO[12][100];
	char Title_ADC_TDC_TOF2_BI[12][100];		char Name_ADC_TDC_TOF2_BI[12][100];

	for(int i=0; i<12; i++){
		sprintf(Title_ADC_TDC_TOF1_UP[i], "TOF1  |  ADC vs TDC  --  Gap %d - UP   (Run %d)", i+1, run_number); 
		sprintf(Name_ADC_TDC_TOF1_UP[i], "TOF1 (ch. %d)", i); 

		sprintf(Title_ADC_TDC_TOF1_DOWN[i], "TOF1  |  ADC vs TDC  --  Gap %d - DOWN   (Run %d)", i+1, run_number); 
		sprintf(Name_ADC_TDC_TOF1_DOWN[i], "TOF1 (ch. %d)", i+12); 

		sprintf(Title_ADC_TDC_TOF2_AO[i], "TOF2  |  ADC vs TDC  --  Gap %d - A OUT   (Run %d)", i+1, run_number); 
		sprintf(Name_ADC_TDC_TOF2_AO[i], "TOF2 (ch. %d)", i); 

		sprintf(Title_ADC_TDC_TOF2_AI[i], "TOF2  |  ADC vs TDC  --  Gap %d - A IN   (Run %d)", i+1, run_number); 
		sprintf(Name_ADC_TDC_TOF2_AI[i], "TOF2 (ch. %d)", i+24); 

		sprintf(Title_ADC_TDC_TOF2_BO[i], "TOF2  |  ADC vs TDC  --  Gap %d - B OUT   (Run %d)", i+1, run_number); 
		sprintf(Name_ADC_TDC_TOF2_BO[i], "TOF2 (ch. %d)", i+12); 

		sprintf(Title_ADC_TDC_TOF2_BI[i], "TOF2  |  ADC vs TDC  --  Gap %d - B IN   (Run %d)", i+1, run_number); 
		sprintf(Name_ADC_TDC_TOF2_BI[i], "TOF2 (ch. %d)", i+36); 
	}


	for (int i=0; i<24; i++) {
		sprintf(Title_ADC_TDC_tof1[i],"ADC vs. TDC (Ch. %d)  --  TOF1",i);
		sprintf(Name_ADC_TDC_tof1[i],"ADC vs. TDC (Ch. %d) - TOF1",i);
		sprintf(Title_ADC_TDC_tof1_peak[i],"ADC vs. TDC (Ch. %d)  --  TOF1 Peak",i);
		sprintf(Name_ADC_TDC_tof1_peak[i],"ADC vs. TDC (Ch. %d) - TOF1 Peak",i);
	}

	for (int i=0; i<48; i++) {
		sprintf(Title_ADC_TDC_tof2[i],"ADC vs. TDC (Ch. %d)  --  TOF2",i);
		sprintf(Name_ADC_TDC_tof2[i],"ADC vs. TDC (Ch. %d) - TOF2",i);
		sprintf(Title_ADC_TDC_tof2_peak[i],"ADC vs. TDC (Ch. %d)  --  TOF2 Peak",i);
		sprintf(Name_ADC_TDC_tof2_peak[i],"ADC vs. TDC (Ch. %d) - TOF2 Peak",i);
	}


	for(int i=0; i<12; i++){
		h_ADC_TDC_tof1[i] = new TH2D(Name_ADC_TDC_TOF1_UP[i], Title_ADC_TDC_TOF1_UP[i], 256, 0, 4000, 128, 0, 4000);
		h_ADC_TDC_tof1[i+12] = new TH2D(Name_ADC_TDC_TOF1_DOWN[i], Title_ADC_TDC_TOF1_DOWN[i], 256, 0, 4000, 128, 0, 4000);

		h_ADC_TDC_tof2[i] = new TH2D(Name_ADC_TDC_TOF2_AO[i], Title_ADC_TDC_TOF2_AO[i], 256, 0, 4000, 128, 0, 4000);
		h_ADC_TDC_tof2[i+12] = new TH2D(Name_ADC_TDC_TOF2_BO[i], Title_ADC_TDC_TOF2_BO[i], 256, 0, 4000, 128, 0, 4000);
		h_ADC_TDC_tof2[i+24] = new TH2D(Name_ADC_TDC_TOF2_AI[i], Title_ADC_TDC_TOF2_AI[i], 256, 0, 4000, 128, 0, 4000);
		h_ADC_TDC_tof2[i+36] = new TH2D(Name_ADC_TDC_TOF2_BI[i], Title_ADC_TDC_TOF2_BI[i], 256, 0, 4000, 128, 0, 4000);
	}

	for (int  i=0; i<24; i++) {
	//h_ADC_TDC_tof1[i] = new TH2D(Name_ADC_TDC_tof1[i],Title_ADC_TDC_tof1[i],1000,0,4000,1000,0,4000);
	h_ADC_TDC_tof1_peak[i] = new TH2D(Name_ADC_TDC_tof1_peak[i],Title_ADC_TDC_tof1_peak[i],256,0,4000,128,0,4000);
	}

	for (int  i=0; i<48; i++) {
	//h_ADC_TDC_tof2[i] = new TH2D(Name_ADC_TDC_tof2[i],Title_ADC_TDC_tof2[i],1000,0,4000,1000,0,4000);
	h_ADC_TDC_tof2_peak[i] = new TH2D(Name_ADC_TDC_tof2_peak[i],Title_ADC_TDC_tof2_peak[i],256,0,4000,128,0,4000);
	}

	//read all entries and fill the histograms
	Int_t nentries = (Int_t)fChain->GetEntries();
	cout <<  "Total Number of Entries :     " <<  nentries << endl;

	cout << "  " << endl;

	//cout << "TDC TOF1 min: " << TDC_TOF1_min << endl;
	//cout << "TDC TOF1 max: " << TDC_TOF1_max << endl;
	//cout << "TDC TOF2 min: " << TDC_TOF2_min << endl;
	//cout << "TDC TOF2 max: " << TDC_TOF2_max << endl;

	cout << "  " << endl;
	cout << "***************  TOF1  ***************" << endl;
	cout << "  " << endl;

	if(flag!=0) nentries=flag;

	for (Int_t i=0; i<nentries; i++) {
		//fChain_TOF1->GetEntry(i);
		fChain->GetEntry(i);
		
		for (Int_t j_TOF1=0; j_TOF1<12; j_TOF1++) {
			ADC_TOF1U[j_TOF1] = ADC_tof1U[j_TOF1]-par_temp_TOF1[1][j_TOF1];
			ADC_TOF1D[j_TOF1] = ADC_tof1D[j_TOF1]-par_temp_TOF1[1][j_TOF1+12];
		}

		for(int i=0; i<12; i++){
			ADC_TOF1[i] = ADC_TOF1U[i];
			ADC_TOF1[i+12] = ADC_TOF1D[i];
		}


		//for (Int_t j_TOF1=0; j_TOF1<24; j_TOF1++) {
	   	//	ADC_TOF1[j_TOF1] = ADC_tof1[j_TOF1]-par_temp_TOF1[1][j_TOF1];
		//}

		if(QUIET==0){
			if(nentries<=30000){
				if(i%1000==1) cout<<"**** "<<i<<" events done"<<endl;
			}
			if(nentries>30000){
				if(i%10000==1) cout<<"**** "<<i<<" events done"<<endl;
			}
		}

		/// Plot ADC vs. TDC for TOF1

		for (int q=0; q<12; q++) {
			h_ADC_TDC_tof1[q]->Fill(TDC_tof1U[q],ADC_TOF1[q]);
			if (TDC_tof1U[q] > TDC_min_TOF1 && TDC_tof1U[q] < TDC_max_TOF1) h_ADC_TDC_tof1_peak[q]->Fill(TDC_tof1U[q],ADC_TOF1[q]);
		}

		for (int q=12; q<24; q++) {
			h_ADC_TDC_tof1[q]->Fill(TDC_tof1D[q-12],ADC_TOF1[q]);
			if (TDC_tof1D[q-12] > TDC_min_TOF1 && TDC_tof1D[q-12] < TDC_max_TOF1) h_ADC_TDC_tof1_peak[q]->Fill(TDC_tof1D[q-12],ADC_TOF1[q]);
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


		//for (Int_t j_TOF2=0; j_TOF2<56; j_TOF2++) {
   		//	ADC_TOF2[j_TOF2] = ADC_tof2[j_TOF2]-par_temp_TOF2[1][j_TOF2];
		//}

		if(QUIET==0){
			if(nentries<=30000){
				if(i%1000==1) cout<<"**** "<<i<<" events done"<<endl;
			}
			if(nentries>30000){
				if(i%10000==1) cout<<"**** "<<i<<" events done"<<endl;
			}
		}

		/// Plot ADC vs. TDC for TOF2

		for (int q=0; q<12; q++) {
			h_ADC_TDC_tof2[q]->Fill(TDC_tof2AO[q],ADC_TOF2[q]);
			if (TDC_tof2AO[q] > TDC_min_TOF2 && TDC_tof2AO[q] < TDC_max_TOF2) h_ADC_TDC_tof2_peak[q]->Fill(TDC_tof2AO[q],ADC_TOF2[q]);
		}

		for (int q=12; q<24; q++) {
			//if (q == 18) {
			//	h_ADC_TDC_tof2[q]->Fill(TDC_tof2BO[q-12],ADC_TOF2[55]);
			//	if (TDC_tof2BO[q-12] > TDC_min_TOF2 && TDC_tof2BO[q-12] < TDC_max_TOF2) h_ADC_TDC_tof2_peak[q]->Fill(TDC_tof2BO[q-12],ADC_TOF2[55]);
			//}
			//else {
				h_ADC_TDC_tof2[q]->Fill(TDC_tof2BO[q-12],ADC_TOF2[q]);
				if (TDC_tof2BO[q-12] > TDC_min_TOF2 && TDC_tof2BO[q-12] < TDC_max_TOF2) h_ADC_TDC_tof2_peak[q]->Fill(TDC_tof2BO[q-12],ADC_TOF2[q]);
			//}
		}

		for (int q=24; q<36; q++) {
			h_ADC_TDC_tof2[q]->Fill(TDC_tof2AI[q-24],ADC_TOF2[q]);
			if (TDC_tof2AI[q-24] > TDC_min_TOF2 && TDC_tof2AI[q-24] < TDC_max_TOF2) h_ADC_TDC_tof2_peak[q]->Fill(TDC_tof2AI[q-24],ADC_TOF2[q]);
		}

		for (int q=36; q<48; q++) {
			h_ADC_TDC_tof2[q]->Fill(TDC_tof2BI[q-36],ADC_TOF2[q]);
			if (TDC_tof2BI[q-36] > TDC_min_TOF2 && TDC_tof2BI[q-36] < TDC_max_TOF2) h_ADC_TDC_tof2_peak[q]->Fill(TDC_tof2BI[q-36],ADC_TOF2[q]);
		}		
	}

	//// Draw TOF1 Histograms

	char Name_Can_ADC_tof1B[100];			char Title_Can_ADC_tof1B[100];

	sprintf(Name_Can_ADC_tof1B,"TOF1 TDC Distributions -- Run %d  (Ch. 0 - 23) -- Mapping Sorted",run_number);
	sprintf(Title_Can_ADC_tof1B,"TOF1 TDC Distributions -- Run %d  (Ch. 0 - 23) -- Mapping Sorted",run_number);

	TCanvas *c7;
	c7 = new TCanvas(Name_Can_ADC_tof1B,Title_Can_ADC_tof1B,900,333); 
	c7->Divide(6,4);
	c7->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=0; ican<12; ican++){
		c7->cd(TOF1_Mapping1[ican]);
		h_ADC_TDC_tof1[ican]->SetAxisRange(TOF1_xmin, TOF1_xmax, "X");
		h_ADC_TDC_tof1[ican]->SetAxisRange(TOF1_ymin, TOF1_ymax, "Y");
		//h_ADC_TDC_tof1_peak[ican]->SetMarkerColor(3);
		//h_ADC_TDC_tof1_peak[ican]->SetMarkerSize(50);
		//h_ADC_TDC_tof1_peak[ican]->Draw();
		h_ADC_TDC_tof1[ican]->Draw("colz");
	}

	for(Int_t ican=0; ican<12; ican++){
		c7->cd(TOF1_Mapping2[ican]);
		h_ADC_TDC_tof1[ican+12]->SetAxisRange(TOF1_xmin, TOF1_xmax, "X");
		h_ADC_TDC_tof1[ican+12]->SetAxisRange(TOF1_ymin, TOF1_ymax, "Y");
		//h_ADC_TDC_tof1_peak[ican+12]->SetMarkerColor(3);
		//h_ADC_TDC_tof1_peak[ican+12]->SetMarkerSize(50);
		//h_ADC_TDC_tof1_peak[ican+12]->Draw();
		h_ADC_TDC_tof1[ican+12]->Draw("colz");
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
		h_ADC_TDC_tof2[ican]->SetAxisRange(TOF2_xmin, TOF2_xmax, "X");
		h_ADC_TDC_tof2[ican]->SetAxisRange(TOF2_ymin, TOF2_ymax, "Y");
		//h_ADC_TDC_tof2_peak[ican]->SetMarkerColor(3);
		//h_ADC_TDC_tof2_peak[ican]->SetMarkerSize(50);
		//h_ADC_TDC_tof2_peak[ican]->Draw();
		h_ADC_TDC_tof2[ican]->Draw("colz");
	}

	for(Int_t ican=0; ican<12; ican++){
		c8->cd(TOF2_Mapping2[ican]);
		h_ADC_TDC_tof2[ican+24]->SetAxisRange(TOF2_xmin, TOF2_xmax, "X");
		h_ADC_TDC_tof2[ican+24]->SetAxisRange(TOF2_ymin, TOF2_ymax, "Y");
		//h_ADC_TDC_tof2_peak[ican+24]->SetMarkerColor(3);
		//h_ADC_TDC_tof2_peak[ican+24]->SetMarkerSize(50);
		//h_ADC_TDC_tof2_peak[ican+24]->Draw();
		h_ADC_TDC_tof2[ican+24]->Draw("colz");
	}
	
	for(Int_t ican=0; ican<12; ican++){
		c8->cd(TOF2_Mapping3[ican]);
		h_ADC_TDC_tof2[ican+12]->SetAxisRange(TOF2_xmin, TOF2_xmax, "X");
		h_ADC_TDC_tof2[ican+12]->SetAxisRange(TOF2_ymin, TOF2_ymax, "Y");
		//h_ADC_TDC_tof2_peak[ican+12]->SetMarkerColor(3);
		//h_ADC_TDC_tof2_peak[ican+12]->SetMarkerSize(50);
		//h_ADC_TDC_tof2_peak[ican+12]->Draw();
		h_ADC_TDC_tof2[ican+12]->Draw("colz");
	}

	for(Int_t ican=0; ican<12; ican++){
		c8->cd(TOF2_Mapping4[ican]);
		h_ADC_TDC_tof2[ican+36]->SetAxisRange(TOF2_xmin, TOF2_xmax, "X");
		h_ADC_TDC_tof2[ican+36]->SetAxisRange(TOF2_ymin, TOF2_ymax, "Y");
		//h_ADC_TDC_tof2_peak[ican+36]->SetMarkerColor(3);
		//h_ADC_TDC_tof2_peak[ican+36]->SetMarkerSize(50);
		//h_ADC_TDC_tof2_peak[ican+36]->Draw();
		h_ADC_TDC_tof2[ican+36]->Draw("colz");
	}
}

