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
#include "TMarker.h"
#include "ANAPATH.h"
#include "mappings.h"
#include "Cuts_and_Windows.h"

#endif


void Tree2Histo(Int_t run_number=5, Int_t flag=0, Int_t QUIET=0) { 

	gStyle->SetOptStat(111111111);

	//Int_t ADC_tof1[24];		Int_t TDC_tof1[24];
	//Int_t ADC_tof2[56];		Int_t TDC_tof2[56];


	//Int_t ADC_tof[96];
	//Int_t TDC_tof[96];

	//Int_t ADC_ac[24];
	//Int_t ADC_ttc[12];
	//Int_t TDC_ttc[12][16];


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// TARGET
	Int_t adc_high_target[256]; 	
	Int_t adc_low_target[256]; 		
	Int_t tdc_le_target[256][16];      	  	
	Int_t tdc_te_target[256][16]; 		

	// SFT
	Int_t ADC_High_sft[128];
	Int_t ADC_Low_sft[128];
	Int_t TDC_LE_sft[128][16];
	Int_t TDC_TE_sft[128][16];

	Int_t HR_tdc[96];
	Int_t TKO_adc[96];
	Int_t V792_adc[256];
	Int_t VT48_tdc[256][16];

	// TOF1
	Int_t ADC_tof1U[12];	Int_t TDC_tof1U[12];
	Int_t ADC_tof1D[12];	Int_t TDC_tof1D[12];

	// TOF2
	Int_t ADC_tof2AO[12];	Int_t TDC_tof2AO[12];
	Int_t ADC_tof2AI[12];	Int_t TDC_tof2AI[12];
	Int_t ADC_tof2BO[12];	Int_t TDC_tof2BO[12];
	Int_t ADC_tof2BI[12];	Int_t TDC_tof2BI[12];

	// AC
	Int_t ADC_acu[12];
	Int_t ADC_acd[12];
	Int_t TDC_acu[12][16];
	Int_t TDC_acd[12][16];

	// GV
	Int_t ADC_gv[12];
	Int_t TDC_gv[12][16];

	// TTC
	Int_t ADC_ttc[12];
	Int_t TDC_ttc[12][16];

	// PGC
	Int_t ADC_pgc[96];
	Int_t TDC_pgc[96][16];

	Int_t ADC_pgc_Gap[12][16];
	Int_t TDC_pgc_Gap[12][8][16];

	// Ck
	Int_t TDC_ck[14][16];
	Int_t TDC_ck_selected[14] = {0};
	double TDC_ck_sum = 0.;


	// Cpi
	Int_t TDC_cpi[14][16];

	// MWPC
	Int_t MwpcADC[512];	
	Int_t ADC_C2XR[56];		Int_t ADC_C2XL[56];		Int_t ADC_C2YR[16];		Int_t ADC_C2YL[16];
	Int_t ADC_C3XR[64];		Int_t ADC_C3XL[64];		Int_t ADC_C3YR[16];		Int_t ADC_C3YL[16];
	Int_t ADC_C4XR[72];		Int_t ADC_C4XL[72];		Int_t ADC_C4YR[16];		Int_t ADC_C4YL[16];

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	char path_input[200];                   //char file_mapping[200];
	sprintf(path_input,"%s",path_merged);          
	//sprintf(file_mapping,path_mapping);
	//sprintf(path_input,"/triumfcs/trshare/trek/E36/Data/April_2015/root");

	char Name_finput[200];
	sprintf(Name_finput,"%s/Run%dMS.root",path_input, run_number);

	char path_output[200];
	sprintf(path_output,"%s",path_histo);

	TFile *foutput;
	char Name_foutput[200];
	sprintf(Name_foutput,"%s/Hist_Run%dM.root",path_output, run_number);
	foutput = new TFile(Name_foutput,"RECREATE");

	//char footer[100];
	//sprintf(footer,"Event_Display_MS.C  --  Run %d ; Event %d",Run_Number,ievt);

	cout << "   " << endl;
	cout << Name_finput << endl;

	TChain *fChain= new TChain("Tree");		
	fChain->Add(Name_finput);		
	fChain->SetMakeClass(1);							

	//fChain->SetBranchAddress("ADC_TOF1",ADC_tof1);
	//fChain->SetBranchAddress("ADC_TOF2",ADC_tof2);
	//fChain->SetBranchAddress("TDC_TOF1",TDC_tof1);
	//fChain->SetBranchAddress("TDC_TOF2",TDC_tof2);


	//fChain->SetBranchAddress("ADC_TOF",ADC_tof);
	//fChain->SetBranchAddress("TDC_TOF",TDC_tof);

	//fChain->SetBranchAddress("AC_ADC",ADC_ac);

	///////////////////////////////////////////////////////


	// TARGET
	fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);		
	fChain->SetBranchAddress("ADC_Low_TARGET",adc_low_target);
	fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);		
	fChain->SetBranchAddress("TDC_TE_TARGET",tdc_te_target);

	// SFT
	fChain->SetBranchAddress("ADC_High_SFT",ADC_High_sft);
	fChain->SetBranchAddress("ADC_Low_SFT",ADC_Low_sft);
	fChain->SetBranchAddress("TDC_LE_SFT",TDC_LE_sft);
	fChain->SetBranchAddress("TDC_TE_SFT",TDC_TE_sft);

	// TOF1
	fChain->SetBranchAddress("ADC_TOF1U",ADC_tof1U);
	fChain->SetBranchAddress("ADC_TOF1D",ADC_tof1D);
	fChain->SetBranchAddress("TDC_TOF1U",TDC_tof1U);
	fChain->SetBranchAddress("TDC_TOF1D",TDC_tof1D);

	// TOF2
	fChain->SetBranchAddress("ADC_TOF2AO",ADC_tof2AO);
	fChain->SetBranchAddress("ADC_TOF2AI",ADC_tof2AI);
	fChain->SetBranchAddress("ADC_TOF2BO",ADC_tof2BO);
	fChain->SetBranchAddress("ADC_TOF2BI",ADC_tof2BI);	
	fChain->SetBranchAddress("TDC_TOF2AO",TDC_tof2AO);
	fChain->SetBranchAddress("TDC_TOF2AI",TDC_tof2AI);
	fChain->SetBranchAddress("TDC_TOF2BO",TDC_tof2BO);
	fChain->SetBranchAddress("TDC_TOF2BI",TDC_tof2BI);	

	// AC
	fChain->SetBranchAddress("ADC_ACU",ADC_acu);
	fChain->SetBranchAddress("ADC_ACD",ADC_acd);
	fChain->SetBranchAddress("TDC_ACU",TDC_acu);
	fChain->SetBranchAddress("TDC_ACD",TDC_acd);

	// GV
	fChain->SetBranchAddress("ADC_GV",ADC_gv);
	fChain->SetBranchAddress("TDC_GV",TDC_gv);

	// TTC
	fChain->SetBranchAddress("ADC_TTC",ADC_ttc);
	fChain->SetBranchAddress("TDC_TTC",TDC_ttc);

	// ALL
	fChain->SetBranchAddress("HR_TDC", HR_tdc);
	fChain->SetBranchAddress("TKO_ADC", TKO_adc);
	fChain->SetBranchAddress("V792_ADC", V792_adc);
	fChain->SetBranchAddress("VT48_TDC", VT48_tdc);

	// PGC
	fChain->SetBranchAddress("ADC_PGC",ADC_pgc);
	fChain->SetBranchAddress("TDC_PGC",TDC_pgc);
	fChain->SetBranchAddress("ADC_PGC_Gap",ADC_pgc_Gap);
	fChain->SetBranchAddress("TDC_PGC_Gap",TDC_pgc_Gap);

	// Ck
	fChain->SetBranchAddress("TDC_Ck",TDC_ck);

	// Cpi
	fChain->SetBranchAddress("TDC_Cpi",TDC_cpi);

	// MWPC
	fChain->SetBranchAddress("MWPCADC",MwpcADC);
	fChain->SetBranchAddress("ADC_C2X_R", ADC_C2XR);
	fChain->SetBranchAddress("ADC_C2X_L", ADC_C2XL);
	fChain->SetBranchAddress("ADC_C2Y_R", ADC_C2YR);
	fChain->SetBranchAddress("ADC_C2Y_L", ADC_C2YL);
	fChain->SetBranchAddress("ADC_C3X_R", ADC_C3XR);
	fChain->SetBranchAddress("ADC_C3X_L", ADC_C3XL);
	fChain->SetBranchAddress("ADC_C3Y_R", ADC_C3YR);
	fChain->SetBranchAddress("ADC_C3Y_L", ADC_C3YL);
	fChain->SetBranchAddress("ADC_C4X_R", ADC_C4XR);
	fChain->SetBranchAddress("ADC_C4X_L", ADC_C4XL);
	fChain->SetBranchAddress("ADC_C4Y_R", ADC_C4YR);
	fChain->SetBranchAddress("ADC_C4Y_L", ADC_C4YL);

	///////////////////////////////////////////////////////////////////////////

	// TARGET
	TH1D *h_ADC_High_TARGET[256];	char Title_ADC_High_TARGET[256][100];	char Name_ADC_High_TARGET[256][100];
	TH1D *h_ADC_Low_TARGET[256];    char Title_ADC_Low_TARGET[256][100];	char Name_ADC_Low_TARGET[256][100];
	TH1D *h_TDC_LE_TARGET[256];     char Title_TDC_LE_TARGET[256][100];	char Name_TDC_LE_TARGET[256][100];     
	TH1D *h_TDC_TE_TARGET[256];		char Title_TDC_TE_TARGET[256][100];	char Name_TDC_TE_TARGET[256][100]; 

	// SFT
	TH1D *h_ADC_High_SFT[128];	char Title_ADC_High_SFT[128][100];	char Name_ADC_High_SFT[128][100];
	TH1D *h_ADC_Low_SFT[128];   char Title_ADC_Low_SFT[128][100];	char Name_ADC_Low_SFT[128][100];
	TH1D *h_TDC_LE_SFT[128];    char Title_TDC_LE_SFT[128][100];	char Name_TDC_LE_SFT[128][100];     
	TH1D *h_TDC_TE_SFT[128];	char Title_TDC_TE_SFT[128][100];	char Name_TDC_TE_SFT[128][100]; 

	// MWPCs
	TH1D *h_MWPCADC[512];	char Title_MWPCADC[512][100];	char Name_MWPCADC[512][100];  
	TH1D *h_ADC_C2XR[56];	char Title_ADC_C2XR[56][100];	char Name_ADC_C2XR[56][100];
	TH1D *h_ADC_C2XL[56];	char Title_ADC_C2XL[56][100];	char Name_ADC_C2XL[56][100];
	TH1D *h_ADC_C2YR[16];	char Title_ADC_C2YR[16][100];	char Name_ADC_C2YR[16][100];
	TH1D *h_ADC_C2YL[16];	char Title_ADC_C2YL[16][100];	char Name_ADC_C2YL[16][100];

	TH1D *h_ADC_C3XR[64];	char Title_ADC_C3XR[64][100];	char Name_ADC_C3XR[64][100];
	TH1D *h_ADC_C3XL[64];	char Title_ADC_C3XL[64][100];	char Name_ADC_C3XL[64][100];
	TH1D *h_ADC_C3YR[16];	char Title_ADC_C3YR[16][100];	char Name_ADC_C3YR[16][100];
	TH1D *h_ADC_C3YL[16];	char Title_ADC_C3YL[16][100];	char Name_ADC_C3YL[16][100];

	TH1D *h_ADC_C4XR[72];	char Title_ADC_C4XR[72][100];	char Name_ADC_C4XR[72][100];
	TH1D *h_ADC_C4XL[72];	char Title_ADC_C4XL[72][100];	char Name_ADC_C4XL[72][100];
	TH1D *h_ADC_C4YR[16];	char Title_ADC_C4YR[16][100];	char Name_ADC_C4YR[16][100];
	TH1D *h_ADC_C4YL[16];	char Title_ADC_C4YL[16][100];	char Name_ADC_C4YL[16][100];

	// TOF1
	TH1D *h_ADC_tof1U[12];   char Title_ADC_tof1U[12][100];	char Name_ADC_tof1U[12][100];
	TH1D *h_ADC_tof1D[12];   char Title_ADC_tof1D[12][100];	char Name_ADC_tof1D[12][100];
	TH1D *h_TDC_tof1U[12];   char Title_TDC_tof1U[12][100];	char Name_TDC_tof1U[12][100];
	TH1D *h_TDC_tof1D[12];   char Title_TDC_tof1D[12][100];	char Name_TDC_tof1D[12][100];

	// TOF2
	TH1D *h_ADC_tof2AO[12];   char Title_ADC_tof2AO[12][100];	char Name_ADC_tof2AO[12][100];
	TH1D *h_ADC_tof2AI[12];   char Title_ADC_tof2AI[12][100];	char Name_ADC_tof2AI[12][100];
	TH1D *h_ADC_tof2BO[12];   char Title_ADC_tof2BO[12][100];	char Name_ADC_tof2BO[12][100];
	TH1D *h_ADC_tof2BI[12];   char Title_ADC_tof2BI[12][100];	char Name_ADC_tof2BI[12][100];
	TH1D *h_TDC_tof2AO[12];   char Title_TDC_tof2AO[12][100];	char Name_TDC_tof2AO[12][100];
	TH1D *h_TDC_tof2AI[12];   char Title_TDC_tof2AI[12][100];	char Name_TDC_tof2AI[12][100];
	TH1D *h_TDC_tof2BO[12];   char Title_TDC_tof2BO[12][100];	char Name_TDC_tof2BO[12][100];
	TH1D *h_TDC_tof2BI[12];   char Title_TDC_tof2BI[12][100];	char Name_TDC_tof2BI[12][100];

	// AC
	TH1D *h_ADC_acu[12];	char Title_ADC_acu[12][100];	char Name_ADC_acu[12][100];
	TH1D *h_ADC_acd[12];	char Title_ADC_acd[12][100];	char Name_ADC_acd[12][100];
	TH1D *h_TDC_acu[12];	char Title_TDC_acu[12][100];	char Name_TDC_acu[12][100];
	TH1D *h_TDC_acd[12];	char Title_TDC_acd[12][100];	char Name_TDC_acd[12][100];

	// GV
	TH1D *h_ADC_gv[12];	char Title_ADC_gv[12][100];	char Name_ADC_gv[12][100];
	TH1D *h_TDC_gv[12];	char Title_TDC_gv[12][100];	char Name_TDC_gv[12][100];

	// TTC
	TH1D *h_ADC_ttc[12];	char Title_TDC_ttc[12][100];	char Name_ADC_ttc[12][100];
	TH1D *h_TDC_ttc[12];	char Title_ADC_ttc[12][100];	char Name_TDC_ttc[12][100];

	// PGC
	//TH1D *h_ADC_pgc[96];	 	char Title_ADC_pgc[96][100];		char Name_ADC_pgc[96][100];
	//TH1D *h_TDC_pgc[96];		char Title_TDC_pgc[96][100];		char Name_TDC_pgc[96][100];

	TH1D *h_ADC_pgc_sorted[12][8];	 	char Title_ADC_pgc_sorted[12][8][100];		char Name_ADC_pgc_sorted[12][8][100];
	TH1D *h_TDC_pgc_sorted[12][8];		char Title_TDC_pgc_sorted[12][8][100];		char Name_TDC_pgc_sorted[12][8][100];

	TH1D *h_ADC_pgc_Gap[12][8];	 		char Title_ADC_pgc_Gap[12][8][100];		char Name_ADC_pgc_Gap[12][8][100];
	TH1D *h_TDC_pgc_Gap[12][8];			char Title_TDC_pgc_Gap[12][8][100];		char Name_TDC_pgc_Gap[12][8][100];

	// Ck
	TH1D *h_TDC_Ck[14];		char Title_TDC_Ck[14][100];		char Name_TDC_Ck[14][100];	
	TH1D *h_TDC_Ck_Avg;     char Title_TDC_Ck_Avg[100];		char Name_TDC_Ck_Avg[100];
	vector <double> vec_Ck_Avg;
	vec_Ck_Avg.clear();

	// Cpi
	TH1D *h_TDC_Cpi[14];		char Title_TDC_Cpi[14][100];	char Name_TDC_Cpi[14][100];

	// ALL
	TH1D *h_HR_tdc[96];		char Title_HR_tdc[96][100];		char Name_HR_tdc[96][100];
	TH1D *h_TKO_adc[96];	char Title_TKO_adc[96][100];	char Name_TKO_adc[96][100];
	TH1D *h_V792_adc[256];	char Title_V792_adc[256][100];	char Name_V792_adc[256][100];
	TH1D *h_VT48_tdc[256];	char Title_VT48_tdc[256][100];	char Name_VT48_tdc[256][100];


	/////////////////////////////////////////////////////////////
	// Setting Histogram Names and Titles
	
	sprintf(Title_TDC_Ck_Avg,"Ck | TDC Average (Run %d)", run_number);
	sprintf(Name_TDC_Ck_Avg,"TDC Ck Average");

	///// 512 Channels
	// MWPCs
	for (int i=0; i<512; i++) {
		sprintf(Title_MWPCADC[i],"Raw ADC (Ch. %d)  --  MWPC  |  Run %d",i, run_number);
		sprintf(Name_MWPCADC[i],"ADC_MWPC (Ch. %d)",i);
	}


	///// 256 Channels
	for(int i=0; i<256; i++){
		// TARGET	
		sprintf(Title_ADC_High_TARGET[i],"Raw ADC High Gain (Ch. %d)  --  TARGET  |  Run %d",i, run_number); 
		sprintf(Title_ADC_Low_TARGET[i],"Raw ADC Low Gain (Ch. %d)  --  TARGET  |  Run %d",i, run_number); 
		sprintf(Title_TDC_LE_TARGET[i],"Raw TDC (LE) (Ch. %d)  --  TARGET  |  Run %d",i, run_number); 
		sprintf(Title_TDC_TE_TARGET[i],"Raw TDC (TE) (Ch. %d)  --  TARGET  |  Run %d",i, run_number); 

		sprintf(Name_ADC_High_TARGET[i],"ADC_High (Ch. %d) - TARGET",i);	
		sprintf(Name_ADC_Low_TARGET[i],"ADC_Low (Ch. %d) - TARGET",i);
		sprintf(Name_TDC_LE_TARGET[i],"TDC_LE (Ch. %d) - TARGET",i);
		sprintf(Name_TDC_TE_TARGET[i],"TDC_TE (Ch. %d) - TARGET",i);

		// V792
		sprintf(Title_V792_adc[i], "Raw ADC  (Ch. %d) -- V792 ADC  |  Run %d", i, run_number);

		sprintf(Name_V792_adc[i], "V792_ADC (Ch. %d)", i);

		// VT48
		sprintf(Title_VT48_tdc[i], "Raw TDC  (Ch. %d) -- VT48 TDC  |  Run %d", i, run_number);

		sprintf(Name_VT48_tdc[i], "VT48_TDC (Ch. %d)", i);
	}

	///// 128 Channels
	for(int i=0; i<128; i++){
		// SFT
		sprintf(Title_ADC_High_SFT[i],"Raw ADC High Gain (Ch. %d)  --  SFT  |  Run %d",i, run_number); 
		sprintf(Title_ADC_Low_SFT[i],"Raw ADC Low Gain (Ch. %d)  --  SFT  |  Run %d",i, run_number); 
		sprintf(Title_TDC_LE_SFT[i],"Raw TDC (LE) (Ch. %d)  --  SFT  |  Run %d",i, run_number); 
		sprintf(Title_TDC_TE_SFT[i],"Raw TDC (TE) (Ch. %d)  --  SFT  |  Run %d",i, run_number); 

		sprintf(Name_ADC_High_SFT[i],"ADC_High (Ch. %d) - SFT",i);	
		sprintf(Name_ADC_Low_SFT[i],"ADC_Low (Ch. %d) - SFT",i);
		sprintf(Name_TDC_LE_SFT[i],"TDC_LE (Ch. %d) - SFT",i);
		sprintf(Name_TDC_TE_SFT[i],"TDC_TE (Ch. %d) - SFT",i);
	}

	///// 96 Channels
	for (int i=0; i<96; i++){
		// TOF
		//sprintf(Title_ADC_tof[i], "Raw ADC (Ch. %d) -- TOF", i);
		//sprintf(Name_ADC_tof[i], "ADC_TOF (Ch. %d)", i);

		//sprintf(Title_TDC_tof[i], "Raw TDC (Ch. %d) -- TOF", i);
		//sprintf(Name_TDC_tof[i], "TDC_TOF (Ch. %d)", i);
	
		// PGC
		//sprintf(Title_ADC_pgc[i], "Raw ADC (Ch. %d) -- PGC", i);
		//sprintf(Name_ADC_pgc[i], "ADC_PGC (Ch. %d)", i);
		//sprintf(Title_TDC_pgc[i], "Raw TDC (Ch. %d) -- PGC", i);
		//sprintf(Name_TDC_pgc[i], "TDC_PGC (Ch. %d)", i);

		// HR TDC
		sprintf(Title_HR_tdc[i], "Raw TDC  (Ch. %d) -- HR TDC  |  Run %d", i, run_number);
		sprintf(Name_HR_tdc[i], "HR_TDC (Ch. %d)", i);

		//TKO ADC
		sprintf(Title_TKO_adc[i], "Raw ADC  (Ch. %d) -- TKO ADC  |  Run %d", i, run_number);
		sprintf(Name_TKO_adc[i], "TKO_ADC (Ch. %d)", i);
	}

	///// 72 Channels
	for(int i=0; i<72; i++){
		sprintf(Title_ADC_C4XR[i] , "Raw ADC (MWPC Ch. %d) -- C4X %d (1 - 6)  |  Run %d", C4XR_ADC[i], i+1, run_number);
		sprintf(Name_ADC_C4XR[i], "ADC_C4X%d_R", i+1);

		sprintf(Title_ADC_C4XL[i] , "Raw ADC (MWPC Ch. %d) -- C4X %d (7 - 12)  |  Run %d", C4XL_ADC[i], i+1, run_number);
		sprintf(Name_ADC_C4XL[i], "ADC_C4X%d_L", i+1);
	}

	///// 64 Channels
	for(int i=0; i<64; i++){
		sprintf(Title_ADC_C3XR[i] , "Raw ADC (MWPC Ch. %d) -- C3X %d (1 - 6)  |  Run %d", C3XR_ADC[i], i+1, run_number);
		sprintf(Name_ADC_C3XR[i], "ADC_C3X%d_R", i+1);

		sprintf(Title_ADC_C3XL[i] , "Raw ADC (MWPC Ch. %d) -- C3X %d (7 - 12)  |  Run %d", C3XL_ADC[i], i+1, run_number);
		sprintf(Name_ADC_C3XL[i], "ADC_C3X%d_L", i+1);
	}

	///// 56 Channels
	for(int i=0; i<56; i++){
		sprintf(Title_ADC_C2XR[i] , "Raw ADC (MWPC Ch. %d) -- C2X %d (1 - 6)  |  Run %d", C2XR_ADC[i], i+1, run_number);
		sprintf(Name_ADC_C2XR[i], "ADC_C2X%d_R", i+1);

		sprintf(Title_ADC_C2XL[i] , "Raw ADC (MWPC Ch. %d) -- C2X %d (7 - 12)  |  Run %d", C2XL_ADC[i], i+1, run_number);
		sprintf(Name_ADC_C2XL[i], "ADC_C2X%d_L", i+1);
	}

	///// 16 Channels
	for(int i=0; i<16; i++){
		sprintf(Title_ADC_C2YR[i] , "Raw ADC (MWPC Ch. %d) -- C2Y %d (1 - 6)  |  Run %d", C2YR_ADC[i], i+1, run_number);
		sprintf(Name_ADC_C2YR[i], "ADC_C2Y%d_R", i+1);

		sprintf(Title_ADC_C2YL[i] , "Raw ADC (MWPC Ch. %d) -- C2Y %d (7 - 12)  |  Run %d", C2YL_ADC[i], i+1, run_number);
		sprintf(Name_ADC_C2YL[i], "ADC_C2Y%d_L", i+1);

		sprintf(Title_ADC_C3YR[i] , "Raw ADC (MWPC Ch. %d) -- C3Y %d (1 - 6)  |  Run %d", C3YR_ADC[i], i+1, run_number);
		sprintf(Name_ADC_C3YR[i], "ADC_C3Y%d_R", i+1);

		sprintf(Title_ADC_C3YL[i] , "Raw ADC (MWPC Ch. %d) -- C3Y %d (7 - 12)  |  Run %d", C3YL_ADC[i], i+1, run_number);
		sprintf(Name_ADC_C3YL[i], "ADC_C3Y%d_L", i+1);

		sprintf(Title_ADC_C4YR[i] , "Raw ADC (MWPC Ch. %d) -- C4Y %d (1 - 6)  |  Run %d", C4YR_ADC[i], i+1, run_number);
		sprintf(Name_ADC_C4YR[i], "ADC_C4Y%d_R", i+1);

		sprintf(Title_ADC_C4YL[i] , "Raw ADC (MWPC Ch. %d) -- C4Y %d (7 - 12)  |  Run %d", C4YL_ADC[i], i+1, run_number);
		sprintf(Name_ADC_C4YL[i], "ADC_C4Y%d_L", i+1);
	}

	///// 56 Channels
	//for (int i=0; i<56; i++) {
		// TOF2
		//sprintf(Title_ADC_tof2[i],"Raw ADC (Ch. %d)  --  TOF2",i);
		//sprintf(Name_ADC_tof2[i],"ADC_TOF2 (Ch. %d)",i);

		//sprintf(Title_TDC_tof2[i],"Raw TDC (Ch. %d)  --  TOF2",i);
		//sprintf(Name_TDC_tof2[i],"TDC_TOF2 (Ch. %d)",i);
	//}

	///// 24 Channels
	//for (int i=0; i<24; i++) {
		// TOF1
	//	sprintf(Title_ADC_tof1[i],"Raw ADC (Ch. %d)  --  TOF1",i);
	//	sprintf(Name_ADC_tof1[i],"ADC_TOF1 (Ch. %d)",i);

	//	sprintf(Title_ADC_ac[i],"Raw ADC (Ch. %d)  --  AC",i);
	//	sprintf(Name_ADC_ac[i],"ADC_AC (Ch. %d)",i);

	//	sprintf(Title_TDC_tof1[i],"Raw TDC (Ch. %d)  --  TOF1",i);
	//	sprintf(Name_TDC_tof1[i],"TDC_TOF1 (Ch. %d)",i);
	//}

	///// 12 Channels
	for (int i=0; i<12; i++) {
		// TOF1
		sprintf(Title_ADC_tof1U[i],"Raw ADC (V792 Ch. %d)  --  TOF1 (%d UP)  |  Run %d",TOF1U_ADC[i], i+1, run_number);
		sprintf(Name_ADC_tof1U[i],"ADC_TOF1 (%d UP)",i+1);
		sprintf(Title_ADC_tof1D[i],"Raw ADC (V792 Ch. %d)  --  TOF1 (%d DOWN)  |  Run %d",TOF1D_ADC[i], i+1, run_number);
		sprintf(Name_ADC_tof1D[i],"ADC_TOF1 (%d DOWN)",i+1);

		sprintf(Title_TDC_tof1U[i],"Raw TDC (HR_TDC Ch. %d)  --  TOF1 (%d UP)  |  Run %d",TOF1U_TDC[i], i+1, run_number);
		sprintf(Name_TDC_tof1U[i],"TDC_TOF1 (%d UP)",i+1);
		sprintf(Title_TDC_tof1D[i],"Raw TDC (HR _TDC Ch. %d)  --  TOF1 (%d DOWN)  |  Run %d",TOF1D_TDC[i], i+1, run_number);
		sprintf(Name_TDC_tof1D[i],"TDC_TOF1 (%d DOWN)",i+1);

		// TOF2
		sprintf(Title_ADC_tof2AO[i],"Raw ADC (V792 Ch. %d)  --  TOF2 (%d A OUT)  |  Run %d",TOF2AO_ADC[i], i+1, run_number);
		sprintf(Name_ADC_tof2AO[i],"ADC_TOF2 (%d AO)",i+1);
		sprintf(Title_ADC_tof2AI[i],"Raw ADC (V792 Ch. %d)  --  TOF2 (%d A IN)  |  Run %d",TOF2AI_ADC[i], i+1, run_number);
		sprintf(Name_ADC_tof2AI[i],"ADC_TOF2 (%d AI)",i+1);
		sprintf(Title_ADC_tof2BO[i],"Raw ADC (V792 Ch. %d)  --  TOF2 (%d B OUT)  |  Run %d",TOF2BO_ADC[i], i+1, run_number);
		sprintf(Name_ADC_tof2BO[i],"ADC_TOF2 (%d BO)",i+1);
		sprintf(Title_ADC_tof2BI[i],"Raw ADC (V792 Ch. %d)  --  TOF2 (%d B IN)  |  Run %d",TOF2BI_ADC[i], i+1, run_number);
		sprintf(Name_ADC_tof2BI[i],"ADC_TOF2 (%d BI)",i+1);

		sprintf(Title_TDC_tof2AO[i],"Raw TDC (HR_TDC Ch. %d)  --  TOF2 (%d A OUT)  |  Run %d",TOF2AO_TDC[i], i+1, run_number);
		sprintf(Name_TDC_tof2AO[i],"TDC_TOF2 (%d AO)",i+1);
		sprintf(Title_TDC_tof2AI[i],"Raw TDC (HR_TDC Ch. %d)  --  TOF2 (%d A IN)  |  Run %d",TOF2AI_TDC[i], i+1, run_number);
		sprintf(Name_TDC_tof2AI[i],"TDC_TOF2 (%d AI)",i+1);
		sprintf(Title_TDC_tof2BO[i],"Raw TDC (HR_TDC Ch. %d)  --  TOF2 (%d B OUT)  |  Run %d",TOF2BO_TDC[i], i+1, run_number);
		sprintf(Name_TDC_tof2BO[i],"TDC_TOF2 (%d BO)",i+1);
		sprintf(Title_TDC_tof2BI[i],"Raw TDC (HR_TDC Ch. %d)  --  TOF2 (%d B IN)  |  Run %d",TOF2BI_TDC[i], i+1, run_number);
		sprintf(Name_TDC_tof2BI[i],"TDC_TOF2 (%d BI)",i+1);

		// AC
		sprintf(Title_ADC_acu[i],"Raw ADC (TKO_ADC Ch. %d)  --  AC (%d UP)  |  Run %d",ACU_ADC[i],i+1, run_number);
		sprintf(Name_ADC_acu[i],"ADC_AC (%d UP)",i+1);
		sprintf(Title_ADC_acd[i],"Raw ADC (TKO_ADC Ch. %d)  --  AC (%d DOWN)  |  Run %d",ACD_ADC[i],i+1, run_number);
		sprintf(Name_ADC_acd[i],"ADC_AC (%d DOWN)",i+1);
	
		sprintf(Title_TDC_acu[i],"Raw TDC (VT48_TDC Ch. %d)  --  AC (%d UP)  |  Run %d",ACU_TDC[i],i+1, run_number);
		sprintf(Name_TDC_acu[i],"TDC_AC (%d UP)",i+1);
		sprintf(Title_TDC_acd[i],"Raw TDC (VT48_TDC Ch. %d)  --  AC (%d DOWN)  |  Run %d",ACD_TDC[i],i+1, run_number);
		sprintf(Name_TDC_acd[i],"TDC_AC (%d DOWN)",i+1);

		// GV
		sprintf(Title_ADC_gv[i],"Raw ADC (V792 Ch. %d)  --  GV (Counter %d)  |  Run %d",GV_ADC[i],i+1, run_number);
		sprintf(Name_ADC_gv[i],"ADC_GV (Counter %d)",i+1);

		sprintf(Title_TDC_gv[i],"Raw TDC (VT48_TDC Ch. %d)  --  GV (Counter %d)  |  Run %d",GV_TDC[i],i+1, run_number);
		sprintf(Name_TDC_gv[i],"TDC_GV (Counter %d)",i+1);
	
		// TTC
		sprintf(Title_ADC_ttc[i],"Raw ADC (TKO_ADC Ch. %d)  --  TTC (Counter %d)  |  Run %d",TTC_ADC[i],i+1, run_number);
		sprintf(Name_ADC_ttc[i],"ADC_TTC (Counter %d)",i+1);

		sprintf(Title_TDC_ttc[i],"Raw TDC (VT48_TDC Ch. %d)  --  TTC (Counter %d)  |  Run %d",TTC_TDC[i],i+1, run_number);
		sprintf(Name_TDC_ttc[i],"TDC_TTC (Counter %d)",i+1);

		// PGC
		for (int j=0; j<8; j++){
			sprintf(Title_ADC_pgc_sorted[i][j],"Raw ADC (V792 Ch. %d)  --  PGC (Counter %d)  |  Run %d",PGC_ADC[8*i+j],i+1, run_number);
			sprintf(Name_ADC_pgc_sorted[i][j],"ADC_PGC (Gap %d-%d)",i+1,j+1);

			sprintf(Title_TDC_pgc_sorted[i][j],"Raw TDC (VT48_TDC Ch. %d)  --  PGC (Counter %d)  |  Run %d",PGC_TDC[8*i+j],i+1, run_number);
			sprintf(Name_TDC_pgc_sorted[i][j],"TDC_PGC (Gap %d-%d)",i+1,j+1);

			sprintf(Title_ADC_pgc_Gap[i][j],"Raw ADC (V792 Ch. %d)  --  PGC_Gap (Counter %d)  |  Run %d",PGC_ADC[8*i+j],i+1, run_number);
			sprintf(Name_ADC_pgc_Gap[i][j],"ADC_PGC_Gap (Gap %d-%d)",i+1,j+1);

			sprintf(Title_TDC_pgc_Gap[i][j],"Raw TDC (VT48_TDC Ch. %d)  --  PGC_Gap (Counter %d)  |  Run %d",PGC_TDC[8*i+j],i+1, run_number);
			sprintf(Name_TDC_pgc_Gap[i][j],"TDC_PGC_Gap (Gap %d-%d)",i+1,j+1);
		}
	}

	for(int i=0; i<14; i++){	
		// Ck
		sprintf(Title_TDC_Ck[i],"Raw TDC (VT48_TDC Ch. %d)  --  Ck (Counter %d)  |  Run %d",VT48Ck[i],i+1, run_number);
		sprintf(Name_TDC_Ck[i],"TDC_Ck (Counter %d)",i+1);

		// Cpi
		sprintf(Title_TDC_Cpi[i],"Raw TDC (VT48_TDC Ch. %d)  --  Cpi (Counter %d)  |  Run %d",VT48Cpi[i],i+1, run_number);
		sprintf(Name_TDC_Cpi[i],"TDC_Cpi (Counter %d)",i+1);
	}


	//////////////////////////////////////////////////////////////////
	// CREATING HISTOGRAMS

	///// 512 Channels
	for (int  i=0; i<512; i++) {
		// MWPCs
		h_MWPCADC[i] = new TH1D(Name_MWPCADC[i],Title_MWPCADC[i],512,0,4096);
	}

	///// 256 Channels
	for (int i=0; i<256; i++){
		// TARGET
		h_ADC_High_TARGET[i] = new TH1D(Name_ADC_High_TARGET[i],Title_ADC_High_TARGET[i],512,0,4096);
		h_ADC_Low_TARGET[i] = new TH1D(Name_ADC_Low_TARGET[i],Title_ADC_Low_TARGET[i],512,0,4096);
		h_TDC_LE_TARGET[i] = new TH1D(Name_TDC_LE_TARGET[i],Title_TDC_LE_TARGET[i],512,0,4096);
		h_TDC_TE_TARGET[i] = new TH1D(Name_TDC_TE_TARGET[i],Title_TDC_TE_TARGET[i],512,0,4096);

		// V792
		h_V792_adc[i] = new TH1D(Name_V792_adc[i], Title_V792_adc[i], 512, 0, 4096);

		// VT48
		h_VT48_tdc[i] = new TH1D(Name_VT48_tdc[i], Title_VT48_tdc[i], 512, 0, 4096);
	}

	///// 128 Channels	
	for (int i=0; i<128; i++){
		// SFT
		h_ADC_High_SFT[i] = new TH1D(Name_ADC_High_SFT[i],Title_ADC_High_SFT[i],512,0,4096);
		h_ADC_Low_SFT[i] = new TH1D(Name_ADC_Low_SFT[i],Title_ADC_Low_SFT[i],512,0,4096);
		h_TDC_LE_SFT[i] = new TH1D(Name_TDC_LE_SFT[i],Title_TDC_LE_SFT[i],512,0,4096);
		h_TDC_TE_SFT[i] = new TH1D(Name_TDC_TE_SFT[i],Title_TDC_TE_SFT[i],512,0,4096);
	}

	///// 96 Channels

		// PGC
	for (int i=0; i<96; i++){
		// TOF
		//h_ADC_tof[i] = new TH1D(Name_ADC_tof[i], Title_ADC_tof[i], 512, 0, 4096);
		//h_TDC_tof[i] = new TH1D(Name_TDC_tof[i], Title_TDC_tof[i], 512, 0, 4096);

		// PGC
		//h_ADC_pgc[i] = new TH1D(Name_ADC_pgc[i], Title_ADC_pgc[i], 1000, 0, 4096);
		//h_TDC_pgc[i] = new TH1D(Name_TDC_pgc[i], Title_TDC_pgc[i], 1000, 0, 4096);

		// HR TDC
		h_HR_tdc[i] = new TH1D(Name_HR_tdc[i], Title_HR_tdc[i], 512, 0, 4096);

		// TKO ADC
		h_TKO_adc[i] = new TH1D(Name_TKO_adc[i], Title_TKO_adc[i], 512, 0, 4096);
	}

	///// 72 Channels	
	for(int i=0; i<72; i++){
		h_ADC_C4XR[i] = new TH1D(Name_ADC_C4XR[i], Title_ADC_C4XR[i], 512, 0, 4096);
		h_ADC_C4XL[i] = new TH1D(Name_ADC_C4XL[i], Title_ADC_C4XL[i], 512, 0, 4096);
	}

	///// 64 Channels	
	for(int i=0; i<64; i++){
		h_ADC_C3XR[i] = new TH1D(Name_ADC_C3XR[i], Title_ADC_C3XR[i], 512, 0, 4096);
		h_ADC_C3XL[i] = new TH1D(Name_ADC_C3XL[i], Title_ADC_C3XL[i], 512, 0, 4096);
	}

	///// 56 Channels	
	for(int i=0; i<56; i++){
		h_ADC_C2XR[i] = new TH1D(Name_ADC_C2XR[i], Title_ADC_C2XR[i], 512, 0, 4096);
		h_ADC_C2XL[i] = new TH1D(Name_ADC_C2XL[i], Title_ADC_C2XL[i], 512, 0, 4096);
	}

	///// 16 Channels	
	for(int i=0; i<16; i++){
    	h_ADC_C2YR[i] = new TH1D(Name_ADC_C2YR[i], Title_ADC_C2YR[i], 512, 0, 4096);
		h_ADC_C2YL[i] = new TH1D(Name_ADC_C2YL[i], Title_ADC_C2YL[i], 512, 0, 4096);

		h_ADC_C3YR[i] = new TH1D(Name_ADC_C3YR[i], Title_ADC_C3YR[i], 512, 0, 4096);
		h_ADC_C3YL[i] = new TH1D(Name_ADC_C3YL[i], Title_ADC_C3YL[i], 512, 0, 4096);

		h_ADC_C4YR[i] = new TH1D(Name_ADC_C4YR[i], Title_ADC_C4YR[i], 512, 0, 4096);
		h_ADC_C4YL[i] = new TH1D(Name_ADC_C4YL[i], Title_ADC_C4YL[i], 512, 0, 4096);
	}

	h_TDC_Ck_Avg = new TH1D(Name_TDC_Ck_Avg, Title_TDC_Ck_Avg, 250, 0, 2000);

	///// 56 Channels
	//for (int  i=0; i<56; i++) {
		// TOF2
	//	h_ADC_tof2[i] = new TH1D(Name_ADC_tof2[i],Title_ADC_tof2[i],512,0,4096);
	//	h_TDC_tof2[i] = new TH1D(Name_TDC_tof2[i],Title_TDC_tof2[i],512,0,4096);
	//}

	///// 24 Channels
	//for (int  i=0; i<24; i++) {
		// TOF1
	//	h_ADC_tof1[i] = new TH1D(Name_ADC_tof1[i],Title_ADC_tof1[i],512,0,4096);
	//	h_TDC_tof1[i] = new TH1D(Name_TDC_tof1[i],Title_TDC_tof1[i],512,0,4096);

		// AC
	//	h_ADC_ac[i] = new TH1D(Name_ADC_ac[i],Title_ADC_ac[i],512,0,4096);
	//}

	///// 12 Channels
	for (int  i=0; i<12; i++) {
		// TOF1
		h_ADC_tof1U[i] = new TH1D(Name_ADC_tof1U[i], Title_ADC_tof1U[i], 512, 0, 4096);
		h_ADC_tof1D[i] = new TH1D(Name_ADC_tof1D[i], Title_ADC_tof1D[i], 512, 0, 4096);
		h_TDC_tof1U[i] = new TH1D(Name_TDC_tof1U[i],Title_TDC_tof1U[i],512,0,4096);
		h_TDC_tof1D[i] = new TH1D(Name_TDC_tof1D[i],Title_TDC_tof1D[i],512,0,4096);

		// TOF2
		h_ADC_tof2AO[i] = new TH1D(Name_ADC_tof2AO[i],Title_ADC_tof2AO[i],512,0,4096);
		h_ADC_tof2AI[i] = new TH1D(Name_ADC_tof2AI[i],Title_ADC_tof2AI[i],512,0,4096);
		h_ADC_tof2BO[i] = new TH1D(Name_ADC_tof2BO[i],Title_ADC_tof2BO[i],512,0,4096);
		h_ADC_tof2BI[i] = new TH1D(Name_ADC_tof2BI[i],Title_ADC_tof2BI[i],512,0,4096);
		h_TDC_tof2AO[i] = new TH1D(Name_TDC_tof2AO[i],Title_TDC_tof2AO[i],512,0,4096);
		h_TDC_tof2AI[i] = new TH1D(Name_TDC_tof2AI[i],Title_TDC_tof2AI[i],512,0,4096);
		h_TDC_tof2BO[i] = new TH1D(Name_TDC_tof2BO[i],Title_TDC_tof2BO[i],512,0,4096);
		h_TDC_tof2BI[i] = new TH1D(Name_TDC_tof2BI[i],Title_TDC_tof2BI[i],512,0,4096);

		// AC
		h_ADC_acu[i] = new TH1D(Name_ADC_acu[i], Title_ADC_acu[i], 512, 0, 4096);
		h_ADC_acd[i] = new TH1D(Name_ADC_acd[i], Title_ADC_acd[i], 512, 0, 4096);
		h_TDC_acu[i] = new TH1D(Name_TDC_acu[i], Title_TDC_acu[i], 512, 0, 4096);
		h_TDC_acd[i] = new TH1D(Name_TDC_acd[i], Title_TDC_acd[i], 512, 0, 4096);

		// GV 
		h_ADC_gv[i] = new TH1D(Name_ADC_gv[i], Title_ADC_gv[i], 512, 0, 4096);
		h_TDC_gv[i] = new TH1D(Name_TDC_gv[i], Title_TDC_gv[i], 512, 0, 4096);

		// TTC
		h_ADC_ttc[i] = new TH1D(Name_ADC_ttc[i], Title_ADC_ttc[i], 512, 0, 4096);
		h_TDC_ttc[i] = new TH1D(Name_TDC_ttc[i], Title_TDC_ttc[i], 512, 0, 4096);

		// PGC
		for (int j=0; j<8; j++){
			h_ADC_pgc_sorted[i][j] = new TH1D(Name_ADC_pgc_sorted[i][j], Title_ADC_pgc_sorted[i][j], 512, 0, 4096);
			h_TDC_pgc_sorted[i][j] = new TH1D(Name_TDC_pgc_sorted[i][j], Title_TDC_pgc_sorted[i][j], 512, 0, 4096);

			h_ADC_pgc_Gap[i][j] = new TH1D(Name_ADC_pgc_Gap[i][j], Title_ADC_pgc_Gap[i][j], 512, 0, 4096);
			h_TDC_pgc_Gap[i][j] = new TH1D(Name_TDC_pgc_Gap[i][j], Title_TDC_pgc_Gap[i][j], 512, 0, 4096);
		}
	}

	for(int i=0; i<14; i++){
		h_TDC_Ck[i] = new TH1D(Name_TDC_Ck[i], Title_TDC_Ck[i], 512, 0, 4096);
		h_TDC_Cpi[i] = new TH1D(Name_TDC_Cpi[i], Title_TDC_Cpi[i], 512, 0, 4096);
	}



	// READ ALL ENTRIES 
	Int_t nentries = (Int_t)fChain->GetEntries();
	cout <<  "Total Number of Entries :     " <<  nentries << endl;

	cout << "  " << endl;
	cout << "***************  Processing....  ***************" << endl;
	cout << "  " << endl;

	if(flag!=0) nentries=flag;
	for (Int_t i=0; i<nentries; i++) {
	//for (Int_t i=0; i<21; i++) {
		fChain->GetEntry(i);

		if(QUIET==0){
			//if(nentries<=30000){
			//	if(i%1000==1) cout<<"**** "<<i<<" events done"<<endl;
			//}
			//if(nentries>30000){
				if(i%10000==1) cout<<"**** "<<i<<" events done"<<endl;
			//}
		}

		///// 512 Channels
		for (Int_t j=0; j<512; j++) {
			// MWPCs
			h_MWPCADC[j]->Fill(MwpcADC[j]);
		}

		///// 256 Channels
		for (Int_t j=0; j<256; j++){
			// TARGET
			h_ADC_High_TARGET[j]->Fill(adc_high_target[j]);
			h_ADC_Low_TARGET[j]->Fill(adc_low_target[j]);
			h_TDC_LE_TARGET[j]->Fill(tdc_le_target[j][0]);
			h_TDC_TE_TARGET[j]->Fill(tdc_te_target[j][0]);

			// V792
			h_V792_adc[j]->Fill(V792_adc[j]);

			// VT48
			if(VT48_tdc[j][0]>=0)  h_VT48_tdc[j]->Fill(VT48_tdc[j][0]);
		}

		///// 128 Channels
		for (Int_t j=0; j<128; j++){
			// SFT
			h_ADC_High_SFT[j]->Fill(ADC_High_sft[j]);
			h_ADC_Low_SFT[j]->Fill(ADC_Low_sft[j]);
			h_TDC_LE_SFT[j]->Fill(TDC_LE_sft[j][0]);
			h_TDC_TE_SFT[j]->Fill(TDC_TE_sft[j][0]);
		}

		///// 96 Channels
		for (Int_t j=0; j<96; j++){
			// TOF
			//h_ADC_tof[j]->Fill(ADC_tof[j]);
			//h_TDC_tof[j]->Fill(TDC_tof[j]);

			// PGC
			//h_ADC_pgc[j]->Fill(ADC_pgc[j]);

			// HR TDC
			h_HR_tdc[j]->Fill(HR_tdc[j]);

			// TKO ADC
			h_TKO_adc[j]->Fill(TKO_adc[j]);
		}
	
		///// 72 Channels
		for (Int_t j=0; j<72; j++){
			h_ADC_C4XR[j]->Fill(ADC_C4XR[j]);
			h_ADC_C4XL[j]->Fill(ADC_C4XL[j]);
		}
	
		///// 64 Channels
		for (Int_t j=0; j<64; j++){
			h_ADC_C3XR[j]->Fill(ADC_C3XR[j]);
			h_ADC_C3XL[j]->Fill(ADC_C3XL[j]);
		}
	
		///// 56 Channels
		for (Int_t j=0; j<56; j++){
			h_ADC_C2XR[j]->Fill(ADC_C2XR[j]);
			h_ADC_C2XL[j]->Fill(ADC_C2XL[j]);
		}

		///// 16 Channels
		for (Int_t j=0; j<16; j++){
			h_ADC_C2YR[j]->Fill(ADC_C2YR[j]);
			h_ADC_C2YL[j]->Fill(ADC_C2YL[j]);
			h_ADC_C3YR[j]->Fill(ADC_C3YR[j]);
			h_ADC_C3YL[j]->Fill(ADC_C3YL[j]);
			h_ADC_C4YR[j]->Fill(ADC_C4YR[j]);
			h_ADC_C4YL[j]->Fill(ADC_C4YL[j]);
		}
		

		///// 84 Channels
		// PGC
		//for (Int_t j=0; j<84; j++)	h_TDC_pgc[j]->Fill(TDC_pgc[j][0]);

		///// 56 Channels
		//for (Int_t j=0; j<56; j++) {
			// TOF2
		//	h_ADC_tof2[j]->Fill(ADC_tof2[j]);
		//	h_TDC_tof2[j]->Fill(TDC_tof2[j]);
		//}

		///// 24 Channels
		//for (Int_t j=0; j<24; j++) {
			// TOF1
		//	h_ADC_tof1[j]->Fill(ADC_tof1[j]);
		//	h_TDC_tof1[j]->Fill(TDC_tof1[j]);
	
			// AC
		//	h_ADC_ac[j]->Fill(ADC_ac[j]);
		//}

		///// 12 Channels
		for (Int_t j=0; j<12; j++) {
			// TOF1
			h_TDC_tof1U[j]->Fill(TDC_tof1U[j]);
			h_TDC_tof1D[j]->Fill(TDC_tof1D[j]);
			h_ADC_tof1U[j]->Fill(ADC_tof1U[j]);
			h_ADC_tof1D[j]->Fill(ADC_tof1D[j]);
	
			// TOF2
			h_ADC_tof2AO[j]->Fill(ADC_tof2AO[j]);
			h_ADC_tof2AI[j]->Fill(ADC_tof2AI[j]);
			h_ADC_tof2BO[j]->Fill(ADC_tof2BO[j]);
			h_ADC_tof2BI[j]->Fill(ADC_tof2BI[j]);
			h_TDC_tof2AO[j]->Fill(TDC_tof2AO[j]);
			h_TDC_tof2AI[j]->Fill(TDC_tof2AI[j]);
			h_TDC_tof2BO[j]->Fill(TDC_tof2BO[j]);
			h_TDC_tof2BI[j]->Fill(TDC_tof2BI[j]);

			// AC
			h_ADC_acu[j]->Fill(ADC_acu[j]);
			h_ADC_acd[j]->Fill(ADC_acd[j]);
			h_TDC_acu[j]->Fill(TDC_acu[j][0]);
			h_TDC_acd[j]->Fill(TDC_acd[j][0]);

			// GV
			h_ADC_gv[j]->Fill(ADC_gv[j]);
			h_TDC_gv[j]->Fill(TDC_gv[j][0]);

			// TTC
			h_ADC_ttc[j]->Fill(ADC_ttc[j]);
			h_TDC_ttc[j]->Fill(TDC_ttc[j][0]);

			// PGC
			for (int k=0; k<8; k++){
				h_ADC_pgc_sorted[j][k]->Fill(ADC_pgc[8*j+k]);
				h_TDC_pgc_sorted[j][k]->Fill(TDC_pgc[8*j+k][0]);

				h_ADC_pgc_Gap[j][k]->Fill(ADC_pgc_Gap[j][k]);
				h_TDC_pgc_Gap[j][k]->Fill(TDC_pgc_Gap[j][k][0]);
			}
		}

		for(int j=0; j<14; j++){
			h_TDC_Ck[j]->Fill(TDC_ck[j][0]);
			h_TDC_Cpi[j]->Fill(TDC_cpi[j][0]);

			for(int k=0; k<8; k++){
				if(TDC_ck[j][k]>=TDC_Ck_min && TDC_ck[j][k]<=TDC_Ck_max){
					TDC_ck_selected[j]=TDC_ck[j][k];
				}
			}
			if(TDC_ck_selected[j]>0) vec_Ck_Avg.push_back(TDC_ck_selected[j]);
			TDC_ck_selected[j]=0;
		}

		for(unsigned j=0; j<vec_Ck_Avg.size(); j++) TDC_ck_sum += vec_Ck_Avg[j];

		h_TDC_Ck_Avg->Fill(TDC_ck_sum/vec_Ck_Avg.size());

		vec_Ck_Avg.clear();
		TDC_ck_sum = 0;

	} // EndLoop over Events


	///////////////////////////////////////////////////////////////////////////////////////////////////
	// CREATE DIRECTORIES

	// SFT
	TDirectory *dir_Raw_SFT = foutput->mkdir("RAW_SFT");
	TDirectory *dir_Raw_SFT_ADC = dir_Raw_SFT->mkdir("ADC");
	TDirectory *sub_Raw_ADC_High_SFT = dir_Raw_SFT_ADC->mkdir("Raw_ADC_High_SFT");
	TDirectory *sub_Raw_ADC_Low_SFT = dir_Raw_SFT_ADC->mkdir("Raw_ADC_Low_SFT");

	TDirectory *dir_Raw_SFT_TDC = dir_Raw_SFT->mkdir("TDC");
	TDirectory *sub_Raw_TDC_LE_SFT = dir_Raw_SFT_TDC->mkdir("Raw_TDC_LE_SFT");
	TDirectory *sub_Raw_TDC_TE_SFT = dir_Raw_SFT_TDC->mkdir("Raw_TDC_TE_SFT");

	// TARGET
	TDirectory *dir_Raw_TARGET = foutput->mkdir("RAW_TARGET");
	TDirectory *dir_Raw_TARGET_ADC = dir_Raw_TARGET->mkdir("ADC");
	TDirectory *sub_Raw_ADC_High_TARGET = dir_Raw_TARGET_ADC->mkdir("Raw_ADC_High_TARGET");
	TDirectory *sub_Raw_ADC_Low_TARGET = dir_Raw_TARGET_ADC->mkdir("Raw_ADC_Low_TARGET");

	TDirectory *dir_Raw_TARGET_TDC = dir_Raw_TARGET->mkdir("TDC");
	TDirectory *sub_Raw_TDC_LE_TARGET = dir_Raw_TARGET_TDC->mkdir("Raw_TDC_LE_TARGET");
	TDirectory *sub_Raw_TDC_TE_TARGET = dir_Raw_TARGET_TDC->mkdir("Raw_TDC_TE_TARGET");

	// MWPCs
	TDirectory *dir_Raw_MWPC = foutput->mkdir("RAW_MWPC");
	TDirectory *sub_Raw_MWPCADC = dir_Raw_MWPC->mkdir("Raw_ADC_MWPC");
	
	TDirectory *sub_Raw_MWPC_C2 = dir_Raw_MWPC->mkdir("Raw_ADC_C2");
	TDirectory *sub_Raw_MWPC_C3 = dir_Raw_MWPC->mkdir("Raw_ADC_C3");
	TDirectory *sub_Raw_MWPC_C4 = dir_Raw_MWPC->mkdir("Raw_ADC_C4");

	TDirectory *sub_Raw_ADC_C2XR = sub_Raw_MWPC_C2->mkdir("Raw_ADC_C2X_R");
	TDirectory *sub_Raw_ADC_C2XL = sub_Raw_MWPC_C2->mkdir("Raw_ADC_C2X_L");
	TDirectory *sub_Raw_ADC_C2YR = sub_Raw_MWPC_C2->mkdir("Raw_ADC_C2Y_R");
	TDirectory *sub_Raw_ADC_C2YL = sub_Raw_MWPC_C2->mkdir("Raw_ADC_C2Y_L");

	TDirectory *sub_Raw_ADC_C3XR = sub_Raw_MWPC_C3->mkdir("Raw_ADC_C3X_R");
	TDirectory *sub_Raw_ADC_C3XL = sub_Raw_MWPC_C3->mkdir("Raw_ADC_C3X_L");
	TDirectory *sub_Raw_ADC_C3YR = sub_Raw_MWPC_C3->mkdir("Raw_ADC_C3Y_R");
	TDirectory *sub_Raw_ADC_C3YL = sub_Raw_MWPC_C3->mkdir("Raw_ADC_C3Y_L");

	TDirectory *sub_Raw_ADC_C4XR = sub_Raw_MWPC_C4->mkdir("Raw_ADC_C4X_R");
	TDirectory *sub_Raw_ADC_C4XL = sub_Raw_MWPC_C4->mkdir("Raw_ADC_C4X_L");
	TDirectory *sub_Raw_ADC_C4YR = sub_Raw_MWPC_C4->mkdir("Raw_ADC_C4Y_R");
	TDirectory *sub_Raw_ADC_C4YL = sub_Raw_MWPC_C4->mkdir("Raw_ADC_C4Y_L");
	
	// TOF1
	TDirectory *dir_Raw_TOF1 = foutput->mkdir("RAW_TOF1");
	TDirectory *dir_Raw_TOF1_ADC = dir_Raw_TOF1->mkdir("ADC");
	TDirectory *sub_Raw_ADC_tof1U = dir_Raw_TOF1_ADC->mkdir("Raw_ADC_tof1U");
	TDirectory *sub_Raw_ADC_tof1D = dir_Raw_TOF1_ADC->mkdir("Raw_ADC_tof1D");

	TDirectory *dir_Raw_TOF1_TDC = dir_Raw_TOF1->mkdir("TDC");
	TDirectory *sub_Raw_TDC_tof1U = dir_Raw_TOF1_TDC->mkdir("Raw_TDC_tof1U");
	TDirectory *sub_Raw_TDC_tof1D = dir_Raw_TOF1_TDC->mkdir("Raw_TDC_tof1D");
	//TDirectory *dir_Raw_TOF1_OLD = dir_Raw_TOF1->mkdir("OLD");

	// TOF2
	TDirectory *dir_Raw_TOF2 = foutput->mkdir("RAW_TOF2");
	TDirectory *dir_Raw_TOF2_ADC = dir_Raw_TOF2->mkdir("ADC");

	TDirectory *sub_Raw_ADC_tof2AO = dir_Raw_TOF2_ADC->mkdir("Raw_ADC_tof2AO");
	TDirectory *sub_Raw_ADC_tof2AI = dir_Raw_TOF2_ADC->mkdir("Raw_ADC_tof2AI");
	TDirectory *sub_Raw_ADC_tof2BO = dir_Raw_TOF2_ADC->mkdir("Raw_ADC_tof2BO");
	TDirectory *sub_Raw_ADC_tof2BI = dir_Raw_TOF2_ADC->mkdir("Raw_ADC_tof2BI");

	TDirectory *dir_Raw_TOF2_TDC = dir_Raw_TOF2->mkdir("TDC");
	TDirectory *sub_Raw_TDC_tof2AO = dir_Raw_TOF2_TDC->mkdir("Raw_TDC_tof2AO");
	TDirectory *sub_Raw_TDC_tof2AI = dir_Raw_TOF2_TDC->mkdir("Raw_TDC_tof2AI");
	TDirectory *sub_Raw_TDC_tof2BO = dir_Raw_TOF2_TDC->mkdir("Raw_TDC_tof2BO");
	TDirectory *sub_Raw_TDC_tof2BI = dir_Raw_TOF2_TDC->mkdir("Raw_TDC_tof2BI");
	//TDirectory *dir_Raw_TOF2_OLD = dir_Raw_TOF2->mkdir("OLD");

	// AC
	TDirectory *dir_Raw_AC = foutput->mkdir("RAW_AC");
	TDirectory *dir_Raw_AC_ADC = dir_Raw_AC->mkdir("ADC");
	TDirectory *sub_Raw_ADC_acu = dir_Raw_AC_ADC->mkdir("Raw_ADC_ACU");
	TDirectory *sub_Raw_ADC_acd = dir_Raw_AC_ADC->mkdir("Raw_ADC_ACD");

	TDirectory *dir_Raw_AC_TDC = dir_Raw_AC->mkdir("TDC");
	TDirectory *sub_Raw_TDC_acu = dir_Raw_AC_TDC->mkdir("Raw_TDC_ACU");
	TDirectory *sub_Raw_TDC_acd = dir_Raw_AC_TDC->mkdir("Raw_TDC_ACD");
	//TDirectory *dir_Raw_AC_OLD = dir_Raw_AC->mkdir("OLD");

	// GV
	TDirectory *dir_Raw_GV = foutput->mkdir("RAW_GV");
	TDirectory *dir_Raw_GV_ADC = dir_Raw_GV->mkdir("ADC");
	TDirectory *sub_Raw_ADC_gv = dir_Raw_GV_ADC->mkdir("Raw_ADC_GV");

	TDirectory *dir_Raw_GV_TDC = dir_Raw_GV->mkdir("TDC");
	TDirectory *sub_Raw_TDC_gv = dir_Raw_GV_TDC->mkdir("Raw_TDC_GV");

	// TTC
	TDirectory *dir_Raw_TTC = foutput->mkdir("RAW_TTC");
	TDirectory *dir_Raw_TTC_ADC = dir_Raw_TTC->mkdir("ADC");
	TDirectory *sub_Raw_ADC_ttc = dir_Raw_TTC_ADC->mkdir("Raw_ADC_TTC");

	TDirectory *dir_Raw_TTC_TDC = dir_Raw_TTC->mkdir("TDC");
	TDirectory *sub_Raw_TDC_ttc = dir_Raw_TTC_TDC->mkdir("Raw_TDC_TTC");

	// Ck
	TDirectory *dir_Raw_Ck = foutput->mkdir("RAW_Ck");
	TDirectory *dir_Raw_Ck_TDC = dir_Raw_Ck->mkdir("TDC");
	TDirectory *sub_Raw_TDC_Ck = dir_Raw_Ck_TDC->mkdir("Raw_TDC_Ck");

	// Cpi
	TDirectory *dir_Raw_Cpi = foutput->mkdir("RAW_Cpi");
	TDirectory *dir_Raw_Cpi_TDC = dir_Raw_Cpi->mkdir("TDC");
	TDirectory *sub_Raw_TDC_Cpi = dir_Raw_Cpi_TDC->mkdir("Raw_TDC_Cpi");

	// PGC
	TDirectory *dir_Raw_PGC_sorted = foutput->mkdir("RAW_PGC");
	TDirectory *dir_Raw_PGC_ADC = dir_Raw_PGC_sorted->mkdir("ADC");
	TDirectory *dir_Raw_PGC_ADC_Gap1 = dir_Raw_PGC_ADC->mkdir("Counter 1");
	TDirectory *dir_Raw_PGC_ADC_Gap2 = dir_Raw_PGC_ADC->mkdir("Counter 2");
	TDirectory *dir_Raw_PGC_ADC_Gap3 = dir_Raw_PGC_ADC->mkdir("Counter 3");
	TDirectory *dir_Raw_PGC_ADC_Gap4 = dir_Raw_PGC_ADC->mkdir("Counter 4");
	TDirectory *dir_Raw_PGC_ADC_Gap5 = dir_Raw_PGC_ADC->mkdir("Counter 5");
	TDirectory *dir_Raw_PGC_ADC_Gap6 = dir_Raw_PGC_ADC->mkdir("Counter 6");
	TDirectory *dir_Raw_PGC_ADC_Gap7 = dir_Raw_PGC_ADC->mkdir("Counter 7");
	TDirectory *dir_Raw_PGC_ADC_Gap8 = dir_Raw_PGC_ADC->mkdir("Counter 8");
	TDirectory *dir_Raw_PGC_ADC_Gap9 = dir_Raw_PGC_ADC->mkdir("Counter 9");
	TDirectory *dir_Raw_PGC_ADC_Gap10 = dir_Raw_PGC_ADC->mkdir("Counter 10");
	TDirectory *dir_Raw_PGC_ADC_Gap11 = dir_Raw_PGC_ADC->mkdir("Counter 11");
	TDirectory *dir_Raw_PGC_ADC_Gap12 = dir_Raw_PGC_ADC->mkdir("Counter 12");

	TDirectory *dir_Raw_PGC_TDC = dir_Raw_PGC_sorted->mkdir("TDC");
	TDirectory *dir_Raw_PGC_TDC_Gap1 = dir_Raw_PGC_TDC->mkdir("Counter 1");
	TDirectory *dir_Raw_PGC_TDC_Gap2 = dir_Raw_PGC_TDC->mkdir("Counter 2");
	TDirectory *dir_Raw_PGC_TDC_Gap3 = dir_Raw_PGC_TDC->mkdir("Counter 3");
	TDirectory *dir_Raw_PGC_TDC_Gap4 = dir_Raw_PGC_TDC->mkdir("Counter 4");
	TDirectory *dir_Raw_PGC_TDC_Gap5 = dir_Raw_PGC_TDC->mkdir("Counter 5");
	TDirectory *dir_Raw_PGC_TDC_Gap6 = dir_Raw_PGC_TDC->mkdir("Counter 6");
	TDirectory *dir_Raw_PGC_TDC_Gap7 = dir_Raw_PGC_TDC->mkdir("Counter 7");
	TDirectory *dir_Raw_PGC_TDC_Gap8 = dir_Raw_PGC_TDC->mkdir("Counter 8");
	TDirectory *dir_Raw_PGC_TDC_Gap9 = dir_Raw_PGC_TDC->mkdir("Counter 9");
	TDirectory *dir_Raw_PGC_TDC_Gap10 = dir_Raw_PGC_TDC->mkdir("Counter 10");
	TDirectory *dir_Raw_PGC_TDC_Gap11 = dir_Raw_PGC_TDC->mkdir("Counter 11");	
	TDirectory *dir_Raw_PGC_TDC_Gap12 = dir_Raw_PGC_TDC->mkdir("Counter 12");


	TDirectory *dir_Raw_PGC_Gap = foutput->mkdir("RAW_PGC_Gap");
	TDirectory *dir_Raw_PGC_ADC_Gap = dir_Raw_PGC_Gap->mkdir("ADC");
	TDirectory *dir_Raw_PGC_ADC_Gap1_Gap = dir_Raw_PGC_ADC_Gap->mkdir("Counter 1");
	TDirectory *dir_Raw_PGC_ADC_Gap2_Gap = dir_Raw_PGC_ADC_Gap->mkdir("Counter 2");
	TDirectory *dir_Raw_PGC_ADC_Gap3_Gap = dir_Raw_PGC_ADC_Gap->mkdir("Counter 3");
	TDirectory *dir_Raw_PGC_ADC_Gap4_Gap = dir_Raw_PGC_ADC_Gap->mkdir("Counter 4");
	TDirectory *dir_Raw_PGC_ADC_Gap5_Gap = dir_Raw_PGC_ADC_Gap->mkdir("Counter 5");
	TDirectory *dir_Raw_PGC_ADC_Gap6_Gap = dir_Raw_PGC_ADC_Gap->mkdir("Counter 6");
	TDirectory *dir_Raw_PGC_ADC_Gap7_Gap = dir_Raw_PGC_ADC_Gap->mkdir("Counter 7");
	TDirectory *dir_Raw_PGC_ADC_Gap8_Gap = dir_Raw_PGC_ADC_Gap->mkdir("Counter 8");
	TDirectory *dir_Raw_PGC_ADC_Gap9_Gap = dir_Raw_PGC_ADC_Gap->mkdir("Counter 9");
	TDirectory *dir_Raw_PGC_ADC_Gap10_Gap = dir_Raw_PGC_ADC_Gap->mkdir("Counter 10");
	TDirectory *dir_Raw_PGC_ADC_Gap11_Gap = dir_Raw_PGC_ADC_Gap->mkdir("Counter 11");
	TDirectory *dir_Raw_PGC_ADC_Gap12_Gap = dir_Raw_PGC_ADC_Gap->mkdir("Counter 12");

	TDirectory *dir_Raw_PGC_TDC_Gap = dir_Raw_PGC_Gap->mkdir("TDC");
	TDirectory *dir_Raw_PGC_TDC_Gap1_Gap = dir_Raw_PGC_TDC_Gap->mkdir("Counter 1");
	TDirectory *dir_Raw_PGC_TDC_Gap2_Gap = dir_Raw_PGC_TDC_Gap->mkdir("Counter 2");	
	TDirectory *dir_Raw_PGC_TDC_Gap3_Gap = dir_Raw_PGC_TDC_Gap->mkdir("Counter 3");
	TDirectory *dir_Raw_PGC_TDC_Gap4_Gap = dir_Raw_PGC_TDC_Gap->mkdir("Counter 4");
	TDirectory *dir_Raw_PGC_TDC_Gap5_Gap = dir_Raw_PGC_TDC_Gap->mkdir("Counter 5");
	TDirectory *dir_Raw_PGC_TDC_Gap6_Gap = dir_Raw_PGC_TDC_Gap->mkdir("Counter 6");
	TDirectory *dir_Raw_PGC_TDC_Gap7_Gap = dir_Raw_PGC_TDC_Gap->mkdir("Counter 7");
	TDirectory *dir_Raw_PGC_TDC_Gap8_Gap = dir_Raw_PGC_TDC_Gap->mkdir("Counter 8");
	TDirectory *dir_Raw_PGC_TDC_Gap9_Gap = dir_Raw_PGC_TDC_Gap->mkdir("Counter 9");
	TDirectory *dir_Raw_PGC_TDC_Gap10_Gap = dir_Raw_PGC_TDC_Gap->mkdir("Counter 10");
	TDirectory *dir_Raw_PGC_TDC_Gap11_Gap = dir_Raw_PGC_TDC_Gap->mkdir("Counter 11");
	TDirectory *dir_Raw_PGC_TDC_Gap12_Gap = dir_Raw_PGC_TDC_Gap->mkdir("Counter 12");


	//TDirectory *dir_Raw_OLD = foutput->mkdir("OLD");
	//TDirectory *sub_Raw_ADC_tof = dir_Raw_OLD->mkdir("Raw_ADC_TOF");
	//TDirectory *sub_Raw_TDC_tof = dir_Raw_OLD->mkdir("Raw_TDC_TOF");
	//TDirectory *sub_Raw_ADC_tof1 = dir_Raw_TOF1_OLD->mkdir("Raw_ADC_tof1");
	//TDirectory *sub_Raw_TDC_tof1 = dir_Raw_TOF1_OLD->mkdir("Raw_TDC_tof1");

	//TDirectory *sub_Raw_ADC_tof2 = dir_Raw_TOF2_OLD->mkdir("Raw_ADC_tof2");
	//TDirectory *sub_Raw_TDC_tof2 = dir_Raw_TOF2_OLD->mkdir("Raw_TDC_tof2");
	//TDirectory *sub_Raw_ADC_ac = dir_Raw_AC_OLD->mkdir("Raw_ADC_AC");


	//TDirectory *dir_Raw_PGC = dir_Raw_OLD->mkdir("RAW (PGC)");
	//TDirectory *sub_Raw_ADC_pgc = dir_Raw_PGC->mkdir("Raw_ADC_PGC");
	//TDirectory *sub_Raw_TDC_pgc = dir_Raw_PGC->mkdir("Raw_TDC_PGC");


	TDirectory *dir_Raw_HRTDC = foutput->mkdir("RAW_HR_TDC");
	TDirectory *dir_Raw_TKO_ADC = foutput->mkdir("RAW_TKO_ADC");
	TDirectory *dir_Raw_V792_ADC = foutput->mkdir("RAW_V792_ADC");
	TDirectory *dir_Raw_VT48_TDC = foutput->mkdir("RAW_VT48_TDC");


	//TDirectory *dir_Raw_TOF1U = dir_Raw_OLD->mkdir("RAW_TOF1_UP");
	//TDirectory *dir_Raw_TOF1D = dir_Raw_OLD->mkdir("RAW_TOF1_DOWN");


	//////////////////////////////////////////////////////////////////////////////////////////////////
	// WRITING HISTOGRAMS

	foutput->cd();

	///// 512 Channels
	for (Int_t ii=0; ii<512; ii++){
		// MWPCs
		sub_Raw_MWPCADC->cd();
		h_MWPCADC[ii]->SetDirectory(sub_Raw_MWPCADC);
		h_MWPCADC[ii]->Write();
	}

	//// 256 Channels

	for (Int_t ii=0; ii<256; ii++){
		// TARGET
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

		// V792
		dir_Raw_V792_ADC->cd();
		h_V792_adc[ii]->SetDirectory(dir_Raw_V792_ADC);
		h_V792_adc[ii]->Write();

		// VT48
		dir_Raw_VT48_TDC->cd();
		h_VT48_tdc[ii]->SetDirectory(dir_Raw_VT48_TDC);
		h_VT48_tdc[ii]->Write();
	}

	///// 128 Channels
	for (Int_t ii=0; ii<128; ii++){
		// SFT
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

	///// 96 Channels
	for (Int_t ii=0; ii<96; ii++){
		// TOF
		//sub_Raw_ADC_tof->cd();
		//h_ADC_tof[ii]->SetDirectory(sub_Raw_ADC_tof);
		//h_ADC_tof[ii]->Write();

		//sub_Raw_TDC_tof->cd();
		//h_TDC_tof[ii]->SetDirectory(sub_Raw_TDC_tof);
		//h_TDC_tof[ii]->Write();

		// HR TDC
		dir_Raw_HRTDC->cd();
		h_HR_tdc[ii]->SetDirectory(dir_Raw_HRTDC);
		h_HR_tdc[ii]->Write();

		// TKO ADC
		dir_Raw_TKO_ADC->cd();
		h_TKO_adc[ii]->SetDirectory(dir_Raw_TKO_ADC);
		h_TKO_adc[ii]->Write();

		// PGC
		//sub_Raw_ADC_pgc->cd();
		//h_ADC_pgc[ii]->SetDirectory(sub_Raw_ADC_pgc);
		//h_ADC_pgc[ii]->Write();
	}
	
	///// 72 Channels
	for (Int_t ii=0; ii<72; ii++){
		sub_Raw_ADC_C4XR->cd();
		h_ADC_C4XR[ii]->SetDirectory(sub_Raw_ADC_C4XR);
		h_ADC_C4XR[ii]->Write();

		sub_Raw_ADC_C4XL->cd();
		h_ADC_C4XL[ii]->SetDirectory(sub_Raw_ADC_C4XL);
		h_ADC_C4XL[ii]->Write();
	}

	///// 64 Channels
	for (Int_t ii=0; ii<64; ii++){
		sub_Raw_ADC_C3XR->cd();
		h_ADC_C3XR[ii]->SetDirectory(sub_Raw_ADC_C3XR);
		h_ADC_C3XR[ii]->Write();

		sub_Raw_ADC_C3XL->cd();
		h_ADC_C3XL[ii]->SetDirectory(sub_Raw_ADC_C3XL);
		h_ADC_C3XL[ii]->Write();
	}

	///// 56 Channels
	for (Int_t ii=0; ii<56; ii++){
		sub_Raw_ADC_C2XR->cd();
		h_ADC_C2XR[ii]->SetDirectory(sub_Raw_ADC_C2XR);
		h_ADC_C2XR[ii]->Write();

		sub_Raw_ADC_C2XL->cd();
		h_ADC_C2XL[ii]->SetDirectory(sub_Raw_ADC_C2XL);
		h_ADC_C2XL[ii]->Write();
	}

	///// 16 Channels
	for (Int_t ii=0; ii<16; ii++){
		sub_Raw_ADC_C2YR->cd();
		h_ADC_C2YR[ii]->SetDirectory(sub_Raw_ADC_C2YR);
		h_ADC_C2YR[ii]->Write();

		sub_Raw_ADC_C2YL->cd();
		h_ADC_C2YL[ii]->SetDirectory(sub_Raw_ADC_C2YL);
		h_ADC_C2YL[ii]->Write();
	
		sub_Raw_ADC_C3YR->cd();
		h_ADC_C3YR[ii]->SetDirectory(sub_Raw_ADC_C3YR);
		h_ADC_C3YR[ii]->Write();

		sub_Raw_ADC_C3YL->cd();
		h_ADC_C3YL[ii]->SetDirectory(sub_Raw_ADC_C3YL);
		h_ADC_C3YL[ii]->Write();

		sub_Raw_ADC_C4YR->cd();
		h_ADC_C4YR[ii]->SetDirectory(sub_Raw_ADC_C4YR);
		h_ADC_C4YR[ii]->Write();

		sub_Raw_ADC_C4YL->cd();
		h_ADC_C4YL[ii]->SetDirectory(sub_Raw_ADC_C4YL);
		h_ADC_C4YL[ii]->Write();
	}	

	///// 84 Channels
	//for (Int_t ii=0; ii<84; ii++){
		// PGC
	//	sub_Raw_TDC_pgc->cd();
	//	h_TDC_pgc[ii]->SetDirectory(sub_Raw_TDC_pgc);
	//	h_TDC_pgc[ii]->Write();
	//}

	///// 56 Channels
	//for (Int_t ii=0; ii<56; ii++){
		// TOF2
	//	sub_Raw_ADC_tof2->cd();	
	//	h_ADC_tof2[ii]->SetDirectory(sub_Raw_ADC_tof2);
	//	h_ADC_tof2[ii]->Write();

	//	sub_Raw_TDC_tof2->cd();
	//	h_TDC_tof2[ii]->SetDirectory(sub_Raw_TDC_tof2);
	//	h_TDC_tof2[ii]->Write();
	//}

	///// 24 Channels
	//for (Int_t ii=0; ii<24; ii++){
		// TOF1
	//	sub_Raw_ADC_tof1->cd();
	//	h_ADC_tof1[ii]->SetDirectory(sub_Raw_ADC_tof1);
	//	h_ADC_tof1[ii]->Write();

	//	sub_Raw_TDC_tof1->cd();
	//	h_TDC_tof1[ii]->SetDirectory(sub_Raw_TDC_tof1);
	//	h_TDC_tof1[ii]->Write();

		// AC
	//	sub_Raw_ADC_ac->cd();
	//	h_ADC_ac[ii]->SetDirectory(sub_Raw_ADC_ac);
	//	h_ADC_ac[ii]->Write();
	//}

	///// 12 Channels
	for (Int_t ii=0; ii<12; ii++){
		// TOF1
		sub_Raw_ADC_tof1U->cd();
		h_ADC_tof1U[ii]->SetDirectory(sub_Raw_ADC_tof1U);
		h_ADC_tof1U[ii]->Write();

		sub_Raw_ADC_tof1D->cd();
		h_ADC_tof1D[ii]->SetDirectory(sub_Raw_ADC_tof1D);
		h_ADC_tof1D[ii]->Write();

		sub_Raw_TDC_tof1U->cd();
		h_TDC_tof1U[ii]->SetDirectory(sub_Raw_TDC_tof1U);
		h_TDC_tof1U[ii]->Write();

		sub_Raw_TDC_tof1D->cd();
		h_TDC_tof1D[ii]->SetDirectory(sub_Raw_TDC_tof1D);
		h_TDC_tof1D[ii]->Write();

		// TOF2
		sub_Raw_ADC_tof2AO->cd();
		h_ADC_tof2AO[ii]->SetDirectory(sub_Raw_ADC_tof2AO);
		h_ADC_tof2AO[ii]->Write();

		sub_Raw_ADC_tof2AI->cd();
		h_ADC_tof2AI[ii]->SetDirectory(sub_Raw_ADC_tof2AI);
		h_ADC_tof2AI[ii]->Write();

		sub_Raw_ADC_tof2BO->cd();
		h_ADC_tof2BO[ii]->SetDirectory(sub_Raw_ADC_tof2BO);
		h_ADC_tof2BO[ii]->Write();

		sub_Raw_ADC_tof2BI->cd();
		h_ADC_tof2BI[ii]->SetDirectory(sub_Raw_ADC_tof2BI);
		h_ADC_tof2BI[ii]->Write();

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

		// AC
		sub_Raw_ADC_acu->cd();
		h_ADC_acu[ii]->SetDirectory(sub_Raw_ADC_acu);
		h_ADC_acu[ii]->Write();

		sub_Raw_ADC_acd->cd();
		h_ADC_acd[ii]->SetDirectory(sub_Raw_ADC_acd);
		h_ADC_acd[ii]->Write();

		sub_Raw_TDC_acu->cd();
		h_TDC_acu[ii]->SetDirectory(sub_Raw_TDC_acu);
		h_TDC_acu[ii]->Write();

		sub_Raw_TDC_acd->cd();
		h_TDC_acd[ii]->SetDirectory(sub_Raw_TDC_acd);
		h_TDC_acd[ii]->Write();

		// GV
		sub_Raw_ADC_gv->cd();
		h_ADC_gv[ii]->SetDirectory(sub_Raw_ADC_gv);
		h_ADC_gv[ii]->Write();

		sub_Raw_TDC_gv->cd();
		h_TDC_gv[ii]->SetDirectory(sub_Raw_TDC_gv);
		h_TDC_gv[ii]->Write();

		// TTC
		sub_Raw_ADC_ttc->cd();
		h_ADC_ttc[ii]->SetDirectory(sub_Raw_ADC_ttc);
		h_ADC_ttc[ii]->Write();

		sub_Raw_TDC_ttc->cd();
		h_TDC_ttc[ii]->SetDirectory(sub_Raw_TDC_ttc);
		h_TDC_ttc[ii]->Write();
	}
	
	sub_Raw_TDC_Ck->cd();
	h_TDC_Ck_Avg->Write();
	for(int ii=0; ii<14; ii++){
		sub_Raw_TDC_Ck->cd();
		h_TDC_Ck[ii]->SetDirectory(sub_Raw_TDC_Ck);
		h_TDC_Ck[ii]->Write();

		sub_Raw_TDC_Cpi->cd();
		h_TDC_Cpi[ii]->SetDirectory(sub_Raw_TDC_Cpi);
		h_TDC_Cpi[ii]->Write();
	}

	for (int ii=0; ii<8; ii++){
		// PGC ADC
		dir_Raw_PGC_ADC_Gap1->cd();
		h_ADC_pgc_sorted[0][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap1);
		h_ADC_pgc_sorted[0][ii]->Write();

		dir_Raw_PGC_ADC_Gap2->cd();
		h_ADC_pgc_sorted[1][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap2);
		h_ADC_pgc_sorted[1][ii]->Write();

		dir_Raw_PGC_ADC_Gap3->cd();
		h_ADC_pgc_sorted[2][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap3);
		h_ADC_pgc_sorted[2][ii]->Write();

		dir_Raw_PGC_ADC_Gap4->cd();
		h_ADC_pgc_sorted[3][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap4);
		h_ADC_pgc_sorted[3][ii]->Write();

		dir_Raw_PGC_ADC_Gap5->cd();
		h_ADC_pgc_sorted[4][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap5);
		h_ADC_pgc_sorted[4][ii]->Write();

		dir_Raw_PGC_ADC_Gap6->cd();
		h_ADC_pgc_sorted[5][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap6);
		h_ADC_pgc_sorted[5][ii]->Write();

		dir_Raw_PGC_ADC_Gap7->cd();
		h_ADC_pgc_sorted[6][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap7);
		h_ADC_pgc_sorted[6][ii]->Write();

		dir_Raw_PGC_ADC_Gap8->cd();
		h_ADC_pgc_sorted[7][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap8);
		h_ADC_pgc_sorted[7][ii]->Write();

		dir_Raw_PGC_ADC_Gap9->cd();
		h_ADC_pgc_sorted[8][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap9);
		h_ADC_pgc_sorted[8][ii]->Write();

		dir_Raw_PGC_ADC_Gap10->cd();
		h_ADC_pgc_sorted[9][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap10);
		h_ADC_pgc_sorted[9][ii]->Write();

		dir_Raw_PGC_ADC_Gap11->cd();
		h_ADC_pgc_sorted[10][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap11);
		h_ADC_pgc_sorted[10][ii]->Write();

		dir_Raw_PGC_ADC_Gap12->cd();
		h_ADC_pgc_sorted[11][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap12);
		h_ADC_pgc_sorted[11][ii]->Write();

		//// PGC TDC
		dir_Raw_PGC_TDC_Gap1->cd();
		h_TDC_pgc_sorted[0][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap1);
		h_TDC_pgc_sorted[0][ii]->Write();

		dir_Raw_PGC_TDC_Gap2->cd();
		h_TDC_pgc_sorted[1][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap2);
		h_TDC_pgc_sorted[1][ii]->Write();

		dir_Raw_PGC_TDC_Gap3->cd();
		h_TDC_pgc_sorted[2][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap3);
		h_TDC_pgc_sorted[2][ii]->Write();

		dir_Raw_PGC_TDC_Gap4->cd();
		h_TDC_pgc_sorted[3][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap4);
		h_TDC_pgc_sorted[3][ii]->Write();

		dir_Raw_PGC_TDC_Gap5->cd();
		h_TDC_pgc_sorted[4][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap5);
		h_TDC_pgc_sorted[4][ii]->Write();

		dir_Raw_PGC_TDC_Gap6->cd();
		h_TDC_pgc_sorted[5][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap6);
		h_TDC_pgc_sorted[5][ii]->Write();

		dir_Raw_PGC_TDC_Gap7->cd();
		h_TDC_pgc_sorted[6][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap7);
		h_TDC_pgc_sorted[6][ii]->Write();

		dir_Raw_PGC_TDC_Gap8->cd();
		h_TDC_pgc_sorted[7][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap8);
		h_TDC_pgc_sorted[7][ii]->Write();

		dir_Raw_PGC_TDC_Gap9->cd();
		h_TDC_pgc_sorted[8][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap9);
		h_TDC_pgc_sorted[8][ii]->Write();

		dir_Raw_PGC_TDC_Gap10->cd();
		h_TDC_pgc_sorted[9][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap10);
		h_TDC_pgc_sorted[9][ii]->Write();

		dir_Raw_PGC_TDC_Gap11->cd();
		h_TDC_pgc_sorted[10][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap11);
		h_TDC_pgc_sorted[10][ii]->Write();

		dir_Raw_PGC_TDC_Gap12->cd();
		h_TDC_pgc_sorted[11][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap12);
		h_TDC_pgc_sorted[11][ii]->Write();
	}


	for (int ii=0; ii<8; ii++){
		// PGC ADC
		dir_Raw_PGC_ADC_Gap1_Gap->cd();
		h_ADC_pgc_Gap[0][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap1_Gap);
		h_ADC_pgc_Gap[0][ii]->Write();

		dir_Raw_PGC_ADC_Gap2_Gap->cd();
		h_ADC_pgc_Gap[1][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap2_Gap);
		h_ADC_pgc_Gap[1][ii]->Write();

		dir_Raw_PGC_ADC_Gap3_Gap->cd();
		h_ADC_pgc_Gap[2][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap3_Gap);
		h_ADC_pgc_Gap[2][ii]->Write();

		dir_Raw_PGC_ADC_Gap4_Gap->cd();
		h_ADC_pgc_Gap[3][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap4_Gap);
		h_ADC_pgc_Gap[3][ii]->Write();

		dir_Raw_PGC_ADC_Gap5_Gap->cd();
		h_ADC_pgc_Gap[4][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap5_Gap);
		h_ADC_pgc_Gap[4][ii]->Write();

		dir_Raw_PGC_ADC_Gap6_Gap->cd();
		h_ADC_pgc_Gap[5][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap6_Gap);
		h_ADC_pgc_Gap[5][ii]->Write();

		dir_Raw_PGC_ADC_Gap7_Gap->cd();
		h_ADC_pgc_Gap[6][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap7_Gap);
		h_ADC_pgc_Gap[6][ii]->Write();

		dir_Raw_PGC_ADC_Gap8_Gap->cd();
		h_ADC_pgc_Gap[7][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap8_Gap);
		h_ADC_pgc_Gap[7][ii]->Write();

		dir_Raw_PGC_ADC_Gap9_Gap->cd();
		h_ADC_pgc_Gap[8][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap9_Gap);
		h_ADC_pgc_Gap[8][ii]->Write();

		dir_Raw_PGC_ADC_Gap10_Gap->cd();
		h_ADC_pgc_Gap[9][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap10_Gap);
		h_ADC_pgc_Gap[9][ii]->Write();

		dir_Raw_PGC_ADC_Gap11_Gap->cd();
		h_ADC_pgc_Gap[10][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap11_Gap);
		h_ADC_pgc_Gap[10][ii]->Write();

		dir_Raw_PGC_ADC_Gap12_Gap->cd();
		h_ADC_pgc_Gap[11][ii]->SetDirectory(dir_Raw_PGC_ADC_Gap12_Gap);
		h_ADC_pgc_Gap[11][ii]->Write();

		// PGC TDC
		dir_Raw_PGC_TDC_Gap1_Gap->cd();
		h_TDC_pgc_Gap[0][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap1_Gap);
		h_TDC_pgc_Gap[0][ii]->Write();

		dir_Raw_PGC_TDC_Gap2_Gap->cd();
		h_TDC_pgc_Gap[1][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap2_Gap);
		h_TDC_pgc_Gap[1][ii]->Write();

		dir_Raw_PGC_TDC_Gap3_Gap->cd();
		h_TDC_pgc_Gap[2][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap3_Gap);
		h_TDC_pgc_Gap[2][ii]->Write();

		dir_Raw_PGC_TDC_Gap4_Gap->cd();
		h_TDC_pgc_Gap[3][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap4_Gap);
		h_TDC_pgc_Gap[3][ii]->Write();

		dir_Raw_PGC_TDC_Gap5_Gap->cd();
		h_TDC_pgc_Gap[4][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap5_Gap);
		h_TDC_pgc_Gap[4][ii]->Write();

		dir_Raw_PGC_TDC_Gap6_Gap->cd();
		h_TDC_pgc_Gap[5][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap6_Gap);
		h_TDC_pgc_Gap[5][ii]->Write();

		dir_Raw_PGC_TDC_Gap7_Gap->cd();
		h_TDC_pgc_Gap[6][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap7_Gap);
		h_TDC_pgc_Gap[6][ii]->Write();

		dir_Raw_PGC_TDC_Gap8_Gap->cd();
		h_TDC_pgc_Gap[7][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap8_Gap);
		h_TDC_pgc_Gap[7][ii]->Write();

		dir_Raw_PGC_TDC_Gap9_Gap->cd();
		h_TDC_pgc_Gap[8][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap9_Gap);
		h_TDC_pgc_Gap[8][ii]->Write();

		dir_Raw_PGC_TDC_Gap10_Gap->cd();
		h_TDC_pgc_Gap[9][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap10_Gap);
		h_TDC_pgc_Gap[9][ii]->Write();

		dir_Raw_PGC_TDC_Gap11_Gap->cd();
		h_TDC_pgc_Gap[10][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap11_Gap);
		h_TDC_pgc_Gap[10][ii]->Write();

		dir_Raw_PGC_TDC_Gap12_Gap->cd();
		h_TDC_pgc_Gap[11][ii]->SetDirectory(dir_Raw_PGC_TDC_Gap12_Gap);
		h_TDC_pgc_Gap[11][ii]->Write();
	}

	cout << "" << endl;
	cout << "" << endl;
	cout << "Histograms have been saved!" << endl;
	cout << "" << endl;
	cout << "" << endl;


	foutput->Close();


	return;
}

