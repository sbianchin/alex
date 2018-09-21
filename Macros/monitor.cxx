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
#include "Thresholds_monitor.h"
#include "CommonParameters.h"
#endif


void monitor(Int_t Run_Number=5, Int_t flag_corr=0, Int_t time_window=0, Int_t flag=0){ 

	// flag_corr = 1    ---->   corrects for TOF1  3d <-> 4d inversion 

	gStyle->Clear();
	TH1::AddDirectory(kFALSE);
	gStyle->SetOptStat(1111111);


 	////////////////////////  Variable Definition  /////////////////////////////////////////////////////////////

	Int_t adc_high_target[256];   		Int_t ADC_High_TARGET[256];    	//Double_t ADC_High_TARGET_corr[256];   	
	Int_t adc_low_target[256]; 			Int_t ADC_Low_TARGET[256]; 	//Double_t ADC_Low_TARGET_corr[256];
	Int_t tdc_le_target[256][16];      	Int_t TDC_LE_TARGET[256];    	//Double_t TDC_LE_TARGET_corr[256];   	
	Int_t tdc_te_target[256][16]; 		Int_t TDC_TE_TARGET[256]; 	//Double_t TDC_TE_TARGET_corr[256];

	Int_t adc_high_sft[128];   			Int_t ADC_High_SFT[128];    	//Double_t ADC_High_SFT_corr[128];   	
	Int_t adc_low_sft[128]; 			Int_t ADC_Low_SFT[128]; 	//Double_t ADC_Low_SFT_corr[128];
	Int_t tdc_le_sft[128][16];      	Int_t TDC_LE_SFT[128];    	//Double_t TDC_LE_SFT_corr[128];   	
	Int_t tdc_te_sft[128][16]; 			Int_t TDC_TE_SFT[128]; 		//Double_t TDC_TE_SFT_corr[128];	

	Int_t ADC_tof1[24];		Int_t ADC_TOF1[24];	
	Int_t ADC_tof2[48];		Int_t ADC_TOF2[48];

	Int_t MwpcADC[512];	Int_t MWPCADC[512];
	//Int_t fiber[128]={-1};

	Int_t TARGET_hit_counter[256] = {0};

	Int_t TDC_thr_low = 100;
	Int_t TDC_thr_high = 1020;
	Int_t ADC_cut_TARGET2 = 2000;

	Int_t kaon_hit_counter[256] = {0};
	Int_t kaon_hit_counter2[256] = {0};
	Int_t path_margin=10;

	Int_t ADC_tof1U[12];  Int_t ADC_TOF1U[12];
	Int_t ADC_tof1D[12];  Int_t ADC_TOF1D[12];

	Int_t TDC_tof1U[12];  Int_t TDC_TOF1U[12];
	Int_t TDC_tof1D[12];  Int_t TDC_TOF1D[12];


	Int_t ADC_tof2AO[12];   Int_t ADC_TOF2AO[12];
	Int_t ADC_tof2BO[12];   Int_t ADC_TOF2BO[12];
	Int_t ADC_tof2AI[12];   Int_t ADC_TOF2AI[12];
	Int_t ADC_tof2BI[12];   Int_t ADC_TOF2BI[12];

	Int_t TDC_tof2AO[12];   Int_t TDC_TOF2AO[12];
	Int_t TDC_tof2BO[12];   Int_t TDC_TOF2BO[12];
	Int_t TDC_tof2AI[12];   Int_t TDC_TOF2AI[12];
	Int_t TDC_tof2BI[12];   Int_t TDC_TOF2BI[12];

	Int_t TOF_Hit_Counter[12] = {0};
	Int_t good_gap_counter = 0;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////


	char path_input[200];                   char file_mapping[200];
	sprintf(path_input,"%s",path_merged);	sprintf(file_mapping,"../Mapping");
	//sprintf(path_input,"/triumfcs/trshare/trek/E36/Data/April_2015/root");

	char Name_finput[200];
	sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);

	////

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

	fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);		//fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
	fChain->SetBranchAddress("ADC_Low_TARGET",adc_low_target);		//fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
	fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);		//fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
	//fChain->SetBranchAddress("TDC_TE_TARGET",tdc_te_target);		fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

	//fChain->SetBranchAddress("ADC_TOF1",ADC_tof1);
	//fChain->SetBranchAddress("ADC_TOF2",ADC_tof2);

	fChain->SetBranchAddress("TDC_TOF1U",TDC_tof1U);
	fChain->SetBranchAddress("TDC_TOF1D",TDC_tof1D);
	fChain->SetBranchAddress("TDC_TOF2AO",TDC_tof2AO);
	fChain->SetBranchAddress("TDC_TOF2AI",TDC_tof2AI);
	fChain->SetBranchAddress("TDC_TOF2BO",TDC_tof2BO);
	fChain->SetBranchAddress("TDC_TOF2BI",TDC_tof2BI); 

	fChain->SetBranchAddress("MWPCADC",MwpcADC);


	fChain->SetBranchAddress("ADC_TOF1U", ADC_tof1U);
	fChain->SetBranchAddress("ADC_TOF1D", ADC_tof1D);

	fChain->SetBranchAddress("TDC_TOF1U", TDC_tof1U);
	fChain->SetBranchAddress("TDC_TOF1D", TDC_tof1D);


	fChain->SetBranchAddress("ADC_TOF2AO", ADC_tof2AO);
	fChain->SetBranchAddress("ADC_TOF2AI", ADC_tof2AI);
	fChain->SetBranchAddress("ADC_TOF2BO", ADC_tof2BO);
	fChain->SetBranchAddress("ADC_TOF2BI", ADC_tof2BI);

	fChain->SetBranchAddress("TDC_TOF2AO", TDC_tof2AO);
	fChain->SetBranchAddress("TDC_TOF2AI", TDC_tof2AI);
	fChain->SetBranchAddress("TDC_TOF2BO", TDC_tof2BO);
	fChain->SetBranchAddress("TDC_TOF2BI", TDC_tof2BI);

	//cout << "TDC TOF1 Thresholds" << endl;
	//cout << TDC_TOF1_min << endl;
	//cout << TDC_TOF1_max << endl;

	//cout << "TDC TOF2 Thresholds" << endl;
	//cout << TDC_TOF2_min << endl;
	//cout << TDC_TOF2_max << endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Histo_Display (Raw Histograms)

	char path_input_Histo[100];
	sprintf(path_input_Histo,"%s",path_histo);

	char Name_finput_Histo[100];
	sprintf(Name_finput_Histo,"%s/Hist_Run%dM.root",path_input_Histo, Run_Number);
	cout << Name_finput_Histo << endl;
	cout << "     " << endl;

	TFile *finput = new TFile(Name_finput_Histo);


	//TH1F *h_ADC_High_Raw[6][64];
	//TH1F *h_ADC_Low_Raw[64];
	//TH1F *h_TDC_LE_Raw[64];
	//TH1F *h_TDC_TE_Raw[64];

	TH1F *h_ADC_tof1_Raw[24];
	TH1F *h_ADC_tof2_Raw[48];

	TH1F *h_TDC_tof1U_Raw[12];   	char Name_TDC_tof1U_Raw[100];
	TH1F *h_TDC_tof1D_Raw[12];   	char Name_TDC_tof1D_Raw[100];
	TH1F *h_TDC_tof2AO_Raw[12];  	char Name_TDC_tof2AO_Raw[100];
	TH1F *h_TDC_tof2AI_Raw[12];  	char Name_TDC_tof2AI_Raw[100];
	TH1F *h_TDC_tof2BO_Raw[12];  	char Name_TDC_tof2BO_Raw[100];
	TH1F *h_TDC_tof2BI_Raw[12];  	char Name_TDC_tof2BI_Raw[100];

	TH1F *h_ADC_tof1U_Raw[12];		char Name_ADC_tof1U_Raw[100];
	TH1F *h_ADC_tof1D_Raw[12];		char Name_ADC_tof1D_Raw[100];

	TH1F *h_ADC_tof2AO_Raw[12];		char Name_ADC_tof2AO_Raw[100];
	TH1F *h_ADC_tof2AI_Raw[12];		char Name_ADC_tof2AI_Raw[100];
	TH1F *h_ADC_tof2BO_Raw[12];		char Name_ADC_tof2BO_Raw[100];
	TH1F *h_ADC_tof2BI_Raw[12];		char Name_ADC_tof2BI_Raw[100];


	//TH1F *h_ADC_MWPC_Raw[512];


	//char Name_ADC_High_Raw[100];
	//char Name_ADC_Low_Raw[100];
	//char Name_TDC_LE_Raw[100];
	//char Name_TDC_TE_Raw[100];
	char Name_ADC_tof1_Raw[100];	
	char Name_ADC_tof2_Raw[100];
	//char Name_ADC_MWPC_Raw[100];

	
	/*
	for(int j=0;j<24;j++){
		sprintf(Name_ADC_tof1_Raw,"ADC_TOF1 (Ch. %d)",j);  // OLD

		h_ADC_tof1_Raw[j]=(TH1F*)finput->FindObjectAny(Name_ADC_tof1_Raw);	// OLD
	}
	*/

	
	for(int j=0; j<12; j++){
		sprintf(Name_ADC_tof1U_Raw,"ADC_TOF1 (%d UP)",j+1);         // NEW
		h_ADC_tof1_Raw[j]=(TH1F*)finput->FindObjectAny(Name_ADC_tof1U_Raw);	// NEW
	}

	for(int j=0; j<12; j++){
		sprintf(Name_ADC_tof1D_Raw,"ADC_TOF1 (%d DOWN)",j+1);         // NEW
		h_ADC_tof1_Raw[j+12]=(TH1F*)finput->FindObjectAny(Name_ADC_tof1D_Raw);	// NEW
	}
	





	
	for(int j=0; j<12; j++){
		//sprintf(Name_TDC_tof1U_Raw,"TDC_TOF1 (Ch. %d)",j);         // OLD
		//sprintf(Name_TDC_tof1D_Raw,"TDC_TOF1 (Ch. %d)",j+12);		 // OLD
		sprintf(Name_TDC_tof1U_Raw,"TDC_TOF1 (%d UP)",j+1);		// NEW
		sprintf(Name_TDC_tof1D_Raw,"TDC_TOF1 (%d DOWN)",j+1);		// NEW

		h_TDC_tof1U_Raw[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof1U_Raw);
		h_TDC_tof1D_Raw[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof1D_Raw);
	}




	//for(int j=0;j<56;j++){
	//	sprintf(Name_ADC_tof2_Raw,"ADC_TOF2 (Ch. %d)",j);	

	//	h_ADC_tof2_Raw[j]=(TH1F*)finput->FindObjectAny(Name_ADC_tof2_Raw);
	//}


	
	for(int j=0; j<12; j++){
		sprintf(Name_ADC_tof2AO_Raw,"ADC_TOF2 (%d AO)",j+1);         // NEW
		h_ADC_tof2_Raw[j]=(TH1F*)finput->FindObjectAny(Name_ADC_tof2AO_Raw);	// NEW
	}

	for(int j=0; j<12; j++){
		sprintf(Name_ADC_tof2AI_Raw,"ADC_TOF2 (%d AI)",j+1);         // NEW
		h_ADC_tof2_Raw[j+12]=(TH1F*)finput->FindObjectAny(Name_ADC_tof2AI_Raw);	// NEW
	}

	for(int j=0; j<12; j++){
		sprintf(Name_ADC_tof2BO_Raw,"ADC_TOF2 (%d BO)",j+1);         // NEW
		h_ADC_tof2_Raw[j+24]=(TH1F*)finput->FindObjectAny(Name_ADC_tof2BO_Raw);	// NEW
	}

	for(int j=0; j<12; j++){
		sprintf(Name_ADC_tof2BI_Raw,"ADC_TOF2 (%d BI)",j+1);         // NEW
		h_ADC_tof2_Raw[j+36]=(TH1F*)finput->FindObjectAny(Name_ADC_tof2BI_Raw);	// NEW
	}
	


	
	for(int j=0;j<12;j++){
		sprintf(Name_TDC_tof2AO_Raw,"TDC_TOF2 (%d AO)",j+1);
		sprintf(Name_TDC_tof2AI_Raw,"TDC_TOF2 (%d AI)",j+1);
		sprintf(Name_TDC_tof2BO_Raw,"TDC_TOF2 (%d BO)",j+1);
		sprintf(Name_TDC_tof2BI_Raw,"TDC_TOF2 (%d BI)",j+1);

		h_TDC_tof2AO_Raw[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof2AO_Raw);
		h_TDC_tof2AI_Raw[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof2AI_Raw);
		h_TDC_tof2BO_Raw[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof2BO_Raw);
		h_TDC_tof2BI_Raw[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof2BI_Raw);
	}
	


	/*
	for(int j=0;j<12;j++){
		sprintf(Name_TDC_tof2AO_Raw,"TDC_TOF2 (Ch. %d)",j);
		sprintf(Name_TDC_tof2AI_Raw,"TDC_TOF2 (Ch. %d)",j+24);
		sprintf(Name_TDC_tof2BO_Raw,"TDC_TOF2 (Ch. %d)",j+12);
		sprintf(Name_TDC_tof2BI_Raw,"TDC_TOF2 (Ch. %d)",j+36);

		h_TDC_tof2AO_Raw[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof2AO_Raw);
		h_TDC_tof2AI_Raw[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof2AI_Raw);
		h_TDC_tof2BO_Raw[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof2BO_Raw);
		h_TDC_tof2BI_Raw[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof2BI_Raw);
	}
	*/

	//for(int j=0;j<512;j++){
	//	sprintf(Name_ADC_MWPC_Raw,"ADC_MWPC (Ch. %d)",j);	

		//h_ADC_MWPC_Raw[j]=(TH1F*)finput->FindObjectAny(Name_ADC_MWPC_Raw);
	//}

	//Canvases

	char Name_Can_ADC_High_b1[100];			char Title_Can_ADC_High_b1[100];
	char Name_Can_ADC_Low_b1[100];			char Title_Can_ADC_Low_b1[100];
	char Name_Can_TDC_LE_b1[100];			char Title_Can_TDC_LE_b1[100];
	char Name_Can_TDC_TE_b1[100];			char Title_Can_TDC_TE_b1[100];
	char Name_Can_ADC_High2_b1[100];		char Title_Can_ADC_High2_b1[100];
	char Name_Can_ADC_Low2_b1[100];			char Title_Can_ADC_Low2_b1[100];
	char Name_Can_TDC_LE2_b1[100];			char Title_Can_TDC_LE2_b1[100];
	char Name_Can_TDC_TE2_b1[100];			char Title_Can_TDC_TE2_b1[100];

	char Name_Can_ADC_High_b2[100];			char Title_Can_ADC_High_b2[100];
	char Name_Can_ADC_Low_b2[100];			char Title_Can_ADC_Low_b2[100];
	char Name_Can_TDC_LE_b2[100];			char Title_Can_TDC_LE_b2[100];
	char Name_Can_TDC_TE_b2[100];			char Title_Can_TDC_TE_b2[100];	
	char Name_Can_ADC_High2_b2[100];		char Title_Can_ADC_High2_b2[100];
	char Name_Can_ADC_Low2_b2[100];			char Title_Can_ADC_Low2_b2[100];
	char Name_Can_TDC_LE2_b2[100];			char Title_Can_TDC_LE2_b2[100];
	char Name_Can_TDC_TE2_b2[100];			char Title_Can_TDC_TE2_b2[100];


	char Name_Can_ADC_High_b3[100];			char Title_Can_ADC_High_b3[100];
	char Name_Can_ADC_Low_b3[100];			char Title_Can_ADC_Low_b3[100];
	char Name_Can_TDC_LE_b3[100];			char Title_Can_TDC_LE_b3[100];
	char Name_Can_TDC_TE_b3[100];			char Title_Can_TDC_TE_b3[100];
	char Name_Can_ADC_High2_b3[100];		char Title_Can_ADC_High2_b3[100];
	char Name_Can_ADC_Low2_b3[100];			char Title_Can_ADC_Low2_b3[100];
	char Name_Can_TDC_LE2_b3[100];			char Title_Can_TDC_LE2_b3[100];
	char Name_Can_TDC_TE2_b3[100];			char Title_Can_TDC_TE2_b3[100];

	char Name_Can_ADC_High_b4[100];			char Title_Can_ADC_High_b4[100];
	char Name_Can_ADC_Low_b4[100];			char Title_Can_ADC_Low_b4[100];
	char Name_Can_TDC_LE_b4[100];			char Title_Can_TDC_LE_b4[100];
	char Name_Can_TDC_TE_b4[100];			char Title_Can_TDC_TE_b4[100];
	char Name_Can_ADC_High2_b4[100];		char Title_Can_ADC_High2_b4[100];
	char Name_Can_ADC_Low2_b4[100];			char Title_Can_ADC_Low2_b4[100];
	char Name_Can_TDC_LE2_b4[100];			char Title_Can_TDC_LE2_b4[100];
	char Name_Can_TDC_TE2_b4[100];			char Title_Can_TDC_TE2_b4[100];

	char Name_Can_ADC_High_b5[100];			char Title_Can_ADC_High_b5[100];
	char Name_Can_ADC_Low_b5[100];			char Title_Can_ADC_Low_b5[100];
	char Name_Can_TDC_LE_b5[100];			char Title_Can_TDC_LE_b5[100];
	char Name_Can_TDC_TE_b5[100];			char Title_Can_TDC_TE_b5[100];
	char Name_Can_ADC_High2_b5[100];		char Title_Can_ADC_High2_b5[100];
	char Name_Can_ADC_Low2_b5[100];			char Title_Can_ADC_Low2_b5[100];	
	char Name_Can_TDC_LE2_b5[100];			char Title_Can_TDC_LE2_b5[100];
	char Name_Can_TDC_TE2_b5[100];			char Title_Can_TDC_TE2_b5[100];

	char Name_Can_ADC_High_b6[100];			char Title_Can_ADC_High_b6[100];
	char Name_Can_ADC_Low_b6[100];			char Title_Can_ADC_Low_b6[100];
	char Name_Can_TDC_LE_b6[100];			char Title_Can_TDC_LE_b6[100];
	char Name_Can_TDC_TE_b6[100];			char Title_Can_TDC_TE_b6[100];
	char Name_Can_ADC_High2_b6[100];		char Title_Can_ADC_High2_b6[100];
	char Name_Can_ADC_Low2_b6[100];			char Title_Can_ADC_Low2_b6[100];
	char Name_Can_TDC_LE2_b6[100];			char Title_Can_TDC_LE2_b6[100];
	char Name_Can_TDC_TE2_b6[100];			char Title_Can_TDC_TE2_b6[100];

	char Name_Can_ADC_tof1_Raw[100];			char Title_Can_ADC_tof1_Raw[100];
	char Name_Can_ADC_tof1B_Raw[100];			char Title_Can_ADC_tof1B_Raw[100];
	char Name_Can_ADC_tof2_Raw[100];			char Title_Can_ADC_tof2_Raw[100];
	char Name_Can_ADC_tof2B_Raw[100];			char Title_Can_ADC_tof2B_Raw[100];
	char Name_Can_ADC_tof2C_Raw[100];			char Title_Can_ADC_tof2C_Raw[100];
	char Name_Can_ADC_tof2D_Raw[100];			char Title_Can_ADC_tof2D_Raw[100];
	char Name_Can_ADC_MWPC_Raw[100];			char Title_Can_ADC_MWPC_Raw[100];
	char Name_Can_ADC_MWPC2_Raw[100];			char Title_Can_ADC_MWPC2_Raw[100];


	sprintf(Name_Can_ADC_High_b1,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board 1",Run_Number);
	sprintf(Title_Can_ADC_High_b1,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board 1",Run_Number);
	sprintf(Name_Can_ADC_High2_b1,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board 1",Run_Number);
	sprintf(Title_Can_ADC_High2_b1,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board 1",Run_Number);

	sprintf(Name_Can_ADC_Low_b1,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board 1",Run_Number);
	sprintf(Title_Can_ADC_Low_b1,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board 1",Run_Number);
	sprintf(Name_Can_ADC_Low2_b1,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board 1",Run_Number);
	sprintf(Title_Can_ADC_Low2_b1,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board 1",Run_Number);

	sprintf(Name_Can_TDC_LE_b1,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board 1",Run_Number);
	sprintf(Title_Can_TDC_LE_b1,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  | EASIROC  Board 1",Run_Number);
	sprintf(Name_Can_TDC_LE2_b1,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board 1",Run_Number);
	sprintf(Title_Can_TDC_LE2_b1,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  | EASIROC  Board 1",Run_Number);

	sprintf(Name_Can_TDC_TE_b1,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board 1",Run_Number);
	sprintf(Title_Can_TDC_TE_b1,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board 1",Run_Number);
	sprintf(Name_Can_TDC_TE2_b1,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board 1",Run_Number);
	sprintf(Title_Can_TDC_TE2_b1,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board 1",Run_Number);


	sprintf(Name_Can_ADC_High_b2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95 (SFT)  |  EASIROC Board 2",Run_Number);
	sprintf(Title_Can_ADC_High_b2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95 (SFT)  |  EASIROC Board 2",Run_Number);
	sprintf(Name_Can_ADC_High2_b2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board 2",Run_Number);
	sprintf(Title_Can_ADC_High2_b2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board 2",Run_Number);

	sprintf(Name_Can_ADC_Low_b2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 64 - 95  (SFT)  |  EASIROC Board 2",Run_Number);
	sprintf(Title_Can_ADC_Low_b2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 64 - 95  (SFT)  |  EASIROC Board 2",Run_Number);
	sprintf(Name_Can_ADC_Low2_b2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board 2",Run_Number);
	sprintf(Title_Can_ADC_Low2_b2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board 2",Run_Number);

	sprintf(Name_Can_TDC_LE_b2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (SFT)  |  EASIROC Board 2",Run_Number);
	sprintf(Title_Can_TDC_LE_b2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (SFT)  |  EASIROC Board 2",Run_Number);
	sprintf(Name_Can_TDC_LE2_b2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board 2",Run_Number);
	sprintf(Title_Can_TDC_LE2_b2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board 2",Run_Number);

	sprintf(Name_Can_TDC_TE_b2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (SFT)  |  EASIROC Board 2",Run_Number);
	sprintf(Title_Can_TDC_TE_b2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (SFT)  |  EASIROC Board 2",Run_Number);
	sprintf(Name_Can_TDC_TE2_b2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board 2",Run_Number);
	sprintf(Title_Can_TDC_TE2_b2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board 2",Run_Number);


	sprintf(Name_Can_ADC_High_b3,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 32  (TARGET)  |  EASIROC Board 3",Run_Number);
	sprintf(Title_Can_ADC_High_b3,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 32  (TARGET)  |  EASIROC Board 3",Run_Number);
	sprintf(Name_Can_ADC_High2_b3,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board 3",Run_Number);
	sprintf(Title_Can_ADC_High2_b3,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board 3",Run_Number);

	sprintf(Name_Can_ADC_Low_b3,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board 3",Run_Number);
	sprintf(Title_Can_ADC_Low_b3,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board 3",Run_Number);
	sprintf(Name_Can_ADC_Low2_b3,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board 3",Run_Number);
	sprintf(Title_Can_ADC_Low2_b3,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board 3",Run_Number);

	sprintf(Name_Can_TDC_LE_b3,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board 3",Run_Number);
	sprintf(Title_Can_TDC_LE_b3,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board 3",Run_Number);
	sprintf(Name_Can_TDC_LE2_b3,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board 3",Run_Number);
	sprintf(Title_Can_TDC_LE2_b3,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board 3",Run_Number);

	sprintf(Name_Can_TDC_TE_b3,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board 3",Run_Number);
	sprintf(Title_Can_TDC_TE_b3,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board 3",Run_Number);
	sprintf(Name_Can_TDC_TE2_b3,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board 3",Run_Number);
	sprintf(Title_Can_TDC_TE2_b3,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board 3",Run_Number);


	sprintf(Name_Can_ADC_High_b4,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board 4",Run_Number);
	sprintf(Title_Can_ADC_High_b4,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board 4",Run_Number);
	sprintf(Name_Can_ADC_High2_b4,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board 4",Run_Number);
	sprintf(Title_Can_ADC_High2_b4,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board 4",Run_Number);

	sprintf(Name_Can_ADC_Low_b4,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board 4",Run_Number);
	sprintf(Title_Can_ADC_Low_b4,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board 4",Run_Number);
	sprintf(Name_Can_ADC_Low2_b4,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board 4",Run_Number);
	sprintf(Title_Can_ADC_Low2_b4,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board 4",Run_Number);

	sprintf(Name_Can_TDC_LE_b4,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board 4",Run_Number);
	sprintf(Title_Can_TDC_LE_b4,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board 4",Run_Number);
	sprintf(Name_Can_TDC_LE2_b4,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board 4",Run_Number);
	sprintf(Title_Can_TDC_LE2_b4,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board 4",Run_Number);

	sprintf(Name_Can_TDC_TE_b4,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board 4",Run_Number);	
	sprintf(Title_Can_TDC_TE_b4,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board 4",Run_Number);
	sprintf(Name_Can_TDC_TE2_b4,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board 4",Run_Number);
	sprintf(Title_Can_TDC_TE2_b4,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board 4",Run_Number);


	sprintf(Name_Can_ADC_High_b5,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board 5",Run_Number);
	sprintf(Title_Can_ADC_High_b5,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board 5",Run_Number);
	sprintf(Name_Can_ADC_High2_b5,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board 5",Run_Number);
	sprintf(Title_Can_ADC_High2_b5,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board 5",Run_Number);

	sprintf(Name_Can_ADC_Low_b5,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board 5",Run_Number);
	sprintf(Title_Can_ADC_Low_b5,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board 5",Run_Number);
	sprintf(Name_Can_ADC_Low2_b5,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board 5",Run_Number);
	sprintf(Title_Can_ADC_Low2_b5,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board 5",Run_Number);

	sprintf(Name_Can_TDC_LE_b5,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board 5",Run_Number);
	sprintf(Title_Can_TDC_LE_b5,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board 5",Run_Number);
	sprintf(Name_Can_TDC_LE2_b5,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board 5",Run_Number);
	sprintf(Title_Can_TDC_LE2_b5,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board 5",Run_Number);

	sprintf(Name_Can_TDC_TE_b5,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board 5",Run_Number);	
	sprintf(Title_Can_TDC_TE_b5,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board 5",Run_Number);
	sprintf(Name_Can_TDC_TE2_b5,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board 5",Run_Number);
	sprintf(Title_Can_TDC_TE2_b5,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board 5",Run_Number);


	sprintf(Name_Can_ADC_High_b6,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board 6",Run_Number);
	sprintf(Title_Can_ADC_High_b6,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board 6",Run_Number);
	sprintf(Name_Can_ADC_High2_b6,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board 6",Run_Number);
	sprintf(Title_Can_ADC_High2_b6,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board 6",Run_Number);

	sprintf(Name_Can_ADC_Low_b6,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board 6",Run_Number);
	sprintf(Title_Can_ADC_Low_b6,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board 6",Run_Number);
	sprintf(Name_Can_ADC_Low2_b6,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board 6",Run_Number);
	sprintf(Title_Can_ADC_Low2_b6,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board 6",Run_Number);

	sprintf(Name_Can_TDC_LE_b6,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board 6",Run_Number);
	sprintf(Title_Can_TDC_LE_b6,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board 6",Run_Number);
	sprintf(Name_Can_TDC_LE2_b6,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board 6",Run_Number);
	sprintf(Title_Can_TDC_LE2_b6,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board 6",Run_Number);

	sprintf(Name_Can_TDC_TE_b6,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board 6",Run_Number);
	sprintf(Title_Can_TDC_TE_b6,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board 6",Run_Number);
	sprintf(Name_Can_TDC_TE2_b6,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board 6",Run_Number);
	sprintf(Title_Can_TDC_TE2_b6,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board 6",Run_Number);


	sprintf(Name_Can_ADC_tof1_Raw,"TOF1 ADCs & TDCs -- Run %d  (Ch. 0 - 23)",Run_Number);
	sprintf(Title_Can_ADC_tof1_Raw,"TOF1 ADCs & TDCs -- Run %d  (Ch. 0 - 23)",Run_Number);

	sprintf(Name_Can_ADC_tof1B_Raw,"TOF1 ADCs & TDCs -- Run %d  (Ch. 0 - 23) -- Mapping Sorted",Run_Number);
	sprintf(Title_Can_ADC_tof1B_Raw,"TOF1 ADCs & TDCs -- Run %d  (Ch. 0 - 23) -- Mapping Sorted",Run_Number);

	sprintf(Name_Can_ADC_tof2_Raw,"TOF2 ADCs -- Run %d  (Ch. 0 - 55)",Run_Number);
	sprintf(Title_Can_ADC_tof2_Raw,"TOF2 ADCs -- Run %d  (Ch. 0 - 55)",Run_Number);
	sprintf(Name_Can_ADC_tof2D_Raw,"TOF2 TDCs -- Run %d  (Ch. 0 - 47)",Run_Number);
	sprintf(Title_Can_ADC_tof2D_Raw,"TOF2 TDCs -- Run %d  (Ch. 0 - 47)",Run_Number);

	sprintf(Name_Can_ADC_tof2B_Raw,"TOF2 ADCs & TDCs -- Run %d  (Ch. 0 - 23) -- Mapping Sorted",Run_Number);
	sprintf(Title_Can_ADC_tof2B_Raw,"TOF2 ADCs & TDCs -- Run %d  (Ch. 0 - 23) -- Mapping Sorted",Run_Number);
	sprintf(Name_Can_ADC_tof2C_Raw,"TOF2 ADCs & TDCs -- Run %d  (Ch. 24 - 47) -- Mapping Sorted",Run_Number);
	sprintf(Title_Can_ADC_tof2C_Raw,"TOF2 ADCs & TDCs -- Run %d  (Ch. 24 - 47) -- Mapping Sorted",Run_Number);

	sprintf(Name_Can_ADC_MWPC_Raw,"MWPC ADCs -- Run %d  (Ch. 0 - 255)",Run_Number);
	sprintf(Title_Can_ADC_MWPC_Raw,"MWPC ADCs -- Run %d  (Ch. 0 - 255)",Run_Number);

	sprintf(Name_Can_ADC_MWPC2_Raw,"MWPC ADCs -- Run %d  (Ch. 256 - 511)",Run_Number);
	sprintf(Title_Can_ADC_MWPC2_Raw,"MWPC ADCs -- Run %d  (Ch. 256 - 511)",Run_Number);


	//TCanvas *c1_b1;
	//c1_b1 = new TCanvas(Name_Can_ADC_High2_b1,Title_Can_ADC_High2_b1,1200,500); 
	//c1_b1->Divide(8,8);
	//c1_b1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	//for(Int_t ican_b1=0; ican_b1<32; ican_b1){
	//	c1_b1->cd(ican_b1+1)->SetLogy();
	//	h_ADC_High_Raw[ican_b1+32]->Draw();

		//c4->cd(ican+1+32);
		//h_TDC_LE[ican+32]->SetAxisRange(0, TDC_Graph_xmax,"X"); 

		//if(TDC_Graph_ymax > 0) {  
		//h_TDC_LE[ican+32]->SetAxisRange(0, TDC_Graph_ymax,"Y");
		//h_TDC_LE[ican+32]->Draw();
	//}


 
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//TDC_TOF_HIT_Distribution

	Int_t TOF1_Mapping1[12] = {1,3,5,7,9,11,13,15,17,19,21,23};
	Int_t TOF1_Mapping2[12] = {2,4,6,8,10,12,14,16,18,20,22,24};

	Int_t TOF2_Mapping1[12] = {1,5,9,13,17,21,25,29,33,37,41,45};
	Int_t TOF2_Mapping2[12] = {2,6,10,14,18,22,26,30,34,38,42,46};
	Int_t TOF2_Mapping3[12] = {3,7,11,15,19,23,27,31,35,39,43,47};
	Int_t TOF2_Mapping4[12] = {4,8,12,16,20,24,28,32,36,40,44,48};	

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

	int TOF1CounterA[24] = {0};
	int TOF1CounterB[24] = {0};
	int TOF1CounterC[24] = {0};
	int TOF1CounterD[24] = {0};

	int TOF2CounterA[48] = {0};
	int TOF2CounterB[48] = {0};
	int TOF2CounterC[48] = {0};
	int TOF2CounterD[48] = {0};

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//MWPC_Event_Display

	//C2 Counters
	Int_t C2YL[16] = {0};
	Int_t C2YR[16] = {0};
	Int_t C2X4L[16] = {0};
	Int_t C2X4R[16] = {0};
	Int_t C2X3L[16] = {0};
	Int_t C2X3R[16] = {0};
	Int_t C2X2L[16] = {0};
	Int_t C2X2R[16] = {0};
	Int_t C2X1L[16] = {0};
	Int_t C2X1R[16] = {0};

	//C3 Counters
	Int_t C3YL[16] = {0};
	Int_t C3YR[16] = {0};
	Int_t C3X4L[16] = {0};
	Int_t C3X4R[16] = {0};
	Int_t C3X3L[16] = {0};
	Int_t C3X3R[16] = {0};
	Int_t C3X2L[16] = {0};
	Int_t C3X2R[16] = {0};
	Int_t C3X1L[16] = {0};
	Int_t C3X1R[16] = {0};

	//C4 Counters
	Int_t C4YL[16] = {0};
	Int_t C4YR[16] = {0};
	Int_t C4X5L[16] = {0};
	Int_t C4X5R[16] = {0};
	Int_t C4X4L[16] = {0};
	Int_t C4X4R[16] = {0};
	Int_t C4X3L[16] = {0};
	Int_t C4X3R[16] = {0};
	Int_t C4X2L[16] = {0};
	Int_t C4X2R[16] = {0};
	Int_t C4X1L[16] = {0};
	Int_t C4X1R[16] = {0};

	Int_t C4X5L_sum = 0;
	Int_t C4X5LE_sum = 0;
	Int_t C4X5R_sum = 0;
	Int_t C4X5RE_sum = 0;
	Int_t C4X5L_count = 0;
	Int_t C4X5LE_count = 0;
	Int_t C4X5R_count = 0;
	Int_t C4X5RE_count = 0;

	Int_t max_C4XL = 0;
	Int_t max_C4XR = 0;
	Int_t max_C3XL = 0;
	Int_t max_C3XR = 0;
	Int_t max_C2XL = 0;
	Int_t max_C2XR = 0;

	TH1D *h_C2YL;	char Title_C2YL[100];	char Name_C2YL[100];
	TH1D *h_C2YR;	char Title_C2YR[100];	char Name_C2YR[100];
	TH1D *h_C2XL[4];	char Title_C2XL[4][100];	char Name_C2XL[4][100];
	TH1D *h_C2XR[4];	char Title_C2XR[4][100];	char Name_C2XR[4][100];

	TH1D *h_C3YL;	char Title_C3YL[100];	char Name_C3YL[100];
	TH1D *h_C3YR;	char Title_C3YR[100];	char Name_C3YR[100];
	TH1D *h_C3XL[4];	char Title_C3XL[4][100];	char Name_C3XL[4][100];
	TH1D *h_C3XR[4];	char Title_C3XR[4][100];	char Name_C3XR[4][100];

	TH1D *h_C4YL;	char Title_C4YL[100];	char Name_C4YL[100];
	TH1D *h_C4YR;	char Title_C4YR[100];	char Name_C4YR[100];
	TH1D *h_C4X5L;	char Title_C4X5L[100];	char Name_C4X5L[100];
	TH1D *h_C4X5R;	char Title_C4X5R[100];	char Name_C4X5R[100];
	TH1D *h_C4XL[4];	char Title_C4XL[4][100];	char Name_C4XL[4][100];
	TH1D *h_C4XR[4];	char Title_C4XR[4][100];	char Name_C4XR[4][100];

	TH1D *h_C2XL_M;	char Title_C2XL_M[100];	char Name_C2XL_M[100];
	TH1D *h_C2XR_M;	char Title_C2XR_M[100];	char Name_C2XR_M[100];
	TH1D *h_C3XL_M;	char Title_C3XL_M[100];	char Name_C3XL_M[100];
	TH1D *h_C3XR_M;	char Title_C3XR_M[100];	char Name_C3XR_M[100];
	TH1D *h_C4XL_M;	char Title_C4XL_M[100];	char Name_C4XL_M[100];
	TH1D *h_C4XR_M;	char Title_C4XR_M[100];	char Name_C4XR_M[100];

	for (int i=1; i<5; i++) {
		sprintf(Title_C2XL[i-1], "MWPC - C2X%dL -- Run %d | ADC > %d (Gaps 7-12)", i, Run_Number, MWPC_thr);
		sprintf(Name_C2XL[i-1], "Run %d", Run_Number);
		sprintf(Title_C2XR[i-1], "MWPC - C2X%dR -- Run %d | ADC > %d (Gaps 1-6)", i, Run_Number, MWPC_thr);
		sprintf(Name_C2XR[i-1], "Run %d", Run_Number);

		sprintf(Title_C3XL[i-1], "MWPC - C3X%dL -- Run %d | ADC > %d (Gaps 7-12)", i, Run_Number, MWPC_thr);
		sprintf(Name_C3XL[i-1], "Run %d", Run_Number);
		sprintf(Title_C3XR[i-1], "MWPC - C3X%dR -- Run %d | ADC > %d (Gaps 1-6)", i, Run_Number, MWPC_thr);
		sprintf(Name_C3XR[i-1], "Run %d", Run_Number);

		sprintf(Title_C4XL[i-1], "MWPC - C4X%dL -- Run %d | ADC > %d (Gaps 7-12)", i, Run_Number, MWPC_thr);
		sprintf(Name_C4XL[i-1], "Run %d", Run_Number);
		sprintf(Title_C4XR[i-1], "MWPC - C4X%dR -- Run %d | ADC > %d (Gaps 1-6)", i, Run_Number, MWPC_thr);
		sprintf(Name_C4XR[i-1], "Run %d", Run_Number);
	}

	sprintf(Title_C2YL, "MWPC - C2YL -- Run %d | ADC > %d (Gaps 7-12)", Run_Number, MWPC_thr);
	sprintf(Name_C2YL, "Run %d", Run_Number);
	sprintf(Title_C2YR, "MWPC - C2YR -- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);
	sprintf(Name_C2YR, "Run %d", Run_Number);

	sprintf(Title_C3YL, "MWPC - C3YL -- Run %d | ADC > %d (Gaps 7-12)", Run_Number, MWPC_thr);
	sprintf(Name_C3YL, "Run %d", Run_Number);
	sprintf(Title_C3YR, "MWPC - C3YR -- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);
	sprintf(Name_C3YR, "Run %d", Run_Number);

	sprintf(Title_C4YL, "MWPC - C4YL -- Run %d | ADC > %d (Gaps 7-12)", Run_Number, MWPC_thr);
	sprintf(Name_C4YL, "Run %d", Run_Number);
	sprintf(Title_C4YR, "MWPC - C4YR -- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);
	sprintf(Name_C4YR, "Run %d", Run_Number);

	sprintf(Title_C4X5L, "MWPC - C4X5L -- Run %d | ADC > %d (Gaps 7-12)" , Run_Number, MWPC_thr);
	sprintf(Name_C4X5L, "Run %d", Run_Number);
	sprintf(Title_C4X5R, "MWPC - C4X5R -- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);
	sprintf(Name_C4X5R, "Run %d", Run_Number);

	sprintf(Title_C2XL_M, "MWPC - C2XL -- Run %d | ADC > %d (Gaps 7-12)", Run_Number, MWPC_thr);
	sprintf(Name_C2XL_M, "Run %d", Run_Number);
	sprintf(Title_C2XR_M, "MWPC - C2XR -- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);
	sprintf(Name_C2XR_M, "Run %d", Run_Number);

	sprintf(Title_C3XL_M, "MWPC - C3XL -- Run %d | ADC > %d (Gaps 7-12)", Run_Number, MWPC_thr);
	sprintf(Name_C3XL_M, "Run %d", Run_Number);
	sprintf(Title_C3XR_M, "MWPC - C3XR -- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);
	sprintf(Name_C3XR_M, "Run %d", Run_Number);

	sprintf(Title_C4XL_M, "MWPC - C4XL -- Run %d | ADC > %d (Gaps 7-12)", Run_Number, MWPC_thr);
	sprintf(Name_C4XL_M, "Run %d", Run_Number);
	sprintf(Title_C4XR_M, "MWPC - C4XR -- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);
	sprintf(Name_C4XR_M, "Run %d", Run_Number);

	for (int i=0; i<4; i++) {
		h_C2XL[i] = new TH1D(Name_C2XL[i], Title_C2XL[i], 18,-1,17);
		h_C2XR[i] = new TH1D(Name_C2XR[i], Title_C2XR[i], 18,-1,17);

		h_C3XL[i] = new TH1D(Name_C3XL[i], Title_C3XL[i], 18,-1,17);
		h_C3XR[i] = new TH1D(Name_C3XL[i], Title_C3XR[i], 18,-1,17);

		h_C4XL[i] = new TH1D(Name_C4XL[i], Title_C4XL[i], 18,-1,17);
		h_C4XR[i] = new TH1D(Name_C4XR[i], Title_C4XR[i], 18,-1,17);
	}

	h_C2YL = new TH1D(Name_C2YL, Title_C2YL, 18,-1,17);
	h_C2YR = new TH1D(Name_C2YR, Title_C2YR, 18,-1,17);

	h_C3YL = new TH1D(Name_C3YL, Title_C3YL, 18,-1,17);
	h_C3YR = new TH1D(Name_C3YR, Title_C3YR, 18,-1,17);

	h_C4YL = new TH1D(Name_C4YL, Title_C4YL, 18,-1,17);
	h_C4YR = new TH1D(Name_C4YR, Title_C4YR, 18,-1,17);

	h_C4X5L = new TH1D(Name_C4X5L, Title_C4X5L, 18,-1,17);
	h_C4X5R = new TH1D(Name_C4X5R, Title_C4X5R, 18,-1,17);

	h_C2XL_M = new TH1D(Name_C2XL_M, Title_C2XL_M, 69,-1,68);
	h_C2XR_M = new TH1D(Name_C2XR_M, Title_C2XR_M, 69,-1,68);
	h_C3XL_M = new TH1D(Name_C3XL_M, Title_C3XL_M, 69,-1,68);
	h_C3XR_M = new TH1D(Name_C3XR_M, Title_C3XR_M, 69,-1,68);
	h_C4XL_M = new TH1D(Name_C4XL_M, Title_C4XL_M, 86,-1,85);
	h_C4XR_M = new TH1D(Name_C4XR_M, Title_C4XR_M, 86,-1,85);



	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Int_t nentries;
	Int_t nentries = (Int_t)fChain->GetEntries();//		Int_t nentries_SFT = (Int_t)fChain->GetEntries();
	//if(nentries_TARGET==nentries_SFT) nentries = nentries_TARGET;
	cout << "  " << endl;
	cout << "****  Number of events: " << nentries << "  **** " <<endl;
	cout << "  " << endl;


	cout << "   " << endl;

	if (flag!=0) nentries=flag;

	for(Int_t i=0; i<nentries; i++){
	//for(Int_t i=106; i<107; i++){
		fChain->GetEntry(i);	//	fChain_SFT->GetEntry(i);

		if(i%10000==1)	cout<<"**** "<<i<<" events done"<<endl;
		//	cout << " " << endl;
		//	cout << "********* Event " << i << endl;
  	
  		for (int j_TARGET=0; j_TARGET<256; j_TARGET++){
  			ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-par_temp_TARGET[1][j_TARGET];
			ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-par_temp_TARGET[1][j_TARGET];
			TDC_LE_TARGET[j_TARGET]=tdc_le_target[j_TARGET][0];
			//TDC_TE_TARGET[j_TARGET]=tdc_te_target[j_TARGET][0];
	  	}	

  
 		for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
			//ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-par_temp_SFT[1][j_SFT];
			//ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-par_temp_SFT[1][j_SFT];
			//TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
			//TDC_TE_SFT[j_SFT]=tdc_te_sft[j_SFT][0];
		}

		///for(int i=0; i<12; i++){
		///	ADC_TOF1[i] = ADC_tof1U[i]-par_temp_TOF1[1][i];
		///	ADC_TOF1[i+12] = ADC_tof1D[i]-par_temp_TOF1[1][i+12];
		///}


		//for (Int_t j_TOF1=0; j_TOF1<24; j_TOF1++) {
		//	ADC_TOF1[j_TOF1] = ADC_tof1[j_TOF1]-par_temp_TOF1[1][j_TOF1];
		//}


		///for(int i=0; i<12; i++){
		///	ADC_TOF2[i] = ADC_tof2AO[i]-par_temp_TOF2[1][i];
		///	ADC_TOF2[i+24] = ADC_tof2AI[i]-par_temp_TOF2[1][i+24];
		///	ADC_TOF2[i+12] = ADC_tof2BO[i]-par_temp_TOF2[1][i+12];
		///	ADC_TOF2[i=36] = ADC_tof2BI[i]-par_temp_TOF2[1][i+36];
		///}




		//for (Int_t j_TOF2=0; j_TOF2<56; j_TOF2++) {
		//	ADC_TOF2[j_TOF2] = ADC_tof2[j_TOF2]-par_temp_TOF2[1][j_TOF2];
		//}

		for (Int_t j_MWPC=0; j_MWPC<512; j_MWPC++) {
			MWPCADC[j_MWPC] = MwpcADC[j_MWPC]-MWPC_thr;
		}
  

		for (int j_TARGET=0; j_TARGET<256; j_TARGET++){
  			if ((ADC_High_TARGET[j_TARGET]>ADC_cut_TARGET2) && (TDC_LE_TARGET[j_TARGET] > TDC_thr_low) && (TDC_LE_TARGET[j_TARGET] < TDC_thr_high)) TARGET_hit_counter[j_TARGET]++;
  		}	


  		//////////////////////////////////////////////////////////////////////////////////

  		int max_index = 0;
		int max_index2 = 0;
		int max_index3 = 0;
		int max_index4 = 0;

		int max_ADC = -1000000000;
		int max_ADC2 = -1000000000;
		int max_ADC3 = -1000000000;
		int max_ADC4 = -1000000000;

		int has_data = 0;

		for(Int_t q=0; q<256; q++){
			if(ADC_Low_TARGET[q]>0) {
			//	cout << "xxxxxx   " << q << "    " << ADC_Low_TARGET[q] << endl;
				has_data++;
			}
		}

		//// Determine four largest ADC LG TARGET hits	

		for(Int_t q=0; q<256; q++){
			if(ADC_Low_TARGET[q]>max_ADC) {
				max_index = q;
				max_ADC = ADC_Low_TARGET[q];
			}
		}

		for(Int_t q=0; q<256; q++){
			if (q == max_index) continue;		
			else {
				if(ADC_Low_TARGET[q]>max_ADC2) {
					max_index2 = q;
					max_ADC2 = ADC_Low_TARGET[q];
				}
			}
		}

		for(Int_t q=0; q<256; q++){
			if ((q == max_index) || (q == max_index2)) continue;		
			else {
				if(ADC_Low_TARGET[q]>max_ADC3) {
					max_index3 = q;
					max_ADC3 = ADC_Low_TARGET[q];
				}
			}
		}

		for(Int_t q=0; q<256; q++){
			if ((q == max_index) || (q == max_index2) || (q == max_index3)) continue;		
			else {
				if(ADC_Low_TARGET[q]>max_ADC4) {
					max_index4 = q;
					max_ADC4 = ADC_Low_TARGET[q];
				}
			}
		}

		/// Remove outlying LG TARGET hits

		double x_cent = Xloc[max_index];
		double y_cent = Yloc[max_index];
	
		double hyp[256] = {-1};

		for(Int_t j=0; j<256; j++){
			hyp[j] = sqrt(pow(x_cent - Xloc[j],2) + pow(y_cent - Yloc[j],2));
		}

		if (has_data > 0) kaon_hit_counter[max_index]++;
		if (has_data > 1 && hyp[max_index2] < path_margin) kaon_hit_counter[max_index2]++;
		if (has_data > 2 && hyp[max_index3] < path_margin) kaon_hit_counter[max_index3]++;
		if (has_data > 3 && hyp[max_index4] < path_margin) kaon_hit_counter[max_index4]++;

		///// Time Cut

		if (time_window!=0) {

			int max_index_TDC = 0;
			int max_index2_TDC = 0;
			int max_index3_TDC = 0;
			int max_index4_TDC = 0;

			int max_ADC_TDC = -1000000000;
			int max_ADC2_TDC = -1000000000;
			int max_ADC3_TDC = -1000000000;
			int max_ADC4_TDC = -1000000000;

			int has_data_TDC = 0;

			for(Int_t q=0; q<256; q++){
				if(ADC_Low_TARGET[q]>0) {
					has_data_TDC++;
				}
			}

			//// Determine four largest ADC LG TARGET hits		

			for(Int_t q=0; q<256; q++){
				if((ADC_Low_TARGET[q]>max_ADC_TDC) && (TDC_LE_TARGET[q] > TDC_min_TARGET) && (TDC_LE_TARGET[q] < TDC_max_TARGET)) {
					max_index_TDC = q;
					max_ADC_TDC = ADC_Low_TARGET[q];
				}
			}

			for(Int_t q=0; q<256; q++){
				if (q == max_index_TDC) continue;		
				else {
					if((ADC_Low_TARGET[q]>max_ADC2_TDC) && (TDC_LE_TARGET[q] > TDC_min_TARGET) && (TDC_LE_TARGET[q] < TDC_max_TARGET)) {
						max_index2_TDC = q;
						max_ADC2_TDC = ADC_Low_TARGET[q];
					}
				}
			}

			for(Int_t q=0; q<256; q++){
				if ((q == max_index_TDC) || (q == max_index2_TDC)) continue;		
				else {
					if((ADC_Low_TARGET[q]>max_ADC3_TDC) && (TDC_LE_TARGET[q] > TDC_min_TARGET) && (TDC_LE_TARGET[q] < TDC_max_TARGET)) {
						max_index3_TDC = q;
						max_ADC3_TDC = ADC_Low_TARGET[q];
					}
				}
			}

			for(Int_t q=0; q<256; q++){
				if ((q == max_index_TDC) || (q == max_index2_TDC) || (q == max_index3_TDC)) continue;		
				else {
					if((ADC_Low_TARGET[q]>max_ADC4_TDC) && (TDC_LE_TARGET[q] > TDC_min_TARGET) && (TDC_LE_TARGET[q] < TDC_max_TARGET)) {
						max_index4_TDC = q;
						max_ADC4_TDC = ADC_Low_TARGET[q];
					}
				}
			}

			/// Remove outlying LG TARGET hits

			double x_cent_TDC = Xloc[max_index_TDC];
			double y_cent_TDC = Yloc[max_index_TDC];
	
			double hyp_TDC[256] = {-1};

			for(Int_t j=0; j<256; j++){
				hyp_TDC[j] = sqrt(pow(x_cent_TDC - Xloc[j],2) + pow(y_cent_TDC - Yloc[j],2));
			}

			if (has_data_TDC > 0) kaon_hit_counter2[max_index_TDC]++;
			if (has_data_TDC > 1 && hyp_TDC[max_index2_TDC] < path_margin) kaon_hit_counter2[max_index2_TDC]++;
			if (has_data_TDC > 2 && hyp_TDC[max_index3_TDC] < path_margin) kaon_hit_counter2[max_index3_TDC]++;
			if (has_data_TDC > 3 && hyp_TDC[max_index4_TDC] < path_margin) kaon_hit_counter2[max_index4_TDC]++;
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//TOF_GAP_Multiplicity.C

		int has_data_TOF = 0;
  	
  		for (int j_TARGET=0; j_TARGET<256; j_TARGET++){
  			ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-par_temp_TARGET[1][j_TARGET];
			ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-par_temp_TARGET[1][j_TARGET];
			TDC_LE_TARGET[j_TARGET]=tdc_le_target[j_TARGET][0];
			TDC_TE_TARGET[j_TARGET]=tdc_te_target[j_TARGET][0];
  		}	


 		for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
			ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-par_temp_SFT[1][j_SFT];
			ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-par_temp_SFT[1][j_SFT];
			TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
			TDC_TE_SFT[j_SFT]=tdc_te_sft[j_SFT][0];
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
		    TDC_TOF2AI[j_TDCTOF] = TDC_tof2AI[j_TDCTOF];
		    TDC_TOF2BO[j_TDCTOF] = TDC_tof2BO[j_TDCTOF];
		    TDC_TOF2BI[j_TDCTOF] = TDC_tof2BI[j_TDCTOF];
  		}


  		/// xxxxxxx
  		if(flag_corr==1){
  			TDC_TOF1D[2] = TDC_tof1D[3];
  			TDC_TOF1D[3] = TDC_tof1D[2];
  		}

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

  		if (time_window == 0 && good_MWPC) {
	  		for (int kk=0; kk<12; kk++) {
				if (kk == 0) {
					if ((ADC_TOF2AO[0]>0) || (ADC_TOF2AI[0]>0) || (ADC_TOF2BO[0]>0) || (ADC_TOF2BI[0]>0)) {
						if((ADC_TOF1U[0] > 0) || (ADC_TOF1D[0] > 0) || (ADC_TOF1U[1] > 0) || (ADC_TOF1D[1] > 0) || (ADC_TOF1U[11] > 0) || (ADC_TOF1D[11] > 0)) {TOF_Hit_Counter[0]++; has_data_TOF++;}
					} 
				}

				/// Replace ADC TOF2[18] with ADC TOF2[55]
				//else if (kk == 6) {
				//	if ((ADC_TOF2[6]>0) || (ADC_TOF2[30]>0) || (ADC_TOF2[55]>0) || (ADC_TOF2[42]>0)) {
				//		if((ADC_TOF1[6] > 0) || (ADC_TOF1[18] > 0) || (ADC_TOF1[7] > 0) || (ADC_TOF1[19] > 0) || (ADC_TOF1[5] > 0) || (ADC_TOF1[17] > 0)) {TOF_Hit_Counter[6]++; has_data_TOF++;}
				//	}
				//}

				else if (kk == 11) {
					if ((ADC_TOF2AO[11]>0) || (ADC_TOF2AI[11]>0) || (ADC_TOF2BO[11]>0) || (ADC_TOF2BI[11]>0)) {
						if((ADC_TOF1U[11] > 0) || (ADC_TOF1D[11] > 0) || (ADC_TOF1U[0] > 0) || (ADC_TOF1D[0] > 0) || (ADC_TOF1U[10] > 0) || (ADC_TOF1D[10] > 0)) {TOF_Hit_Counter[11]++; has_data_TOF++;}
					}
				}
		
				else {
					if ((ADC_TOF2AO[kk]>0) || (ADC_TOF2AI[kk]>0) || (ADC_TOF2BO[kk]>0) || (ADC_TOF2BI[kk]>0))  {
//					if ((ADC_TOF2BI[kk]>0))  {
						if((ADC_TOF1U[kk] > 0) || (ADC_TOF1D[kk] > 0) || (ADC_TOF1U[kk+1] > 0) || (ADC_TOF1D[kk+1] > 0) || (ADC_TOF1U[kk-1] > 0) || (ADC_TOF1D[kk-1] > 0)) {
//							cout << "XXXXXX   " << kk << "   " << "GOTCHA !!" << endl; 
							TOF_Hit_Counter[kk]++; has_data_TOF++;
						}
					}
				}
			}
		}

		
//		cout << "****  " << ADC_TOF2AI[9] << endl;

		bool has_ADC_TOF1_hit[12] = {false};
		bool has_TDC_TOF1_hit[12] = {false};
		bool has_ADC_TOF2_hit[12] = {false};
		bool has_TDC_TOF2_hit[12] = {false};


		/// Determine good gap with time window

		if (time_window != 0 && good_MWPC) {
			for (int kk=0; kk<12; kk++) {
				/// Replace ADC TOF2[18] with ADC TOF2[55]
				//if (kk == 6) {
				//	if ((ADC_TOF2[6]>0) || (ADC_TOF2[30]>0) || (ADC_TOF2[55]>0) || (ADC_TOF2[42]>0)) {has_ADC_TOF2_hit[6]=true;}
				//	if (((TDC_TOF2AO[6]>TDC_TOF2_min && TDC_TOF2AO[6] < TDC_TOF2_max) || (TDC_TOF2AI[6]>TDC_TOF2_min && TDC_TOF2AI[6] < TDC_TOF2_max)) 
				//   || ((TDC_TOF2BO[6]>TDC_TOF2_min && TDC_TOF2BO[6] < TDC_TOF2_max) || (TDC_TOF2BI[6]>TDC_TOF2_min && TDC_TOF2BI[6] < TDC_TOF2_max))) {has_TDC_TOF2_hit[6]=true;}	

				//	if (ADC_TOF1[6]>0 || ADC_TOF1[18]>0) {has_ADC_TOF1_hit[6] = true;}
				// 	if ((TDC_TOF1U[6]>TDC_TOF1_min && TDC_TOF1U[6]<TDC_TOF1_max) || (TDC_TOF1D[6]>TDC_TOF1_min && TDC_TOF1D[6]<TDC_TOF1_max)) {has_TDC_TOF1_hit[6] = true;}
			  	//}

				//else {
				if ((ADC_TOF2AO[kk]>0) || (ADC_TOF2AI[kk]>0) || (ADC_TOF2BO[kk]>0) || (ADC_TOF2BI[kk]>0)) {has_ADC_TOF2_hit[kk]=true;}
				if (((TDC_TOF2AO[kk]>TDC_TOF2_min && TDC_TOF2AO[kk] < TDC_TOF2_max) || (TDC_TOF2AI[kk]>TDC_TOF2_min && TDC_TOF2AI[kk] < TDC_TOF2_max)) 
				   || ((TDC_TOF2BO[kk]>TDC_TOF2_min && TDC_TOF2BO[kk] < TDC_TOF2_max) || (TDC_TOF2BI[kk]>TDC_TOF2_min && TDC_TOF2BI[kk] < TDC_TOF2_max))) {has_TDC_TOF2_hit[kk]=true;}

				if (ADC_TOF1U[kk]>0 || ADC_TOF1D[kk]>0) {has_ADC_TOF1_hit[kk] = true;}
				if ((TDC_TOF1U[kk]>TDC_TOF1_min && TDC_TOF1U[kk]<TDC_TOF1_max) || (TDC_TOF1D[kk]>TDC_TOF1_min && TDC_TOF1D[kk]<TDC_TOF1_max)) {has_TDC_TOF1_hit[kk] = true;}
				//}
			}
	

	
			for (int kk=0; kk<12; kk++) {
				if (kk == 0) {
					if (has_ADC_TOF2_hit[0] && has_TDC_TOF2_hit[0]) {
						if ((has_ADC_TOF1_hit[0] && has_TDC_TOF1_hit[0]) || (has_ADC_TOF1_hit[11] && has_TDC_TOF1_hit[11]) || (has_ADC_TOF1_hit[1] && has_TDC_TOF1_hit[1])) {
							TOF_Hit_Counter[0]++;
							has_data_TOF++;
						}
					}
				}
				else if (kk == 11) {
					if (has_ADC_TOF2_hit[11] && has_TDC_TOF2_hit[11]) {
						if ((has_ADC_TOF1_hit[11] && has_TDC_TOF1_hit[11]) || (has_ADC_TOF1_hit[10] && has_TDC_TOF1_hit[10]) || (has_ADC_TOF1_hit[0] && has_TDC_TOF1_hit[0])) {
							TOF_Hit_Counter[11]++;
							has_data_TOF++;
						}
					}
				}
				else {
					if (has_ADC_TOF2_hit[kk] && has_TDC_TOF2_hit[kk]) {
						if ((has_ADC_TOF1_hit[kk] && has_TDC_TOF1_hit[kk]) || (has_ADC_TOF1_hit[kk-1] && has_TDC_TOF1_hit[kk-1]) || (has_ADC_TOF1_hit[kk+1] && has_TDC_TOF1_hit[kk+1])) {
							TOF_Hit_Counter[kk]++;
							has_data_TOF++;
						}
					}
				}
			}
		}

		if (has_data_TOF > 0) good_gap_counter++;


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//TDC_TOF_HIT_Distribution

		/// Filter TOF1 hits by distribution within peak

		for (int q=0; q<12; q++) {
			if (ADC_TOF1U[q] > 0) {
				if (TDC_TOF1U[q] == -1) TOF1CounterA[q]++;
				else if (TDC_TOF1U[q] > TDC_TOF1U_min[q] && TDC_TOF1U[q] < TDC_TOF1U_max[q]) TOF1CounterB[q]++;
				else if (TDC_TOF1U[q] > 4000) TOF1CounterD[q]++;
				else {TOF1CounterC[q]++;}
			}
		}

		for (int q=0; q<12; q++) {
			if (ADC_TOF1D[q] > 0) {
				if (TDC_TOF1D[q] == -1) TOF1CounterA[q+12]++;
				else if (TDC_TOF1D[q] > TDC_TOF1D_min[q] && TDC_TOF1D[q] < TDC_TOF1D_max[q]) TOF1CounterB[q+12]++;
				else if (TDC_TOF1D[q] > 4000) TOF1CounterD[q+12]++;
				else {TOF1CounterC[q+12]++;}
			}
		}


		/// Filter TOF2 hits by distribution within peak

		for (int q=0; q<12; q++) {
			if (ADC_TOF2AO[q] > 0) {
				if (TDC_tof2AO[q] == -1) TOF2CounterA[q]++;
				else if (TDC_tof2AO[q] > TDC_TOF2AO_min[q] && TDC_tof2AO[q] < TDC_TOF2AO_max[q]) TOF2CounterB[q]++;
				else if (TDC_tof2AO[q] > 4000) TOF2CounterD[q]++;
				else {TOF2CounterC[q]++;}
			}
		}

		for (int q=0; q<12; q++) {
			if (ADC_TOF2BO[q] > 0) {
				if (TDC_tof2BO[q] == -1) TOF2CounterA[q+12]++;
				else if (TDC_tof2BO[q] > TDC_TOF2BO_min[q] && TDC_tof2BO[q] < TDC_TOF2BO_max[q]) TOF2CounterB[q+12]++;
				else if (TDC_tof2BO[q] > 4000) TOF2CounterD[q+12]++;
				else {TOF2CounterC[q+12]++;}
			}
		}

		for (int q=0; q<12; q++) {
			if (ADC_TOF2AI[q] > 0) {
				if (TDC_tof2AI[q] == -1) TOF2CounterA[q+24]++;
				else if (TDC_tof2AI[q] > TDC_TOF2AI_min[q] && TDC_tof2AI[q] < TDC_TOF2AI_max[q]) TOF2CounterB[q+24]++;
				else if (TDC_tof2AI[q] > 4000) TOF2CounterD[q+24]++;
				else {TOF2CounterC[q+24]++;}
			}
		}

		for (int q=0; q<12; q++) {
			if (ADC_TOF2BI[q] > 0) {
				if (TDC_tof2BI[q] == -1) TOF2CounterA[q+36]++;
				else if (TDC_tof2BI[q] > TDC_TOF2BI_min[q] && TDC_tof2BI[q] < TDC_TOF2BI_max[q]) TOF2CounterB[q+36]++;
				else if (TDC_tof2BI[q] > 4000) TOF2CounterD[q+36]++;
				else {TOF2CounterC[q+36]++;}
			}
		}		

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//MWPC_Event_Display

		for (int q=0; q<512; q++) {

			//C3 Counters
			if (q >= 0 && q <= 15) {
				if (MWPCADC[q]>0) C3X3L[q]++;
			}
			if (q >= 16 && q <= 31) {
				if (MWPCADC[q]>0) C3X3R[q-16]++;
			}
			if (q >= 32 && q <= 47) {
				if (MWPCADC[q]>0) C3X2L[q-32]++;
			}
			if (q >= 48 && q <= 63) {
				if (MWPCADC[q]>0) C3X2R[q-48]++;
			}
			if (q >= 64 && q <= 79) {
				if (MWPCADC[q]>0) C3X1L[q-64]++;
			}
			if (q >= 80 && q <= 95) {
				if (MWPCADC[q]>0) C3X1R[q-80]++;
			}

			//C2 Counters
			if (q >= 96 && q <= 111) {
				if (MWPCADC[q]>0) C2YL[q-96]++;
			}
			if (q >= 112 && q <= 127) {
				if (MWPCADC[q]>0) C2YR[q-112]++;
			}
			if (q >= 128 && q <= 143) {
				if (MWPCADC[q]>0) C2X4L[q-128]++;
			}
			if (q >= 144 && q <= 159) {
				if (MWPCADC[q]>0) C2X4R[q-144]++;
			}
			if (q >= 160 && q <= 175) {
				if (MWPCADC[q]>0) C2X3L[q-160]++;
			}
			if (q >= 176 && q <= 191) {
				if (MWPCADC[q]>0) C2X3R[q-176]++;
			}
			if (q >= 192 && q <= 207) {
				if (MWPCADC[q]>0) C2X2L[q-192]++;
			}
			if (q >= 208 && q <= 223) {
				if (MWPCADC[q]>0) C2X2R[q-208]++;
			}
			if (q >= 224 && q <= 239) {
				if (MWPCADC[q]>0) C2X1L[q-224]++;
			}
			if (q >= 240 && q <= 255) {
				if (MWPCADC[q]>0) C2X1R[q-240]++;
			}

			//C4 Counters
			if (q >= 256 && q <= 271) {
				if (MWPCADC[q]>0) C4YL[q-256]++;
			}
			if (q >= 272 && q <= 287) {
				if (MWPCADC[q]>0) C4YR[q-272]++;
			}
			if (q >= 288 && q <= 303) {
				if (MWPCADC[q]>0) {
					if (q >= 288 && q <= 295) {C4X5L[q-288]++; C4X5L_count++; C4X5L_sum = C4X5L_sum + MWPCADC[q];}
					else {C4X5L[q-288]++; C4X5LE_count++; C4X5LE_sum = C4X5LE_sum + MWPCADC[q];}
				}
			}
			if (q >= 304 && q <= 319) {
				if (MWPCADC[q]>0) {
					if (q >= 304 && q <= 311) {C4X5R[q-304]++; C4X5R_count++; C4X5R_sum = C4X5R_sum + MWPCADC[q];}
					else {C4X5R[q-304]++; C4X5RE_count++; C4X5RE_sum = C4X5RE_sum + MWPCADC[q];}
				}
			}
			if (q >= 320 && q <= 335) {
				if (MWPCADC[q]>0) C4X4L[q-320]++;
			}
			if (q >= 336 && q <= 351) {
				if (MWPCADC[q]>0) C4X4R[q-336]++;
			}
			if (q >= 352 && q <= 367) {
				if (MWPCADC[q]>0) C4X3L[q-352]++;
			}
			if (q >= 368 && q <= 383) {
				if (MWPCADC[q]>0) C4X3R[q-368]++;
			}	
			if (q >= 384 && q <= 399) {
				if (MWPCADC[q]>0) C4X2L[q-384]++;
			}
			if (q >= 400 && q <= 415) {
				if (MWPCADC[q]>0) C4X2R[q-400]++;
			}
			if (q >= 416 && q <= 431) {
				if (MWPCADC[q]>0) C4X1L[q-416]++;
			}
			if (q >= 432 && q <= 447) {
				if (MWPCADC[q]>0) C4X1R[q-432]++;
			}

			//C3 Counters
			if (q >= 448 && q <= 463) {
				if (MWPCADC[q]>0) C3YL[q-448]++;
			}
			if (q >= 464 && q <= 479) {
				if (MWPCADC[q]>0) C3YR[q-464]++;
			}
			if (q >= 480 && q <= 495) {
				if (MWPCADC[q]>0) C3X4L[q-480]++;
			}
		}
	
		for (int q=0; q<512; q++) {

			//C3 Counters
			if (q >= 0 && q <= 15) {
				if (MWPCADC[q]>0) C3X3L[q]++;
			}
			if (q >= 16 && q <= 31) {
				if (MWPCADC[q]>0) C3X3R[q-16]++;
			}
			if (q >= 32 && q <= 47) {
				if (MWPCADC[q]>0) C3X2L[q-32]++;
			}
			if (q >= 48 && q <= 63) {
				if (MWPCADC[q]>0) C3X2R[q-48]++;
			}
			if (q >= 64 && q <= 79) {
				if (MWPCADC[q]>0) C3X1L[q-64]++;
			}
			if (q >= 80 && q <= 95) {
				if (MWPCADC[q]>0) C3X1R[q-80]++;
			}

			//C2 Counters
			if (q >= 96 && q <= 111) {
				if (MWPCADC[q]>0) C2YL[q-96]++;
			}
			if (q >= 112 && q <= 127) {
				if (MWPCADC[q]>0) C2YR[q-112]++;
			}
			if (q >= 128 && q <= 143) {
				if (MWPCADC[q]>0) C2X4L[q-128]++;
			}
			if (q >= 144 && q <= 159) {
				if (MWPCADC[q]>0) C2X4R[q-144]++;
			}
			if (q >= 160 && q <= 175) {
				if (MWPCADC[q]>0) C2X3L[q-160]++;
			}
			if (q >= 176 && q <= 191) {
				if (MWPCADC[q]>0) C2X3R[q-176]++;
			}
			if (q >= 192 && q <= 207) {
				if (MWPCADC[q]>0) C2X2L[q-192]++;
			}
			if (q >= 208 && q <= 223) {
				if (MWPCADC[q]>0) C2X2R[q-208]++;
			}
			if (q >= 224 && q <= 239) {
				if (MWPCADC[q]>0) C2X1L[q-224]++;
			}
			if (q >= 240 && q <= 255) {
				if (MWPCADC[q]>0) C2X1R[q-240]++;
			}

			//C4 Counters
			if (q >= 256 && q <= 271) {
				if (MWPCADC[q]>0) C4YL[q-256]++;
			}
			if (q >= 272 && q <= 287) {
				if (MWPCADC[q]>0) C4YR[q-272]++;
			}
			if (q >= 288 && q <= 303) {
				if (MWPCADC[q]>0) {
					if (q >= 288 && q <= 295) {C4X5L[q-288]++; C4X5L_count++; C4X5L_sum = C4X5L_sum + MWPCADC[q];}
					else {C4X5L[q-288]++; C4X5LE_count++; C4X5LE_sum = C4X5LE_sum + MWPCADC[q];}
				}
			}
			if (q >= 304 && q <= 319) {
				if (MWPCADC[q]>0) {
					if (q >= 304 && q <= 311) {C4X5R[q-304]++; C4X5R_count++; C4X5R_sum = C4X5R_sum + MWPCADC[q];}
					else {C4X5R[q-304]++; C4X5RE_count++; C4X5RE_sum = C4X5RE_sum + MWPCADC[q];}
				}
			}
			if (q >= 320 && q <= 335) {
				if (MWPCADC[q]>0) C4X4L[q-320]++;
			}
			if (q >= 336 && q <= 351) {
				if (MWPCADC[q]>0) C4X4R[q-336]++;
			}
			if (q >= 352 && q <= 367) {
				if (MWPCADC[q]>0) C4X3L[q-352]++;
			}
			if (q >= 368 && q <= 383) {
				if (MWPCADC[q]>0) C4X3R[q-368]++;
			}
			if (q >= 384 && q <= 399) {
				if (MWPCADC[q]>0) C4X2L[q-384]++;
			}
			if (q >= 400 && q <= 415) {
				if (MWPCADC[q]>0) C4X2R[q-400]++;
			}
			if (q >= 416 && q <= 431) {
				if (MWPCADC[q]>0) C4X1L[q-416]++;
			}
			if (q >= 432 && q <= 447) {
				if (MWPCADC[q]>0) C4X1R[q-432]++;
			}

			//C3 Counters
			if (q >= 448 && q <= 463) {
				if (MWPCADC[q]>0) C3YL[q-448]++;
			}
			if (q >= 464 && q <= 479) {
				if (MWPCADC[q]>0) C3YR[q-464]++;
			}
			if (q >= 480 && q <= 495) {
				if (MWPCADC[q]>0) C3X4L[q-480]++;
			}
			if (q >= 496 && q <= 511) {
				if (MWPCADC[q]>0) C3X4R[q-496]++;
			}
		}



		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	

	} // EndLoop Over Events

	////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//BeamProfile

	TH2F *hK_HITS;
	TH2F *hK_HITS2;

	char Title_TARGET_h3D[100];	char Name_TARGET_h3D[100];
	char Title_KAON_h3D[100];	char Name_KAON_h3D[100];
	sprintf(Name_TARGET_h3D,"ADC HG TARGET Hit Position - Run %d: ADC > %d", Run_Number, ADC_cut_TARGET); 
	sprintf(Title_TARGET_h3D,"ADC HG TARGET Hit Position -  Run %d", Run_Number); 
	sprintf(Name_KAON_h3D,"ADC LG KAON Hit Position - Run %d: ADC > %d", Run_Number, ADC_cut_TARGET); 
	sprintf(Title_KAON_h3D,"ADC LG KAON Hit Position -  Run %d", Run_Number); 
	hK_HITS = new TH2F(Name_TARGET_h3D,Title_TARGET_h3D,20,-29.45,29.45,20,-29.45,29.45);
	hK_HITS2 = new TH2F(Name_KAON_h3D,Title_KAON_h3D,20,-29.45,29.45,20,-29.45,29.45);

	for(int k=0 ; k<256 ; k++) {
		//cout << "#####  " << k << "   " << kaon_hit_counter[k] << endl;

		hK_HITS->Fill(Xloc[k],Yloc[k],TARGET_hit_counter[k]);
		if (k!=102) hK_HITS2->Fill(Xloc[k],Yloc[k],kaon_hit_counter[k]);
	}


	char Name_Can_BeamProf[100];			char Title_Can_BeamProf[100];
	sprintf(Name_Can_BeamProf,"TARGET Beam Profile (ADC HG & LG) -- Run %d",Run_Number);
	sprintf(Title_Can_BeamProf,"TARGET Beam Profile (ADC HG & LG) -- Run %d",Run_Number);


	TCanvas* c4_TARGET = new TCanvas(Name_Can_BeamProf, Title_Can_BeamProf,1300,700);
	c4_TARGET->Divide(4,2);
	c4_TARGET->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	TH1D *hK_HITSX = hK_HITS->ProjectionX("XProjection");
	TH1D *hK_HITSY = hK_HITS->ProjectionY("YProjection");
	TH1D *hK_HITSX2 = hK_HITS2->ProjectionX("XProjection");
	TH1D *hK_HITSY2 = hK_HITS2->ProjectionY("YProjection");

	c4_TARGET->cd(1);
	gStyle->SetPalette(1);
	hK_HITS->SetMarkerSize(10);
	hK_HITS->SetStats(0);
	hK_HITS->Draw("COLZ");

	c4_TARGET->cd(2);
	hK_HITSX->SetTitle("TARGET ADC Hits - X Projection");
	hK_HITSX->SetStats(0);
	hK_HITSX->Draw();

	c4_TARGET->cd(3);
	hK_HITSY->SetTitle("TARGET ADC Hits - Y Projection");
	hK_HITSY->SetStats(0);
	hK_HITSY->Draw();

	c4_TARGET->cd(4);
	gStyle->SetPalette(1);
	hK_HITS->Draw("Lego2");

	c4_TARGET->cd(5);
	gStyle->SetPalette(1);
	hK_HITS2->SetMarkerSize(10);
	hK_HITS2->SetStats(0);
	hK_HITS2->Draw("COLZ");

	c4_TARGET->cd(6);
	hK_HITSX2->SetTitle("TARGET ADC Hits - X Projection");
	hK_HITSX2->SetStats(0);
	hK_HITSX2->Draw();

	c4_TARGET->cd(7);
	hK_HITSY2->SetTitle("TARGET ADC Hits - Y Projection");
	hK_HITSY2->SetStats(0);
	hK_HITSY2->Draw();

	c4_TARGET->cd(8);
	gStyle->SetPalette(1);
	hK_HITS2->Draw("Lego2");

	char pdf_open[50];
	char pdf_in[50];
	char pdf_closed[50];

	sprintf(pdf_open, "Monitoring_Report_Run%d.pdf(",Run_Number);
	sprintf(pdf_in, "Monitoring_Report_Run%d.pdf",Run_Number);
	sprintf(pdf_closed, "Monitoring_Report_Run%d.pdf)",Run_Number);

	c4_TARGET->Print(pdf_open, "pdf");

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//TOF_GAP_Multiplicity

	//double eff2 = double(double(good_gap_counter)/double(nentries))*100;

	TH1D *h_count;	char Title_h_count[100];	char Name_h_count[100];
	sprintf(Title_h_count,"TOF Gap Multiplicity -- Run %d", Run_Number); 
	sprintf(Name_h_count,"Run %d", Run_Number);
	h_count = new TH1D(Name_h_count,Title_h_count,14,0,14);

	for (int i=0; i<12; i++) {
		//cout << i << 
		h_count->Fill(i+1,TOF_Hit_Counter[i]);
	}

	char Name_Can[100];			char Title_Can[100];
	sprintf(Name_Can,"TOF Gap Multiplicity -- Run %d", Run_Number);
	sprintf(Title_Can,"TOF Gap Multiplicity -- Run %d", Run_Number);

	TCanvas *c111;
	c111 = new TCanvas(Name_Can,Title_Can,800,500); 
	c111->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	c111->cd();
	h_count->Draw(); 

	c111->Print(pdf_in, "pdf");

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//TDC_TOF_HIT_Distribution

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

	sprintf(Name_Can_ADC_tof1B,"TOF1 TDC Distributions -- Run %d",Run_Number);
	sprintf(Title_Can_ADC_tof1B,"TOF1 TDC Distributions -- Run %d",Run_Number);

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
	c7 = new TCanvas(Name_Can_ADC_tof1B,Title_Can_ADC_tof1B,1300,1000); 
	c7->Divide(6,4);
	c7->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=0; ican<12; ican++){
		c7->cd(TOF1_Mapping1[ican]);

		char TOF1Title3[250];
		sprintf(TOF1Title3, "TDC (Ch. %d) - TOF1 | Up-%d : Run %d", ican, ican+1, Run_Number);

		h_TDC_tof1B[ican]->SetTitle(TOF1Title3);
		h_TDC_tof1B[ican]->SetAxisRange(0, TOF1_height_max[ican],"Y");
		h_TDC_tof1B[ican]->SetLineColor(3);
		h_TDC_tof1B[ican]->SetFillStyle(3004);
		h_TDC_tof1B[ican]->SetFillColor(3);
		h_TDC_tof1B[ican]->SetStats(0);	
		h_TDC_tof1B[ican]->Draw();

		h_TDC_tof1A[ican]->SetTitle(TOF1Title3);
		h_TDC_tof1A[ican]->SetAxisRange(0, TOF1_height_max[ican],"Y");
		h_TDC_tof1A[ican]->SetLineColor(2);
		h_TDC_tof1A[ican]->SetFillStyle(3004);
		h_TDC_tof1A[ican]->SetFillColor(2);
		h_TDC_tof1A[ican]->SetStats(0);	
		h_TDC_tof1A[ican]->Draw("same");

		h_TDC_tof1C[ican]->SetTitle(TOF1Title3);
		h_TDC_tof1C[ican]->SetAxisRange(0, TOF1_height_max[ican],"Y");
		h_TDC_tof1C[ican]->SetLineColor(4);
		h_TDC_tof1C[ican]->SetFillStyle(3004);
		h_TDC_tof1C[ican]->SetFillColor(4);
		h_TDC_tof1C[ican]->SetStats(0);	
		h_TDC_tof1C[ican]->Draw("same");

		h_TDC_tof1D[ican]->SetTitle(TOF1Title3);
		h_TDC_tof1D[ican]->SetAxisRange(0, TOF1_height_max[ican],"Y");
		h_TDC_tof1D[ican]->SetLineColor(1);
		h_TDC_tof1D[ican]->SetFillStyle(3004);
		h_TDC_tof1D[ican]->SetFillColor(1);
		h_TDC_tof1D[ican]->SetStats(0);	
		h_TDC_tof1D[ican]->Draw("same");

		l1->Draw("same");
	}

	for(Int_t ican=0; ican<12; ican++){
		c7->cd(TOF1_Mapping2[ican]);

		char TOF1Title4[250];
		sprintf(TOF1Title4, "TDC (Ch. %d) - TOF1 | Down-%d : Run %d", ican+12, ican+1, Run_Number);
		h_TDC_tof1B[ican+12]->SetTitle(TOF1Title4);
		h_TDC_tof1B[ican+12]->SetAxisRange(0, TOF1_height_max[ican+12],"Y");
		h_TDC_tof1B[ican+12]->SetLineColor(3);
		h_TDC_tof1B[ican+12]->SetFillStyle(3004);
		h_TDC_tof1B[ican+12]->SetFillColor(3);
	    h_TDC_tof1B[ican+12]->SetStats(0);	
		h_TDC_tof1B[ican+12]->Draw();

		h_TDC_tof1A[ican+12]->SetTitle(TOF1Title4);
		h_TDC_tof1A[ican+12]->SetAxisRange(0, TOF1_height_max[ican+12],"Y");
		h_TDC_tof1A[ican+12]->SetLineColor(2);
		h_TDC_tof1A[ican+12]->SetFillStyle(3004);
		h_TDC_tof1A[ican+12]->SetFillColor(2);
		h_TDC_tof1A[ican+12]->SetStats(0);	
		h_TDC_tof1A[ican+12]->Draw("same");

		h_TDC_tof1C[ican+12]->SetTitle(TOF1Title4);
		h_TDC_tof1C[ican+12]->SetAxisRange(0, TOF1_height_max[ican+12],"Y");
		h_TDC_tof1C[ican+12]->SetLineColor(4);
		h_TDC_tof1C[ican+12]->SetFillStyle(3004);
		h_TDC_tof1C[ican+12]->SetFillColor(4);
		h_TDC_tof1C[ican+12]->SetStats(0);	
		h_TDC_tof1C[ican+12]->Draw("same");

		h_TDC_tof1D[ican+12]->SetTitle(TOF1Title4);
		h_TDC_tof1D[ican+12]->SetAxisRange(0, TOF1_height_max[ican+12],"Y");
		h_TDC_tof1D[ican+12]->SetLineColor(1);
		h_TDC_tof1D[ican+12]->SetFillStyle(3004);
		h_TDC_tof1D[ican+12]->SetFillColor(1);
		h_TDC_tof1D[ican+12]->SetStats(0);	
		h_TDC_tof1D[ican+12]->Draw("same");

		l1->Draw("same");
	}


	////Draw TOF2 Histograms

	char Name_Can_ADC_tof2B[100];			char Title_Can_ADC_tof2B[100];

	sprintf(Name_Can_ADC_tof2B,"TOF2 TDC Distributions -- Run %d",Run_Number);
	sprintf(Title_Can_ADC_tof2B,"TOF2 TDC Distributions -- Run %d",Run_Number);

	TCanvas *c8;
	c8 = new TCanvas(Name_Can_ADC_tof2B,Title_Can_ADC_tof2B,1300,1000); 
	c8->Divide(8,6);
	c8->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=0; ican<12; ican++){
		c8->cd(TOF2_Mapping1[ican]);

		char TOF2Title5[250];
		sprintf(TOF2Title5, "Raw TDC (Ch. %d) - TOF2 | OutA-%d : Run %d", ican, ican+1, Run_Number);
		h_TDC_tof2B[ican]->SetTitle(TOF2Title5);
		h_TDC_tof2B[ican]->SetAxisRange(0, TOF2_height_max[ican],"Y");
		h_TDC_tof2B[ican]->SetLineColor(3);
		h_TDC_tof2B[ican]->SetFillStyle(3004);
		h_TDC_tof2B[ican]->SetFillColor(3);
		h_TDC_tof2B[ican]->SetStats(0);
		h_TDC_tof2B[ican]->Draw();

		h_TDC_tof2A[ican]->SetTitle(TOF2Title5);
		h_TDC_tof2A[ican]->SetAxisRange(0, TOF2_height_max[ican],"Y");
		h_TDC_tof2A[ican]->SetLineColor(2);
		h_TDC_tof2A[ican]->SetFillStyle(3004);
		h_TDC_tof2A[ican]->SetFillColor(2);
		h_TDC_tof2A[ican]->SetStats(0);
		h_TDC_tof2A[ican]->Draw("same");

		h_TDC_tof2C[ican]->SetTitle(TOF2Title5);
		h_TDC_tof2C[ican]->SetAxisRange(0, TOF2_height_max[ican],"Y");
		h_TDC_tof2C[ican]->SetLineColor(4);
		h_TDC_tof2C[ican]->SetFillStyle(3004);
		h_TDC_tof2C[ican]->SetFillColor(4);
		h_TDC_tof2C[ican]->SetStats(0);
		h_TDC_tof2C[ican]->Draw("same");

		h_TDC_tof2D[ican]->SetTitle(TOF2Title5);
		h_TDC_tof2D[ican]->SetAxisRange(0, TOF2_height_max[ican],"Y");
		h_TDC_tof2D[ican]->SetLineColor(1);
		h_TDC_tof2D[ican]->SetFillStyle(3004);
		h_TDC_tof2D[ican]->SetFillColor(1);
		h_TDC_tof2D[ican]->SetStats(0);
		h_TDC_tof2D[ican]->Draw("same");

		l2->Draw("same");
	}

	for(Int_t ican=0; ican<12; ican++){
		c8->cd(TOF2_Mapping2[ican]);
		
		char TOF2Title6[250];
		sprintf(TOF2Title6, "Raw TDC (Ch. %d) - TOF2 | InA-%d : Run %d", ican+24, ican+1, Run_Number);
		h_TDC_tof2B[ican+24]->SetTitle(TOF2Title6);
		h_TDC_tof2B[ican+24]->SetAxisRange(0, TOF2_height_max[ican+24],"Y");
		h_TDC_tof2B[ican+24]->SetLineColor(3);
		h_TDC_tof2B[ican+24]->SetFillStyle(3004);
		h_TDC_tof2B[ican+24]->SetFillColor(3);
		h_TDC_tof2B[ican+24]->SetStats(0);
		h_TDC_tof2B[ican+24]->Draw();

		h_TDC_tof2A[ican+24]->SetTitle(TOF2Title6);
		h_TDC_tof2A[ican+24]->SetAxisRange(0, TOF2_height_max[ican+24],"Y");
		h_TDC_tof2A[ican+24]->SetLineColor(2);
		h_TDC_tof2A[ican+24]->SetFillStyle(3004);
		h_TDC_tof2A[ican+24]->SetFillColor(2);
		h_TDC_tof2A[ican+24]->SetStats(0);
		h_TDC_tof2A[ican+24]->Draw("same");

		h_TDC_tof2C[ican+24]->SetTitle(TOF2Title6);
		h_TDC_tof2C[ican+24]->SetAxisRange(0, TOF2_height_max[ican+24],"Y");
		h_TDC_tof2C[ican+24]->SetLineColor(4);
		h_TDC_tof2C[ican+24]->SetFillStyle(3004);
		h_TDC_tof2C[ican+24]->SetFillColor(4);
		h_TDC_tof2C[ican+24]->SetStats(0);
		h_TDC_tof2C[ican+24]->Draw("same");

		h_TDC_tof2D[ican+24]->SetTitle(TOF2Title6);
		h_TDC_tof2D[ican+24]->SetAxisRange(0, TOF2_height_max[ican+24],"Y");
		h_TDC_tof2D[ican+24]->SetLineColor(1);
		h_TDC_tof2D[ican+24]->SetFillStyle(3004);
		h_TDC_tof2D[ican+24]->SetFillColor(1);
		h_TDC_tof2D[ican+24]->SetStats(0);
		h_TDC_tof2D[ican+24]->Draw("same");

		l2->Draw("same");
	}
	
	for(Int_t ican=0; ican<12; ican++){
		c8->cd(TOF2_Mapping3[ican]);

		char TOF2Title7[250];
		sprintf(TOF2Title7, "Raw TDC (Ch. %d) - TOF2 | OutB-%d : Run %d", ican+12, ican+1, Run_Number);
		h_TDC_tof2B[ican+12]->SetTitle(TOF2Title7);
		h_TDC_tof2B[ican+12]->SetAxisRange(0, TOF2_height_max[ican+12],"Y");
		h_TDC_tof2B[ican+12]->SetLineColor(3);
		h_TDC_tof2B[ican+12]->SetFillStyle(3004);
		h_TDC_tof2B[ican+12]->SetFillColor(3);
		h_TDC_tof2B[ican+12]->SetStats(0);
		h_TDC_tof2B[ican+12]->Draw();

		h_TDC_tof2A[ican+12]->SetTitle(TOF2Title7);
		h_TDC_tof2A[ican+12]->SetAxisRange(0, TOF2_height_max[ican+12],"Y");
		h_TDC_tof2A[ican+12]->SetLineColor(2);
		h_TDC_tof2A[ican+12]->SetFillStyle(3004);
		h_TDC_tof2A[ican+12]->SetFillColor(2);
		h_TDC_tof2A[ican+12]->SetStats(0);
		h_TDC_tof2A[ican+12]->Draw("same");

		h_TDC_tof2C[ican+12]->SetTitle(TOF2Title7);
		h_TDC_tof2C[ican+12]->SetAxisRange(0, TOF2_height_max[ican+12],"Y");
		h_TDC_tof2C[ican+12]->SetLineColor(4);
		h_TDC_tof2C[ican+12]->SetFillStyle(3004);
		h_TDC_tof2C[ican+12]->SetFillColor(4);
		h_TDC_tof2C[ican+12]->SetStats(0);
		h_TDC_tof2C[ican+12]->Draw("same");

		h_TDC_tof2D[ican+12]->SetTitle(TOF2Title7);
		h_TDC_tof2D[ican+12]->SetAxisRange(0, TOF2_height_max[ican+12],"Y");
		h_TDC_tof2D[ican+12]->SetLineColor(1);
		h_TDC_tof2D[ican+12]->SetFillStyle(3004);
		h_TDC_tof2D[ican+12]->SetFillColor(1);
		h_TDC_tof2D[ican+12]->SetStats(0);
		h_TDC_tof2D[ican+12]->Draw("same");

		l2->Draw("same");
	}

	for(Int_t ican=0; ican<12; ican++){
		c8->cd(TOF2_Mapping4[ican]);

		char TOF2Title8[250];
		sprintf(TOF2Title8, "Raw TDC (Ch. %d) - TOF2 | InB-%d : Run %d", ican+36, ican+1, Run_Number);
		h_TDC_tof2B[ican+36]->SetTitle(TOF2Title8);
		h_TDC_tof2B[ican+36]->SetAxisRange(0, TOF2_height_max[ican+36],"Y");
		h_TDC_tof2B[ican+36]->SetLineColor(3);
		h_TDC_tof2B[ican+36]->SetFillStyle(3004);
		h_TDC_tof2B[ican+36]->SetFillColor(3);
		h_TDC_tof2B[ican+36]->SetStats(0);
		h_TDC_tof2B[ican+36]->Draw();

		h_TDC_tof2A[ican+36]->SetTitle(TOF2Title8);
		h_TDC_tof2A[ican+36]->SetAxisRange(0, TOF2_height_max[ican+36],"Y");
		h_TDC_tof2A[ican+36]->SetLineColor(2);
		h_TDC_tof2A[ican+36]->SetFillStyle(3004);
		h_TDC_tof2A[ican+36]->SetFillColor(2);
		h_TDC_tof2A[ican+36]->SetStats(0);
		h_TDC_tof2A[ican+36]->Draw("same");

		h_TDC_tof2C[ican+36]->SetTitle(TOF2Title8);
		h_TDC_tof2C[ican+36]->SetAxisRange(0, TOF2_height_max[ican+36],"Y");
		h_TDC_tof2C[ican+36]->SetLineColor(4);
		h_TDC_tof2C[ican+36]->SetFillStyle(3004);
		h_TDC_tof2C[ican+36]->SetFillColor(4);
		h_TDC_tof2C[ican+36]->SetStats(0);
		h_TDC_tof2C[ican+36]->Draw("same");

		h_TDC_tof2D[ican+36]->SetTitle(TOF2Title8);
		h_TDC_tof2D[ican+36]->SetAxisRange(0, TOF2_height_max[ican+36],"Y");
		h_TDC_tof2D[ican+36]->SetLineColor(1);
		h_TDC_tof2D[ican+36]->SetFillStyle(3004);
		h_TDC_tof2D[ican+36]->SetFillColor(1);
		h_TDC_tof2D[ican+36]->SetStats(0);
		h_TDC_tof2D[ican+36]->Draw("same");

		l2->Draw("same");
	}


	c7->Print(pdf_in, "pdf");
	c8->Print(pdf_in, "pdf");

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
		
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//MWPC_Event_Display

	for (int i=0; i<16; i++) {
		h_C2YL->Fill(i,C2YL[i]);
		h_C2YR->Fill(i,C2YR[i]);
		h_C3YL->Fill(i,C3YL[i]);
		h_C3YR->Fill(i,C3YR[i]);
		h_C4YL->Fill(i,C4YL[i]);
		h_C4YR->Fill(i,C4YR[i]);

		h_C2XL[0]->Fill(i,C2X1L[i]);
		h_C2XR[0]->Fill(i,C2X1R[i]);
		h_C2XL[1]->Fill(i,C2X2L[i]);
		h_C2XR[1]->Fill(i,C2X2R[i]);
		h_C2XL[2]->Fill(i,C2X3L[i]);
		h_C2XR[2]->Fill(i,C2X3R[i]);
		h_C2XL[3]->Fill(i,C2X4L[i]);
		h_C2XR[3]->Fill(i,C2X4R[i]);

		h_C3XL[0]->Fill(i,C3X1L[i]);
		h_C3XR[0]->Fill(i,C3X1R[i]);
		h_C3XL[1]->Fill(i,C3X2L[i]);
		h_C3XR[1]->Fill(i,C3X2R[i]);
		h_C3XL[2]->Fill(i,C3X3L[i]);
		h_C3XR[2]->Fill(i,C3X3R[i]);
		h_C3XL[3]->Fill(i,C3X4L[i]);
		h_C3XR[3]->Fill(i,C3X4R[i]);

		h_C4XL[0]->Fill(i,C4X1L[i]);
		h_C4XR[0]->Fill(i,C4X1R[i]);
		h_C4XL[1]->Fill(i,C4X2L[i]);
		h_C4XR[1]->Fill(i,C4X2R[i]);
		h_C4XL[2]->Fill(i,C4X3L[i]);
		h_C4XR[2]->Fill(i,C4X3R[i]);
		h_C4XL[3]->Fill(i,C4X4L[i]);
		h_C4XR[3]->Fill(i,C4X4R[i]);
	}

	for(int i=0; i<14; i++){	
		h_C2XL_M->Fill(i,C2X1L[i]);
		h_C2XL_M->Fill(i+14,C2X2L[i]);
		h_C2XL_M->Fill(i+28,C2X3L[i]);
		h_C2XL_M->Fill(i+42,C2X4L[i]);

		h_C2XR_M->Fill(i,C2X1R[i]);
		h_C2XR_M->Fill(i+14,C2X2R[i]);
		h_C2XR_M->Fill(i+28,C2X3R[i]);
		h_C2XR_M->Fill(i+42,C2X4R[i]);
	}

	for(int i=0; i<16; i++){
		h_C3XL_M->Fill(i,C3X1L[i]);
		h_C3XL_M->Fill(i+16,C3X2L[i]);
		h_C3XL_M->Fill(i+32,C3X3L[i]);
		h_C3XL_M->Fill(i+48,C3X4L[i]);

		h_C3XR_M->Fill(i,C3X1R[i]);
		h_C3XR_M->Fill(i+16,C3X2R[i]);
		h_C3XR_M->Fill(i+32,C3X3R[i]);
		h_C3XR_M->Fill(i+48,C3X4R[i]);

		h_C4XL_M->Fill(i,C4X1L[i]);
		h_C4XL_M->Fill(i+16,C4X2L[i]);
		h_C4XL_M->Fill(i+32,C4X3L[i]);
		h_C4XL_M->Fill(i+48,C4X4L[i]);
		h_C4XL_M->Fill(i+64,C4X5L[i]);

		h_C4XR_M->Fill(i,C4X1R[i]);
		h_C4XR_M->Fill(i+16,C4X2R[i]);
		h_C4XR_M->Fill(i+32,C4X3R[i]);
		h_C4XR_M->Fill(i+48,C4X4R[i]);
		h_C4XR_M->Fill(i+64,C4X5R[i]);

		h_C4X5L->Fill(i,C4X5L[i]);
		h_C4X5R->Fill(i,C4X5R[i]);
	}

	for (int i=0; i<16; i++) {
		if (C2X1L[i] > max_C2XL) max_C2XL = C2X1L[i];
		if (C2X1R[i] > max_C2XR) max_C2XR = C2X1R[i];

		if (C3X1L[i] > max_C3XL) max_C3XL = C3X1L[i];
		if (C3X1R[i] > max_C3XR) max_C3XR = C3X1R[i];

		if (C4X1L[i] > max_C4XL) max_C4XL = C4X1L[i];
		if (C4X1R[i] > max_C4XR) max_C4XR = C4X1R[i];
	}

	for (int i=0; i<16; i++) {
		if (C2X2L[i] > max_C2XL) max_C2XL = C2X2L[i];
		if (C2X2R[i] > max_C2XR) max_C2XR = C2X2R[i];

		if (C3X2L[i] > max_C3XL) max_C3XL = C3X2L[i];
		if (C3X2R[i] > max_C3XR) max_C3XR = C3X2R[i];

		if (C4X2L[i] > max_C4XL) max_C4XL = C4X2L[i];
		if (C4X2R[i] > max_C4XR) max_C4XR = C4X2R[i];
	}


	for (int i=0; i<16; i++) {
		if (C2X3L[i] > max_C2XL) max_C2XL = C2X3L[i];
		if (C2X3R[i] > max_C2XR) max_C2XR = C2X3R[i];

		if (C3X3L[i] > max_C3XL) max_C3XL = C3X3L[i];
		if (C3X3R[i] > max_C3XR) max_C3XR = C3X3R[i];

		if (C4X3L[i] > max_C4XL) max_C4XL = C4X3L[i];
		if (C4X3R[i] > max_C4XR) max_C4XR = C4X3R[i];
	}


	for (int i=0; i<16; i++) {
		if (C2X4L[i] > max_C2XL) max_C2XL = C2X4L[i];
		if (C2X4R[i] > max_C2XR) max_C2XR = C2X4R[i];

		if (C3X4L[i] > max_C3XL) max_C3XL = C3X4L[i];
		if (C3X4R[i] > max_C3XR) max_C3XR = C3X4R[i];

		if (C4X4L[i] > max_C4XL) max_C4XL = C4X4L[i];
		if (C4X4R[i] > max_C4XR) max_C4XR = C4X4R[i];
	}

	double Gap_Height_C2XL = double(double(max_C2XL)*0.8);
	double Gap_Height_C2XR = double(double(max_C2XR)*0.8);
	double Gap_Height_C3XL = double(double(max_C3XL)*0.8);
	double Gap_Height_C3XR = double(double(max_C3XR)*0.8);
	double Gap_Height_C4XL = double(double(max_C4XL)*0.8);
	double Gap_Height_C4XR = double(double(max_C4XR)*0.8);

	//double avg_C4X5L = double(double(C4X5L_sum)/double(C4X5L_count));
	//double avg_C4X5LE = double(double(C4X5LE_sum)/double(C4X5LE_count));
	//double avg_C4X5R = double(double(C4X5R_sum)/double(C4X5R_count));
	//double avg_C4X5RE = double(double(C4X5RE_sum)/double(C4X5RE_count));


	char Name_Can_MWPC[100];			char Title_Can_MWPC[100];
	sprintf(Name_Can_MWPC,"MWPC Left Channels -- Run %d",Run_Number);
	sprintf(Title_Can_MWPC,"MWPC Left Channels -- Run %d",Run_Number);

	TCanvas* c1_MWPC = new TCanvas(Name_Can_MWPC, Title_Can_MWPC,1300,1000);
	c1_MWPC->Divide(6,3);
	c1_MWPC->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(int i=0; i<4; i++) {
		c1_MWPC->cd(i+1);
		h_C2XL[i]->SetStats(0);
		h_C2XL[i]->Draw();
	}

	c1_MWPC->cd(5);
	h_C2YL->SetStats(0);
	h_C2YL->Draw();

	for(int i=0; i<4; i++) {
		c1_MWPC->cd(i+7);
		h_C3XL[i]->SetStats(0);
		h_C3XL[i]->Draw();
	}

	c1_MWPC->cd(11);
	h_C3YL->SetStats(0);
	h_C3YL->Draw();;

	for(int i=0; i<4; i++) {
		c1_MWPC->cd(i+13);
		h_C4XL[i]->SetStats(0);
		h_C4XL[i]->Draw();
	}

	c1_MWPC->cd(17);
	h_C4X5L->SetStats(0);
	h_C4X5L->Draw();

	c1_MWPC->cd(18);	
	h_C4YL->SetStats(0);	
	h_C4YL->Draw();

	char Name_Can2[100];			char Title_Can2[100];
	sprintf(Name_Can2,"MWPC Right Channels -- Run %d",Run_Number);
	sprintf(Title_Can2,"MWPC Right Channels -- Run %d",Run_Number);

	TCanvas* c2_MWPC = new TCanvas(Name_Can2, Title_Can2,1300,1000);
	c2_MWPC->Divide(6,3);
	c2_MWPC->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(int i=0; i<4; i++) {
		c2_MWPC->cd(i+1);
		h_C2XR[i]->SetStats(0);
		h_C2XR[i]->Draw();
	}

	c2_MWPC->cd(5);
	h_C2YR->SetStats(0);
	h_C2YR->Draw();;

	for(int i=0; i<4; i++) {
		c2_MWPC->cd(i+7);
		h_C3XR[i]->SetStats(0);
		h_C3XR[i]->Draw();
	}

	c2_MWPC->cd(11);
	h_C3YR->SetStats(0);
	h_C3YR->Draw();

	for(int i=0; i<4; i++) {
		c2_MWPC->cd(i+13);
		h_C4XR[i]->SetStats(0);
		h_C4XR[i]->Draw();
	}

	c2_MWPC->cd(17);
	h_C4X5R->SetStats(0);
	h_C4X5R->Draw();

	c2_MWPC->cd(18);
	h_C4YR->SetStats(0);
	h_C4YR->Draw();

	char Name_Can3[100];			char Title_Can3[100];
	sprintf(Name_Can3,"MWPC Left & Right Channels -- Run %d",Run_Number);
	sprintf(Title_Can3,"MWPC Left & Right Channels -- Run %d",Run_Number);

	TCanvas* c3_MWPC = new TCanvas(Name_Can3, Title_Can3,1300,1000);
	c3_MWPC->Divide(4,3);
	c3_MWPC->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	TH1D *gaps = new TH1D("Gaps1", "Gaps1", 69,-1,68);
	gaps->Fill(16,Gap_Height_C2XL);
	gaps->Fill(33,Gap_Height_C2XL);
	gaps->Fill(50,Gap_Height_C2XL);
	gaps->Fill(67,Gap_Height_C2XL);
	gaps->SetLineColor(2);

	TH1D *gaps2 = new TH1D("Gaps2", "Gaps2", 69,-1,68);
	gaps2->Fill(16,Gap_Height_C2XR);
	gaps2->Fill(33,Gap_Height_C2XR);
	gaps2->Fill(50,Gap_Height_C2XR);
	gaps2->Fill(67,Gap_Height_C2XR);
	gaps2->SetLineColor(2);

	TH1D *gaps3 = new TH1D("Gaps3", "Gaps3", 69,-1,68);
	gaps3->Fill(16,Gap_Height_C3XL);
	gaps3->Fill(33,Gap_Height_C3XL);
	gaps3->Fill(50,Gap_Height_C3XL);
	gaps3->Fill(67,Gap_Height_C3XL);
	gaps3->SetLineColor(2);

	TH1D *gaps4 = new TH1D("Gaps4", "Gaps4", 69,-1,68);
	gaps4->Fill(16,Gap_Height_C3XR);
	gaps4->Fill(33,Gap_Height_C3XR);
	gaps4->Fill(50,Gap_Height_C3XR);
	gaps4->Fill(67,Gap_Height_C3XR);
	gaps4->SetLineColor(2);

	TH1D *gaps5 = new TH1D("Gaps5", "Gaps5", 86,-1,85);
	gaps5->Fill(16,Gap_Height_C4XL);
	gaps5->Fill(33,Gap_Height_C4XL);
	gaps5->Fill(50,Gap_Height_C4XL);
	gaps5->Fill(67,Gap_Height_C4XL);
	gaps5->Fill(84,Gap_Height_C4XL);
	gaps5->SetLineColor(2);

	TH1D *gaps6 = new TH1D("Gaps6", "Gaps6", 86,-1,85);
	gaps6->Fill(16,Gap_Height_C4XR);
	gaps6->Fill(33,Gap_Height_C4XR);
	gaps6->Fill(50,Gap_Height_C4XR);
	gaps6->Fill(67,Gap_Height_C4XR);
	gaps6->Fill(84,Gap_Height_C4XR);
	gaps6->SetLineColor(2);

	c3_MWPC->cd(1);
	h_C2XL_M->SetStats(0);
	h_C2XL_M->Draw();
	gaps->Draw("same");

	c3_MWPC->cd(2);
	h_C2XR_M->SetStats(0);
	h_C2XR_M->Draw();
	gaps2->Draw("same");

	c3_MWPC->cd(3);
	h_C2YL->SetStats(0);
	h_C2YL->Draw();

	c3_MWPC->cd(4);
	h_C2YR->SetStats(0);
	h_C2YR->Draw();

	c3_MWPC->cd(5);
	h_C3XL_M->SetStats(0);
	h_C3XL_M->Draw();
	gaps3->Draw("same");

	c3_MWPC->cd(6);
	h_C3XR_M->SetStats(0);
	h_C3XR_M->Draw();
	gaps4->Draw("same");

	c3_MWPC->cd(7);
	h_C3YL->SetStats(0);
	h_C3YL->Draw();

	c3_MWPC->cd(8);
	h_C3YR->SetStats(0);
	h_C3YR->Draw();

	c3_MWPC->cd(9);
	h_C4XL_M->SetStats(0);
	h_C4XL_M->Draw();
	gaps5->Draw("same");

	c3_MWPC->cd(10);
	h_C4XR_M->SetStats(0);
	h_C4XR_M->Draw();
	gaps6->Draw("same");

	c3_MWPC->cd(11);
	h_C4YL->SetStats(0);
	h_C4YL->Draw();

	c3_MWPC->cd(12);
	h_C4YR->SetStats(0);
	h_C4YR->Draw();

	c1_MWPC->Print(pdf_in, "pdf");
	c2_MWPC->Print(pdf_in, "pdf");
	c3_MWPC->Print(pdf_in, "pdf");

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Histo_Display


	char Can_TOF1_Name[100];			char Can_TOF1_Title[100];
	char Can_TOF2_RIGHT_Name[100];		char Can_TOF2_RIGHT_Title[100];
	char Can_TOF2_LEFT_Name[100];		char Can_TOF2_LEFT_Title[100];

	sprintf(Can_TOF1_Name, "TOF1 ADCs & TDCs  --  Run %d", Run_Number);
	sprintf(Can_TOF1_Title, "TOF1 ADCs & TDCs  --  Run %d", Run_Number);

	sprintf(Can_TOF2_RIGHT_Name, "TOF2 ADCs & TDCs (Gap 1 - 6)  --  Run %d", Run_Number);
	sprintf(Can_TOF2_RIGHT_Title, "TOF2 ADCs & TDCs (Gap 1 - 6)  --  Run %d", Run_Number);

	sprintf(Can_TOF2_LEFT_Name, "TOF2 ADCs & TDCs (Gap 7 - 12)  --  Run %d", Run_Number);
	sprintf(Can_TOF2_LEFT_Title, "TOF2 ADCs & TDCs (Gap 7 - 12)  --  Run %d", Run_Number);

	TCanvas *c1_TOF1_Raw;
	c1_TOF1_Raw = new TCanvas(Can_TOF1_Name,Can_TOF1_Title,1300,1000); 
	c1_TOF1_Raw->Divide(8,6);
	c1_TOF1_Raw->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	TCanvas *c1_TOF2_RIGHT_Raw;
	c1_TOF2_RIGHT_Raw = new TCanvas(Can_TOF2_RIGHT_Name,Can_TOF2_RIGHT_Title,1300,1000);
	c1_TOF2_RIGHT_Raw->Divide(8,6);
	c1_TOF2_RIGHT_Raw->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	TCanvas *c1_TOF2_LEFT_Raw;
	c1_TOF2_LEFT_Raw = new TCanvas(Can_TOF2_LEFT_Name,Can_TOF2_LEFT_Title,1300,1000);
	c1_TOF2_LEFT_Raw->Divide(8,6);
	c1_TOF2_LEFT_Raw->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");


	char TOF1_ADC_UP_Title[100];			char TOF1_TDC_UP_Title[100];
	char TOF1_ADC_DOWN_Title[100];			char TOF1_TDC_DOWN_Title[100];

	char TOF2_ADC_AO_RIGHT_Title[100];		char TOF2_TDC_AO_RIGHT_Title[100];
	char TOF2_ADC_AI_RIGHT_Title[100];		char TOF2_TDC_AI_RIGHT_Title[100];
	char TOF2_ADC_BO_RIGHT_Title[100];		char TOF2_TDC_BO_RIGHT_Title[100];
	char TOF2_ADC_BI_RIGHT_Title[100];		char TOF2_TDC_BI_RIGHT_Title[100];

	char TOF2_ADC_AO_LEFT_Title[100];		char TOF2_TDC_AO_LEFT_Title[100];
	char TOF2_ADC_AI_LEFT_Title[100];		char TOF2_TDC_AI_LEFT_Title[100];
	char TOF2_ADC_BO_LEFT_Title[100];		char TOF2_TDC_BO_LEFT_Title[100];
	char TOF2_ADC_BI_LEFT_Title[100];		char TOF2_TDC_BI_LEFT_Title[100];


	for(int i=0; i<12; i++){
		sprintf(TOF1_ADC_UP_Title, "Raw ADC (Ch. %d) - TOF1 | Up-%d : Run %d", i, i+1, Run_Number);
		sprintf(TOF1_ADC_DOWN_Title, "Raw ADC (Ch. %d) - TOF1 | Down-%d : Run %d", i+12, i+1, Run_Number);


		c1_TOF1_Raw->cd(2*i+1)->SetLogy();
		//h_ADC_tof1_Raw[i]->SetTitle(TOF1_ADC_UP_Title);
		h_ADC_tof1_Raw[i]->Draw();
		c1_TOF1_Raw->cd(2*i+2)->SetLogy();
		//h_ADC_tof1_Raw[i+12]->SetTitle(TOF1_ADC_DOWN_Title);
		h_ADC_tof1_Raw[i+12]->Draw();	

	}
		

	for(int j=0; j<6; j++){

		sprintf(TOF2_ADC_AO_RIGHT_Title, "Raw ADC (Ch. %d) - TOF2 | OutA-%d : Run %d", j, j+1, Run_Number);
		sprintf(TOF2_ADC_AO_LEFT_Title, "Raw ADC (Ch. %d) - TOF2 | OutA-%d : Run %d", j+6, j+7, Run_Number);

		sprintf(TOF2_ADC_AI_RIGHT_Title, "Raw ADC (Ch. %d) - TOF2 | InA-%d : Run %d", j+24, j+1, Run_Number);
		sprintf(TOF2_ADC_AI_LEFT_Title, "Raw ADC (Ch. %d) - TOF2 | InA-%d : Run %d", j+30, j+7, Run_Number);

		sprintf(TOF2_ADC_BO_RIGHT_Title, "Raw ADC (Ch. %d) - TOF2 | OutB-%d : Run %d", j+12, j+1, Run_Number);
		sprintf(TOF2_ADC_BO_LEFT_Title, "Raw ADC (Ch. %d) - TOF2 | OutB-%d : Run %d", j+18, j+7, Run_Number);

		sprintf(TOF2_ADC_BI_RIGHT_Title, "Raw ADC (Ch. %d) - TOF2 | InB-%d : Run %d", j+36, j+1, Run_Number);
		sprintf(TOF2_ADC_BI_LEFT_Title, "Raw ADC (Ch. %d) - TOF2 | InB-%d : Run %d", j+42, j+7, Run_Number);



		c1_TOF2_RIGHT_Raw->cd(4*j+1)->SetLogy();
		//h_ADC_tof2_Raw[j]->SetTitle(TOF2_ADC_AO_RIGHT_Title);
		h_ADC_tof2_Raw[j]->Draw();

		c1_TOF2_RIGHT_Raw->cd(4*j+2)->SetLogy();
		//h_ADC_tof2_Raw[j+12]->SetTitle(TOF2_ADC_AI_RIGHT_Title);
		h_ADC_tof2_Raw[j+12]->Draw();

		c1_TOF2_RIGHT_Raw->cd(4*j+3)->SetLogy();
		//h_ADC_tof2_Raw[j+24]->SetTitle(TOF2_ADC_BO_RIGHT_Title);
		h_ADC_tof2_Raw[j+24]->Draw();

		c1_TOF2_RIGHT_Raw->cd(4*j+4)->SetLogy();
		//h_ADC_tof2_Raw[j+36]->SetTitle(TOF2_ADC_BI_RIGHT_Title);
		h_ADC_tof2_Raw[j+36]->Draw();


		c1_TOF2_LEFT_Raw->cd(4*j+1)->SetLogy();
		//h_ADC_tof2_Raw[j+6]->SetTitle(TOF2_ADC_AO_LEFT_Title);
		h_ADC_tof2_Raw[j+6]->Draw();

		c1_TOF2_LEFT_Raw->cd(4*j+2)->SetLogy();
		//h_ADC_tof2_Raw[j+6+24]->SetTitle(TOF2_ADC_AI_LEFT_Title);
		h_ADC_tof2_Raw[j+6+12]->Draw();

		c1_TOF2_LEFT_Raw->cd(4*j+3)->SetLogy();
		//h_ADC_tof2_Raw[j+6+12]->SetTitle(TOF2_ADC_BO_LEFT_Title);
		h_ADC_tof2_Raw[j+6+24]->Draw();

		c1_TOF2_LEFT_Raw->cd(4*j+4)->SetLogy();
		//h_ADC_tof2_Raw[j+6+36]->SetTitle(TOF2_ADC_BI_LEFT_Title);
		h_ADC_tof2_Raw[j+6+36]->Draw();
	}

	//char TOF2_ADC_BO_LEFT_ch55_Title[100];
	//sprintf(TOF2_ADC_BO_LEFT_ch55_Title, "Raw ADC (Ch. 55) - TOF2 | OutB-7 : Run %d", Run_Number);

	//c1_TOF2_LEFT_Raw->cd(3)->SetLogy();
	//h_ADC_tof2_Raw[55]->SetTitle(TOF2_ADC_BO_LEFT_ch55_Title);
	//h_ADC_tof2_Raw[55]->Draw();


	for(int j=0; j<6; j++){
		sprintf(TOF2_TDC_AO_RIGHT_Title, "Raw TDC (Ch. %d) - TOF2 | OutA-%d : Run %d", j, j+1, Run_Number); 
		sprintf(TOF2_TDC_AO_LEFT_Title, "Raw TDC (Ch. %d) - TOF2 | OutA-%d : Run %d", j+6, j+7, Run_Number); 

		sprintf(TOF2_TDC_AI_RIGHT_Title, "Raw TDC (Ch. %d) - TOF2 | InA-%d : Run %d", j+24, j+1, Run_Number); 
		sprintf(TOF2_TDC_AI_LEFT_Title, "Raw TDC (Ch. %d) - TOF2 | InA-%d : Run %d", j+30, j+7, Run_Number); 

		sprintf(TOF2_TDC_BO_RIGHT_Title, "Raw TDC (Ch. %d) - TOF2 | OutB-%d : Run %d", j+12, j+1, Run_Number); 
		sprintf(TOF2_TDC_BO_LEFT_Title, "Raw TDC (Ch. %d) - TOF2 | OutB-%d : Run %d", j+18, j+7, Run_Number); 

		sprintf(TOF2_TDC_BI_RIGHT_Title, "Raw TDC (Ch. %d) - TOF2 | InB-%d : Run %d", j+36, j+1, Run_Number); 
		sprintf(TOF2_TDC_BI_LEFT_Title, "Raw TDC (Ch. %d) - TOF2 | InB-%d : Run %d", j+42, j+7, Run_Number); 



		c1_TOF2_RIGHT_Raw->cd(4*j+25);
		//h_TDC_tof2AO_Raw[j]->SetTitle(TOF2_TDC_AO_RIGHT_Title);
		h_TDC_tof2AO_Raw[j]->Draw();

		c1_TOF2_RIGHT_Raw->cd(4*j+26);
		//h_TDC_tof2AI_Raw[j]->SetTitle(TOF2_TDC_AI_RIGHT_Title);
		h_TDC_tof2AI_Raw[j]->Draw();

		c1_TOF2_RIGHT_Raw->cd(4*j+27);
		//h_TDC_tof2BO_Raw[j]->SetTitle(TOF2_TDC_BO_RIGHT_Title);
		h_TDC_tof2BO_Raw[j]->Draw();

		c1_TOF2_RIGHT_Raw->cd(4*j+28);
		//h_TDC_tof2BI_Raw[j]->SetTitle(TOF2_TDC_BI_RIGHT_Title);
		h_TDC_tof2BI_Raw[j]->Draw(); 




		c1_TOF2_LEFT_Raw->cd(4*j+25);
		//h_TDC_tof2AO_Raw[j+6]->SetTitle(TOF2_TDC_AO_LEFT_Title);
		h_TDC_tof2AO_Raw[j+6]->Draw();

		c1_TOF2_LEFT_Raw->cd(4*j+26);
		//h_TDC_tof2AI_Raw[j+6]->SetTitle(TOF2_TDC_AI_LEFT_Title);
		h_TDC_tof2AI_Raw[j+6]->Draw();

		c1_TOF2_LEFT_Raw->cd(4*j+27);
		//h_TDC_tof2BO_Raw[j+6]->SetTitle(TOF2_TDC_BO_LEFT_Title);
		h_TDC_tof2BO_Raw[j+6]->Draw();

		c1_TOF2_LEFT_Raw->cd(4*j+28);
		//h_TDC_tof2BI_Raw[j+6]->SetTitle(TOF2_TDC_BI_LEFT_Title);
		h_TDC_tof2BI_Raw[j+6]->Draw();

	}


	for(int j=0; j<12; j++){
		sprintf(TOF1_TDC_UP_Title, "Raw TDC (Ch. %d) - TOF1 | Up-%d : Run %d", j, j+1, Run_Number);
		sprintf(TOF1_TDC_DOWN_Title, "Raw TDC (Ch. %d) - TOF1 | Down-%d : Run %d", j+12, j+1, Run_Number);

		c1_TOF1_Raw->cd(2*j+1+24);
		//h_TDC_tof1U_Raw[j]->SetTitle(TOF1_TDC_UP_Title);	
		h_TDC_tof1U_Raw[j]->Draw();
		c1_TOF1_Raw->cd(2*j+2+24);
		//h_TDC_tof1D_Raw[j]->SetTitle(TOF1_TDC_DOWN_Title);	
	    h_TDC_tof1D_Raw[j]->Draw();
	}



		c1_TOF1_Raw->Print(pdf_in, "pdf");
		c1_TOF2_RIGHT_Raw->Print(pdf_in, "pdf");
		c1_TOF2_LEFT_Raw->Print(pdf_closed, "pdf");

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
} 
