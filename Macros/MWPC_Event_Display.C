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
#include "CommonParameters.h"
#include "C2_Strip_transform.h"
#endif
  
void MWPC_Event_Display(Int_t Run_Number=5, Int_t Event_Number=-1){ 

	gStyle->SetOptStat(0);

	//Int_t MWPC_thr = 300;

	Int_t adc_high_target[256];   		Int_t ADC_High_TARGET[256];    	//Double_t ADC_High_TARGET_corr[256];   	
	Int_t adc_low_target[256]; 			Int_t ADC_Low_TARGET[256]; 	//Double_t ADC_Low_TARGET_corr[256];
	Int_t tdc_le_target[256][16];      	Int_t TDC_LE_TARGET[256];    	//Double_t TDC_LE_TARGET_corr[256];   	
	Int_t tdc_te_target[256][16]; 		Int_t TDC_TE_TARGET[256]; 	//Double_t TDC_TE_TARGET_corr[256];

	Int_t adc_high_sft[128];   		Int_t ADC_High_SFT[128];    	//Double_t ADC_High_SFT_corr[128];   	
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

	Int_t C2XL56[56] = {0};
	Int_t C2XR56[56] = {0};

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

	TH1D *h_C2XR56;
	TH1D *h_C2XL56;
	TH1D *h_Strip_Left;		char Title_Strip_Left[100];		char Name_Strip_Left[100];
	TH1D *h_Strip_Right;	char Title_Strip_Right[100];	char Name_Strip_Right[100];
	TH1D *h_Z_Left;			char Title_Z_Left[100];			char Name_Z_Left[100];
	TH1D *h_Z_Right;		char Title_Z_Right[100];		char Name_Z_Right[100];




	for(int i=1; i<5; i++) {
		sprintf(Title_C2XL[i-1], "MWPC - C2X%dL -- Run %d | ADC > %d (Gaps 7-12)", i, Run_Number, MWPC_thr);
		sprintf(Name_C2XL[i-1], "MWPC - C2X%dL -- Run %d | ADC > %d (Gaps 7-12)", i, Run_Number, MWPC_thr);
		sprintf(Title_C2XR[i-1], "MWPC - C2X%dR -- Run %d | ADC > %d (Gaps 1-6)", i, Run_Number, MWPC_thr);
		sprintf(Name_C2XR[i-1], "MWPC - C2X%dR -- Run %d | ADC > %d (Gaps 1-6)", i, Run_Number, MWPC_thr);

		sprintf(Title_C3XL[i-1], "MWPC - C3X%dL -- Run %d | ADC > %d (Gaps 7-12)", i, Run_Number, MWPC_thr);
		sprintf(Name_C3XL[i-1], "MWPC - C3X%dL -- Run %d | ADC > %d (Gaps 7-12)", i, Run_Number, MWPC_thr);
		sprintf(Title_C3XR[i-1], "MWPC - C3X%dR -- Run %d | ADC > %d (Gaps 1-6)", i, Run_Number, MWPC_thr);
		sprintf(Name_C3XR[i-1], "MWPC - C3X%dR -- Run %d | ADC > %d (Gaps 1-6)", i, Run_Number, MWPC_thr);

		sprintf(Title_C4XL[i-1], "MWPC - C4X%dL -- Run %d | ADC > %d (Gaps 7-12)", i, Run_Number, MWPC_thr);
		sprintf(Name_C4XL[i-1], "MWPC - C4X%dL -- Run %d | ADC > %d (Gaps 7-12)", i, Run_Number, MWPC_thr);
		sprintf(Title_C4XR[i-1], "MWPC - C4X%dR -- Run %d | ADC > %d (Gaps 1-6)", i, Run_Number, MWPC_thr);
		sprintf(Name_C4XR[i-1], "MWPC - C4X%dR -- Run %d | ADC > %d (Gaps 1-6)", i, Run_Number, MWPC_thr);
	}

	sprintf(Title_C2YL, "MWPC - C2YL -- Run %d | ADC > %d (Gaps 7-12)", Run_Number, MWPC_thr);
	sprintf(Name_C2YL, "MWPC - C2YL -- Run %d | ADC > %d (Gaps 7-12)", Run_Number, MWPC_thr);
	sprintf(Title_C2YR, "MWPC - C2YR -- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);
	sprintf(Name_C2YR, "MWPC - C2YR-- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);

	sprintf(Title_C3YL, "MWPC - C3YL -- Run %d | ADC > %d (Gaps 7-12)", Run_Number, MWPC_thr);
	sprintf(Name_C3YL, "MWPC - C3YL -- Run %d | ADC > %d (Gaps 7-12)", Run_Number, MWPC_thr);
	sprintf(Title_C3YR, "MWPC - C3YR -- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);
	sprintf(Name_C3YR, "MWPC - C3YR-- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);

	sprintf(Title_C4YL, "MWPC - C4YL -- Run %d | ADC > %d (Gaps 7-12)", Run_Number, MWPC_thr);
	sprintf(Name_C4YL, "MWPC - C4YL -- Run %d | ADC > %d (Gaps 7-12)", Run_Number, MWPC_thr);
	sprintf(Title_C4YR, "MWPC - C4YR -- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);
	sprintf(Name_C4YR, "MWPC - C4YR-- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);

	sprintf(Title_C4X5L, "MWPC - C4X5L -- Run %d | ADC > %d (Gaps 7-12)" , Run_Number, MWPC_thr);
	sprintf(Name_C4X5L, "MWPC - C4X5L -- Run %d | ADC > %d (Gaps 7-12)", Run_Number, MWPC_thr);
	sprintf(Title_C4X5R, "MWPC - C4X5R -- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);
	sprintf(Name_C4X5R, "MWPC - C4X5R-- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);

	sprintf(Title_C2XL_M, "MWPC - C2XL -- Run %d | ADC > %d (Gaps 7-12)", Run_Number, MWPC_thr);
	sprintf(Name_C2XL_M, "MWPC - C2XL -- Run %d | ADC > %d (Gaps 7-12)", Run_Number, MWPC_thr);
	sprintf(Title_C2XR_M, "MWPC - C2XR -- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);
	sprintf(Name_C2XR_M, "MWPC - C2XR -- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);

	sprintf(Title_C3XL_M, "MWPC - C3XL -- Run %d | ADC > %d (Gaps 7-12)", Run_Number, MWPC_thr);
	sprintf(Name_C3XL_M, "MWPC - C3XL -- Run %d | ADC > %d (Gaps 7-12)", Run_Number, MWPC_thr);
	sprintf(Title_C3XR_M, "MWPC - C3XR -- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);
	sprintf(Name_C3XR_M, "MWPC - C3XR -- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);

	sprintf(Title_C4XL_M, "MWPC - C4XL -- Run %d | ADC > %d (Gaps 7-12)", Run_Number, MWPC_thr);
	sprintf(Name_C4XL_M, "MWPC - C4XL -- Run %d | ADC > %d (Gaps 7-12)", Run_Number, MWPC_thr);
	sprintf(Title_C4XR_M, "MWPC - C4XR -- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);
	sprintf(Name_C4XR_M, "MWPC - C4XR -- Run %d | ADC > %d (Gaps 1-6)", Run_Number, MWPC_thr);

	sprintf(Title_Strip_Left, "STRIPS LEFT");
	sprintf(Name_Strip_Left, "strips left");
	sprintf(Title_Strip_Right, "STRIPS RIGHT");
	sprintf(Name_Strip_Right, "strips right");

	sprintf(Title_Z_Left, "Z LEFT");
	sprintf(Name_Z_Left, "z left");
	sprintf(Title_Z_Right, "Z RIGHT");
	sprintf(Name_Z_Right, "z right");


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

	h_C2XL_M = new TH1D(Name_C2XL_M, Title_C2XL_M, 66,-1,65);
	h_C2XR_M = new TH1D(Name_C2XR_M, Title_C2XR_M, 66,-1,65);
	h_C3XL_M = new TH1D(Name_C3XL_M, Title_C3XL_M, 69,-1,68);
	h_C3XR_M = new TH1D(Name_C3XR_M, Title_C3XR_M, 69,-1,68);	
	h_C4XL_M = new TH1D(Name_C4XL_M, Title_C4XL_M, 86,-1,85);
	h_C4XR_M = new TH1D(Name_C4XR_M, Title_C4XR_M, 86,-1,85);

	h_C2XL56 = new TH1D(Name_C2XL_M, Title_C2XL_M, 66,-1,65);
	h_C2XR56 = new TH1D(Name_C2XR_M, Title_C2XR_M, 66,-1,65);

	h_Strip_Left = new TH1D(Name_Strip_Left, Title_Strip_Left, 58,-1,57);
	h_Strip_Right = new TH1D(Name_Strip_Right, Title_Strip_Right, 58,-1,57);

	h_Z_Left = new TH1D(Name_Z_Left, Title_Z_Left, 80,-40,40);
	h_Z_Right = new TH1D(Name_Z_Right, Title_Z_Right, 80,-40,40);

	char path_input[200];                   char file_mapping[200];
	sprintf(path_input,"%s",path_merged);          sprintf(file_mapping,"../Mapping");
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

	fChain->SetBranchAddress("MWPCADC",MwpcADC);

	//Int_t nentries;
	Int_t nentries = (Int_t)fChain->GetEntries();
	cout << "  " << endl;
	cout << "****  Number of events: " << nentries << "  **** " <<endl;
	cout << "  " << endl;


	cout << "   " << endl;

	Int_t loop_start = 0;
	Int_t loop_end = nentries;

	if (Event_Number >= 0) {
		loop_start = Event_Number;
		loop_end = Event_Number + 1;
	}

	//if (flag!=0) nentries=flag;
	for(Int_t i=loop_start; i<loop_end; i++){
		fChain->GetEntry(i);

		if(i%10000==1)	cout<<"**** "<<i<<" events done"<<endl;
	//	cout << " " << endl;
	//	cout << "********* Event " << i << endl;
  	
	  	for (int j_TARGET=0; j_TARGET<256; j_TARGET++){
	  		ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-par_temp_TARGET[1][j_TARGET];
			ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-par_temp_TARGET[1][j_TARGET];
			TDC_LE_TARGET[j_TARGET]=tdc_le_target[j_TARGET][0];
			TDC_TE_TARGET[j_TARGET]=tdc_te_target[j_TARGET][0];
  	
	//   	cout << j_TARGET << "   " << ADC_Low_TARGET[j_TARGET] << endl;
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
			//if(MWPCADC[j_MWPC]>=0) cout << "xxxxx  " << j_MWPC << "   " <<  MwpcADC[j_MWPC] << endl;
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
	}

	C2X1L[14]=0;	C2X1L[15]=0;	C2X1R[14]=0;	C2X1R[15]=0;
	C2X2L[14]=0;	C2X2L[15]=0;	C2X1R[14]=0;	C2X1R[15]=0;
	C2X3L[14]=0;	C2X3L[15]=0;	C2X3R[14]=0;	C2X3R[15]=0;
	C2X4L[14]=0;	C2X4L[15]=0;	C2X4R[14]=0;	C2X4R[15]=0;


	for (int i=0; i<16; i++) {
		h_C2YL->Fill(i+1,C2YL[i]);
		h_C2YR->Fill(i+1,C2YR[i]);
		h_C3YL->Fill(i+1,C3YL[i]);
		h_C3YR->Fill(i+1,C3YR[i]);
		h_C4YL->Fill(i+1,C4YL[i]);
		h_C4YR->Fill(i+1,C4YR[i]);

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


	for(int i=0; i<16; i++){
		h_C2XL_M->Fill(i,C2X1L[i]);
		h_C2XL_M->Fill(i+16,C2X2L[i]);
		h_C2XL_M->Fill(i+32,C2X3L[i]);
		h_C2XL_M->Fill(i+48,C2X4L[i]);

		h_C2XR_M->Fill(i,C2X1R[i]);
		h_C2XR_M->Fill(i+16,C2X2R[i]);
		h_C2XR_M->Fill(i+32,C2X3R[i]);
		h_C2XR_M->Fill(i+48,C2X4R[i]);
	}
	
	for(int i=0; i<16; i++){
		h_C3XL_M->Fill(i+1,C3X1L[i]);
		h_C3XL_M->Fill(i+17,C3X2L[i]);
		h_C3XL_M->Fill(i+33,C3X3L[i]);	
		h_C3XL_M->Fill(i+49,C3X4L[i]);

		h_C3XR_M->Fill(i+1,C3X1R[i]);
		h_C3XR_M->Fill(i+17,C3X2R[i]);
		h_C3XR_M->Fill(i+33,C3X3R[i]);
		h_C3XR_M->Fill(i+49,C3X4R[i]);

		h_C4XL_M->Fill(i+1,C4X1L[i]);
		h_C4XL_M->Fill(i+17,C4X2L[i]);
		h_C4XL_M->Fill(i+33,C4X3L[i]);
		h_C4XL_M->Fill(i+49,C4X4L[i]);
		h_C4XL_M->Fill(i+65,C4X5L[i]);

		h_C4XR_M->Fill(i+1,C4X1R[i]);
		h_C4XR_M->Fill(i+17,C4X2R[i]);
		h_C4XR_M->Fill(i+33,C4X3R[i]);
		h_C4XR_M->Fill(i+49,C4X4R[i]);
		h_C4XR_M->Fill(i+65,C4X5R[i]);

		h_C4X5L->Fill(i+1,C4X5L[i]);
		h_C4X5R->Fill(i+1,C4X5R[i]);
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

	cout << "" << endl;
	double avg_C4X5L = double(double(C4X5L_sum)/double(C4X5L_count));
	double avg_C4X5LE = double(double(C4X5LE_sum)/double(C4X5LE_count));
	double avg_C4X5R = double(double(C4X5R_sum)/double(C4X5R_count));
	double avg_C4X5RE = double(double(C4X5RE_sum)/double(C4X5RE_count));

	cout << "C4X5L Sum: " << C4X5L_sum << " -- Left Empty Sum: " << C4X5LE_sum << endl;
	cout << "C4X5R Sum: " << C4X5R_sum << " -- Right Empty Sum: " << C4X5RE_sum << endl;
	cout << "" << endl;

	cout << "C4X5L Count: " << C4X5L_count << " -- Left Empty Count: " << C4X5LE_count << endl;
	cout << "C4X5R Count: " << C4X5R_count << " -- Right Empty Count: " << C4X5RE_count << endl;
	cout << "" << endl;

	cout << "C4X5L Average: " << avg_C4X5L << " -- Left Empty Average: " << avg_C4X5LE << endl;
	cout << "C4X5R Average: " << avg_C4X5R << " -- Right Empty Average: " << avg_C4X5RE << endl;

	char Name_Can[100];			char Title_Can[100];
	sprintf(Name_Can,"MWPC Left Channels -- Run %d -- MWPC_Event_Display.C",Run_Number);
	sprintf(Title_Can,"MWPC Left Channels -- Run %d  -- MWPC_Event_Display.C",Run_Number);

	TCanvas* c1_MWPC = new TCanvas(Name_Can, Title_Can, 0,0,1300,600);
	c1_MWPC->Divide(6,3);
	c1_MWPC->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(int i=0; i<4; i++) {
		c1_MWPC->cd(i+1);
		h_C2XL[i]->Draw();
	}

	c1_MWPC->cd(5);
	h_C2YL->Draw();;

	for(int i=0; i<4; i++) {
		c1_MWPC->cd(i+7);
		h_C3XL[i]->Draw();
	}

	c1_MWPC->cd(11);
	h_C3YL->Draw();;

	for(int i=0; i<4; i++) {
		c1_MWPC->cd(i+13);
		h_C4XL[i]->Draw();
	}

	c1_MWPC->cd(17);
	h_C4X5L->Draw();

	c1_MWPC->cd(18);
	h_C4YL->Draw();

	char Name_Can2[100];			char Title_Can2[100];
	sprintf(Name_Can2,"MWPC Right Channels -- Run %d -- MWPC_Event_Display.C",Run_Number);
	sprintf(Title_Can2,"MWPC Right Channels -- Run %d  -- MWPC_Event_Display.C",Run_Number);

	TCanvas* c2_MWPC = new TCanvas(Name_Can2, Title_Can2, 150,300,1300,600);
	c2_MWPC->Divide(6,3);
	c2_MWPC->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(int i=0; i<4; i++) {
		c2_MWPC->cd(i+1);
		h_C2XR[i]->Draw();
	}

	c2_MWPC->cd(5);
	h_C2YR->Draw();;

	for(int i=0; i<4; i++) {
		c2_MWPC->cd(i+7);
		h_C3XR[i]->Draw();
	}

	c2_MWPC->cd(11);
	h_C3YR->Draw();

	for(int i=0; i<4; i++) {
		c2_MWPC->cd(i+13);
		h_C4XR[i]->Draw();
	}

	c2_MWPC->cd(17);
	h_C4X5R->Draw();

	c2_MWPC->cd(18);
	h_C4YR->Draw();

	char Name_Can3[100];			char Title_Can3[100];
	sprintf(Name_Can3,"MWPC Left & Right Channels -- Run %d -- MWPC_Event_Display.C",Run_Number);
	sprintf(Title_Can3,"MWPC Left & Right Channels -- Run %d  -- MWPC_Event_Display.C",Run_Number);

	TCanvas* c3_MWPC = new TCanvas(Name_Can3, Title_Can3, 700,900,1300,600);	
	c3_MWPC->Divide(4,3);
	c3_MWPC->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	TCanvas* c4_MWPC = new TCanvas("Test", "TEST", 700,900,600,900);	
	c4_MWPC->Divide(2,4);
	c4_MWPC->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");


	TH1D *gaps = new TH1D("Gaps1", "Gaps1", 69,-1,68);
	gaps->Fill(14,Gap_Height_C2XL);
	gaps->Fill(15,Gap_Height_C2XL);
	gaps->Fill(30,Gap_Height_C2XL);
	gaps->Fill(31,Gap_Height_C2XL);
	gaps->Fill(46,Gap_Height_C2XL);
	gaps->Fill(47,Gap_Height_C2XL);
	gaps->Fill(62,Gap_Height_C2XL);
	gaps->Fill(63,Gap_Height_C2XL);
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
	h_C2XL_M->Draw();
	gaps->Draw("same");

	c3_MWPC->cd(2);
	h_C2XR_M->Draw();
	gaps->Draw("same");

	c3_MWPC->cd(3);
	h_C2YL->Draw();

	c3_MWPC->cd(4);
	h_C2YR->Draw();

	c3_MWPC->cd(5);
	h_C3XL_M->Draw();
	//gaps3->Draw("same");

	c3_MWPC->cd(6);
	h_C3XR_M->Draw();
	//gaps4->Draw("same");

	c3_MWPC->cd(7);
	h_C3YL->Draw();

	c3_MWPC->cd(8);
	h_C3YR->Draw();

	c3_MWPC->cd(9);
	h_C4XL_M->Draw();
	//gaps5->Draw("same");

	c3_MWPC->cd(10);
	h_C4XR_M->Draw();
	//gaps6->Draw("same");

	c3_MWPC->cd(11);
	h_C4YL->Draw();

	c3_MWPC->cd(12);
	h_C4YR->Draw();




	c4_MWPC->cd(1);
	h_C2XL_M->Draw();
	gaps->Draw("same");

	c4_MWPC->cd(2);
	h_C2XR_M->Draw();
	gaps->Draw("same");

	// Transformation for C2

	Int_t C2XL_Strips[56];
	Int_t C2XR_Strips[56];





	for(int ii=0; ii<14; ii++){
		C2XL56[ii] = C2X1L[ii];
		C2XL56[ii+14] = C2X2L[ii];
		C2XL56[ii+28] = C2X3L[ii];
		C2XL56[ii+42] = C2X4L[ii];

		C2XR56[ii] = C2X1R[ii];
		C2XR56[ii+14] = C2X2R[ii];
		C2XR56[ii+28] = C2X3R[ii];
		C2XR56[ii+42] = C2X4R[ii];
	}

	for(int jj=0; jj<56; jj++){
		h_C2XL56->Fill(jj,C2XL56[jj]);
		h_C2XR56->Fill(jj,C2XR56[jj]);
	}

	for(int kk=0; kk<56; kk++){
		C2XL_Strips[kk] = C2XL56[C2_strip[kk]];
		C2XR_Strips[kk] = C2XR56[C2_strip[kk]];

		h_Strip_Left->Fill(kk, C2XL_Strips[kk]);
		h_Strip_Right->Fill(kk, C2XR_Strips[kk]);

		h_Z_Left->Fill(0.1*C2_ZLoc[kk], C2XL_Strips[kk]);
		h_Z_Right->Fill(0.1*C2_ZLoc[kk], C2XR_Strips[kk]);
	}



	c4_MWPC->cd(3);
	h_C2XL56->Draw();
	//gaps->Draw("same");

	c4_MWPC->cd(4);
	h_C2XR56->Draw();
	//gaps->Draw("same");

	c4_MWPC->cd(5);
	h_Strip_Left->Draw();
	//gaps->Draw("same");

	c4_MWPC->cd(6);
	h_Strip_Right->Draw();
	//gaps->Draw("same");

	c4_MWPC->cd(7);
	h_Z_Left->Draw();
	//gaps->Draw("same");

	c4_MWPC->cd(8);
	h_Z_Right->Draw();



	//for (int i=0; i<56; i++)  cout << i << "   " << C2_strip[i] << "   " <<  C2XL56[C2_strip[i]] << "   " << C2XL_Strips[i] << endl; 

}