#ifndef __CINT__
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include "ANAPATH.h"
#endif

void Histo_Display(Int_t run_number=1, Int_t board_number=0, Int_t iflag=0) {   

	// if iflag=0 ->  ADC_High and TDC_LE ONLY !!!
	// if iflag=1 ->  Both ADCs and both TDCs !

	// board_number = 0  ->  ALL !

	// board_number = 1  ->  SFT (ch. 0 - 63)
	// board_number = 2  ->  SFT (ch. 64 - 127)

	// board_number = 3  ->  TARGET (ch. 0 - 63)
	// board_number = 4  ->  TARGET (ch. 64 - 127)
	// board_number = 5  ->  TARGET (ch. 128 - 191)
	// board_number = 6  ->  TARGET (ch. 192 - 255)

	// board_number = 7 -> TOF1 (Gap 1 - 12)  12 x 8
	// board_number = 77 -> TOF1 (Gap 1 - 12)  8 x 6

	// board_number = 8 -> TOF2 (Gap 1 - 12)

	// board_number = 9 -> MWPC (ch. 0 - 512)

	// board_number = 10 -> AC 
	// board_number = 11 -> GV 
	// board_number = 12 -> TTC 

	// board_number = 13 -> C_k & C_pi (ONLY TDC[0])
	// board_number = 14 -> TDC_Ck Average


	gStyle->Clear();
	TH1::AddDirectory(kFALSE);
	gStyle->SetOptStat(1111111);

	if (!((board_number == 1) || (board_number == 2) || (board_number == 3) || (board_number == 4) || 
		 (board_number == 5) || (board_number == 6) || (board_number == 7) || (board_number == 8) || 
		 (board_number == 9) || (board_number == 10) || (board_number == 11) || (board_number == 12) ||
		 (board_number == 13) || (board_number == 14) || (board_number == 77))) {
		cout << "  " << endl;
		cout << "Please enter a board number between 1 and 12" << endl;
		cout << "  " << endl;
		return;
	}


	Int_t TDC_Graph_ymax = -1;
	Int_t TDC_Graph_xmax = 1500;

	//Int_t TOF1_Graph_ymax = -1;	
	//Int_t TOF1_Graph_xmax = 3000;

	//Int_t TOF2_Graph_ymax = -1;
	//Int_t TOF2_Graph_xmax = 3000;

	TH1F *h_ADC_High[64];
	TH1F *h_ADC_Low[64];
	TH1F *h_TDC_LE[64];
	//TH1F *h_TDC_TE[64];

	//TH1F *h_ADC_tof1[24];
	//TH1F *h_ADC_tof2[56];

	TH1D *h_ADC_tof1U[12];		TH1D *h_ADC_tof1D[12];	
	TH1D *h_TDC_tof1U[12];		TH1D *h_TDC_tof1D[12];

														
	TH1F *h_ADC_tof2AO[12];		TH1F *h_TDC_tof2AO[12];
	TH1F *h_ADC_tof2AI[12];		TH1F *h_TDC_tof2AI[12];
	TH1F *h_ADC_tof2BO[12];		TH1F *h_TDC_tof2BO[12];
	TH1F *h_ADC_tof2BI[12];		TH1F *h_TDC_tof2BI[12];

	TH1F *h_ADC_MWPC[512];	

	TH1D *h_ADC_acu[12];		TH1D *h_TDC_acu[12];
	TH1D *h_ADC_acd[12];		TH1D *h_TDC_acd[12];

	TH1D *h_ADC_gv[12];			TH1D *h_TDC_gv[12];

	TH1D *h_ADC_ttc[12];		TH1D *h_TDC_ttc[12];

	TH1D *h_TDC_Ck[14][2];		TH1D *h_TDC_Cpi[14][2];		

	TH1D *h_TDC_Ck_Avg;	


	char Name_ADC_High[100];
	char Name_ADC_Low[100];
	char Name_TDC_LE[100];
	//char Name_TDC_TE[100];
	char Name_ADC_MWPC[100];
	//char Name_ADC_tof2[100];
	
	char Name_ADC_tof1U[100];	char Name_TDC_tof1U[100];
	char Name_ADC_tof1D[100];	char Name_TDC_tof1D[100];

	char Name_ADC_tof2AO[100];	char Name_TDC_tof2AO[100];
	char Name_ADC_tof2AI[100];	char Name_TDC_tof2AI[100];
	char Name_ADC_tof2BO[100];	char Name_TDC_tof2BO[100];
	char Name_ADC_tof2BI[100];	char Name_TDC_tof2BI[100];

	char Name_ADC_acu[100];		char Name_TDC_acu[100];
	char Name_ADC_acd[100];		char Name_TDC_acd[100];

	char Name_ADC_gv[100];		char Name_TDC_gv[100];

	char Name_ADC_ttc[100];		char Name_TDC_ttc[100];

	char Name_TDC_Ck[100];		char Name_TDC_Cpi[100];

	char Name_TDC_Ck_Avg[100];

	char Name_Can_ADC_High[100];			char Title_Can_ADC_High[100];
	char Name_Can_ADC_Low[100];				char Title_Can_ADC_Low[100];
	char Name_Can_TDC_LE[100];				char Title_Can_TDC_LE[100];
	char Name_Can_TDC_TE[100];				char Title_Can_TDC_TE[100];
	//char Name_Can_ADC_tof1[100];			char Title_Can_ADC_tof1[100];
	//char Name_Can_ADC_tof1B[100];			char Title_Can_ADC_tof1B[100];
	//char Name_Can_ADC_tof2[100];			char Title_Can_ADC_tof2[100];
	//char Name_Can_ADC_tof2B[100];			char Title_Can_ADC_tof2B[100];
	//char Name_Can_ADC_tof2C[100];			char Title_Can_ADC_tof2C[100];	
	//char Name_Can_ADC_tof2D[100];			char Title_Can_ADC_tof2D[100];
	char Name_Can_ADC_MWPC[100];			char Title_Can_ADC_MWPC[100];
	char Name_Can_ADC_MWPC2[100];			char Title_Can_ADC_MWPC2[100];

	char Name_Can_TOF1[100];				char Title_Can_TOF1[100];
	char Name_Can_TOF1_8_6[100];			char Title_Can_TOF1_8_6[100];

	char Name_Can_TOF2_1[100];				char Title_Can_TOF2_1[100];
	char Name_Can_TOF2_2[100];				char Title_Can_TOF2_2[100];

	char Name_Can_AC[100];					char Title_Can_AC[100];

	char Name_Can_GV[100];					char Title_Can_GV[100];

	char Name_Can_TTC[100];					char Title_Can_TTC[100];

	char Name_Can_Ck[100];					char Title_Can_Ck[100];
	char Name_Can_Cpi[100];					char Title_Can_Cpi[100];


	char Name_Can_ADC_High2[100];			char Title_Can_ADC_High2[100];
	char Name_Can_ADC_Low2[100];			char Title_Can_ADC_Low2[100];
	char Name_Can_TDC_LE2[100];				char Title_Can_TDC_LE2[100];
	char Name_Can_TDC_TE2[100];				char Title_Can_TDC_TE2[100];

	char Name_Can_TDC_Ck_Avg[100];			char Title_Can_TDC_Ck_Avg[100];

	Int_t TDC_Max_Can1=0;	Int_t max_Can1=0;
	Int_t TDC_Max_Can2=0;	Int_t max_Can2=0;

	//	Int_t TOF1_Mapping1[12] = {1,3,5,7,9,11,13,15,17,19,21,23};
	//	Int_t TOF1_Mapping2[12] = {2,4,6,8,10,12,14,16,18,20,22,24};

	//	Int_t TOF2_Mapping1[12] = {1,5,9,13,17,21,25,29,33,37,41,45};
	//	Int_t TOF2_Mapping2[12] = {2,6,10,14,18,22,26,30,34,38,42,46};
	//	Int_t TOF2_Mapping3[12] = {3,7,11,15,19,23,27,31,35,39,43,47};
	//	Int_t TOF2_Mapping4[12] = {4,8,12,16,20,24,28,32,36,40,44,48};

	//char path_input[100];
	//sprintf(path_input,"/triumfcs/trshare/trek/E36/Data/April_2015/Histograms");

	char path_input[100];
	sprintf(path_input,"%s",path_histo);

	char Name_finput[100];
	//sprintf(path_input,"/media/bianchin/hdd1/trek/E36/Data/April_2015/Histograms");
	sprintf(Name_finput,"%s/Hist_Run%dM.root",path_input, run_number);
	cout << Name_finput << endl;
	cout << "     " << endl;


	TFile *finput = new TFile(Name_finput);

	if(board_number==1){
		for(int j=0;j<64;j++){
			sprintf(Name_ADC_High,"ADC_High (Ch. %d) - SFT",j);
			sprintf(Name_ADC_Low,"ADC_Low (Ch. %d) - SFT",j);
			sprintf(Name_TDC_LE,"TDC_LE (Ch. %d) - SFT",j);
			//sprintf(Name_TDC_TE,"TDC_TE (Ch. %d) - SFT",j);
	
			h_ADC_High[j]=(TH1F*)finput->FindObjectAny(Name_ADC_High);
			h_ADC_Low[j]=(TH1F*)finput->FindObjectAny(Name_ADC_Low);
			h_TDC_LE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_LE);
			//h_TDC_TE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_TE);
		}
	}

	if(board_number==2){
		for(int j=0;j<64;j++){
			sprintf(Name_ADC_High,"ADC_High (Ch. %d) - SFT",j+64);
			sprintf(Name_ADC_Low,"ADC_Low (Ch. %d) - SFT",j+64);
			sprintf(Name_TDC_LE,"TDC_LE (Ch. %d) - SFT",j+64);
			//sprintf(Name_TDC_TE,"TDC_TE (Ch. %d) - SFT",j+64);
	
			h_ADC_High[j]=(TH1F*)finput->FindObjectAny(Name_ADC_High);
			h_ADC_Low[j]=(TH1F*)finput->FindObjectAny(Name_ADC_Low);
			h_TDC_LE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_LE);
			//h_TDC_TE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_TE);
		}
	}

	if(board_number==3){
		for(int j=0;j<64;j++){
			sprintf(Name_ADC_High,"ADC_High (Ch. %d) - TARGET",j);
			sprintf(Name_ADC_Low,"ADC_Low (Ch. %d) - TARGET",j);
			sprintf(Name_TDC_LE,"TDC_LE (Ch. %d) - TARGET",j);
			//sprintf(Name_TDC_TE,"TDC_TE (Ch. %d) - TARGET",j);
	
			h_ADC_High[j]=(TH1F*)finput->FindObjectAny(Name_ADC_High);
			h_ADC_Low[j]=(TH1F*)finput->FindObjectAny(Name_ADC_Low);
			h_TDC_LE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_LE);
			//h_TDC_TE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_TE);
		}
	}

	if(board_number==4){
		for(int j=0;j<64;j++){
			sprintf(Name_ADC_High,"ADC_High (Ch. %d) - TARGET",j+64);
			sprintf(Name_ADC_Low,"ADC_Low (Ch. %d) - TARGET",j+64);
			sprintf(Name_TDC_LE,"TDC_LE (Ch. %d) - TARGET",j+64);
			//sprintf(Name_TDC_TE,"TDC_TE (Ch. %d) - TARGET",j+64);

			h_ADC_High[j]=(TH1F*)finput->FindObjectAny(Name_ADC_High);
			h_ADC_Low[j]=(TH1F*)finput->FindObjectAny(Name_ADC_Low);
			h_TDC_LE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_LE);
			//h_TDC_TE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_TE);
		}
	}

	if(board_number==5){
		for(int j=0;j<64;j++){
			sprintf(Name_ADC_High,"ADC_High (Ch. %d) - TARGET",j+128);
			sprintf(Name_ADC_Low,"ADC_Low (Ch. %d) - TARGET",j+128);
			sprintf(Name_TDC_LE,"TDC_LE (Ch. %d) - TARGET",j+128);
			//sprintf(Name_TDC_TE,"TDC_TE (Ch. %d) - TARGET",j+128);
		
			h_ADC_High[j]=(TH1F*)finput->FindObjectAny(Name_ADC_High);
			h_ADC_Low[j]=(TH1F*)finput->FindObjectAny(Name_ADC_Low);
			h_TDC_LE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_LE);
			//h_TDC_TE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_TE);
		}
	}

	if(board_number==6){
		for(int j=0;j<64;j++){
			sprintf(Name_ADC_High,"ADC_High (Ch. %d) - TARGET",j+192);
			sprintf(Name_ADC_Low,"ADC_Low (Ch. %d) - TARGET",j+192);
			sprintf(Name_TDC_LE,"TDC_LE (Ch. %d) - TARGET",j+192);
			//sprintf(Name_TDC_TE,"TDC_TE (Ch. %d) - TARGET",j+192);
		
			h_ADC_High[j]=(TH1F*)finput->FindObjectAny(Name_ADC_High);
			h_ADC_Low[j]=(TH1F*)finput->FindObjectAny(Name_ADC_Low);
			h_TDC_LE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_LE);
			//h_TDC_TE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_TE);
		}
	}


	if(board_number==7 || board_number==77){
		for(int j=0; j<12; j++){
			sprintf(Name_ADC_tof1U,"ADC_TOF1 (%d UP)",j+1);
			sprintf(Name_ADC_tof1D,"ADC_TOF1 (%d DOWN)",j+1);
			sprintf(Name_TDC_tof1U,"TDC_TOF1 (%d UP)",j+1);
			sprintf(Name_TDC_tof1D,"TDC_TOF1 (%d DOWN)",j+1);

			h_ADC_tof1U[j]=(TH1D*)finput->FindObjectAny(Name_ADC_tof1U);
			h_ADC_tof1D[j]=(TH1D*)finput->FindObjectAny(Name_ADC_tof1D);
			h_TDC_tof1U[j]=(TH1D*)finput->FindObjectAny(Name_TDC_tof1U);
			h_TDC_tof1D[j]=(TH1D*)finput->FindObjectAny(Name_TDC_tof1D);
		}
	}


	if(board_number==8){
		for(int j=0; j<12; j++){
			sprintf(Name_ADC_tof2AO,"ADC_TOF2 (%d AO)",j+1);	
			sprintf(Name_ADC_tof2AI,"ADC_TOF2 (%d AI)",j+1);	
			sprintf(Name_ADC_tof2BO,"ADC_TOF2 (%d BO)",j+1);	
			sprintf(Name_ADC_tof2BI,"ADC_TOF2 (%d BI)",j+1);	

			sprintf(Name_TDC_tof2AO,"TDC_TOF2 (%d AO)",j+1);	
			sprintf(Name_TDC_tof2AI,"TDC_TOF2 (%d AI)",j+1);	
			sprintf(Name_TDC_tof2BO,"TDC_TOF2 (%d BO)",j+1);	
			sprintf(Name_TDC_tof2BI,"TDC_TOF2 (%d BI)",j+1);	

			h_ADC_tof2AO[j]=(TH1F*)finput->FindObjectAny(Name_ADC_tof2AO);
			h_ADC_tof2AI[j]=(TH1F*)finput->FindObjectAny(Name_ADC_tof2AI);
			h_ADC_tof2BO[j]=(TH1F*)finput->FindObjectAny(Name_ADC_tof2BO);
			h_ADC_tof2BI[j]=(TH1F*)finput->FindObjectAny(Name_ADC_tof2BI);

			h_TDC_tof2AO[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof2AO);
			h_TDC_tof2AI[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof2AI);
			h_TDC_tof2BO[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof2BO);
			h_TDC_tof2BI[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof2BI);
		}
	}


	if(board_number==9){
		for(int j=0; j<512;j++){
			sprintf(Name_ADC_MWPC,"ADC_MWPC (Ch. %d)",j);	

			h_ADC_MWPC[j]=(TH1F*)finput->FindObjectAny(Name_ADC_MWPC);
		}
	}


	if(board_number==10){
		for(int j=0; j<12; j++){
			sprintf(Name_ADC_acu, "ADC_AC (%d UP)", j+1);
			sprintf(Name_ADC_acd, "ADC_AC (%d DOWN)", j+1);
			sprintf(Name_TDC_acu, "TDC_AC (%d UP)", j+1);
			sprintf(Name_TDC_acd, "TDC_AC (%d DOWN)", j+1);

			h_ADC_acu[j]=(TH1D*)finput->FindObjectAny(Name_ADC_acu);
			h_ADC_acd[j]=(TH1D*)finput->FindObjectAny(Name_ADC_acd);
			h_TDC_acu[j]=(TH1D*)finput->FindObjectAny(Name_TDC_acu);
			h_TDC_acd[j]=(TH1D*)finput->FindObjectAny(Name_TDC_acd);
		}
	}


	if(board_number==11){
		for(int j=0; j<12; j++){
			sprintf(Name_ADC_gv, "ADC_GV (Counter %d)", j+1);
			sprintf(Name_TDC_gv, "TDC_GV (Counter %d)", j+1);

			h_ADC_gv[j]=(TH1D*)finput->FindObjectAny(Name_ADC_gv);
			h_TDC_gv[j]=(TH1D*)finput->FindObjectAny(Name_TDC_gv);
		}
	}


	if(board_number==12){
		for(int j=0; j<12; j++){
			sprintf(Name_ADC_ttc, "ADC_TTC (Counter %d)", j+1);
			sprintf(Name_TDC_ttc, "TDC_TTC (Counter %d)", j+1);

			h_ADC_ttc[j]=(TH1D*)finput->FindObjectAny(Name_ADC_ttc);
			h_TDC_ttc[j]=(TH1D*)finput->FindObjectAny(Name_TDC_ttc);
		}
	}

	if(board_number==13){
		for(int j=0; j<14; j++){
			sprintf(Name_TDC_Ck, "TDC_Ck (Counter %d)", j+1);
			sprintf(Name_TDC_Cpi, "TDC_Cpi (Counter %d)", j+1);

			h_TDC_Ck[j][0]=(TH1D*)finput->FindObjectAny(Name_TDC_Ck);
			h_TDC_Cpi[j][0]=(TH1D*)finput->FindObjectAny(Name_TDC_Cpi);
		}
	}

	if(board_number==14){
		sprintf(Name_TDC_Ck_Avg, "TDC Ck Average");
		h_TDC_Ck_Avg=(TH1D*)finput->FindObjectAny(Name_TDC_Ck_Avg);

	}

	//////////////////////////////////////////////////////

	if(board_number==1){
		sprintf(Name_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number,board_number);

		sprintf(Name_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number,board_number);

		sprintf(Name_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  | EASIROC  Board %d",run_number,board_number);
		sprintf(Name_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  | EASIROC  Board %d",run_number,board_number);

		sprintf(Name_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number,board_number);
	}

	if(board_number==2){
		sprintf(Name_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95 (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95 (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number,board_number);

		sprintf(Name_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 64 - 95  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 64 - 95  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number,board_number);

		sprintf(Name_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number,board_number);

		sprintf(Name_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number,board_number);
	}

	if(board_number==3){
		sprintf(Name_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 32  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 32  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number,board_number);

		sprintf(Name_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number,board_number);

		sprintf(Name_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number,board_number);

		sprintf(Name_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number,board_number);
	}

	if(board_number==4){
		sprintf(Name_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number,board_number);

		sprintf(Name_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number,board_number);

		sprintf(Name_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number,board_number);

		sprintf(Name_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number,board_number);
	}


	if(board_number==5){
		sprintf(Name_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number,board_number);

		sprintf(Name_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number,board_number);

		sprintf(Name_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number,board_number);

		sprintf(Name_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number,board_number);
	}

	if(board_number==6){
		sprintf(Name_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number,board_number);

		sprintf(Name_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number,board_number);
					
		sprintf(Name_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number,board_number);

		sprintf(Name_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Name_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number,board_number);
		sprintf(Title_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number,board_number);
	}

	if(board_number==7){
		sprintf(Name_Can_TOF1,"TOF1 ADCs & TDCs -- Run %d  (Gap 1 - 12)",run_number);
		sprintf(Title_Can_TOF1,"TOF1 ADCs & TDCs -- Run %d  (Gap 1 - 12)  |  1U, 2U, 3U, etc...",run_number);
	}

	if(board_number==8){
		sprintf(Name_Can_TOF2_1,"TOF2 ADCs & TDCs -- Run %d  (Gap 1 - 6)",run_number);
		sprintf(Title_Can_TOF2_1,"TOF2 ADCs & TDCs -- Run %d  (Gap 1 -6)  |  1AO, 2AO, 3AO, etc...",run_number);
		sprintf(Name_Can_TOF2_2,"TOF2 ADCs & TDCs -- Run %d  (Gap 7 - 12)",run_number);
		sprintf(Title_Can_TOF2_2,"TOF2 ADCs & TDCs -- Run %d  (Gap 7 - 12)  |  7AO, 8AO, 9AO, etc...",run_number);
	}

	if(board_number==9){
		sprintf(Name_Can_ADC_MWPC,"MWPC ADCs -- Run %d  (Ch. 0 - 255)",run_number);
		sprintf(Title_Can_ADC_MWPC,"MWPC ADCs -- Run %d  (Ch. 0 - 255)",run_number);

		sprintf(Name_Can_ADC_MWPC2,"MWPC ADCs -- Run %d  (Ch. 256 - 511)",run_number);
		sprintf(Title_Can_ADC_MWPC2,"MWPC ADCs -- Run %d  (Ch. 256 - 511)",run_number);
	}

	if(board_number==10){
		sprintf(Name_Can_AC,"Aerogel Counter -- Run %d  (Gap 1 - 12)",run_number);
		sprintf(Title_Can_AC,"Aerogel Counter -- Run %d  (Gap 1 - 12)",run_number);
	}

	if(board_number==11){
		sprintf(Name_Can_GV,"Gap Veto -- Run %d  (Gap 1 - 12)",run_number);
		sprintf(Title_Can_GV,"Gap Veto -- Run %d  (Gap 1 - 12)",run_number);
	}

	if(board_number==12){
		sprintf(Name_Can_TTC,"TTCs -- Run %d  (Gap 1 - 12)",run_number);
		sprintf(Title_Can_TTC,"TTCs -- Run %d  (Gap 1 - 12)",run_number);
	}

	if(board_number==13){
	//	sprintf(Name_Can_Ck,"Ck -- Run %d  (Counter 1 - 14)",run_number);
	//	sprintf(Title_Can_Ck,"Ck -- Run %d  (Counter 1 - 14)",run_number);
		sprintf(Name_Can_Ck,"Ck & Cpi -- Run %d  (Counter 1 - 14)  |  TDC[0]",run_number);
		sprintf(Title_Can_Ck,"Ck & Cpi -- Run %d  (Counter 1 - 14)  |  TDC[0]",run_number);

		sprintf(Name_Can_Cpi,"Cpi -- Run %d  (Counter 1 - 14)",run_number);
		sprintf(Title_Can_Cpi,"Cpi -- Run %d  (Counter 1 - 14)",run_number);
	}

	if(board_number==14){
		sprintf(Name_Can_TDC_Ck_Avg,"TDC Ck Average -- Run %d",run_number);
		sprintf(Title_Can_TDC_Ck_Avg,"TDC Ck Average");
	}

	if(board_number==77){
		sprintf(Name_Can_TOF1_8_6,"TOF1 ADCs & TDCs -- Run %d  (Gap 1 - 12)",run_number);
		sprintf(Title_Can_TOF1_8_6,"TOF1 ADCs & TDCs -- Run %d  (Gap 1 - 12)  |  1U, 1D, 2U, 2D, etc...",run_number);
	}




	if (board_number <= 6) {
		for(Int_t ii=0; ii<32; ii++){
			max_Can2=h_TDC_LE[ii]->GetMaximum();
			max_Can1=h_TDC_LE[ii+32]->GetMaximum();
			if(max_Can2 > TDC_Max_Can2){
				TDC_Max_Can2 = max_Can2;
			}
		
			if(max_Can1 > TDC_Max_Can1){
				TDC_Max_Can1 = max_Can1;
			}
		}
	}




	if ((board_number <= 6) && (iflag==0)) {

		TCanvas *c4;
		c4 = new TCanvas(Name_Can_ADC_High2,Title_Can_ADC_High2,1200,500); 
		c4->Divide(8,8);
		c4->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<32; ican++){
			//if(board_number==2 && (ican==0 || ican==1 || ican==2 || ican==3 ||ican==4 || ican==5 || ican==6 ||  ican==7 || ican==20)) c4->cd(ican+1)->SetFillColor(45);
			//if(board_number==2 && (ican==0 || ican==1 || ican==2 || ican==3 ||ican==4 || ican==5 || ican==6 ||  ican==7 || ican==20)) c4->cd(ican+1+32)->SetFillColor(45);
			//if(run_number>=1680 && board_number==1 && (ican>=0 && ican<=7)) c4->cd(ican+1)->SetFillColor(kYellow-9);
			//if(run_number>=1680 && board_number==1 && (ican>=0 && ican<=7)) c4->cd(ican+1+32)->SetFillColor(kYellow-9);
			//else c4->cd(ican+1);
			c4->cd(ican+1)->SetLogy();  
			h_ADC_High[ican+32]->Draw();
			c4->cd(ican+1+32);
			h_TDC_LE[ican+32]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			if(TDC_Graph_ymax > 0) {  
				h_TDC_LE[ican+32]->SetAxisRange(0, TDC_Graph_ymax,"Y");
			}	
			h_TDC_LE[ican+32]->Draw();
		}
	}


	if ((board_number <= 6) && (iflag==0)){

		TCanvas *c2;
		c2 = new TCanvas(Name_Can_TDC_LE,Title_Can_TDC_LE,0,200,1200,500); 
		c2->Divide(8,8);
		c2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<32; ican++){
			//if(board_number==1 && (ican==23 || ican==24 || ican==25 || ican==26 ||ican==27 || ican==28 || ican==29)) c2->cd(ican+1)->SetFillColor(45);
			//if(board_number==1 && (ican==23 || ican==24 || ican==25 || ican==26 ||ican==27 || ican==28 || ican==29)) c2->cd(ican+1+32)->SetFillColor(45);
			//if(run_number>=1680 && board_number==1 && (ican>=32 && ican<=39)) c2->cd(ican+1)->SetFillColor(kYellow-9);
			//if(run_number>=1680 && board_number==2 && (ican>=8 && ican<=15)) c2->cd(ican+1)->SetFillColor(kYellow-9);
			//if(run_number>=1680 && board_number==2 && (ican>=8 && ican<=15)) c2->cd(ican+1+32)->SetFillColor(kYellow-9);
			//if(run_number>=1680 && board_number==4 && (ican>=24 && ican<=31)) c2->cd(ican+1)->SetFillColor(kYellow-9);
			//if(run_number>=1680 && board_number==4 && (ican>=24 && ican<=31)) c2->cd(ican+1+32)->SetFillColor(kYellow-9);
			//if(run_number>=1680 && board_number==5 && (ican>=8 && ican<=15)) c2->cd(ican+1)->SetFillColor(kYellow-9);
			//if(run_number>=1680 && board_number==5 && (ican>=8 && ican<=15)) c2->cd(ican+1+32)->SetFillColor(kYellow-9);
			//if(run_number>=1680 && board_number==6 && (ican>=16 && ican<=23)) c2->cd(ican+1)->SetFillColor(kYellow-9);
			//if(run_number>=1680 && board_number==6 && (ican>=16 && ican<=23)) c2->cd(ican+1+32)->SetFillColor(kYellow-9);
			//else c2->cd(ican+1);
			c2->cd(ican+1)->SetLogy();
			h_ADC_High[ican]->Draw();
			c2->cd(ican+1+32);
			h_TDC_LE[ican]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			if(TDC_Graph_ymax > 0) {  
				h_TDC_LE[ican]->SetAxisRange(0, TDC_Graph_ymax,"Y");
			}
			h_TDC_LE[ican]->Draw();
		}
	}

	if((board_number <= 6) && (iflag==1)){

		TCanvas *c1;
		c1 = new TCanvas(Name_Can_ADC_Low,Title_Can_ADC_Low,1200,500); 
		c1->Divide(8,8);
		c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<32; ican++){
			c1->cd(ican+1)->SetLogy();  
			h_ADC_Low[ican]->Draw();
			c1->cd(ican+1+32);
			h_TDC_LE[ican]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			h_TDC_LE[ican]->Draw();
		}

		TCanvas *c3;
		c3 = new TCanvas(Name_Can_ADC_Low2,Title_Can_ADC_Low2,1200,500); 
		c3->Divide(8,8);
		c3->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<32; ican++){
			//if(board_number==1 && (ican==23 || ican==24 || ican==25 || ican==26 ||ican==27 || ican==28 || ican==29)) c3->cd(ican+1)->SetFillColor(45);
			//if(board_number==2 && (ican==32 || ican==33 || ican==34 || ican==35 ||ican==36 || ican==37 || ican==38 ||  ican==39 || ican==52)) c3->cd(ican+1)->SetFillColor(45);
			//if(run_number>=1680 && board_number==1 && (ican>=32 && ican<=39)) c3->cd(ican+1)->SetFillColor(kYellow-9);
			//if(run_number>=1680 && board_number==2 && (ican>=8 && ican<=15)) c3->cd(ican+1)->SetFillColor(kYellow-9);
			//if(run_number>=1680 && board_number==4 && (ican>=24 && ican<=31)) c3->cd(ican+1)->SetFillColor(kYellow-9);
			//if(run_number>=1680 && board_number==5 && (ican>=8 && ican<=15)) c3->cd(ican+1)->SetFillColor(kYellow-9);
			//if(run_number>=1680 && board_number==6 && (ican>=16 && ican<=23)) c3->cd(ican+1)->SetFillColor(kYellow-9);
			//else c3->cd(ican+1);
			c3->cd(ican+1)->SetLogy();  
			h_ADC_Low[ican+32]->Draw();
			c3->cd(ican+1+32);
			h_TDC_LE[ican+32]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			h_TDC_LE[ican+32]->Draw();
		}
	}


	if (board_number == 7) {

		TCanvas *c8;
		c8 = new TCanvas(Name_Can_TOF1,Title_Can_TOF1,1200,500); 
		c8->Divide(12,4);
		c8->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<12; ican++){
			c8->cd(ican+1)->SetLogy();
			h_ADC_tof1U[ican]->Draw();
		}

		for(Int_t ican=0; ican<12; ican++){
			c8->cd(ican+1+12);
			h_TDC_tof1U[ican]->Draw();
		}

		for(Int_t ican=0; ican<12; ican++){
			c8->cd(ican+1+24)->SetLogy();
			h_ADC_tof1D[ican]->Draw();
		}

		for(Int_t ican=0; ican<12; ican++){
			c8->cd(ican+1+36);
			h_TDC_tof1D[ican]->Draw();
		}


	}


	if (board_number == 8) {
		// Gap 7 - 12
		TCanvas *c10_2;
		c10_2 = new TCanvas(Name_Can_TOF2_2,Title_Can_TOF2_2,1200,700); 
		c10_2->Divide(6,8);
		c10_2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<6; ican++){
			c10_2->cd(ican+1)->SetLogy();
			h_ADC_tof2AO[ican+6]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c10_2->cd(ican+1+6);
			h_TDC_tof2AO[ican+6]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c10_2->cd(ican+1+12)->SetLogy();
			h_ADC_tof2AI[ican+6]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c10_2->cd(ican+1+18);
			h_TDC_tof2AI[ican+6]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c10_2->cd(ican+1+24)->SetLogy();
			h_ADC_tof2BO[ican+6]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c10_2->cd(ican+1+30);
			h_TDC_tof2BO[ican+6]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c10_2->cd(ican+1+36)->SetLogy();
			h_ADC_tof2BI[ican+6]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c10_2->cd(ican+1+42);
			h_TDC_tof2BI[ican+6]->Draw();
		}


		// Gap 1 - 6
		TCanvas *c10_1;
		c10_1 = new TCanvas(Name_Can_TOF2_1,Title_Can_TOF2_1,1200,700); 
		c10_1->Divide(6,8);
		c10_1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");


		for(Int_t ican=0; ican<6; ican++){
			c10_1->cd(ican+1)->SetLogy();
			h_ADC_tof2AO[ican]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c10_1->cd(ican+1+6);
			h_TDC_tof2AO[ican]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c10_1->cd(ican+1+12)->SetLogy();
			h_ADC_tof2AI[ican]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c10_1->cd(ican+1+18);
			h_TDC_tof2AI[ican]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c10_1->cd(ican+1+24)->SetLogy();
			h_ADC_tof2BO[ican]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c10_1->cd(ican+1+30);
			h_TDC_tof2BO[ican]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c10_1->cd(ican+1+36)->SetLogy();
			h_ADC_tof2BI[ican]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c10_1->cd(ican+1+42);
			h_TDC_tof2BI[ican]->Draw();
		}
	}


	if (board_number == 9) {

		TCanvas *c9;
		c9 = new TCanvas(Name_Can_ADC_MWPC,Title_Can_ADC_MWPC,1400,1000); 
		c9->Divide(16,16);
		c9->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<256; ican++){
			c9->cd(ican+1)->SetLogy();
			h_ADC_MWPC[ican]->Draw();
		}

		TCanvas *c10;
		c10 = new TCanvas(Name_Can_ADC_MWPC2,Title_Can_ADC_MWPC2,0,200,1400,1000); 
		c10->Divide(16,16);
		c10->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<256; ican++){
			c10->cd(ican+1)->SetLogy();
			h_ADC_MWPC[ican+256]->Draw();
		}
	}


	if (board_number == 10) {
	
		TCanvas *c12;
		c12 = new TCanvas(Name_Can_AC,Title_Can_AC,1300,450); 
		c12->Divide(12,4);
		c12->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<12; ican++){
			c12->cd(ican+1)->SetLogy();
			h_ADC_acu[ican]->Draw();
		}

		for(Int_t ican=0; ican<12; ican++){
			c12->cd(ican+1+12);
			h_TDC_acu[ican]->Draw();
		}

		for(Int_t ican=0; ican<12; ican++){
			c12->cd(ican+1+24)->SetLogy();
			h_ADC_acd[ican]->Draw();
		}

		for(Int_t ican=0; ican<12; ican++){
			c12->cd(ican+1+36);
			h_TDC_acd[ican]->Draw();
		}
	}

	
	if (board_number == 11) {
	
		TCanvas *c13;
		c13 = new TCanvas(Name_Can_GV,Title_Can_GV,1000,800); 
		c13->Divide(6,4);
		c13->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<6; ican++){
			c13->cd(ican+1)->SetLogy();
			h_ADC_gv[ican]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c13->cd(ican+1+6);
			h_TDC_gv[ican]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c13->cd(ican+1+12)->SetLogy();
			h_ADC_gv[ican+6]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c13->cd(ican+1+18);
			h_TDC_gv[ican+6]->Draw();
		}
	}

	
	if (board_number == 12) {
	
		TCanvas *c14;
		c14 = new TCanvas(Name_Can_TTC,Title_Can_TTC,1000,800); 
		c14->Divide(6,4);
		c14->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<6; ican++){
			c14->cd(ican+1)->SetLogy();
			h_ADC_ttc[ican]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c14->cd(ican+1+6);
			h_TDC_ttc[ican]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c14->cd(ican+1+12)->SetLogy();
			h_ADC_ttc[ican+6]->Draw();
		}

		for(Int_t ican=0; ican<6; ican++){
			c14->cd(ican+1+18);
			h_TDC_ttc[ican+6]->Draw();
		}
	}

	if (board_number == 13) {

		TCanvas *c15;
		c15 = new TCanvas(Name_Can_Ck,Title_Can_Ck,1200,800); 
//		c15->Divide(4,4);
		c15->Divide(7,4);
		c15->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		//TCanvas *c16;
		//c16 = new TCanvas(Name_Can_Cpi,Title_Can_Cpi,1000,800); 
		//c16->Divide(4,4);
		//c16->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<14; ican++){
			c15->cd(ican+1);
			h_TDC_Ck[ican][0]->Draw();
		}

		for(Int_t ican=0; ican<14; ican++){
			c15->cd(ican+1+14);
			h_TDC_Cpi[ican][0]->Draw();
		}
	}

	if (board_number==14){
		TCanvas *c16;
		c16 = new TCanvas(Name_Can_TDC_Ck_Avg, Title_Can_TDC_Ck_Avg, 800, 600);
		c16->cd();
		h_TDC_Ck_Avg->Draw();
	}


	if (board_number == 77) {

		TCanvas *c77;
		c77 = new TCanvas(Name_Can_TOF1_8_6,Title_Can_TOF1_8_6,1200,500); 
		c77->Divide(8,6);
		c77->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<12; ican++){
			c77->cd(2*ican+1)->SetLogy();
			h_ADC_tof1U[ican]->Draw();
			c77->cd(2*ican+2)->SetLogy();
			h_ADC_tof1D[ican]->Draw();
		}

		for(Int_t ican=0; ican<12; ican++){
			c77->cd(2*ican+1+24);
			h_TDC_tof1U[ican]->Draw();
			c77->cd(2*ican+2+24);
			h_TDC_tof1D[ican]->Draw();
		}
	}


	return;

}

