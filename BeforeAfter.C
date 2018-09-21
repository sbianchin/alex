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

void BeforeaAfter(Int_t run_number_before=2853, Int_t run_number_after=3082, Int_t board_number=0) {   

	gStyle->Clear();
	TH1::AddDirectory(kFALSE);
	gStyle->SetOptStat(1111111);

	Int_t TDC_Graph_ymax = -1;
	Int_t TDC_Graph_xmax = 1500;

	TH1D *h_ADC_High_before[64];		TH1D *h_ADC_High_after[64];
	TH1D *h_ADC_Low_before[64];			TH1D *h_ADC_Low_after[64];
	TH1D *h_TDC_LE_before[64];			TH1D *h_TDC_LE_after[64];
	TH1D *h_TDC_TE_before[64];			TH1D *h_TDC_TE_after[64];

	char Name_ADC_High_before[100];		char Name_ADC_High_after[100];
	char Name_ADC_Low_before[100];		char Name_ADC_Low_after[100];
	char Name_TDC_LE_before[100];		char Name_TDC_LE_after[100];
	char Name_TDC_TE_before[100];		char Name_TDC_TE_after[100];

	char Name_Can_ADC_High[100];			char Title_Can_ADC_High[100];
	char Name_Can_ADC_Low[100];				char Title_Can_ADC_Low[100];
	char Name_Can_TDC_LE[100];				char Title_Can_TDC_LE[100];
	char Name_Can_TDC_TE[100];				char Title_Can_TDC_TE[100];

	char Name_Can_ADC_High2[100];			char Title_Can_ADC_High2[100];
	char Name_Can_ADC_Low2[100];			char Title_Can_ADC_Low2[100];
	char Name_Can_TDC_LE2[100];				char Title_Can_TDC_LE2[100];
	char Name_Can_TDC_TE2[100];				char Title_Can_TDC_TE2[100];

	char path_input_before[100];
	sprintf(path_input_before,"%s",path_histo);

	char path_input_after[100];
	sprintf(path_input_after,"%s",path_histo);

	char Name_finput_before[100];
	char Name_finput_after[100];

	sprintf(Name_finput_before,"%s/Hist_Run%dM.root",path_input_before, run_number_before);
	sprintf(Name_finput_after,"%s/Hist_Run%dM.root",path_input_after, run_number_after);

	cout << "     " << endl;
	cout << Name_finput_before << endl;
	cout << "     " << endl;
	cout << Name_finput_after << endl;
	cout << "     " << endl;

	TFile *finput_before = new TFile(Name_finput_before);

	if(board_number==1){
		for(int j=0;j<64;j++){
			sprintf(Name_ADC_High_before,"ADC_High (Ch. %d) - SFT",j);
			sprintf(Name_ADC_Low_before,"ADC_Low (Ch. %d) - SFT",j);
			sprintf(Name_TDC_LE_before,"TDC_LE (Ch. %d) - SFT",j);
			sprintf(Name_TDC_TE_before,"TDC_TE (Ch. %d) - SFT",j);
	
			h_ADC_High_before[j]=(TH1D*)finput_before->FindObjectAny(Name_ADC_High_before);
			h_ADC_Low_before[j]=(TH1D*)finput_before->FindObjectAny(Name_ADC_Low_before);
			h_TDC_LE_before[j]=(TH1D*)finput_before->FindObjectAny(Name_TDC_LE_before);
			h_TDC_TE_before[j]=(TH1D*)finput_before->FindObjectAny(Name_TDC_TE_before);
		}
	}

	if(board_number==2){
		for(int j=0;j<64;j++){
			sprintf(Name_ADC_High_before,"ADC_High (Ch. %d) - SFT",j+64);
			sprintf(Name_ADC_Low_before,"ADC_Low (Ch. %d) - SFT",j+64);
			sprintf(Name_TDC_LE_before,"TDC_LE (Ch. %d) - SFT",j+64);
			sprintf(Name_TDC_TE_before,"TDC_TE (Ch. %d) - SFT",j+64);
	
			h_ADC_High_before[j]=(TH1D*)finput_before->FindObjectAny(Name_ADC_High_before);
			h_ADC_Low_before[j]=(TH1D*)finput_before->FindObjectAny(Name_ADC_Low_before);
			h_TDC_LE_before[j]=(TH1D*)finput_before->FindObjectAny(Name_TDC_LE_before);
			h_TDC_TE_before[j]=(TH1D*)finput_before->FindObjectAny(Name_TDC_TE_before);
		}
	}

	if(board_number==3){
		for(int j=0;j<64;j++){
			sprintf(Name_ADC_High_before,"ADC_High (Ch. %d) - TARGET",j);
			sprintf(Name_ADC_Low_before,"ADC_Low (Ch. %d) - TARGET",j);
			sprintf(Name_TDC_LE_before,"TDC_LE (Ch. %d) - TARGET",j);
			sprintf(Name_TDC_TE_before,"TDC_TE (Ch. %d) - TARGET",j);
	
			h_ADC_High_before[j]=(TH1D*)finput_before->FindObjectAny(Name_ADC_High_before);
			h_ADC_Low_before[j]=(TH1D*)finput_before->FindObjectAny(Name_ADC_Low_before);
			h_TDC_LE_before[j]=(TH1D*)finput_before->FindObjectAny(Name_TDC_LE_before);
			h_TDC_TE_before[j]=(TH1D*)finput_before->FindObjectAny(Name_TDC_TE_before);
		}
	}

	if(board_number==4){
		for(int j=0;j<64;j++){
			sprintf(Name_ADC_High_before,"ADC_High (Ch. %d) - TARGET",j+64);
			sprintf(Name_ADC_Low_before,"ADC_Low (Ch. %d) - TARGET",j+64);
			sprintf(Name_TDC_LE_before,"TDC_LE (Ch. %d) - TARGET",j+64);
			sprintf(Name_TDC_TE_before,"TDC_TE (Ch. %d) - TARGET",j+64);

			h_ADC_High_before[j]=(TH1D*)finput_before->FindObjectAny(Name_ADC_High_before);
			h_ADC_Low_before[j]=(TH1D*)finput_before->FindObjectAny(Name_ADC_Low_before);
			h_TDC_LE_before[j]=(TH1D*)finput_before->FindObjectAny(Name_TDC_LE_before);
			h_TDC_TE_before[j]=(TH1D*)finput_before->FindObjectAny(Name_TDC_TE_before);
		}
	}

	if(board_number==5){
		for(int j=0;j<64;j++){
			sprintf(Name_ADC_High_before,"ADC_High (Ch. %d) - TARGET",j+128);
			sprintf(Name_ADC_Low_before,"ADC_Low (Ch. %d) - TARGET",j+128);
			sprintf(Name_TDC_LE_before,"TDC_LE (Ch. %d) - TARGET",j+128);
			sprintf(Name_TDC_TE_before,"TDC_TE (Ch. %d) - TARGET",j+128);
		
			h_ADC_High_before[j]=(TH1D*)finput_before->FindObjectAny(Name_ADC_High_before);
			h_ADC_Low_before[j]=(TH1D*)finput_before->FindObjectAny(Name_ADC_Low_before);
			h_TDC_LE_before[j]=(TH1D*)finput_before->FindObjectAny(Name_TDC_LE_before);
			h_TDC_TE_before[j]=(TH1D*)finput_before->FindObjectAny(Name_TDC_TE_before);
		}
	}

	if(board_number==6){
		for(int j=0;j<64;j++){
			sprintf(Name_ADC_High_before,"ADC_High (Ch. %d) - TARGET",j+192);
			sprintf(Name_ADC_Low_before,"ADC_Low (Ch. %d) - TARGET",j+192);
			sprintf(Name_TDC_LE_before,"TDC_LE (Ch. %d) - TARGET",j+192);
			sprintf(Name_TDC_TE_before,"TDC_TE (Ch. %d) - TARGET",j+192);
		
			h_ADC_High_before[j]=(TH1D*)finput_before->FindObjectAny(Name_ADC_High_before);
			h_ADC_Low_before[j]=(TH1D*)finput_before->FindObjectAny(Name_ADC_Low_before);
			h_TDC_LE_before[j]=(TH1D*)finput_before->FindObjectAny(Name_TDC_LE_before);
			h_TDC_TE_before[j]=(TH1D*)finput_before->FindObjectAny(Name_TDC_TE_before);
		}
	}



	TFile *finput_after = new TFile(Name_finput_after);

	if(board_number==1){
		for(int j=0;j<64;j++){
			sprintf(Name_ADC_High_after,"ADC_High (Ch. %d) - SFT",j);
			sprintf(Name_ADC_Low_after,"ADC_Low (Ch. %d) - SFT",j);
			sprintf(Name_TDC_LE_after,"TDC_LE (Ch. %d) - SFT",j);
			sprintf(Name_TDC_TE_after,"TDC_TE (Ch. %d) - SFT",j);
	
			h_ADC_High_after[j]=(TH1D*)finput_after->FindObjectAny(Name_ADC_High_after);
			h_ADC_Low_after[j]=(TH1D*)finput_after->FindObjectAny(Name_ADC_Low_after);
			h_TDC_LE_after[j]=(TH1D*)finput_after->FindObjectAny(Name_TDC_LE_after);
			h_TDC_TE_after[j]=(TH1D*)finput_after->FindObjectAny(Name_TDC_TE_after);
		}
	}

	if(board_number==2){
		for(int j=0;j<64;j++){
			sprintf(Name_ADC_High_after,"ADC_High (Ch. %d) - SFT",j+64);
			sprintf(Name_ADC_Low_after,"ADC_Low (Ch. %d) - SFT",j+64);
			sprintf(Name_TDC_LE_after,"TDC_LE (Ch. %d) - SFT",j+64);
			sprintf(Name_TDC_TE_after,"TDC_TE (Ch. %d) - SFT",j+64);
	
			h_ADC_High_after[j]=(TH1D*)finput_after->FindObjectAny(Name_ADC_High_after);
			h_ADC_Low_after[j]=(TH1D*)finput_after->FindObjectAny(Name_ADC_Low_after);
			h_TDC_LE_after[j]=(TH1D*)finput_after->FindObjectAny(Name_TDC_LE_after);
			h_TDC_TE_after[j]=(TH1D*)finput_after->FindObjectAny(Name_TDC_TE_after);
		}
	}

	if(board_number==3){
		for(int j=0;j<64;j++){
			sprintf(Name_ADC_High_after,"ADC_High (Ch. %d) - TARGET",j);
			sprintf(Name_ADC_Low_after,"ADC_Low (Ch. %d) - TARGET",j);
			sprintf(Name_TDC_LE_after,"TDC_LE (Ch. %d) - TARGET",j);
			sprintf(Name_TDC_TE_after,"TDC_TE (Ch. %d) - TARGET",j);
	
			h_ADC_High_after[j]=(TH1D*)finput_after->FindObjectAny(Name_ADC_High_after);
			h_ADC_Low_after[j]=(TH1D*)finput_after->FindObjectAny(Name_ADC_Low_after);
			h_TDC_LE_after[j]=(TH1D*)finput_after->FindObjectAny(Name_TDC_LE_after);
			h_TDC_TE_after[j]=(TH1D*)finput_after->FindObjectAny(Name_TDC_TE_after);
		}
	}

	if(board_number==4){
		for(int j=0;j<64;j++){
			sprintf(Name_ADC_High_after,"ADC_High (Ch. %d) - TARGET",j+64);
			sprintf(Name_ADC_Low_after,"ADC_Low (Ch. %d) - TARGET",j+64);
			sprintf(Name_TDC_LE_after,"TDC_LE (Ch. %d) - TARGET",j+64);
			sprintf(Name_TDC_TE_after,"TDC_TE (Ch. %d) - TARGET",j+64);

			h_ADC_High_after[j]=(TH1D*)finput_after->FindObjectAny(Name_ADC_High_after);
			h_ADC_Low_after[j]=(TH1D*)finput_after->FindObjectAny(Name_ADC_Low_after);
			h_TDC_LE_after[j]=(TH1D*)finput_after->FindObjectAny(Name_TDC_LE_after);
			h_TDC_TE_after[j]=(TH1D*)finput_after->FindObjectAny(Name_TDC_TE_after);
		}
	}

	if(board_number==5){
		for(int j=0;j<64;j++){
			sprintf(Name_ADC_High_after,"ADC_High (Ch. %d) - TARGET",j+128);
			sprintf(Name_ADC_Low_after,"ADC_Low (Ch. %d) - TARGET",j+128);
			sprintf(Name_TDC_LE_after,"TDC_LE (Ch. %d) - TARGET",j+128);
			sprintf(Name_TDC_TE_after,"TDC_TE (Ch. %d) - TARGET",j+128);
		
			h_ADC_High_after[j]=(TH1D*)finput_after->FindObjectAny(Name_ADC_High_after);
			h_ADC_Low_after[j]=(TH1D*)finput_after->FindObjectAny(Name_ADC_Low_after);
			h_TDC_LE_after[j]=(TH1D*)finput_after->FindObjectAny(Name_TDC_LE_after);
			h_TDC_TE_after[j]=(TH1D*)finput_after->FindObjectAny(Name_TDC_TE_after);
		}
	}

	if(board_number==6){
		for(int j=0;j<64;j++){
			sprintf(Name_ADC_High_after,"ADC_High (Ch. %d) - TARGET",j+192);
			sprintf(Name_ADC_Low_after,"ADC_Low (Ch. %d) - TARGET",j+192);
			sprintf(Name_TDC_LE_after,"TDC_LE (Ch. %d) - TARGET",j+192);
			sprintf(Name_TDC_TE_after,"TDC_TE (Ch. %d) - TARGET",j+192);
		
			h_ADC_High_after[j]=(TH1D*)finput_after->FindObjectAny(Name_ADC_High_after);
			h_ADC_Low_after[j]=(TH1D*)finput_after->FindObjectAny(Name_ADC_Low_after);
			h_TDC_LE_after[j]=(TH1D*)finput_after->FindObjectAny(Name_TDC_LE_after);
			h_TDC_TE_after[j]=(TH1D*)finput_after->FindObjectAny(Name_TDC_TE_after);
		}
	}



	if(board_number==1){
		sprintf(Name_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);

		sprintf(Name_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);

		sprintf(Name_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 0 - 31  (SFT)  | EASIROC  Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 32 - 63  (SFT)  | EASIROC  Board %d",run_number_before,run_number_after,board_number);

		sprintf(Name_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
	}

	if(board_number==2){
		sprintf(Name_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 64 - 95 (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 64 - 95 (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);

		sprintf(Name_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 64 - 95  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 64 - 95  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);

		sprintf(Name_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 64 - 95  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 64 - 95  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);

		sprintf(Name_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 64 - 95  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 64 - 95  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
	}

	if(board_number==3){
		sprintf(Name_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 0 - 32  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 0 - 32  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);

		sprintf(Name_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);

		sprintf(Name_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);

		sprintf(Name_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
	}

	if(board_number==4){
		sprintf(Name_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);

		sprintf(Name_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);

		sprintf(Name_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);

		sprintf(Name_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
	}


	if(board_number==5){
		sprintf(Name_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);

		sprintf(Name_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);

		sprintf(Name_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);

		sprintf(Name_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
	}

	if(board_number==6){
		sprintf(Name_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);

		sprintf(Name_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d - %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
					
		sprintf(Name_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);

		sprintf(Name_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Name_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
		sprintf(Title_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d - %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number_before,run_number_after,board_number);
	}


		TCanvas *c2;
		c2 = new TCanvas(Name_Can_ADC_High2,Title_Can_ADC_High2,1200,500); 
		//c2 = new TCanvas("jnvjnsd","feo",1200,500); 
		c2->Divide(8,8);
		c2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<32; ican++){
			c2->cd(ican+1)->SetLogy();  
			h_ADC_High_before[ican+32]->Draw();
			h_ADC_High_after[ican+32]->SetLineColor(2);
			h_ADC_High_after[ican+32]->Draw("same");
			c2->cd(ican+1+32);
			h_TDC_LE_before[ican+32]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			//if(TDC_Graph_ymax > 0) {  
			//	h_TDC_LE_before[ican+32]->SetAxisRange(0, TDC_Graph_ymax,"Y");
			//}	
			h_TDC_LE_before[ican+32]->Draw();
			h_TDC_LE_after[ican+32]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			h_TDC_LE_after[ican+32]->SetLineColor(2);
			h_TDC_LE_after[ican+32]->Draw("same");

		}


		TCanvas *c1;
		c1 = new TCanvas(Name_Can_ADC_Low,Title_Can_ADC_Low,1200,500); 
		c1->Divide(8,8);
		c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<32; ican++){
			c1->cd(ican+1)->SetLogy();  
			h_ADC_Low_before[ican]->Draw();
			h_ADC_Low_after[ican]->SetLineColor(2);
			h_ADC_Low_after[ican]->Draw("same");
			c1->cd(ican+1+32);
			h_TDC_LE_before[ican]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			h_TDC_LE_before[ican]->Draw();
			h_TDC_LE_after[ican]->SetAxisRange(0, TDC_Graph_xmax,"X");
			h_TDC_LE_after[ican]->SetLineColor(2);
			h_TDC_LE_after[ican]->Draw("same");
		}






}
