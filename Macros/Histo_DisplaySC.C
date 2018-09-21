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

void Histo_DisplaySC(Int_t run_number=1, Int_t board_number=0, Int_t iflag=0) {   

// if iflag=0 ->  ADC_High and TDC_LE ONLY !!!
// if iflag=1 ->  Both ADCs and both TDCs !

// board number = 1  ->  SFT (ch. 0 - 63)
// board number = 2  ->  SFT (ch. 64 - 127)

// board number = 3  ->  TARGET (ch. 0 - 63)
// board number = 4  ->  TARGET (ch. 64 - 127)
// board number = 5  ->  TARGET (ch. 128 - 191)
// board number = 6  ->  TARGET (ch. 192 - 255)

// board number = 7 -> TOF1 (ch. 0 - 23)
// board number = 8 -> TOF1 (ch. 0 - 23) - Sorted
// board number = 9 -> TOF2 (ch. 0 - 47)
// board number = 10 -> TOF2 (ch. 0 - 47) - Sorted
// board number = 11 -> MWPC (ch. 0 -512)

gStyle->SetOptStat(1111111);

for (int i = 0; i<1; i++) {
	if (!((board_number == 1) || (board_number == 2) || (board_number == 3) || (board_number == 4) || (board_number == 5) || (board_number == 6) || (board_number == 7) || (board_number == 8) || (board_number == 9) || (board_number == 10) || (board_number == 11))) {
		cout << "Please enter a board number between 1 and 11" << endl;
		break;
	}
}

Int_t TDC_Graph_ymax = -1;
Int_t TDC_Graph_xmax = 1500;

Int_t TOF1_Graph_ymax = -1;
Int_t TOF1_Graph_xmax = 3000;

Int_t TOF2_Graph_ymax = -1;
Int_t TOF2_Graph_xmax = 3000;

TH1F *h_ADC_High[64];
TH1F *h_ADC_Low[64];
TH1F *h_TDC_LE[64];
TH1F *h_TDC_TE[64];

TH1F *h_ADC_tof1[24];
TH1F *h_ADC_tof2[56];

TH1F *h_TDC_tof1U[12];   	char Name_TDC_tof1U[100];
TH1F *h_TDC_tof1D[12];   	char Name_TDC_tof1D[100];
TH1F *h_TDC_tof2AO[12];  	char Name_TDC_tof2AO[100];
TH1F *h_TDC_tof2AI[12];  	char Name_TDC_tof2AI[100];
TH1F *h_TDC_tof2BO[12];  	char Name_TDC_tof2BO[100];
TH1F *h_TDC_tof2BI[12];  	char Name_TDC_tof2BI[100];

TH1F *h_ADC_MWPC[512];

char Name_ADC_High[100];
char Name_ADC_Low[100];
char Name_TDC_LE[100];
char Name_TDC_TE[100];
char Name_ADC_tof1[100];
char Name_ADC_tof2[100];
char Name_ADC_MWPC[100];

char Name_Can_ADC_High[100];			char Title_Can_ADC_High[100];
char Name_Can_ADC_Low[100];				char Title_Can_ADC_Low[100];
char Name_Can_TDC_LE[100];				char Title_Can_TDC_LE[100];
char Name_Can_TDC_TE[100];			    char Title_Can_TDC_TE[100];
char Name_Can_ADC_tof1[100];			char Title_Can_ADC_tof1[100];
char Name_Can_ADC_tof1B[100];			char Title_Can_ADC_tof1B[100];
char Name_Can_ADC_tof2[100];			char Title_Can_ADC_tof2[100];
char Name_Can_ADC_tof2B[100];			char Title_Can_ADC_tof2B[100];
char Name_Can_ADC_tof2C[100];			char Title_Can_ADC_tof2C[100];
char Name_Can_ADC_tof2D[100];			char Title_Can_ADC_tof2D[100];
char Name_Can_ADC_MWPC[100];			char Title_Can_ADC_MWPC[100];
char Name_Can_ADC_MWPC2[100];			char Title_Can_ADC_MWPC2[100];

//char Name_Can_ADC_High1[100];			char Title_Can_ADC_High1[100];
//char Name_Can_ADC_Low1[100];			char Title_Can_ADC_Low1[100];
//char Name_Can_TDC_LE1[100];				char Title_Can_TDC_LE1[100];
//char Name_Can_TDC_TE1[100];				char Title_Can_TDC_TE1[100];

char Name_Can_ADC_High2[100];			char Title_Can_ADC_High2[100];
char Name_Can_ADC_Low2[100];			char Title_Can_ADC_Low2[100];
char Name_Can_TDC_LE2[100];				char Title_Can_TDC_LE2[100];
char Name_Can_TDC_TE2[100];				char Title_Can_TDC_TE2[100];

Int_t TDC_Max_Can1=0;	Int_t max_Can1=0;
Int_t TDC_Max_Can2=0;	Int_t max_Can2=0;

Int_t TOF1_Mapping1[12] = {1,3,5,7,9,11,13,15,17,19,21,23};
Int_t TOF1_Mapping2[12] = {2,4,6,8,10,12,14,16,18,20,22,24};

Int_t TOF2_Mapping1[12] = {1,5,9,13,17,21,25,29,33,37,41,45};
Int_t TOF2_Mapping2[12] = {2,6,10,14,18,22,26,30,34,38,42,46};
Int_t TOF2_Mapping3[12] = {3,7,11,15,19,23,27,31,35,39,43,47};
Int_t TOF2_Mapping4[12] = {4,8,12,16,20,24,28,32,36,40,44,48};

//char path_input[100];
//sprintf(path_input,"/triumfcs/trshare/trek/E36/Data/April_2015/Histograms");

char path_input[100];
sprintf(path_input,path_histo);

char Name_finput[100];
//sprintf(path_input,"/media/bianchin/hdd1/trek/E36/Data/April_2015/Histograms");
sprintf(Name_finput,"%s/Hist_Run%dMSC.root",path_input, run_number);
cout << Name_finput << endl;
cout << "     " << endl;



TFile *finput = new TFile(Name_finput);

if(board_number==1){
	for(int j=0;j<64;j++){
		sprintf(Name_ADC_High,"ADC_High (Ch. %d) - SFT",j);
		sprintf(Name_ADC_Low,"ADC_Low (Ch. %d) - SFT",j);
		sprintf(Name_TDC_LE,"TDC_LE (Ch. %d) - SFT",j);
		sprintf(Name_TDC_TE,"TDC_TE (Ch. %d) - SFT",j);
	
		h_ADC_High[j]=(TH1F*)finput->FindObjectAny(Name_ADC_High);
		h_ADC_Low[j]=(TH1F*)finput->FindObjectAny(Name_ADC_Low);
		h_TDC_LE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_LE);
		h_TDC_TE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_TE);
	}

}



if(board_number==2){
	for(int j=0;j<64;j++){
		sprintf(Name_ADC_High,"ADC_High (Ch. %d) - SFT",j+64);
		sprintf(Name_ADC_Low,"ADC_Low (Ch. %d) - SFT",j+64);
		sprintf(Name_TDC_LE,"TDC_LE (Ch. %d) - SFT",j+64);
		sprintf(Name_TDC_TE,"TDC_TE (Ch. %d) - SFT",j+64);
	
		h_ADC_High[j]=(TH1F*)finput->FindObjectAny(Name_ADC_High);
		h_ADC_Low[j]=(TH1F*)finput->FindObjectAny(Name_ADC_Low);
		h_TDC_LE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_LE);
		h_TDC_TE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_TE);
	}
}

if(board_number==3){
	for(int j=0;j<64;j++){
		sprintf(Name_ADC_High,"ADC_High (Ch. %d) - TARGET",j);
		sprintf(Name_ADC_Low,"ADC_Low (Ch. %d) - TARGET",j);
		sprintf(Name_TDC_LE,"TDC_LE (Ch. %d) - TARGET",j);
		sprintf(Name_TDC_TE,"TDC_TE (Ch. %d) - TARGET",j);
	
		h_ADC_High[j]=(TH1F*)finput->FindObjectAny(Name_ADC_High);
		h_ADC_Low[j]=(TH1F*)finput->FindObjectAny(Name_ADC_Low);
		h_TDC_LE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_LE);
		h_TDC_TE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_TE);
	}
}

if(board_number==4){
	for(int j=0;j<64;j++){
		sprintf(Name_ADC_High,"ADC_High (Ch. %d) - TARGET",j+64);
		sprintf(Name_ADC_Low,"ADC_Low (Ch. %d) - TARGET",j+64);
		sprintf(Name_TDC_LE,"TDC_LE (Ch. %d) - TARGET",j+64);
		sprintf(Name_TDC_TE,"TDC_TE (Ch. %d) - TARGET",j+64);

		h_ADC_High[j]=(TH1F*)finput->FindObjectAny(Name_ADC_High);
		h_ADC_Low[j]=(TH1F*)finput->FindObjectAny(Name_ADC_Low);
		h_TDC_LE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_LE);
		h_TDC_TE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_TE);
	}
}

if(board_number==5){
	for(int j=0;j<64;j++){
		sprintf(Name_ADC_High,"ADC_High (Ch. %d) - TARGET",j+128);
		sprintf(Name_ADC_Low,"ADC_Low (Ch. %d) - TARGET",j+128);
		sprintf(Name_TDC_LE,"TDC_LE (Ch. %d) - TARGET",j+128);
		sprintf(Name_TDC_TE,"TDC_TE (Ch. %d) - TARGET",j+128);
		
		h_ADC_High[j]=(TH1F*)finput->FindObjectAny(Name_ADC_High);
		h_ADC_Low[j]=(TH1F*)finput->FindObjectAny(Name_ADC_Low);
		h_TDC_LE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_LE);
		h_TDC_TE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_TE);
	}
}

if(board_number==6){
	for(int j=0;j<64;j++){
		sprintf(Name_ADC_High,"ADC_High (Ch. %d) - TARGET",j+192);
		sprintf(Name_ADC_Low,"ADC_Low (Ch. %d) - TARGET",j+192);
		sprintf(Name_TDC_LE,"TDC_LE (Ch. %d) - TARGET",j+192);
		sprintf(Name_TDC_TE,"TDC_TE (Ch. %d) - TARGET",j+192);
		
		h_ADC_High[j]=(TH1F*)finput->FindObjectAny(Name_ADC_High);
		h_ADC_Low[j]=(TH1F*)finput->FindObjectAny(Name_ADC_Low);
		h_TDC_LE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_LE);
		h_TDC_TE[j]=(TH1F*)finput->FindObjectAny(Name_TDC_TE);
	}
}

if(board_number==7 || board_number==8){
	for(int j=0;j<24;j++){
		sprintf(Name_ADC_tof1,"ADC_TOF1 (Ch. %d)",j);

		h_ADC_tof1[j]=(TH1F*)finput->FindObjectAny(Name_ADC_tof1);
	}

	for(int j=0;j<12;j++){
		sprintf(Name_TDC_tof1U,"TDC_TOF1 (Ch. %d)",j);
		sprintf(Name_TDC_tof1D,"TDC_TOF1 (Ch. %d)",j+12);

		h_TDC_tof1U[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof1U);
		h_TDC_tof1D[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof1D);
	}
}

if(board_number==9 || board_number==10){
	for(int j=0;j<56;j++){
		sprintf(Name_ADC_tof2,"ADC_TOF2 (Ch. %d)",j);	

		h_ADC_tof2[j]=(TH1F*)finput->FindObjectAny(Name_ADC_tof2);
	}

	for(int j=0;j<12;j++){
		sprintf(Name_TDC_tof2AO,"TDC_TOF2 (Ch. %d)",j);
		sprintf(Name_TDC_tof2AI,"TDC_TOF2 (Ch. %d)",j+24);
		sprintf(Name_TDC_tof2BO,"TDC_TOF2 (Ch. %d)",j+12);
		sprintf(Name_TDC_tof2BI,"TDC_TOF2 (Ch. %d)",j+36);

		h_TDC_tof2AO[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof2AO);
		h_TDC_tof2AI[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof2AI);
		h_TDC_tof2BO[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof2BO);
		h_TDC_tof2BI[j]=(TH1F*)finput->FindObjectAny(Name_TDC_tof2BI);
	}
}

if(board_number==11){
	for(int j=0;j<512;j++){
		sprintf(Name_ADC_MWPC,"ADC_MWPC (Ch. %d)",j);	

		h_ADC_MWPC[j]=(TH1F*)finput->FindObjectAny(Name_ADC_MWPC);
	}
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
sprintf(Name_Can_ADC_tof1,"TOF1 ADCs & TDCs -- Run %d  (Ch. 0 - 23)",run_number);
sprintf(Title_Can_ADC_tof1,"TOF1 ADCs & TDCs -- Run %d  (Ch. 0 - 23)",run_number);
}

if(board_number==8){
sprintf(Name_Can_ADC_tof1B,"TOF1 ADCs & TDCs -- Run %d  (Gap 1 - 12)  |  1U, 1D, 2U, 2D, etc...",run_number);
sprintf(Title_Can_ADC_tof1B,"TOF1 ADCs & TDCs -- Run %d  (Gap 0 - 12)  |  1U, 1D, 2U, 2D, etc...",run_number);
}

if(board_number==9){
sprintf(Name_Can_ADC_tof2,"TOF2 ADCs -- Run %d  (Ch. 0 - 55)",run_number);
sprintf(Title_Can_ADC_tof2,"TOF2 ADCs -- Run %d  (Ch. 0 - 55)",run_number);
sprintf(Name_Can_ADC_tof2D,"TOF2 TDCs -- Run %d  (Ch. 0 - 47)",run_number);
sprintf(Title_Can_ADC_tof2D,"TOF2 TDCs -- Run %d  (Ch. 0 - 47)",run_number);
}

if(board_number==10){
sprintf(Name_Can_ADC_tof2B,"TOF2 ADCs & TDCs -- Run %d  (Gap 1 - 6)  |  1AO, 1AI, 1BO, 1BI, etc...",run_number);
sprintf(Title_Can_ADC_tof2B,"TOF2 ADCs & TDCs -- Run %d  (Gap 1 - 6)  |  1AO, 1AI, 1BO, 1BI, etc...",run_number);
sprintf(Name_Can_ADC_tof2C,"TOF2 ADCs & TDCs -- Run %d  (Gap 7 - 12)  |  7AO, 7AI, 7BO, 7BI, etc...",run_number);
sprintf(Title_Can_ADC_tof2C,"TOF2 ADCs & TDCs -- Run %d  (Gap 7 - 12)",run_number);
}

if(board_number==11){
sprintf(Name_Can_ADC_MWPC,"MWPC ADCs -- Run %d  (Ch. 0 - 255)",run_number);
sprintf(Title_Can_ADC_MWPC,"MWPC ADCs -- Run %d  (Ch. 0 - 255)",run_number);

sprintf(Name_Can_ADC_MWPC2,"MWPC ADCs -- Run %d  (Ch. 256 - 511)",run_number);
sprintf(Title_Can_ADC_MWPC2,"MWPC ADCs -- Run %d  (Ch. 256 - 511)",run_number);
}


////////////////////////////////////////////////
/*
if (board_number <= 6){
	if(iflag==1){
		TCanvas *c1;
		c1 = new TCanvas(Name_Can_TDC_TE,Title_Can_TDC_TE,1200,500); 
		c1->Divide(8,8);
		c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
	}
}

if (board_number <= 6){
//if(iflag==1){
	TCanvas *c2;
	c2 = new TCanvas(Name_Can_TDC_LE,Title_Can_TDC_LE,1200,500); 
	c2->Divide(8,8);
	c2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
//}
}

if (board_number <= 6){
	if(iflag==1){
		TCanvas *c3;
		c3 = new TCanvas(Name_Can_ADC_Low,Title_Can_ADC_Low,1200,500); 
		c3->Divide(8,8);
		c3->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
	}
}

if (board_number <= 6){
//if(iflag==1){
	TCanvas *c4;
	c4 = new TCanvas(Name_Can_ADC_High,Title_Can_ADC_High,1200,500); 
	c4->Divide(8,8);
	c4->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
//}
}

if(board_number==7){
	TCanvas *c5;
	c5 = new TCanvas(Name_Can_ADC_tof1,Title_Can_ADC_tof1,1200,500); 
	c5->Divide(8,3);
	c5->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
}

if(board_number==8){
	TCanvas *c6;
	c6 = new TCanvas(Name_Can_ADC_tof2,Title_Can_ADC_tof2,1200,500); 
	c6->Divide(8,6);
	c6->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
}
*/

////////////////////////////////////////////////////////////

if (board_number <= 6) {
	for(Int_t ii=0; ii<32; ii++){
		max_Can2=h_TDC_LE[ii]->GetMaximum();
		max_Can1=h_TDC_LE[ii+32]->GetMaximum();
		if(max_Can2 > TDC_Max_Can2){
			TDC_Max_Can2 = max_Can2;
//			cout << ii << "  " << TDC_Max_Can2 << endl;
		}
		
		if(max_Can1 > TDC_Max_Can1){
			TDC_Max_Can1 = max_Can1;
//			cout << ii+32 << "  " << TDC_Max_Can1 << endl;
		}
	}
}
	//cout << "Finally :   " << TDC_Max_Can2 << "  " << TDC_Max_Can1 << endl;



/*
if(iflag==1){

	TCanvas *c1;
	c1 = new TCanvas(Name_Can_TDC_TE,Title_Can_TDC_TE,1200,500); 
	c1->Divide(8,8);
	c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=0; ican<32; ican++){
		if(board_number==1 && (ican==23 || ican==24 || ican==25 || ican==26 ||ican==27 || ican==28 || ican==29)) c1->cd(ican+1)->SetFillColor(45);
		if(board_number==2 && (ican==32 || ican==33 || ican==34 || ican==35 ||ican==36 || ican==37 || ican==38 ||  ican==39 || ican==52)) c1->cd(ican+1)->SetFillColor(45);
		if(run_number>=1680 && board_number==1 && (ican>=32 && ican<=39)) c1->cd(ican+1)->SetFillColor(kYellow-9);
		if(run_number>=1680 && board_number==2 && (ican>=8 && ican<=15)) c1->cd(ican+1)->SetFillColor(kYellow-9);
		if(run_number>=1680 && board_number==4 && (ican>=24 && ican<=31)) c1->cd(ican+1)->SetFillColor(kYellow-9);
		if(run_number>=1680 && board_number==5 && (ican>=8 && ican<=15)) c1->cd(ican+1)->SetFillColor(kYellow-9);
		if(run_number>=1680 && board_number==6 && (ican>=16 && ican<=23)) c1->cd(ican+1)->SetFillColor(kYellow-9);
		else c1->cd(ican+1);
		h_TDC_TE[ican]->SetAxisRange(0, TDC_Graph_xmax,"X");
		if(TDC_Graph_ymax > 0) {  
			h_TDC_TE[ican]->SetAxisRange(0, TDC_Graph_ymax,"Y");
		}
		h_TDC_TE[ican]->Draw();
	}
}
*/

if ((board_number <= 6) && (iflag==0)) {

	TCanvas *c4;
	c4 = new TCanvas(Name_Can_ADC_High2,Title_Can_ADC_High2,1200,500); 
	c4->Divide(8,8);
	c4->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=0; ican<32; ican++){
		//if(board_number==1 && (ican==23 || ican==24 || ican==25 || ican==26 ||ican==27 || ican==28 || ican==29)) c4->cd(ican+1)->SetFillColor(45);
		if(board_number==2 && (ican==0 || ican==1 || ican==2 || ican==3 ||ican==4 || ican==5 || ican==6 ||  ican==7 || ican==20)) c4->cd(ican+1)->SetFillColor(45);
		if(board_number==2 && (ican==0 || ican==1 || ican==2 || ican==3 ||ican==4 || ican==5 || ican==6 ||  ican==7 || ican==20)) c4->cd(ican+1+32)->SetFillColor(45);
		if(run_number>=1680 && board_number==1 && (ican>=0 && ican<=7)) c4->cd(ican+1)->SetFillColor(kYellow-9);
		if(run_number>=1680 && board_number==1 && (ican>=0 && ican<=7)) c4->cd(ican+1+32)->SetFillColor(kYellow-9);
		//if(run_number>=1680 && board_number==2 && (ican>=8 && ican<=15)) c4->cd(ican+1)->SetFillColor(kYellow-9);
		//if(run_number>=1680 && board_number==4 && (ican>=24 && ican<=31)) c4->cd(ican+1)->SetFillColor(kYellow-9);
		//if(run_number>=1680 && board_number==5 && (ican>=8 && ican<=15)) c4->cd(ican+1)->SetFillColor(kYellow-9);
		//if(run_number>1591 && board_number==6 && (ican>=16 && ican<=23)) c4->cd(ican+1)->SetFillColor(kYellow-9);
		else c4->cd(ican+1);
		c4->cd(ican+1)->SetLogy();  
		h_ADC_High[ican+32]->Draw();
		c4->cd(ican+1+32);
		h_TDC_LE[ican+32]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
	//	cout << "TDC_Max_Can2   :  " << TDC_Max_Can2 << endl;
	//	h_TDC_LE[ican+32]->SetMaximum(TDC_Max_Can2);
		if(TDC_Graph_ymax > 0) {  
			h_TDC_LE[ican+32]->SetAxisRange(0, TDC_Graph_ymax,"Y");
		}	
		h_TDC_LE[ican+32]->SetMaximum(1.1*TDC_Max_Can1);
		h_TDC_LE[ican+32]->Draw();
//		c4->cd(ican+1)->SetLogy();  
//		h_ADC_High[ican+32]->Draw();
	}
}


if ((board_number <= 6) && (iflag==0)){

	TCanvas *c2;
	c2 = new TCanvas(Name_Can_TDC_LE,Title_Can_TDC_LE,0,200,1200,500); 
	c2->Divide(8,8);
	c2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=0; ican<32; ican++){
		if(board_number==1 && (ican==23 || ican==24 || ican==25 || ican==26 ||ican==27 || ican==28 || ican==29)) c2->cd(ican+1)->SetFillColor(45);
		if(board_number==1 && (ican==23 || ican==24 || ican==25 || ican==26 ||ican==27 || ican==28 || ican==29)) c2->cd(ican+1+32)->SetFillColor(45);
		//if(board_number==2 && (ican==32 || ican==33 || ican==34 || ican==35 ||ican==36 || ican==37 || ican==38 ||  ican==39 || ican==52)) c2->cd(ican+1)->SetFillColor(45);
		if(run_number>=1680 && board_number==1 && (ican>=32 && ican<=39)) c2->cd(ican+1)->SetFillColor(kYellow-9);
		if(run_number>=1680 && board_number==2 && (ican>=8 && ican<=15)) c2->cd(ican+1)->SetFillColor(kYellow-9);
		if(run_number>=1680 && board_number==2 && (ican>=8 && ican<=15)) c2->cd(ican+1+32)->SetFillColor(kYellow-9);
		if(run_number>=1680 && board_number==4 && (ican>=24 && ican<=31)) c2->cd(ican+1)->SetFillColor(kYellow-9);
		if(run_number>=1680 && board_number==4 && (ican>=24 && ican<=31)) c2->cd(ican+1+32)->SetFillColor(kYellow-9);
		if(run_number>=1680 && board_number==5 && (ican>=8 && ican<=15)) c2->cd(ican+1)->SetFillColor(kYellow-9);
		if(run_number>=1680 && board_number==5 && (ican>=8 && ican<=15)) c2->cd(ican+1+32)->SetFillColor(kYellow-9);
		if(run_number>=1680 && board_number==6 && (ican>=16 && ican<=23)) c2->cd(ican+1)->SetFillColor(kYellow-9);
		if(run_number>=1680 && board_number==6 && (ican>=16 && ican<=23)) c2->cd(ican+1+32)->SetFillColor(kYellow-9);
		else c2->cd(ican+1);
		c2->cd(ican+1)->SetLogy();
		h_ADC_High[ican]->Draw();
		c2->cd(ican+1+32);
		h_TDC_LE[ican]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
//		cout << "TDC_Max_Can1   :  " << TDC_Max_Can1 << endl;
//		h_TDC_LE[ican]->SetMaximum(TDC_Max_Can1);
		if(TDC_Graph_ymax > 0) {  
			h_TDC_LE[ican]->SetAxisRange(0, TDC_Graph_ymax,"Y");
		}
		//if(board_number==1 && (ican==0 || ican==1)) h_TDC_LE[ican]->SetFillColor(2);
		//c2->cd(ican+1+32);
		h_TDC_LE[ican]->SetMaximum(1.1*TDC_Max_Can2);
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
		if(board_number==1 && (ican==23 || ican==24 || ican==25 || ican==26 ||ican==27 || ican==28 || ican==29)) c3->cd(ican+1)->SetFillColor(45);
		if(board_number==2 && (ican==32 || ican==33 || ican==34 || ican==35 ||ican==36 || ican==37 || ican==38 ||  ican==39 || ican==52)) c3->cd(ican+1)->SetFillColor(45);
		if(run_number>=1680 && board_number==1 && (ican>=32 && ican<=39)) c3->cd(ican+1)->SetFillColor(kYellow-9);
		if(run_number>=1680 && board_number==2 && (ican>=8 && ican<=15)) c3->cd(ican+1)->SetFillColor(kYellow-9);
		if(run_number>=1680 && board_number==4 && (ican>=24 && ican<=31)) c3->cd(ican+1)->SetFillColor(kYellow-9);
		if(run_number>=1680 && board_number==5 && (ican>=8 && ican<=15)) c3->cd(ican+1)->SetFillColor(kYellow-9);
		if(run_number>=1680 && board_number==6 && (ican>=16 && ican<=23)) c3->cd(ican+1)->SetFillColor(kYellow-9);
		else c3->cd(ican+1);
		c3->cd(ican+1)->SetLogy();  
		h_ADC_Low[ican+32]->Draw();
		c3->cd(ican+1+32);
		h_TDC_LE[ican+32]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
		h_TDC_LE[ican+32]->Draw();
	}
}

///********************TOF1

if (board_number == 7) {
	
	TCanvas *c5;
	c5 = new TCanvas(Name_Can_ADC_tof1,Title_Can_ADC_tof1,1200,500); 
	c5->Divide(8,6);
	c5->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=0; ican<24; ican++){
		c5->cd(ican+1)->SetLogy();

//		h_ADC_tof1[ican]->SetAxisRange(0, TOF1_Graph_xmax,"X");
//		if(TOF1_Graph_ymax > 0) {  
//			h_ADC_tof1[ican]->SetAxisRange(0, TOF1_Graph_ymax,"Y");
//		}
		h_ADC_tof1[ican]->Draw();
	}

	for(Int_t ican=0; ican<12; ican++){
		c5->cd(ican+25);

//		h_ADC_tof1[ican]->SetAxisRange(0, TOF1_Graph_xmax,"X");
//		if(TOF1_Graph_ymax > 0) {  
//			h_ADC_tof1[ican]->SetAxisRange(0, TOF1_Graph_ymax,"Y");
//		}
		h_TDC_tof1U[ican]->Draw();
	}

	for(Int_t ican=0; ican<12; ican++){
		c5->cd(ican+37);

//		h_ADC_tof1[ican]->SetAxisRange(0, TOF1_Graph_xmax,"X");
//		if(TOF1_Graph_ymax > 0) {  
//			h_ADC_tof1[ican]->SetAxisRange(0, TOF1_Graph_ymax,"Y");
//		}
		h_TDC_tof1D[ican]->Draw();
	}
}

if (board_number == 8) {
	
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
//		if(TOF1_Graph_ymax > 0) {  
//			h_ADC_tof1[ican]->SetAxisRange(0, TOF1_Graph_ymax,"Y");
//		}
		h_TDC_tof1U[ican]->SetTitle(TOF1Title3);
		h_TDC_tof1U[ican]->Draw();
	}

	for(Int_t ican=0; ican<12; ican++){
		c7->cd(TOF1_Mapping2[ican]+24);

		char TOF1Title4[250];
		sprintf(TOF1Title4, "Raw TDC (Ch. %d) - TOF1 | Down-%d : Run %d", ican+12, ican+1, run_number);
//		h_ADC_tof1[ican]->SetAxisRange(0, TOF1_Graph_xmax,"X");
//		if(TOF1_Graph_ymax > 0) {  
//			h_ADC_tof1[ican]->SetAxisRange(0, TOF1_Graph_ymax,"Y");
//		}
		h_TDC_tof1D[ican]->SetTitle(TOF1Title4);
		h_TDC_tof1D[ican]->Draw();
	}
}

///********************TOF2

if (board_number == 9) {

	TCanvas *c6;
	c6 = new TCanvas(Name_Can_ADC_tof2,Title_Can_ADC_tof2,1200,500); 
	c6->Divide(8,7);
	c6->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=0; ican<56; ican++){
		c6->cd(ican+1)->SetLogy();
//		h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_xmax,"X");
//		if(TOF2_Graph_ymax > 0) {  
//			h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_ymax,"Y");
//		}
		h_ADC_tof2[ican]->Draw();
	}

	TCanvas *c11;
	c11 = new TCanvas(Name_Can_ADC_tof2D,Title_Can_ADC_tof2D,1200,500); 
	c11->Divide(8,6);
	c11->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=0; ican<12; ican++){
		c11->cd(ican+1);
//		h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_xmax,"X");
//		if(TOF2_Graph_ymax > 0) {  
//			h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_ymax,"Y");
//		}
		h_TDC_tof2AO[ican]->Draw();
	}

	for(Int_t ican=0; ican<12; ican++){
		c11->cd(ican+25);
//		h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_xmax,"X");
//		if(TOF2_Graph_ymax > 0) {  
//			h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_ymax,"Y");
//		}
		h_TDC_tof2AI[ican]->Draw();
	}

	for(Int_t ican=0; ican<12; ican++){
		c11->cd(ican+13);
//		h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_xmax,"X");
//		if(TOF2_Graph_ymax > 0) {  
//			h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_ymax,"Y");
//		}
		h_TDC_tof2BO[ican]->Draw();
	}

	for(Int_t ican=0; ican<12; ican++){
		c11->cd(ican+37);
//		h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_xmax,"X");
//		if(TOF2_Graph_ymax > 0) {  
//			h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_ymax,"Y");
//		}
		h_TDC_tof2BI[ican]->Draw();
	}



}

if (board_number == 10) {

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
		h_TDC_tof2AO[ican+6]->SetTitle(TOF2Title5);
		h_TDC_tof2AO[ican+6]->Draw();
	}

	for(Int_t ican=0; ican<6; ican++){
		c12->cd(TOF2_Mapping2[ican+6]);
		
		char TOF2Title6[250];
		sprintf(TOF2Title6, "Raw TDC (Ch. %d) - TOF2 | InA-%d : Run %d", ican+24, ican+7, run_number);
		h_TDC_tof2AI[ican+6]->SetTitle(TOF2Title6);
		h_TDC_tof2AI[ican+6]->Draw();
	}
	
	for(Int_t ican=0; ican<6; ican++){
		c12->cd(TOF2_Mapping3[ican+6]);

		char TOF2Title7[250];
		sprintf(TOF2Title7, "Raw TDC (Ch. %d) - TOF2 | OutB-%d : Run %d", ican+12, ican+7, run_number);
		h_TDC_tof2BO[ican+6]->SetTitle(TOF2Title7);
		h_TDC_tof2BO[ican+6]->Draw();
	}

	for(Int_t ican=0; ican<6; ican++){
		c12->cd(TOF2_Mapping4[ican+6]);

		char TOF2Title8[250];
		sprintf(TOF2Title8, "Raw TDC (Ch. %d) - TOF2 | InB-%d : Run %d", ican+36, ican+7, run_number);
		h_TDC_tof2BI[ican+6]->SetTitle(TOF2Title8);
		h_TDC_tof2BI[ican+6]->Draw();
	}


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
		h_TDC_tof2AO[ican]->SetTitle(TOF2Title5);
		h_TDC_tof2AO[ican]->Draw();
	}

	for(Int_t ican=0; ican<6; ican++){
		c8->cd(TOF2_Mapping2[ican+6]);
		
		char TOF2Title6[250];
		sprintf(TOF2Title6, "Raw TDC (Ch. %d) - TOF2 | InA-%d : Run %d", ican+24, ican+1, run_number);
		h_TDC_tof2AI[ican]->SetTitle(TOF2Title6);
		h_TDC_tof2AI[ican]->Draw();
	}
	
	for(Int_t ican=0; ican<6; ican++){
		c8->cd(TOF2_Mapping3[ican+6]);

		char TOF2Title7[250];
		sprintf(TOF2Title7, "Raw TDC (Ch. %d) - TOF2 | OutB-%d : Run %d", ican+12, ican+1, run_number);
		h_TDC_tof2BO[ican]->SetTitle(TOF2Title7);
		h_TDC_tof2BO[ican]->Draw();
	}

	for(Int_t ican=0; ican<6; ican++){
		c8->cd(TOF2_Mapping4[ican+6]);

		char TOF2Title8[250];
		sprintf(TOF2Title8, "Raw TDC (Ch. %d) - TOF2 | InB-%d : Run %d", ican+36, ican+1, run_number);
		h_TDC_tof2BI[ican]->SetTitle(TOF2Title8);
		h_TDC_tof2BI[ican]->Draw();
	}
}

///********************MWPC

if (board_number == 11) {

	TCanvas *c9;
	c9 = new TCanvas(Name_Can_ADC_MWPC,Title_Can_ADC_MWPC,1400,1000); 
	c9->Divide(16,16);
	c9->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=0; ican<256; ican++){
		c9->cd(ican+1)->SetLogy();
//		h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_xmax,"X");
//		if(TOF2_Graph_ymax > 0) {  
//			h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_ymax,"Y");
//		}
		h_ADC_MWPC[ican]->Draw();
	}

	TCanvas *c10;
	c10 = new TCanvas(Name_Can_ADC_MWPC2,Title_Can_ADC_MWPC2,0,200,1400,1000); 
	c10->Divide(16,16);
	c10->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=0; ican<256; ican++){
		c10->cd(ican+1)->SetLogy();
//		h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_xmax,"X");
//		if(TOF2_Graph_ymax > 0) {  
//			h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_ymax,"Y");
//		}
		h_ADC_MWPC[ican+256]->Draw();
	}
}

return;
}