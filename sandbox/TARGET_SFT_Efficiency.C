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
#include <math.h>
//#include "TSpectrum.h"
#include "TMarker.h"
//#include "Event_Display_MS.h"
#include "ANAPATH.h"
#include "Thresholds.h"
#include "CommonParameters.h"
#endif
  
void TARGET_SFT_Efficiency(Int_t Run_Number=5, Int_t Enable_TARGET_Eff=0, Int_t nevents=0) { 

Int_t thr;					thr = ADC_cut_SFT;
Int_t TDCthr1;				TDCthr1 = TDC_min_SFT;
Int_t TDCthr2;				TDCthr2 = TDC_max_SFT;

gStyle->SetOptStat(0);

Int_t ADC_cut_TARGET2=850;
Int_t ADC_cut_TARGET3=1250;

//Int_t on_track_counter[256] = {0};
//Int_t track_hit_counter[256] = {0};
Int_t bar_hit_counter[256] = {0};
Int_t bar_adj_counter[256] = {0};
//Int_t bar_track_counter[256] = {0};
float Efficiency[256] = {0};
Int_t event_counter = 0;
Int_t good_event=0;

Int_t Angle_ADC_cut = 0;	
Int_t path_margin = 5;
Int_t cut_margin = 12;
Int_t height_margin = 5;
Int_t width_margin = 3;

Int_t adc_high_target[256];   		Int_t ADC_High_TARGET[256];    	//Double_t ADC_High_TARGET_corr[256];   	
Int_t adc_low_target[256]; 			Int_t ADC_Low_TARGET[256]; 	//Double_t ADC_Low_TARGET_corr[256];
Int_t tdc_le_target[256][16];      	Int_t TDC_LE_TARGET[256];    	//Double_t TDC_LE_TARGET_corr[256];   	
Int_t tdc_te_target[256][16]; 		Int_t TDC_TE_TARGET[256]; 	//Double_t TDC_TE_TARGET_corr[256];

Int_t adc_high_sft[128];     	 Int_t ADC_High_sft[128];      Double_t ADC_High_SFT_corr[128];    
Int_t adc_low_sft[128];    		 Int_t ADC_Low_SFT[128];   //Double_t ADC_Low_SFT_corr[128];
Int_t tdc_le_sft[128][16];       Int_t TDC_LE_sft[128][16];      //Double_t TDC_LE_SFT_corr[128];    
Int_t tdc_te_sft[128][16];   	 Int_t TDC_TE_sft[128][16];    //Double_t TDC_TE_SFT_corr[128];  

Int_t ADC_tof1[24];	Int_t ADC_TOF1[24];	
Int_t ADC_tof2[56];	Int_t ADC_TOF2[56];

Int_t TDC_tof1U[12];  Int_t TDC_TOF1U[12];
Int_t TDC_tof1D[12];  Int_t TDC_TOF1D[12];

Int_t TDC_tof2AO[12];   Int_t TDC_TOF2AO[12];
Int_t TDC_tof2BO[12];   Int_t TDC_TOF2BO[12];
Int_t TDC_tof2AI[12];   Int_t TDC_TOF2AI[12];
Int_t TDC_tof2BI[12];   Int_t TDC_TOF2BI[12];

Int_t MwpcADC[512];	Int_t MWPCADC[512];

Int_t fiber[128]={-1};
Int_t fiberTDC[128]={-1};

/////////////////////////////////////////////////////////////////////////////////////////////////

int channel0[8] = {5,13,14,15,26,27,40,41};
int channel1[8] = {27,40,41,56,57,73,90,91};
int channel2[8] = {90,91,109,127,145,163,180,181};
int channel3[8] = {180,181,197,212,213,226,227,239};
int channel4[8] = {212,213,238,239,247,248,249,255};
int channel5[8] = {242,247,250,251,252,253,254,255};
int channel6[8] = {214,215,228,229,241,242,243,250};
int channel7[8] = {164,165,182,198,199,214,215,228};
int channel8[8] = {74,75,92,110,128,146,164,165};
int channel9[8] = {16,28,29,42,43,58,74,75};
int channel10[8] = {0,6,7,8,16,17,28,29};
int channel11[8] = {0,1,2,3,4,5,8,13};


////////////////////

Int_t counter1 = 0;
Int_t counter2 = 0;
Int_t counter3 = 0;
Int_t counter4 = 0;
Int_t counter5 = 0;
Int_t counter6 = 0;
Int_t counter7 = 0;
Int_t counter8 = 0;

////////////////////

Int_t counter1a = 0;
Int_t counter2a = 0;
Int_t counter3a = 0;
Int_t counter4a = 0;
Int_t counter5a = 0;
Int_t counter6a = 0;
Int_t counter7a = 0;
Int_t counter8a = 0;

/////////////////////

Int_t helicity1or = 0;
Int_t helicity2or = 0;
Int_t helicity1and = 0;
Int_t helicity2and = 0;

/////////////////////

Int_t helicity1ora = 0;
Int_t helicity2ora = 0;
Int_t helicity1anda = 0;
Int_t helicity2anda = 0;

/////////////////////

Int_t count_or = 0;
Int_t count_and = 0;

Int_t count_hel_or = 0;
Int_t count_hel_and = 0;

///////////////////////

Int_t count_ora = 0;
Int_t count_anda = 0;

Int_t count_hel_ora = 0;
Int_t count_hel_anda = 0;

///////////////////////

Int_t onehitcount = 0;

//////////////////////

Int_t countlayer1or = 0;
Int_t countlayer1and = 0;
Int_t countlayer2or = 0;
Int_t countlayer2and = 0;
Int_t countlayer3or = 0;
Int_t countlayer3and = 0;
Int_t countlayer4or = 0;
Int_t countlayer4and = 0;
Int_t counthelicity1 = 0;
Int_t counthelicity2 = 0;

/////////////////////

Int_t countlayer1_2[15] = {0};
Int_t countlayer3_4[17] = {0};

Int_t counterp1 = 0;
Int_t counterp2 = 0;
Int_t counterp3 = 0;
Int_t counterp4 = 0;

Int_t p1_p2 = 0;
Int_t p3_p4 = 0;

/////////////////////////////////////////////////////////////////////////////////////////////////////

char path_input[200];                   char file_mapping[200];
sprintf(path_input,"%s",path_merged);          sprintf(file_mapping,"../Mapping");
//sprintf(path_input,"/triumfcs/trshare/trek/E36/Data/April_2015/root");

char Name_finput[200];
sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);

char par_map[200];		char par_finput[200];
//sprintf(par_finput,"/media/bianchin/hdd1/trek/E36/Macros/April_2015/Mapping/SFT_Mapping_Oct14.txt");
sprintf(par_map,"../Mapping");
sprintf(par_finput,"%s/SFT_Mapping_Oct14.txt",par_map);

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

cout << "   " << endl;
cout << Name_finput << endl;

TChain *fChain= new TChain("Tree");		
fChain->Add(Name_finput);		
fChain->SetMakeClass(1);							

fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);		fChain->SetBranchAddress("ADC_High_SFT",ADC_High_sft);
fChain->SetBranchAddress("ADC_Low_TARGET",adc_low_target);		fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);		fChain->SetBranchAddress("TDC_LE_SFT",TDC_LE_sft);
fChain->SetBranchAddress("TDC_TE_TARGET",tdc_te_target);		fChain->SetBranchAddress("TDC_TE_SFT",TDC_TE_sft);

fChain->SetBranchAddress("ADC_TOF1",ADC_TOF1);
fChain->SetBranchAddress("ADC_TOF2",ADC_TOF2);

fChain->SetBranchAddress("TDC_TOF1U",TDC_tof1U);
fChain->SetBranchAddress("TDC_TOF1D",TDC_tof1D);
fChain->SetBranchAddress("TDC_TOF2AO",TDC_tof2AO);
fChain->SetBranchAddress("TDC_TOF2AI",TDC_tof2AI);
fChain->SetBranchAddress("TDC_TOF2BO",TDC_tof2BO);
fChain->SetBranchAddress("TDC_TOF2BI",TDC_tof2BI);    

fChain->SetBranchAddress("MWPCADC",MwpcADC);		

Int_t nentries = (Int_t)fChain->GetEntries();
cout << "  " << endl;
cout << "****  Number of events: " << nentries << "  **** " <<endl;
cout << "  " << endl;

cout << "   " << endl;
cout << "Mapping File:   " << par_finput << endl;
cout << "" << endl;

cout << "   " << endl;

///////////////////////////////////////////////////////////

char h_target_ADC_title[100];
sprintf(h_target_ADC_title,"Run %d -- ADC > %d",Run_Number, ADC_cut_TARGET);
TH2D *h_target_ADC = new TH2D(h_target_ADC_title,h_target_ADC_title,700,-35,35,700,-35,35);

char h_target_ADCA_title[100];
sprintf(h_target_ADCA_title,"Run %d -- ADC > %d -- Outlier Cut",Run_Number, ADC_cut_TARGET);
TH2D *h_target_ADCA = new TH2D(h_target_ADCA_title,h_target_ADCA_title,700,-35,35,700,-35,35);

TH1D *h_count;	char Title_h_count[100];	char Name_h_count[100];
TH1D *h_count2;	char Title_h_count2[100];	char Name_h_count2[100];
TH1D *h_count3;	char Title_h_count3[100];	char Name_h_count3[100];
TH1D *h_count4;	char Title_h_count4[100];	char Name_h_count4[100];
TH1D *h_count5;	char Title_h_count5[100];	char Name_h_count5[100];
TH1D *h_count6;	char Title_h_count6[100];	char Name_h_count6[100];
TH1D *h_count7;	char Title_h_count7[100];	char Name_h_count7[100];
TH1D *h_count8;	char Title_h_count8[100];	char Name_h_count8[100];
TH1D *h_count9;	char Title_h_count9[100];	char Name_h_count9[100];
TH1D *h_count10; char Title_h_count10[100];	char Name_h_count10[100];
TH1D *h_count11; char Title_h_count11[100];	char Name_h_count11[100];
TH1D *h_count12; char Title_h_count12[100];	char Name_h_count12[100];
TH1D *h_count13; char Title_h_count13[100];	char Name_h_count13[100];
TH1D *h_count14; char Title_h_count14[100];	char Name_h_count14[100];
TH1D *h_count15; char Title_h_count15[100];	char Name_h_count15[100];
TH1D *h_count16; char Title_h_count16[100];	char Name_h_count16[100];
TH1D *h_count17; char Title_h_count17[100];	char Name_h_count17[100];
TH1D *h_count18; char Title_h_count18[100];	char Name_h_count18[100];
TH1D *h_count19; char Title_h_count19[100];	char Name_h_count19[100];
TH1D *h_count20; char Title_h_count20[100];	char Name_h_count20[100];
TH1D *h_good_events;  char Title_h_good_events[100];  char Name_h_good_events[100];

////////////////////////////////////////////

sprintf(Title_h_count,"SFT Efficiency - Run %d  |  ADC > %d , %d < TDC < %d", Run_Number, thr, TDCthr1, TDCthr2); 
sprintf(Name_h_count,"SFT/Helicity Counter"); 
h_count = new TH1D(Name_h_count,Title_h_count,7,0,7);

sprintf(Title_h_count2,"Run %d", Run_Number); 
sprintf(Name_h_count2,"SFT/Helicity Counter"); 
h_count2 = new TH1D(Name_h_count2,Title_h_count2,7,0,7);

sprintf(Title_h_count3,"Run %d", Run_Number); 
sprintf(Name_h_count3,"SFT/Helicity Counter"); 
h_count3 = new TH1D(Name_h_count3,Title_h_count3,7,0,7);

sprintf(Title_h_count4,"Run %d", Run_Number); 
sprintf(Name_h_count4,"SFT/Helicity Counter"); 
h_count4 = new TH1D(Name_h_count4,Title_h_count4,7,0,7);

sprintf(Title_h_count5,"Run %d", Run_Number); 
sprintf(Name_h_count5,"SFT/Helicity Counter"); 
h_count5 = new TH1D(Name_h_count5,Title_h_count5,7,0,7);

////

sprintf(Title_h_count6,"SFT Efficiency - Run %d  |  ADC > %d , %d < TDC < %d", Run_Number, thr, TDCthr1, TDCthr2); 
sprintf(Name_h_count6,"SFT/Helicity Counter"); 
h_count6 = new TH1D(Name_h_count6,Title_h_count6,7,0,7);

sprintf(Title_h_count7,"Run %d", Run_Number); 
sprintf(Name_h_count7,"SFT/Helicity Counter"); 
h_count7 = new TH1D(Name_h_count7,Title_h_count7,7,0,7);

sprintf(Title_h_count8,"Run %d", Run_Number); 
sprintf(Name_h_count8,"SFT/Helicity Counter"); 
h_count8 = new TH1D(Name_h_count8,Title_h_count8,7,0,7);

sprintf(Title_h_count9,"Run %d", Run_Number); 
sprintf(Name_h_count9,"SFT/Helicity Counter"); 
h_count9 = new TH1D(Name_h_count9,Title_h_count9,7,0,7);

sprintf(Title_h_count10,"Run %d", Run_Number); 
sprintf(Name_h_count10,"SFT/Helicity Counter"); 
h_count10 = new TH1D(Name_h_count10,Title_h_count10,7,0,7);

////

sprintf(Title_h_count11,"SFT Efficiency - Run %d  |  ADC > %d , %d < TDC < %d", Run_Number, thr, TDCthr1, TDCthr2); 
sprintf(Name_h_count11,"SFT/Helicity Counter"); 
h_count11 = new TH1D(Name_h_count11,Title_h_count11,7,0,7);

sprintf(Title_h_count12,"Run %d", Run_Number); 
sprintf(Name_h_count12,"SFT/Helicity Counter"); 
h_count12 = new TH1D(Name_h_count12,Title_h_count12,7,0,7);

sprintf(Title_h_count13,"Run %d", Run_Number); 
sprintf(Name_h_count13,"SFT/Helicity Counter"); 
h_count13 = new TH1D(Name_h_count13,Title_h_count13,7,0,7);

sprintf(Title_h_count14,"Run %d", Run_Number); 
sprintf(Name_h_count14,"SFT/Helicity Counter"); 
h_count14 = new TH1D(Name_h_count14,Title_h_count14,7,0,7);

sprintf(Title_h_count15,"Run %d", Run_Number); 
sprintf(Name_h_count15,"SFT/Helicity Counter"); 
h_count15 = new TH1D(Name_h_count15,Title_h_count15,7,0,7);

////

sprintf(Title_h_count16,"SFT Efficiency - Run %d  |  ADC > %d , %d < TDC < %d", Run_Number, thr, 0, TDCthr2); 
sprintf(Name_h_count16,"SFT/Helicity Counter"); 
h_count16 = new TH1D(Name_h_count16,Title_h_count16,7,0,7);

sprintf(Title_h_count17,"Run %d", Run_Number); 
sprintf(Name_h_count17,"SFT/Helicity Counter"); 
h_count17 = new TH1D(Name_h_count17,Title_h_count17,7,0,7);

sprintf(Title_h_count18,"Run %d", Run_Number); 
sprintf(Name_h_count18,"SFT/Helicity Counter"); 
h_count18 = new TH1D(Name_h_count18,Title_h_count18,7,0,7);

sprintf(Title_h_count19,"Run %d", Run_Number); 
sprintf(Name_h_count19,"SFT/Helicity Counter"); 
h_count19 = new TH1D(Name_h_count19,Title_h_count19,7,0,7);

sprintf(Title_h_count20,"Run %d", Run_Number); 
sprintf(Name_h_count20,"SFT/Helicity Counter"); 
h_count20 = new TH1D(Name_h_count20,Title_h_count20,7,0,7);

sprintf(Title_h_good_events,"SFT Efficiency - Run %d  |  ADC > %d , %d < TDC < %d", Run_Number, thr, TDCthr1, TDCthr2); 
sprintf(Name_h_good_events,"SFT/Helicity Counter"); 
h_good_events = new TH1D(Name_h_good_events,Title_h_good_events,7,0,7);

////

TH1D *hlayer1_2;	char Title_hlayer1_2[100];	char Name_hlayer1_2[100];
TH1D *hlayer3_4;	char Title_hlayer3_4[100];	char Name_hlayer3_4[100];

sprintf(Title_hlayer1_2,"SFT Layer Distances  |  Run %d - Layers 1 & 2", Run_Number); 
sprintf(Name_hlayer1_2,"Layer 1"); 
hlayer1_2 = new TH1D(Name_hlayer1_2,Title_hlayer1_2,15,0,15);

sprintf(Title_hlayer3_4,"SFT Layer Distances  |  Run %d - Layers 3 & 4", Run_Number); 
sprintf(Name_hlayer3_4,"Layer 2"); 
hlayer3_4 = new TH1D(Name_hlayer3_4,Title_hlayer3_4,17,0,17);


////////////////////////////////////////////////////////////

if (nevents!=0) nentries=nevents;
for(Int_t i=0; i<nentries; i++){
	fChain->GetEntry(i);	//	fChain_SFT->GetEntry(i);

	counter1 = 0;
	counter2 = 0;
	counter3 = 0;
	counter4 = 0;
	counter5 = 0;
	counter6 = 0;
	counter7 = 0;
	counter8 = 0;
	helicity1or = 0;
	helicity2or = 0;
	helicity1and = 0;
	helicity2and = 0;

	counter1a = 0;
	counter2a = 0;
	counter3a = 0;
	counter4a = 0;
	counter5a = 0;
	counter6a = 0;
	counter7a = 0;
	counter8a = 0;
	helicity1ora = 0;
	helicity2ora = 0;
	helicity1anda = 0;
	helicity2anda = 0;

	counterp1 = 0;
	counterp2 = 0;
	counterp3 = 0;
	counterp4 = 0;

	p1_p2 = 0;
	p3_p4 = 0;

	Int_t count_target=0;
	
	if (Enable_TARGET_Eff!=0) {
		if(i%2000==1)	cout<<"**** "<<i<<" events done"<<endl;
	}
	else {
		if(i%10000==1)	cout<<"**** "<<i<<" events done"<<endl;
	}

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
		ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-ADC_cut_TARGET2;
		TDC_LE_TARGET[j_TARGET]=tdc_le_target[j_TARGET][0];
		TDC_TE_TARGET[j_TARGET]=tdc_te_target[j_TARGET][0];
  	}

	for (Int_t j_TOF1=0; j_TOF1<24; j_TOF1++) {
		ADC_TOF1[j_TOF1] = ADC_TOF1[j_TOF1]-par_temp_TOF1[1][j_TOF1];
	}

	for (Int_t j_TOF2=0; j_TOF2<56; j_TOF2++) {
		ADC_TOF2[j_TOF2] = ADC_TOF2[j_TOF2]-par_temp_TOF2[1][j_TOF2];
	}

	for (Int_t j_MWPC=0; j_MWPC<512; j_MWPC++) {
		MWPCADC[j_MWPC] = MwpcADC[j_MWPC]-MWPC_thr;
	}

	for (Int_t j_TDCTOF=0; j_TDCTOF<12; j_TDCTOF++) {
	    TDC_TOF1U[j_TDCTOF] = TDC_tof1U[j_TDCTOF];
	    TDC_TOF1D[j_TDCTOF] = TDC_tof1D[j_TDCTOF];
	    TDC_TOF2AO[j_TDCTOF] = TDC_tof2AO[j_TDCTOF];
	    TDC_TOF2BO[j_TDCTOF] = TDC_tof2BO[j_TDCTOF];
	    TDC_TOF2AI[j_TDCTOF] = TDC_tof2AI[j_TDCTOF];
	    TDC_TOF2BI[j_TDCTOF] = TDC_tof2BI[j_TDCTOF];
  	}

	Int_t has_data_counter = 0;
	Int_t has_data_ADCA = 0;
	h_target_ADC->Reset();
	h_target_ADCA->Reset();
	double xcoord = 0;
	Int_t unique_x = 0;
	int has_data_TOF = 0;

	int max_index = 0;
	int max_index2 = 0;
	int max_index3 = 0;
	int max_index4 = 0;

	int max_ADC = 0;
	int max_ADC2 = 0;
	int max_ADC3 = 0;
	int max_ADC4 = 0;

	bool has_TDC_hit[256] = {false};

	for(Int_t p=0; p<256; p++){
		for (Int_t k=0; k<4; k++) {
			if ((tdc_le_target[p][k]>=TDC_min_TARGET) && (tdc_le_target[p][k]<=TDC_max_TARGET)) has_TDC_hit[p] = true;
		}
	}	

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

	bool has_ADC_TOF1_hit[12] = {false};
	bool has_TDC_TOF1_hit[12] = {false};
	bool has_ADC_TOF2_hit[12] = {false};
	bool has_TDC_TOF2_hit[12] = {false};

	for (int kk=0; kk<12; kk++) {
		if (kk == 6) {
			if ((ADC_TOF2[6]>0) || (ADC_TOF2[30]>0) || (ADC_TOF2[55]>0) || (ADC_TOF2[42]>0)) {has_ADC_TOF2_hit[6]=true;}
			if (((TDC_TOF2AO[6]>TDC_TOF2_min && TDC_TOF2AO[6] < TDC_TOF2_max) || (TDC_TOF2AI[6]>TDC_TOF2_min && TDC_TOF2AI[6] < TDC_TOF2_max)) 
		    || ((TDC_TOF2BO[6]>TDC_TOF2_min && TDC_TOF2BO[6] < TDC_TOF2_max) || (TDC_TOF2BI[6]>TDC_TOF2_min && TDC_TOF2BI[6] < TDC_TOF2_max))) {has_TDC_TOF2_hit[6]=true;}

			if (ADC_TOF1[6]>0 || ADC_TOF1[18]>0) {has_ADC_TOF1_hit[6] = true;}
		  	if ((TDC_TOF1U[6]>TDC_TOF1_min && TDC_TOF1U[6]<TDC_TOF1_max) || (TDC_TOF1D[6]>TDC_TOF1_min && TDC_TOF1D[6]<TDC_TOF1_max)) {has_TDC_TOF1_hit[6] = true;}
	  	}

		else {
			if ((ADC_TOF2[kk]>0) || (ADC_TOF2[kk+24]>0) || (ADC_TOF2[kk+12]>0) || (ADC_TOF2[kk+36]>0)) {has_ADC_TOF2_hit[kk]=true;}
			if (((TDC_TOF2AO[kk]>TDC_TOF2_min && TDC_TOF2AO[kk] < TDC_TOF2_max) || (TDC_TOF2AI[kk]>TDC_TOF2_min && TDC_TOF2AI[kk] < TDC_TOF2_max)) 
		    || ((TDC_TOF2BO[kk]>TDC_TOF2_min && TDC_TOF2BO[kk] < TDC_TOF2_max) || (TDC_TOF2BI[kk]>TDC_TOF2_min && TDC_TOF2BI[kk] < TDC_TOF2_max))) {has_TDC_TOF2_hit[kk]=true;}

			if (ADC_TOF1[kk]>0 || ADC_TOF1[kk+12]>0) {has_ADC_TOF1_hit[kk] = true;}
		  	if ((TDC_TOF1U[kk]>TDC_TOF1_min && TDC_TOF1U[kk]<TDC_TOF1_max) || (TDC_TOF1D[kk]>TDC_TOF1_min && TDC_TOF1D[kk]<TDC_TOF1_max)) {has_TDC_TOF1_hit[kk] = true;}
		}
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

	bool good_MWPC = false;

	if ((C2X_hit >=1) && (C2Y_hit >=1) && (C3X_hit >=1) && 
			(C3Y_hit >=1) && (C4X_hit >=1) && (C4Y_hit >=1)) good_MWPC = true;

	int gap_hit[12] = {0};
	  int ADC_TOF1_hit[12] = {0};
	  int ADCTDC_TOF1_hit[12] = {0};
	  int ADC_TOF2_hit[12] = {0};
	  int ADCTDC_TOF2_hit[12] = {0};

	  //// Determine which gaps have ADC & TDC TARGET hits

	  for (int k=0; k<8; k++) {
	    if (ADC_High_TARGET[channel0[k]] > 0 && has_TDC_hit[channel0[k]]) gap_hit[0] = 1;
	    if (ADC_High_TARGET[channel1[k]] > 0 && has_TDC_hit[channel1[k]]) gap_hit[1] = 1;
	    if (ADC_High_TARGET[channel3[k]] > 0 && has_TDC_hit[channel3[k]]) gap_hit[3] = 1;
	    if (ADC_High_TARGET[channel4[k]] > 0 && has_TDC_hit[channel4[k]]) gap_hit[4] = 1;
	    if (ADC_High_TARGET[channel6[k]] > 0 && has_TDC_hit[channel6[k]]) gap_hit[6] = 1;
	    if (ADC_High_TARGET[channel7[k]] > 0 && has_TDC_hit[channel7[k]]) gap_hit[7] = 1;
	    if (ADC_High_TARGET[channel9[k]] > 0 && has_TDC_hit[channel9[k]]) gap_hit[9] = 1;
	    if (ADC_High_TARGET[channel10[k]] > 0 && has_TDC_hit[channel10[k]]) gap_hit[10] = 1;
	  }

	  for (int k=0; k<8; k++) {
	    if (ADC_High_TARGET[channel2[k]] > 0 && has_TDC_hit[channel2[k]]) gap_hit[2] = 1;
	    if (ADC_High_TARGET[channel5[k]] > 0 && has_TDC_hit[channel5[k]]) gap_hit[5] = 1;
	    if (ADC_High_TARGET[channel8[k]] > 0 && has_TDC_hit[channel8[k]]) gap_hit[8] = 1;
	    if (ADC_High_TARGET[channel11[k]] > 0 && has_TDC_hit[channel11[k]]) gap_hit[11] = 1;
	  }

	  for (int k=0; k<12; k++) {
	      if (has_ADC_TOF1_hit[k]) {
	        if (has_TDC_TOF1_hit[k]) {ADCTDC_TOF1_hit[k]++;}
	        else {ADC_TOF1_hit[k]++;}
	      }
	      if (has_ADC_TOF2_hit[k]) {
	        if (has_TDC_TOF2_hit[k]) {ADCTDC_TOF2_hit[k]++;}
	        else {ADC_TOF2_hit[k]++;}
	      }
	  }

	  int gap_counter[12] = {0};

	  //// Select gap for track fitting

	  for (int k=0; k<12; k++) {
	    if ((has_ADC_TOF2_hit[k] && has_TDC_TOF2_hit[k]) && (has_ADC_TOF1_hit[k] && has_TDC_TOF1_hit[k])) gap_counter[k] = gap_counter[k] + 5;
	    if ((has_ADC_TOF2_hit[k] || (has_ADC_TOF2_hit[k] && has_TDC_TOF2_hit[k])) && (has_ADC_TOF1_hit[k] || (has_ADC_TOF1_hit[k] && has_TDC_TOF1_hit[k]))) gap_counter[k] = gap_counter[k] + 2;
	    if (gap_hit[k] > 0) gap_counter[k] = gap_counter[k] + 4;
	    if (ADCTDC_TOF1_hit[k] > 0) gap_counter[k] = gap_counter[k] + 2;
	    if (ADCTDC_TOF2_hit[k] > 0) gap_counter[k] = gap_counter[k] + 2;
	    if (ADC_TOF2_hit[k] > 0) gap_counter[k]++;
	    if (ADC_TOF1_hit[k] > 0) gap_counter[k]++;
	  }

	  /*
	  for (int k=0; k<12; k++) {
	    if (gap_hit[k] > 0) gap_counter[k]++;
	    if (ADC_TOF1_hit[k] > 0) gap_counter[k]++;
	    if (ADCTDC_TOF1_hit[k] > 0) gap_counter[k]++;
	    if (ADC_TOF2_hit[k] > 0) gap_counter[k]++;
	    if (ADCTDC_TOF2_hit[k] > 0) gap_counter[k]++;
	  }
	  */

	  int max_gap[12] = {0};
	  int high_gap_hit = 0;
	  int gap_count[12] = {0};

	  for (int k=0; k<12; k++) {
	    if (gap_counter[k] >= high_gap_hit) {
	      high_gap_hit = gap_counter[k];
	      max_gap[k] = max_gap[k]+gap_counter[k];
	      }
	  }

	///// Select edge bar for track fitting

	int gap_to_fit;

	  for (int k=0; k<12; k++) {
	    if (max_gap[k]>=high_gap_hit) gap_to_fit = k;
	  }

	  double Xloc_gap = 0;
	  double Yloc_gap = 0;

	  if (gap_to_fit == 0) {
	    for (int j=0; j<8; j++) {
	      if ((ADC_High_TARGET[channel0[j]]>0) && (has_TDC_hit[channel0[j]])) {
	        Xloc_gap = Xloc[channel0[j]];
	        Yloc_gap = Yloc[channel0[j]];
	      }
	    }
	    if (Xloc_gap == 0) Xloc_gap = Xloc[channel0[3]];
	    if (Yloc_gap == 0) Yloc_gap = Yloc[channel0[3]];
	  }
	   if (gap_to_fit == 1) {
	    for (int j=0; j<8; j++) {
	      if ((ADC_High_TARGET[channel1[j]]>0) && (has_TDC_hit[channel1[j]])) {
	        Xloc_gap = Xloc[channel1[j]];
	        Yloc_gap = Yloc[channel1[j]];
	      }
	    }
	    if (Xloc_gap == 0) Xloc_gap = Xloc[channel1[3]];
	    if (Yloc_gap == 0) Yloc_gap = Yloc[channel1[3]];
	  }
	   if (gap_to_fit == 2) {
	    for (int j=0; j<8; j++) {
	      if ((ADC_High_TARGET[channel2[j]]>0) && (has_TDC_hit[channel2[j]])) {
	        Xloc_gap = Xloc[channel2[j]];
	        Yloc_gap = Yloc[channel2[j]];
	      }
	    }
	    if (Xloc_gap == 0) Xloc_gap = Xloc[channel2[3]];
	    if (Yloc_gap == 0) Yloc_gap = Yloc[channel2[3]];
	  }
	   if (gap_to_fit == 3) {
	    for (int j=0; j<8; j++) {
	      if ((ADC_High_TARGET[channel3[j]]>0) && (has_TDC_hit[channel3[j]])) {
	        Xloc_gap = Xloc[channel3[j]];
	        Yloc_gap = Yloc[channel3[j]];
	      }
	    }
	    if (Xloc_gap == 0) Xloc_gap = Xloc[channel3[3]];
	    if (Yloc_gap == 0) Yloc_gap = Yloc[channel3[3]];
	  }
	   if (gap_to_fit == 4) {
	    for (int j=0; j<8; j++) {
	      if ((ADC_High_TARGET[channel4[j]]>0) && (has_TDC_hit[channel4[j]])) {
	        Xloc_gap = Xloc[channel4[j]];
	        Yloc_gap = Yloc[channel4[j]];
	      }
	    }
	    if (Xloc_gap == 0) Xloc_gap = Xloc[channel4[3]];
	    if (Yloc_gap == 0) Yloc_gap = Yloc[channel4[3]];
	  }
	   if (gap_to_fit == 5) {
	    for (int j=0; j<8; j++) {
	      if ((ADC_High_TARGET[channel5[j]]>0) && (has_TDC_hit[channel5[j]])) {
	        Xloc_gap = Xloc[channel5[j]];
	        Yloc_gap = Yloc[channel5[j]];
	      }
	    }
	    if (Xloc_gap == 0) Xloc_gap = Xloc[channel5[3]];
	    if (Yloc_gap == 0) Yloc_gap = Yloc[channel5[3]];
	  }
	   if (gap_to_fit == 6) {
	    for (int j=0; j<8; j++) {
	      if ((ADC_High_TARGET[channel6[j]]>0) && (has_TDC_hit[channel6[j]])) {
	        Xloc_gap = Xloc[channel6[j]];
	        Yloc_gap = Yloc[channel6[j]];
	      }
	    }
	    if (Xloc_gap == 0) Xloc_gap = Xloc[channel6[3]];
	    if (Yloc_gap == 0) Yloc_gap = Yloc[channel6[3]];
	  }
	   if (gap_to_fit == 7) {
	    for (int j=0; j<8; j++) {
	      if ((ADC_High_TARGET[channel7[j]]>0) && (has_TDC_hit[channel7[j]])) {
	        Xloc_gap = Xloc[channel7[j]];
	        Yloc_gap = Yloc[channel7[j]];
	      }
	    }
	    if (Xloc_gap == 0) Xloc_gap = Xloc[channel7[3]];
	    if (Yloc_gap == 0) Yloc_gap = Yloc[channel7[3]];
	  }
	   if (gap_to_fit == 8) {
	    for (int j=0; j<8; j++) {
	      if ((ADC_High_TARGET[channel8[j]]>0) && (has_TDC_hit[channel8[j]])) {
	        Xloc_gap = Xloc[channel8[j]];
	        Yloc_gap = Yloc[channel8[j]];
	      }
	    }
	    if (Xloc_gap == 0) Xloc_gap = Xloc[channel8[3]];
	    if (Yloc_gap == 0) Yloc_gap = Yloc[channel8[3]];
	  }
	   if (gap_to_fit == 9) {
	    for (int j=0; j<8; j++) {
	      if ((ADC_High_TARGET[channel9[j]]>0) && (has_TDC_hit[channel9[j]])) {
	        Xloc_gap = Xloc[channel9[j]];
	        Yloc_gap = Yloc[channel9[j]];
	      }
	    }
	    if (Xloc_gap == 0) Xloc_gap = Xloc[channel9[3]];
	    if (Yloc_gap == 0) Yloc_gap = Yloc[channel9[3]];
	  }
	   if (gap_to_fit == 10) {
	    for (int j=0; j<8; j++) {
	      if ((ADC_High_TARGET[channel10[j]]>0) && (has_TDC_hit[channel10[j]])) {
	        Xloc_gap = Xloc[channel10[j]];
	        Yloc_gap = Yloc[channel10[j]];
	      }
	    }
	    if (Xloc_gap == 0) Xloc_gap = Xloc[channel10[3]];
	    if (Yloc_gap == 0) Yloc_gap = Yloc[channel10[3]];
	  }
	   if (gap_to_fit == 11) {
	    for (int j=0; j<8; j++) {
	      if ((ADC_High_TARGET[channel11[j]]>0) && (has_TDC_hit[channel11[j]])) {
	        Xloc_gap = Xloc[channel11[j]];
	        Yloc_gap = Yloc[channel11[j]];
	      }
	    }
	    if (Xloc_gap == 0) Xloc_gap = Xloc[channel11[3]];
	    if (Yloc_gap == 0) Yloc_gap = Yloc[channel11[3]];
	  }

	/// Determine which gap point is closest to edge bar

	double xdistance1 = pow((TOF_Xloc[(gap_to_fit*3)]-Xloc_gap),2);
	double ydistance1 = pow((TOF_Yloc[(gap_to_fit*3)]-Yloc_gap),2);

	double xdistance2 = pow((TOF_Xloc[(gap_to_fit*3)+1]-Xloc_gap),2);
	double ydistance2 = pow((TOF_Yloc[(gap_to_fit*3)+1]-Yloc_gap),2);

	double xdistance3 = pow((TOF_Xloc[(gap_to_fit*3)+2]-Xloc_gap),2);
	double ydistance3 = pow((TOF_Yloc[(gap_to_fit*3)+2]-Yloc_gap),2);

	double xhyp1 = double(sqrt(double(xdistance1) + double(ydistance1)));
	double xhyp2 = double(sqrt(double(xdistance2) + double(ydistance2)));
	double xhyp3 = double(sqrt(double(xdistance3) + double(ydistance3)));

	int closest_gap_point;

	if ((xhyp1 <= xhyp2) && (xhyp1 <= xhyp3)) {
		closest_gap_point = 0;
 	}
	else if (xhyp2 <= xhyp3) {
		closest_gap_point = 1;
	}
	else {
		closest_gap_point = 2;
	}

	double x_cent = Xloc[max_index];
	double y_cent = Yloc[max_index];
	double hyp[256] = {-1};

	for(Int_t j=0; j<256; j++){
		if(ADC_High_TARGET[j]>=Angle_ADC_cut) {
			hyp[j] = sqrt(pow(x_cent - Xloc[j],2) + pow(y_cent - Yloc[j],2));
		}
	}

	h_target_ADC->Fill(Xloc[max_index],Yloc[max_index]);
	h_target_ADC->Fill(TOF_Xloc[(gap_to_fit*3) + closest_gap_point],TOF_Yloc[(gap_to_fit*3) + closest_gap_point]);

	h_target_ADCA->Fill(Xloc[max_index],Yloc[max_index]);

	/// Fill TARGET Histograms

	for(Int_t p=0; p<256; p++){

		if ((p == max_index) || (p == max_index2) || (p == max_index3) || (p == max_index4)) continue;

		if ((has_TDC_hit[p]) && (ADC_High_TARGET[p]>=Angle_ADC_cut) && (adc_low_target[p]-ADC_cut_TARGET2>0)) {
			h_target_ADC->Fill(Xloc[p],Yloc[p]);
			has_data_counter++;
			if (xcoord != Xloc[p]) {
				xcoord = Xloc[p];
				unique_x++;
			}
		}
	}					

	if ((has_data_counter >= 5) && (unique_x > 1) && (has_data_TOF > 0) && (good_MWPC)) {

//		cout << "Calculating efficiency: Event " << i <<endl;
//		cout << i << endl;
		event_counter++;
		
		TF1 *fit_line_ADC;			
		h_target_ADC->Fit("pol1","Q0CM");
		fit_line_ADC = h_target_ADC->GetFunction("pol1");

		 /// Calculate distance of each point from fit_line_ADC

		double distances1[256] = {0};
		double x_distance1[256] = {0};
		double y_distance1[256] = {0};

		for (int q=0; q<256; q++) {
			x_distance1[q] = pow((fit_line_ADC->Eval(Xloc[q],0,0) - Yloc[q]),2);
			y_distance1[q] = pow((fit_line_ADC->GetX(Yloc[q],-30,30) - Xloc[q]),2);
	
			distances1[q] = double(sqrt(double(x_distance1[q]) + double(y_distance1[q])));
		}

		/// Remove outliers from track fit

		for(int j=0; j<256; j++) {
			if ((j == max_index) || (j == max_index2) || (j == max_index3) || (j == max_index4)) continue;
			if (((distances1[j] < path_margin) || (fabs((fit_line_ADC->Eval(Xloc[j],0,0) - Yloc[j])) < height_margin)) && (has_TDC_hit[j]) && (adc_high_target[j]-ADC_cut_TARGET3>=Angle_ADC_cut) /* || (fabs((fit_line_ADC3->GetX(Yloc[j],-30,30) - Xloc[j])) < width_margin)*/){
				h_target_ADCA->Fill(Xloc[j],Yloc[j]);
				has_data_ADCA++;
//				cout << j << " - " << distances[j] << endl;
			}		
		}

		if (has_data_ADCA > 2) {

		    h_target_ADCA->Fill(TOF_Xloc[(gap_to_fit*3) + closest_gap_point],TOF_Yloc[(gap_to_fit*3) + closest_gap_point]);

		    /// Refit histogram with filtered data + closest gap point
		
			TF1 *fit_line_ADCA;
			h_target_ADCA->Fit("pol1", "Q0CM");
			fit_line_ADCA = h_target_ADCA->GetFunction("pol1");
			
			double distances[256] = {0};
			double x_distance[256] = {0};
			double y_distance[256] = {0};
			
			for (int q=0; q<256; q++) {
				x_distance[q] = pow((fit_line_ADCA->Eval(Xloc[q],0,0) - Yloc[q]),2);
				y_distance[q] = pow((fit_line_ADCA->GetX(Yloc[q],-30,30) - Xloc[q]),2);
	
				distances[q] = double(sqrt(double(x_distance[q]) + double(y_distance[q])));
			}

			for(Int_t m=0; m<256; m++){

				if ((m == max_index) || (m == max_index2) || (m == max_index3) || (m == max_index4)) continue;

				if (((distances[m] < path_margin) || (fabs((fit_line_ADC->Eval(Xloc[m],0,0) - Yloc[m])) < height_margin) || (fabs((fit_line_ADC->GetX(Yloc[m],-30,30) - Xloc[m])) < width_margin)) /*&& (hyp[m] > 9)*/) {

					if ((m == 0) || (m == 6) || (m == 16) || (m == 28) || (m == 42) || (m == 58) || (m == 74) || (m == 92) || (m == 110) || (m == 128) || (m == 146) || (m == 164) || (m == 182) || (m == 198) || (m == 214) || (m == 228) || (m == 240) || (m == 250)) {
						if ((ADC_High_TARGET[m+1] > 0)) {
							bar_adj_counter[m]++;
							if (ADC_High_TARGET[m] > 0) {
								bar_hit_counter[m]++;
							}
						}
					}

					else if ((m == 5) || (m == 15) || (m == 27) || (m == 41) || (m == 57) || (m == 73) || (m == 91) || (m == 109) || (m == 127) || (m == 145) || (m == 163) || (m == 181) || (m == 197) || (m == 213) || (m == 227) || (m == 239) || (m == 249) || (m == 255)) {
						if ((ADC_High_TARGET[m-1] > 0)) {
							bar_adj_counter[m]++;
							if (ADC_High_TARGET[m] > 0) {
								bar_hit_counter[m]++;
							}
						}
					}

					else {
						if ((ADC_High_TARGET[m-1] > 0) || (ADC_High_TARGET[m+1] > 0)) {
							bar_adj_counter[m]++;
							if (ADC_High_TARGET[m] > 0) {
								bar_hit_counter[m]++;
							}
						}
					}
				} 
			}
		}
	}

	if ((has_data_counter >= 5) && (unique_x == 1) && (has_data_TOF > 0) && (good_MWPC)) {
//	cout << "Vertical line track" << endl;
//	cout << i << endl;
	event_counter++;

		for(Int_t m=0; m<256; m++) {

			if ((m == max_index2) || (m == max_index3) || (m == max_index4)) continue;

			if (Xloc[m] == xcoord) {
				
				bar_adj_counter[m]++;

				if ((m == 0) || (m == 6) || (m == 16) || (m == 28) || (m == 42) || (m == 58) || (m == 74) || (m == 92) || (m == 110) || (m == 128) || (m == 146) || (m == 164) || (m == 182) || (m == 198) || (m == 214) || (m == 228) || (m == 240) || (m == 250)) {
					if ((ADC_High_TARGET[m] > 0) || (ADC_High_TARGET[m+1] > 0)) {
//						bar_adj_counter[m]++;
//						if (ADC_High_TARGET[m] > 0) {
							bar_hit_counter[m]++;
//						}
					}
				}

				else if ((m == 5) || (m == 15) || (m == 27) || (m == 41) || (m == 57) || (m == 73) || (m == 91) || (m == 109) || (m == 127) || (m == 145) || (m == 163) || (m == 181) || (m == 197) || (m == 213) || (m == 227) || (m == 239) || (m == 249) || (m == 255)) {
					if ((ADC_High_TARGET[m] > 0) || (ADC_High_TARGET[m-1] > 0)) {
//						bar_adj_counter[m]++;
//						if (ADC_High_TARGET[m] > 0) {
							bar_hit_counter[m]++;
//						}
					}
				}

				else {
					if ((ADC_High_TARGET[m] > 0) || (ADC_High_TARGET[m+1] > 0) || (ADC_High_TARGET[m-1] > 0)) {
//						bar_adj_counter[m]++;
//						if (ADC_High_TARGET[m] > 0) {
							bar_hit_counter[m]++;
//						}
					}
				}

			}
		}
	}

	//SFT Efficiency
	if ((has_data_counter > 3) && (has_data_TOF > 0) && (good_MWPC)) {
		good_event++;

	for(Int_t ii=0; ii<128; ii++){
		ADC_High_sft[ii] = ADC_High_sft[ii] - par_temp_SFT[1][ii];
		fiber[ii]=ADC_High_sft[par_temp[1][ii]]-par_temp_SFT[1][ii];
		if(TDC_LE_sft[par_temp[1][ii]][0]>TDC_min_SFT && TDC_LE_sft[par_temp[1][ii]][0]<TDC_max_SFT)
			fiberTDC[ii]=TDC_LE_sft[par_temp[1][ii]][0];
		else if(TDC_LE_sft[par_temp[1][ii]][1]>TDC_min_SFT && TDC_LE_sft[par_temp[1][ii]][1]<TDC_max_SFT)
			fiberTDC[ii]=TDC_LE_sft[par_temp[1][ii]][1];
		else if(TDC_LE_sft[par_temp[1][ii]][2]>TDC_min_SFT && TDC_LE_sft[par_temp[1][ii]][2]<TDC_max_SFT)
			fiberTDC[ii]=TDC_LE_sft[par_temp[1][ii]][2];
		else if(TDC_LE_sft[par_temp[1][ii]][3]>TDC_min_SFT && TDC_LE_sft[par_temp[1][ii]][3]<TDC_max_SFT)
			fiberTDC[ii]=TDC_LE_sft[par_temp[1][ii]][3];
		else if(TDC_LE_sft[par_temp[1][ii]][4]>TDC_min_SFT && TDC_LE_sft[par_temp[1][ii]][4]<TDC_max_SFT)
			fiberTDC[ii]=TDC_LE_sft[par_temp[1][ii]][4];
		else if(TDC_LE_sft[par_temp[1][ii]][5]>TDC_min_SFT && TDC_LE_sft[par_temp[1][ii]][5]<TDC_max_SFT)
			fiberTDC[ii]=TDC_LE_sft[par_temp[1][ii]][5];
		else
			fiberTDC[ii]=TDC_LE_sft[par_temp[1][ii]][6];
	}

	///Layer 1
		for (Int_t w=0; w<15; w++) {
			if (((fiber[w]>0) && (fiberTDC[w]>TDCthr1) && (fiberTDC[w]<TDCthr2)) || ((fiberTDC[w+64]>TDCthr1) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				counter1++;
				if (counterp1 == 0) {
					counterp1 = w;
				}
			}
			if (((fiber[w]>0) && (fiberTDC[w]>TDCthr1) && (fiberTDC[w]<TDCthr2)) && ((fiberTDC[w+64]>TDCthr1) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				counter5++;
			}
		}

	///Layer 2
		for (Int_t w=15; w<30; w++) {
			if (((fiber[w]>0) && (fiberTDC[w]>TDCthr1) && (fiberTDC[w]<TDCthr2)) || ((fiberTDC[w+64]>TDCthr1) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				counter2++;
				if (counterp2 == 0) {
					counterp2 = w;
				}
			}
			if (((fiber[w]>0) && (fiberTDC[w]>TDCthr1) && (fiberTDC[w]<TDCthr2)) && ((fiberTDC[w+64]>TDCthr1) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				counter6++;
			}
		}

	///Layer 3
		for (Int_t w=30; w<47; w++) {
			if (((fiber[w]>0) && (fiberTDC[w]>TDCthr1) && (fiberTDC[w]<TDCthr2)) || ((fiberTDC[w+64]>TDCthr1) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				counter3++;
				if (counterp3 == 0) {
					counterp3 = w;
				}
			}
			if (((fiber[w]>0) && (fiberTDC[w]>TDCthr1) && (fiberTDC[w]<TDCthr2)) && ((fiberTDC[w+64]>TDCthr1) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				counter7++;
			}
		}

	///Layer 4
		for (Int_t w=47; w<64; w++) {
			if (((fiber[w]>0) && (fiberTDC[w]>TDCthr1) && (fiberTDC[w]<TDCthr2)) || ((fiberTDC[w+64]>TDCthr1) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				counter4++;
				if (counterp4 == 0) {
					counterp4 = w;
				}
			}
			if (((fiber[w]>0) && (fiberTDC[w]>TDCthr1) && (fiberTDC[w]<TDCthr2)) && ((fiberTDC[w+64]>TDCthr1) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				counter8++;
			}
		}

///////////////////////////////////////////////////////////////////////// NO TDC CUT Efficiency

		for (Int_t w=0; w<15; w++) {
			if (((fiber[w]>0) && (fiberTDC[w]>0) && (fiberTDC[w]<TDCthr2)) || ((fiberTDC[w+64]>0) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				counter1a++;
			}
			if (((fiber[w]>0) && (fiberTDC[w]>0) && (fiberTDC[w]<TDCthr2)) && ((fiberTDC[w+64]>0) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				counter5a++;
			}
		}

	///Layer 2
		for (Int_t w=15; w<30; w++) {
			if (((fiber[w]>0) && (fiberTDC[w]>0) && (fiberTDC[w]<TDCthr2)) || ((fiberTDC[w+64]>0) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				counter2a++;
			}
			if (((fiber[w]>0) && (fiberTDC[w]>0) && (fiberTDC[w]<TDCthr2)) && ((fiberTDC[w+64]>0) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				counter6a++;
			}
		}

	///Layer 3
		for (Int_t w=30; w<47; w++) {
			if (((fiber[w]>0) && (fiberTDC[w]>0) && (fiberTDC[w]<TDCthr2)) || ((fiberTDC[w+64]>0) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				counter3a++;
			}
			if (((fiber[w]>0) && (fiberTDC[w]>0) && (fiberTDC[w]<TDCthr2)) && ((fiberTDC[w+64]>0) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				counter7a++;
			}
		}

	///Layer 4
		for (Int_t w=47; w<64; w++) {
			if (((fiber[w]>0) && (fiberTDC[w]>0) && (fiberTDC[w]<TDCthr2)) || ((fiberTDC[w+64]>0) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				counter4a++;
			}
			if (((fiber[w]>0) && (fiberTDC[w]>0) && (fiberTDC[w]<TDCthr2)) && ((fiberTDC[w+64]>0) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				counter8a++;
			}
		}

/////////////////////////////////////////////////////////////////////////

		for (Int_t w=0; w<30; w++) {
			if (((fiber[w]>0) && (fiberTDC[w]>TDCthr1) && (fiberTDC[w]<TDCthr2)) || ((fiberTDC[w+64]>TDCthr1) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				helicity1or++;
			}
			if (((fiber[w]>0) && (fiberTDC[w]>TDCthr1) && (fiberTDC[w]<TDCthr2)) && ((fiberTDC[w+64]>TDCthr1) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				helicity1and++;
			}
		}

		for (Int_t w=30; w<63; w++) {
			if (((fiber[w]>0) && (fiberTDC[w]>TDCthr1) && (fiberTDC[w]<TDCthr2)) || ((fiberTDC[w+64]>TDCthr1) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				helicity2or++;
			}
			if (((fiber[w]>0) && (fiberTDC[w]>TDCthr1) && (fiberTDC[w]<TDCthr2)) && ((fiberTDC[w+64]>TDCthr1) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				helicity2and++;
			}
		}

	////////////////////////////////////////////////////////////////////////NO TDC Cut Efficiency
		for (Int_t w=0; w<30; w++) {
			if (((fiber[w]>0) && (fiberTDC[w]>0) && (fiberTDC[w]<TDCthr2)) || ((fiberTDC[w+64]>0) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				helicity1ora++;
			}
			if (((fiber[w]>0) && (fiberTDC[w]>0) && (fiberTDC[w]<TDCthr2)) && ((fiberTDC[w+64]>0) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				helicity1anda++;
			}
		}

		for (Int_t w=30; w<63; w++) {
			if (((fiber[w]>0) && (fiberTDC[w]>0) && (fiberTDC[w]<TDCthr2)) || ((fiberTDC[w+64]>0) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>0))) {
				helicity2ora++;
			}
			if (((fiber[w]>0) && (fiberTDC[w]>0) && (fiberTDC[w]<TDCthr2)) && ((fiberTDC[w+64]>0) && (fiberTDC[w+64]<TDCthr2) && (fiber[w+64]>thr))) {
				helicity2anda++;
			}
		}

	// Individual Layer Counter	
		if (counter1>0) countlayer1or++;
		if (counter5>0) countlayer1and++;
		if (counter2>0) countlayer2or++;
		if (counter6>0) countlayer2and++;
		if (counter3>0) countlayer3or++;
		if (counter7>0) countlayer3and++;
		if (counter4>0) countlayer4or++;
		if (counter8>0) countlayer4and++;
		if (helicity1or>0) counthelicity1++;
		if (helicity2or>0) counthelicity2++;

	//// Layer Counter
		if ((counter1 > 0) && (counter2 > 0) && (counter3 > 0) && (counter4 > 0)) count_or++;;
		if ((counter5 > 0) && (counter6 > 0) && (counter7 > 0) && (counter8 > 0)) count_and++;
	
		if ((counter1 == 1) && (counter2 == 1) && (counter3 == 1) && (counter4 == 1)) {
			onehitcount++;;
	//		cout << "One Hit Per Layer: Event " << j << endl;
			p1_p2 = abs((counterp2-15)-counterp1);
			p3_p4 = abs((counterp4-17)-counterp3);
	//		cout << p1_p2 << endl;
	//		cout << p3_p4 << endl;
			(countlayer1_2[p1_p2])++;
			(countlayer3_4[p3_p4])++;
		}

	//// Helicity Counter
		if ((helicity2or > 0) && (helicity1or > 0)) count_hel_or++;
		if ((helicity2and > 0) && (helicity1and > 0)) count_hel_and++;;

	////////////////NO TDC CUT
		if ((counter1a > 0) && (counter2a > 0) && (counter3a > 0) && (counter4a > 0)) count_ora++;
		if ((counter5a > 0) && (counter6a > 0) && (counter7a > 0) && (counter8a > 0)) count_anda++;
		if ((helicity2ora > 0) && (helicity1ora > 0)) count_hel_ora++;
		if ((helicity2anda > 0) && (helicity1anda > 0)) count_hel_anda++;
	}
}

cout << "" << endl;
	cout << "Number of events calculated: " << good_event << " -- Percentage of total events: " << 100*(double(double(good_event)/double(nentries))) << "%" <<endl;

if (Enable_TARGET_Eff!=0) {	
	cout << "" << endl;
	for (Int_t j=0; j<256; j++) {
		if (bar_adj_counter[j] == 0) {
			Efficiency[j] = 0;
		}
		else {
			Efficiency[j] = double(double(bar_hit_counter[j])/double(bar_adj_counter[j]));
		}
		cout << "Channel "<< j << " bhit: " << bar_hit_counter[j] << " -- badj: " << bar_adj_counter[j] << " -- Eff: " << Efficiency[j] << endl;
	}

	/// TARGET Efficiency Canvas

	TH1D *TARGET_h_count;	char Title_TARGET_h_count[100];	char Name_TARGET_h_count[100];
	sprintf(Title_TARGET_h_count,"TARGET Efficiency - Run %d: ADC > %d", Run_Number, ADC_cut_TARGET); 
	sprintf(Name_TARGET_h_count,"TARGET Efficiency"); 
	TARGET_h_count = new TH1D(Name_TARGET_h_count,Title_TARGET_h_count,258,-1,256);

	for(int k=0 ; k<256 ; k++) {
		TARGET_h_count->Fill(k,Efficiency[k]);
	}

	char Name_Can[100];			char Title_Can[100];
	sprintf(Name_Can,"ADC -- Run %d -- TARGET_Efficiency_NEW.C",Run_Number);
	sprintf(Title_Can,"ADC -- Run %d  -- TARGET_Efficiency_NEW.C",Run_Number);

	TCanvas* c1_TARGET = new TCanvas(Name_Can, Title_Can, 500,300,1300,600);
	c1_TARGET->cd();
	TARGET_h_count->Draw();
}

/// SFT Efficiency Canvas

for(Int_t i=0; i<15; i++) hlayer1_2->Fill(i,countlayer1_2[i]);
for(Int_t j=0; j<17; j++) hlayer3_4->Fill(j,countlayer3_4[j]);

double eff_SFT_or = (double(count_or)/double(good_event))*100;
double eff_SFT_and = (double(count_and)/double(good_event))*100;
double eff_hel_or = (double(count_hel_or)/double(good_event))*100;
double eff_hel_and = (double(count_hel_and)/double(good_event))*100;

double eff_l1_or = (double(countlayer1or)/double(good_event))*100;
double eff_l1_and = (double(countlayer1and)/double(good_event))*100;
double eff_l2_or = (double(countlayer2or)/double(good_event))*100;
double eff_l2_and = (double(countlayer2and)/double(good_event))*100;

double eff_l3_or = (double(countlayer3or)/double(good_event))*100;
double eff_l3_and = (double(countlayer3and)/double(good_event))*100;
double eff_l4_or = (double(countlayer4or)/double(good_event))*100;
double eff_l4_and = (double(countlayer4and)/double(good_event))*100;

double eff_SFT_ora = (double(count_or)/double(nentries))*100;
double eff_SFT_anda = (double(count_and)/double(nentries))*100;
double eff_hel_ora = (double(count_hel_or)/double(nentries))*100;
double eff_hel_anda = (double(count_hel_and)/double(nentries))*100;

double scale = 0.01; 
eff_SFT_or = floor(eff_SFT_or / scale + 0.5) * scale;
eff_SFT_and = floor(eff_SFT_and / scale + 0.5) * scale;
eff_hel_or = floor(eff_hel_or / scale + 0.5) * scale;
eff_hel_and = floor(eff_hel_and / scale + 0.5) * scale;

eff_l1_or = floor(eff_l1_or / scale + 0.5) * scale;
eff_l1_and = floor( eff_l1_and / scale + 0.5) * scale;
eff_l2_or = floor(eff_l2_or / scale + 0.5) * scale;
eff_l2_and = floor(eff_l2_and / scale + 0.5) * scale;

eff_l3_or = floor(eff_l3_or / scale + 0.5) * scale;
eff_l3_and = floor( eff_l3_and / scale + 0.5) * scale;
eff_l4_or = floor(eff_l4_or / scale + 0.5) * scale;
eff_l4_and = floor(eff_l4_and / scale + 0.5) * scale;

eff_SFT_ora = floor(eff_SFT_ora / scale + 0.5) * scale;
eff_SFT_anda = floor(eff_SFT_anda / scale + 0.5) * scale;
eff_hel_ora = floor(eff_hel_ora / scale + 0.5) * scale;
eff_hel_anda = floor(eff_hel_anda / scale + 0.5) * scale;

/////////////////////////////////////////////////////

char str1[100];		char str9[100];
char str2[100];		char str10[100];
char str3[100];		char str11[100];
char str4[100];		char str12[100];
char str5[100];		char str13[100];
char str6[100];		char str14[100];
char str7[100];		char str15[100];
char str8[100];		char str16[100];

sprintf(str1,"%4.2f%%",eff_SFT_or);
sprintf(str2,"%4.2f%%",eff_SFT_and);
sprintf(str3,"%4.2f%%",eff_hel_or);
sprintf(str4,"%4.2f%%",eff_hel_and);
sprintf(str5,"%4.2f%%",eff_l1_or);
sprintf(str6,"%4.2f%%",eff_l1_and);
sprintf(str7,"%4.2f%%",eff_l2_or);
sprintf(str8,"%4.2f%%",eff_l2_and);
sprintf(str9,"%4.2f%%",eff_l3_or);
sprintf(str10,"%4.2f%%",eff_l3_and);
sprintf(str11,"%4.2f%%",eff_l4_or);
sprintf(str12,"%4.2f%%",eff_l4_and);
sprintf(str13,"%4.2f%%",eff_SFT_ora);
sprintf(str14,"%4.2f%%",eff_SFT_anda);
sprintf(str15,"%4.2f%%",eff_hel_ora);
sprintf(str16,"%4.2f%%",eff_hel_anda);

///

string SFTorEff = str1;
string SFTandEff = str2;
string HelorEff = str3;
string HelandEff = str4;

string L1orEff = str5;
string L1andEff = str6;
string L2orEff = str7;
string L2andEff = str8;

string L3orEff = str9;
string L3andEff = str10;
string L4orEff = str11;
string L4andEff = str12;

string SFTorEffa = str13;
string SFTandEffa = str14;
string HelorEffa = str15;
string HelandEffa = str16;

///
h_count->Fill(1,good_event);
h_count2->Fill(3,count_hel_or);
h_count3->Fill(3,count_hel_and);
h_count4->Fill(5,count_or);
h_count5->Fill(5,count_and);

//h_count6->Fill(1,counthelicity1);
h_count6->Fill(1,good_event);
h_count7->Fill(3,countlayer1or);
h_count8->Fill(3,countlayer1and);
h_count9->Fill(5,countlayer2or);
h_count10->Fill(5,countlayer2and);

//h_count11->Fill(1,counthelicity2);
h_count11->Fill(1,good_event);
h_count12->Fill(3,countlayer3or);
h_count13->Fill(3,countlayer3and);
h_count14->Fill(5,countlayer4or);
h_count15->Fill(5,countlayer4and);

/// No TDC Cut

h_count16->Fill(1,nentries);
h_count17->Fill(3,count_hel_or);
h_count18->Fill(3,count_hel_and);
h_count19->Fill(5,count_or);
h_count20->Fill(5,count_and);

///

h_good_events->Fill(1,good_event);

h_count->GetXaxis()->SetLabelOffset(999);
h_count->GetXaxis()->SetLabelSize(0);
h_count2->GetXaxis()->SetLabelOffset(999);
h_count2->GetXaxis()->SetLabelSize(0);
h_count3->GetXaxis()->SetLabelOffset(999);
h_count3->GetXaxis()->SetLabelSize(0);
h_count4->GetXaxis()->SetLabelOffset(999);
h_count4->GetXaxis()->SetLabelSize(0);
h_count5->GetXaxis()->SetLabelOffset(999);
h_count5->GetXaxis()->SetLabelSize(0);

h_count6->GetXaxis()->SetLabelOffset(999);
h_count6->GetXaxis()->SetLabelSize(0);
h_count7->GetXaxis()->SetLabelOffset(999);
h_count7->GetXaxis()->SetLabelSize(0);
h_count8->GetXaxis()->SetLabelOffset(999);
h_count8->GetXaxis()->SetLabelSize(0);
h_count9->GetXaxis()->SetLabelOffset(999);
h_count9->GetXaxis()->SetLabelSize(0);
h_count10->GetXaxis()->SetLabelOffset(999);
h_count10->GetXaxis()->SetLabelSize(0);

h_count11->GetXaxis()->SetLabelOffset(999);
h_count11->GetXaxis()->SetLabelSize(0);
h_count12->GetXaxis()->SetLabelOffset(999);
h_count12->GetXaxis()->SetLabelSize(0);
h_count13->GetXaxis()->SetLabelOffset(999);
h_count13->GetXaxis()->SetLabelSize(0);
h_count14->GetXaxis()->SetLabelOffset(999);
h_count14->GetXaxis()->SetLabelSize(0);
h_count15->GetXaxis()->SetLabelOffset(999);
h_count15->GetXaxis()->SetLabelSize(0);

h_count16->GetXaxis()->SetLabelOffset(999);
h_count16->GetXaxis()->SetLabelSize(0);
h_count17->GetXaxis()->SetLabelOffset(999);
h_count17->GetXaxis()->SetLabelSize(0);
h_count18->GetXaxis()->SetLabelOffset(999);
h_count18->GetXaxis()->SetLabelSize(0);
h_count19->GetXaxis()->SetLabelOffset(999);
h_count19->GetXaxis()->SetLabelSize(0);
h_count20->GetXaxis()->SetLabelOffset(999);
h_count20->GetXaxis()->SetLabelSize(0);
h_good_events->GetXaxis()->SetLabelOffset(999);
h_good_events->GetXaxis()->SetLabelSize(0);

//////////////////
double expn = 0.75;
double expn2 = 0.4;

double scl_hel_or = double(count_hel_or)*expn;
double scl_hel_and = double(count_hel_and)*expn2;
double scl_count_or = double(count_or)*expn;
double scl_count_and = double(count_and)*expn2;

double scl_layer1or = double(countlayer1or)*expn;
double scl_layer1and = double(countlayer1and)*expn2;
double scl_layer2or = double(countlayer2or)*expn;
double scl_layer2and = double(countlayer2and)*expn2;

double scl_layer3or = double(countlayer3or)*expn;
double scl_layer3and = double(countlayer3and)*expn2;
double scl_layer4or = double(countlayer4or)*expn;
double scl_layer4and = double(countlayer4and)*expn2;

double scl_hel_ora = double(count_hel_ora)*expn;
double scl_hel_anda = double(count_hel_anda)*expn2;
double scl_count_ora = double(count_ora)*expn;
double scl_count_anda = double(count_anda)*expn2;

double scl_labels1 = double(good_event)/-10.0;
double scl_labels2 = double(nentries)/-10.0;

////////////////////////////////////////////////////

char Name_Can2[100];			char Title_Can2[100];
sprintf(Name_Can2,"Run %d -- SFT_Efficiency.C",Run_Number);
sprintf(Title_Can2,"Run %d  -- SFT_Efficiency.C",Run_Number);

TCanvas *c1;
c1 = new TCanvas(Name_Can2,Title_Can2,1000,900); 
c1->Divide(2,3);
c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

c1->cd(3);

TLegend *l1;
l1 = new TLegend(0.75,0.75,0.99,0.99);
l1->AddEntry(h_good_events, "\"Good Events\"", "F");
l1->AddEntry(h_count2, "Helicity | Up OR Down", "F");
l1->AddEntry(h_count3, "Helicity | Up AND Down", "F");
l1->AddEntry(h_count4, "Layer | Up OR Down", "F");
l1->AddEntry(h_count5, "Layer | Up AND Down", "F");

h_good_events->SetLineColor(1);
h_good_events->SetFillStyle(3008);
h_good_events->SetFillColor(1);
h_good_events->Draw();

h_count2->SetLineColor(2);
h_count2->SetFillStyle(3004);
h_count2->SetFillColor(2);
h_count2->Draw("same");

h_count3->SetLineColor(7);
h_count3->SetFillStyle(3008);
h_count3->SetFillColor(7);
h_count3->Draw("same");

h_count4->SetLineColor(6);
h_count4->SetFillStyle(3004);
h_count4->SetFillColor(6);
h_count4->Draw("same");

h_count5->SetLineColor(3);
h_count5->SetFillStyle(3008);
h_count5->SetFillColor(3);
h_count5->Draw("same");

l1->Draw("same");
TLatex *tex4 = new TLatex(0.2,scl_labels1,"\"Good Events\"");
tex4->SetLineWidth(2);
tex4->Draw();
TLatex *tex5 = new TLatex(3.087098,scl_labels1,"Helicity");
tex5->SetLineWidth(2);
tex5->Draw();
TLatex *tex6 = new TLatex(5.195532,scl_labels1,"Layer");
tex6->SetLineWidth(2);
tex6->Draw();

TLatex *tex7 = new TLatex(3.1,scl_hel_or,HelorEff.c_str());
tex7->SetLineWidth(2);
tex7->Draw();
TLatex *tex8 = new TLatex(3.1,scl_hel_and,HelandEff.c_str());
tex8->SetLineWidth(2);
tex8->Draw();
TLatex *tex9 = new TLatex(5.1,scl_count_or,SFTorEff.c_str());
tex9->SetLineWidth(2);
tex9->Draw();
TLatex *tex10 = new TLatex(5.1,scl_count_and,SFTandEff.c_str());
tex10->SetLineWidth(2);
tex10->Draw();

///////////////////////////////////////////////////////////////////
c1->cd(1);

TLegend *l2;
l2 = new TLegend(0.75,0.75,0.99,0.99);
l2->AddEntry(h_good_events, "\"Good Events\"", "F");
l2->AddEntry(h_count2, "Layer 1 | Up OR Down", "F");
l2->AddEntry(h_count3, "Layer 1 | Up AND Down", "F");
l2->AddEntry(h_count4, "Layer 2 | Up OR Down", "F");
l2->AddEntry(h_count5, "Layer 2 | Up AND Down", "F");

h_good_events->SetLineColor(1);
h_good_events->SetFillStyle(3008);
h_good_events->SetFillColor(1);
h_good_events->Draw();

h_count7->SetLineColor(2);
h_count7->SetFillStyle(3004);
h_count7->SetFillColor(2);
h_count7->Draw("same");

h_count8->SetLineColor(7);
h_count8->SetFillStyle(3008);
h_count8->SetFillColor(7);
h_count8->Draw("same");

h_count9->SetLineColor(6);
h_count9->SetFillStyle(3004);
h_count9->SetFillColor(6);
h_count9->Draw("same");

h_count10->SetLineColor(3);
h_count10->SetFillStyle(3008);
h_count10->SetFillColor(3);
h_count10->Draw("same");

l2->Draw("same");
TLatex *tex11 = new TLatex(0.2,scl_labels1,"\"Good Events\"");
tex11->SetLineWidth(2);
tex11->Draw();
TLatex *tex12 = new TLatex(3.087098,scl_labels1,"Layer 1");
tex12->SetLineWidth(2);
tex12->Draw();
TLatex *tex13 = new TLatex(5.1,scl_labels1,"Layer 2");
tex13->SetLineWidth(2);
tex13->Draw();

TLatex *tex14 = new TLatex(3.1,scl_layer1or,L1orEff.c_str());
tex14->SetLineWidth(2);
tex14->Draw();
TLatex *tex15 = new TLatex(3.1,scl_layer1and,L1andEff.c_str());
tex15->SetLineWidth(2);
tex15->Draw();
TLatex *tex16 = new TLatex(5.1,scl_layer2or,L2orEff.c_str());
tex16->SetLineWidth(2);
tex16->Draw();
TLatex *tex17 = new TLatex(5.1,scl_layer2and,L2andEff.c_str());
tex17->SetLineWidth(2);
tex17->Draw();

///////////////////////////////////////////////////////
c1->cd(2);

TLegend *l3;
l3 = new TLegend(0.75,0.75,0.99,0.99);
l3->AddEntry(h_good_events, "\"Good Events\"", "F");
l3->AddEntry(h_count2, "Layer 3 | Up OR Down", "F");
l3->AddEntry(h_count3, "Layer 3 | Up AND Down", "F");
l3->AddEntry(h_count4, "Layer 4 | Up OR Down", "F");
l3->AddEntry(h_count5, "Layer 4 | Up AND Down", "F");

h_good_events->SetLineColor(1);
h_good_events->SetFillStyle(3008);
h_good_events->SetFillColor(1);
h_good_events->Draw();

h_count12->SetLineColor(2);
h_count12->SetFillStyle(3004);
h_count12->SetFillColor(2);
h_count12->Draw("same");

h_count13->SetLineColor(7);
h_count13->SetFillStyle(3008);
h_count13->SetFillColor(7);
h_count13->Draw("same");

h_count14->SetLineColor(6);
h_count14->SetFillStyle(3004);
h_count14->SetFillColor(6);
h_count14->Draw("same");

h_count15->SetLineColor(3);
h_count15->SetFillStyle(3008);
h_count15->SetFillColor(3);
h_count15->Draw("same");


l3->Draw("same");
TLatex *tex18 = new TLatex(0.2,scl_labels1,"\"Good Events\"");
tex18->SetLineWidth(2);
tex18->Draw();
TLatex *tex19 = new TLatex(3.087098,scl_labels1,"Layer 3");
tex19->SetLineWidth(2);
tex19->Draw();
TLatex *tex20 = new TLatex(5.1,scl_labels1,"Layer 4");
tex20->SetLineWidth(2);
tex20->Draw();

TLatex *tex21 = new TLatex(3.1,scl_layer3or,L3orEff.c_str());
tex21->SetLineWidth(2);
tex21->Draw();
TLatex *tex22 = new TLatex(3.1,scl_layer3and,L3andEff.c_str());
tex22->SetLineWidth(2);
tex22->Draw();
TLatex *tex23 = new TLatex(5.1,scl_layer4or,L4orEff.c_str());
tex23->SetLineWidth(2);
tex23->Draw();
TLatex *tex24 = new TLatex(5.1,scl_layer4and,L4andEff.c_str());
tex24->SetLineWidth(2);
tex24->Draw();

/////////////////////////////////

c1->cd(4);

TLegend *l4;
l4 = new TLegend(0.75,0.75,0.99,0.99);
l4->AddEntry(h_count16, "Total Number of Events", "F");
l4->AddEntry(h_good_events, "\"Good Events\"", "F");
l4->AddEntry(h_count2, "Helicity | Up OR Down", "F");
l4->AddEntry(h_count3, "Helicity | Up AND Down", "F");
l4->AddEntry(h_count4, "Layer | Up OR Down", "F");
l4->AddEntry(h_count5, "Layer | Up AND Down", "F");

h_count16->SetLineColor(1);
h_count16->SetFillStyle(3004);
h_count16->SetFillColor(1);
h_count16->Draw();

h_good_events->SetLineColor(1);
h_good_events->SetFillStyle(3008);
h_good_events->SetFillColor(1);
h_good_events->Draw("same");


h_count17->SetLineColor(2);
h_count17->SetFillStyle(3004);
h_count17->SetFillColor(2);
h_count17->Draw("same");

h_count18->SetLineColor(7);
h_count18->SetFillStyle(3008);
h_count18->SetFillColor(7);
h_count18->Draw("same");

h_count19->SetLineColor(6);
h_count19->SetFillStyle(3004);
h_count19->SetFillColor(6);
h_count19->Draw("same");

h_count20->SetLineColor(3);
h_count20->SetFillStyle(3008);
h_count20->SetFillColor(3);
h_count20->Draw("same");

l4->Draw("same");


TLatex *tex25 = new TLatex(0.2,scl_labels2,"Total Entries");
tex25->SetLineWidth(2);
tex25->Draw();
TLatex *tex26 = new TLatex(3.087098,scl_labels2,"Helicity");
tex26->SetLineWidth(2);
tex26->Draw();
TLatex *tex27 = new TLatex(5.195532,scl_labels2,"Layer");
tex27->SetLineWidth(2);
tex27->Draw();

TLatex *tex28 = new TLatex(3.1,scl_hel_ora,HelorEffa.c_str());
tex28->SetLineWidth(2);
tex28->Draw();
TLatex *tex29 = new TLatex(3.1,scl_hel_anda,HelandEffa.c_str());
tex29->SetLineWidth(2);
tex29->Draw();
TLatex *tex30 = new TLatex(5.1,scl_count_ora,SFTorEffa.c_str());
tex30->SetLineWidth(2);
tex30->Draw();
TLatex *tex31 = new TLatex(5.1,scl_count_anda,SFTandEffa.c_str());
tex31->SetLineWidth(2);
tex31->Draw();
////////////////

c1->cd(5);
hlayer1_2->Draw();

c1->cd(6);
hlayer3_4->Draw();

return;

} // End void
