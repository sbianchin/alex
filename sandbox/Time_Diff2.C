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
#include "TGaxis.h"
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
#include "TEllipse.h"
#include "TMarker.h"
#include "ANAPATH.h"
#include "Thresholds.h"
#include "CommonParameters.h"
#endif
void Time_Diff2(Int_t Run_Number=5)
{ 
	
	gStyle->Clear();
	TH1::AddDirectory(kFALSE);
	gStyle->SetOptStat(1111111);

	
	int TDC_Mean_min = 2400;
	int TDC_Mean_max = 2700;



 	Int_t adc_high_target[256];       Int_t ADC_High_TARGET[256];    
	Int_t adc_low_target[256];        Int_t ADC_Low_TARGET[256];  
	Int_t tdc_le_target[256][16];     Int_t TDC_LE_TARGET[256];     
	Int_t tdc_te_target[256][16];     //Int_t TDC_TE_TARGET[256];  

	Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];        Double_t ADC_High_SFT_corr[128];    
	Int_t adc_low_sft[128];           //Int_t ADC_Low_SFT[128];   
	Int_t tdc_le_sft[128][16];        Int_t TDC_LE_SFT[128];         
	Int_t tdc_te_sft[128][16];        //Int_t TDC_TE_SFT[128];  

	Int_t ADC_tof1[24];               Int_t ADC_TOF1[24]; 
	Int_t ADC_tof2[56];               Int_t ADC_TOF2[56];

	Int_t TDC_tof1U[12];              Int_t TDC_TOF1U[12];
	Int_t TDC_tof1D[12];              Int_t TDC_TOF1D[12];

	Float_t TDC_TOF1_Mean[12]={0.};

	Int_t TDC_tof2AO[12];             Int_t TDC_TOF2AO[12];
	Int_t TDC_tof2BO[12];             Int_t TDC_TOF2BO[12];
	Int_t TDC_tof2AI[12];             Int_t TDC_TOF2AI[12];
	Int_t TDC_tof2BI[12];             Int_t TDC_TOF2BI[12];

	Int_t MwpcADC[512];               Int_t MWPCADC[512];

	char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!

	char path_input[200];                   char file_mapping[200];
  	sprintf(path_input,"%s",path_merged);          

  	sprintf(file_mapping,"../Mapping");

  	char Name_finput[200];
  	sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);

  	char par_finput[200];
  	sprintf(par_finput,"%s/%s",file_mapping,source_mapping);
  	
  	char par_finput2[200];
  	sprintf(par_finput2,"%s/MWPC_map.txt",file_mapping);

	char par_finput3[200];
	sprintf(par_finput3,"%s/ADC_TARGET_Thresholds.txt",file_mapping);

	char par_finput4[200];
	sprintf(par_finput4,"%s/ADC_SFT_Thresholds.txt",file_mapping);

	char par_finput5[200];
	sprintf(par_finput5,"%s/ADC_TOF1_Thresholds.txt",file_mapping);

	char par_finput6[200];
	sprintf(par_finput6,"%s/ADC_TOF2_Thresholds.txt",file_mapping);


  	cout << "   " << endl;
  	cout << "File opened:  " << Name_finput << endl;
  	cout << "SFT Mapping File:  " << par_finput << endl;
  	cout << "MWPC Mapping File:  " << par_finput2 << endl;

  	//Change variables by applying Mapping
	Int_t par_temp[2][128];
	ifstream fdat(par_finput,ios::in);
	for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
	fdat.close();

	char par_temp2[512][50];
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
  
	// Histogram Declarations
	TH1F *h_TDC_SFT[128];			char h_TDC_SFT_Name[128][200];			char h_TDC_SFT_Title[128][200];
	TH1F *h_TDC_TOF1_UP[12];		char h_TDC_TOF1_UP_Name[12][200];		char h_TDC_TOF1_UP_Title[12][200];
	TH1F *h_TDC_TOF1_DOWN[12];		char h_TDC_TOF1_DOWN_Name[12][200];		char h_TDC_TOF1_DOWN_Title[12][200];
	TH1F *h_TDC_TOF1_MEAN[12];		char h_TDC_TOF1_MEAN_Name[12][200];		char h_TDC_TOF1_MEAN_Title[12][200];
	TH1F *h_Diff[128];				char h_Diff_Name[128][200];				char h_Diff_Title[128][200];

	TH1F *h_diff1[9];				char h_diff1_Name[9][200];				char h_diff1_Title[9][200];
	TH1F *h_diff2[9];				char h_diff2_Name[9][200];				char h_diff2_Title[9][200];

	sprintf(h_diff1_Title[0], "Run %d  |  TDC_TOF1_Mean(11) - TDC_TOF1_Mean(5)", Run_Number);
	sprintf(h_diff1_Name[0], "No TDC Cut");
	sprintf(h_diff1_Title[1], "Run %d  |  TDC_TOF1_Mean(11) - TDC_TOF1_Mean(6)", Run_Number);
	sprintf(h_diff1_Name[1], "No TDC Cut");
	sprintf(h_diff1_Title[2], "Run %d  |  TDC_TOF1_Mean(11) - TDC_TOF1_Mean(7)", Run_Number);
	sprintf(h_diff1_Name[2], "No TDC Cut");
	sprintf(h_diff1_Title[3], "Run %d  |  TDC_TOF1_Mean(12) - TDC_TOF1_Mean(5)", Run_Number);
	sprintf(h_diff1_Name[3], "No TDC Cut");
	sprintf(h_diff1_Title[4], "Run %d  |  TDC_TOF1_Mean(12) - TDC_TOF1_Mean(6)", Run_Number);
	sprintf(h_diff1_Name[4], "No TDC Cut");
	sprintf(h_diff1_Title[5], "Run %d  |  TDC_TOF1_Mean(12) - TDC_TOF1_Mean(7)", Run_Number);
	sprintf(h_diff1_Name[5], "No TDC Cut");
	sprintf(h_diff1_Title[6], "Run %d  |  TDC_TOF1_Mean(1) - TDC_TOF1_Mean(5)", Run_Number);
	sprintf(h_diff1_Name[6], "No TDC Cut");
	sprintf(h_diff1_Title[7], "Run %d  |  TDC_TOF1_Mean(1) - TDC_TOF1_Mean(6)", Run_Number);
	sprintf(h_diff1_Name[7], "No TDC Cut");
	sprintf(h_diff1_Title[8], "Run %d  |  TDC_TOF1_Mean(1) - TDC_TOF1_Mean(7)", Run_Number);
	sprintf(h_diff1_Name[8], "No TDC Cut");

	sprintf(h_diff2_Title[0], "Run %d  |  TDC_TOF1_Mean(11) - TDC_TOF1_Mean(5)", Run_Number);
	sprintf(h_diff2_Name[0], "With TDC Cut");
	sprintf(h_diff2_Title[1], "Run %d  |  TDC_TOF1_Mean(11) - TDC_TOF1_Mean(6)", Run_Number);
	sprintf(h_diff2_Name[1], "With TDC Cut");
	sprintf(h_diff2_Title[2], "Run %d  |  TDC_TOF1_Mean(11) - TDC_TOF1_Mean(7)", Run_Number);
	sprintf(h_diff2_Name[2], "With TDC Cut");
	sprintf(h_diff2_Title[3], "Run %d  |  TDC_TOF1_Mean(12) - TDC_TOF1_Mean(5)", Run_Number);
	sprintf(h_diff2_Name[3], "With TDC Cut");
	sprintf(h_diff2_Title[4], "Run %d  |  TDC_TOF1_Mean(12) - TDC_TOF1_Mean(6)", Run_Number);
	sprintf(h_diff2_Name[4], "With TDC Cut");
	sprintf(h_diff2_Title[5], "Run %d  |  TDC_TOF1_Mean(12) - TDC_TOF1_Mean(7)", Run_Number);
	sprintf(h_diff2_Name[5], "With TDC Cut");
	sprintf(h_diff2_Title[6], "Run %d  |  TDC_TOF1_Mean(1) - TDC_TOF1_Mean(5)", Run_Number);
	sprintf(h_diff2_Name[6], "With TDC Cut");
	sprintf(h_diff2_Title[7], "Run %d  |  TDC_TOF1_Mean(1) - TDC_TOF1_Mean(6)", Run_Number);
	sprintf(h_diff2_Name[7], "With TDC Cut");
	sprintf(h_diff2_Title[8], "Run %d  |  TDC_TOF1_Mean(1) - TDC_TOF1_Mean(7)", Run_Number);
	sprintf(h_diff2_Name[8], "With TDC Cut");

	for(int i=0; i<9; i++){
		h_diff1[i] = new TH1F(h_diff1_Name[i], h_diff1_Title[i], 4000, -2000, 2000);
		h_diff2[i] = new TH1F(h_diff2_Name[i], h_diff2_Title[i], 400, -200, 200);
	}

	for(int i=0; i<128; i++)
	{
		sprintf(h_TDC_SFT_Title[i], "SFT TDC Channel %d  |  Run %d", i, Run_Number);
		sprintf(h_TDC_SFT_Name[i], "Channel %d", i);

		sprintf(h_Diff_Title[i], "TDC(TOF1_Mean) - TDC(SFT)  |  Channel %d  --  Run %d", i, Run_Number);
		sprintf(h_Diff_Name[i], "Channel %d", i);

		h_TDC_SFT[i] = new TH1F(h_TDC_SFT_Name[i], h_TDC_SFT_Title[i], 1200, 0, 1200);
		h_Diff[i] = new TH1F(h_Diff_Name[i], h_Diff_Title[i], 500, 500, 1000);
	}

	for(int j=0; j<12; j++)
	{
		sprintf(h_TDC_TOF1_UP_Name[j], "TOF1 (TDC Ch. %d)", j);
		sprintf(h_TDC_TOF1_UP_Title[j], "Run %d  |  TDC_TOF1-UP  --  GAP %d", Run_Number, j+1);

		sprintf(h_TDC_TOF1_DOWN_Name[j], "TOF1 (TDC Ch. %d)", j+12);
		sprintf(h_TDC_TOF1_DOWN_Title[j], "Run %d  |  TDC_TOF1-DOWN  --  GAP %d", Run_Number, j+1);

		sprintf(h_TDC_TOF1_MEAN_Name[j], "TDC_Mean (Ch. %d)", j);
		sprintf(h_TDC_TOF1_MEAN_Title[j], "Run %d  |  TDC_TOF1-MEAN  --  GAP %d", Run_Number, j+1);

		h_TDC_TOF1_UP[j] = new TH1F(h_TDC_TOF1_UP_Name[j], h_TDC_TOF1_UP_Title[j], 1000, 0, 4100);
		h_TDC_TOF1_DOWN[j] = new TH1F(h_TDC_TOF1_DOWN_Name[j], h_TDC_TOF1_DOWN_Title[j], 1000, 0, 4100);
		h_TDC_TOF1_MEAN[j] = new TH1F(h_TDC_TOF1_MEAN_Name[j], h_TDC_TOF1_MEAN_Title[j], 2500, 1000, 3500);
	}

	int counter = 0;

	TChain *fChain= new TChain("Tree");   
	fChain->Add(Name_finput);   
	fChain->SetMakeClass(1);              

	fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);    fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
	fChain->SetBranchAddress("ADC_Low_TARGET",adc_low_target);      fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
	fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);        fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
	fChain->SetBranchAddress("TDC_TE_TARGET",tdc_te_target);        fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

	fChain->SetBranchAddress("ADC_TOF1",ADC_tof1);
	fChain->SetBranchAddress("ADC_TOF2",ADC_tof2);

	fChain->SetBranchAddress("TDC_TOF1U",TDC_tof1U);
	fChain->SetBranchAddress("TDC_TOF1D",TDC_tof1D);
	fChain->SetBranchAddress("TDC_TOF2AO",TDC_tof2AO);
	fChain->SetBranchAddress("TDC_TOF2AI",TDC_tof2AI);
	fChain->SetBranchAddress("TDC_TOF2BO",TDC_tof2BO);
	fChain->SetBranchAddress("TDC_TOF2BI",TDC_tof2BI);    

	fChain->SetBranchAddress("MWPCADC",MwpcADC);

	Int_t nentries = (Int_t)fChain->GetEntries();
	cout << "Total Number of Events:  " << nentries <<endl;
	cout << "  " << endl;
	cout << "*****************************************************" << endl;
	cout << "  " << endl;

	//for(Int_t ievent=0; ievent<1; ievent++)		// Loop over Events
	for(Int_t ievent=0; ievent<nentries; ievent++)
	{
    	fChain->GetEntry(ievent);  

    	if(ievent%10000==1) cout << "****  " << ievent << "  events done" << endl;

    	for(int init=0; init<12; init++)	TDC_TOF1_Mean[init]=0;

    	for(int j_TARGET=0; j_TARGET<256; j_TARGET++){
      		ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-ADC_cut_TARGET;
      		ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-ADC_cut_TARGET;
      		TDC_LE_TARGET[j_TARGET]=tdc_le_target[j_TARGET][0];
      		//TDC_TE_TARGET[j_TARGET]=tdc_te_target[j_TARGET][0];  
    	}

    	for(Int_t j_SFT=0; j_SFT<128; j_SFT++){
      		ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-par_temp_SFT[1][j_SFT];
      		//ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-par_temp_SFT[1][j_SFT];
      		TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
      		//TDC_TE_SFT[j_SFT]=tdc_te_sft[j_SFT][0]; 
    	}

    	for(Int_t j_TOF1=0; j_TOF1<24; j_TOF1++) {
      		ADC_TOF1[j_TOF1] = ADC_tof1[j_TOF1]-par_temp_TOF1[1][j_TOF1];
    	}

    	for(Int_t j_TOF2=0; j_TOF2<56; j_TOF2++) {
      		ADC_TOF2[j_TOF2] = ADC_tof2[j_TOF2]-par_temp_TOF2[1][j_TOF2];
    	}

    	for(Int_t j_MWPC=0; j_MWPC<512; j_MWPC++) {
      		MWPCADC[j_MWPC] = MwpcADC[j_MWPC]-MWPC_thr;
    	}

    	for(Int_t j_TDCTOF=0; j_TDCTOF<12; j_TDCTOF++) {
      		TDC_TOF1U[j_TDCTOF] = TDC_tof1U[j_TDCTOF];
      		TDC_TOF1D[j_TDCTOF] = TDC_tof1D[j_TDCTOF];
      		TDC_TOF2AO[j_TDCTOF] = TDC_tof2AO[j_TDCTOF];
      		TDC_TOF2BO[j_TDCTOF] = TDC_tof2BO[j_TDCTOF];
      		TDC_TOF2AI[j_TDCTOF] = TDC_tof2AI[j_TDCTOF];
      		TDC_TOF2BI[j_TDCTOF] = TDC_tof2BI[j_TDCTOF];

      		//if((TDC_TOF1U[j_TDCTOF]>=TDC_TOF1_min && TDC_TOF1U[j_TDCTOF]<=TDC_TOF1_max) && 
      		//	(TDC_TOF1D[j_TDCTOF]>=TDC_TOF1_min && TDC_TOF1D[j_TDCTOF]<=TDC_TOF1_max))
      		if(TDC_TOF1U[j_TDCTOF]>0 && TDC_TOF1D[j_TDCTOF]>0)
      		{
      			TDC_TOF1_Mean[j_TDCTOF] = 0.5 * (TDC_TOF1U[j_TDCTOF] + TDC_TOF1D[j_TDCTOF]);
      		}

      		//if((TDC_TOF1U[j_TDCTOF]>=TDC_TOF1_min && TDC_TOF1U[j_TDCTOF]<=TDC_TOF1_max) && 
      		//	(TDC_TOF1D[j_TDCTOF]<TDC_TOF1_min || TDC_TOF1D[j_TDCTOF]>TDC_TOF1_max))
        	if(TDC_TOF1U[j_TDCTOF]>0 && TDC_TOF1D[j_TDCTOF]<=0)

      		{
      			TDC_TOF1_Mean[j_TDCTOF] = TDC_TOF1U[j_TDCTOF];
      		}

      		//if((TDC_TOF1D[j_TDCTOF]>=TDC_TOF1_min && TDC_TOF1D[j_TDCTOF]<=TDC_TOF1_max) && 
      		//	(TDC_TOF1U[j_TDCTOF]<TDC_TOF1_min || TDC_TOF1U[j_TDCTOF]>TDC_TOF1_max))
        	if(TDC_TOF1U[j_TDCTOF]<=0 && TDC_TOF1D[j_TDCTOF]>0)
      		{
      			TDC_TOF1_Mean[j_TDCTOF] = TDC_TOF1D[j_TDCTOF];
      		}   		

        	if(TDC_TOF1U[j_TDCTOF]<=0 && TDC_TOF1D[j_TDCTOF]<=0)
      		{
      			TDC_TOF1_Mean[j_TDCTOF] = -999;
      		}   		
      	}

      		if(TDC_TOF1_Mean[10]>0 && TDC_TOF1_Mean[4]>0)	h_diff1[0]->Fill(TDC_TOF1_Mean[10] - TDC_TOF1_Mean[4]);
       		if(TDC_TOF1_Mean[10]>0 && TDC_TOF1_Mean[5]>0)	h_diff1[1]->Fill(TDC_TOF1_Mean[10] - TDC_TOF1_Mean[5]);
      		if(TDC_TOF1_Mean[10]>0 && TDC_TOF1_Mean[6]>0)	h_diff1[2]->Fill(TDC_TOF1_Mean[10] - TDC_TOF1_Mean[6]);
      		if(TDC_TOF1_Mean[11]>0 && TDC_TOF1_Mean[4]>0)	h_diff1[3]->Fill(TDC_TOF1_Mean[11] - TDC_TOF1_Mean[4]);
      		if(TDC_TOF1_Mean[11]>0 && TDC_TOF1_Mean[5]>0)	h_diff1[4]->Fill(TDC_TOF1_Mean[11] - TDC_TOF1_Mean[5]);
      		if(TDC_TOF1_Mean[11]>0 && TDC_TOF1_Mean[6]>0)	h_diff1[5]->Fill(TDC_TOF1_Mean[11] - TDC_TOF1_Mean[6]);
      		if(TDC_TOF1_Mean[0]>0 && TDC_TOF1_Mean[4]>0)	h_diff1[6]->Fill(TDC_TOF1_Mean[0] - TDC_TOF1_Mean[4]);
      		if(TDC_TOF1_Mean[0]>0 && TDC_TOF1_Mean[5]>0)	h_diff1[7]->Fill(TDC_TOF1_Mean[0] - TDC_TOF1_Mean[5]);
      		if(TDC_TOF1_Mean[0]>0 && TDC_TOF1_Mean[6]>0)	h_diff1[8]->Fill(TDC_TOF1_Mean[0] - TDC_TOF1_Mean[6]);
     		
      		if(TDC_TOF1_Mean[10]>=TDC_Mean_min && TDC_TOF1_Mean[10]<=TDC_Mean_max && TDC_TOF1_Mean[4]>=TDC_Mean_min && TDC_TOF1_Mean[4]<=TDC_Mean_max){
      			h_diff2[0]->Fill(TDC_TOF1_Mean[10] - TDC_TOF1_Mean[4]);
      		}

      		if(TDC_TOF1_Mean[10]>=TDC_Mean_min && TDC_TOF1_Mean[10]<=TDC_Mean_max && TDC_TOF1_Mean[5]>=TDC_Mean_min && TDC_TOF1_Mean[5]<=TDC_Mean_max){
      			h_diff2[1]->Fill(TDC_TOF1_Mean[10] - TDC_TOF1_Mean[5]);
      		}

      		if(TDC_TOF1_Mean[10]>=TDC_Mean_min && TDC_TOF1_Mean[10]<=TDC_Mean_max && TDC_TOF1_Mean[6]>=TDC_Mean_min && TDC_TOF1_Mean[6]<=TDC_Mean_max){
      			h_diff2[2]->Fill(TDC_TOF1_Mean[10] - TDC_TOF1_Mean[6]);
      		}

      		if(TDC_TOF1_Mean[11]>=TDC_Mean_min && TDC_TOF1_Mean[11]<=TDC_Mean_max && TDC_TOF1_Mean[4]>=TDC_Mean_min && TDC_TOF1_Mean[4]<=TDC_Mean_max){
      			h_diff2[3]->Fill(TDC_TOF1_Mean[11] - TDC_TOF1_Mean[4]);
      		}

      		if(TDC_TOF1_Mean[11]>=TDC_Mean_min && TDC_TOF1_Mean[11]<=TDC_Mean_max && TDC_TOF1_Mean[5]>=TDC_Mean_min && TDC_TOF1_Mean[5]<=TDC_Mean_max){
      			h_diff2[4]->Fill(TDC_TOF1_Mean[11] - TDC_TOF1_Mean[5]);
      		}

      		if(TDC_TOF1_Mean[11]>=TDC_Mean_min && TDC_TOF1_Mean[11]<=TDC_Mean_max && TDC_TOF1_Mean[6]>=TDC_Mean_min && TDC_TOF1_Mean[6]<=TDC_Mean_max){
      			h_diff2[5]->Fill(TDC_TOF1_Mean[11] - TDC_TOF1_Mean[6]);
      		}

      		if(TDC_TOF1_Mean[0]>=TDC_Mean_min && TDC_TOF1_Mean[0]<=TDC_Mean_max && TDC_TOF1_Mean[4]>=TDC_Mean_min && TDC_TOF1_Mean[4]<=TDC_Mean_max){
      			h_diff2[6]->Fill(TDC_TOF1_Mean[0] - TDC_TOF1_Mean[4]);
      		}

      		if(TDC_TOF1_Mean[0]>=TDC_Mean_min && TDC_TOF1_Mean[0]<=TDC_Mean_max && TDC_TOF1_Mean[5]>=TDC_Mean_min && TDC_TOF1_Mean[5]<=TDC_Mean_max){
      			h_diff2[7]->Fill(TDC_TOF1_Mean[0] - TDC_TOF1_Mean[5]);
      		}

      		if(TDC_TOF1_Mean[0]>=TDC_Mean_min && TDC_TOF1_Mean[0]<=TDC_Mean_max && TDC_TOF1_Mean[6]>=TDC_Mean_min && TDC_TOF1_Mean[6]<=TDC_Mean_max){
      			h_diff2[8]->Fill(TDC_TOF1_Mean[0] - TDC_TOF1_Mean[6]);
      		}


//      		cout << ievent << "   " << TDC_TOF1U[0] << "    " << TDC_TOF1D[0] << "    " << TDC_TOF1_Mean[0] << endl; 
      		//if(TDC_TOF1_Mean[10]>0 && TDC_TOF1_Mean[4]>0)
      		//{
      		//	counter++;
      		//	cout << counter << "   " << ievent << "   " << TDC_TOF1_Mean[10] << "    " << TDC_TOF1_Mean[4] << "    " << TDC_TOF1_Mean[10] - TDC_TOF1_Mean[4] << endl;
      		//}
      		//cout << ievent << "   " << TDC_TOF1_Mean[10] - TDC_TOF1_Mean[4] << "    " << TDC_TOF1_Mean[10] - TDC_TOF1_Mean[5] << "    " << TDC_TOF1_Mean[10] - TDC_TOF1_Mean[6] << endl; 
          	//	cout << "Event  " << ievent << ":  " << "TDC TOF1 Mean Channel " << j_TDCTOF << " :  " << TDC_TOF1_Mean[j_TDCTOF] << endl;
    	
  	
    	for(int i=0; i<128; i++)
    	{
    		for(int j=0; j<12; j++)
    		{
    			if(TDC_LE_SFT[i]>=TDC_min_SFT && TDC_LE_SFT[i]<=TDC_max_SFT && TDC_TOF1_Mean[j]>0)
    			{
    		//		cout << ievent << "  " << i << "  " << j << "  " << TDC_LE_SFT[i] << "  " << TDC_TOF1_Mean[j] << "   " << float(TDC_LE_SFT[i]) - 0.025*float(TDC_TOF1_Mean[j])<< endl;
    				h_Diff[i]->Fill(float(TDC_LE_SFT[i]) - 0.025*float(TDC_TOF1_Mean[j]));
   			}
    		}
    	}


    	//cout << "  " << endl;

		for(int j=0; j<128; j++)
		{ 
			if(ADC_High_SFT[j]<0)      ADC_High_SFT_corr[j]=0; 
			if(ADC_High_SFT[j]>=0)     ADC_High_SFT_corr[j]=ADC_High_SFT[j]; 
		}

		//for(int k=0; k<128; k++)
		//{
		//	cout << "SFT Channel " << k << " (TDC):  " << TDC_LE_SFT[k] << endl;;
		//}

		//cout << "  " << endl;
		//cout << "///////////////////////////////////////////////" << endl;
		//cout << "  " << endl;

		//for(int k=0; k<12; k++)
		//{
		//	cout << "TOF1-Up " << k << " (TDC):  " << TDC_TOF1U[k] << endl;
		//}

		//cout << "  " << endl;
		//cout << "///////////////////////////////////////////////" << endl;
		//cout << "  " << endl;

		//for(int k=0; k<12; k++)
		//{
		//	cout << "TOF1-Down " << k << " (TDC):  " << TDC_TOF1D[k] << endl;
		//}


		for(int j=0; j<128; j++)
		{
			if(TDC_LE_SFT[j]>=0)
			{
				h_TDC_SFT[j]->Fill(TDC_LE_SFT[j]);
			}
		}

		for(int k=0; k<12; k++)
		{
			//if(TDC_TOF1U[k]>=0)
			//{
				h_TDC_TOF1_UP[k]->Fill(TDC_TOF1U[k]);
			//}

			//if(TDC_TOF1D[k]>=0)
			//{
				h_TDC_TOF1_DOWN[k]->Fill(TDC_TOF1D[k]);
			//}

			//if(TDC_TOF1_Mean[k]>0)
			//{
				h_TDC_TOF1_MEAN[k]->Fill(TDC_TOF1_Mean[k]);
			//}
		}


	} // EndLoop over Events

	// Canvas Declaration

	char c_Diff1_Title[100];
	char c_Diff2_Title[100];

	sprintf(c_Diff1_Title,"TDC(TOF1_Mean) - TDC(SFT)  |  Run %d  --  Channels 0 - 63", Run_Number);
	sprintf(c_Diff2_Title,"TDC(TOF1_Mean) - TDC(SFT)  |  Run %d  --  Channels 63 - 127", Run_Number);
	
	//TCanvas *c_SFT;
	//c_SFT = new TCanvas("SFT","SFT",0,200,1050,700);
	//c_SFT->Divide(16,8);
	//c_SFT->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	//TCanvas *c_TOF1;
	//c_TOF1 = new TCanvas("TOF1","TOF1",0,200,1050,700);
	//c_TOF1->Divide(6,4);
	//c_TOF1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	//TCanvas *c_TOF1_MEAN;
	//c_TOF1_MEAN = new TCanvas("TOF1 MEAN","TOF1 MEAN",0,200,1050,700);
	//c_TOF1_MEAN->Divide(6,4);
	//c_TOF1_MEAN->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	//TCanvas *c_Diff;
	//c_Diff = new TCanvas("Diff","Diff",0,200,1050,700);
	//c_Diff->Divide(16,8);
	//c_Diff->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	TCanvas *c1;
	c1 = new TCanvas("c1","c1",0,200,1200,500);
	c1->Divide(8,6);
	c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	TCanvas *c2;
	c2 = new TCanvas("c2","c2",500,3000,1200,500);
	c2->Divide(7,3);
	c2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");


	//TCanvas *c_Diff1;
	//c_Diff1 = new TCanvas("Diff1",c_Diff1_Title,0,200,1050,700);
	//c_Diff1->Divide(8,8);
	//c_Diff1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	TLatex *tex_condition1[128];
	TLatex *tex_condition2[128];

	char cond1[100];	sprintf(cond1,"%d #leq TDC(TOF1) #leq %d", TDC_TOF1_min, TDC_TOF1_max);
	char cond2[100];	sprintf(cond2,"%d #leq TDC(SFT) #leq %d", TDC_min_SFT, TDC_max_SFT);

	int max=0;



	c1->cd(1);
	h_TDC_TOF1_UP[0]->Draw();
	c1->cd(2);
	h_TDC_TOF1_DOWN[0]->Draw();
	c1->cd(3);
	h_TDC_TOF1_UP[1]->Draw();
	c1->cd(4);
	h_TDC_TOF1_DOWN[1]->Draw();
	c1->cd(5);
	h_TDC_TOF1_UP[2]->Draw();
	c1->cd(6);
	h_TDC_TOF1_DOWN[2]->Draw();
	c1->cd(7);
	h_TDC_TOF1_UP[3]->Draw();
	c1->cd(8);
	h_TDC_TOF1_DOWN[3]->Draw();
	c1->cd(9);
	h_TDC_TOF1_UP[4]->Draw();
	c1->cd(10);
	h_TDC_TOF1_DOWN[4]->Draw();
	c1->cd(11);
	h_TDC_TOF1_UP[5]->Draw();
	c1->cd(12);
	h_TDC_TOF1_DOWN[5]->Draw();

	c1->cd(13);
	h_TDC_TOF1_UP[6]->Draw();
	c1->cd(14);
	h_TDC_TOF1_DOWN[6]->Draw();
	c1->cd(15);
	h_TDC_TOF1_UP[7]->Draw();
	c1->cd(16);
	h_TDC_TOF1_DOWN[7]->Draw();
	c1->cd(17);
	h_TDC_TOF1_UP[8]->Draw();
	c1->cd(18);
	h_TDC_TOF1_DOWN[8]->Draw();
	c1->cd(19);
	h_TDC_TOF1_UP[9]->Draw();
	c1->cd(20);
	h_TDC_TOF1_DOWN[9]->Draw();
	c1->cd(21);
	h_TDC_TOF1_UP[10]->Draw();
	c1->cd(22);
	h_TDC_TOF1_DOWN[10]->Draw();
	c1->cd(23);
	h_TDC_TOF1_UP[11]->Draw();
	c1->cd(24);
	h_TDC_TOF1_DOWN[11]->Draw();

	int max_mean[12] = {0};

	for(int i=0; i<12; i++)  max_mean[i] = h_TDC_TOF1_MEAN[i]->GetMaximum();

	TLine *Min_Line1 = new TLine(2400, 0, 2400, 1.05*max_mean[0]);		TLine *Max_Line1 = new TLine(2700, 0, 2700, 1.05*max_mean[0]);
	TLine *Min_Line2 = new TLine(2400, 0, 2400, 1.05*max_mean[1]);		TLine *Max_Line2 = new TLine(2700, 0, 2700, 1.05*max_mean[1]);
	TLine *Min_Line3 = new TLine(2400, 0, 2400, 1.05*max_mean[2]);		TLine *Max_Line3 = new TLine(2700, 0, 2700, 1.05*max_mean[2]);
	TLine *Min_Line4 = new TLine(2400, 0, 2400, 1.05*max_mean[3]);		TLine *Max_Line4 = new TLine(2700, 0, 2700, 1.05*max_mean[3]);
	TLine *Min_Line5 = new TLine(2400, 0, 2400, 1.05*max_mean[4]);		TLine *Max_Line5 = new TLine(2700, 0, 2700, 1.05*max_mean[4]);
	TLine *Min_Line6 = new TLine(2400, 0, 2400, 1.05*max_mean[5]);		TLine *Max_Line6 = new TLine(2700, 0, 2700, 1.05*max_mean[5]);
	TLine *Min_Line7 = new TLine(2400, 0, 2400, 1.05*max_mean[6]);		TLine *Max_Line7 = new TLine(2700, 0, 2700, 1.05*max_mean[6]);
	TLine *Min_Line8 = new TLine(2400, 0, 2400, 1.05*max_mean[7]);		TLine *Max_Line8 = new TLine(2700, 0, 2700, 1.05*max_mean[7]);
	TLine *Min_Line9 = new TLine(2400, 0, 2400, 1.05*max_mean[8]);		TLine *Max_Line9 = new TLine(2700, 0, 2700, 1.05*max_mean[8]);
	TLine *Min_Line10 = new TLine(2400, 0, 2400, 1.05*max_mean[9]);		TLine *Max_Line10 = new TLine(2700, 0, 2700, 1.05*max_mean[9]);
	TLine *Min_Line11 = new TLine(2400, 0, 2400, 1.05*max_mean[10]);	TLine *Max_Line11 = new TLine(2700, 0, 2700, 1.05*max_mean[10]);
	TLine *Min_Line12 = new TLine(2400, 0, 2400, 1.05*max_mean[11]);	TLine *Max_Line12 = new TLine(2700, 0, 2700, 1.05*max_mean[11]);
	

	Min_Line1->SetLineStyle(2);			Max_Line1->SetLineStyle(2);
	Min_Line1->SetLineWidth(2);			Max_Line1->SetLineWidth(2);
	Min_Line1->SetLineColor(2); 		Max_Line1->SetLineColor(2); 
	
	Min_Line2->SetLineStyle(2);			Max_Line2->SetLineStyle(2);
	Min_Line2->SetLineWidth(2);			Max_Line2->SetLineWidth(2);
	Min_Line2->SetLineColor(2); 		Max_Line2->SetLineColor(2); 
	
	Min_Line3->SetLineStyle(2);			Max_Line3->SetLineStyle(2);
	Min_Line3->SetLineWidth(2);			Max_Line3->SetLineWidth(2);
	Min_Line3->SetLineColor(2); 		Max_Line3->SetLineColor(2); 
	
	Min_Line4->SetLineStyle(2);			Max_Line4->SetLineStyle(2);
	Min_Line4->SetLineWidth(2);			Max_Line4->SetLineWidth(2);
	Min_Line4->SetLineColor(2); 		Max_Line4->SetLineColor(2); 
	
	Min_Line5->SetLineStyle(2);			Max_Line5->SetLineStyle(2);
	Min_Line5->SetLineWidth(2);			Max_Line5->SetLineWidth(2);
	Min_Line5->SetLineColor(2); 		Max_Line5->SetLineColor(2); 
	
	Min_Line6->SetLineStyle(2);			Max_Line6->SetLineStyle(2);
	Min_Line6->SetLineWidth(2);			Max_Line6->SetLineWidth(2);
	Min_Line6->SetLineColor(2); 		Max_Line6->SetLineColor(2); 
	
	Min_Line7->SetLineStyle(2);			Max_Line7->SetLineStyle(2);
	Min_Line7->SetLineWidth(2);			Max_Line7->SetLineWidth(2);
	Min_Line7->SetLineColor(2); 		Max_Line7->SetLineColor(2); 
	
	Min_Line8->SetLineStyle(2);			Max_Line8->SetLineStyle(2);
	Min_Line8->SetLineWidth(2);			Max_Line8->SetLineWidth(2);
	Min_Line8->SetLineColor(2); 		Max_Line8->SetLineColor(2); 
	
	Min_Line9->SetLineStyle(2);			Max_Line9->SetLineStyle(2);
	Min_Line9->SetLineWidth(2);			Max_Line9->SetLineWidth(2);
	Min_Line9->SetLineColor(2); 		Max_Line9->SetLineColor(2); 
	
	Min_Line10->SetLineStyle(2);		Max_Line10->SetLineStyle(2);
	Min_Line10->SetLineWidth(2);		Max_Line10->SetLineWidth(2);
	Min_Line10->SetLineColor(2); 		Max_Line10->SetLineColor(2); 
	
	Min_Line11->SetLineStyle(2);		Max_Line11->SetLineStyle(2);
	Min_Line11->SetLineWidth(2);		Max_Line11->SetLineWidth(2);
	Min_Line11->SetLineColor(2); 		Max_Line11->SetLineColor(2); 
	
	Min_Line12->SetLineStyle(2);		Max_Line12->SetLineStyle(2);
	Min_Line12->SetLineWidth(2);		Max_Line12->SetLineWidth(2);
	Min_Line12->SetLineColor(2); 		Max_Line12->SetLineColor(2); 


	c1->cd(34);
	h_TDC_TOF1_MEAN[0]->Draw();
	Min_Line1->Draw();
	Max_Line1->Draw();

	c1->cd(35);
	h_TDC_TOF1_MEAN[1]->Draw();
	Min_Line2->Draw();
	Max_Line2->Draw();

	c1->cd(36);
	h_TDC_TOF1_MEAN[2]->Draw();
	Min_Line3->Draw();
	Max_Line3->Draw();

	c1->cd(37);
	h_TDC_TOF1_MEAN[3]->Draw();
	Min_Line4->Draw();
	Max_Line4->Draw();

	c1->cd(38);
	h_TDC_TOF1_MEAN[4]->Draw();
	Min_Line5->Draw();
	Max_Line5->Draw();

	c1->cd(39);
	h_TDC_TOF1_MEAN[5]->Draw();
	Min_Line6->Draw();
	Max_Line6->Draw();

	c1->cd(42);
	h_TDC_TOF1_MEAN[6]->Draw();
	Min_Line7->Draw();
	Max_Line7->Draw();

	c1->cd(43);
	h_TDC_TOF1_MEAN[7]->Draw();
	Min_Line8->Draw();
	Max_Line8->Draw();

	c1->cd(44);
	h_TDC_TOF1_MEAN[8]->Draw();
	Min_Line9->Draw();
	Max_Line9->Draw();

	c1->cd(45);
	h_TDC_TOF1_MEAN[9]->Draw();
	Min_Line10->Draw();
	Max_Line10->Draw();

	c1->cd(46);
	h_TDC_TOF1_MEAN[10]->Draw();
	Min_Line11->Draw();
	Max_Line11->Draw();

	c1->cd(47);
	h_TDC_TOF1_MEAN[11]->Draw();
	Min_Line12->Draw();
	Max_Line12->Draw();

	TF1 *fitline_Gap1;      
	h_diff2[0]->Fit("gaus", "QCM");
	fitline_Gap1 = h_diff2[0]->GetFunction("gaus");
	fitline_Gap1->SetLineWidth(2);
	fitline_Gap1->SetLineColor(2);

	
	TF1 *fitline_Gap2;      
	h_diff2[1]->Fit("gaus", "QCM");
	fitline_Gap2 = h_diff2[1]->GetFunction("gaus");
	fitline_Gap2->SetLineWidth(2);
	fitline_Gap2->SetLineColor(2);
	
	TF1 *fitline_Gap3;      
	h_diff2[2]->Fit("gaus", "QCM");
	fitline_Gap3 = h_diff2[2]->GetFunction("gaus");
	fitline_Gap3->SetLineWidth(2);
	fitline_Gap3->SetLineColor(2);
	
	TF1 *fitline_Gap4;      
	h_diff2[3]->Fit("gaus", "QCM");
	fitline_Gap4 = h_diff2[3]->GetFunction("gaus");
	fitline_Gap4->SetLineWidth(2);
	fitline_Gap4->SetLineColor(2);
	
	TF1 *fitline_Gap5;      
	h_diff2[4]->Fit("gaus", "QCM");
	fitline_Gap5 = h_diff2[4]->GetFunction("gaus");
	fitline_Gap5->SetLineWidth(2);
	fitline_Gap5->SetLineColor(2);
	
	TF1 *fitline_Gap6;      
	h_diff2[5]->Fit("gaus", "QCM");
	fitline_Gap6 = h_diff2[5]->GetFunction("gaus");
	fitline_Gap6->SetLineWidth(2);
	fitline_Gap6->SetLineColor(2);
	
	TF1 *fitline_Gap7;      
	h_diff2[6]->Fit("gaus", "QCM");
	fitline_Gap7 = h_diff2[6]->GetFunction("gaus");
	fitline_Gap7->SetLineWidth(2);
	fitline_Gap7->SetLineColor(2);
	
	TF1 *fitline_Gap8;      
	h_diff2[7]->Fit("gaus", "QCM");
	fitline_Gap8 = h_diff2[7]->GetFunction("gaus");
	fitline_Gap8->SetLineWidth(2);
	fitline_Gap8->SetLineColor(2);
	
	TF1 *fitline_Gap9;      
	h_diff2[8]->Fit("gaus", "QCM");
	fitline_Gap9 = h_diff2[8]->GetFunction("gaus");
	fitline_Gap9->SetLineWidth(2);
	fitline_Gap9->SetLineColor(2);

	float sigma_fit[9]={0.};	float mean_fit[9]={0.};

	mean_fit[0] = fitline_Gap1->GetParameter(1);	sigma_fit[0] = fitline_Gap1->GetParameter(2);
	mean_fit[1] = fitline_Gap2->GetParameter(1);	sigma_fit[1] = fitline_Gap2->GetParameter(2);
	mean_fit[2] = fitline_Gap3->GetParameter(1);	sigma_fit[2] = fitline_Gap3->GetParameter(2);
	mean_fit[3] = fitline_Gap4->GetParameter(1);	sigma_fit[3] = fitline_Gap4->GetParameter(2);
	mean_fit[4] = fitline_Gap5->GetParameter(1);	sigma_fit[4] = fitline_Gap5->GetParameter(2);
	mean_fit[5] = fitline_Gap6->GetParameter(1);	sigma_fit[5] = fitline_Gap6->GetParameter(2);
	mean_fit[6] = fitline_Gap7->GetParameter(1);	sigma_fit[6] = fitline_Gap7->GetParameter(2);
	mean_fit[7] = fitline_Gap8->GetParameter(1);	sigma_fit[7] = fitline_Gap8->GetParameter(2);
	mean_fit[8] = fitline_Gap9->GetParameter(1);	sigma_fit[8] = fitline_Gap9->GetParameter(2);


	int max_diff2[9]={0};
	for(int j=0; j<9; j++)	max_diff2[j]=h_diff2[j]->GetMaximum();


	TLatex *Delta_T[9]={0};

	char delta_t[9][100];

	for(int i=0; i<9; i++){
		sprintf(delta_t[i], "#sigma = %3.4f ns", 0.025*sigma_fit[i]/sqrt(2));
		Delta_T[i] = new TLatex(-190, 0.95*max_diff2[i], delta_t[i]);
		Delta_T[i]->SetTextSize(0.035);
	}

	TF1 *f1=new TF1("f1", "x", -50, 50);
    TGaxis *A1 = new TGaxis(-200,0,200,0,"f1",510,"-");


	c2->cd(1);
	h_diff1[0]->Draw();
	c2->cd(2);
	h_diff1[1]->Draw();
	c2->cd(3);
	h_diff1[2]->Draw();

	c2->cd(8);
	h_diff1[3]->Draw();
	c2->cd(9);
	h_diff1[4]->Draw();
	c2->cd(10);
	h_diff1[5]->Draw();

	c2->cd(15);
	h_diff1[6]->Draw();
	c2->cd(16);
	h_diff1[7]->Draw();
	c2->cd(17);
	h_diff1[8]->Draw();

	c2->cd(5);
	h_diff2[0]->Draw();
	Delta_T[0]->Draw("same");

	c2->cd(6);
	h_diff2[1]->Draw();
	Delta_T[1]->Draw("same");

	c2->cd(7);
	h_diff2[2]->Draw();
	Delta_T[2]->Draw("same");

	c2->cd(12);
	h_diff2[3]->Draw();
	Delta_T[3]->Draw("same");

	c2->cd(13);
	h_diff2[4]->Draw();
	Delta_T[4]->Draw("same");

	c2->cd(14);
	h_diff2[5]->Draw();
	Delta_T[5]->Draw("same");


	c2->cd(19);
	h_diff2[6]->Draw();
	Delta_T[6]->Draw("same");

	c2->cd(20);
	h_diff2[7]->Draw();
	Delta_T[7]->Draw("same");

	c2->cd(21);
	h_diff2[8]->Draw();
	Delta_T[8]->Draw("same");


	c1->cd(47);
	h_TDC_TOF1_MEAN[11]->Draw();
	Min_Line12->Draw();
	Max_Line12->Draw();




} //EndVoid