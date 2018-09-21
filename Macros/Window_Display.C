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
#include "ADC_Thresholds.h"
#include "TDC_Windows.h"

#endif


void Window_Display(Int_t run_number=3994, Int_t flag=0, Int_t gap_selected=0, Int_t nevt=0) { 
  
  	gStyle->Clear();
  	TH1::AddDirectory(kFALSE);
	gStyle->SetOptStat(111111111);



	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int gap_counter[12];
 	int high_gap_hit;
  	int gap_to_fit;
  	int gap_TOF2;
  	int list_counter=0;

  	char Event_List_Title[100];
  	sprintf(Event_List_Title,"Event_List_Run%d_TOF1_%d.txt",run_number,gap_selected);

	ofstream Event_List;
	Event_List.open (Event_List_Title);

	// SFT
	Int_t ADC_High_sft[128];
	Int_t ADC_Low_sft[128];
	Int_t TDC_LE_sft[128][16];
	Int_t TDC_TE_sft[128][16];

	// TARGET
	Int_t adc_high_target[256]; 	
	Int_t adc_low_target[256]; 		
	Int_t tdc_le_target[256][16];      	  	
	Int_t tdc_te_target[256][16]; 		

	// TOF1
	Int_t ADC_tof1U[12];	Int_t TDC_tof1U[12];
	Int_t ADC_tof1D[12];	Int_t TDC_tof1D[12];

	// TOF2
	Int_t ADC_tof2AO[12];	Int_t TDC_tof2AO[12];
	Int_t ADC_tof2AI[12];	Int_t TDC_tof2AI[12];
	Int_t ADC_tof2BO[12];	Int_t TDC_tof2BO[12];
	Int_t ADC_tof2BI[12];	Int_t TDC_tof2BI[12];

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	char path_input[200];                   //char file_mapping[200];
	sprintf(path_input,"%s",path_merged);          
	//sprintf(file_mapping,path_mapping);
	//sprintf(path_input,"/triumfcs/trshare/trek/E36/Data/April_2015/root");

	char Name_finput[200];
	sprintf(Name_finput,"%s/Run%dMS.root",path_input, run_number);


	cout << "   " << endl;
	cout << Name_finput << endl;

	TChain *fChain= new TChain("Tree");		
	fChain->Add(Name_finput);		
	fChain->SetMakeClass(1);							


	///////////////////////////////////////////////////////

	// SFT
	if(flag==1 || flag==2){
		fChain->SetBranchAddress("ADC_High_SFT",ADC_High_sft);
		fChain->SetBranchAddress("ADC_Low_SFT",ADC_Low_sft);
		fChain->SetBranchAddress("TDC_LE_SFT",TDC_LE_sft);
		fChain->SetBranchAddress("TDC_TE_SFT",TDC_TE_sft);
	}

	// TARGET
	if(flag==3 || flag==4 || flag==5 || flag==6){
		fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);		
		fChain->SetBranchAddress("ADC_Low_TARGET",adc_low_target);
		fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);		
		fChain->SetBranchAddress("TDC_TE_TARGET",tdc_te_target);
	}

	// TOF1
	//if(flag==7){
		fChain->SetBranchAddress("ADC_TOF1U",ADC_tof1U);
		fChain->SetBranchAddress("ADC_TOF1D",ADC_tof1D);
		fChain->SetBranchAddress("TDC_TOF1U",TDC_tof1U);
		fChain->SetBranchAddress("TDC_TOF1D",TDC_tof1D);
	//}

	// TOF2
	//if(flag==8){
		fChain->SetBranchAddress("ADC_TOF2AO",ADC_tof2AO);
		fChain->SetBranchAddress("ADC_TOF2AI",ADC_tof2AI);
		fChain->SetBranchAddress("ADC_TOF2BO",ADC_tof2BO);
		fChain->SetBranchAddress("ADC_TOF2BI",ADC_tof2BI);	
		fChain->SetBranchAddress("TDC_TOF2AO",TDC_tof2AO);
		fChain->SetBranchAddress("TDC_TOF2AI",TDC_tof2AI);
		fChain->SetBranchAddress("TDC_TOF2BO",TDC_tof2BO);
		fChain->SetBranchAddress("TDC_TOF2BI",TDC_tof2BI);	
	//}


	///////////////////////////////////////////////////////////////////////////
	
	// SFT
		TH1D *h_ADC_High_SFT[128];	char Title_ADC_High_SFT[128][100];	char Name_ADC_High_SFT[128][100];
		TH1D *h_ADC_Low_SFT[128];   char Title_ADC_Low_SFT[128][100];	char Name_ADC_Low_SFT[128][100];
		TH1D *h_TDC_LE_SFT[128];    char Title_TDC_LE_SFT[128][100];	char Name_TDC_LE_SFT[128][100];     
		TH1D *h_TDC_TE_SFT[128];	char Title_TDC_TE_SFT[128][100];	char Name_TDC_TE_SFT[128][100]; 

		TH1D *h_ADC_High_SFT_cut[128];	 char Title_ADC_High_SFT_cut[128][100];		char Name_ADC_High_SFT_cut[128][100];
		TH1D *h_TDC_LE_SFT_cut[128];     char Title_TDC_LE_SFT_cut[128][100];		char Name_TDC_LE_SFT_cut[128][100];     

	// TARGET
		TH1D *h_ADC_High_TARGET[256];	char Title_ADC_High_TARGET[256][100];	char Name_ADC_High_TARGET[256][100];
		TH1D *h_ADC_Low_TARGET[256];    char Title_ADC_Low_TARGET[256][100];	char Name_ADC_Low_TARGET[256][100];
		TH1D *h_TDC_LE_TARGET[256];     char Title_TDC_LE_TARGET[256][100];		char Name_TDC_LE_TARGET[256][100];     
		TH1D *h_TDC_TE_TARGET[256];		char Title_TDC_TE_TARGET[256][100];		char Name_TDC_TE_TARGET[256][100]; 

		TH1D *h_ADC_High_TARGET_cut[256];	char Title_ADC_High_TARGET_cut[256][100];	char Name_ADC_High_TARGET_cut[256][100];
		TH1D *h_TDC_LE_TARGET_cut[256];     char Title_TDC_LE_TARGET_cut[256][100];		char Name_TDC_LE_TARGET_cut[256][100];     
		
	// TOF1
		TH1D *h_ADC_tof1U[12];   char Title_ADC_tof1U[12][100];	char Name_ADC_tof1U[12][100];
		TH1D *h_ADC_tof1D[12];   char Title_ADC_tof1D[12][100];	char Name_ADC_tof1D[12][100];
		TH1D *h_TDC_tof1U[12];   char Title_TDC_tof1U[12][100];	char Name_TDC_tof1U[12][100];
		TH1D *h_TDC_tof1D[12];   char Title_TDC_tof1D[12][100];	char Name_TDC_tof1D[12][100];

		TH1D *h_ADC_tof1U_cut[12];   char Title_ADC_tof1U_cut[12][100];	char Name_ADC_tof1U_cut[100];
		TH1D *h_ADC_tof1D_cut[12];   char Title_ADC_tof1D_cut[12][100];	char Name_ADC_tof1D_cut[100];
		TH1D *h_TDC_tof1U_cut[12];   char Title_TDC_tof1U_cut[12][100];	char Name_TDC_tof1U_cut[100];
		TH1D *h_TDC_tof1D_cut[12];   char Title_TDC_tof1D_cut[12][100];	char Name_TDC_tof1D_cut[100];
		
	// TOF2
		TH1D *h_ADC_tof2AO[12];   char Title_ADC_tof2AO[12][100];	char Name_ADC_tof2AO[12][100];
		TH1D *h_ADC_tof2AI[12];   char Title_ADC_tof2AI[12][100];	char Name_ADC_tof2AI[12][100];
		TH1D *h_ADC_tof2BO[12];   char Title_ADC_tof2BO[12][100];	char Name_ADC_tof2BO[12][100];
		TH1D *h_ADC_tof2BI[12];   char Title_ADC_tof2BI[12][100];	char Name_ADC_tof2BI[12][100];
		TH1D *h_TDC_tof2AO[12];   char Title_TDC_tof2AO[12][100];	char Name_TDC_tof2AO[12][100];
		TH1D *h_TDC_tof2AI[12];   char Title_TDC_tof2AI[12][100];	char Name_TDC_tof2AI[12][100];
		TH1D *h_TDC_tof2BO[12];   char Title_TDC_tof2BO[12][100];	char Name_TDC_tof2BO[12][100];
		TH1D *h_TDC_tof2BI[12];   char Title_TDC_tof2BI[12][100];	char Name_TDC_tof2BI[12][100];

		TH1D *h_ADC_tof2AO_cut[12];   char Title_ADC_tof2AO_cut[12][100];	char Name_ADC_tof2AO_cut[100];
		TH1D *h_ADC_tof2AI_cut[12];   char Title_ADC_tof2AI_cut[12][100];	char Name_ADC_tof2AI_cut[100];
		TH1D *h_ADC_tof2BO_cut[12];   char Title_ADC_tof2BO_cut[12][100];	char Name_ADC_tof2BO_cut[100];
		TH1D *h_ADC_tof2BI_cut[12];   char Title_ADC_tof2BI_cut[12][100];	char Name_ADC_tof2BI_cut[100];
		TH1D *h_TDC_tof2AO_cut[12];   char Title_TDC_tof2AO_cut[12][100];	char Name_TDC_tof2AO_cut[100];
		TH1D *h_TDC_tof2AI_cut[12];   char Title_TDC_tof2AI_cut[12][100];	char Name_TDC_tof2AI_cut[100];
		TH1D *h_TDC_tof2BO_cut[12];   char Title_TDC_tof2BO_cut[12][100];	char Name_TDC_tof2BO_cut[100];
		TH1D *h_TDC_tof2BI_cut[12];   char Title_TDC_tof2BI_cut[12][100];	char Name_TDC_tof2BI_cut[100];


	/////////////////////////////////////////////////////////////
	// Setting Histogram Names and Titles

	// OLD

	///// SFT
	if(flag==1 || flag==2){
		for(int i=0; i<128; i++){
			sprintf(Title_ADC_High_SFT_cut[i],"Raw ADC High Gain (Ch. %d)  --  SFT",i); 
			sprintf(Title_ADC_Low_SFT[i],"Raw ADC Low Gain (Ch. %d)  --  SFT",i); 
			sprintf(Title_TDC_LE_SFT_cut[i],"Raw TDC (LE) (Ch. %d)  --  SFT",i); 
			sprintf(Title_TDC_TE_SFT[i],"Raw TDC (TE) (Ch. %d)  --  SFT",i); 

			sprintf(Title_ADC_High_SFT[i],"Raw ADC High Gain (Ch. %d)  |  %d #leq TDC #leq %d  --  SFT",i,SFT_TDC_min[i],SFT_TDC_max[i]); 
			sprintf(Title_TDC_LE_SFT[i],"Raw TDC (LE) (Ch. %d)  |  %d #leq TDC #leq %d  --  SFT",i,SFT_TDC_min[i],SFT_TDC_max[i]); 

			sprintf(Name_ADC_High_SFT[i],"ADC_High (Ch. %d)",i);	
			sprintf(Name_ADC_Low_SFT[i],"ADC_Low (Ch. %d)",i);
			sprintf(Name_TDC_LE_SFT[i],"TDC_LE (Ch. %d)",i);
			sprintf(Name_TDC_TE_SFT[i],"TDC_TE (Ch. %d)",i);

			sprintf(Name_ADC_High_SFT_cut[i],"ADC_High (Ch. %d) - SFT",i);	
			sprintf(Name_TDC_LE_SFT_cut[i],"TDC_LE (Ch. %d) - SFT",i);
		}
	}

	///// TARGET
	if(flag==3 || flag==4 || flag==5 || flag==6){
		for(int i=0; i<256; i++){
			sprintf(Title_ADC_High_TARGET_cut[i],"Raw ADC High Gain (Ch. %d)  --  TARGET",i); 
			sprintf(Title_ADC_Low_TARGET[i],"Raw ADC Low Gain (Ch. %d)  --  TARGET",i); 
			sprintf(Title_TDC_LE_TARGET_cut[i],"Raw TDC (LE) (Ch. %d)  --  TARGET",i); 
			sprintf(Title_TDC_TE_TARGET[i],"Raw TDC (TE) (Ch. %d)  --  TARGET",i); 

			sprintf(Title_ADC_High_TARGET[i],"Raw ADC High Gain (Ch. %d)  |  %d #leq TDC #leq %d  --  TARGET",i,TARGET_TDC_min[i],TARGET_TDC_max[i]); 
			sprintf(Title_TDC_LE_TARGET[i],"Raw TDC (LE) (Ch. %d)  |  %d #leq TDC #leq %d  --  TARGET",i,TARGET_TDC_min[i],TARGET_TDC_max[i]); 

			sprintf(Name_ADC_High_TARGET[i],"ADC_High (Ch. %d)",i);	
			sprintf(Name_ADC_Low_TARGET[i],"ADC_Low (Ch. %d)",i);
			sprintf(Name_TDC_LE_TARGET[i],"TDC_LE (Ch. %d)",i);
			sprintf(Name_TDC_TE_TARGET[i],"TDC_TE (Ch. %d)",i);

			sprintf(Name_ADC_High_TARGET_cut[i],"ADC_High (Ch. %d) - TARGET",i);	
			sprintf(Name_TDC_LE_TARGET_cut[i],"TDC_LE (Ch. %d) - TARGET",i);
		}
	}	

	///// TOF1
	if(flag==7){		

		sprintf(Name_ADC_tof1U_cut,"ADC_TOF1 UP");
		sprintf(Name_ADC_tof1D_cut,"ADC_TOF1 DOWN");
		sprintf(Name_TDC_tof1U_cut,"TDC_TOF1 UP");
		sprintf(Name_TDC_tof1D_cut,"TDC_TOF1 DOWN");

		for (int i=0; i<12; i++){
			sprintf(Title_ADC_tof1U[i],"Raw ADC (V792 Ch. %d)  |  %d #leq TDC #leq %d  --  TOF1 (%d UP)",TOF1U_ADC[i],TOF1U_TDC_min[i],TOF1U_TDC_max[i],i+1);
			sprintf(Name_ADC_tof1U[i],"ADC_TOF1 (%d UP)",i+1);
			sprintf(Title_ADC_tof1D[i],"Raw ADC (V792 Ch. %d)  |  %d #leq TDC #leq %d  --  TOF1 (%d DOWN)",TOF1D_ADC[i],TOF1D_TDC_min[i],TOF1D_TDC_max[i],i+1);
			sprintf(Name_ADC_tof1D[i],"ADC_TOF1 (%d DOWN)",i+1);

			sprintf(Title_ADC_tof1U_cut[i],"Raw ADC (V792 Ch. %d)  --  TOF1 (%d UP)",TOF1U_ADC[i],i+1);
			sprintf(Title_ADC_tof1D_cut[i],"Raw ADC (V792 Ch. %d)  --  TOF1 (%d DOWN)",TOF1D_ADC[i],i+1);

			sprintf(Title_TDC_tof1U[i],"Raw TDC (HR_TDC Ch. %d)  |  %d #leq TDC #leq %d  --  TOF1 (%d UP)",TOF1U_TDC[i],TOF1U_TDC_min[i],TOF1U_TDC_max[i],i+1);
			sprintf(Name_TDC_tof1U[i],"TDC_TOF1 (%d UP)",i+1);
			sprintf(Title_TDC_tof1D[i],"Raw TDC (HR _TDC Ch. %d)  |  %d #leq TDC #leq %d  --  TOF1 (%d DOWN)",TOF1D_TDC[i],TOF1D_TDC_min[i],TOF1D_TDC_max[i],i+1);
			sprintf(Name_TDC_tof1D[i],"TDC_TOF1 (%d DOWN)",i+1);

			sprintf(Title_TDC_tof1U_cut[i],"Raw TDC (HR_TDC Ch. %d)  --  TOF1 (%d UP)",TOF1U_TDC[i],i+1);
			sprintf(Title_TDC_tof1D_cut[i],"Raw TDC (HR _TDC Ch. %d)  --  TOF1 (%d DOWN)",TOF1D_TDC[i],i+1);
		}
	}
		

	///// TOF2
	if(flag==8){

		sprintf(Name_ADC_tof2AO_cut,"ADC_TOF2 A OUT");
		sprintf(Name_ADC_tof2AI_cut,"ADC_TOF2 A IN");
		sprintf(Name_ADC_tof2BO_cut,"ADC_TOF2 B OUT");
		sprintf(Name_ADC_tof2BI_cut,"ADC_TOF2 B IN");
		sprintf(Name_TDC_tof2AO_cut,"TDC_TOF2 A OUT");
		sprintf(Name_TDC_tof2AI_cut,"TDC_TOF2 A IN");
		sprintf(Name_TDC_tof2BO_cut,"TDC_TOF2 B OUT");
		sprintf(Name_TDC_tof2BI_cut,"TDC_TOF2 B IN");

		for(int i=0; i<12; i++){
			sprintf(Title_ADC_tof2AO[i],"Raw ADC (V792 Ch. %d)  |  %d #leq TDC #leq %d  --  TOF2 (%d A OUT)",TOF2AO_ADC[i],TOF2AO_TDC_min[i],TOF2AO_TDC_max[i],i+1);
			sprintf(Name_ADC_tof2AO[i],"ADC_TOF2 (%d AO)",i+1);
			sprintf(Title_ADC_tof2AI[i],"Raw ADC (V792 Ch. %d)  |  %d #leq TDC #leq %d  --  TOF2 (%d A IN)",TOF2AI_ADC[i],TOF2AI_TDC_min[i],TOF2AI_TDC_max[i],i+1);
			sprintf(Name_ADC_tof2AI[i],"ADC_TOF2 (%d AI)",i+1);
			sprintf(Title_ADC_tof2BO[i],"Raw ADC (V792 Ch. %d)  |  %d #leq TDC #leq %d  --  TOF2 (%d B OUT)",TOF2BO_ADC[i],TOF2BO_TDC_min[i],TOF2BO_TDC_max[i],i+1);
			sprintf(Name_ADC_tof2BO[i],"ADC_TOF2 (%d BO)",i+1);
			sprintf(Title_ADC_tof2BI[i],"Raw ADC (V792 Ch. %d)  |  %d #leq TDC #leq %d  --  TOF2 (%d B IN)",TOF2BI_ADC[i],TOF2BI_TDC_min[i],TOF2BI_TDC_max[i],i+1);
			sprintf(Name_ADC_tof2BI[i],"ADC_TOF2 (%d BI)",i+1);

			sprintf(Title_ADC_tof2AO_cut[i],"Raw ADC (V792 Ch. %d)  --  TOF2 (%d A OUT)",TOF2AO_ADC[i],i+1);
			sprintf(Title_ADC_tof2AI_cut[i],"Raw ADC (V792 Ch. %d)  --  TOF2 (%d A IN)",TOF2AI_ADC[i],i+1);
			sprintf(Title_ADC_tof2BO_cut[i],"Raw ADC (V792 Ch. %d)  --  TOF2 (%d B OUT)",TOF2BO_ADC[i],i+1);
			sprintf(Title_ADC_tof2BI_cut[i],"Raw ADC (V792 Ch. %d)  --  TOF2 (%d B IN)",TOF2BI_ADC[i],i+1);

			sprintf(Title_TDC_tof2AO[i],"Raw TDC (HR_TDC Ch. %d)  |  %d #leq TDC #leq %d  --  TOF2 (%d A OUT)",TOF2AO_TDC[i],TOF2AO_TDC_min[i],TOF2AO_TDC_max[i],i+1);
			sprintf(Name_TDC_tof2AO[i],"TDC_TOF2 (%d AO)",i+1);
			sprintf(Title_TDC_tof2AI[i],"Raw TDC (HR_TDC Ch. %d)  |  %d #leq TDC #leq %d  --  TOF2 (%d A IN)",TOF2AI_TDC[i],TOF2AI_TDC_min[i],TOF2AI_TDC_max[i],i+1);
			sprintf(Name_TDC_tof2AI[i],"TDC_TOF2 (%d AI)",i+1);
			sprintf(Title_TDC_tof2BO[i],"Raw TDC (HR_TDC Ch. %d)  |  %d #leq TDC #leq %d  --  TOF2 (%d B OUT)",TOF2BO_TDC[i],TOF2BO_TDC_min[i],TOF2BO_TDC_max[i],i+1);
			sprintf(Name_TDC_tof2BO[i],"TDC_TOF2 (%d BO)",i+1);
			sprintf(Title_TDC_tof2BI[i],"Raw TDC (HR_TDC Ch. %d)  |  %d #leq TDC #leq %d  --  TOF2 (%d B IN)",TOF2BI_TDC[i],TOF2BI_TDC_min[i],TOF2BI_TDC_max[i],i+1);
			sprintf(Name_TDC_tof2BI[i],"TDC_TOF2 (%d BI)",i+1);

			sprintf(Title_TDC_tof2AO_cut[i],"Raw TDC (HR_TDC Ch. %d)  --  TOF2 (%d A OUT)",TOF2AO_TDC[i],i+1);
			sprintf(Title_TDC_tof2AI_cut[i],"Raw TDC (HR_TDC Ch. %d)  --  TOF2 (%d A IN)",TOF2AI_TDC[i],i+1);
			sprintf(Title_TDC_tof2BO_cut[i],"Raw TDC (HR_TDC Ch. %d)  --  TOF2 (%d B OUT)",TOF2BO_TDC[i],i+1);
			sprintf(Title_TDC_tof2BI_cut[i],"Raw TDC (HR_TDC Ch. %d)  --  TOF2 (%d B In)",TOF2BI_TDC[i],i+1);

		}
	}

	//////////////////////////////////////////////////////////////////
	// CREATING HISTOGRAMS

	///// SFT	
	if(flag==1 || flag==2){
		for (int i=0; i<128; i++){
			h_ADC_High_SFT[i] = new TH1D(Name_ADC_High_SFT[i],Title_ADC_High_SFT[i],512,0,4096);
			h_ADC_Low_SFT[i] = new TH1D(Name_ADC_Low_SFT[i],Title_ADC_Low_SFT[i],512,0,4096);
			h_TDC_LE_SFT[i] = new TH1D(Name_TDC_LE_SFT[i],Title_TDC_LE_SFT[i],512,0,4096);
			h_TDC_TE_SFT[i] = new TH1D(Name_TDC_TE_SFT[i],Title_TDC_TE_SFT[i],512,0,4096);

			h_ADC_High_SFT_cut[i] = new TH1D(Name_ADC_High_SFT_cut[i],Title_ADC_High_SFT_cut[i],512,0,4096);
			h_TDC_LE_SFT_cut[i] = new TH1D(Name_TDC_LE_SFT_cut[i],Title_TDC_LE_SFT_cut[i],512,0,4096);
		}
	}

	///// TARGET
	if(flag==3 || flag==4 || flag==5 || flag==6){
		for (int i=0; i<256; i++){
			h_ADC_High_TARGET[i] = new TH1D(Name_ADC_High_TARGET[i],Title_ADC_High_TARGET[i],512,0,4096);
			h_ADC_Low_TARGET[i] = new TH1D(Name_ADC_Low_TARGET[i],Title_ADC_Low_TARGET[i],512,0,4096);
			h_TDC_LE_TARGET[i] = new TH1D(Name_TDC_LE_TARGET[i],Title_TDC_LE_TARGET[i],512,0,4096);
			h_TDC_TE_TARGET[i] = new TH1D(Name_TDC_TE_TARGET[i],Title_TDC_TE_TARGET[i],512,0,4096);

			h_ADC_High_TARGET_cut[i] = new TH1D(Name_ADC_High_TARGET_cut[i],Title_ADC_High_TARGET_cut[i],512,0,4096);
			h_TDC_LE_TARGET_cut[i] = new TH1D(Name_TDC_LE_TARGET_cut[i],Title_TDC_LE_TARGET_cut[i],512,0,4096);
		}
	}	

	///// TOF1
	if(flag==7){
		for (int  i=0; i<12; i++) {
			h_ADC_tof1U[i] = new TH1D(Name_ADC_tof1U[i], Title_ADC_tof1U[i], 512, 0, 4096);
			h_ADC_tof1D[i] = new TH1D(Name_ADC_tof1D[i], Title_ADC_tof1D[i], 512, 0, 4096);
			h_TDC_tof1U[i] = new TH1D(Name_TDC_tof1U[i], Title_TDC_tof1U[i],512,0,4096);
			h_TDC_tof1D[i] = new TH1D(Name_TDC_tof1D[i], Title_TDC_tof1D[i],512,0,4096);

			h_ADC_tof1U_cut[i] = new TH1D(Name_ADC_tof1U_cut, Title_ADC_tof1U_cut[i], 512, 0, 4096);
			h_ADC_tof1D_cut[i] = new TH1D(Name_ADC_tof1D_cut, Title_ADC_tof1D_cut[i], 512, 0, 4096);
			h_TDC_tof1U_cut[i] = new TH1D(Name_TDC_tof1U_cut, Title_TDC_tof1U_cut[i],512,0,4096);
			h_TDC_tof1D_cut[i] = new TH1D(Name_TDC_tof1D_cut, Title_TDC_tof1D_cut[i],512,0,4096);
		}
	}

	///// TOF2
	if(flag==8){
		for(int i=0; i<12; i++){
			h_ADC_tof2AO[i] = new TH1D(Name_ADC_tof2AO[i],Title_ADC_tof2AO[i],512,0,4096);
			h_ADC_tof2AI[i] = new TH1D(Name_ADC_tof2AI[i],Title_ADC_tof2AI[i],512,0,4096);
			h_ADC_tof2BO[i] = new TH1D(Name_ADC_tof2BO[i],Title_ADC_tof2BO[i],512,0,4096);
			h_ADC_tof2BI[i] = new TH1D(Name_ADC_tof2BI[i],Title_ADC_tof2BI[i],512,0,4096);
			h_TDC_tof2AO[i] = new TH1D(Name_TDC_tof2AO[i],Title_TDC_tof2AO[i],512,0,4096);
			h_TDC_tof2AI[i] = new TH1D(Name_TDC_tof2AI[i],Title_TDC_tof2AI[i],512,0,4096);
			h_TDC_tof2BO[i] = new TH1D(Name_TDC_tof2BO[i],Title_TDC_tof2BO[i],512,0,4096);
			h_TDC_tof2BI[i] = new TH1D(Name_TDC_tof2BI[i],Title_TDC_tof2BI[i],512,0,4096);

			h_ADC_tof2AO_cut[i] = new TH1D(Name_ADC_tof2AO_cut,Title_ADC_tof2AO_cut[i],512,0,4096);
			h_ADC_tof2AI_cut[i] = new TH1D(Name_ADC_tof2AI_cut,Title_ADC_tof2AI_cut[i],512,0,4096);
			h_ADC_tof2BO_cut[i] = new TH1D(Name_ADC_tof2BO_cut,Title_ADC_tof2BO_cut[i],512,0,4096);
			h_ADC_tof2BI_cut[i] = new TH1D(Name_ADC_tof2BI_cut,Title_ADC_tof2BI_cut[i],512,0,4096);
			h_TDC_tof2AO_cut[i] = new TH1D(Name_TDC_tof2AO_cut,Title_TDC_tof2AO_cut[i],512,0,4096);
			h_TDC_tof2AI_cut[i] = new TH1D(Name_TDC_tof2AI_cut,Title_TDC_tof2AI_cut[i],512,0,4096);
			h_TDC_tof2BO_cut[i] = new TH1D(Name_TDC_tof2BO_cut,Title_TDC_tof2BO_cut[i],512,0,4096);
			h_TDC_tof2BI_cut[i] = new TH1D(Name_TDC_tof2BI_cut,Title_TDC_tof2BI_cut[i],512,0,4096);
		}
	}

	////////////////////////////////////////////////////////////////////////

	int n_counter = 0;
	int n_plus_counter = 0;
	int n_minus_counter = 0;
	int miss_counter = 0;
	int counter_test = 0;

	// READ ALL ENTRIES 
	Int_t nentries = (Int_t)fChain->GetEntries();
	cout <<  "Total Number of Entries :     " <<  nentries << endl;

	cout << "  " << endl;
	cout << "***************  Processing....  ***************" << endl;
	cout << "  " << endl;

	if(nevt!=0) nentries=nevt;
	for (Int_t i=0; i<nentries; i++) {
	//for (Int_t i=0; i<2000; i++) {
		fChain->GetEntry(i);

		if(i%10000==1) cout<<"**** "<<i<<" events done"<<endl;

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///// GAP SCORING
		
		high_gap_hit = 0;
  		gap_to_fit = 0;
  		gap_TOF2 = 0;
  		//counter_test = 0;

		for(int j=0; j<12; j++) gap_counter[j] = 0;

	    for(int j=0; j<12; j++){
        //	cout << "XXXXXX   " << j << "    " << ADC_tof1U[j] << endl;
     		if(ADC_tof1U[j] >= TOF1U_ADC_Thr[j]){ 
        		gap_counter[j]++;
        //		cout << " ####   " << i << "    " <<  j << "   " << gap_counter[j] << endl;
 		    }
	    }

    	for(int j=0; j<12; j++){
      		if(ADC_tof1D[j] >= TOF1D_ADC_Thr[j]){
        		gap_counter[j]++;
		    }
	    }

    	for(int j=0; j<12; j++){
    		if(TDC_tof1U[j] >= TOF1U_TDC_min[j] && TDC_tof1U[j] <= TOF1U_TDC_max[j]){
        		gap_counter[j]++; 
      		}
    	}

    	for(int j=0; j<12; j++){
      		if(TDC_tof1D[j] >= TOF1D_TDC_min[j] && TDC_tof1D[j] <= TOF1D_TDC_max[j]){
        		gap_counter[j]++;
     		}
    	}


    	for(int j=0; j<12; j++){
      		if(j!=0 && j!=11 && gap_counter[j]>0){

        		if(ADC_tof2AO[j] >= TOF2AO_ADC_Thr[j]) gap_counter[j]++;
        		if(ADC_tof2AO[j-1] >= TOF2AO_ADC_Thr[j-1]) gap_counter[j]++;
        		if(ADC_tof2AO[j+1] >= TOF2AO_ADC_Thr[j+1]) gap_counter[j]++;
 
	    	    if(ADC_tof2AI[j] >= TOF2AI_ADC_Thr[j]) gap_counter[j]++;
    		    if(ADC_tof2AI[j-1] >= TOF2AI_ADC_Thr[j-1]) gap_counter[j]++;
        		if(ADC_tof2AI[j+1] >= TOF2AI_ADC_Thr[j+1]) gap_counter[j]++;   

        		if(ADC_tof2BO[j] >= TOF2BO_ADC_Thr[j]) gap_counter[j]++;
        		if(ADC_tof2BO[j-1] >= TOF2BO_ADC_Thr[j-1]) gap_counter[j]++;
        		if(ADC_tof2BO[j+1] >= TOF2BO_ADC_Thr[j+1]) gap_counter[j]++;   

        		if(ADC_tof2BI[j] >= TOF2BI_ADC_Thr[j]) gap_counter[j]++;
        		if(ADC_tof2BI[j-1] >= TOF2BI_ADC_Thr[j-1]) gap_counter[j]++;
        		if(ADC_tof2BI[j+1] >= TOF2BI_ADC_Thr[j+1]) gap_counter[j]++;

        		if(TDC_tof2AO[j] >= TOF2AO_TDC_min[j] && TDC_tof2AO[j] <= TOF2AO_TDC_max[j]) gap_counter[j]++;   
        		if(TDC_tof2AO[j-1] >= TOF2AO_TDC_min[j-1] && TDC_tof2AO[j-1] <= TOF2AO_TDC_max[j-1]) gap_counter[j]++;   
        		if(TDC_tof2AO[j+1] >= TOF2AO_TDC_min[j+1] && TDC_tof2AO[j+1] <= TOF2AO_TDC_max[j+1]) gap_counter[j]++;   

        		if(TDC_tof2AI[j] >= TOF2AI_TDC_min[j] && TDC_tof2AI[j] <= TOF2AI_TDC_max[j]) gap_counter[j]++;   
        		if(TDC_tof2AI[j-1] >= TOF2AI_TDC_min[j-1] && TDC_tof2AI[j-1] <= TOF2AI_TDC_max[j-1]) gap_counter[j]++;   
        		if(TDC_tof2AI[j+1] >= TOF2AI_TDC_min[j+1] && TDC_tof2AI[j+1] <= TOF2AI_TDC_max[j+1]) gap_counter[j]++;   

        		if(TDC_tof2BO[j] >= TOF2BO_TDC_min[j] && TDC_tof2BO[j] <= TOF2BO_TDC_max[j]) gap_counter[j]++;   
        		if(TDC_tof2BO[j-1] >= TOF2BO_TDC_min[j-1] && TDC_tof2BO[j-1] <= TOF2BO_TDC_max[j-1]) gap_counter[j]++;   
        		if(TDC_tof2BO[j+1] >= TOF2BO_TDC_min[j+1] && TDC_tof2BO[j+1] <= TOF2BO_TDC_max[j+1]) gap_counter[j]++;   

        		if(TDC_tof2BI[j] >= TOF2BI_TDC_min[j] && TDC_tof2BI[j] <= TOF2BI_TDC_max[j]) gap_counter[j]++;   
        		if(TDC_tof2BI[j-1] >= TOF2BI_TDC_min[j-1] && TDC_tof2BI[j-1] <= TOF2BI_TDC_max[j-1]) gap_counter[j]++;   
        		if(TDC_tof2BI[j+1] >= TOF2BI_TDC_min[j+1] && TDC_tof2BI[j+1] <= TOF2BI_TDC_max[j+1]) gap_counter[j]++;   
      		}
    	}

    	if(gap_counter[0]>0){
      		if(ADC_tof2AO[0] >= TOF2AO_ADC_Thr[0]) gap_counter[0]++;
      		if(ADC_tof2AO[1] >= TOF2AO_ADC_Thr[1]) gap_counter[0]++;
      		if(ADC_tof2AO[11] >= TOF2AO_ADC_Thr[11]) gap_counter[0]++;

      		if(ADC_tof2AI[0] >= TOF2AI_ADC_Thr[0]) gap_counter[0]++;
      		if(ADC_tof2AI[1] >= TOF2AI_ADC_Thr[1]) gap_counter[0]++;
      		if(ADC_tof2AI[11] >= TOF2AI_ADC_Thr[11]) gap_counter[0]++;   

      		if(ADC_tof2BO[0] >= TOF2BO_ADC_Thr[0]) gap_counter[0]++;
      		if(ADC_tof2BO[1] >= TOF2BO_ADC_Thr[1]) gap_counter[0]++;
      		if(ADC_tof2BO[11] >= TOF2BO_ADC_Thr[11]) gap_counter[0]++;   

      		if(ADC_tof2BI[0] >= TOF2BI_ADC_Thr[0]) gap_counter[0]++;
      		if(ADC_tof2BI[1] >= TOF2BI_ADC_Thr[1]) gap_counter[0]++;
      		if(ADC_tof2BI[11] >= TOF2BI_ADC_Thr[11]) gap_counter[0]++;

      		if(TDC_tof2AO[0] >= TOF2AO_TDC_min[0] && TDC_tof2AO[0] <= TOF2AO_TDC_max[0]) gap_counter[0]++;   
      		if(TDC_tof2AO[1] >= TOF2AO_TDC_min[1] && TDC_tof2AO[1] <= TOF2AO_TDC_max[1]) gap_counter[0]++;   
      		if(TDC_tof2AO[11] >= TOF2AO_TDC_min[11] && TDC_tof2AO[11] <= TOF2AO_TDC_max[11]) gap_counter[0]++;   

      		if(TDC_tof2AI[0]>= TOF2AI_TDC_min[0] && TDC_tof2AI[0] <= TOF2AI_TDC_max[0]) gap_counter[0]++;   
      		if(TDC_tof2AI[1]>= TOF2AI_TDC_min[1] && TDC_tof2AI[1] <= TOF2AI_TDC_max[1]) gap_counter[0]++;   
      		if(TDC_tof2AI[11]>= TOF2AI_TDC_min[11] && TDC_tof2AI[11] <= TOF2AI_TDC_max[11]) gap_counter[0]++;   

      		if(TDC_tof2BO[0]>= TOF2BO_TDC_min[0] && TDC_tof2BO[0] <= TOF2BO_TDC_max[0]) gap_counter[0]++;   
      		if(TDC_tof2BO[1]>= TOF2BO_TDC_min[1] && TDC_tof2BO[1] <= TOF2BO_TDC_max[1]) gap_counter[0]++;   
      		if(TDC_tof2BO[11]>= TOF2BO_TDC_min[11] && TDC_tof2BO[11] <= TOF2BO_TDC_max[11]) gap_counter[0]++;   

      		if(TDC_tof2BI[0]>= TOF2BI_TDC_min[0] && TDC_tof2BI[0] <= TOF2BI_TDC_max[0]) gap_counter[0]++;   
      		if(TDC_tof2BI[1]>= TOF2BI_TDC_min[1] && TDC_tof2BI[1] <= TOF2BI_TDC_max[1]) gap_counter[0]++;   
      		if(TDC_tof2BI[11]>= TOF2BI_TDC_min[11] && TDC_tof2BI[11] <= TOF2BI_TDC_max[11]) gap_counter[0]++;   
    	}

    	if(gap_counter[11]>0){
      		if(ADC_tof2AO[11] >= TOF2AO_ADC_Thr[11]) gap_counter[11]++;
      		if(ADC_tof2AO[10] >= TOF2AO_ADC_Thr[10]) gap_counter[11]++;
      		if(ADC_tof2AO[0] >= TOF2AO_ADC_Thr[0]) gap_counter[11]++;

      		if(ADC_tof2AI[11] >= TOF2AI_ADC_Thr[11]) gap_counter[11]++;
      		if(ADC_tof2AI[10] >= TOF2AI_ADC_Thr[10]) gap_counter[11]++;
      		if(ADC_tof2AI[0] >= TOF2AI_ADC_Thr[0]) gap_counter[11]++;   

      		if(ADC_tof2BO[11] >= TOF2BO_ADC_Thr[11]) gap_counter[11]++;
      		if(ADC_tof2BO[10] >= TOF2BO_ADC_Thr[10]) gap_counter[11]++;
      		if(ADC_tof2BO[0] >= TOF2BO_ADC_Thr[0]) gap_counter[11]++;   

      		if(ADC_tof2BI[11] >= TOF2BI_ADC_Thr[11]) gap_counter[11]++;
      		if(ADC_tof2BI[10] >= TOF2BI_ADC_Thr[10]) gap_counter[11]++;
      		if(ADC_tof2BI[0] >= TOF2BI_ADC_Thr[0]) gap_counter[11]++;

      		if(TDC_tof2AO[11] >= TOF2AO_TDC_min[11] && TDC_tof2AO[11] <= TOF2AO_TDC_max[11]) gap_counter[11]++;   
      		if(TDC_tof2AO[10] >= TOF2AO_TDC_min[10] && TDC_tof2AO[10] <= TOF2AO_TDC_max[10]) gap_counter[11]++;   
      		if(TDC_tof2AO[0] >= TOF2AO_TDC_min[0] && TDC_tof2AO[0] <= TOF2AO_TDC_max[0]) gap_counter[11]++;   

      		if(TDC_tof2AI[11]>= TOF2AI_TDC_min[11] && TDC_tof2AI[11] <= TOF2AI_TDC_max[11]) gap_counter[11]++;   
      		if(TDC_tof2AI[10]>= TOF2AI_TDC_min[10] && TDC_tof2AI[10] <= TOF2AI_TDC_max[10]) gap_counter[11]++;   
      		if(TDC_tof2AI[0]>= TOF2AI_TDC_min[0] && TDC_tof2AI[0] <= TOF2AI_TDC_max[0]) gap_counter[11]++;   

      		if(TDC_tof2BO[11]>= TOF2BO_TDC_min[11] && TDC_tof2BO[11] <= TOF2BO_TDC_max[11]) gap_counter[11]++;   
      		if(TDC_tof2BO[10]>= TOF2BO_TDC_min[10] && TDC_tof2BO[10] <= TOF2BO_TDC_max[10]) gap_counter[11]++;   
      		if(TDC_tof2BO[0]>= TOF2BO_TDC_min[0] && TDC_tof2BO[0] <= TOF2BO_TDC_max[0]) gap_counter[11]++;   

      		if(TDC_tof2BI[11]>= TOF2BI_TDC_min[11] && TDC_tof2BI[11] <= TOF2BI_TDC_max[11]) gap_counter[11]++;   
      		if(TDC_tof2BI[10]>= TOF2BI_TDC_min[10] && TDC_tof2BI[10] <= TOF2BI_TDC_max[10]) gap_counter[11]++;   
      		if(TDC_tof2BI[0]>= TOF2BI_TDC_min[0] && TDC_tof2BI[0] <= TOF2BI_TDC_max[0]) gap_counter[11]++;   
    	}



  		for (int k=0; k<12; k++) {
    		if (gap_counter[k] >= high_gap_hit) {
      			high_gap_hit = gap_counter[k];
			    gap_to_fit = k+1;
      		}
  		}
    
    	for(int j=0; j<12; j++){
 		    if((TDC_tof2AO[j] >= TOF2AO_TDC_min[j] && TDC_tof2AO[j] <= TOF2AO_TDC_max[j]) ||
 		       (TDC_tof2AI[j] >= TOF2AI_TDC_min[j] && TDC_tof2AI[j] <= TOF2AI_TDC_max[j]) ||
 		       (TDC_tof2BO[j] >= TOF2BO_TDC_min[j] && TDC_tof2BO[j] <= TOF2BO_TDC_max[j]) ||
 		       (TDC_tof2BI[j] >= TOF2BI_TDC_min[j] && TDC_tof2BI[j] <= TOF2BI_TDC_max[j])){

 		       		gap_TOF2 = j+1;
			}
 		}

 		//if(gap_TOF2==gap_selected){
 		//	counter_test++;
 		//	cout << " ######  " << i << "    " << gap_to_fit << "   " << gap_TOF2 << "     " << gap_TOF2-gap_to_fit << "   " << counter_test << endl; 
 		//}

 		if(gap_TOF2-gap_to_fit==0) n_counter++;
 		else if(gap_TOF2-gap_to_fit==-1) n_plus_counter++;
 		else if(gap_TOF2-gap_to_fit==1)	n_minus_counter++;
 		else if((gap_TOF2-gap_to_fit==-11) && gap_TOF2==1 && gap_to_fit==12){
 			n_minus_counter++;
 		}
 		else if((gap_TOF2-gap_to_fit==-11) && gap_TOF2==1 && gap_to_fit==12){
 			n_minus_counter++;
 		}
 		else if((gap_TOF2-gap_to_fit==11) && gap_TOF2==12 && gap_to_fit==1){
 			n_plus_counter++;
 		}
 		else miss_counter++;

  		//cout << i << "      " << gap_to_fit << "     "  << gap_TOF2 << "    " << gap_TOF2 - gap_to_fit << endl;
  		//cout << n_counter << "    " << n_plus_counter << "    " << n_minus_counter << "    " << miss_counter << endl;
    //	cout << "#####  " << i << "    " << gap_to_fit << endl;

  	//	for(int j=0; j<12; j++){
    //		cout << "#####  " << i << "    " << j+1 << "    " << gap_counter[j] << "   " << gap_to_fit << endl;
    //	}
 
    //for(int i=0; i<12; i++)    !!!! TO CHECK !!!
    //{
    //  if(Good_tof1[i] && Good_tof2[i]) gap_counter[i]++;
    //}
  //}


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		///// SFT
		if(flag==1 || flag==2){
			for (Int_t j=0; j<128; j++){
				h_ADC_High_SFT[j]->Fill(ADC_High_sft[j]);
				h_ADC_Low_SFT[j]->Fill(ADC_Low_sft[j]);
				h_TDC_LE_SFT[j]->Fill(TDC_LE_sft[j][0]);
				h_TDC_TE_SFT[j]->Fill(TDC_TE_sft[j][0]);

				if(TDC_LE_sft[j][0]>=SFT_TDC_min[j] && TDC_LE_sft[j][0]<=SFT_TDC_max[j]){
					h_TDC_LE_SFT_cut[j]->Fill(TDC_LE_sft[j][0]);
				}

				if((TDC_LE_sft[j][0]>=SFT_TDC_min[j] && TDC_LE_sft[j][0]<=SFT_TDC_max[j]) ||
					(TDC_LE_sft[j][1]>=SFT_TDC_min[j] && TDC_LE_sft[j][1]<=SFT_TDC_max[j]) ||
					(TDC_LE_sft[j][2]>=SFT_TDC_min[j] && TDC_LE_sft[j][2]<=SFT_TDC_max[j]) ||
					(TDC_LE_sft[j][3]>=SFT_TDC_min[j] && TDC_LE_sft[j][3]<=SFT_TDC_max[j])){
					h_ADC_High_SFT_cut[j]->Fill(ADC_High_sft[j]);
				}
			}
		}

		///// TARGET
		if(flag==3 || flag==4 || flag==5 || flag==6){

			//if(gap_to_fit == gap_selected){
			//	list_counter++;
//				myfile << "Writing this to a file.\n";
			//	myfile << "Writing this to a file.\n";
			//	Event_List << list_counter << "   " << i << "\n";
				//cout << "#############  " << gap_selected << "     " << i << endl;	
			//}				


			for (Int_t j=0; j<256; j++){
				if(gap_selected == 0){
					h_ADC_High_TARGET[j]->Fill(adc_high_target[j]);
					h_ADC_Low_TARGET[j]->Fill(adc_low_target[j]);
					h_TDC_LE_TARGET[j]->Fill(tdc_le_target[j][0]);
					h_TDC_TE_TARGET[j]->Fill(tdc_te_target[j][0]);

					if(tdc_le_target[j][0]>=TARGET_TDC_min[j] && tdc_le_target[j][0]<=TARGET_TDC_max[j]){
						h_TDC_LE_TARGET_cut[j]->Fill(tdc_le_target[j][0]);
					}

					if((tdc_le_target[j][0]>=TARGET_TDC_min[j] && tdc_le_target[j][0]<=TARGET_TDC_max[j]) ||
						(tdc_le_target[j][1]>=TARGET_TDC_min[j] && tdc_le_target[j][1]<=TARGET_TDC_max[j]) ||
						(tdc_le_target[j][2]>=TARGET_TDC_min[j] && tdc_le_target[j][2]<=TARGET_TDC_max[j]) ||
						(tdc_le_target[j][3]>=TARGET_TDC_min[j] && tdc_le_target[j][3]<=TARGET_TDC_max[j])){
						h_ADC_High_TARGET_cut[j]->Fill(adc_high_target[j]);
					}
				}

				else if(gap_to_fit == gap_selected){
					h_ADC_High_TARGET[j]->Fill(adc_high_target[j]);
					h_ADC_Low_TARGET[j]->Fill(adc_low_target[j]);
					h_TDC_LE_TARGET[j]->Fill(tdc_le_target[j][0]);
					h_TDC_TE_TARGET[j]->Fill(tdc_te_target[j][0]);

					if(tdc_le_target[j][0]>=TARGET_TDC_min[j] && tdc_le_target[j][0]<=TARGET_TDC_max[j]){
						h_TDC_LE_TARGET_cut[j]->Fill(tdc_le_target[j][0]);
					}

					if((tdc_le_target[j][0]>=TARGET_TDC_min[j] && tdc_le_target[j][0]<=TARGET_TDC_max[j]) ||
						(tdc_le_target[j][1]>=TARGET_TDC_min[j] && tdc_le_target[j][1]<=TARGET_TDC_max[j]) ||
						(tdc_le_target[j][2]>=TARGET_TDC_min[j] && tdc_le_target[j][2]<=TARGET_TDC_max[j]) ||
						(tdc_le_target[j][3]>=TARGET_TDC_min[j] && tdc_le_target[j][3]<=TARGET_TDC_max[j])){
						h_ADC_High_TARGET_cut[j]->Fill(adc_high_target[j]);
					}
				}
			}
		}	

		///// TOF1
		if(flag==7){
			for (Int_t j=0; j<12; j++){
				h_TDC_tof1U[j]->Fill(TDC_tof1U[j]);
				h_TDC_tof1D[j]->Fill(TDC_tof1D[j]);
				h_ADC_tof1U[j]->Fill(ADC_tof1U[j]);
				h_ADC_tof1D[j]->Fill(ADC_tof1D[j]);

				if(TDC_tof1U[j]>=TOF1U_TDC_min[j] && TDC_tof1U[j]<=TOF1U_TDC_max[j]){
					h_ADC_tof1U_cut[j]->Fill(ADC_tof1U[j]);
					h_TDC_tof1U_cut[j]->Fill(TDC_tof1U[j]);
				}

				if(TDC_tof1D[j]>=TOF1D_TDC_min[j] && TDC_tof1D[j]<=TOF1D_TDC_max[j]){
					h_ADC_tof1D_cut[j]->Fill(ADC_tof1D[j]);
					h_TDC_tof1D_cut[j]->Fill(TDC_tof1D[j]);
				}
			}
		}
	
		///// TOF2
		if(flag==8){
			for(Int_t j=0; j<12; j++){
				h_ADC_tof2AO[j]->Fill(ADC_tof2AO[j]);
				h_ADC_tof2AI[j]->Fill(ADC_tof2AI[j]);
				h_ADC_tof2BO[j]->Fill(ADC_tof2BO[j]);
				h_ADC_tof2BI[j]->Fill(ADC_tof2BI[j]);
				h_TDC_tof2AO[j]->Fill(TDC_tof2AO[j]);
				h_TDC_tof2AI[j]->Fill(TDC_tof2AI[j]);
				h_TDC_tof2BO[j]->Fill(TDC_tof2BO[j]);
				h_TDC_tof2BI[j]->Fill(TDC_tof2BI[j]);

				if(TDC_tof2AO[j]>=TOF2AO_TDC_min[j] && TDC_tof2AO[j]<=TOF2AO_TDC_max[j]){
					h_ADC_tof2AO_cut[j]->Fill(ADC_tof2AO[j]);
					h_TDC_tof2AO_cut[j]->Fill(TDC_tof2AO[j]);		
				}

				if(TDC_tof2AI[j]>=TOF2AI_TDC_min[j] && TDC_tof2AI[j]<=TOF2AI_TDC_max[j]){
					h_ADC_tof2AI_cut[j]->Fill(ADC_tof2AI[j]);
					h_TDC_tof2AI_cut[j]->Fill(TDC_tof2AI[j]);		
				}

				if(TDC_tof2BO[j]>=TOF2BO_TDC_min[j] && TDC_tof2BO[j]<=TOF2BO_TDC_max[j]){
					h_ADC_tof2BO_cut[j]->Fill(ADC_tof2BO[j]);
					h_TDC_tof2BO_cut[j]->Fill(TDC_tof2BO[j]);		
				}

				if(TDC_tof2BI[j]>=TOF2BI_TDC_min[j] && TDC_tof2BI[j]<=TOF2BI_TDC_max[j]){
					h_ADC_tof2BI_cut[j]->Fill(ADC_tof2BI[j]);
					h_TDC_tof2BI_cut[j]->Fill(TDC_tof2BI[j]);		
				}
			}
		}


	} // EndLoop over Events

	Event_List.close();

	//////////////////////////////////////////////////////////////////////////////////////////////////
	
	cout << "  " << endl;
	cout << n_counter << " ( " << 100*(double(n_counter)/double(nentries)) << "\% ) " << "      " 
		 << n_plus_counter << " ( " << 100*(double(n_plus_counter)/double(nentries)) << "\% ) " <<  "      " 
		 << n_minus_counter << " ( " << 100*(double(n_minus_counter)/double(nentries)) << "\% ) " << "      " 
		 << miss_counter << " ( " << 100*(double(miss_counter)/double(nentries)) << "\% ) " << endl;
	cout << "  " << endl;
	

	// CREATING CANEVAS
	
	// SFT & TARGET
	char Name_Can_ADC_High[100];			char Title_Can_ADC_High[100];
	char Name_Can_ADC_High2[100];			char Title_Can_ADC_High2[100];
	char Name_Can_ADC_Low[100];				char Title_Can_ADC_Low[100];
	char Name_Can_ADC_Low2[100];			char Title_Can_ADC_Low2[100];
	char Name_Can_TDC_LE[100];				char Title_Can_TDC_LE[100];
	char Name_Can_TDC_LE2[100];				char Title_Can_TDC_LE2[100];
	char Name_Can_TDC_TE[100];				char Title_Can_TDC_TE[100];
	char Name_Can_TDC_TE2[100];				char Title_Can_TDC_TE2[100];

	// TOF1
	char Name_Can_TOF1[100];			char Title_Can_TOF1[100];

	// TOF2
	char Name_Can_TOF2_1[100];			char Title_Can_TOF2_1[100];
	char Name_Can_TOF2_2[100];			char Title_Can_TOF2_2[100];


	// SFT
	if(flag==1){
		sprintf(Name_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  | EASIROC  Board %d",run_number,flag);
		sprintf(Name_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  | EASIROC  Board %d",run_number,flag);

		sprintf(Name_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (SFT)  |  EASIROC Board %d",run_number,flag);
	}
	
	if(flag==2){
		sprintf(Name_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95 (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95 (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 64 - 95  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 64 - 95  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (SFT)  |  EASIROC Board %d",run_number,flag);
	}

	// TARGET
	if(flag==3){
		sprintf(Name_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 32  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 32  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 0 - 31  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 32 - 63  (TARGET)  |  EASIROC Board %d",run_number,flag);
	}

	if(flag==4){
		sprintf(Name_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 64 - 95  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 96 - 127  (TARGET)  |  EASIROC Board %d",run_number,flag);
	}


	if(flag==5){
		sprintf(Name_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 128 - 159  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 160 - 191  (TARGET)  |  EASIROC Board %d",run_number,flag);
	}

	if(flag==6){
		sprintf(Name_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_High,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_High2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_Low,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_Low2,"EASIROC ADC(LG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number,flag);
					
		sprintf(Name_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_LE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_LE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_TE,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 192 - 223  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Name_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_TE2,"EASIROC ADC(HG) & TDC(LE) -- Run %d  Ch. 224 - 255  (TARGET)  |  EASIROC Board %d",run_number,flag);
	}

	// TOF1
	if(flag==7){
		sprintf(Name_Can_TOF1,"TOF1 ADCs & TDCs -- Run %d  (Gap 1 - 12)",run_number);
		sprintf(Title_Can_TOF1,"TOF1 ADCs & TDCs -- Run %d  (Gap 1 - 12)",run_number);
	}

	// TOF2
	if(flag==8){
		sprintf(Name_Can_TOF2_1,"TOF2 ADCs & TDCs -- Run %d  (Gap 1 - 6)",run_number);
		sprintf(Title_Can_TOF2_1,"TOF2 ADCs & TDCs -- Run %d  (Gap 1 -6)",run_number);
		sprintf(Name_Can_TOF2_2,"TOF2 ADCs & TDCs -- Run %d  (Gap 7 - 12)",run_number);
		sprintf(Title_Can_TOF2_2,"TOF2 ADCs & TDCs -- Run %d  (Gap 7 - 12)",run_number);
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// DRAWING
	Int_t TDC_Graph_ymax = -1;
	Int_t TDC_Graph_xmax = 1500;

	Int_t TOF1_Graph_ymax = -1;	
	Int_t TOF1_Graph_xmax = 3000;

	Int_t TOF2_Graph_ymax = -1;
	Int_t TOF2_Graph_xmax = 3000;

	// SFT
	if(flag==1 || flag==2){
		for(int i=0; i<128; i++){
			h_ADC_High_SFT[i]->SetFillColor(4);
			h_ADC_High_SFT_cut[i]->SetFillColor(2);

			h_TDC_LE_SFT_cut[i]->SetFillColor(2);
		}
	}

	// TARGET
	if(flag==3 || flag==4 || flag==5 || flag==6){
		for(int i=0; i<256; i++){
			h_ADC_High_TARGET[i]->SetFillColor(4);
			h_ADC_High_TARGET_cut[i]->SetFillColor(2);

			h_TDC_LE_TARGET_cut[i]->SetFillColor(2);
		}
	}

	// TOF1
	if(flag==7){
		for(int i=0; i<12; i++){
			h_ADC_tof1U[i]->SetFillColor(4);
			h_ADC_tof1D[i]->SetFillColor(4);
			h_ADC_tof1U_cut[i]->SetFillColor(2);
			h_ADC_tof1D_cut[i]->SetFillColor(2);

			h_TDC_tof1U_cut[i]->SetFillColor(2);
			h_TDC_tof1D_cut[i]->SetFillColor(2);
		}
	}

	// TOF2
	if(flag==8){
		for(int i=0; i<12; i++){
			h_ADC_tof2AO[i]->SetFillColor(4);
			h_ADC_tof2AI[i]->SetFillColor(4);
			h_ADC_tof2BO[i]->SetFillColor(4);
			h_ADC_tof2BI[i]->SetFillColor(4);
			h_ADC_tof2AO_cut[i]->SetFillColor(2);
			h_ADC_tof2AI_cut[i]->SetFillColor(2);
			h_ADC_tof2BO_cut[i]->SetFillColor(2);
			h_ADC_tof2BI_cut[i]->SetFillColor(2);

			h_TDC_tof2AO_cut[i]->SetFillColor(2);
			h_TDC_tof2AI_cut[i]->SetFillColor(2);
			h_TDC_tof2BO_cut[i]->SetFillColor(2);
			h_TDC_tof2BI_cut[i]->SetFillColor(2);
		}
	}



	// SFT
	if(flag==1){
		TCanvas *c_SFT2;
		c_SFT2 = new TCanvas(Name_Can_ADC_High2,Title_Can_ADC_High2,1200,500); 
		c_SFT2->Divide(8,8);
		c_SFT2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<32; ican++){
			c_SFT2->cd(ican+1)->SetLogy();  
			h_ADC_High_SFT[ican+32]->Draw();
			h_ADC_High_SFT_cut[ican+32]->Draw("same");
			c_SFT2->cd(ican+1+32);
			h_TDC_LE_SFT[ican+32]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			h_TDC_LE_SFT_cut[ican+32]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			if(TDC_Graph_ymax > 0) {  
				h_TDC_LE_SFT[ican+32]->SetAxisRange(0, TDC_Graph_ymax,"Y");
				h_TDC_LE_SFT_cut[ican+32]->SetAxisRange(0, TDC_Graph_ymax,"Y");
			}	
			h_TDC_LE_SFT[ican+32]->Draw();
			h_TDC_LE_SFT_cut[ican+32]->Draw("same");
		}


		TCanvas *c_SFT;
		c_SFT = new TCanvas(Name_Can_TDC_LE,Title_Can_TDC_LE,0,200,1200,500); 
		c_SFT->Divide(8,8);
		c_SFT->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<32; ican++){
			c_SFT->cd(ican+1)->SetLogy();
			h_ADC_High_SFT[ican]->Draw();
			h_ADC_High_SFT_cut[ican]->Draw("same");
			c_SFT->cd(ican+1+32);
			h_TDC_LE_SFT[ican]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			h_TDC_LE_SFT_cut[ican]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			if(TDC_Graph_ymax > 0) {  
				h_TDC_LE_SFT[ican]->SetAxisRange(0, TDC_Graph_ymax,"Y");
				h_TDC_LE_SFT_cut[ican]->SetAxisRange(0, TDC_Graph_ymax,"Y");
			}
			h_TDC_LE_SFT[ican]->Draw();
			h_TDC_LE_SFT_cut[ican]->Draw("same");
		}
	}
	
	if(flag==2){
		TCanvas *c_SFT2;
		c_SFT2 = new TCanvas(Name_Can_ADC_High2,Title_Can_ADC_High2,1200,500); 
		c_SFT2->Divide(8,8);
		c_SFT2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<32; ican++){
			c_SFT2->cd(ican+1)->SetLogy();  
			h_ADC_High_SFT[ican+96]->Draw();
			h_ADC_High_SFT_cut[ican+96]->Draw("same");
			c_SFT2->cd(ican+1+32);
			h_TDC_LE_SFT[ican+96]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			h_TDC_LE_SFT_cut[ican+96]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			if(TDC_Graph_ymax > 0) {  
				h_TDC_LE_SFT[ican+96]->SetAxisRange(0, TDC_Graph_ymax,"Y");
				h_TDC_LE_SFT_cut[ican+96]->SetAxisRange(0, TDC_Graph_ymax,"Y");
			}	
			h_TDC_LE_SFT[ican+96]->Draw();
			h_TDC_LE_SFT_cut[ican+96]->Draw("same");
		}


		TCanvas *c_SFT;
		c_SFT = new TCanvas(Name_Can_TDC_LE,Title_Can_TDC_LE,0,200,1200,500); 
		c_SFT->Divide(8,8);
		c_SFT->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<32; ican++){
			c_SFT->cd(ican+1)->SetLogy();
			h_ADC_High_SFT[ican+64]->Draw();
			h_ADC_High_SFT_cut[ican+64]->Draw("same");
			c_SFT->cd(ican+1+32);
			h_TDC_LE_SFT[ican+64]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			h_TDC_LE_SFT_cut[ican+64]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			if(TDC_Graph_ymax > 0) {  
				h_TDC_LE_SFT[ican+64]->SetAxisRange(0, TDC_Graph_ymax,"Y");
				h_TDC_LE_SFT_cut[ican+64]->SetAxisRange(0, TDC_Graph_ymax,"Y");
			}
			h_TDC_LE_SFT[ican+64]->Draw();
			h_TDC_LE_SFT_cut[ican+64]->Draw("same");
		}
	}

	if(flag==3){
		TCanvas *c_TARGET2;
		c_TARGET2 = new TCanvas(Name_Can_ADC_High2,Title_Can_ADC_High2,1200,500); 
		c_TARGET2->Divide(8,8);
		c_TARGET2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<32; ican++){
			c_TARGET2->cd(ican+1)->SetLogy();  
			h_ADC_High_TARGET[ican+32]->Draw();
			h_ADC_High_TARGET_cut[ican+32]->Draw("same");
			c_TARGET2->cd(ican+1+32);
			h_TDC_LE_TARGET[ican+32]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			h_TDC_LE_TARGET_cut[ican+32]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			if(TDC_Graph_ymax > 0) {  
				h_TDC_LE_TARGET[ican+32]->SetAxisRange(0, TDC_Graph_ymax,"Y");
				h_TDC_LE_TARGET_cut[ican+32]->SetAxisRange(0, TDC_Graph_ymax,"Y");
			}	
			h_TDC_LE_TARGET[ican+32]->Draw();
			h_TDC_LE_TARGET_cut[ican+32]->Draw("same");
		}


		TCanvas *c_TARGET;
		c_TARGET = new TCanvas(Name_Can_TDC_LE,Title_Can_TDC_LE,0,200,1200,500); 
		c_TARGET->Divide(8,8);
		c_TARGET->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<32; ican++){
			c_TARGET->cd(ican+1)->SetLogy();
			h_ADC_High_TARGET[ican]->Draw();
			h_ADC_High_TARGET_cut[ican]->Draw("same");
			c_TARGET->cd(ican+1+32);
			h_TDC_LE_TARGET[ican]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			h_TDC_LE_TARGET_cut[ican]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			if(TDC_Graph_ymax > 0) {  
				h_TDC_LE_TARGET[ican]->SetAxisRange(0, TDC_Graph_ymax,"Y");
				h_TDC_LE_TARGET_cut[ican]->SetAxisRange(0, TDC_Graph_ymax,"Y");
			}
			h_TDC_LE_TARGET[ican]->Draw();
			h_TDC_LE_TARGET_cut[ican]->Draw("same");
		}
	}

	if(flag==4){
		TCanvas *c_TARGET2;
		c_TARGET2 = new TCanvas(Name_Can_ADC_High2,Title_Can_ADC_High2,1200,500); 
		c_TARGET2->Divide(8,8);
		c_TARGET2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<32; ican++){
			c_TARGET2->cd(ican+1)->SetLogy();  
			h_ADC_High_TARGET[ican+96]->Draw();
			h_ADC_High_TARGET_cut[ican+96]->Draw("same");
			c_TARGET2->cd(ican+1+32);
			h_TDC_LE_TARGET[ican+96]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			h_TDC_LE_TARGET_cut[ican+96]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			if(TDC_Graph_ymax > 0) {  
				h_TDC_LE_TARGET[ican+96]->SetAxisRange(0, TDC_Graph_ymax,"Y");
				h_TDC_LE_TARGET_cut[ican+96]->SetAxisRange(0, TDC_Graph_ymax,"Y");
			}	
			h_TDC_LE_TARGET[ican+96]->Draw();
			h_TDC_LE_TARGET_cut[ican+96]->Draw("same");
		}


		TCanvas *c_TARGET;
		c_TARGET = new TCanvas(Name_Can_TDC_LE,Title_Can_TDC_LE,0,200,1200,500); 
		c_TARGET->Divide(8,8);
		c_TARGET->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<32; ican++){
			c_TARGET->cd(ican+1)->SetLogy();
			h_ADC_High_TARGET[ican+64]->Draw();
			h_ADC_High_TARGET_cut[ican+64]->Draw("same");
			c_TARGET->cd(ican+1+32);
			h_TDC_LE_TARGET[ican+64]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			h_TDC_LE_TARGET_cut[ican+64]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			if(TDC_Graph_ymax > 0) {  
				h_TDC_LE_TARGET[ican+64]->SetAxisRange(0, TDC_Graph_ymax,"Y");
				h_TDC_LE_TARGET_cut[ican+64]->SetAxisRange(0, TDC_Graph_ymax,"Y");
			}
			h_TDC_LE_TARGET[ican+64]->Draw();
			h_TDC_LE_TARGET_cut[ican+64]->Draw("same");
		}
	}

	if(flag==5){
		TCanvas *c_TARGET2;
		c_TARGET2 = new TCanvas(Name_Can_ADC_High2,Title_Can_ADC_High2,1200,500); 
		c_TARGET2->Divide(8,8);
		c_TARGET2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<32; ican++){
			c_TARGET2->cd(ican+1)->SetLogy();  
			h_ADC_High_TARGET[ican+160]->Draw();
			h_ADC_High_TARGET_cut[ican+160]->Draw("same");
			c_TARGET2->cd(ican+1+32);
			h_TDC_LE_TARGET[ican+160]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			h_TDC_LE_TARGET_cut[ican+160]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			if(TDC_Graph_ymax > 0) {  
				h_TDC_LE_TARGET[ican+160]->SetAxisRange(0, TDC_Graph_ymax,"Y");
				h_TDC_LE_TARGET_cut[ican+160]->SetAxisRange(0, TDC_Graph_ymax,"Y");
			}	
			h_TDC_LE_TARGET[ican+160]->Draw();
			h_TDC_LE_TARGET_cut[ican+160]->Draw("same");
		}


		TCanvas *c_TARGET;
		c_TARGET = new TCanvas(Name_Can_TDC_LE,Title_Can_TDC_LE,0,200,1200,500); 
		c_TARGET->Divide(8,8);
		c_TARGET->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<32; ican++){
			c_TARGET->cd(ican+1)->SetLogy();
			h_ADC_High_TARGET[ican+128]->Draw();
			h_ADC_High_TARGET_cut[ican+128]->Draw("same");
			c_TARGET->cd(ican+1+32);
			h_TDC_LE_TARGET[ican+128]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			h_TDC_LE_TARGET_cut[ican+128]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			if(TDC_Graph_ymax > 0) {  
				h_TDC_LE_TARGET[ican+128]->SetAxisRange(0, TDC_Graph_ymax,"Y");
				h_TDC_LE_TARGET_cut[ican+128]->SetAxisRange(0, TDC_Graph_ymax,"Y");
			}
			h_TDC_LE_TARGET[ican+128]->Draw();
			h_TDC_LE_TARGET_cut[ican+128]->Draw("same");
		}
	}

	if(flag==6){
		TCanvas *c_TARGET2;
		c_TARGET2 = new TCanvas(Name_Can_ADC_High2,Title_Can_ADC_High2,1200,500); 
		c_TARGET2->Divide(8,8);
		c_TARGET2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<32; ican++){
			c_TARGET2->cd(ican+1)->SetLogy();  
			h_ADC_High_TARGET[ican+224]->Draw();
			h_ADC_High_TARGET_cut[ican+224]->Draw("same");
			c_TARGET2->cd(ican+1+32);
			h_TDC_LE_TARGET[ican+224]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			h_TDC_LE_TARGET_cut[ican+224]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			if(TDC_Graph_ymax > 0) {  
				h_TDC_LE_TARGET[ican+224]->SetAxisRange(0, TDC_Graph_ymax,"Y");
				h_TDC_LE_TARGET_cut[ican+224]->SetAxisRange(0, TDC_Graph_ymax,"Y");
			}	
			h_TDC_LE_TARGET[ican+224]->Draw();
			h_TDC_LE_TARGET_cut[ican+224]->Draw("same");
		}


		TCanvas *c_TARGET;
		c_TARGET = new TCanvas(Name_Can_TDC_LE,Title_Can_TDC_LE,0,200,1200,500); 
		c_TARGET->Divide(8,8);
		c_TARGET->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<32; ican++){
			c_TARGET->cd(ican+1)->SetLogy();
			h_ADC_High_TARGET[ican+192]->Draw();
			h_ADC_High_TARGET_cut[ican+192]->Draw("same");
			c_TARGET->cd(ican+1+32);
			h_TDC_LE_TARGET[ican+192]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			h_TDC_LE_TARGET_cut[ican+192]->SetAxisRange(0, TDC_Graph_xmax,"X"); 
			if(TDC_Graph_ymax > 0) {  
				h_TDC_LE_TARGET[ican+192]->SetAxisRange(0, TDC_Graph_ymax,"Y");
				h_TDC_LE_TARGET_cut[ican+192]->SetAxisRange(0, TDC_Graph_ymax,"Y");
			}
			h_TDC_LE_TARGET[ican+192]->Draw();
			h_TDC_LE_TARGET_cut[ican+192]->Draw("same");
		}
	}

	if (flag==7) {
		TCanvas *c_TOF1;
		c_TOF1 = new TCanvas(Name_Can_TOF1,Title_Can_TOF1,1200,500); 
		c_TOF1->Divide(12,4);
		c_TOF1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<12; ican++){
			c_TOF1->cd(ican+1)->SetLogy();
			h_ADC_tof1U[ican]->Draw();
			h_ADC_tof1U_cut[ican]->Draw("same");
		}

		for(Int_t ican=0; ican<12; ican++){
			c_TOF1->cd(ican+1+12);
			h_TDC_tof1U[ican]->Draw();
			h_TDC_tof1U_cut[ican]->Draw("same");
		}

		for(Int_t ican=0; ican<12; ican++){
			c_TOF1->cd(ican+1+24)->SetLogy();
			h_ADC_tof1D[ican]->Draw();
			h_ADC_tof1D_cut[ican]->Draw("same");
		}

		for(Int_t ican=0; ican<12; ican++){
			c_TOF1->cd(ican+1+36);
			h_TDC_tof1D[ican]->Draw();
			h_TDC_tof1D_cut[ican]->Draw("same");
		}
	}

	if (flag == 8) {
		// Gap 7 - 12
		TCanvas *c_TOF2_2;
		c_TOF2_2 = new TCanvas(Name_Can_TOF2_2,Name_Can_TOF2_2,1200,700); 
		c_TOF2_2->Divide(6,8);
		c_TOF2_2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(Int_t ican=0; ican<6; ican++){
			c_TOF2_2->cd(ican+1)->SetLogy();
			h_ADC_tof2AO[ican+6]->Draw();
			h_ADC_tof2AO_cut[ican+6]->Draw("same");
		}

		for(Int_t ican=0; ican<6; ican++){
			c_TOF2_2->cd(ican+1+6);
			h_TDC_tof2AO[ican+6]->Draw();
			h_TDC_tof2AO_cut[ican+6]->Draw("same");
		}

		for(Int_t ican=0; ican<6; ican++){
			c_TOF2_2->cd(ican+1+12)->SetLogy();
			h_ADC_tof2AI[ican+6]->Draw();
			h_ADC_tof2AI_cut[ican+6]->Draw("same");
		}

		for(Int_t ican=0; ican<6; ican++){
			c_TOF2_2->cd(ican+1+18);
			h_TDC_tof2AI[ican+6]->Draw();
			h_TDC_tof2AI_cut[ican+6]->Draw("same");
		}

		for(Int_t ican=0; ican<6; ican++){
			c_TOF2_2->cd(ican+1+24)->SetLogy();
			h_ADC_tof2BO[ican+6]->Draw();
			h_ADC_tof2BO_cut[ican+6]->Draw("same");
		}

		for(Int_t ican=0; ican<6; ican++){
			c_TOF2_2->cd(ican+1+30);
			h_TDC_tof2BO[ican+6]->Draw();
			h_TDC_tof2BO_cut[ican+6]->Draw("same");
		}

		for(Int_t ican=0; ican<6; ican++){
			c_TOF2_2->cd(ican+1+36)->SetLogy();
			h_ADC_tof2BI[ican+6]->Draw();
			h_ADC_tof2BI_cut[ican+6]->Draw("same");
		}

		for(Int_t ican=0; ican<6; ican++){
			c_TOF2_2->cd(ican+1+42);
			h_TDC_tof2BI[ican+6]->Draw();
			h_TDC_tof2BI_cut[ican+6]->Draw("same");
		}


		// Gap 1 - 6
		TCanvas *c_TOF2_1;
		c_TOF2_1 = new TCanvas(Name_Can_TOF2_1,Name_Can_TOF2_1,1200,700); 
		c_TOF2_1->Divide(6,8);
		c_TOF2_1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");


		for(Int_t ican=0; ican<6; ican++){
			c_TOF2_1->cd(ican+1)->SetLogy();
			h_ADC_tof2AO[ican]->Draw();
			h_ADC_tof2AO_cut[ican]->Draw("same");
		}

		for(Int_t ican=0; ican<6; ican++){
			c_TOF2_1->cd(ican+1+6);
			h_TDC_tof2AO[ican]->Draw();
			h_TDC_tof2AO_cut[ican]->Draw("same");
		}

		for(Int_t ican=0; ican<6; ican++){
			c_TOF2_1->cd(ican+1+12)->SetLogy();
			h_ADC_tof2AI[ican]->Draw();
			h_ADC_tof2AI_cut[ican]->Draw("same");
		}

		for(Int_t ican=0; ican<6; ican++){
			c_TOF2_1->cd(ican+1+18);
			h_TDC_tof2AI[ican]->Draw();
			h_TDC_tof2AI_cut[ican]->Draw("same");
		}

		for(Int_t ican=0; ican<6; ican++){
			c_TOF2_1->cd(ican+1+24)->SetLogy();
			h_ADC_tof2BO[ican]->Draw();
			h_ADC_tof2BO_cut[ican]->Draw("same");
		}

		for(Int_t ican=0; ican<6; ican++){
			c_TOF2_1->cd(ican+1+30);
			h_TDC_tof2BO[ican]->Draw();
			h_TDC_tof2BO_cut[ican]->Draw("same");
		}

		for(Int_t ican=0; ican<6; ican++){
			c_TOF2_1->cd(ican+1+36)->SetLogy();
			h_ADC_tof2BI[ican]->Draw();
			h_ADC_tof2BI_cut[ican]->Draw("same");
		}

		for(Int_t ican=0; ican<6; ican++){
			c_TOF2_1->cd(ican+1+42);
			h_TDC_tof2BI[ican]->Draw();
			h_TDC_tof2BI_cut[ican]->Draw("same");
		}
	}


	n_counter = 0;
	n_plus_counter = 0;
	n_minus_counter = 0;
	miss_counter = 0;
	counter_test = 0;

	//cout << "TEST TEST TEST   " << gap_to_fit << "   " <<  gap_TOF2 << endl;

	return;
}

