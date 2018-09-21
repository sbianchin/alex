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
#include "TPaveStats.h"
#include "ANAPATH.h"
#include "Thresholds.h"
#include "CommonParameters.h"
#include "mappings.h"
#endif

void test_TDC_cut(Int_t run_number=3994, Int_t flag=0, Int_t fflag=0)
{
  	gStyle->Clear();
  	TH1::AddDirectory(kFALSE);
  	gStyle->SetOptStat(1111111);

  	if(flag==0) return;

  	int count=0;

 	// SFT
  	Int_t adc_high_sft[128];		Int_t ADC_High_SFT[128];
  	Int_t tdc_le_sft[128][16];		Int_t TDC_LE_SFT[128];

  	// TARGET
  	Int_t adc_high_target[256];		Int_t ADC_High_TARGET[256];
  	Int_t tdc_le_target[256][16];	Int_t TDC_LE_TARGET[256];

  	// TOF1
  	Int_t adc_tof1u[12];	Int_t tdc_tof1u[12];
  	Int_t adc_tof1d[12];	Int_t tdc_tof1d[12];

  	Int_t ADC_TOF1U[12];	Int_t TDC_TOF1U[12];
  	Int_t ADC_TOF1D[12];	Int_t TDC_TOF1D[12];

  	// TOF2
  	Int_t adc_tof2AO[12];	Int_t tdc_tof2AO[12];
  	Int_t adc_tof2AI[12];	Int_t tdc_tof2AI[12];
  	Int_t adc_tof2BO[12];	Int_t tdc_tof2BO[12];
  	Int_t adc_tof2BI[12];	Int_t tdc_tof2BI[12];

  	Int_t ADC_TOF2AO[12];	Int_t TDC_TOF2AO[12];
  	Int_t ADC_TOF2AI[12];	Int_t TDC_TOF2AI[12];
  	Int_t ADC_TOF2BO[12];	Int_t TDC_TOF2BO[12];
  	Int_t ADC_TOF2BI[12];	Int_t TDC_TOF2BI[12];

  	// AC
  	Int_t adc_acu[12];  	Int_t ADC_ACU[12];	
  	Int_t adc_acd[12];  	Int_t ADC_ACD[12];
	
	Int_t tdc_acu[12][8]; 		Int_t TDC_ACU[12];
	Int_t tdc_acd[12][8];		Int_t TDC_ACD[12];


 	//////////////////////////////////////////////////////////////////////////
  	/// OPEN FILES
	
	char file_mapping[200];
  	
  	char path_input[200];
	sprintf(path_input,"%s",path_merged);          

	char Name_finput[200];
	sprintf(Name_finput,"%s/Run%dMS.root",path_input, run_number);

 	 char path_cuts[200];
 	 sprintf(path_cuts,"%s","../Cuts");

	char Name_SFT_input[200];
	char Name_TARGET_input[200];
	char Name_TOF1_input[200];
	char Name_TOF2_input[200];
	char Name_AC_input[200];

  	Int_t SFT_TDC_min[128];			Int_t SFT_TDC_max[128];
  	Int_t TARGET_TDC_min[256];		Int_t TARGET_TDC_max[256];
  	Int_t TOF1_TDC_min[24];			Int_t TOF1_TDC_max[24];
  	Int_t TOF2_TDC_min[48];			Int_t TOF2_TDC_max[48];
  	Int_t AC_TDC_min[24];			Int_t AC_TDC_max[24];


  	if(flag==1 || flag==2){
  		sprintf(Name_SFT_input,"%s/SFT_TDC_Cuts.txt",path_cuts);

  		Int_t par_temp_SFT[3][128];
  		ifstream fdat_SFT(Name_SFT_input,ios::in);
  		for(int i=0; i<128; i++) fdat_SFT >> par_temp_SFT[0][i] >> par_temp_SFT[1][i] >> par_temp_SFT[2][i];
  		fdat_SFT.close();

	  	for(int j=0; j<128; j++){
  			SFT_TDC_min[j] = par_temp_SFT[1][j];
 			SFT_TDC_max[j] = par_temp_SFT[2][j];
  		}
  	}

  	if(flag==3 || flag==4 || flag==5 || flag==6){
  		sprintf(Name_TARGET_input,"%s/TARGET_TDC_Cuts.txt",path_cuts);
  		
  		Int_t par_temp_TARGET[3][256];
  		ifstream fdat_TARGET(Name_TARGET_input,ios::in);
  		for(int i=0; i<256; i++) fdat_TARGET >> par_temp_TARGET[0][i] >> par_temp_TARGET[1][i] >> par_temp_TARGET[2][i];
  		fdat_TARGET.close();
	
  		for(int j=0; j<256; j++){
  			TARGET_TDC_min[j] = par_temp_TARGET[1][j];
 			TARGET_TDC_max[j] = par_temp_TARGET[2][j];
  		}
  	}

  	if(flag==8){
  		sprintf(Name_TOF1_input,"%s/TOF1_TDC_Cuts.txt",path_cuts);
	
	  	Int_t par_temp_TOF1[3][24];
	  	ifstream fdat_TOF1(Name_TOF1_input,ios::in);
	  	for(int i=0; i<24; i++) fdat_TOF1 >> par_temp_TOF1[0][i] >> par_temp_TOF1[1][i] >> par_temp_TOF1[2][i];
	  	fdat_TOF1.close();

	  	for(int j=0; j<24; j++){
  			TOF1_TDC_min[j] = par_temp_TOF1[1][j];
 			TOF1_TDC_max[j] = par_temp_TOF1[2][j];
  		}
  	}

  	if(flag==9){
  		sprintf(Name_TOF2_input,"%s/TOF2_TDC_Cuts.txt",path_cuts);
  	
  		Int_t par_temp_TOF2[3][48];
  		ifstream fdat_TOF2(Name_TOF2_input,ios::in);
  		for(int i=0; i<48; i++) fdat_TOF2 >> par_temp_TOF2[0][i] >> par_temp_TOF2[1][i] >> par_temp_TOF2[2][i];
  		fdat_TOF2.close();

	  	for(int j=0; j<48; j++){
	  		TOF2_TDC_min[j] = par_temp_TOF2[1][j];
	 		TOF2_TDC_max[j] = par_temp_TOF2[2][j];
	 	}
  	}

  	if(flag==10){
  		sprintf(Name_AC_input,"%s/AC_TDC_Cuts.txt",path_cuts);

	  	Int_t par_temp_AC[3][24];
	  	ifstream fdat_AC(Name_AC_input,ios::in);
	  	for(int i=0; i<24; i++) fdat_AC >> par_temp_AC[0][i] >> par_temp_AC[1][i] >> par_temp_AC[2][i];
	  	fdat_AC.close();

	  	for(int j=0; j<24; j++){
	  		AC_TDC_min[j] = par_temp_AC[1][j];
	 		AC_TDC_max[j] = par_temp_AC[2][j];
	  	}
  	}

	
	cout << "   " << endl;
	cout << "/// ROOT FILE" << endl;
	cout << Name_finput << endl;

	cout << "   " << endl;
	cout << "/// CUT FILE(S)" << endl;
	if(flag==1 || flag==2) cout << Name_SFT_input << endl;
	if(flag==3 || flag==4 || flag==5 || flag==6) cout << Name_TARGET_input << endl;
	if(flag==8) cout << Name_TOF1_input << endl;
	if(flag==9) cout << Name_TOF2_input << endl;
	if(flag==10) cout << Name_AC_input << endl;
	cout << "   " << endl;

 	//////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////
	/// SET BRANCHES

	TChain *fChain= new TChain("Tree");		
	fChain->Add(Name_finput);		
	fChain->SetMakeClass(1);							


	// SFT
	if(flag==1 || flag==2){
		fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
		fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
	}

	// TARGET
	if(flag==3 || flag==4 || flag==5 || flag==6){	
		fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);		
		fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);		
	}

	// TOF1
	if(flag==8){
		fChain->SetBranchAddress("ADC_TOF1U",adc_tof1u);
		fChain->SetBranchAddress("ADC_TOF1D",adc_tof1d);
		fChain->SetBranchAddress("TDC_TOF1U",tdc_tof1u);
		fChain->SetBranchAddress("TDC_TOF1D",tdc_tof1d);
	}

	// TOF2
	if(flag==9){
		fChain->SetBranchAddress("ADC_TOF2AO",adc_tof2AO);
		fChain->SetBranchAddress("ADC_TOF2AI",adc_tof2AI);
		fChain->SetBranchAddress("ADC_TOF2BO",adc_tof2BO);
		fChain->SetBranchAddress("ADC_TOF2BI",adc_tof2BI);	
		fChain->SetBranchAddress("TDC_TOF2AO",tdc_tof2AO);
		fChain->SetBranchAddress("TDC_TOF2AI",tdc_tof2AI);
		fChain->SetBranchAddress("TDC_TOF2BO",tdc_tof2BO);
		fChain->SetBranchAddress("TDC_TOF2BI",tdc_tof2BI);	
	}

	// AC
	if(flag==10){
		fChain->SetBranchAddress("ADC_ACU",adc_acu);
		fChain->SetBranchAddress("ADC_ACD",adc_acd);
		fChain->SetBranchAddress("TDC_ACU",tdc_acu);
		fChain->SetBranchAddress("TDC_ACD",tdc_acd);
	}

	//////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////
	/// HISTOGRAM DECLARATION

	// SFT
	TH1D *h_ADC_High_SFT[128];	char Title_ADC_High_SFT[128][100];	char Name_ADC_High_SFT[128][100];
	TH1D *h_TDC_LE_SFT[128];    char Title_TDC_LE_SFT[128][100];	char Name_TDC_LE_SFT[128][100];     

	TH1D *h_ADC_High_SFT_Raw[128];	char Title_ADC_High_SFT_Raw[128][100];	char Name_ADC_High_SFT_Raw[128][100];
	TH1D *h_TDC_LE_SFT_Raw[128];    char Title_TDC_LE_SFT_Raw[128][100];	char Name_TDC_LE_SFT_Raw[128][100];     

	// TARGET
	TH1D *h_ADC_High_TARGET[256];	char Title_ADC_High_TARGET[256][100];	char Name_ADC_High_TARGET[256][100];
	TH1D *h_TDC_LE_TARGET[256];     char Title_TDC_LE_TARGET[256][100];	char Name_TDC_LE_TARGET[256][100];     

	TH1D *h_ADC_High_TARGET_Raw[256];	char Title_ADC_High_TARGET_Raw[256][100];	char Name_ADC_High_TARGET_Raw[256][100];
	TH1D *h_TDC_LE_TARGET_Raw[256];     char Title_TDC_LE_TARGET_Raw[256][100];		char Name_TDC_LE_TARGET_Raw[256][100];     

	// TOF1
	TH1D *h_ADC_tof1U[12];   char Title_ADC_tof1U[12][100];	char Name_ADC_tof1U[12][100];
	TH1D *h_ADC_tof1D[12];   char Title_ADC_tof1D[12][100];	char Name_ADC_tof1D[12][100];
	TH1D *h_TDC_tof1U[12];   char Title_TDC_tof1U[12][100];	char Name_TDC_tof1U[12][100];
	TH1D *h_TDC_tof1D[12];   char Title_TDC_tof1D[12][100];	char Name_TDC_tof1D[12][100];

	TH1D *h_ADC_tof1U_Raw[12];   char Title_ADC_tof1U_Raw[12][100];	char Name_ADC_tof1U_Raw[12][100];
	TH1D *h_ADC_tof1D_Raw[12];   char Title_ADC_tof1D_Raw[12][100];	char Name_ADC_tof1D_Raw[12][100];
	TH1D *h_TDC_tof1U_Raw[12];   char Title_TDC_tof1U_Raw[12][100];	char Name_TDC_tof1U_Raw[12][100];
	TH1D *h_TDC_tof1D_Raw[12];   char Title_TDC_tof1D_Raw[12][100];	char Name_TDC_tof1D_Raw[12][100];

	// TOF2
	TH1D *h_ADC_tof2AO[12];   char Title_ADC_tof2AO[12][100];	char Name_ADC_tof2AO[12][100];
	TH1D *h_ADC_tof2AI[12];   char Title_ADC_tof2AI[12][100];	char Name_ADC_tof2AI[12][100];
	TH1D *h_ADC_tof2BO[12];   char Title_ADC_tof2BO[12][100];	char Name_ADC_tof2BO[12][100];
	TH1D *h_ADC_tof2BI[12];   char Title_ADC_tof2BI[12][100];	char Name_ADC_tof2BI[12][100];
	TH1D *h_TDC_tof2AO[12];   char Title_TDC_tof2AO[12][100];	char Name_TDC_tof2AO[12][100];
	TH1D *h_TDC_tof2AI[12];   char Title_TDC_tof2AI[12][100];	char Name_TDC_tof2AI[12][100];
	TH1D *h_TDC_tof2BO[12];   char Title_TDC_tof2BO[12][100];	char Name_TDC_tof2BO[12][100];
	TH1D *h_TDC_tof2BI[12];   char Title_TDC_tof2BI[12][100];	char Name_TDC_tof2BI[12][100];

	TH1D *h_ADC_tof2AO_Raw[12];   char Title_ADC_tof2AO_Raw[12][100];	char Name_ADC_tof2AO_Raw[12][100];
	TH1D *h_ADC_tof2AI_Raw[12];   char Title_ADC_tof2AI_Raw[12][100];	char Name_ADC_tof2AI_Raw[12][100];
	TH1D *h_ADC_tof2BO_Raw[12];   char Title_ADC_tof2BO_Raw[12][100];	char Name_ADC_tof2BO_Raw[12][100];
	TH1D *h_ADC_tof2BI_Raw[12];   char Title_ADC_tof2BI_Raw[12][100];	char Name_ADC_tof2BI_Raw[12][100];
	TH1D *h_TDC_tof2AO_Raw[12];   char Title_TDC_tof2AO_Raw[12][100];	char Name_TDC_tof2AO_Raw[12][100];
	TH1D *h_TDC_tof2AI_Raw[12];   char Title_TDC_tof2AI_Raw[12][100];	char Name_TDC_tof2AI_Raw[12][100];
	TH1D *h_TDC_tof2BO_Raw[12];   char Title_TDC_tof2BO_Raw[12][100];	char Name_TDC_tof2BO_Raw[12][100];
	TH1D *h_TDC_tof2BI_Raw[12];   char Title_TDC_tof2BI_Raw[12][100];	char Name_TDC_tof2BI_Raw[12][100];

	// AC
	TH1D *h_ADC_acu[12];	char Title_ADC_acu[12][100];	char Name_ADC_acu[12][100];
	TH1D *h_ADC_acd[12];	char Title_ADC_acd[12][100];	char Name_ADC_acd[12][100];
	TH1D *h_TDC_acu[12];	char Title_TDC_acu[12][100];	char Name_TDC_acu[12][100];
	TH1D *h_TDC_acd[12];	char Title_TDC_acd[12][100];	char Name_TDC_acd[12][100];

	TH1D *h_ADC_acu_Raw[12];	char Title_ADC_acu_Raw[12][100];	char Name_ADC_acu_Raw[12][100];
	TH1D *h_ADC_acd_Raw[12];	char Title_ADC_acd_Raw[12][100];	char Name_ADC_acd_Raw[12][100];
	TH1D *h_TDC_acu_Raw[12];	char Title_TDC_acu_Raw[12][100];	char Name_TDC_acu_Raw[12][100];
	TH1D *h_TDC_acd_Raw[12];	char Title_TDC_acd_Raw[12][100];	char Name_TDC_acd_Raw[12][100];

	//////////////////////////////////////////////////////////////////////////

	
	//////////////////////////////////////////////////////////////////////////
	/// SET HISTOGRAM TITLES & NAMES

	// SFT
	if(flag==1 || flag==2){
		for(int i=0; i<128; i++){
			sprintf(Title_ADC_High_SFT[i],"ADC High Gain (Ch. %d)  --  SFT",i); 
			sprintf(Title_TDC_LE_SFT[i],"TDC (LE) (Ch. %d)  --  SFT",i); 

			sprintf(Name_ADC_High_SFT[i],"ADC_High (Ch. %d) - SFT",i);	
			sprintf(Name_TDC_LE_SFT[i],"TDC_LE (Ch. %d) - SFT",i);

			sprintf(Title_ADC_High_SFT_Raw[i],"Raw ADC High Gain (Ch. %d)  --  SFT",i); 
			sprintf(Title_TDC_LE_SFT_Raw[i],"Raw TDC (LE) (Ch. %d)  --  SFT",i); 

			sprintf(Name_ADC_High_SFT_Raw[i],"ADC_High (Ch. %d) - SFT",i);	
			sprintf(Name_TDC_LE_SFT_Raw[i],"TDC_LE (Ch. %d) - SFT",i);
		}
	}	

	// TARGET
	if(flag==3 || flag==4 || flag==5 || flag==6){
		for(int i=0; i<256; i++){
			sprintf(Title_ADC_High_TARGET[i],"ADC High Gain (Ch. %d)  --  TARGET",i); 
			sprintf(Title_TDC_LE_TARGET[i],"TDC (LE) (Ch. %d)  --  TARGET",i); 

			sprintf(Name_ADC_High_TARGET[i],"ADC_High (Ch. %d) - TARGET",i);	
			sprintf(Name_TDC_LE_TARGET[i],"TDC_LE (Ch. %d) - TARGET",i);

			sprintf(Title_ADC_High_TARGET_Raw[i],"Raw ADC High Gain (Ch. %d)  --  TARGET",i); 
			sprintf(Title_TDC_LE_TARGET_Raw[i],"Raw TDC (LE) (Ch. %d)  --  TARGET",i); 

			sprintf(Name_ADC_High_TARGET_Raw[i],"ADC_High (Ch. %d) - TARGET",i);	
			sprintf(Name_TDC_LE_TARGET_Raw[i],"TDC_LE (Ch. %d) - TARGET",i);
		}
	}	

	// TOF1
	if(flag==8){
		for(int i=0; i<12; i++){
			sprintf(Title_ADC_tof1U[i],"ADC (V792 Ch. %d)  --  TOF1 (%d UP)",TOF1U_ADC[i],i+1);
			sprintf(Title_ADC_tof1D[i],"ADC (V792 Ch. %d)  --  TOF1 (%d DOWN)",TOF1D_ADC[i],i+1);
			sprintf(Title_TDC_tof1U[i],"TDC (HR_TDC Ch. %d)  --  TOF1 (%d UP)",TOF1U_TDC[i],i+1);
			sprintf(Title_TDC_tof1D[i],"TDC (HR _TDC Ch. %d)  --  TOF1 (%d DOWN)",TOF1D_TDC[i],i+1);

			sprintf(Name_ADC_tof1U[i],"ADC_TOF1 (%d UP)",i+1);
			sprintf(Name_ADC_tof1D[i],"ADC_TOF1 (%d DOWN)",i+1);
			sprintf(Name_TDC_tof1U[i],"TDC_TOF1 (%d UP)",i+1);
			sprintf(Name_TDC_tof1D[i],"TDC_TOF1 (%d DOWN)",i+1);

			sprintf(Title_ADC_tof1U_Raw[i],"Raw ADC (V792 Ch. %d)  --  TOF1 (%d UP)",TOF1U_ADC[i],i+1);
			sprintf(Title_ADC_tof1D_Raw[i],"Raw ADC (V792 Ch. %d)  --  TOF1 (%d DOWN)",TOF1D_ADC[i],i+1);
			sprintf(Title_TDC_tof1U_Raw[i],"Raw TDC (HR_TDC Ch. %d)  --  TOF1 (%d UP)",TOF1U_TDC[i],i+1);
			sprintf(Title_TDC_tof1D_Raw[i],"Raw TDC (HR _TDC Ch. %d)  --  TOF1 (%d DOWN)",TOF1D_TDC[i],i+1);

			sprintf(Name_ADC_tof1U_Raw[i],"ADC_TOF1 (%d UP)",i+1);
			sprintf(Name_ADC_tof1D_Raw[i],"ADC_TOF1 (%d DOWN)",i+1);
			sprintf(Name_TDC_tof1U_Raw[i],"TDC_TOF1 (%d UP)",i+1);
			sprintf(Name_TDC_tof1D_Raw[i],"TDC_TOF1 (%d DOWN)",i+1);
		}
	}

	// TOF2
	if(flag==9){
		for(int i=0; i<12; i++){
			sprintf(Title_ADC_tof2AO[i],"ADC (V792 Ch. %d)  --  TOF2 (%d A OUT)",TOF2AO_ADC[i],i+1);
			sprintf(Title_ADC_tof2AI[i],"ADC (V792 Ch. %d)  --  TOF2 (%d A IN)",TOF2AI_ADC[i],i+1);
			sprintf(Title_ADC_tof2BO[i],"ADC (V792 Ch. %d)  --  TOF2 (%d B OUT)",TOF2BO_ADC[i],i+1);
			sprintf(Title_ADC_tof2BI[i],"ADC (V792 Ch. %d)  --  TOF2 (%d B IN)",TOF2BI_ADC[i],i+1);
			sprintf(Title_TDC_tof2AO[i],"TDC (HR_TDC Ch. %d)  --  TOF2 (%d A OUT)",TOF2AO_TDC[i],i+1);
			sprintf(Title_TDC_tof2AI[i],"TDC (HR_TDC Ch. %d)  --  TOF2 (%d A IN)",TOF2AI_TDC[i],i+1);
			sprintf(Title_TDC_tof2BO[i],"TDC (HR_TDC Ch. %d)  --  TOF2 (%d B OUT)",TOF2BO_TDC[i],i+1);
			sprintf(Title_TDC_tof2BI[i],"TDC (HR_TDC Ch. %d)  --  TOF2 (%d B IN)",TOF2BI_TDC[i],i+1);

			sprintf(Name_ADC_tof2AO[i],"ADC_TOF2 (%d AO)",i+1);
			sprintf(Name_ADC_tof2AI[i],"ADC_TOF2 (%d AI)",i+1);
			sprintf(Name_ADC_tof2BO[i],"ADC_TOF2 (%d BO)",i+1);
			sprintf(Name_ADC_tof2BI[i],"ADC_TOF2 (%d BI)",i+1);
			sprintf(Name_TDC_tof2AO[i],"TDC_TOF2 (%d AO)",i+1);
			sprintf(Name_TDC_tof2AI[i],"TDC_TOF2 (%d AI)",i+1);
			sprintf(Name_TDC_tof2BO[i],"TDC_TOF2 (%d BO)",i+1);
			sprintf(Name_TDC_tof2BI[i],"TDC_TOF2 (%d BI)",i+1);

			sprintf(Title_ADC_tof2AO_Raw[i],"Raw ADC (V792 Ch. %d)  --  TOF2 (%d A OUT)",TOF2AO_ADC[i],i+1);
			sprintf(Title_ADC_tof2AI_Raw[i],"Raw ADC (V792 Ch. %d)  --  TOF2 (%d A IN)",TOF2AI_ADC[i],i+1);
			sprintf(Title_ADC_tof2BO_Raw[i],"Raw ADC (V792 Ch. %d)  --  TOF2 (%d B OUT)",TOF2BO_ADC[i],i+1);
			sprintf(Title_ADC_tof2BI_Raw[i],"Raw ADC (V792 Ch. %d)  --  TOF2 (%d B IN)",TOF2BI_ADC[i],i+1);
			sprintf(Title_TDC_tof2AO_Raw[i],"Raw TDC (HR_TDC Ch. %d)  --  TOF2 (%d A OUT)",TOF2AO_TDC[i],i+1);
			sprintf(Title_TDC_tof2AI_Raw[i],"Raw TDC (HR_TDC Ch. %d)  --  TOF2 (%d A IN)",TOF2AI_TDC[i],i+1);
			sprintf(Title_TDC_tof2BO_Raw[i],"Raw TDC (HR_TDC Ch. %d)  --  TOF2 (%d B OUT)",TOF2BO_TDC[i],i+1);
			sprintf(Title_TDC_tof2BI_Raw[i],"Raw TDC (HR_TDC Ch. %d)  --  TOF2 (%d B IN)",TOF2BI_TDC[i],i+1);

			sprintf(Name_ADC_tof2AO_Raw[i],"ADC_TOF2 (%d AO)",i+1);
			sprintf(Name_ADC_tof2AI_Raw[i],"ADC_TOF2 (%d AI)",i+1);
			sprintf(Name_ADC_tof2BO_Raw[i],"ADC_TOF2 (%d BO)",i+1);
			sprintf(Name_ADC_tof2BI_Raw[i],"ADC_TOF2 (%d BI)",i+1);
			sprintf(Name_TDC_tof2AO_Raw[i],"TDC_TOF2 (%d AO)",i+1);
			sprintf(Name_TDC_tof2AI_Raw[i],"TDC_TOF2 (%d AI)",i+1);
			sprintf(Name_TDC_tof2BO_Raw[i],"TDC_TOF2 (%d BO)",i+1);
			sprintf(Name_TDC_tof2BI_Raw[i],"TDC_TOF2 (%d BI)",i+1);
		}
	}

	// AC
	if(flag==10){
		for(int i=0; i<12; i++){
			sprintf(Title_ADC_acu[i],"ADC (TKO_ADC Ch. %d)  --  AC (%d UP)",ACU_ADC[i],i+1);
			sprintf(Title_ADC_acd[i],"ADC (TKO_ADC Ch. %d)  --  AC (%d DOWN)",ACD_ADC[i],i+1);
			sprintf(Title_TDC_acu[i],"TDC (VT48_TDC Ch. %d)  --  AC (%d UP)",ACU_TDC[i],i+1);
			sprintf(Title_TDC_acd[i],"TDC (VT48_TDC Ch. %d)  --  AC (%d DOWN)",ACD_TDC[i],i+1);
	
			sprintf(Name_ADC_acu[i],"ADC_AC (%d UP)",i+1);
			sprintf(Name_ADC_acd[i],"ADC_AC (%d DOWN)",i+1);
			sprintf(Name_TDC_acu[i],"TDC_AC (%d UP)",i+1);
			sprintf(Name_TDC_acd[i],"TDC_AC (%d DOWN)",i+1);

			sprintf(Title_ADC_acu_Raw[i],"Raw ADC (TKO_ADC Ch. %d)  --  AC (%d UP)",ACU_ADC[i],i+1);
			sprintf(Title_ADC_acd_Raw[i],"Raw ADC (TKO_ADC Ch. %d)  --  AC (%d DOWN)",ACD_ADC[i],i+1);
			sprintf(Title_TDC_acu_Raw[i],"Raw TDC (VT48_TDC Ch. %d)  --  AC (%d UP)",ACU_TDC[i],i+1);
			sprintf(Title_TDC_acd_Raw[i],"Raw TDC (VT48_TDC Ch. %d)  --  AC (%d DOWN)",ACD_TDC[i],i+1);
	
			sprintf(Name_ADC_acu_Raw[i],"ADC_AC (%d UP)",i+1);
			sprintf(Name_ADC_acd_Raw[i],"ADC_AC (%d DOWN)",i+1);
			sprintf(Name_TDC_acu_Raw[i],"TDC_AC (%d UP)",i+1);
			sprintf(Name_TDC_acd_Raw[i],"TDC_AC (%d DOWN)",i+1);
		}
	}
		
	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	/// CREATING HISTOGRAMS

	// SFT
	if(flag==1 || flag==2){
		for (int i=0; i<128; i++){
			h_ADC_High_SFT[i] = new TH1D(Name_ADC_High_SFT[i],Title_ADC_High_SFT[i],512,0,4096);
			h_TDC_LE_SFT[i] = new TH1D(Name_TDC_LE_SFT[i],Title_TDC_LE_SFT[i],512,0,4096);

			h_ADC_High_SFT_Raw[i] = new TH1D(Name_ADC_High_SFT_Raw[i],Title_ADC_High_SFT_Raw[i],512,0,4096);
			h_TDC_LE_SFT_Raw[i] = new TH1D(Name_TDC_LE_SFT_Raw[i],Title_TDC_LE_SFT_Raw[i],512,0,4096);

			h_ADC_High_SFT[i]->SetLineWidth(2);
			h_ADC_High_SFT[i]->SetLineColor(4);
		}
	}

	// TARGET
	if(flag==3 || flag==4 || flag==5 || flag==6){
		for (int i=0; i<256; i++){
			h_ADC_High_TARGET[i] = new TH1D(Name_ADC_High_TARGET[i],Title_ADC_High_TARGET[i],512,0,4096);
			h_TDC_LE_TARGET[i] = new TH1D(Name_TDC_LE_TARGET[i],Title_TDC_LE_TARGET[i],512,0,4096);

			h_ADC_High_TARGET_Raw[i] = new TH1D(Name_ADC_High_TARGET_Raw[i],Title_ADC_High_TARGET_Raw[i],512,0,4096);
			h_TDC_LE_TARGET_Raw[i] = new TH1D(Name_TDC_LE_TARGET_Raw[i],Title_TDC_LE_TARGET_Raw[i],512,0,4096);

			h_ADC_High_TARGET[i]->SetLineWidth(2);
			h_ADC_High_TARGET[i]->SetLineColor(4);
		}
	}

	// TOF1
	if(flag==8){
		for(int i=0; i<12; i++){
			h_ADC_tof1U[i] = new TH1D(Name_ADC_tof1U[i], Title_ADC_tof1U[i], 512, 0, 4096);
			h_ADC_tof1D[i] = new TH1D(Name_ADC_tof1D[i], Title_ADC_tof1D[i], 512, 0, 4096);
			h_TDC_tof1U[i] = new TH1D(Name_TDC_tof1U[i],Title_TDC_tof1U[i],512,0,4096);
			h_TDC_tof1D[i] = new TH1D(Name_TDC_tof1D[i],Title_TDC_tof1D[i],512,0,4096);

			h_ADC_tof1U_Raw[i] = new TH1D(Name_ADC_tof1U_Raw[i], Title_ADC_tof1U_Raw[i], 512, 0, 4096);
			h_ADC_tof1D_Raw[i] = new TH1D(Name_ADC_tof1D_Raw[i], Title_ADC_tof1D_Raw[i], 512, 0, 4096);
			h_TDC_tof1U_Raw[i] = new TH1D(Name_TDC_tof1U_Raw[i],Title_TDC_tof1U_Raw[i],512,0,4096);
			h_TDC_tof1D_Raw[i] = new TH1D(Name_TDC_tof1D_Raw[i],Title_TDC_tof1D_Raw[i],512,0,4096);

			h_ADC_tof1U[i]->SetLineWidth(2);
			h_ADC_tof1U[i]->SetLineColor(4);

			h_ADC_tof1D[i]->SetLineWidth(2);
			h_ADC_tof1D[i]->SetLineColor(4);
		}
	}

	// TOF2
	if(flag==9){
		for(int i=0; i<12; i++){
			h_ADC_tof2AO[i] = new TH1D(Name_ADC_tof2AO[i],Title_ADC_tof2AO[i],512,0,4096);
			h_ADC_tof2AI[i] = new TH1D(Name_ADC_tof2AI[i],Title_ADC_tof2AI[i],512,0,4096);
			h_ADC_tof2BO[i] = new TH1D(Name_ADC_tof2BO[i],Title_ADC_tof2BO[i],512,0,4096);
			h_ADC_tof2BI[i] = new TH1D(Name_ADC_tof2BI[i],Title_ADC_tof2BI[i],512,0,4096);
			h_TDC_tof2AO[i] = new TH1D(Name_TDC_tof2AO[i],Title_TDC_tof2AO[i],512,0,4096);
			h_TDC_tof2AI[i] = new TH1D(Name_TDC_tof2AI[i],Title_TDC_tof2AI[i],512,0,4096);
			h_TDC_tof2BO[i] = new TH1D(Name_TDC_tof2BO[i],Title_TDC_tof2BO[i],512,0,4096);
			h_TDC_tof2BI[i] = new TH1D(Name_TDC_tof2BI[i],Title_TDC_tof2BI[i],512,0,4096);

			h_ADC_tof2AO_Raw[i] = new TH1D(Name_ADC_tof2AO_Raw[i],Title_ADC_tof2AO_Raw[i],512,0,4096);
			h_ADC_tof2AI_Raw[i] = new TH1D(Name_ADC_tof2AI_Raw[i],Title_ADC_tof2AI_Raw[i],512,0,4096);
			h_ADC_tof2BO_Raw[i] = new TH1D(Name_ADC_tof2BO_Raw[i],Title_ADC_tof2BO_Raw[i],512,0,4096);
			h_ADC_tof2BI_Raw[i] = new TH1D(Name_ADC_tof2BI_Raw[i],Title_ADC_tof2BI_Raw[i],512,0,4096);
			h_TDC_tof2AO_Raw[i] = new TH1D(Name_TDC_tof2AO_Raw[i],Title_TDC_tof2AO_Raw[i],512,0,4096);
			h_TDC_tof2AI_Raw[i] = new TH1D(Name_TDC_tof2AI_Raw[i],Title_TDC_tof2AI_Raw[i],512,0,4096);
			h_TDC_tof2BO_Raw[i] = new TH1D(Name_TDC_tof2BO_Raw[i],Title_TDC_tof2BO_Raw[i],512,0,4096);
			h_TDC_tof2BI_Raw[i] = new TH1D(Name_TDC_tof2BI_Raw[i],Title_TDC_tof2BI_Raw[i],512,0,4096);

			h_ADC_tof2AO[i]->SetLineWidth(2);
			h_ADC_tof2AO[i]->SetLineColor(4);

			h_ADC_tof2AI[i]->SetLineWidth(2);
			h_ADC_tof2AI[i]->SetLineColor(4);

			h_ADC_tof2BO[i]->SetLineWidth(2);
			h_ADC_tof2BO[i]->SetLineColor(4);

			h_ADC_tof2BI[i]->SetLineWidth(2);
			h_ADC_tof2BI[i]->SetLineColor(4);
		}
	}

	// AC
	if(flag==10){
		for(int i=0; i<12; i++){
			h_ADC_acu[i] = new TH1D(Name_ADC_acu[i], Title_ADC_acu[i], 512, 0, 4096);
			h_ADC_acd[i] = new TH1D(Name_ADC_acd[i], Title_ADC_acd[i], 512, 0, 4096);
			h_TDC_acu[i] = new TH1D(Name_TDC_acu[i], Title_TDC_acu[i], 512, 0, 4096);
			h_TDC_acd[i] = new TH1D(Name_TDC_acd[i], Title_TDC_acd[i], 512, 0, 4096);

			h_ADC_acu_Raw[i] = new TH1D(Name_ADC_acu_Raw[i], Title_ADC_acu_Raw[i], 512, 0, 4096);
			h_ADC_acd_Raw[i] = new TH1D(Name_ADC_acd_Raw[i], Title_ADC_acd_Raw[i], 512, 0, 4096);
			h_TDC_acu_Raw[i] = new TH1D(Name_TDC_acu_Raw[i], Title_TDC_acu_Raw[i], 512, 0, 4096);
			h_TDC_acd_Raw[i] = new TH1D(Name_TDC_acd_Raw[i], Title_TDC_acd_Raw[i], 512, 0, 4096);

			h_ADC_acu[i]->SetLineWidth(2);
			h_ADC_acu[i]->SetLineColor(4);

			h_ADC_acd[i]->SetLineWidth(2);
			h_ADC_acd[i]->SetLineColor(4);
		}
	}

	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	/// CREATE CANVAS
	/// CANVAS NAMES & TITLES 
	
	// SFT
	char Name_Can_ADC_High_SFT[100];	char Title_Can_ADC_High_SFT[100];
	char Name_Can_ADC_High_SFT2[100];	char Title_Can_ADC_High_SFT2[100];
	char Name_Can_TDC_LE_SFT[100];		char Title_Can_TDC_LE_SFT[100];
	char Name_Can_TDC_LE_SFT2[100];		char Title_Can_TDC_LE_SFT2[100];

	// TARGET 
	char Name_Can_ADC_High_TARGET[100];		char Title_Can_ADC_High_TARGET[100];
	char Name_Can_ADC_High_TARGET2[100];	char Title_Can_ADC_High_TARGET2[100];
	char Name_Can_TDC_LE_TARGET[100];		char Title_Can_TDC_LE_TARGET[100];
	char Name_Can_TDC_LE_TARGET2[100];		char Title_Can_TDC_LE_TARGET2[100];

	// TOF1
	char Name_Can_TOF1[100];	char Title_Can_TOF1[100];

	// TOF2
	char Name_Can_TOF2[100];	char Title_Can_TOF2[100];
	char Name_Can_TOF2_2[100];	char Title_Can_TOF2_2[100];

	// AC
	char Name_Can_AC[100];	char Title_Can_AC[100];
	

	// SFT (ch. 0 - 63)
	if(flag==1){
		sprintf(Name_Can_ADC_High_SFT,"EASIROC ADC(HG) -- Run %d  Ch. 0 - 63  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_High_SFT,"EASIROC ADC(HG) -- Run %d  Ch. 0 - 63  (SFT)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_TDC_LE_SFT,"EASIROC TDC(LE) -- Run %d  Ch. 0 - 63  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_LE_SFT,"EASIROC TDC(LE) -- Run %d  Ch. 0 - 63  (SFT)  | EASIROC  Board %d",run_number,flag);
	}

	// SFT (ch. 64 - 127)
	if(flag==2){
		sprintf(Name_Can_ADC_High_SFT,"EASIROC ADC(HG) -- Run %d  Ch. 64 - 127 (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_High_SFT,"EASIROC ADC(HG) -- Run %d  Ch. 64 - 127 (SFT)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_TDC_LE_SFT,"EASIROC TDC(LE) -- Run %d  Ch. 64 - 127  (SFT)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_LE_SFT,"EASIROC TDC(LE) -- Run %d  Ch. 64 - 127  (SFT)  |  EASIROC Board %d",run_number,flag);
	}

	// TARGET (ch. 0 - 63)
	if(flag==3){
		sprintf(Name_Can_ADC_High_TARGET,"EASIROC ADC(HG) -- Run %d  Ch. 0 - 63  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_High_TARGET,"EASIROC ADC(HG) -- Run %d  Ch. 0 - 63  (TARGET)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_TDC_LE_TARGET,"EASIROC TDC(LE) -- Run %d  Ch. 0 - 63  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_LE_TARGET,"EASIROC TDC(LE) -- Run %d  Ch. 0 - 63  (TARGET)  |  EASIROC Board %d",run_number,flag);
	}

	// TARGET (ch. 64 - 127)
	if(flag==4){
		sprintf(Name_Can_ADC_High_TARGET,"EASIROC ADC(HG) -- Run %d  Ch. 64 - 127  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_High_TARGET,"EASIROC ADC(HG) -- Run %d  Ch. 64 - 127  (TARGET)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_TDC_LE_TARGET,"EASIROC TDC(LE) -- Run %d  Ch. 64 - 127  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_LE_TARGET,"EASIROC TDC(LE) -- Run %d  Ch. 64 - 127  (TARGET)  |  EASIROC Board %d",run_number,flag);
	}

	// TARGET (ch. 128 - 191)
	if(flag==5){
		sprintf(Name_Can_ADC_High_TARGET,"EASIROC ADC(HG) -- Run %d  Ch. 128 - 191  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_High_TARGET,"EASIROC ADC(HG) -- Run %d  Ch. 128 - 191  (TARGET)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_TDC_LE_TARGET,"EASIROC TDC(LE) -- Run %d  Ch. 128 - 191  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_LE_TARGET,"EASIROC TDC(LE) -- Run %d  Ch. 128 - 191  (TARGET)  |  EASIROC Board %d",run_number,flag);
	}

	// TARGET (ch. 192 - 255)
	if(flag==6){
		sprintf(Name_Can_ADC_High_TARGET,"EASIROC ADC(HG) -- Run %d  Ch. 192 - 255  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_ADC_High_TARGET,"EASIROC ADC(HG) -- Run %d  Ch. 192 - 255  (TARGET)  |  EASIROC Board %d",run_number,flag);

		sprintf(Name_Can_TDC_LE_TARGET,"EASIROC TDC(LE) -- Run %d  Ch. 192 - 255  (TARGET)  |  EASIROC Board %d",run_number,flag);
		sprintf(Title_Can_TDC_LE_TARGET,"EASIROC TDC(LE) -- Run %d  Ch. 192 - 255  (TARGET)  |  EASIROC Board %d",run_number,flag);
	}

	// TOF1 
	if(flag==8){
		sprintf(Name_Can_TOF1,"TOF1 ADCs & TDCs -- Run %d  (Gap 1 - 12)",run_number);
		sprintf(Title_Can_TOF1,"TOF1 ADCs & TDCs -- Run %d  (Gap 1 - 12)",run_number);
	}

	// TOF2
	if(flag==9){
		sprintf(Name_Can_TOF2,"TOF2 ADCs & TDCs -- Run %d  (Gap 1 - 6)",run_number);
		sprintf(Title_Can_TOF2,"TOF2 ADCs & TDCs -- Run %d  (Gap 1 -6)",run_number);
		sprintf(Name_Can_TOF2_2,"TOF2 ADCs & TDCs -- Run %d  (Gap 7 - 12)",run_number);
		sprintf(Title_Can_TOF2_2,"TOF2 ADCs & TDCs -- Run %d  (Gap 7 - 12)",run_number);
	}

	// AC
	if(flag==10){
		sprintf(Name_Can_AC,"Aerogel Counter -- Run %d  (Gap 1 - 12)",run_number);
		sprintf(Title_Can_AC,"Aerogel Counter -- Run %d  (Gap 1 - 12)",run_number);
	}

	//////////////////////////////////////////////////////////////////////////
 
	/// READ ALL ENTRIES
	Int_t nentries = (Int_t)fChain->GetEntries();
	
	cout <<  "Total Number of Entries :     " <<  nentries << endl;

	cout << "  " << endl;
	cout << "***************  Processing....  ***************" << endl;
	cout << "  " << endl;

	if(fflag==0) count = nentries;
	if(fflag!=0) count = fflag;

	for(Int_t i=0; i<count; i++){
	//for(Int_t i=0; i<10000; i++){
		fChain->GetEntry(i);

		if(nentries<=30000){
			if(i%1000==1) cout<<"**** "<<i<<" events done"<<endl;
		}
		if(nentries>30000){
			if(i%10000==1) cout<<"**** "<<i<<" events done"<<endl;
		}

		//////////////////////////////////////////////////////////////
		/// RAW HISTOGRAMS

		// SFT
		if(flag==1 || flag==2){
			for(int j=0; j<128; j++){
				h_ADC_High_SFT_Raw[j]->Fill(adc_high_sft[j]);
				h_TDC_LE_SFT_Raw[j]->Fill(tdc_le_sft[j][0]);
				h_TDC_LE_SFT_Raw[j]->SetAxisRange(0, 1500, "X");
			}
		}

		// TARGET
		if(flag==3 || flag==4 || flag==5 || flag==6){
			for(int j=0; j<256; j++){
				h_ADC_High_TARGET_Raw[j]->Fill(adc_high_target[j]);
				h_TDC_LE_TARGET_Raw[j]->Fill(tdc_le_target[j][0]);
				h_TDC_LE_TARGET_Raw[j]->SetAxisRange(0, 1500, "X");
			}
		}

		// TOF1
		if(flag==8){
			for(int j=0; j<12; j++){
				h_ADC_tof1U_Raw[j]->Fill(adc_tof1u[j]);
				h_ADC_tof1D_Raw[j]->Fill(adc_tof1d[j]);
				h_TDC_tof1U_Raw[j]->Fill(tdc_tof1u[j]);
				h_TDC_tof1D_Raw[j]->Fill(tdc_tof1d[j]);
				h_TDC_tof1U_Raw[j]->SetAxisRange(0, 1500, "X");
				h_TDC_tof1D_Raw[j]->SetAxisRange(0, 1500, "X");
			}
		}

		// TOF2
		if(flag==9){
			for(int j=0; j<12; j++){
				h_ADC_tof2AO_Raw[j]->Fill(adc_tof2AO[j]);
				h_ADC_tof2AI_Raw[j]->Fill(adc_tof2AI[j]);
				h_ADC_tof2BO_Raw[j]->Fill(adc_tof2BO[j]);
				h_ADC_tof2BI_Raw[j]->Fill(adc_tof2BI[j]);
				h_TDC_tof2AO_Raw[j]->Fill(tdc_tof2AO[j]);
				h_TDC_tof2AI_Raw[j]->Fill(tdc_tof2AI[j]);
				h_TDC_tof2BO_Raw[j]->Fill(tdc_tof2BO[j]);
				h_TDC_tof2BI_Raw[j]->Fill(tdc_tof2BI[j]);
				h_TDC_tof2AO_Raw[j]->SetAxisRange(0, 2500, "X");
				h_TDC_tof2AI_Raw[j]->SetAxisRange(0, 2500, "X");
				h_TDC_tof2BO_Raw[j]->SetAxisRange(0, 2500, "X");
				h_TDC_tof2BI_Raw[j]->SetAxisRange(0, 2500, "X");
			}
		}
	
		if(flag==10){
			for(int j=0; j<12; j++){
				h_ADC_acu_Raw[j]->Fill(adc_acu[j]);
				h_ADC_acd_Raw[j]->Fill(adc_acd[j]);
				h_TDC_acu_Raw[j]->Fill(tdc_acu[j][0]);
				h_TDC_acd_Raw[j]->Fill(tdc_acd[j][0]);
				h_TDC_acu_Raw[j]->SetAxisRange(0, 3000, "X");
				h_TDC_acd_Raw[j]->SetAxisRange(0, 3000, "X");
			}
		}

		//////////////////////////////////////////////////////////////


		//////////////////////////////////////////////////////////////
		/// APPLY CUTS

		// SFT
		if(flag==1 || flag==2){
			for(int k=0; k<128; k++){
				if((tdc_le_sft[k][0]>=SFT_TDC_min[k] && tdc_le_sft[k][0]<=SFT_TDC_max[k]) ||
				   (tdc_le_sft[k][1]>=SFT_TDC_min[k] && tdc_le_sft[k][1]<=SFT_TDC_max[k]) ||
				   (tdc_le_sft[k][2]>=SFT_TDC_min[k] && tdc_le_sft[k][2]<=SFT_TDC_max[k]) ||
				   (tdc_le_sft[k][3]>=SFT_TDC_min[k] && tdc_le_sft[k][3]<=SFT_TDC_max[k]) ||
				   (tdc_le_sft[k][4]>=SFT_TDC_min[k] && tdc_le_sft[k][4]<=SFT_TDC_max[k]) ||
				   (tdc_le_sft[k][5]>=SFT_TDC_min[k] && tdc_le_sft[k][5]<=SFT_TDC_max[k]) ||
				   (tdc_le_sft[k][6]>=SFT_TDC_min[k] && tdc_le_sft[k][6]<=SFT_TDC_max[k]) ||
				   (tdc_le_sft[k][7]>=SFT_TDC_min[k] && tdc_le_sft[k][7]<=SFT_TDC_max[k]) ||
				   (tdc_le_sft[k][8]>=SFT_TDC_min[k] && tdc_le_sft[k][8]<=SFT_TDC_max[k]) ||
				   (tdc_le_sft[k][9]>=SFT_TDC_min[k] && tdc_le_sft[k][9]<=SFT_TDC_max[k]) ||
				   (tdc_le_sft[k][10]>=SFT_TDC_min[k] && tdc_le_sft[k][10]<=SFT_TDC_max[k]) ||
				   (tdc_le_sft[k][11]>=SFT_TDC_min[k] && tdc_le_sft[k][11]<=SFT_TDC_max[k]) ||
				   (tdc_le_sft[k][12]>=SFT_TDC_min[k] && tdc_le_sft[k][12]<=SFT_TDC_max[k]) ||
				   (tdc_le_sft[k][13]>=SFT_TDC_min[k] && tdc_le_sft[k][13]<=SFT_TDC_max[k]) ||
				   (tdc_le_sft[k][14]>=SFT_TDC_min[k] && tdc_le_sft[k][14]<=SFT_TDC_max[k]) ||
				   (tdc_le_sft[k][15]>=SFT_TDC_min[k] && tdc_le_sft[k][15]<=SFT_TDC_max[k])){
						h_ADC_High_SFT[k]->Fill(adc_high_sft[k]);
				}
			}	
		}		

		// TARGET
		if(flag==3 || flag==4 || flag==5 || flag==6){
			for(int k=0; k<256; k++){
				if((tdc_le_target[k][0]>=TARGET_TDC_min[k] && tdc_le_target[k][0]<=TARGET_TDC_max[k]) ||
				   (tdc_le_target[k][1]>=TARGET_TDC_min[k] && tdc_le_target[k][1]<=TARGET_TDC_max[k]) ||
				   (tdc_le_target[k][2]>=TARGET_TDC_min[k] && tdc_le_target[k][2]<=TARGET_TDC_max[k]) ||
				   (tdc_le_target[k][3]>=TARGET_TDC_min[k] && tdc_le_target[k][3]<=TARGET_TDC_max[k]) ||
				   (tdc_le_target[k][4]>=TARGET_TDC_min[k] && tdc_le_target[k][4]<=TARGET_TDC_max[k]) ||
				   (tdc_le_target[k][5]>=TARGET_TDC_min[k] && tdc_le_target[k][5]<=TARGET_TDC_max[k]) ||
				   (tdc_le_target[k][6]>=TARGET_TDC_min[k] && tdc_le_target[k][6]<=TARGET_TDC_max[k]) ||
				   (tdc_le_target[k][7]>=TARGET_TDC_min[k] && tdc_le_target[k][7]<=TARGET_TDC_max[k]) ||
				   (tdc_le_target[k][8]>=TARGET_TDC_min[k] && tdc_le_target[k][8]<=TARGET_TDC_max[k]) ||
				   (tdc_le_target[k][9]>=TARGET_TDC_min[k] && tdc_le_target[k][9]<=TARGET_TDC_max[k]) ||
				   (tdc_le_target[k][10]>=TARGET_TDC_min[k] && tdc_le_target[k][10]<=TARGET_TDC_max[k]) ||
				   (tdc_le_target[k][11]>=TARGET_TDC_min[k] && tdc_le_target[k][11]<=TARGET_TDC_max[k]) ||
				   (tdc_le_target[k][12]>=TARGET_TDC_min[k] && tdc_le_target[k][12]<=TARGET_TDC_max[k]) ||
				   (tdc_le_target[k][13]>=TARGET_TDC_min[k] && tdc_le_target[k][13]<=TARGET_TDC_max[k]) ||
				   (tdc_le_target[k][14]>=TARGET_TDC_min[k] && tdc_le_target[k][14]<=TARGET_TDC_max[k]) ||
				   (tdc_le_target[k][15]>=TARGET_TDC_min[k] && tdc_le_target[k][15]<=TARGET_TDC_max[k])){
						h_ADC_High_TARGET[k]->Fill(adc_high_target[k]);
				}
			}
		}	

		// TOF1
		if(flag==8){
			for(int k=0; k<12; k++){
				if(tdc_tof1u[k]>=TOF1_TDC_min[k] && tdc_tof1u[k]<=TOF1_TDC_max[k]){
					h_ADC_tof1U[k]->Fill(adc_tof1u[k]);
				}

				if(tdc_tof1d[k]>=TOF1_TDC_min[k+12] && tdc_tof1d[k]<=TOF1_TDC_max[k+12]){
					h_ADC_tof1D[k]->Fill(adc_tof1d[k]);
				}
			}
		}

		// TOF2
		if(flag==9){
			for(int k=0; k<12; k++){
				if(tdc_tof2AO[k]>=TOF2_TDC_min[k] && tdc_tof2AO[k]<=TOF2_TDC_max[k]){
					h_ADC_tof2AO[k]->Fill(adc_tof2AO[k]);
				}

				if(tdc_tof2AI[k]>=TOF2_TDC_min[k+12] && tdc_tof2AI[k]<=TOF2_TDC_max[k+12]){
					h_ADC_tof2AI[k]->Fill(adc_tof2AI[k]);
				}

				if(tdc_tof2BO[k]>=TOF2_TDC_min[k+24] && tdc_tof2BO[k]<=TOF2_TDC_max[k+24]){
					h_ADC_tof2BO[k]->Fill(adc_tof2BO[k]);
				}

				if(tdc_tof2BI[k]>=TOF2_TDC_min[k+36] && tdc_tof2BI[k]<=TOF2_TDC_max[k+36]){
					h_ADC_tof2BI[k]->Fill(adc_tof2BI[k]);
				}
			}
		}	

		// AC
		if(flag==10){
			for(int k=0; k<12; k++){
				if((tdc_acu[k][0]>=AC_TDC_min[k] && tdc_acu[k][0]<=AC_TDC_max[k]) ||
				   (tdc_acu[k][1]>=AC_TDC_min[k] && tdc_acu[k][1]<=AC_TDC_max[k]) ||
				   (tdc_acu[k][2]>=AC_TDC_min[k] && tdc_acu[k][2]<=AC_TDC_max[k]) ||
				   (tdc_acu[k][3]>=AC_TDC_min[k] && tdc_acu[k][3]<=AC_TDC_max[k]) ||
				   (tdc_acu[k][4]>=AC_TDC_min[k] && tdc_acu[k][4]<=AC_TDC_max[k]) ||
				   (tdc_acu[k][5]>=AC_TDC_min[k] && tdc_acu[k][5]<=AC_TDC_max[k]) ||
				   (tdc_acu[k][6]>=AC_TDC_min[k] && tdc_acu[k][6]<=AC_TDC_max[k]) ||
				   (tdc_acu[k][7]>=AC_TDC_min[k] && tdc_acu[k][7]<=AC_TDC_max[k])){
						h_ADC_acu[k]->Fill(adc_acu[k]);
						h_TDC_acu[k]->Fill(tdc_acu[k][0]);
				}

				if((tdc_acd[k][0]>=AC_TDC_min[k+12] && tdc_acd[k][0]<=AC_TDC_max[k+12]) ||
				   (tdc_acd[k][1]>=AC_TDC_min[k+12] && tdc_acd[k][1]<=AC_TDC_max[k+12]) ||
				   (tdc_acd[k][2]>=AC_TDC_min[k+12] && tdc_acd[k][2]<=AC_TDC_max[k+12]) ||
				   (tdc_acd[k][3]>=AC_TDC_min[k+12] && tdc_acd[k][3]<=AC_TDC_max[k+12]) ||
				   (tdc_acd[k][4]>=AC_TDC_min[k+12] && tdc_acd[k][4]<=AC_TDC_max[k+12]) ||
				   (tdc_acd[k][5]>=AC_TDC_min[k+12] && tdc_acd[k][5]<=AC_TDC_max[k+12]) ||
				   (tdc_acd[k][6]>=AC_TDC_min[k+12] && tdc_acd[k][6]<=AC_TDC_max[k+12]) ||
				   (tdc_acd[k][7]>=AC_TDC_min[k+12] && tdc_acd[k][7]<=AC_TDC_max[k+12])){
						h_ADC_acd[k]->Fill(adc_acd[k]);
						h_TDC_acd[k]->Fill(tdc_acd[k][0]);
				}
			}
		}

		//////////////////////////////////////////////////////////////
	}   // EndLoop over Events

  	//////////////////////////////////////////////////////////////////////////
	/// TLINES

	TLine *SFT_Line_min[128];		TLine *SFT_Line_max[128];
	TLine *TARGET_Line_min[256];	TLine *TARGET_Line_max[256];
	TLine *TOF1U_Line_min[12];		TLine *TOF1U_Line_max[12];
	TLine *TOF1D_Line_min[12];		TLine *TOF1D_Line_max[12];
	TLine *TOF2AO_Line_min[12];		TLine *TOF2AO_Line_max[12];
	TLine *TOF2AI_Line_min[12];		TLine *TOF2AI_Line_max[12];
	TLine *TOF2BO_Line_min[12];		TLine *TOF2BO_Line_max[12];
	TLine *TOF2BI_Line_min[12];		TLine *TOF2BI_Line_max[12];
	TLine *ACU_Line_min[12];		TLine *ACU_Line_max[12];
	TLine *ACD_Line_min[12];		TLine *ACD_Line_max[12];


	// SFT
	if(flag==1 || flag==2){
		Double_t SFT_Histo_max[128];
		for(int i=0; i<128; i++){
			SFT_Histo_max[i] = h_TDC_LE_SFT_Raw[i]->GetMaximum();

			SFT_Line_min[i] = new TLine(SFT_TDC_min[i], 0, SFT_TDC_min[i], SFT_Histo_max[i]);
			SFT_Line_max[i] = new TLine(SFT_TDC_max[i], 0, SFT_TDC_max[i], SFT_Histo_max[i]);

			SFT_Line_min[i]->SetLineColor(2);
			SFT_Line_max[i]->SetLineColor(2);
		}
	}

	// TARGET
	if(flag==3 || flag==4 || flag==5 || flag==6){
		Double_t TARGET_Histo_max[256];
		for(int i=0; i<256; i++){
			TARGET_Histo_max[i] = h_TDC_LE_TARGET_Raw[i]->GetMaximum();

			TARGET_Line_min[i] = new TLine(TARGET_TDC_min[i], 0, TARGET_TDC_min[i], TARGET_Histo_max[i]);
			TARGET_Line_max[i] = new TLine(TARGET_TDC_max[i], 0, TARGET_TDC_max[i], TARGET_Histo_max[i]);

			TARGET_Line_min[i]->SetLineColor(2);
			TARGET_Line_max[i]->SetLineColor(2);
		}
	}

	// TOF1
	if(flag==8){
		Double_t TOF1U_Histo_max[12];
		Double_t TOF1D_Histo_max[12];
		for(int i=0; i<12; i++){
			TOF1U_Histo_max[i] = h_TDC_tof1U_Raw[i]->GetMaximum();
			TOF1D_Histo_max[i] = h_TDC_tof1D_Raw[i]->GetMaximum();

			TOF1U_Line_min[i] = new TLine(TOF1_TDC_min[i], 0, TOF1_TDC_min[i], TOF1U_Histo_max[i]);
			TOF1U_Line_max[i] = new TLine(TOF1_TDC_max[i], 0, TOF1_TDC_max[i], TOF1U_Histo_max[i]);
			TOF1D_Line_min[i] = new TLine(TOF1_TDC_min[i+12], 0, TOF1_TDC_min[i+12], TOF1D_Histo_max[i]);
			TOF1D_Line_max[i] = new TLine(TOF1_TDC_max[i+12], 0, TOF1_TDC_max[i+12], TOF1D_Histo_max[i]);

			TOF1U_Line_min[i]->SetLineColor(2);
			TOF1U_Line_max[i]->SetLineColor(2);
			TOF1D_Line_min[i]->SetLineColor(2);
			TOF1D_Line_max[i]->SetLineColor(2);
		}
	}

	// TOF2
	if(flag==9){
		Double_t TOF2AO_Histo_max[12];
		Double_t TOF2AI_Histo_max[12];
		Double_t TOF2BO_Histo_max[12];
		Double_t TOF2BI_Histo_max[12];
		for(int i=0; i<12; i++){
			TOF2AO_Histo_max[i] = h_TDC_tof2AO_Raw[i]->GetMaximum();
			TOF2AI_Histo_max[i] = h_TDC_tof2AI_Raw[i]->GetMaximum();
			TOF2BO_Histo_max[i] = h_TDC_tof2BO_Raw[i]->GetMaximum();
			TOF2BI_Histo_max[i] = h_TDC_tof2BI_Raw[i]->GetMaximum();

			TOF2AO_Line_min[i] = new TLine(TOF2_TDC_min[i], 0, TOF2_TDC_min[i], TOF2AO_Histo_max[i]);
			TOF2AO_Line_max[i] = new TLine(TOF2_TDC_max[i], 0, TOF2_TDC_max[i], TOF2AO_Histo_max[i]);
			TOF2AI_Line_min[i] = new TLine(TOF2_TDC_min[i+12], 0, TOF2_TDC_min[i+12], TOF2AI_Histo_max[i]);
			TOF2AI_Line_max[i] = new TLine(TOF2_TDC_max[i+12], 0, TOF2_TDC_max[i+12], TOF2AI_Histo_max[i]);
			TOF2BO_Line_min[i] = new TLine(TOF2_TDC_min[i+24], 0, TOF2_TDC_min[i+24], TOF2BO_Histo_max[i]);
			TOF2BO_Line_max[i] = new TLine(TOF2_TDC_max[i+24], 0, TOF2_TDC_max[i+24], TOF2BO_Histo_max[i]);
			TOF2BI_Line_min[i] = new TLine(TOF2_TDC_min[i+36], 0, TOF2_TDC_min[i+36], TOF2BI_Histo_max[i]);
			TOF2BI_Line_max[i] = new TLine(TOF2_TDC_max[i+36], 0, TOF2_TDC_max[i+36], TOF2BI_Histo_max[i]);

			TOF2AO_Line_min[i]->SetLineColor(2);
			TOF2AO_Line_max[i]->SetLineColor(2);
			TOF2AI_Line_min[i]->SetLineColor(2);
			TOF2AI_Line_max[i]->SetLineColor(2);
			TOF2BO_Line_min[i]->SetLineColor(2);
			TOF2BO_Line_max[i]->SetLineColor(2);
			TOF2BI_Line_min[i]->SetLineColor(2);
			TOF2BI_Line_max[i]->SetLineColor(2);
		}
	}

	// AC
	if(flag==10){
		Double_t ACU_Histo_max[12];
		Double_t ACD_Histo_max[12];
		for(int i=0; i<12; i++){
			ACU_Histo_max[i] = h_TDC_acu_Raw[i]->GetMaximum();
			ACD_Histo_max[i] = h_TDC_acd_Raw[i]->GetMaximum();

			ACU_Line_min[i] = new TLine(AC_TDC_min[i], 0, AC_TDC_min[i], ACU_Histo_max[i]);
			ACU_Line_max[i] = new TLine(AC_TDC_max[i], 0, AC_TDC_max[i], ACU_Histo_max[i]);
			ACD_Line_min[i] = new TLine(AC_TDC_min[i+12], 0, AC_TDC_min[i+12], ACD_Histo_max[i]);
			ACD_Line_max[i] = new TLine(AC_TDC_max[i+12], 0, AC_TDC_max[i+12], ACD_Histo_max[i]);

			ACU_Line_min[i]->SetLineColor(2);
			ACU_Line_max[i]->SetLineColor(2);
			ACD_Line_min[i]->SetLineColor(2);
			ACD_Line_max[i]->SetLineColor(2);
		}
	}

  	//////////////////////////////////////////////////////////////////////////


	////////////////////////////////////////////////////////////////////////
	/// DRAWING
	
	TPaveStats *stats1;
	TPaveStats *stats2;


	// SFT
	if(flag==1 || flag==2){
		TCanvas *c_SFT_TDC;
		c_SFT_TDC = new TCanvas(Name_Can_TDC_LE_SFT, Title_Can_TDC_LE_SFT,1200,500);
		c_SFT_TDC->Divide(8,8);
		c_SFT_TDC->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		if(flag==1){
			for(int i=0; i<64; i++){
				c_SFT_TDC->cd(i+1);
				h_TDC_LE_SFT_Raw[i]->Draw();
				SFT_Line_min[i] ->Draw();
				SFT_Line_max[i]->Draw();
			}
		}

		if(flag==2){
			for(int i=0; i<64; i++){
				c_SFT_TDC->cd(i+1);
				h_TDC_LE_SFT_Raw[i+64]->Draw();
				SFT_Line_min[i+64]->Draw();
				SFT_Line_max[i+64]->Draw();
			}
		}

		TCanvas *c_SFT_ADC;
		c_SFT_ADC = new TCanvas(Name_Can_ADC_High_SFT, Title_Can_ADC_High_SFT,500,-300,1200,500);
		c_SFT_ADC->Divide(8,8);
		c_SFT_ADC->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		
		if(flag==1){
			for(int i=0; i<64; i++){
				c_SFT_ADC->cd(i+1)->SetLogy();
				h_ADC_High_SFT_Raw[i]->Draw();

			    c_SFT_ADC->Update();
        		stats1 = (TPaveStats*)(h_ADC_High_SFT_Raw[i]->GetListOfFunctions()->FindObject("stats"));
        		stats1->SetX1NDC(.55);
        		stats1->SetX2NDC(.75);

				h_ADC_High_SFT[i]->Draw("sames");
			    c_SFT_ADC->Update();
        		stats2 = (TPaveStats*)(h_ADC_High_SFT[i]->GetListOfFunctions()->FindObject("stats"));
  				stats2->SetTextColor(4);
			}
		}	

		if(flag==2){
			for(int i=0; i<64; i++){
				c_SFT_ADC->cd(i+1)->SetLogy();
				h_ADC_High_SFT_Raw[i+64]->Draw();

				c_SFT_ADC->Update();
        		stats1 = (TPaveStats*)(h_ADC_High_SFT_Raw[i+64]->GetListOfFunctions()->FindObject("stats"));
        		stats1->SetX1NDC(.55);
        		stats1->SetX2NDC(.75);

				h_ADC_High_SFT[i+64]->Draw("sames");
				c_SFT_ADC->Update();
        		stats2 = (TPaveStats*)(h_ADC_High_SFT[i+64]->GetListOfFunctions()->FindObject("stats"));
  				stats2->SetTextColor(4);

			}
		}	
	}	


	// TARGET
	if(flag==3 || flag==4 || flag==5 || flag==6){
		TCanvas *c_TARGET_TDC;
		c_TARGET_TDC = new TCanvas(Name_Can_TDC_LE_TARGET, Title_Can_TDC_LE_TARGET,1200,500);
		c_TARGET_TDC->Divide(8,8);
		c_TARGET_TDC->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
	
		if(flag==3){
			for(int i=0; i<64; i++){
				c_TARGET_TDC->cd(i+1);
				h_TDC_LE_TARGET_Raw[i]->Draw();
				TARGET_Line_min[i]->Draw();
				TARGET_Line_max[i]->Draw();
			}
		}

		if(flag==4){
			for(int i=0; i<64; i++){
				c_TARGET_TDC->cd(i+1);
				h_TDC_LE_TARGET_Raw[i+64]->Draw();
				TARGET_Line_min[i+64]->Draw();
				TARGET_Line_max[i+64]->Draw();
			}
		}

		if(flag==5){
			for(int i=0; i<64; i++){
				c_TARGET_TDC->cd(i+1);
				h_TDC_LE_TARGET_Raw[i+128]->Draw();
				TARGET_Line_min[i+128]->Draw();
				TARGET_Line_max[i+128]->Draw();
			}
		}

		if(flag==6){
			for(int i=0; i<64; i++){
				c_TARGET_TDC->cd(i+1);
				h_TDC_LE_TARGET_Raw[i+192]->Draw();
				TARGET_Line_min[i+192]->Draw();
				TARGET_Line_max[i+192]->Draw();
			}
		}


		TCanvas *c_TARGET_ADC;
		c_TARGET_ADC = new TCanvas(Name_Can_ADC_High_TARGET, Title_Can_ADC_High_TARGET,500,-300,1200,500);
		c_TARGET_ADC->Divide(8,8);
		c_TARGET_ADC->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
	
		if(flag==3){
			for(int i=0; i<64; i++){
				c_TARGET_ADC->cd(i+1)->SetLogy();
				h_ADC_High_TARGET_Raw[i]->Draw();

				c_TARGET_ADC->Update();
        		stats1 = (TPaveStats*)(h_ADC_High_TARGET_Raw[i]->GetListOfFunctions()->FindObject("stats"));
        		stats1->SetX1NDC(.55);
        		stats1->SetX2NDC(.75);

				h_ADC_High_TARGET[i]->Draw("sames");
			    c_TARGET_ADC->Update();
        		stats2 = (TPaveStats*)(h_ADC_High_TARGET[i]->GetListOfFunctions()->FindObject("stats"));
  				stats2->SetTextColor(4);
			}
		}

		if(flag==4){
			for(int i=0; i<64; i++){
				c_TARGET_ADC->cd(i+1)->SetLogy();
				h_ADC_High_TARGET_Raw[i+64]->Draw();
	
				c_TARGET_ADC->Update();
        		stats1 = (TPaveStats*)(h_ADC_High_TARGET_Raw[i+64]->GetListOfFunctions()->FindObject("stats"));
        		stats1->SetX1NDC(.55);
        		stats1->SetX2NDC(.75);

				h_ADC_High_TARGET[i+64]->Draw("sames");
			    c_TARGET_ADC->Update();
        		stats2 = (TPaveStats*)(h_ADC_High_TARGET[i+64]->GetListOfFunctions()->FindObject("stats"));
  				stats2->SetTextColor(4);
			}
		}

		if(flag==5){
			for(int i=0; i<64; i++){
				c_TARGET_ADC->cd(i+1)->SetLogy();
				h_ADC_High_TARGET_Raw[i+128]->Draw();
			
				c_TARGET_ADC->Update();
        		stats1 = (TPaveStats*)(h_ADC_High_TARGET_Raw[i+128]->GetListOfFunctions()->FindObject("stats"));
        		stats1->SetX1NDC(.55);
        		stats1->SetX2NDC(.75);

				h_ADC_High_TARGET[i+128]->Draw("sames");
			    c_TARGET_ADC->Update();
        		stats2 = (TPaveStats*)(h_ADC_High_TARGET[i+128]->GetListOfFunctions()->FindObject("stats"));
  				stats2->SetTextColor(4);
			}
		}

		if(flag==6){
			for(int i=0; i<64; i++){
				c_TARGET_ADC->cd(i+1)->SetLogy();
				h_ADC_High_TARGET_Raw[i+192]->Draw();
			
				c_TARGET_ADC->Update();
        		stats1 = (TPaveStats*)(h_ADC_High_TARGET_Raw[i+192]->GetListOfFunctions()->FindObject("stats"));
        		stats1->SetX1NDC(.55);
        		stats1->SetX2NDC(.75);
			
				h_ADC_High_TARGET[i+192]->Draw("sames");
			    c_TARGET_ADC->Update();
        		stats2 = (TPaveStats*)(h_ADC_High_TARGET[i+192]->GetListOfFunctions()->FindObject("stats"));
  				stats2->SetTextColor(4);
			}
		}
	}

	// TOF1
	if(flag==8){
		TCanvas *c_TOF1;
		c_TOF1 = new TCanvas(Name_Can_TOF1, Title_Can_TOF1,1200,500);
		c_TOF1->Divide(12,4);
		c_TOF1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
	
		for(int i=0; i<12; i++){
			c_TOF1->cd(i+1)->SetLogy();
			h_ADC_tof1U_Raw[i]->Draw();
			
			c_TOF1->Update();
        	stats1 = (TPaveStats*)(h_ADC_tof1U_Raw[i]->GetListOfFunctions()->FindObject("stats"));
        	stats1->SetX1NDC(.55);
        	stats1->SetX2NDC(.75);

			h_ADC_tof1U[i]->Draw("sames");
			c_TOF1->Update();
        	stats2 = (TPaveStats*)(h_ADC_tof1U[i]->GetListOfFunctions()->FindObject("stats"));
  			stats2->SetTextColor(4);
		}

		for(int i=0; i<12; i++){
			c_TOF1->cd(i+1+12);
			h_TDC_tof1U_Raw[i]->Draw();
			TOF1U_Line_min[i]->Draw();
			TOF1U_Line_max[i]->Draw();
		}

		for(Int_t i=0; i<12; i++){
			c_TOF1->cd(i+1+24)->SetLogy();
			h_ADC_tof1D_Raw[i]->Draw();
			
			c_TOF1->Update();
        	stats1 = (TPaveStats*)(h_ADC_tof1D_Raw[i]->GetListOfFunctions()->FindObject("stats"));
        	stats1->SetX1NDC(.55);
        	stats1->SetX2NDC(.75);

			h_ADC_tof1D[i]->Draw("sames");
			c_TOF1->Update();
        	stats2 = (TPaveStats*)(h_ADC_tof1D[i]->GetListOfFunctions()->FindObject("stats"));
  			stats2->SetTextColor(4);
		}

		for(Int_t i=0; i<12; i++){
			c_TOF1->cd(i+1+36);
			h_TDC_tof1D_Raw[i]->Draw();
			TOF1D_Line_min[i]->Draw();
			TOF1D_Line_max[i]->Draw();
		}
	}

	// TOF2
	if(flag==9){
		TCanvas *c_TOF2_2;
		c_TOF2_2 = new TCanvas(Name_Can_TOF2_2, Title_Can_TOF2_2,1200,500);
		c_TOF2_2->Divide(6,8);
		c_TOF2_2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(int i=0; i<6; i++){
			c_TOF2_2->cd(i+1)->SetLogy();
			h_ADC_tof2AO_Raw[i+6]->Draw();
			
			c_TOF2_2->Update();
        	stats1 = (TPaveStats*)(h_ADC_tof2AO_Raw[i+6]->GetListOfFunctions()->FindObject("stats"));
        	stats1->SetX1NDC(.55);
        	stats1->SetX2NDC(.75);

			h_ADC_tof2AO[i+6]->Draw("sames");
			c_TOF2_2->Update();
        	stats2 = (TPaveStats*)(h_ADC_tof2AO[i+6]->GetListOfFunctions()->FindObject("stats"));
  			stats2->SetTextColor(4);
		}

		for(int i=0; i<6; i++){
			c_TOF2_2->cd(i+1+6);
			h_TDC_tof2AO_Raw[i+6]->Draw();
			TOF2AO_Line_min[i+6]->Draw();
			TOF2AO_Line_max[i+6]->Draw();
		}

		for(int i=0; i<6; i++){
			c_TOF2_2->cd(i+1+12)->SetLogy();
			h_ADC_tof2AI_Raw[i+6]->Draw();
			
			c_TOF2_2->Update();
        	stats1 = (TPaveStats*)(h_ADC_tof2AI_Raw[i+6]->GetListOfFunctions()->FindObject("stats"));
        	stats1->SetX1NDC(.55);
        	stats1->SetX2NDC(.75);

			h_ADC_tof2AI[i+6]->Draw("sames");
			c_TOF2_2->Update();
        	stats2 = (TPaveStats*)(h_ADC_tof2AI[i+6]->GetListOfFunctions()->FindObject("stats"));
  			stats2->SetTextColor(4);
		}

		for(int i=0; i<6; i++){
			c_TOF2_2->cd(i+1+18);
			h_TDC_tof2AI_Raw[i+6]->Draw();
			TOF2AI_Line_min[i+6]->Draw();
			TOF2AI_Line_max[i+6]->Draw();
		}

		for(int i=0; i<6; i++){
			c_TOF2_2->cd(i+1+24)->SetLogy();
			h_ADC_tof2BO_Raw[i+6]->Draw();
			
			c_TOF2_2->Update();
        	stats1 = (TPaveStats*)(h_ADC_tof2BO_Raw[i+6]->GetListOfFunctions()->FindObject("stats"));
        	stats1->SetX1NDC(.55);
        	stats1->SetX2NDC(.75);

			h_ADC_tof2BO[i+6]->Draw("sames");
			c_TOF2_2->Update();
        	stats2 = (TPaveStats*)(h_ADC_tof2BO[i+6]->GetListOfFunctions()->FindObject("stats"));
  			stats2->SetTextColor(4);
		}

		for(int i=0; i<6; i++){
			c_TOF2_2->cd(i+1+30);
			h_TDC_tof2BO_Raw[i+6]->Draw();
			TOF2BO_Line_min[i+6]->Draw();
			TOF2BO_Line_max[i+6]->Draw();
		}

		for(int i=0; i<6; i++){
			c_TOF2_2->cd(i+1+36)->SetLogy();
			h_ADC_tof2BI_Raw[i+6]->Draw();
			
			c_TOF2_2->Update();
        	stats1 = (TPaveStats*)(h_ADC_tof2BI_Raw[i+6]->GetListOfFunctions()->FindObject("stats"));
        	stats1->SetX1NDC(.55);
        	stats1->SetX2NDC(.75);

			h_ADC_tof2BI[i+6]->Draw("sames");
			c_TOF2_2->Update();
        	stats2 = (TPaveStats*)(h_ADC_tof2BI[i+6]->GetListOfFunctions()->FindObject("stats"));
  			stats2->SetTextColor(4);
		}

		for(int i=0; i<6; i++){
			c_TOF2_2->cd(i+1+42);
			h_TDC_tof2BI_Raw[i+6]->Draw();
			TOF2BI_Line_min[i+6]->Draw();
			TOF2BI_Line_max[i+6]->Draw();
		}

		TCanvas *c_TOF2_1;
		c_TOF2_1 = new TCanvas(Name_Can_TOF2, Title_Can_TOF2,500,-300,1200,500);
		c_TOF2_1->Divide(6,8);
		c_TOF2_1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(int i=0; i<6; i++){
			c_TOF2_1->cd(i+1)->SetLogy();
			h_ADC_tof2AO_Raw[i]->Draw();
			
			c_TOF2_1->Update();
        	stats1 = (TPaveStats*)(h_ADC_tof2AO_Raw[i]->GetListOfFunctions()->FindObject("stats"));
        	stats1->SetX1NDC(.55);
        	stats1->SetX2NDC(.75);

			h_ADC_tof2AO[i]->Draw("sames");
			c_TOF2_1->Update();
        	stats2 = (TPaveStats*)(h_ADC_tof2AO[i]->GetListOfFunctions()->FindObject("stats"));
  			stats2->SetTextColor(4);
		}

		for(int i=0; i<6; i++){
			c_TOF2_1->cd(i+1+6);
			h_TDC_tof2AO_Raw[i]->Draw();
			TOF2AO_Line_min[i]->Draw();
			TOF2AO_Line_max[i]->Draw();
		}

		for(int i=0; i<6; i++){
			c_TOF2_1->cd(i+1+12)->SetLogy();
			h_ADC_tof2AI_Raw[i]->Draw();
			
			c_TOF2_1->Update();
        	stats1 = (TPaveStats*)(h_ADC_tof2AI_Raw[i]->GetListOfFunctions()->FindObject("stats"));
        	stats1->SetX1NDC(.55);
        	stats1->SetX2NDC(.75);

			h_ADC_tof2AI[i]->Draw("sames");
			c_TOF2_1->Update();
        	stats2 = (TPaveStats*)(h_ADC_tof2AI[i]->GetListOfFunctions()->FindObject("stats"));
  			stats2->SetTextColor(4);
		}

		for(int i=0; i<6; i++){
			c_TOF2_1->cd(i+1+18);
			h_TDC_tof2AI_Raw[i]->Draw();
			TOF2AI_Line_min[i]->Draw();
			TOF2AI_Line_max[i]->Draw();
		}

		for(int i=0; i<6; i++){
			c_TOF2_1->cd(i+1+24)->SetLogy();
			h_ADC_tof2BO_Raw[i]->Draw();
			
			c_TOF2_1->Update();
        	stats1 = (TPaveStats*)(h_ADC_tof2BO_Raw[i]->GetListOfFunctions()->FindObject("stats"));
        	stats1->SetX1NDC(.55);
        	stats1->SetX2NDC(.75);

			h_ADC_tof2BO[i]->Draw("sames");
			c_TOF2_1->Update();
        	stats2 = (TPaveStats*)(h_ADC_tof2BO[i]->GetListOfFunctions()->FindObject("stats"));
  			stats2->SetTextColor(4);
		}

		for(int i=0; i<6; i++){
			c_TOF2_1->cd(i+1+30);
			h_TDC_tof2BO_Raw[i]->Draw();
			TOF2BO_Line_min[i]->Draw();
			TOF2BO_Line_max[i]->Draw();
		}

		for(int i=0; i<6; i++){
			c_TOF2_1->cd(i+1+36)->SetLogy();
			h_ADC_tof2BI_Raw[i]->Draw();
			
			c_TOF2_1->Update();
        	stats1 = (TPaveStats*)(h_ADC_tof2BI_Raw[i]->GetListOfFunctions()->FindObject("stats"));
        	stats1->SetX1NDC(.55);
        	stats1->SetX2NDC(.75);

			h_ADC_tof2BI[i]->Draw("sames");
			c_TOF2_1->Update();
        	stats2 = (TPaveStats*)(h_ADC_tof2BI[i]->GetListOfFunctions()->FindObject("stats"));
  			stats2->SetTextColor(4);
		}

		for(int i=0; i<6; i++){
			c_TOF2_1->cd(i+1+42);
			h_TDC_tof2BI_Raw[i]->Draw();
			TOF2BI_Line_min[i]->Draw();
			TOF2BI_Line_max[i]->Draw();
		}
	}

	// AC
	if(flag==10){
		TCanvas *c_AC;
		c_AC = new TCanvas(Name_Can_AC,Title_Can_AC,1300,450); 
		c_AC->Divide(12,4);
		c_AC->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

		for(int i=0; i<12; i++){
			c_AC->cd(i+1)->SetLogy();
			h_ADC_acu_Raw[i]->Draw();
			
			c_AC->Update();
        	stats1 = (TPaveStats*)(h_ADC_acu_Raw[i]->GetListOfFunctions()->FindObject("stats"));
        	stats1->SetX1NDC(.55);
        	stats1->SetX2NDC(.75);

			h_ADC_acu[i]->Draw("sames");
			c_AC->Update();
        	stats2 = (TPaveStats*)(h_ADC_acu[i]->GetListOfFunctions()->FindObject("stats"));
  			stats2->SetTextColor(4);
		}

		for(int i=0; i<12; i++){
			c_AC->cd(i+1+12);
			h_TDC_acu_Raw[i]->Draw();
			ACU_Line_min[i]->Draw();
			ACU_Line_max[i]->Draw();
		//	h_TDC_acu[i]->Draw();
		}

		for(int i=0; i<12; i++){
			c_AC->cd(i+1+24)->SetLogy();
			h_ADC_acd_Raw[i]->Draw();
			
			c_AC->Update();
        	stats1 = (TPaveStats*)(h_ADC_acd_Raw[i]->GetListOfFunctions()->FindObject("stats"));
        	stats1->SetX1NDC(.55);
        	stats1->SetX2NDC(.75);

			h_ADC_acd[i]->Draw("sames");
			c_AC->Update();
        	stats2 = (TPaveStats*)(h_ADC_acd[i]->GetListOfFunctions()->FindObject("stats"));
  			stats2->SetTextColor(4);
		}

		for(int i=0; i<12; i++){
			c_AC->cd(i+1+36);
			h_TDC_acd_Raw[i]->Draw();
			ACD_Line_min[i]->Draw();
			ACD_Line_max[i]->Draw();
		//	h_TDC_acd[i]->Draw();
		}
	}

	////////////////////////////////////////////////////////////////////////
}
