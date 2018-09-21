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
#endif
  
void TARGET_Event_Counter(Int_t Run_Number=5, Int_t flag=0, Int_t threshold = 5){ 

gStyle->SetOptStat(1111);

char source_mapping[]="SFT_Mapping_Oct14.txt";  // Mapping file !!!

//Double_t flag_size_TARGET=2;
//Double_t flag_size_SFT=2.4;
//Double_t flag_size_palette=2.3;

//Double_t flag_size_TARGET=1.35;
//Double_t flag_size_SFT=1.3;
//Double_t flag_size_palette=1.6;

//if(flag!="TARGET" && flag!="SFT"){ flag_size_TARGET=1.35; flag_size_SFT=1.3; flag_size_palette=1.6;}

//Int_t ADC_High_corr_max=0;
//Int_t index_max=0;
//Int_t TDC_min_TARGET=0;
//Int_t TDC_max_TARGET=1050;

//Int_t TDC_min_SFT=750;
//Int_t TDC_max_SFT=1500;

//Int_t ADC_TOF1_thr = 400;

bool has_entries = false;
Int_t no_hit_count = 0;
Int_t yes_hit_count = 0;

Int_t hit_count = 0;

char ADC_cut[100];		sprintf(ADC_cut,"(ADC >= %d)",ADC_cut_SFT);

Int_t adc_high_target[256];   		Int_t ADC_High_TARGET[256];    	//Double_t ADC_High_TARGET_corr[256];   	
Int_t adc_low_target[256]; 		Int_t ADC_Low_TARGET[256]; 	//Double_t ADC_Low_TARGET_corr[256];
//Int_t tdc_le_target[256][16];      	Int_t TDC_LE_TARGET[256];    	Double_t TDC_LE_TARGET_corr[256];   	
//Int_t tdc_te_target[256][16]; 	Int_t TDC_TE_TARGET[256]; 	//Double_t TDC_TE_TARGET_corr[256];

//Int_t adc_high_sft[128];   		Int_t ADC_High_SFT[128];    	Double_t ADC_High_SFT_corr[128];   	
//Int_t adc_low_sft[128]; 		Int_t ADC_Low_SFT[128]; 	//Double_t ADC_Low_SFT_corr[128];
//Int_t tdc_le_sft[128][16];      	Int_t TDC_LE_SFT[128];    	//Double_t TDC_LE_SFT_corr[128];   	
//Int_t tdc_te_sft[128][16]; 		Int_t TDC_TE_SFT[128]; 		//Double_t TDC_TE_SFT_corr[128];


char path_input[200];                   char file_mapping[200];
sprintf(path_input,"%s",path_merged);          sprintf(file_mapping,"../Mapping");

char Name_finput[100];
sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);

char par_finput[200];
//sprintf(par_finput,"/triumfcs/trshare/trek/E36/Macros/April_2015/Mapping/%s",source_mapping);
sprintf(par_finput,"%s/%s",file_mapping,source_mapping);

cout << "   " << endl;
cout << Name_finput << endl;

cout << "   " << endl;
cout << "Mapping File:   " << par_finput << endl;
cout << "" << endl;

Int_t par_temp[2][128];
ifstream fdat(par_finput,ios::in);
for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
fdat.close();

//for(Int_t jj=0; jj<128; jj++){
//	cout << par_temp[0][jj] << "   " << par_temp[1][jj] << endl;
//}

/*

TChain *fChain_TARGET = new TChain("Tree_TARGET");						TChain *fChain_SFT = new TChain("Tree_SFT");
fChain_TARGET->Add(Name_finput);										fChain_SFT->Add(Name_finput);
fChain_TARGET->SetMakeClass(1);											fChain_SFT->SetMakeClass(1);

fChain_TARGET->SetBranchAddress("ADC_High_TARGET",adc_high_target);		fChain_SFT->SetBranchAddress("ADC_High_SFT",adc_high_sft);
fChain_TARGET->SetBranchAddress("ADC_Low_TARGET",adc_low_target);		fChain_SFT->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
fChain_TARGET->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);			fChain_SFT->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
fChain_TARGET->SetBranchAddress("TDC_TE_TARGET",tdc_te_target);			fChain_SFT->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

*/

TChain *fChain= new TChain("Tree");		
fChain->Add(Name_finput);		
fChain->SetMakeClass(1);							

fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);		//fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
fChain->SetBranchAddress("ADC_Low_TARGET",adc_low_target);		//fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
//fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);		//fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
//fChain->SetBranchAddress("TDC_TE_TARGET",tdc_te_target);		//fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

//fChain->SetBranchAddress("ADC_TOF1",ADC_tof1);
//fChain->SetBranchAddress("ADC_TOF2",ADC_tof2);

//Int_t nentries;
Int_t nentries = (Int_t)fChain->GetEntries();//		Int_t nentries_SFT = (Int_t)fChain->GetEntries();
//if(nentries_TARGET==nentries_SFT) nentries = nentries_TARGET;
cout << "  " << endl;
cout << "****  Number of events: " << nentries << "  **** " <<endl;
cout << "  " << endl;


cout << "   " << endl;

if (flag!=0) nentries = flag;
for(Int_t i=0; i<nentries; i++){
	fChain->GetEntry(i);	//	fChain_SFT->GetEntry(i);
	has_entries = false;
	hit_count = 0;

	if(i%10000==1)	cout<<"**** "<<i<<" events done"<<endl;
  	
  	for (int j_TARGET=0; j_TARGET<256; j_TARGET++){
  		ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-ADC_cut_TARGET;
		ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-ADC_cut_TARGET;
  	
    //	cout << j_TARGET << "   " << ADC_High_TARGET[j_TARGET] << "   " << TDC_LE_TARGET[j_TARGET] << endl;
  	}	

	for(int icol=0 ; icol<256 ; icol++){	
		if(ADC_High_TARGET[icol]>0) {
			has_entries = true;
			hit_count++;
		}
	}

	if (hit_count > threshold) {
		yes_hit_count++;
		cout << "Non-Empty: Event " << i << endl;
	}
	else {
		no_hit_count++;
//		cout << "Empty: Event " << i << endl;
	}


}

cout << "" << endl;
cout << "Number of events with hits: " << yes_hit_count << endl;
cout << "Number of events with no hits: " << no_hit_count << endl;

}
