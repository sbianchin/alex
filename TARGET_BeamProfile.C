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
#endif
  
void TARGET_BeamProfile(Int_t Run_Number=5, Int_t flag=0, Int_t ievt=0){ 

gStyle->SetOptStat(111111);

Int_t adc_high_target[256];   		Int_t ADC_High_TARGET[256];    	//Double_t ADC_High_TARGET_corr[256];   	
//Int_t adc_low_target[256]; 		Int_t ADC_Low_TARGET[256]; 	//Double_t ADC_Low_TARGET_corr[256];
Int_t tdc_le_target[256][16];      	Int_t TDC_LE_TARGET[256];    	//Double_t TDC_LE_TARGET_corr[256];   	
//Int_t tdc_te_target[256][16]; 		Int_t TDC_TE_TARGET[256]; 	//Double_t TDC_TE_TARGET_corr[256];

//Int_t adc_high_sft[128];   		Int_t ADC_High_SFT[128];    	//Double_t ADC_High_SFT_corr[128];   	
//Int_t adc_low_sft[128]; 		Int_t ADC_Low_SFT[128]; 	//Double_t ADC_Low_SFT_corr[128];
//Int_t tdc_le_sft[128][16];      	Int_t TDC_LE_SFT[128];    	//Double_t TDC_LE_SFT_corr[128];   	
//Int_t tdc_te_sft[128][16]; 		Int_t TDC_TE_SFT[128]; 		//Double_t TDC_TE_SFT_corr[128];	

//Int_t ADC_tof1[24];	Int_t ADC_TOF1[24];	
//Int_t ADC_tof2[56];	Int_t ADC_TOF2[56];

//Int_t MwpcADC[512];	Int_t MWPCADC[512];
//Int_t fiber[128]={-1};

Int_t TARGET_hit_counter[256] = {0};

Int_t TDC_thr_low = 100;
Int_t TDC_thr_high = 1020;
Int_t ADC_cut_TARGET2 = 2000;

/////////////////////////////////////////////////////////////////////////////////////////////////


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

char footer[100];
sprintf(footer,"Event_Display_MS.C  --  Run %d ; Event %d",Run_Number,ievt);

cout << "   " << endl;
cout << Name_finput << endl;

cout << "   " << endl;
cout << "SFT Mapping File:   " << par_finput << endl;
cout << "" << endl;
cout << "MWPC Mapping File:   " << par_finput2 << endl;
cout << "" << endl;

////

/*
Int_t par_temp[2][128];
ifstream fdat(par_finput,ios::in);
for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
fdat.close();

char par_temp2[512][50];
//char par_temp3[512][50];
ifstream fdat2(par_finput2,ios::in);
for(Int_t ii=0; ii<512; ii++) fdat2 >> par_temp2[ii];
fdat2.close();
*/

Int_t par_temp_TARGET[2][256];
ifstream fdat3(par_finput3,ios::in);
for(Int_t ii=0; ii<256; ii++) fdat3 >> par_temp_TARGET[0][ii] >> par_temp_TARGET[1][ii];
fdat3.close();

/*
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
*/


TChain *fChain= new TChain("Tree");		
fChain->Add(Name_finput);		
fChain->SetMakeClass(1);							

fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);		//fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
//fChain->SetBranchAddress("ADC_Low_TARGET",adc_low_target);		fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);		//fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
//->SetBranchAddress("TDC_TE_TARGET",tdc_te_target);		fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

//fChain->SetBranchAddress("ADC_TOF1",ADC_tof1);
//fChain->SetBranchAddress("ADC_TOF2",ADC_tof2);

//fChain->SetBranchAddress("MWPCADC",MwpcADC);

//Int_t nentries;
Int_t nentries = (Int_t)fChain->GetEntries();//		Int_t nentries_SFT = (Int_t)fChain->GetEntries();
//if(nentries_TARGET==nentries_SFT) nentries = nentries_TARGET;
cout << "  " << endl;
cout << "****  Number of events: " << nentries << "  **** " <<endl;
cout << "  " << endl;


cout << "   " << endl;

if (flag!=0) nentries=flag;
for(Int_t i=0; i<nentries; i++){
//for(Int_t i=0; i<15000; i++){
	fChain->GetEntry(i);	//	fChain_SFT->GetEntry(i);

	if(i%10000==1)	cout<<"**** "<<i<<" events done"<<endl;
//	cout << " " << endl;
//	cout << "********* Event " << i << endl;
  	
  	for (int j_TARGET=0; j_TARGET<256; j_TARGET++){
  		ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-par_temp_TARGET[1][j_TARGET];
		//ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-par_temp_TARGET[1][j_TARGET];
		TDC_LE_TARGET[j_TARGET]=tdc_le_target[j_TARGET][0];
		//TDC_TE_TARGET[j_TARGET]=tdc_te_target[j_TARGET][0];
  	
//   	cout << j_TARGET << "   " << ADC_Low_TARGET[j_TARGET] << endl;
  	}	

  /*
 	for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
		ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-par_temp_SFT[1][j_SFT];
		ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-par_temp_SFT[1][j_SFT];
		TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
		TDC_TE_SFT[j_SFT]=tdc_te_sft[j_SFT][0];
	
	//	cout << j_SFT << "   " << ADC_High_SFT[j_SFT] << "   " << TDC_LE_SFT[j_SFT] << endl;
	}

	for (Int_t j_TOF1=0; j_TOF1<24; j_TOF1++) {
//		cout << j_TOF1 << " -- " << ADC_tof1[j_TOF1] << endl;
		ADC_TOF1[j_TOF1] = ADC_tof1[j_TOF1]-par_temp_TOF1[1][j_TOF1];
	}

	for (Int_t j_TOF2=0; j_TOF2<56; j_TOF2++) {
		ADC_TOF2[j_TOF2] = ADC_tof2[j_TOF2]-par_temp_TOF2[1][j_TOF2];
	}

	for (Int_t j_MWPC=0; j_MWPC<512; j_MWPC++) {
		MWPCADC[j_MWPC] = MwpcADC[j_MWPC]-MWPC_thr;
	}
  */

	for (int j_TARGET=0; j_TARGET<256; j_TARGET++){
  		if ((ADC_High_TARGET[j_TARGET]>ADC_cut_TARGET2) && (TDC_LE_TARGET[j_TARGET] > TDC_thr_low) && (TDC_LE_TARGET[j_TARGET] < TDC_thr_high)) TARGET_hit_counter[j_TARGET]++;
  	}	
}

TH2F *hK_HITS;
char Title_TARGET_h3D[100];	char Name_TARGET_h3D[100];
sprintf(Name_TARGET_h3D,"ADC HG TARGET Hit Positions - Run %d: ADC > %d", Run_Number, ADC_cut_TARGET); 
sprintf(Title_TARGET_h3D,"ADC HG TARGET Hit Positions -  Run %d", Run_Number); 
hK_HITS = new TH2F(Name_TARGET_h3D,Title_TARGET_h3D,20,-29.45,29.45,20,-29.45,29.45);

for(int k=0 ; k<256 ; k++) {
	hK_HITS->Fill(Xloc[k],Yloc[k],TARGET_hit_counter[k]);
}

/*
char Name_Can[100];			char Title_Can[100];
sprintf(Name_Can,"ADC -- Run %d -- Kaon_Centroid_Histo3D.C",Run_Number);
sprintf(Title_Can,"ADC -- Run %d  -- Kaon_Centroid_Histo3D.C",Run_Number);

TCanvas* c1_TARGET = new TCanvas(Name_Can, Title_Can, 500,300,500,500);
c1_TARGET->cd();
gStyle->SetPalette(1);
hK_HITS->Draw("Lego2");
*/

char Name_Can2[100];			char Title_Can2[100];
sprintf(Name_Can2,"ADC -- Run %d -- TARGET_BeamProfile.C",Run_Number);
sprintf(Title_Can2,"ADC -- Run %d  -- TARGET_BeamProfile.C",Run_Number);

TCanvas* c2_TARGET = new TCanvas(Name_Can2, Title_Can2, 700,900,1000,1000);
c2_TARGET->Divide(2,2);
c2_TARGET->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

c2_TARGET->cd(1);
gStyle->SetPalette(1);
hK_HITS->SetMarkerSize(10);
hK_HITS->Draw("COLZ");

c2_TARGET->cd(2);
TH1D *hK_HITSX = hK_HITS->ProjectionX("XProjection");
hK_HITSX->SetTitle("TARGET ADC Hits - X Projection");
hK_HITSX->Draw();

c2_TARGET->cd(3);
TH1D *hK_HITSY = hK_HITS->ProjectionY("YProjection");
hK_HITSY->SetTitle("TARGET ADC Hits - Y Projection");
hK_HITSY->Draw();

c2_TARGET->cd(4);
gStyle->SetPalette(1);
hK_HITS->Draw("Lego2");

}