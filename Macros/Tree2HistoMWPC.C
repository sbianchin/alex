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
  
void Tree2HistoMWPC(Int_t run_number=5, Int_t flag=0, Int_t QUIET=0) { 

gStyle->SetOptStat(11111111);

Int_t adc_high_target[256]; 	
Int_t adc_low_target[256]; 		
Int_t tdc_le_target[256][16];      	  	
Int_t tdc_te_target[256][16]; 		

Int_t ADC_High_sft[128];
Int_t ADC_Low_sft[128];
Int_t TDC_LE_sft[128][16];
Int_t TDC_TE_sft[128][16];

Int_t ADC_TOF1[24];
Int_t ADC_TOF2[56];

Int_t ADC_tof1U[12];  Int_t ADC_TOF1U[12];
Int_t ADC_tof1D[12];  Int_t ADC_TOF1D[12];

Int_t ADC_tof2AO[12];   Int_t ADC_TOF2AO[12];
Int_t ADC_tof2BO[12];   Int_t ADC_TOF2BO[12];
Int_t ADC_tof2AI[12];   Int_t ADC_TOF2AI[12];
Int_t ADC_tof2BI[12];   Int_t ADC_TOF2BI[12];

Int_t MwpcADC[512];

Int_t TOF1_Mapping1[12] = {1,3,5,7,9,11,13,15,17,19,21,23};
Int_t TOF1_Mapping2[12] = {2,4,6,8,10,12,14,16,18,20,22,24};

Int_t TOF2_Mapping1[12] = {1,5,9,13,17,21,25,29,33,37,41,45};
Int_t TOF2_Mapping2[12] = {2,6,10,14,18,22,26,30,34,38,42,46};
Int_t TOF2_Mapping3[12] = {3,7,11,15,19,23,27,31,35,39,43,47};
Int_t TOF2_Mapping4[12] = {4,8,12,16,20,24,28,32,36,40,44,48};	

char path_input[200];                   char file_mapping[200];
sprintf(path_input,"%s",path_merged);          sprintf(file_mapping,"../Mapping");
//sprintf(path_input,"/triumfcs/trshare/trek/E36/Data/April_2015/root");

char Name_finput[200];
sprintf(Name_finput,"%s/Run%dMS.root",path_input, run_number);

char path_output[200];
sprintf(path_output,"%s",path_histo);

char par_finput[200];
sprintf(par_finput,"/data/trek/E36/offline/Macros/Run_Number_%d_MWPC_Events_Cut_%d.txt",run_number,MWPC_thr);

Int_t par_temp[50000] = {-1};
ifstream fdat(par_finput,ios::in);
for(Int_t ii=0; ii<50000; ii++) fdat >> par_temp[ii];
fdat.close();

//for (Int_t ii=0; ii<1000; ii++) cout << par_temp[ii] << endl;

//char footer[100];
//sprintf(footer,"Event_Display_MS.C  --  Run %d ; Event %d",Run_Number,ievt);

cout << "   " << endl;
cout << Name_finput << endl;
cout << "   " << endl;
cout << par_finput << endl;

TChain *fChain= new TChain("Tree");		
fChain->Add(Name_finput);		
fChain->SetMakeClass(1);							

fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);		
fChain->SetBranchAddress("ADC_Low_TARGET",adc_low_target);
fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);		
fChain->SetBranchAddress("TDC_TE_TARGET",tdc_te_target);

fChain->SetBranchAddress("ADC_High_SFT",ADC_High_sft);
fChain->SetBranchAddress("ADC_Low_SFT",ADC_Low_sft);
fChain->SetBranchAddress("TDC_LE_SFT",TDC_LE_sft);
fChain->SetBranchAddress("TDC_TE_SFT",TDC_TE_sft);

//fChain->SetBranchAddress("ADC_TOF1",ADC_TOF1);
//fChain->SetBranchAddress("ADC_TOF2",ADC_TOF2);
fChain->SetBranchAddress("ADC_TOF1U", ADC_tof1U);
fChain->SetBranchAddress("ADC_TOF1D", ADC_tof1D);

fChain->SetBranchAddress("ADC_TOF2AO", ADC_tof2AO);
fChain->SetBranchAddress("ADC_TOF2AI", ADC_tof2AI);
fChain->SetBranchAddress("ADC_TOF2BO", ADC_tof2BO);
fChain->SetBranchAddress("ADC_TOF2BI", ADC_tof2BI);

fChain->SetBranchAddress("MWPCADC",MwpcADC);			

TH1D *h_ADC_tof1[24];   char Title_ADC_tof1[24][100];	char Name_ADC_tof1[24][100];
TH1D *h_ADC_tof2[56];	char Title_ADC_tof2[56][100];	char Name_ADC_tof2[56][100];      


for (int i=0; i<24; i++) {
	sprintf(Title_ADC_tof1[i],"Raw ADC (Ch. %d)  --  TOF1",i);
	sprintf(Name_ADC_tof1[i],"ADC_TOF1 (Ch. %d)",i);
}

for (int i=0; i<56; i++) {
	sprintf(Title_ADC_tof2[i],"Raw ADC (Ch. %d)  --  TOF2",i);
	sprintf(Name_ADC_tof2[i],"ADC_TOF2 (Ch. %d)",i);
}

//create histograms

for (int  i=0; i<24; i++) {
h_ADC_tof1[i] = new TH1D(Name_ADC_tof1[i],Title_ADC_tof1[i],1000,0,4100);
}

for (int  i=0; i<56; i++) {
h_ADC_tof2[i] = new TH1D(Name_ADC_tof2[i],Title_ADC_tof2[i],1000,0,4100);
}

//read all entries and fill the histograms
Int_t nentries = (Int_t)fChain->GetEntries();
cout <<  "Total Number of Entries :     " <<  nentries << endl;
cout << "" << endl;



//if(flag!=0) nentries=flag;
for (Int_t i=0; i<50000; i++) {
fChain->GetEntry(i);
//fChain->GetEntry(par_temp[i]);

	if(i%10000==1) cout<<"**** "<<i<<" events done"<<endl;
	//if ((par_temp[i]>0) || ((par_temp[i]==0) && (i==0))) {
	
		for (Int_t j_TOF1=0; j_TOF1<12; j_TOF1++) {
			ADC_TOF1U[j_TOF1] = ADC_tof1U[j_TOF1];
			ADC_TOF1D[j_TOF1] = ADC_tof1D[j_TOF1];
		}

		for (Int_t j_TOF2=0; j_TOF2<12; j_TOF2++) {
			ADC_TOF2AO[j_TOF2] = ADC_tof2AO[j_TOF2];
			ADC_TOF2AI[j_TOF2] = ADC_tof2AI[j_TOF2];
			ADC_TOF2BO[j_TOF2] = ADC_tof2BO[j_TOF2];
			ADC_TOF2BI[j_TOF2] = ADC_tof2BI[j_TOF2];
		}


		for(int i=0; i<12; i++){
			ADC_TOF1[i] = ADC_TOF1U[i];
			ADC_TOF1[i+12] = ADC_TOF1D[i];
			ADC_TOF2[i] = ADC_TOF2AO[i];
			ADC_TOF2[i+12] = ADC_TOF2BO[i];
			ADC_TOF2[i+24] = ADC_TOF2AI[i];
			ADC_TOF2[i+36] = ADC_TOF2BI[i];
		}


		for (Int_t j=0; j<24; j++) {
			h_ADC_tof1[j]->Fill(ADC_TOF1[j]);
		}

		for (Int_t j=0; j<56; j++) {
			h_ADC_tof2[j]->Fill(ADC_TOF2[j]);
		}
	//}
}

char Name_Can_ADC_tof1[100];			char Title_Can_ADC_tof1[100];
sprintf(Name_Can_ADC_tof1,"TOF1 ADCs -- Run %d  (Ch. 0 - 23) -- Mapping Sorted",run_number);
sprintf(Title_Can_ADC_tof1,"TOF1 ADCs -- Run %d  (Ch. 0 - 23) -- Mapping Sorted",run_number);

char Name_Can_ADC_tof2[100];			char Title_Can_ADC_tof2[100];
sprintf(Name_Can_ADC_tof2,"TOF2 ADCs -- Run %d  (Ch. 0 - 47) -- Mapping Sorted",run_number);
sprintf(Title_Can_ADC_tof2,"TOF2 ADCs -- Run %d  (Ch. 0 - 47) -- Mapping Sorted",run_number);

	TCanvas *c7;
	c7 = new TCanvas(Name_Can_ADC_tof1,Title_Can_ADC_tof1,1200,500); 
	c7->Divide(6,4);
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

	TCanvas *c8;
	c8 = new TCanvas(Name_Can_ADC_tof2,Title_Can_ADC_tof2,1500,800); 
	c8->Divide(8,6);
	c8->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=0; ican<12; ican++){
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

	for(Int_t ican=0; ican<12; ican++){
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
	
	for(Int_t ican=0; ican<12; ican++){
		c8->cd(TOF2_Mapping3[ican])->SetLogy();

		char TOF2Title3[250];
		sprintf(TOF2Title3, "Raw ADC (Ch. %d) - TOF2 | OutB-%d : Run %d", ican+12, ican+1, run_number);
//		h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_xmax,"X");
//		if(TOF2_Graph_ymax > 0) {  
//			h_ADC_tof2[ican]->SetAxisRange(0, TOF2_Graph_ymax,"Y");
//		}
//		if (ican+12 == 18) {
//			sprintf(TOF2Title3, "Raw ADC (Ch. 55) - TOF2 | OutB-7 : Run %d", run_number);
//			h_ADC_tof2[55]->SetTitle(TOF2Title3);
//			h_ADC_tof2[55]->Draw();
//		}
//		else {
			h_ADC_tof2[ican+12]->SetTitle(TOF2Title3);
			h_ADC_tof2[ican+12]->Draw();
//		}
	}

	for(Int_t ican=0; ican<12; ican++){
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

return;

}
