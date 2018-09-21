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
//#include "ListGoodEvents.h"
#endif
  
void Good_Event_Counter(Int_t Run_Number=5, Int_t flag=0, Int_t ievt=0){ 

gStyle->SetOptStat(0);

char* source_mapping="SFT_Mapping_Oct14.txt";  // Mapping file !!!

Int_t adc_high_target[256];   		Int_t ADC_High_TARGET[256];    	Double_t ADC_High_TARGET_corr[256];   	
Int_t adc_low_target[256]; 		Int_t ADC_Low_TARGET[256]; 	//Double_t ADC_Low_TARGET_corr[256];
Int_t tdc_le_target[256][16];      	Int_t TDC_LE_TARGET[256];    	Double_t TDC_LE_TARGET_corr[256];   	
Int_t tdc_te_target[256][16]; 		Int_t TDC_TE_TARGET[256]; 	//Double_t TDC_TE_TARGET_corr[256];

Int_t adc_high_sft[128];   		Int_t ADC_High_SFT[128];    	Double_t ADC_High_SFT_corr[128];   	
Int_t adc_low_sft[128]; 		Int_t ADC_Low_SFT[128]; 	//Double_t ADC_Low_SFT_corr[128];
Int_t tdc_le_sft[128][16];      	Int_t TDC_LE_SFT[128];    	//Double_t TDC_LE_SFT_corr[128];   	
Int_t tdc_te_sft[128][16]; 		Int_t TDC_TE_SFT[128]; 		//Double_t TDC_TE_SFT_corr[128];	

Int_t ADC_TOF1[24];	
Int_t ADC_TOF2[56];

Int_t MwpcADC[512];	Int_t MWPCADC[512];
Int_t fiber[128]={-1};

Int_t notg_counter = 0;
Int_t mwpcok_counter = 0; 
Int_t mwpcng_counter = 0; 
Int_t c2ng_counter = 0;
Int_t c3ng_counter = 0;
Int_t c4ng_counter = 0;


char path_input[200];                   char file_mapping[200];
sprintf(path_input,path_merged);          sprintf(file_mapping,"../Mapping");
//sprintf(path_input,"/triumfcs/trshare/trek/E36/Data/April_2015/root");

char Name_finput[200];
sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);

////

char par_finput[200];
//sprintf(par_finput,"/triumfcs/trshare/trek/E36/Macros/April_2015/Mapping/%s",source_mapping);
sprintf(par_finput,"%s/%s",file_mapping,source_mapping);

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

//for(Int_t jj=0; jj<256; jj++){
//	cout << par_temp3[0][jj] << "   " << par_temp3[1][jj] << endl;
//}

////

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

fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);		fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
fChain->SetBranchAddress("ADC_Low_TARGET",adc_low_target);		fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);		fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
fChain->SetBranchAddress("TDC_TE_TARGET",tdc_te_target);		fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

fChain->SetBranchAddress("ADC_TOF1",ADC_TOF1);
fChain->SetBranchAddress("ADC_TOF2",ADC_TOF2);

fChain->SetBranchAddress("MWPCADC",MwpcADC);

//Int_t nentries;
const Int_t nentries = (Int_t)fChain->GetEntries();//		Int_t nentries_SFT = (Int_t)fChain->GetEntries();
//if(nentries_TARGET==nentries_SFT) nentries = nentries_TARGET;
cout << "  " << endl;
cout << "****  Number of events: " << nentries << "  **** " <<endl;
cout << "  " << endl;


cout << "   " << endl;

//if (flag!=0) nentries = flag;
for(Int_t i=0; i<nentries; i++){
	fChain->GetEntry(i);	//	fChain_SFT->GetEntry(i);

	if(i%10000==1)	cout<<"**** "<<i<<" events done"<<endl;

//	cout << " " << endl;
//	cout << "********* Event " << i << endl;

//	int island[nentries];
//	memset(island,0,nentries*sizeof(int));

	int has_data_TARGET = 0;
	int has_data_SFT = 0;
	int has_data_TOF = 0;
	int has_data_MWPC = 0;

	int C2X_hit = 0;
	int C2Y_hit = 0;
	int C3X_hit = 0;
	int C3Y_hit = 0;
	int C4X_hit = 0;
	int C4Y_hit = 0;
  	
  	for (int j_TARGET=0; j_TARGET<256; j_TARGET++){
  		ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-par_temp_TARGET[1][j_TARGET];
		ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-par_temp_TARGET[1][j_TARGET];
		TDC_LE_TARGET[j_TARGET]=tdc_le_target[j_TARGET][0];
		TDC_TE_TARGET[j_TARGET]=tdc_te_target[j_TARGET][0];
  	
//    	cout << j_TARGET << "   " << ADC_High_TARGET[j_TARGET] << "   " << TDC_LE_TARGET[j_TARGET] << endl;
  	}	


 	for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
		ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-par_temp_SFT[1][j_SFT];
		ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-par_temp_SFT[1][j_SFT];
		TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
		TDC_TE_SFT[j_SFT]=tdc_te_sft[j_SFT][0];
	
	//	cout << j_SFT << "   " << ADC_High_SFT[j_SFT] << "   " << TDC_LE_SFT[j_SFT] << endl;
	}

	for (Int_t j_TOF1=0; j_TOF1<24; j_TOF1++) {
//		cout << j_TOF1 << " -- " << ADC_TOF1[j_TOF1] << endl;
		ADC_TOF1[j_TOF1] = ADC_TOF1[j_TOF1]-par_temp_TOF1[1][j_TOF1];
	}

	for (Int_t j_TOF2=0; j_TOF2<56; j_TOF2++) {
		ADC_TOF2[j_TOF2] = ADC_TOF2[j_TOF2]-par_temp_TOF2[1][j_TOF2];
	}

	for (Int_t j_MWPC=0; j_MWPC<512; j_MWPC++) {
		MWPCADC[j_MWPC] = MwpcADC[j_MWPC]-MWPC_thr;
	}

	///

	for (Int_t p=0; p<256; p++) {
		if (ADC_High_TARGET[p] > 0) has_data_TARGET++;
	}

	for (Int_t p=0; p<128; p++) {
		fiber[p]=ADC_High_SFT[par_temp[1][p]];
		if (fiber[p] > 0) has_data_SFT++;
	}

	for (Int_t p=0; p<512; p++) {
		if (MWPCADC[p] > 0) has_data_MWPC++;
	}

	for (int kk=0; kk<12; kk++) {
		if (kk == 0) {
			if ((ADC_TOF2[0]>0 && ADC_TOF2[24]>0) || (ADC_TOF2[12]>0 && ADC_TOF2[36]>0)) {
				if((ADC_TOF1[0] > 0) || (ADC_TOF1[12] > 0) || (ADC_TOF1[1] > 0) || (ADC_TOF1[13] > 0) || (ADC_TOF1[11] > 0) || (ADC_TOF1[23] > 0)) has_data_TOF++;
			} 
		}

		else if (kk == 6) {
			if ((ADC_TOF2[6]>0 && ADC_TOF2[30]>0) || (ADC_TOF2[55]>0 && ADC_TOF2[42]>0)) {
				if((ADC_TOF1[6] > 0) || (ADC_TOF1[18] > 0) || (ADC_TOF1[7] > 0) || (ADC_TOF1[19] > 0) || (ADC_TOF1[5] > 0) || (ADC_TOF1[17] > 0)) has_data_TOF++;
			}
		}

		else if (kk == 11) {
			if ((ADC_TOF2[11]>0 && ADC_TOF2[35]>0) || (ADC_TOF2[23]>0 && ADC_TOF2[47]>0)) {
				if((ADC_TOF1[11] > 0) || (ADC_TOF1[23] > 0) || (ADC_TOF1[0] > 0) || (ADC_TOF1[12] > 0) || (ADC_TOF1[10] > 0) || (ADC_TOF1[22] > 0)) has_data_TOF++;
			}
		}
	
		else {
			if ((ADC_TOF2[kk]>0 && ADC_TOF2[kk+24]>0) || (ADC_TOF2[kk+12]>0 && ADC_TOF2[kk+36]>0)) {
				if((ADC_TOF1[kk] > 0) || (ADC_TOF1[kk+12] > 0) || (ADC_TOF1[kk+1] > 0) || (ADC_TOF1[kk+13] > 0) || (ADC_TOF1[kk-1] > 0) || (ADC_TOF1[kk+11] > 0)) has_data_TOF++;
			}
		}

	}

	///
	for (Int_t p=0; p<120; p++) {
		if ((p >= 104) && (p <= 111)) continue;
		if (MWPCADC[p] > 0) C2X_hit++;
	}
	
	for (Int_t p=128; p<160; p++) {
		if (MWPCADC[p] > 0) C2Y_hit++;
	}

	for (Int_t p=160; p<288; p++) {
		if (MWPCADC[p] > 0) C3X_hit++;
	}

	for (Int_t p=288; p<320; p++) {
		if (MWPCADC[p] > 0) C3Y_hit++;
	}

	for (Int_t p=321; p<472; p++) {
		if ((p >= 456) && (p <= 463)) continue;
		if (MWPCADC[p] > 0) C4X_hit++;
	}

	for (Int_t p=480; p<512; p++) {
		if (MWPCADC[p] > 0) C4Y_hit++;
	}

//(C2X_hit >= 1 || C2Y_hit >= 1) && (C3X_hit >= 1 || C3Y_hit >= 1) && (C4X_hit >= 1 || C4Y_hit >= 1))

	if (has_data_TOF > 0 && has_data_TARGET == 0) {
		notg_counter++;
		if ((C2X_hit + C2Y_hit >= 1) && (C3X_hit + C3Y_hit >= 1) && (C4X_hit + C4Y_hit >= 1)) {
			mwpcok_counter++;
		}
		else {
			mwpcng_counter++;
		}
		if (C2X_hit == 0 && C2Y_hit == 0) c2ng_counter++;
		if (C3X_hit == 0 && C3Y_hit == 0) c3ng_counter++;
		if (C4X_hit == 0 && C4Y_hit == 0) c4ng_counter++;
	}
		

}

double eff = (double(double(notg_counter)/double(notg_counter)))*100;
double eff2 = (double(double(mwpcok_counter)/double(notg_counter)))*100;
double eff3 = (double(double(mwpcng_counter)/double(notg_counter)))*100;
double eff4 = (double(double(c2ng_counter)/double(notg_counter)))*100;
double eff5 = (double(double(c3ng_counter)/double(notg_counter)))*100;
double eff6 = (double(double(c4ng_counter)/double(notg_counter)))*100;

cout << "" << endl;
cout << "Gap w/Target Empty - notg: " << notg_counter << " -- Percentage of notg events: " << eff << "%" <<endl;
cout << "MWPC C2 C3 & C4 Chambers All Hit: " << mwpcok_counter << " -- Percentage of notg events: " << eff2 << "%" <<endl;
cout << "MWPC Missing Hit in C2 C3 or C4: " << mwpcng_counter << " -- Percentage of notg events: " << eff3 << "%" <<endl;
cout << "No C2 Hits: " << c2ng_counter << " -- Percentage of notg events: " << eff4 << "%" <<endl;
cout << "No C3 Hits: " << c3ng_counter << " -- Percentage of notg events: " << eff5 << "%" <<endl;
cout << "No C4 Hits: " << c4ng_counter << " -- Percentage of notg events: " << eff6 << "%" <<endl;
cout << "" << endl;

TH1D *h_count;	char Title_h_count[100];	char Name_h_count[100];
sprintf(Title_h_count,"Good Event Counter -- Run %d : MWPC > %d", Run_Number, MWPC_thr); 
sprintf(Name_h_count,"Number of Entries"); 
h_count = new TH1D(Name_h_count,Title_h_count,13,0,13);

TH1D *h_count2;	char Title_h_count2[100];	char Name_h_count2[100];
sprintf(Title_h_count2,"Good Event Counter -- Run %d", Run_Number); 
sprintf(Name_h_count2,"Counter 2"); 
h_count2 = new TH1D(Name_h_count2,Title_h_count2,13,0,13);

TH1D *h_count3;	char Title_h_count3[100];	char Name_h_count3[100];
sprintf(Title_h_count3,"Good Event Counter -- Run %d", Run_Number); 
sprintf(Name_h_count3,"Counter 3"); 
h_count3 = new TH1D(Name_h_count3,Title_h_count3,13,0,13);

TH1D *h_count4;	char Title_h_count4[100];	char Name_h_count4[100];
sprintf(Title_h_count4,"Good Event Counter -- Run %d", Run_Number); 
sprintf(Name_h_count4,"Counter 4"); 
h_count4 = new TH1D(Name_h_count4,Title_h_count4,13,0,13);

TH1D *h_count5;	char Title_h_count5[100];	char Name_h_count5[100];
sprintf(Title_h_count5,"Good Event Counter -- Run %d", Run_Number); 
sprintf(Name_h_count5,"Counter 1"); 
h_count5 = new TH1D(Name_h_count5,Title_h_count5,13,0,13);

TH1D *h_count6;	char Title_h_count6[100];	char Name_h_count6[100];
sprintf(Title_h_count6,"Good Event Counter -- Run %d", Run_Number); 
sprintf(Name_h_count6,"Counter 5"); 
h_count6 = new TH1D(Name_h_count6,Title_h_count6,13,0,13);

h_count->Fill(1,notg_counter);
h_count2->Fill(3,mwpcok_counter);
h_count3->Fill(5,mwpcng_counter);
h_count4->Fill(7,c2ng_counter);
h_count5->Fill(9,c3ng_counter);
h_count6->Fill(11,c4ng_counter);

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

char Name_Can2[100];			char Title_Can2[100];
sprintf(Name_Can2,"Good Event Counter -- Run %d", Run_Number);
sprintf(Title_Can2,"Good Event Counter -- Run %d", Run_Number);

TCanvas *c1;
c1 = new TCanvas(Name_Can2,Title_Can2,1000,900); 
c1->Divide(2,3);
c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

TLegend *l1;
l1 = new TLegend(0.75,0.75,0.99,0.99);
l1->AddEntry(h_count, "Gap w/Target Empty - notg", "F");
l1->AddEntry(h_count2, "MWPC C2 C3 & C4 Chambers All Hit", "F");
l1->AddEntry(h_count3, "MWPC Missing Hit in C2 C3 or C4", "F");
l1->AddEntry(h_count4, "No C2 Hits", "F");
l1->AddEntry(h_count5, "No C3 Hits", "F");
l1->AddEntry(h_count6, "No C4 Hits", "F");

h_count->SetLineColor(1);
h_count->SetFillStyle(3004);
h_count->SetFillColor(1);
h_count->Draw();

h_count2->SetLineColor(2);
h_count2->SetFillStyle(3004);
h_count2->SetFillColor(2);
h_count2->Draw("same");

h_count3->SetLineColor(3);
h_count3->SetFillStyle(3004);
h_count3->SetFillColor(3);
h_count3->Draw("same");

h_count4->SetLineColor(4);
h_count4->SetFillStyle(3004);
h_count4->SetFillColor(4);
h_count4->Draw("same");

h_count5->SetLineColor(7);
h_count5->SetFillStyle(3004);
h_count5->SetFillColor(7);
h_count5->Draw("same");

h_count6->SetLineColor(6);
h_count6->SetFillStyle(3004);
h_count6->SetFillColor(6);
h_count6->Draw("same");
l1->Draw("same");

return;

}

