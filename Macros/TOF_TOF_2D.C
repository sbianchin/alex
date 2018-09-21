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
  
void TOF_TOF_2D(Int_t run_number=5, Int_t flag=0, Int_t QUIET=0) { 

gStyle->SetOptStat(1111111);

Int_t TDC_min_TOF1 = 0;
Int_t TDC_max_TOF1 = 4000;
Int_t TDC_min_TOF2 = 0;
Int_t TDC_max_TOF2 = 4000;

///Adjust chart ranges

Int_t TOF1_xmin = 0;
Int_t TOF1_xmax = 400;
Int_t TOF1_ymin = 0;
Int_t TOF1_ymax = 400;

Int_t TOF2_xmin = 0;
Int_t TOF2_xmax = 400;
Int_t TOF2_ymin = 0;
Int_t TOF2_ymax = 400;

char source_mapping[]="SFT_Mapping_Oct14.txt";

Int_t adc_high_target[256]; 	
Int_t adc_low_target[256]; 		
Int_t tdc_le_target[256][16];      	  	
Int_t tdc_te_target[256][16]; 		

Int_t ADC_High_sft[128];
Int_t ADC_Low_sft[128];
Int_t TDC_LE_sft[128][16];
Int_t TDC_TE_sft[128][16];

Int_t ADC_tof1[24];
Int_t ADC_tof2[56];

Int_t ADC_TOF1[24];
Int_t ADC_TOF2[56];

Int_t ADC_tof1U[12];  Int_t ADC_TOF1U[12];
Int_t ADC_tof1D[12];  Int_t ADC_TOF1D[12];

Int_t ADC_tof2AO[12];   Int_t ADC_TOF2AO[12];
Int_t ADC_tof2BO[12];   Int_t ADC_TOF2BO[12];
Int_t ADC_tof2AI[12];   Int_t ADC_TOF2AI[12];
Int_t ADC_tof2BI[12];   Int_t ADC_TOF2BI[12];

Int_t TDC_tof1U[12];
Int_t TDC_tof1D[12];

Int_t TDC_tof2AO[12];
Int_t TDC_tof2BO[12];
Int_t TDC_tof2AI[12];
Int_t TDC_tof2BI[12];

Int_t MwpcADC[512];

//Int_t TOF1_Mapping1[12] = {1,3,5,7,9,11,13,15,17,19,21,23};
//Int_t TOF1_Mapping2[12] = {2,4,6,8,10,12,14,16,18,20,22,24};

//Int_t TOF2_Mapping1[12] = {1,5,9,13,17,21,25,29,33,37,41,45};
//Int_t TOF2_Mapping2[12] = {2,6,10,14,18,22,26,30,34,38,42,46};
// TOF2_Mapping3[12] = {3,7,11,15,19,23,27,31,35,39,43,47};
//Int_t TOF2_Mapping4[12] = {4,8,12,16,20,24,28,32,36,40,44,48};	

char path_input[200];                   char file_mapping[200];
sprintf(path_input,"%s",path_merged);          sprintf(file_mapping,"../Mapping");
//sprintf(path_input,"/triumfcs/trshare/trek/E36/Data/April_2015/root");

char Name_finput[200];
sprintf(Name_finput,"%s/Run%dMS.root",path_input, run_number);

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

char par_finput7[200];
//sprintf(par_finput,"/triumfcs/trshare/trek/E36/Macros/April_2015/Mapping/%s",source_mapping);
sprintf(par_finput7,"%s/TDC_TOF1_Shifts.txt",file_mapping);

char par_finput8[200];
//sprintf(par_finput,"/triumfcs/trshare/trek/E36/Macros/April_2015/Mapping/%s",source_mapping);
sprintf(par_finput8,"%s/TDC_TOF2_Shifts.txt",file_mapping);

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

Int_t par_temp_TDC_TOF1[2][24];
ifstream fdat7(par_finput7,ios::in);
for(Int_t ii=0; ii<24; ii++) fdat7 >> par_temp_TDC_TOF1[0][ii] >> par_temp_TDC_TOF1[1][ii];
fdat7.close();

Int_t par_temp_TDC_TOF2[2][48];
ifstream fdat8(par_finput8,ios::in);
for(Int_t ii=0; ii<48; ii++) fdat8 >> par_temp_TDC_TOF2[0][ii] >> par_temp_TDC_TOF2[1][ii];
fdat8.close();

//char footer[100];
//sprintf(footer,"Event_Display_MS.C  --  Run %d ; Event %d",Run_Number,ievt);


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

//fChain->SetBranchAddress("ADC_TOF1",ADC_tof1);
//fChain->SetBranchAddress("ADC_TOF2",ADC_tof2);
fChain->SetBranchAddress("MWPCADC",MwpcADC);
fChain->SetBranchAddress("ADC_TOF1U", ADC_tof1U);
fChain->SetBranchAddress("ADC_TOF1D", ADC_tof1D);

fChain->SetBranchAddress("ADC_TOF2AO", ADC_tof2AO);
fChain->SetBranchAddress("ADC_TOF2AI", ADC_tof2AI);
fChain->SetBranchAddress("ADC_TOF2BO", ADC_tof2BO);
fChain->SetBranchAddress("ADC_TOF2BI", ADC_tof2BI);

fChain->SetBranchAddress("TDC_TOF1U",TDC_tof1U);
fChain->SetBranchAddress("TDC_TOF1D",TDC_tof1D);
fChain->SetBranchAddress("TDC_TOF2AO",TDC_tof2AO);
fChain->SetBranchAddress("TDC_TOF2AI",TDC_tof2AI);
fChain->SetBranchAddress("TDC_TOF2BO",TDC_tof2BO);
fChain->SetBranchAddress("TDC_TOF2BI",TDC_tof2BI);			

TH2D *h_ADC_TDC_tof1[12];   char Title_ADC_TDC_tof1[12][100];	char Name_ADC_TDC_tof1[12][100];
//TH2D *h_ADC_TDC_tof1_peak[12];   char Title_ADC_TDC_tof1_peak[12][100];	char Name_ADC_TDC_tof1_peak[12][100];

TH2D *h_ADC_TDC_tof2[24];   char Title_ADC_TDC_tof2[24][100];	char Name_ADC_TDC_tof2[24][100];
//TH2D *h_ADC_TDC_tof2_peak[24];   char Title_ADC_TDC_tof2_peak[24][100];	char Name_ADC_TDC_tof2_peak[24][100];

for (int i=0; i<12; i++) {
	sprintf(Title_ADC_TDC_tof1[i],"TDC Down vs. TDC Up (Gap %d)  --  TOF1 | Run %d",i+1, run_number);
//	sprintf(Name_ADC_TDC_tof1[i],"TDC Down vs. TDC Up  (Gap %d) - TOF1",i+1);
	sprintf(Name_ADC_TDC_tof1[i],"TOF1 -- Gap %d  (Run %d)",i+1, run_number);
//	sprintf(Title_ADC_TDC_tof1_peak[i],"ADC vs. TDC (Ch. %d)  --  TOF1 Peak",i);
//	sprintf(Name_ADC_TDC_tof1_peak[i],"ADC vs. TDC (Ch. %d) - TOF1 Peak",i);
}

for (int i=0; i<24; i++) {
	if (i<12) {
		sprintf(Title_ADC_TDC_tof2[i],"TDC Out vs. TDC In  (Gap %dA)  --  TOF2 | Run %d",i+1, run_number);
		sprintf(Name_ADC_TDC_tof2[i],"TOF2 -- Gap %dA  (Run %d)",i+1, run_number);
	}
	else {
		sprintf(Title_ADC_TDC_tof2[i],"TDC Out vs. TDC In (Gap %dB)  --  TOF2 | Run %d",i-11, run_number);
		sprintf(Name_ADC_TDC_tof2[i],"TOF2 -- Gap %dB  (Run %d)",i-11, run_number);
	}
//	sprintf(Title_ADC_TDC_tof2_peak[i],"ADC vs. TDC (Ch. %d)  --  TOF2 Peak",i);
//	sprintf(Name_ADC_TDC_tof2_peak[i],"ADC vs. TDC (Ch. %d) - TOF2 Peak",i);
}


for (int  i=0; i<12; i++) {
h_ADC_TDC_tof1[i] = new TH2D(Name_ADC_TDC_tof1[i],Title_ADC_TDC_tof1[i],256,0,4000,256,0,4000);
//h_ADC_TDC_tof1_peak[i] = new TH2D(Name_ADC_TDC_tof1_peak[i],Title_ADC_TDC_tof1_peak[i],1000,0,4000,1000,0,4000);
}

for (int  i=0; i<24; i++) {
h_ADC_TDC_tof2[i] = new TH2D(Name_ADC_TDC_tof2[i],Title_ADC_TDC_tof2[i],256,0,4000,256,0,4000);
//h_ADC_TDC_tof2_peak[i] = new TH2D(Name_ADC_TDC_tof2_peak[i],Title_ADC_TDC_tof2_peak[i],1000,0,4000,1000,0,4000);
}


cout << "   " << endl;
cout << "File opened: " << Name_finput << endl;
cout << "  "  << endl;

/// TOF1 Print outs !
cout << "   " << endl;
for(int kk = 0; kk < 9; kk++){
	cout << "TOF1 Gap " << kk+1 << ":    " << "UP = " << par_temp_TDC_TOF1[1][kk] << "       "  << "DOWN = "<< par_temp_TDC_TOF1[1][kk+12] << endl;
}
for(int kkk = 9; kkk < 12; kkk++){
	cout << "TOF1 Gap " << kkk+1 << ":   " << "UP = " << par_temp_TDC_TOF1[1][kkk] << "       "  << "DOWN = "<< par_temp_TDC_TOF1[1][kkk+12] << endl;
}
cout << "   " << endl;

/// TOF2 Print outs !
cout << "   " << endl;
for(int jj = 0; jj < 9; jj++){
	cout << "TOF2 Gap " << jj+1 << ":    " << "A_OUT = " <<  par_temp_TDC_TOF2[1][jj] << "        " << "B_OUT = " << par_temp_TDC_TOF2[1][jj+12] << "      " << "A_IN = " << par_temp_TDC_TOF2[1][jj+24] << "       " << "B_IN = " << par_temp_TDC_TOF2[1][jj+36] << endl;
}
for(int jj = 9; jj < 12; jj++){
	cout << "TOF2 Gap " << jj+1 << ":   " << "A_OUT = " <<  par_temp_TDC_TOF2[1][jj] << "        " << "B_OUT = " << par_temp_TDC_TOF2[1][jj+12] << "      " << "A_IN = " << par_temp_TDC_TOF2[1][jj+24] << "       " << "B_IN = " << par_temp_TDC_TOF2[1][jj+36] << endl;
}

cout << "  " << endl;
cout << "  " << endl; 
//read all entries and fill the histograms
Int_t nentries = (Int_t)fChain->GetEntries();
cout <<  "Total Number of Entries :     " <<  nentries << endl;
cout << "  " << endl;

if(flag!=0) nentries=flag;

for (Int_t i=0; i<nentries; i++) {
//fChain_TOF1->GetEntry(i);
fChain->GetEntry(i);

	for (Int_t j_TOF1=0; j_TOF1<12; j_TOF1++) {
		ADC_TOF1U[j_TOF1] = ADC_tof1U[j_TOF1]-par_temp_TOF1[1][j_TOF1];
		ADC_TOF1D[j_TOF1] = ADC_tof1D[j_TOF1]-par_temp_TOF1[1][j_TOF1+12];
	}
	
	for(int i=0; i<12; i++){
		ADC_TOF1[i] = ADC_TOF1U[i];
		ADC_TOF1[i+12] = ADC_TOF1D[i];
	}


//	for (Int_t j_TOF1=0; j_TOF1<24; j_TOF1++) {
//	   ADC_TOF1[j_TOF1] = ADC_tof1[j_TOF1]-par_temp_TOF1[1][j_TOF1];
//	}

	for (Int_t j_TDCTOF=0; j_TDCTOF<12; j_TDCTOF++) {
	   TDC_tof1U[j_TDCTOF] = TDC_tof1U[j_TDCTOF] - par_temp_TDC_TOF1[1][j_TDCTOF];
	   TDC_tof1D[j_TDCTOF] = TDC_tof1D[j_TDCTOF] - par_temp_TDC_TOF1[1][j_TDCTOF+12];;
	}

	if(QUIET==0){
		if(nentries<=30000){
			if(i%1000==1) cout<<"**** "<<i<<" events done"<<endl;
		}
		if(nentries>30000){
			if(i%10000==1) cout<<"**** "<<i<<" events done"<<endl;
		}
	}

	for (int q=0; q<12; q++) {
		//if ((TDC_tof1U[q] > TDC_min_TOF1 && TDC_tof1U[q] < TDC_max_TOF1) && (TDC_tof1D[q] > TDC_min_TOF1 && TDC_tof1D[q] < TDC_max_TOF1))
			h_ADC_TDC_tof1[q]->Fill(TDC_tof1U[q],TDC_tof1D[q]);
	}

	for (Int_t j_TOF2=0; j_TOF2<12; j_TOF2++) {
		ADC_TOF2AO[j_TOF2] = ADC_tof2AO[j_TOF2]-par_temp_TOF2[1][j_TOF2];
		ADC_TOF2AI[j_TOF2] = ADC_tof2AI[j_TOF2]-par_temp_TOF2[1][j_TOF2+24];
		ADC_TOF2BO[j_TOF2] = ADC_tof2BO[j_TOF2]-par_temp_TOF2[1][j_TOF2+12];
		ADC_TOF2BI[j_TOF2] = ADC_tof2BI[j_TOF2]-par_temp_TOF2[1][j_TOF2+36];
	}

	for(int i=0; i<12; i++){
		ADC_TOF2[i] = ADC_TOF2AO[i];
		ADC_TOF2[i+12] = ADC_TOF2BO[i];
		ADC_TOF2[i+24] = ADC_TOF2AI[i];
		ADC_TOF2[i+36] = ADC_TOF2BI[i];
	}


	//for (Int_t j_TOF2=0; j_TOF2<56; j_TOF2++) {
   	//	ADC_TOF2[j_TOF2] = ADC_tof2[j_TOF2]-par_temp_TOF2[1][j_TOF2];
	//}

	for (Int_t j_TDCTOF=0; j_TDCTOF<12; j_TDCTOF++) {
	    TDC_tof2AO[j_TDCTOF] = TDC_tof2AO[j_TDCTOF] - par_temp_TDC_TOF2[1][j_TDCTOF];
	    TDC_tof2BO[j_TDCTOF] = TDC_tof2BO[j_TDCTOF] - par_temp_TDC_TOF2[1][j_TDCTOF+12];
	    TDC_tof2AI[j_TDCTOF] = TDC_tof2AI[j_TDCTOF] - par_temp_TDC_TOF2[1][j_TDCTOF+24];
	    TDC_tof2BI[j_TDCTOF] = TDC_tof2BI[j_TDCTOF] - par_temp_TDC_TOF2[1][j_TDCTOF+36];
  	}

	for (int q=0; q<12; q++) {
		//if ((TDC_tof2AO[q] > TDC_min_TOF2 && TDC_tof2AO[q] < TDC_max_TOF2) && (TDC_tof2AI[q] > TDC_min_TOF2 && TDC_tof2AI[q] < TDC_max_TOF2)) 
			h_ADC_TDC_tof2[q]->Fill(TDC_tof2AI[q],TDC_tof2AO[q]);
	}

	for (int q=0; q<12; q++) {
	//	if ((TDC_tof2BI[q] > TDC_min_TOF2 && TDC_tof2BI[q] < TDC_max_TOF2) && (TDC_tof2BO[q] > TDC_min_TOF2 && TDC_tof2BO[q] < TDC_max_TOF2)) 
			h_ADC_TDC_tof2[q+12]->Fill(TDC_tof2BI[q],TDC_tof2BO[q]);
	}		

}


////Draw TOF1 Histograms

char Name_Can_ADC_tof1B[100];			char Title_Can_ADC_tof1B[100];

sprintf(Name_Can_ADC_tof1B,"TOF1 TDC Up vs. Down -- Run %d  (Ch. 0 - 23) -- Mapping Sorted",run_number);
sprintf(Title_Can_ADC_tof1B,"TOF1 TDC Up vs. Down -- Run %d  (Ch. 0 - 23) -- Mapping Sorted",run_number);

TCanvas *c7;
	c7 = new TCanvas(Name_Can_ADC_tof1B,Title_Can_ADC_tof1B,1000,400); 
	c7->Divide(4,3);
	c7->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=0; ican<12; ican++){
		c7->cd(ican+1);
		h_ADC_TDC_tof1[ican]->SetAxisRange(TOF1_xmin, TOF1_xmax, "X");
		h_ADC_TDC_tof1[ican]->SetAxisRange(TOF1_ymin, TOF1_ymax, "Y");
		//h_ADC_TDC_tof1_peak[ican]->SetMarkerColor(3);
		//h_ADC_TDC_tof1_peak[ican]->SetMarkerSize(50);
		//h_ADC_TDC_tof1_peak[ican]->Draw();
		h_ADC_TDC_tof1[ican]->Draw("colz");
	}

////Draw TOF2 Histograms

char Name_Can_ADC_tof2B[100];			char Title_Can_ADC_tof2B[100];

sprintf(Name_Can_ADC_tof2B,"TOF2 TDC In vs. Out -- Run %d  (Ch. 0 - 47) -- Mapping Sorted",run_number);
sprintf(Title_Can_ADC_tof2B,"TOF2 TDC In vs. Out -- Run %d  (Ch. 0 - 47) -- Mapping Sorted",run_number);

TCanvas *c8;
	c8 = new TCanvas(Name_Can_ADC_tof2B,Title_Can_ADC_tof2B,1200,500); 
	c8->Divide(6,4);
	c8->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=0; ican<12; ican++){
		c8->cd(ican+1);
		h_ADC_TDC_tof2[ican]->SetAxisRange(TOF2_xmin, TOF2_xmax, "X");
		h_ADC_TDC_tof2[ican]->SetAxisRange(TOF2_ymin, TOF2_ymax, "Y");
		//h_ADC_TDC_tof2_peak[ican]->SetMarkerColor(3);
		//h_ADC_TDC_tof2_peak[ican]->SetMarkerSize(50);
		//h_ADC_TDC_tof2_peak[ican]->Draw();
		h_ADC_TDC_tof2[ican]->Draw("colz");
	}

	for(Int_t ican=12; ican<24; ican++){
		c8->cd(ican+1);
		h_ADC_TDC_tof2[ican]->SetAxisRange(TOF2_xmin, TOF2_xmax, "X");
		h_ADC_TDC_tof2[ican]->SetAxisRange(TOF2_ymin, TOF2_ymax, "Y");
		//h_ADC_TDC_tof2_peak[ican+24]->SetMarkerColor(3);
		//h_ADC_TDC_tof2_peak[ican+24]->SetMarkerSize(50);
		//h_ADC_TDC_tof2_peak[ican+24]->Draw();
		h_ADC_TDC_tof2[ican]->Draw("colz");
	}

}

