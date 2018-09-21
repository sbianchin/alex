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
#include "ANAPATH.h"
#include "Thresholds.h"
#include "ADC_Thresholds.h"
#include "Cuts_and_Windows.h"
//#include "Event_Display_MS.h"
#endif

void SFT_Fiber_Efficiency(Int_t run_number=1181, Int_t flag=0) {

gStyle->SetOptStat(111111111);

Int_t count_l1or = 0;
Int_t count_l1and = 0;
Int_t count_l2or = 0;
Int_t count_l2and = 0;
Int_t count_l3or = 0;
Int_t count_l3and = 0;
Int_t count_l4or = 0;
Int_t count_l4and = 0;

Int_t countlayer1or[15] = {0};
Int_t countlayer2or[15] = {0};
Int_t countlayer3or[17] = {0};
Int_t countlayer4or[17] = {0};

Int_t countlayer1and[15] = {0};
Int_t countlayer2and[15] = {0};
Int_t countlayer3and[17] = {0};
Int_t countlayer4and[17] = {0};

Int_t ADC_High_target[256];		//Int_t ADC_High_sft[128];
Int_t ADC_Low_target[256];		Int_t ADC_Low_sft[128];
Int_t TDC_LE_target[256][16];		Int_t TDC_LE_sft[128][16];
Int_t TDC_TE_target[256][16];		Int_t TDC_TE_sft[128][16];

Int_t  ADC_High_sft[128];

Int_t ADC_tof1[24];
Int_t ADC_tof2[48];

Int_t fiber[128]={-1};
//double fiber[128]={-1.};
Int_t fiberTDC[128]={-1};

char Name_Can[100];			char Title_Can[100];
sprintf(Name_Can,"ADC -- Run %d -- SFT_Fiber_Efficiency.C",run_number);
sprintf(Title_Can,"ADC -- Run %d  -- SFT_Fiber_Efficiency.C",run_number);

char par_map[200];		char par_finput[200];
//sprintf(par_finput,"/media/bianchin/hdd1/trek/E36/Macros/April_2015/Mapping/SFT_Mapping_Oct14.txt");
sprintf(par_map,"../Mapping");
sprintf(par_finput,"%s/SFT_Mapping_Oct14.txt",par_map);

char path_input[200];
//sprintf(path_input,"/media/bianchin/hdd1/trek/E36/Data/June_2015/root/Merged");
sprintf(path_input,"%s",path_merged);

char Name_finput[200];
sprintf(Name_finput,"%s/Run%dMS.root",path_input, run_number);

///

char par_finput2[200];
//sprintf(par_finput,"/triumfcs/trshare/trek/E36/Macros/April_2015/Mapping/%s",source_mapping);
sprintf(par_finput2,"%s/MWPC_map.txt",par_map);

char par_finput3[200];
//sprintf(par_finput,"/triumfcs/trshare/trek/E36/Macros/April_2015/Mapping/%s",source_mapping);
sprintf(par_finput3,"%s/ADC_TARGET_Thresholds.txt",par_map);

char par_finput4[200];
//sprintf(par_finput,"/triumfcs/trshare/trek/E36/Macros/April_2015/Mapping/%s",source_mapping);
sprintf(par_finput4,"%s/ADC_SFT_Thresholds.txt",par_map);

///

Int_t par_temp[2][128];
ifstream fdat(par_finput,ios::in);
for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
fdat.close();

Int_t par_temp_TARGET[2][256];
ifstream fdat3(par_finput3,ios::in);
for(Int_t ii=0; ii<256; ii++) fdat3 >> par_temp_TARGET[0][ii] >> par_temp_TARGET[1][ii];
fdat3.close();

Int_t par_temp_SFT[2][128];
ifstream fdat4(par_finput4,ios::in);
for(Int_t ii=0; ii<128; ii++) fdat4 >> par_temp_SFT[0][ii] >> par_temp_SFT[1][ii];
fdat4.close();

/////////////////////// Histograms
TH1D *h_layer1or;	char Title_h_layer1or[100];	char Name_h_layer1or[100];
TH1D *h_layer1and;	char Title_h_layer1and[100];	char Name_h_layer1and[100];
TH1D *h_layer2or;	char Title_h_layer2or[100];	char Name_h_layer2or[100];
TH1D *h_layer2and;	char Title_h_layer2and[100];	char Name_h_layer2and[100];
TH1D *h_layer3or;	char Title_h_layer3or[100];	char Name_h_layer3or[100];
TH1D *h_layer3and;	char Title_h_layer3and[100];	char Name_h_layer3and[100];
TH1D *h_layer4or;	char Title_h_layer4or[100];	char Name_h_layer4or[100];
TH1D *h_layer4and;	char Title_h_layer4and[100];	char Name_h_layer4and[100];

sprintf(Title_h_layer1or,"SFT Fiber Efficiencies Layer 1 - Run %d: ADC > %d , %d < TDC < %d", run_number, SFT_ADC_Thr_HG_Offset, TDC_min_SFT, TDC_max_SFT); 
sprintf(Name_h_layer1or,"Layer 1"); 
h_layer1or = new TH1D(Name_h_layer1or,Title_h_layer1or,17,0,17);

sprintf(Title_h_layer1and,"SFT Fiber Efficiencies Layer 1 - Run %d: ADC > %d , %d < TDC < %d", run_number, SFT_ADC_Thr_HG_Offset, TDC_min_SFT, TDC_max_SFT); 
sprintf(Name_h_layer1and,"Layer 1"); 
h_layer1and = new TH1D(Name_h_layer1and,Title_h_layer1and,17,0,17);

sprintf(Title_h_layer2or,"SFT Fiber Efficiencies Layer 2 - Run %d: ADC > %d , %d < TDC < %d", run_number, SFT_ADC_Thr_HG_Offset, TDC_min_SFT, TDC_max_SFT); 
sprintf(Name_h_layer2or,"Layer 2"); 
h_layer2or = new TH1D(Name_h_layer2or,Title_h_layer2or,17,0,17);

sprintf(Title_h_layer2and,"SFT Fiber Efficiencies Layer 2 - Run %d: ADC > %d , %d < TDC < %d", run_number, SFT_ADC_Thr_HG_Offset, TDC_min_SFT, TDC_max_SFT); 
sprintf(Name_h_layer2and,"Layer 2"); 
h_layer2and = new TH1D(Name_h_layer2and,Title_h_layer2and,17,0,17);

sprintf(Title_h_layer3or,"SFT Fiber Efficiencies Layer 3 - Run %d: ADC > %d , %d < TDC < %d", run_number, SFT_ADC_Thr_HG_Offset, TDC_min_SFT, TDC_max_SFT); 
sprintf(Name_h_layer3or,"Layer 3"); 
h_layer3or = new TH1D(Name_h_layer3or,Title_h_layer3or,19,0,19);

sprintf(Title_h_layer3and,"SFT Fiber Efficiencies Layer 3 - Run %d: ADC > %d , %d < TDC < %d", run_number, SFT_ADC_Thr_HG_Offset, TDC_min_SFT, TDC_max_SFT); 
sprintf(Name_h_layer3and,"Layer 3"); 
h_layer3and = new TH1D(Name_h_layer3and,Title_h_layer3and,19,0,19);

sprintf(Title_h_layer4or,"SFT Fiber Efficiencies Layer 4 - Run %d: ADC > %d , %d < TDC < %d", run_number, SFT_ADC_Thr_HG_Offset, TDC_min_SFT, TDC_max_SFT); 
sprintf(Name_h_layer4or,"Layer 4"); 
h_layer4or = new TH1D(Name_h_layer4or,Title_h_layer4or,19,0,19);

sprintf(Title_h_layer4and,"SFT Fiber Efficiencies Layer 4 - Run %d: ADC > %d , %d < TDC < %d", run_number, SFT_ADC_Thr_HG_Offset, TDC_min_SFT, TDC_max_SFT); 
sprintf(Name_h_layer4and,"Layer 4"); 
h_layer4and = new TH1D(Name_h_layer4and,Title_h_layer4and,19,0,19);

///////////////////////

//TChain *fChain_TARGET = new TChain("Tree_TARGET");
//TChain *fChain_SFT = new TChain("Tree_SFT");

//fChain_TARGET->Add(Name_finput);      
//fChain_SFT->Add(Name_finput);      

cout << "   " << endl;
cout << "Mapping File:   " << par_finput << endl;
cout << "" << endl;

//cout << "    " << endl;

//cout << "I'm in the File !" << endl;
//cout << Name_finput << endl;
//cout << "    " << endl;

TChain *fChain = new TChain("Tree");

fChain->Add(Name_finput);      

cout << "    " << endl;
cout << "I'm in the File !" << endl;
cout << Name_finput << endl;
cout << "    " << endl;

fChain->SetMakeClass(1);
fChain->SetBranchAddress("ADC_High_TARGET",ADC_High_target);
fChain->SetBranchAddress("ADC_Low_TARGET",ADC_Low_target);
fChain->SetBranchAddress("TDC_LE_TARGET",TDC_LE_target);
fChain->SetBranchAddress("TDC_TE_TARGET",TDC_TE_target);

fChain->SetBranchAddress("ADC_High_SFT",ADC_High_sft);
fChain->SetBranchAddress("ADC_Low_SFT",ADC_Low_sft);
fChain->SetBranchAddress("TDC_LE_SFT",TDC_LE_sft);
fChain->SetBranchAddress("TDC_TE_SFT",TDC_TE_sft);

//fChain->SetBranchAddress("ADC_TOF1",ADC_tof1);

//fChain->SetBranchAddress("ADC_TOF2",ADC_tof2);

//fChain_TARGET->SetMakeClass(1);
//fChain_TARGET->SetBranchAddress("ADC_High_TARGET",ADC_High_target);
//fChain_TARGET->SetBranchAddress("ADC_Low_TARGET",ADC_Low_target);
//fChain_TARGET->SetBranchAddress("TDC_LE_TARGET",TDC_LE_target);
//fChain_TARGET->SetBranchAddress("TDC_TE_TARGET",TDC_TE_target);

//fChain_SFT->SetMakeClass(1);
//fChain_SFT->SetBranchAddress("ADC_High_SFT",ADC_High_sft);
//fChain_SFT->SetBranchAddress("ADC_Low_SFT",ADC_Low_sft);
//fChain_SFT->SetBranchAddress("TDC_LE_SFT",TDC_LE_sft);
//fChain_SFT->SetBranchAddress("TDC_TE_SFT",TDC_TE_sft);


Int_t HG_SFT_ADC_Thr[128] = {0};

bool l1_or = false;
bool l1_and = false;
bool l2_or = false;
bool l2_and = false;
bool l3_or = false;
bool l3_and = false;
bool l4_or = false;
bool l4_and = false;



Int_t nentries = (Int_t)fChain->GetEntries();
cout <<  "DEBUG :     " <<  nentries << endl;
 
for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;

if(flag!=0) nentries=flag;
for (Int_t j=0; j<nentries+1; j++) {
//for (Int_t j=0; j<1; j++) {

//	fChain_SFT->GetEntry(j);
	fChain->GetEntry(j);	

	if(j%10000==1)
	cout<<"**** "<< j <<" events done"<<endl;

	l1_or = false;
	l1_and = false;
	l2_or = false;
	l2_and = false;
	l3_or = false;
	l3_and = false;
	l4_or = false;
	l4_and = false;

	// Variable Change
	for(Int_t ii=0; ii<128; ii++){

		//cout << "DEBUG:  " << ii << "  " << ADC_High_sft[ii] << endl;

		//ADC_High_sft[ii] = ADC_High_sft[ii] - par_temp_SFT[1][ii];
		//cout << "BEFORE  " << ADC_High_sft[ii] << "  " << HG_SFT_ADC_Thr[ii]<< endl;
		ADC_High_sft[ii] = ADC_High_sft[ii] - HG_SFT_ADC_Thr[ii];
		fiber[ii]=ADC_High_sft[par_temp[1][ii]];
		fiberTDC[ii]=TDC_LE_sft[par_temp[1][ii]][0];
		//cout << "AFTER  " << ADC_High_sft[ii] << "  " <<  SFT_ADC_Thr_HG[ii] << endl;
	}


///Layer 1
	for (Int_t i=0; i<15; i++) {
		if (((fiber[i]>0) && (fiberTDC[i]>TDC_min_SFT) && (fiberTDC[i]<TDC_max_SFT)) || ((fiberTDC[i+64]>TDC_min_SFT) && (fiberTDC[i+64]<TDC_max_SFT) && (fiber[i+64]>0))) {
			countlayer1or[i]++;
			l1_or = true;
		}
		if (((fiber[i]>0) && (fiberTDC[i]>TDC_min_SFT) && (fiberTDC[i]<TDC_max_SFT)) && ((fiberTDC[i+64]>TDC_min_SFT) && (fiberTDC[i+64]<TDC_max_SFT) && (fiber[i+64]>0))) {
			countlayer1and[i]++;
			l1_and = true;
		}
	}
///Layer 2
	for (Int_t i=15; i<30; i++) {
		if (((fiber[i]>0) && (fiberTDC[i]>TDC_min_SFT) && (fiberTDC[i]<TDC_max_SFT)) || ((fiberTDC[i+64]>TDC_min_SFT) && (fiberTDC[i+64]<TDC_max_SFT) && (fiber[i+64]>0))) {
			countlayer2or[i-15]++;
			l2_or = true;
		}
		if (((fiber[i]>0) && (fiberTDC[i]>TDC_min_SFT) && (fiberTDC[i]<TDC_max_SFT)) && ((fiberTDC[i+64]>TDC_min_SFT) && (fiberTDC[i+64]<TDC_max_SFT) && (fiber[i+64]>0))) {
			countlayer2and[i-15]++;
			l2_and = true;
		}
	}
///Layer 3
	for (Int_t i=30; i<47; i++) {
		if (((fiber[i]>0) && (fiberTDC[i]>TDC_min_SFT) && (fiberTDC[i]<TDC_max_SFT)) || ((fiberTDC[i+64]>TDC_min_SFT) && (fiberTDC[i+64]<TDC_max_SFT) && (fiber[i+64]>0))) {
			countlayer3or[i-30]++;
			l3_or = true;
		}
		if (((fiber[i]>0) && (fiberTDC[i]>TDC_min_SFT) && (fiberTDC[i]<TDC_max_SFT)) && ((fiberTDC[i+64]>TDC_min_SFT) && (fiberTDC[i+64]<TDC_max_SFT) && (fiber[i+64]>0))) {
			countlayer3and[i-30]++;
			l3_and = true;
		}
	}
///Layer 4
	for (Int_t i=47; i<64; i++) {
		if (((fiber[i]>0) && (fiberTDC[i]>TDC_min_SFT) && (fiberTDC[i]<TDC_max_SFT)) || ((fiberTDC[i+64]>TDC_min_SFT) && (fiberTDC[i+64]<TDC_max_SFT) && (fiber[i+64]>0))) {
			countlayer4or[i-47]++;
			l4_or = true;
		}
		if (((fiber[i]>0) && (fiberTDC[i]>TDC_min_SFT) && (fiberTDC[i]<TDC_max_SFT)) && ((fiberTDC[i+64]>TDC_min_SFT) && (fiberTDC[i+64]<TDC_max_SFT) && (fiber[i+64]>0))) {
			countlayer4and[i-47]++;
			l4_and = true;
		}
	}

// Individual Layer Counter	
	if (l1_or) count_l1or++; 
	if (l1_and) count_l1and++; 
	if (l2_or) count_l2or++;
	if (l2_and) count_l2and++; 	
	if (l3_or) count_l3or++; 	
	if (l3_and) count_l3and++; 	
	if (l4_or) count_l4or++; 	
	if (l4_and) count_l4and++; 	
}

//cout << count_l1or++ << endl; 
//cout << count_l1and++ << endl; 
//cout << count_l2or++ << endl; 
//cout << count_l2and++ << endl; 	
//cout << count_l3or++ << endl; 	
//cout << count_l3and++ << endl; 	
//cout << count_l4or++ << endl; 	
//cout << count_l4and++ << endl; 


/*
for(Int_t i=0; i<15; i++) {
	cout << "Layer 1 Fiber " << i << " count: " << countlayer1or[i] << endl;
}

cout << " " << endl;

for(Int_t i=0; i<15; i++) {
	cout << "Layer 2 Fiber " << i << " count: " << countlayer2or[i] << endl;
}

cout << " " << endl;

for(Int_t i=0; i<17; i++) {
	cout << "Layer 3 Fiber " << i << " count: " << countlayer3or[i] << endl;
}

cout << " " << endl;

for(Int_t i=0; i<17; i++) {
	cout << "Layer 4 Fiber " << i << " count: " << countlayer4or[i] << endl;
}	

for(Int_t i=0; i<15; i++) h_layer1or->Fill(i+1,countlayer1or[i]);
for(Int_t i=0; i<15; i++) h_layer1and->Fill(i+1,countlayer1and[i]);
for(Int_t i=0; i<15; i++) h_layer2or->Fill(i+1,countlayer2or[i]);
for(Int_t i=0; i<15; i++) h_layer2and->Fill(i+1,countlayer2and[i]);

for(Int_t j=0; j<17; j++) h_layer3or->Fill(j+1,countlayer3or[j]);
for(Int_t j=0; j<17; j++) h_layer3and->Fill(j+1,countlayer3and[j]);
for(Int_t j=0; j<17; j++) h_layer4or->Fill(j+1,countlayer4or[j]);
for(Int_t j=0; j<17; j++) h_layer4and->Fill(j+1,countlayer4and[j]);
*/

cout << "TOTAL  :  " << count_l4or << endl;
for(int i=0; i<17; i++) cout << i << "    " << countlayer4or[i] << endl;

for(Int_t i=0; i<15; i++) h_layer1or->Fill(i+1,double(countlayer1or[i])/double(count_l1or));
for(Int_t i=0; i<15; i++) h_layer1and->Fill(i+1,double(countlayer1and[i])/double(count_l1or));
for(Int_t i=0; i<15; i++) h_layer2or->Fill(i+1,double(countlayer2or[i])/double(count_l2or));
for(Int_t i=0; i<15; i++) h_layer2and->Fill(i+1,double(countlayer2and[i])/double(count_l2or));

for(Int_t j=0; j<17; j++) h_layer3or->Fill(j+1,double(countlayer3or[j])/double(count_l3or));
for(Int_t j=0; j<17; j++) h_layer3and->Fill(j+1,double(countlayer3and[j])/double(count_l3or));
for(Int_t j=0; j<17; j++) h_layer4or->Fill(j+1,double(countlayer4or[j])/double(count_l4or));
for(Int_t j=0; j<17; j++) h_layer4and->Fill(j+1,double(countlayer4and[j])/double(count_l4or));

TCanvas *c1;
c1 = new TCanvas(Name_Can,Title_Can,1000,900); 
c1->Divide(2,2);
c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

TLegend *l1;
l1 = new TLegend(0.82,0.87,0.99,0.99);
l1->AddEntry(h_layer1or, "Up OR Down", "F");
l1->AddEntry(h_layer1and, "Up AND Down", "F");

c1->cd(1);
h_layer1or->SetLineColor(4);
h_layer1or->SetFillStyle(3004);
h_layer1or->SetFillColor(4);
h_layer1or->Draw();

h_layer1and->SetLineColor(2);
h_layer1and->SetFillStyle(3005);
h_layer1and->SetFillColor(2);
h_layer1and->Draw("same");
l1->Draw("same");

c1->cd(2);
h_layer2or->SetLineColor(4);
h_layer2or->SetFillStyle(3004);
h_layer2or->SetFillColor(4);
h_layer2or->Draw();

h_layer2and->SetLineColor(2);
h_layer2and->SetFillStyle(3005);
h_layer2and->SetFillColor(2);
h_layer2and->Draw("same");
l1->Draw("same");

c1->cd(3);

h_layer3or->SetLineColor(4);
h_layer3or->SetFillStyle(3004);
h_layer3or->SetFillColor(4);
h_layer3or->Draw();

h_layer3and->SetLineColor(2);
h_layer3and->SetFillStyle(3005);
h_layer3and->SetFillColor(2);
h_layer3and->Draw("same");
l1->Draw("same");

c1->cd(4);

h_layer4or->SetLineColor(4);
h_layer4or->SetFillStyle(3004);
h_layer4or->SetFillColor(4);
h_layer4or->Draw();

h_layer4and->SetLineColor(2);
h_layer4and->SetFillStyle(3005);
h_layer4and->SetFillColor(2);
h_layer4and->Draw("same");
l1->Draw("same");


}
