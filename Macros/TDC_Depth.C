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
#include "TTree.h"
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
#include "ANAPATH.h"
#include "Thresholds.h" 
//#include "TSpectrum.h"
#endif

void TDC_Depth(Int_t run_number=1181, Int_t flag=0) {

gStyle->SetOptStat(111111111);

Int_t counter_TDCA[66][4]= {{0}};
Int_t counter_TDCB[66][4]= {{0}};

Int_t counter_TARGET_TDC_A[66][4]= {{0}};
Int_t counter_TARGET_TDC_B[66][4]= {{0}};
Int_t counter_TARGET_TDC_C[66][4]= {{0}};
Int_t counter_TARGET_TDC_D[66][4]= {{0}};

////

Int_t ADC_High_target[256];		Int_t ADC_High_sft[128];
Int_t ADC_Low_target[256];		Int_t ADC_Low_sft[128];
Int_t TDC_LE_target[256][16];		Int_t TDC_LE_sft[128][16];
Int_t TDC_TE_target[256][16];		Int_t TDC_TE_sft[128][16];

//Int_t fiber[128]={-1};

char Name_Can[100];			char Title_Can[100];
sprintf(Name_Can,"TDC -- Run %d -- TDC_LE_SFT_Multiplicity.C",run_number);
sprintf(Title_Can,"TDC -- Run %d  -- TDC_LE_SFT_Multiplicity.C",run_number);

//char path_input[100];
//sprintf(path_input,"/triumfcs/trshare/trek/E36/Data/April_2015/root");

char path_input[100];
sprintf(path_input,"%s",path_merged);

char Name_finput[200];
sprintf(Name_finput,"%s/Run%dMS.root",path_input, run_number);

//// SFT

TH1D *h_TDC_0A;		char Title_h_TDC_0A[100];	char Name_h_TDC_0A[100];
TH1D *h_TDC_0B;		char Title_h_TDC_0B[100];	char Name_h_TDC_0B[100];
TH1D *h_TDC_1A;		char Title_h_TDC_1A[100];	char Name_h_TDC_1A[100];
TH1D *h_TDC_1B;		char Title_h_TDC_1B[100];	char Name_h_TDC_1B[100];

TH1D *h_TDC_2A;		char Title_h_TDC_2A[100];	char Name_h_TDC_2A[100];
TH1D *h_TDC_2B;		char Title_h_TDC_2B[100];	char Name_h_TDC_2B[100];
TH1D *h_TDC_3A;		char Title_h_TDC_3A[100];	char Name_h_TDC_3A[100];
TH1D *h_TDC_3B;		char Title_h_TDC_3B[100];	char Name_h_TDC_3B[100];

//// TARGET

TH1D *h_TDC_TARGET_0A;		char Title_h_TDC_TARGET_0A[100];	char Name_h_TDC_TARGET_0A[100];
TH1D *h_TDC_TARGET_0B;		char Title_h_TDC_TARGET_0B[100];	char Name_h_TDC_TARGET_0B[100];
TH1D *h_TDC_TARGET_1A;		char Title_h_TDC_TARGET_1A[100];	char Name_h_TDC_TARGET_1A[100];
TH1D *h_TDC_TARGET_1B;		char Title_h_TDC_TARGET_1B[100];	char Name_h_TDC_TARGET_1B[100];

TH1D *h_TDC_TARGET_2A;		char Title_h_TDC_TARGET_2A[100];	char Name_h_TDC_TARGET_2A[100];
TH1D *h_TDC_TARGET_2B;		char Title_h_TDC_TARGET_2B[100];	char Name_h_TDC_TARGET_2B[100];
TH1D *h_TDC_TARGET_3A;		char Title_h_TDC_TARGET_3A[100];	char Name_h_TDC_TARGET_3A[100];
TH1D *h_TDC_TARGET_3B;		char Title_h_TDC_TARGET_3B[100];	char Name_h_TDC_TARGET_3B[100];

TH1D *h_TDC_TARGET_0C;		char Title_h_TDC_TARGET_0C[100];	char Name_h_TDC_TARGET_0C[100];
TH1D *h_TDC_TARGET_0D;		char Title_h_TDC_TARGET_0D[100];	char Name_h_TDC_TARGET_0D[100];
TH1D *h_TDC_TARGET_1C;		char Title_h_TDC_TARGET_1C[100];	char Name_h_TDC_TARGET_1C[100];
TH1D *h_TDC_TARGET_1D;		char Title_h_TDC_TARGET_1D[100];	char Name_h_TDC_TARGET_1D[100];

TH1D *h_TDC_TARGET_2C;		char Title_h_TDC_TARGET_2C[100];	char Name_h_TDC_TARGET_2C[100];
TH1D *h_TDC_TARGET_2D;		char Title_h_TDC_TARGET_2D[100];	char Name_h_TDC_TARGET_2D[100];
TH1D *h_TDC_TARGET_3C;		char Title_h_TDC_TARGET_3C[100];	char Name_h_TDC_TARGET_3C[100];
TH1D *h_TDC_TARGET_3D;		char Title_h_TDC_TARGET_3D[100];	char Name_h_TDC_TARGET_3D[100];

////////////////////////////////////////////

sprintf(Title_h_TDC_0A,"SFT | %d < TDC < %d : Run %d - Channels 0-63", TDC_min_SFT, TDC_max_SFT ,run_number); 
sprintf(Name_h_TDC_0A,"TDCA[0]"); 
h_TDC_0A = new TH1D(Name_h_TDC_0A,Title_h_TDC_0A,66,-1,65);

sprintf(Title_h_TDC_0B,"SFT | %d < TDC < %d : Run %d - Channels 64-127", TDC_min_SFT, TDC_max_SFT ,run_number); 
sprintf(Name_h_TDC_0B,"TDCB[0]"); 
h_TDC_0B = new TH1D(Name_h_TDC_0B,Title_h_TDC_0B,66,63,129);

sprintf(Title_h_TDC_1A,"TDC_LE_SFT Fiber Multiplicities: Run %d - TDC[1]", run_number); 
sprintf(Name_h_TDC_1A,"TDCA[1]"); 
h_TDC_1A = new TH1D(Name_h_TDC_1A,Title_h_TDC_1A,66,-1,65);

sprintf(Title_h_TDC_1B,"TDC_LE_SFT Fiber Multiplicities: Run %d - TDC[1]", run_number); 
sprintf(Name_h_TDC_1B,"TDCB[1]"); 
h_TDC_1B = new TH1D(Name_h_TDC_1B,Title_h_TDC_1B,66,63,129);

sprintf(Title_h_TDC_2A,"TDC_LE_SFT Fiber Multiplicities: Run %d - TDC[2]", run_number); 
sprintf(Name_h_TDC_2A,"TDCA[2]"); 
h_TDC_2A = new TH1D(Name_h_TDC_2A,Title_h_TDC_2A,66,-1,65);

sprintf(Title_h_TDC_2B,"TDC_LE_SFT Fiber Multiplicities: Run %d - TDC[2]", run_number); 
sprintf(Name_h_TDC_2B,"TDCB[2]"); 
h_TDC_2B = new TH1D(Name_h_TDC_2B,Title_h_TDC_2B,66,63,129);

sprintf(Title_h_TDC_3A,"TDC_LE_SFT Fiber Multiplicities: Run %d - TDC[3]", run_number); 
sprintf(Name_h_TDC_3A,"TDCA[3]"); 
h_TDC_3A = new TH1D(Name_h_TDC_3A,Title_h_TDC_3A,66,-1,65);

sprintf(Title_h_TDC_3B,"TDC_LE_SFT Fiber Multiplicities: Run %d - TDC[3]", run_number); 
sprintf(Name_h_TDC_3B,"TDCB[3]"); 
h_TDC_3B = new TH1D(Name_h_TDC_3B,Title_h_TDC_3B,66,63,129);

////////////////////////////////////////////////////////////////////////////////////

sprintf(Title_h_TDC_TARGET_0A,"TARGET | %d < TDC < %d : Run %d - Channels 0-63", TDC_min_TARGET, TDC_max_TARGET ,run_number); 
sprintf(Name_h_TDC_TARGET_0A,"TDCTA[0]"); 
h_TDC_TARGET_0A = new TH1D(Name_h_TDC_TARGET_0A,Title_h_TDC_TARGET_0A,66,-1,65);

sprintf(Title_h_TDC_TARGET_0B,"TARGET | %d < TDC < %d : Run %d - Channels 64-127", TDC_min_TARGET, TDC_max_TARGET ,run_number); 
sprintf(Name_h_TDC_TARGET_0B,"TDCTB[0]"); 
h_TDC_TARGET_0B = new TH1D(Name_h_TDC_TARGET_0B,Title_h_TDC_TARGET_0B,66,63,129);

sprintf(Title_h_TDC_TARGET_0C,"TARGET | %d < TDC < %d : Run %d - Channels 128-191", TDC_min_TARGET, TDC_max_TARGET ,run_number); 
sprintf(Name_h_TDC_TARGET_0C,"TDCTC[0]"); 
h_TDC_TARGET_0C = new TH1D(Name_h_TDC_TARGET_0C,Title_h_TDC_TARGET_0C,66,127,193);

sprintf(Title_h_TDC_TARGET_0D,"TARGET | %d < TDC < %d : Run %d - Channels 192-255", TDC_min_TARGET, TDC_max_TARGET ,run_number); 
sprintf(Name_h_TDC_TARGET_0D,"TDCTD[0]"); 
h_TDC_TARGET_0D = new TH1D(Name_h_TDC_TARGET_0D,Title_h_TDC_TARGET_0D,66,191,257);

sprintf(Title_h_TDC_TARGET_1A,"TDC_LE_SFT Fiber Multiplicities: Run %d - TDC[1]", run_number); 
sprintf(Name_h_TDC_TARGET_1A,"TDCTA[1]"); 
h_TDC_TARGET_1A = new TH1D(Name_h_TDC_TARGET_1A,Title_h_TDC_TARGET_1A,66,-1,65);

sprintf(Title_h_TDC_TARGET_1B,"TDC_LE_SFT Fiber Multiplicities: Run %d - TDC[1]", run_number); 
sprintf(Name_h_TDC_TARGET_1B,"TDCTB[1]"); 
h_TDC_TARGET_1B = new TH1D(Name_h_TDC_TARGET_1B,Title_h_TDC_TARGET_1B,66,63,129);

sprintf(Title_h_TDC_TARGET_1C,"TDC_LE_SFT Fiber Multiplicities: Run %d - TDC[1]", run_number); 
sprintf(Name_h_TDC_TARGET_1C,"TDCTC[1]"); 
h_TDC_TARGET_1C = new TH1D(Name_h_TDC_TARGET_1C,Title_h_TDC_TARGET_1C,66,127,193);

sprintf(Title_h_TDC_TARGET_1D,"TDC_LE_SFT Fiber Multiplicities: Run %d - TDC[1]", run_number); 
sprintf(Name_h_TDC_TARGET_1D,"TDCTD[1]"); 
h_TDC_TARGET_1D = new TH1D(Name_h_TDC_TARGET_1D,Title_h_TDC_TARGET_1D,66,191,257);

sprintf(Title_h_TDC_TARGET_2A,"TDC_LE_SFT Fiber Multiplicities: Run %d - TDC[2]", run_number); 
sprintf(Name_h_TDC_TARGET_2A,"TDCTA[2]"); 
h_TDC_TARGET_2A = new TH1D(Name_h_TDC_TARGET_2A,Title_h_TDC_TARGET_2A,66,-1,65);

sprintf(Title_h_TDC_TARGET_2B,"TDC_LE_SFT Fiber Multiplicities: Run %d - TDC[2]", run_number); 
sprintf(Name_h_TDC_TARGET_2B,"TDCTB[2]"); 
h_TDC_TARGET_2B = new TH1D(Name_h_TDC_TARGET_2B,Title_h_TDC_TARGET_2B,66,63,129);

sprintf(Title_h_TDC_TARGET_2C,"TDC_LE_SFT Fiber Multiplicities: Run %d - TDC[2]", run_number); 
sprintf(Name_h_TDC_TARGET_2C,"TDCTC[2]"); 
h_TDC_TARGET_2C = new TH1D(Name_h_TDC_TARGET_2C,Title_h_TDC_TARGET_2C,66,127,193);

sprintf(Title_h_TDC_TARGET_2D,"TDC_LE_SFT Fiber Multiplicities: Run %d - TDC[2]", run_number); 
sprintf(Name_h_TDC_TARGET_2D,"TDCTD[2]"); 
h_TDC_TARGET_2D = new TH1D(Name_h_TDC_TARGET_2D,Title_h_TDC_TARGET_2D,66,191,257);

sprintf(Title_h_TDC_TARGET_3A,"TDC_LE_SFT Fiber Multiplicities: Run %d - TDC[3]", run_number); 
sprintf(Name_h_TDC_TARGET_3A,"TDCTA[3]"); 
h_TDC_TARGET_3A = new TH1D(Name_h_TDC_TARGET_3A,Title_h_TDC_TARGET_3A,66,-1,65);

sprintf(Title_h_TDC_TARGET_3B,"TDC_LE_SFT Fiber Multiplicities: Run %d - TDC[3]", run_number); 
sprintf(Name_h_TDC_TARGET_3B,"TDCTB[3]"); 
h_TDC_TARGET_3B = new TH1D(Name_h_TDC_TARGET_3B,Title_h_TDC_TARGET_3B,66,63,129);

sprintf(Title_h_TDC_TARGET_3C,"TDC_LE_SFT Fiber Multiplicities: Run %d - TDC[3]", run_number); 
sprintf(Name_h_TDC_TARGET_3C,"TDCTC[3]"); 
h_TDC_TARGET_3C = new TH1D(Name_h_TDC_TARGET_3C,Title_h_TDC_TARGET_3C,66,127,193);

sprintf(Title_h_TDC_TARGET_3D,"TDC_LE_SFT Fiber Multiplicities: Run %d - TDC[3]", run_number); 
sprintf(Name_h_TDC_TARGET_3D,"TDCTD[3]"); 
h_TDC_TARGET_3D = new TH1D(Name_h_TDC_TARGET_3D,Title_h_TDC_TARGET_3D,66,191,257);

////////////////////////////////////////////////////////////////////////////////////

//TChain *fChain_TARGET = new TChain("Tree_TARGET");
//TChain *fChain_SFT = new TChain("Tree_SFT");

//fChain_TARGET->Add(Name_finput);      
//fChain_SFT->Add(Name_finput);

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

Int_t nentries = (Int_t)fChain->GetEntries();
cout <<  "DEBUG :     " <<  nentries << endl;

if(flag!=0) nentries=flag;
for (Int_t j=0; j<nentries; j++) {
	
	if(j%10000==1)
	cout<<"**** "<< j <<" events done"<<endl;
//	fChain_SFT->GetEntry(j);
	fChain->GetEntry(j);
	
	for (Int_t i=0; i<64; i++) {
		for (Int_t m=0; m<4; m++) {
			if ((TDC_LE_sft[i][m] > TDC_min_SFT) && (TDC_LE_sft[i][m] < TDC_max_SFT)) {
				counter_TDCA[i][m]++;
			}
		}

	}

	for (Int_t i=64; i<128; i++) {
		for (Int_t m=0; m<4; m++) {
			if ((TDC_LE_sft[i][m] > TDC_min_SFT) && (TDC_LE_sft[i][m] < TDC_max_SFT)) {
				counter_TDCB[i-64][m]++;
			}
		}

	}

///////////////////////////////////////////////////////

	for (Int_t i=0; i<64; i++) {
		for (Int_t m=0; m<4; m++) {
			if ((TDC_LE_target[i][m] > TDC_min_TARGET) && (TDC_LE_target[i][m] < TDC_max_TARGET)) {
				counter_TARGET_TDC_A[i][m]++;
			}
		}

	}

	for (Int_t i=64; i<128; i++) {
		for (Int_t m=0; m<4; m++) {
			if ((TDC_LE_target[i][m] > TDC_min_TARGET) && (TDC_LE_target[i][m] < TDC_max_TARGET)) {
				counter_TARGET_TDC_B[i-64][m]++;
			}
		}

	}

	for (Int_t i=128; i<192; i++) {
		for (Int_t m=0; m<4; m++) {
			if ((TDC_LE_target[i][m] > TDC_min_TARGET) && (TDC_LE_target[i][m] < TDC_max_TARGET)) {
				counter_TARGET_TDC_C[i-128][m]++;
			}
		}

	}

	for (Int_t i=192; i<256; i++) {
		for (Int_t m=0; m<4; m++) {
			if ((TDC_LE_target[i][m] > TDC_min_TARGET) && (TDC_LE_target[i][m] < TDC_max_TARGET)) {
				counter_TARGET_TDC_D[i-192][m]++;
			}
		}

	}
	
		
}

for(Int_t i=0; i<64; i++) h_TDC_0A->Fill(i,counter_TDCA[i][0]);
for(Int_t j=0; j<64; j++) h_TDC_0B->Fill(j+64,counter_TDCB[j][0]);

for(Int_t k=0; k<64; k++) h_TDC_1A->Fill(k,counter_TDCA[k][1]);
for(Int_t m=0; m<64; m++) h_TDC_1B->Fill(m+64,counter_TDCB[m][1]);

for(Int_t i=0; i<64; i++) h_TDC_2A->Fill(i,counter_TDCA[i][2]);
for(Int_t j=0; j<64; j++) h_TDC_2B->Fill(j+64,counter_TDCB[j][2]);

for(Int_t k=0; k<64; k++) h_TDC_3A->Fill(k,counter_TDCA[k][3]);
for(Int_t m=0; m<64; m++) h_TDC_3B->Fill(m+64,counter_TDCB[m][3]);

/////

for(Int_t i=0; i<64; i++) h_TDC_TARGET_0A->Fill(i,counter_TARGET_TDC_A[i][0]);
for(Int_t j=0; j<64; j++) h_TDC_TARGET_0B->Fill(j+64,counter_TARGET_TDC_B[j][0]);
for(Int_t k=0; k<64; k++) h_TDC_TARGET_0C->Fill(k+128,counter_TARGET_TDC_C[k][0]);
for(Int_t m=0; m<64; m++) h_TDC_TARGET_0D->Fill(m+192,counter_TARGET_TDC_D[m][0]);

for(Int_t i=0; i<64; i++) h_TDC_TARGET_1A->Fill(i,counter_TARGET_TDC_A[i][1]);
for(Int_t j=0; j<64; j++) h_TDC_TARGET_1B->Fill(j+64,counter_TARGET_TDC_B[j][1]);
for(Int_t k=0; k<64; k++) h_TDC_TARGET_1C->Fill(k+128,counter_TARGET_TDC_C[k][1]);
for(Int_t m=0; m<64; m++) h_TDC_TARGET_1D->Fill(m+192,counter_TARGET_TDC_D[m][1]);

for(Int_t i=0; i<64; i++) h_TDC_TARGET_2A->Fill(i,counter_TARGET_TDC_A[i][2]);
for(Int_t j=0; j<64; j++) h_TDC_TARGET_2B->Fill(j+64,counter_TARGET_TDC_B[j][2]);
for(Int_t k=0; k<64; k++) h_TDC_TARGET_2C->Fill(k+128,counter_TARGET_TDC_C[k][2]);
for(Int_t m=0; m<64; m++) h_TDC_TARGET_2D->Fill(m+192,counter_TARGET_TDC_D[m][2]);

for(Int_t i=0; i<64; i++) h_TDC_TARGET_3A->Fill(i,counter_TARGET_TDC_A[i][3]);
for(Int_t j=0; j<64; j++) h_TDC_TARGET_3B->Fill(j+64,counter_TARGET_TDC_B[j][3]);
for(Int_t k=0; k<64; k++) h_TDC_TARGET_3C->Fill(k+128,counter_TARGET_TDC_C[k][3]);
for(Int_t m=0; m<64; m++) h_TDC_TARGET_3D->Fill(m+192,counter_TARGET_TDC_D[m][3]);

//for(Int_t i=0; i<16; i++) {
//	cout << hlayer1->GetBinContent(i) << endl;
//}

TCanvas *c1;
c1 = new TCanvas(Name_Can,Title_Can,1400,1000); 
c1->Divide(2,3);
c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

TLegend *l1;
l1 = new TLegend(0.85,0.9,0.99,0.99);
l1->AddEntry(h_TDC_0A, "TDC_LE_SFT[0]", "F");
l1->AddEntry(h_TDC_1A, "TDC_LE_SFT[1]", "F");
l1->AddEntry(h_TDC_2A, "TDC_LE_SFT[2]", "F");
l1->AddEntry(h_TDC_3A, "TDC_LE_SFT[3]", "F");

TLegend *l2;
l2 = new TLegend(0.85,0.9,0.99,0.99);
l2->AddEntry(h_TDC_TARGET_0A, "TDC_LE_TARGET[0]", "F");
l2->AddEntry(h_TDC_TARGET_1A, "TDC_LE_TARGET[1]", "F");
l2->AddEntry(h_TDC_TARGET_2A, "TDC_LE_TARGET[2]", "F");
l2->AddEntry(h_TDC_TARGET_3A, "TDC_LE_TARGET[3]", "F");

c1->cd(1);
h_TDC_0A->SetLineColor(1);
h_TDC_0A->SetFillStyle(3005);
h_TDC_0A->SetFillColor(1);
h_TDC_0A->Draw();

h_TDC_1A->SetLineColor(2);
h_TDC_1A->SetFillStyle(3004);
h_TDC_1A->SetFillColor(2);
h_TDC_1A->Draw("same");

h_TDC_2A->SetLineColor(3);
h_TDC_2A->SetFillStyle(3005);
h_TDC_2A->SetFillColor(3);
h_TDC_2A->Draw("same");

h_TDC_3A->SetLineColor(4);
h_TDC_3A->SetFillStyle(3004);
h_TDC_3A->SetFillColor(4);
h_TDC_3A->Draw("same");
l1->Draw("same");

c1->cd(2);
h_TDC_0B->SetLineColor(1);
h_TDC_0B->SetFillStyle(3005);
h_TDC_0B->SetFillColor(1);
h_TDC_0B->Draw();

h_TDC_1B->SetLineColor(2);
h_TDC_1B->SetFillStyle(3004);
h_TDC_1B->SetFillColor(2);
h_TDC_1B->Draw("same");

h_TDC_2B->SetLineColor(3);
h_TDC_2B->SetFillStyle(3005);
h_TDC_2B->SetFillColor(3);
h_TDC_2B->Draw("same");

h_TDC_3B->SetLineColor(4);
h_TDC_3B->SetFillStyle(3004);
h_TDC_3B->SetFillColor(4);
h_TDC_3B->Draw("same");
l1->Draw("same");

c1->cd(3);
h_TDC_TARGET_0A->SetLineColor(1);
h_TDC_TARGET_0A->SetFillStyle(3005);
h_TDC_TARGET_0A->SetFillColor(1);
h_TDC_TARGET_0A->Draw();

h_TDC_TARGET_1A->SetLineColor(2);
h_TDC_TARGET_1A->SetFillStyle(3004);
h_TDC_TARGET_1A->SetFillColor(2);
h_TDC_TARGET_1A->Draw("same");

h_TDC_TARGET_2A->SetLineColor(3);
h_TDC_TARGET_2A->SetFillStyle(3005);
h_TDC_TARGET_2A->SetFillColor(3);
h_TDC_TARGET_2A->Draw("same");

h_TDC_TARGET_3A->SetLineColor(4);
h_TDC_TARGET_3A->SetFillStyle(3004);
h_TDC_TARGET_3A->SetFillColor(4);
h_TDC_TARGET_3A->Draw("same");
l2->Draw("same");

c1->cd(4);
h_TDC_TARGET_0B->SetLineColor(1);
h_TDC_TARGET_0B->SetFillStyle(3005);
h_TDC_TARGET_0B->SetFillColor(1);
h_TDC_TARGET_0B->Draw();

h_TDC_TARGET_1B->SetLineColor(2);
h_TDC_TARGET_1B->SetFillStyle(3004);
h_TDC_TARGET_1B->SetFillColor(2);
h_TDC_TARGET_1B->Draw("same");

h_TDC_TARGET_2B->SetLineColor(3);
h_TDC_TARGET_2B->SetFillStyle(3005);
h_TDC_TARGET_2B->SetFillColor(3);
h_TDC_TARGET_2B->Draw("same");

h_TDC_TARGET_3B->SetLineColor(4);
h_TDC_TARGET_3B->SetFillStyle(3004);
h_TDC_TARGET_3B->SetFillColor(4);
h_TDC_TARGET_3B->Draw("same");
l2->Draw("same");

c1->cd(5);
h_TDC_TARGET_0C->SetLineColor(1);
h_TDC_TARGET_0C->SetFillStyle(3005);
h_TDC_TARGET_0C->SetFillColor(1);
h_TDC_TARGET_0C->Draw();

h_TDC_TARGET_1C->SetLineColor(2);
h_TDC_TARGET_1C->SetFillStyle(3004);
h_TDC_TARGET_1C->SetFillColor(2);
h_TDC_TARGET_1C->Draw("same");

h_TDC_TARGET_2C->SetLineColor(3);
h_TDC_TARGET_2C->SetFillStyle(3005);
h_TDC_TARGET_2C->SetFillColor(3);
h_TDC_TARGET_2C->Draw("same");

h_TDC_TARGET_3C->SetLineColor(4);
h_TDC_TARGET_3C->SetFillStyle(3004);
h_TDC_TARGET_3C->SetFillColor(4);
h_TDC_TARGET_3C->Draw("same");
l2->Draw("same");

c1->cd(6);
h_TDC_TARGET_0D->SetLineColor(1);
h_TDC_TARGET_0D->SetFillStyle(3005);
h_TDC_TARGET_0D->SetFillColor(1);
h_TDC_TARGET_0D->Draw();

h_TDC_TARGET_1D->SetLineColor(2);
h_TDC_TARGET_1D->SetFillStyle(3004);
h_TDC_TARGET_1D->SetFillColor(2);
h_TDC_TARGET_1D->Draw("same");

h_TDC_TARGET_2D->SetLineColor(3);
h_TDC_TARGET_2D->SetFillStyle(3005);
h_TDC_TARGET_2D->SetFillColor(3);
h_TDC_TARGET_2D->Draw("same");

h_TDC_TARGET_3D->SetLineColor(4);
h_TDC_TARGET_3D->SetFillStyle(3004);
h_TDC_TARGET_3D->SetFillColor(4);
h_TDC_TARGET_3D->Draw("same");
l2->Draw("same");

}
