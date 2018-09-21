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
#include "ANAPATH.h"
#include "Thresholds.h"
#endif
void TDC_Multiplicity2(Int_t run_number=1, Int_t flag=0) {  

gStyle->SetOptStat(0);

Int_t ADC_High_target[256];			Int_t ADC_High_sft[128];
Int_t ADC_Low_target[256];			Int_t ADC_Low_sft[128];
Int_t TDC_LE_target[256][16];		Int_t TDC_LE_sft[128][16];
Int_t TDC_TE_target[256][16];		Int_t TDC_TE_sft[128][16];

Int_t ADC_tof1[24];
Int_t ADC_tof2[48];

Int_t mult_counter[149] = {0};
Int_t mult_counter2[149] = {0};

Int_t TARGET_empty_counter = 0;
Int_t TARGET_counter = 0;

char path_input[200];                   char file_mapping[200];
sprintf(path_input,path_merged);          sprintf(file_mapping,"../Mapping");
//sprintf(path_input,"/triumfcs/trshare/trek/E36/Data/April_2015/root");

char Name_finput[200];
sprintf(Name_finput,"%s/Run%dMS.root",path_input, run_number);

////

char par_finput3[200];
//sprintf(par_finput,"/triumfcs/trshare/trek/E36/Macros/April_2015/Mapping/%s",source_mapping);
sprintf(par_finput3,"%s/ADC_TARGET_Thresholds.txt",file_mapping);

Int_t par_temp_TARGET[2][256];
ifstream fdat3(par_finput3,ios::in);
for(Int_t ii=0; ii<256; ii++) fdat3 >> par_temp_TARGET[0][ii] >> par_temp_TARGET[1][ii];
fdat3.close();


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

fChain->SetBranchAddress("ADC_TOF1",ADC_tof1);

fChain->SetBranchAddress("ADC_TOF2",ADC_tof2);

TH1D *h_TDC_LE_TARGET;		char Title_h_TDC_LE_TARGET[100];	char Name_h_TDC_LE_TARGET[100];
TH1D *h_TDC_LE_TARGET2;		char Title_h_TDC_LE_TARGET2[100];	char Name_h_TDC_LE_TARGET2[100];
TH1D *h_TDC_LE_TARGET3;		char Title_h_TDC_LE_TARGET3[100];	char Name_h_TDC_LE_TARGET3[100];
TH1D *h_TDC_LE_TARGET4;		char Title_h_TDC_LE_TARGET4[100];	char Name_h_TDC_LE_TARGET4[100];
TH1D *h_TDC_LE_TARGET5;		char Title_h_TDC_LE_TARGET5[100];	char Name_h_TDC_LE_TARGET5[100];

sprintf(Title_h_TDC_LE_TARGET,"TDC_LE_TARGET[0] > 0 -- Run %d", run_number); 
sprintf(Name_h_TDC_LE_TARGET,"TDC_LE_TARGET[0] > 0"); 
h_TDC_LE_TARGET = new TH1D(Name_h_TDC_LE_TARGET,Title_h_TDC_LE_TARGET,150,-1,149);

sprintf(Title_h_TDC_LE_TARGET2,"%d < TDC_LE_TARGET[0] < %d -- Run %d",TDC_min_SFT, TDC_max_SFT, run_number); 
sprintf(Name_h_TDC_LE_TARGET2,"%d < TDC_LE_TARGET[0] < %d -- Run %d",TDC_min_SFT, TDC_max_SFT, run_number); 
h_TDC_LE_TARGET2 = new TH1D(Name_h_TDC_LE_TARGET2,Title_h_TDC_LE_TARGET2,150,-1,149);

sprintf(Title_h_TDC_LE_TARGET3,"ADC_High_TARGET > %d -- Run %d", ADC_cut_TARGET, run_number); 
sprintf(Name_h_TDC_LE_TARGET3,"ADC_High_TARGET > %d -- Run %d", ADC_cut_TARGET, run_number); 
h_TDC_LE_TARGET3 = new TH1D(Name_h_TDC_LE_TARGET3,Title_h_TDC_LE_TARGET3,7,0,7);

sprintf(Title_h_TDC_LE_TARGET4,"Empty Target"); 
sprintf(Name_h_TDC_LE_TARGET4,"Empty Target"); 
h_TDC_LE_TARGET4 = new TH1D(Name_h_TDC_LE_TARGET4,Title_h_TDC_LE_TARGET4,7,0,7);

sprintf(Title_h_TDC_LE_TARGET5,"1 to 5 Target Hits"); 
sprintf(Name_h_TDC_LE_TARGET5,"1 to 5 Target Hits"); 
h_TDC_LE_TARGET5 = new TH1D(Name_h_TDC_LE_TARGET5,Title_h_TDC_LE_TARGET5,7,0,7);

//read all entries and fill the histograms
Int_t nentries = (Int_t)fChain->GetEntries();
cout <<  "Total Number of Entries :     " <<  nentries << endl;
cout << "    " << endl;

if(flag!=0) nentries=flag;
for (Int_t i=0; i<nentries; i++) {
	fChain->GetEntry(i);
	
	Int_t TDC_Hit_Counter = 0;
	Int_t TDC_Hit_Counter2 = 0;
	Int_t ADC_Hit_Counter = 0;
	
	if(i%10000==1)
	cout<<"**** "<< i <<" events done"<<endl;
		
	for (Int_t j=0; j<256; j++) {
		if (TDC_LE_target[j][0] > 0) TDC_Hit_Counter++;
		if (TDC_LE_target[j][0] > TDC_min_SFT && TDC_LE_target[j][0] < TDC_max_SFT) TDC_Hit_Counter2++;
		if (ADC_High_target[j] > par_temp_TARGET[1][j]) ADC_Hit_Counter++;
	}

	mult_counter[TDC_Hit_Counter]++;
	mult_counter2[TDC_Hit_Counter2]++;

	if (ADC_Hit_Counter == 0) TARGET_empty_counter++;
	if (ADC_Hit_Counter >= 1 && ADC_Hit_Counter <= 5) TARGET_counter++;
}

for(Int_t q=0; q<300; q++) {
	h_TDC_LE_TARGET->Fill(q,mult_counter[q]);
	h_TDC_LE_TARGET2->Fill(q,mult_counter2[q]);
}

double total_entries = double(double(nentries)/double(nentries));
double emtpy_target = double(double(TARGET_empty_counter)/double(nentries));
double one_five_target = double(double(TARGET_counter)/double(nentries));

h_TDC_LE_TARGET3->Fill(1,total_entries);
h_TDC_LE_TARGET4->Fill(3,emtpy_target);
h_TDC_LE_TARGET5->Fill(5,one_five_target);

char Name_Can[100];			char Title_Can[100];
sprintf(Name_Can,"TDC_LE_TARGET -- Run %d -- TDC_Multiplicity2.C",run_number);
sprintf(Title_Can,"TDC_LE_TARGET -- Run %d  -- TDC_Multiplicity2.C",run_number);

char Name_Legend[100];
sprintf(Name_Legend,"%d < TDC < %d",TDC_min_SFT, TDC_max_SFT);

TCanvas *c1;
c1 = new TCanvas(Name_Can,Title_Can,500,300,1000,1000);
c1->Divide(2,2);
c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

TLegend *l1;
l1 = new TLegend(0.75,0.75,0.99,0.99);
l1->AddEntry(h_TDC_LE_TARGET3, "Total Entries", "F");
l1->AddEntry(h_TDC_LE_TARGET4, "Empty Target", "F");
l1->AddEntry(h_TDC_LE_TARGET5, "1 to 5 Target Hits", "F");

c1->cd(1);
//h_TDC_LE_TARGET->SetLineColor(4);
//h_TDC_LE_TARGET->SetFillStyle(3004);
//h_TDC_LE_TARGET->SetFillColor(4);
h_TDC_LE_TARGET->Draw();

c1->cd(2);
//h_TDC_LE_TARGET2->SetLineColor(2);
//h_TDC_LE_TARGET2->SetFillStyle(3005);
//h_TDC_LE_TARGET2->SetFillColor(2);
h_TDC_LE_TARGET2->Draw();

c1->cd(3);
h_TDC_LE_TARGET3->SetLineColor(4);
h_TDC_LE_TARGET3->SetFillStyle(3004);
h_TDC_LE_TARGET3->SetFillColor(4);
h_TDC_LE_TARGET3->Draw();

h_TDC_LE_TARGET4->SetLineColor(3);
h_TDC_LE_TARGET4->SetFillStyle(3004);
h_TDC_LE_TARGET4->SetFillColor(3);
h_TDC_LE_TARGET4->Draw("same");

h_TDC_LE_TARGET5->SetLineColor(2);
h_TDC_LE_TARGET5->SetFillStyle(3004);
h_TDC_LE_TARGET5->SetFillColor(2);
h_TDC_LE_TARGET5->Draw("same");
l1->Draw("same");

}
