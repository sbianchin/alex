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
//#include "Event_Display_MS.h"
#endif

void SFT_Efficiency_gap(Int_t run_number=1710, Int_t SFT_ADC_cut=950, Int_t SFT_TDC_min=750, Int_t SFT_TDC_max=850, Int_t NHit_TARGET_good_event=5, Int_t TARGET_ADC_cut=1000, Int_t Nb_events=-1, Int_t flag=0) {

gStyle->SetOptStat(0);

//Int_t cond=0;        // EVERY GOOD EVENT + COMBINATORY
Int_t cond=1; 	 // CLEAN EVENTS ONLY



Int_t count_p1_p2=0;

Int_t thr;						thr = SFT_ADC_cut;
Int_t TDCthr1;					TDCthr1 = SFT_TDC_min;
Int_t TDCthr2;					TDCthr2 = SFT_TDC_max;

//char* source_mapping="SFT_Mapping_Oct14.txt";  // Mapping file !!!

Int_t counter1 = 0;
Int_t counter2 = 0;
Int_t counter3 = 0;
Int_t counter4 = 0;
Int_t counter5 = 0;
Int_t counter6 = 0;
Int_t counter7 = 0;
Int_t counter8 = 0;

////////////////////
Int_t counter1a = 0;
Int_t counter2a = 0;
Int_t counter3a = 0;
Int_t counter4a = 0;
Int_t counter5a = 0;
Int_t counter6a = 0;
Int_t counter7a = 0;
Int_t counter8a = 0;
/////////////////////

Int_t helicity1or = 0;
Int_t helicity2or = 0;
Int_t helicity1and = 0;
Int_t helicity2and = 0;

/////////////////////
Int_t helicity1ora = 0;
Int_t helicity2ora = 0;
Int_t helicity1anda = 0;
Int_t helicity2anda = 0;
/////////////////////

Int_t count_or = 0;
Int_t count_and = 0;

Int_t count_hel_or = 0;
Int_t count_hel_and = 0;

///////////////////////

Int_t count_ora = 0;
Int_t count_anda = 0;

Int_t count_hel_ora = 0;
Int_t count_hel_anda = 0;

///////////////////////

Int_t onehitcount = 0;

//////////////////////
Int_t countlayer1or = 0;
Int_t countlayer1and = 0;
Int_t countlayer2or = 0;
Int_t countlayer2and = 0;
Int_t countlayer3or = 0;
Int_t countlayer3and = 0;
Int_t countlayer4or = 0;
Int_t countlayer4and = 0;
Int_t counthelicity1 = 0;
Int_t counthelicity2 = 0;

/////////////////////
Int_t countlayer1_2[15] = {0};
Int_t countlayer3_4[17] = {0};

Int_t counterp1[300] = {-1};		Int_t ip1=0;
Int_t counterp2[300] = {-1};		Int_t ip2=0;
Int_t counterp3[300] = {-1};		Int_t ip3=0;
Int_t counterp4[300] = {-1};		Int_t ip4=0;
Int_t p1p2[300]={-1};				Int_t ip12=0;
Int_t p3p4[300]={-1};				Int_t ip34=0;

Int_t counters1 = 0;
Int_t counters2 = 0;
Int_t counters3 = 0;
Int_t counters4 = 0;


Int_t p1_p2 = 0;
Int_t p3_p4 = 0;

Int_t count5 = 0;
Int_t count6 = 0;
Int_t count9 = 0;
Int_t count10 = 0;
//////////////////////

Int_t ADC_High_target[256];		Int_t ADC_High_sft[128];
Int_t ADC_Low_target[256];		Int_t ADC_Low_sft[128];
Int_t TDC_LE_target[256][16];		Int_t TDC_LE_sft[128][16];
Int_t TDC_TE_target[256][16];		Int_t TDC_TE_sft[128][16];

Int_t fiber[128]={-1};
Int_t fiberTDC[128]={-1};

Int_t ADC_tof1[24];
Int_t ADC_tof2[48];

char Name_Can[100];			char Title_Can[100];
sprintf(Name_Can,"Run %d -- SFT_Efficiency2.C",run_number);
sprintf(Title_Can,"Run %d  -- SFT_Efficiency2.C",run_number);

char par_map[200];		char par_finput[200];
//sprintf(par_finput,"/media/bianchin/hdd1/trek/E36/Macros/April_2015/Mapping/SFT_Mapping_Oct14.txt");
sprintf(par_map,"../Mapping");
sprintf(par_finput,"%s/SFT_Mapping_Oct14.txt",par_map);

char path_input[200];
//sprintf(path_input,"/media/bianchin/hdd1/trek/E36/Data/June_2015/root/Merged");
sprintf(path_input,path_merged);

char Name_finput[200];
sprintf(Name_finput,"%s/Run%dMS.root",path_input, run_number);

////

TH1D *h_count;	char Title_h_count[100];	char Name_h_count[100];
TH1D *h_count2;	char Title_h_count2[100];	char Name_h_count2[100];
TH1D *h_count3;	char Title_h_count3[100];	char Name_h_count3[100];
TH1D *h_count4;	char Title_h_count4[100];	char Name_h_count4[100];
TH1D *h_count5;	char Title_h_count5[100];	char Name_h_count5[100];
TH1D *h_count6;	char Title_h_count6[100];	char Name_h_count6[100];
TH1D *h_count7;	char Title_h_count7[100];	char Name_h_count7[100];
TH1D *h_count8;	char Title_h_count8[100];	char Name_h_count8[100];
TH1D *h_count9;	char Title_h_count9[100];	char Name_h_count9[100];
TH1D *h_count10; char Title_h_count10[100];	char Name_h_count10[100];
TH1D *h_count11; char Title_h_count11[100];	char Name_h_count11[100];
TH1D *h_count12; char Title_h_count12[100];	char Name_h_count12[100];
TH1D *h_count13; char Title_h_count13[100];	char Name_h_count13[100];
TH1D *h_count14; char Title_h_count14[100];	char Name_h_count14[100];
TH1D *h_count15; char Title_h_count15[100];	char Name_h_count15[100];
TH1D *h_count16; char Title_h_count16[100];	char Name_h_count16[100];
TH1D *h_count17; char Title_h_count17[100];	char Name_h_count17[100];
TH1D *h_count18; char Title_h_count18[100];	char Name_h_count18[100];
TH1D *h_count19; char Title_h_count19[100];	char Name_h_count19[100];
TH1D *h_count20; char Title_h_count20[100];	char Name_h_count20[100];
TH1D *h_good_events;  char Title_h_good_events[100];  char Name_h_good_events[100];

////////////////////////////////////////////

sprintf(Title_h_count,"SFT Efficiency - Run %d  |  ADC > %d , %d < TDC < %d", run_number, thr, TDCthr1, TDCthr2); 
sprintf(Name_h_count,"SFT/Helicity Counter"); 
h_count = new TH1D(Name_h_count,Title_h_count,7,0,7);

sprintf(Title_h_count2,"Run %d", run_number); 
sprintf(Name_h_count2,"SFT/Helicity Counter"); 
h_count2 = new TH1D(Name_h_count2,Title_h_count2,7,0,7);

sprintf(Title_h_count3,"Run %d", run_number); 
sprintf(Name_h_count3,"SFT/Helicity Counter"); 
h_count3 = new TH1D(Name_h_count3,Title_h_count3,7,0,7);

sprintf(Title_h_count4,"Run %d", run_number); 
sprintf(Name_h_count4,"SFT/Helicity Counter"); 
h_count4 = new TH1D(Name_h_count4,Title_h_count4,7,0,7);

sprintf(Title_h_count5,"Run %d", run_number); 
sprintf(Name_h_count5,"SFT/Helicity Counter"); 
h_count5 = new TH1D(Name_h_count5,Title_h_count5,7,0,7);

////

sprintf(Title_h_count6,"SFT Efficiency - Run %d  |  ADC > %d , %d < TDC < %d", run_number, thr, TDCthr1, TDCthr2); 
sprintf(Name_h_count6,"SFT/Helicity Counter"); 
h_count6 = new TH1D(Name_h_count6,Title_h_count6,7,0,7);

sprintf(Title_h_count7,"Run %d", run_number); 
sprintf(Name_h_count7,"SFT/Helicity Counter"); 
h_count7 = new TH1D(Name_h_count7,Title_h_count7,7,0,7);

sprintf(Title_h_count8,"Run %d", run_number); 
sprintf(Name_h_count8,"SFT/Helicity Counter"); 
h_count8 = new TH1D(Name_h_count8,Title_h_count8,7,0,7);

sprintf(Title_h_count9,"Run %d", run_number); 
sprintf(Name_h_count9,"SFT/Helicity Counter"); 
h_count9 = new TH1D(Name_h_count9,Title_h_count9,7,0,7);

sprintf(Title_h_count10,"Run %d", run_number); 
sprintf(Name_h_count10,"SFT/Helicity Counter"); 
h_count10 = new TH1D(Name_h_count10,Title_h_count10,7,0,7);

////

sprintf(Title_h_count11,"SFT Efficiency - Run %d  |  ADC > %d , %d < TDC < %d", run_number, thr, TDCthr1, TDCthr2); 
sprintf(Name_h_count11,"SFT/Helicity Counter"); 
h_count11 = new TH1D(Name_h_count11,Title_h_count11,7,0,7);

sprintf(Title_h_count12,"Run %d", run_number); 
sprintf(Name_h_count12,"SFT/Helicity Counter"); 
h_count12 = new TH1D(Name_h_count12,Title_h_count12,7,0,7);

sprintf(Title_h_count13,"Run %d", run_number); 
sprintf(Name_h_count13,"SFT/Helicity Counter"); 
h_count13 = new TH1D(Name_h_count13,Title_h_count13,7,0,7);

sprintf(Title_h_count14,"Run %d", run_number); 
sprintf(Name_h_count14,"SFT/Helicity Counter"); 
h_count14 = new TH1D(Name_h_count14,Title_h_count14,7,0,7);

sprintf(Title_h_count15,"Run %d", run_number); 
sprintf(Name_h_count15,"SFT/Helicity Counter"); 
h_count15 = new TH1D(Name_h_count15,Title_h_count15,7,0,7);

////
sprintf(Title_h_count16,"SFT Efficiency - Run %d  |  ADC > %d , %d < TDC < %d", run_number, thr, TDCthr1, TDCthr2); 
sprintf(Name_h_count16,"SFT/Helicity Counter"); 
h_count16 = new TH1D(Name_h_count16,Title_h_count16,7,0,7);

sprintf(Title_h_count17,"Run %d", run_number); 
sprintf(Name_h_count17,"SFT/Helicity Counter"); 
h_count17 = new TH1D(Name_h_count17,Title_h_count17,7,0,7);

sprintf(Title_h_count18,"Run %d", run_number); 
sprintf(Name_h_count18,"SFT/Helicity Counter"); 
h_count18 = new TH1D(Name_h_count18,Title_h_count18,7,0,7);

sprintf(Title_h_count19,"Run %d", run_number); 
sprintf(Name_h_count19,"SFT/Helicity Counter"); 
h_count19 = new TH1D(Name_h_count19,Title_h_count19,7,0,7);

sprintf(Title_h_count20,"Run %d", run_number); 
sprintf(Name_h_count20,"SFT/Helicity Counter"); 
h_count20 = new TH1D(Name_h_count20,Title_h_count20,7,0,7);

sprintf(Title_h_good_events,"SFT Efficiency - Run %d  |  ADC > %d , %d < TDC < %d", run_number, thr, TDCthr1, TDCthr2); 
sprintf(Name_h_good_events,"SFT/Helicity Counter"); 
h_good_events = new TH1D(Name_h_good_events,Title_h_good_events,7,0,7);
////

TH1D *hlayer1_2[12];	char Title_hlayer1_2[12][100];	char Name_hlayer1_2[12][100];
TH1D *hlayer3_4[12];	char Title_hlayer3_4[12][100];	char Name_hlayer3_4[12][100];

for(Int_t ii=0; i<12; ii++){
sprintf(Title_hlayer1_2[ii],"SFT Layer Distances  |  Run %d - Layers 1 & 2 (Gap %d)", run_number,ii+1); 
sprintf(Name_hlayer1_2[ii],"Layer 1"); 
hlayer1_2[ii] = new TH1D(Name_hlayer1_2[ii],Title_hlayer1_2[ii],15,0,15);

sprintf(Title_hlayer3_4[ii],"SFT Layer Distances  |  Run %d - Layers 3 & 4 (Gap %d)", run_number,ii+1); 
sprintf(Name_hlayer3_4[ii],"Layer 2"); 
hlayer3_4[ii] = new TH1D(Name_hlayer3_4[ii],Title_hlayer3_4[ii],17,0,17);

////

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

cout << "   " << endl;
cout << "Mapping File:   " << par_finput << endl;
cout << "" << endl;

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

Int_t par_temp[2][128];
ifstream fdat(par_finput,ios::in);
for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
fdat.close();

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
Int_t good_event=0;
Int_t nentries = (Int_t)fChain->GetEntries();
cout <<  "Total Number of Events :     " <<  nentries << endl;
cout << " " << endl;
cout << " " << endl;

if(Nb_events!=-1) nentries=Nb_events;
for (Int_t j=0; j<nentries; j++) {

ip1=0;
ip2=0;
ip3=0;
ip4=0;
ip12=0;
ip34=0;


	if(j%10000==1)
	cout<<"**** "<< j <<" events done"<<endl;
	
//	fChain_SFT->GetEntry(j);
	fChain->GetEntry(j);
	counter1 = 0;
	counter2 = 0;
	counter3 = 0;
	counter4 = 0;
	counter5 = 0;
	counter6 = 0;
	counter7 = 0;
	counter8 = 0;
	helicity1or = 0;
	helicity2or = 0;
	helicity1and = 0;
	helicity2and = 0;

	counter1a = 0;
	counter2a = 0;
	counter3a = 0;
	counter4a = 0;
	counter5a = 0;
	counter6a = 0;
	counter7a = 0;
	counter8a = 0;
	helicity1ora = 0;
	helicity2ora = 0;
	helicity1anda = 0;
	helicity2anda = 0;

	for(Int_t iz=0; iz<128; iz++){
		counterp1[iz]=-1;
		counterp2[iz]=-1;
		counterp3[iz]=-1;
		counterp4[iz]=-1;
		p1p2[iz]=-1;
		p3p4[iz]=-1;
	}

	//counterp1 = 0;
	//counterp2 = 0;
	//counterp3 = 0;
	//counterp4 = 0;

	p1_p2 = 0;
	p3_p4 = 0;

	Int_t count_target=0;
	//Int_t good_event=0;

	//Condition on the TARGET
	for(Int_t it=0; it<256; it++){
		if(ADC_High_target[it]>TARGET_ADC_cut) count_target++;

	}

	//cout << j << "   " << count_target << endl;


	if(count_target>=NHit_TARGET_good_event && (ADC_tof1[0]>400 || ADC_tof1[12]>400)){
		good_event++;
	
	// Variable Change
	for(Int_t ii=0; ii<128; ii++){
		fiber[ii]=ADC_High_sft[par_temp[SFT_Efficiency_gap.C1][ii]];
		if(TDC_LE_sft[par_temp[1][ii]][0]>SFT_TDC_min && TDC_LE_sft[par_temp[1][ii]][0]<SFT_TDC_max)
			fiberTDC[ii]=TDC_LE_sft[par_temp[1][ii]][0];
		else if(TDC_LE_sft[par_temp[1][ii]][1]>SFT_TDC_min && TDC_LE_sft[par_temp[1][ii]][1]<SFT_TDC_max)
			fiberTDC[ii]=TDC_LE_sft[par_temp[1][ii]][1];
		else if(TDC_LE_sft[par_temp[1][ii]][2]>SFT_TDC_min && TDC_LE_sft[par_temp[1][ii]][2]<SFT_TDC_max)
			fiberTDC[ii]=TDC_LE_sft[par_temp[1][ii]][2];
		else if(TDC_LE_sft[par_temp[1][ii]][3]>SFT_TDC_min && TDC_LE_sft[par_temp[1][ii]][3]<SFT_TDC_max)
			fiberTDC[ii]=TDC_LE_sft[par_temp[1][ii]][3];
		else if(TDC_LE_sft[par_temp[1][ii]][4]>SFT_TDC_min && TDC_LE_sft[par_temp[1][ii]][4]<SFT_TDC_max)
			fiberTDC[ii]=TDC_LE_sft[par_temp[1][ii]][4];
		else if(TDC_LE_sft[par_temp[1][ii]][5]>SFT_TDC_min && TDC_LE_sft[par_temp[1][ii]][5]<SFT_TDC_max)
			fiberTDC[ii]=TDC_LE_sft[par_temp[1][ii]][5];
		else
			fiberTDC[ii]=TDC_LE_sft[par_temp[1][ii]][6];

		//	cout << j << "  " << count_target << "  " << ii << "   " << TDC_LE_sft[par_temp[1][ii]][0] << "   " << TDC_LE_sft[par_temp[1][ii]][1] << "   " << TDC_LE_sft[par_temp[1][ii]][2] << " ||  " << fiberTDC[ii] << endl;;
	
	}


/*
	// Variable Change
	for(Int_t ii=0; ii<128; ii++){
		fiber[ii]=ADC_High_sft[par_temp[1][ii]];
		fiberTDC[ii]=TDC_LE_sft[par_temp[1][ii]][0];
		cout << ii << "   " << TDC_LE_sft[par_temp[1][ii]][0] << "   " << TDC_LE_sft[par_temp[1][ii]][1] << endl;
	}
*/

	
//	for (Int_t j=0; j<64; j++) {
//	cout << j << " " << ADC_High_sft[j] << " " << ADC_High_sft[j+64] << endl;
//	}


///Layer 1
	for (Int_t i=0; i<15; i++) {
		if (((fiber[i]>thr) && (fiberTDC[i]>TDCthr1) && (fiberTDC[i]<TDCthr2)) || ((fiberTDC[i+64]>TDCthr1) && (fiberTDC[i+64]<TDCthr2) && (fiber[i+64]>thr))) {
			counter1++;
			counters1=i;
			counterp1[ip1] = i+1;
			ip1++;
		}
		if (((fiber[i]>thr) && (fiberTDC[i]>TDCthr1) && (fiberTDC[i]<TDCthr2)) && ((fiberTDC[i+64]>TDCthr1) && (fiberTDC[i+64]<TDCthr2) && (fiber[i+64]>thr))) {
			counter5++;
		}
	}

///Layer 2
	for (Int_t i=15; i<30; i++) {
		if (((fiber[i]>thr) && (fiberTDC[i]>TDCthr1) && (fiberTDC[i]<TDCthr2)) || ((fiberTDC[i+64]>TDCthr1) && (fiberTDC[i+64]<TDCthr2) && (fiber[i+64]>thr))) {
			counter2++;
			counters2=i;
			counterp2[ip2] = (i-15)+1;
			ip2++;
		}
		if (((fiber[i]>thr) && (fiberTDC[i]>TDCthr1) && (fiberTDC[i]<TDCthr2)) && ((fiberTDC[i+64]>TDCthr1) && (fiberTDC[i+64]<TDCthr2) && (fiber[i+64]>thr))) {
			counter6++;
		}
	}
 
///Layer 3
	for (Int_t i=30; i<47; i++) {
		if (((fiber[i]>thr) && (fiberTDC[i]>TDCthr1) && (fiberTDC[i]<TDCthr2)) || ((fiberTDC[i+64]>TDCthr1) && (fiberTDC[i+64]<TDCthr2) && (fiber[i+64]>thr))) {
			counter3++;
			counters3=i;
			counterp3[ip3] = (i-30)+1;
			ip3++;
		}
		if (((fiber[i]>thr) && (fiberTDC[i]>TDCthr1) && (fiberTDC[i]<TDCthr2)) && ((fiberTDC[i+64]>TDCthr1) && (fiberTDC[i+64]<TDCthr2) && (fiber[i+64]>thr))) {
			counter7++;
		}
	}

///Layer 4
	for (Int_t i=47; i<64; i++) {
		if (((fiber[i]>thr) && (fiberTDC[i]>TDCthr1) && (fiberTDC[i]<TDCthr2)) || ((fiberTDC[i+64]>TDCthr1) && (fiberTDC[i+64]<TDCthr2) && (fiber[i+64]>thr))) {
			counter4++;
			counters4=i;
			counterp4[ip4] = (i-47)+1;
			ip4++;
		}
		if (((fiber[i]>thr) && (fiberTDC[i]>TDCthr1) && (fiberTDC[i]<TDCthr2)) && ((fiberTDC[i+64]>TDCthr1) && (fiberTDC[i+64]<TDCthr2) && (fiber[i+64]>thr))) {
			counter8++;
		}
	}

	
	if(cond==1){
		if ((counter1 == 1) && (counter2 == 1) && (counter3 == 1) && (counter4 == 1)) {
			onehitcount++;
//			cout << "One Hit Per Layer: Event " << j << endl;
			p1_p2 = abs((counters2-15)-counters1);
			p3_p4 = abs((counters4-17)-counters3);
			cout << j << "   "  << counters1+1 << "  " << counters2-15+1 << "  " << p1_p2 << "  " << counters3 << "  " << counters4-17+1 << "  " << p3_p4 << endl;
			count_p1_p2++;
			(countlayer1_2[p1_p2])++;
			(countlayer3_4[p3_p4])++;
		}
	}

		// for(Int_t j=0; j<128; j++) cout << j << "     " << counterp1[j] << "  " << counterp2[j] << "  " << counterp3[j] << "  " << counterp4[j] << endl;
		//cout << j << "  " << count_p1_p2 << "   p1_p2  :  " << p1_p2 << endl;

		if(p1_p2==5) count5++;
		if(p1_p2==6) count6++;
		if(p1_p2==9) count9++;
		if(p1_p2==10) count10++;

		////////(p1 - p2) Calculation (with combinatory)
		for(Int_t i1=0; i1<128; i1++){
			if(counterp1[i1]>=0){
				for(Int_t i2=0; i2<128; i2++){
					if(counterp2[i2]>=0){
					p1p2[ip12] = abs(counterp1[i1]-counterp2[i2]);
					if(flag==1) cout << "TEST12: " << counterp1[i1] << "  " << counterp2[i2] << "  " << ip12 << "   " << p1p2[ip12] << endl;
					ip12++;
					}
				}
			}
		}
		////////////////////////////////////////////////////////////////////////////
		if(flag==1) cout << " " << endl;
		////////(p3 - p4) Calculation (with combinatory)
		for(Int_t i3=0; i3<128; i3++){
			if(counterp3[i3]>=0){
				for(Int_t i4=0; i4<128; i4++){
					if(counterp4[i4]>=0){
					p3p4[ip34] = abs(counterp3[i3]-counterp4[i4]);
					if(flag==1) cout << "TEST34: " << counterp3[i3] << "  " << counterp4[i4] << "  " << ip34 << "   " << p3p4[ip34] << endl;
					ip34++;
					//cout << "ip34 " << ip34 << endl;
					}
				}
			}
		}
		////////////////////////////////////////////////////////////////////////////
		if(flag==1) cout << " " << endl;

	if(cond==0){
		for(Int_t jj=0; jj<128; jj++){
			if(p1p2[jj]>=0){
				hlayer1_2->Fill(p1p2[jj]);
				if(flag==1) cout << jj+1 << "   |p1 - p2|   :  " << p1p2[jj] << endl;	
			}
		}
		
		if(flag==1) cout << "  " << endl;

		for(Int_t kk=0; kk<128; kk++){
			if(p3p4[kk]>=0){
				hlayer3_4->Fill(p3p4[kk]);
				if(flag==1) cout << kk+1 << "   |p3 - p4|   :  " << p3p4[kk] << endl;	
			}
		}
	}
		

		//if(cond==1){
		//	for(Int_t i=0; i<15; i++) hlayer1_2->Fill(i,countlayer1_2[i]);
		//	for(Int_t j=0; j<17; j++) hlayer3_4->Fill(j,countlayer3_4[j]);
		//}

		///////////


		//if(p1_p2==9 || p1_p2==10) cout << "p1 - p2 = 9 or 10 |  Evt: " << j << endl;

		//cout << j << "    p1 - p2  =  " << p1_p2 << "   ||   p3 - p4 =  " << p3_p4 << endl;
//		cout << p3_p4 << endl;
//		(countlayer1_2[p1_p2])++;
		//(countlayer3_4[p3_p4])++;
	//}
//	if (or_all) { 
//		count_or++;
//		cout << "SFT-OR " << j << endl;
//	}
//	if (and_all) {
//		count_and++;
//		cout << "SFT-AND " << j << endl;
//	}
//
//// Helicity Counter
	if ((helicity2or > 0) && (helicity1or > 0)) count_hel_or++;
	if ((helicity2and > 0) && (helicity1and > 0)) count_hel_and++;;

//	if (or_helicity_all) { 
//		count_hel_or++;
//		cout << "Hel-OR: " << j << endl;
//	}
	
//	if (and_helicity_all) {
//		count_hel_and++;
//		cout << "Hel-AND: " << j << endl;
//	}

////////////////NO TDC CUT
	if ((counter1a > 0) && (counter2a > 0) && (counter3a > 0) && (counter4a > 0)) count_ora++;
	if ((counter5a > 0) && (counter6a > 0) && (counter7a > 0) && (counter8a > 0)) count_anda++;
	if ((helicity2ora > 0) && (helicity1ora > 0)) count_hel_ora++;
	if ((helicity2anda > 0) && (helicity1anda > 0)) count_hel_anda++;

} // Endif
//cout << "Number of Good Events:  " << NHit_TARGET_good_event << endl;

} // EndLoop over the Event

	if(cond==1){
	//	if(ADC_tof1[0]>400 || ADC_tof1[12]>400){
		for(Int_t ii=0; ii<12; ii++){
			for(Int_t i=0; i<15; i++) hlayer1_2[ii]->Fill(i,countlayer1_2[i]);
			for(Int_t j=0; j<17; j++) hlayer3_4[ii]->Fill(j,countlayer3_4[j]);
		}
	}


cout << " " << endl;
cout << " " << endl;
cout << "Number of Good Events:  " << good_event << " / " << nentries << "   --->  " << 100*(double(good_event)/double(nentries)) << "%" << endl;
cout << " " << endl;
cout << " " << endl;

//cout << "p1 - p2 = 5 : " << count5 << " events" << endl;
//cout << "p1 - p2 = 6 : " << count6 << " events" << endl;
//cout << "p1 - p2 = 9 : " << count9 << " events" << endl;
//cout << "p1 - p2 = 10 : " << count10 << " events" << endl;

//for(Int_t i=0; i<15; i++) hlayer1_2->Fill(i,countlayer1_2[i]);
//for(Int_t j=0; j<17; j++) hlayer3_4->Fill(j,countlayer3_4[j]);
//cout << "Individual SFT counter" << endl;
//cout << "OR  " << "AND" << endl;
//cout << counter1 << "   " << counter5 << endl;
//cout << counter2 << "   " << counter6 << endl;
//cout << counter3 << "   " << counter7 << endl;
//cout << counter4 << "   " << counter8 << endl;

//cout << "Individual Helicity counter" << endl;
//cout << "OR  " << "AND" << endl;
//cout << helicity1or << "   " << helicity1and << endl;
//cout << helicity2or << "   " << helicity2and << endl;

//cout << "Event SFT counter" << endl;
//cout << "OR  " << "AND" << endl;
//cout << count_or << "   " << count_and << endl;

//cout << "Event Helicity counter" << endl;
//cout << "OR  " << "AND" << endl;
//cout << count_hel_or << "   " << count_hel_and << endl;

//cout << "One hit for each layer count: " << onehitcount << endl;

double eff_SFT_or = (double(count_or)/double(good_event))*100;
double eff_SFT_and = (double(count_and)/double(good_event))*100;
double eff_hel_or = (double(count_hel_or)/double(good_event))*100;
double eff_hel_and = (double(count_hel_and)/double(good_event))*100;

//double eff_l1_or = (double(countlayer1or)/double(counthelicity1))*100;
//double eff_l1_and = (double(countlayer1and)/double(counthelicity1))*100;
//double eff_l2_or = (double(countlayer2or)/double(counthelicity1))*100;
//double eff_l2_and = (double(countlayer2and)/double(counthelicity1))*100;

//double eff_l3_or = (double(countlayer3or)/double(counthelicity2))*100;
//double eff_l3_and = (double(countlayer3and)/double(counthelicity2))*100;
//double eff_l4_or = (double(countlayer4or)/double(counthelicity2))*100;
//double eff_l4_and = (double(countlayer4and)/double(counthelicity2))*100;

//cout << "counterhlicity1  :  " << counthelicity1 << endl;
//cout << "nentries  :  " << nentries << endl;

double eff_l1_or = (double(countlayer1or)/double(good_event))*100;
double eff_l1_and = (double(countlayer1and)/double(good_event))*100;
double eff_l2_or = (double(countlayer2or)/double(good_event))*100;
double eff_l2_and = (double(countlayer2and)/double(good_event))*100;

double eff_l3_or = (double(countlayer3or)/double(good_event))*100;
double eff_l3_and = (double(countlayer3and)/double(good_event))*100;
double eff_l4_or = (double(countlayer4or)/double(good_event))*100;
double eff_l4_and = (double(countlayer4and)/double(good_event))*100;


double eff_SFT_ora = (double(count_or)/double(nentries))*100;
double eff_SFT_anda = (double(count_and)/double(nentries))*100;
double eff_hel_ora = (double(count_hel_or)/double(nentries))*100;
double eff_hel_anda = (double(count_hel_and)/double(nentries))*100;

double scale = 0.01; 
eff_SFT_or = floor(eff_SFT_or / scale + 0.5) * scale;
eff_SFT_and = floor(eff_SFT_and / scale + 0.5) * scale;
eff_hel_or = floor(eff_hel_or / scale + 0.5) * scale;
eff_hel_and = floor(eff_hel_and / scale + 0.5) * scale;

eff_l1_or = floor(eff_l1_or / scale + 0.5) * scale;
eff_l1_and = floor( eff_l1_and / scale + 0.5) * scale;
eff_l2_or = floor(eff_l2_or / scale + 0.5) * scale;
eff_l2_and = floor(eff_l2_and / scale + 0.5) * scale;

eff_l3_or = floor(eff_l3_or / scale + 0.5) * scale;
eff_l3_and = floor( eff_l3_and / scale + 0.5) * scale;
eff_l4_or = floor(eff_l4_or / scale + 0.5) * scale;
eff_l4_and = floor(eff_l4_and / scale + 0.5) * scale;

eff_SFT_ora = floor(eff_SFT_ora / scale + 0.5) * scale;
eff_SFT_anda = floor(eff_SFT_anda / scale + 0.5) * scale;
eff_hel_ora = floor(eff_hel_ora / scale + 0.5) * scale;
eff_hel_anda = floor(eff_hel_anda / scale + 0.5) * scale;

/////////////////////////////////////////////////////

char str1[100];		char str9[100];
char str2[100];		char str10[100];
char str3[100];		char str11[100];
char str4[100];		char str12[100];
char str5[100];		char str13[100];
char str6[100];		char str14[100];
char str7[100];		char str15[100];
char str8[100];		char str16[100];

sprintf(str1,"%4.2f%%",eff_SFT_or);
sprintf(str2,"%4.2f%%",eff_SFT_and);
sprintf(str3,"%4.2f%%",eff_hel_or);
sprintf(str4,"%4.2f%%",eff_hel_and);
sprintf(str5,"%4.2f%%",eff_l1_or);
sprintf(str6,"%4.2f%%",eff_l1_and);
sprintf(str7,"%4.2f%%",eff_l2_or);
sprintf(str8,"%4.2f%%",eff_l2_and);
sprintf(str9,"%4.2f%%",eff_l3_or);
sprintf(str10,"%4.2f%%",eff_l3_and);
sprintf(str11,"%4.2f%%",eff_l4_or);
sprintf(str12,"%4.2f%%",eff_l4_and);
sprintf(str13,"%4.2f%%",eff_SFT_ora);
sprintf(str14,"%4.2f%%",eff_SFT_anda);
sprintf(str15,"%4.2f%%",eff_hel_ora);
sprintf(str16,"%4.2f%%",eff_hel_anda);

///

string SFTorEff = str1;
string SFTandEff = str2;
string HelorEff = str3;
string HelandEff = str4;

string L1orEff = str5;
string L1andEff = str6;
string L2orEff = str7;
string L2andEff = str8;

string L3orEff = str9;
string L3andEff = str10;
string L4orEff = str11;
string L4andEff = str12;

string SFTorEffa = str13;
string SFTandEffa = str14;
string HelorEffa = str15;
string HelandEffa = str16;

///
h_count->Fill(1,good_event);
h_count2->Fill(3,count_hel_or);
h_count3->Fill(3,count_hel_and);
h_count4->Fill(5,count_or);
h_count5->Fill(5,count_and);

//h_count6->Fill(1,counthelicity1);
h_count6->Fill(1,good_event);
h_count7->Fill(3,countlayer1or);
h_count8->Fill(3,countlayer1and);
h_count9->Fill(5,countlayer2or);
h_count10->Fill(5,countlayer2and);

//h_count11->Fill(1,counthelicity2);
h_count11->Fill(1,good_event);
h_count12->Fill(3,countlayer3or);
h_count13->Fill(3,countlayer3and);
h_count14->Fill(5,countlayer4or);
h_count15->Fill(5,countlayer4and);

/// No TDC Cut

h_count16->Fill(1,nentries);
h_count17->Fill(3,count_hel_or);
h_count18->Fill(3,count_hel_and);
h_count19->Fill(5,count_or);
h_count20->Fill(5,count_and);

///

h_good_events->Fill(1,good_event);

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
h_count7->GetXaxis()->SetLabelOffset(999);
h_count7->GetXaxis()->SetLabelSize(0);
h_count8->GetXaxis()->SetLabelOffset(999);
h_count8->GetXaxis()->SetLabelSize(0);
h_count9->GetXaxis()->SetLabelOffset(999);
h_count9->GetXaxis()->SetLabelSize(0);
h_count10->GetXaxis()->SetLabelOffset(999);
h_count10->GetXaxis()->SetLabelSize(0);

h_count11->GetXaxis()->SetLabelOffset(999);
h_count11->GetXaxis()->SetLabelSize(0);
h_count12->GetXaxis()->SetLabelOffset(999);
h_count12->GetXaxis()->SetLabelSize(0);
h_count13->GetXaxis()->SetLabelOffset(999);
h_count13->GetXaxis()->SetLabelSize(0);
h_count14->GetXaxis()->SetLabelOffset(999);
h_count14->GetXaxis()->SetLabelSize(0);
h_count15->GetXaxis()->SetLabelOffset(999);
h_count15->GetXaxis()->SetLabelSize(0);

h_count16->GetXaxis()->SetLabelOffset(999);
h_count16->GetXaxis()->SetLabelSize(0);
h_count17->GetXaxis()->SetLabelOffset(999);
h_count17->GetXaxis()->SetLabelSize(0);
h_count18->GetXaxis()->SetLabelOffset(999);
h_count18->GetXaxis()->SetLabelSize(0);
h_count19->GetXaxis()->SetLabelOffset(999);
h_count19->GetXaxis()->SetLabelSize(0);
h_count20->GetXaxis()->SetLabelOffset(999);
h_count20->GetXaxis()->SetLabelSize(0);
h_good_events->GetXaxis()->SetLabelOffset(999);
h_good_events->GetXaxis()->SetLabelSize(0);

//////////////////
double expn = 0.75;
double expn2 = 0.4;

double scl_hel_or = double(count_hel_or)*expn;
double scl_hel_and = double(count_hel_and)*expn2;
double scl_count_or = double(count_or)*expn;
double scl_count_and = double(count_and)*expn2;

double scl_layer1or = double(countlayer1or)*expn;
double scl_layer1and = double(countlayer1and)*expn2;
double scl_layer2or = double(countlayer2or)*expn;
double scl_layer2and = double(countlayer2and)*expn2;

double scl_layer3or = double(countlayer3or)*expn;
double scl_layer3and = double(countlayer3and)*expn2;
double scl_layer4or = double(countlayer4or)*expn;
double scl_layer4and = double(countlayer4and)*expn2;

double scl_hel_ora = double(count_hel_ora)*expn;
double scl_hel_anda = double(count_hel_anda)*expn2;
double scl_count_ora = double(count_ora)*expn;
double scl_count_anda = double(count_anda)*expn2;

double scl_labels1 = double(good_event)/-10.0;
double scl_labels2 = double(nentries)/-10.0;

////////////////////////////////////////////////////

TCanvas *c1;
c1 = new TCanvas(Name_Can,Title_Can,1000,900); 
c1->Divide(2,3);
c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

c1->cd(3);

TLegend *l1;
l1 = new TLegend(0.75,0.75,0.99,0.99);
l1->AddEntry(h_good_events, "\"Good Events\"", "F");
l1->AddEntry(h_count2, "Helicity | Up OR Down", "F");
l1->AddEntry(h_count3, "Helicity | Up AND Down", "F");
l1->AddEntry(h_count4, "Layer | Up OR Down", "F");
l1->AddEntry(h_count5, "Layer | Up AND Down", "F");

h_good_events->SetLineColor(1);
h_good_events->SetFillStyle(3008);
h_good_events->SetFillColor(1);
h_good_events->Draw();

h_count2->SetLineColor(2);
h_count2->SetFillStyle(3004);
h_count2->SetFillColor(2);
h_count2->Draw("same");

h_count3->SetLineColor(7);
h_count3->SetFillStyle(3008);
h_count3->SetFillColor(7);
h_count3->Draw("same");

h_count4->SetLineColor(6);
h_count4->SetFillStyle(3004);
h_count4->SetFillColor(6);
h_count4->Draw("same");

h_count5->SetLineColor(3);
h_count5->SetFillStyle(3008);
h_count5->SetFillColor(3);
h_count5->Draw("same");

l1->Draw("same");
TLatex *tex4 = new TLatex(0.1,scl_labels1,"\"Good Events\"");
tex4->SetLineWidth(2);
tex4->Draw();
TLatex *tex5 = new TLatex(3.087098,scl_labels1,"Helicity");
tex5->SetLineWidth(2);
tex5->Draw();
TLatex *tex6 = new TLatex(5.195532,scl_labels1,"Layer");
tex6->SetLineWidth(2);
tex6->Draw();

TLatex *tex7 = new TLatex(3.1,scl_hel_or,HelorEff.c_str());
tex7->SetLineWidth(2);
tex7->Draw();
TLatex *tex8 = new TLatex(3.1,scl_hel_and,HelandEff.c_str());
tex8->SetLineWidth(2);
tex8->Draw();
TLatex *tex9 = new TLatex(5.1,scl_count_or,SFTorEff.c_str());
tex9->SetLineWidth(2);
tex9->Draw();
TLatex *tex10 = new TLatex(5.1,scl_count_and,SFTandEff.c_str());
tex10->SetLineWidth(2);
tex10->Draw();

///////////////////////////////////////////////////////////////////
c1->cd(1);

TLegend *l2;
l2 = new TLegend(0.75,0.75,0.99,0.99);
l2->AddEntry(h_good_events, "\"Good Events\"", "F");
l2->AddEntry(h_count2, "Layer 1 | Up OR Down", "F");
l2->AddEntry(h_count3, "Layer 1 | Up AND Down", "F");
l2->AddEntry(h_count4, "Layer 2 | Up OR Down", "F");
l2->AddEntry(h_count5, "Layer 2 | Up AND Down", "F");

h_good_events->SetLineColor(1);
h_good_events->SetFillStyle(3008);
h_good_events->SetFillColor(1);
h_good_events->Draw();

h_count7->SetLineColor(2);
h_count7->SetFillStyle(3004);
h_count7->SetFillColor(2);
h_count7->Draw("same");

h_count8->SetLineColor(7);
h_count8->SetFillStyle(3008);
h_count8->SetFillColor(7);
h_count8->Draw("same");

h_count9->SetLineColor(6);
h_count9->SetFillStyle(3004);
h_count9->SetFillColor(6);
h_count9->Draw("same");

h_count10->SetLineColor(3);
h_count10->SetFillStyle(3008);
h_count10->SetFillColor(3);
h_count10->Draw("same");

l2->Draw("same");
TLatex *tex11 = new TLatex(0.1,scl_labels1,"\"Good Events\"");
tex11->SetLineWidth(2);
tex11->Draw();
TLatex *tex12 = new TLatex(3.087098,scl_labels1,"Layer 1");
tex12->SetLineWidth(2);
tex12->Draw();
TLatex *tex13 = new TLatex(5.1,scl_labels1,"Layer 2");
tex13->SetLineWidth(2);
tex13->Draw();

TLatex *tex14 = new TLatex(3.1,scl_layer1or,L1orEff.c_str());
tex14->SetLineWidth(2);
tex14->Draw();
TLatex *tex15 = new TLatex(3.1,scl_layer1and,L1andEff.c_str());
tex15->SetLineWidth(2);
tex15->Draw();
TLatex *tex16 = new TLatex(5.1,scl_layer2or,L2orEff.c_str());
tex16->SetLineWidth(2);
tex16->Draw();
TLatex *tex17 = new TLatex(5.1,scl_layer2and,L2andEff.c_str());
tex17->SetLineWidth(2);
tex17->Draw();

///////////////////////////////////////////////////////
c1->cd(2);

TLegend *l3;
l3 = new TLegend(0.75,0.75,0.99,0.99);
l3->AddEntry(h_good_events, "\"Good Events\"", "F");
l3->AddEntry(h_count2, "Layer 3 | Up OR Down", "F");
l3->AddEntry(h_count3, "Layer 3 | Up AND Down", "F");
l3->AddEntry(h_count4, "Layer 4 | Up OR Down", "F");
l3->AddEntry(h_count5, "Layer 4 | Up AND Down", "F");

h_good_events->SetLineColor(1);
h_good_events->SetFillStyle(3008);
h_good_events->SetFillColor(1);
h_good_events->Draw();

h_count12->SetLineColor(2);
h_count12->SetFillStyle(3004);
h_count12->SetFillColor(2);
h_count12->Draw("same");

h_count13->SetLineColor(7);
h_count13->SetFillStyle(3008);
h_count13->SetFillColor(7);
h_count13->Draw("same");

h_count14->SetLineColor(6);
h_count14->SetFillStyle(3004);
h_count14->SetFillColor(6);
h_count14->Draw("same");

h_count15->SetLineColor(3);
h_count15->SetFillStyle(3008);
h_count15->SetFillColor(3);
h_count15->Draw("same");


l3->Draw("same");
TLatex *tex18 = new TLatex(0.1,scl_labels1,"\"Good Events\"");
tex18->SetLineWidth(2);
tex18->Draw();
TLatex *tex19 = new TLatex(3.087098,scl_labels1,"Layer 3");
tex19->SetLineWidth(2);
tex19->Draw();
TLatex *tex20 = new TLatex(5.1,scl_labels1,"Layer 4");
tex20->SetLineWidth(2);
tex20->Draw();

TLatex *tex21 = new TLatex(3.1,scl_layer3or,L3orEff.c_str());
tex21->SetLineWidth(2);
tex21->Draw();
TLatex *tex22 = new TLatex(3.1,scl_layer3and,L3andEff.c_str());
tex22->SetLineWidth(2);
tex22->Draw();
TLatex *tex23 = new TLatex(5.1,scl_layer4or,L4orEff.c_str());
tex23->SetLineWidth(2);
tex23->Draw();
TLatex *tex24 = new TLatex(5.1,scl_layer4and,L4andEff.c_str());
tex24->SetLineWidth(2);
tex24->Draw();

/////////////////////////////////

c1->cd(4);

TLegend *l4;
l4 = new TLegend(0.75,0.75,0.99,0.99);
l4->AddEntry(h_count16, "Total Number of Events", "F");
l4->AddEntry(h_good_events, "\"Good Events\"", "F");
l4->AddEntry(h_count2, "Helicity | Up OR Down", "F");
l4->AddEntry(h_count3, "Helicity | Up AND Down", "F");
l4->AddEntry(h_count4, "Layer | Up OR Down", "F");
l4->AddEntry(h_count5, "Layer | Up AND Down", "F");

h_count16->SetLineColor(1);
h_count16->SetFillStyle(3004);
h_count16->SetFillColor(1);
h_count16->Draw();

h_good_events->SetLineColor(1);
h_good_events->SetFillStyle(3008);
h_good_events->SetFillColor(1);
h_good_events->Draw("same");


h_count17->SetLineColor(2);
h_count17->SetFillStyle(3004);
h_count17->SetFillColor(2);
h_count17->Draw("same");

h_count18->SetLineColor(7);
h_count18->SetFillStyle(3008);
h_count18->SetFillColor(7);
h_count18->Draw("same");

h_count19->SetLineColor(6);
h_count19->SetFillStyle(3004);
h_count19->SetFillColor(6);
h_count19->Draw("same");

h_count20->SetLineColor(3);
h_count20->SetFillStyle(3008);
h_count20->SetFillColor(3);
h_count20->Draw("same");

l4->Draw("same");


TLatex *tex25 = new TLatex(0.,scl_labels2,"Total Entries");
tex25->SetLineWidth(2);
tex25->Draw();
TLatex *tex26 = new TLatex(3.087098,scl_labels2,"Helicity");
tex26->SetLineWidth(2);
tex26->Draw();
TLatex *tex27 = new TLatex(5.195532,scl_labels2,"Layer");
tex27->SetLineWidth(2);
tex27->Draw();

TLatex *tex28 = new TLatex(3.1,scl_hel_ora,HelorEffa.c_str());
tex28->SetLineWidth(2);
tex28->Draw();
TLatex *tex29 = new TLatex(3.1,scl_hel_anda,HelandEffa.c_str());
tex29->SetLineWidth(2);
tex29->Draw();
TLatex *tex30 = new TLatex(5.1,scl_count_ora,SFTorEffa.c_str());
tex30->SetLineWidth(2);
tex30->Draw();
TLatex *tex31 = new TLatex(5.1,scl_count_anda,SFTandEffa.c_str());
tex31->SetLineWidth(2);
tex31->Draw();
////////////////

c1->cd(5);
hlayer1_2->Draw();

c1->cd(6);
hlayer3_4->Draw();

} 
