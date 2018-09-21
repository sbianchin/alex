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
#include "TMarker.h"
#include "ANAPATH.h"
#include "Thresholds.h"
#include "CommonParameters.h"
#endif


void MWPC_counter(Int_t Run_Number=2008, Int_t thr=300){ 
	
	gStyle->Clear();
	TH1::AddDirectory(kFALSE);
	gStyle->SetOptStat(1111111);

	Int_t MwpcADC[512];	//Int_t MWPCADC[512];
	Int_t counter[512] = {0};
	
	char path_input[200];                   char file_mapping[200];
	sprintf(path_input,"%s",path_merged);	sprintf(file_mapping,"../Mapping");
	//sprintf(path_input,"/triumfcs/trshare/trek/E36/Data/April_2015/root");

	char Name_finput[200];
	sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);

	TChain *fChain= new TChain("Tree");		
	fChain->Add(Name_finput);		
	fChain->SetMakeClass(1);							
	
	fChain->SetBranchAddress("MWPCADC",MwpcADC);
	
	Int_t nentries = (Int_t)fChain->GetEntries();//		Int_t nentries_SFT = (Int_t)fChain->GetEntries();
	//if(nentries_TARGET==nentries_SFT) nentries = nentries_TARGET;
	cout << "  " << endl;
	cout << "****  Number of events: " << nentries << "  **** " <<endl;
	cout << "  " << endl;


	cout << "   " << endl;

	for(int i=0; i<nentries; i++){
			
		fChain->GetEntry(i);	//	fChain_SFT->GetEntry(i);

		if(i%10000==1)	cout<<"**** "<<i<<" events done"<<endl;

		for(int j=0; j<512; j++){
			if(MwpcADC[j] > thr){
				counter[j]++;			
			}		
		}	
		

	}

	for(int k=0; k<512; k++){
		cout << k << "    " << counter[k] << endl; 	
	}


}

