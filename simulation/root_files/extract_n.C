#include <stdio.h>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <string.h>
#include "TSystem.h"
#include "TFile.h"
#include "TProfile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
using namespace std;

void batch_extract_n(int Run_Number) {
	char Name_finput[200];
	sprintf(Name_finput,"G4Run%d.root", Run_Number);
	cout << "File opened:  " << Name_finput << endl;
	int n;
	TChain *fChain= new TChain("Kaon");
	fChain->Add(Name_finput);
	fChain->SetMakeClass(1);
	fChain->SetBranchAddress("n", &n);
	Int_t nentries = (Int_t)fChain->GetEntries();

	char output[100];
	//for (int i = 0; i < more_than_one_track.size(); ++i)
	sprintf(output,"RUN_%d_extract_n.csv",Run_Number);
	ofstream output_file;
	output_file.open(output);
	output_file << "Entry Number, Event Number," << endl;
	for (int i = 0; i < nentries; ++i) {

		fChain->GetEntry(i);
		output_file << i << ", " << n << "," << endl;
	}

}
