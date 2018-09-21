#ifndef __CINT__
#include <stdio.h>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <TROOT.h>
#include <typeinfo>
#include "TSystem.h"
#include "TFile.h"
#include "TProfile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TPolyLine.h"
#include "TLine.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TBranch.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TEllipse.h"
#include "TMarker.h"
#include "ANAPATH.h"
//#include "../../Macros/CommonParameters.h"
#include "ADC_TARGET_Pedestals.h"
#include "TDC_Windows.h"
#include "Cuts_and_Windows.h"
#include "Pedestals.h"
#include "intersect.cxx"
#include "GlobalParameters.h"
#include "G4DataRootApr19.h" //defines the variables in the .root files
#include "TOFsParameters.h"
#include "TargetParameters.h"
#include "G4Parameters.h"
#include "G4HistosApr.h"
#include "Plot_Simulation_Event_Display.h"
#endif

using namespace std;

int extract_tof1n_with_id(int ievt, TChain *fChain) {


  //Get data from tof1N
  fChain->GetEntry(ievt);
  Int_t n_tracks = 0;
  //Check for entries not -999
  for (int i = 0; i < 5; i++) {
    cout << tof1N[i] << " ";
    if (tof1N[i] != -999) { //normal
      n_tracks++;
    }
  }
  return n_tracks;
}

vector <int> extract_tof1n_with_id(int ievt, TChain *fChain, int ManyTracks) {


  //Get data from tof1N
  fChain->GetEntry(ievt);
  vector <int> tracks;

  //Check for entries not -999
  for (int i = 0; i < 5; i++) {
    //cout << tof1N[i] << " ";
    if (tof1N[i] != -999) { //normal
      tracks.push_back(tof1N[i]);
    }
  }
  return tracks;
}

void batch_extract_tof1N(int Run_Number) {
  char Name_finput[200];
  sprintf(Name_finput,"G4Run%d.root", Run_Number);
  cout << "File opened:  " << Name_finput << endl;

  TChain *fChain= new TChain("Kaon");
  fChain->Add(Name_finput);
  fChain->SetMakeClass(1);
  fChain->SetBranchAddress("tof1N[5]", tof1N);

  Int_t nentries = (Int_t)fChain->GetEntries();

  vector <int> more_than_one_track;
  vector <int> less_than_one_track;
  vector <int> one_track;
  vector <int> distinct_tracks;
  //cout << extract_tof1n_with_id(0, fChain);
  //return;
  for (int i = 0; i < nentries; ++i) {
    vector <int> num_track = extract_tof1n_with_id(i, fChain, 1);
    if (num_track.size() > 1){ //Check number of tracks
      more_than_one_track.push_back(i);

      //Find distinct tracks
      set<int> s( num_track.begin(), num_track.end() );

      //If more than 1 distinct tracks
      if (s.size() > 1){
        //for(vector<int>::const_iterator i = num_track.begin(); i != num_track.end(); ++i)
          //cout << *i << ' ';
        //cout << endl;
        distinct_tracks.push_back(i);
      }

    }
    else if (num_track.size() < 1) {
      less_than_one_track.push_back(i);
    }
    else {
      one_track.push_back(i);
    }
  }
  char output[100];
  //for (int i = 0; i < more_than_one_track.size(); ++i)
  sprintf(output,"RUN_%d_extract_tof1N.txt",Run_Number);

  ofstream output_file;
  output_file.open(output);
  output_file << "Distinctively More than 1 track" << endl;
  for(vector<int>::const_iterator i = distinct_tracks.begin(); i != distinct_tracks.end(); ++i) {
    output_file << *i << '\n';
}
  /*
  output_file << "More than 1 track" << endl;
  for(vector<int>::const_iterator i = more_than_one_track.begin(); i != more_than_one_track.end(); ++i) {
    output_file << *i << '\n';
}

  output_file << "Less than 1 track" << endl;
  for(vector<int>::const_iterator i = less_than_one_track.begin(); i != less_than_one_track.end(); ++i) {
    output_file << *i << '\n';
  }
  output_file << "1 track" << endl;
  for(vector<int>::const_iterator i = one_track.begin(); i != one_track.end(); ++i) {
    output_file << *i << '\n';
  }
  */
  output_file.close();
}
