#include <map>
#include "Event_Display_5.3.C"
using namespace std;

void Batch_5_3(int Run_Number, int Read_File = 0) {
  if (Read_File == 0) {
    gROOT->Clear();

    gStyle->Clear();
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(11);

    char path_input[200];
    sprintf(path_input,"%s",path_merged);

    char Name_finput[200];
    sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);
    cout << "File opened:  " << Name_finput << endl;

    TChain *fChain= new TChain("Kaon");
    fChain->Add(Name_finput);
    fChain->SetMakeClass(1);
    Int_t nentries = (Int_t)fChain->GetEntries();

    // for (int i = 0; i < nentries; i++) {
    for (int i = 5000; i < 10000; i++) {
      Event_Display_5_3(Run_Number, i, 0, 0, 0, 0);
    }
  }
  else {
    char events[200];
    sprintf(events,"/Users/alng/Documents/Datascience/TREK/offline/alex/RUN_%d_events.txt",
    Run_Number);
    ifstream fdat_read(events, ios::in);

    while (fdat_read.good()) {
      int Read_ievt;
      fdat_read >> Read_ievt;
      Event_Display_5_3(Run_Number, Read_ievt, 0, 0, 0, 0);
    }
  }
}
