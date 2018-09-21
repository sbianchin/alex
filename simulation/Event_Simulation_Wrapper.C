
#include "Experimental_Event_Display_Simulation_1.C"

void Event_Simulation_Wrapper(int Run_Number) {
  //gROOT->LoadMacro("Experimental_Event_Display_Simulation_1.C");
  char tof1_track[200];
  sprintf(tof1_track,"/Users/alng/Documents/Datascience/TREK/offline/alex/simulation/RUN_%d_extract_tof1N.txt",
  Run_Number);

  ifstream fdat_read(tof1_track, ios::in);

  char Name_finput[200];
  sprintf(Name_finput,"G4Run%d.root", Run_Number);
  cout << "File opened:  " << Name_finput << endl;

  TChain *fChain= new TChain("Kaon");
  fChain->Add(Name_finput);
  fChain->SetMakeClass(1);

  Int_t nentries = (Int_t)fChain->GetEntries();

  while (fdat_read.good()) {
    int Read_ievt;
    fdat_read >> Read_ievt;
    Event_Display_5_1(Run_Number, Read_ievt);
  }
}
