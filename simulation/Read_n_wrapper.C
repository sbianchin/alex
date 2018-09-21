#include <map>
using namespace std;

void Event_Simulation_Wrapper(int Run_Number) {
  //gROOT->LoadMacro("Experimental_Event_Display_Simulation_1.C");
  char event_track[200];
  sprintf(event_track,"/Users/alng/Documents/Datascience/TREK/offline/alex/simulation/RUN_%d_extract_n.txt",
  Run_Number);

  // ifstream fdat_read(tof1_track, ios::in);

  char Name_finput[200];
  sprintf(Name_finput,"G4Run%d.root", Run_Number);
  cout << "File opened:  " << Name_finput << endl;
  int n;
  TChain *fChain= new TChain("Kaon");
  fChain->Add(Name_finput);
  fChain->SetMakeClass(1);
  fChain->SetBranchAddress("n", &n);

  Int_t nentries = (Int_t)fChain->GetEntries();
  map <int, int> event_to_n;
  ofstream n_to_event;
  n_to_event.open(event_track);
  // for(auto& kv : event_to_n) {
  //   event_track << kv.second;
  // }
  n_to_event << "map<int, int> root_entry_to_tree_entry = {" << endl;
  for (int i = 0; i < nentries; i++) {
    fChain->GetEntry(i);
    // event_to_n[i] = n;
    n_to_event << "{"<<n << ", "<< i << "},"<< endl;
  }
  n_to_event << " }" << endl;
  n_to_event.close();
}
