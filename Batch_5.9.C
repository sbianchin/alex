#include <stdio.h>
#include <chrono>
#include "Batch_Variables.h"
#include "Event_Display_5.9.C"


using namespace std;

void Batch_5_9(int Run_Number, int Read_File = 0) {
  if (Read_File == 0) {
    // gROOT->Clear();
    gROOT->SetBatch(1);


    gStyle->Clear();
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(11);

    char path_input[200];
    sprintf(path_input,"%s",path_merged);

    char Name_finput[200];
    sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);
    cout << "File opened:  " << Name_finput << endl;



    TChain *fChain= new TChain("Tree");
    fChain->Add(Name_finput);
    fChain->SetMakeClass(1);
    Int_t nentries = (Int_t)fChain->GetEntries();
    // Int_t nentries = StartChain(fChain, Name_finput);

    // for (int i = 0; i < nentries; i++) {
    int i = 4750;
    // std::chrono::high_resolution_clock::time_point t1, t2;

    while (i < nentries) {
      // if (i % 400 == 0) {
      //   t1 = std::chrono::high_resolution_clock::now();
      // }
      // if (i % 400 == 399) {
      //   t2 = std::chrono::high_resolution_clock::now();
      //   cout << "Time taken for 100 events "
      //               << std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count()
      //               << " microseconds\n";
      // }
      cout << i << endl;
      // Event_Display_5_8(Run_Number, i, 0, 0, 0, 0);
      Event_Display_5_9(Run_Number, i, 0, 0, 0, 0);
      // cout << dist_to_k_stop;
      ++i;
    }
  }
  else {
    char events[200];
    sprintf(events,"/Users/alng/Documents/Datascience/TREK/offline/alex/RUN_%d_read_file.txt",
    Run_Number);
    ifstream fdat_read(events, ios::in);
    int Read_ievt;
    // while (!feof (file))
    // {
    //   printf ("%d ", i);
    //   fscanf (file, "%d", &i);
    // }
    while (fdat_read.good()) {

      fdat_read >> Read_ievt;
      cout << Read_ievt << endl;
      Event_Display_5_9(Run_Number, Read_ievt, 0, 0, 0, 0);
    }
  }
}
