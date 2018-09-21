#include <stdio.h>
#include <chrono>
// #include "Batch_Variables.h"
#include "Event_Display_5.5.C"


using namespace std;


Int_t StartChain(TChain &fChain, char Name_finput[200]) {

    fChain.Add(Name_finput);
    fChain.SetMakeClass(1);
    fChain.SetBranchAddress("VT48_TDC",tdc_vt48);

    fChain.SetBranchAddress("TDC_Trig",tdc_trigger); // tdc_trigger for TDC_diff calculation

    fChain.SetBranchAddress("ADC_High_TARGET",adc_high_target);    fChain.SetBranchAddress("ADC_High_SFT",adc_high_sft);
    fChain.SetBranchAddress("ADC_Low_TARGET",adc_low_target);      fChain.SetBranchAddress("ADC_Low_SFT",adc_low_sft);
    fChain.SetBranchAddress("TDC_LE_TARGET",tdc_le_target);        fChain.SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
    fChain.SetBranchAddress("TDC_TE_TARGET",tdc_te_target);        fChain.SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

    fChain.SetBranchAddress("ADC_TOF1U",ADC_tof1U);
    fChain.SetBranchAddress("ADC_TOF1D",ADC_tof1D);
    fChain.SetBranchAddress("TDC_TOF1U",TDC_tof1U);
    fChain.SetBranchAddress("TDC_TOF1D",TDC_tof1D);
    fChain.SetBranchAddress("ADC_TOF2AO",ADC_tof2AO);
    fChain.SetBranchAddress("ADC_TOF2AI",ADC_tof2AI);
    fChain.SetBranchAddress("ADC_TOF2BO",ADC_tof2BO);
    fChain.SetBranchAddress("ADC_TOF2BI",ADC_tof2BI);
    fChain.SetBranchAddress("TDC_TOF2AO",TDC_tof2AO);
    fChain.SetBranchAddress("TDC_TOF2AI",TDC_tof2AI);
    fChain.SetBranchAddress("TDC_TOF2BO",TDC_tof2BO);
    fChain.SetBranchAddress("TDC_TOF2BI",TDC_tof2BI);
    fChain.SetBranchAddress("MWPCADC",MwpcADC);
    fChain.SetBranchAddress("ADC_C2X_R",adc_c2x_r);
    fChain.SetBranchAddress("ADC_C2X_L",adc_c2x_l);
    fChain.SetBranchAddress("ADC_C2Y_R",adc_c2y_r);
    fChain.SetBranchAddress("ADC_C2Y_L",adc_c2y_l);
    fChain.SetBranchAddress("ADC_C3X_R",adc_c3x_r);
    fChain.SetBranchAddress("ADC_C3X_L",adc_c3x_l);
    fChain.SetBranchAddress("ADC_C3Y_R",adc_c3y_r);
    fChain.SetBranchAddress("ADC_C3Y_L",adc_c3y_l);
    fChain.SetBranchAddress("ADC_C4X_R",adc_c4x_r);
    fChain.SetBranchAddress("ADC_C4X_L",adc_c4x_l);
    fChain.SetBranchAddress("ADC_C4Y_R",adc_c4y_r);
    fChain.SetBranchAddress("ADC_C4Y_L",adc_c4y_l);
    fChain.SetBranchAddress("TDC_Ck", tdc_ck);
    fChain.SetBranchAddress("TDC_Cpi", tdc_cpi);
    fChain.SetBranchAddress("EvFlag", Event_flag);
    fChain.SetBranchAddress("EvTag", Event_tag);

    Int_t nentries = (Int_t)fChain.GetEntries();
    return nentries;

}


void Batch_5_5(int Run_Number, int min_ievt, int max_ievt, int flag) {
  gROOT->SetBatch(1);


  gStyle->Clear();
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(11);
  char batch_output[50];
  sprintf(batch_output, "RUN_%d_Batch_5.5.txt", Run_Number);
  ofstream fout;
  fout.open(batch_output);
  fout << "";
  fout.close();

  char path_input[200];
  sprintf(path_input,"%s",path_merged);
  char Name_finput[200];
  sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);
  cout << "File opened:  " << Name_finput << endl;

  Int_t nentries = StartChain(fChain, Name_finput);
  std::chrono::high_resolution_clock::time_point t1, t2;
  t1 = std::chrono::high_resolution_clock::now();
  if (flag == 0) {
    int i = min_ievt;
    while (i < max_ievt) {
      if (i % 10000 == 9999) {
        cout << i << "events done!" << endl;
      }
      Event_Display_5_5(Run_Number, i, 0, 0, 0, 0, 1, 1);
      ++i;
    }
  }
  else if (flag == 1) {
    int i = min_ievt;
    while (i < max_ievt) {
      Event_Display_5_5(Run_Number, i, 0, 0, 0, 1, 0, 1);
      ++i;
    }
  }
  else if (flag == 2) {
    int i = 0;
    while (i < nentries) {
      if (i % 10000 == 9999) {
        cout << i << "events done!" << endl;
      }
      Event_Display_5_5(Run_Number, i, 0, 0, 0, 0, 1, 1);
      ++i;
    }
  }
  else if (flag == 3){
    char events[200];
    sprintf(events,"/Users/alng/Documents/Datascience/TREK/offline/alex/RUN_%d_read_file.txt",
    Run_Number);
    ifstream fdat_read(events, ios::in);
    int Read_ievt;
    while (fdat_read.good()) {

      fdat_read >> Read_ievt;
      cout << Read_ievt << endl;
      Event_Display_5_5(Run_Number, Read_ievt, 0, 0, 0, 0, 1, 1);
    }
  }
  else cout << "Error: Flag needs to be 0, 1, 2, or 3.";
  t2 = std::chrono::high_resolution_clock::now();
  cout << "Time taken for all events "
                << std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count()
                << " seconds\n";
}
