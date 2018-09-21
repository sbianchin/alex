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
#include "CommonParameters.h"
#include "ADC_TARGET_Pedestals.h"
#include "TDC_Windows.h"
#include "Cuts_and_Windows.h"
#include "MWPC_Thr2.h"
#include "Pedestals.h"
#include "intersect.cxx"
#include "C2_Strip_transform.h"
#include "Channel_to_Strip.h"
#include "SFT_functions.h"
#include "Plot_Event_Display.C"
#endif

using namespace std;
int kaon_fiber(double x_bar, double y_bar);
/*
void FitChiSquare(int vec_size, vector<double> vec_xx, vector<double> vec_yy) {
  double xiyi[vec_size];
  for (int i = 0; i < vec_size; ++i) {
    xiyi[i] = vec_xx[i] * vec_yy[i];
  }
  double sum_xiyi = 0;
  for(auto& num : xiyi) sum_xiyi += num;
  double sum_xi = accumulate(vec_xx.begin(), vec_xx.end(), 0);
  double sum_yi = accumulate(vec_yy.begin(), vec_yy.end(), 0);
  double xixi[vec_size];
  for (int i = 0; i < vec_size; ++i) {
    xixi[i] = vec_xx[i] * vec_xx[i];
  }
  double sum_xixi = 0;
  for(auto& num : xixi) sum_xiyi += num;
  double m = (vec_size * sum_xiyi - sum_xi * sum_yi) / (vec_size * sum_xixi - sum_xi * sum_xi);
  double b = (sum_xixi * sum_yi - sum_xiyi * sum_xi) / (vec_size * sum_xixi - sum_xi * sum_xi);

  double didi[vec_size];
  for (int i = 0; i < vec_size; ++i) {
    didi[i] = (vec_yy[i] - (m * vec_xx[i] + b)) * (vec_yy[i] - (m * vec_xx[i] + b));
  }
  double sum_didi = 0;
  for(auto& num : didi) sum_didi += num;
  dm_sq = (sum_didi / (vec_size * sum_xixi - sum_xi * sum_xi)) * (vec_size / (vec_size - 2));
  dm = pow(dm_sq, 0.5);

  db_sq = (sum_didi * sum_xixi) / (vec_size * sum_xixi - sum_xi * sum_xi) / (vec_size - 2);
  db = pow(db_sq, 0.5);


}
*/
/*
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

    Int_t nentries = (Int_t)fChain.GetEntries();
    return nentries;

}
*/

void Event_Display_5_9(Int_t Run_Number=5, Int_t ievt=0, int walk=0, Int_t enable_cout=0, Int_t display = 0, Int_t graphics = 1){

  int Switch=1; // Displays hit with no HG, but LG (0 = OFF ; 1 = ON)
  int Rotate=1; // When TOF1 is 12 or 6, or when TOF2 is 6 or 12, rotate by -90 deg to fit a horizontal line (0 = OFF ; 1 = ON)
  int MWPC_switch = graphics; // Display wire chamber hits (0 = OFF ; 1 = ON)

  ifstream blacklist;
  blacklist.open("Event_Blacklist_Test.txt");

  if(enable_cout!=0 && enable_cout!=1 && enable_cout!=2 && enable_cout!=9){
    cout << "  " << endl;
    cout << "Flag Error !" << endl;
    cout << "  " << endl;
    return;
  }

  int T_limit = 3;

  // gROOT->Clear();

  gStyle->Clear();
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(11);

  char path_input[200];
  sprintf(path_input,"%s",path_merged);

  char Name_finput[200];
  sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);

  char footer[100];
  sprintf(footer,"Event_Display.C  --  Run %d ; Event %d",Run_Number,ievt);

  char Version[100] = "Version 5.3";
  char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!

  TChain fChain("Tree");
  fChain.Add(Name_finput);
  fChain.SetMakeClass(1);

  Double_t flag_size_TARGET=1.35;
  Double_t flag_size_SFT=1.3;
  Double_t flag_size_palette=1.6;

  Int_t HG_TARGET_ADC_Thr[256] = {0};
  Int_t LG_TARGET_ADC_Thr[256] = {0};
  Int_t HG_SFT_ADC_Thr[128] = {0};
  Int_t LG_SFT_ADC_Thr[128] = {0};

  for(int i=0; i<256; i++)  HG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Ped_HG[i]) + TARGET_ADC_Thr_HG_Offset;
  for(int i=0; i<256; i++)  LG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Ped_LG[i]) + TARGET_ADC_Thr_LG_Offset;
  for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;
  for(int i=0; i<128; i++)  LG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_LG[i]) + SFT_ADC_Thr_LG_Offset;

  Int_t TDC_min_SFT = SFT_TDC_min[0];
  Int_t TDC_max_SFT = SFT_TDC_max[0];

  bool TARGET_High_has_neighbours[256] = {false}; // Array of High gain target hits which have no neighbouring targets hit

  int n_hit = 2;   // ## Number of hits required in the TARGET

  Int_t ADC_High_corr_max=0;

  Int_t ADC_cut_TARGET2 = 850;

  char ADC_cut[100];    sprintf(ADC_cut,"(ADC >= %d)",HG_SFT_ADC_Thr[0]);

  Int_t adc_high_target[256];       Int_t ADC_High_TARGET[256];
  Int_t adc_low_target[256];        Int_t ADC_Low_TARGET[256];
  Int_t tdc_le_target[256][16];     Int_t TDC_LE_TARGET[256];
  Int_t tdc_te_target[256][16];
  Int_t TDC_LE_TARGET_corrected[256][6]={0};
  char TDC_LE_TARGET_corr[256][6][20];

  Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];        Double_t ADC_High_SFT_corr[128];
  Int_t adc_low_sft[128];           Int_t ADC_Low_SFT[128];
  Int_t tdc_le_sft[128][16];        Int_t TDC_LE_SFT[128];
  Int_t tdc_te_sft[128][16];

  Int_t ADC_TOF1[24];

  Int_t ADC_tof1U[12];              Int_t ADC_TOF1U[12];
  Int_t ADC_tof1D[12];              Int_t ADC_TOF1D[12];

  Int_t TDC_tof1U[12];              Int_t TDC_TOF1U[12];
  Int_t TDC_tof1D[12];              Int_t TDC_TOF1D[12];

  Int_t ADC_tof2AO[12];             Int_t ADC_TOF2AO[12];
  Int_t ADC_tof2BO[12];             Int_t ADC_TOF2BO[12];
  Int_t ADC_tof2AI[12];             Int_t ADC_TOF2AI[12];
  Int_t ADC_tof2BI[12];             Int_t ADC_TOF2BI[12];

  Int_t TDC_tof2AO[12];             Int_t TDC_TOF2AO[12];
  Int_t TDC_tof2BO[12];             Int_t TDC_TOF2BO[12];
  Int_t TDC_tof2AI[12];             Int_t TDC_TOF2AI[12];
  Int_t TDC_tof2BI[12];             Int_t TDC_TOF2BI[12];

  Int_t MwpcADC[512];

  Int_t adc_c2x_r[56];              double ADC_C2X_R[56];
  Int_t adc_c2x_l[56];              double ADC_C2X_L[56];
  Int_t adc_c2y_r[16];              double ADC_C2Y_R[16];
  Int_t adc_c2y_l[16];              double ADC_C2Y_L[16];

  Int_t adc_c3x_r[64];              double ADC_C3X_R[64];
  Int_t adc_c3x_l[64];              double ADC_C3X_L[64];
  Int_t adc_c3y_r[16];              double ADC_C3Y_R[16];
  Int_t adc_c3y_l[16];              double ADC_C3Y_L[16];

  Int_t adc_c4x_r[72];              double ADC_C4X_R[72];
  Int_t adc_c4x_l[72];              double ADC_C4X_L[72];
  Int_t adc_c4y_r[16];              double ADC_C4Y_R[16];
  Int_t adc_c4y_l[16];              double ADC_C4Y_L[16];

  Int_t tdc_ck[14][16];
  Int_t tdc_cpi[14][16];

  double TDC_diff = -1;
  int tdc_trigger[2][16];
  double tdc_ck_corr = 0.;


  Int_t has_TDC_SFT_hit[128] = {0};

  Bool_t Event_flag[40] = {false};

  float R_TARGET = 29.0;
  float R_TOF1 = 47.1;
  float R_SFT_L1 = 40.0;

  char run_string[100];            char event_string[100];

  char h_ADC_title[200];
  sprintf(h_ADC_title,"(ADC offset = %d) | (%d < TDC < %d)  --  Run %d (Event %d)",SFT_ADC_Thr_HG_Offset, TDC_min_SFT,TDC_max_SFT,Run_Number,ievt);

  double Z_TOF1[12] = {999.99};

  const double Offset_TOF1[12] = {180., 410., 0., 0., 310., 280., -30., 70., 0., 0., 390., 430.};
  const double Slope_correction_TOF1[12] = {0.70, 0.78, 1., 1., 0.70, 0.60, 0.86, 0.77, 1., 1., 0.67, 0.78};

  char file_mapping[200];
  sprintf(file_mapping,"../Mapping");

  char par_finput[200];
  sprintf(par_finput,"%s/%s",file_mapping,source_mapping);

  char par_finput2[200];
  sprintf(par_finput2,"%s/MWPC_map2.txt",file_mapping);
  Int_t par_temp[2][128];
  ifstream fdat(par_finput,ios::in);
  for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
  fdat.close();

  char par_temp2[512][50];
  ifstream fdat2(par_finput2,ios::in);
  for(Int_t ii=0; ii<512; ii++) fdat2 >> par_temp2[ii];
  fdat2.close();

  vector <double> vec_xx_lepton;            vec_xx_lepton.clear();
  vector <double> vec_yy_lepton;            vec_yy_lepton.clear();
  vector <double> vec_ex_lepton;            vec_ex_lepton.clear();
  vector <double> vec_ey_lepton;            vec_ey_lepton.clear();

  vector <double> vec_xx_lepton_test;       vec_xx_lepton_test.clear();
  vector <double> vec_yy_lepton_test;       vec_yy_lepton_test.clear();

  vector <double> vec_xx_lepton_rotate;     vec_xx_lepton_rotate.clear();
  vector <double> vec_yy_lepton_rotate;     vec_yy_lepton_rotate.clear();
  vector <double> vec_ex_lepton_rotate;     vec_ex_lepton_rotate.clear();
  vector <double> vec_ey_lepton_rotate;     vec_ey_lepton_rotate.clear();

  vector <int> lepton_bars;             lepton_bars.clear();
  vector <int> lepton_bars_rotate;      lepton_bars_rotate.clear();
  vector <int> vec_kaon_bars;           vec_kaon_bars.clear();

  vector <double> vec_lepton_size;          vec_lepton_size.clear();
  vector <double> vec_lepton_rotate_size;   vec_lepton_rotate_size.clear();
  vector <int> vec_bar;                     vec_bar.clear();
  vector <int> vec_bar_rotate;              vec_bar_rotate.clear();
  vector <double> vec_yprime;               vec_yprime.clear();
  vector <double> vec_yprime_rotate;        vec_yprime_rotate.clear();
  vector <double> vec_Dy;                   vec_Dy.clear();
  vector <double> vec_Dy_rotate;            vec_Dy_rotate.clear();

  double sumS = 999.;
  double sumS_rotate = 999.;

  vector <double> vec_xx_kaon;            vec_xx_kaon.clear();
  vector <double> vec_yy_kaon;            vec_yy_kaon.clear();
  vector <double> vec_ex_kaon;            vec_ex_kaon.clear();
  vector <double> vec_ey_kaon;            vec_ey_kaon.clear();

  vector <double> vec_Ck;            vec_Ck.clear();
  vector <double> vec_Cpi;           vec_Cpi.clear();

  vector <double> vec_xx_TOF1_Marker;     vec_xx_TOF1_Marker.clear();
  vector <double> vec_yy_TOF1_Marker;     vec_yy_TOF1_Marker.clear();

  vector <double> vec_xx_TOF1;      vec_xx_TOF1.clear();
  vector <double> vec_yy_TOF1;      vec_yy_TOF1.clear();

  vector <double> vec_xx_TOF1_closest;    vec_xx_TOF1_closest.clear();
  vector <double> vec_yy_TOF1_closest;    vec_yy_TOF1_closest.clear();

  vector <double> vec_xx_int_TDC_Gap_Fibers;   vec_xx_int_TDC_Gap_Fibers.clear();
  vector <double> vec_yy_int_TDC_Gap_Fibers;   vec_yy_int_TDC_Gap_Fibers.clear();

  vector <double> vec_xx_int_TDC_Gap_Fibers_SFT;   vec_xx_int_TDC_Gap_Fibers_SFT.clear();
  vector <double> vec_yy_int_TDC_Gap_Fibers_SFT;   vec_yy_int_TDC_Gap_Fibers_SFT.clear();

  vector <double> vec_xx_int_TDC_TARGET;      vec_xx_int_TDC_TARGET.clear();
  vector <double> vec_yy_int_TDC_TARGET;      vec_yy_int_TDC_TARGET.clear();

  vector <double> vec_xx_kaon_stop;   vec_xx_kaon_stop.clear();
  vector <double> vec_yy_kaon_stop;   vec_yy_kaon_stop.clear();

  Int_t TDC_ck_selected[14] = {0};
  Int_t TDC_cpi_selected[14] = {0};

  Int_t TDC_ck_sum = 0;       double TDC_ck_avg = 0.;     double TDC_ck_sigma = 0.;
  Int_t TDC_cpi_sum = 0;      double TDC_cpi_avg = 0.;    double TDC_cpi_sigma = 0.;

  double TDC_ck_sigma2 = 0.;
  double TDC_cpi_sigma2 = 0.;

  Int_t TDC_ck_sum2 = 0;    double TDC_ck_avg2=0.;    int TDC_ck_counter = 0;
  Int_t TDC_cpi_sum2 = 0;   double TDC_cpi_avg2=0.;   int TDC_cpi_counter = 0;

  // B0 counter
  int tdc_vt48[256][16];
  fChain.SetBranchAddress("VT48_TDC",tdc_vt48);
  // vector<int> vec_tdc_b0_6;  vec_tdc_b0_6.clear();
  // vector<int> vec_tdc_b0_7;  vec_tdc_b0_7.clear();

  vector<int> vec_TARGET_bar_selected;  vec_TARGET_bar_selected.clear();

  double Gap[12][3][2] = {{{0}}};


  /////////////////////   Dave's Time Walk Correction File  ////////////////////////
  double par_in[256][3] = {0.};
  double par_err[356][3] = {0.};
  Int_t ADCmax = 3450;
  double Yfi = 0.;
  double Ani = 0.;
  double Yni = 0.;
  double Tpi = 0.;
  bool NoFile = false;

  char ParsTarg1[100];
  sprintf(ParsTarg1,"TimeWalk%d.dat",Run_Number);

  if(!ifstream(ParsTarg1)) NoFile = true;

  sprintf(run_string,"Run %d ; Event %d",Run_Number,ievt);
  sprintf(event_string,"Run %d ; Event %d",Run_Number,ievt);


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

  // Int_t nentries = (Int_t)fChain->GetEntries();

  bool Good_Event=false;

  for(int ivt=ievt; ivt<ievt+1; ivt++){
    auto t1 = std::chrono::high_resolution_clock::now();

    fChain.GetEntry(ivt);

    auto t2 = std::chrono::high_resolution_clock::now();
    cout << "Time taken "
                << std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count()
                << " microseconds\n";

    for (int j_TARGET=0; j_TARGET<256; j_TARGET++){
      ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-HG_TARGET_ADC_Thr[j_TARGET];
      ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-LG_TARGET_ADC_Thr[j_TARGET];
      TDC_LE_TARGET[j_TARGET]=tdc_le_target[j_TARGET][0];
    }


    ////////////////////////  Dave's Time Walk Correction ///////////////////////

    if(walk==0){
      ifstream parTARGdat(ParsTarg1,ios::in);
      Int_t ij = 0;
      Int_t ik = 0;

      // Read in parameters and their errors. (errors not used)
      for(Int_t ii = 0; ii < nBars; ii++){
        parTARGdat >> ij >> ik >> par_in[ii][0] >> par_err[ii][0];
        parTARGdat >> ij >> ik >> par_in[ii][1] >> par_err[ii][1];
        parTARGdat >> ij >> ik >> par_in[ii][2] >> par_err[ii][2];
      }

      for(Int_t ii = 0; ii<256; ii++){

        Yfi = par_in[ii][0] - par_in[ii][1]/sqrt(ADCmax - par_in[ii][2]);
        Ani = adc_high_target[ii]-TARGET_ADC_Ped_HG[ii]; //SEB

        if((Ani >= TARGET_ADC_Thr_HG_Offset) && (Ani < ADCmax)){ // SEB
          Yni = par_in[ii][0] - par_in[ii][1]/sqrt(Ani - par_in[ii][2]);
          Tpi = Yfi - Yni;
          for(int jj=0; jj<16; jj++){
            if(tdc_le_target[ii][jj]>0) tdc_le_target[ii][jj] = tdc_le_target[ii][jj] + Tpi;
          }
        }
      }
    }

    /////////////////////////////////////////////////////////////////////////////

    for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
      ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-HG_SFT_ADC_Thr[j_SFT];
      TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
      ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-LG_SFT_ADC_Thr[j_SFT];
    }

    for(int i=0; i<12; i++){
      ADC_TOF1[i] = ADC_tof1U[i]-TOF1U_ADC_Thr[i];
      ADC_TOF1[i+12] = ADC_tof1D[i]-TOF1D_ADC_Thr[i];
    }

    for (Int_t j_TOF2=0; j_TOF2<12; j_TOF2++) {
      ADC_TOF2AO[j_TOF2] = ADC_tof2AO[j_TOF2]-TOF2AO_ADC_Thr[j_TOF2];
      ADC_TOF2BO[j_TOF2] = ADC_tof2BO[j_TOF2]-TOF2AI_ADC_Thr[j_TOF2];
      ADC_TOF2AI[j_TOF2] = ADC_tof2AI[j_TOF2]-TOF2BO_ADC_Thr[j_TOF2];
      ADC_TOF2BI[j_TOF2] = ADC_tof2BI[j_TOF2]-TOF2BI_ADC_Thr[j_TOF2];
   }

    for(Int_t j_C2=0; j_C2<56; j_C2++){
      ADC_C2X_R[j_C2] = double(adc_c2x_r[j_C2])-(ADC_C2X_Thr_R[j_C2]+MWPC_Thr_Offset_C2X);
      ADC_C2X_L[j_C2] = double(adc_c2x_l[j_C2])-(ADC_C2X_Thr_L[j_C2]+MWPC_Thr_Offset_C2X);
    }

    for(Int_t j_C3=0; j_C3<64; j_C3++){
      ADC_C3X_R[j_C3] = double(adc_c3x_r[j_C3])-(ADC_C3X_Thr_R[j_C3]+MWPC_Thr_Offset_C3X);
      ADC_C3X_L[j_C3] = double(adc_c3x_l[j_C3])-(ADC_C3X_Thr_L[j_C3]+MWPC_Thr_Offset_C3X);
    }

    for(Int_t j_C4=0; j_C4<72; j_C4++){
      ADC_C4X_R[j_C4] = double(adc_c4x_r[j_C4])-(ADC_C4X_Thr_R[j_C4]+MWPC_Thr_Offset_C4X);
      ADC_C4X_L[j_C4] = double(adc_c4x_l[j_C4])-(ADC_C4X_Thr_L[j_C4]+MWPC_Thr_Offset_C4X);
    }

    for(Int_t j_CY=0; j_CY<16; j_CY++){
      ADC_C2Y_R[j_CY] = double(adc_c2y_r[j_CY])-(ADC_C2Y_Thr_R[j_CY]+MWPC_Thr_Offset_C2Y);
      ADC_C2Y_L[j_CY] = double(adc_c2y_l[j_CY])-(ADC_C2Y_Thr_L[j_CY]+MWPC_Thr_Offset_C2Y);
      ADC_C3Y_R[j_CY] = double(adc_c3y_r[j_CY])-(ADC_C3Y_Thr_R[j_CY]+MWPC_Thr_Offset_C3Y);
      ADC_C3Y_L[j_CY] = double(adc_c3y_l[j_CY])-(ADC_C3Y_Thr_L[j_CY]+MWPC_Thr_Offset_C3Y);
      ADC_C4Y_R[j_CY] = double(adc_c4y_r[j_CY])-(ADC_C4Y_Thr_R[j_CY]+MWPC_Thr_Offset_C4Y);
      ADC_C4Y_L[j_CY] = double(adc_c4y_l[j_CY])-(ADC_C4Y_Thr_L[j_CY]+MWPC_Thr_Offset_C4Y);
    }

    for (Int_t j_TDCTOF=0; j_TDCTOF<12; j_TDCTOF++) {
      TDC_TOF1U[j_TDCTOF] = TDC_tof1U[j_TDCTOF];
      TDC_TOF1D[j_TDCTOF] = TDC_tof1D[j_TDCTOF];
      TDC_TOF2AO[j_TDCTOF] = TDC_tof2AO[j_TDCTOF];
      TDC_TOF2BO[j_TDCTOF] = TDC_tof2BO[j_TDCTOF];
      TDC_TOF2AI[j_TDCTOF] = TDC_tof2AI[j_TDCTOF];
      TDC_TOF2BI[j_TDCTOF] = TDC_tof2BI[j_TDCTOF];
    }

    //********* GOOD TARGET EVENTS
    bool Good_TARGET_Event = false;
    int count_TARGET_evts = 0;
    for(int i=0; i<256; i++)
    {
       if((ADC_High_TARGET[i]>0 && tdc_le_target[i][0]>=TARGET_TDC_min[i] && tdc_le_target[i][0]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]>0 && tdc_le_target[i][1]>=TARGET_TDC_min[i] && tdc_le_target[i][1]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]>0 && tdc_le_target[i][2]>=TARGET_TDC_min[i] && tdc_le_target[i][2]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]>0 && tdc_le_target[i][3]>=TARGET_TDC_min[i] && tdc_le_target[i][3]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]<0  && ADC_Low_TARGET[i]>0)) count_TARGET_evts++;
    }

    if(count_TARGET_evts >= n_hit) Good_TARGET_Event = true;

    if(enable_cout == 9) Good_TARGET_Event = true;

    //********* GOOD TOF EVENTS
    bool Good_TOF1_ADC[12]={false};   bool Good_TOF2_ADC[12]={false};
    bool Good_TOF1_TDC[12]={false};   bool Good_TOF2_TDC[12]={false};
    bool Good_TOF1[12]={false};       bool Good_TOF2[12]={false};
    bool Good_TOFs[12]={false};
    bool Good_TOF_Event = false;

    bool Event_On_Blacklist = false;
    string Current_Event;
    int current_event;

    if(blacklist.fail()){
      cout << "Error: Could not read blacklist file." << endl;
    }
    else{
        while(getline(blacklist,Current_Event)){
      sscanf(Current_Event.c_str(), "%d", &current_event);
        if(current_event == ievt){
          Event_On_Blacklist = true;
          break;
        }
      }
    }

    blacklist.close();

    for(int i =0; i<12; i++){
      ADC_TOF1U[i] = ADC_TOF1[i];
      ADC_TOF1D[i] = ADC_TOF1[i+12];
    }

    for(int i=0; i<12; i++){

      if(ADC_TOF1U[i]>=0 || ADC_TOF1D[i]>=0)  Good_TOF1_ADC[i] = true;

      if((TDC_TOF1U[i]>=TOF1U_TDC_min[i] && TDC_TOF1U[i]<=TOF1U_TDC_max[i]) ||
         (TDC_TOF1D[i]>=TOF1D_TDC_min[i] && TDC_TOF1D[i]<=TOF1D_TDC_max[i]))  Good_TOF1_TDC[i] = true;

      if(Good_TOF1_TDC[i] || Good_TOF1_ADC[i]) Good_TOF1[i] = true;

      if((ADC_TOF2AO[i]>=0 && ADC_TOF2AI[i]>=0) || (ADC_TOF2BO[i]>=0 && ADC_TOF2BI[i]>=0))  Good_TOF2_ADC[i] = true;

      if(((TDC_TOF2AO[i]>=TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<=TOF2AO_TDC_max[i])  &&
          (TDC_TOF2AI[i]>=TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<=TOF2AI_TDC_max[i])) ||
         ((TDC_TOF2BO[i]>=TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<=TOF2BO_TDC_max[i])  &&
          (TDC_TOF2BI[i]>=TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<=TOF2BI_TDC_max[i])))  Good_TOF2_TDC[i] = true;

      if(Good_TOF2_ADC[i] || Good_TOF2_TDC[i]) Good_TOF2[i] = true;
    }

    for(int k=0; k<12; k++){
      if(k!=0 && k!=11){
        if((Good_TOF2[k] && Good_TOF1[k-1]) || (Good_TOF2[k] && Good_TOF1[k]) || (Good_TOF2[k] && Good_TOF1[k+1])) Good_TOFs[k] = true;
      }
    }

    if((Good_TOF2[0] && Good_TOF1[11]) || (Good_TOF2[0] && Good_TOF1[0]) || (Good_TOF2[0] && Good_TOF1[1]))  Good_TOFs[0] = true;

    if((Good_TOF2[11] && Good_TOF1[10]) || (Good_TOF2[11] && Good_TOF1[11]) || (Good_TOF2[11] && Good_TOF1[0]))  Good_TOFs[11] = true;

    for(int kk=0; kk<12; kk++){
      if(Good_TOFs[kk]) Good_TOF_Event = true;
    }

    if(enable_cout == 9) Good_TOF_Event = true;


    //********* GOOD MWPC EVENTS
    int count_C2X = 0;    int count_C2Y = 0;
    int count_C3X = 0;    int count_C3Y = 0;
    int count_C4X = 0;    int count_C4Y = 0;
    bool Good_MWPC_Event = false;

    for(int i=0; i<56; i++){
      if(ADC_C2X_R[i]>0. || ADC_C2X_L[i]>0.) count_C2X++;
    }

    for(int ii=0; ii<16; ii++){
      if(ADC_C2Y_R[ii]>0. || ADC_C2Y_L[ii]>0.) count_C2Y++;
    }

    for(int j=0; j<64; j++){
      if(ADC_C3X_R[j]>0. || ADC_C3X_L[j]>0.) count_C3X++;
    }

    for(int jj=0; jj<16; jj++){
      if(ADC_C3Y_R[jj]>0. || ADC_C3Y_L[jj]>0.) count_C3Y++;
    }

    for(int k=0; k<72; k++){
      if(ADC_C4X_R[k]>0. || ADC_C4X_L[k]>0.) count_C4X++;
    }

    for(int kk=0; kk<16; kk++){
      if(ADC_C4Y_R[kk]>0. || ADC_C4Y_L[kk]>0.) count_C4Y++;
    }

    if(count_C2X>0 && count_C2Y>0 && count_C3X>0 && count_C3Y>0 && count_C4X>0 && count_C4Y>0)  Good_MWPC_Event = true;

    if(enable_cout == 9) Good_MWPC_Event = true;

    if(Good_TARGET_Event && Good_TOF_Event && Good_MWPC_Event && !Event_On_Blacklist)  Good_Event = true;

    //if(Good_Event) cout << "Event: "<< ievt << "   --  GOOD EVENT!" << endl;

    if(Good_TARGET_Event && Good_TOF_Event && !Good_MWPC_Event){
      cout << "Event: "<< ievt << "   --  NO MWPC!" << endl;
      cout << " " << endl;
      cout << ">>>>>  Please, choose another event" << endl;
      cout << " " << endl;
      break;
    }

    if(Good_TARGET_Event && !Good_TOF_Event && Good_MWPC_Event){
      cout << "Event: "<< ievt << "   --  NO TOF!" << endl;
      cout << " " << endl;
      cout << ">>>>>  Please, choose another event" << endl;
      cout << " " << endl;
      break;
    }

    if(!Good_TARGET_Event && Good_TOF_Event && Good_MWPC_Event){
      cout << "Event: "<< ievt << "   --  NO TARGET!" << endl;
      cout << " " << endl;
      cout << ">>>>>  Please, choose another event" << endl;
      cout << " " << endl;
      break;
    }

    if(Good_TARGET_Event && !Good_TOF_Event && !Good_MWPC_Event){
      cout << "Event: "<< ievt << "   --  NO TOF and NO MWPC!" << endl;
      cout << " " << endl;
      cout << ">>>>>  Please, choose another event" << endl;
      cout << " " << endl;
      break;
    }

    if(!Good_TARGET_Event && Good_TOF_Event && !Good_MWPC_Event){
      cout << "Event: "<< ievt << "   --  NO TARGET and NO MWPC!" << endl;
      cout << " " << endl;
      cout << ">>>>>  Please, choose another event" << endl;
      cout << " " << endl;
      break;
    }

    if(!Good_TARGET_Event && !Good_TOF_Event && Good_MWPC_Event){
      cout << "Event: "<< ievt << "   --  NO TARGET and NO TOF!" << endl;
      cout << " " << endl;
      cout << ">>>>>  Please, choose another event" << endl;
      cout << " " << endl;
      break;
    }

    if(!Good_TARGET_Event && !Good_TOF_Event && !Good_MWPC_Event){
      cout << "Event: "<< ievt << "   --  NO TARGET and NO TOF and NO MWPC!" << endl;
      cout << " " << endl;
      cout << ">>>>>  Please, choose another event" << endl;
      cout << " " << endl;
      break;
    }

    if(Event_On_Blacklist){
      cout << "Event: "<< ievt << " is on the blacklist." << endl;
      cout << " " << endl;
      cout << ">>>>>  Please, choose another event" << endl;
      cout << " " << endl;
      break;
    }
  }

  if(!Good_Event) return;

  for(int j=0 ; j<128 ; j++){
    if(ADC_High_SFT[j]<0)      ADC_High_SFT_corr[j]=0;
    if(ADC_High_SFT[j]>=0)     ADC_High_SFT_corr[j]=ADC_High_SFT[j];
  }




  //////////////////////// LEGENDS /////////////////////////////////////////////////

  //####### TARGET!









  //////// DETERMINE FIBER WITH HIGHEST AND SECOND HIGHEST LOW GAIN AMPLITUDE
  Int_t TDC_average = -1;

  int max_index_all[256];
  int max_ADC_all[256];
  int max_index_flag;

  for(int i = 0; i<256; i++){
    max_ADC_all[i] = -100000000;
    max_index_all[i] = -1;
  }

  for(int j = 0; j<256; j++){
    for(int i=0; i<256; i++){
      max_index_flag = 0;

      for(int k = 0; k<256; k++){
        if(i == max_index_all[k]){
        max_index_flag = 1;
        break;
        }
      }
      if (max_index_flag == 1) continue;

      else {
        if(ADC_Low_TARGET[i]>max_ADC_all[j]){
          max_index_all[j] = i;
          max_ADC_all[j] = ADC_Low_TARGET[i];
        }
      }
    }
  }

  //Calculate TDC Average
  int TDC_LG_max = -1;
  int TDC_LG_max2 = -1;
  int index_max1=0;
  int index_max2=0;

  for(int i = 0; i<256; i++){
    for(int j=0; j<6; j++){
      if(TDC_LG_max == -1){
        if(tdc_le_target[max_index_all[i]][j] >= TDC_Thr_min && tdc_le_target[max_index_all[i]][j]  <= TDC_Thr_max)
          TDC_LG_max = tdc_le_target[max_index_all[i]][j];
        index_max1 = max_index_all[i];
      }
      else if(TDC_LG_max2 == -1){
        if(tdc_le_target[max_index_all[i]][j] >= TDC_Thr_min && tdc_le_target[max_index_all[i]][j] <= TDC_Thr_max)
          TDC_LG_max2 = tdc_le_target[max_index_all[i]][j];
        index_max2 = max_index_all[i];
      }
    }
  }

  if(abs(TDC_LG_max - TDC_LG_max2) <= T_limit) TDC_average = (TDC_LG_max + TDC_LG_max2)/2;
  else if(TDC_LG_max > TDC_LG_max2) TDC_average = TDC_LG_max;
  else if(TDC_LG_max < TDC_LG_max2) TDC_average = TDC_LG_max2;

  int kaon_TDC_min = TDC_average + TDC_Avg_Offset_min;
  int kaon_TDC_max = TDC_average + TDC_Avg_Offset_max;
  int TDC_min_Kstop = TDC_average + TDC_Kstop_Avg_Offset_min;
  int TDC_max_Kstop = TDC_average + TDC_Kstop_Avg_Offset_max;

  Int_t TDC_min_TARGET = kaon_TDC_min;

  char ch_ADC_cut_TARGET[100];
  char ch_ADC_and_TDC_cut_Kstop[100];
  char ch_ADC_and_TDC_cut[100];

  bool has_TDC_hit[256] = {false};
  bool has_TDC_hit_Kstop[256] = {false};


  for(Int_t i=0; i<256; i++){
    for (Int_t k=0; k<4; k++) {
      if ((tdc_le_target[i][k]>=TDC_min_Kstop) && (tdc_le_target[i][k]<=TDC_max_Kstop)) has_TDC_hit_Kstop[i] = true;
      if ((tdc_le_target[i][k]>=kaon_TDC_min) && (tdc_le_target[i][k]<=kaon_TDC_max)) has_TDC_hit[i] = true;
    }
  }

  sprintf(ch_ADC_cut_TARGET,"(ADC offset = %d)",TARGET_ADC_Thr_HG_Offset);
  sprintf(ch_ADC_and_TDC_cut_Kstop,"(ADC: HG #geq %d, LG #geq %d | %d #leq TDC K Stop #leq %d)",HG_KAON,LG_KAON,TDC_min_Kstop,TDC_max_Kstop);
  sprintf(ch_ADC_and_TDC_cut,"(ADC offset = %d | %d #leq TDC #leq %d)",TARGET_ADC_Thr_HG_Offset,kaon_TDC_min,kaon_TDC_max);

  char str_lepton_fit1[100];
  char str_lepton_fit2[100];
  char str_lepton_fit3[100];
  char str_kaon_fit[100];
  char str_kaon_fit_final[100];
  char str_lepton_fit_final[100];
  char str_final[100];

  //sprintf(str_lepton_fit1,"Lepton Fit 1  |  Run %d ; Event %d",Run_Number,ievt);
  sprintf(str_lepton_fit1,"Lepton Fit 1  |  Run %d  --  Event %d",Run_Number,ievt);
  sprintf(str_lepton_fit2,"Lepton Fit 2  |  Run %d  --  Event %d",Run_Number,ievt);
  sprintf(str_lepton_fit3,"Lepton Fit 3  |  Run %d  --  Event %d",Run_Number,ievt);
  sprintf(str_kaon_fit,"Kaon Fit 1  |  Run %d  --  Event %d",Run_Number,ievt);
  sprintf(str_kaon_fit_final,"Kaon Fit  |  Run %d  --  Event %d",Run_Number,ievt);
  sprintf(str_lepton_fit_final,"Lepton Fit  |  Run %d  --  Event %d",Run_Number,ievt);
  sprintf(str_final,"Final Fit  |  Run %d  --  Event %d",Run_Number,ievt);

  // TEllipse *ell = new TEllipse(0, 0, R_TOF1, 0);
  // TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0);
  // TEllipse *ell_L1 = new TEllipse(0, 0, R_SFT_L1, 0);

  for(Int_t ii=0; ii<128; ii++){
    if(ADC_High_SFT_corr[ii]>ADC_High_corr_max) ADC_High_corr_max=ADC_High_SFT_corr[ii];
  }

  for(Int_t ii=0; ii<128; ii++){
    for (Int_t qq=0; qq<6; qq++){
      if(tdc_le_sft[ii][qq] > SFT_TDC_min[ii] && tdc_le_sft[ii][qq] < SFT_TDC_max[ii]) has_TDC_SFT_hit[ii]++;
    }
  }


  Int_t Angle_ADC_cut = 0;
  double x_inc = 0;
  double y_inc = 0;
  int hit_count = 0;
  int count = 0;

  int max_index = 0;
  int max_index2 = 0;
  int max_index3 = 0;
  int max_index4 = 0;

  int max_ADC = -100000000;
  int max_ADC2 = -100000000;
  int max_ADC3 = -100000000;
  int max_ADC4 = -100000000;

  for(Int_t q=0; q<256; q++){
    if(ADC_High_TARGET[q]>max_ADC){
      max_index = q;
      max_ADC = ADC_High_TARGET[q];
    }
  }

  for(Int_t q=0; q<256; q++){
    if (q == max_index) continue;
    else{
      if(ADC_High_TARGET[q]>max_ADC2){
        max_index2 = q;
        max_ADC2 = ADC_High_TARGET[q];
      }
    }
  }

  for(Int_t q=0; q<256; q++){
    if ((q == max_index) || (q == max_index2)) continue;
    else{
      if(ADC_High_TARGET[q]>max_ADC3) {
        max_index3 = q;
        max_ADC3 = ADC_High_TARGET[q];
      }
    }
  }

  for(Int_t q=0; q<256; q++){
    if ((q == max_index) || (q == max_index2) || (q == max_index3)) continue;
    else{
      if(ADC_High_TARGET[q]>max_ADC4) {
        max_index4 = q;
        max_ADC4 = ADC_High_TARGET[q];
      }
    }
  }

  double x_cent = Xloc[max_index];
  double y_cent = Yloc[max_index];
  double hyp[256] = {-1};

  for(Int_t j=0; j<256; j++){
    hyp[j] = sqrt(pow(x_cent - Xloc[j],2) + pow(y_cent - Yloc[j],2));
  }

  bool has_ADC_TOF1_hit[12] = {false};
  bool has_TDC_TOF1_hit[12] = {false};
  bool has_ADC_TOF2_hit[12] = {false};
  bool has_TDC_TOF2_hit[12] = {false};
  bool has_both_ADC_TOF1_hit[12] = {false};
  bool has_both_TDC_TOF1_hit[12] = {false};

  ///Set TOF2 Lines
  for(int i = 0; i < 12; i++){
    if ((ADC_TOF2AO[i]>0 && ADC_TOF2AI[i]>0) || (ADC_TOF2BO[i]>0 && ADC_TOF2BI[i]>0)) {has_ADC_TOF2_hit[i]=true;}
    if ((((TDC_TOF2AO[i]>TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<TOF2AO_TDC_max[i]) && (TDC_TOF2AI[i]>TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<TOF2AI_TDC_max[i])))
    ||  (((TDC_TOF2BO[i]>TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<TOF2BO_TDC_max[i]) && (TDC_TOF2BI[i]>TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<TOF2BI_TDC_max[i])))) {has_TDC_TOF2_hit[i]=true;}
  }



  ///Set TOF1 Lines
  for(int i = 0; i < 12; i++){
    if (ADC_TOF1U[i]>0 && ADC_TOF1D[i]>0) {has_both_ADC_TOF1_hit[i] = true;}
    if (ADC_TOF1U[i]>0 || ADC_TOF1D[i]>0) {has_ADC_TOF1_hit[i] = true;}
    if ((TDC_TOF1U[i]>TOF1U_TDC_min[i] && TDC_TOF1U[i]<TOF1U_TDC_max[i]) && (TDC_TOF1D[i]>TOF1D_TDC_min[i] && TDC_TOF1D[i]<TOF1D_TDC_max[i])) {has_both_TDC_TOF1_hit[i] = true;}
    if ((TDC_TOF1U[i]>TOF1U_TDC_min[i] && TDC_TOF1U[i]<TOF1U_TDC_max[i]) || (TDC_TOF1D[i]>TOF1D_TDC_min[i] && TDC_TOF1D[i]<TOF1D_TDC_max[i])) {has_TDC_TOF1_hit[i] = true;}
  }

    //int gap_hit[12] = {0};
  int ADC_TOF1_hit[12] = {0};
  int ADCTDC_TOF1_hit[12] = {0};
  int ADC_TOF2_hit[12] = {0};
  int ADCTDC_TOF2_hit[12] = {0};

  for (int k=0; k<12; k++) {
    if(has_ADC_TOF1_hit[k]){
      if(has_TDC_TOF1_hit[k]) ADCTDC_TOF1_hit[k]++;
      else ADC_TOF1_hit[k]++;
    }
    if(has_ADC_TOF2_hit[k]){
      if(has_TDC_TOF2_hit[k]) ADCTDC_TOF2_hit[k]++;
      else ADC_TOF2_hit[k]++;
    }
  }

  int selected_TOF2 = 0;

  // Determine which TOF2 is hit
  for(int i = 0; i<12; i++){
    if(has_TDC_TOF2_hit[i] && has_ADC_TOF2_hit[i]) selected_TOF2 = i + 1;
  }

  if(selected_TOF2 == 0){
    for(int i = 0; i<12; i++){
      if(has_TDC_TOF2_hit[i] || has_ADC_TOF2_hit[i]) selected_TOF2 = i+1;
    }
  }

  int gap_counter[12] = {0};

  //// GAP SCORING !
  int scoring_type = 2;         // scoring_type = 1  --->  Oscar's Method
                                // scoring_type = 2  --->  Sebastien's Method  TOF1[i]  TOF2[i-1], TOF2[i], TOF2[i+1]
                                // scoring_type = 3  --->  Sebastien's Method  TOF2[i]  TOF1[i-1], TOF1[i], TOF1[i+1]
  if(scoring_type==2){
    for(int i=0; i<12; i++){
      if(ADC_TOF1U[i]>=0) gap_counter[i]++;
    }

    for(int i=0; i<12; i++){
      if(ADC_TOF1D[i]>=0) gap_counter[i]++;
    }

    for(int i=0; i<12; i++){
      if(TDC_TOF1U[i]>=TOF1U_TDC_min[i] && TDC_TOF1U[i]<=TOF1U_TDC_max[i]) gap_counter[i]++;
    }

    for(int i=0; i<12; i++){
      if(TDC_TOF1D[i]>=TOF1D_TDC_min[i] && TDC_TOF1D[i]<=TOF1D_TDC_max[i]) gap_counter[i]++;
    }

    for(int i=0; i<12; i++){
      if(i!=0 && i!=11 && gap_counter[i]>0){

        if(ADC_TOF2AO[i]>=0) gap_counter[i]++;

        if(ADC_TOF2AO[i-1]>=0) gap_counter[i]++;

        if(ADC_TOF2AO[i+1]>=0) gap_counter[i]++;

        if(ADC_TOF2AI[i]>=0) gap_counter[i]++;

        if(ADC_TOF2AI[i-1]>=0) gap_counter[i]++;

        if(ADC_TOF2AI[i+1]>=0) gap_counter[i]++;

        if(ADC_TOF2BO[i]>=0) gap_counter[i]++;

        if(ADC_TOF2BO[i-1]>=0) gap_counter[i]++;

        if(ADC_TOF2BO[i+1]>=0) gap_counter[i]++;

        if(ADC_TOF2BI[i]>=0) gap_counter[i]++;

        if(ADC_TOF2BI[i-1]>=0) gap_counter[i]++;

        if(ADC_TOF2BI[i+1]>=0) gap_counter[i]++;

        if(TDC_TOF2AO[i]>=TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<=TOF2AO_TDC_max[i]) gap_counter[i]++;

        if(TDC_TOF2AO[i-1]>=TOF2AO_TDC_min[i-1] && TDC_TOF2AO[i-1]<=TOF2AO_TDC_max[i-1]) gap_counter[i]++;

        if(TDC_TOF2AO[i+1]>=TOF2AO_TDC_min[i+1] && TDC_TOF2AO[i+1]<=TOF2AO_TDC_max[i+1]) gap_counter[i]++;

        if(TDC_TOF2AI[i]>=TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<=TOF2AI_TDC_max[i]) gap_counter[i]++;

        if(TDC_TOF2AI[i-1]>=TOF2AI_TDC_min[i-1] && TDC_TOF2AI[i-1]<=TOF2AI_TDC_max[i-1]) gap_counter[i]++;

        if(TDC_TOF2AI[i+1]>=TOF2AI_TDC_min[i+1] && TDC_TOF2AI[i+1]<=TOF2AI_TDC_max[i+1]) gap_counter[i]++;

        if(TDC_TOF2BO[i]>=TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<=TOF2BO_TDC_max[i]) gap_counter[i]++;

        if(TDC_TOF2BO[i-1]>=TOF2BO_TDC_min[i-1] && TDC_TOF2BO[i-1]<=TOF2BO_TDC_max[i-1]) gap_counter[i]++;

        if(TDC_TOF2BO[i+1]>=TOF2BO_TDC_min[i+1] && TDC_TOF2BO[i+1]<=TOF2BO_TDC_max[i+1]) gap_counter[i]++;

        if(TDC_TOF2BI[i]>=TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<=TOF2BI_TDC_max[i]) gap_counter[i]++;

        if(TDC_TOF2BI[i-1]>=TOF2BI_TDC_min[i-1] && TDC_TOF2BI[i-1]<=TOF2BI_TDC_max[i-1]) gap_counter[i]++;

        if(TDC_TOF2BI[i+1]>=TOF2BI_TDC_min[i+1] && TDC_TOF2BI[i+1]<=TOF2BI_TDC_max[i+1]) gap_counter[i]++;
      }
    }

    if(gap_counter[0]>0){
      if(ADC_TOF2AO[0]>=0) gap_counter[0]++;

      if(ADC_TOF2AO[1]>=0) gap_counter[0]++;

      if(ADC_TOF2AO[11]>=0) gap_counter[0]++;

      if(ADC_TOF2AI[0]>=0) gap_counter[0]++;

      if(ADC_TOF2AI[1]>=0) gap_counter[0]++;

      if(ADC_TOF2AI[11]>=0) gap_counter[0]++;

      if(ADC_TOF2BO[0]>=0) gap_counter[0]++;

      if(ADC_TOF2BO[1]>=0) gap_counter[0]++;

      if(ADC_TOF2BO[11]>=0) gap_counter[0]++;

      if(ADC_TOF2BI[0]>=0) gap_counter[0]++;

      if(ADC_TOF2BI[1]>=0) gap_counter[0]++;

      if(ADC_TOF2BI[11]>=0) gap_counter[0]++;

      if(TDC_TOF2AO[0]>=TOF2AO_TDC_min[0] && TDC_TOF2AO[0]<=TOF2AO_TDC_max[0]) gap_counter[0]++;

      if(TDC_TOF2AO[1]>=TOF2AO_TDC_min[1] && TDC_TOF2AO[1]<=TOF2AO_TDC_max[1]) gap_counter[0]++;

      if(TDC_TOF2AO[11]>=TOF2AO_TDC_min[11] && TDC_TOF2AO[11]<=TOF2AO_TDC_max[11]) gap_counter[0]++;

      if(TDC_TOF2AI[0]>=TOF2AI_TDC_min[0] && TDC_TOF2AI[0]<=TOF2AI_TDC_max[0]) gap_counter[0]++;

      if(TDC_TOF2AI[1]>=TOF2AI_TDC_min[1] && TDC_TOF2AI[1]<=TOF2AI_TDC_max[1]) gap_counter[0]++;

      if(TDC_TOF2AI[11]>=TOF2AI_TDC_min[11] && TDC_TOF2AI[11]<=TOF2AI_TDC_max[11]) gap_counter[0]++;

      if(TDC_TOF2BO[0]>=TOF2BO_TDC_min[0] && TDC_TOF2BO[0]<=TOF2BO_TDC_max[0]) gap_counter[0]++;

      if(TDC_TOF2BO[1]>=TOF2BO_TDC_min[1] && TDC_TOF2BO[1]<=TOF2BO_TDC_max[1]) gap_counter[0]++;

      if(TDC_TOF2BO[11]>=TOF2BO_TDC_min[11] && TDC_TOF2BO[11]<=TOF2BO_TDC_max[11]) gap_counter[0]++;

      if(TDC_TOF2BI[0]>=TOF2BI_TDC_min[0] && TDC_TOF2BI[0]<=TOF2BI_TDC_max[0]) gap_counter[0]++;

      if(TDC_TOF2BI[1]>=TOF2BI_TDC_min[1] && TDC_TOF2BI[1]<=TOF2BI_TDC_max[1]) gap_counter[0]++;

      if(TDC_TOF2BI[11]>=TOF2BI_TDC_min[11] && TDC_TOF2BI[11]<=TOF2BI_TDC_max[11]) gap_counter[0]++;
    }

    if(gap_counter[11]>0){
      if(ADC_TOF2AO[11]>=0) gap_counter[11]++;
      if(ADC_TOF2AO[10]>=0) gap_counter[11]++;
      if(ADC_TOF2AO[0]>=0) gap_counter[11]++;

      if(ADC_TOF2AI[11]>=0) gap_counter[11]++;
      if(ADC_TOF2AI[10]>=0) gap_counter[11]++;
      if(ADC_TOF2AI[0]>=0) gap_counter[11]++;

      if(ADC_TOF2BO[11]>=0) gap_counter[11]++;
      if(ADC_TOF2BO[10]>=0) gap_counter[11]++;
      if(ADC_TOF2BO[0]>=0) gap_counter[11]++;

      if(ADC_TOF2BI[11]>=0) gap_counter[11]++;
      if(ADC_TOF2BI[10]>=0) gap_counter[11]++;
      if(ADC_TOF2BI[0]>=0) gap_counter[11]++;

      if(TDC_TOF2AO[11]>=TOF2AO_TDC_min[11] && TDC_TOF2AO[11]<=TOF2AO_TDC_max[11]) gap_counter[11]++;
      if(TDC_TOF2AO[10]>=TOF2AO_TDC_min[10] && TDC_TOF2AO[10]<=TOF2AO_TDC_max[10]) gap_counter[11]++;
      if(TDC_TOF2AO[0]>=TOF2AO_TDC_min[0] && TDC_TOF2AO[0]<=TOF2AO_TDC_max[0]) gap_counter[11]++;

      if(TDC_TOF2AI[11]>=TOF2AI_TDC_min[11] && TDC_TOF2AI[11]<=TOF2AI_TDC_max[11]) gap_counter[11]++;
      if(TDC_TOF2AI[10]>=TOF2AI_TDC_min[10] && TDC_TOF2AI[10]<=TOF2AI_TDC_max[10]) gap_counter[11]++;
      if(TDC_TOF2AI[0]>=TOF2AI_TDC_min[0] && TDC_TOF2AI[0]<=TOF2AI_TDC_max[0]) gap_counter[11]++;

      if(TDC_TOF2BO[11]>=TOF2BO_TDC_min[11] && TDC_TOF2BO[11]<=TOF2BO_TDC_max[11]) gap_counter[11]++;
      if(TDC_TOF2BO[10]>=TOF2BO_TDC_min[10] && TDC_TOF2BO[10]<=TOF2BO_TDC_max[10]) gap_counter[11]++;
      if(TDC_TOF2BO[0]>=TOF2BO_TDC_min[0] && TDC_TOF2BO[0]<=TOF2BO_TDC_max[0]) gap_counter[11]++;

      if(TDC_TOF2BI[11]>=TOF2BI_TDC_min[11] && TDC_TOF2BI[11]<=TOF2BI_TDC_max[11]) gap_counter[11]++;
      if(TDC_TOF2BI[10]>=TOF2BI_TDC_min[10] && TDC_TOF2BI[10]<=TOF2BI_TDC_max[10]) gap_counter[11]++;
      if(TDC_TOF2BI[0]>=TOF2BI_TDC_min[0] && TDC_TOF2BI[0]<=TOF2BI_TDC_max[0]) gap_counter[11]++;
    }
  }

  int high_gap_hit = 0;
  int gap_to_fit = 0;
  int score_max = 0;

  vector<int> tof1_ties;

  for (int k=0; k<12; k++) {
    if(gap_counter[k]>=high_gap_hit){
      if(gap_counter[k] == high_gap_hit) tof1_ties.push_back(k);
      else{
        tof1_ties.clear();
        tof1_ties.push_back(k);
      }
      high_gap_hit = gap_counter[k];
    }
  }

  if(tof1_ties.size() > 1){  // only when there are more than 2 elements in tof1_ties!
  	for(unsigned int j=0; j<tof1_ties.size(); j++){
  		for(int jj=0; jj<8; jj++){
  			if(ADC_High_TARGET[channel[tof1_ties[j]][jj]] > 0) gap_counter[tof1_ties[j]]++;
  			if(ADC_Low_TARGET[channel[tof1_ties[j]][jj]] > 0) gap_counter[tof1_ties[j]]++;
  		}
  	}
  }

  for(int k=0; k<12; k++){
  	if(gap_counter[k] > score_max){
  		score_max = gap_counter[k];
  		gap_to_fit = k+1;
  	}
  }

  // Determine the k_stop bars
  bool k_stop_bar[256] = {false};
  bool k_stop_bar_initial[256] = {false};
  vector<int> good_k_stop_bars;

  for(int i = 0; i<256; i++){
    if(ADC_High_TARGET[i]>0 || ADC_Low_TARGET[i]>0) vec_TARGET_bar_selected.push_back(i);
    if(ADC_High_TARGET[i] > HG_KAON && ADC_Low_TARGET[i] > LG_KAON && has_TDC_hit_Kstop[i]){
      good_k_stop_bars.push_back(i);
      k_stop_bar[i] = true;
      k_stop_bar_initial[i] = true;
    }
  }

  // Fill KStop bars
  for(int i = 0; i<256; i++){
    if(k_stop_bar[i]){
      vec_xx_kaon.push_back(Xloc[i]);
      vec_ex_kaon.push_back(TARGET_Errors_X);

      vec_yy_kaon.push_back(Yloc[i]);
      vec_ey_kaon.push_back(TARGET_Errors_Y);

      vec_kaon_bars.push_back(i);
    }
  }

  //delete gr_kaon;               gr_kaon = nullptr;

  // TGraph gr_kaon = TGraphErrors(vec_xx_kaon.size(),&vec_xx_kaon[0],&vec_yy_kaon[0],
  //   &vec_ey_kaon[0],&vec_ey_kaon[0]);
  TGraph gr_kaon = TGraph(vec_xx_kaon.size(),&vec_xx_kaon[0],&vec_yy_kaon[0]);
  gr_kaon.SetMarkerStyle(21);
  gr_kaon.SetMarkerColor(kBlue-6);
  gr_kaon.SetMarkerSize(0.8);
  gr_kaon.GetXaxis()->SetLimits(-50.,50.);
  gr_kaon.GetYaxis()->SetRangeUser(-50.,50.);

  // TGraph *gr_kaon_bk;
  TGraph gr_kaon_bk = TGraph(vec_xx_kaon.size(),&vec_xx_kaon[0],&vec_yy_kaon[0]);

  // TF1 *gr_kaon_fit = new TF1("kaon_fit", "pol1");
  TF1 gr_kaon_fit = TF1("kaon_fit", "pol1");
  TF1 *gr_kaon_fit_ptr;
  double a_fit_kaon = 0.;
  double b_fit_kaon = 0.;
  double Chis_kaon = 0.;
  int ndf_kaon = 99;
  bool kaon_bk = false;




  if(vec_xx_kaon.size()>=5){
    gr_kaon.Fit("kaon_fit","Q");
    gr_kaon_fit_ptr = gr_kaon.GetFunction("kaon_fit");
    gr_kaon_fit_ptr->SetLineWidth(2);
    gr_kaon_fit_ptr->SetLineColor(4);
    a_fit_kaon = gr_kaon_fit_ptr->GetParameter(1);
    b_fit_kaon = gr_kaon_fit_ptr->GetParameter(0);
    Chis_kaon = gr_kaon_fit_ptr->GetChisquare();
    ndf_kaon = gr_kaon_fit_ptr->GetNDF();
  }


  if(a_fit_kaon > 1000){
    kaon_bk = true;

    gr_kaon_bk.SetMarkerStyle(21);
    gr_kaon_bk.SetMarkerColor(kBlue-6);
    gr_kaon_bk.SetMarkerSize(0.8);
    gr_kaon_bk.GetXaxis()->SetLimits(-50.,50.);
    gr_kaon_bk.GetYaxis()->SetRangeUser(-50.,50.);

    gr_kaon_bk.Fit("kaon_fit","Q");
    gr_kaon_fit_ptr = gr_kaon_bk.GetFunction("kaon_fit");
    gr_kaon_fit_ptr->SetLineWidth(2);
    gr_kaon_fit_ptr->SetLineColor(4);
    a_fit_kaon = gr_kaon_fit_ptr->GetParameter(1);
    b_fit_kaon = gr_kaon_fit_ptr->GetParameter(0);
  }


  // Determine if a hit target has any hit neighbours
  for(int i = 0; i<256; i++){
    for(int j=0; j<8; j++){
      if((TARGET_neighbours[i][j] != -1 && ADC_High_TARGET[TARGET_neighbours[i][j]] > Angle_ADC_cut && has_TDC_hit[TARGET_neighbours[i][j]] && !k_stop_bar[TARGET_neighbours[i][j]]) ||
         (TARGET_neighbours[i][j] != -1 && ADC_High_TARGET[TARGET_neighbours[i][j]] <0 && ADC_Low_TARGET[TARGET_neighbours[i][j]]>0 && Switch==1 && !k_stop_bar[TARGET_neighbours[i][j]])){
          TARGET_High_has_neighbours[i] = true;
          break;
      }
    }
  }


  for(Int_t i=0; i<256; i++){
    if(TARGET_High_has_neighbours[i]){
      if(ADC_High_TARGET[i]>Angle_ADC_cut && has_TDC_hit[i]){
        if(!k_stop_bar[i]){

          if(gap_to_fit==6 || gap_to_fit==12){


            vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[i]]);
            vec_ex_lepton.push_back(TARGET_Errors_X);

            vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[i]]);
            vec_ey_lepton.push_back(TARGET_Errors_Y);

            lepton_bars.push_back(i);
          }
          else{
            vec_xx_lepton.push_back(Xloc[i]);
            vec_ex_lepton.push_back(TARGET_Errors_X);

            vec_yy_lepton.push_back(Yloc[i]);
            vec_ey_lepton.push_back(TARGET_Errors_Y);

            lepton_bars.push_back(i);
          }

         vec_lepton_size.push_back(Xloc[i]);
        }
        count++;

        if(IsIn(i,channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],   // Additional weight for fibers close to the edge if hit
                  channel[gap_to_fit-1][2], channel[gap_to_fit-1][3],
                  channel[gap_to_fit-1][4], channel[gap_to_fit-1][5],
                  channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){
          if(!k_stop_bar[i]){

            if(gap_to_fit==6 || gap_to_fit==12){
              vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[i]]);
              vec_ex_lepton.push_back(TARGET_Errors_X);
              vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[i]]);
              vec_ex_lepton.push_back(TARGET_Errors_X);

              vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[i]]);
              vec_ey_lepton.push_back(TARGET_Errors_Y);
              vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[i]]);
              vec_ey_lepton.push_back(TARGET_Errors_Y);
            }
            else{
              vec_xx_lepton.push_back(Xloc[i]);
              vec_ex_lepton.push_back(TARGET_Errors_X);
              vec_xx_lepton.push_back(Xloc[i]);
              vec_ex_lepton.push_back(TARGET_Errors_X);

              vec_yy_lepton.push_back(Yloc[i]);
              vec_ey_lepton.push_back(TARGET_Errors_Y);
              vec_yy_lepton.push_back(Yloc[i]);
              vec_ey_lepton.push_back(TARGET_Errors_Y);
            }

            vec_lepton_size.push_back(Xloc[i]);
          }
        }
      }

      if(ADC_High_TARGET[i]>Angle_ADC_cut){
        x_inc = x_inc + Xloc[i];
        y_inc = y_inc + Yloc[i];
        hit_count++;
      }

      if(ADC_High_TARGET[i]<0 && ADC_Low_TARGET[i]>0 && Switch==1){
        if(!k_stop_bar[i]){

          if(gap_to_fit==6 || gap_to_fit==12){
            vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[i]]);
            vec_ex_lepton.push_back(TARGET_Errors_X);

            vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[i]]);
            vec_ey_lepton.push_back(TARGET_Errors_Y);

            lepton_bars.push_back(i);
          }
          else{
            vec_xx_lepton.push_back(Xloc[i]);
            vec_ex_lepton.push_back(TARGET_Errors_X);

            vec_yy_lepton.push_back(Yloc[i]);
            vec_ey_lepton.push_back(TARGET_Errors_Y);

            lepton_bars.push_back(i);
          }

          vec_lepton_size.push_back(Xloc[i]);
        }
        count++;
      }
    }
  }

  for(int g=0; g<12; g++){
    Gap[g][0][0] = TOF_Xloc[3*g];
    Gap[g][1][0] = TOF_Xloc[3*g+1];
    Gap[g][2][0] = TOF_Xloc[3*g+2];

    Gap[g][0][1] = TOF_Yloc[3*g];
    Gap[g][1][1] = TOF_Yloc[3*g+1];
    Gap[g][2][1] = TOF_Yloc[3*g+2];
  }

  for(int i = 0; i<3; i++){

    if(gap_to_fit == 6){
      vec_xx_lepton.push_back(Gap[2][1][0]);
      vec_ex_lepton.push_back(TOF1_Errors_X[2][1]);

      vec_yy_lepton.push_back(Gap[2][1][1]);
      vec_ey_lepton.push_back(TOF1_Errors_Y[2][1]);
    }
    else if(gap_to_fit == 12){
      vec_xx_lepton.push_back(Gap[8][1][0]);
      vec_ex_lepton.push_back(TOF1_Errors_X[2][1]);

      vec_yy_lepton.push_back(Gap[8][1][1]);
      vec_ey_lepton.push_back(TOF1_Errors_Y[2][1]);
    }
    else{
      vec_xx_lepton.push_back(Gap[gap_to_fit-1][1][0]);
      vec_ex_lepton.push_back(TOF1_Errors_X[gap_to_fit-1][1]);

      vec_yy_lepton.push_back(Gap[gap_to_fit-1][1][1]);
      vec_ey_lepton.push_back(TOF1_Errors_Y[gap_to_fit-1][1]);
    }

  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // TGraph *gr_lepton_1 = new TGraphErrors(vec_xx_lepton.size(),&vec_xx_lepton[0],&vec_yy_lepton[0],
  //   &vec_ex_lepton[0],&vec_ey_lepton[0]);
  //
  // gr_lepton_1->SetMarkerStyle(21);
  // gr_lepton_1->SetMarkerColor(2);
  // gr_lepton_1->SetMarkerSize(0.8);
  // gr_lepton_1->GetXaxis()->SetLimits(-50.,50.);
  // gr_lepton_1->GetYaxis()->SetRangeUser(-50.,50.);
  // TGraph gr_lepton_1 = TGraphErrors(vec_xx_lepton.size(),&vec_xx_lepton[0],&vec_yy_lepton[0],
  //   &vec_ex_lepton[0],&vec_ey_lepton[0]);
  TGraph gr_lepton_1 = TGraph(vec_xx_lepton.size(),&vec_xx_lepton[0],&vec_yy_lepton[0]);

  gr_lepton_1.SetMarkerStyle(21);
  gr_lepton_1.SetMarkerColor(2);
  gr_lepton_1.SetMarkerSize(0.8);
  gr_lepton_1.GetXaxis()->SetLimits(-50.,50.);
  gr_lepton_1.GetYaxis()->SetRangeUser(-50.,50.);

  double a_lepton_fit_1 = 0.;
  double b_lepton_fit_1 = 0.;
  double Chis_lepton_fit_1 = 0.;
  int ndf_lepton_fit_1 = 0;

  // TF1 *func_lepton_fit_1 = new TF1("Lepton_fit_1", "pol1");
  TF1 func_lepton_fit_1 = TF1("Lepton_fit_1", "pol1");

  // if(gap_to_fit==1 || gap_to_fit==2 || gap_to_fit==7 || gap_to_fit==8){
  //   func_lepton_fit_1->SetParameter(0,0);
  //   func_lepton_fit_1->SetParameter(1,1);
  // }
  // if(gap_to_fit==4 || gap_to_fit==5 || gap_to_fit==10 || gap_to_fit==11){
  //   func_lepton_fit_1->SetParameter(0,0);
  //   func_lepton_fit_1->SetParameter(1,-1);
  // }
  // else{
  //   func_lepton_fit_1->SetParameter(0,0);
  //   func_lepton_fit_1->SetParameter(1,1);
  // }
  //
  // func_lepton_fit_1->SetParLimits(0,-50,50);
  // func_lepton_fit_1->SetParLimits(1,-50,50);
  if(gap_to_fit==1 || gap_to_fit==2 || gap_to_fit==7 || gap_to_fit==8){
    func_lepton_fit_1.SetParameter(0,0);
    func_lepton_fit_1.SetParameter(1,1);
  }
  if(gap_to_fit==4 || gap_to_fit==5 || gap_to_fit==10 || gap_to_fit==11){
    func_lepton_fit_1.SetParameter(0,0);
    func_lepton_fit_1.SetParameter(1,-1);
  }
  else{
    func_lepton_fit_1.SetParameter(0,0);
    func_lepton_fit_1.SetParameter(1,1);
  }

  func_lepton_fit_1.SetParLimits(0,-50,50);
  func_lepton_fit_1.SetParLimits(1,-50,50);
  // return;
  gr_lepton_1.Fit("Lepton_fit_1","QB");

  // delete func_lepton_fit_1;

  // TF1 *func_lepton_fit_1_ptr = gr_lepton_1->GetFunction("Lepton_fit_1");
  // func_lepton_fit_1->SetLineWidth(2);
  // func_lepton_fit_1->SetLineColor(2);
  // a_lepton_fit_1 = func_lepton_fit_1->GetParameter(1);
  // b_lepton_fit_1 = func_lepton_fit_1->GetParameter(0);
  // Chis_lepton_fit_1 = func_lepton_fit_1->GetChisquare();
  // ndf_lepton_fit_1 = func_lepton_fit_1->GetNDF();
  TF1 *func_lepton_fit_1_ptr = gr_lepton_1.GetFunction("Lepton_fit_1");
  func_lepton_fit_1_ptr->SetLineWidth(2);
  func_lepton_fit_1_ptr->SetLineColor(2);
  a_lepton_fit_1 = func_lepton_fit_1_ptr->GetParameter(1);
  b_lepton_fit_1 = func_lepton_fit_1_ptr->GetParameter(0);
  Chis_lepton_fit_1 = func_lepton_fit_1_ptr->GetChisquare();
  ndf_lepton_fit_1 = func_lepton_fit_1_ptr->GetNDF();




////////////////////////////////////////////////////////////////////////////////////////////////////////

  for(int k=0; k<256; k++){
    for(int l=0; l<6; l++){
        TDC_LE_TARGET_corrected[k][l] = 9999;
        if(tdc_le_target[k][l]>-1)  TDC_LE_TARGET_corrected[k][l] = tdc_le_target[k][l] - TDC_average;
    }
  }

  for(int k=0; k<256; k++){
    for(int l=0; l<6; l++){
        if(TDC_LE_TARGET_corrected[k][l] == 9999) sprintf(TDC_LE_TARGET_corr[k][l],"----");
        if(TDC_LE_TARGET_corrected[k][l] != 9999) sprintf(TDC_LE_TARGET_corr[k][l],"%4d",TDC_LE_TARGET_corrected[k][l]);
    }
  }

  //// Don't do anything if the event has less than n_hits hits in the TARGET
  if (enable_cout == 9) count = n_hit+1;
  if(count<n_hit){
    gROOT->Reset();
    gROOT->Clear();
    cout << " >>>>  Event "<< ievt << " has fewer than "<< n_hit << " high gain hits in the TARGET (outliers and K stop removed)!" << endl;
    cout << " >>>>  Please, choose another event" << endl;
    cout << " " << endl;
    return;
  }

  // Arrays for finding mwpc clustering
  double C2X_L[56] = {0.};
  double C2X_R[56] = {0.};
  double C2Y_L[16] = {0.};
  double C2Y_R[16] = {0.};
  double C3X_L[64] = {0.};
  double C3X_R[64] = {0.};
  double C3Y_L[16] = {0.};
  double C3Y_R[16] = {0.};
  double C4X_L[72] = {0.};
  double C4X_R[72] = {0.};
  double C4Y_L[16] = {0.};
  double C4Y_R[16] = {0.};

  double C2X_centroid = 0.0;

  //C2 Counters
  for(int i=0; i<56; i++){
    if (Good_Event && ADC_C2X_L[i]>0) C2X_L[i] = ADC_C2X_L[i];
    if (Good_Event && ADC_C2X_R[i]>0) C2X_R[i] = ADC_C2X_R[i];
  }

  for(int i=0; i<64; i++){
    if (Good_Event && ADC_C3X_L[i]>0) C3X_L[i] = ADC_C3X_L[i];
    if (Good_Event && ADC_C3X_R[i]>0) C3X_R[i] = ADC_C3X_R[i];
  }

  for(int i=0; i<72; i++){
    if (Good_Event && ADC_C4X_L[i]>0) C4X_L[i] = ADC_C4X_L[i];
    if (Good_Event && ADC_C4X_R[i]>0) C4X_R[i] = ADC_C4X_R[i];
  }

  for(int i=0; i<16; i++){
    if (Good_Event && ADC_C2Y_L[i]>0) C2Y_L[i] = ADC_C2Y_L[i];
    if (Good_Event && ADC_C2Y_R[i]>0) C2Y_R[i] = ADC_C2Y_R[i];
    if (Good_Event && ADC_C3Y_L[i]>0) C3Y_L[i] = ADC_C3Y_L[i];
    if (Good_Event && ADC_C3Y_R[i]>0) C3Y_R[i] = ADC_C3Y_R[i];
    if (Good_Event && ADC_C4Y_L[i]>0) C4Y_L[i] = ADC_C4Y_L[i];
    if (Good_Event && ADC_C4Y_R[i]>0) C4Y_R[i] = ADC_C4Y_R[i];
  }

  // Find clustering of mwpcs
  ////////////////////////////
  bool first_cluster = true;
  int cluster_spacing = 0;
  int cluster_length_count = 0;

  int C2X_clusters = 0;
  int C2Y_clusters = 0;
  int C3X_clusters = 0;
  int C3Y_clusters = 0;
  int C4X_clusters = 0;
  int C4Y_clusters = 0;

  vector<int> C2X_cluster_index; // Hold the starting index of each cluster.
  vector<int> C2X_cluster_length; // Hold the length of each cluster.
  vector<int> C2Y_cluster_index;
  vector<int> C2Y_cluster_length;

  vector<int> C3X_cluster_index;
  vector<int> C3X_cluster_length;
  vector<int> C3Y_cluster_index;
  vector<int> C3Y_cluster_length;

  vector<int> C4X_cluster_index;
  vector<int> C4X_cluster_length;
  vector<int> C4Y_cluster_index;
  vector<int> C4Y_cluster_length;

  if(selected_TOF2 > 6){  // LEFT
    first_cluster = true;
    cluster_spacing = 0;
    cluster_length_count = 0;

    for(int i = 0; i<56; i++){
      if(C2X_L[i] > 0. && first_cluster){
        cluster_spacing = MWPC_cluster_separation + 1;
        first_cluster = false;
      }

      if(C2X_L[i] > 0.){
        if(cluster_spacing > MWPC_cluster_separation){
          C2X_clusters++;
          C2X_cluster_index.push_back(i);
        }
        cluster_length_count++;
        cluster_spacing = 0;

        if(i+2 < 56 && C2X_L[i+1] <= 0. && C2X_L[i+2] <= 0.) C2X_cluster_length.push_back(cluster_length_count);
        else if(i + 2 == 56 && C2X_L[i+1] <= 0.) C2X_cluster_length.push_back(cluster_length_count);
        else if(i + 1 == 56) C2X_cluster_length.push_back(cluster_length_count);
      }
      else{
        cluster_spacing++;
        if(i != 0 && C2X_L[i-1] <= 0.) cluster_length_count = 0;
        else if ( i != 0) cluster_length_count++;
      }
    }

    // count C2Y clusters
    first_cluster = true;
    cluster_spacing = 0;
    cluster_length_count = 0;

    for(int i = 0; i<16; i++){
      if(C2Y_L[i] > 0 && first_cluster){
        cluster_spacing = MWPC_cluster_separation + 1;
        first_cluster = false;
      }

      if(C2Y_L[i] > 0){
        if(cluster_spacing > MWPC_cluster_separation){
          C2Y_clusters++;
          C2Y_cluster_index.push_back(i);
        }
        cluster_length_count++;
        cluster_spacing = 0;

        if(i+2 < 16 && C2Y_L[i+1] <= 0 && C2Y_L[i+2] <= 0) C2Y_cluster_length.push_back(cluster_length_count);
        else if(i + 2 == 16 && C2Y_L[i+1] <= 0) C2Y_cluster_length.push_back(cluster_length_count);
        else if(i + 1 == 16) C2Y_cluster_length.push_back(cluster_length_count);
      }
      else{
        cluster_spacing++;
        if(i != 0 && C2Y_L[i-1] <= 0) cluster_length_count = 0;
        else if ( i != 0) cluster_length_count++;
      }
    }

    // count C3X clusters
    first_cluster = true;
    cluster_spacing = 0;
    cluster_length_count = 0;

    for(int i = 0; i<64; i++){
      if(C3X_L[i] > 0 && first_cluster){
        cluster_spacing = MWPC_cluster_separation + 1;
        first_cluster = false;
      }

      if(C3X_L[i] > 0){
        if(cluster_spacing > MWPC_cluster_separation){
          C3X_clusters++;
          C3X_cluster_index.push_back(i);
        }
        cluster_length_count++;
        cluster_spacing = 0;

        if(i+2 < 64 && C3X_L[i+1] <= 0 && C3X_L[i+2] <= 0) C3X_cluster_length.push_back(cluster_length_count);
        else if(i + 2 == 64 && C3X_L[i+1] <= 0) C3X_cluster_length.push_back(cluster_length_count);
        else if(i + 1 == 64) C3X_cluster_length.push_back(cluster_length_count);
      }
      else{
        cluster_spacing++;
        if(i != 0 && C3X_L[i-1] <= 0) cluster_length_count = 0;
        else if ( i != 0) cluster_length_count++;
      }
    }

    // count C3Y clusters
    first_cluster = true;
    cluster_spacing = 0;
    cluster_length_count = 0;

    for(int i = 0; i<16; i++){
      if(C3Y_L[i] > 0 && first_cluster){
        cluster_spacing = MWPC_cluster_separation + 1;
        first_cluster = false;
      }

      if(C3Y_L[i] > 0){
        if(cluster_spacing > MWPC_cluster_separation){
          C3Y_clusters++;
          C3Y_cluster_index.push_back(i);
        }
        cluster_length_count++;
        cluster_spacing = 0;

        if(i+2 < 16 && C3Y_L[i+1] <= 0 && C3Y_L[i+2] <= 0) C3Y_cluster_length.push_back(cluster_length_count);
        else if(i + 2 == 16 && C3Y_L[i+1] <= 0) C3Y_cluster_length.push_back(cluster_length_count);
        else if(i + 1 == 16) C3Y_cluster_length.push_back(cluster_length_count);
      }
      else{
        cluster_spacing++;
        if(i != 0 && C3Y_L[i-1] <= 0) cluster_length_count = 0;
        else if ( i != 0) cluster_length_count++;
      }
    }

    // count C4X clusters
    first_cluster = true;
    cluster_spacing = 0;
    cluster_length_count = 0;

    for(int i = 0; i<72; i++){
      if(C4X_L[i] > 0 && first_cluster){
        cluster_spacing = MWPC_cluster_separation + 1;
        first_cluster = false;
      }

      if(C4X_L[i] > 0){
        if(cluster_spacing > MWPC_cluster_separation){
          C4X_clusters++;
          C4X_cluster_index.push_back(i);
        }
        cluster_length_count++;
        cluster_spacing = 0;

        if(i+2 < 72 && C4X_L[i+1] <= 0 && C4X_L[i+2] <= 0) C4X_cluster_length.push_back(cluster_length_count);
        else if(i == 70 && C4X_L[71] <= 0) C4X_cluster_length.push_back(cluster_length_count);
        else if(i == 71) C4X_cluster_length.push_back(cluster_length_count);
      }
      else{
        cluster_spacing++;
        if(i != 0 && C4X_L[i-1] <= 0) cluster_length_count = 0;
        else if ( i != 0) cluster_length_count++;
      }
    }

    // count C4Y clusters
    first_cluster = true;
    cluster_spacing = 0;
    cluster_length_count = 0;

    for(int i = 0; i<16; i++){
      if(C4Y_L[i] > 0 && first_cluster){
        cluster_spacing = MWPC_cluster_separation + 1;
        first_cluster = false;
      }

      if(C4Y_L[i] > 0){
        if(cluster_spacing > MWPC_cluster_separation){
          C4Y_clusters++;
          C4Y_cluster_index.push_back(i);
        }
        cluster_length_count++;
        cluster_spacing = 0;

        if(i+2 < 16 && C4Y_L[i+1] <= 0 && C4Y_L[i+2] <= 0) C4Y_cluster_length.push_back(cluster_length_count);
        else if(i + 2 == 16 && C4Y_L[i+1] <= 0) C4Y_cluster_length.push_back(cluster_length_count);
        else if(i + 1 == 16) C4Y_cluster_length.push_back(cluster_length_count);
      }
      else{
        cluster_spacing++;
        if(i != 0 && C4Y_L[i-1] <= 0) cluster_length_count = 0;
        else if ( i != 0) cluster_length_count++;
      }
    }
  }
  else{  // RIGHT
    // count C2X clusters
    first_cluster = true;
    cluster_spacing = 0;
    cluster_length_count = 0;

    for(int i = 0; i<56; i++){
      if(C2X_R[i] > 0 && first_cluster){
        cluster_spacing = MWPC_cluster_separation + 1;
        first_cluster = false;
      }

      if(C2X_R[i] > 0){
        if(cluster_spacing > MWPC_cluster_separation){
          C2X_clusters++;
          C2X_cluster_index.push_back(i);
        }
        cluster_length_count++;
        cluster_spacing = 0;

        if(i+2 < 56 && C2X_R[i+1] <= 0 && C2X_R[i+2] <= 0) C2X_cluster_length.push_back(cluster_length_count);
        else if(i + 2 == 56 && C2X_R[i+1] <= 0) C2X_cluster_length.push_back(cluster_length_count);
        else if(i + 1 == 56)C2X_cluster_length.push_back(cluster_length_count);
      }
      else{
        cluster_spacing++;
        if(i != 0 && C2X_R[i-1] <= 0) cluster_length_count = 0;
        else if ( i != 0) cluster_length_count++;
      }
    }

    // count C2Y clusters
    first_cluster = true;
    cluster_spacing = 0;
    cluster_length_count = 0;

    for(int i = 0; i<16; i++){
      if(C2Y_R[i] > 0 && first_cluster){
        cluster_spacing = MWPC_cluster_separation + 1;
        first_cluster = false;
      }

      if(C2Y_R[i] > 0){
        if(cluster_spacing > MWPC_cluster_separation){
          C2Y_clusters++;
          C2Y_cluster_index.push_back(i);
        }
        cluster_length_count++;
        cluster_spacing = 0;

        if(i+2 < 16 && C2Y_R[i+1] <= 0 && C2Y_R[i+2] <= 0) C2Y_cluster_length.push_back(cluster_length_count);
        else if(i + 2 == 16 && C2Y_R[i+1] <= 0) C2Y_cluster_length.push_back(cluster_length_count);
        else if(i + 1 == 16) C2Y_cluster_length.push_back(cluster_length_count);
      }
      else{
        cluster_spacing++;
        if(i != 0 && C2Y_R[i-1] <= 0) cluster_length_count = 0;
        else if ( i != 0) cluster_length_count++;
      }
    }

    // count C3X clusters
    first_cluster = true;
    cluster_spacing = 0;
    cluster_length_count = 0;

    for(int i = 0; i<64; i++){
      if(C3X_R[i] > 0 && first_cluster){
        cluster_spacing = MWPC_cluster_separation + 1;
        first_cluster = false;
      }

      if(C3X_R[i] > 0){
        if(cluster_spacing > MWPC_cluster_separation){
          C3X_clusters++;
          C3X_cluster_index.push_back(i);
        }
        cluster_length_count++;
        cluster_spacing = 0;

        if(i+2 < 64 && C3X_R[i+1] <= 0 && C3X_R[i+2] <= 0) C3X_cluster_length.push_back(cluster_length_count);
        else if(i + 2 == 64 && C3X_R[i+1] <= 0) C3X_cluster_length.push_back(cluster_length_count);
        else if(i + 1 == 64) C3X_cluster_length.push_back(cluster_length_count);
      }
      else{
        cluster_spacing++;
        if(i != 0 && C3X_R[i-1] <= 0) cluster_length_count = 0;
        else if (i != 0) cluster_length_count++;
      }
    }

    // count C3Y clusters
    first_cluster = true;
    cluster_spacing = 0;
    cluster_length_count = 0;

    for(int i = 0; i<16; i++){
      if(C3Y_R[i] > 0 && first_cluster){
        cluster_spacing = MWPC_cluster_separation + 1;
        first_cluster = false;
      }

      if(C3Y_R[i] > 0){
        if(cluster_spacing > MWPC_cluster_separation){
          C3Y_clusters++;
          C3Y_cluster_index.push_back(i);
        }
        cluster_length_count++;
        cluster_spacing = 0;

        if(i+2 < 16 && C3Y_R[i+1] <= 0 && C3Y_R[i+2] <= 0) C3Y_cluster_length.push_back(cluster_length_count);
        else if(i + 2 == 16 && C3Y_R[i+1] <= 0) C3Y_cluster_length.push_back(cluster_length_count);
        else if(i + 1 == 16) C3Y_cluster_length.push_back(cluster_length_count);
      }
      else{
        cluster_spacing++;
        if(i != 0 && C3Y_R[i-1] <= 0) cluster_length_count = 0;
        else if ( i != 0) cluster_length_count++;
      }
    }

    // count C4X clusters
    first_cluster = true;
    cluster_spacing = 0;
    cluster_length_count = 0;

    for(int i = 0; i<72; i++){
      if(C4X_R[i] > 0 && first_cluster){
        cluster_spacing = MWPC_cluster_separation + 1;
        first_cluster = false;
      }

      if(C4X_R[i] > 0){
        if(cluster_spacing > MWPC_cluster_separation){
          C4X_clusters++;
          C4X_cluster_index.push_back(i);
        }
        cluster_length_count++;
        cluster_spacing = 0;

        if(i+2 < 72 && C4X_R[i+1] <= 0 && C4X_R[i+2] <= 0) C4X_cluster_length.push_back(cluster_length_count);
        else if(i == 70 && C4X_R[71] <= 0) C4X_cluster_length.push_back(cluster_length_count);
        else if(i == 71) C4X_cluster_length.push_back(cluster_length_count);
      }
      else{
        cluster_spacing++;
        if(i != 0 && C4X_R[i-1] <= 0) cluster_length_count = 0;
        else if ( i != 0) cluster_length_count++;
      }
    }


    // count C4Y clusters
    first_cluster = true;
    cluster_spacing = 0;
    cluster_length_count = 0;

    for(int i = 0; i<16; i++){
      if(C4Y_R[i] > 0 && first_cluster){
        cluster_spacing = MWPC_cluster_separation + 1;
        first_cluster = false;
      }

      if(C4Y_R[i] > 0){
        if(cluster_spacing > MWPC_cluster_separation){
          C4Y_clusters++;
          C4Y_cluster_index.push_back(i);
        }
        cluster_length_count++;
        cluster_spacing = 0;

        if(i+2 < 16 && C4Y_R[i+1] <= 0 && C4Y_R[i+2] <= 0) C4Y_cluster_length.push_back(cluster_length_count);
        else if(i + 2 == 16 && C4Y_R[i+1] <= 0) C4Y_cluster_length.push_back(cluster_length_count);
        else if(i + 1 == 16) C4Y_cluster_length.push_back(cluster_length_count);
      }
      else{
        cluster_spacing++;
        if(i != 0 && C4Y_R[i-1] <= 0) cluster_length_count = 0;
        else if (i != 0) cluster_length_count++;
      }
    }
  }

  double C2_centroid_num = 0;
  double C2_centroid_den = 0;

  if(C2X_clusters == 1){
    if(selected_TOF2>=7){
      for(int i=0; i<56; i++){
        if(C2X_L[i] > 0.){
          C2_centroid_num += C2X_L[i]*C2_ZLoc[i];
          C2_centroid_den += C2X_L[i];
        }
      }
    }
    else{
      for(int i=0; i<56; i++){
        if(C2X_R[i] > 0){
          C2_centroid_num += C2X_R[i]*C2_ZLoc[i];
          C2_centroid_den += C2X_R[i];
        }
      }
    }
    C2X_centroid = double(C2_centroid_num)/C2_centroid_den;
  }

  ////////////////////////////////////////////////////////////////
  /// Draw wire chamber hits.
  // B0 Counter
  // for(int i=0; i<16; i++){
  //   if(tdc_vt48[6][i]>TDC_B0_min && tdc_vt48[6][i]<TDC_B0_max) vec_tdc_b0_6.push_back(tdc_vt48[6][i]);
  //   else vec_tdc_b0_6.push_back(-1);
  //
  //   if(tdc_vt48[7][i]>TDC_B0_min && tdc_vt48[7][i]<TDC_B0_max) vec_tdc_b0_7.push_back(tdc_vt48[7][i]);
  //   else vec_tdc_b0_7.push_back(-1);
  // }
  //
  // sort(vec_tdc_b0_6.begin(), vec_tdc_b0_6.end());
  // sort(vec_tdc_b0_7.begin(), vec_tdc_b0_7.end());


  // ### Display Ck and Cpi infomation

  for(int ic = 0; ic < 14; ic++){
    for(int jc=0; jc<6; jc++){
      if(tdc_ck[ic][jc]>=TDC_Ck_min && tdc_ck[ic][jc]<=TDC_Ck_max) TDC_ck_selected[ic]=tdc_ck[ic][jc];
      if(tdc_cpi[ic][jc]>=TDC_Cpi_min && tdc_cpi[ic][jc]<=TDC_Cpi_max) TDC_cpi_selected[ic]=tdc_cpi[ic][jc];
    }

    if(TDC_ck_selected[ic]>0) vec_Ck.push_back(TDC_ck_selected[ic]);
    if(TDC_cpi_selected[ic]>0) vec_Cpi.push_back(TDC_cpi_selected[ic]);
    TDC_ck_selected[ic]=0;
    TDC_cpi_selected[ic]=0;
  }

  for(unsigned int ik=0; ik<vec_Ck.size(); ik++) TDC_ck_sum += vec_Ck[ik];
  for(unsigned int ip=0; ip<vec_Cpi.size(); ip++) TDC_cpi_sum += vec_Cpi[ip];

  if(vec_Ck.size()>0) TDC_ck_avg = double(TDC_ck_sum)/double(vec_Ck.size());
  else TDC_ck_avg = -1;

  if(vec_Cpi.size()>0) TDC_cpi_avg = double(TDC_cpi_sum)/double(vec_Cpi.size());
  else TDC_cpi_avg = -1;

  for(unsigned int i=0; i<vec_Ck.size(); i++) TDC_ck_sigma2 += pow((vec_Ck[i]-TDC_ck_avg),2);
  for(unsigned int j=0; j<vec_Cpi.size(); j++) TDC_cpi_sigma2 += pow((vec_Cpi[j]-TDC_cpi_avg),2);

  if(vec_Ck.size()>0) TDC_ck_sigma = sqrt(TDC_ck_sigma2/vec_Ck.size());
  else TDC_ck_sigma = -1;

  if(vec_Cpi.size()>0) TDC_cpi_sigma = sqrt(TDC_cpi_sigma2/vec_Cpi.size());
  else TDC_cpi_sigma = -1;

  for(unsigned int i=0; i<vec_Ck.size(); i++){
    if(abs(vec_Ck[i]-TDC_ck_avg) <= 1.4*TDC_ck_sigma){
      TDC_ck_sum2 += vec_Ck[i];
      TDC_ck_counter++;
    }
  }

  for(unsigned int j=0; j<vec_Cpi.size(); j++){
    if(abs(vec_Cpi[j]-TDC_cpi_avg) <= 1.4*TDC_cpi_sigma){
      TDC_cpi_sum2 += vec_Cpi[j];
      TDC_cpi_counter++;
    }
  }

  if(TDC_ck_counter>0) TDC_ck_avg2 = double(TDC_ck_sum2)/double(TDC_ck_counter);
  else TDC_ck_avg2 = -1;

  if(TDC_cpi_counter>0) TDC_cpi_avg2 = double(TDC_cpi_sum2)/double(TDC_cpi_counter);
  else TDC_cpi_avg2 = -1;

  int n_bar_selected;
  n_bar_selected = vec_TARGET_bar_selected.size();
  char TDC_LE_TARGET_sel_string[n_bar_selected][16][20];


  for(unsigned int i=0; i<vec_TARGET_bar_selected.size(); i++){
    for(int j=0; j<16; j++){
      if(tdc_le_target[vec_TARGET_bar_selected[i]][j] > 0){
        sprintf(TDC_LE_TARGET_sel_string[i][j],"%5.1f",550 - (0.91*tdc_le_target[vec_TARGET_bar_selected[i]][j] - 0.625*(tdc_vt48[0][0]-TDC_ck_avg2)));
      }
      else sprintf(TDC_LE_TARGET_sel_string[i][j],"-----");
    }
  }


  // TCanvas *c2;
  // c2 = new TCanvas("Event_Display.C  --  TARGET & SFT","Event_Display.C  --  TARGET & SFT",0,200,1050,700);

  Double_t par0_ADC = 0;  float par0_TDC = 0.;
  Double_t par1_ADC = 0;  float par1_TDC = 0.;

  char angle_string_ADC[100];


  double X_BAR = 999.;
  double Y_BAR = 999.;

  cout << "  " << endl;

  const int gap_to_fit_rotate = gap_to_fit;

  /// Data counters
  int has_data_TDC2 = 0;
  int has_data_ADC2 = 0;
  int has_data_ADC3 = 0;
  int has_data_ADC4 = 0;
  int has_data_ADCA = 0;

  ///// Select edge fiber for track fitting
  double Xloc_gap = 0;
  double Yloc_gap = 0;

  if (gap_to_fit == 0){
    for (int j=0; j<8; j++){
      if ((ADC_High_TARGET[channel[0][j]]>0) && (has_TDC_hit[channel[0][j]])){
        Xloc_gap = Xloc[channel[0][j]];
        Yloc_gap = Yloc[channel[0][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[0][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[0][3]];
  }
  if (gap_to_fit == 1){
    for (int j=0; j<8; j++){
      if ((ADC_High_TARGET[channel[1][j]]>0) && (has_TDC_hit[channel[1][j]])){
        Xloc_gap = Xloc[channel[1][j]];
        Yloc_gap = Yloc[channel[1][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[1][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[1][3]];
  }
  if (gap_to_fit == 2){
    for (int j=0; j<8; j++){
      if ((ADC_High_TARGET[channel[2][j]]>0) && (has_TDC_hit[channel[2][j]])){
        Xloc_gap = Xloc[channel[2][j]];
        Yloc_gap = Yloc[channel[2][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[2][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[2][3]];
  }
  if (gap_to_fit == 3){
    for (int j=0; j<8; j++){
      if ((ADC_High_TARGET[channel[3][j]]>0) && (has_TDC_hit[channel[3][j]])){
        Xloc_gap = Xloc[channel[3][j]];
        Yloc_gap = Yloc[channel[3][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[3][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[3][3]];
  }
  if (gap_to_fit == 4){
    for (int j=0; j<8; j++){
      if ((ADC_High_TARGET[channel[4][j]]>0) && (has_TDC_hit[channel[4][j]])){
        Xloc_gap = Xloc[channel[4][j]];
        Yloc_gap = Yloc[channel[4][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[4][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[4][3]];
  }
  if (gap_to_fit == 5){
    for (int j=0; j<8; j++){
      if ((ADC_High_TARGET[channel[5][j]]>0) && (has_TDC_hit[channel[5][j]])){
        Xloc_gap = Xloc[channel[5][j]];
        Yloc_gap = Yloc[channel[5][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[5][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[5][3]];
  }
  if (gap_to_fit == 6){
    for (int j=0; j<8; j++){
      if ((ADC_High_TARGET[channel[6][j]]>0) && (has_TDC_hit[channel[6][j]])){
        Xloc_gap = Xloc[channel[6][j]];
        Yloc_gap = Yloc[channel[6][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[6][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[6][3]];
  }
  if (gap_to_fit == 7){
    for (int j=0; j<8; j++){
      if ((ADC_High_TARGET[channel[7][j]]>0) && (has_TDC_hit[channel[7][j]])){
        Xloc_gap = Xloc[channel[7][j]];
        Yloc_gap = Yloc[channel[7][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[7][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[7][3]];
  }
  if (gap_to_fit == 8){
    for (int j=0; j<8; j++){
      if ((ADC_High_TARGET[channel[8][j]]>0) && (has_TDC_hit[channel[8][j]])){
        Xloc_gap = Xloc[channel[8][j]];
        Yloc_gap = Yloc[channel[8][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[8][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[8][3]];
  }
  if (gap_to_fit == 9){
    for (int j=0; j<8; j++){
      if ((ADC_High_TARGET[channel[9][j]]>0) && (has_TDC_hit[channel[9][j]])){
        Xloc_gap = Xloc[channel[9][j]];
        Yloc_gap = Yloc[channel[9][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[9][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[9][3]];
  }
  if (gap_to_fit == 10){
    for (int j=0; j<8; j++){
      if ((ADC_High_TARGET[channel[10][j]]>0) && (has_TDC_hit[channel[10][j]])){
        Xloc_gap = Xloc[channel[10][j]];
        Yloc_gap = Yloc[channel[10][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[10][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[10][3]];
  }
  if (gap_to_fit == 11){
    for (int j=0; j<8; j++){
      if ((ADC_High_TARGET[channel[11][j]]>0) && (has_TDC_hit[channel[11][j]])){
        Xloc_gap = Xloc[channel[11][j]];
        Yloc_gap = Yloc[channel[11][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[11][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[11][3]];
  }

  ////Fill TARGET tracking histograms with data
  for(Int_t i=0; i<256; i++){
    if ((i == max_index) || (i == max_index2) || (i == max_index3) || (i == max_index4)) continue;
    if((has_TDC_hit[i]) && (ADC_High_TARGET[i]>Angle_ADC_cut) && (hyp[i] > 12.5)) has_data_ADC2++;
    if((has_TDC_hit[i]) && (ADC_High_TARGET[i]>Angle_ADC_cut) && (adc_low_target[i]-ADC_cut_TARGET2>0)) has_data_TDC2++;
    if((has_TDC_hit[i]) && (ADC_High_TARGET[i]>Angle_ADC_cut) && (adc_low_target[i]-ADC_cut_TARGET2>0)) has_data_ADC3++;
    if((has_TDC_hit[i]) && (ADC_High_TARGET[i]>Angle_ADC_cut) && (hyp[i] < 12.5)) has_data_ADC4++;
  }

  if (has_data_ADC3 > 1){

    double xcoord = 0;
    Int_t unique_x = 0;

    if (has_data_ADCA > 2){

      /// Ring intercept coordinates
      float determinant = 4*(pow(par0_ADC,2))*(pow(par1_ADC,2)) - 4*(pow(par1_ADC,2) + 1)*(pow(par0_ADC,2)-1600);
      float x_circle_int1 = (-2*(par0_ADC)*(par1_ADC) + sqrt(determinant))/(2*(pow(par1_ADC,2) + 1));
      float y_circle_int1 = (par1_ADC)*(x_circle_int1) + par0_ADC;
      float x_circle_int2 = (-2*(par0_ADC)*(par1_ADC) - sqrt(determinant))/(2*(pow(par1_ADC,2) + 1));
      float y_circle_int2 = (par1_ADC)*(x_circle_int2) + par0_ADC;

      if (unique_x == 1){
        x_circle_int1 = xcoord;
        y_circle_int1 = sqrt(1600-(pow(xcoord,2)));
        x_circle_int2 = xcoord;
        y_circle_int2 = sqrt(1600-(pow(xcoord,2)))*-1;
      }

      double SFTxdistance1 = pow((x_circle_int1-Xloc_gap),2);
      double SFTydistance1 = pow((y_circle_int1-Yloc_gap),2);

      double SFTxdistance2 = pow((x_circle_int2-Xloc_gap),2);
      double SFTydistance2 = pow((y_circle_int2-Yloc_gap),2);

      double SFTxhyp1 = double(sqrt(double(SFTxdistance1) + double(SFTydistance1)));
      double SFTxhyp2 = double(sqrt(double(SFTxdistance2) + double(SFTydistance2)));





      if (has_data_ADC4 > 1){

        double x_intercept = double(double(par0_ADC - par0_TDC)/double(par1_TDC - par1_ADC));
        double y_intercept = double(par1_ADC)*double(x_intercept) + par0_ADC;
        double x_distance4[256] = {0};
        double y_distance4[256] = {0};
        int distances4[256];

        for (int q=0; q<256; q++) {
          x_distance4[q] = pow((Xloc[q]-x_intercept),2);
          y_distance4[q] = pow((Yloc[q]-y_intercept),2);
          distances4[q] = double(sqrt(double(x_distance4[q]) + double(y_distance4[q])));
        }

        double min_distance = 10000.0;

        for (int q=0; q<256; q++) {
          x_distance4[q] = pow((Xloc[q]-x_intercept),2);
          y_distance4[q] = pow((Yloc[q]-y_intercept),2);
          distances4[q] = double(sqrt(double(x_distance4[q]) + double(y_distance4[q])));
        }

        for (int q=0; q<256; q++) {
          if (distances4[q] < min_distance) {
            min_distance = distances4[q];
          }
        }
      }
    }
  }
  // else cout << "Histo Fit 4 Is Empty" << endl;

  double a_fit_GoodLG=0.;                     float b_fit_GoodLG=0.;
  double a_fit_GoodLG_weighted=0.;            float b_fit_GoodLG_weighted=0.;
  double a_fit_TDC_selected_weighted=0.;      float b_fit_TDC_selected_weighted=0.;

  for(int i=0; i<12; i++)
  {
    vec_xx_TOF1_Marker.push_back(Gap[i][0][0]);
    vec_xx_TOF1_Marker.push_back(Gap[i][1][0]);
    vec_xx_TOF1_Marker.push_back(Gap[i][2][0]);

    vec_yy_TOF1_Marker.push_back(Gap[i][0][1]);
    vec_yy_TOF1_Marker.push_back(Gap[i][1][1]);
    vec_yy_TOF1_Marker.push_back(Gap[i][2][1]);
  }

  for(int g=0; g<3; g++){
    vec_xx_TOF1.push_back(Gap[gap_to_fit-1][g][0]);
    vec_yy_TOF1.push_back(Gap[gap_to_fit-1][g][1]);
  }

  // ell->SetFillStyle(0);                 ell_Target->SetFillStyle(0);                 ell_L1->SetFillStyle(0);
  // ell->SetLineColor(6);                 ell_Target->SetLineColor(1);                 ell_L1->SetLineColor(4);
  // ell->SetLineWidth(1);                 ell_Target->SetLineWidth(1);                 ell_L1->SetLineWidth(1);
  //
  // Gap1l->SetLineWidth(10);
  // Gap1l->SetLineColor(15);
  //
  // Gap2l->SetLineWidth(10);
  // Gap2l->SetLineColor(15);
  //
  // Gap3l->SetLineWidth(10);
  // Gap3l->SetLineColor(15);
  //
  // Gap4l->SetLineWidth(10);
  // Gap4l->SetLineColor(15);
  //
  // Gap5l->SetLineWidth(10);
  // Gap5l->SetLineColor(15);
  //
  // Gap6l->SetLineWidth(10);
  // Gap6l->SetLineColor(15);
  //
  // Gap7l->SetLineWidth(10);
  // Gap7l->SetLineColor(15);
  //
  // Gap8l->SetLineWidth(10);
  // Gap8l->SetLineColor(15);
  //
  // Gap9l->SetLineWidth(10);
  // Gap9l->SetLineColor(15);
  //
  // Gap10l->SetLineWidth(10);
  // Gap10l->SetLineColor(15);
  //
  // Gap11l->SetLineWidth(10);
  // Gap11l->SetLineColor(15);
  //
  // Gap12l->SetLineWidth(10);
  // Gap12l->SetLineColor(15);
  //
  // // c2->cd(6);

  vec_xx_lepton.clear();            vec_ex_lepton.clear();
  vec_yy_lepton.clear();            vec_ey_lepton.clear();
  vec_xx_lepton_rotate.clear();     vec_ex_lepton_rotate.clear();
  vec_yy_lepton_rotate.clear();     vec_ey_lepton_rotate.clear();
  vec_lepton_size.clear();
  int lepton_counter = 0;

  if(gap_to_fit == 6 || gap_to_fit == 12){
    for(unsigned int i=0; i<lepton_bars.size(); i++){
      if(distance_to_line(Xloc[TARGET_Rotated_index[lepton_bars[i]]],Yloc[TARGET_Rotated_index[lepton_bars[i]]],a_lepton_fit_1,b_lepton_fit_1) <= max_dist
        && TARGET_High_has_neighbours[lepton_bars[i]] && !k_stop_bar[lepton_bars[i]]){
        if(ADC_High_TARGET[lepton_bars[i]]>Angle_ADC_cut && has_TDC_hit[lepton_bars[i]]){

          vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[lepton_bars[i]]]);
          vec_ex_lepton.push_back(TARGET_Errors_X);

          vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[lepton_bars[i]]]);
          vec_ey_lepton.push_back(TARGET_Errors_Y);

          lepton_bars_rotate.push_back(lepton_bars[i]);
          lepton_counter++;

          if(gap_to_fit == 6){
            if(IsIn(TARGET_Rotated_index[lepton_bars[i]],channel[3-1][0], channel[3-1][1],
                      channel[3-1][2], channel[3-1][3],
                      channel[3-1][4], channel[3-1][5],
                      channel[3-1][6], channel[3-1][7])){

              vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[lepton_bars[i]]]);
              vec_ex_lepton.push_back(TARGET_Errors_X);

              vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[lepton_bars[i]]]);
              vec_ey_lepton.push_back(TARGET_Errors_Y);

              lepton_bars_rotate.push_back(lepton_bars[i]);
            }
          }

          if(gap_to_fit == 12){
            if(IsIn(TARGET_Rotated_index[lepton_bars[i]],channel[9-1][0], channel[9-1][1],
                      channel[9-1][2], channel[9-1][3],
                      channel[9-1][4], channel[9-1][5],
                      channel[9-1][6], channel[9-1][7])){

              vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[lepton_bars[i]]]);
              vec_ex_lepton.push_back(TARGET_Errors_X);

              vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[lepton_bars[i]]]);
              vec_ey_lepton.push_back(TARGET_Errors_Y);

              lepton_bars_rotate.push_back(lepton_bars[i]);
            }
          }
        }

        if(ADC_High_TARGET[lepton_bars[i]]<0 && ADC_Low_TARGET[lepton_bars[i]]>0 && Switch==1){
          vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[lepton_bars[i]]]);
          vec_ex_lepton.push_back(TARGET_Errors_X);

          vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[lepton_bars[i]]]);
          vec_ey_lepton.push_back(TARGET_Errors_Y);

          lepton_bars_rotate.push_back(lepton_bars[i]);
          lepton_counter++;
        }
      }
    }
  }
  else{
    for(unsigned int i=0; i<lepton_bars.size(); i++){
      if(distance_to_line(Xloc[lepton_bars[i]],Yloc[lepton_bars[i]],a_lepton_fit_1,b_lepton_fit_1) <= max_dist
        && TARGET_High_has_neighbours[lepton_bars[i]] && !k_stop_bar[lepton_bars[i]]){
        if(ADC_High_TARGET[lepton_bars[i]]>Angle_ADC_cut && has_TDC_hit[lepton_bars[i]]){

          vec_xx_lepton.push_back(Xloc[lepton_bars[i]]);
          vec_ex_lepton.push_back(TARGET_Errors_X);

          vec_yy_lepton.push_back(Yloc[lepton_bars[i]]);
          vec_ey_lepton.push_back(TARGET_Errors_Y);

          lepton_counter++;

          if(IsIn(lepton_bars[i],channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],
                    channel[gap_to_fit-1][2], channel[gap_to_fit-1][3],
                    channel[gap_to_fit-1][4], channel[gap_to_fit-1][5],
                    channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){

            vec_xx_lepton.push_back(Xloc[lepton_bars[i]]);
            vec_ex_lepton.push_back(TARGET_Errors_X);

            vec_yy_lepton.push_back(Yloc[lepton_bars[i]]);
            vec_ey_lepton.push_back(TARGET_Errors_Y);
          }
        }

        if(ADC_High_TARGET[lepton_bars[i]]<0 && ADC_Low_TARGET[lepton_bars[i]]>0 && Switch==1){
          vec_xx_lepton.push_back(Xloc[lepton_bars[i]]);
          vec_ex_lepton.push_back(TARGET_Errors_X);

          vec_yy_lepton.push_back(Yloc[lepton_bars[i]]);
          vec_ey_lepton.push_back(TARGET_Errors_Y);

          lepton_counter++;
        }
      }
    }
  }

  double a_lepton_fit_2 = 0.;
  double b_lepton_fit_2 = 0.;
  double Chis_lepton_fit_2 = 0.;
  int ndf_lepton_fit_2 = 0;

  // TGraph *gr_lepton_2;
  // TGraph gr_lepton_2 = TGraphErrors(vec_xx_lepton.size(),&vec_xx_lepton[0],&vec_yy_lepton[0],
                                //  &vec_ex_lepton[0],&vec_ey_lepton[0]);
  TGraph gr_lepton_2 = TGraph(vec_xx_lepton.size(),&vec_xx_lepton[0],&vec_yy_lepton[0]);


  // TF1 *func_lepton_fit_2 = new TF1("lepton_fit_2", "pol1");
  TF1 func_lepton_fit_2 = TF1("lepton_fit_2", "pol1");

  // if(vec_xx_lepton.size()>0){
  //
  //   gr_lepton_2->SetMarkerStyle(21);
  //   gr_lepton_2->SetMarkerColor(2);
  //   gr_lepton_2->SetMarkerSize(0.8);
  //   gr_lepton_2->GetXaxis()->SetLimits(-50.,50.);
  //   gr_lepton_2->GetYaxis()->SetRangeUser(-50.,50.);
  //
  //   if(gap_to_fit==1 || gap_to_fit==2 || gap_to_fit==7 || gap_to_fit==8){
  //     func_lepton_fit_2->SetParameter(0,0);
  //     func_lepton_fit_2->SetParameter(1,1);
  //   }
  //   if(gap_to_fit==4 || gap_to_fit==5 || gap_to_fit==10 || gap_to_fit==11){
  //     func_lepton_fit_2->SetParameter(0,0);
  //     func_lepton_fit_2->SetParameter(1,-1);
  //   }
  //   else{
  //     func_lepton_fit_2->SetParameter(0,0);
  //     func_lepton_fit_2->SetParameter(1,1);
  //   }
  //
  //
  //   func_lepton_fit_2->SetParLimits(0,-50,50);
  //   func_lepton_fit_2->SetParLimits(1,-50,50);
  if(vec_xx_lepton.size()>0){

    gr_lepton_2.SetMarkerStyle(21);
    gr_lepton_2.SetMarkerColor(2);
    gr_lepton_2.SetMarkerSize(0.8);
    gr_lepton_2.GetXaxis()->SetLimits(-50.,50.);
    gr_lepton_2.GetYaxis()->SetRangeUser(-50.,50.);

    if(gap_to_fit==1 || gap_to_fit==2 || gap_to_fit==7 || gap_to_fit==8){
      func_lepton_fit_2.SetParameter(0,0);
      func_lepton_fit_2.SetParameter(1,1);
    }
    if(gap_to_fit==4 || gap_to_fit==5 || gap_to_fit==10 || gap_to_fit==11){
      func_lepton_fit_2.SetParameter(0,0);
      func_lepton_fit_2.SetParameter(1,-1);
    }
    else{
      func_lepton_fit_2.SetParameter(0,0);
      func_lepton_fit_2.SetParameter(1,1);
    }


    func_lepton_fit_2.SetParLimits(0,-50,50);
    func_lepton_fit_2.SetParLimits(1,-50,50);


    gr_lepton_2.Fit("lepton_fit_2","QB");
    // delete func_lepton_fit_2;
    // func_lepton_fit_2 = gr_lepton_2->GetFunction("lepton_fit_2");
    // func_lepton_fit_2->SetLineWidth(2);
    // func_lepton_fit_2->SetLineColor(2);
    //
    // a_lepton_fit_2 = func_lepton_fit_2->GetParameter(1);
    // b_lepton_fit_2 = func_lepton_fit_2->GetParameter(0);
    // Chis_lepton_fit_2 = func_lepton_fit_2->GetChisquare();
    // ndf_lepton_fit_2 = func_lepton_fit_2->GetNDF();
    TF1 *func_lepton_fit_2_ptr = gr_lepton_2.GetFunction("lepton_fit_2");
    func_lepton_fit_2_ptr->SetLineWidth(2);
    func_lepton_fit_2_ptr->SetLineColor(2);

    a_lepton_fit_2 = func_lepton_fit_2_ptr->GetParameter(1);
    b_lepton_fit_2 = func_lepton_fit_2_ptr->GetParameter(0);
    Chis_lepton_fit_2 = func_lepton_fit_2_ptr->GetChisquare();
    ndf_lepton_fit_2 = func_lepton_fit_2_ptr->GetNDF();
  }
  // else cout << "Empty Fit 2" << endl;

  // TGraph gr2_Leptons_rotate = TGraphErrors(vec_xx_lepton_rotate.size(),&vec_xx_lepton_rotate[0],&vec_yy_lepton_rotate[0],
  //                                               &vec_ex_lepton[0],&vec_ey_lepton[0]);
  // TGraph gr2_Leptons_rotate = TGraphErrors(vec_xx_lepton_rotate.size(),&vec_xx_lepton_rotate[0],&vec_yy_lepton_rotate[0]);
  //
  // gr2_Leptons_rotate.SetMarkerStyle(25);
  // gr2_Leptons_rotate.SetMarkerColor(2);
  // gr2_Leptons_rotate.SetMarkerSize(0.8);
  // gr2_Leptons_rotate.GetXaxis()->SetLimits(-50.,50.);
  // gr2_Leptons_rotate.GetYaxis()->SetRangeUser(-50.,50.);

  float x_int_TDC_selected_weighted[2];         float y_int_TDC_selected_weighted[2];
  float x_int_GoodLG[2];                        float y_int_GoodLG[2];
  float x_int_GoodLG_weighted[2];               float y_int_GoodLG_weighted[2];
  float x_int_TDC_Gap_Fibers[2];                float y_int_TDC_Gap_Fibers[2];
  float x_int_TDC_Gap_Fibers_SFT[2];            float y_int_TDC_Gap_Fibers_SFT[2];
  float x_int_TARGET[2];                        float y_int_TARGET[2];

  x_int_GoodLG[0] = intersectx1(a_fit_GoodLG, b_fit_GoodLG, R_TOF1);
  x_int_GoodLG[1] = intersectx2(a_fit_GoodLG, b_fit_GoodLG, R_TOF1);
  y_int_GoodLG[0] = y1_int(x_int_GoodLG[0], a_fit_GoodLG, b_fit_GoodLG);
  y_int_GoodLG[1] = y2_int(x_int_GoodLG[1], a_fit_GoodLG, b_fit_GoodLG);

  x_int_TDC_Gap_Fibers[0] = intersectx1(a_lepton_fit_2, b_lepton_fit_2, R_TOF1);
  x_int_TDC_Gap_Fibers[1] = intersectx2(a_lepton_fit_2, b_lepton_fit_2, R_TOF1);
  y_int_TDC_Gap_Fibers[0] = y1_int(x_int_TDC_Gap_Fibers[0], a_lepton_fit_2, b_lepton_fit_2);
  y_int_TDC_Gap_Fibers[1] = y2_int(x_int_TDC_Gap_Fibers[1], a_lepton_fit_2, b_lepton_fit_2);

  x_int_TDC_Gap_Fibers_SFT[0] = intersectx1(a_lepton_fit_2, b_lepton_fit_2, R_SFT_L1);
  x_int_TDC_Gap_Fibers_SFT[1] = intersectx2(a_lepton_fit_2, b_lepton_fit_2, R_SFT_L1);
  y_int_TDC_Gap_Fibers_SFT[0] = y1_int(x_int_TDC_Gap_Fibers_SFT[0], a_lepton_fit_2, b_lepton_fit_2);
  y_int_TDC_Gap_Fibers_SFT[1] = y2_int(x_int_TDC_Gap_Fibers_SFT[1], a_lepton_fit_2, b_lepton_fit_2);


  /// Selection of the Good Intersect Coordinates
  //////////////////////////////////////////////////////////////////////
  float x_GoodLG_intersect1=0.;    float y_GoodLG_intersect1=0.;
  float x_TDC_Gap_Fibers=0.;       float y_TDC_Gap_Fibers=0.;

  float dist1_GoodLG[2];
  float dist1_TDC_Gap_Fibers[2];
  float dist1_TDC_Gap_Fibers_SFT[2];

  for(int i=0; i<2; i++){

    if(gap_to_fit == 6){
      dist1_GoodLG[i] = distance(x_int_GoodLG[i], y_int_GoodLG[i], Gap[3-1][1][0], Gap[3-1][1][1]);
      dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[3-1][1][0], Gap[3-1][1][1]);
      dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[3-1][1][0], Gap[3-1][1][1]);
    }
    else if(gap_to_fit == 12){
      dist1_GoodLG[i] = distance(x_int_GoodLG[i], y_int_GoodLG[i], Gap[9-1][1][0], Gap[9-1][1][1]);
      dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[9-1][1][0], Gap[9-1][1][1]);
      dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[9-1][1][0], Gap[9-1][1][1]);
    }
    else{
      dist1_GoodLG[i] = distance(x_int_GoodLG[i], y_int_GoodLG[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
    }
  }

  if(dist1_GoodLG[0] < dist1_GoodLG[1]){
    x_GoodLG_intersect1 = x_int_GoodLG[0];
    y_GoodLG_intersect1 = y_int_GoodLG[0];
  }
  else if(dist1_GoodLG[1] < dist1_GoodLG[0]){
    x_GoodLG_intersect1 = x_int_GoodLG[1];
    y_GoodLG_intersect1 = y_int_GoodLG[1];
  }
  else cout << "ERROR !" << endl;

  if(dist1_TDC_Gap_Fibers[0] < dist1_TDC_Gap_Fibers[1]){
    x_TDC_Gap_Fibers = x_int_TDC_Gap_Fibers[0];
    y_TDC_Gap_Fibers = y_int_TDC_Gap_Fibers[0];
  }
  else if(dist1_TDC_Gap_Fibers[1] < dist1_TDC_Gap_Fibers[0]){
    x_TDC_Gap_Fibers = x_int_TDC_Gap_Fibers[1];
    y_TDC_Gap_Fibers = y_int_TDC_Gap_Fibers[1];
  }
  else cout << "ERROR !" << endl;

  /// Selection of the Good TOF1 Section
  //////////////////////////////////////////////////////////////////////
  float dist2_TDC_selected[3];
  float dist2_GoodLG[3];

  float dist2_TDC_selected_min = 1000.;
  float dist2_GoodLG_min = 1000.;
  int selected_TDC_selected = 0;

  for(int ii=0; ii<3; ii++){

    if(gap_to_fit == 6){
      dist2_TDC_selected[ii] = distance(x_TDC_Gap_Fibers, y_TDC_Gap_Fibers, Gap[3-1][ii][0], Gap[3-1][ii][1]);
      dist2_GoodLG[ii] = distance(x_GoodLG_intersect1, y_GoodLG_intersect1, Gap[3-1][ii][0], Gap[3-1][ii][1]);
    }
    else if(gap_to_fit == 12){
      dist2_TDC_selected[ii] = distance(x_TDC_Gap_Fibers, y_TDC_Gap_Fibers, Gap[9-1][ii][0], Gap[9-1][ii][1]);
      dist2_GoodLG[ii] = distance(x_GoodLG_intersect1, y_GoodLG_intersect1, Gap[9-1][ii][0], Gap[9-1][ii][1]);
    }
    else{
      dist2_TDC_selected[ii] = distance(x_TDC_Gap_Fibers, y_TDC_Gap_Fibers, Gap[gap_to_fit-1][ii][0], Gap[gap_to_fit-1][ii][1]);
      dist2_GoodLG[ii] = distance(x_GoodLG_intersect1, y_GoodLG_intersect1, Gap[gap_to_fit-1][ii][0], Gap[gap_to_fit-1][ii][1]);
    }

    if(dist2_TDC_selected[ii] <= dist2_TDC_selected_min){
      dist2_TDC_selected_min = dist2_TDC_selected[ii];
      selected_TDC_selected = ii;
    }

    if(dist2_GoodLG[ii] <= dist2_GoodLG_min) dist2_GoodLG_min = dist2_GoodLG[ii];
  }

  vec_xx_TOF1_closest.push_back(Gap[gap_to_fit_rotate-1][selected_TDC_selected][0]);
  vec_yy_TOF1_closest.push_back(Gap[gap_to_fit_rotate-1][selected_TDC_selected][1]);

  for(int i = 0; i<3; i++){

    if(gap_to_fit == 6){
      vec_xx_lepton.push_back(Gap[3-1][selected_TDC_selected][0]);
      vec_ex_lepton.push_back(TOF1_Errors_X[3-1][selected_TDC_selected]);

      vec_yy_lepton.push_back(Gap[3-1][selected_TDC_selected][1]);
      vec_ey_lepton.push_back(TOF1_Errors_Y[3-1][selected_TDC_selected]);
    }
    else if(gap_to_fit == 12){
      vec_xx_lepton.push_back(Gap[9-1][selected_TDC_selected][0]);
      vec_ex_lepton.push_back(TOF1_Errors_X[9-1][selected_TDC_selected]);

      vec_yy_lepton.push_back(Gap[9-1][selected_TDC_selected][1]);
      vec_ey_lepton.push_back(TOF1_Errors_Y[9-1][selected_TDC_selected]);
    }
    else{
    vec_xx_lepton.push_back(Gap[gap_to_fit-1][selected_TDC_selected][0]);
    vec_ex_lepton.push_back(TOF1_Errors_X[gap_to_fit-1][selected_TDC_selected]);

    vec_yy_lepton.push_back(Gap[gap_to_fit-1][selected_TDC_selected][1]);
    vec_ey_lepton.push_back(TOF1_Errors_Y[gap_to_fit-1][selected_TDC_selected]);
    }
  }

  double ParError = 999.99;
  double ChiS = 0.0;
  int ndf = 0;

  // TGraph gr_lepton_3 = TGraphErrors(vec_xx_lepton.size(),&vec_xx_lepton[0],&vec_yy_lepton[0],
  //                                        &vec_ex_lepton[0],&vec_ey_lepton[0]);
  TGraph gr_lepton_3 = TGraph(vec_xx_lepton.size(),&vec_xx_lepton[0],&vec_yy_lepton[0]);

  gr_lepton_3.SetMarkerStyle(21);
  gr_lepton_3.SetMarkerColor(2);
  gr_lepton_3.SetMarkerSize(0.8);
  gr_lepton_3.GetXaxis()->SetLimits(-50.,50.);
  gr_lepton_3.GetYaxis()->SetRangeUser(-50.,50.);

  double a_lepton_fit_3 = 0.;
  double b_lepton_fit_3 = 0.;
  double Chis_lepton_fit_3 = 0.;
  int ndf_lepton_fit_3 = 0;

  // TF1 *func_lepton_fit_3 = new TF1("lepton_fit_3", "pol1");
  TF1 func_lepton_fit_3 = TF1("lepton_fit_3", "pol1");

  // if(gap_to_fit==1 || gap_to_fit==2 || gap_to_fit==7 || gap_to_fit==8){
  //   func_lepton_fit_3->SetParameter(0,0);
  //   func_lepton_fit_3->SetParameter(1,1);
  // }
  // if(gap_to_fit==4 || gap_to_fit==5 || gap_to_fit==10 || gap_to_fit==11){
  //   func_lepton_fit_3->SetParameter(0,0);
  //   func_lepton_fit_3->SetParameter(1,-1);
  // }
  // else{
  //     func_lepton_fit_3->SetParameter(0,0);
  //     func_lepton_fit_3->SetParameter(1,1);
  //   }
  //
  // func_lepton_fit_3->SetParLimits(0,-50,50);
  // func_lepton_fit_3->SetParLimits(1,-50,50);
  if(gap_to_fit==1 || gap_to_fit==2 || gap_to_fit==7 || gap_to_fit==8){
    func_lepton_fit_3.SetParameter(0,0);
    func_lepton_fit_3.SetParameter(1,1);
  }
  if(gap_to_fit==4 || gap_to_fit==5 || gap_to_fit==10 || gap_to_fit==11){
    func_lepton_fit_3.SetParameter(0,0);
    func_lepton_fit_3.SetParameter(1,-1);
  }
  else{
      func_lepton_fit_3.SetParameter(0,0);
      func_lepton_fit_3.SetParameter(1,1);
    }

  func_lepton_fit_3.SetParLimits(0,-50,50);
  func_lepton_fit_3.SetParLimits(1,-50,50);

  gr_lepton_3.Fit("lepton_fit_3","QB");



  // delete func_lepton_fit_3;
  // func_lepton_fit_3 = gr_lepton_3->GetFunction("lepton_fit_3");
  // func_lepton_fit_3->SetLineWidth(2);
  // func_lepton_fit_3->SetLineColor(2);
  //
  //
  // a_lepton_fit_3 = func_lepton_fit_3->GetParameter(1);
  // b_lepton_fit_3 = func_lepton_fit_3->GetParameter(0);
  //
  // ParError = func_lepton_fit_3->GetParError(1);
  // Chis_lepton_fit_3 = func_lepton_fit_3->GetChisquare();
  // ndf_lepton_fit_3 = func_lepton_fit_3->GetNDF();
  TF1 *func_lepton_fit_3_ptr = gr_lepton_3.GetFunction("lepton_fit_3");
  func_lepton_fit_3_ptr->SetLineWidth(2);
  func_lepton_fit_3_ptr->SetLineColor(2);


  a_lepton_fit_3 = func_lepton_fit_3_ptr->GetParameter(1);
  b_lepton_fit_3 = func_lepton_fit_3_ptr->GetParameter(0);

  ParError = func_lepton_fit_3_ptr->GetParError(1);
  Chis_lepton_fit_3 = func_lepton_fit_3_ptr->GetChisquare();
  ndf_lepton_fit_3 = func_lepton_fit_3_ptr->GetNDF();

  if(gap_to_fit==6 || gap_to_fit==12){
    for(unsigned int j=0; j<lepton_bars_rotate.size(); j++){

      vec_xx_lepton_rotate.push_back(Xloc[lepton_bars_rotate[j]]);
      vec_ex_lepton_rotate.push_back(TARGET_Errors_X);

      vec_yy_lepton_rotate.push_back(Yloc[lepton_bars_rotate[j]]);
      vec_ey_lepton_rotate.push_back(TARGET_Errors_Y);

    }
  }



  // TGraph gr3_Leptons_rotate = TGraphErrors(vec_xx_lepton_rotate.size(),&vec_xx_lepton_rotate[0],&vec_yy_lepton_rotate[0],
  //                                               &vec_ex_lepton_rotate[0],&vec_ey_lepton_rotate[0]);
  TGraph gr3_Leptons_rotate = TGraph(vec_xx_lepton_rotate.size(),&vec_xx_lepton_rotate[0],&vec_yy_lepton_rotate[0]);

  gr3_Leptons_rotate.SetMarkerStyle(21);
  gr3_Leptons_rotate.SetMarkerColor(2);
  gr3_Leptons_rotate.SetMarkerSize(0.8);
  gr3_Leptons_rotate.GetXaxis()->SetLimits(-50.,50.);
  gr3_Leptons_rotate.GetYaxis()->SetRangeUser(-50.,50.);

  x_int_TDC_Gap_Fibers[0] = intersectx1(a_lepton_fit_3, b_lepton_fit_3, R_TOF1);
  x_int_TDC_Gap_Fibers[1] = intersectx2(a_lepton_fit_3, b_lepton_fit_3, R_TOF1);
  x_int_TDC_Gap_Fibers_SFT[0] = intersectx1(a_lepton_fit_3, b_lepton_fit_3, R_SFT_L1);
  x_int_TDC_Gap_Fibers_SFT[1] = intersectx2(a_lepton_fit_3, b_lepton_fit_3, R_SFT_L1);

  y_int_TDC_Gap_Fibers[0] = y1_int(x_int_TDC_Gap_Fibers[0], a_lepton_fit_3, b_lepton_fit_3);
  y_int_TDC_Gap_Fibers[1] = y2_int(x_int_TDC_Gap_Fibers[1], a_lepton_fit_3, b_lepton_fit_3);
  y_int_TDC_Gap_Fibers_SFT[0] = y1_int(x_int_TDC_Gap_Fibers_SFT[0], a_lepton_fit_3, b_lepton_fit_3);
  y_int_TDC_Gap_Fibers_SFT[1] = y2_int(x_int_TDC_Gap_Fibers_SFT[1], a_lepton_fit_3, b_lepton_fit_3);

  x_int_TARGET[0] = intersectx1(a_lepton_fit_3, b_lepton_fit_3, R_TARGET);
  x_int_TARGET[1] = intersectx2(a_lepton_fit_3, b_lepton_fit_3, R_TARGET);
  y_int_TARGET[0] = y1_int(x_int_TARGET[0], a_lepton_fit_3, b_lepton_fit_3);
  y_int_TARGET[1] = y2_int(x_int_TARGET[1], a_lepton_fit_3, b_lepton_fit_3);

  // target_uncertainty = find_target_uncertainty(x_int_TARGET, a_lepton_fit_3, b_lepton_fit_3, );
  float x_TDC_Gap_Fibers_intersect1=0.;          float y_TDC_Gap_Fibers_intersect1=0.;
  float x_TDC_Gap_Fibers_SFT_intersect1=0.;      float y_TDC_Gap_Fibers_SFT_intersect1=0.;
  float x_TARGET_intersect=0;                    float y_TARGET_intersect=0;
  float x_Arrows=0;                               float y_Arrows=0;

  float dist1_TARGET_intersect[2];

  for(int i=0; i<2; i++){
    if(gap_to_fit==6){
      dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[3-1][1][0], Gap[3-1][1][1]);
      dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[3-1][1][0], Gap[3-1][1][1]);
      dist1_TARGET_intersect[i] = distance(x_int_TARGET[i], y_int_TARGET[i], Gap[3-1][1][0], Gap[3-1][1][1]);
    }
    else if(gap_to_fit==12){
      dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[9-1][1][0], Gap[9-1][1][1]);
      dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[9-1][1][0], Gap[9-1][1][1]);
      dist1_TARGET_intersect[i] = distance(x_int_TARGET[i], y_int_TARGET[i], Gap[9-1][1][0], Gap[9-1][1][1]);
    }
    else{
      dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TARGET_intersect[i] = distance(x_int_TARGET[i], y_int_TARGET[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
    }
  }

  if(dist1_TDC_Gap_Fibers[0] < dist1_TDC_Gap_Fibers[1]){
    x_TDC_Gap_Fibers_intersect1 = x_int_TDC_Gap_Fibers[0];
    y_TDC_Gap_Fibers_intersect1 = y_int_TDC_Gap_Fibers[0];
  }
  else if(dist1_TDC_Gap_Fibers[1] < dist1_TDC_Gap_Fibers[0]){
    x_TDC_Gap_Fibers_intersect1 = x_int_TDC_Gap_Fibers[1];
    y_TDC_Gap_Fibers_intersect1 = y_int_TDC_Gap_Fibers[1];
  }
  else cout << "ERROR !" << endl;

  if(dist1_TDC_Gap_Fibers_SFT[0] < dist1_TDC_Gap_Fibers_SFT[1]){
    x_TDC_Gap_Fibers_SFT_intersect1 = x_int_TDC_Gap_Fibers_SFT[0];
    y_TDC_Gap_Fibers_SFT_intersect1 = y_int_TDC_Gap_Fibers_SFT[0];
  }
  else if(dist1_TDC_Gap_Fibers_SFT[1] < dist1_TDC_Gap_Fibers_SFT[0]){
    x_TDC_Gap_Fibers_SFT_intersect1 = x_int_TDC_Gap_Fibers_SFT[1];
    y_TDC_Gap_Fibers_SFT_intersect1 = y_int_TDC_Gap_Fibers_SFT[1];
  }
  else cout << "ERROR !" << endl;

  if(dist1_TARGET_intersect[0] < dist1_TARGET_intersect[1]){
    x_TARGET_intersect = x_int_TARGET[0];
    y_TARGET_intersect = y_int_TARGET[0];
    x_Arrows = x_int_TARGET[1];
    y_Arrows = y_int_TARGET[1];
  }
  else if(dist1_TARGET_intersect[1] < dist1_TARGET_intersect[0]){
    x_TARGET_intersect = x_int_TARGET[1];
    y_TARGET_intersect = y_int_TARGET[1];
    x_Arrows = x_int_TARGET[0];
    y_Arrows = y_int_TARGET[0];
  }
  else cout << "ERROR !" << endl;

  int Axis_Vector_Length = 10;
  // TArrow *x_guide;
  // TArrow *y_guide;
  // TArrow *x_guide_rotate;
  // TArrow *y_guide_rotate;

  //if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
  if(gap_to_fit==6 || gap_to_fit==12){

    vec_xx_int_TDC_Gap_Fibers.push_back(y_TDC_Gap_Fibers_intersect1);
    vec_yy_int_TDC_Gap_Fibers.push_back(-x_TDC_Gap_Fibers_intersect1);

    vec_xx_int_TDC_Gap_Fibers_SFT.push_back(y_TDC_Gap_Fibers_SFT_intersect1);
    vec_yy_int_TDC_Gap_Fibers_SFT.push_back(-x_TDC_Gap_Fibers_SFT_intersect1);

    vec_xx_int_TDC_TARGET.push_back(y_TARGET_intersect);
    vec_yy_int_TDC_TARGET.push_back(-x_TARGET_intersect);

    // TArrow *x_guide = new TArrow(y_Arrows,-x_Arrows,y_Arrows + Axis_Vector_Length,-x_Arrows, 0.005, "|>");
    // TArrow *y_guide = new TArrow(y_Arrows,-x_Arrows,y_Arrows,-x_Arrows + Axis_Vector_Length, 0.005, "|>");
  }
  else{
    vec_xx_int_TDC_Gap_Fibers.push_back(x_TDC_Gap_Fibers_intersect1);
    vec_yy_int_TDC_Gap_Fibers.push_back(y_TDC_Gap_Fibers_intersect1);

    vec_xx_int_TDC_Gap_Fibers_SFT.push_back(x_TDC_Gap_Fibers_SFT_intersect1);
    vec_yy_int_TDC_Gap_Fibers_SFT.push_back(y_TDC_Gap_Fibers_SFT_intersect1);

    vec_xx_int_TDC_TARGET.push_back(x_TARGET_intersect);
    vec_yy_int_TDC_TARGET.push_back(y_TARGET_intersect);

    // TArrow *x_guide = new TArrow(x_Arrows, y_Arrows, x_Arrows + Axis_Vector_Length, y_Arrows, 0.005, "|>");
    // TArrow *y_guide = new TArrow(x_Arrows, y_Arrows, x_Arrows, y_Arrows + Axis_Vector_Length, 0.005, "|>");
  }

  /// Distance from K-stop bar to best fit line
  double dist_to_k_stop = 0;

  if(X_BAR != -10000){
    if((gap_to_fit_rotate==6 || gap_to_fit_rotate==12 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
      dist_to_k_stop = distance_to_line(Y_BAR,-X_BAR,a_lepton_fit_3,b_lepton_fit_3);
    }
    else{
      dist_to_k_stop = distance_to_line(X_BAR,Y_BAR,a_lepton_fit_3,b_lepton_fit_3);
    }
  }
  else{
    dist_to_k_stop = -10000;
  }

  /// Angle Calculation
  float a_final_guide = 0.;
  float alpha_guide = 0.;
  float tanalpha_guide = 0.;
  float angle_final_guide = 0.;

  if(gap_to_fit==6 || gap_to_fit==12) a_final_guide = a_lepton_fit_3;
  else a_final_guide = a_lepton_fit_3;

  tanalpha_guide = a_final_guide;
  alpha_guide = atan(tanalpha_guide);

  /// Determination of Delta Phi
  double Delta_phi = 999.99;  double Delta_phi_deg = 999.99;
  Delta_phi = sin(alpha_guide)*cos(alpha_guide)*(ParError/a_final_guide);

  Delta_phi_deg = (180/PI)*Delta_phi;

  // x_guide->SetLineWidth(2);     y_guide->SetLineWidth(2);
  // x_guide->SetLineColor(4);     y_guide->SetLineColor(4);

  // x_guide_rotate = new TArrow(y_TARGET_intersect, -x_TARGET_intersect, y_TARGET_intersect + Axis_Vector_Length, -x_TARGET_intersect, 0.005, "|>");
  // y_guide_rotate = new TArrow(y_TARGET_intersect, -x_TARGET_intersect, y_TARGET_intersect, -x_TARGET_intersect + Axis_Vector_Length, 0.005, "|>");

  double x_exit_rotate;
  double y_exit_rotate;

  x_exit_rotate = -y_TARGET_intersect;
  y_exit_rotate = x_TARGET_intersect;

  // x_guide_rotate->SetLineWidth(2);     y_guide_rotate->SetLineWidth(2);
  // x_guide_rotate->SetLineColor(4);     y_guide_rotate->SetLineColor(4);

  if(gap_to_fit==6 || gap_to_fit==12){
    alpha_guide = alpha_guide*(180/PI) + 90;

    if((-x_TDC_Gap_Fibers_intersect1 + x_TARGET_intersect) < 0.) angle_final_guide = 180. + alpha_guide;
    else{

      if((y_TDC_Gap_Fibers_intersect1 - y_TARGET_intersect) >= 0) angle_final_guide = alpha_guide;
      else angle_final_guide = alpha_guide + 360.;
    }
  }
  else{

    if((x_TDC_Gap_Fibers_intersect1 - x_TARGET_intersect) < 0) angle_final_guide = 180. + alpha_guide * (180./PI);
    else{

      if((y_TDC_Gap_Fibers_intersect1 - y_TARGET_intersect) >= 0) angle_final_guide = alpha_guide * (180./PI);
      else angle_final_guide = alpha_guide * (180./PI) + 360.;
    }
  }




  double a_fit_lepton_rotate = 999.;
  double b_fit_lepton_rotate = 999.;
  // cout << "" << endl;


  a_fit_lepton_rotate = -1/a_lepton_fit_3;
  b_fit_lepton_rotate = b_lepton_fit_3/a_lepton_fit_3;





  // if(X_BAR != -10000){
  //   double x_tof1_intersect_1 = 0;
  //   double y_tof1_intersect_1 = 0;
  //   double x_tof1_intersect_2 = 0;
  //   double y_tof1_intersect_2 = 0;
  //   double x_tof1_intersect = 0;
  //   double y_tof1_intersect = 0;
  //
  //   double alpha = angle_final_guide - 90.0;
  //
  //   if((gap_to_fit_rotate==6 || gap_to_fit_rotate==12 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
  //     x_tof1_intersect_1 = y_TDC_Gap_Fibers_intersect1 - dist_to_k_stop*cos(alpha*M_PI/180.0);
  //     x_tof1_intersect_2 = y_TDC_Gap_Fibers_intersect1 + dist_to_k_stop*cos(alpha*M_PI/180.0);
  //     y_tof1_intersect_1 = -x_TDC_Gap_Fibers_intersect1 - dist_to_k_stop*sin(alpha*M_PI/180.0);
  //     y_tof1_intersect_2 = -x_TDC_Gap_Fibers_intersect1 + dist_to_k_stop*sin(alpha*M_PI/180.0);
  //
  //     if(distance(X_BAR,Y_BAR,x_tof1_intersect_1,y_tof1_intersect_1) < distance(X_BAR,Y_BAR,x_tof1_intersect_2,y_tof1_intersect_2)){
  //       x_tof1_intersect = x_tof1_intersect_1;
  //       y_tof1_intersect = y_tof1_intersect_1;
  //     }
  //     else{
  //       x_tof1_intersect = x_tof1_intersect_2;
  //       y_tof1_intersect = y_tof1_intersect_2;
  //     }
  //   }
  //   else{
  //     x_tof1_intersect_1 = x_TDC_Gap_Fibers_intersect1 - dist_to_k_stop*cos(alpha*M_PI/180.0);
  //     x_tof1_intersect_2 = x_TDC_Gap_Fibers_intersect1 + dist_to_k_stop*cos(alpha*M_PI/180.0);
  //     y_tof1_intersect_1 = y_TDC_Gap_Fibers_intersect1 - dist_to_k_stop*sin(alpha*M_PI/180.0);
  //     y_tof1_intersect_2 = y_TDC_Gap_Fibers_intersect1 + dist_to_k_stop*sin(alpha*M_PI/180.0);
  //
  //     if(distance(X_BAR,Y_BAR,x_tof1_intersect_1,y_tof1_intersect_1) < distance(X_BAR,Y_BAR,x_tof1_intersect_2,y_tof1_intersect_2)){
  //       x_tof1_intersect = x_tof1_intersect_1;
  //       y_tof1_intersect = y_tof1_intersect_1;
  //     }
  //     else{
  //       x_tof1_intersect = x_tof1_intersect_2;
  //       y_tof1_intersect = y_tof1_intersect_2;
  //     }
  //   }
  // }
  // else{
  //   cout << "No K-Stop for track length." << endl;
  // }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Calculation of XBAR and YBAR

  double X_weights = 0.;
  double Y_weights = 0.;
  double total_energy = 0.;
  vector<double> vec_kaon_centroid_coordinates;
  vector<double> vec_fit_lines_intersect;
  vector<double> vec_k_stop_coordinates;
  vec_kaon_centroid_coordinates.clear();
  vec_fit_lines_intersect.clear();
  vec_k_stop_coordinates.clear();

  for(unsigned int k=0; k<vec_kaon_bars.size(); k++){
    X_weights += ADC_Low_TARGET[vec_kaon_bars[k]]*Xloc[vec_kaon_bars[k]];
    Y_weights += ADC_Low_TARGET[vec_kaon_bars[k]]*Yloc[vec_kaon_bars[k]];
    total_energy += ADC_Low_TARGET[vec_kaon_bars[k]];
  }

  if(vec_kaon_bars.size()>0 || total_energy!=0){
    vec_kaon_centroid_coordinates.push_back(X_weights/total_energy);
    vec_kaon_centroid_coordinates.push_back(Y_weights/total_energy);
  }
  else{
    vec_kaon_centroid_coordinates.push_back(999.99);
    vec_kaon_centroid_coordinates.push_back(999.99);
  }

  if(vec_kaon_bars.size()>0){
    if(gap_to_fit==6 || gap_to_fit==12){
      vec_fit_lines_intersect = _2lines_intersect(a_fit_lepton_rotate, b_fit_lepton_rotate, a_fit_kaon, b_fit_kaon);
    }
    else{
      vec_fit_lines_intersect = _2lines_intersect(a_lepton_fit_3, b_lepton_fit_3, a_fit_kaon, b_fit_kaon);
    }
  }
  else{
    vec_fit_lines_intersect.push_back(999.99);
    vec_fit_lines_intersect.push_back(999.99);
  }

  if(vec_kaon_bars.size() < 5){
    vec_k_stop_coordinates.push_back(vec_kaon_centroid_coordinates[0]);
    vec_k_stop_coordinates.push_back(vec_kaon_centroid_coordinates[1]);
  }
  else{
    vec_k_stop_coordinates.push_back(vec_fit_lines_intersect[0]);
    vec_k_stop_coordinates.push_back(vec_fit_lines_intersect[1]);
  }

  vec_xx_kaon_stop.push_back(vec_k_stop_coordinates[0]);
  vec_yy_kaon_stop.push_back(vec_k_stop_coordinates[1]);






  int i_kaon_bar = 999;
  i_kaon_bar = kaon_fiber(vec_k_stop_coordinates[0],vec_k_stop_coordinates[1]);
  double length_in_target = 0;

  length_in_target = distance(vec_xx_kaon_stop[0],vec_yy_kaon_stop[0],
                              vec_xx_int_TDC_TARGET[0],vec_yy_int_TDC_TARGET[0]); //CORRECT



  ///////////////////////////////////////////////////////////////////////////////////////////////////////////





    //
    // delete fChain;
    // delete ell; delete ell_Target; delete ell_L1;
    // delete gr_kaon_fit;
    // delete gr_kaon;
    // delete gr_kaon_bk;
    // delete gr_kaon_fit;
    // delete gr_kaon_fit;
    // if(gr_lepton_1 != NULL)    delete gr_lepton_1;
    // return;
    // delete func_lepton_fit_1;
    // delete gr_lepton_1;
    // delete func_lepton_fit_2;
    // delete gr_lepton_2;
    // delete gr2_Leptons_rotate;
    // delete func_lepton_fit_3;
    // delete gr_lepton_3;
    // delete gr3_Leptons_rotate;
    // delete x_guide;
    // delete y_guide;
    // delete x_guide_rotate;
    // delete y_guide_rotate;
    // delete tex_Angle_guide;
    // delete tex_ChiS;
    // delete best_fit_rotate; delete A1;
    // delete A2;
    // delete x_sft;
    // delete y_sft; delete x_tof1;
    // delete y_tof1; delete x_target;
    // delete y_target;
    // delete k_stop_line;
    // delete tex_Kstop_X; delete tex_Kstop_Y;
    // delete tex_Kbar;
    // delete gr_TOF1_Markers;
    // delete gr_Target_Center; delete gr_TOF1;
    // delete gr_TOF1_closest; delete gr_int_TDC_Gap_Fibers;
    // delete gr_int_TDC_Gap_Fibers_SFT; delete gr_int_TDC_TARGET;
    // delete gr_kaon_stop;
    /*
    if (graphics == 0) {
      c2->Close();
      // delete h_ADC_L1_DS;
      // delete h_ADC_L2_DS;
      // delete h_ADC_L3_DS;
      // delete h_ADC_L4_DS;
      // delete h_ADC_L1_US;
      // delete h_ADC_L2_US;
      // delete h_ADC_L3_US;
      // delete h_ADC_L4_US;






      // delete TOF_line1; delete TOF_line2; delete TOF_line3; delete TOF_line4;
      // delete TOF_line5; delete TOF_line6; delete TOF_line7; delete TOF_line8;
      // delete TOF_line9; delete TOF_line10; delete TOF_line11; delete TOF_line12;
      //
      // delete TOF_line13; delete TOF_line14; delete TOF_line15; delete TOF_line16;
      // delete TOF_line17; delete TOF_line18; delete TOF_line19; delete TOF_line20;
      // delete TOF_line21; delete TOF_line22; delete TOF_line23; delete TOF_line24;

    }
    */
    // delete c2;




  // return;



  if(display != 0){
    cout << endl;
    cout << endl;
    cout << " //////////////////   Chi Squared Calculation   ///////////////////////" << endl;
    cout << endl;
    cout << fixed;

    if((gap_to_fit_rotate==6 || gap_to_fit_rotate==12 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
      sumS_rotate = 0.;
      for(unsigned j=0; j<vec_xx_lepton_rotate.size(); j++){
        vec_yprime_rotate.push_back(-a_lepton_fit_3*vec_yy_lepton_rotate[j] - b_lepton_fit_3);
        vec_Dy_rotate.push_back(vec_xx_lepton_rotate[j] - (-a_lepton_fit_3*vec_yy_lepton_rotate[j] - b_lepton_fit_3));
      }

      cout << "y = mx + b  with:  m = ";
      cout << setw(7) << setprecision(3) << a_lepton_fit_3;
      cout << " and b = " << setw(8) << setprecision(3) << b_lepton_fit_3 << endl;
      cout << endl;
      cout << " Bar#      X         Y         X'        DX          DX^2" << endl;
      for(unsigned i=0; i<vec_xx_lepton_rotate.size(); i++){
        cout << " " << setw(3) << vec_bar_rotate[i] << "   ";
        cout << " " << setw(6) << setprecision(2) << vec_xx_lepton_rotate[i] << "   ";
        cout << " " << setw(6) << setprecision(2) << vec_yy_lepton_rotate[i] << "   ";
        cout << " " << setw(7) << setprecision(3) << vec_yprime_rotate[i] << "   ";
        cout << " " << setw(6) << setprecision(3) << vec_Dy_rotate[i] << "   ";
        cout << " " << setw(8) << setprecision(3) << pow(vec_Dy_rotate[i],2) << "   ";
        cout << endl;
        sumS_rotate += pow(vec_Dy_rotate[i],2);
      }
      cout << endl;
      cout << "Sum of DX^2 = " << sumS_rotate << endl;
    }
    else{
      sumS = 0.;
      for(unsigned j=0; j<vec_xx_lepton.size(); j++){
        vec_yprime.push_back(a_lepton_fit_3*vec_xx_lepton[j] + b_lepton_fit_3);
        vec_Dy.push_back(vec_yy_lepton[j] - (a_lepton_fit_3*vec_xx_lepton[j] + b_lepton_fit_3));
      }

      cout << "y = mx + b  with:  m = ";
      cout << setw(7) << setprecision(3) << a_lepton_fit_3;
      cout << " and b = " << setw(8) << setprecision(3) << b_lepton_fit_3 << endl;
      cout << endl;
      cout << " Bar#      X         Y         Y'        DY          DY^2" << endl;
      for(unsigned i=0; i<vec_xx_lepton.size()-3; i++){
        cout << " " << setw(3) << vec_bar[i] << "   ";
        cout << " " << setw(6) << setprecision(2) << vec_xx_lepton[i] << "   ";
        cout << " " << setw(6) << setprecision(2) << vec_yy_lepton[i] << "   ";
        cout << " " << setw(7) << setprecision(3) << vec_yprime[i] << "   ";
        cout << " " << setw(6) << setprecision(3) << vec_Dy[i] << "   ";
        cout << " " << setw(8) << setprecision(3) << pow(vec_Dy[i],2) << "   ";
        cout << endl;
        sumS += pow(vec_Dy[i],2);
      }
      cout << endl;
      cout << "Sum of DY^2 = " << sumS << endl;
      cout << endl;
    }
  }

  // Calculate TDC_diff
  tdc_ck_corr = tdc_trigger[0][0] - TDC_ck_avg2;
  TDC_diff = 550 - (0.91*TDC_average - 0.625*tdc_ck_corr);


  ofstream fout;
  fout.open("output_Batch_5.9.txt", ios::app);

  fout << fixed;
  fout << setw(4) << Run_Number << "  ";
  fout << setw(7) << ievt << "  ";
  fout << setw(2) << gap_to_fit_rotate << "  ";
  fout << setw(2) << selected_TOF2 << "  ";
  fout << setw(7) << setprecision(3) << angle_final_guide << "  ";
  fout << setw(5) << setprecision(2) << Delta_phi_deg << "  ";
  fout << setw(8) << setprecision(2) << Chis_lepton_fit_3 << "  ";
  fout << setw(3) << setprecision(0) << ndf_lepton_fit_3 << "  ";
  fout << setw(7) << setprecision(2) << Chis_lepton_fit_3/ndf_lepton_fit_3 << "  ";
  fout << setw(3) << lepton_counter << "  ";
  fout << setw(7) << setprecision(2) << vec_xx_int_TDC_TARGET[0] << "  ";
  fout << setw(7) << setprecision(2) << vec_yy_int_TDC_TARGET[0] << "  ";
  fout << setw(6) << setprecision(2) << vec_kaon_centroid_coordinates[0] << "  ";
  fout << setw(6) << setprecision(2) << vec_kaon_centroid_coordinates[1] << "  ";
  fout << setw(6) << setprecision(2) << vec_fit_lines_intersect[0] << "  ";
  fout << setw(6) << setprecision(2) << vec_fit_lines_intersect[1] << "  ";
  fout << setw(6) << setprecision(2) << vec_k_stop_coordinates[0] << "  ";
  fout << setw(6) << setprecision(2) << vec_k_stop_coordinates[1] << "  ";
  fout << setw(3) << i_kaon_bar << "  ";
  fout << setw(3) << vec_kaon_bars.size() << "  ";
  fout << setw(6) << setprecision(2) << Chis_kaon/ndf_kaon << "  ";
  fout << setw(3) << vec_Ck.size() << "  ";
  fout << setw(3) << vec_Cpi.size() << "  ";
  fout << setw(8) << setprecision(3) << length_in_target << "  ";
  fout << setw(8) << setprecision(3) << C2X_centroid << "  ";
  fout << setw(6) << setprecision(1) << TDC_diff << "  ";
  fout << endl;
  fout.close();





  return;

} // End void
