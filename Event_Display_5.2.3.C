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

double SFT_Test_CR_Event(int Run_Number, int evt, double phi, bool to_print, double C2X_centroid, double length_in_target);
vector<double> SFT_Test_CR_Event2(int Run_Number, int evt, double phi, double C2X_centroid, double length_in_target);
int kaon_fiber(double x_bar, double y_bar);

void Event_Display_5_3(Int_t Run_Number=5, Int_t ievt=0, int walk=0, Int_t enable_cout=0, Int_t display = 0, Int_t graphics = 0){

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

  gROOT->Clear();

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

  TChain *fChain= new TChain("Tree");
  fChain->Add(Name_finput);
  fChain->SetMakeClass(1);

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

  Int_t tdc_Ck[14][16];
  Int_t tdc_Cpi[14][16];

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

  vector <int> lepton_bars;
  vector <int> lepton_bars_rotate;
  vector <int> kaon_bars;

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
  fChain->SetBranchAddress("VT48_TDC",tdc_vt48);
  vector<int> vec_tdc_b0_6;  vec_tdc_b0_6.clear();
  vector<int> vec_tdc_b0_7;  vec_tdc_b0_7.clear();
  //int tdc_b0_6 = -1; int tdc_b0_7 = -1; int tdc_b0 = -1;

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

  cout << "   " << endl;
  cout << "   " << endl;
  cout << "************************************************************************************************************" << endl;

  cout << "File opened:  " << Name_finput << endl;
  if(!NoFile) cout << "Time Walk Parameter File:  " << ParsTarg1 << endl;
  if(NoFile) cout << "Time Walk Parameter File:  " << "NOT FOUND !" << endl;
  cout << "SFT Mapping File:  " << par_finput << endl;
  cout << "MWPC Mapping File:  " << par_finput2 << endl;

  sprintf(run_string,"Run %d ; Event %d",Run_Number,ievt);
  sprintf(event_string,"Run %d ; Event %d",Run_Number,ievt);

  fChain->SetBranchAddress("TDC_Trig",tdc_trigger); // tdc_trigger for TDC_diff calculation

  fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);    fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
  fChain->SetBranchAddress("ADC_Low_TARGET",adc_low_target);      fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
  fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);        fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
  fChain->SetBranchAddress("TDC_TE_TARGET",tdc_te_target);        fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

  fChain->SetBranchAddress("ADC_TOF1U",ADC_tof1U);
  fChain->SetBranchAddress("ADC_TOF1D",ADC_tof1D);
  fChain->SetBranchAddress("TDC_TOF1U",TDC_tof1U);
  fChain->SetBranchAddress("TDC_TOF1D",TDC_tof1D);

  fChain->SetBranchAddress("ADC_TOF2AO",ADC_tof2AO);
  fChain->SetBranchAddress("ADC_TOF2AI",ADC_tof2AI);
  fChain->SetBranchAddress("ADC_TOF2BO",ADC_tof2BO);
  fChain->SetBranchAddress("ADC_TOF2BI",ADC_tof2BI);
  fChain->SetBranchAddress("TDC_TOF2AO",TDC_tof2AO);
  fChain->SetBranchAddress("TDC_TOF2AI",TDC_tof2AI);
  fChain->SetBranchAddress("TDC_TOF2BO",TDC_tof2BO);
  fChain->SetBranchAddress("TDC_TOF2BI",TDC_tof2BI);

  fChain->SetBranchAddress("MWPCADC",MwpcADC);

  fChain->SetBranchAddress("ADC_C2X_R",adc_c2x_r);
  fChain->SetBranchAddress("ADC_C2X_L",adc_c2x_l);
  fChain->SetBranchAddress("ADC_C2Y_R",adc_c2y_r);
  fChain->SetBranchAddress("ADC_C2Y_L",adc_c2y_l);
  fChain->SetBranchAddress("ADC_C3X_R",adc_c3x_r);
  fChain->SetBranchAddress("ADC_C3X_L",adc_c3x_l);
  fChain->SetBranchAddress("ADC_C3Y_R",adc_c3y_r);
  fChain->SetBranchAddress("ADC_C3Y_L",adc_c3y_l);
  fChain->SetBranchAddress("ADC_C4X_R",adc_c4x_r);
  fChain->SetBranchAddress("ADC_C4X_L",adc_c4x_l);
  fChain->SetBranchAddress("ADC_C4Y_R",adc_c4y_r);
  fChain->SetBranchAddress("ADC_C4Y_L",adc_c4y_l);

  fChain->SetBranchAddress("TDC_Ck", tdc_Ck);
  fChain->SetBranchAddress("TDC_Cpi", tdc_Cpi);

  fChain->SetBranchAddress("EvFlag", Event_flag);

  Int_t nentries = (Int_t)fChain->GetEntries();
  cout << "Total Number of Events:  " << nentries <<endl;
  //cout << "  " << endl;
  cout << "Event_Display.C -- " << Version << endl;
  cout << "************************************************************************************************************" << endl;
  cout << "  " << endl;
  cout << "  " << endl;

  bool Good_Event=false;

  for(int ivt=ievt; ivt<ievt+1; ivt++){
    fChain->GetEntry(ivt);

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

    if(Good_Event) cout << "Event: "<< ievt << "   --  GOOD EVENT!" << endl;

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

  cout << "  " << endl;
  cout << "  " << endl;

  if(!Good_Event) return;

  cout << " **** RAW ADCs" << endl;
  cout << " ///////////   TOF1   ///////////   //////////////////////   TOF2   //////////////////////" << endl;
  for(Int_t j_TOF=0; j_TOF<12; j_TOF++){
    printf(" %2d   UP: %4d  |  DOWN: %4d    ||    %2d   AO: %4d  |  AI: %4d  |  BO: %4d  |  BI:%4d\n",
    j_TOF+1, ADC_tof1U[j_TOF], ADC_tof1D[j_TOF], j_TOF+1, ADC_tof2AO[j_TOF], ADC_tof2AI[j_TOF], ADC_tof2BO[j_TOF], ADC_tof2BI[j_TOF]);
  }

  cout << "  " << endl;
  cout << "  " << endl;

  for(int j=0 ; j<128 ; j++){
    if(ADC_High_SFT[j]<0)      ADC_High_SFT_corr[j]=0;
    if(ADC_High_SFT[j]>=0)     ADC_High_SFT_corr[j]=ADC_High_SFT[j];
  }


  ////////////////////  GEOMETRY !!!  ////////////////////////////////////////////

  // ######### TARGET !
  TLine *hline1 = new TLine(0.38,0.14,0.62,0.14);    hline1->SetLineWidth(2);
  TLine *hline2 = new TLine(0.30,0.18,0.70,0.18);    hline2->SetLineWidth(2);
  TLine *hline3 = new TLine(0.26,0.22,0.74,0.22);    hline3->SetLineWidth(2);
  TLine *hline4 = new TLine(0.22,0.26,0.78,0.26);    hline4->SetLineWidth(2);
  TLine *hline5 = new TLine(0.18,0.30,0.82,0.30);    hline5->SetLineWidth(2);
  TLine *hline6 = new TLine(0.18,0.34,0.82,0.34);    hline6->SetLineWidth(2);
  TLine *hline7 = new TLine(0.14,0.38,0.86,0.38);    hline7->SetLineWidth(2);
  TLine *hline8 = new TLine(0.14,0.42,0.86,0.42);    hline8->SetLineWidth(2);
  TLine *hline9 = new TLine(0.14,0.46,0.86,0.46);    hline9->SetLineWidth(2);
  TLine *hline10 = new TLine(0.14,0.50,0.86,0.50);   hline10->SetLineWidth(2);
  TLine *hline11 = new TLine(0.14,0.54,0.86,0.54);   hline11->SetLineWidth(2);
  TLine *hline12 = new TLine(0.14,0.58,0.86,0.58);   hline12->SetLineWidth(2);
  TLine *hline13 = new TLine(0.14,0.62,0.86,0.62);   hline13->SetLineWidth(2);
  TLine *hline14 = new TLine(0.18,0.66,0.82,0.66);   hline14->SetLineWidth(2);
  TLine *hline15 = new TLine(0.18,0.70,0.82,0.70);   hline15->SetLineWidth(2);
  TLine *hline16 = new TLine(0.22,0.74,0.78,0.74);   hline16->SetLineWidth(2);
  TLine *hline17 = new TLine(0.26,0.78,0.74,0.78);   hline17->SetLineWidth(2);
  TLine *hline18 = new TLine(0.30,0.82,0.70,0.82);   hline18->SetLineWidth(2);
  TLine *hline19 = new TLine(0.38,0.86,0.62,0.86);   hline19->SetLineWidth(2);

  TLine *vline1 = new TLine(0.14,0.38,0.14,0.62);    vline1->SetLineWidth(2);
  TLine *vline2 = new TLine(0.18,0.30,0.18,0.70);    vline2->SetLineWidth(2);
  TLine *vline3 = new TLine(0.22,0.26,0.22,0.74);    vline3->SetLineWidth(2);
  TLine *vline4 = new TLine(0.26,0.22,0.26,0.78);    vline4->SetLineWidth(2);
  TLine *vline5 = new TLine(0.30,0.18,0.30,0.82);    vline5->SetLineWidth(2);
  TLine *vline6 = new TLine(0.34,0.18,0.34,0.82);    vline6->SetLineWidth(2);
  TLine *vline7 = new TLine(0.38,0.14,0.38,0.86);    vline7->SetLineWidth(2);
  TLine *vline8 = new TLine(0.42,0.14,0.42,0.86);    vline8->SetLineWidth(2);
  TLine *vline9 = new TLine(0.46,0.14,0.46,0.86);    vline9->SetLineWidth(2);
  TLine *vline10 = new TLine(0.50,0.14,0.50,0.86);   vline10->SetLineWidth(2);
  TLine *vline11 = new TLine(0.54,0.14,0.54,0.86);   vline11->SetLineWidth(2);
  TLine *vline12 = new TLine(0.58,0.14,0.58,0.86);   vline12->SetLineWidth(2);
  TLine *vline13 = new TLine(0.62,0.14,0.62,0.86);   vline13->SetLineWidth(2);
  TLine *vline14 = new TLine(0.66,0.18,0.66,0.82);   vline14->SetLineWidth(2);
  TLine *vline15 = new TLine(0.70,0.18,0.70,0.82);   vline15->SetLineWidth(2);
  TLine *vline16 = new TLine(0.74,0.22,0.74,0.78);   vline16->SetLineWidth(2);
  TLine *vline17 = new TLine(0.78,0.26,0.78,0.74);   vline17->SetLineWidth(2);
  TLine *vline18 = new TLine(0.82,0.30,0.82,0.70);   vline18->SetLineWidth(2);
  TLine *vline19 = new TLine(0.86,0.38,0.86,0.62);   vline19->SetLineWidth(2);

  TLine *vblue1 = new TLine(0.38,0.70,0.38,0.74);      vblue1->SetLineWidth(5);   vblue1->SetLineColor(kBlue-9);
  TLine *vblue2 = new TLine(0.42,0.62,0.42,0.66);      vblue2->SetLineWidth(5);   vblue2->SetLineColor(kBlue-9);
  TLine *vblue3 = new TLine(0.30,0.54,0.30,0.58);      vblue3->SetLineWidth(5);   vblue3->SetLineColor(kBlue-9);
  TLine *vblue4 = new TLine(0.86,0.50,0.86,0.54);      vblue4->SetLineWidth(5);   vblue4->SetLineColor(kBlue-9);
  TLine *vblue5 = new TLine(0.70,0.42,0.70,0.46);      vblue5->SetLineWidth(5);   vblue5->SetLineColor(kBlue-9);
  TLine *vblue6 = new TLine(0.58,0.34,0.58,0.38);      vblue6->SetLineWidth(5);   vblue6->SetLineColor(kBlue-9);
  TLine *vblue7 = new TLine(0.62,0.26,0.62,0.30);      vblue7->SetLineWidth(5);   vblue7->SetLineColor(kBlue-9);
  TLine *vblue8 = new TLine(0.62,0.14,0.62,0.18);      vblue8->SetLineWidth(5);   vblue8->SetLineColor(kBlue-9);

  ////////////////TARGET TOF
  TLine *TOF_line9 = new TLine(0.052539,0.38,0.052539,0.62);
  TLine *TOF_line8 = new TLine(0.172539,0.172154,0.052539,0.38);
  TLine *TOF_line7 = new TLine(0.172539,0.172154,0.38,0.0521539);
  TLine *TOF_line6 = new TLine(0.38,0.0521539,0.62,0.0521539);
  TLine *TOF_line5 = new TLine(0.62,0.0521539,0.8278461,0.172154);
  TLine *TOF_line4 = new TLine(0.8278461,0.172154,0.9478461,0.38);
  TLine *TOF_line3 = new TLine(0.9478461,0.38,0.9478461,0.62);
  TLine *TOF_line2 = new TLine(0.9478461,0.62,0.8278461,0.8278461);
  TLine *TOF_line1 = new TLine(0.8278461,0.8278461,0.62,0.9478461);
  TLine *TOF_line12 = new TLine(0.62,0.9478461,0.38,0.9478461);
  TLine *TOF_line11 = new TLine(0.38,0.9478461,0.172539,0.8278461);
  TLine *TOF_line10 = new TLine(0.052539,0.62,0.172539,0.8278461);

  TLine *TOF1_line[12];
  TOF1_line[0] = TOF_line1;
  TOF1_line[1] = TOF_line2;
  TOF1_line[2] = TOF_line3;
  TOF1_line[3] = TOF_line4;
  TOF1_line[4] = TOF_line5;
  TOF1_line[5] = TOF_line6;
  TOF1_line[6] = TOF_line7;
  TOF1_line[7] = TOF_line8;
  TOF1_line[8] = TOF_line9;
  TOF1_line[9] = TOF_line10;
  TOF1_line[10] = TOF_line11;
  TOF1_line[11] = TOF_line12;

  for(int i=0; i<12; i++){
  	TOF1_line[i]->SetLineWidth(10);
 	  TOF1_line[i]->SetLineColor(17);
  }

  TOF_line1->SetLineColor(17);
  TOF_line2->SetLineColor(17);
  TOF_line3->SetLineColor(17);
  TOF_line4->SetLineColor(17);
  TOF_line5->SetLineColor(17);
  TOF_line6->SetLineColor(17);
  TOF_line7->SetLineColor(17);
  TOF_line8->SetLineColor(17);
  TOF_line9->SetLineColor(17);
  TOF_line10->SetLineColor(17);
  TOF_line11->SetLineColor(17);
  TOF_line12->SetLineColor(17);

  TLine *TOF_line13 = new TLine(0.8348214,0.8393338,0.625558,0.9620716);        TOF_line13->SetLineWidth(20);
  TLine *TOF_line14 = new TLine(0.9603795,0.627551,0.8394717,0.8369272);        TOF_line14->SetLineWidth(20);
  TLine *TOF_line15 = new TLine(0.961,0.38,0.961,0.62);                         TOF_line15->SetLineWidth(20);
  TLine *TOF_line16 = new TLine(0.8394717,0.1630728,0.9580543,0.372449);        TOF_line16->SetLineWidth(20);
  TLine *TOF_line17 = new TLine(0.6232329,0.040335,0.8324963,0.1606662);        TOF_line17->SetLineWidth(20);
  TLine *TOF_line18 = new TLine(0.38,0.039,0.62,0.039);                         TOF_line18->SetLineWidth(20);
  TLine *TOF_line19 = new TLine(0.1651786,0.1606662,0.3721168,0.040335);        TOF_line19->SetLineWidth(20);
  TLine *TOF_line20 = new TLine(0.1605283,0.1630728,0.03962054,0.372449);       TOF_line20->SetLineWidth(20);
  TLine *TOF_line21 = new TLine(0.04,0.38,0.04,0.62);                           TOF_line21->SetLineWidth(20);
  TLine *TOF_line22 = new TLine(0.03962054,0.6251444,0.1605283,0.8345206);      TOF_line22->SetLineWidth(20);
  TLine *TOF_line23 = new TLine(0.3721168,0.9572584,0.1605283,0.8369272);       TOF_line23->SetLineWidth(20);
  TLine *TOF_line24 = new TLine(0.62,0.96,0.38,0.96);                           TOF_line24->SetLineWidth(20);

  TOF_line13->SetLineColor(15);
  TOF_line14->SetLineColor(15);
  TOF_line15->SetLineColor(15);
  TOF_line16->SetLineColor(15);
  TOF_line17->SetLineColor(15);
  TOF_line18->SetLineColor(15);
  TOF_line19->SetLineColor(15);
  TOF_line20->SetLineColor(15);
  TOF_line21->SetLineColor(15);
  TOF_line22->SetLineColor(15);
  TOF_line23->SetLineColor(15);
  TOF_line24->SetLineColor(15);

  // ######### SFT !

  TLine *hline_DS[8];   TLine *vline_DS[68];
  TLine *hline_US[8];   TLine *vline_US[68];
  TLine *split_line;

  // Horizontal Lines
  hline_DS[0] = new TLine(0.1,0.8,0.7,0.8);   hline_US[0] = new TLine(0.1,0.76,0.7,0.76);
  hline_DS[1] = new TLine(0.1,0.76,0.7,0.76);   hline_US[1] = new TLine(0.1,0.72,0.7,0.72);
  hline_DS[2] = new TLine(0.1,0.66,0.7,0.66);   hline_US[2] = new TLine(0.1,0.62,0.7,0.62);
  hline_DS[3] = new TLine(0.1,0.62,0.7,0.62);   hline_US[3] = new TLine(0.1,0.58,0.7,0.58);
  hline_DS[4] = new TLine(0.1,0.52,0.78,0.52);    hline_US[4] = new TLine(0.1,0.48,0.78,0.48);
  hline_DS[5] = new TLine(0.1,0.48,0.78,0.48);    hline_US[5] = new TLine(0.1,0.44,0.78,0.44);
  hline_DS[6] = new TLine(0.1,0.38,0.78,0.38);    hline_US[6] = new TLine(0.1,0.34,0.78,0.34);
  hline_DS[7] = new TLine(0.1,0.34,0.78,0.34);    hline_US[7] = new TLine(0.1,0.30,0.78,0.30);

  for(Int_t ih_DS=0; ih_DS<8; ih_DS++)  hline_DS[ih_DS]->SetLineWidth(2);
  for(Int_t ih_US=0; ih_US<8; ih_US++)  hline_US[ih_US]->SetLineWidth(2);

  //Vertical Lines
  for(Int_t iv_DS1=0; iv_DS1<16; iv_DS1++)  vline_DS[iv_DS1] = new TLine(0.1+0.04*iv_DS1,0.76,0.1+0.04*iv_DS1,0.80);
  for(Int_t iv_DS2=0; iv_DS2<16; iv_DS2++)  vline_DS[iv_DS2+16] = new TLine(0.1+0.04*iv_DS2,0.62,0.1+0.04*iv_DS2,0.66);
  for(Int_t iv_DS3=0; iv_DS3<18; iv_DS3++)  vline_DS[iv_DS3+32] = new TLine(0.1+0.04*iv_DS3,0.48,0.1+0.04*iv_DS3,0.52);
  for(Int_t iv_DS4=0; iv_DS4<18; iv_DS4++)  vline_DS[iv_DS4+50] = new TLine(0.1+0.04*iv_DS4,0.34,0.1+0.04*iv_DS4,0.38);
  for(Int_t iv_DS=0; iv_DS<68; iv_DS++)   vline_DS[iv_DS]->SetLineWidth(2);

  for(Int_t iv_US1=0; iv_US1<16; iv_US1++)  vline_US[iv_US1] = new TLine(0.1+0.04*iv_US1,0.72,0.1+0.04*iv_US1,0.76);
  for(Int_t iv_US2=0; iv_US2<16; iv_US2++)  vline_US[iv_US2+16] = new TLine(0.1+0.04*iv_US2,0.58,0.1+0.04*iv_US2,0.62);
  for(Int_t iv_US3=0; iv_US3<18; iv_US3++)  vline_US[iv_US3+32] = new TLine(0.1+0.04*iv_US3,0.44,0.1+0.04*iv_US3,0.48);
  for(Int_t iv_US4=0; iv_US4<18; iv_US4++)  vline_US[iv_US4+50] = new TLine(0.1+0.04*iv_US4,0.30,0.1+0.04*iv_US4,0.34);
  for(Int_t iv_US=0; iv_US<68; iv_US++)   vline_US[iv_US]->SetLineWidth(2);

  //////////////////////// LEGENDS /////////////////////////////////////////////////

  //####### TARGET!
  TLatex *tex_Title_ADC_High_TARGET;
  TLatex *tex_Subtitle_ADC_High_TARGET;

  TLatex *tex_run_Title_TARGET;
  TLatex *tex_run_TARGET;
  TLatex *tex_event_Title_TARGET;
  TLatex *tex_event_TARGET;

  TLatex *tex_palette_TARGET[10];
  TLatex *tex_palette_TARGET2[10];
  TLatex *tex_Legend_TARGET[36];

  TLatex *tex_palette_MWPC_C2[10];
  TLatex *tex_palette_MWPC_C3[10];
  TLatex *tex_palette_MWPC_C4[10];

  tex_run_Title_TARGET = new TLatex(0.034,0.15,"Run :");
  tex_run_Title_TARGET->SetTextSize(0.038);
  tex_run_Title_TARGET->SetLineWidth(2);

  tex_run_TARGET = new TLatex(0.034,0.125,run_string);
  tex_run_TARGET->SetTextSize(0.038);
  tex_run_TARGET->SetLineWidth(2);

  tex_event_Title_TARGET = new TLatex(0.034,0.07,"Event :");
  tex_event_Title_TARGET->SetTextSize(0.038);
  tex_event_Title_TARGET->SetLineWidth(2);

  tex_event_TARGET = new TLatex(0.036,0.0,event_string);
  tex_event_TARGET->SetTextSize(0.038);
  tex_event_TARGET->SetLineWidth(2);

  tex_Legend_TARGET[0] = new TLatex(0.36,0.83,"0");     tex_Legend_TARGET[9] = new TLatex(0.09,0.47,"128");
  tex_Legend_TARGET[1] = new TLatex(0.28,0.79,"6");     tex_Legend_TARGET[10] = new TLatex(0.09,0.43,"146");
  tex_Legend_TARGET[2] = new TLatex(0.225,0.75,"16");     tex_Legend_TARGET[11] = new TLatex(0.09,0.39,"164");
  tex_Legend_TARGET[3] = new TLatex(0.185,0.71,"28");     tex_Legend_TARGET[12] = new TLatex(0.13,0.35,"182");
  tex_Legend_TARGET[4] = new TLatex(0.145,0.67,"42");     tex_Legend_TARGET[13] = new TLatex(0.13,0.31,"198");
  tex_Legend_TARGET[5] = new TLatex(0.145,0.63,"58");     tex_Legend_TARGET[14] = new TLatex(0.17,0.27,"214");
  tex_Legend_TARGET[6] = new TLatex(0.1025,0.59,"74");    tex_Legend_TARGET[15] = new TLatex(0.21,0.23,"228");
  tex_Legend_TARGET[7] = new TLatex(0.1025,0.55,"92");    tex_Legend_TARGET[16] = new TLatex(0.25,0.19,"240");
  tex_Legend_TARGET[8] = new TLatex(0.09,0.51,"110");     tex_Legend_TARGET[17] = new TLatex(0.33,0.15,"250");

  tex_Legend_TARGET[18] = new TLatex(0.635,0.83,"5");      tex_Legend_TARGET[27] = new TLatex(0.87,0.47,"145");
  tex_Legend_TARGET[19] = new TLatex(0.71,0.79,"15");     tex_Legend_TARGET[28] = new TLatex(0.87,0.43,"163");
  tex_Legend_TARGET[20] = new TLatex(0.75,0.75,"27");      tex_Legend_TARGET[29] = new TLatex(0.87,0.39,"181");
  tex_Legend_TARGET[21] = new TLatex(0.79,0.71,"41");      tex_Legend_TARGET[30] = new TLatex(0.83,0.35,"197");
  tex_Legend_TARGET[22] = new TLatex(0.83,0.67,"57");      tex_Legend_TARGET[31] = new TLatex(0.83,0.31,"213");
  tex_Legend_TARGET[23] = new TLatex(0.83,0.63,"73");      tex_Legend_TARGET[32] = new TLatex(0.79,0.27,"227");
  tex_Legend_TARGET[24] = new TLatex(0.885,0.59,"91");      tex_Legend_TARGET[33] = new TLatex(0.75,0.23,"239");
  tex_Legend_TARGET[25] = new TLatex(0.87,0.55,"109");     tex_Legend_TARGET[34] = new TLatex(0.71,0.19,"249");
  tex_Legend_TARGET[26] = new TLatex(0.87,0.51,"127");      tex_Legend_TARGET[35] = new TLatex(0.64,0.15,"255");

  for(Int_t jj=0; jj<36; jj++){
    tex_Legend_TARGET[jj]->SetTextSize(0.03);
    tex_Legend_TARGET[jj]->SetLineWidth(2);
  }

  //####### SFT!
  TLatex *tex_event_SFT;
  TLatex *tex_run_Title_SFT;
  TLatex *tex_run_SFT;
  TLatex *tex_event_Title_SFT;

  tex_run_Title_SFT = new TLatex(0.034,0.15,"Run :");
  tex_run_Title_SFT->SetTextSize(0.038);
  tex_run_Title_SFT->SetLineWidth(2);

  tex_run_SFT = new TLatex(0.034,0.125,run_string);
  tex_run_SFT->SetTextSize(0.020);
  tex_run_SFT->SetLineWidth(2);

  tex_event_Title_SFT = new TLatex(0.034,0.07,"Event :");
  tex_event_Title_SFT->SetTextSize(0.038);
  tex_event_Title_SFT->SetLineWidth(2);

  tex_event_SFT = new TLatex(0.036,0.045,event_string);
  tex_event_SFT->SetTextSize(0.02);
  tex_event_SFT->SetLineWidth(2);


  //////////////////////////////////   MARKERS   //////////////////////////////////////////

  // ######### TARGET !
  TMarker *marker_ADC_TARGET[256];
  TMarker *marker_ADCL_TARGET[256];
  TMarker *marker_TDC_TARGET[256];
  TMarker *palette_TARGET[10];
  TMarker *palette_MWPC_C2[10];
  TMarker *palette_MWPC_C3[10];
  TMarker *palette_MWPC_C4[10];

  for(Int_t i1=0; i1<6; i1++) marker_ADC_TARGET[i1] = new TMarker(0.40+0.04*i1,0.84,21);
  for(Int_t i2=0; i2<10; i2++) marker_ADC_TARGET[i2+6] = new TMarker(0.32+0.04*i2,0.80,21);
  for(Int_t i3=0; i3<12; i3++) marker_ADC_TARGET[i3+16] = new TMarker(0.28+0.04*i3,0.76,21);
  for(Int_t i4=0; i4<14; i4++) marker_ADC_TARGET[i4+28] = new TMarker(0.24+0.04*i4,0.72,21);
  for(Int_t i5=0; i5<16; i5++) marker_ADC_TARGET[i5+42] = new TMarker(0.20+0.04*i5,0.68,21);
  for(Int_t i6=0; i6<16; i6++) marker_ADC_TARGET[i6+58] = new TMarker(0.20+0.04*i6,0.64,21);
  for(Int_t i7=0; i7<18; i7++) marker_ADC_TARGET[i7+74] = new TMarker(0.16+0.04*i7,0.60,21);
  for(Int_t i8=0; i8<18; i8++) marker_ADC_TARGET[i8+92] = new TMarker(0.16+0.04*i8,0.56,21);
  for(Int_t i9=0; i9<18; i9++) marker_ADC_TARGET[i9+110] = new TMarker(0.16+0.04*i9,0.52,21);
  for(Int_t i10=0; i10<18; i10++) marker_ADC_TARGET[i10+128] = new TMarker(0.16+0.04*i10,0.48,21);
  for(Int_t i11=0; i11<18; i11++) marker_ADC_TARGET[i11+146] = new TMarker(0.16+0.04*i11,0.44,21);
  for(Int_t i12=0; i12<18; i12++) marker_ADC_TARGET[i12+164] = new TMarker(0.16+0.04*i12,0.40,21);
  for(Int_t i13=0; i13<16; i13++) marker_ADC_TARGET[i13+182] = new TMarker(0.20+0.04*i13,0.36,21);
  for(Int_t i14=0; i14<16; i14++) marker_ADC_TARGET[i14+198] = new TMarker(0.20+0.04*i14,0.32,21);
  for(Int_t i15=0; i15<14; i15++) marker_ADC_TARGET[i15+214] = new TMarker(0.24+0.04*i15,0.28,21);
  for(Int_t i16=0; i16<12; i16++) marker_ADC_TARGET[i16+228] = new TMarker(0.28+0.04*i16,0.24,21);
  for(Int_t i17=0; i17<10; i17++) marker_ADC_TARGET[i17+240] = new TMarker(0.32+0.04*i17,0.20,21);
  for(Int_t i18=0; i18<6; i18++) marker_ADC_TARGET[i18+250] = new TMarker(0.40+0.04*i18,0.16,21);
  for(Int_t iSize=0; iSize<256; iSize++)  marker_ADC_TARGET[iSize]->SetMarkerSize(flag_size_TARGET);

  for(Int_t i1=0; i1<6; i1++) marker_ADCL_TARGET[i1] = new TMarker(0.40+0.04*i1,0.84,21);
  for(Int_t i2=0; i2<10; i2++) marker_ADCL_TARGET[i2+6] = new TMarker(0.32+0.04*i2,0.80,21);
  for(Int_t i3=0; i3<12; i3++) marker_ADCL_TARGET[i3+16] = new TMarker(0.28+0.04*i3,0.76,21);
  for(Int_t i4=0; i4<14; i4++) marker_ADCL_TARGET[i4+28] = new TMarker(0.24+0.04*i4,0.72,21);
  for(Int_t i5=0; i5<16; i5++) marker_ADCL_TARGET[i5+42] = new TMarker(0.20+0.04*i5,0.68,21);
  for(Int_t i6=0; i6<16; i6++) marker_ADCL_TARGET[i6+58] = new TMarker(0.20+0.04*i6,0.64,21);
  for(Int_t i7=0; i7<18; i7++) marker_ADCL_TARGET[i7+74] = new TMarker(0.16+0.04*i7,0.60,21);
  for(Int_t i8=0; i8<18; i8++) marker_ADCL_TARGET[i8+92] = new TMarker(0.16+0.04*i8,0.56,21);
  for(Int_t i9=0; i9<18; i9++) marker_ADCL_TARGET[i9+110] = new TMarker(0.16+0.04*i9,0.52,21);
  for(Int_t i10=0; i10<18; i10++) marker_ADCL_TARGET[i10+128] = new TMarker(0.16+0.04*i10,0.48,21);
  for(Int_t i11=0; i11<18; i11++) marker_ADCL_TARGET[i11+146] = new TMarker(0.16+0.04*i11,0.44,21);
  for(Int_t i12=0; i12<18; i12++) marker_ADCL_TARGET[i12+164] = new TMarker(0.16+0.04*i12,0.40,21);
  for(Int_t i13=0; i13<16; i13++) marker_ADCL_TARGET[i13+182] = new TMarker(0.20+0.04*i13,0.36,21);
  for(Int_t i14=0; i14<16; i14++) marker_ADCL_TARGET[i14+198] = new TMarker(0.20+0.04*i14,0.32,21);
  for(Int_t i15=0; i15<14; i15++) marker_ADCL_TARGET[i15+214] = new TMarker(0.24+0.04*i15,0.28,21);
  for(Int_t i16=0; i16<12; i16++) marker_ADCL_TARGET[i16+228] = new TMarker(0.28+0.04*i16,0.24,21);
  for(Int_t i17=0; i17<10; i17++) marker_ADCL_TARGET[i17+240] = new TMarker(0.32+0.04*i17,0.20,21);
  for(Int_t i18=0; i18<6; i18++) marker_ADCL_TARGET[i18+250] = new TMarker(0.40+0.04*i18,0.16,21);
  for(Int_t iSize=0; iSize<256; iSize++)  marker_ADCL_TARGET[iSize]->SetMarkerSize(flag_size_TARGET);

  for(Int_t i1=0; i1<6; i1++) marker_TDC_TARGET[i1] = new TMarker(0.40+0.04*i1,0.84,21);
  for(Int_t i2=0; i2<10; i2++) marker_TDC_TARGET[i2+6] = new TMarker(0.32+0.04*i2,0.80,21);
  for(Int_t i3=0; i3<12; i3++) marker_TDC_TARGET[i3+16] = new TMarker(0.28+0.04*i3,0.76,21);
  for(Int_t i4=0; i4<14; i4++) marker_TDC_TARGET[i4+28] = new TMarker(0.24+0.04*i4,0.72,21);
  for(Int_t i5=0; i5<16; i5++) marker_TDC_TARGET[i5+42] = new TMarker(0.20+0.04*i5,0.68,21);
  for(Int_t i6=0; i6<16; i6++) marker_TDC_TARGET[i6+58] = new TMarker(0.20+0.04*i6,0.64,21);
  for(Int_t i7=0; i7<18; i7++) marker_TDC_TARGET[i7+74] = new TMarker(0.16+0.04*i7,0.60,21);
  for(Int_t i8=0; i8<18; i8++) marker_TDC_TARGET[i8+92] = new TMarker(0.16+0.04*i8,0.56,21);
  for(Int_t i9=0; i9<18; i9++) marker_TDC_TARGET[i9+110] = new TMarker(0.16+0.04*i9,0.52,21);
  for(Int_t i10=0; i10<18; i10++) marker_TDC_TARGET[i10+128] = new TMarker(0.16+0.04*i10,0.48,21);
  for(Int_t i11=0; i11<18; i11++) marker_TDC_TARGET[i11+146] = new TMarker(0.16+0.04*i11,0.44,21);
  for(Int_t i12=0; i12<18; i12++) marker_TDC_TARGET[i12+164] = new TMarker(0.16+0.04*i12,0.40,21);
  for(Int_t i13=0; i13<16; i13++) marker_TDC_TARGET[i13+182] = new TMarker(0.20+0.04*i13,0.36,21);
  for(Int_t i14=0; i14<16; i14++) marker_TDC_TARGET[i14+198] = new TMarker(0.20+0.04*i14,0.32,21);
  for(Int_t i15=0; i15<14; i15++) marker_TDC_TARGET[i15+214] = new TMarker(0.24+0.04*i15,0.28,21);
  for(Int_t i16=0; i16<12; i16++) marker_TDC_TARGET[i16+228] = new TMarker(0.28+0.04*i16,0.24,21);
  for(Int_t i17=0; i17<10; i17++) marker_TDC_TARGET[i17+240] = new TMarker(0.32+0.04*i17,0.20,21);
  for(Int_t i18=0; i18<6; i18++) marker_TDC_TARGET[i18+250] = new TMarker(0.40+0.04*i18,0.16,21);
  for(Int_t iSize=0; iSize<256; iSize++)  marker_TDC_TARGET[iSize]->SetMarkerSize(flag_size_TARGET);


  // ######### SFT !
  TMarker *marker_DS[64];
  TMarker *marker_US[64];

  for(Int_t i1=0; i1<15; i1++) marker_DS[i1] = new TMarker(0.12+0.04*i1,0.78,21);
  for(Int_t i2=0; i2<15; i2++) marker_DS[i2+15] = new TMarker(0.12+0.04*i2,0.64,21);
  for(Int_t i3=0; i3<17; i3++) marker_DS[i3+30] = new TMarker(0.12+0.04*i3,0.50,21);
  for(Int_t i4=0; i4<17; i4++) marker_DS[i4+47] = new TMarker(0.12+0.04*i4,0.36,21);

  for(Int_t iSize=0; iSize<64; iSize++) marker_DS[iSize]->SetMarkerSize(flag_size_SFT);

  for(Int_t iUS1=0; iUS1<15; iUS1++) marker_US[iUS1] = new TMarker(0.12+0.04*iUS1,0.74,21);
  for(Int_t iUS2=0; iUS2<15; iUS2++) marker_US[iUS2+15] = new TMarker(0.12+0.04*iUS2,0.60,21);
  for(Int_t iUS3=0; iUS3<17; iUS3++) marker_US[iUS3+30] = new TMarker(0.12+0.04*iUS3,0.46,21);
  for(Int_t iUS4=0; iUS4<17; iUS4++) marker_US[iUS4+47] = new TMarker(0.12+0.04*iUS4,0.32,21);
  for(Int_t iSize_US=0; iSize_US<64; iSize_US++)  marker_US[iSize_US]->SetMarkerSize(flag_size_SFT);


  //////////////////////////   Palette  ///////////////////////////
  Int_t ADC_palette_TARGET[10]={0,3,6,9,12,15,18,21,24,27};
  char ADC_palette_string_TARGET[10][100];

  Int_t ADC_palette_TARGET2[10]={0,1,2,3,4,5,6,7,8,9};
  char ADC_palette_string_TARGET2[10][100];

  for(int j=0; j<10; j++) sprintf(ADC_palette_string_TARGET[j],"%d",ADC_palette_TARGET[j]);
  sprintf(ADC_palette_string_TARGET[9],"%d+",ADC_palette_TARGET[9]);

  for(int j=0; j<10; j++) sprintf(ADC_palette_string_TARGET2[j],"%d",ADC_palette_TARGET2[j]);
  sprintf(ADC_palette_string_TARGET2[9],"%d+",ADC_palette_TARGET2[9]);

  palette_TARGET[0] = new TMarker(0.54,0.075,21);     palette_TARGET[0]->SetMarkerColor(kOrange+10);      palette_TARGET[0]->SetMarkerSize(flag_size_palette);
  palette_TARGET[1] = new TMarker(0.58,0.075,21);     palette_TARGET[1]->SetMarkerColor(kOrange+7);     palette_TARGET[1]->SetMarkerSize(flag_size_palette);
  palette_TARGET[2] = new TMarker(0.62,0.075,21);     palette_TARGET[2]->SetMarkerColor(kOrange+1);     palette_TARGET[2]->SetMarkerSize(flag_size_palette);
  palette_TARGET[3] = new TMarker(0.66,0.075,21);     palette_TARGET[3]->SetMarkerColor(kOrange-4);     palette_TARGET[3]->SetMarkerSize(flag_size_palette);
  palette_TARGET[4] = new TMarker(0.70,0.075,21);     palette_TARGET[4]->SetMarkerColor(kYellow-9);     palette_TARGET[4]->SetMarkerSize(flag_size_palette);
  palette_TARGET[5] = new TMarker(0.74,0.075,21);     palette_TARGET[5]->SetMarkerColor(kYellow-7);     palette_TARGET[5]->SetMarkerSize(flag_size_palette);
  palette_TARGET[6] = new TMarker(0.78,0.075,21);     palette_TARGET[6]->SetMarkerColor(kYellow-0);     palette_TARGET[6]->SetMarkerSize(flag_size_palette);
  palette_TARGET[7] = new TMarker(0.82,0.075,21);     palette_TARGET[7]->SetMarkerColor(kSpring-4);     palette_TARGET[7]->SetMarkerSize(flag_size_palette);
  palette_TARGET[8] = new TMarker(0.86,0.075,21);     palette_TARGET[8]->SetMarkerColor(kSpring-2);     palette_TARGET[8]->SetMarkerSize(flag_size_palette);
  palette_TARGET[9] = new TMarker(0.90,0.075,21);     palette_TARGET[9]->SetMarkerColor(kGreen-0);      palette_TARGET[9]->SetMarkerSize(flag_size_palette);
  TLatex *tex_palette_TARGET_scale = new TLatex(0.905,0.04,"x 100");        tex_palette_TARGET_scale->SetTextSize(0.02);

  tex_palette_TARGET[0] = new TLatex(0.510,0.04,ADC_palette_string_TARGET[0]);      tex_palette_TARGET[0]->SetTextSize(0.02);
  tex_palette_TARGET[1] = new TLatex(0.545,0.04,ADC_palette_string_TARGET[1]);      tex_palette_TARGET[1]->SetTextSize(0.02);
  tex_palette_TARGET[2] = new TLatex(0.585,0.04,ADC_palette_string_TARGET[2]);      tex_palette_TARGET[2]->SetTextSize(0.02);
  tex_palette_TARGET[3] = new TLatex(0.625,0.04,ADC_palette_string_TARGET[3]);      tex_palette_TARGET[3]->SetTextSize(0.02);
  tex_palette_TARGET[4] = new TLatex(0.665,0.04,ADC_palette_string_TARGET[4]);      tex_palette_TARGET[4]->SetTextSize(0.02);
  tex_palette_TARGET[5] = new TLatex(0.705,0.04,ADC_palette_string_TARGET[5]);      tex_palette_TARGET[5]->SetTextSize(0.02);
  tex_palette_TARGET[6] = new TLatex(0.745,0.04,ADC_palette_string_TARGET[6]);      tex_palette_TARGET[6]->SetTextSize(0.02);
  tex_palette_TARGET[7] = new TLatex(0.785,0.04,ADC_palette_string_TARGET[7]);      tex_palette_TARGET[7]->SetTextSize(0.02);
  tex_palette_TARGET[8] = new TLatex(0.825,0.04,ADC_palette_string_TARGET[8]);      tex_palette_TARGET[8]->SetTextSize(0.02);
  tex_palette_TARGET[9] = new TLatex(0.865,0.04,ADC_palette_string_TARGET[9]);      tex_palette_TARGET[9]->SetTextSize(0.02);

  tex_palette_TARGET2[0] = new TLatex(0.510,0.04,ADC_palette_string_TARGET2[0]);      tex_palette_TARGET2[0]->SetTextSize(0.02);
  tex_palette_TARGET2[1] = new TLatex(0.545,0.04,ADC_palette_string_TARGET2[1]);      tex_palette_TARGET2[1]->SetTextSize(0.02);
  tex_palette_TARGET2[2] = new TLatex(0.585,0.04,ADC_palette_string_TARGET2[2]);      tex_palette_TARGET2[2]->SetTextSize(0.02);
  tex_palette_TARGET2[3] = new TLatex(0.625,0.04,ADC_palette_string_TARGET2[3]);      tex_palette_TARGET2[3]->SetTextSize(0.02);
  tex_palette_TARGET2[4] = new TLatex(0.665,0.04,ADC_palette_string_TARGET2[4]);      tex_palette_TARGET2[4]->SetTextSize(0.02);
  tex_palette_TARGET2[5] = new TLatex(0.705,0.04,ADC_palette_string_TARGET2[5]);      tex_palette_TARGET2[5]->SetTextSize(0.02);
  tex_palette_TARGET2[6] = new TLatex(0.745,0.04,ADC_palette_string_TARGET2[6]);      tex_palette_TARGET2[6]->SetTextSize(0.02);
  tex_palette_TARGET2[7] = new TLatex(0.785,0.04,ADC_palette_string_TARGET2[7]);      tex_palette_TARGET2[7]->SetTextSize(0.02);
  tex_palette_TARGET2[8] = new TLatex(0.825,0.04,ADC_palette_string_TARGET2[8]);      tex_palette_TARGET2[8]->SetTextSize(0.02);
  tex_palette_TARGET2[9] = new TLatex(0.865,0.04,ADC_palette_string_TARGET2[9]);      tex_palette_TARGET2[9]->SetTextSize(0.02);

  TLatex *tex_DS_SFT;
  TLatex *tex_US_SFT;
  TLatex *tex_footer_SFT;
  TLatex *tex_version;

  TLatex *tex_number1_SFT;
  TLatex *tex_number2_SFT;
  TLatex *tex_number3_SFT;
  TLatex *tex_number4_SFT;
  TLatex *tex_number5_SFT;
  TLatex *tex_number6_SFT;
  TLatex *tex_number7_SFT;
  TLatex *tex_number8_SFT;
  TLatex *tex_number9_SFT;
  TLatex *tex_number10_SFT;
  TLatex *tex_number11_SFT;
  TLatex *tex_number12_SFT;
  TLatex *tex_number13_SFT;
  TLatex *tex_number14_SFT;
  TLatex *tex_number15_SFT;
  TLatex *tex_number16_SFT;
  TLatex *tex_number17_SFT;

  TLatex *tex_number21_SFT; TLatex *tex_number21_bis_SFT;
  TLatex *tex_number22_SFT; TLatex *tex_number22_bis_SFT;
  TLatex *tex_number23_SFT; TLatex *tex_number23_bis_SFT;
  TLatex *tex_number24_SFT; TLatex *tex_number24_bis_SFT;

  char Event_flag_title[100];
  sprintf(Event_flag_title,"Event Flag = %d", Event_flag[0]);

  TLatex *tex_Event_flag;

  tex_Event_flag = new TLatex(0.02,0.2,Event_flag_title);
  tex_Event_flag->SetTextSize(0.04);
  tex_Event_flag->SetLineWidth(2);

  tex_DS_SFT = new TLatex(0.02,0.88,"SFT  --  Downstream & Upstream");
  tex_DS_SFT->SetTextSize(0.05);
  tex_DS_SFT->SetLineWidth(2);

  tex_US_SFT = new TLatex(0.02,0.41,"SFT  --  Upstream");
  tex_US_SFT->SetTextSize(0.05);
  tex_US_SFT->SetLineWidth(2);

  tex_footer_SFT = new TLatex(0.02,0.08,footer);
  tex_footer_SFT->SetTextSize(0.04);
  tex_footer_SFT->SetLineWidth(2);

  tex_version = new TLatex(0.02,0.02,Version);
  tex_version->SetTextSize(0.04);
  tex_version->SetLineWidth(2);

  tex_number1_SFT = new TLatex(0.115,0.81,"1");
  tex_number1_SFT->SetTextSize(0.025);
  tex_number1_SFT->SetLineWidth(2);

  tex_number2_SFT = new TLatex(0.15,0.81,"2");
  tex_number2_SFT->SetTextSize(0.025);
  tex_number2_SFT->SetLineWidth(2);

  tex_number3_SFT = new TLatex(0.19,0.81,"3");
  tex_number3_SFT->SetTextSize(0.025);
  tex_number3_SFT->SetLineWidth(2);

  tex_number4_SFT = new TLatex(0.23,0.81,"4");
  tex_number4_SFT->SetTextSize(0.025);
  tex_number4_SFT->SetLineWidth(2);

  tex_number5_SFT = new TLatex(0.27,0.81,"5");
  tex_number5_SFT->SetTextSize(0.025);
  tex_number5_SFT->SetLineWidth(2);

  tex_number6_SFT = new TLatex(0.31,0.81,"6");
  tex_number6_SFT->SetTextSize(0.025);
  tex_number6_SFT->SetLineWidth(2);

  tex_number7_SFT = new TLatex(0.35,0.81,"7");
  tex_number7_SFT->SetTextSize(0.025);
  tex_number7_SFT->SetLineWidth(2);

  tex_number8_SFT = new TLatex(0.39,0.81,"8");
  tex_number8_SFT->SetTextSize(0.025);
  tex_number8_SFT->SetLineWidth(2);

  tex_number9_SFT = new TLatex(0.43,0.81,"9");
  tex_number9_SFT->SetTextSize(0.025);
  tex_number9_SFT->SetLineWidth(2);

  tex_number10_SFT = new TLatex(0.465,0.81,"10");
  tex_number10_SFT->SetTextSize(0.025);
  tex_number10_SFT->SetLineWidth(2);

  tex_number11_SFT = new TLatex(0.505,0.81,"11");
  tex_number11_SFT->SetTextSize(0.025);
  tex_number11_SFT->SetLineWidth(2);

  tex_number12_SFT = new TLatex(0.545,0.81,"12");
  tex_number12_SFT->SetTextSize(0.025);
  tex_number12_SFT->SetLineWidth(2);

  tex_number13_SFT = new TLatex(0.585,0.81,"13");
  tex_number13_SFT->SetTextSize(0.025);
  tex_number13_SFT->SetLineWidth(2);

  tex_number14_SFT = new TLatex(0.625,0.81,"14");
  tex_number14_SFT->SetTextSize(0.025);
  tex_number14_SFT->SetLineWidth(2);

  tex_number15_SFT = new TLatex(0.665,0.81,"15");
  tex_number15_SFT->SetTextSize(0.025);
  tex_number15_SFT->SetLineWidth(2);

  tex_number16_SFT = new TLatex(0.705,0.81,"16");
  tex_number16_SFT->SetTextSize(0.025);
  tex_number16_SFT->SetLineWidth(2);

  tex_number17_SFT = new TLatex(0.745,0.81,"17");
  tex_number17_SFT->SetTextSize(0.025);
  tex_number17_SFT->SetLineWidth(2);

  tex_number21_SFT = new TLatex(0.03,0.77,"L1-D");  tex_number21_bis_SFT = new TLatex(0.03,0.73,"L1-U");
  tex_number21_SFT->SetTextSize(0.03);      tex_number21_bis_SFT->SetTextSize(0.03);
  tex_number21_SFT->SetLineWidth(2);        tex_number21_bis_SFT->SetLineWidth(2);

  tex_number22_SFT = new TLatex(0.03,0.63,"L2-D");  tex_number22_bis_SFT = new TLatex(0.03,0.59,"L2-U");
  tex_number22_SFT->SetTextSize(0.03);      tex_number22_bis_SFT->SetTextSize(0.03);
  tex_number22_SFT->SetLineWidth(2);        tex_number22_bis_SFT->SetLineWidth(2);

  tex_number23_SFT = new TLatex(0.03,0.49,"L3-D");  tex_number23_bis_SFT = new TLatex(0.03,0.45,"L3-U");
  tex_number23_SFT->SetTextSize(0.03);      tex_number23_bis_SFT->SetTextSize(0.03);
  tex_number23_SFT->SetLineWidth(2);        tex_number23_bis_SFT->SetLineWidth(2);

  tex_number24_SFT = new TLatex(0.03,0.35,"L4-D");  tex_number24_bis_SFT = new TLatex(0.03,0.31,"L4-U");
  tex_number24_SFT->SetTextSize(0.03);      tex_number24_bis_SFT->SetTextSize(0.03);
  tex_number24_SFT->SetLineWidth(2);        tex_number24_bis_SFT->SetLineWidth(2);

  TLatex *tex_parameter1;
  TLatex *tex_parameter2;
  TLatex *tex_parameter3;
  TLatex *tex_parameter4;
  TLatex *tex_parameter5;
  TLatex *tex_parameter6;
  TLatex *tex_parameter7;
  TLatex *tex_parameter8;

  tex_parameter1 = new TLatex(0.06,0.48,"Run #:");
  tex_parameter1->SetTextSize(0.03);
  tex_parameter1->SetLineWidth(2);

  tex_parameter2 = new TLatex(0.06,0.445,"Event #:");
  tex_parameter2->SetTextSize(0.03);
  tex_parameter2->SetLineWidth(2);

  tex_parameter3 = new TLatex(0.06,0.41,"PHI angle:");
  tex_parameter3->SetTextSize(0.03);
  tex_parameter3->SetLineWidth(2);

  tex_parameter4 = new TLatex(0.06,0.30,"Total Number of Events : ");
  tex_parameter4->SetTextSize(0.03);
  tex_parameter4->SetLineWidth(2);

  tex_parameter5 = new TLatex(0.06,0.26,"Number of hits in L1 or L2 (Helicity 1) : ");
  tex_parameter5->SetTextSize(0.03);
  tex_parameter5->SetLineWidth(2);

  tex_parameter6 = new TLatex(0.06,0.22,"Number of hits in L3 or L4 (Helicity 2) :");
  tex_parameter6->SetTextSize(0.03);
  tex_parameter6->SetLineWidth(2);

  tex_parameter7 = new TLatex(0.06,0.18,"Number of hits in Helicity 1 AND Helicity 2 : ");
  tex_parameter7->SetTextSize(0.03);
  tex_parameter7->SetLineWidth(2);

  tex_parameter8 = new TLatex(0.06,0.12,"Total Efficiency : ");
  tex_parameter8->SetTextSize(0.03);
  tex_parameter8->SetLineWidth(2);


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

  printf("                    Fiber  ADC-Ped   TDC\n");
  printf("LG ADC 1st Max.  :  %4d    %4d    %4d\n",index_max1, ADC_Low_TARGET[index_max1] + TARGET_ADC_Thr_LG_Offset, TDC_LG_max);
  printf("LG ADC 2nd Max.  :  %4d    %4d    %4d\n",index_max2, ADC_Low_TARGET[index_max2] + TARGET_ADC_Thr_LG_Offset, TDC_LG_max2);
  printf("TDC Average  :  %4d\n",TDC_average);

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

  TH1F *h_ADC_L1_DS = new TH1F("h_ADC_L1_DS",h_ADC_title,130,-1,129);
  TH1F *h_ADC_L2_DS = new TH1F("h_ADC_L2_DS",h_ADC_title,130,-1,129);
  TH1F *h_ADC_L3_DS = new TH1F("h_ADC_L3_DS",h_ADC_title,130,-1,129);
  TH1F *h_ADC_L4_DS = new TH1F("h_ADC_L4_DS",h_ADC_title,130,-1,129);
  TH1F *h_ADC_L1_US = new TH1F("h_ADC_L1_US",h_ADC_title,66,63,129);
  TH1F *h_ADC_L2_US = new TH1F("h_ADC_L2_US",h_ADC_title,66,63,129);
  TH1F *h_ADC_L3_US = new TH1F("h_ADC_L3_US",h_ADC_title,66,63,129);
  TH1F *h_ADC_L4_US = new TH1F("h_ADC_L4_US",h_ADC_title,66,63,129);

  TH2F *h_EMPTY = new TH2F("Final Fit", "Final Fit", 500, -50, 50, 500, -50, 50);
  TH2F *h_Circle = new TH2F("Test", "TARGET and TOF1", 500, -50, 50, 500, -50, 50);
  TH2F *h_Target = new TH2F("Test Target", "Target", 3000, -50, 50, 3000, -50, 50);
  TH2F *h_Target_Center = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  TH2F *h_kaon = new TH2F("Kaon", "Kaon", 500, -50, 50, 500, -50, 50);
  TH2F *h_kaon_copy = new TH2F("Kaon Copy", "Kaon Copy", 500, -50, 50, 500, -50, 50);
  TH2F *h_max = new TH2F("Max", "Max", 500, -50, 50, 500, -50, 50);
  TH2F *h_max_copy = new TH2F("Max Copy", "Max Copy", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1 = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_closest = new TH2F("h_TOF1_closest", "h_TOF1_closest", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_rotate = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50); // ROTATE_CHANGE
  TH2F *h_int_TDC_Gap_Fibers = new TH2F("h_int_TDC_Gap_Fibers", "h_int_TDC_Gap_Fibers", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_Gap_Fibers_SFT = new TH2F("h_int_TDC_Gap_Fibers_SFT","h_int_TDC_Gap_Fibers_SFT", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_Gap_Fibers_rotate = new TH2F("h_int_TDC_Gap_Fibers", "h_int_TDC_Gap_Fibers", 500, -50, 50, 500, -50, 50);  // ROTATE_CHANGE
  TH2F *h_int_TDC_Gap_Fibers_SFT_rotate = new TH2F("h_int_TDC_Gap_Fibers_SFT","h_int_TDC_Gap_Fibers_SFT", 500, -50, 50, 500, -50, 50);  // ROTATE_CHANGE
  TH2F *h_int_TDC_TARGET = new TH2F("h_int_TDC_TARGET", "h_int_TDC_TARGET", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_TARGET_rotate = new TH2F("h_int_TDC_TARGET_rotate", "h_int_TDC_TARGET_rotate", 500, -50, 50, 500, -50, 50);
  TH2F *h_Centroid = new TH2F("Centroid", "Centroid", 500, -50, 50, 500, -50, 50);
  TH2F *h_2lines_intersect = new TH2F("2lines_intersect", "2lines_intersect", 500, -50, 50, 500, -50, 50);
  TH2F *h_kaon_stop = new TH2F("kaon_stop", "kaon_stop", 500, -50, 50, 500, -50, 50);

  TEllipse *ell = new TEllipse(0, 0, R_TOF1, 0);
  TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0);
  TEllipse *ell_L1 = new TEllipse(0, 0, R_SFT_L1, 0);

  for(Int_t ii=0; ii<128; ii++){
    if(ADC_High_SFT_corr[ii]>ADC_High_corr_max) ADC_High_corr_max=ADC_High_SFT_corr[ii];
  }

  for(Int_t ii=0; ii<128; ii++){
    for (Int_t qq=0; qq<6; qq++){
      if(tdc_le_sft[ii][qq] > SFT_TDC_min[ii] && tdc_le_sft[ii][qq] < SFT_TDC_max[ii]) has_TDC_SFT_hit[ii]++;
    }
  }

  for(Int_t imark1=0; imark1<15; imark1++){
    if((ADC_High_SFT_corr[imark1]!=0) && (has_TDC_SFT_hit[imark1] > 0)){
      h_ADC_L1_DS->Fill(imark1,ADC_High_SFT_corr[imark1]);
      h_ADC_L1_DS->SetFillColor(2);
    }
    if((ADC_High_SFT_corr[imark1+64]!=0) && (has_TDC_SFT_hit[imark1+64] > 0)){
      h_ADC_L1_US->Fill(imark1+64,ADC_High_SFT_corr[imark1+64]);
      h_ADC_L1_US->SetFillColor(2);
    }
  }

  for(Int_t imark2=0; imark2<15; imark2++){
    if((ADC_High_SFT_corr[imark2+15]!=0) && (has_TDC_SFT_hit[imark2+15] > 0)){
      h_ADC_L2_DS->Fill(imark2+15,ADC_High_SFT_corr[imark2+15]);
      h_ADC_L2_DS->SetFillColor(4);
    }
    if((ADC_High_SFT_corr[imark2+79]!=0) && (has_TDC_SFT_hit[imark2+79] > 0)){
      h_ADC_L2_US->Fill(imark2+79,ADC_High_SFT_corr[imark2+79]);
      h_ADC_L2_US->SetFillColor(4);
    }
  }

  for(Int_t imark3=0; imark3<2; imark3++){
    if((ADC_High_SFT_corr[imark3+30]!=0) && (has_TDC_SFT_hit[imark3+30] > 0)){
      h_ADC_L3_DS->Fill(imark3+30,ADC_High_SFT_corr[imark3+30]);
      h_ADC_L3_DS->SetFillColor(3);
    }
    if((ADC_High_SFT_corr[imark3+94]!=0) && (has_TDC_SFT_hit[imark3+94] > 0)){
      h_ADC_L3_US->Fill(imark3+94,ADC_High_SFT_corr[imark3+94]);
      h_ADC_L3_US->SetFillColor(3);
    }
  }

  for(Int_t imark3=2; imark3<17; imark3++){
    if((ADC_High_SFT_corr[imark3+47]!=0) && (has_TDC_SFT_hit[imark3+47] > 0)){
      h_ADC_L3_DS->Fill(imark3+47,ADC_High_SFT_corr[imark3+47]);
      h_ADC_L3_DS->SetFillColor(3);
    }
    if((ADC_High_SFT_corr[imark3+111]!=0) && (has_TDC_SFT_hit[imark3+111] > 0)){
      h_ADC_L3_US->Fill(imark3+111,ADC_High_SFT_corr[imark3+111]);
      h_ADC_L3_US->SetFillColor(3);
    }
  }

  for(Int_t imark4=0; imark4<17; imark4++){
    if((ADC_High_SFT_corr[imark4+32]!=0) && (has_TDC_SFT_hit[imark4+32] > 0)){
      h_ADC_L4_DS->Fill(imark4+32,ADC_High_SFT_corr[imark4+32]);
      h_ADC_L4_DS->SetFillColor(1);
    }
    if((ADC_High_SFT_corr[imark4+96]!=0) && (has_TDC_SFT_hit[imark4+96] > 0)){
      h_ADC_L4_US->Fill(imark4+96,ADC_High_SFT_corr[imark4+96]);
      h_ADC_L4_US->SetFillColor(1);
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

  if (has_TDC_TOF2_hit[0] || has_ADC_TOF2_hit[0]) TOF_line13->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[0] && has_ADC_TOF2_hit[0]) TOF_line13->SetLineColor(3);

  if (has_TDC_TOF2_hit[1] || has_ADC_TOF2_hit[1]) TOF_line14->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[1] && has_ADC_TOF2_hit[1]) TOF_line14->SetLineColor(3);

  if (has_TDC_TOF2_hit[2] || has_ADC_TOF2_hit[2]) TOF_line15->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[2] && has_ADC_TOF2_hit[2]) TOF_line15->SetLineColor(3);

  if (has_TDC_TOF2_hit[3] || has_ADC_TOF2_hit[3]) TOF_line16->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[3] && has_ADC_TOF2_hit[3]) TOF_line16->SetLineColor(3);

  if (has_TDC_TOF2_hit[4] || has_ADC_TOF2_hit[4]) TOF_line17->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[4] && has_ADC_TOF2_hit[4]) TOF_line17->SetLineColor(3);

  if (has_TDC_TOF2_hit[5] || has_ADC_TOF2_hit[5]) TOF_line18->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[5] && has_ADC_TOF2_hit[5]) TOF_line18->SetLineColor(3);

  if (has_TDC_TOF2_hit[6] || has_ADC_TOF2_hit[6]) TOF_line19->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[6] && has_ADC_TOF2_hit[6]) TOF_line19->SetLineColor(3);

  if (has_TDC_TOF2_hit[7] || has_ADC_TOF2_hit[7]) TOF_line20->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[7] && has_ADC_TOF2_hit[7]) TOF_line20->SetLineColor(3);

  if (has_TDC_TOF2_hit[8] || has_ADC_TOF2_hit[8]) TOF_line21->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[8] && has_ADC_TOF2_hit[8]) TOF_line21->SetLineColor(3);

  if (has_TDC_TOF2_hit[9] || has_ADC_TOF2_hit[9]) TOF_line22->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[9] && has_ADC_TOF2_hit[9]) TOF_line22->SetLineColor(3);

  if (has_TDC_TOF2_hit[10] || has_ADC_TOF2_hit[10]) TOF_line23->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[10] && has_ADC_TOF2_hit[10]) TOF_line23->SetLineColor(3);

  if (has_TDC_TOF2_hit[11] || has_ADC_TOF2_hit[11]) TOF_line24->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[11] && has_ADC_TOF2_hit[11]) TOF_line24->SetLineColor(3);


  ///Set TOF1 Lines
  for(int i = 0; i < 12; i++){
    if (ADC_TOF1U[i]>0 && ADC_TOF1D[i]>0) {has_both_ADC_TOF1_hit[i] = true;}
    if (ADC_TOF1U[i]>0 || ADC_TOF1D[i]>0) {has_ADC_TOF1_hit[i] = true;}
    if ((TDC_TOF1U[i]>TOF1U_TDC_min[i] && TDC_TOF1U[i]<TOF1U_TDC_max[i]) && (TDC_TOF1D[i]>TOF1D_TDC_min[i] && TDC_TOF1D[i]<TOF1D_TDC_max[i])) {has_both_TDC_TOF1_hit[i] = true;}
    if ((TDC_TOF1U[i]>TOF1U_TDC_min[i] && TDC_TOF1U[i]<TOF1U_TDC_max[i]) || (TDC_TOF1D[i]>TOF1D_TDC_min[i] && TDC_TOF1D[i]<TOF1D_TDC_max[i])) {has_TDC_TOF1_hit[i] = true;}
  }

  if (has_both_ADC_TOF1_hit[0]) TOF_line1->SetLineColor(kRed);
  if (has_both_TDC_TOF1_hit[0]) TOF_line1->SetLineColor(kYellow);
  if (has_ADC_TOF1_hit[0] || has_TDC_TOF1_hit[0]) TOF_line1->SetLineColor(kOrange);
  if (has_both_TDC_TOF1_hit[0]) TOF_line1->SetLineColor(3);

  if (has_TDC_TOF1_hit[0]) TOF_line1->SetLineColor(kYellow);
  if (has_ADC_TOF1_hit[0] && has_TDC_TOF1_hit[0]) TOF_line1->SetLineColor(3);

  if (has_ADC_TOF1_hit[1]) TOF_line2->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[1]) TOF_line2->SetLineColor(kYellow);
  if (has_ADC_TOF1_hit[1] && has_TDC_TOF1_hit[1]) TOF_line2->SetLineColor(3);

  if (has_ADC_TOF1_hit[2]) TOF_line3->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[2]) TOF_line3->SetLineColor(kYellow);
  if (has_ADC_TOF1_hit[2] && has_TDC_TOF1_hit[2]) TOF_line3->SetLineColor(3);

  if (has_ADC_TOF1_hit[3]) TOF_line4->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[3]) TOF_line4->SetLineColor(kYellow);
  if (has_ADC_TOF1_hit[3] && has_TDC_TOF1_hit[3]) TOF_line4->SetLineColor(3);

  if (has_ADC_TOF1_hit[4]) TOF_line5->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[4]) TOF_line5->SetLineColor(kYellow);
  if (has_ADC_TOF1_hit[4] && has_TDC_TOF1_hit[4]) TOF_line5->SetLineColor(3);

  if (has_ADC_TOF1_hit[5]) TOF_line6->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[5]) TOF_line6->SetLineColor(kYellow);
  if (has_ADC_TOF1_hit[5] && has_TDC_TOF1_hit[5]) TOF_line6->SetLineColor(3);

  if (has_ADC_TOF1_hit[6]) TOF_line7->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[6]) TOF_line7->SetLineColor(kYellow);
  if (has_ADC_TOF1_hit[6] && has_TDC_TOF1_hit[6]) TOF_line7->SetLineColor(3);

  if (has_ADC_TOF1_hit[7]) TOF_line8->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[7]) TOF_line8->SetLineColor(kYellow);
  if (has_ADC_TOF1_hit[7] && has_TDC_TOF1_hit[7]) TOF_line8->SetLineColor(3);

  if (has_ADC_TOF1_hit[8]) TOF_line9->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[8]) TOF_line9->SetLineColor(kYellow);
  if (has_ADC_TOF1_hit[8] && has_TDC_TOF1_hit[8]) TOF_line9->SetLineColor(3);

  if (has_ADC_TOF1_hit[9]) TOF_line10->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[9]) TOF_line10->SetLineColor(kYellow);
  if (has_ADC_TOF1_hit[9] && has_TDC_TOF1_hit[9]) TOF_line10->SetLineColor(3);

  if (has_ADC_TOF1_hit[10]) TOF_line11->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[10]) TOF_line11->SetLineColor(kYellow);
  if (has_ADC_TOF1_hit[10] && has_TDC_TOF1_hit[10]) TOF_line11->SetLineColor(3);

  if (has_ADC_TOF1_hit[11]) TOF_line12->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[11]) TOF_line12->SetLineColor(kYellow);
  if (has_ADC_TOF1_hit[11] && has_TDC_TOF1_hit[11]) TOF_line12->SetLineColor(3);

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


  //double total_energy = 0.0;
  //double X_weights = 0.0;
  //double Y_weights = 0.0;

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
      //cout << "////////////////////////////   DEBUG1 : " << i << endl;
      vec_xx_kaon.push_back(Xloc[i]);
      vec_ex_kaon.push_back(TARGET_Errors_X);

      vec_yy_kaon.push_back(Yloc[i]);
      vec_ey_kaon.push_back(TARGET_Errors_Y);

      kaon_bars.push_back(i);
    }
  }

  //delete gr_kaon;               gr_kaon = nullptr;

  TGraph *gr_kaon = new TGraphErrors(vec_xx_kaon.size(),&vec_xx_kaon[0],&vec_yy_kaon[0],
    &vec_ey_kaon[0],&vec_ey_kaon[0]);
  gr_kaon->SetMarkerStyle(21);
  gr_kaon->SetMarkerColor(kBlue-6);
  gr_kaon->SetMarkerSize(0.8);
  gr_kaon->GetXaxis()->SetLimits(-50.,50.);
  gr_kaon->GetYaxis()->SetRangeUser(-50.,50.);

  TGraph *gr_kaon_bk;
  gr_kaon_bk = new TGraph(vec_xx_kaon.size(),&vec_xx_kaon[0],&vec_yy_kaon[0]);

  TF1 *gr_kaon_fit = new TF1("kaon_fit", "pol1");
  double a_fit_kaon = 0.;
  double b_fit_kaon = 0.;
  bool kaon_bk = false;
  //int a_fit_kaon = 0;
  //int b_fit_kaon = 0;

  //cout << "SEB1 : " << a_fit_kaon << "  " << b_fit_kaon << endl;
/*    cout << endl;
    cout << endl;
    cout << endl;
    cout << endl;
    cout << endl;
    cout << endl;
*/
  if(vec_xx_kaon.size()>=5){
    //gr_kaon->Fit("pol1","QCM");
    //gr_kaon_fit->SetParameter(0,0);
    //gr_kaon_fit->SetParameter(1,1);
    //gr_kaon_fit->SetParLimits(0,-50,50);
    //gr_kaon_fit->SetParLimits(0,-50,50);
    gr_kaon->Fit("kaon_fit","Q");
    gr_kaon_fit = gr_kaon->GetFunction("kaon_fit");
    gr_kaon_fit->SetLineWidth(2);
    gr_kaon_fit->SetLineColor(4);
    a_fit_kaon = gr_kaon_fit->GetParameter(1);
    b_fit_kaon = gr_kaon_fit->GetParameter(0);
    //cout << "SEB3 : " << gr_kaon_fit->GetParameter(1) << "  " << gr_kaon_fit->GetParameter(0) << endl;
  }
  //cout << "SEB2 : " << a_fit_kaon << "  " << b_fit_kaon << endl;

  if(a_fit_kaon > 1000){
    kaon_bk = true;

    gr_kaon_bk->SetMarkerStyle(21);
    gr_kaon_bk->SetMarkerColor(kBlue-6);
    gr_kaon_bk->SetMarkerSize(0.8);
    gr_kaon_bk->GetXaxis()->SetLimits(-50.,50.);
    gr_kaon_bk->GetYaxis()->SetRangeUser(-50.,50.);

    gr_kaon_bk->Fit("kaon_fit","Q");
    gr_kaon_fit = gr_kaon_bk->GetFunction("kaon_fit");
    gr_kaon_fit->SetLineWidth(2);
    gr_kaon_fit->SetLineColor(4);
    a_fit_kaon = gr_kaon_fit->GetParameter(1);
    b_fit_kaon = gr_kaon_fit->GetParameter(0);
  }

  //cout << "SEB2 : " << a_fit_kaon << "  " << b_fit_kaon << endl;


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
    //if(TARGET_High_has_neighbours[i] && distance_to_line(Xloc[i],Yloc[i],a_fit_kaon,b_fit_kaon) >= dist_to_kaon_fit){
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

        //cout << "//// LEPTON : " << vec_xx_lepton.size() << endl;

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

              //lepton_bars.push_back(i);
              //lepton_bars.push_back(i);
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

              //lepton_bars.push_back(i);
              //lepton_bars.push_back(i);
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
          //cout << "//// BLUE BARS : " << vec_xx_lepton.size() << endl;
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

 // cout << "//// TOF1 : " << vec_xx_lepton.size() << endl;

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

  //cout << "//// TOF1 : " << vec_xx_lepton.size() << endl;


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TGraph *gr_lepton_1 = new TGraphErrors(vec_xx_lepton.size(),&vec_xx_lepton[0],&vec_yy_lepton[0],
    &vec_ex_lepton[0],&vec_ey_lepton[0]);

  gr_lepton_1->SetMarkerStyle(21);
  gr_lepton_1->SetMarkerColor(2);
  gr_lepton_1->SetMarkerSize(0.8);
  gr_lepton_1->GetXaxis()->SetLimits(-50.,50.);
  gr_lepton_1->GetYaxis()->SetRangeUser(-50.,50.);

  double a_lepton_fit_1 = 0.;
  double b_lepton_fit_1 = 0.;
  double Chis_lepton_fit_1 = 0.;
  int ndf_lepton_fit_1 = 0;

  TF1 *func_lepton_fit_1 = new TF1("Lepton_fit_1", "pol1");
  //TF1 *gr_test_kaon_fit = new TF1("gr_test_kaon_fit", "pol1");

  //if(vec_xx_lepton.size()>0){
    //gr_test_lepton->Fit("pol1","QCM");

    if(gap_to_fit==1 || gap_to_fit==2 || gap_to_fit==7 || gap_to_fit==8){
      func_lepton_fit_1->SetParameter(0,0);
      func_lepton_fit_1->SetParameter(1,1);
    }
    if(gap_to_fit==4 || gap_to_fit==5 || gap_to_fit==10 || gap_to_fit==11){
      func_lepton_fit_1->SetParameter(0,0);
      func_lepton_fit_1->SetParameter(1,-1);
    }
    else{
      func_lepton_fit_1->SetParameter(0,0);
      func_lepton_fit_1->SetParameter(1,1);
    }

    func_lepton_fit_1->SetParLimits(0,-50,50);
    func_lepton_fit_1->SetParLimits(1,-50,50);

    gr_lepton_1->Fit("Lepton_fit_1","Q");
    func_lepton_fit_1 = gr_lepton_1->GetFunction("Lepton_fit_1");
    func_lepton_fit_1->SetLineWidth(2);
    func_lepton_fit_1->SetLineColor(2);
    a_lepton_fit_1 = func_lepton_fit_1->GetParameter(1);
    b_lepton_fit_1 = func_lepton_fit_1->GetParameter(0);
    Chis_lepton_fit_1 = func_lepton_fit_1->GetChisquare();
    ndf_lepton_fit_1 = func_lepton_fit_1->GetNDF();
  //}

  /*cout << endl;
  cout << "//////////////////////  FIT 1" << endl;
  cout << endl;
  cout << "a_lepton_1 : " << a_lepton_fit_1 << endl;
  cout << "b_lepton_1 : " << b_lepton_fit_1 << endl;
  cout << "ChiS_lepton_1 : " << Chis_lepton_fit_1 << endl;
  cout << "ndf_lepton_1 : " << ndf_lepton_fit_1 << endl;
  cout << endl;
*/

  //if(vec_xx_kaon.size()>0){
    //gr_kaon->Fit("pol1","QCM");
    //gr_kaon->Fit("pol1","QC");
    //gr_kaon_fit = gr_test_kaon->GetFunction("pol1");
    //gr_kaon_fit->SetLineWidth(2);
    //gr_test_kaon_fit->SetLineColor(4);
    //a_test_kaon_fit = gr_test_kaon_fit->GetParameter(1);
  //}



  //cout << "ALBANE INTIALE : " << a_test_lepton_fit << "   " << a_test_kaon_fit << endl;

  //cout << "ALBANE TEST : " << a_lepton_fit_1 << "   " << a_fit_kaon << endl;


////////////////////////////////////////////////////////////////////////////////////////////////////////

/*    vec_xx_lepton.clear();
    vec_ex_lepton.clear();
    vec_yy_lepton.clear();
    vec_ey_lepton.clear();

    //vec_xx_lepton_test.clear();
    //vec_yy_lepton_test.clear();
    vec_lepton_size.clear();



  //if(a_test_lepton_fit/a_test_kaon_fit<0 && vec_lepton_size.size()>5){
  if(a_lepton_fit_1/a_fit_kaon < 0){

    if(gap_to_fit == 6 || gap_to_fit == 12){

      for(Int_t i=0; i<256; i++){
        if(TARGET_High_has_neighbours[i] &&
           distance_to_line(Xloc[TARGET_Rotated_index[i]],Yloc[TARGET_Rotated_index[i]],a_fit_kaon,b_fit_kaon) >= dist_to_kaon_fit){
          if(ADC_High_TARGET[i]>Angle_ADC_cut &&
          has_TDC_hit[i] &&
          !k_stop_bar[i]){

            vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[i]]);
            vec_ex_lepton.push_back(TARGET_Errors_X);

            vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[i]]);
            vec_ey_lepton.push_back(TARGET_Errors_Y);

            count++;

            if(gap_to_fit == 6){

              if(IsIn(TARGET_Rotated_index[i],channel[2][0], channel[2][1],   // Additional weight for fibers close to the edge if hit
                                              channel[2][2], channel[2][3],
                                              channel[2][4], channel[2][5],
                                              channel[2][6], channel[2][7])){

                vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[i]]);
                vec_ex_lepton.push_back(TARGET_Errors_X);
                vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[i]]);
                vec_ex_lepton.push_back(TARGET_Errors_X);

                vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[i]]);
                vec_ey_lepton.push_back(TARGET_Errors_Y);
                vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[i]]);
                vec_ey_lepton.push_back(TARGET_Errors_Y);
              }
            }

            if(gap_to_fit == 12){
              if(IsIn(TARGET_Rotated_index[i],channel[8][0], channel[8][1],   // Additional weight for fibers close to the edge if hit
                                              channel[8][2], channel[8][3],
                                              channel[8][4], channel[8][5],
                                              channel[8][6], channel[8][7])){

                vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[i]]);
                vec_ex_lepton.push_back(TARGET_Errors_X);
                vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[i]]);
                vec_ex_lepton.push_back(TARGET_Errors_X);

                vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[i]]);
                vec_ey_lepton.push_back(TARGET_Errors_Y);
                vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[i]]);
                vec_ey_lepton.push_back(TARGET_Errors_Y);
              }
            }
          }

          if(ADC_High_TARGET[TARGET_Rotated_index[i]]>Angle_ADC_cut){
            x_inc = x_inc + Xloc[TARGET_Rotated_index[i]];
            y_inc = y_inc + Yloc[TARGET_Rotated_index[i]];
            hit_count++;
          }

          if(ADC_High_TARGET[TARGET_Rotated_index[i]]<0 && ADC_Low_TARGET[TARGET_Rotated_index[i]]>0 && Switch==1){
            if(!k_stop_bar[TARGET_Rotated_index[i]]){
              cout << "TEST LOOP 3" << endl;

              vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[i]]);
              vec_ex_lepton.push_back(TARGET_Errors_X);

              vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[i]]);
              vec_ey_lepton.push_back(TARGET_Errors_Y);

            }
            count++;
          }
        }
      }
    }
    else{
      for(Int_t i=0; i<256; i++){

        if(TARGET_High_has_neighbours[i] && distance_to_line(Xloc[i],Yloc[i],a_fit_kaon,b_fit_kaon) >= dist_to_kaon_fit){
        //if(TARGET_High_has_neighbours[i]){
          if(ADC_High_TARGET[i]>Angle_ADC_cut && has_TDC_hit[i]){
            if(!k_stop_bar[i]){

              cout << "TEST LOOP 1" << endl;

              vec_xx_lepton.push_back(Xloc[i]);
              vec_ex_lepton.push_back(TARGET_Errors_X);

              vec_yy_lepton.push_back(Yloc[i]);
              vec_ey_lepton.push_back(TARGET_Errors_Y);

            }
            count++;

            if(IsIn(i,channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],   // Additional weight for fibers close to the edge if hit
                    channel[gap_to_fit-1][2], channel[gap_to_fit-1][3],
                    channel[gap_to_fit-1][4], channel[gap_to_fit-1][5],
                    channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){
              if(!k_stop_bar[i]){

                cout << "TEST LOOP 2" << endl;
                vec_xx_lepton.push_back(Xloc[i]);
                vec_ex_lepton.push_back(TARGET_Errors_X);
                vec_xx_lepton.push_back(Xloc[i]);
                vec_ex_lepton.push_back(TARGET_Errors_X);

                vec_yy_lepton.push_back(Yloc[i]);
                vec_ey_lepton.push_back(TARGET_Errors_Y);
                vec_yy_lepton.push_back(Yloc[i]);
                vec_ey_lepton.push_back(TARGET_Errors_Y);
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
              cout << "TEST LOOP 3" << endl;

              vec_xx_lepton.push_back(Xloc[i]);
              vec_ex_lepton.push_back(TARGET_Errors_X);

              vec_yy_lepton.push_back(Yloc[i]);
              vec_ey_lepton.push_back(TARGET_Errors_Y);
            }
            count++;
          }
        }
      }
    }
  }

  cout << "/////////////////////////////////////////////  SIZE GR2 : " << vec_xx_lepton.size() << endl;

 */ /////////////////////////////////////////////////////////////////////////////////////////////////////////////////




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

  /// Print Values
  if(Good_Event && enable_cout==1){
    if(!NoFile) cout << "//////  TARGET  (WALK CORRECTION APPLIED !)  //////" << endl;
    if(NoFile) cout << "//////  TARGET  //////" << endl;
    cout << "Event: " << ievt << "     ADC Threshold Offset (HG): " << TARGET_ADC_Thr_HG_Offset << "     ADC Threshold Offset (LG): " << TARGET_ADC_Thr_LG_Offset << endl << endl;
    cout << "Fiber  HG-Ped  LG-Ped    TDC[0]   TDC[1]   TDC[2]   TDC[3]   TDC[4]   TDC[5]   Thr(HG)   Thr(LG)   HG/LG" << endl;  // SEBASTIEN
    for(Int_t jj=0; jj<256; jj++){
      if(TDC_LE_TARGET[jj]>-1 || tdc_le_target[jj][1]>-1 || tdc_le_target[jj][2]>-1 || tdc_le_target[jj][3]>-1){
        printf("%3d     %4d    %4d     %4d     %4d     %4d     %4d     %4d     %4d      %4d      %4d     %2.2f\n",
       jj, ADC_High_TARGET[jj] + TARGET_ADC_Thr_HG_Offset, ADC_Low_TARGET[jj] + TARGET_ADC_Thr_LG_Offset,
       TDC_LE_TARGET[jj], tdc_le_target[jj][1], tdc_le_target[jj][2], tdc_le_target[jj][3], tdc_le_target[jj][4], tdc_le_target[jj][5],
       HG_TARGET_ADC_Thr[jj], LG_TARGET_ADC_Thr[jj], float(adc_high_target[jj]-round(TARGET_ADC_Ped_HG[jj]))/float(adc_low_target[jj]-round(TARGET_ADC_Ped_LG[jj])));
      }
    }
  }

  if(Good_Event && (enable_cout==0 || enable_cout==9)){
    if(!NoFile) cout << "//////  TARGET  (WALK CORRECTION APPLIED !)  //////" << endl;
    if(NoFile) cout << "//////  TARGET  //////" << endl;
    cout << "Event: " << ievt << "     ADC Threshold Offset (HG): " << TARGET_ADC_Thr_HG_Offset << "     ADC Threshold Offset (LG): " << TARGET_ADC_Thr_LG_Offset << endl << endl;
    cout << "Fiber  HG-Ped   LG-Ped   T[0]-Av   T[1]-Av   T[2]-Av   T[3]-Av   T[4]-Av   T[5]-Av   Thr(HG)   Thr(LG)   HG/LG" << endl;
      for(Int_t jj=0; jj<256; jj++){
      if(ADC_High_TARGET[jj]>0 || (ADC_Low_TARGET[jj]>0 && Switch==1)){
        printf("%3d     %4d    %4d      %s       %s      %s      %s      %s      %s     %4d      %4d     %2.2f\n",
        jj, ADC_High_TARGET[jj]+TARGET_ADC_Thr_HG_Offset, ADC_Low_TARGET[jj] + TARGET_ADC_Thr_LG_Offset,
        TDC_LE_TARGET_corr[jj][0], TDC_LE_TARGET_corr[jj][1], TDC_LE_TARGET_corr[jj][2], TDC_LE_TARGET_corr[jj][3], TDC_LE_TARGET_corr[jj][4], TDC_LE_TARGET_corr[jj][5],
        HG_TARGET_ADC_Thr[jj], LG_TARGET_ADC_Thr[jj], float(adc_high_target[jj]-round(TARGET_ADC_Ped_HG[jj]))/float(adc_low_target[jj]-round(TARGET_ADC_Ped_LG[jj])));
      }
    }
  }

  if(Good_Event && enable_cout==2){
    if(!NoFile) cout << "//////  TARGET  (WALK CORRECTION APPLIED !)  //////" << endl;
    if(NoFile) cout << "//////  TARGET  //////" << endl;
    cout << "Event: " << ievt << "     ADC Threshold Offset (HG): " << TARGET_ADC_Thr_HG_Offset << "     ADC Threshold Offset (LG): " << TARGET_ADC_Thr_LG_Offset << endl << endl;
    cout << "Fiber  HG-Ped   LG-Ped   T[0]-Av   T[1]-Av   T[2]-Av   T[3]-Av   T[4]-Av   T[5]-Av   Thr(HG)   Thr(LG)   HG/LG" << endl;   // SEBASTIEN
    for(Int_t jj=0; jj<256; jj++){
        printf("%3d     %4d     %4d      %4d      %4d      %4d      %4d      %4d      %4d     %4d      %4d     %2.2f\n",
        jj, ADC_High_TARGET[jj]+TARGET_ADC_Thr_HG_Offset, ADC_Low_TARGET[jj] + TARGET_ADC_Thr_LG_Offset,
        TDC_LE_TARGET_corrected[jj][0], TDC_LE_TARGET_corrected[jj][1], TDC_LE_TARGET_corrected[jj][2], TDC_LE_TARGET_corrected[jj][3], TDC_LE_TARGET_corrected[jj][4], TDC_LE_TARGET_corrected[jj][5],
        HG_TARGET_ADC_Thr[jj], LG_TARGET_ADC_Thr[jj], float(adc_high_target[jj]-round(TARGET_ADC_Ped_HG[jj]))/float(adc_low_target[jj]-round(TARGET_ADC_Ped_LG[jj])));
   }
  }

  if(Good_Event && (enable_cout==0 || enable_cout==1 || enable_cout==9)){
    cout << " " << endl;
    cout << " " << endl;
    cout << "//////  SFT  //////" << endl;
    cout << "ADC Threshold Offset (HG): " << SFT_ADC_Thr_HG_Offset << endl;
    cout << "ADC Threshold Offset (LG): " << SFT_ADC_Thr_LG_Offset << endl << endl;
    cout << "Channel   Layer     Fiber     HG-Thr     LG-Thr   TDC[0]    TDC[1]    TDC[2]    TDC[3]   (HG) Thr   (LG) Thr" << endl;
    for(Int_t jj=0; jj<128; jj++){
      if((ADC_High_SFT_corr[jj] != 0) && (has_TDC_SFT_hit[jj] > 0)){
        printf("%3d       %3d      %3d-%1c      %4d       %4d      %4d     %4d      %4d      %4d       %4d       %4d\n",
        jj, SFT_channel_to_layer[jj], SFT_channel_to_fiber[jj], SFT_channel_to_US_DS[jj], ADC_High_SFT[jj], ADC_Low_SFT[jj],
        TDC_LE_SFT[jj], tdc_le_sft[jj][1], tdc_le_sft[jj][2], tdc_le_sft[jj][3], HG_SFT_ADC_Thr[jj], HG_SFT_ADC_Thr[jj]);
      }
    }
  }

  if(Good_Event && enable_cout==2){
    cout << " " << endl;
    cout << " " << endl;
    cout << "//////  SFT  //////" << endl << endl;
    cout << "ADC Threshold Offset (HG): " << SFT_ADC_Thr_HG_Offset << endl;
    cout << "ADC Threshold Offset (LG): " << SFT_ADC_Thr_LG_Offset << endl << endl;
    cout << "Channel   Layer     Fiber     HG-Thr     LG-Thr   TDC[0]    TDC[1]    TDC[2]    TDC[3]    Thr(HG)    Thr(LG)" << endl;
    for(Int_t jj=0; jj<128; jj++){
      printf("%3d       %3d      %3d-%1c      %4d       %4d      %4d     %4d      %4d      %4d       %4d       %4d\n",
      jj, SFT_channel_to_layer[jj], SFT_channel_to_fiber[jj], SFT_channel_to_US_DS[jj], ADC_High_SFT[jj], ADC_Low_SFT[jj],
      TDC_LE_SFT[jj], tdc_le_sft[jj][1], tdc_le_sft[jj][2], tdc_le_sft[jj][3], HG_SFT_ADC_Thr[jj], HG_SFT_ADC_Thr[jj]);

    }
  }

  cout << "" << endl;
  cout << "" << endl;

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

  if(gap_to_fit >= 7){ //LEFT
    for(int q=0; q<56; q++){ // C2 Counters
      if(Good_Event && ADC_C2X_L[q]>0){
        cout << " ADC MWPC Channel " << C2XL_Channels[q] << " -- " << "C2_X_" << q+1 << "_Left: " << int(ADC_C2X_L[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){
      if(Good_Event && ADC_C2Y_L[q]>0){
        cout << " ADC MWPC Channel " << C2YL_Channels[q] << " -- " << "C2_Y_" << q+1 << "_Left: " << int(ADC_C2Y_L[q]) << endl;
      }
    }
    for(int q=0; q<56; q++){ // C2 Counters
      if(Good_Event && ADC_C2X_R[q]>0){
        cout << " ADC MWPC Channel " << C2XR_Channels[q] << " -- " << "C2_X_" << q+1 << "_Right: " << int(ADC_C2X_R[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){
      if(Good_Event && ADC_C2Y_R[q]>0){
        cout << " ADC MWPC Channel " << C2YR_Channels[q] << " -- " << "C2_Y_" << q+1 << "_Right: " << int(ADC_C2Y_R[q]) << endl;
      }
    }

    cout << " " << endl;

    for(int q=0; q<64; q++){ // C3 Counters
      if(Good_Event && ADC_C3X_L[q]>0){
        cout << " ADC MWPC Channel " << C3XL_Channels[q] << " -- " << "C3_X_" << q+1 << "_Left: " << int(ADC_C3X_L[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){
      if(Good_Event && ADC_C3Y_L[q]>0){
        cout << " ADC MWPC Channel " << C3YL_Channels[q] << " -- " << "C3_Y_" << q+1 << "_Left: " << int(ADC_C3Y_L[q]) << endl;
      }
    }
    for(int q=0; q<64; q++){ // C3 Counters
      if(Good_Event && ADC_C3X_R[q]>0){
        cout << " ADC MWPC Channel " << C3XR_Channels[q] << " -- " << "C3_X_" << q+1 << "_Right: " << int(ADC_C3X_R[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){ // C3 Counters
      if(Good_Event && ADC_C3Y_R[q]>0){
        cout << " ADC MWPC Channel " << C3YR_Channels[q] << " -- " << "C3_Y_" << q+1 << "_Right: " << int(ADC_C3Y_R[q]) << endl;
      }
    }

    cout << " " << endl;

    for(int q=0; q<72; q++){ // C4 Counters
      if(Good_Event && ADC_C4X_L[q]>0){
        cout << " ADC MWPC Channel " << C4XL_Channels[q] << " -- " << "C4_X_" << q+1 << "_Left: " << int(ADC_C4X_L[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){
      if(Good_Event && ADC_C4Y_L[q]>0){
        cout << " ADC MWPC Channel " << C4YL_Channels[q] << " -- " << "C4_Y_" << q+1 << "_Left: " << int(ADC_C4Y_L[q]) << endl;
      }
    }
    for(int q=0; q<72; q++){ // C4 Counters
      if(Good_Event && ADC_C4X_R[q]>0){
        cout << " ADC MWPC Channel " << C4XR_Channels[q] << " -- " << "C4_X_" << q+1 << "_Right: " << int(ADC_C4X_R[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){
      if(Good_Event && ADC_C4Y_R[q]>0){
        cout << " ADC MWPC Channel " << C4YR_Channels[q] << " -- " << "C4_Y_" << q+1 << "_Right: " << int(ADC_C4Y_R[q]) << endl;
      }
    }
  }
  else{
    for(int q=0; q<56; q++){ // C2 Counters
      if(Good_Event && ADC_C2X_R[q]>0){
        cout << " ADC MWPC Channel " << C2XR_Channels[q] << " -- " << "C2_X_" << q+1 << "_Right: " << int(ADC_C2X_R[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){
      if(Good_Event && ADC_C2Y_R[q]>0){
        cout << " ADC MWPC Channel " << C2YR_Channels[q] << " -- " << "C2_Y_" << q+1 << "_Right: " << int(ADC_C2Y_R[q]) << endl;
      }
    }
    for(int q=0; q<56; q++){ // C2 Counters
      if(Good_Event && ADC_C2X_L[q]>0){
        cout << " ADC MWPC Channel " << C2XL_Channels[q] << " -- " << "C2_X_" << q+1 << "_Left: " << int(ADC_C2X_L[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){
      if(Good_Event && ADC_C2Y_L[q]>0){
        cout << " ADC MWPC Channel " << C2YL_Channels[q] << " -- " << "C2_Y_" << q+1 << "_Left: " << int(ADC_C2Y_L[q]) << endl;
      }
    }

    cout << " " << endl;

    for(int q=0; q<64; q++){ // C3 Counters
      if(Good_Event && ADC_C3X_R[q]>0){
        cout << " ADC MWPC Channel " << C3XR_Channels[q] << " -- " << "C3_X_" << q+1 << "_Right: " << int(ADC_C3X_R[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){
      if(Good_Event && ADC_C3Y_R[q]>0){
        cout << " ADC MWPC Channel " << C3YR_Channels[q] << " -- " << "C3_Y_" << q+1 << "_Right: " << int(ADC_C3Y_R[q]) << endl;
      }
    }
    for(int q=0; q<64; q++){ // C3 Counters
      if(Good_Event && ADC_C3X_L[q]>0){
        cout << " ADC MWPC Channel " << C3XL_Channels[q] << " -- " << "C3_X_" << q+1 << "_Left: " << int(ADC_C3X_L[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){
      if(Good_Event && ADC_C3Y_L[q]>0){
        cout << " ADC MWPC Channel " << C3YL_Channels[q] << " -- " << "C3_Y_" << q+1 << "_Left: " << int(ADC_C3Y_L[q]) << endl;
      }
    }

    cout << " " << endl;

    for(int q=0; q<72; q++){ // C4 Counters
      if(Good_Event && ADC_C4X_R[q]>0){
        cout << " ADC MWPC Channel " << C4XR_Channels[q] << " -- " << "C4_X_" << q+1 << "_Right: " << int(ADC_C4X_R[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){
      if(Good_Event && ADC_C4Y_R[q]>0){
        cout << " ADC MWPC Channel " << C4YR_Channels[q] << " -- " << "C4_Y_" << q+1 << "_Right: " << int(ADC_C4Y_R[q]) << endl;
      }
    }
    for(int q=0; q<72; q++){ // C4 Counters
      if(Good_Event && ADC_C4X_L[q]>0){
        cout << " ADC MWPC Channel " << C4XL_Channels[q] << " -- " << "C4_X_" << q+1 << "_Left: " << int(ADC_C4X_L[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){
      if(Good_Event && ADC_C4Y_L[q]>0){
        cout << " ADC MWPC Channel " << C4YL_Channels[q] << " -- " << "C4_Y_" << q+1 << "_Left: " << int(ADC_C4Y_L[q]) << endl;
      }
    }
  }

  cout << endl;
  cout << endl;

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

  // Display MWPC counter information
  cout << endl;
  if(selected_TOF2 > 6) cout << " LEFT" << endl;
  else cout << " RIGHT" << endl;

  cout << " C2X clusters = " << C2X_clusters << endl;
  cout << " C2Y clusters = " << C2Y_clusters << endl;
  cout << " C3X clusters = " << C3X_clusters << endl;
  cout << " C3Y clusters = " << C3Y_clusters << endl;
  cout << " C4X clusters = " << C4X_clusters << endl;
  cout << " C4Y clusters = " << C4Y_clusters << endl;

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
  ////////////////////////////////////////////////////////////////
  if(MWPC_switch==1){

    TCanvas *cMWPC;
    cMWPC = new TCanvas("Wire Chamber","Wire Chamber",50,50,1050,700);
    cMWPC->Divide(1,3);

    int C2X_count = 0;
    int C2Y_count = 0;
    int C3X_count = 0;
    int C3Y_count = 0;
    int C4X_count = 0;
    int C4Y_count = 0;

    vector<int> C2X_wire_numbers;
    vector<int> C2Y_wire_numbers;
    vector<int> C3X_wire_numbers;
    vector<int> C3Y_wire_numbers;
    vector<int> C4X_wire_numbers;
    vector<int> C4Y_wire_numbers;

    // counters for number of wc hits.
    if(selected_TOF2 > 6){
      for(int i=0; i<56; i++){
        if(C2X_L[i] > 0){
          C2X_count++;
          C2X_wire_numbers.push_back(i);
        }
      }

      for(int i=0; i<64; i++){
        if(C3X_L[i] > 0){
          C3X_count++;
          C3X_wire_numbers.push_back(i);
        }
      }

      for(int i=0; i<72; i++){
        if(C4X_L[i] > 0){
          C4X_count++;
          C4X_wire_numbers.push_back(i);
        }
      }

      for(int i=0; i<16; i++){
        if(C2Y_L[i] > 0){
          C2Y_count++;
          C2Y_wire_numbers.push_back(i);
        }
      }

      for(int i=0; i<16; i++){
        if(C3Y_L[i] > 0){
          C3Y_count++;
          C3Y_wire_numbers.push_back(i);
        }
      }

      for(int i=0; i<16; i++){
        if(C4Y_L[i] > 0){
          C4Y_count++;
          C4Y_wire_numbers.push_back(i);
        }
      }
    }
    else{
      for(int i=0; i<56; i++){
        if(C2X_R[i] > 0){
          C2X_count++;
          C2X_wire_numbers.push_back(i);
        }
      }

      for(int i=0; i<64; i++){
        if(C3X_R[i] > 0){
          C3X_count++;
          C3X_wire_numbers.push_back(i);
        }
      }

      for(int i=0; i<72; i++){
        if(C4X_R[i] > 0){
          C4X_count++;
          C4X_wire_numbers.push_back(i);
        }
      }

      for(int i=0; i<16; i++){
        if(C2Y_R[i] > 0){
          C2Y_count++;
          C2Y_wire_numbers.push_back(i);
        }
      }

      for(int i=0; i<16; i++){
        if(C3Y_R[i] > 0){
          C3Y_count++;
          C3Y_wire_numbers.push_back(i);
        }
      }

      for(int i=0; i<16; i++){
        if(C4Y_R[i] > 0){
          C4Y_count++;
          C4Y_wire_numbers.push_back(i);
        }
      }
    }

    TH2F *C2_hist = new TH2F("", "", 500, -280,280,500, -1,18);
    TH2F *C3_hist = new TH2F("", "", 500, -1,64,500, -1,18);
    TH2F *C4_hist = new TH2F("", "", 500, -1,72,500, -1,18);

    C2_hist->SetLabelSize(0.06,"X");
    C2_hist->SetLabelSize(0.06,"Y");
    C3_hist->SetLabelSize(0.06,"X");
    C3_hist->SetLabelSize(0.06,"Y");
    C4_hist->SetLabelSize(0.06,"X");
    C4_hist->SetLabelSize(0.06,"Y");

    C2_hist->SetNdivisions(18,"Y");
    C3_hist->SetNdivisions(18,"Y");
    C4_hist->SetNdivisions(18,"Y");

    TGaxis *C2_top = new TGaxis(-280,18,280,18,"pol1",510,"+U");
    TGaxis *C2_right = new TGaxis(280,-1,280,18,"pol1",510,"-U");
    TGaxis *C3_top = new TGaxis(-1,18,64,18,"pol1",510,"+U");
    TGaxis *C3_right = new TGaxis(64,-1,64,18,"pol1",510,"-U");
    TGaxis *C4_top = new TGaxis(-1,18,72,18,"pol1",510,"+U");
    TGaxis *C4_right = new TGaxis(72,-1,72,18,"pol1",510,"-U");

    TLine *C2X_line[C2X_count];
    TLine *C2Y_line[C2Y_count];
    TLine *C3X_line[C3X_count];
    TLine *C3Y_line[C3Y_count];
    TLine *C4X_line[C4X_count];
    TLine *C4Y_line[C4Y_count];

    const int MWPC_color_thr[10] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450};

    double ADC_palette_MWPC[10];

    for(int i = 0; i < 10; i++){
      ADC_palette_MWPC[i] = MWPC_color_thr[i];
    }

    char ADC_palette_string_MWPC[10][100];

    for(int j=0; j<10; j++) sprintf(ADC_palette_string_MWPC[j],"%2.1f",ADC_palette_MWPC[j]);
    sprintf(ADC_palette_string_MWPC[9],"%2.1f+",ADC_palette_MWPC[9]);

    double flag_size_palette_MWPC = 1;
    double text_size_MWPC = 0.04;

    double C2_palette_loc = -323;
    double C2_tex_palette_loc = -317;
    double C3_palette_loc = -6;
    double C3_tex_palette_loc = -5.3;
    double C4_palette_loc = -7;
    double C4_tex_palette_loc = -6.2;

    // C2 legend
    palette_MWPC_C2[0] = new TMarker(C2_palette_loc,15,21);     palette_MWPC_C2[0]->SetMarkerColor(kOrange+10);    palette_MWPC_C2[0]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C2[1] = new TMarker(C2_palette_loc,14,21);     palette_MWPC_C2[1]->SetMarkerColor(kOrange+7);     palette_MWPC_C2[1]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C2[2] = new TMarker(C2_palette_loc,13,21);     palette_MWPC_C2[2]->SetMarkerColor(kOrange+1);     palette_MWPC_C2[2]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C2[3] = new TMarker(C2_palette_loc,12,21);     palette_MWPC_C2[3]->SetMarkerColor(kOrange-4);     palette_MWPC_C2[3]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C2[4] = new TMarker(C2_palette_loc,11,21);     palette_MWPC_C2[4]->SetMarkerColor(kYellow-9);     palette_MWPC_C2[4]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C2[5] = new TMarker(C2_palette_loc,10,21);     palette_MWPC_C2[5]->SetMarkerColor(kYellow-7);     palette_MWPC_C2[5]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C2[6] = new TMarker(C2_palette_loc,9,21);      palette_MWPC_C2[6]->SetMarkerColor(kYellow-0);     palette_MWPC_C2[6]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C2[7] = new TMarker(C2_palette_loc,8,21);      palette_MWPC_C2[7]->SetMarkerColor(kSpring-4);     palette_MWPC_C2[7]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C2[8] = new TMarker(C2_palette_loc,7,21);      palette_MWPC_C2[8]->SetMarkerColor(kSpring-2);     palette_MWPC_C2[8]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C2[9] = new TMarker(C2_palette_loc,6,21);      palette_MWPC_C2[9]->SetMarkerColor(kGreen-0);      palette_MWPC_C2[9]->SetMarkerSize(flag_size_palette_MWPC);

    tex_palette_MWPC_C2[0] = new TLatex(C2_tex_palette_loc,14.5,ADC_palette_string_MWPC[0]);      tex_palette_MWPC_C2[0]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C2[1] = new TLatex(C2_tex_palette_loc,13.5,ADC_palette_string_MWPC[1]);      tex_palette_MWPC_C2[1]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C2[2] = new TLatex(C2_tex_palette_loc,12.5,ADC_palette_string_MWPC[2]);      tex_palette_MWPC_C2[2]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C2[3] = new TLatex(C2_tex_palette_loc,11.5,ADC_palette_string_MWPC[3]);      tex_palette_MWPC_C2[3]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C2[4] = new TLatex(C2_tex_palette_loc,10.5,ADC_palette_string_MWPC[4]);      tex_palette_MWPC_C2[4]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C2[5] = new TLatex(C2_tex_palette_loc,9.5,ADC_palette_string_MWPC[5]);       tex_palette_MWPC_C2[5]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C2[6] = new TLatex(C2_tex_palette_loc,8.5,ADC_palette_string_MWPC[6]);       tex_palette_MWPC_C2[6]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C2[7] = new TLatex(C2_tex_palette_loc,7.5,ADC_palette_string_MWPC[7]);       tex_palette_MWPC_C2[7]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C2[8] = new TLatex(C2_tex_palette_loc,6.5,ADC_palette_string_MWPC[8]);       tex_palette_MWPC_C2[8]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C2[9] = new TLatex(C2_tex_palette_loc,5.5,ADC_palette_string_MWPC[9]);       tex_palette_MWPC_C2[9]->SetTextSize(text_size_MWPC);

    palette_MWPC_C3[0] = new TMarker(C3_palette_loc,15,21);     palette_MWPC_C3[0]->SetMarkerColor(kOrange+10);    palette_MWPC_C3[0]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C3[1] = new TMarker(C3_palette_loc,14,21);     palette_MWPC_C3[1]->SetMarkerColor(kOrange+7);     palette_MWPC_C3[1]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C3[2] = new TMarker(C3_palette_loc,13,21);     palette_MWPC_C3[2]->SetMarkerColor(kOrange+1);     palette_MWPC_C3[2]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C3[3] = new TMarker(C3_palette_loc,12,21);     palette_MWPC_C3[3]->SetMarkerColor(kOrange-4);     palette_MWPC_C3[3]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C3[4] = new TMarker(C3_palette_loc,11,21);     palette_MWPC_C3[4]->SetMarkerColor(kYellow-9);     palette_MWPC_C3[4]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C3[5] = new TMarker(C3_palette_loc,10,21);     palette_MWPC_C3[5]->SetMarkerColor(kYellow-7);     palette_MWPC_C3[5]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C3[6] = new TMarker(C3_palette_loc,9,21);      palette_MWPC_C3[6]->SetMarkerColor(kYellow-0);     palette_MWPC_C3[6]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C3[7] = new TMarker(C3_palette_loc,8,21);      palette_MWPC_C3[7]->SetMarkerColor(kSpring-4);     palette_MWPC_C3[7]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C3[8] = new TMarker(C3_palette_loc,7,21);      palette_MWPC_C3[8]->SetMarkerColor(kSpring-2);     palette_MWPC_C3[8]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C3[9] = new TMarker(C3_palette_loc,6,21);      palette_MWPC_C3[9]->SetMarkerColor(kGreen-0);      palette_MWPC_C3[9]->SetMarkerSize(flag_size_palette_MWPC);

    TLatex *tex_palette_MWPC_scale = new TLatex(-4,6,"x 100");                 tex_palette_MWPC_scale->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C3[0] = new TLatex(C3_tex_palette_loc,14.5,ADC_palette_string_MWPC[0]);      tex_palette_MWPC_C3[0]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C3[1] = new TLatex(C3_tex_palette_loc,13.5,ADC_palette_string_MWPC[1]);      tex_palette_MWPC_C3[1]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C3[2] = new TLatex(C3_tex_palette_loc,12.5,ADC_palette_string_MWPC[2]);      tex_palette_MWPC_C3[2]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C3[3] = new TLatex(C3_tex_palette_loc,11.5,ADC_palette_string_MWPC[3]);      tex_palette_MWPC_C3[3]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C3[4] = new TLatex(C3_tex_palette_loc,10.5,ADC_palette_string_MWPC[4]);      tex_palette_MWPC_C3[4]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C3[5] = new TLatex(C3_tex_palette_loc,9.5,ADC_palette_string_MWPC[5]);       tex_palette_MWPC_C3[5]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C3[6] = new TLatex(C3_tex_palette_loc,8.5,ADC_palette_string_MWPC[6]);       tex_palette_MWPC_C3[6]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C3[7] = new TLatex(C3_tex_palette_loc,7.5,ADC_palette_string_MWPC[7]);       tex_palette_MWPC_C3[7]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C3[8] = new TLatex(C3_tex_palette_loc,6.5,ADC_palette_string_MWPC[8]);       tex_palette_MWPC_C3[8]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C3[9] = new TLatex(C3_tex_palette_loc,5.5,ADC_palette_string_MWPC[9]);       tex_palette_MWPC_C3[9]->SetTextSize(text_size_MWPC);

    palette_MWPC_C4[0] = new TMarker(C4_palette_loc,15,21);     palette_MWPC_C4[0]->SetMarkerColor(kOrange+10);    palette_MWPC_C4[0]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C4[1] = new TMarker(C4_palette_loc,14,21);     palette_MWPC_C4[1]->SetMarkerColor(kOrange+7);     palette_MWPC_C4[1]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C4[2] = new TMarker(C4_palette_loc,13,21);     palette_MWPC_C4[2]->SetMarkerColor(kOrange+1);     palette_MWPC_C4[2]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C4[3] = new TMarker(C4_palette_loc,12,21);     palette_MWPC_C4[3]->SetMarkerColor(kOrange-4);     palette_MWPC_C4[3]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C4[4] = new TMarker(C4_palette_loc,11,21);     palette_MWPC_C4[4]->SetMarkerColor(kYellow-9);     palette_MWPC_C4[4]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C4[5] = new TMarker(C4_palette_loc,10,21);     palette_MWPC_C4[5]->SetMarkerColor(kYellow-7);     palette_MWPC_C4[5]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C4[6] = new TMarker(C4_palette_loc,9,21);      palette_MWPC_C4[6]->SetMarkerColor(kYellow-0);     palette_MWPC_C4[6]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C4[7] = new TMarker(C4_palette_loc,8,21);      palette_MWPC_C4[7]->SetMarkerColor(kSpring-4);     palette_MWPC_C4[7]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C4[8] = new TMarker(C4_palette_loc,7,21);      palette_MWPC_C4[8]->SetMarkerColor(kSpring-2);     palette_MWPC_C4[8]->SetMarkerSize(flag_size_palette_MWPC);
    palette_MWPC_C4[9] = new TMarker(C4_palette_loc,6,21);      palette_MWPC_C4[9]->SetMarkerColor(kGreen-0);      palette_MWPC_C4[9]->SetMarkerSize(flag_size_palette_MWPC);

    tex_palette_MWPC_C4[0] = new TLatex(C4_tex_palette_loc,14.5,ADC_palette_string_MWPC[0]);      tex_palette_MWPC_C4[0]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C4[1] = new TLatex(C4_tex_palette_loc,13.5,ADC_palette_string_MWPC[1]);      tex_palette_MWPC_C4[1]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C4[2] = new TLatex(C4_tex_palette_loc,12.5,ADC_palette_string_MWPC[2]);      tex_palette_MWPC_C4[2]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C4[3] = new TLatex(C4_tex_palette_loc,11.5,ADC_palette_string_MWPC[3]);      tex_palette_MWPC_C4[3]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C4[4] = new TLatex(C4_tex_palette_loc,10.5,ADC_palette_string_MWPC[4]);      tex_palette_MWPC_C4[4]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C4[5] = new TLatex(C4_tex_palette_loc,9.5,ADC_palette_string_MWPC[5]);       tex_palette_MWPC_C4[5]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C4[6] = new TLatex(C4_tex_palette_loc,8.5,ADC_palette_string_MWPC[6]);       tex_palette_MWPC_C4[6]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C4[7] = new TLatex(C4_tex_palette_loc,7.5,ADC_palette_string_MWPC[7]);       tex_palette_MWPC_C4[7]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C4[8] = new TLatex(C4_tex_palette_loc,6.5,ADC_palette_string_MWPC[8]);       tex_palette_MWPC_C4[8]->SetTextSize(text_size_MWPC);
    tex_palette_MWPC_C4[9] = new TLatex(C4_tex_palette_loc,5.5,ADC_palette_string_MWPC[9]);       tex_palette_MWPC_C4[9]->SetTextSize(text_size_MWPC);

    TLatex *tex_Title_MWPC_C2;
    TLatex *tex_Title_MWPC_C3;
    TLatex *tex_Title_MWPC_C4;

    char c2_title[200];
    char c3_title[200];
    char c4_title[200];

    if(selected_TOF2 > 6){
      sprintf(c2_title,"C2 Left   |   Run Number %d   |    Event %d    |    C2X Thr = %5.2f    |    C2Y Thr = %5.2f   ", Run_Number, ievt, MWPC_Thr_Offset_C2X, MWPC_Thr_Offset_C2Y);
      sprintf(c3_title,"C3 Left   |   Run Number %d   |    Event %d    |    C3X Thr = %5.2f    |    C3Y Thr = %5.2f   ", Run_Number, ievt, MWPC_Thr_Offset_C3X, MWPC_Thr_Offset_C3Y);
      sprintf(c4_title,"C4 Left   |   Run Number %d   |    Event %d    |    C4X Thr = %5.2f    |    C4Y Thr = %5.2f   ", Run_Number, ievt, MWPC_Thr_Offset_C4X, MWPC_Thr_Offset_C4Y);
    }
    else{
      sprintf(c2_title,"C2 Right   |   Run Number %d   |    Event %d    |    C2X Thr = %5.2f    |    C2Y Thr = %5.2f   ", Run_Number, ievt, MWPC_Thr_Offset_C2X, MWPC_Thr_Offset_C2Y);
      sprintf(c3_title,"C3 Right   |   Run Number %d   |    Event %d    |    C3X Thr = %5.2f    |    C3Y Thr = %5.2f   ", Run_Number, ievt, MWPC_Thr_Offset_C3X, MWPC_Thr_Offset_C3Y);
      sprintf(c4_title,"C4 Right   |   Run Number %d   |    Event %d    |    C4X Thr = %5.2f    |    C4Y Thr = %5.2f   ", Run_Number, ievt, MWPC_Thr_Offset_C4X, MWPC_Thr_Offset_C4Y);
    }

    tex_Title_MWPC_C2 = new TLatex(-280,18.5,c2_title);
    tex_Title_MWPC_C3 = new TLatex(0,18.5,c3_title);
    tex_Title_MWPC_C4 = new TLatex(0,18.5,c4_title);

    tex_Title_MWPC_C2->SetTextSize(0.07645875);
    tex_Title_MWPC_C2->SetLineWidth(2);
    tex_Title_MWPC_C2->Draw();

    tex_Title_MWPC_C3->SetTextSize(0.07645875);
    tex_Title_MWPC_C3->SetLineWidth(2);
    tex_Title_MWPC_C3->Draw();

    tex_Title_MWPC_C4->SetTextSize(0.07645875);
    tex_Title_MWPC_C4->SetLineWidth(2);
    tex_Title_MWPC_C4->Draw();

    for(int i = 0; i < C2X_count; i++){
      C2X_line[i] = new TLine(C2_ZLoc[*(C2X_wire_numbers.begin() + i)], -1, C2_ZLoc[*(C2X_wire_numbers.begin() + i)], 18);

      if(selected_TOF2 > 6){
        if(C2X_L[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C2X_L[*(C2X_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C2X_line[i]->SetLineColor(kOrange+10);
        if(C2X_L[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C2X_L[*(C2X_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C2X_line[i]->SetLineColor(kOrange+7);
        if(C2X_L[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C2X_L[*(C2X_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C2X_line[i]->SetLineColor(kOrange+1);
        if(C2X_L[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C2X_L[*(C2X_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C2X_line[i]->SetLineColor(kOrange-4);
        if(C2X_L[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C2X_L[*(C2X_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C2X_line[i]->SetLineColor(kOrange-9);
        if(C2X_L[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C2X_L[*(C2X_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C2X_line[i]->SetLineColor(kOrange-7);
        if(C2X_L[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C2X_L[*(C2X_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C2X_line[i]->SetLineColor(kOrange-0);
        if(C2X_L[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C2X_L[*(C2X_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C2X_line[i]->SetLineColor(kOrange-4);
        if(C2X_L[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C2X_L[*(C2X_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C2X_line[i]->SetLineColor(kOrange-2);
        if(C2X_L[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C2X_line[i]->SetLineColor(kGreen-0);
      }
      else{
        if(C2X_R[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C2X_R[*(C2X_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C2X_line[i]->SetLineColor(kOrange+10);
        if(C2X_R[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C2X_R[*(C2X_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C2X_line[i]->SetLineColor(kOrange+7);
        if(C2X_R[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C2X_R[*(C2X_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C2X_line[i]->SetLineColor(kOrange+1);
        if(C2X_R[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C2X_R[*(C2X_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C2X_line[i]->SetLineColor(kOrange-4);
        if(C2X_R[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C2X_R[*(C2X_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C2X_line[i]->SetLineColor(kOrange-9);
        if(C2X_R[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C2X_R[*(C2X_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C2X_line[i]->SetLineColor(kOrange-7);
        if(C2X_R[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C2X_R[*(C2X_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C2X_line[i]->SetLineColor(kOrange-0);
        if(C2X_R[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C2X_R[*(C2X_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C2X_line[i]->SetLineColor(kOrange-4);
        if(C2X_R[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C2X_R[*(C2X_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C2X_line[i]->SetLineColor(kOrange-2);
        if(C2X_R[*(C2X_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C2X_line[i]->SetLineColor(kGreen-0);
      }
      C2X_line[i]->SetLineWidth(8);
    }

    for(int i = 0; i < C2Y_count; i++){
      C2Y_line[i] = new TLine(-280,*(C2Y_wire_numbers.begin() + i)+1, 280, *(C2Y_wire_numbers.begin() + i)+1);

      if(selected_TOF2 > 6){
        if(C2Y_L[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C2Y_L[*(C2Y_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C2Y_line[i]->SetLineColor(kOrange+10);
        if(C2Y_L[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C2Y_L[*(C2Y_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C2Y_line[i]->SetLineColor(kOrange+7);
        if(C2Y_L[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C2Y_L[*(C2Y_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C2Y_line[i]->SetLineColor(kOrange+1);
        if(C2Y_L[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C2Y_L[*(C2Y_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C2Y_line[i]->SetLineColor(kOrange-4);
        if(C2Y_L[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C2Y_L[*(C2Y_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C2Y_line[i]->SetLineColor(kOrange-9);
        if(C2Y_L[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C2Y_L[*(C2Y_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C2Y_line[i]->SetLineColor(kOrange-7);
        if(C2Y_L[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C2Y_L[*(C2Y_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C2Y_line[i]->SetLineColor(kOrange-0);
        if(C2Y_L[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C2Y_L[*(C2Y_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C2Y_line[i]->SetLineColor(kOrange-4);
        if(C2Y_L[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C2Y_L[*(C2Y_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C2Y_line[i]->SetLineColor(kOrange-2);
        if(C2Y_L[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C2Y_line[i]->SetLineColor(kGreen-0);
      }
      else{
        if(C2Y_R[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C2Y_R[*(C2Y_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C2Y_line[i]->SetLineColor(kOrange+10);
        if(C2Y_R[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C2Y_R[*(C2Y_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C2Y_line[i]->SetLineColor(kOrange+7);
        if(C2Y_R[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C2Y_R[*(C2Y_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C2Y_line[i]->SetLineColor(kOrange+1);
        if(C2Y_R[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C2Y_R[*(C2Y_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C2Y_line[i]->SetLineColor(kOrange-4);
        if(C2Y_R[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C2Y_R[*(C2Y_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C2Y_line[i]->SetLineColor(kOrange-9);
        if(C2Y_R[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C2Y_R[*(C2Y_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C2Y_line[i]->SetLineColor(kOrange-7);
        if(C2Y_R[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C2Y_R[*(C2Y_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C2Y_line[i]->SetLineColor(kOrange-0);
        if(C2Y_R[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C2Y_R[*(C2Y_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C2Y_line[i]->SetLineColor(kOrange-4);
        if(C2Y_R[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C2Y_R[*(C2Y_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C2Y_line[i]->SetLineColor(kOrange-2);
        if(C2Y_R[*(C2Y_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C2Y_line[i]->SetLineColor(kGreen-0);
      }
      C2Y_line[i]->SetLineWidth(8);
    }

    for(int i = 0; i < C3X_count; i++){
      C3X_line[i] = new TLine(*(C3X_wire_numbers.begin() + i)+1, -1, *(C3X_wire_numbers.begin() + i)+1, 18);

      if(selected_TOF2 > 6){
        if(C3X_L[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C3X_L[*(C3X_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C3X_line[i]->SetLineColor(kOrange+10);
        if(C3X_L[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C3X_L[*(C3X_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C3X_line[i]->SetLineColor(kOrange+7);
        if(C3X_L[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C3X_L[*(C3X_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C3X_line[i]->SetLineColor(kOrange+1);
        if(C3X_L[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C3X_L[*(C3X_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C3X_line[i]->SetLineColor(kOrange-4);
        if(C3X_L[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C3X_L[*(C3X_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C3X_line[i]->SetLineColor(kOrange-9);
        if(C3X_L[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C3X_L[*(C3X_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C3X_line[i]->SetLineColor(kOrange-7);
        if(C3X_L[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C3X_L[*(C3X_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C3X_line[i]->SetLineColor(kOrange-0);
        if(C3X_L[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C3X_L[*(C3X_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C3X_line[i]->SetLineColor(kOrange-4);
        if(C3X_L[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C3X_L[*(C3X_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C3X_line[i]->SetLineColor(kOrange-2);
        if(C3X_L[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C3X_line[i]->SetLineColor(kGreen-0);
      }
      else{
        if(C3X_R[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C3X_R[*(C3X_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C3X_line[i]->SetLineColor(kOrange+10);
        if(C3X_R[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C3X_R[*(C3X_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C3X_line[i]->SetLineColor(kOrange+7);
        if(C3X_R[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C3X_R[*(C3X_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C3X_line[i]->SetLineColor(kOrange+1);
        if(C3X_R[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C3X_R[*(C3X_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C3X_line[i]->SetLineColor(kOrange-4);
        if(C3X_R[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C3X_R[*(C3X_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C3X_line[i]->SetLineColor(kOrange-9);
        if(C3X_R[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C3X_R[*(C3X_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C3X_line[i]->SetLineColor(kOrange-7);
        if(C3X_R[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C3X_R[*(C3X_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C3X_line[i]->SetLineColor(kOrange-0);
        if(C3X_R[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C3X_R[*(C3X_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C3X_line[i]->SetLineColor(kOrange-4);
        if(C3X_R[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C3X_R[*(C3X_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C3X_line[i]->SetLineColor(kOrange-2);
        if(C3X_R[*(C3X_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C3X_line[i]->SetLineColor(kGreen-0);
      }
      C3X_line[i]->SetLineWidth(8);
    }

    for(int i = 0; i < C3Y_count; i++){
      C3Y_line[i] = new TLine(-1,*(C3Y_wire_numbers.begin() + i)+1, 64, *(C3Y_wire_numbers.begin() + i)+1);

      if(selected_TOF2 > 6){
        if(C3Y_L[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C3Y_L[*(C3Y_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C3Y_line[i]->SetLineColor(kOrange+10);
        if(C3Y_L[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C3Y_L[*(C3Y_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C3Y_line[i]->SetLineColor(kOrange+7);
        if(C3Y_L[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C3Y_L[*(C3Y_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C3Y_line[i]->SetLineColor(kOrange+1);
        if(C3Y_L[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C3Y_L[*(C3Y_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C3Y_line[i]->SetLineColor(kOrange-4);
        if(C3Y_L[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C3Y_L[*(C3Y_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C3Y_line[i]->SetLineColor(kOrange-9);
        if(C3Y_L[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C3Y_L[*(C3Y_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C3Y_line[i]->SetLineColor(kOrange-7);
        if(C3Y_L[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C3Y_L[*(C3Y_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C3Y_line[i]->SetLineColor(kOrange-0);
        if(C3Y_L[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C3Y_L[*(C3Y_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C3Y_line[i]->SetLineColor(kOrange-4);
        if(C3Y_L[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C3Y_L[*(C3Y_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C3Y_line[i]->SetLineColor(kOrange-2);
        if(C3Y_L[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C3Y_line[i]->SetLineColor(kGreen-0);
      }
      else{
        if(C3Y_R[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C3Y_R[*(C3Y_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C3Y_line[i]->SetLineColor(kOrange+10);
        if(C3Y_R[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C3Y_R[*(C3Y_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C3Y_line[i]->SetLineColor(kOrange+7);
        if(C3Y_R[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C3Y_R[*(C3Y_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C3Y_line[i]->SetLineColor(kOrange+1);
        if(C3Y_R[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C3Y_R[*(C3Y_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C3Y_line[i]->SetLineColor(kOrange-4);
        if(C3Y_R[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C3Y_R[*(C3Y_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C3Y_line[i]->SetLineColor(kOrange-9);
        if(C3Y_R[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C3Y_R[*(C3Y_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C3Y_line[i]->SetLineColor(kOrange-7);
        if(C3Y_R[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C3Y_R[*(C3Y_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C3Y_line[i]->SetLineColor(kOrange-0);
        if(C3Y_R[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C3Y_R[*(C3Y_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C3Y_line[i]->SetLineColor(kOrange-4);
        if(C3Y_R[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C3Y_R[*(C3Y_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C3Y_line[i]->SetLineColor(kOrange-2);
        if(C3Y_R[*(C3Y_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C3Y_line[i]->SetLineColor(kGreen-0);
      }
      C3Y_line[i]->SetLineWidth(8);
    }

    for(int i = 0; i < C4X_count; i++){
      C4X_line[i] = new TLine(*(C4X_wire_numbers.begin() + i)+1, -1, *(C4X_wire_numbers.begin() + i)+1, 18);

      if(selected_TOF2 > 6){
        if(C4X_L[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C4X_L[*(C4X_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C4X_line[i]->SetLineColor(kOrange+10);
        if(C4X_L[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C4X_L[*(C4X_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C4X_line[i]->SetLineColor(kOrange+7);
        if(C4X_L[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C4X_L[*(C4X_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C4X_line[i]->SetLineColor(kOrange+1);
        if(C4X_L[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C4X_L[*(C4X_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C4X_line[i]->SetLineColor(kOrange-4);
        if(C4X_L[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C4X_L[*(C4X_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C4X_line[i]->SetLineColor(kOrange-9);
        if(C4X_L[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C4X_L[*(C4X_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C4X_line[i]->SetLineColor(kOrange-7);
        if(C4X_L[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C4X_L[*(C4X_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C4X_line[i]->SetLineColor(kOrange-0);
        if(C4X_L[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C4X_L[*(C4X_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C4X_line[i]->SetLineColor(kOrange-4);
        if(C4X_L[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C4X_L[*(C4X_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C4X_line[i]->SetLineColor(kOrange-2);
        if(C4X_L[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C4X_line[i]->SetLineColor(kGreen-0);
      }
      else{
        if(C4X_R[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C4X_R[*(C4X_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C4X_line[i]->SetLineColor(kOrange+10);
        if(C4X_R[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C4X_R[*(C4X_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C4X_line[i]->SetLineColor(kOrange+7);
        if(C4X_R[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C4X_R[*(C4X_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C4X_line[i]->SetLineColor(kOrange+1);
        if(C4X_R[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C4X_R[*(C4X_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C4X_line[i]->SetLineColor(kOrange-4);
        if(C4X_R[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C4X_R[*(C4X_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C4X_line[i]->SetLineColor(kOrange-9);
        if(C4X_R[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C4X_R[*(C4X_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C4X_line[i]->SetLineColor(kOrange-7);
        if(C4X_R[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C4X_R[*(C4X_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C4X_line[i]->SetLineColor(kOrange-0);
        if(C4X_R[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C4X_R[*(C4X_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C4X_line[i]->SetLineColor(kOrange-4);
        if(C4X_R[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C4X_R[*(C4X_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C4X_line[i]->SetLineColor(kOrange-2);
        if(C4X_R[*(C4X_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C4X_line[i]->SetLineColor(kGreen-0);
      }
      C4X_line[i]->SetLineWidth(8);
    }

    for(int i = 0; i < C4Y_count; i++){
      C4Y_line[i] = new TLine(-1,*(C4Y_wire_numbers.begin() + i)+1, 72, *(C4Y_wire_numbers.begin() + i)+1);

      if(selected_TOF2 > 6){
        if(C4Y_L[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C4Y_L[*(C4Y_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C4Y_line[i]->SetLineColor(kOrange+10);
        if(C4Y_L[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C4Y_L[*(C4Y_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C4Y_line[i]->SetLineColor(kOrange+7);
        if(C4Y_L[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C4Y_L[*(C4Y_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C4Y_line[i]->SetLineColor(kOrange+1);
        if(C4Y_L[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C4Y_L[*(C4Y_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C4Y_line[i]->SetLineColor(kOrange-4);
        if(C4Y_L[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C4Y_L[*(C4Y_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C4Y_line[i]->SetLineColor(kOrange-9);
        if(C4Y_L[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C4Y_L[*(C4Y_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C4Y_line[i]->SetLineColor(kOrange-7);
        if(C4Y_L[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C4Y_L[*(C4Y_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C4Y_line[i]->SetLineColor(kOrange-0);
        if(C4Y_L[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C4Y_L[*(C4Y_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C4Y_line[i]->SetLineColor(kOrange-4);
        if(C4Y_L[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C4Y_L[*(C4Y_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C4Y_line[i]->SetLineColor(kOrange-2);
        if(C4Y_L[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C4Y_line[i]->SetLineColor(kGreen-0);
      }
      else{
        if(C4Y_R[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C4Y_R[*(C4Y_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C4Y_line[i]->SetLineColor(kOrange+10);
        if(C4Y_R[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C4Y_R[*(C4Y_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C4Y_line[i]->SetLineColor(kOrange+7);
        if(C4Y_R[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C4Y_R[*(C4Y_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C4Y_line[i]->SetLineColor(kOrange+1);
        if(C4Y_R[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C4Y_R[*(C4Y_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C4Y_line[i]->SetLineColor(kOrange-4);
        if(C4Y_R[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C4Y_R[*(C4Y_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C4Y_line[i]->SetLineColor(kOrange-9);
        if(C4Y_R[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C4Y_R[*(C4Y_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C4Y_line[i]->SetLineColor(kOrange-7);
        if(C4Y_R[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C4Y_R[*(C4Y_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C4Y_line[i]->SetLineColor(kOrange-0);
        if(C4Y_R[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C4Y_R[*(C4Y_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C4Y_line[i]->SetLineColor(kOrange-4);
        if(C4Y_R[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C4Y_R[*(C4Y_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C4Y_line[i]->SetLineColor(kOrange-2);
        if(C4Y_R[*(C4Y_wire_numbers.begin() + i)] > MWPC_color_thr[6]                                                               ) C4Y_line[i]->SetLineColor(kGreen-0);
      }
      C4Y_line[i]->SetLineWidth(8);
    }


    gStyle->SetOptStat(0); // Remove legends

    Draw_MWPC(cMWPC,
      C2_hist, C3_hist, C4_hist,
      tex_Title_MWPC_C2,
      tex_Title_MWPC_C3,
      tex_Title_MWPC_C4,
      tex_palette_MWPC_C2,
      tex_palette_MWPC_C3,
      tex_palette_MWPC_C4,
      palette_MWPC_C2,
      palette_MWPC_C3,
      palette_MWPC_C4,
      C2_top, C2_right,
      C3_top, C3_right,
      C4_top, C4_right,
      C2X_count, C2Y_count,
      C3X_count, C3Y_count,
      C4X_count, C4Y_count,
      C2X_line,
      C2Y_line,
      C3X_line,
      C3Y_line,
      C4X_line,
      C4Y_line);
    }

  // B0 Counter
  for(int i=0; i<16; i++){
    if(tdc_vt48[6][i]>TDC_B0_min && tdc_vt48[6][i]<TDC_B0_max) vec_tdc_b0_6.push_back(tdc_vt48[6][i]);
    else vec_tdc_b0_6.push_back(-1);

    if(tdc_vt48[7][i]>TDC_B0_min && tdc_vt48[7][i]<TDC_B0_max) vec_tdc_b0_7.push_back(tdc_vt48[7][i]);
    else vec_tdc_b0_7.push_back(-1);
  }

  sort(vec_tdc_b0_6.begin(), vec_tdc_b0_6.end());
  sort(vec_tdc_b0_7.begin(), vec_tdc_b0_7.end());


  // ### Display Ck and Cpi infomation

  cout << "   " << endl;
  cout << "   " << endl;

  cout << "//////  C_k / C_pi  //////" << endl;
  cout << "Counter       TDC[0]           TDC[1]           TDC[2]           TDC[3]           TDC[4]           TDC[5]" << endl;

  for(int ic = 0; ic < 14; ic++){
    printf("%3d        %4d / %4d      %4d / %4d      %4d / %4d      %4d / %4d      %4d / %4d      %4d / %4d\n",ic,
    tdc_Ck[ic][0], tdc_Cpi[ic][0], tdc_Ck[ic][1], tdc_Cpi[ic][1], tdc_Ck[ic][2], tdc_Cpi[ic][2],
    tdc_Ck[ic][3], tdc_Cpi[ic][3], tdc_Ck[ic][4], tdc_Cpi[ic][4], tdc_Ck[ic][5], tdc_Cpi[ic][5]);

    for(int jc=0; jc<6; jc++){
      if(tdc_Ck[ic][jc]>=TDC_Ck_min && tdc_Ck[ic][jc]<=TDC_Ck_max) TDC_ck_selected[ic]=tdc_Ck[ic][jc];
      if(tdc_Cpi[ic][jc]>=TDC_Cpi_min && tdc_Cpi[ic][jc]<=TDC_Cpi_max) TDC_cpi_selected[ic]=tdc_Cpi[ic][jc];
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

  cout << endl;
  cout << "N(C_k) = " << vec_Ck.size() << "   " << "TDC Ck Avg. = " << TDC_ck_avg << "   " << "Sigma = " << TDC_ck_sigma << endl;
  cout << "N(C_pi) =  " << vec_Cpi.size() << "   " << "TDC Cpi Avg. = " << TDC_cpi_avg << "   " << "Sigma = " << TDC_cpi_sigma << endl;
  cout << "TDC Mean Ck = " << TDC_ck_avg2 << endl;
  cout << "TDC Mean Cpi = " << TDC_cpi_avg2 << endl;
  cout << endl;
  cout << "TDC Mean Ck - TDC Average = " << TDC_ck_avg2 - TDC_average << endl;
  cout << endl;

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

  if(Good_Event && (enable_cout==0 || enable_cout==9)){
    cout << endl;
    cout << "////// TDC TARGET - TDC Ck (corrected)  //////" << endl;
    cout << "Fiber  HG-Ped   LG-Ped     T[0]-Ck    T[1]-Ck    T[2]-Ck    T[3]-Ck    T[4]-Ck    T[5]-Ck" << endl;
    for(unsigned int jj=0; jj<vec_TARGET_bar_selected.size(); jj++){
        printf("%3d     %4d    %4d        %s      %s      %s      %s      %s      %s\n",
        vec_TARGET_bar_selected[jj], ADC_High_TARGET[vec_TARGET_bar_selected[jj]]+TARGET_ADC_Thr_HG_Offset, ADC_Low_TARGET[vec_TARGET_bar_selected[jj]] + TARGET_ADC_Thr_LG_Offset,
        TDC_LE_TARGET_sel_string[jj][0], TDC_LE_TARGET_sel_string[jj][1], TDC_LE_TARGET_sel_string[jj][2], TDC_LE_TARGET_sel_string[jj][3], TDC_LE_TARGET_sel_string[jj][4], TDC_LE_TARGET_sel_string[jj][5]);
    }
    cout << endl;
  }


  double X_BAR = 999.;
  double Y_BAR = 999.;
  int i_kaon_bar;
  TCanvas *c2;
  c2 = new TCanvas("Event_Display.C  --  TARGET & SFT","Event_Display.C  --  TARGET & SFT",0,200,1050,700);
  if (graphics == 1) {


    c2->Divide(3,2);
    c2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

    c2->cd(1);
    Subplot_1(ch_ADC_cut_TARGET, tex_Legend_TARGET, tex_event_TARGET,
      ADC_High_TARGET, marker_ADC_TARGET,
      ADC_Low_TARGET, TDC_min_TARGET,
      Switch, has_TDC_hit, palette_TARGET, tex_palette_TARGET,
      tex_palette_TARGET_scale, max_index,
      max_index2, max_index3, max_index4,
      marker_TDC_TARGET);

    c2->cd(2);
    Subplot_2(ch_ADC_cut_TARGET, tex_Legend_TARGET, tex_event_TARGET,
      ADC_High_TARGET, marker_ADC_TARGET,
      ADC_Low_TARGET, TDC_min_TARGET,
      Switch, has_TDC_hit, palette_TARGET, tex_palette_TARGET,
      tex_palette_TARGET_scale, max_index,
      max_index2, max_index3, max_index4,
      marker_TDC_TARGET);

    c2->cd(3);
    i_kaon_bar = Subplot_3(ch_ADC_cut_TARGET,
      tex_Legend_TARGET, tex_event_TARGET,
      ADC_High_TARGET, marker_ADCL_TARGET,
      ADC_Low_TARGET,
      has_TDC_hit_Kstop,
      palette_TARGET, tex_palette_TARGET,
      tex_palette_TARGET_scale,
      good_k_stop_bars, X_BAR, Y_BAR);
    // double X_BAR = subplot_3_return_value.X_BAR;
    // double Y_BAR = subplot_3_return_value.Y_BAR;
    // int i_kaon_bar = subplot_3_return_value.i_kaon_bar;
    // TLatex *tex_Label_Centroid_x = subplot_3_return_value.tex_Label_Centroid_x;
    // TLatex *tex_Label_Centroid_y = subplot_3_return_value.tex_Label_Centroid_y;
    cout << "hello, display is " << display <<endl;
    cout << "X, YBar" << X_BAR << ", " << Y_BAR << endl;


    SetupGaps();
  }
  else {
    if (graphics == 0) {
      int i_kaon_bar = 999;
      i_kaon_bar = kaon_fiber(X_BAR,Y_BAR);
    }
  }





  Double_t par0_ADC = 0;  float par0_TDC = 0.;
  Double_t par1_ADC = 0;  float par1_TDC = 0.;

  char angle_string_ADC[100];

  TLatex *tex_angle_ADC;

  if (graphics == 1) {
    TCanvas *c_SFT;
    c_SFT = new TCanvas("Event_Display.C  --  SFT","Event_Display.C  -- SFT",0,200,700,350);
    c_SFT->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
    c_SFT->Divide(2,1);

    c_SFT->cd(1);
    SFT_subplot_1(Run_Number, ievt,
      Event_flag,
      hline_DS, vline_DS,
      hline_US, vline_US,
      ADC_High_SFT_corr, has_TDC_SFT_hit,
      ADC_High_SFT, ADC_Low_TARGET,
      ADC_Low_SFT,
      marker_DS, marker_US,
      par_temp, Switch);


    c_SFT->cd(2);
    SFT_subplot_2(h_ADC_L1_DS, h_ADC_L2_DS,
    h_ADC_L3_DS, h_ADC_L4_DS,
    h_ADC_L1_US, h_ADC_L2_US,
    h_ADC_L3_US, h_ADC_L4_US, ADC_High_corr_max);
  }




  //// Fill K Stop bars
  //for(int i = 0; i<256; i++){
  //  if(k_stop_bar[i]){
  //    cout << "TEST KBARS :  " << i << endl;
  //    vec_xx_kaon.push_back(Xloc[i]);
  //   vec_yy_kaon.push_back(Yloc[i]);
  //  }
  //}

  cout << "  " << endl;

  const int gap_to_fit_rotate = gap_to_fit;

/*  if((gap_to_fit==12 || gap_to_fit==6 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
    int ADC_High_TARGET_temp[256];
    int ADC_Low_TARGET_temp[256];
    bool has_TDC_hit_temp[256];
    bool has_TDC_hit_Kstop_temp[256];
    int HG_TARGET_ADC_Thr_temp[256];
    int LG_TARGET_ADC_Thr_temp[256];
    int TARGET_High_has_neighbours_temp[256];
    bool k_stop_bar_temp[256];

    if(gap_to_fit==11) gap_to_fit = 2;
    else if(gap_to_fit==12) gap_to_fit = 3;
    else if(gap_to_fit==1) gap_to_fit = 4;
    else if(gap_to_fit==5) gap_to_fit = 8;
    else if(gap_to_fit==6) gap_to_fit = 9;
    else if(gap_to_fit==7) gap_to_fit = 10;

    vec_xx_lepton.clear();  vec_yy_lepton.clear();
    vec_ex_lepton.clear();  vec_ey_lepton.clear();
    vec_lepton_size.clear();  vec_lepton_rotate_size.clear();

    for(int i = 0; i<256; i++){
      ADC_High_TARGET_temp[i] = ADC_High_TARGET[i];
      ADC_Low_TARGET_temp[i] = ADC_Low_TARGET[i];
      HG_TARGET_ADC_Thr_temp[i] = HG_TARGET_ADC_Thr[i];
      LG_TARGET_ADC_Thr_temp[i] = LG_TARGET_ADC_Thr[i];
      has_TDC_hit_temp[i] = has_TDC_hit[i];
      has_TDC_hit_Kstop_temp[i] = has_TDC_hit_Kstop[i];
      TARGET_High_has_neighbours_temp[i] = TARGET_High_has_neighbours[i];
      k_stop_bar_temp[i] = k_stop_bar[i];
    }

    for(int i=0; i<256; i++){
      ADC_High_TARGET[i] = ADC_High_TARGET_temp[TARGET_Rotated_index[i]];
      ADC_Low_TARGET[i] = ADC_Low_TARGET_temp[TARGET_Rotated_index[i]];
      HG_TARGET_ADC_Thr[i] = HG_TARGET_ADC_Thr_temp[TARGET_Rotated_index[i]];
      LG_TARGET_ADC_Thr[i] = LG_TARGET_ADC_Thr_temp[TARGET_Rotated_index[i]];
      has_TDC_hit[i] = has_TDC_hit_temp[TARGET_Rotated_index[i]];
      has_TDC_hit_Kstop[i] = has_TDC_hit_Kstop_temp[TARGET_Rotated_index[i]];
      adc_low_target[i] = adc_low_target[TARGET_Rotated_index[i]];
      TARGET_High_has_neighbours[i] = TARGET_High_has_neighbours_temp[TARGET_Rotated_index[i]];
      k_stop_bar[i] = k_stop_bar_temp[TARGET_Rotated_index[i]];
    }

    max_index = 0;
    max_index2 = 0;
    max_index3 = 0;
    max_index4 = 0;

    max_ADC = -100000000;
    max_ADC2 = -100000000;
    max_ADC3 = -100000000;
    max_ADC4 = -100000000;

    for(Int_t q=0; q<256; q++){
      if(ADC_High_TARGET[q]>max_ADC){
        max_index = q;
        max_ADC = ADC_Low_TARGET[q];
      }
    }

    for(Int_t q=0; q<256; q++){
      if (q == max_index) continue;
      else {
        if(ADC_High_TARGET[q]>max_ADC2){
          max_index2 = q;
          max_ADC2 = ADC_Low_TARGET[q];
        }
      }
    }

    for(Int_t q=0; q<256; q++){
      if ((q == max_index) || (q == max_index2)) continue;
      else {
        if(ADC_High_TARGET[q]>max_ADC3){
          max_index3 = q;
          max_ADC3 = ADC_Low_TARGET[q];
        }
      }
    }

    for(Int_t q=0; q<256; q++){
      if ((q == max_index) || (q == max_index2) || (q == max_index3)) continue;
      else {
        if(ADC_High_TARGET[q]>max_ADC4) {
          max_index4 = q;
          max_ADC4 = ADC_Low_TARGET[q];
        }
      }
    }

    x_cent = Xloc[max_index];
    y_cent = Yloc[max_index];

    for(Int_t j=0; j<256; j++) hyp[j] = sqrt(pow(x_cent - Xloc[j],2) + pow(y_cent - Yloc[j],2));

    count = 0;
    hit_count = 0;

    for(Int_t i=0; i<256; i++){
      if(TARGET_High_has_neighbours[i]){
        if(ADC_High_TARGET[i]>Angle_ADC_cut && has_TDC_hit[i]){
          if(!k_stop_bar[i]){
            vec_xx_lepton.push_back(Xloc[i]);
            vec_yy_lepton.push_back(Yloc[i]);
            vec_lepton_size.push_back(Xloc[i]);
          }
          count++;
          if(IsIn(i,channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],   // Additional weight for fibers close to the edge if hit
                    channel[gap_to_fit-1][2], channel[gap_to_fit-1][3],
                    channel[gap_to_fit-1][4], channel[gap_to_fit-1][5],
                    channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){
            if(!k_stop_bar[i]){
              vec_xx_lepton.push_back(Xloc[i]);
              vec_xx_lepton.push_back(Xloc[i]);
              vec_yy_lepton.push_back(Yloc[i]);
              vec_yy_lepton.push_back(Yloc[i]);
              //vec_lepton_size.push_back(Xloc[i]);
            }
          }
        }

        if(ADC_High_TARGET[i]>Angle_ADC_cut){
          x_inc += Xloc[i];
          y_inc += Yloc[i];
          hit_count++;
        }

        if(ADC_High_TARGET[i]<0 && ADC_Low_TARGET[i]>0 && Switch==1){
          if(!k_stop_bar[i]){
            vec_xx_lepton.push_back(Xloc[i]);
            vec_yy_lepton.push_back(Yloc[i]);
            vec_lepton_size.push_back(Xloc[i]);
          }
          count++;
        }
      }
    }
  }
*/
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

      double SFT_x_intercept;
      double SFT_y_intercept;

      if (SFTxhyp1 < SFTxhyp2){
        SFT_x_intercept = x_circle_int1;
        SFT_y_intercept = y_circle_int1;
      }
      else {
        SFT_x_intercept = x_circle_int2;
        SFT_y_intercept = y_circle_int2;
      }

      cout << "" << endl;
      cout << "SFT Circle Intercept" << endl;
      cout << "X coordinate: " << SFT_x_intercept << " -- Y coordinate: " << SFT_y_intercept << endl;

      /// Determine angle phi between centre of TARGET and SFT ring intercept
      double SFT_phi;

      if (SFT_x_intercept >0) SFT_phi = (atan2(SFT_x_intercept,SFT_y_intercept))*57.2957795;
      else SFT_phi = 360 + (atan2(SFT_x_intercept,SFT_y_intercept))*57.2957795;

      sprintf(angle_string_ADC,"Phi = %3.4f deg.",SFT_phi);
      tex_angle_ADC=new TLatex(0,-9,angle_string_ADC);
      tex_angle_ADC->SetTextSize(0.05);
      tex_angle_ADC->Draw();

      cout << "" << endl;
      cout << "Angle between SFT Layer Ring Intercept and Centre of TARGET: " << SFT_phi << " deg." << endl;
      cout << "" << endl;

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
  else cout << "Histo Fit 4 Is Empty" << endl;

  cout << "//////  TOF1 SCORING  //////" << endl;
  for (int k=0; k<12; k++) cout << "TOF1 score " << k+1 << ": " << gap_counter[k] << endl;

  cout << "" << endl;
  cout << "TOF1 SELECTED FOR FITTING:  TOF1 " << gap_to_fit_rotate << endl;

  cout << "" << endl;
  cout << "" << endl;
  for (Int_t i = 0; i<12; i++){
    if(Good_Event){
      if (ADC_TOF1[i] > 0 || ADC_TOF1[i+12] > 0){
        if ((TDC_TOF1U[i]>TOF1U_TDC_min[i] && TDC_TOF1U[i]<TOF1U_TDC_max[i]) || (TDC_TOF1D[i]>TOF1D_TDC_min[i] && TDC_TOF1D[i]<TOF1D_TDC_max[i])){
          printf("ADC TOF1 Up-%2d: %5d -- Down-%2d: %5d  |  TDC TOF1 Up-%2d: %5d -- Down-%2d: %5d\n", i+1, ADC_TOF1[i], i+1, ADC_TOF1[i+12], i+1, TDC_TOF1U[i], i+1, TDC_TOF1D[i]);
        }
        else {
          printf("ADC TOF1 Up-%2d: %5d -- Down-%2d: %5d  |\n", i+1, ADC_TOF1[i], i+1, ADC_TOF1[i+12]);
        }
      }
      else {
        if ((TDC_TOF1U[i]>TOF1U_TDC_min[i] && TDC_TOF1U[i]<TOF1U_TDC_max[i]) || (TDC_TOF1D[i]>TOF1D_TDC_min[i] && TDC_TOF1D[i]<TOF1D_TDC_max[i])) {
          printf("                                         |  TDC TOF1 Up-%2d: %5d -- Down-%2d: %5d\n", i+1, TDC_TOF1U[i], i+1, TDC_TOF1D[i]);
        }
      }
    }
  }

  cout << "" << endl;

  for (Int_t i = 0; i<12; i++) {
    if (ADC_TOF2AO[i] > 0 || ADC_TOF2AI[i] > 0){
      if ((TDC_TOF2AO[i]>TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<TOF2AO_TDC_max[i]) || (TDC_TOF2AI[i]>TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<TOF2AI_TDC_max[i])) {
        printf("ADC TOF2 OutA-%2d: %5d -- InA-%2d: %5d  |  TDC TOF2 OutA-%2d: %5d -- InA-%2d: %5d\n", i+1, ADC_TOF2AO[i], i+1, ADC_TOF2AI[i], i+1, TDC_TOF2AO[i], i+1, TDC_TOF2AI[i]);
      }
      else {
        printf("ADC TOF2 OutA-%2d: %5d -- InA-%2d: %5d  |\n", i+1, ADC_TOF2AO[i], i+1, ADC_TOF2AI[i]);
      }
    }
    else {
      if ((TDC_TOF2AO[i]>TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<TOF2AO_TDC_max[i]) || (TDC_TOF2AI[i]>TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<TOF2AI_TDC_max[i])) {
        printf("                                        |  TDC TOF2 OutA-%2d: %5d -- InA-%2d: %5d\n", i+1, TDC_TOF2AO[i], i+1, TDC_TOF2AI[i]);
      }
    }
  }

  for (Int_t i = 0; i<12; i++){
    if (i == 6) {
      if (ADC_TOF2BO[i] > 0 || ADC_TOF2BO[i] > 0) {
        if ((TDC_TOF2BO[6]>TOF2BO_TDC_min[6] && TDC_TOF2BO[6]<TOF2BO_TDC_max[6]) || (TDC_TOF2BI[6]>TOF2BI_TDC_min[6] && TDC_TOF2BI[6]<TOF2BI_TDC_max[6])) {
          printf("ADC TOF2 OutB-%2d: %5d -- InB-%2d: %5d  |  TDC TOF2 OutB-%2d: %5d -- InB-%2d: %5d\n", i+1, ADC_TOF2BO[i], i+1, ADC_TOF2BI[i], i+1, TDC_TOF2BO[i], i+1, TDC_TOF2BI[i]);
        }
        else {
          cout << "ADC TOF2 OutB-" << i+1 << ": " << ADC_TOF2BO[i] << " -- InB-" << i+1 << ": " << ADC_TOF2BI[i] << endl;
          printf("ADC TOF2 OutB-%2d: %5d -- InB-%2d: %5d  |\n", i+1, ADC_TOF2BO[i], i+1, ADC_TOF2BI[i]);
        }
      }
      else {
        if ((TDC_TOF2BO[6]>TOF2BO_TDC_min[6] && TDC_TOF2BO[6]>TOF2BO_TDC_max[6]) || (TDC_TOF2BI[6]>TOF2BI_TDC_min[6] && TDC_TOF2BI[6]>TOF2BI_TDC_max[6])) {
          printf("                                        |  TDC TOF2 OutB-%2d: %5d -- InB-%2d: %5d\n", i+1, TDC_TOF2BO[i], i+1, TDC_TOF2BI[i]);
        }
      }
    }
    else {
      if (ADC_TOF2BO[i] > 0 || ADC_TOF2BI[i] > 0) {
        if ((TDC_TOF2BO[i]>TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<TOF2BO_TDC_max[i]) || (TDC_TOF2BI[i]>TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<TOF2BI_TDC_max[i])) {
          printf("ADC TOF2 OutB-%2d: %5d -- InB-%2d: %5d  |  TDC TOF2 OutB-%2d: %5d -- InB-%2d: %5d\n", i+1, ADC_TOF2BO[i], i+1, ADC_TOF2BI[i], i+1, TDC_TOF2BO[i], i+1, TDC_TOF2BI[i]);
        }
        else {
          printf("ADC TOF2 OutB-%2d: %5d -- InB-%2d: %5d  |\n", i+1, ADC_TOF2BO[i], i+1, ADC_TOF2BI[i]);
        }
      }
      else {
        if ((TDC_TOF2BO[i]>TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<TOF2BO_TDC_max[i]) || (TDC_TOF2BI[i]>TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<TOF2BI_TDC_max[i])) {
         printf("                                        |  TDC TOF2 OutB-%2d: %5d -- InB-%2d: %5d\n", i+1, TDC_TOF2BO[i], i+1, TDC_TOF2BI[i]);
        }
      }
    }
  }

  //double Gap[12][3][2] = {{{0}}};
  double a_fit_GoodLG=0.;                     float b_fit_GoodLG=0.;
  double a_fit_GoodLG_weighted=0.;            float b_fit_GoodLG_weighted=0.;
  double a_fit_TDC_selected_weighted=0.;      float b_fit_TDC_selected_weighted=0.;
  //double a_fit_TDC_Gap_Fibers=0.;             double b_fit_TDC_Gap_Fibers=0.;

  //for(int g=0; g<12; g++)
  //{
  //  Gap[g][0][0] = TOF_Xloc[3*g];
  //  Gap[g][1][0] = TOF_Xloc[3*g+1];
  //  Gap[g][2][0] = TOF_Xloc[3*g+2];
  //
  //  Gap[g][0][1] = TOF_Yloc[3*g];
  //  Gap[g][1][1] = TOF_Yloc[3*g+1];
  //  Gap[g][2][1] = TOF_Yloc[3*g+2];
  //}

  for(int i=0; i<12; i++)
  {
    h_Circle->Fill(Gap[i][0][0], Gap[i][0][1]);
    h_Circle->Fill(Gap[i][1][0], Gap[i][1][1]);
    h_Circle->Fill(Gap[i][2][0], Gap[i][2][1]);
  }

  h_Target_Center->Fill(0., 0.);

  for(int t=0; t<256; t++) h_Target->Fill(Xloc[t],Yloc[t]);

  h_max->Fill(Xloc[max_index], Yloc[max_index]);        h_max_copy->Fill(Xloc[max_index], Yloc[max_index]);
  h_kaon->Fill(Xloc[max_index], Yloc[max_index]);       h_kaon_copy->Fill(Xloc[max_index], Yloc[max_index]);
  h_kaon->Fill(Xloc[max_index2], Yloc[max_index2]);     h_kaon_copy->Fill(Xloc[max_index2], Yloc[max_index2]);
  h_kaon->Fill(Xloc[max_index3], Yloc[max_index3]);     h_kaon_copy->Fill(Xloc[max_index3], Yloc[max_index3]);
  h_kaon->Fill(Xloc[max_index4], Yloc[max_index4]);     h_kaon_copy->Fill(Xloc[max_index4], Yloc[max_index4]);

  // ROTATE_CHANGE
  //if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
  if(gap_to_fit==6 || gap_to_fit==12){
    for(int g=0; g<3; g++) h_TOF1_rotate->Fill(Gap[gap_to_fit_rotate-1][g][0], Gap[gap_to_fit_rotate-1][g][1]);
  }
  else{
    for(int g=0; g<3; g++) h_TOF1->Fill(Gap[gap_to_fit-1][g][0], Gap[gap_to_fit-1][g][1]);
  }

  h_Circle->SetMarkerStyle(5);
  h_Circle->SetMarkerSize(1.2);
  h_Circle->SetLineWidth(2);

  h_Target->SetMarkerStyle(25);         h_Target_Center->SetMarkerStyle(5);
  h_Target->SetMarkerColor(1);          h_Target_Center->SetMarkerColor(1);
  h_Target->SetMarkerSize(0.8);         h_Target_Center->SetMarkerSize(1);

  h_kaon->SetMarkerStyle(21);           h_kaon_copy->SetMarkerStyle(21);
  h_kaon->SetMarkerColor(1);            h_kaon_copy->SetMarkerColor(1);
  h_kaon->SetMarkerSize(0.8);           h_kaon_copy->SetMarkerSize(0.8);

  h_max->SetMarkerStyle(21);            h_max_copy->SetMarkerStyle(21);
  h_max->SetMarkerColor(kViolet+1);     h_max_copy->SetMarkerColor(kViolet+1);
  h_max->SetMarkerSize(1);              h_max_copy->SetMarkerSize(1);

  h_TOF1->SetMarkerStyle(20);           h_TOF1_closest->SetMarkerStyle(20);
  h_TOF1->SetMarkerColor(2);            h_TOF1_closest->SetMarkerColor(3);
  h_TOF1->SetMarkerSize(1.5);           h_TOF1_closest->SetMarkerSize(1.5);

  h_TOF1_rotate->SetMarkerStyle(20);
  h_TOF1_rotate->SetMarkerColor(2);
  h_TOF1_rotate->SetMarkerSize(1.5);

  ell->SetFillStyle(0);                 ell_Target->SetFillStyle(0);                 ell_L1->SetFillStyle(0);
  ell->SetLineColor(6);                 ell_Target->SetLineColor(1);                 ell_L1->SetLineColor(4);
  ell->SetLineWidth(1);                 ell_Target->SetLineWidth(1);                 ell_L1->SetLineWidth(1);

  Gap1l->SetLineWidth(10);
  Gap1l->SetLineColor(15);

  Gap2l->SetLineWidth(10);
  Gap2l->SetLineColor(15);

  Gap3l->SetLineWidth(10);
  Gap3l->SetLineColor(15);

  Gap4l->SetLineWidth(10);
  Gap4l->SetLineColor(15);

  Gap5l->SetLineWidth(10);
  Gap5l->SetLineColor(15);

  Gap6l->SetLineWidth(10);
  Gap6l->SetLineColor(15);

  Gap7l->SetLineWidth(10);
  Gap7l->SetLineColor(15);

  Gap8l->SetLineWidth(10);
  Gap8l->SetLineColor(15);

  Gap9l->SetLineWidth(10);
  Gap9l->SetLineColor(15);

  Gap10l->SetLineWidth(10);
  Gap10l->SetLineColor(15);

  Gap11l->SetLineWidth(10);
  Gap11l->SetLineColor(15);

  Gap12l->SetLineWidth(10);
  Gap12l->SetLineColor(15);

  // Add weight to center point of selected TOF1 (gap_to_fit_rotate)
  //for(int i = 0; i<3; i++){

  //    vec_xx_lepton.push_back(Gap[gap_to_fit-1][1][0]);
  //    vec_ex_lepton.push_back(TOF1_Errors_X[gap_to_fit-1][1]);

  //    vec_yy_lepton.push_back(Gap[gap_to_fit-1][1][1]);
  //    vec_ey_lepton.push_back(TOF1_Errors_Y[gap_to_fit-1][1]);

  //    cout << "FIT2_TOF1 : " << i << endl;
  //}

  //double a_lepton_2 = 0.;
  //double b_lepton_2 = 0.;
  //double Chis_lepton_2 = 0.;
  //double ndf_lepton_2 = 0.;

  //cout << "ALBANE SIZE NEW : " << vec_xx_lepton.size() << endl;

  //TGraph *gr_lepton_2 = new TGraphErrors(vec_xx_lepton.size(),&vec_xx_lepton[0],&vec_yy_lepton[0],
  //  &vec_ex_lepton[0],&vec_ey_lepton[0]);

  //TF1 *func_lepton_fit_2 = new TF1("lepton_fit_2", "pol1");


  //gr_lepton_2->SetLineColor(2);
  //gr_lepton_2->SetLineWidth(4);
  //gr_lepton_2->SetMarkerColor(4);
  //gr_lepton_2->SetMarkerStyle(21);
  //gr_lepton_2->GetXaxis()->SetLimits(-50.,50.);
  //gr_lepton_2->GetYaxis()->SetRangeUser(-50.,50.);

  //gr_lepton_2->Fit("lepton_fit_2","QC");
  //func_lepton_fit_2 = gr_lepton_2->GetFunction("lepton_fit_2");
  //a_lepton_2 = func_lepton_fit_2->GetParameter(1);
  //b_lepton_2 = func_lepton_fit_2->GetParameter(0);
  //Chis_lepton_2 = func_lepton_fit_2->GetChisquare();
  //ndf_lepton_2 = func_lepton_fit_2->GetNDF();


  //a_fit_TDC_Gap_Fibers = a_lepton_2;
  //b_fit_TDC_Gap_Fibers = b_lepton_2;

  //cout << endl;
  //cout << "//////////////  FIT 2  " << endl;
  //cout << endl;
  //cout << "a_lepton_2 : " << a_lepton_2 << endl;
  //cout << "b_lepton_2 : " << b_lepton_2 << endl;
  //cout << "Chis_lepton_2 : " << Chis_lepton_2 << endl;
  //cout << "ndf_lepton_2 : " << ndf_lepton_2 << endl;
  //cout << endl;
  //cout << "////////////////////" << endl;
  //cout << endl;



  //cout << "ALBANE 1 : " << a_lepton_fit_3 << "   " << b_lepton_fit_3 << endl;

  vec_xx_lepton.clear();            vec_ex_lepton.clear();
  vec_yy_lepton.clear();            vec_ey_lepton.clear();
  vec_xx_lepton_rotate.clear();     vec_ex_lepton_rotate.clear();
  vec_yy_lepton_rotate.clear();     vec_ey_lepton_rotate.clear();
  vec_lepton_size.clear();

  //gr_lepton_2->Set(0);

  //vec_xx_kaon.clear();      vec_ex_kaon.clear();
  //vec_yy_kaon.clear();      vec_ey_kaon.clear();


    /*// Fill KStop bars
    for(int i = 0; i<256; i++){
      if(k_stop_bar_initial[i]){
      //if(k_stop_bar[i]){
        cout << "////////////////////////////   DEBUG2 : " << i << endl;
        vec_xx_kaon.push_back(Xloc[i]);
        vec_ex_kaon.push_back(1.6);

        vec_yy_kaon.push_back(Yloc[i]);
        vec_ey_kaon.push_back(1.6);
      }
    } */

    /*TGraph *gr_kaon_final = new TGraphErrors(vec_xx_kaon.size(),&vec_xx_kaon[0],&vec_yy_kaon[0],
      &vec_ex_kaon[0],&vec_ey_kaon[0]);
    gr_kaon_final->SetMarkerStyle(21);
    gr_kaon_final->SetMarkerColor(kBlue-6);
    gr_kaon_final->SetMarkerSize(0.8);
    gr_kaon_final->GetXaxis()->SetLimits(-50.,50.);
    gr_kaon_final->GetYaxis()->SetRangeUser(-50.,50.);
*/

    /*cout << endl;
    cout << endl;
    cout << endl;
    cout << endl;
    cout << endl;
    TF1 *gr_kaon_final_fit = new TF1("test1", "pol1");
    double a_fit_final_kaon = 0.;
    double b_fit_final_kaon = 0.;
    //int a_fit_kaon = 0;
    //int b_fit_kaon = 0;
    cout << endl;
    cout << endl;
    cout << endl;
    cout << endl;
    cout << endl;
    cout << endl;
*/
    //cout << "SEB1 : " << a_fit_kaon << "  " << b_fit_kaon << endl;

    /*if(vec_xx_kaon.size()>0){
      //gr_kaon_final->Fit("pol1","Q");
      gr_kaon_final->Fit("test1","Q");
      gr_kaon_final_fit = gr_kaon_final->GetFunction("test1");
      //gr_kaon_final_fit = gr_kaon_final->GetFunction("pol1");
      gr_kaon_final_fit->SetLineWidth(2);
      gr_kaon_final_fit->SetLineColor(4);
      a_fit_final_kaon = gr_kaon_final_fit->GetParameter(1);
      b_fit_final_kaon = gr_kaon_final_fit->GetParameter(0);
      //cout << "SEB3 : " << gr_kaon_fit->GetParameter(1) << "  " << gr_kaon_fit->GetParameter(0) << endl;
    }*/


  //vec_xx_lepton.push_back(Xloc[49]);    vec_xx_lepton.push_back(Yloc[49]);

  //a_fit_TDC_Gap_Fibers = a_lepton_fit_1;
  //b_fit_TDC_Gap_Fibers = b_lepton_fit_1;

  //cout << "TEST TEST TEST : " << a_fit_TDC_Gap_Fibers << "   " << b_fit_TDC_Gap_Fibers << endl;

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

          //vec_xx_lepton_rotate.push_back(Xloc[TARGET_Rotated_index[i]]);
          //vec_ex_lepton_rotate.push_back(TARGET_Errors_X);

          //vec_yy_lepton_rotate.push_back(Yloc[TARGET_Rotated_index[i]]);
          //vec_ey_lepton_rotate.push_back(TARGET_Errors_Y);

          //vec_bar.push_back(i);
          //vec_bar_rotate.push_back(TARGET_Rotated_index[i]);
          //vec_lepton_size.push_back(Xloc[i]);
          //vec_lepton_rotate_size.push_back(Xloc[i]);

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

              //vec_xx_lepton_rotate.push_back(Xloc[TARGET_Rotated_index[i]]);
              //vec_ex_lepton_rotate.push_back(TARGET_Errors_X);

              //vec_yy_lepton_rotate.push_back(Yloc[TARGET_Rotated_index[i]]);
              //vec_ey_lepton_rotate.push_back(TARGET_Errors_Y);

              //vec_bar.push_back(i);
              //vec_bar_rotate.push_back(TARGET_Rotated_index[i]);
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

              //vec_xx_lepton_rotate.push_back(Xloc[TARGET_Rotated_index[i]]);
              //vec_ex_lepton_rotate.push_back(TARGET_Errors_X);

              //vec_yy_lepton_rotate.push_back(Yloc[TARGET_Rotated_index[i]]);
              //vec_ey_lepton_rotate.push_back(TARGET_Errors_Y);

              //vec_bar.push_back(i);
              //vec_bar_rotate.push_back(TARGET_Rotated_index[i]);
            }
          }
        }

        if(ADC_High_TARGET[lepton_bars[i]]<0 && ADC_Low_TARGET[lepton_bars[i]]>0 && Switch==1){
          vec_xx_lepton.push_back(Xloc[TARGET_Rotated_index[lepton_bars[i]]]);
          vec_ex_lepton.push_back(TARGET_Errors_X);

          vec_yy_lepton.push_back(Yloc[TARGET_Rotated_index[lepton_bars[i]]]);
          vec_ey_lepton.push_back(TARGET_Errors_Y);

          lepton_bars_rotate.push_back(lepton_bars[i]);

          //vec_xx_lepton_rotate.push_back(Xloc[TARGET_Rotated_index[i]]);
          //vec_ex_lepton_rotate.push_back(TARGET_Errors_X);

          //vec_yy_lepton_rotate.push_back(Yloc[TARGET_Rotated_index[i]]);
          //vec_ey_lepton_rotate.push_back(TARGET_Errors_Y);

          //vec_bar.push_back(i);
          //vec_bar_rotate.push_back(TARGET_Rotated_index[i]);
          //vec_lepton_size.push_back(Xloc[i]);
          //vec_lepton_rotate_size.push_back(Xloc[i]);
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

          //vec_xx_lepton_rotate.push_back(Xloc[TARGET_Rotated_index[i]]);
          //vec_ex_lepton_rotate.push_back(TARGET_Errors_X);

          //vec_yy_lepton_rotate.push_back(Yloc[TARGET_Rotated_index[i]]);
          //vec_ey_lepton_rotate.push_back(TARGET_Errors_Y);

          //vec_bar.push_back(i);
          //vec_bar_rotate.push_back(TARGET_Rotated_index[i]);
          //vec_lepton_size.push_back(Xloc[i]);
          //vec_lepton_rotate_size.push_back(Xloc[i]);

          if(IsIn(lepton_bars[i],channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],
                    channel[gap_to_fit-1][2], channel[gap_to_fit-1][3],
                    channel[gap_to_fit-1][4], channel[gap_to_fit-1][5],
                    channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){

            vec_xx_lepton.push_back(Xloc[lepton_bars[i]]);
            vec_ex_lepton.push_back(TARGET_Errors_X);

            vec_yy_lepton.push_back(Yloc[lepton_bars[i]]);
            vec_ey_lepton.push_back(TARGET_Errors_Y);

            //vec_xx_lepton_rotate.push_back(Xloc[TARGET_Rotated_index[i]]);
            //vec_ex_lepton_rotate.push_back(TARGET_Errors_X);

            //vec_yy_lepton_rotate.push_back(Yloc[TARGET_Rotated_index[i]]);
            //vec_ey_lepton_rotate.push_back(TARGET_Errors_Y);

            //vec_bar.push_back(i);
            //vec_bar_rotate.push_back(TARGET_Rotated_index[i]);
          }
        }

        if(ADC_High_TARGET[lepton_bars[i]]<0 && ADC_Low_TARGET[lepton_bars[i]]>0 && Switch==1){
          vec_xx_lepton.push_back(Xloc[lepton_bars[i]]);
          vec_ex_lepton.push_back(TARGET_Errors_X);

          vec_yy_lepton.push_back(Yloc[lepton_bars[i]]);
          vec_ey_lepton.push_back(TARGET_Errors_Y);

          //vec_xx_lepton_rotate.push_back(Xloc[TARGET_Rotated_index[i]]);
          //vec_ex_lepton_rotate.push_back(TARGET_Errors_X);

          //vec_yy_lepton_rotate.push_back(Yloc[TARGET_Rotated_index[i]]);
          //vec_ey_lepton_rotate.push_back(TARGET_Errors_Y);

          //vec_bar.push_back(i);
          //vec_bar_rotate.push_back(TARGET_Rotated_index[i]);
          //vec_lepton_size.push_back(Xloc[i]);
          //vec_lepton_rotate_size.push_back(Xloc[i]);
        }
      }
    }
  }





/*  for(int i = 0; i<256; i++){
    if(distance_to_line(Xloc[i],Yloc[i],a_lepton_fit_1,b_lepton_fit_1) <= max_dist && TARGET_High_has_neighbours[i] && !k_stop_bar[i]){
      if(ADC_High_TARGET[i]>Angle_ADC_cut && has_TDC_hit[i]){

          vec_xx_lepton.push_back(Xloc[i]);
          vec_ex_lepton.push_back(TARGET_Errors_X);

          vec_yy_lepton.push_back(Yloc[i]);
          vec_ey_lepton.push_back(TARGET_Errors_Y);

          vec_xx_lepton_rotate.push_back(Xloc[TARGET_Rotated_index[i]]);
          vec_ex_lepton_rotate.push_back(TARGET_Errors_X);

          vec_yy_lepton_rotate.push_back(Yloc[TARGET_Rotated_index[i]]);
          vec_ey_lepton_rotate.push_back(TARGET_Errors_Y);

          vec_bar.push_back(i);
          vec_bar_rotate.push_back(TARGET_Rotated_index[i]);
          vec_lepton_size.push_back(Xloc[i]);
          vec_lepton_rotate_size.push_back(Xloc[i]);

          if(IsIn(i,channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],
                    channel[gap_to_fit-1][2], channel[gap_to_fit-1][3],
                    channel[gap_to_fit-1][4], channel[gap_to_fit-1][5],
                    channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){
            vec_xx_lepton.push_back(Xloc[i]);
            vec_ex_lepton.push_back(TARGET_Errors_X);

            vec_yy_lepton.push_back(Yloc[i]);
            vec_ey_lepton.push_back(TARGET_Errors_Y);

            vec_xx_lepton_rotate.push_back(Xloc[TARGET_Rotated_index[i]]);
            vec_ex_lepton_rotate.push_back(TARGET_Errors_X);

            vec_yy_lepton_rotate.push_back(Yloc[TARGET_Rotated_index[i]]);
            vec_ey_lepton_rotate.push_back(TARGET_Errors_Y);

            vec_bar.push_back(i);
            vec_bar_rotate.push_back(TARGET_Rotated_index[i]);
          }
      }
      if(ADC_High_TARGET[i]<0 && ADC_Low_TARGET[i]>0 && Switch==1){
        vec_xx_lepton.push_back(Xloc[i]);
        vec_ex_lepton.push_back(TARGET_Errors_X);

        vec_yy_lepton.push_back(Yloc[i]);
        vec_ey_lepton.push_back(TARGET_Errors_Y);

        vec_xx_lepton_rotate.push_back(Xloc[TARGET_Rotated_index[i]]);
        vec_ex_lepton_rotate.push_back(TARGET_Errors_X);

        vec_yy_lepton_rotate.push_back(Yloc[TARGET_Rotated_index[i]]);
        vec_ey_lepton_rotate.push_back(TARGET_Errors_Y);

        vec_bar.push_back(i);
        vec_bar_rotate.push_back(TARGET_Rotated_index[i]);
        vec_lepton_size.push_back(Xloc[i]);
        vec_lepton_rotate_size.push_back(Xloc[i]);
      }
    }
  }

*/





    double a_lepton_fit_2 = 0.;
    double b_lepton_fit_2 = 0.;
    double Chis_lepton_fit_2 = 0.;
    int ndf_lepton_fit_2 = 0;

    TGraph *gr_lepton_2;
    gr_lepton_2 = new TGraphErrors(vec_xx_lepton.size(),&vec_xx_lepton[0],&vec_yy_lepton[0],
      &vec_ex_lepton[0],&vec_ey_lepton[0]);

    TF1 *func_lepton_fit_2;

  if(vec_xx_lepton.size()>0){

    gr_lepton_2->SetMarkerStyle(21);
    gr_lepton_2->SetMarkerColor(2);
    gr_lepton_2->SetMarkerSize(0.8);
    gr_lepton_2->GetXaxis()->SetLimits(-50.,50.);
    gr_lepton_2->GetYaxis()->SetRangeUser(-50.,50.);

    func_lepton_fit_2 = new TF1("lepton_fit_2", "pol1");

    if(gap_to_fit==1 || gap_to_fit==2 || gap_to_fit==7 || gap_to_fit==8){
      func_lepton_fit_2->SetParameter(0,0);
      func_lepton_fit_2->SetParameter(1,1);
    }
    if(gap_to_fit==4 || gap_to_fit==5 || gap_to_fit==10 || gap_to_fit==11){
      func_lepton_fit_2->SetParameter(0,0);
      func_lepton_fit_2->SetParameter(1,-1);
    }
    else{
      func_lepton_fit_2->SetParameter(0,0);
      func_lepton_fit_2->SetParameter(1,1);
    }


    func_lepton_fit_2->SetParLimits(0,-50,50);
    func_lepton_fit_2->SetParLimits(1,-50,50);


    gr_lepton_2->Fit("lepton_fit_2","Q");
    func_lepton_fit_2 = gr_lepton_2->GetFunction("lepton_fit_2");
    func_lepton_fit_2->SetLineWidth(2);
    func_lepton_fit_2->SetLineColor(2);

    a_lepton_fit_2 = func_lepton_fit_2->GetParameter(1);
    b_lepton_fit_2 = func_lepton_fit_2->GetParameter(0);
    Chis_lepton_fit_2 = func_lepton_fit_2->GetChisquare();
    ndf_lepton_fit_2 = func_lepton_fit_2->GetNDF();
  }
  else cout << "Empty Fit 2" << endl;


  /*cout << endl;
  cout << "//////////////  FIT 2  " << endl;
  cout << endl;
  cout << "a_lepton_2 : " << a_lepton_fit_2 << endl;
  cout << "b_lepton_2 : " << b_lepton_fit_2 << endl;
  cout << "Chis_lepton_2 : " << Chis_lepton_fit_2 << endl;
  cout << "ndf_lepton_2 : " << ndf_lepton_fit_2 << endl;
  cout << endl;
  cout << "////////////////////" << endl;
  cout << endl;
*/

  //TGraph *gr2_Leptons_rotate = new TGraph(vec_xx_lepton_rotate.size(),&vec_xx_lepton_rotate[0],&vec_yy_lepton_rotate[0]);
  TGraph *gr2_Leptons_rotate = new TGraphErrors(vec_xx_lepton_rotate.size(),&vec_xx_lepton_rotate[0],&vec_yy_lepton_rotate[0],
    &vec_ex_lepton[0],&vec_ey_lepton[0]);
  gr2_Leptons_rotate->SetMarkerStyle(25);
  gr2_Leptons_rotate->SetMarkerColor(2);
  gr2_Leptons_rotate->SetMarkerSize(0.8);
  gr2_Leptons_rotate->GetXaxis()->SetLimits(-50.,50.);
  gr2_Leptons_rotate->GetYaxis()->SetRangeUser(-50.,50.);

  //a_fit_TDC_Gap_Fibers = a_gr_Leptons;
  //b_fit_TDC_Gap_Fibers = b_gr_Leptons;

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

  //x_int_TDC_Gap_Fibers[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);
  //x_int_TDC_Gap_Fibers[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);
  //y_int_TDC_Gap_Fibers[0] = y1_int(x_int_TDC_Gap_Fibers[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);
  //y_int_TDC_Gap_Fibers[1] = y2_int(x_int_TDC_Gap_Fibers[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);

  x_int_TDC_Gap_Fibers[0] = intersectx1(a_lepton_fit_2, b_lepton_fit_2, R_TOF1);
  x_int_TDC_Gap_Fibers[1] = intersectx2(a_lepton_fit_2, b_lepton_fit_2, R_TOF1);
  y_int_TDC_Gap_Fibers[0] = y1_int(x_int_TDC_Gap_Fibers[0], a_lepton_fit_2, b_lepton_fit_2);
  y_int_TDC_Gap_Fibers[1] = y2_int(x_int_TDC_Gap_Fibers[1], a_lepton_fit_2, b_lepton_fit_2);

  //x_int_TDC_Gap_Fibers_SFT[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_SFT_L1);
  //x_int_TDC_Gap_Fibers_SFT[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_SFT_L1);
  //y_int_TDC_Gap_Fibers_SFT[0] = y1_int(x_int_TDC_Gap_Fibers_SFT[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);
  //y_int_TDC_Gap_Fibers_SFT[1] = y2_int(x_int_TDC_Gap_Fibers_SFT[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);

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














/*  for(int ii=0; ii<3; ii++){
    dist2_TDC_selected[ii] = distance(x_TDC_Gap_Fibers, y_TDC_Gap_Fibers, Gap[gap_to_fit-1][ii][0], Gap[gap_to_fit-1][ii][1]);
    dist2_GoodLG[ii] = distance(x_GoodLG_intersect1, y_GoodLG_intersect1, Gap[gap_to_fit-1][ii][0], Gap[gap_to_fit-1][ii][1]);

    if(dist2_TDC_selected[ii] <= dist2_TDC_selected_min){
      dist2_TDC_selected_min = dist2_TDC_selected[ii];
      selected_TDC_selected = ii;
    }

    if(dist2_GoodLG[ii] <= dist2_GoodLG_min) dist2_GoodLG_min = dist2_GoodLG[ii];
  }
*/
  h_TOF1_closest->Fill(Gap[gap_to_fit_rotate-1][selected_TDC_selected][0], Gap[gap_to_fit_rotate-1][selected_TDC_selected][1]);

  //cout << "RUQUIER : " << gap_to_fit-1 << "    " << selected_TDC_selected << endl;

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


  h_int_TDC_Gap_Fibers->SetMarkerStyle(20);   h_int_TDC_Gap_Fibers_SFT->SetMarkerStyle(20);
  h_int_TDC_Gap_Fibers->SetMarkerColor(4);    h_int_TDC_Gap_Fibers_SFT->SetMarkerColor(4);
  h_int_TDC_Gap_Fibers->SetMarkerSize(0.8);   h_int_TDC_Gap_Fibers_SFT->SetMarkerSize(0.8);

  h_int_TDC_Gap_Fibers_rotate->SetMarkerStyle(20);   h_int_TDC_Gap_Fibers_SFT_rotate->SetMarkerStyle(20);
  h_int_TDC_Gap_Fibers_rotate->SetMarkerColor(4);    h_int_TDC_Gap_Fibers_SFT_rotate->SetMarkerColor(4);
  h_int_TDC_Gap_Fibers_rotate->SetMarkerSize(0.8);   h_int_TDC_Gap_Fibers_SFT_rotate->SetMarkerSize(0.8);

  h_int_TDC_TARGET->SetMarkerStyle(20);       h_int_TDC_TARGET_rotate->SetMarkerStyle(20);
  h_int_TDC_TARGET->SetMarkerColor(1);        h_int_TDC_TARGET_rotate->SetMarkerColor(1);
  h_int_TDC_TARGET->SetMarkerSize(0.8);       h_int_TDC_TARGET_rotate->SetMarkerSize(0.8);

  double ParError = 999.99;
  double ChiS = 0.0;
  int ndf = 0;
  //double a_gr_Leptons = 0.;
  //double b_gr_Leptons = 0.;


  TGraph *gr_lepton_3 = new TGraphErrors(vec_xx_lepton.size(),&vec_xx_lepton[0],&vec_yy_lepton[0],
    &vec_ex_lepton[0],&vec_ey_lepton[0]);

  gr_lepton_3->SetMarkerStyle(21);                    //gr3_Leptons_test->SetMarkerStyle(25);
  gr_lepton_3->SetMarkerColor(2);                     //gr3_Leptons_test->SetMarkerColor(5);
  gr_lepton_3->SetMarkerSize(0.8);                    //gr3_Leptons_test->SetMarkerSize(0.8);
  gr_lepton_3->GetXaxis()->SetLimits(-50.,50.);       //gr3_Leptons_test->GetXaxis()->SetLimits(-50.,50.);
  gr_lepton_3->GetYaxis()->SetRangeUser(-50.,50.);    //gr3_Leptons_test->GetYaxis()->SetRangeUser(-50.,50.);

  double a_lepton_fit_3 = 0.;
  double b_lepton_fit_3 = 0.;
  double Chis_lepton_fit_3 = 0.;
  int ndf_lepton_fit_3 = 0;

  TF1 *func_lepton_fit_3 = new TF1("lepton_fit_3", "pol1");

  if(gap_to_fit==1 || gap_to_fit==2 || gap_to_fit==7 || gap_to_fit==8){
    func_lepton_fit_3->SetParameter(0,0);
    func_lepton_fit_3->SetParameter(1,1);
  }
  if(gap_to_fit==4 || gap_to_fit==5 || gap_to_fit==10 || gap_to_fit==11){
    func_lepton_fit_3->SetParameter(0,0);
    func_lepton_fit_3->SetParameter(1,-1);
  }
  else{
      func_lepton_fit_3->SetParameter(0,0);
      func_lepton_fit_3->SetParameter(1,1);
    }


  func_lepton_fit_3->SetParLimits(0,-50,50);
  func_lepton_fit_3->SetParLimits(1,-50,50);

  gr_lepton_3->Fit("lepton_fit_3","Q");
  func_lepton_fit_3 = gr_lepton_3->GetFunction("lepton_fit_3");
  func_lepton_fit_3->SetLineWidth(2);
  func_lepton_fit_3->SetLineColor(2);


  a_lepton_fit_3 = func_lepton_fit_3->GetParameter(1);
  b_lepton_fit_3 = func_lepton_fit_3->GetParameter(0);

  ParError = func_lepton_fit_3->GetParError(1);
  Chis_lepton_fit_3 = func_lepton_fit_3->GetChisquare();
  ndf_lepton_fit_3 = func_lepton_fit_3->GetNDF();

  /*cout << endl;
  cout << "//////////////  FIT 3  " << endl;
  cout << endl;
  cout << "a_lepton_3 : " << a_lepton_fit_3 << endl;
  cout << "b_lepton_3 : " << b_lepton_fit_3 << endl;
  cout << "Chis_lepton_3 : " << Chis_lepton_fit_3 << endl;
  cout << "ndf_lepton_3 : " << ndf_lepton_fit_3 << endl;
  cout << endl;
  cout << "////////////////////" << endl;
  cout << endl;
*/

  if(gap_to_fit==6 || gap_to_fit==12){
    for(unsigned int j=0; j<lepton_bars_rotate.size(); j++){

      vec_xx_lepton_rotate.push_back(Xloc[lepton_bars_rotate[j]]);
      vec_ex_lepton_rotate.push_back(TARGET_Errors_X);

      vec_yy_lepton_rotate.push_back(Yloc[lepton_bars_rotate[j]]);
      vec_ey_lepton_rotate.push_back(TARGET_Errors_Y);

    }
  }


  TGraph *gr3_Leptons_rotate = new TGraphErrors(vec_xx_lepton_rotate.size(),&vec_xx_lepton_rotate[0],&vec_yy_lepton_rotate[0],
    &vec_ex_lepton_rotate[0],&vec_ey_lepton_rotate[0]);


  gr3_Leptons_rotate->SetMarkerStyle(21);
  gr3_Leptons_rotate->SetMarkerColor(2);
  gr3_Leptons_rotate->SetMarkerSize(0.8);
  gr3_Leptons_rotate->GetXaxis()->SetLimits(-50.,50.);
  gr3_Leptons_rotate->GetYaxis()->SetRangeUser(-50.,50.);

  //a_gr_Leptons = a_lepton_fit_3;
  //b_gr_Leptons = b_lepton_fit_3;

  //a_fit_TDC_Gap_Fibers = a_gr_Leptons;
  //b_fit_TDC_Gap_Fibers = b_gr_Leptons;

  //x_int_TDC_selected_weighted[0] = intersectx1(a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted, R_TOF1);
  //x_int_TDC_selected_weighted[1] = intersectx2(a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted, R_TOF1);

  //y_int_TDC_selected_weighted[0] = y1_int(x_int_TDC_selected_weighted[0], a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted);
  //y_int_TDC_selected_weighted[1] = y2_int(x_int_TDC_selected_weighted[1], a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted);

  //x_int_GoodLG_weighted[0] = intersectx1(a_fit_GoodLG_weighted, b_fit_GoodLG_weighted, R_TOF1);
  //x_int_GoodLG_weighted[1] = intersectx2(a_fit_GoodLG_weighted, b_fit_GoodLG_weighted, R_TOF1);

  //y_int_GoodLG_weighted[0] = y1_int(x_int_GoodLG_weighted[0], a_fit_GoodLG_weighted, b_fit_GoodLG_weighted);
  //y_int_GoodLG_weighted[1] = y2_int(x_int_GoodLG_weighted[1], a_fit_GoodLG_weighted, b_fit_GoodLG_weighted);

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

  //float x_TDC_selected_weighted_intersect1=0.;   float y_TDC_selected_weighted_intersect1=0.;
  //float x_GoodLG_weighted_intersect1=0.;         float y_GoodLG_weighted_intersect1=0.;
  float x_TDC_Gap_Fibers_intersect1=0.;          float y_TDC_Gap_Fibers_intersect1=0.;
  float x_TDC_Gap_Fibers_SFT_intersect1=0.;      float y_TDC_Gap_Fibers_SFT_intersect1=0.;
  float x_TARGET_intersect=0;                    float y_TARGET_intersect=0;

  //float dist1_TDC_selected_weighted[2];
  //float dist1_GoodLG_weighted[2];
  float dist1_TARGET_intersect[2];

  for(int i=0; i<2; i++){
    if(gap_to_fit==6){
      //dist1_TDC_selected_weighted[i] = distance(x_int_TDC_selected_weighted[i], y_int_TDC_selected_weighted[i], Gap[3-1][1][0], Gap[3-1][1][1]);
      //dist1_GoodLG_weighted[i] = distance(x_int_GoodLG_weighted[i], y_int_GoodLG_weighted[i], Gap[3-1][1][0], Gap[3-1][1][1]);
      dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[3-1][1][0], Gap[3-1][1][1]);
      dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[3-1][1][0], Gap[3-1][1][1]);
      dist1_TARGET_intersect[i] = distance(x_int_TARGET[i], y_int_TARGET[i], Gap[3-1][1][0], Gap[3-1][1][1]);
    }
    else if(gap_to_fit==12){
      //dist1_TDC_selected_weighted[i] = distance(x_int_TDC_selected_weighted[i], y_int_TDC_selected_weighted[i], Gap[9-1][1][0], Gap[9-1][1][1]);
      //dist1_GoodLG_weighted[i] = distance(x_int_GoodLG_weighted[i], y_int_GoodLG_weighted[i], Gap[9-1][1][0], Gap[9-1][1][1]);
      dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[9-1][1][0], Gap[9-1][1][1]);
      dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[9-1][1][0], Gap[9-1][1][1]);
      dist1_TARGET_intersect[i] = distance(x_int_TARGET[i], y_int_TARGET[i], Gap[9-1][1][0], Gap[9-1][1][1]);
    }
    else{
      //dist1_TDC_selected_weighted[i] = distance(x_int_TDC_selected_weighted[i], y_int_TDC_selected_weighted[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      //dist1_GoodLG_weighted[i] = distance(x_int_GoodLG_weighted[i], y_int_GoodLG_weighted[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TARGET_intersect[i] = distance(x_int_TARGET[i], y_int_TARGET[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
    }
  }

/*  if(dist1_TDC_selected_weighted[0] < dist1_TDC_selected_weighted[1]){
    x_TDC_selected_weighted_intersect1 = x_int_TDC_selected_weighted[0];
    y_TDC_selected_weighted_intersect1 = y_int_TDC_selected_weighted[0];
  }
  else if(dist1_TDC_selected_weighted[1] < dist1_TDC_selected_weighted[0]){
    x_TDC_selected_weighted_intersect1 = x_int_TDC_selected_weighted[1];
    y_TDC_selected_weighted_intersect1 = y_int_TDC_selected_weighted[1];
  }
  else cout << "ERROR !" << endl;
*/

/*  if(dist1_GoodLG_weighted[0] < dist1_GoodLG_weighted[1]){
    x_GoodLG_weighted_intersect1 = x_int_GoodLG_weighted[0];
    y_GoodLG_weighted_intersect1 = y_int_GoodLG_weighted[0];
  }
  else if(dist1_GoodLG_weighted[1] < dist1_GoodLG_weighted[0]){
    x_GoodLG_weighted_intersect1 = x_int_GoodLG_weighted[1];
    y_GoodLG_weighted_intersect1 = y_int_GoodLG_weighted[1];
  }
  else cout << "ERROR !" << endl;
*/
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
  }
  else if(dist1_TARGET_intersect[1] < dist1_TARGET_intersect[0]){
    x_TARGET_intersect = x_int_TARGET[1];
    y_TARGET_intersect = y_int_TARGET[1];
  }
  else cout << "ERROR !" << endl;


  //if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
  if(gap_to_fit==6 || gap_to_fit==12){
    //h_int_TDC_Gap_Fibers_rotate->Fill(-y_TDC_Gap_Fibers_intersect1, x_TDC_Gap_Fibers_intersect1);
    //h_int_TDC_Gap_Fibers_SFT_rotate->Fill(-y_TDC_Gap_Fibers_SFT_intersect1, x_TDC_Gap_Fibers_SFT_intersect1);
    h_int_TDC_Gap_Fibers_rotate->Fill(y_TDC_Gap_Fibers_intersect1, -x_TDC_Gap_Fibers_intersect1);
    h_int_TDC_Gap_Fibers_SFT_rotate->Fill(y_TDC_Gap_Fibers_SFT_intersect1, -x_TDC_Gap_Fibers_SFT_intersect1);
    //h_int_TDC_Gap_Fibers_rotate->Fill(x_TDC_Gap_Fibers_intersect1, y_TDC_Gap_Fibers_intersect1);
    //h_int_TDC_Gap_Fibers_SFT_rotate->Fill(x_TDC_Gap_Fibers_SFT_intersect1, y_TDC_Gap_Fibers_SFT_intersect1);
  }
  else{
    h_int_TDC_Gap_Fibers->Fill(x_TDC_Gap_Fibers_intersect1, y_TDC_Gap_Fibers_intersect1);
    h_int_TDC_Gap_Fibers_SFT->Fill(x_TDC_Gap_Fibers_SFT_intersect1, y_TDC_Gap_Fibers_SFT_intersect1);
  }
  h_int_TDC_TARGET->Fill(x_TARGET_intersect, y_TARGET_intersect);

  //TLine *CenterLine_TDC = new TLine(0., 0., x_TDC_Gap_Fibers_intersect1, y_TDC_Gap_Fibers_intersect1);
  //TLine *CenterLine_GoodLG = new TLine(0., 0., x_GoodLG_weighted_intersect1, y_GoodLG_weighted_intersect1);


  //CenterLine_TDC->SetLineWidth(2);    CenterLine_GoodLG->SetLineWidth(2);
  //CenterLine_TDC->SetLineColor(3);    CenterLine_GoodLG->SetLineColor(3);

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

  //a_final_TDC = y_TDC_selected_weighted_intersect1 / x_TDC_selected_weighted_intersect1;
  //a_final_GoodLG = y_GoodLG_weighted_intersect1 / x_GoodLG_weighted_intersect1;
  //a_final_Gap_Fibers = y_TDC_Gap_Fibers_intersect1 / x_TDC_Gap_Fibers_intersect1;

  //a_final_guide = (y_TDC_Gap_Fibers_intersect1 - y_TARGET_intersect) / (x_TDC_Gap_Fibers_intersect1 - x_TARGET_intersect);

  if(gap_to_fit==6 || gap_to_fit==12) a_final_guide = a_lepton_fit_3;
  else a_final_guide = a_lepton_fit_3;

  //tanalpha_GoodLG = a_final_GoodLG;       tanalpha_Gap_Fibers = a_final_Gap_Fibers;
  //alpha_GoodLG = atan(tanalpha_GoodLG);   alpha_Gap_Fibers = atan(tanalpha_Gap_Fibers);

  tanalpha_guide = a_final_guide;
  alpha_guide = atan(tanalpha_guide);

  //cout << "COMPARISON : " << "  " << "a_final_guide = " << a_final_guide << endl;
  //cout << "COMPARISON : " << "  " << "a_lepton_fit_3 = " << a_lepton_fit_3 << endl;

  /// Determination of Delta Phi
  double Delta_phi = 999.99;  double Delta_phi_deg = 999.99;
  Delta_phi = sin(alpha_guide)*cos(alpha_guide)*(ParError/a_final_guide);

  Delta_phi_deg = (180/PI)*Delta_phi;

  int Axis_Vector_Length = 10;
  TArrow *x_guide;
  TArrow *y_guide;
  TArrow *x_guide_rotate;
  TArrow *y_guide_rotate;

  x_guide = new TArrow(x_TARGET_intersect, y_TARGET_intersect, x_TARGET_intersect + Axis_Vector_Length, y_TARGET_intersect, 0.005, "|>");
  y_guide = new TArrow(x_TARGET_intersect, y_TARGET_intersect, x_TARGET_intersect, y_TARGET_intersect + Axis_Vector_Length, 0.005, "|>");

  x_guide->SetLineWidth(2);     y_guide->SetLineWidth(2);
  x_guide->SetLineColor(4);     y_guide->SetLineColor(4);

  //x_guide_rotate = new TArrow(-y_TARGET_intersect, x_TARGET_intersect, -y_TARGET_intersect + Axis_Vector_Length, x_TARGET_intersect, 0.005, "|>");
  //y_guide_rotate = new TArrow(-y_TARGET_intersect, x_TARGET_intersect, -y_TARGET_intersect, x_TARGET_intersect + Axis_Vector_Length, 0.005, "|>");

  x_guide_rotate = new TArrow(y_TARGET_intersect, -x_TARGET_intersect, y_TARGET_intersect + Axis_Vector_Length, -x_TARGET_intersect, 0.005, "|>");
  y_guide_rotate = new TArrow(y_TARGET_intersect, -x_TARGET_intersect, y_TARGET_intersect, -x_TARGET_intersect + Axis_Vector_Length, 0.005, "|>");

  double x_exit_rotate;
  double y_exit_rotate;

  x_exit_rotate = -y_TARGET_intersect;
  y_exit_rotate = x_TARGET_intersect;

  x_guide_rotate->SetLineWidth(2);     y_guide_rotate->SetLineWidth(2);
  x_guide_rotate->SetLineColor(4);     y_guide_rotate->SetLineColor(4);



  //if(x_GoodLG_weighted_intersect1 < 0) angle_final_GoodLG = 180. + alpha_GoodLG * (180./PI);

  //if(x_GoodLG_weighted_intersect1 >= 0){
  //  if(y_GoodLG_weighted_intersect1 >= 0) angle_final_GoodLG = alpha_GoodLG * (180./PI);
  //  else angle_final_GoodLG = alpha_GoodLG * (180./PI) + 360.;
  //}

  //if(x_TDC_Gap_Fibers_intersect1 < 0)  angle_final_Gap_Fibers = 180. + alpha_Gap_Fibers * (180./PI);

  //if(x_TDC_Gap_Fibers_intersect1 >= 0){
  //  if(y_TDC_Gap_Fibers_intersect1 >= 0)  angle_final_Gap_Fibers = alpha_Gap_Fibers * (180./PI);
  //  else angle_final_Gap_Fibers = alpha_Gap_Fibers * (180./PI) + 360.;
  //}

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



  char Angle_guide_string[30];
  char ChiS_string[30];
  char ChiS_COS[30];


  //if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1) angle_final_guide += 180.;
  if(angle_final_guide > 360.) angle_final_guide -= 360.;
  //if(angle_final_guide < 0.) angle_final_guide += 360;

  //sprintf(Angle_GoodLD_string,"#phi = %3.2f#circ", angle_final_GoodLG);
  //sprintf(Angle_Gap_Fibers_string,"#phi = %3.2f#circ  (OLD)", angle_final_Gap_Fibers);
  sprintf(Angle_guide_string,"#phi = %3.2f#circ#pm%3.2f", angle_final_guide, Delta_phi_deg);
  //sprintf(ChiS_string,"#chi^{2}/(ndf-2) = %3.2f", ChiS/(ndf-2));
  sprintf(ChiS_string,"#chi^{2}/ndf = %3.2f", Chis_lepton_fit_3/(ndf_lepton_fit_3));

  //TLatex *tex_Angle_GoodLG;
  //TLatex *tex_Angle_Gap_Fibers;
  TLatex *tex_Angle_guide;
  TLatex *tex_ChiS;
  TLatex *tex_ChiS_COS;

  //tex_Angle_GoodLG = new TLatex(-45.,43.,Angle_GoodLD_string);
  //tex_Angle_GoodLG->SetTextSize(0.05);
  //tex_Angle_GoodLG->SetLineWidth(2);

  //tex_Angle_Gap_Fibers = new TLatex(-45.,43.,Angle_Gap_Fibers_string);
  //tex_Angle_Gap_Fibers->SetTextSize(0.05);
  //tex_Angle_Gap_Fibers->SetLineWidth(2);

  tex_Angle_guide = new TLatex(-45.,43.,Angle_guide_string);
  tex_Angle_guide->SetTextSize(0.05);
  tex_Angle_guide->SetLineWidth(2);

  tex_ChiS = new TLatex(-45.,37.,ChiS_string);
  tex_ChiS->SetTextSize(0.05);
  tex_ChiS->SetLineWidth(2);

  cout << "" << endl;
  cout << "" << endl;
  cout << "Final Angle (deg.) = " << angle_final_guide << " +- " << Delta_phi_deg << "  ( GUIDE )" << endl;
  cout << "Chi Square = " << Chis_lepton_fit_3 << endl;
  cout << "NdF = " << ndf_lepton_fit_3 << endl;
  cout << "Reduced Chi Square = " << Chis_lepton_fit_3/ndf_lepton_fit_3 << endl;
  cout << "" << endl;
  cout << "" << endl;

  cout << "Kaon stop :  " << "X_Kstop = " << X_BAR << "   " << "Y_Kstop = " << Y_BAR << "   " << "Bar = " << i_kaon_bar << endl;

  if(dist_to_k_stop != -10000){
    cout << "Distance from K-stop to best fit line = " << dist_to_k_stop << endl << endl;
  }
  else{
    cout << "No k-stop centroid." << endl << endl;
  }

  if((gap_to_fit_rotate==6 || gap_to_fit_rotate==12 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
    cout << "Intersect with SFT Layer 1:  " << "x = " << -y_TDC_Gap_Fibers_SFT_intersect1 << "   " << "y = " << x_TDC_Gap_Fibers_SFT_intersect1 << endl;
    cout << "Intersect with TOF1:  " << "x = " << -y_TDC_Gap_Fibers_intersect1 << "   " << "y = " << x_TDC_Gap_Fibers_intersect1 << endl;
  }
  else{
    cout << "Intersect with SFT Layer 1:  " << "x = " << x_TDC_Gap_Fibers_SFT_intersect1 << "   " << "y = " << y_TDC_Gap_Fibers_SFT_intersect1 << endl;
    cout << "Intersect with TOF1:  " << "x = " << x_TDC_Gap_Fibers_intersect1 << "   " << "y = " << y_TDC_Gap_Fibers_intersect1 << endl;
  }
  double Yth_TOF1=999.;
  double New_ChiS=999.;
  double ChiS_cos=999.;
  double a_fit_lepton_rotate = 999.;
  double b_fit_lepton_rotate = 999.;
  cout << "" << endl;


  a_fit_lepton_rotate = -1/a_lepton_fit_3;
  b_fit_lepton_rotate = b_lepton_fit_3/a_lepton_fit_3;

  TLine *best_fit_rotate = new TLine(50*a_lepton_fit_3 + b_lepton_fit_3,-50,-50*a_lepton_fit_3 + b_lepton_fit_3,50); //ROTATE_CHANGE



  if((gap_to_fit_rotate==6 || gap_to_fit_rotate==12 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
    cout << "X_TOF1 = " << Gap[gap_to_fit-1][selected_TDC_selected][0] << "   ;   Y_TOF1 = " << Gap[gap_to_fit-1][selected_TDC_selected][1] << endl;
    Yth_TOF1 = a_lepton_fit_3*Gap[gap_to_fit-1][selected_TDC_selected][0] + b_lepton_fit_3;
    cout << "Calculated Y_TOF1 = " << Yth_TOF1 << endl;
    cout << "Delta_Y = " << Gap[gap_to_fit-1][selected_TDC_selected][1] - Yth_TOF1 << endl;
    //New_ChiS = ChiS - 3*((Gap[gap_to_fit-1][selected_TDC_selected][1] - Yth_TOF1)*(Gap[gap_to_fit-1][selected_TDC_selected][1] - Yth_TOF1));
    //New_ChiS = ChiS - 3*((Gap[gap_to_fit-1][selected_TDC_selected][1] - Yth_TOF1)*(Gap[gap_to_fit-1][selected_TDC_selected][1] - Yth_TOF1));
    //ChiS_cos = (New_ChiS/(ndf-3))*(pow(cos(angle_final_guide*M_PI/180.),2));
    //cout << "New ChiS = " << New_ChiS << endl;
    //cout << "NDF' = NDF - 3*TOF1 = " << ndf-3 << endl;
    //cout << "New reduced ChiS = " << New_ChiS/(ndf-3) << endl;
  }
  else{
    cout << "X_TOF1 = " << Gap[gap_to_fit-1][selected_TDC_selected][0] << "   ;   Y_TOF1 = " << Gap[gap_to_fit-1][selected_TDC_selected][1] << endl;
    Yth_TOF1 = a_lepton_fit_3*Gap[gap_to_fit-1][selected_TDC_selected][0] + b_lepton_fit_3;
    cout << "Calculated Y_TOF1 = " << Yth_TOF1 << endl;
    cout << "Delta_Y = " << Gap[gap_to_fit-1][selected_TDC_selected][1] - Yth_TOF1 << endl;
    //New_ChiS = ChiS - 3*((Gap[gap_to_fit-1][selected_TDC_selected][1] - Yth_TOF1)*(Gap[gap_to_fit-1][selected_TDC_selected][1] - Yth_TOF1));
    //ChiS_cos = (New_ChiS/(ndf-3))*(pow(cos(angle_final_guide*M_PI/180.),2));
    //cout << "New ChiS = " << New_ChiS << endl;
    //cout << "NDF' = NDF - 3*TOF1 = " << ndf-3 << endl;
    //cout << "New reduced ChiS = " << New_ChiS/(ndf-3) << endl;
  }
  cout << " " << endl;



  //sprintf(ChiS_COS,"#chi^{2}cos^{2}#phi = %3.2f", ChiS_cos);
  //tex_ChiS_COS = new TLatex(-45.,31.,ChiS_COS);
  //tex_ChiS_COS->SetTextSize(0.05);
  //tex_ChiS_COS->SetLineWidth(2);


  best_fit_rotate->SetLineWidth(2); //ROTATE_CHANGE
  best_fit_rotate->SetLineColor(kRed); // ROTATE_CHANGE
  h_Centroid->Fill(X_BAR,Y_BAR);

  h_int_TDC_TARGET_rotate->Fill(x_exit_rotate,y_exit_rotate);

  TGaxis *A1 = new TGaxis(-50,50,50,50,"pol1",510,"-U");
  TGaxis *A2 = new TGaxis(50,-50,50,50,"pol1",510,"+U");

  TLatex *x_sft;
  TLatex *y_sft;
  TLatex *x_tof1;
  TLatex *y_tof1;

  char X_SFT_String[30];
  char Y_SFT_String[30];
  char X_TOF1_String[30];
  char Y_TOF1_String[30];

  if((gap_to_fit_rotate==6 || gap_to_fit_rotate==12 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
    sprintf(X_SFT_String,"X(40) = %4.2f", -y_TDC_Gap_Fibers_SFT_intersect1);
    sprintf(Y_SFT_String,"Y(40) = %4.2f", x_TDC_Gap_Fibers_SFT_intersect1);
    sprintf(X_TOF1_String,"X(47.1) = %4.2f", -y_TDC_Gap_Fibers_intersect1);
    sprintf(Y_TOF1_String,"Y(47.1) = %4.2f", x_TDC_Gap_Fibers_intersect1);
  }
  else{
    sprintf(X_SFT_String,"X(40) = %4.2f", x_TDC_Gap_Fibers_SFT_intersect1);
    sprintf(Y_SFT_String,"Y(40) = %4.2f", y_TDC_Gap_Fibers_SFT_intersect1);
    sprintf(X_TOF1_String,"X(47.1) = %4.2f", x_TDC_Gap_Fibers_intersect1);
    sprintf(Y_TOF1_String,"Y(47.1) = %4.2f", y_TDC_Gap_Fibers_intersect1);
  }

  x_sft = new TLatex(10.,42.,X_SFT_String);
  x_sft->SetTextSize(0.05);
  x_sft->SetLineWidth(2);

  y_sft = new TLatex(10.,37.,Y_SFT_String);
  y_sft->SetTextSize(0.05);
  y_sft->SetLineWidth(2);

  x_tof1 = new TLatex(-45.,-40.,X_TOF1_String);
  x_tof1->SetTextSize(0.05);
  x_tof1->SetLineWidth(2);

  y_tof1 = new TLatex(-45.,-45.,Y_TOF1_String);
  y_tof1->SetTextSize(0.05);
  y_tof1->SetLineWidth(2);

  // K Stop line / length in target
  TLine *k_stop_line;

  if(X_BAR != -10000){
    if(angle_final_guide >= 90 && angle_final_guide <= 270){ // second x point in the left
      k_stop_line = new TLine(X_BAR, Y_BAR, -50, tan(angle_final_guide*M_PI/180.0)*(-50-X_BAR) + Y_BAR);
    }
    else{ // second x point in the right
      k_stop_line = new TLine(X_BAR, Y_BAR, 50, tan(angle_final_guide*M_PI/180.0)*(50-X_BAR) + Y_BAR);
    }
  }
  else{
    k_stop_line = new TLine(0,0,0,0);
  }

  k_stop_line->SetLineColor(3);
  k_stop_line->SetLineWidth(2);

  double length_in_target = 0;

  if(X_BAR != -10000){
    double x_tof1_intersect_1 = 0;
    double y_tof1_intersect_1 = 0;
    double x_tof1_intersect_2 = 0;
    double y_tof1_intersect_2 = 0;
    double x_tof1_intersect = 0;
    double y_tof1_intersect = 0;

    double alpha = angle_final_guide - 90.0;

    if((gap_to_fit_rotate==6 || gap_to_fit_rotate==12 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
      x_tof1_intersect_1 = y_TDC_Gap_Fibers_intersect1 - dist_to_k_stop*cos(alpha*M_PI/180.0);
      x_tof1_intersect_2 = y_TDC_Gap_Fibers_intersect1 + dist_to_k_stop*cos(alpha*M_PI/180.0);
      y_tof1_intersect_1 = -x_TDC_Gap_Fibers_intersect1 - dist_to_k_stop*sin(alpha*M_PI/180.0);
      y_tof1_intersect_2 = -x_TDC_Gap_Fibers_intersect1 + dist_to_k_stop*sin(alpha*M_PI/180.0);

      if(distance(X_BAR,Y_BAR,x_tof1_intersect_1,y_tof1_intersect_1) < distance(X_BAR,Y_BAR,x_tof1_intersect_2,y_tof1_intersect_2)){
        x_tof1_intersect = x_tof1_intersect_1;
        y_tof1_intersect = y_tof1_intersect_1;
      }
      else{
        x_tof1_intersect = x_tof1_intersect_2;
        y_tof1_intersect = y_tof1_intersect_2;
      }
    }
    else{
      x_tof1_intersect_1 = x_TDC_Gap_Fibers_intersect1 - dist_to_k_stop*cos(alpha*M_PI/180.0);
      x_tof1_intersect_2 = x_TDC_Gap_Fibers_intersect1 + dist_to_k_stop*cos(alpha*M_PI/180.0);
      y_tof1_intersect_1 = y_TDC_Gap_Fibers_intersect1 - dist_to_k_stop*sin(alpha*M_PI/180.0);
      y_tof1_intersect_2 = y_TDC_Gap_Fibers_intersect1 + dist_to_k_stop*sin(alpha*M_PI/180.0);

      if(distance(X_BAR,Y_BAR,x_tof1_intersect_1,y_tof1_intersect_1) < distance(X_BAR,Y_BAR,x_tof1_intersect_2,y_tof1_intersect_2)){
        x_tof1_intersect = x_tof1_intersect_1;
        y_tof1_intersect = y_tof1_intersect_1;
      }
      else{
        x_tof1_intersect = x_tof1_intersect_2;
        y_tof1_intersect = y_tof1_intersect_2;
      }
    }

    length_in_target = distance(X_BAR,Y_BAR,x_tof1_intersect,y_tof1_intersect); //CORRECT
    printf("Length of track in target (xy plane) = %5.2f\n", length_in_target);
  }
  else{
    cout << "No K-Stop for track length." << endl;
  }

//  TGraph *gr_kaon = new TGraph(vec_xx_kaon.size(),&vec_xx_kaon[0],&vec_yy_kaon[0]);
//  gr_kaon->SetMarkerStyle(21);
//  gr_kaon->SetMarkerColor(kBlue-6);
//  gr_kaon->SetMarkerSize(0.8);
//  gr_kaon->GetXaxis()->SetLimits(-50.,50.);
//  gr_kaon->GetYaxis()->SetRangeUser(-50.,50.);

//  TF1 *gr_kaon_fit;
//  double a_fit_kaon = 0.;
//  double b_fit_kaon = 0.;

//  if(vec_xx_kaon.size()>0){
//    gr_kaon->Fit("pol1","QCM");
//    gr_kaon_fit = gr_kaon->GetFunction("pol1");
//    gr_kaon_fit->SetLineWidth(2);
//    gr_kaon_fit->SetLineColor(4);
//    a_fit_kaon = gr_kaon_fit->GetParameter(1);
//    b_fit_kaon = gr_kaon_fit->GetParameter(0);
//  }

  vector<double> vec_test;

/*  if((gap_to_fit_rotate==6 || gap_to_fit_rotate==12 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
    a_lepton_fit_3 = a_fit_lepton_rotate;
  }
  else{
    a_lepton_fit_3 = a_lepton_fit_3;
  }

*/
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Insert here !



  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Calculation of XBAR and YBAR

  //double X_BAR = 0.;
  //double Y_BAR = 0.;
  double X_weights = 0.;
  double Y_weights = 0.;
  double total_energy = 0.;
  vector<double> vec_kaon_centroid_coordinates;
  vector<double> vec_fit_lines_intersect;
  vector<double> vec_k_stop_coordinates;
  vec_kaon_centroid_coordinates.clear();
  vec_fit_lines_intersect.clear();
  vec_k_stop_coordinates.clear();

  for(unsigned int k=0; k<kaon_bars.size(); k++){
    X_weights += ADC_Low_TARGET[kaon_bars[k]]*Xloc[kaon_bars[k]];
    Y_weights += ADC_Low_TARGET[kaon_bars[k]]*Yloc[kaon_bars[k]];
    total_energy += ADC_Low_TARGET[kaon_bars[k]];
  }

  if(kaon_bars.size()>0 || total_energy!=0){
    vec_kaon_centroid_coordinates.push_back(X_weights/total_energy);
    vec_kaon_centroid_coordinates.push_back(Y_weights/total_energy);
  }
  else{
    vec_kaon_centroid_coordinates.push_back(999.99);
    vec_kaon_centroid_coordinates.push_back(999.99);
  }

  if(kaon_bars.size()>0){
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


  if(kaon_bars.size() < 5){
    vec_k_stop_coordinates.push_back(vec_kaon_centroid_coordinates[0]);
    vec_k_stop_coordinates.push_back(vec_kaon_centroid_coordinates[1]);
  }
  else{
    vec_k_stop_coordinates.push_back(vec_fit_lines_intersect[0]);
    vec_k_stop_coordinates.push_back(vec_fit_lines_intersect[1]);


/*    if(gap_to_fit==6 || gap_to_fit==12){
      vec_k_stop_coordinates = _2lines_intersect(a_fit_lepton_rotate, b_fit_lepton_rotate, a_fit_kaon, b_fit_kaon);
    }
    else{
      k_stop_coordinates = _2lines_intersect(a_lepton_fit_3, b_lepton_fit_3, a_fit_kaon, b_fit_kaon);
    }
*/}

    h_kaon_stop->Fill(vec_k_stop_coordinates[0],vec_k_stop_coordinates[1]);
    h_kaon_stop->SetMarkerStyle(20);
    h_kaon_stop->SetMarkerColor(1);
    h_kaon_stop->SetMarkerSize(1);



  h_Centroid->SetMarkerStyle(5);         h_2lines_intersect->SetMarkerStyle(20);
  h_Centroid->SetMarkerColor(1);         h_2lines_intersect->SetMarkerColor(1);
  h_Centroid->SetMarkerSize(4);          h_2lines_intersect->SetMarkerSize(0.8);



  //if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){ //ROTATE_CHANGE

  if (graphics == 1) {
    char ch_centroid_x[200];
    char ch_centroid_y[200];

    if(X_BAR != -10000){
      sprintf(ch_centroid_x,"#bar{x} = %5.2f",X_BAR);
      sprintf(ch_centroid_y,"#bar{y} = %5.2f",Y_BAR);
    }
    else{
      sprintf(ch_centroid_x,"#bar{x} = NULL");
      sprintf(ch_centroid_y,"#bar{y} = NULL");
    }

    char ch_Kbar[100];
    sprintf(ch_Kbar, "K_{stop} Bar = %d",i_kaon_bar);

    TLatex *tex_Kbar = new TLatex(0.64,0.9295171,ch_Kbar);
    tex_Kbar->SetTextSize(0.05);
    tex_Kbar->SetLineWidth(2);
    tex_Kbar->Draw();

    TLatex *tex_Label_Centroid_x = new TLatex(22.,-40.,ch_centroid_x);
    tex_Label_Centroid_x->SetTextSize(0.05);
    tex_Label_Centroid_x->SetLineWidth(2);

    TLatex *tex_Label_Centroid_y = new TLatex(22.,-45.,ch_centroid_y);
    tex_Label_Centroid_y->SetTextSize(0.05);
    tex_Label_Centroid_y->SetLineWidth(2);

    c2->cd(6)->Update();
    c2->cd(6)->Range(-50, -50, 50, 50);
    vec_test = Subplot_6(h_EMPTY, A1,
      A2, h_Circle,
      ell, ell_Target,
      ell_L1, h_Target_Center, gap_to_fit,
      h_TOF1, h_TOF1_closest,
      h_int_TDC_Gap_Fibers_rotate, h_int_TDC_Gap_Fibers,
      gr3_Leptons_rotate, gr_lepton_3,
      best_fit_rotate, h_TOF1_rotate,
      h_int_TDC_Gap_Fibers_SFT,
      h_int_TDC_Gap_Fibers_SFT_rotate,
      x_guide, y_guide,
      x_guide_rotate, y_guide_rotate,
      tex_Angle_guide, tex_ChiS,
      x_sft, y_sft, x_tof1, y_tof1,
      tex_Label_Centroid_x, tex_Label_Centroid_y,
      a_fit_lepton_rotate,  b_fit_lepton_rotate,
      a_lepton_fit_3,  b_lepton_fit_3,
      a_fit_kaon, b_fit_kaon,
      gr_kaon, gr_kaon_bk,
      h_kaon_stop,
      vec_xx_kaon.size(), kaon_bk);

    c2->cd(4)->Range(-50, -50, 50, 50);
    Subplot_4(h_EMPTY, A1,
      A2, h_Circle,
      ell, ell_Target,
      ell_L1, h_Target_Center, gap_to_fit,
      h_TOF1, h_TOF1_closest,
      h_int_TDC_Gap_Fibers_rotate, h_int_TDC_Gap_Fibers,
      gr3_Leptons_rotate, gr_lepton_3,
      best_fit_rotate);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////]
    c2->cd(5)->Range(-50, -50, 50, 50);
    Subplot_5(h_EMPTY, A1,
      A2, h_Circle,
      ell, ell_Target,
      ell_L1, h_Target_Center, gap_to_fit,
      h_TOF1, h_TOF1_closest,
      h_int_TDC_Gap_Fibers_rotate, h_int_TDC_Gap_Fibers,
      gr_kaon, gr_kaon_bk,
      vec_xx_kaon.size(), kaon_bk);

    TCanvas *c_TEST;
    c_TEST = new TCanvas("TEST","TEST",0,200,800,800);
    c_TEST->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
    c_TEST->Divide(2,2);

    c_TEST->cd(1)->Range(-50, -50, 50, 50);
    TEST_subplot_1(h_EMPTY, A1,
      A2, h_Circle,
      ell, ell_Target,
      ell_L1, h_Target_Center,
      gr_lepton_1);

    c_TEST->cd(2)->Range(-50, -50, 50, 50);
    TEST_subplot_2(h_EMPTY, A1,
      A2, h_Circle,
      ell, ell_Target,
      ell_L1, h_Target_Center,
      vec_xx_kaon.size(), kaon_bk,
      gr_kaon_bk, gr_kaon);


    c_TEST->cd(3)->Range(-50, -50, 50, 50);
    TEST_subplot_1(h_EMPTY, A1,
      A2, h_Circle,
      ell, ell_Target,
      ell_L1, h_Target_Center,
      gr_lepton_2);



    c_TEST->cd(4)->Range(-50, -50, 50, 50);
    TEST_subplot_1(h_EMPTY, A1,
      A2, h_Circle,
      ell, ell_Target,
      ell_L1, h_Target_Center,
      gr_lepton_3);
  }
  else {
    if (graphics == 0) {
      if(gap_to_fit==6 || gap_to_fit==12){ //ROTATE_CHANGE
        vec_test = _2lines_intersect(a_fit_lepton_rotate, b_fit_lepton_rotate, a_fit_kaon, b_fit_kaon);
      }
      else{
          vec_test = _2lines_intersect(a_lepton_fit_3, b_lepton_fit_3, a_fit_kaon, b_fit_kaon);
      }

      // delete c2;
      c2->Close();
      // gSystem->ProcessEvents();
    }
  }
  // c2->Close();









  vector<double> ZZ = SFT_Test_CR_Event2(Run_Number, ievt, angle_final_guide, C2X_centroid, length_in_target);
  vector<double> Z_selected;
  Z_selected.empty();

  /// TOF1_Z_range
  Z_TOF1[gap_to_fit_rotate-1] =
  Slope_correction_TOF1[gap_to_fit_rotate-1]*((0.025*150*0.5*(TDC_TOF1U[gap_to_fit_rotate-1] - TDC_TOF1D[gap_to_fit_rotate-1])) + Offset_TOF1[gap_to_fit_rotate-1]);

  cout << "Z_TOF1  =  " << Z_TOF1[gap_to_fit_rotate-1] << " mm" << endl;

  for(unsigned i=0; i<ZZ.size(); i++){
      Z_selected.push_back(ZZ[i]);
  }

  /// SORT VECTOR ELEMENTS
  vector<double> Z_selected_sorted;
  Z_selected_sorted.empty();

  sort(Z_selected.begin(), Z_selected.end());

  cout << "Z selected:  ";
  for(unsigned k=0; k<Z_selected.size(); k++){
    cout << Z_selected[k] << "  ";
  }

  cout << endl;
  cout << endl;

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

  cout << endl;
  cout << endl;
  cout << "////////////////////////////////////////////////////////////" << endl;
  cout << endl;


  cout << fixed;
  cout << "DEBUG : " << "  ";
  cout << setw(4) << Run_Number << "  ";
  cout << setw(7) << ievt << "  ";
  ////cout << setw(7) << Read_Run_Number << "  ";
  cout << setw(2) << gap_to_fit_rotate << "  ";
  cout << setw(2) << selected_TOF2 << "  ";
  cout << setw(7) << setprecision(3) << angle_final_guide << "  ";
  cout << setw(5) << setprecision(2) << Delta_phi_deg << "  ";
  cout << setw(8) << setprecision(2) << Chis_lepton_fit_3 << "  ";
  cout << setw(3) << setprecision(0) << ndf_lepton_fit_3 << "  ";
  cout << setw(7) << setprecision(2) << Chis_lepton_fit_3/ndf_lepton_fit_3 << "    ||  ";
  //cout << setw(8) << setprecision(2) << New_ChiS << "  ";
  //cout << setw(3) << setprecision(0) << ndf-3 << "  ";
  //cout << setw(7) << setprecision(2) << New_ChiS/(ndf-3) << "  ";
  //cout << setw(7) << setprecision(2) << (New_ChiS/(ndf-3))*(pow(cos(angle_final_guide*M_PI/180.),2)) << "  ";
  cout << setw(6) << setprecision(2) << vec_kaon_centroid_coordinates[0] << "  ";
  cout << setw(6) << setprecision(2) << vec_kaon_centroid_coordinates[1] << "   ||  ";

  cout << setw(6) << setprecision(2) << vec_fit_lines_intersect[0] << "  ";
  cout << setw(6) << setprecision(2) << vec_fit_lines_intersect[1] << "   ||  ";

  cout << setw(6) << setprecision(2) << vec_k_stop_coordinates[0] << "  ";
  cout << setw(6) << setprecision(2) << vec_k_stop_coordinates[1] << "    ||  ";

  //cout << setw(6) << setprecision(2) << X_BAR << "  ";
  //cout << setw(6) << setprecision(2) << Y_BAR << "  ";
  cout << setw(3) << i_kaon_bar << "  ";
  cout << setw(2) << vec_xx_kaon.size() << "  ";
  ////cout << setw(2) << K_stop_count << "  ";
  cout << setw(7) << setprecision(2) << x_TARGET_intersect << "  ";
  cout << setw(7) << setprecision(2) << y_TARGET_intersect << "  ";
  cout << setw(3) << vec_Ck.size() << "  ";
  cout << setw(3) << vec_Cpi.size() << "  ";
  cout << setw(8) << setprecision(3) << length_in_target << "  ";
  cout << setw(8) << setprecision(3) << C2X_centroid << "  ";
  cout << setw(6) << setprecision(1) << TDC_diff << "  ";
  cout << setw(6) << setprecision(2) << vec_test[0] << "  ";
  cout << setw(6) << setprecision(2) << vec_test[1] << "  ";

  //cout << setw(8) << setprecision(1) << 550 - (0.91*TDC_average - 0.625*tdc_Ck) << "  ";

  cout << endl;


  cout << endl;

  ////gr_Leptons->Delete();
  ////gr_Leptons_fit3->Delete();
  ////gr_Leptons_fit->Delete();
  //gr_Leptons_fit2->Delete();
  //gr2_Leptons->Delete();
  ////gr3_Leptons->Delete();
  //gr_kaon_fit->Delete();
  //gr_Leptons_fit2->Delete();
  ////gr_Leptons_fit3->Delete();
  ////gr_test_lepton->Delete();
  ////gr_test_kaon->Delete();

  /*
  delete gr_kaon_fit;           gr_kaon_fit = nullptr;
  delete gr_test_lepton_fit;    gr_test_lepton_fit = nullptr;
  delete gr_test_kaon_fit;      gr_test_kaon_fit = nullptr;
  delete gr_Leptons_fit;        gr_Leptons_fit = nullptr;
  //delete gr_Leptons_fit2;       gr_Leptons_fit2 = nullptr;
  delete gr_Leptons_fit3;       gr_Leptons_fit3 = nullptr;


  delete gr_kaon;               gr_kaon = nullptr;
  delete gr_test_lepton;        gr_test_lepton = nullptr;
  delete gr_test_kaon;          gr_test_kaon = nullptr;
  delete gr_Leptons;            gr_Leptons = nullptr;
  //delete gr2_Leptons;           gr2_Leptons = nullptr;
  delete gr2_Leptons_rotate;    gr2_Leptons_rotate = nullptr;
  delete gr3_Leptons;           gr3_Leptons = nullptr;
  delete gr3_Leptons_test;      gr3_Leptons_test = nullptr;
  delete gr3_Leptons_rotate;    gr3_Leptons_rotate = nullptr;
  */

  return;


} // End void


double SFT_Test_CR_Event(int Run_Number, int evt, double phi, bool to_print, double C2X_centroid, double length_in_target){
  Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];
  Int_t adc_low_sft[128];
  Int_t tdc_le_sft[128][16];
  Int_t tdc_te_sft[128][16];

  Int_t HG_SFT_ADC_Thr[128] = {0};

  for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;

  Double_t ADC_High_SFT_corr[128];
  Int_t has_TDC_SFT_hit[128] = {0};

  char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!

  char file_mapping[200];
  sprintf(file_mapping,"../Mapping");

  char par_finput[200];
  sprintf(par_finput,"%s/%s",file_mapping,source_mapping);

  Int_t par_temp[2][128];
  ifstream fdat(par_finput,ios::in);
  for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
  fdat.close();

  char path_input[200];
  sprintf(path_input,"%s",path_merged);

  char Name_finput[200];
  sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);

  TChain *fChain= new TChain("Tree");
  fChain->Add(Name_finput);
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
  fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
  fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
  fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

  fChain->GetEntry(evt);

  for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
    ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-HG_SFT_ADC_Thr[j_SFT];
  }

  for(int j=0 ; j<128 ; j++){
    if(ADC_High_SFT[j]<0) ADC_High_SFT_corr[j]=0;
    else ADC_High_SFT_corr[j]=ADC_High_SFT[j];
  }

  for(Int_t ii=0; ii<128; ii++){
    for (Int_t qq=0; qq<6; qq++) {
      if (tdc_le_sft[ii][qq] > SFT_TDC_min[ii] && tdc_le_sft[ii][qq] < SFT_TDC_max[ii]) has_TDC_SFT_hit[ii]++;
    }
  }


  double sft_z_selected = 0.;
  sft_z_selected = SFT_print(ADC_High_SFT_corr, has_TDC_SFT_hit, SFT_channel_to_fiber, evt, phi, to_print, 0, false, C2X_centroid, length_in_target);

  return sft_z_selected;
}

vector<double> SFT_Test_CR_Event2(int Run_Number, int evt, double phi, double C2X_centroid, double length_in_target){
  Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];
  Int_t adc_low_sft[128];
  Int_t tdc_le_sft[128][16];
  Int_t tdc_te_sft[128][16];

  Int_t HG_SFT_ADC_Thr[128] = {0};


  for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;

  Double_t ADC_High_SFT_corr[128];
  Int_t has_TDC_SFT_hit[128] = {0};

  char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!

  char file_mapping[200];
  sprintf(file_mapping,"../Mapping");

  char par_finput[200];
  sprintf(par_finput,"%s/%s",file_mapping,source_mapping);

  Int_t par_temp[2][128];
  ifstream fdat(par_finput,ios::in);
  for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
  fdat.close();

  char path_input[200];
  sprintf(path_input,"%s",path_merged);

  char Name_finput[200];
  sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);


  TChain *fChain= new TChain("Tree");
  fChain->Add(Name_finput);
  fChain->SetMakeClass(1);
  fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
  fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
  fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
  fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

  fChain->GetEntry(evt);

  for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
    ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-HG_SFT_ADC_Thr[j_SFT];
  }

  for(int j=0 ; j<128 ; j++){
    if(ADC_High_SFT[j]<0) ADC_High_SFT_corr[j]=0;
    else ADC_High_SFT_corr[j]=ADC_High_SFT[j];
  }

  for(Int_t ii=0; ii<128; ii++){
    for(Int_t qq=0; qq<6; qq++) {
      if(tdc_le_sft[ii][qq] > SFT_TDC_min[ii] && tdc_le_sft[ii][qq] < SFT_TDC_max[ii]) has_TDC_SFT_hit[ii]++;
    }
  }


  //double sft_z_selected = 0.;
  vector<double> ZZ;
  ZZ = Z_Avg(ADC_High_SFT_corr, has_TDC_SFT_hit, SFT_channel_to_fiber, phi, true, 0, false, C2X_centroid, length_in_target);

  //return sft_z_selected;
  return ZZ;
}
