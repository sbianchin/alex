#ifndef __CINT__
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>
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
#include "ADC_Thresholds_CR.h"
#include "TDC_Windows_CR.h"
#include "Cuts_and_Windows_CR.h"
#include "MWPC_Thr_CR.h"
#endif

#include "intersect.cxx"
#include "C2_Strip_transform.h"
#include "Channel_to_Strip.h"
#include "SFT_functions.h"
//#include "SFT_functions_copy.h"



using namespace std;

static double C2_intersect_SFT_top=999.99;
static double C2_intersect_SFT_bottom=999.99;
static int ana_flag=1;


void SFT_Test_CR_Event(int Run_Number, int evt, double phi, int gap_to_fit);
void SFT_angle_selector(int Run_Number, int ievt, double phi, int gap_to_fit_left, int gap_to_fit_right);

void Event_Display_Cosmic_Rays(Int_t Run_Number=5, Int_t ievt=0, Int_t enable_cout=9, Int_t display_fits = 0){ 
  
  int Switch=1; // Displays hit with no HG, but LG (0 = OFF ; 1 = ON)
  int Rotate=1; // When TOF1 is 12 or 6, or when TOF2 is 6 or 12, rotate by -90 deg to fit a horizontal line (0 = OFF ; 1 = ON)
  int MWPC_switch = 1; // Display canvas wire chamber hits (0 = OFF ; 1 = ON)
  
  // Event blacklist.
  ifstream blacklist;
  blacklist.open("Event_Blacklist_Test.txt");
  

  if(enable_cout!=0 && enable_cout!=1 && enable_cout!=2 && enable_cout!=9){
    cout << "  " << endl;
    cout << "Flag Error !" << endl;
    cout << "  " << endl;
    return;
  }

  int T_limit = 3;

  gStyle->Clear();
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(11);
  
  char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!

  Double_t flag_size_TARGET=1.35;
  Double_t flag_size_SFT=1.3;
  Double_t flag_size_palette=1.6;

  Int_t HG_TARGET_ADC_Thr[256] = {0};
  Int_t LG_TARGET_ADC_Thr[256] = {0};
  Int_t HG_SFT_ADC_Thr[128] = {0};
  Int_t LG_SFT_ADC_Thr[128] = {0};

  
  for(int i=0; i<256; i++)  HG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Thr_HG[i]) + TARGET_ADC_Thr_HG_Offset;
  for(int i=0; i<256; i++)  LG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Thr_LG[i]) + TARGET_ADC_Thr_LG_Offset;
  for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;
  for(int i=0; i<128; i++)  LG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_LG[i]) + SFT_ADC_Thr_LG_Offset;
  
 

  Int_t TDC_min_TARGET = TARGET_TDC_min[0];
  Int_t TDC_max_TARGET = TARGET_TDC_max[0];
  Int_t ADC_TARGET_Thr = HG_TARGET_ADC_Thr[0];
  Int_t TDC_min_SFT = SFT_TDC_min[0];
  Int_t TDC_max_SFT = SFT_TDC_max[0];
  
  bool TARGET_High_has_neighbours[256] = {false}; // Array of High gain target hits which have no neighbouring targets hit   
  

  int n_hit = 2;   // ## Number of hits required in the TARGET

  Int_t ADC_High_corr_max=0;

  Int_t ADC_cut_TARGET2 = 850;

  char ADC_cut[100];    sprintf(ADC_cut,"(ADC >= %d)",HG_SFT_ADC_Thr[0]);

  Int_t adc_high_target[256];       Int_t ADC_High_TARGET[256]={0};    
  Int_t adc_low_target[256];        Int_t ADC_Low_TARGET[256];
  Int_t ADC_Low_TARGET_ped[256];    Int_t ADC_High_TARGET_ped[256];
  Int_t tdc_le_target[256][16];     Int_t TDC_LE_TARGET[256];     
  Int_t tdc_te_target[256][16];     //Int_t TDC_TE_TARGET[256];  
  Int_t TDC_LE_TARGET_corrected[256][6]={0};
  char TDC_LE_TARGET_corr[256][6][20];

  Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];        Double_t ADC_High_SFT_corr[128];    
  Int_t adc_low_sft[128];           Int_t ADC_Low_SFT[128];   
  Int_t tdc_le_sft[128][16];        Int_t TDC_LE_SFT[128];         
  Int_t tdc_te_sft[128][16];        //Int_t TDC_TE_SFT[128];  

  Int_t ADC_TOF1[24];   
  //Int_t ADC_TOF2[48];  // !!!
   
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

  Int_t MwpcADC[512];               Int_t MWPCADC[512];
  Int_t adc_c2x_r[56];              Int_t ADC_C2X_R[56];
  Int_t adc_c2x_l[56];              Int_t ADC_C2X_L[56];
  Int_t adc_c2y_r[16];              Int_t ADC_C2Y_R[16];
  Int_t adc_c2y_l[16];              Int_t ADC_C2Y_L[16];

  Int_t adc_c3x_r[64];              Int_t ADC_C3X_R[64];
  Int_t adc_c3x_l[64];              Int_t ADC_C3X_L[64];
  Int_t adc_c3y_r[16];              Int_t ADC_C3Y_R[16];
  Int_t adc_c3y_l[16];              Int_t ADC_C3Y_L[16];

  Int_t adc_c4x_r[72];              Int_t ADC_C4X_R[72];
  Int_t adc_c4x_l[72];              Int_t ADC_C4X_L[72];
  Int_t adc_c4y_r[16];              Int_t ADC_C4Y_R[16];
  Int_t adc_c4y_l[16];              Int_t ADC_C4Y_L[16];

  Int_t tdc_Ck[14][16];
  Int_t tdc_Cpi[14][16];

  Int_t has_TDC_SFT_hit[128] = {0};

  float R_TARGET = 29.0;
  float R_TOF1 = 47.1;
  float R_SFT_L1 = 40.0;


  char run_string[100];            char event_string[100];

  char h_ADC_title[200];



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

 /////////////////////////////////////////////////////////////////////////////////////////////////////

  char h_target_ADC_title[100];
  char h_target_TDC_title[100];
  //  sprintf(h_target_ADC_title,"Run %d (Event %d) -- ADC #geq %d)",Run_Number, ievt, par_temp_TARGET[1][0]);
  sprintf(h_target_ADC_title,"Run %d (Event %d) -- ADC #geq %d)",Run_Number, ievt, ADC_TARGET_Thr);
  //  sprintf(h_target_TDC_title,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt, par_temp_TARGET[1][0],TDC_min_TARGET,TDC_max_TARGET);
  sprintf(h_target_TDC_title,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt, ADC_TARGET_Thr,TDC_min_TARGET,TDC_max_TARGET);

  char h_target_ADC_title2[100];
  char h_target_TDC_title2[100];
  char h_target_ADC_title3[100];
  char h_target_ADC_title4[100];
  sprintf(h_target_ADC_title2,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt, ADC_TARGET_Thr,TDC_min_TARGET, TDC_max_TARGET);
  sprintf(h_target_TDC_title2,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt, ADC_TARGET_Thr, TDC_min_TARGET,TDC_max_TARGET);
  sprintf(h_target_ADC_title3,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt, ADC_TARGET_Thr, TDC_min_TARGET,TDC_max_TARGET);
  sprintf(h_target_ADC_title4,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt, ADC_TARGET_Thr, TDC_min_TARGET,TDC_max_TARGET);
  
  sprintf(h_ADC_title,"(ADC offset = %d) | (%d < TDC < %d)  --  Run %d (Event %d)",SFT_ADC_Thr_HG_Offset, TDC_min_SFT,TDC_max_SFT,Run_Number,ievt);

  char h_target_TDC_copy_Name[200];     char h_target_TDC_copy_Title[200];
  char h_TDC_selected_Name[200];        char h_TDC_selected_Title[200];
  char h_GoodLG_Name[200];              char h_GoodLG_Title[200];
  char h_TDC_selected2_Name[200];       char h_TDC_selected2_Title[200];
  char h_GoodLG_weighted_Name[200];     char h_GoodLG_weighted_Title[200];
  char h_Gap_Fibers_Name[200];          char h_Gap_Fibers_Title[200];

  sprintf(h_target_TDC_copy_Name,"Event %d (Run %d)", ievt, Run_Number);
  sprintf(h_target_TDC_copy_Title,"ADC_HG #geq %d  |  %d #leq TDC #leq %d", ADC_TARGET_Thr, TDC_min_TARGET, TDC_max_TARGET);
 
  sprintf(h_TDC_selected_Name,"Event %d (Run %d)", ievt, Run_Number);
  sprintf(h_TDC_selected_Title,"ADC_HG #geq %d  |  %d #leq TDC #leq %d", ADC_TARGET_Thr+100, TDC_min_TARGET, TDC_max_TARGET);

  sprintf(h_GoodLG_Name,"Event %d (Run %d)", ievt, Run_Number);
  sprintf(h_GoodLG_Title,"ADC_LG #geq %d  |  %d #leq TDC #leq %d", ADC_TARGET_Thr,TDC_min_TARGET,TDC_max_TARGET);

  sprintf(h_TDC_selected2_Name,"Event %d (Run %d)", ievt, Run_Number);
  sprintf(h_TDC_selected2_Title,"ADC_HG #geq %d  |  %d #leq TDC #leq %d  (WEIGHTED)", ADC_TARGET_Thr+100,TDC_min_TARGET,TDC_max_TARGET);

  sprintf(h_GoodLG_weighted_Name,"Event %d (Run %d)", ievt, Run_Number);
  sprintf(h_GoodLG_weighted_Title,"ADC_LG #geq %d  |  %d #leq TDC #leq %d  (WEIGHTED)", ADC_TARGET_Thr,TDC_min_TARGET,TDC_max_TARGET);

  sprintf(h_Gap_Fibers_Name,"Event %d (Run %d)", ievt, Run_Number);
  sprintf(h_Gap_Fibers_Title,"ADC LG offset =  %d  |  %d #leq TDC #leq %d", TARGET_ADC_Thr_LG_Offset,TDC_min_TARGET,TDC_max_TARGET);


  TH1F *h_ADC_L1_DS = new TH1F("h_ADC_L1_DS",h_ADC_title,130,-1,129);
  TH1F *h_ADC_L2_DS = new TH1F("h_ADC_L2_DS",h_ADC_title,130,-1,129);
  TH1F *h_ADC_L3_DS = new TH1F("h_ADC_L3_DS",h_ADC_title,130,-1,129);
  TH1F *h_ADC_L4_DS = new TH1F("h_ADC_L4_DS",h_ADC_title,130,-1,129);

  TH1F *h_ADC_L1_US = new TH1F("h_ADC_L1_US",h_ADC_title,66,63,129);
  TH1F *h_ADC_L2_US = new TH1F("h_ADC_L2_US",h_ADC_title,66,63,129);
  TH1F *h_ADC_L3_US = new TH1F("h_ADC_L3_US",h_ADC_title,66,63,129);
  TH1F *h_ADC_L4_US = new TH1F("h_ADC_L4_US",h_ADC_title,66,63,129);

  TH2F *h_Circle = new TH2F("Test", "TARGET and TOF1", 500, -50, 50, 500, -50, 50);

  TH2F *h_target_ADC = new TH2F("Histo Fit 0",h_target_ADC_title,3000,-50,50,3000,-50,50);
  TH2F *h_target_TDC = new TH2F("Histo Fit 1",h_target_TDC_title,3000,-50,50,3000,-50,50);
  TH2F *h_target_TDC_copy = new TH2F(h_target_TDC_copy_Name,h_target_TDC_copy_Title,3000,-50,50,3000,-50,50);
  TH2F *h_TDC_selected = new TH2F(h_TDC_selected_Name, h_TDC_selected_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_TDC_selected2 = new TH2F(h_TDC_selected2_Name, h_TDC_selected2_Title, 500, -50, 50, 500, -50, 50);

  TH2F *h_target_ADC2 = new TH2F("Histo Fit 2",h_target_ADC_title2,3000,-50,50,3000,-50,50);
  TH2F *h_target_TDC2 = new TH2F("Histo Fit 3",h_target_TDC_title2,3000,-50,50,3000,-50,50);

  TH2F *h_target_ADC3 = new TH2F("Histo Fit 4",h_target_ADC_title3,3000,-50,50,3000,-50,50);
  TH2F *h_target_ADC4 = new TH2F("Histo Fit 5",h_target_ADC_title4,3000,-50,50,3000,-50,50);
  TH2F *h_target_ADCA = new TH2F("Histo Fit 6",h_target_ADC_title3,3000,-50,50,3000,-50,50);

  TH2F *h_Target = new TH2F("Test Target", "Target", 3000, -50, 50, 3000, -50, 50);
  TH2F *h_Target_Center = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);

  TH2F *h_kaon = new TH2F("Kaon", "Kaon", 500, -50, 50, 500, -50, 50);
  TH2F *h_kaon_copy = new TH2F("Kaon Copy", "Kaon Copy", 500, -50, 50, 500, -50, 50);
    
  TH2F *h_max = new TH2F("Max", "Max", 500, -50, 50, 500, -50, 50);
  TH2F *h_max_copy = new TH2F("Max Copy", "Max Copy", 500, -50, 50, 500, -50, 50);

  TH2F *h_GoodLG = new TH2F(h_GoodLG_Name, h_GoodLG_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_GoodLG_weighted = new TH2F(h_GoodLG_weighted_Name, h_GoodLG_weighted_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_GoodLG_copy = new TH2F("Good LG Copy", "Good LG Copy", 500, -50, 50, 500, -50, 50);

  TH2F *h_TOF1 = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_closest = new TH2F("h_TOF1_closest", "h_TOF1_closest", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_rotate = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50); // ROTATE_CHANGE
  TH2F *h_TOF1_closest_rotate = new TH2F("h_TOF1_closest", "h_TOF1_closest", 500, -50, 50, 500, -50, 50); //ROTATE_CHANGE

  TH2F *h_int_TDC = new TH2F("h_int_TDC", "h_int_TDC", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_selected = new TH2F("h_int_TDC_selected", "h_int_TDC_selected", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_selected_weighted = new TH2F("h_int_TDC_selected_weighted", "h_int_TDC_selected_weighted", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_Gap_Fibers = new TH2F("h_int_TDC_Gap_Fibers", "h_int_TDC_Gap_Fibers", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_Gap_Fibers_SFT = new TH2F("h_int_TDC_Gap_Fibers_SFT","h_int_TDC_Gap_Fibers_SFT", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_Gap_Fibers_rotate = new TH2F("h_int_TDC_Gap_Fibers", "h_int_TDC_Gap_Fibers", 500, -50, 50, 500, -50, 50);  // ROTATE_CHANGE
  TH2F *h_int_TDC_Gap_Fibers_SFT_rotate = new TH2F("h_int_TDC_Gap_Fibers_SFT","h_int_TDC_Gap_Fibers_SFT", 500, -50, 50, 500, -50, 50);  // ROTATE_CHANGE
  TH2F *h_int_TDC_SFT = new TH2F("h_int_TDC_SFT", "h_int_TDC_SFT", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_TARGET = new TH2F("h_int_TDC_TARGET", "h_int_TDC_TARGET", 500, -50, 50, 500, -50, 50);


  TH2F *h_int_GoodLG = new TH2F("h_int_GoodLG", "h_int_GoodLG", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_GoodLG_weighted = new TH2F("h_int_GoodLG_weighted", "h_int_GoodLG_weighted", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_GoodLG_SFT = new TH2F("h_int_GoodLG_SFT", "h_int_GoodLG_SFT", 500, -50, 50, 500, -50, 50);
    
  TH2F *h_TDC_Gap_Fibers = new TH2F(h_Gap_Fibers_Name, h_Gap_Fibers_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_K_Stop_Bars = new TH2F(h_Gap_Fibers_Name, h_Gap_Fibers_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_TDC_Gap_Fibers_copy = new TH2F(h_Gap_Fibers_Name, h_Gap_Fibers_Title, 500, -50, 50, 500, -50, 50); //ROTATE_CHANGE
  TH2F *h_TDC_Gap_Fibers_kaon = new TH2F(h_Gap_Fibers_Name, h_Gap_Fibers_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_TARGET_LG_Blue = new TH2F(h_Gap_Fibers_Name, h_Gap_Fibers_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_Centroid = new TH2F("Centroid", "Centroid", 500, -50, 50, 500, -50, 50);
  TEllipse *ell = new TEllipse(0, 0, R_TOF1, 0);
  TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0);
  TEllipse *ell_L1 = new TEllipse(0, 0, R_SFT_L1, 0);
  

  TLine *HorizontalAxis = new TLine(0., 0., 50., 0.);



  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////***mwpc+target drawings***//////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  char h_C2_target_title[200];
  sprintf(h_C2_target_title,"C2 Chambers and Target  |  Run %d  |  Event %d", Run_Number, ievt);

  TH2F *h_C2_target = new TH2F("C2 Chambers and Target",h_C2_target_title,700,-700,700, 700,-700,700);
  TEllipse *ell_C2 = new TEllipse(0,0,629.4,0);



  char path_input[200];                   
  sprintf(path_input,"%s",path_merged);          


  char Name_finput[200];
  sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);

  char footer[100];
  sprintf(footer,"Event_Display_MS.C  --  Run %d ; Event %d",Run_Number,ievt);

  cout << "   " << endl;
  cout << "   " << endl;
  cout << "************************************************************************************************************" << endl;

  //cout << "   " << endl;
  cout << "File opened:  " << Name_finput << endl;
  cout << "SFT Mapping File:  " << par_finput << endl;
  cout << "MWPC Mapping File:  " << par_finput2 << endl;


  sprintf(run_string,"Run %d ; Event %d",Run_Number,ievt);
  sprintf(event_string,"Run %d ; Event %d",Run_Number,ievt);


  TChain *fChain= new TChain("Tree");   
  fChain->Add(Name_finput);   
  fChain->SetMakeClass(1);              

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
 
  Int_t nentries = (Int_t)fChain->GetEntries();
  cout << "Total Number of Events:  " << nentries <<endl;
  //cout << "  " << endl;
  cout << "************************************************************************************************************" << endl;
  cout << "  " << endl;
  cout << "  " << endl;

  bool Good_Event=false;
  bool Good_tof1[12] = {false};
  bool Good_tof2[12] = {false};

  for(int ivt=ievt; ivt<ievt+1; ivt++){
    fChain->GetEntry(ivt);  
    
    for(int init=0; init<256; init++){
      ADC_High_TARGET[init]=0;
      HG_TARGET_ADC_Thr[init] = 0;
      LG_TARGET_ADC_Thr[init] = 0;
    }
    
    for(int init=0; init<128; init++){
      HG_SFT_ADC_Thr[init] = 0;
      LG_SFT_ADC_Thr[init] = 0;
    }
  
    for(int i=0; i<256; i++)  HG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Thr_HG[i]) + TARGET_ADC_Thr_HG_Offset;
    for(int i=0; i<256; i++)  LG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Thr_LG[i]) + TARGET_ADC_Thr_LG_Offset;
    for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;
    for(int i=0; i<128; i++)  LG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_LG[i]) + SFT_ADC_Thr_LG_Offset;

    for (int j_TARGET=0; j_TARGET<256; j_TARGET++){
      ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-HG_TARGET_ADC_Thr[j_TARGET];
      ADC_High_TARGET_ped[j_TARGET]=adc_high_target[j_TARGET]-HG_TARGET_ADC_Thr[j_TARGET] + TARGET_ADC_Thr_HG_Offset;
      ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-LG_TARGET_ADC_Thr[j_TARGET];
      ADC_Low_TARGET_ped[j_TARGET]=adc_low_target[j_TARGET]-LG_TARGET_ADC_Thr[j_TARGET] + TARGET_ADC_Thr_LG_Offset; 
      TDC_LE_TARGET[j_TARGET]=tdc_le_target[j_TARGET][0];
   } 
  
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

    for (Int_t j_MWPC=0; j_MWPC<512; j_MWPC++) {
      MWPCADC[j_MWPC] = MwpcADC[j_MWPC]-MWPC_ADC_Thr[j_MWPC];
    }

    for(Int_t j_C2=0; j_C2<56; j_C2++){
      ADC_C2X_R[j_C2] = adc_c2x_r[j_C2]-ADC_C2X_Thr[j_C2];
      ADC_C2X_L[j_C2] = adc_c2x_l[j_C2]-ADC_C2X_Thr[j_C2];
    }

    for(Int_t j_C3=0; j_C3<64; j_C3++){
      ADC_C3X_R[j_C3] = adc_c3x_r[j_C3]-ADC_C3X_Thr[j_C3];
      ADC_C3X_L[j_C3] = adc_c3x_l[j_C3]-ADC_C3X_Thr[j_C3];
    }

 
    for(Int_t j_C4=0; j_C4<72; j_C4++){
      ADC_C4X_R[j_C4] = adc_c4x_r[j_C4]-ADC_C4X_Thr[j_C4];
      ADC_C4X_L[j_C4] = adc_c4x_l[j_C4]-ADC_C4X_Thr[j_C4];
    }


    for(Int_t j_CY=0; j_CY<16; j_CY++){
      ADC_C2Y_R[j_CY] = adc_c2y_r[j_CY]-ADC_C2Y_Thr[j_CY];
      ADC_C2Y_L[j_CY] = adc_c2y_l[j_CY]-ADC_C2Y_Thr[j_CY];
      ADC_C3Y_R[j_CY] = adc_c3y_r[j_CY]-ADC_C3Y_Thr[j_CY];
      ADC_C3Y_L[j_CY] = adc_c3y_l[j_CY]-ADC_C3Y_Thr[j_CY];
      ADC_C4Y_R[j_CY] = adc_c4y_r[j_CY]-ADC_C4Y_Thr[j_CY];
      ADC_C4Y_L[j_CY] = adc_c4y_l[j_CY]-ADC_C4Y_Thr[j_CY];
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
       if((ADC_High_TARGET[i]>=0 && tdc_le_target[i][0]>=TARGET_TDC_min[i] && tdc_le_target[i][0]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]>=0 && tdc_le_target[i][1]>=TARGET_TDC_min[i] && tdc_le_target[i][1]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]>=0 && tdc_le_target[i][2]>=TARGET_TDC_min[i] && tdc_le_target[i][2]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]>=0 && tdc_le_target[i][3]>=TARGET_TDC_min[i] && tdc_le_target[i][3]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]<0  && ADC_Low_TARGET[i]>0))
      { 
       count_TARGET_evts++;
      }
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

    for(int i=0; i<12; i++)
    {

      if(ADC_TOF1U[i]>=0 || ADC_TOF1D[i]>=0)  Good_TOF1_ADC[i] = true;
      
      if((TDC_TOF1U[i]>=TOF1U_TDC_min[i] && TDC_TOF1U[i]<=TOF1U_TDC_max[i]) ||
         (TDC_TOF1D[i]>=TOF1D_TDC_min[i] && TDC_TOF1D[i]<=TOF1D_TDC_max[i]))  Good_TOF1_TDC[i] = true;

      if(Good_TOF1_TDC[i] || Good_TOF1_ADC[i]) Good_TOF1[i] = true;

      if((ADC_TOF2AO[i]>=0 && ADC_TOF2AI[i]>=0) || (ADC_TOF2BO[i]>=0 && ADC_TOF2BI[i]>=0))  Good_TOF2_ADC[i] = true;

      if(((TDC_TOF2AO[i]>=TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<=TOF2AO_TDC_max[i])  &&
          (TDC_TOF2AI[i]>=TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<=TOF2AI_TDC_max[i])) ||
         ((TDC_TOF2BO[i]>=TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<=TOF2BO_TDC_max[i])  &&
          (TDC_TOF2BI[i]>=TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<=TOF2BI_TDC_max[i])))  Good_TOF2_TDC[i] = true;

      if(Good_TOF2_ADC[i] || Good_TOF2_TDC[i])
       Good_TOF2[i] = true;
    }


    for(int k=0; k<12; k++)
    {
      if(k!=0 && k!=11)
      {
        if((Good_TOF2[k] && Good_TOF1[k-1]) || (Good_TOF2[k] && Good_TOF1[k]) || (Good_TOF2[k] && Good_TOF1[k+1]))
        {
          Good_TOFs[k] = true;
        }
      }
    }

    if((Good_TOF2[0] && Good_TOF1[11]) || (Good_TOF2[0] && Good_TOF1[0]) || (Good_TOF2[0] && Good_TOF1[1]))  Good_TOFs[0] = true;

    if((Good_TOF2[11] && Good_TOF1[10]) || (Good_TOF2[11] && Good_TOF1[11]) || (Good_TOF2[11] && Good_TOF1[0]))  Good_TOFs[11] = true;

    for(int kk=0; kk<12; kk++){   
      if(Good_TOFs[kk]) Good_TOF_Event = true;
    }     

    for(int i=0; i<12; i++){
      Good_tof1[i] = Good_TOF1[i];
      Good_tof2[i] = Good_TOF2[i];
    }

    if(enable_cout == 9) Good_TOF_Event = true;


    //********* GOOD MWPC EVENTS  
    int count_C2X = 0;    int count_C2Y = 0;
    int count_C3X = 0;    int count_C3Y = 0;
    int count_C4X = 0;    int count_C4Y = 0;
    bool Good_MWPC_Event = false;

    for(int i=128; i<254; i++){
      if(i!=142 && i!=143 && i!=158 && i!=159 && i!=174 && i!=175 && i!=190 && i!=191 && 
         i!=206 && i!=207 && i!=223 && i!=224 && i!=238 && i!=239){
      
        if(MWPCADC[i]>=0)  count_C2X++; 
      } 
    } 
    
    for(int ii=96; ii<128; ii++){
      
      if(MWPCADC[ii]>=0)  count_C2Y++;
    }
    
    for(int j=0; j<96; j++){

      if(MWPCADC[j]>=0) count_C3X++;
    }
    
    for(int j=480; j<511; j++){

      if(MWPCADC[j]>=0) count_C3X++;
    }

    for(int j=448; j<479; j++){

      if(MWPCADC[j]>=0) count_C3Y++;
    }
    
    for(int k=288; k<448; k++){
      if((k>=288 && k<=295) || (k>=304 && k<=311) || (k>=320 && k<=447)){
        if(MWPCADC[k]>=0) count_C4X++;
      }
    }

    for(int kk=256; kk<288; kk++){
      if(MWPCADC[kk]>=0)  count_C4Y++;
    }


    if(count_C2X>0 && count_C2Y>0 && count_C3X>0 && count_C3Y>0 && count_C4X>0 && count_C4Y>0)  Good_MWPC_Event = true;
  
    if(enable_cout == 9) Good_MWPC_Event = true;


    if(Good_TARGET_Event && Good_TOF_Event && Good_MWPC_Event && !Event_On_Blacklist)  Good_Event = true;
    //if(Good_TARGET_Event && Good_TOF_Event && Good_MWPC_Event)  Good_Event = true;

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

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
  TLine *TOF_line9 = new TLine(0.052539,0.38,0.052539,0.62);           TOF_line9->SetLineWidth(10);
  TLine *TOF_line8 = new TLine(0.172539,0.172154,0.052539,0.38);       TOF_line8->SetLineWidth(10);
  TLine *TOF_line7 = new TLine(0.172539,0.172154,0.38,0.0521539);      TOF_line7->SetLineWidth(10);
  TLine *TOF_line6 = new TLine(0.38,0.0521539,0.62,0.0521539);         TOF_line6->SetLineWidth(10);
  TLine *TOF_line5 = new TLine(0.62,0.0521539,0.8278461,0.172154);     TOF_line5->SetLineWidth(10);
  TLine *TOF_line4 = new TLine(0.8278461,0.172154,0.9478461,0.38);     TOF_line4->SetLineWidth(10);
  TLine *TOF_line3 = new TLine(0.9478461,0.38,0.9478461,0.62);         TOF_line3->SetLineWidth(10);
  TLine *TOF_line2 = new TLine(0.9478461,0.62,0.8278461,0.8278461);    TOF_line2->SetLineWidth(10);
  TLine *TOF_line1 = new TLine(0.8278461,0.8278461,0.62,0.9478461);    TOF_line1->SetLineWidth(10);
  TLine *TOF_line12 = new TLine(0.62,0.9478461,0.38,0.9478461);        TOF_line12->SetLineWidth(10);
  TLine *TOF_line11 = new TLine(0.38,0.9478461,0.172539,0.8278461);    TOF_line11->SetLineWidth(10);
  TLine *TOF_line10 = new TLine(0.052539,0.62,0.172539,0.8278461);     TOF_line10->SetLineWidth(10);

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


  ///////////////////////////////////////////////////////////////////////////////////////////////////

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

  TLatex *tex_palette_MWPC[10]; 
 
  tex_run_Title_TARGET = new TLatex(0.034,0.15,"Run :");
  tex_run_Title_TARGET->SetTextSize(0.038);
  tex_run_Title_TARGET->SetLineWidth(2);

  tex_run_TARGET = new TLatex(0.034,0.125,run_string);
  tex_run_TARGET->SetTextSize(0.038);
  tex_run_TARGET->SetLineWidth(2);

  tex_event_Title_TARGET = new TLatex(0.034,0.07,"Event :");
  tex_event_Title_TARGET->SetTextSize(0.038);
  tex_event_Title_TARGET->SetLineWidth(2);

  tex_event_TARGET = new TLatex(0.036,0.06,event_string);
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

  ////////////////////////////////////////////////////////////////////////////////////////



  //////////////////////////////////   MARKERS   //////////////////////////////////////////

  // ######### TARGET !

  TMarker *marker_ADC_TARGET[256];
  TMarker *marker_ADCL_TARGET[256];
  TMarker *marker_TDC_TARGET[256];
  TMarker *palette_TARGET[10];
  TMarker *palette_MWPC[10];

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

  //////////////////////////////////////////////////////////////////

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

 /////////////////////////////////////////////////////////////////


  TLatex *tex_DS_SFT;
  TLatex *tex_US_SFT;
  TLatex *tex_footer_SFT;

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


  TLatex *tex_palette_MWPC_C2[10]; 
  TMarker *palette_MWPC_C2[10];

  tex_DS_SFT = new TLatex(0.02,0.88,"SFT  --  Downstream & Upstream");
  tex_DS_SFT->SetTextSize(0.05);
  tex_DS_SFT->SetLineWidth(2);

  tex_US_SFT = new TLatex(0.02,0.41,"SFT  --  Upstream");
  tex_US_SFT->SetTextSize(0.05);
  tex_US_SFT->SetLineWidth(2);

  tex_footer_SFT = new TLatex(0.02,0.02,footer);
  tex_footer_SFT->SetTextSize(0.025);
  tex_footer_SFT->SetLineWidth(2);

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



  for(Int_t ii=0; ii<128; ii++){
    if(ADC_High_SFT_corr[ii]>ADC_High_corr_max){
      ADC_High_corr_max=ADC_High_SFT_corr[ii];
    }
  }

  for(Int_t ii=0; ii<128; ii++){
    for (Int_t qq=0; qq<6; qq++) {
      if (tdc_le_sft[ii][qq] > SFT_TDC_min[ii] && tdc_le_sft[ii][qq] < SFT_TDC_max[ii]) has_TDC_SFT_hit[ii]++;
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


  bool has_TDC_hit[256] = {false};
  bool has_TDC_hit_Kstop[256] = {false};

  for(Int_t i=0; i<256; i++){
    for (Int_t k=0; k<4; k++) {
      if ((tdc_le_target[i][k]>=TARGET_TDC_min[i]) && (tdc_le_target[i][k]<=TARGET_TDC_max[i])) has_TDC_hit[i] = true;
      if ((tdc_le_target[i][k]>=TDC_min_Kstop) && (tdc_le_target[i][k]<=TDC_max_Kstop)) has_TDC_hit_Kstop[i] = true;
   }
  }


  char ch_ADC_cut_TARGET[100];    sprintf(ch_ADC_cut_TARGET,"(ADC offset = %d)",TARGET_ADC_Thr_HG_Offset);
  char ch_ADC_and_TDC_cut[100];   sprintf(ch_ADC_and_TDC_cut,"(ADC offset = %d | %d #leq TDC #leq %d)",TARGET_ADC_Thr_HG_Offset,TDC_min_TARGET,TDC_max_TARGET);
  char ch_ADC_and_TDC_cut_Kstop[100];   sprintf(ch_ADC_and_TDC_cut_Kstop,"(ADC: HG #geq %d, LG #geq %d | %d #leq TDC K Stop #leq %d)",HG_KAON,LG_KAON,TDC_min_Kstop,TDC_max_Kstop);

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
    if(ADC_High_TARGET[q]>max_ADC) {
      max_index = q;
      max_ADC = ADC_High_TARGET[q];
    }
  }

  for(Int_t q=0; q<256; q++){
    if (q == max_index) continue;   
    else {
      if(ADC_High_TARGET[q]>max_ADC2) {
        max_index2 = q;
        max_ADC2 = ADC_High_TARGET[q];
      }
    }
  }

  for(Int_t q=0; q<256; q++){
    if ((q == max_index) || (q == max_index2)) continue;    
    else {
      if(ADC_High_TARGET[q]>max_ADC3) {
        max_index3 = q;
        max_ADC3 = ADC_High_TARGET[q];
      }
    }
  }

  for(Int_t q=0; q<256; q++){
    if ((q == max_index) || (q == max_index2) || (q == max_index3)) continue;   
    else {
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

  
  for(int i=0; i<12; i++){
    if(ADC_TOF1U[i]>0 && ADC_TOF1D[i]>0)
    has_both_ADC_TOF1_hit[i] = true;  
  }
  
  
    ///Set TOF1 Lines

  for(int i = 0; i < 12; i++){
    if (ADC_TOF1U[i]>0 || ADC_TOF1D[i]>0) {has_ADC_TOF1_hit[i] = true;}
    if ((TDC_TOF1U[i]>TOF1U_TDC_min[i] && TDC_TOF1U[i]<TOF1U_TDC_max[i]) || (TDC_TOF1D[i]>TOF1D_TDC_min[i] && TDC_TOF1D[i]<TOF1D_TDC_max[i])) {has_TDC_TOF1_hit[i] = true;}
  }

  if (has_TDC_TOF1_hit[0]) TOF_line1->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[0] && has_ADC_TOF1_hit[0]) TOF_line1->SetLineColor(3);

  if (has_TDC_TOF1_hit[1]) TOF_line2->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[1] && has_ADC_TOF1_hit[1]) TOF_line2->SetLineColor(3);

  if (has_TDC_TOF1_hit[2]) TOF_line3->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[2] && has_ADC_TOF1_hit[2]) TOF_line3->SetLineColor(3);

  if (has_TDC_TOF1_hit[3]) TOF_line4->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[3] && has_ADC_TOF1_hit[3]) TOF_line4->SetLineColor(3);

  if (has_TDC_TOF1_hit[4]) TOF_line5->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[4] && has_ADC_TOF1_hit[4]) TOF_line5->SetLineColor(3);

  if (has_TDC_TOF1_hit[5]) TOF_line6->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[5] && has_ADC_TOF1_hit[5]) TOF_line6->SetLineColor(3);

  if (has_TDC_TOF1_hit[6]) TOF_line7->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[6] && has_ADC_TOF1_hit[6]) TOF_line7->SetLineColor(3);

  if (has_TDC_TOF1_hit[7]) TOF_line8->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[7] && has_ADC_TOF1_hit[7]) TOF_line8->SetLineColor(3);

  if (has_TDC_TOF1_hit[8]) TOF_line9->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[8] && has_ADC_TOF1_hit[8]) TOF_line9->SetLineColor(3);

  if (has_TDC_TOF1_hit[9]) TOF_line10->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[9] && has_ADC_TOF1_hit[9]) TOF_line10->SetLineColor(3);

  if (has_TDC_TOF1_hit[10]) TOF_line11->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[10] && has_ADC_TOF1_hit[10]) TOF_line11->SetLineColor(3);

  if (has_TDC_TOF1_hit[11]) TOF_line12->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[11] && has_ADC_TOF1_hit[11]) TOF_line12->SetLineColor(3);

  int gap_to_fit_left = 0;
  int gap_to_fit_right = 0;

  int high_adc_u = 0;
  int high_adc_d = 0;

  int tie_count = 0;

  // for(int i = 0; i<6; i++){
  //   if(has_TDC_TOF1_hit[i] && has_ADC_TOF1_hit[i]){
  //     if(ADC_TOF1U[i] > high_adc_u && ADC_TOF1D[i] > high_adc_d){
  //       gap_to_fit_right = i+1;
  //       high_adc_u = ADC_TOF1U[i];
  //       high_adc_d = ADC_TOF1D[i];
  //     }
  //   }
  // }    

  for(int i = 0; i<6; i++){
    if(has_TDC_TOF1_hit[i] && has_ADC_TOF1_hit[i]){
      gap_to_fit_right = i+1;
      tie_count++;
    }
  }

  if(tie_count > 1){
    for(int i = 0; i<6; i++){
      if(has_TDC_TOF1_hit[i] && has_ADC_TOF1_hit[i]){
        if(ADC_TOF1U[i] > high_adc_u && ADC_TOF1D[i] > high_adc_d){
          gap_to_fit_right = i+1;
          high_adc_u = ADC_TOF1U[i];
          high_adc_d = ADC_TOF1D[i];
        }
      }
    }       
  }
  
  high_adc_u = 0;
  high_adc_d = 0;
  tie_count = 0;

  for(int i = 6; i<12; i++){
    if(has_TDC_TOF1_hit[i] && has_ADC_TOF1_hit[i]){
      gap_to_fit_left = i+1;
      tie_count++;
    }
  }  

  if(tie_count>1){
    for(int i = 6; i<12; i++){
      if(has_TDC_TOF1_hit[i] && has_ADC_TOF1_hit[i]){
        if(ADC_TOF1U[i] > high_adc_u && ADC_TOF1D[i] > high_adc_d){
          gap_to_fit_left = i+1;
          high_adc_u = ADC_TOF1U[i];
          high_adc_d = ADC_TOF1D[i];
        }
      }
    } 
  }



  int gap_hit[12] = {0};
  int ADC_TOF1_hit[12] = {0};
  int ADCTDC_TOF1_hit[12] = {0};
  int ADC_TOF2_hit[12] = {0};
  int ADCTDC_TOF2_hit[12] = {0};

 
  for(int j=0; j<12; j++){
    for(int k=0; k<8; k++){
      if (ADC_High_TARGET[channel[j][k]] > 0 && has_TDC_hit[channel[j][k]]) gap_hit[j] = 1;
    }
  }


  for (int k=0; k<12; k++) {
      if (has_ADC_TOF1_hit[k]) {
        if (has_TDC_TOF1_hit[k]) {ADCTDC_TOF1_hit[k]++;}
        else {ADC_TOF1_hit[k]++;}
      }
      if (has_ADC_TOF2_hit[k]) {
        if (has_TDC_TOF2_hit[k]) {ADCTDC_TOF2_hit[k]++;}
        else {ADC_TOF2_hit[k]++;}
      }
  }
  


  int gap_counter[12] = {0};


  //// GAP SCORING !
  
  int scoring_type = 2;         // scoring_type = 1  --->  Oscar's Method
                                // scoring_type = 2  --->  Sebastien's Method  TOF1[i]  TOF2[i-1], TOF2[i], TOF2[i+1]
                                // scoring_type = 3  --->  Sebastien's Method  TOF2[i]  TOF1[i-1], TOF1[i], TOF1[i+1]
  //Sebastien's Method
  //cout << "BIG TEST  !" << endl;

  if(scoring_type==2){
    for(int i=0; i<12; i++){
      if(ADC_TOF1U[i]>=0){
        gap_counter[i]++;
      }
    }

    for(int i=0; i<12; i++){
      if(ADC_TOF1D[i]>=0){
        gap_counter[i]++;
     }
    }

    for(int i=0; i<12; i++){
      if(TDC_TOF1U[i]>=TOF1U_TDC_min[i] && TDC_TOF1U[i]<=TOF1U_TDC_max[i]){
        gap_counter[i]++; 
      }
    }

    for(int i=0; i<12; i++){
      if(TDC_TOF1D[i]>=TOF1D_TDC_min[i] && TDC_TOF1D[i]<=TOF1D_TDC_max[i]){
        gap_counter[i]++;
      }
    }
  }
  

  int high_gap_hit = 0;
  int gap_to_fit = 0;



  vector<int> tof1_ties;
  int adc_tof1u_min = -10000;
  int adc_tof1d_min = -10000;

  for (int k=0; k<12; k++) {
    if ((gap_counter[k]>=high_gap_hit) && (has_TDC_TOF1_hit[k] || has_ADC_TOF1_hit[k])){
      if(high_gap_hit == gap_counter[k])
        tof1_ties.push_back(k);
      else{
        tof1_ties.clear();
        tof1_ties.push_back(k);
      }

      high_gap_hit = gap_counter[k];
      gap_to_fit = k+1;
    }
  }
  
  if(!tof1_ties.empty()){
    for(vector<int>::iterator it = tof1_ties.begin(); it != tof1_ties.end(); it++){
      if(ADC_TOF1U[*it] >= adc_tof1u_min && ADC_TOF1D[*it] >= adc_tof1d_min){
        gap_to_fit = *it + 1;
        adc_tof1u_min = ADC_TOF1U[*it];
        adc_tof1d_min = ADC_TOF1D[*it];
      }
    }
  }


  ////////////////////// ADDITIONAL


  // Determine the k_stop bars  
  bool k_stop_bar[256] = {false};
  vector<int> good_k_stop_bars;
  
  double X_BAR = 0;
  double Y_BAR = 0;
  
    
  double total_energy = 0.0;
  double X_weights = 0.0;
  double Y_weights = 0.0;   
    
  for(int i = 0; i<256; i++){
    if(ADC_High_TARGET_ped[i] > HG_KAON && ADC_Low_TARGET_ped[i] > LG_KAON && has_TDC_hit_Kstop[i]){
    good_k_stop_bars.push_back(i);
    k_stop_bar[i] = true;
    }
  }    

  
  // Determine if a hit target has any hit neighbours
  // bool TARGET_High_has_neighbours[256] = {false};
  for(int i = 0; i<256; i++){
    for(int j=0; j<8; j++){
      if((TARGET_neighbours[i][j] != -1 && ADC_High_TARGET[TARGET_neighbours[i][j]] >= Angle_ADC_cut && has_TDC_hit[TARGET_neighbours[i][j]] && !k_stop_bar[TARGET_neighbours[i][j]]) ||
         (TARGET_neighbours[i][j] != -1 && ADC_High_TARGET[TARGET_neighbours[i][j]] <0 && ADC_Low_TARGET[TARGET_neighbours[i][j]]>=0 && Switch==1 && !k_stop_bar[TARGET_neighbours[i][j]])){
           TARGET_High_has_neighbours[i] = true;
           break;
      }
    }
  }
  
  

  for(Int_t i=0; i<256; i++){
  if(TARGET_High_has_neighbours[i]){
      if(ADC_High_TARGET[i]>=Angle_ADC_cut && has_TDC_hit[i]){
        h_target_TDC->Fill(Xloc[i],Yloc[i]);
        h_target_TDC_copy->Fill(Xloc[i],Yloc[i]);
        
        if(!k_stop_bar[i]){
          h_TDC_Gap_Fibers->Fill(Xloc[i],Yloc[i]);
        }
        count++;
        

        if(IsIn(i,channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],   // Additional weight for fibers close to the edge if hit
                  channel[gap_to_fit-1][2], channel[gap_to_fit-1][3], 
                  channel[gap_to_fit-1][4], channel[gap_to_fit-1][5], 
                  channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){
          if(!k_stop_bar[i]){         
            h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
            h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
          }
          
          
          h_TDC_selected2->Fill(Xloc[i], Yloc[i]);
          h_TDC_selected2->Fill(Xloc[i], Yloc[i]);
  
          h_TDC_selected->Fill(Xloc[i], Yloc[i]);
          h_TDC_selected->Fill(Xloc[i], Yloc[i]);
        }
      }

      if(ADC_High_TARGET[i]>=Angle_ADC_cut){
        h_target_ADC->Fill(Xloc[i],Yloc[i]);
        x_inc = x_inc + Xloc[i];
        y_inc = y_inc + Yloc[i];
        hit_count++;
      }
    
      if(ADC_High_TARGET[i]<0 && ADC_Low_TARGET[i]>=0 && Switch==1){
        if(!k_stop_bar[i])  
          h_TDC_Gap_Fibers->Fill(Xloc[i],Yloc[i]);
        h_TARGET_LG_Blue->Fill(Xloc[i],Yloc[i]);
        count++;
      }
    }
  }
  ///////// DETERMINE FIBER WITH HIGHEST AND SECOND HIGHEST LOW GAIN AMPLITUDE
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
        if(ADC_Low_TARGET[i]>max_ADC_all[j]) {
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
  printf("LG ADC 1st Max.  :  %4d    %4d    %4d\n",index_max1, ADC_Low_TARGET_ped[index_max1], TDC_LG_max);
  printf("LG ADC 2nd Max.  :  %4d    %4d    %4d\n",index_max2, ADC_Low_TARGET_ped[index_max2], TDC_LG_max2);
  printf("TDC Average  :  %4d\n",TDC_average);

  cout << " " << endl;


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
    cout << "//////  TARGET  //////" << endl;
    cout << "Event: " << ievt << endl << endl;
    cout << "ADC Threshold Offset (HG): " << TARGET_ADC_Thr_HG_Offset << endl;
    cout << "ADC Threshold Offset (LG): " << TARGET_ADC_Thr_LG_Offset << endl << endl;
    cout << "Fiber  HG-Ped  LG-Ped    TDC[0]   TDC[1]   TDC[2]   TDC[3]   TDC[4]   TDC[5]  (HG) Thr  (LG) Thr  HG/LG" << endl;  // SEBASTIEN
    for(Int_t jj=0; jj<256; jj++){
      if(TDC_LE_TARGET[jj]>-1 || tdc_le_target[jj][1]>-1 || tdc_le_target[jj][2]>-1 || tdc_le_target[jj][3]>-1){
        printf("%3d     %4d    %4d     %4d     %4d     %4d     %4d     %4d     %4d      %4d      %4d     %2.2f\n",
       jj, ADC_High_TARGET[jj] + TARGET_ADC_Thr_HG_Offset, ADC_Low_TARGET[jj] + TARGET_ADC_Thr_LG_Offset, 
       TDC_LE_TARGET[jj], tdc_le_target[jj][1], tdc_le_target[jj][2], tdc_le_target[jj][3], tdc_le_target[jj][4], tdc_le_target[jj][5],
       HG_TARGET_ADC_Thr[jj], LG_TARGET_ADC_Thr[jj], float(adc_high_target[jj]-round(TARGET_ADC_Thr_HG[jj]))/float(adc_low_target[jj]-round(TARGET_ADC_Thr_LG[jj])));
      }
    }
  }      


  if(Good_Event && (enable_cout==0 || enable_cout==9)){
    cout << "//////  TARGET  //////" << endl;
    cout << "Event: " << ievt << endl << endl;
    cout << "ADC Threshold Offset (HG): " << TARGET_ADC_Thr_HG_Offset << endl;
    cout << "ADC Threshold Offset (LG): " << TARGET_ADC_Thr_LG_Offset << endl << endl;
    cout << "Fiber  HG-Ped   LG-Ped   T[0]-Av   T[1]-Av   T[2]-Av   T[3]-Av   T[4]-Av   T[5]-Av  (HG) Thr  (LG) Thr  HG/LG" << endl;
      for(Int_t jj=0; jj<256; jj++){
      if(ADC_High_TARGET[jj]>=0 || (ADC_Low_TARGET[jj]>=0 && Switch==1)){
        printf("%3d     %4d    %4d      %s       %s      %s      %s      %s      %s     %4d      %4d     %2.2f\n", 
        jj, ADC_High_TARGET[jj]+TARGET_ADC_Thr_HG_Offset, ADC_Low_TARGET[jj] + TARGET_ADC_Thr_LG_Offset, 
        TDC_LE_TARGET_corr[jj][0], TDC_LE_TARGET_corr[jj][1], TDC_LE_TARGET_corr[jj][2], TDC_LE_TARGET_corr[jj][3], TDC_LE_TARGET_corr[jj][4], TDC_LE_TARGET_corr[jj][5],
        HG_TARGET_ADC_Thr[jj], LG_TARGET_ADC_Thr[jj], float(adc_high_target[jj]-round(TARGET_ADC_Thr_HG[jj]))/float(adc_low_target[jj]-round(TARGET_ADC_Thr_LG[jj])));
      }
    }
  }    

  if(Good_Event && enable_cout==2){
    cout << "//////  TARGET  //////" << endl;
    cout << "Event: " << ievt << endl << endl;
    cout << "ADC Threshold Offset (HG): " << TARGET_ADC_Thr_HG_Offset << endl;
    cout << "ADC Threshold Offset (LG): " << TARGET_ADC_Thr_LG_Offset << endl << endl;
    cout << "Fiber  HG-Ped   LG-Ped   T[0]-Av   T[1]-Av   T[2]-Av   T[3]-Av   T[4]-Av   T[5]-Av  (HG) Thr  (LG) Thr  HG/LG" << endl;   // SEBASTIEN  
    for(Int_t jj=0; jj<256; jj++){
        printf("%3d     %4d     %4d      %4d      %4d      %4d      %4d      %4d      %4d     %4d      %4d     %2.2f\n", 
        jj, ADC_High_TARGET[jj]+TARGET_ADC_Thr_HG_Offset, ADC_Low_TARGET[jj] + TARGET_ADC_Thr_LG_Offset, 
        TDC_LE_TARGET_corrected[jj][0], TDC_LE_TARGET_corrected[jj][1], TDC_LE_TARGET_corrected[jj][2], TDC_LE_TARGET_corrected[jj][3], TDC_LE_TARGET_corrected[jj][4], TDC_LE_TARGET_corrected[jj][5],
        HG_TARGET_ADC_Thr[jj], LG_TARGET_ADC_Thr[jj], float(adc_high_target[jj]-round(TARGET_ADC_Thr_HG[jj]))/float(adc_low_target[jj]-round(TARGET_ADC_Thr_LG[jj])));
   }
  }    



  if(Good_Event && (enable_cout==0 || enable_cout==1 || enable_cout==9)){
    cout << " " << endl;
    cout << " " << endl;
    cout << "//////  SFT  //////" << endl << endl;
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
    cout << "Channel   Layer     Fiber     HG-Thr     LG-Thr   TDC[0]    TDC[1]    TDC[2]    TDC[3]   (HG) Thr   (LG) Thr" << endl;
    for(Int_t jj=0; jj<128; jj++){
      printf("%3d       %3d      %3d-%1c      %4d       %4d      %4d     %4d      %4d      %4d       %4d       %4d\n", 
      jj, SFT_channel_to_layer[jj], SFT_channel_to_fiber[jj], SFT_channel_to_US_DS[jj], ADC_High_SFT[jj], ADC_Low_SFT[jj],
      TDC_LE_SFT[jj], tdc_le_sft[jj][1], tdc_le_sft[jj][2], tdc_le_sft[jj][3], HG_SFT_ADC_Thr[jj], HG_SFT_ADC_Thr[jj]);
      
    }
  }


  cout << "" << endl;
  cout << "" << endl;
  for (Int_t i = 0; i<12; i++) {
    if(Good_Event){
      if (ADC_TOF1[i] > 0 || ADC_TOF1[i+12] > 0) {
        if ((TDC_TOF1U[i]>TOF1U_TDC_min[i] && TDC_TOF1U[i]<TOF1U_TDC_max[i]) || (TDC_TOF1D[i]>TOF1D_TDC_min[i] && TDC_TOF1D[i]<TOF1D_TDC_max[i])){
           printf(" ADC TOF1 Up-%2d: %5d -- Down-%2d: %5d  |  TDC TOF1 Up-%2d: %5d -- Down-%2d: %5d\n", i+1, ADC_TOF1[i], i+1, ADC_TOF1[i+12], i+1, TDC_TOF1U[i], i+1, TDC_TOF1D[i]);
        }
        else {
          printf(" ADC TOF1 Up-%2d: %5d -- Down-%2d: %5d  |\n", i+1, ADC_TOF1[i], i+1, ADC_TOF1[i+12]);
        }
      }
      else {
       if ((TDC_TOF1U[i]>TOF1U_TDC_min[i] && TDC_TOF1U[i]<TOF1U_TDC_max[i]) || (TDC_TOF1D[i]>TOF1D_TDC_min[i] && TDC_TOF1D[i]<TOF1D_TDC_max[i])) {
          printf("                                          |  TDC TOF1 Up-%2d: %5d -- Down-%2d: %5d\n", i+1, TDC_TOF1U[i], i+1, TDC_TOF1D[i]);
        }
      }
    }  
  }

  cout << "" << endl;
  cout << "" << endl;


  int C2X_L[56] = {0};
  int C2X_R[56] = {0};
  int C2Y_L[16] = {0};
  int C2Y_R[16] = {0};
  int C3X_L[64] = {0};
  int C3X_R[64] = {0};
  int C3Y_L[16] = {0};
  int C3Y_R[16] = {0};
  int C4X_L[72] = {0};
  int C4X_R[72] = {0};
  int C4Y_L[16] = {0};
  int C4Y_R[16] = {0};


  //C2 Counters

  for(int i=0; i<56; i++){
    if (Good_Event && ADC_C2X_L[i]>0){
      C2X_L[i] = ADC_C2X_L[i];
    }
    if (Good_Event && ADC_C2X_R[i]>0){
      C2X_R[i] = ADC_C2X_R[i];
    }
  }

  for(int i=0; i<64; i++){
    if (Good_Event && ADC_C3X_L[i]>0){
      C3X_L[i] = ADC_C3X_L[i];
    }
    if (Good_Event && ADC_C3X_R[i]>0){
      C3X_R[i] = ADC_C3X_R[i];
    }
  }   

  for(int i=0; i<72; i++){
    if (Good_Event && ADC_C4X_L[i]>0){
      C4X_L[i] = ADC_C4X_L[i];
    }
    if (Good_Event && ADC_C4X_R[i]>0){
      C4X_R[i] = ADC_C4X_R[i];
    }
  } 

  for(int i=0; i<16; i++){
    if (Good_Event && ADC_C2Y_L[i]>0) {
      C2Y_L[i] = ADC_C2Y_L[i];
      }
    if (Good_Event && ADC_C2Y_R[i]>0) {
      C2Y_R[i] = ADC_C2Y_R[i];
      }
    if (Good_Event && ADC_C3Y_L[i]>0) {
      C3Y_L[i] = ADC_C3Y_L[i];
      }
    if (Good_Event && ADC_C3Y_R[i]>0) {
      C3Y_R[i] = ADC_C3Y_R[i];
      }
    if (Good_Event && ADC_C4Y_L[i]>0) {
      C4Y_L[i] = ADC_C4Y_L[i];
      }
    if (Good_Event && ADC_C4Y_R[i]>0) {
      C4Y_R[i] = ADC_C4Y_R[i];
      }                
  }        
     

  if(gap_to_fit >= 7){ //LEFT

    for(int q=0; q<56; q++){ // C2 Counters
      if(Good_Event && ADC_C2X_L[q]>0){
        cout << " ADC MWPC Channel " << C2XL_Channels[q] << " -- " << "C2_X_" << q+1 << "_Left: " << ADC_C2X_L[q] << endl;
      }
    }  
    for(int q=0; q<16; q++){  
      if(Good_Event && ADC_C2Y_L[q]>0){
        cout << " ADC MWPC Channel " << C2YL_Channels[q] << " -- " << "C2_Y_" << q+1 << "_Left: " << ADC_C2Y_L[q] << endl;
      }
    }
    for(int q=0; q<56; q++){ // C2 Counters
      if(Good_Event && ADC_C2X_R[q]>0){
        cout << " ADC MWPC Channel " << C2XR_Channels[q] << " -- " << "C2_X_" << q+1 << "_Right: " << ADC_C2X_R[q] << endl;
      }
    }
    for(int q=0; q<16; q++){  
      if(Good_Event && ADC_C2Y_R[q]>0){
        cout << " ADC MWPC Channel " << C2YR_Channels[q] << " -- " << "C2_Y_" << q+1 << "_Right: " << ADC_C2Y_R[q] << endl;
      }
    }

    cout << " " << endl;

    for(int q=0; q<64; q++){ // C3 Counters
      if(Good_Event && ADC_C3X_L[q]>0){
        cout << " ADC MWPC Channel " << C3XL_Channels[q] << " -- " << "C3_X_" << q+1 << "_Left: " << ADC_C3X_L[q] << endl;
      }
    }
    for(int q=0; q<16; q++){  
      if(Good_Event && ADC_C3Y_L[q]>0){
        cout << " ADC MWPC Channel " << C3YL_Channels[q] << " -- " << "C3_Y_" << q+1 << "_Left: " << ADC_C3Y_L[q] << endl;
      }
    }
    for(int q=0; q<64; q++){ // C3 Counters
      if(Good_Event && ADC_C3X_R[q]>0){
        cout << " ADC MWPC Channel " << C3XR_Channels[q] << " -- " << "C3_X_" << q+1 << "_Right: " << ADC_C3X_R[q] << endl;
      }
    }
    for(int q=0; q<16; q++){ // C3 Counters
      if(Good_Event && ADC_C3Y_R[q]>0){
        cout << " ADC MWPC Channel " << C3YR_Channels[q] << " -- " << "C3_Y_" << q+1 << "_Right: " << ADC_C3Y_R[q] << endl;
      }
    }

    cout << " " << endl;

    for(int q=0; q<72; q++){ // C4 Counters
      if(Good_Event && ADC_C4X_L[q]>0){
        cout << " ADC MWPC Channel " << C4XL_Channels[q] << " -- " << "C4_X_" << q+1 << "_Left: " << ADC_C4X_L[q] << endl;
      }
    }
    for(int q=0; q<16; q++){   
      if(Good_Event && ADC_C4Y_L[q]>0){
        cout << " ADC MWPC Channel " << C4YL_Channels[q] << " -- " << "C4_Y_" << q+1 << "_Left: " << ADC_C4Y_L[q] << endl;
      }
    }
    for(int q=0; q<72; q++){ // C4 Counters
      if(Good_Event && ADC_C4X_R[q]>0){
        cout << " ADC MWPC Channel " << C4XR_Channels[q] << " -- " << "C4_X_" << q+1 << "_Right: " << ADC_C4X_R[q] << endl;
      }
    }
    for(int q=0; q<16; q++){   
      if(Good_Event && ADC_C4Y_R[q]>0){
        cout << " ADC MWPC Channel " << C4YR_Channels[q] << " -- " << "C4_Y_" << q+1 << "_Right: " << ADC_C4Y_R[q] << endl;
      }
    }
  }
  else{


    for(int q=0; q<56; q++){ // C2 Counters
      if(Good_Event && ADC_C2X_R[q]>0){
        cout << " ADC MWPC Channel " << C2XR_Channels[q] << " -- " << "C2_X_" << q+1 << "_Right: " << ADC_C2X_R[q] << endl;
      }
    }
    for(int q=0; q<16; q++){   
      if(Good_Event && ADC_C2Y_R[q]>0){ 
        cout << " ADC MWPC Channel " << C2YR_Channels[q] << " -- " << "C2_Y_" << q+1 << "_Right: " << ADC_C2Y_R[q] << endl;
      }
    }
    for(int q=0; q<56; q++){ // C2 Counters
      if(Good_Event && ADC_C2X_L[q]>0){ 
        cout << " ADC MWPC Channel " << C2XL_Channels[q] << " -- " << "C2_X_" << q+1 << "_Left: " << ADC_C2X_L[q] << endl;
      }
    }
    for(int q=0; q<16; q++){   
      if(Good_Event && ADC_C2Y_L[q]>0){
        cout << " ADC MWPC Channel " << C2YL_Channels[q] << " -- " << "C2_Y_" << q+1 << "_Left: " << ADC_C2Y_L[q] << endl;
      }
    }

    cout << " " << endl;

    for(int q=0; q<64; q++){ // C3 Counters
      if(Good_Event && ADC_C3X_R[q]>0){
        cout << " ADC MWPC Channel " << C3XR_Channels[q] << " -- " << "C3_X_" << q+1 << "_Right: " << ADC_C3X_R[q] << endl;
      }
    }
    for(int q=0; q<16; q++){   
      if(Good_Event && ADC_C3Y_R[q]>0){
        cout << " ADC MWPC Channel " << C3YR_Channels[q] << " -- " << "C3_Y_" << q+1 << "_Right: " << ADC_C3Y_R[q] << endl;
      }
    }
    for(int q=0; q<64; q++){ // C3 Counters
      if(Good_Event && ADC_C3X_L[q]>0){
        cout << " ADC MWPC Channel " << C3XL_Channels[q] << " -- " << "C3_X_" << q+1 << "_Left: " << ADC_C3X_L[q] << endl;
      }
    }
    for(int q=0; q<16; q++){   
      if(Good_Event && ADC_C3Y_L[q]>0){
        cout << " ADC MWPC Channel " << C3YL_Channels[q] << " -- " << "C3_Y_" << q+1 << "_Left: " << ADC_C3Y_L[q] << endl;
      }
    }

    cout << " " << endl;

    for(int q=0; q<72; q++){ // C4 Counters
      if(Good_Event && ADC_C4X_R[q]>0){
        cout << " ADC MWPC Channel " << C4XR_Channels[q] << " -- " << "C4_X_" << q+1 << "_Right: " << ADC_C4X_R[q] << endl;
      }
    }
    for(int q=0; q<16; q++){   
      if(Good_Event && ADC_C4Y_R[q]>0){
        cout << " ADC MWPC Channel " << C4YR_Channels[q] << " -- " << "C4_Y_" << q+1 << "_Right: " << ADC_C4Y_R[q] << endl;
      }
    }
    for(int q=0; q<72; q++){ // C4 Counters
      if(Good_Event && ADC_C4X_L[q]>0){
        cout << " ADC MWPC Channel " << C4XL_Channels[q] << " -- " << "C4_X_" << q+1 << "_Left: " << ADC_C4X_L[q] << endl;
      }
    }
    for(int q=0; q<16; q++){   
      if(Good_Event && ADC_C4Y_L[q]>0){
        cout << " ADC MWPC Channel " << C4YL_Channels[q] << " -- " << "C4_Y_" << q+1 << "_Left: " << ADC_C4Y_L[q] << endl;
      }
    }
  }

  cout << endl;
  cout << endl;

  /*
  for (int q = 224; q < 238; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
    }
  }
  for (int q = 192; q < 206; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
    }
  }
  for (int q = 160; q < 174; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
    }
  }
  for (int q = 128; q < 142; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
    }
  }
  for (int q = 96; q < 112; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
    }
  } 

  for (int q = 240; q < 254; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
    }
  }
  for (int q = 208; q < 222; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
    }
  }
  for (int q = 176; q < 190; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
    }
  }
  for (int q = 144; q < 158; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
    }
  }
  // for (int q = 96; q < 128; q++) {
  //   if (Good_Event && MWPCADC[q]>0) {
  //     cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
  //   }
  // }

  for (int q = 112; q < 128; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
    }
  } 
  cout << "" << endl;

  //C3 Counters

  for (int q = 64; q < 80; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C3X_L[q-64] = MWPCADC[q];
    }
  }
  for (int q = 32; q < 48; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C3X_L[q-16] = MWPCADC[q];
    }
  }
  for (int q = 0; q < 16; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C3X_L[q+32] = MWPCADC[q];
    }
  }
  for (int q = 480; q < 496; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C3X_L[q-432] = MWPCADC[q];
    }
  }
  for (int q = 448; q < 464; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C3Y_L[q-448] = MWPCADC[q];
    }
  }

  for (int q = 80; q < 96; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C3X_R[q-80] = MWPCADC[q];
    }
  }
  for (int q = 48; q < 64; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C3X_R[q-32] = MWPCADC[q];
    }
  }
  for (int q = 16; q < 32; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C3X_R[q+16] = MWPCADC[q];
    }
  }
  for (int q = 496; q < 512; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C3X_R[q-448] = MWPCADC[q];
    }
  }
  for (int q = 464; q < 480; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C3Y_R[q-464] = MWPCADC[q];
    }
  }

  cout << "" << endl;

  //C4 Counters

  for (int q = 416; q < 432; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C4X_L[q-416] = MWPCADC[q];
    }
  }
  for (int q = 384; q < 400; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C4X_L[q-368] = MWPCADC[q];
    }
  }
  for (int q = 352; q < 368; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C4X_L[q-320] = MWPCADC[q];
    }
  }
  for (int q = 320; q < 336; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C4X_L[q-272] = MWPCADC[q];
    }
  }
  for (int q = 288; q < 296; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C4X_L[q-224] = MWPCADC[q];
    }
  }
  for (int q = 256; q < 272; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C4Y_L[q-256] = MWPCADC[q];
    }
  }
  
  for (int q = 432; q < 448; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C4X_R[q-432] = MWPCADC[q];
    }
  }
  for (int q = 400; q < 416; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C4X_R[q-384] = MWPCADC[q];
    }
  }
  for (int q = 368; q < 384; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C4X_R[q-336] = MWPCADC[q];
    }
  }
  for (int q = 336; q < 352; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C4X_R[q-288] = MWPCADC[q];
    }
  }
  for (int q = 304; q < 312; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C4X_R[q-240] = MWPCADC[q];
    }
  }

  for (int q = 272; q < 288; q++) {
    if (Good_Event && MWPCADC[q]>0) {
      cout << " ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
      //C4Y_R[q-272] = MWPCADC[q];
    }
  }
  */


  // Find clustering of mwpcs
  ////////////////////////////
  bool first_cluster = true;
  int cluster_spacing = 0;
  int cluster_length_count = 0;

  vector<int> C2X_L_cluster_index; // Hold the starting index of each cluster.
  vector<int> C2X_L_cluster_length; // Hold the length of each cluster.
  vector<int> C2Y_L_cluster_index;
  vector<int> C2Y_L_cluster_length;

  vector<int> C3X_L_cluster_index;
  vector<int> C3X_L_cluster_length;
  vector<int> C3Y_L_cluster_index;
  vector<int> C3Y_L_cluster_length;    

  vector<int> C4X_L_cluster_index;
  vector<int> C4X_L_cluster_length;  
  vector<int> C4Y_L_cluster_index;
  vector<int> C4Y_L_cluster_length;  



  vector<int> C2X_R_cluster_index; // Hold the starting index of each cluster.
  vector<int> C2X_R_cluster_length; // Hold the length of each cluster.
  vector<int> C2Y_R_cluster_index;
  vector<int> C2Y_R_cluster_length;

  vector<int> C3X_R_cluster_index;
  vector<int> C3X_R_cluster_length;
  vector<int> C3Y_R_cluster_index;
  vector<int> C3Y_R_cluster_length;    

  vector<int> C4X_R_cluster_index;
  vector<int> C4X_R_cluster_length;  
  vector<int> C4Y_R_cluster_index;
  vector<int> C4Y_R_cluster_length;   

  int C2X_L_clusters = 0;
  int C2Y_L_clusters = 0;
  int C3X_L_clusters = 0;
  int C3Y_L_clusters = 0;
  int C4X_L_clusters = 0;
  int C4Y_L_clusters = 0;

  int C2X_R_clusters = 0;
  int C2Y_R_clusters = 0;
  int C3X_R_clusters = 0;
  int C3Y_R_clusters = 0;
  int C4X_R_clusters = 0;
  int C4Y_R_clusters = 0;  



  // count C2X clusters
  first_cluster = true;
  cluster_spacing = 0;
  cluster_length_count = 0;


  for(int i = 0; i<56; i++){
    if(C2X_L[i] > 0 && first_cluster){
      cluster_spacing = MWPC_cluster_separation + 1;
      first_cluster = false;
    }

  //cout << " -> " << i << "  " << cluster_spacing << "  " << MWPC_cluster_separation << "  " << first_cluster << endl;

    if(C2X_L[i] > 0){
      if(cluster_spacing > MWPC_cluster_separation){
        C2X_L_clusters++;
        C2X_L_cluster_index.push_back(i);
      }
      cluster_length_count++;
      cluster_spacing = 0;

      if(i+2 < 56 && C2X_L[i+1] <= 0 && C2X_L[i+2] <= 0){
        C2X_L_cluster_length.push_back(cluster_length_count);
      }
      else if(i + 2 == 56 && C2X_L[i+1] <= 0){
        C2X_L_cluster_length.push_back(cluster_length_count);
      }
      else if(i + 1 == 56){
        C2X_L_cluster_length.push_back(cluster_length_count);
      }

    }
    else{
      cluster_spacing++;
      if(i != 0 && C2X_L[i-1] <= 0){
        cluster_length_count = 0;
      }
      else if ( i != 0){
        cluster_length_count++;
      }
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
        C2Y_L_clusters++;
        C2Y_L_cluster_index.push_back(i);
      }
      cluster_length_count++;
      cluster_spacing = 0;

      if(i+2 < 16 && C2Y_L[i+1] <= 0 && C2Y_L[i+2] <= 0){
        C2Y_L_cluster_length.push_back(cluster_length_count);
      }
      else if(i + 2 == 16 && C2Y_L[i+1] <= 0){
        C2Y_L_cluster_length.push_back(cluster_length_count);
      }
      else if(i + 1 == 16){
        C2Y_L_cluster_length.push_back(cluster_length_count);
      }

    }
    else{
      cluster_spacing++;
      if(i != 0 && C2Y_L[i-1] <= 0){
        cluster_length_count = 0;
      }
      else if ( i != 0){
        cluster_length_count++;
      }
    }
  }      


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
        C2X_R_clusters++;
        C2X_R_cluster_index.push_back(i);
      }
      cluster_length_count++;
      cluster_spacing = 0;

      if(i+2 < 56 && C2X_R[i+1] <= 0 && C2X_R[i+2] <= 0){
        C2X_R_cluster_length.push_back(cluster_length_count);
      }
      else if(i + 2 == 56 && C2X_R[i+1] <= 0){
        C2X_R_cluster_length.push_back(cluster_length_count);
      }
      else if(i + 1 == 56){
        C2X_R_cluster_length.push_back(cluster_length_count);
      }

    }
    else{
      cluster_spacing++;
      if(i != 0 && C2X_R[i-1] <= 0){
        cluster_length_count = 0;
      }
      else if ( i != 0){
        cluster_length_count++;
      }
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
        C2Y_R_clusters++;
        C2Y_R_cluster_index.push_back(i);
      }
      cluster_length_count++;
      cluster_spacing = 0;

      if(i+2 < 16 && C2Y_R[i+1] <= 0 && C2Y_R[i+2] <= 0){
        C2Y_R_cluster_length.push_back(cluster_length_count);
      }
      else if(i + 2 == 16 && C2Y_R[i+1] <= 0){
        C2Y_R_cluster_length.push_back(cluster_length_count);
      }
      else if(i + 1 == 16){
        C2Y_R_cluster_length.push_back(cluster_length_count);
      }

    }
    else{
      cluster_spacing++;
      if(i != 0 && C2Y_R[i-1] <= 0){
        cluster_length_count = 0;
      }
      else if ( i != 0){
        cluster_length_count++;
      }
    }
  }      

 
     



  // Display MWPC counter information
  cout << endl;

  cout << " Left Clustering" << endl;
  cout << " C2X clusters = " << C2X_L_clusters << endl;
  cout << " C2Y clusters = " << C2Y_L_clusters << endl << endl; 

  cout << " Right Clustering" << endl;
  cout << " C2X clusters = " << C2X_R_clusters << endl;
  cout << " C2Y clusters = " << C2Y_R_clusters << endl << endl;



  // Calculate and display centroids for C2X-L, C2X-R, C2Y-L, C2Y-R

  const int max_length = 100;

  int C2X_L_cluster_index_arr[max_length];
  int C2X_R_cluster_index_arr[max_length];
  int C2Y_L_cluster_index_arr[max_length];
  int C2Y_R_cluster_index_arr[max_length];

  int C2X_L_cluster_length_arr[max_length];
  int C2X_R_cluster_length_arr[max_length];
  int C2Y_L_cluster_length_arr[max_length];
  int C2Y_R_cluster_length_arr[max_length];

  double C2X_L_centroids[max_length];
  double C2X_R_centroids[max_length];
  double C2Y_L_centroids[max_length];
  double C2Y_R_centroids[max_length];

  double C2X_L_centroids_rotate[max_length];
  double C2X_R_centroids_rotate[max_length];
  double C2Y_L_centroids_rotate[max_length];
  double C2Y_R_centroids_rotate[max_length]; 


  for(int i = 0; i< max_length; i++){
    C2X_L_cluster_index_arr[i] = -10000;
    C2X_R_cluster_index_arr[i] = -10000;
    C2Y_L_cluster_index_arr[i] = -10000;
    C2Y_R_cluster_index_arr[i] = -10000;

    C2X_L_cluster_length_arr[i] = -1;
    C2X_R_cluster_length_arr[i] = -1;
    C2Y_L_cluster_length_arr[i] = -1;
    C2Y_R_cluster_length_arr[i] = -1;

    C2X_L_centroids[i] = -10000.0;
    C2X_R_centroids[i] = -10000.0;
    C2Y_L_centroids[i] = -10000.0;
    C2Y_R_centroids[i] = -10000.0;

    C2X_L_centroids_rotate[i] = -10000.0;
    C2X_R_centroids_rotate[i] = -10000.0;
    C2Y_L_centroids_rotate[i] = -10000.0;
    C2Y_R_centroids_rotate[i] = -10000.0;      
  } 


  // Store all index/lengths in arrays.
  for(unsigned int i = 0; i < C2X_L_cluster_index.size(); i++){
    C2X_L_cluster_index_arr[i] = *(C2X_L_cluster_index.begin() + i);
  }
  for(unsigned int i = 0; i < C2X_R_cluster_index.size(); i++){
    C2X_R_cluster_index_arr[i] = *(C2X_R_cluster_index.begin() + i);
  }
  for(unsigned int i = 0; i < C2Y_L_cluster_index.size(); i++){
    C2Y_L_cluster_index_arr[i] = *(C2Y_L_cluster_index.begin() + i);
  }
  for(unsigned int i = 0; i < C2Y_R_cluster_index.size(); i++){
    C2Y_R_cluster_index_arr[i] = *(C2Y_R_cluster_index.begin() + i);
  }

  for(unsigned int i = 0; i < C2X_L_cluster_length.size(); i++){
    C2X_L_cluster_length_arr[i] = *(C2X_L_cluster_length.begin() + i);
  }
  for(unsigned int i = 0; i < C2X_R_cluster_length.size(); i++){
    C2X_R_cluster_length_arr[i] = *(C2X_R_cluster_length.begin() + i);
  }
  for(unsigned int i = 0; i < C2Y_L_cluster_length.size(); i++){
    C2Y_L_cluster_length_arr[i] = *(C2Y_L_cluster_length.begin() + i);
  }
  for(unsigned int i = 0; i < C2Y_R_cluster_index.size(); i++){
    C2Y_R_cluster_length_arr[i] = *(C2Y_R_cluster_length.begin() + i);
  }



  // Calculate all centroids for C2-all
  
  double centroid_sum = 0.0;
  double centroid_den = 0.0;


  for(unsigned int i = 0; i < 100; i++){
    centroid_sum = 0.0;
    centroid_den = 0.0;

    if(C2X_L_cluster_index_arr[i] < -1000)
      break;

    for(int j = C2X_L_cluster_index_arr[i]; j<C2X_L_cluster_index_arr[i] + C2X_L_cluster_length_arr[i]; j++){
      if(C2X_L[j] > 0){
        centroid_sum += C2X_L[j]*C2_ZLoc[j];
        centroid_den += C2X_L[j];
      }
    }
   C2X_L_centroids[i] = centroid_sum/centroid_den;
  }



  for(unsigned int i = 0; i < 100; i++){
    centroid_sum = 0.0;
    centroid_den = 0.0;

    if(C2X_R_cluster_index_arr[i] < -1000)
      break;

    for(int j = C2X_R_cluster_index_arr[i]; j<C2X_R_cluster_index_arr[i] + C2X_R_cluster_length_arr[i]; j++){
      if(C2X_R[j] > 0){
        centroid_sum += C2X_R[j]*C2_ZLoc[j];
        centroid_den += C2X_R[j];
      }
    }

    C2X_R_centroids[i] = centroid_sum/centroid_den;
  }



  for(unsigned int i = 0; i < 100; i++){
    centroid_sum = 0.0;
    centroid_den = 0.0;

    if(C2Y_L_cluster_index_arr[i] < -1000)
      break;

    for(int j = C2Y_L_cluster_index_arr[i]; j<C2Y_L_cluster_index_arr[i] + C2Y_L_cluster_length_arr[i]; j++){
      if(C2Y_L[j] > 0){
        centroid_sum += C2Y_L[j]*C2_YLoc[j];
        centroid_den += C2Y_L[j];
      }
    }

    C2Y_L_centroids[i] = centroid_sum/centroid_den;
  }



  for(unsigned int i = 0; i < 100; i++){
    centroid_sum = 0.0;
    centroid_den = 0.0;

    if(C2Y_R_cluster_index_arr[i] < -1000)
      break;

    for(int j = C2Y_R_cluster_index_arr[i]; j<C2Y_R_cluster_index_arr[i] + C2Y_R_cluster_length_arr[i]; j++){
      if(C2Y_R[j] > 0){
        centroid_sum += C2Y_R[j]*C2_YLoc[j];
        centroid_den += C2Y_R[j];
      }
    }

    C2Y_R_centroids[i] = centroid_sum/centroid_den;
  }        





  cout << endl;

  for(int i = 0; i < 100; i++){
    if(C2X_L_centroids[i] < -1000)
      break;

    cout << "C2X_L centroid = " << C2X_L_centroids[i] << endl;
  }

  cout << endl;

  for(int i = 0; i < 100; i++){
    if(C2X_R_centroids[i] < -1000)
      break;

    cout << "C2X_R centroid = " << C2X_R_centroids[i] << endl;
  }

  cout << endl;

  for(int i = 0; i < 100; i++){
    if(C2Y_L_centroids[i] < -1000)
      break;

    cout << "C2Y_L centroid = " << C2Y_L_centroids[i] << endl;
  }

  cout <<endl;

  for(int i = 0; i < 100; i++){
    if(C2Y_R_centroids[i] < -1000)
      break;

    cout << "C2Y_R centroid = " << C2Y_R_centroids[i] << endl;
  }       

  cout << endl << endl;     


    ////////////////////////////////////////////////////////////////
    /// Draw wire chamber hits.
    ////////////////////////////////////////////////////////////////  
    if(MWPC_switch==1){

      TCanvas *cMWPC;
      cMWPC = new TCanvas("Wire Chamber","Wire Chamber",50,50,1050,700);
      cMWPC->Divide(2,3);
      cMWPC->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

    
      int C2X_L_count = 0;
      int C2Y_L_count = 0;
      int C3X_L_count = 0;
      int C3Y_L_count = 0;
      int C4X_L_count = 0;
      int C4Y_L_count = 0;

      vector<int> C2X_L_wire_numbers;
      vector<int> C2Y_L_wire_numbers;
      vector<int> C3X_L_wire_numbers;
      vector<int> C3Y_L_wire_numbers;
      vector<int> C4X_L_wire_numbers;
      vector<int> C4Y_L_wire_numbers;

      int C2X_R_count = 0;
      int C2Y_R_count = 0;
      int C3X_R_count = 0;
      int C3Y_R_count = 0;
      int C4X_R_count = 0;
      int C4Y_R_count = 0;

      vector<int> C2X_R_wire_numbers;
      vector<int> C2Y_R_wire_numbers;
      vector<int> C3X_R_wire_numbers;
      vector<int> C3Y_R_wire_numbers;
      vector<int> C4X_R_wire_numbers;
      vector<int> C4Y_R_wire_numbers;      


      // counters for number of wc hits.
        for(int i=0; i<56; i++){
          if(C2X_L[i] > 0){
            C2X_L_count++;
            C2X_L_wire_numbers.push_back(i);
          }
        }

        for(int i=0; i<64; i++){
          if(C3X_L[i] > 0){
            C3X_L_count++;
            C3X_L_wire_numbers.push_back(i);
          }
        }

        for(int i=0; i<72; i++){
          if(C4X_L[i] > 0){
            C4X_L_count++;
            C4X_L_wire_numbers.push_back(i);
          }
        }

        for(int i=0; i<16; i++){
          if(C2Y_L[i] > 0){
            C2Y_L_count++;
            C2Y_L_wire_numbers.push_back(i);
          }
        }

        for(int i=0; i<16; i++){
          if(C3Y_L[i] > 0){
            C3Y_L_count++;
            C3Y_L_wire_numbers.push_back(i);
          }
        }

        for(int i=0; i<16; i++){
          if(C4Y_L[i] > 0){
            C4Y_L_count++;
            C4Y_L_wire_numbers.push_back(i);
          }
        }
      

        for(int i=0; i<56; i++){
          if(C2X_R[i] > 0){
            C2X_R_count++;
            C2X_R_wire_numbers.push_back(i);
          }
        }

        for(int i=0; i<64; i++){
          if(C3X_R[i] > 0){
            C3X_R_count++;
            C3X_R_wire_numbers.push_back(i);
          }
        }

        for(int i=0; i<72; i++){
          if(C4X_R[i] > 0){
            C4X_R_count++;
            C4X_R_wire_numbers.push_back(i);
          }
        }

        for(int i=0; i<16; i++){
          if(C2Y_R[i] > 0){
            C2Y_R_count++;
            C2Y_R_wire_numbers.push_back(i);
          }
        }

        for(int i=0; i<16; i++){
          if(C3Y_R[i] > 0){
            C3Y_R_count++;
            C3Y_R_wire_numbers.push_back(i);
          }
        }

        for(int i=0; i<16; i++){
          if(C4Y_R[i] > 0){
            C4Y_R_count++;
            C4Y_R_wire_numbers.push_back(i);
          }
        }
      

      TH2F *C2_L_hist = new TH2F("", "", 500, -280,280,500, -1,18);
      TH2F *C3_L_hist = new TH2F("", "", 500, -1,64,500, -1,18); 
      TH2F *C4_L_hist = new TH2F("", "", 500, -1,72,500, -1,18);

      TH2F *C2_R_hist = new TH2F("", "", 500, -280,280,500, -1,18);
      TH2F *C3_R_hist = new TH2F("", "", 500, -1,64,500, -1,18); 
      TH2F *C4_R_hist = new TH2F("", "", 500, -1,72,500, -1,18);      

      C2_L_hist->SetLabelSize(0.06,"X");
      C2_L_hist->SetLabelSize(0.06,"Y");
      C3_L_hist->SetLabelSize(0.06,"X");
      C3_L_hist->SetLabelSize(0.06,"Y");
      C4_L_hist->SetLabelSize(0.06,"X");
      C4_L_hist->SetLabelSize(0.06,"Y");

      C2_L_hist->SetNdivisions(18,"Y");
      C3_L_hist->SetNdivisions(18,"Y");
      C4_L_hist->SetNdivisions(18,"Y");

      C2_R_hist->SetLabelSize(0.06,"X");
      C2_R_hist->SetLabelSize(0.06,"Y");
      C3_R_hist->SetLabelSize(0.06,"X");
      C3_R_hist->SetLabelSize(0.06,"Y");
      C4_R_hist->SetLabelSize(0.06,"X");
      C4_R_hist->SetLabelSize(0.06,"Y");

      C2_R_hist->SetNdivisions(18,"Y");
      C3_R_hist->SetNdivisions(18,"Y");
      C4_R_hist->SetNdivisions(18,"Y"); 


      TGaxis *C2_top = new TGaxis(-280,18,280,18,"pol1",510,"+U");
      TGaxis *C2_right = new TGaxis(280,-1,280,18,"pol1",510,"-U");
      TGaxis *C3_top = new TGaxis(-1,18,64,18,"pol1",510,"+U");
      TGaxis *C3_right = new TGaxis(64,-1,64,18,"pol1",510,"-U");
      TGaxis *C4_top = new TGaxis(-1,18,72,18,"pol1",510,"+U");
      TGaxis *C4_right = new TGaxis(72,-1,72,18,"pol1",510,"-U");   
      

      TLine *C2X_L_line[C2X_L_count];
      TLine *C2Y_L_line[C2Y_L_count];
      TLine *C3X_L_line[C3X_L_count];
      TLine *C3Y_L_line[C3Y_L_count];
      TLine *C4X_L_line[C4X_L_count];
      TLine *C4Y_L_line[C4Y_L_count];

      TLine *C2X_R_line[C2X_R_count];
      TLine *C2Y_R_line[C2Y_R_count];
      TLine *C3X_R_line[C3X_R_count];
      TLine *C3Y_R_line[C3Y_R_count];
      TLine *C4X_R_line[C4X_R_count];
      TLine *C4Y_R_line[C4Y_R_count];      


      const int MWPC_color_thr[10] = {0, 150, 300, 450, 600, 750, 900, 1050, 1200, 1350};

      double ADC_palette_MWPC[10];
      for(int i = 0; i < 10; i++){
        ADC_palette_MWPC[i] = MWPC_color_thr[i]/100.0;  
      }
      char ADC_palette_string_MWPC[10][100];

      for(int j=0; j<10; j++) sprintf(ADC_palette_string_MWPC[j],"%2.1f",ADC_palette_MWPC[j]);
      sprintf(ADC_palette_string_MWPC[9],"%2.1f+",ADC_palette_MWPC[9]);

      double flag_size_palette_MWPC = 1;
      double text_size_MWPC = 0.04;

      double C2_palette_loc = -335;
      double C2_tex_palette_loc = -325;
      // C2 legend
      palette_MWPC_C2[0] = new TMarker(C2_palette_loc,16,21);     palette_MWPC_C2[0]->SetMarkerColor(kOrange+10);    palette_MWPC_C2[0]->SetMarkerSize(flag_size_palette_MWPC);
      palette_MWPC_C2[1] = new TMarker(C2_palette_loc,15,21);     palette_MWPC_C2[1]->SetMarkerColor(kOrange+7);     palette_MWPC_C2[1]->SetMarkerSize(flag_size_palette_MWPC);
      palette_MWPC_C2[2] = new TMarker(C2_palette_loc,14,21);     palette_MWPC_C2[2]->SetMarkerColor(kOrange+1);     palette_MWPC_C2[2]->SetMarkerSize(flag_size_palette_MWPC);
      palette_MWPC_C2[3] = new TMarker(C2_palette_loc,13,21);     palette_MWPC_C2[3]->SetMarkerColor(kOrange-4);     palette_MWPC_C2[3]->SetMarkerSize(flag_size_palette_MWPC);
      palette_MWPC_C2[4] = new TMarker(C2_palette_loc,12,21);     palette_MWPC_C2[4]->SetMarkerColor(kYellow-9);     palette_MWPC_C2[4]->SetMarkerSize(flag_size_palette_MWPC);
      palette_MWPC_C2[5] = new TMarker(C2_palette_loc,11,21);     palette_MWPC_C2[5]->SetMarkerColor(kYellow-7);     palette_MWPC_C2[5]->SetMarkerSize(flag_size_palette_MWPC);
      palette_MWPC_C2[6] = new TMarker(C2_palette_loc,10,21);     palette_MWPC_C2[6]->SetMarkerColor(kYellow-0);     palette_MWPC_C2[6]->SetMarkerSize(flag_size_palette_MWPC);
      palette_MWPC_C2[7] = new TMarker(C2_palette_loc,9,21);      palette_MWPC_C2[7]->SetMarkerColor(kSpring-4);     palette_MWPC_C2[7]->SetMarkerSize(flag_size_palette_MWPC);
      palette_MWPC_C2[8] = new TMarker(C2_palette_loc,8,21);      palette_MWPC_C2[8]->SetMarkerColor(kSpring-2);     palette_MWPC_C2[8]->SetMarkerSize(flag_size_palette_MWPC);
      palette_MWPC_C2[9] = new TMarker(C2_palette_loc,7,21);      palette_MWPC_C2[9]->SetMarkerColor(kGreen-0);      palette_MWPC_C2[9]->SetMarkerSize(flag_size_palette_MWPC);

      TLatex *tex_palette_MWPC_scale_C2 = new TLatex(C2_tex_palette_loc,6,"x 100");               
      tex_palette_MWPC_scale_C2->SetTextSize(text_size_MWPC);

      tex_palette_MWPC_C2[0] = new TLatex(C2_tex_palette_loc,16,ADC_palette_string_MWPC[0]);      tex_palette_MWPC_C2[0]->SetTextSize(text_size_MWPC);
      tex_palette_MWPC_C2[1] = new TLatex(C2_tex_palette_loc,15,ADC_palette_string_MWPC[1]);      tex_palette_MWPC_C2[1]->SetTextSize(text_size_MWPC);
      tex_palette_MWPC_C2[2] = new TLatex(C2_tex_palette_loc,14,ADC_palette_string_MWPC[2]);      tex_palette_MWPC_C2[2]->SetTextSize(text_size_MWPC);
      tex_palette_MWPC_C2[3] = new TLatex(C2_tex_palette_loc,13,ADC_palette_string_MWPC[3]);      tex_palette_MWPC_C2[3]->SetTextSize(text_size_MWPC);
      tex_palette_MWPC_C2[4] = new TLatex(C2_tex_palette_loc,12,ADC_palette_string_MWPC[4]);      tex_palette_MWPC_C2[4]->SetTextSize(text_size_MWPC);
      tex_palette_MWPC_C2[5] = new TLatex(C2_tex_palette_loc,11,ADC_palette_string_MWPC[5]);      tex_palette_MWPC_C2[5]->SetTextSize(text_size_MWPC);
      tex_palette_MWPC_C2[6] = new TLatex(C2_tex_palette_loc,10,ADC_palette_string_MWPC[6]);      tex_palette_MWPC_C2[6]->SetTextSize(text_size_MWPC);
      tex_palette_MWPC_C2[7] = new TLatex(C2_tex_palette_loc,9,ADC_palette_string_MWPC[7]);       tex_palette_MWPC_C2[7]->SetTextSize(text_size_MWPC);
      tex_palette_MWPC_C2[8] = new TLatex(C2_tex_palette_loc,8,ADC_palette_string_MWPC[8]);       tex_palette_MWPC_C2[8]->SetTextSize(text_size_MWPC);
      tex_palette_MWPC_C2[9] = new TLatex(C2_tex_palette_loc,7,ADC_palette_string_MWPC[9]);  



      palette_MWPC[0] = new TMarker(-7,16,21);     palette_MWPC[0]->SetMarkerColor(kOrange+10);    palette_MWPC[0]->SetMarkerSize(flag_size_palette_MWPC);
      palette_MWPC[1] = new TMarker(-7,15,21);     palette_MWPC[1]->SetMarkerColor(kOrange+7);     palette_MWPC[1]->SetMarkerSize(flag_size_palette_MWPC);
      palette_MWPC[2] = new TMarker(-7,14,21);     palette_MWPC[2]->SetMarkerColor(kOrange+1);     palette_MWPC[2]->SetMarkerSize(flag_size_palette_MWPC);
      palette_MWPC[3] = new TMarker(-7,13,21);     palette_MWPC[3]->SetMarkerColor(kOrange-4);     palette_MWPC[3]->SetMarkerSize(flag_size_palette_MWPC);
      palette_MWPC[4] = new TMarker(-7,12,21);     palette_MWPC[4]->SetMarkerColor(kYellow-9);     palette_MWPC[4]->SetMarkerSize(flag_size_palette_MWPC);
      palette_MWPC[5] = new TMarker(-7,11,21);     palette_MWPC[5]->SetMarkerColor(kYellow-7);     palette_MWPC[5]->SetMarkerSize(flag_size_palette_MWPC);
      palette_MWPC[6] = new TMarker(-7,10,21);     palette_MWPC[6]->SetMarkerColor(kYellow-0);     palette_MWPC[6]->SetMarkerSize(flag_size_palette_MWPC);
      palette_MWPC[7] = new TMarker(-7,9,21);      palette_MWPC[7]->SetMarkerColor(kSpring-4);     palette_MWPC[7]->SetMarkerSize(flag_size_palette_MWPC);
      palette_MWPC[8] = new TMarker(-7,8,21);      palette_MWPC[8]->SetMarkerColor(kSpring-2);     palette_MWPC[8]->SetMarkerSize(flag_size_palette_MWPC);
      palette_MWPC[9] = new TMarker(-7,7,21);      palette_MWPC[9]->SetMarkerColor(kGreen-0);      palette_MWPC[9]->SetMarkerSize(flag_size_palette_MWPC);

      TLatex *tex_palette_MWPC_scale = new TLatex(-6,6,"x 100");               tex_palette_MWPC_scale->SetTextSize(text_size_MWPC);
      tex_palette_MWPC[0] = new TLatex(-6,16,ADC_palette_string_MWPC[0]);      tex_palette_MWPC[0]->SetTextSize(text_size_MWPC);
      tex_palette_MWPC[1] = new TLatex(-6,15,ADC_palette_string_MWPC[1]);      tex_palette_MWPC[1]->SetTextSize(text_size_MWPC);
      tex_palette_MWPC[2] = new TLatex(-6,14,ADC_palette_string_MWPC[2]);      tex_palette_MWPC[2]->SetTextSize(text_size_MWPC);
      tex_palette_MWPC[3] = new TLatex(-6,13,ADC_palette_string_MWPC[3]);      tex_palette_MWPC[3]->SetTextSize(text_size_MWPC);
      tex_palette_MWPC[4] = new TLatex(-6,12,ADC_palette_string_MWPC[4]);      tex_palette_MWPC[4]->SetTextSize(text_size_MWPC);
      tex_palette_MWPC[5] = new TLatex(-6,11,ADC_palette_string_MWPC[5]);      tex_palette_MWPC[5]->SetTextSize(text_size_MWPC);
      tex_palette_MWPC[6] = new TLatex(-6,10,ADC_palette_string_MWPC[6]);      tex_palette_MWPC[6]->SetTextSize(text_size_MWPC);
      tex_palette_MWPC[7] = new TLatex(-6,9,ADC_palette_string_MWPC[7]);       tex_palette_MWPC[7]->SetTextSize(text_size_MWPC);
      tex_palette_MWPC[8] = new TLatex(-6,8,ADC_palette_string_MWPC[8]);       tex_palette_MWPC[8]->SetTextSize(text_size_MWPC);
      tex_palette_MWPC[9] = new TLatex(-6,7,ADC_palette_string_MWPC[9]);       tex_palette_MWPC[9]->SetTextSize(text_size_MWPC);



      TLatex *tex_Title_MWPC_C2_L;
      TLatex *tex_Title_MWPC_C3_L;
      TLatex *tex_Title_MWPC_C4_L;    

      TLatex *tex_Title_MWPC_C2_R;
      TLatex *tex_Title_MWPC_C3_R;
      TLatex *tex_Title_MWPC_C4_R;       

      char c2_L_title[200];
      char c3_L_title[200];
      char c4_L_title[200];

      char c2_R_title[200];
      char c3_R_title[200];
      char c4_R_title[200];      


 
      sprintf(c2_L_title,"C2 Left  |  Run Number %d  |   Event %d   |   C2 Thr = %d   ", Run_Number, ievt, MWPC_ADC_Thr[0]);
      sprintf(c3_L_title,"C3 Left  |  Run Number %d  |   Event %d   |   C3 Thr = %d   ", Run_Number, ievt, MWPC_ADC_Thr[0]);
      sprintf(c4_L_title,"C4 Left  |  Run Number %d  |   Event %d   |   C4 Thr = %d   ", Run_Number, ievt, MWPC_ADC_Thr[0]);

      sprintf(c2_R_title,"C2 Right  |  Run Number %d  |   Event %d   |   C2 Thr = %d   ", Run_Number, ievt, MWPC_ADC_Thr[0]);
      sprintf(c3_R_title,"C3 Right  |  Run Number %d  |   Event %d   |   C3 Thr = %d   ", Run_Number, ievt, MWPC_ADC_Thr[0]);
      sprintf(c4_R_title,"C4 Right  |  Run Number %d  |   Event %d   |   C4 Thr = %d   ", Run_Number, ievt, MWPC_ADC_Thr[0]);



      tex_Title_MWPC_C2_L = new TLatex(-280,18.5,c2_L_title);
      tex_Title_MWPC_C3_L = new TLatex(0,18.5,c3_L_title);
      tex_Title_MWPC_C4_L = new TLatex(0,18.5,c4_L_title);

      tex_Title_MWPC_C2_L->SetTextSize(0.07);
      tex_Title_MWPC_C2_L->SetLineWidth(2);
      tex_Title_MWPC_C2_L->Draw();

      tex_Title_MWPC_C3_L->SetTextSize(0.07);
      tex_Title_MWPC_C3_L->SetLineWidth(2);
      tex_Title_MWPC_C3_L->Draw();

      tex_Title_MWPC_C4_L->SetTextSize(0.07);
      tex_Title_MWPC_C4_L->SetLineWidth(2);
      tex_Title_MWPC_C4_L->Draw();



      tex_Title_MWPC_C2_R = new TLatex(-280,18.5,c2_R_title);
      tex_Title_MWPC_C3_R = new TLatex(0,18.5,c3_R_title);
      tex_Title_MWPC_C4_R = new TLatex(0,18.5,c4_R_title);

      tex_Title_MWPC_C2_R->SetTextSize(0.07);
      tex_Title_MWPC_C2_R->SetLineWidth(2);
      tex_Title_MWPC_C2_R->Draw();

      tex_Title_MWPC_C3_R->SetTextSize(0.07);
      tex_Title_MWPC_C3_R->SetLineWidth(2);
      tex_Title_MWPC_C3_R->Draw();

      tex_Title_MWPC_C4_R->SetTextSize(0.07);
      tex_Title_MWPC_C4_R->SetLineWidth(2);
      tex_Title_MWPC_C4_R->Draw();      



      for(int i = 0; i < C2X_L_count; i++){
        C2X_L_line[i] = new TLine(C2_ZLoc[*(C2X_L_wire_numbers.begin() + i)], -1, C2_ZLoc[*(C2X_L_wire_numbers.begin() + i)], 18);

          if(C2X_L[*(C2X_L_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C2X_L[*(C2X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C2X_L_line[i]->SetLineColor(kOrange+10);
          if(C2X_L[*(C2X_L_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C2X_L[*(C2X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C2X_L_line[i]->SetLineColor(kOrange+7);
          if(C2X_L[*(C2X_L_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C2X_L[*(C2X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C2X_L_line[i]->SetLineColor(kOrange+1);
          if(C2X_L[*(C2X_L_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C2X_L[*(C2X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C2X_L_line[i]->SetLineColor(kOrange-4);
          if(C2X_L[*(C2X_L_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C2X_L[*(C2X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C2X_L_line[i]->SetLineColor(kOrange-9);
          if(C2X_L[*(C2X_L_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C2X_L[*(C2X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C2X_L_line[i]->SetLineColor(kOrange-7);
          if(C2X_L[*(C2X_L_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C2X_L[*(C2X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C2X_L_line[i]->SetLineColor(kOrange-0);
          if(C2X_L[*(C2X_L_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C2X_L[*(C2X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C2X_L_line[i]->SetLineColor(kOrange-4);
          if(C2X_L[*(C2X_L_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C2X_L[*(C2X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C2X_L_line[i]->SetLineColor(kOrange-2);
          if(C2X_L[*(C2X_L_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C2X_L_line[i]->SetLineColor(kGreen-0);

        C2X_L_line[i]->SetLineWidth(4.5);
      }

      for(int i = 0; i < C2X_R_count; i++){
        C2X_R_line[i] = new TLine(C2_ZLoc[*(C2X_R_wire_numbers.begin() + i)], -1, C2_ZLoc[*(C2X_R_wire_numbers.begin() + i)], 18);

          if(C2X_R[*(C2X_R_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C2X_R[*(C2X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C2X_R_line[i]->SetLineColor(kOrange+10);
          if(C2X_R[*(C2X_R_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C2X_R[*(C2X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C2X_R_line[i]->SetLineColor(kOrange+7);
          if(C2X_R[*(C2X_R_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C2X_R[*(C2X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C2X_R_line[i]->SetLineColor(kOrange+1);
          if(C2X_R[*(C2X_R_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C2X_R[*(C2X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C2X_R_line[i]->SetLineColor(kOrange-4);
          if(C2X_R[*(C2X_R_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C2X_R[*(C2X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C2X_R_line[i]->SetLineColor(kOrange-9);
          if(C2X_R[*(C2X_R_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C2X_R[*(C2X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C2X_R_line[i]->SetLineColor(kOrange-7);
          if(C2X_R[*(C2X_R_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C2X_R[*(C2X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C2X_R_line[i]->SetLineColor(kOrange-0);
          if(C2X_R[*(C2X_R_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C2X_R[*(C2X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C2X_R_line[i]->SetLineColor(kOrange-4);
          if(C2X_R[*(C2X_R_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C2X_R[*(C2X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C2X_R_line[i]->SetLineColor(kOrange-2);
          if(C2X_R[*(C2X_R_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C2X_R_line[i]->SetLineColor(kGreen-0);
        
        C2X_R_line[i]->SetLineWidth(4.5);
      }



      for(int i = 0; i < C2Y_L_count; i++){
        C2Y_L_line[i] = new TLine(-280,*(C2Y_L_wire_numbers.begin() + i)+1, 280, *(C2Y_L_wire_numbers.begin() + i)+1);

          if(C2Y_L[*(C2Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C2Y_L[*(C2Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C2Y_L_line[i]->SetLineColor(kOrange+10);
          if(C2Y_L[*(C2Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C2Y_L[*(C2Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C2Y_L_line[i]->SetLineColor(kOrange+7);
          if(C2Y_L[*(C2Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C2Y_L[*(C2Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C2Y_L_line[i]->SetLineColor(kOrange+1);
          if(C2Y_L[*(C2Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C2Y_L[*(C2Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C2Y_L_line[i]->SetLineColor(kOrange-4);
          if(C2Y_L[*(C2Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C2Y_L[*(C2Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C2Y_L_line[i]->SetLineColor(kOrange-9);
          if(C2Y_L[*(C2Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C2Y_L[*(C2Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C2Y_L_line[i]->SetLineColor(kOrange-7);
          if(C2Y_L[*(C2Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C2Y_L[*(C2Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C2Y_L_line[i]->SetLineColor(kOrange-0);
          if(C2Y_L[*(C2Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C2Y_L[*(C2Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C2Y_L_line[i]->SetLineColor(kOrange-4);
          if(C2Y_L[*(C2Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C2Y_L[*(C2Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C2Y_L_line[i]->SetLineColor(kOrange-2);
          if(C2Y_L[*(C2Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C2Y_L_line[i]->SetLineColor(kGreen-0);

        C2Y_L_line[i]->SetLineWidth(4.5);
      }


      for(int i = 0; i < C2Y_R_count; i++){
        C2Y_R_line[i] = new TLine(-280,*(C2Y_R_wire_numbers.begin() + i)+1, 280, *(C2Y_R_wire_numbers.begin() + i)+1);

          if(C2Y_R[*(C2Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C2Y_R[*(C2Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C2Y_R_line[i]->SetLineColor(kOrange+10);
          if(C2Y_R[*(C2Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C2Y_R[*(C2Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C2Y_R_line[i]->SetLineColor(kOrange+7);
          if(C2Y_R[*(C2Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C2Y_R[*(C2Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C2Y_R_line[i]->SetLineColor(kOrange+1);
          if(C2Y_R[*(C2Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C2Y_R[*(C2Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C2Y_R_line[i]->SetLineColor(kOrange-4);
          if(C2Y_R[*(C2Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C2Y_R[*(C2Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C2Y_R_line[i]->SetLineColor(kOrange-9);
          if(C2Y_R[*(C2Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C2Y_R[*(C2Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C2Y_R_line[i]->SetLineColor(kOrange-7);
          if(C2Y_R[*(C2Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C2Y_R[*(C2Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C2Y_R_line[i]->SetLineColor(kOrange-0);
          if(C2Y_R[*(C2Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C2Y_R[*(C2Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C2Y_R_line[i]->SetLineColor(kOrange-4);
          if(C2Y_R[*(C2Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C2Y_R[*(C2Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C2Y_R_line[i]->SetLineColor(kOrange-2);
          if(C2Y_R[*(C2Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C2Y_R_line[i]->SetLineColor(kGreen-0);
        

        C2Y_R_line[i]->SetLineWidth(4.5);
      }


      for(int i = 0; i < C3X_L_count; i++){
        C3X_L_line[i] = new TLine(*(C3X_L_wire_numbers.begin() + i)+1, -1, *(C3X_L_wire_numbers.begin() + i)+1, 18);

          if(C3X_L[*(C3X_L_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C3X_L[*(C3X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C3X_L_line[i]->SetLineColor(kOrange+10);
          if(C3X_L[*(C3X_L_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C3X_L[*(C3X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C3X_L_line[i]->SetLineColor(kOrange+7);
          if(C3X_L[*(C3X_L_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C3X_L[*(C3X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C3X_L_line[i]->SetLineColor(kOrange+1);
          if(C3X_L[*(C3X_L_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C3X_L[*(C3X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C3X_L_line[i]->SetLineColor(kOrange-4);
          if(C3X_L[*(C3X_L_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C3X_L[*(C3X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C3X_L_line[i]->SetLineColor(kOrange-9);
          if(C3X_L[*(C3X_L_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C3X_L[*(C3X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C3X_L_line[i]->SetLineColor(kOrange-7);
          if(C3X_L[*(C3X_L_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C3X_L[*(C3X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C3X_L_line[i]->SetLineColor(kOrange-0);
          if(C3X_L[*(C3X_L_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C3X_L[*(C3X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C3X_L_line[i]->SetLineColor(kOrange-4);
          if(C3X_L[*(C3X_L_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C3X_L[*(C3X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C3X_L_line[i]->SetLineColor(kOrange-2);
          if(C3X_L[*(C3X_L_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C3X_L_line[i]->SetLineColor(kGreen-0);


        C3X_L_line[i]->SetLineWidth(4.5);
      }
      for(int i = 0; i < C3X_R_count; i++){
        C3X_R_line[i] = new TLine(*(C3X_R_wire_numbers.begin() + i)+1, -1, *(C3X_R_wire_numbers.begin() + i)+1, 18);

          if(C3X_R[*(C3X_R_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C3X_R[*(C3X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C3X_R_line[i]->SetLineColor(kOrange+10);
          if(C3X_R[*(C3X_R_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C3X_R[*(C3X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C3X_R_line[i]->SetLineColor(kOrange+7);
          if(C3X_R[*(C3X_R_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C3X_R[*(C3X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C3X_R_line[i]->SetLineColor(kOrange+1);
          if(C3X_R[*(C3X_R_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C3X_R[*(C3X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C3X_R_line[i]->SetLineColor(kOrange-4);
          if(C3X_R[*(C3X_R_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C3X_R[*(C3X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C3X_R_line[i]->SetLineColor(kOrange-9);
          if(C3X_R[*(C3X_R_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C3X_R[*(C3X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C3X_R_line[i]->SetLineColor(kOrange-7);
          if(C3X_R[*(C3X_R_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C3X_R[*(C3X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C3X_R_line[i]->SetLineColor(kOrange-0);
          if(C3X_R[*(C3X_R_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C3X_R[*(C3X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C3X_R_line[i]->SetLineColor(kOrange-4);
          if(C3X_R[*(C3X_R_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C3X_R[*(C3X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C3X_R_line[i]->SetLineColor(kOrange-2);
          if(C3X_R[*(C3X_R_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C3X_R_line[i]->SetLineColor(kGreen-0);
        
        C3X_R_line[i]->SetLineWidth(4.5);
      }



      for(int i = 0; i < C3Y_L_count; i++){
        C3Y_L_line[i] = new TLine(-1,*(C3Y_L_wire_numbers.begin() + i)+1, 64, *(C3Y_L_wire_numbers.begin() + i)+1);

          if(C3Y_L[*(C3Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C3Y_L[*(C3Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C3Y_L_line[i]->SetLineColor(kOrange+10);
          if(C3Y_L[*(C3Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C3Y_L[*(C3Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C3Y_L_line[i]->SetLineColor(kOrange+7);
          if(C3Y_L[*(C3Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C3Y_L[*(C3Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C3Y_L_line[i]->SetLineColor(kOrange+1);
          if(C3Y_L[*(C3Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C3Y_L[*(C3Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C3Y_L_line[i]->SetLineColor(kOrange-4);
          if(C3Y_L[*(C3Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C3Y_L[*(C3Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C3Y_L_line[i]->SetLineColor(kOrange-9);
          if(C3Y_L[*(C3Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C3Y_L[*(C3Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C3Y_L_line[i]->SetLineColor(kOrange-7);
          if(C3Y_L[*(C3Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C3Y_L[*(C3Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C3Y_L_line[i]->SetLineColor(kOrange-0);
          if(C3Y_L[*(C3Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C3Y_L[*(C3Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C3Y_L_line[i]->SetLineColor(kOrange-4);
          if(C3Y_L[*(C3Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C3Y_L[*(C3Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C3Y_L_line[i]->SetLineColor(kOrange-2);
          if(C3Y_L[*(C3Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C3Y_L_line[i]->SetLineColor(kGreen-0);

        C3Y_L_line[i]->SetLineWidth(4.5);
      }   

      for(int i = 0; i < C3Y_R_count; i++){
        C3Y_R_line[i] = new TLine(-1,*(C3Y_R_wire_numbers.begin() + i)+1, 64, *(C3Y_R_wire_numbers.begin() + i)+1);

          if(C3Y_R[*(C3Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C3Y_R[*(C3Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C3Y_R_line[i]->SetLineColor(kOrange+10);
          if(C3Y_R[*(C3Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C3Y_R[*(C3Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C3Y_R_line[i]->SetLineColor(kOrange+7);
          if(C3Y_R[*(C3Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C3Y_R[*(C3Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C3Y_R_line[i]->SetLineColor(kOrange+1);
          if(C3Y_R[*(C3Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C3Y_R[*(C3Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C3Y_R_line[i]->SetLineColor(kOrange-4);
          if(C3Y_R[*(C3Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C3Y_R[*(C3Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C3Y_R_line[i]->SetLineColor(kOrange-9);
          if(C3Y_R[*(C3Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C3Y_R[*(C3Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C3Y_R_line[i]->SetLineColor(kOrange-7);
          if(C3Y_R[*(C3Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C3Y_R[*(C3Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C3Y_R_line[i]->SetLineColor(kOrange-0);
          if(C3Y_R[*(C3Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C3Y_R[*(C3Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C3Y_R_line[i]->SetLineColor(kOrange-4);
          if(C3Y_R[*(C3Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C3Y_R[*(C3Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C3Y_R_line[i]->SetLineColor(kOrange-2);
          if(C3Y_R[*(C3Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C3Y_R_line[i]->SetLineColor(kGreen-0);

        C3Y_R_line[i]->SetLineWidth(4.5);
      }  


      for(int i = 0; i < C4X_L_count; i++){
        C4X_L_line[i] = new TLine(*(C4X_L_wire_numbers.begin() + i)+1, -1, *(C4X_L_wire_numbers.begin() + i)+1, 18);

          if(C4X_L[*(C4X_L_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C4X_L[*(C4X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C4X_L_line[i]->SetLineColor(kOrange+10);
          if(C4X_L[*(C4X_L_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C4X_L[*(C4X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C4X_L_line[i]->SetLineColor(kOrange+7);
          if(C4X_L[*(C4X_L_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C4X_L[*(C4X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C4X_L_line[i]->SetLineColor(kOrange+1);
          if(C4X_L[*(C4X_L_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C4X_L[*(C4X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C4X_L_line[i]->SetLineColor(kOrange-4);
          if(C4X_L[*(C4X_L_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C4X_L[*(C4X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C4X_L_line[i]->SetLineColor(kOrange-9);
          if(C4X_L[*(C4X_L_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C4X_L[*(C4X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C4X_L_line[i]->SetLineColor(kOrange-7);
          if(C4X_L[*(C4X_L_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C4X_L[*(C4X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C4X_L_line[i]->SetLineColor(kOrange-0);
          if(C4X_L[*(C4X_L_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C4X_L[*(C4X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C4X_L_line[i]->SetLineColor(kOrange-4);
          if(C4X_L[*(C4X_L_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C4X_L[*(C4X_L_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C4X_L_line[i]->SetLineColor(kOrange-2);
          if(C4X_L[*(C4X_L_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C4X_L_line[i]->SetLineColor(kGreen-0);

        C4X_L_line[i]->SetLineWidth(4.5);
      }

      for(int i = 0; i < C4X_R_count; i++){
        C4X_R_line[i] = new TLine(*(C4X_R_wire_numbers.begin() + i)+1, -1, *(C4X_R_wire_numbers.begin() + i)+1, 18);

          if(C4X_R[*(C4X_R_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C4X_R[*(C4X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C4X_R_line[i]->SetLineColor(kOrange+10);
          if(C4X_R[*(C4X_R_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C4X_R[*(C4X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C4X_R_line[i]->SetLineColor(kOrange+7);
          if(C4X_R[*(C4X_R_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C4X_R[*(C4X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C4X_R_line[i]->SetLineColor(kOrange+1);
          if(C4X_R[*(C4X_R_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C4X_R[*(C4X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C4X_R_line[i]->SetLineColor(kOrange-4);
          if(C4X_R[*(C4X_R_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C4X_R[*(C4X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C4X_R_line[i]->SetLineColor(kOrange-9);
          if(C4X_R[*(C4X_R_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C4X_R[*(C4X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C4X_R_line[i]->SetLineColor(kOrange-7);
          if(C4X_R[*(C4X_R_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C4X_R[*(C4X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C4X_R_line[i]->SetLineColor(kOrange-0);
          if(C4X_R[*(C4X_R_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C4X_R[*(C4X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C4X_R_line[i]->SetLineColor(kOrange-4);
          if(C4X_R[*(C4X_R_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C4X_R[*(C4X_R_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C4X_R_line[i]->SetLineColor(kOrange-2);
          if(C4X_R[*(C4X_R_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C4X_R_line[i]->SetLineColor(kGreen-0);
        
        C4X_R_line[i]->SetLineWidth(4.5);
      }

      for(int i = 0; i < C4Y_L_count; i++){
        C4Y_L_line[i] = new TLine(-1,*(C4Y_L_wire_numbers.begin() + i)+1, 72, *(C4Y_L_wire_numbers.begin() + i)+1);

          if(C4Y_L[*(C4Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C4Y_L[*(C4Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C4Y_L_line[i]->SetLineColor(kOrange+10);
          if(C4Y_L[*(C4Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C4Y_L[*(C4Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C4Y_L_line[i]->SetLineColor(kOrange+7);
          if(C4Y_L[*(C4Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C4Y_L[*(C4Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C4Y_L_line[i]->SetLineColor(kOrange+1);
          if(C4Y_L[*(C4Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C4Y_L[*(C4Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C4Y_L_line[i]->SetLineColor(kOrange-4);
          if(C4Y_L[*(C4Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C4Y_L[*(C4Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C4Y_L_line[i]->SetLineColor(kOrange-9);
          if(C4Y_L[*(C4Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C4Y_L[*(C4Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C4Y_L_line[i]->SetLineColor(kOrange-7);
          if(C4Y_L[*(C4Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C4Y_L[*(C4Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C4Y_L_line[i]->SetLineColor(kOrange-0);
          if(C4Y_L[*(C4Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C4Y_L[*(C4Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C4Y_L_line[i]->SetLineColor(kOrange-4);
          if(C4Y_L[*(C4Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C4Y_L[*(C4Y_L_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C4Y_L_line[i]->SetLineColor(kOrange-2);
          if(C4Y_L[*(C4Y_L_wire_numbers.begin() + i)] > MWPC_color_thr[9]                                                               ) C4Y_L_line[i]->SetLineColor(kGreen-0);

        C4Y_L_line[i]->SetLineWidth(4.5);
      }    


      for(int i = 0; i < C4Y_R_count; i++){
        C4Y_R_line[i] = new TLine(-1,*(C4Y_R_wire_numbers.begin() + i)+1, 72, *(C4Y_R_wire_numbers.begin() + i)+1);

          if(C4Y_R[*(C4Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[0] && C4Y_R[*(C4Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[1]) C4Y_R_line[i]->SetLineColor(kOrange+10);
          if(C4Y_R[*(C4Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[1] && C4Y_R[*(C4Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[2]) C4Y_R_line[i]->SetLineColor(kOrange+7);
          if(C4Y_R[*(C4Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[2] && C4Y_R[*(C4Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[3]) C4Y_R_line[i]->SetLineColor(kOrange+1);
          if(C4Y_R[*(C4Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[3] && C4Y_R[*(C4Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[4]) C4Y_R_line[i]->SetLineColor(kOrange-4);
          if(C4Y_R[*(C4Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[4] && C4Y_R[*(C4Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[5]) C4Y_R_line[i]->SetLineColor(kOrange-9);
          if(C4Y_R[*(C4Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[5] && C4Y_R[*(C4Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[6]) C4Y_R_line[i]->SetLineColor(kOrange-7);
          if(C4Y_R[*(C4Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[6] && C4Y_R[*(C4Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[7]) C4Y_R_line[i]->SetLineColor(kOrange-0);
          if(C4Y_R[*(C4Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[7] && C4Y_R[*(C4Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[8]) C4Y_R_line[i]->SetLineColor(kOrange-4);
          if(C4Y_R[*(C4Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[8] && C4Y_R[*(C4Y_R_wire_numbers.begin() + i)] <= MWPC_color_thr[9]) C4Y_R_line[i]->SetLineColor(kOrange-2);
          if(C4Y_R[*(C4Y_R_wire_numbers.begin() + i)] > MWPC_color_thr[6]                                                               ) C4Y_R_line[i]->SetLineColor(kGreen-0);

        C4Y_R_line[i]->SetLineWidth(4.5);
      }           

      gStyle->SetOptStat(0); // Remove legneds          


      cMWPC->cd(1);

      cMWPC->cd(1)->Update();
      cMWPC->cd(1)->Range(-1,-1,56,16);

      C2_L_hist->Draw();

      tex_Title_MWPC_C2_L->Draw();

      for(Int_t ipal=0; ipal<10; ipal++)  palette_MWPC_C2[ipal]->Draw();
      for(Int_t ileg=0; ileg<10; ileg++)  tex_palette_MWPC_C2[ileg]->Draw();
      tex_palette_MWPC_scale_C2->Draw(); 


      C2_top->Draw("same");
      C2_right->Draw("same");

      for(int i = 0; i < C2X_L_count; i++){
        C2X_L_line[i]->Draw();
      }

      for(int i = 0; i < C2Y_L_count; i++){
        C2Y_L_line[i]->Draw();
      }   



      cMWPC->cd(2);

      cMWPC->cd(2)->Update();
      cMWPC->cd(2)->Range(-1,-1,56,16);

      C2_R_hist->Draw();

      tex_Title_MWPC_C2_R->Draw();

      for(Int_t ipal=0; ipal<10; ipal++)  palette_MWPC_C2[ipal]->Draw();
      for(Int_t ileg=0; ileg<10; ileg++)  tex_palette_MWPC_C2[ileg]->Draw();
      tex_palette_MWPC_scale_C2->Draw(); 


      C2_top->Draw("same");
      C2_right->Draw("same");

      for(int i = 0; i < C2X_R_count; i++){
        C2X_R_line[i]->Draw();
      }

      for(int i = 0; i < C2Y_R_count; i++){
        C2Y_R_line[i]->Draw();
      }   



      cMWPC->cd(3);

      cMWPC->cd(3)->Update();
      cMWPC->cd(3)->Range(-1,-1,64,16);

      C3_L_hist->Draw();
      tex_Title_MWPC_C3_L->Draw();

      for(Int_t ipal=0; ipal<10; ipal++)  palette_MWPC[ipal]->Draw();
      for(Int_t ileg=0; ileg<10; ileg++)  tex_palette_MWPC[ileg]->Draw();
      tex_palette_MWPC_scale->Draw(); 


      C3_top->Draw("same");
      C3_right->Draw("same");

      for(int i = 0; i < C3X_L_count; i++){
        C3X_L_line[i]->Draw();
      }

      for(int i = 0; i < C3Y_L_count; i++){
        C3Y_L_line[i]->Draw();
      }     


      cMWPC->cd(4);

      cMWPC->cd(4)->Update();
      cMWPC->cd(4)->Range(-1,-1,64,16);

      C3_R_hist->Draw();
      tex_Title_MWPC_C3_R->Draw();

      for(Int_t ipal=0; ipal<10; ipal++)  palette_MWPC[ipal]->Draw();
      for(Int_t ileg=0; ileg<10; ileg++)  tex_palette_MWPC[ileg]->Draw();
      tex_palette_MWPC_scale->Draw(); 


      C3_top->Draw("same");
      C3_right->Draw("same");

      for(int i = 0; i < C3X_R_count; i++){
        C3X_R_line[i]->Draw();
      }

      for(int i = 0; i < C3Y_R_count; i++){
        C3Y_R_line[i]->Draw();
      }   


      cMWPC->cd(5);

      cMWPC->cd(5)->Update();
      cMWPC->cd(5)->Range(-1,-1,72,16);

      C4_L_hist->Draw();
      tex_Title_MWPC_C4_L->Draw();

      for(Int_t ipal=0; ipal<10; ipal++)  palette_MWPC[ipal]->Draw();
      for(Int_t ileg=0; ileg<10; ileg++)  tex_palette_MWPC[ileg]->Draw();
      tex_palette_MWPC_scale->Draw();        

      C4_top->Draw("same");
      C4_right->Draw("same");

      for(int i = 0; i < C4X_L_count; i++){
        C4X_L_line[i]->Draw();
      }

      for(int i = 0; i < C4Y_L_count; i++){
        C4Y_L_line[i]->Draw();
      }    
    


      cMWPC->cd(6);

      cMWPC->cd(6)->Update();
      cMWPC->cd(6)->Range(-1,-1,72,16);

      C4_R_hist->Draw();
      tex_Title_MWPC_C4_R->Draw();

      for(Int_t ipal=0; ipal<10; ipal++)  palette_MWPC[ipal]->Draw();
      for(Int_t ileg=0; ileg<10; ileg++)  tex_palette_MWPC[ileg]->Draw();
      tex_palette_MWPC_scale->Draw();        

      C4_top->Draw("same");
      C4_right->Draw("same");

      for(int i = 0; i < C4X_R_count; i++){
        C4X_R_line[i]->Draw();
      }

      for(int i = 0; i < C4Y_R_count; i++){
        C4Y_R_line[i]->Draw();
      }    
    }

    cout << "Event: " << ievt << endl;


  TCanvas *c2;
  c2 = new TCanvas("Event_Display.C  --  TARGET & SFT","Event_Display.C  --  TARGET & SFT",0,200,1050,700);
  c2->Divide(3,2);
  c2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");


  c2->cd(1);
  hline1->Draw();     vline1->Draw();   
  hline2->Draw();     vline2->Draw();
  hline3->Draw();     vline3->Draw();
  hline4->Draw();     vline4->Draw();
  hline5->Draw();     vline5->Draw();
  hline6->Draw();     vline6->Draw();
  hline7->Draw();     vline7->Draw();
  hline8->Draw();     vline8->Draw();
  hline9->Draw();     vline9->Draw(); 
  hline10->Draw();    vline10->Draw();
  hline11->Draw();    vline11->Draw();
  hline12->Draw();    vline12->Draw();
  hline13->Draw();    vline13->Draw();
  hline14->Draw();    vline14->Draw();
  hline15->Draw();    vline15->Draw();
  hline16->Draw();    vline16->Draw();
  hline17->Draw();    vline17->Draw();
  hline18->Draw();    vline18->Draw();
  hline19->Draw();    vline19->Draw();

 
  TOF_line13->Draw();
  TOF_line14->Draw();
  TOF_line15->Draw();
  TOF_line16->Draw();
  TOF_line17->Draw();
  TOF_line18->Draw();
  TOF_line19->Draw();
  TOF_line20->Draw();
  TOF_line21->Draw();
  TOF_line22->Draw();
  TOF_line23->Draw();
  TOF_line24->Draw();

  TOF_line1->Draw();
  TOF_line2->Draw();
  TOF_line3->Draw();
  TOF_line4->Draw();
  TOF_line5->Draw();
  TOF_line6->Draw();
  TOF_line7->Draw();
  TOF_line8->Draw();
  TOF_line9->Draw();
  TOF_line10->Draw();
  TOF_line11->Draw();
  TOF_line12->Draw();

  for(Int_t ileg=0; ileg<36; ileg++)    tex_Legend_TARGET[ileg]->Draw();

  tex_event_TARGET->Draw();

  tex_Title_ADC_High_TARGET = new TLatex(0.01759134,0.9295171,"ADC HG Cut");
  tex_Title_ADC_High_TARGET->SetTextSize(0.07645875);
  tex_Title_ADC_High_TARGET->SetLineWidth(2);
  tex_Title_ADC_High_TARGET->Draw();

  
  tex_Subtitle_ADC_High_TARGET = new TLatex(0.01759134,0.88,ch_ADC_cut_TARGET);
  tex_Subtitle_ADC_High_TARGET->SetTextSize(0.04);
  tex_Subtitle_ADC_High_TARGET->SetLineWidth(2);
  tex_Subtitle_ADC_High_TARGET->Draw();
  



  for(Int_t icol=0; icol<256; icol++){
    if(ADC_High_TARGET[icol]>=0    && ADC_High_TARGET[icol]<300 ){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange+10); marker_ADC_TARGET[icol]->Draw();} 
    if(ADC_High_TARGET[icol]>=300  && ADC_High_TARGET[icol]<600 ){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange+7); marker_ADC_TARGET[icol]->Draw();}    
    if(ADC_High_TARGET[icol]>=600  && ADC_High_TARGET[icol]<900 ){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange+1); marker_ADC_TARGET[icol]->Draw();}    
    if(ADC_High_TARGET[icol]>=900  && ADC_High_TARGET[icol]<1200){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange-4); marker_ADC_TARGET[icol]->Draw();}   
    if(ADC_High_TARGET[icol]>=1200 && ADC_High_TARGET[icol]<1500){ marker_ADC_TARGET[icol]->SetMarkerColor(kYellow-9); marker_ADC_TARGET[icol]->Draw();}    
    if(ADC_High_TARGET[icol]>=1500 && ADC_High_TARGET[icol]<1800){ marker_ADC_TARGET[icol]->SetMarkerColor(kYellow-7); marker_ADC_TARGET[icol]->Draw();}    
    if(ADC_High_TARGET[icol]>=1800 && ADC_High_TARGET[icol]<2100){ marker_ADC_TARGET[icol]->SetMarkerColor(kYellow-0); marker_ADC_TARGET[icol]->Draw();}    
    if(ADC_High_TARGET[icol]>=2100 && ADC_High_TARGET[icol]<2400){ marker_ADC_TARGET[icol]->SetMarkerColor(kSpring-4); marker_ADC_TARGET[icol]->Draw();}    
    if(ADC_High_TARGET[icol]>=2400 && ADC_High_TARGET[icol]<2700){ marker_ADC_TARGET[icol]->SetMarkerColor(kSpring-2); marker_ADC_TARGET[icol]->Draw();}    
    if(ADC_High_TARGET[icol]>=2700){ marker_ADC_TARGET[icol]->SetMarkerColor(kGreen-0); marker_ADC_TARGET[icol]->Draw();} 
  
  if(ADC_High_TARGET[icol]<0 && ADC_Low_TARGET[icol]>0 && TDC_min_TARGET && Switch==1){
    if(has_TDC_hit[icol])
      marker_ADC_TARGET[icol]->SetMarkerColor(11);
    else
      marker_ADC_TARGET[icol]->SetMarkerColor(kBlue);
    
    marker_ADC_TARGET[icol]->Draw();
    }     
  }

  for(Int_t ipal=0; ipal<10; ipal++)  palette_TARGET[ipal]->Draw();
  for(Int_t ileg=0; ileg<10; ileg++)  tex_palette_TARGET[ileg]->Draw();
  tex_palette_TARGET_scale->Draw();
  
  
  if (ADC_High_TARGET[max_index] > 0) {
    marker_TDC_TARGET[max_index]->SetMarkerColor(kViolet+1);
    marker_TDC_TARGET[max_index]->Draw();
  }

  if (ADC_High_TARGET[max_index2] > 0) {
    marker_TDC_TARGET[max_index2]->SetMarkerColor(1);
    marker_TDC_TARGET[max_index2]->Draw();
  }

  if (ADC_High_TARGET[max_index3] > 0) {
    marker_TDC_TARGET[max_index3]->SetMarkerColor(1);
    marker_TDC_TARGET[max_index3]->Draw();
  }

  if (ADC_High_TARGET[max_index4] > 0) {
    marker_TDC_TARGET[max_index4]->SetMarkerColor(1);
    marker_TDC_TARGET[max_index4]->Draw();
  }

  vblue2->Draw();
  vblue4->Draw();
  vblue6->Draw();
  vblue8->Draw();

  c2->cd(2);
  hline1->Draw();     vline1->Draw();   
  hline2->Draw();     vline2->Draw();
  hline3->Draw();     vline3->Draw();
  hline4->Draw();     vline4->Draw();
  hline5->Draw();     vline5->Draw();
  hline6->Draw();     vline6->Draw();
  hline7->Draw();     vline7->Draw();
  hline8->Draw();     vline8->Draw();
  hline9->Draw();     vline9->Draw();
  hline10->Draw();    vline10->Draw();
  hline11->Draw();    vline11->Draw();
  hline12->Draw();    vline12->Draw();
  hline13->Draw();    vline13->Draw();
  hline14->Draw();    vline14->Draw();
  hline15->Draw();    vline15->Draw();
  hline16->Draw();    vline16->Draw();
  hline17->Draw();    vline17->Draw();
  hline18->Draw();    vline18->Draw();
  hline19->Draw();    vline19->Draw();

  TOF_line13->Draw();
  TOF_line14->Draw();
  TOF_line15->Draw();
  TOF_line16->Draw();
  TOF_line17->Draw();
  TOF_line18->Draw();
  TOF_line19->Draw();
  TOF_line20->Draw();
  TOF_line21->Draw();
  TOF_line22->Draw();
  TOF_line23->Draw();
  TOF_line24->Draw();

  TOF_line1->Draw();
  TOF_line2->Draw();
  TOF_line3->Draw();
  TOF_line4->Draw();
  TOF_line5->Draw();
  TOF_line6->Draw();
  TOF_line7->Draw();
  TOF_line8->Draw();
  TOF_line9->Draw();
  TOF_line10->Draw();
  TOF_line11->Draw();
  TOF_line12->Draw();

  for(Int_t ileg=0; ileg<36; ileg++)    tex_Legend_TARGET[ileg]->Draw();

  tex_event_TARGET->Draw();

  tex_Title_ADC_High_TARGET = new TLatex(0.01759134,0.9295171,"ADC HG & TDC Cut");
  tex_Title_ADC_High_TARGET->SetTextSize(0.07645875);
  tex_Title_ADC_High_TARGET->SetLineWidth(2);
  tex_Title_ADC_High_TARGET->Draw();

  tex_Subtitle_ADC_High_TARGET = new TLatex(0.01759134,0.88,ch_ADC_and_TDC_cut);
  tex_Subtitle_ADC_High_TARGET->SetTextSize(0.04);
  tex_Subtitle_ADC_High_TARGET->SetLineWidth(2);
  tex_Subtitle_ADC_High_TARGET->Draw();

  for(Int_t icol=0; icol<256; icol++){
    if (has_TDC_hit[icol]) {
      if(ADC_High_TARGET[icol]>=0 && ADC_High_TARGET[icol]<300){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange+10); marker_ADC_TARGET[icol]->Draw();} 
      if(ADC_High_TARGET[icol]>=300 && ADC_High_TARGET[icol]<600){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange+7); marker_ADC_TARGET[icol]->Draw();}    
      if(ADC_High_TARGET[icol]>=600 && ADC_High_TARGET[icol]<900){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange+1); marker_ADC_TARGET[icol]->Draw();}    
      if(ADC_High_TARGET[icol]>=900 && ADC_High_TARGET[icol]<1200){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange-4); marker_ADC_TARGET[icol]->Draw();}   
      if(ADC_High_TARGET[icol]>=1200 && ADC_High_TARGET[icol]<1500){ marker_ADC_TARGET[icol]->SetMarkerColor(kYellow-9); marker_ADC_TARGET[icol]->Draw();}    
      if(ADC_High_TARGET[icol]>=1500 && ADC_High_TARGET[icol]<1800){ marker_ADC_TARGET[icol]->SetMarkerColor(kYellow-7); marker_ADC_TARGET[icol]->Draw();}    
      if(ADC_High_TARGET[icol]>=1800 && ADC_High_TARGET[icol]<2100){ marker_ADC_TARGET[icol]->SetMarkerColor(kYellow-0); marker_ADC_TARGET[icol]->Draw();}    
      if(ADC_High_TARGET[icol]>=2100 && ADC_High_TARGET[icol]<2400){ marker_ADC_TARGET[icol]->SetMarkerColor(kSpring-4); marker_ADC_TARGET[icol]->Draw();}    
      if(ADC_High_TARGET[icol]>=2400 && ADC_High_TARGET[icol]<2700){ marker_ADC_TARGET[icol]->SetMarkerColor(kSpring-2); marker_ADC_TARGET[icol]->Draw();}    
      if(ADC_High_TARGET[icol]>=2700){ marker_ADC_TARGET[icol]->SetMarkerColor(kGreen-0); marker_ADC_TARGET[icol]->Draw();}
    } 
    
  if(ADC_High_TARGET[icol]<0 && ADC_Low_TARGET[icol]>0 && TDC_min_TARGET && Switch==1){
    if(has_TDC_hit[icol])
    marker_ADC_TARGET[icol]->SetMarkerColor(11);
    else
      marker_ADC_TARGET[icol]->SetMarkerColor(kBlue);
    
    marker_ADC_TARGET[icol]->Draw();
    } 

  }

  for(Int_t ipal=0; ipal<10; ipal++)  palette_TARGET[ipal]->Draw();
  for(Int_t ileg=0; ileg<10; ileg++)  tex_palette_TARGET[ileg]->Draw();
  tex_palette_TARGET_scale->Draw();


  if ((ADC_Low_TARGET[max_index] > 0) && (ADC_High_TARGET[max_index] > 0) && (has_TDC_hit[max_index])) {
    marker_TDC_TARGET[max_index]->SetMarkerColor(kViolet+1);
    marker_TDC_TARGET[max_index]->Draw();
  }

  if ((ADC_Low_TARGET[max_index2] > 0) && (ADC_High_TARGET[max_index2] > 0) && (has_TDC_hit[max_index2])) {
    marker_TDC_TARGET[max_index2]->SetMarkerColor(1);
    marker_TDC_TARGET[max_index2]->Draw();
  }

  if ((ADC_Low_TARGET[max_index3] > 0) && (ADC_High_TARGET[max_index3] > 0) && (has_TDC_hit[max_index3])) {
    marker_TDC_TARGET[max_index3]->SetMarkerColor(1);
    marker_TDC_TARGET[max_index3]->Draw();
  }

  if ((ADC_Low_TARGET[max_index4] > 0) && (ADC_High_TARGET[max_index4] > 0) && (has_TDC_hit[max_index4])) {
    marker_TDC_TARGET[max_index4]->SetMarkerColor(1);
    marker_TDC_TARGET[max_index4]->Draw();
  }

  vblue2->Draw();
  vblue4->Draw(); 
  vblue6->Draw();
  vblue8->Draw();

  Double_t par0_ADC = 0;  float par0_TDC = 0.;
  Double_t par1_ADC = 0;  float par1_TDC = 0.;
  Double_t phi_TDC = 0;

  char angle_string_ADC[100]; char angle_string_TDC[100];
  
  TLatex *tex_angle_ADC;  TLatex *tex_angle_TDC;
  
  c2->cd(3);

  hline1->Draw();     vline1->Draw();   
  hline2->Draw();     vline2->Draw();
  hline3->Draw();     vline3->Draw();
  hline4->Draw();     vline4->Draw();
  hline5->Draw();     vline5->Draw(); 
  hline6->Draw();     vline6->Draw();
  hline7->Draw();     vline7->Draw();
  hline8->Draw();     vline8->Draw();
  hline9->Draw();     vline9->Draw();
  hline10->Draw();    vline10->Draw();
  hline11->Draw();    vline11->Draw();
  hline12->Draw();    vline12->Draw();
  hline13->Draw();    vline13->Draw();
  hline14->Draw();    vline14->Draw();
  hline15->Draw();    vline15->Draw();  
  hline16->Draw();    vline16->Draw();
  hline17->Draw();    vline17->Draw();
  hline18->Draw();    vline18->Draw();
  hline19->Draw();    vline19->Draw();

  TOF_line13->Draw();
  TOF_line14->Draw();
  TOF_line15->Draw();
  TOF_line16->Draw();
  TOF_line17->Draw();
  TOF_line18->Draw();
  TOF_line19->Draw();
  TOF_line20->Draw();
  TOF_line21->Draw();
  TOF_line22->Draw();
  TOF_line23->Draw();
  TOF_line24->Draw();

  TOF_line1->Draw();
  TOF_line2->Draw();
  TOF_line3->Draw();
  TOF_line4->Draw();
  TOF_line5->Draw();
  TOF_line6->Draw();
  TOF_line7->Draw();
  TOF_line8->Draw();
  TOF_line9->Draw();
  TOF_line10->Draw();
  TOF_line11->Draw();
  TOF_line12->Draw();

  for(Int_t ipal=0; ipal<10; ipal++)  palette_TARGET[ipal]->Draw();
  for(Int_t ileg=0; ileg<10; ileg++)  tex_palette_TARGET2[ileg]->Draw();
  tex_palette_TARGET_scale->Draw();
  
  for(Int_t ileg=0; ileg<36; ileg++)    tex_Legend_TARGET[ileg]->Draw();

  tex_event_TARGET->Draw();

  TLatex *tex_Title_ADC_Low_TARGET = new TLatex(0.01759134,0.9295171,"K Stop Cuts");
  tex_Title_ADC_Low_TARGET->SetTextSize(0.07);
  tex_Title_ADC_Low_TARGET->SetLineWidth(2);
  tex_Title_ADC_Low_TARGET->Draw();


  TLatex *tex_Subtitle_ADC_Low_TARGET = new TLatex(0.01759134,0.88,ch_ADC_and_TDC_cut_Kstop);
  tex_Subtitle_ADC_Low_TARGET->SetTextSize(0.04);
  tex_Subtitle_ADC_Low_TARGET->SetLineWidth(2);
  tex_Subtitle_ADC_Low_TARGET->Draw();
  

   
  /// Determine centroid for K Stop

        
  // Compute energy weighted centroids

  for(vector<int>::iterator it = good_k_stop_bars.begin(); it != good_k_stop_bars.end(); it++){
    X_weights += ADC_Low_TARGET_ped[*it]*Xloc[*it];
    Y_weights += ADC_Low_TARGET_ped[*it]*Yloc[*it];
    total_energy += ADC_Low_TARGET_ped[*it];
  }  
  
  X_BAR = X_weights/total_energy;
  Y_BAR = Y_weights/total_energy;    
  
  
  //Determine closest bar to the centroid, and draw in red
  double closest_to_centroid = 1000;
  int closest_to_centroid_index = -1;
  
  for(int i = 0; i<256; i++){
    if(distance(X_BAR,Y_BAR,Xloc[i],Yloc[i]) <= closest_to_centroid){
      closest_to_centroid_index = i;
      closest_to_centroid = distance(X_BAR,Y_BAR,Xloc[i],Yloc[i]);
    }
  }
  
  // 3rd window will show the condition for a kstop event
  for(Int_t icol=0; icol<256; icol++){
    if (has_TDC_hit_Kstop[icol]) {
      if(ADC_Low_TARGET_ped[icol]>=LG_KAON && ADC_High_TARGET_ped[icol] >=HG_KAON){ 
    marker_ADCL_TARGET[icol]->SetMarkerColor(12); marker_ADCL_TARGET[icol]->Draw();
    h_K_Stop_Bars->Fill(Xloc[icol],Yloc[icol]);
    }     
    } 
  }  
  
  if(!good_k_stop_bars.empty()){
    marker_ADCL_TARGET[closest_to_centroid_index]->SetMarkerColor(2); 
    marker_ADCL_TARGET[closest_to_centroid_index]->Draw();
  }

  char ch_centroid[200];    sprintf(ch_centroid,"#bar{x} = %3.2f, #bar{y} = %3.2f",X_BAR, Y_BAR);
  
  TLatex *tex_Label_Centroid = new TLatex(0.55,0.9295171,ch_centroid);
  tex_Label_Centroid->SetTextSize(0.05);
  tex_Label_Centroid->SetLineWidth(2);
  tex_Label_Centroid->Draw();  

  vblue2->Draw();
  vblue4->Draw();
  vblue6->Draw();
  vblue8->Draw();
  
  
  
  
  
  
  

  c2->cd(4);
  for(Int_t hdraw_DS=0; hdraw_DS<8; hdraw_DS++) hline_DS[hdraw_DS]->Draw();
  for(Int_t hdraw_US=0; hdraw_US<8; hdraw_US++) hline_US[hdraw_US]->Draw();

  for(Int_t vdraw_DS=0; vdraw_DS<68; vdraw_DS++)  vline_DS[vdraw_DS]->Draw();
  for(Int_t vdraw_US=0; vdraw_US<68; vdraw_US++)  vline_US[vdraw_US]->Draw();

  tex_DS_SFT->Draw();
  tex_footer_SFT->Draw();

  tex_number1_SFT->Draw();  tex_number10_SFT->Draw();
  tex_number2_SFT->Draw();  tex_number11_SFT->Draw();
  tex_number3_SFT->Draw();  tex_number12_SFT->Draw();
  tex_number4_SFT->Draw();  tex_number13_SFT->Draw();
  tex_number5_SFT->Draw();  tex_number14_SFT->Draw();
  tex_number6_SFT->Draw();  tex_number15_SFT->Draw();
  tex_number7_SFT->Draw();  tex_number16_SFT->Draw();
  tex_number8_SFT->Draw();  tex_number17_SFT->Draw();
  tex_number9_SFT->Draw();

  tex_number21_SFT->Draw();   tex_number21_bis_SFT->Draw();
  tex_number22_SFT->Draw();   tex_number22_bis_SFT->Draw();
  tex_number23_SFT->Draw();   tex_number23_bis_SFT->Draw();
  tex_number24_SFT->Draw();   tex_number24_bis_SFT->Draw();

  for(Int_t imark1=0; imark1<15; imark1++){
    if((ADC_High_SFT_corr[par_temp[1][imark1]]!=0) && (has_TDC_SFT_hit[par_temp[1][imark1]])){
      marker_DS[par_temp[0][imark1]]->SetMarkerColor(2);
      marker_DS[par_temp[0][imark1]]->Draw();
    }

    if((ADC_High_SFT_corr[par_temp[1][imark1+64]]!=0) && (has_TDC_SFT_hit[par_temp[1][imark1+64]])){
      marker_US[par_temp[0][imark1+64]]->SetMarkerColor(2);
      marker_US[par_temp[0][imark1+64]]->Draw();
    }
    
    if(ADC_High_SFT[imark1]<0 && ADC_Low_SFT[imark1]>0 && Switch==1){
      marker_DS[par_temp[0][imark1]]->SetMarkerColor(kGray);
      marker_DS[par_temp[0][imark1]]->Draw();
    } 
    
    if(ADC_High_SFT[imark1+64]<0 && ADC_Low_SFT[imark1+64]>0 && Switch==1){
      marker_US[par_temp[0][imark1+64]]->SetMarkerColor(kGray);
      marker_US[par_temp[0][imark1+64]]->Draw();
    }
  }


  for(Int_t imark2=0; imark2<15; imark2++){
    if((ADC_High_SFT_corr[par_temp[1][imark2+15]]!=0) && (has_TDC_SFT_hit[par_temp[1][imark2+15]])){
      marker_DS[par_temp[0][imark2+15]]->SetMarkerColor(4);
      marker_DS[par_temp[0][imark2+15]]->Draw();
    }

    if((ADC_High_SFT_corr[par_temp[1][imark2+79]]!=0) && (has_TDC_SFT_hit[par_temp[1][imark2+79]])){
      marker_US[par_temp[0][imark2+79]]->SetMarkerColor(4);
      marker_US[par_temp[0][imark2+79]]->Draw();
    }
    
    if(ADC_High_SFT[imark2+15]<0 && ADC_Low_SFT[imark2+15]>0 && Switch==1){
      marker_DS[par_temp[0][imark2+15]]->SetMarkerColor(kGray);
      marker_DS[par_temp[0][imark2+15]]->Draw();
    }
    
    if(ADC_High_SFT[imark2+79]<0 && ADC_Low_SFT[imark2+79]>0 && Switch==1){
      marker_US[par_temp[0][imark2+79]]->SetMarkerColor(kGray);
      marker_US[par_temp[0][imark2+79]]->Draw();
    }
  }


  for(Int_t imark3=0; imark3<17; imark3++){
    if((ADC_High_SFT_corr[par_temp[1][imark3+30]]!=0) && (has_TDC_SFT_hit[par_temp[1][imark3+30]])){
      marker_DS[par_temp[0][imark3+30]]->SetMarkerColor(3);
      marker_DS[par_temp[0][imark3+30]]->Draw();
    }

    if((ADC_High_SFT_corr[par_temp[1][imark3+94]]!=0) && (has_TDC_SFT_hit[par_temp[1][imark3+94]])){
      marker_US[par_temp[0][imark3+94]]->SetMarkerColor(3);
      marker_US[par_temp[0][imark3+94]]->Draw();
    }
    
    if(ADC_High_SFT[imark3+30]<0 && ADC_Low_SFT[imark3+30]>0 && Switch==1){
      marker_DS[par_temp[0][imark3+30]]->SetMarkerColor(kGray);
      marker_DS[par_temp[0][imark3+30]]->Draw();
    }
    
    if(ADC_High_SFT[imark3+94]<0 && ADC_Low_SFT[imark3+94]>0 && Switch==1){
      marker_US[par_temp[0][imark3+94]]->SetMarkerColor(kGray);
      marker_US[par_temp[0][imark3+94]]->Draw();
    }
  }


  for(Int_t imark4=0; imark4<17; imark4++){
    if((ADC_High_SFT_corr[par_temp[1][imark4+47]]!=0) && (has_TDC_SFT_hit[par_temp[1][imark4+47]])){
      marker_DS[par_temp[0][imark4+47]]->SetMarkerColor(1);
      marker_DS[par_temp[0][imark4+47]]->Draw();
    }

    if((ADC_High_SFT_corr[par_temp[1][imark4+111]]!=0) && (has_TDC_SFT_hit[par_temp[1][imark4+111]])){
      marker_US[par_temp[0][imark4+111]]->SetMarkerColor(1);
      marker_US[par_temp[0][imark4+111]]->Draw();
    }
    
    if(ADC_High_SFT[imark4+47]<0 && ADC_Low_SFT[imark4+47]>0 && Switch==1){
      marker_DS[par_temp[0][imark4+47]]->SetMarkerColor(kGray);
      marker_DS[par_temp[0][imark4+47]]->Draw();
    }
    
    if(ADC_High_SFT[imark4+111]<0 && ADC_Low_SFT[imark4+111]>0 && Switch==1){
      marker_US[par_temp[0][imark4+111]]->SetMarkerColor(kGray);
      marker_US[par_temp[0][imark4+111]]->Draw();
    }
  }
  

  c2->cd(5);
  
  /// uncomment for sft graph

  h_ADC_L1_DS->SetMaximum(1.1*ADC_High_corr_max);
  h_ADC_L1_DS->Draw();      h_ADC_L1_US->Draw("same");
  h_ADC_L2_DS->Draw("same");    h_ADC_L2_US->Draw("same");
  h_ADC_L3_DS->Draw("same");    h_ADC_L3_US->Draw("same");
  h_ADC_L4_DS->Draw("same");    h_ADC_L4_US->Draw("same");

  // Split LineADC_High_corr_max
  split_line = new TLine(64,0,64,1.1*ADC_High_corr_max);
  split_line->SetLineWidth(4);  
  split_line->SetLineColor(4);  
  split_line->SetLineStyle(2);
  split_line->Draw();

  ///
  
  
  ///////////////////////////////////////////////////////////////////////

  // Fill K Stop bars
  for(int i = 0; i<256; i++){
    if(k_stop_bar[i])
      h_TDC_Gap_Fibers_kaon->Fill(Xloc[i],Yloc[i]);
  }       
 
  /////
  ///////////////////////////////////////////////////////// K STOP CENTROID END
  
  cout << "  " << endl;


  TF1 *fit_line_ADC4;
  TF1 *fit_line_ADCA;
  
  
        // TARGET ROTATION (90 deg.)
        
  // ROTATE_CHANGE: any line tagged with ROTATE_CHANGE is used when tof1 is 6 or 12.
  // ROTATE_CHANGE
  // ROTATE_CHANGE    
        
  const int gap_to_fit_rotate = gap_to_fit;
         
  if((gap_to_fit==12 || gap_to_fit==6) && Rotate==1){    
    int ADC_High_TARGET_temp[256];
    int ADC_Low_TARGET_temp[256];
    int adc_low_target_temp[256];
    bool has_TDC_hit_temp[256];
    bool has_TDC_hit_Kstop_temp[256];
    int HG_TARGET_ADC_Thr_temp[256];
    int LG_TARGET_ADC_Thr_temp[256];
    int TARGET_High_has_neighbours_temp[256];
    bool k_stop_bar_temp[256];
  

    if(gap_to_fit==11)
      gap_to_fit = 2;
    else if(gap_to_fit==12)
      gap_to_fit = 3;
    else if(gap_to_fit==1)
      gap_to_fit = 4;
    else if(gap_to_fit==5)
      gap_to_fit = 8;
    else if(gap_to_fit==6)
      gap_to_fit = 9;
    else if(gap_to_fit==7)
      gap_to_fit = 10;
    

  
  h_TDC_Gap_Fibers->Reset();
  h_TDC_Gap_Fibers->Reset();
  h_target_TDC->Reset();
  h_target_TDC_copy->Reset();
  h_TDC_selected->Reset();
  h_TDC_selected2->Reset();
  h_TARGET_LG_Blue->Reset();

    
  for(int i = 0; i<256; i++){
    ADC_High_TARGET_temp[i] = ADC_High_TARGET[i];
    ADC_Low_TARGET_temp[i] = ADC_Low_TARGET[i];
    HG_TARGET_ADC_Thr_temp[i] = HG_TARGET_ADC_Thr[i];
    LG_TARGET_ADC_Thr_temp[i] = LG_TARGET_ADC_Thr[i];
    has_TDC_hit_temp[i] = has_TDC_hit[i];
    has_TDC_hit_Kstop_temp[i] = has_TDC_hit_Kstop[i];
    adc_low_target_temp[i] = adc_low_target[i];
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
      if(ADC_High_TARGET[q]>max_ADC) {
        max_index = q;
        max_ADC = ADC_Low_TARGET[q];
      }
    }

    for(Int_t q=0; q<256; q++){
      if (q == max_index) continue;   
      else {
        if(ADC_High_TARGET[q]>max_ADC2) {
          max_index2 = q;
          max_ADC2 = ADC_Low_TARGET[q];
        }
      }
    }

    for(Int_t q=0; q<256; q++){
      if ((q == max_index) || (q == max_index2)) continue;    
      else {
        if(ADC_High_TARGET[q]>max_ADC3) {
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
  
  
    for(Int_t j=0; j<256; j++){
      hyp[j] = sqrt(pow(x_cent - Xloc[j],2) + pow(y_cent - Yloc[j],2));
    }
    
    count = 0;
    hit_count = 0;
    
    for(Int_t i=0; i<256; i++){
      if(TARGET_High_has_neighbours[i]){ 
        if(ADC_High_TARGET[i]>=Angle_ADC_cut && has_TDC_hit[i]){
          h_target_TDC->Fill(Xloc[i],Yloc[i]);
          h_target_TDC_copy->Fill(Xloc[i],Yloc[i]);
          
          if(!k_stop_bar[i]){
            h_TDC_Gap_Fibers->Fill(Xloc[i],Yloc[i]);

          }
            
          count++;
          
          
          if(IsIn(i,channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],   // Additional weight for fibers close to the edge if hit
                    channel[gap_to_fit-1][2], channel[gap_to_fit-1][3], 
                    channel[gap_to_fit-1][4], channel[gap_to_fit-1][5], 
                    channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){
            
            if(!k_stop_bar[i]){     
              h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
              h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
            }
        
            h_TDC_selected2->Fill(Xloc[i], Yloc[i]);
            h_TDC_selected2->Fill(Xloc[i], Yloc[i]);
        
            h_TDC_selected->Fill(Xloc[i], Yloc[i]);
            h_TDC_selected->Fill(Xloc[i], Yloc[i]);
          }
        }

        if(ADC_High_TARGET[i]>=Angle_ADC_cut){
          h_target_ADC->Fill(Xloc[i],Yloc[i]);
          x_inc += Xloc[i];
          y_inc += Yloc[i];
          hit_count++;
        }
    
        if(ADC_High_TARGET[i]<0 && ADC_Low_TARGET[i]>=0 && Switch==1){
          if(!k_stop_bar[i])  
            h_TDC_Gap_Fibers->Fill(Xloc[i],Yloc[i]);{
          }
        h_TARGET_LG_Blue->Fill(Xloc[i],Yloc[i]);
        count++;
        }
      }
    }
  }


  /// Data counters
  int has_data_TDC2 = 0;
  int has_data_ADC2 = 0;
  int has_data_ADC3 = 0;
  int has_data_ADC4 = 0;
  int has_data_ADCA = 0;

  h_target_ADC4->Fill(Xloc[max_index],Yloc[max_index]);
  h_target_ADC4->Fill(Xloc[max_index2],Yloc[max_index2]);
  h_target_ADC4->Fill(Xloc[max_index3],Yloc[max_index3]);
  h_target_ADC4->Fill(Xloc[max_index4],Yloc[max_index4]);

  ///// Select edge fiber for track fitting


  double Xloc_gap = 0;
  double Yloc_gap = 0;

  if (gap_to_fit == 0) {
    for (int j=0; j<8; j++) {
      if ((ADC_High_TARGET[channel[0][j]]>0) && (has_TDC_hit[channel[0][j]])) {
        Xloc_gap = Xloc[channel[0][j]];
        Yloc_gap = Yloc[channel[0][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[0][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[0][3]];
  }
   if (gap_to_fit == 1) {
    for (int j=0; j<8; j++) {
      if ((ADC_High_TARGET[channel[1][j]]>0) && (has_TDC_hit[channel[1][j]])) {
        Xloc_gap = Xloc[channel[1][j]];
        Yloc_gap = Yloc[channel[1][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[1][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[1][3]];
  }
   if (gap_to_fit == 2) {
    for (int j=0; j<8; j++) {
      if ((ADC_High_TARGET[channel[2][j]]>0) && (has_TDC_hit[channel[2][j]])) {
        Xloc_gap = Xloc[channel[2][j]];
        Yloc_gap = Yloc[channel[2][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[2][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[2][3]];
  }
   if (gap_to_fit == 3) {
    for (int j=0; j<8; j++) {
      if ((ADC_High_TARGET[channel[3][j]]>0) && (has_TDC_hit[channel[3][j]])) {
        Xloc_gap = Xloc[channel[3][j]];
        Yloc_gap = Yloc[channel[3][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[3][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[3][3]];
  }
   if (gap_to_fit == 4) {
    for (int j=0; j<8; j++) {
      if ((ADC_High_TARGET[channel[4][j]]>0) && (has_TDC_hit[channel[4][j]])) {
        Xloc_gap = Xloc[channel[4][j]];
        Yloc_gap = Yloc[channel[4][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[4][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[4][3]];
  }
   if (gap_to_fit == 5) {
    for (int j=0; j<8; j++) {
      if ((ADC_High_TARGET[channel[5][j]]>0) && (has_TDC_hit[channel[5][j]])) {
        Xloc_gap = Xloc[channel[5][j]];
        Yloc_gap = Yloc[channel[5][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[5][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[5][3]];
  }
   if (gap_to_fit == 6) {
    for (int j=0; j<8; j++) {
      if ((ADC_High_TARGET[channel[6][j]]>0) && (has_TDC_hit[channel[6][j]])) {
        Xloc_gap = Xloc[channel[6][j]];
        Yloc_gap = Yloc[channel[6][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[6][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[6][3]];
  }
   if (gap_to_fit == 7) {
    for (int j=0; j<8; j++) {
      if ((ADC_High_TARGET[channel[7][j]]>0) && (has_TDC_hit[channel[7][j]])) {
        Xloc_gap = Xloc[channel[7][j]];
        Yloc_gap = Yloc[channel[7][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[7][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[7][3]];
  }
   if (gap_to_fit == 8) {
    for (int j=0; j<8; j++) {
      if ((ADC_High_TARGET[channel[8][j]]>0) && (has_TDC_hit[channel[8][j]])) {
        Xloc_gap = Xloc[channel[8][j]];
        Yloc_gap = Yloc[channel[8][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[8][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[8][3]];
  }
   if (gap_to_fit == 9) {
    for (int j=0; j<8; j++) {
      if ((ADC_High_TARGET[channel[9][j]]>0) && (has_TDC_hit[channel[9][j]])) {
        Xloc_gap = Xloc[channel[9][j]];
        Yloc_gap = Yloc[channel[9][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[9][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[9][3]];
  }
   if (gap_to_fit == 10) {
    for (int j=0; j<8; j++) {
      if ((ADC_High_TARGET[channel[10][j]]>0) && (has_TDC_hit[channel[10][j]])) {
        Xloc_gap = Xloc[channel[10][j]];
        Yloc_gap = Yloc[channel[10][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[10][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[10][3]];
  }
   if (gap_to_fit == 11) {
    for (int j=0; j<8; j++) {
      if ((ADC_High_TARGET[channel[11][j]]>0) && (has_TDC_hit[channel[11][j]])) {
        Xloc_gap = Xloc[channel[11][j]];
        Yloc_gap = Yloc[channel[11][j]];
      }
    }
    if (Xloc_gap == 0) Xloc_gap = Xloc[channel[11][3]];
    if (Yloc_gap == 0) Yloc_gap = Yloc[channel[11][3]];
  }

  /// Determine which gap point is closest to edge fiber

  double xdistance1 = pow((TOF_Xloc[(gap_to_fit*3)]-Xloc_gap),2);
  double ydistance1 = pow((TOF_Yloc[(gap_to_fit*3)]-Yloc_gap),2);

  double xdistance2 = pow((TOF_Xloc[(gap_to_fit*3)+1]-Xloc_gap),2);
  double ydistance2 = pow((TOF_Yloc[(gap_to_fit*3)+1]-Yloc_gap),2);

  double xdistance3 = pow((TOF_Xloc[(gap_to_fit*3)+2]-Xloc_gap),2);
  double ydistance3 = pow((TOF_Yloc[(gap_to_fit*3)+2]-Yloc_gap),2);

  double xhyp1 = double(sqrt(double(xdistance1) + double(ydistance1)));
  double xhyp2 = double(sqrt(double(xdistance2) + double(ydistance2)));
  double xhyp3 = double(sqrt(double(xdistance3) + double(ydistance3)));


  int closest_gap_point;

  if ((xhyp1 <= xhyp2) && (xhyp1 <= xhyp3)) {
    closest_gap_point = 0;
  }
  else if (xhyp2 <= xhyp3) {
    closest_gap_point = 1;
  }
  else {
    closest_gap_point = 2;
  }

  h_target_ADC3->Fill(TOF_Xloc[(gap_to_fit*3) + closest_gap_point],TOF_Yloc[(gap_to_fit*3) + closest_gap_point]);

  ////Fill TARGET tracking histograms with data

  for(Int_t i=0; i<256; i++){

    if ((i == max_index) || (i == max_index2) || (i == max_index3) || (i == max_index4)) continue;

    if((has_TDC_hit[i]) && (ADC_High_TARGET[i]>=Angle_ADC_cut) && (hyp[i] > 12.5)) {
      h_target_ADC2->Fill(Xloc[i],Yloc[i]);
      has_data_ADC2++;
    }
    if((has_TDC_hit[i]) && (ADC_High_TARGET[i]>=Angle_ADC_cut) && (adc_low_target[i]-ADC_cut_TARGET2>0)) {
      h_target_TDC2->Fill(Xloc[i],Yloc[i]);
      has_data_TDC2++;
    }
    if((has_TDC_hit[i]) && (ADC_High_TARGET[i]>=Angle_ADC_cut) && (adc_low_target[i]-ADC_cut_TARGET2>0)) {
      h_target_ADC3->Fill(Xloc[i],Yloc[i]);
      has_data_ADC3++;
    }
    if((has_TDC_hit[i]) && (ADC_High_TARGET[i]>=Angle_ADC_cut) && (hyp[i] < 12.5)) {
      h_target_ADC4->Fill(Xloc[i],Yloc[i]);
      has_data_ADC4++;
    }
  }


  if (display_fits ==1){
    h_target_ADC2->SetMarkerStyle(25);
    h_target_ADC2->SetMarkerColor(4);

    h_target_TDC2->SetMarkerStyle(25);
    h_target_TDC2->SetMarkerColor(4);

    cout << "  " << endl;

    if (has_data_TDC2 > 1){
      phi_TDC=(180/PI)*atan(par1_TDC);
      sprintf(angle_string_TDC,"Phi = %3.4f deg.",phi_TDC);
      tex_angle_TDC=new TLatex(0,-9,angle_string_TDC);
      tex_angle_TDC->SetTextSize(0.05);
    }
    else {
      cout << "Histo Fit 3 Is Empty" << endl;
    }

    h_target_ADC4->SetMarkerStyle(25);
    h_target_ADC4->SetMarkerColor(4);
  }  /// TO CHECK !!

  //// Used for primary track fitting

  h_target_ADC3->SetMarkerStyle(25);
  h_target_ADC3->SetMarkerColor(4);

  if (has_data_ADC3 > 1){
    double xcoord = 0;
    Int_t unique_x = 0;


    if (has_data_ADCA > 2) {


      h_target_ADCA->Fill(TOF_Xloc[(gap_to_fit*3) + closest_gap_point],TOF_Yloc[(gap_to_fit*3) + closest_gap_point]);
      
      /// Refit histogram with filtered data + closest gap point

      h_target_ADCA->SetMarkerStyle(25);
      h_target_ADCA->SetMarkerColor(4);
      h_target_ADCA->Draw();
      h_target_ADCA->Fit("pol1", "QCM");
      fit_line_ADCA = h_target_ADCA->GetFunction("pol1");
      fit_line_ADCA->SetLineWidth(2);
      fit_line_ADCA->SetLineColor(2);

      
      /// Retrieve slope & intercept for refit line
      par0_ADC=fit_line_ADCA->GetParameter(0);
      par1_ADC=fit_line_ADCA->GetParameter(1);



      /// Ring intercept coordinates
      float determinant = 4*(pow(par0_ADC,2))*(pow(par1_ADC,2)) - 4*(pow(par1_ADC,2) + 1)*(pow(par0_ADC,2)-1600);
      float x_circle_int1 = (-2*(par0_ADC)*(par1_ADC) + sqrt(determinant))/(2*(pow(par1_ADC,2) + 1));
      float y_circle_int1 = (par1_ADC)*(x_circle_int1) + par0_ADC;

      float x_circle_int2 = (-2*(par0_ADC)*(par1_ADC) - sqrt(determinant))/(2*(pow(par1_ADC,2) + 1));
      float y_circle_int2 = (par1_ADC)*(x_circle_int2) + par0_ADC;

      if (unique_x == 1) {
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

      if (SFTxhyp1 < SFTxhyp2) {
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

      if (SFT_x_intercept >0) {
        SFT_phi = (atan2(SFT_x_intercept,SFT_y_intercept))*57.2957795;
      }
      else {
        SFT_phi = 360 + (atan2(SFT_x_intercept,SFT_y_intercept))*57.2957795;
      }

      sprintf(angle_string_ADC,"Phi = %3.4f deg.",SFT_phi);
      tex_angle_ADC=new TLatex(0,-9,angle_string_ADC);
      tex_angle_ADC->SetTextSize(0.05);
      tex_angle_ADC->Draw();

      cout << "" << endl;
      cout << "Angle between SFT Layer Ring Intercept and Centre of TARGET: " << SFT_phi << " deg." << endl;
      cout << "" << endl;

      if (has_data_ADC4 > 1) {

        h_target_ADC4->Fit("pol1", "Q0CM");
        fit_line_ADC4 = h_target_ADC4->GetFunction("pol1");
        fit_line_ADC4->SetLineWidth(2);
        fit_line_ADC4->SetLineColor(1);

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
  else {
    cout << "Histo Fit 4 Is Empty" << endl;
  }

  cout << "//////  TOF1 SCORING  //////" << endl;
  for (int k=0; k<12; k++) {
    cout << "TOF1 score " << k+1 << ": " << gap_counter[k] << endl;
  }

  cout << "" << endl;  
  for (int k=0; k<12; k++) {
  }

  cout << "TOF1 SELECTED FOR FITTING:  TOF1 " << gap_to_fit_rotate << endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    TF1 *fit_line_GoodLG;           TF1 *fit_line_GoodLG_weighed;
    TF1 *fit_line_TDC_selected;     TF1 *fit_line_TDC_selected2;
    TF1 *fit_line_TDC_Gap_Fibers;



    float Gap[12][3][2] = {{{0}}};
    
    float a_fit_TDC_selected=0.;               float b_fit_TDC_selected=0.;
    float a_fit_GoodLG=0.;                     float b_fit_GoodLG=0.;
    float a_fit_GoodLG_weighted=0.;            float b_fit_GoodLG_weighted=0.;
    float a_fit_TDC_selected_weighted=0.;      float b_fit_TDC_selected_weighted=0.;
    float a_fit_TDC_Gap_Fibers=0.;             float b_fit_TDC_Gap_Fibers=0.;



    for(int g=0; g<12; g++){
      Gap[g][0][0] = TOF_Xloc[3*g];
      Gap[g][1][0] = TOF_Xloc[3*g+1];
      Gap[g][2][0] = TOF_Xloc[3*g+2];
      
      Gap[g][0][1] = TOF_Yloc[3*g];
      Gap[g][1][1] = TOF_Yloc[3*g+1];
      Gap[g][2][1] = TOF_Yloc[3*g+2];
    }

   
    for(int l=0; l<256; l++){
      if(adc_low_target[l]>=LG_TARGET_ADC_Thr[l] && has_TDC_hit[l])
      {
        h_GoodLG->Fill(Xloc[l], Yloc[l]);
        h_GoodLG_copy->Fill(Xloc[l], Yloc[l]);
        h_GoodLG_weighted->Fill(Xloc[l], Yloc[l]);
      }
    }

    for(int i=0; i<12; i++){
      h_Circle->Fill(Gap[i][0][0], Gap[i][0][1]);
      h_Circle->Fill(Gap[i][1][0], Gap[i][1][1]);
      h_Circle->Fill(Gap[i][2][0], Gap[i][2][1]);
  
    }

    h_Target_Center->Fill(0., 0.);

    for(int t=0; t<256; t++){
      h_Target->Fill(Xloc[t],Yloc[t]);
    }

    h_max->Fill(Xloc[max_index], Yloc[max_index]);        h_max_copy->Fill(Xloc[max_index], Yloc[max_index]);
   
    h_kaon->Fill(Xloc[max_index], Yloc[max_index]);       h_kaon_copy->Fill(Xloc[max_index], Yloc[max_index]);
    h_kaon->Fill(Xloc[max_index2], Yloc[max_index2]);     h_kaon_copy->Fill(Xloc[max_index2], Yloc[max_index2]);
    h_kaon->Fill(Xloc[max_index3], Yloc[max_index3]);     h_kaon_copy->Fill(Xloc[max_index3], Yloc[max_index3]);
    h_kaon->Fill(Xloc[max_index4], Yloc[max_index4]);     h_kaon_copy->Fill(Xloc[max_index4], Yloc[max_index4]);

  // ROTATE_CHANGE
  if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6) && Rotate==1){  
    for(int g=0; g<3; g++){
          h_TOF1_rotate->Fill(Gap[gap_to_fit_rotate-1][g][0], Gap[gap_to_fit_rotate-1][g][1]);
          h_TDC_selected->Fill(Gap[gap_to_fit-1][g][0], Gap[gap_to_fit-1][g][1]);  
        } 
  }
  else{
    for(int g=0; g<3; g++){
        h_TOF1->Fill(Gap[gap_to_fit-1][g][0], Gap[gap_to_fit-1][g][1]);
        h_TDC_selected->Fill(Gap[gap_to_fit-1][g][0], Gap[gap_to_fit-1][g][1]);  
      } 
  }

    h_Circle->SetMarkerStyle(5);
    h_Circle->SetMarkerSize(1.2);
    h_Circle->SetLineWidth(2);

    h_Target->SetMarkerStyle(25);         h_Target_Center->SetMarkerStyle(5);
    h_Target->SetMarkerColor(1);          h_Target_Center->SetMarkerColor(1);
    h_Target->SetMarkerSize(0.8);         h_Target_Center->SetMarkerSize(1);

    h_GoodLG->SetMarkerStyle(21);         h_GoodLG_copy->SetMarkerStyle(21);          h_GoodLG_weighted->SetMarkerStyle(21);
    h_GoodLG->SetMarkerColor(2);          h_GoodLG_copy->SetMarkerColor(2);           h_GoodLG_weighted->SetMarkerColor(2);
    h_GoodLG->SetMarkerSize(0.8);         h_GoodLG_copy->SetMarkerSize(0.8);          h_GoodLG_weighted->SetMarkerSize(0.8);

    h_target_TDC->SetMarkerStyle(25);     h_target_TDC_copy->SetMarkerStyle(25);      h_TDC_selected->SetMarkerStyle(25);   h_TDC_selected2->SetMarkerStyle(25);
    h_target_TDC->SetMarkerColor(2);      h_target_TDC_copy->SetMarkerColor(2);       h_TDC_selected->SetMarkerColor(2);;   h_TDC_selected2->SetMarkerColor(2);
    h_target_TDC->SetMarkerSize(0.8);     h_target_TDC_copy->SetMarkerSize(0.8);      h_TDC_selected->SetMarkerSize(0.8);   h_TDC_selected2->SetMarkerSize(0.8);

    h_kaon->SetMarkerStyle(21);           h_kaon_copy->SetMarkerStyle(21);            h_TDC_Gap_Fibers->SetMarkerStyle(25);
    h_kaon->SetMarkerColor(1);            h_kaon_copy->SetMarkerColor(1);             h_TDC_Gap_Fibers->SetMarkerColor(2);
    h_kaon->SetMarkerSize(0.8);           h_kaon_copy->SetMarkerSize(0.8);            h_TDC_Gap_Fibers->SetMarkerSize(0.8);
    
    h_TDC_Gap_Fibers_kaon->SetMarkerStyle(21);
    h_TDC_Gap_Fibers_kaon->SetMarkerColor(kBlue-6);
    h_TDC_Gap_Fibers_kaon->SetMarkerSize(0.8);
    
    h_max->SetMarkerStyle(21);            h_max_copy->SetMarkerStyle(21);             h_TARGET_LG_Blue->SetMarkerStyle(25);
    h_max->SetMarkerColor(kViolet+1);     h_max_copy->SetMarkerColor(kViolet+1);      h_TARGET_LG_Blue->SetMarkerColor(4);
    h_max->SetMarkerSize(1);              h_max_copy->SetMarkerSize(1);               h_TARGET_LG_Blue->SetMarkerSize(0.8);

    h_TOF1->SetMarkerStyle(20);           h_TOF1_closest->SetMarkerStyle(20);         h_TDC_Gap_Fibers_copy->SetMarkerStyle(25); // ROTATE_CHANGE   
    h_TOF1->SetMarkerColor(2);            h_TOF1_closest->SetMarkerColor(3);          h_TDC_Gap_Fibers_copy->SetMarkerColor(2); // ROTATE_CHANGE
    h_TOF1->SetMarkerSize(1.5);           h_TOF1_closest->SetMarkerSize(1.5);         h_TDC_Gap_Fibers_copy->SetMarkerSize(0.8); // ROTATE_CHANGE
  
    h_TOF1_rotate->SetMarkerStyle(20);    h_TOF1_closest_rotate->SetMarkerStyle(20);  // ROTATE_CHANGE
    h_TOF1_rotate->SetMarkerColor(2);     h_TOF1_closest_rotate->SetMarkerColor(3);   // ROTATE_CHANGE
    h_TOF1_rotate->SetMarkerSize(1.5);    h_TOF1_closest_rotate->SetMarkerSize(1.5);  // ROTATE_CHANGE
    
    h_Centroid->SetMarkerStyle(5);
    h_Centroid->SetMarkerColor(1);         
    h_Centroid->SetMarkerSize(4); 
    
  
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




    c2->cd(6);
    h_TDC_selected->Draw();
    h_TDC_selected->Fit("pol1", "QCM");
    fit_line_TDC_selected = h_TDC_selected->GetFunction("pol1");
    fit_line_TDC_selected->SetLineWidth(2);
    fit_line_TDC_selected->SetLineColor(1);
    a_fit_TDC_selected=fit_line_TDC_selected->GetParameter(1);      
    b_fit_TDC_selected=fit_line_TDC_selected->GetParameter(0);
    HorizontalAxis->Draw("same");

    if(h_GoodLG->GetEntries()!=0){
      h_GoodLG->Fit("pol1", "QCM");
      fit_line_GoodLG = h_GoodLG->GetFunction("pol1");
      fit_line_GoodLG->SetLineWidth(2);
      fit_line_GoodLG->SetLineColor(1);
      a_fit_GoodLG=fit_line_GoodLG->GetParameter(1);      
      b_fit_GoodLG=fit_line_GoodLG->GetParameter(0);
    }
    else{
    cout << "Empty Histogram" << endl;
    }
  


  // Add weight to center point of selected TOF1 (gap_to_fit_rotate)
  for(int i = 0; i<3; i++){
      h_TDC_Gap_Fibers->Fill(Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
  }
  

  // Fit 1: Targets hit and center point of selected TOF1
  h_TDC_Gap_Fibers->Fit("pol1", "QCM");
  fit_line_TDC_Gap_Fibers = h_TDC_Gap_Fibers->GetFunction("pol1");
  fit_line_TDC_Gap_Fibers->SetLineWidth(2);
  fit_line_TDC_Gap_Fibers->SetLineColor(2);
  a_fit_TDC_Gap_Fibers=h_TDC_Gap_Fibers->GetFunction("pol1")->GetParameter(1);      
  b_fit_TDC_Gap_Fibers=h_TDC_Gap_Fibers->GetFunction("pol1")->GetParameter(0);

  // ReFill histogram with points closer than max_dist
  h_TDC_Gap_Fibers->Reset();
  h_TDC_Gap_Fibers_copy->Reset();
  for(int i = 0; i<256; i++){   
    if(distance_to_line(Xloc[i],Yloc[i],a_fit_TDC_Gap_Fibers,b_fit_TDC_Gap_Fibers) <= max_dist && TARGET_High_has_neighbours[i] && !k_stop_bar[i]){   
      if(ADC_High_TARGET[i]>=Angle_ADC_cut && has_TDC_hit[i]){
      h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
      h_TDC_Gap_Fibers_copy->Fill(Xloc[TARGET_Rotated_index[i]],Yloc[TARGET_Rotated_index[i]]); // ROTATE_CHANGE
      
      
        if(IsIn(i,channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],   
                  channel[gap_to_fit-1][2], channel[gap_to_fit-1][3], 
                  channel[gap_to_fit-1][4], channel[gap_to_fit-1][5], 
                  channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){
        h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);     
        h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
        h_TDC_Gap_Fibers_copy->Fill(Xloc[TARGET_Rotated_index[i]],Yloc[TARGET_Rotated_index[i]]); // ROTATE_CHANGE
        h_TDC_Gap_Fibers_copy->Fill(Xloc[TARGET_Rotated_index[i]],Yloc[TARGET_Rotated_index[i]]); // ROTATE_CHANGE      
        }
      }
      if(ADC_High_TARGET[i]<0 && ADC_Low_TARGET[i]>=0 && Switch==1){
        h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
        h_TDC_Gap_Fibers_copy->Fill(Xloc[TARGET_Rotated_index[i]],Yloc[TARGET_Rotated_index[i]]);  // ROTATE_CHANGE
      }
    } 
  }
  
 if(h_TDC_Gap_Fibers->GetEntries() != 0){
    //Fit 2: Targets hit within max_distance of the best fit line from Fit 1
    h_TDC_Gap_Fibers->Fit("pol1", "QCM");
    fit_line_TDC_Gap_Fibers = h_TDC_Gap_Fibers->GetFunction("pol1");
    fit_line_TDC_Gap_Fibers->SetLineWidth(2);
    fit_line_TDC_Gap_Fibers->SetLineColor(2);
    a_fit_TDC_Gap_Fibers=fit_line_TDC_Gap_Fibers->GetParameter(1);      
    b_fit_TDC_Gap_Fibers=fit_line_TDC_Gap_Fibers->GetParameter(0);
    }
  else
    cout << "Empty FIT 2" << endl;
  

    float x_int_TDC[2];                           float y_int_TDC[2];
    float x_int_TDC_selected[2];                  float y_int_TDC_selected[2];
    float x_int_TDC_selected_weighted[2];         float y_int_TDC_selected_weighted[2];
    float x_int_GoodLG[2];                        float y_int_GoodLG[2];
    float x_int_GoodLG_weighted[2];               float y_int_GoodLG_weighted[2];
    float x_int_TDC_Gap_Fibers[2];                float y_int_TDC_Gap_Fibers[2];
    float x_int_TDC_SFT[2];                       float y_int_TDC_SFT[2];
    float x_int_GoodLG_SFT[2];                    float y_int_GoodLG_SFT[2];
    float x_int_TDC_Gap_Fibers_SFT[2];            float y_int_TDC_Gap_Fibers_SFT[2];
    float x_int_TARGET[2];                        float y_int_TARGET[2];

    x_int_TDC_selected[0] = intersectx1(a_fit_TDC_selected, b_fit_TDC_selected, R_TOF1);
    x_int_TDC_selected[1] = intersectx2(a_fit_TDC_selected, b_fit_TDC_selected, R_TOF1);
    y_int_TDC_selected[0] = y1_int(x_int_TDC_selected[0], a_fit_TDC_selected, b_fit_TDC_selected);
    y_int_TDC_selected[1] = y2_int(x_int_TDC_selected[1], a_fit_TDC_selected, b_fit_TDC_selected);
  

    x_int_GoodLG[0] = intersectx1(a_fit_GoodLG, b_fit_GoodLG, R_TOF1);
    x_int_GoodLG[1] = intersectx2(a_fit_GoodLG, b_fit_GoodLG, R_TOF1);
    y_int_GoodLG[0] = y1_int(x_int_GoodLG[0], a_fit_GoodLG, b_fit_GoodLG);
    y_int_GoodLG[1] = y2_int(x_int_GoodLG[1], a_fit_GoodLG, b_fit_GoodLG);


    x_int_TDC_Gap_Fibers[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);
    x_int_TDC_Gap_Fibers[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);   
    y_int_TDC_Gap_Fibers[0] = y1_int(x_int_TDC_Gap_Fibers[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers); 
    y_int_TDC_Gap_Fibers[1] = y2_int(x_int_TDC_Gap_Fibers[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);

    x_int_TDC_Gap_Fibers_SFT[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_SFT_L1);
    x_int_TDC_Gap_Fibers_SFT[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_SFT_L1);   
    y_int_TDC_Gap_Fibers_SFT[0] = y1_int(x_int_TDC_Gap_Fibers_SFT[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers); 
    y_int_TDC_Gap_Fibers_SFT[1] = y2_int(x_int_TDC_Gap_Fibers_SFT[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);


    /// Selection of the Good Intersect Coordinates
    //////////////////////////////////////////////////////////////////////
    
    float x_TDC_sel_intersect1=0.;   float y_TDC_sel_intersect1=0.;
    float x_GoodLG_intersect1=0.;    float y_GoodLG_intersect1=0.;
    float x_TDC_Gap_Fibers=0.;       float y_TDC_Gap_Fibers=0.;
    float x_TDC_Gap_Fibers_SFT=0.;   float y_TDC_Gap_Fibers_SFT=0.;
 
    float dist1_TDC_selected[2];
    float dist1_GoodLG[2];
    float dist1_TDC_Gap_Fibers[2];
    float dist1_TDC_Gap_Fibers_SFT[2];

    for(int i=0; i<2; i++){
      dist1_TDC_selected[i] = distance(x_int_TDC_selected[i], y_int_TDC_selected[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_GoodLG[i] = distance(x_int_GoodLG[i], y_int_GoodLG[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
    }    
    
    if(dist1_TDC_selected[0] < dist1_TDC_selected[1]){
      x_TDC_sel_intersect1 = x_int_TDC_selected[0];
      y_TDC_sel_intersect1 = y_int_TDC_selected[0];
    }
    else if(dist1_TDC_selected[1] < dist1_TDC_selected[0]){
      x_TDC_sel_intersect1 = x_int_TDC_selected[1];
      y_TDC_sel_intersect1 = y_int_TDC_selected[1];
    }
    else cout << "ERROR !" << endl;


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

    if(dist1_TDC_Gap_Fibers_SFT[0] < dist1_TDC_Gap_Fibers_SFT[1]){
      x_TDC_Gap_Fibers_SFT = x_int_TDC_Gap_Fibers_SFT[0];
      y_TDC_Gap_Fibers_SFT = y_int_TDC_Gap_Fibers_SFT[0];
    }
    else if(dist1_TDC_Gap_Fibers_SFT[1] < dist1_TDC_Gap_Fibers_SFT[0]){
      x_TDC_Gap_Fibers_SFT = x_int_TDC_Gap_Fibers_SFT[1];
      y_TDC_Gap_Fibers_SFT = y_int_TDC_Gap_Fibers_SFT[1];
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
      dist2_TDC_selected[ii] = distance(x_TDC_Gap_Fibers, y_TDC_Gap_Fibers, Gap[gap_to_fit-1][ii][0], Gap[gap_to_fit-1][ii][1]);
      dist2_GoodLG[ii] = distance(x_GoodLG_intersect1, y_GoodLG_intersect1, Gap[gap_to_fit-1][ii][0], Gap[gap_to_fit-1][ii][1]);

      if(dist2_TDC_selected[ii] <= dist2_TDC_selected_min){
        dist2_TDC_selected_min = dist2_TDC_selected[ii];
        selected_TDC_selected = ii;
      }

      if(dist2_GoodLG[ii] <= dist2_GoodLG_min){
        dist2_GoodLG_min = dist2_GoodLG[ii];
      }
    }
  

    h_TOF1_closest->Fill(Gap[gap_to_fit_rotate-1][selected_TDC_selected][0], Gap[gap_to_fit_rotate-1][selected_TDC_selected][1]);
  
    

    //////////////////////////////////////////////////////////////////////




  h_int_TDC->Fill(x_int_TDC[0],y_int_TDC[0]);
  h_int_TDC->Fill(x_int_TDC[1],y_int_TDC[1]);


  h_int_TDC_selected->Fill(x_TDC_sel_intersect1, y_TDC_sel_intersect1);

  h_int_GoodLG->Fill(x_GoodLG_intersect1, y_GoodLG_intersect1);

    /// Add Weight on TOF1
    //////////////////////////////////////////////////////////////////
  for(int i = 0; i<5; i++)
      h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected][0], Gap[gap_to_fit-1][selected_TDC_selected][1]);
   
  for(int i = 0; i<5; i++)
    h_GoodLG_weighted->Fill(Gap[gap_to_fit-1][selected_TDC_selected][0], Gap[gap_to_fit-1][selected_TDC_selected][1]);
  
  for(int i = 0; i<3; i++)
    h_TDC_Gap_Fibers->Fill(Gap[gap_to_fit-1][selected_TDC_selected][0], Gap[gap_to_fit-1][selected_TDC_selected][1]);
    //////////////////////////////////////////////////////////////////
 

    h_int_TDC->SetMarkerStyle(20);            h_int_GoodLG->SetMarkerStyle(20);      h_int_GoodLG_weighted->SetMarkerStyle(20);
    h_int_TDC->SetMarkerColor(4);             h_int_GoodLG->SetMarkerColor(4);       h_int_GoodLG_weighted->SetMarkerColor(4);
    h_int_TDC->SetMarkerSize(0.8);            h_int_GoodLG->SetMarkerSize(0.8);      h_int_GoodLG_weighted->SetMarkerSize(0.8);

    h_int_TDC_selected->SetMarkerStyle(20);   h_int_TDC_selected_weighted->SetMarkerStyle(20);
    h_int_TDC_selected->SetMarkerColor(4);    h_int_TDC_selected_weighted->SetMarkerColor(4);
    h_int_TDC_selected->SetMarkerSize(0.8);   h_int_TDC_selected_weighted->SetMarkerSize(0.8);

    h_int_TDC_Gap_Fibers->SetMarkerStyle(20);   h_int_TDC_Gap_Fibers_SFT->SetMarkerStyle(20);
    h_int_TDC_Gap_Fibers->SetMarkerColor(4);    h_int_TDC_Gap_Fibers_SFT->SetMarkerColor(4);
    h_int_TDC_Gap_Fibers->SetMarkerSize(0.8);   h_int_TDC_Gap_Fibers_SFT->SetMarkerSize(0.8);
    
    h_int_TDC_Gap_Fibers_rotate->SetMarkerStyle(20);   h_int_TDC_Gap_Fibers_SFT_rotate->SetMarkerStyle(20);
    h_int_TDC_Gap_Fibers_rotate->SetMarkerColor(4);    h_int_TDC_Gap_Fibers_SFT_rotate->SetMarkerColor(4);
    h_int_TDC_Gap_Fibers_rotate->SetMarkerSize(0.8);   h_int_TDC_Gap_Fibers_SFT_rotate->SetMarkerSize(0.8);

    h_int_TDC_SFT->SetMarkerStyle(20);        h_int_GoodLG_SFT->SetMarkerStyle(20);
    h_int_TDC_SFT->SetMarkerColor(1);         h_int_GoodLG_SFT->SetMarkerColor(1);
    h_int_TDC_SFT->SetMarkerSize(0.8);        h_int_GoodLG_SFT->SetMarkerSize(0.8);

    h_int_TDC_SFT->SetMarkerStyle(20);
    h_int_TDC_SFT->SetMarkerColor(1);
    h_int_TDC_SFT->SetMarkerSize(0.8);

    double ParError = 999.99;

    h_TDC_selected2->Fit("pol1", "QCM");
    fit_line_TDC_selected2 = h_TDC_selected2->GetFunction("pol1");
    fit_line_TDC_selected2->SetLineWidth(2);
    fit_line_TDC_selected2->SetLineColor(2);
    a_fit_TDC_selected_weighted=fit_line_TDC_selected2->GetParameter(1);      
    b_fit_TDC_selected_weighted=fit_line_TDC_selected2->GetParameter(0);

    h_GoodLG_weighted->Fit("pol1", "QCM");
    fit_line_GoodLG_weighed = h_GoodLG_weighted->GetFunction("pol1");
    fit_line_GoodLG_weighed->SetLineWidth(2);
    fit_line_GoodLG_weighed->SetLineColor(1);
    a_fit_GoodLG_weighted=fit_line_GoodLG_weighed->GetParameter(1);      
    b_fit_GoodLG_weighted=fit_line_GoodLG_weighed->GetParameter(0);

    h_TDC_Gap_Fibers->Fit("pol1", "QCM");
    fit_line_TDC_Gap_Fibers = h_TDC_Gap_Fibers->GetFunction("pol1");
    fit_line_TDC_Gap_Fibers->SetLineWidth(2);
    fit_line_TDC_Gap_Fibers->SetLineColor(2);
    a_fit_TDC_Gap_Fibers=fit_line_TDC_Gap_Fibers->GetParameter(1);      
    b_fit_TDC_Gap_Fibers=fit_line_TDC_Gap_Fibers->GetParameter(0);
    ParError=fit_line_TDC_Gap_Fibers->GetParError(1);



    x_int_TDC_selected_weighted[0] = intersectx1(a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted, R_TOF1);
    x_int_TDC_selected_weighted[1] = intersectx2(a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted, R_TOF1); 
    x_int_TDC_SFT[0] = intersectx1(a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted, R_SFT_L1);
    x_int_TDC_SFT[1] = intersectx2(a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted, R_SFT_L1);

    y_int_TDC_selected_weighted[0] = y1_int(x_int_TDC_selected_weighted[0], a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted);
    y_int_TDC_selected_weighted[1] = y2_int(x_int_TDC_selected_weighted[1], a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted);
    y_int_TDC_SFT[0] = y1_int(x_int_TDC_SFT[0], a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted);
    y_int_TDC_SFT[1] = y2_int(x_int_TDC_SFT[1], a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted);

    x_int_GoodLG_weighted[0] = intersectx1(a_fit_GoodLG_weighted, b_fit_GoodLG_weighted, R_TOF1);
    x_int_GoodLG_weighted[1] = intersectx2(a_fit_GoodLG_weighted, b_fit_GoodLG_weighted, R_TOF1);
    x_int_GoodLG_SFT[0] = intersectx1(a_fit_GoodLG_weighted, b_fit_GoodLG_weighted, R_SFT_L1);
    x_int_GoodLG_SFT[1] = intersectx2(a_fit_GoodLG_weighted, b_fit_GoodLG_weighted, R_SFT_L1);

    y_int_GoodLG_weighted[0] = y1_int(x_int_GoodLG_weighted[0], a_fit_GoodLG_weighted, b_fit_GoodLG_weighted);
    y_int_GoodLG_weighted[1] = y2_int(x_int_GoodLG_weighted[1], a_fit_GoodLG_weighted, b_fit_GoodLG_weighted);
    y_int_GoodLG_SFT[0] = y1_int(x_int_GoodLG_SFT[0], a_fit_GoodLG_weighted, b_fit_GoodLG_weighted);
    y_int_GoodLG_SFT[1] = y2_int(x_int_GoodLG_SFT[1], a_fit_GoodLG_weighted, b_fit_GoodLG_weighted);

 
    x_int_TDC_Gap_Fibers[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);
    x_int_TDC_Gap_Fibers[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);   
    x_int_TDC_Gap_Fibers_SFT[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_SFT_L1);
    x_int_TDC_Gap_Fibers_SFT[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_SFT_L1);   

    y_int_TDC_Gap_Fibers[0] = y1_int(x_int_TDC_Gap_Fibers[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers); 
    y_int_TDC_Gap_Fibers[1] = y2_int(x_int_TDC_Gap_Fibers[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);
    y_int_TDC_Gap_Fibers_SFT[0] = y1_int(x_int_TDC_Gap_Fibers_SFT[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers); 
    y_int_TDC_Gap_Fibers_SFT[1] = y2_int(x_int_TDC_Gap_Fibers_SFT[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);

    x_int_TARGET[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TARGET);
    x_int_TARGET[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TARGET);
    y_int_TARGET[0] = y1_int(x_int_TARGET[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);
    y_int_TARGET[1] = y2_int(x_int_TARGET[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);


    float x_TDC_selected_weighted_intersect1=0.;   float y_TDC_selected_weighted_intersect1=0.;
    float x_GoodLG_weighted_intersect1=0.;         float y_GoodLG_weighted_intersect1=0.;
    float x_TDC_SFT_intersect1=0.;                 float y_TDC_SFT_intersect1=0.;
    float x_GoodLG_SFT_intersect1=0.;              float y_GoodLG_SFT_intersect1=0.;
    float x_TDC_Gap_Fibers_intersect1=0.;          float y_TDC_Gap_Fibers_intersect1=0.;
    float x_TDC_Gap_Fibers_SFT_intersect1=0.;      float y_TDC_Gap_Fibers_SFT_intersect1=0.;
    float x_TARGET_intersect=0;                    float y_TARGET_intersect=0;
   
    float dist1_TDC_selected_weighted[2];
    float dist1_TDC_SFT[2];
    float dist1_GoodLG_weighted[2];
    float dist1_GoodLG_SFT[2];
    float dist1_TARGET_intersect[2];
    
  
    for(int i=0; i<2; i++){
      dist1_TDC_selected_weighted[i] = distance(x_int_TDC_selected_weighted[i], y_int_TDC_selected_weighted[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TDC_SFT[i] = distance(x_int_TDC_SFT[i], y_int_TDC_SFT[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_GoodLG_weighted[i] = distance(x_int_GoodLG_weighted[i], y_int_GoodLG_weighted[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_GoodLG_SFT[i] = distance(x_int_GoodLG_SFT[i], y_int_GoodLG_SFT[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TARGET_intersect[i] = distance(x_int_TARGET[i], y_int_TARGET[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
    }    
    
  

    if(dist1_TDC_selected_weighted[0] < dist1_TDC_selected_weighted[1]){
      x_TDC_selected_weighted_intersect1 = x_int_TDC_selected_weighted[0];
      y_TDC_selected_weighted_intersect1 = y_int_TDC_selected_weighted[0];
    }
    else if(dist1_TDC_selected_weighted[1] < dist1_TDC_selected_weighted[0]){
      x_TDC_selected_weighted_intersect1 = x_int_TDC_selected_weighted[1];
      y_TDC_selected_weighted_intersect1 = y_int_TDC_selected_weighted[1];
    }
    else cout << "ERROR !" << endl;


    if(dist1_GoodLG_weighted[0] < dist1_GoodLG_weighted[1]){
      x_GoodLG_weighted_intersect1 = x_int_GoodLG_weighted[0];
      y_GoodLG_weighted_intersect1 = y_int_GoodLG_weighted[0];
    }
    else if(dist1_GoodLG_weighted[1] < dist1_GoodLG_weighted[0]){
      x_GoodLG_weighted_intersect1 = x_int_GoodLG_weighted[1];
      y_GoodLG_weighted_intersect1 = y_int_GoodLG_weighted[1];
    }
    else cout << "ERROR !" << endl;


    if(dist1_TDC_SFT[0] < dist1_TDC_SFT[1]){
      x_TDC_SFT_intersect1 = x_int_TDC_SFT[0];
      y_TDC_SFT_intersect1 = y_int_TDC_SFT[0];
    }
    else if(dist1_TDC_SFT[1] < dist1_TDC_SFT[0]){
      x_TDC_SFT_intersect1 = x_int_TDC_SFT[1];
      y_TDC_SFT_intersect1 = y_int_TDC_SFT[1];
    }
    else cout << "ERROR !" << endl;


    if(dist1_GoodLG_SFT[0] < dist1_GoodLG_SFT[1]){
      x_GoodLG_SFT_intersect1 = x_int_GoodLG_SFT[0];
      y_GoodLG_SFT_intersect1 = y_int_GoodLG_SFT[0];
    }
    else if(dist1_GoodLG_SFT[1] < dist1_GoodLG_SFT[0]){
      x_GoodLG_SFT_intersect1 = x_int_GoodLG_SFT[1];
      y_GoodLG_SFT_intersect1 = y_int_GoodLG_SFT[1];
    }
    else cout << "ERROR !" << endl;
    //    cout << "DEBUG: " << x_TDC_Gap_Fibers_intersect1 << endl;
    //    cout << "DEBUG2 : " << dist1_TDC_Gap_Fibers[0] << "  " << dist1_TDC_Gap_Fibers[1] << endl;

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
      x_TARGET_intersect = x_int_TARGET[1];
      y_TARGET_intersect = y_int_TARGET[1];
    }
    else if(dist1_TARGET_intersect[1] < dist1_TARGET_intersect[0]){
      x_TARGET_intersect = x_int_TARGET[0];
      y_TARGET_intersect = y_int_TARGET[0];
    }
    else cout << "ERROR !" << endl;

   


  h_int_TDC_selected_weighted->Fill(x_TDC_selected_weighted_intersect1, y_TDC_selected_weighted_intersect1);
  h_int_GoodLG_weighted->Fill(x_GoodLG_weighted_intersect1, y_GoodLG_weighted_intersect1);
  h_int_TDC_SFT->Fill(x_TDC_SFT_intersect1, y_TDC_SFT_intersect1);
  h_int_GoodLG_SFT->Fill(x_GoodLG_SFT_intersect1, y_GoodLG_SFT_intersect1);
    
  if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6) && Rotate==1){
    h_int_TDC_Gap_Fibers_rotate->Fill(-y_TDC_Gap_Fibers_intersect1, x_TDC_Gap_Fibers_intersect1);
    h_int_TDC_Gap_Fibers_SFT_rotate->Fill(-y_TDC_Gap_Fibers_SFT_intersect1, x_TDC_Gap_Fibers_SFT_intersect1);   
  }
  else{
      h_int_TDC_Gap_Fibers->Fill(x_TDC_Gap_Fibers_intersect1, y_TDC_Gap_Fibers_intersect1);
      h_int_TDC_Gap_Fibers_SFT->Fill(x_TDC_Gap_Fibers_SFT_intersect1, y_TDC_Gap_Fibers_SFT_intersect1);   
  }
    h_int_TDC_TARGET->Fill(x_TARGET_intersect, y_TARGET_intersect);

    TLine *CenterLine_TDC = new TLine(0., 0., x_TDC_Gap_Fibers_intersect1, y_TDC_Gap_Fibers_intersect1);
    TLine *CenterLine_GoodLG = new TLine(0., 0., x_GoodLG_weighted_intersect1, y_GoodLG_weighted_intersect1);


    CenterLine_TDC->SetLineWidth(2);    CenterLine_GoodLG->SetLineWidth(2);   
    CenterLine_TDC->SetLineColor(3);    CenterLine_GoodLG->SetLineColor(3);
    HorizontalAxis->SetLineWidth(2);    HorizontalAxis->SetLineColor(1);

    /// Angle Calculation
    ///////////////////////////////////////////////////////////////////////////////

    float a_final_TDC = 0.;           float a_final_GoodLG = 0.;             float a_final_Gap_Fibers = 0.;       float a_final_guide = 0.;
    float alpha_TDC = 0.;             float alpha_GoodLG = 0.;               float alpha_Gap_Fibers = 0.;         float alpha_guide = 0.;
    float tanalpha_TDC = 0.;          float tanalpha_GoodLG = 0.;            float tanalpha_Gap_Fibers = 0.;      float tanalpha_guide = 0.;
    float angle_final_TDC = 0.;       float angle_final_GoodLG = 0.;         float angle_final_Gap_Fibers = 0.;   float angle_final_guide = 0.;


    a_final_TDC = y_TDC_selected_weighted_intersect1 / x_TDC_selected_weighted_intersect1;
    a_final_GoodLG = y_GoodLG_weighted_intersect1 / x_GoodLG_weighted_intersect1;
    a_final_Gap_Fibers = y_TDC_Gap_Fibers_intersect1 / x_TDC_Gap_Fibers_intersect1;
    a_final_guide = (y_TDC_Gap_Fibers_intersect1 - y_TARGET_intersect) / (x_TDC_Gap_Fibers_intersect1 - x_TARGET_intersect);

    tanalpha_TDC = a_final_TDC;       tanalpha_GoodLG = a_final_GoodLG;       tanalpha_Gap_Fibers = a_final_Gap_Fibers;
    alpha_TDC = atan(tanalpha_TDC);   alpha_GoodLG = atan(tanalpha_GoodLG);   alpha_Gap_Fibers = atan(tanalpha_Gap_Fibers);

    tanalpha_guide = a_final_guide;
    alpha_guide = atan(tanalpha_guide);


    //Determination of Delta Phi
    double Delta_phi = 999.99;  double Delta_phi_deg = 999.99;
    Delta_phi = sin(alpha_guide)*cos(alpha_guide)*(ParError/a_final_guide);
    Delta_phi_deg = (180/PI)*Delta_phi;
    //cout << "DT:  " << Delta_phi_deg << endl; 


    int Axis_Vector_Length = 10;
    TArrow *x_guide;
    TArrow *y_guide;
    TArrow *x_guide_rotate;
    TArrow *y_guide_rotate;


    x_guide = new TArrow(x_TARGET_intersect, y_TARGET_intersect, x_TARGET_intersect + Axis_Vector_Length, y_TARGET_intersect, 0.005, "|>");
    y_guide = new TArrow(x_TARGET_intersect, y_TARGET_intersect, x_TARGET_intersect, y_TARGET_intersect + Axis_Vector_Length, 0.005, "|>");

    x_guide->SetLineWidth(2);     y_guide->SetLineWidth(2);
    x_guide->SetLineColor(4);     y_guide->SetLineColor(4);
    
    x_guide_rotate = new TArrow(-y_TARGET_intersect, x_TARGET_intersect, -y_TARGET_intersect + Axis_Vector_Length, x_TARGET_intersect, 0.005, "|>");
    y_guide_rotate = new TArrow(-y_TARGET_intersect, x_TARGET_intersect, -y_TARGET_intersect, x_TARGET_intersect + Axis_Vector_Length, 0.005, "|>");

    x_guide_rotate->SetLineWidth(2);     y_guide_rotate->SetLineWidth(2);
    x_guide_rotate->SetLineColor(4);     y_guide_rotate->SetLineColor(4);   


    //### Option 2
    if(x_TDC_selected_weighted_intersect1 < 0) angle_final_TDC = 180. + alpha_TDC * (180./PI);

    if(x_TDC_selected_weighted_intersect1 >= 0){
      if(y_TDC_selected_weighted_intersect1 >= 0) angle_final_TDC = alpha_TDC * (180./PI);
      else angle_final_TDC = alpha_TDC * (180./PI) + 360.;
    }


    if(x_GoodLG_weighted_intersect1 < 0) angle_final_GoodLG = 180. + alpha_GoodLG * (180./PI);

    if(x_GoodLG_weighted_intersect1 >= 0){
      if(y_GoodLG_weighted_intersect1 >= 0) angle_final_GoodLG = alpha_GoodLG * (180./PI);
      else angle_final_GoodLG = alpha_GoodLG * (180./PI) + 360.;
    }


    if(x_TDC_Gap_Fibers_intersect1 < 0)  angle_final_Gap_Fibers = 180. + alpha_Gap_Fibers * (180./PI);

    if(x_TDC_Gap_Fibers_intersect1 >= 0){
      if(y_TDC_Gap_Fibers_intersect1 >= 0)  angle_final_Gap_Fibers = alpha_Gap_Fibers * (180./PI);
      else angle_final_Gap_Fibers = alpha_Gap_Fibers * (180./PI) + 360.;
    }
      //cout << "DEBUG :" << angle_final_guide << endl;
       //cout << "DEBUG1 : " << x_TDC_Gap_Fibers_intersect1 << "  " << x_TARGET_intersect << endl;


    if((x_TDC_Gap_Fibers_intersect1 - x_TARGET_intersect) < 0) angle_final_guide = 180. + alpha_guide * (180./PI);

    if((x_TDC_Gap_Fibers_intersect1 - x_TARGET_intersect) >= 0){
      if((y_TDC_Gap_Fibers_intersect1 - y_TARGET_intersect) >= 0)  angle_final_guide = alpha_guide * (180./PI);
      else angle_final_guide = alpha_guide * (180./PI) + 360.;
    }

    char Angle_TDC_string[30];      char Angle_GoodLD_string[30];      char Angle_Gap_Fibers_string[30];      char Angle_guide_string[30];

    
    if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6) && Rotate==1)
      angle_final_guide += 90.;
    if(angle_final_guide > 360.)
    angle_final_guide -= 360.;
    

    sprintf(Angle_TDC_string,"#phi = %3.2f#circ", angle_final_TDC);
    sprintf(Angle_GoodLD_string,"#phi = %3.2f#circ", angle_final_GoodLG);
    sprintf(Angle_Gap_Fibers_string,"#phi = %3.2f#circ  (OLD)", angle_final_Gap_Fibers);
    //sprintf(Angle_guide_string,"#phi = %3.2f#circ", angle_final_guide);
    sprintf(Angle_guide_string,"#phi = %3.2f#circ#pm%3.2f", angle_final_guide, Delta_phi_deg);



    TLatex *tex_Angle_TDC;
    TLatex *tex_Angle_GoodLG;
    TLatex *tex_Angle_Gap_Fibers;
    TLatex *tex_Angle_guide;

    tex_Angle_TDC = new TLatex(-45.,43.,Angle_TDC_string);
    tex_Angle_TDC->SetTextSize(0.05);
    tex_Angle_TDC->SetLineWidth(2);

    tex_Angle_GoodLG = new TLatex(-45.,43.,Angle_GoodLD_string);
    tex_Angle_GoodLG->SetTextSize(0.05);
    tex_Angle_GoodLG->SetLineWidth(2);

    tex_Angle_Gap_Fibers = new TLatex(-45.,43.,Angle_Gap_Fibers_string);
    tex_Angle_Gap_Fibers->SetTextSize(0.05);
    tex_Angle_Gap_Fibers->SetLineWidth(2);

    tex_Angle_guide = new TLatex(-45.,43.,Angle_guide_string);
    //tex_Angle_guide->SetTextSize(0.05);
    tex_Angle_guide->SetTextSize(0.045);
    tex_Angle_guide->SetLineWidth(2);

  
    cout << "" << endl;
    cout << "" << endl;
    //cout << "Final Angle (deg.) = " << angle_final_guide << " +- " << Delta_phi_deg << "  ( GUIDE )" << endl;
    cout << "" << endl;
    cout << "" << endl;
    
    if((gap_to_fit_rotate==6 || gap_to_fit_rotate==12) && Rotate==1){
      cout << "Intersect with SFT Layer 1:  " << "x = " << -y_TDC_Gap_Fibers_SFT_intersect1 << "   " << "y = " << x_TDC_Gap_Fibers_SFT_intersect1 << endl;
      cout << "Intersect with TOF1:  " << "x = " << -y_TDC_Gap_Fibers_intersect1 << "   " << "y = " << x_TDC_Gap_Fibers_intersect1 << endl;
    }
    else{
      cout << "Intersect with SFT Layer 1:  " << "x = " << x_TDC_Gap_Fibers_SFT_intersect1 << "   " << "y = " << y_TDC_Gap_Fibers_SFT_intersect1 << endl;
      cout << "Intersect with TOF1:  " << "x = " << x_TDC_Gap_Fibers_intersect1 << "   " << "y = " << y_TDC_Gap_Fibers_intersect1 << endl;  
  }
    cout << "" << endl;

    ///////////////////////////////////////////////////////////////////////////////


    TLine *best_fit_rotate = new TLine(50*a_fit_TDC_Gap_Fibers - b_fit_TDC_Gap_Fibers,-50,-50*a_fit_TDC_Gap_Fibers - b_fit_TDC_Gap_Fibers,50); //ROTATE_CHANGE
    best_fit_rotate->SetLineWidth(2); //ROTATE_CHANGE
    best_fit_rotate->SetLineColor(kRed); // ROTATE_CHANGE
    h_Centroid->Fill(X_BAR,Y_BAR);

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
    
    if((gap_to_fit_rotate==6 || gap_to_fit_rotate==12) && Rotate==1){
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



    

    double target_slope = 0;

    if((gap_to_fit_rotate==6 || gap_to_fit_rotate==12) && Rotate==1){
      target_slope = -1.0/a_fit_TDC_Gap_Fibers;
    }  
    else{
      target_slope = a_fit_TDC_Gap_Fibers;
    }

    double c2_centroid_slope = 0.0;
    double c2_centroid_int = 0.0;
    int selected_C2L_index = 0;
    int selected_C2R_index = 0;


    double x1_rotated = 0;
    double y1_rotated = 0;
    double x2_rotated = 0;
    double y2_rotated = 0;      



    double x1_rotated_selected = 0;
    double y1_rotated_selected = 0;
    double x2_rotated_selected = 0;
    double y2_rotated_selected = 0;    
       



    double min_slope_difference = 10000000.0;



    if(gap_to_fit_left!=0 && gap_to_fit_right!=0){


      for(int i = 0; i < max_length; i++){
        if(C2Y_L_centroids[i] < -5000)
          break;

        for(int j = 0; j < max_length; j++){
          if(C2Y_R_centroids[j] < -5000)
            break;



          if(gap_to_fit_left==12 && gap_to_fit_right==6){
            c2_centroid_slope = points_to_slope(C2Y_L_centroids[i], 629.4, C2Y_R_centroids[j], -629.4);
            if(fabs(c2_centroid_slope - target_slope) < min_slope_difference){
              selected_C2L_index = i;
              selected_C2R_index = j;
              min_slope_difference = fabs(c2_centroid_slope - target_slope);

              x1_rotated_selected = C2Y_L_centroids[i];
              x2_rotated_selected = -C2Y_R_centroids[j]; //gap 6
              y1_rotated_selected = 629.4;
              y2_rotated_selected = -629.4;
            }
          }
          else if(gap_to_fit_left==12 && gap_to_fit_right==5){
            x1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,0,'x'); 
            y1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,0,'y'); 
            x2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-150,'x');
            y2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-150,'y'); 

            c2_centroid_slope = points_to_slope(x1_rotated, y1_rotated, x2_rotated, y2_rotated);

            if(fabs(c2_centroid_slope - target_slope) < min_slope_difference){
              selected_C2L_index = i;
              selected_C2R_index = j;
              min_slope_difference = fabs(c2_centroid_slope - target_slope);

              x1_rotated_selected = x1_rotated;
              x2_rotated_selected = x2_rotated;
              y1_rotated_selected = y1_rotated;
              y2_rotated_selected = y2_rotated;  
            }
          }
          else if(gap_to_fit_left==12 && gap_to_fit_right==4){
            x1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,0,'x'); 
            y1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,0,'y'); 
            x2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-120,'x');
            y2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-120,'y'); 

            c2_centroid_slope = points_to_slope(x1_rotated, y1_rotated, x2_rotated, y2_rotated);

            if(fabs(c2_centroid_slope - target_slope) < min_slope_difference){
              selected_C2L_index = i;
              selected_C2R_index = j;
              min_slope_difference = fabs(c2_centroid_slope - target_slope);

              x1_rotated_selected = x1_rotated;
              x2_rotated_selected = x2_rotated;
              y1_rotated_selected = y1_rotated;
              y2_rotated_selected = y2_rotated;  
            }
          }
          else if(gap_to_fit_left==11 && gap_to_fit_right==6){
            x1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,30,'x'); 
            y1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,30,'y'); 
            x2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,180,'x');
            y2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,180,'y'); 

            c2_centroid_slope = points_to_slope(x1_rotated, y1_rotated, x2_rotated, y2_rotated);

            if(fabs(c2_centroid_slope - target_slope) < min_slope_difference){
              selected_C2L_index = i;
              selected_C2R_index = j;
              min_slope_difference = fabs(c2_centroid_slope - target_slope);

              x1_rotated_selected = x1_rotated;
              x2_rotated_selected = x2_rotated;
              y1_rotated_selected = y1_rotated;
              y2_rotated_selected = y2_rotated;  
            }
          }
          else if(gap_to_fit_left==11 && gap_to_fit_right==5){
            x1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,30,'x'); 
            y1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,30,'y'); 
            x2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-150,'x');
            y2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-150,'y'); 

            c2_centroid_slope = points_to_slope(x1_rotated, y1_rotated, x2_rotated, y2_rotated);

            if(fabs(c2_centroid_slope - target_slope) < min_slope_difference){
              selected_C2L_index = i;
              selected_C2R_index = j;
              min_slope_difference = fabs(c2_centroid_slope - target_slope);

              x1_rotated_selected = x1_rotated;
              x2_rotated_selected = x2_rotated;
              y1_rotated_selected = y1_rotated;
              y2_rotated_selected = y2_rotated;  
            }
          }
          else if(gap_to_fit_left==11 && gap_to_fit_right==4){
            x1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,30,'x'); 
            y1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,30,'y'); 
            x2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-120,'x');
            y2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-120,'y'); 

            c2_centroid_slope = points_to_slope(x1_rotated, y1_rotated, x2_rotated, y2_rotated);

            if(fabs(c2_centroid_slope - target_slope) < min_slope_difference){
              selected_C2L_index = i;
              selected_C2R_index = j;
              min_slope_difference = fabs(c2_centroid_slope - target_slope);

              x1_rotated_selected = x1_rotated;
              x2_rotated_selected = x2_rotated;
              y1_rotated_selected = y1_rotated;
              y2_rotated_selected = y2_rotated;  
            }
          }
          else if(gap_to_fit_left==10 && gap_to_fit_right==6){
            x1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,60,'x'); 
            y1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,60,'y'); 
            x2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,180,'x');
            y2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,180,'y'); 

            c2_centroid_slope = points_to_slope(x1_rotated, y1_rotated, x2_rotated, y2_rotated);

            if(fabs(c2_centroid_slope - target_slope) < min_slope_difference){
              selected_C2L_index = i;
              selected_C2R_index = j;
              min_slope_difference = fabs(c2_centroid_slope - target_slope);

              x1_rotated_selected = x1_rotated;
              x2_rotated_selected = x2_rotated;
              y1_rotated_selected = y1_rotated;
              y2_rotated_selected = y2_rotated;  
            }
          }
          else if(gap_to_fit_left==10 && gap_to_fit_right==5){
            x1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,60,'x'); 
            y1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,60,'y'); 
            x2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-150,'x');
            y2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-150,'y'); 

            c2_centroid_slope = points_to_slope(x1_rotated, y1_rotated, x2_rotated, y2_rotated);

            if(fabs(c2_centroid_slope - target_slope) < min_slope_difference){
              selected_C2L_index = i;
              selected_C2R_index = j;
              min_slope_difference = fabs(c2_centroid_slope - target_slope);

              x1_rotated_selected = x1_rotated;
              x2_rotated_selected = x2_rotated;
              y1_rotated_selected = y1_rotated;
              y2_rotated_selected = y2_rotated;  
            }
          }
          else if(gap_to_fit_left==10 && gap_to_fit_right==4){
            x1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,60,'x'); 
            y1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,60,'y'); 
            x2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-120,'x');
            y2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-120,'y'); 

            c2_centroid_slope = points_to_slope(x1_rotated, y1_rotated, x2_rotated, y2_rotated);

            if(fabs(c2_centroid_slope - target_slope) < min_slope_difference){
              selected_C2L_index = i;
              selected_C2R_index = j;
              min_slope_difference = fabs(c2_centroid_slope - target_slope);

              x1_rotated_selected = x1_rotated;
              x2_rotated_selected = x2_rotated;
              y1_rotated_selected = y1_rotated;
              y2_rotated_selected = y2_rotated;  
            }
          }
          else if(gap_to_fit_left==8 && gap_to_fit_right==2){
            x1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,120,'x'); 
            y1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,120,'y'); 
            x2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-60,'x');
            y2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-60,'y'); 

            c2_centroid_slope = points_to_slope(x1_rotated, y1_rotated, x2_rotated, y2_rotated);

            if(fabs(c2_centroid_slope - target_slope) < min_slope_difference){
              selected_C2L_index = i;
              selected_C2R_index = j;
              min_slope_difference = fabs(c2_centroid_slope - target_slope);

              x1_rotated_selected = x1_rotated;
              x2_rotated_selected = x2_rotated;
              y1_rotated_selected = y1_rotated;
              y2_rotated_selected = y2_rotated;  
            }
          }
          else if(gap_to_fit_left==8 && gap_to_fit_right==1){
            x1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,120,'x'); 
            y1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,120,'y'); 
            x2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-30,'x');
            y2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-30,'y'); 

            c2_centroid_slope = points_to_slope(x1_rotated, y1_rotated, x2_rotated, y2_rotated);

            if(fabs(c2_centroid_slope - target_slope) < min_slope_difference){
              selected_C2L_index = i;
              selected_C2R_index = j;
              min_slope_difference = fabs(c2_centroid_slope - target_slope);

              x1_rotated_selected = x1_rotated;
              x2_rotated_selected = x2_rotated;
              y1_rotated_selected = y1_rotated;
              y2_rotated_selected = y2_rotated;  
            }
          }
          else if(gap_to_fit_left==7 && gap_to_fit_right==2){
            x1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,150,'x'); 
            y1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,150,'y'); 
            x2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-60,'x');
            y2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-60,'y'); 

            c2_centroid_slope = points_to_slope(x1_rotated, y1_rotated, x2_rotated, y2_rotated);

            if(fabs(c2_centroid_slope - target_slope) < min_slope_difference){
              selected_C2L_index = i;
              selected_C2R_index = j;
              min_slope_difference = fabs(c2_centroid_slope - target_slope);

              x1_rotated_selected = x1_rotated;
              x2_rotated_selected = x2_rotated;
              y1_rotated_selected = y1_rotated;
              y2_rotated_selected = y2_rotated;  
            }
          }
          else if(gap_to_fit_left==7 && gap_to_fit_right==1){
            x1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,150,'x'); 
            y1_rotated = rotate_by_angle(C2Y_L_centroids[i],629.4,150,'y'); 
            x2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-30,'x');
            y2_rotated = rotate_by_angle(C2Y_R_centroids[j],629.4,-30,'y'); 

            c2_centroid_slope = points_to_slope(x1_rotated, y1_rotated, x2_rotated, y2_rotated);

            if(fabs(c2_centroid_slope - target_slope) < min_slope_difference){
              selected_C2L_index = i;
              selected_C2R_index = j;
              min_slope_difference = fabs(c2_centroid_slope - target_slope);

              x1_rotated_selected = x1_rotated;
              x2_rotated_selected = x2_rotated;
              y1_rotated_selected = y1_rotated;
              y2_rotated_selected = y2_rotated;  
            }
          }
        }
      }
    }


    TLine *C2_line;

    if(min_slope_difference != 10000000.0){
      c2_centroid_slope = points_to_slope(x1_rotated_selected, y1_rotated_selected, x2_rotated_selected, y2_rotated_selected);
      c2_centroid_int = points_to_y_int(x1_rotated_selected, y1_rotated_selected, x2_rotated_selected, y2_rotated_selected);
      C2_line = new TLine(-50, c2_centroid_slope*(-50.0) + c2_centroid_int, 50, c2_centroid_slope*(50.0)+ c2_centroid_int);

        C2_line->SetLineWidth(2);
        C2_line->SetLineColor(kGreen);
    }
    else{
      C2_line = new TLine(0.0, 0.0, 0.0, 50.0);
      C2_line->SetLineWidth(2);
      C2_line->SetLineColor(kGreen);     
    }    

    
  
  if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6) && Rotate==1){ //ROTATE_CHANGE
      
      c2->cd(6)->Update();
      c2->cd(6)->Range(-50, -50, 50, 50); 
      
      h_TDC_Gap_Fibers_copy->Draw();
      h_TDC_Gap_Fibers_kaon->Draw("same");
      best_fit_rotate->Draw("same");

      A1->Draw();
      A2->Draw();
      h_Circle->Draw("same");


      C2_line->Draw("same");
   
      Gap1l->Draw();
      Gap2l->Draw("same");
      Gap3l->Draw("same");
      Gap4l->Draw("same");
      Gap5l->Draw("same");
      Gap6l->Draw("same");
      Gap7l->Draw("same");
      Gap8l->Draw("same");
      Gap9l->Draw("same");
      Gap10l->Draw("same");
      Gap11l->Draw("same");
      Gap12l->Draw("same");
      ell->Draw("same");
      ell_Target->Draw("same");

      ell_L1->Draw("same");
      h_Target_Center->Draw("same");
      h_TOF1_rotate->Draw("same");
      h_TOF1_closest->Draw("same");
      h_int_TDC_Gap_Fibers_rotate->Draw("same");
      h_int_TDC_Gap_Fibers_SFT_rotate->Draw("same"); 

      x_guide_rotate->Draw();
      y_guide_rotate->Draw();


      tex_Angle_guide->Draw("same");

      h_int_TDC_TARGET->Draw("same");

      
      h_Centroid->Draw("same"); //Centroid for k stop
      x_sft->Draw("same");
      y_sft->Draw("same");
      x_tof1->Draw("same");
      y_tof1->Draw("same");      
  }
  else{
    h_TDC_Gap_Fibers->Draw();
    h_TDC_Gap_Fibers_kaon->Draw("same");

    A1->Draw();
    A2->Draw();
    h_Circle->Draw("same");


    C2_line->Draw("same");
   
 
    Gap1l->Draw();
    Gap2l->Draw("same");
    Gap3l->Draw("same");
    Gap4l->Draw("same");
    Gap5l->Draw("same");
    Gap6l->Draw("same");
    Gap7l->Draw("same");
    Gap8l->Draw("same");
    Gap9l->Draw("same");
    Gap10l->Draw("same");
    Gap11l->Draw("same");
    Gap12l->Draw("same");
    ell->Draw("same");
    ell_Target->Draw("same");ell_Target->Draw("same");

    ell_L1->Draw("same");
    h_Target_Center->Draw("same");
    h_TOF1->Draw("same");
    h_TOF1_closest->Draw("same");
    h_int_TDC_Gap_Fibers->Draw("same");
    h_int_TDC_Gap_Fibers_SFT->Draw("same");

    x_guide->Draw();
    y_guide->Draw();

    tex_Angle_guide->Draw("same");

    h_int_TDC_TARGET->Draw("same");


    
    h_Centroid->Draw("same"); 
    x_sft->Draw("same");
    y_sft->Draw("same");
    x_tof1->Draw("same");
    y_tof1->Draw("same");     
  }


  ///////////////////////////////////////////////////////////
   //C2 Target drawing //


  double C2_gap_11_edge_1_x = rotate_by_angle(-80.0, 629.4, 30, 'x');
  double C2_gap_11_edge_2_x = rotate_by_angle(80.0, 629.4, 30, 'x'); 
  double C2_gap_11_edge_1_y = rotate_by_angle(-80.0, 629.4, 30, 'y');
  double C2_gap_11_edge_2_y = rotate_by_angle(80.0, 629.4, 30, 'y');

  double C2_gap_5_edge_1_x = rotate_by_angle(-80.0, 629.4, -150, 'x'); 
  double C2_gap_5_edge_2_x = rotate_by_angle(80.0, 629.4, -150, 'x');
  double C2_gap_5_edge_1_y = rotate_by_angle(-80.0, 629.4, -150, 'y'); 
  double C2_gap_5_edge_2_y = rotate_by_angle(80.0, 629.4, -150, 'y'); 

  double C2_gap_12_edge_1_x = rotate_by_angle(-80.0, 629.4, 0, 'x');
  double C2_gap_12_edge_2_x = rotate_by_angle(80.0, 629.4, 0, 'x');
  double C2_gap_12_edge_1_y = rotate_by_angle(-80.0, 629.4, 0, 'y');
  double C2_gap_12_edge_2_y = rotate_by_angle(80.0, 629.4, 0, 'y');

  double C2_gap_6_edge_1_x = rotate_by_angle(-80.0, 629.4, 180, 'x');
  double C2_gap_6_edge_2_x = rotate_by_angle(80.0, 629.4, 180, 'x');
  double C2_gap_6_edge_1_y = rotate_by_angle(-80.0, 629.4, 180, 'y');
  double C2_gap_6_edge_2_y = rotate_by_angle(80.0, 629.4, 180, 'y');

  double C2_gap_1_edge_1_x = rotate_by_angle(-80.0, 629.4, -30, 'x'); 
  double C2_gap_1_edge_2_x = rotate_by_angle(80.0, 629.4, -30, 'x'); 
  double C2_gap_1_edge_1_y = rotate_by_angle(-80.0, 629.4, -30, 'y'); 
  double C2_gap_1_edge_2_y = rotate_by_angle(80.0, 629.4, -30, 'y'); 

  double C2_gap_7_edge_1_x = rotate_by_angle(-80.0, 629.4, 150, 'x'); 
  double C2_gap_7_edge_2_x = rotate_by_angle(80.0, 629.4, 150, 'x'); 
  double C2_gap_7_edge_1_y = rotate_by_angle(-80.0, 629.4, 150, 'y'); 
  double C2_gap_7_edge_2_y = rotate_by_angle(80.0, 629.4, 150, 'y'); 

  double C2_gap_2_edge_1_x = rotate_by_angle(-80.0, 629.4, -60, 'x'); 
  double C2_gap_2_edge_2_x = rotate_by_angle(80.0, 629.4, -60, 'x'); 
  double C2_gap_2_edge_1_y = rotate_by_angle(-80.0, 629.4, -60, 'y');
  double C2_gap_2_edge_2_y = rotate_by_angle(80.0, 629.4, -60, 'y');

  double C2_gap_4_edge_1_x = rotate_by_angle(-80.0, 629.4, -120, 'x'); 
  double C2_gap_4_edge_2_x = rotate_by_angle(80.0, 629.4, -120, 'x'); 
  double C2_gap_4_edge_1_y = rotate_by_angle(-80.0, 629.4, -120, 'y');
  double C2_gap_4_edge_2_y = rotate_by_angle(80.0, 629.4, -120, 'y');  

  double C2_gap_8_edge_1_x = rotate_by_angle(-80.0, 629.4, 120, 'x'); 
  double C2_gap_8_edge_2_x = rotate_by_angle(80.0, 629.4, 120, 'x');
  double C2_gap_8_edge_1_y = rotate_by_angle(-80.0, 629.4, 120, 'y');  
  double C2_gap_8_edge_2_y = rotate_by_angle(80.0, 629.4, 120, 'y');

  double C2_gap_10_edge_1_x = rotate_by_angle(-80.0, 629.4, 60, 'x');
  double C2_gap_10_edge_2_x = rotate_by_angle(80.0, 629.4, 60, 'x'); 
  double C2_gap_10_edge_1_y = rotate_by_angle(-80.0, 629.4, 60, 'y');
  double C2_gap_10_edge_2_y = rotate_by_angle(80.0, 629.4, 60, 'y'); 

  TLine *C2_line_display;


  if(min_slope_difference != 10000000.0){
    C2_line_display = new TLine(-700, c2_centroid_slope*(-700.0) + c2_centroid_int, 700, c2_centroid_slope*(700.0)+ c2_centroid_int);
   }
  else{
    C2_line_display = new TLine(0, 0, 0, 0);
   }
  


  C2_line_display->SetLineWidth(0.5);
  C2_line_display->SetLineColor(kGreen);



  TLine *C2_gap_11 = new TLine(C2_gap_11_edge_1_x,C2_gap_11_edge_1_y,C2_gap_11_edge_2_x,C2_gap_11_edge_2_y);
  C2_gap_11->SetLineWidth(4);
  C2_gap_11->SetLineColor(kBlack);

  TLine *C2_gap_5 = new TLine(C2_gap_5_edge_1_x,C2_gap_5_edge_1_y,C2_gap_5_edge_2_x,C2_gap_5_edge_2_y);
  C2_gap_5->SetLineWidth(4);
  C2_gap_5->SetLineColor(kBlack);

  TLine *C2_gap_12 = new TLine(C2_gap_12_edge_1_x,C2_gap_12_edge_1_y,C2_gap_12_edge_2_x,C2_gap_12_edge_2_y);
  C2_gap_12->SetLineWidth(4);
  C2_gap_12->SetLineColor(kBlack);

  TLine *C2_gap_6 = new TLine(C2_gap_6_edge_1_x,C2_gap_6_edge_1_y,C2_gap_6_edge_2_x,C2_gap_6_edge_2_y);
  C2_gap_6->SetLineWidth(4);
  C2_gap_6->SetLineColor(kBlack);  

  TLine *C2_gap_1 = new TLine(C2_gap_1_edge_1_x,C2_gap_1_edge_1_y,C2_gap_1_edge_2_x,C2_gap_1_edge_2_y);
  C2_gap_1->SetLineWidth(4);
  C2_gap_1->SetLineColor(kBlack);

  TLine *C2_gap_7 = new TLine(C2_gap_7_edge_1_x,C2_gap_7_edge_1_y,C2_gap_7_edge_2_x,C2_gap_7_edge_2_y);
  C2_gap_7->SetLineWidth(4);
  C2_gap_7->SetLineColor(kBlack);    

  TLine *C2_gap_2 = new TLine(C2_gap_2_edge_1_x,C2_gap_2_edge_1_y,C2_gap_2_edge_2_x,C2_gap_2_edge_2_y);
  C2_gap_2->SetLineWidth(4);
  C2_gap_2->SetLineColor(kBlack);    

  TLine *C2_gap_4 = new TLine(C2_gap_4_edge_1_x,C2_gap_4_edge_1_y,C2_gap_4_edge_2_x,C2_gap_4_edge_2_y);
  C2_gap_4->SetLineWidth(4);
  C2_gap_4->SetLineColor(kBlack);    

  TLine *C2_gap_8 = new TLine(C2_gap_8_edge_1_x,C2_gap_8_edge_1_y,C2_gap_8_edge_2_x,C2_gap_8_edge_2_y);
  C2_gap_8->SetLineWidth(4);
  C2_gap_8->SetLineColor(kBlack);    

  TLine *C2_gap_10 = new TLine(C2_gap_10_edge_1_x,C2_gap_10_edge_1_y,C2_gap_10_edge_2_x,C2_gap_10_edge_2_y);
  C2_gap_10->SetLineWidth(4);
  C2_gap_10->SetLineColor(kBlack);          


  ell_C2->SetLineColor(kBlue);



  /// C2 lines
  //TH2F *h_target_ADC = new TH2F("Histo Fit 0",h_target_ADC_title,3000,-50,50,3000,-50,50);
  char C2_z_histo_title[200];
  sprintf(C2_z_histo_title,"C2, SFT, and Target Side View  | Run Number %d  | Event %d", Run_Number, ievt);

  TH2F *C2_z_histo = new TH2F("C2 Chambers, SFT, and Target Side view", C2_z_histo_title, 500,-300,300,500,-700,700);

  TLine *target_top_line = new TLine(-100,R_TARGET,100,R_TARGET);
  TLine *target_bottom_line = new TLine(-100,-R_TARGET,100,-R_TARGET);
  TLine *target_left_line = new TLine(-100,R_TARGET,-100,-R_TARGET);
  TLine *target_right_line = new TLine(100,R_TARGET,100,-R_TARGET);

  TLine *C2_top_line = new TLine(-280,629.4,280,629.4);
  TLine *C2_bottom_line = new TLine(-280,-629.4,280,-629.4);

  TLine *SFT_top_line = new TLine(-135,R_SFT_L1,135,R_SFT_L1);
  TLine *SFT_bottom_line = new TLine(-135,-R_SFT_L1,135,-R_SFT_L1);

  ///////////////////
  //C2X_L_centroids//
  ///////////////////
  int c2x_num_lines = 0;
  int c2xl_num_lines = 0;
  int c2xr_num_lines = 0;



  for(int i = 0; i < max_length; i++){
    if(C2X_L_centroids[i] < -500)
      break;

    c2xl_num_lines++;
  }

  for(int j = 0; j < max_length; j++){
    if(C2X_R_centroids[j] < -500)
      break;

    c2xr_num_lines++;
  }
  


  TLine *C2X_line[c2xl_num_lines][c2xr_num_lines];

  cout << "C2 Chambers Z values: " << endl;
  for(int i = 0; i < c2xl_num_lines; i++){
    for(int j = 0; j < c2xr_num_lines; j++){
      if(gap_to_fit_right==1 || gap_to_fit_right==2){
        C2X_line[i][j] = new TLine(C2X_L_centroids[i],-629.4,C2X_R_centroids[j],629.4);
      }
      else{
        C2X_line[i][j] = new TLine(C2X_L_centroids[i],629.4,C2X_R_centroids[j],-629.4);
      } 

      C2X_line[i][j]->SetLineColor(kGreen);

      double x1_C2_coord = 0.0;
      double x2_C2_coord = 0.0;
      double y1_C2_coord = 0.0;
      double y2_C2_coord = 0.0;


      if(gap_to_fit_right==1 || gap_to_fit_right==2){
        x1_C2_coord = C2X_R_centroids[j];  
        x2_C2_coord = C2X_L_centroids[i];   
      }
      else{
        x1_C2_coord = C2X_L_centroids[i]; 
        x2_C2_coord = C2X_R_centroids[j];  
      }

      y1_C2_coord = 629.4; 
      y2_C2_coord = -629.4;


      C2_intersect_SFT_top = C2_intersect(points_to_slope(x1_C2_coord,y1_C2_coord,x2_C2_coord,y2_C2_coord),points_to_y_int(x1_C2_coord,y1_C2_coord,x2_C2_coord,y2_C2_coord),R_SFT_L1);
      C2_intersect_SFT_bottom = C2_intersect(points_to_slope(x1_C2_coord,y1_C2_coord,x2_C2_coord,y2_C2_coord),points_to_y_int(x1_C2_coord,y1_C2_coord,x2_C2_coord,y2_C2_coord),-R_SFT_L1);

      if(x1_C2_coord == x2_C2_coord){
      C2_intersect_SFT_top = x1_C2_coord;
      C2_intersect_SFT_bottom = x1_C2_coord;
     }
  
      if(C2_intersect_SFT_top > -135 && C2_intersect_SFT_top < 135) 
        cout << "C2 Z intersect = " << C2_intersect_SFT_top << "  --  Top "<< endl;

      if(C2_intersect_SFT_bottom > -135 && C2_intersect_SFT_bottom < 135)
        cout << "C2 Z intersect = " << C2_intersect_SFT_bottom << "  --  Bottom "<< endl;
    }
  }
    cout << endl;
    cout << "Final Angle (deg.) = " << angle_final_guide << " +- " << Delta_phi_deg << "  ( GUIDE )" << endl;
    cout << endl;

  if(c2xl_num_lines==0 || c2xr_num_lines==0){
    C2_intersect_SFT_top = 999.99;
    C2_intersect_SFT_bottom = 999.99;
  }

  C2_top_line->SetLineWidth(2);
  C2_bottom_line->SetLineWidth(2);

  C2_top_line->SetLineColor(kBlue);
  C2_bottom_line->SetLineColor(kBlue);

  SFT_top_line->SetLineColor(kRed);
  SFT_bottom_line->SetLineColor(kRed);


  // TARGET XY view
    TCanvas *C2_Target_canvas;
    C2_Target_canvas = new TCanvas("C2 with Target","C2 with Target",50,50,1000,500);
    C2_Target_canvas->Divide(2,1);
    C2_Target_canvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");


    C2_Target_canvas->cd(1);

    h_C2_target->Draw();
    ell_C2->Draw("same");
    ell_Target->Draw("same");

    C2_gap_11->Draw("same");
    C2_gap_5->Draw("same");
    C2_gap_12->Draw("same");
    C2_gap_6->Draw("same");  
    C2_gap_1->Draw("same");
    C2_gap_7->Draw("same"); 
    C2_gap_2->Draw("same"); 
    C2_gap_4->Draw("same"); 
    C2_gap_8->Draw("same"); 
    C2_gap_10->Draw("same");         
    C2_line_display->Draw("same");




    //// ADD SECOND CANVAS, showing z-view of c2 chambers and sft. draw closest line.
    //////////////////////////////////////////////////////////////////////////////////

    C2_Target_canvas->cd(2);

    C2_z_histo->Draw();

    target_top_line->Draw("same");
    target_bottom_line->Draw("same");
    target_left_line->Draw("same");
    target_right_line->Draw("same");

    C2_top_line->Draw("same");
    C2_bottom_line->Draw("same");    

    SFT_top_line->Draw("same");
    SFT_bottom_line->Draw("same");


    for(int i = 0; i < c2xl_num_lines; i++){
      for(int j = 0; j < c2xr_num_lines; j++){
        C2X_line[i][j]->Draw("same");
      }
    }    



  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //cout << "DEBUG: " <<  C2_intersect_SFT_top << "  " << C2_intersect_SFT_bottom << endl;

  SFT_angle_selector(Run_Number,ievt,angle_final_guide, gap_to_fit_left, gap_to_fit_right);
  
  if(ana_flag == 1) Ana3(Run_Number, ievt, ADC_High_SFT_corr, has_TDC_SFT_hit, SFT_channel_to_fiber, angle_final_guide, gap_to_fit_left, gap_to_fit_right, C2_intersect_SFT_top, C2_intersect_SFT_bottom);

  
  return;
}

void SFT_angle_selector(int Run_Number, int ievt, double phi, int gap_to_fit_left, int gap_to_fit_right){

  if(gap_to_fit_left==12 && gap_to_fit_right==6 ){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,12);
      SFT_Test_CR_Event(Run_Number,ievt,phi,6);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,12);
      SFT_Test_CR_Event(Run_Number,ievt,phi+180,6);
    }
  }
  if(gap_to_fit_left==12 && gap_to_fit_right==5){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,12);
      SFT_Test_CR_Event(Run_Number,ievt,phi,5);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,12);
      SFT_Test_CR_Event(Run_Number,ievt,phi+180,5);
    }    
  }
  if(gap_to_fit_left==12 && gap_to_fit_right==4 ){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,12);
      SFT_Test_CR_Event(Run_Number,ievt,phi,4);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,12);
      SFT_Test_CR_Event(Run_Number,ievt,phi+180,4);
    }      
  }
  if(gap_to_fit_left==11 && gap_to_fit_right==6 ){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,11);
      SFT_Test_CR_Event(Run_Number,ievt,phi,6);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,11);
      SFT_Test_CR_Event(Run_Number,ievt,phi+180,6);
    }      
  }
  if(gap_to_fit_left==11 && gap_to_fit_right==5 ){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,11);
      SFT_Test_CR_Event(Run_Number,ievt,phi,5);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,11);
      SFT_Test_CR_Event(Run_Number,ievt,phi+180,5);
    }        
  }
  if(gap_to_fit_left==11 && gap_to_fit_right==4 ){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,11);
      SFT_Test_CR_Event(Run_Number,ievt,phi,4);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,11);
      SFT_Test_CR_Event(Run_Number,ievt,phi+180,4);
    }      
  }
  if(gap_to_fit_left==10 && gap_to_fit_right==6 ){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,10);
      SFT_Test_CR_Event(Run_Number,ievt,phi,6);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,10);
      SFT_Test_CR_Event(Run_Number,ievt,phi+180,6);
    }       
  }
  if(gap_to_fit_left==10 && gap_to_fit_right==5 ){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,10);
      SFT_Test_CR_Event(Run_Number,ievt,phi,5);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,10);
      SFT_Test_CR_Event(Run_Number,ievt,phi+180,5);
    }      
  }
  if(gap_to_fit_left==10 && gap_to_fit_right==4 ){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,10);
      SFT_Test_CR_Event(Run_Number,ievt,phi,4);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,10);
      SFT_Test_CR_Event(Run_Number,ievt,phi+180,4);
    }      
  }
  if(gap_to_fit_left==8 && gap_to_fit_right==2 ){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,2);
      SFT_Test_CR_Event(Run_Number,ievt,phi,8);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,2);
      SFT_Test_CR_Event(Run_Number,ievt,phi+180,8);
    }      
  }
  if(gap_to_fit_left==8 && gap_to_fit_right==1 ){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,1);
      SFT_Test_CR_Event(Run_Number,ievt,phi,8);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,1);
      SFT_Test_CR_Event(Run_Number,ievt,phi+180,8);
    }   
  }
  if(gap_to_fit_left==7 && gap_to_fit_right==2 ){
    if(phi > 270){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,2);
      SFT_Test_CR_Event(Run_Number,ievt,phi,7);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,2);
      SFT_Test_CR_Event(Run_Number,ievt,phi+180,7);
    }   
  }
  if(gap_to_fit_left==7 && gap_to_fit_right==1){
    if(phi > 270){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,1);
      SFT_Test_CR_Event(Run_Number,ievt,phi,7);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,1);
      SFT_Test_CR_Event(Run_Number,ievt,phi+180,7);
    }   
  }

 

  return;
}



void SFT_Test_CR_Event(int Run_Number, int evt, double phi, int gap_to_fit){
  Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];        
  Int_t adc_low_sft[128];           Int_t ADC_Low_SFT[128];   
  Int_t tdc_le_sft[128][16];        Int_t TDC_LE_SFT[128];         
  Int_t tdc_te_sft[128][16];        Int_t TDC_TE_SFT[128];  

  Int_t HG_SFT_ADC_Thr[128] = {0};
  Int_t LG_SFT_ADC_Thr[128] = {0};

  
  for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;
  for(int i=0; i<128; i++)  LG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_LG[i]) + SFT_ADC_Thr_LG_Offset;

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
    TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
    ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-LG_SFT_ADC_Thr[j_SFT];
  }

  for(int j=0 ; j<128 ; j++){ 
    if(ADC_High_SFT[j]<0)     
      ADC_High_SFT_corr[j]=0; 
    else
      ADC_High_SFT_corr[j]=ADC_High_SFT[j]; 
  }  

  for(Int_t ii=0; ii<128; ii++){
    for (Int_t qq=0; qq<6; qq++) {
      if (tdc_le_sft[ii][qq] > SFT_TDC_min[ii] && tdc_le_sft[ii][qq] < SFT_TDC_max[ii]) has_TDC_SFT_hit[ii]++;
    }
  }


 if(ana_flag == 2) double z = SFT_print(ADC_High_SFT_corr, has_TDC_SFT_hit, SFT_channel_to_fiber, evt, phi, true, 0, false, 0, 0);
 // Ana2(ADC_High_SFT_corr, has_TDC_SFT_hit, SFT_channel_to_fiber, evt, phi, true, 0, false, 0, 0);

  //cout << endl;
  //cout << endl; 
}




