#ifndef __CINT__
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>
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
#include "TStyle.h"  
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TBranch.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TEllipse.h"
#include "TMarker.h"
#include "ANAPATH.h"
#include "Thresholds.h"
#include "CommonParameters.h"
#endif

#include "intersect.cxx"
  
void Event_Display_dev(Int_t Run_Number=5, Int_t ievt=0, Int_t enable_cout=0, Int_t display_fits = 0)
{ 

  gStyle->Clear();
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(11);
  
  char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!

  Double_t flag_size_TARGET=1.35;
  Double_t flag_size_SFT=1.3;
  Double_t flag_size_palette=1.6;


  Int_t ADC_High_corr_max=0;

  Int_t ADC_cut_TARGET2 = 850;

  char ADC_cut[100];    sprintf(ADC_cut,"(ADC >= %d)",ADC_cut_SFT);

  Int_t adc_high_target[256];       Int_t ADC_High_TARGET[256];    
  Int_t adc_low_target[256];        Int_t ADC_Low_TARGET[256];  
  Int_t tdc_le_target[256][16];     Int_t TDC_LE_TARGET[256];     
  Int_t tdc_te_target[256][16];     //Int_t TDC_TE_TARGET[256];  

  Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];        Double_t ADC_High_SFT_corr[128];    
  Int_t adc_low_sft[128];           //Int_t ADC_Low_SFT[128];   
  Int_t tdc_le_sft[128][16];        Int_t TDC_LE_SFT[128];         
  Int_t tdc_te_sft[128][16];        //Int_t TDC_TE_SFT[128];  

  Int_t ADC_tof1[24];               Int_t ADC_TOF1[24];   Int_t ADC_TOF1U[12];  Int_t ADC_TOF1D[12];
  Int_t ADC_tof2[56];               Int_t ADC_TOF2[56];
  Int_t ADC_TOF2AO[12];             Int_t ADC_TOF2AI[12];
  Int_t ADC_TOF2BO[12];             Int_t ADC_TOF2BI[12];

  Int_t TDC_tof1U[12];              Int_t TDC_TOF1U[12];
  Int_t TDC_tof1D[12];              Int_t TDC_TOF1D[12];

  Int_t TDC_tof2AO[12];             Int_t TDC_TOF2AO[12];
  Int_t TDC_tof2BO[12];             Int_t TDC_TOF2BO[12];
  Int_t TDC_tof2AI[12];             Int_t TDC_TOF2AI[12];
  Int_t TDC_tof2BI[12];             Int_t TDC_TOF2BI[12];

  Int_t MwpcADC[512];               Int_t MWPCADC[512];

  Int_t has_TDC_SFT_hit[128] = {0};

  float R_TOF1 = 44.5;
  float R_SFT_L1 = 40.0;


  char run_string[100];            char event_string[100];

  char h_ADC_title[200];

  sprintf(h_ADC_title,"(ADC #geq %d) | (%d #leq TDC #leq %d)  --  Run %d (Event %d)",ADC_cut_SFT, TDC_min_SFT,TDC_max_SFT,Run_Number,ievt);


 /////////////////////////////////////////////////////////////////////////////////////////////////////

  char h_target_ADC_title[100];
  char h_target_TDC_title[100];
  sprintf(h_target_ADC_title,"Run %d (Event %d) -- ADC #geq %d)",Run_Number, ievt, ADC_cut_TARGET);
  sprintf(h_target_TDC_title,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt, ADC_cut_TARGET,TDC_min_TARGET,TDC_max_TARGET);

  char h_target_ADC_title2[100];
  char h_target_TDC_title2[100];
  char h_target_ADC_title3[100];
  char h_target_ADC_title4[100];
  sprintf(h_target_ADC_title2,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt, ADC_cut_TARGET,TDC_min_TARGET,TDC_max_TARGET);
  sprintf(h_target_TDC_title2,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt, ADC_cut_TARGET,TDC_min_TARGET,TDC_max_TARGET);
  sprintf(h_target_ADC_title3,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt, ADC_cut_TARGET,TDC_min_TARGET,TDC_max_TARGET);
  sprintf(h_target_ADC_title4,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt, ADC_cut_TARGET,TDC_min_TARGET,TDC_max_TARGET);

  char h_target_TDC_copy_Name[200];     char h_target_TDC_copy_Title[200];
  char h_TDC_selected_Name[200];        char h_TDC_selected_Title[200];
  char h_GoodLG_Name[200];              char h_GoodLG_Title[200];
  char h_TDC_selected2_Name[200];       char h_TDC_selected2_Title[200];
  char h_GoodLG_weighted_Name[200];     char h_GoodLG_weighted_Title[200];

  sprintf(h_target_TDC_copy_Name,"Event %d (Run %d)", ievt, Run_Number);
  sprintf(h_target_TDC_copy_Title,"ADC_HG #geq %d  |  %d #leq TDC #leq %d", ADC_cut_TARGET,TDC_min_TARGET,TDC_max_TARGET);
 
  sprintf(h_TDC_selected_Name,"Event %d (Run %d)", ievt, Run_Number);
  sprintf(h_TDC_selected_Title,"ADC_HG #geq %d  |  %d #leq TDC #leq %d", ADC_cut_TARGET+300,TDC_min_TARGET,TDC_max_TARGET);

  sprintf(h_GoodLG_Name,"Event %d (Run %d)", ievt, Run_Number);
  sprintf(h_GoodLG_Title,"ADC_LG #geq %d  |  %d #leq TDC #leq %d", ADC_cut_TARGET,TDC_min_TARGET,TDC_max_TARGET);

  sprintf(h_TDC_selected2_Name,"Event %d (Run %d)", ievt, Run_Number);
  sprintf(h_TDC_selected2_Title,"ADC_HG #geq %d  |  %d #leq TDC #leq %d  (WEIGHTED)", ADC_cut_TARGET+300,TDC_min_TARGET,TDC_max_TARGET);

  sprintf(h_GoodLG_weighted_Name,"Event %d (Run %d)", ievt, Run_Number);
  sprintf(h_GoodLG_weighted_Title,"ADC_LG #geq %d  |  %d #leq TDC #leq %d  (WEIGHTED)", ADC_cut_TARGET,TDC_min_TARGET,TDC_max_TARGET);


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

  //TH2F *h_GoodLG = new TH2F(h_GoodLG_Name, h_GoodLG_Title, 500, -50, 50, 500, -50, 50);
  //TH2F *h_GoodLG_weighted = new TH2F(h_GoodLG_weighted_Name, h_GoodLG_weighted_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_GoodLG_copy = new TH2F("Good LG Copy", "Good LG Copy", 500, -50, 50, 500, -50, 50);

  TH2F *h_TOF1 = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_closest = new TH2F("h_TOF1_closest", "h_TOF1_closest", 500, -50, 50, 500, -50, 50);

  TH2F *h_int_TDC = new TH2F("h_int_TDC", "h_int_TDC", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_selected = new TH2F("h_int_TDC_selected", "h_int_TDC_selected", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_selected_weighted = new TH2F("h_int_TDC_selected_weighted", "h_int_TDC_selected_weighted", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_Gap_Fibers = new TH2F("h_int_TDC_Gap_Fibers", "h_int_TDC_Gap_Fibers", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_SFT = new TH2F("h_int_TDC_SFT", "h_int_TDC_SFT", 500, -50, 50, 500, -50, 50);

  TH2F *h_TDC_Gap_Fibers = new TH2F("h_TDC_Gap_Fibers", "h_TDC_Gap_Fibers", 500, -50, 50, 500, -50, 50);
  TEllipse *ell = new TEllipse(0, 0, R_TOF1, 0);
  TEllipse *ell_Target = new TEllipse(0, 0, 29, 0);
  TEllipse *ell_L1 = new TEllipse(0, 0, R_SFT_L1, 0);

  TLine *VerticalAxis = new TLine(0., 0., 0., 50.);


  char path_input[200];                   char file_mapping[200];
  sprintf(path_input,"%s",path_merged);          


  sprintf(file_mapping,"../Mapping");

  char Name_finput[200];
  sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);



  char par_finput[200];
  sprintf(par_finput,"%s/%s",file_mapping,source_mapping);

  char par_finput2[200];
  sprintf(par_finput2,"%s/MWPC_map.txt",file_mapping);

  char par_finput3[200];
  sprintf(par_finput3,"%s/ADC_TARGET_Thresholds.txt",file_mapping);

  char par_finput4[200];
  sprintf(par_finput4,"%s/ADC_SFT_Thresholds.txt",file_mapping);

  char par_finput5[200];
  sprintf(par_finput5,"%s/ADC_TOF1_Thresholds.txt",file_mapping);

  char par_finput6[200];
  sprintf(par_finput6,"%s/ADC_TOF2_Thresholds.txt",file_mapping);



  char footer[100];
  sprintf(footer,"Event_Display_MS.C  --  Run %d ; Event %d",Run_Number,ievt);

  cout << "   " << endl;
  cout << "   " << endl;
  cout << "*****************************************************" << endl;

  //cout << "   " << endl;
  cout << "File opened:  " << Name_finput << endl;
  cout << "SFT Mapping File:  " << par_finput << endl;
  cout << "MWPC Mapping File:  " << par_finput2 << endl;



  Int_t par_temp[2][128];
  ifstream fdat(par_finput,ios::in);
  for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
  fdat.close();

  char par_temp2[512][50];
  ifstream fdat2(par_finput2,ios::in);
  for(Int_t ii=0; ii<512; ii++) fdat2 >> par_temp2[ii];
  fdat2.close();

  Int_t par_temp_TARGET[2][256];
  ifstream fdat3(par_finput3,ios::in);
  for(Int_t ii=0; ii<256; ii++) fdat3 >> par_temp_TARGET[0][ii] >> par_temp_TARGET[1][ii];
  fdat3.close();

  Int_t par_temp_SFT[2][128];
  ifstream fdat4(par_finput4,ios::in);
  for(Int_t ii=0; ii<128; ii++) fdat4 >> par_temp_SFT[0][ii] >> par_temp_SFT[1][ii];
  fdat4.close();

  Int_t par_temp_TOF1[2][24];
  ifstream fdat5(par_finput5,ios::in);
  for(Int_t ii=0; ii<24; ii++) fdat5 >> par_temp_TOF1[0][ii] >> par_temp_TOF1[1][ii];
  fdat5.close();

  Int_t par_temp_TOF2[2][56];
  ifstream fdat6(par_finput6,ios::in);
  for(Int_t ii=0; ii<56; ii++) fdat6 >> par_temp_TOF2[0][ii] >> par_temp_TOF2[1][ii];
  fdat6.close();


  sprintf(run_string,"Run %d ; Event %d",Run_Number,ievt);
  sprintf(event_string,"Run %d ; Event %d",Run_Number,ievt);


  TChain *fChain= new TChain("Tree");   
  fChain->Add(Name_finput);   
  fChain->SetMakeClass(1);              

  fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);    fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
  fChain->SetBranchAddress("ADC_Low_TARGET",adc_low_target);      fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
  fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);        fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
  fChain->SetBranchAddress("TDC_TE_TARGET",tdc_te_target);        fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

  fChain->SetBranchAddress("ADC_TOF1",ADC_tof1);
  fChain->SetBranchAddress("ADC_TOF2",ADC_tof2);

  fChain->SetBranchAddress("TDC_TOF1U",TDC_tof1U);
  fChain->SetBranchAddress("TDC_TOF1D",TDC_tof1D);
  fChain->SetBranchAddress("TDC_TOF2AO",TDC_tof2AO);
  fChain->SetBranchAddress("TDC_TOF2AI",TDC_tof2AI);
  fChain->SetBranchAddress("TDC_TOF2BO",TDC_tof2BO);
  fChain->SetBranchAddress("TDC_TOF2BI",TDC_tof2BI);    

  fChain->SetBranchAddress("MWPCADC",MwpcADC);

  Int_t nentries = (Int_t)fChain->GetEntries();
  cout << "Total Number of Events:  " << nentries <<endl;
  //cout << "  " << endl;
  cout << "*****************************************************" << endl;
  cout << "  " << endl;
  cout << "  " << endl;

  bool Good_Event=false;
  bool Good_tof1[12] = {false};
  bool Good_tof2[12] = {false};

  for(Int_t i=ievt; i<ievt+1; i++){
    fChain->GetEntry(i);  
   
    for (int j_TARGET=0; j_TARGET<256; j_TARGET++){
      ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-ADC_cut_TARGET;
      ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-ADC_cut_TARGET;
      TDC_LE_TARGET[j_TARGET]=tdc_le_target[j_TARGET][0];
      //TDC_TE_TARGET[j_TARGET]=tdc_te_target[j_TARGET][0];  
    } 


    for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
      ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-par_temp_SFT[1][j_SFT];
      //ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-par_temp_SFT[1][j_SFT];
      TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
      //TDC_TE_SFT[j_SFT]=tdc_te_sft[j_SFT][0]; 
    }

    for (Int_t j_TOF1=0; j_TOF1<24; j_TOF1++) {
      ADC_TOF1[j_TOF1] = ADC_tof1[j_TOF1]-par_temp_TOF1[1][j_TOF1];
    }

    for(Int_t j_TOF1=0; j_TOF1<12; j_TOF1++){
      ADC_TOF1U[j_TOF1] = ADC_TOF1[j_TOF1];
      ADC_TOF1D[j_TOF1] = ADC_TOF1[j_TOF1+12];
    }


    for (Int_t j_TOF2=0; j_TOF2<56; j_TOF2++) {
      ADC_TOF2[j_TOF2] = ADC_tof2[j_TOF2]-par_temp_TOF2[1][j_TOF2];
    }

    for(Int_t j_TOF2=0; j_TOF2<12; j_TOF2++){
      ADC_TOF2AO[j_TOF2] = ADC_TOF2[j_TOF2];
      ADC_TOF2BO[j_TOF2] = ADC_TOF2[j_TOF2+12];
      ADC_TOF2AI[j_TOF2] = ADC_TOF2[j_TOF2+24];
      ADC_TOF2BI[j_TOF2] = ADC_TOF2[j_TOF2+36];
    }
    ADC_TOF2BO[6] = ADC_TOF2[55];


    for (Int_t j_MWPC=0; j_MWPC<512; j_MWPC++) {
      MWPCADC[j_MWPC] = MwpcADC[j_MWPC]-MWPC_thr;
    }

    for (Int_t j_TDCTOF=0; j_TDCTOF<12; j_TDCTOF++) {
      TDC_TOF1U[j_TDCTOF] = TDC_tof1U[j_TDCTOF];
      TDC_TOF1D[j_TDCTOF] = TDC_tof1D[j_TDCTOF];
      TDC_TOF2AO[j_TDCTOF] = TDC_tof2AO[j_TDCTOF];
      TDC_TOF2BO[j_TDCTOF] = TDC_tof2BO[j_TDCTOF];
      TDC_TOF2AI[j_TDCTOF] = TDC_tof2AI[j_TDCTOF];
      TDC_TOF2BI[j_TDCTOF] = TDC_tof2BI[j_TDCTOF];
    }  
  

    /// FILTERS
    //********* GOOD TARGET EVENTS
    bool Good_TARGET_Event = false;
    int count_TARGET_evts = 0; 
    for(int i=0; i<256; i++)
    {
      if((ADC_High_TARGET[i]>=0 && tdc_le_target[i][0]>=TDC_min_TARGET && tdc_le_target[i][0]<=TDC_max_TARGET) ||
         (ADC_High_TARGET[i]>=0 && tdc_le_target[i][1]>=TDC_min_TARGET && tdc_le_target[i][1]<=TDC_max_TARGET) ||
         (ADC_High_TARGET[i]>=0 && tdc_le_target[i][2]>=TDC_min_TARGET && tdc_le_target[i][2]<=TDC_max_TARGET) ||
         (ADC_High_TARGET[i]>=0 && tdc_le_target[i][3]>=TDC_min_TARGET && tdc_le_target[i][3]<=TDC_max_TARGET))
      { 
        count_TARGET_evts++;
      }
    }

    if(count_TARGET_evts >= 5) Good_TARGET_Event = true;

    
    //********* GOOD TOF EVENTS
    bool Good_TOF1_ADC[12]={false};   bool Good_TOF2_ADC[12]={false};
    bool Good_TOF1_TDC[12]={false};   bool Good_TOF2_TDC[12]={false};
    bool Good_TOF1[12]={false};       bool Good_TOF2[12]={false};
    bool Good_TOFs[12]={false};
    bool Good_TOF_Event = false;


    for(int i=0; i<12; i++)
    {
      if(ADC_TOF1U[i]>=0 || ADC_TOF1D[i]>=0)  Good_TOF1_ADC[i] = true;
      
      if((TDC_TOF1U[i]>=TDC_TOF1_min && TDC_TOF1U[i]<=TDC_TOF1_max) ||
         (TDC_TOF1D[i]>=TDC_TOF1_min && TDC_TOF1D[i]<=TDC_TOF1_max))  Good_TOF1_TDC[i] = true;

      if(Good_TOF1_ADC[i] && Good_TOF1_TDC[i]) Good_TOF1[i] = true;


      if(ADC_TOF2AO[i]>=0 || ADC_TOF2AI[i]>=0 || ADC_TOF2BO[i]>=0 || ADC_TOF2BI[i]>=0)  Good_TOF2_ADC[i] = true;

      if((TDC_TOF2AO[i]>=TDC_TOF2_min && TDC_TOF2AO[i]<=TDC_TOF2_max) ||
         (TDC_TOF2AI[i]>=TDC_TOF2_min && TDC_TOF2AI[i]<=TDC_TOF2_max) ||
         (TDC_TOF2BO[i]>=TDC_TOF2_min && TDC_TOF2BO[i]<=TDC_TOF2_max) ||
         (TDC_TOF2BI[i]>=TDC_TOF2_min && TDC_TOF2BI[i]<=TDC_TOF2_max))  Good_TOF2_TDC[i] = true;

      if(Good_TOF2_ADC[i] && Good_TOF2_TDC[i]) Good_TOF2[i] = true;
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

    /*
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
      if((k>=288 && k<=295) || (k>=304 && k<=311) || (k>=321 && k<=447)){
        if(MWPCADC[k]>=0) count_C4X++;
      }
    }

    for(int kk=256; kk<288; kk++){
      if(MWPCADC[kk]>=0)  count_C4Y++;
    }

   
    if(count_C2X>0 && count_C2Y>0 && count_C3X>0 && count_C3Y>0 && count_C4X>0 && count_C4Y>0)  Good_MWPC_Event = true;
    */

    bool Good_MWPC_Event = true;

    if(Good_TARGET_Event && Good_TOF_Event && Good_MWPC_Event)  Good_Event = true;

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
  }

  cout << "  " << endl;
  cout << "  " << endl;

  if(!Good_Event) return;

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


  TLine *vblue1 = new TLine(0.38,0.70,0.38,0.74);      vblue1->SetLineWidth(7);   vblue1->SetLineColor(4);
  TLine *vblue2 = new TLine(0.42,0.62,0.42,0.66);      vblue2->SetLineWidth(7);   vblue2->SetLineColor(4);    
  TLine *vblue3 = new TLine(0.30,0.54,0.30,0.58);      vblue3->SetLineWidth(7);   vblue3->SetLineColor(4);    
  TLine *vblue4 = new TLine(0.86,0.50,0.86,0.54);      vblue4->SetLineWidth(7);   vblue4->SetLineColor(4);    
  TLine *vblue5 = new TLine(0.70,0.42,0.70,0.46);      vblue5->SetLineWidth(7);   vblue5->SetLineColor(4);    
  TLine *vblue6 = new TLine(0.58,0.34,0.58,0.38);      vblue6->SetLineWidth(7);   vblue6->SetLineColor(4);    
  TLine *vblue7 = new TLine(0.62,0.26,0.62,0.30);      vblue7->SetLineWidth(7);   vblue7->SetLineColor(4);    
  TLine *vblue8 = new TLine(0.62,0.14,0.62,0.18);      vblue8->SetLineWidth(7);   vblue8->SetLineColor(4);

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
  TLatex *tex_Legend_TARGET[18];
 
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

  for(Int_t jj=0; jj<18; jj++){
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
    for (Int_t qq=0; qq<4; qq++) {
      if (tdc_le_sft[ii][qq] > TDC_min_SFT && tdc_le_sft[ii][qq] < TDC_max_SFT) has_TDC_SFT_hit[ii]++;
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

  for(Int_t i=0; i<256; i++){
    for (Int_t k=0; k<4; k++) {
      if ((tdc_le_target[i][k]>=TDC_min_TARGET) && (tdc_le_target[i][k]<=TDC_max_TARGET)) has_TDC_hit[i] = true;
    }
  }

  char ch_ADC_cut_TARGET[100];    sprintf(ch_ADC_cut_TARGET,"(ADC #geq %d)",ADC_cut_TARGET);
  char ch_ADC_and_TDC_cut[100];   sprintf(ch_ADC_and_TDC_cut,"(ADC #geq %d | %d #leq TDC #leq %d)",ADC_cut_TARGET,TDC_min_TARGET,TDC_max_TARGET);

  Int_t Angle_ADC_cut = 0;

  bool has_data_TDC = false;

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
    if(ADC_Low_TARGET[q]>max_ADC) {
      max_index = q;
      max_ADC = ADC_Low_TARGET[q];
    }
  }

  for(Int_t q=0; q<256; q++){
    if (q == max_index) continue;   
    else {
      if(ADC_Low_TARGET[q]>max_ADC2) {
        max_index2 = q;
        max_ADC2 = ADC_Low_TARGET[q];
      }
    }
  }

  for(Int_t q=0; q<256; q++){
    if ((q == max_index) || (q == max_index2)) continue;    
    else {
      if(ADC_Low_TARGET[q]>max_ADC3) {
        max_index3 = q;
        max_ADC3 = ADC_Low_TARGET[q];
      }
    }
  }

  for(Int_t q=0; q<256; q++){
    if ((q == max_index) || (q == max_index2) || (q == max_index3)) continue;   
    else {
      if(ADC_Low_TARGET[q]>max_ADC4) {
        max_index4 = q;
        max_ADC4 = ADC_Low_TARGET[q];
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

  ///Set TOF1 Lines

  if (ADC_TOF1[0]>0 || ADC_TOF1[12]>0) {has_ADC_TOF1_hit[0] = true;}
  if ((TDC_TOF1U[0]>TDC_TOF1_min && TDC_TOF1U[0]<TDC_TOF1_max) || (TDC_TOF1D[0]>TDC_TOF1_min && TDC_TOF1D[0]<TDC_TOF1_max)) {has_TDC_TOF1_hit[0] = true;}
  if (has_ADC_TOF1_hit[0]) TOF_line1->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[0] && has_ADC_TOF1_hit[0]) TOF_line1->SetLineColor(3);

  if (ADC_TOF1[1]>0 || ADC_TOF1[13]>0) {has_ADC_TOF1_hit[1] = true;}
  if ((TDC_TOF1U[1]>TDC_TOF1_min && TDC_TOF1U[1]<TDC_TOF1_max) || (TDC_TOF1D[1]>TDC_TOF1_min && TDC_TOF1D[1]<TDC_TOF1_max)) {has_TDC_TOF1_hit[1] = true;}
  if (has_ADC_TOF1_hit[1]) TOF_line2->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[1] && has_ADC_TOF1_hit[1]) TOF_line2->SetLineColor(3);

  if (ADC_TOF1[2]>0 || ADC_TOF1[14]>0) {has_ADC_TOF1_hit[2] = true;}
  if ((TDC_TOF1U[2]>TDC_TOF1_min && TDC_TOF1U[2]<TDC_TOF1_max) || (TDC_TOF1D[2]>TDC_TOF1_min && TDC_TOF1D[2]<TDC_TOF1_max)) {has_TDC_TOF1_hit[2] = true;}
  if (has_ADC_TOF1_hit[2]) TOF_line3->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[2] && has_ADC_TOF1_hit[2]) TOF_line3->SetLineColor(3);

  if (ADC_TOF1[3]>0 || ADC_TOF1[15]>0) {has_ADC_TOF1_hit[3] = true;}
  if ((TDC_TOF1U[3]>TDC_TOF1_min && TDC_TOF1U[3]<TDC_TOF1_max) || (TDC_TOF1D[3]>TDC_TOF1_min && TDC_TOF1D[3]<TDC_TOF1_max)) {has_TDC_TOF1_hit[3] = true;}
  if (has_ADC_TOF1_hit[3]) TOF_line4->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[3] && has_ADC_TOF1_hit[3]) TOF_line4->SetLineColor(3);

  if (ADC_TOF1[4]>0 || ADC_TOF1[16]>0) {has_ADC_TOF1_hit[4] = true;}
  if ((TDC_TOF1U[4]>TDC_TOF1_min && TDC_TOF1U[4]<TDC_TOF1_max) || (TDC_TOF1D[4]>TDC_TOF1_min && TDC_TOF1D[4]<TDC_TOF1_max)) {has_TDC_TOF1_hit[4] = true;}
  if (has_ADC_TOF1_hit[4]) TOF_line5->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[4] && has_ADC_TOF1_hit[4]) TOF_line5->SetLineColor(3);

  if (ADC_TOF1[5]>0 || ADC_TOF1[17]>0) {has_ADC_TOF1_hit[5] = true;}
  if ((TDC_TOF1U[5]>TDC_TOF1_min && TDC_TOF1U[5]<TDC_TOF1_max) || (TDC_TOF1D[5]>TDC_TOF1_min && TDC_TOF1D[5]<TDC_TOF1_max)) {has_TDC_TOF1_hit[5] = true;}
  if (has_ADC_TOF1_hit[5]) TOF_line6->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[5] && has_ADC_TOF1_hit[5]) TOF_line6->SetLineColor(3);

  if (ADC_TOF1[6]>0 || ADC_TOF1[18]>0) {has_ADC_TOF1_hit[6] = true;}
  if ((TDC_TOF1U[6]>TDC_TOF1_min && TDC_TOF1U[6]<TDC_TOF1_max) || (TDC_TOF1D[6]>TDC_TOF1_min && TDC_TOF1D[6]<TDC_TOF1_max)) {has_TDC_TOF1_hit[6] = true;}
  if (has_ADC_TOF1_hit[6]) TOF_line7->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[6] && has_ADC_TOF1_hit[6]) TOF_line7->SetLineColor(3);

  if (ADC_TOF1[7]>0 || ADC_TOF1[19]>0) {has_ADC_TOF1_hit[7] = true;}
  if ((TDC_TOF1U[7]>TDC_TOF1_min && TDC_TOF1U[7]<TDC_TOF1_max) || (TDC_TOF1D[7]>TDC_TOF1_min && TDC_TOF1D[7]<TDC_TOF1_max)) {has_TDC_TOF1_hit[7] = true;}
  if (has_ADC_TOF1_hit[7]) TOF_line8->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[7] && has_ADC_TOF1_hit[7]) TOF_line8->SetLineColor(3);

  if (ADC_TOF1[8]>0 || ADC_TOF1[20]>0) {has_ADC_TOF1_hit[8] = true;}
  if ((TDC_TOF1U[8]>TDC_TOF1_min && TDC_TOF1U[8]<TDC_TOF1_max) || (TDC_TOF1D[8]>TDC_TOF1_min && TDC_TOF1D[8]<TDC_TOF1_max)) {has_TDC_TOF1_hit[8] = true;}
  if (has_ADC_TOF1_hit[8]) TOF_line9->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[8] && has_ADC_TOF1_hit[8]) TOF_line9->SetLineColor(3);

  if (ADC_TOF1[9]>0 || ADC_TOF1[21]>0) {has_ADC_TOF1_hit[9] = true;}
  if ((TDC_TOF1U[9]>TDC_TOF1_min && TDC_TOF1U[9]<TDC_TOF1_max) || (TDC_TOF1D[9]>TDC_TOF1_min && TDC_TOF1D[9]<TDC_TOF1_max)) {has_TDC_TOF1_hit[9] = true;}
  if (has_ADC_TOF1_hit[9]) TOF_line10->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[9] && has_ADC_TOF1_hit[9]) TOF_line10->SetLineColor(3);

  if (ADC_TOF1[10]>0 || ADC_TOF1[22]>0) {has_ADC_TOF1_hit[10] = true;}
  if ((TDC_TOF1U[10]>TDC_TOF1_min && TDC_TOF1U[10]<TDC_TOF1_max) || (TDC_TOF1D[10]>TDC_TOF1_min && TDC_TOF1D[10]<TDC_TOF1_max)) {has_TDC_TOF1_hit[10] = true;}
  if (has_ADC_TOF1_hit[10]) TOF_line11->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[10] && has_ADC_TOF1_hit[10]) TOF_line11->SetLineColor(3);

  if (ADC_TOF1[11]>0 || ADC_TOF1[23]>0) {has_ADC_TOF1_hit[11] = true;}
  if ((TDC_TOF1U[11]>TDC_TOF1_min && TDC_TOF1U[11]<TDC_TOF1_max) || (TDC_TOF1D[11]>TDC_TOF1_min && TDC_TOF1D[11]<TDC_TOF1_max)) {has_TDC_TOF1_hit[11] = true;}
  if (has_ADC_TOF1_hit[11]) TOF_line12->SetLineColor(kOrange+10);
  if (has_TDC_TOF1_hit[11] && has_ADC_TOF1_hit[11]) TOF_line12->SetLineColor(3);

  ///Set TOF2 Lines

  if ((ADC_TOF2[0]>0 || ADC_TOF2[24]>0) || (ADC_TOF2[12]>0 || ADC_TOF2[36]>0)) {has_ADC_TOF2_hit[0]=true;}
  if (((TDC_TOF2AO[0]>TDC_TOF2_min && TDC_TOF2AO[0] < TDC_TOF2_max) || (TDC_TOF2AI[0]>TDC_TOF2_min && TDC_TOF2AI[0] < TDC_TOF2_max)) 
    || ((TDC_TOF2BO[0]>TDC_TOF2_min && TDC_TOF2BO[0] < TDC_TOF2_max) || (TDC_TOF2BI[0]>TDC_TOF2_min && TDC_TOF2BI[0] < TDC_TOF2_max))) {has_TDC_TOF2_hit[0]=true;}
  if (has_ADC_TOF2_hit[0]) TOF_line13->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[0] && has_ADC_TOF2_hit[0]) TOF_line13->SetLineColor(3);

  if ((ADC_TOF2[1]>0 || ADC_TOF2[25]>0) || (ADC_TOF2[13]>0 || ADC_TOF2[37]>0)) {has_ADC_TOF2_hit[1]=true;}
  if (((TDC_TOF2AO[1]>TDC_TOF2_min && TDC_TOF2AO[1] < TDC_TOF2_max) || (TDC_TOF2AI[1]>TDC_TOF2_min && TDC_TOF2AI[1] < TDC_TOF2_max)) 
    || ((TDC_TOF2BO[1]>TDC_TOF2_min && TDC_TOF2BO[1] < TDC_TOF2_max) || (TDC_TOF2BI[1]>TDC_TOF2_min && TDC_TOF2BI[1] < TDC_TOF2_max))) {has_TDC_TOF2_hit[1]=true;}
  if (has_ADC_TOF2_hit[1]) TOF_line14->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[1] && has_ADC_TOF2_hit[1]) TOF_line14->SetLineColor(3);

  if ((ADC_TOF2[2]>0 || ADC_TOF2[26]>0) || (ADC_TOF2[14]>0 || ADC_TOF2[38]>0)) {has_ADC_TOF2_hit[2]=true;}
  if (((TDC_TOF2AO[2]>TDC_TOF2_min && TDC_TOF2AO[2] < TDC_TOF2_max) || (TDC_TOF2AI[2]>TDC_TOF2_min && TDC_TOF2AI[2] < TDC_TOF2_max)) 
    || ((TDC_TOF2BO[2]>TDC_TOF2_min && TDC_TOF2BO[2] < TDC_TOF2_max) || (TDC_TOF2BI[2]>TDC_TOF2_min && TDC_TOF2BI[2] < TDC_TOF2_max))) {has_TDC_TOF2_hit[2]=true;}
  if (has_ADC_TOF2_hit[2]) TOF_line15->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[2] && has_ADC_TOF2_hit[2]) TOF_line15->SetLineColor(3);

  if ((ADC_TOF2[3]>0 || ADC_TOF2[27]>0) || (ADC_TOF2[15]>0 || ADC_TOF2[39]>0)) {has_ADC_TOF2_hit[3]=true;}
  if (((TDC_TOF2AO[3]>TDC_TOF2_min && TDC_TOF2AO[3] < TDC_TOF2_max) || (TDC_TOF2AI[3]>TDC_TOF2_min && TDC_TOF2AI[3] < TDC_TOF2_max)) 
    || ((TDC_TOF2BO[3]>TDC_TOF2_min && TDC_TOF2BO[3] < TDC_TOF2_max) || (TDC_TOF2BI[3]>TDC_TOF2_min && TDC_TOF2BI[3] < TDC_TOF2_max))) {has_TDC_TOF2_hit[3]=true;}
  if (has_ADC_TOF2_hit[3]) TOF_line16->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[3] && has_ADC_TOF2_hit[3]) TOF_line16->SetLineColor(3);

  if ((ADC_TOF2[4]>0 || ADC_TOF2[28]>0) || (ADC_TOF2[16]>0 || ADC_TOF2[40]>0)) {has_ADC_TOF2_hit[4]=true;}
  if (((TDC_TOF2AO[4]>TDC_TOF2_min && TDC_TOF2AO[4] < TDC_TOF2_max) || (TDC_TOF2AI[4]>TDC_TOF2_min && TDC_TOF2AI[4] < TDC_TOF2_max)) 
    || ((TDC_TOF2BO[4]>TDC_TOF2_min && TDC_TOF2BO[4] < TDC_TOF2_max) || (TDC_TOF2BI[4]>TDC_TOF2_min && TDC_TOF2BI[4] < TDC_TOF2_max))) {has_TDC_TOF2_hit[4]=true;}
  if (has_ADC_TOF2_hit[4]) TOF_line17->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[4] && has_ADC_TOF2_hit[4]) TOF_line17->SetLineColor(3);

  if ((ADC_TOF2[5]>0 || ADC_TOF2[29]>0) || (ADC_TOF2[17]>0 || ADC_TOF2[41]>0)) {has_ADC_TOF2_hit[5]=true;}
  if (((TDC_TOF2AO[5]>TDC_TOF2_min && TDC_TOF2AO[5] < TDC_TOF2_max) || (TDC_TOF2AI[5]>TDC_TOF2_min && TDC_TOF2AI[5] < TDC_TOF2_max)) 
    || ((TDC_TOF2BO[5]>TDC_TOF2_min && TDC_TOF2BO[5] < TDC_TOF2_max) || (TDC_TOF2BI[5]>TDC_TOF2_min && TDC_TOF2BI[5] < TDC_TOF2_max))) {has_TDC_TOF2_hit[5]=true;}
  if (has_ADC_TOF2_hit[5]) TOF_line18->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[5] && has_ADC_TOF2_hit[5]) TOF_line18->SetLineColor(3);

  if ((ADC_TOF2[6]>0 || ADC_TOF2[30]>0) || (ADC_TOF2[55]>0 || ADC_TOF2[42]>0)) {has_ADC_TOF2_hit[6]=true;}
  if (((TDC_TOF2AO[6]>TDC_TOF2_min && TDC_TOF2AO[6] < TDC_TOF2_max) || (TDC_TOF2AI[6]>TDC_TOF2_min && TDC_TOF2AI[6] < TDC_TOF2_max)) 
    || ((TDC_TOF2BO[6]>TDC_TOF2_min && TDC_TOF2BO[6] < TDC_TOF2_max) || (TDC_TOF2BI[6]>TDC_TOF2_min && TDC_TOF2BI[6] < TDC_TOF2_max))) {has_TDC_TOF2_hit[6]=true;}
  if (has_ADC_TOF2_hit[6]) TOF_line19->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[6] && has_ADC_TOF2_hit[6]) TOF_line19->SetLineColor(3);

  if ((ADC_TOF2[7]>0 || ADC_TOF2[31]>0) || (ADC_TOF2[19]>0 || ADC_TOF2[43]>0)) {has_ADC_TOF2_hit[7]=true;}
  if (((TDC_TOF2AO[7]>TDC_TOF2_min && TDC_TOF2AO[7] < TDC_TOF2_max) || (TDC_TOF2AI[7]>TDC_TOF2_min && TDC_TOF2AI[7] < TDC_TOF2_max)) 
    || ((TDC_TOF2BO[7]>TDC_TOF2_min && TDC_TOF2BO[7] < TDC_TOF2_max) || (TDC_TOF2BI[7]>TDC_TOF2_min && TDC_TOF2BI[7] < TDC_TOF2_max))) {has_TDC_TOF2_hit[7]=true;}
  if (has_ADC_TOF2_hit[7]) TOF_line20->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[7] && has_ADC_TOF2_hit[7]) TOF_line20->SetLineColor(3);

  if ((ADC_TOF2[8]>0 || ADC_TOF2[32]>0) || (ADC_TOF2[20]>0 || ADC_TOF2[44]>0)) {has_ADC_TOF2_hit[8]=true;}
  if (((TDC_TOF2AO[8]>TDC_TOF2_min && TDC_TOF2AO[8] < TDC_TOF2_max) || (TDC_TOF2AI[8]>TDC_TOF2_min && TDC_TOF2AI[8] < TDC_TOF2_max)) 
    || ((TDC_TOF2BO[8]>TDC_TOF2_min && TDC_TOF2BO[8] < TDC_TOF2_max) || (TDC_TOF2BI[8]>TDC_TOF2_min && TDC_TOF2BI[8] < TDC_TOF2_max))) {has_TDC_TOF2_hit[8]=true;}
  if (has_ADC_TOF2_hit[8]) TOF_line21->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[8] && has_ADC_TOF2_hit[8]) TOF_line21->SetLineColor(3);

  if ((ADC_TOF2[9]>0 || ADC_TOF2[33]>0) || (ADC_TOF2[21]>0 || ADC_TOF2[45]>0)) {has_ADC_TOF2_hit[9]=true;}
  if (((TDC_TOF2AO[9]>TDC_TOF2_min && TDC_TOF2AO[9] < TDC_TOF2_max) || (TDC_TOF2AI[9]>TDC_TOF2_min && TDC_TOF2AI[9] < TDC_TOF2_max)) 
    || ((TDC_TOF2BO[9]>TDC_TOF2_min && TDC_TOF2BO[9] < TDC_TOF2_max) || (TDC_TOF2BI[9]>TDC_TOF2_min && TDC_TOF2BI[9] < TDC_TOF2_max))) {has_TDC_TOF2_hit[9]=true;}
  if (has_ADC_TOF2_hit[9]) TOF_line22->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[9] && has_ADC_TOF2_hit[9]) TOF_line22->SetLineColor(3);

  if ((ADC_TOF2[10]>0 || ADC_TOF2[34]>0) || (ADC_TOF2[22]>0 || ADC_TOF2[46]>0)) {has_ADC_TOF2_hit[10]=true;}
  if (((TDC_TOF2AO[10]>TDC_TOF2_min && TDC_TOF2AO[10] < TDC_TOF2_max) || (TDC_TOF2AI[10]>TDC_TOF2_min && TDC_TOF2AI[10] < TDC_TOF2_max)) 
    || ((TDC_TOF2BO[10]>TDC_TOF2_min && TDC_TOF2BO[10] < TDC_TOF2_max) || (TDC_TOF2BI[10]>TDC_TOF2_min && TDC_TOF2BI[10] < TDC_TOF2_max))) {has_TDC_TOF2_hit[10]=true;}
  if (has_ADC_TOF2_hit[10]) TOF_line23->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[10] && has_ADC_TOF2_hit[10]) TOF_line23->SetLineColor(3);

  if ((ADC_TOF2[11]>0 || ADC_TOF2[35]>0) || (ADC_TOF2[23]>0 || ADC_TOF2[47]>0)) {has_ADC_TOF2_hit[11]=true;}
  if (((TDC_TOF2AO[11]>TDC_TOF2_min && TDC_TOF2AO[11] < TDC_TOF2_max) || (TDC_TOF2AI[11]>TDC_TOF2_min && TDC_TOF2AI[11] < TDC_TOF2_max)) 
    || ((TDC_TOF2BO[11]>TDC_TOF2_min && TDC_TOF2BO[11] < TDC_TOF2_max) || (TDC_TOF2BI[11]>TDC_TOF2_min && TDC_TOF2BI[11] < TDC_TOF2_max))) {has_TDC_TOF2_hit[11]=true;}
  if (has_ADC_TOF2_hit[11]) TOF_line24->SetLineColor(kOrange+10);
  if (has_TDC_TOF2_hit[11] && has_ADC_TOF2_hit[11]) TOF_line24->SetLineColor(3);

  int gap_hit[12] = {0};
  int ADC_TOF1_hit[12] = {0};
  int ADCTDC_TOF1_hit[12] = {0};
  int ADC_TOF2_hit[12] = {0};
  int ADCTDC_TOF2_hit[12] = {0};

  //// Determine which gaps have ADC & TDC TARGET hits

	//int channel[12][8] = {{0}};

  //for(int i=0; i<12; i++){
  //  for(int j=0; j<8; j++){
	//	  channel[i][j]=Channel[j+8*i];
	//  }
  //}		
	 


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


  if(scoring_type==2){
    for(int i=0; i<12; i++){
      if(ADC_TOF1U[i]>=0) gap_counter[i]++;
    }

    for(int i=0; i<12; i++){
      if(ADC_TOF1D[i]>=0) gap_counter[i]++;
    }

    for(int i=0; i<12; i++){
      if(TDC_TOF1U[i]>=TDC_TOF1_min && TDC_TOF1U[i]<=TDC_TOF1_max) gap_counter[i]++;
    }

    for(int i=0; i<12; i++){
      if(TDC_TOF1D[i]>=TDC_TOF1_min && TDC_TOF1D[i]<=TDC_TOF1_max) gap_counter[i]++;
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

        if(TDC_TOF2AO[i]>=TDC_TOF2_min && TDC_TOF2AO[i]<=TDC_TOF2_max) gap_counter[i]++;   
        if(TDC_TOF2AO[i-1]>=TDC_TOF2_min && TDC_TOF2AO[i-1]<=TDC_TOF2_max) gap_counter[i]++;   
        if(TDC_TOF2AO[i+1]>=TDC_TOF2_min && TDC_TOF2AO[i+1]<=TDC_TOF2_max) gap_counter[i]++;   

        if(TDC_TOF2AI[i]>=TDC_TOF2_min && TDC_TOF2AI[i]<=TDC_TOF2_max) gap_counter[i]++;   
        if(TDC_TOF2AI[i-1]>=TDC_TOF2_min && TDC_TOF2AI[i-1]<=TDC_TOF2_max) gap_counter[i]++;   
        if(TDC_TOF2AI[i+1]>=TDC_TOF2_min && TDC_TOF2AI[i+1]<=TDC_TOF2_max) gap_counter[i]++;   

        if(TDC_TOF2BO[i]>=TDC_TOF2_min && TDC_TOF2BO[i]<=TDC_TOF2_max) gap_counter[i]++;   
        if(TDC_TOF2BO[i-1]>=TDC_TOF2_min && TDC_TOF2BO[i-1]<=TDC_TOF2_max) gap_counter[i]++;   
        if(TDC_TOF2BO[i+1]>=TDC_TOF2_min && TDC_TOF2BO[i+1]<=TDC_TOF2_max) gap_counter[i]++;   

        if(TDC_TOF2BI[i]>=TDC_TOF2_min && TDC_TOF2BI[i]<=TDC_TOF2_max) gap_counter[i]++;   
        if(TDC_TOF2BI[i-1]>=TDC_TOF2_min && TDC_TOF2BI[i-1]<=TDC_TOF2_max) gap_counter[i]++;   
        if(TDC_TOF2BI[i+1]>=TDC_TOF2_min && TDC_TOF2BI[i+1]<=TDC_TOF2_max) gap_counter[i]++;   
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

      if(TDC_TOF2AO[0]>=TDC_TOF2_min && TDC_TOF2AO[0]<=TDC_TOF2_max) gap_counter[0]++;   
      if(TDC_TOF2AO[1]>=TDC_TOF2_min && TDC_TOF2AO[1]<=TDC_TOF2_max) gap_counter[0]++;   
      if(TDC_TOF2AO[11]>=TDC_TOF2_min && TDC_TOF2AO[11]<=TDC_TOF2_max) gap_counter[0]++;   

      if(TDC_TOF2AI[0]>=TDC_TOF2_min && TDC_TOF2AI[0]<=TDC_TOF2_max) gap_counter[0]++;   
      if(TDC_TOF2AI[1]>=TDC_TOF2_min && TDC_TOF2AI[1]<=TDC_TOF2_max) gap_counter[0]++;   
      if(TDC_TOF2AI[11]>=TDC_TOF2_min && TDC_TOF2AI[11]<=TDC_TOF2_max) gap_counter[0]++;   

      if(TDC_TOF2BO[0]>=TDC_TOF2_min && TDC_TOF2BO[0]<=TDC_TOF2_max) gap_counter[0]++;   
      if(TDC_TOF2BO[1]>=TDC_TOF2_min && TDC_TOF2BO[1]<=TDC_TOF2_max) gap_counter[0]++;   
      if(TDC_TOF2BO[11]>=TDC_TOF2_min && TDC_TOF2BO[11]<=TDC_TOF2_max) gap_counter[0]++;   

      if(TDC_TOF2BI[0]>=TDC_TOF2_min && TDC_TOF2BI[0]<=TDC_TOF2_max) gap_counter[0]++;   
      if(TDC_TOF2BI[1]>=TDC_TOF2_min && TDC_TOF2BI[1]<=TDC_TOF2_max) gap_counter[0]++;   
      if(TDC_TOF2BI[11]>=TDC_TOF2_min && TDC_TOF2BI[11]<=TDC_TOF2_max) gap_counter[0]++;   
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

      if(TDC_TOF2AO[11]>=TDC_TOF2_min && TDC_TOF2AO[11]<=TDC_TOF2_max) gap_counter[11]++;   
      if(TDC_TOF2AO[10]>=TDC_TOF2_min && TDC_TOF2AO[10]<=TDC_TOF2_max) gap_counter[11]++;   
      if(TDC_TOF2AO[0]>=TDC_TOF2_min && TDC_TOF2AO[0]<=TDC_TOF2_max) gap_counter[11]++;   

      if(TDC_TOF2AI[11]>=TDC_TOF2_min && TDC_TOF2AI[11]<=TDC_TOF2_max) gap_counter[11]++;   
      if(TDC_TOF2AI[10]>=TDC_TOF2_min && TDC_TOF2AI[10]<=TDC_TOF2_max) gap_counter[11]++;   
      if(TDC_TOF2AI[0]>=TDC_TOF2_min && TDC_TOF2AI[0]<=TDC_TOF2_max) gap_counter[11]++;   

      if(TDC_TOF2BO[11]>=TDC_TOF2_min && TDC_TOF2BO[11]<=TDC_TOF2_max) gap_counter[11]++;   
      if(TDC_TOF2BO[10]>=TDC_TOF2_min && TDC_TOF2BO[10]<=TDC_TOF2_max) gap_counter[11]++;   
      if(TDC_TOF2BO[0]>=TDC_TOF2_min && TDC_TOF2BO[0]<=TDC_TOF2_max) gap_counter[11]++;   

      if(TDC_TOF2BI[11]>=TDC_TOF2_min && TDC_TOF2BI[11]<=TDC_TOF2_max) gap_counter[11]++;   
      if(TDC_TOF2BI[10]>=TDC_TOF2_min && TDC_TOF2BI[10]<=TDC_TOF2_max) gap_counter[11]++;   
      if(TDC_TOF2BI[0]>=TDC_TOF2_min && TDC_TOF2BI[0]<=TDC_TOF2_max) gap_counter[11]++;   
    } 

    for(int i=0; i<12; i++)
    {
      if(Good_tof1[i] && Good_tof2[i]) gap_counter[i]++;
    }
  }

  //Sebastien's Method
  if(scoring_type==3){
    for(int i=0; i<12; i++){
      if(ADC_TOF2AO[i]>=0) gap_counter[i]++;
    }

    for(int i=0; i<12; i++){
      if(ADC_TOF2AI[i]>=0) gap_counter[i]++;
    }

    for(int i=0; i<12; i++){
      if(ADC_TOF2BO[i]>=0) gap_counter[i]++;
    }

    for(int i=0; i<12; i++){
      if(ADC_TOF2BI[i]>=0) gap_counter[i]++;
    }


    for(int i=0; i<12; i++){
      if(TDC_TOF2AO[i]>=TDC_TOF2_min && TDC_TOF2AO[i]<=TDC_TOF2_max) gap_counter[i]++;
    }

    for(int i=0; i<12; i++){
      if(TDC_TOF2AI[i]>=TDC_TOF2_min && TDC_TOF2AI[i]<=TDC_TOF2_max) gap_counter[i]++;
    }

    for(int i=0; i<12; i++){
      if(TDC_TOF2BO[i]>=TDC_TOF2_min && TDC_TOF2BO[i]<=TDC_TOF2_max) gap_counter[i]++;
    }

    for(int i=0; i<12; i++){
      if(TDC_TOF2BI[i]>=TDC_TOF2_min && TDC_TOF2BI[i]<=TDC_TOF2_max) gap_counter[i]++;
    }


    for(int i=0; i<12; i++){
      if(i!=0 && i!=11 && gap_counter[i]>0){


      if(ADC_TOF1U[i]>=0) gap_counter[i]++;
      if(ADC_TOF1U[i-1]>=0) gap_counter[i]++;
      if(ADC_TOF1U[i+1]>=0) gap_counter[i]++;

      if(ADC_TOF1D[i]>=0) gap_counter[i]++;
      if(ADC_TOF1D[i-1]>=0) gap_counter[i]++;
      if(ADC_TOF1D[i+1]>=0) gap_counter[i]++;   

      if(TDC_TOF1U[i]>=TDC_TOF1_min && TDC_TOF1U[i]<=TDC_TOF1_max) gap_counter[i]++;   
      if(TDC_TOF1U[i-1]>=TDC_TOF1_min && TDC_TOF1U[i-1]<=TDC_TOF1_max) gap_counter[i]++;   
      if(TDC_TOF1U[i+1]>=TDC_TOF1_min && TDC_TOF1U[i+1]<=TDC_TOF1_max) gap_counter[i]++;   

      if(TDC_TOF1D[i]>=TDC_TOF1_min && TDC_TOF1D[i]<=TDC_TOF1_max) gap_counter[i]++;   
      if(TDC_TOF1D[i-1]>=TDC_TOF1_min && TDC_TOF1D[i-1]<=TDC_TOF1_max) gap_counter[i]++;   
      if(TDC_TOF1D[i+1]>=TDC_TOF1_min && TDC_TOF1D[i+1]<=TDC_TOF1_max) gap_counter[i]++;   
   }
  }

  if(gap_counter[0]>0){
    if(ADC_TOF1U[0]>=0) gap_counter[0]++;
    if(ADC_TOF1U[1]>=0) gap_counter[0]++;
    if(ADC_TOF1U[11]>=0) gap_counter[0]++;

    if(ADC_TOF1D[0]>=0) gap_counter[0]++;
    if(ADC_TOF1D[1]>=0) gap_counter[0]++;
    if(ADC_TOF1D[11]>=0) gap_counter[0]++;   

    if(TDC_TOF1U[0]>=TDC_TOF1_min && TDC_TOF1U[0]<=TDC_TOF1_max) gap_counter[0]++;   
    if(TDC_TOF1U[1]>=TDC_TOF1_min && TDC_TOF1U[1]<=TDC_TOF1_max) gap_counter[0]++;   
    if(TDC_TOF1U[11]>=TDC_TOF1_min && TDC_TOF1U[11]<=TDC_TOF1_max) gap_counter[0]++;   

    if(TDC_TOF1D[0]>=TDC_TOF1_min && TDC_TOF1D[0]<=TDC_TOF1_max) gap_counter[0]++;   
    if(TDC_TOF1D[1]>=TDC_TOF1_min && TDC_TOF1D[1]<=TDC_TOF1_max) gap_counter[0]++;   
    if(TDC_TOF1D[11]>=TDC_TOF1_min && TDC_TOF1D[11]<=TDC_TOF1_max) gap_counter[0]++;   
  }

  if(gap_counter[11]>0){
    if(ADC_TOF1U[11]>=0) gap_counter[11]++;
    if(ADC_TOF1U[10]>=0) gap_counter[11]++;
    if(ADC_TOF1U[0]>=0) gap_counter[11]++;

    if(ADC_TOF1D[11]>=0) gap_counter[11]++;
    if(ADC_TOF1D[10]>=0) gap_counter[11]++;
    if(ADC_TOF1D[0]>=0) gap_counter[11]++;   

    if(TDC_TOF1U[11]>=TDC_TOF1_min && TDC_TOF1U[11]<=TDC_TOF1_max) gap_counter[11]++;   
    if(TDC_TOF1U[10]>=TDC_TOF1_min && TDC_TOF1U[10]<=TDC_TOF1_max) gap_counter[11]++;   
    if(TDC_TOF1U[0]>=TDC_TOF1_min && TDC_TOF1U[0]<=TDC_TOF1_max) gap_counter[11]++;   

    if(TDC_TOF1D[11]>=TDC_TOF1_min && TDC_TOF1D[11]<=TDC_TOF1_max) gap_counter[11]++;   
    if(TDC_TOF1D[10]>=TDC_TOF1_min && TDC_TOF1D[10]<=TDC_TOF1_max) gap_counter[11]++;   
    if(TDC_TOF1D[0]>=TDC_TOF1_min && TDC_TOF1D[0]<=TDC_TOF1_max) gap_counter[11]++;   
  }
  }



  
  //Oscar's method
  if(scoring_type==1){
  for (int k=0; k<12; k++) {
    if ((has_ADC_TOF2_hit[k] && has_TDC_TOF2_hit[k]) && (has_ADC_TOF1_hit[k] && has_TDC_TOF1_hit[k])) gap_counter[k] = gap_counter[k] + 5;
    if ((has_ADC_TOF2_hit[k] || (has_ADC_TOF2_hit[k] && has_TDC_TOF2_hit[k])) && (has_ADC_TOF1_hit[k] || (has_ADC_TOF1_hit[k] && has_TDC_TOF1_hit[k]))) gap_counter[k] = gap_counter[k] + 2;
    if (gap_hit[k] > 0) gap_counter[k] = gap_counter[k] + 4;
    if (ADCTDC_TOF1_hit[k] > 0) gap_counter[k] = gap_counter[k] + 2;
    if (ADCTDC_TOF2_hit[k] > 0) gap_counter[k] = gap_counter[k] + 2;
    if (ADC_TOF2_hit[k] > 0) gap_counter[k]++;
    if (ADC_TOF1_hit[k] > 0) gap_counter[k]++;
  }
  }


  int max_gap[12] = {0};
  int high_gap_hit = 0;
  int gap_to_fit = 0;

  for (int k=0; k<12; k++) {
    if (gap_counter[k] >= high_gap_hit) {
      high_gap_hit = gap_counter[k];
      max_gap[k] = max_gap[k]+gap_counter[k];
      gap_to_fit = k+1;
      }
  }

  ////////////////////// ADDITIONAL

  for(int i=0; i<8; i++){
    cout << "TEST TEST TEST:  " << channel[gap_to_fit-1][i] << endl;
    //h_TDC_Gap_Fibers->Fill(Xloc[channel[gap_to_fit-1][i]], Yloc[channel[gap_to_fit-1][i]]);
//    h_TDC_selected2->Fill(Xloc[channel[gap_to_fit-1][i]], Yloc[channel[gap_to_fit-1][i]]);
  }

  for(Int_t i=0; i<256; i++){
    if(ADC_High_TARGET[i]>=Angle_ADC_cut && has_TDC_hit[i]){
      h_target_TDC->Fill(Xloc[i],Yloc[i]);
      h_target_TDC_copy->Fill(Xloc[i],Yloc[i]);
      count++;
      has_data_TDC = true;
      if(ADC_High_TARGET[i]>=300){
        h_TDC_selected->Fill(Xloc[i],Yloc[i]);
        h_TDC_selected2->Fill(Xloc[i],Yloc[i]);
        //h_TDC_Gap_Fibers->Fill(Xloc[i],Yloc[i]);


        if(IsIn(i,channel[gap_to_fit-1][0], channel[gap_to_fit-1][1], 
                  channel[gap_to_fit-1][2], channel[gap_to_fit-1][3], 
                  channel[gap_to_fit-1][4], channel[gap_to_fit-1][5], 
                  channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){
            //h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
            //h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);

            h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);


            h_TDC_selected2->Fill(Xloc[i], Yloc[i]);
            h_TDC_selected2->Fill(Xloc[i], Yloc[i]);
  
            h_TDC_selected->Fill(Xloc[i], Yloc[i]);
            h_TDC_selected->Fill(Xloc[i], Yloc[i]);
        }
      }
    }  

    if(ADC_High_TARGET[i]>=Angle_ADC_cut) {
      h_target_ADC->Fill(Xloc[i],Yloc[i]);
      x_inc = x_inc + Xloc[i];
      y_inc = y_inc + Yloc[i];
      hit_count++;
    }
  }

  //// Don't do anything if the event has less than 5 hits in the TARGET
  if(count<5){
    gROOT->Reset();
    gROOT->Clear();
    cout << " >>>>  Event "<< ievt << " has less than 5 hits in the TARGET !" << endl;
    cout << " >>>>  Please, choose another event" << endl;
    cout << " " << endl;
    return; 
  }


  /// Print Values
  if(Good_Event && (enable_cout==0 || enable_cout==1)){
    cout << "//////  TARGET  //////" << endl;
    cout << "Fiber  HG-1000  LG-1000      TDC[0]      [1]      [2]      [3]" << endl;
    for(Int_t jj=0; jj<256; jj++){
    printf("%3d    %4d      %4d       %4d       %4d     %4d     %4d\n", jj, ADC_High_TARGET[jj], ADC_Low_TARGET[jj], TDC_LE_TARGET[jj], tdc_le_target[jj][1], tdc_le_target[jj][2], tdc_le_target[jj][3]);
    }
  }  
  
  cout << " " << endl;

  if(Good_Event && enable_cout==1){
    cout << " " << endl;
    cout << " " << endl;
    cout << "//////  SFT  //////" << endl;
    cout << "Fiber  HG-950      TDC[0]      [1]      [2]      [3]" << endl;
    for(Int_t jj=0; jj<128; jj++){
      printf("%3d    %4d       %4d       %4d     %4d     %4d\n", jj, ADC_High_SFT[jj], TDC_LE_SFT[jj], tdc_le_sft[jj][1], tdc_le_sft[jj][2], tdc_le_sft[jj][3]);
    }
  }

  if(Good_Event && enable_cout==2){
    cout << "//////  TARGET  //////" << endl;
    cout << "Fiber  HG-1000  LG-1000      TDC[0]      [1]      [2]      [3]" << endl;
    for(Int_t jj=0; jj<256; jj++){
      if((TDC_LE_TARGET[jj] >= TDC_min_TARGET && TDC_LE_TARGET[jj] <= TDC_max_TARGET) || 
        (tdc_le_target[jj][1] >= TDC_min_TARGET && tdc_le_target[jj][1] <= TDC_max_TARGET) || 
        (tdc_le_target[jj][2] >= TDC_min_TARGET && tdc_le_target[jj][2] <= TDC_max_TARGET) || 
        (tdc_le_target[jj][3] >= TDC_min_TARGET && tdc_le_target[jj][3] <= TDC_max_TARGET)){
        printf("%3d    %4d      %4d       %4d       %4d     %4d     %4d\n", jj, ADC_High_TARGET[jj], ADC_Low_TARGET[jj], TDC_LE_TARGET[jj], tdc_le_target[jj][1], tdc_le_target[jj][2], tdc_le_target[jj][3]);
      }
    }

    cout << " " << endl;
    cout << " " << endl;
    cout << "//////  SFT  //////" << endl;
    cout << "Fiber  HG-950      TDC[0]      [1]      [2]      [3]" << endl;
    for(Int_t jj=0; jj<128; jj++){
      if((TDC_LE_SFT[jj] >= TDC_min_SFT && TDC_LE_SFT[jj] <= TDC_max_SFT) || 
        (tdc_le_sft[jj][1] >= TDC_min_SFT && tdc_le_sft[jj][1] <= TDC_max_SFT) ||
        (tdc_le_sft[jj][2] >= TDC_min_SFT && tdc_le_sft[jj][2] <= TDC_max_SFT) || 
        (tdc_le_sft[jj][3] >= TDC_min_SFT && tdc_le_sft[jj][3] <= TDC_max_SFT)){
        printf("%3d    %4d       %4d       %4d     %4d     %4d\n", jj, ADC_High_SFT[jj], TDC_LE_SFT[jj], tdc_le_sft[jj][1], tdc_le_sft[jj][2], tdc_le_sft[jj][3]);
      }
    }
  }

  if(Good_Event && (enable_cout!=0 && enable_cout!=1 && enable_cout!=2)){
    cout << "Flag Error !" << endl;
  }

  cout << "" << endl;

  for (Int_t i = 0; i<12; i++) {
    if(Good_Event){
      if (ADC_TOF1[i] > 0 || ADC_TOF1[i+12] > 0) {
        if ((TDC_TOF1U[i]>TDC_TOF1_min && TDC_TOF1U[i]<TDC_TOF1_max) || (TDC_TOF1D[i]>TDC_TOF1_min && TDC_TOF1D[i]<TDC_TOF1_max)) {
          cout << "ADC TOF1 Up-" << i+1 << ": " << ADC_TOF1[i] << " -- Down-" << i+1 << ": " << ADC_TOF1[i+12] << " | TDC TOF1 Up-" << i+1 << ": " << TDC_TOF1U[i] << " -- Down-" << i+1 << ": " << TDC_TOF1D[i]<< endl;
        }
        else {cout << "ADC TOF1 Up-" << i+1 << ": " << ADC_TOF1[i] << " -- Down-" << i+1 << ": " << ADC_TOF1[i+12]<< endl;}
      }
      else {
        if ((TDC_TOF1U[i]>TDC_TOF1_min && TDC_TOF1U[i]<TDC_TOF1_max) || (TDC_TOF1D[i]>TDC_TOF1_min && TDC_TOF1D[i]<TDC_TOF1_max)) {
          cout << "TDC TOF1 Up-" << i+1 << ": " << TDC_TOF1U[i] << " -- Down-" << i+1 << ": " << TDC_TOF1D[i] << endl;
        }
      }
    }  
  }

  cout << "" << endl;

  for (Int_t i = 0; i<12; i++) {
    if (ADC_TOF2[i] > 0 || ADC_TOF2[i+24] > 0) {
      if ((TDC_TOF2AO[i] > TDC_TOF2_min && TDC_TOF2AO[i] < TDC_TOF2_max) || (TDC_TOF2AI[i] > TDC_TOF2_min && TDC_TOF2AI[i] < TDC_TOF2_max)) {
      cout << "ADC TOF2 OutA-" << i+1 << ": " << ADC_TOF2[i] << " -- InA-" << i+1 << ": " << ADC_TOF2[i+24]<< " | TDC TOF2 OutA-" << i+1 << ": " << TDC_TOF2AO[i] << " -- InA-" << i+1 << ": " << TDC_TOF2AI[i]<< endl;
      }
      else {cout << "ADC TOF2 OutA-" << i+1 << ": " << ADC_TOF2[i] << " -- InA-" << i+1 << ": " << ADC_TOF2[i+24]<< endl;}
    }
    else {
      if ((TDC_TOF2AO[i] > TDC_TOF2_min && TDC_TOF2AO[i] < TDC_TOF2_max) || (TDC_TOF2AI[i] > TDC_TOF2_min && TDC_TOF2AI[i] < TDC_TOF2_max)) {
      cout << "TDC TOF2 OutA-" << i+1 << ": " << TDC_TOF2AO[i] << " -- InA-" << i+1 << ": " << TDC_TOF2AI[i]<< endl;
      }
    }
  }


  for (Int_t i = 0; i<12; i++) {
    if (i == 6) {
      if (ADC_TOF2[55] > 0 || ADC_TOF2[i+36] > 0) {
        if ((TDC_TOF2BO[6] > TDC_TOF2_min && TDC_TOF2BO[6] < TDC_TOF2_max) || (TDC_TOF2BI[6] > TDC_TOF2_min && TDC_TOF2BI[6] < TDC_TOF2_max)) {
        cout << "ADC TOF2 OutB-" << i+1 << ": " << ADC_TOF2[55] << " -- InB-" << i+1 << ": " << ADC_TOF2[i+36]<< " | TDC TOF2 OutB-" << i+1 << ": " << TDC_TOF2BO[i] << " -- InB-" << i+1 << ": " << TDC_TOF2BI[i]<< endl;
        }
        else {cout << "ADC TOF2 OutB-" << i+1 << ": " << ADC_TOF2[55] << " -- InB-" << i+1 << ": " << ADC_TOF2[i+36]<< endl;}
      }
      else {
        if ((TDC_TOF2BO[6] > TDC_TOF2_min && TDC_TOF2BO[6] > TDC_TOF2_max) || (TDC_TOF2BI[6] > TDC_TOF2_min && TDC_TOF2BI[6] > TDC_TOF2_max)) {
        cout << "TDC TOF2 OutB-" << i+1 << ": " << TDC_TOF2BO[i] << " -- InB-" << i+1 << ": " << TDC_TOF2BI[i]<< endl;
        }
      }
    }
    else {
      if (ADC_TOF2[i+12] > 0 || ADC_TOF2[i+36] > 0) {
        if ((TDC_TOF2BO[i] > TDC_TOF2_min && TDC_TOF2BO[i] < TDC_TOF2_max) || (TDC_TOF2BI[i] > TDC_TOF2_min && TDC_TOF2BI[i] < TDC_TOF2_max)) {
        cout << "ADC TOF2 OutB-" << i+1 << ": " << ADC_TOF2[i+12] << " -- InB-" << i+1 << ": " << ADC_TOF2[i+36]<< " | TDC TOF2 OutB-" << i+1 << ": " << TDC_TOF2BO[i] << " -- InB-" << i+1 << ": " << TDC_TOF2BI[i]<< endl;
        }
        else {cout << "ADC TOF2 OutB-" << i+1 << ": " << ADC_TOF2[i+12] << " -- InB-" << i+1 << ": " << ADC_TOF2[i+36]<< endl;}
      }
      else {
        if ((TDC_TOF2BO[i] > TDC_TOF2_min && TDC_TOF2BO[i] < TDC_TOF2_max) || (TDC_TOF2BI[i] > TDC_TOF2_min && TDC_TOF2BI[i] < TDC_TOF2_max)) {
        cout << "TDC TOF2 OutB-" << i+1 << ": " << TDC_TOF2BO[i] << " -- InB-" << i+1 << ": " << TDC_TOF2BI[i]<< endl;
        }
      }
    }
  }

  cout << "" << endl;

      //C2 Counters

      for (int q = 224; q < 238; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 192; q < 206; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 160; q < 174; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 128; q < 142; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      ///

      for (int q = 240; q < 254; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 208; q < 222; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 176; q < 190; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 144; q < 158; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }
    
      for (int q = 96; q < 128; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      cout << "" << endl;

      //C3 Counters

      for (int q = 64; q < 80; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 32; q < 48; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 0; q < 16; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 480; q < 496; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 80; q < 96; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 48; q < 64; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 16; q < 32; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 496; q < 512; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 448; q < 464; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }
      for (int q = 464; q < 480; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      cout << "" << endl;

      //C4 Counters

      for (int q = 416; q < 432; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 384; q < 400; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 352; q < 368; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 320; q < 336; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 288; q < 296; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      ////

      for (int q = 432; q < 448; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 400; q < 416; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 368; q < 384; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 336; q < 352; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 304; q < 312; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      ///

      for (int q = 256; q < 272; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }

      for (int q = 272; q < 288; q++) {
        if (Good_Event && MWPCADC[q]>0) {
          cout << "ADC MWPC Channel " << q << " -- " << par_temp2[q] << ": " << MWPCADC[q] << endl;
        }
      }


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

  for(Int_t ileg=0; ileg<18; ileg++)    tex_Legend_TARGET[ileg]->Draw();

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

  for(Int_t ipal=0; ipal<10; ipal++)  palette_TARGET[ipal]->Draw();
  for(Int_t ileg=0; ileg<10; ileg++)  tex_palette_TARGET[ileg]->Draw();
  tex_palette_TARGET_scale->Draw();
  
  
  if ((ADC_Low_TARGET[max_index] > 0) && (ADC_High_TARGET[max_index] > 0)) {
    marker_TDC_TARGET[max_index]->SetMarkerColor(kViolet+1);
    marker_TDC_TARGET[max_index]->Draw();
  }

  if ((ADC_Low_TARGET[max_index2] > 0) && (ADC_High_TARGET[max_index2] > 0)) {
    marker_TDC_TARGET[max_index2]->SetMarkerColor(1);
    marker_TDC_TARGET[max_index2]->Draw();
  }

  if ((ADC_Low_TARGET[max_index3] > 0) && (ADC_High_TARGET[max_index3] > 0)) {
    marker_TDC_TARGET[max_index3]->SetMarkerColor(1);
    marker_TDC_TARGET[max_index3]->Draw();
  }

  if ((ADC_Low_TARGET[max_index4] > 0) && (ADC_High_TARGET[max_index4] > 0)) {
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

  for(Int_t ileg=0; ileg<18; ileg++)    tex_Legend_TARGET[ileg]->Draw();

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
  Double_t Pi=3.141592653589793238462643383279502884;

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
  
  for(Int_t ileg=0; ileg<18; ileg++)    tex_Legend_TARGET[ileg]->Draw();

  tex_event_TARGET->Draw();

  TLatex *tex_Title_ADC_Low_TARGET = new TLatex(0.01759134,0.9295171,"ADC LG & TDC Cut");
  tex_Title_ADC_Low_TARGET->SetTextSize(0.07);
  tex_Title_ADC_Low_TARGET->SetLineWidth(2);
  tex_Title_ADC_Low_TARGET->Draw();


  TLatex *tex_Subtitle_ADC_Low_TARGET = new TLatex(0.01759134,0.88,ch_ADC_and_TDC_cut);
  tex_Subtitle_ADC_Low_TARGET->SetTextSize(0.04);
  tex_Subtitle_ADC_Low_TARGET->SetLineWidth(2);
  tex_Subtitle_ADC_Low_TARGET->Draw();


  for(Int_t icol=0; icol<256; icol++){
    if (has_TDC_hit[icol]) {
      if(ADC_Low_TARGET[icol]>=0 && ADC_Low_TARGET[icol]<100){ marker_ADCL_TARGET[icol]->SetMarkerColor(kOrange+10); marker_ADCL_TARGET[icol]->Draw();} 
      if(ADC_Low_TARGET[icol]>=100 && ADC_Low_TARGET[icol]<200){ marker_ADCL_TARGET[icol]->SetMarkerColor(kOrange+7); marker_ADCL_TARGET[icol]->Draw();}    
      if(ADC_Low_TARGET[icol]>=200 && ADC_Low_TARGET[icol]<300){ marker_ADCL_TARGET[icol]->SetMarkerColor(kOrange+1); marker_ADCL_TARGET[icol]->Draw();}    
      if(ADC_Low_TARGET[icol]>=300 && ADC_Low_TARGET[icol]<400){ marker_ADCL_TARGET[icol]->SetMarkerColor(kOrange-4); marker_ADCL_TARGET[icol]->Draw();}    
      if(ADC_Low_TARGET[icol]>=400 && ADC_Low_TARGET[icol]<500){ marker_ADCL_TARGET[icol]->SetMarkerColor(kYellow-9); marker_ADCL_TARGET[icol]->Draw();}    
      if(ADC_Low_TARGET[icol]>=500 && ADC_Low_TARGET[icol]<600){ marker_ADCL_TARGET[icol]->SetMarkerColor(kYellow-7); marker_ADCL_TARGET[icol]->Draw();}    
      if(ADC_Low_TARGET[icol]>=600 && ADC_Low_TARGET[icol]<700){ marker_ADCL_TARGET[icol]->SetMarkerColor(kYellow-0); marker_ADCL_TARGET[icol]->Draw();}    
      if(ADC_Low_TARGET[icol]>=700 && ADC_Low_TARGET[icol]<800){ marker_ADCL_TARGET[icol]->SetMarkerColor(kSpring-4); marker_ADCL_TARGET[icol]->Draw();}    
      if(ADC_Low_TARGET[icol]>=800 && ADC_Low_TARGET[icol]<900){ marker_ADCL_TARGET[icol]->SetMarkerColor(kSpring-2); marker_ADCL_TARGET[icol]->Draw();}    
      if(ADC_Low_TARGET[icol]>=900){ marker_ADCL_TARGET[icol]->SetMarkerColor(kGreen-0); marker_ADCL_TARGET[icol]->Draw();}
    } 
  }

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
  }
  

  c2->cd(5);
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

  

  c2->cd(6);

  cout << "  " << endl;


  TF1 *fit_line_ADC4;
  TF1 *fit_line_ADCA;



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


  if (display_fits ==1) {
  h_target_ADC2->SetMarkerStyle(25);
  h_target_ADC2->SetMarkerColor(4);


  h_target_TDC2->SetMarkerStyle(25);
  h_target_TDC2->SetMarkerColor(4);

  cout << "  " << endl;

  if (has_data_TDC2 > 1) {

    phi_TDC=(180/Pi)*atan(par1_TDC);
    sprintf(angle_string_TDC,"Phi = %3.4f deg.",phi_TDC);
    tex_angle_TDC=new TLatex(0,-9,angle_string_TDC);
    tex_angle_TDC->SetTextSize(0.05);
    //tex_angle_TDC->Draw();
  }
  else {
    cout << "Histo Fit 3 Is Empty" << endl;
  }

  h_target_ADC4->SetMarkerStyle(25);
  h_target_ADC4->SetMarkerColor(4);



  }  /// TO CHECK !!

  //// Used for primary track fitting

  c2->cd(6);
  h_target_ADC3->SetMarkerStyle(25);
  h_target_ADC3->SetMarkerColor(4);

  if (has_data_ADC3 > 1) {
    

    for(Int_t j=0; j<256; j++) {
    }


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
        fit_line_ADC4->SetLineColor(2);

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

  cout << "//////  GAP SCORING  //////" << endl;
  for (int k=0; k<12; k++) {
    cout << "Gap score " << k+1 << ": " << gap_counter[k] << endl;
  }

  cout << "" << endl;  
  for (int k=0; k<12; k++) {
  //if (max_gap[k]>=high_gap_hit) cout << "Most probable gap hit: " << k+1 << endl;
  }

  cout << "GAP SELECTED FOR FITTING:  Gap " << gap_to_fit << endl;




    TF1 *fit_line_GoodLG;           TF1 *fit_line_GoodLG_weighed;
    TF1 *fit_line_TDC_selected;     TF1 *fit_line_TDC_selected2;
    TF1 *fit_line_TDC_Gap_Fibers;


    float Gap[12][3][2] = {0};

    float a_fit_TDC_selected;               float b_fit_TDC_selected;
    float a_fit_GoodLG;                     float b_fit_GoodLG;
    float a_fit_GoodLG_weighted;            float b_fit_GoodLG_weighted;
    float a_fit_TDC_selected_weighted;      float b_fit_TDC_selected_weighted;
    float a_fit_TDC_Gap_Fibers;             float b_fit_TDC_Gap_Fibers;


    for(int g=0; g<12; g++)
    {
      Gap[g][0][0] = TOF_Xloc[3*g];
      Gap[g][1][0] = TOF_Xloc[3*g+1];
      Gap[g][2][0] = TOF_Xloc[3*g+2];
      
      Gap[g][0][1] = TOF_Yloc[3*g];
      Gap[g][1][1] = TOF_Yloc[3*g+1];
      Gap[g][2][1] = TOF_Yloc[3*g+2];
    }

   
    for(int l=0; l<256; l++)
    {
      if(adc_low_target[l]>=ADC_cut_TARGET && has_TDC_hit[l])
      {
    //    h_GoodLG->Fill(Xloc[l], Yloc[l]);
        h_GoodLG_copy->Fill(Xloc[l], Yloc[l]);
        h_TDC_selected2->Fill(Xloc[l], Yloc[l]);
    //    h_GoodLG_weighted->Fill(Xloc[l], Yloc[l]);
      }
    }

    for(int i=0; i<12; i++)
    {
      h_Circle->Fill(Gap[i][0][0], Gap[i][0][1]);
      h_Circle->Fill(Gap[i][1][0], Gap[i][1][1]);
      h_Circle->Fill(Gap[i][2][0], Gap[i][2][1]);
  
    }

    h_Target_Center->Fill(0., 0.);

    for(int t=0; t<256; t++)
    {
      h_Target->Fill(Xloc[t],Yloc[t]);
    }

    h_max->Fill(Xloc[max_index], Yloc[max_index]);        h_max_copy->Fill(Xloc[max_index], Yloc[max_index]);
   
    h_kaon->Fill(Xloc[max_index], Yloc[max_index]);       h_kaon_copy->Fill(Xloc[max_index], Yloc[max_index]);
    h_kaon->Fill(Xloc[max_index2], Yloc[max_index2]);     h_kaon_copy->Fill(Xloc[max_index2], Yloc[max_index2]);
    h_kaon->Fill(Xloc[max_index3], Yloc[max_index3]);     h_kaon_copy->Fill(Xloc[max_index3], Yloc[max_index3]);
    h_kaon->Fill(Xloc[max_index4], Yloc[max_index4]);     h_kaon_copy->Fill(Xloc[max_index4], Yloc[max_index4]);


    for(int g=0; g<3; g++){
      h_TOF1->Fill(Gap[gap_to_fit-1][g][0], Gap[gap_to_fit-1][g][1]);
      h_TDC_selected->Fill(Gap[gap_to_fit-1][g][0], Gap[gap_to_fit-1][g][1]); 
      h_TDC_selected2->Fill(Gap[gap_to_fit-1][g][0], Gap[gap_to_fit-1][g][1]);     
    }



    h_Circle->SetMarkerStyle(5);
    h_Circle->SetMarkerSize(1.2);
    h_Circle->SetLineWidth(2);

    h_Target->SetMarkerStyle(25);         h_Target_Center->SetMarkerStyle(5);
    h_Target->SetMarkerColor(1);          h_Target_Center->SetMarkerColor(1);
    h_Target->SetMarkerSize(0.8);         h_Target_Center->SetMarkerSize(1);

    h_GoodLG_copy->SetMarkerStyle(21);
    h_GoodLG_copy->SetMarkerColor(2);
    h_GoodLG_copy->SetMarkerSize(0.8);

    h_target_TDC->SetMarkerStyle(25);     h_target_TDC_copy->SetMarkerStyle(25);      h_TDC_selected->SetMarkerStyle(25);   h_TDC_selected2->SetMarkerStyle(25);
    h_target_TDC->SetMarkerColor(2);      h_target_TDC_copy->SetMarkerColor(2);       h_TDC_selected->SetMarkerColor(2);;   h_TDC_selected2->SetMarkerColor(2);
    h_target_TDC->SetMarkerSize(0.8);     h_target_TDC_copy->SetMarkerSize(0.8);      h_TDC_selected->SetMarkerSize(0.8);   h_TDC_selected2->SetMarkerSize(0.8);

    h_kaon->SetMarkerStyle(21);           h_kaon_copy->SetMarkerStyle(21);            h_TDC_Gap_Fibers->SetMarkerStyle(25);
    h_kaon->SetMarkerColor(1);            h_kaon_copy->SetMarkerColor(1);             h_TDC_Gap_Fibers->SetMarkerColor(4);
    h_kaon->SetMarkerSize(0.8);           h_kaon_copy->SetMarkerSize(0.8);            h_TDC_Gap_Fibers->SetMarkerSize(0.8);

    h_max->SetMarkerStyle(21);            h_max_copy->SetMarkerStyle(21);
    h_max->SetMarkerColor(kViolet+1);     h_max_copy->SetMarkerColor(kViolet+1);
    h_max->SetMarkerSize(1);              h_max_copy->SetMarkerSize(1);

    h_TOF1->SetMarkerStyle(20);           h_TOF1_closest->SetMarkerStyle(20);
    h_TOF1->SetMarkerColor(2);            h_TOF1_closest->SetMarkerColor(3);
    h_TOF1->SetMarkerSize(1.5);           h_TOF1_closest->SetMarkerSize(1.5);
  
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
    VerticalAxis->Draw("same");
    //h_Target_Center->Draw("same");
    h_GoodLG_copy->Draw("same");
    //h_kaon_copy->Draw("same");
    //h_max_copy->Draw("same");
    //h_TOF1->Draw("same");


     

    //h_TDC_Gap_Fibers->Draw();
    //h_TDC_Gap_Fibers->Fit("pol1", "QCM");
    //fit_line_TDC_Gap_Fibers = h_TDC_Gap_Fibers->GetFunction("pol1");
    //fit_line_TDC_Gap_Fibers->SetLineWidth(2);
    //fit_line_TDC_Gap_Fibers->SetLineColor(1);
    //a_fit_TDC_Gap_Fibers=fit_line_TDC_Gap_Fibers->GetParameter(1);      
    //b_fit_TDC_Gap_Fibers=fit_line_TDC_Gap_Fibers->GetParameter(0);


 

    float x_int_TDC[2];                           float y_int_TDC[2];
    float x_int_TDC_selected[2];                  float y_int_TDC_selected[2];
    float x_int_TDC_selected_weighted[2];         float y_int_TDC_selected_weighted[2];
    float x_int_TDC_Gap_Fibers[2];                float y_int_TDC_Gap_Fibers[2];
    float x_int_TDC_SFT[2];                       float y_int_TDC_SFT[2];
    float x_int_GoodLG_SFT[2];                    float y_int_GoodLG_SFT[2];


    x_int_TDC_selected[0] = intersectx1(a_fit_TDC_selected, b_fit_TDC_selected, R_TOF1);
    x_int_TDC_selected[1] = intersectx2(a_fit_TDC_selected, b_fit_TDC_selected, R_TOF1);
    y_int_TDC_selected[0] = y1_int(x_int_TDC_selected[0], a_fit_TDC_selected, b_fit_TDC_selected);
    y_int_TDC_selected[1] = y2_int(x_int_TDC_selected[1], a_fit_TDC_selected, b_fit_TDC_selected);
  



    //x_int_TDC_Gap_Fibers[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);
    //x_int_TDC_Gap_Fibers[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);   
    //y_int_TDC_Gap_Fibers[0] = y1_int(x_int_TDC_Gap_Fibers[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers); 
    //y_int_TDC_Gap_Fibers[1] = y2_int(x_int_TDC_Gap_Fibers[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);

    /// Selection of the Good Intersect Coordinates
    //////////////////////////////////////////////////////////////////////
    
    float x_TDC_sel_intersect1=0.;   float y_TDC_sel_intersect1=0.;
    float x_TDC_Gap_Fibers=0.;       float y_TDC_Gap_Fibers=0.;
  
    float dist1_TDC_selected[2];
    float dist1_GoodLG[2];
    float dist1_TDC_Gap_Fibers[2];

    for(int i=0; i<2; i++)
    {
      dist1_TDC_selected[i] = distance(x_int_TDC_selected[i], y_int_TDC_selected[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
    //  dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
    }    
    
    if(dist1_TDC_selected[0] < dist1_TDC_selected[1])
    {
      x_TDC_sel_intersect1 = x_int_TDC_selected[0];
      y_TDC_sel_intersect1 = y_int_TDC_selected[0];
    }
    else if(dist1_TDC_selected[1] < dist1_TDC_selected[0])
    {
      x_TDC_sel_intersect1 = x_int_TDC_selected[1];
      y_TDC_sel_intersect1 = y_int_TDC_selected[1];
    }
    else cout << "ERROR !" << endl;


    
    //if(dist1_TDC_Gap_Fibers[0] < dist1_TDC_Gap_Fibers[1])
    //{
    //  x_TDC_Gap_Fibers = x_int_TDC_Gap_Fibers[0];
    //  y_TDC_Gap_Fibers = y_int_TDC_Gap_Fibers[0];
    //}
    //else if(dist1_TDC_Gap_Fibers[1] < dist1_TDC_Gap_Fibers[0])
    //{
    //  x_TDC_Gap_Fibers = x_int_TDC_Gap_Fibers[1];
    //  y_TDC_Gap_Fibers = y_int_TDC_Gap_Fibers[1];
    //}
    //else cout << "ERROR !" << endl;


    /// Selection of the Good TOF1 Section
    //////////////////////////////////////////////////////////////////////
    float dist2_TDC_selected[3];

    float dist2_TDC_selected_min = 1000.;
    int selected_TDC_selected = 0;

    for(int ii=0; ii<3; ii++)
    {
      dist2_TDC_selected[ii] = distance(x_TDC_sel_intersect1, y_TDC_sel_intersect1, Gap[gap_to_fit-1][ii][0], Gap[gap_to_fit-1][ii][1]);

      if(dist2_TDC_selected[ii] <= dist2_TDC_selected_min)
      {
        dist2_TDC_selected_min = dist2_TDC_selected[ii];
        selected_TDC_selected = ii;
      }


    }

    h_TOF1_closest->Fill(Gap[gap_to_fit-1][selected_TDC_selected][0], Gap[gap_to_fit-1][selected_TDC_selected][1]);

    //////////////////////////////////////////////////////////////////////




    h_int_TDC->Fill(x_int_TDC[0],y_int_TDC[0]);
    h_int_TDC->Fill(x_int_TDC[1],y_int_TDC[1]);


    h_int_TDC_selected->Fill(x_TDC_sel_intersect1, y_TDC_sel_intersect1);
    //h_int_TDC_Gap_Fibers->Fill(x_TDC_Gap_Fibers, y_TDC_Gap_Fibers);
 

    /// Add Weight on TOF1
    //////////////////////////////////////////////////////////////////
    
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected][0], Gap[gap_to_fit-1][selected_TDC_selected][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected][0], Gap[gap_to_fit-1][selected_TDC_selected][1]);
   
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);

    /*
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    h_TDC_selected2->Fill(Gap[gap_to_fit-1][selected_TDC_selected+1][0], Gap[gap_to_fit-1][selected_TDC_selected+1][1]);
    */






    //////////////////////////////////////////////////////////////////


    h_int_TDC->SetMarkerStyle(20);            
    h_int_TDC->SetMarkerColor(4);             
    h_int_TDC->SetMarkerSize(0.8);            

    h_int_TDC_selected->SetMarkerStyle(20);   h_int_TDC_selected_weighted->SetMarkerStyle(20);
    h_int_TDC_selected->SetMarkerColor(4);    h_int_TDC_selected_weighted->SetMarkerColor(4);
    h_int_TDC_selected->SetMarkerSize(0.8);   h_int_TDC_selected_weighted->SetMarkerSize(0.8);

    h_int_TDC_Gap_Fibers->SetMarkerStyle(20);
    h_int_TDC_Gap_Fibers->SetMarkerColor(4);
    h_int_TDC_Gap_Fibers->SetMarkerSize(0.8);

    h_int_TDC_SFT->SetMarkerStyle(20);        
    h_int_TDC_SFT->SetMarkerColor(1);         
    h_int_TDC_SFT->SetMarkerSize(0.8);        


    h_TDC_selected2->Draw();
    h_TDC_selected2->Fit("pol1", "QCM");
    fit_line_TDC_selected2 = h_TDC_selected2->GetFunction("pol1");
    fit_line_TDC_selected2->SetLineWidth(2);
    fit_line_TDC_selected2->SetLineColor(2);
    a_fit_TDC_selected_weighted=fit_line_TDC_selected2->GetParameter(1);      
    b_fit_TDC_selected_weighted=fit_line_TDC_selected2->GetParameter(0);



    x_int_TDC_selected_weighted[0] = intersectx1(a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted, R_TOF1);
    x_int_TDC_selected_weighted[1] = intersectx2(a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted, R_TOF1); 
    x_int_TDC_SFT[0] = intersectx1(a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted, R_SFT_L1);
    x_int_TDC_SFT[1] = intersectx2(a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted, R_SFT_L1);

    y_int_TDC_selected_weighted[0] = y1_int(x_int_TDC_selected_weighted[0], a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted);
    y_int_TDC_selected_weighted[1] = y2_int(x_int_TDC_selected_weighted[1], a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted);
    y_int_TDC_SFT[0] = y1_int(x_int_TDC_SFT[0], a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted);
    y_int_TDC_SFT[1] = y2_int(x_int_TDC_SFT[1], a_fit_TDC_selected_weighted, b_fit_TDC_selected_weighted);


    float x_TDC_selected_weighted_intersect1=0.;   float y_TDC_selected_weighted_intersect1=0.;
    float x_TDC_SFT_intersect1=0.;                 float y_TDC_SFT_intersect1=0.;
  
    float dist1_TDC_selected_weighted[2];
    float dist1_TDC_SFT[2];

    for(int i=0; i<2; i++)
    {
      dist1_TDC_selected_weighted[i] = distance(x_int_TDC_selected_weighted[i], y_int_TDC_selected_weighted[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TDC_SFT[i] = distance(x_int_TDC_SFT[i], y_int_TDC_SFT[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
   }    
    
    if(dist1_TDC_selected_weighted[0] < dist1_TDC_selected_weighted[1])
    {
      x_TDC_selected_weighted_intersect1 = x_int_TDC_selected_weighted[0];
      y_TDC_selected_weighted_intersect1 = y_int_TDC_selected_weighted[0];
    }
    else if(dist1_TDC_selected_weighted[1] < dist1_TDC_selected_weighted[0])
    {
      x_TDC_selected_weighted_intersect1 = x_int_TDC_selected_weighted[1];
      y_TDC_selected_weighted_intersect1 = y_int_TDC_selected_weighted[1];
    }
    else cout << "ERROR !" << endl;




    if(dist1_TDC_SFT[0] < dist1_TDC_SFT[1])
    {
      x_TDC_SFT_intersect1 = x_int_TDC_SFT[0];
      y_TDC_SFT_intersect1 = y_int_TDC_SFT[0];
    }
    else if(dist1_TDC_SFT[1] < dist1_TDC_SFT[0])
    {
      x_TDC_SFT_intersect1 = x_int_TDC_SFT[1];
      y_TDC_SFT_intersect1 = y_int_TDC_SFT[1];
    }
    else cout << "ERROR !" << endl;



    h_int_TDC_selected_weighted->Fill(x_TDC_selected_weighted_intersect1, y_TDC_selected_weighted_intersect1);
    h_int_TDC_SFT->Fill(x_TDC_SFT_intersect1, y_TDC_SFT_intersect1);

    TLine *CenterLine_TDC = new TLine(0., 0., x_TDC_selected_weighted_intersect1, y_TDC_selected_weighted_intersect1);


    CenterLine_TDC->SetLineWidth(2);      
    CenterLine_TDC->SetLineColor(3);    
    VerticalAxis->SetLineWidth(2);      VerticalAxis->SetLineColor(1);

    /// Angle Calculation
    ///////////////////////////////////////////////////////////////////////////////

    float a_final_TDC = 0.;
    float alpha_TDC = 0.;
    float tanalpha_TDC = 0.;
    float angle_final_TDC = 0.;

    a_final_TDC = y_TDC_selected_weighted_intersect1 / x_TDC_selected_weighted_intersect1;

    tanalpha_TDC = -a_final_TDC;
    alpha_TDC = atan(tanalpha_TDC);

    if(x_TDC_selected_weighted_intersect1 >= 0) angle_final_TDC = 90. + (alpha_TDC * (180./Pi));
    else angle_final_TDC = 180. + (90. + (alpha_TDC * (180./Pi)));  


    
    char Angle_TDC_string[30];      //char Angle_GoodLD_string[30];

    sprintf(Angle_TDC_string,"#phi = %3.2f", angle_final_TDC);


    TLatex *tex_Angle_TDC;

    tex_Angle_TDC = new TLatex(-45.,43.,Angle_TDC_string);
    tex_Angle_TDC->SetTextSize(0.05);
    tex_Angle_TDC->SetLineWidth(2);

  
    cout << "" << endl;
    cout << "" << endl;
    cout << "Final Angle (deg.) = " << angle_final_TDC << "  (TDC_weighted)" << endl;
    cout << "" << endl;
    cout << "Intersect with SFT Layer 1:  " << "x = " << x_TDC_SFT_intersect1 << "   " << "y = " << y_TDC_SFT_intersect1 << endl;
    cout << "" << endl;


    TF1 *f1=new TF1("f1", "x", -50, 50);
    TGaxis *A1 = new TGaxis(-50,50,50,50,"f1",510,"-U");
    TGaxis *A2 = new TGaxis(50,-50,50,50,"f1",510,"+U");



    h_TDC_selected2->Draw();

    A1->Draw();
    A2->Draw();
    h_GoodLG_copy->Draw("same");
    h_Circle->Draw("same");
    h_kaon_copy->Draw("same");
    h_max_copy->Draw("same");
   
    Gap1l->Draw("same");
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
    h_TOF1->Draw("same");
    h_TOF1_closest->Draw("same");
    VerticalAxis->Draw("same");

    CenterLine_TDC->Draw("same");
    h_int_TDC_selected_weighted->Draw("same");
    h_int_TDC_SFT->Draw("same");
    tex_Angle_TDC->Draw("same");
 


  if (display_fits == 1)
  {  

    TCanvas *cc1;
    cc1 = new TCanvas("Fitting Steps","Fitting Steps",50,50,700,700);
    cc1->Divide(2,2);
    cc1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");


    cc1->cd(1);
    h_target_TDC_copy->Draw();
    A1->Draw();
    A2->Draw();
    h_Circle->Draw("same");
    h_Target->Draw("same");
    h_target_TDC_copy->Draw("same");
    h_GoodLG_copy->Draw("same");
    h_kaon_copy->Draw("same");
    h_max_copy->Draw("same");
    Gap1l->Draw("same");
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
    ell_L1->Draw("same");
    ell_Target->Draw("same");
    h_TOF1->Draw("same");


    cc1->cd(2);
    h_TDC_selected->Draw();
    A1->Draw();
    A2->Draw();
    h_Circle->Draw("same");
    h_GoodLG_copy->Draw("same");
    h_kaon_copy->Draw("same");
    h_max_copy->Draw("same");
    Gap1l->Draw("same");
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
    h_TOF1->Draw("same");
    h_TOF1_closest->Draw("same");
    h_int_TDC_selected->Draw("same");

    
    cc1->cd(3);
    h_TDC_Gap_Fibers->Draw();
    A1->Draw();
    A2->Draw();
    h_Circle->Draw("same");
    //h_kaon_copy->Draw("same");
    //h_max_copy->Draw("same");
    //h_GoodLG_copy->Draw("same");
    Gap1l->Draw("same");
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
    h_TOF1->Draw("same");
    h_TOF1_closest->Draw("same");
    h_int_TDC_Gap_Fibers->Draw("same");


    cc1->cd(4);
    h_TDC_selected2->Draw();
    A1->Draw();
    A2->Draw();
    h_Circle->Draw("same");
    Gap1l->Draw("same");
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
    h_GoodLG_copy->Draw("same");
    h_kaon_copy->Draw("same");
    h_max_copy->Draw("same");
    h_TOF1->Draw("same");
    h_TOF1_closest->Draw("same");
    h_int_TDC_selected_weighted->Draw("same");
    h_int_TDC_SFT->Draw("same");
    CenterLine_TDC->Draw("same");
    VerticalAxis->Draw("same");
    tex_Angle_TDC->Draw("same");
  }


  return;
} // End void













