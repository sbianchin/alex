#ifndef __CINT__
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
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


void Event_Batch_CR(Int_t Run_Number=5, Int_t first_event=0, Int_t last_event=10)
{ 
  int Switch=1; // Displays hit with no HG, but LG (0 = OFF ; 1 = ON)
  int Rotate=1; // When TOF1 is 12 or 6, rotate by -90 deg to fit a horizontal line (0 = OFF ; 1 = ON)


  int T_limit = 3;
  
  
  char path_input[200];                   
  sprintf(path_input,"%s",path_merged); 
  
  char Name_finput[200];
  sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);
  
  

  Int_t adc_high_target[256];       Int_t ADC_High_TARGET[256];    
  Int_t adc_low_target[256];        Int_t ADC_Low_TARGET[256];  
  Int_t tdc_le_target[256][16];     Int_t TDC_LE_TARGET[256];     
  Int_t tdc_te_target[256][16];    

  Int_t ADC_tof1[24];               Int_t ADC_TOF1[24];   


  Int_t ADC_tof1U[12];              Int_t ADC_TOF1U[12];
  Int_t ADC_tof1D[12];              Int_t ADC_TOF1D[12];

  Int_t TDC_tof1U[12];              Int_t TDC_TOF1U[12];
  Int_t TDC_tof1D[12];              Int_t TDC_TOF1D[12];

  Int_t tdc_Ck[14][16];
  Int_t tdc_Cpi[14][16];


  Int_t MwpcADC[512];               Int_t MWPCADC[512];
  
  Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128]; 
  Int_t tdc_le_sft[128][16];        Int_t TDC_LE_SFT[128];  
  Int_t adc_low_sft[128];           Int_t ADC_Low_SFT[128];  


  TChain *fChain= new TChain("Tree");   
  fChain->Add(Name_finput);   
  fChain->SetMakeClass(1);              

  fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);    
  fChain->SetBranchAddress("ADC_Low_TARGET",adc_low_target);      
  fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);        
  fChain->SetBranchAddress("TDC_TE_TARGET",tdc_te_target);       

  fChain->SetBranchAddress("ADC_TOF1U",ADC_tof1U);
  fChain->SetBranchAddress("ADC_TOF1D",ADC_tof1D);
  fChain->SetBranchAddress("TDC_TOF1U",TDC_tof1U);
  fChain->SetBranchAddress("TDC_TOF1D",TDC_tof1D); 

  fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
  fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
  fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);

  fChain->SetBranchAddress("MWPCADC",MwpcADC);
   
    

  char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!
  
  Int_t HG_TARGET_ADC_Thr[256] = {0};
  Int_t LG_TARGET_ADC_Thr[256] = {0};

  Int_t ADC_Low_TARGET_ped[256];    
  Int_t ADC_High_TARGET_ped[256];
  Int_t LG_SFT_ADC_Thr[128] = {0};
  Int_t HG_SFT_ADC_Thr[128] = {0};
  
  for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;
  for(int i=0; i<256; i++)  HG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Thr_HG[i]) + TARGET_ADC_Thr_HG_Offset;
  for(int i=0; i<256; i++)  LG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Thr_LG[i]) + TARGET_ADC_Thr_LG_Offset;
  for(int i=0; i<128; i++)  LG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_LG[i]) + SFT_ADC_Thr_LG_Offset;

  Int_t TDC_min_TARGET = TARGET_TDC_min[0];
  Int_t TDC_max_TARGET = TARGET_TDC_max[0];
  Int_t ADC_TARGET_Thr = HG_TARGET_ADC_Thr[0];

  float R_TARGET = 29.0;
  float R_TOF1 = 47.1;
  float R_SFT_L1 = 40.0;


  /// Good Event Variables
  ////////////////////////
  bool Good_Event=false;


  /// Good Target Event Variables
  ///////////////////////////////
  bool Good_TARGET_Event = false;
  int count_TARGET_evts = 0; 


  /// Good TOF Variables
  //////////////////////
  bool Good_TOF_Event = false;


  /// Good MWPC Events Variables
  //////////////////////////////
  int count_C2X = 0;    int count_C2Y = 0;
  bool Good_MWPC_Event = false;


  /// Gap Hit Variables
  /////////////////////
  int high_gap_hit = 0;
  int gap_to_fit = 0;


  /// Target Rotation Variables
  /////////////////////////////
  int ADC_High_TARGET_temp[256];
  int ADC_Low_TARGET_temp[256];
  int adc_low_target_temp[256];
  bool has_TDC_hit_temp[256];
  bool k_stop_bar_temp[256];
  int TARGET_High_has_neighbours_temp[256]; 


  /// Gap Location Variables
  //////////////////////////
  float Gap[12][3][2] = {{{0}}};

  for(int g=0; g<12; g++){
    Gap[g][0][0] = TOF_Xloc[3*g];
    Gap[g][1][0] = TOF_Xloc[3*g+1];
    Gap[g][2][0] = TOF_Xloc[3*g+2];

    Gap[g][0][1] = TOF_Yloc[3*g];
    Gap[g][1][1] = TOF_Yloc[3*g+1];
    Gap[g][2][1] = TOF_Yloc[3*g+2];
  }  


  /// Intersection Coordinate Variables
  /////////////////////////////////////
  float x_TDC_Gap_Fibers=0.;       float y_TDC_Gap_Fibers=0.;
  float dist1_TDC_Gap_Fibers[2];   
  float x_int_TDC_Gap_Fibers_SFT[2];  float y_int_TDC_Gap_Fibers_SFT[2];

  /// Selection of good TOF1 Variables
  ////////////////////////////////////
  float dist2_TDC_selected[3];
  float dist2_TDC_selected_min = 1000.;
  int selected_TDC_selected = 0; 


  vector<int> tof1_ties;
  int adc_tof1u_min = -10000;
  int adc_tof1d_min = -10000;


  /// Target Intersect Variables
  //////////////////////////////
  float x_TDC_Gap_Fibers_intersect1=0.;          float y_TDC_Gap_Fibers_intersect1=0.;
  float x_TARGET_intersect=0;                    float y_TARGET_intersect=0;

  float dist1_TARGET_intersect[2];             float dist1_TDC_Gap_Fibers_SFT[2];


  /// Centroid Variables
  //////////////////////
  double X_weights = 0.0;
  double Y_weights = 0.0;
  double total_energy = 0.0;
  double X_BAR;
  double Y_BAR;


  double closest_to_centroid = 1000;
  int closest_to_centroid_index = -1;



  /// Intersect Variables
  ///////////////////////
  float x_TDC_Gap_Fibers_SFT_intersect1=0.;      float y_TDC_Gap_Fibers_SFT_intersect1=0.;



  /// Angle Calculation Variables
  ///////////////////////////////
  float alpha_guide = 0.;
  float angle_final_guide = 0.;


  /// MWPC clustering variables
  bool first_cluster = true;
  int cluster_spacing = 0;
  int C2X_clusters = 0;
  int C2Y_clusters = 0;
  int C3X_clusters = 0;
  int C3Y_clusters = 0;
  int C4X_clusters = 0;
  int C4Y_clusters = 0;    

  char MWPC_L_R = ' ';  


  /// Other Variables
  ///////////////////
  const int n_hit = 2;

  const Int_t Angle_ADC_cut = 0;
  int count = 0; 

  int Event_Type = 0;

  bool target_near_tof1 = false;



  /// MWPC Mapping histogram variables
  ////////////////////////////////////
  int C2_histo_display = 1;

    TH1F *h_C2X_L = new TH1F("C2X_L","C2X_L",56,0,56);
    TH1F *h_C2X_R = new TH1F("C2X_R","C2X_R",56,0,56);
    TH1F *h_C2Y_L = new TH1F("C2Y_L","C2Y_L",16,0,16);
    TH1F *h_C2Y_R = new TH1F("C2Y_R","C2Y_R",16,0,16);

   



  // Run through events  
  for(int ievt=first_event; ievt <= last_event; ievt++){

    if(ievt%10000 == 0) cout << "On event: " << ievt << endl; 
    Event_Type = 0;


    bool TARGET_High_has_neighbours[256] = {false}; // Array of High gain target hits which have neighbouring targets hit   
    
   /////////////////////////////////////////////////////////////////////////////////////////////////////

    Good_Event=false;
    bool Good_tof1[12] = {false};

    for(int ivt=ievt; ivt<ievt+1; ivt++){
      fChain->GetEntry(ivt);  


      for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
        ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-HG_SFT_ADC_Thr[j_SFT];
        TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
        ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-LG_SFT_ADC_Thr[j_SFT];
       
      }

      for (int j_TARGET=0; j_TARGET<256; j_TARGET++){
        ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-HG_TARGET_ADC_Thr[j_TARGET];
        ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-LG_TARGET_ADC_Thr[j_TARGET];
        TDC_LE_TARGET[j_TARGET]=tdc_le_target[j_TARGET][0];

        ADC_High_TARGET_ped[j_TARGET]=adc_high_target[j_TARGET]-HG_TARGET_ADC_Thr[j_TARGET] + TARGET_ADC_Thr_HG_Offset;
        ADC_Low_TARGET_ped[j_TARGET]=adc_low_target[j_TARGET]-LG_TARGET_ADC_Thr[j_TARGET] + TARGET_ADC_Thr_LG_Offset; 
      } 

      for(int i=0; i<12; i++){
        ADC_TOF1[i] = ADC_tof1U[i]-TOF1U_ADC_Thr[i];
        ADC_TOF1[i+12] = ADC_tof1D[i]-TOF1D_ADC_Thr[i];
      }


      for (Int_t j_MWPC=0; j_MWPC<512; j_MWPC++) {
        MWPCADC[j_MWPC] = MwpcADC[j_MWPC]-MWPC_ADC_Thr[j_MWPC];
      }

      for (Int_t j_TDCTOF=0; j_TDCTOF<12; j_TDCTOF++) {
        TDC_TOF1U[j_TDCTOF] = TDC_tof1U[j_TDCTOF];
        TDC_TOF1D[j_TDCTOF] = TDC_tof1D[j_TDCTOF];
      }  


      //********* GOOD TARGET EVENTS
      Good_TARGET_Event = false;
      count_TARGET_evts = 0; 
      for(int i=0; i<256; i++){
        if((ADC_High_TARGET[i]>=0 && tdc_le_target[i][0]>=TARGET_TDC_min[i] && tdc_le_target[i][0]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]>=0 && tdc_le_target[i][1]>=TARGET_TDC_min[i] && tdc_le_target[i][1]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]>=0 && tdc_le_target[i][2]>=TARGET_TDC_min[i] && tdc_le_target[i][2]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]>=0 && tdc_le_target[i][3]>=TARGET_TDC_min[i] && tdc_le_target[i][3]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]<0  && ADC_Low_TARGET[i]>0))
        {  
          count_TARGET_evts++;
        }
      }

     if(count_TARGET_evts >= 10) Good_TARGET_Event = true;


      //********* GOOD TOF EVENTS
      bool Good_TOF1_ADC[12]={false};
      bool Good_TOF1_TDC[12]={false};
      bool Good_TOF1[12]={false}; 
      Good_TOF_Event = false;

  
    
      for(int i = 0; i < 12; i++){
        ADC_TOF1U[i] = ADC_TOF1[i];
        ADC_TOF1D[i] = ADC_TOF1[i+12];
      } 
    
    

      for(int i=0; i<12; i++){
        if(ADC_TOF1U[i]>=0 || ADC_TOF1D[i]>=0)  Good_TOF1_ADC[i] = true;

        if((TDC_TOF1U[i]>=TOF1U_TDC_min[i] && TDC_TOF1U[i]<=TOF1U_TDC_max[i]) ||
         (TDC_TOF1D[i]>=TOF1D_TDC_min[i] && TDC_TOF1D[i]<=TOF1D_TDC_max[i]))  Good_TOF1_TDC[i] = true;

        if(Good_TOF1_TDC[i] || Good_TOF1_ADC[i]) Good_TOF1[i] = true;
      }


      for(int i=0; i<12; i++){
        Good_tof1[i] = Good_TOF1[i];
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


  bool has_ADC_TOF1_hit[12] = {false};
  bool has_TDC_TOF1_hit[12] = {false};


    ///Set TOF1 Lines

  for(int i = 0; i<12; i++){
    if (ADC_TOF1U[i]>=500 || ADC_TOF1D[i]>=500) {has_ADC_TOF1_hit[i] = true;}
    if ((TDC_TOF1U[i]>=TOF1U_TDC_min[i] && TDC_TOF1U[i]<=TOF1U_TDC_max[i]) || (TDC_TOF1D[i]>=TOF1D_TDC_min[i] && TDC_TOF1D[i]<=TOF1D_TDC_max[i])) {has_TDC_TOF1_hit[i] = true;}
  }   



  int gap_hit[12] = {0};
  int gap_counter[12] = {0};


      //// GAP SCORING
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


  high_gap_hit = 0;
  gap_to_fit = 0;

  tof1_ties.clear();
  adc_tof1u_min = -10000;
  adc_tof1d_min = -10000;

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

  Double_t ADC_High_SFT_corr[128];    
  Int_t ADC_High_corr_max=0;


  for(int j=0 ; j<128 ; j++){ 
    if(ADC_High_SFT[j]<0)      ADC_High_SFT_corr[j]=0; 
    if(ADC_High_SFT[j]>=0)     ADC_High_SFT_corr[j]=ADC_High_SFT[j]; 
  }

  for(Int_t ii=0; ii<128; ii++){
    if(ADC_High_SFT_corr[ii]>ADC_High_corr_max){
      ADC_High_corr_max=ADC_High_SFT_corr[ii];
    }
  }



  /// SFT
  Int_t has_TDC_SFT_hit[128] = {0};

  for(Int_t ii=0; ii<128; ii++){
    for (Int_t qq=0; qq<6; qq++){
      if (tdc_le_sft[ii][qq] > SFT_TDC_min[ii] && tdc_le_sft[ii][qq] < SFT_TDC_max[ii]) has_TDC_SFT_hit[ii]++;
    }
  }  

  vector<int> layer_1_p_fibers;
  vector<int> layer_2_p_fibers;
  vector<int> layer_3_p_fibers;
  vector<int> layer_4_p_fibers;

  bool L1_sft_flag = false;
  bool L2_sft_flag = false;
  bool L3_sft_flag = false;
  bool L4_sft_flag = false;

  int L1_sft_count = 0;
  int L2_sft_count = 0;
  int L3_sft_count = 0;
  int L4_sft_count = 0;

  int L1_sft_double_count = 0; 
  int L2_sft_double_count = 0;  
  int L3_sft_double_count = 0;   
  int L4_sft_double_count = 0;  


  for(int i = 0; i <128; i++){
    if(ADC_High_SFT_corr[i] != 0 && has_TDC_SFT_hit[i]){
      if(SFT_channel_to_layer[i] == 1){

        L1_sft_flag = false;

        for(vector<int>::iterator it = layer_1_p_fibers.begin(); it != layer_1_p_fibers.end(); it++){
          if(*it == SFT_channel_to_fiber[i]){
            L1_sft_flag = true;
            L1_sft_double_count++;
            break;
          }
        }

        if(!L1_sft_flag){  
          layer_1_p_fibers.push_back(SFT_channel_to_fiber[i]);
          L1_sft_count++;
        }
      }
      else if(SFT_channel_to_layer[i] == 2){

      L2_sft_flag = false;

        for(vector<int>::iterator it = layer_2_p_fibers.begin(); it != layer_2_p_fibers.end(); it++){
          if(*it == SFT_channel_to_fiber[i]){
            L2_sft_flag = true;
            L2_sft_double_count++;
            break;
          }
        }

        if(!L2_sft_flag){  
          layer_2_p_fibers.push_back(SFT_channel_to_fiber[i]);
          L2_sft_count++;
        }
      }
      else if(SFT_channel_to_layer[i] == 3){

        L3_sft_flag = false;

        for(vector<int>::iterator it = layer_3_p_fibers.begin(); it != layer_3_p_fibers.end(); it++){
          if(*it == SFT_channel_to_fiber[i]){
            L3_sft_flag = true;
            L3_sft_double_count++;
            break;
          }
        }

        if(!L3_sft_flag){  
          layer_3_p_fibers.push_back(SFT_channel_to_fiber[i]);
          L3_sft_count++;
        }
      }
      else if(SFT_channel_to_layer[i] == 4){

        L4_sft_flag = false;

        for(vector<int>::iterator it = layer_4_p_fibers.begin(); it != layer_4_p_fibers.end(); it++){
          if(*it == SFT_channel_to_fiber[i]){
            L4_sft_flag = true;
            L4_sft_double_count++;
            break;
          }
        }

        if(!L4_sft_flag){  
          layer_4_p_fibers.push_back(SFT_channel_to_fiber[i]);
          L4_sft_count++;
        }
      }
    }
  }   




  // MWPC clustering 

  // Arrays for finding mwpc clustering
  bool C2X_L_temp[56] = {false};
  bool C2X_R_temp[56] = {false};

  bool C2X_L[56] = {false};
  bool C2X_R[56] = {false};
  bool C2Y_L[16] = {false};
  bool C2Y_R[16] = {false};

  int C2Y_L_hits = 0;
  int C2Y_R_hits = 0;
  int C2X_L_hits = 0;
  int C2X_R_hits = 0;


  //C2 Counters
  for (int q = 224; q < 238; q++) {
    if (MWPCADC[q]>0) {
       C2X_L_temp[q-224] = true;
    }
  }

  for (int q = 192; q < 206; q++) {
    if (MWPCADC[q]>0) {
      C2X_L_temp[q-178] = true;
    }
  }

  for (int q = 160; q < 174; q++) {
    if (MWPCADC[q]>0) {
      C2X_L_temp[q-132] = true;
    }
  }

  for (int q = 128; q < 142; q++) {
    if (MWPCADC[q]>0) {
      C2X_L_temp[q-86] = true;
    }
  }

  for (int q = 240; q < 254; q++) {
    if (MWPCADC[q]>0) {
      C2X_R_temp[q-240] = true;
    }
  }

  for (int q = 208; q < 222; q++) {
    if (MWPCADC[q]>0) {
      C2X_R_temp[q-194] = true;
    }
  }

  for (int q = 176; q < 190; q++) {
    if (MWPCADC[q]>0) {
      C2X_R_temp[q-148] = true;
    }
  }

  for (int q = 144; q < 158; q++) {
    if (MWPCADC[q]>0) {
      C2X_R_temp[q-102] = true;
    }
  }

  // Fill C2Y bool array

  for(int q = 96; q < 112; q++){
    if (MWPCADC[q]>0) {
      C2Y_L[q-96] = true;
      C2Y_L_hits++;
    }
  }

  for(int q = 112; q < 128; q++){
    if (MWPCADC[q]>0) {
      C2Y_R[q-112] = true;
      C2Y_R_hits++;
    }
  }    


  for(int i = 0; i < 56; i++){
    C2X_L[i] = C2X_L_temp[C2_strip[i]];
    C2X_R[i] = C2X_R_temp[C2_strip[i]];
  }

  for(int i = 0; i < 56; i++){
    if(C2X_L[i]>0){
      C2X_L_hits++;
    }
    if(C2X_R[i]>0){
      C2X_R_hits++;
    }    
  }  




  if(C2_histo_display == 1){
    for(int i = 0; i<56; i++){
      if(C2X_L[i])   
        h_C2X_L->Fill(i+0.5);
      if(C2X_R[i])   
        h_C2X_R->Fill(i+0.5);
    }  

    for(int i = 0; i<16; i++){
      if(C2Y_L[i])
        h_C2Y_L->Fill(i+0.5);
      if(C2Y_R[i])
        h_C2Y_R->Fill(i+0.5);  
    }   
  }



  // Find clustering of MWPCs
  ////////////////////////////=
  int C2X_L_clusters = 0;
  int C2Y_L_clusters = 0;
  int C2X_R_clusters = 0;
  int C2Y_R_clusters = 0;



  // count C2X clusters

  first_cluster = true;
  cluster_spacing = 0;

  for(int i = 0; i<56; i++){
    if(C2X_L[i] && first_cluster){
      cluster_spacing = MWPC_cluster_separation + 1;
      first_cluster = false;
    }

    if(C2X_L[i]){
      if(cluster_spacing > MWPC_cluster_separation){
        C2X_L_clusters++;
      }
      cluster_spacing = 0;
    }
    else{
      cluster_spacing++;
    }
  }

  // count C2Y clusters
  first_cluster = true;
  cluster_spacing = 0;

  for(int i = 0; i<16; i++){
    if(C2Y_L[i] && first_cluster){
      cluster_spacing = MWPC_cluster_separation + 1;
      first_cluster = false;
    }

    if(C2Y_L[i]){
      if(cluster_spacing > MWPC_cluster_separation){
        C2Y_L_clusters++;
      }
      cluster_spacing = 0;
    }
    else{
      cluster_spacing++;
    }
  }    


  first_cluster = true;
  cluster_spacing = 0;

  // count C2X clusters
  for(int i = 0; i<56; i++){
    if(C2X_R[i] && first_cluster){
      cluster_spacing = MWPC_cluster_separation + 1;
      first_cluster = false;
    }

    if(C2X_R[i]){
      if(cluster_spacing > MWPC_cluster_separation){
        C2X_R_clusters++;
      }
      cluster_spacing = 0;
    }
    else{
      cluster_spacing++;
    }
  }


  
  // count C2Y clusters
  first_cluster = true;
  cluster_spacing = 0;

  for(int i = 0; i<16; i++){
    if(C2Y_R[i] && first_cluster){
      cluster_spacing = MWPC_cluster_separation + 1;
      first_cluster = false;
    }

    if(C2Y_R[i]){
      if(cluster_spacing > MWPC_cluster_separation){
        C2Y_R_clusters++;
      }
      cluster_spacing = 0;
    }
    else{
      cluster_spacing++;
    }
  }    

      
// TOF1 selected. TOF1 goes from 1 to 12.      
int selected_tof1_top = 12;
int selected_tof1_bottom = 6;



if(has_ADC_TOF1_hit[selected_tof1_top-1] && has_ADC_TOF1_hit[selected_tof1_bottom-1] && 
   has_TDC_TOF1_hit[selected_tof1_top-1] && has_TDC_TOF1_hit[selected_tof1_bottom-1] &&
   Good_TARGET_Event &&
   //C2Y_R_clusters >0 && C2Y_L_clusters >0 &&
   //C2Y_L_hits > 1 && C2Y_R_hits > 1 && C2X_L_hits > 1 && C2X_R_hits > 1 &&
   C2X_R_clusters ==1 && C2Y_R_clusters ==1 && C2X_L_clusters ==1 && C2Y_L_clusters ==1 &&
   //C2X_R_clusters >0 && C2Y_R_clusters >0 && C2X_L_clusters >0 && C2Y_L_clusters >0 &&  
   L1_sft_count > 0 && L2_sft_count > 0 && L3_sft_count > 0 && L4_sft_count > 0 && 
   L1_sft_count < 4 && L2_sft_count < 4 && L3_sft_count < 4 && L4_sft_count < 4)
     cout << ievt << endl;


  }//End for loop for entries;



  // TCanvas *C2_Target_canvas;
  // C2_Target_canvas = new TCanvas("C2 with Target","C2 with Target",50,50,1000,500);
  // C2_Target_canvas->Divide(2,2);
  // C2_Target_canvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

  // C2_Target_canvas->cd(1);
  // h_C2X_L->Draw();

  // C2_Target_canvas->cd(2);
  // h_C2X_R->Draw();

  // C2_Target_canvas->cd(3);
  // h_C2Y_L->Draw();

  // C2_Target_canvas->cd(4);
  // h_C2Y_R->Draw();  


  delete fChain;

  return;
} // End void


