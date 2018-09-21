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
#include "SFT_functions_CR.h"


void SFT_Test_CR_Event(int Run_Number, int evt, double phi_top, double phi_bottom,
  double C2_top, double C2_bottom, double SFT_Array[2],
  double ADC_High_SFT_corr[128], int has_TDC_SFT_hit[128]);
void SFT_angle_selector(int Run_Number, int ievt, double phi, int gap_to_fit_left, int gap_to_fit_right, 
  double C2_intersect_SFT_top, double C2_intersect_SFT_bottom, double SFT_Array[2], double ADC_High_SFT_corr[128], int has_TDC_SFT_hit[128]);

void C2_SFT_Histo(Int_t Run_Number=4558, Int_t num_events=0)
{ 
  int Switch=1; // Displays hit with no HG, but LG (0 = OFF ; 1 = ON)
  int Rotate=1; // When TOF1 is 12 or 6, rotate by -90 deg to fit a horizontal line (0 = OFF ; 1 = ON)


  int T_limit = 3;

  gStyle->SetOptStat(1111111);
  
  
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

  Int_t adc_c2x_r[56];              Int_t ADC_C2X_R[56];
  Int_t adc_c2x_l[56];              Int_t ADC_C2X_L[56];
  Int_t adc_c2y_r[16];              Int_t ADC_C2Y_R[16];
  Int_t adc_c2y_l[16];              Int_t ADC_C2Y_L[16];

  
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

  fChain->SetBranchAddress("ADC_C2X_R",adc_c2x_r);
  fChain->SetBranchAddress("ADC_C2X_L",adc_c2x_l);
  fChain->SetBranchAddress("ADC_C2Y_R",adc_c2y_r);
  fChain->SetBranchAddress("ADC_C2Y_L",adc_c2y_l);

  double left_edge = -100.0;
  double right_edge = 100.0;
  int num_bins = 40;


  char C2_Target_Theta_histo_name[200];
  sprintf(C2_Target_Theta_histo_name,"C2 - Target #theta | Run %d", Run_Number);


  TH2F *h_TDC_Gap_Fibers = new TH2F("Target Histo", "Target Histo", 500, -50, 50, 500, -50, 50);

  TH1F *h_TOF1Position = new TH1F("TOF1Position","TOF1 Position",100,-250,250);

  TH1F *h_gaps_10_4 = new TH1F("Top = 10, Bottom = 4","Top = 10, Bottom = 4", num_bins,left_edge,right_edge);
  TH1F *h_gaps_11_5 = new TH1F("Top = 11, Bottom = 5","Top = 11, Bottom = 5", num_bins,left_edge,right_edge);
  TH1F *h_gaps_12_6 = new TH1F("Top = 12, Bottom = 6","Top = 12, Bottom = 6", num_bins,left_edge,right_edge);
  TH1F *h_gaps_1_7 = new TH1F("Top = 1, Bottom = 7","Top = 1, Bottom = 7", num_bins,left_edge,right_edge);
  TH1F *h_gaps_2_8 = new TH1F("Top = 2, Bottom = 8","Top = 2, Bottom = 8", num_bins,left_edge,right_edge);
  TH1F *h_gaps_all = new TH1F("All ","ALl", num_bins,left_edge,right_edge);

  TH1F *h_delta_theta = new TH1F(C2_Target_Theta_histo_name,C2_Target_Theta_histo_name,61,-30,30);



  char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!

  Int_t nentries = (Int_t)fChain->GetEntries(); 
  
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


  /// INPUT FILE FOR ANGLE
  ////////////////////////
  ////////////////////////

  if(num_events == 0){
    num_events = nentries;
  }
  cout << num_events << endl;
  char event_angles_title[200];
  sprintf(event_angles_title,"Run%d_Event_Angles_%d.txt", Run_Number,num_events);

  ifstream event_angles;
  event_angles.open(event_angles_title);

  vector<double> all_event_angles;
  string Current_Event;
  float current_event;

  if(event_angles.fail()){
    cout << "Error: Could not read angle file." << endl;
    return;
  }
  else{
    while(getline(event_angles, Current_Event)){
      sscanf(Current_Event.c_str(), "%f", &current_event);  
      all_event_angles.push_back(current_event);
      //cout << Current_Event << endl;
      //cout << current_event << endl;
    }
  }



  event_angles.close();



  ofstream angles;
  angles.open("Run_4560_gap6_gap12");


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

  double a_fit_TDC_Gap_Fibers = 0;
  double b_fit_TDC_Gap_Fibers = 0;


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


  /// TOF - Z variables
  /////////////////////
  Int_t TOF1PosCl[12] = { -15, -39, 15, -22, -83, -63, -62, -100, -63, -64, -94, -91};

  double TOF1_pos_top = 0.0;
  double TOF1_pos_bot = 0.0;



  /// Angle Calculation Variables
  ///////////////////////////////
  float alpha_guide = 0.;
  float angle_final_guide = 0.;


  /// MWPC clustering variables
  bool first_cluster = true;
  int cluster_spacing = 0;
  int C2X_clusters = 0;
  int C2Y_clusters = 0;

  char MWPC_L_R = ' ';  


  /// Other Variables
  ///////////////////
  const int n_hit = 2;

  const Int_t Angle_ADC_cut = 0;
  int count = 0; 

  int Event_Type = 0;

  bool target_near_tof1 = false;


  // Run through events  
  for(int ievt=0; ievt < num_events; ievt++){

    if(ievt%10000  == 0) cout << "On event: " << ievt << endl; 
    Event_Type = 0;


    bool TARGET_High_has_neighbours[256] = {false}; // Array of High gain target hits which have neighbouring targets hit   
    
   /////////////////////////////////////////////////////////////////////////////////////////////////////


    Int_t TOF1Pos[12] = {0};

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


      // for (Int_t j_MWPC=0; j_MWPC<512; j_MWPC++) {
      //   MWPCADC[j_MWPC] = MwpcADC[j_MWPC]-MWPC_ADC_Thr[j_MWPC];
      // }

      for(Int_t j_C2=0; j_C2<56; j_C2++){
        ADC_C2X_R[j_C2] = adc_c2x_r[j_C2]-ADC_C2X_Thr[j_C2];
        ADC_C2X_L[j_C2] = adc_c2x_l[j_C2]-ADC_C2X_Thr[j_C2];
      }

      for(Int_t j_CY=0; j_CY<16; j_CY++){
        ADC_C2Y_R[j_CY] = adc_c2y_r[j_CY]-ADC_C2Y_Thr[j_CY];
        ADC_C2Y_L[j_CY] = adc_c2y_l[j_CY]-ADC_C2Y_Thr[j_CY];
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

        if(Good_TOF1_TDC[i]){
          TOF1Pos[i] = ((TDC_TOF1U[i] - TDC_TOF1D[i])/2) -TOF1PosCl[i];
          TOF1Pos[i] = 1.8*TOF1Pos[i];
          h_TOF1Position->Fill(TOF1Pos[i]);
        }
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


    int gap_to_fit_left = 0;
    int gap_to_fit_right = 0;

    int high_adc_u = 0;
    int high_adc_d = 0;

    for(int i = 0; i<6; i++){
      if(has_TDC_TOF1_hit[i] && has_ADC_TOF1_hit[i]){
        if(ADC_TOF1U[i] > high_adc_u && ADC_TOF1D[i] > high_adc_d){
          gap_to_fit_right = i+1;
          high_adc_u = ADC_TOF1U[i];
          high_adc_d = ADC_TOF1D[i];
        }
      }
    }    

    high_adc_u = 0;
    high_adc_d = 0;

    for(int i = 6; i<12; i++){
      if(has_TDC_TOF1_hit[i] && has_ADC_TOF1_hit[i]){
        if(ADC_TOF1U[i] > high_adc_u && ADC_TOF1D[i] > high_adc_d){
          gap_to_fit_left = i+1;
          high_adc_u = ADC_TOF1U[i];
          high_adc_d = ADC_TOF1D[i];
        }
      }
    }  

    // cout << "gap_to_fit_left tof = " << TOF1Pos[gap_to_fit_left-1] << endl;
    // cout << "gap_to_fit_right tof = " << TOF1Pos[gap_to_fit_right-1] << endl;


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


    for(int i=128; i<254; i++){
      if(i!=142 && i!=143 && i!=158 && i!=159 && i!=174 && i!=175 && i!=190 && i!=191 && 
       i!=206 && i!=207 && i!=223 && i!=224 && i!=238 && i!=239){

        if(MWPCADC[i]>=0)  count_C2X++; 
    } 
  } 

  for(int ii=96; ii<128; ii++){

    if(MWPCADC[ii]>=0)  count_C2Y++;
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
  // bool C2X_L_temp[56] = {false};
  // bool C2X_R_temp[56] = {false};

  int C2X_L[56] = {0};
  int C2X_R[56] = {0};
  int C2Y_L[16] = {0};
  int C2Y_R[16] = {0};

  int C2Y_L_hits = 0;
  int C2Y_R_hits = 0;
  int C2X_L_hits = 0;
  int C2X_R_hits = 0;


  for(int i=0; i<56; i++){
    if (ADC_C2X_L[i]>=0){
      C2X_L[i] = ADC_C2X_L[i];
      C2X_L_hits++;
    }
    if (ADC_C2X_R[i]>=0){
      C2X_R[i] = ADC_C2X_R[i];
      C2X_R_hits++;
    }
  }

  for(int i=0; i<16; i++){
    if (ADC_C2Y_L[i]>=0) {
      C2Y_L[i] = ADC_C2Y_L[i];
      C2Y_L_hits++;
    }
    if (ADC_C2Y_R[i]>=0) {
      C2Y_R[i] = ADC_C2Y_R[i];
      C2Y_R_hits++;
    }
  }  



  // Find clustering of MWPCs
  ////////////////////////////

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


  bool good_gap_to_fit = false;
  if((gap_to_fit_left == 10 && gap_to_fit_right == 4)||
    (gap_to_fit_left == 11 && gap_to_fit_right == 5) ||
    (gap_to_fit_left == 12 && gap_to_fit_right == 6) ||
    (gap_to_fit_left == 7 && gap_to_fit_right == 1)  ||
    (gap_to_fit_left == 8 && gap_to_fit_right == 2))
      good_gap_to_fit = true;

  // ADD CODE TO READ ANGLES INTO angle_final_guide

  //if(all_event_angles[ievt])
  angle_final_guide = all_event_angles[ievt];
  a_fit_TDC_Gap_Fibers = tan(angle_final_guide*M_PI/180.0);



  double target_slope = 0;
  target_slope = a_fit_TDC_Gap_Fibers;


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


    double C2_Z_Top = 0;
    double C2_Z_Bottom = 0;
    for(int i = 0; i < c2xl_num_lines; i++){
      for(int j = 0; j < c2xr_num_lines; j++){

        double x1_C2_coord = 0.0;
        double x2_C2_coord = 0.0;
        double y1_C2_coord = 0.0;
        double y2_C2_coord = 0.0;

        double x1_C2_coord_global = 0.0;
        double x2_C2_coord_global = 0.0;
        double y1_C2_coord_global = 0.0;
        double y2_C2_coord_global = 0.0;


        if(gap_to_fit_right==1 || gap_to_fit_right==2){
          x1_C2_coord = C2X_R_centroids[j];  
          x2_C2_coord = C2X_L_centroids[i];   
        }
        else{
          x1_C2_coord = C2X_L_centroids[i]; 
          x2_C2_coord = C2X_R_centroids[j];  
        }

        y1_C2_coord = 629.4; 
        y2_C2_coord = 629.4;


        double C2_intersect_SFT_top = 0;
        double C2_intersect_SFT_bottom = 0;


        C2_intersect_SFT_top = C2_intersect(points_to_slope(x1_C2_coord,y1_C2_coord,x2_C2_coord,y2_C2_coord),points_to_y_int(x1_C2_coord,y1_C2_coord,x2_C2_coord,y2_C2_coord),R_SFT_L1);
        C2_intersect_SFT_bottom = C2_intersect(points_to_slope(x1_C2_coord,y1_C2_coord,x2_C2_coord,y2_C2_coord),points_to_y_int(x1_C2_coord,y1_C2_coord,x2_C2_coord,y2_C2_coord),-R_SFT_L1);

        C2_Z_Top = C2_intersect_SFT_top;
        C2_Z_Bottom = C2_intersect_SFT_bottom;

        //cout << "C2_Z_Top = " << C2_Z_Top << endl;

        // if(C2_intersect_SFT_top > -135 && C2_intersect_SFT_top < 135)
        //   cout << "C2 Z intersect = " << C2_intersect_SFT_top << "  --  Top "<< endl;

        // if(C2_intersect_SFT_bottom > -135 && C2_intersect_SFT_bottom < 135)
        //   cout << "C2 Z intersect = " << C2_intersect_SFT_bottom << "  --  Bottom "<< endl;

        if(gap_to_fit_left == 10 && gap_to_fit_right == 4){
          x1_C2_coord_global = rotate_by_angle(x1_C2_coord,y1_C2_coord,60,'x');
          x2_C2_coord_global = rotate_by_angle(x2_C2_coord,y2_C2_coord,-120,'x');
          y1_C2_coord_global = rotate_by_angle(x1_C2_coord,y1_C2_coord,60,'y');
          y2_C2_coord_global = rotate_by_angle(x2_C2_coord,y2_C2_coord,-120,'y');
        }
        else if(gap_to_fit_left == 11 && gap_to_fit_right == 5){
          x1_C2_coord_global = rotate_by_angle(x1_C2_coord,y1_C2_coord,30,'x');
          x2_C2_coord_global = rotate_by_angle(x2_C2_coord,y2_C2_coord,-150,'x');
          y1_C2_coord_global = rotate_by_angle(x1_C2_coord,y1_C2_coord,30,'y');
          y2_C2_coord_global = rotate_by_angle(x2_C2_coord,y2_C2_coord,-150,'y');
        }
        else if(gap_to_fit_left == 12 && gap_to_fit_right == 6){
          x1_C2_coord_global = rotate_by_angle(x1_C2_coord,y1_C2_coord,0,'x');
          x2_C2_coord_global = rotate_by_angle(x2_C2_coord,y2_C2_coord,180,'x');
          y1_C2_coord_global = rotate_by_angle(x1_C2_coord,y1_C2_coord,0,'y');
          y2_C2_coord_global = rotate_by_angle(x2_C2_coord,y2_C2_coord,180,'y');
        }
        else if(gap_to_fit_left == 7 && gap_to_fit_right == 1){
          x1_C2_coord_global = rotate_by_angle(x1_C2_coord,y1_C2_coord,-30,'x');
          x2_C2_coord_global = rotate_by_angle(x2_C2_coord,y2_C2_coord,150,'x');
          y1_C2_coord_global = rotate_by_angle(x1_C2_coord,y1_C2_coord,-30,'y');
          y2_C2_coord_global = rotate_by_angle(x2_C2_coord,y2_C2_coord,150,'y');
        }
        else if(gap_to_fit_left == 8 && gap_to_fit_right == 2){
          x1_C2_coord_global = rotate_by_angle(x1_C2_coord,y1_C2_coord,-60,'x');
          x2_C2_coord_global = rotate_by_angle(x2_C2_coord,y2_C2_coord,120,'x');
          y1_C2_coord_global = rotate_by_angle(x1_C2_coord,y1_C2_coord,-60,'y');
          y2_C2_coord_global = rotate_by_angle(x2_C2_coord,y2_C2_coord,120,'y');
        }



        if(Good_TARGET_Event && angle_final_guide >= 0 && good_gap_to_fit && 
          C2X_R_clusters ==1 && C2Y_R_clusters ==1 && C2X_L_clusters ==1 && C2Y_L_clusters ==1){
          double C2_angle = atan(points_to_slope(x1_C2_coord_global,y1_C2_coord_global,x2_C2_coord_global,y2_C2_coord_global))*180.0/M_PI;
          if(C2_angle < 0)
            C2_angle += 360;

          if(angle_final_guide < 180 && C2_angle > 180){
            C2_angle -= 180;
          }
          else if(angle_final_guide > 180 && C2_angle < 180){
            C2_angle += 180;
          }

          h_delta_theta->Fill(C2_angle - angle_final_guide);

          // if(fabs(C2_angle - angle_final_guide) >= 8)
          //   cout << "Event: " << ievt << ", d(theta) = " << C2_angle - angle_final_guide << endl;


          //cout << C2_angle << endl;
        } 
      }
    }


    double SFT_Z_top = 0;
    double SFT_Z_bottom = 0;

    if(gap_to_fit_right == 1 || gap_to_fit_right == 2){
      TOF1_pos_top = TOF1Pos[gap_to_fit_right - 1];
      TOF1_pos_bot = TOF1Pos[gap_to_fit_left - 1];
    }
    else{
      TOF1_pos_top = TOF1Pos[gap_to_fit_left - 1];
      TOF1_pos_bot = TOF1Pos[gap_to_fit_right - 1];
    }  

    double SFT_Array[2] = {0.0};

    if(Good_TARGET_Event && /*C2_Z_Top != 0 && */ angle_final_guide >= 0 &&
      C2X_R_clusters ==1 && C2Y_R_clusters ==1 && C2X_L_clusters ==1 && C2Y_L_clusters ==1 &&
      L1_sft_count > 0 && L2_sft_count > 0 && L3_sft_count > 0 && L4_sft_count > 0 && 
      L1_sft_count < 4 && L2_sft_count < 4 && L3_sft_count < 4 && L4_sft_count < 4){
      //cout << "PRINT" << endl;
        SFT_angle_selector(Run_Number, ievt, angle_final_guide, gap_to_fit_left, gap_to_fit_right, C2_Z_Top, C2_Z_Bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
      }

    SFT_Z_top = SFT_Array[0];
    SFT_Z_bottom= SFT_Array[1];  

    //cout << C2_Z_Top << endl;

    if((angle_final_guide > 85 && angle_final_guide < 90) ||
       (angle_final_guide > 265 && angle_final_guide < 275)){
      angles << ievt << endl;
    }



    if(Good_TARGET_Event && C2_Z_Top != 0 && angle_final_guide >= 0 && SFT_Z_top > 0.1 && 
      C2X_R_clusters ==1 && C2Y_R_clusters ==1 && C2X_L_clusters ==1 && C2Y_L_clusters ==1 &&
      L1_sft_count > 0 && L2_sft_count > 0 && L3_sft_count > 0 && L4_sft_count > 0 && 
      L1_sft_count < 4 && L2_sft_count < 4 && L3_sft_count < 4 && L4_sft_count < 4){
        // cout << "-   " << C2_Z_Top-SFT_Z_top << endl;
        // cout << "C2  " << C2_Z_Top << endl;
        // cout << "SFT " << SFT_Z_top << endl;

        if(gap_to_fit_left == 10 && gap_to_fit_right == 4){
          h_gaps_10_4->Fill(C2_Z_Top-SFT_Z_top);
          h_gaps_10_4->Fill(C2_Z_Bottom-SFT_Z_bottom);
          good_gap_to_fit = true;
        }
        else if(gap_to_fit_left == 11 && gap_to_fit_right == 5){
          h_gaps_11_5->Fill(C2_Z_Top-SFT_Z_top);
          h_gaps_11_5->Fill(C2_Z_Bottom-SFT_Z_bottom);
          good_gap_to_fit = true;
        }
        else if(gap_to_fit_left == 12 && gap_to_fit_right == 6){
          h_gaps_12_6->Fill(C2_Z_Top-SFT_Z_top);
          h_gaps_12_6->Fill(C2_Z_Bottom-SFT_Z_bottom);
          good_gap_to_fit = true;
        }
        else if(gap_to_fit_left == 7 && gap_to_fit_right == 1){
          h_gaps_1_7->Fill(C2_Z_Top-SFT_Z_top);
          h_gaps_1_7->Fill(C2_Z_Bottom-SFT_Z_bottom);
          good_gap_to_fit = true;
        }
        else if(gap_to_fit_left == 8 && gap_to_fit_right == 2){
          h_gaps_2_8->Fill(C2_Z_Top-SFT_Z_top);
          h_gaps_2_8->Fill(C2_Z_Bottom-SFT_Z_bottom);
          good_gap_to_fit = true;
        }

        if(good_gap_to_fit){
          h_gaps_all->Fill(C2_Z_Top-SFT_Z_top);
          h_gaps_all->Fill(C2_Z_Bottom-SFT_Z_bottom);
        }
      
    }



    h_TDC_Gap_Fibers->Reset();

  }//End for loop for entries;


  angles.close();

  // TCanvas *C2_Target_canvas;
  // C2_Target_canvas = new TCanvas("C2 - SFT","C2 - SFT",50,50,1000,500);
  // C2_Target_canvas->Divide(3,2);
  // C2_Target_canvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

  // C2_Target_canvas->cd(1);
  // h_gaps_10_4->Draw();

  // C2_Target_canvas->cd(2);
  // h_gaps_11_5->Draw();

  // C2_Target_canvas->cd(3);
  // h_gaps_12_6->Draw();

  // C2_Target_canvas->cd(4);
  // h_gaps_1_7->Draw();

  // C2_Target_canvas->cd(5);
  // h_gaps_2_8->Draw();

  // C2_Target_canvas->cd(6);
  // h_TOF1Position->Draw();
  // // h_gaps_all->Draw();


  double theta_histo_a = 0;
  double theta_histo_b = 0;
  double theta_histo_c = 0;

  TF1 *theta_gaus_fit;

  h_delta_theta->Fit("gaus", "QCM");
  theta_gaus_fit = h_delta_theta->GetFunction("gaus");
  theta_histo_a = theta_gaus_fit->GetParameter(0);
  theta_histo_b = theta_gaus_fit->GetParameter(1);
  theta_histo_c = theta_gaus_fit->GetParameter(2);

  // cout << "PARAM A = " << theta_histo_a << endl;
  // cout << "PARAM B = " << theta_histo_b << endl;
  // cout << "PARAM C = " << theta_histo_c << endl;

  char fit_func_string[200];
  sprintf(fit_func_string, "Fit Function = %3.1fexp#left(-0.5#frac{(x-%3.1f)^{2}}{%3.1f^{2}}#right)",theta_histo_a,theta_histo_b,theta_histo_c);

  TLatex *fit_func = new TLatex(10,theta_histo_a*(2.0/3.0),fit_func_string);
  fit_func->SetTextSize(0.03);
  fit_func->SetLineWidth(2);

  TCanvas *Delta_Theta_Canvas;
  Delta_Theta_Canvas = new TCanvas("C2 - Target theta","C2 - Target theta",50,50,1000,500);
  Delta_Theta_Canvas->cd(1);
  h_delta_theta->Draw("same");
  theta_gaus_fit->Draw("same");
  fit_func->Draw("same");


  delete fChain;

  return;
} // End void


void SFT_angle_selector(int Run_Number, int ievt, double phi, int gap_to_fit_left, int gap_to_fit_right, 
  double C2_intersect_SFT_top, double C2_intersect_SFT_bottom, double SFT_Array[2], double ADC_High_SFT_corr[128], int has_TDC_SFT_hit[128]){

  if(gap_to_fit_left==12 && gap_to_fit_right==6){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,phi, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,phi+180,C2_intersect_SFT_top, C2_intersect_SFT_bottom,  SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }
  }
  if(gap_to_fit_left==12 && gap_to_fit_right==5){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,phi, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,phi+180, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }    
  }
  if(gap_to_fit_left==12 && gap_to_fit_right==4){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,phi, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array,ADC_High_SFT_corr, has_TDC_SFT_hit);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,phi+180, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }      
  }
  if(gap_to_fit_left==11 && gap_to_fit_right==6){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,phi, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,phi+180, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }      
  }
  if(gap_to_fit_left==11 && gap_to_fit_right==5){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,phi, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,phi+180, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }        
  }
  if(gap_to_fit_left==11 && gap_to_fit_right==4){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,phi, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,phi+180, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }      
  }
  if(gap_to_fit_left==10 && gap_to_fit_right==6){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,phi, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,phi+180, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }       
  }
  if(gap_to_fit_left==10 && gap_to_fit_right==5){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,phi, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,phi+180, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }      
  }
  if(gap_to_fit_left==10 && gap_to_fit_right==4){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,phi, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
     }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,phi+180, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }      
  }
  if(gap_to_fit_left==8 && gap_to_fit_right==2){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,phi, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,phi+180, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }      
  }
  if(gap_to_fit_left==8 && gap_to_fit_right==1){
    if(phi > 180){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,phi, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,phi+180, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }   
  }
  if(gap_to_fit_left==7 && gap_to_fit_right==2){
    if(phi > 270){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,phi, C2_intersect_SFT_top, C2_intersect_SFT_bottom,  SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,phi+180, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }   
  }
  if(gap_to_fit_left==7 && gap_to_fit_right==1){
    if(phi > 270){
      SFT_Test_CR_Event(Run_Number,ievt,phi-180,phi, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }
    else{
      SFT_Test_CR_Event(Run_Number,ievt,phi,phi+180, C2_intersect_SFT_top, C2_intersect_SFT_bottom, SFT_Array, ADC_High_SFT_corr, has_TDC_SFT_hit);
    }   
  }

  return;
}

void SFT_Test_CR_Event(int Run_Number, int evt, double phi_top, double phi_bottom,
  double C2_top, double C2_bottom, double SFT_Array[2],
  double ADC_High_SFT_corr[128], int has_TDC_SFT_hit[128]){

  SFT_print(ADC_High_SFT_corr, has_TDC_SFT_hit, SFT_channel_to_fiber, 
    evt, phi_top, phi_bottom, C2_top,C2_bottom, SFT_Array, false);

  return;
}
