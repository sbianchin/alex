#ifndef __CINT__
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>
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
#include "ADC_Thresholds.h"
#include "TDC_Windows.h"
#include "Cuts_and_Windows.h"
#include "MWPC_Thr.h"
#endif

#include "intersect.cxx"
#include "C2_Strip_transform.h"


void HK_Batch(Int_t Run_Number=5, Int_t first_event=0, Int_t last_event=10)
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
  Int_t ADC_tof2[48]; 

  Int_t ADC_tof1U[12];              Int_t ADC_TOF1U[12];
  Int_t ADC_tof1D[12];              Int_t ADC_TOF1D[12];

  Int_t TDC_tof1U[12];              Int_t TDC_TOF1U[12];
  Int_t TDC_tof1D[12];              Int_t TDC_TOF1D[12];

  Int_t tdc_Ck[14][16];
  Int_t tdc_Cpi[14][16];


  Int_t ADC_tof2AO[12];             Int_t ADC_TOF2AO[12];
  Int_t ADC_tof2BO[12];             Int_t ADC_TOF2BO[12];
  Int_t ADC_tof2AI[12];             Int_t ADC_TOF2AI[12];
  Int_t ADC_tof2BI[12];             Int_t ADC_TOF2BI[12];

  Int_t TDC_tof2AO[12];             Int_t TDC_TOF2AO[12];
  Int_t TDC_tof2BO[12];             Int_t TDC_TOF2BO[12];
  Int_t TDC_tof2AI[12];             Int_t TDC_TOF2AI[12];
  Int_t TDC_tof2BI[12];             Int_t TDC_TOF2BI[12];

  Int_t MwpcADC[512];               Int_t MWPCADC[512];
  
  double ParError = 999.99;
  double ChiS = 0.0;
  double ndf = 0.0;
  double prob = 0.0;

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
  
  fChain->SetBranchAddress("ADC_TOF2AO",ADC_tof2AO);
  fChain->SetBranchAddress("ADC_TOF2AI",ADC_tof2AI);
  fChain->SetBranchAddress("ADC_TOF2BO",ADC_tof2BO);
  fChain->SetBranchAddress("ADC_TOF2BI",ADC_tof2BI);    
  fChain->SetBranchAddress("TDC_TOF2AO",TDC_tof2AO);
  fChain->SetBranchAddress("TDC_TOF2AI",TDC_tof2AI);
  fChain->SetBranchAddress("TDC_TOF2BO",TDC_tof2BO);
  fChain->SetBranchAddress("TDC_TOF2BI",TDC_tof2BI);   

  fChain->SetBranchAddress("TDC_Ck", tdc_Ck);
  fChain->SetBranchAddress("TDC_Cpi", tdc_Cpi); 

  fChain->SetBranchAddress("MWPCADC",MwpcADC);
  
  TH2F *h_TDC_Gap_Fibers = new TH2F("h_Gap_Fibers_Name", "h_Gap_Fibers_Title", 500, -50, 50, 500, -50, 50);
    
    
  /// OUTPUT FILES
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////  
  char Event_Angle_Title[100];
  sprintf(Event_Angle_Title,"Event_Angle_Run%d__Event%d_to_Event%d.txt",Run_Number, first_event, last_event);

  FILE *Event_Angle;
  Event_Angle = fopen(Event_Angle_Title, "w");

  fprintf(Event_Angle, "Run Number %d, Event %d to Event %d.\n", Run_Number, first_event, last_event);
  fprintf(Event_Angle, "Event Types: \n");
  fprintf(Event_Angle, "            0 - Good Event\n");
  fprintf(Event_Angle, "            1 - No Target\n");
  fprintf(Event_Angle, "            2 - No TOF\n");
  fprintf(Event_Angle, "            3 - No MWPC\n");
  fprintf(Event_Angle, "            4 - No Target, No TOF\n");
  fprintf(Event_Angle, "            5 - No Target, No MWPC\n");
  fprintf(Event_Angle, "            6 - No TOF, No MWPC\n");
  fprintf(Event_Angle, "            7 - No Target, No TOF, No MWPC\n");
  fprintf(Event_Angle, "            8 - No Target fibers to fit\n");
  fprintf(Event_Angle, "            9 - Event on blacklist\n\n");
  fprintf(Event_Angle, "Event     Event Type     Angle     X-K Stop  Y-K Stop  K-Stop Bar    X(40)     Y(40)    X(47.1)   Y(47.1)    N(C_k)    N(C_pi)    N(C2X)    N(C2Y)    N(C3X)    N(C3Y)    N(C4X)    N(C4Y)     MWPC    \n");


  // File for writing csv, no header.
  char Event_Angle_csv_Title[100];
  sprintf(Event_Angle_csv_Title,"Event_Angle_Run%d__Event%d_to_Event%d.csv",Run_Number, first_event, last_event);

  FILE *Event_Angle_csv;
  Event_Angle_csv = fopen(Event_Angle_csv_Title, "w");



  char Bad_Event_Statistics_Title[100];
  sprintf(Bad_Event_Statistics_Title, "Bad_Event_Statistics_Run%d__Event%d_to_Event%d.txt",Run_Number, first_event, last_event);

  ofstream Bad_Event_Statistics;
  Bad_Event_Statistics.open(Bad_Event_Statistics_Title);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  /// END OUTPUT FILES

  char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!
  
  Int_t HG_TARGET_ADC_Thr[256] = {0};
  Int_t LG_TARGET_ADC_Thr[256] = {0};

  Int_t ADC_Low_TARGET_ped[256];    
  Int_t ADC_High_TARGET_ped[256];
  
  
  for(int i=0; i<256; i++)  HG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Thr_HG[i]) + TARGET_ADC_Thr_HG_Offset;
  for(int i=0; i<256; i++)  LG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Thr_LG[i]) + TARGET_ADC_Thr_LG_Offset;


  Int_t TDC_min_TARGET = TARGET_TDC_min[0];
  Int_t TDC_max_TARGET = TARGET_TDC_max[0];
  Int_t ADC_TARGET_Thr = HG_TARGET_ADC_Thr[0];

  float R_TARGET = 29.0;
  float R_TOF1 = 47.1;
  float R_SFT_L1 = 40.0;



  // Event blacklist. Store all black list events in a vector.
  ifstream blacklist;
  blacklist.open("Event_Blacklist_Test.txt"); 

  vector<int> blacklist_events;

  bool Event_On_Blacklist = false;
  string Current_Event;
  int current_event;

  if(blacklist.fail()){
    cout << "Error: Could not read blacklist file." << endl;
  }
  else{
    while(getline(blacklist,Current_Event)){
      sscanf(Current_Event.c_str(), "%d", &current_event);  
      blacklist_events.push_back(current_event);
    }
  }

  blacklist.close();  


  /// Bad Event Statistics variables
  //////////////////////////////////
  int no_mwpc_count = 0;
  int no_tof_count = 0;
  int no_target_count = 0;
  int no_tof_no_mwpc_count = 0;
  int no_target_no_mwpc_count = 0;
  int no_target_no_tof_count = 0;
  int no_target_no_tof_no_mwpc_count = 0;
  int no_target_near_tof1 = 0;
  int blacklist_event_count = 0;


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
  int count_C3X = 0;    int count_C3Y = 0;
  int count_C4X = 0;    int count_C4Y = 0;
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


  /// Best Fit Line Variables
  ///////////////////////////
  float a_fit_TDC_Gap_Fibers=0.;             float b_fit_TDC_Gap_Fibers=0.;


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

  int selected_TOF2 = 0;

  const Int_t Angle_ADC_cut = 0;
  int count = 0; 

  int Event_Type = 0;

  bool target_near_tof1 = false;


  // N(C_k) and N(C_pi) counters
  ///////////////////////////////
  int TDC_Ck_counter = 0;
  int TDC_Cpi_counter = 0;


  // Run through events  
  for(int ievt=first_event; ievt <= last_event; ievt++){

    if(ievt%10000 == 0) cout << "On event: " << ievt << endl; 
    Event_Type = 0;


    bool TARGET_High_has_neighbours[256] = {false}; // Array of High gain target hits which have neighbouring targets hit   
    
   /////////////////////////////////////////////////////////////////////////////////////////////////////

    Good_Event=false;
    bool Good_tof1[12] = {false};
    bool Good_tof2[12] = {false};

    for(int ivt=ievt; ivt<ievt+1; ivt++){
      fChain->GetEntry(ivt);  

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


      for (Int_t j_TOF2=0; j_TOF2<12; j_TOF2++) {
        ADC_TOF2AO[j_TOF2] = ADC_tof2AO[j_TOF2]-TOF2AO_ADC_Thr[j_TOF2];
        ADC_TOF2BO[j_TOF2] = ADC_tof2BO[j_TOF2]-TOF2AI_ADC_Thr[j_TOF2];
        ADC_TOF2AI[j_TOF2] = ADC_tof2AI[j_TOF2]-TOF2BO_ADC_Thr[j_TOF2];
        ADC_TOF2BI[j_TOF2] = ADC_tof2BI[j_TOF2]-TOF2BI_ADC_Thr[j_TOF2];
      }

      for (Int_t j_MWPC=0; j_MWPC<512; j_MWPC++) {
        MWPCADC[j_MWPC] = MwpcADC[j_MWPC]-MWPC_ADC_Thr[j_MWPC];
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

     if(count_TARGET_evts >= n_hit) Good_TARGET_Event = true;


      //********* GOOD TOF EVENTS
      bool Good_TOF1_ADC[12]={false};   bool Good_TOF2_ADC[12]={false};
      bool Good_TOF1_TDC[12]={false};   bool Good_TOF2_TDC[12]={false};
      bool Good_TOF1[12]={false};       bool Good_TOF2[12]={false};
      bool Good_TOFs[12]={false};
      Good_TOF_Event = false;




      // Check if current event is on the blacklist.
      Event_On_Blacklist = false;

      for(vector<int>::iterator it = blacklist_events.begin(); it != blacklist_events.end(); it++){
        if(ievt == *it){
          Event_On_Blacklist = true;
          break;
        }
      }  
    
      for(int i = 0; i < 12; i++){
        ADC_TOF1U[i] = ADC_TOF1[i];
        ADC_TOF1D[i] = ADC_TOF1[i+12];
      } 
    
    

      for(int i=0; i<12; i++){
        if(ADC_TOF1U[i]>=0 || ADC_TOF1D[i]>=0)  Good_TOF1_ADC[i] = true;

        if((TDC_TOF1U[i]>=TOF1U_TDC_min[i] && TDC_TOF1U[i]<=TOF1U_TDC_max[i]) ||
         (TDC_TOF1D[i]>=TOF1D_TDC_min[i] && TDC_TOF1D[i]<=TOF1D_TDC_max[i]))  Good_TOF1_TDC[i] = true;

        if(Good_TOF1_TDC[i] || Good_TOF1_ADC[i]) Good_TOF1[i] = true;


        if((ADC_TOF2AO[i]>=0 && ADC_TOF2AI[i]>=0) || (ADC_TOF2BO[i]>=0 && ADC_TOF2BI[i]>=0)) Good_TOF2_ADC[i] = true;

        if(((TDC_TOF2AO[i]>=TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<=TOF2AO_TDC_max[i])  &&
          (TDC_TOF2AI[i]>=TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<=TOF2AI_TDC_max[i])) ||
         ((TDC_TOF2BO[i]>=TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<=TOF2BO_TDC_max[i])  &&
          (TDC_TOF2BI[i]>=TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<=TOF2BI_TDC_max[i])))  Good_TOF2_TDC[i] = true;

        if(Good_TOF2_ADC[i] || Good_TOF2_TDC[i]) Good_TOF2[i] = true;
      }


      for(int k=0; k<12; k++){
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



      //********* GOOD MWPC EVENTS  
      count_C2X = 0;     count_C2Y = 0;
      count_C3X = 0;     count_C3Y = 0;
      count_C4X = 0;     count_C4Y = 0;
      Good_MWPC_Event = false;

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
      
      
      if(Good_TARGET_Event && Good_TOF_Event && Good_MWPC_Event && !Event_On_Blacklist)  Good_Event = true;

      
      if(Good_TARGET_Event && Good_TOF_Event && !Good_MWPC_Event){
        no_mwpc_count++;
        Event_Type = 3;
        break;
      }
    
      if(Good_TARGET_Event && !Good_TOF_Event && Good_MWPC_Event){
        no_tof_count++;
        Event_Type = 2;
        break;
      }

      if(!Good_TARGET_Event && Good_TOF_Event && Good_MWPC_Event){
        no_target_count++;
        Event_Type = 1;
        break;
      }
      
      if(Good_TARGET_Event && !Good_TOF_Event && !Good_MWPC_Event){
        Event_Type = 6;
        break;
      }
      
      if(!Good_TARGET_Event && Good_TOF_Event && !Good_MWPC_Event){
        no_target_no_mwpc_count++;
        Event_Type = 5;
        break;
      }
      
      if(!Good_TARGET_Event && !Good_TOF_Event && Good_MWPC_Event){
        no_target_no_tof_count++;
        Event_Type = 4;
        break;
      }
      
      if(!Good_TARGET_Event && !Good_TOF_Event && !Good_MWPC_Event){
        no_target_no_tof_no_mwpc_count++;
        Event_Type = 7;
        break;
      }

      if(Event_On_Blacklist){
        blacklist_event_count++;
        Event_Type = 9;
        break;
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
    bool has_ADC_TOF2_hit[12] = {false};
    bool has_TDC_TOF2_hit[12] = {false};


    ///Set TOF1 Lines

for(int i = 0; i<12; i++){
  if (ADC_TOF1U[i]>=0 || ADC_TOF1D[i]>=0) {has_ADC_TOF1_hit[i] = true;}
  if ((TDC_TOF1U[i]>=TOF1U_TDC_min[i] && TDC_TOF1U[i]<=TOF1U_TDC_max[i]) || (TDC_TOF1D[i]>=TOF1D_TDC_min[i] && TDC_TOF1D[i]<=TOF1D_TDC_max[i])) {has_TDC_TOF1_hit[i] = true;}
}   

    ///Set TOF2 Lines
for(int i = 0; i<12;i++){
  if ((ADC_TOF2AO[i]>=0 && ADC_TOF2AI[i]>=0) || (ADC_TOF2BO[i]>=0 && ADC_TOF2BI[i]>=0)) {has_ADC_TOF2_hit[i]=true;}
  if (((TDC_TOF2AO[i]>=TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<=TOF2AO_TDC_max[i]) && (TDC_TOF2AI[i]>=TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<=TOF2AI_TDC_max[i])) 
    || ((TDC_TOF2BO[i]>=TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<=TOF2BO_TDC_max[i]) && (TDC_TOF2BI[i]>=TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<=TOF2BI_TDC_max[i]))) {has_TDC_TOF2_hit[i]=true;}
}   


selected_TOF2 = 0;

      // Determine which TOF2 is hit
for(int i = 0; i<12; i++){
  if(has_TDC_TOF2_hit[i] && has_ADC_TOF2_hit[i])
    selected_TOF2 = i + 1;  
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


  if(selected_TOF2 == 0)
    selected_TOF2 = gap_to_fit;
// MWPC clustering 

    // Arrays for finding mwpc clustering

    bool C2X_L_temp[56] = {false};
    bool C2X_R_temp[56] = {false};

    bool C2X_L[56] = {false};
    bool C2X_R[56] = {false};
    bool C2Y_L[16] = {false};
    bool C2Y_R[16] = {false};
    bool C3X_L[64] = {false};
    bool C3X_R[64] = {false};
    bool C3Y_L[16] = {false};
    bool C3Y_R[16] = {false};
    bool C4X_L[72] = {false};
    bool C4X_R[72] = {false};
    bool C4Y_L[16] = {false};
    bool C4Y_R[16] = {false};


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




  for(int i = 0; i < 56; i++){
    C2X_L[i] = C2X_L_temp[C2_strip[i]];
    C2X_R[i] = C2X_R_temp[C2_strip[i]];
  }



    // Fill C2Y bool array

    for(int q = 96; q < 112; q++){
      if (MWPCADC[q]>0) {
        C2Y_L[q-96] = true;
      }
    }

    for(int q = 112; q < 128; q++){
      if (MWPCADC[q]>0) {
        C2Y_R[q-112] = true;
      }
    }      


    //C3 Counters

    for (int q = 64; q < 80; q++) {
      if (MWPCADC[q]>0) {
        C3X_L[q-64] = true;
      }
    }

    for (int q = 32; q < 48; q++) {
      if (MWPCADC[q]>0) {
        C3X_L[q-16] = true;
      }
    }

    for (int q = 0; q < 16; q++) {
      if (MWPCADC[q]>0) {
        C3X_L[q+32] = true;
      }
    }

    for (int q = 480; q < 496; q++) {
      if (MWPCADC[q]>0) {
        C3X_L[q-432] = true;
      }
    }

    for (int q = 80; q < 96; q++) {
      if (MWPCADC[q]>0) {
        C3X_R[q-80] = true;
      }
    }

    for (int q = 48; q < 64; q++) {
      if (MWPCADC[q]>0) {
        C3X_R[q-32] = true;
      }
    }

    for (int q = 16; q < 32; q++) {
      if (MWPCADC[q]>0) {
        C3X_R[q+16] = true;
      }
    }

    for (int q = 496; q < 512; q++) {
      if (MWPCADC[q]>0) {
        C3X_R[q-448] = true;
      }
    }

    for (int q = 448; q < 464; q++) {
      if (MWPCADC[q]>0) {
        C3Y_L[q-448] = true;
      }
    }
    for (int q = 464; q < 480; q++) {
      if (MWPCADC[q]>0) {
        C3Y_R[q-464] = true;
      }
    }



    //C4 Counters

    for (int q = 416; q < 432; q++) {
      if (MWPCADC[q]>0) {
        C4X_L[q-416] = true;
      }
    }

    for (int q = 384; q < 400; q++) {
      if (MWPCADC[q]>0) {
        C4X_L[q-368] = true;
      }
    }

    for (int q = 352; q < 368; q++) {
      if (MWPCADC[q]>0) {
        C4X_L[q-320] = true;
      }
    }

    for (int q = 320; q < 336; q++) {
      if (MWPCADC[q]>0) {
        C4X_L[q-272] = true;
      }
    }

    for (int q = 288; q < 296; q++) {
      if (MWPCADC[q]>0) {
        C4X_L[q-224] = true;
      }
    }

    for (int q = 432; q < 448; q++) {
      if (MWPCADC[q]>0) {
        C4X_R[q-432] = true;
      }
    }

    for (int q = 400; q < 416; q++) {
      if (MWPCADC[q]>0) {
        C4X_R[q-384] = true;
      }
    }

    for (int q = 368; q < 384; q++) {
      if (MWPCADC[q]>0) {
        C4X_R[q-356] = true;
      }
    }

    for (int q = 336; q < 352; q++) {
      if (MWPCADC[q]>0) {
        C4X_R[q-288] = true;
      }
    }

    for (int q = 304; q < 312; q++) {
      if (MWPCADC[q]>0) {
        C4X_R[q-240] = true;
      }
    }

    for (int q = 256; q < 272; q++) {
      if (MWPCADC[q]>0) {
        C4Y_L[q-256] = true;
      }
    }

    for (int q = 272; q < 288; q++) {
      if (MWPCADC[q]>0) {
        C4Y_R[q-272] = true;
      }
    }


// Find clustering of MWPCs
////////////////////////////=
C2X_clusters = 0;
C2Y_clusters = 0;
C3X_clusters = 0;
C3Y_clusters = 0;
C4X_clusters = 0;
C4Y_clusters = 0;

if(selected_TOF2 > 6){  // LEFT
  MWPC_L_R = 'L';

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
        C2X_clusters++;
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
        C2Y_clusters++;
      }
      cluster_spacing = 0;
    }
    else{
      cluster_spacing++;
    }
  }    

  // count C3X clusters
  first_cluster = true;
  cluster_spacing = 0;

  for(int i = 0; i<64; i++){
    if(C3X_L[i] && first_cluster){
      cluster_spacing = MWPC_cluster_separation + 1;
      first_cluster = false;
    }

    if(C3X_L[i]){
      if(cluster_spacing > MWPC_cluster_separation){
        C3X_clusters++;
      }
      cluster_spacing = 0;
    }
    else{
      cluster_spacing++;
    }
  } 

  // count C3Y clusters
  first_cluster = true;
  cluster_spacing = 0;

  for(int i = 0; i<16; i++){
    if(C3Y_L[i] && first_cluster){
      cluster_spacing = MWPC_cluster_separation + 1;
      first_cluster = false;
    }

    if(C3Y_L[i]){
      if(cluster_spacing > MWPC_cluster_separation){
        C3Y_clusters++;
      }
      cluster_spacing = 0;
    }
    else{
      cluster_spacing++;
    }
  }    

  // count C4X clusters
  first_cluster = true;
  cluster_spacing = 0;

  for(int i = 0; i<72; i++){
    if(C4X_L[i] && first_cluster){
      cluster_spacing = MWPC_cluster_separation + 1;
      first_cluster = false;
    }

    if(C4X_L[i]){
      if(cluster_spacing > MWPC_cluster_separation){
        C4X_clusters++;
      }
      cluster_spacing = 0;
    }
    else{
      cluster_spacing++;
    }
  } 

  // count C4Y clusters
  first_cluster = true;
  cluster_spacing = 0;

  for(int i = 0; i<16; i++){
    if(C4Y_L[i] && first_cluster){
      cluster_spacing = MWPC_cluster_separation + 1;
      first_cluster = false;
    }

    if(C4Y_L[i]){
      if(cluster_spacing > MWPC_cluster_separation){
        C4Y_clusters++;
      }
      cluster_spacing = 0;
    }
    else{
      cluster_spacing++;
    }
  }    
}


else{                // RIGHT
MWPC_L_R = 'R';
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
        C2X_clusters++;
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
        C2Y_clusters++;
      }
      cluster_spacing = 0;
    }
    else{
      cluster_spacing++;
    }
  }    

  // count C3X clusters
  first_cluster = true;
  cluster_spacing = 0;

  for(int i = 0; i<64; i++){
    if(C3X_R[i] && first_cluster){
      cluster_spacing = MWPC_cluster_separation + 1;
      first_cluster = false;
    }

    if(C3X_R[i]){
      if(cluster_spacing > MWPC_cluster_separation){
        C3X_clusters++;
      }
      cluster_spacing = 0;
    }
    else{
      cluster_spacing++;
    }
  } 

  // count C3Y clusters
  first_cluster = true;
  cluster_spacing = 0;

  for(int i = 0; i<16; i++){
    if(C3Y_R[i] && first_cluster){
      cluster_spacing = MWPC_cluster_separation + 1;
      first_cluster = false;
    }

    if(C3Y_R[i]){
      if(cluster_spacing > MWPC_cluster_separation){
        C3Y_clusters++;
      }
      cluster_spacing = 0;
    }
    else{
      cluster_spacing++;
    }
  }    

  // count C4X clusters
  first_cluster = true;
  cluster_spacing = 0;

  for(int i = 0; i<72; i++){
    if(C4X_R[i] && first_cluster){
      cluster_spacing = MWPC_cluster_separation + 1;
      first_cluster = false;
    }

    if(C4X_R[i]){
      if(cluster_spacing > MWPC_cluster_separation){
        C4X_clusters++;
      }
      cluster_spacing = 0;
    }
    else{
      cluster_spacing++;
    }
  } 

  // count C4Y clusters
  first_cluster = true;
  cluster_spacing = 0;

  for(int i = 0; i<16; i++){
    if(C4Y_R[i] && first_cluster){
      cluster_spacing = MWPC_cluster_separation + 1;
      first_cluster = false;
    }

    if(C4Y_R[i]){
      if(cluster_spacing > MWPC_cluster_separation){
        C4Y_clusters++;
      }
      cluster_spacing = 0;
    }
    else{
      cluster_spacing++;
    }
  }
}     




bool k_stop_bar[256] = {false};
X_weights = 0.0;
Y_weights = 0.0;
total_energy = 0.0;

for(int i = 0; i<256; i++){
  if(ADC_High_TARGET_ped[i] > HG_KAON && ADC_Low_TARGET_ped[i] > LG_KAON && has_TDC_hit_Kstop[i]){
    k_stop_bar[i] = true;
  }
}  

for(int i = 0; i<256; i++){
  if(k_stop_bar[i]){
    X_weights += ADC_Low_TARGET_ped[i]*Xloc[i];
    Y_weights += ADC_Low_TARGET_ped[i]*Yloc[i];
    total_energy += ADC_Low_TARGET_ped[i];
  }
}  

X_BAR = X_weights/total_energy;
Y_BAR = Y_weights/total_energy;  



      //Determine closest bar to the centroid, and draw in red
closest_to_centroid = 1000;
closest_to_centroid_index = -1;

for(int i = 0; i<256; i++){
  if(distance(X_BAR,Y_BAR,Xloc[i],Yloc[i]) <= closest_to_centroid){
    closest_to_centroid_index = i;
    closest_to_centroid = distance(X_BAR,Y_BAR,Xloc[i],Yloc[i]);
  }
}


    // Determine if a hit target has any hit neighbours which have high gain
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


      ////////////////////// N(C_k) and N(C_pi) 
TDC_Ck_counter = 0;
TDC_Cpi_counter = 0;

for(int ic = 0; ic < 14; ic++){
  if(((tdc_Ck[ic][0] >= TDC_Ck_min) && (tdc_Ck[ic][0] <= TDC_Ck_max))||
   ((tdc_Ck[ic][1] >= TDC_Ck_min) && (tdc_Ck[ic][1] <= TDC_Ck_max))||
   ((tdc_Ck[ic][2] >= TDC_Ck_min) && (tdc_Ck[ic][2] <= TDC_Ck_max))||
   ((tdc_Ck[ic][3] >= TDC_Ck_min) && (tdc_Ck[ic][3] <= TDC_Ck_max))||
   ((tdc_Ck[ic][4] >= TDC_Ck_min) && (tdc_Ck[ic][4] <= TDC_Ck_max))||
   ((tdc_Ck[ic][5] >= TDC_Ck_min) && (tdc_Ck[ic][5] <= TDC_Ck_max)))
   TDC_Ck_counter++;

 if(((tdc_Cpi[ic][0] >= TDC_Cpi_min) && (tdc_Cpi[ic][0] <= TDC_Cpi_max))||
   ((tdc_Cpi[ic][1] >= TDC_Cpi_min) && (tdc_Cpi[ic][1] <= TDC_Cpi_max))||
   ((tdc_Cpi[ic][2] >= TDC_Cpi_min) && (tdc_Cpi[ic][2] <= TDC_Cpi_max))||
   ((tdc_Cpi[ic][3] >= TDC_Cpi_min) && (tdc_Cpi[ic][3] <= TDC_Cpi_max))||
   ((tdc_Cpi[ic][4] >= TDC_Cpi_min) && (tdc_Cpi[ic][4] <= TDC_Cpi_max))||
   ((tdc_Cpi[ic][5] >= TDC_Cpi_min) && (tdc_Cpi[ic][5] <= TDC_Cpi_max)))
   TDC_Cpi_counter++; 
}


    ////////////////////// First Histogram Fill

count = 0;

    // TARGET ROTATION (90 deg.)

    // ROTATE_CHANGE: any line tagged with ROTATE_CHANGE is used when tof1 is 6 or 12.

const int gap_to_fit_rotate = gap_to_fit;

if((gap_to_fit==12 || gap_to_fit==6 || selected_TOF2==6 || selected_TOF2==12) && Rotate==1){    

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

  for(int i = 0; i<256; i++){
    ADC_High_TARGET_temp[i] = ADC_High_TARGET[i];
    ADC_Low_TARGET_temp[i] = ADC_Low_TARGET[i];
    has_TDC_hit_temp[i] = has_TDC_hit[i];
    adc_low_target_temp[i] = adc_low_target[i];
    TARGET_High_has_neighbours_temp[i] = TARGET_High_has_neighbours[i];
    k_stop_bar_temp[i] = k_stop_bar[i];
  }

  for(int i=0; i<256; i++){
    ADC_High_TARGET[i] = ADC_High_TARGET_temp[TARGET_Rotated_index[i]];
    ADC_Low_TARGET[i] = ADC_Low_TARGET_temp[TARGET_Rotated_index[i]];
    has_TDC_hit[i] = has_TDC_hit_temp[TARGET_Rotated_index[i]];
    adc_low_target[i] = adc_low_target[TARGET_Rotated_index[i]];
    TARGET_High_has_neighbours[i] = TARGET_High_has_neighbours_temp[TARGET_Rotated_index[i]];
    k_stop_bar[i] = k_stop_bar_temp[TARGET_Rotated_index[i]];
  }
}   



for(Int_t i=0; i<256; i++){
  if(TARGET_High_has_neighbours[i]){
    if(ADC_High_TARGET[i]>=Angle_ADC_cut && has_TDC_hit[i]){
      if(!k_stop_bar[i]){
        h_TDC_Gap_Fibers->Fill(Xloc[i],Yloc[i]);
        count++;
      }

          if(IsIn(i,channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],   // Additional weight for fibers close to the edge if hit
            channel[gap_to_fit-1][2], channel[gap_to_fit-1][3], 
            channel[gap_to_fit-1][4], channel[gap_to_fit-1][5], 
            channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){
            if(!k_stop_bar[i]){
              h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
              h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
            }
          }
        }

        if(ADC_High_TARGET[i]<0 && ADC_Low_TARGET[i]>=0 && Switch==1){
          if(!k_stop_bar[i]){
            h_TDC_Gap_Fibers->Fill(Xloc[i],Yloc[i]);
            count++;
          }
        }
      }
    }

    if(count == 0)
      Event_Type = 8;

    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    a_fit_TDC_Gap_Fibers=0.;             b_fit_TDC_Gap_Fibers=0.;


    // Add weight to center point of selected TOF1 (gap_to_fit_rotate)
    for(int i = 0; i<3; i++){
      h_TDC_Gap_Fibers->Fill(Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
    }


    if(h_TDC_Gap_Fibers->GetEntries() != 0){
      // Fit 1: Targets hit and center point of selected TOF1
      h_TDC_Gap_Fibers->Fit("pol1", "QCM0");
      a_fit_TDC_Gap_Fibers=h_TDC_Gap_Fibers->GetFunction("pol1")->GetParameter(1);      
      b_fit_TDC_Gap_Fibers=h_TDC_Gap_Fibers->GetFunction("pol1")->GetParameter(0);    
    }
    else 
          cout << "Empty Histo 1" << endl; // Dont attempt to fit if histo is empty, will crash.

        
    // ReFill histogram with points closer than max_dist
        h_TDC_Gap_Fibers->Reset();
        for(int i = 0; i<256; i++){   
          if(distance_to_line(Xloc[i],Yloc[i],a_fit_TDC_Gap_Fibers,b_fit_TDC_Gap_Fibers) <= max_dist && TARGET_High_has_neighbours[i] && !k_stop_bar[i]){ 

            if(ADC_High_TARGET[i]>=Angle_ADC_cut && has_TDC_hit[i]){
              h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);


              if(IsIn(i,channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],   
                channel[gap_to_fit-1][2], channel[gap_to_fit-1][3], 
                channel[gap_to_fit-1][4], channel[gap_to_fit-1][5], 
                channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){
                h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);     
              h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]); 
            }
          }
          if(ADC_High_TARGET[i]<0 && ADC_Low_TARGET[i]>=0 && Switch==1){
            h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
          }
        } 
      }

    //Fit 2: Targets hit within max_distance of the best fit line from Fit 1
      if(h_TDC_Gap_Fibers->GetEntries() != 0){
        h_TDC_Gap_Fibers->Fit("pol1", "QCM0");
        a_fit_TDC_Gap_Fibers=h_TDC_Gap_Fibers->GetFunction("pol1")->GetParameter(1);      
        b_fit_TDC_Gap_Fibers=h_TDC_Gap_Fibers->GetFunction("pol1")->GetParameter(0);
      }
      


      float x_int_TDC[2];                           float y_int_TDC[2];
      float x_int_TDC_Gap_Fibers[2];                float y_int_TDC_Gap_Fibers[2];
      float x_int_TARGET[2];                        float y_int_TARGET[2];

      x_int_TDC_Gap_Fibers[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);
      x_int_TDC_Gap_Fibers[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);   
      y_int_TDC_Gap_Fibers[0] = y1_int(x_int_TDC_Gap_Fibers[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers); 
      y_int_TDC_Gap_Fibers[1] = y2_int(x_int_TDC_Gap_Fibers[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);
      
      

      /// Selection of the Good Intersect Coordinates
      //////////////////////////////////////////////////////////////////////
      for(int i=0; i<2; i++)
        dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);

      
      if(dist1_TDC_Gap_Fibers[0] < dist1_TDC_Gap_Fibers[1])
      {
        x_TDC_Gap_Fibers = x_int_TDC_Gap_Fibers[0];
        y_TDC_Gap_Fibers = y_int_TDC_Gap_Fibers[0];
      }
      else
      {
        x_TDC_Gap_Fibers = x_int_TDC_Gap_Fibers[1];
        y_TDC_Gap_Fibers = y_int_TDC_Gap_Fibers[1];
      }

      /// Selection of the Good TOF1 Section
      //////////////////////////////////////////////////////////////////////
      dist2_TDC_selected_min = 1000.;

      for(int ii=0; ii<3; ii++)
      {
        dist2_TDC_selected[ii] = distance(x_TDC_Gap_Fibers, y_TDC_Gap_Fibers, Gap[gap_to_fit-1][ii][0], Gap[gap_to_fit-1][ii][1]);
        
        if(dist2_TDC_selected[ii] <= dist2_TDC_selected_min)
        {
          dist2_TDC_selected_min = dist2_TDC_selected[ii];
          selected_TDC_selected = ii;
        }
      }



      /// Add Weight on TOF1
      //////////////////////////////////////////////////////////////////
      for(int i = 0; i<3; i++)
        h_TDC_Gap_Fibers->Fill(Gap[gap_to_fit-1][selected_TDC_selected][0], Gap[gap_to_fit-1][selected_TDC_selected][1]);
      //////////////////////////////////////////////////////////////////

      if(h_TDC_Gap_Fibers->GetEntries() != 0){
        h_TDC_Gap_Fibers->Fit("pol1", "QCM0");
        a_fit_TDC_Gap_Fibers = h_TDC_Gap_Fibers->GetFunction("pol1")->GetParameter(1);      
        b_fit_TDC_Gap_Fibers = h_TDC_Gap_Fibers->GetFunction("pol1")->GetParameter(0);
        ParError= h_TDC_Gap_Fibers->GetFunction("pol1")-> GetParError(1);
        ChiS = h_TDC_Gap_Fibers->GetFunction("pol1")-> GetChisquare();
        ndf = h_TDC_Gap_Fibers->GetFunction("pol1")-> GetNDF();
        prob = h_TDC_Gap_Fibers->GetFunction("pol1")-> GetProb();
      }
      else
        cout << "Empty Histo 3" << endl; // Dont attempt to fit if histo is empty, will crash.
      
      
      x_int_TDC_Gap_Fibers_SFT[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_SFT_L1);
      x_int_TDC_Gap_Fibers_SFT[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_SFT_L1);   
      y_int_TDC_Gap_Fibers_SFT[0] = y1_int(x_int_TDC_Gap_Fibers_SFT[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers); 
      y_int_TDC_Gap_Fibers_SFT[1] = y2_int(x_int_TDC_Gap_Fibers_SFT[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);  

      x_TDC_Gap_Fibers_intersect1=0.;          y_TDC_Gap_Fibers_intersect1=0.;
      x_TDC_Gap_Fibers_SFT_intersect1=0.;      y_TDC_Gap_Fibers_SFT_intersect1=0.;


      x_int_TDC_Gap_Fibers[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);
      x_int_TDC_Gap_Fibers[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);  

      y_int_TDC_Gap_Fibers[0] = y1_int(x_int_TDC_Gap_Fibers[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers); 
      y_int_TDC_Gap_Fibers[1] = y2_int(x_int_TDC_Gap_Fibers[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);

      x_int_TARGET[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TARGET);
      x_int_TARGET[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TARGET);
      y_int_TARGET[0] = y1_int(x_int_TARGET[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);
      y_int_TARGET[1] = y2_int(x_int_TARGET[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);




      for(int i=0; i<2; i++){ 
        dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
        dist1_TARGET_intersect[i] = distance(x_int_TARGET[i], y_int_TARGET[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
        dist1_TDC_Gap_Fibers_SFT[i] = distance(x_int_TDC_Gap_Fibers_SFT[i], y_int_TDC_Gap_Fibers_SFT[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      }

      if(dist1_TDC_Gap_Fibers[0] < dist1_TDC_Gap_Fibers[1]){
        x_TDC_Gap_Fibers_intersect1 = x_int_TDC_Gap_Fibers[0];
        y_TDC_Gap_Fibers_intersect1 = y_int_TDC_Gap_Fibers[0];
      }
      else{
        x_TDC_Gap_Fibers_intersect1 = x_int_TDC_Gap_Fibers[1];
        y_TDC_Gap_Fibers_intersect1 = y_int_TDC_Gap_Fibers[1];
      }

      if(dist1_TARGET_intersect[0] < dist1_TARGET_intersect[1]){
        x_TARGET_intersect = x_int_TARGET[1];
        y_TARGET_intersect = y_int_TARGET[1];
      }
      else{
        x_TARGET_intersect = x_int_TARGET[0];
        y_TARGET_intersect = y_int_TARGET[0];
      }


      if(dist1_TDC_Gap_Fibers_SFT[0] < dist1_TDC_Gap_Fibers_SFT[1]){
        x_TDC_Gap_Fibers_SFT_intersect1 = x_int_TDC_Gap_Fibers_SFT[0];
        y_TDC_Gap_Fibers_SFT_intersect1 = y_int_TDC_Gap_Fibers_SFT[0];
      }
      else if(dist1_TDC_Gap_Fibers_SFT[1] < dist1_TDC_Gap_Fibers_SFT[0]){
        x_TDC_Gap_Fibers_SFT_intersect1 = x_int_TDC_Gap_Fibers_SFT[1];
        y_TDC_Gap_Fibers_SFT_intersect1 = y_int_TDC_Gap_Fibers_SFT[1];
      }


      /// Angle Calculation
      ///////////////////////////////////////////////////////////////////////////////
      alpha_guide = atan((y_TDC_Gap_Fibers_intersect1 - y_TARGET_intersect) / (x_TDC_Gap_Fibers_intersect1 - x_TARGET_intersect));


      /////////////////////### Option 2
      if((x_TDC_Gap_Fibers_intersect1 - x_TARGET_intersect) < 0) angle_final_guide = 180. + alpha_guide * (180./PI);

      if((x_TDC_Gap_Fibers_intersect1 - x_TARGET_intersect) >= 0){
        if((y_TDC_Gap_Fibers_intersect1 - y_TARGET_intersect) >= 0)  angle_final_guide = alpha_guide * (180./PI);
        else angle_final_guide = alpha_guide * (180./PI) + 360.;
      }


      if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6 || selected_TOF2==6 || selected_TOF2==12) && Rotate==1)
        angle_final_guide += 90.;
      if(angle_final_guide > 360.)
        angle_final_guide -= 360.;
      
      // Readable file
      if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6 || selected_TOF2==6 || selected_TOF2==12) && Rotate==1 && count != 0)
        fprintf(Event_Angle,"%5d %10d     %10.2f%10.2f%10.2f%10d%13.2f%10.2f%10.2f%10.2f%10d%10d%10d%10d%10d%10d%10d%10d%10c\n",
          ievt, Event_Type,  angle_final_guide,X_BAR, Y_BAR, closest_to_centroid_index,-y_TDC_Gap_Fibers_SFT_intersect1,
          x_TDC_Gap_Fibers_SFT_intersect1,-y_TDC_Gap_Fibers_intersect1,x_TDC_Gap_Fibers_intersect1, TDC_Ck_counter, TDC_Cpi_counter,
          C2X_clusters, C2Y_clusters, C3X_clusters, C3Y_clusters, C4X_clusters, C4Y_clusters, MWPC_L_R);
      else if(count != 0)
        fprintf(Event_Angle,"%5d %10d     %10.2f%10.2f%10.2f%10d%13.2f%10.2f%10.2f%10.2f%10d%10d%10d%10d%10d%10d%10d%10d%10c\n",
          ievt, Event_Type,  angle_final_guide,X_BAR,Y_BAR, closest_to_centroid_index, x_TDC_Gap_Fibers_SFT_intersect1,y_TDC_Gap_Fibers_SFT_intersect1,
          x_TDC_Gap_Fibers_intersect1,y_TDC_Gap_Fibers_intersect1, TDC_Ck_counter, TDC_Cpi_counter,
          C2X_clusters, C2Y_clusters, C3X_clusters, C3Y_clusters, C4X_clusters, C4Y_clusters, MWPC_L_R);  
      else
        fprintf(Event_Angle,"%5d %10d     %10.2s%10.2f%10.2f%10d%13.2s%10.2s%10.2s%10.2s%10d%10d%10d%10d%10d%10d%10d%10d%10c\n",
          ievt, Event_Type,  "-",X_BAR,Y_BAR, closest_to_centroid_index, "-","-","-","-", TDC_Ck_counter, TDC_Cpi_counter,
          C2X_clusters, C2Y_clusters, C3X_clusters, C3Y_clusters, C4X_clusters, C4Y_clusters, MWPC_L_R);  
      
      // CSV file
      if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6 || selected_TOF2==6 || selected_TOF2==12) && Rotate==1 && count != 0)
        fprintf(Event_Angle_csv,"%5d,%10d,%10.2f,%10.2f,%10.2f,%10d,%13.2f,%10.2f,%10.2f,%10.2f,%10d,%10d,%c,%10d,%10d,%10d,%10d,%10d,%10d\n",
          ievt, Event_Type,  angle_final_guide,X_BAR, Y_BAR, closest_to_centroid_index,-y_TDC_Gap_Fibers_SFT_intersect1,
          x_TDC_Gap_Fibers_SFT_intersect1,-y_TDC_Gap_Fibers_intersect1,x_TDC_Gap_Fibers_intersect1, TDC_Ck_counter, TDC_Cpi_counter, MWPC_L_R,
          C2X_clusters, C2Y_clusters, C3X_clusters, C3Y_clusters, C4X_clusters, C4Y_clusters);
      else if(count != 0)
        fprintf(Event_Angle_csv,"%5d,%10d,%10.2f,%10.2f,%10.2f,%10d,%13.2f,%10.2f,%10.2f,%10.2f,%10d,%10d,%c,%10d,%10d,%10d,%10d,%10d,%10d\n",
          ievt, Event_Type,  angle_final_guide,X_BAR,Y_BAR, closest_to_centroid_index, x_TDC_Gap_Fibers_SFT_intersect1,y_TDC_Gap_Fibers_SFT_intersect1,
          x_TDC_Gap_Fibers_intersect1,y_TDC_Gap_Fibers_intersect1, TDC_Ck_counter, TDC_Cpi_counter, MWPC_L_R,
          C2X_clusters, C2Y_clusters, C3X_clusters, C3Y_clusters, C4X_clusters, C4Y_clusters);  
      else
        fprintf(Event_Angle_csv,"%5d,%10d,%10.2s,%10.2f,%10.2f,%10d,%13.2s,%10.2s,%10.2s,%10.2s,%10d,%10d,%c,%10d,%10d,%10d,%10d,%10d,%10d\n",
          ievt, Event_Type,  "-",X_BAR,Y_BAR, closest_to_centroid_index, "-","-","-","-", TDC_Ck_counter, TDC_Cpi_counter, MWPC_L_R,
          C2X_clusters, C2Y_clusters, C3X_clusters, C3Y_clusters, C4X_clusters, C4Y_clusters);  

      //if(Switch_Display>0){ 
      //cout << fixed;
      //cout << setw(4) << Run_Number << "  ";
      //cout << setw(7) << ievt << "  ";
      //cout << setw(2) << gap_to_fit_rotate << "  ";
      //cout << setw(7) << setprecision(3) << angle_final_guide << "  ";
      ////cout << setw(6) << setprecision(3) << Delta_phi_deg << "  ";
      //cout << setw(8) << setprecision(3) << ChiS/ndf << "  ";
      //cout << setw(7) << setprecision(3) << X_BAR << "  ";
      //cout << setw(7) << setprecision(3) << Y_BAR << "  ";
      //cout << endl;
      ////cout << setw(9) << setprecision(3) << Z_TOF1[gap_to_fit_rotate-1] << "  ";


      //cout << Run_Number << "  " 
      //     << ievt << "  "
      //     << gap_to_fit_rotate << "  "
      //     << angle_final_guide << "  "
      //     << ChiS/ndf << "  "
      //     << X_BAR << "  "
      //     << Y_BAR << "  "
      //     << closest_to_centroid_index << "  "
      //     << endl;  

      h_TDC_Gap_Fibers->Reset();
  }//End for loop for entries;
  
  
  
  Bad_Event_Statistics << "No MWPC count: " <<  no_mwpc_count << endl;
  Bad_Event_Statistics << "No TOF count: " << no_tof_count << endl;
  Bad_Event_Statistics << "No Target count: " << no_target_count << endl;
  Bad_Event_Statistics << "No TOF, No MWPC count: " << no_tof_no_mwpc_count << endl;
  Bad_Event_Statistics << "No Target, No MWPC count: " << no_target_no_mwpc_count << endl;
  Bad_Event_Statistics << "No Target, No TOF count: " << no_target_no_tof_count << endl;
  Bad_Event_Statistics << "No Target, No TOF, No MWPC count: " << no_target_no_tof_no_mwpc_count << endl;
  Bad_Event_Statistics << "Number of events on blacklist: " << blacklist_event_count << endl;

  
  Bad_Event_Statistics.close(); 
  fclose(Event_Angle);  

  delete h_TDC_Gap_Fibers;
  delete fChain;

  return;
} // End void


