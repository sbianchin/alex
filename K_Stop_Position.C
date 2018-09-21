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
#include "ADC_Thresholds.h"
#include "TDC_Windows.h"
#include "Cuts_and_Windows.h"
#include "MWPC_Thr.h"
#endif
 
#include "intersect.cxx"
  

void K_Stop_Position(Int_t Run_Number=5, Int_t first_event=0, Int_t last_event=10)
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
	  Int_t TDC_LE_TARGET_corrected[256][6]={0};
	  char TDC_LE_TARGET_corr[256][6][20];
	
	  Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];  
	  Int_t adc_low_sft[128];           Int_t ADC_Low_SFT[128];   
	  Int_t tdc_le_sft[128][16];        Int_t TDC_LE_SFT[128];         
	  Int_t tdc_te_sft[128][16];    
	
	  Int_t ADC_tof1[24];               Int_t ADC_TOF1[24];   
	  Int_t ADC_tof2[48]; 
	   
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
	  
    
    /// OUTPUT FILES
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////  
  	char Event_K_Stop_Title[100];
  	sprintf(Event_K_Stop_Title,"Event_K_Stop_Position_Run%d__Event%d_to_Event%d.txt",Run_Number, first_event, last_event);
  	
  	FILE *Event_K_Stop_Position;

	Event_K_Stop_Position = fopen(Event_K_Stop_Title, "w");
	 
	fprintf(Event_K_Stop_Position,"Event       Xbar        Ybar         ADC Low Max bar\n");
	
	
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
	/// END OUTPUT FILES
	
	char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!

	
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
	
	
	
	/// Other Variables
	///////////////////
	const int n_hit = 2;
    const Int_t Angle_ADC_cut = 0;
	int count = 0; 
	
	
	/// X-Y K_stop Centroid Variables
	int total_energy = 0;
	double X_weights = 0.0;
	double Y_weights = 0.0;
	
	
	// Run through events  
	for(int ievt=first_event; ievt <= last_event; ievt++){
		
	  if(ievt%1000 == 0) cout << "On event: " << ievt << endl;	
	  
	 /////////////////////////////////////////////////////////////////////////////////////////////////////

	
	  Good_Event=false;
	  bool Good_tof1[12] = {false};
	  bool Good_tof2[12] = {false};
	  
	  
	  // Get new event data, and determine if the event is a good event.
	  for(int ivt=ievt; ivt<ievt+1; ivt++){
	    fChain->GetEntry(ivt);  
	  
	  for (int j_TARGET=0; j_TARGET<256; j_TARGET++){
	    ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-HG_TARGET_ADC_Thr[j_TARGET];
	    ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-LG_TARGET_ADC_Thr[j_TARGET];
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
	        (ADC_High_TARGET[i]>=0 && tdc_le_target[i][3]>=TARGET_TDC_min[i] && tdc_le_target[i][3]<=TARGET_TDC_max[i]))
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
	
   
	    
	    
	  // Check if current event is on the blacklist/.
	  Event_On_Blacklist = false;
	    
      for(vector<int>::iterator it = blacklist_events.begin(); it != blacklist_events.end(); it++){
		if(ievt == *it){
	      Event_On_Blacklist = true;
		  break;
		}
	  }  
		
		
		
		
	
	  for(int i=0; i<12; i++){
	
	    if(ADC_TOF1U[i]>=0 || ADC_TOF1D[i]>=0)  Good_TOF1_ADC[i] = true;
	      
	    if((TDC_TOF1U[i]>=TOF1U_TDC_min[i] && TDC_TOF1U[i]<=TOF1U_TDC_max[i]) ||
	       (TDC_TOF1D[i]>=TOF1D_TDC_min[i] && TDC_TOF1D[i]<=TOF1D_TDC_max[i]))  Good_TOF1_TDC[i] = true;
	
	    if(Good_TOF1_TDC[i]) Good_TOF1[i] = true;
	
	
	    if(ADC_TOF2AO[i]>=0 || ADC_TOF2AI[i]>=0 || ADC_TOF2BO[i]>=0 || ADC_TOF2BI[i]>=0)  Good_TOF2_ADC[i] = true;
	
	    if(((TDC_TOF2AO[i]>=TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<=TOF2AO_TDC_max[i])  &&
	        (TDC_TOF2AI[i]>=TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<=TOF2AI_TDC_max[i])) ||
	       ((TDC_TOF2BO[i]>=TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<=TOF2BO_TDC_max[i])  &&
	        (TDC_TOF2BI[i]>=TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<=TOF2BI_TDC_max[i])))  Good_TOF2_TDC[i] = true;
	
	    if(Good_TOF2_ADC[i] && Good_TOF2_TDC[i]) Good_TOF2[i] = true;
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
	      fprintf(Event_K_Stop_Position, "%-12d NO MWPC\n", ievt);
	      break;
	    }
	
	    if(Good_TARGET_Event && !Good_TOF_Event && Good_MWPC_Event){
	      fprintf(Event_K_Stop_Position,"%-12d  NO TOF\n",ievt);
	      break;
	    }
	    
	    if(!Good_TARGET_Event && Good_TOF_Event && Good_MWPC_Event){
	      fprintf(Event_K_Stop_Position,"%-12d  NO TARGET\n",ievt);
	      break;
	    }
	
	    if(Good_TARGET_Event && !Good_TOF_Event && !Good_MWPC_Event){
	      fprintf(Event_K_Stop_Position,"%-12d  NO TOF, NO MWPC\n",ievt);
	      break;
	    }
	
	    if(!Good_TARGET_Event && Good_TOF_Event && !Good_MWPC_Event){
	      fprintf(Event_K_Stop_Position,"%-12d  NO TARGET, NO MWPC\n",ievt);
	      break;
	    }
	
	    if(!Good_TARGET_Event && !Good_TOF_Event && Good_MWPC_Event){
	      fprintf(Event_K_Stop_Position,"%-12d  NO TARGET, NO TOF\n",ievt);
	      break;
	    }
	
	    if(!Good_TARGET_Event && !Good_TOF_Event && !Good_MWPC_Event){
	      fprintf(Event_K_Stop_Position,"%-12d  NO TOF, NO TARGET, NO MWPC\n",ievt);
	      break;
	    }
	    
	    if(Event_On_Blacklist){
		  fprintf(Event_K_Stop_Position,"%-12d  Blacklist Event\n",ievt);
	      break;
		}

	    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  }
	
	
	  if(!Good_Event){
		continue;  
	  } 
	  
	  count = 0;
	  
	  
	  bool has_TDC_hit[256] = {false};
	  for(Int_t i=0; i<256; i++){
	    for (Int_t k=0; k<4; k++) {
	      if ((tdc_le_target[i][k]>=TARGET_TDC_min[i]) && (tdc_le_target[i][k]<=TARGET_TDC_max[i])) has_TDC_hit[i] = true;
	   }
	  }	  
	  
	  
	  
	  // Determine if a hit target has any hit neighbours
	  bool TARGET_High_has_neighbours[256] = {false};
      for(int i = 0; i<256; i++){
	    for(int j=0; j<8; j++){
	      if((TARGET_neighbours[i][j] != -1 && ADC_High_TARGET[TARGET_neighbours[i][j]] >= Angle_ADC_cut && has_TDC_hit[TARGET_neighbours[i][j]]) ||
		     (TARGET_neighbours[i][j] != -1 && ADC_High_TARGET[TARGET_neighbours[i][j]] <0 && ADC_Low_TARGET[TARGET_neighbours[i][j]]>=0 && Switch==1)){
              TARGET_High_has_neighbours[i] = true;
	          break;
	      }
	    }
      }	  
	 
	  
	  
	  for(int i = 0; i<256; i++){
	    if(TARGET_High_has_neighbours[i]){
	      if(ADC_High_TARGET[i]>=Angle_ADC_cut && has_TDC_hit[i])
	        count++;
	    }	  
      }
	  
	  if(count<n_hit){
	    fprintf(Event_K_Stop_Position,"%-12d has less than %d hits in the TARGET\n", ievt, n_hit);
        continue;
      }	  
	  
	  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	  /// Determine bar with maximum low gain adc
	  int max_energy = -100000;
	  int low_gain_bar_max = -1;
	  
	  for(int i = 0; i<256; i++){
		if(ADC_Low_TARGET[i] >= max_energy){
		  low_gain_bar_max = i;  
		  max_energy = ADC_Low_TARGET[i];
	    }
	  }
	    	  
	  /// Determine K Stop Position
	  
	  
	  bool has_TDC_hit_Kstop[256] = {false};

      for(Int_t i=0; i<256; i++){
        for (Int_t k=0; k<4; k++) {
          if ((tdc_le_target[i][k]>=TDC_min_Kstop) && (tdc_le_target[i][k]<=TDC_max_Kstop)) has_TDC_hit_Kstop[i] = true;
        }
      }
	  
	  vector<int> good_k_stop_bars;
	  
	  for(int i = 0; i<256; i++){
		if(ADC_High_TARGET[i] > HG_KAON && ADC_Low_TARGET[i] > LG_KAON && has_TDC_hit_Kstop[i])
		  good_k_stop_bars.push_back(i);
	  }
	  
	  
	  total_energy = 0;
	  X_weights = 0.0;
	  Y_weights = 0.0;
	  
	  // Compute energy weighted centroids
	  if(good_k_stop_bars.empty()){
	    fprintf(Event_K_Stop_Position,"%-12d No Good K Stop Bars\n", ievt);
	  }
	  else{
		for(vector<int>::iterator it = good_k_stop_bars.begin(); it != good_k_stop_bars.end(); it++){
		  X_weights += ADC_Low_TARGET[*it]*Xloc[*it];
		  Y_weights += ADC_Low_TARGET[*it]*Yloc[*it];
		  total_energy += ADC_Low_TARGET[*it];
	    }  
	  //                            "Event       Xbar       Ybar       ADC Low Max bar"  
	  fprintf(Event_K_Stop_Position,"%-6d  %8.3f     %8.3f             %-4d\n",ievt, X_weights/double(total_energy), Y_weights/double(total_energy), low_gain_bar_max);	    
	  }
	  
	  
	  

	  
	  
	  
	  good_k_stop_bars.clear();

	}//End for loop for entries;
	
	
  fclose(Event_K_Stop_Position);


  delete fChain;
	    
  return;
} // End void
	
	
