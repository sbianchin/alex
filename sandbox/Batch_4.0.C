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
//#include "ADC_Thresholds.h"
#include "ADC_TARGET_Pedestals.h"
#include "TDC_Windows.h"
#include "Cuts_and_Windows.h"
//#include "MWPC_Thr.h"
#include "MWPC_Thr2.h"
#endif

#include "intersect.cxx"
#include "C2_Strip_transform.h"
#include "Channel_to_Strip.h"
//#include "SFT_functions.h"
#include "Batch_Variables.h"
#include "Batch_Event_Display_functions.h"

using namespace std;
 
int kaon_fiber(double x_bar, double y_bar);
void Batch_4_0(Int_t Run_Number=5, Int_t ievt_min=2, Int_t ievt_max=2, Int_t Switch_Display=0){ 
  

  gROOT->SetBatch(1);

  gStyle->Clear();
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(11);

  bool flag = false;

  vector<int> nan_evts;   //nan events
  nan_evts.clear();
 

  if(Switch_Display==0 || Switch_Display==1 || Switch_Display==2){
    flag = true;
  }
  if(!flag){
    cout << endl;
    cout << "Flag Error! " << endl;
    cout << "Please, choose another flag" << endl;
    cout << "0: Screen Printout (Reduced Statistics)" << endl;
    cout << "1: File Printout (Reduced Statistics) " << endl;
    cout << "2: File Printout (Full Statistics)" << endl;
    cout << endl;
  }

  /////////////////////   Dave's Time Walk Correction File  ////////////////////////

  //Int_t nBars = 256;
  double par_in[256][3] = {0.};
  double par_err[356][3] = {0.};
  Int_t ADCmax = 3450;
  //Int_t HGpedOffset = 10;
  double Yfi = 0.;
  double Ani = 0.;
  double Yni = 0.;
  double Tpi = 0.;
  //Int_t TDCnew[256][16] = {-1};
  bool NoFile = false;

  char ParsTarg1[100];
  sprintf(ParsTarg1,"TimeWalk3994.dat");
  //sprintf(ParsTarg1,"TimeWalk%d.dat",Run_Number);

  if(!ifstream(ParsTarg1)) NoFile = true;
 
  ///////////////////////////////////////////////////////////////////////////////////

  Int_t TDC_min_TARGET = TARGET_TDC_min[0];
  Int_t TDC_max_TARGET = TARGET_TDC_max[0];

  sprintf(ADC_cut,"(ADC >= %d)",HG_SFT_ADC_Thr[0]);

  sprintf(file_mapping,"../Mapping");

  sprintf(par_finput,"%s/%s",file_mapping,source_mapping);

  sprintf(par_finput2,"%s/MWPC_map2.txt",file_mapping);
  ifstream fdat(par_finput,ios::in);
  for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
  fdat.close();

  ifstream fdat2(par_finput2,ios::in);
  for(Int_t ii=0; ii<512; ii++) fdat2 >> par_temp2[ii];
  fdat2.close();
  
  sprintf(h_Gap_Fibers_Name,"Event %d (Run %d)", ievt_min, Run_Number);
  sprintf(h_Gap_Fibers_Title,"ADC LG offset =  %d  |  %d #leq TDC #leq %d", TARGET_ADC_Thr_LG_Offset,TDC_min_TARGET,TDC_max_TARGET);

   
  TH2F *h_TDC_Gap_Fibers = new TH2F(h_Gap_Fibers_Name, h_Gap_Fibers_Title, 500, -50, 50, 500, -50, 50);

  sprintf(path_input,"%s",path_merged);          
  sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);
  sprintf(footer,"Event_Display_MS.C  --  Run %d ; Event %d",Run_Number,ievt_min);
  sprintf(run_string,"Run %d ; Event %d",Run_Number,ievt_min);
  sprintf(event_string,"Run %d ; Event %d",Run_Number,ievt_min);
  
  if(Switch_Display==0 || Switch_Display==1){
    cout << endl;
    cout << "File opened:  " << Name_finput << endl;
    if(!NoFile) cout << "Time Walk Parameter File:  " << ParsTarg1 << endl;
    if(NoFile) cout << "Time Walk Parameter File:  " << "NOT FOUND !" << endl;
    cout << endl;
  } 

  TChain *fChain= new TChain("Tree");   
  fChain->Add(Name_finput);   
  fChain->SetMakeClass(1);              
  
  fChain->SetBranchAddress("TDC_Trig",tdc_trigger);

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

  fChain->SetBranchAddress("TDC_Ck", tdc_ck);
  fChain->SetBranchAddress("TDC_Cpi", tdc_cpi);

  fChain->SetBranchAddress("EvFlag", Event_flag);
 
  if(Switch_Display==2){
    ievt_max = fChain->GetEntries()-1;
    //cout << "DEBUG:  " << ievt_max << endl;
  }
  
  int count = 0; 

  x_TDC_Gap_Fibers_intersect1 = 999.99;   x_TARGET_intersect = 999.99;
  y_TDC_Gap_Fibers_intersect1 = 999.99;   y_TARGET_intersect = 999.99;

  ZZ.clear();
  DZ.clear();
  Track_Length.clear();
  Z_selected.clear();
  Z_selected_sorted.clear();

  vec_tdc_ck.clear();
  vec_tdc_ck_col.clear();

  sprintf(output,"RUN_%d_DATA.txt",Run_Number);
  ofstream fout;
  fout.open(output);

  sprintf(output2,"Run_%d_Nan_Events.txt",Run_Number);
  ofstream fout2;
  fout2.open(output2);

  cout << endl;



  for(int ivt=ievt_min; ivt<ievt_max+1; ivt++){ // LOOP OVER EVENTS
    ZZ.clear();
    fChain->GetEntry(ivt);  
    K_stop_count = 0;

    vec_xx_lepton.clear();
    vec_yy_lepton.clear();

    if((ivt-ievt_min)%10000==0 && (ivt-ievt_min)!=0) cout << "***  " << (ivt-ievt_min) << " events done!" << endl;

    for(int i=0; i<256; i++)  HG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Ped_HG[i]) + TARGET_ADC_Thr_HG_Offset;
    for(int i=0; i<256; i++)  LG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Ped_LG[i]) + TARGET_ADC_Thr_LG_Offset;
    for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;
    
    for (int j_TARGET=0; j_TARGET<256; j_TARGET++){
      ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-HG_TARGET_ADC_Thr[j_TARGET];
      ADC_High_TARGET_ped[j_TARGET]=adc_high_target[j_TARGET]-TARGET_ADC_Ped_HG[j_TARGET];
      ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-LG_TARGET_ADC_Thr[j_TARGET];
      ADC_Low_TARGET_ped[j_TARGET]=adc_low_target[j_TARGET]-TARGET_ADC_Ped_LG[j_TARGET]; 
    } 
  
      ////////////////////////  Dave's Time Walk Correction ///////////////////////

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

      //if((Ani >= TARGET_ADC_Ped_HG[ii]) && (Ani < ADCmax)){ // SEB  /// <----- here !
      if((Ani >= TARGET_ADC_Thr_HG_Offset) && (Ani < ADCmax)){ // SEB
        Yni = par_in[ii][0] - par_in[ii][1]/sqrt(Ani - par_in[ii][2]);
        Tpi = Yfi - Yni;
        for(int jj=0; jj<16; jj++){
          if(tdc_le_target[ii][jj]>0) tdc_le_target[ii][jj] = tdc_le_target[ii][jj] + Tpi;
        }
      }
    }

    /////////////////////////////////////////////////////////////////////////////
    // Ck and C_pi
    TDC_Ck_counter = 0;
    TDC_Cpi_counter = 0;

    vec_tdc_ck.clear();
    TDC_ck_avg = 0.;   TDC_ck_avg2 =0.;
    TDC_ck_sum = 0;    TDC_ck_sum2 = 0;
    TDC_ck_sigma = 0.; TDC_ck_sigma2 = 0.;
    TDC_ck_counter = 0;
    ck_mean = -1;
    TDC_diff = -1; 

    for(int ic = 0; ic < 14; ic++){
      if(((tdc_ck[ic][0] >= TDC_Ck_min) && (tdc_ck[ic][0] <= TDC_Ck_max))||
         ((tdc_ck[ic][1] >= TDC_Ck_min) && (tdc_ck[ic][1] <= TDC_Ck_max))||
         ((tdc_ck[ic][2] >= TDC_Ck_min) && (tdc_ck[ic][2] <= TDC_Ck_max))||
         ((tdc_ck[ic][3] >= TDC_Ck_min) && (tdc_ck[ic][3] <= TDC_Ck_max))||
         ((tdc_ck[ic][4] >= TDC_Ck_min) && (tdc_ck[ic][4] <= TDC_Ck_max))||
         ((tdc_ck[ic][5] >= TDC_Ck_min) && (tdc_ck[ic][5] <= TDC_Ck_max))) TDC_Ck_counter++;
         
      if(((tdc_cpi[ic][0] >= TDC_Cpi_min) && (tdc_cpi[ic][0] <= TDC_Cpi_max))||
         ((tdc_cpi[ic][1] >= TDC_Cpi_min) && (tdc_cpi[ic][1] <= TDC_Cpi_max))||
         ((tdc_cpi[ic][2] >= TDC_Cpi_min) && (tdc_cpi[ic][2] <= TDC_Cpi_max))||
         ((tdc_cpi[ic][3] >= TDC_Cpi_min) && (tdc_cpi[ic][3] <= TDC_Cpi_max))||
         ((tdc_cpi[ic][4] >= TDC_Cpi_min) && (tdc_cpi[ic][4] <= TDC_Cpi_max))||
         ((tdc_cpi[ic][5] >= TDC_Cpi_min) && (tdc_cpi[ic][5] <= TDC_Cpi_max))) TDC_Cpi_counter++;
    
      for(int jc=0; jc<16; jc++){
        if(tdc_ck[ic][jc]>=TDC_Ck_min && tdc_ck[ic][jc]<=TDC_Ck_max){
          vec_tdc_ck_col.push_back(tdc_ck[ic][jc]);
        }
      }
      if(vec_tdc_ck_col.size()>0) vec_tdc_ck.push_back(vec_tdc_ck_col[vec_tdc_ck_col.size()-1]);
      vec_tdc_ck_col.clear();
    }

     for(unsigned int ik=0; ik<vec_tdc_ck.size(); ik++){
      TDC_ck_sum += vec_tdc_ck[ik];
    }

    if(vec_tdc_ck.size()>0) TDC_ck_avg = double(TDC_ck_sum) / double(vec_tdc_ck.size());
    else TDC_ck_avg = -1; 

    for(unsigned int i=0; i<vec_tdc_ck.size(); i++) TDC_ck_sigma2 += pow((vec_tdc_ck[i]-TDC_ck_avg),2);
  
    //cout << "TEST : " << ivt << "  "
    //     << vec_tdc_ck.size() << "  " 
    //     << TDC_ck_sigma2 << endl;


    if(vec_tdc_ck.size()>0) TDC_ck_sigma = sqrt(TDC_ck_sigma2/vec_tdc_ck.size());
    else TDC_ck_sigma = -1;
      
    for(unsigned int i=0; i<vec_tdc_ck.size(); i++){
      if(abs(vec_tdc_ck[i]-TDC_ck_avg) <= 1.4*TDC_ck_sigma){
        TDC_ck_sum2 += vec_tdc_ck[i];
        TDC_ck_counter++;
      }
    }

    if(TDC_ck_counter>0) TDC_ck_avg2 = double(TDC_ck_sum2)/double(TDC_ck_counter);
    else TDC_ck_avg2 = -1;

    if(vec_tdc_ck.size()>0) ck_mean = double(TDC_ck_sum) / double(vec_tdc_ck.size());
    else ck_mean = -1;
    //tdc_ck_corr = tdc_trigger[0][0] - ck_mean;
    tdc_ck_corr = tdc_trigger[0][0] - TDC_ck_avg2;



    /////////////////////////////////////////////////////////////////////////////
  
    for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
      ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-HG_SFT_ADC_Thr[j_SFT];
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
    Good_TARGET_Event = false;
    count_TARGET_evts = 0; 
    for(int i=0; i<256; i++)
    {
       if((ADC_High_TARGET[i]>0 && tdc_le_target[i][0]>=TARGET_TDC_min[i] && tdc_le_target[i][0]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]>0 && tdc_le_target[i][1]>=TARGET_TDC_min[i] && tdc_le_target[i][1]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]>0 && tdc_le_target[i][2]>=TARGET_TDC_min[i] && tdc_le_target[i][2]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]>0 && tdc_le_target[i][3]>=TARGET_TDC_min[i] && tdc_le_target[i][3]<=TARGET_TDC_max[i]) ||
         (ADC_High_TARGET[i]<0  && ADC_Low_TARGET[i]>0))
      { 
       count_TARGET_evts++;
      }
    }

    if(count_TARGET_evts >= n_hit) Good_TARGET_Event = true;

    
    //********* GOOD TOF EVENTS
    for(int i=0; i<12; i++){
      Good_TOF1_ADC[i] = false;   
      Good_TOF2_ADC[i] = false;
      Good_TOF1_TDC[i] = false;   
      Good_TOF2_TDC[i] = false;
      Good_TOF1[i] = false;       
      Good_TOF2[i] = false;
      Good_TOFs[i] = false;
    }
    
    Good_TOF_Event = false;

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

    //********* GOOD MWPC EVENTS  
    count_C2X = 0;    count_C2Y = 0;
    count_C3X = 0;    count_C3Y = 0;
    count_C4X = 0;    count_C4Y = 0;
    Good_MWPC_Event = false;

    for(int i=0; i<56; i++){
      if(ADC_C2X_R[i]>0. || ADC_C2X_L[i]>0.) count_C2X++;
    }

    for(int ii=0; ii<16; ii++){
      if(ADC_C2Y_R[ii]>0. || ADC_C2Y_L[ii]>0.) count_C2Y++;
    }

    for(int j=0; j<64; j++){
      if(ADC_C3X_R[j]>0. || ADC_C3X_L[j]>0.){
        count_C3X++;
      }
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

    if(Good_TARGET_Event && Good_TOF_Event && Good_MWPC_Event && !Event_On_Blacklist)  Good_Event = true;

    if(Good_TARGET_Event && Good_TOF_Event && !Good_MWPC_Event) continue;
    if(Good_TARGET_Event && !Good_TOF_Event && Good_MWPC_Event) continue;
    if(!Good_TARGET_Event && Good_TOF_Event && Good_MWPC_Event) continue;
    if(Good_TARGET_Event && !Good_TOF_Event && !Good_MWPC_Event) continue;
    if(!Good_TARGET_Event && Good_TOF_Event && !Good_MWPC_Event) continue;
    if(!Good_TARGET_Event && !Good_TOF_Event && Good_MWPC_Event) continue;
    if(!Good_TARGET_Event && !Good_TOF_Event && !Good_MWPC_Event) continue;
    
    if(!Good_Event) return;

    for(int j=0 ; j<128 ; j++){ 
      if(ADC_High_SFT[j]<0)      ADC_High_SFT_corr[j]=0; 
      if(ADC_High_SFT[j]>=0)     ADC_High_SFT_corr[j]=ADC_High_SFT[j]; 
    }

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

    //////// DETERMINE FIBER WITH HIGHEST AND SECOND HIGHEST LOW GAIN AMPLITUDE
    Int_t TDC_average = -1;

    int max_index_all[256];
    int max_ADC_all[256];
    int max_index_flag;
  
    for(int i = 0; i<256; i++){
      max_ADC_all[i] = -100000000;
      max_index_all[i] = -1;
    }
  
    for(int j=0; j<256; j++){
      for(int i=0; i<256; i++){
        max_index_flag = 0;  
  
        for(int k = 0; k<256; k++){
          if(i == max_index_all[k]){
            max_index_flag = 1;
            break;
          }
        }  
        if(max_index_flag == 1) continue;     
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

    for(int i = 0; i<256; i++){
      for(int j=0; j<6; j++){  
        if(TDC_LG_max == -1){ 
          if(tdc_le_target[max_index_all[i]][j] >= TDC_Thr_min && tdc_le_target[max_index_all[i]][j]  <= TDC_Thr_max)
            TDC_LG_max = tdc_le_target[max_index_all[i]][j];
        }
        else if(TDC_LG_max2 == -1){
          if(tdc_le_target[max_index_all[i]][j] >= TDC_Thr_min && tdc_le_target[max_index_all[i]][j] <= TDC_Thr_max)
            TDC_LG_max2 = tdc_le_target[max_index_all[i]][j];
        }
      }     
    } 
  
    if(abs(TDC_LG_max - TDC_LG_max2) <= T_limit) TDC_average = (TDC_LG_max + TDC_LG_max2)/2;
    else if(TDC_LG_max > TDC_LG_max2) TDC_average = TDC_LG_max;
    else if(TDC_LG_max < TDC_LG_max2) TDC_average = TDC_LG_max2;

    TDC_diff = 550 - (0.91*TDC_average - 0.625*tdc_ck_corr);

    //cout << "TEST : " << ivt << "  " << TDC_diff << endl;


    int kaon_TDC_min = TDC_average + TDC_Avg_Offset_min;
    int kaon_TDC_max = TDC_average + TDC_Avg_Offset_max;
    int TDC_min_Kstop = TDC_average + TDC_Kstop_Avg_Offset_min;
    int TDC_max_Kstop = TDC_average + TDC_Kstop_Avg_Offset_max;
  
    TDC_min_TARGET = kaon_TDC_min;
    TDC_max_TARGET = kaon_TDC_max;

    for(int i=0; i<256; i++){
      has_TDC_hit[i] = false;
      has_TDC_hit_Kstop[i] = false;
    }

    for(Int_t i=0; i<256; i++){
      for (Int_t k=0; k<4; k++) {
        if ((tdc_le_target[i][k]>=TDC_min_Kstop) && (tdc_le_target[i][k]<=TDC_max_Kstop)) has_TDC_hit_Kstop[i] = true;
        if ((tdc_le_target[i][k]>=kaon_TDC_min) && (tdc_le_target[i][k]<=kaon_TDC_max)) has_TDC_hit[i] = true;
      }
    }

    sprintf(ch_ADC_cut_TARGET,"(ADC offset = %d)",TARGET_ADC_Thr_HG_Offset);
    sprintf(ch_ADC_and_TDC_cut,"(ADC offset = %d | %d #leq TDC #leq %d)",TARGET_ADC_Thr_HG_Offset,TDC_min_TARGET,TDC_max_TARGET);
    sprintf(ch_ADC_and_TDC_cut_Kstop,"(ADC: HG #geq %d, LG #geq %d | %d #leq TDC K Stop #leq %d)",HG_KAON,LG_KAON,TDC_min_Kstop,TDC_max_Kstop);

    Angle_ADC_cut = 0;

    x_inc = 0;
    y_inc = 0;
    hit_count = 0;
    count = 0; 

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

    x_cent = Xloc[max_index];   y_cent = Yloc[max_index];
  
    for(Int_t j=0; j<256; j++){
      hyp[j] = -1;
      hyp[j] = sqrt(pow(x_cent - Xloc[j],2) + pow(y_cent - Yloc[j],2));
    }

    for(int i=0; i<12; i++){
      has_ADC_TOF1_hit[i] = false;
      has_TDC_TOF1_hit[i] = false;
      has_ADC_TOF2_hit[i] = false;
      has_TDC_TOF2_hit[i] = false;
    }

    ///Set TOF2 Lines
    for(int i = 0; i < 12; i++){
      if ((ADC_TOF2AO[i]>0 && ADC_TOF2AI[i]>0) || (ADC_TOF2BO[i]>0 && ADC_TOF2BI[i]>0)) {has_ADC_TOF2_hit[i]=true;}
      if ((((TDC_TOF2AO[i]>TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<TOF2AO_TDC_max[i]) && (TDC_TOF2AI[i]>TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<TOF2AI_TDC_max[i]))) 
      ||  (((TDC_TOF2BO[i]>TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<TOF2BO_TDC_max[i]) && (TDC_TOF2BI[i]>TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<TOF2BI_TDC_max[i])))) {has_TDC_TOF2_hit[i]=true;}
    }

    ///Set TOF1 Lines
    for(int i = 0; i < 12; i++){
      if (ADC_TOF1U[i]>0 || ADC_TOF1D[i]>0) {has_ADC_TOF1_hit[i] = true;}
      if ((TDC_TOF1U[i]>TOF1U_TDC_min[i] && TDC_TOF1U[i]<TOF1U_TDC_max[i]) || (TDC_TOF1D[i]>TOF1D_TDC_min[i] && TDC_TOF1D[i]<TOF1D_TDC_max[i])) {has_TDC_TOF1_hit[i] = true;}
    }
  
    for(int i=0; i<12; i++){
      ADC_TOF1_hit[i] = 0;
      ADCTDC_TOF1_hit[i] = 0;
      ADC_TOF2_hit[i] = 0;
      ADCTDC_TOF2_hit[i] = 0;
    }
 
    for (int k=0; k<12; k++) {
      if (has_ADC_TOF1_hit[k]) {
        if (has_TDC_TOF1_hit[k]) ADCTDC_TOF1_hit[k]++;
        else {ADC_TOF1_hit[k]++;}
      }
      if (has_ADC_TOF2_hit[k]) {
        if (has_TDC_TOF2_hit[k]) ADCTDC_TOF2_hit[k]++;
        else {ADC_TOF2_hit[k]++;}
      }
    }
  
    selected_TOF2 = 0;
  
    // Determine which TOF2 is hit
    for(int i = 0; i<12; i++){
      if(has_TDC_TOF2_hit[i] && has_ADC_TOF2_hit[i]) selected_TOF2 = i + 1;  
    }

    if(selected_TOF2 == 0){
      for(int i = 0; i<12; i++){
        if(has_TDC_TOF2_hit[i] || has_ADC_TOF2_hit[i]) selected_TOF2 = i+1;
      }
    }

  for(int i=0; i<12; i++) gap_counter[i] = 0;

  
  ////////// Clustering (C2X_centroid)
  for(int i=0; i<56; i++){
    C2X_L[i] = 0.;
    C2X_R[i] = 0.;
  }

  for(int j=0; j<64; j++){
    C3X_L[j] = 0.;
    C3X_R[j] = 0.;
  }

  for(int k=0; k<72; k++){
    C4X_L[k] = 0.;
    C4X_R[k] = 0.;
  }

  for(int l=0; l<16; l++){
    C2Y_L[l] = 0.;
    C2Y_R[l] = 0.;
    C3Y_L[l] = 0.;
    C3Y_R[l] = 0.;
    C4Y_L[l] = 0.;
    C4Y_R[l] = 0.;
  }

  C2X_centroid = 0.0;

  
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
        //cout << " ADC MWPC Channel " << C2XL_Channels[q] << " -- " << "C2_X_" << q+1 << "_Left: " << int(ADC_C2X_L[q]) << endl;
      }
    }  
    for(int q=0; q<16; q++){  
      if(Good_Event && ADC_C2Y_L[q]>0){
        //cout << " ADC MWPC Channel " << C2YL_Channels[q] << " -- " << "C2_Y_" << q+1 << "_Left: " << int(ADC_C2Y_L[q]) << endl;
      }
    }
    for(int q=0; q<56; q++){ // C2 Counters
      if(Good_Event && ADC_C2X_R[q]>0){
        //cout << " ADC MWPC Channel " << C2XR_Channels[q] << " -- " << "C2_X_" << q+1 << "_Right: " << int(ADC_C2X_R[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){  
      if(Good_Event && ADC_C2Y_R[q]>0){
        //cout << " ADC MWPC Channel " << C2YR_Channels[q] << " -- " << "C2_Y_" << q+1 << "_Right: " << int(ADC_C2Y_R[q]) << endl;
      }
    }

    //cout << " " << endl;

    for(int q=0; q<64; q++){ // C3 Counters
      if(Good_Event && ADC_C3X_L[q]>0){
        //cout << " ADC MWPC Channel " << C3XL_Channels[q] << " -- " << "C3_X_" << q+1 << "_Left: " << int(ADC_C3X_L[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){  
      if(Good_Event && ADC_C3Y_L[q]>0){
        //cout << " ADC MWPC Channel " << C3YL_Channels[q] << " -- " << "C3_Y_" << q+1 << "_Left: " << int(ADC_C3Y_L[q]) << endl;
      }
    }
    for(int q=0; q<64; q++){ // C3 Counters
      if(Good_Event && ADC_C3X_R[q]>0){
        //cout << " ADC MWPC Channel " << C3XR_Channels[q] << " -- " << "C3_X_" << q+1 << "_Right: " << int(ADC_C3X_R[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){ // C3 Counters
      if(Good_Event && ADC_C3Y_R[q]>0){
        //cout << " ADC MWPC Channel " << C3YR_Channels[q] << " -- " << "C3_Y_" << q+1 << "_Right: " << int(ADC_C3Y_R[q]) << endl;
      }
    }

    //cout << " " << endl;

    for(int q=0; q<72; q++){ // C4 Counters
      if(Good_Event && ADC_C4X_L[q]>0){
        //cout << " ADC MWPC Channel " << C4XL_Channels[q] << " -- " << "C4_X_" << q+1 << "_Left: " << int(ADC_C4X_L[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){   
      if(Good_Event && ADC_C4Y_L[q]>0){
        //cout << " ADC MWPC Channel " << C4YL_Channels[q] << " -- " << "C4_Y_" << q+1 << "_Left: " << int(ADC_C4Y_L[q]) << endl;
      }
    }
    for(int q=0; q<72; q++){ // C4 Counters
      if(Good_Event && ADC_C4X_R[q]>0){
        //cout << " ADC MWPC Channel " << C4XR_Channels[q] << " -- " << "C4_X_" << q+1 << "_Right: " << int(ADC_C4X_R[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){   
      if(Good_Event && ADC_C4Y_R[q]>0){
        //cout << " ADC MWPC Channel " << C4YR_Channels[q] << " -- " << "C4_Y_" << q+1 << "_Right: " << int(ADC_C4Y_R[q]) << endl;
      }
    }
  }
  else{
    for(int q=0; q<56; q++){ // C2 Counters
      if(Good_Event && ADC_C2X_R[q]>0){
        //cout << " ADC MWPC Channel " << C2XR_Channels[q] << " -- " << "C2_X_" << q+1 << "_Right: " << int(ADC_C2X_R[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){   
      if(Good_Event && ADC_C2Y_R[q]>0){ 
        //cout << " ADC MWPC Channel " << C2YR_Channels[q] << " -- " << "C2_Y_" << q+1 << "_Right: " << int(ADC_C2Y_R[q]) << endl;
      }
    }
    for(int q=0; q<56; q++){ // C2 Counters
      if(Good_Event && ADC_C2X_L[q]>0){ 
        //cout << " ADC MWPC Channel " << C2XL_Channels[q] << " -- " << "C2_X_" << q+1 << "_Left: " << int(ADC_C2X_L[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){   
      if(Good_Event && ADC_C2Y_L[q]>0){
        //cout << " ADC MWPC Channel " << C2YL_Channels[q] << " -- " << "C2_Y_" << q+1 << "_Left: " << int(ADC_C2Y_L[q]) << endl;
      }
    }

    //cout << " " << endl;

    for(int q=0; q<64; q++){ // C3 Counters
      if(Good_Event && ADC_C3X_R[q]>0){
        //cout << " ADC MWPC Channel " << C3XR_Channels[q] << " -- " << "C3_X_" << q+1 << "_Right: " << int(ADC_C3X_R[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){   
      if(Good_Event && ADC_C3Y_R[q]>0){
        //cout << " ADC MWPC Channel " << C3YR_Channels[q] << " -- " << "C3_Y_" << q+1 << "_Right: " << int(ADC_C3Y_R[q]) << endl;
      }
    }
    for(int q=0; q<64; q++){ // C3 Counters
      if(Good_Event && ADC_C3X_L[q]>0){
        //cout << " ADC MWPC Channel " << C3XL_Channels[q] << " -- " << "C3_X_" << q+1 << "_Left: " << int(ADC_C3X_L[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){   
      if(Good_Event && ADC_C3Y_L[q]>0){
        //cout << " ADC MWPC Channel " << C3YL_Channels[q] << " -- " << "C3_Y_" << q+1 << "_Left: " << int(ADC_C3Y_L[q]) << endl;
      }
    }

    //cout << " " << endl;

    for(int q=0; q<72; q++){ // C4 Counters
      if(Good_Event && ADC_C4X_R[q]>0){
        //cout << " ADC MWPC Channel " << C4XR_Channels[q] << " -- " << "C4_X_" << q+1 << "_Right: " << int(ADC_C4X_R[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){   
      if(Good_Event && ADC_C4Y_R[q]>0){
        //cout << " ADC MWPC Channel " << C4YR_Channels[q] << " -- " << "C4_Y_" << q+1 << "_Right: " << int(ADC_C4Y_R[q]) << endl;
      }
    }
    for(int q=0; q<72; q++){ // C4 Counters
      if(Good_Event && ADC_C4X_L[q]>0){
        //cout << " ADC MWPC Channel " << C4XL_Channels[q] << " -- " << "C4_X_" << q+1 << "_Left: " << int(ADC_C4X_L[q]) << endl;
      }
    }
    for(int q=0; q<16; q++){   
      if(Good_Event && ADC_C4Y_L[q]>0){
        //cout << " ADC MWPC Channel " << C4YL_Channels[q] << " -- " << "C4_Y_" << q+1 << "_Left: " << int(ADC_C4Y_L[q]) << endl;
      }
    }
  }

  
  //cout << endl;
  //cout << endl;

  // Find clustering of mwpcs
  ////////////////////////////
  first_cluster = true;
  cluster_spacing = 0;
  cluster_length_count = 0;

  C2X_clusters = 0;
  C2Y_clusters = 0;
  C3X_clusters = 0;
  C3Y_clusters = 0;
  C4X_clusters = 0;
  C4Y_clusters = 0;

  C2X_cluster_index.clear(); // Hold the starting index of each cluster.
  C2X_cluster_length.clear(); // Hold the length of each cluster.
  C2Y_cluster_index.clear();
  C2Y_cluster_length.clear();

  C3X_cluster_index.clear();
  C3X_cluster_length.clear();
  C3Y_cluster_index.clear();
  C3Y_cluster_length.clear();    

  C4X_cluster_index.clear();
  C4X_cluster_length.clear();  
  C4Y_cluster_index.clear();
  C4Y_cluster_length.clear();  

  
  if(selected_TOF2 > 6){  // LEFT
    // count C2X clusters
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

        if(i+2 < 56 && C2X_L[i+1] <= 0. && C2X_L[i+2] <= 0.){
          C2X_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 2 == 56 && C2X_L[i+1] <= 0.){
          C2X_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 1 == 56){
          C2X_cluster_length.push_back(cluster_length_count);
        }
      }
      else{
        cluster_spacing++;
        if(i != 0 && C2X_L[i-1] <= 0.){
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
          C2Y_clusters++;
          C2Y_cluster_index.push_back(i);
        }
        cluster_length_count++;
        cluster_spacing = 0;

        if(i+2 < 16 && C2Y_L[i+1] <= 0 && C2Y_L[i+2] <= 0){
          C2Y_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 2 == 16 && C2Y_L[i+1] <= 0){
          C2Y_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 1 == 16){
          C2Y_cluster_length.push_back(cluster_length_count);
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

        if(i+2 < 64 && C3X_L[i+1] <= 0 && C3X_L[i+2] <= 0){
          C3X_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 2 == 64 && C3X_L[i+1] <= 0){
          C3X_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 1 == 64){
          C3X_cluster_length.push_back(cluster_length_count);
        }
      }
      else{
        cluster_spacing++;
        if(i != 0 && C3X_L[i-1] <= 0){
          cluster_length_count = 0;
        }
        else if ( i != 0){
          cluster_length_count++;
        }
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

        if(i+2 < 16 && C3Y_L[i+1] <= 0 && C3Y_L[i+2] <= 0){
          C3Y_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 2 == 16 && C3Y_L[i+1] <= 0){
          C3Y_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 1 == 16){
          C3Y_cluster_length.push_back(cluster_length_count);
        }
      }
      else{
        cluster_spacing++;
        if(i != 0 && C3Y_L[i-1] <= 0){
          cluster_length_count = 0;
        }
        else if ( i != 0){
          cluster_length_count++;
        }
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

        if(i+2 < 72 && C4X_L[i+1] <= 0 && C4X_L[i+2] <= 0){
          C4X_cluster_length.push_back(cluster_length_count);
        }
        else if(i == 70 && C4X_L[71] <= 0){
          C4X_cluster_length.push_back(cluster_length_count);
        }
        else if(i == 71){
          C4X_cluster_length.push_back(cluster_length_count);
        }
      }
      else{
        cluster_spacing++;
        if(i != 0 && C4X_L[i-1] <= 0){
          cluster_length_count = 0;
        }
        else if ( i != 0){
          cluster_length_count++;
        }
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

        if(i+2 < 16 && C4Y_L[i+1] <= 0 && C4Y_L[i+2] <= 0){
          C4Y_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 2 == 16 && C4Y_L[i+1] <= 0){
          C4Y_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 1 == 16){
          C4Y_cluster_length.push_back(cluster_length_count);
        }
      }
      else{
        cluster_spacing++;
        if(i != 0 && C4Y_L[i-1] <= 0){
          cluster_length_count = 0;
        }
        else if ( i != 0){
          cluster_length_count++;
        }
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

        if(i+2 < 56 && C2X_R[i+1] <= 0 && C2X_R[i+2] <= 0){
          C2X_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 2 == 56 && C2X_R[i+1] <= 0){
          C2X_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 1 == 56){
          C2X_cluster_length.push_back(cluster_length_count);
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
          C2Y_clusters++;
          C2Y_cluster_index.push_back(i);
        }
        cluster_length_count++;
        cluster_spacing = 0;

        if(i+2 < 16 && C2Y_R[i+1] <= 0 && C2Y_R[i+2] <= 0){
          C2Y_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 2 == 16 && C2Y_R[i+1] <= 0){
          C2Y_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 1 == 16){
          C2Y_cluster_length.push_back(cluster_length_count);
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

        if(i+2 < 64 && C3X_R[i+1] <= 0 && C3X_R[i+2] <= 0){
          C3X_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 2 == 64 && C3X_R[i+1] <= 0){
          C3X_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 1 == 64){
          C3X_cluster_length.push_back(cluster_length_count);
        }
      }
      else{
        cluster_spacing++;
        if(i != 0 && C3X_R[i-1] <= 0){
          cluster_length_count = 0;
        }
        else if ( i != 0){
          cluster_length_count++;
        }
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

        if(i+2 < 16 && C3Y_R[i+1] <= 0 && C3Y_R[i+2] <= 0){
          C3Y_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 2 == 16 && C3Y_R[i+1] <= 0){
          C3Y_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 1 == 16){
          C3Y_cluster_length.push_back(cluster_length_count);
        }
      }
      else{
        cluster_spacing++;
        if(i != 0 && C3Y_R[i-1] <= 0){
          cluster_length_count = 0;
        }
        else if ( i != 0){
          cluster_length_count++;
        }
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

        if(i+2 < 72 && C4X_R[i+1] <= 0 && C4X_R[i+2] <= 0){
          C4X_cluster_length.push_back(cluster_length_count);
        }
        else if(i == 70 && C4X_R[71] <= 0){
          C4X_cluster_length.push_back(cluster_length_count);
        }
        else if(i == 71){
          C4X_cluster_length.push_back(cluster_length_count);
        }
      }
      else{
        cluster_spacing++;
        if(i != 0 && C4X_R[i-1] <= 0){
          cluster_length_count = 0;
        }
        else if ( i != 0){
          cluster_length_count++;
        }
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

        if(i+2 < 16 && C4Y_R[i+1] <= 0 && C4Y_R[i+2] <= 0){
          C4Y_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 2 == 16 && C4Y_R[i+1] <= 0){
          C4Y_cluster_length.push_back(cluster_length_count);
        }
        else if(i + 1 == 16){
          C4Y_cluster_length.push_back(cluster_length_count);
        }
      }
      else{
        cluster_spacing++;
        if(i != 0 && C4Y_R[i-1] <= 0){
          cluster_length_count = 0;
        }
        else if ( i != 0){
          cluster_length_count++;
        }
      }
    }   
  }  

  // Display MWPC counter information
  //cout << endl;
  //if(selected_TOF2 > 6) cout << " LEFT" << endl;
  //else cout << " RIGHT" << endl;

  //cout << " C2X clusters = " << C2X_clusters << endl;
  //cout << " C2Y clusters = " << C2Y_clusters << endl;
  //cout << " C3X clusters = " << C3X_clusters << endl;
  //cout << " C3Y clusters = " << C3Y_clusters << endl;
  //cout << " C4X clusters = " << C4X_clusters << endl;
  //cout << " C4Y clusters = " << C4Y_clusters << endl;

  C2_centroid_num = 0;
  C2_centroid_den = 0;

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


  // LENGTH IN TARGET
  
  //length_in_target = 0;
  

  /*
  if(X_BAR != -10000){
    x_tof1_intersect_1 = 0;
    y_tof1_intersect_1 = 0;
    x_tof1_intersect_2 = 0;
    y_tof1_intersect_2 = 0;
    x_tof1_intersect = 0;
    y_tof1_intersect = 0;

    alpha = angle_final_guide - 90.0;

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
    //printf("Length of track in target (xy plane) = %5.2f\n", length_in_target);
  }
  //else{
    //cout << "No K-Stop for track length." << endl;
  //}
  */























 
    //// GAP SCORING !
    scoring_type = 2;           // scoring_type = 1  --->  Oscar's Method

    //Sebastien's Method
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
  
    high_gap_hit = 0;
    gap_to_fit = 0;
    score_max = 0;

    tof1_ties.clear();

    for (int k=0; k<12; k++) {
      //if((gap_counter[k]>=high_gap_hit) && (has_TDC_TOF1_hit[k] || has_ADC_TOF1_hit[k])){
      if(gap_counter[k]>=high_gap_hit){
        if(gap_counter[k] == high_gap_hit) tof1_ties.push_back(k); 
        else{
          tof1_ties.clear();
          tof1_ties.push_back(k);
        }
        high_gap_hit = gap_counter[k];
        //gap_to_fit = k+1;
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
    for(int i=0; i<256; i++) k_stop_bar[i] = false;
    good_k_stop_bars.clear();
  
    X_BAR = 0;
    Y_BAR = 0;
   
    total_energy = 0.0;
    X_weights = 0.0;
    Y_weights = 0.0;   
    
    for(int i = 0; i<256; i++){
      if(ADC_High_TARGET[i] > HG_KAON && ADC_Low_TARGET[i] > LG_KAON && has_TDC_hit_Kstop[i]){
      good_k_stop_bars.push_back(i);
      k_stop_bar[i] = true;
      }
    }    

    // Determine if a hit target has any hit neighbours
   for(int i=0; i<256; i++) TARGET_High_has_neighbours[i] = false;

    // Determine if a hit target has any hit neighbours
    for(int i = 0; i<256; i++){
      for(int j=0; j<8; j++){
        if((TARGET_neighbours[i][j] != -1 && ADC_High_TARGET[TARGET_neighbours[i][j]] > Angle_ADC_cut && has_TDC_hit[TARGET_neighbours[i][j]] && !k_stop_bar[TARGET_neighbours[i][j]]) ||
          (TARGET_neighbours[i][j] != -1 && ADC_High_TARGET[TARGET_neighbours[i][j]] <0 && ADC_Low_TARGET[TARGET_neighbours[i][j]]>0 && !k_stop_bar[TARGET_neighbours[i][j]])){
          TARGET_High_has_neighbours[i] = true;
          break;
        }
      }
    }

    for(Int_t i=0; i<256; i++){
      if(k_stop_bar[i]) K_stop_count++;
      if(TARGET_High_has_neighbours[i]){
        if(ADC_High_TARGET[i]>Angle_ADC_cut && has_TDC_hit[i]){
        
          if(!k_stop_bar[i]){
            h_TDC_Gap_Fibers->Fill(Xloc[i],Yloc[i]);
            vec_xx_lepton.push_back(Xloc[i]);
            vec_yy_lepton.push_back(Yloc[i]);
          }
          count++;

          if(IsIn(i,channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],   // Additional weight for fibers close to the edge if hit
                    channel[gap_to_fit-1][2], channel[gap_to_fit-1][3], 
                    channel[gap_to_fit-1][4], channel[gap_to_fit-1][5], 
                    channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){
            if(!k_stop_bar[i]){         
              h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
              h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
              vec_xx_lepton.push_back(Xloc[i]);
              vec_xx_lepton.push_back(Xloc[i]);
              vec_yy_lepton.push_back(Yloc[i]);
              vec_yy_lepton.push_back(Yloc[i]);
            }        
          }
        }
    
        if(ADC_High_TARGET[i]<0 && ADC_Low_TARGET[i]>0){
          if(!k_stop_bar[i]){  
            h_TDC_Gap_Fibers->Fill(Xloc[i],Yloc[i]);
            vec_xx_lepton.push_back(Xloc[i]);
            vec_yy_lepton.push_back(Yloc[i]);
            count++;
          }
        }
      }
    }

    ///////// DETERMINE FIBER WITH HIGHEST AND SECOND HIGHEST LOW GAIN AMPLITUDE
    TDC_average = -1;

    for(int i=0; i<256; i++){
      max_index_all[i] = 999;
      max_ADC_all[i] = 999;
    }

    max_index_flag = 999;
  
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
      
        if(max_index_flag == 1) continue;   
    
        else {
          if(ADC_Low_TARGET[i]>max_ADC_all[j]){
            max_index_all[j] = i;
            max_ADC_all[j] = ADC_Low_TARGET[i];
          }
        }
      } 
    } 

    //Calculate TDC Average
    TDC_LG_max = -1;
    TDC_LG_max2 = -1;
  
    for(int i = 0; i<256; i++){
      for(int j=0; j<6; j++){  
        if(TDC_LG_max == -1){ 
          if(tdc_le_target[max_index_all[i]][j] >= TDC_Thr_min && tdc_le_target[max_index_all[i]][j]  <= TDC_Thr_max)
            TDC_LG_max = tdc_le_target[max_index_all[i]][j];
        }
        else if(TDC_LG_max2 == -1){
          if(tdc_le_target[max_index_all[i]][j] >= TDC_Thr_min && tdc_le_target[max_index_all[i]][j] <= TDC_Thr_max)
            TDC_LG_max2 = tdc_le_target[max_index_all[i]][j];
        }
      }     
    } 
  
    if(abs(TDC_LG_max - TDC_LG_max2) <= T_limit) TDC_average = (TDC_LG_max + TDC_LG_max2)/2;
    else if(TDC_LG_max > TDC_LG_max2) TDC_average = TDC_LG_max;
    else if(TDC_LG_max < TDC_LG_max2) TDC_average = TDC_LG_max2;
 
    //cout << "TEST : " << ivt << "  " << TDC_average << endl;

    //// Don't do anything if the event has less than n_hits hits in the TARGET
    if(count<n_hit){
      gROOT->Reset();
      gROOT->Clear();
      continue;
    }
      
    // Compute energy weighted centroids
    for(vector<int>::iterator it = good_k_stop_bars.begin(); it != good_k_stop_bars.end(); it++){
      X_weights += ADC_Low_TARGET[*it]*Xloc[*it];
      Y_weights += ADC_Low_TARGET[*it]*Yloc[*it];
      total_energy += ADC_Low_TARGET[*it];
    }  
  
    if(total_energy != 0){
      X_BAR = X_weights/total_energy;
      Y_BAR = Y_weights/total_energy;    
    }
    else{
      X_BAR = 999.99;
      Y_BAR = 999.99;
    }
  
    // TARGET ROTATION (90 deg.)      
    // ROTATE_CHANGE: any line tagged with ROTATE_CHANGE is used when tof1 is 6 or 12.
    gap_to_fit_rotate = gap_to_fit;
         
    if((gap_to_fit==12 || gap_to_fit==6 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){    
      for(int i=0; i<256; i++){
        ADC_High_TARGET_temp[i] = 9999;
        ADC_Low_TARGET_temp[i] = 9999;
        has_TDC_hit_temp[i] = false;
        has_TDC_hit_Kstop_temp[i] = false;
        HG_TARGET_ADC_Thr_temp[i] = 9999;
        LG_TARGET_ADC_Thr_temp[i] = 9999;
        TARGET_High_has_neighbours_temp[i] = 9999;
        k_stop_bar_temp[i] = false;
      }
  
      if(gap_to_fit==11) gap_to_fit = 2;
      else if(gap_to_fit==12) gap_to_fit = 3;
      else if(gap_to_fit==1) gap_to_fit = 4;
      else if(gap_to_fit==5) gap_to_fit = 8;
      else if(gap_to_fit==6) gap_to_fit = 9;
      else if(gap_to_fit==7) gap_to_fit = 10;

      h_TDC_Gap_Fibers->Reset();
      vec_xx_lepton.clear();
      vec_yy_lepton.clear();
    
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

    
      count = 0;
      hit_count = 0;
    
      for(Int_t i=0; i<256; i++){
        if(TARGET_High_has_neighbours[i]){ 
          if(ADC_High_TARGET[i]>Angle_ADC_cut && has_TDC_hit[i]){        
            if(!k_stop_bar[i]){
              h_TDC_Gap_Fibers->Fill(Xloc[i],Yloc[i]);
              vec_xx_lepton.push_back(Xloc[i]);
              vec_yy_lepton.push_back(Yloc[i]);
            }
            
            count++;
                  
            if(IsIn(i,channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],   // Additional weight for fibers close to the edge if hit
                      channel[gap_to_fit-1][2], channel[gap_to_fit-1][3], 
                      channel[gap_to_fit-1][4], channel[gap_to_fit-1][5], 
                      channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){          
              if(!k_stop_bar[i]){     
                h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
                h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
                vec_xx_lepton.push_back(Xloc[i]);
                vec_xx_lepton.push_back(Xloc[i]);
                vec_yy_lepton.push_back(Yloc[i]);
                vec_yy_lepton.push_back(Yloc[i]);
              }
            }
          }   
          if(ADC_High_TARGET[i]<0 && ADC_Low_TARGET[i]>0){
            if(!k_stop_bar[i])  
              h_TDC_Gap_Fibers->Fill(Xloc[i],Yloc[i]);
              vec_xx_lepton.push_back(Xloc[i]);
              vec_yy_lepton.push_back(Yloc[i]);
              count++;
          }
        }
      }
    }

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(int i=0; i<12; i++){
      for(int j=0; j<3; j++){
        for(int k=0; k<2; k++){
          Gap[i][j][k] = 0;
        }
      }
    }
    a_fit_TDC_Gap_Fibers = 0.;   
    b_fit_TDC_Gap_Fibers = 0.;

    for(int g=0; g<12; g++){
      Gap[g][0][0] = TOF_Xloc[3*g];
      Gap[g][1][0] = TOF_Xloc[3*g+1];
      Gap[g][2][0] = TOF_Xloc[3*g+2];
    
      Gap[g][0][1] = TOF_Yloc[3*g];
      Gap[g][1][1] = TOF_Yloc[3*g+1];
      Gap[g][2][1] = TOF_Yloc[3*g+2];
    }

    // Add weight to center point of selected TOF1 (gap_to_fit_rotate)
    for(int i = 0; i<3; i++){
      h_TDC_Gap_Fibers->Fill(Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      vec_xx_lepton.push_back(Gap[gap_to_fit-1][1][0]);
      vec_yy_lepton.push_back(Gap[gap_to_fit-1][1][1]);
    }
  
    TGraph *gr_Leptons = new TGraph(vec_xx_lepton.size(),&vec_xx_lepton[0],&vec_yy_lepton[0]);
    gr_Leptons->Fit("pol1","QCM");
    a_gr_Leptons = gr_Leptons->GetFunction("pol1")->GetParameter(1);
    b_gr_Leptons = gr_Leptons->GetFunction("pol1")->GetParameter(0);
    a_fit_TDC_Gap_Fibers = a_gr_Leptons;
    b_fit_TDC_Gap_Fibers = b_gr_Leptons;

    // ReFill histogram with points closer than max_dist
    h_TDC_Gap_Fibers->Reset();
    vec_xx_lepton.clear();
    vec_yy_lepton.clear();
    gr_Leptons->Set(0);

    for(int i = 0; i<256; i++){   
      if(distance_to_line(Xloc[i],Yloc[i],a_fit_TDC_Gap_Fibers,b_fit_TDC_Gap_Fibers) <= max_dist && TARGET_High_has_neighbours[i] && !k_stop_bar[i]){   
        if(ADC_High_TARGET[i]>Angle_ADC_cut && has_TDC_hit[i]){
          h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
           vec_xx_lepton.push_back(Xloc[i]);
           vec_yy_lepton.push_back(Yloc[i]);
        if(IsIn(i,channel[gap_to_fit-1][0], channel[gap_to_fit-1][1],   
                channel[gap_to_fit-1][2], channel[gap_to_fit-1][3], 
                channel[gap_to_fit-1][4], channel[gap_to_fit-1][5], 
                channel[gap_to_fit-1][6], channel[gap_to_fit-1][7])){
            h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);     
            h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
            vec_xx_lepton.push_back(Xloc[i]);
            vec_yy_lepton.push_back(Yloc[i]);
          }
        }
        if(ADC_High_TARGET[i]<0 && ADC_Low_TARGET[i]>0){
          h_TDC_Gap_Fibers->Fill(Xloc[i], Yloc[i]);
          vec_xx_lepton.push_back(Xloc[i]);
          vec_yy_lepton.push_back(Yloc[i]);
        } 
      } 
    }
  
    if(vec_xx_lepton.size() != 0){
      TGraph *gr2_Leptons = new TGraph(vec_xx_lepton.size(),&vec_xx_lepton[0],&vec_yy_lepton[0]);
      gr2_Leptons->Fit("pol1","QCM");
      a_gr_Leptons = gr2_Leptons->GetFunction("pol1")->GetParameter(1);
      b_gr_Leptons = gr2_Leptons->GetFunction("pol1")->GetParameter(0);
    }
    else continue;

    a_fit_TDC_Gap_Fibers = a_gr_Leptons;
    b_fit_TDC_Gap_Fibers = b_gr_Leptons;

    for(int i=0; i<2; i++){
      x_int_TDC_Gap_Fibers[i] = 999.;     
      y_int_TDC_Gap_Fibers[i] = 999.;
    }

    x_int_TDC_Gap_Fibers[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);
    x_int_TDC_Gap_Fibers[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);   
    y_int_TDC_Gap_Fibers[0] = y1_int(x_int_TDC_Gap_Fibers[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers); 
    y_int_TDC_Gap_Fibers[1] = y2_int(x_int_TDC_Gap_Fibers[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);

    /// Selection of the Good Intersect Coordinates
    //////////////////////////////////////////////////////////////////////
    x_TDC_Gap_Fibers = 0.;       
    y_TDC_Gap_Fibers = 0.;

    for(int i=0; i<2; i++) dist1_TDC_Gap_Fibers[i] = 0.;

    for(int i=0; i<2; i++){
      dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
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
  
    if(dist1_TDC_Gap_Fibers[0] < dist1_TDC_Gap_Fibers[1])
    {
      x_TDC_Gap_Fibers = x_int_TDC_Gap_Fibers[0];
      y_TDC_Gap_Fibers = y_int_TDC_Gap_Fibers[0];
    }
    else if(dist1_TDC_Gap_Fibers[1] < dist1_TDC_Gap_Fibers[0])
    {
      x_TDC_Gap_Fibers = x_int_TDC_Gap_Fibers[1];
      y_TDC_Gap_Fibers = y_int_TDC_Gap_Fibers[1];
    }

    //////////////////////////////////////////////////////////////////////
    for(int i=0; i<3; i++){
      dist2_TDC_selected[i] = 0.;
    }

    dist2_TDC_selected_min = 1000.;
    selected_TDC_selected = 0;

    for(int ii=0; ii<3; ii++)
    {
      dist2_TDC_selected[ii] = distance(x_TDC_Gap_Fibers, y_TDC_Gap_Fibers, Gap[gap_to_fit-1][ii][0], Gap[gap_to_fit-1][ii][1]);

      if(dist2_TDC_selected[ii] <= dist2_TDC_selected_min)
      {
        dist2_TDC_selected_min = dist2_TDC_selected[ii];
        selected_TDC_selected = ii;
      }
    }

    //////////////////////////////////////////////////////////////////////
    for(int i = 0; i<3; i++){
      h_TDC_Gap_Fibers->Fill(Gap[gap_to_fit-1][selected_TDC_selected][0], Gap[gap_to_fit-1][selected_TDC_selected][1]);
      vec_xx_lepton.push_back(Gap[gap_to_fit-1][selected_TDC_selected][0]);
      vec_yy_lepton.push_back(Gap[gap_to_fit-1][selected_TDC_selected][1]);
    }
    //////////////////////////////////////////////////////////////////

    ParError = 999.99;

    TGraph *gr3_Leptons = new TGraph(vec_xx_lepton.size(),&vec_xx_lepton[0],&vec_yy_lepton[0]);
    gr3_Leptons->Fit("pol1","QCM");
    a_gr_Leptons = gr3_Leptons->GetFunction("pol1")->GetParameter(1);
    b_gr_Leptons = gr3_Leptons->GetFunction("pol1")->GetParameter(0);

    ParError = gr3_Leptons->GetFunction("pol1")->GetParError(1);
    ChiS = gr3_Leptons->GetFunction("pol1")->GetChisquare();
    ndf = gr3_Leptons->GetFunction("pol1")->GetNDF();
    prob = gr3_Leptons->GetFunction("pol1")->GetProb();

    a_fit_TDC_Gap_Fibers = a_gr_Leptons;
    b_fit_TDC_Gap_Fibers = b_gr_Leptons;

    x_int_TDC_Gap_Fibers[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);
    x_int_TDC_Gap_Fibers[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TOF1);   

    y_int_TDC_Gap_Fibers[0] = y1_int(x_int_TDC_Gap_Fibers[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers); 
    y_int_TDC_Gap_Fibers[1] = y2_int(x_int_TDC_Gap_Fibers[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);

    x_int_TARGET[0] = intersectx1(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TARGET);
    x_int_TARGET[1] = intersectx2(a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, R_TARGET);
    y_int_TARGET[0] = y1_int(x_int_TARGET[0], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);
    y_int_TARGET[1] = y2_int(x_int_TARGET[1], a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers);

    x_TDC_Gap_Fibers_intersect1 = 0.;          y_TDC_Gap_Fibers_intersect1 = 0.;
    x_TARGET_intersect = 0;                    y_TARGET_intersect = 0;

    for(int i=0; i<2; i++) dist1_TARGET_intersect[i] = 0.;

    for(int i=0; i<2; i++)
    {
      dist1_TDC_Gap_Fibers[i] = distance(x_int_TDC_Gap_Fibers[i], y_int_TDC_Gap_Fibers[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
      dist1_TARGET_intersect[i] = distance(x_int_TARGET[i], y_int_TARGET[i], Gap[gap_to_fit-1][1][0], Gap[gap_to_fit-1][1][1]);
    }    
 
    if(dist1_TDC_Gap_Fibers[0] < dist1_TDC_Gap_Fibers[1])
    {
      x_TDC_Gap_Fibers_intersect1 = x_int_TDC_Gap_Fibers[0];
      y_TDC_Gap_Fibers_intersect1 = y_int_TDC_Gap_Fibers[0];
    }
    else if(dist1_TDC_Gap_Fibers[1] < dist1_TDC_Gap_Fibers[0])
    {
      x_TDC_Gap_Fibers_intersect1 = x_int_TDC_Gap_Fibers[1];
      y_TDC_Gap_Fibers_intersect1 = y_int_TDC_Gap_Fibers[1];
    }

    if(dist1_TARGET_intersect[0] < dist1_TARGET_intersect[1])
    {
      x_TARGET_intersect = x_int_TARGET[0];
      y_TARGET_intersect = y_int_TARGET[0];
    }
    else if(dist1_TARGET_intersect[1] < dist1_TARGET_intersect[0])
    {
      x_TARGET_intersect = x_int_TARGET[1];
      y_TARGET_intersect = y_int_TARGET[1];
    }

  ///////////////////////////////////////////////////////////////////////////////

  dist_to_k_stop = 0.;
  //cout << "BERNARD:  " 
  //     << X_BAR << "  " << Y_BAR << "  "
  //     << a_fit_TDC_Gap_Fibers << "  " << b_fit_TDC_Gap_Fibers << endl;

  if(X_BAR != -10000){
    if((gap_to_fit_rotate==6 || gap_to_fit_rotate==12 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
      dist_to_k_stop = distance_to_line(Y_BAR,-X_BAR,a_fit_TDC_Gap_Fibers,b_fit_TDC_Gap_Fibers);
    }
    else{
      dist_to_k_stop = distance_to_line(X_BAR,Y_BAR,a_fit_TDC_Gap_Fibers,b_fit_TDC_Gap_Fibers);
    }
  }
  else{
    dist_to_k_stop = -10000;
  }

  //cout << "PIERRE :  " << dist_to_k_stop << endl;

  /// Angle Calculation
  ///////////////////////////////////////////////////////////////////////////////

    a_final_guide = 0.;
    alpha_guide = 0.;
    tanalpha_guide = 0.;
    angle_final_guide = 0.;
    
    a_final_guide = (y_TDC_Gap_Fibers_intersect1 - y_TARGET_intersect) / (x_TDC_Gap_Fibers_intersect1 - x_TARGET_intersect);
 
    tanalpha_guide = a_final_guide;
    alpha_guide = atan(tanalpha_guide);

    //Determination of Delta Phi
    Delta_phi = 999.99;  Delta_phi_deg = 999.99;
    //cout << "ALB: " << ivt << "  " << alpha_guide << "  " << ParError << "  " << a_final_guide << "  " << Delta_phi << endl;
    Delta_phi = sin(alpha_guide)*cos(alpha_guide)*(ParError/a_final_guide);
    Delta_phi_deg = (180/PI)*Delta_phi;

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
  
    if((x_TDC_Gap_Fibers_intersect1 - x_TARGET_intersect) < 0) angle_final_guide = 180. + alpha_guide * (180./PI);

    if((x_TDC_Gap_Fibers_intersect1 - x_TARGET_intersect) >= 0){
      if((y_TDC_Gap_Fibers_intersect1 - y_TARGET_intersect) >= 0)  angle_final_guide = alpha_guide * (180./PI);
      else angle_final_guide = alpha_guide * (180./PI) + 360.;
    }
    
    if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1) angle_final_guide += 90.;
    if(angle_final_guide > 360.) angle_final_guide -= 360.;

    ///////////////////////////////////////////////////////////////////////////////
    TLatex *x_sft;
    TLatex *y_sft;
    TLatex *x_tof1;
    TLatex *y_tof1;    
   
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

    length_in_target = 0.;
      

    if(X_BAR != -10000){
      x_tof1_intersect_1 = 0.;
      y_tof1_intersect_1 = 0.;
      x_tof1_intersect_2 = 0.;
      y_tof1_intersect_2 = 0.;
      x_tof1_intersect = 0.;
      y_tof1_intersect = 0.;

      alpha = angle_final_guide - 90.0;

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
      //cout << "VALENTIN :  " << x_TDC_Gap_Fibers_intersect1 << "  " << y_TDC_Gap_Fibers_intersect1 << "  "
      //     << dist_to_k_stop << "  " << alpha << endl;

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
         
         //cout << "PAULINE :  " << X_BAR << "  "
         //<< Y_BAR << "  "
         //<< x_tof1_intersect << "  "
         //<< y_tof1_intersect << "  "
         //<< endl; 

        }
      }


    length_in_target = distance(X_BAR,Y_BAR,x_tof1_intersect,y_tof1_intersect); //CORRECT
    }

    Z_selected.clear();
     
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    i_kaon_bar = kaon_fiber(X_BAR,Y_BAR);

    Z_TOF1[gap_to_fit_rotate-1] = 
    Slope_correction_TOF1[gap_to_fit_rotate-1]*((0.025*150*0.5*(TDC_TOF1U[gap_to_fit_rotate-1] - TDC_TOF1D[gap_to_fit_rotate-1])) + Offset_TOF1[gap_to_fit_rotate-1]);

    Yth_TOF1 = a_fit_TDC_Gap_Fibers*Gap[gap_to_fit-1][selected_TDC_selected][0] + b_fit_TDC_Gap_Fibers;
    New_ChiS = ChiS - 3*((Gap[gap_to_fit-1][selected_TDC_selected][1] - Yth_TOF1)*(Gap[gap_to_fit-1][selected_TDC_selected][1] - Yth_TOF1));

    if(Switch_Display==0){ 
      if(ndf!=3){
        if(!isnan(Delta_phi_deg)){
          cout << fixed;
          cout << setw(4) << Run_Number << "  ";
          cout << setw(7) << ivt << "  ";
          cout << setw(2) << gap_to_fit_rotate << "  ";
          cout << setw(2) << selected_TOF2 << "  ";
          cout << setw(7) << setprecision(3) << angle_final_guide << "  ";
          cout << setw(5) << setprecision(2) << Delta_phi_deg << "  ";
          //cout << setw(8) << setprecision(2) << ChiS << "  ";
          //cout << setw(3) << setprecision(0) << ndf-2 << "  ";
          //cout << setw(7) << setprecision(2) << ChiS/(ndf-2) << "  ";
          cout << setw(8) << setprecision(2) << New_ChiS << "  ";
          cout << setw(3) << setprecision(0) << ndf-3 << "  ";
          cout << setw(7) << setprecision(2) << New_ChiS/(ndf-3) << "  ";
          cout << setw(7) << setprecision(2) << (New_ChiS/(ndf-3))*(pow(cos(angle_final_guide*M_PI/180.),2)) << "  ";
          cout << setw(6) << setprecision(2) << X_BAR << "  ";
          cout << setw(6) << setprecision(2) << Y_BAR << "  ";
          cout << setw(3) << i_kaon_bar << "  ";
          cout << setw(2) << K_stop_count << "  ";
          cout << setw(7) << setprecision(2) << x_TARGET_intersect << "  ";
          cout << setw(7) << setprecision(2) << y_TARGET_intersect << "  ";
          cout << setw(3) << TDC_Ck_counter << "  ";
          cout << setw(3) << TDC_Cpi_counter << "  ";
          cout << setw(8) << setprecision(3) << length_in_target << "  ";
          cout << setw(8) << setprecision(3) << C2X_centroid << "  ";
          cout << setw(6) << setprecision(1) << TDC_diff << "  ";
          cout << endl;
        }
        else nan_evts.push_back(ivt);
      }
    }

    if(Switch_Display==1 || Switch_Display==2){
      if(ndf!= 3){
        if(!isnan(Delta_phi_deg)){
          fout << fixed;
          fout << setw(4) << Run_Number << "  ";
          fout << setw(7) << ivt << "  ";
          fout << setw(2) << gap_to_fit_rotate << "  ";
          fout << setw(2) << selected_TOF2 << "  ";
          fout << setw(7) << setprecision(3) << angle_final_guide << "  ";
          fout << setw(5) << setprecision(2) << Delta_phi_deg << "  ";
          //fout << setw(8) << setprecision(2) << ChiS << "  ";
          //fout << setw(3) << setprecision(0) << ndf-2 << "  ";
          //fout << setw(7) << setprecision(2) << ChiS/(ndf-2) << "  ";
          fout << setw(8) << setprecision(2) << New_ChiS << "  ";
          fout << setw(3) << setprecision(0) << ndf-3 << "  ";
          fout << setw(7) << setprecision(2) << New_ChiS/(ndf-3) << "  ";
          fout << setw(7) << setprecision(2) << (New_ChiS/(ndf-3))*(pow(cos(angle_final_guide*M_PI/180.),2)) << "  ";
          fout << setw(6) << setprecision(2) << X_BAR << "  ";
          fout << setw(6) << setprecision(2) << Y_BAR << "  ";
          fout << setw(3) << i_kaon_bar << "  ";
          fout << setw(2) << K_stop_count << "  ";
          fout << setw(7) << setprecision(2) << x_TARGET_intersect << "  ";
          fout << setw(7) << setprecision(2) << y_TARGET_intersect << "  ";
          fout << setw(3) << TDC_Ck_counter << "  ";
          fout << setw(3) << TDC_Cpi_counter << "  ";
          fout << setw(8) << setprecision(3) << length_in_target << "  ";
          fout << setw(8) << setprecision(3) << C2X_centroid << "  ";
          fout << setw(6) << setprecision(1) << TDC_diff << "  ";
          fout << endl;
        }
        else nan_evts.push_back(ivt);
      }
    }
    
    h_TDC_Gap_Fibers->Reset();
    gr_Leptons->Set(0);

  } // ENDLOOP OVER EVENTS! (BATCH LOOP)
  
  cout << endl;
  cout << endl;
  
  cout << nan_evts.size() << " NaN events removed !" << endl;
  cout << endl;

  for(unsigned int i=0; i<nan_evts.size(); i++){
    fout2 << i+1 << "  " << nan_evts[i] << endl;
  }
 
  cout << endl;

  fout.close();
  fout2.close();

  if(!gROOT->IsBatch()) return;

  delete h_TDC_Gap_Fibers;

} // End void

int kaon_fiber(double x_bar, double y_bar){

  double dist_min = 999.;
  int i_min = 999;

  for(int i=0; i<256; i++){
    if(sqrt(pow(abs(Xloc[i] - x_bar),2) + pow(abs(Yloc[i] - y_bar),2)) < dist_min){
      dist_min = sqrt(pow(abs(Xloc[i] - x_bar),2) + pow(abs(Yloc[i] - y_bar),2));
      i_min = i;
    }
  }

  return i_min;  
}




