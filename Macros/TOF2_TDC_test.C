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
#include "Cuts_and_Windows.h"
#include "MWPC_Thr.h"
#endif

#include "intersect.cxx"

// view_type: 1 =  AO or AI or BO or BI
//            2 = (AO and AI) or (BO and BI)
//            3 = Exactly one of AO,AI,BO,BI
//            4 = Atleast one of (AO or AI)   AND   atleast one (BO or BI)

void TOF2_TDC_test(Int_t Run_Number = 5, Int_t first_event=0, Int_t last_event=10){
	
	
  Int_t ADC_tof2AO[12];             Int_t ADC_TOF2AO[12];
  Int_t ADC_tof2BO[12];             Int_t ADC_TOF2BO[12];
  Int_t ADC_tof2AI[12];             Int_t ADC_TOF2AI[12];
  Int_t ADC_tof2BI[12];             Int_t ADC_TOF2BI[12];
	
  Int_t TDC_tof2AO[12];             Int_t TDC_TOF2AO[12];
  Int_t TDC_tof2BO[12];             Int_t TDC_TOF2BO[12];
  Int_t TDC_tof2AI[12];             Int_t TDC_TOF2AI[12];
  Int_t TDC_tof2BI[12];             Int_t TDC_TOF2BI[12];	
	
	
	
  char path_input[200];                   
  sprintf(path_input,"%s",path_merged);    
       	
  char Name_finput[200];
  sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);
  
  TChain *fChain= new TChain("Tree");   
  fChain->Add(Name_finput);   
  fChain->SetMakeClass(1);  
  
  fChain->SetBranchAddress("ADC_TOF2AO",ADC_tof2AO);
  fChain->SetBranchAddress("ADC_TOF2AI",ADC_tof2AI);
  fChain->SetBranchAddress("ADC_TOF2BO",ADC_tof2BO);
  fChain->SetBranchAddress("ADC_TOF2BI",ADC_tof2BI);    
  fChain->SetBranchAddress("TDC_TOF2AO",TDC_tof2AO);
  fChain->SetBranchAddress("TDC_TOF2AI",TDC_tof2AI);
  fChain->SetBranchAddress("TDC_TOF2BO",TDC_tof2BO);
  fChain->SetBranchAddress("TDC_TOF2BI",TDC_tof2BI);   
  
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////    
  
  
  const Double_t w = 800;
  const Double_t h = 400;    
    
  const float bin_edges[8] = {0.0,0.5,1.5,2.5,3.5,4.5,5.5,6.5};	
  const float bin_edges_bool[3] = {0.0,0.5,1.5};
  TH1F *TDC_TOF2_1 = new TH1F("TDC_TOF2_1", "TDC_TOF2_1", 7, bin_edges);
  TH1F *TDC_TOF2_2 = new TH1F("TDC_TOF2_2", "TDC_TOF2_2", 7, bin_edges);
  TH1F *TDC_TOF2_3 = new TH1F("TDC_TOF2_3", "TDC_TOF2_3", 2, bin_edges_bool);
  TH1F *TDC_TOF2_4 = new TH1F("TDC_TOF2_4", "TDC_TOF2_4", 2, bin_edges_bool);
  
  
  
    
  const char *view_type_1 = "AO or  AI or BO or BI";
  const char *view_type_2 = "(AO and AI) or (BO and BI)";
  const char *view_type_3 = "Exactly one of AO,AI,BO,BI";
  const char *view_type_4 = "Atleast one of (AO or AI)   AND   atleast one (BO or BI)";
    
  int Good_TOF2_counter_1 = 0;
  int Good_TOF2_counter_2 = 0;
  int Good_TOF2_counter_3 = 0;
  int Good_TOF2_counter_4 = 0;  
	
  for(int ievt = first_event; ievt < last_event; ievt++){	
	  
	// Output every 10000 events  
    if(ievt%10000 == 0)
      cout << "Event: " << ievt << endl;


    fChain->GetEntry(ievt);
    
    int Good_TOF2_TDC_1[12] = {false};
    int Good_TOF2_TDC_2[12] = {false};
    int Good_TOF2_TDC_3[12] = {false};
    int Good_TOF2_TDC_4[12] = {false};    
	
	
	
	for (Int_t j_TDCTOF=0; j_TDCTOF<12; j_TDCTOF++) {
	  TDC_TOF2AO[j_TDCTOF] = TDC_tof2AO[j_TDCTOF];
	  TDC_TOF2BO[j_TDCTOF] = TDC_tof2BO[j_TDCTOF];
	  TDC_TOF2AI[j_TDCTOF] = TDC_tof2AI[j_TDCTOF];
	  TDC_TOF2BI[j_TDCTOF] = TDC_tof2BI[j_TDCTOF];
	}  
	
    
    
    
    // Count number of Good TOF2s in one event
    //////////////////////////////////////////
	// AO or  AI or BO or BI	
    for(int i=0; i<12; i++){
	    if((TDC_TOF2AO[i]>=TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<=TOF2AO_TDC_max[i]) ||
	       (TDC_TOF2AI[i]>=TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<=TOF2AI_TDC_max[i]) ||
	       (TDC_TOF2BO[i]>=TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<=TOF2BO_TDC_max[i]) ||
	       (TDC_TOF2BI[i]>=TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<=TOF2BI_TDC_max[i]))    Good_TOF2_TDC_1[i] = true;
	}
	// (AO and AI) or (BO and BI)
    for(int i=0; i<12; i++){
	    if(((TDC_TOF2AO[i]>=TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<=TOF2AO_TDC_max[i]) &&
	        (TDC_TOF2AI[i]>=TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<=TOF2AI_TDC_max[i])) ||
	       ((TDC_TOF2BO[i]>=TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<=TOF2BO_TDC_max[i]) &&
	        (TDC_TOF2BI[i]>=TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<=TOF2BI_TDC_max[i])))  Good_TOF2_TDC_2[i] = true;
    }
    // Exactly one of AO,AI,BO,BI
    for(int i=0; i<12; i++){
      if(((TDC_TOF2AO[i]>=TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<=TOF2AO_TDC_max[i]) &&
	     !(TDC_TOF2AI[i]>=TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<=TOF2AI_TDC_max[i]) &&
	     !(TDC_TOF2BO[i]>=TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<=TOF2BO_TDC_max[i]) &&
	     !(TDC_TOF2BI[i]>=TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<=TOF2BI_TDC_max[i])) ||
	     
	     (!(TDC_TOF2AO[i]>=TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<=TOF2AO_TDC_max[i]) &&
	       (TDC_TOF2AI[i]>=TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<=TOF2AI_TDC_max[i]) &&
	      !(TDC_TOF2BO[i]>=TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<=TOF2BO_TDC_max[i]) &&
	      !(TDC_TOF2BI[i]>=TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<=TOF2BI_TDC_max[i])) ||
	      
	     (!(TDC_TOF2AO[i]>=TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<=TOF2AO_TDC_max[i]) &&
	      !(TDC_TOF2AI[i]>=TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<=TOF2AI_TDC_max[i]) &&
	       (TDC_TOF2BO[i]>=TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<=TOF2BO_TDC_max[i]) &&
	      !(TDC_TOF2BI[i]>=TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<=TOF2BI_TDC_max[i])) ||
	  
	     (!(TDC_TOF2AO[i]>=TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<=TOF2AO_TDC_max[i]) &&
	      !(TDC_TOF2AI[i]>=TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<=TOF2AI_TDC_max[i]) &&
	      !(TDC_TOF2BO[i]>=TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<=TOF2BO_TDC_max[i]) &&
	       (TDC_TOF2BI[i]>=TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<=TOF2BI_TDC_max[i])))  Good_TOF2_TDC_3[i] = true;
	}
	// Atleast one of (AO or AI)   AND   atleast one (BO or BI)
	for(int i=0; i<12; i++){
	  if(((TDC_TOF2AO[i]>=TOF2AO_TDC_min[i] && TDC_TOF2AO[i]<=TOF2AO_TDC_max[i]) ||
	      (TDC_TOF2AI[i]>=TOF2AI_TDC_min[i] && TDC_TOF2AI[i]<=TOF2AI_TDC_max[i])) &&
         ((TDC_TOF2BO[i]>=TOF2BO_TDC_min[i] && TDC_TOF2BO[i]<=TOF2BO_TDC_max[i])  ||
	      (TDC_TOF2BI[i]>=TOF2BI_TDC_min[i] && TDC_TOF2BI[i]<=TOF2BI_TDC_max[i])))   Good_TOF2_TDC_4[i] = true;	      	
	}
	
	
	
	// Find total number of Good TOF2s in one event
	///////////////////////////////////////////////
    Good_TOF2_counter_1 = 0;
    Good_TOF2_counter_2 = 0;
    Good_TOF2_counter_3 = 0;
    Good_TOF2_counter_4 = 0;    
	 
	for(int i=0; i<12; i++){
	   if(Good_TOF2_TDC_1[i])
	     Good_TOF2_counter_1++;
	   if(Good_TOF2_TDC_2[i])
	     Good_TOF2_counter_2++;
	   if(Good_TOF2_TDC_3[i])
	     Good_TOF2_counter_3++;
	   if(Good_TOF2_TDC_4[i])
	     Good_TOF2_counter_4++;	     
	}
    
    
    //Set name of histogram based on view type
    //////////////////////////////////////////
    TDC_TOF2_1->SetNameTitle("TDC_TOF2_1", view_type_1);
    TDC_TOF2_2->SetNameTitle("TDC_TOF2_2", view_type_2);
    TDC_TOF2_3->SetNameTitle("TDC_TOF2_3", view_type_3);
    TDC_TOF2_4->SetNameTitle("TDC_TOF2_4", view_type_4);    
            
      
	TDC_TOF2_1->Fill(double(Good_TOF2_counter_1) + 0.4);
	TDC_TOF2_2->Fill(double(Good_TOF2_counter_2) + 0.4);
	TDC_TOF2_3->Fill(double(Good_TOF2_counter_3) + 0.4);
	TDC_TOF2_4->Fill(double(Good_TOF2_counter_4) + 0.4);	
	
  }
  
  // SETTING CANVAS
  
  // 2 Seperate canvases
  //TCanvas *c1 = new TCanvas("c","c", w, h);
  //c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
  //c1->Divide(2,1);
  //TCanvas *c2 = new TCanvas("c2","c2", w, h);
  //c2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
  //c2->Divide(2,1);  
  
  // 1 canvas for all histograms
  TCanvas *c3 = new TCanvas("c3","c3", w, w);
  c3->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
  c3->Divide(2,2);   
  
  
  
  
  // DRAW CANVAS AND HSTOGRAMS
  
  // 2 Seperate canvases
  //c1->cd(1);
  //TDC_TOF2_1->Draw();
  //c1->cd(2);
  //TDC_TOF2_2->Draw();
  
  //c2->cd(1);
  //TDC_TOF2_3->Draw();
  //c2->cd(2);
  //TDC_TOF2_4->Draw();  
  
  
  // All Histograms on 1 canvas
  
  TDC_TOF2_3->SetNdivisions(2);
  TDC_TOF2_4->SetNdivisions(2);  
  
  
  c3->cd(1);
  TDC_TOF2_1->Draw();
  c3->cd(2);
  TDC_TOF2_2->Draw();
  c3->cd(3);
  TDC_TOF2_3->Draw();
  c3->cd(4);
  TDC_TOF2_4->Draw();    


return;

}
