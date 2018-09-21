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
#include "TMarker.h"
#include "ANAPATH.h"
#include "mappings.h"
#include "Cuts_and_Windows.h"
#include "ADC_TARGET_Pedestals.h"
#include "CommonParameters.h"
#endif


void TDC_Ck_Avg(Int_t run_number=5){ 

	gStyle->SetOptStat(111111111);

	// Ck
	Int_t TDC_ck[14][16];
	Int_t TDC_ck_selected[14] = {0};
	//double TDC_ck_sum = 0.;

	//TARGET
    Int_t adc_high_target[256];       Int_t ADC_High_TARGET[256];    
	Int_t tdc_le_target[256][16];     Int_t TDC_LE_TARGET[256];     
	Int_t adc_low_target[256];        Int_t ADC_Low_TARGET[256];
	Int_t LG_TARGET_ADC_Thr[256] = {0};
	int T_limit = 3;

	char path_input[200];                   //char file_mapping[200];
	sprintf(path_input,"%s",path_merged);          

	char Name_finput[200];
	sprintf(Name_finput,"%s/Run%dMS.root",path_input, run_number);

	cout << "   " << endl;
	cout << Name_finput << endl;

	TChain *fChain= new TChain("Tree");		
	fChain->Add(Name_finput);		
	fChain->SetMakeClass(1);							

	// Ck
	fChain->SetBranchAddress("TDC_Ck",TDC_ck);
  
    // TARGET
    fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);
  	fChain->SetBranchAddress("ADC_Low_TARGET",adc_low_target);
    fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);

	TH1D *h_TDC_Ck_Avg;     char Title_TDC_Ck_Avg[100];		char Name_TDC_Ck_Avg[100];
	TH1D *h_TDC_Average;    char Title_TDC_Average[100];	char Name_TDC_Average[100];
	TH1D *h_TDC_Diff;     	char Title_TDC_Diff[100];		char Name_TDC_Diff[100];
	vector <double> vec_Ck;	
	vec_Ck.clear();

	Int_t TDC_ck_sum = 0;       double TDC_ck_avg = 0.;     double TDC_ck_sigma = 0.;
    double TDC_ck_sigma2 = 0.;
    Int_t TDC_ck_sum2 = 0;    double TDC_ck_avg2=0.;    int TDC_ck_counter = 0;
    int TDC_Ck_counter =0;

  	Int_t TDC_average = -1;
  	int max_index_all[256];
  	int max_ADC_all[256];
  	int max_index_flag;
   	int TDC_LG_max = -1;
  	int TDC_LG_max2 = -1;
  	int index_max1=0;
  	int index_max2=0;
 
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
  	sprintf(ParsTarg1,"TimeWalk%d.dat",run_number);

  	if(!ifstream(ParsTarg1)) NoFile = true;
 
  	///////////////////////////////////////////////////////////////////////////////////
  	cout << "   " << endl;
  	cout << "   " << endl;
  	cout << "************************************************************************************************************" << endl;

  	//cout << "   " << endl;
  	cout << "File opened:  " << Name_finput << endl;
  	if(!NoFile) cout << "Time Walk Parameter File:  " << ParsTarg1 << endl;
  	if(NoFile) cout << "Time Walk Parameter File:  " << "NOT FOUND !" << endl;



	/////////////////////////////////////////////////////////////
	// Setting Histogram Names and Titles
	
	sprintf(Title_TDC_Ck_Avg,"Ck | TDC Average (Run %d)", run_number);
	sprintf(Name_TDC_Ck_Avg,"TDC Ck Average");

	sprintf(Title_TDC_Average,"Ck | Average (Run %d)", run_number);
	sprintf(Name_TDC_Average,"TDC Average");

	sprintf(Title_TDC_Diff,"Ck | Diff (Run %d)", run_number);
	sprintf(Name_TDC_Diff,"TDC Diff");


	h_TDC_Ck_Avg = new TH1D(Name_TDC_Ck_Avg, Title_TDC_Ck_Avg, 500, 1300, 1800);
	h_TDC_Average = new TH1D(Name_TDC_Average, Title_TDC_Average, 500, 600, 1100);
	h_TDC_Diff = new TH1D(Name_TDC_Diff, Title_TDC_Diff, 500, 450, 950);

	// READ ALL ENTRIES 
	Int_t nentries = (Int_t)fChain->GetEntries();
	cout <<  "Total Number of Entries :     " <<  nentries << endl;

	cout << "  " << endl;
	cout << "***************  Processing....  ***************" << endl;
	cout << "  " << endl;

	//for (Int_t ivt=0; ivt<nentries; ivt++) {
	for (Int_t ivt=0; ivt<10000; ivt++) {

		fChain->GetEntry(ivt);
		if(ivt%10000==1) cout << "**** " << ivt << " events done" << endl;

  		for(int i = 0; i<256; i++){  
  		    LG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Ped_LG[i]) + TARGET_ADC_Thr_LG_Offset;    
  			TDC_LE_TARGET[i]=tdc_le_target[i][0];
  			ADC_Low_TARGET[i]=adc_low_target[i]-LG_TARGET_ADC_Thr[i];
  		}
 

 	 	////////////////////////  Dave's Time Walk Correction ///////////////////////
  		for(int i=0; i<256; i++){
  			for(int j=0; j<3; j++){
  				par_in[i][j]=0.;
  				par_err[i][j]=0.;
  			}
  		}
  		ADCmax = 3450;
  		Yfi = 0.;
  		Ani = 0.;
  		Yni = 0.;
  		Tpi = 0.;
  		NoFile = false;


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

        	//if((Ani >= TARGET_ADC_Ped_HG[ii]) && (Ani < ADCmax)){ // SEB
        	if((Ani >= TARGET_ADC_Thr_HG_Offset) && (Ani < ADCmax)){ // SEB
          		Yni = par_in[ii][0] - par_in[ii][1]/sqrt(Ani - par_in[ii][2]);
          		Tpi = Yfi - Yni;
          		for(int jj=0; jj<16; jj++){
            		if(tdc_le_target[ii][jj]>0) tdc_le_target[ii][jj] = tdc_le_target[ii][jj] + Tpi;
          		}
        	}
      	}

    	/////////////////////////////////////////////////////////////////////////////




		//////////////////////////////////////////////////////////////////////////////////////////////
  		TDC_average = -1;
  
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


  		TDC_LG_max = -1;
  		TDC_LG_max2 = -1;
  		index_max1=0;
  		index_max2=0;
  
  		for(int i = 0; i<256; i++){
  			for(int j=0; j<6; j++){  
    			if(TDC_LG_max == -1){ 
     				if(tdc_le_target[max_index_all[i]][j] >= TDC_Thr_min && tdc_le_target[max_index_all[i]][j]  <= TDC_Thr_max)
        				TDC_LG_max = tdc_le_target[max_index_all[i]][j];
        				index_max1 = max_index_all[i];
        				//cout << "TEST1 : " << tdc_le_target[max_index_all[i]][j] << endl;
    			}
    			else if(TDC_LG_max2 == -1){
      					if(tdc_le_target[max_index_all[i]][j] >= TDC_Thr_min && tdc_le_target[max_index_all[i]][j] <= TDC_Thr_max)
        					TDC_LG_max2 = tdc_le_target[max_index_all[i]][j];
       						index_max2 = max_index_all[i];
        				    //cout << "TEST2 : " << index_max2 << endl;
      			}
  			}     
  		} 
  
  //cout << "TEST: " << 
  //cout << "TEST : " << 	TDC_LG_max << "  " << TDC_LG_max2 << "  " << T_limit << endl;

  if(abs(TDC_LG_max - TDC_LG_max2) <= T_limit) TDC_average = (TDC_LG_max + TDC_LG_max2)/2;
  else if(TDC_LG_max > TDC_LG_max2) TDC_average = TDC_LG_max;
  else if(TDC_LG_max < TDC_LG_max2) TDC_average = TDC_LG_max2;

  //cout << "PAULINE : " << ivt << "   " << TDC_average << endl;

////////////////////////////////////////////////////////////////////////////////////////////////////////////

 
		for(int j=0; j<14; j++){
			for(int k=0; k<8; k++){
				if(TDC_ck[j][k]>=TDC_Ck_min && TDC_ck[j][k]<=TDC_Ck_max) TDC_ck_selected[j]=TDC_ck[j][k];
			}
			if(TDC_ck_selected[j]>0) vec_Ck.push_back(TDC_ck_selected[j]);
			TDC_ck_selected[j]=0;
		}

		for(unsigned j=0; j<vec_Ck.size(); j++) TDC_ck_sum += vec_Ck[j];

	    if(vec_Ck.size()>0) TDC_ck_avg = double(TDC_ck_sum)/double(vec_Ck.size()); 
	    else TDC_ck_avg = -1;

	    for(unsigned int i=0; i<vec_Ck.size(); i++) TDC_ck_sigma2 += pow((vec_Ck[i]-TDC_ck_avg),2);

        if(vec_Ck.size()>0) TDC_ck_sigma = sqrt(TDC_ck_sigma2/vec_Ck.size());
        else TDC_ck_sigma = -1;

	    for(unsigned int i=0; i<vec_Ck.size(); i++){
        	if(abs(vec_Ck[i]-TDC_ck_avg) < 1.5*TDC_ck_sigma){
      			TDC_ck_sum2 += vec_Ck[i];
      			TDC_ck_counter++;
    		}
  		}

  		if(TDC_ck_counter>0) TDC_ck_avg2 = TDC_ck_sum2/TDC_ck_counter;
  		else TDC_ck_avg2 = -1;

  		//cout << "ALBANE : " << ivt << "   " << TDC_ck_avg2 << endl;	

		h_TDC_Ck_Avg->Fill(TDC_ck_avg2);
		h_TDC_Average->Fill(TDC_average);
		h_TDC_Diff->Fill(TDC_ck_avg2-TDC_average);

		//cout << "TEST : " << ivt << "  " << TDC_ck_avg2-TDC_average << endl;
 
		vec_Ck.clear();
		TDC_ck_sum = 0;
		TDC_ck_sum2 = 0;
		TDC_ck_avg = 0;
		TDC_ck_avg2 = 0;
		TDC_ck_sigma = 0;
		TDC_ck_sigma2 = 0;
		TDC_ck_counter =0;



	} // EndLoop over Events

	
	char Name_Can_TDC_Ck_Avg[100];
	char Title_Can_TDC_Ck_Avg[100];

	sprintf(Name_Can_TDC_Ck_Avg,"TDC Ck Average -- Run %d",run_number);
	sprintf(Title_Can_TDC_Ck_Avg,"TDC Ck Average");

	TCanvas *c1;
	c1 = new TCanvas(Name_Can_TDC_Ck_Avg, Title_Can_TDC_Ck_Avg, 600, 1200);
	c1->Divide(1,3);
	c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
	
	c1->cd(1);
	h_TDC_Ck_Avg->Draw();

	c1->cd(2);
	h_TDC_Average->Draw();

	c1->cd(3);
	h_TDC_Diff->Draw();


	cout << "" << endl;
    cout << "" << endl;

	return;
}

