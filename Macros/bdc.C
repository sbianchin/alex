#ifndef __CINT__
#include <stdio.h>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <string.h>
#include <algorithm>
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

using namespace std;

void bdc(int run_number = 3994){

	gStyle->SetOptStat(1111111111);
	TH1::AddDirectory(kFALSE);

	char path_input[200];
	sprintf(path_input,"%s",path_merged);

	char Name_finput[200];
	sprintf(Name_finput,"%s/Run%dMS.root",path_input,run_number);

	cout << "   " << endl;
	cout << Name_finput << endl;

	TChain *fChain = new TChain("Tree");
	fChain->Add(Name_finput);
	fChain->SetMakeClass(1);

	// Trigger
	int tdc_trigger[2][16];
	fChain->SetBranchAddress("TDC_Trig",tdc_trigger);
	TH1D *h_TDC_Trigger;  char Title_TDC_Trigger[100];	char Name_TDC_Trigger[100];
	TH1D *h_TDC_Trig_VT;  char Title_TDC_Trig_VT[100];	char Name_TDC_Trig_VT[100];
	sprintf(Title_TDC_Trigger,"TDC Trigger");
	sprintf(Name_TDC_Trigger,"h_TDC_Trigger (Run %d)",run_number);
	sprintf(Title_TDC_Trig_VT,"TDC Trigger VT84");
	sprintf(Name_TDC_Trig_VT,"h_TDC_Trig_VT (Run %d)",run_number);
	h_TDC_Trigger = new TH1D(Name_TDC_Trigger, Title_TDC_Trigger, 75, 1800, 2100);
	h_TDC_Trig_VT = new TH1D(Name_TDC_Trig_VT, Title_TDC_Trig_VT, 75, 1800, 2100);

	// VT48
	int tdc_vt48[256][16];
	fChain->SetBranchAddress("VT48_TDC",tdc_vt48);

	// Beam Defining Counter
	int tdc_bdc[16] = {-1};
	vector<int> vec_tdc_bdc;  vec_tdc_bdc.clear();
	int tdc_bdc_selected = -1;
	vector<int> vec_tdc_b0_6;  vec_tdc_b0_6.clear();
	vector<int> vec_tdc_b0_7;  vec_tdc_b0_7.clear();
	int tdc_b0_6 = -1; int tdc_b0_7 = -1; int tdc_b0 = -1;
	TH1D *h_TDC_BDC_Raw;   char Title_TDC_BDC_Raw[100];	    char Name_TDC_BDC_Raw[100];
	TH1D *h_TDC_BDC_Diff;  char Title_TDC_BDC_Diff[100];	char Name_TDC_BDC_Diff[100];
	TH1D *h_TDC_BDC_NULL;  char Title_TDC_BDC_NULL[100];	char Name_TDC_BDC_NULL[100];
	TH1D *h_TDC_B0;  char Title_TDC_B0[100];	char Name_TDC_B0[100];
	TH1D *h_TDC_B0_NULL;  char Title_TDC_B0_NULL[100];	char Name_TDC_B0_NULL[100];
	sprintf(Title_TDC_BDC_Raw,"Beam Defining Counter");
	sprintf(Name_TDC_BDC_Raw,"h_TDC_BDC_Raw (Run %d)", run_number);
	sprintf(Title_TDC_BDC_Diff,"Beam Defining Counter Corrected");
	sprintf(Name_TDC_BDC_Diff,"h_TDC_BDC_Diff (Run %d)", run_number);
	sprintf(Title_TDC_BDC_NULL,"Beam Defining Counter NULL");
	sprintf(Name_TDC_BDC_NULL,"h_TDC_BDC_NULL (Run %d)", run_number);
	sprintf(Title_TDC_B0,"B0 Counter");
	sprintf(Name_TDC_B0,"h_TDC_B0 (Run %d)", run_number);
	sprintf(Title_TDC_B0_NULL,"B0 Counter NULL");
	sprintf(Name_TDC_B0_NULL,"h_TDC_B0_NULL (Run %d)", run_number);
	h_TDC_BDC_Raw = new TH1D(Name_TDC_BDC_Raw, Title_TDC_BDC_Raw, 75, 1450, 1750);
	h_TDC_BDC_Diff = new TH1D(Name_TDC_BDC_Diff, Title_TDC_BDC_Diff, 75, 200, 500);
	h_TDC_BDC_NULL = new TH1D(Name_TDC_BDC_NULL, Title_TDC_BDC_NULL, 75, -10, 10);
	h_TDC_B0 = new TH1D(Name_TDC_B0, Title_TDC_B0, 250, 1000, 2000);
	h_TDC_B0_NULL = new TH1D(Name_TDC_B0_NULL, Title_TDC_B0_NULL, 75, -10, 10);

	// Ck
	int tdc_ck[14][16] = {-1};
  	int TDC_ck_sum = 0;       double TDC_ck_avg = 0.;     double TDC_ck_sigma = 0.;
  	int TDC_cpi_sum = 0;      double TDC_cpi_avg = 0.;    double TDC_cpi_sigma = 0.;
  	double TDC_ck_sigma2 = 0.;
  	double TDC_cpi_sigma2 = 0.;
  	int TDC_ck_sum2 = 0;    double TDC_ck_avg2=0.;    int TDC_ck_counter = 0;
  	int TDC_cpi_sum2 = 0;   double TDC_cpi_avg2=0.;   int TDC_cpi_counter = 0;   
	double ck_mean = 0.;
	double tdc_ck_corr = 0.;
	vector<int> vec_tdc_ck;      vec_tdc_ck.clear();
	vector<int> vec_tdc_ck_col;  vec_tdc_ck_col.clear();
	fChain->SetBranchAddress("TDC_Ck",tdc_ck);
	TH1D *h_TDC_CK;  char Title_TDC_CK[100];	char Name_TDC_CK[100];
	TH1D *h_TDC_CK_corr;  char Title_TDC_CK_corr[100];	char Name_TDC_CK_corr[100];
	TH1D *h_TDC_CK_BDC_Diff;  char Title_TDC_CK_BDC_Diff[100];	char Name_TDC_CK_BDC_Diff[100];
	TH1D *h_TDC_CK_BDC_Diff2;  char Title_TDC_CK_BDC_Diff2[100];	char Name_TDC_CK_BDC_Diff2[100];
	sprintf(Title_TDC_CK,"TDC Mean Ck");
	sprintf(Name_TDC_CK,"h_TDC_CK (Run %d)", run_number);
	sprintf(Title_TDC_CK_corr,"TDC Mean Ck Corrected");
	sprintf(Name_TDC_CK_corr,"h_TDC_CK_corr (Run %d)",run_number);
	sprintf(Title_TDC_CK_BDC_Diff,"TDC Mean Ck - TDC BDC");
	sprintf(Name_TDC_CK_BDC_Diff,"h_TDC_CK_BDC_Diff (Run %d)", run_number);
	sprintf(Title_TDC_CK_BDC_Diff2,"TDC Mean Ck corr. - TDC BDC corr.");
	sprintf(Name_TDC_CK_BDC_Diff2,"Name_TDC_CK_BDC_Diff2 (Run %d)", run_number);
	h_TDC_CK = new TH1D(Name_TDC_CK, Title_TDC_CK, 75, 1400, 1700);
	h_TDC_CK_corr = new TH1D(Name_TDC_CK_corr, Title_TDC_CK_corr, 75, 200, 500);
	h_TDC_CK_BDC_Diff = new TH1D(Name_TDC_CK_BDC_Diff, Title_TDC_CK_BDC_Diff, 75, 1050, 1350);
	h_TDC_CK_BDC_Diff2 = new TH1D(Name_TDC_CK_BDC_Diff2, Title_TDC_CK_BDC_Diff2, 200, -50, 50);

	// TDC TARGET
	int index_max[2]={-1};
	int ADC_Low_max[2];	int ADC_Low_max_bar[2]={0};
	vector <int> vec_TDC_indexes;   vec_TDC_indexes.clear();
	vector <int> vec_tdc_LG_TARGET; vec_tdc_LG_TARGET.clear();
	vector <int> vec_tdc_LG_TARGET_sorted; vec_tdc_LG_TARGET_sorted.clear();
	vector <int> vec_index_max;  vec_index_max.clear();
    int adc_high_target[256];       int ADC_High_TARGET[256];    
	int tdc_le_target[256][16];     int TDC_LE_TARGET[256];     
	int adc_low_target[256];        int ADC_Low_TARGET[256];
	int LG_TARGET_ADC_Thr[256] = {0};
	int max_ADC_all[256];   int max_index_all[256];
	int max_index_flag;
  	int TDC_average = -1;
   	int TDC_LG_max1 = -1;    int TDC_LG_max2 = -1;
   	double TDC_diff = -1;
	fChain->SetBranchAddress("ADC_High_TARGET", adc_high_target);
	fChain->SetBranchAddress("TDC_LE_TARGET", tdc_le_target);
	fChain->SetBranchAddress("ADC_Low_TARGET" ,adc_low_target);
	TH1D *h_TDC_Average; char Title_TDC_Average[100]; char Name_TDC_Average[100];
	TH1D *h_TDC_Diff;    char Title_TDC_Diff[100];    char Name_TDC_Diff[100];
	sprintf(Title_TDC_Average,"TDC Average");
	sprintf(Name_TDC_Average,"h_TDC_Average (Run %d)", run_number);
	sprintf(Title_TDC_Diff,"TDC Diff");
	sprintf(Name_TDC_Diff,"h_TDC_Diff (Run %d)", run_number);
	h_TDC_Average = new TH1D(Name_TDC_Average, Title_TDC_Average, 300, 700,1000);
	h_TDC_Diff = new TH1D(Name_TDC_Diff, Title_TDC_Diff, 4000, -2000, 2000);

	// Time Walk Correction
	double par_in[256][3] = {0.};
	double par_err[256][3] = {0.};
	int ADCmax = 3450;
	double Yfi = 0.;
	double Ani = 0.;
	double Yni = 0.;
	double Tpi = 0.;
	bool NoFile = false;
	int T_limit = 3;

	int max1 = 0;
	int max2 = 0;
    int index_max1 = 0;
    int index_max2 = 0;

	char ParsTarg1[100];
	sprintf(ParsTarg1, "TimeWalk%d.dat", run_number);

	if(!ifstream(ParsTarg1)) NoFile = true;
	cout << endl; 
	cout << endl;
	cout << "***********************************************" << endl;
	cout << "File opened:  " << Name_finput << endl;
	if(!NoFile) cout << "Time Walk Parameter File:  " << ParsTarg1 << endl;
  	if(NoFile) cout << "Time Walk Parameter File:  " << "NOT FOUND !" << endl;

	// Get the total number of entries
	int nentries = (int)fChain->GetEntries();
	cout << "Total Number of Entries :  " << nentries << endl;
	cout << endl;

	// Read in Time Walk parameters
	ifstream parTARGdat(ParsTarg1,ios::in);
	int ij = 0; int ik = 0; 

	for(int ii=0; ii<256; ii++){
		parTARGdat >> ij >> ik >> par_in[ii][0] >> par_err[ii][0];
		parTARGdat >> ij >> ik >> par_in[ii][1] >> par_err[ii][1];
		parTARGdat >> ij >> ik >> par_in[ii][2] >> par_err[ii][2];
	}

	// Loop over events
	for(int ivt=0; ivt<nentries; ivt++){
	//for(int ivt=0; ivt<21; ivt++){
		vec_tdc_bdc.clear();
		tdc_bdc_selected = -1;
		int temp_max = -1;

		ADCmax = 3450;
		Yfi = 0.; Ani = 0.; Yni = 0.; Tpi = 0.; NoFile = false;
		ij = 0; ik = 0;

		TDC_LG_max1 = -1;
  		TDC_LG_max2 = -1;
  		index_max1=0;
  		index_max2=0;
  		TDC_average = -1;

		fChain->GetEntry(ivt);
		if(ivt%10000==1) cout << "**** " << ivt << " events done !" << endl;

		h_TDC_Trigger->Fill(tdc_trigger[0][0]);
		h_TDC_Trig_VT->Fill(tdc_vt48[0][0]);

		for(int i=0; i<256; i++){
			LG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Ped_LG[i]) + TARGET_ADC_Thr_LG_Offset;
  			TDC_LE_TARGET[i]=tdc_le_target[i][0];
  			ADC_Low_TARGET[i]=adc_low_target[i]-LG_TARGET_ADC_Thr[i];
  			vec_tdc_LG_TARGET.push_back(adc_low_target[i]-round(TARGET_ADC_Ped_LG[i]));
  			vec_tdc_LG_TARGET_sorted.push_back(adc_low_target[i]-round(TARGET_ADC_Ped_LG[i]));
  			
  			max_ADC_all[i] = -1000000;
  			max_index_all[i] = -1;

  			Yfi = par_in[i][0] - par_in[i][1]/sqrt(ADCmax - par_in[i][2]);
        	Ani = adc_high_target[i]-TARGET_ADC_Ped_HG[i]; //SEB

        	if((Ani >= TARGET_ADC_Thr_HG_Offset) && (Ani < ADCmax)){ // SEB
         		Yni = par_in[i][0] - par_in[i][1]/sqrt(Ani - par_in[i][2]);
          		Tpi = Yfi - Yni;
          		for(int jj=0; jj<16; jj++){
            		if(tdc_le_target[i][jj]>0) tdc_le_target[i][jj] = tdc_le_target[i][jj] + Tpi;
          		}
        	}

        	for(int j=0; j<16; j++){
        		if(tdc_le_target[i][j]>=TDC_Thr_min && tdc_le_target[i][j]<=TDC_Thr_max){
        			vec_TDC_indexes.push_back(i);
        		}
        	}

        	max1 = -9999;
        	index_max1 = 0;
        	max2 = -9999;
        	index_max2 = 0;
       	}
        	
        for(unsigned int k=0; k<vec_TDC_indexes.size(); k++){
        	if(ADC_Low_TARGET[vec_TDC_indexes[k]]>max1){
         		max1 = ADC_Low_TARGET[vec_TDC_indexes[k]];
        		index_max1 = vec_TDC_indexes[k]; 
        	}
        }
        
        for(unsigned int k=0; k<vec_TDC_indexes.size(); k++){
        	if(ADC_Low_TARGET[vec_TDC_indexes[k]]>max2 && vec_TDC_indexes[k]!=index_max1){
         		max2 = ADC_Low_TARGET[vec_TDC_indexes[k]];
        		index_max2 = vec_TDC_indexes[k]; 
        	}
        }

        for(int j=0; j<16; j++){
        	if(tdc_le_target[index_max1][j]>=TDC_Thr_min && tdc_le_target[index_max2][j]<=TDC_Thr_max){
        		TDC_LG_max1 = tdc_le_target[index_max1][j];
        	}
        	if(vec_TDC_indexes.size()>1 && tdc_le_target[index_max2][j]>=TDC_Thr_min && tdc_le_target[index_max2][j]<=TDC_Thr_max){
        		TDC_LG_max2 = tdc_le_target[index_max2][j];
        	}      	
        	if(vec_TDC_indexes.size()==0) TDC_LG_max1 = -1;
         	if(vec_TDC_indexes.size()<=0) TDC_LG_max2 = -1;
       }

    	if(vec_TDC_indexes.size()>0){
  			if(abs(TDC_LG_max1 - TDC_LG_max2) <= T_limit) TDC_average = (TDC_LG_max1 + TDC_LG_max2)/2;
        	else if(TDC_LG_max1 > TDC_LG_max2) TDC_average = TDC_LG_max1;
  			else if(TDC_LG_max1 < TDC_LG_max2) TDC_average = TDC_LG_max2;
  		}
  		else TDC_average = -1;

		sort(vec_tdc_LG_TARGET_sorted.begin(), vec_tdc_LG_TARGET_sorted.end());
		ADC_Low_max[0] = vec_tdc_LG_TARGET_sorted[255];
		ADC_Low_max[1] = vec_tdc_LG_TARGET_sorted[254];

		for(unsigned i=0; i<vec_tdc_LG_TARGET.size(); i++){
			for(int j=0; j<2; j++){
				if(vec_tdc_LG_TARGET[i] == ADC_Low_max[j]) ADC_Low_max_bar[j] = i;
			}
		}

		for(int i=0; i<14; i++){
			for(int j=0; j<16; j++){
				if(tdc_ck[i][j]>=TDC_Ck_min && tdc_ck[i][j]<=TDC_Ck_max){
					vec_tdc_ck_col.push_back(tdc_ck[i][j]);
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
  
  		if(vec_tdc_ck.size()>0) TDC_ck_sigma = sqrt(TDC_ck_sigma2/vec_tdc_ck.size());
  		else TDC_ck_sigma = -1;
 
  		//cout << "PAULINE : " << vec_tdc_ck.size() << endl;
  		for(unsigned int i=0; i<vec_tdc_ck.size(); i++){
    		if(abs(vec_tdc_ck[i]-TDC_ck_avg) <= 1.4*TDC_ck_sigma){
      			TDC_ck_sum2 += vec_tdc_ck[i];
      			TDC_ck_counter++;
    		}
  		}

  		if(TDC_ck_counter>0) TDC_ck_avg2 = double(TDC_ck_sum2)/double(TDC_ck_counter);
  		else TDC_ck_avg2 = -1;

		ck_mean = double(TDC_ck_sum) / double(vec_tdc_ck.size());
		//tdc_ck_corr = tdc_trigger[0][0] - ck_mean;
		tdc_ck_corr = tdc_trigger[0][0] - TDC_ck_avg2;

		TDC_diff = -1.;
		TDC_diff = 550 - (0.91*TDC_average - 0.625*tdc_ck_corr);
  		
		h_TDC_Average->Fill(TDC_average);
		h_TDC_Diff->Fill(TDC_diff);
		h_TDC_CK->Fill(ck_mean);
		h_TDC_CK_corr->Fill(tdc_ck_corr);

		
		for(int i=0; i<16; i++){
			if(tdc_vt48[4][i]>TDC_Ck_min && tdc_vt48[4][i]<TDC_Ck_max){
				vec_tdc_bdc.push_back(tdc_vt48[4][i]);
			}

			if(tdc_vt48[6][i]>=TDC_Ck_min && tdc_vt48[6][i]<=TDC_Ck_max){
				vec_tdc_b0_6.push_back(tdc_vt48[6][i]);
			}
			else vec_tdc_b0_6.push_back(-1);

            if(tdc_vt48[7][i]>=TDC_Ck_min && tdc_vt48[7][i]<=TDC_Ck_max){
				vec_tdc_b0_7.push_back(tdc_vt48[7][i]);
            }
            else vec_tdc_b0_7.push_back(-1);
		}

		sort(vec_tdc_b0_6.begin(), vec_tdc_b0_6.end());
		sort(vec_tdc_b0_7.begin(), vec_tdc_b0_7.end());

		tdc_b0_6 = vec_tdc_b0_6[vec_tdc_b0_6.size()-1];
		tdc_b0_7 = vec_tdc_b0_7[vec_tdc_b0_7.size()-1];

		if(tdc_b0_6>0 && tdc_b0_7>0) tdc_b0 = 0.5*(tdc_b0_6+tdc_b0_7);
		if(tdc_b0_6 == -1) tdc_b0 = tdc_b0_7;
		if(tdc_b0_7 == -1) tdc_b0 = tdc_b0_6;
		if(tdc_b0_6==-1 && tdc_b0_7==-1) tdc_b0 = -1;

		if(tdc_b0<0) h_TDC_B0_NULL->Fill(tdc_b0);
		else h_TDC_B0->Fill(tdc_b0);

		for(unsigned int v=0; v<vec_TDC_indexes.size(); v++){
 			vec_tdc_LG_TARGET.push_back(adc_low_target[vec_TDC_indexes[v]]-round(TARGET_ADC_Ped_LG[vec_TDC_indexes[v]]));
  			vec_tdc_LG_TARGET_sorted.push_back(adc_low_target[vec_TDC_indexes[v]]-round(TARGET_ADC_Ped_LG[vec_TDC_indexes[v]]));
		}

		sort(vec_tdc_LG_TARGET_sorted.begin(), vec_tdc_LG_TARGET_sorted.end(), greater<int>());

		for(unsigned int r=0; r<vec_tdc_LG_TARGET_sorted.size(); r++){
			if(vec_tdc_LG_TARGET[r] == vec_tdc_LG_TARGET_sorted[0]) index_max[0] = r;
			if(vec_tdc_LG_TARGET[r] == vec_tdc_LG_TARGET_sorted[1]) index_max[1] = r;
		}

		ADC_Low_max_bar[0] = vec_TDC_indexes[index_max[0]];
		ADC_Low_max_bar[1] = vec_TDC_indexes[index_max[1]];

		if(vec_tdc_bdc.size()>0){
			temp_max = vec_tdc_bdc[vec_tdc_bdc.size()-1];
			tdc_bdc_selected = tdc_trigger[0][0] - temp_max;
			h_TDC_BDC_Raw->Fill(temp_max);
			h_TDC_BDC_Diff->Fill(tdc_trigger[0][0] - temp_max);
			h_TDC_CK_BDC_Diff->Fill(ck_mean - tdc_bdc_selected);
			h_TDC_CK_BDC_Diff2->Fill(tdc_ck_corr - tdc_bdc_selected);
		}
		else h_TDC_BDC_NULL->Fill(temp_max);

		vec_tdc_LG_TARGET.clear();
		vec_tdc_LG_TARGET_sorted.clear();
		vec_TDC_indexes.clear();
		vec_tdc_b0_6.clear();
		vec_tdc_b0_7.clear();

		vec_tdc_ck.clear();
		TDC_ck_sum = 0;
		TDC_ck_sum2 = 0;
		TDC_ck_counter = 0; 
		TDC_ck_avg = 0;
		TDC_ck_avg2 = 0;
		TDC_ck_sigma = 0.;
		TDC_ck_sigma2 = 0.;
		ck_mean = 0.;

	} // EndLoop over events

	TCanvas *c1;
  	c1 = new TCanvas("Title","Title",0,200,700,800);
  	c1->Divide(3,5);
  	c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

  	c1->cd(1);
	h_TDC_Trigger->Draw();

	c1->cd(2);
	h_TDC_Trig_VT->Draw();

	c1->cd(4);
	h_TDC_BDC_Raw->Draw();

	c1->cd(5);
	h_TDC_BDC_NULL->Draw();

	c1->cd(6);
	h_TDC_BDC_Diff->Draw();

	c1->cd(7);
	h_TDC_CK->Draw();

	c1->cd(8);
	h_TDC_CK_corr->Draw();

	c1->cd(10);
	h_TDC_CK_BDC_Diff->Draw();

	c1->cd(11);
	h_TDC_CK_BDC_Diff2->Draw();

	c1->cd(13);
	h_TDC_Average->Draw();

	c1->cd(14);
	h_TDC_CK_corr->Draw();

	c1->cd(15);
	h_TDC_Diff->Draw();

	TCanvas *c2;
  	c2 = new TCanvas("Title2","Title2",0,200,700,800);
  	c2->Divide(2,2);
  	c2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

 	c2->cd(1);
	h_TDC_BDC_Raw->Draw();

	c2->cd(2);
	h_TDC_BDC_NULL->Draw();

	c2->cd(3);
	h_TDC_B0->Draw();

	c2->cd(4);
	h_TDC_B0_NULL->Draw();
 	

	cout << endl;

	return;

}