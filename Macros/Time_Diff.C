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
#include "TEllipse.h"
#include "TMarker.h"
#include "ANAPATH.h"
#include "Thresholds.h"
#include "CommonParameters.h"
#endif
void Time_Diff(Int_t Run_Number=5)
{ 

	gStyle->Clear();
	TH1::AddDirectory(kFALSE);
	gStyle->SetOptStat(1111);


 	Int_t adc_high_target[256];       Int_t ADC_High_TARGET[256];    
	Int_t adc_low_target[256];        Int_t ADC_Low_TARGET[256];  
	Int_t tdc_le_target[256][16];     Int_t TDC_LE_TARGET[256];     
	Int_t tdc_te_target[256][16];     //Int_t TDC_TE_TARGET[256];  

	Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];        Double_t ADC_High_SFT_corr[128];    
	Int_t adc_low_sft[128];           //Int_t ADC_Low_SFT[128];   
	Int_t tdc_le_sft[128][16];        Int_t TDC_LE_SFT[128];         
	Int_t tdc_te_sft[128][16];        //Int_t TDC_TE_SFT[128];  

	Int_t ADC_tof1[24];               Int_t ADC_TOF1[24]; 
	Int_t ADC_tof2[56];               Int_t ADC_TOF2[56];

	Int_t TDC_tof1U[12];              Int_t TDC_TOF1U[12];
	Int_t TDC_tof1D[12];              Int_t TDC_TOF1D[12];

	Float_t TDC_TOF1_Mean[12]={0.};

	Int_t TDC_tof2AO[12];             Int_t TDC_TOF2AO[12];
	Int_t TDC_tof2BO[12];             Int_t TDC_TOF2BO[12];
	Int_t TDC_tof2AI[12];             Int_t TDC_TOF2AI[12];
	Int_t TDC_tof2BI[12];             Int_t TDC_TOF2BI[12];

	Int_t MwpcADC[512];               Int_t MWPCADC[512];

	char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!

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


  	cout << "   " << endl;
  	cout << "File opened:  " << Name_finput << endl;
  	cout << "SFT Mapping File:  " << par_finput << endl;
  	cout << "MWPC Mapping File:  " << par_finput2 << endl;

  	//Change variables by applying Mapping
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
  
	// Histogram Declarations
	TH1F *h_TDC_SFT[128];			char h_TDC_SFT_Name[128][200];			char h_TDC_SFT_Title[128][200];
	TH1F *h_TDC_TOF1_UP[12];		char h_TDC_TOF1_UP_Name[12][200];		char h_TDC_TOF1_UP_Title[12][200];
	TH1F *h_TDC_TOF1_DOWN[12];		char h_TDC_TOF1_DOWN_Name[12][200];		char h_TDC_TOF1_DOWN_Title[12][200];
	TH1F *h_TDC_TOF1_MEAN[12];		char h_TDC_TOF1_MEAN_Name[12][200];		char h_TDC_TOF1_MEAN_Title[12][200];
	TH1F *h_Diff[128];				char h_Diff_Name[128][200];				char h_Diff_Title[128][200];

	for(int i=0; i<128; i++)
	{
		sprintf(h_TDC_SFT_Title[i], "SFT TDC Channel %d  |  Run %d", i, Run_Number);
		sprintf(h_TDC_SFT_Name[i], "Channel %d", i);

		sprintf(h_Diff_Title[i], "TDC(TOF1_Mean) - TDC(SFT)  |  Channel %d  --  Run %d", i, Run_Number);
		sprintf(h_Diff_Name[i], "Channel %d", i);

		h_TDC_SFT[i] = new TH1F(h_TDC_SFT_Name[i], h_TDC_SFT_Title[i], 1200, 0, 1200);
		h_Diff[i] = new TH1F(h_Diff_Name[i], h_Diff_Title[i], 500, 500, 1000);
	}

	for(int j=0; j<12; j++)
	{
		sprintf(h_TDC_TOF1_UP_Name[j], "Channel %d", j);
		sprintf(h_TDC_TOF1_UP_Title[j], "TOF1-UP TDC Channel %d  |  Run %d", j, Run_Number);

		sprintf(h_TDC_TOF1_DOWN_Name[j], "Channel %d", j);
		sprintf(h_TDC_TOF1_DOWN_Title[j], "TOF1-DOWN TDC Channel %d  |  Run %d", j, Run_Number);

		sprintf(h_TDC_TOF1_MEAN_Name[j], "Channel %d", j);
		sprintf(h_TDC_TOF1_MEAN_Title[j], "TOF1-MEAN TDC Channel %d  |  Run %d", j, Run_Number);

		h_TDC_TOF1_UP[j] = new TH1F(h_TDC_TOF1_UP_Name[j], h_TDC_TOF1_UP_Title[j], 3000, 0, 3000);
		h_TDC_TOF1_DOWN[j] = new TH1F(h_TDC_TOF1_UP_Name[j], h_TDC_TOF1_DOWN_Title[j], 3000, 0, 3000);
		h_TDC_TOF1_MEAN[j] = new TH1F(h_TDC_TOF1_MEAN_Name[j], h_TDC_TOF1_MEAN_Title[j], 3000, 0, 3000);
	}


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
	cout << "  " << endl;
	cout << "*****************************************************" << endl;
	cout << "  " << endl;

	//for(Int_t ievent=0; ievent<1; ievent++)		// Loop over Events
	for(Int_t ievent=0; ievent<nentries; ievent++)
	{
    	fChain->GetEntry(ievent);  

    	if(ievent%10000==1) cout << "****  " << ievent << "  events done" << endl;

    	for(int init=0; init<12; init++)	TDC_TOF1_Mean[init]=0;

    	for(int j_TARGET=0; j_TARGET<256; j_TARGET++){
      		ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-ADC_cut_TARGET;
      		ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-ADC_cut_TARGET;
      		TDC_LE_TARGET[j_TARGET]=tdc_le_target[j_TARGET][0];
      		//TDC_TE_TARGET[j_TARGET]=tdc_te_target[j_TARGET][0];  
    	}

    	for(Int_t j_SFT=0; j_SFT<128; j_SFT++){
      		ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-par_temp_SFT[1][j_SFT];
      		//ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-par_temp_SFT[1][j_SFT];
      		TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
      		//TDC_TE_SFT[j_SFT]=tdc_te_sft[j_SFT][0]; 
    	}

    	for(Int_t j_TOF1=0; j_TOF1<24; j_TOF1++) {
      		ADC_TOF1[j_TOF1] = ADC_tof1[j_TOF1]-par_temp_TOF1[1][j_TOF1];
    	}

    	for(Int_t j_TOF2=0; j_TOF2<56; j_TOF2++) {
      		ADC_TOF2[j_TOF2] = ADC_tof2[j_TOF2]-par_temp_TOF2[1][j_TOF2];
    	}

    	for(Int_t j_MWPC=0; j_MWPC<512; j_MWPC++) {
      		MWPCADC[j_MWPC] = MwpcADC[j_MWPC]-MWPC_thr;
    	}

    	for(Int_t j_TDCTOF=0; j_TDCTOF<12; j_TDCTOF++) {
      		TDC_TOF1U[j_TDCTOF] = TDC_tof1U[j_TDCTOF];
      		TDC_TOF1D[j_TDCTOF] = TDC_tof1D[j_TDCTOF];
      		TDC_TOF2AO[j_TDCTOF] = TDC_tof2AO[j_TDCTOF];
      		TDC_TOF2BO[j_TDCTOF] = TDC_tof2BO[j_TDCTOF];
      		TDC_TOF2AI[j_TDCTOF] = TDC_tof2AI[j_TDCTOF];
      		TDC_TOF2BI[j_TDCTOF] = TDC_tof2BI[j_TDCTOF];

      		if((TDC_TOF1U[j_TDCTOF]>=TDC_TOF1_min && TDC_TOF1U[j_TDCTOF]<=TDC_TOF1_max) && 
      			(TDC_TOF1D[j_TDCTOF]>=TDC_TOF1_min && TDC_TOF1D[j_TDCTOF]<=TDC_TOF1_max))
      		{
      			TDC_TOF1_Mean[j_TDCTOF] = 0.5 * (TDC_TOF1U[j_TDCTOF] + TDC_TOF1D[j_TDCTOF]);
      		}

      		if((TDC_TOF1U[j_TDCTOF]>=TDC_TOF1_min && TDC_TOF1U[j_TDCTOF]<=TDC_TOF1_max) && 
      			(TDC_TOF1D[j_TDCTOF]<TDC_TOF1_min || TDC_TOF1D[j_TDCTOF]>TDC_TOF1_max))
      		{
      			TDC_TOF1_Mean[j_TDCTOF] = TDC_TOF1U[j_TDCTOF];
      		}

      		if((TDC_TOF1D[j_TDCTOF]>=TDC_TOF1_min && TDC_TOF1D[j_TDCTOF]<=TDC_TOF1_max) && 
      			(TDC_TOF1U[j_TDCTOF]<TDC_TOF1_min || TDC_TOF1U[j_TDCTOF]>TDC_TOF1_max))
      		{
      			TDC_TOF1_Mean[j_TDCTOF] = TDC_TOF1D[j_TDCTOF];
      		}   		

          	//	cout << "Event  " << ievent << ":  " << "TDC TOF1 Mean Channel " << j_TDCTOF << " :  " << TDC_TOF1_Mean[j_TDCTOF] << endl;
    	}
  	
    	for(int i=0; i<128; i++)
    	{
    		for(int j=0; j<12; j++)
    		{
    			if(TDC_LE_SFT[i]>=TDC_min_SFT && TDC_LE_SFT[i]<=TDC_max_SFT && TDC_TOF1_Mean[j]>0)
    			{
    		//		cout << ievent << "  " << i << "  " << j << "  " << TDC_LE_SFT[i] << "  " << TDC_TOF1_Mean[j] << "   " << float(TDC_LE_SFT[i]) - 0.025*float(TDC_TOF1_Mean[j])<< endl;
    				h_Diff[i]->Fill(float(TDC_LE_SFT[i]) - 0.025*float(TDC_TOF1_Mean[j]));
   			}
    		}
    	}


    	//cout << "  " << endl;

		for(int j=0; j<128; j++)
		{ 
			if(ADC_High_SFT[j]<0)      ADC_High_SFT_corr[j]=0; 
			if(ADC_High_SFT[j]>=0)     ADC_High_SFT_corr[j]=ADC_High_SFT[j]; 
		}

		//for(int k=0; k<128; k++)
		//{
		//	cout << "SFT Channel " << k << " (TDC):  " << TDC_LE_SFT[k] << endl;;
		//}

		//cout << "  " << endl;
		//cout << "///////////////////////////////////////////////" << endl;
		//cout << "  " << endl;

		//for(int k=0; k<12; k++)
		//{
		//	cout << "TOF1-Up " << k << " (TDC):  " << TDC_TOF1U[k] << endl;
		//}

		//cout << "  " << endl;
		//cout << "///////////////////////////////////////////////" << endl;
		//cout << "  " << endl;

		//for(int k=0; k<12; k++)
		//{
		//	cout << "TOF1-Down " << k << " (TDC):  " << TDC_TOF1D[k] << endl;
		//}


		for(int j=0; j<128; j++)
		{
			if(TDC_LE_SFT[j]>=0)
			{
				h_TDC_SFT[j]->Fill(TDC_LE_SFT[j]);
			}
		}

		for(int k=0; k<12; k++)
		{
			if(TDC_TOF1U[k]>=0)
			{
				h_TDC_TOF1_UP[k]->Fill(TDC_TOF1U[k]);
			}

			if(TDC_TOF1D[k]>=0)
			{
				h_TDC_TOF1_DOWN[k]->Fill(TDC_TOF1D[k]);
			}

			if(TDC_TOF1_Mean[k]>0)
			{
				h_TDC_TOF1_MEAN[k]->Fill(TDC_TOF1_Mean[k]);
			}
		}


	} // EndLoop over Events

	// Canvas Declaration

	char c_Diff1_Title[100];
	char c_Diff2_Title[100];

	sprintf(c_Diff1_Title,"TDC(TOF1_Mean) - TDC(SFT)  |  Run %d  --  Channels 0 - 63", Run_Number);
	sprintf(c_Diff2_Title,"TDC(TOF1_Mean) - TDC(SFT)  |  Run %d  --  Channels 63 - 127", Run_Number);
	
	//TCanvas *c_SFT;
	//c_SFT = new TCanvas("SFT","SFT",0,200,1050,700);
	//c_SFT->Divide(16,8);
	//c_SFT->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	//TCanvas *c_TOF1;
	//c_TOF1 = new TCanvas("TOF1","TOF1",0,200,1050,700);
	//c_TOF1->Divide(6,4);
	//c_TOF1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	//TCanvas *c_TOF1_MEAN;
	//c_TOF1_MEAN = new TCanvas("TOF1 MEAN","TOF1 MEAN",0,200,1050,700);
	//c_TOF1_MEAN->Divide(6,4);
	//c_TOF1_MEAN->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	//TCanvas *c_Diff;
	//c_Diff = new TCanvas("Diff","Diff",0,200,1050,700);
	//c_Diff->Divide(16,8);
	//c_Diff->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	TCanvas *c_Diff2;
	c_Diff2 = new TCanvas("Diff2",c_Diff2_Title,0,200,1050,700);
	c_Diff2->Divide(8,8);
	c_Diff2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	TCanvas *c_Diff1;
	c_Diff1 = new TCanvas("Diff1",c_Diff1_Title,0,200,1050,700);
	c_Diff1->Divide(8,8);
	c_Diff1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	TLatex *tex_condition1[128];
	TLatex *tex_condition2[128];

	char cond1[100];	sprintf(cond1,"%d #leq TDC(TOF1) #leq %d", TDC_TOF1_min, TDC_TOF1_max);
	char cond2[100];	sprintf(cond2,"%d #leq TDC(SFT) #leq %d", TDC_min_SFT, TDC_max_SFT);

	int max=0;

	//tex_condition1 = new TLatex(-744.,0.99, cond1);
  	//tex_condition1->SetTextSize(0.03);
  	//tex_condition1->SetLineWidth(2);

  	//tex_condition2 = new TLatex(-744.,0.94, cond2);
  	//tex_condition2->SetTextSize(0.03);
  	//tex_condition2->SetLineWidth(2);
	
	for(int ii=0; ii<128; ii++)
	{
		max = h_Diff[ii]->GetMaximum();
		tex_condition1[ii] = new TLatex(520., 0.99*max, cond1);
	  	tex_condition1[ii]->SetTextSize(0.03);
 		tex_condition1[ii]->SetLineWidth(2);
	  	
	  	tex_condition2[ii] = new TLatex(520.,0.94*max, cond2);
  		tex_condition2[ii]->SetTextSize(0.03);
  		tex_condition2[ii]->SetLineWidth(2);


	}

	/*
	for(int ican1=0; ican1<128; ican1++)
	{
		c_SFT->cd(ican1+1);
		h_TDC_SFT[ican1]->Draw();
	}

	for(int ican2=0; ican2<12; ican2++)
	{
		c_TOF1->cd(ican2+1);
		h_TDC_TOF1_UP[ican2]->Draw();

		c_TOF1->cd(ican2+13);
		h_TDC_TOF1_DOWN[ican2]->Draw();
	}

	for(int ican3=0; ican3<12; ican3++)
	{
		c_TOF1_MEAN->cd(ican3+1);
		h_TDC_TOF1_MEAN[ican3]->Draw();
	}

	for(int ican4=0; ican4<128; ican4++)
	{
		c_Diff->cd(ican4+1);
		h_Diff[ican4]->Draw();
	}
	*/
	
	for(int ican5=0; ican5<64; ican5++)
	{
		c_Diff1->cd(ican5+1);
		h_Diff[ican5]->Draw();
		tex_condition1[ican5]->Draw();
		tex_condition2[ican5]->Draw();
	}
	
	for(int ican6=0; ican6<64; ican6++)
	{
		c_Diff2->cd(ican6+1);
		h_Diff[ican6+64]->Draw();
		tex_condition1[ican6+64]->Draw();
		tex_condition2[ican6+64]->Draw();
	}


} //EndVoid