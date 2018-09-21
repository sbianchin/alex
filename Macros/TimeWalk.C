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
#include "ANAPATH.h"
#include "Pedestals.h"
#include "Batch_Variables.h"
#endif

using namespace std;


void TimeWalk(Int_t Run_Number=5, Int_t ievt_min=0, Int_t ievt_max=10){

    Int_t nBars = 256;
    double par_in[256][3] = {0.};
    double par_err[356][3] = {0.};
    Int_t ADCmax = 3450;
    Int_t HGpedOffset = 10;
    double Yfi = 0.;
    double Ani = 0.;
    double Yni = 0.;
    double Tpi = 0.;
    Int_t TDCnew[256] = {-1};

    bool NoFile = false;



    // Get Target TDC vs ADC Fit Parameters
    char ParsTarg1[100];
    //sprintf(ParsTarg1,"%sTimeWalk%d.dat",anal_save_path,RunNum);
    sprintf(ParsTarg1,"/media/bianchin/hdd1/trek/E36/TEST/offline_new/sandbox/TimeWalk%d.dat",Run_Number);
    sprintf(path_input,"%s",path_merged);          
    sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);

    cout << endl;
    cout << "Opened File:  " << endl;
    cout << Name_finput << endl;
    cout << "" << endl;
    cout << "TimeWalk fit parameters file:   " << endl;

    TH2F *h_TDC_vs_ADC_beforeTW[256];    char Title_TDC_vs_ADC_beforeTW[256][100];   char Name_TDC_vs_ADC_beforeTW[256][100];
    TH2F *h_TDC_vs_ADC_afterTW[256];     char Title_TDC_vs_ADC_afterTW[256][100];    char Name_TDC_vs_ADC_afterTW[256][100];

    for(int i=0; i<256; i++){
        sprintf(Title_TDC_vs_ADC_beforeTW[i],"TDC vs ADC (Ch. %d)  --  BEFORE Time Walk Correction  |  Run %d",i, Run_Number); 
        sprintf(Name_TDC_vs_ADC_beforeTW[i]," TDC vs ADC(Ch. %d) - BEFORE Time Walk Correction",i);
        h_TDC_vs_ADC_beforeTW[i] = new TH2F(Name_TDC_vs_ADC_beforeTW[i],Title_TDC_vs_ADC_beforeTW[i],512,0,3000,512,0,3000);

        sprintf(Title_TDC_vs_ADC_afterTW[i],"TDC vs ADC (Ch. %d)  --  AFTER Time Walk Correction  |  Run %d",i, Run_Number); 
        sprintf(Name_TDC_vs_ADC_afterTW[i]," TDC vs ADC(Ch. %d) - AFTER Time Walk Correction",i);
        h_TDC_vs_ADC_afterTW[i] = new TH2F(Name_TDC_vs_ADC_afterTW[i],Title_TDC_vs_ADC_afterTW[i],512,0,3000,512,0,3000);
    }
 
    TChain *fChain= new TChain("Tree");   
    //fChain->Add(Name_finput);   
    fChain->Add("/media/bianchin/hdd1/trek/E36/Data/Oct_2015/root/Merged/Run3994MS.root");   
    fChain->SetMakeClass(1);              
  
    fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);
    fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);

    if(std::ifstream(ParsTarg1)){
        cout << ParsTarg1 << "   FOUND! " << endl;
        cout << endl;
    }
    else{
        cout << ParsTarg1 << "   NOT FOUND! " << std::endl;
        cout << endl;
        NoFile = true;
        //break;
    }
    std::ifstream parTARGdat(ParsTarg1,std::ios::in);
    Int_t ii = 0;
    Int_t ij = 0;
    Int_t ik = 0;
    
    // Read in parameters and their errors. (errors not used)
    for(Int_t ii = 0; ii < nBars; ii++){
        parTARGdat >> ij >> ik >> par_in[ii][0] >> par_err[ii][0];
        //cout << ij << " " << ik << " " << par_in[ii][0] <<endl;
        parTARGdat >> ij >> ik >> par_in[ii][1] >> par_err[ii][1];
        //cout << ij << " " << ik << " " << par_in[ii][1] <<endl;
        parTARGdat >> ij >> ik >> par_in[ii][2] >> par_err[ii][2];
        //cout << ij << " " << ik << " " << par_in[ii][2] <<endl;
    }


    for(int ivt=ievt_min; ivt<ievt_max+1; ivt++){
        fChain->GetEntry(ivt);  
        if(ivt%10000==0 && ivt !=0) cout << "####  " << ivt << " events done!" << endl;

        for(Int_t ii = 0; ii<nBars; ii++){
            Yfi = par_in[ii][0] - par_in[ii][1]/sqrt(ADCmax - par_in[ii][2]);
            Ani = adc_high_target[ii]-Target_HG_ped[ii];

            //cout << "DEBUG:  " << ii << "  " << Yfi << "   " << Ani << "   " << adc_high_target[ii] << "   " << Target_HG_ped[ii] << endl;

            if((Ani >= Target_HG_ped[ii]) && (Ani < ADCmax)){
                Yni = par_in[ii][0] - par_in[ii][1]/sqrt(Ani - par_in[ii][2]);
                Tpi = Yfi - Yni;
                TDCnew[ii] = tdc_le_target[ii][0] + Tpi;
            }
            if((tdc_le_target[ii][0] > 700) && (tdc_le_target[ii][0] < 1200)) {
                //cout << "OK!  " << ivt << "  " << ii << "  " << endl;
                if((adc_high_target[ii] > (Target_HG_ped[ii] + HGpedOffset)) && (adc_high_target[ii] < ADCmax)){
                    //cout << "OK2!  " << ivt << "  " << ii << "  " << endl;
                h_TDC_vs_ADC_beforeTW[ii]->Fill(Ani,tdc_le_target[ii][0]);
                h_TDC_vs_ADC_afterTW[ii]->Fill(Ani,TDCnew[ii]);
                }
            }

        }
    }

    TFile *foutput;
    char Name_foutput[200];
    sprintf(Name_foutput,"TEST_Hist_Run%dM.root",Run_Number);
    foutput = new TFile(Name_foutput,"RECREATE");

    TDirectory *dir_before_TW = foutput->mkdir("BEFORE_TIME_WALK");
    TDirectory *dir_after_TW = foutput->mkdir("AFTER_TIME_WALK");

    foutput->cd();

    for (Int_t ii=0; ii<256; ii++){
        dir_before_TW->cd();
        h_TDC_vs_ADC_beforeTW[ii]->SetDirectory(dir_before_TW);
        h_TDC_vs_ADC_beforeTW[ii]->Write();

        dir_after_TW->cd();
        h_TDC_vs_ADC_afterTW[ii]->SetDirectory(dir_after_TW);
        h_TDC_vs_ADC_afterTW[ii]->Write();
     }

    cout << endl;

    return;
}