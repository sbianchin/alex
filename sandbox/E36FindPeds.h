//////////////////////////////////////////////////////////
// E36FindPeds
// Finds pedestals and calculates pedestal removal cut
// Finds TDC peak and calculates time cuts (upper and lower)
// Results output to Target and SFT pedestal files.
//////////////////////////////////////////////////////////

#ifndef E36FindPeds_h
#define E36FindPeds_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <stdio.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include <TObjArray.h>
#include <TVector3.h>
#include <TVector.h>
#include <TStyle.h>
#include "TF1.h"
#include "TH2.h"
#include "TH1.h"
#include "TStyle.h"
#include "TMath.h"
#include "CommonParameters.h"

//Some Pedestal Finder Global Variables
   TFile *hfout;
//Pedestal results Target
Double_t adcPeakT[nBars], adcRMST[nBars];
Double_t ladcPeakT[nBars], ladcRMST[nBars];
Double_t tdcPeakT[nBars], tdcRMST[nBars];
Int_t BarHitsT[nBars], BarProblemT[nBars], DeadBarIfLess = 5, MinHits = 5;
//Pedestal results SFT
Double_t adcPeakF[nFibrs], adcRMSF[nFibrs];
Double_t ladcPeakF[nFibrs], ladcRMSF[nFibrs];
Double_t tdcPeakF[nFibrs], tdcRMSF[nFibrs];
Int_t BarHitsF[nFibrs], BarProblemF[nFibrs];

Int_t nChi = 10;
Int_t nEvG=0; // number of events to go through

// Pedestal output files
char Pedft[100];
char Pedff[100];
// filename is (path)file containing input trees.
char filename[100];
Int_t NoFile = 1;  // Didn't find the file if < 0
//histos output root file
  char histo_savefile[100];

// Upper and lower level cuts to concentrate on pedestals. Input via Cuts().
Int_t HGcutT = 1000, LGcutT = 1000, MinCutT = 100;
Int_t HGcutF = 1000, LGcutF = 1000, MinCutF = 100;
// Cuts for finding means before doing fits
Int_t T_LowMean = 500, T_HighMean = 1000;
Int_t F_LowMean = 500, F_HighMean = 1000;
// Minimum number of counts in tdc
Int_t T_mintdc = 50, F_mintdc = 50;
// Minimum number of counts in adc
Int_t T_minhighadc = 1000, F_minhighadc =1000;
Int_t T_minlowadc = 1000, F_minlowadc =1000;
//Number of RMS to add to peak mean to calculate cuts
Double_t NPsig = 5, NTsig = 2; //.

// Define Histograms
TH1F *htdcT[nBars];
TH1F *hadcT[nBars];
TH1F *hadcLT[nBars];
TH1F *htdcF[nFibrs];
TH1F *hadcF[nFibrs];
TH1F *hadcLF[nFibrs];
TH1F *hTARGvsSFT;
TH1F *hTARG_HITS;
TH1F *hSFT_HITS;
TH1F *hMissBarT;
TH1F *hMissBarF;
TH2F *hMeanBarT; // high gain Monitor Target vs bar number
TH2F *hMeanBarF; // high gain Monitor SFT vs fibre number
TH2F *hRMSbarT; // hign gain RMS vs bar number
TH2F *hRMSbarF; // hign gain RMS vs fibre number
TH2F *hMeanBarTL; // low gain MonitorTarget vs bar number
TH2F *hMeanBarFL; // low gain Monitor SFT vs fibre number
TH2F *hRMSbarTL; // low gain RMS vs bar number
TH2F *hRMSbarFL; // low gain RMS vs fibre number
TH2F *hMeanTDCT; // tdc MonitorTarget vs bar number
TH2F *hMeanTDCF; // tdc Monitor SFT vs fibre number
TH2F *hRMSTDCT; // tdc RMS vs bar number
TH2F *hRMSTDCF; // tdc RMS vs fibre number

class E36FindPeds {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
    Int_t           ntrack_TARGET;
    Int_t           ADC_High_TARGET[256];
    Int_t           ADC_Low_TARGET[256];
    Int_t           TDC_LE_TARGET[256][16];
    Int_t           TDC_TE_TARGET[256][16];
    //SFT
    Int_t           ntrack_SFT;
    Int_t           ADC_High_SFT[128];
    Int_t           ADC_Low_SFT[128];
    Int_t           TDC_LE_SFT[128][16];
    Int_t           TDC_TE_SFT[128][16];
    //TOF
    Int_t           ntrack_TOF1;
    Int_t           ADC_TOF1[24];   //[ntrack_TOF1]
    Int_t           ntrack_TOF2;
    Int_t           ADC_TOF2[48];   //[ntrack_TOF2]

    // List of branches
    TBranch        *b_ntrack_TARGET;   //!
    TBranch        *b_ADC_High_TARGET;   //!
    TBranch        *b_ADC_Low_TARGET;   //!
    TBranch        *b_TDC_LE_TARGET;   //!
    TBranch        *b_TDC_TE_TARGET;   //!
    //SFT
    TBranch        *b_ntrack_SFT;   //!
    TBranch        *b_ADC_High_SFT;   //!
    TBranch        *b_ADC_Low_SFT;   //!
    TBranch        *b_TDC_LE_SFT;   //!
    TBranch        *b_TDC_TE_SFT;   //!
    //TOF
    TBranch        *b_ntrack_TOF1;   //!
    TBranch        *b_ADC_TOF1;   //!
    TBranch        *b_ntrack_TOF2;   //!
    TBranch        *b_ADC_TOF2;   //!
    
   E36FindPeds(TTree *tree=0);
   virtual ~E36FindPeds();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Int_t FrstEv);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
    virtual void     TMakeHistos();
    virtual void     TFindHits();
    virtual void     TFindPeds();
};

#endif

#ifdef E36FindPeds_cxx
E36FindPeds::E36FindPeds(TTree *tree) : fChain(0)
{
// tree file looked for here

       printf("Run %s \n", filename);
   if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
       if (!f || !f->IsOpen()) {
           f = new TFile(filename);
       }
       f->GetObject("Tree",tree);
   }
   Init(tree);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
E36FindPeds::~E36FindPeds()
{
  if (!fChain) return;
   delete fChain->GetCurrentFile();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t E36FindPeds::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Long64_t E36FindPeds::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void E36FindPeds::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

//Target
    fChain->SetBranchAddress("ntrack_TARGET", &ntrack_TARGET, &b_ntrack_TARGET);
    fChain->SetBranchAddress("ADC_High_TARGET", ADC_High_TARGET, &b_ADC_High_TARGET);
    fChain->SetBranchAddress("ADC_Low_TARGET", ADC_Low_TARGET, &b_ADC_Low_TARGET);
    fChain->SetBranchAddress("TDC_LE_TARGET", TDC_LE_TARGET, &b_TDC_LE_TARGET);
    fChain->SetBranchAddress("TDC_TE_TARGET", TDC_TE_TARGET, &b_TDC_TE_TARGET);
//SFT
    fChain->SetBranchAddress("ntrack_SFT", &ntrack_SFT, &b_ntrack_SFT);
    fChain->SetBranchAddress("ADC_High_SFT", ADC_High_SFT, &b_ADC_High_SFT);
    fChain->SetBranchAddress("ADC_Low_SFT", ADC_Low_SFT, &b_ADC_Low_SFT);
    fChain->SetBranchAddress("TDC_LE_SFT", TDC_LE_SFT, &b_TDC_LE_SFT);
    fChain->SetBranchAddress("TDC_TE_SFT", TDC_TE_SFT, &b_TDC_TE_SFT);
//TOF
    fChain->SetBranchAddress("ntrack_TOF1", &ntrack_TOF1, &b_ntrack_TOF1);
    fChain->SetBranchAddress("ADC_TOF1", ADC_TOF1, &b_ADC_TOF1);
    fChain->SetBranchAddress("ntrack_TOF2", &ntrack_TOF2, &b_ntrack_TOF2);
    fChain->SetBranchAddress("ADC_TOF2", ADC_TOF2, &b_ADC_TOF2);
   Notify();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t E36FindPeds::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void E36FindPeds::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t E36FindPeds::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Make Histograms
void E36FindPeds::TMakeHistos() {
    //set up histos Target
    hTARGvsSFT = new TH1F("EvRatio","Event Ratio",200,-100,100);
    hTARG_HITS = new TH1F("Bars Hit","Bars Hit",260,0,260);
    hMissBarT = new TH1F("MissingBars","Missing Bars",260,0,260);
    hMeanBarT = new TH2F("HighMeansT","High ADC MeansT",4000,0,4000,256,0,256); //monitor gains
    hMeanBarT->SetMarkerStyle(29);
    hMeanBarTL = new TH2F("LowMeansT","Low ADC MeansT",4000,0,4000,256,0,256); //monitor gains
    hMeanBarTL->SetMarkerStyle(29);
    hRMSbarT = new TH2F("HighRMST","High ADC RMST",1000,0,1000,256,0,256); //monitor resolution
    hRMSbarT->SetMarkerStyle(29);
    hRMSbarTL = new TH2F("LowRMST","Low ADC RMST",1000,0,1000,256,0,256); //monitor resolution
    hRMSbarTL->SetMarkerStyle(29);
    hMeanTDCT = new TH2F("TDCMeansT","TDC MeansT",1024,0,1024,256,0,256); //monitor gains
    hMeanTDCT->SetMarkerStyle(29);
    hRMSTDCT = new TH2F("TDCRMST","TDC RMST",1000,0,1000,256,0,256); //monitor resolution
    hRMSTDCT->SetMarkerStyle(29);
    //set up histos SFT
    hSFT_HITS = new TH1F("Fibres Hit","Fifres Hit",130,0,130);
    hMissBarF = new TH1F("MissingFibres","Missing Fibres",130,0,130);
    hMeanBarF = new TH2F("HighMeansF","High ADC MeansF",4000,0,4000,130,0,130); //monitor gains
    hMeanBarF->SetMarkerStyle(29);
    hMeanBarFL = new TH2F("LowMeansF","Low ADC MeansF",4000,0,4000,130,0,130); //monitor gains
    hMeanBarFL->SetMarkerStyle(29);
    hRMSbarF = new TH2F("HighRMSF","High ADC RMSF",1000,0,1000,130,0,130); //monitor resolution
    hRMSbarF->SetMarkerStyle(29);
    hRMSbarFL = new TH2F("LowRMSF","Low ADC RMSF",1000,0,1000,130,0,130); //monitor resolution
    hRMSbarFL->SetMarkerStyle(29);
    hMeanTDCF = new TH2F("TDCMeansF","TDC MeansF",1024,0,1024,130,0,130); //monitor gains
    hMeanTDCF->SetMarkerStyle(29);
    hRMSTDCF = new TH2F("TDCRMSF","TDC RMSF",1000,0,1000,130,0,130); //monitor resolution
    hRMSTDCF->SetMarkerStyle(29);
   for(int ch = 0; ch<nBars; ++ch){
        std::ostringstream _name;
        _name << "adcT" << ch;
        std::string ahist_name = _name.str();
        hadcT[ch] = new TH1F(ahist_name.c_str(), ahist_name.c_str(),4096,0,4096);
        //-------------------------
        std::ostringstream _nameLT;
        _nameLT << "adcLT" << ch;
        std::string ahist_nameLT = _nameLT.str();
        hadcLT[ch] = new TH1F(ahist_nameLT.c_str(), ahist_nameLT.c_str(),4096,0,4096);
        //-------------------------
       std::ostringstream _namet;
       _namet << "tdcT" << ch;
       std::string ahist_namet = _namet.str();
       htdcT[ch] = new TH1F(ahist_namet.c_str(), ahist_namet.c_str(),1050,-10,1040);
       //-------------------------
   }
   for(int ch = 0; ch<nFibrs; ++ch){
        std::ostringstream _nameF;
        _nameF << "adcF" << ch;
        std::string ahist_nameF = _nameF.str();
        hadcF[ch] = new TH1F(ahist_nameF.c_str(), ahist_nameF.c_str(),4096,0,4096);
        //-------------------------
        std::ostringstream _nameLF;
        _nameLF << "adcLF" << ch;
        std::string ahist_nameLF = _nameLF.str();
        hadcLF[ch] = new TH1F(ahist_nameLF.c_str(), ahist_nameLF.c_str(),4096,0,4096);
        //-------------------------
       std::ostringstream _nameft;
       _nameft << "tdcF" << ch;
       std::string ahist_nameft = _nameft.str();
       htdcF[ch] = new TH1F(ahist_nameft.c_str(), ahist_nameft.c_str(),1050,-10,1040);
       //-------------------------
   }
   return;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Find Hits routine
void E36FindPeds::TFindHits() {
	//printf("Event %d \n", nEvent);
    // 1st fill hits histo.
    // Set counters etc to 0
    nhits =0;
    nhits_t =0;
    nhitsL = 0;
    nhitsF =0;
    nhitsLF = 0;
    nhits_Ft =0;
// Loop over input channels and look for hits.
    // #### For now just the 1st tdc bin?!?!?! ####
// Target
 	for(Int_t ii = 0; ii<nBars; ++ii){
        if( ( ADC_High_TARGET[ii] > MinCutT && ADC_High_TARGET[ii] <  HGcutT )){
             // found a hit - fill histo etc.
            hadcT[ii]->Fill(ADC_High_TARGET[ii]);
            nhits++;
        }
	// now look for "valid" tdc times - collect tdc info where there is a "good" adc
        if(ADC_High_TARGET[ii] > HGcutT && ADC_High_TARGET[ii] < 4000){
            hTARG_HITS->Fill(ii);
            htdcT[ii]->Fill(TDC_LE_TARGET[ii][0]);
            nhits_t++;
        }

        if(ADC_Low_TARGET[ii] < LGcutT){
            hadcLT[ii]->Fill(ADC_Low_TARGET[ii]);
            nhitsL++;
        }
    }
    // SFT
 	for(Int_t ii = 0; ii<nFibrs; ++ii){
        if( ( ADC_High_SFT[ii] > MinCutF && ADC_High_SFT[ii] <  HGcutF )){
             // found a hit - fill histo etc.
            hadcF[ii]->Fill(ADC_High_SFT[ii]);
            nhitsF++;
        }
        if(ADC_High_SFT[ii] > HGcutF && ADC_High_SFT[ii] < 4000){
            hSFT_HITS->Fill(ii);
            htdcF[ii]->Fill(TDC_LE_SFT[ii][0]);
            nhits_Ft++;
        }

        if(ADC_Low_SFT[ii] < LGcutF){
            hadcLF[ii]->Fill(ADC_Low_SFT[ii]);
            nhitsLF++;
        }
    }
    if(nhits > 0 ){
    hTARGvsSFT->Fill(25);
    }
    if(nhitsF > 0 ){
    hTARGvsSFT->Fill(-25);
    }
    return;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Find Peds routine
void E36FindPeds::TFindPeds() {
    // Loop over Target histograms.
    for(Int_t ii = 0; ii<nBars; ++ii){
        Double_t H_Mean = hadcT[ii]->GetMaximumBin(); //GetMean();
        Int_t H_LowCut = H_Mean - 100.0;
        Int_t H_HighCut = H_Mean + 100.0;
        hadcT[ii]->GetXaxis()->SetRangeUser(H_LowCut,H_HighCut);
        Int_t adcH = hadcT[ii]->GetEntries();
        if(adcH > T_minhighadc){
        TF1 *fitH = new TF1("fitH", "gaus", H_LowCut,H_HighCut);
        hadcT[ii]->Fit("fitH","R");
        adcPeakT[ii] = fitH->GetParameter(1);
        adcRMST[ii] = fitH->GetParameter(2);
   //cout << "bar  " << ii << " Ok H " << adcPeakT[ii] << "  " << adcRMST[ii] << endl;
        hMeanBarT->Fill(adcPeakT[ii], ii);
        hRMSbarT->Fill(adcRMST[ii], ii);
        }
	else {
	adcPeakT[ii] = H_Mean; 
        adcRMST[ii] = hadcT[ii]->GetRMS();
        hMeanBarT->Fill(adcPeakT[ii], ii);
        hRMSbarT->Fill(adcRMST[ii], ii);
	}
        Double_t L_Mean = hadcLT[ii]->GetMaximumBin(); //GetMean();
        Int_t L_LowCut = L_Mean - 100.0;
        Int_t L_HighCut = L_Mean + 100.0;
        hadcLT[ii]->GetXaxis()->SetRangeUser(L_LowCut,L_HighCut);
            Int_t adcL = hadcLT[ii]->GetEntries();
            if(adcL > T_minlowadc){
        TF1 *fitL = new TF1("fitL", "gaus", L_LowCut,L_HighCut);
        hadcLT[ii]->Fit("fitL","R");
        ladcRMST[ii] = fitL->GetParameter(2);
        ladcPeakT[ii] =  fitL->GetParameter(1);
   //cout << "bar  " << ii << " Ok L " << ladcPeakT[ii] << "  " << ladcRMST[ii] << endl;
        hMeanBarTL->Fill(ladcPeakT[ii], ii);
        hRMSbarTL->Fill(ladcRMST[ii], ii);
            }
	else {
	ladcPeakT[ii] = L_Mean; 
        ladcRMST[ii] = hadcLT[ii]->GetRMS();
        hMeanBarTL->Fill(ladcPeakT[ii], ii);
        hRMSbarTL->Fill(ladcRMST[ii], ii);
	}
        htdcT[ii]->GetXaxis()->SetRangeUser(T_LowMean,T_HighMean);
        Double_t T_Mean = htdcT[ii]->GetMaximumBin(); //GetMean();
        Int_t T_LowCut = T_Mean - 100.0;
        Int_t T_HighCut = T_Mean + 100.0;
        htdcT[ii]->GetXaxis()->SetRangeUser(T_LowCut,T_HighCut);
        Int_t tdcC = htdcT[ii]->GetEntries();
        if(tdcC > T_mintdc){
        TF1 *fitT = new TF1("fitT", "gaus", T_LowCut,T_HighCut);
        htdcT[ii]->Fit("fitT","R");
        tdcRMST[ii]= fitT->GetParameter(2);
        tdcPeakT[ii]= fitT->GetParameter(1);
   cout << "bar  " << ii << " Ok T " << tdcPeakT[ii] << "  " << tdcRMST[ii] << endl;
   cout << "bar  " << ii << " Ok T " << tdcC << " " << T_LowCut << " " << T_HighCut << endl;
        hMeanTDCT->Fill(tdcPeakT[ii], ii);
        hRMSTDCT->Fill(tdcRMST[ii], ii);
        }
	else {
	tdcPeakT[ii] = T_Mean; 
        tdcRMST[ii] = htdcT[ii]->GetRMS();
        hMeanTDCT->Fill(tdcPeakT[ii], ii);
        hRMSTDCT->Fill(tdcRMST[ii], ii);
	}
	BarHitsT[ii] = hTARG_HITS->GetBinContent(ii+1);
        if(BarHitsT[ii] < MinHits){
            hMissBarT->Fill(ii,10);
        }
    }
    // Loop over SFT histograms.
    for(Int_t ii = 0; ii<nFibrs; ++ii){
        Double_t HF_Mean = hadcF[ii]->GetMaximumBin(); //GetMean();
        Int_t H_LowCut = HF_Mean - 100.0;
        Int_t H_HighCut = HF_Mean + 100.0;
        hadcF[ii]->GetXaxis()->SetRangeUser(H_LowCut,H_HighCut);
        Int_t adcH = hadcF[ii]->GetEntries();
        if(adcH > F_minhighadc){
        TF1 *fitH = new TF1("fitH", "gaus", H_LowCut,H_HighCut);
        hadcF[ii]->Fit("fitH","R");
  cout << "fibre  " << ii << " Ok H " << adcH << endl;
        adcPeakF[ii] = fitH->GetParameter(1);
        adcRMSF[ii] = fitH->GetParameter(2);
        hMeanBarF->Fill(adcPeakF[ii], ii);
        hRMSbarF->Fill(adcRMSF[ii], ii);
        }
	else {
	adcPeakF[ii] = HF_Mean; 
        adcRMSF[ii] = hadcF[ii]->GetRMS();
        hMeanBarF->Fill(adcPeakF[ii], ii);
        hRMSbarF->Fill(adcRMSF[ii], ii);
	}
        Double_t LF_Mean = hadcLF[ii]->GetMaximumBin(); //GetMean();
        Int_t L_LowCut = LF_Mean - 100.0;
        Int_t L_HighCut = LF_Mean + 100.0;
        hadcLF[ii]->GetXaxis()->SetRangeUser(L_LowCut,L_HighCut);
            Int_t adcL = hadcLF[ii]->GetEntries();
            if(adcL > F_minlowadc){
        TF1 *fitL = new TF1("fitL", "gaus", L_LowCut,L_HighCut);
        hadcLF[ii]->Fit("fitL","R");
  cout << "fibre  " << ii << " Ok L " << adcL << endl;
        ladcRMSF[ii] = fitL->GetParameter(2);
        ladcPeakF[ii] =  fitL->GetParameter(1);
        hMeanBarFL->Fill(ladcPeakF[ii], ii);
        hRMSbarFL->Fill(ladcRMSF[ii], ii);
            }
	else {
	ladcPeakF[ii] = LF_Mean; 
        ladcRMSF[ii] = hadcLF[ii]->GetRMS();
        hMeanBarFL->Fill(ladcPeakF[ii], ii);
        hRMSbarFL->Fill(ladcRMSF[ii], ii);
	}
        htdcF[ii]->GetXaxis()->SetRangeUser(F_LowMean,F_HighMean);
        Double_t Ft_Mean = htdcF[ii]->GetMaximumBin(); //GetMean();
        Int_t T_LowCut = Ft_Mean - 100.0;
        Int_t T_HighCut = Ft_Mean + 100.0;
        htdcF[ii]->GetXaxis()->SetRangeUser(T_LowCut,T_HighCut);
        Int_t tdcC = htdcF[ii]->GetEntries();
  cout << "fibre  " << ii << " Ok T " << tdcC << " " << T_LowCut << " " << T_HighCut << endl;
        if(tdcC > F_mintdc){
        TF1 *fitT = new TF1("fitT", "gaus", T_LowCut,T_HighCut);
        htdcF[ii]->Fit("fitT","R");
  //cout << "fibre  " << ii << " Ok T " << tdcC << " " << T_LowCut << " " << T_HighCut << endl;
        tdcRMSF[ii]= fitT->GetParameter(2);
        tdcPeakF[ii]= fitT->GetParameter(1);
        hMeanTDCF->Fill(tdcPeakT[ii], ii);
        hRMSTDCF->Fill(tdcRMST[ii], ii);
        }
	else {
	tdcPeakF[ii] = Ft_Mean; 
        tdcRMSF[ii] = htdcF[ii]->GetRMS();
        hMeanTDCF->Fill(tdcPeakT[ii], ii);
        hRMSTDCF->Fill(tdcRMST[ii], ii);
	}
        BarHitsF[ii] = hSFT_HITS->GetBinContent(ii+1);
        if(BarHitsF[ii] < MinHits){
            hMissBarF->Fill(ii,10);
        }
    }
}
#endif // #ifdef E36FindPeds_cxx
