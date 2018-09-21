//////////////////////////////////////////////////////////
// E36Monitor
// Compares new run with "ideal" run
// outputs a warnings file for changes beyond "limits"
// DRG May 2015
//////////////////////////////////////////////////////////

#ifndef E36Monitor_h 
#define E36Monitor_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
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

//Some Global Variables
TFile *hfout;
// Monitoring limits &
Double_t  MaxAdcShift[2]={300,300}, MaxTdcShift[2]={5,5}, MaxRateShift[2]={100,100};Double_t MaxBeamShift[2] = {3.0,3.0};
Int_t MinHits[2] = {5,5}, DeadBarIfLess[2] = {5,5};
Double_t  RunRatio, RunEvents[2] = {0,0};
//Monitoring results
  //Target Monitor
Double_t adcPeakT[nBars][2], adcRMST[nBars][2];
Double_t ladcPeakT[nBars][2], ladcRMST[nBars][2];
Double_t tdcPeakT[nBars][2], tdcRMST[nBars][2];
Int_t BarHitsT[nBars][2], BarProblemT[nBars];
Double_t adcShiftT[nBars], tdcShiftT[nBars];
//Beam monitor
Double_t KxPos[2], KxWid[2], KyPos[2], KyWid[2];
Double_t BeamShiftX, BeamShiftY;
///SFT Monitor
Double_t adcPeakF[nFibrs][2], adcRMSF[nFibrs][2];
Double_t ladcPeakF[nFibrs][2], ladcRMSF[nFibrs][2];
Double_t tdcPeakF[nFibrs][2], tdcRMSF[nFibrs][2];
//Double_t adcHpedF[nFibrs], adcLpedF[nFibrs], tdcLcutF[nFibrs], tdcHcutF[nFibrs];
Int_t BarHitsF[nFibrs][2], BarProblemF[nFibrs];
Double_t adcShiftF[nBars], tdcShiftF[nBars];
// Temporary saves of tdc info
Int_t TDC_LE_Save[nBars];
// Keep ADC & TDC for nhits cut filling of histos
Int_t ADC_SAVE[nBars], TDC_SAVE[nBars], HitBar[nBars];
 	Int_t ip = 0;

// run_file(filename) is file containing input trees.
   char filename[100];
   char run_file;
//histos output root file
  char histo_savefile[100];
   Int_t file_num;
   Int_t NoFile = 1;

//TObjArray *Hlist(0); // create an array of Histograms
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
TH1F *hK_Xhits;
TH1F *hK_Yhits;
TH1F *htdcdepth; // histo of depth of tdc hit
TH1F *hMissBarT;
TH1F *hMissBarF;
TH2F *hMeanBarT; // high gain Monitor Target vs bar number
TH2F *hMeanBarF; // high gain Monitor SFT vs fibre number
TH2F *hRMSbarT; // hign gain RMS vs bar number
TH2F *hRMSbarF; // hign gain RMS vs fibre number
TH2F *hMeanBarTL; // low gain Monitor Target vs bar number
TH2F *hMeanBarFL; // low gain Monitor SFT vs fibre number
TH2F *hRMSbarTL; // low gain RMS vs bar number
TH2F *hRMSbarFL; // low gain RMS vs fibre number
TH2F *hMeanTDCT; // tdc Monitor Target vs bar number
TH2F *hMeanTDCF; // tdc Monitor SFT vs fibre number
TH2F *hRMSTDCT; // tdc RMS vs bar number
TH2F *hRMSTDCF; // tdc RMS vs fibre number*/
TH2F *hHGT_HITS; // Target High Gain Hit Positions
TH2F *hLGT_HITS; // Target Low Gain Hit Positions
TH2F *hK_HITS; // Target Kaon Hit Positions
TH2F *hK_HITS_LR; // Low Res Target Kaon Hit Positions

class E36Monitor {
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
    
   E36Monitor(TTree *tree=0);
   virtual ~E36Monitor();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Analyser();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void    TMakeHistos();
   virtual void    TFindHits();
};

#endif

#ifdef E36Monitor_cxx
E36Monitor::E36Monitor(TTree *tree) : fChain(0)
{
// tree file looked for here

   if (tree == 0) {
       printf("Run %s \n", filename);  
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
       if (!f || !f->IsOpen()) {
           f = new TFile(filename);
       }
       f->GetObject("Tree",tree);
   }
   Init(tree);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
E36Monitor::~E36Monitor()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t E36Monitor::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Long64_t E36Monitor::LoadTree(Long64_t entry)
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
void E36Monitor::Init(TTree *tree)
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
Bool_t E36Monitor::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void E36Monitor::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t E36Monitor::Cut(Long64_t entry)
{
// This function may be called from Analyser.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Make Histograms
 void E36Monitor::TMakeHistos() {
//set up histos Target
    hHGT_HITS = new TH2F("HITS_HG_T","Target HG Hit Positions",700,-35.,35.,700,-35.,35.);
    //hList->Add(hHGT_HITS);
    hLGT_HITS = new TH2F("HITS_LG_T","Target LG Hit Positions",700,-35.,35.,700,-35.,35.);
    hK_HITS = new TH2F("KaonHitPositions","Kaon Hit Positions",700,-35.,35.,700,-35.,35.);
    hK_HITS->SetMarkerStyle(2);
     hK_HITS_LR = new TH2F("LRKaonHitPositions","LR Kaon Hit Positions",70,-35.,35.,70,-35.,35.);
     hK_Xhits = new TH1F("KaonXHit","Kaon X Positions",700,-35.,35.);
     hK_Yhits = new TH1F("KaonYHit","Kaon Y Positions",700,-35.,35.);
    hTARGvsSFT = new TH1F("EvRatio","Event Ratio",200,-100,100);
    hTARG_HITS = new TH1F("BarsHit","Bars Hit",260,0,260);
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
       htdcT[ch] = new TH1F(ahist_namet.c_str(), ahist_namet.c_str(),1050,0,1050);
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
       htdcF[ch] = new TH1F(ahist_nameft.c_str(), ahist_nameft.c_str(),1050,0,1050);
       //-------------------------
   }
  return;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Find Hits routine
 void E36Monitor::TFindHits() {

// Set counters etc to 0
    nhits =0;
    nhits_t =0;
    nhitsL = 0;
    nhitsF =0;
    nhitsLF = 0;
    nhits_Ft =0;
     nKhits =0;
	ip = 0;
// loop over adc & tdc info to look for hits
    // #### For now just the 1st tdc bin?!?!?! ####
// Target
 	for(Int_t ii = 0; ii<nBars; ++ii){
	/*Int_t HGcut = HGped[ii];
	Int_t TDCL = TcutL[ii];
	Int_t TDCH = TcutH[ii];*/
    //printf(" Bar %d Adc %d Tdc %d \n", ii, ADC_High_TARGET[ii], TDC_LE_TARGET[ii][0]);
        if( ADC_High_TARGET[ii] >  HGpedT[ii] && TDC_LE_TARGET[ii][0] > TcutLT[ii] &&
	TDC_LE_TARGET[ii][0] < TcutHT[ii]){
             // found a hit - fill histo etc.
            hadcT[ii]->Fill(ADC_High_TARGET[ii]);
            hTARG_HITS->Fill(ii);
               Xhit[nhits] = Xloc[ii];
               Yhit[nhits] = Yloc[ii];
            HitBar[nhits] = ii;
            hHGT_HITS->Fill(Xhit[ii],Yhit[ii]);
            htdcT[ii]->Fill(TDC_LE_TARGET[ii][0]);
            nhits++;
            nhits_t++;
        }
        // search K hits
        if ( (ADC_High_TARGET[ii] > 3000) &&
             (TDC_LE_TARGET[ii][0] >  TcutLT[ii]  && TDC_LE_TARGET[ii][0] < TcutHT[ii] )) {
            XKhit[nKhits] = Xloc[ii];
            YKhit[nKhits] = Yloc[ii];
            hK_HITS->Fill(XKhit[nKhits],YKhit[nKhits]);
            hK_HITS_LR->Fill(XKhit[nKhits],YKhit[nKhits]);
            hK_Xhits->Fill(XKhit[nKhits]);
            hK_Yhits->Fill(YKhit[nKhits]);
            nKhits++;
        }
        // search low gain hits
        if(ADC_Low_TARGET[ii] > LGpedT[ii] && TDC_LE_TARGET[ii][0] > TcutLT[ii] &&
	TDC_LE_TARGET[ii][0] < TcutHT[ii]){
           hLGT_HITS->Fill(Xhit[ii],Yhit[ii]);
            hadcLT[ii]->Fill(ADC_Low_TARGET[ii]);
            nhitsL++;
        }
    }
    /// SFT
 	for(Int_t ii = 0; ii<nFibrs; ++ii){
        //printf(" Fiber %d Adc %d Tdc %d \n", ii, ADC_High_SFT[ii], TDC_LE_SFT[ii][0]);
        if( ( ADC_High_SFT[ii] > HGpedF[ii] && TDC_LE_SFT[ii][0] >  TcutLF[ii]  &&
	TDC_LE_SFT[ii][0] < TcutHF[ii] )){
             // found a hit - fill histo etc.
            hadcF[ii]->Fill(ADC_High_SFT[ii]);
            nhitsF++;
            hSFT_HITS->Fill(ii);
            htdcF[ii]->Fill(TDC_LE_SFT[ii][0]);
            nhits_Ft++;
        }

        if(ADC_Low_SFT[ii] > LGpedF[ii]){
            hadcLF[ii]->Fill(ADC_Low_SFT[ii]);
            nhitsLF++;
        }
    }
    if(nhits > 0 ){
    hTARGvsSFT->Fill(25);
    RunEvents[file_num]++;
    }
    if(nhitsF > 0 ){
    hTARGvsSFT->Fill(-25);
    }
   return;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif // #ifdef E36Monitor_cxx
