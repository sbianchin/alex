#define E36FindPeds_cxx
//////////////////////////////////////////////////////////
// liear Track Fit for 96 Bars in 256.
// Uses tree that is combined from 3 X 32 channel boards. Combine.C!!
// (DRG)
/////////////////////////////////////////////////////////

#include <stdio.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<string.h>
#include<cmath>
#include "E36FindPeds.h"
#include "TF1.h"
#include "TH1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TFile.h"
#include "TTree.h"
#include "ANAPATH.h"

void E36FindPeds::Loop(Int_t FrstEv) {
    //Set up histos root file
    hfout = new TFile(histo_savefile,"recreate");
    
    //Set up Histos
    cout << "Make Histos" << endl;
    E36FindPeds::TMakeHistos();
    
    // Now go to data.
    printf("Going to Data \n");
    if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
	Int_t TotalEntries = nentries;
	cout << "Total Entries  " << TotalEntries << endl;
   Long64_t nbytes = 0, nb = 0;
   nEvG =0;
   for (Long64_t jentry=FrstEv; jentry<nentries;jentry++) {
       nEvent = jentry;
       //cout << "Event  " << nEvent << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

// Go find the hits and fill histos
           E36FindPeds::TFindHits();
   }
// Go find pedestal info
               E36FindPeds::TFindPeds();
    // Calculate pedestal cut positions
    // And output them and any WARNINGS!!
        ofstream poutt;
        poutt.open (Pedft);
        ofstream poutf;
        poutf.open (Pedff);
    Int_t HPedOutT, LPedOutT, LToutT, HToutT;
        for(Int_t ii = 0; ii<nBars; ++ii){
            // calculate for adc pedestals
            HGpedT[ii] = adcPeakT[ii] + NPsig*adcRMST[ii];
            HPedOutT = HGpedT[ii];
            LGpedT[ii] = ladcPeakT[ii] + NPsig*ladcRMST[ii];
            LPedOutT = LGpedT[ii];
            TcutLT[ii] = tdcPeakT[ii] - NTsig*tdcRMST[ii];
            LToutT = TcutLT[ii];
            TcutHT[ii] = tdcPeakT[ii] + NTsig*tdcRMST[ii];
            HToutT = TcutHT[ii];
            poutt << ii << " "<< HPedOutT << " " << LPedOutT << " " << LToutT << " " << HToutT << endl;
       }
    Int_t HPedOutF, LPedOutF, LToutF, HToutF;
        for(Int_t ii = 0; ii<nFibrs; ++ii){
            // calculate for adc pedestals
            HGpedF[ii] = adcPeakF[ii] + NPsig*adcRMSF[ii];
            HPedOutF = HGpedF[ii];
            LGpedF[ii] = ladcPeakF[ii] + NPsig*ladcRMSF[ii];
            LPedOutF = LGpedF[ii];
            TcutLF[ii] = tdcPeakF[ii] - NTsig*tdcRMSF[ii];
            LToutF = TcutLF[ii];
            TcutHF[ii] = tdcPeakF[ii] + NTsig*tdcRMSF[ii];
            HToutF = TcutHF[ii];
            poutf << ii << " "<< HPedOutF << " " << LPedOutF << " " << LToutF << " " << HToutF << endl;
        }
 printf("Close Out Time \n");
    //fout.close (doutfile);
    hfout->Write();

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t InOutFiles(Int_t RunNum){
// data input file.
    sprintf(filename,"%sRun%dMS.root",data_path,RunNum);
    if(std::ifstream(filename)) {
    cout << filename << " FOUND!! " << endl;
    }
    else{
    cout << filename << "  NOT FOUND!! " << endl;
    NoFile = -1;
    return -1;
    }
//Set up histos save file
    char save_name[] = "Pedestals";
    sprintf(histo_savefile,"%s%s%dOut.root",histo_save_path,save_name,RunNum);
    cout << "Histogram Save File  "  << histo_savefile << endl;	
// Set up save files for pedestals
    sprintf(Pedft,"%s%sT%d.dat",anal_save_path,save_name,RunNum);
    cout << "Target Pedestal File:   " << Pedft << endl;
    sprintf(Pedff,"%s%sF%d.dat",anal_save_path,save_name,RunNum);
    cout << "SFT Pedestal File:   " << Pedff << endl;
    return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t Cuts(Int_t HT, Int_t LT, Int_t minT, Int_t HF, Int_t LF, Int_t minF, Double_t NPW, Double_t NTW){
    HGcutT = HT;
    LGcutT = LT;
    MinCutT = minT;
    HGcutF = HF;
    LGcutF = LF;
    MinCutF = minF;
    NPsig = NPW;
    NTsig = NTW;
    cout << "Cuts changed " << HT << " " << LT << " " << minT << " " << HF << " " << LF << " " << minF << " " << NPW << " " << NTW << endl;
   return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// .L E36FindPeds.C then use as follows to start analysis:
Int_t GetPeds(Int_t Rnum) {
    NoFile = 1;
    InOutFiles(Rnum);
    if(NoFile < 0 ){
    return 1;
    }
    Int_t FirstEv = 1;
    E36FindPeds pedslook;
    pedslook.Loop(FirstEv);
    return 0;
}

