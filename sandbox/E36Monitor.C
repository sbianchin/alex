#define E36Monitor_cxx
//////////////////////////////////////////////////////////
// E36Monitor
// Target Analysis (DRG)
/////////////////////////////////////////////////////////

#include <stdio.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<string.h>
#include<cmath>
#include "E36Monitor.h"
#include "TF1.h"
#include "TH2.h"
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "ANAPATH.h"
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void E36Monitor::Analyser(){
// Set up Histos file and trees
    //OpenHistoFile();

     hfout = new TFile(histo_savefile,"recreate");

//Set up Histos
   E36Monitor::TMakeHistos();

// Now go to data.
    printf("Going to Data \n");
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
    Int_t NENTRY = nentries;
    RunEvents[file_num]= 0;
    printf("Nentries %d \n", NENTRY);
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
// Find the hits
   E36Monitor::TFindHits();
}
return;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t InOutFiles(char* savename, Int_t RunNum){
    sprintf(filename,"/data/trek/E36/Data/June_2015/root/Merged/Run%dMS.root",RunNum);
    cout << "Data File  "  << filename << endl;	
    if(std::ifstream(filename)) {
    cout << filename << " FOUND!! " << endl;
    }
    else{
    cout << filename << "  NOT FOUND!! " << endl;
    NoFile = -1;
    return -1;
    }
//    cout << "Input run  " << filename << endl;
//Set up histos save file name
    sprintf(histo_savefile,"/data/trek/E36/Data/June_2015/root/%s%dOut.root",savename,RunNum);
    cout << "Histogram Save File  "  << histo_savefile << endl;	
    return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Int_t MonitorValues() {
// Loop over Target histograms.
    for(Int_t ii = 0; ii<nBars; ++ii){
    hadcT[ii]->GetXaxis()->SetRangeUser(HGpedT[ii],4095);
        adcPeakT[ii][file_num] = hadcT[ii]->GetMean();
        adcRMST[ii][file_num] = hadcT[ii]->GetRMS();
	Int_t PeakMean = adcPeakT[ii][file_num];
        hMeanBarT->Fill(PeakMean,ii);
        hRMSbarT->Fill(adcRMST[ii][file_num], ii);
        hadcLT[ii]->GetXaxis()->SetRangeUser(LGpedT[ii],4095);
        ladcRMST[ii][file_num]=hadcLT[ii]->GetRMS();
        ladcPeakT[ii][file_num]=hadcLT[ii]->GetMean();
        hMeanBarTL->Fill(ladcPeakT[ii][file_num], ii);
        hRMSbarTL->Fill(ladcRMST[ii][file_num], ii);
        htdcT[ii]->GetXaxis()->SetRangeUser(500,1000);
        tdcRMST[ii][file_num]= htdcT[ii]->GetRMS();
        tdcPeakT[ii][file_num]= htdcT[ii]->GetMean();
        hMeanTDCT->Fill(tdcPeakT[ii][file_num], ii);
        hRMSTDCT->Fill(tdcRMST[ii][file_num], ii);
        BarHitsT[ii][file_num] = hTARG_HITS->GetBinContent(ii+1);
        if(BarHitsT[ii][file_num] < MinHits[0]){
            hMissBarT->Fill(ii,10);
        }
    }
     KxPos[file_num]=hK_Xhits->GetMean();
     KxWid[file_num]=hK_Xhits->GetRMS();
     KyPos[file_num]=hK_Yhits->GetMean();
     KyWid[file_num]=hK_Yhits->GetRMS();

// Loop over SFT histograms.
    for(Int_t ii = 0; ii<nFibrs; ++ii){
    //hadcF[ii]->GetXaxis()->SetRangeUser(HGpedF[ii],4095);
        adcPeakF[ii][file_num] = hadcF[ii]->GetMean();
        adcRMSF[ii][file_num] = hadcF[ii]->GetRMS();
        hMeanBarF->Fill(adcPeakF[ii][file_num], ii);
        hRMSbarF->Fill(adcRMSF[ii][file_num], ii);
        hadcLF[ii]->GetXaxis()->SetRangeUser(LGpedF[ii],4095);
        ladcRMSF[ii][file_num]=hadcLF[ii]->GetRMS();
        ladcPeakF[ii][file_num]=hadcLF[ii]->GetMean();
        hMeanBarFL->Fill(ladcPeakF[ii][file_num], ii);
        hRMSbarFL->Fill(ladcRMSF[ii][file_num], ii);
        htdcF[ii]->GetXaxis()->SetRangeUser(500,1000);
        tdcRMSF[ii][file_num]= htdcF[ii]->GetRMS();
        tdcPeakF[ii][file_num]= htdcF[ii]->GetMean();
        hMeanTDCF->Fill(tdcPeakF[ii][file_num], ii);
        hRMSTDCF->Fill(tdcRMSF[ii][file_num], ii);
 
        BarHitsF[ii][file_num] = hSFT_HITS->GetBinContent(ii+1);
        if(BarHitsF[ii][file_num] < MinHits[0]){
            hMissBarF->Fill(ii,10);
        }
    }
     hfout->Write();
   printf("Closing Time %d\n", file_num);
    return 0;
	}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t MonitorSummary(Int_t runnum){
	char warnings[100];
	sprintf(warnings,"/data/trek/E36/Data/June_2015/root/TargetSFTWarnings%d.dat",runnum);
       ofstream dout;
        dout.open (warnings);
        RunRatio = RunEvents[0]/RunEvents[1];  // NOTE: this is from count of bars only!!
        printf("Events %f %f %f \n", RunEvents[0], RunEvents[1], RunRatio);
// Target loop
    dout << "R1Events " << RunEvents[0] << " R2Events " << RunEvents[1] << " Ratio " << RunRatio << endl;
    for(Int_t ii = 0; ii<nBars; ++ii){
        // look for adc shifts
        adcShiftT[ii] = adcPeakT[ii][0] - adcPeakT[ii][1];
	      if(fabs(adcShiftT[ii]) > MaxAdcShift[0] ){
	      dout << "Bar# " << ii << " ADC SHIFT " << adcShiftT[ii] << endl;
	      }
          // look for tdc shifts
	    tdcShiftT[ii] = tdcPeakT[ii][0] - tdcPeakT[ii][1];
	      if(fabs(tdcShiftT[ii]) > MaxTdcShift[0] ){
	      dout << "Bar# " << ii << " TDC SHIFT " << tdcShiftT[ii] << " " << MaxTdcShift << endl;
	      }
          // look for dead bars
	    if(BarHitsT[ii][0] < DeadBarIfLess[0] || BarHitsT[ii][1] < DeadBarIfLess[0] ){
	      dout << "Dead Bar# " << ii << endl;
	      }
          // look for "Problem Bars, i.e. has count rate in bar changed"
	    BarProblemT[ii] = BarHitsT[ii][0] - RunRatio*BarHitsT[ii][1];
	    if(fabs(BarProblemT[ii]) > MaxRateShift[0] ){
	    dout << "Bad Bar# " << ii << " Rate Change " << BarProblemT[ii] << endl;
	    }
        }
// Beam problems
    BeamShiftX = KxPos[0] - KxPos[1];
    BeamShiftY = KyPos[0] - KyPos[1];
    if(fabs(BeamShiftX) > MaxBeamShift[0] || fabs(BeamShiftY) > MaxBeamShift[0] ){
        dout << "Beam Shifted X " << BeamShiftX << " Y " << BeamShiftY << endl;
    }
//SFT loop
      for(Int_t ii = 0; ii<nFibrs; ++ii){
        // look for adc shifts
        adcShiftF[ii] = adcPeakF[ii][0] - adcPeakF[ii][1];
	      if(fabs(adcShiftF[ii]) > MaxAdcShift[1] ){
	      dout << "SFT# " << ii << " ADC SHIFT " << adcShiftF[ii] << endl;
	      }
          // look for tdc shifts
	    tdcShiftF[ii] = tdcPeakF[ii][0] - tdcPeakF[ii][1];
	      if(fabs(tdcShiftF[ii]) > MaxTdcShift[1] ){
	      dout << "SFT# " << ii << " TDC SHIFT " << tdcShiftF[ii] << endl;
	      }
          // look for dead fibres
	    if(BarHitsF[ii][0] < DeadBarIfLess[1] || BarHitsF[ii][1] < DeadBarIfLess[1] ){
	      dout << "Dead SFT# " << ii << endl;
	      }
          // look for "Problem Fibres, i.e. has count rate in fibre changed"
	    BarProblemF[ii] = BarHitsF[ii][0] - RunRatio*BarHitsF[ii][1];
	    if(fabs(BarProblemF[ii]) > MaxRateShift[1] ){
	    dout << "Bad SFT# " << ii << " Rate Change " << BarProblemF[ii] << endl;
	    }
        }
	return 0;
    }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t PedsIn(char* TargFile, char* SFTFile){
    char PedfT[100];
    sprintf(PedfT,"/data/trek/E36/Data/June_2015/root/%s",TargFile);
    cout << "" << endl;
    cout << "Target Pedestals File:   " << PedfT << endl;
    cout << "" << endl;
    if(std::ifstream(PedfT)) {
    cout << PedfT << " FOUND!! " << endl;
    }
    else{
    cout << PedfT << "  NOT FOUND!! " << endl;
    NoFile = -1;
    return -1;
    }
   ifstream pedTdat(PedfT,ios::in);
    Int_t ii = 0;
    Int_t ij = 0;
    while(!pedTdat.eof()){
        pedTdat >> ij >> HGpedT[ii] >> LGpedT[ii] >> TcutLT[ii] >> TcutHT[ii];
        ii++;
    }
    char PedfF[100];
    sprintf(PedfF,"/data/trek/E36/Data/June_2015/root/%s",SFTFile);
    cout << "" << endl;
    cout << "SFT Pedestals File:   " << PedfF << endl;
    cout << "" << endl;
    if(std::ifstream(PedfF)) {
    cout << PedfF << " FOUND!! " << endl;
    }
    else{
    cout << PedfF << "  NOT FOUND!! " << endl;
    NoFile = -1;
    return -1;
    }
    ifstream pedFdat(PedfF,ios::in);
    Int_t il = 0;
    Int_t ik = 0;
    while(!pedFdat.eof()){
        pedFdat >> ik >> HGpedF[il] >> LGpedF[il] >> TcutLF[il] >> TcutHF[il];
        il++;
    }
    return 0;
    }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t Limits(){
    char Lfile[] = "Limits";
    char limitsfile[100];
sprintf(limitsfile,"/data/trek/E36/Data/June_2015/root/%s.dat",Lfile);
cout << "" << endl;
cout << "Limits File:   " << limitsfile << endl;
cout << "" << endl;
if(std::ifstream(limitsfile)) {
    cout << limitsfile << " FOUND!! " << endl;
}
else{
    cout << limitsfile << "  NOT FOUND!! " << endl;
    NoFile = -1;
    return -1;
}
ifstream limitsdat(limitsfile,ios::in);
    Int_t il = 0;
    while(!limitsdat.eof()){
    limitsdat >> MaxAdcShift[il] >> MaxTdcShift[il] >> MaxRateShift[il] >> MinHits[il] >> DeadBarIfLess[il] >> MaxBeamShift[il];
        if(il==0){
            cout << "Target " << endl;
        }
        if(il==1){
            cout << "SFT " << endl;
        }
        cout << " MaxAdcShift " << "MaxTdcShift " << "MaxRateShift " << "MinHits " << "DeadBarIfLess " << "MaxBeamShift " << endl;
    cout << MaxAdcShift[il] << " " << MaxTdcShift[il] << " " << MaxRateShift[il] << " " << MinHits[il] << " " << DeadBarIfLess[il] << " " << MaxBeamShift[il] << endl;
        il++;
    }
/*  MaxAdcShift If ADC mean shifts more than this number raise warning.
    MaxTdcShift  If TDC mean shifts more than this number raise warning.
    MaxRateShift If Bar count rate shifts more than this number raise warning.
    MinHits  bar histoed if nhits not greater than this number.
    DeadBarIfLess  If Bar counts less than this declare dead and raise warning
    MaxBeamShift  If Kaon beam shiftsin either x or y raise warning.
    */
    return 0;
    }
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // NoBars routine/info is not used but may be in the future.
Int_t NoBars(char* nb_file){
    char NoBarf[100];
    sprintf(NoBarf,"%s",nb_file);
    cout << "" << endl;
    cout << "No Bar File:   " << nb_file << endl;
    cout << "" << endl;
    ifstream bdat(NoBarf,ios::in);
    Int_t ii = 0;
    while(!bdat.eof()){
        bdat >> NoBar[ii];
        cout << "NB   " << ii << "    " << NoBar[ii] << endl;
        ii++;
    }
    nNoBars = ii-1;
    cout << "nNB   " << nNoBars << endl;
    bdat.close();
    return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// .L E36Monitor.C then use as follows to start analysis:
// First file is treated as "IDEAL" run and second file is compared to it!!
Int_t DataCompare(Int_t Rnum0, Int_t Rnum1) {
    NoFile = 1;
// Get the Pedestal info
    char Tfile[] = "PedestalsT";
    char Ffile[] = "PedestalsF";
    char targpedsfile0[100];
    char sftpedsfile0[100];
    sprintf(targpedsfile0,"%s%d.dat",Tfile,Rnum0);
    cout << "Target Pedestal File " << targpedsfile0 << endl;
    sprintf(sftpedsfile0,"%s%d.dat",Ffile,Rnum0);
    cout << "SFT Pedestal File " << sftpedsfile0 << endl;
    PedsIn(targpedsfile0,sftpedsfile0);
    if(NoFile < 0) {
    return 1;
    }
    Limits(); // get limits before doing comparison
    if(NoFile < 0) {
        return 1;
    }
    char save_name[] = "Monitor";
    file_num = 0;  // analyse 1st run
    InOutFiles(save_name, Rnum0);
    if(NoFile < 0) {
    return 1;
    }
    E36Monitor means0;
    means0.Analyser();
    MonitorValues(); // extract the 1st run values of the parameters being monitored 
    NoFile = 1;
    file_num = 1;  // analyse 2nd run
    InOutFiles(save_name, Rnum1);
    if(NoFile < 0) {
    return 1;
    }
    // Get 2nd file Pedestal info
    //char Tfile[] = "PedestalsT";
    //char Ffile[] = "PedestalsF";
    char targpedsfile1[100];
    char sftpedsfile1[100];
    sprintf(targpedsfile1,"%s%d.dat",Tfile,Rnum1);
    cout << "Target Pedestal File " << targpedsfile1 << endl;
    sprintf(sftpedsfile1,"%s%d.dat",Ffile,Rnum1);
    cout << "SFT Pedestal File " << sftpedsfile1 << endl;
    PedsIn(targpedsfile1,sftpedsfile1);
    if(NoFile < 0) {
        return 1;
    }
    E36Monitor means1;
    means1.Analyser();
    MonitorValues(); // extract the 2nd run values of the parameters being monitored
    MonitorSummary(Rnum1); // compare the two sets of parameters and output results.
   return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// .L E36Monitor.C then use as follows to start analysis:
// First file is treated as "IDEAL" run and second file is compared to it!!
Int_t CompareRuns(Int_t Rnum0, Int_t Rnum1) {
    NoFile = 1;
    Limits(); // get limits before doing comparison
    if(NoFile < 0) {
        return 1;
    }
    char ideal_name[] = "IdealRunvalues";
    char idealfile0[100];
    sprintf(idealfile0,"/data/trek/E36/Data/June_2015/root/%s%d.dat",ideal_name,Rnum0);
    if(std::ifstream(idealfile0)) {
        cout << idealfile0 << " FOUND!! " << endl;
    }
    else{
        cout << idealfile0 << "  NOT FOUND!! " << endl;
        NoFile = -1;
        return -1;
    }
    ifstream idealdat(idealfile0,ios::in);
    Int_t ij = 0;
    idealdat >> RunEvents[0] >> KxPos[0] >> KyPos[0];
    //cout << RunEvents[0] << " " << KxPos[0] << " " << KyPos[0] << endl;
    for(Int_t ii = 0; ii<nBars; ++ii){
        idealdat >> ij >> adcPeakT[ii][0] >> tdcPeakT[ii][0] >> BarHitsT[ii][0];
        //cout << ij << " " << adcPeakT[ii][0] << " " << tdcPeakT[ii][0] << " " << BarHitsT[ii][0] << endl;
    }
    ij = 0;
    for(Int_t ii = 0; ii<nFibrs; ++ii){
        idealdat >> ij >> adcPeakF[ii][0] >> tdcPeakF[ii][0] >> BarHitsF[ii][0];
        //cout << ii << " " << adcPeakF[ii][0] << " " << tdcPeakF[ii][0] << " " << BarHitsF[ii][0] << endl;
    }
    // Get the Pedestal info
    char Tfile[] = "PedestalsT";
    char Ffile[] = "PedestalsF";
    char targpedsfile1[100];
    char sftpedsfile1[100];
    sprintf(targpedsfile1,"%s%d.dat",Tfile,Rnum1);
    cout << "Target Pedestal File " << targpedsfile1 << endl;
    sprintf(sftpedsfile1,"%s%d.dat",Ffile,Rnum1);
    cout << "SFT Pedestal File " << sftpedsfile1 << endl;
    PedsIn(targpedsfile1,sftpedsfile1);
    if(NoFile < 0) {
        return 1;
    }
    char save_name[] = "Monitor";
    NoFile = 1;
    file_num = 1;  // analyse 2nd run
    InOutFiles(save_name, Rnum1);
    if(NoFile < 0) {
        return 1;
    }
    E36Monitor means1;
    means1.Analyser();
    MonitorValues(); // extract the 2nd run values of the parameters being monitored
    MonitorSummary(Rnum1); // compare the two sets of parameters and output results.
    return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t QuickLook(Int_t Rnum0) {
    NoFile = 1;
    char save_name[] = "QuickCk";
    file_num = 0;
    InOutFiles(save_name, Rnum0);
    if(NoFile < 0) {
    return 1;
    }
    for(Int_t ii = 0; ii < nBars; ii++){
        HGpedT[ii] = 100;
        LGpedT[ii] = 100;
        TcutLT[ii] = 100;
        TcutHT[ii] = 1020;
    }
    for(Int_t il = 0; il < nFibrs; il++){
        HGpedF[il] = 100;
        LGpedF[il] = 100;
        TcutLF[il] = 100;
        TcutHF[il] = 1020;
    }
    E36Monitor quick;
    quick.Analyser();
     hfout->Write();
   printf("Closing Time %d\n", file_num);
   return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t IdealRun(Int_t Rnum0) {
    NoFile = 1;
    char save_name[] = "IdealRun";
    // Get the Pedestal info
    char Tfile[] = "PedestalsT";
    char Ffile[] = "PedestalsF";
    char targpedsfile0[100];
    char sftpedsfile0[100];
    sprintf(targpedsfile0,"%s%d.dat",Tfile,Rnum0);
    cout << "Target Pedestal File " << targpedsfile0 << endl;
    sprintf(sftpedsfile0,"%s%d.dat",Ffile,Rnum0);
    cout << "SFT Pedestal File " << sftpedsfile0 << endl;
    PedsIn(targpedsfile0,sftpedsfile0);
    if(NoFile < 0) {
        return 1;
    }
    file_num = 0;
    InOutFiles(save_name, Rnum0);
    if(NoFile < 0) {
        return 1;
    }
    E36Monitor ideal;
    ideal.Analyser();
    MonitorValues(); // extract the 1st run values of the parameters being monitored
// Output for later use
    char idealvalues[100];
    sprintf(idealvalues,"/data/trek/E36/Data/June_2015/root/%svalues%d.dat",save_name,Rnum0);
    ofstream idealout;
    idealout.open (idealvalues);
    idealout << RunEvents[0] << " " << KxPos[0] << " " << KyPos[0] << endl;
    for(Int_t ii = 0; ii<nBars; ++ii){
        idealout << ii << " " << adcPeakT[ii][0] << " " << tdcPeakT[ii][0] << " " << BarHitsT[ii][0] << endl;
    }
    for(Int_t ii = 0; ii<nFibrs; ++ii){
        idealout << ii << " " << adcPeakF[ii][0] << " " << tdcPeakF[ii][0] << " " << BarHitsF[ii][0] << endl;
    }
   
    //hfout->Write();
    //printf("Closing Time %d\n", file_num);
    return 0;
}
