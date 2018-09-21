#ifndef GLOBAL_PARAMETERS_H
#define GLOBAL_PARAMETERS_H

TFile *hfout;
//Some Global Variables
double PI = 3.14159265;
double Radion = 180.0 / PI;
char data_path[] = "/Users/davidgill/Documents/Trek/E36/Data/root/";
char anal_save_path[] = "Results/";
char histo_save_path[] = "RootFiles/";
char runinput_path[] = "../Data/root/";
const int nGaps = 12;
const int nBars = 256;
const int nFibrs = 128;
const int nEasiTDC = 8;
const int nTOF1 = 24;
const int nTOF1U = 12;
const int nTOF1D = 12;
const int nTOF2 = 48;
const int nTOF2AO = 12;
const int nTOF2AI = 12;
const int nTOF2BO = 12;
const int nTOF2BI = 12;
const int nPGCadc = 84;
const int nPGCtdc = 84;
const int nCk = 14;
const int nCpi = 14;
const int nGapVadc = 12;
const int nGapVtdc = 12;
const int nTTC = 12;
const int nHRTDC = 96;
const int nMWPCADC = 512;
const int nACU_ADC = 12;
const int nACD_ADC = 12;
const int nACUtdc = 12;
const int nACDtdc = 12;
const int nV792ADC = 192;
const int nV792TDC = 192;
const int nVT48TDC = 240;
const int nEvFlags = 40;
const int nTDCdepth = 16;
// Event number
int nEvent;
int RunNum[100];
int Rnum;
Int_t RunFileNum;
bool NoFile;  // Didn't find the file if false
int doRuns;
int TotalEvents = 0;
Int_t EvType = -1;
int FrstEv;
Int_t LastEv = 10;
Int_t Event_Out = 0;
Int_t Event_Num = 0;
int FirstRun;
int LastRun;
Int_t nBadevents = 0;
// Bools
bool FrstTime = true;
bool BadEvent = false;
// filename is (path)file containing input trees.
char filename[100];
//histos output root file
char histo_savefile[100];
char RunsFile[100];
// Depth of hit in multi hit TDCs
//Find TDC depths
Int_t TDmax;
Int_t iTargDepth[nBars];
Int_t iSFTDepth[nBars];
Int_t iCkDepth[nBars];
Int_t iCpiDepth[nBars];

#endif
