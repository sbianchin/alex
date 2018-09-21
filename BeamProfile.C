#ifndef __CINT__
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include "ANAPATH.h"
#include "E36Monitor.h"
#endif

void BeamProfile(Int_t run_number=1) {   

gStyle->SetPalette(1);

TH2F *BeamProfile;
TH1F *hK_Xhits;
TH1F *hK_Yhits;

char path_input[100];
sprintf(path_input,path_root);

char Name_finput[100];
sprintf(Name_finput,"/data/trek/E36/Data/June_2015/root/QuickCk%dOut.root",run_number);
cout << Name_finput << endl;
cout << "     " << endl;


char Name_Can[100];
char Title_Can[100];
sprintf(Name_Can,"Beam Profile  |  Run %d",run_number);
sprintf(Title_Can,"Beam Profile  |  Run %d",run_number);


TCanvas *c1;
c1 = new TCanvas(Name_Can,Title_Can,800,600); 
c1->Divide(2,2);
c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");


TFile *finput = new TFile(Name_finput);


//sprintf(Name_SFT_ADC_with_TDC_cut,"ADC High Gain (Ch. %d) -- SFT",j);
BeamProfile=(TH2F*)finput->FindObjectAny("LRKaonHitPositions");
hK_Xhits=(TH1F*)finput->FindObjectAny("KaonXHit");
hK_Yhits=(TH1F*)finput->FindObjectAny("KaonYHit");



c1->cd(1);
BeamProfile->Draw("Lego2");

c1->cd(2);
BeamProfile->Draw("colz");

c1->cd(3);
hK_Xhits->Draw();

c1->cd(4);
hK_Yhits->Draw();

//cout << "TEST" << endl;
		
/*

if(board_number==2){
	for(int j=0;j<64;j++){
		sprintf(Name_SFT_ADC_with_TDC_cut,"ADC High Gain (Ch. %d) -- SFT",j+64);

		hK_Xhits=(TH1F*)finput->FindObjectAny("KaonXHit");
	}
}



//////////////////////////////////////////////////////

if(board_number==1){
sprintf(Name_Can_SFT_ADC_with_TDC_cut,"ADC HG with TDC cut -- Run %d  Ch. 0 - 63  (SFT)  |  EASIROC Board %d",run_number,board_number);
sprintf(Title_Can_SFT_ADC_with_TDC_cut,"ADC HG with TDC cut -- Run %d  Ch. 0 - 63  (SFT)  |  EASIROC Board %d",run_number,board_number);
}

if(board_number==2){
sprintf(Name_Can_SFT_ADC_with_TDC_cut,"ADC HG with TDC cut -- Run %d  Ch. 64 - 127  (SFT)  |  EASIROC Board %d",run_number,board_number);
sprintf(Title_Can_SFT_ADC_with_TDC_cut,"ADC HG with TDC cut -- Run %d  Ch. 64 - 127  (SFT)  |  EASIROC Board %d",run_number,board_number);
}



if (board_number <= 2){

	TCanvas *c2;
	c2 = new TCanvas(Name_Can_SFT_ADC_with_TDC_cut,Name_Can_SFT_ADC_with_TDC_cut,1200,500); 
	c2->Divide(8,8);
	c2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

	for(Int_t ican=0; ican<64; ican++){
		if(board_number==1 && (ican==23 || ican==24 || ican==25 || ican==26 ||ican==27 || ican==28 || ican==29)) c2->cd(ican+1)->SetFillColor(45);
		if(board_number==2 && (ican==32 || ican==33 || ican==34 || ican==35 ||ican==36 || ican==37 || ican==38 ||  ican==39 || ican==52)) c2->cd(ican+1)->SetFillColor(45);
	 	c2->cd(ican+1)->SetLogy();
		h_SFT_ADC_with_TDC_cut[ican]->Draw();
	}
}

*/

//return;
}













