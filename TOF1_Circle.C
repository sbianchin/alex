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
//#include "TSpectrum.h"
#include "TMarker.h"
//#include "Event_Display_MS.h"
#include "ANAPATH.h"
#include "Thresholds.h"
#include "CommonParameters.h"
#endif

void TOF1_Circle()
{

	float Gap1[3][2];		float Gap7[3][2];
	float Gap2[3][2];		float Gap8[3][2];
	float Gap3[3][2];		float Gap9[3][2];
	float Gap4[3][2];		float Gap10[3][2];
	float Gap5[3][2];		float Gap11[3][2];
	float Gap6[3][2];		float Gap12[3][2];

	TH2F *h_Circle = new TH2F("Test", "TARGET and TOF1", 500, -50, 50, 500, -50, 50);
	TH2F *h_Target = new TH2F("Test Target", "Target", 3000, -50, 50, 3000, -50, 50);

	TEllipse *ell = new TEllipse(0, 0, 44.5,0);

	TCanvas *c1;
  	c1 = new TCanvas("cc1","cc1",0,200,500,500);
  	//c1->Divide(3,2);
  	c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

//  	TLine *Gap1l = new TLine(28.9952, 34.7881, 15.7548, 42.2881);
//  	TLine *Gap2l = new TLine(42.2881, 15.7548, 34.7881, 28.7452);
//  	TLine *Gap3l = new TLine(44.5, -7.5, 44.5, 7.5);
// 	    TLine *Gap4l = new TLine(42.2881, -15.7548, 34.7881, -28.7452);
//  	TLine *Gap5l = new TLine(28.7452, -34.7881, 15.7548, -42.2881);
//  	TLine *Gap6l = new TLine(-7.5, -44.5, 7.5, -44.5);
//  	TLine *Gap7l = new TLine(-28.7452, -34.7881, -15.7548, -42.2881);
//  	TLine *Gap8l = new TLine(-42.2881, -15.7548, -34.7881, -28.7452);
//  	TLine *Gap9l = new TLine(-44.5, -7.5, -44.5, 7.5);
//   	TLine *Gap10l = new TLine(-42.2881, 15.7548, -34.7881, 28.7452);
// 	    TLine *Gap11l = new TLine(-28.9952, 34.7881, -15.7548, 42.2881);
//  	TLine *Gap12l = new TLine(-7.5, 44.5, 7.5, 44.5);

  	for(int t=0; t<256; t++)
  	{
  		h_Target->Fill(Xloc[t],Yloc[t]);
  	}

	for(int i=0; i<3; i++)
	{
		Gap1[i][0] = TOF_Gap1XLoc[i];	Gap7[i][0] = TOF_Gap7XLoc[i];
		Gap2[i][0] = TOF_Gap2XLoc[i];	Gap8[i][0] = TOF_Gap8XLoc[i];
		Gap3[i][0] = TOF_Gap3XLoc[i];	Gap9[i][0] = TOF_Gap9XLoc[i];
		Gap4[i][0] = TOF_Gap4XLoc[i];	Gap10[i][0] = TOF_Gap10XLoc[i];
		Gap5[i][0] = TOF_Gap5XLoc[i];	Gap11[i][0] = TOF_Gap11XLoc[i];
		Gap6[i][0] = TOF_Gap6XLoc[i];	Gap12[i][0] = TOF_Gap12XLoc[i];
	}

	for(int j=0; j<3; j++)
	{
		Gap1[j][1] = TOF_Gap1YLoc[j];	Gap7[j][1] = TOF_Gap7YLoc[j];
		Gap2[j][1] = TOF_Gap2YLoc[j];	Gap8[j][1] = TOF_Gap8YLoc[j];
		Gap3[j][1] = TOF_Gap3YLoc[j];	Gap9[j][1] = TOF_Gap9YLoc[j];
		Gap4[j][1] = TOF_Gap4YLoc[j];	Gap10[j][1] = TOF_Gap10YLoc[j];
		Gap5[j][1] = TOF_Gap5YLoc[j];	Gap11[j][1] = TOF_Gap11YLoc[j];
		Gap6[j][1] = TOF_Gap6YLoc[j];	Gap12[j][1] = TOF_Gap12YLoc[j];
	}

	cout << "  " << endl;
	cout << "  " << endl;

	//for(int k=0; k<3; k++)
	//{
	//	cout << Gap12[k][0] << "  " << Gap12[k][1] << endl;
	//}

	for(int ii=0; ii<3; ii++)
	{
		h_Circle->Fill(Gap1[ii][0],Gap1[ii][1]);		h_Circle->Fill(Gap7[ii][0],Gap7[ii][1]);
		h_Circle->Fill(Gap2[ii][0],Gap2[ii][1]);		h_Circle->Fill(Gap8[ii][0],Gap8[ii][1]);
		h_Circle->Fill(Gap3[ii][0],Gap3[ii][1]);		h_Circle->Fill(Gap9[ii][0],Gap9[ii][1]);
		h_Circle->Fill(Gap4[ii][0],Gap4[ii][1]);		h_Circle->Fill(Gap10[ii][0],Gap10[ii][1]);
		h_Circle->Fill(Gap5[ii][0],Gap5[ii][1]);		h_Circle->Fill(Gap11[ii][0],Gap11[ii][1]);
		h_Circle->Fill(Gap6[ii][0],Gap6[ii][1]);		h_Circle->Fill(Gap12[ii][0],Gap12[ii][1]);

	}
		
	c1->cd();
	h_Circle->SetMarkerStyle(5);
	h_Circle->SetMarkerSize(1.2);
	h_Circle->SetLineWidth(2);
	h_Circle->Draw();
	h_Target->SetMarkerStyle(25);
	h_Target->SetMarkerColor(4);
	h_Target->SetMarkerSize(1.2);
	h_Target->Draw("same");
	Gap1l->Draw("same");
	//Gap2l->Draw("same");
	//Gap3l->Draw("same");
	//Gap4l->Draw("same");
	//Gap5l->Draw("same");
	//Gap6l->Draw("same");
	//Gap7l->Draw("same");
	//Gap8l->Draw("same");
	//Gap9l->Draw("same");
	//Gap10l->Draw("same");
	//Gap11l->Draw("same");
	Gap12l->Draw("same");
	ell->SetFillStyle(0);
	ell->SetLineColor(6);
	ell->Draw("same");

}