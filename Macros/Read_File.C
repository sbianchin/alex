#ifndef __CINT__
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
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
#include "TGaxis.h"
#include "TRandom.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TPolyLine.h"
#include "TLine.h"
#include "TArrow.h"
#include "TStyle.h"  
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TBranch.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TEllipse.h"
#include "TMarker.h"
#include "ANAPATH.h"
#endif

using namespace std;


void Read_File()
{

	gStyle->Clear();
  	TH1::AddDirectory(kFALSE);
  	gStyle->SetOptStat(0);

// set N_bins, x_min, x_max

	TH1F *h_var1 = new TH1F("#Delta#Phi", "#Phi_C2 - #Phi_Tgt", 44, -11, 11);
	TH1F *h_var2 = new TH1F("#DeltaZ_TOP", "Z_C2 - Z_SFT  |  TOP", 44, -11, 11);
	TH1F *h_var3 = new TH1F("#DeltaZ_BOTTOM", "Z_C2 - Z_SFT  |  BOTTOM", 44, -11, 11);

	double var1_Mean, var1_Sigma, var1_max;
	double var2_Mean, var2_Sigma, var2_max;
	double var3_Mean, var3_Sigma, var3_max;

	char var1_Mean_string[50], var1_Sigma_string[50];
	char var2_Mean_string[50], var2_Sigma_string[50];
	char var3_Mean_string[50], var3_Sigma_string[50];

	TLatex *tex_var1_Mean;
	TLatex *tex_var1_Sigma;

	TLatex *tex_var2_Mean;
	TLatex *tex_var2_Sigma;

	TLatex *tex_var3_Mean;
	TLatex *tex_var3_Sigma;


    ifstream inputFile("file_4143.txt");
    string line;

    int counter = 0;

    while (getline(inputFile, line))
    {
          istringstream ss(line);

        double var1, var2, var3;

        ss >> var1 >> var2 >> var3;
        //counter++;
        //cout << counter << "  " << var1 << "  " << var2 << "   " << var3 << endl;

        h_var1->Fill(var1);
        h_var2->Fill(var2);
        h_var3->Fill(var3);
    }

    TF1 *fit_line_var1;
    TF1 *fit_line_var2;
    TF1 *fit_line_var3;

    //set range for fit_line_var3

    TF1 *g1 = new TF1("g1", "gaus", -4,5);
    g1->SetLineColor(2);
    TF1 *g2 = new TF1("g2", "gaus", -4,4);
    g2->SetLineColor(2);
    TF1 *g3 = new TF1("g3", "gaus", -4,4);
    g3->SetLineColor(2);

    TCanvas *cc1;
    cc1 = new TCanvas("test","test",50,50,450,700);
    cc1->Divide(1,3);
    cc1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

    cc1->cd(1);
    h_var1->Fit(g1, "R");
    var1_max = h_var1->GetMaximum();
    var1_Mean = g1->GetParameter(1);
    var1_Sigma = g1->GetParameter(2);

	sprintf(var1_Mean_string, "#Delta#Phi = %2.3f", var1_Mean);
	sprintf(var1_Sigma_string, "#sigma = %2.3f", var1_Sigma);
	tex_var1_Mean = new TLatex(2.9,0.8*var1_max,var1_Mean_string);
	tex_var1_Sigma = new TLatex(3.6,0.8*var1_max-17,var1_Sigma_string);
	tex_var1_Mean->SetTextSize(0.04);
	tex_var1_Sigma->SetTextSize(0.04);

    h_var1->Draw();
	tex_var1_Mean->Draw();
	tex_var1_Sigma->Draw();

    cout << endl << endl;

    
    cc1->cd(2);
    h_var2->Fit(g2, "R");
    var2_max = h_var2->GetMaximum();
    var2_Mean = g2->GetParameter(1);
    var2_Sigma = g2->GetParameter(2);

	sprintf(var2_Mean_string, "#DeltaZ = %2.3f", var2_Mean);
	sprintf(var2_Sigma_string, "#sigma = %2.3f", var2_Sigma);
	tex_var2_Mean = new TLatex(2.9,0.8*var2_max,var2_Mean_string);
	tex_var2_Sigma = new TLatex(3.45,0.8*var2_max-17,var2_Sigma_string);
	tex_var2_Mean->SetTextSize(0.04);
	tex_var2_Sigma->SetTextSize(0.04);

	cout << "TEST:  " << var2_Mean_string << "  " << var2_Sigma_string << endl;

    h_var2->Draw();
	tex_var2_Mean->Draw();
	tex_var2_Sigma->Draw();

    cout << endl << endl;

    cc1->cd(3);
    h_var3->Fit(g3, "R");
    var3_max = h_var3->GetMaximum();
    var3_Mean = g3->GetParameter(1);
    var3_Sigma = g3->GetParameter(2);

	sprintf(var3_Mean_string, "#DeltaZ = %2.3f", var3_Mean);
	sprintf(var3_Sigma_string, "#sigma = %2.3f", var3_Sigma);
	tex_var3_Mean = new TLatex(2.9,0.8*var3_max,var3_Mean_string);
	tex_var3_Sigma = new TLatex(3.45,0.8*var3_max-17,var3_Sigma_string);
	tex_var3_Mean->SetTextSize(0.04);
	tex_var3_Sigma->SetTextSize(0.04);

    h_var3->Draw();
	tex_var3_Mean->Draw();
	tex_var3_Sigma->Draw();

	cout << "TEST   " << var1_max << "  " << var2_max << "  " << var3_max << endl;

}