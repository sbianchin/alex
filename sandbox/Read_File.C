/*#ifndef __CINT__
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
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

//vector<double> par1;
//vector<double> par2;
//vector<double> par3;

double par_temp1 = 0; 
double par_temp2 = 0;
double par_temp3 = 0;

string evt;

ifstream finput;
finput.open("file_test.dat");

if(finput.fail()){
	cout << "Error: Could not read any file." << endl;
}
else{
	while(getline(finput,evt)){
		//ifstream fdat(finput,ios::in);
		//fdat >> par_temp1 >> par_temp2 >> par_temp2;
		//cout << "TEST:  " << par_temp1 << "  " << par_temp2 << "  " << par_temp3 << endl;
		cout << "TEST  " << evt << endl;
		ifstream fdat(finput, ios::in);
		fdat >> par_temp1 >> par_temp2 >> par_temp3;
		cout << "TEST2 :  " << par_temp1 << "  " << par_temp2 << "  " << par_temp3 << endl;
		fdat.close();
		}
}

    
finput.close();



} 
*/

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
  	gStyle->SetOptStat(1111);

	TH1F *h_var1 = new TH1F("#Delta#Phi", "#Phi_C2 - #Phi_Tgt", 22, -11, 11);
	TH1F *h_var2 = new TH1F("#DeltaZ_TOP", "Z_C2 - Z_SFT  |  TOP", 22, -11, 11);
	TH1F *h_var3 = new TH1F("#DeltaZ_BOTTOM", "Z_C2 - Z_SFT  |  BOTTOM", 22, -11, 11);


    ifstream inputFile("file_4143.dat");
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

    TCanvas *cc1;
    cc1 = new TCanvas("test","test",50,50,450,700);
    cc1->Divide(1,3);
    cc1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

    cc1->cd(1);
    h_var1->Fit("gaus");
    fit_line_var1 = h_var1->GetFunction("gaus");
    fit_line_var1->SetLineColor(2);
    h_var1->Draw();

    cout << endl << endl;

    cc1->cd(2);
    h_var2->Fit("gaus");
    fit_line_var2 = h_var2->GetFunction("gaus");
    fit_line_var2->SetLineColor(2);
    h_var2->Draw();

    cout << endl << endl;

    cc1->cd(3);
    h_var3->Fit("gaus");
    fit_line_var3 = h_var3->GetFunction("gaus");
    fit_line_var3->SetLineColor(2);
    h_var3->Draw();

}