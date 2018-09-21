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
#include <typeinfo>
#include <dlib/optimization.h>
#include "TSystem.h"
#include "TFile.h"
#include "TProfile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TImage.h"
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
//#include "ANAPATH.h"
//#include "../../Macros/CommonParameters.h"
//#include "ADC_TARGET_Pedestals.h"
//#include "TDC_Windows.h"
//#include "Pedestals.h"
#include "intersect.cxx"
#include "GlobalParameters.h"
#include "G4DataRootApr19.h" //defines the variables in the .root files
#include "TOFsParameters.h"
#include "TargetParameters.h"
#include "G4Parameters.h"
#include "G4HistosApr.h"
typedef dlib::matrix<double,4,1> parameter_vector;
#include "Plot_Simulation_Event_Display.h"
#endif

using namespace std;

Double_t MinE = 0.03;
Double_t Tstep = 0.2;

//Defining some data structures
struct Lepton {
vector<double> vec_xx_lepton;
vector<double> vec_yy_lepton;
vector <double> vec_xx_lepton_rotate;
vector <double> vec_yy_lepton_rotate;
vector <double> vec_xx_kaon;
vector <double> vec_yy_kaon;

};



Lepton FindKaons(Int_t Run_Number, vector<double> vec_xx_lepton, vector<double> vec_yy_lepton, vector <double> vec_xx_lepton_rotate,
  vector <double> vec_yy_lepton_rotate, vector <double> vec_xx_kaon, vector <double> vec_yy_kaon){

  Double_t Emax[4] = {-1.0};
  Double_t tEmax[4] = {-1.0};
  Double_t Tmin[4] = {-1000.0};
  Double_t Tmax[4] = {-1000.0};
  Int_t idEmax[4] = {-1};
  Int_t idTmin[4] = {-1};
  Int_t idTmax[4] = {-1};
  // Find minimum time
  Tmin[0] = Thit[0];


  //Tmax[0] = Thit[nhits-1];
  idTmin[0] = 0;
  idTmax[0] = nhits - 1;


  for(Int_t ii = 0; ii<nhits; ++ii){
     //For each entry in Ehit,
     //ensure that the minimum T is in 0th position in Tmin
     //ensure that the maximum T is in 0th position in Tmax
      if((Ehit[ii] > MinE)) { // ignore "pedestal"
         // cout<<ii<<" Thit "<<Thit[ii]<<endl;
          if(Tmin[0] > Thit[ii]){
            //if there's a time that's smaller, change it
              Tmin[0] = Thit[ii];
              idTmin[0] = ii;
          }
          if(Thit[ii] > Tmax[0]){
            //add maximum time
              Tmax[0] = Thit[ii];
              idTmax[0] = ii;
          }
        //  cout<<ii<<" Tmin0 "<<Tmin[0]<<" "<<idTmin[0]<<endl;
        //  cout<<ii<<" Tmax0 "<<Tmax[0]<<" "<<idTmax[0]<<endl;
      }
  }
  //1st element in Tmin is 10 after 0th position of Tmin
  Tmin[1] = Tmin[0] + 10.0;

  for(Int_t ii = 0; ii<nhits; ++ii){
    //ensure that the 2nd smallest time is in 1st position in Tmin
    //ensure that the 2nd biggest time is in 1st position in Tmax
      if((Ehit[ii] > MinE) && (ii != idTmin[0])) { // ignore "Tmin[0]"
         // cout<<ii<<" Thit1 "<<Thit[ii]<<endl;
          if(Tmin[1] > Thit[ii]){
              Tmin[1] = Thit[ii];
              idTmin[1] = ii;
          }
          if((ii != idTmax[0])) { // ignore "Tmin[0]"
          if(Thit[ii] > Tmax[1]){
              Tmax[1] = Thit[ii];
              idTmax[1] = ii;
          }
          }
        //  cout<<ii<<" Tmin1 "<<Tmin[1]<<" "<<idTmin[1]<<endl;
        //  cout<<ii<<" Tmax1 "<<Tmax[1]<<" "<<idTmax[1]<<endl;
      }
  }
  // Assume max energy is a kaon hit!

  cout << "Assume max energy is a kaon hit!" << endl;
  cout << "Initially, the max energy is" << Emax[0] << endl;
  cout << "Initially, the time of max energy is" << tEmax[0] << endl;
  for(Int_t ii = 0; ii<nhits; ++ii){
      if((Ehit[ii] > MinE)) { // ignore "pedestal"
      //pedestal is already ignored in FindTrack
          if(Ehit[ii] > Emax[0]){
              Emax[0] = Ehit[ii];
              tEmax[0] = Thit[ii];
              idEmax[0] = ii;
          }
      }
  }
  cout << "Currently the max energy is" << Emax[0] << endl;
  cout << "Currently, the time of max energy is" << tEmax[0] << endl;
  cout << "Currently, the id of max energy is" << idEmax[0] << endl;

  for(Int_t ii = 0; ii<nhits; ++ii){
      if((Ehit[ii] > MinE) && (ii != idEmax[0])) { // ignore "pedestal" and id0
          if(Ehit[ii] > Emax[1]){
              Emax[1] = Ehit[ii];
              tEmax[1] = Thit[ii];
              idEmax[1] = ii;
          }
      }
  }
  cout << "Current, the second max energy is" << Emax[1] << endl;
  for (int i = 0; i != 4; i++) cout << endl << Tmin[i] ;
  for (int i = 0; i != 4; i++) cout << endl << Tmax[i];
  Double_t minMean = (Tmin[0] + Tmin[1])/2.0;
  Double_t maxMean = (Tmax[0] + Tmax[1])/2.0;
  Tstep = (maxMean - minMean)/2.0;
  //Do not consider events with low time separation for Run Number 524
  bool time_separation = true;
  if (to_string(Run_Number)[0] < 5 || Run_Number == 524) {
    if (Tstep < 2) {
      cout << "Error: Not enough time separation for lepton/kaon discrimination";
      time_separation = false;
    }
  }

  cout << endl << Tstep;

 // cout<<"MinTimes "<<Tmin[0]<<" "<<Tmin[1]<<" "<<minMean<<endl;
 // cout<<"MaxTimes "<<Tmax[0]<<" "<<Tmax[1]<<" "<<maxMean<<endl;
 // cout<<" Time step "<<Tstep<<endl;
  bool changeTstep = false;
  if(Tmin[1]  > (Tmin[0] + Tstep)){
      minMean = Tmin[0];
      changeTstep = true;
  }
  if(Tmax[0]  > (Tmax[1] + Tstep)){
      maxMean = Tmax[0];
      changeTstep = true;
  }
  if(changeTstep){
      Tstep = (maxMean - minMean)/2.0;
    //  cout<<"Min Max "<<minMean<<" "<<maxMean<<" Time step "<<Tstep<<endl;
  }

  cout << "changeTstep "<< changeTstep  << endl;
  cout << endl << Tstep << endl;

  //Discriminate between Lepton and Kaon
  for(Int_t ii = 0; ii<nhits; ++ii){

      if (time_separation) { //if there is large enough separation
        Double_t minDiff = Thit[ii] - minMean;
        //cout << minDiff;
       // cout<<"Time Diff "<<minDiff<<endl;
        if((Ehit[ii] > MinE)) { // ignore "pedestal"
            if((minDiff > Tstep) && (minDiff < 70.0)) {
              //cout << "if((minDiff > Tstep) && (minDiff < 70.0)) {" << endl;
              //printf("%f,%f", Xhit[ii],Yhit[ii]);
              //printf("%f,%f,%f", minDiff, Ehit[ii],Thit[ii]);
              //cout << endl;
               // cout<<"Target dE dt "<<ii<<" "<<Ehit[ii]<<" "<<Thit[ii]<<" "<<Xhit[ii]<<" "<<Yhit[ii]<<endl;

                vec_xx_lepton.push_back(Xhit[ii]);
                vec_yy_lepton.push_back(Yhit[ii]);
                vec_xx_lepton_rotate.push_back(Xloc[TARGET_Rotated_index[hit_bars[ii]]]);
                vec_yy_lepton_rotate.push_back(Yloc[TARGET_Rotated_index[hit_bars[ii]]]);
                //hXY->Fill(Xhit[ii],Yhit[ii]); //fill tracking (black triangles)
                //hTvsE->Fill(Ehit[ii],Thit[ii]);
            }
            if((minDiff > -1.0) && (minDiff <= Tstep)) {
              //cout << "if((minDiff > -1.0) && (minDiff <= Tstep)) {" << endl;
              //printf("%f,%f", Xhit[ii],Yhit[ii]);
              //printf("%f,%f,%f", minDiff, Ehit[ii],Thit[ii]);
              if (Run_Number == 5042) {
                vec_xx_lepton.push_back(Xhit[ii]); //In May 4, everything is a lepton
                vec_yy_lepton.push_back(Yhit[ii]);
                vec_xx_lepton_rotate.push_back(Xloc[TARGET_Rotated_index[hit_bars[ii]]]);
                vec_yy_lepton_rotate.push_back(Yloc[TARGET_Rotated_index[hit_bars[ii]]]);
              }
              else {
                vec_xx_kaon.push_back(Xhit[ii]);
                vec_yy_kaon.push_back(Yhit[ii]);
              }

              //cout << endl;
              //  cout<<"Target dE low dt "<<ii<<" "<<Ehit[ii]<<" "<<Thit[ii]<<" "<<Xhit[ii]<<" "<<Yhit[ii]<<endl;
                //hXYlowT->Fill(Xhit[ii],Yhit[ii]); //fill kaon (red triangles)
                //hTvsElowT->Fill(Ehit[ii],Thit[ii]);
            }
        }
      }
      else { //Add all the bars as leptons if there isn't
        vec_xx_lepton.push_back(Xhit[ii]);
        vec_yy_lepton.push_back(Yhit[ii]);
        vec_xx_lepton_rotate.push_back(Xloc[TARGET_Rotated_index[hit_bars[ii]]]);
        vec_yy_lepton_rotate.push_back(Yloc[TARGET_Rotated_index[hit_bars[ii]]]);
      }
  }

  cout << "This is in FindKaons" << endl;
  for (vector<double>::const_iterator i = vec_xx_kaon.begin(); i != vec_xx_kaon.end(); ++i)
    cout << *i << ' ';
  cout << endl;
  for (vector<double>::const_iterator i = vec_yy_kaon.begin(); i != vec_yy_kaon.end(); ++i)
    cout << *i << ' ';
  cout << "End of FindKaons" << endl;

  Lepton vector_lepton;
  vector_lepton.vec_xx_lepton = vec_xx_lepton;
  vector_lepton.vec_yy_lepton = vec_yy_lepton;
  vector_lepton.vec_xx_lepton_rotate = vec_xx_lepton_rotate;
  vector_lepton.vec_yy_lepton_rotate = vec_yy_lepton_rotate;
  vector_lepton.vec_xx_kaon = vec_xx_kaon;
  vector_lepton.vec_yy_kaon = vec_yy_kaon;
    return vector_lepton;
}





void Event_Display_5_1(Int_t Run_Number=5, Int_t ievt=0, int walk=0, Int_t enable_cout=0, Int_t display = 0){
  //write to file output


  int Switch=1; // Displays hit with no HG, but LG (0 = OFF ; 1 = ON)
  int Rotate=1; // When TOF1 is 12 or 6, or when TOF2 is 6 or 12, rotate by -90 deg to fit a horizontal line (0 = OFF ; 1 = ON)
  int T_limit = 3;
  Double_t MinE = 0.03;
  Double_t Tstep = 0.2;

  gStyle->Clear();
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(11);

  char Name_finput[200];
  sprintf(Name_finput,"G4Run%d.root", Run_Number);
  cout << "File opened:  " << Name_finput << endl;

  char footer[100];
  sprintf(footer,"Event_Display_Simulation.C  --  Run %d ; Event %d",Run_Number,ievt);

  char Version[100] = "Version 5.0";
  char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!

  TChain *fChain= new TChain("Kaon");
  fChain->Add(Name_finput);
  fChain->SetMakeClass(1);
  fChain->SetBranchAddress("n", &n);
  fChain->SetBranchAddress("s_hit", &s_hit);
  fChain->SetBranchAddress("x_hit", &x_hit);
  fChain->SetBranchAddress("tof1charge[5]", tof1charge);
  fChain->SetBranchAddress("tof1N[5]", tof1N);
  fChain->SetBranchAddress("px[10]", px);
  fChain->SetBranchAddress("py[10]", py);
  fChain->SetBranchAddress("pz[10]", pz);
  fChain->SetBranchAddress("u_x", &u_x);
  fChain->SetBranchAddress("u_y", &u_y);
  fChain->SetBranchAddress("u_z", &u_z);
  fChain->SetBranchAddress("genEvt", &genEvt);
  fChain->SetBranchAddress("genEne[10]", genEne);
  fChain->SetBranchAddress("cosTh", &cosTh);
  fChain->SetBranchAddress("sinTh", &sinTh);
  fChain->SetBranchAddress("Phi", &Phi);
  fChain->SetBranchAddress("dv_x", &dv_x);
  fChain->SetBranchAddress("dv_y", &dv_y);
  fChain->SetBranchAddress("dv_z", &dv_z);
  fChain->SetBranchAddress("targetE[5]", targetE);
  fChain->SetBranchAddress("targetdE[256]", targetdE);
  fChain->SetBranchAddress("targetdt[256]", targetdt);
  fChain->SetBranchAddress("targL", &targL);
  fChain->SetBranchAddress("t1pL", &t1pL);
  fChain->SetBranchAddress("tof1E[5]", tof1E);
  fChain->SetBranchAddress("tof1wpos[5][3]", tof1wpos);
  fChain->SetBranchAddress("tof1dt", &tof1dt);
  fChain->SetBranchAddress("secondP[5][100]", secondP);
  fChain->SetBranchAddress("t1", &t1);
  fChain->SetBranchAddress("t1px", &t1px);
  fChain->SetBranchAddress("t1py", &t1py);
  fChain->SetBranchAddress("t1pz", &t1pz);
  fChain->SetBranchAddress("SFTwpos1[5][3]", SFTwpos1);
  fChain->SetBranchAddress("SFTwpos2[5][3]", SFTwpos2);
  fChain->SetBranchAddress("SFTpx", &SFTpx);
  fChain->SetBranchAddress("SFTpy", &SFTpy);
  fChain->SetBranchAddress("SFTpz", &SFTpz);






  vector <double> vec_xx_lepton;           vec_xx_lepton.clear();
  vector <double> vec_yy_lepton;           vec_yy_lepton.clear();
  vector <double> vec_xx_lepton_test;      vec_xx_lepton_test.clear();
  vector <double> vec_yy_lepton_test;      vec_yy_lepton_test.clear();
  vector <double> vec_xx_lepton_rotate;    vec_xx_lepton_rotate.clear();
  vector <double> vec_yy_lepton_rotate;    vec_yy_lepton_rotate.clear();
  vector <double> vec_lepton_size;         vec_lepton_size.clear();
  vector <double> vec_lepton_rotate_size;  vec_lepton_rotate_size.clear();
  vector <int> vec_bar;                    vec_bar.clear();
  vector <int> vec_bar_rotate;             vec_bar_rotate.clear();
  vector <double> vec_yprime;              vec_yprime.clear();
  vector <double> vec_yprime_rotate;       vec_yprime_rotate.clear();
  vector <double> vec_Dy;                  vec_Dy.clear();
  vector <double> vec_Dy_rotate;           vec_Dy_rotate.clear();
  vector<double> vec_test;



  vector <double> vec_xx_kaon;           vec_xx_kaon.clear();
  vector <double> vec_yy_kaon;           vec_yy_kaon.clear();


  cout << "   " << endl;
  cout << "   " << endl;
  cout << "************************************************************************************************************" << endl;

  cout << "File opened:  " << Name_finput << endl;
  char run_string[100];
  sprintf(run_string,"Run %d ; Event %d",Run_Number,ievt);



  Int_t nentries = (Int_t)fChain->GetEntries();
  cout << "Total Number of Events:  " << nentries <<endl;
  //cout << "  " << endl;
  cout << "Event_Display.C -- " << Version << endl;
  cout << "************************************************************************************************************" << endl;
  cout << "  " << endl;
  cout << "  " << endl;
  //printf("Now Target \n");



  fChain->GetEntry(ievt);

  //Prints out Data from Event
  cout << "Here" << endl;
  cout << Phi << endl;
  //for (int i = 0; i != 5; i++) cout << tof1N[i] << " ";
  cout << "Tof1N" << endl;
  for (int i = 0; i != 5; i++) cout << tof1N[i] << " ";

  //for (int i = 0; i != 5; i++) cout << tof1wpos[i][0] << " ";
  //for (int i = 0; i != 5; i++) cout << tof1wpos[i][1] << " ";
  //for (int i = 0; i != 5; i++) cout << tof1wpos[i][2] << " ";
  cout << endl;
  for(Int_t ii = 0; ii<nBars; ++ii){
  if((targetdE[ii] > 0.0)) {
     //  cout<<"Target dE dt "<<ii<<" "<<targetdE[ii]<<" "<<targetdt[ii]<<" "<<Xloc[ii]<<" "<<Yloc[ii]<<endl;
      cout << ii << " ";

    }

      }
  cout << endl;
  for(Int_t ii = 0; ii<nBars; ++ii){
  if(targetdt[ii] > 0.0) {
     //  cout<<"Target dE dt "<<ii<<" "<<targetdE[ii]<<" "<<targetdt[ii]<<" "<<Xloc[ii]<<" "<<Yloc[ii]<<endl;
      cout << ii << " ";

    }
      }
  Xloc_TOF_flipped = -1*tof1wpos[0][1];
  Yloc_TOF_flipped = -1*tof1wpos[0][0];
  Yloc_TOF = tof1wpos[0][1];
  Xloc_TOF = tof1wpos[0][0];
  cout << "Xloc_TOF, Yloc_TOF" << tof1wpos[0][0] << ", " << Yloc_TOF << endl;

  double k_stop[2];
  k_stop[0] = dv_x;
  k_stop[1] = dv_y;
  AGAIN:
    Tmax = -1;
    Emax = -1;
    iMaxE = -1;
    iMaxT = -1;
    nhits =0;
    nLGhits =0;
    nKhits =0;
    Tdiffmax = 5;
    HGpedOffset = 5;
    LGpedOffset = 10;
    Double_t AngleTOF1 = 0.0;
    Int_t nTOF1s = 0;
    Int_t dTOF1s[5] = {-999};
    Int_t pTOF1s[5] = {-999};
    //for(Int_t id = 0; id<5; id++) cout << tof1N[id] << endl;
    //return;
    for(Int_t id = 0; id<5; id++){
        if(tof1N[id] >= 0){
            dTOF1s[nTOF1s] = id;
            pTOF1s[nTOF1s] = tof1N[id];
            nTOF1s++;
       cout<< " tof1N_n "<<id<<" "<<tof1N[id]<<endl;
        }

    }
    Int_t irt; //stores TOF1 Location
    for(Int_t it = 0; it <12; it++){
        for(Int_t id = 0; id<nTOF1s; id++){
            if(it == tof1N[id]){
                irt = tof1N[id];
                cout<<" TOF1 Location "<<irt<<" "<<TOF_Xloc[irt]<<" "<<TOF_Yloc[irt]<<endl;
            }
        }
    }
    // reset hit stuff
    for(Int_t ii = 0; ii<nBars; ++ii){
        Xhit[ii] = -1000.0;
        Yhit[ii] = -1000.0;
        Ehit[ii] = -1000.0;
        Thit[ii] = -1000.0;
    }
    //Find hits
        for(Int_t ii = 0; ii<nBars; ++ii){
        if((targetdE[ii] > MinE) && (targetdt[ii] > 0.01)) {
           //  cout<<"Target dE dt "<<ii<<" "<<targetdE[ii]<<" "<<targetdt[ii]<<" "<<Xloc[ii]<<" "<<Yloc[ii]<<endl;
            Xhit[nhits] = Xloc[ii];
            Yhit[nhits] = Yloc[ii];
            Ehit[nhits] = targetdE[ii];
            Thit[nhits] = targetdt[ii];
            nhits++;

          }
            }

    //cout << MinE << endl;
    //for (int i = 0; i != 256; i++) cout << targetdE[i] << " " ;
    //for (vector<double>::const_iterator i = vec_xx_lepton.begin(); i != vec_xx_lepton.end(); ++i)
    //  cout << *i << ' ';
    cout<<"Found nhits "<<nhits<<endl;

    //TCanvas *c2;
    c2 = new TCanvas("Event_Display.C  --  TARGET & SFT","Event_Display.C  --  TARGET & SFT",0,200,1050,350);
    //Divide the canvas
    c2->Divide(3,1);
    //c2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

    //Define Empty Plots and Axis


    Lepton vector_lepton_kaon = FindKaons(Run_Number, vec_xx_lepton, vec_yy_lepton, vec_xx_lepton_rotate, vec_yy_lepton_rotate, vec_xx_kaon, vec_yy_kaon);

    vec_xx_lepton = vector_lepton_kaon.vec_xx_lepton;
    vec_yy_lepton = vector_lepton_kaon.vec_yy_lepton;
    vec_xx_lepton_rotate = vector_lepton_kaon.vec_xx_lepton_rotate;
    vec_yy_lepton_rotate = vector_lepton_kaon.vec_yy_lepton_rotate;
    vec_xx_kaon = vector_lepton_kaon.vec_xx_kaon;
    vec_yy_kaon = vector_lepton_kaon.vec_yy_kaon;
    //for (vector<double>::const_iterator i = vec_xx_lepton.begin(); i != vec_xx_lepton.end(); ++i)
      //cout << *i << ' ';


    Graphs_lepton_kaon graphs_lepton_kaon = Prepare_Plot(Run_Number, ievt, vec_xx_lepton, vec_yy_lepton, vec_xx_lepton_rotate,
      vec_yy_lepton_rotate, vec_xx_kaon, vec_yy_kaon);
    gr3_Leptons_rotate = graphs_lepton_kaon.gr3_Leptons_rotate;
    gr3_Leptons = graphs_lepton_kaon.gr3_Leptons;
    gr_kaon = graphs_lepton_kaon.gr_kaon;
    a_fit_TDC_Gap_Fibers = graphs_lepton_kaon.a_fit_TDC_Gap_Fibers;
    b_fit_TDC_Gap_Fibers = graphs_lepton_kaon.b_fit_TDC_Gap_Fibers;
    a_gr_Leptons = graphs_lepton_kaon.a_gr_Leptons;
    b_gr_Leptons = graphs_lepton_kaon.b_gr_Leptons;
    a_fit_kaon = graphs_lepton_kaon.a_fit_kaon;
    b_fit_kaon = graphs_lepton_kaon.b_fit_kaon;
    /*
    Subplot_4(vec_xx_lepton, vec_yy_lepton, vec_xx_lepton_rotate,
      vec_yy_lepton_rotate, vec_xx_kaon, vec_yy_kaon);
    Subplot_5(vec_xx_lepton, vec_yy_lepton, vec_xx_lepton_rotate,
      vec_yy_lepton_rotate, vec_xx_kaon, vec_yy_kaon);
      */
      Subplot_1(Run_Number, ievt, gr3_Leptons_rotate, gr3_Leptons,
        a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers,
        Xloc_TOF, Yloc_TOF, Xloc_TOF_flipped, Yloc_TOF_flipped, dv_x, dv_y, irt);
      Subplot_2(gr_kaon, vec_xx_kaon, a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, irt);

      Subplot_3(gr3_Leptons_rotate, gr3_Leptons, gr_kaon,
        vec_xx_kaon, a_gr_Leptons, b_gr_Leptons,
      a_fit_TDC_Gap_Fibers, b_fit_TDC_Gap_Fibers, a_fit_kaon,
      b_fit_kaon, dv_x, dv_y, vec_test, 3);
    //Draw Plot 2

    //TImage *img = TImage::Create();
    //img->FromPad(c2);
    /*
    char image_name[100];
    sprintf(image_name, "./images/RUN_%d_event_%d.png", Run_Number, ievt);
    c2->SaveAs(image_name);
    */

  return;
}
