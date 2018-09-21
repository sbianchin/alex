#include <stdio.h>
#include<iostream>  //Header file for cin & cout
#include<fstream>
#include<sstream>
#include "TChain.h"
#include "TF1.h"
#include "TLine.h"
#include "TH2.h"
#include "TCanvas.h"
#include "GlobalParameters.h"
#include "G4DataRootApr19.h" //defines the variables in the .root files
#include "TOFsParameters.h"
#include "TargetParameters.h"
#include "G4Parameters.h"
#include "G4HistosApr.h"
using namespace std; //calling the standard directory
TF1 *fitfunc1;
TLine *fitfunc2;
Double_t MinE = 0.03;
Double_t Tstep = 0.2;
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void FindKaons(){
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
    for (int i = 0; i != nhits; i++) cout << endl << Thit[i];
    //for (int i = 0; i != nhits; i++) cout << endl << Ehit[i];
    cout << endl << Tstep << endl;
    for(Int_t ii = 0; ii<nhits; ++ii){
        Double_t minDiff = Thit[ii] - minMean;
       // cout<<"Time Diff "<<minDiff<<endl;
        if((Ehit[ii] > MinE)) { // ignore "pedestal"
            if((minDiff > Tstep) && (minDiff < 70.0)) {
              cout << "if((minDiff > Tstep) && (minDiff < 70.0)) {" << endl;
              //printf("%f,%f", Xhit[ii],Yhit[ii]);
              //printf("%f,%f,%f", minDiff, Ehit[ii],Thit[ii]);
              cout << endl;
               // cout<<"Target dE dt "<<ii<<" "<<Ehit[ii]<<" "<<Thit[ii]<<" "<<Xhit[ii]<<" "<<Yhit[ii]<<endl;
                hXY->Fill(Xhit[ii],Yhit[ii]); //fill tracking (black triangles)
                hTvsE->Fill(Ehit[ii],Thit[ii]);
            }
            if((minDiff > -1.0) && (minDiff <= Tstep)) {
              cout << "if((minDiff > -1.0) && (minDiff <= Tstep)) {" << endl;
              printf("%f,%f", Xhit[ii],Yhit[ii]);
              //printf("%f,%f,%f", minDiff, Ehit[ii],Thit[ii]);
              cout << endl;
              //  cout<<"Target dE low dt "<<ii<<" "<<Ehit[ii]<<" "<<Thit[ii]<<" "<<Xhit[ii]<<" "<<Yhit[ii]<<endl;
                hXYlowT->Fill(Xhit[ii],Yhit[ii]); //fill kaon (red triangles)
                hTvsElowT->Fill(Ehit[ii],Thit[ii]);
            }
        }
    }
    return;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int FindTrack(Int_t Run_Number,Int_t Event_Num) //main()  //Main Program
{
    char Title1[100];
    c1 = new TCanvas("Fit1", "Fit 1", 700,900,1300,600);
    c1->Divide(2,1);
    //c2 = new TCanvas("Fill2", "Fill 2", 700,900,1300,600);
    //c2->Divide(3,1);
    cout<<Run_Number<<" FirstEvent "<<FirstEvent<<" "<<Event_Num<<endl;
    if(FirstEvent == 1){
        FirstEvent = -1;
        MakeFindTrackHistos();
    }
    char Name_finput[200];
    sprintf(Name_finput,"G4Run%d.root", Run_Number);
    cout << "File opened:  " << Name_finput << endl;
    //NewValid:
    hXY->Reset();
    hXYTOF1->Reset();
    hXYlowT->Reset();
    //hTOF1->Reset();
    hTvsE->Reset();
    hTvsElowT->Reset();
        int NextEv;

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

    cout << "*****************************************************" << endl;
    //cout << "File opened:  " << Name_finput << endl;
    cout <<  "Event " << Event_Num << endl;

///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Find Hits routine
    //AGAIN:
   Int_t nentries = (Int_t)fChain->GetEntries();

    cout << "Total Number of Events:  " << nentries <<endl;
    Long64_t ivt=Event_Num;
    Event_Out = Event_Num;
    fChain->GetEntry(ivt);

    Int_t Jevent = n;
    cout<< "New Event: " << ivt << " "<< Jevent <<endl;
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       /* for(Int_t ii = 0; ii < nGaps; ii++) {
            TOF1X[ii] = -1000;
            TOF1Y[ii] = -1000;
        }*/
// Target
///+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
    //printf("Now Target \n");
    // Set counters etc to 0
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
    for(Int_t id = 0; id<5; id++){
        if(tof1N[id] >= 0){
            dTOF1s[nTOF1s] = id;
            pTOF1s[nTOF1s] = tof1N[id];
            nTOF1s++;
            cout<< " tof1N_n "<<id<<" "<<tof1N[id]<<endl;
        }

    }

    for(Int_t it = 0; it <12; it++){
       //cout<< " tof1N "<<it<<" "<<tof1N[0]<<" "<<tof1E[0]<<" "<<tof1dt<<endl;
        for(Int_t id = 0; id<nTOF1s; id++){
            if(it == tof1N[id]){
                Int_t irt = tof1N[id];
                cout<<" TOF1 Location "<<irt<<" "<<TOF_Xloc[irt]<<" "<<TOF_Yloc[irt]<<endl;
                hXYTOF1->Fill(TOF_Xloc[irt],TOF_Yloc[irt]);
            }
        }
    }
    // reset hit stuff
    cout << "if((targetdE[ii] > " << MinE << " ) && (targetdt[ii] > 0.01)) {" << endl;
    cout << "ii tells you at what bar the targetdE is" << endl;
    cout << "Add the x- and y-coordinate of the bar into Xhit and YHit" << endl;
    cout << "Add the energy and time of the bar into Ehit and Thit" << endl;
    for(Int_t ii = 0; ii<nBars; ++ii){
        Xhit[ii] = -1000.0;
        Yhit[ii] = -1000.0;
        Ehit[ii] = -1000.0;
        Thit[ii] = -1000.0;
    }
    //Find hits
        for(Int_t ii = 0; ii<nBars; ++ii){
        if((targetdE[ii] > MinE) && (targetdt[ii] > 0.01)) {
            //targetdE (Energy) above Threshold
            //targetdt (time) more than 0.01 (pedestal?)
            //cout<<"Target dE dt "<<ii<<" "<<targetdE[ii]<<" "<<targetdt[ii]<<" "<<Xloc[ii]<<" "<<Yloc[ii]<<endl;
            Xhit[nhits] = Xloc[ii];
            Yhit[nhits] = Yloc[ii];
            Ehit[nhits] = targetdE[ii];
            Thit[nhits] = targetdt[ii];
            nhits++;
            }
        }
    cout<<"Found nhits "<<nhits<<endl;
    FindKaons();
    c1->cd(1);
    //Draw Energy vs. Time
    sprintf(Title1,"Event %d/%d",Event_Num,Jevent);
    hTvsE->SetXTitle(Title1);
    hTvsE->Draw(); //Draw track (low energy, high time i.e. late event)
    hTvsElowT->Draw("same"); //Draw Kaon (high energy, low time i.e. early event)
    //Drawing X vs. Y (position)
    c1->cd(2);
    hXY->Draw(); //Draws track positions (black triangles)
    hXYlowT->Draw("same"); //Draws Kaon Positions (red triangle)
    hXYTOF1->Draw("same"); //Draws XYTOF1 (pink triangle)

    cout<< "Next Event: " << endl;

     return 0;
}
