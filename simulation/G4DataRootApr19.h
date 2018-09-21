#ifndef G4_DATA_ROOT_APR_19_H
#define G4_DATA_ROOT_APR_19_H

#include "TChain.h"

// Declaration of leaf types
// change lines 98, 100  Min_E_k and two_tracks_min_Chi_sq   MH  Oct/Nov/2017
Int_t           n;
Int_t           s_hit;
Int_t           x_hit;
Double_t        tof1charge[5];
Int_t           tof1N[5];
Double_t        px[10];
Double_t        py[10];
Double_t        pz[10];
Double_t        u_x;
Double_t        u_y;
Double_t        u_z;
Double_t        uto_x; // ⎫ 
Double_t        uto_y; // ⎬ momentum unit vector of lepton at exit of detector
Double_t        uto_z; // ⎭

Double_t        delta_x;
Double_t        delta_y;
Double_t        delta_z;
Double_t        delta_stop_x;
Double_t        delta_stop_y;
Double_t        delta_stop_z;
Double_t        delta_ux;
Double_t        delta_uy;
Double_t        delta_uz;
Double_t        delta_energy;
Double_t        delta_length;
Double_t        delta_length_xy;
Double_t        delta_time;



Double_t        genEvt;
Double_t        genEne[10];
Double_t        cosTh;
Double_t        sinTh;
Double_t        Phi;
Double_t        dv_x;
Double_t        dv_y;
Double_t        dv_z;
Double_t        targetE[5];
Double_t        targetdE[256];
Double_t        targetdt[256];
Double_t        targL;
Double_t        t1pL;
Double_t        tof1E[5];
Double_t        tof1wpos[5][3];
Double_t        tof1dt;
Char_t          secondP[5][100];
Double_t        t1;
Double_t        t1px;
Double_t        t1py;
Double_t        t1pz;
Double_t        SFTwpos1[5][3];
Double_t        SFTwpos2[5][3];
Double_t        SFTpx;
Double_t        SFTpy;
Double_t        SFTpz;

TChain fChain("Kaon");

//Defining Lepton data structure
struct Lepton {
vector<double> vec_xx_lepton;
vector<double> vec_yy_lepton;
vector <double> vec_ex_lepton;
vector <double> vec_ey_lepton;
vector <double> vec_xx_lepton_rotate;
vector <double> vec_yy_lepton_rotate;
vector <double> vec_xx_kaon;
vector <double> vec_yy_kaon;
vector <double> vec_ex_kaon;
vector <double> vec_ey_kaon;
vector <int> lepton_hit_bars;
vector <int> kaon_hit_bars;
vector <double> vec_xx_lepton_low_energy;
vector <double> vec_yy_lepton_low_energy;
vector <double> vec_xx_kaon_low_energy;
vector <double> vec_yy_kaon_low_energy;
vector <double> vec_ex_lepton_low_energy;
vector <double> vec_ey_lepton_low_energy;
vector <double> vec_ex_kaon_low_energy;
vector <double> vec_ey_kaon_low_energy;
Double_t dv_x;
Double_t dv_y;
Int_t tof1N[5];
double reduced_loss = 0;
double reduced_ChiS = 0;
int ndf = 0;
double two_track_reduced_loss = 0;
double first_track_ChiS;
double second_track_ChiS;
bool y_event, high_chi_square_high_bar;
double k_stop_x,k_stop_y;
double angle_lepton_all;
double angle_lepton_1;
double angle_lepton_2;
bool no_leptons = false;
bool no_kaons = false;
map<int, vector<double>> kaon_time_energy_map;
};

//Define cut parameters
Double_t MinE = 0.0;
Double_t Tstep = 0.2;
double min_E_valid_hit = 0.0;
// modified Nov 8 MH  2.0 --> 1.0 MeV
double Min_E_k = 1; double Min_T_k = 3;
// modified Oct 23 MH  2.5 --> 1.0
double two_tracks_min_Chi_sq = 1.0;

#endif
