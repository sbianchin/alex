#ifndef FIND_KAONS_H
#define FIND_KAONS_H

#include "G4DataRootApr19.h"
#include "Many_Tracks_Fit_2.2.h" // Coordinate type

Int_t StartChain(TChain &fChain, char Name_finput[200]);

double first_fit_residual(vector <double> vec_xx, vector <double> vec_yy, double m, double b);

Coordinate KStopByIntersection(std::vector<Int_t> kaon_hit_bars, std::vector<Double_t> vec_xx_kaon, std::vector<Double_t> vec_yy_kaon, vector<Double_t> vec_xx_lepton, vector<Double_t> vec_yy_lepton);
Coordinate KStopByCentroid(std::vector<Int_t> kaon_hit_bars, map<int, vector<double>> kaon_time_energy_map);

Lepton FindKaons(Int_t Switch_Output);

Lepton MC_Event_Display(Int_t Run_Number=5, Int_t ievt=0, Int_t Switch_Display = 1, Int_t Switch_Output = 1, Int_t batch = 0);

#endif

