#ifndef MANY_TRACKS_2_2_H
#define MANY_TRACKS_2_2_H

#include "G4DataRootApr19.h"


struct TrackVariables {
	vector<int> primary_track;
	double first_track_ChiS;
	double ndf_1;
	double angle_primary;

	vector<int> secondary_track;
	double second_track_ChiS;
	double ndf_2;
	double angle_secondary;

	double angle_between;
	double E_positron;
	double primary_bars_energy;
	vector<int> vec_l1;
	vector<int> vec_l2;
};

class Coordinate {
	public:
		Coordinate(Double_t x_, Double_t y_, Double_t ex_=0, Double_t ey_=0) : x(x_), y(y_), ex(ex_), ey(ey_) {}


		Double_t x;
		Double_t y;
		Double_t ex;
		Double_t ey;
};

struct FitParams {
	Double_t m;
	Double_t b;
	Double_t ChiS;
	Double_t ndf;
};


pair<vector<int>, vector<int>> stop_at_primary(map <int, vector<int>> lepton_connections_graph, vector<int> vec_l1, vector<int> vec_l2);

vector< vector<int> >  DeterminePrimaryTrack(double tof1_hit_x, double tof1_hit_y, int gap, vector< vector<int> > bar_tracks);

FitParams GetFitParamsFromBars(vector<int> bar_track, int tof1_weight = 0, double tof1_x = std::nan(""), double tof1_y = std::nan(""), int gap = -1);

Double_t LineToDisplayAngle(Double_t slope, Double_t tof1_x, Double_t tof1_y);

Double_t PointToDisplayAngle(Coordinate point);

Double_t DifferenceDisplayAngle(Double_t angle_1, Double_t angle_2, Bool_t only_sharp = true);

std::pair<Int_t, Coordinate> MaxEnergyTOF1Hit(void);

TrackVariables FindTrack(vector<int> lepton_hit_bars, double dv_x, double dv_y, double tof1_hit_x, double tof1_hit_y, int Gap);

Lepton FindTwoTracks(Int_t Run_Number=5, Int_t ievt=0, int savefile = 0, Int_t Switch_Display = 1, Int_t Switch_Output = 1, Int_t batch = 0);

#endif

