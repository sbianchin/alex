#ifndef __CINT__
#include <stdio.h>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>
#include <queue>
#include <algorithm>
#include <TROOT.h>
#include <typeinfo>
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
#include "TMath.h"
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
#include "intersect.cxx"
#include "GlobalParameters.h"
#include "G4DataRootApr19.h" //defines the variables in the .root files
#include "TweakParameters.h"
#include "TOFsParameters.h"
#include "TargetParameters.h"
#include "G4Parameters.h"
#include "G4HistosApr.h"
#include "TPaveText.h"
using namespace std;

#include "Many_Tracks_Fit_2.2.h"
#include "compute_largest_triangle.C"
#include "breadth_first_search.C"

#include "EventPlot.C"
#endif


// global EventPlots to be plottet into by functions and drawn later
TCanvas *c2;
EventPlot* track1_event_plot;
EventPlot* track2_event_plot;
EventPlot* both_tracks_event_plot;

#include "FindKaons.C" // needs the c2 canvas

char canvas_directory[50] = "run816";

pair<vector<int>, vector<int>> stop_at_primary(map <int, vector<int>> lepton_connections_graph, vector<int> vec_l1, vector<int> vec_l2) {
  vector<int> primary_track, secondary_track;
  primary_track = vec_l1;

  int border_primary_leptons = 0;
  // Stop the delta when it hits the primary
  // Currently the delta makes it back to k-stop
  for (unsigned int i = 0; i < vec_l2.size(); ++i) {
      // Does current bar border any of the primary track bars?
      bool border_primary = false;
      for (unsigned int l1_bar = 0; l1_bar < vec_l1.size(); ++l1_bar) {
          vector<int> border_leptons = lepton_connections_graph[vec_l1[l1_bar]];
          if (find(border_leptons.begin(), border_leptons.end(),
            vec_l2[i]) != border_leptons.end()) {
              border_primary = true;
              break;
            }
      }
      if (border_primary) {
        border_primary_leptons += 1;
        //Stops back-tracking when there are 2 bars that borders the primary
        if (border_primary_leptons < 3) {
          secondary_track.push_back(vec_l2[i]);
        }
      }
      else {
        secondary_track.push_back(vec_l2[i]);
      }
  }

  return {primary_track, secondary_track};

}



FitParams GetFitParams(vector<double> vec_x, vector<double> vec_y, vector<double> vec_ex, vector<double> vec_ey) {
    TGraphErrors* track = new TGraphErrors(vec_x.size(),&vec_x[0],
          &vec_y[0], &vec_ex[0],
          &vec_ey[0]);
    TF1 *track_fit = new TF1("track", "pol1");
    FitParams fit_params;
    if (vec_x.size() > 0) {
      track->Fit("track","QW");
      track->Fit("track","QC");
      delete track_fit;
      track_fit = track->GetFunction("track");


      fit_params.ChiS = track_fit->GetChisquare();
      fit_params.ndf = track_fit->GetNDF();
      fit_params.m = track_fit->GetParameter(1);
      fit_params.b = track_fit->GetParameter(0);

    }
    delete track;

    return fit_params;
  }


/*
 * Determine the primary track out of an list of separate tracks. For every
 * track a fit with and without the TOF1 hit point with highest energy is
 * computed. The track with the minimal angle between these fits is promoted to
 * be the primary track and is moved to the first element of the returned list
 * of tracks. The TOF1 hit is weighted 3.
 *
 * This method was chosen over using the track with the minimal chi^2/ndf of a
 * fit with the TOF1, because that method was prone to misidentification if the
 * real primary track scattered a lot, the secondart didn't and the angle
 * difference was small. Other possible methods include looking of tracks that
 * extend up until the outer bound of the detector and align well with the
 * TOF1.
 *
 * double tof1_hit_x      x coordinate of TOF1 hit point
 * double tof1_hit_y      y coordinate of TOF1 hit point
 * int gap                index of the cap the TOF1 was located in
 * vector< vector<int> >  a list of tracks, themselves lists of bar numbers
 *
 * returns: vector< vector<int> >  a list of tracks, themselves lists of bar
 *                                 numbers. The primary track is located at
                                   .first()
 */
vector< vector<int> >  DeterminePrimaryTrack(double tof1_hit_x, double tof1_hit_y, int gap, const int tof1_weight, vector< vector<int> > bar_tracks) {
	double min_angle_diff = std::numeric_limits<double>::infinity();
	vector< vector<int> >::iterator min_angle_diff_track_it = bar_tracks.begin();

	// can't use range-based loop because we need to access the iterator
	for(auto bar_track = bar_tracks.begin(); bar_track != bar_tracks.end(); ++bar_track) {
		if (bar_track->size() > 0) {
			// map bar track to track of spacial coordiantes
			vector<double> space_track_x;
			vector<double> space_track_y;
			vector<double> space_track_ex;
			vector<double> space_track_ey;
			for (vector<int>::const_iterator i = bar_track->begin(); i != bar_track->end(); ++i) {
				int bar = *i;
				space_track_x.push_back(Xloc[bar]);
				space_track_y.push_back(Yloc[bar]);
				space_track_ex.push_back(TARGET_Errors_X);
				space_track_ey.push_back(TARGET_Errors_Y);
			}

			FitParams fit_params_notof = GetFitParams(space_track_x, space_track_y, space_track_ex, space_track_ey);

			// add in highest energy TOF1 hit, weighted by 3
			for (int i = 0; i < tof1_weight; ++i) {
				space_track_x.push_back(tof1_hit_x);
				space_track_y.push_back(tof1_hit_y);
				space_track_ex.push_back(TOF1_Errors_X[gap][0]);
				space_track_ey.push_back(TOF1_Errors_Y[gap][0]);
			}

			FitParams fit_params_tof = GetFitParams(space_track_x, space_track_y, space_track_ex, space_track_ey);

			const double angle_diff = TMath::Abs(TMath::ATan(fit_params_notof.m) - TMath::ATan(fit_params_tof.m));
			if (angle_diff < min_angle_diff) {
				min_angle_diff = angle_diff;
				min_angle_diff_track_it = bar_track;
			}
		}
	}

	// move track with minimum angle difference to the front of the vector
	std::iter_swap(min_angle_diff_track_it, bar_tracks.begin());
	return bar_tracks;
}


/*
 * Maps a bar track (vector<int>) to spacial coordianates and performs a linear
 * fit with respect to uncertainties. A TOF1 hit can be included if the
 * optional parameter tof1_weight is set to a value > 0 (and the consecutive
 * arguments are also given).
 *
 * vector<int> bar_track  list of bar IDs of a track
 * int tof1_weight        optional, number of times the TOF1 hit is added to
 *                          the track before fitting
 * double tof1_x          optional, x coordinate of TOF1 hit point
 * double tof1_y          optional, y coordinate of TOF1 hit point
 * int gap                optional, index of the gap the TOF1 hitpoint is
 *                          located in.
 *
 * returns: FitParams  object containing various fit paramets (chi^2, ndf, m, b)
 */
FitParams GetFitParamsFromBars(vector<int> bar_track, int tof1_weight, double tof1_x, double tof1_y, int gap) {
	vector<double> track_x;
	vector<double> track_y;
	vector<double> track_ex;
	vector<double> track_ey;

	for (vector<int>::const_iterator i = bar_track.begin(); i != bar_track.end(); ++i) {
		int bar = *i;
		track_x.push_back(Xloc[bar]);
		track_y.push_back(Yloc[bar]);
		track_ex.push_back(TARGET_Errors_X);
		track_ey.push_back(TARGET_Errors_Y);
	}

	for (int i = 0; i < tof1_weight; i++) {
		track_x.push_back(tof1_x);
		track_y.push_back(tof1_y);
		track_ex.push_back(TOF1_Errors_X[gap][0]);
		track_ey.push_back(TOF1_Errors_Y[gap][0]);
	}

	return  GetFitParams(track_x, track_y, track_ex, track_ey);
}


/*
 * Calculate an oriented angle for a line fit. The angle returned is in
 * [0, 360) degrees, with 0 at 3 o'clock. In order to find the orientation of
 * the line, the TOF1 hit point is used.
 *
 * Double_t slope      the slope of the line (fit parameter m)
 * Double_t tof1_x     x coordinate of the TOF1 hit point
 * Double_t tof1_y     y coordinate of the TOF1 hit point
 *
 * returns: Double_t   angle in degrees
 */
Double_t LineToDisplayAngle(Double_t slope, Double_t tof1_x, Double_t tof1_y) {
	Double_t angle = TMath::ATan(slope) * TMath::RadToDeg();

	if (tof1_x < 0) {
		angle += 180.;
	} else if (tof1_y < 0) {
		angle += 360.;
	}

	if (angle > 360.) {
		angle -= 360.;
	}

	return angle;
}


/*
 * Calculated an oriented angle for a point in the plane. The angle returned is
 * in [0, 360) degrees, with 0 at 3 o'clock.
 *
 * Coordinate point  the x,y coordinates of the point
 *
 * returns: Double_t  angle in degrees
 */
Double_t PointToDisplayAngle(Coordinate point) {
	Double_t angle = TMath::ATan2(point.y, point.x) * TMath::RadToDeg();
	// convert from [-180, 180) to [0, 360)
	angle = std::fmod(angle + 360, 360);
	return angle;
}


/* Calculate the difference of oriented angles. This has to take care of the
 * wrap-around at 0/360. The optinal argument limits the return value to sharp
 * angles (< 90) and defaults to true.
 *
 * Double_t angle_1
 * Double_t angle_2   the angles to calculate the difference of, in degrees and
 *                      in the range [0, 360)
 * Bool_t only_sharp  optional, defaults to true, wether angles > 90 deg are
 *                      returned or not.
 *
 * returns: Double_t  the difference of the angles in dregrees.
 */
Double_t DifferenceDisplayAngle(Double_t angle_1, Double_t angle_2, Bool_t only_sharp) {
	Double_t diff = 180. - TMath::Abs( TMath::Abs(angle_1 - angle_2) - 180.);

	if (only_sharp && diff > 90.) {
		diff = 180. - diff;
	}

	return diff;
}


/*
 * Find the TOF1 hit point with maximum energy and the corresponding gap.
 *
 * returns: std::pair<Int_t, Coordinate>
 *    .first   the index of the gap the hitpoint is located in. Use with the
 *             array tof1N[]
 *    .second  a Coordinate pair with x and y location of the hit point
 */
std::pair<Int_t, Coordinate> MaxEnergyTOF1Hit(void) {
  Int_t max_TOF_energy_track = -1;
  Double_t max_TOF_energy = -1;

  for (int i = 0; i < 5; ++i) {
    if (max_TOF_energy < tof1E[i]) {
      max_TOF_energy_track = i;
      max_TOF_energy = tof1E[i];
    }
  }

  const Coordinate tof1_hit(tof1wpos[max_TOF_energy_track][0], tof1wpos[max_TOF_energy_track][1]);

  return std::pair<Int_t, Coordinate>(max_TOF_energy_track, tof1_hit);
}


TrackVariables FindTrack(vector<int> lepton_hit_bars, double dv_x, double dv_y, double tof1_hit_x, double tof1_hit_y, int Gap, Bool_t Switch_Display) {

  //Determine the 3 points determining the largest triangle_area

  vector<int> largest_triangle = compute_largest_triangle(lepton_hit_bars);

  //Determine the lepton bar that lies closest to k-stop
  double closest_lepton_dist_to_kstop = 100000; 
  int closest_lepton_to_kstop = 999.;

  // if (Switch_Output == 1) {
  //   cout << "\nBars determining largest triangle:\n";
  //   for (int i = 0; i < 3; ++i) {
  //     cout << largest_triangle[i] << ", " << Xloc[largest_triangle[i]] << ", " << Yloc[largest_triangle[i]] << endl;
  //   }
  //   cout << endl;
  // }

  for (int i = 0; i < 3; ++i) {
    double dist_to_k_stop = distance(Xloc[largest_triangle[i]],
                  Yloc[largest_triangle[i]], dv_x*10, dv_y*10);
    // Determine bar in triangle closest to K-stop
    if (dist_to_k_stop < closest_lepton_dist_to_kstop) {
      closest_lepton_dist_to_kstop = dist_to_k_stop;
      closest_lepton_to_kstop = largest_triangle[i];
    }
  }

  //For an event to be a y-event, make a shortest-path breadth-first-traversal
  //The track to get from one point to another should take a v-shape
  vector<int> lepton_ends;
  for (int i = 0; i < 3; ++i) {
    if (closest_lepton_to_kstop != largest_triangle[i]) {
      lepton_ends.push_back(largest_triangle[i]);
    }
  }

  //Determine bar in lepton_hit_bars closest to k-stop
  closest_lepton_dist_to_kstop = 100000;
  for (unsigned int i = 0; i < lepton_hit_bars.size(); ++i) {
    double dist_to_k_stop = distance(Xloc[lepton_hit_bars[i]],
                  Yloc[lepton_hit_bars[i]], dv_x*10, dv_y*10);
    if (dist_to_k_stop < closest_lepton_dist_to_kstop) {
      closest_lepton_dist_to_kstop = dist_to_k_stop;
      closest_lepton_to_kstop = lepton_hit_bars[i];
    }
  }

  vector<int> vec_l1, vec_l2;
  if (PATH_TRAVERSIAL_USE_ALL) {
    //First track (bar numbers stored in vec_l1) is one endpoint back to k-stop bar
	vec_l1 = FindShortestPath(lepton_hit_bars, lepton_ends[0], closest_lepton_to_kstop);
    //Second track (bar numbers stored in vec_l2) is the other endpoint back to k-stop bar
    vec_l2 = FindShortestPath(lepton_hit_bars, lepton_ends[1], closest_lepton_to_kstop);
  } else {
    path_length traversal = dijkstra(lepton_hit_bars, lepton_ends[0], closest_lepton_to_kstop);
    vec_l1 = traversal.path;
    traversal = dijkstra(lepton_hit_bars, lepton_ends[1], closest_lepton_to_kstop);
	vec_l2 = traversal.path;
  }

  //Reassign the rest of the bars
  map <int, vector<int>> lepton_connections_graph;
  map<int, bool> already_assigned;
  queue<int> track_assign;
  //Build graph of lepton bars and their edges
  for (unsigned int i = 0; i < lepton_hit_bars.size(); ++i) {
    already_assigned[lepton_hit_bars[i]] = false;
    for (int neighbor = 0; neighbor < 8; ++neighbor) {
      int adj_node = TARGET_neighbours[lepton_hit_bars[i]][neighbor];
      if (find(lepton_hit_bars.begin(), lepton_hit_bars.end(),
        adj_node) != lepton_hit_bars.end()) {
          lepton_connections_graph[lepton_hit_bars[i]].push_back(adj_node);
        }
    }

  }

  //If a lepton bar borders one track but doesn't border the other track,
  // then it gets added to that track
  for (unsigned int i = 0; i < vec_l1.size(); ++i) already_assigned[vec_l1[i]] = true;
  for (unsigned int i = 0; i < vec_l2.size(); ++i) already_assigned[vec_l2[i]] = true;
  for (unsigned int i = 0; i < vec_l1.size(); ++i) {

      vector<int> neighbors = lepton_connections_graph[vec_l1[i]];
      for (unsigned int neighbor = 0; neighbor < neighbors.size(); ++neighbor) {
        if (already_assigned[neighbors[neighbor]] == false) {
              track_assign.push(neighbors[neighbor]);
              //Check if the adjacent node only borders vec_l1
              vector<int> adj_node_neighbors = lepton_connections_graph[neighbors[neighbor]];
              bool only_borders_l1 = true;
              for (unsigned int j = 0; j < adj_node_neighbors.size(); ++j) {
                  if (find(vec_l2.begin(), vec_l2.end(),
                        adj_node_neighbors[j]) != vec_l2.end()

                  ) {
                    only_borders_l1 = false;

                  }
              }
              if (only_borders_l1) {
                vec_l1.push_back(neighbors[neighbor]);
                already_assigned[neighbors[neighbor]] = true;
              }
          }

      }
  }
  for (unsigned int i = 0; i < vec_l2.size(); ++i) {
      // cout << "\nCurrent bar is: " << vec_l2[i] << "\n";
      // cout << "The neighbors are: ";

      vector<int> neighbors = lepton_connections_graph[vec_l2[i]];
      // for (int j = 0; j < neighbors.size(); ++j) cout << neighbors[j] << ", ";
      for (unsigned int neighbor = 0; neighbor < neighbors.size(); ++neighbor) {
          //If not already assigned, then push them into track assign queue
          if (already_assigned[neighbors[neighbor]] == false) {
              track_assign.push(neighbors[neighbor]);
              //Check if the adjacent node only borders vec_l1
              vector<int> adj_node_neighbors = lepton_connections_graph[neighbors[neighbor]];
              bool only_borders_l2 = true;
              for (unsigned int j = 0; j < adj_node_neighbors.size(); ++j) {
                  if (find(vec_l1.begin(), vec_l1.end(),
                        adj_node_neighbors[j]) != vec_l1.end()

                  ) {
                    only_borders_l2 = false;

                  }
              }
              if (only_borders_l2) {
                vec_l2.push_back(neighbors[neighbor]);
                already_assigned[neighbors[neighbor]] = true;
              }




          }

      }
  }

  //What are the furthest points in tracks 1 and 2 from one lepton-end?
  double farthest_lepton_from_k_stop = 0.0;
  for (unsigned int i = 0; i < vec_l1.size(); ++i) {
    // cout << vec_l1[i] << ", " << Xloc[vec_l1[i]] << ", " << Yloc[vec_l1[i]] << endl;
    double dist_to_k_stop = distance(Xloc[vec_l1[i]],
                  Yloc[vec_l1[i]], Xloc[closest_lepton_to_kstop],
                                Yloc[closest_lepton_to_kstop]);
    // cout << "Distance to from lepton hit to lepton closest to k-stop: " << dist_to_k_stop << endl;
    if (dist_to_k_stop > farthest_lepton_from_k_stop) {
      farthest_lepton_from_k_stop = dist_to_k_stop;
      lepton_ends[0] = vec_l1[i];
    }
  }
  farthest_lepton_from_k_stop = 0.0;
  for (unsigned int i = 0; i < vec_l2.size(); ++i) {
    // cout << vec_l2[i] << ", " << Xloc[vec_l2[i]] << ", " << Yloc[vec_l2[i]] << endl;
    double dist_to_k_stop = distance(Xloc[vec_l2[i]],
                  Yloc[vec_l2[i]], Xloc[closest_lepton_to_kstop],
                                Yloc[closest_lepton_to_kstop]);
    // cout << "Distance to from lepton hit to lepton closest to k-stop: " << dist_to_k_stop << endl;
    if (dist_to_k_stop < farthest_lepton_from_k_stop) {
      farthest_lepton_from_k_stop = dist_to_k_stop;
      lepton_ends[1] = vec_l2[i];
    }
  }



  // Is there a path from one end point to another?
  path_length traversal = dijkstra(lepton_hit_bars, lepton_ends[0], lepton_ends[1]);
  vector<int> end_to_end = traversal.path;
  // cout << "Path has length: " << traversal.length << endl;
  // cout << "Traversing from one end to another: ";
  // for (int i = 0; i < traversal.path.size(); ++i) {
  //  cout << traversal.path[i] << ", ";
  // }
  // cout << endl;

  //If there's no path, then not a y-event
  bool y_event = true;
  if (end_to_end.size() == 0) {
	  //cout << "Not a Y-event!";
	  y_event = false;
  }
  if (vec_l1.size() == 0 || vec_l2.size() == 0) {
	  //cout << "Not a Y-event!";
	  y_event = false;
  }

  /*
  //Check for v_shape
  if (vec_l1.size() != 0 && vec_l2.size() != 0 && end_to_end.size() != 0) {
    bool v_shape = true;
  
    int array_length = end_to_end.size() - 1;
    int horizontal_shape[array_length];
    int vertical_shape[array_length];
    cout << "\nHorizontal shape: ";
    for (int i = 0; i < array_length; ++i) cout << horizontal_shape[i] << ", ";
    cout << "\nVertical  shape: ";
    for (int i = 0; i < array_length; ++i) cout << vertical_shape[i] << ", ";
  
  
    for (int i = 0; i < array_length; ++i) {
      horizontal_shape[i] = 0; vertical_shape[i] = 0;
      if (Xloc[end_to_end[i+1]] - Xloc[end_to_end[i]] < 0) horizontal_shape[i] = -1;
      if (0 < Xloc[end_to_end[i+1]] - Xloc[end_to_end[i]]) horizontal_shape[i] = 1;
      if (Xloc[end_to_end[i+1]] - Xloc[end_to_end[i]] == 0) {
        if (i > 0 && horizontal_shape[i - 1] != 0) horizontal_shape[i] = horizontal_shape[i - 1];
      }
  
      if (Yloc[end_to_end[i+1]] - Yloc[end_to_end[i]] < 0) vertical_shape[i] = -1;
      if (0 < Yloc[end_to_end[i+1]] - Yloc[end_to_end[i]]) vertical_shape[i] = 1;
      if (Yloc[end_to_end[i+1]] - Yloc[end_to_end[i]] == 0) {
        if (i > 0 && vertical_shape[i - 1] != 0) vertical_shape[i] = vertical_shape[i - 1];
      }
    }
  
    int shape = 0;
    bool all_zeros_horizontal = false;
    while (horizontal_shape[shape] == 0 && shape < array_length - 1) {
      if (shape == array_length - 2 && horizontal_shape[shape + 1] == 0) {
        all_zeros_horizontal = true;
      }
      if (horizontal_shape[shape + 1] != 0) {
        for (;shape >= 0; --shape) horizontal_shape[shape] = horizontal_shape[shape + 1];
        break;
      }
      ++shape;
    }
  
    shape = 0;
    bool all_zeros_vertical = false;
    cout << "\narray_length: " << array_length << "\n";
  
    while (vertical_shape[shape] == 0 && shape < array_length - 1) {
      if (shape == array_length - 2 && vertical_shape[shape + 1] == 0) {
        all_zeros_vertical = true;
      }
      if (vertical_shape[shape + 1] != 0) {
  
        for (;shape >= 0; --shape) vertical_shape[shape] = vertical_shape[shape + 1];
  
        break;
      }
      ++shape;
    }
  
    cout << "\nHorizontal shape: ";
    for (int i = 0; i < array_length; ++i) cout << horizontal_shape[i] << ", ";
    cout << "\nVertical  shape: ";
    for (int i = 0; i < array_length; ++i) cout << vertical_shape[i] << ", ";
  
    //If v-shaped, then is y-event
  
    int horizontal_total = 0; int vertical_total = 0;
    for (int i = 0; i < array_length; ++i) {
      horizontal_total += horizontal_shape[i];
      vertical_total += vertical_shape[i];
    }
    if ((abs(horizontal_total) == array_length || all_zeros_horizontal)
      && (abs(vertical_total) == array_length || all_zeros_vertical)) {
      cout << "Not a y-event";
      v_shape = false;
      y_event = false;
      // How to Handle vertical and horizontal events?
    }
  
    //Determine the density of the path between two ends
    if (!v_shape) {
        vector<int> density, border;
        // Determine the highest and lowest X, Y in path endpoints
        double max_X, min_X, max_Y, min_Y;
  
        if (Xloc[lepton_ends[0]] <= Xloc[lepton_ends[1]]) {
          max_X = Xloc[lepton_ends[1]]; min_X = Xloc[lepton_ends[0]];
        }
        else {
          max_X = Xloc[lepton_ends[0]]; min_X = Xloc[lepton_ends[1]];
        }
  
        if (Yloc[lepton_ends[0]] <= Yloc[lepton_ends[1]]) {
          max_Y = Yloc[lepton_ends[1]]; min_Y = Yloc[lepton_ends[0]];
        }
        else {
          max_Y = Yloc[lepton_ends[0]]; min_Y = Yloc[lepton_ends[1]];
        }
  
  
        for (int i = 0; i < 256; ++i) {
          // If lepton bars lie within the box defined by highest, lowest X, Y
          if (min_X <= Xloc[i] <= max_X && min_Y <= Yloc[i] <= max_Y) {
            density.push_back(i);
          }
          if (min_X <= Xloc[i] <= max_X && min_Y <= Yloc[i] <= max_Y) {
            density.push_back(i);
          }
          // Check lepton bars on the borders
          if (horizontal_total > 0 && vertical_total > 0
              && Xloc[i] == min_X && Yloc[i] == max_Y) {
              border.push_back(i);
            }
          if (horizontal_total > 0 && vertical_total < 0
              && Xloc[i] == min_X && Yloc[i] == max_Y) {
            border.push_back(i);
            }
          if (horizontal_total < 0 && vertical_total > 0
              && Xloc[i] == min_X && Yloc[i] == max_Y) {
            border.push_back(i);
            }
        }
    }
  }
  vector_lepton_kaon.y_event = y_event;
  */

  const int tof1_weight = FIT_TOF1_WEIGHT;

  // determine primary track
  vector< vector<int> > bar_tracks {vec_l1, vec_l2};
  vector< vector<int> > sorted_tracks = DeterminePrimaryTrack(tof1_hit_x, tof1_hit_y, Gap, tof1_weight, bar_tracks);
  vector<int> primary_track_bars = sorted_tracks[0];
  vector<int> secondary_track_bars = sorted_tracks[1];

  // only keep bars directly aligned with the fits on the tracks
  pair<vector<int>, vector<int>> truncated_track_bars = stop_at_primary(lepton_connections_graph, primary_track_bars, secondary_track_bars);
  vector<int> primary_track = truncated_track_bars.first;
  vector<int> secondary_track = truncated_track_bars.second;
  if (primary_track.size() == 0) {
    primary_track = secondary_track;
    secondary_track = {};
  }

  // calculate angles, energy and plot fits
  // TOF1 hit is included with tof1_weight only in primary track
  FitParams fit_params_primary = GetFitParamsFromBars(primary_track, tof1_weight, tof1_hit_x, tof1_hit_y, Gap);
  FitParams fit_params_secondary = GetFitParamsFromBars(secondary_track);

  Double_t angle_primary = LineToDisplayAngle(fit_params_primary.m, tof1_hit_x, tof1_hit_y);
  Double_t angle_secondary = LineToDisplayAngle(fit_params_secondary.m, tof1_hit_x, tof1_hit_y);
  Double_t angle_between = DifferenceDisplayAngle(angle_primary,angle_secondary);

  TF1* fit_primary = new TF1("fit_primary", "[0]*x+[1]", -50, 50);
  fit_primary->SetParameter(0, fit_params_primary.m);
  fit_primary->SetParameter(1, fit_params_primary.b);
  if(Switch_Display == 1) {
	  track1_event_plot->ShowTF1(fit_primary, kRed);
	  both_tracks_event_plot->ShowTF1(fit_primary, kRed);
  }

  TF1* fit_secondary = new TF1("fit_secondary", "[0]*x+[1]", -50, 50);
  fit_secondary->SetParameter(0, fit_params_secondary.m);
  fit_secondary->SetParameter(1, fit_params_secondary.b);
  if(Switch_Display == 1) {
	  track2_event_plot->ShowTF1(fit_secondary, kGreen);
	  both_tracks_event_plot->ShowTF1(fit_secondary, kGreen);
  }

  double primary_bars_energy = 0;
  for (unsigned int i = 0; i < primary_track.size(); ++i) {
    primary_bars_energy += targetdE[primary_track[i]];
  }
  double E_positron = tof1E[0] + delta_energy + primary_bars_energy;

  TrackVariables track_variables = {
    .primary_track       = primary_track,
    .first_track_ChiS    = fit_params_primary.ChiS,
    .ndf_1               = fit_params_primary.ndf,
    .angle_primary       = angle_primary,

    .secondary_track     = secondary_track,
    .second_track_ChiS   = fit_params_secondary.ChiS,
    .ndf_2               = fit_params_secondary.ndf,
    .angle_secondary     = angle_secondary,

    .angle_between       = angle_between,
    .E_positron          = E_positron,
    .primary_bars_energy = primary_bars_energy,
    .vec_l1              = vec_l1,
    .vec_l2              = vec_l2
  };

  return track_variables;
}


Lepton FindTwoTracks(Int_t Run_Number, Int_t ievt, int savefile, Int_t Switch_Display, Int_t Switch_Output, Int_t batch){
  //Open Batch Output file
  char batch_output[50];
  sprintf(batch_output, "RUN_%d_Find_Two_Tracks.txt", Run_Number);
  ofstream fout;
  fout.open(batch_output, ios::app);

  if (Switch_Display == 1) {
	  c2 = new TCanvas("Event_Display.C  --  TARGET & SFT","Event_Display.C  --  TARGET & SFT",500,500,1050,700);
	  c2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
	  c2->Divide(3,2);

	  track1_event_plot = new EventPlot("Track 1");
	  track2_event_plot = new EventPlot("Track 2");
	  both_tracks_event_plot = new EventPlot("Both Tracks");
  }

  //Perform a Single Track Fit
  Lepton vector_lepton_kaon = MC_Event_Display(Run_Number, ievt, Switch_Display, Switch_Output, batch);


  //Unpack Lepton and Kaon Tracks
  vector<double> vec_xx_lepton = vector_lepton_kaon.vec_xx_lepton;
  vector<double> vec_yy_lepton = vector_lepton_kaon.vec_yy_lepton;
  vector<double> vec_xx_lepton_rotate = vector_lepton_kaon.vec_xx_lepton_rotate;
  vector<double> vec_yy_lepton_rotate = vector_lepton_kaon.vec_yy_lepton_rotate;
  vector<double> vec_xx_kaon = vector_lepton_kaon.vec_xx_kaon;
  vector<double> vec_yy_kaon = vector_lepton_kaon.vec_yy_kaon;
  vector<int> lepton_hit_bars = vector_lepton_kaon.lepton_hit_bars;
  vector<int> kaon_hit_bars = vector_lepton_kaon.kaon_hit_bars;
  double k_stop_x = vector_lepton_kaon.k_stop_x;
  double k_stop_y = vector_lepton_kaon.k_stop_y;
  Int_t tof1N[5];
  copy(vector_lepton_kaon.tof1N, vector_lepton_kaon.tof1N+5, tof1N);
  double reduced_loss = vector_lepton_kaon.reduced_loss;
  double reduced_ChiS = vector_lepton_kaon.reduced_ChiS;
  int ndf = vector_lepton_kaon.ndf;

  const std::pair<Int_t, Coordinate> tof1_hit = MaxEnergyTOF1Hit();
  const Int_t max_TOF_energy_track = tof1_hit.first;
  const Double_t tof1_hit_x = tof1_hit.second.x;
  const Double_t tof1_hit_y = tof1_hit.second.y;

  const Double_t k_stop_error = distance(k_stop_x, k_stop_y, dv_x*10, dv_y*10);
  const Double_t k_stop_radius = TMath::Sqrt( TMath::Sq(k_stop_x) + TMath::Sq(k_stop_y) );

  // use exit angle if provided by simulation data. Otherwise use angle of
  // line through k-stop and TOF1
	Double_t angle_primary_sim = std::nan("");
  if (fChain.GetListOfBranches()->FindObject("uto_x")) {
	  angle_primary_sim = PointToDisplayAngle(Coordinate(uto_x, uto_y));
  } else {
	  const Double_t slope_primary_sim = (dv_y*10 - tof1_hit_y) / (dv_x*10 - tof1_hit_x);
	 angle_primary_sim = LineToDisplayAngle(slope_primary_sim, tof1_hit_x, tof1_hit_y);
  }

  // Only Do 2 track reconstruction if one-track fit is bad and
  if (reduced_loss > TWO_TRACK_MIN_CHISQ && lepton_hit_bars.size() > TWO_TRACK_MIN_LEPTON_BARS) {
     TrackVariables two_tracks = FindTrack(lepton_hit_bars, dv_x, dv_y, tof1_hit_x, tof1_hit_y, tof1N[max_TOF_energy_track], Switch_Display);

     double first_track_ChiS = two_tracks.first_track_ChiS;
     double ndf_1 = two_tracks.ndf_1;
     double second_track_ChiS = two_tracks.second_track_ChiS;
     double ndf_2 = two_tracks.ndf_2;
     double E_positron = two_tracks.E_positron;
     double primary_bars_energy = two_tracks.primary_bars_energy;
     double angle_between = two_tracks.angle_between;
	   double angle_primary = two_tracks.angle_primary;
	   double angle_secondary = two_tracks.angle_secondary;
     vector<int> vec_l1 = two_tracks.vec_l1;
     vector<int> vec_l2 = two_tracks.vec_l2;
     vector<int> primary_track = two_tracks.primary_track;
     vector<int> secondary_track = two_tracks.secondary_track;

	 const Double_t delta_angle_primary = DifferenceDisplayAngle(angle_primary, vector_lepton_kaon.angle_lepton_all);
	 const Double_t delta_angle_secondary = DifferenceDisplayAngle(angle_secondary, vector_lepton_kaon.angle_lepton_all);
	 const Double_t angle_primary_error = DifferenceDisplayAngle(angle_primary, angle_primary_sim);

      if (Switch_Output == 1) {
        cout << "\nFirst track ChiS: " << first_track_ChiS << ", ndf: " << ndf_1 << ", ChiS/ndf: "<< first_track_ChiS/ndf_1 << "\n";
        cout << "\nSecond track ChiS: " << second_track_ChiS << ", ndf: " << ndf_2 << ", ChiS/ndf: "<< second_track_ChiS/ndf_2 << "\n";

        cout << "\nE_TOF + E_DELTA + sum(dE_leptonbars) = E_positron\n";
        cout << tof1E[0] << " + " << delta_energy << " + " << primary_bars_energy << " = "
          << E_positron << "\n";
      }

      if (Switch_Display == 1) {
		track1_event_plot->ShowBars(primary_track, kRed);
		track1_event_plot->ShowPoint(tof1_hit_x, tof1_hit_y, kCyan, kFullSquare, "TOF1 Hit")->SetMarkerSize(1);
		track1_event_plot->ShowValue("#theta_1 (deg)", angle_primary);
		track1_event_plot->ShowValue("#Delta#theta_1 (deg)", delta_angle_primary);
		track1_event_plot->Draw(c2, 4);

		track2_event_plot->ShowBars(secondary_track, kGreen);
		track2_event_plot->ShowPoint(tof1_hit_x, tof1_hit_y, kCyan, kFullSquare, "TOF1 Hit")->SetMarkerSize(1);
		track2_event_plot->ShowValue("#theta_2 (deg)", angle_secondary);
		track2_event_plot->ShowValue("#Delta#theta_2 (deg)", delta_angle_secondary);
		track2_event_plot->Draw(c2, 5);

		both_tracks_event_plot->ShowBars(primary_track, kRed);
		both_tracks_event_plot->ShowBars(secondary_track, kGreen);
		both_tracks_event_plot->ShowPoint(tof1_hit_x, tof1_hit_y, kCyan, kFullSquare, "TOF1 Hit")->SetMarkerSize(1);
		both_tracks_event_plot->ShowPoint(k_stop_x, k_stop_y, kMagenta, kStar, "K-Stop (determined)");
		both_tracks_event_plot->ShowPoint(dv_x*10, dv_y*10, kBlue, kPlus, "K-Stop (simulation)");
		both_tracks_event_plot->ShowValue("Delta Energy", delta_energy);
		both_tracks_event_plot->ShowValue("inner angle (deg)", angle_between);
		both_tracks_event_plot->ShowValue("#theta_sim", angle_primary_sim);
		both_tracks_event_plot->ShowValue("Primary Error", angle_primary_error);
		both_tracks_event_plot->Draw(c2, 6);

		/*
		// testing code using the old plotting API
		// need to #include "Plot_Simulation_Event_Display.h"
		// needs cTEST global canvas + canvas creation in FindKaons.C
        TCanvas *cTEST = new TCanvas("Testing","Testing",0,200,1050,350);
        cTEST->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
        cTEST->Divide(3,1);
        const char *trackA_test = "Track A";
        const char *trackB_test = "Track B";
        const char *trackA_andB_test = "Track A+B";

        Subplot_4_test(cTEST, vec_l1, 1, trackA_test, 2);
        Subplot_4_test(cTEST, vec_l2, 2, trackB_test, 8);
        Subplot_6_test(cTEST, vec_l1, vec_l2, trackA_andB_test, 3, 2, 8,
          dv_x, dv_y, 0, 0);
        if (savefile == 1) {
          char before_reassignment[100];
          sprintf(before_reassignment, "./%s/RUN_%d_event_%d_before_reassignment.png", canvas_directory, Run_Number, ievt);
          cTEST->SaveAs(before_reassignment);
          // delete cTEST;
        }
		*/

      }

      if (Switch_Display == 1 && savefile == 1) {

        if (batch == 1) {
          char plot_output[50];
          sprintf(plot_output, "RUN_%d_FindTwoTracks.pdf", Run_Number);
          c2->Update();
          c2->Print(plot_output);
          // delete c2;
        }
        else {
          char image_name[100];
          sprintf(image_name, "./%s/RUN_%d_event_%d_many_tracks.png",
          canvas_directory, Run_Number, ievt);
          c2->SaveAs(image_name);
        }

      }

      if (savefile == 1) {
        fout << fixed;
        fout        /* << setw(4) */                    << Run_Number;
        fout << "," /* << setw(5) */                    << ievt;
        fout << "," /* << setw(7) */ << setprecision(2) << reduced_ChiS;
        fout << "," /* << setw(4) */                    << ndf;
        fout << "," /* << setw(7) */ << setprecision(2) << first_track_ChiS/ndf_1;
        fout << "," /* << setw(4) */ << setprecision(0) << ndf_1;
        fout << "," /* << setw(7) */ << setprecision(2) << second_track_ChiS/ndf_2;
        fout << "," /* << setw(4) */ << setprecision(0) << ndf_2;
        fout << "," /* << setw(7) */ << setprecision(2) << E_positron;
        fout << "," /* << setw(7) */ << setprecision(2) << delta_energy;
        fout << "," /* << setw(7) */ << setprecision(2) << delta_length*10; //Raw delta length given in cm
        fout << "," /* << setw(7) */ << setprecision(2) << delta_length_xy; //Calculated xy length in mm
        fout << "," /* << setw(7) */ << setprecision(2) << angle_between;
		fout << "," /*            */                    << "1"; // two track fit
		fout << "," /* << setw(7) */ << setprecision(2) << angle_primary_sim;
		fout << "," /* << setw(7) */ << setprecision(2) << vector_lepton_kaon.angle_lepton_all; // single track fit angle
		fout << "," /* << setw(7) */ << setprecision(2) << angle_primary;
		fout << "," /* << setw(7) */ << setprecision(2) << delta_angle_primary;
		fout << "," /* << setw(7) */ << setprecision(2) << angle_secondary;
		fout << "," /* << setw(7) */ << setprecision(2) << delta_angle_secondary;
		fout << "," /* << setw(7) */ << setprecision(2) << k_stop_error;
		fout << "," /* << setw(7) */ << setprecision(2) << angle_primary_error;
		fout << "," /* << setw(7) */ << setprecision(2) << targL; // track length in target (sim)
		fout << ","                                     << vector_lepton_kaon.no_leptons;
		fout << ","                                     << vector_lepton_kaon.no_kaons;
		fout << "," /* << setw(7) */ << setprecision(2) << k_stop_radius;
		fout << "," /* << setw(7) */ << setprecision(2) << TWO_TRACK_MIN_CHISQ;
		fout << "," /* << setw(7) */ << setprecision(2) << TWO_TRACK_MIN_LEPTON_BARS;
		fout << "," /* << setw(7) */ << setprecision(2) << FIT_TOF1_WEIGHT;
		fout << "," /* << setw(7) */ << setprecision(2) << K_STOP_CENTROID_THRESH;
		fout << "," /* << setw(7) */ << setprecision(2) << PATH_TRAVERSIAL_USE_ALL;
		fout << "," /* << setw(7) */ << setprecision(2) << PATH_TRAVERSIAL_DIJKSTRA_JUMP_RADIUS;
		fout << "," /* << setw(7) */ << setprecision(2) << PATH_TRAVERSIAL_ALL_PENALTY;
		fout << "\n";
      }

  } else {
  if (savefile == 1) {
    double primary_bars_energy = 0;
    for (unsigned int i = 0; i < lepton_hit_bars.size(); ++i) {
      primary_bars_energy += targetdE[lepton_hit_bars[i]];
    }
    double E_positron = tof1E[0] + delta_energy + primary_bars_energy;

	 const Double_t angle_primary_error = DifferenceDisplayAngle(vector_lepton_kaon.angle_lepton_all, angle_primary_sim);

    fout << fixed;
    fout        /* << setw(4) */                    << Run_Number;
    fout << "," /* << setw(5) */                    << ievt;
    fout << "," /* << setw(7) */ << setprecision(2) << reduced_ChiS;
    fout << "," /* << setw(4) */                    << ndf;
    fout << "," /* << setw(9) */ ;                  // first_track_ChiS / ndf_1
    fout << "," /* << setw(6) */ ;                  // ndf_1
    fout << "," /* << setw(9) */ ;                  // second_track_ChisS / ndf_2
    fout << "," /* << setw(6) */ ;                  // ndf_2
    fout << "," /* << setw(7) */ << setprecision(2) << E_positron;
    fout << "," /* << setw(7) */ << setprecision(2) << delta_energy;
    fout << "," /* << setw(7) */ << setprecision(2) << delta_length*10; //Raw delta length given in cm
    fout << "," /* << setw(7) */ << setprecision(2) << delta_length_xy;
    fout << "," /* << setw(9) */ ;                  // angle_between
    fout << "," /*            */                    << "0"; // one track fit
	fout << "," /* << setw(7) */ << setprecision(2) << angle_primary_sim;
	fout << "," /* << setw(7) */ << setprecision(2) << vector_lepton_kaon.angle_lepton_all; // single track fit angle
    fout << "," /* << setw(9) */ ;                  // angle_primary
    fout << "," /* << setw(9) */ ;                  // delta_angle_primary
    fout << "," /* << setw(9) */ ;                  // angle_secondary
    fout << "," /* << setw(9) */ ;                  // delta_angle_secondary
	fout << "," /* << setw(7) */ << setprecision(2) << k_stop_error;
	fout << "," /* << setw(7) */ << setprecision(2) << angle_primary_error;
	fout << "," /* << setw(7) */ << setprecision(2) << targL; // track length in target (sim)
	fout << ","                                     << vector_lepton_kaon.no_leptons;
	fout << ","                                     << vector_lepton_kaon.no_kaons;
	fout << "," /* << setw(7) */ << setprecision(2) << k_stop_radius;
	fout << "," /* << setw(7) */ << setprecision(2) << TWO_TRACK_MIN_CHISQ;
	fout << "," /* << setw(7) */ << setprecision(2) << TWO_TRACK_MIN_LEPTON_BARS;
	fout << "," /* << setw(7) */ << setprecision(2) << FIT_TOF1_WEIGHT;
	fout << "," /* << setw(7) */ << setprecision(2) << K_STOP_CENTROID_THRESH;
	fout << "," /* << setw(7) */ << setprecision(2) << PATH_TRAVERSIAL_USE_ALL;
	fout << "," /* << setw(7) */ << setprecision(2) << PATH_TRAVERSIAL_DIJKSTRA_JUMP_RADIUS;
	fout << "," /* << setw(7) */ << setprecision(2) << PATH_TRAVERSIAL_ALL_PENALTY;
	fout << "\n";
    }
  }

  fout.close();

  return vector_lepton_kaon;
}
