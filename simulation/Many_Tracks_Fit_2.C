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
typedef dlib::matrix<double,4,1> two_track_parameter_vector;
#include "Plot_Simulation_Event_Display.h"

#include "Experimental_Event_Display_Simulation_1.h"
#endif
using namespace std;







//Defining some data structures
struct Track {
  std::vector<int> lepton_track;
  int kaon_neighbor = 0;
};

struct Coordinate {
  double x;
  double y;
};




double Distance_from_Bar_to_Line(int bar, Coordinate endpoint_1, Coordinate endpoint_2) {
  double v_first = (endpoint_2.x-endpoint_1.x)*(endpoint_1.y-Yloc[bar]);
  double v_second = (endpoint_1.x-Xloc[bar])*(endpoint_2.y-endpoint_1.y);
  double v = abs(v_first - v_second);
  double line_length_squared = pow(endpoint_2.x-endpoint_1.x, 2)
    + pow(endpoint_2.y-endpoint_1.y, 2);
  double line_length = pow(line_length_squared, 0.5);
  return v/line_length;
}


double calculate_double_residual(const vector<double> data,
  const two_track_parameter_vector& params)
  {
  double m1 = params(0); double b1 = params(1); double m2 = params(2); double b2 = params(3);
  double x = data[0]; double y = data[1];
  // cout << "X, Y of Point" << vec_xx_lepton[point] << ", " << vec_yy_lepton[point] << ", ";
  // cout << "Residual 1, 2: " << m1*vec_xx_lepton[point] + b1 - vec_yy_lepton[point] << ", ";
  // cout << m2*vec_xx_lepton[point] + b2 - vec_yy_lepton[point] << ", ";
  double ssr_1 = abs(m1*x + b1 - y); //first SSR
  double ssr_2 = abs(m2*x + b2 - y); //first SSR
  // cout << "ssr1, ssr2: ";
  // cout << ssr_1 << ", " << ssr_2 << endl;
  if (ssr_1 <= ssr_2) {
    return pow(ssr_1, 2);
  }
  else {
    return pow(ssr_2, 2);
  }
}

double calculate_2_track_residual(const vector<double> data,
  const two_track_parameter_vector& params)
  {
  double m1 = params(0); double b1 = params(1); double m2 = params(2); double b2 = params(3);
  double x = data[0]; double y = data[1];

  // if (x > 50) {
  //   x = x - 100; y = y - 100;
  //   double ssr_1 = abs(m1*x + b1 - y);
  //   double ssr_2 = abs(m2*x + b2 - y);
  //   if (ssr_1 <= ssr_2) {
  //     return ssr_1 * 0.25;
  //   }
  //   else return ssr_2 * 0.25;
  //
  // }
  // cout << "X, Y of Point" << vec_xx_lepton[point] << ", " << vec_yy_lepton[point] << ", ";
  // cout << "Residual 1, 2: " << m1*vec_xx_lepton[point] + b1 - vec_yy_lepton[point] << ", ";
  // cout << m2*vec_xx_lepton[point] + b2 - vec_yy_lepton[point] << ", ";
  double ssr_1 = abs(m1*x + b1 - y); //first SSR
  double ssr_2 = abs(m2*x + b2 - y); //first SSR
  // cout << endl << "ssr 1, 2: " << ssr_1 << ", " << ssr_2 << endl;
  // double ssr_1 = distance_to_line(x, y, m1, b1); //first SSR
  // double ssr_2 = distance_to_line(x, y, m2, b2); //first SSR

  // cout << "ssr1, ssr2: ";
  // cout << ssr_1 << ", " << ssr_2 << endl;
  vector<double> intersect = _2lines_intersect(m1, b1, m2, b2);
  double dist = distance(intersect[0], intersect[1], x, y);
  if (ssr_1 <= ssr_2) {
    if (!isnan(dist)) {
        return pow(ssr_1, 1) + dist*0.15;
    }
    else return pow(ssr_1, 1);
  }
  else {
    if (!isnan(dist)) {
        return pow(ssr_2, 1) + dist*0.15;
    }
    else return pow(ssr_2, 1);
  }


}

double fit_to_bar_numbers(vector <int> bars_hit, double m, double b) {
    double loss = 0;
    for (int i = 0; i < bars_hit.size(); ++i) {

      loss += pow((Yloc[bars_hit[i]] - m*Xloc[bars_hit[i]] - b),2);
    }
    return loss;
}



Lepton FindTrack(Int_t Run_Number=5, Int_t ievt=0, int savefile = 0, int current_tof1 = 5,
  Int_t Switch_Display = 1, Int_t Switch_Output = 1){
  //Find Lepton and Kaon Tracks
  Lepton vector_lepton_kaon = Event_Display_5_1(Run_Number, ievt, Switch_Display, Switch_Output);

  //Unpack Lepton and Kaon Tracks
  vector<double> vec_xx_lepton = vector_lepton_kaon.vec_xx_lepton;
  vector<double> vec_yy_lepton = vector_lepton_kaon.vec_yy_lepton;
  vector<double> vec_xx_lepton_rotate = vector_lepton_kaon.vec_xx_lepton_rotate;
  vector<double> vec_yy_lepton_rotate = vector_lepton_kaon.vec_yy_lepton_rotate;
  vector<double> vec_xx_kaon = vector_lepton_kaon.vec_xx_kaon;
  vector<double> vec_yy_kaon = vector_lepton_kaon.vec_yy_kaon;
  vector<int> lepton_hit_bars = vector_lepton_kaon.lepton_hit_bars;
  vector<int> kaon_hit_bars = vector_lepton_kaon.kaon_hit_bars;

  Int_t tof1N[5];
  copy(vector_lepton_kaon.tof1N, vector_lepton_kaon.tof1N+5, tof1N);
  double reduced_loss = vector_lepton_kaon.reduced_loss;
  int ndf = vector_lepton_kaon.ndf;
  cout << "Kaon npts: " << vec_xx_kaon.size() << endl;
  for (vector<int>::const_iterator i = kaon_hit_bars.begin();
    i != kaon_hit_bars.end(); ++i)
    cout << *i << ", ";
  cout << endl << "Lepton npts: " << vec_xx_lepton.size() << endl;
  for (vector<int>::const_iterator i = lepton_hit_bars.begin();
    i != lepton_hit_bars.end(); ++i)
    cout << *i << ", ";
  // cout << "Chi2: " << ChiS << endl;
  // cout << "ndf: " << ndf << endl;
  // cout << "Reduced Chi2: " << ChiS/ndf << endl;

  // Only Do 2 track reconstruction of one-track fit is bad
  if (reduced_loss>0.0) {
    for (vector<int>::const_iterator i = lepton_hit_bars.begin();
      i != lepton_hit_bars.end(); ++i)
      cout << *i << ' ';


    vector<double> k_stop;
    //Calculate K-stop
    //K-stop is where the kaon bars touch the lepton bars
    map<int, vector<int>> intersection_points;
    for (int i = 0; i < kaon_hit_bars.size(); ++i) {
      cout << "Neighbor: ";
      for (int neighbor = 0; neighbor < 8; ++neighbor) {
        cout << TARGET_neighbours[kaon_hit_bars[i]][neighbor] << " ";
        if (find(lepton_hit_bars.begin(), lepton_hit_bars.end(),
          TARGET_neighbours[kaon_hit_bars[i]][neighbor]) != lepton_hit_bars.end()) {
            intersection_points[TARGET_neighbours[kaon_hit_bars[i]][neighbor]].push_back(kaon_hit_bars[i]);
            cout << endl << "kaon intersection: " << kaon_hit_bars[i] << " ";
            cout << "lepton at intersection: " << TARGET_neighbours[kaon_hit_bars[i]][neighbor] << "\n";

          }
        cout << "\n";
      }
    }

    for (map<int, vector<int>>::iterator junction_lepton = intersection_points.begin();
      junction_lepton != intersection_points.end(); ++junction_lepton) {
        // For each junction lepton, find the neighboring kaons

        //candidate_stop.push_back(junction_lepton);
        cout << "Lepton at intersection: " << junction_lepton->first << endl;
        vector<int> candidate_stop = junction_lepton->second;
        cout << "Kaon intersecting: ";
        for (int i = 0; i < candidate_stop.size(); ++i) cout << candidate_stop[i] << " ";
      }
      cout << endl;



    // Plot of Lepton Graph
    TGraph *gr_Leptons = new TGraph(vec_xx_lepton.size(),&vec_xx_lepton[0],&vec_yy_lepton[0]);

    //Check which bars neighboring hit TOF are hit
    vector <int> lepton_ends;
    vector<int> tof1_value (tof1N,tof1N+5); //Remove duplicate TOF1 values
    std::vector<int>::iterator as = unique(tof1_value.begin(), tof1_value.end());
    tof1_value.resize( distance(tof1_value.begin(),as));

    for (int current_tof1 = 0; current_tof1 < 12; ++current_tof1) {
      for (int neighbor_bar = 0; neighbor_bar < 8; ++neighbor_bar) {
        for (int bar_hit = 0; bar_hit < lepton_hit_bars.size(); ++bar_hit) {
          if (channel[current_tof1][neighbor_bar] == lepton_hit_bars[bar_hit] &&
            find(lepton_ends.begin(), lepton_ends.end(),
              lepton_hit_bars[bar_hit]) == lepton_ends.end()) {
            lepton_ends.push_back(lepton_hit_bars[bar_hit]);
          }
        }
      }
    }

    cout << "Lepton ends: ";
    for (int i = 0; i < lepton_ends.size(); ++i) cout << lepton_ends[i] << " " ;
    // return vector_lepton_kaon;

    //Calculate Distances between Lepton Bars to determine End Points
    if (lepton_ends.size() >= 2) {
      cout << endl << "More than 2 leptons. Finding end leptons of 2 tracks" << endl;
      string bitmask(2, 1); // K leading 1's
      bitmask.resize(lepton_ends.size(), 0); // N-K trailing 0's

        // For all pairs of lepton hits that intersect kaon bars
      map<double, vector<int>> lepton_endpoint_distances;
        do {
            vector<int> lepton_distance;
            for (int i = 0; i < lepton_ends.size(); ++i) // [0..N-1] integers
            {
                if (bitmask[i]) lepton_distance.push_back(lepton_ends[i]);
            }
            // Calculate their distance
            lepton_endpoint_distances[distance(Xloc[lepton_distance[0]],
              Yloc[lepton_distance[0]],
              Xloc[lepton_distance[1]], Yloc[lepton_distance[1]])] = lepton_distance;

        } while (prev_permutation(bitmask.begin(), bitmask.end()));


    } else {
      cout << "Error! Not a typical Y event\n\n\n\n\n\n\n";
      cout << "Event Number: " << ievt << endl;
      cout << "reduced_loss: " << reduced_loss;
    }

    map<double, pair<double, double> > k_stop_candidates;



    cout << endl << "Testing k-stop Candidates" << endl;
    double k_stop_x = -dv_y*10;
    double k_stop_y = dv_x*10;
    if (intersection_points.size() > 0) {
      cout << "# of intersection points: " << intersection_points.size() << endl;
      // return vector_lepton_kaon;
      for (map<int, vector<int>>::iterator junction_lepton = intersection_points.begin();
        junction_lepton != intersection_points.end(); ++junction_lepton) {
          // For each junction lepton, find the neighboring kaons

          //candidate_stop.push_back(junction_lepton);
          vector<int> candidate_stop = junction_lepton->second;

          candidate_stop.push_back(junction_lepton->first);

          double k_stop_x_candidate = 0; double k_stop_y_candidate = 0;
          for (int bar = 0; bar < candidate_stop.size(); ++bar) {
            cout << bar << ", " << Xloc[candidate_stop[bar]] << ", " << Yloc[candidate_stop[bar]] << endl;
            k_stop_x_candidate = k_stop_x_candidate + Xloc[candidate_stop[bar]];
            k_stop_y_candidate = k_stop_y_candidate + Yloc[candidate_stop[bar]];
          }

          // Average these to get k-stop points
          k_stop_x_candidate = k_stop_x_candidate / candidate_stop.size();
          k_stop_y_candidate = k_stop_y_candidate / candidate_stop.size();
          cout << "This candidate kstop is: " << k_stop_x_candidate << ", " << k_stop_y_candidate << endl;
          double k_stop_chiS = 0.0;

          for (int end = 0; end < lepton_ends.size(); ++end) {
            // Draw lines between candidate K-stop and end-points
            double lepton_end_x = Xloc[lepton_ends[end]];
            double lepton_end_y = Yloc[lepton_ends[end]];
            double y_int = points_to_y_int(k_stop_x_candidate, k_stop_y_candidate, lepton_end_x, lepton_end_y);
            double slope = points_to_slope(k_stop_x_candidate, k_stop_y_candidate, lepton_end_x, lepton_end_y);
            cout << y_int << ", " << slope << endl;
            TF1 *k_to_lepton = new TF1("k_to_lepton", "[0]*x+[1]", -50., 50.);
            k_to_lepton->SetParameter(0, slope);
            k_to_lepton->SetParameter(1, y_int);
            //Find ChiSquare of these lines
            k_stop_chiS = k_stop_chiS + gr_Leptons->Chisquare(k_to_lepton);
          }
          pair<double, double> k_candidate_coordinates = make_pair(k_stop_x_candidate, k_stop_y_candidate);
          k_stop_candidates[k_stop_chiS] = k_candidate_coordinates;
          cout << "The kstop candidates are " << endl;
          for(auto it = k_stop_candidates.cbegin(); it != k_stop_candidates.cend(); ++it)
            {
                cout << it->second.first << " " << it->second.second << "\n";
            }
          }
          cout << endl;
          k_stop_x = (k_stop_candidates.cbegin()->second).first;
          k_stop_y = (k_stop_candidates.cbegin()->second).second;



    }



    cout << "The k stop are: " << k_stop_x << ", " << k_stop_y << endl;


    double m1 = 1; double b1 = 0; double m2 = -1; double b2 = 0;
    double first_lepton_x = vec_xx_lepton[0];
    double first_lepton_y = vec_yy_lepton[0];
    double second_lepton_x = vec_xx_lepton.back();
    //all events analyzed for 2 tracks have non-empty leptons
    double second_lepton_y = vec_yy_lepton.back();
    if (lepton_ends.size() > 1) {
      // lepton_ends[0] = 42;
      // lepton_ends[1] = 74;
      first_lepton_x = Xloc[lepton_ends[0]];
      first_lepton_y = Yloc[lepton_ends[0]];
      second_lepton_x = Xloc[lepton_ends[1]];
      second_lepton_y = Yloc[lepton_ends[1]];
      m1 = points_to_slope(k_stop_x, k_stop_y, first_lepton_x, first_lepton_y);
      m2 = points_to_slope(k_stop_x, k_stop_y, second_lepton_x, second_lepton_y);
      b1 = points_to_y_int(k_stop_x, k_stop_y, first_lepton_x, first_lepton_y);
      b2 = points_to_y_int(k_stop_x, k_stop_y, second_lepton_x, second_lepton_y);

    }



    // return vector_lepton_kaon;


    vector<vector<double>> lepton_coordinates;

    for (int i = 0; i < vec_xx_lepton.size(); ++i)
      {
          // save the pair
          vector<double> input;
          input.push_back(vec_xx_lepton[i]);
          input.push_back(vec_yy_lepton[i]);
          lepton_coordinates.push_back(input);
      }


    two_track_parameter_vector params;
    params = m1, b1, m2, b2;
    double ssr = 0.0;
    for (int i = 0; i < lepton_coordinates.size(); ++i) {
      ssr = ssr + pow(calculate_2_track_residual(lepton_coordinates[i],params),2);
    }
    // return vector_lepton_kaon;

    vector<int> vec_l1_initial; vec_l1_initial.clear();
    vector <int> vec_l2_initial; vec_l2_initial.clear();

    for (int i = 0; i < lepton_hit_bars.size(); ++i) {
        int bar = lepton_hit_bars[i];
        double ssr_1 = abs(m1*Xloc[bar] + b1 - Yloc[bar]); //first SSR
        double ssr_2 = abs(m2*Xloc[bar] + b2 - Yloc[bar]); //first SSR
        // double distance_to_l1 = distance_to_line(Xloc[bar], Yloc[bar],  params(0), params(1));
        // double distance_to_l2 = distance_to_line(Xloc[bar], Yloc[bar],  params(2), params(3));

        if (ssr_1 <= ssr_2) {
          vec_l1_initial.push_back(bar);
        }
        else {
          vec_l2_initial.push_back(bar);
        }

      }

    cout << "Plotting initial guess" << endl;



    if (savefile==1) {
      TCanvas *c_initial_params = new TCanvas("Initial Track","Initial Track",0,200,350,350);
      c_initial_params->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
      c_initial_params->Divide(1,1);
      const char *track1_initial = "Track 1";
      const char *track2_initial = "Track 2";
      Subplot_4_test(c_initial_params, vec_l1_initial, 1, track1_initial, 2, m1, b1);
      Subplot_4_test(c_initial_params, vec_l2_initial, 2, track2_initial, 8, m2, b2);
      Subplot_6_test(c_initial_params, vec_l1_initial, vec_l2_initial, "Both tracks", 3, 2, 8,
        m1, b1, m2, b2, dv_x, dv_y, k_stop_x, k_stop_y, intersection_points.size());

      char initial_guess[100];
      sprintf(initial_guess, "./many_tracks_many_bars_initial_guess_images/RUN_%d_event_%d_initial_guess.png", Run_Number, ievt);
      c_initial_params->SaveAs(initial_guess);
    }

    dlib::solve_least_squares(dlib::objective_delta_stop_strategy(1e-7),
                               calculate_2_track_residual,
                               dlib::derivative(calculate_2_track_residual),
                               lepton_coordinates,
                               params);

    vector<int> vec_l1; vec_l1.clear();
    vector <int> vec_l2; vec_l2.clear();

    // //Determine Primary Track
    // // k_stop_x = dv_x;
    // // k_stop_y = dv_y;

    double loss = 0.0;
    for (int i = 0; i < lepton_coordinates.size(); ++i) {
      loss = loss + pow(calculate_2_track_residual(lepton_coordinates[i],params),2);
    }
    cout << "Total ssr for final pars is " << loss << endl;
    vector_lepton_kaon.two_track_reduced_loss = loss/ndf;
    cout << "Dividing vectors into what lines they are close to" << endl;

    m1 = params(0); b1 = params(1); m2 = params(2); b2 = params(3);
    for (int i = 0; i < lepton_hit_bars.size(); ++i) {
        int bar = lepton_hit_bars[i];
        double ssr_1 = abs(m1*Xloc[bar] + b1 - Yloc[bar]); //first SSR
        double ssr_2 = abs(m2*Xloc[bar] + b2 - Yloc[bar]); //first SSR
        // double distance_to_l1 = distance_to_line(Xloc[bar], Yloc[bar],  params(0), params(1));
        // double distance_to_l2 = distance_to_line(Xloc[bar], Yloc[bar],  params(2), params(3));

        if (ssr_1 <= ssr_2) {
          vec_l1.push_back(bar);
        }
        else {
          vec_l2.push_back(bar);
        }

      }

    vector_lepton_kaon.first_track_ChiS = fit_to_bar_numbers(vec_l1, m1, b1)/(vec_l1.size()-2);
    vector_lepton_kaon.second_track_ChiS = fit_to_bar_numbers(vec_l2, m2, b2)/(vec_l2.size()-2);

    TCanvas *cTEST;
    if (savefile == 1) {
      cTEST = new TCanvas("Testing","Testing",0,200,1050,350);
      cTEST->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");
      cTEST->Divide(3,1);
      const char *track1_test = "Track 1";
      const char *track2_test = "Track 2";

      Subplot_4_test(cTEST, vec_l1, 1, track1_test, 2, m1, b1);
      Subplot_4_test(cTEST, vec_l2, 2, track2_test, 8, m2, b2);
      cout << "# of Intersection Points: " << intersection_points.size() << endl;
      Subplot_6_test(cTEST, vec_l1, vec_l2, track2_test, 3, 2, 8,
        m1, b1, m2, b2, dv_x, dv_y, k_stop_x, k_stop_y, intersection_points.size());

      char before_reassignment[100];
      sprintf(before_reassignment, "./many_tracks_many_bars_initial_guess_images/RUN_%d_event_%d_before_reassignment.png", Run_Number, ievt);
      cTEST->SaveAs(before_reassignment);
    }

    // Fit 2nd track
    // Track reassignment
    for (int i = 0; i < vec_l1.size(); ++i) {
      cout << vec_l1[i] << " ";
    }
    cout << endl;
    vector<int> second_track;
    int second_track_num;
    if (vec_l1.size() <= vec_l2.size()) {

       if (vec_l2.size() - vec_l1.size() >= 3) {
         cout << "Primary track l2 is 3+ bars longer than l1" << endl;
         for (auto i = vec_l1.begin(); i != vec_l1.end();) {

           int bar = *i;

           double distance_to_l2 = distance_to_line(Xloc[bar], Yloc[bar],  m2, b2);
           if (distance_to_l2 < 5) {
             cout << "Removing bar " << bar << " from l1" << endl;
             vec_l2.push_back(bar);

             vec_l1.erase(i);

           }
           else ++i;

         }

       }
       else if (vec_l2.size() - vec_l1.size() > 0) {
         cout << "Primary track l2 is 0+ bars longer than l1" << endl;
         for (auto i = vec_l1.begin(); i != vec_l1.end();) {
           int bar = *i;
           double distance_to_l2 = distance_to_line(Xloc[bar], Yloc[bar],  m2, b2);
           if (distance_to_l2 < 2) {
             cout << "Removing bar " << bar << " from l1" << endl;
             vec_l2.push_back(bar);

             vec_l1.erase(i);

           }
           else ++i;

         }

       }
       second_track_num = 1;
        second_track = vec_l1;
        // cout << endl << "Secondary track is line 1";
     }
    else {

      if (vec_l1.size() - vec_l2.size() >= 3) {
        cout << "Primary track l1 is 3+ bars longer than l2" << endl;
        for (auto i = vec_l2.begin(); i != vec_l2.end();) {

          int bar = *i;
          // cout << "Considering bar " << bar << endl;
          double distance_to_l1 = distance_to_line(Xloc[bar], Yloc[bar],  m1, b1);
          // cout << "This bar is " << distance_to_l1 << " mm away from l1" << endl;
          if (distance_to_l1 < 5) {
            cout << "Removing bar " << bar << " from l2" << endl;
            vec_l1.push_back(bar);
            vec_l2.erase(i);

          }
          else ++i;
        }
      }
      else if (vec_l1.size() - vec_l2.size() > 0) {
        cout << "Primary track l1 is 0+ bars longer than l2" << endl;
        for (auto i = vec_l2.begin(); i != vec_l2.end();) {
          int bar = *i;
          double distance_to_l1 = distance_to_line(Xloc[bar], Yloc[bar],  m1, b1);
          if (distance_to_l1 < 2) {
            cout << "Removing bar " << bar << " from l1" << endl;
            vec_l1.push_back(bar);
            vec_l2.erase(i);

          }
          else ++i;
        }
      }
      second_track_num = 2;
      second_track = vec_l2;
      // cout << "Secondary track is line 2";
    }
    // return;
    cout << "Vector 1 size, Vector 2 size:";
    cout << endl << vec_l1.size() << ", " << vec_l2.size() << endl;
    // // return;
    // cout << "VecL1 ";
    // for (int i = 0; i != vec_l1.size(); ++i) {
    //     cout << vec_l1[i] << " " << endl;
    //   }
    // cout << endl << "VecL2 ";
    //   for (int i = 0; i != vec_l2.size(); ++i) {
    //       cout << vec_l2[i] << " " << endl;
    //     }


    // if (vec_l1.size() <= vec_l2.size()) {
    //   m1 = gr_Leptons_2nd_track_fit->GetParameter(1);
    //   b1 = gr_Leptons_2nd_track_fit->GetParameter(0);
    // }
    // else {
    //   double m_fit_track2 = gr_Leptons_2nd_track_fit->GetParameter(1);
    //   double b_fit_track2  = gr_Leptons_2nd_track_fit->GetParameter(0);
    // }

    // m1 = 3.91795; b1 =  -48.2186;
    // m2 = 0.573396; b2 = -13.0612;


    //How many bars are near the


    if (Switch_Display == 1) {

      const char *track1 = "Track 1";
      const char *track2 = "Track 2";
      const char *both_tracks = "Both Tracks";

      Subplot_4(c2, vec_l1, 4, track1, 2);
      Subplot_4(c2, vec_l2, 5, track2, 8);
      Subplot_6(vec_l1, vec_l2, both_tracks, 6, 2, 8, m1, b1, m2, b2,
        second_track_num, dv_x, dv_y, k_stop_x, k_stop_y,  intersection_points.size());



      if (savefile == 1) {
        char image_name[100];
        sprintf(image_name, "./many_tracks_many_bars_initial_guess_images/RUN_%d_event_%d_many_tracks.png", Run_Number, ievt);
        c2->SaveAs(image_name);

      }
    }

  }

  return vector_lepton_kaon;
}

//void FindManyTracks(Int_t Run_Number, Int_t ievt)

  //Plot reconstructed tracks





  //Plot reconstructed tracks
