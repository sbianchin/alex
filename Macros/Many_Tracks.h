#include "breadth_first_search.C"
#include "compute_largest_triangle.C"

struct TrackVariables {
  double first_track_ChiS, ndf_1, second_track_ChiS, ndf_2, E_positron, between_angle, primary_bars_energy;
  vector<int> primary_track, secondary_track, vec_l1, vec_l2;
};

struct FitParams {

double m, b, ChiS, ndf;

};

FitParams GetFitParams(vector<double> vec_x, vector<double> vec_y,
  vector<double> vec_ex, vector<double> vec_ey) {
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


  pair<vector<int>, vector<int>> stop_at_primary(map <int, vector<int>> lepton_connections_graph,
    vector<int> vec_l1, vector<int> vec_l2) {
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


//Find 2 tracks function
// Input (in order): lepton bars, k_stop x coordinate, k_stop y coordinate,
// tof1_hit x coordinate, tof1_hit y coordinate, hit Gap (necessary for TGraphErrors fit with tof1)
TrackVariables FindTrack(vector<int> lepton_hit_bars, double dv_x, double dv_y,
  double tof1_hit_x, double tof1_hit_y, int Gap) {

  //Determine the 3 points determining the largest triangle_area

  vector<int> largest_triangle = compute_largest_triangle(lepton_hit_bars);

  //Determine the lepton bar that lies closest to k-stop
  double closest_lepton_dist_to_kstop = 100000; int closest_lepton_to_kstop;
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
  for (int i = 0; i < lepton_hit_bars.size(); ++i) {
    double dist_to_k_stop = distance(Xloc[lepton_hit_bars[i]],
                  Yloc[lepton_hit_bars[i]], dv_x*10, dv_y*10);
    if (dist_to_k_stop < closest_lepton_dist_to_kstop) {
      closest_lepton_dist_to_kstop = dist_to_k_stop;
      closest_lepton_to_kstop = lepton_hit_bars[i];
    }
  }



  //First track (bar numbers stored in vec_l1) is one endpoint back to k-stop bar
  path_length traversal = dijkstra(lepton_hit_bars, lepton_ends[0], closest_lepton_to_kstop);
  vector<int> vec_l1 = traversal.path;


  //Second track (bar numbers stored in vec_l2) is the other endpoint back to k-stop bar
  traversal = dijkstra(lepton_hit_bars, lepton_ends[1], closest_lepton_to_kstop);
  vector<int> vec_l2 = traversal.path;


  //Reassign the rest of the bars
  map <int, vector<int>> lepton_connections_graph;
  map<int, bool> already_assigned;
  queue<int> track_assign;
  //Build graph of lepton bars and their edges
  for (int i = 0; i < lepton_hit_bars.size(); ++i) {
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
  for (int i = 0; i < vec_l1.size(); ++i) already_assigned[vec_l1[i]] = true;
  for (int i = 0; i < vec_l2.size(); ++i) already_assigned[vec_l2[i]] = true;
  for (int i = 0; i < vec_l1.size(); ++i) {

      vector<int> neighbors = lepton_connections_graph[vec_l1[i]];
      for (int neighbor = 0; neighbor < neighbors.size(); ++neighbor) {
        if (already_assigned[neighbors[neighbor]] == false) {
              track_assign.push(neighbors[neighbor]);
              //Check if the adjacent node only borders vec_l1
              vector<int> adj_node_neighbors = lepton_connections_graph[neighbors[neighbor]];
              bool only_borders_l1 = true;
              for (int j = 0; j < adj_node_neighbors.size(); ++j) {
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
  for (int i = 0; i < vec_l2.size(); ++i) {
      // cout << "\nCurrent bar is: " << vec_l2[i] << "\n";
      // cout << "The neighbors are: ";

      vector<int> neighbors = lepton_connections_graph[vec_l2[i]];
      // for (int j = 0; j < neighbors.size(); ++j) cout << neighbors[j] << ", ";
      for (int neighbor = 0; neighbor < neighbors.size(); ++neighbor) {
          //If not already assigned, then push them into track assign queue
          if (already_assigned[neighbors[neighbor]] == false) {
              track_assign.push(neighbors[neighbor]);
              //Check if the adjacent node only borders vec_l1
              vector<int> adj_node_neighbors = lepton_connections_graph[neighbors[neighbor]];
              bool only_borders_l2 = true;
              for (int j = 0; j < adj_node_neighbors.size(); ++j) {
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
  for (int i = 0; i < vec_l1.size(); ++i) {
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
  for (int i = 0; i < vec_l2.size(); ++i) {
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
  traversal = dijkstra(lepton_hit_bars, lepton_ends[0], lepton_ends[1]);
  vector<int> end_to_end = traversal.path;
  // cout << "Path has length: " << traversal.length << endl;
  // cout << "Traversing from one end to another: ";
  // for (int i = 0; i < traversal.path.size(); ++i) {
  //  cout << traversal.path[i] << ", ";
  // }
  cout << endl;

  //If there's no path, then not a y-event
  bool y_event = true;
  if (end_to_end.size() == 0) { cout << "Not a Y-event!"; y_event = false;}
  if (vec_l1.size() == 0 || vec_l2.size() == 0) { cout << "Not a Y-event!"; y_event = false;}

  //Check for v_shape
  // if (vec_l1.size() != 0 && vec_l2.size() != 0 && end_to_end.size() != 0) {
  //   bool v_shape = true;
  //
  //   int array_length = end_to_end.size() - 1;
  //   int horizontal_shape[array_length];
  //   int vertical_shape[array_length];
  //   cout << "\nHorizontal shape: ";
  //   for (int i = 0; i < array_length; ++i) cout << horizontal_shape[i] << ", ";
  //   cout << "\nVertical  shape: ";
  //   for (int i = 0; i < array_length; ++i) cout << vertical_shape[i] << ", ";
  //
  //
  //   for (int i = 0; i < array_length; ++i) {
  //     horizontal_shape[i] = 0; vertical_shape[i] = 0;
  //     if (Xloc[end_to_end[i+1]] - Xloc[end_to_end[i]] < 0) horizontal_shape[i] = -1;
  //     if (0 < Xloc[end_to_end[i+1]] - Xloc[end_to_end[i]]) horizontal_shape[i] = 1;
  //     if (Xloc[end_to_end[i+1]] - Xloc[end_to_end[i]] == 0) {
  //       if (i > 0 && horizontal_shape[i - 1] != 0) horizontal_shape[i] = horizontal_shape[i - 1];
  //     }
  //
  //     if (Yloc[end_to_end[i+1]] - Yloc[end_to_end[i]] < 0) vertical_shape[i] = -1;
  //     if (0 < Yloc[end_to_end[i+1]] - Yloc[end_to_end[i]]) vertical_shape[i] = 1;
  //     if (Yloc[end_to_end[i+1]] - Yloc[end_to_end[i]] == 0) {
  //       if (i > 0 && vertical_shape[i - 1] != 0) vertical_shape[i] = vertical_shape[i - 1];
  //     }
  //   }
  //
  //   int shape = 0;
  //   bool all_zeros_horizontal = false;
  //   while (horizontal_shape[shape] == 0 && shape < array_length - 1) {
  //     if (shape == array_length - 2 && horizontal_shape[shape + 1] == 0) {
  //       all_zeros_horizontal = true;
  //     }
  //     if (horizontal_shape[shape + 1] != 0) {
  //       for (;shape >= 0; --shape) horizontal_shape[shape] = horizontal_shape[shape + 1];
  //       break;
  //     }
  //     ++shape;
  //   }
  //
  //   shape = 0;
  //   bool all_zeros_vertical = false;
  //   cout << "\narray_length: " << array_length << "\n";
  //
  //   while (vertical_shape[shape] == 0 && shape < array_length - 1) {
  //     if (shape == array_length - 2 && vertical_shape[shape + 1] == 0) {
  //       all_zeros_vertical = true;
  //     }
  //     if (vertical_shape[shape + 1] != 0) {
  //
  //       for (;shape >= 0; --shape) vertical_shape[shape] = vertical_shape[shape + 1];
  //
  //       break;
  //     }
  //     ++shape;
  //   }
  //
  //   cout << "\nHorizontal shape: ";
  //   for (int i = 0; i < array_length; ++i) cout << horizontal_shape[i] << ", ";
  //   cout << "\nVertical  shape: ";
  //   for (int i = 0; i < array_length; ++i) cout << vertical_shape[i] << ", ";
  //
  //   //If v-shaped, then is y-event
  //
  //   int horizontal_total = 0; int vertical_total = 0;
  //   for (int i = 0; i < array_length; ++i) {
  //     horizontal_total += horizontal_shape[i];
  //     vertical_total += vertical_shape[i];
  //   }
  //   if ((abs(horizontal_total) == array_length || all_zeros_horizontal)
  //     && (abs(vertical_total) == array_length || all_zeros_vertical)) {
  //     cout << "Not a y-event";
  //     v_shape = false;
  //     y_event = false;
  //     // How to Handle vertical and horizontal events?
  //   }
  //
  //   //Determine the density of the path between two ends
  //   if (!v_shape) {
  //       vector<int> density, border;
  //       // Determine the highest and lowest X, Y in path endpoints
  //       double max_X, min_X, max_Y, min_Y;
  //
  //       if (Xloc[lepton_ends[0]] <= Xloc[lepton_ends[1]]) {
  //         max_X = Xloc[lepton_ends[1]]; min_X = Xloc[lepton_ends[0]];
  //       }
  //       else {
  //         max_X = Xloc[lepton_ends[0]]; min_X = Xloc[lepton_ends[1]];
  //       }
  //
  //       if (Yloc[lepton_ends[0]] <= Yloc[lepton_ends[1]]) {
  //         max_Y = Yloc[lepton_ends[1]]; min_Y = Yloc[lepton_ends[0]];
  //       }
  //       else {
  //         max_Y = Yloc[lepton_ends[0]]; min_Y = Yloc[lepton_ends[1]];
  //       }
  //
  //
  //       for (int i = 0; i < 256; ++i) {
  //         // If lepton bars lie within the box defined by highest, lowest X, Y
  //         if (min_X <= Xloc[i] <= max_X && min_Y <= Yloc[i] <= max_Y) {
  //           density.push_back(i);
  //         }
  //         if (min_X <= Xloc[i] <= max_X && min_Y <= Yloc[i] <= max_Y) {
  //           density.push_back(i);
  //         }
  //         // Check lepton bars on the borders
  //         if (horizontal_total > 0 && vertical_total > 0
  //             && Xloc[i] == min_X && Yloc[i] == max_Y) {
  //             border.push_back(i);
  //           }
  //         if (horizontal_total > 0 && vertical_total < 0
  //             && Xloc[i] == min_X && Yloc[i] == max_Y) {
  //           border.push_back(i);
  //           }
  //         if (horizontal_total < 0 && vertical_total > 0
  //             && Xloc[i] == min_X && Yloc[i] == max_Y) {
  //           border.push_back(i);
  //           }
  //       }
  //   }
  // }
  // vector_lepton_kaon.y_event = y_event;












  // Distinguish between primary and Secondary Track

  vector<int> primary_track, secondary_track;
  map<int, bool> visited;
  double second_track_ChiS, first_track_ChiS;
  int ndf_1, ndf_2;



  //Determine Chi-sq of modified first track found
  if (vec_l1.size() > 0) {
    vector<double> vec_l1_leptons_x, vec_l1_leptons_y;
    vector<double> vec_l1_leptons_ex, vec_l1_leptons_ey;
    for (vector<int>::const_iterator i = vec_l1.begin();
      i != vec_l1.end(); ++i) {
        int bar = *i;
        vec_l1_leptons_x.push_back(Xloc[bar]);
        vec_l1_leptons_y.push_back(Yloc[bar]);
        vec_l1_leptons_ex.push_back(TARGET_Errors_X);
        vec_l1_leptons_ey.push_back(TARGET_Errors_Y);
      }


    //Add in highest energy TOF1 hit, weighted by 3
    for (int i = 0; i < 3; ++i) {
      vec_l1_leptons_x.push_back(tof1_hit_x);
      vec_l1_leptons_y.push_back(tof1_hit_y);
      vec_l1_leptons_ex.push_back(TOF1_Errors_X[Gap][0]);
      vec_l1_leptons_ey.push_back(TOF1_Errors_Y[Gap][0]);
    }

    FitParams fit_params_l1 = GetFitParams(vec_l1_leptons_x, vec_l1_leptons_y, vec_l1_leptons_ex, vec_l1_leptons_ey);
    first_track_ChiS = fit_params_l1.ChiS;
    ndf_1 = fit_params_l1.ndf;


  }

  //Determine Chi-sq of modified second track found
  if (vec_l2.size() > 0) {
    vector<double> vec_l2_leptons_x, vec_l2_leptons_y;
    vector<double> vec_l2_leptons_ex, vec_l2_leptons_ey;

    for (vector<int>::const_iterator i = vec_l2.begin();
      i != vec_l2.end(); ++i) {
        int bar = *i;
        vec_l2_leptons_x.push_back(Xloc[bar]);
        vec_l2_leptons_y.push_back(Yloc[bar]);
        vec_l2_leptons_ex.push_back(TARGET_Errors_X);
        vec_l2_leptons_ey.push_back(TARGET_Errors_Y);
      }
    //Add in highest energy TOF1 hit, weighted by 3
    for (int i = 0; i < 3; ++i) {
      vec_l2_leptons_x.push_back(tof1_hit_x);
      vec_l2_leptons_y.push_back(tof1_hit_y);
      vec_l2_leptons_ex.push_back(TOF1_Errors_X[Gap][0]);
      vec_l2_leptons_ey.push_back(TOF1_Errors_Y[Gap][0]);
    }


    FitParams fit_params_l2 = GetFitParams(vec_l2_leptons_x, vec_l2_leptons_y, vec_l2_leptons_ex, vec_l2_leptons_ey);
    second_track_ChiS = fit_params_l2.ChiS;
    ndf_2 = fit_params_l2.ndf;


  }




  double m1, m2;
  // The track with lower reduced Chi-sq is the primary track
  if (first_track_ChiS/ndf_1  < second_track_ChiS/ndf_2) {


    pair<vector<int>, vector<int>> truncated_secondary = stop_at_primary(lepton_connections_graph, vec_l1, vec_l2);
    primary_track = truncated_secondary.first;
    secondary_track = truncated_secondary.second;


    // Determine ChiSq of both tracks
    vector<double> leptons_x1, leptons_y1;
    vector<double> leptons_ex1, leptons_ey1;

    for (vector<int>::const_iterator i = primary_track.begin();
      i != primary_track.end(); ++i) {
        int bar = *i;
        leptons_x1.push_back(Xloc[bar]);
        leptons_y1.push_back(Yloc[bar]);
        leptons_ex1.push_back(TARGET_Errors_X);
        leptons_ey1.push_back(TARGET_Errors_Y);
      }


    FitParams fit_params_primary = GetFitParams(leptons_x1, leptons_y1, leptons_ex1, leptons_ey1);
    first_track_ChiS = fit_params_primary.ChiS;
    ndf_1 = fit_params_primary.ndf;
    m1 = fit_params_primary.m;

    vector<double> leptons_x2, leptons_y2;
    vector<double> leptons_ex2, leptons_ey2;

    for (vector<int>::const_iterator i = secondary_track.begin();
      i != secondary_track.end(); ++i) {
        int bar = *i;
        leptons_x2.push_back(Xloc[bar]);
        leptons_y2.push_back(Yloc[bar]);
        leptons_ex2.push_back(TARGET_Errors_X);
        leptons_ey2.push_back(TARGET_Errors_Y);
      }

    FitParams fit_params_secondary = GetFitParams(leptons_x2, leptons_y2, leptons_ex2, leptons_ey2);
    second_track_ChiS = fit_params_secondary.ChiS;
    ndf_2 = fit_params_secondary.ndf;
    m2 = fit_params_secondary.m;


  }
  else {

    pair<vector<int>, vector<int>> truncated_secondary = stop_at_primary(lepton_connections_graph, vec_l2, vec_l1);
    primary_track = truncated_secondary.first;
    secondary_track = truncated_secondary.second;


    // Determine ChiSq of both tracks
    vector<double> leptons_x1, leptons_y1;
    vector<double> leptons_ex1, leptons_ey1;

    for (vector<int>::const_iterator i = primary_track.begin();
      i != primary_track.end(); ++i) {
        int bar = *i;
        leptons_x1.push_back(Xloc[bar]);
        leptons_y1.push_back(Yloc[bar]);
        leptons_ex1.push_back(TARGET_Errors_X);
        leptons_ey1.push_back(TARGET_Errors_Y);
      }


    FitParams fit_params_primary = GetFitParams(leptons_x1, leptons_y1, leptons_ex1, leptons_ey1);
    first_track_ChiS = fit_params_primary.ChiS;
    ndf_1 = fit_params_primary.ndf;
    m1 = fit_params_primary.m;

    vector<double> leptons_x2, leptons_y2;
    vector<double> leptons_ex2, leptons_ey2;

    for (vector<int>::const_iterator i = secondary_track.begin();
      i != secondary_track.end(); ++i) {
        int bar = *i;
        leptons_x2.push_back(Xloc[bar]);
        leptons_y2.push_back(Yloc[bar]);
        leptons_ex2.push_back(TARGET_Errors_X);
        leptons_ey2.push_back(TARGET_Errors_Y);
      }

    FitParams fit_params_secondary = GetFitParams(leptons_x2, leptons_y2, leptons_ex2, leptons_ey2);
    second_track_ChiS = fit_params_secondary.ChiS;
    ndf_2 = fit_params_secondary.ndf;
    m2 = fit_params_secondary.m;


  }

  if (primary_track.size() == 0) {
    primary_track = secondary_track;
    secondary_track = {};
  }

  // Determine angle between lines of best fit of 2 tracks
  double between_angle = abs(atan(m1) - atan(m2));
  if (between_angle > TMath::Pi()/2) between_angle = TMath::Pi() - between_angle;
  between_angle = between_angle * 180 / TMath::Pi(); // Convert to degrees





  // double primary_bars_energy = 0;
  // for (int i = 0; i < primary_track.size(); ++i) {
  //   primary_bars_energy += targetdE[primary_track[i]];
  // }
  // double E_positron = tof1E[0] + delta_energy + primary_bars_energy;



  TrackVariables two_tracks;
  two_tracks.first_track_ChiS = first_track_ChiS;
  two_tracks.ndf_1 = ndf_1;
  two_tracks.second_track_ChiS = second_track_ChiS;
  two_tracks.ndf_2 = ndf_2;
  // two_tracks.E_positron = E_positron;
  // two_tracks.primary_bars_energy = primary_bars_energy;
  two_tracks.between_angle = between_angle;
  two_tracks.vec_l1 = vec_l1;
  two_tracks.vec_l2 = vec_l2;
  two_tracks.primary_track = primary_track;
  two_tracks.secondary_track = secondary_track;


  return two_tracks;



}
