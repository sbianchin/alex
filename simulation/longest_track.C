// #include "GlobalParameters.h"
// #include "G4Parameters.h"
// #include <stdio.h>

vector<int> longest_nonborder_track_aux(vector<int> lepton_bars, int current_vertex,
  vector<int> max_len_path = {}, vector<int> path_so_far = {}) {

    // cout << "Considering current vertex: " << current_vertex << endl;
    path_so_far.push_back(current_vertex);


    if (path_so_far.size() > max_len_path.size()) max_len_path = path_so_far;
    for (int neighbor = 0; neighbor < 8; ++neighbor) {
      vector<int> new_branch = path_so_far;


      if (find(lepton_bars.begin(), lepton_bars.end(),
        TARGET_neighbours[current_vertex][neighbor]) != lepton_bars.end() &&
          find(path_so_far.begin(), path_so_far.end(),
            TARGET_neighbours[current_vertex][neighbor]) == path_so_far.end()
      ) {
          int new_vertex = TARGET_neighbours[current_vertex][neighbor];
          max_len_path = longest_nonborder_track_aux(lepton_bars, new_vertex, max_len_path, new_branch);
        }

    }
    return max_len_path;
  }

vector<int> longest_nonborder_track(vector<int> lepton_bars) {
  vector<int> max_len_path = {};
  for (int i = 0; i < lepton_bars.size(); ++i) {
      max_len_path = longest_nonborder_track_aux(lepton_bars, lepton_bars[i], max_len_path);

  }
  return max_len_path;
}

// int longest_track() {
//   vector<int> lepton_bars = {0,1,8,9,6,7,17,18,30,31,167,168,179,255};
//   vector<int> max_len_path = longest_nonborder_track(lepton_bars);
//   for (int i = 0; i < max_len_path.size(); ++i) cout << max_len_path[i] << ", ";
//   return 0;
// }
