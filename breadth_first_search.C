#include <queue>

struct path_length { vector<int> path = {}; double length = DBL_MAX; };



path_length dijkstra(const vector<int> &lepton_bars, int source, int target) {
    map<int, path_length> min_distance;
    path_length source_path;
    source_path.path = {source}; source_path.length = 0;
    min_distance[ source ] = source_path;
    queue< int > active_vertices;
    active_vertices.push( source );
    map<int, bool> visited;

    // cout << "Source is: " << source << "; Target is: " << target << endl;



    // return source_path;

    while (!active_vertices.empty()) {
      int where = active_vertices.front();
      active_vertices.pop();
      visited[where] = true;
      path_length path_to_here = min_distance[where];
      if (where == target) return path_to_here;

      bool no_neighbor = true;

      for (int neighbor = 0; neighbor < 8; ++neighbor) {
        int adj_node = TARGET_neighbours[where][neighbor];
        // cout << "Considering neighbor: " << adj_node << " ";
        // cout << endl << "Is it true that" << adj_node << " is in bars hit? " << (find(lepton_bars.begin(), lepton_bars.end(),
        //   adj_node) != lepton_bars.end()) << endl;
        // cout << "Is it visited? " << visited[adj_node] << endl;
        if (find(lepton_bars.begin(), lepton_bars.end(),
          adj_node) != lepton_bars.end() &&
            visited[adj_node] == false
        ) {
            no_neighbor = false;
            double distance_to_adj_node = distance(Xloc[where], Yloc[where],
              Xloc[adj_node], Yloc[adj_node]);

            if (min_distance[adj_node].length > min_distance[where].length + distance_to_adj_node) {
                // active_vertices.erase( adj_node );
                vector<int> path_to_adj_node = path_to_here.path;
                path_to_adj_node.push_back(adj_node);
                min_distance[adj_node].path = path_to_adj_node;
                min_distance[adj_node].length = min_distance[where].length + distance_to_adj_node;
                active_vertices.push( adj_node );
                }
          }
        }
        // skipping bars if immediate neighbors not found
        // cout << endl << "Node " << where << " has no neighbor? " << no_neighbor << endl;
        if (no_neighbor) {
          for (unsigned int bar = 0; bar < lepton_bars.size(); ++bar) {
            double distance_to_other_bar = distance(Xloc[where], Yloc[where], Xloc[lepton_bars[bar]], Yloc[lepton_bars[bar]]);
            int adj_node = lepton_bars[bar];
            if (0 < distance_to_other_bar && distance_to_other_bar < 9 &&
                min_distance[adj_node].length > min_distance[where].length + distance_to_other_bar
              && visited[adj_node] == false) {
                  // active_vertices.insert( adj_node );
                  vector<int> path_to_adj_node = path_to_here.path;
                  path_to_adj_node.push_back(adj_node);
                  min_distance[adj_node].path = path_to_adj_node;
                  min_distance[adj_node].length = min_distance[where].length + distance_to_other_bar;
                  active_vertices.push( adj_node );
                  }
            }
        }
    }

    //if no path found, return empty vector
    source_path.path = {};
    return source_path;


}
