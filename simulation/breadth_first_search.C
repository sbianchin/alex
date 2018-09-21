struct path_length { vector<int> path = {}; double length = DBL_MAX; };


path_length dijkstra(const vector<int> &lepton_bars, int source, int target) {
    map<int, path_length> min_distance;
    path_length source_path;
    source_path.path = {source};
	source_path.length = 0;
    min_distance[source] = source_path;
    queue<int> active_vertices;
    active_vertices.push(source);
    map<int, bool> visited;

    while (!active_vertices.empty()) {
      int where = active_vertices.front();
      active_vertices.pop();
      visited[where] = true;
      path_length path_to_here = min_distance[where];

      if (where == target) {
		  return path_to_here;
	  }

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
            double distance_to_adj_node = distance(Xloc[where], Yloc[where], Xloc[adj_node], Yloc[adj_node]);

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

	    // if no neighbors are found, look at all bars in a 12mm (~3 bars
	    // orthogonal, ~2 bars diagonal) radius
	    // this allows the algorithm to jump over gaps in the lepton_bars
        if (no_neighbor) {
          for (unsigned int bar = 0; bar < lepton_bars.size(); ++bar) {
            double distance_to_other_bar = distance(Xloc[where], Yloc[where], Xloc[lepton_bars[bar]], Yloc[lepton_bars[bar]]);
            int adj_node = lepton_bars[bar];
            if (0 < distance_to_other_bar && distance_to_other_bar < PATH_TRAVERSIAL_DIJKSTRA_JUMP_RADIUS
                && min_distance[adj_node].length > min_distance[where].length + distance_to_other_bar
				&& visited[adj_node] == false
			) {
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


std::vector<int> FindShortestPath(const vector<int> &lepton_bars, int source, int target) {
	map<int, path_length> min_distance;
	path_length source_path;
	source_path.path = {source};
	source_path.length = 0;
	min_distance[source] = source_path;
	queue<int> active_vertices;
	active_vertices.push(source);
	map<int, bool> visited;

	while (!active_vertices.empty()) {
		int where = active_vertices.front();
		active_vertices.pop();
		visited[where] = true;
		path_length path_to_here = min_distance[where];

		if (where == target) {
			std::vector<int> lepton_path;
			for (int bar : path_to_here.path) {
				if (find(lepton_bars.begin(), lepton_bars.end(), bar) != lepton_bars.end()) {
					lepton_path.push_back(bar);
				}
			}
			return lepton_path;
		}

		for (int neighbor = 0; neighbor < 8; ++neighbor) {
			int adj_node = TARGET_neighbours[where][neighbor];

			// neigbor field is out of bounds
			if (adj_node < 0 || adj_node > nBars) {
				continue;
			}

			if (visited[adj_node] == false) {
				double distance_to_adj_node = distance(Xloc[where], Yloc[where], Xloc[adj_node], Yloc[adj_node]);

				// if adj_node is not in lepton_bars, penalize the edge
				// this way, the algorithm can jump over gaps in the track of
				// lepton_bars, but will snap back onto the lepton_bars if
				// possible.
				if (find(lepton_bars.begin(), lepton_bars.end(), adj_node) == lepton_bars.end()) {
					distance_to_adj_node *= PATH_TRAVERSIAL_ALL_PENALTY;
				}

				double new_path_length = min_distance[where].length + distance_to_adj_node;
				if (min_distance[adj_node].length > new_path_length) {
					vector<int> path_to_adj_node = path_to_here.path;
					path_to_adj_node.push_back(adj_node);
					min_distance[adj_node].path = path_to_adj_node;
					min_distance[adj_node].length = new_path_length;
					active_vertices.push(adj_node);
				}
			}
		}
	}

	//if no path found, return empty vector
	return std::vector<int>({});
}
