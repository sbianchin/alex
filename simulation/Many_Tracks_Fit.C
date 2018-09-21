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
//#include "../../MAcros/CommonParameters.h"
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
#include "Plot_Simulation_Event_Display.h"
#include "Experimental_Event_Display_Simulation_1.h"
#endif

//Defining some data structures
struct Track {
  vector<int> lepton_track;
  int kaon_neighbor = 0;
};



Track depth_first_visit_neighbors(vector<int> lepton_track, vector<int> lepton_hit_bars,
  vector<int> kaon_hit_bars, int bar, vector<int> already_visited) {

  cout << "Visiting " << lepton_track[bar] << endl;

  if (find(already_visited.begin(), already_visited.end(), //avoid
  lepton_track[bar]) != already_visited.end()) {
  Track track;
  track.lepton_track = lepton_track;
  return track;
  }


  for (int neighbor_bar = 0; neighbor_bar < 8; ++neighbor_bar) { //for all neighbors of these
      //Check

      cout << "Visiting neighbor " << TARGET_neighbours[lepton_track[bar]][neighbor_bar] << " ";
      if (find(kaon_hit_bars.begin(), kaon_hit_bars.end(), //avoid
      TARGET_neighbours[lepton_track[bar]][neighbor_bar]) != kaon_hit_bars.end()) {
      cout << "Neighbor is a kaon. Ending traversal." << endl;
      //return track;
      Track track;
      track.kaon_neighbor = 1;
      track.lepton_track = lepton_track;
      return track;
      }

      if (find(lepton_track.begin(), lepton_track.end(), //avoid adding in a vertex already visited
        TARGET_neighbours[lepton_track[bar]][neighbor_bar]) == lepton_track.end() &&
        find(lepton_hit_bars.begin(), lepton_hit_bars.end(),
      TARGET_neighbours[lepton_track[bar]][neighbor_bar]) != lepton_hit_bars.end())
        {
        cout << "Whoa! Adding bar " << TARGET_neighbours[lepton_track[bar]][neighbor_bar] << " to track" << endl;
        lepton_track.push_back(TARGET_neighbours[lepton_track[bar]][neighbor_bar]);
        already_visited.push_back(lepton_track[bar]);
        cout << "Current track is ";
        for (vector<int>::const_iterator i = lepton_track.begin();
          i != lepton_track.end(); ++i)
          cout << *i << ' ';
        cout << endl;
        Track track = depth_first_visit_neighbors(lepton_track, lepton_hit_bars, kaon_hit_bars,
          lepton_track.size() - 1, already_visited);
        if (track.kaon_neighbor == 1) {
          return track;
        }
        cout << endl;
        cout << "Back to visiting " << lepton_track[bar] << endl;
      }
    }
  //return lepton_track;
  Track track;
  track.lepton_track = lepton_track;
  return track;
}

Track breadth_first_visit_neighbors(vector<int> lepton_track, vector<int> lepton_hit_bars,
  vector<int> kaon_hit_bars, int bar, vector<int> already_visited) {

  cout << "Visiting " << lepton_track[bar] << endl;

  if (find(already_visited.begin(), already_visited.end(), //avoid
  lepton_track[bar]) != already_visited.end()) {
  Track track;
  track.lepton_track = lepton_track;
  return track;
  }


  for (int neighbor_bar = 0; neighbor_bar < 8; ++neighbor_bar) { //for all neighbors of these
      //Check

      cout << "Visiting neighbor " << TARGET_neighbours[lepton_track[bar]][neighbor_bar] << " ";
      if (find(kaon_hit_bars.begin(), kaon_hit_bars.end(), //avoid
      TARGET_neighbours[lepton_track[bar]][neighbor_bar]) != kaon_hit_bars.end()) {
      cout << "Neighbor is a kaon. Ending traversal." << endl;
      //return track;
      Track track;
      track.kaon_neighbor = 1;
      track.lepton_track = lepton_track;
      return track;
      }

      if (find(lepton_track.begin(), lepton_track.end(), //avoid adding in a vertex already visited
        TARGET_neighbours[lepton_track[bar]][neighbor_bar]) == lepton_track.end() &&
        find(lepton_hit_bars.begin(), lepton_hit_bars.end(),
      TARGET_neighbours[lepton_track[bar]][neighbor_bar]) != lepton_hit_bars.end())
        {
        cout << "Whoa! Adding bar " << TARGET_neighbours[lepton_track[bar]][neighbor_bar] << " to track" << endl;
        lepton_track.push_back(TARGET_neighbours[lepton_track[bar]][neighbor_bar]);
        already_visited.push_back(lepton_track[bar]);
        cout << "Current track is ";
        for (vector<int>::const_iterator i = lepton_track.begin();
          i != lepton_track.end(); ++i)
          cout << *i << ' ';
        cout << endl;
        //Track track = visit_neighbors(lepton_track, lepton_hit_bars, kaon_hit_bars,
        //  lepton_track.size() - 1, already_visited);
        // if (track.kaon_neighbor == 1) {
          // return track;
        // }
        cout << endl;
        cout << "Back to visiting " << lepton_track[bar] << endl;
      }
    }
  //return lepton_track;
  Track track;
  track.lepton_track = lepton_track;
  return track;
}


void FindTrack(Int_t Run_Number=5, Int_t ievt=0, int current_tof1 = 5,
  int walk=0, Int_t enable_cout=0, Int_t display = 0) {
  //Find Lepton and Kaon Tracks
  Lepton vector_lepton_kaon = Event_Display_5_1(Run_Number, ievt,
    walk, enable_cout, display);

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
  cout << "Lepton Bars Hit" ;
  for (vector<double>::const_iterator i = vec_xx_lepton.begin();
    i != vec_xx_lepton.end(); ++i)
    cout << *i << ' ';
  cout << endl << "Y coords";
  for (vector<double>::const_iterator i = vec_yy_lepton.begin();
    i != vec_yy_lepton.end(); ++i)
    cout << *i << ' ';

  //Define some variables and containers
  //int current_tof1 = ;
  vector < vector<int> > lepton_tracks;
  vector<int> already_visited;

  //Check if bars neighboring hit TOF are hit
  for (int neighbor_bar = 0; neighbor_bar < 8; ++neighbor_bar) {
    for (int bar_hit = 0; bar_hit < lepton_hit_bars.size(); ++bar_hit) {
      //cout << neighbor_bar << " " << bar_hit << endl;
      //cout << channel[current_tof1][neighbor_bar] << " " << lepton_hit_bars[bar_hit] << endl;
      if (channel[current_tof1][neighbor_bar] == lepton_hit_bars[bar_hit]) {
        vector<int> a_track = {lepton_hit_bars[bar_hit]};
        lepton_tracks.push_back(a_track);

      }
    }
  }
  //Reconstruct tracks
  //bool build_track = true;
  int build_track = 0;
  int tracking_lepton = 0;
  cout << "Before the track was build the lepton track is ";
  for (int i = 0; i < lepton_tracks.size(); i++)
    {
        for (int j = 0; j < lepton_tracks[i].size(); j++)
        {
            cout << lepton_tracks[i][j];
        }
    }
  cout << endl;
  cout << "The bars that were leptons are ";
  for (vector<int>::const_iterator i = lepton_hit_bars.begin();
    i != lepton_hit_bars.end(); ++i)
    cout << *i << ' ';
  cout << endl;
  cout << "Before the track was built the already visited vector is ";
  for (vector<int>::const_iterator i = already_visited.begin();
    i != already_visited.end(); ++i)
    cout << *i << ' ';
  cout << endl;

  cout << "Track building" << endl;


  for (int track_number = 0; track_number < lepton_tracks.size(); ++track_number) {

    Track track = depth_first_visit_neighbors(lepton_tracks[track_number], lepton_hit_bars,
      kaon_hit_bars, 0, already_visited);
    cout << "This track is ";
    for (vector<int>::const_iterator i = track.lepton_track.begin();
      i != track.lepton_track.end(); ++i)
      cout << *i << ' ';
    cout << endl;
    cout << "\n\n\n\n\n\n\n" << "OK Going to next track" << endl;
    lepton_tracks[track_number] = track.lepton_track;
    /*for (int bar = tracking_lepton; bar < lepton_track.size(); ++bar) { //for each hit TOF neighbor bar


      if (track.visited == 1) {
        cout << "Bar already visited" << endl;
        continue;
      }
      lepton_track = track.lepton_track;
        cout << "Tracking lepton " << tracking_lepton << endl;
        ++tracking_lepton;
        cout << "Lepton track built so far is ";
        for (vector<int>::const_iterator i = lepton_track.begin();
          i != lepton_track.end(); ++i)
          cout << *i << ' ';
        cout << endl;


        already_visited.push_back(lepton_track[bar]);
        cout << "Already visited ";
        for (vector<int>::const_iterator i = already_visited.begin();
          i != already_visited.end(); ++i)
          cout << *i << ' ';
        cout << endl;
      }
      // build_track = false;
      //++build_track;
    //}
    */
  }
    cout << "After the track was built the lepton tracks are ";

    for (int i = 0; i < lepton_tracks.size(); i++)
      {
          for (int j = 0; j < lepton_tracks[i].size(); j++)
          {
              cout << lepton_tracks[i][j] << " ";
          }
          cout << endl;
      }
    cout << endl;
    /*cout << "Lepton track is ";
    for (vector<int>::const_iterator i = lepton_track.begin();
      i != lepton_track.end(); ++i)
      cout << *i << ' ';
    cout << endl;
  }*/
}

// void FindManyTracks(Int_t Run_Number, Int_t ievt)

  //Plot reconstructed tracks





  //Plot reconstructed tracks
