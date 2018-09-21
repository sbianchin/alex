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
typedef dlib::matrix<double,4,1> parameter_vector;
#include "Plot_Simulation_Event_Display.h"
#include "Experimental_Event_Display_Simulation_1.h"
#endif

//Defining some data structures
struct Track {
  vector<int> lepton_track;
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

void FindTrack(Int_t Run_Number=5, Int_t ievt=0, int savefile = 0; int current_tof1 = 5,
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
  Double_t dv_x = vector_lepton_kaon.dv_x;
  Double_t dv_y = vector_lepton_kaon.dv_y;
  Int_t tof1N[5];
  copy(vector_lepton_kaon.tof1N, vector_lepton_kaon.tof1N+5, tof1N);
  for (vector<int>::const_iterator i = lepton_hit_bars.begin();
    i != lepton_hit_bars.end(); ++i)
    cout << *i << ' ';


  //Calculate K-stop
  //calculate_K_stop(kaon_hit_bars, lepton_hit_bars)

  // Draw 2 tracks
  Coordinate endpoint_1; endpoint_1.x = Xloc[102]; endpoint_1.y = Yloc[102]; // hard-coded
  Coordinate endpoint_2; endpoint_2.x = Xloc[241]; endpoint_2.y = Yloc[241]; // hard-coded
  Coordinate endpoint_3; endpoint_3.x = Xloc[252]; endpoint_3.y = Yloc[252]; // hard-coded

  vector<int> vec_l1;
  vector <int> vec_l2;

  cout << "Dividing vectors into what lines they are close to";

  //Assign bars to each track
  for (vector<int>::const_iterator i = lepton_hit_bars.begin();
    i != lepton_hit_bars.end(); ++i) {
      cout << "Considering bar " <<*i << ' ' << endl;
      int bar = *i;
      cout << "The bar is " << bar;
      double distance_to_l1 = Distance_from_Bar_to_Line(bar, endpoint_1, endpoint_2);
      double distance_to_l2 = Distance_from_Bar_to_Line(bar, endpoint_1, endpoint_3);

      cout << "Distance to Line l1: " << distance_to_l1 << " ";
      cout << "Distance to Line l2: " << distance_to_l2;
      if (distance_to_l1 <= distance_to_l2) { //Bar gets assigned to closer track
        vec_l1.push_back(bar);
        cout << "Pushing to vec_l1";
      }
      else {
        vec_l2.push_back(bar);
        cout << "Pushing to vec_l2";
      }

    }

    // Cheating: Erasing bars for better fit
    cout << endl << "Before erasure vec_l2 is ";
    for (vector<int>::const_iterator i = vec_l2.begin();
      i != vec_l2.end(); ++i) {
        cout << *i << " ";
      }
    // vec_l2.erase(remove(vec_l2.begin(), vec_l2.end(), 120), vec_l2.end());
    // vec_l2.erase(remove(vec_l2.begin(), vec_l2.end(), 155), vec_l2.end());
    // vec_l2.erase(remove(vec_l2.begin(), vec_l2.end(), 244), vec_l2.end());
    // vec_l2.erase(remove(vec_l2.begin(), vec_l2.end(), 252), vec_l2.end());
    // vec_l2.erase(remove(vec_l2.begin(), vec_l2.end(), 233), vec_l2.end());

    cout << endl << "After erasure vec_l2 is ";
    for (vector<int>::const_iterator i = vec_l2.begin();
      i != vec_l2.end(); ++i) {
        cout << *i << " ";
      }

  //Plot tracks
  Subplot_4(vec_l1, 4, 2);
  Subplot_4(vec_l2, 5, 8);
  Subplot_6(vec_l1, vec_l2, 6, 2, 8);

  if (savefile == 1) {
    
  }

}
