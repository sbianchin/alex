#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include "Many_Tracks_Fit_2.C"
using namespace std;

vector<int> y_events_525 = {742, 1125, 2645, 4470, 6638, 12865, 13277,
13824, 17058, 17589, 21475, 26300, 26533, 31759, 32619, 36079, 39247,
39347, 42139, 43789, 49618};
void Y_Events_Extract(int Run_Number) {
  //gROOT->LoadMacro("Experimental_Event_Display_Simulation_1.C");
  char event_track[200];
  sprintf(event_track,"/Users/alng/Documents/Datascience/TREK/offline/alex/simulation/RUN_%d_extract_n.txt",
  Run_Number);
  cout << "File opened:  " << event_track << endl;

  //Read Run_525 File
  ifstream fdat_read(event_track, ios::in);
  int y_event = 0;
  char y_records[200];
  sprintf(y_records,"/Users/alng/Documents/Datascience/TREK/offline/alex/simulation/RUN_%d_y_events.txt",
  Run_Number);
  cout << "File saved:  " << y_records << endl;
  int Read_n, Read_ievt;
  ofstream y_event_extract;
  y_event_extract.open(y_records);
  while (fdat_read.good()) {
    fdat_read >> Read_n >> Read_ievt;
    if (Read_n == y_events_525[y_event]) {
      cout << Read_ievt << endl;
      y_event_extract << Read_ievt << endl;
      ++y_event;

    }



  }
  y_event_extract.close();
}



void Many_Tracks_Fit_Wrapper(int Run_Number) {
  char y_records[200];
  sprintf(y_records,"/Users/alng/Documents/Datascience/TREK/offline/alex/simulation/RUN_%d_y_events.txt",
  Run_Number);
  cout << "File opened:  " << y_records << endl;

  //Read Run_525 File
  ifstream fdat_read(y_records, ios::in);
  int Read_ievt;
  while (fdat_read.good()) {
    fdat_read >> Read_ievt;
    FindTrack(Run_Number, Read_ievt, 1);
  }
}
