#include "Many_Tracks_Fit_2.2.C"

void Batch_Distance_Fit(int Run_Number, int flag, int min_evt = 0, int max_evt = 0) {

  ofstream fout;


  char Name_finput[200];
  sprintf(Name_finput,"G4Run%d.root", Run_Number);
  cout << "File opened:  " << Name_finput << endl;
  TChain fChain("Kaon");
  fChain.Add(Name_finput);
  Int_t nentries = (Int_t)fChain.GetEntries();

  if (flag == 0) {
    char batch_output[50];
    sprintf(batch_output, "RUN_%d_Batch_Distance_Fit.txt", Run_Number);
    fout.open(batch_output);
    int i = 0;
    while (i < nentries) {

      Lepton vector_lepton_kaon = MC_Event_Display(Run_Number, i, 1);
      fout << i << ", " << vector_lepton_kaon.reduced_loss << ", "
        << vector_lepton_kaon.reduced_ChiS << ", "<< endl;
      ++i;
    }
  }
  else if (flag == 1) {
    char batch_output[50];
    sprintf(batch_output, "RUN_%d_Batch_Distance_Fit_Y_events.txt", Run_Number);
    fout.open(batch_output);
    char events[200];
    sprintf(events,
      "/Users/alng/Documents/Datascience/TREK/offline/alex/simulation/RUN_%d_read_distance_fit.txt",
    Run_Number);
    ifstream fdat_read(events, ios::in);
    int i;
    while (fdat_read.good()) {

      fdat_read >> i;
      cout << i << endl;
      Lepton vector_lepton_kaon = FindTwoTracks(Run_Number, i, 0, 5, 0, 0);
      // cout << i << ", " << vector_lepton_kaon.reduced_loss << ", "
      //   << vector_lepton_kaon.reduced_ChiS << ", "
      //   << vector_lepton_kaon.two_track_reduced_loss << ", "
      //   << vector_lepton_kaon.vec_xx_lepton.size() << ", "
      //   << vector_lepton_kaon.first_track_ChiS << ", "
      //   << vector_lepton_kaon.second_track_ChiS << ", " << endl;
      // return;
      fout << i << ", " << vector_lepton_kaon.reduced_loss << ", "
        << vector_lepton_kaon.reduced_ChiS << ", "
        // << vector_lepton_kaon.two_track_reduced_loss << ", "
        << vector_lepton_kaon.vec_xx_lepton.size() << ", " << endl;
        // << vector_lepton_kaon.first_track_ChiS << ", "
        // << vector_lepton_kaon.second_track_ChiS << ", " << endl;
      // else fout << i << ", " << vector_lepton_kaon.reduced_loss << ", "
      //   << vector_lepton_kaon.reduced_ChiS << ", "
      //   << vector_lepton_kaon.two_track_reduced_loss << ", "
      //   << vector_lepton_kaon.vec_xx_lepton.size() << ", "
      //   << vector_lepton_kaon.second_track_ChiS << ", "
      //   << vector_lepton_kaon.first_track_ChiS << ", " << endl;
    }
  }
  else if (flag == 2) {
    char batch_output[50];
    // sprintf(batch_output, "RUN_%d_Batch_Distance_Fit_Large_Distance.txt", Run_Number);
    sprintf(batch_output, "RUN_%d_Batch_read_Large_Distance.txt", Run_Number);
    fout.open(batch_output);
    int i = 0;
    while (i < nentries) {

      Lepton vector_lepton_kaon = MC_Event_Display(Run_Number, i, 0, 0);
      // fout << i << ", " << vector_lepton_kaon.reduced_loss << ", "
      //   << vector_lepton_kaon.reduced_ChiS << ", " << vector_lepton_kaon.vec_xx_lepton.size()<< endl;
      // return;
      if (vector_lepton_kaon.reduced_loss > 3 && !isinf(vector_lepton_kaon.reduced_loss)
          && vector_lepton_kaon.lepton_hit_bars.size() > 12) {
        // fout << i << ", " << vector_lepton_kaon.reduced_loss << ", "
        //   << vector_lepton_kaon.reduced_ChiS << ", " << vector_lepton_kaon.vec_xx_lepton.size()<< endl;
        fout << i << endl;

      }
      ++i;
    }
  }
  else if (flag == 3) {
    char batch_output[50];
    sprintf(batch_output, "RUN_%d_Batch_two_tracks.txt", Run_Number);
    fout.open(batch_output);
    char events[200];
    sprintf(events,"/Users/alng/Documents/Datascience/TREK/offline/alex/simulation/RUN_%d_Batch_read_Large_Distance.txt",
    Run_Number);
    ifstream fdat_read(events, ios::in);
    int i;
    while (fdat_read.good()) {

      fdat_read >> i;
      cout << endl << "Reading event: " << i << endl;
      Lepton vector_lepton_kaon = FindTwoTracks(Run_Number, i, 0, 5, 0, 0);
      // cout <<  i << ", " << vector_lepton_kaon.reduced_loss << ", "
      //   << vector_lepton_kaon.reduced_ChiS << ", "
      //   << vector_lepton_kaon.two_track_reduced_loss << ", "
      //   << vector_lepton_kaon.vec_xx_lepton.size() << ", "
      //   << vector_lepton_kaon.first_track_ChiS << ", "
      //   << vector_lepton_kaon.second_track_ChiS << ", " << endl;
      // return;
      if (vector_lepton_kaon.first_track_ChiS < vector_lepton_kaon.second_track_ChiS) fout << i << ", " << vector_lepton_kaon.reduced_loss << ", "
        << vector_lepton_kaon.reduced_ChiS << ", "
        << vector_lepton_kaon.two_track_reduced_loss << ", "
        << vector_lepton_kaon.vec_xx_lepton.size() << ", "
        << vector_lepton_kaon.first_track_ChiS << ", "
        << vector_lepton_kaon.second_track_ChiS << ", " << endl;
      else fout << i << ", " << vector_lepton_kaon.reduced_loss << ", "
        << vector_lepton_kaon.reduced_ChiS << ", "
        << vector_lepton_kaon.two_track_reduced_loss << ", "
        << vector_lepton_kaon.vec_xx_lepton.size() << ", "
        << vector_lepton_kaon.second_track_ChiS << ", "
        << vector_lepton_kaon.first_track_ChiS << ", " << endl;

    }
  }
  else if (flag == 4) {
    // char batch_output[50];
    // sprintf(batch_output, "RUN_%d_Batch_two_tracks.txt", Run_Number);
    // fout.open(batch_output);
    char events[200];
    sprintf(events,"/Users/alng/Documents/Datascience/TREK/offline/alex/simulation/RUN_%d_Batch_read_Large_Distance.txt",
    Run_Number);
    ifstream fdat_read(events, ios::in);
    int i;
    while (fdat_read.good()) {

      fdat_read >> i;
      cout << endl << "Reading event: " << i << endl;
      Lepton vector_lepton_kaon = FindTwoTracks(Run_Number, i, 1);

    }
  }
  else if (flag == 5) {
    // char batch_output[50];
    // sprintf(batch_output, "RUN_%d_Batch_two_tracks.txt", Run_Number);
    // fout.open(batch_output);

    int i = 0;
    while (i < nentries) {

      cout << endl << "Reading event: " << i << endl;
      Lepton vector_lepton_kaon = FindTwoTracks(Run_Number, i, 1, 0, 0, 0);
      ++i;
    }
  }
  else if (flag == 6) {
    char events[200];
    sprintf(events,"/Users/alng/Documents/Datascience/TREK/offline/alex/simulation/RUN_%d_read_distance_fit.txt",
    Run_Number);
    ifstream fdat_read(events, ios::in);
    int i;
    while (fdat_read.good()) {

      fdat_read >> i;
      cout << i << endl;
      Lepton vector_lepton_kaon = FindTwoTracks(Run_Number, i, 1);
    }
  }
  else if (flag == 7) {
    char batch_output[50];
    sprintf(batch_output, "RUN_%d_y_event_determination.txt", Run_Number);
    fout.open(batch_output);
    char events[200];
    sprintf(events,"/Users/alng/Documents/Datascience/TREK/offline/alex/simulation/RUN_%d_Batch_read_Large_Distance.txt",
    Run_Number);
    ifstream fdat_read(events, ios::in);
    int i;
    int stop = 0;
    while (fdat_read.good() && stop < 100000) {

      fdat_read >> i;
      cout << endl << "Reading event: " << i << endl;
      Lepton vector_lepton_kaon = FindTwoTracks(Run_Number, i, 0, 5, 0, 0);

      fout << i << ", " << vector_lepton_kaon.y_event << "\n";
      cout << i;
      ++stop;
    }
  }
  else if (flag == 8) {
    char batch_output[50];
    sprintf(batch_output, "RUN_%d_Find_Two_Tracks.txt", Run_Number);
    fout.open(batch_output);
    fout << "";
    fout.close();
    if (max_evt == 0) {
      max_evt = nentries;
    }
    int i = min_evt;
    // while (i < nentries) {
    while (i < max_evt) {


      Lepton vector_lepton_kaon = FindTwoTracks(Run_Number, i, 1);
      ++i;
    }
  }


  fout.close();
}
