#include "Many_Tracks_Fit_2.2.C"
#include <sys/stat.h>
#include <cstdarg>

// flag = 0: batch mode without display, dump to Run_#_FindTwoTracks.txt
// flag = 1: display mode, all canvases saved to pdf (<1000 events allowed)
// min_evt default at 0 (start at beginning)
// if max_evt omitted or set to 0, then proceeds until the end of root file
void Batch_Two_Tracks_Fit(int Run_Number, int flag, int min_evt = 0, int max_evt = 0) {

  ofstream fout;
  char batch_output[50];
  sprintf(batch_output, "RUN_%d_Find_Two_Tracks.txt", Run_Number);

  // check if file exists
  struct stat buf;
  bool fout_existed = false;
  if (stat(batch_output, &buf) != -1) {
	fout_existed = true;
  }

  // if fout didn't exist ie is empty, print header
  //fout.open(batch_output, ios::app);
  fout.open(batch_output, ios::trunc);
//  if (false == fout_existed) {
	  fout << "#Run_Number,ievt,reduced_ChiS,ndf,first_track_ChiS/ndf_1,ndf_1,second_track_ChiS/ndf_2,ndf_2,E_positron,delta_energy,delta_length*10,delta_length_xy,angle_between,two track fit,angle_primary_sim,angle_lepton_all,angle_primary,delta_angle_primary,angle_secondary,delta_angle_secondary,k_stop_error,angle_primary_error,targL,no_leptons,no_kaons,k_stop_radius,TWO_TRACK_MIN_CHISQ,TWO_TRACK_MIN_LEPTON_BARS,FIT_TOF1_WEIGHT,K_STOP_CENTROID_THRESH,PATH_TRAVERSIAL_USE_ALL,PATH_TRAVERSIAL_DIJKSTRA_JUMP_RADIUS,PATH_TRAVERSIAL_ALL_PENALTY\n";
//  }
  fout.close();

  char Name_finput[200];
  sprintf(Name_finput,"G4Run%d.root", Run_Number);
  cout << "File opened:  " << Name_finput << endl;


  Int_t nentries = StartChain(fChain, Name_finput);

  if (max_evt == 0) {
    max_evt = nentries;
  }
  if (flag == 0) { // analyze all events, no plots 
    for (int i = min_evt; i < max_evt; ++i) {
		//std::cout << "\r";
		std::cout << std::setw(8) << i << "/" << (max_evt-1) << " (" << std::setw(4) << std::setprecision(3) << float(i)/float(max_evt-1)*100 << "%)";
		std::cout << std::endl;
        FindTwoTracks(Run_Number, i, 1, 0, 0, 1);
    }
	cout << std::endl;
  }
  else if (flag == 1) { //analyze all events, plots to pdf
    char plot_output[50];
    sprintf(plot_output, "RUN_%d_FindTwoTracks.pdf[", Run_Number);
    c2 = new TCanvas("Event_Display.C  --  TARGET & SFT","Event_Display.C  --  TARGET & SFT",0,200,1050,700);
    c2->Print(plot_output);
    delete c2;
    for (int i = min_evt; i < max_evt; ++i) {
      cout << "Reading: " << i << "\n";
        FindTwoTracks(Run_Number, i, 1, 1, 0, 1);
    }
    sprintf(plot_output, "RUN_%d_FindTwoTracks.pdf]", Run_Number);
    c2->Print(plot_output);

  }




}

/*
 * Run Batch_Two_Tracks(n, 0) on multiple runs. The function takes a variable
 * number of run_numbers. The first argument is the number of following run
 * numbers.
 *
 * const int arg_count  number of runs to analyze
 * ...                  one or more integer run numbers to analyze
 */
void Multi_Batch(const int arg_count, ...) {
	va_list args;
	va_start(args, arg_count);

	for (int i = 0; i < arg_count; i++) {
		int run_number = va_arg(args, int);
		Batch_Two_Tracks_Fit(run_number, 0);
	}

	va_end(args);
}
