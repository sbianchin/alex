#include "FindKaons.h"
#include "EventPlot.h"
#include "TweakParameters.h"
#include "/mnt/hdd1/trek/E36/TEST/offline_new/sandbox/ADC_TARGET_Pedestals.h"
#include "/mnt/hdd1/trek/E36/TEST/offline_new/sandbox/TDC_Windows.h"

Int_t StartChain(TChain &fChain, char Name_finput[200]){
    fChain.Reset();

    fChain.Add(Name_finput);
    fChain.SetMakeClass(1);
    fChain.SetBranchAddress("n", &n);  //
    fChain.SetBranchAddress("s_hit", &s_hit); // 
    fChain.SetBranchAddress("x_hit", &x_hit); //
    fChain.SetBranchAddress("tof1charge[5]", tof1charge); //
    fChain.SetBranchAddress("tof1N[5]", tof1N);
    fChain.SetBranchAddress("px[10]", px); //
    fChain.SetBranchAddress("py[10]", py); //
    fChain.SetBranchAddress("pz[10]", pz); //
    fChain.SetBranchAddress("u_x", &u_x); //
    fChain.SetBranchAddress("u_y", &u_y); //
    fChain.SetBranchAddress("u_z", &u_z); //
    fChain.SetBranchAddress("uto_x", &uto_x);  //
    fChain.SetBranchAddress("uto_y", &uto_y); //
    fChain.SetBranchAddress("uto_z", &uto_z); //

    fChain.SetBranchAddress("delta_x", &delta_x);
    fChain.SetBranchAddress("delta_y", &delta_y);
    fChain.SetBranchAddress("delta_z", &delta_z);
    fChain.SetBranchAddress("delta_stop_x", &delta_stop_x);
    fChain.SetBranchAddress("delta_stop_y", &delta_stop_y);
    fChain.SetBranchAddress("delta_stop_z", &delta_stop_z);
    fChain.SetBranchAddress("delta_ux", &delta_ux);
    fChain.SetBranchAddress("delta_uy", &delta_uy);
    fChain.SetBranchAddress("delta_uz", &delta_uz);
    fChain.SetBranchAddress("delta_energy", &delta_energy);
    fChain.SetBranchAddress("delta_length", &delta_length);
    fChain.SetBranchAddress("delta_time", &delta_time);

    fChain.SetBranchAddress("genEvt", &genEvt); //
    fChain.SetBranchAddress("genEne[10]", genEne); //
    fChain.SetBranchAddress("cosTh", &cosTh); //
    fChain.SetBranchAddress("sinTh", &sinTh); //
    fChain.SetBranchAddress("Phi", &Phi); 
    fChain.SetBranchAddress("dv_x", &dv_x); 
    fChain.SetBranchAddress("dv_y", &dv_y); 
    fChain.SetBranchAddress("dv_z", &dv_z); //
    fChain.SetBranchAddress("targetE[5]", targetE);  //
    fChain.SetBranchAddress("targetdE[256]", targetdE);
    fChain.SetBranchAddress("targetdt[256]", targetdt);
    fChain.SetBranchAddress("targL", &targL); //
    fChain.SetBranchAddress("t1pL", &t1pL); //
    fChain.SetBranchAddress("tof1E[5]", tof1E); 
    fChain.SetBranchAddress("tof1wpos[5][3]", tof1wpos);
    fChain.SetBranchAddress("tof1dt", &tof1dt); //
    fChain.SetBranchAddress("secondP[5][100]", secondP); //
    fChain.SetBranchAddress("t1", &t1); //
    fChain.SetBranchAddress("t1px", &t1px); //
    fChain.SetBranchAddress("t1py", &t1py); // 
    fChain.SetBranchAddress("t1pz", &t1pz); //
    fChain.SetBranchAddress("SFTwpos1[5][3]", SFTwpos1); //
    fChain.SetBranchAddress("SFTwpos2[5][3]", SFTwpos2); //
    fChain.SetBranchAddress("SFTpx", &SFTpx); //
    fChain.SetBranchAddress("SFTpy", &SFTpy); //
    fChain.SetBranchAddress("SFTpz", &SFTpz); //

    Int_t nentries = 0;
    nentries = (Int_t)fChain.GetEntries();
    return nentries;
}

Int_t StartChain_Exp(TChain &eChain, char Name_finput[200]){
    eChain.Reset();
  
    eChain.Add(Name_finput);
    eChain.SetMakeClass(1);

    eChain.SetBranchAddress("TDC_Trig",tdc_trigger); // tdc_trigger for TDC_diff calculation
    eChain.SetBranchAddress("ADC_High_TARGET",adc_high_target);
    eChain.SetBranchAddress("ADC_Low_TARGET",adc_low_target);
    eChain.SetBranchAddress("TDC_LE_TARGET",tdc_le_target);
    eChain.SetBranchAddress("TDC_TE_TARGET",tdc_te_target);

    eChain.SetBranchAddress("ADC_High_SFT",adc_high_sft);
    eChain.SetBranchAddress("ADC_Low_SFT",adc_low_sft);
    eChain.SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
    eChain.SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

    eChain.SetBranchAddress("ADC_TOF1U",ADC_tof1U);
    eChain.SetBranchAddress("ADC_TOF1D",ADC_tof1D);
    eChain.SetBranchAddress("TDC_TOF1U",TDC_tof1U);
    eChain.SetBranchAddress("TDC_TOF1D",TDC_tof1D);

    eChain.SetBranchAddress("ADC_TOF2AO",ADC_tof2AO);
    eChain.SetBranchAddress("ADC_TOF2AI",ADC_tof2AI);
    eChain.SetBranchAddress("ADC_TOF2BO",ADC_tof2BO);
    eChain.SetBranchAddress("ADC_TOF2BI",ADC_tof2BI);
    eChain.SetBranchAddress("TDC_TOF2AO",TDC_tof2AO);
    eChain.SetBranchAddress("TDC_TOF2AI",TDC_tof2AI);
    eChain.SetBranchAddress("TDC_TOF2BO",TDC_tof2BO);
    eChain.SetBranchAddress("TDC_TOF2BI",TDC_tof2BI);

    eChain.SetBranchAddress("MWPCADC",MwpcADC);

    eChain.SetBranchAddress("ADC_C2X_R",adc_c2x_r);
    eChain.SetBranchAddress("ADC_C2X_L",adc_c2x_l);
    eChain.SetBranchAddress("ADC_C2Y_R",adc_c2y_r);
    eChain.SetBranchAddress("ADC_C2Y_L",adc_c2y_l);
    eChain.SetBranchAddress("ADC_C3X_R",adc_c3x_r);
    eChain.SetBranchAddress("ADC_C3X_L",adc_c3x_l);
    eChain.SetBranchAddress("ADC_C3Y_R",adc_c3y_r);
    eChain.SetBranchAddress("ADC_C3Y_L",adc_c3y_l);
    eChain.SetBranchAddress("ADC_C4X_R",adc_c4x_r);
    eChain.SetBranchAddress("ADC_C4X_L",adc_c4x_l);
    eChain.SetBranchAddress("ADC_C4Y_R",adc_c4y_r);
    eChain.SetBranchAddress("ADC_C4Y_L",adc_c4y_l);

    eChain.SetBranchAddress("TDC_Ck", tdc_ck);
    eChain.SetBranchAddress("TDC_Cpi", tdc_cpi);

    eChain.SetBranchAddress("VT48_TDC",tdc_vt48);

    eChain.SetBranchAddress("EvFlag", Event_flag);
    eChain.SetBranchAddress("EvTag", Event_tag);

    Int_t nentries = 0;
    nentries = (Int_t)eChain.GetEntries();
    return nentries;
}


double first_fit_residual(vector <double> vec_xx, vector <double> vec_yy, double m, double b) {
    double loss = 0;
    for (unsigned int i = 0; i < vec_xx.size(); ++i) {
      loss += pow(distance_to_line(vec_xx[i], vec_yy[i], m, b),2);
    }
    return loss;


}


/*
 * Find the K-stop by intersecting linear fits for all lepton and all kaon bars.
 * The TOF1 hit point is weighted 3 in the lepton fit.
 *
 * std::vector<Int_t> kaon_hit_bars
 * std::vector<Double_t> vec_xx_kaon    Vector of spacial x components of all
 *                                      kaon bars
 * std::vector<Double_t> vec_yy_kaon    Vector of spacial y components of all
 *                                      kaon bars
 * std::vector<Double_t> vec_xx_lepton  Vector of spacial x components of all
 *                                      lepton bars
 * std::vector<Double_t> vec_yy_lepton  Vector of spacial y components of all
 *                                      lepton bars
 *
 * returns: Coordinate  the coordinate of the K-stop
 */
Coordinate KStopByIntersection(TGraphErrors* kaon_graph, TGraphErrors* lepton_graph) {
    TF1* kaon_fit = kaon_graph->GetFunction("kaon_fit");
    Double_t m_kaon = kaon_fit->GetParameter(1);
    Double_t b_kaon = kaon_fit->GetParameter(0);

    TF1* lepton_fit = lepton_graph->GetFunction("lepton_fit_3");
    Double_t m_lepton = lepton_fit->GetParameter(1);
    Double_t b_lepton = lepton_fit->GetParameter(0);

    vector<double> k_stop = _2lines_intersect(m_lepton, b_lepton, m_kaon, b_kaon);
    Double_t k_stop_x = k_stop[0];
    Double_t k_stop_y = k_stop[1];

	return Coordinate(k_stop_x, k_stop_y);
}


/*
 * Find K-stop by taking the centroid of all kaon bars with certain energy and
 * below a time threshold.
 *
 * std::vector<Int_t> kaon_hit_bars  a vector with kaon bar indices
 * map<int, vector<double>> kaon_time_energy_map
 *                                   a map that maps kaon bar indices to a
 *                                   vector with energy as the the first element
 *                                   and the time as the second element.
 *
 * return: Coordinate  coordinates of the K-stop
 */
Coordinate KStopByCentroid(std::vector<Int_t> kaon_hit_bars, map<int, vector<double>> kaon_time_energy_map) {
	Coordinate k_stop(0,0);
    vector<int> kaon_stop_determination;
    // Find kaon bars with minimum energy and time
    for (unsigned int i = 0; i < kaon_hit_bars.size(); ++i) {
      double kaon_energy = kaon_time_energy_map[kaon_hit_bars[i]][0];
      double kaon_time = kaon_time_energy_map[kaon_hit_bars[i]][1];
      if (Min_E_k < kaon_energy && kaon_time < Min_T_k) {
        kaon_stop_determination.push_back(kaon_hit_bars[i]);
      }
    }
    if (kaon_stop_determination.size() == 0) {
		kaon_stop_determination = kaon_hit_bars;
	}

	for (unsigned int i = 0; i < kaon_stop_determination.size(); ++i) {
		k_stop.x += Xloc[kaon_stop_determination[i]];
		k_stop.y += Yloc[kaon_stop_determination[i]];
	}
	k_stop.x /= kaon_stop_determination.size();
	k_stop.y /= kaon_stop_determination.size();

	return k_stop;
}


Coordinate DetermineKStop(TGraphErrors* kaon_graph, std::vector<Int_t> kaon_hit_bars, map<int, vector<double>> kaon_time_energy_map, TGraphErrors* lepton_graph) {
  Coordinate k_stop(std::nan(""), std::nan(""));

  //Determine x,y-range of kaon
  Double_t min_kaon_x = std::numeric_limits<Double_t>::infinity();
  Double_t min_kaon_y = std::numeric_limits<Double_t>::infinity();
  Double_t max_kaon_x = - std::numeric_limits<Double_t>::infinity();
  Double_t max_kaon_y = - std::numeric_limits<Double_t>::infinity();

  std::vector<Double_t> vec_xx_kaon;
  std::vector<Double_t> vec_yy_kaon;
  for (Int_t bar : kaon_hit_bars) {
     vec_xx_kaon.push_back(Xloc[bar]);
     vec_yy_kaon.push_back(Yloc[bar]);
  }

  for (unsigned int i = 0; i < kaon_hit_bars.size(); ++i) {
    if (min_kaon_x > vec_xx_kaon[i]) {
      min_kaon_x = vec_xx_kaon[i];
    }
    if (min_kaon_y > vec_yy_kaon[i]) {
      min_kaon_y = vec_yy_kaon[i];
    }
    if (max_kaon_x < vec_xx_kaon[i]) {
      max_kaon_x = vec_xx_kaon[i];
    }
    if (max_kaon_y < vec_yy_kaon[i]) {
      max_kaon_y = vec_yy_kaon[i];
    }
  }
  Double_t kaon_x_range = max_kaon_x - min_kaon_x;
  Double_t kaon_y_range = max_kaon_y - min_kaon_y;

  /*
   * If a minimum number of Kaon bars was hit and they cover a certain range,
   * try using the intersection method for the K-stop, else use the centroid.
   * If the resulting point is too far away from any kaon bars or lies between
   * the target and the SFT then also use the centroid.
   */
  if (kaon_hit_bars.size() > 2 && (kaon_x_range > 6.2 || kaon_y_range > 6.2)) {
	  k_stop = KStopByIntersection(kaon_graph, lepton_graph);

	  Double_t dist_from_kaon_track = std::numeric_limits<Double_t>::infinity();
	  for (size_t i = 0; i < vec_xx_kaon.size(); ++i) {
		  Double_t dist = TMath::Sq(vec_xx_kaon[i] - k_stop.x) + TMath::Sq(vec_yy_kaon[i] - k_stop.y);
		  if (dist < dist_from_kaon_track) {
			  dist_from_kaon_track = dist;
		  }
	  }
	  dist_from_kaon_track = TMath::Sqrt(dist_from_kaon_track);

	  const Double_t k_stop_dist_threshold = K_STOP_CENTROID_THRESH /* mm */;
	  const Double_t radius_target = 30 /* mm */;
	  const Double_t radius_SFT = 40 /* mm */;
	  const Double_t radius_k_stop = TMath::Sqrt( TMath::Sq(k_stop.x) + TMath::Sq(k_stop.y));

	  if (dist_from_kaon_track > k_stop_dist_threshold) {
		  k_stop = KStopByCentroid(kaon_hit_bars, kaon_time_energy_map);
	  } else if (radius_k_stop > radius_target && radius_k_stop < radius_SFT) {
		  k_stop = KStopByCentroid(kaon_hit_bars, kaon_time_energy_map);
	  }
  } else {
	  k_stop = KStopByCentroid(kaon_hit_bars, kaon_time_energy_map);
  }

  return k_stop;
}


Lepton FindKaons(Int_t Switch_Output){

  vector <double> vec_xx_lepton;              vec_xx_lepton.clear();
  vector <double> vec_xx_lepton_low_energy;   vec_xx_lepton_low_energy.clear();
  vector <double> vec_yy_lepton_low_energy;   vec_yy_lepton_low_energy.clear();
  vector <double> vec_xx_kaon_low_energy;     vec_xx_kaon_low_energy.clear();
  vector <double> vec_yy_kaon_low_energy;     vec_yy_kaon_low_energy.clear();
  vector <double> vec_ex_lepton_low_energy;   vec_ex_lepton_low_energy.clear();
  vector <double> vec_ey_lepton_low_energy;   vec_ey_lepton_low_energy.clear();
  vector <double> vec_ex_kaon_low_energy;     vec_ex_kaon_low_energy.clear();
  vector <double> vec_ey_kaon_low_energy;     vec_ey_kaon_low_energy.clear();
  vector <double> vec_yy_lepton;              vec_yy_lepton.clear();
  vector <double> vec_ex_lepton;              vec_ex_lepton.clear();
  vector <double> vec_ey_lepton;              vec_ey_lepton.clear();
  vector <double> vec_xx_lepton_test;         vec_xx_lepton_test.clear();
  vector <double> vec_yy_lepton_test;         vec_yy_lepton_test.clear();
  vector <double> vec_xx_lepton_rotate;       vec_xx_lepton_rotate.clear();
  vector <double> vec_yy_lepton_rotate;       vec_yy_lepton_rotate.clear();
  vector <double> vec_ex_kaon;                vec_ex_kaon.clear();
  vector <double> vec_ey_kaon;                vec_ey_kaon.clear();


  vector <double> vec_xx_kaon;                vec_xx_kaon.clear();
  vector <double> vec_yy_kaon;                vec_yy_kaon.clear();

  vector <int> lepton_hit_bars;               lepton_hit_bars.clear();
  vector <int> kaon_hit_bars;                 kaon_hit_bars.clear();
  Double_t Emax[4] = {-1.0};
  Double_t tEmax[4] = {-1.0};
  Double_t Tmin[4] = {-1000.0};
  Double_t Tmax[4] = {-1000.0};
  Int_t idEmax[4] = {-1};
  Int_t idTmin[4] = {-1};
  Int_t idTmax[4] = {-1};
  // Find minimum time
  Tmin[0] = Thit[0];


  //Tmax[0] = Thit[nhits-1];
  idTmin[0] = 0;
  idTmax[0] = nhits - 1;

  for(Int_t ii = 0; ii<nhits; ++ii){
     //For each entry in Ehit,
     //ensure that the minimum T is in 0th position in Tmin
     //ensure that the maximum T is in 0th position in Tmax
      if((Ehit[ii] > MinE)) { // ignore "pedestal"
         // cout<<ii<<" Thit "<<Thit[ii]<<endl;
          if(Tmin[0] > Thit[ii]){
            //if there's a time that's smaller, change it
              Tmin[0] = Thit[ii];
              idTmin[0] = ii;
          }
          if(Thit[ii] > Tmax[0]){
            //add maximum time
              Tmax[0] = Thit[ii];
              idTmax[0] = ii;
          }
        //  cout<<ii<<" Tmin0 "<<Tmin[0]<<" "<<idTmin[0]<<endl;
        //  cout<<ii<<" Tmax0 "<<Tmax[0]<<" "<<idTmax[0]<<endl;
      }
  }
  //1st element in Tmin is 10 after 0th position of Tmin
  Tmin[1] = Tmin[0] + 10.0;

  for(Int_t ii = 0; ii<nhits; ++ii){
    //ensure that the 2nd smallest time is in 1st position in Tmin
    //ensure that the 2nd biggest time is in 1st position in Tmax
      if((Ehit[ii] > MinE) && (ii != idTmin[0])) { // ignore "Tmin[0]"
         // cout<<ii<<" Thit1 "<<Thit[ii]<<endl;
          if(Tmin[1] > Thit[ii]){
              Tmin[1] = Thit[ii];
              idTmin[1] = ii;
          }
          if((ii != idTmax[0])) { // ignore "Tmin[0]"
          if(Thit[ii] > Tmax[1]){
              Tmax[1] = Thit[ii];
              idTmax[1] = ii;
          }
          }
        //  cout<<ii<<" Tmin1 "<<Tmin[1]<<" "<<idTmin[1]<<endl;
        //  cout<<ii<<" Tmax1 "<<Tmax[1]<<" "<<idTmax[1]<<endl;
      }
  }
  // Assume max energy is a kaon hit!


  for(Int_t ii = 0; ii<nhits; ++ii){
      if((Ehit[ii] > MinE)) { // ignore "pedestal"
      //pedestal is already ignored in FindTrack
          if(Ehit[ii] > Emax[0]){
              Emax[0] = Ehit[ii];
              tEmax[0] = Thit[ii];
              idEmax[0] = ii;
          }
      }
  }


  for(Int_t ii = 0; ii<nhits; ++ii){
      if((Ehit[ii] > MinE) && (ii != idEmax[0])) { // ignore "pedestal" and id0
          if(Ehit[ii] > Emax[1]){
              Emax[1] = Ehit[ii];
              tEmax[1] = Thit[ii];
              idEmax[1] = ii;
          }
      }
  }

  //Calculate Time separation between Lepton and Kaons
  bool time_separation = true;
  Double_t minMean = (Tmin[0] + Tmin[1])/2.0;
  Double_t maxMean = (Tmax[0] + Tmax[1])/2.0;
  // Tstep = (maxMean - minMean)/2.0;
  Tstep = (maxMean - minMean)/2.0;
  // if (to_string(Run_Number)[0] < 5 || Run_Number == 524 || Run_Number == 525 |) {
    if (Tstep < 0.5) {
      if (Switch_Output == 1) {
        cout << "Error: Not enough time separation for lepton/kaon discrimination" << endl;
      }
      time_separation = false;
    }
  // }

  bool changeTstep = false;
  if(Tmin[1]  > (Tmin[0] + Tstep)){
      minMean = Tmin[0];
      changeTstep = true;
  }
  if(Tmax[0]  > (Tmax[1] + Tstep)){
      maxMean = Tmax[0];
      changeTstep = true;
  }
  if(changeTstep){
      Tstep = (maxMean - minMean)/2.0;
  }

  // cout << endl << "Tstep: " << Tstep << endl;


  double time_kaon = 0.;
  double time_lepton = 0.;
  int n_leptons_small_energy = 0; // Number of bars that do not make the cut

  // cout << "Time of hit: " ;
  //Discriminate between Lepton and Kaon
  map<int, vector<double>> lepton_time_energy_map;
  map<int, vector<double>> kaon_time_energy_map;
  for(Int_t ii = 0; ii<nhits; ++ii){

      if (time_separation) { //if there is large enough separation
        Double_t minDiff = Thit[ii] - minMean;
        //cout << minDiff;
       // cout<<"Time Diff "<<minDiff<<endl;
        if((minDiff > Tstep) && (minDiff < 70.0)) {
          vector<double> energy_time = {Ehit[ii], Thit[ii]};
          lepton_time_energy_map[hit_bars[ii]] = energy_time;
          if((Ehit[ii] > min_E_valid_hit)) { // ignore "pedestal"

                vec_xx_lepton.push_back(Xhit[ii]);
                vec_yy_lepton.push_back(Yhit[ii]);
                vec_ex_lepton.push_back(TARGET_Errors_X);
                vec_ey_lepton.push_back(TARGET_Errors_Y);
                vec_xx_lepton_rotate.push_back(Xloc[TARGET_Rotated_index[hit_bars[ii]]]);
                vec_yy_lepton_rotate.push_back(Yloc[TARGET_Rotated_index[hit_bars[ii]]]);
                lepton_hit_bars.push_back(hit_bars[ii]);
                time_lepton = time_lepton + Thit[ii];


            }
            else {
              //By tuning the min_E_valid_hit parameter, we don't have to consider
              //bars with low dE, but still be able to view them in our displays
              ++n_leptons_small_energy;
              vec_xx_lepton_low_energy.push_back(Xhit[ii]);
              vec_yy_lepton_low_energy.push_back(Yhit[ii]);
              vec_ex_lepton_low_energy.push_back(TARGET_Errors_X);
              vec_ey_lepton_low_energy.push_back(TARGET_Errors_Y);
            }
          }
          else if((minDiff > -1.0) && (minDiff <= Tstep)) {
              vector<double> energy_time = {Ehit[ii], Thit[ii]};
              kaon_time_energy_map[hit_bars[ii]] = energy_time;
              if((Ehit[ii] > min_E_valid_hit)) { // ignore "pedestal"               
                  vec_xx_kaon.push_back(Xhit[ii]);
                  vec_yy_kaon.push_back(Yhit[ii]);
                  kaon_hit_bars.push_back(hit_bars[ii]);
                  vec_ex_kaon.push_back(TARGET_Errors_X);
                  vec_ey_kaon.push_back(TARGET_Errors_Y);
                  time_kaon = time_kaon + Thit[ii];

              }
              else {
                vec_xx_kaon_low_energy.push_back(Xhit[ii]);
                vec_yy_kaon_low_energy.push_back(Yhit[ii]);
                vec_ex_lepton_low_energy.push_back(TARGET_Errors_X);
                vec_ey_lepton_low_energy.push_back(TARGET_Errors_Y);
              }

        }
      }
      else { //Add all the bars as leptons if there isn't enough separation
      vector<double> energy_time = {Ehit[ii], Thit[ii]};
      lepton_time_energy_map[hit_bars[ii]] = energy_time;
      if((Ehit[ii] > min_E_valid_hit)) {
        vec_xx_lepton.push_back(Xhit[ii]);
        vec_yy_lepton.push_back(Yhit[ii]);
        vec_ex_lepton.push_back(TARGET_Errors_X);
        vec_ey_lepton.push_back(TARGET_Errors_Y);
        vec_xx_lepton_rotate.push_back(Xloc[TARGET_Rotated_index[hit_bars[ii]]]);
        vec_yy_lepton_rotate.push_back(Yloc[TARGET_Rotated_index[hit_bars[ii]]]);
        lepton_hit_bars.push_back(hit_bars[ii]);
        time_lepton = time_lepton + Thit[ii];

      } else {
        ++n_leptons_small_energy;
        vec_xx_lepton_low_energy.push_back(Xhit[ii]);
        vec_yy_lepton_low_energy.push_back(Yhit[ii]);
        vec_ex_lepton_low_energy.push_back(TARGET_Errors_X);
        vec_ey_lepton_low_energy.push_back(TARGET_Errors_Y);
      }
    }
  }

  // //Add TOF1 hits for fit
  // for (int i = 0; i < 5; ++i) {
  //   if (tof1wpos[i][0] > 10) {
  //     vec_xx_lepton.push_back(-tof1wpos[i][1]);
  //     vec_yy_lepton.push_back(tof1wpos[i][0]);
  //     vec_ex_lepton.push_back(TOF1_Errors_X[(tof1N[i]+8)%12][0]);
  //     vec_ey_lepton.push_back(TOF1_Errors_Y[(tof1N[i]+8)%12][0]);
  //
  //   }
  //
  // }


  if(Switch_Output == 1) {
    if (lepton_hit_bars.size() > 0) time_lepton = time_lepton / lepton_hit_bars.size();
    if (kaon_hit_bars.size() > 0) time_kaon = time_kaon / kaon_hit_bars.size();
    cout << endl << "T_Kaon_avg (ns): " << time_kaon;
    cout << endl << "T_Lepton_avg (ns): " << time_lepton;
    cout << endl << "n_leptons_small_energy: " << n_leptons_small_energy;
    cout << "\nDelta_time (ns): " << delta_time << "\n";
    cout << endl << "Lepton Bars Printed in Order (Bar #): Energy, Time" << endl;
    for ( const auto & item : lepton_time_energy_map )
        std::cout << item.first << ": " << (item.second)[0] << ", "
        << (item.second)[1] << endl;

    cout << endl << "Kaon Bars Printed in Order (Bar #): Energy, Time" << endl;
    for ( const auto & item : kaon_time_energy_map )
        std::cout << item.first << ": " << (item.second)[0] << ", "
        << (item.second)[1] << endl;
        vector<pair<double,vector<double>>> lepton_energy_time_vector( lepton_time_energy_map.begin(), lepton_time_energy_map.end() );

        // Sort the vector according to the word count in descending order.
        std::sort( lepton_energy_time_vector.begin(), lepton_energy_time_vector.end(),
                   []( const pair<double,vector<double>> & lhs, const pair<double,vector<double>> & rhs )
                   { return (lhs.second)[1] < (rhs.second)[1]; } );

        // Print out the vector.
        cout << endl << "Lepton Bar (Time Order): Energy, Time" << endl;
        for ( const auto & item : lepton_energy_time_vector )
            std::cout << item.first << ": " << (item.second)[0] << ", "
            << (item.second)[1] << endl;

        vector<pair<double,vector<double>>> kaon_energy_time_vector( kaon_time_energy_map.begin(), kaon_time_energy_map.end() );

        // Sort the vector according to the word count in descending order.
        std::sort( kaon_energy_time_vector.begin(), kaon_energy_time_vector.end(),
                   []( const pair<double,vector<double>> & lhs, const pair<double,vector<double>> & rhs )
                   { return (lhs.second)[1] < (rhs.second)[1]; } );

        // Print out the vector.
        cout << endl << "Kaon Bar (Time Order): Energy, Time" << endl;
        for ( const auto & item : kaon_energy_time_vector )
            std::cout << item.first << ": " << (item.second)[0] << ", "
            << (item.second)[1] << endl;

  }



  Lepton vector_lepton_kaon;
  vector_lepton_kaon.vec_xx_lepton = vec_xx_lepton;
  vector_lepton_kaon.vec_yy_lepton = vec_yy_lepton;
  vector_lepton_kaon.vec_ex_lepton = vec_ex_lepton;
  vector_lepton_kaon.vec_ey_lepton = vec_ey_lepton;
  vector_lepton_kaon.vec_xx_lepton_rotate = vec_xx_lepton_rotate;
  vector_lepton_kaon.vec_yy_lepton_rotate = vec_yy_lepton_rotate;
  vector_lepton_kaon.vec_xx_kaon = vec_xx_kaon;
  vector_lepton_kaon.vec_yy_kaon = vec_yy_kaon;
  vector_lepton_kaon.vec_ex_kaon = vec_ex_kaon;
  vector_lepton_kaon.vec_ey_kaon = vec_ey_kaon;
  vector_lepton_kaon.vec_xx_lepton_low_energy = vec_xx_lepton_low_energy;
  vector_lepton_kaon.vec_yy_lepton_low_energy = vec_yy_lepton_low_energy;
  vector_lepton_kaon.vec_xx_kaon_low_energy = vec_xx_kaon_low_energy;
  vector_lepton_kaon.vec_yy_kaon_low_energy = vec_yy_kaon_low_energy;
  vector_lepton_kaon.vec_ex_lepton_low_energy = vec_ex_lepton_low_energy;
  vector_lepton_kaon.vec_ey_lepton_low_energy = vec_ey_lepton_low_energy;
  vector_lepton_kaon.vec_ex_kaon_low_energy = vec_ex_kaon_low_energy;
  vector_lepton_kaon.vec_ey_kaon_low_energy = vec_ey_kaon_low_energy;

  vector_lepton_kaon.kaon_time_energy_map = kaon_time_energy_map;
  vector_lepton_kaon.lepton_hit_bars = lepton_hit_bars;
  vector_lepton_kaon.kaon_hit_bars = kaon_hit_bars;

  return vector_lepton_kaon;
}


Lepton MC_Event_Display(Int_t Run_Number, Int_t ievt, Int_t Switch_Display, Int_t Switch_Output, Int_t batch){
  //write to file output


  int Switch=1; // Displays hit with no HG, but LG (0 = OFF ; 1 = ON)
  int Rotate=1; // When TOF1 is 12 or 6, or when TOF2 is 6 or 12, rotate by -90 deg to fit a horizontal line (0 = OFF ; 1 = ON)
  int T_limit = 3;

  Double_t Tstep = 0.2;

  gStyle->Clear();
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(11);


  char Name_finput[200] = "";
  sprintf(Name_finput,"G4Run%d.root", Run_Number);


  if (batch == 0) {
	StartChain(fChain, Name_finput);
  }

  fChain.GetEntry(ievt);


  double Xloc_TOF_flipped = -tof1wpos[0][1];
  double Yloc_TOF_flipped = tof1wpos[0][0];
  double Yloc_TOF = tof1wpos[0][1];
  double Xloc_TOF = tof1wpos[0][0];


  AGAIN:
    Tmax = -1;
    Emax = -1;
    iMaxE = -1;
    iMaxT = -1;
    nhits =0;
    nLGhits =0;
    nKhits =0;
    Tdiffmax = 5;
    HGpedOffset = 5;
    LGpedOffset = 10;
    Double_t AngleTOF1 = 0.0;
    Int_t nTOF1s = 0;
    Int_t dTOF1s[5] = {-999};
    Int_t pTOF1s[5] = {-999};
    //for(Int_t id = 0; id<5; id++) cout << tof1N[id] << endl;
    //return;
    for(Int_t id = 0; id<5; id++){
        if(tof1N[id] >= 0){
            dTOF1s[nTOF1s] = id;
            pTOF1s[nTOF1s] = tof1N[id];
            nTOF1s++;
      //  cout<< " tof1N_n "<<id<<" "<<tof1N[id]<<endl;
        }

    }
    Int_t irt; //stores TOF1 Location
    for(Int_t it = 0; it <12; it++){
        for(Int_t id = 0; id<nTOF1s; id++){
            if(it == tof1N[id]){
                irt = tof1N[id];
                // cout<<" TOF1 Location "<<irt<<" "<<TOF_Xloc[irt]<<" "<<TOF_Yloc[irt]<<endl;
            }
        }
    }
    //Check if we need to rotate for fit
    if ((find(begin(tof1N), end(tof1N), 5) != end(tof1N)) ||
      (find(begin(tof1N), end(tof1N), 11) != end(tof1N))) {
        irt = 1;
      }
    // reset hit stuff
    for(Int_t ii = 0; ii<nBars; ++ii){
        Xhit[ii] = -1000.0;
        Yhit[ii] = -1000.0;
        Ehit[ii] = -1000.0;
        Thit[ii] = -1000.0;
        hit_bars[ii] = -1;
    }
    //Find hits
    // cout << "Energy of hit: ";



    for(Int_t ii = 0; ii<nBars; ++ii){
    if((targetdE[ii] > MinE) && (targetdt[ii] > 0.0)) {
        // cout<<"Target dE dt "<<ii<<" "<<targetdE[ii]<<" "<<targetdt[ii]<<" "<<Xloc[ii]<<" "<<Yloc[ii]<<endl;
        // Xhit[nhits] = Xloc[kek2Triumf[ii]];
        // Yhit[nhits] = Yloc[kek2Triumf[ii]];
        // Ehit[nhits] = targetdE[ii];
        // cout << targetdE[ii] << " ";
        // cout << ii << " ";
        // Thit[nhits] = targetdt[ii];
        // hit_bars[nhits] = kek2Triumf[ii];
        // nhits++;
        // bar_energy.push_back(targetdE[ii]);
        // bar_time.push_back(targetdt[ii]);
        Xhit[nhits] = Xloc[ii];
        Yhit[nhits] = Yloc[ii];
        Ehit[nhits] = targetdE[ii];
        Thit[nhits] = targetdt[ii];
        hit_bars[nhits] = ii;
        nhits++;


      }
        }


    //Calculate the xy-length of delta (defined in G4DataRootApr19.h)
    delta_length_xy = distance(delta_x*10, delta_y*10, delta_stop_x*10, delta_stop_y*10 );


    if (Switch_Output == 1) {
      cout << "File opened:  " << Name_finput << endl;
      cout << "Bars hit: ";
      vector<double> bar_energy, bar_time;
      for(Int_t ii = 0; ii<nBars; ++ii){
      if((targetdE[ii] > MinE) && (targetdt[ii] > 0.0)) {
        bar_energy.push_back(targetdE[ii]);
        bar_time.push_back(targetdt[ii]);
        }
      }
      cout << endl;
      sort (bar_energy.begin(), bar_energy.end());
      sort (bar_time.begin(), bar_time.end());

      cout << endl << "Bar energies, from low to high: ";
      for (unsigned int i = 0; i < bar_energy.size(); ++i) cout << bar_energy[i] << " ";
      cout << endl;
      cout << endl << "Bar time, from low to high: ";
      for (unsigned int i = 0; i < bar_time.size(); ++i) cout << bar_time[i] << " ";
      cout << endl;

      cout << "   " << endl;
      cout << "   " << endl;
      cout << "************************************************************************************************************" << endl;

      cout << "File opened:  " << Name_finput << endl;
      char run_string[100];
      sprintf(run_string,"Run %d ; Event %d",Run_Number,ievt);


      char footer[100];
      sprintf(footer,"Event_Display_Simulation.C  --  Run %d ; Event %d",Run_Number,ievt);

      char Version[100] = "Version 5.0";
      cout << "Analyzing event: " << n << " (Jun's numer)/ " << ievt << " (Root number)\n";
      char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!
      Int_t nentries = (Int_t)fChain.GetEntries();
      cout << "Total Number of Events:  " << nentries <<endl;
      //cout << "  " << endl;
      cout << "Event_Display.C -- " << Version << endl;
      cout << "************************************************************************************************************" << endl;
      cout << "  " << endl;
      cout << "  " << endl;
      //printf("Now Target \n");
      // Set counters etc to 0
      double calculated_phi_degrees = atan2( (tof1wpos[0][1] - dv_y), (tof1wpos[0][0] - dv_x)) * 180 / M_PI;
      if (calculated_phi_degrees < 0) {
        cout << "Calculated Phi (Degrees)= " << calculated_phi_degrees + 360<< endl;
        cout << "Calculated Phi (Radians)= " << atan2( (tof1wpos[0][1] - dv_y), (tof1wpos[0][0] - dv_x)) + 2*M_PI << endl;

      }
      else {
        cout << "Calculated Phi (Degrees)= " << calculated_phi_degrees << endl;
        cout << "Calculated Phi (Radians)= " << atan2( (tof1wpos[0][1] - dv_y), (tof1wpos[0][0] - dv_x)) << endl;


      }
      cout << "Tof1E[i]: ";
      for (int i = 0; i != 5; i++) cout << tof1E[i] << " ";
      //for (int i = 0; i != 5; i++) cout << tof1N[i] << " ";
      cout << "\nTof1 Gap[i]: " << endl;
      for (int i = 0; i != 5; i++) cout << tof1N[i] << " ";
      cout << endl;
      cout << "TOF1X[i]: ";
      for (int i = 0; i != 5; i++) cout << tof1wpos[i][0] << " ";
      cout << endl;
      cout << "TOF1Y[i]: ";
      for (int i = 0; i != 5; i++) cout << tof1wpos[i][1] << " ";
      cout << endl;

      //for (int i = 0; i != 5; i++) cout << tof1wpos[i][0] << " ";
      //for (int i = 0; i != 5; i++) cout << tof1wpos[i][1] << " ";
      //for (int i = 0; i != 5; i++) cout << tof1wpos[i][2] << " ";
      cout << endl;
      // for(Int_t ii = 0; ii<nBars; ++ii){
      // if((targetdE[ii] > 0.0)) {
      //    //  cout<<"Target dE dt "<<ii<<" "<<targetdE[ii]<<" "<<targetdt[ii]<<" "<<Xloc[ii]<<" "<<Yloc[ii]<<endl;
      //     cout << kek2Triumf[ii] << " ";
      //
      //   }
      //
      //     }
      // cout << endl;
      // for(Int_t ii = 0; ii<nBars; ++ii){
      // if(targetdt[ii] > 0.0) {
      //    //  cout<<"Target dE dt "<<ii<<" "<<targetdE[ii]<<" "<<targetdt[ii]<<" "<<Xloc[ii]<<" "<<Yloc[ii]<<endl;
      //     cout << kek2Triumf[ii] << " ";
      //
      //   }
      //     }
      double closest_distance_to_delta = 100000; 
      int closest_bar_to_delta = 999;
      
      for (int i = 0; i < 256; ++i) {
        if (closest_distance_to_delta > distance(Xloc[i],Yloc[i],delta_x*10,delta_y*10)) {
          closest_bar_to_delta = i;
          closest_distance_to_delta = distance(Xloc[i],Yloc[i],delta_x*10,delta_y*10);
        }
      }

      double closest_distance_to_delta_stop = 100000; 
      int closest_bar_to_delta_stop = 999;
      
      for (int i = 0; i < 256; ++i) {
        if (closest_distance_to_delta_stop > distance(Xloc[i],Yloc[i],delta_stop_x*10,delta_stop_y*10)) {
          closest_bar_to_delta_stop = i;
          closest_distance_to_delta_stop = distance(Xloc[i],Yloc[i],delta_stop_x*10,delta_stop_y*10);
        }
      }

      cout << "\nRaw delta_x,y,z (cm): " << delta_x << ", " << delta_y << ", " << delta_z << "\n";
      cout << "Raw track direction cosines: " << delta_ux << ", " << delta_uy << ", " << delta_uz << "\n";
      cout << "\nInitial delta_x,y,z (mm): " << delta_x*10 << ", " << delta_y*10 << ", " << delta_z*10 << "\n";
      cout << "Final delta_x,y,z (mm): " << delta_stop_x*10 << ", " << delta_stop_y*10 << ", " << delta_stop_z*10 << "\n";
      cout << "Bar closest to delta production: " << closest_bar_to_delta << "\n";
      cout << "Bar closest to delta stop: " << closest_bar_to_delta_stop << "\n";

      cout << "Delta Energy (MeV): " << delta_energy << "\n";

      cout << "Delta Length--xyz, xy (mm): " << delta_length*10 << ", " << delta_length_xy << "\n";

    }





    //Lepton vector_lepton_kaon = FindKaons(vec_xx_lepton, vec_yy_lepton, vec_xx_lepton_rotate, vec_yy_lepton_rotate, vec_xx_kaon, vec_yy_kaon);
    Lepton vector_lepton_kaon = FindKaons(Switch_Output);

    vector<double> vec_xx_lepton = vector_lepton_kaon.vec_xx_lepton;
    vector<double> vec_yy_lepton = vector_lepton_kaon.vec_yy_lepton;
    vector<double> vec_ex_lepton = vector_lepton_kaon.vec_ex_lepton;
    vector<double> vec_ey_lepton = vector_lepton_kaon.vec_ey_lepton;
    vector<double> vec_xx_lepton_rotate = vector_lepton_kaon.vec_xx_lepton_rotate;
    vector<double> vec_yy_lepton_rotate = vector_lepton_kaon.vec_yy_lepton_rotate;
    vector<double> vec_xx_kaon = vector_lepton_kaon.vec_xx_kaon;
    vector<double> vec_yy_kaon = vector_lepton_kaon.vec_yy_kaon;
    vector<double> vec_ex_kaon = vector_lepton_kaon.vec_ex_kaon;
    vector<double> vec_ey_kaon = vector_lepton_kaon.vec_ey_kaon;
    vector<double> vec_xx_lepton_low_energy = vector_lepton_kaon.vec_xx_lepton_low_energy;
    vector<double> vec_yy_lepton_low_energy = vector_lepton_kaon.vec_yy_lepton_low_energy;
    vector<double> vec_xx_kaon_low_energy = vector_lepton_kaon.vec_xx_kaon_low_energy;
    vector<double> vec_yy_kaon_low_energy = vector_lepton_kaon.vec_yy_kaon_low_energy;
    vector<double> vec_ex_lepton_low_energy = vector_lepton_kaon.vec_ex_lepton_low_energy;
    vector<double> vec_ey_lepton_low_energy = vector_lepton_kaon.vec_ey_lepton_low_energy;
    vector<double> vec_ex_kaon_low_energy = vector_lepton_kaon.vec_ex_kaon_low_energy;
    vector<double> vec_ey_kaon_low_energy = vector_lepton_kaon.vec_ey_kaon_low_energy;
    vector <int> lepton_hit_bars = vector_lepton_kaon.lepton_hit_bars;
    vector <int> kaon_hit_bars = vector_lepton_kaon.kaon_hit_bars;
    vector_lepton_kaon.dv_x = dv_x;
    vector_lepton_kaon.dv_y = dv_y;

    copy(tof1N, tof1N+5, vector_lepton_kaon.tof1N);


    // Graphs_lepton_kaon graphs_lepton_kaon = Prepare_Plot(Run_Number, ievt,
    // // Prepare_Plot(Run_Number, ievt,
    //   vec_xx_lepton, vec_yy_lepton,
    //   vec_ex_lepton, vec_ey_lepton,
    //   vec_xx_lepton_rotate,
    //   vec_yy_lepton_rotate,
    //   vec_xx_kaon,
    //   vec_yy_kaon, vec_ex_kaon,
    //   vec_ey_kaon);
    //
    // TGraphErrors *gr3_Leptons_rotate = graphs_lepton_kaon.gr3_Leptons_rotate;
    //Fill Graph with Leptons

    vector<double> vec_xx_lepton_plot = vec_xx_lepton;
    vector<double> vec_yy_lepton_plot = vec_yy_lepton;
    vector<double> vec_ex_lepton_plot = vec_ex_lepton;
    vector<double> vec_ey_lepton_plot = vec_ey_lepton;

    for (int i = 0; i < FIT_TOF1_WEIGHT; ++i) {
      vec_xx_lepton_plot.push_back(tof1wpos[0][0]);
      vec_yy_lepton_plot.push_back(tof1wpos[0][1]);
      vec_ex_lepton_plot.push_back(TOF1_Errors_X[tof1N[0]][0]);
      vec_ey_lepton_plot.push_back(TOF1_Errors_Y[tof1N[0]][0]);
    }


    TGraphErrors *gr3_Leptons = new TGraphErrors(vec_xx_lepton_plot.size(),&vec_xx_lepton_plot[0],&vec_yy_lepton_plot[0],
      &vec_ex_lepton_plot[0],&vec_ey_lepton_plot[0]);


    gr3_Leptons->SetMarkerStyle(21);
    gr3_Leptons->SetMarkerColor(kRed);
    gr3_Leptons->SetMarkerSize(0.8);
    gr3_Leptons->GetXaxis()->SetLimits(-50.,50.);
    gr3_Leptons->GetYaxis()->SetRangeUser(-50.,50.);

    TF1 *func_lepton_fit_3 = new TF1("lepton_fit_3", "pol1");
    double a_gr_Leptons = std::nan("");
	  double b_gr_Leptons = std::nan("");
    double ChiS = 999.;
    double loss = 999.;
    int ndf = 999;
    if (vec_xx_lepton.size() > 0) {
      // if(tof1N[0]==1 || tof1N[0]==2 || tof1N[0]==7 || tof1N[0]==8){
      //   func_lepton_fit_3->SetParameter(0,0);
      //   func_lepton_fit_3->SetParameter(1,1);
      // }
      // if(tof1N[0]==4 || tof1N[0]==5 || tof1N[0]==10 || tof1N[0]==11){
      //   func_lepton_fit_3->SetParameter(0,0);
      //   func_lepton_fit_3->SetParameter(1,-1);
      // }
      // else{
      //   func_lepton_fit_3->SetParameter(0,0);
      //   func_lepton_fit_3->SetParameter(1,-1);
      // }
      // func_lepton_fit_3->SetParLimits(0,-50,50);
      // func_lepton_fit_3->SetParLimits(1,-50,50);
      gr3_Leptons->Fit("lepton_fit_3","QW");
      gr3_Leptons->Fit("lepton_fit_3","QC");
      TF1 *func_lepton_fit_3_ptr = gr3_Leptons->GetFunction("lepton_fit_3");
	    func_lepton_fit_3_ptr->SetLineColor(kRed);
	    func_lepton_fit_3_ptr->SetLineWidth(2);
      a_gr_Leptons = func_lepton_fit_3_ptr->GetParameter(1);
      b_gr_Leptons = func_lepton_fit_3_ptr->GetParameter(0);
      ChiS = gr3_Leptons->GetFunction("lepton_fit_3")->GetChisquare();
      ndf = gr3_Leptons->GetFunction("lepton_fit_3")->GetNDF();
      loss = first_fit_residual(vec_xx_lepton, vec_yy_lepton, a_gr_Leptons, b_gr_Leptons);


      vector_lepton_kaon.reduced_loss = loss/ndf;
      vector_lepton_kaon.reduced_ChiS = ChiS/ndf;

      vector_lepton_kaon.ndf = ndf;
    }
    else {
	  vector_lepton_kaon.no_leptons = true;
      if (Switch_Output == 1) {
        cout << "No leptons. No fit attempted.";
      }
      if (Switch_Display == 0) {
        return vector_lepton_kaon;
      }
    }

    if (Switch_Output == 1) {
      cout << "Parameters of Single-Track lepton fit: " << a_gr_Leptons << ", " << b_gr_Leptons << endl;
      cout << endl;

      cout << "Single track ChiS: " << ChiS << ", ndf: " << ndf << ", ChiS/ndf: "<< ChiS/ndf << "\n";

      cout << "Sum Perpendicular Squares to Line: " << loss << endl;
      cout << "Sum Perpendicular Squares to Line/ndf: " << loss/ndf << endl;
    }

      //Prepare Kaon Graph
      TGraphErrors *gr_kaon = new TGraphErrors(vec_xx_kaon.size(),&vec_xx_kaon[0],&vec_yy_kaon[0],
                                         &vec_ey_kaon[0],&vec_ey_kaon[0]);
      gr_kaon->SetMarkerStyle(21);
      gr_kaon->SetMarkerColor(kBlue-6);
      gr_kaon->SetMarkerSize(0.8);
      gr_kaon->GetXaxis()->SetLimits(-50.,50.);
      gr_kaon->GetYaxis()->SetRangeUser(-50.,50.);
      TF1 *gr_kaon_fit = new TF1("kaon_fit", "pol1");
      double a_fit_kaon = 0.;
      double b_fit_kaon = 0.;

      //Fit Kaon

      gr_kaon_fit->SetParameter(0,0);
      gr_kaon_fit->SetParameter(1,1);
      gr_kaon_fit->SetParLimits(0,-50,50);
      gr_kaon_fit->SetParLimits(1,-50,50);

	  if (vec_xx_kaon.size() == 0) {
		  vector_lepton_kaon.no_kaons = true;
		  return vector_lepton_kaon;
	  }

      //if(vec_xx_kaon.size()>2){
        gr_kaon->Fit("kaon_fit","Q");
        gr_kaon_fit = gr_kaon->GetFunction("kaon_fit");
        a_fit_kaon = gr_kaon_fit->GetParameter(1);
        b_fit_kaon = gr_kaon_fit->GetParameter(0);
        gr_kaon_fit->SetLineWidth(2);
        gr_kaon_fit->SetLineColor(9);
      //}

  /*
      vector<double> k_stop2 = _2lines_intersect(a_gr_Leptons, b_gr_Leptons, a_fit_kaon, b_fit_kaon);
      k_stop_x = k_stop2[0];
      k_stop_y = k_stop2[1];

	  vector_lepton_kaon.k_stop_x = k_stop_x;
	  vector_lepton_kaon.k_stop_y = k_stop_y;
  */
      irt = 0;


	Coordinate k_stop = DetermineKStop(gr_kaon, kaon_hit_bars, vector_lepton_kaon.kaon_time_energy_map, gr3_Leptons);
    double k_stop_x = k_stop.x;
    double k_stop_y = k_stop.y;
	vector_lepton_kaon.k_stop_x = k_stop_x;
	vector_lepton_kaon.k_stop_y = k_stop_y;

	const std::pair<Int_t, Coordinate> tof1_hit = MaxEnergyTOF1Hit();
    const Double_t tof1_hit_x = tof1_hit.second.x;
    const Double_t tof1_hit_y = tof1_hit.second.y;
	const Double_t lepton_angle = LineToDisplayAngle(a_gr_Leptons, tof1_hit_x, tof1_hit_y);
	vector_lepton_kaon.angle_lepton_all = lepton_angle;

	const Double_t k_stop_error = distance(k_stop_x, k_stop_y, dv_x*10, dv_y*10);

    if (Switch_Display == 1) {
		TArrow *delta_u_arrow = new TArrow(delta_x*10,delta_y*10,delta_stop_x*10,delta_stop_y*10,0.05,"|>");
		TArrow *delta_u_arrow_init = new TArrow(delta_x*10, delta_y*10, delta_x*10+delta_ux*20, delta_y*10+delta_uy*20, 0.05, "|>");

		EventPlot* lepton_event_plot = new EventPlot("Lepton Fit");
		lepton_event_plot->ShowTGraphErrors(gr3_Leptons);
		lepton_event_plot->ShowPoint(delta_x*10, delta_y*10, 30, kFullSquare, "Short MC Delta Origin")->SetMarkerSize(0.8);
		lepton_event_plot->ShowTArrow(delta_u_arrow);
		lepton_event_plot->ShowTArrow(delta_u_arrow_init)->SetLineColor(13);
		lepton_event_plot->ShowValue("#theta_Lepton (deg)", lepton_angle);
		lepton_event_plot->Draw(c2, 1);

		EventPlot* kaon_event_plot = new EventPlot("Kaon Fit");
		kaon_event_plot->ShowTGraphErrors(gr_kaon);
		kaon_event_plot->ShowValue("Run", Run_Number);
		kaon_event_plot->ShowValue("Sim Event", n);
		kaon_event_plot->ShowValue("Event", ievt);
		kaon_event_plot->Draw(c2, 2);

		EventPlot* k_stop_event_plot = new EventPlot("K-Stop Determination");
		k_stop_event_plot->ShowTGraphErrors(gr_kaon);
		k_stop_event_plot->ShowTGraphErrors(gr3_Leptons);
		k_stop_event_plot->ShowPoint(k_stop_x, k_stop_y, 6, kStar, "K-Stop (determined)");
		k_stop_event_plot->ShowPoint(dv_x*10, dv_y*10, 4, kPlus, "K-Stop (simulation)");
		k_stop_event_plot->ShowValue("Delta Energy", delta_energy);
		k_stop_event_plot->ShowValue("K-Stop x", k_stop_x);
		k_stop_event_plot->ShowValue("K-Stop y", k_stop_y);
		k_stop_event_plot->ShowValue("K-Stop Error", k_stop_error);
		k_stop_event_plot->ShowPoint(delta_x*10, delta_y*10, 30, kFullSquare, "Short MC Delta Origin")->SetMarkerSize(0.8);
		k_stop_event_plot->ShowTArrow(delta_u_arrow);
		k_stop_event_plot->ShowTArrow(delta_u_arrow_init)->SetLineColor(13);
		k_stop_event_plot->Draw(c2, 3);
    } else {
      delete func_lepton_fit_3;
      delete gr3_Leptons;
    }

    return vector_lepton_kaon;
}


//Lepton MC_Event_Display_Exp(Int_t Run_Number, Int_t ievt, Int_t Switch_Display, Int_t Switch_Output, Int_t batch){
void MC_Event_Display_Exp(Int_t Run_Number, Int_t ievt, Int_t Switch_Display, Int_t Switch_Output, Int_t batch){
  //write to file output


  int Switch=1; // Displays hit with no HG, but LG (0 = OFF ; 1 = ON)
  int Rotate=1; // When TOF1 is 12 or 6, or when TOF2 is 6 or 12, rotate by -90 deg to fit a horizontal line (0 = OFF ; 1 = ON)
  int T_limit = 3;
  int nentries;

  Double_t Tstep = 0.2;

  gStyle->Clear();
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(11);


  char Name_finput[200] = "";
  //sprintf(Name_finput,"G4Run%d.root", Run_Number);
  sprintf(Name_finput,"/mnt/hdd1/trek/E36/Data/Oct_2015/root/Merged/Run%dMS.root", Run_Number);


  if (batch == 0) {
  StartChain_Exp(eChain, Name_finput);
  }

  eChain.GetEntry(ievt);

  for(int i=0; i<256; i++){
    cout << i << "  " << adc_high_target[i] << endl;
  }



  //double Xloc_TOF_flipped = -tof1wpos[0][1];
  //double Yloc_TOF_flipped = tof1wpos[0][0];
  //double Yloc_TOF = tof1wpos[0][1];
  //double Xloc_TOF = tof1wpos[0][0];


  //AGAIN:
  //  Tmax = -1;
  //  Emax = -1;
  //  iMaxE = -1;
  //  iMaxT = -1;
  //  nhits =0;
  //  nLGhits =0;
  //  nKhits =0;
  //  Tdiffmax = 5;
  //  HGpedOffset = 5;
  //  LGpedOffset = 10;
  //  Double_t AngleTOF1 = 0.0;
  //  Int_t nTOF1s = 0;
  //  Int_t dTOF1s[5] = {-999};
  //  Int_t pTOF1s[5] = {-999};
    //for(Int_t id = 0; id<5; id++) cout << tof1N[id] << endl;
    //return;
    
    //for(Int_t id = 0; id<5; id++){
    //    if(tof1N[id] >= 0){
    //        dTOF1s[nTOF1s] = id;
    //        pTOF1s[nTOF1s] = tof1N[id];
    //        nTOF1s++;
      //  cout<< " tof1N_n "<<id<<" "<<tof1N[id]<<endl;
    //    }

    //}
    Int_t irt; //stores TOF1 Location
    //for(Int_t it = 0; it <12; it++){
    //    for(Int_t id = 0; id<nTOF1s; id++){
    //        if(it == tof1N[id]){
    //            irt = tof1N[id];
    //            // cout<<" TOF1 Location "<<irt<<" "<<TOF_Xloc[irt]<<" "<<TOF_Yloc[irt]<<endl;
    //        }
    //    }
    //}

    //Check if we need to rotate for fit
    //if ((find(begin(tof1N), end(tof1N), 5) != end(tof1N)) ||
    //  (find(begin(tof1N), end(tof1N), 11) != end(tof1N))) {
    //    irt = 1;
    //  }

    // reset hit stuff
    for(Int_t ii = 0; ii<nBars; ++ii){
        Xhit[ii] = -1000.0;
        Yhit[ii] = -1000.0;
        Ehit[ii] = -1000.0;
        Thit[ii] = -1000.0;
        hit_bars[ii] = -1;
    }
    //Find hits
    // cout << "Energy of hit: ";

    Int_t ADC_High_TARGET[256];
    Int_t ADC_Low_TARGET[256];
    Int_t TDC_LE_TARGET[256];
    Int_t HG_TARGET_ADC_Thr[256];
    Int_t LG_TARGET_ADC_Thr[256];
    //TARGET_ADC_Ped_HG
    //TARGET_ADC_Ped_LG
    Int_t TARGET_ADC_Thr_HG_Offset = 20;
    Int_t TARGET_ADC_Thr_LG_Offset = 20;

    for(int i=0; i<256; i++) HG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Ped_HG[i]) + TARGET_ADC_Thr_HG_Offset;
    for(int i=0; i<256; i++)  LG_TARGET_ADC_Thr[i] = round(TARGET_ADC_Ped_LG[i]) + TARGET_ADC_Thr_LG_Offset;

    for (int j_TARGET=0; j_TARGET<256; j_TARGET++){
      ADC_High_TARGET[j_TARGET]=adc_high_target[j_TARGET]-HG_TARGET_ADC_Thr[j_TARGET];
      ADC_Low_TARGET[j_TARGET]=adc_low_target[j_TARGET]-LG_TARGET_ADC_Thr[j_TARGET];
      TDC_LE_TARGET[j_TARGET]=tdc_le_target[j_TARGET][0];
    }

    for(int k=0; k<256; k++){
    //  cout << k << "  " << ADC_High_TARGET[k] << "  " << ADC_Low_TARGET[k] << "  " << TDC_LE_TARGET[k] << endl;
       if((ADC_High_TARGET[k]>0 && tdc_le_target[k][0]>=TARGET_TDC_min[k] && tdc_le_target[k][0]<=TARGET_TDC_max[k]) ||
         (ADC_High_TARGET[k]>0 && tdc_le_target[k][1]>=TARGET_TDC_min[k] && tdc_le_target[k][1]<=TARGET_TDC_max[k]) ||
         (ADC_High_TARGET[k]>0 && tdc_le_target[k][2]>=TARGET_TDC_min[k] && tdc_le_target[k][2]<=TARGET_TDC_max[k]) ||
         (ADC_High_TARGET[k]>0 && tdc_le_target[k][3]>=TARGET_TDC_min[k] && tdc_le_target[k][3]<=TARGET_TDC_max[k]) ||
         (ADC_High_TARGET[k]<0  && ADC_Low_TARGET[k]>0)) cout << "ALBANE !" << "  " << k << endl;
    }

    return;

    for(Int_t ii = 0; ii<nBars; ++ii){
      cout << "targetdE : " << " " << ii << "  " << targetdE[ii] << "  " << MinE << endl; 
     if((targetdE[ii] > MinE) && (targetdt[ii] > 0.0)) {
     //if(adc_high_target){
       // cout<<"Target dE dt "<<ii<<" "<<targetdE[ii]<<" "<<targetdt[ii]<<" "<<Xloc[ii]<<" "<<Yloc[ii]<<endl;
        // Xhit[nhits] = Xloc[kek2Triumf[ii]];
        // Yhit[nhits] = Yloc[kek2Triumf[ii]];
        // Ehit[nhits] = targetdE[ii];
        // cout << targetdE[ii] << " ";
        // cout << ii << " ";
        // Thit[nhits] = targetdt[ii];
        // hit_bars[nhits] = kek2Triumf[ii];
        // nhits++;
        // bar_energy.push_back(targetdE[ii]);
        // bar_time.push_back(targetdt[ii]);
        Xhit[nhits] = Xloc[ii];
        Yhit[nhits] = Yloc[ii];
        Ehit[nhits] = targetdE[ii];
        Thit[nhits] = targetdt[ii];
        hit_bars[nhits] = ii;
        nhits++;


      }
        }



    //Calculate the xy-length of delta (defined in G4DataRootApr19.h)
    delta_length_xy = distance(delta_x*10, delta_y*10, delta_stop_x*10, delta_stop_y*10 );


    if (Switch_Output == 1) {
      cout << "File opened:  " << Name_finput << endl;
      cout << "Bars hit: ";
      vector<double> bar_energy, bar_time;
      for(Int_t ii = 0; ii<nBars; ++ii){
      if((targetdE[ii] > MinE) && (targetdt[ii] > 0.0)) {
        bar_energy.push_back(targetdE[ii]);
        bar_time.push_back(targetdt[ii]);
        }
      }
      cout << endl;
      sort (bar_energy.begin(), bar_energy.end());
      sort (bar_time.begin(), bar_time.end());

      cout << endl << "Bar energies, from low to high: ";
      for (unsigned int i = 0; i < bar_energy.size(); ++i) cout << bar_energy[i] << " ";
      cout << endl;
      cout << endl << "Bar time, from low to high: ";
      for (unsigned int i = 0; i < bar_time.size(); ++i) cout << bar_time[i] << " ";
      cout << endl;

      cout << "   " << endl;
      cout << "   " << endl;
      cout << "************************************************************************************************************" << endl;

      cout << "File opened:  " << Name_finput << endl;
      char run_string[100];
      sprintf(run_string,"Run %d ; Event %d",Run_Number,ievt);


      char footer[100];
      sprintf(footer,"Event_Display_Simulation.C  --  Run %d ; Event %d",Run_Number,ievt);

      char Version[100] = "Version 5.0";
      cout << "Analyzing event: " << n << " (Jun's numer)/ " << ievt << " (Root number)\n";
      char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!
      Int_t nentries = (Int_t)fChain.GetEntries();
      cout << "Total Number of Events:  " << nentries <<endl;
      //cout << "  " << endl;
      cout << "Event_Display.C -- " << Version << endl;
      cout << "************************************************************************************************************" << endl;
      cout << "  " << endl;
      cout << "  " << endl;
      //printf("Now Target \n");
      // Set counters etc to 0
      double calculated_phi_degrees = atan2( (tof1wpos[0][1] - dv_y), (tof1wpos[0][0] - dv_x)) * 180 / M_PI;
      if (calculated_phi_degrees < 0) {
        cout << "Calculated Phi (Degrees)= " << calculated_phi_degrees + 360<< endl;
        cout << "Calculated Phi (Radians)= " << atan2( (tof1wpos[0][1] - dv_y), (tof1wpos[0][0] - dv_x)) + 2*M_PI << endl;

      }
      else {
        cout << "Calculated Phi (Degrees)= " << calculated_phi_degrees << endl;
        cout << "Calculated Phi (Radians)= " << atan2( (tof1wpos[0][1] - dv_y), (tof1wpos[0][0] - dv_x)) << endl;


      }
      cout << "Tof1E[i]: ";
      for (int i = 0; i != 5; i++) cout << tof1E[i] << " ";
      //for (int i = 0; i != 5; i++) cout << tof1N[i] << " ";
      cout << "\nTof1 Gap[i]: " << endl;
      for (int i = 0; i != 5; i++) cout << tof1N[i] << " ";
      cout << endl;
      cout << "TOF1X[i]: ";
      for (int i = 0; i != 5; i++) cout << tof1wpos[i][0] << " ";
      cout << endl;
      cout << "TOF1Y[i]: ";
      for (int i = 0; i != 5; i++) cout << tof1wpos[i][1] << " ";
      cout << endl;

      //for (int i = 0; i != 5; i++) cout << tof1wpos[i][0] << " ";
      //for (int i = 0; i != 5; i++) cout << tof1wpos[i][1] << " ";
      //for (int i = 0; i != 5; i++) cout << tof1wpos[i][2] << " ";
      cout << endl;
      // for(Int_t ii = 0; ii<nBars; ++ii){
      // if((targetdE[ii] > 0.0)) {
      //    //  cout<<"Target dE dt "<<ii<<" "<<targetdE[ii]<<" "<<targetdt[ii]<<" "<<Xloc[ii]<<" "<<Yloc[ii]<<endl;
      //     cout << kek2Triumf[ii] << " ";
      //
      //   }
      //
      //     }
      // cout << endl;
      // for(Int_t ii = 0; ii<nBars; ++ii){
      // if(targetdt[ii] > 0.0) {
      //    //  cout<<"Target dE dt "<<ii<<" "<<targetdE[ii]<<" "<<targetdt[ii]<<" "<<Xloc[ii]<<" "<<Yloc[ii]<<endl;
      //     cout << kek2Triumf[ii] << " ";
      //
      //   }
      //     }
      double closest_distance_to_delta = 100000; 
      int closest_bar_to_delta = 999;
      
      for (int i = 0; i < 256; ++i) {
        if (closest_distance_to_delta > distance(Xloc[i],Yloc[i],delta_x*10,delta_y*10)) {
          closest_bar_to_delta = i;
          closest_distance_to_delta = distance(Xloc[i],Yloc[i],delta_x*10,delta_y*10);
        }
      }

      double closest_distance_to_delta_stop = 100000; 
      int closest_bar_to_delta_stop = 999;
      
      for (int i = 0; i < 256; ++i) {
        if (closest_distance_to_delta_stop > distance(Xloc[i],Yloc[i],delta_stop_x*10,delta_stop_y*10)) {
          closest_bar_to_delta_stop = i;
          closest_distance_to_delta_stop = distance(Xloc[i],Yloc[i],delta_stop_x*10,delta_stop_y*10);
        }
      }

      cout << "\nRaw delta_x,y,z (cm): " << delta_x << ", " << delta_y << ", " << delta_z << "\n";
      cout << "Raw track direction cosines: " << delta_ux << ", " << delta_uy << ", " << delta_uz << "\n";
      cout << "\nInitial delta_x,y,z (mm): " << delta_x*10 << ", " << delta_y*10 << ", " << delta_z*10 << "\n";
      cout << "Final delta_x,y,z (mm): " << delta_stop_x*10 << ", " << delta_stop_y*10 << ", " << delta_stop_z*10 << "\n";
      cout << "Bar closest to delta production: " << closest_bar_to_delta << "\n";
      cout << "Bar closest to delta stop: " << closest_bar_to_delta_stop << "\n";

      cout << "Delta Energy (MeV): " << delta_energy << "\n";

      cout << "Delta Length--xyz, xy (mm): " << delta_length*10 << ", " << delta_length_xy << "\n";

    }





    //Lepton vector_lepton_kaon = FindKaons(vec_xx_lepton, vec_yy_lepton, vec_xx_lepton_rotate, vec_yy_lepton_rotate, vec_xx_kaon, vec_yy_kaon);
    Lepton vector_lepton_kaon = FindKaons(Switch_Output);

    vector<double> vec_xx_lepton = vector_lepton_kaon.vec_xx_lepton;
    vector<double> vec_yy_lepton = vector_lepton_kaon.vec_yy_lepton;
    vector<double> vec_ex_lepton = vector_lepton_kaon.vec_ex_lepton;
    vector<double> vec_ey_lepton = vector_lepton_kaon.vec_ey_lepton;
    vector<double> vec_xx_lepton_rotate = vector_lepton_kaon.vec_xx_lepton_rotate;
    vector<double> vec_yy_lepton_rotate = vector_lepton_kaon.vec_yy_lepton_rotate;
    vector<double> vec_xx_kaon = vector_lepton_kaon.vec_xx_kaon;
    vector<double> vec_yy_kaon = vector_lepton_kaon.vec_yy_kaon;
    vector<double> vec_ex_kaon = vector_lepton_kaon.vec_ex_kaon;
    vector<double> vec_ey_kaon = vector_lepton_kaon.vec_ey_kaon;
    vector<double> vec_xx_lepton_low_energy = vector_lepton_kaon.vec_xx_lepton_low_energy;
    vector<double> vec_yy_lepton_low_energy = vector_lepton_kaon.vec_yy_lepton_low_energy;
    vector<double> vec_xx_kaon_low_energy = vector_lepton_kaon.vec_xx_kaon_low_energy;
    vector<double> vec_yy_kaon_low_energy = vector_lepton_kaon.vec_yy_kaon_low_energy;
    vector<double> vec_ex_lepton_low_energy = vector_lepton_kaon.vec_ex_lepton_low_energy;
    vector<double> vec_ey_lepton_low_energy = vector_lepton_kaon.vec_ey_lepton_low_energy;
    vector<double> vec_ex_kaon_low_energy = vector_lepton_kaon.vec_ex_kaon_low_energy;
    vector<double> vec_ey_kaon_low_energy = vector_lepton_kaon.vec_ey_kaon_low_energy;
    vector <int> lepton_hit_bars = vector_lepton_kaon.lepton_hit_bars;
    vector <int> kaon_hit_bars = vector_lepton_kaon.kaon_hit_bars;
    vector_lepton_kaon.dv_x = dv_x;
    vector_lepton_kaon.dv_y = dv_y;

    copy(tof1N, tof1N+5, vector_lepton_kaon.tof1N);


    // Graphs_lepton_kaon graphs_lepton_kaon = Prepare_Plot(Run_Number, ievt,
    // // Prepare_Plot(Run_Number, ievt,
    //   vec_xx_lepton, vec_yy_lepton,
    //   vec_ex_lepton, vec_ey_lepton,
    //   vec_xx_lepton_rotate,
    //   vec_yy_lepton_rotate,
    //   vec_xx_kaon,
    //   vec_yy_kaon, vec_ex_kaon,
    //   vec_ey_kaon);
    //
    // TGraphErrors *gr3_Leptons_rotate = graphs_lepton_kaon.gr3_Leptons_rotate;
    //Fill Graph with Leptons

    vector<double> vec_xx_lepton_plot = vec_xx_lepton;
    vector<double> vec_yy_lepton_plot = vec_yy_lepton;
    vector<double> vec_ex_lepton_plot = vec_ex_lepton;
    vector<double> vec_ey_lepton_plot = vec_ey_lepton;

    for (int i = 0; i < FIT_TOF1_WEIGHT; ++i) {
      vec_xx_lepton_plot.push_back(tof1wpos[0][0]);
      vec_yy_lepton_plot.push_back(tof1wpos[0][1]);
      vec_ex_lepton_plot.push_back(TOF1_Errors_X[tof1N[0]][0]);
      vec_ey_lepton_plot.push_back(TOF1_Errors_Y[tof1N[0]][0]);
    }


    TGraphErrors *gr3_Leptons = new TGraphErrors(vec_xx_lepton_plot.size(),&vec_xx_lepton_plot[0],&vec_yy_lepton_plot[0],
      &vec_ex_lepton_plot[0],&vec_ey_lepton_plot[0]);


    gr3_Leptons->SetMarkerStyle(21);
    gr3_Leptons->SetMarkerColor(kRed);
    gr3_Leptons->SetMarkerSize(0.8);
    gr3_Leptons->GetXaxis()->SetLimits(-50.,50.);
    gr3_Leptons->GetYaxis()->SetRangeUser(-50.,50.);

    TF1 *func_lepton_fit_3 = new TF1("lepton_fit_3", "pol1");
    double a_gr_Leptons = std::nan("");
    double b_gr_Leptons = std::nan("");
    double ChiS = 999.;
    double loss = 999.;
    int ndf = 999;
    if (vec_xx_lepton.size() > 0) {
      // if(tof1N[0]==1 || tof1N[0]==2 || tof1N[0]==7 || tof1N[0]==8){
      //   func_lepton_fit_3->SetParameter(0,0);
      //   func_lepton_fit_3->SetParameter(1,1);
      // }
      // if(tof1N[0]==4 || tof1N[0]==5 || tof1N[0]==10 || tof1N[0]==11){
      //   func_lepton_fit_3->SetParameter(0,0);
      //   func_lepton_fit_3->SetParameter(1,-1);
      // }
      // else{
      //   func_lepton_fit_3->SetParameter(0,0);
      //   func_lepton_fit_3->SetParameter(1,-1);
      // }
      // func_lepton_fit_3->SetParLimits(0,-50,50);
      // func_lepton_fit_3->SetParLimits(1,-50,50);
      gr3_Leptons->Fit("lepton_fit_3","QW");
      gr3_Leptons->Fit("lepton_fit_3","QC");
      TF1 *func_lepton_fit_3_ptr = gr3_Leptons->GetFunction("lepton_fit_3");
      func_lepton_fit_3_ptr->SetLineColor(kRed);
      func_lepton_fit_3_ptr->SetLineWidth(2);
      a_gr_Leptons = func_lepton_fit_3_ptr->GetParameter(1);
      b_gr_Leptons = func_lepton_fit_3_ptr->GetParameter(0);
      ChiS = gr3_Leptons->GetFunction("lepton_fit_3")->GetChisquare();
      ndf = gr3_Leptons->GetFunction("lepton_fit_3")->GetNDF();
      loss = first_fit_residual(vec_xx_lepton, vec_yy_lepton, a_gr_Leptons, b_gr_Leptons);


      vector_lepton_kaon.reduced_loss = loss/ndf;
      vector_lepton_kaon.reduced_ChiS = ChiS/ndf;

      vector_lepton_kaon.ndf = ndf;
    }
    else {
    vector_lepton_kaon.no_leptons = true;
      if (Switch_Output == 1) {
        cout << "No leptons. No fit attempted.";
      }
      if (Switch_Display == 0) {
        return vector_lepton_kaon;
      }
    }

    if (Switch_Output == 1) {
      cout << "Parameters of Single-Track lepton fit: " << a_gr_Leptons << ", " << b_gr_Leptons << endl;
      cout << endl;

      cout << "Single track ChiS: " << ChiS << ", ndf: " << ndf << ", ChiS/ndf: "<< ChiS/ndf << "\n";

      cout << "Sum Perpendicular Squares to Line: " << loss << endl;
      cout << "Sum Perpendicular Squares to Line/ndf: " << loss/ndf << endl;
    }

      //Prepare Kaon Graph
      TGraphErrors *gr_kaon = new TGraphErrors(vec_xx_kaon.size(),&vec_xx_kaon[0],&vec_yy_kaon[0],
                                         &vec_ey_kaon[0],&vec_ey_kaon[0]);
      gr_kaon->SetMarkerStyle(21);
      gr_kaon->SetMarkerColor(kBlue-6);
      gr_kaon->SetMarkerSize(0.8);
      gr_kaon->GetXaxis()->SetLimits(-50.,50.);
      gr_kaon->GetYaxis()->SetRangeUser(-50.,50.);
      TF1 *gr_kaon_fit = new TF1("kaon_fit", "pol1");
      double a_fit_kaon = 0.;
      double b_fit_kaon = 0.;

      //Fit Kaon

      gr_kaon_fit->SetParameter(0,0);
      gr_kaon_fit->SetParameter(1,1);
      gr_kaon_fit->SetParLimits(0,-50,50);
      gr_kaon_fit->SetParLimits(1,-50,50);

    if (vec_xx_kaon.size() == 0) {
      vector_lepton_kaon.no_kaons = true;
      return vector_lepton_kaon;
    }

      //if(vec_xx_kaon.size()>2){
        gr_kaon->Fit("kaon_fit","Q");
        gr_kaon_fit = gr_kaon->GetFunction("kaon_fit");
        a_fit_kaon = gr_kaon_fit->GetParameter(1);
        b_fit_kaon = gr_kaon_fit->GetParameter(0);
        gr_kaon_fit->SetLineWidth(2);
        gr_kaon_fit->SetLineColor(9);
      //}

  /*
      vector<double> k_stop2 = _2lines_intersect(a_gr_Leptons, b_gr_Leptons, a_fit_kaon, b_fit_kaon);
      k_stop_x = k_stop2[0];
      k_stop_y = k_stop2[1];

    vector_lepton_kaon.k_stop_x = k_stop_x;
    vector_lepton_kaon.k_stop_y = k_stop_y;
  */
      irt = 0;


  Coordinate k_stop = DetermineKStop(gr_kaon, kaon_hit_bars, vector_lepton_kaon.kaon_time_energy_map, gr3_Leptons);
    double k_stop_x = k_stop.x;
    double k_stop_y = k_stop.y;
  vector_lepton_kaon.k_stop_x = k_stop_x;
  vector_lepton_kaon.k_stop_y = k_stop_y;

  const std::pair<Int_t, Coordinate> tof1_hit = MaxEnergyTOF1Hit();
    const Double_t tof1_hit_x = tof1_hit.second.x;
    const Double_t tof1_hit_y = tof1_hit.second.y;
  const Double_t lepton_angle = LineToDisplayAngle(a_gr_Leptons, tof1_hit_x, tof1_hit_y);
  vector_lepton_kaon.angle_lepton_all = lepton_angle;

  const Double_t k_stop_error = distance(k_stop_x, k_stop_y, dv_x*10, dv_y*10);

    if (Switch_Display == 1) {
    TArrow *delta_u_arrow = new TArrow(delta_x*10,delta_y*10,delta_stop_x*10,delta_stop_y*10,0.05,"|>");
    TArrow *delta_u_arrow_init = new TArrow(delta_x*10, delta_y*10, delta_x*10+delta_ux*20, delta_y*10+delta_uy*20, 0.05, "|>");

    EventPlot* lepton_event_plot = new EventPlot("Lepton Fit");
    lepton_event_plot->ShowTGraphErrors(gr3_Leptons);
    lepton_event_plot->ShowPoint(delta_x*10, delta_y*10, 30, kFullSquare, "Short MC Delta Origin")->SetMarkerSize(0.8);
    lepton_event_plot->ShowTArrow(delta_u_arrow);
    lepton_event_plot->ShowTArrow(delta_u_arrow_init)->SetLineColor(13);
    lepton_event_plot->ShowValue("#theta_Lepton (deg)", lepton_angle);
    lepton_event_plot->Draw(c2, 1);

    EventPlot* kaon_event_plot = new EventPlot("Kaon Fit");
    kaon_event_plot->ShowTGraphErrors(gr_kaon);
    kaon_event_plot->ShowValue("Run", Run_Number);
    kaon_event_plot->ShowValue("Sim Event", n);
    kaon_event_plot->ShowValue("Event", ievt);
    kaon_event_plot->Draw(c2, 2);

    EventPlot* k_stop_event_plot = new EventPlot("K-Stop Determination");
    k_stop_event_plot->ShowTGraphErrors(gr_kaon);
    k_stop_event_plot->ShowTGraphErrors(gr3_Leptons);
    k_stop_event_plot->ShowPoint(k_stop_x, k_stop_y, 6, kStar, "K-Stop (determined)");
    k_stop_event_plot->ShowPoint(dv_x*10, dv_y*10, 4, kPlus, "K-Stop (simulation)");
    k_stop_event_plot->ShowValue("Delta Energy", delta_energy);
    k_stop_event_plot->ShowValue("K-Stop x", k_stop_x);
    k_stop_event_plot->ShowValue("K-Stop y", k_stop_y);
    k_stop_event_plot->ShowValue("K-Stop Error", k_stop_error);
    k_stop_event_plot->ShowPoint(delta_x*10, delta_y*10, 30, kFullSquare, "Short MC Delta Origin")->SetMarkerSize(0.8);
    k_stop_event_plot->ShowTArrow(delta_u_arrow);
    k_stop_event_plot->ShowTArrow(delta_u_arrow_init)->SetLineColor(13);
    k_stop_event_plot->Draw(c2, 3);
    } else {
      delete func_lepton_fit_3;
      delete gr3_Leptons;
    }

    return vector_lepton_kaon;
}
