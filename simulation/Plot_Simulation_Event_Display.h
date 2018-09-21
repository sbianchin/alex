

//const int gap_to_fit = 3;
//const int gap_to_fit_rotate = gap_to_fit;
//const int selected_TOF2 = 0;

//Defining data structure holding graph objects
struct Graphs_lepton_kaon {
// TGraphErrors *gr_kaon;
TGraphErrors *gr3_Leptons_rotate;
// TGraphErrors *gr3_Leptons;
double a_fit_TDC_Gap_Fibers=0.;             double b_fit_TDC_Gap_Fibers=0.;
double a_gr_Leptons = 999.99;
double b_gr_Leptons = 999.99;
double a_fit_kaon = 0.;
double b_fit_kaon = 0.;
};

struct FitParams {

double m, b, ChiS, ndf;

};

TCanvas *cTEST;
TCanvas *c2;
TGaxis *A1;
TGaxis *A2;
TH2F *h_EMPTY;
TGraph *gr_kaon;
TH2F *h_TOF1;
TH2F *h_TOF1_closest;
TH2F *h_TOF1_rotate;
TH2F *h_2lines_intersect;
TH2F *h_Target_Center;
TLine *best_fit_rotate;
//TGraph *gr_kaon;
//TGraph *gr3_Leptons_rotate;
//TGraph *gr3_Leptons;

//c2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",0,0, "ProcEvent(Int_t,Int_t,Int_t,TObject*)");

//Plotting parameters
float R_TOF1 = 47.1;
float R_TARGET = 29.0;
float R_SFT_L1 = 40.0;

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


void SetupGaps() {
  Gap1l->SetLineWidth(10);
  Gap1l->SetLineColor(15);

  Gap2l->SetLineWidth(10);
  Gap2l->SetLineColor(15);

  Gap3l->SetLineWidth(10);
  Gap3l->SetLineColor(15);

  Gap4l->SetLineWidth(10);
  Gap4l->SetLineColor(15);

  Gap5l->SetLineWidth(10);
  Gap5l->SetLineColor(15);

  Gap6l->SetLineWidth(10);
  Gap6l->SetLineColor(15);

  Gap7l->SetLineWidth(10);
  Gap7l->SetLineColor(15);

  Gap8l->SetLineWidth(10);
  Gap8l->SetLineColor(15);

  Gap9l->SetLineWidth(10);
  Gap9l->SetLineColor(15);

  Gap10l->SetLineWidth(10);
  Gap10l->SetLineColor(15);

  Gap11l->SetLineWidth(10);
  Gap11l->SetLineColor(15);

  Gap12l->SetLineWidth(10);
  Gap12l->SetLineColor(15);

}

// Graphs_lepton_kaon Prepare_Plot(Int_t Run_Number, Int_t ievt,
// // void Prepare_Plot(Int_t Run_Number, Int_t ievt,
//   vector<double> vec_xx_lepton,
//   vector<double> vec_yy_lepton, vector<double> vec_ex_lepton,
//   vector<double> vec_ey_lepton,
//   vector <double> vec_xx_lepton_rotate,
//   vector <double> vec_yy_lepton_rotate,
//   vector <double> vec_xx_kaon, vector <double> vec_yy_kaon,
//   vector <double> vec_ex_kaon, vector <double> vec_ey_kaon) {
//
//
//
//   //Define the Horizontal and Vertical Lines
//   TLine *hline1 = new TLine(0.38,0.14,0.62,0.14);    hline1->SetLineWidth(2);
//   TLine *hline2 = new TLine(0.30,0.18,0.70,0.18);    hline2->SetLineWidth(2);
//   TLine *hline3 = new TLine(0.26,0.22,0.74,0.22);    hline3->SetLineWidth(2);
//   TLine *hline4 = new TLine(0.22,0.26,0.78,0.26);    hline4->SetLineWidth(2);
//   TLine *hline5 = new TLine(0.18,0.30,0.82,0.30);    hline5->SetLineWidth(2);
//   TLine *hline6 = new TLine(0.18,0.34,0.82,0.34);    hline6->SetLineWidth(2);
//   TLine *hline7 = new TLine(0.14,0.38,0.86,0.38);    hline7->SetLineWidth(2);
//   TLine *hline8 = new TLine(0.14,0.42,0.86,0.42);    hline8->SetLineWidth(2);
//   TLine *hline9 = new TLine(0.14,0.46,0.86,0.46);    hline9->SetLineWidth(2);
//   TLine *hline10 = new TLine(0.14,0.50,0.86,0.50);   hline10->SetLineWidth(2);
//   TLine *hline11 = new TLine(0.14,0.54,0.86,0.54);   hline11->SetLineWidth(2);
//   TLine *hline12 = new TLine(0.14,0.58,0.86,0.58);   hline12->SetLineWidth(2);
//   TLine *hline13 = new TLine(0.14,0.62,0.86,0.62);   hline13->SetLineWidth(2);
//   TLine *hline14 = new TLine(0.18,0.66,0.82,0.66);   hline14->SetLineWidth(2);
//   TLine *hline15 = new TLine(0.18,0.70,0.82,0.70);   hline15->SetLineWidth(2);
//   TLine *hline16 = new TLine(0.22,0.74,0.78,0.74);   hline16->SetLineWidth(2);
//   TLine *hline17 = new TLine(0.26,0.78,0.74,0.78);   hline17->SetLineWidth(2);
//   TLine *hline18 = new TLine(0.30,0.82,0.70,0.82);   hline18->SetLineWidth(2);
//   TLine *hline19 = new TLine(0.38,0.86,0.62,0.86);   hline19->SetLineWidth(2);
//
//   TLine *vline1 = new TLine(0.14,0.38,0.14,0.62);    vline1->SetLineWidth(2);
//   TLine *vline2 = new TLine(0.18,0.30,0.18,0.70);    vline2->SetLineWidth(2);
//   TLine *vline3 = new TLine(0.22,0.26,0.22,0.74);    vline3->SetLineWidth(2);
//   TLine *vline4 = new TLine(0.26,0.22,0.26,0.78);    vline4->SetLineWidth(2);
//   TLine *vline5 = new TLine(0.30,0.18,0.30,0.82);    vline5->SetLineWidth(2);
//   TLine *vline6 = new TLine(0.34,0.18,0.34,0.82);    vline6->SetLineWidth(2);
//   TLine *vline7 = new TLine(0.38,0.14,0.38,0.86);    vline7->SetLineWidth(2);
//   TLine *vline8 = new TLine(0.42,0.14,0.42,0.86);    vline8->SetLineWidth(2);
//   TLine *vline9 = new TLine(0.46,0.14,0.46,0.86);    vline9->SetLineWidth(2);
//   TLine *vline10 = new TLine(0.50,0.14,0.50,0.86);   vline10->SetLineWidth(2);
//   TLine *vline11 = new TLine(0.54,0.14,0.54,0.86);   vline11->SetLineWidth(2);
//   TLine *vline12 = new TLine(0.58,0.14,0.58,0.86);   vline12->SetLineWidth(2);
//   TLine *vline13 = new TLine(0.62,0.14,0.62,0.86);   vline13->SetLineWidth(2);
//   TLine *vline14 = new TLine(0.66,0.18,0.66,0.82);   vline14->SetLineWidth(2);
//   TLine *vline15 = new TLine(0.70,0.18,0.70,0.82);   vline15->SetLineWidth(2);
//   TLine *vline16 = new TLine(0.74,0.22,0.74,0.78);   vline16->SetLineWidth(2);
//   TLine *vline17 = new TLine(0.78,0.26,0.78,0.74);   vline17->SetLineWidth(2);
//   TLine *vline18 = new TLine(0.82,0.30,0.82,0.70);   vline18->SetLineWidth(2);
//   TLine *vline19 = new TLine(0.86,0.38,0.86,0.62);   vline19->SetLineWidth(2);
//
//   ////////////////TARGET TOF
//   TLine *TOF_line9 = new TLine(0.052539,0.38,0.052539,0.62);
//   TLine *TOF_line8 = new TLine(0.172539,0.172154,0.052539,0.38);
//   TLine *TOF_line7 = new TLine(0.172539,0.172154,0.38,0.0521539);
//   TLine *TOF_line6 = new TLine(0.38,0.0521539,0.62,0.0521539);
//   TLine *TOF_line5 = new TLine(0.62,0.0521539,0.8278461,0.172154);
//   TLine *TOF_line4 = new TLine(0.8278461,0.172154,0.9478461,0.38);
//   TLine *TOF_line3 = new TLine(0.9478461,0.38,0.9478461,0.62);
//   TLine *TOF_line2 = new TLine(0.9478461,0.62,0.8278461,0.8278461);
//   TLine *TOF_line1 = new TLine(0.8278461,0.8278461,0.62,0.9478461);
//   TLine *TOF_line12 = new TLine(0.62,0.9478461,0.38,0.9478461);
//   TLine *TOF_line11 = new TLine(0.38,0.9478461,0.172539,0.8278461);
//   TLine *TOF_line10 = new TLine(0.052539,0.62,0.172539,0.8278461);
//   TLine *TOF_line13 = new TLine(0.8348214,0.8393338,0.625558,0.9620716);        TOF_line13->SetLineWidth(20);
//   TLine *TOF_line14 = new TLine(0.9603795,0.627551,0.8394717,0.8369272);        TOF_line14->SetLineWidth(20);
//   TLine *TOF_line15 = new TLine(0.961,0.38,0.961,0.62);                         TOF_line15->SetLineWidth(20);
//   TLine *TOF_line16 = new TLine(0.8394717,0.1630728,0.9580543,0.372449);        TOF_line16->SetLineWidth(20);
//   TLine *TOF_line17 = new TLine(0.6232329,0.040335,0.8324963,0.1606662);        TOF_line17->SetLineWidth(20);
//   TLine *TOF_line18 = new TLine(0.38,0.039,0.62,0.039);                         TOF_line18->SetLineWidth(20);
//   TLine *TOF_line19 = new TLine(0.1651786,0.1606662,0.3721168,0.040335);        TOF_line19->SetLineWidth(20);
//   TLine *TOF_line20 = new TLine(0.1605283,0.1630728,0.03962054,0.372449);       TOF_line20->SetLineWidth(20);
//   TLine *TOF_line21 = new TLine(0.04,0.38,0.04,0.62);                           TOF_line21->SetLineWidth(20);
//   TLine *TOF_line22 = new TLine(0.03962054,0.6251444,0.1605283,0.8345206);      TOF_line22->SetLineWidth(20);
//   TLine *TOF_line23 = new TLine(0.3721168,0.9572584,0.1605283,0.8369272);       TOF_line23->SetLineWidth(20);
//   TLine *TOF_line24 = new TLine(0.62,0.96,0.38,0.96);
//
//   //Define the Text
//   TLatex *tex_Title_ADC_High_TARGET;
//   TLatex *tex_Subtitle_ADC_High_TARGET;
//
//   char ch_ADC_cut_TARGET[100];
//   sprintf(ch_ADC_cut_TARGET,"(ADC offset)");
//
//   char event_string[100];
//   sprintf(event_string,"Run %d ; Event %d",Run_Number,ievt);
//
//   TLatex *tex_event_TARGET;
//   tex_event_TARGET = new TLatex(0.036,0.0,event_string);
//   tex_event_TARGET->SetTextSize(0.038);
//   tex_event_TARGET->SetLineWidth(2);
//
//   //Define the Legend Target
//   TLatex *tex_Legend_TARGET[36];
//   tex_Legend_TARGET[0] = new TLatex(0.36,0.83,"0");     tex_Legend_TARGET[9] = new TLatex(0.09,0.47,"128");
//   tex_Legend_TARGET[1] = new TLatex(0.28,0.79,"6");     tex_Legend_TARGET[10] = new TLatex(0.09,0.43,"146");
//   tex_Legend_TARGET[2] = new TLatex(0.225,0.75,"16");     tex_Legend_TARGET[11] = new TLatex(0.09,0.39,"164");
//   tex_Legend_TARGET[3] = new TLatex(0.185,0.71,"28");     tex_Legend_TARGET[12] = new TLatex(0.13,0.35,"182");
//   tex_Legend_TARGET[4] = new TLatex(0.145,0.67,"42");     tex_Legend_TARGET[13] = new TLatex(0.13,0.31,"198");
//   tex_Legend_TARGET[5] = new TLatex(0.145,0.63,"58");     tex_Legend_TARGET[14] = new TLatex(0.17,0.27,"214");
//   tex_Legend_TARGET[6] = new TLatex(0.1025,0.59,"74");    tex_Legend_TARGET[15] = new TLatex(0.21,0.23,"228");
//   tex_Legend_TARGET[7] = new TLatex(0.1025,0.55,"92");    tex_Legend_TARGET[16] = new TLatex(0.25,0.19,"240");
//   tex_Legend_TARGET[8] = new TLatex(0.09,0.51,"110");     tex_Legend_TARGET[17] = new TLatex(0.33,0.15,"250");
//
//   tex_Legend_TARGET[18] = new TLatex(0.635,0.83,"5");      tex_Legend_TARGET[27] = new TLatex(0.87,0.47,"145");
//   tex_Legend_TARGET[19] = new TLatex(0.71,0.79,"15");     tex_Legend_TARGET[28] = new TLatex(0.87,0.43,"163");
//   tex_Legend_TARGET[20] = new TLatex(0.75,0.75,"27");      tex_Legend_TARGET[29] = new TLatex(0.87,0.39,"181");
//   tex_Legend_TARGET[21] = new TLatex(0.79,0.71,"41");      tex_Legend_TARGET[30] = new TLatex(0.83,0.35,"197");
//   tex_Legend_TARGET[22] = new TLatex(0.83,0.67,"57");      tex_Legend_TARGET[31] = new TLatex(0.83,0.31,"213");
//   tex_Legend_TARGET[23] = new TLatex(0.83,0.63,"73");      tex_Legend_TARGET[32] = new TLatex(0.79,0.27,"227");
//   tex_Legend_TARGET[24] = new TLatex(0.885,0.59,"91");      tex_Legend_TARGET[33] = new TLatex(0.75,0.23,"239");
//   tex_Legend_TARGET[25] = new TLatex(0.87,0.55,"109");     tex_Legend_TARGET[34] = new TLatex(0.71,0.19,"249");
//   tex_Legend_TARGET[26] = new TLatex(0.87,0.51,"127");      tex_Legend_TARGET[35] = new TLatex(0.64,0.15,"255");
//
//   //TARGET variables
//   TMarker *marker_ADC_TARGET[256];
//
//   //Set Gap style
//   SetupGaps();
//
//
//
//   //Prepare Lepton Fit
//   double a_gr_Leptons = 999.99;
//   double b_gr_Leptons = 999.99;
//   double a_fit_TDC_Gap_Fibers=0.;             double b_fit_TDC_Gap_Fibers=0.;
//
//
//
//
//
//
//
//   //Prepare Rotated Leptons Graph
//   TGraphErrors *gr3_Leptons_rotate;
//   if (vec_xx_lepton.size() > 0) {
//     gr3_Leptons_rotate = new TGraphErrors(vec_xx_lepton.size(),&vec_xx_lepton[0],&vec_yy_lepton[0],
//     &vec_ex_lepton[0],&vec_ey_lepton[0]);
//     TF1 *gr_Leptons_fit3_rotate;
//
//
//     gr3_Leptons_rotate->SetMarkerStyle(21);
//     gr3_Leptons_rotate->SetMarkerColor(2);
//     gr3_Leptons_rotate->SetMarkerSize(0.8);
//     gr3_Leptons_rotate->GetXaxis()->SetLimits(-50.,50.);
//     gr3_Leptons_rotate->GetYaxis()->SetRangeUser(-50.,50.);
//
//     TF1 *func_lepton_fit_3_rotate = new TF1("lepton_fit_3_rotate", "pol1");
//     gr3_Leptons_rotate->Fit("lepton_fit_3_rotate","QC");
//
//     TF1 *func_lepton_fit_3_rotate_ptr = gr3_Leptons_rotate->GetFunction("lepton_fit_3_rotate");
//     func_lepton_fit_3_rotate_ptr->SetLineWidth(2);
//     func_lepton_fit_3_rotate_ptr->SetLineColor(2);
//   }
//
//   Graphs_lepton_kaon graphs_lepton_kaon;
//   graphs_lepton_kaon.gr3_Leptons_rotate = gr3_Leptons_rotate;
//   // graphs_lepton_kaon.gr3_Leptons = gr3_Leptons;
//
//   return graphs_lepton_kaon;
//   //return;
//
// }



/*
void Subplot_4(vector<double> vec_xx_lepton ,
  vector<double> vec_yy_lepton ,
  vector <double> vec_xx_lepton_rotate ,
  vector <double> vec_yy_lepton_rotate ,
  vector <double> vec_xx_kaon , vector <double> vec_yy_kaon ,
  */
void Subplot_1(Int_t Run_Number, Int_t ievt, TGraphErrors *gr3_Leptons,
  vector <double> vec_xx_lepton_low_energy,
  vector <double> vec_yy_lepton_low_energy,
  vector <double> vec_ex_lepton_low_energy,
  vector <double> vec_ey_lepton_low_energy)

{

  char Event_Display[200];
  // Add in Entry Number
  // sprintf(Event_Display,"Lepton Fit - Run Number: %d, Event Number: %d, Entry Number: %d",
  // Run_Number, ievt);
  sprintf(Event_Display,"Lepton Fit-Run #: %d, Event #: %d/%d",
  Run_Number, n, ievt);

  TH2F *h_EMPTY = new TH2F(Event_Display, Event_Display, 500, -50, 50, 500, -50, 50);
  TGaxis *A1 = new TGaxis(-50,50,50,50,"pol1",510,"-U");
  TGaxis *A2 = new TGaxis(50,-50,50,50,"pol1",510,"+U");
  TH2F *h_Target_Center = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  h_Target_Center->Fill(0., 0.);
  // h_Target_Center->Fill(Xloc_TOF_flipped, Yloc_TOF_flipped);
  // h_Target_Center->Fill(tof1wpos[0][0], tof1wpos[0][1]);

  h_Target_Center->SetMarkerStyle(5);
  h_Target_Center->SetMarkerColor(1);
  h_Target_Center->SetMarkerSize(2);



  TGraphErrors *gr3_Leptons_low_energy = new TGraphErrors(
    vec_xx_lepton_low_energy.size(),&vec_xx_lepton_low_energy[0],
    &vec_yy_lepton_low_energy[0],
    &vec_ex_lepton_low_energy[0],&vec_ey_lepton_low_energy[0]);



  gr3_Leptons_low_energy->SetMarkerStyle(21);
  gr3_Leptons_low_energy->SetMarkerColor(kBlack);
  gr3_Leptons_low_energy->SetMarkerSize(0.8);
  gr3_Leptons_low_energy->GetXaxis()->SetLimits(-50.,50.);
  gr3_Leptons_low_energy->GetYaxis()->SetRangeUser(-50.,50.);


  TGraph *delta_bar = new TGraph(1);
  delta_bar->SetPoint(0, delta_x*10, delta_y*10);
  delta_bar->SetMarkerStyle(21);
  delta_bar->SetMarkerColor(30);
  delta_bar->SetMarkerSize(0.8);






  TArrow *delta_u_arrow = new TArrow(delta_x*10,delta_y*10,delta_stop_x*10,delta_stop_y*10,0.05,"|>");
  //  ar2->SetAngle(40);
   delta_u_arrow->SetLineWidth(2);



  TH2F *k_stop_point = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  // k_stop_point->Fill(dv_x*10, dv_y*10);
  k_stop_point->SetMarkerStyle(2);
  k_stop_point->SetMarkerColor(4);
  k_stop_point->SetMarkerSize(2.5);



  TH2F *h_TOF1 = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_closest = new TH2F("h_TOF1_closest", "h_TOF1_closest", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_rotate = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50); // ROTATE_CHANGE
  TH2F *h_2lines_intersect = new TH2F("2lines_intersect", "2lines_intersect", 500, -50, 50, 500, -50, 50);



  TH2F *h_Circle = new TH2F("Test", "TARGET and TOF1", 500, -50, 50, 500, -50, 50);
  double Gap[12][3][2] = {{{0}}};
  for(int g=0; g<12; g++){
    Gap[g][0][0] = TOF_Xloc[3*g];
    Gap[g][1][0] = TOF_Xloc[3*g+1];
    Gap[g][2][0] = TOF_Xloc[3*g+2];

    Gap[g][0][1] = TOF_Yloc[3*g];
    Gap[g][1][1] = TOF_Yloc[3*g+1];
    Gap[g][2][1] = TOF_Yloc[3*g+2];
  }

  h_Circle->SetMarkerStyle(5);
  h_Circle->SetMarkerSize(1.2);
  h_Circle->SetLineWidth(2);



  h_TOF1->SetMarkerStyle(20);           h_TOF1_closest->SetMarkerStyle(20);
  h_TOF1->SetMarkerColor(2);            h_TOF1_closest->SetMarkerColor(3);
  h_TOF1->SetMarkerSize(1.5);           h_TOF1_closest->SetMarkerSize(1.5);

  h_TOF1_rotate->SetMarkerStyle(20);
  h_TOF1_rotate->SetMarkerColor(2);
  h_TOF1_rotate->SetMarkerSize(1.5);

  TEllipse *ell = new TEllipse(0, 0, R_TOF1, 0);
  TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0);
  TEllipse *ell_L1 = new TEllipse(0, 0, R_SFT_L1, 0);
  ell->SetFillStyle(0);                 ell_Target->SetFillStyle(0);                 ell_L1->SetFillStyle(0);
  ell->SetLineColor(6);                 ell_Target->SetLineColor(1);                 ell_L1->SetLineColor(4);
  ell->SetLineWidth(1);                 ell_Target->SetLineWidth(1);                 ell_L1->SetLineWidth(1);







  TLine* Gaps[12] = {Gap1l,
  Gap2l, Gap3l, Gap4l, Gap5l,
  Gap6l, Gap7l, Gap8l, Gap9l,Gap10l, Gap11l, Gap12l};

  for (int i = 0; i < 5; ++i ) {
    if (tof1N[i] != -999) Gaps[tof1N[i]-1]->SetLineColor(42);
    // if (tof1N[i] != -999) Gaps[tof1N[i]]->SetLineColor(42);
  }
  //Draw Plot 1

  c2->cd(1)->Range(-50, -50, 50, 50);
  h_EMPTY->Draw();

  Gap1l->Draw();
  Gap2l->Draw("same");
  Gap3l->Draw("same");
  Gap4l->Draw("same");
  Gap5l->Draw("same");
  Gap6l->Draw("same");
  Gap7l->Draw("same");
  Gap8l->Draw("same");
  Gap9l->Draw("same");
  Gap10l->Draw("same");
  Gap11l->Draw("same");
  Gap12l->Draw("Xsame");

  // best_fit_distance->Draw("same");
  //gr3_Leptons_test->Draw("P");
  A1->Draw();
  A2->Draw();

  h_Circle->Draw("same");

  // hFLIP->Draw("same");
  k_stop_point->Draw("same");
  //if(vec_xx_kaon.size()>0) gr_kaon->Draw("P");

  ell->Draw("same");

  ell_Target->Draw("same");

  ell_L1->Draw("same");

  h_Target_Center->Draw("same");

  //h_TOF1->Draw("same");
  //h_TOF1_closest->Draw("same");

  gr3_Leptons->Draw("P");
  gr3_Leptons_low_energy->Draw("P");
  delta_bar->Draw("P");
  delta_u_arrow->Draw("same");


}

void Subplot_2 (TGraphErrors *gr_kaon,
  vector <double> vec_xx_kaon_low_energy,
  vector <double> vec_yy_kaon_low_energy,
  vector <double> vec_ex_kaon_low_energy,
  vector <double> vec_ey_kaon_low_energy,
  Int_t irt)
{


  TH2F *h_EMPTY = new TH2F("Kaon Fit", "Kaon Fit", 500, -50, 50, 500, -50, 50);
  TGaxis *A1 = new TGaxis(-50,50,50,50,"pol1",510,"-U");
  TGaxis *A2 = new TGaxis(50,-50,50,50,"pol1",510,"+U");
  TH2F *h_Target_Center = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  h_Target_Center->Fill(0., 0.);
  h_Target_Center->SetMarkerStyle(5);
  h_Target_Center->SetMarkerColor(1);
  h_Target_Center->SetMarkerSize(1);
  TH2F *h_TOF1 = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_closest = new TH2F("h_TOF1_closest", "h_TOF1_closest", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_rotate = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50); // ROTATE_CHANGE
  TH2F *h_2lines_intersect = new TH2F("2lines_intersect", "2lines_intersect", 500, -50, 50, 500, -50, 50);


  TGraphErrors *gr_kaon_low_energy = new TGraphErrors(
    vec_xx_kaon_low_energy.size(),&vec_xx_kaon_low_energy[0],
    &vec_yy_kaon_low_energy[0],
    &vec_ex_kaon_low_energy[0],&vec_ey_kaon_low_energy[0]);


  gr_kaon_low_energy->SetMarkerStyle(21);
  gr_kaon_low_energy->SetMarkerColor(kOrange);
  gr_kaon_low_energy->SetMarkerSize(0.8);
  gr_kaon_low_energy->GetXaxis()->SetLimits(-50.,50.);
  gr_kaon_low_energy->GetYaxis()->SetRangeUser(-50.,50.);

  TH2F *h_Circle = new TH2F("Test", "TARGET and TOF1", 500, -50, 50, 500, -50, 50);
  double Gap[12][3][2] = {{{0}}};
  for(int g=0; g<12; g++){
    Gap[g][0][0] = TOF_Xloc[3*g];
    Gap[g][1][0] = TOF_Xloc[3*g+1];
    Gap[g][2][0] = TOF_Xloc[3*g+2];

    Gap[g][0][1] = TOF_Yloc[3*g];
    Gap[g][1][1] = TOF_Yloc[3*g+1];
    Gap[g][2][1] = TOF_Yloc[3*g+2];
  }


  h_Circle->SetMarkerStyle(5);
  h_Circle->SetMarkerSize(1.2);
  h_Circle->SetLineWidth(2);



  h_TOF1->SetMarkerStyle(20);           h_TOF1_closest->SetMarkerStyle(20);
  h_TOF1->SetMarkerColor(2);            h_TOF1_closest->SetMarkerColor(3);
  h_TOF1->SetMarkerSize(1.5);           h_TOF1_closest->SetMarkerSize(1.5);

  h_TOF1_rotate->SetMarkerStyle(20);
  h_TOF1_rotate->SetMarkerColor(2);
  h_TOF1_rotate->SetMarkerSize(1.5);

  TEllipse *ell = new TEllipse(0, 0, R_TOF1, 0);
  TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0);
  TEllipse *ell_L1 = new TEllipse(0, 0, R_SFT_L1, 0);
  ell->SetFillStyle(0);                 ell_Target->SetFillStyle(0);                 ell_L1->SetFillStyle(0);
  ell->SetLineColor(6);                 ell_Target->SetLineColor(1);                 ell_L1->SetLineColor(4);
  ell->SetLineWidth(1);                 ell_Target->SetLineWidth(1);                 ell_L1->SetLineWidth(1);

  const int gap_to_fit = irt;
  const int gap_to_fit_rotate = gap_to_fit;
  const int selected_TOF2 = 0;

  c2->cd(2)->Range(-50, -50, 50, 50);

  h_EMPTY->Draw(); //x

  Gap1l->Draw();
  Gap2l->Draw("same");
  Gap3l->Draw("same");
  Gap4l->Draw("same");
  Gap5l->Draw("same");
  Gap6l->Draw("same");
  Gap7l->Draw("same");
  Gap8l->Draw("same");
  Gap9l->Draw("same");
  Gap10l->Draw("same");
  Gap11l->Draw("same");
  Gap12l->Draw("Xsame");

  A1->Draw(); //x
  A2->Draw(); //x

  h_Circle->Draw("same"); //x
  //gr3_Leptons_rotate->Draw("P"); //x
  gr_kaon->Draw("P"); //x
  gr_kaon_low_energy->Draw("P");
  //best_fit_rotate->Draw("same");

  ell->Draw("same");
  ell_Target->Draw("same");

  ell_L1->Draw("same");

  h_Target_Center->Draw("same"); //x
}

void Subplot_3(TGraphErrors *gr3_Leptons,
  TGraphErrors *gr_kaon,
  double dv_x, double dv_y,
  double k_stop_x, double k_stop_y,
  vector <double> vec_xx_lepton_low_energy,
  vector <double> vec_yy_lepton_low_energy,
  vector <double> vec_xx_kaon_low_energy,
  vector <double> vec_yy_kaon_low_energy,
  vector <double> vec_ex_lepton_low_energy,
  vector <double> vec_ey_lepton_low_energy,
  vector <double> vec_ex_kaon_low_energy,
  vector <double> vec_ey_kaon_low_energy){
  
  c2->cd(3);
  char plot_title[100];
  sprintf(plot_title, "Final Fit: (Delta Energy: %0.2f)", delta_energy);
  TH2F *h_EMPTY = new TH2F(plot_title, plot_title, 500, -50, 50, 500, -50, 50);
  TGaxis *A1 = new TGaxis(-50,50,50,50,"pol1",510,"-U");
  TGaxis *A2 = new TGaxis(50,-50,50,50,"pol1",510,"+U");
  TH2F *h_Target_Center = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  h_Target_Center->Fill(0., 0.);
  // h_Target_Center->Fill(tof1wpos[0][0], tof1wpos[0][1]);
  h_Target_Center->SetMarkerStyle(5);
  h_Target_Center->SetMarkerColor(1);
  h_Target_Center->SetMarkerSize(2);
  TH2F *h_TOF1 = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_closest = new TH2F("h_TOF1_closest", "h_TOF1_closest", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_rotate = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50); // ROTATE_CHANGE
  TH2F *h_2lines_intersect = new TH2F("2lines_intersect", "2lines_intersect", 500, -50, 50, 500, -50, 50);

  TH2F *k_stop_point = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  k_stop_point->Fill(dv_x*10, dv_y*10);
  k_stop_point->SetMarkerStyle(2);
  k_stop_point->SetMarkerColor(4);
  k_stop_point->SetMarkerSize(4.5);

  TH2F *k_stop_experimental = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  k_stop_experimental->Fill(k_stop_x,k_stop_y);
  k_stop_experimental->SetMarkerStyle(3);
  k_stop_experimental->SetMarkerColor(6);
  k_stop_experimental->SetMarkerSize(2.5);


  TGraphErrors *gr3_Leptons_low_energy = new TGraphErrors(
    vec_xx_lepton_low_energy.size(),&vec_xx_lepton_low_energy[0],
    &vec_yy_lepton_low_energy[0],
    &vec_ex_lepton_low_energy[0],&vec_ey_lepton_low_energy[0]);


  TGraph *delta_bar = new TGraph(1);
  delta_bar->SetPoint(0, delta_x*10, delta_y*10);
  delta_bar->SetMarkerStyle(21);
  delta_bar->SetMarkerColor(30);
  delta_bar->SetMarkerSize(0.8);
  TArrow *delta_u_arrow = new TArrow(delta_x*10,delta_y*10,delta_stop_x*10,delta_stop_y*10,0.05,"|>");
  //  ar2->SetAngle(40);
   delta_u_arrow->SetLineWidth(2);

  gr3_Leptons_low_energy->SetMarkerStyle(21);
  gr3_Leptons_low_energy->SetMarkerColor(kBlack);
  gr3_Leptons_low_energy->SetMarkerSize(0.8);
  gr3_Leptons_low_energy->GetXaxis()->SetLimits(-50.,50.);
  gr3_Leptons_low_energy->GetYaxis()->SetRangeUser(-50.,50.);

  TGraphErrors *gr_kaon_low_energy = new TGraphErrors(
    vec_xx_kaon_low_energy.size(),&vec_xx_kaon_low_energy[0],
    &vec_yy_kaon_low_energy[0],
    &vec_ex_kaon_low_energy[0],&vec_ey_kaon_low_energy[0]);


  gr_kaon_low_energy->SetMarkerStyle(21);
  gr_kaon_low_energy->SetMarkerColor(kOrange);
  gr_kaon_low_energy->SetMarkerSize(0.8);
  gr_kaon_low_energy->GetXaxis()->SetLimits(-50.,50.);
  gr_kaon_low_energy->GetYaxis()->SetRangeUser(-50.,50.);

  TH2F *h_Circle = new TH2F("Test", "TARGET and TOF1", 500, -50, 50, 500, -50, 50);
  double Gap[12][3][2] = {{{0}}};


  h_Circle->SetMarkerStyle(5);
  h_Circle->SetMarkerSize(1.2);
  h_Circle->SetLineWidth(2);



  h_TOF1->SetMarkerStyle(20);           h_TOF1_closest->SetMarkerStyle(20);
  h_TOF1->SetMarkerColor(2);            h_TOF1_closest->SetMarkerColor(3);
  h_TOF1->SetMarkerSize(1.5);           h_TOF1_closest->SetMarkerSize(1.5);

  h_TOF1_rotate->SetMarkerStyle(20);
  h_TOF1_rotate->SetMarkerColor(2);
  h_TOF1_rotate->SetMarkerSize(1.5);

  TEllipse *ell = new TEllipse(0, 0, R_TOF1, 0);
  TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0);
  TEllipse *ell_L1 = new TEllipse(0, 0, R_SFT_L1, 0);
  ell->SetFillStyle(0);                 ell_Target->SetFillStyle(0);                 ell_L1->SetFillStyle(0);
  ell->SetLineColor(6);                 ell_Target->SetLineColor(1);                 ell_L1->SetLineColor(4);
  ell->SetLineWidth(1);                 ell_Target->SetLineWidth(1);                 ell_L1->SetLineWidth(1);






  //h_2lines_intersect->Fill(vec_test[0],vec_test[1]);
  h_EMPTY->Draw();

  Gap1l->Draw();
  Gap2l->Draw("same");
  Gap3l->Draw("same");
  Gap4l->Draw("same");
  Gap5l->Draw("same");
  Gap6l->Draw("same");
  Gap7l->Draw("same");
  Gap8l->Draw("same");
  Gap9l->Draw("same");
  Gap10l->Draw("same");
  Gap11l->Draw("same");
  Gap12l->Draw("Xsame");

  //gr_Leptons->Draw("P");
  gr_kaon->Draw("P");
  A1->Draw();
  A2->Draw();
  h_Circle->Draw("same");

  ell->Draw("same");
  ell_Target->Draw("same");

  ell_L1->Draw("same");
  h_Target_Center->Draw("same");
  h_TOF1->Draw("same");
  h_TOF1_closest->Draw("same");
  //h_int_TDC_Gap_Fibers->Draw("same");
  //h_int_TDC_Gap_Fibers_SFT->Draw("same");

  //x_guide->Draw();
  //y_guide->Draw();

  //tex_Angle_guide->Draw("same");
  //tex_ChiS->Draw("same");
  //h_int_TDC_TARGET->Draw("same");

  //h_Centroid->Draw("same");
  //h_2lines_intersect->Draw("same");
  //x_sft->Draw("same");
  //y_sft->Draw("same");
  //x_tof1->Draw("same");
  //y_tof1->Draw("same");
  //tex_Label_Centroid_x->Draw("same");
  //tex_Label_Centroid_y->Draw("same");

  //k_stop_line->Draw("same");
  k_stop_point->Draw("same");
  k_stop_experimental->Draw("same");
  gr3_Leptons->Draw("P");
  gr3_Leptons_low_energy->Draw("P");
  gr_kaon_low_energy->Draw("P");
  delta_bar->Draw("P");
  delta_u_arrow->Draw("same");


 }

  vector<vector<double>> fit_rotate(vector<int> lepton_bars) {
  vector<double> lepton_x_coords, lepton_y_coords;
  vector<double> lepton_ex, lepton_ey;
  //rotate CCW 90 degrees
  for (unsigned int bar = 0; bar < lepton_bars.size(); ++bar) {
      lepton_x_coords.push_back(-Yloc[lepton_bars[bar]]);
      lepton_y_coords.push_back(Xloc[lepton_bars[bar]]);
      lepton_ex.push_back(TARGET_Errors_X);
      lepton_ey.push_back(TARGET_Errors_Y);

    }
  TGraphErrors* lepton_track_rotate = new TGraphErrors(lepton_x_coords.size(),&lepton_x_coords[0],&lepton_y_coords[0]
                                          , &lepton_ex[0], &lepton_ey[0]);
  TF1 *lepton_track_rotate_fit = new TF1("lepton_track_rotate", "pol1");

  lepton_track_rotate->Fit("lepton_track_rotate","QW");
  lepton_track_rotate->Fit("lepton_track_rotate","QC");
  //
  lepton_track_rotate_fit = lepton_track_rotate->GetFunction("lepton_track_rotate");

  double m  = lepton_track_rotate_fit->GetParameter(1);
  double b = lepton_track_rotate_fit->GetParameter(0);

  //Return 2 points
  return {{-50*m+b,50}, {50*m+b,-50}};

}


void Subplot_4(TCanvas *canvas, vector<int> lepton_bars, int Subplot,
  const char *plot_title, Color_t color, int gap) //Plot tracks of Leptons
{


  canvas->cd(Subplot);
  // char *p =;
  // sprintf(plot_title, title);
  TH2F *h_EMPTY = new TH2F(plot_title, plot_title, 500, -50, 50, 500, -50, 50);
  TGaxis *A1 = new TGaxis(-50,50,50,50,"pol1",510,"-U");
  TGaxis *A2 = new TGaxis(50,-50,50,50,"pol1",510,"+U");
  TH2F *h_Target_Center = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  h_Target_Center->Fill(0., 0.);

  //h_Target_Center->Fill(tof1wpos[0][0], tof1wpos[0][1]);
  h_Target_Center->SetMarkerStyle(5);
  h_Target_Center->SetMarkerColor(1);
  h_Target_Center->SetMarkerSize(2);
  TH2F *h_TOF1 = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_closest = new TH2F("h_TOF1_closest", "h_TOF1_closest", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_rotate = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50); // ROTATE_CHANGE
  TH2F *h_2lines_intersect = new TH2F("2lines_intersect", "2lines_intersect", 500, -50, 50, 500, -50, 50);
  // TLine *best_fit_rotate = new TLine(50*a_fit_TDC_Gap_Fibers - b_fit_TDC_Gap_Fibers,-50,-50*a_fit_TDC_Gap_Fibers - b_fit_TDC_Gap_Fibers,50); //ROTATE_CHANGE

  vector<double> vec_X_TOF1, vec_Y_TOF1;
  vector<double> vec_EX_TOF1, vec_EY_TOF1;

  vec_X_TOF1.push_back(tof1wpos[0][0]);
  vec_Y_TOF1.push_back(tof1wpos[0][1]);
  vec_EX_TOF1.push_back(TOF1_Errors_X[gap][0]);
  vec_EY_TOF1.push_back(TOF1_Errors_Y[gap][0]);

  TGraphErrors *gr_TOF1 = new TGraphErrors(vec_X_TOF1.size(),&vec_X_TOF1[0],&vec_Y_TOF1[0],&vec_EX_TOF1[0],&vec_EY_TOF1[0]);
  gr_TOF1->SetMarkerStyle(21);
  gr_TOF1->SetMarkerColor(2);
  gr_TOF1->SetMarkerSize(0.8);


  vector<double> lepton_x_coords, lepton_y_coords;
  vector<double> lepton_ex, lepton_ey;
  for (unsigned int bar = 0; bar < lepton_bars.size(); ++bar) {
      lepton_x_coords.push_back(Xloc[lepton_bars[bar]]);
      lepton_y_coords.push_back(Yloc[lepton_bars[bar]]);
      lepton_ex.push_back(TARGET_Errors_X);
      lepton_ey.push_back(TARGET_Errors_Y);

    }
    // //Add TOF1 hits for fit
    // for (int i = 0; i < 5; ++i) {
    //   if (tof1wpos[i][0] != -999) {
    //     double ssr_1 = abs(m*(-tof1wpos[i][1]) + b - tof1wpos[i][0]); //first SSR
    //     double ssr_2 = abs(m2*(-tof1wpos[i][1]) + b2 - tof1wpos[i][0]);
    //     if (ssr_1 <= ssr_2) {
    //       lepton_x_coords.push_back(-tof1wpos[i][1]);
    //       lepton_y_coords.push_back(tof1wpos[i][0]);
    //       lepton_ex.push_back(TOF1_Errors_X[(tof1N[i]+8)%12][0]);
    //       lepton_ey.push_back(TOF1_Errors_Y[(tof1N[i]+8)%12][0]);
    //     }
    //
    //
    //   }
    //
    // }

  TLine *best_fit_rotate;
  TGraphErrors *lepton_track = new TGraphErrors(lepton_x_coords.size(),&lepton_x_coords[0],&lepton_y_coords[0]
                                            , &lepton_ex[0], &lepton_ey[0]);

  Int_t fitStatus = 999;
  if (lepton_bars.size() > 0) {
  //    lepton_track = new TGraphErrors(lepton_x_coords.size(),&lepton_x_coords[0],&lepton_y_coords[0]
  //                                            , &lepton_ex[0], &lepton_ey[0]);
    TF1 *lepton_track_fit = new TF1("lepton_track", "pol1");

    lepton_track->SetMarkerStyle(21);
    lepton_track->SetMarkerColor(color);
    lepton_track->SetMarkerSize(0.8);
    lepton_track->GetXaxis()->SetLimits(-50.,50.);
    lepton_track->GetYaxis()->SetRangeUser(-50.,50.);
    fitStatus = lepton_track->Fit("lepton_track","QW");
    if (fitStatus == 0) { //if fit works well
      lepton_track->Fit("lepton_track","QC");
      //
      TF1 *lepton_track_fit_ptr = lepton_track->GetFunction("lepton_track");
      lepton_track_fit_ptr->SetLineWidth(2);
      lepton_track_fit_ptr->SetLineColor(color);

    }
    else { //if fit fails (rotation necessary)
      vector<vector<double>> after_rotation = fit_rotate(lepton_bars);
      double x1 = after_rotation[0][0]; double y1 = after_rotation[0][1];
      double x2 = after_rotation[1][0]; double y2 = after_rotation[1][1];
      best_fit_rotate = new TLine(x1,y1,x2,y2);
      best_fit_rotate->SetLineWidth(2); //ROTATE_CHANGE
      best_fit_rotate->SetLineColor(color); // ROTATE_CHANGE
      TF1 *lepton_track_fit_ptr = lepton_track->GetFunction("lepton_track");
      lepton_track_fit_ptr->SetLineColor(0);
    }

  }

  TH2F *h_Circle = new TH2F("Test", "TARGET and TOF1", 500, -50, 50, 500, -50, 50);
  double Gap[12][3][2] = {{{0}}};
  /*
  for(int g=0; g<12; g++){
    Gap[g][0][0] = TOF_Xloc[3*g];
    Gap[g][1][0] = TOF_Xloc[3*g+1];
    Gap[g][2][0] = TOF_Xloc[3*g+2];

    Gap[g][0][1] = TOF_Yloc[3*g];
    Gap[g][1][1] = TOF_Yloc[3*g+1];
    Gap[g][2][1] = TOF_Yloc[3*g+2];
  }



  for(int i=0; i<12; i++)
  {
    h_Circle->Fill(Gap[i][0][0], Gap[i][0][1]);
    h_Circle->Fill(Gap[i][1][0], Gap[i][1][1]);
    h_Circle->Fill(Gap[i][2][0], Gap[i][2][1]);
  }
  */

  //h_Circle->SetMarkerStyle(5);
  //h_Circle->SetMarkerSize(1.2);
  //h_Circle->SetLineWidth(2);



  h_TOF1->SetMarkerStyle(20);           h_TOF1_closest->SetMarkerStyle(20);
  h_TOF1->SetMarkerColor(2);            h_TOF1_closest->SetMarkerColor(3);
  h_TOF1->SetMarkerSize(1.5);           h_TOF1_closest->SetMarkerSize(1.5);

  h_TOF1_rotate->SetMarkerStyle(20);
  h_TOF1_rotate->SetMarkerColor(2);
  h_TOF1_rotate->SetMarkerSize(1.5);

  TEllipse *ell = new TEllipse(0, 0, R_TOF1, 0);
  TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0);
  TEllipse *ell_L1 = new TEllipse(0, 0, R_SFT_L1, 0);

  ell->SetFillStyle(0);                 ell_Target->SetFillStyle(0);                 ell_L1->SetFillStyle(0);
  ell->SetLineColor(6);                 ell_Target->SetLineColor(1);                 ell_L1->SetLineColor(4);
  ell->SetLineWidth(1);                 ell_Target->SetLineWidth(1);                 ell_L1->SetLineWidth(1);

  const int gap_to_fit = 0;
  const int gap_to_fit_rotate = gap_to_fit;
  const int selected_TOF2 = 0;

  canvas->cd(Subplot)->Range(-50, -50, 50, 50);

  h_EMPTY->Draw(); //x

  Gap1l->Draw();
  Gap2l->Draw("same");
  Gap3l->Draw("same");
  Gap4l->Draw("same");
  Gap5l->Draw("same");
  Gap6l->Draw("same");
  Gap7l->Draw("same");
  Gap8l->Draw("same");
  Gap9l->Draw("same");
  Gap10l->Draw("same");
  Gap11l->Draw("same");
  Gap12l->Draw("Xsame");

  A1->Draw(); //x
  A2->Draw(); //x

  //h_Circle->Draw("same"); //x
  if(lepton_bars.size()>0) lepton_track->Draw("P"); //x
  if (fitStatus != 0) best_fit_rotate->Draw("same");
  //best_fit_rotate->Draw("same");

  ell->Draw("same");
  ell_Target->Draw("same");

  ell_L1->Draw("same");

  h_Target_Center->Draw("same");   
  if(gap != 999) gr_TOF1->Draw("P");
}




void Subplot_4_test(TCanvas *canvas, vector<int> lepton_bars, int Subplot,
  const char *plot_title, Color_t color){ //Plot tracks of Leptons


  canvas->cd(Subplot);
  // char *p =;
  // sprintf(plot_title, title);
  TH2F *h_EMPTY = new TH2F(plot_title, plot_title, 500, -50, 50, 500, -50, 50);
  TGaxis *A1 = new TGaxis(-50,50,50,50,"pol1",510,"-U");
  TGaxis *A2 = new TGaxis(50,-50,50,50,"pol1",510,"+U");
  TH2F *h_Target_Center = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  h_Target_Center->Fill(0., 0.);
  h_Target_Center->SetMarkerStyle(5);
  h_Target_Center->SetMarkerColor(1);
  h_Target_Center->SetMarkerSize(1);
  TH2F *h_TOF1 = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_closest = new TH2F("h_TOF1_closest", "h_TOF1_closest", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_rotate = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50); // ROTATE_CHANGE
  TH2F *h_2lines_intersect = new TH2F("2lines_intersect", "2lines_intersect", 500, -50, 50, 500, -50, 50);
  // TLine *best_fit_rotate = new TLine(50*a_fit_TDC_Gap_Fibers - b_fit_TDC_Gap_Fibers,-50,-50*a_fit_TDC_Gap_Fibers - b_fit_TDC_Gap_Fibers,50); //ROTATE_CHANGE

  vector<double> lepton_x_coords, lepton_y_coords;
  vector<double> lepton_ex, lepton_ey;
  for (vector<int>::const_iterator i = lepton_bars.begin();
    i != lepton_bars.end(); ++i) {
      int bar = *i;
      lepton_x_coords.push_back(Xloc[bar]);
      lepton_y_coords.push_back(Yloc[bar]);
      lepton_ex.push_back(TARGET_Errors_X);
      lepton_ey.push_back(TARGET_Errors_Y);

    }
  // //Add TOF1 hits for fit
  // for (int i = 0; i < 5; ++i) {
  //   if (tof1wpos[i][0] != -999) {
  //     double ssr_1 = abs(m*(-tof1wpos[i][1]) + b - tof1wpos[i][0]); //first SSR
  //     double ssr_2 = abs(m2*(-tof1wpos[i][1]) + b2 - tof1wpos[i][0]);
  //     if (ssr_1 < ssr_2) {
  //       lepton_x_coords.push_back(-tof1wpos[i][1]);
  //       lepton_y_coords.push_back(tof1wpos[i][0]);
  //       lepton_ex.push_back(TOF1_Errors_X[(tof1N[i]+8)%12][0]);
  //       lepton_ey.push_back(TOF1_Errors_Y[(tof1N[i]+8)%12][0]);
  //     }
  //
  //
  //   }
  //
  // }
  TGraphErrors *lepton_track = new TGraphErrors(lepton_x_coords.size(),&lepton_x_coords[0],&lepton_y_coords[0]
                                            , &lepton_ex[0], &lepton_ey[0]);

  if (lepton_bars.size() > 0) {
  //    lepton_track = new TGraphErrors(lepton_x_coords.size(),&lepton_x_coords[0],&lepton_y_coords[0]
  //                                            , &lepton_ex[0], &lepton_ey[0]);
    // TF1 *lepton_track_fit = new TF1("lepton_track", "pol1");

    lepton_track->SetMarkerStyle(21);
    lepton_track->SetMarkerColor(color);
    lepton_track->SetMarkerSize(0.8);
    lepton_track->GetXaxis()->SetLimits(-50.,50.);
    lepton_track->GetYaxis()->SetRangeUser(-50.,50.);
    // TF1 *lepton_track_fit = new TF1("track1", "pol1");
    // lepton_track->Fit("track1","QW");
    // lepton_track->Fit("track1","QC");
    // lepton_track_fit->SetLineWidth(2);
    // lepton_track_fit->SetLineColor(color);

  }

  TH2F *h_Circle = new TH2F("Test", "TARGET and TOF1", 500, -50, 50, 500, -50, 50);
  double Gap[12][3][2] = {{{0}}};
  /*
  for(int g=0; g<12; g++){
    Gap[g][0][0] = TOF_Xloc[3*g];
    Gap[g][1][0] = TOF_Xloc[3*g+1];
    Gap[g][2][0] = TOF_Xloc[3*g+2];

    Gap[g][0][1] = TOF_Yloc[3*g];
    Gap[g][1][1] = TOF_Yloc[3*g+1];
    Gap[g][2][1] = TOF_Yloc[3*g+2];
  }



  for(int i=0; i<12; i++)
  {
    h_Circle->Fill(Gap[i][0][0], Gap[i][0][1]);
    h_Circle->Fill(Gap[i][1][0], Gap[i][1][1]);
    h_Circle->Fill(Gap[i][2][0], Gap[i][2][1]);
  }
  */

  h_Circle->SetMarkerStyle(5);
  h_Circle->SetMarkerSize(1.2);
  h_Circle->SetLineWidth(2);



  h_TOF1->SetMarkerStyle(20);           h_TOF1_closest->SetMarkerStyle(20);
  h_TOF1->SetMarkerColor(2);            h_TOF1_closest->SetMarkerColor(3);
  h_TOF1->SetMarkerSize(1.5);           h_TOF1_closest->SetMarkerSize(1.5);

  h_TOF1_rotate->SetMarkerStyle(20);
  h_TOF1_rotate->SetMarkerColor(2);
  h_TOF1_rotate->SetMarkerSize(1.5);

  TEllipse *ell = new TEllipse(0, 0, R_TOF1, 0);
  TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0);
  TEllipse *ell_L1 = new TEllipse(0, 0, R_SFT_L1, 0);
  ell->SetFillStyle(0);                 ell_Target->SetFillStyle(0);                 ell_L1->SetFillStyle(0);
  ell->SetLineColor(6);                 ell_Target->SetLineColor(1);                 ell_L1->SetLineColor(4);
  ell->SetLineWidth(1);                 ell_Target->SetLineWidth(1);                 ell_L1->SetLineWidth(1);

  const int gap_to_fit = 0;
  const int gap_to_fit_rotate = gap_to_fit;
  const int selected_TOF2 = 0;

  canvas->cd(Subplot)->Range(-50, -50, 50, 50);

  h_EMPTY->Draw(); //x

  Gap1l->Draw();
  Gap2l->Draw("same");
  Gap3l->Draw("same");
  Gap4l->Draw("same");
  Gap5l->Draw("same");
  Gap6l->Draw("same");
  Gap7l->Draw("same");
  Gap8l->Draw("same");
  Gap9l->Draw("same");
  Gap10l->Draw("same");
  Gap11l->Draw("same");
  Gap12l->Draw("Xsame");

  A1->Draw(); //x
  A2->Draw(); //x

  h_Circle->Draw("same"); //x
  //gr3_Leptons_rotate->Draw("P"); //x
  if(lepton_bars.size()>0) lepton_track->Draw("P"); //x
  //best_fit_rotate->Draw("same");

  ell->Draw("same");
  ell_Target->Draw("same");

  ell_L1->Draw("same");

  h_Target_Center->Draw("same"); //x
}

void Subplot_6(vector<int> lepton_bars_l1, vector<int> lepton_bars_l2,
  const char *plot_title, int Subplot, Color_t color_1, Color_t color_2,
  Double_t dv_x,
  Double_t dv_y, double k_stop_x,  double k_stop_y, int gap) //Plot tracks of Leptons
{
  c2->cd(Subplot);

  TH2F *h_EMPTY = new TH2F(plot_title, plot_title, 500, -50, 50, 500, -50, 50);
  TGaxis *A1 = new TGaxis(-50,50,50,50,"pol1",510,"-U");
  TGaxis *A2 = new TGaxis(50,-50,50,50,"pol1",510,"+U");
  TH2F *h_Target_Center = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  h_Target_Center->Fill(0., 0.);
  //h_Target_Center->Fill(tof1wpos[0][0], tof1wpos[0][1]);
  h_Target_Center->SetMarkerStyle(5);
  h_Target_Center->SetMarkerColor(1);
  h_Target_Center->SetMarkerSize(2);
  TH2F *h_TOF1 = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_closest = new TH2F("h_TOF1_closest", "h_TOF1_closest", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_rotate = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50); // ROTATE_CHANGE
  TH2F *h_2lines_intersect = new TH2F("2lines_intersect", "2lines_intersect", 500, -50, 50, 500, -50, 50);
  //TLine *best_fit_rotate = new TLine(50*a_fit_TDC_Gap_Fibers - b_fit_TDC_Gap_Fibers,-50,-50*a_fit_TDC_Gap_Fibers - b_fit_TDC_Gap_Fibers,50); //ROTATE_CHANGE

  vector<double> vec_X_TOF1, vec_Y_TOF1;
  vector<double> vec_EX_TOF1, vec_EY_TOF1;

  vec_X_TOF1.push_back(tof1wpos[0][0]);
  vec_Y_TOF1.push_back(tof1wpos[0][1]);
  vec_EX_TOF1.push_back(TOF1_Errors_X[gap][0]);
  vec_EY_TOF1.push_back(TOF1_Errors_Y[gap][0]);

  TGraphErrors *gr_TOF1 = new TGraphErrors(vec_X_TOF1.size(),&vec_X_TOF1[0],&vec_Y_TOF1[0],&vec_EX_TOF1[0],&vec_EY_TOF1[0]);
  gr_TOF1->SetMarkerStyle(21);
  gr_TOF1->SetMarkerColor(2);
  gr_TOF1->SetMarkerSize(0.8);

  //Fill 2 tracks with x and y coordinates of hits
  vector<double> lepton_x_coords_l1, lepton_y_coords_l1;
  vector<double> lepton_ex_l1, lepton_ey_l1;
  for (vector<int>::const_iterator i = lepton_bars_l1.begin();
    i != lepton_bars_l1.end(); ++i) {
      int bar = *i;
      lepton_x_coords_l1.push_back(Xloc[bar]);
      lepton_y_coords_l1.push_back(Yloc[bar]);
      lepton_ex_l1.push_back(TARGET_Errors_X);
      lepton_ey_l1.push_back(TARGET_Errors_Y);
    }

  vector<double> lepton_x_coords_l2, lepton_y_coords_l2;
  vector<double> lepton_ex_l2, lepton_ey_l2;
  for (vector<int>::const_iterator i = lepton_bars_l2.begin();
    i != lepton_bars_l2.end(); ++i) {
      int bar = *i;
      lepton_x_coords_l2.push_back(Xloc[bar]);
      lepton_y_coords_l2.push_back(Yloc[bar]);
      lepton_ex_l2.push_back(TARGET_Errors_X);
      lepton_ey_l2.push_back(TARGET_Errors_Y);
    }

  // //Add TOF1 hits for fit
  // for (int i = 0; i < 5; ++i) {
  //   if (tof1wpos[i][0] != -999) {
  //     double ssr_1 = abs(m1*(-tof1wpos[i][1]) + b1 - tof1wpos[i][0]); //first SSR
  //     double ssr_2 = abs(m2*(-tof1wpos[i][1]) + b2 - tof1wpos[i][0]);
  //     if (ssr_1 < ssr_2) {
  //       lepton_x_coords_l1.push_back(-tof1wpos[i][1]);
  //       lepton_y_coords_l1.push_back(tof1wpos[i][0]);
  //       lepton_ex_l1.push_back(TOF1_Errors_X[(tof1N[i]+8)%12][0]);
  //       lepton_ey_l1.push_back(TOF1_Errors_Y[(tof1N[i]+8)%12][0]);
  //     }
  //     else {
  //       lepton_x_coords_l2.push_back(-tof1wpos[i][1]);
  //       lepton_y_coords_l2.push_back(tof1wpos[i][0]);
  //       lepton_ex_l2.push_back(TOF1_Errors_X[(tof1N[i]+8)%12][0]);
  //       lepton_ey_l2.push_back(TOF1_Errors_Y[(tof1N[i]+8)%12][0]);
  //     }
  //
  //
  //   }
  //
  // }
  //Plot of 2 tracks
  TGraphErrors *lepton_track_1, *lepton_track_2;
  Int_t fitStatus1 = 0;
  Int_t fitStatus2 = 0;
  TLine *best_fit_rotate1, *best_fit_rotate2;
  if (lepton_bars_l1.size() > 0) {

    lepton_track_1 = new TGraphErrors(lepton_x_coords_l1.size(),&lepton_x_coords_l1[0],
                &lepton_y_coords_l1[0], &lepton_ex_l1[0],
                &lepton_ey_l1[0] );
    TF1 *lepton_track_fit_1 = new TF1("lepton_track_1", "pol1");
    lepton_track_1->SetMarkerStyle(21);
    lepton_track_1->SetMarkerColor(color_1);
    lepton_track_1->SetMarkerSize(0.8);
    lepton_track_1->GetXaxis()->SetLimits(-50.,50.);
    lepton_track_1->GetYaxis()->SetRangeUser(-50.,50.);
    fitStatus1 = lepton_track_1->Fit("lepton_track_1","QW");
    if (fitStatus1 == 0) { //if fit works well
      lepton_track_1->Fit("lepton_track_1","QC");
      //
      TF1 *lepton_track_fit1_ptr = lepton_track_1->GetFunction("lepton_track_1");
      lepton_track_fit1_ptr->SetLineWidth(2);
      lepton_track_fit1_ptr->SetLineColor(color_1);

    }
    else { //if fit fails (rotation necessary)
      vector<vector<double>> after_rotation1 = fit_rotate(lepton_bars_l1);
      double x1 = after_rotation1[0][0]; double y1 = after_rotation1[0][1];
      double x2 = after_rotation1[1][0]; double y2 = after_rotation1[1][1];
      best_fit_rotate1 = new TLine(x1,y1,x2,y2);
      best_fit_rotate1->SetLineWidth(2); //ROTATE_CHANGE
      best_fit_rotate1->SetLineColor(color_1); // ROTATE_CHANGE
      TF1 *lepton_track_fit1_ptr = lepton_track_1->GetFunction("lepton_track_1");
      lepton_track_fit1_ptr->SetLineColor(0);
    }
  }

  if (lepton_bars_l2.size() > 0) {
    lepton_track_2 = new TGraphErrors(lepton_x_coords_l2.size(),&lepton_x_coords_l2[0],
          &lepton_y_coords_l2[0], &lepton_ex_l2[0],
          &lepton_ey_l2[0]);
    TF1 *lepton_track_fit_2 = new TF1("lepton_track_2", "pol1");

    lepton_track_2->SetMarkerStyle(21);
    lepton_track_2->SetMarkerColor(color_2);
    lepton_track_2->SetMarkerSize(0.8);
    lepton_track_2->GetXaxis()->SetLimits(-50.,50.);
    lepton_track_2->GetYaxis()->SetRangeUser(-50.,50.);
    fitStatus2 = lepton_track_2->Fit("lepton_track_2","QW");
    if (fitStatus2 == 0) { //if fit works well
      lepton_track_2->Fit("lepton_track_2","QC");
      //
      TF1 *lepton_track_fit2_ptr = lepton_track_2->GetFunction("lepton_track_2");
      lepton_track_fit2_ptr->SetLineWidth(2);
      lepton_track_fit2_ptr->SetLineColor(color_2);

    }
    else { //if fit fails (rotation necessary)
      vector<vector<double>> after_rotation2 = fit_rotate(lepton_bars_l2);
      double x1 = after_rotation2[0][0]; double y1 = after_rotation2[0][1];
      double x2 = after_rotation2[1][0]; double y2 = after_rotation2[1][1];
      best_fit_rotate2 = new TLine(x1,y1,x2,y2);
      best_fit_rotate2->SetLineWidth(2); //ROTATE_CHANGE
      best_fit_rotate2->SetLineColor(color_2); // ROTATE_CHANGE
      TF1 *lepton_track_fit2_ptr = lepton_track_2->GetFunction("lepton_track_2");
      lepton_track_fit2_ptr->SetLineColor(0);
    }

  }

  TH2F *h_Circle = new TH2F("Test", "TARGET and TOF1", 500, -50, 50, 500, -50, 50);
  double Gap[12][3][2] = {{{0}}};
  /*
  for(int g=0; g<12; g++){
    Gap[g][0][0] = TOF_Xloc[3*g];
    Gap[g][1][0] = TOF_Xloc[3*g+1];
    Gap[g][2][0] = TOF_Xloc[3*g+2];

    Gap[g][0][1] = TOF_Yloc[3*g];
    Gap[g][1][1] = TOF_Yloc[3*g+1];
    Gap[g][2][1] = TOF_Yloc[3*g+2];
  }



  for(int i=0; i<12; i++)
  {
    h_Circle->Fill(Gap[i][0][0], Gap[i][0][1]);
    h_Circle->Fill(Gap[i][1][0], Gap[i][1][1]);
    h_Circle->Fill(Gap[i][2][0], Gap[i][2][1]);
  }
  */

  h_Circle->SetMarkerStyle(5);
  h_Circle->SetMarkerSize(1.2);
  h_Circle->SetLineWidth(2);



  h_TOF1->SetMarkerStyle(20);           h_TOF1_closest->SetMarkerStyle(20);
  h_TOF1->SetMarkerColor(2);            h_TOF1_closest->SetMarkerColor(3);
  h_TOF1->SetMarkerSize(1.5);           h_TOF1_closest->SetMarkerSize(1.5);

  h_TOF1_rotate->SetMarkerStyle(20);
  h_TOF1_rotate->SetMarkerColor(2);
  h_TOF1_rotate->SetMarkerSize(1.5);

  TEllipse *ell = new TEllipse(0, 0, R_TOF1, 0);
  TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0);
  TEllipse *ell_L1 = new TEllipse(0, 0, R_SFT_L1, 0);
  ell->SetFillStyle(0);                 ell_Target->SetFillStyle(0);                 ell_L1->SetFillStyle(0);
  ell->SetLineColor(6);                 ell_Target->SetLineColor(1);                 ell_L1->SetLineColor(4);
  ell->SetLineWidth(1);                 ell_Target->SetLineWidth(1);                 ell_L1->SetLineWidth(1);

  const int gap_to_fit = 0;
  const int gap_to_fit_rotate = gap_to_fit;
  const int selected_TOF2 = 0;

  c2->cd(Subplot)->Range(-50, -50, 50, 50);


  h_EMPTY->Draw(); //x

  Gap1l->Draw();
  Gap2l->Draw("same");
  Gap3l->Draw("same");
  Gap4l->Draw("same");
  Gap5l->Draw("same");
  Gap6l->Draw("same");
  Gap7l->Draw("same");
  Gap8l->Draw("same");
  Gap9l->Draw("same");
  Gap10l->Draw("same");
  Gap11l->Draw("same");
  Gap12l->Draw("Xsame");

  A1->Draw(); //x
  A2->Draw(); //x

  h_Circle->Draw("same"); //x
  //gr3_Leptons_rotate->Draw("P"); //x
  if(lepton_bars_l1.size()>0) lepton_track_1->Draw("P"); //x
  if(lepton_bars_l2.size()>0) lepton_track_2->Draw("P"); //x
  if (fitStatus1 != 0) best_fit_rotate1->Draw("same");
  if (fitStatus2 != 0) best_fit_rotate2->Draw("same");
  //best_fit_rotate->Draw("same");

  ell->Draw("same");
  ell_Target->Draw("same");

  ell_L1->Draw("same");

  h_Target_Center->Draw("same"); 
  if(gap != 999) gr_TOF1->Draw("P");






  TH2F *k_stop_experimental = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  k_stop_experimental->Fill(k_stop_x,k_stop_y);
  k_stop_experimental->SetMarkerStyle(3);
  k_stop_experimental->SetMarkerColor(6);
  k_stop_experimental->SetMarkerSize(2.5);
  k_stop_experimental->Draw("same");

  TH2F *k_stop_point = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  k_stop_point->Fill(dv_x*10,dv_y*10);
  k_stop_point->SetMarkerStyle(2);
  k_stop_point->SetMarkerColor(4);
  k_stop_point->SetMarkerSize(2.5);
  k_stop_point->Draw("same");
}

void Subplot_6_test(TCanvas *canvas, vector<int> lepton_bars_l1, vector<int> lepton_bars_l2,
  const char *plot_title, int Subplot, Color_t color_1, Color_t color_2,
  Double_t dv_x,
  Double_t dv_y, double k_stop_x,  double k_stop_y) //Plot tracks of Leptons
{
  canvas->cd(Subplot);
  TH2F *h_EMPTY = new TH2F(plot_title, plot_title, 500, -50, 50, 500, -50, 50);
  TGaxis *A1 = new TGaxis(-50,50,50,50,"pol1",510,"-U");
  TGaxis *A2 = new TGaxis(50,-50,50,50,"pol1",510,"+U");
  TH2F *h_Target_Center = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  h_Target_Center->Fill(0., 0.);
  h_Target_Center->SetMarkerStyle(5);
  h_Target_Center->SetMarkerColor(1);
  h_Target_Center->SetMarkerSize(1);
  TH2F *h_TOF1 = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_closest = new TH2F("h_TOF1_closest", "h_TOF1_closest", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_rotate = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50); // ROTATE_CHANGE
  TH2F *h_2lines_intersect = new TH2F("2lines_intersect", "2lines_intersect", 500, -50, 50, 500, -50, 50);
  //TLine *best_fit_rotate = new TLine(50*a_fit_TDC_Gap_Fibers - b_fit_TDC_Gap_Fibers,-50,-50*a_fit_TDC_Gap_Fibers - b_fit_TDC_Gap_Fibers,50); //ROTATE_CHANGE

  //Fill 2 tracks with x and y coordinates of hits
  vector<double> lepton_x_coords_l1, lepton_y_coords_l1;
  vector<double> lepton_ex_l1, lepton_ey_l1;
  for (vector<int>::const_iterator i = lepton_bars_l1.begin();
    i != lepton_bars_l1.end(); ++i) {
      int bar = *i;
      lepton_x_coords_l1.push_back(Xloc[bar]);
      lepton_y_coords_l1.push_back(Yloc[bar]);
      lepton_ex_l1.push_back(TARGET_Errors_X);
      lepton_ey_l1.push_back(TARGET_Errors_Y);
    }

  vector<double> lepton_x_coords_l2, lepton_y_coords_l2;
  vector<double> lepton_ex_l2, lepton_ey_l2;
  for (vector<int>::const_iterator i = lepton_bars_l2.begin();
    i != lepton_bars_l2.end(); ++i) {
      int bar = *i;
      lepton_x_coords_l2.push_back(Xloc[bar]);
      lepton_y_coords_l2.push_back(Yloc[bar]);
      lepton_ex_l2.push_back(TARGET_Errors_X);
      lepton_ey_l2.push_back(TARGET_Errors_Y);
    }



  //Plot of 2 tracks
  TGraphErrors *lepton_track_1;

  if (lepton_bars_l1.size() > 0) {

    lepton_track_1 = new TGraphErrors(lepton_x_coords_l1.size(),&lepton_x_coords_l1[0],
                &lepton_y_coords_l1[0], &lepton_ex_l1[0],
                &lepton_ey_l1[0] );
    // TF1 *lepton_track_fit_1 = new TF1("lepton_track_1", "pol1");
    lepton_track_1->SetMarkerStyle(21);
    lepton_track_1->SetMarkerColor(color_1);
    lepton_track_1->SetMarkerSize(0.8);
    lepton_track_1->GetXaxis()->SetLimits(-50.,50.);
    lepton_track_1->GetYaxis()->SetRangeUser(-50.,50.);
    TF1 *lepton_track_fit_1 = new TF1("lepton_track_1", "pol1");
    lepton_track_1->Fit("lepton_track_1","QW");
    lepton_track_1->Fit("lepton_track_1","QC");
    lepton_track_fit_1 = lepton_track_1->GetFunction("lepton_track_1");
    lepton_track_fit_1->SetLineWidth(2);
    lepton_track_fit_1->SetLineColor(color_1);

  }

  TGraphErrors *lepton_track_2;
  TF1 *lepton_track_fit_2;
  if (lepton_bars_l2.size() > 0) {
    lepton_track_2 = new TGraphErrors(lepton_x_coords_l2.size(),&lepton_x_coords_l2[0],
          &lepton_y_coords_l2[0], &lepton_ex_l2[0],
          &lepton_ey_l2[0]);
    // TF1 *lepton_track_fit_2 = new TF1("lepton_track_2", "pol1");

    lepton_track_2->SetMarkerStyle(21);
    lepton_track_2->SetMarkerColor(color_2);
    lepton_track_2->SetMarkerSize(0.8);
    lepton_track_2->GetXaxis()->SetLimits(-50.,50.);
    lepton_track_2->GetYaxis()->SetRangeUser(-50.,50.);
    TF1 *lepton_track_fit_2 = new TF1("lepton_track_2", "pol1");
    lepton_track_2->Fit("lepton_track_2","QW");
    lepton_track_2->Fit("lepton_track_2","QC");
    lepton_track_fit_2 = lepton_track_2->GetFunction("lepton_track_2");
    lepton_track_fit_2->SetLineWidth(2);
    lepton_track_fit_2->SetLineColor(color_1);

  }

  TH2F *h_Circle = new TH2F("Test", "TARGET and TOF1", 500, -50, 50, 500, -50, 50);
  double Gap[12][3][2] = {{{0}}};
  /*
  for(int g=0; g<12; g++){
    Gap[g][0][0] = TOF_Xloc[3*g];
    Gap[g][1][0] = TOF_Xloc[3*g+1];
    Gap[g][2][0] = TOF_Xloc[3*g+2];

    Gap[g][0][1] = TOF_Yloc[3*g];
    Gap[g][1][1] = TOF_Yloc[3*g+1];
    Gap[g][2][1] = TOF_Yloc[3*g+2];
  }



  for(int i=0; i<12; i++)
  {
    h_Circle->Fill(Gap[i][0][0], Gap[i][0][1]);
    h_Circle->Fill(Gap[i][1][0], Gap[i][1][1]);
    h_Circle->Fill(Gap[i][2][0], Gap[i][2][1]);
  }
  */

  h_Circle->SetMarkerStyle(5);
  h_Circle->SetMarkerSize(1.2);
  h_Circle->SetLineWidth(2);



  h_TOF1->SetMarkerStyle(20);           h_TOF1_closest->SetMarkerStyle(20);
  h_TOF1->SetMarkerColor(2);            h_TOF1_closest->SetMarkerColor(3);
  h_TOF1->SetMarkerSize(1.5);           h_TOF1_closest->SetMarkerSize(1.5);

  h_TOF1_rotate->SetMarkerStyle(20);
  h_TOF1_rotate->SetMarkerColor(2);
  h_TOF1_rotate->SetMarkerSize(1.5);

  TEllipse *ell = new TEllipse(0, 0, R_TOF1, 0);
  TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0);
  TEllipse *ell_L1 = new TEllipse(0, 0, R_SFT_L1, 0);
  ell->SetFillStyle(0);                 ell_Target->SetFillStyle(0);                 ell_L1->SetFillStyle(0);
  ell->SetLineColor(6);                 ell_Target->SetLineColor(1);                 ell_L1->SetLineColor(4);
  ell->SetLineWidth(1);                 ell_Target->SetLineWidth(1);                 ell_L1->SetLineWidth(1);

  const int gap_to_fit = 0;
  const int gap_to_fit_rotate = gap_to_fit;
  const int selected_TOF2 = 0;

  canvas->cd(Subplot)->Range(-50, -50, 50, 50);

  h_EMPTY->Draw(); //x

  Gap1l->Draw();
  Gap2l->Draw("same");
  Gap3l->Draw("same");
  Gap4l->Draw("same");
  Gap5l->Draw("same");
  Gap6l->Draw("same");
  Gap7l->Draw("same");
  Gap8l->Draw("same");
  Gap9l->Draw("same");
  Gap10l->Draw("same");
  Gap11l->Draw("same");
  Gap12l->Draw("Xsame");

  A1->Draw(); //x
  A2->Draw(); //x

  h_Circle->Draw("same"); //x
  //gr3_Leptons_rotate->Draw("P"); //x
  if(lepton_bars_l1.size()>0) {
    lepton_track_1->Draw("P"); //x
  }
  if(lepton_bars_l2.size()>0) {
    lepton_track_2->Draw("P");
  } //x
  //best_fit_rotate->Draw("same");

  ell->Draw("same");
  ell_Target->Draw("same");

  ell_L1->Draw("same");

  h_Target_Center->Draw("same"); //x

  vector<int> second_track;





  TH2F *k_stop_point = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  k_stop_point->Fill(dv_x*10, dv_y*10);
  k_stop_point->SetMarkerStyle(2);
  k_stop_point->SetMarkerColor(4);
  k_stop_point->SetMarkerSize(2.5);
  k_stop_point->Draw("same");

  TH2F *k_stop_experimental = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  k_stop_experimental->Fill(k_stop_x,k_stop_y);
  k_stop_experimental->SetMarkerStyle(3);
  k_stop_experimental->SetMarkerColor(6);
  k_stop_experimental->SetMarkerSize(2.5);
  k_stop_experimental->Draw("same");
}
