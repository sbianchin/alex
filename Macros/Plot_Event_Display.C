//some type definitions
struct SUBPLOT_3_RETURN {
  double X_BAR;
  double Y_BAR;
  int i_kaon_bar;
  TLatex *tex_Label_Centroid_x;
  TLatex *tex_Label_Centroid_y;
};

int kaon_fiber(double x_bar, double y_bar){

  double dist_min = 999.;
  int i_min = 999;


  for(int i=0; i<256; i++){
    if(sqrt(pow(abs(Xloc[i] - x_bar),2) + pow(abs(Yloc[i] - y_bar),2)) < dist_min){
      dist_min = sqrt(pow(abs(Xloc[i] - x_bar),2) + pow(abs(Yloc[i] - y_bar),2));
      i_min = i;
    }
  }

  return i_min;
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


void DrawTargetFrames(char ch_ADC_cut_TARGET[100],
  TLatex *tex_Legend_TARGET[36], TLatex *tex_event_TARGET,
  TMarker *palette_TARGET[10], TLatex *tex_palette_TARGET[10],
  TLatex *tex_palette_TARGET_scale,
  TLine *TOF_line1, TLine *TOF_line2,
  TLine *TOF_line3, TLine *TOF_line4,
  TLine *TOF_line5, TLine *TOF_line6,
  TLine *TOF_line7, TLine *TOF_line8,
  TLine *TOF_line9, TLine *TOF_line10,
  TLine *TOF_line11, TLine *TOF_line12,
  TLine *TOF_line13, TLine *TOF_line14,
  TLine *TOF_line15, TLine *TOF_line16,
  TLine *TOF_line17, TLine *TOF_line18,
  TLine *TOF_line19, TLine *TOF_line20,
  TLine *TOF_line21, TLine *TOF_line22,
  TLine *TOF_line23, TLine *TOF_line24) {
  TLine *hline1 = new TLine(0.38,0.14,0.62,0.14);    hline1->SetLineWidth(2);
  TLine *hline2 = new TLine(0.30,0.18,0.70,0.18);    hline2->SetLineWidth(2);
  TLine *hline3 = new TLine(0.26,0.22,0.74,0.22);    hline3->SetLineWidth(2);
  TLine *hline4 = new TLine(0.22,0.26,0.78,0.26);    hline4->SetLineWidth(2);
  TLine *hline5 = new TLine(0.18,0.30,0.82,0.30);    hline5->SetLineWidth(2);
  TLine *hline6 = new TLine(0.18,0.34,0.82,0.34);    hline6->SetLineWidth(2);
  TLine *hline7 = new TLine(0.14,0.38,0.86,0.38);    hline7->SetLineWidth(2);
  TLine *hline8 = new TLine(0.14,0.42,0.86,0.42);    hline8->SetLineWidth(2);
  TLine *hline9 = new TLine(0.14,0.46,0.86,0.46);    hline9->SetLineWidth(2);
  TLine *hline10 = new TLine(0.14,0.50,0.86,0.50);   hline10->SetLineWidth(2);
  TLine *hline11 = new TLine(0.14,0.54,0.86,0.54);   hline11->SetLineWidth(2);
  TLine *hline12 = new TLine(0.14,0.58,0.86,0.58);   hline12->SetLineWidth(2);
  TLine *hline13 = new TLine(0.14,0.62,0.86,0.62);   hline13->SetLineWidth(2);
  TLine *hline14 = new TLine(0.18,0.66,0.82,0.66);   hline14->SetLineWidth(2);
  TLine *hline15 = new TLine(0.18,0.70,0.82,0.70);   hline15->SetLineWidth(2);
  TLine *hline16 = new TLine(0.22,0.74,0.78,0.74);   hline16->SetLineWidth(2);
  TLine *hline17 = new TLine(0.26,0.78,0.74,0.78);   hline17->SetLineWidth(2);
  TLine *hline18 = new TLine(0.30,0.82,0.70,0.82);   hline18->SetLineWidth(2);
  TLine *hline19 = new TLine(0.38,0.86,0.62,0.86);   hline19->SetLineWidth(2);

  TLine *vline1 = new TLine(0.14,0.38,0.14,0.62);    vline1->SetLineWidth(2);
  TLine *vline2 = new TLine(0.18,0.30,0.18,0.70);    vline2->SetLineWidth(2);
  TLine *vline3 = new TLine(0.22,0.26,0.22,0.74);    vline3->SetLineWidth(2);
  TLine *vline4 = new TLine(0.26,0.22,0.26,0.78);    vline4->SetLineWidth(2);
  TLine *vline5 = new TLine(0.30,0.18,0.30,0.82);    vline5->SetLineWidth(2);
  TLine *vline6 = new TLine(0.34,0.18,0.34,0.82);    vline6->SetLineWidth(2);
  TLine *vline7 = new TLine(0.38,0.14,0.38,0.86);    vline7->SetLineWidth(2);
  TLine *vline8 = new TLine(0.42,0.14,0.42,0.86);    vline8->SetLineWidth(2);
  TLine *vline9 = new TLine(0.46,0.14,0.46,0.86);    vline9->SetLineWidth(2);
  TLine *vline10 = new TLine(0.50,0.14,0.50,0.86);   vline10->SetLineWidth(2);
  TLine *vline11 = new TLine(0.54,0.14,0.54,0.86);   vline11->SetLineWidth(2);
  TLine *vline12 = new TLine(0.58,0.14,0.58,0.86);   vline12->SetLineWidth(2);
  TLine *vline13 = new TLine(0.62,0.14,0.62,0.86);   vline13->SetLineWidth(2);
  TLine *vline14 = new TLine(0.66,0.18,0.66,0.82);   vline14->SetLineWidth(2);
  TLine *vline15 = new TLine(0.70,0.18,0.70,0.82);   vline15->SetLineWidth(2);
  TLine *vline16 = new TLine(0.74,0.22,0.74,0.78);   vline16->SetLineWidth(2);
  TLine *vline17 = new TLine(0.78,0.26,0.78,0.74);   vline17->SetLineWidth(2);
  TLine *vline18 = new TLine(0.82,0.30,0.82,0.70);   vline18->SetLineWidth(2);
  TLine *vline19 = new TLine(0.86,0.38,0.86,0.62);   vline19->SetLineWidth(2);

  TLine *vblue1 = new TLine(0.38,0.70,0.38,0.74);      vblue1->SetLineWidth(5);   vblue1->SetLineColor(kBlue-9);
  TLine *vblue2 = new TLine(0.42,0.62,0.42,0.66);      vblue2->SetLineWidth(5);   vblue2->SetLineColor(kBlue-9);
  TLine *vblue3 = new TLine(0.30,0.54,0.30,0.58);      vblue3->SetLineWidth(5);   vblue3->SetLineColor(kBlue-9);
  TLine *vblue4 = new TLine(0.86,0.50,0.86,0.54);      vblue4->SetLineWidth(5);   vblue4->SetLineColor(kBlue-9);
  TLine *vblue5 = new TLine(0.70,0.42,0.70,0.46);      vblue5->SetLineWidth(5);   vblue5->SetLineColor(kBlue-9);
  TLine *vblue6 = new TLine(0.58,0.34,0.58,0.38);      vblue6->SetLineWidth(5);   vblue6->SetLineColor(kBlue-9);
  TLine *vblue7 = new TLine(0.62,0.26,0.62,0.30);      vblue7->SetLineWidth(5);   vblue7->SetLineColor(kBlue-9);
  TLine *vblue8 = new TLine(0.62,0.14,0.62,0.18);      vblue8->SetLineWidth(5);   vblue8->SetLineColor(kBlue-9);

  hline1->Draw();     vline1->Draw();   TOF_line13->Draw();   TOF_line1->Draw();
  hline2->Draw();     vline2->Draw();   TOF_line14->Draw();   TOF_line2->Draw();
  hline3->Draw();     vline3->Draw();   TOF_line15->Draw();   TOF_line3->Draw();
  hline4->Draw();     vline4->Draw();   TOF_line16->Draw();   TOF_line4->Draw();
  hline5->Draw();     vline5->Draw();   TOF_line17->Draw();   TOF_line5->Draw();
  hline6->Draw();     vline6->Draw();   TOF_line18->Draw();   TOF_line6->Draw();
  hline7->Draw();     vline7->Draw();   TOF_line19->Draw();   TOF_line7->Draw();
  hline8->Draw();     vline8->Draw();   TOF_line20->Draw();   TOF_line8->Draw();
  hline9->Draw();     vline9->Draw();   TOF_line21->Draw();   TOF_line9->Draw();
  hline10->Draw();    vline10->Draw();  TOF_line22->Draw();   TOF_line10->Draw();
  hline11->Draw();    vline11->Draw();  TOF_line23->Draw();   TOF_line11->Draw();
  hline12->Draw();    vline12->Draw();  TOF_line24->Draw();   TOF_line12->Draw();
  hline13->Draw();    vline13->Draw();
  hline14->Draw();    vline14->Draw();
  hline15->Draw();    vline15->Draw();
  hline16->Draw();    vline16->Draw();
  hline17->Draw();    vline17->Draw();
  hline18->Draw();    vline18->Draw();
  hline19->Draw();    vline19->Draw();

  vblue2->Draw();
  vblue4->Draw();
  vblue6->Draw();
  vblue8->Draw();

  for(Int_t ileg=0; ileg<36; ileg++)    tex_Legend_TARGET[ileg]->Draw();

  tex_event_TARGET->Draw();

  TLatex *tex_Title_ADC_High_TARGET;
  tex_Title_ADC_High_TARGET = new TLatex(0.01759134,0.9295171,"ADC HG Cut");
  tex_Title_ADC_High_TARGET->SetTextSize(0.07645875);
  tex_Title_ADC_High_TARGET->SetLineWidth(2);
  tex_Title_ADC_High_TARGET->Draw();

  TLatex *tex_Subtitle_ADC_High_TARGET;
  tex_Subtitle_ADC_High_TARGET = new TLatex(0.01759134,0.88,ch_ADC_cut_TARGET);
  tex_Subtitle_ADC_High_TARGET->SetTextSize(0.04);
  tex_Subtitle_ADC_High_TARGET->SetLineWidth(2);
  tex_Subtitle_ADC_High_TARGET->Draw();

  for(Int_t ipal=0; ipal<10; ipal++) palette_TARGET[ipal]->Draw();
  for(Int_t ileg=0; ileg<10; ileg++) tex_palette_TARGET[ileg]->Draw();
  tex_palette_TARGET_scale->Draw();
}

void Subplot_1(char ch_ADC_cut_TARGET[100],
  TLatex *tex_Legend_TARGET[36], TLatex *tex_event_TARGET,
  Int_t ADC_High_TARGET[256], TMarker *marker_ADC_TARGET[256],
  Int_t ADC_Low_TARGET[256], Int_t TDC_min_TARGET,
  Int_t Switch, bool has_TDC_hit[256],
  TMarker *palette_TARGET[10], TLatex *tex_palette_TARGET[10],
  TLatex *tex_palette_TARGET_scale,
  int max_index, int max_index2, int max_index3, int max_index4,
  TMarker *marker_TDC_TARGET[256], TLine *TOF_line1, TLine *TOF_line2,
  TLine *TOF_line3, TLine *TOF_line4,
  TLine *TOF_line5, TLine *TOF_line6,
  TLine *TOF_line7, TLine *TOF_line8,
  TLine *TOF_line9, TLine *TOF_line10,
  TLine *TOF_line11, TLine *TOF_line12,
  TLine *TOF_line13, TLine *TOF_line14,
  TLine *TOF_line15, TLine *TOF_line16,
  TLine *TOF_line17, TLine *TOF_line18,
  TLine *TOF_line19, TLine *TOF_line20,
  TLine *TOF_line21, TLine *TOF_line22,
  TLine *TOF_line23, TLine *TOF_line24) {

  DrawTargetFrames(ch_ADC_cut_TARGET, tex_Legend_TARGET, tex_event_TARGET,
    palette_TARGET, tex_palette_TARGET,
    tex_palette_TARGET_scale, TOF_line1, TOF_line2,
    TOF_line3, TOF_line4,
    TOF_line5, TOF_line6,
    TOF_line7, TOF_line8,
    TOF_line9, TOF_line10,
    TOF_line11, TOF_line12,
    TOF_line13, TOF_line14,
    TOF_line15, TOF_line16,
    TOF_line17, TOF_line18,
    TOF_line19, TOF_line20,
    TOF_line21, TOF_line22,
    TOF_line23, TOF_line24);





  for(Int_t icol=0; icol<256; icol++){
    if(ADC_High_TARGET[icol]>=0 && ADC_High_TARGET[icol]<300){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange+10); marker_ADC_TARGET[icol]->Draw();}
    if(ADC_High_TARGET[icol]>=300 && ADC_High_TARGET[icol]<600){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange+7); marker_ADC_TARGET[icol]->Draw();}
    if(ADC_High_TARGET[icol]>=600 && ADC_High_TARGET[icol]<900){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange+1); marker_ADC_TARGET[icol]->Draw();}
    if(ADC_High_TARGET[icol]>=900 && ADC_High_TARGET[icol]<1200){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange-4); marker_ADC_TARGET[icol]->Draw();}
    if(ADC_High_TARGET[icol]>=1200 && ADC_High_TARGET[icol]<1500){ marker_ADC_TARGET[icol]->SetMarkerColor(kYellow-9); marker_ADC_TARGET[icol]->Draw();}
    if(ADC_High_TARGET[icol]>=1500 && ADC_High_TARGET[icol]<1800){ marker_ADC_TARGET[icol]->SetMarkerColor(kYellow-7); marker_ADC_TARGET[icol]->Draw();}
    if(ADC_High_TARGET[icol]>=1800 && ADC_High_TARGET[icol]<2100){ marker_ADC_TARGET[icol]->SetMarkerColor(kYellow-0); marker_ADC_TARGET[icol]->Draw();}
    if(ADC_High_TARGET[icol]>=2100 && ADC_High_TARGET[icol]<2400){ marker_ADC_TARGET[icol]->SetMarkerColor(kSpring-4); marker_ADC_TARGET[icol]->Draw();}
    if(ADC_High_TARGET[icol]>=2400 && ADC_High_TARGET[icol]<2700){ marker_ADC_TARGET[icol]->SetMarkerColor(kSpring-2); marker_ADC_TARGET[icol]->Draw();}
    if(ADC_High_TARGET[icol]>=2700){ marker_ADC_TARGET[icol]->SetMarkerColor(kGreen-0); marker_ADC_TARGET[icol]->Draw();}

    if(ADC_High_TARGET[icol]<0 && ADC_Low_TARGET[icol]>0 && TDC_min_TARGET && Switch==1){
      if(has_TDC_hit[icol]) marker_ADC_TARGET[icol]->SetMarkerColor(11);
      else marker_ADC_TARGET[icol]->SetMarkerColor(kBlue);
      marker_ADC_TARGET[icol]->Draw();
    }
  }

  for(Int_t ipal=0; ipal<10; ipal++) palette_TARGET[ipal]->Draw();
  for(Int_t ileg=0; ileg<10; ileg++) tex_palette_TARGET[ileg]->Draw();
  tex_palette_TARGET_scale->Draw();

  if(ADC_High_TARGET[max_index] > 0){
    marker_TDC_TARGET[max_index]->SetMarkerColor(kViolet+1);
    marker_TDC_TARGET[max_index]->Draw();
  }

  if(ADC_High_TARGET[max_index2] > 0){
    marker_TDC_TARGET[max_index2]->SetMarkerColor(1);
    marker_TDC_TARGET[max_index2]->Draw();
  }

  if(ADC_High_TARGET[max_index3] > 0){
    marker_TDC_TARGET[max_index3]->SetMarkerColor(1);
    marker_TDC_TARGET[max_index3]->Draw();
  }

  if(ADC_High_TARGET[max_index4] > 0){
    marker_TDC_TARGET[max_index4]->SetMarkerColor(1);
    marker_TDC_TARGET[max_index4]->Draw();
  }
}

void Subplot_2(char ch_ADC_cut_TARGET[100],
  TLatex *tex_Legend_TARGET[36], TLatex *tex_event_TARGET,
  Int_t ADC_High_TARGET[256], TMarker *marker_ADC_TARGET[256],
  Int_t ADC_Low_TARGET[256], Int_t TDC_min_TARGET,
  Int_t Switch, bool has_TDC_hit[256],
  TMarker *palette_TARGET[10], TLatex *tex_palette_TARGET[10],
  TLatex *tex_palette_TARGET_scale,
  int max_index, int max_index2, int max_index3, int max_index4,
  TMarker *marker_TDC_TARGET[256], TLine *TOF_line1, TLine *TOF_line2,
  TLine *TOF_line3, TLine *TOF_line4,
  TLine *TOF_line5, TLine *TOF_line6,
  TLine *TOF_line7, TLine *TOF_line8,
  TLine *TOF_line9, TLine *TOF_line10,
  TLine *TOF_line11, TLine *TOF_line12,
  TLine *TOF_line13, TLine *TOF_line14,
  TLine *TOF_line15, TLine *TOF_line16,
  TLine *TOF_line17, TLine *TOF_line18,
  TLine *TOF_line19, TLine *TOF_line20,
  TLine *TOF_line21, TLine *TOF_line22,
  TLine *TOF_line23, TLine *TOF_line24) {

  DrawTargetFrames(ch_ADC_cut_TARGET, tex_Legend_TARGET, tex_event_TARGET,
    palette_TARGET, tex_palette_TARGET,
    tex_palette_TARGET_scale, TOF_line1, TOF_line2,
    TOF_line3, TOF_line4,
    TOF_line5, TOF_line6,
    TOF_line7, TOF_line8,
    TOF_line9, TOF_line10,
    TOF_line11, TOF_line12,
    TOF_line13, TOF_line14,
    TOF_line15, TOF_line16,
    TOF_line17, TOF_line18,
    TOF_line19, TOF_line20,
    TOF_line21, TOF_line22,
    TOF_line23, TOF_line24);



  for(Int_t icol=0; icol<256; icol++){
    if(has_TDC_hit[icol]) {
      if(ADC_High_TARGET[icol]>=0 && ADC_High_TARGET[icol]<300){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange+10); marker_ADC_TARGET[icol]->Draw();}
      if(ADC_High_TARGET[icol]>=300 && ADC_High_TARGET[icol]<600){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange+7); marker_ADC_TARGET[icol]->Draw();}
      if(ADC_High_TARGET[icol]>=600 && ADC_High_TARGET[icol]<900){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange+1); marker_ADC_TARGET[icol]->Draw();}
      if(ADC_High_TARGET[icol]>=900 && ADC_High_TARGET[icol]<1200){ marker_ADC_TARGET[icol]->SetMarkerColor(kOrange-4); marker_ADC_TARGET[icol]->Draw();}
      if(ADC_High_TARGET[icol]>=1200 && ADC_High_TARGET[icol]<1500){ marker_ADC_TARGET[icol]->SetMarkerColor(kYellow-9); marker_ADC_TARGET[icol]->Draw();}
      if(ADC_High_TARGET[icol]>=1500 && ADC_High_TARGET[icol]<1800){ marker_ADC_TARGET[icol]->SetMarkerColor(kYellow-7); marker_ADC_TARGET[icol]->Draw();}
      if(ADC_High_TARGET[icol]>=1800 && ADC_High_TARGET[icol]<2100){ marker_ADC_TARGET[icol]->SetMarkerColor(kYellow-0); marker_ADC_TARGET[icol]->Draw();}
      if(ADC_High_TARGET[icol]>=2100 && ADC_High_TARGET[icol]<2400){ marker_ADC_TARGET[icol]->SetMarkerColor(kSpring-4); marker_ADC_TARGET[icol]->Draw();}
      if(ADC_High_TARGET[icol]>=2400 && ADC_High_TARGET[icol]<2700){ marker_ADC_TARGET[icol]->SetMarkerColor(kSpring-2); marker_ADC_TARGET[icol]->Draw();}
      if(ADC_High_TARGET[icol]>=2700){ marker_ADC_TARGET[icol]->SetMarkerColor(kGreen-0); marker_ADC_TARGET[icol]->Draw();}
    }

    if(ADC_High_TARGET[icol]<0 && ADC_Low_TARGET[icol]>0 && TDC_min_TARGET && Switch==1){
      if(has_TDC_hit[icol]) marker_ADC_TARGET[icol]->SetMarkerColor(11);
      else marker_ADC_TARGET[icol]->SetMarkerColor(kBlue);
      marker_ADC_TARGET[icol]->Draw();
    }
  }



  if((ADC_Low_TARGET[max_index] > 0) && (ADC_High_TARGET[max_index] > 0) && (has_TDC_hit[max_index])){
   	marker_TDC_TARGET[max_index]->SetMarkerColor(kViolet+1);
    marker_TDC_TARGET[max_index]->Draw();
  }

  if((ADC_Low_TARGET[max_index2] > 0) && (ADC_High_TARGET[max_index2] > 0) && (has_TDC_hit[max_index2])){
    marker_TDC_TARGET[max_index2]->SetMarkerColor(1);
    marker_TDC_TARGET[max_index2]->Draw();
  }

  if((ADC_Low_TARGET[max_index3] > 0) && (ADC_High_TARGET[max_index3] > 0) && (has_TDC_hit[max_index3])){
    marker_TDC_TARGET[max_index3]->SetMarkerColor(1);
    marker_TDC_TARGET[max_index3]->Draw();
  }

  if((ADC_Low_TARGET[max_index4] > 0) && (ADC_High_TARGET[max_index4] > 0) && (has_TDC_hit[max_index4])){
    marker_TDC_TARGET[max_index4]->SetMarkerColor(1);
    marker_TDC_TARGET[max_index4]->Draw();
  }
}

void Subplot_3(
  TLatex *tex_Legend_TARGET[36], TLatex *tex_event_TARGET,
  char ch_ADC_and_TDC_cut_Kstop[100],
  Int_t ADC_High_TARGET[256], TMarker *marker_ADCL_TARGET[256],
  Int_t ADC_Low_TARGET[256],
  bool has_TDC_hit_Kstop[256],
  TLatex *tex_version2,
  TLine *TOF_line1, TLine *TOF_line2,
  TLine *TOF_line3, TLine *TOF_line4,
  TLine *TOF_line5, TLine *TOF_line6,
  TLine *TOF_line7, TLine *TOF_line8,
  TLine *TOF_line9, TLine *TOF_line10,
  TLine *TOF_line11, TLine *TOF_line12,
  TLine *TOF_line13, TLine *TOF_line14,
  TLine *TOF_line15, TLine *TOF_line16,
  TLine *TOF_line17, TLine *TOF_line18,
  TLine *TOF_line19, TLine *TOF_line20,
  TLine *TOF_line21, TLine *TOF_line22,
  TLine *TOF_line23, TLine *TOF_line24) {

  TLine *hline1 = new TLine(0.38,0.14,0.62,0.14);    hline1->SetLineWidth(2);
  TLine *hline2 = new TLine(0.30,0.18,0.70,0.18);    hline2->SetLineWidth(2);
  TLine *hline3 = new TLine(0.26,0.22,0.74,0.22);    hline3->SetLineWidth(2);
  TLine *hline4 = new TLine(0.22,0.26,0.78,0.26);    hline4->SetLineWidth(2);
  TLine *hline5 = new TLine(0.18,0.30,0.82,0.30);    hline5->SetLineWidth(2);
  TLine *hline6 = new TLine(0.18,0.34,0.82,0.34);    hline6->SetLineWidth(2);
  TLine *hline7 = new TLine(0.14,0.38,0.86,0.38);    hline7->SetLineWidth(2);
  TLine *hline8 = new TLine(0.14,0.42,0.86,0.42);    hline8->SetLineWidth(2);
  TLine *hline9 = new TLine(0.14,0.46,0.86,0.46);    hline9->SetLineWidth(2);
  TLine *hline10 = new TLine(0.14,0.50,0.86,0.50);   hline10->SetLineWidth(2);
  TLine *hline11 = new TLine(0.14,0.54,0.86,0.54);   hline11->SetLineWidth(2);
  TLine *hline12 = new TLine(0.14,0.58,0.86,0.58);   hline12->SetLineWidth(2);
  TLine *hline13 = new TLine(0.14,0.62,0.86,0.62);   hline13->SetLineWidth(2);
  TLine *hline14 = new TLine(0.18,0.66,0.82,0.66);   hline14->SetLineWidth(2);
  TLine *hline15 = new TLine(0.18,0.70,0.82,0.70);   hline15->SetLineWidth(2);
  TLine *hline16 = new TLine(0.22,0.74,0.78,0.74);   hline16->SetLineWidth(2);
  TLine *hline17 = new TLine(0.26,0.78,0.74,0.78);   hline17->SetLineWidth(2);
  TLine *hline18 = new TLine(0.30,0.82,0.70,0.82);   hline18->SetLineWidth(2);
  TLine *hline19 = new TLine(0.38,0.86,0.62,0.86);   hline19->SetLineWidth(2);

  TLine *vline1 = new TLine(0.14,0.38,0.14,0.62);    vline1->SetLineWidth(2);
  TLine *vline2 = new TLine(0.18,0.30,0.18,0.70);    vline2->SetLineWidth(2);
  TLine *vline3 = new TLine(0.22,0.26,0.22,0.74);    vline3->SetLineWidth(2);
  TLine *vline4 = new TLine(0.26,0.22,0.26,0.78);    vline4->SetLineWidth(2);
  TLine *vline5 = new TLine(0.30,0.18,0.30,0.82);    vline5->SetLineWidth(2);
  TLine *vline6 = new TLine(0.34,0.18,0.34,0.82);    vline6->SetLineWidth(2);
  TLine *vline7 = new TLine(0.38,0.14,0.38,0.86);    vline7->SetLineWidth(2);
  TLine *vline8 = new TLine(0.42,0.14,0.42,0.86);    vline8->SetLineWidth(2);
  TLine *vline9 = new TLine(0.46,0.14,0.46,0.86);    vline9->SetLineWidth(2);
  TLine *vline10 = new TLine(0.50,0.14,0.50,0.86);   vline10->SetLineWidth(2);
  TLine *vline11 = new TLine(0.54,0.14,0.54,0.86);   vline11->SetLineWidth(2);
  TLine *vline12 = new TLine(0.58,0.14,0.58,0.86);   vline12->SetLineWidth(2);
  TLine *vline13 = new TLine(0.62,0.14,0.62,0.86);   vline13->SetLineWidth(2);
  TLine *vline14 = new TLine(0.66,0.18,0.66,0.82);   vline14->SetLineWidth(2);
  TLine *vline15 = new TLine(0.70,0.18,0.70,0.82);   vline15->SetLineWidth(2);
  TLine *vline16 = new TLine(0.74,0.22,0.74,0.78);   vline16->SetLineWidth(2);
  TLine *vline17 = new TLine(0.78,0.26,0.78,0.74);   vline17->SetLineWidth(2);
  TLine *vline18 = new TLine(0.82,0.30,0.82,0.70);   vline18->SetLineWidth(2);
  TLine *vline19 = new TLine(0.86,0.38,0.86,0.62);   vline19->SetLineWidth(2);

  TLine *vblue1 = new TLine(0.38,0.70,0.38,0.74);      vblue1->SetLineWidth(5);   vblue1->SetLineColor(kBlue-9);
  TLine *vblue2 = new TLine(0.42,0.62,0.42,0.66);      vblue2->SetLineWidth(5);   vblue2->SetLineColor(kBlue-9);
  TLine *vblue3 = new TLine(0.30,0.54,0.30,0.58);      vblue3->SetLineWidth(5);   vblue3->SetLineColor(kBlue-9);
  TLine *vblue4 = new TLine(0.86,0.50,0.86,0.54);      vblue4->SetLineWidth(5);   vblue4->SetLineColor(kBlue-9);
  TLine *vblue5 = new TLine(0.70,0.42,0.70,0.46);      vblue5->SetLineWidth(5);   vblue5->SetLineColor(kBlue-9);
  TLine *vblue6 = new TLine(0.58,0.34,0.58,0.38);      vblue6->SetLineWidth(5);   vblue6->SetLineColor(kBlue-9);
  TLine *vblue7 = new TLine(0.62,0.26,0.62,0.30);      vblue7->SetLineWidth(5);   vblue7->SetLineColor(kBlue-9);
  TLine *vblue8 = new TLine(0.62,0.14,0.62,0.18);      vblue8->SetLineWidth(5);   vblue8->SetLineColor(kBlue-9);

  hline1->Draw();     vline1->Draw();   TOF_line13->Draw();   TOF_line1->Draw();
  hline2->Draw();     vline2->Draw();   TOF_line14->Draw();   TOF_line2->Draw();
  hline3->Draw();     vline3->Draw();   TOF_line15->Draw();   TOF_line3->Draw();
  hline4->Draw();     vline4->Draw();   TOF_line16->Draw();   TOF_line4->Draw();
  hline5->Draw();     vline5->Draw();   TOF_line17->Draw();   TOF_line5->Draw();
  hline6->Draw();     vline6->Draw();   TOF_line18->Draw();   TOF_line6->Draw();
  hline7->Draw();     vline7->Draw();   TOF_line19->Draw();   TOF_line7->Draw();
  hline8->Draw();     vline8->Draw();   TOF_line20->Draw();   TOF_line8->Draw();
  hline9->Draw();     vline9->Draw();   TOF_line21->Draw();   TOF_line9->Draw();
  hline10->Draw();    vline10->Draw();  TOF_line22->Draw();   TOF_line10->Draw();
  hline11->Draw();    vline11->Draw();  TOF_line23->Draw();   TOF_line11->Draw();
  hline12->Draw();    vline12->Draw();  TOF_line24->Draw();   TOF_line12->Draw();
  hline13->Draw();    vline13->Draw();
  hline14->Draw();    vline14->Draw();
  hline15->Draw();    vline15->Draw();
  hline16->Draw();    vline16->Draw();
  hline17->Draw();    vline17->Draw();
  hline18->Draw();    vline18->Draw();
  hline19->Draw();    vline19->Draw();

  vblue2->Draw();
  vblue4->Draw();
  vblue6->Draw();
  vblue8->Draw();

  tex_version2->Draw();

  for(Int_t ileg=0; ileg<36; ileg++) tex_Legend_TARGET[ileg]->Draw();
  tex_event_TARGET->Draw();

  TLatex *tex_Title_ADC_Low_TARGET = new TLatex(0.01759134,0.9295171,"K Stop Cuts");
  tex_Title_ADC_Low_TARGET->SetTextSize(0.07);
  tex_Title_ADC_Low_TARGET->SetLineWidth(2);
  tex_Title_ADC_Low_TARGET->Draw();

  TLatex *tex_Subtitle_ADC_Low_TARGET = new TLatex(0.01759134,0.88,ch_ADC_and_TDC_cut_Kstop);
  tex_Subtitle_ADC_Low_TARGET->SetTextSize(0.04);
  tex_Subtitle_ADC_Low_TARGET->SetLineWidth(2);
  tex_Subtitle_ADC_Low_TARGET->Draw();

  //Determine closest bar to the centroid, and draw in red
  //Not used elsewhere
  double X_BAR = 999.;
  double Y_BAR = 999.;
  double closest_to_centroid = 1000;
  int closest_to_centroid_index = -1;

  for(int i = 0; i<256; i++){
  if(distance(X_BAR,Y_BAR,Xloc[i],Yloc[i]) <= closest_to_centroid){
    closest_to_centroid_index = i;
    closest_to_centroid = distance(X_BAR,Y_BAR,Xloc[i],Yloc[i]);
    }
  }

  // 3rd window will show the condition for a kstop event
  for(Int_t icol=0; icol<256; icol++){
    if (has_TDC_hit_Kstop[icol]){
      if(ADC_Low_TARGET[icol]>LG_KAON && ADC_High_TARGET[icol]>HG_KAON){
        marker_ADCL_TARGET[icol]->SetMarkerColor(12); marker_ADCL_TARGET[icol]->Draw();
      }
    }
  }


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



// void Subplot_4(TH2F *h_lepton_fit_final, TGaxis *A1,
//   TGaxis *A2, TGraph *gr_Target_Center,
//   TEllipse *ell, TEllipse *ell_Target,
//   TEllipse *ell_L1, int gap_to_fit,
//   TGraph gr3_Leptons_rotate, TGraph gr_lepton_3,
//   TLine *best_fit_rotate,
//   TGraph *gr_TOF1, TGraph *gr_TOF1_closest) {
//   h_lepton_fit_final->Draw();
//   A1->Draw();
//   A2->Draw();
//   gr_Target_Center->Draw("P");
//
//   //gr_TOF1_Markers->Draw("P");
//
//   Gap1l->Draw();
//   Gap2l->Draw("same");
//   Gap3l->Draw("same");
//   Gap4l->Draw("same");
//   Gap5l->Draw("same");
//   Gap6l->Draw("same");
//   Gap7l->Draw("same");
//   Gap8l->Draw("same");
//   Gap9l->Draw("same");
//   Gap10l->Draw("same");
//   Gap11l->Draw("same");
//   Gap12l->Draw("Xsame");
//
//   ell->Draw("same");
//   ell_Target->Draw("same");
//   ell_L1->Draw("same");
//   gr_TOF1->Draw("P");
//   gr_TOF1_closest->Draw("P");
//   //if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){ //ROTATE_CHANGE
//   if(gap_to_fit==6 || gap_to_fit==12){
//     gr3_Leptons_rotate.Draw("Psame");
//     best_fit_rotate->Draw("same");
//   }
//   else gr_lepton_3.Draw("P");
//
//
// }
//
//
//
//
// void Subplot_5(TH2F *h_kaon_fit_final, TGaxis *A1,
//   TGaxis *A2, TGraph *gr_Target_Center,
//   TEllipse *ell, TEllipse *ell_Target,
//   TEllipse *ell_L1,
//   TGraph gr_kaon, TGraph gr_kaon_bk,
//   int vec_xx_kaon_size, bool kaon_bk) {
//
//     h_kaon_fit_final->Draw();
//     A1->Draw();
//     A2->Draw();
//     gr_Target_Center->Draw("P");
//
//     //gr_TOF1_Markers->Draw("P");
//
//     Gap1l->Draw();
//     Gap2l->Draw("same");
//     Gap3l->Draw("same");
//     Gap4l->Draw("same");
//     Gap5l->Draw("same");
//     Gap6l->Draw("same");
//     Gap7l->Draw("same");
//     Gap8l->Draw("same");
//     Gap9l->Draw("same");
//     Gap10l->Draw("same");
//     Gap11l->Draw("same");
//     Gap12l->Draw("Xsame");
//
//     ell->Draw("same");
//     ell_Target->Draw("same");
//     ell_L1->Draw("same");
//
//     if(vec_xx_kaon_size >0){
//       if(kaon_bk) gr_kaon_bk.Draw("P");
//       else gr_kaon.Draw("P");
//     }
// }
//
// void Subplot_6(TH2F *h_final, TGaxis *A1,
//   TGaxis *A2, TGraph *gr_Target_Center,
//   TEllipse *ell, TEllipse *ell_Target,
//   TEllipse *ell_L1,
//   int gap_to_fit,
//   TGraph *gr_TOF1, TGraph *gr_TOF1_closest,
//   TGraph *gr_int_TDC_TARGET, TGraph *gr_int_TDC_Gap_Fibers_SFT,
//   TGraph *gr_kaon_stop,
//   TArrow *x_guide, TArrow *y_guide,
//   TLatex *tex_Angle_guide, TLatex *tex_ChiS,
//   TLatex *tex_Kstop_X, TLatex *tex_Kstop_Y,
//   TLatex *x_sft, TLatex *y_sft, TLatex *x_target, TLatex *y_target,
//   TGraph gr3_Leptons_rotate, TGraph gr_lepton_3,
//   TLine *best_fit_rotate,
//   TGraph gr_kaon, TGraph gr_kaon_bk,
//   int vec_xx_kaon_size, bool kaon_bk) {
//
//     h_final->Draw();
//     A1->Draw();
//     A2->Draw();
//     gr_Target_Center->Draw("P");
//
//     //gr_TOF1_Markers->Draw("P");
//
//     Gap1l->Draw();
//     Gap2l->Draw("same");
//     Gap3l->Draw("same");
//     Gap4l->Draw("same");
//     Gap5l->Draw("same");
//     Gap6l->Draw("same");
//     Gap7l->Draw("same");
//     Gap8l->Draw("same");
//     Gap9l->Draw("same");
//     Gap10l->Draw("same");
//     Gap11l->Draw("same");
//     Gap12l->Draw("same");
//
//     ell->Draw("same");
//     ell_Target->Draw("same");
//     ell_L1->Draw("same");
//
//     x_guide->Draw();
//     y_guide->Draw();
//
//     //if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
//     if(gap_to_fit==6 || gap_to_fit==12){
//       gr3_Leptons_rotate.Draw("P");
//       best_fit_rotate->Draw("same");
//     }
//     else{
//       gr_lepton_3.Draw("P");
//     }
//
//     if(vec_xx_kaon_size >0){
//       if(kaon_bk) gr_kaon_bk.Draw("P");
//       else gr_kaon.Draw("P");
//     }
//
//     gr_TOF1->Draw("P");
//     gr_TOF1_closest->Draw("P");
//
//     gr_int_TDC_TARGET->Draw("P");
//     gr_int_TDC_Gap_Fibers_SFT->Draw("P");
//     //gr_int_TDC_Gap_Fibers->Draw("P");
//     gr_kaon_stop->Draw("P");
//
//     tex_Angle_guide->Draw("same");
//     tex_ChiS->Draw("same");
//     tex_Kstop_X->Draw("same");
//     tex_Kstop_Y->Draw("same");
//     x_sft->Draw("same");
//     y_sft->Draw("same");
//     x_target->Draw("same");
//     y_target->Draw("same");
// }

void Subplot_4(TH2F *h_lepton_fit_final, TGaxis *A1,
  TGaxis *A2, TGraph *gr_Target_Center,
  TEllipse *ell, TEllipse *ell_Target,
  TEllipse *ell_L1, int gap_to_fit,
  TGraph *gr3_Leptons_rotate, TGraph *gr_lepton_3,
  TLine *best_fit_rotate,
  TGraph *gr_TOF1, TGraph *gr_TOF1_closest) {
  h_lepton_fit_final->Draw();
  A1->Draw();
  A2->Draw();
  gr_Target_Center->Draw("P");

  //gr_TOF1_Markers->Draw("P");

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

  ell->Draw("same");
  ell_Target->Draw("same");
  ell_L1->Draw("same");
  gr_TOF1->Draw("P");
  gr_TOF1_closest->Draw("P");
  //if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){ //ROTATE_CHANGE
  if(gap_to_fit==6 || gap_to_fit==12){
    gr3_Leptons_rotate->Draw("Psame");
    best_fit_rotate->Draw("same");
  }
  else gr_lepton_3->Draw("P");


}




void Subplot_5(TH2F *h_kaon_fit_final, TGaxis *A1,
  TGaxis *A2, TGraph *gr_Target_Center,
  TEllipse *ell, TEllipse *ell_Target,
  TEllipse *ell_L1,
  TGraph *gr_kaon, TGraph *gr_kaon_bk,
  vector<double> vec_xx_kaon, bool kaon_bk,
  bool vertical_kaon) {

    h_kaon_fit_final->Draw();
    A1->Draw();
    A2->Draw();
    gr_Target_Center->Draw("P");

    //gr_TOF1_Markers->Draw("P");

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

    ell->Draw("same");
    ell_Target->Draw("same");
    ell_L1->Draw("same");

    if(vec_xx_kaon.size() >0){
      if (vertical_kaon) {
        gr_kaon->Draw("P");
        TLine *kaon_fit = new TLine(vec_xx_kaon[0], 50., vec_xx_kaon[0], -50.);
        kaon_fit->SetLineWidth(2);
        kaon_fit->SetLineColor(4);
        kaon_fit->Draw();
      }
      else {
        if(kaon_bk) gr_kaon_bk->Draw("P");
        else gr_kaon->Draw("P");
      }
    }
}

void Subplot_6(TH2F *h_final, TGaxis *A1,
  TGaxis *A2, TGraph *gr_Target_Center,
  TEllipse *ell, TEllipse *ell_Target,
  TEllipse *ell_L1,
  int gap_to_fit,
  TGraph *gr_TOF1, TGraph *gr_TOF1_closest,
  TGraph *gr_int_TDC_TARGET, TGraph *gr_int_TDC_Gap_Fibers_SFT,
  TGraph *gr_kaon_stop,
  TArrow *x_guide, TArrow *y_guide,
  TLatex *tex_Angle_guide, TLatex *tex_ChiS,
  TLatex *tex_Kstop_X, TLatex *tex_Kstop_Y,
  TLatex *x_sft, TLatex *y_sft, TLatex *x_target, TLatex *y_target,
  TGraph *gr3_Leptons_rotate, TGraph *gr_lepton_3,
  TLine *best_fit_rotate,
  TGraph *gr_kaon, TGraph *gr_kaon_bk,
  vector<double> vec_xx_kaon, bool kaon_bk,
  bool vertical_kaon) {

    h_final->Draw();
    A1->Draw();
    A2->Draw();
    gr_Target_Center->Draw("P");

    //gr_TOF1_Markers->Draw("P");

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
    Gap12l->Draw("same");

    ell->Draw("same");
    ell_Target->Draw("same");
    ell_L1->Draw("same");

    x_guide->Draw();
    y_guide->Draw();

    //if((gap_to_fit_rotate==12 || gap_to_fit_rotate==6 || selected_TOF2==12 || selected_TOF2==6) && Rotate==1){
    if(gap_to_fit==6 || gap_to_fit==12){
      gr3_Leptons_rotate->Draw("P");
      best_fit_rotate->Draw("same");
    }
    else{
      gr_lepton_3->Draw("P");
    }

    if(vec_xx_kaon.size() >0){
      if (vertical_kaon) {
        gr_kaon->Draw("P");
        TLine *kaon_fit = new TLine(vec_xx_kaon[0], 50., vec_xx_kaon[0], -50.);
        kaon_fit->SetLineWidth(2);
        kaon_fit->SetLineColor(4);
        kaon_fit->Draw();
      }
      else {
        if(kaon_bk) gr_kaon_bk->Draw("P");
        else gr_kaon->Draw("P");
      }
    }

    gr_TOF1->Draw("P");
    gr_TOF1_closest->Draw("P");

    gr_int_TDC_TARGET->Draw("P");
    gr_int_TDC_Gap_Fibers_SFT->Draw("P");
    //gr_int_TDC_Gap_Fibers->Draw("P");
    gr_kaon_stop->Draw("P");

    tex_Angle_guide->Draw("same");
    tex_ChiS->Draw("same");
    tex_Kstop_X->Draw("same");
    tex_Kstop_Y->Draw("same");
    x_sft->Draw("same");
    y_sft->Draw("same");
    x_target->Draw("same");
    y_target->Draw("same");
}

void Two_tracks_Subplot_1(TCanvas *canvas,
  Int_t Run_Number, Int_t ievt, TGraph *gr3_Leptons, Int_t gap_to_fit)

{

  char Event_Display[200];
  // Add in Entry Number
  // sprintf(Event_Display,"Lepton Fit - Run Number: %d, Event Number: %d, Entry Number: %d",
  // Run_Number, ievt);
  sprintf(Event_Display,"Lepton Fit-Run #: %d, Event #: %d",
  Run_Number, ievt);

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


  // TGraph *delta_bar = new TGraph(1);
  // delta_bar->SetPoint(0, delta_x*10, delta_y*10);
  // delta_bar->SetMarkerStyle(21);
  // delta_bar->SetMarkerColor(30);
  // delta_bar->SetMarkerSize(0.8);






  // TArrow *delta_u_arrow = new TArrow(delta_x*10,delta_y*10,delta_stop_x*10,delta_stop_y*10,0.05,"|>");
  // //  ar2->SetAngle(40);
  //  delta_u_arrow->SetLineWidth(2);



  // TH2F *k_stop_point = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  // k_stop_point->Fill(dv_x*10, dv_y*10);
  // k_stop_point->SetMarkerStyle(2);
  // k_stop_point->SetMarkerColor(4);
  // k_stop_point->SetMarkerSize(2.5);



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

 Gaps[gap_to_fit-1]->SetLineColor(42);
    // if (tof1N[i] != -999) Gaps[tof1N[i]]->SetLineColor(42);
  //Draw Plot 1

  canvas->cd(1)->Range(-50, -50, 50, 50);
  h_EMPTY->Draw();
  // best_fit_distance->Draw("same");
  //gr3_Leptons_test->Draw("P");
  A1->Draw();
  A2->Draw();

  h_Circle->Draw("same");

  // hFLIP->Draw("same");
  //if(vec_xx_kaon.size()>0) gr_kaon->Draw("P");

  Gap1l->Draw("same");
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
  Gap12l->Draw("same");

  ell->Draw("same");

  ell_Target->Draw("same");

  ell_L1->Draw("same");

  h_Target_Center->Draw("same");

  //h_TOF1->Draw("same");
  //h_TOF1_closest->Draw("same");

  gr3_Leptons->Draw("P");
  // gr3_Leptons_low_energy->Draw("P");
  // delta_bar->Draw("P");
  // delta_u_arrow->Draw("same");


}


void Two_tracks_Subplot_2 (TCanvas *canvas, TGraph *gr_kaon, TGraph *gr_kaon_bk,
vector<double> vec_xx_kaon, bool kaon_bk,
bool vertical_kaon)
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

  // const int gap_to_fit = irt;
  // const int gap_to_fit_rotate = gap_to_fit;
  // const int selected_TOF2 = 0;

  canvas->cd(2)->Range(-50, -50, 50, 50);

  h_EMPTY->Draw(); //x

  A1->Draw(); //x
  A2->Draw(); //x

  h_Circle->Draw("same"); //x
  if(vec_xx_kaon.size() >0){
    if (vertical_kaon) {
      gr_kaon->Draw("P");
      TLine *kaon_fit = new TLine(vec_xx_kaon[0], 50., vec_xx_kaon[0], -50.);
      kaon_fit->SetLineWidth(2);
      kaon_fit->SetLineColor(4);
      kaon_fit->Draw();
    }
    else {
      if(kaon_bk) gr_kaon_bk->Draw("P");
      else gr_kaon->Draw("P");
    }
  }
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

  ell->Draw("same");
  ell_Target->Draw("same");

  ell_L1->Draw("same");

  h_Target_Center->Draw("same"); //x
}

void Two_tracks_Subplot_3(TCanvas *canvas, TGraph *gr3_Leptons,
  TGraph *gr_kaon, TGraph *gr_kaon_bk,
  vector<double> vec_xx_kaon, bool kaon_bk,
  bool vertical_kaon, TGraph *gr_kaon_stop)
{
  canvas->cd(3);
  char plot_title[100];
  sprintf(plot_title, "Leptons and Kaons");
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

  // TH2F *k_stop_point = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  // k_stop_point->Fill(dv_x*10, dv_y*10);
  // k_stop_point->SetMarkerStyle(2);
  // k_stop_point->SetMarkerColor(4);
  // k_stop_point->SetMarkerSize(4.5);

  // TH2F *k_stop_experimental = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  // k_stop_experimental->Fill(k_stop_x,k_stop_y);
  // k_stop_experimental->SetMarkerStyle(3);
  // k_stop_experimental->SetMarkerColor(6);
  // k_stop_experimental->SetMarkerSize(2.5);




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

  if(vec_xx_kaon.size() >0){
    if (vertical_kaon) {
      gr_kaon->Draw("P");
      TLine *kaon_fit = new TLine(vec_xx_kaon[0], 50., vec_xx_kaon[0], -50.);
      kaon_fit->SetLineWidth(2);
      kaon_fit->SetLineColor(4);
      kaon_fit->Draw();
    }
    else {
      if(kaon_bk) gr_kaon_bk->Draw("P");
      else gr_kaon->Draw("P");
    }
  }
  A1->Draw();
  A2->Draw();
  h_Circle->Draw("same");


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
  Gap12l->Draw("same");
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
  // k_stop_point->Draw("same");
  // k_stop_experimental->Draw("same");
  gr3_Leptons->Draw("P");
  gr_kaon_stop->Draw("P");

}


void Two_tracks_Subplot_4(TCanvas *canvas,
  vector<int> lepton_bars,
  TGraph *gr_TOF1, TGraph *gr_TOF1_closest,
  int Subplot, const char *plot_title,
  Color_t color) //Plot tracks of Leptons
{


  canvas->cd(Subplot);
  // char *p =;
  // sprintf(plot_title, title);
  TH2F *h_EMPTY = new TH2F(plot_title, plot_title, 500, -50, 50, 500, -50, 50);
  TGaxis *A1 = new TGaxis(-50,50,50,50,"pol1",510,"-U");
  TGaxis *A2 = new TGaxis(50,-50,50,50,"pol1",510,"+U");
  TH2F *h_Target_Center = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  h_Target_Center->Fill(0., 0.);
  // h_Target_Center->Fill(tof1wpos[0][0], tof1wpos[0][1]);
  h_Target_Center->SetMarkerStyle(5);
  h_Target_Center->SetMarkerColor(1);
  h_Target_Center->SetMarkerSize(2);

  vector<double> lepton_x_coords, lepton_y_coords;
  vector<double> lepton_ex, lepton_ey;
  for (unsigned int bar = 0; bar < lepton_bars.size(); ++bar) {
      lepton_x_coords.push_back(Xloc[lepton_bars[bar]]);
      lepton_y_coords.push_back(Yloc[lepton_bars[bar]]);
      lepton_ex.push_back(TARGET_Errors_X);
      lepton_ey.push_back(TARGET_Errors_Y);

    }

  TLine *best_fit_rotate;
  TGraphErrors *lepton_track;
  Int_t fitStatus;
  if (lepton_bars.size() > 0) {
    lepton_track = new TGraphErrors(lepton_x_coords.size(),&lepton_x_coords[0],&lepton_y_coords[0]
                                            , &lepton_ex[0], &lepton_ey[0]);
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

  A1->Draw(); //x
  A2->Draw(); //x

  h_Circle->Draw("same"); //x
  if(lepton_bars.size()>0) lepton_track->Draw("P"); //x
  if (fitStatus != 0) best_fit_rotate->Draw("same");
  //best_fit_rotate->Draw("same");

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

  ell->Draw("same");
  ell_Target->Draw("same");

  ell_L1->Draw("same");

  h_Target_Center->Draw("same"); //x
  gr_TOF1->Draw("P");
  gr_TOF1_closest->Draw("P");

}

void Two_tracks_Subplot_6(TCanvas *canvas, vector<int> lepton_bars_l1, vector<int> lepton_bars_l2,
  TGraph *gr_TOF1, TGraph *gr_TOF1_closest,
  TGraph *gr_kaon_stop,
  const char *plot_title, int Subplot,
  Color_t color_1, Color_t color_2) //Plot tracks of Leptons
{
  canvas->cd(Subplot);
  TH2F *h_EMPTY = new TH2F(plot_title, plot_title, 500, -50, 50, 500, -50, 50);
  TGaxis *A1 = new TGaxis(-50,50,50,50,"pol1",510,"-U");
  TGaxis *A2 = new TGaxis(50,-50,50,50,"pol1",510,"+U");
  TH2F *h_Target_Center = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);
  h_Target_Center->Fill(0., 0.);
  // h_Target_Center->Fill(tof1wpos[0][0], tof1wpos[0][1]);
  h_Target_Center->SetMarkerStyle(5);
  h_Target_Center->SetMarkerColor(1);
  h_Target_Center->SetMarkerSize(2);

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
  return;

  TH2F *h_Circle = new TH2F("Test", "TARGET and TOF1", 500, -50, 50, 500, -50, 50);
  double Gap[12][3][2] = {{{0}}};


  h_Circle->SetMarkerStyle(5);
  h_Circle->SetMarkerSize(1.2);
  h_Circle->SetLineWidth(2);

  TEllipse *ell = new TEllipse(0, 0, R_TOF1, 0);
  TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0);
  TEllipse *ell_L1 = new TEllipse(0, 0, R_SFT_L1, 0);
  ell->SetFillStyle(0);                 ell_Target->SetFillStyle(0);                 ell_L1->SetFillStyle(0);
  ell->SetLineColor(6);                 ell_Target->SetLineColor(1);                 ell_L1->SetLineColor(4);
  ell->SetLineWidth(1);                 ell_Target->SetLineWidth(1);                 ell_L1->SetLineWidth(1);


  canvas->cd(Subplot)->Range(-50, -50, 50, 50);

  h_EMPTY->Draw(); //x

  A1->Draw(); //x
  A2->Draw(); //x

  h_Circle->Draw("same"); //x
  //gr3_Leptons_rotate->Draw("P"); //x
  if(lepton_bars_l1.size()>0) lepton_track_1->Draw("P"); //x
  if(lepton_bars_l2.size()>0) lepton_track_2->Draw("P"); //x
  if (fitStatus1 != 0) best_fit_rotate1->Draw("same");
  if (fitStatus2 != 0) best_fit_rotate2->Draw("same");
  //best_fit_rotate->Draw("same");

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

  ell->Draw("same");
  ell_Target->Draw("same");

  ell_L1->Draw("same");

  h_Target_Center->Draw("same"); //x
  gr_TOF1->Draw("P");
  gr_TOF1_closest->Draw("P");
  gr_kaon_stop->Draw("P");







}


void TEST_subplot_1(TH2F *h_Lepton_fit, TGaxis *A1,
  TGaxis *A2, TGraph *gr_Target_Center,
  TEllipse *ell, TEllipse *ell_Target,
  TEllipse *ell_L1,
  TGraph *gr_lepton) {
    h_Lepton_fit->Draw();
    A1->Draw();
    A2->Draw();
    gr_Target_Center->Draw("P");

    //gr_TOF1_Markers->Draw("P");

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

    ell->Draw("same");
    ell_Target->Draw("same");
    ell_L1->Draw("same");

    gr_lepton->Draw("P");
}

void TEST_subplot_2(TH2F *h_kaon_fit, TGaxis *A1,
  TGaxis *A2, TGraph *gr_Target_Center,
  TEllipse *ell, TEllipse *ell_Target,
  TEllipse *ell_L1,
  int vec_xx_kaon_size, bool kaon_bk,
  TGraph *gr_kaon_bk, TGraph *gr_kaon,
  vector<double> vec_xx_kaon, bool vertical_kaon) {
    h_kaon_fit->Draw();
    A1->Draw();
    A2->Draw();
    gr_Target_Center->Draw("P");

    //gr_TOF1_Markers->Draw("P");

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

    ell->Draw("same");
    ell_Target->Draw("same");
    ell_L1->Draw("same");

    //if(vec_xx_kaon_size>0){
    //  if(kaon_bk) gr_kaon_bk->Draw("P");
    //  else gr_kaon->Draw("P");
    //}

    if(vec_xx_kaon.size() >0){
      if (vertical_kaon) {
        gr_kaon->Draw("P");
        TLine *kaon_fit = new TLine(vec_xx_kaon[0], 50., vec_xx_kaon[0], -50.);
        kaon_fit->SetLineWidth(2);
        kaon_fit->SetLineColor(4);
        kaon_fit->Draw();
      }
      else {
        if(kaon_bk) gr_kaon_bk->Draw("P");
        else gr_kaon->Draw("P");
      }
    }



}

/*void Draw_MWPC(TCanvas *cMWPC,
  TH2F *C2_hist, TH2F *C3_hist, TH2F *C4_hist,
  TLatex *tex_Title_MWPC_C2,
  TLatex *tex_Title_MWPC_C3,
  TLatex *tex_Title_MWPC_C4,
  TLatex *tex_palette_MWPC_C2[10],
  TLatex *tex_palette_MWPC_C3[10],
  TLatex *tex_palette_MWPC_C4[10],
  TMarker *palette_MWPC_C2[10],
  TMarker *palette_MWPC_C3[10],
  TMarker *palette_MWPC_C4[10],
  TGaxis *C2_top, TGaxis *C2_right,
  TGaxis *C3_top, TGaxis *C3_right,
  TGaxis *C4_top, TGaxis *C4_right,
  int C2X_count, int C2Y_count,
  int C3X_count, int C3Y_count,
  int C4X_count, int C4Y_count,
  TLine *C2X_line[C2X_count],
  TLine *C2Y_line[C2Y_count],
  TLine *C3X_line[C3X_count],
  TLine *C3Y_line[C3Y_count],
  TLine *C4X_line[C4X_count],
  TLine *C4Y_line[C4Y_count]
  ){
  cMWPC->cd(1);
  cMWPC->cd(1)->Update();
  cMWPC->cd(1)->Range(-280,-1,280,16);

  C2_hist->Draw();

  tex_Title_MWPC_C2->Draw();

  for(Int_t ipal=0; ipal<10; ipal++)  palette_MWPC_C2[ipal]->Draw();
  for(Int_t ileg=0; ileg<10; ileg++)  tex_palette_MWPC_C2[ileg]->Draw();

  C2_top->Draw("same");
  C2_right->Draw("same");

  for(int i = 0; i < C2X_count; i++) C2X_line[i]->Draw();
  for(int i = 0; i < C2Y_count; i++) C2Y_line[i]->Draw();

  cMWPC->cd(2);
  cMWPC->cd(2)->Update();
  cMWPC->cd(2)->Range(-1,-1,64,16);

  C3_hist->Draw();
  tex_Title_MWPC_C3->Draw();

  for(Int_t ipal=0; ipal<10; ipal++)  palette_MWPC_C3[ipal]->Draw();
  for(Int_t ileg=0; ileg<10; ileg++)  tex_palette_MWPC_C3[ileg]->Draw();

  C3_top->Draw("same");
  C3_right->Draw("same");

  for(int i = 0; i < C3X_count; i++) C3X_line[i]->Draw();
  for(int i = 0; i < C3Y_count; i++) C3Y_line[i]->Draw();

  cMWPC->cd(3);
  cMWPC->cd(3)->Update();
  cMWPC->cd(3)->Range(-1,-1,72,16);

  C4_hist->Draw();
  tex_Title_MWPC_C4->Draw();

  for(Int_t ipal=0; ipal<10; ipal++)  palette_MWPC_C4[ipal]->Draw();
  for(Int_t ileg=0; ileg<10; ileg++)  tex_palette_MWPC_C4[ileg]->Draw();

  C4_top->Draw("same");
  C4_right->Draw("same");

  for(int i = 0; i < C4X_count; i++) C4X_line[i]->Draw();
  for(int i = 0; i < C4Y_count; i++) C4Y_line[i]->Draw();

}
*/
void SFT_subplot_1(Int_t Run_Number, Int_t ievt,
  Bool_t Event_flag[40],
  TLine *hline_DS[8], TLine *vline_DS[68],
  TLine *hline_US[8], TLine *vline_US[68],
  Double_t ADC_High_SFT_corr[128], Int_t has_TDC_SFT_hit[128],
  Int_t ADC_High_SFT[128],
  Int_t ADC_Low_SFT[128],
  TMarker *marker_DS[64], TMarker *marker_US[64],
  Int_t par_temp[2][128], int Switch) {
  for(Int_t hdraw_DS=0; hdraw_DS<8; hdraw_DS++) hline_DS[hdraw_DS]->Draw();
  for(Int_t hdraw_US=0; hdraw_US<8; hdraw_US++) hline_US[hdraw_US]->Draw();
  for(Int_t vdraw_DS=0; vdraw_DS<68; vdraw_DS++)  vline_DS[vdraw_DS]->Draw();
  for(Int_t vdraw_US=0; vdraw_US<68; vdraw_US++)  vline_US[vdraw_US]->Draw();

  TLatex *tex_DS_SFT;
  TLatex *tex_footer_SFT;
  TLatex *tex_version;
  TLatex *tex_Event_flag;

  tex_DS_SFT = new TLatex(0.02,0.88,"SFT  --  Downstream & Upstream");
  tex_DS_SFT->SetTextSize(0.05);
  tex_DS_SFT->SetLineWidth(2);

  char Event_flag_title[100];
  sprintf(Event_flag_title,"Event Flag = %d", Event_flag[0]);
  tex_Event_flag = new TLatex(0.02,0.2,Event_flag_title);
  tex_Event_flag->SetTextSize(0.04);
  tex_Event_flag->SetLineWidth(2);

  char footer[100];
  sprintf(footer,"Event_Display.C  --  Run %d ; Event %d",Run_Number,ievt);
  tex_footer_SFT = new TLatex(0.02,0.08,footer);
  tex_footer_SFT->SetTextSize(0.04);
  tex_footer_SFT->SetLineWidth(2);

  char Version[100] = "Version 5.3";
  char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!
  tex_version = new TLatex(0.02,0.02,Version);
  tex_version->SetTextSize(0.04);
  tex_version->SetLineWidth(2);

  tex_DS_SFT->Draw();
  tex_Event_flag->Draw();
  tex_footer_SFT->Draw();
  tex_version->Draw();

  TLatex *tex_number1_SFT;
  TLatex *tex_number2_SFT;
  TLatex *tex_number3_SFT;
  TLatex *tex_number4_SFT;
  TLatex *tex_number5_SFT;
  TLatex *tex_number6_SFT;
  TLatex *tex_number7_SFT;
  TLatex *tex_number8_SFT;
  TLatex *tex_number9_SFT;
  TLatex *tex_number10_SFT;
  TLatex *tex_number11_SFT;
  TLatex *tex_number12_SFT;
  TLatex *tex_number13_SFT;
  TLatex *tex_number14_SFT;
  TLatex *tex_number15_SFT;
  TLatex *tex_number16_SFT;
  TLatex *tex_number17_SFT;

  TLatex *tex_number21_SFT; TLatex *tex_number21_bis_SFT;
  TLatex *tex_number22_SFT; TLatex *tex_number22_bis_SFT;
  TLatex *tex_number23_SFT; TLatex *tex_number23_bis_SFT;
  TLatex *tex_number24_SFT; TLatex *tex_number24_bis_SFT;

  tex_number1_SFT = new TLatex(0.115,0.81,"1");
  tex_number1_SFT->SetTextSize(0.025);
  tex_number1_SFT->SetLineWidth(2);

  tex_number2_SFT = new TLatex(0.15,0.81,"2");
  tex_number2_SFT->SetTextSize(0.025);
  tex_number2_SFT->SetLineWidth(2);

  tex_number3_SFT = new TLatex(0.19,0.81,"3");
  tex_number3_SFT->SetTextSize(0.025);
  tex_number3_SFT->SetLineWidth(2);

  tex_number4_SFT = new TLatex(0.23,0.81,"4");
  tex_number4_SFT->SetTextSize(0.025);
  tex_number4_SFT->SetLineWidth(2);

  tex_number5_SFT = new TLatex(0.27,0.81,"5");
  tex_number5_SFT->SetTextSize(0.025);
  tex_number5_SFT->SetLineWidth(2);

  tex_number6_SFT = new TLatex(0.31,0.81,"6");
  tex_number6_SFT->SetTextSize(0.025);
  tex_number6_SFT->SetLineWidth(2);

  tex_number7_SFT = new TLatex(0.35,0.81,"7");
  tex_number7_SFT->SetTextSize(0.025);
  tex_number7_SFT->SetLineWidth(2);

  tex_number8_SFT = new TLatex(0.39,0.81,"8");
  tex_number8_SFT->SetTextSize(0.025);
  tex_number8_SFT->SetLineWidth(2);

  tex_number9_SFT = new TLatex(0.43,0.81,"9");
  tex_number9_SFT->SetTextSize(0.025);
  tex_number9_SFT->SetLineWidth(2);

  tex_number10_SFT = new TLatex(0.465,0.81,"10");
  tex_number10_SFT->SetTextSize(0.025);
  tex_number10_SFT->SetLineWidth(2);

  tex_number11_SFT = new TLatex(0.505,0.81,"11");
  tex_number11_SFT->SetTextSize(0.025);
  tex_number11_SFT->SetLineWidth(2);

  tex_number12_SFT = new TLatex(0.545,0.81,"12");
  tex_number12_SFT->SetTextSize(0.025);
  tex_number12_SFT->SetLineWidth(2);

  tex_number13_SFT = new TLatex(0.585,0.81,"13");
  tex_number13_SFT->SetTextSize(0.025);
  tex_number13_SFT->SetLineWidth(2);

  tex_number14_SFT = new TLatex(0.625,0.81,"14");
  tex_number14_SFT->SetTextSize(0.025);
  tex_number14_SFT->SetLineWidth(2);

  tex_number15_SFT = new TLatex(0.665,0.81,"15");
  tex_number15_SFT->SetTextSize(0.025);
  tex_number15_SFT->SetLineWidth(2);

  tex_number16_SFT = new TLatex(0.705,0.81,"16");
  tex_number16_SFT->SetTextSize(0.025);
  tex_number16_SFT->SetLineWidth(2);

  tex_number17_SFT = new TLatex(0.745,0.81,"17");
  tex_number17_SFT->SetTextSize(0.025);
  tex_number17_SFT->SetLineWidth(2);

  tex_number21_SFT = new TLatex(0.03,0.77,"L1-D");  tex_number21_bis_SFT = new TLatex(0.03,0.73,"L1-U");
  tex_number21_SFT->SetTextSize(0.03);      tex_number21_bis_SFT->SetTextSize(0.03);
  tex_number21_SFT->SetLineWidth(2);        tex_number21_bis_SFT->SetLineWidth(2);

  tex_number22_SFT = new TLatex(0.03,0.63,"L2-D");  tex_number22_bis_SFT = new TLatex(0.03,0.59,"L2-U");
  tex_number22_SFT->SetTextSize(0.03);      tex_number22_bis_SFT->SetTextSize(0.03);
  tex_number22_SFT->SetLineWidth(2);        tex_number22_bis_SFT->SetLineWidth(2);

  tex_number23_SFT = new TLatex(0.03,0.49,"L3-D");  tex_number23_bis_SFT = new TLatex(0.03,0.45,"L3-U");
  tex_number23_SFT->SetTextSize(0.03);      tex_number23_bis_SFT->SetTextSize(0.03);
  tex_number23_SFT->SetLineWidth(2);        tex_number23_bis_SFT->SetLineWidth(2);

  tex_number24_SFT = new TLatex(0.03,0.35,"L4-D");  tex_number24_bis_SFT = new TLatex(0.03,0.31,"L4-U");
  tex_number24_SFT->SetTextSize(0.03);      tex_number24_bis_SFT->SetTextSize(0.03);
  tex_number24_SFT->SetLineWidth(2);        tex_number24_bis_SFT->SetLineWidth(2);


  tex_number1_SFT->Draw();  tex_number10_SFT->Draw();
  tex_number2_SFT->Draw();  tex_number11_SFT->Draw();
  tex_number3_SFT->Draw();  tex_number12_SFT->Draw();
  tex_number4_SFT->Draw();  tex_number13_SFT->Draw();
  tex_number5_SFT->Draw();  tex_number14_SFT->Draw();
  tex_number6_SFT->Draw();  tex_number15_SFT->Draw();
  tex_number7_SFT->Draw();  tex_number16_SFT->Draw();
  tex_number8_SFT->Draw();  tex_number17_SFT->Draw();
  tex_number9_SFT->Draw();

  tex_number21_SFT->Draw();   tex_number21_bis_SFT->Draw();
  tex_number22_SFT->Draw();   tex_number22_bis_SFT->Draw();
  tex_number23_SFT->Draw();   tex_number23_bis_SFT->Draw();
  tex_number24_SFT->Draw();   tex_number24_bis_SFT->Draw();

  for(Int_t imark1=0; imark1<15; imark1++){
    if((ADC_High_SFT_corr[par_temp[1][imark1]]!=0) && (has_TDC_SFT_hit[par_temp[1][imark1]])){
      marker_DS[par_temp[0][imark1]]->SetMarkerColor(2);
      marker_DS[par_temp[0][imark1]]->Draw();
    }

    if((ADC_High_SFT_corr[par_temp[1][imark1+64]]!=0) && (has_TDC_SFT_hit[par_temp[1][imark1+64]])){
      marker_US[par_temp[0][imark1+64]]->SetMarkerColor(2);
      marker_US[par_temp[0][imark1+64]]->Draw();
    }

    if(ADC_High_SFT[imark1]<0 && ADC_Low_SFT[imark1]>0 && Switch==1){
      marker_DS[par_temp[0][imark1]]->SetMarkerColor(kGray);
      marker_DS[par_temp[0][imark1]]->Draw();
    }

    if(ADC_High_SFT[imark1+64]<0 && ADC_Low_SFT[imark1+64]>0 && Switch==1){
      marker_US[par_temp[0][imark1+64]]->SetMarkerColor(kGray);
      marker_US[par_temp[0][imark1+64]]->Draw();
    }
  }

  for(Int_t imark2=0; imark2<15; imark2++){
    if((ADC_High_SFT_corr[par_temp[1][imark2+15]]!=0) && (has_TDC_SFT_hit[par_temp[1][imark2+15]])){
      marker_DS[par_temp[0][imark2+15]]->SetMarkerColor(4);
      marker_DS[par_temp[0][imark2+15]]->Draw();
    }

    if((ADC_High_SFT_corr[par_temp[1][imark2+79]]!=0) && (has_TDC_SFT_hit[par_temp[1][imark2+79]])){
      marker_US[par_temp[0][imark2+79]]->SetMarkerColor(4);
      marker_US[par_temp[0][imark2+79]]->Draw();
    }

    if(ADC_High_SFT[imark2+15]<0 && ADC_Low_SFT[imark2+15]>0 && Switch==1){
      marker_DS[par_temp[0][imark2+15]]->SetMarkerColor(kGray);
      marker_DS[par_temp[0][imark2+15]]->Draw();
    }

    if(ADC_High_SFT[imark2+79]<0 && ADC_Low_SFT[imark2+79]>0 && Switch==1){
      marker_US[par_temp[0][imark2+79]]->SetMarkerColor(kGray);
      marker_US[par_temp[0][imark2+79]]->Draw();
    }
  }

  for(Int_t imark3=0; imark3<17; imark3++){
    if((ADC_High_SFT_corr[par_temp[1][imark3+30]]!=0) && (has_TDC_SFT_hit[par_temp[1][imark3+30]])){
      marker_DS[par_temp[0][imark3+30]]->SetMarkerColor(3);
      marker_DS[par_temp[0][imark3+30]]->Draw();
    }

    if((ADC_High_SFT_corr[par_temp[1][imark3+94]]!=0) && (has_TDC_SFT_hit[par_temp[1][imark3+94]])){
      marker_US[par_temp[0][imark3+94]]->SetMarkerColor(3);
      marker_US[par_temp[0][imark3+94]]->Draw();
    }

    if(ADC_High_SFT[imark3+30]<0 && ADC_Low_SFT[imark3+30]>0 && Switch==1){
      marker_DS[par_temp[0][imark3+30]]->SetMarkerColor(kGray);
      marker_DS[par_temp[0][imark3+30]]->Draw();
    }

    if(ADC_High_SFT[imark3+94]<0 && ADC_Low_SFT[imark3+94]>0 && Switch==1){
      marker_US[par_temp[0][imark3+94]]->SetMarkerColor(kGray);
      marker_US[par_temp[0][imark3+94]]->Draw();
    }
  }

  for(Int_t imark4=0; imark4<17; imark4++){
    if((ADC_High_SFT_corr[par_temp[1][imark4+47]]!=0) && (has_TDC_SFT_hit[par_temp[1][imark4+47]])){
      marker_DS[par_temp[0][imark4+47]]->SetMarkerColor(1);
      marker_DS[par_temp[0][imark4+47]]->Draw();
    }

    if((ADC_High_SFT_corr[par_temp[1][imark4+111]]!=0) && (has_TDC_SFT_hit[par_temp[1][imark4+111]])){
      marker_US[par_temp[0][imark4+111]]->SetMarkerColor(1);
      marker_US[par_temp[0][imark4+111]]->Draw();
    }

    if(ADC_High_SFT[imark4+47]<0 && ADC_Low_SFT[imark4+47]>0 && Switch==1){
      marker_DS[par_temp[0][imark4+47]]->SetMarkerColor(kGray);
      marker_DS[par_temp[0][imark4+47]]->Draw();
    }

    if(ADC_High_SFT[imark4+111]<0 && ADC_Low_SFT[imark4+111]>0 && Switch==1){
      marker_US[par_temp[0][imark4+111]]->SetMarkerColor(kGray);
      marker_US[par_temp[0][imark4+111]]->Draw();
    }
  }
}

void SFT_subplot_2(TH1F *h_ADC_L1_DS, TH1F *h_ADC_L2_DS,
TH1F *h_ADC_L3_DS, TH1F *h_ADC_L4_DS,
TH1F *h_ADC_L1_US, TH1F *h_ADC_L2_US,
TH1F *h_ADC_L3_US, TH1F *h_ADC_L4_US, Int_t ADC_High_corr_max) {
  h_ADC_L1_DS->SetMaximum(1.1*ADC_High_corr_max);
  h_ADC_L1_DS->Draw();          h_ADC_L1_US->Draw("same");
  h_ADC_L2_DS->Draw("same");    h_ADC_L2_US->Draw("same");
  h_ADC_L3_DS->Draw("same");    h_ADC_L3_US->Draw("same");
  h_ADC_L4_DS->Draw("same");    h_ADC_L4_US->Draw("same");

  // Split LineADC_High_corr_max
  TLine *split_line;
  split_line = new TLine(64,0,64,1.1*ADC_High_corr_max);
  split_line->SetLineWidth(4);
  split_line->SetLineColor(4);
  split_line->SetLineStyle(2);
  split_line->Draw();
}
