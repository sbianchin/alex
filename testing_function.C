void testing_function(){
  TH2F *h_lepton_fit_final = new TH2F("Lepton Fit", "str_lepton_fit_final", 500, -50, 50, 500, -50, 50);

  // TEllipse *ell = new TEllipse(0, 0, R_TOF1, 0);
  // TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0);
  // TEllipse *ell_L1 = new TEllipse(0, 0, R_SFT_L1, 0);
  //
  // TGaxis *A1 = new TGaxis(-50,50,50,50,"pol1",510,"-U");
  // TGaxis *A2 = new TGaxis(50,-50,50,50,"pol1",510,"+U");
  //
  // vector<double> vec_xx_Target_Center;  vec_xx_Target_Center.clear();
  // vector<double> vec_yy_Target_Center;  vec_yy_Target_Center.clear();
  // vec_xx_Target_Center.push_back(0.);
  // vec_yy_Target_Center.push_back(0.);
  //
  // TGraph *gr_Target_Center = new TGraph(vec_xx_Target_Center.size(),&vec_xx_Target_Center[0],&vec_yy_Target_Center[0]);
  // gr_Target_Center->SetMarkerStyle(5);
  // gr_Target_Center->SetMarkerColor(1);
  // gr_Target_Center->SetMarkerSize(1);

  TCanvas *a_test;
  a_test = new TCanvas("My test","My test",0,100,700,350);
  a_test->Divide(1,1);
  a_test->cd(1);

  h_lepton_fit_final->Draw();
  // A1->Draw();
  // A2->Draw();
  // gr_Target_Center->Draw("P");
  //
  // //gr_TOF1_Markers->Draw("P");
  //
  // Gap1l->Draw();
  // Gap2l->Draw("same");
  // Gap3l->Draw("same");
  // Gap4l->Draw("same");
  // Gap5l->Draw("same");
  // Gap6l->Draw("same");
  // Gap7l->Draw("same");
  // Gap8l->Draw("same");
  // Gap9l->Draw("same");
  // Gap10l->Draw("same");
  // Gap11l->Draw("same");
  // Gap12l->Draw("Xsame");
  //
  // TGraph *gr_TOF1 = new TGraph(vec_xx_TOF1.size(),&vec_xx_TOF1[0],&vec_yy_TOF1[0]);
  // gr_TOF1->SetMarkerStyle(20);
  // gr_TOF1->SetMarkerColor(2);
  // gr_TOF1->SetMarkerSize(1.5);
  //
  //
  // TGraph *gr_TOF1_closest = new TGraph(vec_xx_TOF1_closest.size(),&vec_xx_TOF1_closest[0],&vec_yy_TOF1_closest[0]);
  // gr_TOF1_closest->SetMarkerStyle(20);
  // gr_TOF1_closest->SetMarkerColor(3);
  // gr_TOF1_closest->SetMarkerSize(1.5);
  //
  // ell->Draw("same");
  // ell_Target->Draw("same");
  // ell_L1->Draw("same");
  // gr_TOF1->Draw("P");
  // gr_TOF1_closest->Draw("P");
  double x[] = {0, 1, 2, 3, 4};
  double y[] = {0, 2, 4, 1, 3};
  double ex[] = {0.1, 0.2, 0.3, 0.4, 0.5};
  double ey[] = {1, 0.5, 1, 0.5, 1};
  // TGraphErrors *ge = new TGraphErrors(5, x, y, ex, ey);
  // ge->Draw("same");
  TGraphErrors ge = TGraphErrors(5, x, y, ex, ey);
  ge.Draw("same");
  // ge.Paint();

  // // gr_lepton_3.Draw("psame");
  // a_test->Update();
  // gPad->Update();
  gPad->Modified();
  gPad->Update();
}
