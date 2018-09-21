
void histogram_xy(int Run_Number, int flag, int min_evt = 0, int max_evt = 0) {

  ofstream fout;


  char Name_finput[200];
  sprintf(Name_finput,"G4Run%d.root", Run_Number);
  cout << "File opened:  " << Name_finput << endl;
  TChain fChain("Kaon");
  fChain.Add(Name_finput);
  Int_t nentries = (Int_t)fChain.GetEntries();
  double dv_x, dv_y, u_x, u_y;
  fChain.SetBranchAddress("dv_x", &dv_x);
  fChain.SetBranchAddress("dv_y", &dv_y);
  fChain.SetBranchAddress("u_x",&u_x);
	fChain.SetBranchAddress("u_y",&u_y);

  TCanvas *c1 = new TCanvas("c1", "c1",300,300);
  gStyle->SetOptStat(0);
  c1->Divide(1,1);
  c1->cd(1);
  TH2F *h2 = new TH2F("h2","",100,-50,50,100,-50,50);
  for (int i = 0; i < nentries; ++i) {
    fChain.GetEntry(i);
    h2->Fill(dv_x*10,dv_y*10);
    TMarker* Kstp=new TMarker(u_x/3.1+11,u_y/3.1+10,24);
  	Kstp->Draw("same");



  }
  gStyle->SetPalette(1);
  // h2->GetZaxis()->SetRangeUser(0,6);
  char name[50];
  sprintf(name,"Run %d K-stop Histogram",Run_Number); //JIIIIIIII
	h2->SetTitle(name);
  h2->SetXTitle("X-axis (mm)");
	h2->SetYTitle("Y-axis (mm)");
	h2->Draw("colz");







}

// {
//    TCanvas *c1 = new TCanvas("c1", "c1",900,900);
//    gStyle->SetOptStat(0);
//
//    // Create the three pads
//    TPad *center_pad = new TPad("center_pad", "center_pad",0.0,0.0,0.6,0.6);
//    center_pad->Draw();
//
//    right_pad = new TPad("right_pad", "right_pad",0.55,0.0,1.0,0.6);
//    right_pad->Draw();
//
//    top_pad = new TPad("top_pad", "top_pad",0.0,0.55,0.6,1.0);
//    top_pad->Draw();
//
//    // Create, fill and project a 2D histogram.
//    TH2F *h2 = new TH2F("h2","",40,-4,4,40,-20,20);
//    Float_t px, py;
//    for (Int_t i = 0; i < 25000; i++) {
//       gRandom->Rannor(px,py);
//       h2->Fill(px,5*py);
//    }
//    TH1D * projh2X = h2->ProjectionX();
//    TH1D * projh2Y = h2->ProjectionY();
//
//    // Drawing
//    center_pad->cd();
//    gStyle->SetPalette(1);
//    h2->Draw("COL");
//
//    top_pad->cd();
//    projh2X->SetFillColor(kBlue+1);
//    projh2X->Draw("bar");
//
//    right_pad->cd();
//    projh2Y->SetFillColor(kBlue-2);
//    projh2Y->Draw("hbar");
//
//    c1->cd();
//    TLatex *t = new TLatex();
//    t->SetTextFont(42);
//    t->SetTextSize(0.02);
//    t->DrawLatex(0.6,0.88,"This example demonstrate how to display");
//    t->DrawLatex(0.6,0.85,"a histogram and its two projections.");
// }
