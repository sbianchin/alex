TH2F *hXY; // position - satisfies E and T cuts
TH2F *hXYTOF1; // TOF1 position
TH2F *hXYlowT; // position - satisfies E and T cuts
//TH2F *hXYTOF1time; // the valid gap TOF1 position { empty circle}
TH2F *hTvsE; // time vs energy
TH2F *hTvsElowT; // time vs energy - satisfy T cut
//TH2F *hTvsETOF1; // time vs energy of highest above time hits peak
TCanvas* c1;
//TCanvas* c2;
//TF1 *fitfunc1;
//TLine *fitfunc2;
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t MakeFindTrackHistos(){
    hXY = new TH2F("YvsX", "Y vs X ",100,-50,50,100,-50,50);
    hXY->SetMarkerStyle(22);
    hXY->SetMarkerColor(1);  //Black
    hXYTOF1 = new TH2F("YvsXTOF1", "Y vs X TOF1",100,-50,50,100,-50,50);
    hXYTOF1->SetMarkerStyle(22);
    hXYTOF1->SetMarkerColor(6);  //Mauve
    hXYlowT = new TH2F("YvsXlowT", "Y vs X low T",100,-50,50,100,-50,50);
    hXYlowT->SetMarkerStyle(23);
    hXYlowT->SetMarkerColor(2);  //Red
    //hTOF1 = new TH2F("TOF1", "TOF1",100,-50,50,100,-50,50);
    //hTOF1->SetMarkerStyle(4);
    // hTvsE = new TH2F("TLvsE", "T vs E ",300,0.0,10.0,300,0.0,75.0);
    hTvsE = new TH2F("TLvsE", "T vs E ",300,0.0,10.0,300, 0.0, 40);
    hTvsE->SetMarkerStyle(22);
    hTvsE->SetMarkerColor(1);  //Black
    // hTvsElowT = new TH2F("TLvsElowT", "T vs E Low T",300,0.0,10.0,300,0.0,75.0);
    hTvsElowT = new TH2F("TLvsElowT", "T vs E Low T",300,0.0,10.0,300,0.0, 40);
    hTvsElowT->SetMarkerStyle(23);
    hTvsElowT->SetMarkerColor(2);  //Red
    return 0;
}
