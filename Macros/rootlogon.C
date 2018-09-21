{
//TColor *color;
gROOT->SetStyle("Plain");
//Int_t colors[]={51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76};
//Int_t colory[26];
//Float_t r[]={0.0,0.0,0.0,0.0,0.0,0.4,0.6,1.0,1.0,1.0,1.0,1.0,1.0,0.8,0.6,0.4};
//Float_t g[]={0.4,0.6,0.8,1.0,1.0,1.0,1.0,1.0,1.0,0.8,0.6,0.4,0.2,0.0,0.0,0.0}; 
//Float_t b[]={1.0,1.0,1.0,0.8,0.4,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

//for(int i=0;i<16;i++){
//  gROOT->GetListOfColors()->Remove(gROOT->GetColor(colors[i]));
  //gROOT->GetListOfColors()->Remove(gROOT->GetColor(colors[i]+100));
  //gROOT->GetListOfColors()->Remove(gROOT->GetColor(colors[i]+150));
//  color = new TColor(colors[i],r[i],g[i],b[i]);
//}
   
//gStyle->SetPalette(16,colors);

gROOT->LoadMacro("ProcEvent.C");
gSystem->Load("libTree"); 
gSystem->Load("libGpad"); 

}
