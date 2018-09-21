#include "TH1F.h"
#ifndef __CINT__

#include "TAttMarker.h"
#include "Buttons.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFrame.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TMatrix.h"
#include "TPad.h"
#include "TPostScript.h"
#include "TStopwatch.h" 
#include "TStyle.h"
#include "TSystem.h"
#include "TText.h"
#include "TUnixSystem.h"
//#include "RQ_OBJECT.h"

#endif
 

//enum EEventType {
//   kNoEvent       =  0,
//   kButton1Down   =  1, kButton2Down   =  2, kButton3Down   =  3, kKeyDown  =  4,
//   kButton1Up     = 11, kButton2Up     = 12, kButton3Up     = 13, kKeyUp    = 14,
//   kButton1Motion = 21, kButton2Motion = 22, kButton3Motion = 23, kKeyPress = 24,
//   kButton1Locate = 41, kButton2Locate = 42, kButton3Locate = 43,
//   kMouseMotion   = 51, kMouseEnter    = 52, kMouseLeave    = 53,
//   kButton1Double = 61, kButton2Double = 62, kButton3Double = 63
//};

//___________________________________________________________________

void ProcEvent(Int_t event, Int_t px, Int_t py, TObject *sel)
{
   //  print event type and current cursor position


   TCanvas *c = (TCanvas *) gTQSender;
   TPad *pad = (TPad *) c->GetSelectedPad();
   
    if(!pad) return;
    gROOT->SetEditHistograms(kFALSE);
    //printf("event=%d, px=%d, py=%d\n", event, px, py);
    //Float_t x = pad->AbsPixeltoX(px);
    //Float_t y = pad->AbsPixeltoY(py);
    //x = pad->PadtoX(x);
    //y = pad->PadtoY(y);
    //printf("x=%.3g, y=%.3g\n",x,y);
   if(event==kButton1Double){ 
     pad->Pop();
     //printf("%s %d\n",pad->GetName(),pad->GetNumber());
     pad->cd();
     TCanvas *c_blow;
     TIter next(pad->GetListOfPrimitives());
     if((TCanvas*)gROOT->GetListOfCanvases()->FindObject("c_blow")){
       c_blow=(TCanvas*)gROOT->GetListOfCanvases()->FindObject("c_blow");
       c_blow->Clear();
       //printf("show %d\n",c_blow->GetUniqueID());
       //c_blow->GetCanvasImp()->Show();
       //c_blow->RaiseWindow();
       //c_blow->GetCanvasImp()->RaiseWindow();
       //c_blow->Flush();
       //gVirtualX->RaiseWindow(c_blow->GetUniqueID());
       //c_blow->Show();

     }
     else{
 //      c_blow = new TCanvas("c_blow","blowup",750,0,743,525);
       c_blow = new TCanvas("c_blow","blowup",100,50,600,600);
       c_blow->SetFillColor(10); //white
       //c_blow->ToggleEventStatus();
       //c_blow->SetCrosshair();
       c_blow->Draw();
     }
     c_blow->cd();

     gROOT->SetSelectedPad((TPad*)c_blow);
     TPad *clone =(TPad*)pad->Clone();
     clone->SetPad(0,0,1,1);
     clone->Draw();
     clone->Modified();
     clone->Update();
     c_blow->Show();
  }
}
//___________________________________________________________________

