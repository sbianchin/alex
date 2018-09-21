//Some TOF1 & TOF2 Variables
//Control
//TFile *hfout;
int TOF1SUM[nGaps];
int TOF1MT[nGaps];
// Test numbers for TOF1 neighbours
int GLnum[nGaps] = {11,0,1,2,3,4,5,6,7,8,9,10};
int GHnum[nGaps] = {1,2,3,4,5,6,7,8,9,10,11,0};
//TOF1 tdc cuts
int TOF1UtdcLcut[nGaps], TOF1UtdcHcut[nGaps];
int TOF1DtdcLcut[nGaps], TOF1DtdcHcut[nGaps];
// TOF2 tdc cuts
int TOF2AItdcLcut[nGaps], TOF2AItdcHcut[nGaps];
int TOF2AOtdcLcut[nGaps], TOF2AOtdcHcut[nGaps];
int TOF2BItdcLcut[nGaps], TOF2BItdcHcut[nGaps];
int TOF2BOtdcLcut[nGaps], TOF2BOtdcHcut[nGaps];
// TOF1 ADC Sum cuts
//int TOF1SumLowCut[nGaps] = {240, 500, 100, 300, 400, 450, 550, 520, 350, 400, 520, 420};
int TOF1SumLowCut[nGaps] = {240, 400, 150, 300, 400, 450, 450, 350, 350, 400, 350, 400};
// TOF1 TDC means
int TOF1MeanT[nGaps];
int TOF1NewMT[nGaps];
int TOF1MTlowCut = 500;  //400;
int TOF1MThghCut = 750;  //1125;
//TOF1 positions and cuts
int TOF1PosMin = -50;
int TOF1PosMax =  50;
int TOF1Pos[nGaps];
int TOF1NewPos[nGaps];
int TOF1Diff[nGaps];
//int TOF1PosLcut[nGaps] = {-60,-90,-80,-70,-80,-60,-60,-100,-60,-70,-80,-80};
int TOF1PosCl[nGaps] ={ -14, -38, 14, -22, -89, -66, -65, -110, -71, -71, -96, -96};
// TOF2
int TOF2APos[nGaps];
int TOF2BPos[nGaps];
int TOF2AMeanT[nGaps];
int TOF2BMeanT[nGaps];
int TOF2MeanT[nGaps];
int TOF2ASUM[nGaps];
int TOF2BSUM[nGaps];
//TOF2 ADC cuts
int TOF2upper = 225;
int TOF2adcMin[nTOF2] = {200,200,300,350,500,500,500,500,550,450,500,500,450,500,500,500,
    500,500,500,500,500,450,400,500,500,500,500,500,500,500,500,500,500,500,500,500,500,
    500,500,500,500,450,500,500,500,500,500,500};
 //TOF2 ADC Sum Cuts
 int TOF2ASCut = 700;
 int TOF2BSCut = 700;
//TOF2A positions(Diff) and cuts
int TOF2APosUedge[nGaps] = {300,325,160,220,120,240,260,270,340,110,270,200};
//TOF2A positions and cuts
//int TOF2BPosLcut[nGaps] = {-100,-130,-130,-100,-40,-150,-100,-20,0,-70,-70,-120};
int TOF2BPosUedge[nGaps] = {70,145,100,60,60,120,150,-10,-165,-20,110,160};
// TOF2 TDC mean & cuts
int TOF2MTlowCut = 800;
int TOF2MThghCut = 2200;
// parameters & values for new tdc calculations
Double_t Yfu, Ynu;
Double_t Yfd, Ynd;
Double_t An, Tp;
Int_t TOF1D_TDCnewL[nGaps], TOF1U_TDCnewL[nGaps];
Int_t TOF1D_TDCnewC[nGaps], TOF1U_TDCnewC[nGaps];
Int_t TOF1D_TDCnewH[nGaps], TOF1U_TDCnewH[nGaps];
Int_t TOF1UNewTDC = -1;
Int_t TOF1DNewTDC = -1;
// Fit break points
// - NOTE these will change if POSs converted to mm
Int_t PosCut1D[nGaps] = {-9,-16,-28,-18,24,31,18,21,26,25,31,24};
Int_t PosCut1U[nGaps] = {-7,-14,-28,-16,24,31,23,21,27,29,34,25};
Int_t Pos1Dmax = 55;
Int_t Pos1Umin = -55;
int TOF1_Gap[nTOF1] = {1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,7,8,9,10,11,12};
int TOF1_UD[nTOF1] = {1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2};
int TOF2_Gap[nTOF2] = {1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,7,8,9,10,11,12,
    1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,7,8,9,10,11,12};
int TOF2_UD[nTOF2] = {1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,
    3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4};
//TOF angle limits
Double_t alimL[12] = {15,45,75,105,135,165,195,225,255,285,315,345};
Double_t alimH[12] = {45,75,105,135,165,195,225,255,285,315,345,15};
// TOF1 centers
float TOF_Xloc[12] = {22.25,38.538,44.5,38.538,22.25,0,-22.25,-38.538,-44.5,-38.538,-22.25,0};
float TOF_Yloc[12] = {38.538,22.25,0,-22.25,-38.538,-44.5,-38.538,-22.25,0,22.25,38.538,44.5};
bool TOF1x_Gap[12] = {false};
bool TOF1y_Gap[12] = {false};
// fit parameters
const Int_t n1par = 4;
const Int_t n2par = 2;
Double_t par1dL[nGaps][n1par];
Double_t par1dL_err[nGaps][n1par];
Double_t par1uL[nGaps][n1par];
Double_t par1uL_err[nGaps][n1par];
Double_t par1dC[nGaps][n1par];
Double_t par1dC_err[nGaps][n1par];
Double_t par1uC[nGaps][n1par];
Double_t par1uC_err[nGaps][n1par];
Double_t par1dH[nGaps][n1par];
Double_t par1dH_err[nGaps][n1par];
Double_t par1uH[nGaps][n1par];
Double_t par1uH_err[nGaps][n1par];
Double_t par2AI[nGaps][n2par];
Double_t par2AI_err[nGaps][n2par];
Double_t par2AO[nGaps][n2par];
Double_t par2AO_err[nGaps][n2par];
Double_t par2BI[nGaps][n2par];
Double_t par2BI_err[nGaps][n2par];
Double_t par2BO[nGaps][n2par];
Double_t par2BO_err[nGaps][n2par];

/*double TOF1_Errors_X[12][3] = {(sqrt(pow(8.22,2)/12))*sin(60*PI/180), (sqrt(pow(8.22,2)/12))*sin(60*PI/180), (sqrt(pow(8.22,2)/12))*sin(60*PI/180),
                               (sqrt(pow(8.22,2)/12))*sin(30*PI/180), (sqrt(pow(8.22,2)/12))*sin(30*PI/180), (sqrt(pow(8.22,2)/12))*sin(30*PI/180),
                               0., 0., 0.,
                               (sqrt(pow(8.22,2)/12))*sin(30*PI/180), (sqrt(pow(8.22,2)/12))*sin(30*PI/180), (sqrt(pow(8.22,2)/12))*sin(30*PI/180),
                               (sqrt(pow(8.22,2)/12))*sin(60*PI/180), (sqrt(pow(8.22,2)/12))*sin(60*PI/180), (sqrt(pow(8.22,2)/12))*sin(60*PI/180),
                               (sqrt(pow(8.22,2)/12)), (sqrt(pow(8.22,2)/12)), (sqrt(pow(8.22,2)/12)),
                               (sqrt(pow(8.22,2)/12))*sin(60*PI/180), (sqrt(pow(8.22,2)/12))*sin(60*PI/180), (sqrt(pow(8.22,2)/12))*sin(60*PI/180),
                               (sqrt(pow(8.22,2)/12))*sin(30*PI/180), (sqrt(pow(8.22,2)/12))*sin(30*PI/180), (sqrt(pow(8.22,2)/12))*sin(30*PI/180),
                               0., 0., 0.,
                               (sqrt(pow(8.22,2)/12))*sin(30*PI/180), (sqrt(pow(8.22,2)/12))*sin(30*PI/180), (sqrt(pow(8.22,2)/12))*sin(30*PI/180),
                               (sqrt(pow(8.22,2)/12))*sin(60*PI/180), (sqrt(pow(8.22,2)/12))*sin(60*PI/180), (sqrt(pow(8.22,2)/12))*sin(60*PI/180),
                               (sqrt(pow(8.22,2)/12)), (sqrt(pow(8.22,2)/12)), (sqrt(pow(8.22,2)/12))};


double TOF1_Errors_Y[12][3] = {(sqrt(pow(8.22,2)/12))*cos(60*PI/180), (sqrt(pow(8.22,2)/12))*cos(60*PI/180), (sqrt(pow(8.22,2)/12))*cos(60*PI/180),
                               (sqrt(pow(8.22,2)/12))*cos(30*PI/180), (sqrt(pow(8.22,2)/12))*cos(30*PI/180), (sqrt(pow(8.22,2)/12))*cos(30*PI/180),
                               (sqrt(pow(8.22,2)/12)), (sqrt(pow(8.22,2))/12), (sqrt(pow(8.22,2))/12),
                               (sqrt(pow(8.22,2)/12))*cos(30*PI/180), (sqrt(pow(8.22,2)/12))*cos(30*PI/180), (sqrt(pow(8.22,2)/12))*cos(30*PI/180),
                               (sqrt(pow(8.22,2)/12))*cos(60*PI/180), (sqrt(pow(8.22,2)/12))*cos(60*PI/180), (sqrt(pow(8.22,2)/12))*cos(60*PI/180),
                               0., 0., 0.,
                               (sqrt(pow(8.22,2)/12))*cos(60*PI/180), (sqrt(pow(8.22,2)/12))*cos(60*PI/180), (sqrt(pow(8.22,2)/12))*cos(60*PI/180),
                               (sqrt(pow(8.22,2)/12))*cos(30*PI/180), (sqrt(pow(8.22,2)/12))*cos(30*PI/180), (sqrt(pow(8.22,2)/12))*cos(30*PI/180),
                               (sqrt(pow(8.22,2)/12)), (sqrt(pow(8.22,2)/12)), (sqrt(pow(8.22,2)/12)),
                               (sqrt(pow(8.22,2)/12))*cos(30*PI/180), (sqrt(pow(8.22,2)/12))*cos(30*PI/180), (sqrt(pow(8.22,2)/12))*cos(30*PI/180),
                               (sqrt(pow(8.22,2)/12))*cos(60*PI/180), (sqrt(pow(8.22,2)/12))*cos(60*PI/180), (sqrt(pow(8.22,2)/12))*cos(60*PI/180),
                               0., 0., 0.};
*/       

double sigma_par = sqrt(pow(8.22,2)/12);   // sigma_parallel
double sigma_perp = sqrt(pow(5.,2)/12);    // sigma perpendicular                 

double TOF1_Errors_X[12][3] = {sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), 
                               sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180),
                               sigma_par*cos(90*PI/180) + sigma_perp*cos(00*PI/180), sigma_par*cos(90*PI/180) + sigma_perp*cos(00*PI/180), sigma_par*cos(90*PI/180) + sigma_perp*cos(00*PI/180),
                               sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180),
                               sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180),
                               sigma_par*cos(00*PI/180) + sigma_perp*cos(90*PI/180), sigma_par*cos(00*PI/180) + sigma_perp*cos(90*PI/180), sigma_par*cos(00*PI/180) + sigma_perp*cos(90*PI/180),
                               sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180),
                               sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180),
                               sigma_par*cos(90*PI/180) + sigma_perp*cos(00*PI/180), sigma_par*cos(90*PI/180) + sigma_perp*cos(00*PI/180), sigma_par*cos(90*PI/180) + sigma_perp*cos(00*PI/180),
                               sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180), sigma_par*cos(60*PI/180) + sigma_perp*cos(30*PI/180),
                               sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180), sigma_par*cos(30*PI/180) + sigma_perp*cos(60*PI/180),
                               sigma_par*cos(00*PI/180) + sigma_perp*cos(90*PI/180), sigma_par*cos(00*PI/180) + sigma_perp*cos(90*PI/180), sigma_par*cos(00*PI/180) + sigma_perp*cos(90*PI/180)};


double TOF1_Errors_Y[12][3] = {sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), 
                               sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180),
                               sigma_par*sin(90*PI/180) + sigma_perp*sin(00*PI/180), sigma_par*sin(90*PI/180) + sigma_perp*sin(00*PI/180), sigma_par*sin(90*PI/180) + sigma_perp*sin(00*PI/180),
                               sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180),
                               sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180),
                               sigma_par*sin(00*PI/180) + sigma_perp*sin(90*PI/180), sigma_par*sin(00*PI/180) + sigma_perp*sin(90*PI/180), sigma_par*sin(00*PI/180) + sigma_perp*sin(90*PI/180),
                               sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180),
                               sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180),
                               sigma_par*sin(90*PI/180) + sigma_perp*sin(00*PI/180), sigma_par*sin(90*PI/180) + sigma_perp*sin(00*PI/180), sigma_par*sin(90*PI/180) + sigma_perp*sin(00*PI/180),
                               sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180), sigma_par*sin(60*PI/180) + sigma_perp*sin(30*PI/180),
                               sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180), sigma_par*sin(30*PI/180) + sigma_perp*sin(60*PI/180),
                               sigma_par*sin(00*PI/180) + sigma_perp*sin(90*PI/180), sigma_par*sin(00*PI/180) + sigma_perp*sin(90*PI/180), sigma_par*sin(00*PI/180) + sigma_perp*sin(90*PI/180)};
