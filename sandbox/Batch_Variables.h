  //#include "TDC_Windows.h"
  
  using namespace std;
  
  double Yth_TOF1=999.;
  double New_ChiS=999.;

  //int C2X_clusters = 0;
  //int C2Y_clusters = 0;

  //vector <int> C2X_cluster_index;
  //vector <int> C2X_cluster_length;
  //vector <int> C2Y_cluster_index;
  //vector <int> C2Y_cluster_length;

  double C2X_L[56] = {0.};
  double C2X_R[56] = {0.};
  double C2Y_L[16] = {0.};
  double C2Y_R[16] = {0.};
  
  double C3X_L[64] = {0.};
  double C3X_R[64] = {0.};
  double C3Y_L[16] = {0.};
  double C3Y_R[16] = {0.};
  
  double C4X_L[72] = {0.};
  double C4X_R[72] = {0.};
  double C4Y_L[16] = {0.};
  double C4Y_R[16] = {0.};

  //bool first_cluster = true;
  //int cluster_spacing = 0;
  //int cluster_length_count = 0;

  double C2X_centroid = 999;

  vector <double> vec_xx_lepton;
  vector <double> vec_yy_lepton;
  vector <double> vec_xx_lepton_rotate;
  vector <double> vec_yy_lepton_rotate;

  vector <double> vec_xx_kaon;
  vector <double> vec_yy_kaon;
  
  vector<int> vec_tdc_ck;
  vector<int> vec_tdc_ck_col;


  double a_gr_Leptons = 999.99;
  double b_gr_Leptons = 999.99;

  int Rotate=1; // When TOF1 is 12 or 6, or when TOF2 is 6 or 12, rotate by -90 deg to fit a horizontal line (0 = OFF ; 1 = ON)
  int i_kaon_bar = 999;
  int K_stop_count = 0; 
  int T_limit = 3;
  
  char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!

  //Double_t flag_size_TARGET=1.35;
  //Double_t flag_size_SFT=1.3;
  //Double_t flag_size_palette=1.6;

  Int_t HG_TARGET_ADC_Thr[256] = {0};
  Int_t LG_TARGET_ADC_Thr[256] = {0};
  Int_t HG_SFT_ADC_Thr[128] = {0};
  Int_t LG_SFT_ADC_Thr[128] = {0};

  //Int_t TDC_min_TARGET = TARGET_TDC_min[0];
  //Int_t TDC_max_TARGET = TARGET_TDC_max[0];
  //Int_t ADC_TARGET_Thr = HG_TARGET_ADC_Thr[0];
  //Int_t TDC_min_SFT = SFT_TDC_min[0];
  //Int_t TDC_max_SFT = SFT_TDC_max[0];
  
  bool TARGET_High_has_neighbours[256] = {false}; // Array of High gain target hits which have no neighbouring targets hit   
  
  int n_hit = 2;   // ## Number of hits required in the TARGET

  Int_t ADC_High_corr_max=0;

  Int_t ADC_cut_TARGET2 = 850;

  char ADC_cut[100];    //sprintf(ADC_cut,"(ADC >= %d)",HG_SFT_ADC_Thr[0]);
  
  int tdc_trigger[2][16];
  double tdc_ck_corr = 0.;
  double TDC_ck_avg =0.;
  double TDC_ck_avg2 = 0.;
  int TDC_ck_sum = 0;  int TDC_ck_counter = 0;
  int TDC_ck_sum2 = 0;
  double TDC_ck_sigma = 0.;  double TDC_ck_sigma2 = 0.;
  double ck_mean = 1;
  double TDC_diff = -1;


  Int_t adc_high_target[256];       Int_t ADC_High_TARGET[256];    
  Int_t adc_low_target[256];        Int_t ADC_Low_TARGET[256];
  Int_t ADC_Low_TARGET_ped[256];    Int_t ADC_High_TARGET_ped[256];
  Int_t tdc_le_target[256][16];     //Int_t TDC_LE_TARGET[256];     
  Int_t tdc_te_target[256][16];     //Int_t TDC_TE_TARGET[256];  
  Int_t TDC_LE_TARGET_corrected[256][6]={0};
  char TDC_LE_TARGET_corr[256][6][20];

  Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];        Double_t ADC_High_SFT_corr[128];    
  Int_t adc_low_sft[128];           Int_t ADC_Low_SFT[128];   
  Int_t tdc_le_sft[128][16];        Int_t TDC_LE_SFT[128];         
  Int_t tdc_te_sft[128][16];        Int_t TDC_TE_SFT[128];  

  Int_t ADC_TOF1[24];   
  //Int_t ADC_TOF2[48];  // !!!
   
  Int_t ADC_tof1U[12];              Int_t ADC_TOF1U[12];
  Int_t ADC_tof1D[12];              Int_t ADC_TOF1D[12];

  Int_t TDC_tof1U[12];              Int_t TDC_TOF1U[12];
  Int_t TDC_tof1D[12];              Int_t TDC_TOF1D[12];

  Int_t ADC_tof2AO[12];             Int_t ADC_TOF2AO[12];
  Int_t ADC_tof2BO[12];             Int_t ADC_TOF2BO[12];
  Int_t ADC_tof2AI[12];             Int_t ADC_TOF2AI[12];
  Int_t ADC_tof2BI[12];             Int_t ADC_TOF2BI[12];

  Int_t TDC_tof2AO[12];             Int_t TDC_TOF2AO[12];
  Int_t TDC_tof2BO[12];             Int_t TDC_TOF2BO[12];
  Int_t TDC_tof2AI[12];             Int_t TDC_TOF2AI[12];
  Int_t TDC_tof2BI[12];             Int_t TDC_TOF2BI[12];

  Int_t MwpcADC[512];               Int_t MWPCADC[512];

  Int_t adc_c2x_r[56];              double ADC_C2X_R[56];
  Int_t adc_c2x_l[56];              double ADC_C2X_L[56];
  Int_t adc_c2y_r[16];              double ADC_C2Y_R[16];
  Int_t adc_c2y_l[16];              double ADC_C2Y_L[16];

  Int_t adc_c3x_r[64];              double ADC_C3X_R[64];
  Int_t adc_c3x_l[64];              double ADC_C3X_L[64];
  Int_t adc_c3y_r[16];              double ADC_C3Y_R[16];
  Int_t adc_c3y_l[16];              double ADC_C3Y_L[16];

  Int_t adc_c4x_r[72];              double ADC_C4X_R[72];
  Int_t adc_c4x_l[72];              double ADC_C4X_L[72];
  Int_t adc_c4y_r[16];              double ADC_C4Y_R[16];
  Int_t adc_c4y_l[16];              double ADC_C4Y_L[16];  

  Int_t tdc_ck[14][16];
  Int_t tdc_cpi[14][16];

  Int_t has_TDC_SFT_hit[128] = {0};

  Bool_t Event_flag[40] = {false};
  
  float R_TARGET = 29.0;
  float R_TOF1 = 47.1;
  float R_SFT_L1 = 40.0;
  //float R_C2 = 629.4;

  char run_string[100];            char event_string[100];

  char h_ADC_title[200];


  double Z_TOF1[12] = {999.99};

  const double Offset_TOF1[12] = {180., 410., 0., 0., 310., 280., -30., 70., 0., 0., 390., 430.};
  const double Slope_correction_TOF1[12] = {0.70, 0.78, 1., 1., 0.70, 0.60, 0.86, 0.77, 1., 1., 0.67, 0.78};
  

  char file_mapping[200];
  //sprintf(file_mapping,"../Mapping");

  char par_finput[200];
  //sprintf(par_finput,"%s/%s",file_mapping,source_mapping);

  char par_finput2[200];
  //sprintf(par_finput2,"%s/MWPC_map2.txt",file_mapping);
  Int_t par_temp[2][128];
  //ifstream fdat(par_finput,ios::in);
  //for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
  //fdat.close();

  char par_temp2[512][50];
  //ifstream fdat2(par_finput2,ios::in);
  //for(Int_t ii=0; ii<512; ii++) fdat2 >> par_temp2[ii];
  //fdat2.close();
  
 /////////////////////////////////////////////////////////////////////////////////////////////////////

  char h_target_ADC_title[100];
  char h_target_TDC_title[100];
  //sprintf(h_target_ADC_title,"Run %d (Event %d) -- ADC #geq %d)",Run_Number, ievt, par_temp_TARGET[1][0]);
  //sprintf(h_target_ADC_title,"Run %d (Event %d) -- ADC #geq %d)",Run_Number, ievt_min, ADC_TARGET_Thr);
  //sprintf(h_target_TDC_title,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt, par_temp_TARGET[1][0],TDC_min_TARGET,TDC_max_TARGET);
  //sprintf(h_target_TDC_title,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt_min, ADC_TARGET_Thr,TDC_min_TARGET,TDC_max_TARGET);

  char h_target_ADC_title2[100];
  char h_target_TDC_title2[100];
  char h_target_ADC_title3[100];
  char h_target_ADC_title4[100];
  //sprintf(h_target_ADC_title2,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt_min, ADC_TARGET_Thr,TDC_min_TARGET, TDC_max_TARGET);
  //sprintf(h_target_TDC_title2,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt_min, ADC_TARGET_Thr, TDC_min_TARGET,TDC_max_TARGET);
  //sprintf(h_target_ADC_title3,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt_min, ADC_TARGET_Thr, TDC_min_TARGET,TDC_max_TARGET);
  //sprintf(h_target_ADC_title4,"Run %d (Event %d) -- ADC #geq %d | %d #leq TDC #leq %d",Run_Number, ievt_min, ADC_TARGET_Thr, TDC_min_TARGET,TDC_max_TARGET);
  
  //sprintf(h_ADC_title,"(ADC offset = %d) | (%d < TDC < %d)  --  Run %d (Event %d)",SFT_ADC_Thr_HG_Offset, TDC_min_SFT,TDC_max_SFT,Run_Number,ievt_min);

  char h_target_TDC_copy_Name[200];     char h_target_TDC_copy_Title[200];
  char h_TDC_selected_Name[200];        char h_TDC_selected_Title[200];
  char h_GoodLG_Name[200];              char h_GoodLG_Title[200];
  char h_TDC_selected2_Name[200];       char h_TDC_selected2_Title[200];
  char h_GoodLG_weighted_Name[200];     char h_GoodLG_weighted_Title[200];
  char h_Gap_Fibers_Name[200];          char h_Gap_Fibers_Title[200];

  //sprintf(h_target_TDC_copy_Name,"Event %d (Run %d)", ievt_min, Run_Number);
  //sprintf(h_target_TDC_copy_Title,"ADC_HG #geq %d  |  %d #leq TDC #leq %d", ADC_TARGET_Thr, TDC_min_TARGET, TDC_max_TARGET);
 
  //sprintf(h_TDC_selected_Name,"Event %d (Run %d)", ievt_min, Run_Number);
  //sprintf(h_TDC_selected_Title,"ADC_HG #geq %d  |  %d #leq TDC #leq %d", ADC_TARGET_Thr+100, TDC_min_TARGET, TDC_max_TARGET);

  //sprintf(h_GoodLG_Name,"Event %d (Run %d)", ievt_min, Run_Number);
  //sprintf(h_GoodLG_Title,"ADC_LG #geq %d  |  %d #leq TDC #leq %d", ADC_TARGET_Thr,TDC_min_TARGET,TDC_max_TARGET);

  //sprintf(h_TDC_selected2_Name,"Event %d (Run %d)", ievt_min, Run_Number);
  //sprintf(h_TDC_selected2_Title,"ADC_HG #geq %d  |  %d #leq TDC #leq %d  (WEIGHTED)", ADC_TARGET_Thr+100,TDC_min_TARGET,TDC_max_TARGET);

  //sprintf(h_GoodLG_weighted_Name,"Event %d (Run %d)", ievt_min, Run_Number);
  //sprintf(h_GoodLG_weighted_Title,"ADC_LG #geq %d  |  %d #leq TDC #leq %d  (WEIGHTED)", ADC_TARGET_Thr,TDC_min_TARGET,TDC_max_TARGET);

  //sprintf(h_Gap_Fibers_Name,"Event %d (Run %d)", ievt_min, Run_Number);
  //sprintf(h_Gap_Fibers_Title,"ADC LG offset =  %d  |  %d #leq TDC #leq %d", TARGET_ADC_Thr_LG_Offset,TDC_min_TARGET,TDC_max_TARGET);
  

  //TH1F *h_ADC_L1_DS = new TH1F("h_ADC_L1_DS",h_ADC_title,130,-1,129);
  //TH1F *h_ADC_L2_DS = new TH1F("h_ADC_L2_DS",h_ADC_title,130,-1,129);
  //TH1F *h_ADC_L3_DS = new TH1F("h_ADC_L3_DS",h_ADC_title,130,-1,129);
  //TH1F *h_ADC_L4_DS = new TH1F("h_ADC_L4_DS",h_ADC_title,130,-1,129);

  //TH1F *h_ADC_L1_US = new TH1F("h_ADC_L1_US",h_ADC_title,66,63,129);
  //TH1F *h_ADC_L2_US = new TH1F("h_ADC_L2_US",h_ADC_title,66,63,129);
  //TH1F *h_ADC_L3_US = new TH1F("h_ADC_L3_US",h_ADC_title,66,63,129);
  //TH1F *h_ADC_L4_US = new TH1F("h_ADC_L4_US",h_ADC_title,66,63,129);

  /*
  TH2F *h_Circle = new TH2F("Test", "TARGET and TOF1", 500, -50, 50, 500, -50, 50);

  TH2F *h_target_ADC = new TH2F("Histo Fit 0",h_target_ADC_title,3000,-50,50,3000,-50,50);
  TH2F *h_target_TDC = new TH2F("Histo Fit 1",h_target_TDC_title,3000,-50,50,3000,-50,50);
  TH2F *h_target_TDC_copy = new TH2F(h_target_TDC_copy_Name,h_target_TDC_copy_Title,3000,-50,50,3000,-50,50);
  TH2F *h_TDC_selected = new TH2F(h_TDC_selected_Name, h_TDC_selected_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_TDC_selected2 = new TH2F(h_TDC_selected2_Name, h_TDC_selected2_Title, 500, -50, 50, 500, -50, 50);

  TH2F *h_target_ADC2 = new TH2F("Histo Fit 2",h_target_ADC_title2,3000,-50,50,3000,-50,50);
  TH2F *h_target_TDC2 = new TH2F("Histo Fit 3",h_target_TDC_title2,3000,-50,50,3000,-50,50);

  TH2F *h_target_ADC3 = new TH2F("Histo Fit 4",h_target_ADC_title3,3000,-50,50,3000,-50,50);
  TH2F *h_target_ADC4 = new TH2F("Histo Fit 5",h_target_ADC_title4,3000,-50,50,3000,-50,50);
  TH2F *h_target_ADCA = new TH2F("Histo Fit 6",h_target_ADC_title3,3000,-50,50,3000,-50,50);

  TH2F *h_Target = new TH2F("Test Target", "Target", 3000, -50, 50, 3000, -50, 50);
  TH2F *h_Target_Center = new TH2F("Target Center", "Target Center", 500, -50, 50, 500, -50, 50);

  TH2F *h_kaon = new TH2F("Kaon", "Kaon", 500, -50, 50, 500, -50, 50);
  TH2F *h_kaon_copy = new TH2F("Kaon Copy", "Kaon Copy", 500, -50, 50, 500, -50, 50);
    
  TH2F *h_max = new TH2F("Max", "Max", 500, -50, 50, 500, -50, 50);
  TH2F *h_max_copy = new TH2F("Max Copy", "Max Copy", 500, -50, 50, 500, -50, 50);

  TH2F *h_GoodLG = new TH2F(h_GoodLG_Name, h_GoodLG_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_GoodLG_weighted = new TH2F(h_GoodLG_weighted_Name, h_GoodLG_weighted_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_GoodLG_copy = new TH2F("Good LG Copy", "Good LG Copy", 500, -50, 50, 500, -50, 50);

  TH2F *h_TOF1 = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_closest = new TH2F("h_TOF1_closest", "h_TOF1_closest", 500, -50, 50, 500, -50, 50);
  TH2F *h_TOF1_rotate = new TH2F("h_TOF1", "h_TOF1", 500, -50, 50, 500, -50, 50); // ROTATE_CHANGE
  TH2F *h_TOF1_closest_rotate = new TH2F("h_TOF1_closest", "h_TOF1_closest", 500, -50, 50, 500, -50, 50); //ROTATE_CHANGE

  TH2F *h_int_TDC = new TH2F("h_int_TDC", "h_int_TDC", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_selected = new TH2F("h_int_TDC_selected", "h_int_TDC_selected", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_selected_weighted = new TH2F("h_int_TDC_selected_weighted", "h_int_TDC_selected_weighted", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_Gap_Fibers = new TH2F("h_int_TDC_Gap_Fibers", "h_int_TDC_Gap_Fibers", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_Gap_Fibers_SFT = new TH2F("h_int_TDC_Gap_Fibers_SFT","h_int_TDC_Gap_Fibers_SFT", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_Gap_Fibers_rotate = new TH2F("h_int_TDC_Gap_Fibers", "h_int_TDC_Gap_Fibers", 500, -50, 50, 500, -50, 50);  // ROTATE_CHANGE
  TH2F *h_int_TDC_Gap_Fibers_SFT_rotate = new TH2F("h_int_TDC_Gap_Fibers_SFT","h_int_TDC_Gap_Fibers_SFT", 500, -50, 50, 500, -50, 50);  // ROTATE_CHANGE
  TH2F *h_int_TDC_SFT = new TH2F("h_int_TDC_SFT", "h_int_TDC_SFT", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_TDC_TARGET = new TH2F("h_int_TDC_TARGET", "h_int_TDC_TARGET", 500, -50, 50, 500, -50, 50);


  TH2F *h_int_GoodLG = new TH2F("h_int_GoodLG", "h_int_GoodLG", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_GoodLG_weighted = new TH2F("h_int_GoodLG_weighted", "h_int_GoodLG_weighted", 500, -50, 50, 500, -50, 50);
  TH2F *h_int_GoodLG_SFT = new TH2F("h_int_GoodLG_SFT", "h_int_GoodLG_SFT", 500, -50, 50, 500, -50, 50);
    
  TH2F *h_TDC_Gap_Fibers = new TH2F(h_Gap_Fibers_Name, h_Gap_Fibers_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_K_Stop_Bars = new TH2F(h_Gap_Fibers_Name, h_Gap_Fibers_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_TDC_Gap_Fibers_copy = new TH2F(h_Gap_Fibers_Name, h_Gap_Fibers_Title, 500, -50, 50, 500, -50, 50); //ROTATE_CHANGE
  TH2F *h_TDC_Gap_Fibers_kaon = new TH2F(h_Gap_Fibers_Name, h_Gap_Fibers_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_TARGET_LG_Blue = new TH2F(h_Gap_Fibers_Name, h_Gap_Fibers_Title, 500, -50, 50, 500, -50, 50);
  TH2F *h_Centroid = new TH2F("Centroid", "Centroid", 500, -50, 50, 500, -50, 50);
  TEllipse *ell = new TEllipse(0, 0, R_TOF1, 0);
  TEllipse *ell_Target = new TEllipse(0, 0, R_TARGET, 0);
  TEllipse *ell_L1 = new TEllipse(0, 0, R_SFT_L1, 0);
  

  TLine *HorizontalAxis = new TLine(0., 0., 50., 0.);
  */


  char path_input[200];                   
  //sprintf(path_input,"%s",path_merged);          


  char Name_finput[200];
  //sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);

  char footer[100];
  //sprintf(footer,"Event_Display_MS.C  --  Run %d ; Event %d",Run_Number,ievt_min);


  //sprintf(run_string,"Run %d ; Event %d",Run_Number,ievt_min);
  //sprintf(event_string,"Run %d ; Event %d",Run_Number,ievt_min);

  /*
  TChain *fChain= new TChain("Tree");   
  fChain->Add(Name_finput);   
  fChain->SetMakeClass(1);              

  fChain->SetBranchAddress("ADC_High_TARGET",adc_high_target);    fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
  fChain->SetBranchAddress("ADC_Low_TARGET",adc_low_target);      fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
  fChain->SetBranchAddress("TDC_LE_TARGET",tdc_le_target);        fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
  fChain->SetBranchAddress("TDC_TE_TARGET",tdc_te_target);        fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);

  fChain->SetBranchAddress("ADC_TOF1U",ADC_tof1U);
  fChain->SetBranchAddress("ADC_TOF1D",ADC_tof1D);
  fChain->SetBranchAddress("TDC_TOF1U",TDC_tof1U);
  fChain->SetBranchAddress("TDC_TOF1D",TDC_tof1D); 
  
  fChain->SetBranchAddress("ADC_TOF2AO",ADC_tof2AO);
  fChain->SetBranchAddress("ADC_TOF2AI",ADC_tof2AI);
  fChain->SetBranchAddress("ADC_TOF2BO",ADC_tof2BO);
  fChain->SetBranchAddress("ADC_TOF2BI",ADC_tof2BI);    
  fChain->SetBranchAddress("TDC_TOF2AO",TDC_tof2AO);
  fChain->SetBranchAddress("TDC_TOF2AI",TDC_tof2AI);
  fChain->SetBranchAddress("TDC_TOF2BO",TDC_tof2BO);
  fChain->SetBranchAddress("TDC_TOF2BI",TDC_tof2BI);    

  fChain->SetBranchAddress("MWPCADC",MwpcADC);

  fChain->SetBranchAddress("ADC_C2X_R",adc_c2x_r);
  fChain->SetBranchAddress("ADC_C2X_L",adc_c2x_l);
  fChain->SetBranchAddress("ADC_C2Y_R",adc_c2y_r);
  fChain->SetBranchAddress("ADC_C2Y_L",adc_c2y_l);
  fChain->SetBranchAddress("ADC_C3X_R",adc_c3x_r);
  fChain->SetBranchAddress("ADC_C3X_L",adc_c3x_l);
  fChain->SetBranchAddress("ADC_C3Y_R",adc_c3y_r);
  fChain->SetBranchAddress("ADC_C3Y_L",adc_c3y_l);
  fChain->SetBranchAddress("ADC_C4X_R",adc_c4x_r);
  fChain->SetBranchAddress("ADC_C4X_L",adc_c4x_l);
  fChain->SetBranchAddress("ADC_C4Y_R",adc_c4y_r);
  fChain->SetBranchAddress("ADC_C4Y_L",adc_c4y_l);  

  fChain->SetBranchAddress("TDC_Ck", tdc_Ck);
  fChain->SetBranchAddress("TDC_Cpi", tdc_Cpi);

  fChain->SetBranchAddress("EvFlag", Event_flag);
 
  //Int_t nentries = (Int_t)fChain->GetEntries();
  */


  
  bool Good_Event=false;
  //bool Good_tof1[12] = {false};
  //bool Good_tof2[12] = {false};

  bool Good_TARGET_Event = false;
  int count_TARGET_evts = 0; 

  bool Good_TOF1_ADC[12]={false};   bool Good_TOF2_ADC[12]={false};
  bool Good_TOF1_TDC[12]={false};   bool Good_TOF2_TDC[12]={false};
  bool Good_TOF1[12]={false};       bool Good_TOF2[12]={false};
  bool Good_TOFs[12]={false};
  bool Good_TOF_Event = false;

  bool Event_On_Blacklist = false;
  string Current_Event;
  //int current_event;

  int count_C2X = 0;    int count_C2Y = 0;
  int count_C3X = 0;    int count_C3Y = 0;
  int count_C4X = 0;    int count_C4Y = 0;
  bool Good_MWPC_Event = false;

  bool has_TDC_hit[256] = {false};
  bool has_TDC_hit_Kstop[256] = {false};

  char ch_ADC_cut_TARGET[100];;
  char ch_ADC_and_TDC_cut[100];
  char ch_ADC_and_TDC_cut_Kstop[100];

  Int_t Angle_ADC_cut = 0;

  double x_inc = 0;
  double y_inc = 0;
  int hit_count = 0;
  //int count = 0; 

  int max_index = 0;
  int max_index2 = 0;
  int max_index3 = 0;
  int max_index4 = 0;

  int max_ADC = -100000000;
  int max_ADC2 = -100000000;
  int max_ADC3 = -100000000;
  int max_ADC4 = -100000000;  

  double x_cent;
  double y_cent;
  
  double hyp[256] = {-1};
  
  bool has_ADC_TOF1_hit[12] = {false};
  bool has_TDC_TOF1_hit[12] = {false};
  bool has_ADC_TOF2_hit[12] = {false};
  bool has_TDC_TOF2_hit[12] = {false};
  //bool has_both_ADC_TOF1_hit[12] = {false};

  //int gap_hit[12] = {0};
  int ADC_TOF1_hit[12] = {0};
  int ADCTDC_TOF1_hit[12] = {0};
  int ADC_TOF2_hit[12] = {0};
  int ADCTDC_TOF2_hit[12] = {0};

  int selected_TOF2 = 0;

  int gap_counter[12] = {0};
  int scoring_type = 2;         // scoring_type = 1  --->  Oscar's Method

  int high_gap_hit = 0;
  int gap_to_fit = 0;
  int score_max = 0;
  
  vector<int> tof1_ties;
  int tdc_tof1u_min = -10000;
  int tdc_tof1d_min = -10000;
  
  bool k_stop_bar[256] = {false};
  vector<int> good_k_stop_bars;
  
  double X_BAR = 0;
  double Y_BAR = 0;
  
  double total_energy = 0.0;
  double X_weights = 0.0;
  double Y_weights = 0.0;   

  Int_t TDC_average = -1;

  int max_index_all[256];
  int max_ADC_all[256];
  int max_index_flag;

  int TDC_LG_max = -1;
  int TDC_LG_max2 = -1;
  //int index_max1=0;
  //int index_max2=0;
  
  //int C2X_L[56] = {0};
  //int C2X_R[56] = {0};
  //int C2Y_L[16] = {0};
  //int C2Y_R[16] = {0};
  //int C3X_L[64] = {0};
  //int C3X_R[64] = {0};
  //int C3Y_L[16] = {0};
  //int C3Y_R[16] = {0};
  //int C4X_L[72] = {0};
  //int C4X_R[72] = {0};
  //int C4Y_L[16] = {0};
  //int C4Y_R[16] = {0};

  //double C2X_centroid = 0.0;

  bool first_cluster = true;
  int cluster_spacing = 0;
  int cluster_length_count = 0;

  int C2X_clusters = 0;
  int C2Y_clusters = 0;
  int C3X_clusters = 0;
  int C3Y_clusters = 0;
  int C4X_clusters = 0;
  int C4Y_clusters = 0;
  
  vector<int> C2X_cluster_index; // Hold the starting index of each cluster.
  vector<int> C2X_cluster_length; // Hold the length of each cluster.
  vector<int> C2Y_cluster_index;
  vector<int> C2Y_cluster_length;

  vector<int> C3X_cluster_index;
  vector<int> C3X_cluster_length;
  vector<int> C3Y_cluster_index;
  vector<int> C3Y_cluster_length;    

  vector<int> C4X_cluster_index;
  vector<int> C4X_cluster_length;  
  vector<int> C4Y_cluster_index;
  vector<int> C4Y_cluster_length;  
  
  double C2_centroid_num = 0;
  double C2_centroid_den = 0;
  double C2X_cenroid = 999.;

  int TDC_Ck_counter = 0;
  int TDC_Cpi_counter = 0;

  Double_t par0_ADC = 0;  float par0_TDC = 0.;
  Double_t par1_ADC = 0;  float par1_TDC = 0.;
  //Double_t phi_TDC = 0;

  //char angle_string_ADC[100]; char angle_string_TDC[100];
  
  //TLatex *tex_angle_ADC;  TLatex *tex_angle_TDC;

  double closest_to_centroid = 1000;
  //int closest_to_centroid_index = -1;

  //char ch_centroid[200];    

  int gap_to_fit_rotate = 0;
  
  int ADC_High_TARGET_temp[256];
  int ADC_Low_TARGET_temp[256];
  //int adc_low_target_temp[256];
  bool has_TDC_hit_temp[256];
  bool has_TDC_hit_Kstop_temp[256];
  int HG_TARGET_ADC_Thr_temp[256];
  int LG_TARGET_ADC_Thr_temp[256];
  int TARGET_High_has_neighbours_temp[256];
  bool k_stop_bar_temp[256];

  int has_data_TDC2 = 0;
  int has_data_ADC2 = 0;
  int has_data_ADC3 = 0;
  int has_data_ADC4 = 0;
  int has_data_ADCA = 0;

  double Xloc_gap = 0;
  double Yloc_gap = 0;

  double xdistance1;
  double ydistance1;

  double xdistance2;
  double ydistance2;

  double xdistance3;
  double ydistance3;

  double xhyp1;
  double xhyp2;
  double xhyp3;

  int closest_gap_point = 9999;
  
  //double xcoord = 0;
  //Int_t unique_x = 0;

  //float determinant;
  //float x_circle_int1;
  //float y_circle_int1;

  //float x_circle_int2;
  //float y_circle_int2;

  //double SFTxdistance1;
  //double SFTydistance1;

  //double SFTxdistance2;
  //double SFTydistance2;

  //double SFTxhyp1;
  //double SFTxhyp2;

  //double SFT_x_intercept = 999.;
  //double SFT_y_intercept = 999.;

  //double SFT_phi = 999.;

  double x_intercept = 999.;
  double y_intercept = 999.;

  double x_distance4[256] = {0};
  double y_distance4[256] = {0};
  int distances4[256];

  double min_distance = 9999.;

  float Gap[12][3][2] = {{{0}}};
  
  float a_fit_TDC_selected=0.;               float b_fit_TDC_selected=0.;
  float a_fit_GoodLG=0.;                     float b_fit_GoodLG=0.;
  float a_fit_GoodLG_weighted=0.;            float b_fit_GoodLG_weighted=0.;
  float a_fit_TDC_selected_weighted=0.;      float b_fit_TDC_selected_weighted=0.;
  float a_fit_TDC_Gap_Fibers=0.;             float b_fit_TDC_Gap_Fibers=0.;

  float x_int_TDC[2];                           float y_int_TDC[2];
  float x_int_TDC_selected[2];                  float y_int_TDC_selected[2];
  float x_int_TDC_selected_weighted[2];         float y_int_TDC_selected_weighted[2];
  float x_int_GoodLG[2];                        float y_int_GoodLG[2];
  float x_int_GoodLG_weighted[2];               float y_int_GoodLG_weighted[2];
  float x_int_TDC_Gap_Fibers[2];                float y_int_TDC_Gap_Fibers[2];
  float x_int_TDC_SFT[2];                       float y_int_TDC_SFT[2];
  float x_int_GoodLG_SFT[2];                    float y_int_GoodLG_SFT[2];
  float x_int_TDC_Gap_Fibers_SFT[2];            float y_int_TDC_Gap_Fibers_SFT[2];
  float x_int_TARGET[2];                        float y_int_TARGET[2];

  float x_TDC_sel_intersect1=0.;   float y_TDC_sel_intersect1=0.;
  float x_GoodLG_intersect1=0.;    float y_GoodLG_intersect1=0.;
  float x_TDC_Gap_Fibers=0.;       float y_TDC_Gap_Fibers=0.;
  //float x_TDC_Gap_Fibers_SFT=0.;   float y_TDC_Gap_Fibers_SFT=0.;

  float dist1_TDC_selected[2];
  float dist1_GoodLG[2];
  float dist1_TDC_Gap_Fibers[2];
  float dist1_TDC_Gap_Fibers_SFT[2];

  float dist2_TDC_selected[3];
  float dist2_GoodLG[3];

  float dist2_TDC_selected_min = 1000.;
  float dist2_GoodLG_min = 1000.;
  int selected_TDC_selected = 0;

  double ParError = 999.99;
  double ChiS = 0.0;
  double ndf = 0.0;
  double prob = 0.0;
  
  float x_TDC_selected_weighted_intersect1=0.;   float y_TDC_selected_weighted_intersect1=0.;
  float x_GoodLG_weighted_intersect1=0.;         float y_GoodLG_weighted_intersect1=0.;
  float x_TDC_SFT_intersect1=0.;                 float y_TDC_SFT_intersect1=0.;
  float x_GoodLG_SFT_intersect1=0.;              float y_GoodLG_SFT_intersect1=0.;
  float x_TDC_Gap_Fibers_intersect1=0.;          float y_TDC_Gap_Fibers_intersect1=0.;
  float x_TDC_Gap_Fibers_SFT_intersect1=0.;      float y_TDC_Gap_Fibers_SFT_intersect1=0.;
  float x_TARGET_intersect=0;                    float y_TARGET_intersect=0;
 
  float dist1_TDC_selected_weighted[2];
  float dist1_TDC_SFT[2];
  float dist1_GoodLG_weighted[2];
  float dist1_GoodLG_SFT[2];
  float dist1_TARGET_intersect[2];

  double dist_to_k_stop = 0.;

  float a_final_TDC = 0.;           float a_final_GoodLG = 0.;             float a_final_Gap_Fibers = 0.;       float a_final_guide = 0.;
  float alpha_TDC = 0.;             float alpha_GoodLG = 0.;               float alpha_Gap_Fibers = 0.;         float alpha_guide = 0.;
  float tanalpha_TDC = 0.;          float tanalpha_GoodLG = 0.;            float tanalpha_Gap_Fibers = 0.;      float tanalpha_guide = 0.;
  float angle_final_TDC = 0.;       float angle_final_GoodLG = 0.;         float angle_final_Gap_Fibers = 0.;   float angle_final_guide = 0.;

  double Delta_phi = 999.99;  double Delta_phi_deg = 999.99;

  int Axis_Vector_Length = 10;

  float x_Kstop = 0.;     float x_SFT1 = 0.;
  float y_Kstop = 0.;     float y_SFT1 = 0.;
  float x2_Keito = 0.;    float y2_Keito = 0.;
  float x3_Keito = 0.;    float y3_Keito = 0.;

  float length = 20;
  float angle_Keito = 0.;

  //float a_final_Keito = 0.;

  char Angle_TDC_string[30];      char Angle_GoodLD_string[30];      char Angle_Gap_Fibers_string[30];      char Angle_guide_string[30];
  char Angle_Keito_string[30];

  char X_SFT_String[30];
  char Y_SFT_String[30];
  char X_TOF1_String[30];
  char Y_TOF1_String[30];

  double length_in_target = 0;

  double x_tof1_intersect_1 = 0;
  double y_tof1_intersect_1 = 0;
  double x_tof1_intersect_2 = 0;
  double y_tof1_intersect_2 = 0;
  double x_tof1_intersect = 0;
  double y_tof1_intersect = 0;

  double alpha = 0.;
  
  vector<double> ZZ;  //ZZ.clear();
  vector<double> DZ;  //DZ.clear();
  vector<double> Track_Length;  //Track_Length.clear();
  vector<double> Z_selected;  //Z_selected.clear();
  
  //double TOF1_Z_range = 25.; 
  //double TOF1_Z_min = 999.;
  //double TOF1_Z_max = 999.;

  vector<double> Z_selected_sorted;   //Z_selected_sorted.clear();
  
  char output[100];
  //sprintf(output,"RUN_%d_DATA.dat",Run_Number);
  //ofstream fout;
  //fout.open(output);

  char output2[100];
  //sprintf(output2,"Run_%d_Zlist.dat",Run_Number);
  //ofstream fout2;
  //fout2.open(output2);
