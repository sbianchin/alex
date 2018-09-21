vector<double> ZZ;
vector<double> Z_selected;


  bool To_rotate = false;
  int error_flag = 1;
  int T_limit = 3;
  char path_input[200];
  char Name_finput[200];
  char Version[100] = "Version 5.5";  // TO MOVE
  char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!

  Int_t HG_TARGET_ADC_Thr[256] = {0};
  Int_t LG_TARGET_ADC_Thr[256] = {0};
  Int_t HG_SFT_ADC_Thr[128] = {0};
  Int_t LG_SFT_ADC_Thr[128] = {0};

  Int_t TDC_min_SFT;
  Int_t TDC_max_SFT;

  bool TARGET_High_has_neighbours[256] = {false};

  Int_t ADC_High_corr_max = 0;

  Int_t ADC_cut_TARGET2 = 850;

  Int_t adc_high_target[256];       Int_t ADC_High_TARGET[256];
  Int_t adc_low_target[256];        Int_t ADC_Low_TARGET[256];
  Int_t tdc_le_target[256][16];     Int_t TDC_LE_TARGET[256];
  Int_t tdc_te_target[256][16];
  Int_t TDC_LE_TARGET_corrected[256][6]={0};
  char TDC_LE_TARGET_corr[256][6][20];  // TO MOVE 

  Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];        Double_t ADC_High_SFT_corr[128];
  Int_t adc_low_sft[128];           Int_t ADC_Low_SFT[128];
  Int_t tdc_le_sft[128][16];        Int_t TDC_LE_SFT[128];
  Int_t tdc_te_sft[128][16];

  Int_t ADC_TOF1[24];

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

  Int_t MwpcADC[512];

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

  double TDC_diff = -1;
  int tdc_trigger[2][16];
  double tdc_ck_corr = 0.;

  Int_t has_TDC_SFT_hit[128] = {0};

  Bool_t Event_flag[40] = {false};
  Int_t Event_tag[14] = {999};

  //const float R_TARGET = 29.0;
  //const float R_TOF1 = 47.1;
  //const float R_SFT_L1 = 40.0;

  double Z_TOF1[12] = {999.99};

  const double Offset_TOF1[12] = {180., 410., 0., 0., 310., 280., -30., 70., 0., 0., 390., 430.};
  const double Slope_correction_TOF1[12] = {0.70, 0.78, 1., 1., 0.70, 0.60, 0.86, 0.77, 1., 1., 0.67, 0.78};

  char file_mapping[200];
  char par_finput[200];
  char par_finput2[200];

  Int_t par_temp[2][128];
  char par_temp2[512][50];

  vector <double> vec_xx_lepton;
  vector <double> vec_yy_lepton;
  vector <double> vec_ex_lepton;
  vector <double> vec_ey_lepton;

  vector <double> vec_xx_lepton_test;
  vector <double> vec_yy_lepton_test;

  vector <double> vec_xx_lepton_rotate;
  vector <double> vec_yy_lepton_rotate;
  vector <double> vec_ex_lepton_rotate;
  vector <double> vec_ey_lepton_rotate;

  vector <int> vec_lepton_bars;
  vector <int> vec_lepton_bars_final;
  vector <int> vec_lepton_bars_rotate;
  vector <int> vec_kaon_bars;

  vector <double> vec_lepton_size;
  vector <double> vec_lepton_rotate_size;
  vector <int> vec_bar;
  vector <int> vec_bar_rotate;
  vector <double> vec_yprime;
  vector <double> vec_yprime_rotate;
  vector <double> vec_Dy;
  vector <double> vec_Dy_rotate;

  double sumS = 999.;
  double sumS_rotate = 999.;

  vector <double> vec_xx_kaon;
  vector <double> vec_yy_kaon;
  vector <double> vec_ex_kaon;
  vector <double> vec_ey_kaon;

  vector <double> vec_Ck;
  vector <double> vec_Cpi;

  vector <double> vec_xx_TOF1_Marker;
  vector <double> vec_yy_TOF1_Marker;

  vector <double> vec_xx_TOF1;
  vector <double> vec_yy_TOF1;

  vector <double> vec_xx_TOF1_closest;
  vector <double> vec_yy_TOF1_closest;

  vector <double> vec_xx_int_TDC_Gap_Fibers;
  vector <double> vec_yy_int_TDC_Gap_Fibers;

  vector <double> vec_xx_int_TDC_Gap_Fibers_SFT;
  vector <double> vec_yy_int_TDC_Gap_Fibers_SFT;

  vector <double> vec_xx_int_TDC_TARGET;
  vector <double> vec_yy_int_TDC_TARGET;

  vector <double> vec_xx_kaon_stop;
  vector <double> vec_yy_kaon_stop;

  Int_t TDC_ck_selected[14] = {0};
  Int_t TDC_cpi_selected[14] = {0};

  Int_t TDC_ck_sum = 0;       double TDC_ck_avg = 0.;     double TDC_ck_sigma = 0.;
  Int_t TDC_cpi_sum = 0;      double TDC_cpi_avg = 0.;    double TDC_cpi_sigma = 0.;

  double TDC_ck_sigma2 = 0.;
  double TDC_cpi_sigma2 = 0.;

  Int_t TDC_ck_sum2 = 0;    double TDC_ck_avg2=0.;    int TDC_ck_counter = 0;
  Int_t TDC_cpi_sum2 = 0;   double TDC_cpi_avg2=0.;   int TDC_cpi_counter = 0;

  int tdc_vt48[256][16];
  vector<int> vec_tdc_b0_6;
  vector<int> vec_tdc_b0_7;

  vector<int> vec_TARGET_bar_selected;

  double Gap[12][3][2] = {{{0}}};

  bool Good_Event=false;

  bool Good_TARGET_Event = false;
  int count_TARGET_evts = 0;

  bool Good_TOF1_ADC[12]={false};   bool Good_TOF2_ADC[12]={false};
  bool Good_TOF1_TDC[12]={false};   bool Good_TOF2_TDC[12]={false};
  bool Good_TOF1[12]={false};       bool Good_TOF2[12]={false};
  bool Good_TOFs[12]={false};
  bool Good_TOF_Event = false;
  bool Event_On_Blacklist = false;
  string Current_Event;
  int current_event;

  int count_C2X = 0;    int count_C2Y = 0;
  int count_C3X = 0;    int count_C3Y = 0;
  int count_C4X = 0;    int count_C4Y = 0;
  bool Good_MWPC_Event = false;

  Int_t TDC_average = -1;

  int max_index_all[256];
  int max_ADC_all[256];
  int max_index_flag;

  int TDC_LG_max = -1;
  int TDC_LG_max2 = -1;
  int index_max1=0;
  int index_max2=0;

  int kaon_TDC_min;
  int kaon_TDC_max;
  //int TDC_min_Kstop;
  //int TDC_max_Kstop;

  Int_t TDC_min_TARGET;

  bool has_TDC_hit[256] = {false};
  bool has_TDC_hit_Kstop[256] = {false};

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
  bool has_both_ADC_TOF1_hit[12] = {false};
  bool has_both_TDC_TOF1_hit[12] = {false};

  int ADC_TOF1_hit[12] = {0};
  int ADCTDC_TOF1_hit[12] = {0};
  int ADC_TOF2_hit[12] = {0};
  int ADCTDC_TOF2_hit[12] = {0};

  int selected_TOF2 = 0;
  int gap_counter[12] = {0};

  int high_gap_hit = 0;
  int gap_to_fit = 0;
  int score_max = 0;

  vector<int> tof1_ties;

  bool k_stop_bar[256] = {false};
  bool k_stop_bar_initial[256] = {false};
  vector<int> good_k_stop_bars;

  TGraph *gr_kaon;
  TGraph *gr_kaon_bk;
  TF1 *gr_kaon_fit;

  double a_fit_kaon = 0.;
  double b_fit_kaon = 0.;
  double Chis_kaon = 0.;
  int ndf_kaon = 99;
  bool kaon_bk = false;

  TGraph *gr_lepton_1;
  TF1 *func_lepton_fit_1;

  double a_lepton_fit_1 = 0.;
  double b_lepton_fit_1 = 0.;
  double Chis_lepton_fit_1 = 0.;
  int ndf_lepton_fit_1 = 0;

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

  double C2X_centroid = 0.0;

  bool first_cluster = true;
  int cluster_spacing = 0;
  int cluster_length_count = 0;

  int C2X_clusters = 0;
  int C2Y_clusters = 0;
  int C3X_clusters = 0;
  int C3Y_clusters = 0;
  int C4X_clusters = 0;
  int C4Y_clusters = 0;

  vector<int> C2X_cluster_index;
  vector<int> C2X_cluster_length;
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

  int n_bar_selected=999;
  //char TDC_LE_TARGET_sel_string[n_bar_selected][16][20];

  Double_t par0_ADC = 0;  float par0_TDC = 0.;
  Double_t par1_ADC = 0;  float par1_TDC = 0.;

  double X_BAR = 999.;
  double Y_BAR = 999.;


  int gap_to_fit_rotate = 99;

  /// Data counters
  int has_data_TDC2 = 0;
  int has_data_ADC2 = 0;
  int has_data_ADC3 = 0;
  int has_data_ADC4 = 0;
  int has_data_ADCA = 0;

  ///// Select edge fiber for track fitting
  double Xloc_gap = 0;
  double Yloc_gap = 0;

  double xcoord = 0;
  Int_t unique_x = 0;

  double determinant;
  double x_circle_int1;
  double y_circle_int1;
  double x_circle_int2;
  double y_circle_int2;

  double SFTxdistance1;
  double SFTydistance1;

  double SFTxdistance2;
  double SFTydistance2;

  double SFTxhyp1;
  double SFTxhyp2;

  double SFT_x_intercept;
  double SFT_y_intercept;

  double SFT_phi;

  double x_intercept;
  double y_intercept;
  double x_distance4[256] = {0};
  double y_distance4[256] = {0};
  int distances4[256];

  double min_distance = 10000.0;

  double a_fit_GoodLG=0.;                     float b_fit_GoodLG=0.;
  double a_fit_GoodLG_weighted=0.;            float b_fit_GoodLG_weighted=0.;
  double a_fit_TDC_selected_weighted=0.;      float b_fit_TDC_selected_weighted=0.;

  int lepton_counter = 0;

  double a_lepton_fit_2 = 0.;
  double b_lepton_fit_2 = 0.;
  double Chis_lepton_fit_2 = 0.;
  int ndf_lepton_fit_2 = 0;

  TGraph *gr_lepton_2;
  TF1 *func_lepton_fit_2;
  TGraph *gr2_Leptons_rotate;
 
  double x_int_TDC_selected_weighted[2];         double y_int_TDC_selected_weighted[2];
  double x_int_GoodLG[2];                        double y_int_GoodLG[2];
  double x_int_GoodLG_weighted[2];               double y_int_GoodLG_weighted[2];
  double x_int_TDC_Gap_Fibers[2];                double y_int_TDC_Gap_Fibers[2];
  double x_int_TDC_Gap_Fibers_SFT[2];            double y_int_TDC_Gap_Fibers_SFT[2];
  double x_int_TARGET[2];                        double y_int_TARGET[2];

  double x_GoodLG_intersect1=0.;    double y_GoodLG_intersect1=0.;
  double x_TDC_Gap_Fibers=0.;       double y_TDC_Gap_Fibers=0.;

  double dist1_GoodLG[2];
  double dist1_TDC_Gap_Fibers[2];
  double dist1_TDC_Gap_Fibers_SFT[2];

  float dist2_TDC_selected[3];
  float dist2_GoodLG[3];

  float dist2_TDC_selected_min = 1000.;
  float dist2_GoodLG_min = 1000.;
  int selected_TDC_selected = 0;

  double ParError = 999.;
  double ChiS = 0.0;
  int ndf = 0;

  double a_lepton_fit_3 = 0.;
  double b_lepton_fit_3 = 0.;
  double Chis_lepton_fit_3 = 0.;
  int ndf_lepton_fit_3 = 0;

  TGraph *gr_lepton_3;
  TF1 *func_lepton_fit_3;
  TGraph *gr3_Leptons_rotate;

  double x_TDC_Gap_Fibers_intersect1=0.;          double y_TDC_Gap_Fibers_intersect1=0.;
  double x_TDC_Gap_Fibers_SFT_intersect1=0.;      double y_TDC_Gap_Fibers_SFT_intersect1=0.;
  double x_TARGET_intersect=0;                    double y_TARGET_intersect=0;
  double x_Arrows=0;                              double y_Arrows=0;

  double dist1_TARGET_intersect[2];

  double dist_to_k_stop = 0;

  double a_final_guide = 0.;
  double alpha_guide = 0.;
  double tanalpha_guide = 0.;
  double angle_final_guide = 0.;
  double Delta_phi = 999.99;  
  double Delta_phi_deg = 999.99;
  double Delta_X = 999.99;
  double Delta_Y = 999.99;


  double x_exit_rotate;
  double y_exit_rotate;

  double Yth_TOF1=999.;
  double New_ChiS=999.;
  double ChiS_cos=999.;
  double a_fit_lepton_rotate = 999.;
  double b_fit_lepton_rotate = 999.;

  double x_tof1_intersect_1 = 0;
  double y_tof1_intersect_1 = 0;
  double x_tof1_intersect_2 = 0;
  double y_tof1_intersect_2 = 0;
  double x_tof1_intersect = 0;
  double y_tof1_intersect = 0;

  double alpha;

  double X_weights = 0.;
  double Y_weights = 0.;
  double total_energy = 0.;
  vector<double> vec_kaon_centroid_coordinates;
  vector<double> vec_fit_lines_intersect;
  vector<double> vec_k_stop_coordinates;

  int i_kaon_bar = 999;
  double length_in_target = 0.;

  double sum_TDC_lepton = 0.;
  double sum_TDC_kaon = 0.;

  int counter_TDC_lepton = 0.;
  int counter_TDC_kaon = 0.;

  double Average_TDC_lepton = 0.;
  double Average_TDC_kaon = 0.;

  double distance_min_lepton = 999.;
  int lepton_out_bar = 999;

  int sum_ADC_HG_lepton = 0;

