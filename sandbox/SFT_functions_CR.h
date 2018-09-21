// Functions to use with SFT hit Z-values

using namespace std;
#include <vector>



double calculate_Z_position(double phi = 90., int N = 7 , double p = 4.0, int layer = 4){
  double Z_position = 0.0;

  double PI = 3.1415;
  
  phi *= PI/180.0;

  // Z_0 CONSTANT OFFSET FOR EACH LAYER
  /////////////////////////////////////////////////////////////
  //const double z_offset_const[4] = {-131.0, -136.5, -138.68, -135.83};
  const double z_offset_const[4] = {-125.298, -130.764, -138.13, -134.949};

  // Number of fibers in each ribbon/layer
  const int fibers_per_ribbon[4] = {15, 15, 17, 17};

  // Angle of wrapping
  const double alpha[4] = {3.57*(PI/180.0), 3.5*(PI/180.0), 3.87*(PI/180.0), 3.79*(PI/180.0)};
                  // OLD  3.6                 3.6               3.9                3.9
  // Radius of layer 
  const double layer_radius[4] = {80.5/2, 82.2/2, 84.2/2, 85.9/2};

  // Fiber widths and spacing between each fiber
  const double fiber_diameter = 1.0;
  const double fiber_spacing = 0.05;

  // horizontal width of each ribbon
  const double ribbon_width[4] = {15.75/cos(alpha[0]), 15.75/cos(alpha[1]), 
                                  17.85/cos(alpha[2]), 17.85/cos(alpha[3])};


  // Calculate the z position
  if(layer == 1){
    Z_position = z_offset_const[0] + layer_radius[0]*phi*tan(alpha[0]) 
               + (N-1)*ribbon_width[0]
               + (p-1.0)*ribbon_width[0]/fibers_per_ribbon[0]
               + 0.5*fiber_diameter/cos(alpha[0]);   
  }
  else if(layer == 2){
    Z_position = z_offset_const[1] + layer_radius[1]*phi*tan(alpha[1]) 
               + (N-1)*ribbon_width[1]
               + (p-1.0)*ribbon_width[1]/fibers_per_ribbon[1]
               + 0.5*fiber_diameter/cos(alpha[1]);     
  }
  else if(layer == 3){
    Z_position = z_offset_const[2] + layer_radius[2]*(2*PI-phi)*tan(alpha[2]) 
               + (N-1)*ribbon_width[2]
               + (fibers_per_ribbon[2]-p)*ribbon_width[2]/fibers_per_ribbon[2]
               + 0.5*fiber_diameter/cos(alpha[2]);   
  }
  else if(layer == 4){
    Z_position = z_offset_const[3] + layer_radius[3]*(2*PI-phi)*tan(alpha[3]) 
               + (N-1)*ribbon_width[3]
               + (fibers_per_ribbon[3]-p)*ribbon_width[3]/fibers_per_ribbon[3]
               + 0.5*fiber_diameter/cos(alpha[3]);   
  }

  if(p <= 0)
    Z_position = -1;
  

  return Z_position;
}
                                          
double *calculate_min_Z_position(double phi, 
  vector<double> layer_1_p_fibers, vector<double> layer_2_p_fibers, vector<double> layer_3_p_fibers, vector<double> layer_4_p_fibers, 
  double C2_Z, double return_array[24]){

  const double C2_SFT_delta_Z_thr = 7.0;


  for(int i = 0; i<24; i++)
    return_array[i] = -1;

  // Set p1_p2 / p3_p4 flags to pass value
  return_array[14] = 2;
  return_array[15] = 2;

  // Array to store all Z positions
  double Z_positions[4][20] = {0};

  // Number of turns per layer
  const int num_turns[4] = {17, 17, 15, 15};

  // Number of fibers in each ribbon/layer
  const int fibers_per_ribbon[4] = {15, 15, 17, 17};  


  bool no_p1_flag = false;
  bool no_p2_flag = false;
  bool no_p3_flag = false;
  bool no_p4_flag = false;

/////////////////////////////////////////////////////////////////////
// ADD NEW SFT Z ALGORITHM HERE.

  vector<double> good_p1_fiber;
  vector<double> good_p2_fiber;
  vector<double> good_p3_fiber;
  vector<double> good_p4_fiber;

  vector<int> p1_turn_number;
  vector<int> p2_turn_number;
  vector<int> p3_turn_number;
  vector<int> p4_turn_number;

  vector<double> p1_z_value;
  vector<double> p2_z_value;
  vector<double> p3_z_value;
  vector<double> p4_z_value;

  vector<double> C2_p1_delta_z;
  vector<double> C2_p2_delta_z;
  vector<double> C2_p3_delta_z;
  vector<double> C2_p4_delta_z;

  double current_sft_z = 0;

  vector<double> z_values;

  int count = 0;


  for(vector<double>::iterator it = layer_1_p_fibers.begin(); it != layer_1_p_fibers.end(); it++){
    for(int N=0; N<num_turns[0]; N++){
      current_sft_z = calculate_Z_position(phi,N,*it,1);

      if(fabs(current_sft_z - C2_Z) < C2_SFT_delta_Z_thr){
        good_p1_fiber.push_back(*it);
        p1_turn_number.push_back(N);
        p1_z_value.push_back(current_sft_z);
        C2_p1_delta_z.push_back(C2_Z-current_sft_z);
      }
    }  
  }


  for(vector<double>::iterator it = layer_2_p_fibers.begin(); it != layer_2_p_fibers.end(); it++){
    for(int N=0; N<num_turns[1]; N++){
      current_sft_z = calculate_Z_position(phi,N,*it,2);

      if(fabs(current_sft_z - C2_Z) < C2_SFT_delta_Z_thr){
        good_p2_fiber.push_back(*it);
        p2_turn_number.push_back(N);
        p2_z_value.push_back(current_sft_z);
        C2_p2_delta_z.push_back(C2_Z-current_sft_z);
      }
    }
  }

  for(vector<double>::iterator it = layer_3_p_fibers.begin(); it != layer_3_p_fibers.end(); it++){
    for(int N=0; N<num_turns[2]; N++){
      current_sft_z = calculate_Z_position(phi,N,*it,3);

      if(fabs(current_sft_z - C2_Z) < C2_SFT_delta_Z_thr){
        good_p3_fiber.push_back(*it);
        p3_turn_number.push_back(N);
        p3_z_value.push_back(current_sft_z);
        C2_p3_delta_z.push_back(C2_Z-current_sft_z);
      }
    }
  }

  for(vector<double>::iterator it = layer_4_p_fibers.begin(); it != layer_4_p_fibers.end(); it++){
    for(int N=0; N<num_turns[3]; N++){
      current_sft_z = calculate_Z_position(phi,N,*it,4);

      if(fabs(current_sft_z - C2_Z) < C2_SFT_delta_Z_thr){
        good_p4_fiber.push_back(*it);
        p4_turn_number.push_back(N);
        p4_z_value.push_back(current_sft_z);
        C2_p4_delta_z.push_back(C2_Z-current_sft_z);
      }
    }
  } 


  // Select best p-value (corresponding to smallest C2_px_delta_z)
  int vector_lengths = C2_p1_delta_z.size();
  double smallest_delta_z = 100;
  int best_z_value_index_1 = -1;

  int good_p_fiber_count = 0;

  if(vector_lengths != 0){
    for(int i=0; i<vector_lengths; i++){
      if(fabs(C2_p1_delta_z[i]) < smallest_delta_z){
        smallest_delta_z = fabs(C2_p1_delta_z[i]);
        best_z_value_index_1 = i;
      }
    }
    good_p_fiber_count++;
  }


  vector_lengths = C2_p2_delta_z.size();
  smallest_delta_z = 100;
  int best_z_value_index_2 = -1;  

  if(vector_lengths != 0){
    for(int i=0; i<vector_lengths; i++){
      if(fabs(C2_p2_delta_z[i]) < smallest_delta_z){
        smallest_delta_z = fabs(C2_p2_delta_z[i]);
        best_z_value_index_2 = i;
      }
    }
    good_p_fiber_count++;
  }


  vector_lengths = C2_p3_delta_z.size();
  smallest_delta_z = 100;
  int best_z_value_index_3 = -1;    

  if(vector_lengths != 0){
    for(int i=0; i<vector_lengths; i++){
      if(fabs(C2_p3_delta_z[i]) < smallest_delta_z){
        smallest_delta_z = fabs(C2_p3_delta_z[i]);
        best_z_value_index_3 = i;
      }
    }
    good_p_fiber_count++;
  }


  vector_lengths = C2_p4_delta_z.size();
  smallest_delta_z = 100;
  int best_z_value_index_4 = -1;    

  if(vector_lengths != 0){
    for(int i=0; i<vector_lengths; i++){
      if(fabs(C2_p4_delta_z[i]) < smallest_delta_z){
        smallest_delta_z = fabs(C2_p4_delta_z[i]);
        best_z_value_index_4 = i;
      }
    }
    good_p_fiber_count++;
  }      


  double z_avg_numerator = 0.0;
  double z_avg_denominator = 0.0;

  if(best_z_value_index_1 != -1){
    z_avg_numerator += p1_z_value[best_z_value_index_1];
    z_avg_denominator++;
    z_values.push_back(p1_z_value[best_z_value_index_1]);
    return_array[2] = p1_turn_number[best_z_value_index_1];
    return_array[6] = good_p1_fiber[best_z_value_index_1];
    return_array[10] = p1_z_value[best_z_value_index_1];
  }

  if(best_z_value_index_2 != -1){
    z_avg_numerator += p2_z_value[best_z_value_index_2];
    z_avg_denominator++;
    z_values.push_back(p2_z_value[best_z_value_index_2]);
    return_array[3] = p2_turn_number[best_z_value_index_2];
    return_array[7] = good_p2_fiber[best_z_value_index_2];
    return_array[11] = p2_z_value[best_z_value_index_2];
  }

  if(best_z_value_index_3 != -1){
    z_avg_numerator += p3_z_value[best_z_value_index_3];
    z_avg_denominator++;
    z_values.push_back(p3_z_value[best_z_value_index_3]);
    return_array[4] = p3_turn_number[best_z_value_index_3];
    return_array[8] = good_p3_fiber[best_z_value_index_3];
    return_array[12] = p3_z_value[best_z_value_index_3];
  }

  if(best_z_value_index_4 != -1){
    z_avg_numerator += p4_z_value[best_z_value_index_4];
    z_avg_denominator++;
    z_values.push_back(p4_z_value[best_z_value_index_4]);
    return_array[5] = p4_turn_number[best_z_value_index_4];
    return_array[9] = good_p4_fiber[best_z_value_index_4];
    return_array[13] = p4_z_value[best_z_value_index_4];
  }    

  return_array[0] = z_avg_numerator/z_avg_denominator;



  double largest_z_val = -1000;
  double smallest_z_val = 1000;

  for(vector<double>::iterator it = z_values.begin(); it != z_values.end(); it++){
    if(*it > largest_z_val)
      largest_z_val = *it;
    if(*it < smallest_z_val)
      smallest_z_val = *it;
  }

  return_array[1] = largest_z_val - smallest_z_val;



  // z_avgs_top = *return_array_top;
  // delta_z_top = *(return_array_top + 1);
  // N1_store_top = *(return_array_top + 2);
  // N2_store_top = *(return_array_top + 3);
  // N3_store_top = *(return_array_top + 4);
  // N4_store_top = *(return_array_top + 5);
  // p1_store_top = *(return_array_top + 6);
  // p2_store_top = *(return_array_top + 7);
  // p3_store_top = *(return_array_top + 8);
  // p4_store_top = *(return_array_top + 9);
  // z1_store_top = *(return_array_top + 10);
  // z2_store_top = *(return_array_top + 11);
  // z3_store_top = *(return_array_top + 12);
  // z4_store_top = *(return_array_top + 13);
  // p1_p2_flag_top = *(return_array_top + 14);
  // p3_p4_flag_top = *(return_array_top + 15);
  // z1_z2_flag_top = *(return_array_top + 16);
  // z3_z4_flag_top = *(return_array_top + 17);
  // z1234_flag_top = *(return_array_top + 18);

  return return_array;
}

vector<double> Layer_averaging(vector<int> p_fibers, bool *bad_sft_flag, int layer){
  vector<double> p_fibers_d;

  int layer_length;

  if(layer==1 || layer==2)
    layer_length = 15;
  else
    layer_length = 17;

  int p_hits[layer_length];
  for(int i = 0; i<layer_length; i++){
    p_hits[i]= -1;
  }

  for(vector<int>::iterator it = p_fibers.begin(); it != p_fibers.end(); it++){
    p_hits[*it - 1] = 1;
  }


  int p_cluster_separation = 0;

  // Find clustering of p fibers
  ////////////////////////////
  bool first_cluster = true;
  int cluster_spacing = 0;
  int cluster_length_count = 0;

  vector<int> p_cluster_index; // Hold the starting index of each cluster.
  vector<int> p_cluster_length; // Hold the length of each cluster. 


  for(int i=0; i<layer_length; i++){
    if(i!=0){
      if(p_hits[i]>0 && p_hits[i-1] <= 0){
        cluster_length_count++;
        p_cluster_index.push_back(i);
      }
      else if(p_hits[i]>0 && p_hits[i-1] > 0){
        cluster_length_count++;
      }
      else if(p_hits[i]<0 && p_hits[i-1]>0){
        p_cluster_length.push_back(cluster_length_count);
        cluster_length_count = 0;
      }
      else if(p_hits[i]<0 && p_hits[i-1] < 0){
        cluster_length_count = 0;
      }
    }
    else{
      if(p_hits[i]>0){
        cluster_length_count++;
        p_cluster_index.push_back(i);
      }
    }

    if(cluster_length_count!=0)
      p_cluster_length.push_back(cluster_length_count);
  }




  // Calculate and display centroids p fibers

  const int max_length = 100;

  int p_cluster_index_arr[max_length];
  int p_cluster_length_arr[max_length];

  double p_centroids[max_length];


  for(int i = 0; i< max_length; i++){
    p_cluster_index_arr[i] = -10000;

    p_cluster_length_arr[i] = -1;

    p_centroids[i] = -10000.0;  
  } 


  // Store all index/lengths in arrays.
  for(unsigned int i = 0; i < p_cluster_index.size(); i++){
    p_cluster_index_arr[i] = *(p_cluster_index.begin() + i);
  }

  for(unsigned int i = 0; i < p_cluster_length.size(); i++){
    p_cluster_length_arr[i] = *(p_cluster_length.begin() + i);
  }

  // Calculate all centroids for C2-all
  
  int centroid_sum = 0;
  int centroid_den = 0;


  for(unsigned int i = 0; i < 100; i++){
    if(p_cluster_index_arr[i] < -1000)
      break;

    p_centroids[i] = p_cluster_index[i] +(p_cluster_length[i] - 1)/2.0 + 1;
    p_fibers_d.push_back(p_centroids[i]);
  }  






  // for(vector<int>::iterator it = p_fibers.begin(); it != p_fibers.end(); it++){
  //   p_fibers_d.push_back(double(*it));
  // }

  // uncomment for flagging bad sft 
  // *bad_sft_flag = true;

  return p_fibers_d;
}


void SFT_print(double ADC_High_SFT_corr[128], int has_TDC_SFT_hit[128], int SFT_channel_to_fiber[128], 
  int evt, double phi_top, double phi_bottom, double C2_top, double C2_bottom, double SFT_Array[2], bool to_print){

  const int sft_hit_count_max = 3;
  const int NUM_ITERATIONS = 100;

  float R_TARGET = 29.0;
  float R_TOF1 = 47.1;
  float R_SFT_L1 = 40.0;
  float R_C2 = 629.4;

  int L1_sft_count = 0;
  int L2_sft_count = 0;
  int L3_sft_count = 0;
  int L4_sft_count = 0;


  vector<int> layer_1_p_fibers;
  vector<int> layer_2_p_fibers;
  vector<int> layer_3_p_fibers;
  vector<int> layer_4_p_fibers;

  vector<double> layer_1_adc;
  vector<double> layer_2_adc;
  vector<double> layer_3_adc;
  vector<double> layer_4_adc;


  bool L1_sft_flag = false;
  bool L2_sft_flag = false;
  bool L3_sft_flag = false;
  bool L4_sft_flag = false;


  for(int i = 0; i <128; i++){
    if(ADC_High_SFT_corr[i] != 0 && has_TDC_SFT_hit[i]){
      if(SFT_channel_to_layer[i] == 1){

        L1_sft_flag = false;

        for(vector<int>::iterator it = layer_1_p_fibers.begin(); it != layer_1_p_fibers.end(); it++){
          if(*it == SFT_channel_to_fiber[i]){
            L1_sft_flag = true;
            break;
          }
        }

        if(!L1_sft_flag){  
          layer_1_p_fibers.push_back(SFT_channel_to_fiber[i]);
          layer_1_adc.push_back(ADC_High_SFT_corr[i]);
          L1_sft_count++;
        }


      }
      else if(SFT_channel_to_layer[i] == 2){

        L2_sft_flag = false;

        for(vector<int>::iterator it = layer_2_p_fibers.begin(); it != layer_2_p_fibers.end(); it++){
          if(*it == SFT_channel_to_fiber[i]){
            L2_sft_flag = true;
            break;
          }
        }

        if(!L2_sft_flag){  
          layer_2_p_fibers.push_back(SFT_channel_to_fiber[i]);
          layer_2_adc.push_back(ADC_High_SFT_corr[i]);
          L2_sft_count++;
        }

      }
      else if(SFT_channel_to_layer[i] == 3){

        L3_sft_flag = false;

        for(vector<int>::iterator it = layer_3_p_fibers.begin(); it != layer_3_p_fibers.end(); it++){
          if(*it == SFT_channel_to_fiber[i]){
            L3_sft_flag = true;
            break;
          }
        }

        if(!L3_sft_flag){  
          layer_3_p_fibers.push_back(SFT_channel_to_fiber[i]);
          layer_3_adc.push_back(ADC_High_SFT_corr[i]);
          L3_sft_count++;
        }

      }
      else if(SFT_channel_to_layer[i] == 4){

        L4_sft_flag = false;

        for(vector<int>::iterator it = layer_4_p_fibers.begin(); it != layer_4_p_fibers.end(); it++){
          if(*it == SFT_channel_to_fiber[i]){
            L4_sft_flag = true;
            break;
          }
        }

        if(!L4_sft_flag){  
          layer_4_p_fibers.push_back(SFT_channel_to_fiber[i]);
          layer_4_adc.push_back(ADC_High_SFT_corr[i]);
          L4_sft_count++;
        }
      }
    }
  }   


  // pass vectors to new p-averaging function. return vector<double> of averaged p-fibers.
  // vector<double> Layer_averaging(vector<int> p_fibers)
  vector<double> layer_1_p_fibers_d;
  vector<double> layer_2_p_fibers_d;
  vector<double> layer_3_p_fibers_d;
  vector<double> layer_4_p_fibers_d;
  
  bool bad_sft = false;
  bool *bad_sft_flag = &bad_sft;

  layer_1_p_fibers_d = Layer_averaging(layer_1_p_fibers, bad_sft_flag, 1);
  layer_2_p_fibers_d = Layer_averaging(layer_2_p_fibers, bad_sft_flag, 2);
  layer_3_p_fibers_d = Layer_averaging(layer_3_p_fibers, bad_sft_flag, 3);
  layer_4_p_fibers_d = Layer_averaging(layer_4_p_fibers, bad_sft_flag, 4);

     

  // Test if more than 1 layer is empty or if any layer has more than 3 p fibers
  if((layer_1_p_fibers_d.size() == 0 && layer_2_p_fibers_d.size() == 0) ||  
   (layer_3_p_fibers_d.size() == 0 && layer_4_p_fibers_d.size() == 0)){
    if(to_print)
      cout << "Need 1 hit in each helicity." << endl;

    SFT_Array[0] = {-1000};
    SFT_Array[1] = {-1000};
    return;
  }

  bool flag = false;
  if(layer_1_p_fibers_d.size() > sft_hit_count_max){
    if(to_print)
      cout << "More than " << sft_hit_count_max << " SFT hits in layer 1" << endl;
    flag = true;
  }
  if(layer_2_p_fibers_d.size() > sft_hit_count_max){
    if(to_print)
      cout << "More than " << sft_hit_count_max << " SFT hits in layer 2" << endl;
    flag = true;
  }
  if(layer_3_p_fibers_d.size() > sft_hit_count_max){
    if(to_print)
      cout << "More than " << sft_hit_count_max << " SFT hits in layer 3" << endl;
    flag = true;
  }
  if(layer_4_p_fibers_d.size() > sft_hit_count_max){
    if(to_print)
      cout << "More than " << sft_hit_count_max << " SFT hits in layer 4" << endl;
    flag = true;
  }  

  if(flag){
    SFT_Array[0] = {-1000};
    SFT_Array[1] = {-1000};
    return;
  }

  int iteration_counter = 0;

  // Storage vectors
  double z_avgs_top;
  double delta_z_top;
  double p1_store_top = -1.0;
  double p2_store_top = -1.0;
  double p3_store_top = -1.0;
  double p4_store_top = -1.0;

  int N1_store_top;
  int N2_store_top;
  int N3_store_top;
  int N4_store_top;

  double z1_store_top;
  double z2_store_top;
  double z3_store_top;
  double z4_store_top;  

  int p1_p2_flag_top;
  int p3_p4_flag_top;

  int z1_z2_flag_top;
  int z3_z4_flag_top;
  int z1234_flag_top;


  double z_avgs_bottom;
  double delta_z_bottom;
  double p1_store_bottom = -1.0;
  double p2_store_bottom = -1.0;
  double p3_store_bottom = -1.0;
  double p4_store_bottom = -1.0;

  int N1_store_bottom;
  int N2_store_bottom;
  int N3_store_bottom;
  int N4_store_bottom;

  double z1_store_bottom;
  double z2_store_bottom;
  double z3_store_bottom;
  double z4_store_bottom;  

  int p1_p2_flag_bottom;
  int p3_p4_flag_bottom;

  int z1_z2_flag_bottom;
  int z3_z4_flag_bottom;
  int z1234_flag_bottom;



  double smallest_delta_z[NUM_ITERATIONS] = {0};

  double *return_array_top;
  double *return_array_bottom;

  // Case 1: all layers have atleast one hit
  double array_to_return_top[24] = {0};
  double array_to_return_bottom[24] = {0};

  // flag_array[14] == p1-p2 flag
  // flag_array[15] == p3-p4 flag
  // flag_array[16] == z12 out of threshold
  // flag_array[17] == z34 out of thresholf
  // flag_array[18] == z1234 out of threshold

  return_array_top = calculate_min_Z_position(phi_top, layer_1_p_fibers_d,layer_2_p_fibers_d, layer_3_p_fibers_d, layer_4_p_fibers_d, 
     C2_top, array_to_return_top);


  z_avgs_top = *return_array_top;
  delta_z_top = *(return_array_top + 1);
  N1_store_top = *(return_array_top + 2);
  N2_store_top = *(return_array_top + 3);
  N3_store_top = *(return_array_top + 4);
  N4_store_top = *(return_array_top + 5);
  p1_store_top = *(return_array_top + 6);
  p2_store_top = *(return_array_top + 7);
  p3_store_top = *(return_array_top + 8);
  p4_store_top = *(return_array_top + 9);
  z1_store_top = *(return_array_top + 10);
  z2_store_top = *(return_array_top + 11);
  z3_store_top = *(return_array_top + 12);
  z4_store_top = *(return_array_top + 13);
  p1_p2_flag_top = *(return_array_top + 14);
  p3_p4_flag_top = *(return_array_top + 15);
  z1_z2_flag_top = *(return_array_top + 16);
  z3_z4_flag_top = *(return_array_top + 17);
  z1234_flag_top = *(return_array_top + 18);


  int it_1 = 0; 
  int it_2 = 0; 
  int it_3 = 0; 
  int it_4 = 0; 

  // Remove p-fibers used in top z value;
  for(vector<double>::iterator it = layer_1_p_fibers_d.begin(); it != layer_1_p_fibers_d.end(); it++){
    if(fabs(*it - p1_store_top) < 0.1){
      it_1++;
      //layer_1_p_fibers_d.erase(it);
    }
  }
  for(vector<double>::iterator it = layer_2_p_fibers_d.begin(); it != layer_2_p_fibers_d.end(); it++){
    if(fabs(*it - p2_store_top) < 0.1){
      it_2++;
      //layer_2_p_fibers_d.erase(it);
    }
  }
  for(vector<double>::iterator it = layer_3_p_fibers_d.begin(); it != layer_3_p_fibers_d.end(); it++){
    if(fabs(*it - p3_store_top) < 0.1){
      it_3++;
      //layer_3_p_fibers_d.erase(it);
    }
  }
  for(vector<double>::iterator it = layer_4_p_fibers_d.begin(); it != layer_4_p_fibers_d.end(); it++){
    if(fabs(*it - p4_store_top) < 0.1){
      it_4++;
      //layer_4_p_fibers_d.erase(it);
    }
  }

  // cout << "T1" << endl;
  // layer_1_p_fibers_d.erase(layer_1_p_fibers_d.begin() + it_1);
  // cout << "T2" << endl;
  // layer_2_p_fibers_d.erase(layer_2_p_fibers_d.begin() + it_2);
  // cout << "T3" << endl;
  // layer_3_p_fibers_d.erase(layer_3_p_fibers_d.begin() + it_3);
  // cout << "T4" << endl;
  // layer_4_p_fibers_d.erase(layer_4_p_fibers_d.begin() + it_4);
  // cout << "T5" << endl;


  return_array_bottom = calculate_min_Z_position(phi_bottom, layer_1_p_fibers_d,layer_2_p_fibers_d, layer_3_p_fibers_d, layer_4_p_fibers_d, 
    C2_bottom, array_to_return_bottom);


  z_avgs_bottom = *return_array_bottom;
  delta_z_bottom = *(return_array_bottom + 1);
  N1_store_bottom = *(return_array_bottom + 2);
  N2_store_bottom = *(return_array_bottom + 3);
  N3_store_bottom = *(return_array_bottom + 4);
  N4_store_bottom = *(return_array_bottom + 5);
  p1_store_bottom = *(return_array_bottom + 6);
  p2_store_bottom = *(return_array_bottom + 7);
  p3_store_bottom = *(return_array_bottom + 8);
  p4_store_bottom = *(return_array_bottom + 9);
  z1_store_bottom = *(return_array_bottom + 10);
  z2_store_bottom = *(return_array_bottom + 11);
  z3_store_bottom = *(return_array_bottom + 12);
  z4_store_bottom = *(return_array_bottom + 13);
  p1_p2_flag_bottom = *(return_array_bottom + 14);
  p3_p4_flag_bottom = *(return_array_bottom + 15);
  z1_z2_flag_bottom = *(return_array_bottom + 16);
  z3_z4_flag_bottom = *(return_array_bottom + 17);
  z1234_flag_bottom = *(return_array_bottom + 18);

  SFT_Array[0] = z_avgs_top;
  SFT_Array[1] = z_avgs_bottom;

  //cout << SFT_Array[0] << endl;


  if(to_print)
    cout<<endl;

  if(to_print){
    cout << "Event #     Angle      L1        L2        L3        L4    ||   Z Avg      Delta Z" << endl;
    cout << "===========================================================||======================" << endl;
  }
    
  // cout << "p1_store_top = " << p1_store_top << endl;
  // cout << "p2_store_top = " << p2_store_top << endl;
  // cout << "p3_store_top = " << p3_store_top << endl;
  // cout << "p4_store_top = " << p4_store_top << endl;
  // cout << "p1_store_bottom = " << p1_store_bottom << endl;
  // cout << "p2_store_bottom = " << p2_store_bottom << endl;
  // cout << "p3_store_bottom = " << p3_store_bottom << endl;
  // cout << "p4_store_bottom = " << p4_store_bottom << endl;


  if(to_print){
    if(p1_store_top >0 && p2_store_top > 0 && p3_store_top > 0 && p4_store_top > 0){
      printf("%-5d   %10.2f   %5.1f     %5.1f     %5.1f     %5.1f   ||\n", evt, phi_top, p1_store_top,p2_store_top,
        p3_store_top, p4_store_top);
      printf("N(i)               %5d     %5d     %5d     %5d     ||\n", N1_store_top,N2_store_top,N3_store_top,N4_store_top);
      printf("Z(i)               %7.2f   %7.2f    %7.2f   %7.2f  || %7.2f      %5.3f\n", z1_store_top,z2_store_top,z3_store_top,z4_store_top,
       z_avgs_top,delta_z_top);

      cout << "===========================================================||====================== " << endl;
    }
    else if(p1_store_top <=0 && p2_store_top > 0 && p3_store_top > 0 && p4_store_top > 0){
      printf("%-5d   %10.2f %5c       %5.1f     %5.1f     %5.1f   ||\n", evt, phi_top, '-',p2_store_top,
        p3_store_top, p4_store_top);
      printf("N(i)               %5c     %5d     %5d     %5d     ||\n", '-',N2_store_top,N3_store_top,N4_store_top);
      printf("Z(i)             %7c     %7.2f    %7.2f   %7.2f  || %7.2f      %5.1f \n", '-',z2_store_top,z3_store_top,z4_store_top,
       z_avgs_top,delta_z_top);

      cout << "===========================================================||====================== " << endl; 
    }
    else if(p1_store_top >0 && p2_store_top <= 0 && p3_store_top > 0 && p4_store_top > 0){
      printf("%-5d   %10.2f   %5.1f   %5c       %5.1f     %5.1f   ||\n", evt, phi_top, p1_store_top,'-',
        p3_store_top, p4_store_top);
      printf("N(i)               %5d     %5c     %5d     %5d     ||\n", N1_store_top,'-',N3_store_top,N4_store_top);
      printf("Z(i)               %7.2f %7c      %7.2f   %7.2f  || %7.2f      %5.1f     \n", z1_store_top,'-',z3_store_top,z4_store_top,
       z_avgs_top,delta_z_top);

      cout << "===========================================================||====================== " << endl;
    }
    else if(p1_store_top >0 && p2_store_top > 0 && p3_store_top <= 0 && p4_store_top > 0){
      printf("%-5d   %10.2f   %5.1f     %5.1f   %5c       %5.1f   ||\n", evt, phi_top, p1_store_top,p2_store_top,
        '-', p4_store_top);
      printf("N(i)               %5d     %5d     %5c     %5d     ||\n", N1_store_top,N2_store_top,'-',N4_store_top);
      printf("Z(i)               %7.2f   %7.2f %7c      %7.2f  || %7.2f      %5.1f     \n", z1_store_top,z2_store_top,'-',z4_store_top,
       z_avgs_top,delta_z_top);

      cout << "===========================================================||====================== " << endl;
    }
    else if(p1_store_top >0 && p2_store_top > 0 && p3_store_top > 0 && p4_store_top <= 0){
      printf("%-5d   %10.2f   %5.1f     %5.1f     %5.1f   %5c     ||\n", evt, phi_top, p1_store_top,p2_store_top,
        p3_store_top, '-');
      printf("N(i)               %5d     %5d     %5d     %5c     ||\n", N1_store_top,N2_store_top,N3_store_top,'-');
      printf("Z(i)               %7.2f   %7.2f    %7.2f%7c     || %7.2f      %5.1f        \n", z1_store_top,z2_store_top,z3_store_top,'-',
       z_avgs_top,delta_z_top);

      cout << "===========================================================||====================== " << endl;
    }
    else if(p1_store_top <= 0 && p2_store_top >0 && p3_store_top <= 0 && p4_store_top >0){
      printf("%-5d   %10.2f %5c       %5.1f   %5c       %5.1f   ||\n", evt, phi_top, '-',p2_store_top,
        '-', p4_store_top);
      printf("N(i)               %5c     %5d     %5c     %5d     ||\n",'-',N2_store_top,'-',N4_store_top);
      printf("Z(i)             %7c     %7.2f %7c      %7.2f  || %7.2f      %5.1f     \n", '-',z2_store_top,'-',z4_store_top,
       z_avgs_top,delta_z_top);

      cout << "===========================================================||====================== " << endl;
    }
    else if(p1_store_top <= 0 && p2_store_top >0 && p3_store_top >0 && p4_store_top <= 0){
      printf("%-5d   %10.2f %5c       %5.1f     %5.1f   %5c     ||\n", evt, phi_top, '-',p2_store_top,
        p3_store_top, '-');
      printf("N(i)               %5c     %5d     %5d     %5c     ||\n", '-',N2_store_top,N3_store_top,'-');
      printf("Z(i)             %7c     %7.2f    %7.2f%7c     || %7.2f      %5.1f        \n", '-',z2_store_top,z3_store_top,'-',
       z_avgs_top,delta_z_top);

      cout << "===========================================================||====================== " << endl;
    }
    else if(p1_store_top >0 && p2_store_top <= 0 && p3_store_top <= 0 && p4_store_top >0){
      printf("%-5d   %10.2f   %5.1f   %5c     %5c       %5.1f   ||\n", evt, phi_top, p1_store_top,'-',
        '-', p4_store_top);
      printf("N(i)               %5d     %5c     %5c     %5d     ||\n", N1_store_top,'-','-',N4_store_top);
      printf("Z(i)               %7.2f %7c   %7c      %7.2f  || %7.2f      %5.1f    \n", z1_store_top,'-','-',z4_store_top,
       z_avgs_top,delta_z_top);

      cout << "===========================================================||====================== " << endl;
    }
    else if(p1_store_top >0 && p2_store_top <= 0 && p3_store_top >0 && p4_store_top <= 0){
      printf("%-5d   %10.2f   %5.1f   %5c       %5.1f   %5c     ||\n", evt, phi_top, p1_store_top,'-',
        p3_store_top, '-');
      printf("N(i)               %5d     %5c     %5d     %5c     ||\n", N1_store_top,'-',N3_store_top,'-');
      printf("Z(i)               %7.2f %7c      %7.2f%7c     || %7.2f      %5.1f     \n", z1_store_top,'-',z3_store_top,'-',
       z_avgs_top,delta_z_top);

      cout << "===========================================================||====================== " << endl;
    }


    if(p1_store_bottom >0 && p2_store_bottom > 0 && p3_store_bottom > 0 && p4_store_bottom > 0){
      printf("%-5d   %10.2f   %5.1f     %5.1f     %5.1f     %5.1f   ||\n", evt, phi_bottom, p1_store_bottom,p2_store_bottom,
        p3_store_bottom, p4_store_bottom);
      printf("N(i)               %5d     %5d     %5d     %5d     ||\n", N1_store_bottom,N2_store_bottom,N3_store_bottom,N4_store_bottom);
      printf("Z(i)               %7.2f   %7.2f    %7.2f   %7.2f  || %7.2f      %5.3f\n", z1_store_bottom,z2_store_bottom,z3_store_bottom,z4_store_bottom,
       z_avgs_bottom,delta_z_bottom);

      cout << "===========================================================||====================== " << endl;
    }
    else if(p1_store_bottom <=0 && p2_store_bottom > 0 && p3_store_bottom > 0 && p4_store_bottom > 0){
      printf("%-5d   %10.2f %5c       %5.1f     %5.1f     %5.1f   ||\n", evt, phi_bottom, '-',p2_store_bottom,
        p3_store_bottom, p4_store_bottom);
      printf("N(i)               %5c     %5d     %5d     %5d     ||\n", '-',N2_store_bottom,N3_store_bottom,N4_store_bottom);
      printf("Z(i)             %7c     %7.2f    %7.2f   %7.2f  || %7.2f      %5.1f \n", '-',z2_store_bottom,z3_store_bottom,z4_store_bottom,
       z_avgs_bottom,delta_z_bottom);

      cout << "===========================================================||====================== " << endl; 
    }
    else if(p1_store_bottom >0 && p2_store_bottom <= 0 && p3_store_bottom > 0 && p4_store_bottom > 0){
      printf("%-5d   %10.2f   %5.1f   %5c       %5.1f     %5.1f   ||\n", evt, phi_bottom, p1_store_bottom,'-',
        p3_store_bottom, p4_store_bottom);
      printf("N(i)               %5d     %5c     %5d     %5d     ||\n", N1_store_bottom,'-',N3_store_bottom,N4_store_bottom);
      printf("Z(i)               %7.2f %7c      %7.2f   %7.2f  || %7.2f      %5.1f     \n", z1_store_bottom,'-',z3_store_bottom,z4_store_bottom,
       z_avgs_bottom,delta_z_bottom);

      cout << "===========================================================||====================== " << endl;
    }
    else if(p1_store_bottom >0 && p2_store_bottom > 0 && p3_store_bottom <= 0 && p4_store_bottom > 0){
      printf("%-5d   %10.2f   %5.1f     %5.1f   %5c       %5.1f   ||\n", evt, phi_bottom, p1_store_bottom,p2_store_bottom,
        '-', p4_store_bottom);
      printf("N(i)               %5d     %5d     %5c     %5d     ||\n", N1_store_bottom,N2_store_bottom,'-',N4_store_bottom);
      printf("Z(i)               %7.2f   %7.2f %7c      %7.2f  || %7.2f      %5.1f     \n", z1_store_bottom,z2_store_bottom,'-',z4_store_bottom,
       z_avgs_bottom,delta_z_bottom);

      cout << "===========================================================||====================== " << endl;
    }
    else if(p1_store_bottom >0 && p2_store_bottom > 0 && p3_store_bottom > 0 && p4_store_bottom <= 0){
      printf("%-5d   %10.2f   %5.1f     %5.1f     %5.1f   %5c     ||\n", evt, phi_bottom, p1_store_bottom,p2_store_bottom,
        p3_store_bottom, '-');
      printf("N(i)               %5d     %5d     %5d     %5c     ||\n", N1_store_bottom,N2_store_bottom,N3_store_bottom,'-');
      printf("Z(i)               %7.2f   %7.2f    %7.2f%7c     || %7.2f      %5.1f        \n", z1_store_bottom,z2_store_bottom,z3_store_bottom,'-',
       z_avgs_bottom,delta_z_bottom);

      cout << "===========================================================||====================== " << endl;
    }
    else if(p1_store_bottom <= 0 && p2_store_bottom >0 && p3_store_bottom <= 0 && p4_store_bottom >0){
      printf("%-5d   %10.2f %5c       %5.1f   %5c       %5.1f   ||\n", evt, phi_bottom, '-',p2_store_bottom,
        '-', p4_store_bottom);
      printf("N(i)               %5c     %5d     %5c     %5d     ||\n",'-',N2_store_bottom,'-',N4_store_bottom);
      printf("Z(i)             %7c     %7.2f %7c      %7.2f  || %7.2f      %5.1f     \n", '-',z2_store_bottom,'-',z4_store_bottom,
       z_avgs_bottom,delta_z_bottom);

      cout << "===========================================================||====================== " << endl;
    }
    else if(p1_store_bottom <= 0 && p2_store_bottom >0 && p3_store_bottom >0 && p4_store_bottom <= 0){
      printf("%-5d   %10.2f %5c       %5.1f     %5.1f   %5c     ||\n", evt, phi_bottom, '-',p2_store_bottom,
        p3_store_bottom, '-');
      printf("N(i)               %5c     %5d     %5d     %5c     ||\n", '-',N2_store_bottom,N3_store_bottom,'-');
      printf("Z(i)             %7c     %7.2f    %7.2f%7c     || %7.2f      %5.1f        \n", '-',z2_store_bottom,z3_store_bottom,'-',
       z_avgs_bottom,delta_z_bottom);

      cout << "===========================================================||====================== " << endl;
    }
    else if(p1_store_bottom >0 && p2_store_bottom <= 0 && p3_store_bottom <= 0 && p4_store_bottom >0){
      printf("%-5d   %10.2f   %5.1f   %5c     %5c       %5.1f   ||\n", evt, phi_bottom, p1_store_bottom,'-',
        '-', p4_store_bottom);
      printf("N(i)               %5d     %5c     %5c     %5d     ||\n", N1_store_bottom,'-','-',N4_store_bottom);
      printf("Z(i)               %7.2f %7c   %7c      %7.2f  || %7.2f      %5.1f    \n", z1_store_bottom,'-','-',z4_store_bottom,
       z_avgs_bottom,delta_z_bottom);

      cout << "===========================================================||====================== " << endl;
    }
    else if(p1_store_bottom >0 && p2_store_bottom <= 0 && p3_store_bottom >0 && p4_store_bottom <= 0){
      printf("%-5d   %10.2f   %5.1f   %5c       %5.1f   %5c     ||\n", evt, phi_bottom, p1_store_bottom,'-',
        p3_store_bottom, '-');
      printf("N(i)               %5d     %5c     %5d     %5c     ||\n", N1_store_bottom,'-',N3_store_bottom,'-');
      printf("Z(i)               %7.2f %7c      %7.2f%7c     || %7.2f      %5.1f     \n", z1_store_bottom,'-',z3_store_bottom,'-',
       z_avgs_bottom,delta_z_bottom);

      cout << "===========================================================||====================== " << endl;
    }    
  }



  //return sft_z_selected;

  return;
}   