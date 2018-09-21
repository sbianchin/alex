// Functions to use with SFT hit Z-values



double calculate_Z_position(double phi = 90., int N = 7 , double p = 4.0, int layer = 4){
  double Z_position = 0.0;

  double PI = 3.1415;
  
  phi *= PI/180.0;

  // Z_0 CONSTANT OFFSET FOR EACH LAYER
  /////////////////////////////////////////////////////////////
  //const double z_offset_const[4] = {-131.0, -136.5, -138.68, -135.83};  old
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
               + fiber_diameter/cos(alpha[0]);   
  }
  else if(layer == 2){
    Z_position = z_offset_const[1] + layer_radius[1]*phi*tan(alpha[1]) 
               + (N-1)*ribbon_width[1]
               + (p-1.0)*ribbon_width[1]/fibers_per_ribbon[1]
               + fiber_diameter/cos(alpha[1]);     
  }
  else if(layer == 3){
    Z_position = z_offset_const[2] + layer_radius[2]*(2*PI-phi)*tan(alpha[2]) 
               + (N-1)*ribbon_width[2]
               + (fibers_per_ribbon[2]-p)*ribbon_width[2]/fibers_per_ribbon[2]
               + fiber_diameter/cos(alpha[2]); 
  }
  else if(layer == 4){
    Z_position = z_offset_const[3] + layer_radius[3]*(2*PI-phi)*tan(alpha[3]) 
               + (N-1)*ribbon_width[3]
               + (fibers_per_ribbon[3]-p)*ribbon_width[3]/fibers_per_ribbon[3]
               + fiber_diameter/cos(alpha[3]);
  }

  if(p <= 0)
    Z_position = -1;
  

  return Z_position;
}

double *calculate_min_Z_position(double phi, double p1, double p2, double p3, double p4, double return_array[20], int iteration_counter, double delta_z[100]){


  for(int i = 0; i<20; i++)
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

  // selected p's in an array
  const double sft_p_fiber[4] = {p1, p2, p3, p4};


  // Variables for calculating Z position from all four layers.
  double smallest_dist = 100.0;
  int layer_1_low_index = 0;
  int layer_2_low_index = 0;
  int layer_3_low_index = 0;
  int layer_4_low_index = 0;
  double n[6] = {0.0};
  double smallest_dist_temp = 0.;


  if(p1 > 0 && p2 > 0 && p3 > 0 && p4 > 0){ // Case 1: all p-fibers have a value

    //test for closeness of p1-p2, p3-p4
    //cout << "    p1 ********* " << p1 << endl;

    // p2-p1 = 6, so
    // Reject if p2-p1 != 4,5,6,7
    if(p2-p1 >= 0){
      if(p2-p1 < 4 || p2-p1 > 7){
        //Reject
        return_array[14] = 0;
      }
    }
    else{ // Case if p1 < p2, occurs when p1 and p2 are in different turn numbers.
      if((p2+15)-p1 < 4 || (p2+15)-p1 > 7){
        // Reject
        return_array[14] = 0;
      }
    }

    // p3-p4 = 1,2,3
    if(p4-p3 >= 0){
      if(p4-p3 < 1 || p4-p3 > 4){
        //Reject
        return_array[15] = 0;
      }
    }
    else{ // Case if p1 < p2, occurs when p1 and p2 are in different turn numbers.
      if((p4+17)-p3 < 1 || (p4+17)-p3 > 4){
        // Reject
        return_array[15] = 0;
      }
    }  
  }
  else if(p1 < 0 && p2 > 0 && p3 > 0 && p4 > 0){ // Case 2: p1 has no value

    //test for closeness of p1-p2, p3-p4

    return_array[14] = 2; // 2 == flag for no value in fiber
 

    // p3-p4 = 1,2,3
    if(p4-p3 >= 0){
      if(p4-p3 < 1 || p4-p3 > 4){
        //Reject
        return_array[15] = 0;
      }
    }
    else{ // Case if p1 < p2, occurs when p1 and p2 are in different turn numbers.
      if((p4+17)-p3 < 1 || (p4+17)-p3 > 4){
        // Reject
        return_array[15] = 0;
      }
    }  
  }
  else if(p1 > 0 && p2 < 0 && p3 > 0 && p4 > 0){ // Case 3: p2 has no value

    //test for closeness of p1-p2, p3-p4

    return_array[14] = 2; // 2 == flag for no value in fiber
 

    // p3-p4 = 1,2,3
    if(p4-p3 >= 0){
      if(p4-p3 < 1 || p4-p3 > 4){
        //Reject
        return_array[15] = 0;
      }
    }
    else{ // Case if p1 < p2, occurs when p1 and p2 are in different turn numbers.
      if((p4+17)-p3 < 1 || (p4+17)-p3 > 4){
        // Reject
        return_array[15] = 0;
      }
    }  
  }  
  else if(p1 > 0 && p2 > 0 && p3 < 0 && p4 > 0){ // Case 4: p3 has no value

    //test for closeness of p1-p2, p3-p4

    // p2-p1 = 6, so
    // Reject if p2-p1 != 5,6,7
    if(p2-p1 >= 0){
      if(p2-p1 < 4 || p2-p1 > 7){
        //Reject
        return_array[14] = 0;
      }
    }
    else{ // Case if p1 < p2, occurs when p1 and p2 are in different turn numbers.
      if((p2+15)-p1 < 4 || (p2+15)-p1 > 7){
        // Reject
        return_array[14] = 0;
      }
    }

    return_array[15] = 2;
  }  
  else if(p1 > 0 && p2 > 0 && p3 > 0 && p4 < 0){ // Case 4: p4 has no value

    //test for closeness of p1-p2, p3-p4

    // p2-p1 = 6, so
    // Reject if p2-p1 != 5,6,7
    if(p2-p1 >= 0){
      if(p2-p1 < 4 || p2-p1 > 7){
        //Reject
        return_array[14] = 0;
      }
    }
    else{ // Case if p1 < p2, occurs when p1 and p2 are in different turn numbers.
      if((p2+15)-p1 < 4 || (p2+15)-p1 > 7){
        // Reject
        return_array[14] = 0;
      }
    }

    return_array[15] = 2;
  }  
  else{ //one hit in each helicity

    return_array[14] = 2;
    return_array[15] = 2;
  }


  // Find Z positions for layer 1
  for(int N = 0; N < num_turns[0]; N++){        
    Z_positions[0][N] = calculate_Z_position(phi, N, p1, 1);
  }

  //Find Z positions for layer 2;
  for(int N = 0; N < num_turns[1]; N++){
    Z_positions[1][N] = calculate_Z_position(phi, N, p2, 2);
  }

  // Find Z positions for layer 3;
  for(int N = 0; N < num_turns[2]; N++){
    Z_positions[2][N] = calculate_Z_position(phi, N, p3, 3); 
  }

  // Find Z positions for layer 4;
  for(int N = 0; N < num_turns[3]; N++){
    Z_positions[3][N] = calculate_Z_position(phi, N, p4, 4);
  }


  bool no_p1_flag = false;
  bool no_p2_flag = false;
  bool no_p3_flag = false;
  bool no_p4_flag = false;

  bool same_turn_flag = false;

  const double smallest_dist_thr = 0.001;


  for(int h = 0; h<num_turns[0]; h++){
    if(p1 <= 0)
      no_p1_flag=true;

    for(int i = 0; i<num_turns[1]; i++){
      if(p2 <= 0)
        no_p2_flag=true;

      for(int j = 0; j<num_turns[2]; j++){
        if(p3 <= 0)
          no_p3_flag=true;

        for(int k = 0; k<num_turns[3]; k++){
          if(p4 <= 0)
            no_p4_flag=true;

            if(!(h==i && i==j && j==k)){

              n[0] = fabs(Z_positions[0][h] - Z_positions[1][i]);
              n[1] = fabs(Z_positions[0][h] - Z_positions[2][j]);
              n[2] = fabs(Z_positions[0][h] - Z_positions[3][k]);
              n[3] = fabs(Z_positions[1][i] - Z_positions[2][j]);
              n[4] = fabs(Z_positions[1][i] - Z_positions[3][k]);
              n[5] = fabs(Z_positions[2][j] - Z_positions[3][k]);



              // Determine smallest distance
              ///////////////////////////////
              ///////////////////////////////

              if(!no_p1_flag && !no_p2_flag && !no_p3_flag && !no_p4_flag){  // has p1,p2,p3,p4
                if(n[0] < smallest_dist &&
                   n[1] < smallest_dist && 
                   n[2] < smallest_dist && 
                   n[3] < smallest_dist && 
                   n[4] < smallest_dist && 
                   n[5] < smallest_dist){
                                                        ///////////////////// remove layer_1_low_index, etc
                   smallest_dist_temp = 0.; 

                   for(int t = 0; t<6; t++){
                      if(n[t] > smallest_dist_temp)
                        smallest_dist_temp = n[t];
                   }

    ///////////////////////////////////////////////////////////////////////////////////////////////////              
                   int smallest_dist_flag = false;    
                                                             
                   for(int i_sd = 0; i_sd < iteration_counter; i_sd++){ 
                     if(fabs(smallest_dist_temp - delta_z[i_sd]) < smallest_dist_thr) 
                       smallest_dist_flag = true;             
                   }                                          
                                                                
                   if(!smallest_dist_flag){                       
                      smallest_dist = smallest_dist_temp;     
                      layer_1_low_index = h;                 
                      layer_2_low_index = i;                
                      layer_3_low_index = j;                 
                      layer_4_low_index = k;             
                                                                     
                      return_array[0] = (Z_positions[0][layer_1_low_index] 
                                       + Z_positions[1][layer_2_low_index] 
                                       + Z_positions[2][layer_3_low_index] 
                                       + Z_positions[3][layer_4_low_index])/4.0; 
                   }

   ///////////////////////////////////////////////////////////////////////////////////////////////////
                }               
              }
              else if(no_p1_flag && !no_p2_flag && !no_p3_flag && !no_p4_flag){ // has p2,p3,p4
                if(n[3] < smallest_dist && 
                   n[4] < smallest_dist && 
                   n[5] < smallest_dist){

                  smallest_dist_temp = 0.; 


                  if(n[3] > smallest_dist_temp)
                    smallest_dist_temp = n[3];
                  if(n[4] > smallest_dist_temp)
                    smallest_dist_temp = n[4];
                  if(n[5] > smallest_dist_temp)
                    smallest_dist_temp = n[5];


                   int smallest_dist_flag = false;    
                                                             
                   for(int i_sd = 0; i_sd < iteration_counter; i_sd++){ 
                     if(fabs(smallest_dist_temp - delta_z[i_sd]) < smallest_dist_thr) 
                       smallest_dist_flag = true;             
                   }                                          
                                                                
                   if(!smallest_dist_flag){                       
                      smallest_dist = smallest_dist_temp;     
                      layer_1_low_index = -1;                 
                      layer_2_low_index = i;                
                      layer_3_low_index = j;                 
                      layer_4_low_index = k;             
                                                                     
                      return_array[0] = (Z_positions[1][layer_2_low_index] 
                                  + Z_positions[2][layer_3_low_index] 
                                  + Z_positions[3][layer_4_low_index])/3.0;  
                   }
                }                  
              }
              else if(!no_p1_flag && no_p2_flag && !no_p3_flag && !no_p4_flag){ // has p1,p3,p4
                if(n[1] < smallest_dist && 
                   n[2] < smallest_dist && 
                   n[5] < smallest_dist){

                   
                   smallest_dist_temp = 0.; 

                   if(n[1] > smallest_dist_temp)
                     smallest_dist_temp = n[1];
                   if(n[2] > smallest_dist_temp)
                     smallest_dist_temp = n[2];
                   if(n[5] > smallest_dist_temp)
                     smallest_dist_temp = n[5];
           

                   int smallest_dist_flag = false;    
                                                             
                   for(int i_sd = 0; i_sd < iteration_counter; i_sd++){ 
                     if(fabs(smallest_dist_temp - delta_z[i_sd]) < smallest_dist_thr) 
                       smallest_dist_flag = true;             
                   }                                          
                                                                
                   if(!smallest_dist_flag){                       
                      smallest_dist = smallest_dist_temp;     
                      layer_1_low_index = h;
                      layer_2_low_index = -1;
                      layer_3_low_index = j;
                      layer_4_low_index = k;             
                                                                     
                      return_array[0] = (Z_positions[0][layer_1_low_index] 
                                    + Z_positions[2][layer_3_low_index] 
                                    + Z_positions[3][layer_4_low_index])/3.0; 
                   }
                } 
              }
              else if(!no_p1_flag && !no_p2_flag && no_p3_flag && !no_p4_flag){ // has p1,p2,p4
                if(n[0] < smallest_dist &&
                   n[2] < smallest_dist && 
                   n[4] < smallest_dist){

                   
                   
                   smallest_dist_temp = 0.; 

                   if(n[0] > smallest_dist_temp)
                     smallest_dist_temp = n[0];
                   if(n[2] > smallest_dist_temp)
                     smallest_dist_temp = n[2];
                   if(n[4] > smallest_dist_temp)
                     smallest_dist_temp = n[4];

                  int smallest_dist_flag = false;    
                                                             
                   for(int i_sd = 0; i_sd < iteration_counter; i_sd++){ 
                     if(fabs(smallest_dist_temp - delta_z[i_sd]) < smallest_dist_thr) 
                       smallest_dist_flag = true;             
                   }                                          
                                                                
                   if(!smallest_dist_flag){                       
                      smallest_dist = smallest_dist_temp;     
                      layer_1_low_index = h;
                      layer_2_low_index = i;
                      layer_3_low_index = -1;
                      layer_4_low_index = k;          
                                                                     
                      return_array[0] = (Z_positions[0][layer_1_low_index] 
                                    + Z_positions[1][layer_2_low_index] 
                                    + Z_positions[3][layer_4_low_index])/3.0; 
                   }


                } 
              }
              else if(!no_p1_flag && !no_p2_flag && !no_p3_flag && no_p4_flag){ // has p1,p2,p3
                if(n[0] < smallest_dist &&
                   n[1] < smallest_dist && 
                   n[3] < smallest_dist){

                   layer_1_low_index = h;
                   layer_2_low_index = i;
                   layer_3_low_index = j;
                   layer_4_low_index = -1;
                   
                   smallest_dist_temp = 0.; 

                   if(n[0] > smallest_dist_temp)
                     smallest_dist_temp = n[0];
                   if(n[1] > smallest_dist_temp)
                     smallest_dist_temp = n[1];
                   if(n[3] > smallest_dist_temp)
                     smallest_dist_temp = n[3];


                   int smallest_dist_flag = false;    
                                                             
                   for(int i_sd = 0; i_sd < iteration_counter; i_sd++){ 
                     if(fabs(smallest_dist_temp - delta_z[i_sd]) < smallest_dist_thr) 
                       smallest_dist_flag = true;             
                   }                                          
                                                                
                   if(!smallest_dist_flag){                       
                      smallest_dist = smallest_dist_temp;     
                      layer_1_low_index = h;
                      layer_2_low_index = i;
                      layer_3_low_index = j;
                      layer_4_low_index = -1;          
                                                                     
                      return_array[0] = (Z_positions[0][layer_1_low_index] 
                                    + Z_positions[1][layer_2_low_index] 
                                    + Z_positions[2][layer_3_low_index])/3.0; 
                   }

                }
              }
              else if(no_p1_flag && !no_p2_flag && no_p3_flag && !no_p4_flag){ // has p2,p4
                if(n[4] < smallest_dist){

                  smallest_dist_temp = n[4];

                    int smallest_dist_flag = false;    
                                                             
                   for(int i_sd = 0; i_sd < iteration_counter; i_sd++){ 
                     if(fabs(smallest_dist_temp - delta_z[i_sd]) < smallest_dist_thr) 
                       smallest_dist_flag = true;             
                   }                                          
                                                                
                   if(!smallest_dist_flag){                       
                      smallest_dist = smallest_dist_temp;     
                      layer_1_low_index = -1;
                      layer_2_low_index = i;
                      layer_3_low_index = -1;
                      layer_4_low_index = k;           
                                                                     
                      return_array[0] = (Z_positions[1][layer_2_low_index] 
                                    + Z_positions[3][layer_4_low_index])/2.0;
                   }

                }   
              }
              else if(no_p1_flag && !no_p2_flag && !no_p3_flag && no_p4_flag){ // has p2,p3
                if(n[3] < smallest_dist){

                  smallest_dist_temp = n[3];

                   int smallest_dist_flag = false;    
                                                             
                   for(int i_sd = 0; i_sd < iteration_counter; i_sd++){ 
                     if(fabs(smallest_dist_temp - delta_z[i_sd]) < smallest_dist_thr) 
                       smallest_dist_flag = true;             
                   }                                          
                                                                
                   if(!smallest_dist_flag){                       
                      smallest_dist = smallest_dist_temp;     
                      layer_1_low_index = -1;
                      layer_2_low_index = i;
                      layer_3_low_index = j;
                      layer_4_low_index = -1;            
                                                                     
                      return_array[0] = (Z_positions[1][layer_2_low_index] 
                                    + Z_positions[2][layer_3_low_index])/2.0;  
                   }
                }
              }
              else if(!no_p1_flag && no_p2_flag && no_p3_flag && !no_p4_flag){ // has p1,p4
                if(n[2] < smallest_dist){

                 smallest_dist_temp = n[2];

                   
                  int smallest_dist_flag = false;    
                                                             
                   for(int i_sd = 0; i_sd < iteration_counter; i_sd++){ 
                     if(fabs(smallest_dist_temp - delta_z[i_sd]) < smallest_dist_thr) 
                       smallest_dist_flag = true;             
                   }                                          
                                                                
                   if(!smallest_dist_flag){                       
                      smallest_dist = smallest_dist_temp;     
                      layer_1_low_index = h;
                      layer_2_low_index = -1;
                      layer_3_low_index = -1;
                      layer_4_low_index = k;             
                                                                     
                      return_array[0] = (Z_positions[0][layer_1_low_index] 
                                    + Z_positions[3][layer_4_low_index])/2.0;  
                   }


                }
              }
              else if(!no_p1_flag && no_p2_flag && !no_p3_flag && no_p4_flag){ // has p1,p3
                if(n[1] < smallest_dist){
                  smallest_dist_temp = n[1];
                   
                   int smallest_dist_flag = false;    
                                                             
                   for(int i_sd = 0; i_sd < iteration_counter; i_sd++){ 
                     if(fabs(smallest_dist_temp - delta_z[i_sd]) < smallest_dist_thr) 
                       smallest_dist_flag = true;             
                   }                                          
                                                                
                   if(!smallest_dist_flag){                       
                      smallest_dist = smallest_dist_temp;;     
                      layer_1_low_index = h;
                      layer_2_low_index = -1;
                      layer_3_low_index = j;
                      layer_4_low_index = -1;             
                                                                     
                      return_array[0] = (Z_positions[0][layer_1_low_index] 
                                    + Z_positions[3][layer_4_low_index])/2.0;  
                   }
                }
              }




              return_array[1] = smallest_dist;
              return_array[2] = layer_1_low_index;    
              return_array[3] = layer_2_low_index;    
              return_array[4] = layer_3_low_index;    
              return_array[5] = layer_4_low_index;             
              return_array[6] = p1;    
              return_array[7] = p2;    
              return_array[8] = p3;    
              return_array[9] = p4;    
              return_array[10] = Z_positions[0][layer_1_low_index];
              return_array[11] = Z_positions[1][layer_2_low_index];
              return_array[12] = Z_positions[2][layer_3_low_index];
              return_array[13] = Z_positions[3][layer_4_low_index];
          }
        }
      }
    }
  }


  double z12_thr = 5.0;
  double z34_thr = 5.0;
  double z_total_thr = 5.0;

  if(fabs(return_array[10] - return_array[11]) > z12_thr){
    return_array[16] = 0;
  }
  if(fabs(return_array[12] - return_array[13]) > z34_thr){
    return_array[17] = 0;
  }


  // This removes the flags on the p2-p1 and p4-p3 conditions.
   // return_array[14] = 2;
   // return_array[15] = 2;

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

  cout << "Layer " << layer << ": ";
  for(unsigned int i = 0; i < 100; i++){
    if(p_cluster_index_arr[i] < -1000)
      break;

    p_centroids[i] = p_cluster_index[i] +(p_cluster_length[i] - 1)/2.0 + 1;
    cout << p_centroids[i] << ", ";
    p_fibers_d.push_back(p_centroids[i]);
  }  

  cout << endl;






  // for(vector<int>::iterator it = p_fibers.begin(); it != p_fibers.end(); it++){
  //   p_fibers_d.push_back(double(*it));
  // }

  // uncomment for flagging bad sft 
  // *bad_sft_flag = true;

  return p_fibers_d;
}


double SFT_print(double ADC_High_SFT_corr[128], int has_TDC_SFT_hit[128], int SFT_channel_to_fiber[128],
 int evt, double phi, bool to_print, double TOF1_pos, bool use_TOF1, double C2X_centroid, double len_in_target){

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
    return -10000;
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

  if(flag)
    return -10000; 

  int iteration_counter = 0;

  // Storage vectors
  vector<double> z_avgs_arr;
  vector<double> delta_z_arr;
  vector<double> p1_store_arr;
  vector<double> p2_store_arr;
  vector<double> p3_store_arr;
  vector<double> p4_store_arr;

  vector<int> N1_store_arr;
  vector<int> N2_store_arr;
  vector<int> N3_store_arr;
  vector<int> N4_store_arr;

  vector<double> z1_store_arr;
  vector<double> z2_store_arr;
  vector<double> z3_store_arr;
  vector<double> z4_store_arr;  

  vector<int> p1_p2_flag;
  vector<int> p3_p4_flag;

  vector<int> z1_z2_flag;
  vector<int> z3_z4_flag;
  vector<int> z1234_flag;

  double smallest_delta_z[NUM_ITERATIONS] = {0};

  double *return_array;

  // Case 1: all layers have atleast one hit
  if(layer_1_p_fibers_d.size() != 0 && layer_2_p_fibers_d.size() != 0 && layer_3_p_fibers_d.size() != 0 && layer_4_p_fibers_d.size() != 0){
    for(vector<double>::iterator it_1 = layer_1_p_fibers_d.begin(); it_1 != layer_1_p_fibers_d.end(); it_1++){
      for(vector<double>::iterator it_2 = layer_2_p_fibers_d.begin(); it_2 != layer_2_p_fibers_d.end(); it_2++){
        for(vector<double>::iterator it_3 = layer_3_p_fibers_d.begin(); it_3 != layer_3_p_fibers_d.end(); it_3++){
          for(vector<double>::iterator it_4 = layer_4_p_fibers_d.begin(); it_4 != layer_4_p_fibers_d.end(); it_4++){ // 

            double array_to_return[20] = {1};
            
            // flag_array[14] == p1-p2 flag
            // flag_array[15] == p3-p4 flag
            // flag_array[16] == z12 out of threshold
            // flag_array[17] == z34 out of thresholf
            // flag_array[18] == z1234 out of threshold

            for(int i=0; i<NUM_ITERATIONS; i++){
              return_array = calculate_min_Z_position(phi, *it_1, *it_2, *it_3, *it_4, array_to_return, i, smallest_delta_z);
              smallest_delta_z[i] = *(return_array+1);

              z_avgs_arr.push_back(*return_array);
              delta_z_arr.push_back(*(return_array + 1));
              N1_store_arr.push_back(*(return_array + 2));
              N2_store_arr.push_back(*(return_array + 3));
              N3_store_arr.push_back(*(return_array + 4));
              N4_store_arr.push_back(*(return_array + 5));
              p1_store_arr.push_back(*(return_array + 6));
              p2_store_arr.push_back(*(return_array + 7));
              p3_store_arr.push_back(*(return_array + 8));
              p4_store_arr.push_back(*(return_array + 9));
              z1_store_arr.push_back(*(return_array + 10));
              z2_store_arr.push_back(*(return_array + 11));
              z3_store_arr.push_back(*(return_array + 12));
              z4_store_arr.push_back(*(return_array + 13));
              p1_p2_flag.push_back(*(return_array + 14));
              p3_p4_flag.push_back(*(return_array + 15));
              z1_z2_flag.push_back(*(return_array + 16));
              z3_z4_flag.push_back(*(return_array + 17));
              z1234_flag.push_back(*(return_array + 18));

              iteration_counter++;
            }
          }
        }
      }
    }
  }

  // Case 2: layer 1 has 0 hits
  else if(layer_1_p_fibers_d.size() == 0 && layer_2_p_fibers_d.size() != 0 && layer_3_p_fibers_d.size() != 0 && layer_4_p_fibers_d.size() != 0){
    for(vector<double>::iterator it_2 = layer_2_p_fibers_d.begin(); it_2 != layer_2_p_fibers_d.end(); it_2++){
      for(vector<double>::iterator it_3 = layer_3_p_fibers_d.begin(); it_3 != layer_3_p_fibers_d.end(); it_3++){
        for(vector<double>::iterator it_4 = layer_4_p_fibers_d.begin(); it_4 != layer_4_p_fibers_d.end(); it_4++){ // 

          double array_to_return[20] = {1};
          double *return_array;
          // flag_array[14] == p1-p2 flag
          // flag_array[15] == p3-p4 flag
          // flag_array[16] == z12 out of threshold
          // flag_array[17] == z34 out of thresholf
          // flag_array[18] == z1234 out of threshold

          for(int i=0; i<NUM_ITERATIONS; i++){
              return_array = calculate_min_Z_position(phi, -1, *it_2, *it_3, *it_4, array_to_return, i, smallest_delta_z);
              smallest_delta_z[i] = *(return_array+1);

              z_avgs_arr.push_back(*return_array);
              delta_z_arr.push_back(*(return_array + 1));
              N1_store_arr.push_back(*(return_array + 2));
              N2_store_arr.push_back(*(return_array + 3));
              N3_store_arr.push_back(*(return_array + 4));
              N4_store_arr.push_back(*(return_array + 5));
              p1_store_arr.push_back(*(return_array + 6));
              p2_store_arr.push_back(*(return_array + 7));
              p3_store_arr.push_back(*(return_array + 8));
              p4_store_arr.push_back(*(return_array + 9));
              z1_store_arr.push_back(*(return_array + 10));
              z2_store_arr.push_back(*(return_array + 11));
              z3_store_arr.push_back(*(return_array + 12));
              z4_store_arr.push_back(*(return_array + 13));
              p1_p2_flag.push_back(*(return_array + 14));
              p3_p4_flag.push_back(*(return_array + 15));
              z1_z2_flag.push_back(*(return_array + 16));
              z3_z4_flag.push_back(*(return_array + 17));
              z1234_flag.push_back(*(return_array + 18));

              iteration_counter++;
            }
        }
      }
    }
  }  

  // Case 3: layer 2 has 0 hits
  else if(layer_1_p_fibers_d.size() != 0 && layer_2_p_fibers_d.size() == 0 && layer_3_p_fibers_d.size() != 0 && layer_4_p_fibers_d.size() != 0){
    for(vector<double>::iterator it_1 = layer_1_p_fibers_d.begin(); it_1 != layer_1_p_fibers_d.end(); it_1++){
      for(vector<double>::iterator it_3 = layer_3_p_fibers_d.begin(); it_3 != layer_3_p_fibers_d.end(); it_3++){
        for(vector<double>::iterator it_4 = layer_4_p_fibers_d.begin(); it_4 != layer_4_p_fibers_d.end(); it_4++){ // 

          double array_to_return[20] = {1};
          double *return_array;
          // flag_array[14] == p1-p2 flag
          // flag_array[15] == p3-p4 flag
          // flag_array[16] == z12 out of threshold
          // flag_array[17] == z34 out of thresholf
          // flag_array[18] == z1234 out of threshold

          for(int i=0; i<NUM_ITERATIONS; i++){
              return_array = calculate_min_Z_position(phi, *it_1, -1, *it_3, *it_4, array_to_return, i, smallest_delta_z);
              smallest_delta_z[i] = *(return_array+1);

              z_avgs_arr.push_back(*return_array);
              delta_z_arr.push_back(*(return_array + 1));
              N1_store_arr.push_back(*(return_array + 2));
              N2_store_arr.push_back(*(return_array + 3));
              N3_store_arr.push_back(*(return_array + 4));
              N4_store_arr.push_back(*(return_array + 5));
              p1_store_arr.push_back(*(return_array + 6));
              p2_store_arr.push_back(*(return_array + 7));
              p3_store_arr.push_back(*(return_array + 8));
              p4_store_arr.push_back(*(return_array + 9));
              z1_store_arr.push_back(*(return_array + 10));
              z2_store_arr.push_back(*(return_array + 11));
              z3_store_arr.push_back(*(return_array + 12));
              z4_store_arr.push_back(*(return_array + 13));
              p1_p2_flag.push_back(*(return_array + 14));
              p3_p4_flag.push_back(*(return_array + 15));
              z1_z2_flag.push_back(*(return_array + 16));
              z3_z4_flag.push_back(*(return_array + 17));
              z1234_flag.push_back(*(return_array + 18));

              iteration_counter++;
            }
        }
      }
    }
  }  

  // Case 4: layer 3 has 0 hits
  else if(layer_1_p_fibers_d.size() != 0 && layer_2_p_fibers_d.size() != 0 && layer_3_p_fibers_d.size() == 0 && layer_4_p_fibers_d.size() != 0){
    for(vector<double>::iterator it_1 = layer_1_p_fibers_d.begin(); it_1 != layer_1_p_fibers_d.end(); it_1++){
      for(vector<double>::iterator it_2 = layer_2_p_fibers_d.begin(); it_2 != layer_2_p_fibers_d.end(); it_2++){
        for(vector<double>::iterator it_4 = layer_4_p_fibers_d.begin(); it_4 != layer_4_p_fibers_d.end(); it_4++){ // 

          double array_to_return[20] = {1};
          double *return_array;
          // flag_array[14] == p1-p2 flag
          // flag_array[15] == p3-p4 flag
          // flag_array[16] == z12 out of threshold
          // flag_array[17] == z34 out of thresholf
          // flag_array[18] == z1234 out of threshold

          for(int i=0; i<NUM_ITERATIONS; i++){
              return_array = calculate_min_Z_position(phi, *it_1, *it_2, -1, *it_4, array_to_return, i, smallest_delta_z);
              smallest_delta_z[i] = *(return_array+1);

              z_avgs_arr.push_back(*return_array);
              delta_z_arr.push_back(*(return_array + 1));
              N1_store_arr.push_back(*(return_array + 2));
              N2_store_arr.push_back(*(return_array + 3));
              N3_store_arr.push_back(*(return_array + 4));
              N4_store_arr.push_back(*(return_array + 5));
              p1_store_arr.push_back(*(return_array + 6));
              p2_store_arr.push_back(*(return_array + 7));
              p3_store_arr.push_back(*(return_array + 8));
              p4_store_arr.push_back(*(return_array + 9));
              z1_store_arr.push_back(*(return_array + 10));
              z2_store_arr.push_back(*(return_array + 11));
              z3_store_arr.push_back(*(return_array + 12));
              z4_store_arr.push_back(*(return_array + 13));
              p1_p2_flag.push_back(*(return_array + 14));
              p3_p4_flag.push_back(*(return_array + 15));
              z1_z2_flag.push_back(*(return_array + 16));
              z3_z4_flag.push_back(*(return_array + 17));
              z1234_flag.push_back(*(return_array + 18));

              iteration_counter++;
            }
        }
      }
    }
  }  
  
  // Case 5: layer 4 has 0 hits
  else if(layer_1_p_fibers_d.size() != 0 && layer_2_p_fibers_d.size() != 0 && layer_3_p_fibers_d.size() != 0 && layer_4_p_fibers_d.size() == 0){
    for(vector<double>::iterator it_1 = layer_1_p_fibers_d.begin(); it_1 != layer_1_p_fibers_d.end(); it_1++){
      for(vector<double>::iterator it_2 = layer_2_p_fibers_d.begin(); it_2 != layer_2_p_fibers_d.end(); it_2++){
        for(vector<double>::iterator it_3 = layer_3_p_fibers_d.begin(); it_3 != layer_3_p_fibers_d.end(); it_3++){

          double array_to_return[20] = {1};
          double *return_array;
          // flag_array[14] == p1-p2 flag
          // flag_array[15] == p3-p4 flag
          // flag_array[16] == z12 out of threshold
          // flag_array[17] == z34 out of thresholf
          // flag_array[18] == z1234 out of threshold

          for(int i=0; i<NUM_ITERATIONS; i++){
              return_array = calculate_min_Z_position(phi, *it_1, *it_2, *it_3, -1, array_to_return, i, smallest_delta_z);
              smallest_delta_z[i] = *(return_array+1);

              z_avgs_arr.push_back(*return_array);
              delta_z_arr.push_back(*(return_array + 1));
              N1_store_arr.push_back(*(return_array + 2));
              N2_store_arr.push_back(*(return_array + 3));
              N3_store_arr.push_back(*(return_array + 4));
              N4_store_arr.push_back(*(return_array + 5));
              p1_store_arr.push_back(*(return_array + 6));
              p2_store_arr.push_back(*(return_array + 7));
              p3_store_arr.push_back(*(return_array + 8));
              p4_store_arr.push_back(*(return_array + 9));
              z1_store_arr.push_back(*(return_array + 10));
              z2_store_arr.push_back(*(return_array + 11));
              z3_store_arr.push_back(*(return_array + 12));
              z4_store_arr.push_back(*(return_array + 13));
              p1_p2_flag.push_back(*(return_array + 14));
              p3_p4_flag.push_back(*(return_array + 15));
              z1_z2_flag.push_back(*(return_array + 16));
              z3_z4_flag.push_back(*(return_array + 17));
              z1234_flag.push_back(*(return_array + 18));

              iteration_counter++;
            }
        }
      }
    }
  }  

  //        1) 0 1 0 1
  else if(layer_1_p_fibers_d.size() == 0 && layer_2_p_fibers_d.size() != 0 && layer_3_p_fibers_d.size() == 0 && layer_4_p_fibers_d.size() != 0){
    for(vector<double>::iterator it_2 = layer_2_p_fibers_d.begin(); it_2 != layer_2_p_fibers_d.end(); it_2++){
      for(vector<double>::iterator it_4 = layer_4_p_fibers_d.begin(); it_4 != layer_4_p_fibers_d.end(); it_4++){ // 

        double array_to_return[20] = {1};
        double *return_array;
        // flag_array[14] == p1-p2 flag
        // flag_array[15] == p3-p4 flag
        // flag_array[16] == z12 out of threshold
        // flag_array[17] == z34 out of thresholf
        // flag_array[18] == z1234 out of threshold

        for(int i=0; i<NUM_ITERATIONS; i++){
              return_array = calculate_min_Z_position(phi, -1, *it_2, -1, *it_4, array_to_return, i, smallest_delta_z);
              smallest_delta_z[i] = *(return_array+1);

              z_avgs_arr.push_back(*return_array);
              delta_z_arr.push_back(*(return_array + 1));
              N1_store_arr.push_back(*(return_array + 2));
              N2_store_arr.push_back(*(return_array + 3));
              N3_store_arr.push_back(*(return_array + 4));
              N4_store_arr.push_back(*(return_array + 5));
              p1_store_arr.push_back(*(return_array + 6));
              p2_store_arr.push_back(*(return_array + 7));
              p3_store_arr.push_back(*(return_array + 8));
              p4_store_arr.push_back(*(return_array + 9));
              z1_store_arr.push_back(*(return_array + 10));
              z2_store_arr.push_back(*(return_array + 11));
              z3_store_arr.push_back(*(return_array + 12));
              z4_store_arr.push_back(*(return_array + 13));
              p1_p2_flag.push_back(*(return_array + 14));
              p3_p4_flag.push_back(*(return_array + 15));
              z1_z2_flag.push_back(*(return_array + 16));
              z3_z4_flag.push_back(*(return_array + 17));
              z1234_flag.push_back(*(return_array + 18));

              iteration_counter++;
            }
      }
    }
  }

  else if(layer_1_p_fibers_d.size() == 0 && layer_2_p_fibers_d.size() != 0 && layer_3_p_fibers_d.size() != 0 && layer_4_p_fibers_d.size() == 0){
    for(vector<double>::iterator it_2 = layer_2_p_fibers_d.begin(); it_2 != layer_2_p_fibers_d.end(); it_2++){
      for(vector<double>::iterator it_3 = layer_3_p_fibers_d.begin(); it_3 != layer_3_p_fibers_d.end(); it_3++){
          
        double array_to_return[20] = {1};
        double *return_array;
        // flag_array[14] == p1-p2 flag
        // flag_array[15] == p3-p4 flag
        // flag_array[16] == z12 out of threshold
        // flag_array[17] == z34 out of thresholf
        // flag_array[18] == z1234 out of threshold

        for(int i=0; i<NUM_ITERATIONS; i++){
              return_array = calculate_min_Z_position(phi, -1, *it_2, *it_3, -1, array_to_return, i, smallest_delta_z);
              smallest_delta_z[i] = *(return_array+1);

              z_avgs_arr.push_back(*return_array);
              delta_z_arr.push_back(*(return_array + 1));
              N1_store_arr.push_back(*(return_array + 2));
              N2_store_arr.push_back(*(return_array + 3));
              N3_store_arr.push_back(*(return_array + 4));
              N4_store_arr.push_back(*(return_array + 5));
              p1_store_arr.push_back(*(return_array + 6));
              p2_store_arr.push_back(*(return_array + 7));
              p3_store_arr.push_back(*(return_array + 8));
              p4_store_arr.push_back(*(return_array + 9));
              z1_store_arr.push_back(*(return_array + 10));
              z2_store_arr.push_back(*(return_array + 11));
              z3_store_arr.push_back(*(return_array + 12));
              z4_store_arr.push_back(*(return_array + 13));
              p1_p2_flag.push_back(*(return_array + 14));
              p3_p4_flag.push_back(*(return_array + 15));
              z1_z2_flag.push_back(*(return_array + 16));
              z3_z4_flag.push_back(*(return_array + 17));
              z1234_flag.push_back(*(return_array + 18));

              iteration_counter++;
            }
      }
    }
  }

  else if(layer_1_p_fibers_d.size() != 0 && layer_2_p_fibers_d.size() == 0 && layer_3_p_fibers_d.size() == 0 && layer_4_p_fibers_d.size() != 0){
    for(vector<double>::iterator it_1 = layer_1_p_fibers_d.begin(); it_1 != layer_1_p_fibers_d.end(); it_1++){
      for(vector<double>::iterator it_4 = layer_4_p_fibers_d.begin(); it_4 != layer_4_p_fibers_d.end(); it_4++){ // 

        double array_to_return[20] = {1};
        double *return_array;
        // flag_array[14] == p1-p2 flag
        // flag_array[15] == p3-p4 flag
        // flag_array[16] == z12 out of threshold
        // flag_array[17] == z34 out of thresholf
        // flag_array[18] == z1234 out of threshold


        for(int i=0; i<NUM_ITERATIONS; i++){
              return_array = calculate_min_Z_position(phi, *it_1, -1, -1, *it_4, array_to_return, i, smallest_delta_z);
              smallest_delta_z[i] = *(return_array+1);

              z_avgs_arr.push_back(*return_array);
              delta_z_arr.push_back(*(return_array + 1));
              N1_store_arr.push_back(*(return_array + 2));
              N2_store_arr.push_back(*(return_array + 3));
              N3_store_arr.push_back(*(return_array + 4));
              N4_store_arr.push_back(*(return_array + 5));
              p1_store_arr.push_back(*(return_array + 6));
              p2_store_arr.push_back(*(return_array + 7));
              p3_store_arr.push_back(*(return_array + 8));
              p4_store_arr.push_back(*(return_array + 9));
              z1_store_arr.push_back(*(return_array + 10));
              z2_store_arr.push_back(*(return_array + 11));
              z3_store_arr.push_back(*(return_array + 12));
              z4_store_arr.push_back(*(return_array + 13));
              p1_p2_flag.push_back(*(return_array + 14));
              p3_p4_flag.push_back(*(return_array + 15));
              z1_z2_flag.push_back(*(return_array + 16));
              z3_z4_flag.push_back(*(return_array + 17));
              z1234_flag.push_back(*(return_array + 18));

              iteration_counter++;
            }
      }
    }
  }

  else if(layer_1_p_fibers_d.size() != 0 && layer_2_p_fibers_d.size() == 0 && layer_3_p_fibers_d.size() != 0 && layer_4_p_fibers_d.size() == 0){
    for(vector<double>::iterator it_1 = layer_1_p_fibers_d.begin(); it_1 != layer_1_p_fibers_d.end(); it_1++){
      for(vector<double>::iterator it_3 = layer_3_p_fibers_d.begin(); it_3 != layer_3_p_fibers_d.end(); it_3++){
      
        double array_to_return[20] = {1};
        double *return_array;
        // flag_array[14] == p1-p2 flag
        // flag_array[15] == p3-p4 flag
        // flag_array[16] == z12 out of threshold
        // flag_array[17] == z34 out of thresholf
        // flag_array[18] == z1234 out of threshold

        for(int i=0; i<NUM_ITERATIONS; i++){
              return_array = calculate_min_Z_position(phi, *it_1, -1, *it_3, -1, array_to_return, i, smallest_delta_z);
              smallest_delta_z[i] = *(return_array+1);

              z_avgs_arr.push_back(*return_array);
              delta_z_arr.push_back(*(return_array + 1));
              N1_store_arr.push_back(*(return_array + 2));
              N2_store_arr.push_back(*(return_array + 3));
              N3_store_arr.push_back(*(return_array + 4));
              N4_store_arr.push_back(*(return_array + 5));
              p1_store_arr.push_back(*(return_array + 6));
              p2_store_arr.push_back(*(return_array + 7));
              p3_store_arr.push_back(*(return_array + 8));
              p4_store_arr.push_back(*(return_array + 9));
              z1_store_arr.push_back(*(return_array + 10));
              z2_store_arr.push_back(*(return_array + 11));
              z3_store_arr.push_back(*(return_array + 12));
              z4_store_arr.push_back(*(return_array + 13));
              p1_p2_flag.push_back(*(return_array + 14));
              p3_p4_flag.push_back(*(return_array + 15));
              z1_z2_flag.push_back(*(return_array + 16));
              z3_z4_flag.push_back(*(return_array + 17));
              z1234_flag.push_back(*(return_array + 18));

              iteration_counter++;
            }
      }
    }
  }


  // Select z_avgs corresponding to delta z < thr

  double z_avg_min = 0;
  double delta_z_min = 1000;
  int delta_z_min_index = 0;

  vector<int> good_delta_z_index;

  vector<double> z_avgs_store;
  bool z_printed = false;


  double delta_z_thr = 5.0;

  double sft_delta_z_selected = 1000;
  double sft_z_selected = 0;

  for(int i = 0; i < iteration_counter; i++){
    if(delta_z_arr[i] < delta_z_thr){
      good_delta_z_index.push_back(i);
    }
  }

  if(to_print)
    cout<<endl;

  if(to_print){
    cout << "Event #     Angle      L1        L2        L3        L4    ||   Z Avg      Delta Z      Track Length" << endl;
    cout << "===========================================================||========================================" << endl;
  }

  for(vector<int>::iterator it = good_delta_z_index.begin(); it != good_delta_z_index.end(); it++){
  

    

    double z_avg = 0.0;
    int z_avg_count = 0;

    if(N1_store_arr[*it] != -1){
      z_avg += z1_store_arr[*it];
      z_avg_count++;
    }

    if(N2_store_arr[*it] != -1){
      z_avg += z2_store_arr[*it];
      z_avg_count++;
    }

    if(N3_store_arr[*it] != -1){
      z_avg += z3_store_arr[*it];
      z_avg_count++;
    }

    if(N4_store_arr[*it] != -1){
      z_avg += z4_store_arr[*it];
      z_avg_count++;
    }
    
    z_printed = false;


    double sft_angle = atan(fabs(C2X_centroid - z_avg/double(z_avg_count))/(R_C2-R_SFT_L1));
    double total_length = len_in_target/cos(sft_angle);

    if(delta_z_arr[*it] <= sft_delta_z_selected && fabs(z_avg/double(z_avg_count)) <=100){
      if(use_TOF1){
        if(fabs(z_avg/double(z_avg_count) - TOF1_pos) < TOF1_Z_cut){
          sft_z_selected = z_avg/double(z_avg_count);
          sft_delta_z_selected = delta_z_arr[*it];
        }
      }
      else{
        sft_z_selected = z_avg/double(z_avg_count);
        sft_delta_z_selected = delta_z_arr[*it];
      }
    }

    if(p1_p2_flag[*it] > 1 && p3_p4_flag[*it] > 1){
      if(z_avgs_store.empty()){
        z_printed = false;
      }
      else{
        for(vector<double>::iterator db = z_avgs_store.begin(); db != z_avgs_store.end(); db++){
          if(fabs(*db - z_avg/double(z_avg_count)) < 0.01){
            z_printed = true;
          }
        }
      }

      if(to_print){
        if(fabs(z_avg/double(z_avg_count)) <=100 && !z_printed){

          if(p1_store_arr[*it] >0 && p2_store_arr[*it] > 0 && p3_store_arr[*it] > 0 && p4_store_arr[*it] > 0){
            printf("%-5d   %10.2f   %5.1f     %5.1f     %5.1f     %5.1f   ||\n", evt, phi, p1_store_arr[*it],p2_store_arr[*it],
              p3_store_arr[*it], p4_store_arr[*it]);
            printf("N(i)               %5d     %5d     %5d     %5d     ||\n", N1_store_arr[*it],N2_store_arr[*it],N3_store_arr[*it],N4_store_arr[*it]);
            printf("Z(i)               %7.2f   %7.2f    %7.2f   %7.2f  || %7.2f      %5.3f            %5.2f\n", z1_store_arr[*it],z2_store_arr[*it],z3_store_arr[*it],z4_store_arr[*it],
             z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            cout << "===========================================================||======================================== " << endl;
          }
          else if(p1_store_arr[*it] <=0 && p2_store_arr[*it] > 0 && p3_store_arr[*it] > 0 && p4_store_arr[*it] > 0){
            printf("%-5d   %10.2f %5c       %5.1f     %5.1f     %5.1f   ||\n", evt, phi, '-',p2_store_arr[*it],
              p3_store_arr[*it], p4_store_arr[*it]);
            printf("N(i)               %5c     %5d     %5d     %5d     ||\n", '-',N2_store_arr[*it],N3_store_arr[*it],N4_store_arr[*it]);
            printf("Z(i)             %7c     %7.2f    %7.2f   %7.2f  || %7.2f      %5.1f            %5.2f\n", '-',z2_store_arr[*it],z3_store_arr[*it],z4_store_arr[*it],
             z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            cout << "===========================================================||======================================== " << endl; 
          }
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] <= 0 && p3_store_arr[*it] > 0 && p4_store_arr[*it] > 0){
            printf("%-5d   %10.2f   %5.1f   %5c       %5.1f     %5.1f   ||\n", evt, phi, p1_store_arr[*it],'-',
              p3_store_arr[*it], p4_store_arr[*it]);
            printf("N(i)               %5d     %5c     %5d     %5d     ||\n", N1_store_arr[*it],'-',N3_store_arr[*it],N4_store_arr[*it]);
            printf("Z(i)               %7.2f %7c      %7.2f   %7.2f  || %7.2f      %5.1f            %5.2f\n", z1_store_arr[*it],'-',z3_store_arr[*it],z4_store_arr[*it],
             z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            cout << "===========================================================||======================================== " << endl;
          }
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] > 0 && p3_store_arr[*it] <= 0 && p4_store_arr[*it] > 0){
            printf("%-5d   %10.2f   %5.1f     %5.1f   %5c       %5.1f   ||\n", evt, phi, p1_store_arr[*it],p2_store_arr[*it],
              '-', p4_store_arr[*it]);
            printf("N(i)               %5d     %5d     %5c     %5d     ||\n", N1_store_arr[*it],N2_store_arr[*it],'-',N4_store_arr[*it]);
            printf("Z(i)               %7.2f   %7.2f %7c      %7.2f  || %7.2f      %5.1f            %5.2f\n", z1_store_arr[*it],z2_store_arr[*it],'-',z4_store_arr[*it],
             z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            cout << "===========================================================||======================================== " << endl;
          }
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] > 0 && p3_store_arr[*it] > 0 && p4_store_arr[*it] <= 0){
            printf("%-5d   %10.2f   %5.1f     %5.1f     %5.1f   %5c     ||\n", evt, phi, p1_store_arr[*it],p2_store_arr[*it],
              p3_store_arr[*it], '-');
            printf("N(i)               %5d     %5d     %5d     %5c     ||\n", N1_store_arr[*it],N2_store_arr[*it],N3_store_arr[*it],'-');
            printf("Z(i)               %7.2f   %7.2f    %7.2f%7c     || %7.2f      %5.1f            %5.2f\n", z1_store_arr[*it],z2_store_arr[*it],z3_store_arr[*it],'-',
             z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            cout << "===========================================================||======================================== " << endl;
          }
          else if(p1_store_arr[*it] <= 0 && p2_store_arr[*it] >0 && p3_store_arr[*it] <= 0 && p4_store_arr[*it] >0){
            printf("%-5d   %10.2f %5c       %5.1f   %5c       %5.1f   ||\n", evt, phi, '-',p2_store_arr[*it],
              '-', p4_store_arr[*it]);
            printf("N(i)               %5c     %5d     %5c     %5d     ||\n",'-',N2_store_arr[*it],'-',N4_store_arr[*it]);
            printf("Z(i)             %7c     %7.2f %7c      %7.2f  || %7.2f      %5.1f            %5.2f\n", '-',z2_store_arr[*it],'-',z4_store_arr[*it],
             z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            cout << "===========================================================||======================================== " << endl;
          }
          else if(p1_store_arr[*it] <= 0 && p2_store_arr[*it] >0 && p3_store_arr[*it] >0 && p4_store_arr[*it] <= 0){
            printf("%-5d   %10.2f %5c       %5.1f     %5.1f   %5c     ||\n", evt, phi, '-',p2_store_arr[*it],
              p3_store_arr[*it], '-');
            printf("N(i)               %5c     %5d     %5d     %5c     ||\n", '-',N2_store_arr[*it],N3_store_arr[*it],'-');
            printf("Z(i)             %7c     %7.2f    %7.2f%7c     || %7.2f      %5.1f            %5.2f\n", '-',z2_store_arr[*it],z3_store_arr[*it],'-',
             z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            cout << "===========================================================||======================================== " << endl;
          }
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] <= 0 && p3_store_arr[*it] <= 0 && p4_store_arr[*it] >0){
            printf("%-5d   %10.2f   %5.1f   %5c     %5c       %5.1f   ||\n", evt, phi, p1_store_arr[*it],'-',
              '-', p4_store_arr[*it]);
            printf("N(i)               %5d     %5c     %5c     %5d     ||\n", N1_store_arr[*it],'-','-',N4_store_arr[*it]);
            printf("Z(i)               %7.2f %7c   %7c      %7.2f  || %7.2f      %5.1f            %5.2f\n", z1_store_arr[*it],'-','-',z4_store_arr[*it],
             z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            cout << "===========================================================||======================================== " << endl;
          }
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] <= 0 && p3_store_arr[*it] >0 && p4_store_arr[*it] <= 0){
            printf("%-5d   %10.2f   %5.1f   %5c       %5.1f   %5c     ||\n", evt, phi, p1_store_arr[*it],'-',
              p3_store_arr[*it], '-');
            printf("N(i)               %5d     %5c     %5d     %5c     ||\n", N1_store_arr[*it],'-',N3_store_arr[*it],'-');
            printf("Z(i)               %7.2f %7c      %7.2f%7c     || %7.2f      %5.1f            %5.2f\n", z1_store_arr[*it],'-',z3_store_arr[*it],'-',
             z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            cout << "===========================================================||======================================== " << endl;
            
          }
        }   
      }
   
    z_avgs_store.push_back(z_avg/double(z_avg_count));
    }
  }


  return sft_z_selected;


}