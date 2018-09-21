// Functions to use with SFT hit Z-values

//void fillVector(vector<int>& newMyVector){
//
//  cout << "Type in a list of numbers (-1 to stop): ";
//  int input;
//  cin >> input;
//
//  while(input != -1){
//    newMyVector.push_back(input);
//    cin >> input;
//  }
//  cout << endl;
// }

void printVector(const vector<int>& newMyVector){

  cout << "Vector: " ;

  for(unsigned int i=0; i<newMyVector.size(); i++){
    cout << newMyVector[i] << " ";
  }

  cout << endl;
 }

void Z_all(double phi=90., double p1=8., double p2=13., double p3=3., double p4=6.){
  double Z1[16]={0.};
  double Z2[16]={0.};
  double Z3[18]={0.};
  double Z4[18]={0.};

  //double DZ_12[16]={0.};
  //double DZ_34[18]={0.};
  //double DZ[18]={0.};

  double PI = 3.14159265359;
  phi *= PI/180.0;

   // Z_0 CONSTANT OFFSET FOR EACH LAYER
  /////////////////////////////////////////////////////////////
  const double z_offset_const[4] = {-125.298, -130.764, -138.13, -134.949};

  // Number of fibers in each ribbon/layer
  const int fibers_per_ribbon[4] = {15, 15, 17, 17};

  // Angle of wrapping
  const double alpha[4] = {3.57*(PI/180.0), 3.5*(PI/180.0), 3.87*(PI/180.0), 3.79*(PI/180.0)};
                  // OLD  3.6                 3.6               3.9                3.9
  // Radius of layer 
  const double layer_radius[4] = {80.5/2, 82.2/2, 84.2/2, 85.9/2};
  
   // Fiber widths and spacing between each fiber
  //const double fiber_diameter = 1.0;
  //const double fiber_spacing = 0.05;

  // horizontal width of each ribbon
  const double ribbon_width[4] = {15.75/cos(alpha[0]), 15.75/cos(alpha[1]), 
                                  17.85/cos(alpha[2]), 17.85/cos(alpha[3])};


  //if(layer == 1){
    for(int N=1; N<16; N++){
        Z1[N] = z_offset_const[0] + layer_radius[0]*phi*tan(alpha[0]) + (N-1)*ribbon_width[0] 
    //            + (p1-1.0)*ribbon_width[0]/fibers_per_ribbon[0] + fiber_diameter/cos(alpha[1]);  
                + (p1-1.0)*ribbon_width[0]/fibers_per_ribbon[0] + 0.5*1.05/cos(alpha[0]);  

        Z2[N] = z_offset_const[1] + layer_radius[1]*phi*tan(alpha[1]) + (N-1)*ribbon_width[1]
    //            + (p2-1.0)*ribbon_width[1]/fibers_per_ribbon[1] + fiber_diameter/cos(alpha[1]);  
                + (p2-1.0)*ribbon_width[1]/fibers_per_ribbon[1] + 0.5*1.05/cos(alpha[1]);  

        //if(Z1[N]-Z2[N]<0) DZ_12[N] = Z2[N]-Z1[N];
        //else DZ_12[N] = Z1[N]-Z2[N];  
      }

    //cout << "" << endl;
    //for(int N=1; N<16; N++){
    //  cout << "Z1(n=" << N << ") = " << Z1[N] << endl;
    //}
  //} 

  //  if(layer == 2){
  //    for(int N=1; N<16; N++){
  //        Z2[N] = z_offset_const[1] + layer_radius[1]*phi*tan(alpha[1]) + (N-1)*ribbon_width[1]
  //                + (p2-1.0)*ribbon_width[1]/fibers_per_ribbon[1] + fiber_diameter/cos(alpha[1]);     
  //    } 

  //    cout << "" << endl;
  //    for(int N=1; N<16; N++){
  //      cout << "Z2(n=" << N << ") = " << Z2[N] << endl;
  //    }
  //  } 

  //  if(layer == 3){
    for(int N=1; N<18; N++){
        Z3[N] = z_offset_const[2] + layer_radius[2]*(2*PI-phi)*tan(alpha[2]) + (N-1)*ribbon_width[2]
  //                + (fibers_per_ribbon[2]-p3)*ribbon_width[2]/fibers_per_ribbon[2] + fiber_diameter/cos(alpha[2]); 
                + (fibers_per_ribbon[2]-p3)*ribbon_width[2]/fibers_per_ribbon[2] + 0.5*1.05/cos(alpha[2]); 
 
        Z4[N] = z_offset_const[3] + layer_radius[3]*(2*PI-phi)*tan(alpha[3]) + (N-1)*ribbon_width[3]
  //                + (fibers_per_ribbon[3]-p4)*ribbon_width[3]/fibers_per_ribbon[3] + fiber_diameter/cos(alpha[3]);
                + (fibers_per_ribbon[3]-p4)*ribbon_width[3]/fibers_per_ribbon[3] + 0.5*1.05/cos(alpha[3]);

        //if(Z3[N]-Z4[N]<0) DZ_34[N] = Z4[N]-Z3[N];
        //else DZ_34[N] = Z3[N]-Z4[N];  
    }

  //    cout << "" << endl;
  //    for(int N=1; N<18; N++){
  //      cout << "Z3(n=" << N << ") = " << Z3[N] << endl;
  //    }
  //  } 

  //  if(layer == 4){
  //    for(int N=1; N<18; N++){
  //    Z4[N] = z_offset_const[3] + layer_radius[3]*(2*PI-phi)*tan(alpha[3]) + (N-1)*ribbon_width[3]
  //            + (fibers_per_ribbon[3]-p4)*ribbon_width[3]/fibers_per_ribbon[3] + fiber_diameter/cos(alpha[3]);
  //    }

  //    cout << "" << endl;
  //    for(int N=1; N<18; N++){
  //      cout << "Z4(n=" << N << ") = " << Z4[N] << endl;
  //    }
  //  } 

  cout << "" << endl;
  cout << " " << endl;
  cout << " =================================================" << endl;
  printf(" p1 = %3.1f,   p2 = %3.1f,   p3 = %3.1f,   p4 = %3.1f                  \n", p1, p2, p3, p4);
  cout << " =================================================" << endl;
  //  cout << "         |                    TOP                     |                   BOTTOM                  "<< endl;
  printf( "         |             phi = %3.3f\n",phi*180./PI);
  cout << "    N    |========================================" << endl;
  cout << "         |      Z1        Z2        Z3        Z4      " << endl;
  cout << " ========|========================================" << endl;
  for(int i=1; i<16; i++){
  printf(" n = %2d  |  %8.3f  %8.3f  %8.3f  %8.3f\n", i,Z1[i],Z2[i],Z3[i],Z4[i]);


  //  cout << "" << endl;
  //  cout << "    N            Z1        Z2        Z3        Z4" << endl;
  //  for(int i=1; i<16; i++){
  //  printf(" n = %2d  *  %8.3f  %8.3f  %8.3f  %8.3f\n", i,Z1[i],Z2[i],Z3[i],Z4[i]);
  //cout << " n = " << i << "  ||   "  << Z1[i] << "  " << Z2[i] << "  " << Z3[i] << "  " << Z4[i] << endl;
  //cout << "" << endl;
  }

  for(int j=16; j<18; j++){
  printf(" n = %2d  *  --------  --------  %8.3f  %8.3f\n", j,Z3[j],Z4[j]);
  //    cout << " n = " << j << "  ||   " << "---" << "  " << "---" << "  " << Z3[j] << "  " << Z4[j] << endl;
  }
  cout << "" << endl;

  }

void Z_all_CR(double phi_top=90., double phi_bot=270., double p1=8., double p2=13., double p3=3., double p4=6.){
 
  double Z1_top[16]={0.};  double Z1_bot[16]={0.};
  double Z2_top[16]={0.};  double Z2_bot[16]={0.};
  double Z3_top[18]={0.};  double Z3_bot[18]={0.};
  double Z4_top[18]={0.};  double Z4_bot[18]={0.};

  //double DZ_12[16]={0.};
  //double DZ_34[18]={0.};
  //double DZ[18]={0.};

  double PI = 3.14159265359;
  phi_top *= PI/180.0;
  phi_bot *= PI/180.0;

   // Z_0 CONSTANT OFFSET FOR EACH LAYER
  /////////////////////////////////////////////////////////////
  const double z_offset_const[4] = {-125.298, -130.764, -138.13, -134.949};

  // Number of fibers in each ribbon/layer
  const int fibers_per_ribbon[4] = {15, 15, 17, 17};

  // Angle of wrapping
  const double alpha[4] = {3.57*(PI/180.0), 3.5*(PI/180.0), 3.87*(PI/180.0), 3.79*(PI/180.0)};
                  // OLD  3.6                 3.6               3.9                3.9
  // Radius of layer 
  const double layer_radius[4] = {80.5/2, 82.2/2, 84.2/2, 85.9/2};
  
   // Fiber widths and spacing between each fiber
  //const double fiber_diameter = 1.0;
  //const double fiber_spacing = 0.05;

  // horizontal width of each ribbon
  const double ribbon_width[4] = {15.75/cos(alpha[0]), 15.75/cos(alpha[1]), 
                                  17.85/cos(alpha[2]), 17.85/cos(alpha[3])};


  //  if(layer == 1){
    for(int N=1; N<16; N++){
        Z1_top[N] = z_offset_const[0] + layer_radius[0]*phi_top*tan(alpha[0]) + (N-1)*ribbon_width[0] 
                + (p1-1.0)*ribbon_width[0]/fibers_per_ribbon[0] + 0.5*1.05/cos(alpha[0]);  

        Z1_bot[N] = z_offset_const[0] + layer_radius[0]*phi_bot*tan(alpha[0]) + (N-1)*ribbon_width[0] 
                + (p1-1.0)*ribbon_width[0]/fibers_per_ribbon[0] + 0.5*1.05/cos(alpha[0]);  

        Z2_top[N] = z_offset_const[1] + layer_radius[1]*phi_top*tan(alpha[1]) + (N-1)*ribbon_width[1]
                + (p2-1.0)*ribbon_width[1]/fibers_per_ribbon[1] + 0.5*1.05/cos(alpha[1]);  

        Z2_bot[N] = z_offset_const[1] + layer_radius[1]*phi_bot*tan(alpha[1]) + (N-1)*ribbon_width[1]
                + (p2-1.0)*ribbon_width[1]/fibers_per_ribbon[1] + 0.5*1.05/cos(alpha[1]);  

    //    if(Z1[N]-Z2[N]<0) DZ_12[N] = Z2[N]-Z1[N];
    //    else DZ_12[N] = Z1[N]-Z2[N];  
    }

  //    cout << "" << endl;
  //    for(int N=1; N<16; N++){
  //      cout << "Z1(n=" << N << ") = " << Z1[N] << endl;
  //    }
  //  } 

  //  if(layer == 2){
  //    for(int N=1; N<16; N++){
  //        Z2[N] = z_offset_const[1] + layer_radius[1]*phi*tan(alpha[1]) + (N-1)*ribbon_width[1]
  //                + (p2-1.0)*ribbon_width[1]/fibers_per_ribbon[1] + fiber_diameter/cos(alpha[1]);     
  //    } 

  //    cout << "" << endl;
  //    for(int N=1; N<16; N++){
  //      cout << "Z2(n=" << N << ") = " << Z2[N] << endl;
  //    }
  //  } 

  //  if(layer == 3){
    for(int N=1; N<18; N++){
        Z3_top[N] = z_offset_const[2] + layer_radius[2]*(2*PI-phi_top)*tan(alpha[2]) + (N-1)*ribbon_width[2]
                + (fibers_per_ribbon[2]-p3)*ribbon_width[2]/fibers_per_ribbon[2] + 0.5*1.05/cos(alpha[2]); 

        Z3_bot[N] = z_offset_const[2] + layer_radius[2]*(2*PI-phi_bot)*tan(alpha[2]) + (N-1)*ribbon_width[2]
                + (fibers_per_ribbon[2]-p3)*ribbon_width[2]/fibers_per_ribbon[2] + 0.5*1.05/cos(alpha[2]); 
 
        Z4_top[N] = z_offset_const[3] + layer_radius[3]*(2*PI-phi_top)*tan(alpha[3]) + (N-1)*ribbon_width[3]
                + (fibers_per_ribbon[3]-p4)*ribbon_width[3]/fibers_per_ribbon[3] + 0.5*1.05/cos(alpha[3]);

        Z4_bot[N] = z_offset_const[3] + layer_radius[3]*(2*PI-phi_bot)*tan(alpha[3]) + (N-1)*ribbon_width[3]
                + (fibers_per_ribbon[3]-p4)*ribbon_width[3]/fibers_per_ribbon[3] + 0.5*1.05/cos(alpha[3]);

   //     if(Z3[N]-Z4[N]<0) DZ_34[N] = Z4[N]-Z3[N];
   //     else DZ_34[N] = Z3[N]-Z4[N];  
    }

  //    cout << "" << endl;
  //    for(int N=1; N<18; N++){
  //      cout << "Z3(n=" << N << ") = " << Z3[N] << endl;
  //    }
  //  } 

  //  if(layer == 4){
  //    for(int N=1; N<18; N++){
  //    Z4[N] = z_offset_const[3] + layer_radius[3]*(2*PI-phi)*tan(alpha[3]) + (N-1)*ribbon_width[3]
  //            + (fibers_per_ribbon[3]-p4)*ribbon_width[3]/fibers_per_ribbon[3] + fiber_diameter/cos(alpha[3]);
  //    }

  //    cout << "" << endl;
  //    for(int N=1; N<18; N++){
  //      cout << "Z4(n=" << N << ") = " << Z4[N] << endl;
  //    }
  //  } 

  cout << "" << endl;
  cout << " " << endl;
  cout << " ======================================================" << endl;
  printf(" p1 = %3.1f,   p2 = %3.1f,   p3 = %3.1f,   p4 = %3.1f            |      \n", p1, p2, p3, p4);
  cout << " =====================================================|===========================================" << endl;
  //  cout << "         |                    TOP                     |                   BOTTOM                  "<< endl;
  printf( "         |             TOP (phi = %3.3f)             |            BOTTOM (phi = %3.3f)\n",phi_top, phi_bot);
  cout << "    N    |============================================|===========================================" << endl;
  cout << "         |      Z1        Z2        Z3        Z4      |        Z1        Z2        Z3        Z4   " << endl;
  cout << " ========|============================================|===========================================" << endl;
  for(int i=1; i<16; i++){
  printf(" n = %2d  |  %8.3f  %8.3f  %8.3f  %8.3f    |    %8.3f  %8.3f  %8.3f  %8.3f\n", i,Z1_top[i],Z2_top[i],Z3_top[i],Z4_top[i],Z1_bot[i],Z2_bot[i],Z3_bot[i],Z4_bot[i]);
  //cout << " n = " << i << "  ||   "  << Z1[i] << "  " << Z2[i] << "  " << Z3[i] << "  " << Z4[i] << endl;
  //cout << "" << endl;
  }

  for(int j=16; j<18; j++){
  printf(" n = %2d  |  --------  --------  %8.3f  %8.3f    |    --------  --------  %8.3f  %8.3f\n", j,Z3_top[j],Z4_top[j], Z3_bot[j],Z4_bot[j]);
  //    cout << " n = " << j << "  ||   " << "---" << "  " << "---" << "  " << Z3[j] << "  " << Z4[j] << endl;
  }
  cout << "" << endl;


  }

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
               + 0.5*(fiber_diameter+fiber_spacing)/cos(alpha[0]);   
  }
  else if(layer == 2){
    Z_position = z_offset_const[1] + layer_radius[1]*phi*tan(alpha[1]) 
               + (N-1)*ribbon_width[1]
               + (p-1.0)*ribbon_width[1]/fibers_per_ribbon[1]
               + 0.5*(fiber_diameter+fiber_spacing)/cos(alpha[1]);     
  }
  else if(layer == 3){
    Z_position = z_offset_const[2] + layer_radius[2]*(2*PI-phi)*tan(alpha[2]) 
               + (N-1)*ribbon_width[2]
               + (fibers_per_ribbon[2]-p)*ribbon_width[2]/fibers_per_ribbon[2]
               + 0.5*(fiber_diameter+fiber_spacing)/cos(alpha[2]); 
  }
  else if(layer == 4){
    Z_position = z_offset_const[3] + layer_radius[3]*(2*PI-phi)*tan(alpha[3]) 
               + (N-1)*ribbon_width[3]
               + (fibers_per_ribbon[3]-p)*ribbon_width[3]/fibers_per_ribbon[3]
               + 0.5*(fiber_diameter+fiber_spacing)/cos(alpha[3]);
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
  //const int fibers_per_ribbon[4] = {15, 15, 17, 17};  

  // selected p's in an array
  //const double sft_p_fiber[4] = {p1, p2, p3, p4};


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

  //bool same_turn_flag = false;

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
                                    + Z_positions[3][layer_3_low_index])/2.0;  
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
  //double z_total_thr = 5.0;

  if(fabs(return_array[10] - return_array[11]) > z12_thr){
    return_array[16] = 0;
  }
  if(fabs(return_array[12] - return_array[13]) > z34_thr){
    return_array[17] = 0;
  }


  // This removes the flags on the p2-p1 and p4-p3 conditions.
   // return_array[14] = 2;
   // return_array[15] = 2;

  return return_array;}

vector<double> Layer_averaging(vector<int> p_fibers, int layer){
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


  //int p_cluster_separation = 0;

  // Find clustering of p fibers
  ////////////////////////////
  //bool first_cluster = true;
  //int cluster_spacing = 0;
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
  //int p_cluster_length_arr[max_length];

  double p_centroids[max_length];


  for(int i = 0; i< max_length; i++){
    p_cluster_index_arr[i] = -10000;

    //p_cluster_length_arr[i] = -1;

    p_centroids[i] = -10000.0;  
  } 


  // Store all index/lengths in arrays.
  for(unsigned int i = 0; i < p_cluster_index.size(); i++){
    p_cluster_index_arr[i] = *(p_cluster_index.begin() + i);
  }

  for(unsigned int i = 0; i < p_cluster_length.size(); i++){
    //p_cluster_length_arr[i] = *(p_cluster_length.begin() + i);
  }

  // Calculate all centroids for C2-all
  
  //int centroid_sum = 0;
  //int centroid_den = 0;

  
  //cout << "  Layer " << layer << ": ";
  for(unsigned int i = 0; i < 100; i++){
    if(p_cluster_index_arr[i] < -1000)
      break;

    p_centroids[i] = p_cluster_index[i] +(p_cluster_length[i] - 1)/2.0 + 1;
    //cout << p_centroids[i] << "  ";
    p_fibers_d.push_back(p_centroids[i]);
  }  

  //cout << endl;
  





  // for(vector<int>::iterator it = p_fibers.begin(); it != p_fibers.end(); it++){
  //   p_fibers_d.push_back(double(*it));
  // }

  // uncomment for flagging bad sft 
  // *bad_sft_flag = true;

  return p_fibers_d;
 }

//vector<double> Batch_Layer_averaging(vector<int> p_fibers, bool *bad_sft_flag, int layer){
vector<double> Batch_Layer_averaging(vector<int> p_fibers, int layer){
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


  //int p_cluster_separation = 0;

  // Find clustering of p fibers
  ////////////////////////////
  //bool first_cluster = true;
  //int cluster_spacing = 0;
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
  //int p_cluster_length_arr[max_length];

  double p_centroids[max_length];


  for(int i = 0; i< max_length; i++){
    p_cluster_index_arr[i] = -10000;

    //p_cluster_length_arr[i] = -1;

    p_centroids[i] = -10000.0;  
  } 


  // Store all index/lengths in arrays.
  for(unsigned int i = 0; i < p_cluster_index.size(); i++){
    p_cluster_index_arr[i] = *(p_cluster_index.begin() + i);
  }

  for(unsigned int i = 0; i < p_cluster_length.size(); i++){
    //p_cluster_length_arr[i] = *(p_cluster_length.begin() + i);
  }

  // Calculate all centroids for C2-all
  
  //int centroid_sum = 0;
  //int centroid_den = 0;

  //cout << "  Layer " << layer << ": ";
  for(unsigned int i = 0; i < 100; i++){
    if(p_cluster_index_arr[i] < -1000)
      break;

    p_centroids[i] = p_cluster_index[i] +(p_cluster_length[i] - 1)/2.0 + 1;
    //cout << p_centroids[i] << "  ";
    p_fibers_d.push_back(p_centroids[i]);
  }  







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

  //float R_TARGET = 29.0;
  //float R_TOF1 = 47.1;
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

  /*  cout << endl;
  cout << endl;
  cout << "####### DEBUG SFT print  ########" << endl;
  cout << "p1 = ";
  for (unsigned i1=0; i1<layer_1_p_fibers.size(); i1++) cout << layer_1_p_fibers[i1] << "  ";
  cout << endl; 
  cout << "p2 = "; 
  for (unsigned i2=0; i2<layer_2_p_fibers.size(); i2++) cout << layer_2_p_fibers[i2] << "  ";
  cout << endl;
  cout << "p3 = ";
  for (unsigned i3=0; i3<layer_3_p_fibers.size(); i3++) cout << layer_3_p_fibers[i3] << "  ";
  cout << endl;
  cout << "p4 = ";
  for (unsigned i4=0; i4<layer_4_p_fibers.size(); i4++) cout << layer_4_p_fibers[i4] << "  ";
  cout << endl;
  cout << "###########################################" << endl;
  cout << endl;
  cout << endl;
  */

  // pass vectors to new p-averaging function. return vector<double> of averaged p-fibers.
  // vector<double> Layer_averaging(vector<int> p_fibers)
  vector<double> layer_1_p_fibers_d;
  vector<double> layer_2_p_fibers_d;
  vector<double> layer_3_p_fibers_d;
  vector<double> layer_4_p_fibers_d;
  
  //bool bad_sft = false;
  //bool *bad_sft_flag = &bad_sft;

  layer_1_p_fibers_d = Layer_averaging(layer_1_p_fibers, 1);
  layer_2_p_fibers_d = Layer_averaging(layer_2_p_fibers, 2);
  layer_3_p_fibers_d = Layer_averaging(layer_3_p_fibers, 3);
  layer_4_p_fibers_d = Layer_averaging(layer_4_p_fibers, 4);
  
  cout << endl;
  cout << "  Layer 1: ";
  for(unsigned l1=0; l1<layer_1_p_fibers_d.size(); l1++) cout << layer_1_p_fibers_d[l1] << "  ";
  cout << endl;
  
  cout << "  Layer 2: ";
  for(unsigned l2=0; l2<layer_2_p_fibers_d.size(); l2++) cout << layer_2_p_fibers_d[l2] << "  ";
  cout << endl;
  
  cout << "  Layer 3: ";
  for(unsigned l3=0; l3<layer_3_p_fibers_d.size(); l3++) cout << layer_3_p_fibers_d[l3] << "  ";
  cout << endl;
  
  cout << "  Layer 4: ";
  for(unsigned l4=0; l4<layer_4_p_fibers_d.size(); l4++) cout << layer_4_p_fibers_d[l4] << "  ";
  cout << endl;

 
  // Test if more than 1 layer is empty or if any layer has more than 3 p fibers
  if((layer_1_p_fibers_d.size() == 0 && layer_2_p_fibers_d.size() == 0) ||  
   (layer_3_p_fibers_d.size() == 0 && layer_4_p_fibers_d.size() == 0)){
    //if(to_print) cout << "Need 1 hit in each helicity." << endl;
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

  //double z_avg_min = 0;
  //double delta_z_min = 1000;
  //int delta_z_min_index = 0;

  vector<int> good_delta_z_index;

  vector<double> z_avgs_store;
  bool z_printed = false;


  double delta_z_thr = 10.0;  // Delta_Z

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
          z_avgs_store.push_back(z_avg/double(z_avg_count));    
        }   
      }
    }
  }

  return sft_z_selected;
 }


//vector<double> Z_Avg(double ADC_High_SFT_corr[128], int has_TDC_SFT_hit[128], int SFT_channel_to_fiber[128],
// int evt, double phi, bool to_print, double TOF1_pos, bool use_TOF1, double C2X_centroid, double len_in_target){
vector<double> Z_Avg(double ADC_High_SFT_corr[128], int has_TDC_SFT_hit[128], int SFT_channel_to_fiber[128],
 double phi, bool to_print, double TOF1_pos, bool use_TOF1, double C2X_centroid, double len_in_target){

  vector<double> ZZ;
  ZZ.clear();

  vector<double> DZ;
  DZ.clear();

  vector<double> Track_Length;
  Track_Length.clear();

  const int sft_hit_count_max = 3;
  const int NUM_ITERATIONS = 100;

  //float R_TARGET = 29.0;
  //float R_TOF1 = 47.1;
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

  /*  cout << endl;
  cout << endl;
  cout << "####### DEBUG SFT print  ########" << endl;
  cout << "p1 = ";
  for (unsigned i1=0; i1<layer_1_p_fibers.size(); i1++) cout << layer_1_p_fibers[i1] << "  ";
  cout << endl; 
  cout << "p2 = "; 
  for (unsigned i2=0; i2<layer_2_p_fibers.size(); i2++) cout << layer_2_p_fibers[i2] << "  ";
  cout << endl;
  cout << "p3 = ";
  for (unsigned i3=0; i3<layer_3_p_fibers.size(); i3++) cout << layer_3_p_fibers[i3] << "  ";
  cout << endl;
  cout << "p4 = ";
  for (unsigned i4=0; i4<layer_4_p_fibers.size(); i4++) cout << layer_4_p_fibers[i4] << "  ";
  cout << endl;
  cout << "###########################################" << endl;
  cout << endl;
  cout << endl;
  */

  // pass vectors to new p-averaging function. return vector<double> of averaged p-fibers.
  // vector<double> Layer_averaging(vector<int> p_fibers)
  vector<double> layer_1_p_fibers_d;
  vector<double> layer_2_p_fibers_d;
  vector<double> layer_3_p_fibers_d;
  vector<double> layer_4_p_fibers_d;
  
  //bool bad_sft = false;
  //bool *bad_sft_flag = &bad_sft;

  layer_1_p_fibers_d = Layer_averaging(layer_1_p_fibers, 1);
  layer_2_p_fibers_d = Layer_averaging(layer_2_p_fibers, 2);
  layer_3_p_fibers_d = Layer_averaging(layer_3_p_fibers, 3);
  layer_4_p_fibers_d = Layer_averaging(layer_4_p_fibers, 4);
  
  //cout << endl;
  //cout << "  Layer 1: ";
  //for(unsigned l1=0; l1<layer_1_p_fibers_d.size(); l1++) cout << layer_1_p_fibers_d[l1] << "  ";
  //cout << endl;
  
  //cout << "  Layer 2: ";
  //for(unsigned l2=0; l2<layer_2_p_fibers_d.size(); l2++) cout << layer_2_p_fibers_d[l2] << "  ";
  //cout << endl;
  
  //cout << "  Layer 3: ";
  //for(unsigned l3=0; l3<layer_3_p_fibers_d.size(); l3++) cout << layer_3_p_fibers_d[l3] << "  ";
  //cout << endl;
  
  //cout << "  Layer 4: ";
  //for(unsigned l4=0; l4<layer_4_p_fibers_d.size(); l4++) cout << layer_4_p_fibers_d[l4] << "  ";
  //cout << endl;

 
  // Test if more than 1 layer is empty or if any layer has more than 3 p fibers
  if((layer_1_p_fibers_d.size() == 0 && layer_2_p_fibers_d.size() == 0) ||  
   (layer_3_p_fibers_d.size() == 0 && layer_4_p_fibers_d.size() == 0)){
    // if(to_print) cout << "Need 1 hit in each helicity." << endl;
    return ZZ;
    //return -10000;
  }

  bool flag = false;
  if(layer_1_p_fibers_d.size() > sft_hit_count_max){
    if(to_print)
      //cout << "More than " << sft_hit_count_max << " SFT hits in layer 1" << endl;
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
    return ZZ;
    //return -10000; 

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

  //double z_avg_min = 0;
  //double delta_z_min = 1000;
  //int delta_z_min_index = 0;

  vector<int> good_delta_z_index;

  bool z_printed = false;


  double delta_z_thr = 10.0;  // Delta_Z

  double sft_delta_z_selected = 1000;
  //double sft_z_selected = 0;

  for(int i = 0; i < iteration_counter; i++){
    if(delta_z_arr[i] < delta_z_thr){
      good_delta_z_index.push_back(i);
    }
  }

  //if(to_print)
    //cout<<endl;

  if(to_print){
    //cout << "Event #     Angle      L1        L2        L3        L4    ||   Z Avg      Delta Z      Track Length" << endl;
    //cout << "===========================================================||========================================" << endl;
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
          //sft_z_selected = z_avg/double(z_avg_count);
          sft_delta_z_selected = delta_z_arr[*it];
        }
      }
      else{
        //sft_z_selected = z_avg/double(z_avg_count);
        sft_delta_z_selected = delta_z_arr[*it];
      }
    }

    if(p1_p2_flag[*it] > 1 && p3_p4_flag[*it] > 1){
      if(ZZ.empty()){
        z_printed = false;
      }
      else{
        for(vector<double>::iterator db = ZZ.begin(); db != ZZ.end(); db++){
          if(fabs(*db - z_avg/double(z_avg_count)) < 0.01){
            z_printed = true;
          }
        }
      }

      if(to_print){
        if(fabs(z_avg/double(z_avg_count)) <=100 && !z_printed){

          if(p1_store_arr[*it] >0 && p2_store_arr[*it] > 0 && p3_store_arr[*it] > 0 && p4_store_arr[*it] > 0){
            //printf("%-5d   %10.2f   %5.1f     %5.1f     %5.1f     %5.1f   ||\n", evt, phi, p1_store_arr[*it],p2_store_arr[*it],
            //  p3_store_arr[*it], p4_store_arr[*it]);
            //printf("N(i)               %5d     %5d     %5d     %5d     ||\n", N1_store_arr[*it],N2_store_arr[*it],N3_store_arr[*it],N4_store_arr[*it]);
            //printf("Z(i)               %7.2f   %7.2f    %7.2f   %7.2f  || %7.2f      %5.3f            %5.2f\n", z1_store_arr[*it],z2_store_arr[*it],z3_store_arr[*it],z4_store_arr[*it],
             //z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            //cout << "===========================================================||======================================== " << endl;
            DZ.push_back(delta_z_arr[*it]);
            Track_Length.push_back(total_length);
          }
          else if(p1_store_arr[*it] <=0 && p2_store_arr[*it] > 0 && p3_store_arr[*it] > 0 && p4_store_arr[*it] > 0){
            //printf("%-5d   %10.2f %5c       %5.1f     %5.1f     %5.1f   ||\n", evt, phi, '-',p2_store_arr[*it],
            //  p3_store_arr[*it], p4_store_arr[*it]);
            //printf("N(i)               %5c     %5d     %5d     %5d     ||\n", '-',N2_store_arr[*it],N3_store_arr[*it],N4_store_arr[*it]);
            //printf("Z(i)             %7c     %7.2f    %7.2f   %7.2f  || %7.2f      %5.1f            %5.2f\n", '-',z2_store_arr[*it],z3_store_arr[*it],z4_store_arr[*it],
            // z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            //cout << "===========================================================||======================================== " << endl; 
            DZ.push_back(delta_z_arr[*it]);
            Track_Length.push_back(total_length);
          }
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] <= 0 && p3_store_arr[*it] > 0 && p4_store_arr[*it] > 0){
            //printf("%-5d   %10.2f   %5.1f   %5c       %5.1f     %5.1f   ||\n", evt, phi, p1_store_arr[*it],'-',
            //  p3_store_arr[*it], p4_store_arr[*it]);
            //printf("N(i)               %5d     %5c     %5d     %5d     ||\n", N1_store_arr[*it],'-',N3_store_arr[*it],N4_store_arr[*it]);
            //printf("Z(i)               %7.2f %7c      %7.2f   %7.2f  || %7.2f      %5.1f            %5.2f\n", z1_store_arr[*it],'-',z3_store_arr[*it],z4_store_arr[*it],
            // z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            //cout << "===========================================================||======================================== " << endl;
            DZ.push_back(delta_z_arr[*it]);
            Track_Length.push_back(total_length);
          }
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] > 0 && p3_store_arr[*it] <= 0 && p4_store_arr[*it] > 0){
            //printf("%-5d   %10.2f   %5.1f     %5.1f   %5c       %5.1f   ||\n", evt, phi, p1_store_arr[*it],p2_store_arr[*it],
            //  '-', p4_store_arr[*it]);
            //printf("N(i)               %5d     %5d     %5c     %5d     ||\n", N1_store_arr[*it],N2_store_arr[*it],'-',N4_store_arr[*it]);
            //printf("Z(i)               %7.2f   %7.2f %7c      %7.2f  || %7.2f      %5.1f            %5.2f\n", z1_store_arr[*it],z2_store_arr[*it],'-',z4_store_arr[*it],
            // z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            //cout << "===========================================================||======================================== " << endl;
            DZ.push_back(delta_z_arr[*it]);
            Track_Length.push_back(total_length);
          }
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] > 0 && p3_store_arr[*it] > 0 && p4_store_arr[*it] <= 0){
            //printf("%-5d   %10.2f   %5.1f     %5.1f     %5.1f   %5c     ||\n", evt, phi, p1_store_arr[*it],p2_store_arr[*it],
            //  p3_store_arr[*it], '-');
            //printf("N(i)               %5d     %5d     %5d     %5c     ||\n", N1_store_arr[*it],N2_store_arr[*it],N3_store_arr[*it],'-');
            //printf("Z(i)               %7.2f   %7.2f    %7.2f%7c     || %7.2f      %5.1f            %5.2f\n", z1_store_arr[*it],z2_store_arr[*it],z3_store_arr[*it],'-',
            // z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            //cout << "===========================================================||======================================== " << endl;
            DZ.push_back(delta_z_arr[*it]);
            Track_Length.push_back(total_length);
          }
          else if(p1_store_arr[*it] <= 0 && p2_store_arr[*it] >0 && p3_store_arr[*it] <= 0 && p4_store_arr[*it] >0){
            //printf("%-5d   %10.2f %5c       %5.1f   %5c       %5.1f   ||\n", evt, phi, '-',p2_store_arr[*it],
            //  '-', p4_store_arr[*it]);
            //printf("N(i)               %5c     %5d     %5c     %5d     ||\n",'-',N2_store_arr[*it],'-',N4_store_arr[*it]);
            //printf("Z(i)             %7c     %7.2f %7c      %7.2f  || %7.2f      %5.1f            %5.2f\n", '-',z2_store_arr[*it],'-',z4_store_arr[*it],
            // z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            //cout << "===========================================================||======================================== " << endl;
            DZ.push_back(delta_z_arr[*it]);
            Track_Length.push_back(total_length);
          }
          else if(p1_store_arr[*it] <= 0 && p2_store_arr[*it] >0 && p3_store_arr[*it] >0 && p4_store_arr[*it] <= 0){
            //printf("%-5d   %10.2f %5c       %5.1f     %5.1f   %5c     ||\n", evt, phi, '-',p2_store_arr[*it],
            //  p3_store_arr[*it], '-');
            //printf("N(i)               %5c     %5d     %5d     %5c     ||\n", '-',N2_store_arr[*it],N3_store_arr[*it],'-');
            //printf("Z(i)             %7c     %7.2f    %7.2f%7c     || %7.2f      %5.1f            %5.2f\n", '-',z2_store_arr[*it],z3_store_arr[*it],'-',
            // z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            //cout << "===========================================================||======================================== " << endl;
            DZ.push_back(delta_z_arr[*it]);
            Track_Length.push_back(total_length);
          }
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] <= 0 && p3_store_arr[*it] <= 0 && p4_store_arr[*it] >0){
            //printf("%-5d   %10.2f   %5.1f   %5c     %5c       %5.1f   ||\n", evt, phi, p1_store_arr[*it],'-',
            //  '-', p4_store_arr[*it]);
            //printf("N(i)               %5d     %5c     %5c     %5d     ||\n", N1_store_arr[*it],'-','-',N4_store_arr[*it]);
            //printf("Z(i)               %7.2f %7c   %7c      %7.2f  || %7.2f      %5.1f            %5.2f\n", z1_store_arr[*it],'-','-',z4_store_arr[*it],
             //z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            //cout << "===========================================================||======================================== " << endl;
            DZ.push_back(delta_z_arr[*it]);
            Track_Length.push_back(total_length);
          }
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] <= 0 && p3_store_arr[*it] >0 && p4_store_arr[*it] <= 0){
            //printf("%-5d   %10.2f   %5.1f   %5c       %5.1f   %5c     ||\n", evt, phi, p1_store_arr[*it],'-',
            //  p3_store_arr[*it], '-');
            //printf("N(i)               %5d     %5c     %5d     %5c     ||\n", N1_store_arr[*it],'-',N3_store_arr[*it],'-');
            //printf("Z(i)               %7.2f %7c      %7.2f%7c     || %7.2f      %5.1f            %5.2f\n", z1_store_arr[*it],'-',z3_store_arr[*it],'-',
            // z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            //cout << "===========================================================||======================================== " << endl;
            DZ.push_back(delta_z_arr[*it]);         
            Track_Length.push_back(total_length);
          }
          ZZ.push_back(z_avg/double(z_avg_count));
        }   
      }
    }
  }

  //for(unsigned i=0; i<ZZ.size(); i++){
    //cout << "UUUU :   " << ZZ[i] << "   " << DZ[i] << "    " << Track_Length[i] << endl; 
    //cout << "UUUU :   " << ZZ[i] << "   " << DZ[i] << endl; 
  //}

  //cout << endl;

  return ZZ;
 }


//vector<double> Delta_Z(double ADC_High_SFT_corr[128], int has_TDC_SFT_hit[128], int SFT_channel_to_fiber[128],
// int evt, double phi, bool to_print, double TOF1_pos, bool use_TOF1, double C2X_centroid, double len_in_target){
vector<double> Delta_Z(double ADC_High_SFT_corr[128], int has_TDC_SFT_hit[128], int SFT_channel_to_fiber[128],
 double phi, bool to_print, double TOF1_pos, bool use_TOF1){

  vector<double> ZZ;
  ZZ.clear();

  vector<double> DZ;
  DZ.clear();

  const int sft_hit_count_max = 3;
  const int NUM_ITERATIONS = 100;

  //float R_TARGET = 29.0;
  //float R_TOF1 = 47.1;
  //float R_SFT_L1 = 40.0;
  //float R_C2 = 629.4;

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

  /*  cout << endl;
  cout << endl;
  cout << "####### DEBUG SFT print  ########" << endl;
  cout << "p1 = ";
  for (unsigned i1=0; i1<layer_1_p_fibers.size(); i1++) cout << layer_1_p_fibers[i1] << "  ";
  cout << endl; 
  cout << "p2 = "; 
  for (unsigned i2=0; i2<layer_2_p_fibers.size(); i2++) cout << layer_2_p_fibers[i2] << "  ";
  cout << endl;
  cout << "p3 = ";
  for (unsigned i3=0; i3<layer_3_p_fibers.size(); i3++) cout << layer_3_p_fibers[i3] << "  ";
  cout << endl;
  cout << "p4 = ";
  for (unsigned i4=0; i4<layer_4_p_fibers.size(); i4++) cout << layer_4_p_fibers[i4] << "  ";
  cout << endl;
  cout << "###########################################" << endl;
  cout << endl;
  cout << endl;
  */

  // pass vectors to new p-averaging function. return vector<double> of averaged p-fibers.
  // vector<double> Layer_averaging(vector<int> p_fibers)
  vector<double> layer_1_p_fibers_d;
  vector<double> layer_2_p_fibers_d;
  vector<double> layer_3_p_fibers_d;
  vector<double> layer_4_p_fibers_d;
  
  //bool bad_sft = false;
  //bool *bad_sft_flag = &bad_sft;

  layer_1_p_fibers_d = Layer_averaging(layer_1_p_fibers, 1);
  layer_2_p_fibers_d = Layer_averaging(layer_2_p_fibers, 2);
  layer_3_p_fibers_d = Layer_averaging(layer_3_p_fibers, 3);
  layer_4_p_fibers_d = Layer_averaging(layer_4_p_fibers, 4);
  
  //cout << endl;
  //cout << "  Layer 1: ";
  //for(unsigned l1=0; l1<layer_1_p_fibers_d.size(); l1++) cout << layer_1_p_fibers_d[l1] << "  ";
  //cout << endl;
  
  //cout << "  Layer 2: ";
  //for(unsigned l2=0; l2<layer_2_p_fibers_d.size(); l2++) cout << layer_2_p_fibers_d[l2] << "  ";
  //cout << endl;
  
  //cout << "  Layer 3: ";
  //for(unsigned l3=0; l3<layer_3_p_fibers_d.size(); l3++) cout << layer_3_p_fibers_d[l3] << "  ";
  //cout << endl;
  
  //cout << "  Layer 4: ";
  //for(unsigned l4=0; l4<layer_4_p_fibers_d.size(); l4++) cout << layer_4_p_fibers_d[l4] << "  ";
  //cout << endl;

 
  // Test if more than 1 layer is empty or if any layer has more than 3 p fibers
  if((layer_1_p_fibers_d.size() == 0 && layer_2_p_fibers_d.size() == 0) ||  
   (layer_3_p_fibers_d.size() == 0 && layer_4_p_fibers_d.size() == 0)){
    // if(to_print) cout << "Need 1 hit in each helicity." << endl;
    return ZZ;
    //return -10000;
  }

  bool flag = false;
  if(layer_1_p_fibers_d.size() > sft_hit_count_max){
    if(to_print)
      //cout << "More than " << sft_hit_count_max << " SFT hits in layer 1" << endl;
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
    return ZZ;
    //return -10000; 

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

  //double z_avg_min = 0;
  //double delta_z_min = 1000;
  //int delta_z_min_index = 0;

  vector<int> good_delta_z_index;

  bool z_printed = false;


  double delta_z_thr = 10.0;  // Delta_Z

  double sft_delta_z_selected = 1000;
  //double sft_z_selected = 0;

  for(int i = 0; i < iteration_counter; i++){
    if(delta_z_arr[i] < delta_z_thr){
      good_delta_z_index.push_back(i);
    }
  }

  //if(to_print)
    //cout<<endl;

  if(to_print){
    //cout << "Event #     Angle      L1        L2        L3        L4    ||   Z Avg      Delta Z      Track Length" << endl;
    //cout << "===========================================================||========================================" << endl;
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


    //double sft_angle = atan(fabs(C2X_centroid - z_avg/double(z_avg_count))/(R_C2-R_SFT_L1));
    //double total_length = len_in_target/cos(sft_angle);

    if(delta_z_arr[*it] <= sft_delta_z_selected && fabs(z_avg/double(z_avg_count)) <=100){
      if(use_TOF1){
        if(fabs(z_avg/double(z_avg_count) - TOF1_pos) < TOF1_Z_cut){
          //sft_z_selected = z_avg/double(z_avg_count);
          sft_delta_z_selected = delta_z_arr[*it];
        }
      }
      else{
        //sft_z_selected = z_avg/double(z_avg_count);
        sft_delta_z_selected = delta_z_arr[*it];
      }
    }

    if(p1_p2_flag[*it] > 1 && p3_p4_flag[*it] > 1){
      if(ZZ.empty()){
        z_printed = false;
      }
      else{
        for(vector<double>::iterator db = ZZ.begin(); db != ZZ.end(); db++){
          if(fabs(*db - z_avg/double(z_avg_count)) < 0.01){
            z_printed = true;
          }
        }
      }

      if(to_print){
        if(fabs(z_avg/double(z_avg_count)) <=100 && !z_printed){

          if(p1_store_arr[*it] >0 && p2_store_arr[*it] > 0 && p3_store_arr[*it] > 0 && p4_store_arr[*it] > 0){
            //printf("%-5d   %10.2f   %5.1f     %5.1f     %5.1f     %5.1f   ||\n", evt, phi, p1_store_arr[*it],p2_store_arr[*it],
            //  p3_store_arr[*it], p4_store_arr[*it]);
            //printf("N(i)               %5d     %5d     %5d     %5d     ||\n", N1_store_arr[*it],N2_store_arr[*it],N3_store_arr[*it],N4_store_arr[*it]);
            //printf("Z(i)               %7.2f   %7.2f    %7.2f   %7.2f  || %7.2f      %5.3f            %5.2f\n", z1_store_arr[*it],z2_store_arr[*it],z3_store_arr[*it],z4_store_arr[*it],
             //z_avg/double(z_avg_count),delta_z_arr[*it], total_length);
            //DZ.push_back(delta_z_arr[*it]);
            //cout << "===========================================================||======================================== " << endl;
          }
          else if(p1_store_arr[*it] <=0 && p2_store_arr[*it] > 0 && p3_store_arr[*it] > 0 && p4_store_arr[*it] > 0){
            //printf("%-5d   %10.2f %5c       %5.1f     %5.1f     %5.1f   ||\n", evt, phi, '-',p2_store_arr[*it],
            //  p3_store_arr[*it], p4_store_arr[*it]);
            //printf("N(i)               %5c     %5d     %5d     %5d     ||\n", '-',N2_store_arr[*it],N3_store_arr[*it],N4_store_arr[*it]);
            //printf("Z(i)             %7c     %7.2f    %7.2f   %7.2f  || %7.2f      %5.1f            %5.2f\n", '-',z2_store_arr[*it],z3_store_arr[*it],z4_store_arr[*it],
            // z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            //cout << "===========================================================||======================================== " << endl; 
          }
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] <= 0 && p3_store_arr[*it] > 0 && p4_store_arr[*it] > 0){
            //printf("%-5d   %10.2f   %5.1f   %5c       %5.1f     %5.1f   ||\n", evt, phi, p1_store_arr[*it],'-',
            //  p3_store_arr[*it], p4_store_arr[*it]);
            //printf("N(i)               %5d     %5c     %5d     %5d     ||\n", N1_store_arr[*it],'-',N3_store_arr[*it],N4_store_arr[*it]);
            //printf("Z(i)               %7.2f %7c      %7.2f   %7.2f  || %7.2f      %5.1f            %5.2f\n", z1_store_arr[*it],'-',z3_store_arr[*it],z4_store_arr[*it],
            // z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            //cout << "===========================================================||======================================== " << endl;
          }
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] > 0 && p3_store_arr[*it] <= 0 && p4_store_arr[*it] > 0){
            //printf("%-5d   %10.2f   %5.1f     %5.1f   %5c       %5.1f   ||\n", evt, phi, p1_store_arr[*it],p2_store_arr[*it],
            //  '-', p4_store_arr[*it]);
            //printf("N(i)               %5d     %5d     %5c     %5d     ||\n", N1_store_arr[*it],N2_store_arr[*it],'-',N4_store_arr[*it]);
            //printf("Z(i)               %7.2f   %7.2f %7c      %7.2f  || %7.2f      %5.1f            %5.2f\n", z1_store_arr[*it],z2_store_arr[*it],'-',z4_store_arr[*it],
            // z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            //cout << "===========================================================||======================================== " << endl;
          }
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] > 0 && p3_store_arr[*it] > 0 && p4_store_arr[*it] <= 0){
            //printf("%-5d   %10.2f   %5.1f     %5.1f     %5.1f   %5c     ||\n", evt, phi, p1_store_arr[*it],p2_store_arr[*it],
            //  p3_store_arr[*it], '-');
            //printf("N(i)               %5d     %5d     %5d     %5c     ||\n", N1_store_arr[*it],N2_store_arr[*it],N3_store_arr[*it],'-');
            //printf("Z(i)               %7.2f   %7.2f    %7.2f%7c     || %7.2f      %5.1f            %5.2f\n", z1_store_arr[*it],z2_store_arr[*it],z3_store_arr[*it],'-',
            // z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            //cout << "===========================================================||======================================== " << endl;
          }
          else if(p1_store_arr[*it] <= 0 && p2_store_arr[*it] >0 && p3_store_arr[*it] <= 0 && p4_store_arr[*it] >0){
            //printf("%-5d   %10.2f %5c       %5.1f   %5c       %5.1f   ||\n", evt, phi, '-',p2_store_arr[*it],
            //  '-', p4_store_arr[*it]);
            //printf("N(i)               %5c     %5d     %5c     %5d     ||\n",'-',N2_store_arr[*it],'-',N4_store_arr[*it]);
            //printf("Z(i)             %7c     %7.2f %7c      %7.2f  || %7.2f      %5.1f            %5.2f\n", '-',z2_store_arr[*it],'-',z4_store_arr[*it],
            // z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            //cout << "===========================================================||======================================== " << endl;
          }
          else if(p1_store_arr[*it] <= 0 && p2_store_arr[*it] >0 && p3_store_arr[*it] >0 && p4_store_arr[*it] <= 0){
            //printf("%-5d   %10.2f %5c       %5.1f     %5.1f   %5c     ||\n", evt, phi, '-',p2_store_arr[*it],
            //  p3_store_arr[*it], '-');
            //printf("N(i)               %5c     %5d     %5d     %5c     ||\n", '-',N2_store_arr[*it],N3_store_arr[*it],'-');
            //printf("Z(i)             %7c     %7.2f    %7.2f%7c     || %7.2f      %5.1f            %5.2f\n", '-',z2_store_arr[*it],z3_store_arr[*it],'-',
            // z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            //cout << "===========================================================||======================================== " << endl;
          }
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] <= 0 && p3_store_arr[*it] <= 0 && p4_store_arr[*it] >0){
            //printf("%-5d   %10.2f   %5.1f   %5c     %5c       %5.1f   ||\n", evt, phi, p1_store_arr[*it],'-',
            //  '-', p4_store_arr[*it]);
            //printf("N(i)               %5d     %5c     %5c     %5d     ||\n", N1_store_arr[*it],'-','-',N4_store_arr[*it]);
            //printf("Z(i)               %7.2f %7c   %7c      %7.2f  || %7.2f      %5.1f            %5.2f\n", z1_store_arr[*it],'-','-',z4_store_arr[*it],
             //z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            //cout << "===========================================================||======================================== " << endl;
          }
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] <= 0 && p3_store_arr[*it] >0 && p4_store_arr[*it] <= 0){
            //printf("%-5d   %10.2f   %5.1f   %5c       %5.1f   %5c     ||\n", evt, phi, p1_store_arr[*it],'-',
            //  p3_store_arr[*it], '-');
            //printf("N(i)               %5d     %5c     %5d     %5c     ||\n", N1_store_arr[*it],'-',N3_store_arr[*it],'-');
            //printf("Z(i)               %7.2f %7c      %7.2f%7c     || %7.2f      %5.1f            %5.2f\n", z1_store_arr[*it],'-',z3_store_arr[*it],'-',
            // z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            //cout << "===========================================================||======================================== " << endl;
            
          }
          ZZ.push_back(z_avg/double(z_avg_count));
        }   
      }
    }
  }

  //for(unsigned i=0; i<z_avgs_store.size(); i++){
  //  cout << "UUUU :   " << z_avgs_store[i] << endl; 
  //}

  //cout << endl;
  return ZZ;
 }





























vector<double> Z_list2(double ADC_High_SFT_corr[128], int has_TDC_SFT_hit[128], int SFT_channel_to_fiber[128],
 int evt, double phi, bool to_print, double TOF1_pos, bool use_TOF1, double C2X_centroid, double len_in_target){

  //cout << endl;
  //cout << endl;
  //cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
  //cout << endl;
  //cout << endl;

  double z_avg = 0.0;

  const int sft_hit_count_max = 3;
  const int NUM_ITERATIONS = 100;

  //float R_TARGET = 29.0;
  //float R_TOF1 = 47.1;
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

  /*  cout << endl;
  cout << endl;
  cout << "####### DEBUG SFT print  ########" << endl;
  cout << "p1 = ";
  for (unsigned i1=0; i1<layer_1_p_fibers.size(); i1++) cout << layer_1_p_fibers[i1] << "  ";
  cout << endl; 
  cout << "p2 = "; 
  for (unsigned i2=0; i2<layer_2_p_fibers.size(); i2++) cout << layer_2_p_fibers[i2] << "  ";
  cout << endl;
  cout << "p3 = ";
  for (unsigned i3=0; i3<layer_3_p_fibers.size(); i3++) cout << layer_3_p_fibers[i3] << "  ";
  cout << endl;
  cout << "p4 = ";
  for (unsigned i4=0; i4<layer_4_p_fibers.size(); i4++) cout << layer_4_p_fibers[i4] << "  ";
  cout << endl;
  cout << "###########################################" << endl;
  cout << endl;
  cout << endl;
  */

  // pass vectors to new p-averaging function. return vector<double> of averaged p-fibers.
  // vector<double> Layer_averaging(vector<int> p_fibers)
  vector<double> layer_1_p_fibers_d;
  vector<double> layer_2_p_fibers_d;
  vector<double> layer_3_p_fibers_d;
  vector<double> layer_4_p_fibers_d;
  
  //bool bad_sft = false;
  //bool *bad_sft_flag = &bad_sft;

  layer_1_p_fibers_d = Layer_averaging(layer_1_p_fibers, 1);
  layer_2_p_fibers_d = Layer_averaging(layer_2_p_fibers, 2);
  layer_3_p_fibers_d = Layer_averaging(layer_3_p_fibers, 3);
  layer_4_p_fibers_d = Layer_averaging(layer_4_p_fibers, 4);
     

  // Test if more than 1 layer is empty or if any layer has more than 3 p fibers
  if((layer_1_p_fibers_d.size() == 0 && layer_2_p_fibers_d.size() == 0) ||  
   (layer_3_p_fibers_d.size() == 0 && layer_4_p_fibers_d.size() == 0)){
    if(to_print){
      //cout << "Need 1 hit in each helicity." << endl;
      z_avg = -10000;
      //return -1000
    }
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
    z_avg = -10000;
    //return -10000; 

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

  //double z_avg_min = 0;
  //double delta_z_min = 1000;
  //int delta_z_min_index = 0;

  vector<int> good_delta_z_index;

  vector<double> z_avgs_store;
  vector<double> DeltaZ;
  bool z_printed = false;


  double delta_z_thr = 10.0;  // Delta_Z

  double sft_delta_z_selected = 1000;
  //double sft_z_selected = 0;

  for(int i = 0; i < iteration_counter; i++){
    if(delta_z_arr[i] < delta_z_thr){
      good_delta_z_index.push_back(i);
    }
  }

  if(to_print)
    cout<<endl;

  //if(to_print){
  //  cout << "Event #     Angle      L1        L2        L3        L4    ||   Z Avg      Delta Z      Track Length" << endl;
  //  cout << "===========================================================||========================================" << endl;
  //}

  for(vector<int>::iterator it = good_delta_z_index.begin(); it != good_delta_z_index.end(); it++){
  

    
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
          //sft_z_selected = z_avg/double(z_avg_count);
          sft_delta_z_selected = delta_z_arr[*it];
        }
      }
      else{
        //sft_z_selected = z_avg/double(z_avg_count);
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

      cout << "ALBANE :  " << z_avg << endl;


      if(to_print){
        if(fabs(z_avg/double(z_avg_count)) <=100 && !z_printed){

          if(p1_store_arr[*it] >0 && p2_store_arr[*it] > 0 && p3_store_arr[*it] > 0 && p4_store_arr[*it] > 0){
            printf("%-5d   %10.2f   %5.1f     %5.1f     %5.1f     %5.1f   ||\n", evt, phi, p1_store_arr[*it],p2_store_arr[*it],
              p3_store_arr[*it], p4_store_arr[*it]);
            printf("N(i)               %5d     %5d     %5d     %5d     ||\n", N1_store_arr[*it],N2_store_arr[*it],N3_store_arr[*it],N4_store_arr[*it]);
            printf("Z(i)               %7.2f   %7.2f    %7.2f   %7.2f  || %7.2f      %5.3f            %5.2f\n", z1_store_arr[*it],z2_store_arr[*it],z3_store_arr[*it],z4_store_arr[*it],
             z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            cout << "===========================================================||======================================== " << endl;
            //cout << "INTER:  " << z_avg/double(z_avg_count) << endl;
            //cout << "SEB: " << z_avg/double(z_avg_count) << endl;
            //z_avgs_store.push_back(z_avg/double(z_avg_count));
            //DeltaZ.push_back(delta_z_arr[*it]);
          }
          
          else if(p1_store_arr[*it] <=0 && p2_store_arr[*it] > 0 && p3_store_arr[*it] > 0 && p4_store_arr[*it] > 0){
            printf("%-5d   %10.2f %5c       %5.1f     %5.1f     %5.1f   ||\n", evt, phi, '-',p2_store_arr[*it],
              p3_store_arr[*it], p4_store_arr[*it]);
            printf("N(i)               %5c     %5d     %5d     %5d     ||\n", '-',N2_store_arr[*it],N3_store_arr[*it],N4_store_arr[*it]);
            printf("Z(i)             %7c     %7.2f    %7.2f   %7.2f  || %7.2f      %5.1f            %5.2f\n", '-',z2_store_arr[*it],z3_store_arr[*it],z4_store_arr[*it],
             z_avg/double(z_avg_count),delta_z_arr[*it], total_length);
            cout << "===========================================================||======================================== " << endl; 
            //cout << "INTER:  " << z_avg/double(z_avg_count) << endl;
            //cout << "SEB: " << z_avg/double(z_avg_count) << endl;
            //z_avgs_store.push_back(z_avg/double(z_avg_count));
            //DeltaZ.push_back(delta_z_arr[*it]);
          }
          
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] <= 0 && p3_store_arr[*it] > 0 && p4_store_arr[*it] > 0){
            printf("%-5d   %10.2f   %5.1f   %5c       %5.1f     %5.1f   ||\n", evt, phi, p1_store_arr[*it],'-',
              p3_store_arr[*it], p4_store_arr[*it]);
            printf("N(i)               %5d     %5c     %5d     %5d     ||\n", N1_store_arr[*it],'-',N3_store_arr[*it],N4_store_arr[*it]);
            printf("Z(i)               %7.2f %7c      %7.2f   %7.2f  || %7.2f      %5.1f            %5.2f\n", z1_store_arr[*it],'-',z3_store_arr[*it],z4_store_arr[*it],
             z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            cout << "===========================================================||======================================== " << endl;
            //cout << "INTER:  " << z_avg/double(z_avg_count) << endl;
            //cout << "SEB: " << z_avg/double(z_avg_count) << endl;
            //z_avgs_store.push_back(z_avg/double(z_avg_count));
            //DeltaZ.push_back(delta_z_arr[*it]);
          }
          
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] > 0 && p3_store_arr[*it] <= 0 && p4_store_arr[*it] > 0){
            printf("%-5d   %10.2f   %5.1f     %5.1f   %5c       %5.1f   ||\n", evt, phi, p1_store_arr[*it],p2_store_arr[*it],
              '-', p4_store_arr[*it]);
            printf("N(i)               %5d     %5d     %5c     %5d     ||\n", N1_store_arr[*it],N2_store_arr[*it],'-',N4_store_arr[*it]);
            printf("Z(i)               %7.2f   %7.2f %7c      %7.2f  || %7.2f      %5.1f            %5.2f\n", z1_store_arr[*it],z2_store_arr[*it],'-',z4_store_arr[*it],
             z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            cout << "===========================================================||======================================== " << endl;
            //cout << "INTER:  " << z_avg/double(z_avg_count) << endl;
            //cout << "SEB: " << z_avg/double(z_avg_count) << endl;
            //z_avgs_store.push_back(z_avg/double(z_avg_count));
            //DeltaZ.push_back(delta_z_arr[*it]);
          }
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] > 0 && p3_store_arr[*it] > 0 && p4_store_arr[*it] <= 0){
          
            printf("%-5d   %10.2f   %5.1f     %5.1f     %5.1f   %5c     ||\n", evt, phi, p1_store_arr[*it],p2_store_arr[*it],
              p3_store_arr[*it], '-');
            printf("N(i)               %5d     %5d     %5d     %5c     ||\n", N1_store_arr[*it],N2_store_arr[*it],N3_store_arr[*it],'-');
            printf("Z(i)               %7.2f   %7.2f    %7.2f%7c     || %7.2f      %5.1f            %5.2f\n", z1_store_arr[*it],z2_store_arr[*it],z3_store_arr[*it],'-',
             z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            cout << "===========================================================||======================================== " << endl;
            //cout << "INTER:  " << z_avg/double(z_avg_count) << endl;
            //cout << "SEB: " << z_avg/double(z_avg_count) << endl;
            //z_avgs_store.push_back(z_avg/double(z_avg_count));
            //DeltaZ.push_back(delta_z_arr[*it]);
          }
          
          else if(p1_store_arr[*it] <= 0 && p2_store_arr[*it] >0 && p3_store_arr[*it] <= 0 && p4_store_arr[*it] >0){
            printf("%-5d   %10.2f %5c       %5.1f   %5c       %5.1f   ||\n", evt, phi, '-',p2_store_arr[*it],
              '-', p4_store_arr[*it]);
            printf("N(i)               %5c     %5d     %5c     %5d     ||\n",'-',N2_store_arr[*it],'-',N4_store_arr[*it]);
            printf("Z(i)             %7c     %7.2f %7c      %7.2f  || %7.2f      %5.1f            %5.2f\n", '-',z2_store_arr[*it],'-',z4_store_arr[*it],
             z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            cout << "===========================================================||======================================== " << endl;
            //cout << "INTER:  " << z_avg/double(z_avg_count) << endl;
            //cout << "SEB: " << z_avg/double(z_avg_count) << endl;
            //z_avgs_store.push_back(z_avg/double(z_avg_count));
            //DeltaZ.push_back(delta_z_arr[*it]);
          }
          
          else if(p1_store_arr[*it] <= 0 && p2_store_arr[*it] >0 && p3_store_arr[*it] >0 && p4_store_arr[*it] <= 0){
            printf("%-5d   %10.2f %5c       %5.1f     %5.1f   %5c     ||\n", evt, phi, '-',p2_store_arr[*it],
              p3_store_arr[*it], '-');
            printf("N(i)               %5c     %5d     %5d     %5c     ||\n", '-',N2_store_arr[*it],N3_store_arr[*it],'-');
            printf("Z(i)             %7c     %7.2f    %7.2f%7c     || %7.2f      %5.1f            %5.2f\n", '-',z2_store_arr[*it],z3_store_arr[*it],'-',
             z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            cout << "===========================================================||======================================== " << endl;
            //cout << "INTER:  " << z_avg/double(z_avg_count) << endl;
            //cout << "SEB: " << z_avg/double(z_avg_count) << endl;
            //z_avgs_store.push_back(z_avg/double(z_avg_count));
            //DeltaZ.push_back(delta_z_arr[*it]);
          }
          
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] <= 0 && p3_store_arr[*it] <= 0 && p4_store_arr[*it] >0){
            printf("%-5d   %10.2f   %5.1f   %5c     %5c       %5.1f   ||\n", evt, phi, p1_store_arr[*it],'-',
              '-', p4_store_arr[*it]);
            printf("N(i)               %5d     %5c     %5c     %5d     ||\n", N1_store_arr[*it],'-','-',N4_store_arr[*it]);
            printf("Z(i)               %7.2f %7c   %7c      %7.2f  || %7.2f      %5.1f            %5.2f\n", z1_store_arr[*it],'-','-',z4_store_arr[*it],
             z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            cout << "===========================================================||======================================== " << endl;
            //cout << "INTER:  " << z_avg/double(z_avg_count) << endl;
            //cout << "SEB: " << z_avg/double(z_avg_count) << endl;
            //z_avgs_store.push_back(z_avg/double(z_avg_count));
            //DeltaZ.push_back(delta_z_arr[*it]);
          }
          
          else if(p1_store_arr[*it] >0 && p2_store_arr[*it] <= 0 && p3_store_arr[*it] >0 && p4_store_arr[*it] <= 0){
            printf("%-5d   %10.2f   %5.1f   %5c       %5.1f   %5c     ||\n", evt, phi, p1_store_arr[*it],'-',
              p3_store_arr[*it], '-');
            printf("N(i)               %5d     %5c     %5d     %5c     ||\n", N1_store_arr[*it],'-',N3_store_arr[*it],'-');
            printf("Z(i)               %7.2f %7c      %7.2f%7c     || %7.2f  I    %5.1f            %5.2f\n", z1_store_arr[*it],'-',z3_store_arr[*it],'-',
             z_avg/double(z_avg_count),delta_z_arr[*it], total_length);

            cout << "===========================================================||======================================== " << endl;        
            //cout << "INTER:  " << z_avg/double(z_avg_count) << endl;
            //cout << "SEB: " << z_avg/double(z_avg_count) << endl;
            //z_avgs_store.push_back(z_avg/double(z_avg_count));
            //DeltaZ.push_back(delta_z_arr[*it]);
          }
          z_avgs_store.push_back(z_avg/double(z_avg_count));

        }   
      }
    }
  }


  //sft_z_selected = z_avgs_store[3];

  for(unsigned i=0; i<z_avgs_store.size(); i++){
    cout << "SSSSS:  " << z_avgs_store[i] << endl;
  }

  return z_avgs_store;}


double Ana2(double ADC_High_SFT_corr[128], int has_TDC_SFT_hit[128], int SFT_channel_to_fiber[128],
 int evt, double phi, bool to_print, double TOF1_pos, bool use_TOF1, double C2X_centroid, double len_in_target){

  const int sft_hit_count_max = 3;
  const int NUM_ITERATIONS = 100;

  //float R_TARGET = 29.0;
  //float R_TOF1 = 47.1;
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

  cout << endl;
  cout << endl;
  cout << "####### DEBUG SFT print  ########" << endl;
  cout << "p1 = ";
  for (unsigned i1=0; i1<layer_1_p_fibers.size(); i1++) cout << layer_1_p_fibers[i1] << "  ";
  cout << endl; 
  cout << "p2 = "; 
  for (unsigned i2=0; i2<layer_2_p_fibers.size(); i2++) cout << layer_2_p_fibers[i2] << "  ";
  cout << endl;
  cout << "p3 = ";
  for (unsigned i3=0; i3<layer_3_p_fibers.size(); i3++) cout << layer_3_p_fibers[i3] << "  ";
  cout << endl;
  cout << "p4 = ";
  for (unsigned i4=0; i4<layer_4_p_fibers.size(); i4++) cout << layer_4_p_fibers[i4] << "  ";
  cout << endl;
  cout << "###########################################" << endl;
  cout << endl;
  cout << endl;
  // pass vectors to new p-averaging function. return vector<double> of averaged p-fibers.
  // vector<double> Layer_averaging(vector<int> p_fibers)
  vector<double> layer_1_p_fibers_d;
  vector<double> layer_2_p_fibers_d;
  vector<double> layer_3_p_fibers_d;
  vector<double> layer_4_p_fibers_d;
  
  //bool bad_sft = false;
  //bool *bad_sft_flag = &bad_sft;

  layer_1_p_fibers_d = Layer_averaging(layer_1_p_fibers, 1);
  layer_2_p_fibers_d = Layer_averaging(layer_2_p_fibers, 2);
  layer_3_p_fibers_d = Layer_averaging(layer_3_p_fibers, 3);
  layer_4_p_fibers_d = Layer_averaging(layer_4_p_fibers, 4);
     

  // Test if more than 1 layer is empty or if any layer has more than 3 p fibers
  if((layer_1_p_fibers_d.size() == 0 && layer_2_p_fibers_d.size() == 0) ||  
   (layer_3_p_fibers_d.size() == 0 && layer_4_p_fibers_d.size() == 0)){
    if(to_print)
      //cout << "Need 1 hit in each helicity." << endl;
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

  //double z_avg_min = 0;
  //double delta_z_min = 1000;
  //int delta_z_min_index = 0;

  vector<int> good_delta_z_index;

  vector<double> z_avgs_store;
  bool z_printed = false;


  double delta_z_thr = 10.0;  // Delta_Z

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

void Ana3(int Run_Number, int evt, double ADC_High_SFT_corr[128], int has_TDC_SFT_hit[128], int SFT_channel_to_fiber[128], double phi, int gap_to_fit_left, int gap_to_fit_right, double C2_Z_top, double C2_Z_bottom){

  cout << "###############################    C2_Z-SFT_Z ANALYSIS    ############################### " << endl;
  cout << endl;
  //cout << "TEST :  " << phi <<  endl;
  //cout << "DEBUG :  " << C2_Z_top << "  " << C2_Z_bottom << endl;
  //cout << "  " << endl;
  //for(int i=0; i<128; i++){
  //  if(ADC_High_SFT_corr[i]!=0) cout << "DEBUG99:  " << evt << "  " << i << "  " << ADC_High_SFT_corr[i] << endl;
  //}


  //const int sft_hit_count_max = 3;
  //const int NUM_ITERATIONS = 100;

  //float R_TARGET = 29.0;
  //float R_TOF1 = 47.1;
  //float R_SFT_L1 = 40.0;
  //float R_C2 = 629.4;

  int L1_sft_count = 0;
  int L2_sft_count = 0;
  int L3_sft_count = 0;
  int L4_sft_count = 0;

  double phi_top = 0.;
  double phi_bottom = 0.;

  vector<int> layer_1_p_fibers;   layer_1_p_fibers.clear();
  vector<int> layer_2_p_fibers;   layer_2_p_fibers.clear();
  vector<int> layer_3_p_fibers;   layer_3_p_fibers.clear();
  vector<int> layer_4_p_fibers;   layer_4_p_fibers.clear();

  vector<double> layer_1_adc;   layer_1_adc.clear();
  vector<double> layer_2_adc;   layer_2_adc.clear();
  vector<double> layer_3_adc;   layer_3_adc.clear();
  vector<double> layer_4_adc;   layer_4_adc.clear();

  vector<double> layer_1_p_fibers_d;  layer_1_p_fibers_d.clear();
  vector<double> layer_2_p_fibers_d;  layer_2_p_fibers_d.clear();
  vector<double> layer_3_p_fibers_d;  layer_3_p_fibers_d.clear();
  vector<double> layer_4_p_fibers_d;  layer_4_p_fibers_d.clear();
  
  vector< vector<double> > SFT_Z1_top;  vector<double> SFT_Z1_top_col;  
  vector< vector<double> > SFT_Z2_top;  vector<double> SFT_Z2_top_col;  
  vector< vector<double> > SFT_Z3_top;  vector<double> SFT_Z3_top_col;  
  vector< vector<double> > SFT_Z4_top;  vector<double> SFT_Z4_top_col;  
  SFT_Z1_top.clear();  SFT_Z1_top_col.clear();  
  SFT_Z2_top.clear();  SFT_Z2_top_col.clear();  
  SFT_Z3_top.clear();  SFT_Z3_top_col.clear();  
  SFT_Z4_top.clear();  SFT_Z4_top_col.clear();  

  vector< vector<double> > SFT_Z1_bottom;  vector<double> SFT_Z1_bottom_col;  
  vector< vector<double> > SFT_Z2_bottom;  vector<double> SFT_Z2_bottom_col;  
  vector< vector<double> > SFT_Z3_bottom;  vector<double> SFT_Z3_bottom_col;  
  vector< vector<double> > SFT_Z4_bottom;  vector<double> SFT_Z4_bottom_col;  
  SFT_Z1_bottom.clear();  SFT_Z1_bottom_col.clear();  
  SFT_Z2_bottom.clear();  SFT_Z2_bottom_col.clear();  
  SFT_Z3_bottom.clear();  SFT_Z3_bottom_col.clear();  
  SFT_Z4_bottom.clear();  SFT_Z4_bottom_col.clear();  

  vector< vector<double> > DZ1_top;  vector<double> DZ1_top_col;  
  vector< vector<double> > DZ2_top;  vector<double> DZ2_top_col;  
  vector< vector<double> > DZ3_top;  vector<double> DZ3_top_col;  
  vector< vector<double> > DZ4_top;  vector<double> DZ4_top_col;  
  DZ1_top.clear();  DZ1_top_col.clear();  
  DZ2_top.clear();  DZ2_top_col.clear();  
  DZ3_top.clear();  DZ3_top_col.clear();  
  DZ4_top.clear();  DZ4_top_col.clear();  

  vector< vector<double> > DZ1_bottom;  vector<double> DZ1_bottom_col;  
  vector< vector<double> > DZ2_bottom;  vector<double> DZ2_bottom_col;  
  vector< vector<double> > DZ3_bottom;  vector<double> DZ3_bottom_col;  
  vector< vector<double> > DZ4_bottom;  vector<double> DZ4_bottom_col;  
  DZ1_bottom.clear();  DZ1_bottom_col.clear();  
  DZ2_bottom.clear();  DZ2_bottom_col.clear();  
  DZ3_bottom.clear();  DZ3_bottom_col.clear();  
  DZ4_bottom.clear();  DZ4_bottom_col.clear();  

  vector<double> DZ1_lowest_top;  vector<double> DZ1_lowest_bottom;
  vector<double> DZ2_lowest_top;  vector<double> DZ2_lowest_bottom;
  vector<double> DZ3_lowest_top;  vector<double> DZ3_lowest_bottom;
  vector<double> DZ4_lowest_top;  vector<double> DZ4_lowest_bottom;
  DZ1_lowest_top.clear();  DZ1_lowest_bottom.clear();
  DZ2_lowest_top.clear();  DZ2_lowest_bottom.clear();
  DZ3_lowest_top.clear();  DZ3_lowest_bottom.clear();
  DZ4_lowest_top.clear();  DZ4_lowest_bottom.clear();

  vector<int> n1_lowest_top;      vector<int> n1_lowest_bottom;
  vector<int> n2_lowest_top;      vector<int> n2_lowest_bottom;
  vector<int> n3_lowest_top;      vector<int> n3_lowest_bottom;
  vector<int> n4_lowest_top;      vector<int> n4_lowest_bottom;
  n1_lowest_top.clear();      n1_lowest_bottom.clear();
  n2_lowest_top.clear();      n2_lowest_bottom.clear();
  n3_lowest_top.clear();      n3_lowest_bottom.clear();
  n4_lowest_top.clear();      n4_lowest_bottom.clear();

  vector<int> p1_lowest_top;      vector<int> p1_lowest_bottom;
  vector<int> p2_lowest_top;      vector<int> p2_lowest_bottom;
  vector<int> p3_lowest_top;      vector<int> p3_lowest_bottom;
  vector<int> p4_lowest_top;      vector<int> p4_lowest_bottom;
  p1_lowest_top.clear();      p1_lowest_bottom.clear();
  p2_lowest_top.clear();      p2_lowest_bottom.clear();
  p3_lowest_top.clear();      p3_lowest_bottom.clear();
  p4_lowest_top.clear();      p4_lowest_bottom.clear();

  vector<double> Z1_lowest_top;      vector<double> Z1_lowest_bottom;
  vector<double> Z2_lowest_top;      vector<double> Z2_lowest_bottom;
  vector<double> Z3_lowest_top;      vector<double> Z3_lowest_bottom;
  vector<double> Z4_lowest_top;      vector<double> Z4_lowest_bottom;
  Z1_lowest_top.clear();      Z1_lowest_bottom.clear();
  Z2_lowest_top.clear();      Z2_lowest_bottom.clear();
  Z3_lowest_top.clear();      Z3_lowest_bottom.clear();
  Z4_lowest_top.clear();      Z4_lowest_bottom.clear();


  double lowest_top = 999.;   double lowest_bottom = 999.;
  int n_lowest_top = 99;      int n_lowest_bottom = 99;
  double z_lowest_top = 999.; double z_lowest_bottom = 999.;
  //int n_lowest_top[4]={0};    int n_lowest_bottom[4] = {0};

  bool L1_sft_flag = false;
  bool L2_sft_flag = false;
  bool L3_sft_flag = false;
  bool L4_sft_flag = false;

  layer_1_p_fibers.clear();
  layer_2_p_fibers.clear();
  layer_3_p_fibers.clear();
  layer_4_p_fibers.clear();

  //cout << "DEBUG :   Gap_to_fit_LEFT : " <<  gap_to_fit_left << "  gap_to_fit_RIGHT : " << gap_to_fit_right << endl;
  //cout << endl;

  //  for(int i=Run_Number; i<Run_Number+1; i++){
  //    if(gap_to_fit_left!=12){
  //    cout << "Nope !" << endl;
  //    break;
      
  //  }

    //if(gap_to_fit_left==11 && gap_to_fit_right!=5) break; 
    //if(gap_to_fit_left==10 && gap_to_fit_right!=4) break;   
    //if(gap_to_fit_left==8 && gap_to_fit_right!=2) break;   
    //if(gap_to_fit_left==7 && gap_to_fit_right!=1) break;   

    if((gap_to_fit_left>=7 && gap_to_fit_left<=12) && (gap_to_fit_right>=1 && gap_to_fit_right<=6)){
      if(phi > 180){
        phi_top = phi-180;
        phi_bottom = phi; 
      }
      else{
        phi_top = phi;
        phi_bottom = phi+180;
      }
    }
    /*
    if(gap_to_fit_left==7 && gap_to_fit_right==2 ){
      if(phi > 270){
        phi_top = phi-180;
        phi_bottom = phi;
      }
      else{
        phi_top = phi;
        phi_bottom = phi+180;
      }   
    }
    if(gap_to_fit_left==7 && gap_to_fit_right==1){
      if(phi > 270){
        phi_top = phi-180;
        phi_bottom = phi;
      }
      else{
        phi_top = phi;
        phi_bottom = phi+180;
      }   
    }
    */ 

    //cout << "TEST: " << phi_top << "  " << phi_bottom << endl; 

    for(int i = 0; i <128; i++){
      if(ADC_High_SFT_corr[i] != 0 && has_TDC_SFT_hit[i]){
      //if(ADC_High_SFT_corr[i] != 0){

        //cout << "DEBUG2:  " << evt << "  " << i << "  " << ADC_High_SFT_corr[i] << endl;
        //cout << "DEBUG2:  " << evt << "  " << i << "  " << ADC_High_SFT_corr[i] << endl;

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

  for(unsigned int i=0; i<layer_2_p_fibers.size(); i++){
  //cout << "DEBUG1:" << evt << "  " << i << "   " << layer_2_p_fibers[i] << endl;
  //cout << "DEBUG1:" << evt << "  " << i < "   " << layer_2_p_fibers[i] << endl;
  }

  //bool bad_sft = false;
  //bool *bad_sft_flag = &bad_sft;

  layer_1_p_fibers_d.clear();
  layer_2_p_fibers_d.clear();
  layer_3_p_fibers_d.clear();
  layer_4_p_fibers_d.clear();

  //layer_1_p_fibers_d = Layer_averaging(layer_1_p_fibers, bad_sft_flag, 1);
  //layer_2_p_fibers_d = Layer_averaging(layer_2_p_fibers, bad_sft_flag, 2);
  //layer_3_p_fibers_d = Layer_averaging(layer_3_p_fibers, bad_sft_flag, 3);
  //layer_4_p_fibers_d = Layer_averaging(layer_4_p_fibers, bad_sft_flag, 4);

  for(unsigned int i1=0; i1<layer_1_p_fibers.size(); i1++) layer_1_p_fibers_d.push_back(double(layer_1_p_fibers[i1]));
  for(unsigned int i2=0; i2<layer_2_p_fibers.size(); i2++) layer_2_p_fibers_d.push_back(double(layer_2_p_fibers[i2]));
  for(unsigned int i3=0; i3<layer_3_p_fibers.size(); i3++) layer_3_p_fibers_d.push_back(double(layer_3_p_fibers[i3]));
  for(unsigned int i4=0; i4<layer_4_p_fibers.size(); i4++) layer_4_p_fibers_d.push_back(double(layer_4_p_fibers[i4]));

  //for(unsigned int i2=0; i2<layer_2_p_fibers.size(); i2++){
  //  layer_2_p_fibers_d[i2].push_back(double(layer_2_p_fibers[i2]);
  //}
  //for(unsigned int i3=0; i3<layer_3_p_fibers.size(); i3++){
  //  layer_3_p_fibers_d[i3].push_back(double(layer_3_p_fibers[i3]);
  //}
  //for(unsigned int i4=0; i4<layer_4_p_fibers.size(); i4++){
  //  layer_4_p_fibers_d[i4].push_back(double(layer_4_p_fibers[i4]);
  //}

  /*
  cout << endl;
  cout << endl;
  cout << "####### Ana 3  ########" << endl;
  cout << "phi (TOP) = " << phi_top << " (deg.)  ,  phi (BOTTOM) = " << phi_bottom << " (deg.)" << endl;
  cout << "C2_Z (TOP) = " << C2_Z_top << "  " << "C2_Z (BOTTOM) = " << C2_Z_bottom << endl;
  cout << "p1 = ";
  for (unsigned int i1=0; i1<layer_1_p_fibers_d.size(); i1++) cout << layer_1_p_fibers_d[i1] << "  ";
  cout << endl; 
  cout << "p2 = "; 
  for (unsigned int i2=0; i2<layer_2_p_fibers_d.size(); i2++) cout << layer_2_p_fibers_d[i2] << "  ";
  cout << endl;
  cout << "p3 = ";
  for (unsigned int i3=0; i3<layer_3_p_fibers_d.size(); i3++) cout << layer_3_p_fibers_d[i3] << "  ";
  cout << endl;
  cout << "p4 = ";
  for (unsigned int i4=0; i4<layer_4_p_fibers_d.size(); i4++) cout << layer_4_p_fibers_d[i4] << "  ";
  cout << endl;
  cout << "###########################################" << endl;
  cout << endl;
  cout << endl;
  */

  for(unsigned int i1=0; i1<layer_1_p_fibers_d.size(); i1++){
    for(int j1=0; j1<15; j1++){
      SFT_Z1_top_col.push_back(calculate_Z_position(phi_top, j1+1,layer_1_p_fibers_d[i1], 1));
      SFT_Z1_bottom_col.push_back(calculate_Z_position(phi_bottom, j1+1,layer_1_p_fibers_d[i1], 1));
      DZ1_top_col.push_back(C2_Z_top - calculate_Z_position(phi_top, j1+1,layer_1_p_fibers_d[i1], 1));
      DZ1_bottom_col.push_back(C2_Z_bottom - calculate_Z_position(phi_bottom, j1+1,layer_1_p_fibers_d[i1], 1));
    } 
    SFT_Z1_top.push_back(SFT_Z1_top_col);
    SFT_Z1_bottom.push_back(SFT_Z1_bottom_col);
    DZ1_top.push_back(DZ1_top_col);
    DZ1_bottom.push_back(DZ1_bottom_col);
    SFT_Z1_top_col.clear();
    SFT_Z1_bottom_col.clear();
    DZ1_top_col.clear();
    DZ1_bottom_col.clear();
  }

  for(unsigned int i2=0; i2<layer_2_p_fibers_d.size(); i2++){
    for(int j2=0; j2<15; j2++){
      SFT_Z2_top_col.push_back(calculate_Z_position(phi_top, j2+1,layer_2_p_fibers_d[i2], 2));
      SFT_Z2_bottom_col.push_back(calculate_Z_position(phi_bottom, j2+1,layer_2_p_fibers_d[i2], 2));
      DZ2_top_col.push_back(C2_Z_top - calculate_Z_position(phi_top, j2+1,layer_2_p_fibers_d[i2], 2));
      DZ2_bottom_col.push_back(C2_Z_bottom - calculate_Z_position(phi_bottom, j2+1,layer_2_p_fibers_d[i2], 2));
     } 
    SFT_Z2_top.push_back(SFT_Z2_top_col);
    SFT_Z2_bottom.push_back(SFT_Z2_bottom_col);
    DZ2_top.push_back(DZ2_top_col);
    DZ2_bottom.push_back(DZ2_bottom_col);
    SFT_Z2_top_col.clear();
    SFT_Z2_bottom_col.clear();
    DZ2_top_col.clear();
    DZ2_bottom_col.clear();
  }

  for(unsigned int i3=0; i3<layer_3_p_fibers_d.size(); i3++){
    for(int j3=0; j3<17; j3++){
      SFT_Z3_top_col.push_back(calculate_Z_position(phi_top, j3+1,layer_3_p_fibers_d[i3], 3));
      SFT_Z3_bottom_col.push_back(calculate_Z_position(phi_bottom, j3+1,layer_3_p_fibers_d[i3], 3));
      DZ3_top_col.push_back(C2_Z_top - calculate_Z_position(phi_top, j3+1,layer_3_p_fibers_d[i3], 3));
      DZ3_bottom_col.push_back(C2_Z_bottom - calculate_Z_position(phi_bottom, j3+1,layer_3_p_fibers_d[i3], 3));
     } 
    SFT_Z3_top.push_back(SFT_Z3_top_col);
    SFT_Z3_bottom.push_back(SFT_Z3_bottom_col);
    DZ3_top.push_back(DZ3_top_col);
    DZ3_bottom.push_back(DZ3_bottom_col);
    SFT_Z3_top_col.clear();
    SFT_Z3_bottom_col.clear();
    DZ3_top_col.clear();
    DZ3_bottom_col.clear();
  }

  for(unsigned int i4=0; i4<layer_4_p_fibers_d.size(); i4++){
    for(int j4=0; j4<17; j4++){
      SFT_Z4_top_col.push_back(calculate_Z_position(phi_top, j4+1,layer_4_p_fibers_d[i4], 4));
      SFT_Z4_bottom_col.push_back(calculate_Z_position(phi_bottom, j4+1,layer_4_p_fibers_d[i4], 4));
      DZ4_top_col.push_back(C2_Z_top - calculate_Z_position(phi_top, j4+1,layer_4_p_fibers_d[i4], 4));
      DZ4_bottom_col.push_back(C2_Z_bottom - calculate_Z_position(phi_bottom, j4+1,layer_4_p_fibers_d[i4], 4));
     } 
    SFT_Z4_top.push_back(SFT_Z4_top_col);
    SFT_Z4_bottom.push_back(SFT_Z4_bottom_col);
    DZ4_top.push_back(DZ4_top_col);
    DZ4_bottom.push_back(DZ4_bottom_col);
    SFT_Z4_top_col.clear();
    SFT_Z4_bottom_col.clear();
    DZ4_top_col.clear();
    DZ4_bottom_col.clear();
  }

  /*
  cout << endl;
  cout << "LAYER 1 " << endl; 
  for(unsigned int i1=0; i1<SFT_Z1_top.size(); i1++){
    for(unsigned int j1=0; j1<SFT_Z1_top[i1].size(); j1++){
        cout << "****************  " << i1 << "   " << j1 << "   " << SFT_Z1_top[i1][j1] << "   " << SFT_Z1_bottom[i1][j1] << "   " << C2_Z_top << "   " << C2_Z_bottom << "   " << DZ1_top[i1][j1] << "   " << DZ1_bottom[i1][j1] << endl;
     }
  }
 
  cout << endl;
  cout << "LAYER 2 " << endl; 
  for(unsigned int i2=0; i2<SFT_Z2_top.size(); i2++){
    for(unsigned int j2=0; j2<SFT_Z2_top[i2].size(); j2++){
        cout << "****************  " << i2 << "   " << j2 << "   " << SFT_Z2_top[i2][j2] << "   " << SFT_Z2_bottom[i2][j2] << "   " << C2_Z_top << "   " << C2_Z_bottom << "   "  << DZ2_top[i2][j2] << "   " << DZ2_bottom[i2][j2] <<  endl;
    }
  }

  cout << endl;
  cout << "LAYER 3 " << endl; 
  for(unsigned int i3=0; i3<SFT_Z3_top.size(); i3++){
    for(unsigned int j3=0; j3<SFT_Z3_top[i3].size(); j3++){
        cout << "****************  " << i3 << "   " << j3 << "   " << SFT_Z3_top[i3][j3] << "   " << SFT_Z3_bottom[i3][j3] << "   " << C2_Z_top << "   " << C2_Z_bottom << "   "  << DZ3_top[i3][j3] << "   " << DZ3_bottom[i3][j3] <<   endl;
    }
  }

  cout << endl;
  cout << "LAYER 4 " << endl; 
  for(unsigned int i4=0; i4<SFT_Z4_top.size(); i4++){
    for(unsigned int j4=0; j4<SFT_Z4_top[i4].size(); j4++){
        cout << "****************  " << i4 << "   " << j4 << "   " << SFT_Z4_top[i4][j4] << "   " << SFT_Z4_bottom[i4][j4] << "   " << C2_Z_top << "   " << C2_Z_bottom << "   "  << DZ4_top[i4][j4] << "   " << DZ4_bottom[i4][j4] <<   endl;
    }
  }
  cout << endl;
  cout << endl;
  */

  for(unsigned int i1=0; i1<SFT_Z1_top.size(); i1++){
    for(unsigned int j1=0; j1<SFT_Z1_top[i1].size(); j1++){
      if(abs(DZ1_top[i1][j1]) < lowest_top){
        lowest_top = abs(DZ1_top[i1][j1]);
        n_lowest_top = j1+1;
        z_lowest_top = SFT_Z1_top[i1][j1];
      }
      if(abs(DZ1_bottom[i1][j1]) < lowest_bottom){
        lowest_bottom = abs(DZ1_bottom[i1][j1]);
        n_lowest_bottom = j1+1;
        z_lowest_bottom = SFT_Z1_bottom[i1][j1];
      }
    }
    n1_lowest_top.push_back(n_lowest_top);
    n1_lowest_bottom.push_back(n_lowest_bottom);
    Z1_lowest_top.push_back(z_lowest_top);
    Z1_lowest_bottom.push_back(z_lowest_bottom);
    DZ1_lowest_top.push_back(lowest_top);
    DZ1_lowest_bottom.push_back(lowest_bottom);
    lowest_top = 999.;
    lowest_bottom = 999.;
    z_lowest_top = 999.;
    z_lowest_bottom = 999.;
    n_lowest_top = 99;
    n_lowest_bottom = 99;
  }

  for(unsigned int i2=0; i2<SFT_Z2_top.size(); i2++){
    for(unsigned int j2=0; j2<SFT_Z2_top[i2].size(); j2++){
      if(abs(DZ2_top[i2][j2]) < lowest_top){
        lowest_top = abs(DZ2_top[i2][j2]);
        n_lowest_top = j2+1;
        z_lowest_top = SFT_Z2_top[i2][j2];
     }
      if(abs(DZ2_bottom[i2][j2]) < lowest_bottom){
        lowest_bottom = abs(DZ2_bottom[i2][j2]);
        n_lowest_bottom = j2+1;
        z_lowest_bottom = SFT_Z2_bottom[i2][j2];
      }
   }
    n2_lowest_top.push_back(n_lowest_top);
    n2_lowest_bottom.push_back(n_lowest_bottom);
    Z2_lowest_top.push_back(z_lowest_top);
    Z2_lowest_bottom.push_back(z_lowest_bottom);
    DZ2_lowest_top.push_back(lowest_top);
    DZ2_lowest_bottom.push_back(lowest_bottom);
    lowest_top = 999.;
    lowest_bottom = 999.;
    z_lowest_top = 999.;
    z_lowest_bottom = 999.;
    n_lowest_top = 99;
    n_lowest_bottom = 99;
  }

  for(unsigned int i3=0; i3<SFT_Z3_top.size(); i3++){
    for(unsigned int j3=0; j3<SFT_Z3_top[i3].size(); j3++){
      if(abs(DZ3_top[i3][j3]) < lowest_top){
        lowest_top = abs(DZ3_top[i3][j3]);
        n_lowest_top = j3+1;
        z_lowest_top = SFT_Z3_top[i3][j3];
      }
      if(abs(DZ3_bottom[i3][j3]) < lowest_bottom){
        lowest_bottom = abs(DZ3_bottom[i3][j3]);
        n_lowest_bottom = j3+1;
        z_lowest_bottom = SFT_Z3_bottom[i3][j3];
      }
   }
    n3_lowest_top.push_back(n_lowest_top);
    n3_lowest_bottom.push_back(n_lowest_bottom);
    Z3_lowest_top.push_back(z_lowest_top);
    Z3_lowest_bottom.push_back(z_lowest_bottom);
    DZ3_lowest_top.push_back(lowest_top);
    DZ3_lowest_bottom.push_back(lowest_bottom);
    lowest_top = 999.;
    lowest_bottom = 999.;
    z_lowest_top = 999.;
    z_lowest_bottom = 999.;
    n_lowest_top = 99;
    n_lowest_bottom = 99;
  }

  for(unsigned int i4=0; i4<SFT_Z4_top.size(); i4++){
    for(unsigned int j4=0; j4<SFT_Z4_top[i4].size(); j4++){
      if(abs(DZ4_top[i4][j4]) < lowest_top){
        lowest_top = abs(DZ4_top[i4][j4]);
        n_lowest_top = j4+1;
        z_lowest_top = SFT_Z4_top[i4][j4];
      }  
      if(abs(DZ4_bottom[i4][j4]) < lowest_bottom){
        lowest_bottom = abs(DZ4_bottom[i4][j4]);
        n_lowest_bottom = j4+1;
        z_lowest_bottom = SFT_Z4_bottom[i4][j4];
      }
   }
    n4_lowest_top.push_back(n_lowest_top);
    n4_lowest_bottom.push_back(n_lowest_bottom);
    Z4_lowest_top.push_back(z_lowest_top);
    Z4_lowest_bottom.push_back(z_lowest_bottom);
    DZ4_lowest_top.push_back(lowest_top);
    DZ4_lowest_bottom.push_back(lowest_bottom);
    lowest_top = 999.;
    lowest_bottom = 999.;
    z_lowest_top = 999.;
    z_lowest_bottom = 999.;
    n_lowest_top = 99;
    n_lowest_bottom = 99;
  }

  /*
  cout << "LAYER 1:" << endl;
  for(unsigned int i1=0; i1<DZ1_lowest_top.size(); i1++){
    cout << n1_lowest_top[i1] << " " << layer_1_p_fibers_d[i1] << "   " << Z1_lowest_top[i1] << "  " << C2_Z_top << "  " << DZ1_lowest_top[i1] << "   ";
  }
  cout << endl;

  for(unsigned int i1=0; i1<DZ1_lowest_bottom.size(); i1++){
    cout << n1_lowest_bottom[i1] << " " << layer_1_p_fibers_d[i1] << "   " << Z1_lowest_bottom[i1] << "  " << C2_Z_bottom << "  " << DZ1_lowest_bottom[i1] << "   ";
  }

  cout << endl;
  cout << endl;

  cout << "LAYER 2:" << endl;
  for(unsigned int i2=0; i2<DZ2_lowest_top.size(); i2++){
    cout << DZ2_lowest_top[i2] << "   ";
  }
  cout << endl;

  for(unsigned int i2=0; i2<DZ2_lowest_bottom.size(); i2++){
    cout << DZ2_lowest_bottom[i2] << "   ";
  }
  cout << endl;
  cout << endl;

  cout << "LAYER 3:" << endl;
  for(unsigned int i3=0; i3<DZ3_lowest_top.size(); i3++){
    cout << DZ3_lowest_top[i3] << "   ";
  }
  cout << endl;

  for(unsigned int i3=0; i3<DZ3_lowest_bottom.size(); i3++){
    cout << DZ3_lowest_bottom[i3] << "   ";
  }
  cout << endl;
  cout << endl;

  cout << "LAYER 4:" << endl;
  for(unsigned int i4=0; i4<DZ4_lowest_top.size(); i4++){
    cout << DZ4_lowest_top[i4] << "   ";
  }
  cout << endl;

  for(unsigned int i4=0; i4<DZ4_lowest_bottom.size(); i4++){
    cout << DZ4_lowest_bottom[i4] << "   ";
  }
  cout << endl;
  cout << endl;
  */
  
  cout << endl;
  cout << " |===========|==========================================|==========================================| " << endl;
  cout << " |   Run#    |                   TOP                    |                  BOTTOM                  |" << endl;
  printf(" |   %4d    |                                          |                                          |\n",Run_Number);
  //  cout << " |   4560    |                                          |                                          | " << endl;
  cout << " |===========|==========================================|==========================================| " << endl;
  printf(" |           |            phi = %8.3f deg.           |             phi = %8.3f deg.          |\n",phi_top,phi_bottom);
  printf(" |   Evt#    |           C2_Z = %8.3f   mm           |            C2_Z = %8.3f   mm          |\n",C2_Z_top,C2_Z_bottom);
  printf(" |   %4d    |==========================================|==========================================|\n",evt);
  //cout << " |           |==========================================|==========================================| " << endl;
  //  cout << " |           |       p       n       Z         DZ       |       p       n       Z         DZ       |" << endl;
  printf(" |           |       p       n       Z         DZ       |       p       n       Z         DZ       |\n");
 cout << " |===========|==========================================|==========================================| " << endl;

  if(layer_1_p_fibers_d.size()==0){
    printf(" |  LAYER 1  |                -----------               |               -----------                |\n");
  }
  else{
    printf(" |  LAYER 1  |     %5.2f    %2d   %8.3f  %8.3f     |     %5.2f    %2d   %8.3f  %8.3f     |\n", 
    layer_1_p_fibers_d[0], n1_lowest_top[0], Z1_lowest_top[0], DZ1_lowest_top[0], layer_1_p_fibers_d[0], n1_lowest_bottom[0], Z1_lowest_bottom[0], DZ1_lowest_bottom[0]);
    for(unsigned int i1=0; i1<DZ1_lowest_top.size()-1; i1++){
    printf(" |           |     %5.2f    %2d   %8.3f  %8.3f     |     %5.2f    %2d   %8.3f  %8.3f     |\n", 
    layer_1_p_fibers_d[i1+1], n1_lowest_top[i1+1], Z1_lowest_top[i1+1], DZ1_lowest_top[i1+1], layer_1_p_fibers_d[i1+1], n1_lowest_bottom[i1+1], Z1_lowest_bottom[i1+1], DZ1_lowest_bottom[i1+1]);
    }
  }

  cout << " |===========|==========================================|==========================================| " << endl;

  if(layer_2_p_fibers_d.size()==0){
    printf(" |  LAYER 2  |               -----------                |               -----------                |\n");
  }
  else{
    printf(" |  LAYER 2  |     %5.2f    %2d   %8.3f  %8.3f     |     %5.2f    %2d   %8.3f  %8.3f     |\n", 
    layer_2_p_fibers_d[0], n2_lowest_top[0], Z2_lowest_top[0], DZ2_lowest_top[0], layer_2_p_fibers_d[0], n2_lowest_bottom[0], Z2_lowest_bottom[0], DZ2_lowest_bottom[0]);
    for(unsigned int i2=0; i2<DZ2_lowest_top.size()-1; i2++){
    printf(" |           |     %5.2f    %2d   %8.3f  %8.3f     |     %5.2f    %2d   %8.3f  %8.3f     |\n", 
    layer_2_p_fibers_d[i2+1], n2_lowest_top[i2+1], Z2_lowest_top[i2+1], DZ2_lowest_top[i2+1], layer_2_p_fibers_d[i2+1], n2_lowest_bottom[i2+1], Z2_lowest_bottom[i2+1], DZ2_lowest_bottom[i2+1]);
    }
  }

  cout << " |===========|==========================================|==========================================| " << endl;

  if(layer_3_p_fibers_d.size()==0){
    printf(" |  LAYER 3  |                -----------                |                -----------                |\n");
  }
  else{
    printf(" |  LAYER 3  |     %5.2f    %2d   %8.3f  %8.3f     |     %5.2f    %2d   %8.3f  %8.3f     |\n", 
    layer_3_p_fibers_d[0], n3_lowest_top[0], Z3_lowest_top[0], DZ3_lowest_top[0], layer_3_p_fibers_d[0], n3_lowest_bottom[0], Z3_lowest_bottom[0], DZ3_lowest_bottom[0]);
    for(unsigned int i3=0; i3<DZ3_lowest_top.size()-1; i3++){
    printf(" |           |     %5.2f    %2d   %8.3f  %8.3f     |     %5.2f    %2d   %8.3f  %8.3f     |\n", 
    layer_3_p_fibers_d[i3+1], n3_lowest_top[i3+1], Z3_lowest_top[i3+1], DZ3_lowest_top[i3+1], layer_3_p_fibers_d[i3+1], n3_lowest_bottom[i3+1], Z3_lowest_bottom[i3+1], DZ3_lowest_bottom[i3+1]);
    }
  }

  cout << " |===========|==========================================|==========================================| " << endl;

  if(layer_4_p_fibers_d.size()==0){
    printf(" |  LAYER 4  |                -----------                |                -----------                |\n");
  }
  else{
    printf(" |  LAYER 4  |     %5.2f    %2d   %8.3f  %8.3f     |     %5.2f    %2d   %8.3f  %8.3f     |\n", 
    layer_4_p_fibers_d[0], n4_lowest_top[0], Z4_lowest_top[0], DZ4_lowest_top[0], layer_4_p_fibers_d[0], n4_lowest_bottom[0], Z4_lowest_bottom[0], DZ4_lowest_bottom[0]);
    for(unsigned int i4=0; i4<DZ4_lowest_top.size()-1; i4++){
    printf(" |           |     %5.2f    %2d   %8.3f  %8.3f     |     %5.2f    %2d   %8.3f  %8.3f     |\n", 
    layer_4_p_fibers_d[i4+1], n4_lowest_top[i4+1], Z4_lowest_top[i4+1], DZ4_lowest_top[i4+1], layer_4_p_fibers_d[i4+1], n4_lowest_bottom[i4+1], Z4_lowest_bottom[i4+1], DZ4_lowest_bottom[i4+1]);
    }
  }

  cout << " |===========|==========================================|==========================================| " << endl;
  cout << endl;
  cout << endl;

  vector<double> Z1_top;  vector<double> Z1_bottom;
  vector<double> Z2_top;  vector<double> Z2_bottom;
  vector<double> Z3_top;  vector<double> Z3_bottom;
  vector<double> Z4_top;  vector<double> Z4_bottom;

  vector<double> p1_top;  vector<double> p1_bottom;
  vector<double> p2_top;  vector<double> p2_bottom;
  vector<double> p3_top;  vector<double> p3_bottom;
  vector<double> p4_top;  vector<double> p4_bottom;

  vector<int> n1_top;  vector<int> n1_bottom;
  vector<int> n2_top;  vector<int> n2_bottom;
  vector<int> n3_top;  vector<int> n3_bottom;
  vector<int> n4_top;  vector<int> n4_bottom;

  vector<double> D_Z1_top;  vector<double> D_Z1_bottom;
  vector<double> D_Z2_top;  vector<double> D_Z2_bottom;
  vector<double> D_Z3_top;  vector<double> D_Z3_bottom;
  vector<double> D_Z4_top;  vector<double> D_Z4_bottom;


  // LAYER 1
  for(unsigned int i1=0; i1<layer_1_p_fibers_d.size(); i1++){
    if(DZ1_lowest_top[i1]<= DZ1_lowest_bottom[i1]){
      n1_top.push_back(n1_lowest_top[i1]);
      p1_top.push_back(layer_1_p_fibers_d[i1]);
      Z1_top.push_back(Z1_lowest_top[i1]);
      D_Z1_top.push_back(DZ1_lowest_top[i1]);
    }  
    else{
      n1_bottom.push_back(n1_lowest_bottom[i1]);
      p1_bottom.push_back(layer_1_p_fibers_d[i1]);
      Z1_bottom.push_back(Z1_lowest_bottom[i1]);
      D_Z1_bottom.push_back(DZ1_lowest_bottom[i1]);
    }
  }

  // LAYER 2
  for(unsigned int i2=0; i2<layer_2_p_fibers_d.size(); i2++){
    if(DZ2_lowest_top[i2]<= DZ2_lowest_bottom[i2]){
      n2_top.push_back(n2_lowest_top[i2]);
      p2_top.push_back(layer_2_p_fibers_d[i2]);
      Z2_top.push_back(Z2_lowest_top[i2]);
      D_Z2_top.push_back(DZ2_lowest_top[i2]);
    }  
    else{
      n2_bottom.push_back(n2_lowest_bottom[i2]);
      p2_bottom.push_back(layer_2_p_fibers_d[i2]);
      Z2_bottom.push_back(Z2_lowest_bottom[i2]);
      D_Z2_bottom.push_back(DZ2_lowest_bottom[i2]);
    }
  } 

  // LAYER 3
  for(unsigned int i3=0; i3<layer_3_p_fibers_d.size(); i3++){
    if(DZ3_lowest_top[i3]<= DZ3_lowest_bottom[i3]){
      n3_top.push_back(n3_lowest_top[i3]);
      p3_top.push_back(layer_3_p_fibers_d[i3]);
      Z3_top.push_back(Z3_lowest_top[i3]);
      D_Z3_top.push_back(DZ3_lowest_top[i3]);
    }  
    else{
      n3_bottom.push_back(n3_lowest_bottom[i3]);
      p3_bottom.push_back(layer_3_p_fibers_d[i3]);
      Z3_bottom.push_back(Z3_lowest_bottom[i3]);
      D_Z3_bottom.push_back(DZ3_lowest_bottom[i3]);
    }
  } 

  // LAYER 4
  for(unsigned int i4=0; i4<layer_4_p_fibers_d.size(); i4++){
    if(DZ4_lowest_top[i4]<= DZ4_lowest_bottom[i4]){
      n4_top.push_back(n4_lowest_top[i4]);
      p4_top.push_back(layer_4_p_fibers_d[i4]);
      Z4_top.push_back(Z4_lowest_top[i4]);
      D_Z4_top.push_back(DZ4_lowest_top[i4]);
    }  
    else{
      n4_bottom.push_back(n4_lowest_bottom[i4]);
      p4_bottom.push_back(layer_4_p_fibers_d[i4]);
      Z4_bottom.push_back(Z4_lowest_bottom[i4]);
      D_Z4_bottom.push_back(DZ4_lowest_bottom[i4]);
    }
  } 

  /*
  for(unsigned int i1=0; i1<Z1_top.size(); i1++){
    cout << "Z1_top(" << n1_top[i1] << "," << p1_top[i1] << ")= " << Z1_top[i1] << ",  DZ1_top= " << D_Z1_top[i1] << endl;
  }
  for(unsigned int j1=0; j1<Z1_bottom.size(); j1++){
    cout << "Z1_bottom(" << n1_bottom[j1] << "," << p1_bottom[j1] << ")= " << Z1_bottom[j1] << ",  DZ1_bottom= " << D_Z1_bottom[j1] << endl;
  }
  cout << endl;

  for(unsigned int i2=0; i2<Z2_top.size(); i2++){
    cout << "Z2_top(" << n2_top[i2] << "," << p2_top[i2] << ")= " << Z2_top[i2] << ",  DZ2_top= " << D_Z2_top[i2] << endl;
  }
  for(unsigned int j2=0; j2<Z2_bottom.size(); j2++){
    cout << "Z2_bottom(" << n2_bottom[j2] << "," << p2_bottom[j2] << ")= " << Z2_bottom[j2] << ",  DZ2_bottom= " << D_Z2_bottom[j2] << endl;
  }
  cout << endl;

  for(unsigned int i3=0; i3<Z3_top.size(); i3++){
    cout << "Z3_top(" << n3_top[i3] << "," << p3_top[i3] << ")= " << Z3_top[i3] << ",  DZ3_top= " << D_Z3_top[i3] << endl;
  }
  for(unsigned int j3=0; j3<Z3_bottom.size(); j3++){
    cout << "Z3_bottom(" << n3_bottom[j3] << "," << p3_bottom[j3] << ")= " << Z3_bottom[j3] << ",  DZ3_bottom= " << D_Z3_bottom[j3] << endl;
  }
  cout << endl;

  for(unsigned int i4=0; i4<Z4_top.size(); i4++){
    cout << "Z4_top(" << n4_top[i4] << "," << p4_top[i4] << ")= " << Z4_top[i4] << ",  DZ4_top= " << D_Z4_top[i4] << endl;
  }
  for(unsigned int j4=0; j4<Z4_bottom.size(); j4++){
    cout << "Z4_bottom(" << n4_bottom[j4] << "," << p4_bottom[j4] << ")= " << Z4_bottom[j4] << ",  DZ4_bottom= " << D_Z4_bottom[j4] << endl;
  }
  cout << endl;
  */

  vector<double> Z_top_final;   vector<double> Z_bottom_final;
  double DZ1_top_min = 999.;   double DZ1_bottom_min = 999.;  
  double DZ2_top_min = 999.;   double DZ2_bottom_min = 999.;  
  double DZ3_top_min = 999.;   double DZ3_bottom_min = 999.;  
  double DZ4_top_min = 999.;   double DZ4_bottom_min = 999.;  

  int i1_top_min = 99;  int i1_bot_min = 99;
  int i2_top_min = 99;  int i2_bot_min = 99;
  int i3_top_min = 99;  int i3_bot_min = 99;
  int i4_top_min = 99;  int i4_bot_min = 99;

  
  // LAYER 1
  for(unsigned int i1=0; i1<Z1_top.size(); i1++){
    if(D_Z1_top[i1] < DZ1_top_min){
      DZ1_top_min = D_Z1_top[i1];
      i1_top_min = i1;
    }
  }
  for(unsigned int j1=0; j1<Z1_bottom.size(); j1++){
    if(D_Z1_bottom[j1] < DZ1_bottom_min){
      DZ1_bottom_min = D_Z1_bottom[j1];
      i1_bot_min = j1;
    }
  }
  if(Z1_top.size()!=0) Z_top_final.push_back(Z1_top[i1_top_min]); 
  if(Z1_bottom.size()!=0) Z_bottom_final.push_back(Z1_bottom[i1_bot_min]); 

  
  //LAYER 2
  for(unsigned int i2=0; i2<Z2_top.size(); i2++){
    if(D_Z2_top[i2] < DZ2_top_min){
      DZ2_top_min = D_Z2_top[i2];
      i2_top_min = i2;
    }
  }
  for(unsigned int j2=0; j2<Z2_bottom.size(); j2++){
    if(D_Z2_bottom[j2] < DZ2_bottom_min){
      DZ2_bottom_min = D_Z2_bottom[j2];
      i2_bot_min = j2;
    }
  }
  if(Z2_top.size()!=0) Z_top_final.push_back(Z2_top[i2_top_min]); 
  if(Z2_bottom.size()!=0) Z_bottom_final.push_back(Z2_bottom[i2_bot_min]); 

  
  
  //LAYER 3
  for(unsigned int i3=0; i3<Z3_top.size(); i3++){
    if(D_Z3_top[i3] < DZ3_top_min){
      DZ3_top_min = D_Z3_top[i3];
      i3_top_min = i3;
    }
  }
  for(unsigned int j3=0; j3<Z3_bottom.size(); j3++){
    if(D_Z3_bottom[j3] < DZ3_bottom_min){
      DZ3_bottom_min = D_Z3_bottom[j3];
      i3_bot_min = j3;
    }
  }
  if(Z3_top.size()!=0) Z_top_final.push_back(Z3_top[i3_top_min]);
  if(Z3_bottom.size()!=0) Z_bottom_final.push_back(Z3_bottom[i3_bot_min]);
 
  //LAYER 4
  for(unsigned int i4=0; i4<Z4_top.size(); i4++){
    if(D_Z4_top[i4] < DZ4_top_min){
      DZ4_top_min = D_Z4_top[i4];
      i4_top_min = i4;
    }
  }
  for(unsigned int j4=0; j4<Z4_bottom.size(); j4++){
    if(D_Z4_bottom[j4] < DZ4_bottom_min){
      DZ4_bottom_min = D_Z4_bottom[j4];
      i4_bot_min = j4;
    }
  }
  if(Z4_top.size()!=0) Z_top_final.push_back(Z4_top[i4_top_min]);
  if(Z4_bottom.size()!=0) Z_bottom_final.push_back(Z4_bottom[i4_bot_min]);


  //cout << "TEST : " << i1_top_min << "  " << i1_bot_min << endl;
  //cout << "TEST : " << Z1_top[i1_top_min] << "  " << Z1_bottom[i1_bot_min] << endl;
  //cout << "TEST : " << i1_bot_min << "  " << Z1_bottom[0] << "  " << Z1_bottom[i1_bot_min] << endl;
  //Z_top_final.push_back(Z3_top[i3_top_min]);
  //Z_bottom_final.push_back(Z1_bottom[i1_bot_min]);



  //Z_top_final.push_back(Z3_top[i_top_min]);
  //Z_bottom_final.push_back(Z3_bottom[i_bot_min]);

  //cout << "TEST:  " << i_top_min << "   " << Z1_top[i_top_min] << "   " << Z_top_final.size() << "   " << Z_top_final[0] << endl;
  if(Z_top_final.size()>0){
    cout << " Z_top :    ";
    for(unsigned int i=0; i<Z_top_final.size(); i++) cout << Z_top_final[i] << "  ";
    cout << endl;
  }
  if(Z_bottom_final.size()>0){
    cout << " Z_bottom : ";
    for(unsigned int j=0; j<Z_bottom_final.size(); j++) cout << Z_bottom_final[j] << "  ";
    cout << endl;
  }
  
  cout << endl;
  cout << endl;

  double Z_top_FINAL = 999.;    double Z_bottom_FINAL = 999.;
  double DZ_top_FINAL = 999.;   double DZ_bottom_FINAL = 999.;
  double sum_top = 0.;          double sum_bottom = 0.;

  if(Z_top_final.size()!=0){
    for(unsigned int i=0; i<Z_top_final.size(); i++){
      sum_top += Z_top_final[i];
    }
    Z_top_FINAL =  sum_top/Z_top_final.size();
    DZ_top_FINAL = Z_top_FINAL - C2_Z_top;
  }

  if(Z_bottom_final.size()!=0){
    for(unsigned int i=0; i<Z_bottom_final.size(); i++){
      sum_bottom += Z_bottom_final[i];
    }
    Z_bottom_FINAL =  sum_bottom/Z_bottom_final.size(); 
    DZ_bottom_FINAL = Z_bottom_FINAL - C2_Z_bottom;
  }

  //if(Z_top_FINAL>500) cout << "Z_TOP_FINAL =   -------  " << endl;
  cout << " ===============================================" << endl;
  cout << endl;
  //if(Z_top_FINAL<500) cout << " TOP:      Z = " << Z_top_FINAL << "  ;  Delta_Z = " << DZ_top_FINAL << endl;
  if(Z_top_FINAL<500) printf("   TOP:     Z = %8.3f  ;  Delta_Z = %7.3f\n",Z_top_FINAL,DZ_top_FINAL);

  //if(Z_bottom_FINAL>500) cout << "Z_BOTTOM_FINAL =   -------  " << endl;
  //if(Z_bottom_FINAL<500) cout << " BOTTOM:   Z = " << Z_bottom_FINAL << "  ;  Delta_Z = " << DZ_bottom_FINAL << endl;
  if(Z_bottom_FINAL<500) printf("   BOTTOM:  Z = %8.3f  ;  Delta_Z = %7.3f\n",Z_bottom_FINAL,DZ_bottom_FINAL);
  cout << endl;
  cout << " ===============================================" << endl;

  cout << endl;
  cout << endl;

 } // End of Ana3

//vector<double> Batch_Ana3(int Run_Number, int evt, double ADC_High_SFT_corr[128], int has_TDC_SFT_hit[128], int SFT_channel_to_fiber[128], double phi, int gap_to_fit_left, int gap_to_fit_right, double C2_Z_top, double C2_Z_bottom){
vector<double> Batch_Ana3(double ADC_High_SFT_corr[128], int has_TDC_SFT_hit[128], int SFT_channel_to_fiber[128], double phi, int gap_to_fit_left, int gap_to_fit_right, double C2_Z_top, double C2_Z_bottom){
  vector<double> vec_output;

  //cout << "###############################    C2_Z-SFT_Z ANALYSIS    ############################### " << endl;
  //cout << endl;
  //cout << "Event " << evt << endl; 
  //cout << "DEBUG : " <<  << endl;
 //  cout << "TEST :  " << gap_to_fit_left << "   " << gap_to_fit_right << endl;
 //  cout << "DEBUG :  " << evt << endl;
  //cout << endl;

  //cout << "INPUTS " << endl;
  //for(int i=0;i<128;i++){
  //  cout << i << "   " << ADC_High_SFT_corr[i] << endl;
  //  cout << i << "   " << SFT_channel_to_fiber[i] << endl;
  //}

  //const int sft_hit_count_max = 3;
  //const int NUM_ITERATIONS = 100;

  //float R_TARGET = 29.0;
  //float R_TOF1 = 47.1;
  //float R_SFT_L1 = 40.0;
  //float R_C2 = 629.4;

  int L1_sft_count = 0;
  int L2_sft_count = 0;
  int L3_sft_count = 0;
  int L4_sft_count = 0;

  double phi_top = 0.;
  double phi_bottom = 0.;

  vector<int> layer_1_p_fibers;
  vector<int> layer_2_p_fibers;
  vector<int> layer_3_p_fibers;
  vector<int> layer_4_p_fibers;

  vector<double> layer_1_adc;
  vector<double> layer_2_adc;
  vector<double> layer_3_adc;
  vector<double> layer_4_adc;

  vector<double> layer_1_p_fibers_d;
  vector<double> layer_2_p_fibers_d;
  vector<double> layer_3_p_fibers_d;
  vector<double> layer_4_p_fibers_d;
  
  vector< vector<double> > SFT_Z1_top;  vector<double> SFT_Z1_top_col;  
  vector< vector<double> > SFT_Z2_top;  vector<double> SFT_Z2_top_col;  
  vector< vector<double> > SFT_Z3_top;  vector<double> SFT_Z3_top_col;  
  vector< vector<double> > SFT_Z4_top;  vector<double> SFT_Z4_top_col;  

  vector< vector<double> > SFT_Z1_bottom;  vector<double> SFT_Z1_bottom_col;  
  vector< vector<double> > SFT_Z2_bottom;  vector<double> SFT_Z2_bottom_col;  
  vector< vector<double> > SFT_Z3_bottom;  vector<double> SFT_Z3_bottom_col;  
  vector< vector<double> > SFT_Z4_bottom;  vector<double> SFT_Z4_bottom_col;  

  vector< vector<double> > DZ1_top;  vector<double> DZ1_top_col;  
  vector< vector<double> > DZ2_top;  vector<double> DZ2_top_col;  
  vector< vector<double> > DZ3_top;  vector<double> DZ3_top_col;  
  vector< vector<double> > DZ4_top;  vector<double> DZ4_top_col;  

  vector< vector<double> > DZ1_bottom;  vector<double> DZ1_bottom_col;  
  vector< vector<double> > DZ2_bottom;  vector<double> DZ2_bottom_col;  
  vector< vector<double> > DZ3_bottom;  vector<double> DZ3_bottom_col;  
  vector< vector<double> > DZ4_bottom;  vector<double> DZ4_bottom_col;  

  vector<double> DZ1_lowest_top;  vector<double> DZ1_lowest_bottom;
  vector<double> DZ2_lowest_top;  vector<double> DZ2_lowest_bottom;
  vector<double> DZ3_lowest_top;  vector<double> DZ3_lowest_bottom;
  vector<double> DZ4_lowest_top;  vector<double> DZ4_lowest_bottom;

  vector<int> n1_lowest_top;      vector<int> n1_lowest_bottom;
  vector<int> n2_lowest_top;      vector<int> n2_lowest_bottom;
  vector<int> n3_lowest_top;      vector<int> n3_lowest_bottom;
  vector<int> n4_lowest_top;      vector<int> n4_lowest_bottom;

  vector<int> p1_lowest_top;      vector<int> p1_lowest_bottom;
  vector<int> p2_lowest_top;      vector<int> p2_lowest_bottom;
  vector<int> p3_lowest_top;      vector<int> p3_lowest_bottom;
  vector<int> p4_lowest_top;      vector<int> p4_lowest_bottom;

  vector<double> Z1_lowest_top;      vector<double> Z1_lowest_bottom;
  vector<double> Z2_lowest_top;      vector<double> Z2_lowest_bottom;
  vector<double> Z3_lowest_top;      vector<double> Z3_lowest_bottom;
  vector<double> Z4_lowest_top;      vector<double> Z4_lowest_bottom;


  double lowest_top = 999.;     double lowest_bottom = 999.;
  int n_lowest_top = 99;        int n_lowest_bottom = 99;
  double z_lowest_top = 999.;   double z_lowest_bottom = 999.;
  //int n_lowest_top[4]={0};    int n_lowest_bottom[4] = {0};

  bool L1_sft_flag = false;
  bool L2_sft_flag = false;
  bool L3_sft_flag = false;
  bool L4_sft_flag = false;

  double Z_top_FINAL = 999.99;    double Z_bottom_FINAL = 999.99;
  double DZ_top_FINAL = 999.99;   double DZ_bottom_FINAL = 999.99;
  double sum_top = 0.;             double sum_bottom = 0.;


 
  //cout << "DEBUG :   Gap_to_fit_LEFT : " <<  gap_to_fit_left << "  gap_to_fit_RIGHT : " << gap_to_fit_right << endl;
  //cout << endl;

  //  for(int i=Run_Number; i<Run_Number+1; i++){
  //    if(gap_to_fit_left!=12){
  //    cout << "Nope !" << endl;
  //    break;
      
  //  }

    //if(gap_to_fit_left==11 && gap_to_fit_right!=5) break; 
    //if(gap_to_fit_left==10 && gap_to_fit_right!=4) break;   
    //if(gap_to_fit_left==8 && gap_to_fit_right!=2) break;   
    //if(gap_to_fit_left==7 && gap_to_fit_right!=1) break;   

    if((gap_to_fit_left>=7 && gap_to_fit_left<=12) && (gap_to_fit_right>=1 && gap_to_fit_right<=6)){
      if(phi > 180){
        phi_top = phi-180;
        phi_bottom = phi; 
      }
      else{
        phi_top = phi;
        phi_bottom = phi+180;
      }
    }
    /*
    if(gap_to_fit_left==7 && gap_to_fit_right==2 ){
      if(phi > 270){
        phi_top = phi-180;
        phi_bottom = phi;
      }
      else{
        phi_top = phi;
        phi_bottom = phi+180;
      }   
    }
    if(gap_to_fit_left==7 && gap_to_fit_right==1){
      if(phi > 270){
        phi_top = phi-180;
        phi_bottom = phi;
      }
      else{
        phi_top = phi;
        phi_bottom = phi+180;
      }   
    }
    */ 


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

  //bool bad_sft = false;
  //bool *bad_sft_flag = &bad_sft;
  
  //layer_1_p_fibers_d.clear();
  //layer_2_p_fibers_d.clear();
  //layer_3_p_fibers_d.clear();
  //layer_4_p_fibers_d.clear();

  //layer_1_p_fibers_d = Batch_Layer_averaging(layer_1_p_fibers, bad_sft_flag, 1);
  //layer_2_p_fibers_d = Batch_Layer_averaging(layer_2_p_fibers, bad_sft_flag, 2);
  //layer_3_p_fibers_d = Batch_Layer_averaging(layer_3_p_fibers, bad_sft_flag, 3);
  //layer_4_p_fibers_d = Batch_Layer_averaging(layer_4_p_fibers, bad_sft_flag, 4);

  for(unsigned int i1=0; i1<layer_1_p_fibers.size(); i1++) layer_1_p_fibers_d.push_back(double(layer_1_p_fibers[i1]));
  for(unsigned int i2=0; i2<layer_2_p_fibers.size(); i2++) layer_2_p_fibers_d.push_back(double(layer_2_p_fibers[i2]));
  for(unsigned int i3=0; i3<layer_3_p_fibers.size(); i3++) layer_3_p_fibers_d.push_back(double(layer_3_p_fibers[i3]));
  for(unsigned int i4=0; i4<layer_4_p_fibers.size(); i4++) layer_4_p_fibers_d.push_back(double(layer_4_p_fibers[i4]));


  //for(unsigned int i2=0; i2<layer_2_p_fibers.size(); i2++){
  //  layer_2_p_fibers_d[i2].push_back(double(layer_2_p_fibers[i2]);
  //}
  //for(unsigned int i3=0; i3<layer_3_p_fibers.size(); i3++){
  //  layer_3_p_fibers_d[i3].push_back(double(layer_3_p_fibers[i3]);
  //}
  //for(unsigned int i4=0; i4<layer_4_p_fibers.size(); i4++){
  //  layer_4_p_fibers_d[i4].push_back(double(layer_4_p_fibers[i4]);
  //}

  /*
  cout << endl;
  cout << endl;
  cout << "####### Ana 3  ########" << endl;
  cout << "phi (TOP) = " << phi_top << " (deg.)  ,  phi (BOTTOM) = " << phi_bottom << " (deg.)" << endl;
  cout << "C2_Z (TOP) = " << C2_Z_top << "  " << "C2_Z (BOTTOM) = " << C2_Z_bottom << endl;
  cout << "p1 = ";
  for (unsigned int i1=0; i1<layer_1_p_fibers_d.size(); i1++) cout << layer_1_p_fibers_d[i1] << "  ";
  cout << endl; 
  cout << "p2 = "; 
  for (unsigned int i2=0; i2<layer_2_p_fibers_d.size(); i2++) cout << layer_2_p_fibers_d[i2] << "  ";
  cout << endl;
  cout << "p3 = ";
  for (unsigned int i3=0; i3<layer_3_p_fibers_d.size(); i3++) cout << layer_3_p_fibers_d[i3] << "  ";
  cout << endl;
  cout << "p4 = ";
  for (unsigned int i4=0; i4<layer_4_p_fibers_d.size(); i4++) cout << layer_4_p_fibers_d[i4] << "  ";
  cout << endl;
  cout << "###########################################" << endl;
  cout << endl;
  cout << endl;
  */

  for(unsigned int i1=0; i1<layer_1_p_fibers_d.size(); i1++){
    for(int j1=0; j1<15; j1++){
      SFT_Z1_top_col.push_back(calculate_Z_position(phi_top, j1+1,layer_1_p_fibers_d[i1], 1));
      SFT_Z1_bottom_col.push_back(calculate_Z_position(phi_bottom, j1+1,layer_1_p_fibers_d[i1], 1));
      DZ1_top_col.push_back(C2_Z_top - calculate_Z_position(phi_top, j1+1,layer_1_p_fibers_d[i1], 1));
      DZ1_bottom_col.push_back(C2_Z_bottom - calculate_Z_position(phi_bottom, j1+1,layer_1_p_fibers_d[i1], 1));
    } 
    SFT_Z1_top.push_back(SFT_Z1_top_col);
    SFT_Z1_bottom.push_back(SFT_Z1_bottom_col);
    DZ1_top.push_back(DZ1_top_col);
    DZ1_bottom.push_back(DZ1_bottom_col);
    SFT_Z1_top_col.clear();
    SFT_Z1_bottom_col.clear();
    DZ1_top_col.clear();
    DZ1_bottom_col.clear();
  }

  for(unsigned int i2=0; i2<layer_2_p_fibers_d.size(); i2++){
    for(int j2=0; j2<15; j2++){
      SFT_Z2_top_col.push_back(calculate_Z_position(phi_top, j2+1,layer_2_p_fibers_d[i2], 2));
      SFT_Z2_bottom_col.push_back(calculate_Z_position(phi_bottom, j2+1,layer_2_p_fibers_d[i2], 2));
      DZ2_top_col.push_back(C2_Z_top - calculate_Z_position(phi_top, j2+1,layer_2_p_fibers_d[i2], 2));
      DZ2_bottom_col.push_back(C2_Z_bottom - calculate_Z_position(phi_bottom, j2+1,layer_2_p_fibers_d[i2], 2));
     } 
    SFT_Z2_top.push_back(SFT_Z2_top_col);
    SFT_Z2_bottom.push_back(SFT_Z2_bottom_col);
    DZ2_top.push_back(DZ2_top_col);
    DZ2_bottom.push_back(DZ2_bottom_col);
    SFT_Z2_top_col.clear();
    SFT_Z2_bottom_col.clear();
    DZ2_top_col.clear();
    DZ2_bottom_col.clear();
  }

  for(unsigned int i3=0; i3<layer_3_p_fibers_d.size(); i3++){
    for(int j3=0; j3<17; j3++){
      SFT_Z3_top_col.push_back(calculate_Z_position(phi_top, j3+1,layer_3_p_fibers_d[i3], 3));
      SFT_Z3_bottom_col.push_back(calculate_Z_position(phi_bottom, j3+1,layer_3_p_fibers_d[i3], 3));
      DZ3_top_col.push_back(C2_Z_top - calculate_Z_position(phi_top, j3+1,layer_3_p_fibers_d[i3], 3));
      DZ3_bottom_col.push_back(C2_Z_bottom - calculate_Z_position(phi_bottom, j3+1,layer_3_p_fibers_d[i3], 3));
     } 
    SFT_Z3_top.push_back(SFT_Z3_top_col);
    SFT_Z3_bottom.push_back(SFT_Z3_bottom_col);
    DZ3_top.push_back(DZ3_top_col);
    DZ3_bottom.push_back(DZ3_bottom_col);
    SFT_Z3_top_col.clear();
    SFT_Z3_bottom_col.clear();
    DZ3_top_col.clear();
    DZ3_bottom_col.clear();
  }

  for(unsigned int i4=0; i4<layer_4_p_fibers_d.size(); i4++){
    for(int j4=0; j4<17; j4++){
      SFT_Z4_top_col.push_back(calculate_Z_position(phi_top, j4+1,layer_4_p_fibers_d[i4], 4));
      SFT_Z4_bottom_col.push_back(calculate_Z_position(phi_bottom, j4+1,layer_4_p_fibers_d[i4], 4));
      DZ4_top_col.push_back(C2_Z_top - calculate_Z_position(phi_top, j4+1,layer_4_p_fibers_d[i4], 4));
      DZ4_bottom_col.push_back(C2_Z_bottom - calculate_Z_position(phi_bottom, j4+1,layer_4_p_fibers_d[i4], 4));
     } 
    SFT_Z4_top.push_back(SFT_Z4_top_col);
    SFT_Z4_bottom.push_back(SFT_Z4_bottom_col);
    DZ4_top.push_back(DZ4_top_col);
    DZ4_bottom.push_back(DZ4_bottom_col);
    SFT_Z4_top_col.clear();
    SFT_Z4_bottom_col.clear();
    DZ4_top_col.clear();
    DZ4_bottom_col.clear();
  }

  /*
  cout << endl;
  cout << "LAYER 1 " << endl; 
  for(unsigned int i1=0; i1<SFT_Z1_top.size(); i1++){
    for(unsigned int j1=0; j1<SFT_Z1_top[i1].size(); j1++){
        cout << "****************  " << i1 << "   " << j1 << "   " << SFT_Z1_top[i1][j1] << "   " << SFT_Z1_bottom[i1][j1] << "   " << C2_Z_top << "   " << C2_Z_bottom << "   " << DZ1_top[i1][j1] << "   " << DZ1_bottom[i1][j1] << endl;
     }
  }
 
  cout << endl;
  cout << "LAYER 2 " << endl; 
  for(unsigned int i2=0; i2<SFT_Z2_top.size(); i2++){
    for(unsigned int j2=0; j2<SFT_Z2_top[i2].size(); j2++){
        cout << "****************  " << i2 << "   " << j2 << "   " << SFT_Z2_top[i2][j2] << "   " << SFT_Z2_bottom[i2][j2] << "   " << C2_Z_top << "   " << C2_Z_bottom << "   "  << DZ2_top[i2][j2] << "   " << DZ2_bottom[i2][j2] <<  endl;
    }
  }

  cout << endl;
  cout << "LAYER 3 " << endl; 
  for(unsigned int i3=0; i3<SFT_Z3_top.size(); i3++){
    for(unsigned int j3=0; j3<SFT_Z3_top[i3].size(); j3++){
        cout << "****************  " << i3 << "   " << j3 << "   " << SFT_Z3_top[i3][j3] << "   " << SFT_Z3_bottom[i3][j3] << "   " << C2_Z_top << "   " << C2_Z_bottom << "   "  << DZ3_top[i3][j3] << "   " << DZ3_bottom[i3][j3] <<   endl;
    }
  }

  cout << endl;
  cout << "LAYER 4 " << endl; 
  for(unsigned int i4=0; i4<SFT_Z4_top.size(); i4++){
    for(unsigned int j4=0; j4<SFT_Z4_top[i4].size(); j4++){
        cout << "****************  " << i4 << "   " << j4 << "   " << SFT_Z4_top[i4][j4] << "   " << SFT_Z4_bottom[i4][j4] << "   " << C2_Z_top << "   " << C2_Z_bottom << "   "  << DZ4_top[i4][j4] << "   " << DZ4_bottom[i4][j4] <<   endl;
    }
  }
  cout << endl;
  cout << endl;
  */

  for(unsigned int i1=0; i1<SFT_Z1_top.size(); i1++){
    for(unsigned int j1=0; j1<SFT_Z1_top[i1].size(); j1++){
      if(abs(DZ1_top[i1][j1]) < lowest_top){
        lowest_top = abs(DZ1_top[i1][j1]);
        n_lowest_top = j1+1;
        z_lowest_top = SFT_Z1_top[i1][j1];
      }
      if(abs(DZ1_bottom[i1][j1]) < lowest_bottom){
        lowest_bottom = abs(DZ1_bottom[i1][j1]);
        n_lowest_bottom = j1+1;
        z_lowest_bottom = SFT_Z1_bottom[i1][j1];
      }
    }
    n1_lowest_top.push_back(n_lowest_top);
    n1_lowest_bottom.push_back(n_lowest_bottom);
    Z1_lowest_top.push_back(z_lowest_top);
    Z1_lowest_bottom.push_back(z_lowest_bottom);
    DZ1_lowest_top.push_back(lowest_top);
    DZ1_lowest_bottom.push_back(lowest_bottom);
    lowest_top = 999.;
    lowest_bottom = 999.;
    z_lowest_top = 999.;
    z_lowest_bottom = 999.;
    n_lowest_top = 99;
    n_lowest_bottom = 99;
  }

  for(unsigned int i2=0; i2<SFT_Z2_top.size(); i2++){
    for(unsigned int j2=0; j2<SFT_Z2_top[i2].size(); j2++){
      if(abs(DZ2_top[i2][j2]) < lowest_top){
        lowest_top = abs(DZ2_top[i2][j2]);
        n_lowest_top = j2+1;
        z_lowest_top = SFT_Z2_top[i2][j2];
     }
      if(abs(DZ2_bottom[i2][j2]) < lowest_bottom){
        lowest_bottom = abs(DZ2_bottom[i2][j2]);
        n_lowest_bottom = j2+1;
        z_lowest_bottom = SFT_Z2_bottom[i2][j2];
      }
   }
    n2_lowest_top.push_back(n_lowest_top);
    n2_lowest_bottom.push_back(n_lowest_bottom);
    Z2_lowest_top.push_back(z_lowest_top);
    Z2_lowest_bottom.push_back(z_lowest_bottom);
    DZ2_lowest_top.push_back(lowest_top);
    DZ2_lowest_bottom.push_back(lowest_bottom);
    lowest_top = 999.;
    lowest_bottom = 999.;
    z_lowest_top = 999.;
    z_lowest_bottom = 999.;
    n_lowest_top = 99;
    n_lowest_bottom = 99;
  }

  for(unsigned int i3=0; i3<SFT_Z3_top.size(); i3++){
    for(unsigned int j3=0; j3<SFT_Z3_top[i3].size(); j3++){
      if(abs(DZ3_top[i3][j3]) < lowest_top){
        lowest_top = abs(DZ3_top[i3][j3]);
        n_lowest_top = j3+1;
        z_lowest_top = SFT_Z3_top[i3][j3];
      }
      if(abs(DZ3_bottom[i3][j3]) < lowest_bottom){
        lowest_bottom = abs(DZ3_bottom[i3][j3]);
        n_lowest_bottom = j3+1;
        z_lowest_bottom = SFT_Z3_bottom[i3][j3];
      }
   }
    n3_lowest_top.push_back(n_lowest_top);
    n3_lowest_bottom.push_back(n_lowest_bottom);
    Z3_lowest_top.push_back(z_lowest_top);
    Z3_lowest_bottom.push_back(z_lowest_bottom);
    DZ3_lowest_top.push_back(lowest_top);
    DZ3_lowest_bottom.push_back(lowest_bottom);
    lowest_top = 999.;
    lowest_bottom = 999.;
    z_lowest_top = 999.;
    z_lowest_bottom = 999.;
    n_lowest_top = 99;
    n_lowest_bottom = 99;
  }

  for(unsigned int i4=0; i4<SFT_Z4_top.size(); i4++){
    for(unsigned int j4=0; j4<SFT_Z4_top[i4].size(); j4++){
      if(abs(DZ4_top[i4][j4]) < lowest_top){
        lowest_top = abs(DZ4_top[i4][j4]);
        n_lowest_top = j4+1;
        z_lowest_top = SFT_Z4_top[i4][j4];
      }  
      if(abs(DZ4_bottom[i4][j4]) < lowest_bottom){
        lowest_bottom = abs(DZ4_bottom[i4][j4]);
        n_lowest_bottom = j4+1;
        z_lowest_bottom = SFT_Z4_bottom[i4][j4];
      }
   }
    n4_lowest_top.push_back(n_lowest_top);
    n4_lowest_bottom.push_back(n_lowest_bottom);
    Z4_lowest_top.push_back(z_lowest_top);
    Z4_lowest_bottom.push_back(z_lowest_bottom);
    DZ4_lowest_top.push_back(lowest_top);
    DZ4_lowest_bottom.push_back(lowest_bottom);
    lowest_top = 999.;
    lowest_bottom = 999.;
    z_lowest_top = 999.;
    z_lowest_bottom = 999.;
    n_lowest_top = 99;
    n_lowest_bottom = 99;
  }

  /*
  cout << "LAYER 1:" << endl;
  for(unsigned int i1=0; i1<DZ1_lowest_top.size(); i1++){
    cout << n1_lowest_top[i1] << " " << layer_1_p_fibers_d[i1] << "   " << Z1_lowest_top[i1] << "  " << C2_Z_top << "  " << DZ1_lowest_top[i1] << "   ";
  }
  cout << endl;

  for(unsigned int i1=0; i1<DZ1_lowest_bottom.size(); i1++){
    cout << n1_lowest_bottom[i1] << " " << layer_1_p_fibers_d[i1] << "   " << Z1_lowest_bottom[i1] << "  " << C2_Z_bottom << "  " << DZ1_lowest_bottom[i1] << "   ";
  }

  cout << endl;
  cout << endl;

  cout << "LAYER 2:" << endl;
  for(unsigned int i2=0; i2<DZ2_lowest_top.size(); i2++){
    cout << DZ2_lowest_top[i2] << "   ";
  }
  cout << endl;

  for(unsigned int i2=0; i2<DZ2_lowest_bottom.size(); i2++){
    cout << DZ2_lowest_bottom[i2] << "   ";
  }
  cout << endl;
  cout << endl;

  cout << "LAYER 3:" << endl;
  for(unsigned int i3=0; i3<DZ3_lowest_top.size(); i3++){
    cout << DZ3_lowest_top[i3] << "   ";
  }
  cout << endl;

  for(unsigned int i3=0; i3<DZ3_lowest_bottom.size(); i3++){
    cout << DZ3_lowest_bottom[i3] << "   ";
  }
  cout << endl;
  cout << endl;

  cout << "LAYER 4:" << endl;
  for(unsigned int i4=0; i4<DZ4_lowest_top.size(); i4++){
    cout << DZ4_lowest_top[i4] << "   ";
  }
  cout << endl;

  for(unsigned int i4=0; i4<DZ4_lowest_bottom.size(); i4++){
    cout << DZ4_lowest_bottom[i4] << "   ";
  }
  cout << endl;
  cout << endl;
  */
  //cout << endl;
  //cout << " |===========|==========================================|==========================================| " << endl;
  //cout << " |   Run#    |                   TOP                    |                  BOTTOM                  |" << endl;
  //printf(" |   %4d    |                                          |                                          |\n",Run_Number);
  //  cout << " |   4560    |                                          |                                          | " << endl;
  //cout << " |===========|==========================================|==========================================| " << endl;
  //printf(" |           |            phi = %8.3f deg.           |             phi = %8.3f deg.          |\n",phi_top,phi_bottom);
  //printf(" |   Evt#    |           C2_Z = %8.3f   mm           |            C2_Z = %8.3f   mm          |\n",C2_Z_top,C2_Z_bottom);
  //printf(" |   %4d    |==========================================|==========================================|\n",evt);
  //cout << " |           |==========================================|==========================================| " << endl;
  //  cout << " |           |       p       n       Z         DZ       |       p       n       Z         DZ       |" << endl;
  //printf(" |           |       p       n       Z         DZ       |       p       n       Z         DZ       |\n");
 //cout << " |===========|==========================================|==========================================| " << endl;

  //if(layer_1_p_fibers_d.size()==0){
  //  printf(" |  LAYER 1  |                -----------               |               -----------                |\n");
  //}
  //else{
  //  printf(" |  LAYER 1  |     %5.2f    %2d   %8.3f  %8.3f     |     %5.2f    %2d   %8.3f  %8.3f     |\n", 
  //  layer_1_p_fibers_d[0], n1_lowest_top[0], Z1_lowest_top[0], DZ1_lowest_top[0], layer_1_p_fibers_d[0], n1_lowest_bottom[0], Z1_lowest_bottom[0], DZ1_lowest_bottom[0]);
  //  for(unsigned int i1=0; i1<DZ1_lowest_top.size()-1; i1++){
  //  printf(" |           |     %5.2f    %2d   %8.3f  %8.3f     |     %5.2f    %2d   %8.3f  %8.3f     |\n", 
  //  layer_1_p_fibers_d[i1+1], n1_lowest_top[i1+1], Z1_lowest_top[i1+1], DZ1_lowest_top[i1+1], layer_1_p_fibers_d[i1+1], n1_lowest_bottom[i1+1], Z1_lowest_bottom[i1+1], DZ1_lowest_bottom[i1+1]);
  //  }
  //}

  //cout << " |===========|==========================================|==========================================| " << endl;

  /*
  if(layer_2_p_fibers_d.size()==0){
    printf(" |  LAYER 2  |               -----------                |               -----------                |\n");
  }
  else{
    printf(" |  LAYER 2  |     %5.2f    %2d   %8.3f  %8.3f     |     %5.2f    %2d   %8.3f  %8.3f     |\n", 
    layer_2_p_fibers_d[0], n2_lowest_top[0], Z2_lowest_top[0], DZ2_lowest_top[0], layer_2_p_fibers_d[0], n2_lowest_bottom[0], Z2_lowest_bottom[0], DZ2_lowest_bottom[0]);
    for(unsigned int i2=0; i2<DZ2_lowest_top.size()-1; i2++){
    printf(" |           |     %5.2f    %2d   %8.3f  %8.3f     |     %5.2f    %2d   %8.3f  %8.3f     |\n", 
    layer_2_p_fibers_d[i2+1], n2_lowest_top[i2+1], Z2_lowest_top[i2+1], DZ2_lowest_top[i2+1], layer_2_p_fibers_d[i2+1], n2_lowest_bottom[i2+1], Z2_lowest_bottom[i2+1], DZ2_lowest_bottom[i2+1]);
    }
  }

  cout << " |===========|==========================================|==========================================| " << endl;
  */
  /*
  if(layer_3_p_fibers_d.size()==0){
    printf(" |  LAYER 3  |                -----------                |                -----------                |\n");
  }
  else{
    printf(" |  LAYER 3  |     %5.2f    %2d   %8.3f  %8.3f     |     %5.2f    %2d   %8.3f  %8.3f     |\n", 
    layer_3_p_fibers_d[0], n3_lowest_top[0], Z3_lowest_top[0], DZ3_lowest_top[0], layer_3_p_fibers_d[0], n3_lowest_bottom[0], Z3_lowest_bottom[0], DZ3_lowest_bottom[0]);
    for(unsigned int i3=0; i3<DZ3_lowest_top.size()-1; i3++){
    printf(" |           |     %5.2f    %2d   %8.3f  %8.3f     |     %5.2f    %2d   %8.3f  %8.3f     |\n", 
    layer_3_p_fibers_d[i3+1], n3_lowest_top[i3+1], Z3_lowest_top[i3+1], DZ3_lowest_top[i3+1], layer_3_p_fibers_d[i3+1], n3_lowest_bottom[i3+1], Z3_lowest_bottom[i3+1], DZ3_lowest_bottom[i3+1]);
    }
  }

  cout << " |===========|==========================================|==========================================| " << endl;
  */
  /*
  if(layer_4_p_fibers_d.size()==0){
    printf(" |  LAYER 4  |                -----------                |                -----------                |\n");
  }
  else{
    printf(" |  LAYER 4  |     %5.2f    %2d   %8.3f  %8.3f     |     %5.2f    %2d   %8.3f  %8.3f     |\n", 
    layer_4_p_fibers_d[0], n4_lowest_top[0], Z4_lowest_top[0], DZ4_lowest_top[0], layer_4_p_fibers_d[0], n4_lowest_bottom[0], Z4_lowest_bottom[0], DZ4_lowest_bottom[0]);
    for(unsigned int i4=0; i4<DZ4_lowest_top.size()-1; i4++){
    printf(" |           |     %5.2f    %2d   %8.3f  %8.3f     |     %5.2f    %2d   %8.3f  %8.3f     |\n", 
    layer_4_p_fibers_d[i4+1], n4_lowest_top[i4+1], Z4_lowest_top[i4+1], DZ4_lowest_top[i4+1], layer_4_p_fibers_d[i4+1], n4_lowest_bottom[i4+1], Z4_lowest_bottom[i4+1], DZ4_lowest_bottom[i4+1]);
    }
  }

  cout << " |===========|==========================================|==========================================| " << endl;
  cout << endl;
  cout << endl;
  */

  vector<double> Z1_top;  vector<double> Z1_bottom;
  vector<double> Z2_top;  vector<double> Z2_bottom;
  vector<double> Z3_top;  vector<double> Z3_bottom;
  vector<double> Z4_top;  vector<double> Z4_bottom;

  vector<double> p1_top;  vector<double> p1_bottom;
  vector<double> p2_top;  vector<double> p2_bottom;
  vector<double> p3_top;  vector<double> p3_bottom;
  vector<double> p4_top;  vector<double> p4_bottom;

  vector<int> n1_top;  vector<int> n1_bottom;
  vector<int> n2_top;  vector<int> n2_bottom;
  vector<int> n3_top;  vector<int> n3_bottom;
  vector<int> n4_top;  vector<int> n4_bottom;

  vector<double> D_Z1_top;  vector<double> D_Z1_bottom;
  vector<double> D_Z2_top;  vector<double> D_Z2_bottom;
  vector<double> D_Z3_top;  vector<double> D_Z3_bottom;
  vector<double> D_Z4_top;  vector<double> D_Z4_bottom;


  // LAYER 1
  for(unsigned int i1=0; i1<layer_1_p_fibers_d.size(); i1++){
    if(DZ1_lowest_top[i1]<= DZ1_lowest_bottom[i1]){
      n1_top.push_back(n1_lowest_top[i1]);
      p1_top.push_back(layer_1_p_fibers_d[i1]);
      Z1_top.push_back(Z1_lowest_top[i1]);
      D_Z1_top.push_back(DZ1_lowest_top[i1]);
    }  
    else{
      n1_bottom.push_back(n1_lowest_bottom[i1]);
      p1_bottom.push_back(layer_1_p_fibers_d[i1]);
      Z1_bottom.push_back(Z1_lowest_bottom[i1]);
      D_Z1_bottom.push_back(DZ1_lowest_bottom[i1]);
    }
  }

  // LAYER 2
  for(unsigned int i2=0; i2<layer_2_p_fibers_d.size(); i2++){
    if(DZ2_lowest_top[i2]<= DZ2_lowest_bottom[i2]){
      n2_top.push_back(n2_lowest_top[i2]);
      p2_top.push_back(layer_2_p_fibers_d[i2]);
      Z2_top.push_back(Z2_lowest_top[i2]);
      D_Z2_top.push_back(DZ2_lowest_top[i2]);
    }  
    else{
      n2_bottom.push_back(n2_lowest_bottom[i2]);
      p2_bottom.push_back(layer_2_p_fibers_d[i2]);
      Z2_bottom.push_back(Z2_lowest_bottom[i2]);
      D_Z2_bottom.push_back(DZ2_lowest_bottom[i2]);
    }
  } 

  // LAYER 3
  for(unsigned int i3=0; i3<layer_3_p_fibers_d.size(); i3++){
    if(DZ3_lowest_top[i3]<= DZ3_lowest_bottom[i3]){
      n3_top.push_back(n3_lowest_top[i3]);
      p3_top.push_back(layer_3_p_fibers_d[i3]);
      Z3_top.push_back(Z3_lowest_top[i3]);
      D_Z3_top.push_back(DZ3_lowest_top[i3]);
    }  
    else{
      n3_bottom.push_back(n3_lowest_bottom[i3]);
      p3_bottom.push_back(layer_3_p_fibers_d[i3]);
      Z3_bottom.push_back(Z3_lowest_bottom[i3]);
      D_Z3_bottom.push_back(DZ3_lowest_bottom[i3]);
    }
  } 


  // LAYER 4
  for(unsigned int i4=0; i4<layer_4_p_fibers_d.size(); i4++){
    if(DZ4_lowest_top[i4]<= DZ4_lowest_bottom[i4]){
      n4_top.push_back(n4_lowest_top[i4]);
      p4_top.push_back(layer_4_p_fibers_d[i4]);
      Z4_top.push_back(Z4_lowest_top[i4]);
      D_Z4_top.push_back(DZ4_lowest_top[i4]);
    }  
    else{
      n4_bottom.push_back(n4_lowest_bottom[i4]);
      p4_bottom.push_back(layer_4_p_fibers_d[i4]);
      Z4_bottom.push_back(Z4_lowest_bottom[i4]);
      D_Z4_bottom.push_back(DZ4_lowest_bottom[i4]);
    }
  } 

  /*
  for(unsigned int i1=0; i1<Z1_top.size(); i1++){
    cout << "Z1_top(" << n1_top[i1] << "," << p1_top[i1] << ")= " << Z1_top[i1] << ",  DZ1_top= " << D_Z1_top[i1] << endl;
  }
  for(unsigned int j1=0; j1<Z1_bottom.size(); j1++){
    cout << "Z1_bottom(" << n1_bottom[j1] << "," << p1_bottom[j1] << ")= " << Z1_bottom[j1] << ",  DZ1_bottom= " << D_Z1_bottom[j1] << endl;
  }
  cout << endl;

  for(unsigned int i2=0; i2<Z2_top.size(); i2++){
    cout << "Z2_top(" << n2_top[i2] << "," << p2_top[i2] << ")= " << Z2_top[i2] << ",  DZ2_top= " << D_Z2_top[i2] << endl;
  }
  for(unsigned int j2=0; j2<Z2_bottom.size(); j2++){
    cout << "Z2_bottom(" << n2_bottom[j2] << "," << p2_bottom[j2] << ")= " << Z2_bottom[j2] << ",  DZ2_bottom= " << D_Z2_bottom[j2] << endl;
  }
  cout << endl;

  for(unsigned int i3=0; i3<Z3_top.size(); i3++){
    cout << "Z3_top(" << n3_top[i3] << "," << p3_top[i3] << ")= " << Z3_top[i3] << ",  DZ3_top= " << D_Z3_top[i3] << endl;
  }
  for(unsigned int j3=0; j3<Z3_bottom.size(); j3++){
    cout << "Z3_bottom(" << n3_bottom[j3] << "," << p3_bottom[j3] << ")= " << Z3_bottom[j3] << ",  DZ3_bottom= " << D_Z3_bottom[j3] << endl;
  }
  cout << endl;

  for(unsigned int i4=0; i4<Z4_top.size(); i4++){
    cout << "Z4_top(" << n4_top[i4] << "," << p4_top[i4] << ")= " << Z4_top[i4] << ",  DZ4_top= " << D_Z4_top[i4] << endl;
  }
  for(unsigned int j4=0; j4<Z4_bottom.size(); j4++){
    cout << "Z4_bottom(" << n4_bottom[j4] << "," << p4_bottom[j4] << ")= " << Z4_bottom[j4] << ",  DZ4_bottom= " << D_Z4_bottom[j4] << endl;
  }
  cout << endl;
  */

  vector<double> Z_top_final;   vector<double> Z_bottom_final;
  double DZ1_top_min = 999.99;   double DZ1_bottom_min = 999.99;  
  double DZ2_top_min = 999.99;   double DZ2_bottom_min = 999.99;  
  double DZ3_top_min = 999.99;   double DZ3_bottom_min = 999.99;  
  double DZ4_top_min = 999.99;   double DZ4_bottom_min = 999.99;  

  int i1_top_min = 99;  int i1_bot_min = 99;
  int i2_top_min = 99;  int i2_bot_min = 99;
  int i3_top_min = 99;  int i3_bot_min = 99;
  int i4_top_min = 99;  int i4_bot_min = 99;

  Z_top_final.clear();
  Z_bottom_final.clear();


  /*
  // LAYER 1
  for(unsigned int i1=0; i1<Z1_top.size(); i1++){
    if(D_Z1_top[i1] < DZ1_top_min){
      DZ1_top_min = D_Z1_top[i1];
      i1_top_min = i1;
    }
  }
  for(unsigned int j1=0; j1<Z1_bottom.size(); j1++){
    if(D_Z1_bottom[j1] < DZ1_bottom_min){
      DZ1_bottom_min = D_Z1_bottom[j1];
      i1_bot_min = j1;
    }
  }
  if(Z1_top.size()!=0) Z_top_final.push_back(Z1_top[i1_top_min]); 
  if(Z1_bottom.size()!=0) Z_bottom_final.push_back(Z1_bottom[i1_bot_min]); 

  
  //LAYER 2
  for(unsigned int i2=0; i2<Z2_top.size(); i2++){
    if(D_Z2_top[i2] < DZ2_top_min){
      DZ2_top_min = D_Z2_top[i2];
      i2_top_min = i2;
    }
  }
  for(unsigned int j2=0; j2<Z2_bottom.size(); j2++){
    if(D_Z2_bottom[j2] < DZ2_bottom_min){
      DZ2_bottom_min = D_Z2_bottom[j2];
      i2_bot_min = j2;
    }
  }
  if(Z2_top.size()!=0) Z_top_final.push_back(Z2_top[i2_top_min]); 
  if(Z2_bottom.size()!=0) Z_bottom_final.push_back(Z2_bottom[i2_bot_min]); 

  
  
  //LAYER 3
  for(unsigned int i3=0; i3<Z3_top.size(); i3++){
    //cout << "D_Z3 :"  << D_Z3_top[i3] << "  " << DZ3_top_min << endl;

    if(D_Z3_top[i3] < DZ3_top_min){
      DZ3_top_min = D_Z3_top[i3];
      i3_top_min = i3;
    }
  }
  //cout << "IIII : " << i3_top_min << endl;

  for(unsigned int j3=0; j3<Z3_bottom.size(); j3++){
    if(D_Z3_bottom[j3] < DZ3_bottom_min){
      DZ3_bottom_min = D_Z3_bottom[j3];
      i3_bot_min = j3;
    }
  }
  if(Z3_top.size()!=0) Z_top_final.push_back(Z3_top[i3_top_min]);
  if(Z3_bottom.size()!=0) Z_bottom_final.push_back(Z3_bottom[i3_bot_min]);
  
  for(unsigned int i=0; i<Z_top_final.size(); i++){
    //cout << "Z_top_final: " << Z_top_final[i] << endl;
  }

  //LAYER 4
  for(unsigned int i4=0; i4<Z4_top.size(); i4++){
    if(D_Z4_top[i4] < DZ4_top_min){
      DZ4_top_min = D_Z4_top[i4];
      i4_top_min = i4;
    }
  }
  for(unsigned int j4=0; j4<Z4_bottom.size(); j4++){
    if(D_Z4_bottom[j4] < DZ4_bottom_min){
      DZ4_bottom_min = D_Z4_bottom[j4];
      i4_bot_min = j4;
    }
  }
  if(Z4_top.size()!=0) Z_top_final.push_back(Z4_top[i4_top_min]);
  if(Z4_bottom.size()!=0) Z_bottom_final.push_back(Z4_bottom[i4_bot_min]);


  //cout << "TEST : " << i1_top_min << "  " << i1_bot_min << endl;
  //cout << "TEST : " << Z1_top[i1_top_min] << "  " << Z1_bottom[i1_bot_min] << endl;
  //cout << "TEST : " << i1_bot_min << "  " << Z1_bottom[0] << "  " << Z1_bottom[i1_bot_min] << endl;
  //Z_top_final.push_back(Z3_top[i3_top_min]);
  //Z_bottom_final.push_back(Z1_bottom[i1_bot_min]);
  */

    // LAYER 1
  for(unsigned int i1=0; i1<Z1_top.size(); i1++){
    if(D_Z1_top[i1] < DZ1_top_min){
      DZ1_top_min = D_Z1_top[i1];
      i1_top_min = i1;
    }
  }
  for(unsigned int j1=0; j1<Z1_bottom.size(); j1++){
    if(D_Z1_bottom[j1] < DZ1_bottom_min){
      DZ1_bottom_min = D_Z1_bottom[j1];
      i1_bot_min = j1;
    }
  }
  if(Z1_top.size()!=0) Z_top_final.push_back(Z1_top[i1_top_min]); 
  if(Z1_bottom.size()!=0) Z_bottom_final.push_back(Z1_bottom[i1_bot_min]); 

  
  //LAYER 2
  for(unsigned int i2=0; i2<Z2_top.size(); i2++){
    if(D_Z2_top[i2] < DZ2_top_min){
      DZ2_top_min = D_Z2_top[i2];
      i2_top_min = i2;
    }
  }
  for(unsigned int j2=0; j2<Z2_bottom.size(); j2++){
    if(D_Z2_bottom[j2] < DZ2_bottom_min){
      DZ2_bottom_min = D_Z2_bottom[j2];
      i2_bot_min = j2;
    }
  }
  if(Z2_top.size()!=0) Z_top_final.push_back(Z2_top[i2_top_min]); 
  if(Z2_bottom.size()!=0) Z_bottom_final.push_back(Z2_bottom[i2_bot_min]); 

  
  
  //LAYER 3
  for(unsigned int i3=0; i3<Z3_top.size(); i3++){
    if(D_Z3_top[i3] < DZ3_top_min){
      DZ3_top_min = D_Z3_top[i3];
      i3_top_min = i3;
    }
  }
  for(unsigned int j3=0; j3<Z3_bottom.size(); j3++){
    if(D_Z3_bottom[j3] < DZ3_bottom_min){
      DZ3_bottom_min = D_Z3_bottom[j3];
      i3_bot_min = j3;
    }
  }
  if(Z3_top.size()!=0) Z_top_final.push_back(Z3_top[i3_top_min]);
  if(Z3_bottom.size()!=0) Z_bottom_final.push_back(Z3_bottom[i3_bot_min]);

  //LAYER 4
  for(unsigned int i4=0; i4<Z4_top.size(); i4++){
    if(D_Z4_top[i4] < DZ4_top_min){
      DZ4_top_min = D_Z4_top[i4];
      i4_top_min = i4;
    }
  }
  for(unsigned int j4=0; j4<Z4_bottom.size(); j4++){
    if(D_Z4_bottom[j4] < DZ4_bottom_min){
      DZ4_bottom_min = D_Z4_bottom[j4];
      i4_bot_min = j4;
    }
  }
  if(Z4_top.size()!=0) Z_top_final.push_back(Z4_top[i4_top_min]);
  if(Z4_bottom.size()!=0) Z_bottom_final.push_back(Z4_bottom[i4_bot_min]);



  //Z_top_final.push_back(Z3_top[i_top_min]);
  //Z_bottom_final.push_back(Z3_bottom[i_bot_min]);

  //cout << "TEST:  " << i_top_min << "   " << Z1_top[i_top_min] << "   " << Z_top_final.size() << "   " << Z_top_final[0] << endl;
  //if(Z_top_final.size()>0){
  //  cout << " Z_top :    ";
  //  for(unsigned int i=0; i<Z_top_final.size(); i++) cout << Z_top_final[i] << "  ";
  //  cout << endl;
  //}
  //if(Z_bottom_final.size()>0){
  //  cout << " Z_bottom : ";
  //  for(unsigned int j=0; j<Z_bottom_final.size(); j++) cout << Z_bottom_final[j] << "  ";
  //  cout << endl;
  //}
  
  //cout << endl;
  //cout << endl;



  //cout << "SEBB: " << Z_top_FINAL << endl;
  //cout << "SIZE : " <<  Z_top_final.size() << "  " << Z_bottom_final.size() << endl;

  if(Z_top_final.size()!=0){
    for(unsigned int i=0; i<Z_top_final.size(); i++){
      sum_top += Z_top_final[i];
    }
    Z_top_FINAL =  sum_top/Z_top_final.size();
    DZ_top_FINAL = Z_top_FINAL - C2_Z_top;
  }

  if(Z_bottom_final.size()!=0){
    for(unsigned int i=0; i<Z_bottom_final.size(); i++){
      sum_bottom += Z_bottom_final[i];
    }
    Z_bottom_FINAL =  sum_bottom/Z_bottom_final.size(); 
    DZ_bottom_FINAL = Z_bottom_FINAL - C2_Z_bottom;
  }


  vec_output.clear();
  vec_output.push_back(Z_top_FINAL);
  vec_output.push_back(Z_bottom_FINAL);
  vec_output.push_back(abs(DZ_top_FINAL));
  vec_output.push_back(abs(DZ_bottom_FINAL));

  //cout << "SEB : " << Z_top_FINAL << endl;

  //if(Z_top_FINAL>500) cout << "Z_TOP_FINAL =   -------  " << endl;
  //cout << " ===============================================" << endl;
  //cout << endl;
  //if(Z_top_FINAL<500) cout << " TOP:      Z = " << Z_top_FINAL << "  ;  Delta_Z = " << DZ_top_FINAL << endl;
  //if(Z_top_FINAL<500) printf("   TOP:     Z = %8.3f  ;  Delta_Z = %7.3f\n",Z_top_FINAL,DZ_top_FINAL);

  //if(Z_bottom_FINAL>500) cout << "Z_BOTTOM_FINAL =   -------  " << endl;
  //if(Z_bottom_FINAL<500) cout << " BOTTOM:   Z = " << Z_bottom_FINAL << "  ;  Delta_Z = " << DZ_bottom_FINAL << endl;
  //if(Z_bottom_FINAL<500) printf("   BOTTOM:  Z = %8.3f  ;  Delta_Z = %7.3f\n",Z_bottom_FINAL,DZ_bottom_FINAL);
  //cout << endl;
  //cout << " ===============================================" << endl;

  //cout << endl;
  //cout << endl;
  //for(unsigned int i=0; i<vec_output.size(); i++){
  //  cout << "TEST : " << vec_output[i] << endl;
  //}

  return vec_output;

 } // End of Ana3
