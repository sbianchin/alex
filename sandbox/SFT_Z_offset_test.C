#include <stdio.h>
#include <math.h>
#include <iostream>

using namespace std;

double calculate_Z_position(double phi = 90., int N = 7 , int p = 4, int layer = 4, bool z_function_test = true){
  double Z_position = 0;
  
  if(z_function_test)
    phi *= M_PI/180.;

  // Z_0 CONSTANT OFFSET FOR EACH LAYER
  /////////////////////////////////////////////////////////////
  const double z_offset_const[4] = {-131.0, -136.5, -138.68, -135.83};

  // Number of fibers in each ribbon/layer
  const int fibers_per_ribbon[4] = {15, 15, 17, 17};

  // Number of turns per layer
  const int num_turns[4] = {17, 17, 15, 15};

  // Angle of wrapping
  const double alpha[4] = {3.6*(M_PI/180.0), 3.6*(M_PI/180.0), 3.9*(M_PI/180.0), 3.9*(M_PI/180.0)};

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
               + (p-1)*ribbon_width[0]/fibers_per_ribbon[0]
               + fiber_diameter/cos(alpha[1]);   
  }
  else if(layer == 2){
    Z_position = z_offset_const[1] + layer_radius[1]*phi*tan(alpha[1]) 
               + (N-1)*ribbon_width[1]
               + (p-1)*ribbon_width[1]/fibers_per_ribbon[1]
               + fiber_diameter/cos(alpha[1]);     

  }
  else if(layer == 3){
    Z_position = z_offset_const[2] + layer_radius[2]*(2*M_PI-phi)*tan(alpha[2]) 
               + (N-1)*ribbon_width[2]
               + (fibers_per_ribbon[2]-p)*ribbon_width[2]/fibers_per_ribbon[2]
               + fiber_diameter/cos(alpha[2]); 
  }
  else if(layer == 4){
    Z_position = z_offset_const[3] + layer_radius[3]*(2*M_PI-phi)*tan(alpha[3]) 
               + (N-1)*ribbon_width[3]
               + (fibers_per_ribbon[3]-p)*ribbon_width[3]/fibers_per_ribbon[3]
               + fiber_diameter/cos(alpha[3]);
  }

  if(z_function_test){
    cout << "Layer = " << layer << endl;
    cout << "N = " << N << endl;
    cout << "p = " << p << endl;
    cout << "Z = " << Z_position << endl;
  }
  

  return Z_position;
}


double calculate_C_position(double phi = 90., int N = 7 , int p = 4, int layer = 4, double z_pos = 0){
  double Z_position = 0;
  double z_const = 0;

  // Z_0 CONSTANT OFFSET FOR EACH LAYER
  /////////////////////////////////////////////////////////////
  //const double z_offset_const[4] = {-89.15, -79.76, -121.88-16.8, -135.83};



  // Number of fibers in each ribbon/layer
  const int fibers_per_ribbon[4] = {15, 15, 17, 17};

  // Number of turns per layer
  const int num_turns[4] = {17, 17, 15, 15};

  // Angle of wrapping
  const double alpha[4] = {3.6*(M_PI/180.0), 3.6*(M_PI/180.0), 3.9*(M_PI/180.0), 3.9*(M_PI/180.0)};

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
    // Z_position = z_offset_const[0] + layer_radius[0]*phi*tan(alpha[0]) 
    //            + (N-1)*ribbon_width[0]
    //            + (p-1)*ribbon_width[0]/fibers_per_ribbon[0]
    //            + fiber_diameter/cos(alpha[1]); 

    z_const = z_pos -(layer_radius[0]*phi*tan(alpha[0]) 
                      + (N-1)*ribbon_width[0]
                      + (p-1)*ribbon_width[0]/fibers_per_ribbon[0]
                      + fiber_diameter/cos(alpha[1]));             
  }
  else if(layer == 2){
    // Z_position = z_offset_const[1] + layer_radius[1]*phi*tan(alpha[1]) 
    //            + (N-1)*ribbon_width[1]
    //            + (p-1)*ribbon_width[1]/fibers_per_ribbon[1]
    //            + fiber_diameter/cos(alpha[1]);     

    z_const = z_pos -(layer_radius[1]*phi*tan(alpha[1]) 
                      + (N-1)*ribbon_width[1]
                      + (p-1)*ribbon_width[1]/fibers_per_ribbon[1]
                      + fiber_diameter/cos(alpha[1]));         

  }

  

  return z_const;
}

void SFT_Z_offset_test(double phi, int p1, int p2, int p3, int p4){


  // Number of turns per layer
  const int num_turns[4] = {17, 17, 15, 13};



  // Convert phi to radians:
  phi *= M_PI/180;
 
  // Array to store all Z positions
  double Z_positions[4][20] = {0};

  // selected p's in an array
  const int sft_p_fiber[4] = {p1, p2, p3, p4};                                  




  for(int N = 0; N < num_turns[0]; N++){        
    Z_positions[0][N] = calculate_Z_position(phi, N, p1, 1, false);
  }

  //Find Z positions for layer 2;
  for(int N = 0; N < num_turns[1]; N++){
    Z_positions[1][N] = calculate_Z_position(phi, N, p2, 2, false);
  }

  // Find Z positions for layer 3;
  for(int N = 0; N < num_turns[2]; N++){
    Z_positions[2][N] = calculate_Z_position(phi, N, p3, 3, false); 
  }

  // Find Z positions for layer 4;
  for(int N = 0; N < num_turns[3]; N++){
    Z_positions[3][N] = calculate_Z_position(phi, N, p4, 4, false);
  }
  

  double Z_34_avg[20] = {0.0};
  double c1[20] = {0.0};
  double c2[20] = {0.0};

  for(int i = 0; i < 13; i++){
    c1[i] = calculate_C_position(phi, i, p1, 1, Z_34_avg[i]);
    c2[i] = calculate_C_position(phi, i, p2, 2, Z_34_avg[i]);
  }


  for(int i = 0; i<13; i++){
    Z_34_avg[i] = (Z_positions[2][i] + Z_positions[3][i])/2.0;
  }

  cout << " N      Layer 1      Layer 2     Layer 3    Layer 4     Z34 avg       C1          C2"<< endl;

  for(int N = 0; N<20; N++){
    printf("%2d     %6.1f       %6.1f       %6.1f     %6.1f      %6.1f     %6.1f      %6.1f\n", N, Z_positions[0][N],
                                                               Z_positions[1][N],
                                                               Z_positions[2][N],
                                                               Z_positions[3][N],
                                                               Z_34_avg[N],
                                                               c1[N],
                                                               c2[N]);
  }
}