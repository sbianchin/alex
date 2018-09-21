#ifndef __CINT__
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>
#include <TROOT.h>
#include "TSystem.h"
#include "TFile.h"
#include "TProfile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TPolyLine.h"
#include "TLine.h"
#include "TArrow.h"
#include "TStyle.h"  
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TBranch.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TEllipse.h"
#include "TMarker.h"
#include "ANAPATH.h"
#include "CommonParameters.h"
#include "ADC_Thresholds.h"
#include "TDC_Windows.h"
#include "Cuts_and_Windows.h"
#include "MWPC_Thr.h"
#endif

#include "intersect.cxx"
#include <algorithm>

using namespace std;


// All lengths in millimeters

double calculate_Z_position(double phi = 90., int N = 7 , double p = 4.0, int layer = 4){
  double Z_position = 0.0;
  
  phi *= M_PI/180.0;

  // Z_0 CONSTANT OFFSET FOR EACH LAYER
  /////////////////////////////////////////////////////////////
  const float M_PI = 3.14195;

  const double z_offset_const[4] = {-131.0, -136.5, -138.68, -135.83};

  // Number of fibers in each ribbon/layer
  const int fibers_per_ribbon[4] = {15, 15, 17, 17};

  // Angle of wrapping
  const double alpha[4] = {3.57*(M_PI/180.0), 3.5*(M_PI/180.0), 3.87*(M_PI/180.0), 3.79*(M_PI/180.0)};
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
               + fiber_diameter/cos(alpha[1]);   
  }
  else if(layer == 2){
    Z_position = z_offset_const[1] + layer_radius[1]*phi*tan(alpha[1]) 
               + (N-1)*ribbon_width[1]
               + (p-1.0)*ribbon_width[1]/fibers_per_ribbon[1]
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

  if(p <= 0)
    Z_position = -1;
  

  return Z_position;
}

double *calculate_min_Z_position(double phi, double p1, double p2, double p3, double p4, double return_array[20]){


  // Convert phi to radians:
  //phi *= M_PI/180;

  // Set p1_p2 / p3_p4 flags == 1
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


  
  // // Compute the Z position with smallest range over all layers

  if(p1 > 0 && p2 > 0 && p3 > 0 && p4 > 0){ // Case 1: all p-fibers have a value
    for(int i = 0; i<num_turns[0];i++){
      for(int j = 0; j<num_turns[1];j++){
        for(int k = 0; k<num_turns[2];k++){
          for(int h = 0; h<num_turns[3];h++){

            n[0] = fabs(Z_positions[0][i] - Z_positions[1][j]);
            n[1] = fabs(Z_positions[0][i] - Z_positions[2][k]);
            n[2] = fabs(Z_positions[0][i] - Z_positions[3][h]);
            n[3] = fabs(Z_positions[1][j] - Z_positions[2][k]);
            n[4] = fabs(Z_positions[1][j] - Z_positions[3][h]);
            n[5] = fabs(Z_positions[2][k] - Z_positions[3][h]);

            if(n[0] < smallest_dist &&
               n[1] < smallest_dist && 
               n[2] < smallest_dist && 
               n[3] < smallest_dist && 
               n[4] < smallest_dist && 
               n[5] < smallest_dist){

               layer_1_low_index = i;
               layer_2_low_index = j;
               layer_3_low_index = k;
               layer_4_low_index = h;
               
               smallest_dist_temp = 0.; 

               for(int t = 0; t<6; t++){
                  if(n[t] > smallest_dist_temp)
                    smallest_dist_temp = n[t];
               }

               smallest_dist = smallest_dist_temp;
               return_array[0] = (Z_positions[0][layer_1_low_index] 
                                + Z_positions[1][layer_2_low_index] 
                                + Z_positions[2][layer_3_low_index] 
                                + Z_positions[3][layer_4_low_index])/4.0;    


               return_array[1] = smallest_dist;
               return_array[2] = i;    
               return_array[3] = j;    
               return_array[4] = k;    
               return_array[5] = h;    
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
  }

  else if(p1 < 0 && p2 > 0 && p3 > 0 && p4 > 0){ // Case 2: p1 has no value
    for(int j = 0; j<num_turns[1];j++){
      for(int k = 0; k<num_turns[2];k++){
        for(int h = 0; h<num_turns[3];h++){

          n[3] = fabs(Z_positions[1][j] - Z_positions[2][k]);
          n[4] = fabs(Z_positions[1][j] - Z_positions[3][h]);
          n[5] = fabs(Z_positions[2][k] - Z_positions[3][h]);

          if(n[3] < smallest_dist && 
             n[4] < smallest_dist && 
             n[5] < smallest_dist){

             layer_1_low_index = -1;
             layer_2_low_index = j;
             layer_3_low_index = k;
             layer_4_low_index = h;
             
             smallest_dist_temp = 0.; 

             
             if(n[3] > smallest_dist_temp)
               smallest_dist_temp = n[3];
             if(n[4] > smallest_dist_temp)
               smallest_dist_temp = n[4];
             if(n[5] > smallest_dist_temp)
               smallest_dist_temp = n[5];

             smallest_dist = smallest_dist_temp;
             return_array[0] = (Z_positions[1][layer_2_low_index] 
                              + Z_positions[2][layer_3_low_index] 
                              + Z_positions[3][layer_4_low_index])/3.0;    
             return_array[1] = smallest_dist;
             return_array[2] = -1;    
             return_array[3] = j;    
             return_array[4] = k;    
             return_array[5] = h;    
             return_array[6] = 0;    
             return_array[7] = p2;    
             return_array[8] = p3;    
             return_array[9] = p4;    
             return_array[10] = -1;
             return_array[11] = Z_positions[1][layer_2_low_index];
             return_array[12] = Z_positions[2][layer_3_low_index];
             return_array[13] = Z_positions[3][layer_4_low_index];
          }
        }
      }
    }
  }

  else if(p1 > 0 && p2 < 0 && p3 > 0 && p4 > 0){ // Case 3: p2 has no value
    for(int i = 0; i<num_turns[0];i++){
        for(int k = 0; k<num_turns[2];k++){
          for(int h = 0; h<num_turns[3];h++){

            n[1] = fabs(Z_positions[0][i] - Z_positions[2][k]);
            n[2] = fabs(Z_positions[0][i] - Z_positions[3][h]);
            n[5] = fabs(Z_positions[2][k] - Z_positions[3][h]);

            if(n[1] < smallest_dist && 
               n[2] < smallest_dist && 
               n[5] < smallest_dist){

               layer_1_low_index = i;
               layer_2_low_index = -1;
               layer_3_low_index = k;
               layer_4_low_index = h;
               
               smallest_dist_temp = 0.; 

               if(n[1] > smallest_dist_temp)
                 smallest_dist_temp = n[1];
               if(n[2] > smallest_dist_temp)
                 smallest_dist_temp = n[2];
               if(n[5] > smallest_dist_temp)
                 smallest_dist_temp = n[5];

               smallest_dist = smallest_dist_temp;
               return_array[0] = (Z_positions[0][layer_1_low_index] 
                                + Z_positions[2][layer_3_low_index] 
                                + Z_positions[3][layer_4_low_index])/3.0;    
               return_array[1] = smallest_dist;
               return_array[2] = i;    
               return_array[3] = -1;    
               return_array[4] = k;    
               return_array[5] = h;    
               return_array[6] = p1;    
               return_array[7] = 0;    
               return_array[8] = p3;    
               return_array[9] = p4;    
               return_array[10] = Z_positions[0][layer_1_low_index];
               return_array[11] = -1;
               return_array[12] = Z_positions[2][layer_3_low_index];
               return_array[13] = Z_positions[3][layer_4_low_index];
            }
          }
        }
      
    }
  }

  else if(p1 > 0 && p2 > 0 && p3 < 0 && p4 > 0){ // Case 4: p3 has no value
    for(int i = 0; i<num_turns[0];i++){
      for(int j = 0; j<num_turns[1];j++){
          for(int h = 0; h<num_turns[3];h++){

            n[0] = fabs(Z_positions[0][i] - Z_positions[1][j]);
            n[2] = fabs(Z_positions[0][i] - Z_positions[3][h]);
            n[4] = fabs(Z_positions[1][j] - Z_positions[3][h]);

            if(n[0] < smallest_dist &&
               n[2] < smallest_dist && 
               n[4] < smallest_dist){

               layer_1_low_index = i;
               layer_2_low_index = j;
               layer_3_low_index = -1;
               layer_4_low_index = h;
               
               smallest_dist_temp = 0.; 

               if(n[0] > smallest_dist_temp)
                 smallest_dist_temp = n[0];
               if(n[2] > smallest_dist_temp)
                 smallest_dist_temp = n[2];
               if(n[4] > smallest_dist_temp)
                 smallest_dist_temp = n[4];

               smallest_dist = smallest_dist_temp;
               return_array[0] = (Z_positions[0][layer_1_low_index] 
                                + Z_positions[1][layer_2_low_index] 
                                + Z_positions[3][layer_4_low_index])/3.0;    

               return_array[1] = smallest_dist;
               return_array[2] = i;    
               return_array[3] = j;    
               return_array[4] = -1;    
               return_array[5] = h;    
               return_array[6] = p1;    
               return_array[7] = p2;    
               return_array[8] = -1;    
               return_array[9] = p4;    
               return_array[10] = Z_positions[0][layer_1_low_index];
               return_array[11] = Z_positions[1][layer_2_low_index];
               return_array[12] = -1;
               return_array[13] = Z_positions[3][layer_4_low_index];
            }
          }
        
      }
    }
  }

  else if(p1 > 0 && p2 > 0 && p3 > 0 && p4 < 0){ // Case 4: p4 has no value
    for(int i = 0; i<num_turns[0];i++){
      for(int j = 0; j<num_turns[1];j++){
        for(int k = 0; k<num_turns[2];k++){

          n[0] = fabs(Z_positions[0][i] - Z_positions[1][j]);
          n[1] = fabs(Z_positions[0][i] - Z_positions[2][k]);
          n[3] = fabs(Z_positions[1][j] - Z_positions[2][k]);

          if(n[0] < smallest_dist &&
             n[1] < smallest_dist && 
             n[3] < smallest_dist){

             layer_1_low_index = i;
             layer_2_low_index = j;
             layer_3_low_index = k;
             layer_4_low_index = -1;
             
             smallest_dist_temp = 0.; 

             if(n[0] > smallest_dist_temp)
               smallest_dist_temp = n[0];
             if(n[1] > smallest_dist_temp)
               smallest_dist_temp = n[1];
             if(n[3] > smallest_dist_temp)
               smallest_dist_temp = n[3];

             smallest_dist = smallest_dist_temp;
             return_array[0] = (Z_positions[0][layer_1_low_index] 
                              + Z_positions[1][layer_2_low_index] 
                              + Z_positions[2][layer_3_low_index])/3.0;    
             return_array[1] = smallest_dist;
             return_array[2] = i;    
             return_array[3] = j;    
             return_array[4] = k;    
             return_array[5] = -1;    
             return_array[6] = p1;    
             return_array[7] = p2;    
             return_array[8] = p3;    
             return_array[9] = -1;    
             return_array[10] = Z_positions[0][layer_1_low_index];
             return_array[11] = Z_positions[1][layer_2_low_index];
             return_array[12] = Z_positions[2][layer_3_low_index];
             return_array[13] = -1;
          }
        }
      }
    }
  }   

  // Cases: 1 hit in each helicity
  //        1) 0 1 0 1
  //        2) 0 1 1 0
  //        3) 1 0 0 1
  //        4) 1 0 1 0


  else if(p1 < 0 && p2 > 0 && p3 < 0 && p4 > 0){ // 1)
    for(int j = 0; j<num_turns[1];j++){
      for(int h = 0; h<num_turns[3];h++){

        n[4] = fabs(Z_positions[1][j] - Z_positions[3][h]);

        if(n[4] < smallest_dist){

           layer_1_low_index = -1;
           layer_2_low_index = j;
           layer_3_low_index = -1;
           layer_4_low_index = h;
          

           smallest_dist = n[4];

           return_array[0] = (Z_positions[1][layer_2_low_index] 
                            + Z_positions[3][layer_4_low_index])/2.0;    

           return_array[1] = smallest_dist;
           return_array[2] = -1;    
           return_array[3] = j;    
           return_array[4] = -1;    
           return_array[5] = h;    
           return_array[6] = -1;    
           return_array[7] = p2;    
           return_array[8] = -1;    
           return_array[9] = p4;    
           return_array[10] = -1;
           return_array[11] = Z_positions[1][layer_2_low_index];
           return_array[12] = -1;
           return_array[13] = Z_positions[3][layer_4_low_index];
        }
      }
    }
  }


  else if(p1 < 0 && p2 > 0 && p3 > 0 && p4 < 0){ //  2) 0 1 1 0
    for(int j = 0; j<num_turns[1];j++){
      for(int k = 0; k<num_turns[2];k++){
      
        n[3] = fabs(Z_positions[1][j] - Z_positions[2][k]);

        if(n[3] < smallest_dist){

           layer_1_low_index = -1;
           layer_2_low_index = j;
           layer_3_low_index = k;
           layer_4_low_index = -1;

           smallest_dist = n[3];

           return_array[0] = (Z_positions[1][layer_2_low_index] 
                            + Z_positions[2][layer_3_low_index])/2.0;    

           return_array[1] = smallest_dist;
           return_array[2] = -1;    
           return_array[3] = j;    
           return_array[4] = k;    
           return_array[5] = -1;    
           return_array[6] = -1;    
           return_array[7] = p2;    
           return_array[8] = p3;    
           return_array[9] = -1;    
           return_array[10] = -1;
           return_array[11] = Z_positions[1][layer_2_low_index];
           return_array[12] = Z_positions[2][layer_3_low_index];
           return_array[13] = -1;
        }
      }
    }
  }


  //        3) 1 0 0 1

  else if(p1 > 0 && p2 < 0 && p3 < 0 && p4 > 0){ // Case 1: all p-fibers have a value
    for(int i = 0; i<num_turns[0];i++){
      for(int h = 0; h<num_turns[3];h++){

        n[2] = fabs(Z_positions[0][i] - Z_positions[3][h]);

        if(n[2] < smallest_dist){

           layer_1_low_index = i;
           layer_2_low_index = -1;
           layer_3_low_index = -1;
           layer_4_low_index = h;


           smallest_dist = n[2];

           return_array[0] = (Z_positions[0][layer_1_low_index] 
                            + Z_positions[3][layer_4_low_index])/2.0;  

           return_array[1] = smallest_dist;
           return_array[2] = i;    
           return_array[3] = -1;    
           return_array[4] = -1;    
           return_array[5] = h;    
           return_array[6] = p1;    
           return_array[7] = -1;    
           return_array[8] = -1;    
           return_array[9] = p4;    
           return_array[10] = Z_positions[0][layer_1_low_index];
           return_array[11] = -1;
           return_array[12] = -1;
           return_array[13] = Z_positions[3][layer_4_low_index];
        }
      }
    }
  }
 


  //        4) 1 0 1 0

  else if(p1 > 0 && p2 < 0 && p3 > 0 && p4 < 0){ // Case 1: all p-fibers have a value
    for(int i = 0; i<num_turns[0];i++){
      for(int k = 0; k<num_turns[2];k++){
    
        n[1] = fabs(Z_positions[0][i] - Z_positions[2][k]);

        if(n[1] < smallest_dist){

           layer_1_low_index = i;
           layer_2_low_index = -1;
           layer_3_low_index = k;
           layer_4_low_index = -1;

           smallest_dist = n[1];

           return_array[0] = (Z_positions[0][layer_1_low_index] 
                            + Z_positions[3][layer_4_low_index])/2.0;    
           return_array[1] = smallest_dist;
           return_array[2] = i;    
           return_array[3] = -1;    
           return_array[4] = k;    
           return_array[5] = -1;    
           return_array[6] = p1;    
           return_array[7] = -1;    
           return_array[8] = p3;    
           return_array[9] = -1;    
           return_array[10] = Z_positions[0][layer_1_low_index];
           return_array[11] = -1;
           return_array[12] = Z_positions[2][layer_3_low_index];
           return_array[13] = -1;
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
  return_array[14] = 2;
  return_array[15] = 2;

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
    p_hits[i]= 0;
  }

  for(vector<int>::iterator it = p_fibers.begin(); it != p_fibers.end(); it++){
    p_hits[*it - 1]++;
  }


  int p_cluster_separation = 1;

  // Find clustering of p fibers
  ////////////////////////////
  bool first_cluster = true;
  int cluster_spacing = 0;
  int cluster_length_count = 0;

  vector<int> p_cluster_index; // Hold the starting index of each cluster.
  vector<int> p_cluster_length; // Hold the length of each cluster. 

  int p_clusters = 0;

  for(int i = 0; i<layer_length; i++){
    if(p_hits[i] > 0 && first_cluster){
      cluster_spacing = p_cluster_separation + 1;
      first_cluster = false;
    }

    if(p_hits[i] > 0){
      if(cluster_spacing > p_cluster_separation){
        p_clusters++;
        p_cluster_index.push_back(i);
      }
      cluster_length_count++;
      cluster_spacing = 0;

      if(i+2 < layer_length && p_hits[i+1] <= 0 && p_hits[i+2] <= 0){
        p_cluster_length.push_back(cluster_length_count);
      }
      else if(i + 2 == layer_length && p_hits[i+1] <= 0){
        p_cluster_length.push_back(cluster_length_count);
      }
      else if(i + 1 == layer_length){
        p_cluster_length.push_back(cluster_length_count);
      }

    }
    else{
      cluster_spacing++;
      if(i != 0 && p_hits[i-1] <= 0){
        cluster_length_count = 0;
      }
      else if ( i != 0){
        cluster_length_count++;
      }
    }
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
    centroid_sum = 0;
    centroid_den = 0;

    if(p_cluster_index_arr[i] < -1000)
      break;

    for(int j = p_cluster_index_arr[i]; j<p_cluster_index_arr[i] + p_cluster_length_arr[i]; j++){
      if(p_hits[j] > 0){
        centroid_sum += p_hits[j];
        centroid_den++;
      }
    }

    p_centroids[i] = centroid_sum/centroid_den;
  }



  for(int i = 0; i < layer_length; i++){
    if(p_centroids[i] < -1000)
      break;

    cout << "p centroid = " << p_centroids[i] << endl;
  }


  for(vector<int>::iterator it = p_fibers.begin(); it != p_fibers.end(); it++){
    p_fibers_d.push_back(double(*it));
  }

  // uncomment for flagging bad sft 
  // *bad_sft_flag = true;

  return p_fibers_d;
}

//                                          degrees
void SFT_Test_T(int Run_Number, int evt, double phi){

  const int sft_hit_count_max = 3; // REDEFINE IN COMMON PARAMATERS  


  int L1_sft_count = 0;
  int L2_sft_count = 0;
  int L3_sft_count = 0;
  int L4_sft_count = 0;



  Int_t adc_high_sft[128];          Int_t ADC_High_SFT[128];        
  Int_t adc_low_sft[128];           Int_t ADC_Low_SFT[128];   
  Int_t tdc_le_sft[128][16];        Int_t TDC_LE_SFT[128];         
  Int_t tdc_te_sft[128][16];        Int_t TDC_TE_SFT[128];  

  Int_t HG_SFT_ADC_Thr[128] = {0};
  Int_t LG_SFT_ADC_Thr[128] = {0};

  
  for(int i=0; i<128; i++)  HG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_HG[i]) + SFT_ADC_Thr_HG_Offset;
  for(int i=0; i<128; i++)  LG_SFT_ADC_Thr[i] = round(SFT_ADC_Thr_LG[i]) + SFT_ADC_Thr_LG_Offset;

  Double_t ADC_High_SFT_corr[128]; 

  Int_t has_TDC_SFT_hit[128] = {0};


  char source_mapping[] = "SFT_Mapping_Oct14.txt";  // Mapping file !!!

  char file_mapping[200];
  sprintf(file_mapping,"../Mapping");

  char par_finput[200];
  sprintf(par_finput,"%s/%s",file_mapping,source_mapping);

  Int_t par_temp[2][128];
  ifstream fdat(par_finput,ios::in);
  for(Int_t ii=0; ii<128; ii++) fdat >> par_temp[0][ii] >> par_temp[1][ii];
  fdat.close();

  char path_input[200];                   
  sprintf(path_input,"%s",path_merged);          


  char Name_finput[200];
  sprintf(Name_finput,"%s/Run%dMS.root",path_input, Run_Number);


  TChain *fChain= new TChain("Tree");   
  fChain->Add(Name_finput);   
  fChain->SetMakeClass(1);              

  fChain->SetBranchAddress("ADC_High_SFT",adc_high_sft);
  fChain->SetBranchAddress("ADC_Low_SFT",adc_low_sft);
  fChain->SetBranchAddress("TDC_LE_SFT",tdc_le_sft);
  fChain->SetBranchAddress("TDC_TE_SFT",tdc_te_sft);


  fChain->GetEntry(evt);  

  for (Int_t j_SFT=0; j_SFT<128; j_SFT++){
    ADC_High_SFT[j_SFT]=adc_high_sft[j_SFT]-HG_SFT_ADC_Thr[j_SFT];
    TDC_LE_SFT[j_SFT]=tdc_le_sft[j_SFT][0];
    ADC_Low_SFT[j_SFT]=adc_low_sft[j_SFT]-LG_SFT_ADC_Thr[j_SFT];
  }

  for(int j=0 ; j<128 ; j++){ 
    if(ADC_High_SFT[j]<0)     
      ADC_High_SFT_corr[j]=0; 
    else
      ADC_High_SFT_corr[j]=ADC_High_SFT[j]; 
  }  

  for(Int_t ii=0; ii<128; ii++){
    for (Int_t qq=0; qq<6; qq++) {
      if (tdc_le_sft[ii][qq] > SFT_TDC_min[ii] && tdc_le_sft[ii][qq] < SFT_TDC_max[ii]) has_TDC_SFT_hit[ii]++;
    }
  }



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

  // if(*bad_sft_flag)
  //   return;
  // else
  //   cout << "*** Good SFT" << endl << endl;
     

  // Test if more than 1 layer is empty or if any layer has more than 3 p fibers
  if((layer_1_p_fibers_d.size() == 0 && layer_2_p_fibers_d.size() == 0) ||  
     (layer_3_p_fibers_d.size() == 0 && layer_4_p_fibers_d.size() == 0)){
    cout << "Need 1 hit in each helicity." << endl;
    return;
  }

  int flag = 0;
  if(layer_1_p_fibers_d.size() > sft_hit_count_max){
    cout << "More than " << sft_hit_count_max << " SFT hits in layer 1" << endl;
    flag = 1;
  }
  if(layer_2_p_fibers_d.size() > sft_hit_count_max){
    cout << "More than " << sft_hit_count_max << " SFT hits in layer 2" << endl;
    flag = 1;
  }
  if(layer_3_p_fibers_d.size() > sft_hit_count_max){
    cout << "More than " << sft_hit_count_max << " SFT hits in layer 3" << endl;
    flag = 1;
  }
  if(layer_4_p_fibers_d.size() > sft_hit_count_max){
    cout << "More than " << sft_hit_count_max << " SFT hits in layer 4" << endl;
    flag = 1;
  }  

  if(flag==1)
    return; 


  int iteration_counter = 0;

  // Storage arrays
  double z_avgs_arr[100];
  double delta_z_arr[100] = {-1};
  double p1_store_arr[100] = {-1};
  double p2_store_arr[100] = {-1};
  double p3_store_arr[100] = {-1};
  double p4_store_arr[100] = {-1};

  int N1_store_arr[100] = {-1};
  int N2_store_arr[100] = {-1};
  int N3_store_arr[100] = {-1};
  int N4_store_arr[100] = {-1};

  double z1_store_arr[100] = {-1000};
  double z2_store_arr[100] = {-1000};
  double z3_store_arr[100] = {-1000};
  double z4_store_arr[100] = {-1000};  

  int p1_p2_flag[100] = {1};
  int p3_p4_flag[100] = {1};

  int z1_z2_flag[100] = {1};
  int z3_z4_flag[100] = {1};
  int z1234_flag[100] = {1};

  // Case 1: all layers have atleast one hit
  if(layer_1_p_fibers_d.size() != 0 && layer_2_p_fibers_d.size() != 0 && layer_3_p_fibers_d.size() != 0 && layer_4_p_fibers_d.size() != 0){
    for(vector<double>::iterator it_1 = layer_1_p_fibers_d.begin(); it_1 != layer_1_p_fibers_d.end(); it_1++){
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

            return_array = calculate_min_Z_position(phi, *it_1, *it_2, *it_3, *it_4, array_to_return);

            z_avgs_arr[iteration_counter] = *return_array;
            delta_z_arr[iteration_counter] = *(return_array + 1);
            N1_store_arr[iteration_counter] = *(return_array + 2);
            N2_store_arr[iteration_counter] = *(return_array + 3);
            N3_store_arr[iteration_counter] = *(return_array + 4);
            N4_store_arr[iteration_counter] = *(return_array + 5);
            p1_store_arr[iteration_counter] = *(return_array + 6);
            p2_store_arr[iteration_counter] = *(return_array + 7);
            p3_store_arr[iteration_counter] = *(return_array + 8);
            p4_store_arr[iteration_counter] = *(return_array + 9);
            z1_store_arr[iteration_counter] = *(return_array + 10);
            z2_store_arr[iteration_counter] = *(return_array + 11);
            z3_store_arr[iteration_counter] = *(return_array + 12);
            z4_store_arr[iteration_counter] = *(return_array + 13);
            p1_p2_flag[iteration_counter] = *(return_array + 14);
            p3_p4_flag[iteration_counter] = *(return_array + 15);
            z1_z2_flag[iteration_counter] = *(return_array + 16);
            z3_z4_flag[iteration_counter] = *(return_array + 17);
            z1234_flag[iteration_counter] = *(return_array + 18);

            iteration_counter++;
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

          return_array = calculate_min_Z_position(phi, -1, *it_2, *it_3, *it_4, array_to_return);

          z_avgs_arr[iteration_counter] = *return_array;
          delta_z_arr[iteration_counter] = *(return_array + 1);
          N1_store_arr[iteration_counter] = *(return_array + 2);
          N2_store_arr[iteration_counter] = *(return_array + 3);
          N3_store_arr[iteration_counter] = *(return_array + 4);
          N4_store_arr[iteration_counter] = *(return_array + 5);
          p1_store_arr[iteration_counter] = *(return_array + 6);
          p2_store_arr[iteration_counter] = *(return_array + 7);
          p3_store_arr[iteration_counter] = *(return_array + 8);
          p4_store_arr[iteration_counter] = *(return_array + 9);
          z1_store_arr[iteration_counter] = *(return_array + 10);
          z2_store_arr[iteration_counter] = *(return_array + 11);
          z3_store_arr[iteration_counter] = *(return_array + 12);
          z4_store_arr[iteration_counter] = *(return_array + 13);
          p1_p2_flag[iteration_counter] = *(return_array + 14);
          p3_p4_flag[iteration_counter] = *(return_array + 15);
          z1_z2_flag[iteration_counter] = *(return_array + 16);
          z3_z4_flag[iteration_counter] = *(return_array + 17);
          z1234_flag[iteration_counter] = *(return_array + 18);

          iteration_counter++;
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

          return_array = calculate_min_Z_position(phi, *it_1, -1, *it_3, *it_4, array_to_return);

          z_avgs_arr[iteration_counter] = *return_array;
          delta_z_arr[iteration_counter] = *(return_array + 1);
          N1_store_arr[iteration_counter] = *(return_array + 2);
          N2_store_arr[iteration_counter] = *(return_array + 3);
          N3_store_arr[iteration_counter] = *(return_array + 4);
          N4_store_arr[iteration_counter] = *(return_array + 5);
          p1_store_arr[iteration_counter] = *(return_array + 6);
          p2_store_arr[iteration_counter] = *(return_array + 7);
          p3_store_arr[iteration_counter] = *(return_array + 8);
          p4_store_arr[iteration_counter] = *(return_array + 9);
          z1_store_arr[iteration_counter] = *(return_array + 10);
          z2_store_arr[iteration_counter] = *(return_array + 11);
          z3_store_arr[iteration_counter] = *(return_array + 12);
          z4_store_arr[iteration_counter] = *(return_array + 13);
          p1_p2_flag[iteration_counter] = *(return_array + 14);
          p3_p4_flag[iteration_counter] = *(return_array + 15);
          z1_z2_flag[iteration_counter] = *(return_array + 16);
          z3_z4_flag[iteration_counter] = *(return_array + 17);
          z1234_flag[iteration_counter] = *(return_array + 18);

          iteration_counter++;
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

          return_array = calculate_min_Z_position(phi, *it_1, *it_2, -1, *it_4, array_to_return);

          z_avgs_arr[iteration_counter] = *return_array;
          delta_z_arr[iteration_counter] = *(return_array + 1);
          N1_store_arr[iteration_counter] = *(return_array + 2);
          N2_store_arr[iteration_counter] = *(return_array + 3);
          N3_store_arr[iteration_counter] = *(return_array + 4);
          N4_store_arr[iteration_counter] = *(return_array + 5);
          p1_store_arr[iteration_counter] = *(return_array + 6);
          p2_store_arr[iteration_counter] = *(return_array + 7);
          p3_store_arr[iteration_counter] = *(return_array + 8);
          p4_store_arr[iteration_counter] = *(return_array + 9);
          z1_store_arr[iteration_counter] = *(return_array + 10);
          z2_store_arr[iteration_counter] = *(return_array + 11);
          z3_store_arr[iteration_counter] = *(return_array + 12);
          z4_store_arr[iteration_counter] = *(return_array + 13);
          p1_p2_flag[iteration_counter] = *(return_array + 14);
          p3_p4_flag[iteration_counter] = *(return_array + 15);
          z1_z2_flag[iteration_counter] = *(return_array + 16);
          z3_z4_flag[iteration_counter] = *(return_array + 17);
          z1234_flag[iteration_counter] = *(return_array + 18);

          iteration_counter++;
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

          return_array = calculate_min_Z_position(phi, *it_1, *it_2, *it_3, -1, array_to_return);

          z_avgs_arr[iteration_counter] = *return_array;
          delta_z_arr[iteration_counter] = *(return_array + 1);
          N1_store_arr[iteration_counter] = *(return_array + 2);
          N2_store_arr[iteration_counter] = *(return_array + 3);
          N3_store_arr[iteration_counter] = *(return_array + 4);
          N4_store_arr[iteration_counter] = *(return_array + 5);
          p1_store_arr[iteration_counter] = *(return_array + 6);
          p2_store_arr[iteration_counter] = *(return_array + 7);
          p3_store_arr[iteration_counter] = *(return_array + 8);
          p4_store_arr[iteration_counter] = *(return_array + 9);
          z1_store_arr[iteration_counter] = *(return_array + 10);
          z2_store_arr[iteration_counter] = *(return_array + 11);
          z3_store_arr[iteration_counter] = *(return_array + 12);
          z4_store_arr[iteration_counter] = *(return_array + 13);
          p1_p2_flag[iteration_counter] = *(return_array + 14);
          p3_p4_flag[iteration_counter] = *(return_array + 15);
          z1_z2_flag[iteration_counter] = *(return_array + 16);
          z3_z4_flag[iteration_counter] = *(return_array + 17);
          z1234_flag[iteration_counter] = *(return_array + 18);

          iteration_counter++;
        }
      }
    }
  }  

  // Cases: 1 hit in each helicity
  //        1) 0 1 0 1
  //        2) 0 1 1 0
  //        3) 1 0 0 1
  //        4) 1 0 1 0


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

        return_array = calculate_min_Z_position(phi, -1, *it_2, -1, *it_4, array_to_return);

        z_avgs_arr[iteration_counter] = *return_array;
        delta_z_arr[iteration_counter] = *(return_array + 1);
        N1_store_arr[iteration_counter] = *(return_array + 2);
        N2_store_arr[iteration_counter] = *(return_array + 3);
        N3_store_arr[iteration_counter] = *(return_array + 4);
        N4_store_arr[iteration_counter] = *(return_array + 5);
        p1_store_arr[iteration_counter] = *(return_array + 6);
        p2_store_arr[iteration_counter] = *(return_array + 7);
        p3_store_arr[iteration_counter] = *(return_array + 8);
        p4_store_arr[iteration_counter] = *(return_array + 9);
        z1_store_arr[iteration_counter] = *(return_array + 10);
        z2_store_arr[iteration_counter] = *(return_array + 11);
        z3_store_arr[iteration_counter] = *(return_array + 12);
        z4_store_arr[iteration_counter] = *(return_array + 13);
        p1_p2_flag[iteration_counter] = *(return_array + 14);
        p3_p4_flag[iteration_counter] = *(return_array + 15);
        z1_z2_flag[iteration_counter] = *(return_array + 16);
        z3_z4_flag[iteration_counter] = *(return_array + 17);
        z1234_flag[iteration_counter] = *(return_array + 18);

        iteration_counter++;
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

        return_array = calculate_min_Z_position(phi, -1, *it_2, *it_3, -1, array_to_return);

        z_avgs_arr[iteration_counter] = *return_array;
        delta_z_arr[iteration_counter] = *(return_array + 1);
        N1_store_arr[iteration_counter] = *(return_array + 2);
        N2_store_arr[iteration_counter] = *(return_array + 3);
        N3_store_arr[iteration_counter] = *(return_array + 4);
        N4_store_arr[iteration_counter] = *(return_array + 5);
        p1_store_arr[iteration_counter] = *(return_array + 6);
        p2_store_arr[iteration_counter] = *(return_array + 7);
        p3_store_arr[iteration_counter] = *(return_array + 8);
        p4_store_arr[iteration_counter] = *(return_array + 9);
        z1_store_arr[iteration_counter] = *(return_array + 10);
        z2_store_arr[iteration_counter] = *(return_array + 11);
        z3_store_arr[iteration_counter] = *(return_array + 12);
        z4_store_arr[iteration_counter] = *(return_array + 13);
        p1_p2_flag[iteration_counter] = *(return_array + 14);
        p3_p4_flag[iteration_counter] = *(return_array + 15);
        z1_z2_flag[iteration_counter] = *(return_array + 16);
        z3_z4_flag[iteration_counter] = *(return_array + 17);
        z1234_flag[iteration_counter] = *(return_array + 18);

        iteration_counter++;
      }
    }
  }




  if(layer_1_p_fibers_d.size() != 0 && layer_2_p_fibers_d.size() == 0 && layer_3_p_fibers_d.size() == 0 && layer_4_p_fibers_d.size() != 0){
    for(vector<double>::iterator it_1 = layer_1_p_fibers_d.begin(); it_1 != layer_1_p_fibers_d.end(); it_1++){
      for(vector<double>::iterator it_4 = layer_4_p_fibers_d.begin(); it_4 != layer_4_p_fibers_d.end(); it_4++){ // 

        double array_to_return[20] = {1};
        double *return_array;
        // flag_array[14] == p1-p2 flag
        // flag_array[15] == p3-p4 flag
        // flag_array[16] == z12 out of threshold
        // flag_array[17] == z34 out of thresholf
        // flag_array[18] == z1234 out of threshold

        return_array = calculate_min_Z_position(phi, *it_1, -1, -1, *it_4, array_to_return);

        z_avgs_arr[iteration_counter] = *return_array;
        delta_z_arr[iteration_counter] = *(return_array + 1);
        N1_store_arr[iteration_counter] = *(return_array + 2);
        N2_store_arr[iteration_counter] = *(return_array + 3);
        N3_store_arr[iteration_counter] = *(return_array + 4);
        N4_store_arr[iteration_counter] = *(return_array + 5);
        p1_store_arr[iteration_counter] = *(return_array + 6);
        p2_store_arr[iteration_counter] = *(return_array + 7);
        p3_store_arr[iteration_counter] = *(return_array + 8);
        p4_store_arr[iteration_counter] = *(return_array + 9);
        z1_store_arr[iteration_counter] = *(return_array + 10);
        z2_store_arr[iteration_counter] = *(return_array + 11);
        z3_store_arr[iteration_counter] = *(return_array + 12);
        z4_store_arr[iteration_counter] = *(return_array + 13);
        p1_p2_flag[iteration_counter] = *(return_array + 14);
        p3_p4_flag[iteration_counter] = *(return_array + 15);
        z1_z2_flag[iteration_counter] = *(return_array + 16);
        z3_z4_flag[iteration_counter] = *(return_array + 17);
        z1234_flag[iteration_counter] = *(return_array + 18);

        iteration_counter++;
      }
    }
  }


  if(layer_1_p_fibers_d.size() != 0 && layer_2_p_fibers_d.size() == 0 && layer_3_p_fibers_d.size() != 0 && layer_4_p_fibers_d.size() == 0){
    for(vector<double>::iterator it_1 = layer_1_p_fibers_d.begin(); it_1 != layer_1_p_fibers_d.end(); it_1++){
      for(vector<double>::iterator it_3 = layer_3_p_fibers_d.begin(); it_3 != layer_3_p_fibers_d.end(); it_3++){
      
        double array_to_return[20] = {1};
        double *return_array;
        // flag_array[14] == p1-p2 flag
        // flag_array[15] == p3-p4 flag
        // flag_array[16] == z12 out of threshold
        // flag_array[17] == z34 out of thresholf
        // flag_array[18] == z1234 out of threshold

        return_array = calculate_min_Z_position(phi, *it_1, -1, *it_3, -1, array_to_return);

        z_avgs_arr[iteration_counter] = *return_array;
        delta_z_arr[iteration_counter] = *(return_array + 1);
        N1_store_arr[iteration_counter] = *(return_array + 2);
        N2_store_arr[iteration_counter] = *(return_array + 3);
        N3_store_arr[iteration_counter] = *(return_array + 4);
        N4_store_arr[iteration_counter] = *(return_array + 5);
        p1_store_arr[iteration_counter] = *(return_array + 6);
        p2_store_arr[iteration_counter] = *(return_array + 7);
        p3_store_arr[iteration_counter] = *(return_array + 8);
        p4_store_arr[iteration_counter] = *(return_array + 9);
        z1_store_arr[iteration_counter] = *(return_array + 10);
        z2_store_arr[iteration_counter] = *(return_array + 11);
        z3_store_arr[iteration_counter] = *(return_array + 12);
        z4_store_arr[iteration_counter] = *(return_array + 13);
        p1_p2_flag[iteration_counter] = *(return_array + 14);
        p3_p4_flag[iteration_counter] = *(return_array + 15);
        z1_z2_flag[iteration_counter] = *(return_array + 16);
        z3_z4_flag[iteration_counter] = *(return_array + 17);
        z1234_flag[iteration_counter] = *(return_array + 18);

        iteration_counter++;
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

  for(int i = 0; i < iteration_counter; i++){
    if(delta_z_arr[i] < delta_z_thr){
      good_delta_z_index.push_back(i);
    }
  }

  cout << "SFT Z Values -- -- Angle Used = " << phi << endl << endl;

 
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

    if(p1_p2_flag[*it] > 1 && p3_p4_flag[*it] > 1){
      if(z_avgs_store.empty()){
        z_printed = false;
      }
      else{
        for(vector<double>::iterator db = z_avgs_store.begin(); db != z_avgs_store.end(); db++){
          if(fabs(*db - z_avg/double(z_avg_count)) < 1){
            z_printed = true;
          }
        }
      }

      if(1/*!z_printed*/){
        cout << "Layer    Turn Number          P               Z" << endl;
        printf("  1          %2d              %3.1f           %3.2f    \n" ,
        N1_store_arr[*it],p1_store_arr[*it],z1_store_arr[*it]);
        printf("  2          %2d              %3.1f           %3.2f    \n" ,
        N2_store_arr[*it],p2_store_arr[*it],z2_store_arr[*it]);
        printf("  3          %2d              %3.1f           %3.2f    \n" ,
        N3_store_arr[*it],p3_store_arr[*it],z3_store_arr[*it]);
        printf("  4          %2d              %3.1f           %3.2f    \n" ,
        N4_store_arr[*it],p4_store_arr[*it],z4_store_arr[*it]);

        cout << "Z Avg = " << z_avg/double(z_avg_count) << endl;
        cout << "Delta Z = " << delta_z_arr[*it] << endl << endl;
      }
    }
   

    z_avgs_store.push_back(z_avg/double(z_avg_count));
  }

  cout << endl;
}
