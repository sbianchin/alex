void calculate_Z_position(double phi = 90., int N = 7 , double p = 4.0, int layer = 4){
  double Z_position = 0.0;

  double PI = 3.1415;
  
  phi *= PI/180.0;

  // Z_0 CONSTANT OFFSET FOR EACH LAYER
  /////////////////////////////////////////////////////////////
  const double z_offset_const[4] = {-131.0, -136.5, -138.68, -135.83};

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
               + 0.5*fiber_diameter;///cos(alpha[1]);     
  }
  else if(layer == 3){
    Z_position = z_offset_const[2] + layer_radius[2]*(2*PI-phi)*tan(alpha[2]) 
               + (N-1)*ribbon_width[2]
               + (fibers_per_ribbon[2]-p)*ribbon_width[2]/fibers_per_ribbon[2]
               + 0.5*fiber_diameter;///cos(alpha[2]);   
  }
  else if(layer == 4){
    Z_position = z_offset_const[3] + layer_radius[3]*(2*PI-phi)*tan(alpha[3]) 
               + (N-1)*ribbon_width[3]
               + (fibers_per_ribbon[3]-p)*ribbon_width[3]/fibers_per_ribbon[3]
               + 0.5*fiber_diameter/cos(alpha[3]);   
  }

  if(p <= 0)
    Z_position = -1;


  double C4_const =  (layer_radius[3]*(2*PI-PI/2.0)*tan(alpha[3]) 
                      + (7-1)*ribbon_width[3]
                      + (fibers_per_ribbon[3]-4)*ribbon_width[3]/fibers_per_ribbon[3]
                     + 0.5*fiber_diameter/cos(alpha[3])); 

  cout << "C4 constant = "  << C4_const << endl;


  return;
}