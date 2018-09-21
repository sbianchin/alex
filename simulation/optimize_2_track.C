#include <stdio.h>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <TROOT.h>
#include <typeinfo>
#include <dlib/optimization.h>
#include "GlobalParameters.h"
#include "TargetParameters.h"
typedef dlib::matrix<double,4,1> two_track_parameter_vector;

two_track_parameter_vector derivative_calculate_2_track_residual(const vector<double> data,
  const two_track_parameter_vector& params)
  {
  double m1 = params(0); double b1 = params(1); double m2 = params(2); double b2 = params(3);
  double x = data[0]; double y = data[1];

  // if (x > 50) {
  //   x = x - 100; y = y - 100;
  //   double ssr_1 = abs(m1*x + b1 - y);
  //   double ssr_2 = abs(m2*x + b2 - y);
  //   if (ssr_1 <= ssr_2) {
  //     return ssr_1 * 0.25;
  //   }
  //   else return ssr_2 * 0.25;
  //
  // }
  // cout << "X, Y of Point" << vec_xx_lepton[point] << ", " << vec_yy_lepton[point] << ", ";
  // cout << "Residual 1, 2: " << m1*vec_xx_lepton[point] + b1 - vec_yy_lepton[point] << ", ";
  // cout << m2*vec_xx_lepton[point] + b2 - vec_yy_lepton[point] << ", ";
  double ssr_1 = abs(m1*x + b1 - y); //first SSR
  double ssr_2 = abs(m2*x + b2 - y); //first SSR
  // cout << endl << "ssr 1, 2: " << ssr_1 << ", " << ssr_2 << endl;
  // double ssr_1 = distance_to_line(x, y, m1, b1); //first SSR
  // double ssr_2 = distance_to_line(x, y, m2, b2); //first SSR

  // cout << "ssr1, ssr2: ";
  // cout << ssr_1 << ", " << ssr_2 << endl;


  // double square_ssr_1 = pow(ssr_1, 2);
  // double square_ssr_2 = pow(ssr_2, 2);
  two_track_parameter_vector derivative;
  if (ssr_1 < ssr_2) {
    double dm = x;
    double db = 1;
    derivative = dm, db, 0, 0;
  }
  else {
    double dm = x;
    double db = 1;
    derivative = 0, 0, dm, db;
  }
  return derivative;


}

double calculate_2_track_residual(const vector<double> data,
  const two_track_parameter_vector& params)
  {
  double m1 = params(0); double b1 = params(1); double m2 = params(2); double b2 = params(3);
  double x = data[0]; double y = data[1];

  // if (x > 50) {
  //   x = x - 100; y = y - 100;
  //   double ssr_1 = abs(m1*x + b1 - y);
  //   double ssr_2 = abs(m2*x + b2 - y);
  //   if (ssr_1 <= ssr_2) {
  //     return ssr_1 * 0.25;
  //   }
  //   else return ssr_2 * 0.25;
  //
  // }
  // cout << "X, Y of Point" << vec_xx_lepton[point] << ", " << vec_yy_lepton[point] << ", ";
  // cout << "Residual 1, 2: " << m1*vec_xx_lepton[point] + b1 - vec_yy_lepton[point] << ", ";
  // cout << m2*vec_xx_lepton[point] + b2 - vec_yy_lepton[point] << ", ";
  double ssr_1 = abs(m1*x + b1 - y); //first SSR
  double ssr_2 = abs(m2*x + b2 - y); //first SSR
  // cout << endl << "ssr 1, 2: " << ssr_1 << ", " << ssr_2 << endl;
  // double ssr_1 = distance_to_line(x, y, m1, b1); //first SSR
  // double ssr_2 = distance_to_line(x, y, m2, b2); //first SSR

  // cout << "ssr1, ssr2: ";
  // cout << ssr_1 << ", " << ssr_2 << endl;
  // vector<double> intersect = _2lines_intersect(m1, b1, m2, b2);
  // double dist = distance(intersect[0], intersect[1], x, y);
  if (ssr_1 <= ssr_2) {
    // if (!isnan(dist)) {
        return pow(ssr_1, 1);// + dist*0.15;
    // }
    // else return pow(ssr_1, 1);
  }
  else {
    // if (!isnan(dist)) {
        return pow(ssr_2, 1);// + dist*0.15;
    // }
    // else return pow(ssr_2, 1);
  }


}

int optimize_2_track() {
  vector<vector<double>> lepton_coordinates;
  vector<int> lepton_bars = {16, 17, 18, 19, 20, 21, 22, 23, 24,
  32, 48, 65, 83, 102, 121, 140};
  vector<double> vec_xx_lepton, vec_yy_lepton;
  for (int i = 0; i < lepton_bars.size(); ++i)
    {
        // save the pair
        vector<double> input;
        input.push_back(Xloc[lepton_bars[i]]);
        input.push_back(Yloc[lepton_bars[i]]);
        lepton_coordinates.push_back(input);
    }
  two_track_parameter_vector params, params2;
  // params = 1.1, 1.1, -1.1, -1.1;
  params = 3, 3, -3.1, -3.1;
  // params2 = 1.1, 1.1, -1.1, -1.1;
  params2 = 3, 3, -3.1, -3.1;

  cout << "derivative analytic: " << trans(derivative_calculate_2_track_residual(lepton_coordinates[0], params)) << endl;
  cout << "approximate derivative: " << trans(dlib::derivative(calculate_2_track_residual)(lepton_coordinates[0], params)) << endl;
  cout << "derivative error: " << length(derivative_calculate_2_track_residual(lepton_coordinates[0], params) -
                                               dlib::derivative(calculate_2_track_residual)(lepton_coordinates[0], params) ) << endl;


  dlib::solve_least_squares_lm(dlib::objective_delta_stop_strategy(1e-7).be_verbose(),
                             calculate_2_track_residual,
                             derivative_calculate_2_track_residual,
                             lepton_coordinates,
                             params);
   cout << "inferred parameters: "<< trans(params) << endl;
  //  cout << "solution error:      "<< length(params - params) << endl;

   dlib::solve_least_squares_lm(dlib::objective_delta_stop_strategy(1e-7).be_verbose(),
                              calculate_2_track_residual,
                              dlib::derivative(calculate_2_track_residual),
                              lepton_coordinates,
                              params2);
    cout << "inferred parameters: "<< trans(params2) << endl;
    // cout << "solution error:      "<< length(params - params) << endl;
  //  cout << endl;
   return 1;

}
