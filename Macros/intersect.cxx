#ifndef __CINT__
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <TROOT.h>
#include "TSystem.h"
#include "TFile.h"
#include "TProfile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TRandom.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TPolyLine.h"
#include "TLine.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TBranch.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TEllipse.h"
#include "TMarker.h"
#endif

#include <vector>

float intersectx1(float a=0., float b=0., float R=44.5)
{
	float x1 = 0.;

	x1 = (-2*a*b - sqrt(4*a*a*b*b - 4*(a*a+1)*(b*b-R*R))) / (2*(a*a+1));


	return x1;
}

float intersectx2(float a=0., float b=0., float R=44.5)
{
	float x2 = 0.;

	x2 = (-2*a*b + sqrt(4*a*a*b*b - 4*(a*a+1)*(b*b-R*R))) / (2*(a*a+1));


	return x2;
}

float y1_int(float x1=0., float a=0., float b=0.)
{
	float y1 = 0.;
	y1 = a * x1 + b;

	return y1;
}

float y2_int(float x2=0., float a=0., float b=0.)
{
	float y2 = 0.;
	y2 = a * x2 + b;

	return y2;
}


float y_int2(float x=0.,float a=0., float b=0.)
{
	float y = 0.;

	y = a * x + b;

	return y;
}

float distance(float x1=0., float y1=0., float x2=0., float y2=0.)
{
	float dist = 0;
	dist = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));

	return dist;

}


float distance_to_line(float x=0., float y=0., float a=0., float b=0.) // Line is given by equation y=a*x + b;
{
	float dist = 0.;
	  
	if((-a*x+y-b) < 0) dist = -(-a*x+y-b)/sqrt(a*a+1);
	else dist = (-a*x+y-b)/sqrt(a*a+1);
	
	return dist;
}	


double points_to_slope(double x1 = 0.0, double y1 = 0.0, double x2 = 1.0, double y2 = 1.0){
	double slope = 0.0;

	slope = (y2-y1)/(x2-x1);

	return slope;
}
	
double points_to_y_int(double x1 = 0.0, double y1 = 0.0, double x2 = 1.0, double y2 = 1.0){
	double y_int = 0;

	y_int = y1 - x1*(y2-y1)/(x2-x1);

	return y_int;
}


double C2_intersect(double a = 0.0, double b = 0.0, double R = 40.0){
	return (R-b)/a;
}
	


bool IsIn(int a=0, int l0=0, int l1=0, int l2=0, int l3=0, int l4=0, int l5=0, int l6=0, int l7=0)
{

	bool output = false;

		if(a==l0 || a==l1 || a==l2 || a==l3 || a==l4 || a==l5 || a==l6 || a==l7){
			output = true;
		}

	return output;
}


double rotate_by_angle(double x=0, double y=0, int theta=0, char coord='x'){

	double new_coordinate = 0;

	if(!(theta==0 || theta==30 || theta==60 || theta==90 || theta==120 || theta==150 || theta==180 ||
		theta==-30 || theta==-60 || theta==-90 || theta==-120 || theta==-150)){
		return x;
	}


	if(coord=='x'){

		new_coordinate = cos(theta*M_PI/180.0)*x - sin(theta*M_PI/180.0)*y;

	}
	else if(coord=='y'){

		new_coordinate = sin(theta*M_PI/180.0)*x + cos(theta*M_PI/180.0)*y;

	}

	return new_coordinate;
}

std::vector<double> _2lines_intersect(double a1=99., double b1=99., double a2=99., double b2=99.){

  std::vector<double> vec_2lines_intersect;

	vec_2lines_intersect.push_back((b2 - b1) / (a1 - a2));
	vec_2lines_intersect.push_back((a1*b2 - a2*b1) / (a1 - a2));

	return vec_2lines_intersect;
}

//double _2lines_intersect(double a1=99., double b1=99., double a2=99., double b2=99.){
//
//	double xc = 99;
//	double yc = 99;
//	double _2lines_int_coordinates[2] = 99.;
//
//	_2lines_int_coordinates[0] = (b2 - b1) / (a1 - a2);
//	_2lines_int_coordinates[1] = (a1*b2 - a2*b1) / (a1 - a2);
//
//	return
//}
