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
//#include "TSpectrum.h"
#include "TMarker.h"
//#include "Event_Display_MS.h"
#include "ANAPATH.h"
#include "Thresholds.h"
//#include "ListGoodEvents.h"
#endif

void Z_Position_Simulator(float phi = 45, float p1 = 0, float p2 = 0, float p3 = 0, float p4 = 0){

	cout << "" << endl;

	float thr = 2.;

	cout << "p1: " << p1 << endl;
	cout << "p2: " << p2 << endl;
	cout << "p3: " << p3 << endl;
	cout << "p4: " << p4 << endl;

	float pi = 3.14195;

	//Constants
	float z1 = 141.;
	float z2 = 135.;
	float z3 = 157.2;
	float z4 = 168.;

	float R1 = 40.25;
	float R2 = 41.1;
	float R3 = 42.1;
	float R4 = 42.95;

	float alpha1 = 0.0623082017;
	float alpha2 = 0.0610864722;
	float alpha3 = 0.067544185;
	float alpha4 = 0.0661479228;

	float phi1 = 165.;
	float phi2 = 15.;
	float phi3 = 75.;
	float phi4 = 285.;

	float Zpos1[18] = {0};
	float Zpos2[18] = {0};
	float Zpos3[18] = {0};
	float Zpos4[18] = {0};

	for (int n=0; n<18; n++) {

		Zpos1[n] = z1 + (phi)*(pi/180)*R1*tan(alpha1) + (n-1)*15.7 + (p1-1)*1.05 + 0.5;
		Zpos2[n] = z2 + (phi)*(pi/180)*R2*tan(alpha2) + (n-1)*15.7 + (p2-1)*1.05 + 0.5;

		if (n<16) Zpos3[n] = z3 + (360-phi)*(pi/180)*(R3*tan(alpha3)) + (n-1)*17.8 + (17-p3)*1.05 + 0.5;
		if (n<14) Zpos4[n] = z4 + (360-phi)*(pi/180)*(R4*tan(alpha4)) + (n-1)*17.8 + (17-p4)*1.05 + 0.5;

		
		//Zpos1[n] = z1 + (phi - phi1)*(pi/180)*(R1*tan(alpha1)) + (n-1)*15.7 + (15-p1)*1.05 + 0.5;
		//Zpos2[n] = z2 + (phi - phi2)*(pi/180)*(R2*tan(alpha2)) + (n-1)*15.7 + (15-p2)*1.05 + 0.5;
		//if (n<16) Zpos3[n] = z3 + (360 - phi + phi3)*(pi/180)*(R3*tan(alpha3)) + (n-1)*17.8 + (p3-1)*1.05 + 0.5;
		//if (n<14) Zpos4[n] = z4 + (360 - phi + phi4)*(pi/180)*(R4*tan(alpha4)) + (n-1)*17.8 + (p4-1)*1.05 + 0.5;
	}

	 for (int n=0; n<18; n++) {
	 	cout << "" << endl;
	 	cout << "N = " << n << endl;
	 	cout << "Z1 with n = " << n << ": " << Zpos1[n] << endl;
	 	cout << "Z2 with n = " << n << ": " << Zpos2[n] << endl;
	 	cout << "Z3 with n = " << n << ": " << Zpos3[n] << endl;
	 	cout << "Z4 with n = " << n << ": " << Zpos4[n] << endl;
	 }

	cout << "" << endl;

	for (int i=0; i<18; i++) {
		for (int j=0; j<18; j++) {
			if (fabs(Zpos1[i]-Zpos2[j]) <= thr) {
				for (int k=0; k<16; k++) {
					if (fabs(Zpos1[i]-Zpos3[k]) <= thr) {
						for (int h=0; h<14; h++) {
							if (fabs(Zpos1[i]-Zpos4[h]) <= thr) {
								cout << "Z1 -- N = " << i << ": " << Zpos1[i] << " | Z2 -- N = " << j << ": " << Zpos2[j] << " | Z3 -- N = " << k << ": " 
								<< Zpos3[k] << " | Z4 -- N = " << h << ": " << Zpos4[h] << endl;
							}
						}
					}
				}
			}
		}
	}

	
} 