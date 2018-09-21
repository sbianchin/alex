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


void _NEED_TO_ROTATE(vector<int> vec_myVector){

	vector <double> vec_xx_vector;	vec_xx_vector.clear();
	vector <double> vec_yy_vector;	vec_yy_vector.clear();

	double sum = 0.;
	double mean = 0.;
	int counter = 0;

	for(unsigned int i=0; i<vec_myVector.size(); i++){
		vec_xx_vector.push_back(Xloc[vec_myVector[i]]);
		vec_yy_vector.push_back(Yloc[vec_myVector[i]]);
		cout << "TEST : " << vec_xx_vector[i] << "  " << vec_yy_vector[i] << endl;
	} 

	for(unsigned int j=0; j<vec_myVector.size(); j++){
		for(unsigned int k=0; k<vec_myVector.size(); k++){
			if(j<k){
				sum += abs(vec_xx_vector[j] - vec_xx_vector[k]);
				counter++;
				cout << "ALBANE : " << j << "  " << abs(vec_xx_vector[j] - vec_xx_vector[k]) << endl;
			}
		}
		cout << endl;
	}

	mean = sum/counter;

	cout << endl;
	cout << "SUM : " << sum << endl;
	cout << "Counter : " << counter << endl;
	cout << "Mean : " << sum/counter << "  " << mean << endl;


	return;
}

void _ADC_HG_SUM(vector<int> vec_myVector){

	int ADC_High_TARGET[256];

	for(unsigned int i=0; i<vec_myVector.size(); i++){
		cout << "RUQUIER : " << ADC_High_TARGET[vec_myVector[i]] << endl;
		//cout << "RUQUIER : " << vec_myVector[i] << endl;
	}

	cout << endl;

	return;
}
