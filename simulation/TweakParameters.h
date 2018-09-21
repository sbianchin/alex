
#ifndef TWEAK_PARAMETERS_H
#define TWEAK_PARAMETERS_H

const Double_t TWO_TRACK_MIN_CHISQ = 1;
const Double_t TWO_TRACK_MIN_LEPTON_BARS = 5;

const Double_t FIT_TOF1_WEIGHT = 3;

const Double_t K_STOP_CENTROID_THRESH = 0.0; // mm

const Bool_t PATH_TRAVERSIAL_USE_ALL = true; // use FindShortestPath() instead of dijkstra()
const Double_t PATH_TRAVERSIAL_DIJKSTRA_JUMP_RADIUS = 9; // mm
const Double_t PATH_TRAVERSIAL_ALL_PENALTY = 1000;

#endif
                            