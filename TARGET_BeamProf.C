void TARGET_BeamProf(Int_t Run_Number=1592)
{
  		
  char beamprofile[100];

  sprintf(beamprofile,"TARGET_BeamProfile(%d)",Run_Number);
	
  gROOT->ProcessLine(".L TARGET_BeamProfile.C+");

  gROOT->ProcessLine(beamprofile);


}