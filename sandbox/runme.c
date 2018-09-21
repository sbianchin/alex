void runme(int Run_Number=3994){
	char Name[100];
	sprintf(Name, "Batch_3_0(%d,0,10,2)", Run_Number);
	gROOT->ProcessLine(".L Batch_3.0.C+");
	gROOT->ProcessLine(Name);
	gSystem->Exit(0);
}
	

