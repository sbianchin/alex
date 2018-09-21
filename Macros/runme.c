void runme(int Run_Number=3994){
	char Name[100];
	//sprintf(Name, "Batch_5_5(%d,0,10,2)", Run_Number);
	sprintf(Name, "Batch_5_5(%d,0,200,1)", Run_Number);
	gROOT->ProcessLine(".L Batch_5.5.C+");
	gROOT->ProcessLine(Name);
	gSystem->Exit(0);
}
	

