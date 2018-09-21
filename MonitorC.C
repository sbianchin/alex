void MonitorC(Int_t Run_Number=1800)
{
	char Name[100];
	sprintf(Name, "monitor(%d, 1)", Run_Number);
	gROOT->ProcessLine(".L monitor.cxx+");
	gROOT->ProcessLine(Name);
//	gROOT->ProcessLine(".q");
}
