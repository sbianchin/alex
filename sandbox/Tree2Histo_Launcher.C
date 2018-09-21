void Tree2Histo_Launcher(Int_t Run_Number=1800)
{
	char Name[100];
	sprintf(Name, "Tree2Histo(%d)", Run_Number);
	gROOT->ProcessLine(".L Tree2Histo.C+");
	gROOT->ProcessLine(Name);
	gROOT->ProcessLine(".q");
}