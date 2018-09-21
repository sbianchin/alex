
void moniTarget(){
 	c1=new TCanvas("c1","",0,0,800,400);
	c1->Divide(2,1);
	gStyle->SetOptStat(0);
	target=new TH2F("target","target",20,0,20,20,0,20);
	target->SetXTitle("X-axis");
	target->SetYTitle("Y-axis");
	time=new TH2F("time","time",20, 0, 20, 20, 0, 20);
	time->SetXTitle("X-axis");
	time->SetYTitle("Y-axis");
      
        char output[100]="target.pdf";
        // char output[100]="target2.pdf";
   
	char name[100];
	double targetdE[256],targetdt[256],u_x,u_y,u_z,px[10],py[10],pz[10];
	int n;
	int targX[1256];
	int targY[1256];
	int m=0;
	/*
	for(int i=0;i<6;i++){  targX[m]=7+i; targY[m]=1; m++;}
	for(int i=0;i<10;i++){ targX[m]=5+i; targY[m]=2; m++;}
	for(int i=0;i<12;i++){ targX[m]=4+i; targY[m]=3; m++;}
	for(int i=0;i<14;i++){ targX[m]=3+i; targY[m]=4; m++;}
	for(int i=0;i<16;i++){ targX[m]=2+i; targY[m]=5; m++;}
	for(int i=0;i<16;i++){ targX[m]=2+i; targY[m]=6; m++;}
	for(int i=0;i<18;i++){ targX[m]=1+i; targY[m]=7; m++;}
	for(int i=0;i<18;i++){ targX[m]=1+i; targY[m]=8; m++;}
	for(int i=0;i<18;i++){ targX[m]=1+i; targY[m]=9; m++;} 
	for(int i=0;i<18;i++){ targX[m]=1+i; targY[m]=10; m++;}
	for(int i=0;i<18;i++){ targX[m]=1+i; targY[m]=11; m++;}
	for(int i=0;i<18;i++){ targX[m]=1+i; targY[m]=12; m++;}
	for(int i=0;i<16;i++){ targX[m]=2+i; targY[m]=13; m++;}
	for(int i=0;i<16;i++){ targX[m]=2+i; targY[m]=14; m++;}
	for(int i=0;i<14;i++){ targX[m]=3+i; targY[m]=15; m++;}
	for(int i=0;i<12;i++){ targX[m]=4+i; targY[m]=16; m++;}
	for(int i=0;i<10;i++){ targX[m]=5+i; targY[m]=17; m++;}
	for(int i=0;i<6;i++){  targX[m]=7+i; targY[m]=18; m++;}
	*/
        for(int i=0;i<6;i++){  targX[m]=7+i; targY[m]=18; m++;}
	for(int i=0;i<10;i++){ targX[m]=5+i; targY[m]=17; m++;}
	for(int i=0;i<12;i++){ targX[m]=4+i; targY[m]=16; m++;}
	for(int i=0;i<14;i++){ targX[m]=3+i; targY[m]=15; m++;}
	for(int i=0;i<16;i++){ targX[m]=2+i; targY[m]=14; m++;}
	for(int i=0;i<16;i++){ targX[m]=2+i; targY[m]=13; m++;}
	for(int i=0;i<18;i++){ targX[m]=1+i; targY[m]=12; m++;}
	for(int i=0;i<18;i++){ targX[m]=1+i; targY[m]=11; m++;}
	for(int i=0;i<18;i++){ targX[m]=1+i; targY[m]=10; m++;} 
	for(int i=0;i<18;i++){ targX[m]=1+i; targY[m]=9; m++;}
	for(int i=0;i<18;i++){ targX[m]=1+i; targY[m]=8; m++;}
	for(int i=0;i<18;i++){ targX[m]=1+i; targY[m]=7; m++;}
	for(int i=0;i<16;i++){ targX[m]=2+i; targY[m]=6; m++;}
	for(int i=0;i<16;i++){ targX[m]=2+i; targY[m]=5; m++;}
	for(int i=0;i<14;i++){ targX[m]=3+i; targY[m]=4; m++;}
	for(int i=0;i<12;i++){ targX[m]=4+i; targY[m]=3; m++;}
	for(int i=0;i<10;i++){ targX[m]=5+i; targY[m]=2; m++;}
	for(int i=0;i<6;i++){  targX[m]=7+i; targY[m]=1; m++;}
	/*
	for(int i=0;i<6;i++){  targX[m]=12-i; targY[m]=18; m++;}
	for(int i=0;i<10;i++){ targX[m]=14-i; targY[m]=17; m++;}
	for(int i=0;i<12;i++){ targX[m]=15-i; targY[m]=16; m++;}
	for(int i=0;i<14;i++){ targX[m]=16-i; targY[m]=15; m++;}
	for(int i=0;i<16;i++){ targX[m]=17-i; targY[m]=14; m++;}
	for(int i=0;i<16;i++){ targX[m]=17-i; targY[m]=13; m++;}
	for(int i=0;i<18;i++){ targX[m]=18-i; targY[m]=12; m++;}
	for(int i=0;i<18;i++){ targX[m]=18-i; targY[m]=11; m++;}
	for(int i=0;i<18;i++){ targX[m]=18-i; targY[m]=10; m++;} 
	for(int i=0;i<18;i++){ targX[m]=18-i; targY[m]=9; m++;}
	for(int i=0;i<18;i++){ targX[m]=18-i; targY[m]=8; m++;}
	for(int i=0;i<18;i++){ targX[m]=18-i; targY[m]=7; m++;}
	for(int i=0;i<16;i++){ targX[m]=17-i; targY[m]=6; m++;}
	for(int i=0;i<16;i++){ targX[m]=17-i; targY[m]=5; m++;}
	for(int i=0;i<14;i++){ targX[m]=16-i; targY[m]=4; m++;}
	for(int i=0;i<12;i++){ targX[m]=15-i; targY[m]=3; m++;}
	for(int i=0;i<10;i++){ targX[m]=14-i; targY[m]=2; m++;}
	for(int i=0;i<6;i++){  targX[m]=12-i; targY[m]=1; m++;}
        */
	file=new TFile("testing.root");
	t1=(TTree*)file->Get("Kaon");
	t1->SetBranchAddress("n",&n);
	//	t1->SetBranchAddress("itrack)", &itrack); //JIIIIIIIIIII
	t1->SetBranchAddress("targetdE[256]",&targetdE);
	t1->SetBranchAddress("targetdt[256]",&targetdt);
	t1->SetBranchAddress("u_x",&u_x);
	t1->SetBranchAddress("u_y",&u_y);	t1->SetBranchAddress("u_z",&u_z);
	t1->SetBranchAddress("px[10]",&px);
	t1->SetBranchAddress("py[10]",&py);
	t1->SetBranchAddress("pz[10]",&pz);

	sprintf(name,"%s[",output);
	c1->Print(name);

	for(int j=1;j<t1->GetEntries();j++){
	t1->GetEntry(j);
	for(int i=0;i<256;i++){
		target->Fill(targX[i],targY[i],targetdE[i]);
		time->Fill(targX[i],targY[i],targetdt[i]);
	}
	c1->cd(1);
	sprintf(name,"Target dE (Event: %d)",n); //JIIIIIIII
	target->SetTitle(name);
	target->GetZaxis()->SetRangeUser(0,2);
	target->Draw("colz");
	DrawTargetFrame();

	Kstp=new TMarker(u_x/3.1+11,u_y/3.1+10,24);
	Kstp->Draw("same");
	line=new TLine(u_x/3.1+11,u_y/3.1+10,u_x/3.1+11+px[0], u_y/3.1+10+py[0]);
	line->Draw("same");

	sprintf(name,"(%0.1f,%0.1f)",u_x,u_y);
	text=new TLatex(u_x/3.1+9,u_y/3.1+11,name);
	text->SetTextSize(0.03);
	text->Draw("same");

	c1->cd(2);
	sprintf(name,"Target dt (Event: %d)",n /* ,itrack*/); //JIIIIIIII
	time->SetTitle(name);
	//	time->GetZaxis()->SetRangeUser(0,0.1);
        //      time->GetZaxis()->SetLabelOffset(-0.02);
	//	time->GetZaxis()->SetRangeUser(0.,0.2);
	time->GetZaxis()->SetRangeUser(0.,40.);
	time->GetZaxis()->SetLabelOffset(-0.1);
	time->Draw("colz");
	DrawTargetFrame();

	Kstp=new TMarker(u_x/3.1+11,u_y/3.1+10,24);
	Kstp->Draw("same");
	line=new TLine(u_x/3.1+11,u_y/3.1+10,u_x/3.1+11+px[0], u_y/3.1+10+py[0]);
	line->Draw("same");

	sprintf(name,"(%0.1f,%0.1f)",u_x,u_y);
	text=new TLatex(u_x/3.1+9,u_y/3.1+11,name);
	text->SetTextSize(0.03);
	text->Draw("same");
	cout<<u_x<<" "<<u_y<<" "<<px[0]<<" "<<py[0]<<" "<<flush;
//return;

	c1->Update();
	c1->Print(output);
	target->Reset();
	time->Reset();

	if(n>1e5)break;
	}

	Kaon->Draw("targetE");
	c1->Print(output);

	sprintf(name,"%s]",output);
	c1->Print(name);

  return;
}



//
// Draw Box Line
//
void DrawTargetFrame(){
	for(int i=0;i<6;i++) DrawBox(7+i,18,8+i,19 , 1,0);
	for(int i=0;i<10;i++) DrawBox(5+i,17,6+i,18 , 1,0);
	for(int i=0;i<12;i++) DrawBox(4+i,16,5+i,17 , 1,0);
	for(int i=0;i<14;i++) DrawBox(3+i,15,4+i,16 , 1,0);
	for(int i=0;i<16;i++) DrawBox(2+i,14,3+i,15 , 1,0);
	for(int i=0;i<16;i++) DrawBox(2+i,13,3+i,14 , 1,0);
	for(int i=0;i<18;i++) DrawBox(1+i,12,2+i,13 , 1,0);
	for(int i=0;i<18;i++) DrawBox(1+i,11,2+i,12 , 1,0);
	for(int i=0;i<18;i++) DrawBox(1+i,10,2+i,11 , 1,0);
	for(int i=0;i<18;i++) DrawBox(1+i,9,2+i,10 , 1,0);
	for(int i=0;i<18;i++) DrawBox(1+i,8,2+i,9 , 1,0);
	for(int i=0;i<18;i++) DrawBox(1+i,7,2+i,8 , 1,0);
	for(int i=0;i<16;i++) DrawBox(2+i,6,3+i,7 , 1,0);
	for(int i=0;i<16;i++) DrawBox(2+i,5,3+i,6 , 1,0);
	for(int i=0;i<14;i++) DrawBox(3+i,4,4+i,5 , 1,0);
	for(int i=0;i<12;i++) DrawBox(4+i,3,5+i,4 , 1,0);
	for(int i=0;i<10;i++) DrawBox(5+i,2,6+i,3 , 1,0);
	for(int i=0;i<6;i++) DrawBox(7+i,1,8+i,2 , 1,0);
return;
}

void DrawBox(int x1, int y1, int x2, int y2, int color, int style){
    line=new TLine(x1,y1,x2,y1);
    line->SetLineColor(color);
    line->SetLineStyle(style);
    line->Draw("same");
    line=new TLine(x1,y2,x2,y2);
    line->SetLineColor(color);
    line->SetLineStyle(style);
    line->Draw("same");
    line=new TLine(x1,y1,x1,y2);
    line->SetLineColor(color);
    line->SetLineStyle(style);
    line->Draw("same");
    line=new TLine(x2,y1,x2,y2);
    line->SetLineColor(color);
    line->SetLineStyle(style);
    line->Draw("same");
}


