//########################################################################################################
// 
// sbn_sens.C -- ROOT macro for generating the SBN Program sensitivity 
//               plots for nue appearance and numu disappearance 
//
//########################################################################################################

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TCanvas.h"
#include "TGraph.h"


void loadStyle(){

    
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.15);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0000);

  gStyle->SetOptTitle(0);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTitleFont(42,"x");
  gStyle->SetTitleFont(42,"y");
  gStyle->SetTitleFont(42,"z");
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetTitleSize(0.06,"z");
  gStyle->SetTitleOffset(1.15,"x");
  gStyle->SetTitleOffset(1.2,"y");
  gStyle->SetLabelSize(0.045,"x");
  gStyle->SetLabelSize(0.045,"y");
  gStyle->SetHistLineStyle(1);

  gStyle->SetFrameLineWidth(2);
  gStyle->SetLineWidth(1);
  gStyle->SetHistLineWidth(2);

}

//========================================================================================================
//
// Simple helper funtion to add text to a Canvas
//
//========================================================================================================
void add_plot_label( char* label, double x, double y, 
		     double size = 0.035, int color = 1, int font = 62, int align = 12 ){

  TLatex *latex = new TLatex( x, y, label );
  latex->SetNDC();
  latex->SetTextSize(size);
  latex->SetTextColor(color);
  latex->SetTextFont(font);
  latex->SetTextAlign(align);
  latex->Draw();

}

//========================================================================================================
//
// Read in LSND data files and add allowed regions to the provided Canvas and Legend 
//
//========================================================================================================
void lsnd_plot( TCanvas* c, TLegend *leg, bool draw_bf = false ){

  c->cd();
  const char* data_dir = "lsnd_data/";
  Double_t  dm2BF[] = {1.2};
  Double_t sin22thBF[] = {0.003};

  const Int_t NDATAFILES = 11;
  const char * file_list[NDATAFILES] = {"llreg_608_1.vec",
					"llreg_608_2.vec",
					"llreg_608_3.vec",
					"llreg_607_1.vec",
					"llreg_607_2.vec",
					"llreg_607_3.vec",
					"llreg_607_4.vec",
					"llreg_607_5.vec",
					"llreg_607_6.vec",
					"llreg_607_7.vec",
					"llreg_607_8.vec"};
  
  Int_t graph_color[NDATAFILES] = {29, 29, 29, 38, 38, 38, 38, 38, 38, 38, 38};

  const Int_t NDATAPOINTS = 500;
  Int_t    nlines;
  Double_t x[NDATAPOINTS], y[NDATAPOINTS];
  Double_t dummy, dummy_old;
  TGraph* gr[NDATAFILES];
  char  filename[100];
  ifstream datafile;
  
  for (Int_t ifile = 0; ifile<NDATAFILES; ifile++) {
    
    nlines = 0;
    for (Int_t i=0;i<NDATAPOINTS;i++){x[i]=0.0;y[i]=0.0;}
    
    strcpy(filename, data_dir);
    strcat(filename, file_list[ifile]);
    datafile.open(filename, ios_base::in);

    //check if the file is open: 
    if (!datafile.is_open() ) {std::cerr << "lsnd_plot.C: file not opened" <<std::endl; return;}
    else {std::cout << "Successfully opened " << filename << std::endl;}
    
    while (!datafile.eof()) {
      datafile >> dummy; 
      datafile >> dummy; 
      datafile >> x[nlines]; 
      datafile >> y[nlines];
      nlines++;
      if (dummy == dummy_old) nlines--; //if last row was empty
      dummy_old = dummy;
    }

    gr[ifile] = new TGraph(nlines,x,y);
    datafile.close();
  }
  std::cout << "Finished reading data files" << std::endl;

  // Draw contours
  for (Int_t ifile = 0; ifile<NDATAFILES; ifile++) {
    gr[ifile]->SetFillColor(graph_color[ifile]);
    gr[ifile]->Draw("LF");
  }
 
  // Add the best fit point;
  TGraph * bfPoint = new TGraph(1, sin22thBF, dm2BF);
  bfPoint->SetLineColor(2);
  bfPoint->SetMarkerStyle(5);
  bfPoint->SetMarkerSize(0.6);
  bfPoint->SetMarkerColor(1);
  if( draw_bf) bfPoint->Draw("P");

  // Add legend items
  leg->AddEntry(gr[NDATAFILES-1],"LSND 90%","f");
  leg->AddEntry(gr[0],"LSND 99%","f");
  if( draw_bf ) leg->AddEntry(bfPoint,"LSND best fit","p");

  return;
}

//======================================================================================================
// 
// Add C. Giunti et. al. global fit contours
//    channel = nue appearance
//
//====================================================================================================== 
void giunti_global_nue(TCanvas* c, TLegend *leg){

  c->cd();
  const char* data_dir = "giunti_global_fits/global/";
  Double_t  dm2BF[] = {1.7};
  Double_t sin22thBF[] = {0.0011};

  const Int_t NDATAFILES = 6;
  const char * file_list[NDATAFILES] = {"nue-app-3sigma-1.dat",
					"nue-app-3sigma-2.dat",
					"nue-app-3sigma-3.dat",
					"nue-app-2sigma-1.dat",
					"nue-app-2sigma-2.dat",
					"nue-app-1sigma.dat"};
  
  Int_t graph_color[NDATAFILES] = {kGreen+2, kGreen+2, kGreen+2, kGreen, kGreen, kYellow};

  const Int_t NDATAPOINTS = 500;
  Int_t    nlines;
  Double_t x[NDATAPOINTS], y[NDATAPOINTS], test;
  TGraph* gr[NDATAFILES];
  char  filename[100];
  ifstream datafile;

  for( Int_t ifile = 0; ifile < NDATAFILES; ifile++ ){

    nlines = 0;
    for (Int_t i=0;i<NDATAPOINTS;i++){x[i]=0.0;y[i]=0.0;} 

    strcpy(filename, data_dir);
    strcat(filename, file_list[ifile]); 
    datafile.open(filename, ios_base::in);
 
    //confirm the file is open: 
    if (!datafile.is_open() ) {std::cerr << filename << " not opened" << std::endl; return;}
    else {std::cout << "Successfully opened " << filename << std::endl;}
    
    while (!datafile.eof()) {
      datafile >> x[nlines]; 
      datafile >> y[nlines];
      test = x[nlines];
      if( test != 0 )
	nlines++;
    }

    gr[ifile] = new TGraph(nlines,x,y);
    datafile.close();
  }

  // Set graph colors and draw
  for( Int_t ifile = 0; ifile < NDATAFILES; ifile++){
    gr[ifile]->SetFillColor(graph_color[ifile]);
    gr[ifile]->Draw("LF");
  }

  // Add the best fit point
  TGraph * bfPoint = new TGraph(1, sin22thBF, dm2BF);
  bfPoint -> SetMarkerStyle(34);
  bfPoint -> SetMarkerColor(1);
  bfPoint -> Draw("P");

  // Add to legend
  leg->AddEntry(gr[NDATAFILES-1],"Global 2017 1#sigma","f");
  leg->AddEntry(gr[NDATAFILES-2],"Global 2017 2#sigma","f");
  leg->AddEntry(gr[0],"Global 2017 3#sigma","f");
  leg->AddEntry(bfPoint,"Global 2017 best fit","p");
  
  return;
}

//======================================================================================================
// 
// Add C. Giunti et. al. global fit contours
//    channel = nue appearance
//
//====================================================================================================== 
void giunti_global_numu(TCanvas* c, TLegend *leg){

  c->cd();
  const char* data_dir = "giunti_global_fits/global/";
  Double_t  dm2BF[] = {1.7};
  Double_t sin22thBF[] = {0.059};

  const Int_t NDATAFILES = 6;
  const char * file_list[NDATAFILES] = {"numu-dis-3sigma-1.dat",
					"numu-dis-3sigma-2.dat",
					"numu-dis-3sigma-3.dat",
					"numu-dis-2sigma-1.dat",
					"numu-dis-2sigma-2.dat",
					"numu-dis-1sigma.dat"};
  
  Int_t graph_color[NDATAFILES] = {kGreen+2, kGreen+2, kGreen+2, kGreen, kGreen, kYellow};

  const Int_t NDATAPOINTS = 500;
  Int_t    nlines;
  Double_t x[NDATAPOINTS], y[NDATAPOINTS], test;
  TGraph* gr[NDATAFILES];
  char  filename[100];
  ifstream datafile;

  for( Int_t ifile = 0; ifile < NDATAFILES; ifile++ ){

    nlines = 0;
    for (Int_t i=0;i<NDATAPOINTS;i++){x[i]=0.0;y[i]=0.0;} 

    strcpy(filename, data_dir);
    strcat(filename, file_list[ifile]); 
    datafile.open(filename, ios_base::in);
 
    //confirm the file is open: 
    if (!datafile.is_open() ) {std::cerr << filename << " not opened" << std::endl; return;}
    else {std::cout << "Successfully opened " << filename << std::endl;}
    
    while (!datafile.eof()) {
      datafile >> x[nlines]; 
      datafile >> y[nlines];
      test = x[nlines];
      if( test != 0 )
	nlines++;
    }

    gr[ifile] = new TGraph(nlines,x,y);
    datafile.close();
  }

  // Set graph colors and draw
  for( Int_t ifile = 0; ifile < NDATAFILES; ifile++){
      gr[ifile]->SetFillColor(graph_color[ifile]);
      gr[ifile]->Draw("LF");
  }

  // Add the best fit point
  TGraph * bfPoint = new TGraph(1, sin22thBF, dm2BF);
  bfPoint -> SetMarkerStyle(34);
  bfPoint -> SetMarkerColor(1);
  bfPoint -> Draw("P");

  // Add to legend
  leg->AddEntry(gr[NDATAFILES-1],"Global 2017 1#sigma","f");
  leg->AddEntry(gr[NDATAFILES-2],"Global 2017 2#sigma","f");
  leg->AddEntry(gr[0],"Global 2017 3#sigma","f");
  leg->AddEntry(bfPoint,"Global 2017 best fit","p");

  
  return;
}



//======================================================================================================
//
// Add MINOS/MINOS+ 90%, 3 sigma, 5 sigma contours from May 2018
//    using published data release & assuming a 2 sided, 2 dof chi2 cut
//    channel = numu appearance
//
//======================================================================================================

void minos_sens( TCanvas* c, TLegend *leg, bool minos90, bool minos3s, bool minos5s ) {
    TFile* minos_data = new TFile("dataRelease.root", "READ");
    std::cout << "Successfully opened dataRelease.root" << std::endl;
    
    TH2D* dm241vsth24 = (TH2D*)minos_data->Get("dm241vsth24");
    
    // SET CONTOURS //
    double x2min = round(dm241vsth24->GetMinimum()*100.0)/100.0;
    double contours[3];
    contours[0] = x2min + 4.61; //90%
    contours[1] = x2min + 9.00; //3sigma
    contours[2] = x2min + 24.37; //5sigma
    dm241vsth24->SetContour(3, contours);
    
    TCanvas* c1_minos = new TCanvas();
    dm241vsth24->Draw("CONT Z LIST");
    c1_minos->Update();
    
    
    // X-AXIS UNIT CONVERSTION //
    TObjArray* conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    TList* contLevel = NULL;
    TGraph* curv = NULL;
    TGraph* gc = NULL;
    double xval, yval, zval;
    
    c1_minos->Close();
    
    //vector of vectors: top level vector is the contour level
    //each vector inside is for each graph within that contour level
    
    vector<vector<double>> x90;
    vector<vector<double>> y90;
    
    vector<vector<double>> x3sig;
    vector<vector<double>> y3sig;
    
    vector<vector<double>> x5sig;
    vector<vector<double>> y5sig;
    
    if (conts == NULL) {
        std::cout << "No contours extracted!" << std::endl;
    }
    
    
    
    else {
        for (int i=0; i<conts->GetSize(); i++) {
            contLevel = (TList*)conts->At(i);
            zval = contours[i];
            
            //get first graph in contour level
            curv = (TGraph*)contLevel->First();
            for (int j=0; j<contLevel->GetSize(); j++) { //loop over graph j in contour level
                vector<double> xpts;
                vector<double> ypts;
                //loop over each point in the graph
                for (int k=0; k<curv->GetN(); k++) {
                    curv->GetPoint(k, xval, yval);
                    xpts.push_back(TMath::Sin(2*xval)*TMath::Sin(2*xval));
                    ypts.push_back(yval);

                }
                
                //add to appropriate top level vector
                if (i==0) {
                    x90.push_back(xpts);
                    y90.push_back(ypts);
                }
                if (i==1) {
                    x3sig.push_back(xpts);
                    y3sig.push_back(ypts);
                }
                else if (i==2) {
                    x5sig.push_back(xpts);
                    y5sig.push_back(ypts);
                }
                
                curv=(TGraph*)contLevel->After(curv); //get the next graph in contour level
                
            }
        }
        
        // PLOTS //
        
        TGraph* cl90_0 = new TGraph(x90[0].size(), &x90[0][0], &y90[0][0]);
        TGraph* cl90_1 = new TGraph(x90[1].size(), &x90[1][0], &y90[1][0]);
        TGraph* cl90_2 = new TGraph(x90[2].size(), &x90[2][0], &y90[2][0]);
        
        TGraph* cl3sig_0 = new TGraph(x3sig[0].size(), &x3sig[0][0], &y3sig[0][0]);
        TGraph* cl3sig_1 = new TGraph(x3sig[1].size(), &x3sig[1][0], &y3sig[1][0]);
        
        TGraph* cl5sig_0 = new TGraph(x5sig[0].size(), &x5sig[0][0], &y5sig[0][0]);
        
        c->cd();
        
        if (minos90 == true) {
            cl90_0->SetLineColor(kBlack);
            cl90_0->SetLineStyle(2);
            cl90_0->SetLineWidth(1);
            cl90_0->Draw("C");
            leg->AddEntry(cl90_0, "MINOS/MINOS+ 90%");
        }
        
        if (minos3s == true) {
            cl3sig_0->SetLineColor(kBlue-10);
            cl3sig_0->SetLineStyle(1);
            cl3sig_0->SetLineWidth(2);
            cl3sig_0->Draw("C");
            leg->AddEntry(cl3sig_0, "MINOS/MINOS+ 3#sigma");
        }
        
        if (minos5s == true) {
            cl5sig_0->SetLineColor(kBlue-10);
            cl5sig_0->SetLineStyle(2);
            cl5sig_0->SetLineWidth(2);
            cl5sig_0->Draw("C");
            leg->AddEntry(cl5sig_0, "MINOS/MINOS+ 5#sigma");
        }
        
        
        
    }
    
    
    return;


}


//======================================================================================================
//
// Generate the SBN nue appearance sensitivity plot
//
//======================================================================================================
void sbn_nue_plot( bool w90 = false, bool w3s = true, bool w5s = true ){

  loadStyle();

  TCanvas * d = new TCanvas("SBN Nue Sensitivity", "SBN Nue Sensitivity", 600, 470);
  d->SetLogx();
  d->SetLogy();
  
  TH2D* hr1 = new TH2D("hr1","hr1",500,0.0001,1,500,0.02,50);
  hr1->Reset();
  hr1->SetFillColor(0);
  hr1->SetTitle(";sin^{2}2#theta_{#mue};#Deltam^{2} (eV^{2})");
  hr1->GetXaxis()->SetTitleOffset(1.5);
  hr1->GetYaxis()->SetTitleOffset(1.0);
  hr1->GetXaxis()->SetTitleSize(0.05);
  hr1->GetYaxis()->SetTitleSize(0.05);
  hr1->GetXaxis()->CenterTitle();
  hr1->GetYaxis()->CenterTitle();
  hr1->SetStats(kFALSE);
  hr1->Draw();
  
  TLegend *legend = new TLegend(0.59,0.52,0.82,0.845);
  legend->SetTextSize(0.022);


  lsnd_plot(d, legend, false);
  giunti_global_nue(d, legend);


  Int_t    nlines = 0;
  Double_t x90[2100],y90[2100];
  Double_t x3[2000],y3[2000];
  Double_t x5[500],y5[500];
  ifstream datafile;
  char firstline[500];
  int c;

  //-----------------------
  //Get the 90% line
  //-----------------------
  datafile.open("sbn_contours/nue90.txt", ios_base::in);

  //check if the file is open: 
  if (!datafile.is_open() ) {std::cerr << "file not opened" <<std::endl; return;}
  else {std::cout << "Successfully opened sbn_contours/nue3s.txt" << std::endl;}
  
  datafile.getline( firstline, 500 );
  //std::cout << firstline << std::endl;
  
  nlines = 0;
  while (!datafile.eof()) {
    datafile >> x90[nlines];
    c = datafile.get();
    datafile >> y90[nlines];
    nlines++;
  }
  datafile.close();
  nlines--; //throw away the last, and likely corrupt entry

  std::cout << "Number entries 90%: " << nlines << std::endl;

  TGraph *nue90_curve = new TGraph(nlines,x90,y90);

  //-----------------------
  //Get the 3 sigma line
  //-----------------------
  datafile.open("sbn_contours/nue3s.txt", ios_base::in);

  //check if the file is open: 
  if (!datafile.is_open() ) {std::cerr << "file not opened" <<std::endl; return;}
  else {std::cout << "Successfully opened sbn_contours/nue3s.txt" << std::endl;}
  
  datafile.getline( firstline, 500 );
  //std::cout << firstline << std::endl;

  nlines = 0;  
  while (!datafile.eof()) {
    datafile >> x3[nlines];
    c = datafile.get();
    datafile >> y3[nlines];
    nlines++;
  }
  datafile.close();
  nlines--; //throw away the last, and likely corrupt entry

  std::cout << "Number entries 3 sigma: " << nlines << std::endl;

  TGraph *nue3s_curve = new TGraph(nlines,x3,y3);
  
  //-----------------------
  //Get the 5 sigma line
  //-----------------------
  datafile.open("sbn_contours/nue5s.txt", ios_base::in);

  //check if the file is open: 
  if (!datafile.is_open() ) {std::cerr << "file not opened" <<std::endl; return;}
  else {std::cout << "Successfully opened sbn_contours/nue5s.txt" << std::endl;}
  
  datafile.getline( firstline, 500 );
  //std::cout << firstline << std::endl;
  
  nlines = 0;
  while (!datafile.eof()) {
    datafile >> x5[nlines];
    c = datafile.get();
    datafile >> y5[nlines];
    nlines++;
  }
  datafile.close();
  nlines--; //throw away the last, and likely corrupt entry

  std::cout << "Number entries 5 sigma: " << nlines << std::endl;

  TGraph *nue5s_curve = new TGraph(nlines,x5,y5);

  //-------------------------------
  // Draw the sensitivity contours
  //-------------------------------
  nue90_curve->SetLineColor(kBlack);
  nue90_curve->SetLineStyle(1);
  nue90_curve->SetLineWidth(1);
  if( w90 ) nue90_curve->Draw("C");

  nue3s_curve->SetLineColor(kRed);
  nue3s_curve->SetLineStyle(1);
  nue3s_curve->SetLineWidth(2);
  if( w3s) nue3s_curve->Draw("C");

  nue5s_curve->SetLineColor(kRed);
  nue5s_curve->SetLineStyle(2);
  nue5s_curve->SetLineWidth(3);
  if( w5s ) nue5s_curve->Draw("C");

  if( w90 ) legend->AddEntry(nue90_curve,"SBN 90%","l");
  if( w3s ) legend->AddEntry(nue3s_curve,"SBN 3#sigma","l");
  if( w5s ) legend->AddEntry(nue5s_curve,"SBN 5#sigma","l");
  
  legend->Draw();

  char label[200];

  if( w90 || w3s || w5s ){ 
    sprintf( label, "SBN sensitivities assume exposures of:");
    add_plot_label( label, 0.185, 0.31, 0.024, 1, 42, 12 );
    sprintf( label, "6.60#times10^{20} protons on target in ICARUS and SBND");
    add_plot_label( label, 0.195, 0.28, 0.024, 1, 42, 12 );
    sprintf( label, "13.2#times10^{20} protons on target in MicroBooNE");
    add_plot_label( label, 0.195, 0.25, 0.024, 1, 42, 12 );
    sprintf( label, "Global 2017: S. Gariazzo et al., arXiv:1703.00860 [hep-ph]");
    add_plot_label( label, 0.185, 0.20, 0.024, 1, 42, 12 );
  }
  else{
    sprintf( label, "Global 2017:");
    add_plot_label( label, 0.185, 0.21, 0.027, 1, 62, 11 );
    sprintf( label, "S. Gariazzo et al., arXiv:1703.00860 [hep-ph]");
    add_plot_label( label, 0.3, 0.21, 0.027, 1, 42, 11 );
  }
  
  sprintf( label, "#nu_{#mu} #rightarrow #nu_{e}  appearance");
  add_plot_label( label, 0.7, 0.87, 0.03, 1, 62, 22 );

  d->RedrawAxis();
  TLine l;
  l.SetLineWidth(2);
  l.SetLineColor(kBlack);
  l.DrawLine(0,hr1->GetYaxis()->GetBinUpEdge(hr1->GetNbinsY()),hr1->GetXaxis()->GetBinUpEdge(hr1->GetNbinsX()),hr1->GetYaxis()->GetBinUpEdge(hr1->GetNbinsY()));
  
  char filename[100];
  strcpy(filename, "SBN_nue");
  if( w90 || w3s || w5s ) strcat(filename, "_sensitivity");
  if( w90 ) strcat(filename, "_90");
  if( w3s ) strcat(filename, "_3s");
  if( w5s ) strcat(filename, "_5s");
  strcat(filename, ".pdf");
  d->Print(filename);
  
  return;
}

//====================================================================================================
//
// Generate the SBN numu appearance sensitivity plot
//
//====================================================================================================
void sbn_numu_plot( bool w90 = false, bool w3s = true, bool w5s = false, bool minos = true ){

  loadStyle();

  TCanvas * d = new TCanvas("SBN Numu Sensitivity", "SBN Numu Sensitivity", 600, 470);
  d->SetLogx();
  d->SetLogy();
  
  TH2D* hr1 = new TH2D("hr1","hr1",500,0.01,1,500,0.02,50);
  hr1->Reset();
  hr1->SetFillColor(0);
  hr1->SetTitle(";sin^{2}2#theta_{#mu#mu};#Deltam^{2} (eV^{2})");
  hr1->GetXaxis()->SetTitleOffset(1.5);
  hr1->GetYaxis()->SetTitleOffset(1.0);
  hr1->GetXaxis()->SetTitleSize(0.05);
  hr1->GetYaxis()->SetTitleSize(0.05);
  hr1->GetXaxis()->CenterTitle();
  hr1->GetYaxis()->CenterTitle();
  hr1->SetStats(kFALSE);
  hr1->Draw();

  TLegend *legend = new TLegend(0.59,0.59,0.82,0.83);
  legend->SetTextSize(0.022);


  giunti_global_numu(d, legend);


  Int_t nlines = 0;
  Double_t x90[3000],y90[3000];
  Double_t x3[2000],y3[2000];
  Double_t x5[500],y5[500];
  ifstream datafile;
  char firstline[500];
  int c;

  //-------------------------------
  //Get the 90% line
  //-------------------------------
  datafile.open("sbn_contours/numu90.txt", ios_base::in);

  //check if the file is open: 
  if (!datafile.is_open() ) {std::cerr << "file not opened" <<std::endl; return;}
  else {std::cout << "Successfully opened sbn_contours/numu90.txt" << std::endl;}
  
  datafile.getline( firstline, 500 );
 
  nlines = 0;
  while (!datafile.eof()) {
    datafile >> x90[nlines];
    c = datafile.get();
    datafile >> y90[nlines];
    nlines++;
  }
  datafile.close();
  nlines--; //throw away the last, and likely corrupt entry

  std::cout << "Number entries 90%: " << nlines << std::endl;

  TGraph *numu90_curve = new TGraph(nlines,x90,y90);

  //-------------------------------
  //Get the 3 sigma line
  //-------------------------------
  datafile.open("sbn_contours/numu3s.txt", ios_base::in);

  //check if the file is open: 
  if (!datafile.is_open() ) {std::cerr << "file not opened" <<std::endl; return;}
  else {std::cout << "Successfully opened sbn_contours/numu3s.txt" << std::endl;}
  
  datafile.getline( firstline, 500 );
 
  nlines = 0;
  while (!datafile.eof()) {
    datafile >> x3[nlines];
    c = datafile.get();
    datafile >> y3[nlines];
    nlines++;
  }
  datafile.close();
  nlines--; //throw away the last, and likely corrupt entry

  std::cout << "Number entries 3 sigma: " << nlines << std::endl;

  TGraph *numu3s_curve = new TGraph(nlines,x3,y3);
  
  //-----------------------------
  //Get the 5 sigma line
  //-----------------------------
  datafile.open("sbn_contours/numu5s.txt", ios_base::in);

  //check if the file is open: 
  if (!datafile.is_open() ) {std::cerr << "file not opened" <<std::endl; return;}
  else {std::cout << "Successfully opened sbn_contours/numu5s.txt" << std::endl;}
  
  datafile.getline( firstline, 500 );
  
  nlines = 0;
  while (!datafile.eof()) {
    datafile >> x5[nlines];
    c = datafile.get();
    datafile >> y5[nlines];
    nlines++;
  }
  datafile.close();
  nlines--; //throw away the last, and likely corrupt entry

  std::cout << "Number entries 5 sigma: " << nlines << std::endl;

  TGraph *numu5s_curve = new TGraph(nlines,x5,y5);
    

  //-------------------------------
  // Draw the sensitivity contours
  //-------------------------------
  numu90_curve->SetLineColor(kBlack);
  numu90_curve->SetLineStyle(1);
  numu90_curve->SetLineWidth(1);
  if( w90 ) numu90_curve->Draw("C");

  numu3s_curve->SetLineColor(kRed);
  numu3s_curve->SetLineStyle(1);
  numu3s_curve->SetLineWidth(2);
  if( w3s ) numu3s_curve->Draw("C");

  numu5s_curve->SetLineColor(kRed);
  numu5s_curve->SetLineStyle(2);
  numu5s_curve->SetLineWidth(3);
  if( w5s ) numu5s_curve->Draw("C");
    
  if( w90 ) legend->AddEntry(numu90_curve,"SBN 90%","l");
  if( w3s ) legend->AddEntry(numu3s_curve,"SBN 3#sigma","l");
  if( w5s ) legend->AddEntry(numu5s_curve,"SBN 5#sigma","l");
    
    
    
  // MINOS/MINOS+ SENSITIVITY //
  if ( minos ) {
      minos_sens(d, legend, w90, w3s, w5s);
  }
    
    
    
  legend->Draw();

  char label[200];

  if( (w90 || w3s || w5s) && (minos == false) ){
    sprintf( label, "SBN sensitivities assume exposures of:");
    add_plot_label( label, 0.185, 0.31, 0.024, 1, 42, 12 );
    sprintf( label, "6.60#times10^{20} protons on target in ICARUS and SBND");
    add_plot_label( label, 0.195, 0.28, 0.024, 1, 42, 12 );
    sprintf( label, "13.2#times10^{20} protons on target in MicroBooNE");
    add_plot_label( label, 0.195, 0.25, 0.024, 1, 42, 12 );
    sprintf( label, "Global 2017: S. Gariazzo et al., arXiv:1703.00860 [hep-ph]");
    add_plot_label( label, 0.185, 0.20, 0.024, 1, 42, 12 );
  }
  else{
    sprintf( label, "Global 2017:");
    add_plot_label( label, 0.185, 0.21, 0.027, 1, 62, 11 );
    sprintf( label, "S. Gariazzo et al., arXiv:1703.00860 [hep-ph]");
    add_plot_label( label, 0.3, 0.21, 0.027, 1, 42, 11 );
  }
  
  sprintf( label, "#nu_{#mu} disappearance");
  add_plot_label( label, 0.7, 0.85, 0.03, 1, 62, 22 );

  gPad->RedrawAxis();
  TLine l;
  l.SetLineWidth(2);
  l.SetLineColor(kBlack);
  l.DrawLine(0,hr1->GetYaxis()->GetBinUpEdge(hr1->GetNbinsY()),hr1->GetXaxis()->GetBinUpEdge(hr1->GetNbinsX()),hr1->GetYaxis()->GetBinUpEdge(hr1->GetNbinsY()));;

  char filename[100];
  strcpy(filename, "SBN_numu");
  if( w90 || w3s || w5s ) strcat(filename, "_sensitivity");
  if( w90 ) strcat(filename, "_90");
  if( w3s ) strcat(filename, "_3s");
  if( w5s ) strcat(filename, "_5s");
  strcat(filename, ".pdf");
  d->Print(filename);
  
  return;
}

//#################################################################################################
//
// END sbn_sens.C
//
//#################################################################################################
