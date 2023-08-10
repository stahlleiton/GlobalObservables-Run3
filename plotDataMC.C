#include "TSystem.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "TClonesArray.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <TPaveText.h>
#include "RooCategory.h"
#include "RooNumIntConfig.h"
#include "RooPlotable.h"
#include <TUnfold.h>
#include <TLorentzVector.h>
#include <vector>
#include <TRandom.h>
#include <TF1.h>
#include <TObjArray.h>
#include <TEfficiency.h>
#include "TLegend.h"
#include <fstream>
#include "TLatex.h"

using namespace std;

int color(int i) {
  if (i==0) return kMagenta+1;
  else if (i==1) return kBlue+1;
  else if (i==2) return kRed+1;
  else if (i==3) return kCyan+1;
  else if (i==4) return kGreen+1;
  else if (i==5) return kOrange+1;
  else if (i==6) return 28;
  else if (i==7) return kBlack;
  else if (i==8) return kAzure+1;
  else if (i==9) return kYellow+1;
  else return kBlack;
}


void plotDataMC() {

  TFile * file = new TFile("CentralityTable_HFtowers200_DataPbPb_usingMC_d20230809_v1_TestRun_xSF0p86_ySFwithAllFilter_ThE4GeV_officiaMC2022__Threshold100__NOMINAL__Normalisation1000_4000.root");
  gStyle->SetOptStat(0);

  double nbin = 100;
  double axisMax = 6000;


  TCanvas* can = new TCanvas ("can","",900,1000);
  TH1F* axisHist = new TH1F("axisHist", "", nbin, 0, axisMax);
  TH1F* pullHist = new TH1F("pullHist","",nbin, 0, axisMax);

  //axisHist->GetYaxis()->SetRangeUser(0, 1.2);
  axisHist->GetYaxis()->SetTitle("Events");
  axisHist->GetXaxis()->SetTitle("");
  axisHist->GetXaxis()->SetNdivisions(505);
  axisHist->GetXaxis()->CenterTitle(true);
  axisHist->GetYaxis()->CenterTitle(true);
  axisHist->GetYaxis()->SetTitleSize(25);
  axisHist->GetYaxis()->SetTitleFont(43);
  axisHist->GetYaxis()->SetTitleOffset(1.4);
  axisHist->GetYaxis()->SetLabelFont(43);
  axisHist->GetYaxis()->SetLabelSize(20);

  //pullHist->GetYaxis()->SetRangeUser(0.68, 1.32);
  pullHist->GetYaxis()->SetTitle("data/MC");
  pullHist->GetYaxis()->SetNdivisions(505);
  pullHist->GetYaxis()->CenterTitle(true);
  pullHist->GetYaxis()->SetNdivisions(505);
  pullHist->GetYaxis()->SetTitleSize(25);
  pullHist->GetYaxis()->SetTitleFont(43);
  pullHist->GetYaxis()->SetTitleOffset(1.4);
  pullHist->GetYaxis()->SetLabelFont(43); 
  pullHist->GetYaxis()->SetLabelSize(20);

  pullHist->GetXaxis()->CenterTitle(true);
  pullHist->GetXaxis()->SetTitle("hiHF");
  pullHist->GetXaxis()->SetNdivisions(510);
  pullHist->GetXaxis()->SetTitleSize(25);
  pullHist->GetXaxis()->SetTitleFont(43);
  pullHist->GetXaxis()->SetTitleOffset(2.5);
  pullHist->GetXaxis()->SetLabelFont(43); 
  pullHist->GetXaxis()->SetLabelSize(25);
  pullHist->SetMarkerColor(kCyan+2);
  pullHist->SetMarkerStyle(kFullCircle);
  pullHist->SetMarkerSize(1);
  pullHist->SetLineColor(kCyan+2);
  pullHist->SetLineWidth(2);

  TLegend* leg = new TLegend(0.15,0.35,0.45,0.5);//(0.62,0.3,0.88,0.45);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  TLine * ly1 = new TLine(0, 1, axisMax, 1);
  ly1->SetLineColor(kRed);
  ly1->SetLineStyle(2);
  ly1->SetLineWidth(2);

  TLatex *  text = new TLatex(0.72 ,0.82,"CMS");
  text->SetNDC();
  text->SetTextFont(61);
  text->SetTextSize(0.075);
  text->SetLineWidth(8);
  
  TLatex *  text1 = new TLatex(0.7 ,0.77,"Internal");
  text1->SetNDC();
  text1->SetTextFont(52);
  text1->SetTextSize(0.055);
  text1->SetLineWidth(2);

  TLatex *  text2 = new TLatex(0.15 , 0.16, "#splitline{pprimaryVertexFilter &&}{#splitline{pclusterCompatibilityFilter &&}{HFTowerTh4 >=2}}");
  text2->SetNDC();
  text2->SetTextFont(42);
  text2->SetTextSize(0.04);
  text2->SetLineWidth(2);


  //cout<<"[INFO] looking for the file"<<endl;
  //TFile* f = TFile::Open(Form("%s",inputName.c_str()));
  //cout<<"[INFO] found the file"<<endl;

  TH1F* dataHist = (TH1F*) file->Get("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1205x02_offline/hfData");
  cout<<"[INFO] found the data hist"<<endl;
  TH1F* mcHist = (TH1F*) file->Get("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1205x02_offline/hfMc");
  cout<<"[INFO] found MC hist"<<endl;

  dataHist->GetXaxis()->SetRangeUser(0, axisMax);
  //dataHist->SetBins(nbin, 0, axisMax);
  dataHist->SetMarkerColor(kCyan+2);
  dataHist->SetMarkerStyle(kFullCircle);
  dataHist->SetMarkerSize(1);
  dataHist->SetLineColor(kCyan+2);
  dataHist->SetLineWidth(1);
  leg->AddEntry(dataHist, "Test Run Data 2022 ", "lp");

  mcHist->GetXaxis()->SetRangeUser(0, axisMax);
  //mcHist->SetBins(nbin, 0, axisMax);
  mcHist->SetMarkerColor(kMagenta+2);
  //mcHist->SetMarkerStyle(kFullCircle);
  //mcHist->SetMarkerSize(1);
  mcHist->SetLineColor(kMagenta+2);
  mcHist->SetLineWidth(2);
  leg->AddEntry(mcHist, "HYDJET Drum5F x 0.86", "lp");

  axisHist->GetYaxis()->SetRangeUser(0.009, 10* mcHist->GetMaximum());
  TLine * lx40 = new TLine(100, 0, 100, 10* mcHist->GetMaximum());
  lx40->SetLineColor(kRed);
  lx40->SetLineStyle(2);
  lx40->SetLineWidth(2);

  can->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.01); 
  pad1->Draw();
  pad1->cd();
  axisHist->SetStats(0);
  axisHist->Draw();
  mcHist->Draw("same hist");
  dataHist->Draw("same");
  leg->Draw("same");
  lx40->Draw("same");
  text->Draw("same");
  text1->Draw("same");
  text2->Draw("same");
  pad1->SetLogy();

  for (int i = 0; i<nbin; i++) {
    double hiHF = (axisMax/nbin)*i;
    double nom = dataHist->GetBinContent(dataHist->FindBin(hiHF));
    double den = mcHist->GetBinContent(mcHist->FindBin(hiHF));
    if (den!=0) {
      pullHist->SetBinContent(pullHist->FindBin(hiHF),nom*1.0/den);
      pullHist->SetBinError(pullHist->FindBin(hiHF),0.00001);
    }
  }

  can->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetTopMargin(0.01);
  //pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();
  pullHist->Draw();
  ly1->Draw("same");

  can->SaveAs("dataMcComparaison.pdf");
}
