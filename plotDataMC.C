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

  TFile * file = new TFile("CentralityTable_HFtowers200_DataPbPb_usingMC_v1_TestRun_xSF0p86_ySFwithAllFilter_ThE5GeV_officiaMC2022__Threshold100__NOMINAL__Normalisation1000_4000__GT_Aug10_NEW.root");
  gStyle->SetOptStat(0);

  double nbin = 25;
  double axisMax = 500;


  TCanvas* can = new TCanvas ("can","",900,1000);
  TH1F* axisHist = new TH1F("axisHist", "", nbin, 0, axisMax);
  TH1F* pullHist1 = new TH1F("pullHist1","",nbin, 0, axisMax);
  TH1F* pullHist2 = new TH1F("pullHist2","",nbin, 0, axisMax);

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

  for (auto& pullHist : {pullHist1, pullHist2}) {
  pullHist->GetYaxis()->SetRangeUser(0.7, 1.3);
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
  pullHist->SetMarkerColor(kRed+2);
  pullHist->SetMarkerStyle(kFullCircle);
  pullHist->SetMarkerSize(1);
  pullHist->SetLineColor(kRed+2);
  pullHist->SetLineWidth(2);
  }
  pullHist2->SetMarkerColor(kGreen+2);
  pullHist2->SetLineColor(kGreen+2);

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

  TLatex *  text2 = new TLatex(0.15 , 0.16, "#splitline{HLT_HIMinimumBias &&}{#splitline{pprimaryVertexFilter &&}{#splitline{pclusterCompatibilityFilter &&}{HFTowerTh4 >=2}}}");
  text2->SetNDC();
  text2->SetTextFont(42);
  text2->SetTextSize(0.04);
  text2->SetLineWidth(2);


  //cout<<"[INFO] looking for the file"<<endl;
  //TFile* f = TFile::Open(Form("%s",inputName.c_str()));
  //cout<<"[INFO] found the file"<<endl;

  TH1F* dataHist1 = (TH1F*) file->Get("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1205x02_offline/hfData1");
  if (!dataHist1) return;
  TH1F* dataHist2 = (TH1F*) file->Get("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1205x02_offline/hfData2");
  if (!dataHist2) return;
  cout<<"[INFO] found the data hist"<<endl;
  TH1F* mcHist1 = (TH1F*) file->Get("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1205x02_offline/hfMc1");
  TH1F* mcHist2 = (TH1F*) file->Get("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1205x02_offline/hfMc1");
  cout<<"[INFO] found MC hist"<<endl;
  TH1F* combHist = (TH1F*) file->Get("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1205x02_offline/hfCombined");
  cout<<"[INFO] found combined hist"<<endl;

  const auto nEntries1 = dataHist1->Integral(dataHist1->FindBin(1000), dataHist1->FindBin(4000));
  const auto nEntries2 = dataHist2->Integral(dataHist2->FindBin(1000), dataHist2->FindBin(4000));
  const auto nEntries3 = mcHist2->Integral(mcHist2->FindBin(1000), mcHist2->FindBin(4000));
  dataHist2->Scale(nEntries1/nEntries2);
  mcHist2->Scale(nEntries1/nEntries3);

  for (auto& dataHist : {dataHist1, dataHist2}) {
  dataHist->GetXaxis()->SetRangeUser(0, axisMax);
  //dataHist->SetBins(nbin, 0, axisMax);
  dataHist->SetMarkerColor(kRed+2);
  dataHist->SetMarkerStyle(kFullCircle);
  dataHist->SetMarkerSize(1);
  dataHist->SetLineColor(kRed+2);
  dataHist->SetLineWidth(1);
  }
  dataHist2->SetMarkerColor(kGreen+2);
  dataHist2->SetLineColor(kGreen+2);
  leg->AddEntry(dataHist1, "Test Run Data 2022 run <= 363320", "lp");
  leg->AddEntry(dataHist2, "Test Run Data 2022 run > 363320", "lp");

  for (auto& mcHist : {mcHist1, mcHist2}) {
  mcHist->GetXaxis()->SetRangeUser(0, axisMax);
  //mcHist->SetBins(nbin, 0, axisMax);
  mcHist->SetMarkerColor(kMagenta+2);
  //mcHist->SetMarkerStyle(kFullCircle);
  //mcHist->SetMarkerSize(1);
  mcHist->SetLineColor(kMagenta+2);
  mcHist->SetLineWidth(2);
  }
  leg->AddEntry(mcHist1, "HYDJET Drum5F x 0.86", "lp");

  combHist->GetXaxis()->SetRangeUser(0, axisMax);
  //combHist->SetBins(nbin, 0, axisMax);
  combHist->SetMarkerColor(kBlue+2);
  //combHist->SetMarkerStyle(kFullCircle);
  //combHist->SetMarkerSize(1);
  combHist->SetLineColor(kBlue+2);
  combHist->SetLineWidth(2);
  leg->AddEntry(combHist, "Combined", "lp");

  axisHist->GetYaxis()->SetRangeUser(0.009, 10* mcHist1->GetMaximum());
  TLine * lx40 = new TLine(100, 0, 100, 10* mcHist1->GetMaximum());
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
  mcHist1->Draw("same hist");
  dataHist1->Draw("same");
  dataHist2->Draw("same");
  combHist->Draw("same");
  leg->Draw("same");
  lx40->Draw("same");
  text->Draw("same");
  text1->Draw("same");
  text2->Draw("same");
  pad1->SetLogy();

  for (int i = 0; i<nbin; i++) {
    double hiHF = (axisMax/nbin)*i;
    double nom1 = dataHist1->GetBinContent(dataHist1->FindBin(hiHF));
    double nom2 = dataHist2->GetBinContent(dataHist2->FindBin(hiHF));
    double den = mcHist2->GetBinContent(mcHist2->FindBin(hiHF));
    if (den!=0) {
      pullHist1->SetBinContent(pullHist1->FindBin(hiHF),nom1*1.0/den);
      pullHist1->SetBinError(pullHist1->FindBin(hiHF),0.00001);
      pullHist2->SetBinContent(pullHist2->FindBin(hiHF),nom2*1.0/den);
      pullHist2->SetBinError(pullHist2->FindBin(hiHF),0.00001);
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
  pullHist2->Draw();
  pullHist1->Draw("same");
  ly1->Draw("same");

  can->SaveAs("dataMcComparaison.pdf");
}
