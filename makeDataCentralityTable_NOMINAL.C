#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TParameter.h"
#include "TMath.h"
#include "TString.h"
#include "TChain.h"
#include "TEfficiency.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <numeric>


void makeDataCentralityTable_NOMINAL(const int hfCoinN = 2,
                                     const int hfCoinThr = 5,
                                     const double threshold = 100.0,
                                     const size_t nbins = 200,
                                     const std::string tag = "CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1205x02_offline")
{
  // Constant parameters
  const auto mcXscale = 0.86;
  const auto threshold_norm = 1000.0;
  const auto thresholdMax = 4000.0;

  // Process data
  const std::string inFileName = "/eos/cms/store/group/phys_heavyions/anstahll/GO2023/HIForest_HIMinimumBias_HIRun2022A_20230811_SKIM.root";
  TFile inFile(inFileName.c_str(), "READ");
  if (!inFile.IsOpen()) throw std::logic_error("Data file was not found!"); 
  const auto& t = inFile.Get<TTree>("hiEvtAnalyzer/HiTree");
  const auto& tskimanalysis = inFile.Get<TTree>("skimanalysis/HltTree");
  const auto& thltanalysis = inFile.Get<TTree>("hltanalysis/HltTree");
  t->AddFriend(tskimanalysis);
  t->AddFriend(thltanalysis);

  UInt_t run;
  t->SetBranchAddress("run", &run);
  std::map<std::string, int> varI, mcVarI;
  const char* numMinHFTowerLbl = Form("numMinHFTower%d", hfCoinThr);
  for (const auto& p : {"HLT_HIMinimumBias_v2", "HLT_HIZeroBias_v4", "pprimaryVertexFilter", "pclusterCompatibilityFilter", numMinHFTowerLbl, "hiBin"})
    t->SetBranchAddress(p, &(varI[p]));
  std::map<std::string, float> varF, mcVarF;
  for (const auto& p : {"hiHF", "vz"})
    t->SetBranchAddress(p, &(varF[p]));
  t->SetBranchStatus("*", 0);
  for (const auto& p : {"run", "hiHF", "HLT_HIMinimumBias_v2", "HLT_HIZeroBias_v4", "pprimaryVertexFilter", "pclusterCompatibilityFilter", numMinHFTowerLbl, "vz", "hiBin"})
    t->SetBranchStatus(p, 1);

  std::vector<std::pair<float, bool>> values, hfdata;
  std::vector<std::pair<int, bool>> hibin;
  TH1D::SetDefaultSumw2();
  TH1F hfData1("hfData1","hf data run <= 363320", 2000,0, 10000);
  TH1F hfData2("hfData2","hf data run > 363320", 2000,0, 10000);
  TEfficiency dataEff1("dataEff1", "dataEff1", 1000, 0, 1000);
  TEfficiency dataEff2("dataEff2", "dataEff2", 1000, 0, 1000);
  std::array<std::array<size_t, 500>, 4> numMinHFTowerV{};

  auto Nevents = t->GetEntries();
  double mcYscale_data(0);
  for(Long64_t iev = 0; iev < Nevents; iev++) {
    if(iev%1000000 == 0) cout<<"Processing data event: " << iev << " / " << Nevents << endl;
    t->GetEntry(iev);
    const auto& parameter = varF.at("hiHF");
    const auto& numMinHFTower = varI.at(numMinHFTowerLbl);
    const bool pass = (varI.at("HLT_HIMinimumBias_v2")>0 && varI.at("pprimaryVertexFilter")>0 && varI.at("pclusterCompatibilityFilter")>0 && numMinHFTower>=hfCoinN);
    if (pass) {
      hfdata.push_back({parameter, (run <= 362320)});
      hibin.push_back({varI.at("hiBin"), (run <= 362320)});
      if (run <= 362320) {
        hfData1.Fill(parameter);
        if (parameter > threshold)
          values.push_back({parameter, false});
        if (parameter>threshold_norm && parameter<thresholdMax)
          mcYscale_data += 1;
      }
      else if (run > 362320)
        hfData2.Fill(parameter);
    }
    if (pass) {
      const bool passTight = (pass && varI.at(numMinHFTowerLbl)>=(hfCoinN+1));
      ((run <= 362320) ? dataEff1 : dataEff2).Fill(passTight, parameter);
    }
    if (varI.at("HLT_HIMinimumBias_v2")>0 && varI.at("pprimaryVertexFilter")>0 && varI.at("pclusterCompatibilityFilter")>0) {
      size_t idx = (parameter > 1000.)*2 + (run <= 362320);
      if ((parameter > 1000.) || (parameter < 140.))//(parameter > 25. && parameter < 60.))//(parameter > 60. && parameter < 140.))
        numMinHFTowerV[idx][numMinHFTower] += 1;
    }
  } //data events loop
  inFile.Close();

  // Compute noise study
  TH1F hfNoise1("hfNoise1","hf noise data run <= 363320", 60, 0, 60);
  TH1F hfNoise2("hfNoise2","hf noise data run > 363320", 60, 0, 60);
  std::array<std::array<double, 200>, 4> nEvt{};
  for (size_t i=0; i<101; i++)
    for (size_t j=0; j<4; j++)
      nEvt[j][i] = std::accumulate(numMinHFTowerV[j].begin()+i, numMinHFTowerV[j].end(), 0);
  for (size_t i=0; i<60; i++) {
    hfNoise1.SetBinContent(i+1, (nEvt[0][i] / nEvt[1][i])*(nEvt[3][i] / nEvt[2][i]));
    hfNoise2.SetBinContent(i+1, nEvt[1][i] / nEvt[1][20]);
  }

  // Process MC
  TFile inputMCfile("/eos/cms/store/group/phys_heavyions/anstahll/GO2023/HiForestMiniAOD_HYDJET_June7_test_SKIM.root", "READ");
  if (!inputMCfile.IsOpen()) throw std::logic_error("MC file was not found!");
  const auto& tmc = inputMCfile.Get<TTree>("hiEvtAnalyzer/HiTree");
  const auto& tskimanalysismc = inputMCfile.Get<TTree>("skimanalysis/HltTree");
  const auto& thltanalysismc = inputMCfile.Get<TTree>("hltanalysis/HltTree");
  tmc->AddFriend(tskimanalysismc);
  tmc->AddFriend(thltanalysismc);

  for (const auto& p : {"HLT_HIMinimumBias_v2", "pprimaryVertexFilter", "pclusterCompatibilityFilter", numMinHFTowerLbl})
    tmc->SetBranchAddress(p, &(mcVarI[p]));
  for (const auto& p : {"hiHF", "vz"})
    tmc->SetBranchAddress(p, &(mcVarF[p]));
  tmc->SetBranchStatus("*", 0);
  for (const auto& p : {"hiHF", "HLT_HIMinimumBias_v2", "pprimaryVertexFilter", "pclusterCompatibilityFilter", "numMinHFTower4", "numMinHFTower5", "vz"})
    tmc->SetBranchStatus(p, 1);

  TH1F hfMc1("hfMc1","hf mc", 2000, 0, 10000);
  TH1F hfMc2("hfMc2","hf mc", 2000, 0, 10000);
  TEfficiency mcEff("mcEff", "mcEff", 1000, 0, 1000);

  Nevents = tmc->GetEntries();
  double mcYscale_mc(0);
  for(Long64_t iev = 0; iev < Nevents; iev++) {
    if(iev%500000 == 0) cout<<"Processing mc event: " << iev << " / " << Nevents << endl;
    tmc->GetEntry(iev);
    const auto parameter = mcVarF.at("hiHF") * mcXscale;
    const bool pass = (mcVarI.at("HLT_HIMinimumBias_v2")>0 && mcVarI.at("pprimaryVertexFilter")>0 && mcVarI.at("pclusterCompatibilityFilter")>0 && mcVarI.at(numMinHFTowerLbl)>=hfCoinN);
    if (pass) {
      hfMc1.Fill(parameter);
      if (parameter>threshold_norm && parameter<thresholdMax)
        mcYscale_mc += 1;
    }
    if (parameter <= threshold)
      values.push_back({parameter, true});
    hfMc2.Fill(parameter);
    mcEff.Fill(pass, parameter);
  } //end of mc loop
  inputMCfile.Close();

  TNtuple nt("nt","","value");
  TH1F hfCombined("hfCombined","hf mc", 2000,0, 10000);

  // Scale MC
  const auto mcYscale = mcYscale_data / mcYscale_mc;
  std::cout<<"[INFO] mcYscale = "<< mcYscale << std::endl;
  hfMc1.Scale(mcYscale);
  hfMc2.Scale(mcYscale);
  for (const auto& v : values) {
    nt.Fill(v.first, v.second ? mcYscale : 1.);
    hfCombined.Fill(v.first, v.second ? mcYscale : 1.);
  }
  const auto totEff = hfData1.Integral() / hfCombined.Integral();

  // Store histograms
  const std::string outputTag = Form("v1_TestRun_xSF0p86_ySFwithAllFilter_%dThE%dGeV_officiaMC2022__Threshold%.0f__NOMINAL__Normalisation%.0f_%.0f__GT_Aug10_NEW", hfCoinN, hfCoinThr, threshold, threshold_norm, thresholdMax);
  TFile outFile(Form("CentralityTable_HFtowers200_DataPbPb_usingMC_%s.root", outputTag.c_str()),"recreate");
  const auto& dir = outFile.mkdir(tag.c_str());
  dir->cd();
  nt.Write();
  hfData1.Write();
  hfData2.Write();
  hfMc1.Write();
  hfMc2.Write();
  mcEff.Write();
  dataEff1.Write();
  dataEff2.Write();
  hfNoise1.Write();
  hfNoise2.Write();
  hfCombined.Write();
  outFile.Close();

  const auto passed = values.size();
  double totalXsec(0);
  for(const auto& v : values)
    totalXsec += v.second ? mcYscale : 1.;

  std::cout << std::endl;
  std::cout << "Selected events = " << passed << std::endl;
  std::cout << "Selected weighed events = " << totalXsec << std::endl;
  std::cout << "Total efficiency = " << totEff << std::endl;

  // Create text file
  ofstream txtfile(Form("output_DataPbPb_usingMC_hiHF_%s.txt", outputTag.c_str()));
  txtfile << "Input tree: " << inFileName << endl;
  txtfile << "Tag name: " << tag << endl;
  txtfile << "Number of events = " << passed << std::endl;
  txtfile << "Number of weighted events = " << totalXsec << std::endl;
  txtfile << "Total efficiency = " << totEff << std::endl;
  txtfile << std::endl;
  txtfile << "-------------------------------------" << std::endl;
  txtfile << "Using MC to correct for peripheral events. Threshold = " << threshold<<". mc X scale factor = "<< mcXscale<< std::endl;
  txtfile << std::endl;
  txtfile << "-------------------------------------" << std::endl;
  txtfile << "hiHF based cuts are: " << std::endl;
  txtfile << "(";

  // Store bin boundaries
  const auto size = values.size();
  std::sort(values.begin(), values.end());
  std::vector<double> binboundaries(nbins+1);
  binboundaries[0] = 0.;
  binboundaries[nbins] = values[size-1].first;
  std::cout << "Events per bin = " << totalXsec/nbins << std::endl;
  double integral(0);
  size_t currentbin(1);
  for(const auto& v : values) {
    const auto& val = v.first;
    integral += v.second ? mcYscale : 1.;
    const auto sum = currentbin*(totalXsec/nbins);
	if(integral > sum) {
      //std::cout << "current bin = " << currentbin << " ; integral = " << integral << " ; sum = " << sum << std::endl;
	  binboundaries[currentbin] = val >= 0 ? val : 0;
	  currentbin++;
	}
  }
  std::cout << currentbin << std::endl;
  for(size_t i = 0; i < nbins; i++)
    txtfile << binboundaries[i] << ", ";
  txtfile << binboundaries[nbins] << ")" << std::endl;
  txtfile << std::endl;
  txtfile<<"-------------------------------------"<<endl;
  txtfile.close();

  // Check bin boundaries
  int newbin, oldbin;
  TFile outf(Form("compare_centralitybins_%s.root", outputTag.c_str()),"recreate");
  TTree t1("anaCentrality_Before362320","analysis level centrality");
  TTree t2("anaCentrality_After362320","analysis level centrality");
  for (auto& t : {&t1, &t2}) {
    t->Branch("newBin",&newbin,"newBin/I");
    t->Branch("oldBin",&oldbin,"oldBin/I");
  }
  for (size_t i=0; i<hfdata.size(); i++) {
    newbin = 199;
    for(size_t b = 0; b < 200; ++b){
      if(hfdata[i].first >= binboundaries[199-b]){
        newbin = b;
        break;
      }
    }
    oldbin = hibin[i].first;
    (hfdata[i].second ? t1 : t2).Fill();
  }
  t1.Write();
  t2.Write();
  outf.Close();
}
