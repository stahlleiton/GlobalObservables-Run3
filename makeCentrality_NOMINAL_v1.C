#include "TFile.h"
#include "TTree.h"
#include <iostream>

void makeCentrality_NOMINAL_v1()
{
  // Constant parameters
  const int  hfCoinThr = 4;

  //Threshold = 100 (Nominal)
  const float bounds[201] = {0, 10.3936, 11.139, 11.8473, 12.5633, 13.2989, 14.0502, 14.8044, 15.5867, 16.3663, 17.17, 18.0086, 18.8808, 19.7601, 20.6558, 21.6331, 22.6134, 23.6261, 24.6832, 25.7968, 26.9304, 28.1161, 29.36, 30.6509, 32.0071, 33.4047, 34.8694, 36.3867, 37.9862, 39.6281, 41.4218, 43.286, 45.2203, 47.236, 49.3532, 51.5665, 53.8534, 56.3585, 58.8423, 61.458, 64.1841, 67.0791, 70.0369, 73.0162, 76.1905, 79.4435, 82.8564, 86.3925, 90.0308, 93.7451, 97.7624, 101.712, 105.745, 110.025, 114.499, 119.139, 124.059, 129.035, 134.178, 139.549, 145.126, 150.904, 156.839, 162.91, 169.239, 175.622, 182.399, 189.278, 196.452, 203.883, 211.603, 219.357, 227.408, 235.656, 244.272, 252.903, 261.808, 271.378, 280.853, 290.514, 300.82, 311.1, 322.012, 333.029, 344.39, 356.01, 367.987, 380.265, 392.762, 405.627, 419.015, 432.48, 446.543, 460.754, 475.059, 490.29, 505.592, 521.326, 537.077, 553.2, 570.138, 587.503, 605.092, 622.876, 640.819, 658.965, 677.708, 697.443, 717.242, 737.287, 757.71, 778.643, 799.681, 821.531, 843.182, 865.345, 888.858, 911.826, 935.783, 960.571, 985.701, 1010.95, 1036.59, 1063.06, 1089.78, 1116.68, 1144.23, 1172.47, 1201, 1229.83, 1259.37, 1290.21, 1321.26, 1352.78, 1383.84, 1416.33, 1448.65, 1481.42, 1515.49, 1549.22, 1584.66, 1621.11, 1657.5, 1694.54, 1732.75, 1771.03, 1809.55, 1849.23, 1888.76, 1928.87, 1970.08, 2012.39, 2055.13, 2098.15, 2142.24, 2185.82, 2231.89, 2277.81, 2324.66, 2372.53, 2421.23, 2470.56, 2520.4, 2572.21, 2623.61, 2675.65, 2729.34, 2784.23, 2838.71, 2895.29, 2953.37, 3011.77, 3070.92, 3131.24, 3192.27, 3255.13, 3319.95, 3385.98, 3452.52, 3520.57, 3592.11, 3663.32, 3735.27, 3810, 3886.55, 3965.05, 4045.09, 4126.84, 4209.26, 4295.86, 4384.6, 4472.23, 4566.55, 4662.89, 4760.9, 4864.15, 4965.96, 5072.55, 5190.35, 5337.33, 8031.6};

  // Process data
  const std::string inFileName = "/eos/cms/store/group/phys_heavyions/anstahll/GO2023/HIForest_HIMinimumBias_HIRun2022A_20230811_SKIM.root";
  TFile inFile(inFileName.c_str(), "READ");
  if (!inFile.IsOpen()) throw std::logic_error("Data file was not found!");
  const auto& tref = inFile.Get<TTree>("hiEvtAnalyzer/HiTree");
  const auto& tskimanalysis = inFile.Get<TTree>("skimanalysis/HltTree");
  const auto& thltanalysis = inFile.Get<TTree>("hltanalysis/HltTree");
  tref->AddFriend(tskimanalysis);
  tref->AddFriend(thltanalysis);

  UInt_t run;
  tref->SetBranchAddress("run", &run);
  std::map<std::string, int> varI;
  const char* numMinHFTowerLbl = Form("numMinHFTower%d", hfCoinThr);
  for (const auto& p : {"HLT_HIMinimumBias_v2", "pprimaryVertexFilter", "pclusterCompatibilityFilter", numMinHFTowerLbl, "hiBin"})
    tref->SetBranchAddress(p, &(varI[p]));
  std::map<std::string, float> varF;
  for (const auto& p : {"hiHF", "vz"})
    tref->SetBranchAddress(p, &(varF[p]));
  tref->SetBranchStatus("*", 0);
  for (const auto& p : {"run", "hiHF", "HLT_HIMinimumBias_v2", "pprimaryVertexFilter", "pclusterCompatibilityFilter", numMinHFTowerLbl, "vz", "hiBin"})
    tref->SetBranchStatus(p, 1);

  int newbin;
  TFile outf("compare_centralitybins_Th100_Withofficial2022MC_NOMINAL_Aug11_v1.root","recreate");
  TTree t1("anaCentrality_Before362320","analysis level centrality");
  TTree t2("anaCentrality_After362320","analysis level centrality");
  for (auto& t : {&t1, &t2}) {
    t->Branch("newBin",&newbin,"newBin/I");
    t->Branch("oldBin",&(varI.at("hiBin")),"oldBin/I");
  }

  const auto& Nevents = tref->GetEntries();
  for(Long64_t iev = 0; iev < Nevents; iev++) {
    if(iev%1000000 == 0) cout<<"Processing data event: " << iev << " / " << Nevents << endl;
    tref->GetEntry(iev);
    const bool pass = (varI.at("HLT_HIMinimumBias_v2")>0 && varI.at("pprimaryVertexFilter")>0 && varI.at("pclusterCompatibilityFilter")>0 && varI.at(numMinHFTowerLbl)>=2);
    if (!pass) continue;
    newbin = 199;
    for(size_t b = 0; b < 200; ++b){
      if(varF.at("hiHF") >= bounds[199-b]){
        newbin = b;
        break;
      }
    }
    ((run <= 362320) ? t1 : t2).Fill();
  }
  t1.Write();
  t2.Write();
  outf.Close();
}
