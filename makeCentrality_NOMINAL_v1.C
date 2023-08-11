#include "TFile.h"
#include "TTree.h"
#include <iostream>
using namespace std;

void makeCentrality_NOMINAL_v1(){

  TFile* infData = TFile::Open("/eos/cms/store/group/phys_heavyions/anstahll/GO2023/HiForestAODZB_TestRun_Data_ZB1.root");

  TFile* outf = new TFile("compare_centralitybins_Th100_Withofficial2022MC_NOMINAL_Aug11_v1.root","recreate");
  //Threshold = 100 (Nominal)
  float bounds[201] = {0, 13.3145, 14.4357, 15.3788, 16.266, 17.0986, 17.952, 18.8266, 19.7014, 20.5792, 21.5344, 22.4929, 23.4772, 24.5017, 25.572, 26.6768, 27.8178, 29.0134, 30.2532, 31.5644, 32.9105, 34.305, 35.7672, 37.275, 38.8872, 40.5392, 42.3589, 44.1401, 46.0941, 48.1159, 50.1956, 52.3762, 54.691, 57.0964, 59.5681, 62.1371, 64.8221, 67.6734, 70.5275, 73.4911, 76.5685, 79.7387, 83.1046, 86.5429, 90.0977, 93.7198, 97.6043, 101.455, 105.557, 109.343, 113.422, 117.411, 121.893, 126.179, 130.787,135.444, 140.796, 145.707, 151.592, 157.613, 163.95, 169.501, 176.017, 182.681, 189.878, 196.804, 203.462, 211.277, 219.342, 227.236, 235.32, 243.558, 252.57, 262.022, 271.51, 280.539, 288.88, 298.446, 309.011, 318.898, 329.053, 340.171, 350.878, 362.12, 373.652, 385.575, 397.356, 409.668, 423.177, 436.935, 450.249, 463.773, 478.832, 494.756, 510.134, 524.757, 540.658, 557.165, 572.879, 590.542, 607.877, 626.532, 643.586, 662.722, 678.161, 699.226, 719.336, 741.804, 763.89, 784.438, 803.974, 826.603, 847.498, 868.371, 889.216, 912.972, 940.382, 965.47, 992.524, 1014.86, 1037.85, 1062.6, 1088.12, 1114.68, 1142.1, 1171.45, 1196.72, 1223.67, 1254.69, 1285.67, 1315.25, 1344.15, 1373.36, 1402.61, 1436.13, 1468.57, 1503.51, 1536.4, 1573.15, 1605.26, 1637.94, 1669.47, 1704.82, 1742.03, 1776.44, 1816.34, 1852.98, 1895.98, 1936.18, 1973.47, 2017.17, 2053.69, 2095.91, 2141.64, 2187.64, 2235.13, 2286.23, 2330.54, 2369.8, 2412.87, 2460, 2514.63, 2566.73, 2617.75,2667.91, 2723.12, 2772.9, 2828.47, 2883.42, 2938.53, 2997.39, 3056.38, 3112.22, 3169.48, 3227.21, 3297.18, 3368.91, 3429.5, 3491.42, 3552.9, 3625.8, 3693.3, 3764.38, 3836.9, 3918, 4002.25, 4074.25, 4152.45, 4231.82, 4314.37, 4402.17, 4486.14, 4573.65, 4674.61, 4776.59, 4871.57, 4972.08, 5066.03, 5181.6, 5340.05, 5920.4};

  float hf; 
  int hiBin, numMinHFTower4;
  int phfCoincFilter3, pprimaryVertexFilter, pclusterCompatibilityFilter;
  int hlt_mb;

  TTree* tref = (TTree*)infData->Get("hiEvtAnalyzer/HiTree");
  TTree *tskim = (TTree*)infData->Get("skimanalysis/HltTree");
  TTree *thlt = (TTree*)infData->Get("hltanalysis/HltTree");

  tref->AddFriend(tskim);
  tref->AddFriend(thlt);

  tref->SetBranchAddress("hiHF", &hf);
  tref->SetBranchAddress("hiBin", &hiBin);
  tref->SetBranchAddress("HLT_HIMinimumBias_v2", &hlt_mb);
  tref->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter);
  tref->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);
  //tskim->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3);
  tref->SetBranchAddress("numMinHFTower4", &numMinHFTower4);

  tref->SetBranchStatus("*", 0);
  for (const auto& p : {"hiHF", "hiBin", "HLT_HIMinimumBias_v2", "pprimaryVertexFilter", "pclusterCompatibilityFilter", "numMinHFTower4"})
    tref->SetBranchStatus(p, 1);

  int newbin;
  TTree* t = new TTree("anaCentrality","analysis level centrality");
  t->Branch("newBin",&newbin,"newBin/I");
  t->Branch("oldBin",&hiBin,"newBin/I");

  int N = tref->GetEntries();
  for(int i = 0; i < N; ++i){
    tref->GetEntry(i);
    //thlt->GetEntry(i);
    //tskim->GetEntry(i);
    
    if(i % 500000 == 0) cout<<"processing event : "<<i<<endl;

    if(hlt_mb != 1) continue;
    if(pprimaryVertexFilter != 1) continue;
    if(pclusterCompatibilityFilter != 1) continue;
    if(numMinHFTower4 < 2)continue;

    newbin = 199; 
    for(int b = 0; b < 200; ++b){
      if(hf >= bounds[199-b]){
	newbin = b;
	break;	    
      }
    }

    t->Fill();
  }

  t->Write();
  outf->Write();

}


