#include <TChain.h>
#include <ROOT/RDataFrame.hxx>
#include <iostream>
#include <memory>

void skimHIForest_data(const std::string& dirin="/eos/cms/store/group/phys_heavyions/anstahll/GO2023/", const std::string& dirout="/eos/cms/store/group/phys_heavyions/anstahll/GO2023/SKIM/")
{
  // prepare multi-threading
  ROOT::EnableImplicitMT(14);

  const std::map<std::string, std::string> treeInfo({{"hiEvtAnalyzer", "HiTree"}, {"skimanalysis", "HltTree"}, {"hltanalysis", "HltTree"}});
  const std::vector<std::string> treeDirV({"hiEvtAnalyzer", "skimanalysis", "hltanalysis"});
  std::map<size_t, std::vector<std::string>> varsUI({ {0, {"run"}} });
  std::map<size_t, std::vector<std::string>> varsI({ {0, {"numMinHFTower4", "numMinHFTower5", "hiBin"}},
                                                     {1, {"pprimaryVertexFilter", "pclusterCompatibilityFilter"}},
                                                     {2, {"HLT_HIMinimumBias_v2"}} });
  std::map<size_t, std::vector<std::string>> varsF({ {0, {"hiHF", "vz"}} });
  std::vector<std::string> columns;
  for (const auto& c : varsUI)
    for (const auto& v : c.second)
      columns.push_back(v);
  for (const auto& c : varsI)
    for (const auto& v : c.second)
      columns.push_back(v);
  for (const auto& c : varsF)
    for (const auto& v : c.second)
      columns.push_back(v);

  for (size_t iF=0; iF<24; iF++) {
    const auto filein = dirin + Form("HIForest_HIMinimumBias%lu_HIRun2022A_20230811.root", iF);
    auto fileout = dirout + filein.substr(filein.rfind("/")+1);
    fileout = fileout.substr(0, fileout.rfind(".root")) + "_SKIM.root";

    std::cout << "Skimming: " << filein << std::endl; 

    // Open trees
    std::vector<std::unique_ptr<TChain>> treeV;
    for (const auto& c : {"hiEvtAnalyzer/HiTree", "skimanalysis/HltTree", "hltanalysis/HltTree"}) {
      treeV.emplace_back(new TChain(c));
      treeV.back()->AddFile(filein.c_str());
    }
    auto& tree = *treeV[0];
    tree.AddFriend(treeV[1].get());
    tree.AddFriend(treeV[2].get());
    tree.SetBranchStatus("*", 0);
    for (const auto& c : columns)
      tree.SetBranchStatus(c.c_str(), 1);

    // Skim trees
    ROOT::RDataFrame df(tree, columns);
    df.Filter("HLT_HIMinimumBias_v2 == 1 && pprimaryVertexFilter == 1").Snapshot("skim", dirout+"skim.root", columns);

    std::cout << "Changing to HiForest format: " << filein << std::endl;

    // Reformat the HiForest
    TFile fskim((dirout+"skim.root").c_str(), "READ");
    const auto& tskim = fskim.Get<TTree>("skim");
    treeV[0].reset(new TChain("hiEvtAnalyzer/HiTree"));
    treeV[0]->AddFile(filein.c_str());

    std::map<std::string, TTree*> newTreeM;
    TFile fout(fileout.c_str(), "RECREATE");
    for (size_t iD=0; iD<3; iD++) {
      fout.cd();
      const auto& tdir = fout.mkdir(treeDirV[iD].c_str());
      tdir->cd();
      treeV[iD]->SetBranchStatus("*", 0);
      for (const auto& c : columns)
        if (treeV[iD]->GetBranch(c.c_str()) != nullptr)
          treeV[iD]->SetBranchStatus(c.c_str(), 1);
      newTreeM[treeDirV[iD]] = treeV[iD]->CloneTree(0);
    }

    std::map<std::string, UInt_t> oldUI, newUI;
    std::map<std::string, int> oldI, newI;
    std::map<std::string, float> oldF, newF;
    for (size_t iD=0; iD<3; iD++) {
      for (const auto& v : varsUI[iD]) {
        tskim->SetBranchAddress(v.c_str(), &(oldUI[v]));
        treeV[iD]->SetBranchAddress(v.c_str(), &(newUI[v]));
      }
      for (const auto& v : varsI[iD]) {
        tskim->SetBranchAddress(v.c_str(), &(oldI[v]));
        treeV[iD]->SetBranchAddress(v.c_str(), &(newI[v]));
      }
      for (const auto& v : varsF[iD]) {
        tskim->SetBranchAddress(v.c_str(), &(oldF[v]));
        treeV[iD]->SetBranchAddress(v.c_str(), &(newF[v]));
      }
    }

    for (Long64_t iEntry=0; iEntry<tskim->GetEntries(); iEntry++) {
      tskim->GetEntry(iEntry);
      for (const auto& v : oldUI)
        newUI.at(v.first) = v.second;
      for (const auto& v : oldI)
        newI.at(v.first) = v.second;
      for (const auto& v : oldF)
        newF.at(v.first) = v.second;
      for (auto& t : newTreeM)
        t.second->Fill();
    }

    // Store trees
    fout.Write();
    fout.Close();
    fskim.Close();
  }
}
