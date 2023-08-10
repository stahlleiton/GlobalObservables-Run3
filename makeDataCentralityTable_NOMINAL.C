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
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <map>

#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"

TDatime* date = new TDatime();
using namespace std;

bool descend(float i,float j) { return (i>j); }

double dEffErrF(double* x, double* par) // Use "err"
{
  // Error function
  
  return 0.5*(1+TMath::Erf((x[0]-par[0])/(TMath::Sqrt(x[0])*par[1])));
}

double dEffErrFPW(double* x, double* par) // Use "errPW"
{
  // Piecewise function changing at par[0]
  // Note that in this function par[2] and par[3] must be correlated in this way: par[3] = par[2]/par[0], to ensure the continuity
  
  if (x[0] <= par[0]) return 0.5*(1+TMath::Erf((x[0]-par[1])/par[2]));
  else return 0.5*(1+TMath::Erf((x[0]-par[1])/(par[3]*x[0])));
}

double dEffStepF(double* x, double* par) // Use "step"
{
  // Step function changing to 1
  // par[0] is the efficiency for x <= par[0]
  
  if (x[0] <= par[0]) return par[1];
  else return 1.0;
}

TF1* fEffStepF = new TF1("hfEfficiency",dEffStepF,0.,5000.,2);
TF1* fEffErrF = new TF1("hfEfficiency",dEffErrF,0.,5000.,2);
TF1* fEffErrFPW = new TF1("hfEfficiency",dEffErrFPW,0.,5000.,3);

const std::map< std::string , TF1* > effFfunc = {
  {"step",fEffStepF},
  {"err",fEffErrF},
  {"errPW",fEffErrFPW}
};

const std::map< std::string , std::vector<double> > effFpars = {
  {"step",{25.0,0.877}},
  {"err",{13.439,0.811}},
  //{"errPW",{11.36,13.5298,2.9839,0.2626}}
  {"errPW",{25.0,22.6610,8.0034,0.3201}} //for hiHF
  //{"errPW",{25.0,14.4182,7.6966,0.3079}} //for hiHFECut
};

const std::map< std::string , std::string > effFstring = {
  {"step","if (x <= par[0]) eff = par[1] ; else eff= 1.0"},
  {"err","eff = 0.5*(1+TMath::Erf((x-par[0])/(TMath::Sqrt(x)*par[1])))"},
  {"errPW","if (x <= par[0]) eff = 0.5*(1+TMath::Erf((x-par[1])/par[2])) ; else eff = 0.5*(1+TMath::Erf((x-par[1])/(par[3]*x)))"}
};

void makeDataCentralityTable_NOMINAL(int nbins = 200, const string label = "hiHF",
                             const char * tag = "CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1205x02_offline",
                             bool useEffFunc = false, const char* effName = "errPW", bool useMC = true, double intEff = -1., 
			     bool geomInfo = false){

  if (useEffFunc && useMC) { cout <<"[ERROR] can't use both eff function and MC methods"; return; }
  
  TH1D::SetDefaultSumw2();

  const char* outputTag = "v1_TestRun_xSF0p86_ySFwithAllFilter_ThE4GeV_officiaMC2022__Threshold100__NOMINAL__Normalisation100_5000__GT_Aug10_NEW";
  TFile *outFile = new TFile(Form("CentralityTable_HFtowers200_DataPbPb_usingMC_d%d_%s.root",date->GetDate(),outputTag),"recreate");

 
  //TString inFileName = "/eos/cms/store/group/phys_heavyions/nsaha/GO2023/HiForest_TestRun_data/HITestRaw/HiForestAODZB_TestRun_Data_ZB9_merged.root";
  TString inFileName = "/eos/cms/store/group/phys_heavyions/nsaha/GO2023/HiForest_TestRun_data_originalGT/HITestRaw1/HiForest_TestRun_data_ZB1_GT/230809_194339/0000/HiForestAODZB_TestRun_Data_ZB1_merged_small.root";
  //TString inFileName = "/eos/cms/store/group/phys_heavyions/nsaha/GO2023/HiForest_TestRun_data/HITestRaw/HiForestAODZB_TestRun_Data_ZB1_merged_small.root";

  TFile *inFile = TFile::Open(inFileName);
  TTree *t = (TTree*)inFile->Get("hiEvtAnalyzer/HiTree");
  TTree * tskimanalysis = (TTree*)inFile->Get("skimanalysis/HltTree");
  TTree * thltanalysis = (TTree*)inFile->Get("hltanalysis/HltTree");

  t->AddFriend(tskimanalysis);
  t->AddFriend(thltanalysis);
  
  
  ofstream txtfile(Form("output_DataPbPb_%s_d%d_%s_%s.txt", useEffFunc?effName:(useMC?"usingMC":"intEff"),date->GetDate(), label.data(),outputTag));
  txtfile << "Input tree: " << inFileName << endl;
  txtfile << "Tag name: " << tag << endl;
  
  TDirectory *dir = outFile->mkdir(tag);
  dir->cd();
  TNtuple * nt = new TNtuple("nt","","value");
  TH1F* hfData = new TH1F("hfData","hf data", 100,0, 6000);
  TH1F* hfMc = new TH1F("hfMc","hf mc", 100,0, 6000);
  TH1F* hfCombined = new TH1F("hfCombined","hf mc", 100,0, 6000);

  const int runNum = 1;
  CentralityBins * bins = new CentralityBins(Form("run%d",runNum), tag, nbins);
  bins->table_.reserve(nbins);
  
  // Determine the binning to be used from input
  bool binHF = label.compare("hiHF") == 0;
  bool binNtrks = label.compare("hiNtracks") == 0;
  
  // Deffine efficiency to weight events
  TF1* fEff(0x0);
  TParameter<double>* gEff(0x0);
  double effPrime = -1;
  if (useEffFunc)
  {
    fEff = effFfunc.at(effName);
    if (!fEff)
    {
      cout << "[ERROR] No efficiency function could be defined " << endl;
      return;
    }
    else cout << "[INFO] Using efficiency function ..." << endl;

    int parS = effFpars.at(effName).size();
    for (int i = 0 ; i < parS ; i++)
    {
      fEff->FixParameter(i,effFpars.at(effName).at(i));
    }
    
    if (!strcmp(effName,"step"))
    {
      cout << "Computing average efficiency for step function" << endl;
      
      const char* selection = "pprimaryVertexFilter && pclusterCompatibilityFilter && phfCoincFilter2Th4";
      
      double ntot = (1./effFpars.at("step").at(1))*t->GetEntries(selection);
      double nlow = t->GetEntries(Form("%s && %s<=%f",selection,label.data(),effFpars.at("step").at(0)));
      double nhigh = t->GetEntries(Form("%s && %s>%f",selection,label.data(),effFpars.at("step").at(0)));
      effPrime = nlow / (ntot - nhigh);
      
      cout << "Average efficiency for  x <= " << effFpars.at("step").at(0) << " = " << effPrime << endl;
      
      fEff->FixParameter(1,effPrime);
    }
  }
  else if (!useMC)
  {
    if (intEff>0 && intEff<=1) gEff = new TParameter<double>("eff",intEff);
    else
    {
      cout << "[ERROR] Input integrated efficiency is not correct, must be (0,1]" << fEff << endl;
      return;
    }
    cout << "Using global efficiency of " << gEff->GetVal() << endl;
  }
  
  
  //Here we need the default Glauber for 2.76 or 5 TeV
  //TFile * inputMCfile = TFile::Open("/eos/cms/store/group/phys_heavyions/nsaha/GO2023/HiForest_HYDJET_official/MinBias_Drum5F_5p36TeV_hydjet/HiForest_TestRun2022_MC_official/230803_133309/0000/HiForestMiniAOD_OfficialMC_Hyd2022_merged.root");

  TFile * inputMCfile = TFile::Open("/eos/cms/store/group/phys_heavyions/nsaha/GO2023/HiForest_HYDJET_official_modifiedGT/MinBias_Drum5F_5p36TeV_hydjet/HiForest_TestRun2022_MC_official_GT/230809_191838/0000/HiForestMiniAOD_OfficialMC_Hyd2022_merged.root");

  TChain * tmc = new TChain("hiEvtAnalyzer/HiTree","");
  TChain * tskimanalysismc = new TChain("skimanalysis/HltTree","");
  TChain * thltanalysismc = new TChain("hltanalysis/HltTree","");

  tmc->Add(inputMCfile->GetName());
  tskimanalysismc->Add(inputMCfile->GetName());
  thltanalysismc->Add(inputMCfile->GetName());

  tmc->AddFriend(tskimanalysismc);
  tmc->AddFriend(thltanalysismc);
  
  double binboundaries[nbins+1];
  vector<float> values;

  UInt_t run;
  float hf, hfecut, hfplus, hfpluseta4, hfminuseta4, hfminus, hfhit, ee, eb, zdc, zdcplus, zdcminus;
  int lumi, npix, npixtrks, ntrks ,pprimaryVertexFilter, pclusterCompatibilityFilter, phfCoincFilter2Th4, pclusterfilter, pvfilter, HLT_HIMinimumBias_v2, numMinHFTower4, numMinHFTower5;
  t->SetBranchAddress("run",    &run);
  //t->SetBranchAddress("lumi",   &lumi);
  t->SetBranchAddress("hiHF",           &hf);
  
  t->SetBranchAddress("HLT_HIMinimumBias_v2", &HLT_HIMinimumBias_v2);
  t->SetBranchAddress("pprimaryVertexFilterHI", &pvfilter);
  t->SetBranchAddress("pclusterCompatibilityFilter", &pclusterfilter);
  t->SetBranchAddress("numMinHFTower4",&numMinHFTower4);
  //t->SetBranchAddress("numMinHFTower5",&numMinHFTower5);

  

  float hf_mc, hfecut_mc;
  int ntrks_mc ,pprimaryVertexFilter_mc, pclusterCompatibilityFilter_mc, phfCoincFilter2Th4_mc, numMinHFTower4_mc, numMinHFTower5_mc, HLT_HIMinimumBias_v2_mc;
  tmc->SetBranchAddress("hiHF", &hf_mc);
  tmc->SetBranchAddress("HLT_HIMinimumBias_v2", &HLT_HIMinimumBias_v2_mc);
  tmc->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter_mc);
  tmc->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter_mc);
  tmc->SetBranchAddress("numMinHFTower4",&numMinHFTower4_mc);
  //tmc->SetBranchAddress("numMinHFTower5",&numMinHFTower5_mc);

  unsigned int Nevents = t->GetEntries();
  //  std::cout<<"Selected Nevents ="<<Nevents<<std::endl;

  double threshold = 0.0;
  double threshold_norm = 0.0;
  double thresholdMax = 999999.0;
  if (useMC) { threshold = 100.0; threshold_norm = 1000; thresholdMax = 5000;}
  double mcXscale = 0.86;
  double mcYscale = 1;
  if (useMC)
    mcYscale = 1.0*t->GetEntries(Form("hiHF>%f && hiHF<%f && run <= 362320 && HLT_HIMinimumBias_v2 ==1 && pprimaryVertexFilterHI ==1 && pclusterCompatibilityFilter ==1 && numMinHFTower4 >=2", threshold, thresholdMax))/tmc->GetEntries(Form("hiHF>%f && hiHF<%f && HLT_HIMinimumBias_v2 ==1 && pprimaryVertexFilter ==1 && pclusterCompatibilityFilter ==1 && numMinHFTower4 >=2", threshold/mcXscale, thresholdMax/mcXscale));

    //if (useMC)
    //mcYscale = 1.0*data_entries/tmc->GetEntries(Form("hiHF>%f && hiHF<%f", threshold/mcXscale, thresholdMax/mcXscale));
  
  cout<<"[INFO] mcYscale = "<< mcYscale <<endl;
  double totalXsec(0);
  double totalXsec_data(0);
  double totalXsec_mc(0);
  int passed(0);
  for(unsigned int iev = 0; iev < Nevents; iev++) {
    
    if(iev%1000000 == 0) cout<<"Processing data event: " << iev << " / " << Nevents << endl;
    t->GetEntry(iev);
    
    if (run <= 362320 && HLT_HIMinimumBias_v2 ==1 && pvfilter ==1 && pclusterfilter ==1 && numMinHFTower4 >=2 )
    {
      float parameter = -1;
      if(binHF) parameter = hf;
      //if(binHFECut) parameter = hfecut;
      //if(binNtrks) parameter = ntrks;
    
	hfData->Fill(parameter);
      if (parameter > threshold) {
	passed++;
	values.push_back(parameter);
	nt->Fill(parameter);
	hfCombined->Fill(parameter, 1.0);
      }
      if (useEffFunc)
	{
	  double eff = fEff->Eval(parameter);
	  //cout << "eff(" << parameter << ") = " << eff << endl;
	  if (eff <= 1 && eff > 0.) totalXsec_data += 1./eff;
	  else totalXsec_data += 1.;
	} else if (parameter > threshold){
	totalXsec_data+=1.;
      }
    }
  } //data events loop
  
  if (!useEffFunc && !useMC)
    totalXsec_data = totalXsec_data/(gEff->GetVal());
  
  if (useMC){
    Nevents = tmc->GetEntries();
    //Nevents = 20000;
    for(unsigned int iev = 0; iev < Nevents; iev++) {
      
      if(iev%5000 == 0) cout<<"Processing mc event: " << iev << " / " << Nevents << endl;
      tmc->GetEntry(iev);

      if (HLT_HIMinimumBias_v2_mc ==1 && pprimaryVertexFilter_mc ==1 && pclusterCompatibilityFilter_mc ==1 && numMinHFTower4_mc >=2 )
      {
	  float parameter = -1;
	  if(binHF) parameter = hf_mc*mcXscale;
	  //if(binHFECut) parameter = hfecut_mc*mcXscale;
	  //if(binNtrks) parameter = ntrks_mc*mcXscale;
	  
	  hfMc->Fill(parameter, mcYscale);
	  
	  if (parameter <= threshold){
	    passed++;
	    values.push_back(parameter);
	    nt->Fill(parameter, mcYscale);
	    hfCombined->Fill(parameter, mcYscale);
	    totalXsec_mc+=1.;
	  }
      }
    } //end of mc loop
  } //end of if mc
  
  totalXsec_mc = totalXsec_mc*mcYscale;
  totalXsec = totalXsec_data + totalXsec_mc;

  cout << endl;
  cout << "Selected events = " << passed << endl;
  cout << "Selected weighed events = " << totalXsec << endl;
  if (useMC)
    cout << "Total efficiency = " << (double) 1.0*hfData->Integral()/hfCombined->Integral() << endl;
  else 
    cout << "Total efficiency = " << (double)(passed/totalXsec) << endl;
  
  sort(values.begin(),values.end());
  
  txtfile << "Number of events = " << passed << endl;
  txtfile << "Number of weighted events = " << totalXsec << endl;
  if (useMC)
    txtfile << "Total efficiency = " << (double) 1.0*hfData->Integral()/hfCombined->Integral()<< endl;
  else
    txtfile << "Total efficiency = " << (double)(passed/totalXsec) << endl;
  txtfile << endl;
  txtfile << "-------------------------------------" << endl;
  if (useEffFunc)
  {
    txtfile << "EFficiency shape is " << effName << " :" << endl;
    txtfile << effFstring.at(effName) << endl;
    txtfile << "with pars : " << endl;
    for (int i = 0 ; i < effFpars.at(effName).size() ; i++)
    {
      txtfile << Form("par[%d] = %f",i,fEff->GetParameter(i)) << endl;
    }
    
    if (effPrime>0) txtfile << "Average efficiency for  x <= " << effFpars.at("step").at(0) << " = " << effPrime << endl;
  }
  if (useMC) { 
    txtfile << "Using MC to correct for peripheral events. Threshold = " << threshold<<". mc X scale factor = "<< mcXscale<< endl;
  }
  txtfile << endl;
  txtfile << "-------------------------------------" << endl;
  txtfile << label.data() << " based cuts are: " << endl;
  txtfile << "(";
  
  unsigned int size = values.size();
  binboundaries[nbins] = values[size-1];
  cout << "Events per bin = " << totalXsec/(double)(nbins) << endl;
  if (useEffFunc)
  {
    binboundaries[0] = 0.;
    
    double integral = 0;
    int currentbin = 1;
    for(unsigned int iev = 0; iev < size; iev++)
    {
      double val = values[iev];
      double eff = fEff->Eval(val);
      if (eff <= 1 && eff > 0.) integral += 1. / eff;
      else integral += 1.;
      
      if(integral > ((double)(currentbin))*(totalXsec/(double)(nbins)))
      {
        cout << "current bin = " << currentbin << " ; integral = " << integral << " ; sum = " << ((double)(currentbin))*(totalXsec/(double)(nbins)) << endl;
        binboundaries[currentbin] = val;
        currentbin++;
      }
    }
    cout << currentbin << endl;
  }
  else if (!useMC)// This way assumes all inefficiency is in the most peripheral bins
  {
    double EFF = gEff->GetVal();
    for(int i = 0; i < nbins; i++) {
      int entry = (int)( (float)(i)*((float)(size)/EFF/(float)(nbins)) - (float)(size)*(1. - EFF)/EFF );
      if(entry < 0 || i == 0) binboundaries[i] = 0;
      else binboundaries[i] = values[entry];
      if(binboundaries[i] < 0) { binboundaries[i] = 0; cout << "*"; }
    }
  }
  else {
      binboundaries[0] = 0.;

      double integral = 0;
      int currentbin = 1;
      for(unsigned int iev = 0; iev < size; iev++)
	{
	  double val = values[iev];
	  if (val<=threshold) integral += 1.*mcYscale;
	  else integral += 1.;

	  if(integral > ((double)(currentbin))*(totalXsec/(double)(nbins)))
	    {
	      cout << "current bin = " << currentbin << " ; integral = " << integral << " ; sum = " << ((double)(currentbin))*(totalXsec/(double)(nbins)) << endl;
	      binboundaries[currentbin] = val;
	      currentbin++;
	    }
	}
      cout << currentbin << endl;
  }
  for(int i = 0; i < nbins; i++) {
    if(binboundaries[i] < 0) binboundaries[i] = 0;
    txtfile << binboundaries[i] << ", ";
  }
  txtfile << binboundaries[nbins] << ")" << endl;
  txtfile << endl;
  
  txtfile<<"-------------------------------------"<<endl;
  txtfile<<"# Bin NpartMean NpartSigma NcollMean NcollSigma bMean bSigma BinEdge"<<endl;
  for(int i = 0; i < nbins; i++){
    int ii = nbins-i;
    
    //if (inputMCtable)
    //{
    //bins->table_[i].n_part_mean = inputMCtable->NpartMeanOfBin(i);
    //bins->table_[i].n_part_var = inputMCtable->NpartSigmaOfBin(i);
    //bins->table_[i].n_coll_mean = inputMCtable->NcollMeanOfBin(i);
    //bins->table_[i].n_coll_var = inputMCtable->NcollSigmaOfBin(i);
    //bins->table_[i].b_mean = inputMCtable->bMeanOfBin(i);
    //bins->table_[i].b_var = inputMCtable->bSigmaOfBin(i);
    //bins->table_[i].n_hard_mean = inputMCtable->NhardMeanOfBin(i);
    //bins->table_[i].n_hard_var = inputMCtable->NhardSigmaOfBin(i);
    //bins->table_[i].ecc2_mean  = inputMCtable->eccentricityMeanOfBin(i);
    //bins->table_[i].ecc2_var = inputMCtable->eccentricitySigmaOfBin(i);
    //}
    //else
    //{
      bins->table_[i].n_part_mean = -99;
      bins->table_[i].n_part_var = -99;
      bins->table_[i].n_coll_mean = -99;
      bins->table_[i].n_coll_var = -99;
      bins->table_[i].b_mean = -99;
      bins->table_[i].b_var = -99;
      bins->table_[i].n_hard_mean = -99;
      bins->table_[i].n_hard_var = -99;
      bins->table_[i].ecc2_mean  = -99;
      bins->table_[i].ecc2_var = -99;
      //}
    
      bins->table_[i].bin_edge = binboundaries[ii-1];
    
    txtfile << i << " " << bins->table_[i].n_part_mean << " " << bins->table_[i].n_part_var << " " << bins->table_[i].n_coll_mean << " " << bins->table_[i].n_coll_var << " " <<bins->table_[i].b_mean << " " << bins->table_[i].b_var << " " << bins->table_[i].n_hard_mean << " " << bins->table_[i].n_hard_var << " " << bins->table_[i].bin_edge << " " << endl;
  }
  txtfile << endl;
  txtfile<<"-------------------------------------"<<endl;
  
  outFile->cd();
  dir->cd();
  bins->Write();
  nt->Write();
  hfData->Write();
  hfMc->Write();
  //  bins->Delete();
  outFile->Write();
  outFile->Close();
  txtfile.close();
}
