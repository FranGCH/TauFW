/*
 * @short: provide pileup weights at drawing level
 * @author: Izaak Neutelings (October 2018)
 * @edit: Francisco Casalinho (January 2026)
 *
 */

#include "TROOT.h"
#include "TFile.h"
#include "TH2.h"
#include "TH2F.h"
#include <iostream>
#include <algorithm>
#include "TSystem.h" // for gSystem
using namespace std;
TString _datadir_pu = gSystem->ExpandPathName("$CMSSW_BASE/src/TauFW/PicoProducer/data/pileup");
TString _fname_pu_data = _datadir_pu+"/Data_PileUp_2022_postEE.root";
TString _fname_pu_mc = _datadir_pu+"/MC_PileUp_2024.root";


vector<float> _pu_data_weights;
vector<float> _pu_mc_weights;



vector<float> getPUWeightsFromHist(TString filename, TString histname="pileup"){
  std::cout << ">>> opening " << filename << std::endl;
  TFile *file = new TFile(filename);
  TH1F* hist = (TH1F*) file->Get(histname);
  hist->SetDirectory(0);
  hist->Scale(1./hist->Integral());
  
  // Convert histogram to vector
  vector<float> weights;
  for(int i = 1; i <= hist->GetNbinsX(); ++i){
    weights.push_back(hist->GetBinContent(i));
  }
  
  file->Close();
  return weights;
}


void readPUFile(TString filename_data=_fname_pu_data, TString filename_mc=_fname_pu_mc){
  _pu_data_weights = getPUWeightsFromHist(filename_data);
  _pu_mc_weights = getPUWeightsFromHist(filename_mc);
}


Float_t getPUWeight(Int_t npu){
  // Bounds checking - ensure npu is within vector range
  if(npu < 0 || npu >= (Int_t)_pu_data_weights.size() || npu >= (Int_t)_pu_mc_weights.size()){
    return 1.0;
  }
  
  float data = _pu_data_weights[npu];
  float mc   = _pu_mc_weights[npu];
  
  if(mc > 0.){
    return data/mc;
  }
  return 1.0;
}


void pileup(){
  std::cout << ">>> initializing pileup.C ... " << std::endl;
}
