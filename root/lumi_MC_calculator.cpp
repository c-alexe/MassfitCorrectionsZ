#include <ROOT/RDataFrame.hxx>
#include "TFile.h"

using namespace std;
using namespace ROOT;
using ROOT::RDF::RNode;

void lumi_MC_calculator() {
    ROOT::EnableImplicitMT();

    // Zmumu cross section in fb from Run 2
    // TODO confirm we can use same xsec as run2  xsec_ZmmPostVFP = 2001.9
    double xsec = 2001.9e+03;  

    // input files relevant
    vector<string> in_files = {};
  
	    in_files = {
			"/scratch/wmass/y2018v6/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv6/240718_131855/0000/Nano*.root",
			"/scratch/wmass/y2018v6/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv6/240718_131855/0001/Nano*.root",
			"/scratch/wmass/y2018v6/DYJetsToMuMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv6/240903_095439/0000/Nano*.root",
			"/scratch/wmass/y2018v6/DYJetsToMuMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv6/240903_095439/0001/Nano*.root",
			"/scratch/wmass/y2018v6/DYJetsToMuMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv6/240903_095439/0002/Nano*.root"
	    };
   
    //TFile* fout = TFile::Open("./inoutfiles/results/lumi_MC.root", "RECREATE");
    ROOT::RDataFrame d( "Events", in_files );
    
    auto dlast = std::make_unique<RNode>(d);
    std::cout <<"Total initial entries count is " << *(dlast->Count()) << std::endl;

    // Define MC weight
    dlast = std::make_unique<RNode>(dlast->Define("weight", [](float weight) -> float
	{
	  return std::copysign(1.0, weight);
    }, {"Generator_weight"} ));   

    // Book weight histo
    auto myHist = dlast->Histo1D({"h_genweight", "h_genweight", 2, -1., 1.1}, "weight");
    std::cout << "Equivalent integrated lumi of MC in fb^-1 = number of weighted events / xsec " << ( myHist->GetBinContent(2) - myHist->GetBinContent(1) ) / xsec << std::endl;
    // myHist->Write();

}