// Fits the Z mass in 4D bins (muon eta+, pt+, eta-, pt-) for data and MC extracting the mass scale and resolution biases.
// It can be used iteratively with massfit.cpp and resolfit.cpp to correct muon momenta in MC to get mass distributions closer to the data.
// Authors: Cristina Alexe, Lorenzo Bianchini

#include <ROOT/RDataFrame.hxx>
#include "TFile.h"
#include "TRandom3.h"
#include "TVector.h"
#include "TVectorT.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include "RooRealVar.h"
#include "RooDerivative.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooCrystalBall.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "TAxis.h"
#include "RooPlot.h"
#include <TMatrixD.h>
#include <TMatrixDSymfwd.h>
#include <TStopwatch.h>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <boost/program_options.hpp>
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/FCNGradientBase.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

//#include <Eigen/Core>
//#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;
using namespace ROOT;
using namespace RooFit;
using namespace ROOT::Minuit2;

typedef ROOT::VecOps::RVec<double> RVecD;
typedef ROOT::VecOps::RVec<unsigned int> RVecUI;
typedef ROOT::VecOps::RVec<int> RVecI;
typedef ROOT::VecOps::RVec<float> RVecF;
typedef ROOT::VecOps::RVec<bool> RVecB;
using ROOT::RDF::RNode;

using namespace boost::program_options;

//fb^-1
constexpr double lumiData2016 = 16.811; 
constexpr double lumiData2017 = 41.556;
constexpr double lumiData2018 = 59.688;

// with cross sections in fb at sqrt(s)=13Tev
constexpr double lumiMC2016 = 3.9375e+08/2001.9e+03;
constexpr double lumiMC2017 = 4.9803e+07/2001.9e+03;
constexpr double lumiMC2018 = 7.99137e+07/2001.9e+03; 

int main(int argc, char* argv[]) {

  TStopwatch sw;
  sw.Start();

  //ROOT::EnableImplicitMT();

  variables_map vm;
  try {
    options_description desc{"Options"};
    desc.add_options()
	  ("help,h", "Help screen")
	  ("minNumEvents",       value<int>()->default_value(100), "min number of events for a histogram to be accepted")
	  ("minNumEventsPerBin", value<int>()->default_value(10), "min number of events for a bin of a histogram to be accepted")
	  ("tag",                value<std::string>()->default_value("closure"), "run type, type of data used")
	  ("run",                value<std::string>()->default_value("closure"), "number of iteration")
	  ("saveMassFitHistos",  bool_switch()->default_value(false), "save pre and postfit mass distribution in 4D bin")
	  ("skipUnsmearedReco",  bool_switch()->default_value(true), "do not fill histograms of mass or jacobian weights per event for reco without smearing")
	  ("firstIter",          value<int>()->default_value(-1), "firstIter (-1: fills data histos, 0: fills MC histos, 1(needs 0): fills jacobians, 2(needs -1,0,1): fits for scale/resolution bias)")
	  ("lastIter",           value<int>()->default_value(-1), "lastIter  (-1: fills data histos, 0: fills MC histos, 1(needs 0): fills jacobians, 2(needs -1,0,1): fits for scale/resolution bias)")
	  ("nRMSforGausFit",     value<float>()->default_value(-1.), "number of RMS for Gaus mass difference fit")
	  ("minNumMassBins",     value<int>()->default_value(4), "min number of mass bins for a histogram to be accepted")
	  ("maxRMS",             value<float>()->default_value(-1.), "max RMS of Gaus mass difference fit for a 4D bin to be included in the mass fit")
	  ("rebin",              value<int>()->default_value(1), "rebin 4D bin histograms before scale fit")
	  ("fitWidth",           bool_switch()->default_value(false), "compute resolution bias")
	  ("fitNorm",            bool_switch()->default_value(false), "compute difference in normalisation in 4D bin")
	  ("usePrevMassFit",     bool_switch()->default_value(false), "use previous mass fit")
	  ("tagPrevMassFit",     value<std::string>()->default_value("closure"), "run type, type of data used")
	  ("runPrevMassFit",     value<std::string>()->default_value("closure"), "number of iteration")
	  ("usePrevResolFit",    bool_switch()->default_value(false), "use previous resolution fit")
	  ("tagPrevResolFit",    value<std::string>()->default_value("closure"), "run type, type of data used")
	  ("runPrevResolFit",    value<std::string>()->default_value("closure"), "number of iteration")
	  ("useKf",              bool_switch()->default_value(false), "use track input from Kalman Filter instead of CVH")
	  ("useCB",              bool_switch()->default_value(false), "under development")
	  ("scaleToData",        bool_switch()->default_value(false), "scale MC to data in 4D bin")
	  ("y2016",              bool_switch()->default_value(false), "use 2016 data")
	  ("y2017",              bool_switch()->default_value(false), "")
	  ("y2018",              bool_switch()->default_value(false), "");

    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);
    if (vm.count("help")) {
	  std::cout << desc << '\n';
	  return 0;
    }
    if (vm.count("tag")) std::cout << "Tag: " << vm["tag"].as<std::string>() << '\n';
    if (vm.count("run")) std::cout << "Run: " << vm["run"].as<std::string>() << '\n';
  }
  catch (const error &ex) {
    std::cerr << ex.what() << '\n';
  }
  
  int minNumEvents            = vm["minNumEvents"].as<int>();
  float nRMSforGausFit        = vm["nRMSforGausFit"].as<float>();
  std::string tag             = vm["tag"].as<std::string>();
  std::string run             = vm["run"].as<std::string>();
  int minNumEventsPerBin      = vm["minNumEventsPerBin"].as<int>();
  int minNumMassBins          = vm["minNumMassBins"].as<int>();
  int rebin                   = vm["rebin"].as<int>();
  int firstIter               = vm["firstIter"].as<int>();
  int lastIter                = vm["lastIter"].as<int>();
  bool saveMassFitHistos      = vm["saveMassFitHistos"].as<bool>();
  bool skipUnsmearedReco      = vm["skipUnsmearedReco"].as<bool>();
  bool fitWidth               = vm["fitWidth"].as<bool>();
  bool fitNorm                = vm["fitNorm"].as<bool>();
  bool usePrevMassFit         = vm["usePrevMassFit"].as<bool>();
  bool usePrevResolFit        = vm["usePrevResolFit"].as<bool>();
  bool useKf                  = vm["useKf"].as<bool>();
  bool useCB                  = vm["useCB"].as<bool>();
  bool y2016                  = vm["y2016"].as<bool>();
  bool y2017                  = vm["y2017"].as<bool>();
  bool y2018                  = vm["y2018"].as<bool>();
  std::string tagPrevMassFit  = vm["tagPrevMassFit"].as<std::string>();
  std::string runPrevMassFit  = vm["runPrevMassFit"].as<std::string>();
  std::string tagPrevResolFit = vm["tagPrevResolFit"].as<std::string>();
  std::string runPrevResolFit = vm["runPrevResolFit"].as<std::string>();
  bool scaleToData            = vm["scaleToData"].as<bool>();
  float maxRMS                = vm["maxRMS"].as<float>();
  
  assert( firstIter>=-1 && lastIter<=2 && firstIter<=lastIter );
  assert( y2016 || y2017 || y2018 );

  double lumi = lumiData2016;
  if(y2017)      lumi = lumiData2017;
  else if(y2018) lumi = lumiData2018; 

  vector<float> pt_edges  = {25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0}; 
  vector<float> eta_edges = {-2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0,
                             0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4};

  TH1F* h_pt_edges  = new TH1F("h_pt_edges", "",  pt_edges.size()-1, pt_edges.data());
  TH1F* h_eta_edges = new TH1F("h_eta_edges", "", eta_edges.size()-1, eta_edges.data());
  
  unsigned int n_pt_bins  = pt_edges.size()-1;
  unsigned int n_eta_bins = eta_edges.size()-1;
  // Number of 4D bins in muon kinematics (eta+, pt+, eta-, pt-)
  int n_bins = n_pt_bins*n_pt_bins*n_eta_bins*n_eta_bins; 

  // Bins in mass
  const int x_nbins   = 40;
  const double x_low  = 70.0;
  const double x_high = 110.0;

  // Bins in mass - gen mass
  const int dm_bins    = 24;
  const double dm_low  = -6.0;
  const double dm_high = 6.0;

  // _nom histograms with AeM set to 0 (remnant from masscales.cpp, needed for massfit.cpp to run for both data and toys) 
  TH1F* h_A_vals_nom = new TH1F("h_A_vals_nom", "", n_eta_bins, 0, n_eta_bins );
  TH1F* h_e_vals_nom = new TH1F("h_e_vals_nom", "", n_eta_bins, 0, n_eta_bins );
  TH1F* h_M_vals_nom = new TH1F("h_M_vals_nom", "", n_eta_bins, 0, n_eta_bins );
  // _prevfit histograms for the sum of the pT scale bias parameters A, e or M from all the previous iterations 
  TH1F* h_A_vals_prevfit = new TH1F("h_A_vals_prevfit", "", n_eta_bins, 0, n_eta_bins );
  TH1F* h_e_vals_prevfit = new TH1F("h_e_vals_prevfit", "", n_eta_bins, 0, n_eta_bins );
  TH1F* h_M_vals_prevfit = new TH1F("h_M_vals_prevfit", "", n_eta_bins, 0, n_eta_bins );
  // _prevfit histograms for the sum of the curvature resolution bias parameters c or d from all the previous iterations
  TH1F* h_c_vals_prevfit = new TH1F("h_c_vals_prevfit", "", n_eta_bins, 0, n_eta_bins );
  TH1F* h_d_vals_prevfit = new TH1F("h_d_vals_prevfit", "", n_eta_bins, 0, n_eta_bins );
  
  // Store curvature (not pT) biases
  VectorXd A_vals_fit( n_eta_bins );
  VectorXd e_vals_fit( n_eta_bins );
  VectorXd M_vals_fit( n_eta_bins );
  VectorXd c_vals_fit( n_eta_bins );
  VectorXd d_vals_fit( n_eta_bins );

  // Initialize curvature bias parameters A,e,M,c,d = 0, they will remain 0 if usePrevMassFit/usePrevResolFit are false
  for(unsigned int i=0; i<n_eta_bins; i++) {
    h_A_vals_nom->SetBinContent(i+1, 0.0);
    h_e_vals_nom->SetBinContent(i+1, 0.0);
    h_M_vals_nom->SetBinContent(i+1, 0.0);
    h_A_vals_prevfit->SetBinContent(i+1, 0.0);
    h_e_vals_prevfit->SetBinContent(i+1, 0.0);
    h_M_vals_prevfit->SetBinContent(i+1, 0.0);
    h_c_vals_prevfit->SetBinContent(i+1, 0.0);
    h_d_vals_prevfit->SetBinContent(i+1, 0.0);
    A_vals_fit(i) = 0.0;
    e_vals_fit(i) = 0.0;
    M_vals_fit(i) = 0.0;
	c_vals_fit(i) = 0.0;
	d_vals_fit(i) = 0.0;
  }           
  
  // Option to keep track of MC and jacobians before(reco) and after(smear0) applying curvature corrections to get closer to data
  // If skipUnsmearedReco == false, we will save mass and jacobian histograms for both reco and smear0 in the dataframe
  // In the scale fit, only smear0 jacobians are used
  std::vector<string> recos = {"reco", "smear0"};

  // Map to histograms containing information for each 4D bin
  std::map<string, TH1D*> h_map;
  for(unsigned int r = 0; r<recos.size(); r++) {
	if(skipUnsmearedReco && recos[r]=="reco") continue;
	// Gaussian mean and rms of the mass - gen mass distribution in a 4D bin
    h_map.insert( std::make_pair<string, TH1D* >("mean_"+recos[r], 0 ) );
    h_map.insert( std::make_pair<string, TH1D* >("rms_"+recos[r],  0 ) );
	// 1/0 if keeping(ignoring) a 4D bin in the fit
    h_map.insert( std::make_pair<string, TH1D* >("mask_"+recos[r],  0 ) );
  }

  // Map to 2D histograms for Crystal Ball jacobians
  std::map<string, TH2D*> h_jac_map;
  for(unsigned int r = 0; r<recos.size(); r++) {
	if(skipUnsmearedReco && recos[r]=="reco") continue;
	// Crystal Ball jacobian event weights in a 4D bin, in a mass - gen mass bin
    h_jac_map.insert( std::make_pair<string, TH2D* >("jscale_cb_per_evt_"+recos[r], 0 ) );
    h_jac_map.insert( std::make_pair<string, TH2D* >("jwidth_cb_per_evt_"+recos[r], 0 ) );
  }

  // Map to positions in RVecF "masses", position 0 is gen (also used slightly differently for weights_jac)
  std::map<string, unsigned int> idx_map;
  idx_map.insert( std::make_pair<string, unsigned int >("reco",   1 ) );
  idx_map.insert( std::make_pair<string, unsigned int >("smear0", 2 ) );

  if(usePrevMassFit) {
    TFile* ffit = TFile::Open(("./massfit_"+tagPrevMassFit+"_"+runPrevMassFit+".root").c_str(), "READ");
    if(ffit!=0) {    
      cout << "Using fit results from " <<  std::string(ffit->GetName()) << " as new nominal for smear0" << endl;
      // Read the sum of the pT scale bias parameters A, e or M from all the previous iterations
      TH1D* h_A_vals_prevfit_in = (TH1D*)ffit->Get("h_A_vals_prevfit");
      TH1D* h_e_vals_prevfit_in = (TH1D*)ffit->Get("h_e_vals_prevfit");
      TH1D* h_M_vals_prevfit_in = (TH1D*)ffit->Get("h_M_vals_prevfit");
      for(unsigned int i=0; i<n_eta_bins; i++) {
        // We will use curvature scale biases, which are (-1) * (pT scale biases)
	    A_vals_fit(i) = -h_A_vals_prevfit_in->GetBinContent(i+1);
	    e_vals_fit(i) = -h_e_vals_prevfit_in->GetBinContent(i+1);
	    M_vals_fit(i) = -h_M_vals_prevfit_in->GetBinContent(i+1);
      }
      // Save the content of h_ _vals_prevfit_in to be passed to massfit.cpp without further changes
      h_A_vals_prevfit->Add(h_A_vals_prevfit_in, 1.0);
      h_e_vals_prevfit->Add(h_e_vals_prevfit_in, 1.0);
      h_M_vals_prevfit->Add(h_M_vals_prevfit_in, 1.0);
      ffit->Close();
    }
    else {
      cout << "No mass fit file!" << endl;
    }
  }

  if(usePrevResolFit) {
    TFile* ffit = TFile::Open(("./resolfit_"+tagPrevResolFit+"_"+runPrevResolFit+".root").c_str(), "READ");
    if(ffit!=0) {    
      cout << "Using fit results from " <<  std::string(ffit->GetName()) << " as MC smear" << endl;
	  // Read the sum of the resolution biases c or d from all the previous iterations
      TH1D* h_c_vals_prevfit_in = (TH1D*)ffit->Get("h_c_vals_prevfit");
      TH1D* h_d_vals_prevfit_in = (TH1D*)ffit->Get("h_d_vals_prevfit");
      for(unsigned int i=0; i<n_eta_bins; i++){
	    c_vals_fit(i) = h_c_vals_prevfit_in->GetBinContent(i+1);
	    d_vals_fit(i) = h_d_vals_prevfit_in->GetBinContent(i+1);
      }
	  // Save the content of h_ _vals_prevfit_in to be passed to resolfit.cpp without further changes
      h_c_vals_prevfit->Add(h_c_vals_prevfit_in, +1.0);
      h_d_vals_prevfit->Add(h_d_vals_prevfit_in, +1.0);
      ffit->Close();
    }
    else {
      cout << "No smear fit file!" << endl;
    }    
  }  

  // Define a single output file, we will write to and read from it at the different iterations 
  // If firstIter = 2, update an existing output file with iter -1,0 and 1 to (over)write iter 2 (the mass fit results)
  TFile* fout = TFile::Open(("./massscales_"+tag+"_"+run+".root").c_str(), firstIter<2 ? "RECREATE" : "UPDATE");
  
  // iter -1 -> data mass histos
  // iter  0 -> MC mass histos + calculation of jacobian terms per event
  // iter  1 -> MC jacobian histos
  // iter  2 -> fit for beta(= mass scale bias) [,alpha(= mass resolution bias), nu(= mass normalization bias)]
  
  // By the end the output file will contain: 
  // - histograms resulted from the event-loop at different iterations
  // - TTree and histograms resulted from the mass fits
  // - OPTIONAL, in the postfit/ folder, pre and postfit mass distribution for the 4D bins

  for(int iter=-1; iter<3; iter++) {

    if( !(iter>=firstIter && iter<=lastIter) ) continue;
    cout << "Doing iter " << iter << endl;

    // Read the input files relevant to the current iteration
    vector<string> in_files = {};
    if(iter>=0) { // MC
      if(y2016) {
        in_files = {
	      "/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_040854/0000/NanoV9MCPostVFP_*.root",
	      "/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_041233/0000/NanoV9MCPostVFP_*.root",
	      "/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_041233/0001/NanoV9MCPostVFP_*.root",
	      "/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_041233/0002/NanoV9MCPostVFP_*.root"
	    };
      }
      else if(y2017) {
	    in_files = {
	      "/scratch/wmass/y2017/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2017_TrackFitV722_NanoProdv3/NanoV9MC2017_*.root"
	    };
      }
      else if(y2018) {
	    in_files = {
	      "/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_1.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_10.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_100.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_101.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_102.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_103.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_104.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_105.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_107.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_108.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_109.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_11.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_110.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_111.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_112.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_113.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_114.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_115.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_116.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_117.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_118.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_119.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_12.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_120.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_121.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_122.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_123.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_125.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_126.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_128.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_129.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_13.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_130.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_131.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_132.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_133.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_135.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_138.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_14.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_142.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_145.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_146.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_147.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_148.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_149.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_15.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_151.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_153.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_154.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_155.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_156.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_157.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_158.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_159.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_16.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_160.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_161.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_162.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_163.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_164.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_165.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_166.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_167.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_168.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_169.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_17.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_170.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_171.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_172.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_173.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_174.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_175.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_176.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_177.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_178.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_179.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_18.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_180.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_182.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_183.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_184.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_185.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_186.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_187.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_188.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_189.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_19.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_190.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_194.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_196.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_197.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_199.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_2.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_20.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_200.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_202.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_204.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_205.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_206.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_207.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_208.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_209.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_21.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_210.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_211.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_212.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_213.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_215.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_217.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_218.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_22.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_220.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_221.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_222.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_223.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_224.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_225.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_226.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_227.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_228.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_229.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_23.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_230.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_231.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_232.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_233.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_234.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_235.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_236.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_237.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_239.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_24.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_240.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_241.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_242.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_243.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_244.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_245.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_246.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_247.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_248.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_249.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_25.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_250.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_251.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_253.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_254.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_255.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_256.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_257.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_258.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_259.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_26.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_260.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_261.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_263.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_264.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_265.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_266.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_267.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_268.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_269.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_27.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_270.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_271.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_272.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_273.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_274.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_275.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_276.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_277.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_278.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_279.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_28.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_280.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_281.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_282.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_283.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_284.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_285.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_286.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_288.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_289.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_29.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_290.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_291.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_293.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_294.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_295.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_296.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_297.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_298.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_299.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_3.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_30.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_300.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_301.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_302.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_303.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_304.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_305.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_306.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_307.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_308.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_309.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_31.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_310.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_311.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_312.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_313.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_314.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_315.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_316.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_317.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_318.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_319.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_32.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_320.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_321.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_322.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_323.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_324.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_325.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_326.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_327.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_329.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_33.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_330.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_332.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_334.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_335.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_336.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_337.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_338.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_339.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_34.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_340.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_341.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_342.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_343.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_344.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_345.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_346.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_347.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_348.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_349.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_35.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_350.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_351.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_352.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_353.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_354.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_355.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_356.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_357.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_358.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_359.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_36.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_360.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_361.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_362.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_363.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_364.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_365.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_366.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_367.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_368.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_369.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_37.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_370.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_371.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_372.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_373.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_374.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_375.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_376.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_377.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_379.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_38.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_380.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_381.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_382.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_383.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_384.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_385.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_386.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_387.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_388.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_389.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_39.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_391.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_392.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_393.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_394.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_395.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_396.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_397.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_398.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_399.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_4.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_401.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_402.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_403.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_404.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_405.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_406.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_407.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_408.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_409.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_410.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_411.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_412.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_413.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_414.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_415.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_416.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_418.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_419.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_42.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_420.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_421.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_422.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_423.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_424.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_425.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_426.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_427.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_428.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_429.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_43.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_430.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_431.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_432.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_433.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_435.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_436.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_437.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_438.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_439.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_44.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_440.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_441.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_442.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_445.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_446.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_447.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_448.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_449.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_45.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_450.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_451.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_452.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_453.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_454.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_455.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_456.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_457.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_458.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_459.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_46.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_461.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_463.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_464.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_465.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_467.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_468.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_470.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_471.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_472.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_473.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_475.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_477.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_479.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_48.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_481.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_483.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_484.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_485.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_486.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_487.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_488.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_489.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_490.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_491.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_492.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_493.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_494.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_495.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_496.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_497.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_498.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_499.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_5.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_500.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_501.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_502.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_503.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_504.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_505.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_506.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_507.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_508.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_509.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_51.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_510.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_511.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_512.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_513.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_514.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_515.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_516.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_517.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_518.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_519.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_52.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_520.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_521.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_522.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_523.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_524.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_525.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_526.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_527.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_528.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_529.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_53.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_530.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_531.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_532.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_533.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_534.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_535.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_536.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_537.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_538.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_539.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_54.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_541.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_542.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_543.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_544.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_545.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_546.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_547.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_548.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_549.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_55.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_550.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_551.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_552.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_553.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_554.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_555.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_556.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_557.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_558.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_559.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_56.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_560.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_561.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_562.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_563.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_564.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_565.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_566.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_567.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_568.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_569.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_57.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_570.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_571.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_572.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_573.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_574.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_575.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_576.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_577.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_578.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_579.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_58.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_580.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_581.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_582.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_583.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_584.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_585.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_586.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_587.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_588.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_589.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_59.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_590.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_591.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_592.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_593.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_594.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_595.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_596.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_597.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_598.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_599.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_6.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_60.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_600.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_601.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_602.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_603.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_604.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_605.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_606.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_607.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_608.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_609.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_61.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_610.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_612.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_614.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_615.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_616.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_617.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_618.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_619.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_62.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_620.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_621.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_622.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_623.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_624.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_625.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_626.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_627.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_628.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_629.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_630.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_631.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_632.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_633.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_634.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_635.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_636.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_637.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_638.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_639.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_64.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_640.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_641.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_642.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_643.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_644.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_645.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_646.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_647.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_648.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_649.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_65.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_650.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_651.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_652.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_653.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_654.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_655.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_656.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_657.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_658.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_659.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_66.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_660.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_661.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_662.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_663.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_665.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_666.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_667.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_668.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_669.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_67.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_670.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_671.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_672.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_673.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_674.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_675.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_676.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_677.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_678.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_679.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_68.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_680.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_681.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_682.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_683.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_684.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_685.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_686.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_687.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_688.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_689.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_69.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_690.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_691.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_692.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_693.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_694.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_695.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_696.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_697.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_698.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_699.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_7.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_70.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_700.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_701.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_702.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_703.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_704.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_705.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_706.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_707.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_708.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_709.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_71.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_711.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_712.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_713.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_714.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_715.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_716.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_717.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_718.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_719.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_72.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_720.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_721.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_722.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_723.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_724.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_725.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_726.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_727.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_728.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_729.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_73.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_730.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_731.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_732.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_733.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_734.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_735.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_736.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_737.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_738.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_739.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_74.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_740.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_741.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_742.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_743.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_744.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_745.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_746.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_747.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_748.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_749.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_750.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_751.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_752.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_753.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_754.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_755.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_756.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_757.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_758.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_759.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_76.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_760.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_761.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_764.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_765.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_766.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_767.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_768.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_769.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_77.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_770.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_771.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_772.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_773.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_774.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_775.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_776.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_777.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_778.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_779.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_78.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_780.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_781.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_782.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_783.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_785.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_786.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_787.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_788.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_789.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_79.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_790.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_791.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_792.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_793.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_794.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_795.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_796.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_799.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_8.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_800.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_801.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_802.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_803.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_804.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_805.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_806.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_807.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_808.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_809.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_81.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_810.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_812.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_813.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_814.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_815.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_816.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_817.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_818.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_82.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_820.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_822.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_823.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_825.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_826.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_827.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_828.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_829.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_830.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_831.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_832.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_833.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_835.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_836.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_837.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_838.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_839.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_84.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_841.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_843.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_844.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_845.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_846.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_847.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_848.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_849.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_85.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_850.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_851.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_852.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_853.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_854.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_855.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_856.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_857.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_858.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_859.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_86.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_860.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_861.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_862.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_863.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_864.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_865.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_866.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_868.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_869.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_87.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_871.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_872.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_873.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_874.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_875.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_876.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_877.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_878.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_879.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_880.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_881.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_883.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_884.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_885.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_886.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_887.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_888.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_89.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_890.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_891.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_892.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_893.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_894.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_895.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_896.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_897.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_898.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_9.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_90.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_900.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_901.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_902.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_903.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_904.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_906.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_907.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_908.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_91.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_910.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_911.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_912.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_913.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_914.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_915.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_916.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_917.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_918.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_919.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_920.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_921.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_922.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_923.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_924.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_925.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_926.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_927.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_928.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_929.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_93.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_930.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_931.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_932.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_933.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_934.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_935.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_936.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_937.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_938.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_939.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_94.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_940.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_942.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_943.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_944.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_945.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_946.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_947.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_948.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_949.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_95.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_950.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_952.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_953.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_954.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_955.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_956.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_957.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_958.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_959.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_960.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_962.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_963.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_964.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_965.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_966.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_967.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_968.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_969.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_97.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_970.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_971.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_972.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_973.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_974.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_975.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_976.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_977.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_979.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_98.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_980.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_981.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_982.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_983.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_984.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_985.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_986.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_987.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_988.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_997.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_998.root",
"/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_999.root",
	      "/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0001/NanoV9MC2018_*.root"
	    };
      }      
    }
    else { // data
      if(y2016) {
	    in_files = {
	      "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016FDataPostVFP_TrackFitV722_NanoProdv6/240509_051502/0000/NanoV9DataPostVFP_*.root",
	      "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016GDataPostVFP_TrackFitV722_NanoProdv6/240509_051653/0000/NanoV9DataPostVFP_*.root",
	      "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016GDataPostVFP_TrackFitV722_NanoProdv6/240509_051653/0001/NanoV9DataPostVFP_*.root",
	      "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016GDataPostVFP_TrackFitV722_NanoProdv6/240509_051653/0002/NanoV9DataPostVFP_*.root",
	      "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016GDataPostVFP_TrackFitV722_NanoProdv6/240509_051653/0003/NanoV9DataPostVFP_*.root",
	      "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016GDataPostVFP_TrackFitV722_NanoProdv6/240509_051653/0004/NanoV9DataPostVFP_*.root",
	      "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016HDataPostVFP_TrackFitV722_NanoProdv6/240509_051807/0000/NanoV9DataPostVFP_*.root",
	      "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016HDataPostVFP_TrackFitV722_NanoProdv6/240509_051807/0001/NanoV9DataPostVFP_*.root",
	      "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016HDataPostVFP_TrackFitV722_NanoProdv6/240509_051807/0002/NanoV9DataPostVFP_*.root",
	      "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016HDataPostVFP_TrackFitV722_NanoProdv6/240509_051807/0003/NanoV9DataPostVFP_*.root",
	      "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016HDataPostVFP_TrackFitV722_NanoProdv6/240509_051807/0004/NanoV9DataPostVFP_*.root",
	      "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016HDataPostVFP_TrackFitV722_NanoProdv6/240509_051807/0005/NanoV9DataPostVFP_*.root",
	      "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016HDataPostVFP_TrackFitV722_NanoProdv6/240509_051807/0006/NanoV9DataPostVFP_*.root"
	    };
      }
      else if(y2017) {
	    in_files = {
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017B_TrackFitV722_NanoProdv3/240127_110915/0000/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017B_TrackFitV722_NanoProdv3/240127_110915/0001/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017B_TrackFitV722_NanoProdv3/240127_110915/0002/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017B_TrackFitV722_NanoProdv3/240127_110915/0003/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017C_TrackFitV722_NanoProdv3/240127_115941/0000/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017C_TrackFitV722_NanoProdv3/240127_115941/0001/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017C_TrackFitV722_NanoProdv3/240127_115941/0002/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017C_TrackFitV722_NanoProdv3/240127_115941/0003/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017C_TrackFitV722_NanoProdv3/240127_115941/0004/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017C_TrackFitV722_NanoProdv3/240127_115941/0005/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017D_TrackFitV722_NanoProdv3/240127_120137/0000/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017D_TrackFitV722_NanoProdv3/240127_120137/0001/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017D_TrackFitV722_NanoProdv3/240127_120137/0002/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017E_TrackFitV722_NanoProdv3/240127_121346/0000/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017E_TrackFitV722_NanoProdv3/240127_121346/0001/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017E_TrackFitV722_NanoProdv3/240127_121346/0002/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017E_TrackFitV722_NanoProdv3/240127_121346/0003/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017E_TrackFitV722_NanoProdv3/240127_121346/0004/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017F_TrackFitV722_NanoProdv3/240127_122701/0000/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017F_TrackFitV722_NanoProdv3/240127_122701/0001/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017F_TrackFitV722_NanoProdv3/240127_122701/0002/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017F_TrackFitV722_NanoProdv3/240127_122701/0003/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017F_TrackFitV722_NanoProdv3/240127_122701/0004/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017F_TrackFitV722_NanoProdv3/240127_122701/0005/NanoV9Data2017_*.root",
	      "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017F_TrackFitV722_NanoProdv3/240127_122701/0006/NanoV9Data2017_*.root"
	    };
      }
      else if(y2018) {
	    in_files = {
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018A_TrackFitV722_NanoProdv3/231102_185937/0000/NanoV9Data2018_*.root",
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018A_TrackFitV722_NanoProdv3/231102_185937/0001/NanoV9Data2018_*.root",
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018A_TrackFitV722_NanoProdv3/231102_185937/0002/NanoV9Data2018_*.root",
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018A_TrackFitV722_NanoProdv3/231102_185937/0003/NanoV9Data2018_*.root",
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018A_TrackFitV722_NanoProdv3/231102_185937/0004/NanoV9Data2018_*.root",
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018A_TrackFitV722_NanoProdv3/231102_185937/0005/NanoV9Data2018_*.root",
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018A_TrackFitV722_NanoProdv3/231102_185937/0006/NanoV9Data2018_*.root",
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018B_TrackFitV722_NanoProdv3/231103_093816/0000/NanoV9Data2018_*.root",
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018B_TrackFitV722_NanoProdv3/231103_093816/0001/NanoV9Data2018_*.root",
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018B_TrackFitV722_NanoProdv3/231103_093816/0002/NanoV9Data2018_*.root",
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018C_TrackFitV722_NanoProdv3/231103_101410/0000/NanoV9Data2018_*.root",
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018C_TrackFitV722_NanoProdv3/231103_101410/0001/NanoV9Data2018_*.root",
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018C_TrackFitV722_NanoProdv3/231103_101410/0002/NanoV9Data2018_*.root",
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018D_TrackFitV722_NanoProdv3/231107_134901/0000/NanoV9Data2018_*.root",
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018D_TrackFitV722_NanoProdv3/231107_134901/0001/NanoV9Data2018_*.root",
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018D_TrackFitV722_NanoProdv3/231107_134901/0002/NanoV9Data2018_*.root",
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018D_TrackFitV722_NanoProdv3/231107_134901/0003/NanoV9Data2018_*.root",
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018D_TrackFitV722_NanoProdv3/231107_134901/0004/NanoV9Data2018_*.root",
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018D_TrackFitV722_NanoProdv3/231107_134901/0005/NanoV9Data2018_*.root",
	      "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018D_TrackFitV722_NanoProdv3/231107_134901/0006/NanoV9Data2018_*.root"
	    }; 
      }
    }
    
	// Define dataframe for the input files relevant to the current iteration 
    ROOT::RDataFrame d( "Events", in_files );
    auto dlast = std::make_unique<RNode>(d);
	double number_input_events = *(dlast->Count());  
    std::cout << " Input file read. Initial event count is " << number_input_events << std::endl;
        
    if(iter>=0) { // MC

      // Define the indices of individual muons passing selection criteria
      dlast = std::make_unique<RNode>(dlast->Define("idxs", [&](UInt_t nMuon, RVecB Muon_looseId, RVecF Muon_dxybs, RVecB Muon_isGlobal, 
								       RVecB Muon_highPurity, RVecB Muon_mediumId, RVecF Muon_pfRelIso04_all, RVecF Muon_pt, RVecF Muon_eta) -> RVecUI 
      {
	    RVecUI out;
	    for(unsigned int i = 0; i < nMuon; i++){
	      if( Muon_looseId[i] && TMath::Abs(Muon_dxybs[i]) < 0.05 && Muon_isGlobal[i] && Muon_highPurity[i] && Muon_mediumId[i] && Muon_pfRelIso04_all[i]<0.15 &&
	        Muon_pt[i] >= pt_edges[0] && Muon_pt[i] < pt_edges[ n_pt_bins ]  && Muon_eta[i]>=eta_edges[0] && Muon_eta[i]<=eta_edges[ n_eta_bins ] ) out.emplace_back(i);
		}
	    return out;
      }, {"nMuon", "Muon_looseId", "Muon_dxybs", "Muon_isGlobal", "Muon_highPurity","Muon_mediumId", "Muon_pfRelIso04_all",
		 useKf ? "Muon_pt" : "Muon_cvhidealPt", useKf ? "Muon_eta" : "Muon_cvhidealEta"} ));

      // Filter to keep only events with exactly 2 oppositely charged, selected muons
      dlast = std::make_unique<RNode>(dlast->Filter( [](RVecUI idxs, RVecI Muon_charge, bool HLT_IsoMu24 )
	  {
	    if( idxs.size()!=2 || !HLT_IsoMu24) return false;
	    if( Muon_charge[idxs[0]]*Muon_charge[idxs[1]] > 0 ) return false;
	    return true;
      }, {"idxs", "Muon_charge", "HLT_IsoMu24"} ));
      
      // Define MC weight
      dlast = std::make_unique<RNode>(dlast->Define("weight", [](float weight) -> float
	  {
	    return std::copysign(1.0, weight);
      }, {"Generator_weight"} ));          
      
      // Define pos and neg curvature k smeared according to the curvature biases A,e,M,c,d computed in previous iterations
      dlast = std::make_unique<RNode>(dlast->Define("Muon_ksmear", [&](RVecUI idxs,
									   RVecF Muon_pt, RVecF Muon_eta, RVecF Muon_phi, RVecF Muon_mass, RVecI Muon_charge,
									   UInt_t nGenPart, RVecI GenPart_status, RVecI GenPart_statusFlags, RVecI GenPart_pdgId,
									   RVecF GenPart_pt, RVecF GenPart_eta, RVecF GenPart_phi, RVecF GenPart_mass) -> RVecF
	  {
	    RVecF out;
	    unsigned int idxP = Muon_charge[idxs[0]]>0 ? idxs[0] : idxs[1];
	    unsigned int idxM = Muon_charge[idxs[0]]>0 ? idxs[1] : idxs[0];
	    // reco muon
	    ROOT::Math::PtEtaPhiMVector muP( Muon_pt[ idxP ], Muon_eta[ idxP ], Muon_phi[ idxP ], Muon_mass[ idxP ] );
	    ROOT::Math::PtEtaPhiMVector muM( Muon_pt[ idxM ], Muon_eta[ idxM ], Muon_phi[ idxM ], Muon_mass[ idxM ] );
	    // gen matching
	    ROOT::Math::PtEtaPhiMVector gmuP( 0., 0., 0., 0. );
	    ROOT::Math::PtEtaPhiMVector gmuM( 0., 0., 0., 0. );
	    for(unsigned int i = 0; i < nGenPart; i++) {
	      bool isGoodGenPart = (GenPart_status[i]==1 && (GenPart_statusFlags[i] & 1 || (GenPart_statusFlags[i] & (1<<5))) && TMath::Abs(GenPart_pdgId[i])==13);
	      if(!isGoodGenPart) continue;
	      ROOT::Math::PtEtaPhiMVector gen(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
	      if( ROOT::Math::VectorUtil::DeltaR(gen, muP) < 0.1 && ROOT::Math::VectorUtil::DeltaR(gen, muM) > 0.1 ) gmuP = gen;
	      else if( ROOT::Math::VectorUtil::DeltaR(gen, muP) > 0.1 && ROOT::Math::VectorUtil::DeltaR(gen, muM) < 0.1 ) gmuM = gen;
	    }
	
	    // gen pT cut
	    if( gmuP.Pt()>10. && gmuM.Pt()>10. ) {
	  
	      float kmuP = 1./muP.Pt();
	      float kmuM = 1./muM.Pt();
	      float kgmuP = 1./gmuP.Pt();
	      float kgmuM = 1./gmuM.Pt();
	  
	      float scale_smear0P = 1.0;
	      float scale_smear0M = 1.0;
	      float resol_smear0P = 0.0;
	      float resol_smear0M = 0.0;
	  
	      // AeMcd are eta dependent 
	      unsigned int ietaP = n_eta_bins;
	      for(unsigned int ieta_p = 0; ieta_p<n_eta_bins; ieta_p++) {
	        float eta_p_low = eta_edges[ieta_p];
	        float eta_p_up  = eta_edges[ieta_p+1];
	        if( muP.Eta()>=eta_p_low && muP.Eta()<eta_p_up) ietaP = ieta_p;	  
	      }
	      unsigned int ietaM = n_eta_bins;
	      for(unsigned int ieta_m = 0; ieta_m<n_eta_bins; ieta_m++) {
	        float eta_m_low = eta_edges[ieta_m];
	        float eta_m_up  = eta_edges[ieta_m+1];
	        if( muM.Eta()>=eta_m_low && muM.Eta()<eta_m_up) ietaM = ieta_m;	  
	      }

	      if(ietaP<n_eta_bins && ietaM<n_eta_bins) {
            // Correct the MC curvature with the curvature scale biases derived in previous iterations (which are respectively equal to (-1)* sum of the pT scale biases from previous iterations)
	        // if usePrevMassFit is false, A,e,M are 0
			scale_smear0P = (1. + A_vals_fit(ietaP) - e_vals_fit(ietaP)*kmuP + M_vals_fit(ietaP)/kmuP);
	        scale_smear0M = (1. + A_vals_fit(ietaM) - e_vals_fit(ietaM)*kmuM - M_vals_fit(ietaM)/kmuM);
	        //cout << "smear0:" << scale_smear0P << ": " << 1 << " + " << A_vals_fit(ietaP) << " - " << e_vals_fit(ietaP)*kmuP << " + " << M_vals_fit(ietaP)/kmuP << endl;
	        cout << "scale_smear0P: " << scale_smear0P << "for ietap"<< ietaP <<" scale_smear0M: "<<scale_smear0M << "for ietam"<< ietaM << endl;
			if(usePrevResolFit) {
			  // Correct the MC curvature with the resolution biases derived in previous iterations (which are respectively equal to the sum of the resolution biases from previous iterations)  	
	          resol_smear0P = TMath::Sqrt( TMath::Max( 1.0 + c_vals_fit(ietaP) + d_vals_fit(ietaP)*kmuP, 0.0)  ) - 1.0;
	          resol_smear0M = TMath::Sqrt( TMath::Max( 1.0 + c_vals_fit(ietaM) + d_vals_fit(ietaM)*kmuM, 0.0)  ) - 1.0;
	          //cout << "resol_smear0P: sqrt( max(1.0 + " << c_vals_fit(ietaP)  << " + " << d_vals_fit(ietaP)*kmuP << ")) - 1.0 = " << resol_smear0P << endl;
	          //cout << "resol_smear0M: sqrt( max(1.0 + " << c_vals_fit(ietaM)  << " + " << d_vals_fit(ietaM)*kmuM << ")) - 1.0 = " << resol_smear0M << endl;  
	        }
	        float kmuPsmear0 = (kgmuP + (kmuP - kgmuP)*(1.0 + resol_smear0P))*scale_smear0P; // if A,e,M,c,d = 0, kmuPsmear0 = kmuP
	        float kmuMsmear0 = (kgmuM + (kmuM - kgmuM)*(1.0 + resol_smear0M))*scale_smear0M; // if A,e,M,c,d = 0, kmuMsmear0 = kmuM	  
	        out.emplace_back( kmuPsmear0 );
	        out.emplace_back( kmuMsmear0 );
		  }	  
		  else {
	        out.emplace_back(0.0);
	        out.emplace_back(0.0);
	      }	      
	    }
	    else {
	      out.emplace_back(0.0);
	      out.emplace_back(0.0);
	    }
	    return out;
	  }, {"idxs", useKf ? "Muon_pt" : "Muon_cvhidealPt", useKf ? "Muon_eta" : "Muon_cvhidealEta", useKf ? "Muon_phi" : "Muon_cvhidealPhi", "Muon_mass",
		 "Muon_charge", "nGenPart", "GenPart_status", "GenPart_statusFlags", "GenPart_pdgId", "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass"} ));
      
	  // Define eta+, pt+, eta-, pt- indexes for each muon pair that passed selection. The 1st entry in "indexes" is for reco, the 2nd for smear0
	  dlast = std::make_unique<RNode>(dlast->Define("indexes", [&](RVecUI idxs, RVecF Muon_pt, RVecF Muon_eta, RVecI Muon_charge, RVecF Muon_ksmear) -> RVecUI
	  {
	    unsigned int idxP = Muon_charge[idxs[0]]>0 ? idxs[0] : idxs[1];
	    unsigned int idxM = Muon_charge[idxs[0]]>0 ? idxs[1] : idxs[0];
	    float ptP  = Muon_pt[idxP];
    	float ptM  = Muon_pt[idxM];
    	float ksmear0P = Muon_ksmear[0]>0. ? Muon_ksmear[0] : 1./(pt_edges[0]-0.01);
    	float ksmear0M = Muon_ksmear[1]>0. ? Muon_ksmear[1] : 1./(pt_edges[0]-0.01);
    	float etaP = Muon_eta[idxP];
    	float etaM = Muon_eta[idxM];
    	RVecUI out;
    	out.emplace_back(n_bins);
    	out.emplace_back(n_bins);
	
    	unsigned int ibin = 0;
    	for(unsigned int ieta_p = 0; ieta_p<n_eta_bins; ieta_p++){
    	  float eta_p_low = eta_edges[ieta_p];
    	  float eta_p_up  = eta_edges[ieta_p+1];      
    	  for(unsigned int ipt_p = 0; ipt_p<n_pt_bins; ipt_p++){
    	    float pt_p_low = pt_edges[ipt_p];
    	    float pt_p_up  = pt_edges[ipt_p+1];
    	    for(unsigned int ieta_m = 0; ieta_m<n_eta_bins; ieta_m++){
    	      float eta_m_low = eta_edges[ieta_m];
    	      float eta_m_up  = eta_edges[ieta_m+1];      
    	      for(unsigned int ipt_m = 0; ipt_m<n_pt_bins; ipt_m++){
    			float pt_m_low = pt_edges[ipt_m];
    			float pt_m_up  = pt_edges[ipt_m+1];
    	      	if( etaP>=eta_p_low && etaP<eta_p_up &&
    		  	etaM>=eta_m_low && etaM<eta_m_up &&
    		  	ptP>=pt_p_low   && ptP<pt_p_up &&
    		  	ptM>=pt_m_low   && ptM<pt_m_up 
    		  	) out[0] = ibin;
    	      	if( etaP>=eta_p_low && etaP<eta_p_up &&
    		  	etaM>=eta_m_low && etaM<eta_m_up &&
    		  	1./ksmear0P>=pt_p_low   && 1./ksmear0P<pt_p_up &&
    		  	1./ksmear0M>=pt_m_low   && 1./ksmear0M<pt_m_up 
    		  	) out[1] = ibin;
    	      	ibin++;
    	      }
    	    }
    	  }
    	}
    	return out;
	  }, {"idxs", useKf ? "Muon_pt" : "Muon_cvhidealPt", useKf ? "Muon_eta" : "Muon_cvhidealEta", "Muon_charge", "Muon_ksmear"} ));
      
	  for(unsigned int r = 0 ; r<recos.size(); r++) {
        dlast = std::make_unique<RNode>(dlast->Define( TString(("index_"+recos[r]).c_str()), [r](RVecUI indexes) 
		{
	  	  return indexes.at(r);
		}, {"indexes"} ));
      }
      
	  // Define gen, reco, smear0 mass per muon pair
	  dlast = std::make_unique<RNode>(dlast->Define("masses", [&](RVecUI idxs,
								       RVecF Muon_pt, RVecF Muon_eta, RVecF Muon_phi, RVecF Muon_mass, RVecI Muon_charge,
								       UInt_t nGenPart, RVecI GenPart_status, RVecI GenPart_statusFlags, RVecI GenPart_pdgId,
								       RVecF GenPart_pt, RVecF GenPart_eta, RVecF GenPart_phi, RVecF GenPart_mass,
								       RVecF Muon_ksmear) -> RVecF
	  {
	    RVecF out;
	    unsigned int idxP = Muon_charge[idxs[0]]>0 ? idxs[0] : idxs[1];
	    unsigned int idxM = Muon_charge[idxs[0]]>0 ? idxs[1] : idxs[0];
	    ROOT::Math::PtEtaPhiMVector muP( Muon_pt[ idxP ], Muon_eta[ idxP ], Muon_phi[ idxP ], Muon_mass[ idxP ] );
	    ROOT::Math::PtEtaPhiMVector muM( Muon_pt[ idxM ], Muon_eta[ idxM ], Muon_phi[ idxM ], Muon_mass[ idxM ] );
	    // gen matching
		ROOT::Math::PtEtaPhiMVector gmuP( 0., 0., 0., 0. );
	    ROOT::Math::PtEtaPhiMVector gmuM( 0., 0., 0., 0. );
	    for(unsigned int i = 0; i < nGenPart; i++) {
	      bool isGoodGenPart = (GenPart_status[i]==1 && (GenPart_statusFlags[i] & 1 || (GenPart_statusFlags[i] & (1<<5))) && TMath::Abs(GenPart_pdgId[i])==13);
	      if(!isGoodGenPart) continue;
	      ROOT::Math::PtEtaPhiMVector gen(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
	      if( ROOT::Math::VectorUtil::DeltaR(gen, muP) < 0.1 && ROOT::Math::VectorUtil::DeltaR(gen, muM) > 0.1) gmuP = gen;
	      else if( ROOT::Math::VectorUtil::DeltaR(gen, muP) > 0.1 && ROOT::Math::VectorUtil::DeltaR(gen, muM) < 0.1) gmuM = gen;
	    }
	
	    if( gmuP.Pt()>10. && gmuM.Pt()>10.) {
	      out.emplace_back( (gmuP + gmuM).M() );
	      out.emplace_back( (muP + muM).M() );
    	  float ksmear0P = Muon_ksmear[0]>0. ? Muon_ksmear[0] : 1./(pt_edges[0]-0.01);
	      float ksmear0M = Muon_ksmear[1]>0. ? Muon_ksmear[1] : 1./(pt_edges[0]-0.01);	  
    	  ROOT::Math::PtEtaPhiMVector muP_smear0( 1./ksmear0P, Muon_eta[ idxP ], Muon_phi[ idxP ], Muon_mass[ idxP ] );
	      ROOT::Math::PtEtaPhiMVector muM_smear0( 1./ksmear0M, Muon_eta[ idxM ], Muon_phi[ idxM ], Muon_mass[ idxM ] );      
	      out.emplace_back( (muP_smear0 + muM_smear0).M() );	
	    } 
	
	    return out;
      }, {"idxs", useKf ? "Muon_pt" : "Muon_cvhidealPt", useKf ? "Muon_eta" : "Muon_cvhidealEta", useKf ? "Muon_phi" : "Muon_cvhidealPhi",
	   "Muon_mass", "Muon_charge", "nGenPart", "GenPart_status", "GenPart_statusFlags", "GenPart_pdgId", "GenPart_pt", "GenPart_eta",
	   "GenPart_phi", "GenPart_mass", "Muon_ksmear"} ));

      for(unsigned int r = 0 ; r<recos.size(); r++) {
		if(skipUnsmearedReco && recos[r]=="reco") continue;

	    unsigned int mpos = idx_map.at(recos[r]);

	    dlast = std::make_unique<RNode>(dlast->Define(TString( (recos[r]+"_m").c_str() ), [mpos](RVecF masses)
		{
	      return masses.size()>0 ? masses.at( mpos ) : -99.;
	    }, {"masses"} ));

        // Define mass - gen mass
	    dlast = std::make_unique<RNode>(dlast->Define(TString( (recos[r]+"_dm").c_str() ), [mpos](RVecF masses)
		{
	      return masses.size()>0 ? masses.at( mpos ) - masses.at(0) : -99.;
	    }, {"masses"} ));

	    dlast = std::make_unique<RNode>(dlast->Define(TString( (recos[r]+"_gm").c_str() ), [mpos](RVecF masses) 
		{
          return masses.size()>0 ? masses.at(0) : -99.;
        }, {"masses"} ));
      }
      
      // Define jacobian weights per event
      dlast = std::make_unique<RNode>(dlast->Define("weights_jac", [skipUnsmearedReco,n_bins,recos,h_map,h_jac_map,idx_map](RVecF masses, RVecUI indexes) -> RVecF
	  {
	    RVecF out;
	    if(masses.size()==0) {
	      for(unsigned int r = 0 ; r<recos.size(); r++) {
	        out.emplace_back(0.0);
	        out.emplace_back(0.0);
	        out.emplace_back(0.0);
	        out.emplace_back(0.0);
	      }
	      return out;
	    }
	
	    float gm  = masses.at(0);
	    for(unsigned int r = 0 ; r<recos.size(); r++) {
			if(skipUnsmearedReco && recos[r]=="reco") {
				out.emplace_back(0.0);
	            out.emplace_back(0.0);
	            out.emplace_back(0.0);
	            out.emplace_back(0.0);
				continue;
			}
	      unsigned int rpos = idx_map.at(recos[r]);
		  // Gaussian mean and rms of the mass - gen mass distribution in a 4D bin
	      TH1D* h_mean = h_map.at("mean_"+recos[r]);
	      TH1D* h_rms  = h_map.at("rms_"+recos[r]);
		  // Crystal Ball jacobian event weight in a 4D bin, in a mass - gen mass bin
	      TH2D* h_jac_scale = h_jac_map.at("jscale_cb_per_evt_"+recos[r]);
	      TH2D* h_jac_width = h_jac_map.at("jwidth_cb_per_evt_"+recos[r]);

	      float m = masses.at( rpos );
	      float dm = m - gm;
	      int ijac_dm = (h_jac_scale->GetYaxis()->FindBin(dm)>0 && h_jac_scale->GetYaxis()->FindBin(dm) < h_jac_scale->GetYaxis()->GetNbins()+1) ? h_jac_scale->GetYaxis()->FindBin(dm) : -99;
	  
	      float delta = 0.;
          float sigma = 0.;
	      if(indexes[r]<n_bins) {
	        delta = h_mean->GetBinContent(indexes[r]+1);
	        sigma = h_rms->GetBinContent(indexes[r]+1);
	      }
	      float jscale = sigma>0. ? +(m - (gm+delta) )*(gm+delta)/sigma/sigma : 0.0;
	      float jwidth = sigma>0. ? +(m - (gm+delta) )*(m - (gm+delta) )/sigma/sigma - 1.0 : 0.0;
	      float jscale_cb = (sigma>0. && ijac_dm>0) ? h_jac_scale->GetBinContent(indexes[r]+1, ijac_dm) : 0.0;
	      float jwidth_cb = (sigma>0. && ijac_dm>0) ? h_jac_width->GetBinContent(indexes[r]+1, ijac_dm) : 0.0;
	      out.emplace_back(jscale);
	      out.emplace_back(jwidth);
	      out.emplace_back(jscale_cb);
	      out.emplace_back(jwidth_cb);
	    }
	    return out;
      }, {"masses", "indexes"} ));

      for(unsigned int r = 0 ; r<recos.size(); r++) {
		if(skipUnsmearedReco && recos[r]=="reco") continue;

        unsigned int jpos = (idx_map.at(recos[r])-1)*4;

	    dlast = std::make_unique<RNode>(dlast->Define( TString((recos[r]+"_jscale_weight").c_str()), [jpos](RVecF weights_jac, float weight) -> float
		{
	      return weights_jac.at( jpos )*weight;
	    }, {"weights_jac", "weight"} ));

	    dlast = std::make_unique<RNode>(dlast->Define( TString((recos[r]+"_jwidth_weight").c_str()), [jpos](RVecF weights_jac, float weight) -> float
		{
	      return weights_jac.at( jpos+1 )*weight;
	    }, {"weights_jac", "weight"} ));

	    dlast = std::make_unique<RNode>(dlast->Define( TString((recos[r]+"_jscale_cb_weight").c_str()), [jpos](RVecF weights_jac, float weight) -> float
		{
          return weights_jac.at( jpos+2 )*weight;
        }, {"weights_jac", "weight"} ));
	
	    dlast = std::make_unique<RNode>(dlast->Define( TString((recos[r]+"_jwidth_cb_weight").c_str()), [jpos](RVecF weights_jac, float weight) -> float
		{
          return weights_jac.at( jpos+3 )*weight;
        }, {"weights_jac", "weight"} ));
      }
      
    }
    
    else { // data
	  // Define indices of individual muons that pass the selection
	  dlast = std::make_unique<RNode>(dlast->Define("idxs", [&](UInt_t nMuon, RVecB Muon_looseId, RVecF Muon_dxybs, RVecB Muon_isGlobal,
								RVecB Muon_highPurity, RVecB Muon_mediumId, RVecF Muon_pfRelIso04_all,
								RVecF Muon_pt, RVecF Muon_eta) -> RVecUI
	  {
	    RVecUI out;
	    for(unsigned int i = 0; i < nMuon; i++) {
	      if( Muon_looseId[i] && TMath::Abs(Muon_dxybs[i]) < 0.05 && Muon_isGlobal[i] && Muon_highPurity[i] && Muon_mediumId[i] && Muon_pfRelIso04_all[i]<0.15 &&
	      Muon_pt[i] >= pt_edges[0] && Muon_pt[i] < pt_edges[ n_pt_bins ]  && Muon_eta[i]>=eta_edges[0] && Muon_eta[i]<=eta_edges[ n_eta_bins ] ) out.emplace_back(i);
	    }
	    return out;
	  }, {"nMuon", "Muon_looseId", "Muon_dxybs", "Muon_isGlobal", "Muon_highPurity", "Muon_mediumId", "Muon_pfRelIso04_all",
	  useKf ? "Muon_pt" : "Muon_cvhPt", useKf ? "Muon_eta" : "Muon_cvhEta"} ));
      
      // Filter for muon pairs
      dlast = std::make_unique<RNode>(dlast->Filter( [](RVecUI idxs, RVecI Muon_charge, bool HLT_IsoMu24 )
	  {
	    if( idxs.size()!=2 || !HLT_IsoMu24) return false;
	    if( Muon_charge[idxs[0]]*Muon_charge[idxs[1]] > 0 ) return false;
	    return true;
      }, {"idxs", "Muon_charge", "HLT_IsoMu24"} ));      
	  
      // Define data weight = 1.0
      dlast = std::make_unique<RNode>(dlast->Define("weight", []()->float{ return 1.0; }, {} ));          

	  // Define eta+, pt+, eta-, pt- indexes for each muon pair that passed selection   
	  dlast = std::make_unique<RNode>(dlast->Define("index_data", [&](RVecUI idxs, RVecF Muon_pt, RVecF Muon_eta, RVecI Muon_charge) -> unsigned int
	  {
	    unsigned int idxP = Muon_charge[idxs[0]]>0 ? idxs[0] : idxs[1];
	    unsigned int idxM = Muon_charge[idxs[0]]>0 ? idxs[1] : idxs[0];
	    float ptP  = Muon_pt[idxP];
	    float ptM  = Muon_pt[idxM];
	    float etaP = Muon_eta[idxP];
	    float etaM = Muon_eta[idxM];
	    unsigned int out = n_bins;	
	    unsigned int ibin = 0;
	    for(unsigned int ieta_p = 0; ieta_p<n_eta_bins; ieta_p++){
	      float eta_p_low = eta_edges[ieta_p];
	      float eta_p_up  = eta_edges[ieta_p+1];      
	      for(unsigned int ipt_p = 0; ipt_p<n_pt_bins; ipt_p++){
	        float pt_p_low = pt_edges[ipt_p];
	        float pt_p_up  = pt_edges[ipt_p+1];
	        for(unsigned int ieta_m = 0; ieta_m<n_eta_bins; ieta_m++){
	          float eta_m_low = eta_edges[ieta_m];
	          float eta_m_up  = eta_edges[ieta_m+1];      
	          for(unsigned int ipt_m = 0; ipt_m<n_pt_bins; ipt_m++){
		        float pt_m_low = pt_edges[ipt_m];
		        float pt_m_up  = pt_edges[ipt_m+1];
	            if( etaP>=eta_p_low && etaP<eta_p_up &&
		        etaM>=eta_m_low && etaM<eta_m_up &&
		        ptP>=pt_p_low   && ptP<pt_p_up &&
		        ptM>=pt_m_low   && ptM<pt_m_up 
		        ) out = ibin;
	            ibin++;
	          }
	        }  
	      }
	    }
	    return out;
	  }, {"idxs", useKf ? "Muon_pt" : "Muon_cvhPt", useKf ? "Muon_eta" : "Muon_cvhEta", "Muon_charge"} ));

	  // Define mass in data    
	  dlast = std::make_unique<RNode>(dlast->Define("data_m", [&](RVecUI idxs,
								       RVecF Muon_pt, RVecF Muon_eta, RVecF Muon_phi, RVecF Muon_mass, RVecI Muon_charge) -> float
	  {
	    float out = 0.0;
	    unsigned int idxP = Muon_charge[idxs[0]]>0 ? idxs[0] : idxs[1];
	    unsigned int idxM = Muon_charge[idxs[0]]>0 ? idxs[1] : idxs[0];
	    ROOT::Math::PtEtaPhiMVector muP( Muon_pt[ idxP ], Muon_eta[ idxP ], Muon_phi[ idxP ], Muon_mass[ idxP ] );
	    ROOT::Math::PtEtaPhiMVector muM( Muon_pt[ idxM ], Muon_eta[ idxM ], Muon_phi[ idxM ], Muon_mass[ idxM ] );
    	out = (muP + muM).M();	  
    	return out;
	  }, {"idxs", useKf ? "Muon_pt" : "Muon_cvhPt", useKf ? "Muon_eta" : "Muon_cvhEta", useKf ? "Muon_phi" : "Muon_cvhPhi", "Muon_mass", "Muon_charge"} ));           
    }
    
    // Vector of pointers to histograms output by the dataframe
    std::vector< ROOT::RDF::RResultPtr<TH1D> > df_histos1D;
    std::vector< ROOT::RDF::RResultPtr<TH2D> > df_histos2D;
    std::vector< ROOT::RDF::RResultPtr<TH3D> > df_histos3D;
  
    if(iter==-1) { // Book data histogram
	  // x-axis: 4D bin index, y-axis: data mass, weight = 1
      df_histos2D.emplace_back(dlast->Histo2D({ "h_data_bin_m", "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index_data", "data_m", "weight" ));
      auto colNames = dlast->GetColumnNames();
      double total = *(dlast->Count());  
      std::cout << colNames.size() << " columns created. Total event count is " << total  << std::endl;
    }
    else if(iter==0) { // Book MC histograms
      //df_histos1D.emplace_back(dlast->Histo1D({"h_gen_m", "nominal", x_nbins, x_low, x_high}, "gen_m", "weight"));
      //df_histos1D.emplace_back(dlast->Histo1D({"h_reco_m", "nominal", x_nbins, x_low, x_high}, "reco_m", "weight"));
      //df_histos1D.emplace_back(dlast->Histo1D({"h_smear_m", "nominal", x_nbins, x_low, x_high}, "smear0_m", "weight"));
      for(unsigned int r = 0 ; r<recos.size(); r++) {
		if(skipUnsmearedReco && recos[r]=="reco") continue;
		// x-axis: 4D bin index, y-axis: MC mass, weight = MC weight
        df_histos2D.emplace_back(dlast->Histo2D({ "h_"+TString(recos[r].c_str())+"_bin_m",    "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high},   "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_m", "weight" ));
    	// x-axis: 4D bin index, y-axis: MC mass - gen mass, weight = MC weight
		df_histos2D.emplace_back(dlast->Histo2D({ "h_"+TString(recos[r].c_str())+"_bin_dm",   "nominal", n_bins, 0, double(n_bins), dm_bins, dm_low, dm_high}, "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_dm", "weight"));
    	//df_histos3D.emplace_back(dlast->Histo3D({ "h_"+TString(recos[r].c_str())+"_bin_gm_dm", "nominal", n_bins, 0, double(n_bins),  x_nbins, x_low, x_high, dm_bins, dm_low, dm_high}, "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_gm", TString(recos[r].c_str())+"_dm", "weight"));
    	//df_histos3D.emplace_back(dlast->Histo3D({ "h_"+TString(recos[r].c_str())+"_bin_gm_m", "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high, x_nbins, x_low, x_high},     "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_gm", TString(recos[r].c_str())+"_m", "weight"));
      }
      auto colNames = dlast->GetColumnNames();
      double total = *(dlast->Count());  
      std::cout << colNames.size() << " columns created. Total event count is " << total  << std::endl;
    }
    else if(iter==1) { // Book jac histograms only for smear0
      for(unsigned int r = 0 ; r<recos.size(); r++){
	    if(skipUnsmearedReco && recos[r]=="reco") continue;
		// x-axis: 4D bin index, y-axis: MC mass, weight = gaussian scale jacobian event weight
    	df_histos2D.emplace_back(dlast->Histo2D({"h_"+TString(recos[r].c_str())+"_bin_jac_scale", "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_m", TString(recos[r].c_str())+"_jscale_weight"));
    	// x-axis: 4D bin index, y-axis: MC mass, weight = gaussian width jacobian event weight
		df_histos2D.emplace_back(dlast->Histo2D({"h_"+TString(recos[r].c_str())+"_bin_jac_width", "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_m", TString(recos[r].c_str())+"_jwidth_weight"));
    	// x-axis: 4D bin index, y-axis: MC mass, weight = Crystal Ball scale jacobian event weight
		df_histos2D.emplace_back(dlast->Histo2D({"h_"+TString(recos[r].c_str())+"_bin_jac_scale_cb", "cb", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_m", TString(recos[r].c_str())+"_jscale_cb_weight"));
    	// x-axis: 4D bin index, y-axis: MC mass, weight = Crystal Ball width jacobian event weight
		df_histos2D.emplace_back(dlast->Histo2D({"h_"+TString(recos[r].c_str())+"_bin_jac_width_cb", "cb", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_m", TString(recos[r].c_str())+"_jwidth_cb_weight"));
      }
    }

    // Write dataframe histograms
    if(iter<2) {
	  fout->cd();
	  std::cout << "Writing histos..." << std::endl;
	  
	  // Scale MC to luminosity in data
	  double lumiMC = lumiMC2016;
	  if(y2017)      lumiMC = lumiMC2017;
	  else if(y2018) lumiMC = lumiMC2018;
	  
	  double sf = lumi>0. ? lumi/lumiMC : 1.0; //double(lumi)/double(minNumEvents);
	  
	  for(auto h : df_histos1D) {
		if(iter>=0) h->Scale(sf); // scale only for MC
		h->Write();
	  }
	  for(auto h : df_histos2D) {
		if(iter>=0) h->Scale(sf); // scale only for MC
		string h_name = std::string(h->GetName());
		std::cout << "Total number of events in 2D histo " << h_name << ": " << h->GetEntries() << std::endl;
		h->Write();
      }
	  for(auto h : df_histos3D) {
       	if(iter>=0) h->Scale(sf); // scale only for MC
       	string h_name = std::string(h->GetName());
       	std::cout << "Total number of events in 3D histo " << h_name << ": " << h->GetEntries() << std::endl;
       	h->Write();
      }
    }
    
    if(iter==0) { 

      cout << "Writing aux files" << endl;
      h_pt_edges->Write();
      h_eta_edges->Write();
      h_A_vals_nom->Write();
      h_e_vals_nom->Write();
      h_M_vals_nom->Write();
      h_A_vals_prevfit->Write();
      h_e_vals_prevfit->Write();
      h_M_vals_prevfit->Write();
      h_c_vals_prevfit->Write();
      h_d_vals_prevfit->Write();

	  // Fill histograms using the results from the dataframe
      
      RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
      gErrorIgnoreLevel = 6001;

      for(unsigned int r = 0 ; r<recos.size(); r++) {

		if(skipUnsmearedReco && recos[r]=="reco") continue;
	
	    TH2D* h_reco_dm = (TH2D*)fout->Get(TString( ("h_"+recos[r]+"_bin_dm").c_str()) );
	    TH2D* h_reco_m  = (TH2D*)fout->Get(TString( ("h_"+recos[r]+"_bin_m").c_str()) );
	    if( h_reco_dm==0 || h_reco_m==0 ) {
	      cout << "h_reco_dm/h_reco_m NOT FOUND" << endl;
	      continue;
	    }
		// Gaussian mean and rms of the mass - gen mass distribution in a 4D bin
	    h_map["mean_"+recos[r]] = new TH1D( TString( ("h_mean_"+recos[r]+"_bin_dm").c_str() ),"", n_bins, 0, double(n_bins));
	    h_map["rms_"+recos[r]]  = new TH1D( TString( ("h_rms_"+recos[r]+"_bin_dm").c_str() ),"", n_bins, 0, double(n_bins));
	    // 1/0 if keeping(ignoring) a 4D bin in the fit
		h_map["mask_"+recos[r]] = new TH1D( TString( ("h_mask_"+recos[r]+"_bin_dm").c_str() ),"", n_bins, 0, double(n_bins));
	
	    // Crystal Ball jacobian event weights in a 4D bin, in a mass - gen mass bin
	    h_jac_map["jscale_cb_per_evt_"+recos[r]] = new TH2D("h_"+TString(recos[r].c_str())+"_bin_jscale_cb_per_evt", "cb", n_bins, 0, double(n_bins), dm_bins, dm_low, dm_high); 
	    h_jac_map["jwidth_cb_per_evt_"+recos[r]] = new TH2D("h_"+TString(recos[r].c_str())+"_bin_jwidth_cb_per_evt", "cb", n_bins, 0, double(n_bins), dm_bins, dm_low, dm_high);
	
	    for(unsigned int i = 0; i<n_bins; i++ ) {
	      if(i%1000==0) cout << "Doing gaus fit for 4D bin " << i << " / " << n_bins << endl;
	      TString projname(Form("bin_%d_", i));
	      projname += TString( recos[r].c_str() );
	      TH1D* hi   = (TH1D*)h_reco_dm->ProjectionY( projname+"_dm", i+1, i+1 );
	      TH1D* hi_m = (TH1D*)h_reco_m->ProjectionY( projname+"_m", i+1, i+1 );
	      double mean_i = 0.0;
	      double meanerr_i = 0.0;
	      double rms_i = 0.0;
	      double rmserr_i = 0.0;
	      //cout << hi_m->Integral() << ", " << hi->Integral() << ", " << hi_m->GetMean() << endl;
		  
		  // 4D bin selection cuts
	      if( hi_m->Integral() > minNumEvents && hi->Integral() > minNumEvents  &&  hi_m->GetMean()>( x_low + 5.0 ) && hi_m->GetMean()<( x_high - 5.0 ) ) { //TODO make this 5.0 an input parameter
	        h_map.at("mask_"+recos[r])->SetBinContent(i+1, 1);

            // Gaus fit
	        TF1* gf = new TF1("gf","[0]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp( -0.5*(x-[1])*(x-[1])/[2]/[2] )",
			      hi->GetXaxis()->GetBinLowEdge(1), hi->GetXaxis()->GetBinUpEdge( hi->GetXaxis()->GetNbins() ));      
	        gf->SetParameter(0, hi->Integral());
	        gf->SetParameter(1, hi->GetMean());
	        gf->SetParameter(2, hi->GetRMS() );
	        float m_min = nRMSforGausFit>0. ? TMath::Max(-nRMSforGausFit*hi->GetRMS(), dm_low) : dm_low;
	        float m_max = nRMSforGausFit>0. ? TMath::Min(+nRMSforGausFit*hi->GetRMS(), dm_high) : dm_high;
	        hi->Fit("gf", "QR", "", m_min, m_max );
	        mean_i    = gf->GetParameter(1);
	        meanerr_i = gf->GetParError(1);
	        rms_i     = TMath::Abs(gf->GetParameter(2));
	        rmserr_i  = gf->GetParError(2);
			if(maxRMS>0. && rms_i>maxRMS) h_map.at("mask_"+recos[r])->SetBinContent(i+1, 0);
	        //cout << "Fit " << mean_i << endl;
	        delete gf;
            
	        // Crystal Ball fit
	        RooRealVar x0("x0", "", mean_i, dm_low, dm_high); 
			RooRealVar mass("mass", "", dm_low, dm_high); 
	        mass.setRange("r1", dm_low, dm_high); 
	        RooRealVar alphaL("alphaL", "", 1.0, 0.2, +10 );
	        RooRealVar alphaR("alphaR", "", 1.0, 0.2, +10 );
	        RooRealVar nL("nL", "", 2, 1, 100 );
	        RooRealVar nR("nR", "", 2, 1, 100 );
	        RooRealVar sigmaL("sigmaL", "", rms_i, rms_i*0.5, rms_i*2 );
	        RooRealVar sigmaR("sigmaR", "", rms_i, rms_i*0.5, rms_i*2 );
	    				 
	        RooDataHist data("data", "", RooArgList(mass), hi );
	        RooCrystalBall pdf("pdf", "", mass, x0, sigmaL, sigmaR, alphaL, nL, alphaR, nR);
	    
	        std::unique_ptr<RooFitResult> res{pdf.fitTo(data,
					      InitialHesse(true),
					      Minimizer("Minuit2"),
					      Range("r1"),
					      Save(), SumW2Error(true),
						  PrintLevel(-1),
						  Verbose(false) )};
	    
	        TH1D* h_der = new TH1D("h_der", "", hi->GetXaxis()->GetNbins()*2, dm_low, dm_high);
	        h_der->Reset();
	        RooDerivative* der = pdf.derivative( mass, 1, 0.001 );
	    
		    // Crystal Ball jacobian event weights in a 4D bin, in a mass - gen mass bin
	        for(int ib=1; ib<=h_der->GetXaxis()->GetNbins();ib++) {
	          double x = h_der->GetXaxis()->GetBinCenter(ib);
	          mass.setVal( x );
	          double fprime = der->getVal();
	          double f = pdf.getVal();
	          h_jac_map.at("jscale_cb_per_evt_"+recos[r])->SetBinContent(i+1, ib, -fprime/f * hi_m->GetMean());
	          h_jac_map.at("jwidth_cb_per_evt_"+recos[r])->SetBinContent(i+1, ib, -(1+x*fprime/f));
	        }
            
	      }
	      else {
	        h_map.at("mask_"+recos[r])->SetBinContent(i+1, 0);
	      }
	      h_map.at("mean_"+recos[r])->SetBinContent(i+1, mean_i);
	      h_map.at("mean_"+recos[r])->SetBinError(i+1, meanerr_i);
	      h_map.at("rms_"+recos[r])->SetBinContent(i+1, rms_i);
	      h_map.at("rms_"+recos[r])->SetBinError(i+1, rmserr_i);
	    } 
      }
      
      fout->cd();
    
      for(unsigned int r = 0 ; r<recos.size(); r++) {	 
		if(skipUnsmearedReco && recos[r]=="reco") continue;

	    h_map["mean_"+recos[r]]->Write();
	    h_map["rms_"+recos[r]]->Write();
	    h_map["mask_"+recos[r]]->Write();
      }
    
	}

    else if(iter==2) {
      if(saveMassFitHistos && fout->GetDirectory("postfit")==0) fout->mkdir("postfit");

      // Make tree with quantities relevant to the fit and the results
      TTree* treescales = new TTree("treescales","treescales");
      int inmassbins, indof, ibinIdx  ;
      float inevents, ibeta, ibetaErr, ialpha, ialphaErr, inu, inuErr, iprob, ichi2old, ichi2new; 
      treescales->Branch("nevents",&inevents,"nevents/F");
      treescales->Branch("beta",&ibeta,"beta/F");
      treescales->Branch("betaErr",&ibetaErr,"betaErr/F");
      treescales->Branch("alpha",&ialpha,"alpha/F");
      treescales->Branch("alphaErr",&ialphaErr,"alphaErr/F");
      treescales->Branch("nu",&inu,"nu/F");
      treescales->Branch("nuErr",&inuErr,"nuErr/F");
      treescales->Branch("prob",&iprob,"prob/F");
      treescales->Branch("chi2old",&ichi2old,"chi2old/F"); // prefit agreement between data and MC
      treescales->Branch("chi2new",&ichi2new,"chi2new/F"); // postfit agreement between data and MC
      treescales->Branch("nmassbins",&inmassbins,"nmassbins/I");
      treescales->Branch("ndof",&indof,"ndof/I");
      treescales->Branch("binIdx",&ibinIdx,"binIdx/I");
      
	  // Define histograms to save results, x-axis is 4D bin index
      TH1D* h_scales  = new TH1D("h_scales", "", n_bins, 0, double(n_bins));
      TH1D* h_widths  = new TH1D("h_widths", "", n_bins, 0, double(n_bins));
      TH1D* h_norms   = new TH1D("h_norms", "", n_bins, 0, double(n_bins));
      TH1D* h_probs   = new TH1D("h_probs", "", n_bins, 0, double(n_bins));
      TH1D* h_masks   = new TH1D("h_masks", "", n_bins, 0, double(n_bins));

      // Get histograms needed for the mass fit
      TH2D* h_data_2D   = (TH2D*)fout->Get("h_data_bin_m");
      TH2D* h_nom_2D    = (TH2D*)fout->Get("h_smear0_bin_m");
      TH1D* h_nom_mask  = (TH1D*)fout->Get("h_mask_smear0_bin_dm");
      TH2D* h_jscale_2D = useCB ? (TH2D*)fout->Get("h_smear0_bin_jac_scale_cb") : (TH2D*)fout->Get("h_smear0_bin_jac_scale");
      TH2D* h_jwidth_2D = useCB ? (TH2D*)fout->Get("h_smear0_bin_jac_width_cb") : (TH2D*)fout->Get("h_smear0_bin_jac_width");
      
      for(unsigned int ibin=0; ibin<n_bins; ibin++) { // Loop over 4D bins

	    // skip empty bins
	    if( h_nom_mask->GetBinContent(ibin+1)<0.5 ) continue;

	    ibinIdx = ibin;
	
	    TH1D* h_data_i    = (TH1D*)h_data_2D->ProjectionY( Form("h_data_i_%d",ibin ),     ibin+1, ibin+1 );
	    TH1D* h_nom_i     = (TH1D*)h_nom_2D->ProjectionY( Form("h_nom_i_%d", ibin),       ibin+1, ibin+1 );
	    TH1D* h_jscale_i  = (TH1D*)h_jscale_2D->ProjectionY( Form("h_jscale_i_%d", ibin), ibin+1, ibin+1 );
	    TH1D* h_jwidth_i  = (TH1D*)h_jwidth_2D->ProjectionY( Form("h_jwidth_i_%d", ibin), ibin+1, ibin+1 );	

	    if(scaleToData) {
	      double data_norm_i = h_data_i->Integral();
	      double mc_norm_i = h_nom_i->Integral();
	      if(data_norm_i>0. && mc_norm_i>0.) {
	        double sf_i = data_norm_i/mc_norm_i;
	        h_nom_i->Scale(sf_i);
	        h_jscale_i->Scale(sf_i);
	        h_jwidth_i->Scale(sf_i);
	      }
	    }
	
	    if(rebin>1) {
	      h_data_i->Rebin(rebin);
	      h_nom_i->Rebin(rebin);
	      h_jscale_i->Rebin(rebin);
	      h_jwidth_i->Rebin(rebin);
	    }
	
	    unsigned int n_mass_bins = 0;

	    // Skip 4D bins with less than minNumMassBins high-stat (>minNumEventsPerBin) data mass bins
	    for(int im = 1 ; im<=h_data_i->GetXaxis()->GetNbins(); im++) {
	      if( h_data_i->GetBinContent(im)>minNumEventsPerBin ) n_mass_bins++;
	    }
	    if( n_mass_bins < minNumMassBins ) {
	      h_scales->SetBinContent(ibin+1, 0.0);
	      h_widths->SetBinContent(ibin+1, 0.0);
	      h_norms->SetBinContent(ibin+1, 0.0);
	      h_probs->SetBinContent(ibin+1, 0.0);
	      h_masks->SetBinContent(ibin+1, 0.0);
	      continue;
	    }

	    inevents = h_data_i->Integral();
	    inmassbins = n_mass_bins;
	
	    // Get mass fit terms
        MatrixXd inv_sqrtV(n_mass_bins,n_mass_bins);
	    MatrixXd inv_V(n_mass_bins,n_mass_bins);
	    for(unsigned int ibm = 0; ibm<n_mass_bins; ibm++ ) {
	      for(unsigned int jbm = 0; jbm<n_mass_bins; jbm++ ) {
	        inv_sqrtV(ibm,jbm) = 0.;
	        inv_V(ibm,jbm) = 0.;
	      }
	    }
	    VectorXd y(n_mass_bins);
	    VectorXd y0(n_mass_bins);
	    VectorXd jscale(n_mass_bins);
	    VectorXd jwidth(n_mass_bins);
	    unsigned int bin_counter = 0;
	    for(int im = 0 ; im<h_data_i->GetXaxis()->GetNbins(); im++) {
	      if( h_data_i->GetBinContent(im+1)>minNumEventsPerBin ) {
	        y(bin_counter)  = h_data_i->GetBinContent(im+1);
	        y0(bin_counter) = h_nom_i->GetBinContent(im+1);	    
	        jscale(bin_counter) = h_jscale_i->GetBinContent(im+1);
	        jwidth(bin_counter) = h_jwidth_i->GetBinContent(im+1);  
	        double mcErr_im = h_nom_i->GetBinError(im+1);
	        inv_V(bin_counter,bin_counter) = lumi>0. ?
	        1./(y(bin_counter)  + mcErr_im*mcErr_im ) :
	        1./(2*mcErr_im*mcErr_im);
	        //cout << TMath::Sqrt(y(bin_counter)) << " (+) " << h_nom_i->GetBinError(im+1) << endl;
	        inv_sqrtV(bin_counter,bin_counter) = TMath::Sqrt( inv_V(bin_counter,bin_counter) );
	        bin_counter++;
	      }
	    }

	    unsigned int n_fit_params = 3;
	    if(!fitWidth) n_fit_params--;
	    if(!fitNorm)  n_fit_params--;
	    MatrixXd jac(n_mass_bins, n_fit_params);
	    for(unsigned int ib=0; ib<n_mass_bins;ib++){
	      jac(ib, 0) = jscale(ib);
	      if(fitWidth) jac(ib, 1) = jwidth(ib);
	      if(fitNorm)  jac(ib, 2) = y0(ib);
	    }

        // Mass fit
	    MatrixXd A = inv_sqrtV*jac;
	    VectorXd b = inv_sqrtV*(y-y0);
	    VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
	    MatrixXd C = (jac.transpose()*inv_V*jac).inverse();
	    MatrixXd rho( C.rows(), C.rows() ) ;
	    for(unsigned int ir = 0; ir<C.rows(); ir++) {
	      for(unsigned int ic = 0; ic<C.rows(); ic++) {
	        rho(ir,ic) = C(ir,ic)/TMath::Sqrt(C(ir,ir)*C(ic,ic));
	      }
	    }
	    MatrixXd chi2old = b.transpose()*b;
	    MatrixXd chi2new = ((b - A*x).transpose())*(b-A*x);
	    int ndof = n_mass_bins-n_fit_params;
	    double chi2norm_old = chi2old(0,0)/(n_mass_bins);
	    double chi2norm_new = chi2new(0,0)/ndof;
	    double prob = TMath::Prob(chi2norm_new*ndof, ndof );
	    //cout << "Bin: " << ibin << ": mass fits with " << n_mass_bins << " mass bins: " << chi2norm_old << " ---> " << chi2norm_new << " (prob = " << prob << ")" << endl;

	    indof = ndof;
	    ibeta = x(0);
	    ibetaErr = TMath::Sqrt(C(0,0));
	    ialpha = fitWidth ? x(1) : 0.;
	    ialphaErr = fitWidth ? TMath::Sqrt(C(1,1)) : 0.;
	    inu = fitNorm ? x(2) : 0.;
	    inuErr = fitNorm ? TMath::Sqrt(C(2,2)) : 0.;
 	    ichi2old = chi2norm_old;
	    ichi2new = chi2norm_new;
	    iprob =  prob;	
	    treescales->Fill();
	    //cout << "Filling tree" << endl;
	
        h_scales->SetBinContent(ibin+1, ibeta+1.0);
        h_scales->SetBinError(ibin+1, ibetaErr);
	    h_norms->SetBinContent(ibin+1, inu+1.0);
	    h_norms->SetBinError(ibin+1, inuErr);
	    h_widths->SetBinContent(ibin+1, ialpha+1.0);
	    h_widths->SetBinError(ibin+1, ialphaErr);
	    h_probs->SetBinContent(ibin+1, prob);
	    h_probs->SetBinError(ibin+1, 0.);
	    h_masks->SetBinContent(ibin+1, 1.0);

        // Optional: save pre and postfit mass distribution in 4D bin
	    if(saveMassFitHistos) {
	      TH1D* h_pre_i   = (TH1D*)h_nom_i->Clone(Form("h_prefit_%d", ibin));
	      TH1D* h_post_i  = (TH1D*)h_nom_i->Clone(Form("h_postfit_%d", ibin));
	      unsigned int bin_counter = 0;
	      for(int im = 0 ; im<h_post_i->GetXaxis()->GetNbins(); im++) {	  
	        if( h_data_i->GetBinContent(im+1)>minNumEventsPerBin ) {
	          h_post_i->SetBinContent( im+1, y0(bin_counter)+(jac*x)(bin_counter) );
	          bin_counter++;
	        }
	      }
	      fout->cd("postfit/");
	      h_data_i->Write(Form("h_data_%d", ibin) ,TObject::kOverwrite);
	      h_pre_i->Write(TString(h_pre_i->GetName()) ,TObject::kOverwrite);
	      h_post_i->Write(TString(h_post_i->GetName()),TObject::kOverwrite);
	    }	
      }
      
      fout->cd();
      h_scales->Write(0,TObject::kOverwrite);
      h_norms->Write(0,TObject::kOverwrite);
      h_widths->Write(0,TObject::kOverwrite);
      h_probs->Write(0,TObject::kOverwrite);
      h_masks->Write(0,TObject::kOverwrite);
      treescales->Write(0,TObject::kOverwrite);
	
      cout << h_masks->Integral() << " scales have been computed" << endl;
    }

  }
  
  sw.Stop();

  std::cout << "Real time: " << sw.RealTime()/60. << " mins " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;

  fout->Close(); 
  
  return 0;
}
