#include <ROOT/RDataFrame.hxx>
#include "TFile.h"
#include "TRandom3.h"
#include "TVector.h"
#include "TVectorT.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraphErrors.h"
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
using namespace ROOT::Minuit2;

typedef ROOT::VecOps::RVec<double> RVecD;
typedef ROOT::VecOps::RVec<unsigned int> RVecUI;
typedef ROOT::VecOps::RVec<int> RVecI;
typedef ROOT::VecOps::RVec<float> RVecF;
typedef ROOT::VecOps::RVec<bool> RVecB;
using ROOT::RDF::RNode;


using namespace boost::program_options;

constexpr double lumiMC = 3.33369e+08/2001.9e+03;
  
int main(int argc, char* argv[]) {

  TStopwatch sw;
  sw.Start();

  ROOT::EnableImplicitMT();

  variables_map vm;
  try {
    options_description desc{"Options"};
    desc.add_options()
	    ("help,h", "Help screen")
	    ("minNumEvents",       value<int>()->default_value(100), "min number of events for a histogram to be accepted")
	    ("minNumEventsPerBin", value<int>()->default_value(10), "min number of events for a bin of a histogram to be accepted")
	    ("lumi",               value<float>()->default_value(16.1), "number of events in data")
	    ("tag",                value<std::string>()->default_value("closure"), "run type, type of toy used")
	    ("run",                value<std::string>()->default_value("closure"), "number of iteration")
	    ("skipUnsmearedReco",  bool_switch()->default_value(true), "do not fill histograms of mass or jacobian weights per event for reco without smearing")
	    ("saveMassFitHistos",  bool_switch()->default_value(false), "save pre and postfit mass distribution in 4D bin")
	    ("firstIter",          value<int>()->default_value(0), "firstIter (0: fills MC and pseudodata histos, 1(needs 0): fills jacobians, 2(needs 0,1): fits for scale/resolution bias)")
	    ("lastIter",           value<int>()->default_value(2), "lastIter  (0: fills MC and pseudodata histos, 1(needs 0): fills jacobians, 2(needs 0,1): fits for scale/resolution bias)")
	    ("biasResolution",     value<float>()->default_value(-1.), "input resolution bias randomly sampled from -biasResolution to biasResolution")
	    ("nRMSforGausFit",     value<float>()->default_value(-1.), "number of RMS for Gaus mass difference fit")
	    ("minNumMassBins",     value<int>()->default_value(4), "min number of mass bins for a histogram to be accepted")
	    ("maxRMS",             value<float>()->default_value(-1.), "max RMS of Gaus mass difference fit for a 4D bin to be included in the mass fit")
	    ("rebin",              value<int>()->default_value(2), "rebin before fit")
	    ("fitWidth",           bool_switch()->default_value(false), "compute resolution bias")
	    ("fitNorm",            bool_switch()->default_value(false), "compute difference in normalisation in 4D bin")
	    ("usePrevMassFit",     bool_switch()->default_value(false), "use previous mass fit")
	    ("tagPrevMassFit",     value<std::string>()->default_value("closure"), "run type, type of toy used")
	    ("runPrevMassFit",     value<std::string>()->default_value("closure"), "number of iteration")
	    ("usePrevResolFit",    bool_switch()->default_value(false), "use previous resolution fit")
	    ("tagPrevResolFit",    value<std::string>()->default_value("closure"), "run type, type of toy used")
	    ("runPrevResolFit",    value<std::string>()->default_value("closure"), "number of iteration")
	    ("useMCasData",        bool_switch()->default_value(false), "use MC as data")
	    ("seed",               value<int>()->default_value(4357), "seed for random numbers");

    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);
    if (vm.count("help")) {
	    std::cout << desc << '\n';
	    return 0;
    }
    if (vm.count("tag"))        std::cout << "Tag: " << vm["tag"].as<std::string>() << '\n';
    if (vm.count("run"))        std::cout << "Run: " << vm["run"].as<std::string>() << '\n';
  }
  catch (const error &ex) {
    std::cerr << ex.what() << '\n';
  }
  
  int minNumEvents       = vm["minNumEvents"].as<int>();
  float lumi             = vm["lumi"].as<float>();
  float nRMSforGausFit   = vm["nRMSforGausFit"].as<float>();
  float biasResolution   = vm["biasResolution"].as<float>();
  std::string tag        = vm["tag"].as<std::string>();
  std::string run        = vm["run"].as<std::string>();
  int seed               = vm["seed"].as<int>();
  int minNumEventsPerBin = vm["minNumEventsPerBin"].as<int>();
  int minNumMassBins     = vm["minNumMassBins"].as<int>();
  int rebin              = vm["rebin"].as<int>();
  int firstIter          = vm["firstIter"].as<int>();
  int lastIter           = vm["lastIter"].as<int>();
  bool skipUnsmearedReco = vm["skipUnsmearedReco"].as<bool>();
  bool saveMassFitHistos = vm["saveMassFitHistos"].as<bool>();
  bool fitWidth          = vm["fitWidth"].as<bool>();
  bool fitNorm           = vm["fitNorm"].as<bool>();
  bool usePrevMassFit    = vm["usePrevMassFit"].as<bool>();
  bool usePrevResolFit   = vm["usePrevResolFit"].as<bool>();
  std::string tagPrevMassFit  = vm["tagPrevMassFit"].as<std::string>();
  std::string runPrevMassFit  = vm["runPrevMassFit"].as<std::string>();
  std::string tagPrevResolFit = vm["tagPrevResolFit"].as<std::string>();
  std::string runPrevResolFit = vm["runPrevResolFit"].as<std::string>();
  bool useMCasData            = vm["useMCasData"].as<bool>();
  float maxRMS                = vm["maxRMS"].as<float>();
    
  TRandom3* ran0 = new TRandom3(seed);

  // Sample random bias in resolution
  if( biasResolution < 0. )
    biasResolution = 0.0;
  else
    biasResolution = ran0->Uniform(-biasResolution,+biasResolution);
  cout << "Bias in resolution: " << biasResolution << endl;
  
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

  // _nom histograms with the input curvature scale bias parameters A, e or M for the current toy
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
  
  float kmean_val = 0.5*( 1./pt_edges[0] + 1./pt_edges[ pt_edges.size()-1] );
  // Store input curvature biases
  VectorXd A_vals_nom( n_eta_bins );
  VectorXd e_vals_nom( n_eta_bins );
  VectorXd M_vals_nom( n_eta_bins );
  // Store curvature (not pT) biases
  VectorXd A_vals_fit( n_eta_bins );
  VectorXd e_vals_fit( n_eta_bins );
  VectorXd M_vals_fit( n_eta_bins );
  VectorXd c_vals_fit( n_eta_bins );
  VectorXd d_vals_fit( n_eta_bins );

  // Generate random curvature scale bias parameters A,e,M
  // Initialize curvature bias parameters from previous iteration A,e,M,c,d = 0, they will remain 0 if usePrevMassFit/usePrevResolFit are false
  
  // Bias for A out
  double val_A_rnd = ran0->Uniform(-0.001, 0.001);
  for(unsigned int i=0; i<n_eta_bins; i++) {
    //double val = ran0->Uniform(-0.001, 0.001);
    // Random, eta dependent input curvature scale bias parameter A
    float y_max = h_eta_edges->GetXaxis()->GetXmax();
    float y_i = h_eta_edges->GetXaxis()->GetBinCenter(i+1); 
    double val = val_A_rnd*( 1 + (y_i/y_max)*(y_i/y_max) );
    A_vals_nom(i) = val;
    h_A_vals_nom->SetBinContent(i+1, val);
    // Prevfit
    h_A_vals_prevfit->SetBinContent(i+1, 0.0);
    h_c_vals_prevfit->SetBinContent(i+1, 0.0);
    h_d_vals_prevfit->SetBinContent(i+1, 0.0);
    A_vals_fit(i) = 0.0;
  }
  // Bias for e out
  for(unsigned int i=0; i<n_eta_bins; i++) {
    // Random input curvature scale bias parameter e
    double val = ran0->Uniform(-0.0001/kmean_val, 0.0001/kmean_val);
    //double val = 0.0;
    e_vals_nom(i) = val;
    h_e_vals_nom->SetBinContent(i+1, val);
    // Prevfit
    h_e_vals_prevfit->SetBinContent(i+1, 0.0);
    e_vals_fit(i) = 0.0;
  }
  // Bias for M out
  for(unsigned int i=0; i<n_eta_bins; i++) {
    // Random input curvature scale bias parameter M
    double val = ran0->Uniform(-0.001*kmean_val, 0.001*kmean_val);
    //double val = 0.0;
    M_vals_nom(i) = val;
    h_M_vals_nom->SetBinContent(i+1, val);
    // Prevfit
    h_M_vals_prevfit->SetBinContent(i+1, 0.0);
    M_vals_fit(i) = 0.0;
  }            

  // Keep track of MC before(reco) and after(smear0) applying curvature corrections from previous iterations
  // smear1 is the pseudodata
  // If skipUnsmearedReco == false, we will save mass and jacobian histograms for both reco and smear0 in the dataframe
  // In the scale fit, only smear0 jacobians are used
  std::vector<string> recos = {"reco", "smear0", "smear1"}; 

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

  // Map to positions in RVecF "masses", position 0 is gen
  std::map<string, unsigned int> idx_map;
  idx_map.insert( std::make_pair<string, unsigned int >("reco",   1 ) );
  idx_map.insert( std::make_pair<string, unsigned int >("smear0", 2 ) );
  idx_map.insert( std::make_pair<string, unsigned int >("smear1", 3 ) );

  // Work out the biased resolution to smear MC curvature to get pseudodata
  TH1D* histobudget = 0;
  TH1D* histohitres = 0;
  TFile* faux = TFile::Open("root/coefficients2016ptfrom20forscaleptfrom20to70forres.root", "READ");
  if(faux!=0) {
    histobudget = (TH1D*)faux->Get("histobudget");
    histohitres = (TH1D*)faux->Get("histohitres");
  } else {
    cout << "No nominal resolution file found! Will quit" << endl;
    return;
  }

  auto resolution = [histobudget, histohitres](float k, float eta, float bias) -> float
  {
    float out = 0.02*k*(1.0 + bias); // Output if input histos are empty
    if(histobudget!=0 && histohitres!=0) {
      int eta_bin = histobudget->FindBin(eta);
      if(eta_bin<1) eta_bin = 1;
      else if(eta_bin>histobudget->GetNbinsX()) eta_bin = histobudget->GetNbinsX();
      float budget = histobudget->GetBinContent( eta_bin );
      float hitres = histohitres->GetBinContent( eta_bin );
      float budget2 = budget*budget;
      float hitres2 = hitres*hitres;
      // Nominal (pT resolution divided by pT)*k
      out = TMath::Sqrt( budget2 + hitres2/k/k )*k;
      // Bias
      out *= (1.0 + bias);
      return out;
    }
    return out;
  };

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
	      // Update the nom histograms
	      h_A_vals_nom->SetBinContent(i+1, h_A_vals_nom->GetBinContent(i+1) - A_vals_fit(i) );
	      h_e_vals_nom->SetBinContent(i+1, h_e_vals_nom->GetBinContent(i+1) - e_vals_fit(i) );
	      h_M_vals_nom->SetBinContent(i+1, h_M_vals_nom->GetBinContent(i+1) - M_vals_fit(i) );
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
      for(unsigned int i=0; i<n_eta_bins; i++) {
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
  // If firstIter = 2, update an existing output file with iter 0 and 1 to (over)write iter 2 (the mass fit results)
  TFile* fout = TFile::Open(("./massscales_"+tag+"_"+run+".root").c_str(), firstIter<2 ? "RECREATE" : "UPDATE");
  
  // iter  0 -> MC and pseudodata mass histos + calculation of jacobian terms per event
  // iter  1 -> MC jacobian histos
  // iter  2 -> fit for beta(= mass scale bias) [,alpha(= mass resolution bias), nu(= mass normalization bias)]

  // By the end the output file will contain: 
  // - histograms resulted from the event-loop at different iterations
  // - TTree and histograms resulted from the mass fits
  // - OPTIONAL, in the postfit/ folder, pre and postfit mass distribution for the 4D bins

  for(unsigned int iter=0; iter<3; iter++) {

    if( !(iter>=firstIter && iter<=lastIter) ) continue;
    cout << "Doing iter " << iter << endl;

    // Define dataframe from the MC files
    ROOT::RDataFrame d( "Events",
			{//"/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_040854/0000/NanoV9MCPostVFP_1.root",
			  "/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_040854/0000/NanoV9MCPostVFP_*.root",
			  "/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_041233/0000/NanoV9MCPostVFP_*.root",
			  "/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_041233/0001/NanoV9MCPostVFP_*.root",
			  "/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_041233/0002/NanoV9MCPostVFP_*.root"
			} );

    // Define vector of different TRandom variables to be used by different threads
    unsigned int nslots = d.GetNSlots();
    std::vector<TRandom3*> rans = {};
    for(unsigned int i = 0; i < nslots; i++){
      rans.emplace_back( new TRandom3(seed + i*10 + iter) );
    }

    auto dlast = std::make_unique<RNode>(d);
    
    // Define MC reco weight
    dlast = std::make_unique<RNode>(dlast->Define("weight_reco", [](float weight) -> float
    {
      return std::copysign(1.0, weight);
    }, {"Generator_weight"} ));

    // Define MC smear0 weight -> even events are used as MC
    dlast = std::make_unique<RNode>(dlast->Define("weight_smear0", [](ULong64_t rdfentry, float weight) -> float
    {
      return std::copysign(1.0, weight)*(rdfentry%2==0);
    }, {"rdfentry_", "Generator_weight"} ));

    // Define MC smear1 weight -> odd events are used as pseudodata
    dlast = std::make_unique<RNode>(dlast->Define("weight_smear1", [](ULong64_t rdfentry, float weight) -> float
    {
      return std::copysign(1.0, weight)*(rdfentry%2==1);
    }, {"rdfentry_", "Generator_weight"} ));
    
    // Define the indices of individual muons passing selection criteria
    dlast = std::make_unique<RNode>(dlast->Define("idxs", [&](UInt_t nMuon, RVecB Muon_looseId, RVecF Muon_dxybs, RVecB Muon_isGlobal, RVecB Muon_highPurity,
                                    RVecB Muon_mediumId, RVecF Muon_pfRelIso04_all, RVecF Muon_pt, RVecF Muon_eta) -> RVecUI
    {
      RVecUI out;
      for(unsigned int i = 0; i < nMuon; i++) {
	      if( Muon_looseId[i] && TMath::Abs(Muon_dxybs[i]) < 0.05 && Muon_isGlobal[i] && Muon_highPurity[i] && Muon_mediumId[i] && Muon_pfRelIso04_all[i]<0.15 &&
	          Muon_pt[i] >= pt_edges[0] && Muon_pt[i] < pt_edges[ n_pt_bins ]  && Muon_eta[i]>=eta_edges[0] && Muon_eta[i]<=eta_edges[ n_eta_bins ]  ) {
	        out.emplace_back(i);
	      }
      }
      return out;
    }, {"nMuon", "Muon_looseId", "Muon_dxybs", "Muon_isGlobal", "Muon_highPurity","Muon_mediumId", "Muon_pfRelIso04_all",	"Muon_cvhPt", "Muon_cvhEta" } ));

    // Filter to keep only events with exactly 2 oppositely charged, selected muons
    dlast = std::make_unique<RNode>(dlast->Filter( [](RVecUI idxs, RVecI Muon_charge, bool HLT_IsoMu24 )
    {
      if( idxs.size()!=2 || !HLT_IsoMu24) return false;
      if( Muon_charge[idxs[0]]*Muon_charge[idxs[1]] > 0 ) return false;
      return true;
    }, {"idxs", "Muon_charge", "HLT_IsoMu24"} ));

    // Define smear0 -> used as MC: pos and neg curvature k smeared according to the curvature biases A,e,M,c,d computed in previous iterations
    // Define smear1 -> used as pseudodata: pos and neg curvature k smeared according to the input curvature biases A,e,M,c,d 
    dlast = std::make_unique<RNode>(dlast->DefineSlot("Muon_ksmear", [&](unsigned int nslot, RVecUI idxs, RVecF Muon_pt, RVecF Muon_eta, RVecF Muon_phi,
                                    RVecF Muon_mass, RVecI Muon_charge, UInt_t nGenPart, RVecI GenPart_status, RVecI GenPart_statusFlags, RVecI GenPart_pdgId,
									                  RVecF GenPart_pt, RVecF GenPart_eta, RVecF GenPart_phi, RVecF GenPart_mass) -> RVecF
    {
      RVecF out;
      unsigned int idxP = Muon_charge[idxs[0]]>0 ? idxs[0] : idxs[1];
      unsigned int idxM = Muon_charge[idxs[0]]>0 ? idxs[1] : idxs[0];
      // reco muon
      ROOT::Math::PtEtaPhiMVector muP( Muon_pt[ idxP ], Muon_eta[ idxP ], Muon_phi[ idxP ], Muon_mass[ idxP ] );
      ROOT::Math::PtEtaPhiMVector muM( Muon_pt[ idxM ], Muon_eta[ idxM ], Muon_phi[ idxM ], Muon_mass[ idxM ] );
      // gen matching
      ROOT::Math::PtEtaPhiMVector gmuP;
      ROOT::Math::PtEtaPhiMVector gmuM;
      for(unsigned int i = 0; i < nGenPart; i++) {
	      bool isGoodGenPart = (GenPart_status[i]==1 && (GenPart_statusFlags[i] & 1 || (GenPart_statusFlags[i] & (1<<5))) && TMath::Abs(GenPart_pdgId[i])==13);
	      if(!isGoodGenPart) continue;
	      ROOT::Math::PtEtaPhiMVector gen(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
	      if( ROOT::Math::VectorUtil::DeltaR(gen, muP) < 0.1 && ROOT::Math::VectorUtil::DeltaR(gen, muM) > 0.1) gmuP = gen;
	      else if( ROOT::Math::VectorUtil::DeltaR(gen, muP) > 0.1 && ROOT::Math::VectorUtil::DeltaR(gen, muM) < 0.1) gmuM = gen;
      }
      
      // gen pT cut
      if( gmuP.Pt()>10. && gmuM.Pt()>10.) {

        // gen curvature
	      float kmuP = 1./gmuP.Pt();
	      float kmuM = 1./gmuM.Pt();

	      float scale_smear0P = 1.0;
	      float scale_smear0M = 1.0;
	      float resol_smear0P = 0.0;
	      float resol_smear0M = 0.0;
	      float scale_smear1P = 1.0;
	      float scale_smear1M = 1.0;

        // AeMcd are eta dependent
	      unsigned int ietaP = n_eta_bins;
	      for(unsigned int ieta_p = 0; ieta_p<n_eta_bins; ieta_p++) {
	        float eta_p_low = eta_edges[ieta_p];
	        float eta_p_up  = eta_edges[ieta_p+1];
	        if( gmuP.Eta()>=eta_p_low && gmuP.Eta()<eta_p_up) ietaP = ieta_p;	  
	      }
	      unsigned int ietaM = n_eta_bins;
	      for(unsigned int ieta_m = 0; ieta_m<n_eta_bins; ieta_m++) {
	        float eta_m_low = eta_edges[ieta_m];
	        float eta_m_up  = eta_edges[ieta_m+1];
	        if( gmuM.Eta()>=eta_m_low && gmuM.Eta()<eta_m_up) ietaM = ieta_m;	  
	      }	

	      if(ietaP<n_eta_bins && ietaM<n_eta_bins) {
          // Correct the MC curvature with the curvature scale biases derived in previous iterations (which are respectively equal to (-1)* sum of the pT scale biases from previous iterations)
	        // if usePrevMassFit is false, A,e,M are 0
          scale_smear0P = (1. + A_vals_fit(ietaP) + e_vals_fit(ietaP)*kmuP - M_vals_fit(ietaP)/kmuP);
	        scale_smear0M = (1. + A_vals_fit(ietaM) + e_vals_fit(ietaM)*kmuM + M_vals_fit(ietaM)/kmuM);
          // Generate pseudodata curvature as gen k smeared according to the input curvature biases A,e,M
	        scale_smear1P = (1. + A_vals_nom(ietaP) + e_vals_nom(ietaP)*kmuP - M_vals_nom(ietaP)/kmuP);
	        scale_smear1M = (1. + A_vals_nom(ietaM) + e_vals_nom(ietaM)*kmuM + M_vals_nom(ietaM)/kmuM);
	        //cout << "smear0:" << scale_smear0P << ": " << 1 << " + " << A_vals_fit(ietaP) << " + " << e_vals_fit(ietaP)*kmuP << " - " << M_vals_fit(ietaP)/kmuP << endl;
	        //cout << "smear1:" << scale_smear1P << ": " << 1 << " + " << A_vals_nom(ietaP) << " + " << e_vals_nom(ietaP)*kmuP << " - " << M_vals_nom(ietaP)/kmuP << endl;

	        if(usePrevResolFit) {
            // Correct the MC curvature with the resolution biases derived in previous iterations (which are respectively equal to the sum of the resolution biases from previous iterations)
	          resol_smear0P = TMath::Sqrt( TMath::Max( 1.0 + c_vals_fit(ietaP) + d_vals_fit(ietaP)*kmuP, 0.0)  ) - 1.0;
	          resol_smear0M = TMath::Sqrt( TMath::Max( 1.0 + c_vals_fit(ietaM) + d_vals_fit(ietaM)*kmuM, 0.0)  ) - 1.0;
	          //cout << "resol_smear0P: sqrt( max(1.0 + " << c_vals_fit(ietaP)  << " + " << d_vals_fit(ietaP)*kmuP << ")) - 1.0 = " << resol_smear0P << endl;
	          //cout << "resol_smear0M: sqrt( max(1.0 + " << c_vals_fit(ietaM)  << " + " << d_vals_fit(ietaM)*kmuM << ")) - 1.0 = " << resol_smear0M << endl;  
	        }
	      }

	      float resol0P = resolution(kmuP, gmuP.Eta(), resol_smear0P);
	      float resol0M = resolution(kmuM, gmuM.Eta(), resol_smear0M);
	      float resol1P = resolution(kmuP, gmuP.Eta(), biasResolution);
	      float resol1M = resolution(kmuM, gmuM.Eta(), biasResolution);

	      if(useMCasData) {
          //smear0 
	        out.emplace_back( (1./gmuP.Pt() + (1./muP.Pt() - 1./gmuP.Pt())*(1.0 + resol_smear0P))*scale_smear0P  );
	        out.emplace_back( (1./gmuM.Pt() + (1./muM.Pt() - 1./gmuM.Pt())*(1.0 + resol_smear0M))*scale_smear0M  );
	        //smear1 
          out.emplace_back( (1./gmuP.Pt() + (1./muP.Pt() - 1./gmuP.Pt())*(1.0 + biasResolution))*scale_smear1P  );
	        out.emplace_back( (1./gmuM.Pt() + (1./muM.Pt() - 1./gmuM.Pt())*(1.0 + biasResolution))*scale_smear1M  );
	      } else {
	        //smear0 -> used as MC
          out.emplace_back( rans[nslot]->Gaus(kmuP*scale_smear0P, resol0P) );
	        out.emplace_back( rans[nslot]->Gaus(kmuM*scale_smear0M, resol0M) );
	        //smear1 -> used as pseudodata
          out.emplace_back( rans[nslot]->Gaus(kmuP*scale_smear1P, resol1P) );
	        out.emplace_back( rans[nslot]->Gaus(kmuM*scale_smear1M, resol1M) );
	      }
      } else {
	      //smear0
	      out.emplace_back(0.0);
	      out.emplace_back(0.0);
	      //smear1
	      out.emplace_back(0.0);
	      out.emplace_back(0.0);
      }
      return out;
    }, {"idxs",	"Muon_cvhPt", "Muon_cvhEta", "Muon_cvhPhi", "Muon_mass", "Muon_charge", "nGenPart", "GenPart_status",
     "GenPart_statusFlags", "GenPart_pdgId", "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass"} ));
    
    // Define eta+, pt+, eta-, pt- indexes for each muon pair that passed selection. The 1st entry in "indexes" is for reco, the 2nd for smear0, 3rd smear 1
    dlast = std::make_unique<RNode>(dlast->Define("indexes", [&](RVecUI idxs, RVecF Muon_pt, RVecF Muon_eta, RVecI Muon_charge, RVecF Muon_ksmear) -> RVecUI
    {
      unsigned int idxP = Muon_charge[idxs[0]]>0 ? idxs[0] : idxs[1];
      unsigned int idxM = Muon_charge[idxs[0]]>0 ? idxs[1] : idxs[0];
      float ptP  = Muon_pt[idxP];
      float ptM  = Muon_pt[idxM];
      float ksmear0P = Muon_ksmear[0]>0. ? Muon_ksmear[0] : 1./(pt_edges[0]-0.01);
      float ksmear0M = Muon_ksmear[1]>0. ? Muon_ksmear[1] : 1./(pt_edges[0]-0.01);
      float ksmear1P = Muon_ksmear[2]>0. ? Muon_ksmear[2] : 1./(pt_edges[0]-0.01);
      float ksmear1M = Muon_ksmear[3]>0. ? Muon_ksmear[3] : 1./(pt_edges[0]-0.01);
      float etaP = Muon_eta[idxP];
      float etaM = Muon_eta[idxM];
      RVecUI out;
      out.emplace_back(n_bins);
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
	            if( etaP>=eta_p_low && etaP<eta_p_up &&
		          etaM>=eta_m_low && etaM<eta_m_up &&
		          1./ksmear1P>=pt_p_low   && 1./ksmear1P<pt_p_up &&
		          1./ksmear1M>=pt_m_low   && 1./ksmear1M<pt_m_up 
		          ) out[2] = ibin;	      
	            ibin++;
	          }
	        }
	      }
      }
      return out;
    }, {"idxs", "Muon_cvhPt", "Muon_cvhEta", "Muon_charge", "Muon_ksmear"} ));

    for(unsigned int r = 0 ; r<recos.size(); r++) {
      dlast = std::make_unique<RNode>(dlast->Define( TString(("index_"+recos[r]).c_str()), [r](RVecUI indexes)
      {
	      return indexes.at(r);
      }, {"indexes"} ));
    }
    
    // Define gen, reco, smear0, smear1 mass per muon pair
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
	      float ksmear1P = Muon_ksmear[2]>0. ? Muon_ksmear[2] : 1./(pt_edges[0]-0.01);
	      float ksmear1M = Muon_ksmear[3]>0. ? Muon_ksmear[3] : 1./(pt_edges[0]-0.01);

	      ROOT::Math::PtEtaPhiMVector muP_smear0( 1./ksmear0P, Muon_eta[ idxP ], Muon_phi[ idxP ], Muon_mass[ idxP ] );
	      ROOT::Math::PtEtaPhiMVector muM_smear0( 1./ksmear0M, Muon_eta[ idxM ], Muon_phi[ idxM ], Muon_mass[ idxM ] );      
	      out.emplace_back( (muP_smear0 + muM_smear0).M() );	

	      ROOT::Math::PtEtaPhiMVector muP_smear1( 1./ksmear1P, Muon_eta[ idxP ], Muon_phi[ idxP ], Muon_mass[ idxP ] );
	      ROOT::Math::PtEtaPhiMVector muM_smear1( 1./ksmear1M, Muon_eta[ idxM ], Muon_phi[ idxM ], Muon_mass[ idxM ] );      
	      out.emplace_back( (muP_smear1 + muM_smear1).M() );
      }
      
      return out;
    }, {"idxs",	"Muon_cvhPt", "Muon_cvhEta", "Muon_cvhPhi", "Muon_mass", "Muon_charge", "nGenPart", "GenPart_status",
     "GenPart_statusFlags", "GenPart_pdgId", "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass",	"Muon_ksmear"} ));


    dlast = std::make_unique<RNode>(dlast->Define("gen_m", [](RVecF masses)
    {
      return masses.size()>0 ? masses.at(0) : -99.;
    }, {"masses"} ));

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
    }
    
    // Define jacobian weights per event
    dlast = std::make_unique<RNode>(dlast->Define("weights_jac", [n_bins,recos,h_map,idx_map](RVecF masses, RVecUI indexes) -> RVecF
    {
      RVecF out;
      if(masses.size()==0) {
	      for(unsigned int r = 0 ; r<recos.size(); r++) {
	        out.emplace_back(0.0);
	        out.emplace_back(0.0);
	      }
	      return out;
      }
      
      float gen_m  = masses.at(0);
      for(unsigned int r = 0 ; r<recos.size(); r++) {
        if(skipUnsmearedReco && recos[r]=="reco") {
				  out.emplace_back(0.0);
	        out.emplace_back(0.0);
				  continue;
			  }
	      unsigned int rpos = idx_map.at(recos[r]);
	      // Gaussian mean and rms of the mass - gen mass distribution in a 4D bin
        TH1D* h_mean = h_map.at("mean_"+recos[r]);
	      TH1D* h_rms  = h_map.at("rms_"+recos[r]);

	      float reco_m = masses.at( rpos );	
	      float reco_delta = 0.;
	      float reco_sigma = 0.;
	      if(indexes[r]<n_bins) {
	        reco_delta = h_mean->GetBinContent(indexes[r]+1);
	        reco_sigma = h_rms->GetBinContent(indexes[r]+1);
	      }
	      float reco_jscale = reco_sigma>0. ? +(reco_m - (gen_m+reco_delta) )*(gen_m+reco_delta)/reco_sigma/reco_sigma : 0.0;
	      float reco_jwidth = reco_sigma>0. ? +(reco_m - (gen_m+reco_delta) )*(reco_m - (gen_m+reco_delta) )/reco_sigma/reco_sigma - 1.0 : 0.0;
	      out.emplace_back(reco_jscale);
	      out.emplace_back(reco_jwidth);
      }
      return out;
    }, {"masses", "indexes"} ));

    for(unsigned int r = 0 ; r<recos.size(); r++) {
      if(skipUnsmearedReco && recos[r]=="reco") continue;

      unsigned int jpos = (idx_map.at(recos[r])-1)*2;
      
      dlast = std::make_unique<RNode>(dlast->Define( TString((recos[r]+"_jscale_weight").c_str()), [jpos](RVecF weights_jac, float weight) -> float
      {
	      //unisgned int jpos = (idx_map[recos[r]]-1)*2;
	      return weights_jac.at( jpos )*weight;
      }, {"weights_jac", "weight_"+recos[r] } ));
      
      dlast = std::make_unique<RNode>(dlast->Define( TString((recos[r]+"_jwidth_weight").c_str()), [jpos](RVecF weights_jac, float weight) -> float 
      {
	      //unisgned int jpos = (idx_map[recos[r]]-1)*2 + 1;
	      return weights_jac.at( jpos+1 )*weight;
      }, {"weights_jac", "weight_"+recos[r] } ));
    }
    
    // Vector of pointers to histograms output by the dataframe
    std::vector<ROOT::RDF::RResultPtr<TH1D> > histos1D;
    std::vector<ROOT::RDF::RResultPtr<TH2D> > histos2D;
  
    //histos1D.emplace_back(dlast->Histo1D({"h_gen_m", "nominal", x_nbins, x_low, x_high}, "gen_m", "weight"));      
    //histos1D.emplace_back(dlast->Histo1D({"h_reco_m", "nominal", x_nbins, x_low, x_high}, "reco_m", "weight"));
    //histos1D.emplace_back(dlast->Histo1D({"h_smear_m", "nominal", x_nbins, x_low, x_high}, "smear_m", "weight"));

    if(iter==0) { // Book MC and pseudodata histograms
      if(!skipUnsmearedReco)	 
        histos2D.emplace_back(dlast->Histo2D({"h_gen_bin_m", "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index_reco", "gen_m", "weight_reco"));
      for(unsigned int r = 0 ; r<recos.size(); r++) {
	      if(skipUnsmearedReco && recos[r]=="reco") continue;
	      histos2D.emplace_back(dlast->Histo2D({ "h_"+TString(recos[r].c_str())+"_bin_m",    "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_m", "weight_"+recos[r] ));
	      histos2D.emplace_back(dlast->Histo2D({ "h_"+TString(recos[r].c_str())+"_bin_dm",   "nominal", n_bins, 0, double(n_bins), 24, -6.0, 6.0},          "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_dm", "weight_"+recos[r]));
      }
      auto colNames = dlast->GetColumnNames();
      double total = *(dlast->Count());  
      std::cout << colNames.size() << " columns created. Total event count is " << total  << std::endl;
    }
    else if(iter==1) { // Book jac histograms
      for(unsigned int r = 0 ; r<recos.size(); r++) {
	      if(skipUnsmearedReco && recos[r]=="reco") continue;
	      histos2D.emplace_back(dlast->Histo2D({"h_"+TString(recos[r].c_str())+"_bin_jac_scale", "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_m", TString(recos[r].c_str())+"_jscale_weight"));
	      histos2D.emplace_back(dlast->Histo2D({"h_"+TString(recos[r].c_str())+"_bin_jac_width", "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_m", TString(recos[r].c_str())+"_jwidth_weight"));
      }
    }
    
    // Write dataframe histograms
    fout->cd();
    std::cout << "Writing histos..." << std::endl;
    
    // Scale MC to luminosity in data
    double sf = lumi>0. ? lumi/lumiMC : 1.0; //double(lumi)/double(minNumEvents);
    
    for(auto h : histos1D) {
      h->Scale(sf);
      string h_name = std::string(h->GetName());
      h->Write();
    }
    for(auto h : histos2D) {
      h->Scale(sf);
      string h_name = std::string(h->GetName());
      std::cout << "Total number of events in 2D histo " << h_name << ": " << h->GetEntries() << std::endl;
      h->Write();
    }    
    std::cout << "Total slots: " << dlast->GetNSlots() << std::endl;

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
          if( hi_m->Integral() > minNumEvents && hi->Integral() > minNumEvents  &&  hi_m->GetMean()>75. && hi_m->GetMean()<105. ){
	          h_map.at("mask_"+recos[r])->SetBinContent(i+1, 1);
	    
            // Gaus fit
            TF1* gf = new TF1("gf","[0]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp( -0.5*(x-[1])*(x-[1])/[2]/[2] )",
			            hi->GetXaxis()->GetBinLowEdge(1), hi->GetXaxis()->GetBinUpEdge( hi->GetXaxis()->GetNbins() ));      
	          gf->SetParameter(0, hi->Integral());
	          gf->SetParameter(1, hi->GetMean());
	          gf->SetParameter(2, hi->GetRMS() );
	          float m_min = nRMSforGausFit>0. ? TMath::Max(-nRMSforGausFit*hi->GetRMS(), -6.0) : -6.0;
	          float m_max = nRMSforGausFit>0. ? TMath::Min(+nRMSforGausFit*hi->GetRMS(), +6.0) : +6.0;
	          hi->Fit("gf", "QR", "", m_min, m_max );
	          mean_i    = gf->GetParameter(1);
	          meanerr_i = gf->GetParError(1);
	          rms_i     = TMath::Abs(gf->GetParameter(2));
	          rmserr_i  = gf->GetParError(2);
	          if(maxRMS>0. && rms_i>maxRMS) h_map.at("mask_"+recos[r])->SetBinContent(i+1, 0);
	          //cout << "Fit " << mean_i << endl;
	          delete gf;
	        } else {
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

    if(iter==2) {
      if(saveMassFitHistos &&GetDirectory("postfit")==0) fout->mkdir("postfit");
      
      // Make tree with quantities relevant to the fit and the results
      TTree* treescales = new TTree("treescales","treescales");
      int inmassbins, indof, ibinIdx  ;
      float inevents, ibeta, ibetaErr, ialpha, ialphaErr, inu, inuErr, iprob, ichi2old, ichi2new; 
      float imean, irms;
      treescales->Branch("nevents",&inevents,"nevents/F");
      treescales->Branch("beta",&ibeta,"beta/F");
      treescales->Branch("betaErr",&ibetaErr,"betaErr/F");
      treescales->Branch("alpha",&ialpha,"alpha/F");
      treescales->Branch("alphaErr",&ialphaErr,"alphaErr/F");
      treescales->Branch("nu",&inu,"nu/F");
      treescales->Branch("nuErr",&inuErr,"nuErr/F");
      treescales->Branch("prob",&iprob,"prob/F");
      treescales->Branch("chi2old",&ichi2old,"chi2old/F");
      treescales->Branch("chi2new",&ichi2new,"chi2new/F");
      treescales->Branch("nmassbins",&inmassbins,"nmassbins/I");
      treescales->Branch("ndof",&indof,"ndof/I");
      treescales->Branch("binIdx",&ibinIdx,"binIdx/I");
      treescales->Branch("mean",&imean,"mean/F");
      treescales->Branch("rms", &irms,"rms/F");
      
      // Define histograms to save results, x-axis is 4D bin index
      TH1D* h_scales  = new TH1D("h_scales", "", n_bins, 0, double(n_bins));
      TH1D* h_widths  = new TH1D("h_widths", "", n_bins, 0, double(n_bins));
      TH1D* h_norms   = new TH1D("h_norms", "", n_bins, 0, double(n_bins));
      TH1D* h_probs   = new TH1D("h_probs", "", n_bins, 0, double(n_bins));
      TH1D* h_masks   = new TH1D("h_masks", "", n_bins, 0, double(n_bins));

      // Get histograms needed for the mass fit
      TH2D* h_data_2D   = (TH2D*)fout->Get("h_smear1_bin_m");
      TH2D* h_nom_2D    = (TH2D*)fout->Get("h_smear0_bin_m");
      TH2D* h_nomdm_2D  = (TH2D*)fout->Get("h_smear0_bin_dm");
      TH1D* h_nom_mask  = (TH1D*)fout->Get("h_mask_smear0_bin_dm");
      TH2D* h_jscale_2D = (TH2D*)fout->Get("h_smear0_bin_jac_scale");
      TH2D* h_jwidth_2D = (TH2D*)fout->Get("h_smear0_bin_jac_width");

      for(unsigned int ibin=0; ibin<n_bins; ibin++) { // Loop over 4D bins

	      // skip empty bins
	      if( h_nom_mask->GetBinContent(ibin+1)<0.5 ) continue;

	      ibinIdx = ibin;
	
	      TH1D* h_data_i    = (TH1D*)h_data_2D->ProjectionY( Form("h_data_i_%d",ibin ),     ibin+1, ibin+1 );
	      TH1D* h_nom_i     = (TH1D*)h_nom_2D->ProjectionY( Form("h_nom_i_%d", ibin),       ibin+1, ibin+1 );
	      TH1D* h_jscale_i  = (TH1D*)h_jscale_2D->ProjectionY( Form("h_jscale_i_%d", ibin), ibin+1, ibin+1 );
	      TH1D* h_jwidth_i  = (TH1D*)h_jwidth_2D->ProjectionY( Form("h_jwidth_i_%d", ibin), ibin+1, ibin+1 );
	      TH1D* h_nomdm_i   = (TH1D*)h_nomdm_2D->ProjectionY( Form("h_nomdm_i_%d", ibin),   ibin+1, ibin+1 );

	      imean = h_nomdm_i->GetMean();
	      irms  = h_nomdm_i->GetRMS();
	
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
	      if( n_mass_bins<minNumMassBins ) {
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
	      for(unsigned int ib=0; ib<n_mass_bins;ib++) {
	        jac(ib, 0) = jscale(ib);
	        if(fitWidth) jac(ib, 1) = jwidth(ib);
	        if(fitNorm) jac(ib, 2) = y0(ib);
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
	      //cout << "Bin: " << ibin << ": mass fits with " << n_mass_bins << " mass bins: " << chi2norm_old << " ---> " << chi2norm_new << " (prob = " << TMath::Prob(chi2norm_new*ndof, ndof ) << ")" << endl;

	      indof = ndof;
	      ibeta = x(0);
	      ibetaErr = TMath::Sqrt(C(0,0));
	      ialpha = fitWidth ? x(1) : 0.;
	      ialphaErr = fitWidth ? TMath::Sqrt(C(1,1)) : 0.;
	      inu = fitNorm ? x(2) : 0.;
	      inuErr = fitNorm ? TMath::Sqrt(C(2,2)) : 0.;
 	      ichi2old = chi2norm_old;
	      ichi2new = chi2norm_new;
	      iprob =  TMath::Prob(chi2norm_new*ndof, ndof );	
	      treescales->Fill();
	      //cout << "Filling tree" << endl;
	
	      h_scales->SetBinContent(ibin+1, ibeta+1.0);
	      h_scales->SetBinError(ibin+1, ibetaErr);
	      h_norms->SetBinContent(ibin+1, inu+1.0);
	      h_norms->SetBinError(ibin+1, inuErr);
	      h_widths->SetBinContent(ibin+1, ialpha+1.0);
	      h_widths->SetBinError(ibin+1, ialphaErr);
	      h_probs->SetBinContent(ibin+1, TMath::Prob(chi2norm_new*ndof, ndof ));
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

    for(auto r : rans) delete r;
  }
  
  sw.Stop();

  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;

  fout->Close(); 

  delete ran0;
  
  return 0;
}
