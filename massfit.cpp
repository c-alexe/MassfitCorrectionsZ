// It can be ran in data mode ->  takes the mass scale biases per 4D bin and fits for the pT scale biases parameters A,e,M per eta bin 
// OR toys mode (closure test) -> generates mass scale biases from dummy AeM biases and fits for AeM from them
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
#include <TMatrixD.h>
#include <TMatrixDSymfwd.h>
#include <TStopwatch.h>
#include <ROOT/RVec.hxx>
#include <iostream>
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
using ROOT::RDF::RNode; 

using namespace boost::program_options;

constexpr double MZ = 91.;
constexpr double GW = 2.5;

class TheoryFcn : public FCNGradientBase {
//class TheoryFcn : public FCNBase {

public:
  TheoryFcn(const int& debug, const int& seed, const int& bias, string fname)
    : errorDef_(1.0), debug_(debug), seed_(seed), bias_(bias)
  {

    ran_ = new TRandom3(seed_);

    // pT and eta binning
    if(bias_==-1) { // data mode, matches binning in masscales_data.cpp
      TFile* fin = TFile::Open(fname.c_str(), "READ");
      if(fin==0) {
	      cout << "No data file found! Will quit" << endl;
  	    return;
      }
      else {
	      cout << string(fin->GetName()) << " found!" << endl;
      }
      TH1F* h_pt_edges = (TH1F*)fin->Get("h_pt_edges");
      TH1F* h_eta_edges = (TH1F*)fin->Get("h_eta_edges");
      unsigned int pt_edges_size  = h_pt_edges->GetXaxis()->GetNbins()+1;
      unsigned int eta_edges_size = h_eta_edges->GetXaxis()->GetNbins()+1;
      
      pt_edges_.reserve( pt_edges_size );
      for(unsigned int i=0; i<h_pt_edges->GetXaxis()->GetNbins(); i++) 
        pt_edges_.push_back( h_pt_edges->GetXaxis()->GetBinLowEdge(i+1) );
      pt_edges_.push_back( h_pt_edges->GetXaxis()->GetBinUpEdge( h_pt_edges->GetXaxis()->GetNbins() ));

      eta_edges_.reserve( eta_edges_size );
      for(unsigned int i=0; i<h_eta_edges->GetXaxis()->GetNbins(); i++)
        eta_edges_.push_back( h_eta_edges->GetXaxis()->GetBinLowEdge(i+1) );
      eta_edges_.push_back( h_eta_edges->GetXaxis()->GetBinUpEdge( h_eta_edges->GetXaxis()->GetNbins() ));
      fin->Close();
    }
    else { // toys mode 
      pt_edges_  = {25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0}; 
      eta_edges_ = {-2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0,
                    0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4};
    }

    n_pt_bins_  = pt_edges_.size()-1;
    n_eta_bins_ = eta_edges_.size()-1;

    n_pars_ = 3*n_eta_bins_;
    
    for(auto p : pt_edges_) k_edges_.emplace_back(1./p); // curvature k
    for(unsigned int i = 0; i < n_pt_bins_; i++) kmean_vals_.emplace_back( 0.5*(k_edges_[i]+k_edges_[i+1]) );
    kmean_val_ = 0.5*(kmean_vals_[n_pt_bins_-1] + kmean_vals_[0]);

    n_data_ = n_eta_bins_*n_eta_bins_*n_pt_bins_*n_pt_bins_; // number of 4D bins
    n_dof_ = 0;
    
    // Prepare storage for fit inputs and results
    scales2_.reserve(n_data_); // mass scale bias -> (beta + 1.0)^2
    scales2Err_.reserve(n_data_);
    masks_.reserve(n_data_); // 1/0 if keeping(ignoring) a 4D bin in the fit
    for(unsigned int idata = 0; idata<n_data_; idata++) {
      scales2_.push_back( 0.0 );
      scales2Err_.push_back( 0.0 );
      masks_.push_back(1);
    }

    // Toys mode: input pT scale biases A,e,M
    x_vals_ = VectorXd(n_pars_); // Holds all input A,e,M biases in a single vector in this order
    A_vals_ = VectorXd(n_eta_bins_);
    e_vals_ = VectorXd(n_eta_bins_);
    M_vals_ = VectorXd(n_eta_bins_);
    // Sum of the pT scale biases A, e or M from all the previous iterations
    A_vals_prevfit_ = VectorXd(n_eta_bins_);
    e_vals_prevfit_ = VectorXd(n_eta_bins_);
    M_vals_prevfit_ = VectorXd(n_eta_bins_);
    for(unsigned int i=0; i<n_eta_bins_; i++) {
      A_vals_(i) = 0.0; 
      e_vals_(i) = 0.0;
      M_vals_(i) = 0.0;
      A_vals_prevfit_(i) = 0.0;
      e_vals_prevfit_(i) = 0.0;
      M_vals_prevfit_(i) = 0.0;
    }
    for(unsigned int i=0; i<n_pars_; i++){
      x_vals_(i) = 0.0;
    }

    if(bias_>0) { // Toys mode only 
      // bias for A out
      for(unsigned int i=0; i<n_eta_bins_; i++) {
	      double val = ran_->Uniform(-0.001, 0.001);
	      if (bias_== 2) {
	        double mid_point = double(n_eta_bins_)*0.5;
	        val = (i-mid_point)*(i-mid_point)/mid_point/mid_point*0.001;
	      }
	      A_vals_(i) = val;
	      x_vals_(i) = val;
      }
      // bias for e out
      for(unsigned int i=0; i<n_eta_bins_; i++) {
	      double val = ran_->Uniform(-0.0001/kmean_val_, 0.0001/kmean_val_);
	      if (bias_== 2) {
	        double mid_point = double(n_eta_bins_)*0.5;
	        val = -(i-mid_point)*(i-mid_point)/mid_point/mid_point*0.0001;
	      }
	      e_vals_(i) = val;
	      x_vals_(i+n_eta_bins_) = val;
      }
      // bias for M out
      for(unsigned int i=0; i<n_eta_bins_; i++) {
	      double val = ran_->Uniform(-0.001*kmean_val_, 0.001*kmean_val_);
	      if (bias_== 2) {
	        double mid_point = double(n_eta_bins_)*0.5;
	        val = (i-mid_point)/mid_point*0.001;
	      }
	      M_vals_(i) = val;
	      x_vals_(i+2*n_eta_bins_) = val;
      }        
      n_dof_ = n_data_ - n_pars_;    
    }
    else if(bias_== -1) { // Data mode
      // Read scales and masks per 4D bin from input file
      TFile* fin = TFile::Open(fname.c_str(), "READ");
      if(fin==0) {
	      cout << "No data file found! Will quit" << endl;
  	    return;
      }
      TH1D* h_scales = (TH1D*)fin->Get("h_scales"); // mass scale bias -> beta + 1.0
      TH1D* h_masks = (TH1D*)fin->Get("h_masks");
      assert( h_scales->GetXaxis()->GetNbins() == n_data_);

      unsigned int n_unmasked_bins = 0;  
      for(unsigned int ibin=0;ibin<h_scales->GetXaxis()->GetNbins(); ibin++) {
	      scales2_[ibin]    = h_scales->GetBinContent(ibin+1)*h_scales->GetBinContent(ibin+1);
	      scales2Err_[ibin] = 2*TMath::Abs(h_scales->GetBinContent(ibin+1))*h_scales->GetBinError(ibin+1);
	      masks_[ibin]      = h_masks->GetBinContent(ibin+1);
	      if( masks_[ibin]>0.5 ) n_unmasked_bins++;
      }

      n_dof_ = n_unmasked_bins - n_pars_;
      n_data_ = n_unmasked_bins;
	
      // AeM values used to generate toy in massscales.cpp OR 0 in massscales_data.cpp
      TH1D* h_A_vals = (TH1D*)fin->Get("h_A_vals_nom");
      TH1D* h_e_vals = (TH1D*)fin->Get("h_e_vals_nom");
      TH1D* h_M_vals = (TH1D*)fin->Get("h_M_vals_nom");
      // Read from massscales_data.cpp the sum of pT scale biases A,e or M obtained from all the previous iterations
      TH1D* h_A_vals_prevfit = (TH1D*)fin->Get("h_A_vals_prevfit");
      TH1D* h_e_vals_prevfit = (TH1D*)fin->Get("h_e_vals_prevfit");
      TH1D* h_M_vals_prevfit = (TH1D*)fin->Get("h_M_vals_prevfit");
      
      assert( h_A_vals->GetXaxis()->GetNbins() == n_eta_bins_ );
      assert( h_e_vals->GetXaxis()->GetNbins() == n_eta_bins_ );
      assert( h_M_vals->GetXaxis()->GetNbins() == n_eta_bins_ );
	
      for(unsigned int i=0; i<n_eta_bins_; i++) {
	      A_vals_(i) = -h_A_vals->GetBinContent(i+1);
	      A_vals_prevfit_(i) = h_A_vals_prevfit->GetBinContent(i+1);
	      x_vals_(i) = A_vals_(i);
      }
      for(unsigned int i=0; i<n_eta_bins_; i++) {
	      e_vals_(i) = -h_e_vals->GetBinContent(i+1);
	      e_vals_prevfit_(i) = h_e_vals_prevfit->GetBinContent(i+1);
	      x_vals_(i+n_eta_bins_) = e_vals_(i);
      }
      for(unsigned int i=0; i<n_eta_bins_; i++) {
	      M_vals_(i) = -h_M_vals->GetBinContent(i+1);
	      M_vals_prevfit_(i) = h_M_vals_prevfit->GetBinContent(i+1);
	      x_vals_(i+2*n_eta_bins_) = M_vals_(i);
      }
      fin->Close();
    }

    // Transformation of external to internal parameters
    U_ = MatrixXd(n_pars_,n_pars_);
    for(unsigned int i=0; i<n_pars_; i++) {
      for(unsigned int j=0; j<n_pars_; j++) {
	      // block(A,A)
	      if(i<n_eta_bins_ && j<n_eta_bins_)
	        U_(i,j) = i==j ? 1.0 : 0.0;
	      // block(A,e)
	      else if(i<n_eta_bins_ && (j>=n_eta_bins_ && j<2*n_eta_bins_) )
	        U_(i,j) = i==(j-n_eta_bins_) ? -1.0*kmean_val_ : 0.0;
	      // block(e,e)
	      else if(i>=n_eta_bins_ && i<2*n_eta_bins_ && j>=n_eta_bins_ && j<2*n_eta_bins_)
	        U_(i,j) = i==j ? kmean_val_ : 0.0;
	      // block(M,M)
	      else if(i>=2*n_eta_bins_  && j>=2*n_eta_bins_)
	        U_(i,j) = i==j ? 1.0/kmean_val_ : 0.0;
	      else U_(i,j) = 0.0;
      }
    }
    //cout << U_ << endl;
    
  }
  
  ~TheoryFcn() { delete ran_;}

  // In toy mode, function to generate mass scale bias^2 values from given AeM
  void generate_data();

  // Function to set the seed value for random numbers
  void set_seed(const int& seed){ ran_->SetSeed(seed);}

  // Function to get external or internal true parameter values from index (for toys)
  double get_true_params(const unsigned int& i, const bool& external) {
    if(external)
      return x_vals_(i);
    else
      return (U_*x_vals_)(i);
  }

  // Get the sum of pT scale biases A,e or M obtained from all the previous iterations
  double get_A_prevfit(const unsigned int& i) {
    return A_vals_prevfit_(i);
  }
  double get_e_prevfit(const unsigned int& i) {
    return e_vals_prevfit_(i);
  }
  double get_M_prevfit(const unsigned int& i) {
    return M_vals_prevfit_(i);
  }

  unsigned int get_n_params(){ return n_pars_;}
  unsigned int get_n_data(){ return n_data_;}
  unsigned int get_n_dof(){ return n_dof_;}

  double get_first_pt_edge(){ return pt_edges_.at(0) ;}
  double get_last_pt_edge(){ return pt_edges_.at(n_pt_bins_) ;} 

  // Get the transformation of external to internal parameters
  double get_U(const unsigned int& i, const unsigned int& j) {
    return U_(i,j);
  }
  
  virtual double Up() const {return errorDef_;}
  virtual void SetErrorDef(double def) {errorDef_ = def;}

  // Function to be minimised from mass scale bias ^2 values, errors and pT scale biases parameters AeM
  virtual double operator()(const vector<double>&) const;
  // Function for analytical gradient of the function to be minimised
  virtual vector<double> Gradient(const vector<double>& ) const;
  virtual bool CheckGradient() const {return true;} 

private:

  vector<double> scales2_;
  vector<double> scales2Err_;
  vector<int> masks_;
  vector<float> pt_edges_;
  vector<double> k_edges_;
  vector<double> kmean_vals_;
  VectorXd A_vals_;
  VectorXd e_vals_;
  VectorXd M_vals_;
  VectorXd x_vals_;
  VectorXd A_vals_prevfit_;
  VectorXd e_vals_prevfit_;
  VectorXd M_vals_prevfit_;
  double kmean_val_;
  vector<float> eta_edges_;
  unsigned int n_pt_bins_;
  unsigned int n_eta_bins_;
  unsigned int n_data_;
  unsigned int n_pars_;
  unsigned int n_dof_;
  int debug_;
  int seed_;
  int bias_;
  double errorDef_;
  MatrixXd U_;
  TRandom3* ran_;
};

// In toy mode, function to generate mass scale bias ^2 values and errors from given pT scale bias AeM via Gaussian sampling
void TheoryFcn::generate_data() {
  double chi2_start = 0.;
  unsigned int ibin = 0;
  for(unsigned int ieta_p = 0; ieta_p<n_eta_bins_; ieta_p++) {
    for(unsigned int ipt_p = 0; ipt_p<n_pt_bins_; ipt_p++) {
      double k_p = kmean_vals_[ipt_p];
      for(unsigned int ieta_m = 0; ieta_m<n_eta_bins_; ieta_m++) {
	      for(unsigned int ipt_m = 0; ipt_m<n_pt_bins_; ipt_m++) {
	        double k_m = kmean_vals_[ipt_m];

	        // Draw error on mass scale bias ^2 centered around ierr2_nom as function of eta
	        double ierr2_nom = 0.0001*(1+double(ieta_p)/n_eta_bins_)*(1+double(ieta_m)/n_eta_bins_);
	        //*(2-0.1*double(ipt_p)/n_pt_bins_)*(2-0.1*double(ipt_m)/n_pt_bins_);
	        double ierr2 = ran_->Gaus(ierr2_nom,  ierr2_nom*0.1);
	        while(ierr2<=0.) 
	          ierr2 = ran_->Gaus(ierr2_nom,  ierr2_nom*0.1);
	        
          // Draw scale^2 centered around iscale2_bias with width ierr2
          // Note: the model is correct as we read - _nom histos  
	        double iscale2_bias =
	          (1.0 + A_vals_(ieta_p) + e_vals_(ieta_p)*k_p - M_vals_(ieta_p)/k_p)*
	          (1.0 + A_vals_(ieta_m) + e_vals_(ieta_m)*k_m + M_vals_(ieta_m)/k_m);
	        double iscale2 = ran_->Gaus(iscale2_bias, ierr2);

	        //if(ibin<3) cout << iscale2 << endl;
	        scales2_[ibin]    = iscale2 ;
	        scales2Err_[ibin] =  ierr2 ;
	        double dchi2 = (scales2_[ibin]-1.0)/scales2Err_[ibin];
	        //cout << dchi2*dchi2 << endl;
	        chi2_start += dchi2*dchi2 ;
	        ibin++;
	      }
      }
    }
  }
  cout << "[Toy mode] Initial chi2 = " << chi2_start << " / " << n_data_ << " ndof has prob " << TMath::Prob(chi2_start, n_data_ ) <<  endl;
  return;
}

// Define function to be minimised from mass scale bias ^2 values and errors and pT scale biases parameters AeM -> will obtain AeM
double TheoryFcn::operator()(const vector<double>& par) const {

  double val = 0.0;
  const unsigned int npars = par.size();

  // Keep track of 4D bins (the data points)
  unsigned int ibin = 0;
  // +ve muon term
  for(unsigned int ieta_p = 0; ieta_p < n_eta_bins_; ieta_p++) { // AeM are eta dependent
    double A_p = par[ieta_p];
    double e_p = par[ieta_p+n_eta_bins_];
    double M_p = par[ieta_p+2*n_eta_bins_]; 
    for(unsigned int ipt_p = 0; ipt_p < n_pt_bins_; ipt_p++) {   
      double k_p = kmean_vals_[ipt_p];
      double p_term = (1.0 + A_p - e_p*(k_p-kmean_val_)/kmean_val_ + M_p/k_p*kmean_val_ );
      // -ve muon term
      for(unsigned int ieta_m = 0; ieta_m < n_eta_bins_; ieta_m++) { // AeM are eta dependent
	      double A_m = par[ieta_m];
	      double e_m = par[ieta_m+n_eta_bins_];
	      double M_m = par[ieta_m+2*n_eta_bins_];
	      for(unsigned int ipt_m = 0; ipt_m < n_pt_bins_; ipt_m++) {	  
	        double k_m = kmean_vals_[ipt_m];
	        double m_term = (1.0 + A_m - e_m*(k_m-kmean_val_)/kmean_val_ - M_m/k_m*kmean_val_);
          // Function to be minimized is ( chi2/ndf - 1 )
	        double ival = (scales2_[ibin] - p_term*m_term)/scales2Err_[ibin];
	        double ival2 = ival*ival;
	        if(masks_[ibin]) // check if accepting or ignoring the 4D bin
	          val += ival2;
	        ibin++;
	      }
      }
    }
  }
  // Function to be minimized is ( chi2/ndf - 1 )
  val /= n_dof_;
  val -= 1.0;
  
  return val;
}

vector<double> TheoryFcn::Gradient(const vector<double> &par ) const {

  //cout << "Using gradient" << endl; 
  vector<double> grad(par.size(), 0.0);

  // Loop over all AeM parameters
  for(unsigned int ipar = 0; ipar < par.size(); ipar++) {
    unsigned int ieta     = ipar % n_eta_bins_;
    unsigned int par_type = ipar / n_eta_bins_; // A, e or M    
    //cout << "ipar " << ipar << ": " << ieta << ", " << par_type << endl;

    // Gradient of function to be minimized wrt to a given parameter in a given 4D bin
    double grad_i = 0.0;    
    // Keep track of 4D bins (the data points)
    unsigned int ibin = 0;
    
    // For the given parameter, check all the 4D bins (data points)
    // +ve muon term
    for(unsigned int ieta_p = 0; ieta_p < n_eta_bins_; ieta_p++) { // AeM are eta dependent
      double A_p = par[ieta_p];
      double e_p = par[ieta_p+n_eta_bins_];
      double M_p = par[ieta_p+2*n_eta_bins_]; 
      for(unsigned int ipt_p = 0; ipt_p < n_pt_bins_; ipt_p++) { 
	      double k_p = kmean_vals_[ipt_p];
	      double p_term = 0.;
	      if(ieta_p != ieta)
          p_term = (1.0 + A_p - e_p*(k_p-kmean_val_)/kmean_val_ + M_p/k_p*kmean_val_);
	      else {
	        if(par_type==0)       p_term = 1.0;
	        else if( par_type==1) p_term = -(k_p-kmean_val_)/kmean_val_;
	        else                  p_term = 1./k_p*kmean_val_;
	      }
        // -ve muon term
	      for(unsigned int ieta_m = 0; ieta_m < n_eta_bins_; ieta_m++) { // AeM are eta dependent
	        double A_m = par[ieta_m];
	        double e_m = par[ieta_m+n_eta_bins_];
	        double M_m = par[ieta_m+2*n_eta_bins_];
	        for(unsigned int ipt_m = 0; ipt_m < n_pt_bins_; ipt_m++) {	  
	          double k_m = kmean_vals_[ipt_m];
	          double m_term = 0.;
	          if(ieta_m != ieta)
              m_term = (1.0 + A_m - e_m*(k_m-kmean_val_)/kmean_val_ - M_m/k_m*kmean_val_);
	          else {
	            if(par_type==0)       m_term = 1.0;
	            else if( par_type==1) m_term = -(k_m-kmean_val_)/kmean_val_;
	            else                  m_term = -1./k_m*kmean_val_;
	          }

	          double ival = -2*(scales2_[ibin] - (1.0 + A_p - e_p*(k_p-kmean_val_)/kmean_val_ + M_p/k_p*kmean_val_)*
			        (1.0 + A_m - e_m*(k_m-kmean_val_)/kmean_val_ - M_m/k_m*kmean_val_))
	            /scales2Err_[ibin]/scales2Err_[ibin];

	          double term = 0.0;
	          if(ieta_p==ieta || ieta_m==ieta) { // if this 4D bin contains the relevant parameter 
	            if(ieta_p!=ieta_m) term = p_term*m_term; // if it contains it only in the +ve or -ve term
	            else { // if it contains it in both the +ve and the -ve term
		            if(par_type==0) // A
		              term = 1.0*(1.0 + A_m - e_m*(k_m-kmean_val_)/kmean_val_ - M_m/k_m*kmean_val_) +
		                (1.0 + A_p - e_p*(k_p-kmean_val_)/kmean_val_ + M_p/k_p*kmean_val_)*1.0;
		            else if(par_type==1) // e 
		              term = -(k_p-kmean_val_)/kmean_val_ * (1.0 + A_m - e_m*(k_m-kmean_val_)/kmean_val_ - M_m/k_m*kmean_val_) -
		                (1.0 + A_p - e_p*(k_p-kmean_val_)/kmean_val_ + M_p/k_p*kmean_val_) * (k_m-kmean_val_)/kmean_val_;
		            else // M
		              term = 1.0/k_p*kmean_val_ * (1.0 + A_m - e_m*(k_m-kmean_val_)/kmean_val_ - M_m/k_m*kmean_val_) -
		                (1.0 + A_p - e_p*(k_p-kmean_val_)/kmean_val_ + M_p/k_p*kmean_val_) * 1.0/k_m*kmean_val_;
	            }
	          }
	          //cout << "ival " << ival << "," << term << endl; 
	          double ig = ival*term;
            // Function to be minimized is ( chi2/ndf - 1 )
	          ig /= n_dof_;
	          //cout << "ibin " << ibin << " += " << ig << endl;
	          if(masks_[ibin]) grad_i += ig; // check if accepting or ignoring the 4D bin
	          ibin++;
	        }
	      }
      }
    }
    //cout << "\t" << ipar << ": " << grad_i << endl;
    grad[ipar] = grad_i; 
  }

  return grad; 
}


  
int main(int argc, char* argv[]) {

  TStopwatch sw;
  sw.Start();

  //ROOT::EnableImplicitMT();

  variables_map vm;
  try {
    options_description desc{"Options"};
    desc.add_options()
	    ("help,h", "Help screen")
	    ("ntoys",  value<long>()->default_value(1), "number of toys, should be 1 to use data")
	    ("tag",    value<std::string>()->default_value("closure"), "tag of input data")
	    ("run",    value<std::string>()->default_value("closure"), "run of input data")
	    ("bias",   value<int>()->default_value(0), "bias [-1 for data, >0 for toys: 1 for uniform random bias, 2 for eta dependent bias]")
	    ("infile", value<std::string>()->default_value("massscales"), "type of input data")
	    ("seed",   value<int>()->default_value(4357), "seed for random toys with different AeM bias");

    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);
    if (vm.count("help")) {
	    std::cout << desc << '\n';
	    return 0;
    }
    if (vm.count("ntoys"))    std::cout << "Number of toys: " << vm["ntoys"].as<long>() << '\n';
    if (vm.count("tag"))      std::cout << "Tag: " << vm["tag"].as<std::string>() << '\n';
    if (vm.count("run"))      std::cout << "Run: " << vm["run"].as<std::string>() << '\n';
  }
  catch (const error &ex) {
    std::cerr << ex.what() << '\n';
  }

  long ntoys         = vm["ntoys"].as<long>();
  std::string tag    = vm["tag"].as<std::string>();
  std::string infile = vm["infile"].as<std::string>();
  std::string run    = vm["run"].as<std::string>();
  int bias           = vm["bias"].as<int>();
  int seed           = vm["seed"].as<int>();

  TFile* fout = TFile::Open(("./massfit_"+tag+"_"+run+".root").c_str(), "RECREATE");
  
  // Define TTree with fit qualities
  TTree* tree = new TTree("tree", "tree");
  double edm, fmin, prob;
  int isvalid, hasAccurateCovar, hasPosDefCovar;
  tree->Branch("edm", &edm, "edm/D");
  tree->Branch("fmin", &fmin, "fmin/D");
  tree->Branch("prob", &prob, "prob/D");
  tree->Branch("isvalid", &isvalid, "isvalid/I");
  tree->Branch("hasAccurateCovar", &hasAccurateCovar, "hasAccurateCovar/I");
  tree->Branch("hasPosDefCovar", &hasPosDefCovar, "hasPosDefCovar/I");

  // Initialize function to be minimized ( chi2/ndf - 1 )
  int debug = 0;
  string infname = infile+"_"+tag+"_"+run+".root";
  TheoryFcn* fFCN = new TheoryFcn(debug, seed, bias, infname);  
  fFCN->SetErrorDef(1.0 / fFCN->get_n_dof());
  unsigned int n_parameters = fFCN->get_n_params();
  // Get the transformation of external to internal parameters
  MatrixXd U(n_parameters,n_parameters);
  for (int i=0; i<n_parameters; i++) {
    for (int j=0; j<n_parameters; j++) {
      U(i,j) = fFCN->get_U(i,j);
    }
  }
  MatrixXd Uinv = U.inverse(); // transformation of internal to external parameters
  
  // Single vectors to store all pT scale bias parameters AeM
  vector<double> tparIn0(n_parameters);  // internal, true (toys)
  vector<double> tparIn(n_parameters);   // internal, fitted
  vector<double> tparInErr(n_parameters);
  vector<double> tparOut0(n_parameters); // external, true (toys)
  vector<double> tparOut(n_parameters);  // external, fitted
  vector<double> tparOutErr(n_parameters);

  // Fit parameters TTree branches 
  for (int i=0; i<n_parameters/3; i++){
    tree->Branch(Form("A%d",i),       &tparOut[i],    Form("A%d/D",i));
    tree->Branch(Form("A%d_true",i),  &tparOut0[i],   Form("A%d_true/D",i));
    tree->Branch(Form("A%d_err",i),   &tparOutErr[i], Form("A%d_err/D",i));
    tree->Branch(Form("A%d_in",i),    &tparIn[i],    Form("A%d_in/D",i));
    tree->Branch(Form("A%d_intrue",i),&tparIn0[i],    Form("A%d_intrue/D",i));
    tree->Branch(Form("A%d_inerr",i), &tparInErr[i], Form("A%d_inerr/D",i));
  }
  for (int i=0; i<n_parameters/3; i++){
    tree->Branch(Form("e%d",i),        &tparOut[i+n_parameters/3],    Form("e%d/D",i));
    tree->Branch(Form("e%d_true",i),   &tparOut0[i+n_parameters/3],   Form("e%d_true/D",i));
    tree->Branch(Form("e%d_err",i),    &tparOutErr[i+n_parameters/3], Form("e%d_err/D",i));
    tree->Branch(Form("e%d_in",i),     &tparIn[i+n_parameters/3],     Form("e%d_in/D",i));
    tree->Branch(Form("e%d_intrue",i), &tparIn0[i+n_parameters/3],    Form("e%d_intrue/D",i));
    tree->Branch(Form("e%d_inerr",i),  &tparInErr[i+n_parameters/3],  Form("e%d_inerr/D",i));
  }
  for (int i=0; i<n_parameters/3; i++){
    tree->Branch(Form("M%d",i),        &tparOut[i+2*n_parameters/3],    Form("M%d/D",i));
    tree->Branch(Form("M%d_true",i),   &tparOut0[i+2*n_parameters/3],   Form("M%d_true/D",i));
    tree->Branch(Form("M%d_err",i),    &tparOutErr[i+2*n_parameters/3], Form("M%d_err/D",i));
    tree->Branch(Form("M%d_in",i),     &tparIn[i+2*n_parameters/3],     Form("M%d_in/D",i));
    tree->Branch(Form("M%d_intrue",i), &tparIn0[i+2*n_parameters/3],    Form("M%d_in/D",i));
    tree->Branch(Form("M%d_inerr",i),  &tparInErr[i+2*n_parameters/3],  Form("M%d_inerr/D",i));
  }

  // Nom histograms needed for toy mode
  TH1D* h_A_vals_nom  = new TH1D("h_A_vals_nom", "A nominal", n_parameters/3, 0, n_parameters/3);
  TH1D* h_e_vals_nom  = new TH1D("h_e_vals_nom", "e nominal", n_parameters/3, 0, n_parameters/3);
  TH1D* h_M_vals_nom  = new TH1D("h_M_vals_nom", "M nominal", n_parameters/3, 0, n_parameters/3);
  TH1D* h_Ain_vals_nom  = new TH1D("h_Ain_vals_nom", "(A-?e#bar{k})", n_parameters/3, 0, n_parameters/3);
  TH1D* h_ein_vals_nom  = new TH1D("h_ein_vals_nom", "e/#bar{k} nominal", n_parameters/3, 0, n_parameters/3);
  TH1D* h_Min_vals_nom  = new TH1D("h_Min_vals_nom", "M#bar{k} nominal", n_parameters/3, 0, n_parameters/3);

  // Save external and internal fitted parameters representing pT scale biases A,e or M obtained in the current iteration   
  TH1D* h_A_vals_fit  = new TH1D("h_A_vals_fit", "#hat{A}", n_parameters/3, 0, n_parameters/3);
  TH1D* h_e_vals_fit  = new TH1D("h_e_vals_fit", "#hat{e}", n_parameters/3, 0, n_parameters/3);
  TH1D* h_M_vals_fit  = new TH1D("h_M_vals_fit", "#hat{M}", n_parameters/3, 0, n_parameters/3);
  TH1D* h_Ain_vals_fit  = new TH1D("h_Ain_vals_fit", "(#hat{A}-#hat{e}#bar{k})", n_parameters/3, 0, n_parameters/3);
  TH1D* h_ein_vals_fit  = new TH1D("h_ein_vals_fit", "#hat{e}#bar{k}", n_parameters/3, 0, n_parameters/3);
  TH1D* h_Min_vals_fit  = new TH1D("h_Min_vals_fit", "#hat{M}/#bar{k}", n_parameters/3, 0, n_parameters/3);
  
  // We will overwrite the _prevfit histograms to contain the sum of pT scale biases A,e or M obtained in all the previous iterations + the ones obtained in the current iteration
  TH1D* h_A_vals_prevfit  = new TH1D("h_A_vals_prevfit", "#hat{A}", n_parameters/3, 0, n_parameters/3);
  TH1D* h_e_vals_prevfit  = new TH1D("h_e_vals_prevfit", "#hat{e}", n_parameters/3, 0, n_parameters/3);
  TH1D* h_M_vals_prevfit  = new TH1D("h_M_vals_prevfit", "#hat{M}", n_parameters/3, 0, n_parameters/3);

  // Toy mode: mass scale biases from input pT scale bias AeM
  TH2D* h_scales_nom_plus   = new TH2D("h_scales_nom_plus", "scales nominal plus; #eta bin", n_parameters/3, 0, n_parameters/3,
				       50, fFCN->get_first_pt_edge(), fFCN->get_last_pt_edge() );
  // Data mode: scales from the sum of pT scale biases A,e or M obtained in all the previous iterations + the ones obtained in the current iteration
  TH2D* h_scales_fit_plus   = new TH2D("h_scales_fit_plus", "scales plus; #eta bin", n_parameters/3, 0, n_parameters/3,
				       50, fFCN->get_first_pt_edge(), fFCN->get_last_pt_edge() );
  TH2D* h_scales_nom_minus  = new TH2D("h_scales_nom_minus", "scales nominal minus; #eta bin", n_parameters/3, 0, n_parameters/3,
				       50, fFCN->get_first_pt_edge(), fFCN->get_last_pt_edge() );
  TH2D* h_scales_fit_minus  = new TH2D("h_scales_fit_minus", "scales minus; #eta bin", n_parameters/3, 0, n_parameters/3,
				       50, fFCN->get_first_pt_edge(), fFCN->get_last_pt_edge() );

  unsigned int maxfcn(numeric_limits<unsigned int>::max());
  double tolerance(0.001);
  int verbosity = int(ntoys<2); 
  ROOT::Minuit2::MnPrint::SetGlobalLevel(verbosity);
  
  if(bias<0) assert( ntoys == 1); // use ntoys==1 for data
  for(unsigned int itoy=0; itoy<ntoys; itoy++) {

    if(itoy%10==0) cout << "Toy " << itoy << " / " << ntoys << endl;

    //fFCN->set_seed(seed);
    if(bias>=0) fFCN->generate_data();
    
    // Define minimization parameters
    MnUserParameters upar;
    double start=0.0, par_error=0.01;
    for (int i=0; i<n_parameters/3; i++) upar.Add(Form("A%d",i), start, par_error);
    for (int i=0; i<n_parameters/3; i++) upar.Add(Form("e%d",i), start, par_error);
    for (int i=0; i<n_parameters/3; i++) upar.Add(Form("M%d",i), start, par_error);      

    // Minimize
    MnMigrad migrad(*fFCN, upar, 1);    
    if(itoy<1) cout << "\tMigrad..." << endl;
    FunctionMinimum min = migrad(maxfcn, tolerance);

    // Fit properties
    edm = double(min.Edm());
    fmin = double(min.Fval());
    prob = TMath::Prob((min.Fval()+1)*fFCN->get_n_dof(), fFCN->get_n_dof() );
    isvalid = int(min.IsValid());
    hasAccurateCovar = int(min.HasAccurateCovar());
    hasPosDefCovar = int(min.HasPosDefCovar());
    
    if(itoy<1) cout << "\tHesse..." << endl;
    MnHesse hesse(1);
    hesse(*fFCN, min);

    if(itoy<1) cout << "\t => final chi2/ndf: " << min.Fval()+1 << " (prob: " << TMath::Prob((min.Fval()+1)*fFCN->get_n_dof(), fFCN->get_n_dof() ) << ")" <<  endl;;
    
    // Internal covariance matrix
    MatrixXd Vin(n_parameters,n_parameters);    
    for(unsigned int i = 0 ; i<n_parameters; i++) {    
      for(unsigned int j = 0 ; j<n_parameters; j++) {
	      Vin(i,j) = i>j ?
	        min.UserState().Covariance().Data()[j+ i*(i+1)/2] :
	        min.UserState().Covariance().Data()[i+ j*(j+1)/2];
      }
    }
    // External covariance matrix
    MatrixXd Vout = Uinv*Vin*Uinv.transpose();

    // Internal fitted parameters
    VectorXd xin(n_parameters);
    VectorXd xinErr(n_parameters);     
    for(unsigned int i = 0 ; i<n_parameters; i++) {
      xin(i)    = min.UserState().Value(i) ;
      xinErr(i) = min.UserState().Error(i) ;
    }
    // External fitted parameters
    VectorXd x = Uinv*xin;
    VectorXd xErr(n_parameters);
    for(unsigned int i = 0 ; i<n_parameters; i++) {
      xErr(i) = TMath::Sqrt(Vout(i,i));
    }    

    // Save scale histograms for first toy / data
    // Note this is incorrect for toys for now due to the AeM signs and nom histograms
    if(itoy<1) {
      for(unsigned int ib = 0 ; ib<n_parameters/3; ib++) {
        Eigen::Vector3d xi;
        xi <<
	        x(ib) + fFCN->get_A_prevfit(ib),
	        x(ib + n_parameters/3) + fFCN->get_e_prevfit(ib),
	        x(ib + 2*n_parameters/3) + fFCN->get_M_prevfit(ib); 
        Eigen::Vector3d xnomi;
        xnomi <<
	        fFCN->get_true_params(ib, true),
	        fFCN->get_true_params(ib + n_parameters/3, true),
	        fFCN->get_true_params(ib + 2*n_parameters/3, true);
        for(unsigned int jb = 0 ; jb<h_scales_nom_plus->GetYaxis()->GetNbins(); jb++) {
	        double kj = 1./h_scales_nom_plus->GetYaxis()->GetBinCenter(jb+1);
	        Eigen::Vector3d ajp;
	        Eigen::Vector3d ajm;
	        ajp << 1.0, -1.*kj, +1./kj;
	        ajm << 1.0, -1.*kj, -1./kj;
	        Eigen::Matrix3d Vj;
	        Vj <<
	          Vout(ib, ib),                  Vout(ib, ib+n_parameters/3),                  Vout(ib, ib+2*n_parameters/3),
	          Vout(ib+n_parameters/3, ib),   Vout(ib+n_parameters/3, ib+n_parameters/3) ,  Vout(ib+n_parameters/3, ib+2*n_parameters/3),
	          Vout(ib+2*n_parameters/3, ib), Vout(ib+2*n_parameters/3, ib+n_parameters/3), Vout(ib+2*n_parameters/3, ib+2*n_parameters/3);
	        MatrixXd scalejp    = ajp.transpose()*xi;
	        MatrixXd scalenomjp = ajp.transpose()*xnomi;
	        MatrixXd Vscalejp   = ajp.transpose()*Vj*ajp;
	        MatrixXd scalejm    = ajm.transpose()*xi;
	        MatrixXd scalenomjm = ajm.transpose()*xnomi;
	        MatrixXd Vscalejm   = ajm.transpose()*Vj*ajm;
	        h_scales_nom_plus->SetBinContent(ib+1, jb+1, 1.0 + scalenomjp(0,0) );
	        h_scales_fit_plus->SetBinContent(ib+1, jb+1, 1.0 + scalejp(0,0) );
	        h_scales_fit_plus->SetBinError(ib+1, jb+1, TMath::Sqrt(Vscalejp(0,0)) );
	        h_scales_nom_minus->SetBinContent(ib+1, jb+1, 1.0 + scalenomjm(0,0) );
	        h_scales_fit_minus->SetBinContent(ib+1, jb+1, 1.0 + scalejm(0,0) );
	        h_scales_fit_minus->SetBinError(ib+1, jb+1, TMath::Sqrt(Vscalejm(0,0)) );
        }      
      }
    }
    
    for(unsigned int i = 0 ; i<n_parameters; i++) {
      // Save parameter values
      tparIn[i]     = xin(i);
      tparIn0[i]    = fFCN->get_true_params(i, false) ;
      tparInErr[i]  = xinErr(i);
      tparOut[i]    = x(i);
      tparOut0[i]   = fFCN->get_true_params(i, true) ;
      tparOutErr[i] = xErr(i);
      //cout << "Param " << i << ": " << x(i) << " +/- " << xErr(i) << ". True value is " << fFCN->get_true_params(i, true) << endl;

      // Save AeM histograms for first toy / data
      if(itoy<1) {
        int ip = i%(n_parameters/3);
        if(i<n_parameters/3) {
	        h_A_vals_fit->SetBinContent(ip+1, x(i));
	        h_A_vals_fit->SetBinError(ip+1, xErr(i));
	        h_A_vals_nom->SetBinContent(ip+1, fFCN->get_true_params(i, true)); // toys
	        h_Ain_vals_fit->SetBinContent(ip+1, xin(i));
	        h_Ain_vals_fit->SetBinError(ip+1, xinErr(i));
	        h_Ain_vals_nom->SetBinContent(ip+1, fFCN->get_true_params(i, false)); // toys
	        h_A_vals_prevfit->SetBinContent(ip+1, fFCN->get_A_prevfit(ip) + x(i));
        }
        else if(i>=n_parameters/3 && i<2*n_parameters/3) {
	        h_e_vals_fit->SetBinContent(ip+1, x(i));
	        h_e_vals_fit->SetBinError(ip+1, xErr(i));
	        h_e_vals_nom->SetBinContent(ip+1, fFCN->get_true_params(i, true)); // toys
	        h_ein_vals_fit->SetBinContent(ip+1, xin(i));
	        h_ein_vals_fit->SetBinError(ip+1, xinErr(i));
	        h_ein_vals_nom->SetBinContent(ip+1, fFCN->get_true_params(i, false)); // toys
	        h_e_vals_prevfit->SetBinContent(ip+1, fFCN->get_e_prevfit(ip) + x(i));
        }
        else {
	        h_M_vals_fit->SetBinContent(ip+1, x(i));
	        h_M_vals_fit->SetBinError(ip+1, xErr(i));
	        h_M_vals_nom->SetBinContent(ip+1, fFCN->get_true_params(i, true)); // toys
	        h_Min_vals_fit->SetBinContent(ip+1, xin(i));
      	  h_Min_vals_fit->SetBinError(ip+1, xinErr(i));
	        h_Min_vals_nom->SetBinContent(ip+1, fFCN->get_true_params(i, false)); // toys
	        h_M_vals_prevfit->SetBinContent(ip+1, fFCN->get_M_prevfit(ip) + x(i));
        }  
      } 
    }

    tree->Fill();

    // Save covariance and correlation matrices for first toy / data
    if(itoy<1) {   
      TH2D* hcov = new TH2D(Form("hcov_%d", itoy), "", n_parameters, 0, n_parameters, n_parameters, 0, n_parameters);
      TH2D* hcor = new TH2D(Form("hcor_%d", itoy), "", n_parameters, 0, n_parameters, n_parameters, 0, n_parameters);  
      TH2D* hcovin = new TH2D(Form("hcovin_%d", itoy), "", n_parameters, 0, n_parameters, n_parameters, 0, n_parameters);
      TH2D* hcorin = new TH2D(Form("hcorin_%d", itoy), "", n_parameters, 0, n_parameters, n_parameters, 0, n_parameters);  

      for(unsigned int i = 0 ; i<n_parameters; i++) {    
        hcov->GetXaxis()->SetBinLabel(i+1, TString(upar.GetName(i).c_str()) );
        hcor->GetXaxis()->SetBinLabel(i+1, TString(upar.GetName(i).c_str()) );
        hcovin->GetXaxis()->SetBinLabel(i+1, TString(upar.GetName(i).c_str()) );
        hcorin->GetXaxis()->SetBinLabel(i+1, TString(upar.GetName(i).c_str()) );
        for(unsigned int j = 0 ; j<n_parameters; j++) {
	        double covin_ij = Vin(i,j);
	        double corin_ij = Vin(i,j)/TMath::Sqrt(Vin(i,i)*Vin(j,j)); 
	        double cov_ij = Vout(i,j);
	        double cor_ij = Vout(i,j)/TMath::Sqrt(Vout(i,i)*Vout(j,j)); 
	        hcov->GetYaxis()->SetBinLabel(j+1, TString(upar.GetName(j).c_str()) );
	        hcor->GetYaxis()->SetBinLabel(j+1, TString(upar.GetName(j).c_str()) );
      	  hcovin->GetYaxis()->SetBinLabel(j+1, TString(upar.GetName(j).c_str()) );
	        hcorin->GetYaxis()->SetBinLabel(j+1, TString(upar.GetName(j).c_str()) );
	        hcovin->SetBinContent(i+1, j+1, covin_ij);
	        hcorin->SetBinContent(i+1, j+1, corin_ij);
	        hcov->SetBinContent(i+1, j+1, cov_ij);
	        hcor->SetBinContent(i+1, j+1, cor_ij);
        }
      }

      hcor->SetMinimum(-1.0);
      hcor->SetMaximum(+1.0);
      hcorin->SetMinimum(-1.0);
      hcorin->SetMaximum(+1.0);
    
      fout->cd();
    
      hcor->Write();
      hcov->Write();
      hcorin->Write();
      hcovin->Write();
    }

    if(verbosity) {
      cout << "Data points: " << fFCN->get_n_data() << endl;
      cout << "Number of parameters: " << fFCN->get_n_params() << endl;
      cout << "chi2/ndf: " << min.Fval()+1 << " (prob: " << TMath::Prob((min.Fval()+1)*fFCN->get_n_dof(), fFCN->get_n_dof() ) << ")" <<  endl;;
      cout << "min is valid: " << min.IsValid() << std::endl;
      cout << "HesseFailed: " << min.HesseFailed() << std::endl;
      cout << "HasCovariance: " << min.HasCovariance() << std::endl;
      cout << "HasValidCovariance: " << min.HasValidCovariance() << std::endl;
      cout << "HasValidParameters: " << min.HasValidParameters() << std::endl;
      cout << "IsAboveMaxEdm: " << min.IsAboveMaxEdm() << std::endl;
      cout << "HasReachedCallLimit: " << min.HasReachedCallLimit() << std::endl;
      cout << "HasAccurateCovar: " << min.HasAccurateCovar() << std::endl;
      cout << "HasPosDefCovar : " << min.HasPosDefCovar() << std::endl;
      cout << "HasMadePosDefCovar : " << min.HasMadePosDefCovar() << std::endl;
    }
  }

  fout->cd();
  tree->Write();

  // For toys, pull distributions of the fitted vs true AeM parameters
  TH1D* hpulls = new TH1D("hpulls", "", n_parameters, 0, n_parameters);
  TH1D* hsigma = new TH1D("hsigma", "", n_parameters, 0, n_parameters);
  for (int i=0; i<n_parameters; i++) {
    TH1D* h = new TH1D(Form("h%d", i), "", 100,-3,3);
    int ip = i%(n_parameters/3);
    if(i<n_parameters/3) {
      tree->Draw(Form("(A%d - A%d_true)/A%d_err>>h%d", ip, ip, ip, i), "", "");
      hpulls->GetXaxis()->SetBinLabel(i+1, Form("A%d", ip));
      h_A_vals_fit->GetXaxis()->SetBinLabel(ip+1, Form("A%d", ip));
      h_Ain_vals_fit->GetXaxis()->SetBinLabel(ip+1, Form("Ain%d", ip));
      h_A_vals_nom->GetXaxis()->SetBinLabel(ip+1, Form("A%d", ip));
      h_Ain_vals_nom->GetXaxis()->SetBinLabel(ip+1, Form("Ain%d", ip));
    }
    else if(i>=n_parameters/3 && i<2*n_parameters/3) {
      tree->Draw(Form("(e%d - e%d_true)/e%d_err>>h%d", ip, ip, ip, i), "", "");
      hpulls->GetXaxis()->SetBinLabel(i+1, Form("e%d", ip));
      h_e_vals_fit->GetXaxis()->SetBinLabel(ip+1, Form("e%d", ip));
      h_ein_vals_fit->GetXaxis()->SetBinLabel(ip+1, Form("ein%d", ip));
      h_e_vals_nom->GetXaxis()->SetBinLabel(ip+1, Form("e%d", ip));
      h_ein_vals_nom->GetXaxis()->SetBinLabel(ip+1, Form("ein%d", ip));
    }
    else {
      tree->Draw(Form("(M%d - M%d_true)/M%d_err>>h%d", ip, ip, ip, i), "", "");
      hpulls->GetXaxis()->SetBinLabel(i+1, Form("M%d", ip));
      h_M_vals_fit->GetXaxis()->SetBinLabel(ip+1, Form("M%d", ip));
      h_Min_vals_fit->GetXaxis()->SetBinLabel(ip+1, Form("Min%d", ip));
      h_M_vals_nom->GetXaxis()->SetBinLabel(ip+1, Form("M%d", ip));
      h_Min_vals_nom->GetXaxis()->SetBinLabel(ip+1, Form("Min%d", ip));
    }
    //cout << i << "-->" << h->GetMean() << endl;
    float pull_i = h->GetMean();
    float pull_i_err = 0.;
    float sigma_i = 0.;
    float sigma_i_err = 0.;
    if(h->GetEntries()>10) {
      h->Fit("gaus", "Q");
      TF1* gaus = (TF1*)h->GetFunction("gaus");
      if(gaus==0) {
	      cout << "no func" << endl;
	      continue;
      }
      pull_i = gaus->GetParameter(1);
      pull_i_err = gaus->GetParError(1);
      sigma_i = gaus->GetParameter(2);
      sigma_i_err = gaus->GetParError(2);
    }
    hpulls->SetBinContent(i+1, pull_i);
    hpulls->SetBinError(i+1, pull_i_err);
    hsigma->SetBinContent(i+1, sigma_i);
    hsigma->SetBinError(i+1, sigma_i_err);
    delete h;
  }
  hpulls->Write();
  hsigma->Write();

  h_A_vals_fit->Write();
  h_e_vals_fit->Write();
  h_M_vals_fit->Write();
  h_A_vals_prevfit->Write();
  h_e_vals_prevfit->Write();
  h_M_vals_prevfit->Write();
  h_Ain_vals_fit->Write();
  h_ein_vals_fit->Write();
  h_Min_vals_fit->Write();
  h_A_vals_nom->Write();
  h_e_vals_nom->Write();
  h_M_vals_nom->Write();
  h_Ain_vals_nom->Write();
  h_ein_vals_nom->Write();
  h_Min_vals_nom->Write();
  h_scales_nom_plus->Write();
  h_scales_fit_plus->Write();
  h_scales_nom_minus->Write();
  h_scales_fit_minus->Write();
  
  sw.Stop();
  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;
  fout->Close(); 

  return 0;
}
