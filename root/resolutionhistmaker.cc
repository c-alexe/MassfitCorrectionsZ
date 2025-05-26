#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
#include <string>

#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string.hpp>

using namespace ROOT;
using namespace ROOT::VecOps;

float deltaPhi(float phi1, float phi2)
{                                                        
  float result = phi1 - phi2;
  while (result > float(M_PI)) result -= float(2*M_PI);
  while (result <= -float(M_PI)) result += float(2*M_PI);
  return result;
}

float deltaR2(float eta1, float phi1, float eta2, float phi2)
{
  float deta = std::abs(eta1-eta2);
  float dphi = deltaPhi(phi1,phi2);
  return deta*deta + dphi*dphi;
}

float deltaR(float eta1, float phi1, float eta2, float phi2)
{
  return std::sqrt(deltaR2(eta1,phi1,eta2,phi2));
}

RVec<Int_t> hasGenMatchIdx(RVec<Float_t> &GenPart_eta,
			 RVec<Float_t> &GenPart_phi, RVec<Float_t> &Cand_eta, 
			 RVec<Float_t> &Cand_phi, double dR_angle=0.1)
{
  RVec<Int_t> isGenMatched;  
  for(unsigned int iGen=0; iGen<GenPart_eta.size(); iGen++){
    float mcmatch_tmp_dr = 999.;
    int idx=-1;
    for(int iCand=0;iCand<Cand_eta.size();iCand++){
      float tmpDR = deltaR(Cand_eta[iCand], Cand_phi[iCand], GenPart_eta[iGen], GenPart_phi[iGen]);
      if (tmpDR < mcmatch_tmp_dr){
          mcmatch_tmp_dr = tmpDR;
          idx=iCand;
      } 
    }
    if (mcmatch_tmp_dr < dR_angle) isGenMatched.push_back(idx);
    else isGenMatched.push_back(-1);
  }
  return isGenMatched;
}

RVec<Float_t> makeResolution(RVec<Float_t> &GenMuonBare_pt, RVec<Float_t> &goodMuons_pt, RVec<Int_t> &idx)
{
  RVec<Float_t> res;
  for (unsigned int i=0; i!=idx.size(); i++) {
    if (idx[i]!=-1) {
      res.push_back(-(goodMuons_pt[idx[i]]-GenMuonBare_pt[i])/goodMuons_pt[idx[i]]);
    }
  }
  return res;
}

RVec<Int_t> makecharge(RVec<Int_t> &pdgid)
{
  RVec<Int_t> charge;
  for (unsigned int i=0; i!=pdgid.size(); i++) {
    Int_t Charge = pdgid[i] > 0 ? -1 : 1;
    charge.push_back(Charge);
  }
  return charge;
}
