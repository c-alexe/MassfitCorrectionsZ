#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "RooMsgService.h"

using namespace ROOT;
using namespace ROOT::VecOps;
using ROOT::RDF::RNode;

int quick_plotter(){

  gErrorIgnoreLevel = 6001;
   
  /*
  //ROOT::EnableImplicitMT();
  RDataFrame df("tree", "cb.root");

  auto h_edm = df.Histo1D({"edm", "edm", 100, -1, 20},"edm");
  auto h_status = df.Histo1D({"status", "status", 20, -10, 10},"status");
  auto h_flag = df.Histo1D({"flag", "flag", 4, 0, 4},"flag");
  
  std::unique_ptr<TFile> fout( TFile::Open("cb_tree_hist.root", "RECREATE") ); 
  h_edm->Write();
  h_status->Write();
  h_flag->Write();
  fout->Close();
  */
  
  TFile* f = TFile::Open("massscales_PostVFP_Iter0.root");
  std::unique_ptr<TFile> fout2( TFile::Open("jac_superimposed.root", "RECREATE") ); 
  // gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_jac.pdf -dBATCH jac_*.pdf

  TH2D* h_dm = (TH2D*)f->Get("h_reco_bin_dm");
  
  TH2D* h_scale_g = (TH2D*)f->Get("h_smear0_bin_jac_scale");
  TH2D* h_scale_cb = (TH2D*)f->Get("h_smear0_bin_jac_scale_cb");
  TH2D* h_width_g = (TH2D*)f->Get("h_smear0_bin_jac_width");
  TH2D* h_width_cb = (TH2D*)f->Get("h_smear0_bin_jac_width_cb");

  TCanvas* c = new TCanvas("c", "canvas", 1200, 400);
  c->Divide(2,1);

  auto leg1 = new TLegend(0.58, 0.68, 0.90, 0.90);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetFillColor(10);
  leg1->SetNColumns(1);
  leg1->SetHeader("");

  double scale_max, scale_min, width_max, width_min;
  double entries;
  
  for(int ibin=0; ibin<h_scale_g->GetXaxis()->GetNbins(); ibin++) {
    TH1D* h_dm_proj = (TH1D*)h_dm->ProjectionY("h_reco_bin_dm",ibin,ibin);
    entries = h_dm_proj->Integral();
    TH1D* h_scale_g_proj = (TH1D*)h_scale_g->ProjectionY("h_smear0_bin_jac_scale",ibin,ibin);
    if (h_scale_g_proj->Integral() == 0.) continue;
    TH1D* h_scale_cb_proj = (TH1D*)h_scale_cb->ProjectionY("h_smear0_bin_jac_scale_cb",ibin,ibin);
    scale_max = h_scale_g_proj->GetBinContent(h_scale_g_proj->GetMaximumBin()) > h_scale_cb_proj->GetBinContent(h_scale_cb_proj->GetMaximumBin()) ? h_scale_g_proj->GetBinContent(h_scale_g_proj->GetMaximumBin()) : h_scale_cb_proj->GetBinContent(h_scale_cb_proj->GetMaximumBin());
    scale_min = h_scale_g_proj->GetBinContent(h_scale_g_proj->GetMinimumBin()) < h_scale_cb_proj->GetBinContent(h_scale_cb_proj->GetMinimumBin()) ? h_scale_g_proj->GetBinContent(h_scale_g_proj->GetMinimumBin()) : h_scale_cb_proj->GetBinContent(h_scale_cb_proj->GetMinimumBin());
    TH1D* h_width_g_proj = (TH1D*)h_width_g->ProjectionY("h_smear0_bin_jac_width",ibin,ibin);
    TH1D* h_width_cb_proj = (TH1D*)h_width_cb->ProjectionY("h_smear0_bin_jac_width_cb",ibin,ibin);
    width_max = h_width_g_proj->GetBinContent(h_width_g_proj->GetMaximumBin()) > h_width_cb_proj->GetBinContent(h_width_cb_proj->GetMaximumBin()) ? h_width_g_proj->GetBinContent(h_width_g_proj->GetMaximumBin()) : h_width_cb_proj->GetBinContent(h_width_cb_proj->GetMaximumBin());
    width_min = h_width_g_proj->GetBinContent(h_width_g_proj->GetMinimumBin()) < h_width_cb_proj->GetBinContent(h_width_cb_proj->GetMinimumBin()) ? h_width_g_proj->GetBinContent(h_width_g_proj->GetMinimumBin()) : h_width_cb_proj->GetBinContent(h_width_cb_proj->GetMinimumBin());

    c->cd(1);
    h_scale_g_proj->SetTitle(TString(("bin "+to_string(ibin)).c_str()));
    h_scale_g_proj->SetLineColor(kBlack);
    h_scale_g_proj->SetMaximum(scale_max+100.);
    h_scale_g_proj->SetMinimum(scale_min-100.);
    h_scale_g_proj->Draw("HIST");
    h_scale_cb_proj->SetLineColor(kGreen);
    h_scale_cb_proj->Draw("HIST SAME");

    leg1->AddEntry(h_scale_g_proj, "Gaus", "l");
    leg1->AddEntry(h_scale_cb_proj, "CB", "l");
    leg1->AddEntry((TObject*)0, TString((to_string(entries)+" events").c_str()), "");
    leg1->Draw("SAME");
    
    c->cd(2);
    h_width_g_proj->SetTitle(TString(("bin "+to_string(ibin)).c_str()));
    h_width_g_proj->SetLineColor(kBlack);
    h_width_g_proj->SetMaximum(width_max+100.);
    h_width_g_proj->SetMinimum(width_min-100.);
    h_width_g_proj->Draw("HIST");
    h_width_cb_proj->SetLineColor(kGreen);
    h_width_cb_proj->Draw("HIST SAME");
    c->cd();

    c->Draw();
    //fout2->WriteObject(c, TString(("jac_"+to_string(ibin)).c_str()));
    c->SaveAs(Form("jac_plots/jac_%d.pdf", ibin));
    
    leg1->Clear();
  }
  
  /*
  gStyle->SetOptStat(0);
  TFile* f = TFile::Open("massscales_PostVFP_Iter0.root");
  
  TH3D* h_dm = (TH3D*)f->Get("h_reco_bin_gm_dm");
  TH3D* h_m = (TH3D*)f->Get("h_reco_bin_gm_m");

  double entries;
  
  TCanvas* c = new TCanvas("c", "canvas", 1200, 400);
  c->Divide(2,1);

  auto leg1 = new TLegend(0.58, 0.68, 0.90, 0.90);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetFillColor(10);
  leg1->SetNColumns(1);
  leg1->SetHeader("");
  
  for(int ibin=0; ibin<h_dm->GetXaxis()->GetNbins(); ibin++){
    h_dm->GetXaxis()->SetRange(ibin,ibin);
    TH2D* h_dm_proj = (TH2D*)h_dm->Project3D("zy");
    
    entries = h_dm_proj->Integral();
    if (entries < 800.) continue;
    
    h_m->GetXaxis()->SetRange(ibin,ibin);
    TH2D* h_m_proj = (TH2D*)h_m->Project3D("zy");
   
    c->cd(1);
    h_dm_proj->SetTitle(TString(("bin "+to_string(ibin)).c_str()));
    h_dm_proj->GetXaxis()->SetTitle("m_gen");
    h_dm_proj->GetYaxis()->SetTitle("m_diff");
    h_dm_proj->Draw("COLZ");

    leg1->AddEntry((TObject*)0, TString((to_string(entries)+" events").c_str()), "");
    leg1->Draw("SAME");
    
    c->cd(2);
    h_m_proj->SetTitle(TString(("bin "+to_string(ibin)).c_str()));
    h_m_proj->GetXaxis()->SetTitle("m_gen");
    h_m_proj->GetYaxis()->SetTitle("m_reco");
    h_m_proj->Draw("COLZ");
    c->cd();

    c->Draw();
    c->SaveAs(Form("kernel_plots/kernel_%d.pdf", ibin));

    leg1->Clear();
  }
  */
  
  return 0;
}
