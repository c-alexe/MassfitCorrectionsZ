//---------------------------------------------------------------------------
// Plot Gaussian and Crystal Ball jacobians superimposed for scale and width
//---------------------------------------------------------------------------

void plot_jacobians() {

  gErrorIgnoreLevel = 6001;
  
  TFile* f = TFile::Open("massscales_PostVFP_Iter0.root");
  //std::unique_ptr<TFile> fout( TFile::Open("jac_superimposed.root", "RECREATE") ); 

  // X-axis is 4D bin, Y-axis is m_reco-m_gen
  TH2D* h_dm = (TH2D*)f->Get("h_smear0_bin_dm");
  
  // X-axis is 4D bin, Y-axis is smear0 (reco mass corrected according to previous iterations) weighted by gaussian or crystal ball jacobian weights
  TH2D* h_scale_g = (TH2D*)f->Get("h_smear0_bin_jac_scale");
  TH2D* h_scale_cb = (TH2D*)f->Get("h_smear0_bin_jac_scale_cb");
  TH2D* h_width_g = (TH2D*)f->Get("h_smear0_bin_jac_width");
  TH2D* h_width_cb = (TH2D*)f->Get("h_smear0_bin_jac_width_cb");

  TCanvas* c = new TCanvas("c", "canvas", 1200, 400);
  // Panels for scale or width jacobians
  c->Divide(2,1);

  auto leg1 = new TLegend(0.58, 0.68, 0.90, 0.90);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetFillColor(10);
  leg1->SetNColumns(1);
  leg1->SetHeader("");

  // Loop over 4D bins
  for(int ibin=0; ibin<h_scale_g->GetXaxis()->GetNbins(); ibin++) {
    // Get number of entries in reco mass difference histogram for this 4D bin
    TH1D* h_dm_proj = (TH1D*)h_dm->ProjectionY("h_smear0_bin_dm",ibin,ibin);
    double entries = h_dm_proj->Integral();

    // Project the gaussian and cb jacobians in 1D histograms
    // scale jacobians
    TH1D* h_scale_g_proj = (TH1D*)h_scale_g->ProjectionY("h_smear0_bin_jac_scale",ibin,ibin);
    if (h_scale_g_proj->Integral() == 0.) continue;
    TH1D* h_scale_cb_proj = (TH1D*)h_scale_cb->ProjectionY("h_smear0_bin_jac_scale_cb",ibin,ibin);
    double scale_max = h_scale_g_proj->GetBinContent(h_scale_g_proj->GetMaximumBin()) > h_scale_cb_proj->GetBinContent(h_scale_cb_proj->GetMaximumBin()) ? h_scale_g_proj->GetBinContent(h_scale_g_proj->GetMaximumBin()) : h_scale_cb_proj->GetBinContent(h_scale_cb_proj->GetMaximumBin());
    double scale_min = h_scale_g_proj->GetBinContent(h_scale_g_proj->GetMinimumBin()) < h_scale_cb_proj->GetBinContent(h_scale_cb_proj->GetMinimumBin()) ? h_scale_g_proj->GetBinContent(h_scale_g_proj->GetMinimumBin()) : h_scale_cb_proj->GetBinContent(h_scale_cb_proj->GetMinimumBin());
    // width jacobians
    TH1D* h_width_g_proj = (TH1D*)h_width_g->ProjectionY("h_smear0_bin_jac_width",ibin,ibin);
    TH1D* h_width_cb_proj = (TH1D*)h_width_cb->ProjectionY("h_smear0_bin_jac_width_cb",ibin,ibin);
    double width_max = h_width_g_proj->GetBinContent(h_width_g_proj->GetMaximumBin()) > h_width_cb_proj->GetBinContent(h_width_cb_proj->GetMaximumBin()) ? h_width_g_proj->GetBinContent(h_width_g_proj->GetMaximumBin()) : h_width_cb_proj->GetBinContent(h_width_cb_proj->GetMaximumBin());
    double width_min = h_width_g_proj->GetBinContent(h_width_g_proj->GetMinimumBin()) < h_width_cb_proj->GetBinContent(h_width_cb_proj->GetMinimumBin()) ? h_width_g_proj->GetBinContent(h_width_g_proj->GetMinimumBin()) : h_width_cb_proj->GetBinContent(h_width_cb_proj->GetMinimumBin());

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
    //fout->WriteObject(c, TString(("jac_"+to_string(ibin)).c_str()));
    c->SaveAs(Form("jac_plots/jac_%d.pdf", ibin));
    
    leg1->Clear();
  } // end loop over 4D bins

  // Execute outside of ROOT
  // cd jac_plots/; gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_jac.pdf -dBATCH jac_*.pdf; rm jac_*.pdf; cd ..
}

//---------------------------------------------------------------------------
// Plot m_smear0 or m_smear0 - m_gen vs m_gen in a 4D bin
//---------------------------------------------------------------------------

void plot_kernels() { 
  
  gStyle->SetOptStat(0);
  // Because this is Iter0, smear0=reco e.g. there are no corrections from previous iterations
  TFile* f = TFile::Open("massscales_PostVFP_Iter0.root");
  
  // X-axis is 4D bin, Y-axis is gen mass, Z-axis is m_smear0-m_gen
  TH3D* h_dm = (TH3D*)f->Get("h_smear0_bin_gm_dm");
  // X-axis is 4D bin, Y-axis is gen mass, Z-axis is m_smear0
  TH3D* h_m = (TH3D*)f->Get("h_smear0_bin_gm_m");
  
  TCanvas* c = new TCanvas("c", "canvas", 1200, 400);
  c->Divide(2,1);

  auto leg1 = new TLegend(0.58, 0.68, 0.90, 0.90);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetFillColor(10);
  leg1->SetNColumns(1);
  leg1->SetHeader("");
  
  // Loop over 4D bins
  for(int ibin=0; ibin<h_dm->GetXaxis()->GetNbins(); ibin++) {
    // Consider a single 4D bin
    h_dm->GetXaxis()->SetRange(ibin,ibin);
    // X-axis is m_gen, Y-axis is m_smear0-m_gen
    TH2D* h_dm_proj = (TH2D*)h_dm->Project3D("zy");
    
    // Skip 4D bins with low statistics
    double entries = h_dm_proj->Integral();
    if (entries < 800.) continue;
    
    h_m->GetXaxis()->SetRange(ibin,ibin);
    // X-axis is m_gen, Y-axis is m_smear0
    TH2D* h_m_proj = (TH2D*)h_m->Project3D("zy");
   
    c->cd(1);
    h_dm_proj->SetTitle(TString(("bin "+to_string(ibin)).c_str()));
    h_dm_proj->GetXaxis()->SetTitle("m_gen");
    h_dm_proj->GetYaxis()->SetTitle("m_smear0 - m_gen");
    h_dm_proj->Draw("COLZ");

    leg1->AddEntry((TObject*)0, TString((to_string(entries)+" events").c_str()), "");
    leg1->Draw("SAME");
    
    c->cd(2);
    h_m_proj->SetTitle(TString(("bin "+to_string(ibin)).c_str()));
    h_m_proj->GetXaxis()->SetTitle("m_gen");
    h_m_proj->GetYaxis()->SetTitle("m_smear0");
    h_m_proj->Draw("COLZ");
    c->cd();

    c->Draw();
    c->SaveAs(Form("kernel_plots/kernel_%d.pdf", ibin));

    leg1->Clear();
  } // end loop over 4D bins

  // Execute outside of ROOT
  // cd kernel_plots/; gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combined_kernel.pdf -dBATCH kernel_*.pdf; rm kernel_*.pdf; cd ..
  
}