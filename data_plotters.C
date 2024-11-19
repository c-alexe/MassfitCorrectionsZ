// -------------------------------------------------------------------------------------------------
// Plot data to corrected MC mass ratio after a given iteration, with eta, pT cuts for the muons
// -------------------------------------------------------------------------------------------------

void plot_data_mc_ratio( TString tag = "PostVFP", TString run = "Iter0", TString selection = "NONE" ) {
// run IterX means the sum of corrections derived in Iters 0 to (X-1) were applied

  TString plotname = TString("ratio_") + tag + TString("_") + run + TString("_") + selection + TString(".png");
  
  TCanvas* c = new TCanvas("c", "canvas", 600, 600);
  
  // Input mass files for data at a given iteration
  TFile* fsIter = TFile::Open("massscales_"+tag+"_"+run+".root", "READ");
  TString mc_name    = "h_smear0_bin_m"; // MC with corrections from previous iterations
  TString data_name  = "h_data_bin_m";
  // X-axis are 4D bins, Y-axis are masses; get MC and data
  TH2D* h2mass_mc = (TH2D*)fsIter->Get(mc_name);
  TH2D* h2mass_data = (TH2D*)fsIter->Get(data_name);

  // Muon eta pT binning
  TH1D* eta_edges = (TH1D*)fsIter->Get("h_eta_edges");
  TH1D* pt_edges = (TH1D*)fsIter->Get("h_pt_edges");

  unsigned int n_eta_bins = eta_edges->GetXaxis()->GetNbins();
  unsigned int n_pt_bins  = pt_edges->GetXaxis()->GetNbins();
  
  unsigned int ibin = 0; // Keep track of 4D bins
  for(unsigned int ieta_p = 0; ieta_p<n_eta_bins; ieta_p++) {
    for(unsigned int ipt_p = 0; ipt_p<n_pt_bins; ipt_p++) {
      for(unsigned int ieta_m = 0; ieta_m<n_eta_bins; ieta_m++) {
	      for(unsigned int ipt_m = 0; ipt_m<n_pt_bins; ipt_m++) {
	        float etap = TMath::Abs(eta_edges->GetXaxis()->GetBinCenter(ieta_p+1));
	        float etam = TMath::Abs(eta_edges->GetXaxis()->GetBinCenter(ieta_m+1));
	        float ptp = TMath::Abs(pt_edges->GetXaxis()->GetBinCenter(ipt_p+1));
	        float ptm = TMath::Abs(pt_edges->GetXaxis()->GetBinCenter(ipt_m+1));
	        bool accept = true;
          // Muon eta or pT selection cuts
	        if(selection=="CC")
	          accept = etap<1.5 && etam<1.5;
	        else if(selection=="FC")
	          accept = (etap<1.5 && etam>1.5) || (etap>1.5 && etam<1.5);
	        else if(selection=="FF")
	          accept = etap>1.5 && etam>1.5;
	        else if(selection=="LL")
	          accept =  (ptp<35 && ptm<35);
	        else if(selection=="LH")
	          accept =  (ptp<35 && ptm>35) || (ptp>35 && ptm<35);
	        else if(selection=="HH")
	          accept =  (ptp>35 && ptm>35);
	        if( !accept ) {
	          for(unsigned int iy=0; iy<h2mass_mc->GetYaxis()->GetNbins(); iy++) {
	            h2mass_mc->SetBinContent(ibin+1, iy+1, 0.0);
	            h2mass_data->SetBinContent(ibin+1, iy+1, 0.0);
	          }
	        }
	        ibin++;
	      }
      }
    }
  }

  // Get mass inclusive in all 4D bins for MC and data
  TH1D* hmass_mc = h2mass_mc->ProjectionY("mc");
  TH1D* hmass_data = h2mass_data->ProjectionY("data");

  // Scale MC to data
  hmass_mc->Scale(hmass_data->Integral()/hmass_mc->Integral());
  // Draw data to MC mass ratio
  hmass_data->Divide(hmass_mc);
  hmass_data->SetLineColor(kBlack);
  hmass_data->SetMaximum(1.1);
  hmass_data->SetMinimum(0.9);
  hmass_data->SetStats(0);
  hmass_data->SetTitle(tag+", "+run+", "+selection);

  c->cd();
  hmass_data->Draw("HISTE");

  c->SaveAs(plotname);
  
}

// -------------------------------------------------------------------------------------------------------------------
// Produce a TTree with input and the sum of fitted bias parameters from different iterations for a given toy or data 
// Produce a canvas with 4 panels for multiple iterations of a single toy or of data
// Panel 1 -> (pseudo)data to corrected MC mass ratio after the latest iteration
// Panels 2-4 -> sum of fitted A, e or M parameters at different iterations for the same toy or data
// -------------------------------------------------------------------------------------------------------------------

void merge_massloop(TString tag = "SmearRealistic_toy0", TString name = "massloop_in_iter2", bool savePng=true, bool isData=false) {
// run IterX means the sum of corrections derived in Iters 0 to (X-1) were applied

  TString plotname = name + TString("_") + tag;
  
  TCanvas* c = new TCanvas("c", "canvas", 1600, 400);
  c->Divide(4,1);

  // Output file and TTree
  TFile* fout = TFile::Open("merge_massloop_"+tag+"_"+name+".root", "RECREATE");
  TTree* treeout = new TTree("tree","");
  
  // Draw first canvas panel -> (pseudo)data to corrected MC mass ratio at the latest iteration

  // Input mass files for the same toy or data, different iterations
  TFile* fsIter0 = TFile::Open("massscales_"+tag+"_Iter0.root", "READ");
  TFile* fsIter1 = TFile::Open("massscales_"+tag+"_Iter1.root", "READ");
  TFile* fsIter2 = TFile::Open("massscales_"+tag+"_Iter2.root", "READ");
  TFile* fsIter3 = TFile::Open("massscales_"+tag+"_Iter3.root", "READ");
  TFile* fsIter4 = TFile::Open("massscales_"+tag+"_Iter4.root", "READ");
  TFile* fsIter5 = TFile::Open("massscales_"+tag+"_Iter5.root", "READ");
  TFile* fsIter6 = TFile::Open("massscales_"+tag+"_Iter6.root", "READ");
  TFile* fsIter7 = TFile::Open("massscales_"+tag+"_Iter7.root", "READ");
  TFile* fsIter8 = TFile::Open("massscales_"+tag+"_Iter8.root", "READ");
  
  TH1D* hmass_smear0 = 0; // use for MC with corrections from previous iterations
  TH1D* hmass_smear1 = 0; // use for (pseudo)data
  TString data_name = isData ? "h_data_bin_m" : "h_smear1_bin_m"; // name 
  TString data_namep = isData ? "data" : "smear1"; // title

  // Get mass inclusive in all 4D bins for MC and (pseudo)data for the latest iteration
  if(string(plotname.Data()).find("iter0")!=string::npos) {
    hmass_smear0 = ((TH2D*)fsIter0->Get("h_smear0_bin_m"))->ProjectionY("smear0"); // get MC
    hmass_smear1 = ((TH2D*)fsIter0->Get(data_name))->ProjectionY(data_namep); // get (pseudo)data
  } else if(string(plotname.Data()).find("iter1")!=string::npos) {
    hmass_smear0 = ((TH2D*)fsIter1->Get("h_smear0_bin_m"))->ProjectionY("smear0");
    hmass_smear1 = ((TH2D*)fsIter1->Get(data_name))->ProjectionY(data_namep);
  } else if(string(plotname.Data()).find("iter2")!=string::npos) {
    hmass_smear0 = ((TH2D*)fsIter2->Get("h_smear0_bin_m"))->ProjectionY("smear0");
    hmass_smear1 = ((TH2D*)fsIter2->Get(data_name))->ProjectionY(data_namep);
  } else if(string(plotname.Data()).find("iter3")!=string::npos && fsIter3!=0) {
    hmass_smear0 = ((TH2D*)fsIter3->Get("h_smear0_bin_m"))->ProjectionY("smear0");
    hmass_smear1 = ((TH2D*)fsIter3->Get(data_name))->ProjectionY(data_namep);
  } else if(string(plotname.Data()).find("iter4")!=string::npos && fsIter4!=0) {
    hmass_smear0 = ((TH2D*)fsIter4->Get("h_smear0_bin_m"))->ProjectionY("smear0");
    hmass_smear1 = ((TH2D*)fsIter4->Get(data_name))->ProjectionY(data_namep);
  } else if(string(plotname.Data()).find("iter5")!=string::npos && fsIter5!=0) {
    hmass_smear0 = ((TH2D*)fsIter5->Get("h_smear0_bin_m"))->ProjectionY("smear0");
    hmass_smear1 = ((TH2D*)fsIter5->Get(data_name))->ProjectionY(data_namep);
  } else if(string(plotname.Data()).find("iter6")!=string::npos && fsIter6!=0) {
    hmass_smear0 = ((TH2D*)fsIter6->Get("h_smear0_bin_m"))->ProjectionY("smear0");
    hmass_smear1 = ((TH2D*)fsIter6->Get(data_name))->ProjectionY(data_namep);
  } else if(string(plotname.Data()).find("iter7")!=string::npos && fsIter7!=0) {
    hmass_smear0 = ((TH2D*)fsIter7->Get("h_smear0_bin_m"))->ProjectionY("smear0");
    hmass_smear1 = ((TH2D*)fsIter7->Get(data_name))->ProjectionY(data_namep);
  } else if(string(plotname.Data()).find("iter8")!=string::npos && fsIter8!=0) {
    hmass_smear0 = ((TH2D*)fsIter8->Get("h_smear0_bin_m"))->ProjectionY("smear0");
    hmass_smear1 = ((TH2D*)fsIter8->Get(data_name))->ProjectionY(data_namep);
  }
  // Scale MC to (pseudo)data
  hmass_smear0->Scale(hmass_smear1->Integral()/hmass_smear0->Integral());
  // Draw (pseudo)data to MC mass ratio
  hmass_smear1->Divide(hmass_smear0);
  hmass_smear1->SetLineColor(kBlack);
  hmass_smear1->SetMaximum(1.1);
  hmass_smear1->SetMinimum(0.9);
  c->cd(1);
  hmass_smear1->SetStats(0);
  hmass_smear1->SetTitle("data / nominal");
  hmass_smear1->Draw("HISTE");
  
  // Draw remaining canvas panels -> sum of fitted A, e or M parameters at different iterations for the same toy or data 
  
  // Input massfit files for the same toy or data, different iterations
  TFile* fIter0 = TFile::Open("massfit_"+tag+"_Iter0.root", "READ");
  TFile* fIter1 = TFile::Open("massfit_"+tag+"_Iter1.root", "READ");
  TFile* fIter2 = TFile::Open("massfit_"+tag+"_Iter2.root", "READ");
  TFile* fIter3 = TFile::Open("massfit_"+tag+"_Iter3.root", "READ");
  TFile* fIter4 = TFile::Open("massfit_"+tag+"_Iter4.root", "READ");
  TFile* fIter5 = TFile::Open("massfit_"+tag+"_Iter5.root", "READ");
  TFile* fIter6 = TFile::Open("massfit_"+tag+"_Iter6.root", "READ");
  TFile* fIter7 = TFile::Open("massfit_"+tag+"_Iter7.root", "READ");
  TFile* fIter8 = TFile::Open("massfit_"+tag+"_Iter8.root", "READ");

  // Use either external or internal parameters
  vector<TString> params = {"A", "e", "M"};
  if(string(name.Data()).find("in")!=string::npos) params = {"Ain", "ein", "Min"};

  // Input massfit TTrees with fitted parameters at different iterations for the same toy or data   
  TTree* t0 = (TTree*) fIter0->Get("tree");
  TTree* t1 = (TTree*) fIter1->Get("tree");
  TTree* t2 = (TTree*) fIter2->Get("tree");
  TTree* t3 = 0;
  if(fIter3!=0) t3 = (TTree*) fIter3->Get("tree");
  TTree* t4 = 0;
  if(fIter4!=0) t4 = (TTree*) fIter4->Get("tree");
  TTree* t5 = 0;
  if(fIter5!=0) t5 = (TTree*) fIter5->Get("tree");
  TTree* t6 = 0;
  if(fIter6!=0) t6 = (TTree*) fIter6->Get("tree");
  TTree* t7 = 0;
  if(fIter7!=0) t7 = (TTree*) fIter7->Get("tree");
  TTree* t8 = 0;
  if(fIter8!=0) t8 = (TTree*) fIter8->Get("tree");
  
  // Read from TTree
  double fmin0, prob0;
  double fmin1, prob1;
  double fmin2, prob2;
  double fmin3, prob3;
  double fmin4, prob4;
  double fmin5, prob5;
  double fmin6, prob6;
  double fmin7, prob7;
  double fmin8, prob8;
  t0->SetBranchAddress("fmin", &fmin0);
  t0->SetBranchAddress("prob", &prob0);
  t1->SetBranchAddress("fmin", &fmin1);
  t1->SetBranchAddress("prob", &prob1);
  t2->SetBranchAddress("fmin", &fmin2);
  t2->SetBranchAddress("prob", &prob2);
  t0->GetEntry(0);
  t1->GetEntry(0);
  t2->GetEntry(0);
  if(fIter3!=0) {
    t3->SetBranchAddress("fmin", &fmin3);
    t3->SetBranchAddress("prob", &prob3);
    t3->GetEntry(0);
  }
  if(fIter4!=0) {
    t4->SetBranchAddress("fmin", &fmin4);
    t4->SetBranchAddress("prob", &prob4);
    t4->GetEntry(0);
  }
  if(fIter5!=0) {
    t5->SetBranchAddress("fmin", &fmin5);
    t5->SetBranchAddress("prob", &prob5);
    t5->GetEntry(0);
  }
  if(fIter6!=0) {
    t6->SetBranchAddress("fmin", &fmin6);
    t6->SetBranchAddress("prob", &prob6);
    t6->GetEntry(0);
  }
  if(fIter7!=0) {
    t7->SetBranchAddress("fmin", &fmin7);
    t7->SetBranchAddress("prob", &prob7);
    t7->GetEntry(0);
  }
  if(fIter8!=0) {
    t8->SetBranchAddress("fmin", &fmin8);
    t8->SetBranchAddress("prob", &prob8);
    t8->GetEntry(0);
  }

  // Histogram with the input values of the biases to use as format template
  TH1D* hp_nom_template  = (TH1D*) fIter0->Get("h_A_vals_nom");
  double A_val_fits[hp_nom_template->GetNbinsX()];
  double A_val_noms[hp_nom_template->GetNbinsX()];
  double A_val_errs[hp_nom_template->GetNbinsX()];
  double e_val_fits[hp_nom_template->GetNbinsX()];
  double e_val_noms[hp_nom_template->GetNbinsX()];
  double e_val_errs[hp_nom_template->GetNbinsX()];
  double M_val_fits[hp_nom_template->GetNbinsX()];
  double M_val_noms[hp_nom_template->GetNbinsX()];
  double M_val_errs[hp_nom_template->GetNbinsX()];

  // For each parameter type A,e,M
  for(unsigned int p = 0; p < params.size(); p++) {
    // Read nominal and fitted bias parameters from different iterations for the same toy or data 
    TH1D* hp_nom  = (TH1D*) fIter0->Get("h_"+params[p]+"_vals_nom");
    //hp_nom->Scale(0.);
    TH1D* hp_fit0 = (TH1D*) fIter0->Get("h_"+params[p]+"_vals_fit");
    TH1D* hp_fit1 = (TH1D*) fIter1->Get("h_"+params[p]+"_vals_fit");
    TH1D* hp_fit2 = (TH1D*) fIter2->Get("h_"+params[p]+"_vals_fit");
    TH1D* hp_fit3 = 0;
    if(fIter3!=0) hp_fit3 = (TH1D*) fIter3->Get("h_"+params[p]+"_vals_fit");
    TH1D* hp_fit4 = 0;
    if(fIter4!=0) hp_fit4 = (TH1D*) fIter4->Get("h_"+params[p]+"_vals_fit");
    TH1D* hp_fit5 = 0;
    if(fIter5!=0) hp_fit5 = (TH1D*) fIter5->Get("h_"+params[p]+"_vals_fit");
    TH1D* hp_fit6 = 0;
    if(fIter6!=0) hp_fit6 = (TH1D*) fIter6->Get("h_"+params[p]+"_vals_fit");
    TH1D* hp_fit7 = 0;
    if(fIter7!=0) hp_fit7 = (TH1D*) fIter7->Get("h_"+params[p]+"_vals_fit");
    TH1D* hp_fit8 = 0;
    if(fIter8!=0) hp_fit8 = (TH1D*) fIter8->Get("h_"+params[p]+"_vals_fit");
    
    // Add results for the fitted parameter of interest per eta bin from different iterations for the same toy or data 
    // The addition of fitted bias parameters makes sense since we correct the MC to better match (pseudo)data at each iter
    // The error is the error on the fitted parameter from the most recently added massfit

    hp_fit0->SetTitle(Form("#chi^{2}/ndof = %.2f (prob=%.2f)", 1+fmin0, prob0));

    if( string(plotname.Data()).find("iter1")!=string::npos ) {
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit0->SetBinError(ib, hp_fit1->GetBinError(ib) );
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit1->SetBinError(ib, 0.);
      hp_fit0->Add(hp_fit1);
      // The title shows the chi^2/ndof of the most recently added massfit
      hp_fit0->SetTitle(Form("#chi^{2}/ndof = %.2f (prob=%.2f)", 1+fmin1, prob1));
    } else if( string(plotname.Data()).find("iter2")!=string::npos ) {
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit0->SetBinError(ib, hp_fit2->GetBinError(ib) );
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit1->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit2->GetXaxis()->GetNbins();ib++) hp_fit2->SetBinError(ib, 0.);
      hp_fit0->Add(hp_fit1);
      hp_fit0->Add(hp_fit2);
      hp_fit0->SetTitle(Form("#chi^{2}/ndof = %.2f (prob=%.2f)", 1+fmin2, prob2));
    } else if( string(plotname.Data()).find("iter3")!=string::npos && fIter3!=0 ) {
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit0->SetBinError(ib, hp_fit3->GetBinError(ib) );
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit1->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit2->GetXaxis()->GetNbins();ib++) hp_fit2->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit2->GetXaxis()->GetNbins();ib++) hp_fit3->SetBinError(ib, 0.);
      hp_fit0->Add(hp_fit1);
      hp_fit0->Add(hp_fit2);
      hp_fit0->Add(hp_fit3);
      hp_fit0->SetTitle(Form("#chi^{2}/ndof = %.2f (prob=%.2f)", 1+fmin3, prob3));
    } else if( string(plotname.Data()).find("iter4")!=string::npos && fIter4!=0 ) {
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit0->SetBinError(ib, hp_fit4->GetBinError(ib) );
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit1->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit2->GetXaxis()->GetNbins();ib++) hp_fit2->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit3->GetXaxis()->GetNbins();ib++) hp_fit3->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit4->GetXaxis()->GetNbins();ib++) hp_fit4->SetBinError(ib, 0.);
      hp_fit0->Add(hp_fit1);
      hp_fit0->Add(hp_fit2);
      hp_fit0->Add(hp_fit3);
      hp_fit0->Add(hp_fit4);
      hp_fit0->SetTitle(Form("#chi^{2}/ndof = %.2f (prob=%.2f)", 1+fmin4, prob4));
    } else if( string(plotname.Data()).find("iter5")!=string::npos && fIter5!=0 ) {
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit0->SetBinError(ib, hp_fit5->GetBinError(ib) );
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit1->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit2->GetXaxis()->GetNbins();ib++) hp_fit2->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit3->GetXaxis()->GetNbins();ib++) hp_fit3->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit4->GetXaxis()->GetNbins();ib++) hp_fit4->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit4->GetXaxis()->GetNbins();ib++) hp_fit5->SetBinError(ib, 0.);
      hp_fit0->Add(hp_fit1);
      hp_fit0->Add(hp_fit2);
      hp_fit0->Add(hp_fit3);
      hp_fit0->Add(hp_fit4);
      hp_fit0->Add(hp_fit5);
      hp_fit0->SetTitle(Form("#chi^{2}/ndof = %.2f (prob=%.2f)", 1+fmin5, prob5));
    } else if( string(plotname.Data()).find("iter6")!=string::npos && fIter6!=0 ) {
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit0->SetBinError(ib, hp_fit6->GetBinError(ib) );
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit1->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit2->GetXaxis()->GetNbins();ib++) hp_fit2->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit3->GetXaxis()->GetNbins();ib++) hp_fit3->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit4->GetXaxis()->GetNbins();ib++) hp_fit4->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit4->GetXaxis()->GetNbins();ib++) hp_fit5->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit4->GetXaxis()->GetNbins();ib++) hp_fit6->SetBinError(ib, 0.);
      hp_fit0->Add(hp_fit1);
      hp_fit0->Add(hp_fit2);
      hp_fit0->Add(hp_fit3);
      hp_fit0->Add(hp_fit4);
      hp_fit0->Add(hp_fit5);
      hp_fit0->Add(hp_fit6);
      hp_fit0->SetTitle(Form("#chi^{2}/ndof = %.2f (prob=%.2f)", 1+fmin6, prob6));
    } else if( string(plotname.Data()).find("iter7")!=string::npos && fIter7!=0 ) {
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit0->SetBinError(ib, hp_fit7->GetBinError(ib) );
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit1->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit2->GetXaxis()->GetNbins();ib++) hp_fit2->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit3->GetXaxis()->GetNbins();ib++) hp_fit3->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit4->GetXaxis()->GetNbins();ib++) hp_fit4->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit4->GetXaxis()->GetNbins();ib++) hp_fit5->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit4->GetXaxis()->GetNbins();ib++) hp_fit6->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit4->GetXaxis()->GetNbins();ib++) hp_fit7->SetBinError(ib, 0.);
      hp_fit0->Add(hp_fit1);
      hp_fit0->Add(hp_fit2);
      hp_fit0->Add(hp_fit3);
      hp_fit0->Add(hp_fit4);
      hp_fit0->Add(hp_fit5);
      hp_fit0->Add(hp_fit6);
      hp_fit0->Add(hp_fit7);
      hp_fit0->SetTitle(Form("#chi^{2}/ndof = %.2f (prob=%.2f)", 1+fmin7, prob7));
    } else if( string(plotname.Data()).find("iter8")!=string::npos && fIter8!=0 ) {
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit0->SetBinError(ib, hp_fit8->GetBinError(ib) );
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit1->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit2->GetXaxis()->GetNbins();ib++) hp_fit2->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit3->GetXaxis()->GetNbins();ib++) hp_fit3->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit4->GetXaxis()->GetNbins();ib++) hp_fit4->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit4->GetXaxis()->GetNbins();ib++) hp_fit5->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit4->GetXaxis()->GetNbins();ib++) hp_fit6->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit4->GetXaxis()->GetNbins();ib++) hp_fit7->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit4->GetXaxis()->GetNbins();ib++) hp_fit8->SetBinError(ib, 0.);
      hp_fit0->Add(hp_fit1);
      hp_fit0->Add(hp_fit2);
      hp_fit0->Add(hp_fit3);
      hp_fit0->Add(hp_fit4);
      hp_fit0->Add(hp_fit5);
      hp_fit0->Add(hp_fit6);
      hp_fit0->Add(hp_fit7);
      hp_fit0->Add(hp_fit8);
      hp_fit0->SetTitle(Form("#chi^{2}/ndof = %.2f (prob=%.2f)", 1+fmin8, prob8));
    } // Finish if over how many iterations to add
    
    // Draw superimposed the nominal and the sum of fitted bias parameters from many iterations for the same toy or data
    hp_nom->SetLineColor(kBlue);
    hp_nom->SetLineWidth(3);
    hp_fit0->SetLineColor(kBlack);
    hp_fit0->SetMarkerColor(kBlack);
    hp_fit0->SetMarkerStyle(kFullCircle);
    c->cd(p+2);
    hp_nom->SetStats(0);
    //hp_nom->SetMaximum( hp_nom->GetMaximum()*1.5);
    //hp_nom->SetMinimum( hp_nom->GetMinimum()*0.5);
    hp_fit0->SetStats(0);
    //hp_fit0->SetMaximum( hp_fit0->GetMaximum()*1.5);
    //hp_fit0->SetMinimum( hp_fit0->GetMinimum()*0.5);

    //if(p>=1) {
    //  hp_fit0->SetMaximum( + (hp_fit0->GetMaximum()+hp_fit0->GetBinError(1))*1.1 );
    //  hp_fit0->SetMinimum( - (hp_fit0->GetMaximum()+hp_fit0->GetBinError(1))*1.1 );
    //}

    if(p==0 && savePng) { // A    
      hp_fit0->SetMaximum( +0.002 );
      hp_fit0->SetMinimum( -0.002 );      
    }
    if(p==1 && savePng) { // e
      hp_fit0->SetMaximum( +0.0025 );
      hp_fit0->SetMinimum( -0.0025 );
    }
    if(p==2 && savePng) { // M
      hp_fit0->SetMaximum( +0.0025 );
      hp_fit0->SetMinimum( -0.0025 );
    }          
    hp_fit0->Draw("HISTPE");
    hp_nom->Draw("HISTSAME");

    // Fill output TTrees with nominal, fitted and errors of parameters after adding multiple iterations for the same toy or data
    // Loop over eta bins
    for(int ib = 0; ib<hp_nom->GetNbinsX(); ib++) {
      // Branch names e.g. A0_intrue
      // Input bias (0 for data)
      TString b_nom_name(Form("%c%d%s", params[p][0], ib, string(params[p].Data()).find("in")!=string::npos ? Form("_intrue") : Form("_true") ));
      // Fitted bias, added from multiple iterations
      TString b_fit_name(Form("%c%d%s", params[p][0], ib, string(params[p].Data()).find("in")!=string::npos ? Form("_in") : Form("") ));
      // Fitted bias error, from the most recently added iteration
      TString b_err_name(Form("%c%d%s", params[p][0], ib, string(params[p].Data()).find("in")!=string::npos ? Form("_inerr") : Form("_err") )); 
      if(p==0) { // A
	      treeout->Branch( b_nom_name.Data(), &(A_val_noms[ib]), (b_nom_name+TString("/D")).Data() );
	      treeout->Branch( b_fit_name.Data(), &(A_val_fits[ib]), (b_fit_name+TString("/D")).Data() );
	      treeout->Branch( b_err_name.Data(), &(A_val_errs[ib]), (b_err_name+TString("/D")).Data() );
      } else if(p==1) { // e
	      treeout->Branch( b_nom_name.Data(), &(e_val_noms[ib]), (b_nom_name+TString("/D")).Data() );
	      treeout->Branch( b_fit_name.Data(), &(e_val_fits[ib]), (b_fit_name+TString("/D")).Data() );
	      treeout->Branch( b_err_name.Data(), &(e_val_errs[ib]), (b_err_name+TString("/D")).Data() );
      } else if(p==2) { // M
	      treeout->Branch( b_nom_name.Data(), &(M_val_noms[ib]), (b_nom_name+TString("/D")).Data() );
	      treeout->Branch( b_fit_name.Data(), &(M_val_fits[ib]), (b_fit_name+TString("/D")).Data() );
	      treeout->Branch( b_err_name.Data(), &(M_val_errs[ib]), (b_err_name+TString("/D")).Data() );
      }
    }

    // Loop over eta bins
    for(int ib = 0; ib<hp_nom->GetNbinsX(); ib++) { 
      if(p==0) { // A
	      A_val_noms[ib] = hp_nom->GetBinContent(ib+1);
	      A_val_fits[ib] = hp_fit0->GetBinContent(ib+1);
	      A_val_errs[ib] = hp_fit0->GetBinError(ib+1);
      } else if(p==1) { // e
 	      e_val_noms[ib] = hp_nom->GetBinContent(ib+1);
	      e_val_fits[ib] = hp_fit0->GetBinContent(ib+1);
	      e_val_errs[ib] = hp_fit0->GetBinError(ib+1);
      } else if(p==2) { // M
	      M_val_noms[ib] = hp_nom->GetBinContent(ib+1);
	      M_val_fits[ib] = hp_fit0->GetBinContent(ib+1);
 	      M_val_errs[ib] = hp_fit0->GetBinError(ib+1);
      }
    }

    fout->cd();
    hp_nom->Write(Form("h_%c_vals_nom", params[p][0]));
  } // Finish loop over parameter types

  c->Update();
  c->Draw();
  if(savePng) c->SaveAs(plotname+".png");
  
  fIter0->Close();
  fIter1->Close();
  fIter2->Close();
  if(fIter3!=0) fIter3->Close();
  if(fIter4!=0) fIter4->Close();
  if(fIter5!=0) fIter5->Close();
  if(fIter6!=0) fIter6->Close();
  if(fIter7!=0) fIter7->Close();
  if(fIter8!=0) fIter8->Close();
  delete c;

  fout->cd();
  treeout->Fill();
  treeout->Write();
  fout->Close();
}