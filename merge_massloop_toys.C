void merge_massloop_toys(TString tag = "SmearRealisticRnd_merged", TString name = "massloop_merged_in", bool batch=false, bool savePng=false,
			 bool plotMean=true){
  
  TString plotname = name + TString("_") + tag + TString(plotMean ? "_mean" : "_sigma");

  TH1D* h_pIter0 = 0;
  TH1D* h_pIter1 = 0;
  TH1D* h_pIter2 = 0;

  TCanvas* c = new TCanvas("c", "canvas", 1600, 400);
  c->Divide(3,1);

  TLegend* leg1 = new TLegend(0.45, 0.65, 0.85, 0.90, "","brNDC");
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.045);
  leg1->SetFillColor(10);
  leg1->SetHeader("");

  TFile* fIter0 = TFile::Open("massfit_"+tag+"_Iter0.root", "READ");
  TFile* fIter1 = TFile::Open("massfit_"+tag+"_Iter1.root", "READ");
  TFile* fIter2 = TFile::Open("massfit_"+tag+"_Iter2.root", "READ");

  vector<TString> params = {"A", "e", "M"};
  
  for(unsigned int p = 0; p < params.size(); p++){


    TH1D* h_template = (TH1D*)fIter0->Get(Form("h_%s_vals_nom", params[p].Data()));
    int nparams = h_template->GetNbinsX();
    c->cd(p+1);
    h_pIter0 = new TH1D(Form("h_%s_Iter0", params[p].Data()), "", nparams, 0, nparams);   
    h_pIter1 = new TH1D(Form("h_%s_Iter1", params[p].Data()), "", nparams, 0, nparams);
    h_pIter2 = new TH1D(Form("h_%s_Iter2", params[p].Data()), "", nparams, 0, nparams);   

    for(int ip = 0; ip<nparams; ip++){
      h_pIter0->GetXaxis()->SetBinLabel(ip+1, h_template->GetXaxis()->GetBinLabel(ip+1));
      h_pIter1->GetXaxis()->SetBinLabel(ip+1, h_template->GetXaxis()->GetBinLabel(ip+1));
      h_pIter2->GetXaxis()->SetBinLabel(ip+1, h_template->GetXaxis()->GetBinLabel(ip+1));
    }
    
    TTree* t0 = (TTree*) fIter0->Get("tree");
    TTree* t1 = (TTree*) fIter1->Get("tree");
    TTree* t2 = (TTree*) fIter2->Get("tree");

    TH1D* hpulls = new TH1D(Form("hpulls%d",p), "", 40, -8, 8);
    for(int ip = 0; ip<nparams; ip++){
      
      TF1* gf = new TF1("gf","[0]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp( -0.5*(x-[1])*(x-[1])/[2]/[2] )", -8, 8);
      gf->SetNpx(10000);
      gf->SetParameter(0, 100.);
      gf->SetParameter(1, 0.);
      gf->SetParameter(2, 1.0 );      

      TString formula(Form("(%s%d_in - %s%d_intrue)/%s%d_inerr>>hpulls%d", params[p].Data(), ip, params[p].Data(), ip, params[p].Data(), ip, p));
      if( string(name.Data()).find("in")==string::npos )
	formula = TString(Form("(%s%d - %s%d_true)/%s%d_err>>hpulls%d", params[p].Data(), ip, params[p].Data(), ip, params[p].Data(), ip, p));
      float mean = 0.;
      float mean_err = 0.;
      float rms = 0.;
      float rms_err = 0.;
      t0->Draw(formula.Data(), "", "");
      hpulls->Fit("gf", "QR", "", -6, 6);
      mean = gf->GetParameter(1);
      mean_err = gf->GetParError(1);
      rms = TMath::Abs(gf->GetParameter(2));
      rms_err = gf->GetParError(2);
      //rms = hpulls->GetRMS();
      //rms_err = hpulls->GetRMSError();
      h_pIter0->SetBinContent(ip+1, plotMean ? mean : rms);
      h_pIter0->SetBinError(ip+1, plotMean ? mean_err : rms_err);
      hpulls->Reset();
	
      t1->Draw(formula.Data(), "", "");
      gf->SetParameter(0, 100.);
      gf->SetParameter(1, 0.);
      gf->SetParameter(2, 1.0 );      
      hpulls->Fit("gf", "QR", "", -3, 3);
      mean = gf->GetParameter(1);
      mean_err = gf->GetParError(1);
      rms = TMath::Abs(gf->GetParameter(2));
      rms_err = gf->GetParError(2);
      //rms = hpulls->GetRMS();
      //rms_err = hpulls->GetRMSError();
      h_pIter1->SetBinContent(ip+1, plotMean ? mean : rms);
      h_pIter1->SetBinError(ip+1, plotMean ? mean_err : rms_err);
      hpulls->Reset();

      t2->Draw(formula.Data(), "", "");
      gf->SetParameter(0, 100.);
      gf->SetParameter(1, 0.);
      gf->SetParameter(2, 1.0 );      
      hpulls->Fit("gf", "QR", "", -3, 3);
      mean = gf->GetParameter(1);
      mean_err = gf->GetParError(1);
      rms = TMath::Abs(gf->GetParameter(2));
      rms_err = gf->GetParError(2);
      //rms = hpulls->GetRMS();
      //rms_err = hpulls->GetRMSError();
      h_pIter2->SetBinContent(ip+1, plotMean ? mean : rms);
      h_pIter2->SetBinError(ip+1, plotMean ? mean_err : rms_err);
      hpulls->Reset();
      delete gf;
    }

    c->cd(p+1);
    h_pIter0->SetLineColor(kBlue);
    h_pIter1->SetLineColor(kGreen);
    h_pIter2->SetLineColor(kRed);
    h_pIter0->SetMarkerStyle(kFullSquare);
    h_pIter1->SetMarkerStyle(kFullTriangleUp);
    h_pIter2->SetMarkerStyle(kFullCircle);
    h_pIter0->SetMarkerColor(kBlue);
    h_pIter1->SetMarkerColor(kGreen);
    h_pIter2->SetMarkerColor(kRed);
    h_pIter0->SetMarkerSize(1.3);
    h_pIter1->SetMarkerSize(1.3);
    h_pIter2->SetMarkerSize(1.3);

    
    h_pIter0->SetTitle(Form("%s of pulls for param %s (%s)",
			    plotMean ? Form("Mean") : Form("#sigma"),
			    params[p].Data(),
			    string(name.Data()).find("in")==string::npos ? Form("out") : Form("in") ));    
    h_pIter0->SetStats(0);
    h_pIter0->SetMinimum(plotMean ? -2.0 : 0.25);
    h_pIter0->SetMaximum(plotMean ? +2.0 : 3.0);
    h_pIter0->Draw("E");
    h_pIter1->Draw("ESAME");
    h_pIter2->Draw("ESAME");
    if(p==0){
      leg1->AddEntry(h_pIter0, "Iter 0", "PL");
      leg1->AddEntry(h_pIter1, "Iter 1", "PL");
      leg1->AddEntry(h_pIter2, "Iter 2", "PL");
    }
    leg1->Draw();
    TF1* line = new TF1("line", plotMean ? "0.0" : "1.0", 0, nparams);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->Draw("same");
  }
    
  c->Update();
  c->Draw();
  if(savePng) c->SaveAs(plotname+".png");
  if(batch){
    fIter0->Close();
    fIter1->Close();
    fIter2->Close();
    delete c;
  }
}
