Bool_t rescorr=false;

// pT range relevant for mass fits 
Double_t innercut=20.; 
Double_t outercut=60.; 

Double_t scalemodel(Double_t *x, Double_t *par)
{
    if ((abs(x[0])<innercut)||(abs(x[0])>outercut)) {
      TF1::RejectPoint();
      return -999.;
   }
   return par[0]+par[1]/abs(x[0])+par[2]*x[0];
}

Double_t resmodel_withcorr(Double_t *x, Double_t *par)
{
    if ((abs(x[0])<innercut)||(abs(x[0])>outercut)) {
      TF1::RejectPoint();
      return -999.;
   }
   return sqrt(par[0]*par[0]+par[1]*par[1]*x[0]*x[0]+par[2]*par[2]/(1+par[3]*par[3]/(x[0]*x[0])));
}

Double_t resmodel(Double_t *x, Double_t *par)
{
    if ((abs(x[0])<innercut)||(abs(x[0])>outercut)) {
      TF1::RejectPoint();
      return -999.;
   }
   return sqrt(par[0]*par[0]+par[1]*par[1]*x[0]*x[0]); // par[0] = a, par[1] = c, x is pT
}

void resolutionfitter() {
	TFile* file=new TFile("nominal_resolution_histograms_2018.root");
	TH3D* histo=(TH3D*)file->Get("histo");
	TFile* output=new TFile("nominal_resolution_coefficients_2018.root","RECREATE");
	output->cd();
	TH1D* scalea = new TH1D("scalea","",24,-2.4,2.4);
	TH1D* scalem = new TH1D("scalem","",24,-2.4,2.4);
	TH1D* scaleeps = new TH1D("scaleeps","",24,-2.4,2.4);
	TH1D* resa = new TH1D("resa","",24,-2.4,2.4);
	TH1D* resb = new TH1D("resb","",24,-2.4,2.4);
	TH1D* resc = new TH1D("resc","",24,-2.4,2.4);
	TH1D* resd = new TH1D("resd","",24,-2.4,2.4);
	for (unsigned int i=0; i!=histo->GetZaxis()->GetNbins(); i++) {
		histo->GetZaxis()->SetRange(i+1,i+1);
		float mineta=histo->GetZaxis()->GetBinLowEdge(i+1), maxeta=histo->GetZaxis()->GetBinUpEdge(i+1);
		TH1D* Histo=new TH1D((std::string("Histo")+std::to_string(i)).c_str(),(std::to_string(mineta)+std::string("<#eta<")+std::to_string(maxeta)).c_str(),histo->GetYaxis()->GetNbins(),histo->GetYaxis()->GetBinLowEdge(1),histo->GetYaxis()->GetBinUpEdge(histo->GetYaxis()->GetNbins()));
		TH1D* Histores=new TH1D((std::string("Histores")+std::to_string(i)).c_str(),(std::to_string(mineta)+std::string("<#eta<")+std::to_string(maxeta)).c_str(),histo->GetYaxis()->GetNbins(),histo->GetYaxis()->GetBinLowEdge(1),histo->GetYaxis()->GetBinUpEdge(histo->GetYaxis()->GetNbins()));
		for (unsigned int j=0; j!=histo->GetYaxis()->GetNbins(); j++) {
			histo->GetYaxis()->SetRange(j+1,j+1);
			float minpt=histo->GetYaxis()->GetBinLowEdge(j+1), maxpt=histo->GetYaxis()->GetBinUpEdge(j+1);
			std::string histname("histo_");
			histname+=std::to_string(i)+std::string("_")+std::to_string(j);
			TH1D* histo1=(TH1D*)histo->Project3D("x")->Clone(histname.c_str());
			if (histo1->GetEntries()<2000) histo1->Rebin(10);
			std::string histtitle;
			histtitle+=std::to_string(mineta)+std::string("<#eta<")+std::to_string(maxeta)+std::string(" ")+std::to_string(minpt)+std::string("<#pt<")+std::to_string(maxpt);
			histo1->SetTitle(histtitle.c_str());
			std::string model("((abs(x-[1])<=[3]*abs([2]))*[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2])))");
			model+=std::string("+((abs(x-[1])>[3]*abs([2]))*[0]*exp([3]*[3]/2)*exp(-[3]/abs([2])*abs(x-[1])))+[4]");
			int maxbin=histo1->GetMaximumBin();
			float hwhm;
			for (unsigned int h=maxbin; h!=histo1->GetXaxis()->GetNbins(); h++) {
				if (histo1->GetBinContent(h+1)>0.5*histo1->GetBinContent(maxbin)) hwhm=histo1->GetBinCenter(h+1)-histo1->GetBinCenter(maxbin);
			}
			float kfit=4, kfit2=4.5; // range for resolution fits
			auto fa1 = new TF1("fa1",model.c_str(),histo1->GetBinCenter(maxbin)-kfit*hwhm,histo1->GetBinCenter(maxbin)+kfit*hwhm);
			fa1->SetParameter(0,histo1->GetBinContent(maxbin));
			fa1->SetParameter(1,histo1->GetBinCenter(maxbin));
			fa1->SetParameter(2,hwhm);
			fa1->SetParameter(3,3); // exponential tails
			fa1->SetParameter(4,0);
			auto r1 = histo1->Fit(fa1, "LSQ", "", histo1->GetBinCenter(maxbin)-kfit*hwhm,histo1->GetBinCenter(maxbin)+kfit*hwhm);
			auto fa2 = new TF1("fa2",model.c_str(),histo1->GetBinCenter(maxbin)-kfit2*abs(r1->Parameter(2)),histo1->GetBinCenter(maxbin)+kfit2*abs(r1->Parameter(2)));
			fa2->SetParameter(0,r1->Parameter(0));
			fa2->SetParameter(1,r1->Parameter(1));
			fa2->SetParameter(2,abs(r1->Parameter(2)));
			fa2->SetParameter(3,r1->Parameter(3));
			fa2->SetParameter(4,r1->Parameter(4));
			auto r2 = histo1->Fit(fa2, "LSQ", "", histo1->GetBinCenter(maxbin)-kfit2*abs(r1->Parameter(2)),histo1->GetBinCenter(maxbin)+kfit2*abs(r1->Parameter(2)));
			Int_t fitStatus = r2;
			if ((fitStatus==0)&&(r2->ParError(1)<0.1)&&(abs(r2->Parameter(1)-histo1->GetBinCenter(maxbin))<0.05)) {
				//histo1->Write();
				Histo->SetBinContent(j+1,r2->Parameter(1));
				Histo->SetBinError(j+1,r2->ParError(1));
				Histores->SetBinContent(j+1,abs(r2->Parameter(2)));
				Histores->SetBinError(j+1,r2->ParError(2));
			}
		}
		Histo->SetMinimum(-0.01);
		Histo->SetMaximum(0.01);
		TF1 *Scalemodel = new TF1("scalemodel",scalemodel,histo->GetYaxis()->GetBinLowEdge(1),histo->GetYaxis()->GetBinUpEdge(histo->GetYaxis()->GetNbins()),3);
		TF1 *Resmodel;
		if (rescorr) Resmodel = new TF1("resmodel",resmodel_withcorr,histo->GetYaxis()->GetBinLowEdge(1),histo->GetYaxis()->GetBinUpEdge(histo->GetYaxis()->GetNbins()),4);
		else Resmodel = new TF1("resmodel",resmodel,histo->GetYaxis()->GetBinLowEdge(1),histo->GetYaxis()->GetBinUpEdge(histo->GetYaxis()->GetNbins()),2);
		Scalemodel->SetParameters(0,0.001);
		Scalemodel->SetParameters(1,0.001);
		Scalemodel->SetParameters(2,0.001);
		Resmodel->SetParameters(0,0.01);
		Resmodel->SetParameters(1,0.001);
		if (rescorr) {
			Resmodel->SetParameters(2,0.01);
			Resmodel->SetParameters(3,0.001);
		}
		auto Scaleresult = Histo->Fit("scalemodel","SQ");
		scalea->SetBinContent(i+1,Scaleresult->Parameter(0));
		scalea->SetBinError(i+1,Scaleresult->ParError(0));
		scaleeps->SetBinContent(i+1,Scaleresult->Parameter(1));
		scaleeps->SetBinError(i+1,Scaleresult->ParError(1));
		scalem->SetBinContent(i+1,Scaleresult->Parameter(2));
		scalem->SetBinError(i+1,Scaleresult->ParError(2));
		auto Resresult = Histores->Fit("resmodel","SQ");
		resa->SetBinContent(i+1,abs(Resresult->Parameter(0)));
		resa->SetBinError(i+1,Resresult->ParError(0));
		resc->SetBinContent(i+1,abs(Resresult->Parameter(1)));
		resc->SetBinError(i+1,Resresult->ParError(1));
		if (rescorr) {
			resb->SetBinContent(i+1,abs(Resresult->Parameter(2)));
			resb->SetBinError(i+1,Resresult->ParError(2));
			resd->SetBinContent(i+1,abs(Resresult->Parameter(3)));
			resd->SetBinError(i+1,Resresult->ParError(3));
		}
		Histo->Write();
		Histores->Write();
	}
	scalea->Write();
	scaleeps->Write();
	scalem->Write();
	resa->Write();
	resc->Write();
	if (rescorr) {
		resb->Write();
		resd->Write();
	}
	output->Close();
}
