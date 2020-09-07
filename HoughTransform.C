void HoughTransform() {
	TFile *f = new TFile("data.root", "READ");
	TTree *tree = (TTree *)f->Get("tree");

	Int_t eventId;
	TLorentzVector *position = new TLorentzVector(0., 0., 0., 0.);
	Double_t momentumMag;
	std::string *particleName = new std::string("");
	tree->SetBranchAddress("eventId", &eventId);
	tree->SetBranchAddress("position", &position);
	tree->SetBranchAddress("momentumMag", &momentumMag);
	tree->SetBranchAddress("particleName", &particleName);

	TH3I *h = new TH3I("h", "h", 200, 6500, 8500, 200, -1000, 1000, 50, 50, 250);
	TH2I *hr = new TH2I("h", "h", 200, 6500, 8500, 50, 50, 250);
	TH2I *h2 = new TH2I("h2", "h2", 200, 6500, 8500, 200, -1000, 1000);
	TH2I *h3 = new TH2I("h3", "h3", 200, 6500, 8500, 200, -1000, 1000);
	for(int j=0; j < tree->GetEntries(); ++j) {
		tree->GetEntry(j);
		if(momentumMag < 25) continue; 
//		if(momentumMag < 40) continue;
//		if(momentumMag > 45) continue; 
		if(particleName->compare("e-") !=0) continue;
		Float_t x = position->Z();
		Float_t y = position->Y();
		h3->Fill(x, y);
		for(int r=90; r < 160; ++r ) {
			for(int p=0; p < 400; ++p) {
				Float_t angle = TMath::Pi()*(-1+p/200.);
				Float_t m = x + r * TMath::Cos(angle);
				Float_t n = y + r * TMath::Sin(angle);
				if(abs(m-x)>20 && abs(n-y)>20) {
					h->Fill(m, n, r);
					h2->Fill(m, n);
				}
			}
		}
	}
	Int_t threshold = 1500;
	for(int c=0; c < h->GetSize(); ++c) {
		if(h2->GetBinContent(c) < threshold) {
			h2->SetBinContent(c, 0);
		}
	}
	TSpectrum2 *s = new TSpectrum2();
	TCanvas *c = new TCanvas("c", "c");
	c->Divide(2, 2);
	c->cd(1);
	//TODO: How to decide search threshold sigma
	Int_t nfound = s->Search(h2, 5, "colz");
	//searching good and ghost peaks (approximation)
	Double_t *xpeaks = s->GetPositionX();
	Double_t *ypeaks = s->GetPositionY();
	for (Int_t pf=0;pf<nfound;pf++) {
		std::cout << xpeaks[pf] << " " << ypeaks[pf] << std::endl;
	}
	h3->Draw("SAME");
	c->cd(2);
	h->Project3D("yz")->Draw("colz");
	c->cd(3);
	h->Project3D("zx")->Draw("colz");
	c->cd(4);
	c->Draw();
	
}
