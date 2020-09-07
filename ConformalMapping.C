void FindCenter(std::vector<Float_t> &trackX, std::vector<Float_t> &trackY) {
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

	TH2F *h = new TH2F("h", "h;u;v", 100, -2000, 2000, 100, -2000, 2000);
	TH2F *circ = new TH2F("circ", "circ;x;y", 100, -2000, 2000, 100, -2000, 2000);
	for(int j=0; j < tree->GetEntries(); ++j) {
		tree->GetEntry(j);
		//if entry is found in vector, skip
		if(std::find(trackX.begin(), trackX.end(), position->Z() - 7650 ) != trackX.end() && std::find(trackY.begin(), trackY.end(), position->Y() ) != trackY.end()) continue;
		if(momentumMag < 380) continue; 
		if(momentumMag > 415) continue; 
		if(particleName->compare("proton") !=0) continue;
		Float_t x = position->Z() - 7650;
		Float_t y = position->Y();
		circ->Fill(x, y);
		for(Int_t u=-2000; u < 2000; u = u+40) {
			//Hough transform
			Float_t v = -(x/y)*u + (x*x + y*y)/(2*y);
			//remove points within inner radius
			if(u*u +v*v < 496*496 ) continue;
			h->Fill(u, v);
		}
	}
	
	Int_t threshold = 10;
	for(int c=0; c < h->GetSize(); ++c) {
		if(h->GetBinContent(c) < threshold) {
			h->SetBinContent(c, 0);
		}
	}
	TSpectrum *s = new TSpectrum(1);
	TSpectrum *xs = new TSpectrum(1);
	Int_t nfound = s->Search(h->ProjectionY(), 2, "", 0.10);
	Double_t *ypeaks = s->GetPositionX();

	Int_t Xc, Yc;
	for(Int_t p=0; p < nfound; ++p) {
		Double_t xp = ypeaks[p];
		TH1D *h_px = h->ProjectionX("h_px", h->ProjectionX()->GetXaxis()->FindBin(xp), h->ProjectionX()->GetXaxis()->FindBin(xp));
		h_px->Rebin(2);
		xs->Search(h_px, 2, "", 0.10);
		Double_t *xpeaks = xs->GetPositionX();
		Yc = ypeaks[p];
		Xc = xpeaks[0];
		delete h_px;
	}

	TEllipse *innerwall = new TEllipse(0, 0, 496); innerwall->SetFillStyle(0); innerwall->SetLineWidth(1);
	TEllipse *outerwall = new TEllipse(0, 0, 840); outerwall->SetFillStyle(0); outerwall->SetLineWidth(1);
	TEllipse *target= new TEllipse(0, 0, 100); target->SetFillStyle(0); target->SetLineWidth(1);
	TEllipse *track = new TEllipse(Xc, Yc, TMath::Sqrt(Xc*Xc + Yc*Yc), 0);

	TCanvas *c = new TCanvas("c", "c", 1000, 1000);
	h->Draw("colz");
	circ->SetMarkerStyle(7);
	circ->Draw("SAME");

	innerwall->Draw("SAME");
	outerwall->Draw("SAME");
	target->Draw("SAME");

	track->SetFillStyle(0);
	track->SetLineWidth(1);
	track->SetLineColor(kRed);
	track->SetNoEdges(kTRUE);
	track->Draw("SAME");
	TMarker *mark = new TMarker(Xc, Yc, kFullCircle);
	mark->Draw("SAME");
	c->Draw();
	std::cout << "center: (" << Xc << ", " << Yc << ")" << std::endl;

	for(int j=0; j < tree->GetEntries(); ++j) {
		tree->GetEntry(j);
		if(momentumMag < 380) continue; 
		if(momentumMag > 415) continue; 
		if(particleName->compare("proton") !=0) continue;
		Float_t x = position->Z() - 7650;
		Float_t y = position->Y();
		Float_t deltaR = 5;
		if(abs((x-Xc)*(x-Xc) + (y-Yc)*(y-Yc) - (Xc-deltaR)*(Xc-deltaR) + (Yc-deltaR)*(Yc-deltaR) ) < 10e5 ) {
			trackX.push_back(x);
			trackY.push_back(y);
		}
	}
}
void ConformalMapping() {
	std::vector<Float_t> trackX, trackY;
	FindCenter(trackX, trackY);
//	FindCenter(trackX, trackY);
}
