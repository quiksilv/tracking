////
Double_t lowMomLimit = 370;
Double_t highMomLimit = 415;
Int_t findTrackLimit = 25;
Float_t deltaRsquared = 5e4;
Int_t conformalMappingThreshold = 20;
Double_t search_sigma = 1.2;
TH1F *hDelta = new TH1F("hDelta", "hDelta", 1000, -1e5, 1e5);
////
void FindCenter(std::vector<Float_t> &trackX, std::vector<Float_t> &trackY, std::vector<Float_t> &vecXc, std::vector<Float_t> &vecYc) {
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
		if(momentumMag < lowMomLimit) continue; 
		if(momentumMag > highMomLimit) continue; 
		if(particleName->compare("proton") !=0) continue;
		Float_t x = position->Z() - 7650;
		Float_t y = position->Y();
		//Float_t z = position->X() - 6450;
		circ->Fill(x, y);
		for(Int_t u=-2000; u < 2000; u = u+40) {
			//Hough transform
			Float_t v = -(x/y)*u + (x*x + y*y)/(2*y);
			//remove points within inner radius
			if(u*u + v*v < 496*496 ) continue;
			h->Fill(u, v);
		}
	}
	
	for(int c=0; c < h->GetSize(); ++c) {
		if(h->GetBinContent(c) < conformalMappingThreshold) {
			h->SetBinContent(c, 0);
		}
	}

	TSpectrum *s = new TSpectrum(10);
	TSpectrum *xs = new TSpectrum(10);
	Int_t nfound = s->Search(h->ProjectionY(), search_sigma, "", 0.05);
	Double_t *ypeaks = s->GetPositionX();
	Double_t xp = ypeaks[0];
	TH1D *h_px = h->ProjectionX("h_px", h->ProjectionX()->GetXaxis()->FindBin(xp), h->ProjectionX()->GetXaxis()->FindBin(xp) );
	xs->Search(h_px, search_sigma, "", 0.05);
	Double_t *xpeaks = xs->GetPositionX();
	Int_t Yc = ypeaks[0];
	Int_t Xc = xpeaks[0];
	//display evolution of the conformal mapping projection
//	TCanvas *temp = new TCanvas("temp", "temp");
//	temp->Divide(2, 2);
//	temp->cd(1);
//	h->Draw("colz");
//	temp->cd(2);
//	h_px->Draw();
//	temp->cd(3);
//	h->ProjectionY()->Draw();
//	temp->SaveAs("htemp.pdf");
	delete h_px;
	delete ypeaks;
	delete xpeaks;
	std::cout << "center: (" << Xc << ", " << Yc << ") radius: " << TMath::Sqrt(Xc*Xc + Yc*Yc) << std::endl;
	if(Xc!=0 && Yc!=0) {
		//push to vector and remove track hits only if nonzero track center is found
		vecXc.push_back(Xc);
		vecYc.push_back(Yc);

		TVector3 v2(Xc, Yc, 0);
		for(int j=0; j < tree->GetEntries(); ++j) {
			tree->GetEntry(j);
			if(momentumMag < lowMomLimit) continue; 
			if(momentumMag > highMomLimit) continue; 
			if(particleName->compare("proton") !=0) continue;
			Float_t x = position->Z() - 7650;
			Float_t y = position->Y();
			TVector3 v1(x, y, 0);
			TVector3 cprod = v1.Cross(v2);
			//check vicinity and direction of curl (positive or negative charged particles)
			hDelta->Fill((x-Xc)*(x-Xc) + (y-Yc)*(y-Yc) - Xc*Xc - Yc*Yc);
			if(abs((x-Xc)*(x-Xc) + (y-Yc)*(y-Yc) - Xc*Xc - Yc*Yc) < deltaRsquared && cprod.Z() >= 0) {
				//push points in the vicinity of the circle into the vector to be excluded in the next iteration
				trackX.push_back(x);
				trackY.push_back(y);
			}
		}
	}
}
void ConformalMapping() {
	std::vector<Float_t> trackX, trackY;
	std::vector<Float_t> Xc, Yc;

	Int_t counter=0;
	Int_t prev_track_size = 0;
	while(counter < findTrackLimit) {
		FindCenter(trackX, trackY, Xc, Yc);
		++counter;
	}

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
	TH2F *found= new TH2F("found", "found;x;y", 100, -2000, 2000, 100, -2000, 2000);
	for(int j=0; j < tree->GetEntries(); ++j) {
		tree->GetEntry(j);
		if(momentumMag < lowMomLimit) continue; 
		if(momentumMag > highMomLimit) continue; 
		if(particleName->compare("proton") !=0) continue;
		Float_t x = position->Z() - 7650;
		Float_t y = position->Y();
		if(std::find(trackX.begin(), trackX.end(), position->Z() - 7650 ) != trackX.end() && std::find(trackY.begin(), trackY.end(), position->Y() ) != trackY.end()) {
			found->Fill(x, y);
		} else {
			circ->Fill(x, y);
		}
	}

	TEllipse *innerwall = new TEllipse(0, 0, 496); innerwall->SetFillStyle(0); innerwall->SetLineWidth(1);
	TEllipse *outerwall = new TEllipse(0, 0, 840); outerwall->SetFillStyle(0); outerwall->SetLineWidth(1);
	TEllipse *target= new TEllipse(0, 0, 100); target->SetFillStyle(0); target->SetLineWidth(1);

	TCanvas *c = new TCanvas("c", "c", 500, 500);
	h->Draw("colz");
	circ->SetMarkerStyle(7);
	circ->Draw("SAME");
	found->SetMarkerColor(kBlue);
	found->SetMarkerStyle(7);
	found->Draw("SAME");

	innerwall->Draw("SAME");
	outerwall->Draw("SAME");
	target->Draw("SAME");
	for(int z=0; z < Xc.size(); ++z) {
		TEllipse *track = new TEllipse(Xc[z], Yc[z], TMath::Sqrt(Xc[z]*Xc[z] + Yc[z]*Yc[z]), 0);
		track->SetFillStyle(0);
		track->SetLineWidth(1);
		track->SetLineColor(kRed);
		track->SetNoEdges(kTRUE);
		track->Draw("SAME");

		TMarker *mark = new TMarker(Xc[z], Yc[z], kFullCircle);
		mark->SetMarkerColor(kRed);
		mark->Draw("SAME");

		TLatex *numbering = new TLatex(Xc[z], Yc[z], Form("%d", z) );
		numbering->Draw("SAME");
	}
	c->Draw();
	c->SaveAs("results.pdf");

	TCanvas *delta = new TCanvas("delta", "delta");
	hDelta->Draw();
	delta->Draw();
}
