////
Double_t lowMomLimit = .004; //artificial limiters for testing the code
Double_t highMomLimit = .1; //artificial
Int_t findTrackLimit = 20; //limit number of tracks to find
Double_t search_sigma = 12; //TSpectrum sigma
TH1F *hDelta = new TH1F("hDelta", "hDelta; #Delta r^{2} [mm^{2}]", 500, -5e4, 5e4);
TH1F *hEstMom = new TH1F("hEstMom", "hEstMom;p [MeV]", 100, 0, 1000);
TH1I *hHitsPerTrack = new TH1I("hHitsPerTrack", "hHitsPerTrack;HitsPerTrack;Count", 100, 0, 100);
////
Int_t tId = 1;
void FindCenter(
	const char * filename,
	int conformalMappingThreshold, //number will change if you change the histogram binning, this is to remove noise in the u-v space
	std::vector<Int_t> &trackId, 
	std::vector<Float_t> &trackX, 
	std::vector<Float_t> &trackY, 
	std::vector<Float_t> &trackZ, 
	std::vector<Float_t> &vecXc, 
	std::vector<Float_t> &vecYc
) {
	TFile *f = new TFile(filename, "READ");
	TTree *tree = (TTree *)f->Get("cdc");

	Int_t eventId;
	TLorentzVector *position = new TLorentzVector(0., 0., 0., 0.);
	Double_t E;
	tree->SetBranchAddress("pos", &position);
	tree->SetBranchAddress("E", &E);

	TH2F *h = new TH2F("h", "h;u;v", 2000, -2000, 2000, 2000, -2000, 2000);
	for(int j=0; j < tree->GetEntries(); ++j) {
		tree->GetEntry(j);
		//if entry is found in vector, skip
		if(std::find(trackX.begin(), trackX.end(), position->Z() - 7650 ) != trackX.end() && std::find(trackY.begin(), trackY.end(), position->Y() ) != trackY.end()) continue;
		if(E < lowMomLimit) continue; 
		if(E > highMomLimit) continue; 
		Float_t x = position->Z() - 7650;
		Float_t y = position->Y();
		for(Int_t u=-2000; u < 2000; ++u) {
			//Hough transform
			Float_t v = -(x/y)*u + (x*x + y*y)/(2*y);
			//remove points within inner radius, however, this removes low momentum tracks
			//if(u*u + v*v < 496*496 ) continue;
			h->Fill(u, v);
		}
	}
	
	for(int c=0; c < h->GetSize(); ++c) {
		if(h->GetBinContent(c) < conformalMappingThreshold) {
			h->SetBinContent(c, 0);
		}
	}

	TSpectrum *s = new TSpectrum(1);
	TSpectrum *xs = new TSpectrum(1);
	Int_t nfound = s->Search(h->ProjectionY(), search_sigma, "", 0.05);
	Double_t *ypeaks = s->GetPositionX();
	Double_t xp = ypeaks[0];
	TH1D *h_px = h->ProjectionX("h_px", h->ProjectionX()->GetXaxis()->FindBin(xp), h->ProjectionX()->GetXaxis()->FindBin(xp) );
	xs->Search(h_px, search_sigma, "", 0.05);
	Double_t *xpeaks = xs->GetPositionX();
	Int_t Yc = ypeaks[0];
	Int_t Xc = xpeaks[0];
	//display evolution of the conformal mapping projection
	TCanvas *temp = new TCanvas("temp", "temp");
	temp->Divide(2, 2);
	temp->cd(1);
	h->Draw("colz");
	temp->cd(2);
	h_px->Draw();
	temp->cd(3);
	h->ProjectionY()->Draw();
	temp->SaveAs("htemp.pdf");
	delete h_px;
	delete ypeaks;
	delete xpeaks;
	std::cout << "center: (" << Xc << ", " << Yc << ") radius: " << TMath::Sqrt(Xc*Xc + Yc*Yc) << "mapthres: " << conformalMappingThreshold << std::endl;
	int _bHitCount = trackId.size();
	if(Xc!=0 && Yc!=0) {
		hEstMom->Fill(0.3 * TMath::Sqrt(Xc*Xc + Yc*Yc) );
		//push to vector and remove track hits only if nonzero track center is found
		vecXc.push_back(Xc);
		vecYc.push_back(Yc);

		TVector3 v2(Xc, Yc, 0);
		for(int j=0; j < tree->GetEntries(); ++j) {
			tree->GetEntry(j);
			if(E < lowMomLimit) continue; 
			if(E > highMomLimit) continue; 
			Float_t x = position->Z() - 7650;
			Float_t y = position->Y();
			Float_t z = position->X() - 6450;
			TVector3 v1(x, y, 0);
			TVector3 cprod = v1.Cross(v2);
			//check vicinity and direction of curl (positive or negative charged particles)
			hDelta->Fill((x-Xc)*(x-Xc) + (y-Yc)*(y-Yc) - Xc*Xc - Yc*Yc);
		}
		TF1 *fit = new TF1("fit", "gaus", -5e3, 5e3);
		hDelta->Fit("fit", "SR");
		Float_t meanRsquared = fit->GetParameter(1);
		Float_t deltaRsquared = fit->GetParameter(2) * 3;
		for(int j=0; j < tree->GetEntries(); ++j) {
			tree->GetEntry(j);
			if(E < lowMomLimit) continue; 
			if(E > highMomLimit) continue; 
			Float_t x = position->Z() - 7650;
			Float_t y = position->Y();
			Float_t z = position->X() - 6450;
			TVector3 v1(x, y, 0);
			TVector3 cprod = v1.Cross(v2);
			//check vicinity and direction of curl (positive or negative charged particles)
			if(abs((x-Xc)*(x-Xc) + (y-Yc)*(y-Yc) - Xc*Xc - Yc*Yc - meanRsquared) < deltaRsquared) {
				//push points in the vicinity of the circle into the vector to be excluded in the next iteration
				if(Xc*Xc + Yc*Yc < 496*496 )  {
					trackId.push_back(tId);
					trackX.push_back(x);
					trackY.push_back(y);
					trackZ.push_back(z);
				} else {
					if(cprod.Z() > 0) {
						trackId.push_back(tId);
						trackX.push_back(x);
						trackY.push_back(y);
						trackZ.push_back(z);
					}
				}
				
			}
		}
	}
	hHitsPerTrack->Fill(trackId.size() - _bHitCount);
	++tId;
}
void ConformalMapping() {
	std::vector<Int_t> trackId;
	std::vector<Float_t> trackX, trackY, trackZ;
	std::vector<Float_t> Xc, Yc;
	//const char * filename = "m01002507.root";
	const char * filename = "test1.root";

	Int_t counter=0;
	Int_t prev_track_size = 0;
	Float_t prevX, prevY;
	Int_t conformalMappingThreshold = 50;
	while(counter < findTrackLimit) {
//	while(true) {
		FindCenter(filename, conformalMappingThreshold, trackId, trackX, trackY, trackZ, Xc, Yc);
		++counter;
		if(abs(prevX-Xc[Xc.size()-1])<1 && abs(prevY-Yc[Yc.size()-1])<1) {
			conformalMappingThreshold -= 5; 
			if(conformalMappingThreshold < 0) break;
		}
		prevX = Xc[Xc.size()-1];
		prevY = Yc[Yc.size()-1];
	}

	TFile *f = new TFile(filename, "READ"); //8e6 protons
	TTree *tree = (TTree *)f->Get("cdc");

	Int_t eventId;
	TLorentzVector *position = new TLorentzVector(0., 0., 0., 0.);
	Double_t E;
	tree->SetBranchAddress("pos", &position);
	tree->SetBranchAddress("E", &E);

	TH2F *h = new TH2F("h", "h;u;v", 100, -2000, 2000, 100, -2000, 2000);
	TH2F *circ = new TH2F("circ", "circ;x;y", 1000, -8000, 8000, 1000, -8000, 8000);
	TH2F *found= new TH2F("found", "found;x;y", 1000, -8000, 8000, 1000, -8000, 8000);
	for(int j=0; j < tree->GetEntries(); ++j) {
		tree->GetEntry(j);
		if(E < lowMomLimit) continue; 
		if(E > highMomLimit) continue; 
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
	c->SaveAs("results.pdf");

	TCanvas *cTemp= new TCanvas("diag", "diag");
	cTemp->Divide(2, 2);
	cTemp->cd(1);
	hEstMom->Draw();
	cTemp->cd(2);
	hDelta->Draw();
	cTemp->cd(3);
	hHitsPerTrack->Draw();
	cTemp->SaveAs("diagnostics.pdf");

	/////////// output to ROOT file for GENFIT ///////////////

	TFile * fOutput = new TFile("analysis.root", "RECREATE");
	Int_t iId;
	Float_t fX;
	Float_t fY;
	Float_t fZ;
	Float_t fPhi;
	Float_t fTheta;
	Float_t fMag;
	TTree * cdc = new TTree("cdc", "cdc");
	cdc->Branch("trackId", &iId);
	cdc->Branch("x", &fX);
	cdc->Branch("y", &fY);
	cdc->Branch("z", &fZ);
	cdc->Branch("phi", &fPhi);
	cdc->Branch("theta", &fTheta);
	cdc->Branch("mag", &fMag);

	for(Int_t i = 0; i < trackId.size(); ++i) {
		iId = trackId[i];
		fX = trackX[i];
		fY = trackY[i];
		fZ = trackZ[i];
		TVector3 v(fX, fY, fZ);
		fPhi = v.Phi();
		fTheta = v.Theta();
		fMag = 0.3e-3 * TMath::Sqrt(fX*fX + fY*fY); //MeV
		cdc->Fill();
	}
	fOutput->Write();
}
