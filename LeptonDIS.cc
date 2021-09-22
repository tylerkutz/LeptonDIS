#include "LeptonDIS.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <HepMC/GenEvent.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>

#include <TFile.h>
#include <TTree.h>

using namespace std;

LeptonDIS::LeptonDIS(const std::string &name) : SubsysReco(name) {

	mOutputFileName = "dis_reconstruction_output.root";

	mElectronE = 10.;
	mHadronE = 100.;

	mCrossingAngle = 25.e-3; // 25 msr

}


LeptonDIS::~LeptonDIS() {

}



void LeptonDIS::SetOutputFile(TString outfilename) {

	mOutputFileName = outfilename;

}

void LeptonDIS::SetBeamEnergies(double electronE, double hadronE) {
	
	mElectronE = electronE;
	mHadronE = hadronE;

}

void LeptonDIS::SetCrossingAngle(double xingAngle) {

	mCrossingAngle = xingAngle;

}

int LeptonDIS::Init(PHCompositeNode* topNode) {

	mEventCounter = 0;

	std::cout << "Setting beam configuration to " << mElectronE << "x" << mHadronE << " GeV" << std::endl;
	m4Ve = TLorentzVector(0.,0.,-sqrt(mElectronE * mElectronE - Me*Me), mElectronE);
	m4Vh = TLorentzVector(0.,0.,sqrt(mHadronE*mHadronE - Mp*Mp), mHadronE);

	mFile = new TFile(mOutputFileName, "RECREATE");
	InitializeTree();

	return Fun4AllReturnCodes::EVENT_OK;

}

int LeptonDIS::process_event(PHCompositeNode* topNode) {

	mEventCounter++;
	if(mEventCounter % 1000 == 0) {
		cout << "Event " << mEventCounter << endl;
	}

	GetNodes(topNode);

	ResetVariables();	

	MCEventInfo();
	LeptonRecon();

	mTree->Fill();

	return Fun4AllReturnCodes::EVENT_OK;

}

int LeptonDIS::End(PHCompositeNode* topNode) {

	mFile->cd();
	mTree->Write("T", TObject::kOverwrite);
	mFile->Close();

	return Fun4AllReturnCodes::EVENT_OK;


}


int LeptonDIS::GetNodes(PHCompositeNode* topNode) {

	mTruth = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");	
	if(!mTruth) {
		cout << PHWHERE << " PHG4TruthInfoContainer node not found!" << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	mTrackMap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
	if(!mTrackMap) {
		cout << PHWHERE << " SvtxTrackMap node not found!" << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	mHepMCMap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
	if(!mHepMCMap) {
		cout << PHWHERE << " MCGenEventMap node not found!" << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}
	

	mVertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
	if(!mVertexMap) {
		cout << PHWHERE << " SvtxVertexMap node not found!" << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	return 1;	

}

void LeptonDIS::MCEventInfo() {

	for (PHHepMCGenEventMap::ConstIter eventIter = mHepMCMap->begin(); eventIter != mHepMCMap->end(); ++eventIter) {
	
		PHHepMCGenEvent* mcEvent = eventIter->second;

		if(mcEvent) {

			HepMC::GenEvent* truthEvent = mcEvent->getEvent();

			// First get truth event info 
			// (should correspond to leptontic truth)
			HepMC::PdfInfo* pdfInfo = truthEvent->pdf_info();
			hepmc_xB = pdfInfo->x2();
			hepmc_Q2 = pdfInfo->scalePDF();	

			// Next loop over truth particles to find virtual photon
			// (should correspond to hadronic/Born truth)
			for (HepMC::GenEvent::particle_const_iterator hepmcPart = truthEvent->particles_begin(); 
				hepmcPart != truthEvent->particles_end(); ++hepmcPart) {

				if ((*hepmcPart)->status() == 3 && (*hepmcPart)->pdg_id() == 23) {
					double pxg = (*hepmcPart)->momentum().px();
					double pyg = (*hepmcPart)->momentum().py();
					double pzg = (*hepmcPart)->momentum().pz();
					double Eg = (*hepmcPart)->momentum().e();
					m4Vg.SetPxPyPzE(pxg, pyg, pzg, Eg);
					CalculateDISKinematics(m4Ve, m4Ve - m4Vg, m4Vh, born_xB, born_Q2, born_W2, born_y, born_eta);				
				}
				
			} // HepMC particle loop
			
		}

	} // EventMap loop

}

void LeptonDIS::LeptonRecon() {

	PHG4TruthInfoContainer::ConstRange range = mTruth->GetPrimaryParticleRange();

	for(PHG4TruthInfoContainer::ConstIterator truthItr = range.first; truthItr != range.second; ++truthItr) {
		
		PHG4Particle* truthParticle = truthItr->second;

		// Reconstructing with electron
		// Only care about matching to truth electron track
		if(truthParticle->get_pid() != 11 || truthParticle->get_parent_id() != 0) continue;
	
		SvtxTrack* recTrack = nullptr;
		
		for(SvtxTrackMap::ConstIter trackItr = mTrackMap->begin(); trackItr != mTrackMap->end(); trackItr++) {

			SvtxTrack* tempTrack = dynamic_cast<SvtxTrack*>(trackItr->second);

			if(!tempTrack) {
				continue;
			}	
			
			if( (tempTrack->get_truth_track_id() - truthParticle->get_track_id()) == 0) {
				recTrack = tempTrack;
			}			


		} // TrackMap iterator

		g4_tr_px = truthParticle->get_px();
		g4_tr_py = truthParticle->get_py();
		g4_tr_pz = truthParticle->get_pz();
		g4_tr_p  = sqrt(g4_tr_px*g4_tr_px + g4_tr_py*g4_tr_py + g4_tr_pz*g4_tr_pz);
		g4_tr_th = acos(g4_tr_pz / g4_tr_p);
		g4_tr_ph = atan2(g4_tr_py, g4_tr_px);		

		TLorentzVector g4_eprime = CorrectCrossingAngle(TLorentzVector(g4_tr_px, g4_tr_py, g4_tr_pz, sqrt(g4_tr_p*g4_tr_p + Me*Me)));

		CalculateDISKinematics(m4Ve, g4_eprime, m4Vh, g4_tr_xB, g4_tr_Q2, g4_tr_W2, g4_tr_y, g4_tr_eta);	

		if(recTrack) {

			mRecTrack = 1;
			rec_tr_px = recTrack->get_px();
			rec_tr_py = recTrack->get_py();
			rec_tr_pz = recTrack->get_pz();
			rec_tr_p  = sqrt(rec_tr_px*rec_tr_px + rec_tr_py*rec_tr_py + rec_tr_pz*rec_tr_pz);
			rec_tr_th = acos(rec_tr_pz / rec_tr_p);
			rec_tr_ph = atan2(rec_tr_py, rec_tr_px);		
		
			TLorentzVector rec_eprime = CorrectCrossingAngle(TLorentzVector(rec_tr_px, rec_tr_py, rec_tr_pz, sqrt(rec_tr_p*rec_tr_p + Me*Me)));

			CalculateDISKinematics(m4Ve, rec_eprime, m4Vh, rec_tr_xB, rec_tr_Q2, rec_tr_W2, rec_tr_y, rec_tr_eta);	

		}

	} // TruthInfo iterator

		
}

void LeptonDIS::InitializeTree() {

	mTree = new TTree("T", "Tree with DIS kinematic variables");

	mTree->Branch("event",		&mEventCounter,	"event/I");

	mTree->Branch("rec_track",	&mRecTrack,	"rec_track/I");

	mTree->Branch("hepmc_xB",	&hepmc_xB,	"hepmc_xB/D");
	mTree->Branch("hepmc_Q2",	&hepmc_Q2,	"hepmc_Q2/D");

	mTree->Branch("born_y",		&born_y,	"born_y/D");
	mTree->Branch("born_Q2",	&born_Q2,	"born_Q2/D");
	mTree->Branch("born_xB",	&born_xB,	"born_xB/D");
	mTree->Branch("born_W2",	&born_W2,	"born_W2/D");
	mTree->Branch("born_eta",	&born_eta,	"born_eta/D");

	mTree->Branch("rec_tr_px",	&rec_tr_px,	"rec_tr_px/D");	
	mTree->Branch("rec_tr_py",	&rec_tr_py,	"rec_tr_py/D");	
	mTree->Branch("rec_tr_pz",	&rec_tr_pz,	"rec_tr_pz/D");	
	mTree->Branch("rec_tr_p",	&rec_tr_p,	"rec_tr_p/D");	
	mTree->Branch("rec_tr_th",	&rec_tr_th,	"rec_tr_th/D");	
	mTree->Branch("rec_tr_ph",	&rec_tr_ph,	"rec_tr_ph/D");	

	mTree->Branch("g4_tr_px",	&g4_tr_px,	"g4_tr_px/D");	
	mTree->Branch("g4_tr_py",	&g4_tr_py,	"g4_tr_py/D");	
	mTree->Branch("g4_tr_pz",	&g4_tr_pz,	"g4_tr_pz/D");	
	mTree->Branch("g4_tr_p",	&g4_tr_p,	"g4_tr_p/D");	
	mTree->Branch("g4_tr_th",	&g4_tr_th,	"g4_tr_th/D");		
	mTree->Branch("g4_tr_ph",	&g4_tr_ph,	"g4_tr_ph/D");	

	mTree->Branch("rec_tr_y",	&rec_tr_y,	"rec_tr_y/D");
	mTree->Branch("rec_tr_Q2",	&rec_tr_Q2,	"rec_tr_Q2/D");
	mTree->Branch("rec_tr_xB",	&rec_tr_xB,	"rec_tr_xB/D");
	mTree->Branch("rec_tr_W2",	&rec_tr_W2,	"rec_tr_W2/D");
	mTree->Branch("rec_tr_eta",	&rec_tr_eta,	"rec_tr_eta/D");

	mTree->Branch("g4_tr_y",	&g4_tr_y,	"g4_tr_y/D");
	mTree->Branch("g4_tr_Q2",	&g4_tr_Q2,	"g4_tr_Q2/D");
	mTree->Branch("g4_tr_xB",	&g4_tr_xB,	"g4_tr_xB/D");
	mTree->Branch("g4_tr_W2",	&g4_tr_W2,	"g4_tr_W2/D");
	mTree->Branch("g4_tr_eta",	&g4_tr_eta,	"g4_tr_eta/D");

}

void LeptonDIS::ResetVariables() {
 
	mRecTrack = 0;

	hepmc_xB = -99.;
	hepmc_Q2 = -99.;

	born_y = -99.; 
	born_Q2 = -99.; 
	born_xB = -99.; 
	born_W2 = -99.; 
	born_eta = -99.;

	rec_tr_px = -99.;
	rec_tr_py = -99.;
	rec_tr_pz = -99.;
	rec_tr_p = -99.;
	rec_tr_th = -99.;
	rec_tr_ph = -99.;
		
	g4_tr_px = -99.;
	g4_tr_py = -99.;
	g4_tr_p  = -99.;
	g4_tr_pz = -99.;
	g4_tr_th = -99.;
	g4_tr_ph = -99.;

	rec_tr_y = -99.;
	rec_tr_Q2 = -99.;
	rec_tr_xB = -99.;
	rec_tr_W2 = -99.;
	rec_tr_eta = -99.;
		  
	g4_tr_y = -99.;
	g4_tr_Q2 = -99.;
	g4_tr_xB = -99.;
	g4_tr_W2 = -99.;
	g4_tr_eta = -99.;	

}

TLorentzVector LeptonDIS::CorrectCrossingAngle(TLorentzVector vecIn) {

	TLorentzVector vecOut = vecIn;
	vecOut.RotateY(mCrossingAngle/2.);
	vecOut.Boost(sin(mCrossingAngle/2.),0,0);
	return vecOut;

}

void LeptonDIS::CalculateDISKinematics(TLorentzVector l, TLorentzVector lp, TLorentzVector P, double& xB, double& Q2, double& W2, double& y, double& eta) {

	TLorentzVector q = l - lp;
	Q2 = -q*q;
	xB = Q2/(2. * (P*q));
	W2 = (P+q)*(P+q);
	y = (q*P)/(l*P);
	eta = lp.PseudoRapidity();

}


