#ifndef LEPTONDIS_H__
#define LEPTONDIS_H__

#include <fun4all/SubsysReco.h>

#include <TLorentzVector.h>
#include <TString.h>

// ROOT classes
class TFile;
class TTree;

// fun4all classes
class SvtxTrackMap;
class SvtxVertexMap;
class GlobalVertexMap;
class PHG4TruthInfoContainer;
class PHHepMCGenEventMap;

/// Definition of this analysis module class
class LeptonDIS : public SubsysReco {

public:

	/// Constructor
	LeptonDIS(const std::string &name = "LeptonDIS");

	// Destructor
	virtual ~LeptonDIS();

	/// SubsysReco initialize processing method
	int Init(PHCompositeNode *);

	/// SubsysReco event processing method
	int process_event(PHCompositeNode *);

	/// SubsysReco end processing method
	int End(PHCompositeNode *);

	void SetBeamEnergies(double, double);
	void SetOutputFile(TString);
	void SetCrossingAngle(double);

private:

	double Me = 0.000510999;	// GeV
	double Mp = 0.938272088;  	// GeV

	double mElectronE;
	double mHadronE;
	double mCrossingAngle;
	
	TLorentzVector m4Ve;	// electron beam
	TLorentzVector m4Vh;	// hadron beam
	TLorentzVector m4Vg;	// virtual photon

	int GetNodes(PHCompositeNode *);

	void MCEventInfo();
	void LeptonRecon();

	TLorentzVector CorrectCrossingAngle(TLorentzVector);

	SvtxTrackMap* mTrackMap;
	SvtxVertexMap* mVertexMap;
	PHG4TruthInfoContainer* mTruth;	

	PHHepMCGenEventMap* mHepMCMap;

	TString mOutputFileName;
	TFile* mFile;
	TTree* mTree;

	void InitializeTree();
	void ResetVariables();

	int mEventCounter;
	int mRecTrack;

	double hepmc_xB;
	double hepmc_Q2;

	double born_y;
	double born_Q2;
	double born_xB;
	double born_W2;
	double born_eta;

	double rec_tr_px;	
	double rec_tr_py;	
	double rec_tr_pz;	
	double rec_tr_p;	
	double rec_tr_th;
	double rec_tr_ph;
	
	double g4_tr_px;	
	double g4_tr_py;	
	double g4_tr_pz;	
	double g4_tr_p;	
	double g4_tr_th;	
	double g4_tr_ph;	
		      
	double rec_tr_y;
	double rec_tr_Q2;
	double rec_tr_xB;
	double rec_tr_W2;
	double rec_tr_eta;

	double g4_tr_y;
	double g4_tr_Q2;
	double g4_tr_xB;
	double g4_tr_W2;
	double g4_tr_eta;

	void CalculateDISKinematics(TLorentzVector, TLorentzVector, TLorentzVector, double&, double&, double&, double&, double&);

};

#endif
