// -*- C++ -*-

#include "PineNut.h"

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief Study of prompt dilepton kinematics sensitive to boson polarizations
  class MC_DILEPTON_PN : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_DILEPTON_PN);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      auto bookFunc = [this](Histo1DPtr & hist, const std::string & histname, const std::vector<double> & bins) {
        book(hist, histname, bins);
      };

      PineNut::init({
	  0, 2, 0, 0,
	  1, 2, 0, 0,
	  1, 2, 1, 0,
	  1, 2, 0, 1,
	}, "NNPDF40MC_nlo_as_01180",
	this);
      
      // set FS cuts from input options
      const double etacut = getOption<double>("ABSETALMAX", 5.);
      const double ptcut = getOption<double>("PTLMIN", 10.);
      
      // Initialise and register projections
      declare(PromptFinalState((Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON)
                               && Cuts::abseta < etacut && Cuts::pT > ptcut*GeV), "Leptons");

      // Book histograms
      PineNut::book_logspace(bookFunc, "lep1_pt", 40, 10, 400);
      PineNut::book(bookFunc, "lep1_costheta", 25, -1, 1);
      PineNut::book(bookFunc, "lep1_ppara", 40, -50, 350);
      PineNut::book(bookFunc, "lep1_pperp", 25, 0, 100);
      //
      PineNut::book_logspace(bookFunc, "lep2_pt", 40, 10, 400);
      PineNut::book(bookFunc, "lep2_costheta", 25, -1, 1);
      PineNut::book(bookFunc, "lep2_ppara", 40, -50, 350);
      PineNut::book(bookFunc, "lep2_pperp", 25, 0, 100);
      //
      PineNut::book(bookFunc, "com_costheta_l1", 25, -1, 1);
      PineNut::book(bookFunc, "com_costheta_l2", 25, -1, 1);
      PineNut::book(bookFunc, "com_ppara_l1", 25, -50, 50);
      PineNut::book(bookFunc, "com_ppara_l2", 25, -50, 50);
      //
      PineNut::book(bookFunc, "com_costheta", 25, -1, 1);
      PineNut::book(bookFunc, "com_ppara", 25, -50, 50);
      PineNut::book(bookFunc, "com_pperp", 25, 0, 100);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const Particles& leptons = apply<FinalState>(event, "Leptons").particlesByPt();
      if (leptons.size() != 2) vetoEvent;
      const Particle& l1 = leptons[0];
      const Particle& l2 = leptons[1];
      
      PineNut::fill(event, "lep1_pt", l1.pT()/GeV);
      PineNut::fill(event, "lep2_pt", l2.pT()/GeV);

      const FourMomentum pcom = l1.mom() + l2.mom();
      const Vector3 betacom = pcom.betaVec();
      const Vector3 unitboostvec = betacom.unit();
      const LorentzTransform comboost = LorentzTransform::mkFrameTransformFromBeta(betacom);

      const double l1_costheta = cos(l1.p3().angle(unitboostvec));
      const double l1_ppara = l2.p3().dot(unitboostvec);
      const double l1_pperp = l2.p3().cross(unitboostvec).mod();

      PineNut::fill(event, "lep1_costheta", l1_costheta);
      PineNut::fill(event, "lep1_ppara", l1_ppara);
      PineNut::fill(event, "lep1_pperp", l1_pperp);

      const double l2_costheta = cos(l2.p3().angle(unitboostvec));
      const double l2_ppara = l2.p3().dot(unitboostvec);
      const double l2_pperp = l2.p3().cross(unitboostvec).mod();

      PineNut::fill(event, "lep2_costheta", l2_costheta);
      PineNut::fill(event, "lep2_ppara", l2_ppara);
      PineNut::fill(event, "lep2_pperp", l2_pperp);

      const FourMomentum p1com = comboost.transform(l1.mom());
      const FourMomentum p2com = comboost.transform(l2.mom());
      const double com_costheta1 = cos(p1com.p3().angle(unitboostvec));
      const double com_costheta2 = cos(p2com.p3().angle(unitboostvec));
      MSG_DEBUG("CoM cos(th)s: " << com_costheta1 << ", " << com_costheta2);
      // assert(com_costheta1 == com_costheta2 && "CoM cos(th)s differ");
      const double com_ppara1 = p1com.p3().dot(unitboostvec);
      const double com_ppara2 = p2com.p3().dot(unitboostvec);
      MSG_DEBUG("CoM p_paras: " << com_ppara1 << ", " << com_ppara2);
      // assert(com_ppara1 == com_ppara2 && "CoM p_paras differ");
      const double com_pperp1 = p1com.p3().cross(unitboostvec).mod();
      const double com_pperp2 = p2com.p3().cross(unitboostvec).mod();
      MSG_DEBUG("CoM p_pperps: " << com_pperp1 << ", " << com_pperp2);
      // assert(com_pperp1 == com_pperp2 && "CoM p_perps differ");

      PineNut::fill(event, "com_costheta_l1", com_costheta1);
      PineNut::fill(event, "com_costheta_l2", com_costheta2);
      PineNut::fill(event, "com_costheta", com_costheta1);
      PineNut::fill(event, "com_costheta", com_costheta2);

      PineNut::fill(event, "com_ppara_l1", com_ppara1);
      PineNut::fill(event, "com_ppara_l2", com_ppara2);
      PineNut::fill(event, "com_ppara", com_ppara1);
      PineNut::fill(event, "com_ppara", com_ppara2);

      PineNut::fill(event, "com_pperp", com_pperp1);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection()/picobarn/sumOfWeights();

      PineNut::scale(sf);

      PineNut::scale("com_costheta", 0.5);
      PineNut::scale("com_ppara", 0.5);
      
      PineNut::writegrids();
    }

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_DILEPTON_PN);


}
