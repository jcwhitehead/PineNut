// -*- C++ -*-

#include "PineNut.h"

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {



  /// @brief MC validation analysis for Z events
  class MC_ZINC_PN : public Analysis {
  public:

    /// Default constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_ZINC_PN);
    
    /// @name Analysis methods
    //@{

    /// Book histograms
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

      
      _dR=0.2;
      if (getOption("SCHEME") == "BARE")  _dR = 0.0;
		  _lepton=PID::ELECTRON;
      if (getOption("LMODE") == "MU")  _lepton = PID::MUON;

      // set FS cuts from input options
      const double etacut = getOption<double>("ABSETALMAX", 3.5);
      const double ptcut = getOption<double>("PTLMIN", 25.);
      
      FinalState fs;
      Cut cut = Cuts::abseta < etacut && Cuts::pT > ptcut*GeV;

      ZFinder zfinder(fs, cut, _lepton, 66.0*GeV, 116.0*GeV, _dR, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::YES);
      declare(zfinder, "ZFinder");

      PineNut::book(bookFunc, "Z_mass", 50, 66.0, 116.0);
      PineNut::book_logspace(bookFunc ,"Z_pT", 100, 1.0, 0.5*(sqrtS()>0.?sqrtS():14000.)/GeV);
      PineNut::book(bookFunc, "Z_pT_peak", 25, 0.0, 25.0);
      PineNut::book(bookFunc, "Z_y", 40, -4.0, 4.0);
      PineNut::book(bookFunc, "Z_phi", 25, 0.0, TWOPI);
      PineNut::book_logspace(bookFunc ,"lepton_pT", 100, 10.0, 0.25*(sqrtS()>0.?sqrtS():14000.)/GeV);
      PineNut::book(bookFunc, "lepton_eta", 40, -4.0, 4.0);

    }


    /// Do the analysis
    void analyze(const Event & event) {

      // cout << "printing weight names inside Rivet: " << endl;
      // for ( auto & en : event.originalGenEvent()->weight_names() ) {
      // 	cout << en << endl;
      // }
      // cout << "printing Rivet weights inside Rivet: " << endl;
      // for ( auto w : event.weights() ) {
      // 	cout << w << endl;
      // }
      // cout << "printing HepMC weights inside Rivet: " << endl;
      // for ( auto w : event.originalGenEvent()->weights() ) {
      // 	cout << w << endl;
      // }
      
      const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");
      if (zfinder.bosons().size() != 1) vetoEvent;

      FourMomentum zmom(zfinder.bosons()[0].momentum());

      PineNut::fill(event, "Z_mass", zmom.mass()/GeV);
      PineNut::fill(event, "Z_pT", zmom.pT()/GeV);
      PineNut::fill(event, "Z_pT_peak", zmom.pT()/GeV);
      PineNut::fill(event, "Z_y", zmom.rapidity());
      PineNut::fill(event, "Z_phi", zmom.phi());
      for (const Particle& l : zfinder.constituents()) {
	PineNut::fill(event, "lepton_pT", l.pT()/GeV);
	PineNut::fill(event, "lepton_eta", l.eta());
      }
    }


    /// Finalize
    void finalize() {
      const double sf = crossSection()/picobarn/sumOfWeights();

      PineNut::scale(sf);
      PineNut::writegrids();
    }

    //@}


  protected:

    /// @name Parameters for specialised e/mu and dressed/bare subclassing
    //@{
    double _dR;
    PdgId _lepton;
    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_Z_mass;
    Histo1DPtr _h_Z_pT;
    Histo1DPtr _h_Z_pT_peak;
    Histo1DPtr _h_Z_y;
    Histo1DPtr _h_Z_phi;
    Histo1DPtr _h_lepton_pT;
    Histo1DPtr _h_lepton_eta;
    //@}

  };

  // The hooks for the plugin system
  RIVET_DECLARE_PLUGIN(MC_ZINC_PN);
}
