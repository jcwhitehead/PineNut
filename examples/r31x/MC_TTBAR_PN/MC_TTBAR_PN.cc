#include "PineNut.h"

// #include "Rivet/Analysis.hh"
// #include "Rivet/AnalysisHandler.hh"

// #ifndef RIVET_ENABLE_HEPMC_3
// #include "HepMC/HepMCDefs.h"
// #include "HepMC/GenEvent.h"
// #include "HepMC/Attribute.h"
// #endif

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/AnalysisLoader.hh"


namespace Rivet {

  class MC_TTBAR_PN : public Analysis {
    
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_TTBAR_PN);

    // template<typename... Args>
    // void mybook(Args&& ... args) {
    //   book(std::forward<Args>(args) ...);
    // }
    
    std::map<std::string, std::vector<double> > bins;

    void init() {      
      auto bookFunc = [this](Histo1DPtr & hist, const std::string & histname, const std::vector<double> & bins) {
        book(hist, histname, bins);
      };

      // format:
      // order in alphas, order in alphaEM, power of log(muF²), power of log(muR²)
      PineNut::init({
	  2, 0, 0, 0,
	  3, 0, 0, 0,
	  3, 0, 1, 0,
	  3, 0, 0, 1,
	}, "NNPDF40MC_nlo_as_01180",
	this);
      
      FinalState fs;

      bins["tt_absrap"] = {0.0, 0.3, 0.6, 0.9, 1.3, 2.5};
      bins["tt_m"]      = {345, 400, 470, 550, 650, 800, 1100, 1600};
      bins["t_absrap"]  = {0, 0.4, 0.8, 1.2, 1.6, 2.5};
      bins["tt_m_vs_tt_absrap"] = {345, 400, 470, 550, 650, 800, 1100, 1600,      //0-0.3
	1945, 2000, 2070, 2150, 2250, 2400, 2700, 3200, //0.3-0.6
	3545, 3600, 3670, 3750, 3850, 4000, 4300, 4800, //0.6-0.9
	5145, 5200, 5270, 5350, 5450, 5600, 5900, 6400, //0.9-1.3
	6745, 6800, 6870, 6950, 7050, 7200, 7500, 8000}; //1.3-2.5
      bins["t_pt_vs_t_absrap"] = {0,    60,   100,  150,  200,  260,  320,  400,  500, 
	560,   600,  650,  700,  760,  820,  900, 1000, 
	1060,  1100, 1150, 1200, 1260, 1320, 1400, 1500, 
	1560,  1600, 1650, 1700, 1760, 1820, 1900, 2000, 
	2060,  2100, 2150, 2200, 2260, 2320, 2400, 2500};
      bins["tt_pt"] = {0, 60, 100, 150, 200, 260, 320, 400, 500};
      bins["tt_dR"] = {0.  , 0.25, 0.5 , 0.75, 1.  , 1.25, 1.5 , 1.75, 2.  , 2.25, 2.5 ,
	2.75, 3.  , 3.25, 3.5 , 3.75, 4.  , 4.25, 4.5 , 4.75, 5.  , 5.25,
	5.5 , 5.75, 6.  , 6.25, 6.5 , 6.75, 7.  , 7.25, 7.5 , 7.75, 8.  ,
	8.25, 8.5 , 8.75, 9.  , 9.25, 9.5 , 9.75, 10.};

      bins["t_pt"]      = {0, 60, 100, 150, 200, 260, 320, 400, 500};
      bins["t_phi"] = {0.    , 0.1565, 0.313 , 0.4695, 0.626 , 0.7825, 0.939 , 1.0955,
	1.252 , 1.4085, 1.565 , 1.7215, 1.878 , 2.0345, 2.191 , 2.3475,
	2.504 , 2.6605, 2.817 , 2.9735, 3.13  , 3.2865, 3.443 , 3.5995,
	3.756 , 3.9125, 4.069 , 4.2255, 4.382 , 4.5385, 4.695 , 4.8515,
	5.008 , 5.1645, 5.321 , 5.4775, 5.634 , 5.7905, 5.947 , 6.1035,
	6.26 , 6.28, 6.30, 6.32, 6.34, 6.36, 6.38, 6.4};

      bins["top_pt"] = bins["t_pt"];
      bins["top_rap"] = {-9, -3.5, -2.5, -1.3, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.3, 2.5, 3.5, 9};    
      bins["top_phi"] = bins["t_phi"];
      bins["top_e"] = {0, 50, 100, 200, 300, 400, 500, 600, 700, 800, 1100, 1600, 4000};
      
      bins["antitop_pt"] = bins["t_pt"];
      bins["antitop_rap"] = bins["top_rap"];
      bins["antitop_phi"] = bins["t_phi"];
      bins["antitop_e"] = bins["top_e"];

      bins["XS"] = {0, 1};
      bins["pmXS"] = {-1, 0, 1};

      for ( auto hist : bins) {
	cout << "booking " << hist.first << endl;
	cout.flush();
	PineNut::book(bookFunc, hist.first, hist.second);
      };
    }

    void analyze(const Event & event) {
      PineNut::fill(event, "pmXS", 0.5 * event.weights()[0] / abs(event.weights()[0]));
      PineNut::fill(event, "XS", 0.5);
      
      VetoedFinalState fs(FinalState(Cuts::pT > 0.*GeV));
      fs.project(event);
      vector<PseudoJet> particles;
      particles.reserve(fs.particles().size());
      
      bool FoundTop = false;
      bool FoundAntiTop = false;
      FourMomentum r_top;
      FourMomentum r_antitop;

      // truth-level top identification
      for ( long unsigned int i = 0; i < fs.particles().size(); i++ ) {
	const Particle &p = fs.particles().at(i);
	if ( p.pid() == 6 ){
	  r_top = p.momentum();
	  FoundTop = true;
	}
	if ( p.pid() == -6 ){
	  r_antitop = p.momentum();
	  FoundAntiTop = true;
	}
      }

      if ( !FoundTop || !FoundAntiTop ) return;
      
      FourMomentum r_ttbar = r_top + r_antitop;

      if ( ( abs(r_ttbar.rapidity()) < 2.5 )
	   && ( r_ttbar.mass() < 1600 ) ) {

	PineNut::fill(event, "tt_absrap", abs(r_ttbar.rapidity()));
	PineNut::fill(event, "tt_m", abs(r_ttbar.mass()));

	double shift = 0.0;
	if (abs(r_ttbar.rapidity()) < 0.3) {
	  shift = 0;
	} else if (abs(r_ttbar.rapidity()) < 0.6) {
	  shift = 1600;
	} else if (abs(r_ttbar.rapidity()) < 0.9) {
	  shift = 3200;
	} else if (abs(r_ttbar.rapidity()) < 1.3) {
	  shift = 4800;
	} else if (abs(r_ttbar.rapidity()) < 2.5) {
	  shift = 6400;
	}
	PineNut::fill(event, "tt_m_vs_tt_absrap", r_ttbar.mass() + shift);
      }

      std::map<std::string, double> obs;

      obs["tt_pt"] = r_ttbar.pt();
      obs["tt_dR"] = deltaR(r_top, r_antitop);
      
      obs["top_pt"] = r_top.pt();
      obs["top_rap"] = r_top.rapidity();
      obs["top_phi"] = r_top.phi();
      obs["top_e"] = r_top.E();
      
      obs["antitop_pt"] = r_antitop.pt();
      obs["antitop_rap"] = r_antitop.rapidity();
      obs["antitop_phi"] = r_antitop.phi();
      obs["antitop_e"] = r_antitop.E();

      if ( (abs(r_top.rapidity()) < 2.5) && (r_top.pt() < 500)
	   && (abs(r_antitop.rapidity()) < 2.5) && (r_antitop.pt() < 500) ) {

	double shift = 0.0;
	if (abs(r_top.rapidity()) < 0.4) {
	  shift = 0;
	} else if (abs(r_top.rapidity()) < 0.8) {
	  shift = 500;
	} else if (abs(r_top.rapidity()) < 1.2) {
	  shift = 1000;
	} else if (abs(r_top.rapidity()) < 1.6) {
	  shift = 1500;
	} else if (abs(r_top.rapidity()) < 2.5) {
	  shift = 2000;
	}
	obs["t_pt_vs_t_absrap"] = r_top.pt() + shift;

	shift = 0.0;
	// For the anti-top quark, repeat the same logic
	if (abs(r_antitop.rapidity()) < 0.4) {
	  shift = 0 ;
	} else if (abs(r_antitop.rapidity()) < 0.8) {
	  shift = 500;
	} else if (abs(r_antitop.rapidity()) < 1.2) {
	  shift = 1000;
	} else if (abs(r_antitop.rapidity()) < 1.6) {
	  shift = 1500;
	} else if (abs(r_antitop.rapidity()) < 2.5) {
	  shift = 2000;
	}
	obs["t_pt_vs_t_absrap"] = r_antitop.pt() + shift;

	for ( auto ob : obs ) {
	  PineNut::fill(event, ob.first, ob.second);
	}

	PineNut::fill(event, "t_absrap", abs(r_top.rapidity()));
	PineNut::fill(event, "t_absrap", abs(r_antitop.rapidity()));

	PineNut::fill(event, "t_pt", r_top.pt());
	PineNut::fill(event, "t_pt", r_antitop.pt());
      }
    }

    void finalize() {
      const double sf = crossSection() / sumOfWeights();

      PineNut::scale(sf);

      PineNut::scale("t_absrap", 0.5);
      PineNut::scale("t_pt", 0.5);
      PineNut::scale("t_pt_vs_t_absrap", 0.5);
      
      PineNut::writegrids();
    }

  };

  RIVET_DECLARE_PLUGIN(MC_TTBAR_PN);

}
