#include<cmath>

#include "pineappl_capi.h"
#include "PineAPPL.hpp"

#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"

#ifndef RIVET_ENABLE_HEPMC_3
#include "HepMC/HepMCDefs.h"
#include "HepMC/GenEvent.h"
#include "HepMC/Attribute.h"
#endif

#include "LHAPDF/LHAPDF.h"

// PineNut:
// A header-only tool to make pineAPPL grids with Rivet analyses.
// Tested with:
//   - Rivet v3.1.x
//   - pineAPPL v0.8.x

namespace PineNut {
  
  // LHAPDF PDF object matching PDF from run (optional?)
  LHAPDF::PDF* runpdf;

  // Rivet histograms
  std::map<string, Rivet::Histo1DPtr> _h;
  // PineAPPL grids
  std::map<string, pineappl_grid*> gridPointers;

  // PineAPPL variables
  std::vector<uint32_t> orders;
  std::map<std::pair<int,int>,int> pineappl_channelids;
  std::string pineappl_prefix = "pineappl";

  // pointer to Rivet analysis
  Rivet::Analysis * analysispointer;
    
    
  // Initialise PineNut
  // (call inside Rivet init())
  void init(std::vector<uint32_t> orders_in,
	    string pdfname,
	    Rivet::Analysis * rana) {
    orders = orders_in;
    runpdf = LHAPDF::mkPDF(pdfname, 0);
    analysispointer = rana;
  };


  // Book histograms within Rivet and PineAPPL
  // (call within init() for each histogram)
  void book(std::function<void(Rivet::Histo1DPtr & hist, const std::string &, const std::vector<double> &)> bookFunc,
	    const std::string & histname,
	    const std::vector<double> & histbins) {
    // book Rivet histogram
    bookFunc(_h[histname], histname, histbins);
    
    // book PineAppl grid
    pineappl_lumi * channels = pineappl_lumi_new();
    
    // for now we use separate grids for each non-zero partonic channel
    // [NB: we are not generating a photon-initiated channel here]
    int cid = 0;
    for ( int a = -5; a < 6; a++ ) {
      for ( int b = -5; b < 6; b++ ) {
	int32_t pids[] = {a, b};
	pineappl_lumi_add(channels, 1, pids, nullptr);
	pineappl_channelids[{a,b}] = cid++;
      }
    }
    
    pineappl_keyval* keyval = pineappl_keyval_new();
    pineappl_grid* grid = pineappl_grid_new(channels,
					    orders.size() / 4,
					    orders.data(),
					    histbins.size() - 1,
					    histbins.data(),
					    keyval);

    // pineappl_grid_set_key_value(grid, "arxiv", "NONE");
    pineappl_grid_set_key_value(grid, "inspire", analysispointer->info().inspireId().c_str());

    pineappl_grid_set_key_value(grid, "rivet_description", analysispointer->info().description().c_str());
    pineappl_grid_set_key_value(grid, "rivet_name", analysispointer->name().c_str());
    // pineappl_grid_set_key_value(grid, "rivet_path", analysispointer->histoPath(histname).c_str());

    // could also set PDF ids, hadron IDs
    
    gridPointers[histname] = grid;
    
    pineappl_keyval_delete(keyval);
    pineappl_lumi_delete(channels);
  };


  // Fill histograms within Rivet and PineAPPL simultaneously, with a single call
  // (call within analyze() for each histogram fill)
  void fill(const Rivet::Event & event,
	    string histname,
	    double histvalue) {

    int pinepid1 = 0, pinepid2 = 0;
    double x1 = 0.0, x2 = 0.0;
    double xf1 = 0.0, xf2 = 0.0;
    double muf = 0.0, mur = 0.0;
    double alphaqcd = 0.0, alphaqed = 0.0;

    std::cout << "trying to fill " << histname << std::endl;
    
    const Rivet::GenEvent* evt = event.genEvent();
    if ( evt->pdf_info()->pdf_id[0] != 0
      && evt->pdf_info()->pdf_id[1] != 0 ) {
      
      pinepid1 = evt->pdf_info()->parton_id[0];
      pinepid2 = evt->pdf_info()->parton_id[1];
      
      if ( pinepid1 == 21 ) pinepid1 = 0;
      if ( pinepid2 == 21 ) pinepid2 = 0;	

      x1 = evt->pdf_info()->x[0];
      x2 = evt->pdf_info()->x[1];

      xf1 = evt->pdf_info()->xf[0];
      xf2 = evt->pdf_info()->xf[1];

      muf = evt->pdf_info()->scale;
    } else {
      std::cout << "PDF info not found (throw error)" << std::endl;
    }
    
    mur = evt->attribute<HepMC3::DoubleAttribute>("event_scale")->value();
    alphaqcd = evt->attribute<HepMC3::DoubleAttribute>("alphaQCD")->value();
    alphaqed = evt->attribute<HepMC3::DoubleAttribute>("alphaQED")->value();

    std::cout << "scales : muf = " << muf << " and mur = " << mur << std::endl;
    double alphaqcd_lhapdf = runpdf->alphasQ(mur);
    std::cout << "alphas from lhapdf = " << alphaqcd_lhapdf
	 << "; alphas from hepmc = " << alphaqcd << std::endl;

    double pineappl_weight = ( xf1 * xf2 == 0.0 ) ? 0.0 : event.weights()[0] * (x1 * x2) / (pow(alphaqcd, orders[0]) * xf1 * xf2);

    std::cout << "about to fill rivet " << histname << std::endl;
    
    // not sure how to do multiple weights...
    _h[histname]->fill(histvalue);

    std::cout << "filled rivet " << histname << " with " << histvalue << std::endl;
    
    std::cout << "pineids: " << pinepid1 << " , " << pinepid2 << " ; channel ID = " << pineappl_channelids[{pinepid1, pinepid2}] << std::endl;
    std::cout << "calling pineappl_grid_fill with "
         << " x1 = " << x1 << " , "
         << " x2 = " << x2 << " , "
         << " muf = " << muf << " , "
         << " orders = " << orders[0] << " , "
         << " histvalue = " << histvalue << " , "
         << " weight = " << pineappl_weight << std::endl;
    pineappl_grid_fill(gridPointers[histname],
		       x1, x2,
		       muf,
		       0,
		       histvalue,
		       pineappl_channelids[{pinepid1, pinepid2}],
		       pineappl_weight);
    std::cout << "filled papl " << histname << " with " << histvalue << std::endl;

  };


  // Scale Rivet histogram and PineAPPL grids simultaneously, with a single call
  // (call within finalize())
  void scale(string histname,
	     double scalefactor) {
    analysispointer->scale(_h[histname], scalefactor);
    pineappl_grid_scale(gridPointers[histname], scalefactor);
  };


  // Write PineAPPL grids to file
  // (call within finalize())
  void writegrids( ) {
    for ( auto const & hgrid : gridPointers ){
      if ( hgrid.second != nullptr ) {
	string pineappl_filename = pineappl_prefix + hgrid.first + "_PineNuttesting.lz4";
	pineappl_grid_write(hgrid.second, pineappl_filename.c_str());
      }
    }
  } 
}
