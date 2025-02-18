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

  // weights less than this value will be ignored
  const double zerotol = 1.e-20;
    
  // Initialise PineNut
  // (call inside Rivet init())
  void init(std::vector<uint32_t> orders_in,
	    string pdfname,
	    Rivet::Analysis * rana) {
    orders = orders_in;
    runpdf = LHAPDF::mkPDF(pdfname, 0);
    analysispointer = rana;
  };
  
  std::vector<double> linspace_bins(double start, double stop, int nbins) {
    if (nbins <= 0) {
      throw std::invalid_argument("Number of bins must be positive");
    }
    if ( stop <= start ) {
      throw std::invalid_argument("Upper limit must exceed lower limit");
    }
    std::vector<double> result(nbins + 1, 0.0);
    double step = (stop - start) / nbins;
     for (int i = 0; i <= nbins; ++i) {
        result[i] = start + i * step;
    }
    return result;
  }

  std::vector<double> logspace_bins(double start, double stop, int nbins) {
    if (nbins <= 0) {
      throw std::invalid_argument("Number of bins must be positive");
    }
    if ( stop <= start ) {
      throw std::invalid_argument("Upper limit must exceed lower limit");
    }
    std::vector<double> result(nbins + 1, 0.0);
    double step = (stop - start) / nbins;
     for (int i = 0; i <= nbins; ++i) {
       result[i] = std::pow(10.0, start + i * step);
    }
    return result;
  }
  
  std::vector<int> extractOrders(const std::string & input) {
    std::vector<int> strorders(4, 0);
    std::vector<std::string> keys = {"as", "ae", "lR", "lF"};
    for (size_t i = 0; i < keys.size(); ++i) {
      size_t pos = input.find(keys[i]);
      if (pos != std::string::npos) {
	pos += keys[i].length();
	size_t end = input.find('_', pos);
	std::string numStr = input.substr(pos, end - pos);
	strorders[i] = std::stoi(numStr);
      }
    }
    return strorders;
  }
  
  int lookUpOrder(const std::vector<int> & order) {
    if ( order.size() != 4 ) {
      throw std::invalid_argument("PineAPPL orders must contain four integers.");
    }
    for ( size_t i = 0; i < orders.size() / 4 ; ++i ) {
      bool match = true;
      for ( size_t j = 0; j < 4; ++j ) {
	if ( order[j] != static_cast<int>(orders[4*i + j]) ) {
	  match = false;
	  break;
	}
      }
      if ( match ) return i;
    }
    std::stringstream orderstr;
    std::copy(order.begin(), order.end(), std::ostream_iterator<int>(orderstr, ","));
    throw std::invalid_argument("No such order registered with PineAPPL. Please consider adding " + orderstr.str() + " to init()");
  }

  
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

  void book(std::function<void(Rivet::Histo1DPtr & hist, const std::string &, const std::vector<double> &)> bookFunc,
	    const std::string & histname,
	    const int nbins,
	    const double xmin,
	    const double xmax) {
    std::vector<double> histbins = linspace_bins(xmin, xmax, nbins);
    book(bookFunc, histname, histbins);
  }

  void book_logspace(std::function<void(Rivet::Histo1DPtr & hist, const std::string &, const std::vector<double> &)> bookFunc,
	    const std::string & histname,
	    const int nbins,
	    const double xmin,
	    const double xmax) {
    if ( xmin <= 0.0 ) {
      throw std::invalid_argument("Lower edge of histogram must be positive for log-scale (not " + std::to_string(xmin) + ")");
    }
    if ( xmin <= 0.0 ) {
      throw std::invalid_argument("Upper edge of histogram must be positive for log-scale (not " + std::to_string(xmax) + ")");
    }
    std::vector<double> histbins = logspace_bins(log10(xmin), log10(xmax), nbins);
    book(bookFunc, histname, histbins);
  }
  

  // Fill histograms within Rivet and PineAPPL simultaneously, with a single call
  // (call within analyze() for each histogram fill)
  void fill(const Rivet::Event & event,
	    string histname,
	    double histvalue) {

    std::cout << "about to fill rivet " << histname << std::endl;    
    _h[histname]->fill(histvalue);
    std::cout << "filled rivet " << histname << " with " << histvalue << std::endl;
    
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
    std::cout << "alphaem from hepmc = " << alphaqed << std::endl;
    std::cout << "alphas from lhapdf = " << alphaqcd_lhapdf
	 << "; alphas from hepmc = " << alphaqcd << std::endl;

    if ( alphaqcd < 0 ) {
      std::cout << "using alphas from lhapdf" << std::endl;
      alphaqcd = alphaqcd_lhapdf;
    }

    // std::cout << "pdf values : " << xf1 << " , " << xf2 << " vs from LHAPDF: " << runpdf->xfxQ(event.genEvent()->pdf_info()->parton_id[0], x1, muf) << " , " << runpdf->xfxQ(event.genEvent()->pdf_info()->parton_id[1], x2, muf) << std::endl;

    // double xf1_lhapdf = runpdf->xfxQ(event.genEvent()->pdf_info()->parton_id[0], x1, muf);
    // double xf2_lhapdf = runpdf->xfxQ(event.genEvent()->pdf_info()->parton_id[1], x2, muf);
    
    // if ( ( abs(xf1 / xf1_lhapdf) - 1 ) > 1e-8 ) {
    //   std::cout << "WARNING PDF 1 differs : " << xf1 << " vs from LHAPDF: " << xf1_lhapdf << std::endl;  
    // }
    // if ( ( abs(xf2 / xf2_lhapdf) - 1 ) > 1e-8 ) {
    //   std::cout << "WARNING PDF 2 differs : " << xf2 << " vs from LHAPDF: " << xf2_lhapdf << std::endl;  
    // }
    
    double pineappl_weight = ( xf1 * xf2 == 0.0 ) ? 0.0 : event.originalGenEvent()->weights()[0] * (x1 * x2) / (xf1 * xf2);

    size_t i_weight = 0;
    for ( auto & wname : event.originalGenEvent()->weight_names() ) {
      std::cout << "considering weight: " << i_weight << ", " << wname << " : " << event.originalGenEvent()->weights()[i_weight] << std::endl;

      // skip nominal and numerically very small weights
      if ( i_weight == 0
	   || abs(event.originalGenEvent()->weights()[i_weight]) < zerotol ) {
	i_weight += 1;
	continue;
      }
      
      // possibly Herwig-specific? perhaps try to do something smarter than this
      std::vector<int> evtorders = extractOrders(wname);
      uint32_t evtorder_idx = lookUpOrder(evtorders);

      // now we have to work out the contribution type
      std::cout << "pineids: " << pinepid1 << " , " << pinepid2 << " ; channel ID = " << pineappl_channelids[{pinepid1, pinepid2}] << std::endl;

      std::cout << "calling pineappl_grid_fill with "
		<< " x1 = " << x1 << " , "
		<< " x2 = " << x2 << " , "
		<< " muf = " << muf << " , "
		<< " evtorder_idx = " << evtorder_idx
		<< " orders = " << evtorders[0] << " , " << evtorders[1] << " , " << evtorders[2] << " , "<< evtorders[3]
		<< " histvalue = " << histvalue << " , "
		<< " weight = " << pineappl_weight * event.originalGenEvent()->weights()[i_weight] / (pow(alphaqcd, evtorders[0])) << std::endl;
      pineappl_grid_fill(gridPointers[histname],
			 x1, x2,
			 muf*muf,
			 evtorder_idx,
			 histvalue,
			 pineappl_channelids[{pinepid1, pinepid2}],
			 pineappl_weight * event.originalGenEvent()->weights()[i_weight] / pow(alphaqcd, evtorders[0]));
      std::cout << "filled papl " << histname << " with " << histvalue << std::endl;
      i_weight += 1;
    }
  };


  // Scale Rivet histogram and PineAPPL grids simultaneously, with a single call
  // (call within finalize())
  void scale(string const histname,
	     double const scalefactor) {
    std::cout << "scaling by scalefactor " << scalefactor << std::endl;
    analysispointer->scale(_h[histname], scalefactor);
    pineappl_grid_scale(gridPointers[histname], scalefactor);
  };


  void scale(double const scalefactor) {
    std::cout << "scaling " << analysispointer->name().c_str()
	      << " by scalefactor " << scalefactor << std::endl;
    for ( auto const & hist : _h ) {
      analysispointer->scale(hist.second, scalefactor);
    }
    for ( auto const & grid : gridPointers ) {
      pineappl_grid_scale(grid.second, scalefactor);
    }
  }
  
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
