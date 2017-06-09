#include <string>
#include <map>
#include <set>
#include <iostream>
#include <utility>
#include <vector>
#include <cstdlib>
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Observation.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/BinByBin.h"

#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooNumConvPdf.h>
#include <RooMsgService.h>
#include <RooWorkspace.h>
#include <RooBreitWigner.h>
#include <RooCBShape.h>
#include <RooExponential.h>
#include <RooGenericPdf.h>

using namespace std;

int histograms_main() {
	//! [part1]
	// First define the location of the "auxiliaries" directory where we can
	// source the input files containing the datacard shapes
	string aux_shapes = "/home/andres/ttHAnalysis/2016/histosCF_mu_summer/datacards/charge_flip/";

	// Create an empty CombineHarvester instance that will hold all of the
	// datacard configuration and histograms etc.
	ch::CombineHarvester cb;
	typedef vector<pair<int, string>> Categories;
	typedef vector<string> VString;

	// Here we will just define two categories for an 8TeV analysis. Each entry in
	// the vector below specifies a bin name and corresponding bin_id.
	VString chns = {"SS", "OS"};

    map<string, VString> bkg_procs;
    bkg_procs["SS"] = {"DY_fake", "WJets", "Singletop", "Diboson", "TTbar"};
    bkg_procs["OS"] = {"DY_fake", "WJets", "Singletop", "Diboson", "TTbar"};
    map<string, VString> sig_procs;
    sig_procs["SS"] = {"DY"};
    sig_procs["OS"] = {"DY"};

    map<string, Categories> cats;
    cats["SS_13TeV"] = {
	    {0,"total"}};
    cats["OS_13TeV"] = {
	    {0,"total"}};

    VString shape_systs = {};
    //VString shape_systs = {"CMS_ttHl_muonESBarrel1", "CMS_ttHl_muonESBarrel2",
    //    "CMS_ttHl_muonESEndcap2", "CMS_ttHl_muonESEndcap2"};

    VString shape_systs_signal = {/*"CMS_ttHl_muonER"*/
    };

    cout << ">> Creating processes and observations...\n";
    for (string era : {"13TeV"}) {
	    for (auto chn : chns) {
		    cb.AddObservations(
				    {"*"}, {"htt"}, {era}, {chn}, cats[chn+"_"+era]);
		    cb.AddProcesses(
				    {"*"}, {"htt"}, {era}, {chn}, bkg_procs[chn], cats[chn+"_"+era], false);
		    cb.AddProcesses(
				    {"*"}, {"htt"}, {era}, {chn}, sig_procs[chn], cats[chn+"_"+era], true);
	    }
    }

    //Some of the code for this is in a nested namespace, so
    // we'll make some using declarations first to simplify things a bit.
    using ch::syst::SystMap;
    using ch::syst::era;
    using ch::syst::bin_id;
    using ch::syst::process;


  //syst on luminosity
  cb.cp().channel({"SS"}).signals()
	  .AddSyst(cb, "lumi", "lnN", SystMap<>::init(1.05));
  cb.cp().channel({"SS"}).backgrounds()
	  .AddSyst(cb, "lumi", "lnN", SystMap<>::init(1.05));

  //syst on normalization
  cb.cp().channel({"SS"}).signals()
	  .AddSyst(cb, "DY_norm", "lnN", SystMap<>::init(5.0));
	cb.cp().channel({"SS"}).process({"DY_fake"})
	  .AddSyst(cb, "fake_norm", "lnN", SystMap<>::init(5.));
	cb.cp().channel({"SS"}).process({"WJets"})
	  .AddSyst(cb, "wjets_norm", "lnN", SystMap<>::init(5.));
	cb.cp().channel({"SS"}).process({"Singletop"})
	  .AddSyst(cb, "singletop_norm", "lnN", SystMap<>::init(5.));
	cb.cp().channel({"SS"}).process({"Diboson"})
	  .AddSyst(cb, "diboson_norm", "lnN", SystMap<>::init(5.));
  cb.cp().channel({"SS"}).process({"TTbar"})
	  .AddSyst(cb, "ttbar_norm", "lnN", SystMap<>::init(5.));
	
  cb.cp().channel({"OS"}).signals()
	  .AddSyst(cb, "lumi", "lnN", SystMap<>::init(1.05));
  cb.cp().channel({"OS"}).backgrounds()
	  .AddSyst(cb, "lumi", "lnN", SystMap<>::init(1.05));

  //syst on normalization
	cb.cp().channel({"OS"}).signals()
	  .AddSyst(cb, "DY_norm", "lnN", SystMap<>::init(5.0));
	cb.cp().channel({"OS"}).process({"DY_fake"})
	  .AddSyst(cb, "fake_norm", "lnN", SystMap<>::init(5.));
	cb.cp().channel({"OS"}).process({"WJets"})
	  .AddSyst(cb, "wjets_norm", "lnN", SystMap<>::init(5.));
	cb.cp().channel({"OS"}).process({"Singletop"})
	  .AddSyst(cb, "singletop_norm", "lnN", SystMap<>::init(5.));
	cb.cp().channel({"OS"}).process({"Diboson"})
	  .AddSyst(cb, "diboson_norm", "lnN", SystMap<>::init(5.));
    cb.cp().channel({"OS"}).process({"TTbar"})
	  .AddSyst(cb, "ttbar_norm", "lnN", SystMap<>::init(5.));
	/*cb.cp().channel({"OS"}).process({"QCD"})
	  .AddSyst(cb, "QCD_norm", "lnN", SystMap<>::init(1.10));
	*/


    auto bins_ss = cb.cp().channel({"SS"}).bin_set();
    /*for (auto bin : bins_ss) {
      std::cout << "bin " << bin << std::endl;
      for (auto shape_syst : shape_systs) {
          cb.cp().channel({"SS"}).bin({bin}).signals().AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
      }
      for (auto shape_syst_sig : shape_systs_signal) {
          cb.cp().channel({"SS"}).bin({bin}).signals().AddSyst(cb, shape_syst_sig, "shape", SystMap<>::init(1.00));
      }        
    }*/
    //Backgrounds for all
    /*for (auto shape_syst : shape_systs) {
        cb.cp().channel({"SS"}).bin({bin}).backgrounds().AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
    }*/        
    //}
    
    auto bins_os = cb.cp().channel({"OS"}).bin_set();
    for (auto bin : bins_os) {
        for (auto shape_syst : shape_systs) {
            cb.cp().channel({"OS"}).bin({bin}).signals().AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
            cb.cp().channel({"OS"}).bin({bin}).backgrounds().AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
        }
        for (auto shape_syst_sig : shape_systs_signal) {
            cb.cp().channel({"OS"}).bin({bin}).signals().AddSyst(cb, shape_syst_sig, "shape", SystMap<>::init(1.00));
        }
    }

	cout << ">> Extracting histograms from input root files...\n";
	for (string era : {"13TeV"}) {
		for (string chn : chns) {
			string file = aux_shapes+ "/prepareDatacards_pseudodata_charge_flip_mass_ll.root";
			cb.cp().channel({chn}).backgrounds().ExtractShapes(
					file, "ttH_charge_flip_" + chn+"_$BIN/x_$PROCESS", "ttH_charge_flip_" + chn+"_$BIN/x_$PROCESS_$SYSTEMATIC");
			cb.cp().channel({chn}).signals().ExtractShapes(
					file, "ttH_charge_flip_" + chn +"_$BIN/x_$PROCESS", "ttH_charge_flip_" + chn+"_$BIN/x_$PROCESS_$SYSTEMATIC");
		}
	}


	auto bbb = ch::BinByBinFactory()
  	.SetAddThreshold(0.05)
    .SetMergeThreshold(0.5)
	  .SetFixNorm(true);

	bbb.AddBinByBin(cb.cp().backgrounds(), cb);
	bbb.AddBinByBin(cb.cp().signals(), cb);

	// This function modifies every entry to have a standardised bin name of
	// the form: {analysis}_{channel}_{bin_id}_{era}
	// which is commonly used in the htt analyses
	ch::SetStandardBinNames(cb);
	
	// First we generate a set of bin names:
	set<string> bins = cb.bin_set();
	// This method will produce a set of unique bin names by considering all
	// Observation, Process and Systematic entries in the CombineHarvester
	// instance.



    for (string chn : chns) {
                string folder = ("/home/andres/tth/chargeFlip/CMSSW_7_4_7/src/tthAnalysis/ChargeFlipEstimation/bin/output_pseudodata_mu_summer_May10/cards/"+chn+"cards/").c_str();
                boost::filesystem::create_directories(folder);
                boost::filesystem::create_directories(folder + "/common");
	    TFile output((folder + "/common/htt_" + chn + ".input.root").c_str(),
			    "RECREATE");
	    auto bins = cb.cp().channel({chn}).bin_set();
	    for (auto b : bins) {
		    cout << ">> Writing datacard for bin: " << b << "\r" << flush;
		    cb.cp().channel({chn}).bin({b}).WriteDatacard(
				    folder + "/" + b + ".txt", output);
	    }
     output.Close();
    }
    return 0;
}


int main() {
    histograms_main();
}

