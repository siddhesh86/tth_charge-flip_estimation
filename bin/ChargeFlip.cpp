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

using namespace std;

int main() {
	//! [part1]
	// First define the location of the "auxiliaries" directory where we can
	// source the input files containing the datacard shapes
	string aux_shapes = "/home/andres/tth/histograms/histosCF_data_eleESER2/datacards";

	// Create an empty CombineHarvester instance that will hold all of the
	// datacard configuration and histograms etc.
	ch::CombineHarvester cb;
	typedef vector<pair<int, string>> Categories;
	typedef vector<string> VString;

	// Here we will just define two categories for an 8TeV analysis. Each entry in
	// the vector below specifies a bin name and corresponding bin_id.
	VString chns = {"SS", "OS"};

map<string, VString> bkg_procs;
//bkg_procs["SS"] = {"TTZ", "WZ", "TTW", "signal", "Rares", "additional_signal_overlap", "background_data_singletop","background_data_WJets","background_data_WW","background_data_ZZ","background_data_TTJets"};
//bkg_procs["OS"] = {"TTZ", "WZ", "TTW", "signal", "Rares", "additional_signal_overlap", "background_data_singletop","background_data_WJets","background_data_WW","background_data_ZZ","background_data_TTJets"};
bkg_procs["SS"] = {"DY_fake", "WJets", "Singletop", "Diboson", "TTbar"};
bkg_procs["OS"] = {"DY_fake", "WJets", "Singletop", "Diboson", "TTbar"};
map<string, VString> sig_procs;
sig_procs["SS"] = {"DY"};
sig_procs["OS"] = {"DY"};

map<string, Categories> cats;
cats["SS_13TeV"] = {
	{0,"BB_LL"},{1,"BB_ML"},{2,"BB_MM"},{3,"BB_HL"},{4,"BB_HM"},{5,"BB_HH"},
  {6,"EE_LL"},{7,"EE_ML"},{8,"EE_MM"},{9,"EE_HL"},{10,"EE_HM"},{11,"EE_HH"},
  {12,"BE_LL"},{13,"BE_ML"},{14,"EB_ML"},{15,"BE_MM"},{16,"BE_HL"},{17,"EB_HL"},
  {18,"BE_HM"},{19,"EB_HM"},{20,"BE_HH"},};
cats["OS_13TeV"] = {
	//{0,"mt_inclusive"},
	{0,"BB_LL"},{1,"BB_ML"},{2,"BB_MM"},{3,"BB_HL"},{4,"BB_HM"},{5,"BB_HH"},
  {6,"EE_LL"},{7,"EE_ML"},{8,"EE_MM"},{9,"EE_HL"},{10,"EE_HM"},{11,"EE_HH"},
  {12,"BE_LL"},{13,"BE_ML"},{14,"EB_ML"},{15,"BE_MM"},{16,"BE_HL"},{17,"EB_HL"},
  {18,"BE_HM"},{19,"EB_HM"},{20,"BE_HH"},};
/*cats["mm_13TeV"] = {
	//{0,"mm_inclusive"},
	{0,"mm_0jet"}, {1,"mm_1jet_zpt_loose"},
	{2,"mm_1jet_zpt_medium"},{3,"mm_1jet_zpt_tight"},
	{4,"mm_2jet_cp"}, {5,"mm_vbf"},
	{6,"mm_1bjet"}, {7,"mm_2bjet"}, {8,"mm_MSSM_btag"}};
*/

/*VString shape_systs = {"CMS_ttHl_btag_HFUp", 
       "CMS_ttHl_btag_HFDown",	
       "CMS_ttHl_btag_HFStats1Up", 
       "CMS_ttHl_btag_HFStats1Down",
       "CMS_ttHl_btag_HFStats2Up", 
       "CMS_ttHl_btag_HFStats2Down",
       "CMS_ttHl_btag_LFUp", 
       "CMS_ttHl_btag_LFDown",	
       "CMS_ttHl_btag_LFStats1Up", 
       "CMS_ttHl_btag_LFStats1Down",
       "CMS_ttHl_btag_LFStats2Up", 
       "CMS_ttHl_btag_LFStats2Down",
       "CMS_ttHl_btag_cErr1Up",
       "CMS_ttHl_btag_cErr1Down",
       "CMS_ttHl_btag_cErr2Up",
       "CMS_ttHl_btag_cErr2Down",
       "CMS_ttHl_JESUp",
       "CMS_ttHl_JESDown",
       "CMS_ttHl_tauESUp",
       "CMS_ttHl_tauESDown"
};*/

VString shape_systs = {"CMS_ttHl_electronESBarrel",
    "CMS_ttHl_electronESEndcap"};

VString shape_systs_signal = {"CMS_ttHl_electronER"
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
	  .AddSyst(cb, "DY_norm", "lnN", SystMap<>::init(1.5));
	cb.cp().channel({"SS"}).process({"DY_fake"})
	  .AddSyst(cb, "fake_norm", "lnN", SystMap<>::init(1.5));
	cb.cp().channel({"SS"}).process({"WJets"})
	  .AddSyst(cb, "wjets_norm", "lnN", SystMap<>::init(1.5));
	cb.cp().channel({"SS"}).process({"Singletop"})
	  .AddSyst(cb, "singletop_norm", "lnN", SystMap<>::init(1.5));
	cb.cp().channel({"SS"}).process({"Diboson"})
	  .AddSyst(cb, "diboson_norm", "lnN", SystMap<>::init(1.5));
  cb.cp().channel({"SS"}).process({"TTbar"})
	  .AddSyst(cb, "ttbar_norm", "lnN", SystMap<>::init(1.5));
	
  cb.cp().channel({"OS"}).signals()
	  .AddSyst(cb, "lumi", "lnN", SystMap<>::init(1.05));
  cb.cp().channel({"OS"}).backgrounds()
	  .AddSyst(cb, "lumi", "lnN", SystMap<>::init(1.05));

  //syst on normalization
	cb.cp().channel({"OS"}).signals()
	  .AddSyst(cb, "DY_norm", "lnN", SystMap<>::init(1.5));
	cb.cp().channel({"OS"}).process({"DY_fake"})
	  .AddSyst(cb, "fake_norm", "lnN", SystMap<>::init(1.5));
	cb.cp().channel({"OS"}).process({"WJets"})
	  .AddSyst(cb, "wjets_norm", "lnN", SystMap<>::init(1.5));
	cb.cp().channel({"OS"}).process({"Singletop"})
	  .AddSyst(cb, "singletop_norm", "lnN", SystMap<>::init(1.5));
	cb.cp().channel({"OS"}).process({"Diboson"})
	  .AddSyst(cb, "diboson_norm", "lnN", SystMap<>::init(1.5));
  cb.cp().channel({"OS"}).process({"TTbar"})
	  .AddSyst(cb, "ttbar_norm", "lnN", SystMap<>::init(1.5));
	/*cb.cp().channel({"OS"}).process({"QCD"})
	  .AddSyst(cb, "QCD_norm", "lnN", SystMap<>::init(1.10));
	*/


  auto bins_ss = cb.cp().channel({"SS"}).bin_set();
  for (auto bin : bins_ss) {
     std::cout << "bin " << bin << std::endl;
     if (bin == "EE_LL"){
        /*for (auto shape_syst : shape_systs) {
            cb.cp().channel({"SS"}).bin({bin}).backgrounds().AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
        }*/
     }
     /*else if (bin == "BB_HH" || bin == "BB_ML" || bin == "BE_HL" || bin == "BE_HM" || bin == "BE_LL" || bin == "EB_HL" || bin == "EB_HM" || bin == "EB_ML" || bin == "EE_MM"){   //No W+jets
        for (auto shape_syst : shape_systs) {
            cb.cp().channel({"SS"}).bin({bin}).signals().AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
            cb.cp().channel({"SS"}).bin({bin}).process({"DY_fake", "Singletop", "Diboson", "TTbar"}).AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
        }
        for (auto shape_syst_sig : shape_systs_signal) {
            cb.cp().channel({"SS"}).bin({bin}).signals().AddSyst(cb, shape_syst_sig, "shape", SystMap<>::init(1.00));
        }
     }
     else if (bin == "BE_HH"){   //No W+jets or singletop
        for (auto shape_syst : shape_systs) {
            cb.cp().channel({"SS"}).bin({bin}).signals().AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
            cb.cp().channel({"SS"}).bin({bin}).process({"DY_fake", "Diboson", "TTbar"}).AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
        }
        for (auto shape_syst_sig : shape_systs_signal) {
            cb.cp().channel({"SS"}).bin({bin}).signals().AddSyst(cb, shape_syst_sig, "shape", SystMap<>::init(1.00));
        }
     }
     else if (bin == "EE_ML"){   //No W+jets or ttbar
        for (auto shape_syst : shape_systs) {
            cb.cp().channel({"SS"}).bin({bin}).signals().AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
            cb.cp().channel({"SS"}).bin({bin}).process({"DY_fake", "Singletop", "Diboson"}).AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
        }
        for (auto shape_syst_sig : shape_systs_signal) {
            cb.cp().channel({"SS"}).bin({bin}).signals().AddSyst(cb, shape_syst_sig, "shape", SystMap<>::init(1.00));
        }
     }
     else if (bin == "EE_HH"){   //No W+jets or ttbar?
        for (auto shape_syst : shape_systs) {
            cb.cp().channel({"SS"}).bin({bin}).signals().AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
            cb.cp().channel({"SS"}).bin({bin}).process({"DY_fake", "Singletop", "Diboson"}).AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
        }
        //TODO: cb.cp().channel({"SS"}).bin({bin}).process({"TTbar"}).AddSyst(cb, "CMS_ttHl_electronESEndcap, "shape", SystMap<>::init(1.00));
        for (auto shape_syst_sig : shape_systs_signal) {
            cb.cp().channel({"SS"}).bin({bin}).signals().AddSyst(cb, shape_syst_sig, "shape", SystMap<>::init(1.00));
        }
     }
     else if (bin == "BB_LL"){   //No various stuff
        for (auto shape_syst : shape_systs) {
            cb.cp().channel({"SS"}).bin({bin}).signals().AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
            cb.cp().channel({"SS"}).bin({bin}).process({"DY_fake", "Diboson"}).AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
        }
        //TODO: cb.cp().channel({"SS"}).bin({bin}).process({"Singletop"}).AddSyst(cb, "CMS_ttHl_electronESBarrel", "shape", SystMap<>::init(1.00));
        for (auto shape_syst_sig : shape_systs_signal) {
            cb.cp().channel({"SS"}).bin({bin}).signals().AddSyst(cb, shape_syst_sig, "shape", SystMap<>::init(1.00));
        }
     }*/
     
     else{  //Normal bins
        for (auto shape_syst : shape_systs) {
            cb.cp().channel({"SS"}).bin({bin}).signals().AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
            cb.cp().channel({"SS"}).bin({bin}).backgrounds().AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
        }
        for (auto shape_syst_sig : shape_systs_signal) {
            cb.cp().channel({"SS"}).bin({bin}).signals().AddSyst(cb, shape_syst_sig, "shape", SystMap<>::init(1.00));
        }
        
     }
  }

  
  auto bins_os = cb.cp().channel({"OS"}).bin_set();
  for (auto bin : bins_os) {
     std::cout << "os bin " << bin << std::endl;
     /*if (bin == "EE_LL"){
        /*cb.cp().process({"ZEE"}).bin({bin}).AddSyst(
           cb, "CMS_htt_$CHANNEL_z$CHANNELShape_" + clip + "_mass" + i + "_$ERA",
           "shape", SystMap<>::init(1.0));*/
           ; 
     //}*/
     /*if (bin == "EE_LL" || bin == "EE_ML" || bin == "EE_MM" || bin == "BB_LL" || bin == "BB_HH" || 
            bin == "BE_HH" || bin == "BE_HL" || bin == "BE_HM" || bin == "BE_ML" || 
            bin == "EB_HM" || bin == "EE_HH" || bin == "EE_HL" || bin == "EE_HM"){   //No W+jets*/
     if(false){
        for (auto shape_syst : shape_systs) {
            cb.cp().channel({"OS"}).bin({bin}).signals().AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
            cb.cp().channel({"OS"}).bin({bin}).process({"DY_fake", "Singletop", "Diboson", "TTbar"}).AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
        }
        for (auto shape_syst_sig : shape_systs_signal) {
            cb.cp().channel({"OS"}).bin({bin}).signals().AddSyst(cb, shape_syst_sig, "shape", SystMap<>::init(1.00));
        }
     }
     
     else{  //Normal bins
        for (auto shape_syst : shape_systs) {
            cb.cp().channel({"OS"}).bin({bin}).signals().AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
            cb.cp().channel({"OS"}).bin({bin}).backgrounds().AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
        }
        for (auto shape_syst_sig : shape_systs_signal) {
            cb.cp().channel({"OS"}).bin({bin}).signals().AddSyst(cb, shape_syst_sig, "shape", SystMap<>::init(1.00));
        }
        
     }
  }

	cout << ">> Extracting histograms from input root files...\n";
	for (string era : {"13TeV"}) {
		for (string chn : chns) {
			string file = aux_shapes+ "/prepareDatacards_data_charge_flip_mass_ll.root";
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
	//! [part8]

	//! [part9]
	// First we generate a set of bin names:
	set<string> bins = cb.bin_set();
	// This method will produce a set of unique bin names by considering all
	// Observation, Process and Systematic entries in the CombineHarvester
	// instance.


	for (string chn : chns) {
                string folder = ("/home/andres/tth/chargeFlip/CMSSW_7_4_7/src/tthAnalysis/ChargeFlipEstimation/bin/output_data_eleESER2/cards/"+chn+"cards/").c_str();
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


}
