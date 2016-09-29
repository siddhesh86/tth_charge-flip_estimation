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
	    {0,"BB_LL"},{1,"BB_ML"},{2,"BB_MM"},{3,"BB_HL"},{4,"BB_HM"},{5,"BB_HH"},
      {6,"EE_LL"},{7,"EE_ML"},{8,"EE_MM"},{9,"EE_HL"},{10,"EE_HM"},{11,"EE_HH"},
      {12,"BE_LL"},{13,"BE_ML"},{14,"EB_ML"},{15,"BE_MM"},{16,"BE_HL"},{17,"EB_HL"},
      {18,"BE_HM"},{19,"EB_HM"},{20,"BE_HH"},};


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
            //Don't add any signal shape uncertainties in this case (all are empty & the fit doesn't work)
         }
         else{  //Normal bins do have signal
            for (auto shape_syst : shape_systs) {
                cb.cp().channel({"SS"}).bin({bin}).signals().AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
            }
            for (auto shape_syst_sig : shape_systs_signal) {
                cb.cp().channel({"SS"}).bin({bin}).signals().AddSyst(cb, shape_syst_sig, "shape", SystMap<>::init(1.00));
            }        
        }
        //Backgrounds for all
        for (auto shape_syst : shape_systs) {
            cb.cp().channel({"SS"}).bin({bin}).backgrounds().AddSyst(cb, shape_syst, "shape", SystMap<>::init(1.00));
        }        
    }

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
    return 0;
}


int shapes_main(bool do_parametric = true) {
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
    //bkg_procs["SS"] = {"DY_fake", "WJets", "Singletop", "Diboson", "TTbar"};
    bkg_procs["SS"] = {"bkg"};
    bkg_procs["OS"] = {"bkg"};
    //bkg_procs["OS"] = {"DY_fake", "WJets", "Singletop", "Diboson", "TTbar"};
    map<string, VString> sig_procs;
    sig_procs["SS"] = {"DY"};
    sig_procs["OS"] = {"DY"};

    map<string, Categories> cats;
    cats["SS_13TeV"] = {
	    {0,"BB_LL"}};/*,{1,"BB_ML"},{2,"BB_MM"},{3,"BB_HL"},{4,"BB_HM"},{5,"BB_HH"},
      {6,"EE_LL"},{7,"EE_ML"},{8,"EE_MM"},{9,"EE_HL"},{10,"EE_HM"},{11,"EE_HH"},
      {12,"BE_LL"},{13,"BE_ML"},{14,"EB_ML"},{15,"BE_MM"},{16,"BE_HL"},{17,"EB_HL"},
      {18,"BE_HM"},{19,"EB_HM"},{20,"BE_HH"},};
    cats["OS_13TeV"] = {
	    {0,"BB_LL"},{1,"BB_ML"},{2,"BB_MM"},{3,"BB_HL"},{4,"BB_HM"},{5,"BB_HH"},
      {6,"EE_LL"},{7,"EE_ML"},{8,"EE_MM"},{9,"EE_HL"},{10,"EE_HM"},{11,"EE_HH"},
      {12,"BE_LL"},{13,"BE_ML"},{14,"EB_ML"},{15,"BE_MM"},{16,"BE_HL"},{17,"EB_HL"},
      {18,"BE_HM"},{19,"EB_HM"},{20,"BE_HH"},};
    */

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
  cb.cp().channel({"SS", "OS"}).AddSyst(cb, "lumi", "lnN", SystMap<>::init(1.05));
  /*cb.cp().channel({"SS"}).backgrounds()
	  .AddSyst(cb, "lumi", "lnN", SystMap<>::init(1.05));
    */
  //syst on normalization
  //cb.cp().channel({"SS"}).signals()
//	  .AddSyst(cb, "DY_norm", "lnN", SystMap<>::init(1.5));
	/*cb.cp().channel({"SS"}).process({"DY_fake"})
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
    */
  //syst on normalization
	//cb.cp().channel({"OS"}).signals()
	//  .AddSyst(cb, "DY_norm", "lnN", SystMap<>::init(1.5));
	/*cb.cp().channel({"OS"}).process({"DY_fake"})
	  .AddSyst(cb, "fake_norm", "lnN", SystMap<>::init(1.5));
	cb.cp().channel({"OS"}).process({"WJets"})
	  .AddSyst(cb, "wjets_norm", "lnN", SystMap<>::init(1.5));
	cb.cp().channel({"OS"}).process({"Singletop"})
	  .AddSyst(cb, "singletop_norm", "lnN", SystMap<>::init(1.5));
	cb.cp().channel({"OS"}).process({"Diboson"})
	  .AddSyst(cb, "diboson_norm", "lnN", SystMap<>::init(1.5));
    cb.cp().channel({"OS"}).process({"TTbar"})
	  .AddSyst(cb, "ttbar_norm", "lnN", SystMap<>::init(1.5));*/
	/*cb.cp().channel({"OS"}).process({"QCD"})
	  .AddSyst(cb, "QCD_norm", "lnN", SystMap<>::init(1.10));
	*/


    // This function modifies every entry to have a standardised bin name of
	// the form: {analysis}_{channel}_{bin_id}_{era}
	// which is commonly used in the htt analyses
	ch::SetStandardBinNames(cb);
	
	// First we generate a set of bin names:
	set<string> bins = cb.bin_set();
	//set<string> hm_bins = cb_hm.SetFromObs(mem_fn(&ch::Observation::bin));
	
	// This method will produce a set of unique bin names by considering all
	// Observation, Process and Systematic entries in the CombineHarvester
	// instance.
    /*if (create_asimov) {
        for (auto const& b : lm_bins) {
          ch::CombineHarvester tmp = std::move(
              cb.cp().bin({b}).backgrounds().syst_type({"shape"}, false));
          TH1F tot_bkg = tmp.GetShape();
          // double bkg_rate = tot_bkg.Integral();
          tot_bkg.Scale(tmp.GetObservedRate()/tot_bkg.Integral());
          // tot_bkg.Add(&sig, 0.0);
          // for (int i = 1; i <= tot_bkg.GetNbinsX(); ++i) {
          //   tot_bkg.SetBinContent(i, std::floor(tot_bkg.GetBinContent(i) + 0.5));
          // }
          tot_bkg.Scale(1.0/tot_bkg.Integral());
          tmp.ForEachObs([&](ch::Observation * obs) {
            obs->set_shape(ch::make_unique<TH1F>(tot_bkg), false);
          });
        }
    } */

    /*ch::CombineHarvester cb_hm = cb;
    
    cb_hm.ForEachObs([&](ch::Observation* in) {
        in->set_bin(in->bin() + "_hm");
    });
    cb_hm.ForEachProc([&](ch::Process* in) {
        in->set_bin(in->bin() + "_hm");
    });
    cb_hm.ForEachSyst([&](ch::Systematic* in) {
        in->set_bin(in->bin() + "_hm");
    });*/

  //set<string> hm_bins = cb_hm.SetFromObs(mem_fn(&ch::Observation::bin));

  /*cb.cp().bin_id({8}).VariableRebin(
    {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120,
      130, 140, 150, 160, 170, 180, 190, 200,225,250,275,300});
  cb.cp().bin_id({9}).VariableRebin(
      {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 250, 300});*/

  // Drop shape uncerts on signal in the hm for now,
  // the fb versions aren't in the datacards and
  // we want to keep the comparison fair
  /*cb_hm.FilterSysts([&](ch::Systematic const* p) {
    return p->type() == "shape" && p->signal();
  });*/

  /*std::cout << "Doing bbb for the low mass...\n";
  cb.cp().process({"W", "QCD", "ZTT", "ZL", "ZJ", "TT", "VV"})
      .MergeBinErrors(0.05, 0.4);
  cb.cp().process({"W", "QCD", "ZTT", "ZL", "ZJ", "TT", "VV"})
      .AddBinByBin(0.05, true, &cb);
  std::cout << "...done\n";
    */


  if (do_parametric) {
    for (string chn : chns) {
        /*cb.FilterSysts([&](ch::Systematic const* p) {
          return p->type() == "shape";
        });*/

        /*std::vector<double> bins_hm;
        double x = 60.;
        while (x < 121.) {
          bins_hm.push_back(x);
          x += 1.;
        }*/
        //cb_hm.VariableRebin(bins_hm);
        cb.PrintAll();

        RooWorkspace ws(Form("htt_%s",chn.data()), Form("htt_%s",chn.data()));
        //for (auto const& b : hm_bins) {
        auto bins = cb.cp().channel({chn}).bin_set();
        for (auto const& b : bins) {
          ch::CombineHarvester tmp = std::move(cb.cp().bin({b}).backgrounds());
          TH1F tot_bkg = tmp.GetShape();
          double bkg_error = 0.;
          double bkg_rate = tot_bkg.IntegralAndError(1, tot_bkg.GetNbinsX(), bkg_error);
          //double bkg_uncert = tmp.GetUncertainty();
           // tot_bkg.Integral();
          std::cout << "bkg_rate: " << bkg_rate << "\t" << bkg_error << "\t" << tmp.GetUncertainty() <<  "\n";
          // Here we just play a small trick to reduce the background to one process,
          // DY_fake, which we then change to the new bkg pdf
          tmp.process({"DY_fake"});
          tmp.ForEachProc([&](ch::Process *proc) {
            proc->set_process("bkg");
            proc->set_rate(1.0);
            proc->set_shape(nullptr, false);
          });


          RooRealVar mtt("CMS_th1x", "CMS_th1x", 0,
                         static_cast<float>(tot_bkg.GetNbinsX()));
          RooRealVar lA("bkgpdf_a", "", 50, 0.01, 1000);
          RooRealVar lB("bkgpdf_b", "", 50, -10500,
                        10500);
          RooRealVar lC("bkgpdf_c", "", 50, 0,
                        10500);
          std::string cond = "((bkgpdf_a+bkgpdf_b*0.001*CMS_th1x) > 0) * ";
          std::string fn = cond += "exp(-CMS_th1x/(bkgpdf_a+bkgpdf_b*0.001*CMS_th1x))";
          std::cout << "fn = " << fn << std::endl;
          RooGenericPdf bkg_pdf("bkgpdf", fn.c_str(),
                                RooArgList(mtt, lA, lB));
          // RooRealVar bkg_norm((b + "_bkgpdf_norm").c_str(), "", bkg_rate, 0., bkg_rate*10.);
          RooRealVar bkg_norm("bkgpdf_norm", "", bkg_rate);
          /*tmp.cp().process({"bkg"})
              .AddSyst(cb, "CMS_htt_norm", "lnN", ch::syst::SystMap<>::init
              (1.0 + (bkg_uncert/bkg_rate)));*/
          bkg_norm.setConstant();
          ws.import(bkg_pdf);
          ws.import(bkg_norm);
        }
        //cb.process({"ggH", "bbH", "bkg"});

        cb.AddWorkspace(ws);
        cb.cp().backgrounds().ExtractPdfs(cb, "htt_"+chn, "$CHANNEL_bkgpdf");
        // cb_hm.PrintAll();
     }
  }
    for (string chn : chns) {
                string folder = ("/home/andres/tth/chargeFlip/CMSSW_7_4_7/src/tthAnalysis/ChargeFlipEstimation/bin/output_shapes_pseudodata_eleESER_mva_0_6_notrig/cards/"+chn+"cards/").c_str();
                boost::filesystem::create_directories(folder);
                boost::filesystem::create_directories(folder + "/common");
	    TFile output((folder + "/common/htt_" + chn + ".input.root").c_str(),
			    "RECREATE");
	    auto bins = cb.cp().channel({chn}).bin_set();
	    for (auto b : bins) {
		    cout << ">> Writing datacard for bin: " << b << "\r" << endl;
		    cout << "____________" << endl;
		    cout << folder + "/" + b + ".txt" << endl;
		    cb.PrintAll();
		    //cout << cb.cp().channel({chn}) << endl;
		    //cout << cb.cp().channel({chn}).bin({b}) << endl;
		    cb.cp().channel({chn}).bin({b}).WriteDatacard(
				    folder + "/" + b + ".txt", output);
	    }
        output.Close();
    }
    return 0;
}



int main() {
    histograms_main();
    //return shapes_main();    
}

