import os
import errno    
from ROOT import TFile, TH1D, TCanvas
import ROOT
from utils import read_category_ratios, bin_names_composite, bin_names_single, mkdir_p, get_bin_name_single, get_bin_name
from plot_pulls import contained_in, get_other_component, get_all_containers

def get_bin_nr_composite(cat):
  return bin_names_composite.index(cat)

def get_bin_nr_single(cat):
  return bin_names_single.index(cat)

def readMisIDRatiosGen(infile, processes = ["DY"]):
  ratios = {}
  f = TFile(infile)
  for p in processes:
    effs = f.Get("gen_ratio/pt_eta_%s" % p)
    print "gen_ratio/pt_eta_%s" % p
    totalHisto = effs.GetTotalHistogram()
    for bin_eta in range(1, totalHisto.GetNbinsY()+1):
      for bin_pt in range(1, totalHisto.GetNbinsX()+1):
        bin = effs.GetGlobalBin(bin_pt, bin_eta)
        eff = effs.GetEfficiency(bin)
        effErrLo = effs.GetEfficiencyErrorLow(bin)
        effErrHi = effs.GetEfficiencyErrorUp(bin)
        ratios[get_bin_name_single(bin_eta, bin_pt)] = (eff, effErrLo, effErrHi)
        #print "Bin (%d, %d): Eff = %f + %f - %f" % (bin_eta, bin_pt, eff * 100, effErrHi * 100, effErrLo*100)
  return ratios

def readCategoryRatiosGen(infile):
  f = TFile(infile) 
  cats = ["BB_LL", "BB_ML", "BB_MM", "BB_HL", "BB_HM", "BB_HH", "EE_LL", "EE_ML", "EE_MM", "EE_HL", "EE_HM", "EE_HH", "BE_LL", "BE_ML", "EB_ML", "BE_MM", "BE_HL", "EB_HL", "BE_HM", "EB_HM", "BE_HH"]
  i = 0
  os_err = ROOT.Double()
  ss_err = ROOT.Double()
  ratios = {}
  for cat in cats:
    histo_OS = f.Get("gen/OS/%s/mass_ll" % cat)
    histo_SS = f.Get("gen/SS/%s/mass_ll" % cat)
    os_count = histo_OS.IntegralAndError(0, histo_OS.GetNbinsX()+2, os_err)
    ss_count = histo_SS.IntegralAndError(0, histo_SS.GetNbinsX()+2, ss_err)
    if os_count > 0:
      ratio = ss_count / (ss_count + os_count)
      err = (ss_count + ss_err) / (ss_count + ss_err + os_count - os_err) - ratio
    else: 
      ratio = 1.
      err = 1.
    #print "%d, %f, %f, %f" % (i, ratio, err, err)
    ratios[get_bin_name(i)] = (ratio, err, err)
    i+=1
  return ratios


def make_pull_plot_gen(category, misIDRatios, catRatios):
  pull_plot = TH1D(category, category, 6, 0, 6 );
  others_plot = TH1D(category+"others", category+"others", 6, 0, 6 );
  bin_names = get_all_containers(category)
  for b in range(1, len(bin_names)+1):
    pull_plot.GetXaxis().SetBinLabel(b,bin_names[b-1])
    others_plot.GetXaxis().SetBinLabel(b,bin_names[b-1])
    (value, err, err_plus) = catRatios[bin_names[b-1]]
    pull_plot.SetBinContent(b, value)
    pull_plot.SetBinError(b, err)

    other = get_other_component(category, bin_names[b-1])
    (valueO, errO, err0_plus) = misIDRatios[other]
    others_plot.SetBinContent(b, valueO)
    others_plot.SetBinError(b, errO)
    #print bin_names[b-1], value, valueO  
  pull_plot.Add(others_plot, -1)
  c = TCanvas("Plot", "Plot", 800,600)
  ROOT.gStyle.SetOptStat(0)
  pull_plot.Draw()
  mydir = "pull_plots_gen/"
  mkdir_p(mydir)
  c.SaveAs("%s/%s_pulls.pdf" % (mydir, category))
  c.SaveAs("%s/%s_pulls.png" % (mydir, category))

if __name__ == "__main__":
  infile = "/hdfs/local/ttH_2tau/andres/ttHAnalysis/2016/histosCF_summer2/histograms/charge_flip/histograms_harvested_stage2_charge_flip_Tight.root"  
        
  misIDRatios = readMisIDRatiosGen(infile)
  catRatios = readCategoryRatiosGen(infile)

  #for bin in bin_names_composite:
  for bin in bin_names_single:
    make_pull_plot_gen(bin, misIDRatios, catRatios)


