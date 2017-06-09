import os
import errno    
from ROOT import TFile, TH1D, TCanvas
import ROOT
from utils import read_category_ratios, readMisIDRatios, bin_names_composite, bin_names_composite_nice, read_exclude_bins, bin_names_single
#, mkdir_p, get_bin_name_single, get_bin_name, 
#from plot_pulls import 
from plot_pulls_gen import readMisIDRatiosGen, readCategoryRatiosGen, make_pull_plot_21, makeCatRatiosFrom6
from matrix_solver import calculate_solution, print_ratios_latex

NSIGMAS = 1.2815515 #Corresponds to p-value of 0.1
EXCLUDED_FILE = "../data/excluded_categories.txt"


def get_bin_nr_composite(cat):
  return bin_names_composite.index(cat)

def get_bin_nr_single(cat):
  return bin_names_single.index(cat)

def select_categories(chi2s):
  f = open(EXCLUDED_FILE, "w")
  for (k,v) in chi2s.items():
    #print k, v > NSIGMAS
    if v > NSIGMAS:
      f.write("%s\n" % k)
  f.close()

if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)
  infile = "/hdfs/local/ttH_2tau/andres/ttHAnalysis/2016/histosCF_summer_June6/histograms/charge_flip/histograms_harvested_stage2_charge_flip_Tight.root"
  FITNAME = "summer_June6"
   
  (misIDRatiosNum, misIDRatios) = readMisIDRatiosGen(infile)
  catRatiosNum, catRatios = readCategoryRatiosGen(infile)
  print_ratios_latex(misIDRatiosNum, "gen")   
  """for k in bin_names_composite_nice:
    print k, catRatios[k]
  for k in bin_names_single:
    print k, misIDRatios[k]"""
  chi2s = make_pull_plot_21(misIDRatios, catRatios, mydir = "pull_plots_all", y_range = (-0.001, 0.011))
  
  (misIDRatiosNum, misIDRatios) = readMisIDRatiosGen(infile, rec = "_rec")
  catRatiosNum, catRatios = readCategoryRatiosGen(infile, gen = "gen_rec")
  print_ratios_latex(misIDRatiosNum, "gen_rec")   
  """for k in bin_names_composite_nice:
    print k, catRatios[k]
  for k in bin_names_single:
    print k, misIDRatios[k]"""
  chi2s = make_pull_plot_21(misIDRatios, catRatios, mydir = "pull_plots_all", y_range = (-0.001, 0.011), name = "gen_rec")
  
  #closure test:
  (exclude_bins, exclude_bins_num) = read_exclude_bins(EXCLUDED_FILE)
  
  catRatiosNum, catRatios = makeCatRatiosFrom6(misIDRatios, exclude_bins)
  calculate_solution(catRatiosNum, exclude_bins_num, FITNAME, "closure", "pseudodata")
  #file_misId = "fit_output_pseudodata_%s/fit_res%s.root" % (FITNAME, fittypestring)
  
  
  #select_categories(chi2s)
  for exclude in [False, True]:
    fittypestring = "_gen_rec"
    file_misId = "fit_output_pseudodata_%s/fit_res%s.root" % (FITNAME, fittypestring)
    name = "gen_rec_fit"
    exclude_bins, exclude_bins_num = [], []
    if exclude:
      (exclude_bins, exclude_bins_num) = read_exclude_bins(EXCLUDED_FILE)
      name += "_exclusions"
      fittypestring += "_exclusions"
    #file_misId = "fit_output_pseudodata_%s/fit_res%s.root" % (FITNAME, fittypestring)
    catRatiosNum, catRatios = readCategoryRatiosGen(infile, exclude_bins, gen = "gen_rec")
    #print catRatiosNum
    #print ":::"
    #print exclude_bins
    calculate_solution(catRatiosNum, exclude_bins_num, FITNAME, fittypestring, "pseudodata")
    print file_misId, fittypestring, exclude
    misIDRatios = readMisIDRatios(file_misId)
    make_pull_plot_21(misIDRatios, catRatios, mydir = "pull_plots_all", name = name, y_range = (-0.001, 0.011), excluded = exclude_bins)
    
  FITTYPE = "" #can use also "shapes" or "hybrid" here if the fit results exist
  for datastring in ["pseudodata", "data"]:
    fittypestring = FITTYPE
    for exclude in [False, True]:
      if len(FITTYPE) > 0: fittypestring = "_"+FITTYPE
      file_cats = "fit_output_%s_%s/results_cat%s.txt" % (datastring, FITNAME, fittypestring)
      name = datastring
      exclude_bins, exclude_bins_num = [], []
      if exclude:
        (exclude_bins, exclude_bins_num) = read_exclude_bins(EXCLUDED_FILE)
        name += "_exclusions"
        fittypestring += "_exclusions"
      file_misId = "fit_output_%s_%s/fit_res%s.root" % (datastring, FITNAME, fittypestring)
      catRatiosNum, catRatios = read_category_ratios(file_cats, exclude_bins)
      calculate_solution(catRatiosNum, exclude_bins_num, FITNAME, fittypestring, datastring)
      misIDRatios = readMisIDRatios(file_misId)
      
      make_pull_plot_21(misIDRatios, catRatios, mydir = "pull_plots_all", name = name, y_range = (-0.001, 0.011), excluded = exclude_bins)      

