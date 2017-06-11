import os
import errno    
from ROOT import TFile, TH1D, TCanvas
import ROOT
from utils import read_category_ratios, readMisIDRatios, bin_names_composite, bin_names_composite_nice, read_exclude_bins, bin_names_single
#, mkdir_p, get_bin_name_single, get_bin_name, 
#from plot_pulls import 
from plot_pulls_gen import readMisIDRatiosGen, readCategoryRatiosGen, make_pull_plot_21, makeCatRatiosFrom6
from matrix_solver import calculate_solution, print_ratios_latex

"""@file docstring
Script for plotting pulls comparing result from 21 and 6 categories
Selects which categories of 21 to drop because of correlations and solves the equations for different cases

@author Andres Tiko <andres.tiko@cern.ch>
"""

#Number of sigmas difference to consider fit results not compatible
NSIGMAS = 1.2815515 #Corresponds to p-value of 0.1
#File to store list of excluded categories
EXCLUDED_FILE = "../data/excluded_categories.txt"

#Writes to the specified file list of categories to exclude
def select_categories(chi2s):
  f = open(EXCLUDED_FILE, "w")
  for (k,v) in chi2s.items():
    #print k, v > NSIGMAS
    if v > NSIGMAS:
      f.write("%s\n" % k)
  f.close()

if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)
  #Input file with ntuples
  infile = "/hdfs/local/ttH_2tau/andres/ttHAnalysis/2016/histosCF_summer_June6/histograms/charge_flip/histograms_harvested_stage2_charge_flip_Tight.root"
  #Name of fits to be used  
  FITNAME = "summer_June6"
   
  #Read generator-level ratios (misIDRatios: 6 for single electrons, catRatios for 21 categories of double electrons)
  #Numeric values are used for matrix solver
  (misIDRatiosNum, misIDRatios) = readMisIDRatiosGen(infile)
  catRatiosNum, catRatios = readCategoryRatiosGen(infile)
  #Print latex results  
  print_ratios_latex(misIDRatiosNum, "gen")   
  #Makes a pull plot comparing the 21 numbers to sums of respective ones from 6
  print "The ratios for gen-level electrons with gen-level pT and eta"
  chi2s = make_pull_plot_21(misIDRatios, catRatios, mydir = "pull_plots_all", y_range = (-0.001, 0.011))
  
  print "\nThe ratios for gen-level electrons with reconstructed pT and eta"
  (misIDRatiosNum, misIDRatios) = readMisIDRatiosGen(infile, rec = "_rec")
  catRatiosNum, catRatios = readCategoryRatiosGen(infile, gen = "gen_rec")
  print_ratios_latex(misIDRatiosNum, "gen_rec")   
  chi2s = make_pull_plot_21(misIDRatios, catRatios, mydir = "pull_plots_all", y_range = (-0.001, 0.011), name = "gen_rec")
  
  print "\nClosure test to see if we get 6 number back if we construct the 21 from the 6 and then fit"
  print "Turns out this underestimates uncertainty (due to correlations"
  (exclude_bins, exclude_bins_num) = read_exclude_bins(EXCLUDED_FILE)
  catRatiosNum, catRatios = makeCatRatiosFrom6(misIDRatios, exclude_bins)
  calculate_solution(catRatiosNum, exclude_bins_num, FITNAME, "closure", "pseudodata")
  
  #Selects categories to drop and writes them to file
  #Comment out the following line and edit the file manually to specify the categories to be dropped yourself
  select_categories(chi2s)

  print "\nMake pull plots for both cases of not dropping and dropping categories"
  for exclude in [False, True]:
    fittypestring = "_gen_rec"
    file_misId = "fit_output_pseudodata_%s/fit_res%s.root" % (FITNAME, fittypestring)
    name = "gen_rec_fit"
    exclude_bins, exclude_bins_num = [], []
    if exclude:
      (exclude_bins, exclude_bins_num) = read_exclude_bins(EXCLUDED_FILE)
      name += "_exclusions"
      fittypestring += "_exclusions"
    catRatiosNum, catRatios = readCategoryRatiosGen(infile, exclude_bins, gen = "gen_rec")
    calculate_solution(catRatiosNum, exclude_bins_num, FITNAME, fittypestring, "pseudodata")
    misIDRatios = readMisIDRatios(file_misId)
    make_pull_plot_21(misIDRatios, catRatios, mydir = "pull_plots_all", name = name, y_range = (-0.001, 0.011), excluded = exclude_bins)
  

  print "\nFit results for pseudodata and data first without and then with excluding some categories"
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

