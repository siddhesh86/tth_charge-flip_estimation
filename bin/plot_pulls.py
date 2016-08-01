import os
import errno    
from ROOT import TFile, TH1D, TCanvas
import ROOT


bin_names_composite = ["BB_LL", "BB_ML", "BB_MM", "BB_HL", "BB_HM", "BB_HH",
      "EE_LL", "EE_ML", "EE_MM", "EE_HL", "EE_HM", "EE_HH",
      "BE_LL", "BE_ML", "EB_ML", "BE_MM", "BE_HL", "EB_HL",
      "BE_HM", "EB_HM", "BE_HH"]
bin_names_single = ["BL", "BM", "BH", "EL", "EM", "EH"]


def get_bin_nr_composite(cat):
  return bin_names_composite.index(cat)

def get_bin_nr_single(cat):
  return bin_names_single.index(cat)

def contained_in(composite_name):
  comp1 = composite_name[0]+composite_name[3]
  comp2 = composite_name[1]+composite_name[4]
  return (comp1, comp2)


def get_all_containers(component):
  containers = []  
  for c in bin_names_composite:
    (c1, c2) = contained_in(c)
    if c1 == component or c2 == component: containers.append(c)
  return containers




def readMisIdRatios(file_misId):
  ROOT.gROOT.SetBatch(True)
  f = TFile(file_misId)
  misIdHisto = f.Get("chargeMisId")
  ratios = []
  for etaBin in range(1, misIdHisto.GetNbinsY()+1):
    for ptBin in range(1, misIdHisto.GetNbinsX()+1):
      ratios.append((misIdHisto.GetBinContent(ptBin, etaBin), misIdHisto.GetBinError(ptBin, etaBin)))
      #print "MisID (%d, %d): %f" % (etaBin, ptBin, ratio*100)
  return ratios

def readCategoryRatios(file_cats):
  f = open(file_cats)  
  bins = []
  ratios = []
  for line in f.readlines():
    spl = line.split(",")
    if float(spl[2]) > 0:
      ratios.append((float(spl[1]), float(spl[2])))
    else:
      ratios.append((float(spl[1]), 0.01))
  return ratios
    
def get_other_component(category, composite_category):
  comps = contained_in(composite_category)
  if comps[0] == category: return comps[1]
  else: return comps[0]

def make_pull_plot(category, misIDRatios, catRatios):
  pull_plot = TH1D(category, category, 6, 0, 6 );
  others_plot = TH1D(category+"others", category+"others", 6, 0, 6 );
  bin_names = get_all_containers(category)
  for b in range(1, len(bin_names)+1):
    pull_plot.GetXaxis().SetBinLabel(b,bin_names[b-1])
    others_plot.GetXaxis().SetBinLabel(b,bin_names[b-1])
    (value, err) = catRatios[get_bin_nr_composite(bin_names[b-1])]
    pull_plot.SetBinContent(b, value)
    pull_plot.SetBinError(b, err)

    other = get_other_component(category, bin_names[b-1])
    (valueO, errO) = misIDRatios[get_bin_nr_single(other)]
    others_plot.SetBinContent(b, valueO)
    others_plot.SetBinError(b, errO)
    #print bin_names[b-1], value, valueO  
  pull_plot.Add(others_plot, -1)
  c = TCanvas("Plot", "Plot", 800,600)
  ROOT.gStyle.SetOptStat(0)
  pull_plot.Draw()
  c.SaveAs("pull_plots/"+category+"_pulls_data.pdf")
  c.SaveAs("pull_plots/"+category+"_pulls_data.png")

if __name__ == "__main__":
  file_misId = "fit_output_data_errfix/fit_res_data.root"
  file_cats = "fit_output_data_errfix/results_cat.txt"
  misIDRatios = readMisIdRatios(file_misId)
  catRatios = readCategoryRatios(file_cats)

  for bin in bin_names_single:
    make_pull_plot(bin, misIDRatios, catRatios)


