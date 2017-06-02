import os
import errno    
from ROOT import TFile, TH1D, TCanvas
import ROOT
from utils import read_category_ratios, bin_names_composite, bin_names_composite_nice, bin_names_single, mkdir_p, get_bin_name_single, get_bin_name, get_bin_name_nice, get_component_cats
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
    #print "gen_ratio/pt_eta_%s" % p
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

def readCategoryRatiosGen(infile, exclude_bins = []):
  f = TFile(infile) 
  i = 0
  os_err = ROOT.Double()
  ss_err = ROOT.Double()
  ratios = {}
  ratios_num = []
  for cat in bin_names_composite:
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
    if err == 0. : err = 1e-8
    #print "%d, %f, %f, %f" % (i, ratio, err, err)
    #print get_bin_name_nice(i), get_bin_name_nice(i) in exclude_bins
    if get_bin_name_nice(i) in exclude_bins: 
      i += 1
      continue
    else:    
      ratios[get_bin_name_nice(i)] = (ratio, err, err)
      ratios_num.append((ratio, err, err))
      i+=1
  return (ratios_num, ratios)


def make_pull_plot_gen(category, misIDRatios, catRatios):
  pull_plot = TH1D(category, category, 6, 0, 6 )
  others_plot = TH1D(category+"others", category+"others", 6, 0, 6 )
  true_plot = TH1D(category+"true", category+"true", 6, 0, 6 )
  bin_names = get_all_containers(category)
  for b in range(1, len(bin_names)+1):
    pull_plot.GetXaxis().SetBinLabel(b,bin_names[b-1])
    #pull_plot.SetAxisRange(-0.006, 0.006,"Y")
    #others_plot.SetAxisRange(-0.006, 0.006,"Y")
    others_plot.GetXaxis().SetBinLabel(b,bin_names[b-1])
    true_plot.GetXaxis().SetBinLabel(b,bin_names[b-1])
    (value, err, err_plus) = catRatios[bin_names[b-1]]
    pull_plot.SetBinContent(b, value)
    pull_plot.SetBinError(b, err)

    other = get_other_component(category, bin_names[b-1])
    (valueO, errO, err0_plus) = misIDRatios[other]
    others_plot.SetBinContent(b, valueO)
    others_plot.SetBinError(b, errO)

    true_plot.SetBinContent(b, misIDRatios[category][0])
    #print bin_names[b-1], value, valueO  
  pull_plot.Add(others_plot, -1)
  c = TCanvas("Plot", "Plot", 1920,1080)
  ROOT.gStyle.SetOptStat(0)
  true_plot.SetLineColor(ROOT.kRed)
  true_plot.SetLineWidth(3)
  true_plot.GetYaxis().SetRangeUser(-0.006, 0.006)
  true_plot.Draw()
  pull_plot.SetLineWidth(3)
  pull_plot.Draw("SAME")
  mydir = "pull_plots_gen/"
  mkdir_p(mydir)
  c.SaveAs("%s/%s_pulls.pdf" % (mydir, category))
  c.SaveAs("%s/%s_pulls.png" % (mydir, category))


def make_pull_plot_21(misIDRatios, catRatios, name = "gen", mydir = "pull_plots_21", y_range = None, excluded = []):
  pull_plots = []
  sum_plots = []
  sum_plots_2 = []
  chi2s = {}
  sum_plot = TH1D("sum_plot", "", 21, 0, 21 )
  gen_plot = TH1D("gen_plot", "", 21, 0, 21 )
  c = TCanvas("Plot", "Plot", 1920,1080)
  ROOT.gStyle.SetOptStat(0)
  test1 = TH1D("test1", "test1", 1, 0, 1 )
  test2 = TH1D("test2", "test2", 1, 0, 1 )
  
  for b in range(1, len(bin_names_composite_nice)+1):
    pull_plots.append(TH1D("cats%d"%b, "", 21, 0, 21 ))
    sum_plots.append(TH1D("sums%d"%b, "sums%d"%b, 21, 0, 21 ))
    sum_plots_2.append(TH1D("sums2_%d"%b, "sums2_%d"%b, 21, 0, 21 ))    
    gen_plot.GetXaxis().SetBinLabel(b,bin_names_composite_nice[b-1])
    sum_plot.GetXaxis().SetBinLabel(b,bin_names_composite_nice[b-1])

  for b in range(1, len(bin_names_composite_nice)+1):
    for i in range(1, len(bin_names_composite_nice)+1):
      pull_plots[b-1].GetXaxis().SetBinLabel(i,bin_names_composite_nice[i-1])
      sum_plots[b-1].GetXaxis().SetBinLabel(i,bin_names_composite_nice[i-1])
      sum_plots_2[b-1].GetXaxis().SetBinLabel(i,bin_names_composite_nice[i-1])
    
    if bin_names_composite_nice[b-1] in excluded: continue

    (cat1, cat2) = get_component_cats(bin_names_composite_nice[b-1])
    (value_gen, err, err_plus) = catRatios[bin_names_composite_nice[b-1]]
    pull_plots[b-1].SetBinContent(b, value_gen)
    pull_plots[b-1].SetBinError(b, err)
        
  
    (value, err, err_plus) = misIDRatios[cat1]
    sum_plots[b-1].SetBinContent(b, value)
    sum_plots[b-1].SetBinError(b, err)
    
    (value, err, err_plus) = misIDRatios[cat2]
    sum_plots_2[b-1].SetBinContent(b, value)
    sum_plots_2[b-1].SetBinError(b, err)

    #print pull_plots[b-1].GetBinContent(b), sum_plots[b-1].GetBinContent(b), sum_plots_2[b-1].GetBinContent(b)
    sum_plots[b-1].Add(sum_plots_2[b-1])
    
    test1.SetBinContent(1, pull_plots[b-1].GetBinContent(b))
    test1.SetBinError(1, pull_plots[b-1].GetBinError(b))
    test2.SetBinContent(1, sum_plots[b-1].GetBinContent(b))
    test2.SetBinError(1, sum_plots[b-1].GetBinError(b))
    #print test1.Chi2Test(test2, "WW")
    #chi2s.append(test1.Chi2Test(test2, "WW"))
    chi2s[bin_names_composite_nice[b-1]] = abs(test1.GetBinContent(1) - test2.GetBinContent(1)) / (test1.GetBinError(1) + test2.GetBinError(1))
    #if value_gen > 0 and False: 
    #  sum_plots[b-1].Scale(1/pull_plots[b-1].Integral())
    #  pull_plots[b-1].Scale(1/pull_plots[b-1].Integral())
    
    gen_plot.Add(pull_plots[b-1])
    sum_plot.Add(sum_plots[b-1])

  if y_range:
    gen_plot.SetAxisRange(y_range[0], y_range[1],"Y")
  #others_plot.SetAxisRange(-0.006, 0.006,"Y")
  #others_plot.GetXaxis().SetBinLabel(b,bin_names[b-1])
  #true_plot.GetXaxis().SetBinLabel(b,bin_names[b-1])
    
  gen_plot.SetLineColor(ROOT.kRed)
  gen_plot.SetLineWidth(2)
  sum_plot.SetLineWidth(2)
  gen_plot.Draw("e1")
  sum_plot.Draw("e1 same")
  mkdir_p(mydir)
  c.SaveAs("%s/pulls_%s.pdf" % (mydir, name))
  c.SaveAs("%s/pulls_%s.png" % (mydir, name))
  return chi2s


if __name__ == "__main__":
  infile = "/hdfs/local/ttH_2tau/andres/ttHAnalysis/2016/histosCF_summer2/histograms/charge_flip/histograms_harvested_stage2_charge_flip_Tight.root"  
        
  misIDRatios = readMisIDRatiosGen(infile)
  catRatios = readCategoryRatiosGen(infile)
  #for bin in bin_names_composite:
  #for bin in bin_names_single:
  #  make_pull_plot_gen(bin, misIDRatios, catRatios)
  make_pull_plot_21(misIDRatios, catRatios)


