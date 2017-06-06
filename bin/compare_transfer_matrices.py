from ROOT import TFile, TH2D, TCanvas
import ROOT
from utils import mkdir_p, bin_names_single


def column_sum(histo, column):
  col_sum = 0.
  for row in range(0, histo.GetNbinsY()+2):
    col_sum += histo.GetBinContent(column, row)
  return col_sum


def column_scaler(histo):
  for col in range(0, histo.GetNbinsX()+2):
    col_sum = column_sum(histo, col)
    if col_sum == 0: continue
    for row in range(0, histo.GetNbinsX()+2):
      print col, row, histo.GetBinContent(col, row)
      histo.SetBinContent(col, row, histo.GetBinContent(col, row) / col_sum)
      histo.SetBinError(col, row, histo.GetBinError(col, row) / col_sum)

def set_axes(histo):
  for i in range(1,7):
    histo.GetXaxis().SetBinLabel(i,bin_names_single[i-1])
    histo.GetYaxis().SetBinLabel(i,bin_names_single[i-1])

def plot_transfer_matrix(histo, title):
  c = TCanvas("Plot", "Plot", 1920,1080)
  ROOT.gStyle.SetOptStat(0)
  #true_plot.SetLineColor(ROOT.kRed)
  #true_plot.SetLineWidth(3)
  #true_plot.GetYaxis().SetRangeUser(-0.006, 0.006)
  set_axes(histo)
  histo.GetXaxis().SetTitle("gen")
  histo.GetYaxis().SetTitle("rec")
  histo.SetNameTitle(title, title)
  histo.Draw("colz")
  histo.Draw("text same")
  mydir = "../plots/transfer_matrices/"
  mkdir_p(mydir)
  c.SaveAs("%s/transfer_matrix_%s.pdf" % (mydir, title))
  c.SaveAs("%s/transfer_matrix_%s.png" % (mydir, title))


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)
  infile = "/hdfs/local/ttH_2tau/andres/ttHAnalysis/2016/histosCF_summer_June5/histograms/charge_flip/histograms_harvested_stage2_charge_flip_Tight.root"
  f = TFile(infile)   
  transfer_matrix_gen_SS = f.Get("gen_ratio/transfer_matrix_noflip")
  transfer_matrix_gen_OS = f.Get("gen_ratio/transfer_matrix_flip")

  column_scaler(transfer_matrix_gen_SS)
  column_scaler(transfer_matrix_gen_OS)

  plot_transfer_matrix(transfer_matrix_gen_SS, "SS")
  plot_transfer_matrix(transfer_matrix_gen_OS, "OS")

