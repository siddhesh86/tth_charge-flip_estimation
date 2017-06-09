from subprocess import call, check_output
import os
import errno    
import math
from ROOT import TFile, TH1D
import ROOT
from make_fits import read_fit_result, mkdir_p

datacard_dir = "output_data_mu_testNew"


def combine_cards():
  dc_dir = "%s/cards" % datacard_dir
  ws_dir = "%s/workspaces" % datacard_dir
  mkdir_p(dc_dir)
  mkdir_p(ws_dir)
  fit_results = []
  for bin in range(1):
    this_card = "%s/card_%d.txt" % (dc_dir, bin)
    this_ws = "%s/workspace_%d.root" %(ws_dir, bin)
    print "combineCards.py pass=%s/SScards/htt_SS_%d_13TeV.txt fail=%s/cards/OScards/htt_OS_%d_13TeV.txt > %s" % (dc_dir, bin, datacard_dir, bin, this_card)
    if bin == 6 and "data" in dc_dir: #SS is totally empty
      fit_results.append((0,1,1))
      continue
      #if bin == 6 and "_data" in dc_dir: #Exclude shape systematics when all are zero
      #call("combineCards.py -s pass=%s/SScards/htt_SS_%d_13TeV.txt fail=%s/cards/OScards/htt_OS_%d_13TeV.txt > %s" % (dc_dir, bin, datacard_dir, bin, this_card), shell = True)
    else:
      call("combineCards.py pass=%s/SScards/htt_SS_%d_13TeV.txt fail=%s/cards/OScards/htt_OS_%d_13TeV.txt > %s" % (dc_dir, bin, datacard_dir, bin, this_card), shell = True)
    
    #Hack to prevent PostFitShapesFromWorkspace messing up directory structure
    call("cp %s ." % (this_card), shell = True)
    print "text2workspace.py %s -o %s -P HiggsAnalysis.CombinedLimit.TagAndProbeModel:tagAndProbe" % (this_card, this_ws)
    call("text2workspace.py %s -o %s -P HiggsAnalysis.CombinedLimit.TagAndProbeModel:tagAndProbe" % (this_card, this_ws), shell = True)
    fit_dir = "./fit_%s/bin%d" % (datacard_dir, bin)
    mkdir_p(fit_dir)
    call("combine -v0 -M MaxLikelihoodFit %s --out %s --plots --saveNormalizations --skipBOnlyFit --saveShapes --saveWithUncertainties --robustFit 1" % (this_ws, fit_dir), shell = True)
    print "combine -v0 -M MaxLikelihoodFit %s --out %s --plots --saveNormalizations --skipBOnlyFit --saveShapes --saveWithUncertainties --robustFit 1" % (this_ws, fit_dir)
    print "PostFitShapesFromWorkspace -d %s -w %s -o %s/output_postfit.root -f %s/mlfit.root:fit_s --postfit --sampling --print" % ("card_%d.txt" % (bin), this_ws, fit_dir, fit_dir)
    call("PostFitShapesFromWorkspace -d %s -w %s -o %s/output_postfit.root -f %s/mlfit.root:fit_s --postfit --sampling" % ("card_%d.txt" % (bin), this_ws, fit_dir, fit_dir), shell=True)
    fit_results.append(read_fit_result("%s/mlfit.root" % fit_dir, "%s/output_postfit.root" % fit_dir))
    #fit_results_mc.append(read_fit_result("%s/mlfit.root" % fit_dir, "%s/output_postfit.root" % fit_dir, isMC = True))
    #print "result", read_fit_result("%s/mlfit.root" % fit_dir, "%s/output_postfit.root" % fit_dir)
    
  
  call("rm card*.txt", shell = True)
  i = 0
  f = open("./fit_%s/results_cat.txt" % (datacard_dir), "w")
  #print fit_results
  for fr in fit_results:
    print "RES: %d %.8f + %.8f - %.8f" % (i, fr[0]/2, fr[1]/2, fr[2]/2)
    f.write("%d, %.8f, %.8f, %.8f\n" % (i, fr[0]/2, fr[1]/2, fr[2]/2))
    i += 1
  f.close()
  
  
combine_cards()
