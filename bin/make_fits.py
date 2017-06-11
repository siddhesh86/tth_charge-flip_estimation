from subprocess import call, check_output
import os
import errno    
import math
from ROOT import TFile, TH1D
import ROOT
from utils import mkdir_p

"""@file docstring
Script for running electron charge flip estimation fit

@author Andres Tiko <andres.tiko@cern.ch>
"""


def read_fit_result(fit_file, postfit_file):
  ROOT.gROOT.SetBatch(True)
  f2 = TFile(postfit_file)
  #POSTFIT SCALED TO SIGNAL r = 1! -> make changes
  fail_h = f2.Get("fail_prefit/DY")
  pass_h = f2.Get("pass_prefit/DY")

  f = TFile(fit_file)
  tree = f.Get("tree_fit_sb")
  
  tree.Draw("mu>>hist")
  hist = ROOT.gDirectory.Get("hist")
  tree.Draw("muErr>>histErr")
  histErr = ROOT.gDirectory.Get("histErr")
  #tree.Draw("muLoErr>>histLo")
  tree.Draw("muErr>>histLo")
  histLo = ROOT.gDirectory.Get("histLo")
  #tree.Draw("muHiErr>>histHi")
  tree.Draw("muErr>>histHi")
  histHi = ROOT.gDirectory.Get("histHi")
  mu = hist.GetMean()
  muErr = histErr.GetMean()
  muLoErr = histLo.GetMean()
  muHiErr = histHi.GetMean()
  
  try:
    #postFit distributions are scaled to scale factor 1, need to multiply by fitted number
    bestFit = mu * pass_h.Integral() / (fail_h.Integral() + pass_h.Integral())
  except AttributeError:
    if fail_h.Integral() > 0:
      bestFit = 0.
    else: raise
  print fail_h.Integral(), fail_h.GetEntries()  
  #print fail_h.Integral(), pass_h.Integral(), bestFit, pass_h.Integral() / (fail_h.Integral() + pass_h.Integral()), f2.Get("fail_postfit"), f2.Get("pass_postfit"), 
    
  """print bestFit - mu * f2.Get("pass_prefit/DY").Integral() / (f2.Get("fail_prefit/DY").Integral() + f2.Get("pass_prefit/DY").Integral()), bestFit, mu * f2.Get("pass_prefit/DY").Integral() / (f2.Get("fail_prefit/DY").Integral() + f2.Get("pass_prefit/DY").Integral()), mu, f2.Get("pass_prefit/DY").Integral() / (f2.Get("fail_prefit/DY").Integral() + f2.Get("pass_prefit/DY").Integral())
  assert(abs(bestFit - mu * f2.Get("pass_prefit/DY").Integral() / (f2.Get("fail_prefit/DY").Integral() + f2.Get("pass_prefit/DY").Integral())) < 0.0001)"""

  print mu, muErr, muHiErr, muLoErr, bestFit
  fitHiErr = max(muHiErr, muErr)/mu * bestFit
  fitLoErr = max(muLoErr, muErr)/mu * bestFit
  if fitHiErr > 1:
    fitHiErr = max(1., mu)
  if fitLoErr > 1:
    fitLoErr = max(1., mu)
  #Use Poisson mean of getting 0 events for uncertainty if no observed events
  if fitHiErr == 0. :
    lambda_poisson0 = -math.log(0.32)
    fitHiErr = lambda_poisson0 / (lambda_poisson0 + fail_h.Integral())
    fitLoErr = fitHiErr
  return (bestFit, fitHiErr, fitLoErr)

#TODO: parallelize this for faster fitting
def make_fits(datacard_dir):
  dc_dir = "%s/cards" % datacard_dir
  ws_dir = "%s/workspaces" % datacard_dir
  mkdir_p(dc_dir)
  mkdir_p(ws_dir)
  fit_results = []
  for bin in range(21):
    this_card = "%s/card_%d.txt" % (dc_dir, bin)
    this_ws = "%s/workspace_%d.root" %(ws_dir, bin)
    #print "combineCards.py pass=%s/SScards/htt_SS_%d_13TeV.txt fail=%s/cards/OScards/htt_OS_%d_13TeV.txt > %s" % (dc_dir, bin, datacard_dir, bin, this_card)
    
    #1. step: combine SS and OS datacatds    
    if bin == 6 and "data" in dc_dir: #SS is totally empty
      fit_results.append((0,1,1))
      continue
    else:
      call("combineCards.py pass=%s/SScards/htt_SS_%d_13TeV.txt fail=%s/cards/OScards/htt_OS_%d_13TeV.txt > %s" % (dc_dir, bin, datacard_dir, bin, this_card), shell = True)
    
    #Hack to prevent PostFitShapesFromWorkspace messing up directory structure - copy datacard to current directory
    call("cp %s ." % (this_card), shell = True)
    #print "text2workspace.py %s -o %s -P HiggsAnalysis.CombinedLimit.TagAndProbeModel:tagAndProbe" % (this_card, this_ws)
    #2. Make Roofit workspace from datacard
    call("text2workspace.py %s -o %s -P HiggsAnalysis.CombinedLimit.TagAndProbeModel:tagAndProbe" % (this_card, this_ws), shell = True)
    #Specify outpud directory for fit results - TODO: move out of bin-directory 
    fit_dir = "./fit_%s/bin%d" % (datacard_dir, bin)
    mkdir_p(fit_dir)

    #3. Perform fit with combine
    #Default fit settings specified in else clause
    #But always do not give convergence - settings for specific fits defined here, adjust as necessary:
    if ((bin in [0,17]) and "_pseudodata" in fit_dir and "summer" in fit_dir):
        specific_settings = "--minimizerStrategy 0"
    elif ((bin in [0,13]) and "_pseudodata" in fit_dir and "summer" in fit_dir):
        #Fit does not work at all
        fit_results.append((0,1,1))
        continue
    elif ((bin in [14]) and "_data" in fit_dir and "summer" in fit_dir):
        specific_settings = "--minimizerStrategy 0"
    else:
        specific_settings = ""
        
    command = "combine -v0 -M MaxLikelihoodFit %s --out %s --plots --saveNormalizations --skipBOnlyFit --saveShapes --saveWithUncertainties %s" % (this_ws, fit_dir, specific_settings)
    #Call combine
    call(command, shell = True)

    #print "PostFitShapesFromWorkspace -d %s -w %s -o %s/output_postfit.root -f %s/mlfit.root:fit_s --postfit --sampling --print" % ("card_%d.txt" % (bin), this_ws, fit_dir, fit_dir)
    #4. Create postfit plots
    call("PostFitShapesFromWorkspace -d %s -w %s -o %s/output_postfit.root -f %s/mlfit.root:fit_s --postfit --sampling" % ("card_%d.txt" % (bin), this_ws, fit_dir, fit_dir), shell=True)
    #5. Add to list of fit results
    fit_results.append(read_fit_result("%s/mlfit.root" % fit_dir, "%s/output_postfit.root" % fit_dir))
    
  #Clean up
  call("rm card*.txt", shell = True)
  i = 0
  f = open("./fit_%s/results_cat.txt" % (datacard_dir), "w")
  #Output fit results
  for fr in fit_results:
    print "RES: %d %.8f + %.8f - %.8f" % (i, fr[0], fr[1], fr[2])
    f.write("%d, %.8f, %.8f, %.8f\n" % (i, fr[0], fr[1], fr[2]))
    i += 1
  f.close()
  
if __name__ == "__main__":
  #Make fits for pseudodata and data
  make_fits(datacard_dir = "output_pseudodata_summer_June6")
  make_fits(datacard_dir = "output_data_summer_June6")
