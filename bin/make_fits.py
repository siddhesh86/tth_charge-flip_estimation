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


def read_fit_result(fit_file, postfit_file, bin):
  ROOT.gROOT.SetBatch(True)
  f2 = TFile(postfit_file)
  fail_h = f2.Get("fail_prefit/DY")
  pass_h = f2.Get("pass_prefit/DY")

  f = TFile(fit_file)
  tree = f.Get("tree_fit_sb")
  
  tree.Draw("fit_status>>histStatus")
  fit_status = ROOT.gDirectory.Get("histStatus").GetMean()
  if fit_status != 0:
    raise AssertionError("Input fit %d did not converge! Fit status %d." % (bin, fit_status))
  

  tree.Draw("mu>>hist")
  hist = ROOT.gDirectory.Get("hist")
  tree.Draw("muErr>>histErr")
  histErr = ROOT.gDirectory.Get("histErr")
  tree.Draw("muLoErr>>histLo")
  #tree.Draw("muErr>>histLo")
  histLo = ROOT.gDirectory.Get("histLo")
  tree.Draw("muHiErr>>histHi")
  #tree.Draw("muErr>>histHi")
  histHi = ROOT.gDirectory.Get("histHi")
  mu = hist.GetMean()
  muErr = histErr.GetMean()
  muLoErr = histLo.GetMean()
  muHiErr = histHi.GetMean()
  #print mu, muErr, muLoErr, muHiErr, fit_file, postfit_file
  try:
    #postFit distributions are scaled to scale factor 1, need to multiply by fitted number
    bestFit = mu * pass_h.Integral() / (fail_h.Integral() + pass_h.Integral())
  except AttributeError:
    if fail_h.Integral() > 0:
      bestFit = 0.
    else: raise
  
  #print mu, muErr, muHiErr, muLoErr, bestFit
  if abs(muHiErr / muLoErr) > 2 or abs(muHiErr / muLoErr < 1./2):   #In this case probably failed to find crossing, use symmetric error
    #All cases manually verified
    #print fit_file
    if False:
      pass    
      #if ((bin in [7]) and "_pseudodata" in fit_file and "summer" in fit_file): #\
      #or ((bin in [7]) and "_data" in fit_file and "summer" in fit_file):
      #  muHiErr = muErr
      #  muLoErr = muErr
    else:
      raise AssertionError("Strange asymmetric errors! Probably failed to find crossing.")
  #fitHiErr = max(muHiErr, muErr)/mu * bestFit
  #fitLoErr = max(muLoErr, muErr)/mu * bestFit
  
  fitHiErr = muHiErr/mu * bestFit
  fitLoErr = muLoErr/mu * bestFit
  #if fitHiErr > 1:
  #  fitHiErr = max(1., mu)
  #if fitLoErr > 1:
  #  fitLoErr = max(1., mu)
  #Use Poisson mean of getting 0 events for uncertainty if no observed events
  if fitHiErr == 0. :
    lambda_poisson0 = -math.log(0.32)
    fitHiErr = lambda_poisson0 / (lambda_poisson0 + fail_h.Integral())
    fitLoErr = fitHiErr
  return (bestFit, fitHiErr, fitLoErr)


def failed_result():
  return (float('nan'), float('nan'), float('nan'))

#TODO: parallelize this for faster fitting
def make_fits(datacard_dir, skip_bins = []):
  dc_dir = "%s/cards" % datacard_dir
  ws_dir = "%s/workspaces" % datacard_dir
  mkdir_p(dc_dir)
  mkdir_p(ws_dir)
  fit_results = []
  for bin in range(21):
    this_card = "%s/card_%d.txt" % (dc_dir, bin)
    this_ws = "%s/workspace_%d.root" %(ws_dir, bin)
    print "combineCards.py pass=%s/SScards/htt_SS_%d_13TeV.txt fail=%s/cards/OScards/htt_OS_%d_13TeV.txt > %s" % (dc_dir, bin, datacard_dir, bin, this_card)
    
    #1. step: combine SS and OS datacatds    
    if bin in skip_bins:
      fit_results.append(failed_result())
      continue
    else:
      call("combineCards.py pass=%s/SScards/htt_SS_%d_13TeV.txt fail=%s/cards/OScards/htt_OS_%d_13TeV.txt > %s" % (dc_dir, bin, datacard_dir, bin, this_card), shell = True)
    
    #Hack to prevent PostFitShapesFromWorkspace messing up directory structure - copy datacard to current directory
    call("cp %s ." % (this_card), shell = True)
    print "text2workspace.py %s -o %s -P HiggsAnalysis.CombinedLimit.TagAndProbeModel:tagAndProbe" % (this_card, this_ws)
    #2. Make Roofit workspace from datacard
    call("text2workspace.py %s -o %s -P HiggsAnalysis.CombinedLimit.TagAndProbeModel:tagAndProbe" % (this_card, this_ws), shell = True)
    #Specify output directory for fit results - TODO: move out of bin-directory 
    fit_dir = "./fit_%s/bin%d" % (datacard_dir, bin)
    mkdir_p(fit_dir)

    #3. Perform fit with combine
    #Default fit settings specified in else clause
    #But always do not give convergence - settings for specific fits defined here, adjust as necessary:
    specific_settings = "--robustFit 1 "
    if ((bin in [2, 5, 7, 9, 15, 16, 18]) and "_pseudodata" in fit_dir and "summer" in fit_dir):
        specific_settings += "--minimizerStrategy 0"
    elif ((bin in [15]) and "_data" in fit_dir and "summer" in fit_dir):
        specific_settings = "--minimizerStrategy 0"
    elif ((bin in [1, 4, 7, 8, 14, 17, 18]) and "_data" in fit_dir and "summer" in fit_dir):
        specific_settings += "--minimizerStrategy 0"
        #elif ((bin in [1]) and "_pseudodata" in fit_dir and "summer" in fit_dir):
        #specific_settings = "--minimizerStrategy 2"
    elif ((bin in [10, 13]) and "_data" in fit_dir and "summer" in fit_dir):
        specific_settings += "--minimizerAlgo Minuit"
    else:
        specific_settings += ""
        
    command = "combine -v0 -M MaxLikelihoodFit %s --out %s --plots --saveNormalizations --skipBOnlyFit --saveShapes --saveWithUncertainties --maxFailedSteps 20 %s" % (this_ws, fit_dir, specific_settings)
    print command    
    #Call combine
    call(command, shell = True)

    #4. Create postfit plots
    #print "PostFitShapesFromWorkspace -d %s -w %s -o %s/output_postfit.root -f %s/mlfit.root:fit_s --postfit --sampling --print" % ("card_%d.txt" % (bin), this_ws, fit_dir, fit_dir)
    call("PostFitShapesFromWorkspace -d %s -w %s -o %s/output_postfit.root -f %s/mlfit.root:fit_s --postfit --sampling" % ("card_%d.txt" % (bin), this_ws, fit_dir, fit_dir), shell=True)
    
    #5. Add to list of fit results
    fit_results.append(read_fit_result("%s/mlfit.root" % fit_dir, "%s/output_postfit.root" % fit_dir, bin = bin))
    
    
  #Clean up
  call("rm card*.txt", shell = True)
  f = open("./fit_%s/results_cat.txt" % (datacard_dir), "w")
  #Output fit results
  for i, fr in enumerate(fit_results):
    print "RES: %d %.8f + %.8f - %.8f" % (i, fr[0], fr[1], fr[2])
    f.write("%d, %.8f, %.8f, %.8f\n" % (i, fr[0], fr[1], fr[2]))
  f.close()
  
if __name__ == "__main__":
  skip_bins = [
    6, # data totally empty for this bin, cannot fit
    7, # pseduodata fit does not converge
  ]
  #Make fits for pseudodata and data
  print "Pseudodata:"
  make_fits(datacard_dir = "output_pseudodata_summer_Aug25", skip_bins = skip_bins)
  print "Data:"
  make_fits(datacard_dir = "output_data_summer_Aug25", skip_bins = skip_bins)
