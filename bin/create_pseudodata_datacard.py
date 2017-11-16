from ROOT import TFile, TDirectory
import ROOT
import numpy as np
import math
from utils import bin_names_composite

"""@file docstring
Script for creating pseudodata and datacards, which are suitable as input to charge flip estimation

#See: <https://en.wikipedia.org/wiki/Least_squares#Weighted_least_squares>
for notation and background

@author Andres Tiko <andres.tiko@cern.ch>
"""

#Samples to use
samples = ["data_obs",
	    "DY",
      "DY_fake", 
	    #"WJets", 
	    "Singletop", 
	    "Diboson", 
	    "TTbar"
  ]

#21 categories of lepton pairs by pT and eta
categories = {}
categories["ele"] = bin_names_composite
#For muons all in the same cateogry
categories["mu"] = ["total"]
prefix = "ttH_charge_flip"
charges = ["OS", "SS"]

#Systematic uncertainties to add to datacard
systematics = {}
systematics["ele"] = ["CMS_ttHl_electronERUp",
       "CMS_ttHl_electronERDown",
       "CMS_ttHl_electronESEndcapUp",
       "CMS_ttHl_electronESEndcapDown",
       "CMS_ttHl_electronESBarrelUp",
       "CMS_ttHl_electronESBarrelDown"]
#Leave the muon uncertainties out for now for faster testing
systematics["mu"] = [#"CMS_ttHl_muonERUp",
       #"CMS_ttHl_muonERDown",
       #"CMS_ttHl_muonESEndcap1Up",
       #"CMS_ttHl_muonESEndcap1Down",
       #"CMS_ttHl_muonESEndcap2Up",
       #"CMS_ttHl_muonESEndcap2Down",
       #"CMS_ttHl_muonESBarrel1Up",
       #"CMS_ttHl_muonESBarrel1Down",
       #"CMS_ttHl_muonESBarrel2Up",
       #"CMS_ttHl_muonESBarrel2Down"
      ]

#function to create datacards
#Also includes some processing
def create_pseudodata(infile, outfile_data, outfile_pseudodata, channel, rebin=1):
    f = TFile(infile)
    fout = TFile(outfile_data, "RECREATE")
    fout_pseudo = TFile(outfile_pseudodata, "RECREATE")
    
    for charge in charges:
        for cat in categories[channel]:
            first = True
            dirname = "%s_%s_%s" % (prefix, charge, cat)
            cd = fout.mkdir(dirname)
            cd_pseudo = fout_pseudo.mkdir(dirname)
            data_histo = ''
            for sample in samples:
                cd.cd()
                #Read nominal histograms                
                #histo_nominal = f.Get("%s/x_%s"  % (dirname, sample))
                histoname = "%s/%s"  % (dirname, sample)
                print("%s"  % histoname)
                histo_nominal = ''
                #histo_nominal = f.Get("%s/%s"  % (dirname, sample))
                histo_nominal = f.Get(histoname)
                #print "Histo:", histo_nominal
                #if histo_nominal is None:
                if not histo_nominal :
                    print ("Histogram %s could not fetch" % histoname)
                    continue
                print("\t\t running for %s histogram" % histoname)
                histo_nominal.Rebin(rebin)
                #print "%s/x_%s"  % (dirname, sample)
                for b in range(1, histo_nominal.GetNbinsX()+1):
                  content = histo_nominal.GetBinContent(b)
                  #If bin less than zero, set to zero with uncertainty
                  if content < 0:
                    histo_nominal.SetBinContent(b, 0.0)
                    histo_nominal.SetBinError(b, math.sqrt(histo_nominal.GetBinError(b)**2 + histo_nominal.GetBinContent(b)) )
                histo_nominal.Write()
                
                #Don't add data to pseudodata
                if sample == "data_obs": continue

                cd_pseudo.cd()
                histo_nominal.Write()
                
                #Add different MCs together as pseudodata
                if first == True:
                    data_histo = histo_nominal.Clone()
                    first = False
                else:
                    data_histo.Add(histo_nominal)
                
                #Add systematics
                for syst in systematics[channel]:
                  if syst.startswith("CMS_ttHl_electronER") and not sample == "DY": continue
                  if syst.startswith("CMS_ttHl_muonER") and not sample == "DY": continue
                  #histo = f.Get("%s/x_%s_%s"  % (dirname, sample, syst))
                  histoname = "%s/%s_%s"  % (dirname, sample, syst)
                  histo = f.Get(histoname)
                  print("\t\t\t\t Systematic %s"  % histoname)
                  histo.Rebin(rebin)
                  cd.cd() 
                  #print histo, dirname, sample, syst
                  for b in range(1, histo.GetNbinsX()+1):
                    content = histo.GetBinContent(b)
                    if content < 0:
                      histo.SetBinContent(b, 0.0)
                      histo.SetBinError(b, math.sqrt(histo.GetBinError(b)**2 + histo.GetBinContent(b)) )
                    
                  histo.Write()
                  cd_pseudo.cd() 
                  histo.Write()

            if not data_histo:
                print("data_histo is null.  'continue'")
                continue
            #data_histo.SetNameTitle("x_data_obs", "x_data_obs")
            data_histo.SetNameTitle("data_obs", "data_obs")
            #Generate poisson yields from MC expectation for pseudodata
            for b in range(1, data_histo.GetNbinsX()+1):
                #print data_histo.GetBinContent(b)
                data_histo.SetBinContent(b, np.random.poisson(data_histo.GetBinContent(b)))
                data_histo.Sumw2(ROOT.kFALSE)
                #print data_histo.GetBinContent(b), data_histo.GetBinError(b)
                #data_histo.SetBinError(b, max(1, math.sqrt(data_histo.GetBinContent(b))))
            data_histo.Write()            
    f.Close()
    fout.Close()
    fout_pseudo.Close()
            

if __name__ == "__main__":
  np.random.seed(123)
  #indir = "/home/andres/ttHAnalysis/2016/histosCF_summer_Aug25_noMassScaling/datacards/charge_flip/"
  indir = "/home/ssawant/ttHAnalysis/2016/histosCF_summer_June6/datacards/charge_flip/"
  infile = "prepareDatacards_charge_flip_mass_ll.root"
  datafile = "prepareDatacards_data_charge_flip_mass_ll.root"
  pseudodatafile = "prepareDatacards_pseudodata_charge_flip_mass_ll.root"
  
  #Create datacards for electrons
  create_pseudodata(indir + infile, 
    indir + datafile, 
    indir + pseudodatafile, 
    "ele",
    rebin = 2    
  )
  
  #Also for muons - enable as needed  
    
  #indir = "/home/andres/ttHAnalysis/2016/histosCF_mu_summer/datacards/charge_flip/"
  #create_pseudodata(indir + infile, 
  #  indir + datafile, 
  #  indir + pseudodatafile, "mu")

