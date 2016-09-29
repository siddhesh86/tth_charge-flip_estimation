from ROOT import TH2F, TFile
from array import array

bin_names_composite = ["BB_LL", "BB_ML", "BB_MM", "BB_HL", "BB_HM", "BB_HH",
      "EE_LL", "EE_ML", "EE_MM", "EE_HL", "EE_HM", "EE_HH",
      "BE_LL", "BE_ML", "EB_ML", "BE_MM", "BE_HL", "EB_HL",
      "BE_HM", "EB_HM", "BE_HH"]
bin_names_single = ["BL", "BM", "BH", "EL", "EM", "EH"]




def read_category_ratios(file_cats, exclude_bins = []):
  f = open(file_cats)  
  bins = []
  ratios = []
  for line in f.readlines():
    spl = line.split(",")
    bin_nr = int(spl[0])
    #Skip some bins
    if bin_nr in exclude_bins: continue
    if float(spl[2]) > 0:
      ratios.append((float(spl[1]), float(spl[2])))
    else:
      ratios.append((float(spl[1]), 0.01))
  return ratios
  
def get_bin_name(bin_nr):
  return bin_names_composite[bin_nr]
    
def bin_names_to_numbers(bin_names):
  bin_nrs = []
  for b in bin_names:
      bin_nrs.append(bin_names_composite.index(b))
  return bin_nrs
  

def fit_results_to_file(rates, uncs, fittype, fitname, datastring):
    fname = "fit_output_%s_%s/fit_res%s.root" % (datastring, fitname, fittype)
    f = TFile(fname,"recreate")
    f.cd()

    binsPt = [15, 25, 50, 1000]
    binsEta = [0, 1.479, 2.5]
    NbinsPt = len(binsPt) - 1
    NbinsEta = len(binsEta) - 1

    h = TH2F("chargeMisId","chargeMisId;p_{T}(e) [GeV];#eta(e)", NbinsPt, array('d',binsPt), NbinsEta, array('d',binsEta))
    for i in range(len(rates)):
        #print "%d %f %f %f" %(i, (i%NbinsPt)+1, (i/NbinsPt)+1, rates[i])
        h.SetBinContent((i % NbinsPt)+1, (i / NbinsPt)+1, rates[i])
        h.SetBinError((i % NbinsPt)+1, (i / NbinsPt)+1, uncs[i])

    f.Write()
    f.Close()
