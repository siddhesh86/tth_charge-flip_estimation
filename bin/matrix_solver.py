import numpy as np
import math
from plot_pulls import bin_names_composite, bin_names_single, readMisIdRatios
from utils import read_category_ratios, bin_names_to_numbers, fit_results_to_file


def make_coefficient_matrix(exclude_bins=[]):
 coeffs = [[ 2,  0,  0,  0,  0,  0,],
 [ 1,  1,  0,  0,  0,  0],
 [ 0,  2,  0,  0,  0,  0],
 [ 1,  0,  1,  0,  0,  0],
 [ 0,  1,  1,  0,  0,  0],
 [ 0,  0,  2,  0,  0,  0],
 [ 0,  0,  0,  2,  0,  0],
 [ 0,  0,  0,  1,  1,  0],
 [ 0,  0,  0,  0,  2,  0],
 [ 0,  0,  0,  1,  0,  1],
 [ 0,  0,  0,  0,  1,  1],
 [ 0,  0,  0,  0,  0,  2],
 [ 1,  0,  0,  1,  0,  0],
 [ 0,  1,  0,  1,  0,  0],
 [ 1,  0,  0,  0,  1,  0],
 [ 0,  1,  0,  0,  1,  0],
 [ 0,  0,  1,  1,  0,  0],
 [ 1,  0,  0,  0,  0,  1],
 [ 0,  0,  1,  0,  1,  0],
 [ 0,  1,  0,  0,  0,  1],
 [ 0,  0,  1,  0,  0,  1]]
 used_coeffs = []
 for i in range(len(coeffs)):
   if not i in exclude_bins:
     used_coeffs.append(coeffs[i]) 
 return np.array(used_coeffs)

def gen_level_probs():
  return [(0.000035157, 0.000033635), (0.000221524, 0.000014054), (0.000267707, 0.000035297), (0.000829932, 0.000155612), (0.001864478, 0.000074103), (0.002657734, 0.000203647)]


def make_cat_ratios_from_gen():
  gen_probs = gen_level_probs()
  indices = [[0,0], [0,1], [1,1], [2,0], [2,1], [2,2], [3,3], [4,3], [4,4], [5,3], [5,4], [5,5],
      [0,3], [1,3], [4,0], [1,4], [2,3], [5,0], [2,4], [5,1], [2,5]]
  cat_ratios = []
  for r in range(len(indices)):
    #print gen_probs[indices[r][0]]
    val = gen_probs[indices[r][0]][0]+gen_probs[indices[r][1]][0]
    err = math.sqrt(gen_probs[indices[r][0]][1]**2+gen_probs[indices[r][1]][1]**2)
    #if err == 0: err = 0.01
    print r, val, err
    cat_ratios.append((val, err))
  print cat_ratios
  return cat_ratios

def make_category_matrix(catRatios, weighted = True):
  b = np.array(catRatios)
  if weighted:
    return (b[:,0], b[:,1])
  else:
    return (b[:,0], b[:,1]/b[:,1])

def solve_matrix(catRatios):
  A = make_coefficient_matrix()
  #print A
  (b, W) = make_category_matrix(catRatios)
  print W
  #print 1/W
  W = np.sqrt(np.diag(1/W))
  #print W
  Aw = np.dot(W,A)
  Bw = np.dot(b,W)
  #print W
  x_lstsq = np.linalg.lstsq(Aw,Bw)[0] # computing the numpy solution

  Q,R = np.linalg.qr(Aw) # qr decomposition of A
  Qb = np.dot(Q.T,Bw) # computing Q^T*b (project b onto the range of A)
  x_qr = np.linalg.solve(R,Qb) # solving R*x = Q^T*b

  print 'qr solution'
  print x_qr
  print 'lstqs solution'
  print x_lstsq
  return x_qr
  
def print_solution(x):
  for x, value in np.ndenumerate(x):
    print x,value*100

def make_gen_check():
  x = solve_matrix(make_cat_ratios_from_gen())
  print_solution(x)


def calculate_M(Aw):
  #print W
  M = np.dot(np.linalg.inv( np.dot(np.transpose(Aw), Aw) ) , np.transpose(Aw))
  #print M
  return M

def calculate_rates(M, b):
  return np.dot(M, b)

def calculate_uncertainties(M, deltab):
  uncs = [] 
  for i in range(6):    
    uncs.append(0)
    for k in range(len(M)):
      uncs[i] += (M[i, k] * deltab[k])**2
    uncs[i] = math.sqrt(uncs[i])
  return np.array(uncs)

def print_solution_with_uncertainties(x, uncs):
  for i in range(len(x)):
    print x[i]*100, "+-", uncs[i]*100
    
def print_latex_header():
   print """\\begingroup\setlength{\\fboxsep}{0pt}
    \colorbox{cyan}{%
    \\begin{tabular}{llccc}
	    \hline
	    ID & & $10\leq\pt<25\GeV$ & $25\leq\pt<50\GeV$ & $50\GeV\leq\pt$ \\\\
	    \hline	    
    \end{tabular}%
}\endgroup"""
    
def print_solution_latex(x, uncs, fittype, datastring):
  latex = """
	    \multirow{2}{*}{%s}   & $0\leq\eta<1.479$    & %.4f $\pm$ %.4f & %.4f $\pm$ %.4f & %.4f $\pm$ %.4f  \\\\
	                            & $1.479\leq\eta<2.5$  & %.4f $\pm$ %.4f & %.4f $\pm$ %.4f & %.4f $\pm$ %.4f  \\\\
	    \hline""" % (datastring, x[0]*100, uncs[0]*100, x[1]*100, uncs[1]*100, x[2]*100, uncs[2]*100, x[3]*100, uncs[3]*100, x[4]*100, uncs[4]*100, x[5]*100, uncs[5]*100)
  print latex  

def calculate(catRatios, exclude_bins = [], weighted = True):
  A = make_coefficient_matrix(exclude_bins)
  (b, W) = make_category_matrix(catRatios, weighted)
  w = np.diag(1/W)
  Aw = np.dot(w, A)  
  M = calculate_M(Aw)
  bw = np.dot(w, b)
  rates = calculate_rates(M, bw)
  (b, err) = make_category_matrix(catRatios)
  uncs = calculate_uncertainties(M, np.dot(w,err))
  return (rates.tolist(), uncs.tolist())
  

if __name__ == "__main__":
  exclude_bin_names = [#"BB_LL", "BB_ML", "BB_HL", "BB_HH", "EE_LL", "EE_ML", "EE_HL", "EE_HH", "BE_HL", "BE_HH",
    #"EE_HM", "BE_HM", "EB_HM",
    #"BE_LL", "BE_ML", "EB_ML"
    ]
  exclude_bins = bin_names_to_numbers(exclude_bin_names)
  for datastring in ["data", "pseudodata"]:
    for fitname in ["eleESER_mva_0_6_notrig", "eleESER2"]:
      for fittype in ["", "shapes", "hybrid"]:
        fittypestring = fittype
        if len(fittype) > 0: fittypestring = "_"+fittype
        file_cats = "fit_output_%s_%s/results_cat%s.txt" % (datastring, fitname, fittypestring)
        print file_cats    
        categoryRatios = read_category_ratios(file_cats, exclude_bins)
        #print categoryRatios  
        #x = solve_matrix(categoryRatios)
        #print_solution(x)
        #print "_"*80
        (rates, uncs) = calculate(categoryRatios, exclude_bins)
        #print_solution_with_uncertainties(rates, uncs)
        #make_gen_check()
        fit_results_to_file(rates, uncs, fittypestring, fitname, datastring)
      
        print_solution_latex(rates, uncs, fittype, datastring)   
  print_latex_header()     
