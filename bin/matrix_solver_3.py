import numpy as np
import math
from plot_pulls import bin_names_composite, bin_names_single
from utils import bin_names_to_numbers, fit_results_to_file


def make_coefficient_matrix(exclude_bins=[]):
 coeffs = [[ 2,  0],
 [ 1,  1],
 [ 0,  2]]
 used_coeffs = []
 for i in range(len(coeffs)):
   if not i in exclude_bins:
     used_coeffs.append(coeffs[i]) 
 return np.array(used_coeffs)

#def gen_level_probs():
#  return [(0.000035157, 0.000033635), (0.000221524, 0.000014054), (0.000267707, 0.000035297), (0.000829932, 0.000155612), (0.001864478, 0.000074103), (0.002657734, 0.000203647)]


def make_cat_ratios_from_gen():
  cat_ratios = [(0.0049654093207578565, 0.0001462259953501589),
        (0.0023022851031497876, 0.000041019420446535004), 
        (0.0004516881846517407, 0.000012535791157441956)]
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
  ##print 1/W
  W = np.sqrt(np.diag(1/W**2))
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
  for i in range(2):    
    uncs.append(0)
    print "LENM", len(M)
    print M
    for k in range(len(M[0])):
      uncs[i] += (M[i, k] * deltab[k]) ** 2
      print i, uncs[i], M[i, k], deltab[k]
    uncs[i] = math.sqrt(uncs[i])
    print "unc", k, uncs[i]
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
    
def print_solution_latex(x, uncs, datastring):
  latex = """
	    \multirow{2}{*}{%s}   & $0\leq\eta<1.479$    & %.4f $\pm$ %.8f   \\\\
	                            & $1.479\leq\eta<2.5$  & %.4f $\pm$ %.8f   \\\\
	    \hline""" % (datastring, x[0]*100, uncs[0]*100, x[1]*100, uncs[1]*100)
  print latex

def print_ratios_latex(ratios, datastring):
  latex = """
	    \multirow{2}{*}{%s}   & $0\leq\eta<1.479$    & %.4f $\pm$ %.4f & %.4f $\pm$ %.4f & %.4f $\pm$ %.4f  \\\\
	                            & $1.479\leq\eta<2.5$  & %.4f $\pm$ %.4f & %.4f $\pm$ %.4f & %.4f $\pm$ %.4f  \\\\
	    \hline""" % (datastring, ratios[0][0]*100, ratios[0][1]*100, ratios[1][0]*100, ratios[1][1]*100, ratios[2][0]*100, ratios[2][1]*100, ratios[3][0]*100, ratios[3][1]*100, ratios[4][0]*100, ratios[4][1]*100, ratios[5][0]*100, ratios[5][1]*100)
  print latex


def calculate(catRatios, exclude_bins = [], weighted = True):
  A = make_coefficient_matrix(exclude_bins)
  (b, W) = make_category_matrix(catRatios, weighted)
  #print W
  w = np.diag(1/W**2)
  Aw = np.dot(w, A)  
  M = calculate_M(Aw)
  bw = np.dot(w, b)
  rates = calculate_rates(M, bw)
  (b, err) = make_category_matrix(catRatios)
  print "B", b
  print "bw", bw
  #print "err", err
  #print "w", w
  #print np.dot(w,err)
  print "M", M
  print "MW", np.dot(M,w)
  print "korrutis", np.dot(w,err)
  print b, np.sqrt(1/W)
  print "rates", rates
  #uncs = calculate_uncertainties(M, np.dot(w,err))
  uncs = calculate_uncertainties(np.dot(M,w), W)
  return (rates.tolist(), uncs.tolist())


def calculate_solution(categoryRatios, exclude_bins, fitname, fittypestring, datastring):
  (rates, uncs) = calculate(categoryRatios, exclude_bins)
  #fit_results_to_file(rates, uncs, fittypestring, fitname, datastring)
  print_solution_latex(rates, uncs, datastring)   

if __name__ == "__main__":
  exclude_bins = []
  #print file_cats    
  categoryRatios = make_cat_ratios_from_gen()
  calculate_solution(categoryRatios, [], "asd", "asd", "gencheck")
  #print categoryRatios  
  #x = solve_matrix(categoryRatios)
  #print_solution(x)
  #print "_"*80
  #print_solution_with_uncertainties(rates, uncs)
  #make_gen_check()
  #print_latex_header()     
