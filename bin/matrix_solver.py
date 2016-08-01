import numpy as np
import math
from plot_pulls import bin_names_composite, bin_names_single, readMisIdRatios, readCategoryRatios


def make_coefficient_matrix():
  return np.array([[ 2,  0,  0,  0,  0,  0,],
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
 [ 0,  0,  1,  0,  0,  1]])

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

def make_category_matrix(catRatios):
  b = np.array(catRatios)
  return (b[:,0], b[:,1])

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

if __name__ == "__main__":
  for file_cats in ["fit_output_data_oldDY_notrig/results_cat.txt", "fit_output_data_newDY_notrig/results_cat.txt", "fit_output_data_newDY/results_cat.txt",
        "fit_output_pseudodata_oldDY_notrig/results_cat.txt", "fit_output_pseudodata_newDY_notrig/results_cat.txt", "fit_output_pseudodata_newDY/results_cat.txt"]:
    print file_cats
    categoryRatios = readCategoryRatios(file_cats)
    print categoryRatios  
    x = solve_matrix(categoryRatios)
    print_solution(x)
    #print "_"*80
    #make_gen_check()
