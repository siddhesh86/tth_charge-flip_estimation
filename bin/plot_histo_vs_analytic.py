import numpy as np
import math
import sys
import os
from plot_pulls import bin_names_composite, bin_names_single, readMisIdRatios
from utils import read_category_ratios, get_bin_name



def make_category_matrix(catRatios, weighted = True):
  b = np.array(catRatios)
  if weighted:
    return (b[:,0], b[:,1])
  else:
    return (b[:,0], b[:,1]/b[:,1])

def make_header():
    latex = ""
    latex += "\documentclass{beamer}\n"
    latex += "\usepackage[utf8]{inputenc}\n"
    latex += "\usepackage{default}\n"
    latex += "\usepackage{graphbox}\n"
    latex += "\\begin{document}\n"
    latex += "\\beamertemplatenavigationsymbolsempty\n"
    #\usepackage[T1]{fontenc}
    return latex

   

def make_footer():
    return """
    \end{document}
    """


def make_latex(indir, ratios_hist, ratios_analytic):
    latex = make_header()
    for b in range(21):
        textcolor = "black"
        if ratios_analytic[b][0] + ratios_analytic[b][1] < ratios_hist[b][0] - ratios_hist[b][1] or \
           ratios_analytic[b][0] - ratios_analytic[b][1] > ratios_hist[b][0] + ratios_hist[b][1]:
            textcolor = "red"
                
    
        latex += "\\begin{frame}{%s}\n" % get_bin_name(b).replace("_", " ")
        latex += "\\begin{columns}[T,onlytextwidth]\n"
        latex += "\column{0.55\\textwidth}\n"
        latex += "With analytic background:\\\\"
        latex += "Pass: \includegraphics[align=c,width=0.8\\textwidth]{%s/bin%d/pass_fit_s_hybrid.png}\\\\ \n" % (indir, b)
        latex += "$ \\textcolor{%s}{%.5f \pm %.5f} $  \\\\ \n" % (textcolor, ratios_analytic[b][0], ratios_analytic[b][1])
        latex += "Fail: \includegraphics[align=c,width=0.8\\textwidth]{%s/bin%d/fail_fit_s_hybrid.png}\\\\ \n" % (indir, b)
        latex += "\column{0.45\\textwidth}\n"
        latex += "With histograms:\\\\"
        latex += "\includegraphics[width=\\textwidth]{%s/bin%d/pass_fit_s.png}\\\\ \n" % (indir, b)
        latex += "$ \\textcolor{%s}{%.5f \pm %.5f} $ \n" % (textcolor, ratios_hist[b][0], ratios_hist[b][1])
        latex += "\includegraphics[width=\\textwidth]{%s/bin%d/fail_fit_s.png}\\\\ \n" % (indir, b)
        latex += "\end{columns}\n"
        latex += "\end{frame}\n"
    latex += make_footer()
    return latex
    


if __name__ == "__main__":
  indir = "fit_output_pseudodata_eleESER_mva_0_6_notrig"
  #"fit_output_pseudodata_eleESER_mva_0_6_notrig/results_cat_shapes.txt", 
  #"fit_output_data_eleESER_mva_0_6_notrig/results_cat_shapes.txt"]:
  outfile = "latexed_hybrid_pseudodata.tex"
  f = open(outfile, "w")
  ratios_hist = read_category_ratios(indir + "/results_cat.txt")
  ratios_analytic = read_category_ratios(indir + "/results_cat_hybrid.txt")
  f.write(make_latex(indir, ratios_hist, ratios_analytic))
  
  
  """print file_cats    
    categoryRatios = read_category_ratios(file_cats, exclude_bins)
    #print categoryRatios  
    #x = solve_matrix(categoryRatios)
    #print_solution(x)
    #print "_"*80
    calculate(categoryRatios, exclude_bins)
    #make_gen_check()"""
