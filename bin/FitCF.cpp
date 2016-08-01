#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <map>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include <TMinuit.h>
#include <Math/Functor.h>
#include <Fit/Fitter.h>

using namespace std;

#define N_CATEGORIES 21

struct binSt{
    binSt() {};
    binSt(vector<float> pt, vector<float> eta) {
        binsPt = pt;
        binsEta = eta;

        nBPt = binsPt.size();
        nBEta = binsEta.size();
    }

    vector<float> binsPt;
    vector<float> binsEta;

    int nBPt;
    int nBEta;

    int getNProb() {return nBPt*nBEta;}
    int getNBPt() {return nBPt;}
    int getNBEta() {return nBEta;}

};

void getSingleEleBinNumbers(int composite_bin, int &bin1, int &bin2){
    /* From the composite bins:
      {0,"BB_LL"},{1,"BB_ML"},{2,"BB_MM"},{3,"BB_HL"},{4,"BB_HM"},{5,"BB_HH"},
      {6,"EE_LL"},{7,"EE_ML"},{8,"EE_MM"},{9,"EE_HL"},{10,"EE_HM"},{11,"EE_HH"},
      {12,"BE_LL"},{13,"BE_ML"},{14,"EB_ML"},{15,"BE_MM"},{16,"BE_HL"},{17,"EB_HL"},
      {18,"BE_HM"},{19,"EB_HM"},{20,"BE_HH"},};
      extract the 2 single bins: 
      {0, "BL"}, {1, "BM"}, {2, "BH"}, {3, "EL"}, {4, "EM"}, {5, "EH}
    */
    if(composite_bin == 0 || composite_bin == 1 || composite_bin == 3 || composite_bin == 12 || composite_bin == 14 || composite_bin == 17)
      bin1 = 0;
    else if(composite_bin == 2 || composite_bin == 4 || composite_bin == 13 || composite_bin == 15 || composite_bin == 19)
      bin1 = 1;
    else if(composite_bin == 5 || composite_bin == 16 || composite_bin == 18 || composite_bin == 20)
      bin1 = 2;
    else if(composite_bin == 6 || composite_bin == 7 || composite_bin == 9)
      bin1 = 3;
    else if(composite_bin == 8 || composite_bin == 10)
      bin1 = 4;
    else if(composite_bin == 11)
      bin1 = 5;

    if(composite_bin == 0)
      bin2 = 0;
    else if(composite_bin == 1 || composite_bin == 2)
      bin2 = 1;
    else if(composite_bin <= 5)
      bin2 = 2;
    else if(composite_bin == 6 || composite_bin == 12 || composite_bin == 13 || composite_bin == 16)
      bin2 = 3;
    else if(composite_bin == 7 || composite_bin == 8 || composite_bin == 14 || composite_bin == 15 || composite_bin == 18)
      bin2 = 4;
    else if(composite_bin == 9 || composite_bin == 10 || composite_bin == 11 || composite_bin == 17 || composite_bin == 19 || composite_bin == 20)
      bin2 = 5;
}


// function Object to be minimized
struct Chi2 {
    vector<std::pair<vector<float>, vector<int> > > _vals;

    void setPoint(float val, float eval, float p1, float p2) {
        vector<float> vals(2,0);
        vector<int> bins(2,0);

        vals[0] = val;
        vals[1] = eval;
        bins[0] = p1;
        bins[1] = p2;

        std::pair<vector<float>, vector<int> > p(vals,bins);
        _vals.push_back(p);
    }

    // implementation of the function to be minimized
    double operator() (const double * param) {
        double chi2 = 0;

        float val,eval;
        int p1,p2;
        for(unsigned int ip = 0;ip<_vals.size();ip++) {
            val = _vals[ip].first[0];
            eval = _vals[ip].first[1];
            p1 = _vals[ip].second[0];
            p2 = _vals[ip].second[1];
            //if(eval == 0) eval = val;
            //if(eval == 0) eval = 0.01;  //If also fitted value was 0
            chi2 += pow( val-(param[p1]+param[p2]), 2)/pow(eval,2);
        }
        //cout << " chi2: " << chi2 << endl;
        return chi2;
    }
};


int read_results(string infile, Chi2 &chi2){
  FILE *fp;
  fp = fopen(infile.data(), "r");
  //char str1[10], str2[10], str3[10], str4[10];
  int bin = -1;
  double mean = 0, up_error = 0, down_error = 0;
  if(NULL == fp)
  {
      printf("\nError in opening file.");
      return -1;
  }

  int i = 0;
  while (EOF != fscanf(fp, "%d, %lf, %lf, %lf", &bin, &mean, &up_error, &down_error)){
    printf("%d %f %f %f\n", bin, mean, up_error, down_error);
    int p1, p2;
    getSingleEleBinNumbers(bin, p1, p2);
    chi2.setPoint( mean, max(up_error, down_error), p1, p2);
    i++;
  }  
  

  //int i = 0;
  /*while(EOF != fscanf(fp, " %[^,], %[^,], %[^,], %d, %f, %f, %f ", bins[i], means[i], up_errors[i], down_errors[i]))
  {
      printf("\n%d, %f, %f, %f ", bins[i], means[i], up_errors[i], down_errors[i]));
      i++;
  }*/
  fclose(fp);
  return 0;
}



binSt setBinSizes(string file) {
    vector<float> binsPt;
    vector<float> binsEta;

    map<string, float> yields;
    map<string, float> eyields;
    map<string, vector<float> > bins;

    binsPt.push_back(10);
    binsPt.push_back(25);
    binsPt.push_back(50);
    //binsPt.push_back(...)
    binsEta.push_back(0);
    binsEta.push_back(1.479);
    //binsEta.push_back(2.5);
    
    binSt binstruct(binsPt, binsEta);
    return binstruct;
}



int main(int argc, char* argv[]) {
    string file;
    bool isData = true;
    string singleCateg = "";
    
    char c;
    while( (c = getopt(argc, argv, "f:d:s:D:a:n:h")) != -1 ){
        switch (c) {
        case 'f': { file = string(optarg); break;}
        case 's': { singleCateg = string(optarg); break;}
        case 'D': { isData = bool(atoi(optarg)); break;}
        default : {
            cout << "configuration options:\n "
                 << "-f : file to read (root or ASCII) \n "
                 << "-s <categ> perform a fit over a single Z category.\n "
                 << "-D run on data (0 per default). \n "
                 <<"-h help \n" << endl;
            return 0; }
        }
    }

    if (isData) ; //temp

    /*int catbins[N_CATEGORIES]; 
    double means[N_CATEGORIES];
    double up_errors[N_CATEGORIES];
    double down_errors[N_CATEGORIES];
    */
    
    //==============================================

    Chi2 chi2;
    binSt bins = setBinSizes(file);
    string which_fit = "pseudodata_testnewconf";
    read_results("fit_output_"+which_fit+"/results_cat.txt", chi2);


    // perform the final fit ====================
    int nvars = bins.getNProb();

    ROOT::Fit::Fitter  fitter;
    ROOT::Math::Functor fcn(chi2,nvars);

    // bloody ROOT and lack of vector handling
    double* vars = new double[nvars];
    fitter.SetFCN(fcn,vars);

    // set step sizes and limits
    for (int i=0; i<nvars; ++i) {
        fitter.Config().ParSettings(i).SetStepSize(0.000001);
        fitter.Config().ParSettings(i).SetLimits(0,0.2);
    }

    bool ok = fitter.FitFCN();
    if (!ok) {
        cout << "The final fit did not converge properly!" << endl;
        return 1;
    }

    fitter.CalculateMinosErrors();
    fitter.CalculateHessErrors();
    const ROOT::Fit::FitResult & result = fitter.Result();
    result.Print(std::cout);

    const double * parFit = result.GetParams();
    const double * parErrs = result.GetErrors();
    cout << "probabilities (in percent) (i/etabin/ptbin) and fit result" << endl;
    for(int i=0; i<nvars; i++)
        cout << i << "    "
             << i/bins.nBPt << "   "
             << i%bins.nBPt << " ==> "
             << parFit[i] * 100 << " +- " << parErrs[i] * 100 << endl;

    //size_t p = file.find(".root");

    string tag = "fit_output_"+which_fit+"/fit_res";
    
    string fname = tag+(isData?"_data":"_MC")+".root";
    TFile* f = new TFile(fname.c_str(),"recreate");
    f->cd();

    // adding the last bin boundaries, bloody root...
    double* binsPt = new double[ bins.nBPt+1 ];
    double* binsEta = new double[ bins.nBEta+1 ];

    for(int ib = 0; ib<max( bins.nBPt, bins.nBEta ); ib++) {
        if(ib<bins.nBPt) binsPt[ib] = bins.binsPt[ib];
        if(ib<bins.nBEta) binsEta[ib] = bins.binsEta[ib];
    }
    binsPt[bins.nBPt] = 1000;
    binsEta[bins.nBEta] = 2.5;

    TH2F* h = new TH2F("chargeMisId","chargeMisId;p_{T}(e) [GeV];#eta(e)", bins.nBPt, binsPt, bins.nBEta, binsEta);
    for(int i=0; i<nvars; i++){
        cout << i << " " << bins.nBPt << " " << (i/bins.nBPt)+1 << " " << (i%bins.nBPt)+1 << " " << parFit[i] << endl;
        h->SetBinContent((i%bins.nBPt)+1, (i/bins.nBPt)+1, parFit[i]);
        h->SetBinError((i%bins.nBPt)+1, (i/bins.nBPt)+1, parErrs[i]);
    }

    f->Write();
    f->Close();
    
}
