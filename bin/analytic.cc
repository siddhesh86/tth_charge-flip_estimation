#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <map>
#include <utility>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include <TGraphAsymmErrors.h>
//#include <TLorentzVector.h>
#include <TMinuit.h>
#include "TKey.h"

#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooArgSet.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooDataHist.h>
#include <RooSimultaneous.h>
#include <RooNumConvPdf.h>
#include <RooMsgService.h>

#include <RooBreitWigner.h>
#include <RooCBShape.h>
#include <RooExponential.h>

#include <Math/Functor.h>
#include <Fit/Fitter.h>



using namespace std;
using namespace RooFit;

typedef vector<string> VString;


float getErr(float m1, float m2, float e1, float e2) {
  /*float e=pow(e1/m1,2)+pow(e2/m2,2);
  if(m1==0) return 1;
  return (m1/m2)*sqrt(e);*/
  float e=pow(e1/m1,2)+pow((e1+e2)/(m1+m2),2);
  if(m1==0) return 1;
  return (m1/(m1+m2))*sqrt(e);
}

RooNumConvPdf* shapeZ(string tag, RooRealVar* x) {

  RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::DataHandling);

  RooRealVar* mZ0=new RooRealVar( ("m_Z0_"+tag).c_str(),"Z0 mass", 91.188, "GeV/c^{2}" );
  RooRealVar* gammaZ0=new RooRealVar( ("gamma_Z0_"+tag).c_str(), "Z0 width",2.4952, "GeV/c^{2}" );
  RooBreitWigner* bw0=new RooBreitWigner( ("bw0"+tag).c_str(),"true BW",*x, *mZ0, *gammaZ0);

  RooRealVar* cb_bias=new RooRealVar( ("cbb_"+tag).c_str(), "bias",0.07, -3.0, 3.0 );
  RooRealVar* cb_width=new RooRealVar( ("cbw_"+tag).c_str(),"width", 1.,0.0,5 );
  RooRealVar* cb_alpha=new RooRealVar( ("cba_"+tag).c_str(),"alpha", 1.2,0.03,2.0 );
  RooRealVar* cb_power=new RooRealVar( ("cbn_"+tag).c_str(),"power", 5 );

  RooCBShape* cb_pdf=new RooCBShape( ("cb_pdf_"+tag).c_str(), "CB shape",
				     *x,*cb_bias, *cb_width, *cb_alpha, *cb_power );

  RooNumConvPdf* bw=new RooNumConvPdf( ("bw_"+tag).c_str(),"Convolution", *x, *cb_pdf, *bw0 );

  return bw;
}


RooAbsPdf* shapeSB(string tag, RooRealVar* x, RooRealVar* nSig, RooRealVar* nBkg) {


  RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::DataHandling);

  RooNumConvPdf* bw = shapeZ(tag, x);

  RooRealVar* exp_tau=new RooRealVar( ("expt_"+tag).c_str(), "tau", -0.05, -40., -0.04);
  RooExponential* exp_pdf=new RooExponential( ("exp_pdf_"+tag).c_str(), "bkg shape", *x, *exp_tau );

  //RooRealVar* n_Z=new RooRealVar( ("N_{sig} "+tag).c_str(),"n Z events",501.0, 500., 10000000.);
  // RooRealVar* n_bkg=new RooRealVar( ("N_{bkg} "+tag).c_str(),"n bkg events", 3., 0., 600000.);
  RooArgList* listPdf=new RooArgList( *exp_pdf, *bw );
  RooArgList* listPdfVal=new RooArgList( *nBkg, *nSig );
  RooAddPdf* bw_tot=new RooAddPdf( "bw_EBEB_MC", "PDF ee", *listPdf, *listPdfVal );

  return bw_tot;

}


vector<float> doSingleFit(TH1* histo, bool isData, string category, string channel, string plotDir) {

  RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval); // 1 for INFO
  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::DataHandling);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Fitting);

  RooRealVar mass("mass","m_{ll}",60,120,"GeV");
  RooDataHist data("hist","hist",mass,histo);

  vector<float> v(4,0);

  if(data.sumEntries() < 0.001 || data.numEntries() <10 ) {
    v[0]=histo->Integral(histo->GetXaxis()->FindBin(70), histo->GetXaxis()->FindBin(110) );
    v[1]=sqrt( v[0] ); //not true for MC
    v[2]=histo->Integral()-v[0];
    v[3]=sqrt( v[2] );
    cout<<"result:\t"<<histo->GetName()<<"\t"<<v[0]<<"\t"<<v[1]<<endl;
    return v;
  }

  RooRealVar nSig("nSig","nSig",data.sumEntries(),0,10000000);
  RooRealVar nBkg("nBkg","nBkg",1.,0,10000000);

  // ostringstream os;
  // os<<categ;
  //RooNumConvPdf bw;
  //RooExponential exp_pdf;
  RooAbsPdf* shape=shapeSB( ("s"+category+channel),&mass, &nSig, &nBkg);

  RooFitResult* result;
  result = shape->fitTo(data, RooFit::SumW2Error(kFALSE), RooFit::Save(kTRUE), RooFit::PrintLevel(-1) );
  
  double N=nSig.getVal();
  double eN=nSig.getError();

  double NB=nBkg.getVal();
  double eNB=nBkg.getError();

  cout<<"result:\t"<<histo->GetName()<<"\t"<<N<<"\t"<<eN<< " ... " << result << endl;

  TCanvas* c=new TCanvas( ("c"+category+channel).c_str(),("c"+category+channel).c_str());
  TCanvas* clog=new TCanvas( ("clog"+category+channel).c_str(),("clog"+category+channel).c_str());
  clog->SetLogy();
  c->cd();
  RooPlot* frame=mass.frame();
  data.plotOn(frame);
  shape->plotOn(frame);  
  //std::cout << bw << "   AAAAAAAAAAAAA      " << exp_pdf << std::endl;
  
  string bw_name = "bw_s"+category+channel;
  string exp_name = "exp_pdf_s"+category+channel;
  shape->plotOn(frame,Components(bw_name.data()),LineColor(kGreen)) ;
  shape->plotOn(frame,Components(exp_name.data()),LineColor(kRed)) ;
  frame->Draw();
  
  
  string passtag = "pass";
  if(channel == "OS")
    passtag = "fail";
  string name=Form("%s/%s_fit_s_shapes.png", plotDir.data(), passtag.data());
  c->SaveAs(name.data());

  clog->cd();
  frame->Draw();
  name=Form("%s/%s_log_fit_s_shapes.png", plotDir.data(), passtag.data());
  clog->SaveAs(name.data());

  delete frame;
  delete c;
  delete clog;

  delete shape;
  //delete data;


  v[0]=N;
  v[1]=eN;
  v[2]=NB;
  v[3]=eNB;

  return v;
}


map<string, vector<float> > doFits(string tag, string file, bool isData, string singleCateg) {
  TFile* f=new TFile(file.c_str(), "read");

  vector<float> vs;
  map<string, vector<float> > vals_ss;
  map<string, vector<float> > vals_os;
  
  VString chns = {"SS", "OS"};

  VString cats = {"BB_LL","BB_ML","BB_MM","BB_HL","BB_HM","BB_HH",
      "EE_LL","EE_ML","EE_MM","EE_HL","EE_HM","EE_HH",
      "BE_LL","BE_ML","EB_ML","BE_MM","BE_HL","EB_HL",
      "BE_HM","EB_HM","BE_HH"};

  string datatag = "data";
  if(!isData)
    datatag = "pseudodata";
  string mydir = Form("output_%s_%s", datatag.data(), tag.data());
  FILE *test=fopen( mydir.data(), "r" );
  if( test==0 ) system( Form("mkdir %s", mydir.data()));
  else fclose( test );
  string myFitDir = Form("fit_%s", mydir.data());
  FILE *test2=fopen( myFitDir.data(), "r" );
  if( test2==0 ) system( Form("mkdir %s", myFitDir.data()));
  else fclose( test2 );
  
  
  ofstream myfile;
  myfile.open(Form("%s/results_cat_shapes.txt", mydir.data()));
  for (int i=0; i<1; i++) {
    string cat = cats[i];
	for (auto channel : chns) {
	    std::cout << Form("ttH_charge_flip_%s_%s", cat.data(), channel.data()) << std::endl;
        TDirectory* dir = (TDirectory*)f->Get(Form("ttH_charge_flip_%s_%s", channel.data(), cat.data())) ;
        std::cout << dir << endl;
        TH1* histo = (TH1*)dir->Get("x_data_obs");
        
        string name=Form("%s_%s", channel.data(), cat.data());

        //if(singleCateg!="" && name.find(singleCateg)==string::npos) continue;
        string plotDir = Form("%s/bin%d", myFitDir.data(), i);
        FILE *plottest=fopen( plotDir.data(), "r" );
        if( plottest==0 ) system( Form("mkdir %s", plotDir.data()));
        else fclose( plottest );

        vs = doSingleFit(histo, isData, cat, channel, plotDir);
        if (channel == "SS")
            vals_ss[ cat ] = vs;
        else
            vals_os[ cat ] = vs;
    }
    double ratio = vals_ss[cat][0] / (vals_ss[cat][0] + vals_os[cat][0]);
    std::cout << vals_ss[cat][1] << " " << vals_ss[cat][0] << " " << vals_ss[cat][1] / vals_ss[cat][0];
    double error = getErr(vals_ss[cat][0], vals_os[cat][0], vals_ss[cat][1], vals_os[cat][1]);
    myfile << i << ", " << ratio << ", " << error << ", " << error << "\n";
    
  }
  myfile.close();
  return vals_ss;

}



int main(int argc, char* argv[]) {

  string file;
  bool isData=false;
  string singleCateg="";
  string tag = "noname";
  char c;

  while ((c = getopt(argc, argv, "f:s:D:t:h")) != -1 ) {
    switch (c) {
      //case 'd': { file=optarg; break;}
    case 'f': { file=string(optarg); break;}
    case 's': { singleCateg=string(optarg); break;}
    case 'D': { isData=bool(atoi(optarg)); break; }
    case 't': { tag=string(optarg); break; }
    default : {
      cout<<"configuration options:\n -f : file to read (root) \n -s <categ> perform a fit over a single Z category. \n -D run on data (0 per default). \n -t tag \n -h help \n"<<endl;
      return 0; }
    }
  }

  //==============================================


  map<string, vector<float> > vals=doFits(tag, file, isData, singleCateg);
}

