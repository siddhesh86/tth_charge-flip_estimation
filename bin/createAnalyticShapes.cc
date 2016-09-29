#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooNumConvPdf.h>
#include <RooMsgService.h>
#include <RooWorkspace.h>
#include <RooBreitWigner.h>
#include <RooCBShape.h>
#include <RooExponential.h>

#include "RooCMSShape.h"

using namespace std;
using namespace RooFit;


RooNumConvPdf* shapeZ(string tag, RooRealVar* x) {
    RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);
    RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
    RooMsgService::instance().getStream(1).removeTopic(RooFit::DataHandling);

    RooRealVar* mZ0 = new RooRealVar( ("m_Z0_"+tag).c_str(),"Z0 mass", 91.188, "GeV/c^{2}" );
    RooRealVar* gammaZ0 = new RooRealVar( ("gamma_Z0_"+tag).c_str(), "Z0 width",2.4952, "GeV/c^{2}" );
    RooBreitWigner* bw0 = new RooBreitWigner( ("bw0"+tag).c_str(),"true BW",*x, *mZ0, *gammaZ0);

    RooRealVar* cb_bias  = new RooRealVar( ("cbb_"+tag).c_str(), "bias",  0.07, -3.00, 3.00 );
    RooRealVar* cb_width = new RooRealVar( ("cbw_"+tag).c_str(), "width", 1.00,  0.00, 5.00 );
    RooRealVar* cb_alpha = new RooRealVar( ("cba_"+tag).c_str(), "alpha", 1.20,  0.03, 2.00 );
    RooRealVar* cb_power = new RooRealVar( ("cbn_"+tag).c_str(), "power", 5.00 );

    RooCBShape* cb_pdf = new RooCBShape( ("cb_pdf_"+tag).c_str(), "CB shape",
                                         *x,*cb_bias, *cb_width, *cb_alpha, *cb_power );

    RooNumConvPdf* bw = new RooNumConvPdf( ("bw_"+tag).c_str(),"Convolution", *x, *cb_pdf, *bw0 );
    return bw;
}


RooAddPdf* shapeModel(string tag, RooRealVar* x, RooRealVar* nSig, RooRealVar* nBkg) {
    RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);
    RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
    RooMsgService::instance().getStream(1).removeTopic(RooFit::DataHandling);

    RooNumConvPdf* bw = shapeZ(tag, x);
    
    RooRealVar* exp_tau     = new RooRealVar( ("expt_"+tag).c_str(), "tau", -0.05, -40., -0.04);
    RooExponential* exp_pdf = new RooExponential( ("exp_pdf_"+tag).c_str(), "bkg shape", *x, *exp_tau );

    /*RooRealVar* exp_alpha = new RooRealVar( ("expa_"+tag).c_str(), "alpha", 40.0, 20.0, 160.0);
    RooRealVar* exp_beta  = new RooRealVar( ("expb_"+tag).c_str(), "beta",  0.05, 0.0, 2.0);
    RooRealVar* exp_gamma = new RooRealVar( ("expg_"+tag).c_str(), "gamma", 0.02, 0.0, 0.1);
    RooRealVar* exp_peak  = new RooRealVar( ("expp_"+tag).c_str(), "peak",  91.2);*/
    //RooCMSShape* exp_pdf = new RooCMSShape( ("exp_pdf_"+tag).c_str(), string("bkg shape").c_str(),
    //                                        *x, *exp_alpha, *exp_beta, *exp_gamma, *exp_peak);

    //RooRealVar* nSig   = new RooRealVar( ("N_{sig} "+tag).c_str(),"n Z events",501.0, 500., 10000000.);
    //RooRealVar* nBkg = new RooRealVar( ("N_{bkg} "+tag).c_str(),"n bkg events", 3., 0., 600000.);
    RooArgList* listPdf = new RooArgList( *exp_pdf, *bw );
    RooArgList* listPdfVal = new RooArgList( *nBkg, *nSig );
    RooAddPdf* bw_tot = new RooAddPdf( "bw_EBEB_MC", "PDF_ee", *listPdf, *listPdfVal );

    return bw_tot;
    //return exp_pdf;
}


int main(int argc, char* argv[]) {
    RooRealVar mass("mass","m_{ll}",60,120,"GeV");
    //RooDataHist data("hist","hist",mass,histo);
    //vector<float> values(4,0);

    RooRealVar nSig("nSig","nSig",/*data.sumEntries()*/ 10,0,10000000);
    RooRealVar nBkg("nBkg","nBkg",1.,0,10000000);

    //string os = (string)histo->GetName();
    
    //RooNumConvPdf* sig = shapeZ("sig", &mass);
    //RooAbsPdf* bg = shapeBG("bg", &mass);
    RooAddPdf* model = shapeModel("model", &mass, &nSig, &nBkg);

    RooWorkspace *w = new RooWorkspace("w","workspace") ;
    w->import(*model);
    //w->import(sig);



    w->Print() ;
    w->writeToFile("analytic_workspace.root") ;

    // Workspace will remain in memory after macro finishes
    //gDirectory->Add(w) ;
}
