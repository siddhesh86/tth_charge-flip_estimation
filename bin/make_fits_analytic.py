vector<float> doSingleFit(TH1* histo, bool isData) {

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
  string os=(string)histo->GetName();
  RooAbsPdf* shape=shapeSB( ("s"+os),&mass, &nSig, &nBkg);

  RooFitResult* result;
  result = shape->fitTo( data ,RooFit::SumW2Error(kFALSE), RooFit::Save(kTRUE), RooFit::PrintLevel(-1) );
  
  double N=nSig.getVal();
  double eN=nSig.getError();

  double NB=nBkg.getVal();
  double eNB=nBkg.getError();

  cout<<"result:\t"<<histo->GetName()<<"\t"<<N<<"\t"<<eN<< " ... " << result << endl;

  TCanvas* c=new TCanvas( ("c"+os).c_str(),("c"+os).c_str());
  RooPlot* frame=mass.frame();
  data.plotOn(frame);
  shape->plotOn(frame);
  frame->Draw();

  FILE *test=fopen( "plots", "r" );
  if( test==0 ) system( "mkdir plots");
  else fclose( test );

  string name="plots/fitData_";
  if(!isData) name="plots/fitMC_";
  c->SaveAs( (name+os+".png").c_str() );

  delete frame;
  delete c;

  delete shape;
  //delete data;


  v[0]=N;
  v[1]=eN;
  v[2]=NB;
  v[3]=eNB;

  return v;
}
