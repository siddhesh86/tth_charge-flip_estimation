cp $CMSSW_BASE/src/tthAnalysis/ChargeFlipEstimation/macros/MaxLikelihoodFit.cc $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/src/
cd $CMSSW_BASE/src
scram b -j8
cd -
