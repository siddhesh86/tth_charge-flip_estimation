cp $CMSSW_BASE/src/tthAnalysis/ChargeFlipEstimation/macros/PostFitShapesFromWorkspace.cpp $CMSSW_BASE/src/CombineHarvester/CombineTools/bin/PostFitShapesFromWorkspace.cpp
cd $CMSSW_BASE/src
scram b -j8
cd -
