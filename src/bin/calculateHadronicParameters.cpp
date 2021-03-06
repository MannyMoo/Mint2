#include <Mint/HadronicParameters.h>
#include <Mint/NamedParameter.h>
#include <Mint/NamedParameterBase.h>
#include <Mint/DalitzEvent.h>
#include <Mint/DalitzEventPattern.h>
#include <TRandom3.h>
#include <Mint/FitAmpSum.h>
#include <iostream>
#include <fstream>
#include <ctime>

using MINT::NamedParameter ;
using namespace std ;

int calculateHadronicParameters(const string& config = string()) {
  if(config.size() > 0){
    cout << "Using config file " << config << endl ;
    MINT::NamedParameterBase::setDefaultInputFile(config) ;
  }

  NamedParameter<int> EventPattern("Event Pattern", 421, -321, 211, 211, -211);
  DalitzEventPattern pat(EventPattern.getVector());
  DalitzEventPattern cpPat(pat) ;
  cpPat[0].antiThis() ;
  
  HadronicParameters::ModelPtr model(new FitAmpSum(pat)) ;
  HadronicParameters::ModelPtr cpmodel(new FitAmpSum(cpPat)) ;
  
  TRandom3 ranLux;
  NamedParameter<int> RandomSeed("RandomSeed", 0);
  ranLux.SetSeed((int)RandomSeed);
  gRandom = &ranLux;

  NamedParameter<int> nBinsPhase("nBinsPhase", 8) ;

  NamedParameter<string> binningType("binning_type", string("model"), (char*)0) ;
  HadronicParameters::BinningPtr binning(0) ;
  if(string(binningType) == string("model3body"))
    binning = HadronicParameters::BinningPtr(new HadronicParameters::ModelBinning3Body(model, cpmodel, nBinsPhase)) ;
  else
    binning = HadronicParameters::BinningPtr(new HadronicParameters::ModelPhaseBinning(model, cpmodel, nBinsPhase)) ;
  
  NamedParameter<string> parsName("parsName", string("hadronicPars"), (char*)0) ;
  HadronicParameters pars(parsName, binning) ;

  NamedParameter<int>  Nevents("Nevents", 10000);
  cout << "Calculating parameters with random seed " << int(RandomSeed) << " for " 
       << int(Nevents) << " events." << endl ;
  int startTime(time(0)) ;
  pars.add(pat, ranLux, int(Nevents)) ;
  int endTime(time(0)) ;
  cout << "Finished calculating. Took " << (endTime - startTime) << " s, " 
       << float(endTime - startTime)/Nevents << " s/event." << endl ;

  cout << "Before finalisation:" << endl;
  pars.Print();
  
  pars.finaliseSum() ;

  cout << "Calculated parameters:" << endl ;
  pars.Print() ;

  NamedParameter<string> outputFile("outputFile", string("hadronicParameters.txt"), (char*)0) ;
  pars.write(outputFile) ;

  return 0 ;
}

int main(int argc, char** argv) {
  if(argc > 1)
    return calculateHadronicParameters(string(argv[1])) ;
  return calculateHadronicParameters() ;
}
