#include <Mint/HadronicParameters.h>
#include <Mint/NamedParameter.h>
#include <Mint/NamedParameterBase.h>
#include <Mint/DalitzEvent.h>
#include <Mint/DalitzEventPattern.h>
#include <TRandom3.h>
#include <Mint/FitAmpSum.h>
#include <iostream>
#include <fstream>

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

  MINT::counted_ptr<HadronicParameters::IBinSign> binSignPlus(new HadronicParameters::BinSign(model)) ;
  MINT::counted_ptr<HadronicParameters::IBinSign> binSignMinus(new HadronicParameters::BinSign(cpmodel)) ;
  
  HadronicParameters parsPlus(model, cpmodel, nBinsPhase, binSignPlus) ;
  HadronicParameters parsMinus(cpmodel, model, nBinsPhase, binSignMinus) ;

  NamedParameter<int>  Nevents("Nevents", 10000);
  for(unsigned i = 0 ; i < int(Nevents) ; ++i){
    DalitzEvent evt(pat, &ranLux) ;
    parsPlus.add(evt) ;
    parsMinus.add(evt) ;
  }

  parsPlus.normalise() ;
  parsMinus.normalise() ;
  cout << "Plus parameters (" << pat << ") :" << endl ;
  parsPlus.Print("plusPars") ;
  cout << endl ;
  cout << "Minus parameters (" << cpPat << ") :" << endl ;
  parsMinus.Print("minusPars") ;

  NamedParameter<string> outputFile("outputFile", string("hadronicParameters.txt"), (char*)0) ;
  ofstream fout ;
  fout.open(outputFile) ;
  parsPlus.Print("parsPlus", fout) ;
  parsMinus.Print("parsMinus", fout) ;
  fout.close() ;

  return 0 ;
}

int main(int argc, char** argv) {
  if(argc > 1)
    return calculateHadronicParameters(string(argv[1])) ;
  return calculateHadronicParameters() ;
}
