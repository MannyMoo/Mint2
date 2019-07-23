// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:18:01 GMT
#include "Mint/NamedParameter.h"
#include "Mint/DalitzEventList.h"
#include "Mint/DalitzEvent.h"
#include "TFile.h"
#include "TRandom3.h"
#include <ctime>
#include <iostream>
#include <memory>

#include <Mint/SplineGenerator.h>
#include <Mint/TimeDependentGenerator.h>
#include <Mint/Eff3piSymmetric.h>

using namespace std;
using namespace MINT;


// ===========
// Main method
// ===========
int genTimeDependent(){
  time_t startTime = time(0);

  TRandom3 ranLux;
  NamedParameter<int> RandomSeed("RandomSeed", 0);
  ranLux.SetSeed((int)RandomSeed);
  gRandom = &ranLux;

  FitAmplitude::AutogenerateFitFile();

  NamedParameter<int>  Nevents("Nevents", 10000);

  NamedParameter<int> EventPattern("Event Pattern", 421, -321, 211, 211, -211);
  DalitzEventPattern pat(EventPattern.getVector());
  DalitzEventPattern cpPat(pat) ;
  cpPat[0].antiThis() ;

  SignalGenerator genD0(pat);
  SignalGenerator genD0bar(cpPat);

  NamedParameter<int> saveEvents("SaveEvents", 1);
  NamedParameter<int> genTimeDependent("genTimeDependent", 0);
  NamedParameter<double> lifetime("lifetime", 0.4101) ;
  double width = 1./lifetime ;
  NamedParameter<double> x("x", 0.0039) ;
  double deltam = x * width ;
  NamedParameter<double> y("y", 0.0065) ;
  double deltagamma = y * width * 2. ;
  NamedParameter<double> qoverp("qoverp", 1.) ;
  NamedParameter<double> phi("phi", 0.) ;

  NamedParameter<string> efficiencyFile("efficiencyFile", string("/home/ppe/n/nmchugh/SummerProject/DaVinciDev_v44r10p1/AGammaD0Tohhpi0/scripts/mint/h_efficiency.root")) ;
  NamedParameter<string> h_efficiencyName( "h_efficiencyName", string("h_efficiency") ) ;
  
  NamedParameter<int> addExpEffects("addExpEffects", 0) ;
  NamedParameter <float> resWidth("resWidth", 0.05) ;

  TH1F* h_efficiency = NULL ; 
  Eff3piSymmetric* sEfficiency = NULL;
  if((bool)addExpEffects){
    TFile* eff_infile = TFile::Open( ((string) efficiencyFile).c_str() ) ;
    eff_infile->GetObject(((string) h_efficiencyName).c_str(), h_efficiency) ;
    sEfficiency = new Eff3piSymmetric();
  }

  
  cout << " got event pattern: " << pat << endl;

  unique_ptr<TimeDependentGenerator> timedepgen ;
  if(genTimeDependent){
    int startinit(time(0)) ;
    timedepgen.reset(new TimeDependentGenerator(pat,
						width, deltam, deltagamma, qoverp, phi, &ranLux, 
						h_efficiency, resWidth, (bool)addExpEffects, sEfficiency)) ;
    cout << "Initialise TimeDependentGenerator took " << time(0) - startinit << " s" << endl ;
  }

  if(!(bool)saveEvents)
    return 0 ;

  DalitzEventList eventList1 ;

  cout << "Generating " << Nevents << " signal events (1)." << endl;
  int startTimeGen(time(0)) ;
  for(int i = 0 ; i < Nevents ; ++i){

    if(i%1000 == 0)
      cout << "Generating candidate " << i << " (" << (time(0)-startTimeGen)/float(i)
	   << " s per candidate)" << endl ;

    MINT::counted_ptr<IDalitzEvent> evt(0) ;
    if(genTimeDependent){
      evt = timedepgen->generate_event() ;
    }
    else {
      // Decide if it's a D0 or D0bar that's being generated.
      SignalGenerator* gen(&genD0) ;
      int tag = 1 ;
      if(ranLux.Rndm() > 0.5){
	gen = &genD0bar ;
	tag = -1 ;
      }

      // Generate the decay time of the candidate.
      double decaytime = ranLux.Exp(lifetime) ;

      evt = gen->newEvent() ;
      evt->setValueInVector(TimeDependentGenerator::GenTimeEvent::ITAG, tag) ;
      evt->setValueInVector(TimeDependentGenerator::GenTimeEvent::IDECAYTIME, decaytime) ;
      evt->setValueInVector(TimeDependentGenerator::GenTimeEvent::ISMEAREDDECAYTIME, -999.) ;
    }
    eventList1.Add(*evt) ;
  }
  
  if(genTimeDependent){
    cout << "Generator efficiency: " << timedepgen->get_gen_efficiency() << endl ;
    cout << "Generator scale: " << timedepgen->get_scale() << endl ;
  }
  
  // Make sure the first event in the list is a D0 so the naming scheme is consistent.
  if(eventList1.begin()->eventPattern() != pat){
    auto ievt = eventList1.begin() ;
    ++ievt ;
    for( ; ievt->eventPattern() != pat && ievt != eventList1.end() ; ++ievt)
      continue ;
    if(ievt != eventList1.end()){
      DalitzEvent evt = *ievt ;
      eventList1.erase(ievt) ;
      eventList1.theVector().insert(eventList1.begin(), evt) ;
    }
  }

  NamedParameter<string> outputfname("outputFileName", string("pipipi0_1.root"), (char*)0) ;
  eventList1.save(outputfname);
  string fnamestr(outputfname) ;
  TFile tuplefile(fnamestr.c_str(), "update") ;
  TNtupleD* ntuple = (TNtupleD*)tuplefile.Get("DalitzEventList") ;
  int tag(0) ;
  double decaytime(0.) ;
  double smeareddecaytime(0.) ;
  TBranch* tagbranch = ntuple->Branch("tag", &tag, "tag/I") ;
  TBranch* decaytimebranch = ntuple->Branch("decaytime", &decaytime, "decaytime/D") ;
  TBranch* smeareddecaytimebranch = ntuple->Branch("smeareddecaytime", &smeareddecaytime, "smeareddecaytime/D") ;
  for(const auto& ievt : eventList1){
    tag = ievt.getValueFromVector(TimeDependentGenerator::GenTimeEvent::infoNames["tag"]) ;
    decaytime = ievt.getValueFromVector(TimeDependentGenerator::GenTimeEvent::infoNames["decaytime"]) ;
    smeareddecaytime = ievt.getValueFromVector(TimeDependentGenerator::GenTimeEvent
					       ::infoNames["smeareddecaytime"]) ;
    tagbranch->Fill() ;
    decaytimebranch->Fill() ;
    smeareddecaytimebranch->Fill() ;
  }
  ntuple->Write(ntuple->GetName(), TObject::kWriteDelete) ;
  tuplefile.Close() ;
  
  DalitzHistoSet datH = eventList1.histoSet();
  datH.save("plotsFromEventList.root");

  cout << " genTimeDependent done. Took " << (time(0) - startTime)/60. 
       << " min. Returning 0." << endl;

  return 0;
}


int main(){

  return genTimeDependent();

}
//
