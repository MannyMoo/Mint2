// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:18:01 GMT
#include "Mint/NamedParameter.h"
#include "Mint/DiskResidentEventList.h"
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

  NamedParameter<double> s13("s13", 0., NamedParameterBase::QUIET) ;
  NamedParameter<double> s23("s23", 0., NamedParameterBase::QUIET) ;

  NamedParameter<string> efficiencyFile("efficiencyFile", string("/home/ppe/n/nmchugh/SummerProject/DaVinciDev_v44r10p1/AGammaD0Tohhpi0/scripts/mint/h_efficiency.root")) ;
  NamedParameter<string> h_efficiencyName( "h_efficiencyName", string("h_efficiency") ) ;
  
  NamedParameter<int> addExpEffects("addExpEffects", 0) ;
  NamedParameter <float> resWidth("resWidth", 0.05) ;

  NamedParameter<string> outputfname("outputFileName", string("pipipi0_1.root"), (char*)0) ;

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

  DiskResidentEventList eventList1(pat, outputfname) ;
  TNtupleD* evttuple = (TNtupleD*)gDirectory->Get("DalitzEventList") ;
  TTree* ntuple = new TTree("TimeEventList", "TimeEventList") ;
  evttuple->AddFriend(ntuple) ;
  int tag(0) ;
  double decaytime(0.) ;
  double smeareddecaytime(0.) ;
  TBranch* tagbranch = ntuple->Branch("tag", &tag, "tag/I") ;
  TBranch* decaytimebranch = ntuple->Branch("decaytime", &decaytime, "decaytime/D") ;
  TBranch* smeareddecaytimebranch = ntuple->Branch("smeareddecaytime", &smeareddecaytime, "smeareddecaytime/D") ;

  cout << "Generating " << Nevents << " signal events (1)." << endl;
  int startTimeGen(time(0)) ;
  for(int i = 0 ; i < Nevents ; ++i){

    if(i%1000 == 0)
      cout << "Generating candidate " << i << " (" << (time(0)-startTimeGen)/float(i)
	   << " s per candidate)" << endl ;

    MINT::counted_ptr<IDalitzEvent> evt(0) ;
    if(genTimeDependent){
      if(s13 && s23)
	evt = timedepgen->generate_event(s13, s23) ;
      else 
	evt = timedepgen->generate_event() ;
    }
    else {
      // Decide if it's a D0 or D0bar that's being generated.
      SignalGenerator* gen(&genD0) ;
      int _tag = 1 ;
      if(ranLux.Rndm() > 0.5){
	gen = &genD0bar ;
	_tag = -1 ;
      }

      // Generate the decay time of the candidate.
      double _decaytime = ranLux.Exp(lifetime) ;

      evt = gen->newEvent() ;
      evt->setValueInVector(TimeDependentGenerator::GenTimeEvent::ITAG, _tag) ;
      evt->setValueInVector(TimeDependentGenerator::GenTimeEvent::IDECAYTIME, _decaytime) ;
      evt->setValueInVector(TimeDependentGenerator::GenTimeEvent::ISMEAREDDECAYTIME, _decaytime) ;
    }
    eventList1.Add(*evt) ;
    tag = evt->getValueFromVector(TimeDependentGenerator::GenTimeEvent::ITAG) ;
    decaytime = evt->getValueFromVector(TimeDependentGenerator::GenTimeEvent::IDECAYTIME) ;
    smeareddecaytime = evt->getValueFromVector(TimeDependentGenerator::GenTimeEvent::ISMEAREDDECAYTIME) ;
    ntuple->Fill() ;
  }
  
  if(genTimeDependent){
    cout << "Generator efficiency: " << timedepgen->get_gen_efficiency() << endl ;
    cout << "Generator scale: " << timedepgen->get_scale() << endl ;
  }
  
  cout << "Save data" << endl ;
  ntuple->Write() ;
  eventList1.Close() ;
  
  cout << "Save Dalitz plots." << endl ;
  DalitzHistoSet datH = eventList1.histoSet();
  datH.save("plotsFromEventList.root");

  cout << " genTimeDependent done. Took " << (time(0) - startTime)/60. 
       << " min. Returning 0." << endl;

  // For debugging, plot PDF vs time at given point in phase space for each tag.
  if(s13 && s23){
    TFile fplots("plotsVsTime.root", "recreate") ;
    unsigned nbins = 1000 ;
    float tmin = 0. ;
    float tmax = 8.2 ;
    DalitzEvent evt(cpPat, s13, s23) ;
    for(int tag = -1 ; tag < 2 ; tag += 2){
      evt.setValueInVector(TimeDependentGenerator::GenTimeEvent::ITAG, tag) ;
      std::ostringstream tagstr ; 
      tagstr << tag ;
      TH1F pdf = timedepgen->draw_pdf_vs_time(evt, nbins, tmin, tmax,
					      "pdf_vs_time_tag" + tagstr.str()) ;
      TH1F env = timedepgen->draw_envelope_vs_time(evt, nbins, tmin, tmax,
						   "env_vs_time_tag" + tagstr.str()) ;
      pdf.Write() ;
      env.Write() ;
      evt = DalitzEvent(pat, s13, s23) ;
    }
    fplots.Close() ;
  }

  return 0;
}


int main(){

  return genTimeDependent();

}
//
