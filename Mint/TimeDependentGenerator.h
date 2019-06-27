#ifndef __TIMEDEPENDENTGENERATOR_H__
#define __TIMEDEPENDENTGENERATOR_H__

#include "Mint/FitAmpSum.h"
#include "Mint/DalitzEvent.h"
#include "Mint/SignalGenerator.h"
#include <memory>
#include <list>
#include <map>
#include <complex>
#include <string>
#include <Mint/SplineGenerator.h>

class TRandom3 ;

// Class to generate time dependent phase space events.
class TimeDependentGenerator {

public :
  // Class to hold the generator at a given time point.
  class GenTimePoint {
  public :
    GenTimePoint(const double _decaytime, FitAmpSum* _model,
		 const double _integral, SignalGenerator* _generator) ;

    const double decaytime ;
    const double integral ;
    std::unique_ptr<FitAmpSum> model ;
    std::unique_ptr<SignalGenerator> generator ;
  } ;

  typedef std::list<GenTimePoint> GenList ;
  typedef std::map<int, GenList> GenMap ;
  typedef std::pair<std::complex<double>, std::complex<double> > AmpPair ;

  // Class to hold the tag (production flavour), decay time and
  // Dalitz event that's generated.
  struct GenTimeEvent {
    int tag ;
    double decaytime ;
    MINT::counted_ptr<IDalitzEvent> evt ;
  } ;

  // Take the CP conjugate of the head of the decay pattern.
  static DalitzEventPattern anti(DalitzEventPattern pat) ;

  /* Constructor, takes:
   name : the name of the generator and the directory in which the integrators will be saved.
   overwrite : whether to overwrite the existing integrator files (if they exist).
   rndm : The random number generator to use.
   precision : The precision to which the integrals must be calculated.
   pattern : The event pattern to be used (the CP conjugate will automatically be added).
   width : the decay width in 1/ps.
   deltam : the delta-mass in 1/ps.
   deltagamma : the delta-gamma in 1/ps.
   qoverp : the magnitude of q/p.
   phi : the phase of q/p.
   tmax : the maximum decay time that'll be generated.
   ntimepoints : the number of points to sample between 0 and tmax when building the generators.
   h_efficiency : (optional) histogram to which efficiency plot will be fitted
  */
  TimeDependentGenerator(const std::string& name, const bool overwrite, TRandom3* rndm, double precision,
			 const DalitzEventPattern& pattern, double width, double deltam,
			 double deltagamma,
			 double qoverp, double phi, double tmax, int ntimepoints,
			 const bool saveIntegEvents = true, double tmin = 0., TH1F* h_efficiency = NULL) ;

  // Get the coefficients of the amplitudes for the produced flavour and the mixed flavour
  // given the tag and decay time.
  AmpPair amplitude_coefficients(const int tag, const double decaytime) ;

  // Generate a flavour.
  int generate_tag() const ;

  // Generate a decay time for the given flavour.
  double generate_decay_time(const int tag) const ;

  // Generate a Dalitz event for the given flavour and decay time.
  MINT::counted_ptr<IDalitzEvent> generate_dalitz_event(const int tag, const double decaytime) const ;

  // Generate a flavour, decay time and Dalitz event.
  GenTimeEvent generate_event() const ;

  // Get the decay time generators.
  const std::map<int, SplineGenerator> time_generators() const ;
  
private :
  const std::string m_name ;
  TRandom3* m_rndm ;
  const DalitzEventPattern m_pattern ;
  const DalitzEventPattern m_cppattern ;
  
  const double m_width ;
  const double m_deltam ;
  const double m_deltagamma ;
  const double m_qoverp ;
  const double m_phi ;

  const double m_tmax ;
  const double m_tmin ;
  const int m_ntimepoints ;

  GenMap m_genmap ;

  std::map<int, SplineGenerator> m_timegenerators ;

  double m_tagintegralfrac ;
  double m_precision ;

  TH1F* m_h_efficiency;
} ;

#endif
