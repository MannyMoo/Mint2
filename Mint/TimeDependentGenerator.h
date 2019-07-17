#ifndef __TIMEDEPENDENTGENERATOR_H__
#define __TIMEDEPENDENTGENERATOR_H__

#include <Mint/FitAmpSum.h>
#include <Mint/DalitzEvent.h>
#include <Mint/SignalGenerator.h>
#include <Mint/FitAmpIncoherentSum.h>
#include <memory>
#include <list>
#include <map>
#include <complex>
#include <string>
#include <TSpline.h>

class TRandom3 ;

// Class to generate time dependent phase space events.
class TimeDependentGenerator {

public :
  // Class to hold the tag (production flavour), decay time and
  // Dalitz event that's generated.
  class GenTimeEvent : public DalitzEvent {
  public :
    GenTimeEvent(const IDalitzEvent&, const int, const double,
		 const double smeareddecaytime = -999.) ;
    GenTimeEvent(const IDalitzEvent&) ;
    int getTag() const ;
    void setTag(int) ;
    double getDecayTime() const ;
    void setDecayTime(double) ;
    double getSmearedDecayTime() const ;
    void setSmearedDecayTime(double) ;
    enum InfoIndex{ITAG = 0, IDECAYTIME, ISMEAREDDECAYTIME} ;
    static std::map<std::string, unsigned> infoNames ;
  private :
    static std::map<std::string, unsigned> makeInfoNames() ;
  } ;

  // Take the CP conjugate of the head of the decay pattern.
  static DalitzEventPattern anti(DalitzEventPattern pat) ;

  /* Constructor, takes:
   pattern : The event pattern to be used (the CP conjugate will automatically be added).
   width : the decay width in 1/ps.
   deltam : the delta-mass in 1/ps.
   deltagamma : the delta-gamma in 1/ps.
   qoverp : the magnitude of q/p.
   phi : the phase of q/p.
   rndm : The random number generator to use.
   h_efficiency : (optional) histogram to which efficiency plot will be fitted
  */
  TimeDependentGenerator(const DalitzEventPattern& pattern, double width, double deltam,
			 double deltagamma, double qoverp, double phi,
			 TRandom3* rndm, TH1F* h_efficiency = NULL,
                         float resWidth = 0.05, bool addExpEffects = false) ;
  TimeDependentGenerator(MINT::counted_ptr<FitAmpSum> model, MINT::counted_ptr<FitAmpSum> cpmodel,
			 double width, double deltam,
			 double deltagamma, double qoverp, double phi,
			 TRandom3*, TH1F* h_efficiency = NULL,
                         float resWidth = 0.05, bool addExpEffects = false) ;

  // Generate a decay time for the given flavour.
  std::pair<double, double> generate_decay_time() const ;

  // Generate a flavour, decay time and Dalitz event.
  MINT::counted_ptr<IDalitzEvent> generate_event() ;

  double get_scale() const ;
  float get_gen_efficiency() const ;

  double pdf_value(int, double, IDalitzEvent&) ;
  double pdf_value(IDalitzEvent&) ;

  typedef std::pair<std::complex<double>, std::complex<double> > AmpPair ;
  /** Get the coefficients of the amplitudes for the produced flavour and the mixed flavour
      given the tag and decay time. */
  AmpPair amplitude_coefficients(const int tag, const double decaytime) ;

private :
  TRandom3* m_rndm ;
  const DalitzEventPattern m_pattern ;
  const DalitzEventPattern m_cppattern ;
  MINT::counted_ptr<FitAmpSum> m_model ;
  MINT::counted_ptr<FitAmpSum> m_cpmodel ;
  FitAmpIncoherentSum m_bothmodel ;
  SignalGenerator m_generator ;

  const double m_width ;
  const double m_deltam ;
  const double m_deltagamma ;
  const std::complex<double> m_qoverp ;
  double m_scale ;
  unsigned m_ngen ;
  unsigned m_naccept ;

  TH1F* m_h_efficiency;
  TSpline3 m_efficiencyFit;

  float m_resWidth;
  bool m_addExpEffects;
} ;

#endif
