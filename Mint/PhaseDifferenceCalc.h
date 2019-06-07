#ifndef __PHASEDIFFERENCECALC_H__
#define __PHASEDIFFERENCECALC_H__

#include "Mint/TimeDependentGenerator.h"
#include <memory>
#include "Mint/DalitzEvent.h"
#include "Mint/FitAmpSum.h"
#include <complex>

class PhaseDifferenceCalc {
 public :
  PhaseDifferenceCalc(const DalitzEventPattern&, const char* fname=0) ;
  
  // Can't be made const as IReturnComplexForEvent::ComplexVal is not const.
  // Returns the phase of A/Abar for the given point in phase space.
  double phase_difference(IDalitzEvent&) ;

  // Returns the cross term A* x Abar fro the given point in phase space.
  std::complex<double> cross_term(IDalitzEvent&) ;
 private :
  FitAmpSum m_model ;
  FitAmpSum m_cpmodel ;

} ;

#endif
