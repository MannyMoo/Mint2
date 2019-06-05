#ifndef __PHASEDIFFERENCECALC_H__
#define __PHASEDIFFERENCECALC_H__

#include "Mint/TimeDependentGenerator.h"
#include <memory>
#include "Mint/DalitzEvent.h"
#include "Mint/FitAmpSum.h"

class PhaseDifferenceCalc {
 public :
  PhaseDifferenceCalc(const DalitzEventPattern&, const char* fname=0) ;
  
  // Can't be made const as IReturnComplexForEvent::ComplexVal is not const.
  double phase_difference(IDalitzEvent&) ;
 private :
  FitAmpSum m_model ;
  FitAmpSum m_cpmodel ;

} ;

#endif
