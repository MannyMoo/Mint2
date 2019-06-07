#include <Mint/PhaseDifferenceCalc.h>

PhaseDifferenceCalc::PhaseDifferenceCalc(const DalitzEventPattern& pat, const char* fname) :
  m_model(pat, fname),
  m_cpmodel(TimeDependentGenerator::anti(pat), fname)
{}

// Can't be made const as IReturnComplexForEvent::ComplexVal is not const.
// Returns the phase of A/Abar for the given point in phase space.
double PhaseDifferenceCalc::phase_difference(IDalitzEvent& evt) {
  std::complex<double> val = m_model.ComplexVal(evt) ;
  std::complex<double> cpval = m_cpmodel.ComplexVal(evt) ;
  val /= cpval ;
  return std::arg(val) ;
}

// Returns the cross term A* x Abar fro the given point in phase space.
std::complex<double> PhaseDifferenceCalc::cross_term(IDalitzEvent& evt) {
  std::complex<double> val = m_model.ComplexVal(evt) ;
  std::complex<double> cpval = m_cpmodel.ComplexVal(evt) ;
  return std::conj(val) * cpval ;
}
