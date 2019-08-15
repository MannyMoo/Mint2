#ifndef __GAUSSIANCONSTRAINTCHI_H__
#define __GAUSSIANCONSTRAINTCHI_H__

#include <Mint/Minimisable.h>

class GaussianConstraintChi2 : public MINT::Minimisable {
 private :
  MINT::MinuitParameterSet::CovMatrix m_covMatrixInv ;
 public :
  GaussianConstraintChi2(MINT::MinuitParameterSet*) ;
  virtual ~GaussianConstraintChi2() {} ;
  virtual double getVal() override ;
} ;

#endif
