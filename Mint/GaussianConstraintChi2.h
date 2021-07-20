#ifndef __GAUSSIANCONSTRAINTCHI_H__
#define __GAUSSIANCONSTRAINTCHI_H__

#include <Mint/Minimisable.h>
#include <TVectorT.h>
#include <vector>

namespace MINT{
  class FitParameter;
}

class GaussianConstraintChi2 : public MINT::Minimisable {
 private :
  typedef std::vector<MINT::FitParameter*> ParVect;
  typedef MINT::MinuitParameterSet::CovMatrix CovMatrix;
  ParVect m_pars;
  CovMatrix m_covMatrixInv ;
  virtual TVectorT<double> getDiffs() ;
 public :
  GaussianConstraintChi2(MINT::MinuitParameterSet*,
			 const ParVect& pars = ParVect(),
			 const CovMatrix& covMatrix = CovMatrix()) ;
  virtual ~GaussianConstraintChi2() {} ;
  virtual double getVal() override ;
} ;

#endif
