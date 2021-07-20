#include <Mint/GaussianConstraintChi2.h>
#include <Mint/FitParameter.h>
#include <stdexcept>

using MINT::MinuitParameterSet ;
using MINT::FitParameter;
using std::vector;

GaussianConstraintChi2::GaussianConstraintChi2(MinuitParameterSet* parset,
					       const ParVect& pars,
					       const CovMatrix& covMatrix) :
  Minimisable(parset),
  m_pars(pars),
  m_covMatrixInv(covMatrix)
{
  // No pars given, take all floating pars from the parset
  if(m_pars.size() == 0)
    for(unsigned i = 0 ; i < parset->size() ; ++i){
      FitParameter* par = (FitParameter*)parset->getParPtr(i);
      if(par->iFixInit() == FitParameter::FIT)
	m_pars.push_back(par);
    }

  // No cov matrix given
  if(m_covMatrixInv.GetNcols() == 0){
    // If the parset has a defined cov matrix, use that
    if(parset->covMatrix().GetNcols() > 0){
      m_covMatrixInv = parset->covMatrix() ;
    }
    // Otherwise, just a diagonal cov matrix of the errors^2
    else {
      m_covMatrixInv = CovMatrix(m_pars.size());
      for(unsigned i = 0 ; i < m_pars.size() ; ++i)
	m_covMatrixInv[i][i] = pow(m_pars[i]->stepInit(), 2) ;
    }
  }

  if(m_pars.size() != m_covMatrixInv.GetNcols())
    throw std::length_error("GaussianConstraintChi2: Incompatible number of\
 parameters and size of covariance matrix!");
  
  m_covMatrixInv.Invert() ;
}

TVectorT<double> GaussianConstraintChi2::getDiffs() {
  TVectorT<double> diffs(getParSet()->size()) ;
  for(unsigned i = 0 ; i < getParSet()->size() ; ++i){
    auto* par = getParSet()->getParPtr(i) ;
    diffs[i] = par->mean() - par->meanInit() ;
  }
  return diffs;
}

double GaussianConstraintChi2::getVal() {
  TVectorT<double> diffs = getDiffs();
  return m_covMatrixInv.Similarity(diffs) ;
}
