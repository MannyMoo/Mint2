#include <Mint/GaussianConstraintChi2.h>
#include <TVectorT.h>

using MINT::MinuitParameterSet ;
GaussianConstraintChi2::GaussianConstraintChi2(MinuitParameterSet* parset) :
  Minimisable(parset),
  m_covMatrixInv(parset->size())
{
  if(parset->covMatrix().GetNcols() > 0){
    m_covMatrixInv = parset->covMatrix() ;
  }
  else {
    for(unsigned i = 0 ; i < parset->size() ; ++i)
      m_covMatrixInv[i][i] = pow(parset->getParPtr(i)->stepInit(), 2) ;
  }
  m_covMatrixInv.Invert() ;
}

double GaussianConstraintChi2::getVal() {
  TVectorT<double> diffs(getParSet()->size()) ;
  for(unsigned i = 0 ; i < getParSet()->size() ; ++i){
    auto* par = getParSet()->getParPtr(i) ;
    diffs[i] = par->mean() - par->meanInit() ;
  }
  return m_covMatrixInv.Similarity(diffs) ;
}
