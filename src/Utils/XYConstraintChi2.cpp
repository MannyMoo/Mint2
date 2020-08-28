#include <Mint/XYConstraintChi2.h>
#include <Mint/FitParameter.h>
#include <Mint/BinFlipChi2.h>
#include <complex>
#include <TMath.h>

using MINT::MinuitParameterSet ;
using MINT::FitParameter;
using std::complex;

MinuitParameterSet* XYConstraintChi2::makeParSet(BinFlipParSet* binflipParSet){
  MinuitParameterSet* parset = new MinuitParameterSet();
  parset->add(new FitParameter("x", 0, 0.004, 0.001));
  parset->add(new FitParameter("y", 0, 0.006, 0.001));
  parset->add(new FitParameter("qoverpm1", 0, 0., 0.001));
  parset->add(new FitParameter("phi", 0, 0., 0.001));
  parset->setCovMatrix(binflipParSet->covMatrix());
  return parset;
}

XYConstraintChi2::XYConstraintChi2(BinFlipParSet* binflipParSet) :
  GaussianConstraintChi2(XYConstraintChi2::makeParSet(binflipParSet)),
  m_binflipParSet(binflipParSet)
{}

XYConstraintChi2::~XYConstraintChi2(){
  for(unsigned i = 0 ; i < getParSet()->size(); ++i)
    delete getParSet()->getParPtr(i);
  delete getParSet();
}

TVectorT<double> XYConstraintChi2::getDiffs() {
  auto vals = BinFlipParSet::fromXY(getParSet()->getParPtr(0)->mean(), getParSet()->getParPtr(1)->mean(),
				    getParSet()->getParPtr(2)->mean()+1,
				    getParSet()->getParPtr(3)->mean()*TMath::Pi()/180.);
  complex<double> zdiff = vals.first - m_binflipParSet->zcp();
  complex<double> dzdiff = vals.second - m_binflipParSet->deltaz();
  TVectorT<double> diffs(4);
  diffs[0] = zdiff.real();
  diffs[1] = zdiff.imag();
  diffs[2] = dzdiff.real();
  diffs[3] = dzdiff.imag();
  return diffs;
}

