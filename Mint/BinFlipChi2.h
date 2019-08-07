#ifndef __BINFLIPCHI2_H__
#define __BINFLIPCHI2_H__

#include <Mint/Minimisable.h>
#include <Mint/HadronicParameters.h>
#include <Mint/TimeBinning.h>
#include <string>
#include <Mint/MinuitParameterSet.h>
#include <Mint/counted_ptr.h>
#include <complex>
#include <Mint/FitParameter.h>
#include <utility>

class BinFlipParSet : public MINT::MinuitParameterSet {
 public :
  BinFlipParSet(double, double, double, double,
		double, double, double, double) ;
  BinFlipParSet(const std::string&) ;
  MINT::FitParameter zcp_Re ;
  MINT::FitParameter zcp_Im ;
  MINT::FitParameter deltaz_Re ;
  MINT::FitParameter deltaz_Im ;
  void fixZeroCPV() ;
  std::complex<double> zcp() const ;
  std::complex<double> deltaz() const ;
  static std::pair<std::complex<double>, std::complex<double> > fromXY(double, double, double, double) ;
} ;

class BinFlipChi2 : public MINT::Minimisable {
 public :
  BinFlipChi2(BinFlipParSet*, const HadronicParameters&, const TimeBinning&) ;
  BinFlipChi2(BinFlipParSet*, const std::string&, const std::string&, const std::string& fname = "") ;

  virtual double getVal() override ;
  virtual void parametersChanged() override ;

 private :
  BinFlipParSet* m_fitPars ;
  HadronicParameters m_hadronicPars ;
  TimeBinning m_timeBinning ;
  std::complex<double> m_zcp ;
  std::complex<double> m_dz ;

  void checkPhaseBinning() const ;
  void init() ;
} ;

#endif
