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
#include <TH2F.h>

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

class BinFlipChi2Base : public MINT::Minimisable {
 public :
  BinFlipChi2Base(BinFlipParSet*) ;
  virtual void parametersChanged() override ;
  TH2F scan2D(unsigned, unsigned, float nSigmaRange = 5., unsigned nBins = 300,
	      bool zeroCentre = true, float scale = 1.) ;
  BinFlipParSet* getBinFlipParSet() ;
 protected :
  BinFlipParSet* m_fitPars ;
  std::complex<double> m_zcp ;
  std::complex<double> m_dz ;
} ;  

class BinFlipChi2 : public BinFlipChi2Base {
 public :
  BinFlipChi2(BinFlipParSet*, const HadronicParameters&, const TimeBinning&) ;
  BinFlipChi2(BinFlipParSet*, const std::string&, const std::string&, const std::string& fname = "") ;

  virtual double getVal() override ;

 private :
  HadronicParameters m_hadronicPars ;
  TimeBinning m_timeBinning ;

  void checkPhaseBinning() const ;
  void init() ;
} ;

class BinFlipChi2Simul : public BinFlipChi2Base {
 public :
  typedef std::deque<BinFlipChi2*> BinFlippers ;
  BinFlipChi2Simul(BinFlipChi2*, BinFlipChi2*) ;
  BinFlipChi2Simul(const BinFlippers&) ;

  virtual double getVal() override ;
  virtual void parametersChanged() override ;
 private :
  BinFlippers m_flippers ;
} ;

#endif
