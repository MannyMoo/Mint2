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
#include <TH2D.h>

class BinFlipChi2Base;

class BinFlipParSet : public MINT::MinuitParameterSet {
  friend class BinFlipChi2Base;
 public :
  BinFlipParSet(double, double, double, double,
		double, double, double, double,
		unsigned long blindingSeed = 0, double zBlindRange = 0., double deltazBlindRange = 0.) ;
  BinFlipParSet(const std::string&) ;
  MINT::FitParameter zcp_Re ;
  MINT::FitParameter zcp_Im ;
  MINT::FitParameter deltaz_Re ;
  MINT::FitParameter deltaz_Im ;

  void fixZeroCPV() ;
  void floatOnlyDeltaY();

  std::complex<double> zcp() const ;
  std::complex<double> deltaz() const ;
  static std::pair<std::complex<double>, std::complex<double> > fromXY(double, double, double, double) ;
 private:
  std::complex<double> unblindZcp() const ;
  std::complex<double> unblindDeltaz() const ;
} ;

class BinFlipChi2Base : public MINT::Minimisable {
 public :
  BinFlipChi2Base(BinFlipParSet*) ;
  virtual void parametersChanged() override ;
  /// Perform a 2D scan of the chi2, while keeping the other two variables fixed.
  TH2D scan2D(unsigned, unsigned, float nSigmaRange = 5., unsigned nBins = 300,
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
  std::deque<std::deque<TH1F> > plotsVsTime(const std::string&) const ;
  void savePlotsVsTime(const std::string&, TFile&) const ;

 protected:
  /// Get the chi2 for the given time & phase bin.
  virtual double binChiSquared(unsigned, unsigned) const;
  HadronicParameters m_hadronicPars ;
  TimeBinning m_timeBinning ;
  double m_lifetime ;

  void checkPhaseBinning() const ;
 private:
  void init() ;
} ;

class AsymmetryChi2 : public BinFlipChi2 {
 public :
  AsymmetryChi2(BinFlipParSet*, const HadronicParameters&, const TimeBinning&) ;
  AsymmetryChi2(BinFlipParSet*, const std::string&, const std::string&, const std::string& fname = "") ;
protected:
  /// Get the chi2 for the given time & phase bin.
  virtual double binChiSquared(unsigned, unsigned) const override;
};

class BinFlipChi2Simul : public BinFlipChi2Base {
 public :
  typedef std::deque<BinFlipChi2Base*> BinFlippers ;
  BinFlipChi2Simul(BinFlipChi2Base*, BinFlipChi2Base*) ;
  BinFlipChi2Simul(const BinFlippers&) ;

  virtual double getVal() override ;
  virtual void parametersChanged() override ;
 private :
  BinFlippers m_flippers ;
} ;

#endif
