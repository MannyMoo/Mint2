#ifndef __TIMEBINNING_H__
#define __TIMEBINNING_H__

#include <deque>
#include <string>
#include <ostream>
#include <Mint/HadronicParameters.h>
#include <TH1F.h>

class TFile ;

class TimeBinning {
 public :
  class Bin {
  public :
    Bin() ;
    Bin(const std::string&, unsigned, const std::string& fname = "") ;
    double t() const ;
    double t2() const ;
    double nPlus() const ;
    double nPlusErr() const ;
    double nMinus() const ;
    double nMinusErr() const ;
    bool usePoissonErrs() const ;
    void add(double, bool, double weight = 1.) ;
    void Print(const std::string&, unsigned, std::ostream& os = std::cout) const ;
    static std::string getName(const std::string&, unsigned) ;
    double chiSquared(double) const ;
  private :
    double m_t ;
    double m_t2 ;
    double m_sumw ;
    double m_sumw2 ;
    double m_sumwplus ;
    double m_sumw2plus ;
    double m_sumwminus ;
    double m_sumw2minus ;
  } ;
  typedef std::deque<Bin> Bins ;
  typedef std::deque<Bins> Bins2D ;

  TimeBinning(const std::vector<double>&, HadronicParameters::BinningPtr) ;
  TimeBinning(const std::string&, const std::string& fname = "") ;

  void add(IDalitzEvent&, int, double, double weight = 1.) ;
  HadronicParameters::BinningPtr phaseBinning() const ;
  int timeBin(double) const ;
  void Print(const std::string&, std::ostream& os = std::cout) const ;
  void write(const std::string&, const std::string&) const ;

  double chiSquared(unsigned, unsigned, double, double) const ;
  unsigned nBinsTime() const ;
  unsigned nBinsPhase() const ;
  const Bin& integratedBin(unsigned) const ;
  const Bin& bin(unsigned, unsigned) const ;
  const Bin& binBar(unsigned, unsigned) const ;

  std::deque<TH1F> plotVsTime(const std::string&, unsigned, int) const ;
  void savePlotsVsTime(const std::string&, TFile&) const ;
 private :
  std::vector<double> m_timeBins ;
  Bins2D m_bins ;
  Bins2D m_binsBar ;
  Bins m_binsInt ;
  HadronicParameters::BinningPtr m_phaseBinning ;

  Bin& _integratedBin(unsigned) ;
  Bin& _bin(unsigned, unsigned)  ;
  Bin& _binBar(unsigned, unsigned)  ;
} ;

#endif
