#ifndef __TIMEBINNING_H__
#define __TIMEBINNING_H__

#include <deque>
#include <string>
#include <ostream>
#include <Mint/HadronicParameters.h>
#include <TH1F.h>
#include <Mint/counted_ptr.h>

class TFile ;
class TSpline3;

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
    static void poissonErrors2(const double, double&, double&);
    void errors2(double&, double&, double&, double&) const ;
    double chiSquared(double) const ;
    Bin operator+(Bin) const;
    Bin& operator+=(const Bin&);
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

  TimeBinning(const std::vector<double>&, HadronicParameters::BinningPtr, double,
	      const TH1* hefficiency = nullptr) ;
  TimeBinning(const std::string&, const std::string& fname = "") ;
  virtual ~TimeBinning() {};
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

  /// Get the asymmetry of D0 to D0bar in the given time & phase bin.
  std::pair<double, double> asymmetry(unsigned, unsigned) const;
  /// Get the chi-squared of the asymmetry in the given time & phase bin, given the expected value
  double asymmetryChiSquared(unsigned, unsigned, double) const;

  std::deque<TH1F> plotVsTime(const std::string&, unsigned, int) const ;
  std::deque<std::deque<TH1F> > plotsVsTime(const std::string&) const ;
  void savePlotsVsTime(const std::string&, TFile&) const ;

  double meanUnmixedTime(unsigned) const ;
  double meanUnmixedTime2(unsigned) const ;
  void setLifetime(double) ;
  double getLifetime() const ;
  bool isConsistent(const TimeBinning&) const;
  TimeBinning operator+(TimeBinning) const;
  TimeBinning& operator+=(const TimeBinning&);
 private :
  std::vector<double> m_timeBins ;
  std::vector<double> m_meant ;
  std::vector<double> m_meant2 ;
  double m_lifetime ;
  Bins2D m_bins ;
  Bins2D m_binsBar ;
  Bins m_binsInt ;
  HadronicParameters::BinningPtr m_phaseBinning ;

  Bin& _integratedBin(unsigned) ;
  Bin& _bin(unsigned, unsigned)  ;
  Bin& _binBar(unsigned, unsigned)  ;

  MINT::counted_ptr<TSpline3> m_efficiencySpline;

  virtual double unmixedTimeMoment(unsigned, double, int) const ;
} ;

#endif
