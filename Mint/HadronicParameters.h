#ifndef __HADRONICPARAMETERS_H__
#define __HADRONICPARAMETERS_H__

#include <Mint/DalitzEvent.h>
#include <Mint/counted_ptr.h>
#include <Mint/FitAmpSum.h>
#include <complex>
#include <deque>
#include <utility>
#include <iostream>

class TRandom3 ;

/** Class describing the hadronic parameters in bins of phase difference.
 */
class HadronicParameters {
 public :
  typedef MINT::counted_ptr<FitAmpSum> ModelPtr ;

  /// Info on the bin number of an Event - output of IPhaseBin::binNumber. Caches as much as possible for efficiency.
  class EventBinNumber {
  public :
    // Should I use a counted_ptr?
    IDalitzEvent* evt ;
    std::complex<double> amp ;
    std::complex<double> ampBar ;
    double F ;
    double Fbar ;
    int binNumber ;

    /// Get the bin sign, +1 or -1.
    int binSign() const ;
    /// Check if the event is in a +ve bin.
    bool isPlus() const ;
  } ;
  /// Interface class for determining which bin an event lives in.
  class PhaseBinningBase {
  public :
    PhaseBinningBase(unsigned nBins) ;
    virtual EventBinNumber binNumber(IDalitzEvent&) const = 0 ;
    virtual std::string type() const = 0 ;
    virtual void Print(const std::string&, std::ostream& os = std::cout) const = 0 ;
    virtual ~PhaseBinningBase() {} ;
    const unsigned nBins ;
  } ;
  /// Class for determining if an event lives in a +ve or -ve bin.
  class ModelPhaseBinning : public PhaseBinningBase {
  public :
    ModelPhaseBinning(ModelPtr, ModelPtr, unsigned) ;
    virtual ~ModelPhaseBinning() {} ;
    virtual EventBinNumber binNumber(IDalitzEvent&) const override ;
    virtual std::string type() const override ;
    virtual void Print(const std::string&, std::ostream& os = std::cout) const override ;
    static MINT::counted_ptr<PhaseBinningBase> 
      fromConfig(const std::string&, const std::string& fname = "") ;
  private :
    ModelPtr m_model ;
    ModelPtr m_cpmodel ;
  } ;

  /// Hadronic parameters in a bin of phase space.
  class Bin {
  private :
    /// Magnitude sq. in the favoured region.
    double m_Fplus ;
    /// Magnitude sq. in the suppressed region.
    double m_Fminus ;
    /// Cross term.
    std::complex<double> m_X ;
    /// Magnitude sq. in the favoured region, for the CP-conjugate decay.
    double m_Fbarplus ;
    /// Magnitude sq. in the suppressed region, for the CP-conjugate decay.
    double m_Fbarminus ;
    /// Cross term, for the CP-conjugate decay.
    std::complex<double> m_Xbar ;
    /// Sum of weights.
    double m_sumw ;
    /// Sum of weights sq.
    double m_sumw2 ;
    /// The normalisation scale. - Are there potentially issues with having different normalisation for D0 and D0bar?
    double m_norm ;
    /// The normalisation scale for the CP-conjugate decay.
    double m_normBar ;

    /// Get the expected ratio of events (suppressed)/(favoured) at the given time for the given mixing parameters.
    double _R(double, double, const std::complex<double>&, const std::complex<double>&,
	      double, double, const std::complex<double>&, const std::complex<double>&) const ;
  public :
    /// Initialise from predetermined parameters.
    Bin(double, double, const std::complex<double>&, double, double, const std::complex<double>&,
	double norm = 1., double normbar = 1., double sumw = 1., double sumw2 = 1.) ;
    /// Initialise an empty bin.
    Bin() ;
    /// Initialise from a config file.
    Bin(const std::string&, unsigned, const std::string& fname = "") ;

    /// Add a DalitzEvent and its conjugate.
    void add(const EventBinNumber&, const EventBinNumber&, double weight = 1.) ;
    /// Get the magnitude sq in the favoured region.
    double Fplus() const ;
    /// Get the magnitude sq in the suppressed region.
    double Fminus() const ;
    /// Get the cross term for the favoured region.
    std::complex<double> Xplus() const ;
    /// Get the cross term for the suppressed region.
    std::complex<double> Xminus() const ;
    /// Get the magnitude sq in the favoured region, for the CP-conjugate decay.
    double Fbarplus() const ;
    /// Get the magnitude sq in the suppressed region, for the CP-conjugate decay.
    double Fbarminus() const ;
    /// Get the cross term for the favoured region, for the CP-conjugate decay.
    std::complex<double> Xbarplus() const ;
    /// Get the cross term for the suppressed region, for the CP-conjugate decay.
    std::complex<double> Xbarminus() const ;
    /// Get the normalisation.
    double getNorm() const ;
    /// Set the normalisation.
    void setNorm(double) ;
    /// Get the normalisation, for the CP-conjugate decay.
    double getNormBar() const ;
    /// Set the normalisation, for the CP-conjugate decay.
    void setNormBar(double) ;
    /// Print the parameters.
    void Print(const std::string& name, unsigned, std::ostream& os=std::cout) const ;
    /// Get the name string.
    static std::string getName(const std::string&, unsigned) ;

    /// Get expected ratio of events (suppressed)/(favoured) at the given time for the given mixing parameters.
    double R(double, double, const std::complex<double>&, const std::complex<double>&) const ;
    /** Get expected ratio of events (suppressed)/(favoured) at the given time for the given mixing parameters,
	for the CP-conjugate decay.*/
    double Rbar(double, double, const std::complex<double>&, const std::complex<double>&) const ;
    
  } ;

  /// A set of bins.
  typedef std::deque<Bin> Bins ;
  /// Pointer to the binning scheme.
  typedef MINT::counted_ptr<PhaseBinningBase> BinningPtr ;

  /// Initialise from a predetermined set of bins.
  HadronicParameters(const Bins&, BinningPtr) ;
  /// Initialise from a binning scheme.
  HadronicParameters(BinningPtr) ;
  /// Initialise from a config file.
  HadronicParameters(const std::string&, const std::string& fname = "") ;

  /// Get the binning scheme.
  const PhaseBinningBase& binning() const ;

  /// Get the bin number for a DalitzEvent.
  int binNumber(IDalitzEvent&) const ;
  /// Get the bin for a DalitzEvent.
  const Bin& bin(IDalitzEvent&) const ;

  /// Add a DalitzEvent.
  void add(IDalitzEvent&, double weight = 1.) ;
  /// Add a number of DalitzEvents generated from flat phase space.
  void add(const DalitzEventPattern&, TRandom3&, unsigned) ;
  /// Add a DalitzEvent generated from flat phase space.
  void add(const DalitzEventPattern&, TRandom3&) ;

  /// Get the integral over phase space.
  std::pair<double, double> integral() const ;

  /// Normalise the parameters.
  std::pair<double, double> normalise(double norm = 1., double normBar = 1.) ;
  /// Print the parameters.
  void Print(const std::string&, std::ostream& os = std::cout) const ;

  /// Get the PhaseBinning type.
  static BinningPtr getPhaseBinning(const std::string&, const std::string& fname = "") ;

 private :
  Bins m_bins ;
  BinningPtr m_phaseBinning ;
} ;

#endif
