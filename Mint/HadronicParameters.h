#ifndef __HADRONICPARAMETERS_H__
#define __HADRONICPARAMETERS_H__

#include <Mint/DalitzEvent.h>
#include <Mint/counted_ptr.h>
#include <Mint/FitAmpSum.h>
#include <complex>
#include <deque>
#include <utility>

class TRandom3 ;

/** Class describing the hadronic parameters in bins of phase difference.
 */
class HadronicParameters {
 public :
  typedef MINT::counted_ptr<FitAmpSum> ModelPtr ;

  /// Event and its CP conjugate - output of BinSign::classifyEvent.
  class EventPair {
  public :
    MINT::counted_ptr<IDalitzEvent> evt ;
    MINT::counted_ptr<IDalitzEvent> cpevt ;
    MINT::counted_ptr<IDalitzEvent> evtPlus ;
    MINT::counted_ptr<IDalitzEvent> evtMinus ;
    double Fplus ;
    double Fminus ;

    /// Get the bin sign, +1 or -1.
    int binSign() const ;
    /// Check if the event is in a +ve bin.
    bool isPlus() const ;
  } ;
  /// Interface class for determining if an event lives in a +ve or -ve bin.
  class IBinSign {
  public :
    virtual EventPair classifyEvent(const IDalitzEvent&) const = 0 ;
    virtual ~IBinSign() {} ;
  } ;
  /// Class for determining if an event lives in a +ve or -ve bin.
  class BinSign : public IBinSign {
  public :
    BinSign(ModelPtr) ;
    virtual ~BinSign() {} ;
    virtual EventPair classifyEvent(const IDalitzEvent&) const override ;
  private :
    ModelPtr m_model ;
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
    /// Sum of weights.
    double m_sumw ;
    /// Sum of weights sq.
    double m_sumw2 ;
    /// The normalisation scale.
    double m_norm ;

    /// The model.
    ModelPtr m_model ;
    /// The conjugate model.
    ModelPtr m_cpmodel ;
  public :
    /// Initialise from predetermined parameters.
    Bin(double, double, const std::complex<double>&, double sumw = 1., double sumw2 = 1., double norm = 1.) ;
    /// Initialise from a model and its conjugate.
    Bin(ModelPtr, ModelPtr) ;

    /// Add a DalitzEvent and its conjugate.
    void add(const HadronicParameters::EventPair&, double weight = 1.) ;
    /// Get the magnitude sq in the favoured region.
    double Fplus() const ;
    /// Get the magnitude sq in the suppressed region.
    double Fminus() const ;
    /// Get the cross term for the favoured region.
    std::complex<double> Xplus() const ;
    /// Get the cross term for the suppressed region.
    std::complex<double> Xminus() const ;
    /// Get the normalisation.
    double norm() const ;
    /// Set the normalisation.
    void setNorm(double) ;
    /// Print the parameters.
    void Print() const ;
  } ;

  /// A set of bins.
  typedef std::deque<Bin> Bins ;
  /// Initialise from a predetermined set of bins.
  HadronicParameters(const DalitzEventPattern&, const DalitzEventPattern&,
		     const Bins&, MINT::counted_ptr<IBinSign>) ;
  /// Initialise from a model and its conjugate, and a number of bins.
  HadronicParameters(const DalitzEventPattern&, const DalitzEventPattern&,
		     ModelPtr, ModelPtr, unsigned, MINT::counted_ptr<IBinSign>) ;

  /// Get the number of bins.
  unsigned nBins() const ;

  /// Get the bin number for a DalitzEvent.
  int binNumber(const IDalitzEvent&) const ;
  /// Get the bin number for an EventPair.
  int binNumber(const EventPair&) const ;
  /// Get the bin for a DalitzEvent.
  Bin& bin(IDalitzEvent&) const ;
  /// Get the bin for an EventPair.
  Bin& bin(const EventPair&) const ;

  /// Add a DalitzEvent.
  void add(const IDalitzEvent&, double weight = 1.) ;
  /// Add an EventPair.
  void add(const EventPair&, double weight = 1.) ;
  /// Add a number of DalitzEvents generated from flat phase space.
  void add(TRandom3&, unsigned) ;
  /// Add a DalitzEvent generated from flat phase space.
  void add(TRandom3&) ;

  /// Normalise the parameters.
  double normalise() ;
  /// Print the parameters.
  void Print() const ;

 private :
  DalitzEventPattern m_pat ;
  DalitzEventPattern m_cppat ;
  Bins m_bins ;
  ModelPtr m_model ;
  ModelPtr m_cpmodel ;
  MINT::counted_ptr<IBinSign> m_binSign ;
} ;

#endif
