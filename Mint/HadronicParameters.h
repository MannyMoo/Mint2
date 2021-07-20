#ifndef __HADRONICPARAMETERS_H__
#define __HADRONICPARAMETERS_H__

#include <Mint/DalitzEvent.h>
#include <Mint/counted_ptr.h>
#include <Mint/FitAmpSum.h>
#include <Mint/FitParameter.h>
#include <complex>
#include <deque>
#include <utility>
#include <iostream>

class TRandom3 ;

namespace MINT{
  class MinuitParameterSet;
  class IMinimisable;
}

/** Class describing the hadronic parameters in bins of phase difference.
 */
class HadronicParameters {
 public :
  typedef MINT::counted_ptr<FitAmpSum> ModelPtr ;
  typedef MINT::counted_ptr<MINT::FitParameter> FitParPtr;
  
  /// Info on the bin of an Event - output of IPhaseBin::binInfo. Caches as much as possible for efficiency.
  class EventBinInfo {
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
    virtual EventBinInfo binInfo(IDalitzEvent&) const = 0 ;
    virtual std::string type() const = 0 ;
    virtual void Print(const std::string&, std::ostream& os = std::cout) const = 0 ;
    virtual ~PhaseBinningBase() {} ;
    const unsigned nBins ;
    int binNumber(IDalitzEvent&) const ;
    virtual bool operator==(const PhaseBinningBase&) const ;
    bool operator!=(const PhaseBinningBase&) const ;
  } ;
  /// Class for determining if an event lives in a +ve or -ve bin.
  class ModelPhaseBinning : public PhaseBinningBase {
  public :
    ModelPhaseBinning(ModelPtr, ModelPtr, unsigned) ;
    virtual ~ModelPhaseBinning() {} ;
    virtual EventBinInfo binInfo(IDalitzEvent&) const override ;
    virtual std::string type() const override ;
    virtual void Print(const std::string&, std::ostream& os = std::cout) const override ;
    static MINT::counted_ptr<PhaseBinningBase> 
      fromConfig(const std::string&, const std::string& fname = "") ;
    // TODO: implement this - not really sure how best to check that the models are the same.
    //virtual bool operator==(const PhaseBinningBase&) const override ;
  protected :
    virtual bool isFavoured(const double, const double, IDalitzEvent&) const ;
    ModelPtr m_model ;
    ModelPtr m_cpmodel ;
  } ;
  /// Phase binning class for 3-body decays using the line s13=s23 to determine favoured/suppressed.
  class ModelBinning3Body : public ModelPhaseBinning {
  public :
    ModelBinning3Body(ModelPtr, ModelPtr, unsigned) ;
    virtual ~ModelBinning3Body() {} ;
    virtual std::string type() const override ;
    static MINT::counted_ptr<PhaseBinningBase> 
      fromConfig(const std::string&, const std::string& fname = "") ;
  protected :
    virtual bool isFavoured(const double, const double, IDalitzEvent&) const override ;
  } ;
  /// Hadronic parameters in a bin of phase space.
  class Bin {
  private :
    /// Name of the bin
    std::string m_name;
    /// Magnitude sq. in the favoured region.
    FitParPtr m_Fplus = nullptr ;
    /// Magnitude sq. in the suppressed region.
    FitParPtr m_Fminus = nullptr ;
    /// Cross term.
    FitParPtr m_C = nullptr;
    FitParPtr m_S = nullptr;
    /// Magnitude sq. in the favoured region, for the CP-conjugate decay.
    FitParPtr m_Fbarplus = nullptr ;
    /// Magnitude sq. in the suppressed region, for the CP-conjugate decay.
    FitParPtr m_Fbarminus = nullptr ;
    /// Cross term, for the CP-conjugate decay.
    FitParPtr m_Cbar = nullptr;
    FitParPtr m_Sbar = nullptr;
    /// Sum of weights.
    double m_sumw = 0. ;
    /// Sum of weights sq.
    double m_sumw2 = 0.;
    /// The normalisation scale. - Are there potentially issues with having different normalisation for D0 and D0bar?
    double m_norm = 1.;
    /// The normalisation scale for the CP-conjugate decay.
    double m_normBar = 1.;

    /// Get zcp + dz
    std::complex<double> sumz(const std::complex<double>& z, const std::complex<double>& dz) const;
    /// Get zcp - dz
    std::complex<double> sumzBar(const std::complex<double>& z, const std::complex<double>& dz) const;

    /// Get the number of suppressed & favoured decays at the given time for the given mixing parameters.
    std::pair<double, double> _N(double, double, double, const std::complex<double>&, const std::complex<double>&,
				 double, double, const std::complex<double>&, const std::complex<double>&) const ;

    /// Get the expected ratio of events (suppressed)/(favoured) at the given time for the given mixing parameters.
    double _R(double, double, double, const std::complex<double>&, const std::complex<double>&,
	      double, double, const std::complex<double>&, const std::complex<double>&) const ;

    /// Initialise a fit parameter from an input file
    void initFitPar(FitParPtr&, FitParPtr, const std::string& fname,
		    const std::string& parName, const std::string& altParName = "");
    
  public :
    /// Initialise from predetermined parameters.
    Bin(const std::string&, unsigned, double, double, const std::complex<double>&,
	double, double, const std::complex<double>&,
	double norm = 1., double normbar = 1., double sumw = 1., double sumw2 = 1.) ;
    /// Initialise an empty bin.
    Bin(const std::string&, unsigned) ;
    /// Initialise from a config file.
    Bin(const std::string&, unsigned, const std::string& fname) ;

    /// Set Abar == A
    void setNoCPV();
    /// Check if CPV is allowed
    bool allowsCPV() const;

    /// Get fit parameters
    std::vector<MINT::FitParameter*> getFitPars();
    /// Get fit parameters
    std::vector<const MINT::FitParameter*> getFitPars() const;
    
    /// Set the parameter set
    void setParSet(MINT::MinuitParameterSet*);
    
    /// Add a DalitzEvent and its conjugate.
    void add(const EventBinInfo&, const EventBinInfo&, double weight = 1.) ;
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
    /// Set the normalisation for both the decay and CP-conjugate decay.
    void normalise(double norm = 1., double normBar = 1);
    /// Call when the sum over MC events is finished when calculating pars
    void finaliseSum();
    /// Print the parameters.
    void Print(const std::string& name, unsigned, std::ostream& os=std::cout) const ;
    /// Get the name string.
    static std::string getName(const std::string&, unsigned) ;

    /// Get expected ratio of events (suppressed)/(favoured) at the given time for the given mixing parameters.
    double R(double, double, double, const std::complex<double>&, const std::complex<double>&) const ;
    /** Get expected ratio of events (suppressed)/(favoured) at the given time for the given mixing parameters,
	for the CP-conjugate decay.*/
    double Rbar(double, double, double, const std::complex<double>&, const std::complex<double>&) const ;

    /// Get the expected number of suppressed & favoured decays at the given time for the given mixing parameters.
    std::pair<double, double> N(double, double, double, const std::complex<double>&,
				const std::complex<double>&) const ;
    /** Get the expected number of suppressed & favoured decays at the given time for the given mixing parameters
	for the CP-conjugate decay.*/
    std::pair<double, double> Nbar(double, double, double, const std::complex<double>&,
				   const std::complex<double>&) const ;
    
    /** Get the asymmetry in the number of decays between D0 & D0bar at the given time for the given mixing
	parameters.*/
    double asymmetry(double, double, double, const std::complex<double>&,
		     const std::complex<double>&) const ;
  } ;

  /// A set of bins.
  typedef std::deque<Bin> Bins ;
  /// Pointer to the binning scheme.
  typedef MINT::counted_ptr<PhaseBinningBase> BinningPtr ;
  typedef MINT::MinuitParameterSet::CovMatrix CovMatrix;
  
  /// Initialise from a predetermined set of bins.
  HadronicParameters(const std::string&, const Bins&, BinningPtr) ;
  /// Initialise from a binning scheme.
  HadronicParameters(const std::string&, BinningPtr) ;
  /// Initialise from a config file.
  HadronicParameters(const std::string&, const std::string& fname = "") ;

  /// Get the binning scheme.
  const PhaseBinningBase& binning() const ;
  /// Get a pointer to the binning scheme.
  HadronicParameters::BinningPtr binningPtr() const ;

  /// Get the bin for a DalitzEvent.
  const Bin& bin(IDalitzEvent&) const ;
  /// Get a bin by number.
  const Bin& bin(unsigned) const ;

  /// Add a DalitzEvent.
  void add(IDalitzEvent&, double weight = 1.) ;
  /// Add a number of DalitzEvents generated from flat phase space.
  void add(const DalitzEventPattern&, TRandom3&, unsigned) ;
  /// Add a DalitzEvent generated from flat phase space.
  void add(const DalitzEventPattern&, TRandom3&) ;

  /// Get the integral over phase space.
  std::pair<double, double> integral() const ;

  /// Call when the sum to generate hadronic pars from MC is finished
  void finaliseSum();
  /// Normalise the parameters.
  std::pair<double, double> normalise(double norm = 1., double normBar = 1.) ;
  /// Print the parameters.
  void Print(std::ostream& os = std::cout) const ;
  /// Write to a file.
  void write(const std::string&) const ;

  /// Get the PhaseBinning type.
  static BinningPtr getPhaseBinning(const std::string&, const std::string& fname = "") ;

  /// Set the parameter set
  void setParSet(MINT::MinuitParameterSet*);

  /// Check if CPV is allowed
  bool allowsCPV() const;

  /// Get the Gaussian constraints and add the parameters to the parset
  MINT::IMinimisable* getConstraints(MINT::MinuitParameterSet*);

  /// Set the covariance matrix
  void setCovMatrix(const CovMatrix& m);

  /// Get floating parameters
  std::vector<const MINT::FitParameter*> getFloatingPars() const;
  
 private :
  std::string m_name;
  Bins m_bins ;
  BinningPtr m_phaseBinning ;
  CovMatrix m_covMatrix = CovMatrix();

  /// Get floating parameters
  std::vector<MINT::FitParameter*> getFloatingPars();

} ;

#endif
