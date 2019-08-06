#include <Mint/HadronicParameters.h>
#include <TMath.h>
#include <Mint/NamedParameter.h>

using namespace std ;
using MINT::NamedParameter ;

int HadronicParameters::EventBinNumber::binSign() const {
  return binNumber > 0 ? 1 : -1 ;
}

bool HadronicParameters::EventBinNumber::isPlus() const {
  return binNumber > 0 ;
}

HadronicParameters::PhaseBinningBase::PhaseBinningBase(unsigned nBins) :
  nBins(nBins)
{}

HadronicParameters::ModelPhaseBinning::ModelPhaseBinning(HadronicParameters::ModelPtr model,
							 HadronicParameters::ModelPtr cpmodel,
							 unsigned nBins) :
  PhaseBinningBase(nBins),
  m_model(model),
  m_cpmodel(cpmodel)
{}

string HadronicParameters::ModelPhaseBinning::type() const {
  return "model" ;
}

HadronicParameters::EventBinNumber HadronicParameters::ModelPhaseBinning::binNumber(IDalitzEvent& evt) const {
  EventBinNumber binNo ;
  binNo.evt = &evt ;
  binNo.amp = m_model->ComplexVal(evt) ;
  binNo.ampBar = m_cpmodel->ComplexVal(evt) ;
  binNo.F = norm(binNo.amp) ;
  binNo.Fbar = norm(binNo.ampBar) ;
  double phasediff = arg(binNo.amp/binNo.ampBar) ;
  if(phasediff < 0.)
    phasediff += 2. * TMath::Pi() ;
  binNo.binNumber = (int(phasediff * nBins/2./TMath::Pi() + 0.5) % nBins + 1) ;
  // Possible issue if F == Fbar ?
  if(binNo.Fbar > binNo.F)
    binNo.binNumber *= -1 ;
  return binNo ;
}

void HadronicParameters::ModelPhaseBinning::Print(const std::string& _name, std::ostream& os) const {
  string name = _name + "_binning_" ;
  os << name << "type\t" << type() << endl ;
  os << name << "nBins\t" << nBins << endl ;
  os << name << "eventPattern" ;
  for(const auto& ipat : m_model->getAmpPtr(0)->getTreePattern().getVectorOfInts())
    os << "\t" << ipat ;
  os << endl ;
  os << name << "cpEventPattern" ;
  for(const auto& ipat : m_cpmodel->getAmpPtr(0)->getTreePattern().getVectorOfInts())
    os << "\t" << ipat ;
  os << endl ;
}

HadronicParameters::BinningPtr HadronicParameters::ModelPhaseBinning::fromConfig(const string& _name,
										 const string& fname) {
  string name = _name + "_binning_ " ;
  if(string(NamedParameter<string>(name + "type", fname.c_str())) != "model")
    return BinningPtr(0) ;
  NamedParameter<int> nBins(name + "nBins", fname.c_str()) ;
  NamedParameter<int> pat(name + "eventPattern", fname.c_str()) ;
  NamedParameter<int> cpPat(name + "cpEventPattern", fname.c_str()) ;
  ModelPtr model(new FitAmpSum(DalitzEventPattern(pat), fname)) ;
  ModelPtr cpModel(new FitAmpSum(DalitzEventPattern(cpPat), fname)) ;
  return BinningPtr(new ModelPhaseBinning(model, cpModel, nBins)) ;
}

HadronicParameters::Bin::Bin(double _Fplus, double _Fminus, const complex<double>& X,
			     double _Fbarplus, double _Fbarminus, const complex<double>& Xbar,
			     double norm, double normBar, double sumw, double sumw2) :
  m_Fplus(_Fplus * sumw * norm),
  m_Fminus(_Fminus * sumw * norm),
  m_X(X),
  m_Fbarplus(_Fbarplus * sumw * normBar),
  m_Fbarminus(_Fbarminus * sumw * normBar),
  m_Xbar(Xbar),
  m_sumw(sumw),
  m_sumw2(sumw2),
  m_norm(norm),
  m_normBar(normBar)
{
  m_X *= sqrt(Fplus() * Fminus()) * sumw * norm ;
  m_Xbar *= sqrt(Fbarplus() * Fbarminus()) * sumw * normBar ;
}

HadronicParameters::Bin::Bin() :
  m_Fplus(0.),
  m_Fminus(0.),
  m_X(0., 0.),
  m_Fbarplus(0.),
  m_Fbarminus(0.),
  m_Xbar(0., 0.),
  m_sumw(0.),
  m_sumw2(0.),
  m_norm(1.),
  m_normBar(1.)
{}

HadronicParameters::Bin::Bin(const string& _name, unsigned number, const string& fname) :
  m_Fplus(0.),
  m_Fminus(0.),
  m_X(0., 0.),
  m_Fbarplus(0.),
  m_Fbarminus(0.),
  m_Xbar(0., 0.),
  m_sumw(0.),
  m_sumw2(0.),
  m_norm(0.),
  m_normBar(1.)
{
  string name = getName(_name, number) ;
  NamedParameter<double> _Fplus(name + "_Fplus", fname.c_str()) ;
  NamedParameter<double> _Fminus(name + "_Fminus", fname.c_str()) ;
  NamedParameter<double> Xplus_Re(name + "_Xplus_Re", fname.c_str()) ;
  NamedParameter<double> Xplus_Im(name + "_Xplus_Im", fname.c_str()) ;

  NamedParameter<double> _Fbarplus(name + "_Fbarplus", fname.c_str()) ;
  NamedParameter<double> _Fbarminus(name + "_Fbarminus", fname.c_str()) ;
  NamedParameter<double> Xbarplus_Re(name + "_Xbarplus_Re", fname.c_str()) ;
  NamedParameter<double> Xbarplus_Im(name + "_Xbarplus_Im", fname.c_str()) ;

  NamedParameter<double> sumw(name + "_sumw", fname.c_str()) ;
  NamedParameter<double> sumw2(name + "_sumw2", fname.c_str()) ;
  NamedParameter<double> norm(name + "_norm", fname.c_str()) ;
  NamedParameter<double> normBar(name + "_normBar", fname.c_str()) ;

  m_norm = norm ;
  m_normBar = normBar ;
  m_sumw = sumw ;
  m_sumw2 = sumw2 ;

  m_Fplus = _Fplus * m_sumw * m_norm ;
  m_Fminus = _Fminus * m_sumw * m_norm ;
  m_X = complex<double>(Xplus_Re, Xplus_Im) * sqrt(Fplus() * Fminus()) * m_sumw * m_norm ;

  m_Fbarplus = _Fbarplus * m_sumw * m_normBar ;
  m_Fbarminus = _Fbarminus * m_sumw * m_normBar ;
  m_Xbar = complex<double>(Xbarplus_Re, Xbarplus_Im) * sqrt(Fbarplus() * Fbarminus()) * m_sumw * m_normBar ;
}


void HadronicParameters::Bin::add(const HadronicParameters::EventBinNumber& evtPlus,
				  const HadronicParameters::EventBinNumber& evtMinus,
				  double weight) {
  m_Fplus += evtPlus.F * weight ;
  m_Fminus += evtMinus.F * weight ;
  m_Fbarplus += evtMinus.Fbar * weight ;
  m_Fbarminus += evtPlus.Fbar * weight ;
  m_X += conj(evtPlus.amp) * evtPlus.ampBar * weight ;
  m_Xbar += conj(evtMinus.ampBar) * evtMinus.amp * weight ;
  m_sumw += weight ;
  m_sumw2 += weight*weight ;
}

double HadronicParameters::Bin::Fplus() const {
  return m_Fplus / m_sumw / m_norm ;
}

double HadronicParameters::Bin::Fminus() const {
  return m_Fminus / m_sumw / m_norm ;
}

complex<double> HadronicParameters::Bin::Xplus() const {
  return m_X / m_sumw / sqrt(Fplus() * Fbarminus()) / m_norm ;
}

complex<double> HadronicParameters::Bin::Xminus() const {
  return conj(Xplus()) ;
}

double HadronicParameters::Bin::Fbarplus() const {
  return m_Fbarplus / m_sumw / m_normBar ;
}

double HadronicParameters::Bin::Fbarminus() const {
  return m_Fbarminus / m_sumw / m_normBar ;
}

complex<double> HadronicParameters::Bin::Xbarplus() const {
  return m_Xbar / m_sumw / sqrt(Fbarplus() * Fminus()) / m_normBar ;
}

complex<double> HadronicParameters::Bin::Xbarminus() const {
  return conj(Xbarplus()) ;
}

double HadronicParameters::Bin::norm() const {
  return m_norm ;
}

void HadronicParameters::Bin::setNorm(double norm) {
  m_norm = norm ;
}

double HadronicParameters::Bin::normBar() const {
  return m_normBar ;
}

void HadronicParameters::Bin::setNormBar(double norm) {
  m_normBar = norm ;
}

void HadronicParameters::Bin::Print(const string& _name, unsigned number, ostream& os) const {
  string name = getName(_name, number) ;
  complex<double> Xp = Xplus() ;
  complex<double> Xbarp = Xbarplus() ;
  os.precision(16) ;
  os << name << "_Fplus\t" << Fplus() << endl ;
  os << name << "_Fminus\t" << Fminus() << endl ;
  os << name << "_Xplus_Re\t" << Xp.real() << endl ;
  os << name << "_Xplus_Im\t" << Xp.imag() << endl ;
  os << name << "_Fbarplus\t" << Fbarplus() << endl ;
  os << name << "_Fbarminus\t" << Fbarminus() << endl ;
  os << name << "_Xbarplus_Re\t" << Xbarp.real() << endl ;
  os << name << "_Xbarplus_Im\t" << Xbarp.imag() << endl ;
  os << name << "_sumw\t" << m_sumw << endl ;
  os << name << "_sumw2\t" << m_sumw2 << endl ;
  os << name << "_norm\t" << m_norm << endl ;
  os << name << "_normBar\t" << m_normBar << endl ;
}

string HadronicParameters::Bin::getName(const string& name, unsigned number) {
  ostringstream sstr ;
  sstr << name << "_bin_" << number ;
  return sstr.str() ;
}

HadronicParameters::HadronicParameters(const HadronicParameters::Bins& bins, 
				       BinningPtr phaseBinning) :
  m_bins(bins),
  m_phaseBinning(phaseBinning)
{}

HadronicParameters::HadronicParameters(BinningPtr phaseBinning) :
  m_bins(phaseBinning->nBins, HadronicParameters::Bin()),
  m_phaseBinning(phaseBinning)
{}

HadronicParameters::HadronicParameters(const string& name, const string& fname) :
  m_bins(),
  m_phaseBinning(HadronicParameters::getPhaseBinning(name, fname))
{
  for(int i = 1 ; i < m_phaseBinning->nBins + 1 ; ++i)
    m_bins.push_back(Bin(name, i, fname)) ;  
}

int HadronicParameters::binNumber(IDalitzEvent& evt) const {
  return m_phaseBinning->binNumber(evt).binNumber ;
}

const HadronicParameters::Bin& HadronicParameters::bin(IDalitzEvent& evt) const {
  int binNo = abs(binNumber(evt)) ;
  return m_bins.at(binNo-1) ;
}

void HadronicParameters::add(IDalitzEvent& evt, double weight) {
  EventBinNumber binNo = m_phaseBinning->binNumber(evt) ;
  DalitzEvent cpEvt(evt) ;
  cpEvt.CP_conjugateYourself() ;
  EventBinNumber cpBinNo = m_phaseBinning->binNumber(cpEvt) ;
  if(binNo.isPlus())
    m_bins.at(binNo.binNumber-1).add(binNo, cpBinNo, weight) ;
  else 
    m_bins.at(cpBinNo.binNumber-1).add(cpBinNo, binNo, weight) ;
}

void HadronicParameters::add(const DalitzEventPattern& pat, TRandom3& rndm, unsigned nevt) {
  for(unsigned i = 0 ; i < nevt ; ++i)
    add(pat, rndm) ;
}

void HadronicParameters::add(const DalitzEventPattern& pat, TRandom3& rndm) {
  DalitzEvent evt(pat, (TRandom*)&rndm) ;
  add(evt) ;
}

pair<double, double> HadronicParameters::integral() const {
  double norm = 0. ;
  double normBar = 0. ;
  for(const auto& bin : m_bins){
    norm += bin.Fplus() + bin.Fminus() ;
    normBar += bin.Fbarplus() + bin.Fbarminus() ;
  }
  return pair<double, double>(norm, normBar) ;
}

pair<double, double> HadronicParameters::normalise(double _norm, double _normBar) {
  for(auto& bin : m_bins){
    bin.setNorm(1.) ;
    bin.setNormBar(1.) ;
  }
  pair<double, double> norm = integral() ;
  norm.first /= _norm ;
  norm.second /= _normBar ;
  for(auto& bin : m_bins){
    bin.setNorm(norm.first) ;
    bin.setNormBar(norm.second) ;
  }
  return norm ;
}

void HadronicParameters::Print(const string& name, ostream& os) const {
  m_phaseBinning->Print(name, os) ;
  unsigned i = 1 ;
  for(const auto& bin : m_bins){
    bin.Print(name, i, os) ;
    i += 1 ;
  }
}
  
HadronicParameters::BinningPtr
HadronicParameters::getPhaseBinning(const string& name, const string& fname) {
  NamedParameter<string> type(name + "_binning_type", fname.c_str()) ;
  string strtype(type) ;
  if(strtype == string("model"))
    return BinningPtr(ModelPhaseBinning::fromConfig(name, fname)) ;
  throw invalid_argument("Unknown binning type: " + strtype) ;
}
