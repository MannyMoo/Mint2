#include <Mint/HadronicParameters.h>
#include <TMath.h>
#include <Mint/NamedParameter.h>

using namespace std ;
using MINT::NamedParameter ;

int HadronicParameters::EventPair::binSign() const {
  return isPlus() ? 1 : -1 ;
}

bool HadronicParameters::EventPair::isPlus() const {
  return evt == evtPlus ;
}

HadronicParameters::BinSign::BinSign(HadronicParameters::ModelPtr model) :
  m_model(model)
{}

string HadronicParameters::BinSign::type() const {
  return "std" ;
}

HadronicParameters::EventPair HadronicParameters::BinSign::classifyEvent(const IDalitzEvent& evt) const {
  HadronicParameters::EventPair evtPair ;
  evtPair.evt = MINT::counted_ptr<IDalitzEvent>(new DalitzEvent(evt)) ;
  DalitzEvent* cpevt = new DalitzEvent(evt) ;
  cpevt->CP_conjugateYourself() ;
  evtPair.cpevt = MINT::counted_ptr<IDalitzEvent>(cpevt) ;
  double f = m_model->RealVal(*evtPair.evt) ;
  double fcp = m_model->RealVal(*evtPair.cpevt) ;
  if(f > fcp){
    evtPair.evtPlus = evtPair.evt ;
    evtPair.evtMinus = evtPair.cpevt ;
    evtPair.Fplus = f ;
    evtPair.Fminus = fcp ;
  }
  else{
    evtPair.evtPlus = evtPair.cpevt ;
    evtPair.evtMinus = evtPair.evt ;
    evtPair.Fplus = fcp ;
    evtPair.Fminus = f ;
  }
  return evtPair ;
}

HadronicParameters::Bin::Bin(double _Fplus, double _Fminus, const complex<double>& X,
			     double sumw, double sumw2, double norm) :
  m_Fplus(_Fplus * sumw * norm),
  m_Fminus(_Fminus * sumw * norm),
  m_X(X),
  m_sumw(sumw),
  m_sumw2(sumw2),
  m_norm(norm),
  m_model(0),
  m_cpmodel(0)
{
  m_X *= sqrt(Fplus() * Fminus()) * sumw * norm ;
}

HadronicParameters::Bin::Bin(ModelPtr model, ModelPtr cpmodel) :
  m_Fplus(0.),
  m_Fminus(0.),
  m_X(0., 0.),
  m_sumw(0.),
  m_sumw2(0.),
  m_norm(1.),
  m_model(model),
  m_cpmodel(cpmodel)
{}

HadronicParameters::Bin::Bin(const string& _name, unsigned number, const string& fname) :
  m_Fplus(0.),
  m_Fminus(0.),
  m_X(0., 0.),
  m_sumw(0.),
  m_sumw2(0.),
  m_norm(0.),
  m_model(0),
  m_cpmodel(0)
{
  string name = getName(_name, number) ;
  NamedParameter<double> _Fplus(name + "_Fplus", fname.c_str()) ;
  NamedParameter<double> _Fminus(name + "_Fminus", fname.c_str()) ;
  NamedParameter<double> Xplus_Re(name + "_Xplus_Re", fname.c_str()) ;
  NamedParameter<double> Xplus_Im(name + "_Xplus_Im", fname.c_str()) ;
  NamedParameter<double> sumw(name + "_sumw", fname.c_str()) ;
  NamedParameter<double> sumw2(name + "_sumw2", fname.c_str()) ;
  NamedParameter<double> norm(name + "_norm", fname.c_str()) ;

  m_norm = norm ;
  m_sumw = sumw ;
  m_sumw2 = sumw2 ;
  m_Fplus = _Fplus * m_sumw * m_norm ;
  m_Fminus = _Fminus * m_sumw * m_norm ;
  m_X = complex<double>(Xplus_Re, Xplus_Im) * sqrt(Fplus() * Fminus()) * m_sumw * m_norm ;
}


void HadronicParameters::Bin::add(const HadronicParameters::EventPair& evtpair,
				  double weight) {
  m_Fplus += evtpair.Fplus ;
  m_Fminus += evtpair.Fminus ;
  complex<double> Aplus = m_model->ComplexVal(*evtpair.evtPlus) ;
  complex<double> Aminus = m_cpmodel->ComplexVal(*evtpair.evtPlus) ;
  m_X += conj(Aplus) * Aminus ;
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
  return m_X / m_sumw / sqrt(Fplus() * Fminus()) / m_norm ;
}

complex<double> HadronicParameters::Bin::Xminus() const {
  return conj(Xplus()) ;
}

double HadronicParameters::Bin::norm() const {
  return m_norm ;
}

void HadronicParameters::Bin::setNorm(double norm) {
  m_norm = norm ;
}

void HadronicParameters::Bin::Print(const string& _name, unsigned number, ostream& os) const {
  string name = getName(_name, number) ;
  complex<double> Xp = Xplus() ;
  complex<double> Xm = Xminus() ;
  os.precision(16) ;
  os << name << "_Fplus\t" << Fplus() << endl ;
  os << name << "_Fminus\t" << Fminus() << endl ;
  os << name << "_Xplus_Re\t" << Xp.real() << endl ;
  os << name << "_Xplus_Im\t" << Xp.imag() << endl ;
  os << name << "_Xminus_Re\t" << Xm.real() << endl ;
  os << name << "_Xminus_Im\t" << Xm.imag() << endl ;
  os << name << "_sumw\t" << m_sumw << endl ;
  os << name << "_sumw2\t" << m_sumw2 << endl ;
  os << name << "_norm\t" << m_norm << endl ;
}

string HadronicParameters::Bin::getName(const string& name, unsigned number) {
  ostringstream sstr ;
  sstr << name << "_bin_" << number ;
  return sstr.str() ;
}

HadronicParameters::HadronicParameters(const HadronicParameters::Bins& bins, 
				       MINT::counted_ptr<IBinSign> binSign) :
  m_bins(bins),
  m_model(0),
  m_cpmodel(0),
  m_binSign(binSign)
{}

HadronicParameters::HadronicParameters(HadronicParameters::ModelPtr model,
				       HadronicParameters::ModelPtr cpmodel,
				       unsigned nbins,
				       MINT::counted_ptr<IBinSign> binSign) :
  m_bins(nbins, HadronicParameters::Bin(model, cpmodel)),
  m_model(model),
  m_cpmodel(cpmodel),
  m_binSign(binSign)
{}

HadronicParameters::HadronicParameters(const string& name, const string& fname) :
  m_bins(),
  m_model(0),
  m_cpmodel(0),
  m_binSign(0)
{
  NamedParameter<int> nbins(name + "_nBins", fname.c_str()) ;
  for(int i = 1 ; i < nbins + 1 ; ++i)
    m_bins.push_back(Bin(name, i, fname)) ;
  
  NamedParameter<int> evtPattern(name + "_eventPattern", fname.c_str()) ;
  m_model = ModelPtr(new FitAmpSum(DalitzEventPattern(evtPattern.getVector()), fname)) ;

  NamedParameter<string> binType(name + "_binType", fname.c_str()) ;
  m_binSign = getBinSign(m_model, binType) ;
}

unsigned HadronicParameters::nBins() const {
  return m_bins.size() ;
}

int HadronicParameters::binNumber(const HadronicParameters::EventPair& evtPair) const {
  complex<double> Amp = m_model->ComplexVal(*evtPair.evtPlus) ;
  complex<double> Ampbar = m_cpmodel->ComplexVal(*evtPair.evtPlus) ;
  double phasediff = arg(Amp/Ampbar) ;
  if(phasediff < 0.)
    phasediff += 2. * TMath::Pi() ;
  return (int(phasediff * nBins()/2./TMath::Pi() + 0.5) % nBins() + 1) * evtPair.binSign() ;
}

int HadronicParameters::binNumber(const IDalitzEvent& evt) const {
  return binNumber(m_binSign->classifyEvent(evt)) ;
}

void HadronicParameters::add(const HadronicParameters::EventPair& evtPair, double weight) {
  int binno = abs(binNumber(evtPair)) - 1;
  m_bins[binno].add(evtPair, weight) ;
}

void HadronicParameters::add(const IDalitzEvent& evt, double weight) {
  return add(m_binSign->classifyEvent(evt), weight) ;
}

void HadronicParameters::add(TRandom3& rndm, unsigned nevt) {
  for(unsigned i = 0 ; i < nevt ; ++i)
    add(rndm) ;
}

void HadronicParameters::add(TRandom3& rndm) {
  DalitzEvent evt(m_model->getAmpPtr(0)->getTreePattern(), (TRandom*)&rndm) ;
  add(evt) ;
}

double HadronicParameters::normalise() {
  for(auto& bin : m_bins)
    bin.setNorm(1.) ;
  double norm = 0. ;
  for(const auto& bin : m_bins)
    norm += bin.Fplus() + bin.Fminus() ;
  for(auto& bin : m_bins)
    bin.setNorm(norm) ;
  return norm ;
}

void HadronicParameters::Print(const string& name, ostream& os) const {
  unsigned i = 1 ;
  os << name << "_eventPattern" ;
  for(const auto& ipat : m_model->getAmpPtr(0)->getTreePattern().getVectorOfInts())
    os << "\t" << ipat ;
  os << endl ;
  os << name << "_nBins\t" << nBins() << endl ;
  os << name << "_binType\t" << m_binSign->type() << endl ;
  for(const auto& bin : m_bins){
    bin.Print(name, i, os) ;
    i += 1 ;
  }
}
  
MINT::counted_ptr<HadronicParameters::IBinSign> 
HadronicParameters::getBinSign(HadronicParameters::ModelPtr model, const string& type) {
  if(type == string("std"))
    return MINT::counted_ptr<IBinSign>(new HadronicParameters::BinSign(model)) ;
  throw invalid_argument("Unknown binning type: " + type) ;
}
