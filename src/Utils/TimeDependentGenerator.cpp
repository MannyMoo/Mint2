#include <Mint/TimeDependentGenerator.h>
#include <Mint/SplineGenerator.h>
#include <TRandom3.h>
#include <iostream>
#include <sys/stat.h>
#include <Mint/DalitzPdfSaveInteg.h>

using namespace std;
using namespace MINT;

// check if a file exists.
bool exists(const string& fname) {
  struct stat statinfo ;
  return stat(fname.c_str(), &statinfo) == 0 ;
}

TimeDependentGenerator::GenTimePoint::GenTimePoint(const double _decaytime, FitAmpSum* _model,
						   const double _integral, SignalGenerator* _generator) :
  decaytime(_decaytime),
  integral(_integral),
  model(_model),
  generator(_generator)
{}

// Take the CP conjugate of the head of the decay pattern.
DalitzEventPattern TimeDependentGenerator::anti(DalitzEventPattern pat) {
  pat[0].antiThis() ;
  return pat ;
}

/* Constructor, takes:
   name : the name of the generator and the directory in which the integrators will be saved.
   overwrite : whether to overwrite the existing integrator files (if they exist).
   rndm : The random number generator to use.
   precision : The precision to which the integrals must be calculated.
   pattern : The event pattern to be used (the CP conjugate will automatically be added).
   width : the decay width in 1/ps.
   deltam : the delta-mass in 1/ps.
   deltagamma : the delta-gamma in 1/ps.
   qoverp : the magnitude of q/p.
   phi : the phase of q/p.
   tmax : the maximum decay time that'll be generated.
   ntimepoints : the number of points to sample between 0 and tmax when building the generators.
   h_efficiency : (optional) histogram to which efficiency plot will be fitted
*/
TimeDependentGenerator::TimeDependentGenerator(const string& name, const bool overwrite, TRandom3* rndm,
					       double precision,
					       const DalitzEventPattern& pattern, double width, double deltam,
					       double deltagamma,
					       double qoverp, double phi, double tmax, int ntimepoints,
					       const bool saveIntegEvents, const double tmin, TH1F* h_efficiency) :
  m_name(name),
  m_rndm(rndm),
  m_pattern(pattern),
  m_cppattern(anti(pattern)),
  m_width(width),
  m_deltam(deltam),
  m_deltagamma(deltagamma),
  m_qoverp(qoverp),
  m_phi(phi),
  m_tmax(tmax),
  m_tmin(tmin),
  m_ntimepoints(ntimepoints),
  m_genmap(),
  m_timegenerators(),
  m_tagintegralfrac(0.),
  m_precision(precision),
  m_h_efficiency(h_efficiency)
{
  // If overwrite is true and the integrators directory exists, delete it.
  if(overwrite && exists(name)){
    cout << "Deleting previous integrators in directory " << name << endl ;
    string cmd("rm -rf " + name) ;
    system(cmd.c_str()) ;
  }

  if(!exists(name)){
    string cmd("mkdir -p " + name) ;
    system(cmd.c_str()) ;
  }
    
  const DalitzEventPattern* patterns[] = {&m_cppattern, &m_pattern} ;
  double sampleinterval = m_ntimepoints > 1 ? (m_tmax - m_tmin)/(m_ntimepoints-1) : 0. ;
  // Loop over flavours.
  for(int tag = -1 ; tag <= 1 ; tag += 2) {
    m_genmap[tag] = GenList() ;
    const DalitzEventPattern* evtpat = patterns[(tag+1)/2] ;
    const DalitzEventPattern* antipat = patterns[((tag+1)/2 + 1) % 2] ;
    // cout << "evtpat " ;
    // evtpat->print() ;
    // cout << " antipat " ;
    // antipat->print() ;
    // cout << endl ;
    vector<double> times ;
    vector<double> integrals ;
    // Loop over decay time sample points.
    for(int i = 0 ; i < m_ntimepoints ; ++i){
      double decaytime = m_tmin + i * sampleinterval ;
      AmpPair amps = amplitude_coefficients(tag, decaytime) ;
      FitAmpSum* model(new FitAmpSum(*evtpat)) ;
      *model *= amps.first ;
      if(amps.second != complex<double>(0., 0.)){
	FitAmpSum antimodel(*antipat) ;
	antimodel *= amps.second ;
	model->add(antimodel) ;
      }
      SignalGenerator* generator = new SignalGenerator(*evtpat, model) ;
      //cout << "Time point " << i << endl ;
      //model->print() ;
      //auto evt = generator->newEvent() ;
      //model->printAllAmps(*evt) ;
      ostringstream fname ;
      fname << m_name << "/tag_" << tag << "_decaytime_" << decaytime ;
      double integral(0.) ;
      // Calculate the integral if necessary.
      if(!exists(fname.str())){
	const string eventsFile(fname.str() + "_events.root") ;
	DalitzPdfSaveInteg dalitz(*evtpat, model, m_precision, fname.str(),
				  eventsFile, "topUp", fname.str()) ;
	integral = dalitz.getIntegralValue() ;
	dalitz.saveIntegrator(fname.str()) ;
	if(!saveIntegEvents){
	  string cmd("rm " + eventsFile) ;
	  system(cmd.c_str()) ;
	}
      }
      // Else retrive the integral from a file.
      else {
	auto intcalc = model->makeIntegrationCalculator() ;
	intcalc->retrieve(fname.str()) ;
	integral = intcalc->integral() ;
      }
      cout << "Make generator with tag " << tag << ", decay time " << decaytime
	   << ", coeffprod " << amps.first.real()
	   << " + " << amps.first.imag() << " j, "
	   << " coeffmix " << amps.second.real() << " + " << amps.second.imag()
	   << " j, integral " << integral << endl ;
      if(integral == 0.){
	complex<double> coeff = amps.first + amps.second ;
	integral = sqrt(coeff.real() * coeff.real() + coeff.imag() * coeff.imag()) ;
      }
      m_genmap[tag].push_back(GenTimePoint(decaytime, model, integral, generator)) ;
      times.push_back(decaytime) ;
      integrals.push_back(integral) ;
    }
    // Make a spline interpolator of the decay time distributions, to be used to generate
    // the decay times.
    ostringstream splinename ;
    splinename << "timespline_tag_" << tag ;
    string splinenamestr = splinename.str() ;
    TSpline3 timespline(splinenamestr.c_str(), &times[0], &integrals[0], times.size()) ;
    timespline.SetName(splinenamestr.c_str()) ;
    m_timegenerators.insert(make_pair(tag, SplineGenerator(rndm, timespline))) ;
  }
  // Calculate the integrated CP asymmetry.
  double integminus = m_timegenerators.find(-1)->second.integral() ;
  double integplus = m_timegenerators.find(1)->second.integral() ;
  cout << "Integrated CP asymmetry is " << (integplus - integminus)/(integplus + integminus) << endl ;
  m_tagintegralfrac = integminus/(integminus + integplus) ;
  double tauminus = m_timegenerators.find(-1)->second.mean() ;
  double tauplus = m_timegenerators.find(1)->second.mean() ;
  cout << "Tau minus: " << tauminus << endl ;
  cout << "Tau plus: " << tauplus << endl ;
  cout << "AGamma: " << (tauminus - tauplus)/(tauminus + tauplus) << endl ;

  TSpline3 efficiencyFit = TSpline3(m_h_efficiency) ;
}

// Get the coefficients of the amplitudes for the produced flavour and the mixed flavour
// given the tag and decay time.
TimeDependentGenerator::AmpPair 
TimeDependentGenerator::amplitude_coefficients(const int tag, const double decaytime) {
  double coeff = exp(-decaytime * 0.5 * (m_width + 0.5 * m_deltagamma)) ;
  complex<double> expterm = exp(complex<double>(0.5 * m_deltagamma * decaytime, m_deltam * decaytime)) ;
  complex<double> plusterm = 1. + expterm ;
  complex<double> minusterm = 1. - expterm ;
  complex<double> coeffprod = coeff * plusterm ;
  complex<double> coeffmix = pow(polar(m_qoverp, m_phi), tag) * coeff * minusterm ;
  return AmpPair(coeffprod, coeffmix) ;
}

// Generate a flavour.
int TimeDependentGenerator::generate_tag() const {
  double rndm = m_rndm->Rndm() ;
  if(rndm < m_tagintegralfrac)
    return -1 ;
  return 1 ;
}

// Generate a decay time for the given flavour.
double TimeDependentGenerator::generate_decay_time(const int tag) const {
  double decaytime = m_tmax + 1. ;
  while(decaytime > m_tmax)
    decaytime = m_timegenerators.find(tag)->second.gen_random() ;
 
  if(m_h_efficiency == NULL){ 
    return decaytime ; 
  }
  else{
    srand(1) ;

    int i = 0 ;
    float r_num = 0, efficiency = 0 ;
    int maxiter = 100000 ;
    while(true){
      efficiency = efficiencyFit.Eval(decaytime) ;
      r_num = rand() / RAND_MAX ;

      if(r_num < efficiency){
        return decaytime ;
      }
 
      decaytime =  m_timegenerators.find(tag)->second.gen_random() ;
      i += 1 ;
      if( i > maxiter ){ 
        cout << "WARNING: Decay time generation limit exceeded." << endl ; 
        return decaytime ;
      }

    }
  }
  return decaytime ; 
}

// Generate a Dalitz event for the given flavour and decay time.
MINT::counted_ptr<IDalitzEvent> TimeDependentGenerator::generate_dalitz_event(const int tag, const double decaytime) const {
  const GenList& genlist = m_genmap.find(tag)->second ;
  GenList::const_iterator igen = genlist.begin() ;
  while(igen->decaytime < decaytime && igen != genlist.end())
    ++igen ;
  if(igen == genlist.end()){
    cerr << "TimeDependentGenerator::generate_dalitz_event: ERROR: Got impossible decay time: "
	 << decaytime << " (tmax = " << m_tmax << ")" << endl ;
    return MINT::counted_ptr<IDalitzEvent>(0) ;
  }
  // Unlikely, but best to check.
  if(decaytime == igen->decaytime)
    return igen->generator->newEvent() ;
  GenList::const_iterator igenprev(igen) ;
  --igen ;
  // Pick between the generators either side of the decay time according to how close they are to it.
  double gensel = m_rndm->Rndm() ;
  if(gensel < 1. - (decaytime - igenprev->decaytime)/(igen->decaytime - igenprev->decaytime))
    return igenprev->generator->newEvent() ;
  return igen->generator->newEvent() ;
}

// Generate a flavour, decay time and Dalitz event.
TimeDependentGenerator::GenTimeEvent TimeDependentGenerator::generate_event() const {
  GenTimeEvent evt ;
  evt.tag = generate_tag() ;
  evt.decaytime = generate_decay_time(evt.tag) ;
  evt.evt = generate_dalitz_event(evt.tag, evt.decaytime) ;
  return evt ;
}

// Get the decay time generators.
const map<int, SplineGenerator> TimeDependentGenerator::time_generators() const {
  return m_timegenerators;
}
