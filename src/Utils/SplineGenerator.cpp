#include <Mint/SplineGenerator.h>
#include <TRandom3.h>
#include <iostream>

using namespace std ;

SplineGenerator::SplineGenerator(TRandom3* rndm, const TSpline3& spline) :
  m_rndm(rndm),
  m_spline(spline),
  m_integral(0.),
  m_boxintegral(0.),
  m_mean(0.),
  m_bins(spline.GetNp())
{
  BinInfo ibin ;
  double xmin(0.), a(0.), b(0.), c(0.), d(0.) ;
  int i = m_spline.GetNp() - 1 ;
  // Get info on the highest knot in the spline.
  m_spline.GetCoeff(i, xmin, a, b, c, d) ;
  ibin.xmin = xmin ;
  ibin.xmax = 1e30 ;
  ibin.integral = 0. ;
  ibin.ymin = 0. ;
  ibin.ymax = 1e30 ;
  m_bins[i] = ibin ;
  --i ;
  // Loop backwards over knots, so we can use the x value of the knot above
  // to set the range.
  for( ; i >= 0 ; --i){
    m_spline.GetCoeff(i, xmin, a, b, c, d) ;
    ibin.xmin = xmin ;
    ibin.xmax = m_bins[i+1].xmin ;
    // cout << "ibin " << i << " xmin " << ibin.xmin << " xmax " << ibin.xmax << endl ;
    // cout << "a " << a << " b " << b << " c " << c << " d " << d << endl ;
    ibin.integral = integral(i, ibin.xmin, ibin.xmax) ;
    // cout << "integral " << ibin.integral << endl ;
    m_integral += ibin.integral ;
    pair<double, double> tps = turning_points(i) ;
    // cout << "tps.first " << tps.first << " tps.second " << tps.second << endl ;
    tps.first = max(min(tps.first, ibin.xmax), ibin.xmin) ;
    tps.second = max(min(tps.second, ibin.xmax), ibin.xmin) ;
    // cout << "tps.first " << tps.first << " tps.second " << tps.second << endl ;
    ibin.ymax = max(
		    max(
			max(m_spline.Eval(ibin.xmin),
			    m_spline.Eval(ibin.xmax)),
			m_spline.Eval(tps.first)),
		    m_spline.Eval(tps.second)) ;
    ibin.ymin = min(
		    min(
			min(m_spline.Eval(ibin.xmin),
			    m_spline.Eval(ibin.xmax)),
			m_spline.Eval(tps.first)),
		    m_spline.Eval(tps.second)) ;
    if(ibin.ymin < 0.){
      cerr << "SplineGenerator ERROR: ymin < 0. (" << ibin.ymin << ") for bin " << i << endl ;
    }
    ibin.boxintegral = ibin.ymax * (ibin.xmax - ibin.xmin) ;
    // cout << "ymax " << ibin.ymax << " ymin " << ibin.ymin << " boxintegral " << ibin.boxintegral << endl ;
    m_boxintegral += ibin.boxintegral ;
    m_bins[i] = ibin ;
    m_mean += mean(i, ibin.xmin, ibin.xmax) ;
  }
  if(m_integral != 0.)
    m_mean /= m_integral ;
}

// Total integral between the minimum and maximum knots.
double SplineGenerator::integral() const {
  return m_integral ;
}

// Mean of the PDF.
double SplineGenerator::mean() const {
  return m_mean ;
}
  
// Generate a random number from the spline.
double SplineGenerator::gen_random() const {
  while(true){
    double boxsel = m_rndm->Rndm() * m_boxintegral ;
    double boxsum(0.) ;
    vector<BinInfo>::const_iterator ibin = m_bins.begin() ;
    for(; ibin != m_bins.end() ; ++ibin){
      boxsum += ibin->boxintegral ;
      if(boxsum >= boxsel)
	break ;
    }
    double x = m_rndm->Rndm() * (ibin->xmax - ibin->xmin) + ibin->xmin ;
    if(m_rndm->Rndm() * ibin->ymax < m_spline.Eval(x))
      return x ;
  }
}

// Get the integral of the spline between xmin and xmax.
double SplineGenerator::integral(double xmin, double xmax) {
  int istart = m_spline.FindX(xmin) ;
  int iend = m_spline.FindX(xmax) ;
  double _integral = integral(istart, xmin, m_bins[istart].xmax)
    + integral(iend, m_bins[iend].xmin, xmax) ;
  for(int ibin = istart + 1 ; ibin != iend ; ++ibin)
    _integral += m_bins[ibin].integral ;
  return _integral ;
}

// Get the spline.
const TSpline3& SplineGenerator::spline() const {
  return m_spline ;
}

// Partial integral for knot i at x.
double SplineGenerator::partial_integral(int i, double x) {
  double xmin(0.), a(0.), b(0.), c(0.), d(0.) ;
  m_spline.GetCoeff(i, xmin, a, b, c, d) ;
  x -= xmin ;
  return x * (a + x * (b/2. + x * (c/3. + x * d/4.))) ;
}

// Integral for knot i between xmin and xmax.
double SplineGenerator::integral(int i, double xmin, double xmax) {
  return partial_integral(i, xmax) - partial_integral(i, xmin) ;
}

// Get the turning points of knot i.
pair<double, double> SplineGenerator::turning_points(int i) const {
  // Note the reversal of coefficient labelling, as the standard
  // quadratic formula is a x^2 + b x + c while TSplinePoly
  // uses a + b x + c x^2 + d x^3.
  double xmin(0.), d(0.), c(0.), b(0.), a(0.) ;
  m_spline.GetCoeff(i, xmin, d, c, b, a) ;
  b *= 2. ;
  a *= 3. ;
  double arg = b*b - 4. * a * c ;
  if(arg < 0.)
    return pair<double, double>(xmin-1e30, xmin-1e30) ;
  return pair<double, double>(xmin + (-b - sqrt(arg))/2./a, xmin + (-b + sqrt(arg))/2./a) ;
}

// The partial integral of x * f(x) used for calculating the mean.
double SplineGenerator::mean_part_integral(int i, double x) {
  double xmin(0.), a(0.), b(0.), c(0.), d(0.) ;
  m_spline.GetCoeff(i, xmin, a, b, c, d) ;
  double xshift = x - xmin ;
  return x * xshift * (a + xshift * (b/2. + xshift * (c/3. + xshift * d/4.)))
    - xshift * xshift * (a/2. + xshift * (b/6. + xshift * (c/12. + xshift * d/20.))) ;
}

// Get the mean between xmin and xmax for knot i.
double SplineGenerator::mean(int i, double xmin, double xmax) {
  return mean_part_integral(i, xmax) - mean_part_integral(i, xmin) ;
}
