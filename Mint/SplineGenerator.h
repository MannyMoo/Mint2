#ifndef __SPLINEGENERATOR_H__
#define __SPLINEGENERATOR_H__

#include <TSpline.h>
#include <vector>
class TRandom3 ;

// Class to generate random numbers according to a spline interpolator.
class SplineGenerator {
 private :
  struct BinInfo {
    double xmin ;
    double xmax ;
    double integral ;
    double ymin ;
    double ymax ;
    double boxintegral ;
  } ;
  
  TRandom3* m_rndm ;
  mutable TSpline3 m_spline ; // Cause ROOT is crap at const correctness.
  double m_integral ;
  double m_boxintegral ;
  double m_mean ;
  std::vector<BinInfo> m_bins ;
  
 public :
  SplineGenerator(TRandom3* rndm, const TSpline3& spline) ;
  
  // Total integral between the minimum and maximum knots.
  double integral() const ;

  // Mean of the PDF.
  double mean() const ;

  // Generate a random number from the spline.
  double gen_random() const ;

  // Get the integral of the spline between xmin and xmax.
  double integral(double xmin, double xmax) ;

  // Get the spline.
  const TSpline3& spline() const ;

 private :
  // Partial integral for knot i at x.
  double partial_integral(int i, double x) ;

  // Integral for knot i between xmin and xmax.
  double integral(int i, double xmin, double xmax) ;
 
  // Get the turning points of knot i.
  std::pair<double, double> turning_points(int i) const ;
  
  // The partial integral of x * f(x) used for calculating the mean.
  double mean_part_integral(int i, double x) ;

  // Get the mean between xmin and xmax for knot i.
  double mean(int i, double xmin, double xmax) ;
} ;

#endif
