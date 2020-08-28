#ifndef __XYCONSTRAINTCHI_H__
#define __XYCONSTRAINTCHI_H__

#include <Mint/GaussianConstraintChi2.h>

class BinFlipParSet;

class XYConstraintChi2 : public GaussianConstraintChi2 {
 private:
  virtual TVectorT<double> getDiffs() override;
  BinFlipParSet* m_binflipParSet;
  static MINT::MinuitParameterSet* makeParSet(BinFlipParSet*);
 public:
  XYConstraintChi2(BinFlipParSet*);
  virtual ~XYConstraintChi2();
};

#endif
