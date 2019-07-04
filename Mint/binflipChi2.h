#include "Mint/IMinimisable.h"
#include "Mint/MinuitParameterSet.h"
#include "Mint/Minimisable.h"
#include "TH2F.h"
#include "TGraph.h"
#include "cmath"
#include "complex"

using namespace MINT;
using namespace std;


class binflipChi2 : public Minimisable{
    int m_nbinsPhase;
    int m_nbinsTime;
    vector<complex<float> > m_X;
    float* m_r;
    float* m_tAv;
    float* m_tSqAv;
    TH2F m_pHistD0;
    TH2F m_pHistD0bar;
    TH2F m_nHistD0;
    TH2F m_nHistD0bar;

  public:
    binflipChi2(int nbinsPhase, int nbinsTime, float* Xreal, float* Ximag, float* r, float* tAv, float* tSqAv,
		      TH2F pHistD0, TH2F pHistD0bar, TH2F nHistD0, TH2F nHistD0bar);
    double getVal();
    vector<vector<TGraph> > getFits(complex<float> zcp, complex<float> deltaz);
};

