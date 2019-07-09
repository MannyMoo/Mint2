#include "Mint/IMinimisable.h"
#include "Mint/MinuitParameterSet.h"
#include "Mint/Minimisable.h"
#include "Mint/FitParameter.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "cmath"
#include "complex"

using namespace MINT;
using namespace std;


class binflipChi2 : public Minimisable{
    unsigned int m_nbinsPhase;
    unsigned int m_nbinsTime;
    vector<complex<float> > m_X;
    vector<float> m_r;
    vector<float> m_tAv;
    vector<float> m_tSqAv;
    TH2F m_pHistD0;
    TH2F m_pHistD0bar;
    TH2F m_nHistD0;
    TH2F m_nHistD0bar;
    FitParameter m_ReZcp;
    FitParameter m_ImZcp;
    FitParameter m_ReDz;
    FitParameter m_ImDz;
    bool m_fakeData;

  public:
    binflipChi2(vector<complex<float> > X, vector<float> r, vector<float> tAv, vector<float> tSqAv, 
                       TH2F pHistD0, TH2F pHistD0bar, TH2F nHistD0, TH2F nHistD0bar, float ReZcp, float ImZcp, 
                       float ReDz, float ImDz, float stepSize, bool fakeData = false, vector<float> Fm = vector<float>(), 
                       vector<float> Fp = vector<float>());
    ~binflipChi2();
    double getVal();
    vector<vector<TGraph> > getFits(complex<float> zcp, complex<float> deltaz);
    void genFakeData(vector<float> Fm, vector<float> Fp);
};

