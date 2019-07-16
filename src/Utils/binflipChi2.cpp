#include "Mint/binflipChi2.h"

using namespace MINT;
using namespace std;

binflipChi2::binflipChi2(vector<complex<float> > X, vector<float> r, vector<float> tAv, vector<float> tSqAv, TH2F pHistD0, 
			 TH2F pHistD0bar, TH2F nHistD0, TH2F nHistD0bar, float ReZcp, float ImZcp, float ReDz, float ImDz, float stepSize,
                         int fakeData, vector<float> Fm, vector<float> Fp):

  Minimisable(new MinuitParameterSet),
  m_X(X),
  m_r(r),
  m_tAv(tAv),
  m_tSqAv(tSqAv),
  m_pHistD0(pHistD0),
  m_pHistD0bar(pHistD0bar),
  m_nHistD0(nHistD0),
  m_nHistD0bar(nHistD0bar),
  m_ReZcp("ReZcp", FitParameter::FIT, ReZcp, stepSize, 0, 0, getParSet()),
  m_ImZcp("ImZcp", FitParameter::FIT, ImZcp, stepSize, 0, 0, getParSet()),
  m_ReDz("ReDz", FitParameter::FIT, ReDz, stepSize, 0, 0, getParSet()),
  m_ImDz("ImDz", FitParameter::FIT, ImDz, stepSize, 0, 0, getParSet()),
  m_fakeData(fakeData),
  m_Fm(Fm),
  m_Fp(Fp)
{
    m_nbinsPhase = m_r.size();
    m_nbinsTime = m_tAv.size();

    if (( m_fakeData > 0 ) && ( Fm.size() == m_nbinsPhase ) && ( Fp.size() == m_nbinsPhase )){
        genFakeData(); 
    }
}


binflipChi2::~binflipChi2(){
    delete getParSet();
}


double binflipChi2::getChi2(const int tag, const int j, const int b, const TH2F& nHist, const TH2F& pHist,
			    const complex<float>& zcp, const complex<float>& deltaz) const {
  if(nHist.GetBinContent(j, b) == 0 || pHist.GetBinContent(j, b) == 0)
    return 0. ;
  const double R = getRVal(tag, j, b, zcp, deltaz) ;
  double den =  pow(nHist.GetBinError(j,b), 2) + pow(pHist.GetBinError(j,b) * R, 2); 
  if(den == 0.)
    return 0. ;

  double num = pow(nHist.GetBinContent(j,b) - pHist.GetBinContent(j,b) * R, 2);
  return num/den ;
} 

double binflipChi2::getRVal(const int tag, const int j, const int b, const complex<float>& zcp,
			    const complex<float>& deltaz) const {
  double numerator = 0., denominator = 0.;
  
  int i = tag == -1 ? 1 : 0 ;

  numerator = m_r[b-1] * ( 1 + 0.25 * m_tSqAv[j-1] * ( pow(zcp,2) - pow(deltaz,2) ).real() );
  numerator += 0.25 * m_tSqAv[j-1] * pow(abs(zcp + (float)pow(-1, i) * deltaz), 2);
  numerator += sqrt(m_r[b-1]) * m_tAv[j-1] * ( conj(m_X[b-1]) * (zcp + (float)pow(-1, i) * deltaz) ).real();
  
  denominator = 1 + 0.25 * m_tSqAv[j-1] * ( pow(zcp,2) - pow(deltaz,2) ).real();
  denominator += m_r[b-1] * 0.25 * m_tSqAv[j-1] * pow(abs(zcp + (float)pow(-1, i) * deltaz), 2);
  denominator += sqrt(m_r[b-1]) * m_tAv[j-1] * ( m_X[b-1] * (zcp + (float)pow(-1, i) * deltaz) ).real();
  
  double Rval = numerator/denominator;
  return Rval ;
}

double binflipChi2::getVal(){    

    double chi2 = 0.;

    complex<float> zcp(m_ReZcp.getCurrentFitVal(), m_ImZcp.getCurrentFitVal());
    complex<float> deltaz(m_ReDz.getCurrentFitVal(), m_ImDz.getCurrentFitVal());

    for(int b = 1; b <= m_nbinsPhase; b++){
    
        for(int j = 1; j <= m_nbinsTime; j++){

	  double D0_term = getChi2(1, j, b, m_nHistD0, m_pHistD0, zcp, deltaz) ;
	  double D0bar_term = getChi2(-1, j, b, m_nHistD0bar, m_pHistD0bar, zcp, deltaz) ;
   
	  chi2 += D0_term + D0bar_term;
        }
    }
    return chi2;
}


vector<vector<TGraph> > binflipChi2::getFits(){

    vector<vector<TGraph> > fits(2);

    complex<float> zcp(m_ReZcp.getCurrentFitVal(), m_ImZcp.getCurrentFitVal());
    complex<float> deltaz(m_ReDz.getCurrentFitVal(), m_ImDz.getCurrentFitVal());
    
    for(int i : {0,1}){

        for(int b = 1; b <= m_nbinsPhase; b++){
	    TGraph temp = TGraph(m_nbinsTime);
	    fits[i].push_back( temp );

            for(int j = 1; j <= m_nbinsTime; j++){
	      double Rval = getRVal(i == 0 ? 1 : -1, j, b, zcp, deltaz) ;
	      fits[i][b-1].SetPoint(j-1, m_tAv[j-1], Rval);
            }
        }
    }
    
    return fits;
}


void binflipChi2::genFakeData(){
    TRandom3 rndm = TRandom3(m_fakeData);
    complex<float> zcp(m_ReZcp.getCurrentFitVal(), m_ImZcp.getCurrentFitVal());
    complex<float> deltaz(m_ReDz.getCurrentFitVal(), m_ImDz.getCurrentFitVal());

    complex<float> qoverp = sqrt( (zcp + deltaz) / (zcp - deltaz) );
    complex<float> z = sqrt( pow(zcp,2) - pow(deltaz,2) );
        
    for(int i : {0,1}){

        complex<float> pqterm;
        if ( i == 0 ){
	    pqterm = ((float)-1)*qoverp;
	}
	else{
	    pqterm = ((float)-1)*pow(qoverp, -1); 
        }


        for(int b = 1; b <= m_nbinsPhase; b++){

     	    for(int j = 1; j <= m_nbinsTime; j++){


                double pval = 0., nval = 0., rNum = 0., err = 0.;

                pval = m_Fp[b-1] * ( 1 + 0.25 * m_tSqAv[j-1] * ( pow(z, 2) ).real() );
                pval += 0.25 * m_tSqAv[j-1] * pow(abs(z), 2) * pow(abs(pqterm), 2) * m_Fm[b-1];
                pval += m_tAv[j-1] * sqrt(m_Fm[b-1] * m_Fp[b-1]) * (pqterm * m_X[b-1] * z).real();

                nval = m_Fm[b-1] * ( 1 + 0.25 * m_tSqAv[j-1] * ( pow(z, 2) ).real() );
                nval += 0.25 * m_tSqAv[j-1] * pow(abs(z), 2) * pow(abs(pqterm), 2) * m_Fp[b-1];
                nval += m_tAv[j-1] * sqrt(m_Fm[b-1] * m_Fp[b-1]) * (pqterm * conj(m_X[b-1]) * z).real();

                if( i == 0 ){
		    if( (m_pHistD0.GetBinContent(j,b) != 0) && (m_pHistD0.GetBinError(j, b) != 0) ){
			err = pval * m_pHistD0.GetBinError(j, b) / m_pHistD0.GetBinContent(j,b);
		    }
                    else{
		        err = sqrt(pval);
		        //err = pval * 0.01;
		    }
                    err = sqrt(pval);//testing
                    if (m_pHistD0.GetBinContent(j,b) != 0){
			rNum = rndm.Gaus(0, err);
			m_pHistD0.SetBinContent(j, b, pval + rNum);
			m_pHistD0.SetBinError(j, b, err);
		    }
                         
		    if( (m_nHistD0.GetBinContent(j,b) != 0) && (m_nHistD0.GetBinError(j, b) != 0) ){
			err = nval * m_nHistD0.GetBinError(j, b) / m_nHistD0.GetBinContent(j,b);
		    }
                    else{
		        err = sqrt(nval);
                        //err = nval * 0.01;
		    }
                    err = sqrt(nval);//testing
                    if (m_nHistD0.GetBinContent(j,b) != 0){
			rNum = rndm.Gaus(0, err);
			m_nHistD0.SetBinContent(j, b, nval + rNum);
			m_nHistD0.SetBinError(j, b, err);
		    }
		}
                else{
		    if( (m_pHistD0bar.GetBinContent(j,b) != 0) && (m_pHistD0bar.GetBinError(j, b) != 0) ){
			err = pval * m_pHistD0bar.GetBinError(j, b) / m_pHistD0bar.GetBinContent(j,b);
		    }
                    else{
		        err = sqrt(pval);
                        //err = pval * 0.01;
		    }
                    err = sqrt(pval);//testing
                    if( m_pHistD0bar.GetBinContent(j,b) != 0 ){
			rNum = rndm.Gaus(0, err);
			m_pHistD0bar.SetBinContent(j, b, pval + rNum);
			m_pHistD0bar.SetBinError(j, b, err);
		    }
                         
		    if( (m_nHistD0bar.GetBinContent(j,b) != 0) && (m_nHistD0bar.GetBinError(j, b) != 0) ){
			err = nval * m_nHistD0bar.GetBinError(j, b) / m_nHistD0bar.GetBinContent(j,b);
		    }
                    else{
		        err = sqrt(nval);
                        //err = nval * 0.01;
		    }
                    err = sqrt(nval);//testing
                    if( m_nHistD0bar.GetBinContent(j,b) != 0 ){
			rNum = rndm.Gaus(0, err);
			m_nHistD0bar.SetBinContent(j, b, nval + rNum);
			m_nHistD0bar.SetBinError(j, b, err);
		    }
		}
	    }
	}
    }       
}

