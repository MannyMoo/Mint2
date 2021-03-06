#include "Mint/binflipChi2.h"

using namespace MINT;
using namespace std;

binflipChi2::binflipChi2(vector<complex<double> > X, vector<double> r, vector<double> tAv, vector<double> tSqAv, TH2F pHistD0, 
			 TH2F pHistD0bar, TH2F nHistD0, TH2F nHistD0bar, double ReZcp, double ImZcp, double ReDz, double ImDz, double stepSize,
                         int fakeData, vector<double> Fm, vector<double> Fp, int verbosity):

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
  m_Fp(Fp),
  m_verbosity(verbosity)
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


double binflipChi2::getVal(){    

    double chi2 = 0.;
    int count = 0;
    vector<vector<TGraph> > fits(2);

    fits = getFits();

    for(int b = 1; b <= m_nbinsPhase; b++){
    
        double* Rvals_pl = fits[0][b-1].GetY();
        double* Rvals_mi = fits[1][b-1].GetY();  
    
        for(int j = 1; j <= m_nbinsTime; j++){

	    double R_pl = Rvals_pl[j-1];
            double R_mi = Rvals_mi[j-1]; 
            double D0_term = 0, D0bar_term = 0;
            double D0_num = 0, D0_den = 0, D0bar_num = 0, D0bar_den = 0; 

            if(( m_nHistD0.GetBinContent(j,b) != 0 )&&( m_pHistD0.GetBinContent(j,b) != 0 )){
		D0_num = pow(m_nHistD0.GetBinContent(j,b) - m_pHistD0.GetBinContent(j,b) * R_pl, 2);
		D0_den =  pow(m_nHistD0.GetBinError(j,b), 2) + pow(m_pHistD0.GetBinError(j,b) * R_pl, 2); 
		if(D0_den != 0){
		    D0_term = D0_num/D0_den;
		}
	    } 
     
            if(( m_nHistD0bar.GetBinContent(j,b) != 0 )&&( m_pHistD0bar.GetBinContent(j,b) != 0 )){
		D0bar_num = pow(m_nHistD0bar.GetBinContent(j,b) - m_pHistD0bar.GetBinContent(j,b) * R_mi, 2);
		D0bar_den =  pow(m_nHistD0bar.GetBinError(j,b), 2) + pow(m_pHistD0bar.GetBinError(j,b) * R_mi, 2); 
		if(D0bar_den != 0){
		    D0bar_term = D0bar_num/D0bar_den;
		}
	    }
   
            chi2 += D0_term + D0bar_term; 
       
            float tol = 0.01;
            if(D0_term > tol){
	        count += 1;
	    }
            if(D0bar_term > tol){
	        count += 1;
	    }
        }
    }

    if( m_verbosity > 0 ){
        cout << "Number of terms contributing to chi2 is: " << count << endl;
    }

    return chi2;
}


vector<vector<TGraph> > binflipChi2::getFits(){

    vector<vector<TGraph> > fits(2);

    complex<double> zcp(m_ReZcp.getCurrentFitVal(), m_ImZcp.getCurrentFitVal());
    complex<double> deltaz(m_ReDz.getCurrentFitVal(), m_ImDz.getCurrentFitVal());
    
    for(int i : {0,1}){

        for(int b = 1; b <= m_nbinsPhase; b++){
	    TGraph temp = TGraph(m_nbinsTime);
	    fits[i].push_back( temp );

            for(int j = 1; j <= m_nbinsTime; j++){
	        double numerator = 0., denominator = 0.;
     
                numerator = m_r[b-1] * ( 1 + 0.25 * m_tSqAv[j-1] * ( pow(zcp,2) - pow(deltaz,2) ).real() );
                numerator += 0.25 * m_tSqAv[j-1] * pow(abs(zcp + (double)pow(-1, i) * deltaz), 2);
    	        numerator += sqrt(m_r[b-1]) * m_tAv[j-1] * ( conj(m_X[b-1]) * (zcp + (double)pow(-1, i) * deltaz) ).real();

                denominator = 1 + 0.25 * m_tSqAv[j-1] * ( pow(zcp,2) - pow(deltaz,2) ).real();
                denominator += m_r[b-1] * 0.25 * m_tSqAv[j-1] * pow(abs(zcp + (double)pow(-1, i) * deltaz), 2);
    	        denominator += sqrt(m_r[b-1]) * m_tAv[j-1] * ( m_X[b-1] * (zcp + (double)pow(-1, i) * deltaz) ).real();

                double Rval = numerator/denominator;
    	        fits[i][b-1].SetPoint(j-1, m_tAv[j-1], Rval);

            }
        }
    }
    
    return fits;
}


void binflipChi2::genFakeData(){
    TRandom3 rndm = TRandom3(m_fakeData);
    complex<double> zcp(m_ReZcp.getCurrentFitVal(), m_ImZcp.getCurrentFitVal());
    complex<double> deltaz(m_ReDz.getCurrentFitVal(), m_ImDz.getCurrentFitVal());

    complex<double> qoverp = sqrt( (zcp + deltaz) / (zcp - deltaz) );
    complex<double> z = sqrt( pow(zcp,2) - pow(deltaz,2) );
        
    for(int i : {0,1}){

        complex<double> pqterm;
        if ( i == 0 ){
	    pqterm = ((double)-1)*qoverp;
	}
	else{
	    pqterm = ((double)-1)*pow(qoverp, -1); 
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
		    }
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
		    }
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
		    }
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
		    }
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


TGraph binflipChi2::getFit(int i, int b){

    TGraph fit = TGraph(m_nbinsTime);
    fit = getFits()[i][b];
    return fit;

}
