#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "haploid.h"

#include <iostream>
using namespace std;

HapInd::HapInd(void)
{
	m_phi = NULL;
	m_psi = NULL;
	m_mgt = NULL; 
	snp_dstr = NULL; 
	m_ploid = 1; 
	m_admix = NULL; 
	m_ps = NULL; 
	m_pk = NULL; 
}
                 
HapInd::HapInd(int n)
{
	nLoci = n; 
	m_phi = NULL;
	m_psi = NULL;
	m_mgt = NULL; 
	snp_dstr = NULL; 
	m_ploid = 1; 
	m_admix = NULL; 
	m_ps = NULL; 
	m_pk = NULL; 
}

HapInd::~HapInd(void)
{
	if(m_phi) {Free2DMatrix(m_phi);    m_phi = NULL; }
	if(m_psi) {Free2DMatrix(m_psi);       m_psi = NULL;}
}

double HapInd::prG(double t, short snp)
{
	double res = 0.0; 
	double tt = (m_mu + (1.0 - 2.0 * m_mu) * t); 
	switch(snp)
	{
		case 0:
			res = 1.0 - tt;
			break;
		case 1:
			res = tt;
			break;
		default:
			res = 1.0;
			break;
	}
	return res; 
}

void HapInd::clean_after_em(int deep)
{
	delete[] m_pjk; m_pjk = NULL; 
	delete[] m_pjs; m_pjs = NULL; 
	Free2DMatrix(m_phi); m_phi = NULL; 
	Free2DMatrix(m_psi); m_psi = NULL; 
	delete[] phiScale; phiScale = NULL; 
	delete[] psiScale; psiScale = NULL; 
	if(deep == 1) 
	{
		Free2DMatrix(m_pk); m_pk = NULL; 
		Free2DMatrix(m_ps); m_ps = NULL; 
	}
}

void HapInd::compute_ps(const int nS, const int nK, ModelParam * pMP, const int mode)
{
	int nsk = nS * nK; 
	double * prz = new double[nsk]; 
	if(m_ps == NULL) 
		m_ps = Allocate2DMatrix(nLoci, nS); 

    CalcAll(nS, nK, pMP, NULL, NULL, NULL, NULL, 1, 1); 
	for (int m = 0; m < nLoci; m++)
	{
		for (int j = 0; j < nsk; j++)
			prz[j] = m_phi[m][j] * m_psi[m][j]; 
		normalize(prz, nsk); 

		for (int j = 0, s = 0; s < nS; s++)
		{
			m_ps[m][s] = 0; 
			for (int k = 0; k < nK; k++, j++)
                m_ps[m][s] += prz[j]; 
		}

	}
	delete[] prz; 
}

void HapInd::compute_ps(const int nS, const int nK, ModelParam * pMP, const int mode, const int b, vector<int> pbreak)
{
	int nsk = nS * nK; 
	double * prz = new double[nsk]; 
	int len = pbreak.at(b+1) - pbreak.at(b); 
	if(m_ps != NULL) {
		Free2DMatrix(m_ps); m_ps = NULL; 
	}
	m_ps = Allocate2DMatrix(len, nS); 

	for (int m = 0; m < len; m++)
		for (int s = 0; s < nS; s++)
			m_ps[m][s] = 0; 
   	forward_backward_block(nS, nK, pMP, NULL, NULL, NULL, NULL, 1, b, pbreak); 
	for (int m = 0; m < len; m++)
	{
		for (int j = 0; j < nsk; j++)
			prz[j] = m_phi[m][j] * m_psi[m][j]; 
		normalize(prz, nsk); 

		for (int s = 0; s < nS; s++)
			m_ps[m][s] = 0; 
		for (int j = 0, s = 0; s < nS; s++)
			for (int k = 0; k < nK; k++, j++)
                m_ps[m][s] += prz[j]; 
	}
	clean_after_em(0); 
	delete[] prz; 
}

void HapInd::compute_pk(const int nS, const int nK, ModelParam * pMP, const int mode)
{
	int nsk = nS * nK; 
	double * prz = new double[nsk]; 
	if(m_pk == NULL) 
		m_pk = Allocate2DMatrix(nLoci, nK);
	
    CalcAll(nS, nK, pMP, NULL, NULL, NULL, NULL, 1, 1); 
	for (int m = 0; m < nLoci; m++)
	{
		for (int j = 0; j < nsk; j++)
			prz[j] = m_phi[m][j] * m_psi[m][j]; 
		normalize(prz, nsk); 

		for (int k = 0; k < nK; k++)
			m_pk[m][k] = 0; 
		for (int j = 0, s = 0; s < nS; s++)
			for (int k = 0; k < nK; k++, j++)
                m_pk[m][k] += prz[j]; 
	}
	delete[] prz; 
}

void HapInd::compute_pk(const int nS, const int nK, ModelParam * pMP, const int mode, const int b, vector<int> pbreak)
{
	int nsk = nS * nK; 
	double * prz = new double[nsk]; 
	int len = pbreak.at(b+1) - pbreak.at(b); 
	if(m_pk != NULL) {
		Free2DMatrix(m_pk); m_pk = NULL; 
	}
	m_pk = Allocate2DMatrix(len, nK); 

	for (int m = 0; m < len; m++)
		for (int k = 0; k < nK; k++)
			m_pk[m][k] = 0; 
   	forward_backward_block(nS, nK, pMP, NULL, NULL, NULL, NULL, 1, b, pbreak); 
	for (int m = 0; m < len; m++)
	{
		for (int j = 0; j < nsk; j++)
			prz[j] = m_phi[m][j] * m_psi[m][j]; 
		normalize(prz, nsk); 

		for (int k = 0; k < nK; k++)
			m_pk[m][k] = 0; 
		for (int j = 0, s = 0; s < nS; s++)
			for (int k = 0; k < nK; k++, j++)
                m_pk[m][k] += prz[j]; 
	}
	clean_after_em(0); 
	delete[] prz; 
}

void HapInd::CalcAll(const int nS, const int nK, ModelParam * pMD, double ** stop, double ** sbot, double ** sjmsk, double ** sjms, int n1, int n2)
{
	if(m_admix == NULL) 
	{
		m_admix = new double[nS]; 
		for (int s = 0; s < nS; s++)
			m_admix[s] = 1.0 / nS; 
	}
	
	double * admix = m_admix; 

	double * rk = pMD->Getrk(); 
	double * rs = pMD->Getrs(); 
	double ** theta = pMD->Gettheta(); 
	double *** beta = pMD->Getbeta(); 

	int nsk = nS * nK; 
	double ts = 0; 
	double * tSum = new double[nS] ; 
		
/////////////////////////////////////////////////////////////////////////////////////

	if(m_psi == NULL)
		m_psi = Allocate2DMatrix(nLoci, nsk); 
	
	//at the rightmost locus M; 
	for (int j = 0; j < nsk; j++) (m_psi)[nLoci-1][j] = 1.0;

	//do the recursive calc backwards; 
	for (int m = nLoci - 1; m > 0; m--)
	{
		ts = 0; 
		for (int j = 0, s = 0; s < nS; s++)
		{
			tSum[s] = 0.0;
			for (int k = 0; k < nK; k++, j++)
				tSum[s] += prG(theta[m][k], GetsnpGT(m)) * m_psi[m][j] * beta[m][s][k]; 
			ts += tSum[s] * admix[s];
		}

		for (int j = 0, s = 0; s < nS; s++)
			for (int k = 0; k < nK; k++, j++)
			{
				double temp =(1.0 - rk[m]) * prG(theta[m][k], GetsnpGT(m)) * m_psi[m][j] + rk[m] * tSum[s]; 
				temp *= (1.0 - rs[m]);
				m_psi[m-1][j] = rs[m] * ts + temp; 
			}
		for (int j = 0; j < nsk; j++)
			m_psi[m-1][j] += 1e-20;

		ts = 0; 
		for (int j = 0; j < nsk; j++)
			ts += m_psi[m-1][j];
		for (int j= 0; j < nsk; j++)
				m_psi[m-1][j] /= ts;
	}
	
	///////////////////////////////////////////////////////////////////////////////

	if(m_phi == NULL)
		m_phi = Allocate2DMatrix(nLoci, nsk); 
	//at locus 0 -- leftmost; 
	for (int j = 0, s = 0; s < nS; s++)
		for (int k = 0; k < nK; k++, j++)
	   		m_phi[0][j] = admix[s] * beta[0][s][k] * prG(theta[0][k], GetsnpGT(0));
	
	double phiscale = 0; 
	for (int m = 0; m < nLoci-1; m++)
	{
		ts = 0; 
		for (int j = 0, s = 0; s < nS; s++)
		{
			tSum[s] = 0.0;
			for (int k = 0; k < nK; k++, j++)
				tSum[s] += m_phi[m][j];
			ts += tSum[s]; 
		}

		for (int j = 0, s = 0; s < nS; s++)
			for (int k = 0; k < nK; k++, j++)
			{   
				double temp = (1.0 - rk[m+1]) * m_phi[m][j] + rk[m+1] * beta[m+1][s][k] * tSum[s]; 
				temp *= (1.0 - rs[m+1]); 
				if(nS > 1) 
					temp += rs[m+1] * admix[s] * beta[m+1][s][k] * ts; 
				m_phi[m+1][j] = temp * prG(theta[m+1][k], GetsnpGT(m+1));
			}   
		for (int j = 0; j < nsk; j++)
		   	m_phi[m+1][j] += 1e-20;

		ts = 0; 
		for (int j = 0; j < nsk; j++)
		   	ts += m_phi[m+1][j];
		for (int j = 0; j < nsk; j++)
			m_phi[m+1][j] /= ts;
		phiscale += log(ts);  
	}                                                          

	m_loglike = phiscale / log(10.0);
	
	///////////////////////////////////////////////////////////////////////////////////
	double * prz2 = new double[nsk]; 
	double * prz = new double[nK]; 
	double * prs = new double[nS]; 
	double * tadmix = new double[nS]; 
	for (int s = 0; s < nS; s++)
		tadmix[s] = 0; 
	for (int m = 0; m < nLoci; m++)
	{
		for (int j = 0; j < nsk; j++)
			prz2[j] = m_phi[m][j] * m_psi[m][j]; 
		normalize(prz2, nsk); 

		for (int k = 0; k < nK; k++)
			prz[k] = 0; 
		for (int j = 0, s = 0; s < nS; s++)
			for (int k = 0; k < nK; k++, j++)
			{
				prz[k] += prz2[j];
				tadmix[s] += prz2[j];
			}

		if(sbot != NULL && stop != NULL) 
		{
			for (int k = 0; k < nK; k++)
				prz[k] *= get_weight(); 
			switch(GetsnpGT(m))
			{
				case 0: 
					for (int j = 0; j < nK; j++)
						sbot[m][j] += prz[j];
					break;
				case 1:
					for (int j = 0; j < nK; j++)
					{
						stop[m][j] += prz[j];
						sbot[m][j] += prz[j];
					}
					break;
				default:
					break;
			}
		}
	}

	////////////////////////////////////////////////////////////////////////
	if(sjmsk != NULL) 
	{
		for (int m = 1; m < nLoci; m++)
		{
			ts = 0; 
			for (int j = 0, s = 0; s < nS; s++)
			{
				tSum[s] = 0.0;
				for (int k = 0; k < nK; k++, j++)
					tSum[s] += m_phi[m-1][j];
				ts += tSum[s]; 
			}

			double tj = 0; 
			for (int j = 0, s = 0; s < nS; s++)
			{
				prs[s] = 0; 
				for (int k = 0; k < nK; k++, j++)
				{                                          
					double bj = prG(theta[m][k], GetsnpGT(m)) * m_psi[m][j];   
					double s0k0 = m_phi[m-1][j] * (1.0 - rs[m]) * (1.0 - rk[m]) * bj;        
					double s0k1 = tSum[s] * (1.0 - rs[m]) * rk[m] * beta[m][s][k] * bj;
					double s1 = ts * rs[m] * admix[s] * beta[m][s][k] * bj;        
					prz2[j] = s0k1; 
					prs[s] += s1; 
					tj += s0k0 + s0k1 + s1; 
				}
			}
			for (int j = 0; j < nsk; j++)
		   		sjmsk[m][j] += prz2[j] / tj * get_weight(); 
			for (int s = 0; s < nS; s++)
				sjms[m][s] += prs[s] / tj * get_weight(); 
		}

		for (int j = 0; j < nsk; j++)
	   		prz2[j] = m_phi[0][j] * m_psi[0][j];
		normalize(prz2, nsk); 
		for (int j = 0; j < nsk; j++)
	   		prz2[j] *= get_weight();

		for (int j = 0, s = 0; s < nS; s++)
			for (int k = 0; k < nK; k++, j++)
			{
				sjmsk[0][j] += prz2[j]; 
				sjms[0][s] += prz2[j]; 
			}
	}

	normalize(tadmix, nS); 
	for (int s = 0; s < nS; s++)
		admix[s] = tadmix[s]; 

	delete[] tSum; 
	delete[] prz2; 
	delete[] prz; 
	delete[] prs; 
	delete[] tadmix; 
}

void HapInd::forward_backward_block(const int nS, const int nK, ModelParam * pMP, double ** stop, double ** sbot, double ** sjmsk, double ** sjms, int which, int b, vector<int> pbreak) 
{
	//[pbreak.at(0), pbreak.at(1)-1]
	//[pbreak.at(1), pbreak.at(2)-1]
	//[pbreak.at(2), pbreak.at(3)-1]
	//
	if(m_admix == NULL) {
		m_admix = new double[nS]; 
		for (int s = 0; s < nS; s++)
			m_admix[s] = 1.0 / nS;
	}
	double * admix = m_admix; 

	double * rs = pMP->Getrs(); 
	double * rk = pMP->Getrk(); 
	double ** theta = pMP->Gettheta(); 
	double *** beta = pMP->Getbeta(); 

	int nsk = nS * nK; 
	double ts = 0; 
	double * prz2 = new double[nsk]; 
	double * prz = new double[nK]; 
	double * prGk = new double[nK]; 
	double * tSum = new double[nS]; 

	/////////////////////////////////////////////////////////////////////////////////////
	if(b == 0) //first time around, compute backward all the way and record at the end points of each interval; 
	{
		double * psi = new double[nsk]; 

		vector<vector<double> > ().swap(vvfwdp); 
		vector<vector<double> > ().swap(vvfwd); 
		vector<vector<double> > ().swap(vvbwd); 
		vector<double > tv1; 
		for (int j = 0; j < nsk; j++)
			tv1.push_back(0); 
		for (unsigned i = 0; i < pbreak.size()-1; i++)
		{
			vvfwdp.push_back(tv1); 
			vvfwd.push_back(tv1); 
			vvbwd.push_back(tv1); 
		}   //initialize vvfwd and vvbwd; 
		vector<double > ().swap(tv1); 

		for (int j = 0; j < nsk; j++) 
			psi[j] = 1.0; 
		for (int j = 0; j < nsk; j++)
			vvbwd.at(pbreak.size()-2).at(j) = psi[j]; 
		for (unsigned f = pbreak.size()-2; f > 0; f--)
		{
			int len = pbreak.at(f+1) - pbreak.at(f); 

			for (int p = len-1; p >= 0; p--)
			{
				int m = pbreak.at(f) + p; 
				for (int k = 0; k < nK; k++)
					prGk[k] = prG(theta[m][k], GetsnpGT(m));

				ts = 0; 
				for (int j = 0, s = 0; s < nS; s++)
				{
					tSum[s] = 0.0;
					for (int k = 0; k < nK; k++, j++)
						tSum[s] +=  prGk[k] * psi[j] * beta[m][s][k]; 
					ts += tSum[s] * admix[s]; 
				}

				for (int j = 0, s = 0; s < nS; s++)
					for (int k = 0; k < nK; k++, j++)
					{
						double temp = (1.0 - rk[m]) * prGk[k] * psi[j] + rk[m] * tSum[s]; 
						temp *= (1.0 - rs[m]);
						psi[j] = temp + rs[m] * ts; 
					}
				for (int j = 0; j < nsk; j++)
					psi[j] += 1e-20; 
				ts = 0; 
				for (int j = 0; j < nsk; j++)
					ts += psi[j]; 

				normalize(psi, nsk); 
			}
			for (int j = 0; j < nsk; j++)
				vvbwd.at(f-1).at(j) = psi[j]; 
		}
		delete[] psi; 
	}
	//hello

	int len = pbreak.at(b+1) - pbreak.at(b); 
	if(m_psi == NULL)
		m_psi = Allocate2DMatrix(len, nsk); 

	for (int j = 0; j < nsk; j++)
		m_psi[len-1][j] = vvbwd.at(b).at(j); 
	//do the recursive calc backwards; 
	for (int p = len-1; p > 0; p--)
	{
		int m = pbreak.at(b) + p; 
		for (int k = 0; k < nK; k++)
			prGk[k] = prG(theta[m][k], GetsnpGT(m));

		ts = 0; 
		for (int j = 0, s = 0; s < nS; s++)
		{
			tSum[s] = 0.0;
			for (int k = 0; k < nK; k++, j++)
				tSum[s] +=  prGk[k] * m_psi[p][j] * beta[m][s][k]; 
			ts += tSum[s] * admix[s]; 
		}

		for (int j = 0, s = 0; s < nS; s++)
			for (int k = 0; k < nK; k++, j++)
			{
				double temp = (1.0 - rk[m]) * prGk[k] * m_psi[p][j] + rk[m] * tSum[s]; 
				temp *= (1.0 - rs[m]);
				m_psi[p-1][j] = temp + rs[m] * ts; 
			}
		for (int j = 0; j < nsk; j++)
			m_psi[p-1][j] += 1e-20; 
		ts = 0; 
		for (int j = 0; j < nsk; j++)
			ts += m_psi[p-1][j]; 
		normalize(m_psi[p-1], nsk); 
	}
	
	///////////////////////////////////////////////////////////////////////////////
	if(m_phi == NULL)
		m_phi = Allocate2DMatrix(len, nsk); 

	//at locus 0 -- leftmost; 

	if(b == 0) {
		for (int k = 0; k < nK; k++)
			prGk[k] = prG(theta[0][k], GetsnpGT(0));
		for (int j = 0, s = 0; s < nS; s++)
			for (int k = 0; k < nK; k++, j++)
				m_phi[0][j] = admix[s] * beta[0][s][k] * prGk[k];
		for (int j = 0; j < nsk; j++)
			vvfwdp.at(0).at(j) = 1; 
		for (int j = 0; j < nsk; j++)
			vvfwd.at(0).at(j) = m_phi[0][j]; 
	}
	else {
		for (int j = 0; j < nsk; j++)
			m_phi[0][j] = vvfwd.at(b).at(j); 
	}
	
	ts = 0; 
	for (int j = 0; j < nsk; j++)
		ts += m_phi[0][j]; 
	double scale = log(ts); 
	for (int p = 0; p < len-1; p++)
	{
		int m = pbreak.at(b) + p; 
		for (int k = 0; k < nK; k++)             	
			prGk[k] = prG(theta[m+1][k], GetsnpGT(m+1)); 

		ts = 0; 
		for (int j = 0, s = 0; s < nS; s++)
		{
			tSum[s] = 0.0;
			for (int k = 0; k < nK; k++, j++)
				tSum[s] += m_phi[p][j];
			ts += tSum[s]; 
		}

		for (int j = 0, s = 0; s < nS; s++)
			for (int k = 0; k < nK; k++, j++)
			{   
				double temp = (1.0 - rk[m+1]) * m_phi[p][j] + rk[m+1] * beta[m+1][s][k] * tSum[s]; 
				temp *= (1.0 - rs[m+1]); 
				temp += rs[m+1] * (admix[s] * beta[m+1][s][k]) * ts; 
				m_phi[p+1][j] = temp * prGk[k]; 
			}   
		for (int j = 0; j < nsk; j++)
			m_phi[p+1][j] += 1e-20; 
		double ts = 0; 
		for (int j = 0; j < nsk; j++)
			ts += m_phi[p+1][j]; 
		for (int j = 0; j < nsk; j++)
			m_phi[p+1][j] /= ts; 
		scale += log(ts); 
	}                                                          
	m_loglike += scale / log(10.0);
	//do one more step to obtain the starting fwd of the next block; 
	if(b < (int) pbreak.size()-2)   //not last block
	{
		for (int j = 0; j < nsk; j++)
			vvfwdp.at(b+1).at(j) = m_phi[len-1][j]; 
		//this is needed to compute expected jumps at the block ends; 
		double * phi1 = new double[nsk]; 
		int p = len-1; 
		int m = pbreak.at(b+1)-1; 
		for (int k = 0; k < nK; k++)
			prGk[k] = prG(theta[m+1][k], GetsnpGT(m+1)); 

		ts = 0; 
		for (int j = 0, s = 0; s < nS; s++)
		{
			tSum[s] = 0.0;
			for (int k = 0; k < nK; k++, j++)
				tSum[s] += m_phi[p][j];
			ts += tSum[s]; 
		}

		for (int j = 0, s = 0; s < nS; s++)
			for (int k = 0; k < nK; k++, j++)
			{   
				double temp = (1.0 - rk[m+1]) * m_phi[p][j] + rk[m+1] * beta[m+1][s][k] * tSum[s]; 
				temp *= (1.0 - rs[m+1]); 
				temp += rs[m+1] * (admix[s] * beta[m+1][s][k]) * ts; 
				phi1[j] = temp * prGk[k]; 
			}   
		normalize(phi1, nsk);
		for (int j = 0; j < nsk; j++)
			vvfwd.at(b+1).at(j) = phi1[j]; 
		delete[] phi1; 
	}                                                          

//	{
//		cout << " b = " << b << endl; 
//		for (int j = 0; j < nsk; j++)
//			cout << vvfwdp.at(b).at(j) << " "; 
//		cout << endl; 
//		for (int j = 0; j < nsk; j++)
//			cout << m_phi[len-1][j] << " "; 
//		cout << endl; 
//	}

	if(sjmsk != NULL) 
	{                         
		for (int p = 0; p < len; p++)
		{                                       
			int m = pbreak.at(b) + p; 
			for (int k = 0; k < nK; k++)
				prGk[k] = prG(theta[m][k], GetsnpGT(m)); 

			ts = 0; 
			for (int j = 0, s = 0; s < nS; s++)
			{
				tSum[s] = 0.0;
				if(p == 0) {
					for (int k = 0; k < nK; k++, j++)
						tSum[s] += vvfwdp.at(b).at(j); 
				}
				else {
					for (int k = 0; k < nK; k++, j++)
						tSum[s] += m_phi[p-1][j];
				}
				ts += tSum[s]; 
			}

			double tj = 0; 
			for (int j = 0, s = 0; s < nS; s++)
			{
				prz[s] = 0; 
				for (int k = 0; k < nK; k++, j++)
				{                    
					double bj = prGk[k] * m_psi[p][j]; 
					double s0k0 = 0; 
					if(p == 0) 
						s0k0 = vvfwdp.at(b).at(j) * (1.0 - rs[m]) * (1.0 - rk[m]) * bj; 
					else 
						s0k0 = m_phi[p-1][j] * (1.0 - rs[m]) * (1.0 - rk[m]) * bj; 
					double s0k1 = tSum[s] * (1.0 - rs[m]) * rk[m] * beta[m][s][k] * bj;   
					double s1 = ts * rs[m] * admix[s] * beta[m][s][k] * bj;  
					prz2[j] = s0k1; 
					prz[s] += s1; 
					tj += s0k0 + s0k1 + s1; 
				}   
			}

			if(m > 0) {
				for (int j = 0; j < nsk; j++)
					sjmsk[m][j] += prz2[j] / tj * get_weight(); 
				for (int j = 0; j < nS; j++)
					sjms[m][j] += prz[j] / tj * get_weight(); 
			}
			else {
				for (int j = 0; j < nsk; j++)
					prz2[j] = m_phi[0][j] * m_psi[0][j];
				normalize(prz2, nsk); 
				for (int j = 0; j < nsk; j++)
					sjmsk[0][j] += prz2[j] * get_weight(); 
			}
		//the first locus, it's important for beta. 
		}
	}

	////////////////////////////////////////////////////////////////////////////
	if(stop != NULL && sbot != NULL) 
	{
		for (int p = 0; p < len; p++)
		{                                       
			int m = pbreak.at(b) + p; 
			for (int j = 0; j < nsk; j++)
				prz2[j] = (m_phi[p][j]) * (m_psi[p][j]); 
			normalize(prz2, nsk); 

			for (int k = 0; k < nK; k++)
				prz[k]= 0; 
			for (int j = 0, s = 0; s < nS; s++)
				for (int k = 0; k < nK; k++, j++)
					prz[k] += prz2[j];

			for (int k = 0; k < nK; k++)
				prz[k] *= get_weight(); 
			int gt = GetsnpGT(m); 
			if(gt == 0) 
			{
				for (int k = 0; k < nK; k++)
					sbot[m][k] += prz[k];
			}
			else //if (gt == 1)
			{
				for (int k = 0; k < nK; k++)
				{
					stop[m][k] += prz[k];      
					sbot[m][k] += prz[k]; 
				}
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////
	delete[] prz; 
	delete[] prz2; 
	delete[] tSum; 
	delete[] prGk; 
}
