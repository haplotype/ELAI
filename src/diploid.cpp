#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "indiv.h"
#include "fpmath.h"
#include "diploid.h"
#include "haploid.h"
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

DipInd::DipInd(void)
{
	m_phi = NULL;
	m_psi = NULL;

	snp_dstr = NULL; 
	m_mgt = NULL; 
	m_ploid = 2; 
	m_ps = NULL; 
	m_pk = NULL; 
	m_beta = NULL; 
	m_maf = NULL; 
	m_admix = NULL; 
	m_admix1 = NULL; 
}

DipInd::DipInd(int n)
{
	nLoci = n; 
	m_phi = NULL;
	m_psi = NULL;

	snp_dstr = NULL; 
	m_mgt = NULL; 
	m_ploid = 2; 
	m_ps = NULL; 
	m_pk = NULL; 
	m_beta = NULL; 
	m_maf = NULL; 
	m_admix = NULL; 
	m_admix1 = NULL; 
}

DipInd::~DipInd(void)
{
 	Free2DMatrix(m_phi);        m_phi = NULL;
	Free2DMatrix(m_psi);      	m_psi = NULL;
} 

double DipInd::prG1(double h, double t2, short gt)
{
    double res = 1.0; 
	if(gt == 0) 
		res = (1.0 - t2);
	else if (gt == 2) 
		res = t2;
	else if (gt == 1)
		res = 0.5;
	return res; 
}

double DipInd::prG(double t1, double t2, short gt)
{
    double res = 1.0; 
	if(gt == 0) 
		res = (1.0 - t2);
//		res = (1.0 - t1) * (1.0 - t2);
	else if (gt == 2) 
		res = t2;
//		res = t1 * t2;
	else if (gt == 1)
		res = t1 * (1-t2) + (1-t1) * t2;
	return res; 
}


double DipInd::prG2(double t1, double t2, short gt)
{
    double res = 1.0; 
	if(gt == 0)
		res = (1.0 - t1) * (1.0 - t2);
	else if (gt == 2) 
		res = t1 * t2;
	else if (gt == 1) 
		res = t1 * (1 - t2) + (1 - t1) * t2; 
	return (res);  
}

void DipInd::clean_after_em(int deep)
{
	if(m_phi){ Free2DMatrix(m_phi); m_phi = NULL;} 
	if(m_psi) {Free2DMatrix(m_psi); m_psi = NULL;} 
	if(phiScale) {delete[] phiScale; phiScale = NULL;} 
	if(psiScale) {delete[] psiScale; psiScale = NULL;} 
	if(m_pjk) {delete[] m_pjk; m_pjk = NULL;} 
	if(m_pjs) {delete[] m_pjs; m_pjs = NULL;} 
	if(deep == 1) 
	{
		Free2DMatrix(m_ps); m_ps = NULL; 
		Free2DMatrix(m_pk); m_pk = NULL; 
		if(m_maf) {delete[] m_maf; m_maf = NULL;}
	}
}

void DipInd::compute_ps(const int nS, const int nK, ModelParam * pMP, const int comp_mode)
{

//	cout << "in dip:compute_ps " << comp_mode << endl; 
//	CalcAll(nS, nK, pMP, NULL, NULL, NULL, NULL, comp_mode, 11); 
	if(comp_mode == 1) 
	{
		if(m_maf == NULL) 
		{
			m_maf = new double[nLoci]; 
			for (int m = 0; m < nLoci; m++)
				m_maf[m] = 0.5; 
		}
//		forward_backward_haploid(nS, nK, pMP, NULL, NULL, NULL, NULL, 0); 
//		for (int j = 0; j < 5; j++)
//		{
//			forward_backward_haploid(nS, nK, pMP, NULL, NULL, NULL, NULL, 1); 
//			forward_backward_haploid(nS, nK, pMP, NULL, NULL, NULL, NULL, 2); 
//		}
//		if(m_ps == NULL) m_ps = Allocate2DMatrix(nLoci, nS); 
//		for (int m = 0; m < nLoci; m++)
//			for (int s = 0; s < nS; s++)
//				m_ps[m][s] = 0; 
//		forward_backward_haploid(nS, nK, pMP, NULL, NULL, NULL, NULL, 1); 
//		forward_backward_haploid(nS, nK, pMP, NULL, NULL, NULL, NULL, 2); 
//		clean_after_em(0); 

		int unit = 2000; 
		int nfile = ceil(nLoci / (double) unit); 
		vector<int> pbreak; 
		for (int i = 0; i < nfile; i++)
			pbreak.push_back(i * unit); 
		pbreak.push_back(nLoci); 
//		for (unsigned b = 0; b < pbreak.size() - 1; b++) 
//			forward_backward_block(nS, nK, pMP, NULL, NULL, NULL, NULL, 0, b, pbreak); 
//		for (int j = 0; j < 3; j++) {
//			for (unsigned b = 0; b < pbreak.size() - 1; b++) 
//				forward_backward_block(nS, nK, pMP, NULL, NULL, NULL, NULL, 1, b, pbreak); 
//			for (unsigned b = 0; b < pbreak.size() - 1; b++) 
//				forward_backward_block(nS, nK, pMP, NULL, NULL, NULL, NULL, 2, b, pbreak); 
//		}
		for (unsigned b = 0; b < pbreak.size() - 1; b++) 
			forward_backward_block1(nS, nK, pMP, 0, b, pbreak); 
		if(m_ps == NULL) m_ps = Allocate2DMatrix(nLoci, nS); 
		for (int m = 0; m < nLoci; m++)
			for (int s = 0; s < nS; s++)
				m_ps[m][s] = 0; 
		for (unsigned b = 0; b < pbreak.size() - 1; b++) 
			forward_backward_block1(nS, nK, pMP, 1, b, pbreak); 
		for (unsigned b = 0; b < pbreak.size() - 1; b++) 
			forward_backward_block1(nS, nK, pMP, 2, b, pbreak); 
		clean_after_em(0); 
	}
	else {
		if(m_ps == NULL) m_ps = Allocate2DMatrix(nLoci, nS*nS); 
		for (int m = 0; m < nLoci; m++)
			for (int s = 0; s < nS*nS; s++)
				m_ps[m][s] = 0; 
//		CalcAll(nS, nK, pMP, NULL, NULL, NULL, NULL, comp_mode, 11); 
        forward_backward_diploid(nS, nK, pMP, NULL, NULL, NULL, NULL, 11); 
	}
}

void DipInd::compute_ps(const int nS, const int nK, ModelParam * pMP, int comp_mode, int b, vector<int> pbreak)
{
	int len = pbreak.at(b+1) - pbreak.at(b); 
	if(m_ps != NULL) {
		Free2DMatrix(m_ps); m_ps = NULL; 
	}
	m_ps = Allocate2DMatrix(len, nS*nS); 
	
	for (int m = 0; m < len; m++)
		for (int s = 0; s < nS*nS; s++)
			m_ps[m][s] = 0; 
   	forward_backward_block2(nS, nK, pMP, 2, b, pbreak); 
	clean_after_em(0); 
}

void DipInd::compute_pk(const int nS, const int nK, ModelParam * pMP, const int comp_mode)
{
	if(m_pk == NULL) m_pk = Allocate2DMatrix(nLoci, nK); 
	for (int m = 0; m < nLoci; m++)
		for (int k = 0; k < nK; k++)
			m_pk[m][k] = 0; 
	CalcAll(nS, nK, pMP, NULL, NULL, NULL, NULL, comp_mode, 11); 
	clean_after_em(0); 
}


void DipInd::compute_pk(const int nS, const int nK, ModelParam * pMP, int comp_mode, int b, vector<int> pbreak)
{
	int len = pbreak.at(b+1) - pbreak.at(b); 
	if(m_pk != NULL) {
		Free2DMatrix(m_pk); m_pk = NULL; 
	}
	m_pk = Allocate2DMatrix(len, nK); 

	for (int m = 0; m < len; m++)
		for (int k = 0; k < nK; k++)
			m_pk[m][k] = 0; 
	if(comp_mode == 2) {
		forward_backward_block2(nS, nK, pMP, 2, b, pbreak); 
		clean_after_em(0); 
	} else {
		if(m_maf == NULL) 
		{
			m_maf = new double[nLoci]; 
			for (int m = 0; m < nLoci; m++)
				m_maf[m] = 0.5; 
		}
		forward_backward_block1(nS, nK, pMP, 0, b, pbreak); 
		forward_backward_block1(nS, nK, pMP, 1, b, pbreak); 
		forward_backward_block1(nS, nK, pMP, 2, b, pbreak); 
		clean_after_em(0); 
	}

}

void DipInd::forward_backward_diploid(const int nS, const int nK, ModelParam * pMP, double ** stop, double ** sbot, double ** sjmsk, double ** sjms, int which) 
{
	if(m_admix == NULL) 
	{
		m_admix = new double[nS]; 
		for (int s = 0; s < nS; s++)
			m_admix[s] = 1.0 / nS; 
	}
	
	double * admix = m_admix; 

	double * rk = pMP->Getrk();
	double * rs = pMP->Getrs();
	double ** theta = pMP->Gettheta(); 
	double *** beta = pMP->Getbeta(); 
	
	int nsk = nS * nK; 
	double ** prz2 = Allocate2DMatrix(nsk, nsk); 
	double ** prz = Allocate2DMatrix(nK, nK); 
	double ** prGk = Allocate2DMatrix(nK, nK); 
	double * tSums = new double[nS];  
	double ts = 0;

	///////////////////////////////////////////////////////////////////////////////
	if(m_psi == NULL)
		m_psi = Allocate2DMatrix(nLoci, nsk*nsk); 

	for (int j1 = 0; j1 < nsk; j1++)             
	for (int j2 = 0; j2 < nsk; j2++)
	   	m_psi[nLoci-1][j1*nsk+j2] = 1.0;

	for (int m = nLoci - 1; m > 0; m--)
	{
		for (int k1 = 0; k1 < nK; k1++)
		for (int k2 = 0; k2 < nK; k2++)
		   	prGk[k1][k2] = prG2(theta[m][k1],theta[m][k2], GetsnpGT(m)); 

		for (int j1 = 0, s1 = 0; s1 < nS; s1++)
		for (int k1 = 0; k1 < nK; k1++, j1++)
		{
			double ts = 0; 
			for (int j2 = 0, s2 = 0; s2 < nS; s2++)
			{
				tSums[s2] = 0; 
				for (int k2 = 0; k2 < nK; k2++, j2++)
					tSums[s2] += prGk[k1][k2] * m_psi[m][j1*nsk+j2] * beta[m][s2][k2]; 
				ts += tSums[s2] * admix[s2];
			}
   			for (int j2 = 0, s2 = 0; s2 < nS; s2++)
   			for (int k2 = 0; k2 < nK; k2++, j2++)
   			{ 
				double t2 = (1.0 - rk[m]) * prGk[k1][k2] * m_psi[m][j1*nsk+j2] + rk[m] * tSums[s2]; 
				t2 *= (1.0 - rs[m]); 
				t2 += rs[m] * ts;
				prz2[j1][j2] = t2;
			}
		}   

		for (int j2 = 0, s2 = 0; s2 < nS; s2++)
		for (int k2 = 0; k2 < nK; k2++, j2++)
		{
			double ts = 0; 
			for (int j1 = 0, s1 = 0; s1 < nS; s1++)
			{	
				tSums[s1] = 0; 
				for (int k1 = 0; k1 < nK; k1++, j1++)
					tSums[s1] += prz2[j1][j2] * beta[m][s1][k1]; 
				ts += tSums[s1] * admix[s1]; 
			}
			for (int j1 = 0, s1 = 0; s1 < nS; s1++)
			for (int k1 = 0; k1 < nK; k1++, j1++)  
			{
				double t1 = (1.0 - rk[m]) * prz2[j1][j2] + rk[m] * tSums[s1];
				t1 *= (1.0 - rs[m]); 
				t1 += rs[m] * ts; 
				m_psi[m-1][j1*nsk+j2] = t1; 
			}
		}

		ts = 0;
		for (int j1 = 0; j1 < nsk; j1++)
			for (int j2 = 0; j2 < nsk; j2++)
				ts += m_psi[m-1][j1*nsk+j2];
		for (int j1 = 0; j1 < nsk; j1++)
			for (int j2 = 0; j2 < nsk; j2++)
				m_psi[m-1][j1*nsk+j2] /= ts;
	}                         
	///////////////////////////////////////////////////////////////////////////////
	
	double phiscale = 0; 
	if(m_phi == NULL)
		m_phi = Allocate2DMatrix(nLoci, nsk*nsk); 

	//at locus 0 -- leftmost; 
	for (int k1 = 0; k1 < nK; k1++)
	for (int k2 = 0; k2 < nK; k2++)
	   	prGk[k1][k2] = prG2(theta[0][k1],theta[0][k2], GetsnpGT(0)); 
	for (int j1 = 0, s1 = 0; s1 < nS; s1++)
	for (int k1 = 0; k1 < nK; k1++, j1++)  
		for (int j2 = 0, s2 = 0; s2 < nS; s2++)
		for (int k2 = 0; k2 < nK; k2++, j2++)
			m_phi[0][j1*nsk+j2] = admix[s1] * admix[s2] * beta[0][s1][k1] * beta[0][s2][k2] * prGk[k1][k2];
		
	for (int m = 0; m < nLoci-1; m++)
	{
		for (int k1 = 0; k1 < nK; k1++)
		for (int k2 = 0; k2 < nK; k2++)
			prGk[k1][k2] = prG2(theta[m+1][k1],theta[m+1][k2], GetsnpGT(m+1)); 

		for (int j1 = 0, s1 = 0; s1 < nS; s1++)
		for (int k1 = 0; k1 < nK; k1++, j1++)
		{
			double ts = 0; 
			for (int j2 = 0, s2 = 0; s2 < nS; s2++)
			{
				tSums[s2] = 0; 
				for (int k2 = 0; k2 < nK; k2++, j2++)
					tSums[s2] += m_phi[m][j1*nsk+j2]; 
				ts += tSums[s2];
			}
			for (int j2 = 0, s2 = 0; s2 < nS; s2++)
			for (int k2 = 0; k2 < nK; k2++, j2++)
			{ 
				double t2 = (1.0 - rk[m+1]) * m_phi[m][j1*nsk+j2] + rk[m+1] * beta[m+1][s2][k2] * tSums[s2];
				t2 *= (1.0 - rs[m+1]); 
                t2 += rs[m+1] * ts * admix[s2] * beta[m+1][s2][k2]; 
				prz2[j1][j2] = t2; 
			} 
		}

		for (int j2 = 0, s2 = 0; s2 < nS; s2++)
		for (int k2 = 0; k2 < nK; k2++, j2++)
		{
			double ts = 0; 
			for (int j1 = 0, s1 = 0; s1 < nS; s1++)
			{	
				tSums[s1] = 0; 
				for (int k1 = 0; k1 < nK; k1++, j1++)
					tSums[s1] += prz2[j1][j2]; 
				ts += tSums[s1]; 
			}
			for (int j1 = 0, s1 = 0; s1 < nS; s1++)
			for (int k1 = 0; k1 < nK; k1++, j1++)  
			{
				double t1 = (1.0 - rk[m+1]) * prz2[j1][j2] + rk[m+1] * beta[m+1][s1][k1] * tSums[s1];
				t1 *= (1.0 - rs[m+1]); 
				t1 += rs[m+1] * ts * admix[s1] * beta[m+1][s1][k1]; 
				m_phi[m+1][j1*nsk+j2] = t1 * prGk[k1][k2]; 
			}
		}

		ts = 0; 
		for (int j1 = 0; j1 < nsk; j1++)
			for (int j2 = 0; j2 < nsk; j2++)
				ts += m_phi[m+1][j1*nsk+j2];
		for (int j1 = 0; j1 < nsk; j1++)
			for (int j2 = 0; j2 < nsk ; j2++)
				m_phi[m+1][j1*nsk+j2] /= ts;
		phiscale += log(ts);  
	}                                                          

	m_loglike = phiscale / log(10.0);

	///////////////////////////////////////////////////////////////////////////////////

	double * tadmix = new double[nS]; 
	for (int s = 0; s < nS; s++)
		tadmix[s] = 0; 
	double * dinc = new double[nK]; 
	for (int m = 0; m < nLoci; m++)
	{
		for (int j1 = 0; j1 < nsk; j1++)
		for (int j2 = 0; j2 < nsk; j2++)
		   	prz2[j1][j2] = m_phi[m][j1*nsk+j2] * m_psi[m][j1*nsk+j2]; 
		normalize(&prz2[0][0], nsk * nsk); 

		for (int k1 = 0; k1 < nK; k1++)
			for (int k2 = 0; k2 < nK; k2++)
				prz[k1][k2] = 0; 
		for (int j1 = 0, s1 = 0; s1 < nS; s1++)
			for (int k1 = 0; k1 < nK; k1++, j1++)
				for (int j2 = 0, s2 = 0; s2 < nS; s2++)
					for (int k2 = 0; k2 < nK; k2++, j2++)
					{
						prz[k1][k2] += prz2[j1][j2]; 
						tadmix[s1] += prz2[j1][j2]; 
					}

		if(m_ps != NULL) 
		{
			if(which == 11) {
				for (int j1 = 0, s1 = 0; s1 < nS; s1++)
					for (int k1 = 0; k1 < nK; k1++, j1++)
						for (int j2 = 0, s2 = 0; s2 < nS; s2++)
							for (int k2 = 0; k2 < nK; k2++, j2++)
								m_ps[m][s1*nS+s2] += prz2[j1][j2]; 
			}
			else {
				for (int j1 = 0, s1 = 0; s1 < nS; s1++)
					for (int k1 = 0; k1 < nK; k1++, j1++)
						for (int j2 = 0, s2 = 0; s2 < nS; s2++)
							for (int k2 = 0; k2 < nK; k2++, j2++)
								m_ps[m][s1] += 2*prz2[j1][j2]; 
			}
		}
		if(m_pk != NULL) 
		{
			for (int j1 = 0, s1 = 0; s1 < nS; s1++)
				for (int k1 = 0; k1 < nK; k1++, j1++)
					for (int j2 = 0, s2 = 0; s2 < nS; s2++)
						for (int k2 = 0; k2 < nK; k2++, j2++)
							m_pk[m][k1] += 2*prz2[j1][j2]; 
		}
                                      
		if(stop != NULL && sbot != NULL) 
		{           
			short gt = GetsnpGT(m);
			for (int j1 = 0; j1 < nK; j1++)
			{
				dinc[j1] = 0; 
				for (int j2 = 0; j2 < nK; j2++)
					dinc[j1] += (j2 == j1 ? (2.0 * prz[j1][j2]) : prz[j1][j2]); 
				dinc[j1] *= get_weight(); 
			}
			if(gt == 0) 
			{
				for (int j1 = 0; j1 < nK; j1++)
					sbot[m][j1] += dinc[j1];
			}
			else if(gt == 2) 
			{
				for (int j1 = 0; j1 < nK; j1++)
				{
					stop[m][j1] += dinc[j1]; 
					sbot[m][j1] += dinc[j1];
				}
			}
			else if(gt == 1) 
			{
				for (int j1 = 0; j1 < nK; j1++)
				{
					sbot[m][j1] += dinc[j1]; 
					for (int j2 = 0; j2 < nK; j2++)
					{
						if(j2 == j1) 
							stop[m][j1] += prz[j1][j2] * weight;    
						else 
						{
							double t1 = theta[m][j1] * (1. - theta[m][j2]);
							double t2 = theta[m][j2] * (1. - theta[m][j1]);					
							stop[m][j1] += prz[j1][j2] * t1 / (t1 + t2) * get_weight();    
						}
					}
				}
			}
		}
	}
	delete[] dinc; 

	//////////////////////////////////////////////////////////////////////////////////
//	if(sjmsk != NULL)
//	{
//		double * jump2j = new double[nsk]; 
//		double * jump2s = new double[nsk]; 
//		double * prf1 = new double[nsk];
//		double * prb1 = new double[nsk];
//		for (int m = 1; m < nLoci; m++)
//		{
//			for (int j = 0; j < nsk; j++)
//				jump2j[j] = jump2s[j] = prf1[j] = prb1[j] = 0; 
//
//			for (int k1 = 0; k1 < nK; k1++)
//			for (int k2 = 0; k2 < nK; k2++)
//				prGk[k1][k2] = prG2(theta[m][k1],theta[m][k2], GetsnpGT(m)); 
//
//			for (int j1 = 0, s1 = 0; s1 < nS; s1++)
//			for (int k1 = 0; k1 < nK; k1++, j1++)
//			{
//				double ts = 0; 
//				for (int j2 = 0, s2 = 0; s2 < nS; s2++)
//				{
//					tSums[s2] = 0; 
//					for (int k2 = 0; k2 < nK; k2++, j2++)
//						tSums[s2] += m_phi[m-1][j1*nsk+j2]; 
//					ts += tSums[s2];
//				}
//
//				for (int j2 = 0, s2 = 0; s2 < nS; s2++)
//				for (int k2 = 0; k2 < nK; k2++, j2++)
//				{ 
////					double t2 = (1.0 - rk[m]) * m_phi[m-1][j1*nsk+j2] + rk[m] * beta[m][s2][k2] * tSums[s2];
////					t2 *= (1.0 - rs[m]); 
////					t2 += rs[m] * ts * admix[s2] * beta[m][s2][k2]; 
////					prz2[j1][j2] = t2; 
//					double jk0 = (1.0 - rs[m]) * (1.0 - rk[m]) * m_phi[m-1][j1*nsk+j2]; 
//					double jk1 = (1.0 - rs[m]) * rk[m] * beta[m][s2][k2] * tSums[s2];
//					double js1 = rs[m] * ts * admix[s2] * beta[m][s2][k2]; 
//					prf1[j2] += jk0 + jk1 + js1; 
//					jump2j[j2] += jk1; 
//					jump2s[j2] += js1; 
//				} 
//			}
//
//			for (int j2 = 0, s2 = 0; s2 < nS; s2++)
//			for (int k2 = 0; k2 < nK; k2++, j2++)
//			{
//				for (int j1 = 0, s1 = 0; s1 < nS; s1++)
//				for (int k1 = 0; k1 < nK; k1++, j1++)  
//				{                                
//					double bj = prGk[k1][k2] * m_psi[m][j1*nsk+j2]; 
//					prb1[j2] += bj; 
//				}
//			}
//
//			double tj = 0; 
//			for (int j = 0, s = 0; s < nS; s++)
//			{
//				for (int k = 0; k < nK; k++, j++)
//				{
//					tj += prf1[j] * prb1[j]; 
//					jump2j[j] *= prb1[j]; 
//					jump2s[j] *= prb1[j]; 
//				}
//			}
//
//			if(tj < 1e-100) continue; 
//			for (int j = 0; j < nsk; j++)
//  		    	sjmsk[m][j] += 2*jump2j[j] / tj; 
//			for (int j = 0, s = 0; s < nS; s++)
//			{
//				double ts = 0; 
//				for (int k = 0; k < nK; k++, j++)
//					ts += jump2s[j]; 
//				sjms[m][s] += 2*ts / tj; 
//			}
//		}
//		delete[] jump2j;
//		delete[] jump2s;
//		delete[] prf1;
//		delete[] prb1;
//
//		for (int j1 = 0; j1 < nsk; j1++)
//			for (int j2 = 0; j2 < nsk; j2++)
//				prz2[j1][j2] = m_phi[0][j1*nsk+j2] * m_psi[0][j1*nsk+j2]; 
//		normalize(&prz2[0][0], nsk * nsk); 
//
//		for (int j1 = 0, s1 = 0; s1 < nS; s1++)
//			for (int k1 = 0; k1 < nK; k1++, j1++)
//				for (int j2 = 0, s2 = 0; s2 < nS; s2++)
//					for (int k2 = 0; k2 < nK; k2++, j2++)
//					{
//						sjmsk[0][j1] += prz2[j1][j2]; 
//						sjms[0][s1] += prz2[j1][j2]; 
//					}
//	}

	if(sjmsk != NULL)
	{

		double ** tmarginal = Allocate2DMatrix(nS, nK); 
		double *** tSumk = Allocate3DMatrix(nS, nK, nS); 
		double ** tDoubleSum = Allocate2DMatrix(nS, nS); 
		double * jump2j = new double[nsk]; 
		double * jump2s = new double[nS]; 
		for (int m = 1; m < nLoci; m++)
		{
			for (int k1 = 0; k1 < nK; k1++)
			for (int k2 = 0; k2 < nK; k2++)
				prGk[k1][k2] = prG2(theta[m][k1],theta[m][k2], GetsnpGT(m)); 

			for (int j1 = 0, s1 = 0; s1 < nS; s1++)
				for (int k1 = 0; k1 < nK; k1++, j1++)
				{
					tmarginal[s1][k1] = 0; 
					for (int j2 = 0, s2 = 0; s2 < nS; s2++)
					{
						tSumk[s1][k1][s2] = 0; 
						for (int k2 = 0; k2 < nK; k2++, j2++)
							tSumk[s1][k1][s2] += m_phi[m-1][j1*nsk+j2]; 
						tmarginal[s1][k1] += tSumk[s1][k1][s2];
					}
				}   //marginal sum; //

			double ts = 0;
			for (int s1 = 0; s1 < nS; s1++)
			{
				tSums[s1] = 0.0; 
				for (int s2 = 0; s2 < nS; s2++)
				{
					tDoubleSum[s1][s2] = 0.0; 
					for (int k1 = 0; k1 < nK; k1++)
						tDoubleSum[s1][s2] += tSumk[s1][k1][s2];
					tSums[s1] += tDoubleSum[s1][s2]; 
				}
				ts += tSums[s1]; 
			}

///////////////////////////////////////////////////////////////
			double tj = 0; 
			for (int j1 = 0, s1 = 0; s1 < nS; s1++)
			{ 
				jump2s[s1] = 0; 
				for (int k1 = 0; k1 < nK; k1++, j1++)  
				{
					jump2j[j1] = 0; 
					for (int j2 = 0, s2 = 0; s2 < nS; s2++)
					for (int k2 = 0; k2 < nK; k2++, j2++)
					{                                                
						double bj = prGk[k1][k2] * m_psi[m][j1*nsk+j2]; 
						double nsns = (1.0 - rs[m]) * (1.0 - rs[m]) * (1.0 - rk[m]) * bj; 

						double jk0 = nsns * (1.0 - rk[m]) * m_phi[m-1][j1*nsk+j2];
						double jk1j1 = nsns * rk[m] * beta[m][s1][k1] * tSumk[s2][k2][s1];
						double jk1j2 = nsns * rk[m] * beta[m][s2][k2] * tSumk[s1][k1][s2]; 
						double jk2 = nsns / (1.0 - rk[m]) * rk[m] * rk[m] * beta[m][s1][k1] * beta[m][s2][k2] * tDoubleSum[s1][s2];
						double t0 = jk0 + jk1j1 + jk1j2 + jk2; 

						double tt = rs[m] * (1.0 - rs[m]) * bj;
						double t1j1 = (rk[m] * beta[m][s1][k1] * tSums[s1] + (1.0 - rk[m]) * tmarginal[s1][k1]) * admix[s2] * beta[m][s2][k2] * tt; 
						double t1j2 = (rk[m] * beta[m][s2][k2] * tSums[s2] + (1.0 - rk[m]) * tmarginal[s2][k2]) * admix[s1] * beta[m][s1][k1] * tt; 
						double t1 = t1j1 + t1j2; 
						double t2 = ts * (rs[m] * admix[s1] * beta[m][s1][k1]) * (rs[m] * admix[s2] * beta[m][s2][k2]) * bj;  
						tj += t0 + t1 + t2;
						jump2j[j1] += 2 * (jk1j1 + jk2 + t1j1); 
						jump2s[s1] += 2 * (t1j1 + t2); 
					} 
				}
			}

			for (int j = 0; j < nsk; j++)
  		    	sjmsk[m][j] += jump2j[j] / tj * get_weight(); 
			for (int s = 0; s < nS; s++)
				sjms[m][s] += jump2s[s] / tj * get_weight(); 
		}
		delete[] jump2j;
		delete[] jump2s;
		Free3DMatrix(tSumk); 
		Free2DMatrix(tDoubleSum); 
		Free2DMatrix(tmarginal); 

		for (int j1 = 0; j1 < nsk; j1++)
			for (int j2 = 0; j2 <= j1; j2++)
			{
				prz2[j1][j2] = m_phi[0][j1*nsk+j2] * m_psi[0][j1*nsk+j2]; 
				prz2[j2][j1] = prz2[j1][j2]; 
			}
		normalize(&prz2[0][0], nsk * nsk); 

		for (int j1 = 0, s1 = 0; s1 < nS; s1++)
			for (int k1 = 0; k1 < nK; k1++, j1++)
				for (int j2 = 0, s2 = 0; s2 < nS; s2++)
					for (int k2 = 0; k2 < nK; k2++, j2++)
						sjmsk[0][j1] += prz2[j1][j2] * get_weight(); 
	}
	normalize(tadmix, nS); 
	for (int s = 0; s < nS; s++)
		admix[s] = tadmix[s]; 
//	for (int s = 0; s < nS; s++)
//		cout << admix[s] << " "; 
//	cout << endl; 

	delete[] tadmix; 
	Free2DMatrix(prGk); 
	Free2DMatrix(prz2); 
	Free2DMatrix(prz);             
	delete[] tSums;  
}	

void DipInd::forward_backward_block2(const int nS, const int nK, ModelParam * pMP, int which, int b, vector<int> pbreak)  
{
	int	mbeg = pbreak.at(b); 
	int	mend = pbreak.at(b+1); 
	int nsize = mend - mbeg; 
		
	if(m_admix == NULL) 
	{
		m_admix = new double[nS]; 
		for (int s = 0; s < nS; s++)
			m_admix[s] = 1.0 / nS; 
	}
	double * admix = m_admix; 

	double * rk = pMP->Getrk();
	double * rs = pMP->Getrs();
	double ** theta = pMP->Gettheta(); 
	double *** beta = pMP->Getbeta(); 
	
	int nsk = nS * nK; 
	int nsk2 = nsk * nsk; 
	double ** prz2 = Allocate2DMatrix(nsk, nsk); 
	double ** prz = Allocate2DMatrix(nK, nK); 
	double ** prGk = Allocate2DMatrix(nK, nK); 
	double * tSums = new double[nS];  

	/////////////////////////////////////////////////////////////////////////////////////////////
	if(b == 0) {
		vector<vector<double> > ().swap(vvfwd); 
		vector<vector<double> > ().swap(vvbwd); 
		vector<double > tv; 
		for (int j = 0; j < nsk2; j++)
			tv.push_back(-1); 
		for (unsigned i = 0; i < pbreak.size()-1; i++)
		{
			vvfwd.push_back(tv); 
			vvbwd.push_back(tv); 
		}
		vector<double > ().swap(tv); 
		//intialize the block boundary; 
		
		double ** psi = Allocate2DMatrix(nsk, nsk); 
		for (int j1 = 0; j1 < nsk; j1++) 
			for (int j2 = 0; j2 < nsk; j2++) 
				psi[j1][j2] = 1.0;

		for (unsigned i = pbreak.size()-1; i > 1; i--)
		{
			for (int m = pbreak.at(i)-1; m >= pbreak.at(i-1); m--)
			{
//				cout <<  m << " "; 
				for (int k1 = 0; k1 < nK; k1++)
				for (int k2 = 0; k2 < nK; k2++)
					prGk[k1][k2] = prG2(theta[m][k1],theta[m][k2], GetsnpGT(m)); 

				for (int j1 = 0, s1 = 0; s1 < nS; s1++)
				for (int k1 = 0; k1 < nK; k1++, j1++)
				{
					double ts = 0; 
					for (int j2 = 0, s2 = 0; s2 < nS; s2++)
					{
						tSums[s2] = 0; 
						for (int k2 = 0; k2 < nK; k2++, j2++)
							tSums[s2] += prGk[k1][k2] * psi[j1][j2] * beta[m][s2][k2]; 
						ts += tSums[s2] * admix[s2];
					}
					for (int j2 = 0, s2 = 0; s2 < nS; s2++)
					for (int k2 = 0; k2 < nK; k2++, j2++)
					{ 
						double t2 = (1.0 - rk[m]) * prGk[k1][k2] * psi[j1][j2] + rk[m] * tSums[s2]; 
						t2 *= (1.0 - rs[m]); 
						t2 += rs[m] * ts;
						prz2[j1][j2] = t2;
					}
				}   

				for (int j2 = 0, s2 = 0; s2 < nS; s2++)
				for (int k2 = 0; k2 < nK; k2++, j2++)
				{
					double ts = 0; 
					for (int j1 = 0, s1 = 0; s1 < nS; s1++)
					{	
						tSums[s1] = 0; 
						for (int k1 = 0; k1 < nK; k1++, j1++)
							tSums[s1] += prz2[j1][j2] * beta[m][s1][k1]; 
						ts += tSums[s1] * admix[s1]; 
					}
					for (int j1 = 0, s1 = 0; s1 < nS; s1++)
					for (int k1 = 0; k1 < nK; k1++, j1++)  
					{
						double t1 = (1.0 - rk[m]) * prz2[j1][j2] + rk[m] * tSums[s1];
						t1 *= (1.0 - rs[m]); 
						t1 += rs[m] * ts; 
						psi[j1][j2] = t1; 
					}
				}
				normalize(&psi[0][0], nsk2); 
			}   //it's a quite surprise that backward recursion can be done in position.                       

			for (int j1 = 0; j1 < nsk; j1++)
			for (int j2 = 0; j2 < nsk; j2++)
				vvbwd.at(i-2).at(j1*nsk+j2) = psi[j1][j2]; 
		}
		Free2DMatrix(psi); 
//		for (int b = 0; b < pbreak.size()-1; b++)
//		{
//			cout << b << " ==== "; 
//			for (int j = 0; j < nsk; j++)
//				cout << vvbwd.at(b).at(j) << " "; 
//			cout << endl; 
//		}
//		exit(0); 
	}
	///////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////
	if(m_psi == NULL)
		m_psi = Allocate2DMatrix(nsize, nsk*nsk); 

	if(mend == nLoci) {
		for (int j = 0; j < nsk2; j++) 
			m_psi[nsize-1][j] = 1.0;
	}
	else {
		for (int j = 0; j < nsk2; j++)
			m_psi[nsize-1][j] = vvbwd.at(b).at(j); 
	}

	for (int n = nsize - 1; n > 0; n--)
	{
		int m = n + mbeg; 
		for (int k1 = 0; k1 < nK; k1++)
		for (int k2 = 0; k2 < nK; k2++)
		   	prGk[k1][k2] = prG2(theta[m][k1],theta[m][k2], GetsnpGT(m)); 

		for (int j1 = 0, s1 = 0; s1 < nS; s1++)
		for (int k1 = 0; k1 < nK; k1++, j1++)
		{
			double ts = 0; 
			for (int j2 = 0, s2 = 0; s2 < nS; s2++)
			{
				tSums[s2] = 0; 
				for (int k2 = 0; k2 < nK; k2++, j2++)
					tSums[s2] += prGk[k1][k2] * m_psi[n][j1*nsk+j2] * beta[m][s2][k2]; 
				ts += tSums[s2] * admix[s2];
			}
   			for (int j2 = 0, s2 = 0; s2 < nS; s2++)
   			for (int k2 = 0; k2 < nK; k2++, j2++)
   			{ 
				double t2 = (1.0 - rk[m]) * prGk[k1][k2] * m_psi[n][j1*nsk+j2] + rk[m] * tSums[s2]; 
				t2 *= (1.0 - rs[m]); 
				t2 += rs[m] * ts;
				prz2[j1][j2] = t2;
			}
		}   

		for (int j2 = 0, s2 = 0; s2 < nS; s2++)
		for (int k2 = 0; k2 < nK; k2++, j2++)
		{
			double ts = 0; 
			for (int j1 = 0, s1 = 0; s1 < nS; s1++)
			{	
				tSums[s1] = 0; 
				for (int k1 = 0; k1 < nK; k1++, j1++)
					tSums[s1] += prz2[j1][j2] * beta[m][s1][k1]; 
				ts += tSums[s1] * admix[s1]; 
			}
			for (int j1 = 0, s1 = 0; s1 < nS; s1++)
			for (int k1 = 0; k1 < nK; k1++, j1++)  
			{
				double t1 = (1.0 - rk[m]) * prz2[j1][j2] + rk[m] * tSums[s1];
				t1 *= (1.0 - rs[m]); 
				t1 += rs[m] * ts; 
				m_psi[n-1][j1*nsk+j2] = t1; 
			}
		}

		normalize(m_psi[n-1], nsk2); 
	}                         
	
	///////////////////////////////////////////////////////////////////////////////
	
	if(m_phi == NULL)
		m_phi = Allocate2DMatrix(nsize, nsk*nsk); 

	if(b == 0) {
		for (int k1 = 0; k1 < nK; k1++)
		for (int k2 = 0; k2 < nK; k2++)
			prGk[k1][k2] = prG2(theta[0][k1],theta[0][k2], GetsnpGT(0)); 
		for (int j1 = 0, s1 = 0; s1 < nS; s1++)
		for (int k1 = 0; k1 < nK; k1++, j1++)  
			for (int j2 = 0, s2 = 0; s2 < nS; s2++)
			for (int k2 = 0; k2 < nK; k2++, j2++)
				m_phi[0][j1*nsk+j2] = admix[s1] * admix[s2] * beta[0][s1][k1] * beta[0][s2][k2] * prGk[k1][k2];
	}
	else {
		for (int j = 0; j < nsk2; j++)
			m_phi[0][j] = vvfwd.at(b).at(j); 
	}
	//at locus 0 -- leftmost; 
		
	for (int n = 0; n < nsize-1; n++)
	{
		int m = n + mbeg; 
//		cout << m << " "; 
		for (int k1 = 0; k1 < nK; k1++)
		for (int k2 = 0; k2 < nK; k2++)
			prGk[k1][k2] = prG2(theta[m+1][k1],theta[m+1][k2], GetsnpGT(m+1)); 

		for (int j1 = 0, s1 = 0; s1 < nS; s1++)
		for (int k1 = 0; k1 < nK; k1++, j1++)
		{
			double ts = 0; 
			for (int j2 = 0, s2 = 0; s2 < nS; s2++)
			{
				tSums[s2] = 0; 
				for (int k2 = 0; k2 < nK; k2++, j2++)
					tSums[s2] += m_phi[n][j1*nsk+j2]; 
				ts += tSums[s2];
			}
			for (int j2 = 0, s2 = 0; s2 < nS; s2++)
			for (int k2 = 0; k2 < nK; k2++, j2++)
			{ 
				double t2 = (1.0 - rk[m+1]) * m_phi[n][j1*nsk+j2] + rk[m+1] * beta[m+1][s2][k2] * tSums[s2];
				t2 *= (1.0 - rs[m+1]); 
                t2 += rs[m+1] * ts * admix[s2] * beta[m+1][s2][k2]; 
				prz2[j1][j2] = t2; 
			} 
		}

		for (int j2 = 0, s2 = 0; s2 < nS; s2++)
		for (int k2 = 0; k2 < nK; k2++, j2++)
		{
			double ts = 0; 
			for (int j1 = 0, s1 = 0; s1 < nS; s1++)
			{	
				tSums[s1] = 0; 
				for (int k1 = 0; k1 < nK; k1++, j1++)
					tSums[s1] += prz2[j1][j2]; 
				ts += tSums[s1]; 
			}
			for (int j1 = 0, s1 = 0; s1 < nS; s1++)
			for (int k1 = 0; k1 < nK; k1++, j1++)  
			{
				double t1 = (1.0 - rk[m+1]) * prz2[j1][j2] + rk[m+1] * beta[m+1][s1][k1] * tSums[s1];
				t1 *= (1.0 - rs[m+1]); 
				t1 += rs[m+1] * ts * admix[s1] * beta[m+1][s1][k1]; 
				m_phi[n+1][j1*nsk+j2] = t1 * prGk[k1][k2]; 
			}
		}

		normalize(m_phi[n+1], nsk2);
	}                                                          
	if(mend < nLoci)  // do one more step; 
	{
		double * phi1 = new double[nsk2]; 
		int n = nsize-1; 
		int m = n + mbeg; 

//		cout << " + " << m << " + " << b+1 << endl; 
		for (int k1 = 0; k1 < nK; k1++)
		for (int k2 = 0; k2 < nK; k2++)
			prGk[k1][k2] = prG2(theta[m+1][k1],theta[m+1][k2], GetsnpGT(m+1)); 

		for (int j1 = 0, s1 = 0; s1 < nS; s1++)
		for (int k1 = 0; k1 < nK; k1++, j1++)
		{
			double ts = 0; 
			for (int j2 = 0, s2 = 0; s2 < nS; s2++)
			{
				tSums[s2] = 0; 
				for (int k2 = 0; k2 < nK; k2++, j2++)
					tSums[s2] += m_phi[n][j1*nsk+j2]; 
				ts += tSums[s2];
			}
			for (int j2 = 0, s2 = 0; s2 < nS; s2++)
			for (int k2 = 0; k2 < nK; k2++, j2++)
			{ 
				double t2 = (1.0 - rk[m+1]) * m_phi[n][j1*nsk+j2] + rk[m+1] * beta[m+1][s2][k2] * tSums[s2];
				t2 *= (1.0 - rs[m+1]); 
                t2 += rs[m+1] * ts * admix[s2] * beta[m+1][s2][k2]; 
				prz2[j1][j2] = t2; 
			} 
		}

		for (int j2 = 0, s2 = 0; s2 < nS; s2++)
		for (int k2 = 0; k2 < nK; k2++, j2++)
		{
			double ts = 0; 
			for (int j1 = 0, s1 = 0; s1 < nS; s1++)
			{	
				tSums[s1] = 0; 
				for (int k1 = 0; k1 < nK; k1++, j1++)
					tSums[s1] += prz2[j1][j2]; 
				ts += tSums[s1]; 
			}
			for (int j1 = 0, s1 = 0; s1 < nS; s1++)
			for (int k1 = 0; k1 < nK; k1++, j1++)  
			{
				double t1 = (1.0 - rk[m+1]) * prz2[j1][j2] + rk[m+1] * beta[m+1][s1][k1] * tSums[s1];
				t1 *= (1.0 - rs[m+1]); 
				t1 += rs[m+1] * ts * admix[s1] * beta[m+1][s1][k1]; 
				phi1[j1*nsk+j2] = t1 * prGk[k1][k2]; 
			}
		}
		normalize(phi1, nsk2);

		for (int j = 0; j < nsk2; j++) 
			vvfwd.at(b+1).at(j) = phi1[j]; 
//		cout << b+1 << "+++++++ " << endl; 
//		for (int j = 0; j < nsk; j++)
//			cout << vvfwd.at(b+1).at(j) << " "; 
//		cout << endl; 
//		for (int j = 0; j < nsk; j++)
//			cout << vvbwd.at(b+1).at(j) << " "; 
//		cout << endl; 
		delete[] phi1; 
	} 
	///////////////////////////////////////////////////////////////////////////////////
	for (int n = 0; n < nsize; n++)
	{
		for (int j1 = 0; j1 < nsk; j1++)
		for (int j2 = 0; j2 < nsk; j2++)
		   	prz2[j1][j2] = m_phi[n][j1*nsk+j2] * m_psi[n][j1*nsk+j2]; 
		normalize(&prz2[0][0], nsk * nsk); 

		if(m_ps != NULL) {
			for (int j1 = 0, s1 = 0; s1 < nS; s1++)
				for (int k1 = 0; k1 < nK; k1++, j1++)
					for (int j2 = 0, s2 = 0; s2 < nS; s2++)
						for (int k2 = 0; k2 < nK; k2++, j2++)
							m_ps[n][s1*nS+s2] += prz2[j1][j2]; 
		}
		if(m_pk != NULL) {
			for (int j1 = 0, s1 = 0; s1 < nS; s1++)
				for (int k1 = 0; k1 < nK; k1++, j1++)
					for (int j2 = 0, s2 = 0; s2 < nS; s2++)
						for (int k2 = 0; k2 < nK; k2++, j2++)
							m_pk[n][k1] += 2*prz2[j1][j2]; 
		}
	}

	Free2DMatrix(prGk); 
	Free2DMatrix(prz2); 
	Free2DMatrix(prz);             
	delete[] tSums;  
}	

void DipInd::forward_backward_block1(const int nS, const int nK, ModelParam * pMP, int which, int b, vector<int> pbreak)  
{
	int	mbeg = pbreak.at(b); 
	int	mend = pbreak.at(b+1); 
	int len = mend - mbeg; 

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

	if(m_admix == NULL) 
	{
		m_admix = new double[nS]; 
		for (int s = 0; s < nS; s++)
			m_admix[s] = 1.0 / nS; 
	}
	if(m_admix1 == NULL) 
	{
		m_admix1 = new double[nS]; 
		for (int s = 0; s < nS; s++)
			m_admix1[s] = 0; 
	}
	double * admix = m_admix; 
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
				if(which == 0) {
					for (int k = 0; k < nK; k++)
						prGk[k] = prG1(m_maf[m], theta[m][k], GetsnpGT(m));
				} else {
					for (int k = 0; k < nK; k++)
						prGk[k] = prG(m_maf[m], theta[m][k], GetsnpGT(m));
				}

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

	if(m_psi == NULL)
		m_psi = Allocate2DMatrix(len, nsk); 

	for (int j = 0; j < nsk; j++)
		m_psi[len-1][j] = vvbwd.at(b).at(j); 
	//do the recursive calc backwards; 
	for (int p = len-1; p > 0; p--)
	{
		int m = pbreak.at(b) + p; 
		if(which == 0) {
			for (int k = 0; k < nK; k++)
				prGk[k] = prG1(m_maf[m], theta[m][k], GetsnpGT(m));
		} else {
			for (int k = 0; k < nK; k++)
				prGk[k] = prG(m_maf[m], theta[m][k], GetsnpGT(m));
		}

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
		if(which == 0) {
			for (int k = 0; k < nK; k++)
				prGk[k] = prG1(m_maf[0], theta[0][k], GetsnpGT(0));
		} else {
			for (int k = 0; k < nK; k++)
				prGk[k] = prG(m_maf[0], theta[0][k], GetsnpGT(0));
		}

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
	for (int j = 0; j < nsk; j++)
		m_phi[0][j] /= ts; 

	for (int p = 0; p < len-1; p++)
	{
		int m = pbreak.at(b) + p; 
		if(which == 0) {
			for (int k = 0; k < nK; k++)             	
				prGk[k] = prG1(m_maf[m+1], theta[m+1][k], GetsnpGT(m+1)); 
		} else {
			for (int k = 0; k < nK; k++)             	
				prGk[k] = prG(m_maf[m+1], theta[m+1][k], GetsnpGT(m+1)); 
		}

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
		ts = 0; 
		for (int j = 0; j < nsk; j++)
			ts += m_phi[p+1][j]; 
		for (int j = 0; j < nsk; j++)
			m_phi[p+1][j] /= ts; 
	}                                                          

	//do one more step to obtain the starting fwd of the next block; 
	if(b < (int) pbreak.size()-2)   //not last block
	{
		for (int j = 0; j < nsk; j++)
			vvfwdp.at(b+1).at(j) = m_phi[len-1][j]; 
		//this is needed to compute expected jumps at the block ends; 
		double * phi1 = new double[nsk]; 
		int p = len-1; 
		int m = pbreak.at(b+1)-1; 
		if(which == 0) {
			for (int k = 0; k < nK; k++)
				prGk[k] = prG1(m_maf[m+1], theta[m+1][k], GetsnpGT(m+1)); 
		}
		else {
			for (int k = 0; k < nK; k++)
				prGk[k] = prG(m_maf[m+1], theta[m+1][k], GetsnpGT(m+1)); 
		}

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
//		normalize(phi1, nsk);
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


	////////////////////////////////////////////////////////////////////////////
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
			{
				prz[k] += prz2[j];
				if(which > 0) m_admix1[s] += prz2[j]; 
			}

		if(which > 0 && m_ps != NULL) 
		{
			for (int j = 0, s = 0; s < nS; s++)
				for (int k = 0; k < nK; k++, j++)
					m_ps[m][s] += prz2[j];
		}

		if(which > 0 && m_pk != NULL) 
		{
			for (int j = 0, s = 0; s < nS; s++)
				for (int k = 0; k < nK; k++, j++)
					m_pk[m][k] += prz2[j];
		}

//			if(which == 1) 
		{
			m_maf[m] = 0; 
			for (int k = 0; k < nK; k++)
				m_maf[m] += prz[k] * theta[m][k]; 
		}
	}

	///////////////////////////////////////////////////////////////////////////
	delete[] prz; 
	delete[] prz2; 
	delete[] tSum; 
	delete[] prGk; 
}

void DipInd::forward_backward_block(const int nS, const int nK, ModelParam * pMP, double ** stop, double ** sbot, double ** sjmsk, double ** sjms, int which, int b, vector<int> pbreak) 
{
	//[pbreak.at(0), pbreak.at(1)-1]
	//[pbreak.at(1), pbreak.at(2)-1]
	//[pbreak.at(2), pbreak.at(3)-1]
	//

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

	double * admix = NULL; 
	if(which % 2 == 0) 
		admix = m_admix; 
	else 
		admix = m_admix + nS; 
	/////////////////////////////////////////////////////////////////////////////////////
//	if(which == 0 && b == 0) 
//	{
//		double * pmar1 = new double[nsk]; //current; 
//		double * pmar2 = new double[nsk]; //next; 
//
//		for (int j = 0, s = 0; s < nS; s++)
//			for (int k = 0; k < nK; k++, j++)
//				pmar1[j] = admix[s] * beta[0][s][k]; 
//		normalize(pmar1, nsk); 
//		for (int k = 0; k < nK; k++)
//			prz[k] = 0; 
//		for (int j = 0, s = 0; s < nS; s++)
//			for (int k = 0; k < nK; k++, j++)
//				prz[k] += pmar1[j]; 
//		m_maf[0] = 0; 
//		for (int k = 0; k < nK; k++)
//			m_maf[0] += prz[k] * theta[0][k];
//
//		for (int m = 0; m < nLoci-1; m++)
//		{
//			double ts = 0; 
//			for (int j = 0, s = 0; s < nS; s++)
//			{
//				tSum[s] = 0; 
//				for (int k = 0; k < nK; k++, j++)
//					tSum[s] += pmar1[j];
//				ts += tSum[s]; 
//			}
//			for (int j = 0, s = 0; s < nS; s++)
//				for (int k = 0; k < nK; k++, j++)
//				{
//					pmar2[j] = rs[m+1] * ts * admix[s] * beta[m+1][s][k];
//					pmar2[j] += (1 - rs[m+1]) * (rk[m+1] * beta[m+1][s][k] * tSum[s] + (1 - rk[m+1]) * pmar1[j]); 
//				}
//			for (int j = 0; j < nsk; j++)
//				pmar1[j] = pmar2[j]; 
//			normalize(pmar1, nsk); 
//
//			for (int k = 0; k < nK; k++)
//				prz[k] = 0; 
//			for (int j = 0, s = 0; s < nS; s++)
//				for (int k = 0; k < nK; k++, j++)
//					prz[k] += pmar1[j]; 
//			m_maf[m+1] = 0; 
//			for (int k = 0; k < nK; k++)
//				m_maf[m+1] += prz[k] * theta[m][k];
//		}
//		delete[] pmar1; 
//		delete[] pmar2; 
//	}

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
				if(which == 0) {
					for (int k = 0; k < nK; k++)
						prGk[k] = prG1(m_maf[m], theta[m][k], GetsnpGT(m));
				} else {
					for (int k = 0; k < nK; k++)
						prGk[k] = prG(m_maf[m], theta[m][k], GetsnpGT(m));
				}

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
		if(which == 0) {
			for (int k = 0; k < nK; k++)
				prGk[k] = prG1(m_maf[m], theta[m][k], GetsnpGT(m));
		} else {
			for (int k = 0; k < nK; k++)
				prGk[k] = prG(m_maf[m], theta[m][k], GetsnpGT(m));
		}

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
		if(which == 0) {
			for (int k = 0; k < nK; k++)
				prGk[k] = prG1(m_maf[0], theta[0][k], GetsnpGT(0));
		} else {
			for (int k = 0; k < nK; k++)
				prGk[k] = prG(m_maf[0], theta[0][k], GetsnpGT(0));
		}

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
	for (int j = 0; j < nsk; j++)
		m_phi[0][j] /= ts; 
	double scale = log(ts); 
	for (int p = 0; p < len-1; p++)
	{
		int m = pbreak.at(b) + p; 
		if(which == 0) {
			for (int k = 0; k < nK; k++)             	
				prGk[k] = prG1(m_maf[m+1], theta[m+1][k], GetsnpGT(m+1)); 
		} else {
			for (int k = 0; k < nK; k++)             	
				prGk[k] = prG(m_maf[m+1], theta[m+1][k], GetsnpGT(m+1)); 
		}

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
		ts = 0; 
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
		if(which == 0) {
			for (int k = 0; k < nK; k++)
				prGk[k] = prG1(m_maf[m+1], theta[m+1][k], GetsnpGT(m+1)); 
		}
		else {
			for (int k = 0; k < nK; k++)
				prGk[k] = prG(m_maf[m+1], theta[m+1][k], GetsnpGT(m+1)); 
		}

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
//		normalize(phi1, nsk);
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
			if(which == 0) {
				for (int k = 0; k < nK; k++)
					prGk[k] = prG1(m_maf[m], theta[m][k], GetsnpGT(m)); 
			} else {
				for (int k = 0; k < nK; k++)
					prGk[k] = prG(m_maf[m], theta[m][k], GetsnpGT(m)); 
			}

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

			if(m == 0) {
				for (int j = 0; j < nsk; j++)
					prz2[j] = m_phi[0][j] * m_psi[0][j];
				normalize(prz2, nsk); 
				for (int j = 0, s = 0; s < nS; s++)
					for (int k = 0; k < nK; k++, j++)
						sjmsk[0][j] += prz2[j] * get_weight(); 
				}
			else {
				for (int j = 0; j < nsk; j++)
					sjmsk[m][j] += prz2[j] / tj * get_weight(); 
				for (int j = 0; j < nS; j++)
					sjms[m][j] += prz[j] / tj * get_weight(); 
			}
		//the first locus, it's important for beta. 
		}
	}

	////////////////////////////////////////////////////////////////////////////
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
			{
				prz[k] += prz2[j];
				if(which > 0) m_admix1[s] += prz2[j]; 
			}

		if(m_ps != NULL) 
		{
			for (int j = 0, s = 0; s < nS; s++)
				for (int k = 0; k < nK; k++, j++)
					m_ps[m][s] += prz2[j];
		}

		if(m_pk != NULL) 
		{
			for (int j = 0, s = 0; s < nS; s++)
				for (int k = 0; k < nK; k++, j++)
					m_pk[m][k] += prz2[j];
		}
//		if(which == 1) 
		{
			m_maf[m] = 0; 
			for (int k = 0; k < nK; k++)
				m_maf[m] += prz[k] * theta[m][k]; 
		}

		if(stop != NULL && sbot != NULL) 
		{
			for (int k = 0; k < nK; k++)
				prz[k] *= get_weight(); 
			int gt = GetsnpGT(m); 
			if(gt == 0) 
			{
				for (int k = 0; k < nK; k++)
					sbot[m][k] += prz[k];
			}
			else if (gt == 2)
			{
				for (int k = 0; k < nK; k++)
				{
					stop[m][k] += prz[k];      
					sbot[m][k] += prz[k]; 
				}
			}
			else if (gt == 1) 
			{
				for (int k = 0; k < nK; k++)
				{
					double t1 = (1 - m_maf[m]) * theta[m][k];
					double t2 = m_maf[m] * (1.0 - theta[m][k]);					
					stop[m][k] += prz[k] * t1 / (t1 + t2); 
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

void DipInd::forward_backward_haploid(const int nS, const int nK, ModelParam * pMP, double ** stop, double ** sbot, double ** sjmsk, double ** sjms, int which) 
{
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

	if(m_admix == NULL) {
		m_admix = new double[nS]; 
		for (int s = 0; s < nS; s++)
			m_admix[s] = 1.0 / nS;
	}
	if(m_admix1 == NULL) {
		m_admix1 = new double[nS]; 
		for (int s = 0; s < nS; s++)
			m_admix1[s] = 0;
	}

	double * admix = m_admix; 
	double * pmaf = m_maf; 
	/////////////////////////////////////////////////////////////////////////////////////

	if(m_psi == NULL)
		m_psi = Allocate2DMatrix(nLoci, nsk); 
	//at the rightmost locus M; 
	for (int j = 0; j < nsk; j++) 
		m_psi[nLoci-1][j] = 1.0;
	//do the recursive calc backwards; 
	for (int m = nLoci - 1; m > 0; m--)
	{
		if(which == 0) {
			for (int k = 0; k < nK; k++)
				prGk[k] = prG1(pmaf[m], theta[m][k], GetsnpGT(m));
		}
		else {
			for (int k = 0; k < nK; k++)
				prGk[k] = prG(pmaf[m], theta[m][k], GetsnpGT(m));
		}


		ts = 0; 
		for (int j = 0, s = 0; s < nS; s++)
		{
			tSum[s] = 0.0;
			for (int k = 0; k < nK; k++, j++)
				tSum[s] +=  prGk[k] * m_psi[m][j] * beta[m][s][k]; 
			ts += tSum[s] * admix[s]; 
		}

		for (int j = 0, s = 0; s < nS; s++)
			for (int k = 0; k < nK; k++, j++)
			{
				double temp = (1.0 - rk[m]) * prGk[k] * m_psi[m][j] + rk[m] * tSum[s]; 
				temp *= (1.0 - rs[m]);
				m_psi[m-1][j] = temp + rs[m] * ts; 
			}
//		for (int j = 0; j < nsk; j++)
//			m_psi[m-1][j] += 1e-20; 
		ts = 0; 
		for (int j = 0; j < nsk; j++)
			ts += m_psi[m-1][j]; 
		normalize(m_psi[m-1], nsk); 
		for (int j = 0; j < nsk; j++)
		{
			if(isnan(m_psi[m-1][j])) {
				for (int j = 0; j < nsk; j++)
					cout << m_psi[m][j] << " "; 
				cout << endl; 
				cout << rk[m] << endl; 
				cout << "m_psi is nan " << endl ; exit(0); 
			}
		}
	}
	
	///////////////////////////////////////////////////////////////////////////////
	if(m_phi == NULL)
		m_phi = Allocate2DMatrix(nLoci, nsk); 

	//at locus 0 -- leftmost; 
	if(which == 0) {
		for (int k = 0; k < nK; k++)
			prGk[k] = prG1(pmaf[0], theta[0][k], GetsnpGT(0));
	}
	else {
		for (int k = 0; k < nK; k++)
			prGk[k] = prG(pmaf[0], theta[0][k], GetsnpGT(0));
	}

	for (int j = 0, s = 0; s < nS; s++)
		for (int k = 0; k < nK; k++, j++)
	   		m_phi[0][j] = admix[s] * beta[0][s][k] * prGk[k];
	
	ts = 0; 
	for (int j = 0; j < nsk; j++)
		ts += m_phi[0][j]; 
	for (int j = 0; j < nsk; j++)
		m_phi[0][j] /= ts;
	double scale = log(ts); 
	for (int m = 0; m < nLoci-1; m++)
	{
		if(which == 0) {
			for (int k = 0; k < nK; k++)
				prGk[k] = prG1(pmaf[m+1], theta[m+1][k], GetsnpGT(m+1)); 
		}
		else {
			for (int k = 0; k < nK; k++)
				prGk[k] = prG(pmaf[m+1], theta[m+1][k], GetsnpGT(m+1)); 
		}

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
				temp += rs[m+1] * (admix[s] * beta[m+1][s][k]) * ts; 
				m_phi[m+1][j] = temp * prGk[k]; 
			}   

//		for (int j = 0; j < nsk; j++)
//		   	m_phi[m+1][j] += 1e-20;
		ts = 0; 
		for (int j = 0; j < nsk; j++)
		   	ts += m_phi[m+1][j];
		for (int j = 0; j < nsk; j++)
			m_phi[m+1][j] /= ts;
		scale += log(ts); 
	}                                                          
	m_loglike += scale / log(10.0);

	{
		if(sjmsk != NULL) 
		{                         
			for (int m = 1; m < nLoci; m++)
			{
				if(which == 0) {
					for (int k = 0; k < nK; k++)
						prGk[k] = prG1(pmaf[m], theta[m][k], GetsnpGT(m)); 
				} else {
					for (int k = 0; k < nK; k++)
						prGk[k] = prG(pmaf[m], theta[m][k], GetsnpGT(m)); 
				}
				

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
					prz[s] = 0; 
					for (int k = 0; k < nK; k++, j++)
					{                    
						double bj = prGk[k] * m_psi[m][j]; 
						double s0k0 =  m_phi[m-1][j] * (1.0 - rs[m]) * (1.0 - rk[m]) * bj; 
						double s0k1 = tSum[s] * (1.0 - rs[m]) * rk[m] * beta[m][s][k] * bj;   
						double s1 = ts * rs[m] * admix[s] * beta[m][s][k] * bj;  
						prz2[j] = s0k1; 
						prz[s] += s1; 
						tj += s0k0 + s0k1 + s1; 
					}   
				}

				for (int j = 0; j < nsk; j++)
					sjmsk[m][j] += prz2[j] / tj * get_weight(); 
				for (int j = 0; j < nS; j++)
					sjms[m][j] += prz[j] / tj * get_weight(); 
			}

			for (int j = 0; j < nsk; j++)
				prz2[j] = m_phi[0][j] * m_psi[0][j];
			normalize(prz2, nsk); 
			for (int j = 0, s = 0; s < nS; s++)
				for (int k = 0; k < nK; k++, j++)
					sjmsk[0][j] += prz2[j] * get_weight(); 
			//the first locus, it's important for beta. 
		}
		////////////////////////////////////////////////////////////////////////////
		double * tadmix = new double[nS]; 
		for (int s = 0; s < nS; s++)
			tadmix[s] = 0; 
		for (int m = 0; m < nLoci; m++)
		{                                       
			for (int j = 0; j < nsk; j++)
				prz2[j] = (m_phi[m][j]) * (m_psi[m][j]); 
			normalize(prz2, nsk); 

			for (int k = 0; k < nK; k++)
				prz[k]= 0; 
			for (int j = 0, s = 0; s < nS; s++)
				for (int k = 0; k < nK; k++, j++)
				{
					prz[k] += prz2[j];
					tadmix[s] += prz2[j]; 
				}

			if(which > 0 && m_ps != NULL) 
			{
				for (int j = 0, s = 0; s < nS; s++)
					for (int k = 0; k < nK; k++, j++)
						m_ps[m][s] += prz2[j];
			}

			if(which > 0 && m_pk != NULL) 
			{
				for (int j = 0, s = 0; s < nS; s++)
					for (int k = 0; k < nK; k++, j++)
						m_pk[m][k] += prz2[j];
			}

//			if(which == 0) 
			{
				pmaf[m] = 0; 
				for (int k = 0; k < nK; k++)
					pmaf[m] += prz[k] * theta[m][k]; 
			}

			if(stop != NULL && sbot != NULL) 
			{
				for (int k = 0; k < nK; k++)
					prz[k] *= get_weight(); 
				int gt = GetsnpGT(m); 
				if(gt == 0) 
				{
					for (int k = 0; k < nK; k++)
						sbot[m][k] += prz[k];
				}
				else if (gt == 2)
				{
					for (int k = 0; k < nK; k++)
					{
						stop[m][k] += prz[k];      
						sbot[m][k] += prz[k]; 
					}
				}
				else if (gt == 1) 
				{
					for (int k = 0; k < nK; k++)
					{
						double t1 = (1 - pmaf[m]) * theta[m][k];
						double t2 = pmaf[m] * (1.0 - theta[m][k]);					
						stop[m][k] += prz[k] * t1 / (t1 + t2); 
						sbot[m][k] += prz[k]; 
					}
				}
			}

		}

		///////////////////////////////////////////////////////////////////////

		if(which > 0) {
			normalize(tadmix, nS); 
			for (int s = 0; s < nS; s++)
				m_admix1[s] += tadmix[s]; 
		}

		delete[] tadmix; 
	}
	///////////////////////////////////////////////////////////////////////////
	delete[] prz; 
	delete[] prz2; 
	delete[] tSum; 
	delete[] prGk; 
}

void DipInd::CalcAll(const int nS, const int nK, ModelParam * pMP, double ** stop, double ** sbot, double ** sjmsk, double ** sjms, const int mode, int which)
{
	int unit = 2000; 
	int nfile = ceil(nLoci / (double) unit); 
	vector<int> pbreak; 
	for (int i = 0; i < nfile; i++)
		pbreak.push_back(i * unit); 
	pbreak.push_back(nLoci); 
//	cout << nfile << " --- "; 
//	for (unsigned i = 0; i < pbreak.size(); i++)
//		cout << pbreak.at(i) << " " ; 
//	cout << endl; 

	if(m_admix == NULL) {
		m_admix = new double[nS]; 
		for (int s = 0; s < nS; s++)
			m_admix[s] = 1.0 / nS;
	}
   	if(m_admix1 == NULL) 
		m_admix1 = new double[nS]; 
	if(mode == 1) 
	{
		for (int s = 0; s < nS; s++)
			m_admix[s] += 0.0; 
		normalize(m_admix, nS); 

		if(m_maf == NULL) 
		{
			m_maf = new double[nLoci]; 
			for (int m = 0; m < nLoci; m++)
				m_maf[m] =  gsl_ran_beta(gsl_r, 1000, 1000); 
		}

		m_loglike = 0; 
//		if(m_admix1 == NULL) m_admix1 = new double[nS]; 
////		for (int s = 0; s < nS; s++)
////			m_admix1[s] = 0; 
////		for (unsigned b = 0; b < pbreak.size() - 1; b++) 
////			forward_backward_block(nS, nK, pMP, NULL, NULL, NULL, NULL, 0, b, pbreak); 
////		normalize(m_admix1, nS); 
////		for (int s = 0; s < nS; s++)
////			m_admix[s] = m_admix1[s]; 
//
//		for (int s = 0; s < nS; s++)
//			m_admix1[s] = 0; 
//		for (unsigned b = 0; b < pbreak.size() - 1; b++) 
//			forward_backward_block(nS, nK, pMP, stop, sbot, sjmsk, sjms, 0, b, pbreak); 
////		normalize(m_admix1, nS); 
////		for (int s = 0; s < nS; s++)
////			m_admix[s] = m_admix1[s]; 
////		for (int s = 0; s < nS; s++)
////			m_admix1[s] = 0; 
//		for (unsigned b = 0; b < pbreak.size() - 1; b++) 
//			forward_backward_block(nS, nK, pMP, stop, sbot, sjmsk, sjms, 1, b, pbreak); 
//		normalize(m_admix1, nS); 
//		for (int s = 0; s < nS; s++)
//			m_admix[s] = m_admix[nS+s] = m_admix1[s]; 

		forward_backward_haploid(nS, nK, pMP, NULL, NULL, NULL, NULL, 0); 
		m_loglike = 0; 
		for (int s = 0; s < nS; s++)
			m_admix1[s] = 0; 
		forward_backward_haploid(nS, nK, pMP, stop, sbot, sjmsk, sjms, 1); 
		forward_backward_haploid(nS, nK, pMP, stop, sbot, sjmsk, sjms, 2); 
		for (int s = 0; s < nS; s++)
			m_admix[s] = m_admix1[s]; 
		normalize(m_admix, nS); 
	} 
	else 
	{
//		for (int s = 0; s < nS; s++)
//			m_admix[s] += 0.01; 
//		normalize(m_admix, nS); 
        forward_backward_diploid(nS, nK, pMP, stop, sbot, sjmsk, sjms, which); 
	}
}
