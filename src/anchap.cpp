#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "anchap.h"
#include "model.h"
#include <iostream>
using namespace std;

AncHap::AncHap(void)
{
	m_phi = NULL;
	m_psi = NULL;
	m_top = NULL;
	m_bot = NULL; 
	m_lod = NULL;
	m_lps = NULL; 
	m_ploid = 1; 
	m_jms = NULL; 
	m_admix = NULL; 
}

AncHap::~AncHap(void)
{
	if(m_phi) {Free2DMatrix(m_phi);    m_phi = NULL; }
	if(m_psi) {Free2DMatrix(m_psi);       m_psi = NULL;}
	if(m_top) {Free2DMatrix(m_top); m_top = NULL;}
	if(m_bot) {Free2DMatrix(m_bot); m_bot = NULL;}
}

double prG(double th, double af )
{
//	 return ((1.0 - af) * (1.0 - th) + af * th);
//	if(res > 0.999) res = 0.999; 
//	else if (res < 0.001) res = 0.001; 
	double res = gsl_ran_beta_pdf(af, global_Ne * th, global_Ne * (1.0-th)); 
	return(res); 
}

void AncHap::CalcAll(const int nLoci, const int nS, ModelParam * pMD, double ** stop, double ** sbot, double ** slod, double ** slps, int dk, int n2)
{
	double * r = pMD->Getru(); 
	double ** eta = pMD->Geteta();
	if(m_admix == NULL)  
	{
		m_admix = new double[nS]; 
		for (int s = 0; s < nS; s++)
			m_admix[s] = 1.0 / nS; 
	}
	double * admix = m_admix; 
//	for (int s = 0; s < nS; s++)
//	    admix[s] = 1.0 / nS; 
//    normalize(admix, nS); 
	
	if(m_psi == NULL)
		m_psi = Allocate2DMatrix(nLoci, nS); 
	for (int s = 0; s < nS; s++)
		(m_psi)[nLoci-1][s] = 1.0;
	for (int m = nLoci - 1; m > 0; m--)
	{
		double ts = 0; 
		for (int s = 0; s < nS; s++)
            ts += prG(eta[m][s], get_af(m)) * m_psi[m][s] * admix[s];
		for (int s = 0; s < nS; s++)
			m_psi[m-1][s] = (1.0 - r[m]) * prG(eta[m][s], get_af(m)) * m_psi[m][s] + r[m] * ts;
		
		ts = 0; 
		for (int s = 0; s < nS; s++)
			ts += m_psi[m-1][s];
		for (int s = 0; s < nS; s++)
			 m_psi[m-1][s] /= ts;
	}
	
	if(m_phi == NULL)
		m_phi = Allocate2DMatrix(nLoci, nS); 
	for (int s = 0; s < nS; s++)
		m_phi[0][s] = admix[s] * prG(eta[0][s], get_af(0));
	for (int m = 0; m < nLoci-1; m++)
	{
		double ts = 0; 
		for (int s = 0; s < nS; s++)
			ts += m_phi[m][s]; 
		for (int s = 0; s < nS; s++)
			m_phi[m+1][s] = ((1.0 - r[m+1]) * m_phi[m][s] + r[m+1] * admix[s] * ts) * prG(eta[m+1][s], get_af(m+1));
		
		ts = 0; 
		for (int s = 0; s < nS; s++)
			ts += m_phi[m+1][s];
		for (int s = 0; s < nS; s++)
			m_phi[m+1][s] /= ts;
	 }                                                          

	///////////////////////////////////////////////////////////////////////
	double * prs = new double[nS]; 
	double * tadmix = new double[nS]; 
	for (int s = 0; s < nS; s++)
		tadmix[s] = 0; 
	for (int m = 0; m < nLoci; m++)
	{
		for (int s = 0; s < nS; s++)
			prs[s] = m_phi[m][s] * m_psi[m][s]; 
		normalize(prs, nS); 
		
		for (int s = 0; s < nS; s++)
        	tadmix[s] += prs[s]; 

//		if(slod != NULL && slps != NULL) 
		{
			double af = get_af(m); 
			double lod = log(af/(1.0 - af));
			for (int s = 0; s < nS; s++)
			{
				slod[m][s] += lod * prs[s]; 
				slps[m][s] += prs[s]; 
			}
		}
//		if(sbot != NULL && stop != NULL) 
		{
			double top = 0; 
			double bot = 0; 
			for (int s = 0; s < nS; s++)
			{
				top += (global_Ne * eta[m][s] - 1.0) * prs[s]; 
				bot += (global_Ne - 2.0) * prs[s];
			}    
			stop[m][dk] += top; 
			sbot[m][dk] += bot; 
		}
	}

	//////////////////////////////////////////////////////////////////////
	if(m_pjs == NULL) 
		m_pjs = new double[nLoci]; 
	for (int m = 1; m < nLoci; m++)
	{
		double ts = 0; 
		for (int s = 0; s < nS; s++)
			ts += m_phi[m-1][s];

		double t0, t1; 
		t0 = t1 = 0; 
		for (int s = 0; s < nS; s++)
		{                                        
			double bs = prG(eta[m][s], get_af(m)) * m_psi[m][s];
			t0 += m_phi[m-1][s] * (1 - r[m]) * bs;
			t1 += ts * r[m] * admix[s] * bs;
		}
    	m_pjs[m] = t1 / (t1 + t0); 
	}

	normalize(tadmix, nS); 
	for (int s = 0; s < nS; s++)
	{
		if(tadmix[s] < 0.01) tadmix[s] = 0.01; 
		else if(tadmix[s] > 0.99) tadmix[s] = 0.99; 
		admix[s] = tadmix[s];
	}
	normalize(admix, nS); 

	delete[] prs; 
	delete[] tadmix; 
}
                              
void AncHap::clean_after_em(int deep)
{
	if(m_pjs) {delete[] m_pjs; m_pjs = NULL;} 
//	if(phiScale) {delete[] phiScale; phiScale = NULL;}
//	if(psiScale) {delete[] psiScale; psiScale = NULL;}
//	Free2DMatrix(m_lps); 	m_lps = NULL;
//	Free2DMatrix(m_lod); 	m_lod = NULL;
//    Free2DMatrix(m_top); 		m_top = NULL;
//    Free2DMatrix(m_bot); 		m_bot = NULL;
//    Free2DMatrix(m_jms); 		m_jms = NULL;
	if(deep) {
		Free2DMatrix(m_phi);    	m_phi = NULL; 
		Free2DMatrix(m_psi);       m_psi = NULL;
	}
}




void AncHap::forward_backward_block(const int nS, const int nK, ModelParam * pMP, double ** stop, double ** sbot, double ** sjmsk, double ** sjms, int which, int b, vector<int> pbreak) 
{
	; 
}
