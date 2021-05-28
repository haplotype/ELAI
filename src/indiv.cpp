#include <stdio.h>
#include <stdlib.h>

#include "indiv.h"
#include <iostream>
using namespace std;

Individual::Individual()
{
	snpGT = NULL; 
	phiScale = NULL; 
	psiScale = NULL; 
	nMissingGT = 0; 
	snp_dstr = NULL; 
	m_mgt = NULL; 
	m_ps = NULL; 
	m_pk = NULL; 
	m_lod = NULL; 
	m_lps = NULL; 
	m_top = NULL; 
	m_bot = NULL; 
	m_mu = 0.0; 
	m_admix = NULL; 
	m_pjk = NULL; 
	m_pjs = NULL; 
	m_temp = 1.0; 
	mask = 0; 
	weight=1.0; 
}


Individual::Individual(int n)
{               
	nLoci = n; 
	snpGT = NULL; 
	phiScale = NULL; 
	psiScale = NULL; 
	nMissingGT = 0; 
	snp_dstr = NULL; 
	m_mgt = NULL; 
	m_ps = NULL; 
	m_pk = NULL; 
	m_lod = NULL; 
	m_lps = NULL; 
	m_top = NULL; 
	m_bot = NULL; 
	m_mu = 0.0; 
	m_admix = NULL; 
	m_pjk = NULL; 
	m_pjs = NULL; 
	m_temp = 1.0; 
	mask = 0; 
	weight=1.0; 
}

Individual::~Individual()
{  
	huskyID.assign("\0"); 
	if(snpGT) {Free1D((void*) snpGT); snpGT = NULL;}
	if(phiScale) {delete[] phiScale;       phiScale = NULL;}
	if(psiScale) {delete[] psiScale;     psiScale = NULL;} 
}

void Individual::set_admix(int a, int nS)
{
	//this only works for panel; 
	if(m_admix == NULL) 
		m_admix = new double[nS*2]; 
	if(a >= 0  && a < nS) 
	{
		for (int s = 0; s < nS*2; s++)
			m_admix[s] = 0.0; 
		m_admix[a] = m_admix[nS+a] = 1.0; 
		normalize(m_admix, nS); 
		normalize(m_admix+nS, nS); 
	}
}

void Individual::set_admix(int nS)
{
	//this only works for panel; 
	if(m_admix == NULL) 
		m_admix = new double[nS*2]; 
	for (int s = 0; s < nS*2; s++)
		m_admix[s] = 1.0 / nS; 
}     

void Individual::set_admix(int s, double v, int ns)
{
	if(m_admix == NULL) 
		m_admix=new double[ns*2]; 
	m_admix[s]= m_admix[ns+s] = v;
}

