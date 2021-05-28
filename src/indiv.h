#ifndef __INDIV_H__
#define __INDIV_H__

#include <iostream>
#include <string>
#include <vector> 
#include "fpmath.h"
using namespace std;

#define 	QQ		(9)

//the parent class of diploid and haploid; 
class Individual
{
protected:
	int nLoci; 
	int   m_ploid; 
	int   m_isPanel; 
	char * snpGT; 					    // 1 dimension array;
	double * m_mgt; 
	double * snp_dstr; 
	string		huskyID;           		// name or ID of the individual;
    double * 		phiScale;               // the scale to prevent phi from underflow; 
 	double * 		psiScale;              // the scale to prevent psi from underflow; 
	int         nMissingGT; 

    double 		m_loglike; 			// likelihood of genotype = log(m_phi(nLoci-1)); 
//	double * m_maf; 
	double m_mu; 
	double ** m_ps; 
	double ** m_pk; 
	double ** m_lod; 
	double ** m_lps; 
	double ** m_top; 
	double ** m_bot; 
	double * m_admix; 
	double * m_pjk; 
	double * m_pjs; 

	double m_temp; 
	int mask; 

	double * m_beta; 
	double weight; //this is to weight rare alleles in HLA allele imputation; 
public:
	Individual();					   	// construct function;
	Individual(int);					   	// construct function;
	virtual ~Individual();              // destruct function; 
	
	inline double * get_pjk(void) {return m_pjk;}
	inline double * get_pjs(void) {return m_pjs;}
    inline void SetID(char * str){huskyID.assign(str);}
	// set husky ID for each individual if available;
	inline void AllocatesnpGT(int len)
		{snpGT = (char*) Allocate1D(sizeof(char), len);}
	// allocate snpGT for pointer; 
    inline void SetsnpGT(int pos, char val) {snpGT[pos] = val;}
   	// for readdipdata;
 	inline string GetID(void){return huskyID;}
	inline short GetsnpGT(int pos) {return short(snpGT[pos]-'0');}
	inline void GetsnpGT(int nLoci, short * snp_gt) 
	{ 
		for (int m = 0; m < nLoci; m++)
			snp_gt[m] = (short) (snpGT[m] - '0'); 
	}
	//		memcpy(snp_gt, snpGT, sizeof(char) * nLoci);}
	inline double GetLikelihood(void){return m_loglike;}


	inline void Setploid(int s) {m_ploid = s;}
	inline void SetisPanel(int s) {m_isPanel = s;}
	inline int Getploid(void) {return m_ploid;}
	inline int GetisPanel(void) {return m_isPanel;}
             
	inline void set_mtemp(double s) {m_temp = s;}
	inline double * get_mbeta(void) {return m_beta;}
	void set_mu(double s) {m_mu = s;}
	double ** get_ps(void) {return m_ps;}
	double ** get_pk(void) {return m_pk;}
	double ** get_lod(void) {return m_lod;}
	double ** get_lps(void) {return m_lps;}
	double ** get_top(void) {return m_top;}
	double ** get_bot(void) {return m_bot;}

	inline void set_mask(void) {mask = 1;}
	inline int get_mask(void) {return mask;}
	inline void set_weight(double s) {weight = s;}
	inline double get_weight(void) {return weight;}
	// interface; 
	
	void set_admix(int a, int nS);
	void set_admix(int nS);
	void set_admix(int s, double v, int ns); 
	double * get_admix(void) {return m_admix;};

	virtual void CalcAll(const int, const int, class ModelParam *, double **, double **, double **, double **, int, int) = 0;	
	virtual void clean_after_em(int) = 0; 

	virtual void compute_ps(const int nS, const int nK, ModelParam * pMD, const int) = 0;
	virtual void compute_ps(const int nS, const int nK, ModelParam * pMD, const int, const int, vector<int>) = 0;
	virtual void compute_pk(const int nS, const int nK, ModelParam * pMD, const int) = 0;
	virtual void compute_pk(const int nS, const int nK, ModelParam * pMD, const int, const int, vector<int>) = 0;
	virtual void forward_backward_block(const int, const int, class ModelParam *, double **, double **, double **, double **, int, int, vector<int>) = 0;	
};
#endif
