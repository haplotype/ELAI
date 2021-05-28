#ifndef __ANCHAP_H__
#define __ANCHAP_H__

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "fp.h"
#include "indiv.h"
#include "model.h"

using namespace std;

class AncHap : public Individual
{
private:
	double ** m_phi; 						// forward MxK
	double ** m_psi; 					// backward MxK
	double *  snp_dstr; 	             //nx2; 
	double * af; 
	double ** m_jms; 						
public:
	AncHap();
	~AncHap();
	
	double * get_admix(void) {return m_admix;}
	double ** get_jms(void) {return m_jms;}
	void allocate_af(int s) {af = new double[s];}
	void set_af(int m, double s) {af[m] = s;}
	double get_af(int m) {return af[m];}

	virtual void CalcAll(const int, const int, class ModelParam *, double **, double **, double **, double **, int, int);
	virtual void clean_after_em(int);   
	virtual void compute_ps(const int, const int, ModelParam *, const int){;};
	virtual void compute_ps(int, int, ModelParam*, int, int, vector<int>){;};
	virtual void compute_pk(int, int, ModelParam*, int){;};
	virtual void compute_pk(int, int, ModelParam*, int, int, vector<int>){;};
	virtual void forward_backward_block(const int, const int, class ModelParam *, double **, double **, double **, double **, int, int, vector<int>);	
	
};

#endif 
