#ifndef __HAPLOID_H__
#define __HAPLOID_H__

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "fp.h"
#include "indiv.h"
#include "model.h"

using namespace std;

class HapInd : public Individual
{
private:
	double ** m_phi; 						// forward MxK
	double ** m_psi; 					// backward MxK

	vector<vector<double > > vvfwdp; 
	vector<vector<double > > vvfwd; 
	vector<vector<double > > vvbwd; 
public:
	HapInd();
	HapInd(int);
    ~HapInd();
	
	double prG(double, short); 
	virtual void CalcAll(const int, const int, class ModelParam *, double **, double **, double **, double **, int, int);

	virtual void compute_ps(const int, const int, ModelParam *, const int);
	virtual void compute_ps(const int, const int, ModelParam * pMD, const int, const int, vector<int>);
	virtual void compute_pk(const int, const int, ModelParam *, const int);
	virtual void compute_pk(const int, const int, ModelParam * pMD, const int, const int, vector<int>);
	virtual void clean_after_em(int); 
	virtual void forward_backward_block(const int, const int, class ModelParam *, double **, double **, double **, double **, int, int, vector<int>);	
};

#endif 
