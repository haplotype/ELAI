#ifndef __DIPLOID_H__
#define __DIPLOID_H__

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "fp.h"
#include "indiv.h"
#include "model.h"

using namespace std;


class DipInd : public Individual
{
private:
    double ** 	m_phi;  		// forward probability, for each invidual, this is an MxKxK matrix; 
	double ** 	m_psi; 		// backward probability.  
	double * hapadmix; 

	double * m_maf; 
	double * m_admix1; 
	vector<vector<double > > vvfwdp; 
	vector<vector<double > > vvfwd; 
	vector<vector<double > > vvbwd; 
//	class Haploid * pHap1, pHap2; 
public:
	DipInd();
	DipInd(int n);
	~DipInd();

	double prG(double, double, short); 
	double prG1(double, double); 
	double prG1(double, double, short); 
	double prG2(double, short); 
	double prG2(double, double, short); 

	void phase(const int, const int, class ModelParam *, double **, double **, double **, double **, int);	
	void impute(const int, const int, class ModelParam *, double **, double **, double **, double **, int);	
	void forward_backward_diploid(const int, const int, class ModelParam *, double **, double **, double **, double **, int);	
	void forward_backward_haploid(const int, const int, class ModelParam *, double **, double **, double **, double **, int);	
//	void forward_backward_block(const int, const int, class ModelParam *, double **, double **, double **, double **, int, int, vector<int>);	
	void forward_backward_block1(const int, const int, class ModelParam *, int, int, vector<int> );	
	void forward_backward_block2(const int, const int, class ModelParam *, int, int, vector<int> );	

	virtual void CalcAll(const int, const int, class ModelParam *, double **, double **, double **, double **, int, int);	

	virtual void compute_ps(const int, const int, ModelParam * pMD, const int);
	virtual void compute_ps(const int, const int, ModelParam * pMD, const int, const int, vector<int>);
	virtual void compute_pk(const int, const int, ModelParam * pMD, const int);
	virtual void compute_pk(const int, const int, ModelParam * pMD, const int, const int, vector<int>);
	virtual void clean_after_em(int); 
	virtual void forward_backward_block(const int, const int, class ModelParam *, double **, double **, double **, double **, int, int, vector<int>);	
};

#endif
