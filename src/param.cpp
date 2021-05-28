#include <iostream>
#include <string.h>
#include <time.h>
#include "fpmath.h"
#include "param.h"

using namespace std;
ModelParam::ModelParam()
{
	ru = NULL; 
	rs = NULL; 
	eta = NULL; 
	alpha = NULL; 
	rk = NULL; 
	theta = NULL; 
	beta = NULL; 
	kappa = NULL; 
}

ModelParam::~ModelParam(void)
{
	Free2DMatrix(kappa); kappa = NULL; 
	Free2DMatrix(eta);  eta = NULL;
	Free2DMatrix(alpha); alpha = NULL; 
	Free2DMatrix(theta);  theta = NULL;
	Free3DMatrix(beta); beta = NULL; 
	delete[] rk;  rk = NULL;
	delete[] rs; rs = NULL; 

}

void ModelParam::Init(int nLoci, int nS, int nK, double rhok, double rhos, double rhou, int minit)
{
	ru = new double[nLoci];
	rs = new double[nLoci];
	rk = new double[nLoci];
	
	for (int m = 1; m < nLoci; m++)
	{
		ru[m] = rhou; 
		rs[m] = rhos; 
		rk[m] = rhok; 
	}
	rs[0] = ru[0] = 1.0; 
	rk[0] = 0.0; 

	if(nS == 1) 
	{
		for (int m = 1; m < nLoci; m++)
	   		rs[m] = 0.0; 
	}

	eta = Allocate2DMatrix(nLoci, nS);
	alpha = Allocate2DMatrix(nLoci, nS);
	theta = Allocate2DMatrix(nLoci, nK);
	beta = Allocate3DMatrix(nLoci, nS, nK); 

	int div = nK / nS; 
	int rem = nK % nS; 
	int share = nS + rem; 
	int uniq = (div - 1); 
	nB = uniq + share; 
	kappa = Allocate2DMatrix(nS, nB); 
	for (int s = 0; s < nS; s++)
	{
		int beg = share + s * uniq; 
		int end = share + (s+1) * uniq; 
		int wh = 0; 
		for (int i = 0; i < share; i++)
			kappa[s][wh++] = i; 
		for (int i = beg; i < end; i++)
			kappa[s][wh++] = i; 
	}
//	for (int s = 0; s < nS; s++)
//	{
//		for (int b = 0; b < nB; b++)
//			cout << kappa[s][b] << " "; 
//		cout << endl; 
//	}

	for (int m = 0; m < nLoci; m++)
	{           
//		kappa[m][0] = gsl_ran_beta(gsl_r, 1000, 1000); 
//		kappa[m][1] = 1 - kappa[m][0]; 
//		kappa[m][2] = gsl_ran_beta(gsl_r, 1000, 1000); 
//		kappa[m][3] = 1 - kappa[m][2]; 
//		for (int i = 0; i < 4; i++)
//			kappa[m][i] = 0.5; 

		for (int s = 0; s < nS; s++)
			for (int k = 0; k < nK; k++)
   			beta[m][s][k] = 1.0 / nK; 
//			for (int b = 0; b < nB; b++)
//			{
//				int k = kappa[s][b]; 
//   			beta[m][s][k] = 1.0 / nB; 
//   		}
			
		for (int s = 0; s < nS; s++)
			eta[m][s] = gsl_ran_beta(gsl_r, 1000, 1000); 

		for (int j = 0; j < nK; j++)
			theta[m][j] = gsl_ran_beta(gsl_r, 1000,1000);; 
	}

	double * bs = new double[nS]; 
	for (int s = 0; s < nS; s++)
		bs[s] = 1000.0; 
	double * bk = new double[nK]; 
	for (int k = 0; k < nK; k++)
		bk[k] = 100.0; 
	if(minit == 1)
	{
    	int m = gsl_rng_uniform_int(gsl_r,(long unsigned int)(nLoci*0.6)) + (long unsigned int) (nLoci * 0.2); 
//		cout << m << endl; 
   		for (int s = 0; s < nS; s++)
   			gsl_ran_dirichlet(gsl_r, nK, bk, beta[m][s]); 
	}
	else 
	{
//		for (int m = 0; m < nLoci; m++)
//			gsl_ran_dirichlet(gsl_r, nS, bs, alpha[m]); 

		for (int m = 0; m < nLoci; m++)
			for (int s = 0; s < nS; s++)
				gsl_ran_dirichlet(gsl_r, nK, bk, beta[m][s]); 
	}
	delete[] bk; 
	delete[] bs; 
}

