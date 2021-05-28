#ifndef __PARAM_H__
#define __PARAM_H__

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "diploid.h"
#include "haploid.h"
#include "fpmath.h"
#include <iostream>
using namespace std;
		 
//differenet populations share theta while having different alpha and r;
class ModelParam
{                                    
private:
	int nB; 
	double * ru; 
	double * rs; 
	double ** eta;  			
	double ** alpha; 
	double * rk; 
	double ** theta;  			
	double *** beta; 
	double ** kappa;  			

public:
	ModelParam();
	~ModelParam(); 
	void Init(int, int, int, double, double, double, int); 
	inline void Setru(int m, double s) {ru[m] = s;}
	inline void Setrs(int m, double s) {rs[m] = s;}
	inline void Seteta(int m, int k, double s) {eta[m][k] = s;}
	inline void Setalpha(int m, int k, double a) {alpha[m][k] = a;}
	inline void Setkappa(int m, int k, double a) {kappa[m][k] = a;}

	inline void Setrk(int m, double s) {rk[m] = s;}
	inline void Settheta(int m, int k, double s) {theta[m][k] = s;}
	inline void Setbeta(int m, int s, int k, double a) {beta[m][s][k] = a;}
	//set is used to update parameters, one by one; 
	//while get functions return pointers;
	inline double * Getru(void) {return ru;}
	inline double * Getrs(void) {return rs;}
	inline double ** Geteta(void) {return eta;}            
	inline double ** Getalpha(void) {return alpha;}
	inline double ** Getkappa(void) {return kappa;}

	inline double * Getrk(void) {return rk;}
	inline double ** Gettheta(void) {return theta;}            
	inline double *** Getbeta(void) {return beta;}

	inline int GetnB(void) {return nB;}
};

#endif
