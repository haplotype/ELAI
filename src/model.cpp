#include <stdio.h>
#include <stdlib.h> 
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include "model.h"
#include "anchap.h"
#include "fpmath.h"
#include "control.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_permutation.h"
#include "gsl/gsl_permute.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_statistics.h"
extern "C"{
#include "gsl/gsl_randist.h"
}
#include "gsl/gsl_rng.h"
#include "gsl/gsl_eigen.h"

using namespace std;

double global_Ne; 

void ModelnData::print_progress_bar(int last, char* str, int p, int total)
{
	if(m_silence == 1) return; 
	int progress = (int) (100.0 * p / total); 
	int barsize = (int) (progress / 2.0); 
	char bar[100];
	memset(bar, '\0', 100); 
	if(!last) {
		for (int i = 0; i < barsize; i++)
			bar[i] = '>'; 
		printf("%s [%-50s] %d%%\r", str, bar, progress); 
		fflush(stdout); 	
	} else {
		for (int i = 0; i < barsize; i++)
			bar[i] = '>'; 
		printf("%s [%-50s] 100%%\n", str, bar); 
	}
}

ModelnData::ModelnData(void)
{
	morgan = NULL; 
	snpmgt = NULL; 
	snpIn01 = NULL; 
	m_silence = 0; 
	m_exclude_maf = -0.0; 
	m_exclude_miss = 1.00; 
	m_exclude_nopos = 1; 
	m_exclude_miss1 = 0; 
	m_allele_coding_mode = 0;  //0 minor is reference; 1 major is referece; genotype is count of reference alleles. 
	m_num = 10;      //this is beta distribution to link between top and bottom layer. 
	nK = 10; 
	nS = 2; 
	randSeed = 0;
	nEMRuns = 1;
	nMaxSteps = 20;
	nWarmSteps = 0; 
	vPARCrs.clear(); 
	vPARCpos.clear(); 
	
	global_Ne = 1; 
	nLoci = 0; 
	
	fnOutput.assign("\0");
	fnRsPos.assign("\0");
	
	nGeneFlank = 0; 
	nDip = 0;
	nHap = 0;
	pMP = NULL;

	m_npath = 2; 

	m_upper_switches = 1; 
	m_lower_switches = 1; 
	m_rrt = 1; 
	m_morgan = -1; 
	m_mixgen = 20; 

	m_verbose = 1; 
	m_impute = 0; 
	m_init = 2; 

	hpa = 0; 
	psa = 0; 
	lai = 0; 
	m_mu = 0.0;
}

ModelnData::~ModelnData()
{
	if(pMP) {delete[] pMP;  pMP = NULL;} 
	vGin.clear();
	vPin.clear(); 
	vFileIndiv.clear(); 
	vFilePloid.clear(); 
	vsRsnum.clear(); 
	nLoci = 0; 
	
	nDip = 0;
	nHap = 0;
}

void ModelnData::InitModelParam(void)
{
	if (pMP == NULL)
		pMP = new ModelParam[nEMRuns];
	
	//li and stephens used 4Nc = 20, 200. 
	//the unit is centi-morgan per mega-basepair. 
	//we assume 1 centi-morgan per mega-basepair.
	//combine together we have 50/nK cross-over events per mega-basepair; 
	//then divide by nK because the branch length is 1/nK. 
    //m_morgan has unit 100 Mpb, so we have 5000 / nK. 

//    global_Ne = m_num; 

	m_upper_switches = m_morgan * m_mixgen; 
	m_lower_switches = m_morgan * 1000.0;    //default can be thought as mixing 1000 generations agao. 
	double rhok = m_lower_switches / nLoci;  
	double rhos = m_upper_switches / nLoci; 
	double rhou = m_rrt / nLoci; 
if(rhok > 0.999) rhok = 0.999; 
if(rhos > 0.999) rhos = 0.999; 
if(rhou > 0.999) rhou = 0.999; 
//cout << "rhok = " << rhok << endl; 
//cout << "rhos = " << rhos << endl; 
//cout << "rhou = " << rhou << endl; 
	fplog << " ### estimated total genetic distance " << m_morgan << endl; 
	fplog << " ### constrained upper layer switches " << m_upper_switches << endl; 
	fplog << " ### constrained lower layer switches " << m_lower_switches << endl; 
	fplog << " ### constrained ancillary switches " << rhou * nLoci << endl; 

    for (int i = 0; i < nEMRuns; i++)
		pMP[i].Init(nLoci,nS, nK, rhok, rhos, rhou, m_init);
}



void ModelnData::assign_pop_label(void)
{               
    for (int i = 0; i < nIndiv; i++)
	{
		int ph = m_phval.at(i); 
		if(ph < 10) continue; 
	
		if(ph - 10 >= nS) {
			fplog << " ### -p <ph>, where ph < 10+nS "; 
			exit(0); 
		}
	   	pIndiv[i]->set_admix(ph-10, nS); 
	}
}

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

struct params
{
	double b; 
	double Ne; 
};

double digamma(double x)
{
	double x1=1.0/x; 
	double x2=x1/x; 
	double x3=x2/x;
	double x4=x3/x;
	double x5=x4/x;

   	double res = x1 + 0.5*x2 + 5.0/24.0 * x3 + 1.0/16.0 * x4 + 47.0 / 48.0 / 120.0 * x5; 
	return (-log(res)); 
}

double f (double x, void * par) 
{               
	struct params * p = (struct params * ) par; 
	double b =  p->b; 
	double Ne =  p->Ne; 
	double nx = Ne * x; 
	double res = digamma(nx + 6.0) - digamma(6.0 + Ne - nx); 
	res -= 1.0 / nx + 1.0 / (nx + 1.0) + 1.0 /(nx + 2.0) + 1.0 / (nx + 3.0) + 1.0 / (nx + 4.0) + 1.0 / (nx + 5.0); 
	res += 1.0 / (Ne - nx) + 1.0 / (Ne - nx + 1.0) + 1.0 / (Ne - nx + 2.0) + 1.0 / (Ne - nx + 3.0) + 1.0 / (Ne -nx + 4.0) + 1.0 / (Ne - nx + 5.0); 
	res -= b; 
	return(res); 
}

void ModelnData::EM_joint(int mode)
{
	int totalsteps = 0; 
	if(mode == 1) totalsteps = nWarmSteps; 
	else if(mode == 2) totalsteps = nMaxSteps; 
	if(totalsteps == 0) return; 
	//mode == 1; iterative marginalization;
	//mode == 2; quadratic;

	int nsk = nS * nK; 

	double ** sum_jmsk = Allocate2DMatrix(nLoci, nsk); 
	double ** sum_jms = Allocate2DMatrix(nLoci, nS); 
	double ** sum_top = Allocate2DMatrix(nLoci, nK); 
	double ** sum_bot = Allocate2DMatrix(nLoci, nK); 

	double ** sum_lod = Allocate2DMatrix(nLoci, nS); 
	double ** sum_lps = Allocate2DMatrix(nLoci, nS); 

	double * sum_ujs = new double[nLoci]; 

	for (int i = 0; i < nIndiv; i++)
		pIndiv[i]->clean_after_em(1); 
//	for (int i = 0; i < nIndiv; i++)
//		cout << pIndiv[i]->get_weight() << " "; 
//	cout << endl; 

	for (int e = 0; e < nEMRuns; e++)
	{
		double lastlike = 0; 

  		assign_pop_label(); 
		for (int steps = 0; steps < totalsteps; steps++)
		{
			double totalloglike = 0; 
			
			for (int m = 0; m < nLoci; m++)              
			{
				sum_ujs[m] = 0.0; 
				for (int j = 0; j < nsk; j++)
					sum_jmsk[m][j] = 0; 
				for (int s = 0; s < nS; s++)
					sum_jms[m][s] = 0; 
				for (int j = 0; j < nK; j++)
					sum_top[m][j] = sum_bot[m][j] = 0.0; 
				for (int s = 0; s < nS; s++)
					sum_lod[m][s] = sum_lps[m][s] = 0; 
			}

			///////////////////////////////////////////////
			double ** th = (pMP+e)->Gettheta(); 
			for (int i = 0; i < nK; i++)
			{
				for (int m = 0; m < nLoci; m++)
				{       if(th[m][i] > 0.999) th[m][i] = 0.999; 
                                        else if(th[m][i] < 0.001) th[m][i] = 0.001; 
					pHap[i]->set_af(m, th[m][i]); 
                                }
				pHap[i]->CalcAll(nLoci, nS, pMP+e, sum_top, sum_bot, sum_lod, sum_lps, i, 1); 
				double * js = pHap[i]->get_pjs(); 
				for (int m = 0; m < nLoci; m++)
					sum_ujs[m] += js[m]; 
				pHap[i]->clean_after_em(1); 
			}
			/////////////////////////////////////////////

			for (int m = 1; m < nLoci; m++)
			{
				double tru = sum_ujs[m] / nK;
				if(tru > 1) tru = 1; 
				(pMP+e)->Setru(m, tru); 
			}
			//update eta; 
			/////////////////////////////////////////////////////////////////////
			if(fabs(global_Ne - 1) > 1e-6) 
			{
				int iter = 0, max_iter = 100;
				const gsl_root_fsolver_type *T;
				gsl_root_fsolver *solver;
				gsl_function FDF;
				struct params para;
				para.Ne = global_Ne; 

				FDF.function = &f;
	//			FDF.df = &df;
	//			FDF.fdf = &fdf;
				FDF.params = &para;

				T = gsl_root_fsolver_brent;
				solver = gsl_root_fsolver_alloc (T);

				for (int m = 0; m < nLoci; m++)
				{
					for (int s = 0; s < nS; s++)
					{
						double x = 0.5;
						double ss = sum_lod[m][s] / sum_lps[m][s]; 
						if(ss > 1e10) 
							x = 0.999; 
						else 
						{
							double x_lo = 0.001, x_hi = 0.999; 
							para.b = ss; 
							int status;
							iter = 0;
							gsl_root_fsolver_set (solver, &FDF, x_lo, x_hi);
							do {
								iter++;
								status = gsl_root_fsolver_iterate (solver);
								x = gsl_root_fsolver_root (solver);
								x_lo = gsl_root_fsolver_x_lower (solver);
								x_hi = gsl_root_fsolver_x_upper (solver);
								status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);
							 }
							while (status == GSL_CONTINUE && iter < max_iter);
						}
						(pMP+e)->Seteta(m, s, x); 
					}   //eta; 
				}
				gsl_root_fsolver_free (solver);
			}
			else
			{
				for (int m = 0; m < nLoci; m++)
					for (int s = 0; s < nS; s++)
					{
						double tt = 0.5; 
						double lod = sum_lod[m][s] / sum_lps[m][s]; 
						if(fabs(lod) > 1e-20) 
						{
							lod /= -3.1415926;  
							tt = 0.5 - atan(lod) / 3.1415925;   
							if(tt > 0.999) tt = 0.999; 
							else if(tt < 0.001) tt = 0.001; 
						}
						(pMP+e)->Seteta(m, s, tt); 
					}   //eta; 
			}
			///////////////////////////////////////////////////////////

			assign_pop_label();                          
		    fplog << "### "; 
			for (int i = 0; i < nIndiv; i++)
			{
				if(pIndiv[i]->get_mask() == 1) continue;
				if(mode == 1) //-w 
					pIndiv[i]->CalcAll(nS, nK, pMP+e, sum_top, sum_bot, sum_jmsk, sum_jms, 1, 0); 
				else          //-s
					pIndiv[i]->CalcAll(nS, nK, pMP+e, sum_top, sum_bot, sum_jmsk, sum_jms, 2, 0); 
				totalloglike += pIndiv[i]->GetLikelihood(); 
				pIndiv[i]->clean_after_em(1); 
			}

			//update theta; 
//			for (int m = 0; m < nLoci; m++)
//			{
//				for (int k = 0; k < nK; k++)
//				{                      
//					double tt = 0.5; 
//					double lod = sum_top[m][k] / sum_bot[m][k]; 
//					if(fabs(lod) > 1e-20) 
//					{
//						lod /= -3.1415926;  
//						tt = 0.5 - atan(lod) / 3.1415925;   
//						if(tt > 0.999) tt = 0.999; 
//						else if(tt < 0.001) tt = 0.001; 
//					}
//					//bound theta; 
//					if(isnan(tt)) {
//						cout << "tt is nan in theta unpdate " << endl; 
//				   		(pMP+e)->Settheta(m, k, 0.5); 
//					}
//					else
//				   		(pMP+e)->Settheta(m, k, tt); 
//				}   //share the lower clusters; 
//			}
			for (int m = 0; m < nLoci; m++)
			{
				for (int k = 0; k < nK; k++)
				{                      
					double tt = (sum_top[m][k]) / (sum_bot[m][k]); 
					if(tt > 0.999) tt = 0.999; 
					else if(tt < 0.001) tt = 0.001; 
					//bound theta; 
//					if(isnan(tt)) {
//						fplog << "tt is nan in theta unpdate " << m << " " << sum_top[m][k] << " " << sum_bot[m][k] << endl; 
//						cout << "tt is nan in theta unpdate " << m << " " << sum_top[m][k] << " " << sum_bot[m][k] << endl; 
////						exit(0); 
//				   		(pMP+e)->Settheta(m, k, 0.5); 
//					}
//					else
				   		(pMP+e)->Settheta(m, k, tt); 
				}   //share the lower clusters; 
			}

			double ploids = 0;   //weighted ploidy
			for (int i = 0; i < nIndiv; i++)
			{
				if(pIndiv[i]->get_mask() == 1) continue;
			   	ploids += pIndiv[i]->Getploid() * pIndiv[i]->get_weight(); 
			}

			//update beta; 
			for (int m = 0; m < nLoci; m++)
	   		{                 
				for (int s = 0; s < nS; s++)
				{
					double total = 0; 
					for (int j = s * nK, k = 0; k < nK; k++, j++)
						total += sum_jmsk[m][j]; 
//					if(total < 1e-100) { cout << "total tiny " << m << " " << s << " " << total << endl; continue;} 
					for (int j = s * nK, k = 0; k < nK; k++, j++)
					{
						double ta = (sum_jmsk[m][j])/ total;
//						if(isnan(ta)) { 
//                            fplog << "ta is nan in beta update " << endl; 
//                            cout << "ta is nan in beta update " << endl; 
//							(pMP+e)->Setbeta(m, s, k, 1.0 / nK);
//						}
//						else
							(pMP+e)->Setbeta(m, s, k, ta);
					}
				}

				double jk = 0; 
				for (int j = 0; j < nsk; j++)
					jk += sum_jmsk[m][j]; 
				double trk = jk / ploids; //coupling; 
//				if(isnan(trk)) {
//					fplog << "trk is nan in rk upate " << endl; 
//					cout << "trk is nan in rk upate " << endl; 
//					(pMP+e)->Setrk(m, 1e-20);
//				}
//				else if(trk > 1) { //cout << "trk too big " << trk << " " << m << " " << ploids << endl; 
//					trk = 1;}
//				else 
				if(trk > 0.999) trk = 0.999; 
					(pMP+e)->Setrk(m, trk); 
				//set rk; 

				double js = 0; 
				for (int s = 0; s < nS; s++)
					js += sum_jms[m][s]; 
				double trs = js / ploids; //coupling; 
//				if(isnan(trs)) {           
//					fplog << "trs is nan in rs upate " << endl; 
//					cout << "trs is nan in rs upate " << endl; 
//					(pMP+e)->Setrs(m, 1e-20);
//				}
//				else if(trs > 1) { cout << "trs too big " << trs << " " << m << endl; trs = 1;}
//				else 
				if(trs > 1) trs = 1; 
					(pMP+e)->Setrs(m, trs); 
				//set rs; 
			}
			(pMP+e)->Setrs(0, 1); 
			(pMP+e)->Setrk(0, 0); 
			
			double ktotal = 0; 
			double * rk = (pMP+e)->Getrk(); 
			for (int m = 1; m < nLoci; m++)
   				ktotal += rk[m]; 
//			if(ktotal > m_lower_switches) 
//			{
//				for (int m = 1; m < nLoci; m++)
//				{
//					double tr = rk[m]; 
//					(pMP+e)->Setrk(m, 1-exp(log(1-tr) / ktotal * m_lower_switches)); 
//				}
//			}
			
			double utotal = 0; 
			double * ru = (pMP+e)->Getru();
			for (int m = 1; m < nLoci; m++)
				utotal += ru[m]; 

			if(utotal > m_rrt) 
			{
				for (int m = 1; m < nLoci; m++)
				{
					double tr = ru[m]; 
					double tt = 1-exp(log(1-tr) / utotal * m_rrt); 
					(pMP+e)->Setru(m, tt); 
				}
			}

			double stotal = 0; 
			double * rs = (pMP+e)->Getrs();
			for (int m = 1; m < nLoci; m++)
				stotal += rs[m]; 
			
			if(stotal > m_upper_switches) 
			{
				for (int m = 1; m < nLoci; m++)
				{
					double tr = rs[m]; 
					double tt = 1-exp(log(1-tr) / stotal * m_upper_switches);
					(pMP+e)->Setrs(m, tt); 
				}
			}

			///////////////////////////////////////////////
			double likeinc = totalloglike - lastlike; 
			lastlike = totalloglike; 

			char buf[100]; 
//			sprintf(buf, "%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f", totalloglike, likeinc, ktotal, stotal, utotal);    
			sprintf(buf, "%8.3f\t%8.3f\t", totalloglike, likeinc);    
			if(mode == 1)  cout << " ### Approximate update";   
			else cout << " ### Exact update "; 
			cout << e << " " << steps << "\t " << buf << endl; 

			if(mode == 1)  fplog << " ### Approximate update";   
			else fplog << " ### Exact upate "; 
			fplog << e << " " << steps << "\t " << buf << endl; 
		}

	}

	Free2DMatrix(sum_top); 
	Free2DMatrix(sum_bot); 
	Free2DMatrix(sum_jmsk); 
	Free2DMatrix(sum_jms); 
	Free2DMatrix(sum_lod); 
	Free2DMatrix(sum_lps); 
	delete[] sum_ujs; 
}                   


void ModelnData::compute_ps(int out_mode, int comp_mode)
{
	comp_mode = 2; 
	int unit = 1000; 
	int nfile = ceil(nLoci / (double) unit); 
	vector<int> pbreak; 
	for (int i = 0; i < nfile; i++)
		pbreak.push_back(i * unit); 
	pbreak.push_back(nLoci); 

//	cout << nfile << " --- "; 
//	for (unsigned i = 0; i < pbreak.size(); i++)
//		cout << pbreak.at(i) << " " ; 
//	cout << endl; 

	fstream outfile; 
	string fn("output/");
	fn.append(fnOutput);
	char buf[100];
	sprintf(buf, ".ps%d%d.txt", comp_mode, out_mode); 
	fn.append(buf);

	outfile.open(fn.c_str(), ios::out);
	if(!outfile.good()) 
	{
		fplog << "ERROR: failed to open file to write." << endl;
		cout << "-elai: failed to open file to write" << endl;
		return; 
	}

	if(unit > nLoci) 
	{
		for (int i = 0; i < nIndiv; i++)
		{
			if(pIndiv[i]->GetisPanel()) continue; 
			if(pIndiv[i]->get_mask() == 1) continue; 
			pIndiv[i]->compute_ps(nS, nK, pMP, comp_mode); 
			double ** ps = pIndiv[i]->get_ps(); 

			char buf[100]; 
			if(pIndiv[i]->Getploid() == 1) 
			{
				for (int m = 0; m < nLoci; m++)
				{
					for (int s = 0; s < nS; s++)
					{
						double	temp = ps[m][s]; 
						sprintf(buf, "%5.3f ", temp); 
						outfile << buf; 
					}
				}
			}
			else {
				for (int m = 0; m < nLoci; m++)
				{
	//					for (int s = 0; s < nS; s++)
	//					{
	//						sprintf(buf, "%5.3f ", ps[m][s]); 
	//						outfile << buf; 
	//					}
	//					outfile << " "; 
					if(out_mode == 1) 
					{
						for (int s = 0; s < nS; s++)
						{
							double temp = 0; 
							for (int s1 = 0; s1 < nS; s1++)
								temp += 2*ps[m][nS*s+s1]; 
							sprintf(buf, "%5.3f ", temp); 
							outfile << buf; 
						}
					}
					else //out_mode == 2; 
					{
						for (int s = 0; s < nS; s++)
							for (int s1 = 0; s1 < nS; s1++)
							{
								double temp = ps[m][nS*s+s1]; 
								sprintf(buf, "%5.3f ", temp); 
								outfile << buf; 
							}
					}
				}
			}
			outfile << endl; 

			pIndiv[i]->clean_after_em(1); 
		}
	}
	else 
	{
		for (int i = 0; i < nIndiv; i++)
		{
			if(pIndiv[i]->GetisPanel()) continue; 
			if(pIndiv[i]->get_mask() == 1) continue; 
			for (int b = 0; b < nfile; b++)
			{
				pIndiv[i]->compute_ps(nS, nK, pMP, comp_mode, b, pbreak); 
				double ** ps = pIndiv[i]->get_ps(); 

				char buf[100]; 
				if(pIndiv[i]->Getploid()==1) 
				{
					for (int p = 0; p < pbreak.at(b+1)-pbreak.at(b); p++)
					{
						for (int s = 0; s < nS; s++)
						{
							double temp = ps[p][s]; 
							sprintf(buf, "%5.3f ", temp); 
							outfile << buf; 
						}
					}
				}
				else {
					for (int p = 0; p < pbreak.at(b+1)-pbreak.at(b); p++)
					{
								   
						if(out_mode == 1) 
						{
							for (int s = 0; s < nS; s++)
							{
								double temp = 0; 
								for (int s1 = 0; s1 < nS; s1++)
									temp += 2*ps[p][nS*s+s1]; 
								sprintf(buf, "%5.3f ", temp); 
								outfile << buf; 
							}
						}
						else //out_mode = 2; 
						{
							for (int s = 0; s < nS; s++)
								for (int s1 = 0; s1 < nS; s1++)
								{
									double temp = ps[p][nS*s+s1]; 
									sprintf(buf, "%5.3f ", temp); 
									outfile << buf; 
								}

						}
					}
				}
			}
			outfile << endl; 

			pIndiv[i]->clean_after_em(1); 
		}
	}
	outfile.close(); 
}

void ModelnData::compute_pk(int out_mode, int comp_mode)
{
	fstream outfile; 
	string fn("output/");
	fn.append(fnOutput);
	char buf[100];
	sprintf(buf, ".pk%d.txt", out_mode); 
	fn.append(buf);

	outfile.open(fn.c_str(), ios::out);
	if(!outfile.good()) 
	{
		fplog << "ERROR: failed to open file to write." << endl;
		cout << "-elai: failed to open file to write" << endl;
		return; 
	}

	if(comp_mode == 1) 
	{
		for (int e = 0; e < nEMRuns; e++)
		{
			for (int i = 0; i < nIndiv; i++)
			{
//				cout << i << " --- " << pIndiv[i]->GetisPanel() << " "; 
				if(pIndiv[i]->GetisPanel()) continue; 
				if(pIndiv[i]->get_mask() == 1) continue; 
				pIndiv[i]->compute_pk(nS, nK, pMP+e, 1); 
				double ** pk = pIndiv[i]->get_pk(); 
				if(pk == NULL) continue; 

				char buf[100]; 
				for (int m = 0; m < nLoci; m++)
				{
					for (int k = 0; k < nK; k++)
					{
						sprintf(buf, "%5.3f ", pk[m][k]); 
						outfile << buf; 
					}
					outfile << " "; 
				}
				outfile << endl; 

				pIndiv[i]->clean_after_em(1); 
//				cout << i << endl; 
			}
			outfile << endl; 
		}
	}
	else //compute mode == 2 && diploid; 
	{
		for (int i = 0; i < nIndiv; i++)
		{
			if(pIndiv[i]->GetisPanel()) continue; 
			if(pIndiv[i]->get_mask() == 1) continue; 
			pIndiv[i]->compute_pk(nS, nK, pMP, 2); 
			double ** pk = pIndiv[i]->get_pk(); 

			char buf[100]; 
			for (int m = 0; m < nLoci; m++)
			{
				for (int k = 0; k < nK; k++)
				{
					double mk = pk[m][k]; 
					sprintf(buf, "%5.3f ", mk); 
					outfile << buf; 
				}
				outfile << " "; 
			}
			outfile << endl; 

			pIndiv[i]->clean_after_em(1); 
		}
	}

	outfile.close(); 
}


void ModelnData::get_het(void)
{
	; 
}

void ModelnData::open_log(void)
{
	string sfn("output/");
	sfn.append(fnOutput);
	sfn.append(".log.txt");
	fplog.open(sfn.c_str(), ios::out);
	if(!fplog.is_open()) 
	{
		cout << "-elai: cannot open log file" << endl;
		safe_exit(); 
	}
}

void ModelnData::close_log(void)
{
	fplog.close();
}

void ModelnData::play_ground(void)
{
	cout << "in play ground" << endl; 
 

//	int abeg[5] = {1599, 2908, 2812, 3871, 3900};
//	int aend[5] = {1601, 2911, 2813, 3872, 3901}; 
//	string hla[5] = {"a","b","c","drb1","dqb1"};
	

	nDip = (int) m_num; 
	nLoci = nEMRuns; 
	int beg = nWarmSteps; 
	int end = nMaxSteps;
//	cout << nLoci << endl; 
//	cout << beg << endl; 
//	cout << end << endl; 
	double ** ps = Allocate2DMatrix(nDip, nS*nLoci); 
//	double ** dist = Allocate2DMatrix(nDip, nDip); 

	string fn; 
	fn.assign(fnFILE); 
	ifstream infile; 
	streambuf * pbuf;
	char delimit[] = ";, \t";
	infile.open(fn.c_str(), ios::in);
	if(!infile.is_open()) 
	{
		cout << "-elai: cannot open ps1 file: " << endl; 
		cout << fn << endl; 
		exit(0); 
	} 
	pbuf = infile.rdbuf();
	// open file; 

	string line; 
	line.assign(getline(pbuf)); 
		
	int ni = 0; 
	while(!line.empty())
	{   
		char *res = strtok((char*)line.c_str(), delimit); 
		if(res == NULL) break; 
		int pos = 0; 
		ps[ni][pos++] = atof(res); 
		while(1)
		{          
			string tmp;
			res = strtok(NULL, delimit); 
			if(res == NULL) break; 
			ps[ni][pos++] = atof(res); 
		}   // snp summary in a2c;    //A C G T ? + -
		line.assign(getline(pbuf));
		ni++; 
	}

//	for (int r = 0; r < 1; r++)
	{

//		for (int i = 0; i < nDip; i++)
//		{
//			for (int j = 0; j <= i; j++)
//			{
//				double dd = 0; 
//				for (int m = beg-1; m < end; m++)
//				{       
//					double temp = 0.0; 
//					for(int s = 0; s < nS; s++)
//						temp += ps[i][m*nS+s] * ps[j][m*nS+s];
//					dd += temp; 
//				}
//				dist[i][j] = dist[j][i] = dd / 4.0 / (1.0 + end - beg); 
//			}
//		}

		fstream outfile; 
		fn.assign(fnDOC); 
//		fn.append(hla[r]); 
		fn.append(".txt"); 
		cout << " write file name " << fn << endl; 

		outfile.open(fn.c_str(), ios::out);
		if(!outfile.good()) 
		{
			cout << "-elai: failed to open file to write" << endl;
			exit(0);  
		}
		for (int i = 0; i < nDip; i++)
		{       
//			for (int j = 0; j < nDip; j++)
			for (int j = beg-1; j < end; j++)
			{
				char buf[100]; 
//				sprintf(buf, "%5.4f ", dist[i][j]); 
				for (int s =0; s < nS; s++)
				{
					sprintf(buf, "%5.4f ", ps[i][j*nS+s]); 
					outfile << buf; 
				}
			}
			outfile << endl; 
		}
		outfile.close(); 
	}

	cout << "finishing writing in play ground" << endl; 

}

void ModelnData::simulate_admixed(void)
{
//	{
//	   cout << "in maf-diff block " << endl;  
//		string fn; 
//		fstream outfile; 
//
//		fn.append(fnOutput);
//		fn.append(".mafdiff.txt");
//		outfile.open(fn.c_str(), ios::out);
//		if(!outfile.good()) 
//		{
//			fplog << "ERROR: failed to open file to write." << endl;
//			cout << "-elai: failed to open file to write" << endl;
//			return; 
//		}
//
//		outfile << "###maf-diff.txt" << endl; 
//		for (int m = 0; m < nLoci; m++)
//		{
//			double af1 = 0; 
//			for (int i = 0; i < 90; i++)
//			{
//				short gt = pIndiv[i]->GetsnpGT(m);
//				af1 += gt; 
//			}
//			af1 /= 90.0; 
//
//			double af2 = 0; 
//			for (int i = 91; i < 180; i++)
//			{
//				short gt = pIndiv[i]->GetsnpGT(m);
//				af2 += gt; 
//			}
//			af2 /= 90.0; 
//
//			if(fabs(af1 - af2) > 0.001) 
//			{
//				outfile << af1 << " \t" << af2 << endl; 
////				cout << af1 << " \t" << af2 << endl; 
//			}
//		}
//		outfile.close(); 
//	}
//	return; 

	short ** snpgt = Allocate2DShortMatrix(nIndiv, nLoci); 
    for (int i = 0; i < nIndiv; i++)
	{
		for (int m = 0; m < nLoci; m++)
		{
			short gt = pIndiv[i]->GetsnpGT(m); 
			snpgt[i][m] = gt; 
		}
	}

	int nsnp = nMaxSteps;  
	cout << "nsnp = " << nsnp << endl; 
	cout << "nloci = " << nLoci << endl; 

	short ** c1 = Allocate2DShortMatrix(60, nLoci); 
	short ** c2 = Allocate2DShortMatrix(60, nLoci); 

	for (int k = 0; k < 10; k++)
	{
		short * pf = &snpgt[120-k*2-1][0]; 
		short * pm = &snpgt[240-k*2-1][0]; 
		short * pc = &snpgt[420-k*2-1][0]; 
		int cur = 0; 
		int rand = 0; 
		while (cur < nLoci)
		{
//			int len = gsl_ran_exponential(gsl_r, (double)nsnp); 
			int len = nsnp; 
//			double rand =  gsl_rng_uniform(gsl_r);
			rand++; 
//			if( rand < 0.47) 
			if(rand % 3 == 0) 
			{
				for (int m = cur; m < min(nLoci, cur+len); m++)
				{
					c1[k][m] = pf[m]; 
					fplog << " 0";
				}
			}
//			else if(rand < 0.95) 
			else if (rand % 3 == 1) 
			{
				for (int m = cur; m < min(nLoci, cur+len); m++)
				{
					c1[k][m] = pm[m]; 
					fplog << " 1"; 
				}
			}
			else
			{
				for (int m = cur; m < min(nLoci, cur+len); m++)
				{
					c1[k][m] = pc[m]; 
					fplog << " 2"; 
				}

			}
			cur += len; 
		}
		fplog << endl; 

		pf = &snpgt[120-k*2-2][0]; 
		pm = &snpgt[240-k*2-2][0]; 
		pc = &snpgt[420-k*2-2][0]; 
		cur = 0; 
//		int rand = 1; 
		while (cur < nLoci)
		{
			int len = gsl_ran_exponential(gsl_r, (double)nsnp); 
//			double rand =  gsl_rng_uniform(gsl_r);
			rand++; 
//			if( rand < 0.6) 
			if(rand % 3 == 0) 
			{
				for (int m = cur; m < min(nLoci, cur+len); m++)
				{
					c2[k][m] = pf[m]; 
					fplog << " 0";
				}
			}
//			else if(rand < 0.8) 
			if(rand % 3 == 1) 
			{
				for (int m = cur; m < min(nLoci, cur+len); m++)
				{
					c2[k][m] = pm[m]; 
					fplog << " 1"; 
				}
			}
			else
			{
				for (int m = cur; m < min(nLoci, cur+len); m++)
				{
					c2[k][m] = pc[m]; 
					fplog << " 2"; 
				}

			}
			cur += len; 
		}
		fplog << endl; 
	}

    cout << "i am here " << endl; 

	string fn; 
	fstream outfile; 

	fn.append(fnOutput);
	fn.append(".admixed.txt");
	outfile.open(fn.c_str(), ios::out);
	if(!outfile.good()) 
	{
		fplog << "ERROR: failed to open file to write." << endl;
		cout << "-elai: failed to open file to write" << endl;
		return; 
	}
	
	outfile << nDip << endl; 
	outfile << nLoci << endl; 
	for (int m = 0; m < nLoci; m++)
	{
		string rs(vsRsnum.at(m)); 
		outfile << rs << " "; 
//		for (int i = 0; i < 50; i++)
//		{
//			short gt = pIndiv[2*i]->GetsnpGT(m);
//			gt += pIndiv[2*i+1]->GetsnpGT(m);
//			if(gt == 0)  outfile << " " << mapRs2mm[rs].first << mapRs2mm[rs].first; 
//			else if(gt == 1)  outfile << " " << mapRs2mm[rs].first << mapRs2mm[rs].second; 
//			else  outfile << " " << mapRs2mm[rs].second << mapRs2mm[rs].second; 
//		}
//		for (int i = 0; i < 50; i++)
//		{
//			short gt = pIndiv[120+2*i]->GetsnpGT(m);
//			gt += pIndiv[120+2*i+1]->GetsnpGT(m);
//			if(gt == 0)  outfile << " " << mapRs2mm[rs].first << mapRs2mm[rs].first; 
//			else if(gt == 1)  outfile << " " << mapRs2mm[rs].first << mapRs2mm[rs].second; 
//			else  outfile << " " << mapRs2mm[rs].second << mapRs2mm[rs].second; 
//		}

        for(int k = 0; k < 10; k++)
		{
			short gt = c1[k][m] + c2[k][m]; 
			if(gt == 0)  outfile << " " << mapRs2mm[rs].first << mapRs2mm[rs].first; 
			else if(gt == 1)  outfile << " " << mapRs2mm[rs].first << mapRs2mm[rs].second; 
			else  outfile << " " << mapRs2mm[rs].second << mapRs2mm[rs].second; 
		}
		outfile << endl; 
	}
	outfile.close(); 

	delete[] c1; 
	delete[] c2; 
}

void ModelnData::compute_lhs(int comp_mode)  //this compute pairwise lhs. 
{
	int ni = 0; 
	for (int i = 0; i < nIndiv; i++)
	{
		if(pIndiv[i]->get_mask() == 1) continue; 
		ni++; 
	}

	double ** lhs = Allocate2DMatrix(ni, ni); 
	for (int i = 0; i < ni; i++)
		for (int j = 0; j < ni; j++)
			lhs[i][j] = 0; 

	cout << endl << " num loci = " << nLoci; 
	int unit = 1000; 
	if(nLoci > unit) 
	{ 
		int nfile = ceil(nLoci / (double) unit); 
		vector<int> pbreak; 
		for (int i = 0; i < nfile; i++)
			pbreak.push_back(i * unit); 
		pbreak.push_back(nLoci); 
		cout << nfile << " --- "; 
		for (unsigned i = 0; i < pbreak.size(); i++)
			cout << pbreak.at(i) << " " ; 
		cout << endl; 
		//bpreak has 11 elements; 

		double ** gt = Allocate2DMatrix(ni, unit*nK); 
		for (int e = 0; e < nEMRuns; e++)
		{
			for (int f = 0; f < nfile; f++)
			{
				cout << f << " of " << nfile << " many blocks" << endl; 
				int ii = 0; 
				for (int i = 0; i < nIndiv; i++)
				{
					if(pIndiv[i]->get_mask() == 1) continue; 
					pIndiv[i]->compute_pk(nS, nK, (pMP+e), comp_mode, f, pbreak); 
					double ** pk = pIndiv[i]->get_pk(); 
   					for (int s = 0; s < pbreak[f+1] - pbreak[f]; s++)
					{
						for (int k = 0; k < nK; k++)
							gt[ii][s*nK+k] = pk[s][k]; 
					}
					ii++;
					pIndiv[i]->clean_after_em(1); 
				}

				for (int i = 0; i < ni; i++)
					for (int j = i; j < ni; j++)
					{
						for (int s = 0; s < pbreak[f+1] - pbreak[f]; s++) 
							for (int k = 0; k < nK; k++)
						    	lhs[i][j] += gt[i][s*nK+k] * gt[j][s*nK+k]; 
					}
			}

		}
		Free2DMatrix(gt); 
		vector<int> ().swap(pbreak); 
	}
	else 
	{
		double ** gt = Allocate2DMatrix(ni, nLoci*nK); 
		for (int e = 0; e < nEMRuns; e++)
		{
			int ii = 0; 
			for (int i = 0; i < nIndiv; i++)
			{
				if(pIndiv[i]->get_mask() == 1) continue; 
				pIndiv[i]->compute_pk(nS, nK, (pMP+e), comp_mode); 
				double ** pk = pIndiv[i]->get_pk(); 
				for (int m = 0; m < nLoci; m++)
				{
					for (int k = 0; k < nK; k++)
						gt[ii][m*nK+k] = pk[m][k]; 
				}
				ii++;
				pIndiv[i]->clean_after_em(1); 
			}

			for (int i = 0; i < ni; i++)
				for (int j = i; j < ni; j++)
					for (int m = 0; m < nLoci; m++) 
						for (int k = 0; k < nK; k++)
							lhs[i][j] += gt[i][m*nK+k] * gt[j][m*nK+k]; 

		}
		Free2DMatrix(gt); 
	}

	for (int i = 0; i < ni; i++)
		for (int j = i; j < ni; j++)
			lhs[i][j] /= nLoci; 

	fstream fo; 
	string fn("output/");
	fn.append(fnOutput);
	fn.append(".lhs.txt");

	fo.open(fn.c_str(), ios::out);
	if(!fo.good()) 
	{
		fplog << "ERROR: failed to open file to write." << endl;
		cout << "-elai: failed to open file to write" << endl;
		return; 
	}
    char buf[100]; 
	for (int i = 0; i < ni; i++)
	{
		for (int j = 0; j < ni; j++)
		{
			double temp = 0; 
			if(j < i) temp = lhs[j][i]; 
			else temp = lhs[i][j]; 
			sprintf(buf, "%5.3f ", temp); 
			fo << buf; 
		}
		fo << endl; 
	}

	fo.close(); 

	Free2DMatrix(lhs); 
}
