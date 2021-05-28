#ifndef __MODEL_H__                
#define __MODEL_H__

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>

#include "gio.h" 
#include "indiv.h" 
#include "diploid.h"
#include "haploid.h"
#include "fpmath.h"
#include "param.h"

#if defined (MPI_ENABLED)
#include "mpi.h"
#endif

using namespace std;

extern double global_Ne; 

class ModelnData
{
public: 
	fstream fplog; // log;
	
private:
	short ** snpIn01; 
	double ** snpmgt; 
	
private:
	int   nS;                   //number of upper clusters; 
    int 	nK; 				//number of lower clusters;  
    class 	ModelParam * pMP;   //nEMruns dimension. 
	double m_rrt;  //rate ratio for training data; 
	double m_upper_switches;
	double m_lower_switches;
	double m_mixgen; 
	double m_morgan;           //total morgan; 
	double m_Ne;               // 4 x effective diploid poplation size / nK; 
	int hpa;                    //haplotype-phenotype association; 
	int psa;                    //population structure analysis; 
	int lai;                    //local ancestry inference; 
	int m_npath;
	
	int m_verbose;              // 1 or 3; 
	int m_impute;               // 0 or 1; different param init; 
	int m_init;               // how to initialize the markers;
	vector<int> m_phval; 
/***************************************************************************  
  above for model parameters                                                                    
***************************************************************************/ 
private:
	int randSeed;               //user may specify random seed;
	int nEMRuns;                //how many different em runs;
	int nMaxSteps;              //maximum steps of each em run;
	int nWarmSteps; 
	int bPanelOnly; 
	double total_likelihood; 
	double logLikelihood;      	//loglikelihood of pop. 
	
/***************************************************************************  
  above for EM runs;
***************************************************************************/ 
private: 
	int nDip;    				//number of diploid individuals;
	int nHap;    				//number of haploid individuals; 
	int nIndiv; 
	int nCohort; 				//nIndiv = nCohort + nPanel; 
	int nPanel; 
    int nLoci;                  //number of snp loci, same for everyone;
	vector<string> indivID;     // individual IDs; 
	class Individual ** pIndiv; 
	class AncHap ** pHap; 
	
	fstream logfile; 
	vector<string> vsRsnum; 
	vector<string> vGin;
	vector<string> vPin; 
	vector<string> vPos; 
	string fnOutput;  			//output prefix; 
	string fnRsPos; 
	vector<int> vFileIndiv;       //num of diploids in genotype files; 
	vector<int> vFilePloid;     //ploid of input genotype files. 
	  
    string fnFILE;    //a dummy filename for flexible use; 
	string fnDOC; 
	
	map<string, pair<char, char> > mapRs2mm;  //rs to major, minor alleles.  
	map<string, double> mapRs2maf;              //rs to minor allele freq. 
	map<string, double> mapRs2var;            //rs to var = 2*maf*(1-maf). 
	map<string, long> mapRs2pos; 
	map<string, int> mapRs2chr; 
	map<string, double> mapRs2gm;     //genetic map; 

	double * morgan; 
	double m_mu; //to smooth the output; 0.05 -- 0.10; 

	string fnGene; 
	int nGeneFlank; 

	double m_num; 
	map<string, int> mPanelRsQ;   

private:
	int m_silence; 
	int m_allele_coding_mode; 
	double m_exclude_maf; 
	double m_exclude_miss; 
	int m_exclude_miss1; 
	int m_exclude_nopos; 
	int m_disable_hap; 
public: 
	inline void set_mu(double s){m_mu =s;}
	inline void Setsilence(int s) {m_silence = s;}
	inline void Setverbose(void) {m_verbose = 3;}
	inline void Setexcludemaf(double s) {m_exclude_maf = s;}
	inline void Setexcludemiss(double s) {m_exclude_miss = s;}
	inline void Setexcludemiss1(void) {m_exclude_miss1 = 1;}
	inline void Setexcludenopos(void) {m_exclude_nopos = 1;}
	inline void Setallelecoding(int s){m_allele_coding_mode = s;} 
 	void print_progress_bar(int last, char* str, int p, int total);     
	vector<string> vPARCrs; 
	vector<int> vPARCpos; 
//#if defined (IMPUTATION)
//private:
//	double percentMasked;       // nLoci * this = nMasked;
//	vector<string> vMaskedSNPrs; 
//	int  parcLoci; 
//	int  nMasked; 
//public:
//	inline void Setmask(double s){percentMasked = s;};
//	inline double Getmask(void) {return percentMasked;}; 
//	void MaskSNP(void);					//randomly mask genotype. substitute 0,1's with ?'s.
//	void MaskRand(void); 
//	void ImputeMasked(void);             //missing genotype '?' imputation. 
//#endif
	
public:
	ModelnData();
	~ModelnData();
	
    
	inline void set_disable_hap(void) {m_disable_hap=1;}
	inline void Setnum(double s) {m_num = s;}

	inline void SetfnFILE(string s) {fnFILE.assign(s);}
	inline void SetfnDOC(string s) {fnDOC.assign(s);}
	inline int  GetnPanel(void) {return nPanel;}
	inline int  GetnCohort(void) {return nCohort;}
	inline void SetvGin(string s) {vGin.push_back(s);}
	inline void SetvPin(string s) {vPin.push_back(s);}
	inline void SetvPos(string s) {vPos.push_back(s);}
	inline void SetfnOutput(string s) {fnOutput.assign(s);}
	inline string GetfnOutput(void) {return fnOutput;}
	inline void SetfnRsPos(string s) {fnRsPos.assign(s);}

	inline void SetnK(int s) {nK = s;}//only in CommandLine(...)
	inline void SetnS(int s) {nS = s;}//only in CommandLine(...)
	inline void SetrandSeed(long int s){randSeed = s;}
	inline void SetnEMRuns(int s){nEMRuns = s;}
	inline void SetnWarmSteps(int s){nWarmSteps = s;}
	inline void SetnMaxSteps(int s){nMaxSteps = s;}
	inline void SetnDip(int s){nDip = s;}
	inline void SetnHap(int s){nHap = s;}
	inline void SetnLoci(int s){nLoci = s;}

	inline int GetnMaxSteps(void) {return nMaxSteps;}
	inline int GetnK(void) {return nK;}
	inline int GetnS(void) {return nS;}
	inline int Getnum_Gin(void) {return (int) vGin.size();}
	inline int Getnum_Pin(void) {return (int) vPin.size();}
	inline long int GetrandSeed(void){return randSeed;}
	inline int GetnLoci(void){return nLoci;}   
										//used by fileio.cpp to allocate memory for snpGT. 
	inline int GetnDip(void){return nDip;}
	inline int GetnHap(void){return nHap;}
	inline int GetnEMRuns(void){return nEMRuns;}
	inline int GetnWarmSteps(void) {return nWarmSteps;}

	//interface
    void InitModelParam(void); 			//initialize parameters.

	void set_minit(int s) {m_init = s;}
	void set_npath(int s) {m_npath = s;}
	int get_npath() {return m_npath;}
	void play_ground(void); 
	void simulate_admixed(void); 
	void set_impute(int s) {m_impute = s;}
	void set_upper_switches(double s) {m_upper_switches=s;}
	void set_mixgen(double s) {m_mixgen=s;}
	void set_rrt(double s) {m_rrt=s; }
	void set_morgan(double s) {m_morgan=s; }
	void assign_pop_label(void); 
	       
	void EM_joint(int);     
	void compute_snp_dstr(void); 
	void compute_ps(int, int); 
	void compute_pk(int, int); 
	void write_admix(void); 
	void read_admix(string); 

	void write_em(int); 
	int  read_em(string); 
	int  read_rs2pos(string, long, long, long);  
	int  read_rs2map(string);  
	int  read_bimbam_genotype(long, long);
	int  read_bimbam_phenotype(void);

    int merger(std::string, int **, int, int&, vector<char>&,  int[]);   
		
	void write_genotype(int, int);
	void open_log(void); 
	void close_log(void); 
	void get_het(void); 
	double calc_bf_hpass(int, int, double*, double**, double*, int, double**);
	double calc_bf_hpassQR(int, int, double*, double**, double*, int, double**);
	double logistic_hpass(int, int ns, double * vph, double ** mgt, double * beta, int nc, double ** cov);  
	double logistic_hpassQR(int, int ns, double * vph, double ** mgt, double * beta, int nc, double ** cov, double&);  
	void leading_eigenvec(int, int, double *, double **, double *); 
	double mixmodel_hpass(int ni, int ns, double * vph, double ** mgt, double * beta, int nc, double ** cov);  
	void hpass_pk(int); 
	void hpass_pkbf(int); 
	void hpass_mgt(int); 
	void compute_lhs(int); 
	void combine_lhs(int); 
	void detect_error(int); 
	void selection(int); 
};    

#endif

