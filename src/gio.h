#ifndef __GIO_H__                
#define __GIO_H__

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


class mmSNP 
{
public:
	std::string rs; 
	int ac[7];    //?ATCG+-
public:
	mmSNP(){;};
	~mmSNP(){;}; 
	void assign(class mmSNP s); 
	double get_missing_rate(void); 
	double get_maf(void); 
	char major_allele(void);  
	char minor_allele(void); 
	int missing_count(void) { return(ac[0]); }; 
	int major_count(void); 
	int minor_count(void);  
	void add(class mmSNP s); 
	int allele_count(void);  
	void flip_strand(void); 
	void init(std::string trs, int * tac); 
	void print(void); 
};   //for merge SNPs; 

//void write_admix(void); 
//void write_em(int); 
//int  read_em(std::string); 
//int  read_rs2pos(std::string, long, long, long);  
//int  read_rs2map(std::string);  
//int  read_bimbam_genotype(long, long);
//int  read_bimbam_phenotype(void);
//void write_genotype(int, int);
//int merger(std::string, int **, int, int&, std::vector<char>&,  int);   

#endif
