#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <stdio.h>    
#include <stdlib.h>
#include "control.h"
#include "model.h"
#include "indiv.h"
#include "diploid.h"
#include "haploid.h"
#include "fpmath.h"
#include "fp.h" 

using namespace std;
#define com_empty 	0
#define com_geno	1
#define com_pheno	2
#define com_out		3
#define com_pos		4
#define com_rand 	8
#define com_em 		11
#define com_step 	12
#define com_warm 	13
#define com_lower_cluster 14
#define com_upper_cluster 15
#define com_panel_only 16
#define com_minit 19

#define com_upper_switches 22
#define com_read_map 23
#define com_morgan 24
#define com_mu 26
#define com_rr_training 29
#define com_num_path 30 
#define com_mixgen 31 

#define com_rem     35
#define com_sem     36

#define com_version 99
#define com_help    100
#define com_wmg		101
#define com_wbg     102   //write best guess; 
#define com_wgd     103   //write genotype distribution; 
#define com_weg     104

#define com_ps1 110
#define com_ps2 111
#define com_pk1 113
#define com_pk2 114
#define com_lhs 115
#define com_bf 116
#define com_error 117
#define com_selection 118


#define com_file    3002
#define com_doc     3005
#define com_num     3003
#define com_note    3300
#define com_allele_coding 4000
#define com_exclude_maf 4001
#define com_exclude_miss_in_one_file 4002
#define com_exclude_miss 4005
#define com_exclude_nopos 4006

#define com_silence 4007
#define com_verbose 4008
#define com_disable_hap 4009


CtrlParam::CtrlParam(void)
{
	pMD = NULL;
	hcom[""] = com_empty;
	hcom["-g"] = com_geno; 
	hcom["-gen"] = com_geno;
	hcom["-p"] = com_pheno; 
	hcom["-phe"] = com_pheno; 
	hcom["-o"] = com_out;
	hcom["-out"] = com_out;
	hcom["-pos"] = com_pos; // rspos file name;
	hcom["-P"] = com_pos; // rspos file name;
	hcom["-r"] = com_rand;
	hcom["-R"] = com_rand;
	hcom["-rand"] = com_rand;
//1* for em related. 	
	hcom["-e"] = com_em;
	hcom["-em"] = com_em;
	hcom["-s"] = com_step;
	hcom["-ss"] = com_step;
    hcom["-w"] = com_warm;
	hcom["-c"] = com_lower_cluster;
	hcom["-C"] = com_upper_cluster; 
	hcom["-k"] = com_lower_cluster; 
	hcom["-K"] = com_upper_cluster; 
	hcom["-rr"] = com_upper_switches; 
	hcom["-rrt"] = com_rr_training; 
	hcom["-rm"] = com_read_map; 
	hcom["-morgan"] = com_morgan; 
	hcom["-mixgen"] = com_mixgen; 
	hcom["-mg"] = com_mixgen; 
	hcom["-npath"] = com_num_path; 
	hcom["-np"] = com_num_path; 

	hcom["--ps1"] = com_ps1; 
	hcom["--ps2"] = com_ps2; 
	hcom["--pk1"] = com_pk1; 
	hcom["--pk2"] = com_pk2; 
	hcom["--lhs"] = com_lhs; 
	hcom["--bf"] = com_bf; 
	hcom["--error"] = com_error; 
	hcom["--sel"] = com_selection; 

	hcom["-mu"] = com_mu; 
	hcom["-rem"] = com_rem; 
	hcom["--rem"] = com_rem; 
	hcom["-sem"] = com_sem; 
	hcom["--sem"] = com_sem; 
	hcom["-minit"] = com_minit; 
	hcom["--minit"] = com_minit; 
//help;     
	hcom["-v"] = com_version; 
	hcom["-ver"] = com_version; 
	hcom["-h"] = com_help;
	hcom["-help"] = com_help; 
	hcom["--help"] = com_help; 
	hcom["-wmg"] = com_wmg; 	//write mean genotype
	hcom["--wmg"] = com_wmg; 	//wriet mean genotype
	hcom["-wbg"] = com_wbg;
	hcom["--wbg"] = com_wbg; 
	hcom["-wgd"] = com_wgd; 
	hcom["--wgd"] = com_wgd; 
	hcom["--weg"] = com_weg; 
	hcom["-weg"] = com_weg; 
	hcom["-file"] = com_file; 
	hcom["-FILE"] = com_file; 
	hcom["-doc"] = com_doc; 
	hcom["-DOC"] = com_doc; 
	
	hcom["-num"] = com_num;
	hcom["-note"] = com_note;
	hcom["--note"] = com_note;
	hcom["-ac"] = com_allele_coding; 
	hcom["-allele-coding"] = com_allele_coding; 
	hcom["-exclude-maf"] = com_exclude_maf; 
	hcom["--exclude-maf"] = com_exclude_maf; 
	hcom["-exclude-miss"] = com_exclude_miss; 
	hcom["--exclude-miss"] = com_exclude_miss; 
	hcom["-exclude-miss1"] = com_exclude_miss_in_one_file; 
	hcom["--exclude-miss1"] = com_exclude_miss_in_one_file; 
	hcom["-exclude-nopos"] = com_exclude_nopos; 
	hcom["--exclude-nopos"] = com_exclude_nopos; 
	hcom["--silence"] = com_silence; 
	hcom["-silence"] = com_silence; 
	hcom["--verbose"] = com_verbose; 
	hcom["-verbose"] = com_verbose; 
	hcom["--disable-hap"] = com_disable_hap; 
} 

CtrlParam::~CtrlParam(void)
{
	;
} 

void CtrlParam::OnHelp(void)
{
	PrintHeader(); 
	cout << " " << endl;
	cout << " FILE I/O RELATED OPTIONS." << endl;
	cout << " -g <filename>    " << " repeatable for multiple input files. must pair with -p" << endl;
	cout << " -p <0/1/10/11/12/...>  " << " must pair with -g, <1> pairing genotypes are cohort" << endl; 
	cout << "    <0> pairing genotypes are panel but no labelling, <10/11/12/...> pairing genotype are labelled panels." << endl;
	cout << " -pos <file>      " << " repeatable for multiple input files" << endl; 
    cout << " -o <prefix>      " << " prefix of all output files, random seeds will be used by default" << endl;  
	cout << " " << endl; 
	cout << " EM RELATED OPTIONS." << endl;
	cout << " -w <num>   " << " specify steps of EM run, default 0, using fast linear approximation." << endl; 
	cout << " -s <num>   " << " specify steps of EM run, default 0, using quadratic algorithm." << endl; 
	cout << " -C <num>         " << " specify number of upper clusters, default 2" << endl; 
	cout << " -c <num>         " << " specify number of lower clusters, default 10" << endl; 
	cout << " -mg <num>         " << " specify number of generations the admixture occured" << endl; 
	cout << " -R <num>         " << " specify random seed, system time by default" << endl; 
	cout << " -sem <num>      " << " to save EM results to prefix.em.txt, if 0 save after warm up EM, if 1 save after joint EM." << endl; 
	cout << " -rem <file>      " << " to read EM results" << endl; 
	cout << " --minit <1/2>      " << " 2: initalize at all markers, 1: (default) only initialize a random marker" << endl; 
	cout << " " << endl; 
//	cout << " --ps1       " << " output marginal p(s|data) for diploid, mean local ancestry. default enabled."  << endl; 
//	cout << " --ps2       " << " output joint p(s1,s2|data) for diploid."  << endl; 
//	cout << " --pk1       " << " output marginal p(k|data) for diploid individual."  << endl; 
//	cout << " --pk2       " << " linear algorithm for model fitting, viterbi path to phase and compute p(k|haplotype)"  << endl; 
//	cout << " " << endl;
//	cout << " " << endl;
	cout << " OTHER OPTIONS." << endl;
	cout << " -v(ver)          " << "print version and citation" << endl;
	cout << " -h(help)         " << "print this help" << endl;
//	cout << " -gene <file>     " << " to read gene file that specify regions of interests" << endl; 
//	cout << " -gf(GF) <num>    " << " pair with -gene to specify gene flanking region in kb" << endl; 
	cout << " -exclude-maf <num> " << endl; 
	cout << "          exclude SNPs whose maf < num, default 0.01" << endl; 
	cout << " -exclude-miss <num> " << endl; 
	cout << "          exclude SNPs whose missing rate > num, default 1" << endl; 
	cout << " -exclude-miss1 " << endl; 
	cout << "          exclude SNPs that are missing in one file" << endl; 
	cout << " -exclude-nopos " << endl; 
	cout << "          exclude SNPs that has no position information" << endl; 
	cout << " --silence        " << "no terminal output" <<  endl; 
	cout << " --disable-hap    " << "ignore phase information (=) in genotype input file" <<  endl; 
	cout << " " << endl; 
	cout << " " << endl; 
}

void CtrlParam::PrintHeader(void)
{
	cout << endl; 
	cout << " ELAI version 1.01, visit http://www.haplotype.org for possible update." << endl;
	cout << " Developed by Yongtao Guan ytguan@gmail.com. All Rights Reserved." << endl; 
	cout << " Last update: 13 May 2021." << endl; 
	cout << endl; 
}

void CtrlParam::BatchRun(int argc, char ** argv)
{
	pMD = new class ModelnData;
                             
	int npath = 2; 
	int silence = 0; 
	int yes_weg = 0; 
	int yes_wmg = 0; 
	int yes_wbg = 0;   
	int yes_wgd = 0; 
	int all_weg = 0; 
	int all_wmg = 1; 
	int all_wbg = 0;   //best guess default only write cohort snps; 
	int all_wgd = 1; 
	int file = 0; 
	int rem = 0; 
	int sem = 1; 
	int ps1 = 1; 
	int ps2 = 0; 
	int pk1 = 0; 
	int pk2 = 0; 
	int lhs = 0; 
	int bf = 0; 
	int error = 0; 
	int sel = 0; 
	string fnEM; 
	string note; 
    string fnmap; 

	for(int i = 1; i < argc; i++) 
	{   
		string str;  
		string opt; 
		if(argv[i][0] != '-') 
			continue;
		
		str.assign(argv[i]);
		opt.assign(str, 0, str.length());

		map<string, int> :: iterator iter; 
		iter = hcom.find(opt); 
		if(iter == hcom.end())
		{
			cout << "-elai: unknown option: " << opt << endl; 
			safe_exit(); 
		}
		
		switch (hcom[opt]) 
		{			
			case com_silence:
				silence = 1; 
				pMD->Setsilence(1);
				break; 
			case com_minit:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				pMD->set_minit(atoi(argv[i+1])); 
				break; 
			case com_verbose:
				pMD->Setverbose();
				break; 
			case com_exclude_maf:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') 
					continue;
				{
					double exmaf = atof(argv[i+1]); 
					if(exmaf == 0)
						pMD->Setexcludemaf(-100);
					else 
						pMD->Setexcludemaf(exmaf);
				}
				break; 
			case com_exclude_miss:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				pMD->Setexcludemiss(atof(argv[i+1]));
				break; 

			case com_exclude_miss_in_one_file:
				pMD->Setexcludemiss1();
				break; 

			case com_exclude_nopos:
				pMD->Setexcludenopos();
				break; 
			case com_allele_coding:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				pMD->Setallelecoding(atoi(argv[i+1]));
				break; 
			case com_note:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				note.assign(argv[i+1]);
				break;
			case com_file:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
#ifndef WINDOWS
				str.clear();
#else
				str.resize(0); 
#endif
				str.assign(argv[i+1]);
				pMD->SetfnFILE(str);
				file = 1; 
				break;
			case com_doc:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				str.clear();
				str.assign(argv[i+1]);
				pMD->SetfnDOC(str);
				break;
			case com_geno:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
#ifndef WINDOWS
				str.clear();
#else
				str.resize(0); 
#endif
				str.clear();
				str.assign(argv[i+1]);
				pMD->SetvGin(str);
				break;
			case com_pheno:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				str.assign(argv[i+1]);
				pMD->SetvPin(str);
				break;
			case com_out:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				str.assign(argv[i+1]);
				pMD->SetfnOutput(str);
				break;
			case com_pos:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				str.clear(); 
				str.assign(argv[i+1]);
				pMD->SetvPos(str); 
				break;
			case com_rand: // seed 
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->SetrandSeed(atoi(argv[i+1]));
				break;
			case com_em:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->SetnEMRuns(atoi(argv[i+1]));
				break;
            case com_num:
                if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
                if(!isdigit(argv[i+1][0]))
                {
                    cout << "wrong augument after option." << endl;
                    exit(0);
                }
                pMD->Setnum(atof(argv[i+1]));
                break;

			case com_step:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->SetnMaxSteps(atoi(argv[i+1]));
				break;			
			case com_num_path:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				npath = atoi(argv[i+1]);
				pMD->set_npath(npath);
				break;			
			case com_warm:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->SetnWarmSteps(atoi(argv[i+1]));
				break;
			case com_rem:
				rem = 1;
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				fnEM.assign(argv[i+1]);
				break;
			case com_sem:
				sem = 1;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				sem += atoi(argv[i+1]);
				break;
			case com_lower_cluster:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->SetnK(atoi(argv[i+1]));
				break;
			case com_upper_cluster:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->SetnS(atoi(argv[i+1]));
				break;
			case com_upper_switches:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->set_upper_switches(atof(argv[i+1]));
				break;
			case com_rr_training:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->set_rrt(atoi(argv[i+1]));
				break;
			case com_mixgen:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->set_mixgen(atof(argv[i+1]));
				break;
			case com_morgan:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->set_morgan(atof(argv[i+1]));
				break;
			case com_mu:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				if(!isdigit(argv[i+1][0]))
				{
					cout << "wrong argument after option." << endl; 
					exit(0); 
				}
				pMD->set_mu(atof(argv[i+1]));
				break;
			case com_read_map:
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				{
					fnmap.clear(); 
					fnmap.assign(argv[i+1]);
				}
				break;
				
			case com_help:
				OnHelp();
				safe_exit();
				break;

			case com_version:
				PrintHeader(); 
				safe_exit(); 
				break; 
				
			case com_weg:
				yes_weg = 1; 
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				all_weg = atoi(argv[i+1]); 
				break;
			case com_wmg:
				yes_wmg = 1; 
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				all_wmg = atoi(argv[i+1]); 
				break;
			case com_wbg:
				yes_wbg = 1; 
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				all_wbg = atoi(argv[i+1]); 
				break;
			case com_wgd:
				yes_wgd = 1; 
				if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
				all_wgd = atoi(argv[i+1]); 
				break;

			case com_ps1:
				ps1 = 1; 
				break;
			case com_ps2:
				ps2 = 1; 
				break;
			case com_pk1:
				pk1 = 1; 
				break;
			case com_pk2:
				pk2 = 1; 
				break;

			case com_lhs:
				lhs = 1; 
				break;
			case com_bf:
				bf = 1; 
				break;

			case com_error:
				error = 1; 
				break;

			case com_selection:
				sel = 1; 
				break;
			case com_disable_hap:
				pMD->set_disable_hap(); 
				break;
   			
			default:
				fprintf(stderr,"Bad option %s\n", argv[i]);
				OnHelp();
				exit(0);
				break; 
		}              
	}
    

	streambuf * cout_strbuf(cout.rdbuf()); 
	//save current cout into a streambuf; 
	ostringstream output; 
	if(silence == 1) 
		cout.rdbuf(output.rdbuf()); 
	//redirect cout to a streambuf output; 


	int randseed = pMD->GetrandSeed();  // default returns -1;
	if (randseed <= 0)
		randseed = (unsigned) time(NULL);
	gsl_rng_set(gsl_r, randseed);


	string prefix(pMD->GetfnOutput()); 
	if(prefix.compare("\0") == 0) {
		char tmp[100];
		memset(tmp, '\0', 100);
		sprintf(tmp, "%d", randseed); 
		prefix.assign(tmp); 
		pMD->SetfnOutput(prefix); 
	}
	//if no random seed provided, use system clock;    
	//if no prefix provided, use randseed; 

	pMD->open_log(); 
	(pMD->fplog) << "## " << note << endl;
	(pMD->fplog) << "## COMMAND: "; 
	for (int i = 0; i < argc; i++)
		(pMD->fplog) << argv[i] << " "; 
	(pMD->fplog) << endl; 
	(pMD->fplog) << "## randseed = " << randseed << endl; 
	cout << "## randseed = " << randseed << endl; 
	///////////////////////////////
	
	int	read_dip = pMD->read_bimbam_genotype(0, -1); 
	if(read_dip == 0)
	{
		cout << "-elai: no valid input." << endl; 
		exit(0); 
	}
	///////////////////////////////

	time_t sec_beg, sec_end; 
	sec_beg = time(NULL); 
	
	{
		if(fnmap.size() > 0)
			pMD->read_rs2map(fnmap); 

		if(yes_weg == 1) 
			pMD->write_genotype(0, all_weg); 
                
//                cout << "rem = " << rem << endl; 

		if(rem == 0)
			pMD->InitModelParam(); 
		else 
			pMD->read_em(fnEM); 

		pMD->EM_joint(1);  
		pMD->EM_joint(2);  
		pMD->write_em(1); 
		pMD->write_admix(); 
				  
//		if(ps1 == 1) 
			pMD->compute_ps(1, 2);  //output mean, always use diploid to compute;  
		if(ps2 == 1) 
			pMD->compute_ps(2, 2);  //output prob. always use diploid to compute; 
		//out_mode; comp_mode; 
	
		if(pk1 == 1) 
			pMD->compute_pk(1, 2);                
		if(lhs == 1) 
			pMD->compute_lhs(2); 
	}

	sec_end = time(NULL); 
	(pMD->fplog) << "## EM seconds used = " << sec_end - sec_beg << endl; 
	
	
	

	{
		pMD->fplog << "## ELAI generate following files in the output directory." << endl; 
		pMD->fplog << "## " << prefix << ".snpdata.txt" << endl; 
		if(yes_weg == 1) pMD->fplog << "## "<< prefix << ".exact.genotype.txt" << endl;
		if(lhs == 1) pMD->fplog << "## "<< prefix << ".lhs.txt" << endl; 
		
		(pMD->fplog) << "## random seed = " << randseed << endl;
		pMD->close_log();
		cout << "-elai: finished, for details see log: "<<  prefix << ".log.txt" << endl; 
	} 


	if(silence == 1)
		cout.rdbuf(cout_strbuf); 
	//restore the original cout; 
//	delete pMD; 
}

