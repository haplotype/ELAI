#include "gio.h"
#include "model.h"
#include "anchap.h"
#include <algorithm> 

void mmSNP::init(std::string trs, int * tac) {
	rs.assign(trs); 
	for (int j = 0; j < 7; j++)
		ac[j] = tac[j]; 
}

void mmSNP::assign(class mmSNP s) {
	rs.assign(s.rs);
	for (int j = 0; j < 7; j++) 
		ac[j] = s.ac[j];
}

char mmSNP::major_allele(void) {
	int tmax = -1; 
	int index = 0; 
	for (int j = 1; j < 7; j++) 
	{
		if(ac[j] > tmax) {
			tmax = ac[j]; 
			index = j; 
		}
	} 
	char code[7] = {'?', 'A', 'T', 'C', 'G', '+', '-'}; 
	return(code[index]); 
}

//note: in major allele the comparison is strictly larger; 
// in minor allele the comparison is <=. 
// this happens to workout with maf = 0.5. 

char mmSNP::minor_allele(void) {
	char ret = '?'; 
	if(allele_count() == 2) {
		double tmin = 1e10; 
		int index = 0; 
		for (int j = 1; j < 7; j++) 
		{
			if(ac[j] <= tmin && ac[j] > 0) {
				tmin = (double) ac[j]; 
				index = j; 
			}
		} 

		char code[7] = {'?', 'A', 'T', 'C', 'G', '+', '-'}; 
		ret = code[index]; 
	}
	return(ret); 
}

int mmSNP::major_count(void) {
	int tmax = 0; 
	for (int j = 1; j < 7; j++) 
	{
		if(ac[j] > tmax) {
			tmax = ac[j]; 
		}
	} 
	return(tmax); 
}

int mmSNP::minor_count(void) {
	int ret = 0; 
	if(allele_count() == 2) {
		double tmin = 1e10; 
		for (int j = 1; j < 7; j++) 
		{
			if(ac[j] <= tmin && ac[j] > 0) {
				tmin = (double) ac[j]; 
			}
		} 
		ret = (int) tmin; 
	}
	return(ret); 
}

void mmSNP::add(class mmSNP s) {
	if(rs.compare(s.rs) == 0) {
		for (int j = 0; j < 7; j++)
			ac[j] += s.ac[j]; 
	}
	else {
		cout << "-warning: try to combine alleles of different SNPs" << endl; 
		cout << rs << " " << s.rs << endl; 
		exit(0); 
	}
}

int mmSNP::allele_count(void) {
	int count = 0; 
	for (int j = 1; j < 7; j++)
	{
		if(ac[j] > 0) count ++; 
	}
	return count; 
}

void mmSNP::flip_strand(void) {
	for (int m = 0; m < 3; m++)
	{
		int a = ac[2*m+1]; 
		ac[2*m+1] = ac[2*m+2]; 
		ac[2*m+2] = a; 
	}
}

double mmSNP::get_missing_rate(void) {
	double sum=0; 
	for (int j = 0; j < 7; j++) 
		sum += (double)ac[j]; 
	return(ac[0]/sum); 
}

double mmSNP::get_maf(void) {
	double ret = 0; 
	if(allele_count() == 2) {
		char a = minor_allele(); 
		char code[7] = {'?', 'A', 'T', 'C', 'G', '+', '-'}; 
		int wh = 0; 
		for (int j = 1; j < 7; j++)
		{
			if(code[j] == a)  {
				wh = j; 
				break; 
			}
		}
		
		double sum = 0; 
		for (int j = 1; j < 7; j++)
			sum += ac[j]; 
		ret = ac[wh] / sum; 
	}
	return ret; 
}

void mmSNP::print(void) {
	cout << rs << ": "; 
	for (int j = 0; j < 7; j++) 
		cout << ac[j] << " "; 
	cout << endl; 
}

bool rscomp(pair<string, pair<long, int> > a, pair<string, pair<long, int> > b)
{
	if(a.second.second == b.second.second)
		return (a.second.first < b.second.first); 
	else 
		return(a.second.second < b.second.second); 
}

int ModelnData::merger(std::string rs, int ** mac, int nf, int& missing_in_one_file, vector<char>& allele_coding, int fac[7])  
{
	int ret = 1; 
	missing_in_one_file = 0; 
	int * flip = new int[nf]; 
	for (int i = 0; i < nf; i++)
	{
		allele_coding[i] = 'X'; 
		flip[i] = 0; 
	}

//	if(rs.compare("rs656127") == 0) 
//	{
//		for (int i = 0; i < nf; i++) 
//		{
//			for (int j = 0; j < 7; j++)
//				cout << mac[i][j] << " "; 
//			cout << endl; 
//		}
//	}
	int * vac = new int[nf]; //vector of allele count; 
	for (int i = 0; i < nf; i++)
	{
		class mmSNP ts; 
		ts.init(rs, mac[i]); 
		vac[i] = ts.allele_count(); 
		if(vac[i] == 0) missing_in_one_file ++; 
	}

	for (int j = 0; j < 7; j++)
	{
		int sum = 0; 
		for (int i = 0; i < nf; i++)
			sum += mac[i][j]; 
		fac[j] = sum; 
	}
	class mmSNP fmm; 
	fmm.init(rs, fac); 
	int ac = fmm.allele_count();  
	if(ac == 2 || ac == 1)
	{
		for (int i = 0; i < nf; i++) 
			allele_coding[i] = fmm.major_allele(); 
	}
	else if (ac == 0) 
	{
		ret = 0; 
	}
	else { //ac == 3 or 4; 
		class mmSNP * pmm = new class mmSNP[nf]; 
		vector<int> index1; 
		vector<int> index2; 
		for (int i = 0; i < nf; i++)
		{
			pmm[i].init(rs, mac[i]); 
			vac[i] = pmm[i].allele_count(); 
			if(vac[i] == 1) index1.push_back(i); 
			else if(vac[i] == 2) index2.push_back(i); 
			else if(vac[i] > 2) ret = 0;  
			// a file has three allele, snp merge unsuccessful!
		}

		class mmSNP tmm; 
		int dummy[7] = {0}; 
		tmm.init(rs, dummy); 
		for (unsigned i = 0; i < index2.size(); i++)
		{
			int ii = index2.at(i); 
			class mmSNP t1; 
			t1.assign(pmm[ii]); 
			t1.add(tmm); 
			int ac1 = t1.allele_count(); 

			class mmSNP t2; 
			t2.assign(pmm[ii]); 
			t2.flip_strand(); 
			t2.add(tmm); 
			int ac2 = t2.allele_count(); 

            if(i == 0) 
				tmm.assign(t1); 
			else {
				if((ac1 == 2 && ac2 != 2)) 
					tmm.assign(t1); 
				else if(ac1 != 2 && ac2 == 2) 
				{
					flip[ii] = 1; 
					tmm.assign(t2); 
				}
				else if(ac1 == 2 && ac2 == 2) 
				{
					cout << rs << " is an AT/CG SNP, keep it as is" << endl; 
					fplog << rs << " is an AT/CG SNP, keep it as is" << endl; 
				}
				else { 
					ret = 0; 
					break; 
				}
			}
		}
		//grant;
//			if(rs.compare("rs656127") == 0) 
//			{
//				for (int i = 0; i < nf; i++)
//					cout << flip[i] << " "; 
//				cout << endl; 
//			}

		for (unsigned i = 0; i < index1.size(); i++)
		{
			int ii = index1.at(i); 
			class mmSNP t1; 
			t1.assign(pmm[ii]); 
			t1.add(tmm); 
			int ac1 = t1.allele_count(); 

			class mmSNP t2; 
			t2.assign(pmm[ii]); 
			t2.flip_strand(); 
			t2.add(tmm); 
			int ac2 = t2.allele_count(); 

			if(ac1 == 2 && ac2 != 2) 
				tmm.assign(t1); 
			else if(ac2 == 2 && ac1 != 2) 
			{
				flip[ii] = 1; 
				tmm.assign(t2); 
			}
			else if (ac1 == 2 && ac2 == 2) 
			{
				cout << rs << " is flip compatiable, keep it as is; but please check" << endl; 
				fplog << rs << " is flip compatiable, keep it as is; but please check" << endl; 
			}
		}

		char ref = tmm.major_allele(); 
	    map<char, char> cc;
		cc['A'] = 'T'; cc['T'] = 'A';
		cc['C'] = 'G'; cc['G'] = 'C';
		cc['N'] = 'N'; cc['?'] = '?';
		cc['+'] = '-'; cc['-'] = '+';

		//this is combined reference major allele; 
		if(ret == 1) {
			for (int i = 0; i < nf; i++)
			{
				if(flip[i] == 0) 
					allele_coding[i] = ref; 
				else 
					allele_coding[i] = cc[ref]; 
			}
			for (int j = 0; j < 7; j++)
				fac[j] = tmm.ac[j]; 
		}
		delete[] pmm; 
	}
	delete[] vac; 
	delete[] flip; 
	return ret; 
}

int ModelnData::read_bimbam_genotype(long start_pos, long end_pos)
{   
//	cout << "m_exclude_miss1 = " << m_exclude_miss1 << endl; 
//	cout << "exclude_miss = " << m_exclude_miss << endl; 
	if(vGin.size() == 0 || vGin.size() != vPin.size()) 
	{
		fplog << "## warning: genotype phenotype files not match" << endl; 
		cout << "-warning: genotype phenotype files not match" << endl; 
		return 0; 
	}

	int which = -1; 
	for (unsigned i = 0; i < vPin.size(); i++)
	{
		if(vPin.at(i).compare("0") == 0) //equal; 
		{
            which = i; 
			break; 
		}
	}
    if(which > 0) 
	{
		string ts; 
		ts.assign(vPin.at(0)); 
		vPin.at(0).assign(vPin.at(which)); 
		vPin.at(which).assign(ts); 
		ts.assign(vGin.at(0)); 
		vGin.at(0).assign(vGin.at(which)); 
		vGin.at(which).assign(ts); 
	}
	//exchange file names so that the panel is at the beginning; 

	fplog << "## warning: number of position files = " << vPos.size() << endl; 
	cout << "-warning: number of position files = " << vPos.size() << endl; 
	unsigned readpos = 0; 
	for(unsigned f = 0; f < vPos.size(); f++)
	{
		fnRsPos.assign(vPos.at(f)); 
		readpos += read_rs2pos(fnRsPos, start_pos, end_pos, nGeneFlank); 
	}
	if(readpos == 0 && vGin.size() > 1) 
	{
		cout << readpos << endl; 
		fplog << "## warning: requires position file(s), all markers without position info will be excluded." << endl; 
		cout << "-warning: requires position file(s), all markers without position info will be excluded." << endl; 
		safe_exit(); 
	}   //initialize mapRs2pos; 
	if(readpos < vPos.size()) 
	{
		fplog << "## warning: one of position files failed" << endl; 
		cout << "-warning: one of position files failed" << endl; 
	}
	fplog << "## warning: position files contain " << mapRs2pos.size() << " records." << endl; 
	cout << "-warning: position files contain " << mapRs2pos.size() << " records." << endl; 

	//first pass, go over all the files collecting summary for each SNP. 
	//allele type, allele count, etc. treating phased panel as unphased (but get the correct nHap, nDip);  
	//merge SNPs based on summary. generate a hash keyed by rs whose content has reference allele
	//for different files (so that genotypes are the counts of reference alleles). If the reference
	//allele is '?', the whole SNP for that file were to set 'QQ'; 
	//sort the key according to the physical position, to get vsRsnum; 
	//Allocate memory for each individual, i.e., we have nIndiv x nLoci; 
	 
	nHap = 0; 
	nDip = 0; 
	nLoci = 0; 
	map<string, int> mrs2yes; 
	map<int, map<string, class mmSNP> > mfile2map; 
	vector<int> vindploid; 

	for (unsigned nf = 0; nf < vGin.size(); nf++)
	{   
		int filerows = 0; 
		int filecols = 0; 
		int filecols_diff = 0; 
		int isPanel = 0; 
		if (vPin.at(nf).compare("0") == 0)  
			isPanel = 1; 
		map<string, class mmSNP> rs2ss; //rs to snp summary; 
		
		ifstream infile; 
		streambuf * pbuf;
		char delimit[] = ";, \t";
		infile.open(vGin.at(nf).c_str(), ios::in);
		if(!infile.is_open()) 
		{
			fplog << "## warning: cannot open genotype file: " << vGin.at(nf) << endl; 
			cout << "-warning: cannot open genotype file: " << vGin.at(nf) << endl; 
			safe_exit(); 
		} 
		pbuf = infile.rdbuf();
		// open file; 

		string line; 
		int ni = 0; 
		{
			line.assign(getline(pbuf)); 
			char *res = strtok((char*)line.c_str(), delimit); 
			if(res != NULL) ni = atoi(res);
			
			string::size_type loc; 
			loc = line.find(".", 0);
			if (loc != string::npos || ni < 1)
			{
				fplog << "## warning: confused -p with -g ? " << endl; 
				cout << "-warning: confused -p with -g ? " << endl; 
				safe_exit(); 
			}
			
			loc = line.find("=", 0);
			if( loc != string::npos && m_disable_hap == 0) 
				vFilePloid.push_back(1); 
			else 
				vFilePloid.push_back(2); 
		}   //number of individuals;
		
		int ns = 0; 
		{
			line.assign(getline(pbuf)); 
			char *res = strtok((char*)line.c_str(), delimit); 
			if(res != NULL) ns = atoi(res);
		}	//number of snps;
				
		string tline(getline(pbuf)); 
		line.assign(tline); 
		char * res = strtok((char*)tline.c_str(), delimit); 
		string head(res);
		if(head.compare("IND") == 0 || head.compare("rsID") == 0)
		{
			line.assign(getline(pbuf)); 
		}   //individual ids
			
		while(!line.empty())
		{   
			//printf("%s\n", line.c_str()); 
			filerows ++; 
			char *res = strtok((char*)line.c_str(), delimit); 
			if(res == NULL) break; 
			string rs;
			rs.assign(res);  

			map<string, long> :: iterator iter; 
			iter = mapRs2pos.find(rs); 
			if(iter == mapRs2pos.end() || iter->second < 0) 
			{
				line.assign(getline(pbuf));
				continue; 
			}   //only process those SNPs that have position info. 

			mrs2yes[rs] = 1;
			//mark the snp; 
			
			filecols = 0; 
			int a2c[7] = {0};   //strictly in the following order  //? A T C G + - 
			for (int j = 0; j < ni; j++)    
			{          
				string tmp;
				res = strtok(NULL, delimit); 
				if(res == NULL) break; 
				tmp.assign(res);
				for (int i = 0; i < 2; i++)
				{
					if (tmp.at(i) == 'N' || tmp.at(i) == '0' || tmp.at(i) == '?') 
						a2c[0]++; 
					else if (tmp.at(i) == 'A') 
						a2c[1]++; 
					else if (tmp.at(i) == 'T') 
						a2c[2]++; 
					else if (tmp.at(i) == 'C') 
						a2c[3]++; 
					else if (tmp.at(i) == 'G') 
						a2c[4]++; 
					else if (tmp.at(i) == '+' || tmp.at(i) == 'I') 
						a2c[5]++; 
					else if (tmp.at(i) == '-' || tmp.at(i) == 'D') 
						a2c[6]++; 
					else {
						cout << "-warning: undefined letters; only ?N0ATCG+- are valid." << endl;  
						fplog << "## warning: undefined letters; only ?N0ATCG+- are valid." << endl;  
						exit(0); 
					}
				}                   
				filecols++; 
			}   // snp summary in a2c;    //? A T C G + -
			if(filecols != ni) filecols_diff++; 

			mmSNP ss; 
			ss.init(rs, a2c); 

			map<string, class mmSNP> :: iterator mmiter; 
			mmiter = rs2ss.find(rs); 
			if (mmiter != rs2ss.end())
			{
				fplog << "## WARNING: Two SNPs have the same rsnum " << mmiter->first << endl; 
				cout << "-warning: two snp have same rsnum " << mmiter->first << endl; 
			}
			
			rs2ss.insert(pair<string, class mmSNP> (rs, ss)); 

			line.assign(getline(pbuf));
		}

		infile.close();
		mfile2map[nf] = rs2ss; 
		
		if(filerows != ns) 
		{
			fplog << "## WARNING: Number of SNPs in file is different from the second line" << endl; 
			cout << "-warning: actual snp number differ from the second line" << endl; 
		}
		if(filecols != ni) 
		{
			fplog << "## WARNING: Number of ind's is different from the first line " << filecols << " " << ni << endl; 
			cout << "-warning: actual indiv number differ from the first line " << filecols << " " << ni << endl; 
		}
		if(filecols_diff > 0)
		{
			fplog << "## WARNING: " << filecols_diff << " SNPs has different # of ind's from the first line" << endl; 
			cout << "-warning: in " << filecols_diff << " snps, individual differ from the first line" << endl; 
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////

		if(vFilePloid.at(nf) == 1)
		{
			vFileIndiv.push_back(2 * ni);
			nHap += 2 * ni; 
			for(int i = 0; i < 2 * ni; i++)
				vindploid.push_back(1); 
		}
		else 
		{
			vFileIndiv.push_back(ni);
			nDip += ni; 
			for(int i = 0; i < ni; i++)
				vindploid.push_back(2); 
		}
		fplog << "## warning: File " << nf << " has " << vFileIndiv.at(nf) << " ind's and " \
			<< rs2ss.size() << " SNPs" << endl; 
		cout << "-warning: file " << nf << " has " << vFileIndiv.at(nf) << " individual and " \
			<< rs2ss.size() << " snps" << endl; 
		rs2ss.clear(); 
	}   // SNP summary data;
	
	////////////////////////////////////////////////////////////////////////////
	int count_nomatch = 0;   //failed to match allele types across multiple files;  
	int count_miss1file = 0; //missing in one file; 
	int count_miss = 0;      //too much missing at random; 
	int count_maf = 0;       //failed due to maf too extreme; 
	map<string, int> ::iterator m_iter; 
	map<string, vector<char> > mrs2coding;   //rs to reference allele coding for each file;  
	int ** mac = Allocate2DIntMatrix((int) vGin.size(), 7);  
	for (m_iter = mrs2yes.begin(); m_iter != mrs2yes.end(); ++m_iter)
	{
		string rs; 
		rs.assign(m_iter->first); 

		vector<char> allele_coding; 
		for (unsigned nf = 0; nf < vGin.size(); nf++)
		{    
			for (int j = 0; j < 7; j++)
				mac[nf][j] = 0; 
			map<string, class mmSNP> :: iterator iter; 
			iter = mfile2map[nf].find(rs); 
			if(iter == mfile2map[nf].end())  
			{
				mac[nf][0] = vFileIndiv.at(nf) * vFilePloid.at(nf); 
			}   //if it can't be found then it is missing; 
			else 
			{
				for (int j = 0; j < 7; j++)
					mac[nf][j] = iter->second.ac[j]; 
			}
			allele_coding.push_back('X'); 
		}
		
		int missing_in_one_file = 0; 
		int fac[7] = {0}; 
		int success = merger(rs, mac, (int) vGin.size(), missing_in_one_file, allele_coding, &fac[0]);  
		
		mrs2coding[rs] = allele_coding;                  
		if(success == 0) 
		{
			mrs2yes[rs] = 0; 
			count_nomatch ++; 
		}
		else {
			class mmSNP tss; 
			tss.init(rs, fac); 
			if(m_exclude_miss1 && missing_in_one_file) //missing test; 
			{
				count_miss1file ++; 
				mrs2yes[rs] = -1; 
			}
			else if(tss.get_missing_rate() > m_exclude_miss)   
			{
				count_miss ++; 
				mrs2yes[rs] = -2; 
			}
			else if(tss.get_maf() < m_exclude_maf) 
			{
				count_maf ++;    
				mrs2yes[rs] = -3; 
			}
			else 
				mapRs2mm[rs] = pair<char, char> (tss.major_allele(), tss.minor_allele()); 
		}
		
		vector<char> ().swap(allele_coding); 
	}
	if(count_nomatch > 0) 
	{
		fplog << "## warning: Exclude " << count_nomatch << " SNPs due to failure to match betweeen files." << endl;
		cout << "-warning: exclude " << count_nomatch << " snps due to failure to match betweeen files." << endl;
	}
	if(count_miss1file > 0) 
	{
		fplog << "## warning: Exclude " << count_miss1file << " SNPs due to miss in at least one file " << endl;
		cout << "-warning: exclude " << count_miss1file << " snps due to miss in at least one file " << endl;
	}
	if(count_miss > 0) 
	{
		fplog << "## warning: Exclude " << count_miss << " SNPs due to miss proportion > " << m_exclude_miss << endl;
		cout << "-warning: exclude " << count_miss << " snps due to miss proportion > " << m_exclude_miss << endl;
	}
	if(count_maf > 0) 
	{
		fplog << "## warning: Exclude " << count_maf << " SNPs due to maf < " << m_exclude_maf << endl;
		cout << "-warning: exclude " << count_maf << " snps due to maf < " << m_exclude_maf << endl;
	}
	
	////////////////////////////////////////////////////////////////////////////////
	vector<pair<string, pair<long, int> > > vp;
	map<string, long> :: iterator pos_iter; 
	map<string, int> :: iterator chr_iter; 
	for (m_iter = mrs2yes.begin(); m_iter != mrs2yes.end(); m_iter++)
	{
		string rs(m_iter->first);
		if(m_iter->second == 1) {
			pair<string, pair<long, int> > tmp;
			tmp.first.assign(m_iter->first);
			pos_iter = mapRs2pos.find(m_iter->first);
			if(pos_iter != mapRs2pos.end())
				tmp.second.first = pos_iter->second;
			else 
				tmp.second.first = 0; 

			chr_iter = mapRs2chr.find(tmp.first); 
			if(chr_iter != mapRs2chr.end())
				tmp.second.second = chr_iter->second; 
			else
				tmp.second.second = 0; 
									   
			vp.push_back(tmp); 
		}
	}
	mrs2yes.clear(); 

	stable_sort(vp.begin(), vp.end(), rscomp); 
	//sort the rs pos vector in order of chr && pos;
	vsRsnum.clear();   //snp id in phyiscal order; 
	map<string, int> mrs2index; //snp index; 
	for(int i = 0; i < (int) vp.size(); i++)
	{
		vsRsnum.push_back(vp.at(i).first);
		mrs2index[vp.at(i).first] = i; 
	}
	vector<pair<string, pair<long, int> > >().swap(vp); 

	nLoci = (int) vsRsnum.size(); 
	nIndiv = nDip + nHap; 
	cout << "number of individuals = " << nIndiv << endl; 
	cout << "number of SNPs = " << nLoci << endl; 

	if(nLoci == 0 || nIndiv == 0) 
	{
		cout << "-warning: no valid snp or individual." <<  endl; 
		exit(0); 
	}
	////////////////////////////////////////////////////////////////////////////////////

	if(m_morgan < 0) 
	{
		m_morgan = 0; 
		for (int m = 1; m < nLoci; m++)
		{
			string rs0 = vsRsnum.at(m-1); 
			string rs1 = vsRsnum.at(m); 
			long pos0 = mapRs2pos[rs0]; 
			long pos1 = mapRs2pos[rs1]; 
			int chr0 = mapRs2chr[rs0]; 
			int chr1 = mapRs2chr[rs1]; 

			if(chr0 == chr1) 
				m_morgan += (double) (pos1 - pos0); 
			else 
				m_morgan += 1e8; 
		}           
		m_morgan /= 1e8; 
	}
    cout << "### m_morgan = " << m_morgan << endl; 
    fplog << "### m_morgan = " << m_morgan << endl; 
	////////////////////////////////////////////////////////////////////////////////////
	
	pIndiv = new Individual * [nIndiv]; 
	for (int i = 0; i < nIndiv; i++)
	{
		if (vindploid.at(i) == 1) {
			pIndiv[i] = new HapInd(nLoci); 
			pIndiv[i]->AllocatesnpGT(nLoci); 
			for (int m = 0; m < nLoci; m++)
				pIndiv[i]->SetsnpGT(m, char('0' + QQ)); 
		} else {
			pIndiv[i] = new DipInd(nLoci);
			pIndiv[i]->AllocatesnpGT(nLoci); 
			for (int m = 0; m < nLoci; m++)
				pIndiv[i]->SetsnpGT(m, char('0' + QQ)); 
		}
	}

	pHap = new AncHap * [nK];  
	for (int i = 0; i < nK; i++)
	{
		pHap[i] = new AncHap; 
	   	pHap[i]->allocate_af(nLoci); 
	}

	
////////////////////////////////////////////////////////////////////////////////////////////////////
	// go over genotype files again; 
	// for each diploid file, scan each SNP, fill in the genotype table. 
	int start_ni = 0; 
	for (int nf = 0; nf < (int)vGin.size(); nf++)
	{       
		streambuf * pbuf;
		ifstream infile; 
		char delimit[] = ";, \t";
		infile.open(vGin.at(nf).c_str(), ios::in);
		if(!infile.is_open()) 
		{
			cout << "-warning: cannot open genotype file: " << vGin.at(nf) << endl; 
			safe_exit(); 
		} 
		pbuf = infile.rdbuf();
		 
		string line; 
		line.assign(getline(pbuf));   //first line; 
		line.assign(getline(pbuf));   //second line; 
		
		string tline(getline(pbuf));   //third line; 
		line.assign(tline); 
		char * res = strtok((char*)tline.c_str(), delimit); 
		string head(res);
		if(head.compare("IND") == 0 || head.compare("rsID") == 0)
			line.assign(getline(pbuf)); 
		//possible individual ids
			
		while(line.size() > 0)
		{   
			char * res = strtok((char*)line.c_str(), delimit); 
			if(res == NULL) break; 
			string rs; 
			rs.assign(res);  
			if(end_pos > 0) 
			{
				map<string, long> :: iterator iter; 
				iter = mapRs2pos.find(rs); 
				if(iter == mapRs2pos.end() || iter->second < 0) 
				{
					line.assign(getline(pbuf));
					continue; 
				}
			}    
			int m = -1; 
			map<string, int> :: iterator s2i; 
			s2i = mrs2index.find(rs); 
			if(s2i == mrs2index.end())
			{
				line.assign(getline(pbuf));
				continue; 
			}
			else 
				m = s2i->second; ; 
			// rs number; 
			char ref = mrs2coding[rs].at(nf); //reference allele; 
			int ni = start_ni; 
			for (int i = 0; i < vFileIndiv.at(nf); i++)
			{                              
				char * res = strtok(NULL, delimit); 
				if(res == NULL) break; 
				string tmp;
				tmp.assign(res);
				if(vFilePloid.at(nf)  == 2) 
				{
					if (ref == 'X' || tmp.at(0) == 'N' || tmp.at(0) == '?' || tmp.at(0) == '0')
					{
						pIndiv[ni]->SetsnpGT(m, char('0'+QQ)); 
						//grant; QQ->2; 
					}
					else 
					{
						int gt = 0; 
						if (tmp.at(0) == ref) gt++;
						if (tmp.at(1) == ref) gt++;
						gt = 2 - gt; 
						//because of coding allele is major, so we flip to minor here; 
						pIndiv[ni]->SetsnpGT(m, char('0'+gt)); 
					}           
					ni++; 
				}
				else if (vFilePloid.at(nf) == 1) 
				{
					if (ref == 'X' || tmp.at(0) == 'N' || tmp.at(0) == '?' || tmp.at(0) == '0')
					{
						 pIndiv[ni]->SetsnpGT(m, char('0'+QQ)); 
						 //grant; QQ->1
					}
					else 
					{
						int gt = 0; 
						if (tmp.at(0) == ref) gt++;
						gt = 1 - gt; 
						pIndiv[ni]->SetsnpGT(m, char('0'+gt)); 
					}                                        
					i++;
					ni++; 
					if (ref == 'X' || tmp.at(1) == 'N' || tmp.at(1) == '?' || tmp.at(1) == '0')
					{
						pIndiv[ni]->SetsnpGT(m, char('0'+QQ)); 
						 //grant; QQ->1
					}
					else 
					{
						int gt = 0; 
						if (tmp.at(1) == ref) gt++;
						gt = 1 - gt; 
						pIndiv[ni]->SetsnpGT(m, char('0'+gt)); 
					}     
					ni++;
				}
			}
			line.assign(getline(pbuf));
		}
		infile.close();
		cout << "-warning: read file " << nf << " again " << endl; 
		start_ni += vFileIndiv.at(nf); 
	}

	////////////////////////////
	read_bimbam_phenotype(); 
	
	int ni = 0; 
	for (int i = 0; i < nIndiv; i++)
	{
		if(pIndiv[i]->GetisPanel()==0 && pIndiv[i]->get_mask()==0)
			ni++; 
	}
//	vector<double> tph; 
//	{
//		std::fstream in;
//		in.open(fnFILE.c_str());
//		if(in.is_open()) {
//			double temp = 0; 
//			while (in >> temp) {
//				tph.push_back(temp);
//			}
//			in.close();
//		}
//		else std::cout << "unable to open file" << std::endl;
//	}
////	int ncase, nctrl; 
////	ncase = nctrl = 0; 
////	for (unsigned i = 0; i < tph.size(); i++)
////	{
////		if(tph.at(i) == 1) ncase++; 
////		else if(tph.at(i) == 0) nctrl++; 
////	}
//
//    cout << ni << " === " << tph.size() << endl; 
	int * tgt = new int[ni]; 

	fstream outfile; 
	string sfn("output/");
	sfn.append(fnOutput);
	sfn.append(".snpinfo.txt");
	outfile.open(sfn.c_str(), ios::out);
	if(outfile.is_open()) 
	{
		outfile << "rs\t pos \t chr \t major\t minor\t maf \t miss \t n0 \t n1 \t n2" << endl;  
		for (int m = 0; m < nLoci; m++)
		{
			string rs = vsRsnum.at(m); 

			char buf[100]; 
			sprintf(buf, "%-s\t", rs.c_str()); 
			outfile << buf; 
			outfile << mapRs2pos[rs] << "\t"; 
			outfile << mapRs2chr[rs] << "\t"; 
			outfile << mapRs2mm[rs].second << "\t" << mapRs2mm[rs].first << "\t";

			int ii = 0; 
			for (int i = 0; i < nIndiv; i++)
			{
				if(pIndiv[i]->GetisPanel() == 0 && pIndiv[i]->get_mask()==0)
					tgt[ii++] = pIndiv[i]->GetsnpGT(m);
			}
			int mi, n0, n1, n2; 
			mi = n0 = n1 = n2 = 0; 
			for (int i = 0; i < ni; i++)
			{
				if(tgt[i] > 2) {
					mi++; 
				}
				else if(tgt[i] == 0) n0++; 
				else if(tgt[i] == 1) n1++; 
				else n2++; 
			}
			double maf = (n2 * 2.0 + n1) / 2.0 / (n2 + n1 + n0); 
			sprintf(buf, "%.3f\t", maf); 
			outfile << buf;
			outfile << mi << "\t" << n0 << "\t" << n1 << "\t" << n2 << endl; 

		}
		outfile.close(); 
	}
	else 
		cout << "-warning: skip writing snpinfo" << endl;

	delete[] tgt; 
/////////////////////////////////////////////////////////////////////////////////////////////////////
	
	fplog << "## number genotype files = " << vGin.size() << endl; 
	fplog << "## number phenotype files = " << vPin.size() << endl; 
	fplog << "## number of diploid = " << nDip << endl; 
	fplog << "## number of haploid = " << nHap << endl; 
	fplog << "## number of individuals = " << nIndiv  << endl;
	fplog << "## number of snp = " << nLoci  << endl;

	return 1; 
}  

int ModelnData::read_bimbam_phenotype(void)
{
	if(vPin.size() == 0 || vPin.size() != vGin.size()) 
		return 0; 
	{
		fstream infile; 
		streambuf * pbuf; 
		string sfn; 
		char delimit[] = ",; :\t";
			   
		for (int nf = 0; nf < (int) vPin.size(); nf++)
		{   
			int ni = vFileIndiv.at(nf);
			sfn.assign(vPin.at(nf)); 
			
//			cout << "length of file name " << sfn.length() << endl; 
            if(sfn.length() < 5 && isdigit(sfn.c_str()[0])) 
			{
				int ph = atoi(sfn.c_str()); 
				for (int i = 0; i < ni; i++)
					m_phval.push_back(ph);
			}
			else 
			{
				infile.open(sfn.c_str(), ios::in);
				if(!infile.is_open())
				{
					fplog << "ERROR: cannot open phenotype file " << nf << endl; 
					cout << "ERROR: cannot open phenotype file " << nf << endl; 
					safe_exit(); 
				}
				pbuf = infile.rdbuf();
				string line;
				
//				cout << " ... in read_bimbam_phenotype >> " << ni << endl; 
				for (int i = 0; i < ni; i++)
				{
					line.assign(getline(pbuf)); 
					char * res = strtok((char*)line.c_str(), delimit); 
					if(res == NULL) break; 
					m_phval.push_back(atoi(res));
					if(vFilePloid.at(nf) == 1) 
						m_phval.push_back(atoi(res));
				}
				infile.close(); 
			}
		}
	}
	
	cout << "-warning: total number of individuals = " << m_phval.size() << endl; 
	fplog << "### total number of individuals " << m_phval.size() << endl; 
//	fplog << "### " ; 
//	for (unsigned i = 0; i < m_phval.size(); i++)
//		fplog << m_phval.at(i) << " "; 
//	fplog << endl; 

	{
		nCohort = 0; 
		nPanel = 0; 
	   	for (int i = 0; i < nIndiv; i++)
		{
			int ph = m_phval.at(i); 
			if(ph >= 10 || ph == 0) 
			{
				nPanel++; 
				pIndiv[i]->SetisPanel(1); 
			}
			else if(ph != 9) 
			{
				pIndiv[i]->SetisPanel(0); 
				nCohort++; 
			}

			if(ph == 9) 
	   			pIndiv[i]->set_mask(); 
		}
	}
        
	double * vn = new double[nS]; 
	for (int s = 0; s < nS; s++)
		vn[s] = 0; 
    for (int i = 0; i < nIndiv; i++)
	{
		int ph = (int) (m_phval.at(i)); 
		if(ph >= 10) 
			vn[ph-10] += 1.0; 
	}               

	int has0 = 0; 
	for (int s = 0; s < nS; s++)
	{
	    if((int) vn[s] == 0) 
			has0 = 1;
	}


    for (int i = 0; i < nIndiv; i++)
	{
        if(pIndiv[i]->GetisPanel() == 1) 
		{
			int ph = (int) (m_phval.at(i)); 
			if(ph >= 10) {
				if(has0) 
					pIndiv[i]->set_weight(1./vn[ph-10] * nCohort ); 
				else 
					pIndiv[i]->set_weight(25./vn[ph-10] * nCohort ); 
			}
			if(ph == 0) {
				pIndiv[i]->set_weight(10./nPanel * nCohort); }

		}
		else 
		   	pIndiv[i]->set_weight(1.0); 

		if(pIndiv[i]->get_mask()==1)
			pIndiv[i]->set_weight(0); 
	}

    for (int i = 0; i < nIndiv; i++)
		cout <<  pIndiv[i]->get_weight() << " "; 
	cout << endl; 

	cout << "nPanel = " << nPanel << " ## " << "nCohort = " << nCohort << endl; 
	fplog << "## number of panel individuals = " << nPanel << endl; 
	fplog << "## number of cohort = " << nCohort << endl; 
	
	delete[] vn; 
	return 1; 
}

int ModelnData::read_rs2pos(string fn, long start_pos, long end_pos, long flank) 
{
	ifstream infile; 
	char delimit[] = ";, \t";
	streambuf * pbuf;
	int snpcount = 0; 
	int count = 0; 

    map<string, long> :: iterator liter; 
    map<string, int> :: iterator iter; 
	
	//first, read in the rs pos file; and put into vRsPos;
	infile.open(fn.c_str(), ios::in); 
	if(!infile.is_open()) 
	{
		return 0; 
	}
	pbuf = infile.rdbuf(); 	
	string line; 
	line.assign(getline(pbuf)); 
	while (!line.empty()) 
	{
		count++; 
//		if( count % 1 == 0) 
		{
			line.append(" 0 "); 
			char * res = strtok((char *)line.c_str(), delimit); 
			string rs(res); 
			res = strtok(NULL, delimit); 
			long pos = atol(res);
			res = strtok(NULL, delimit); 
			int chr = atoi(res);
			
			if (end_pos > 0) 
			{
				if(pos > start_pos - flank && pos < end_pos + flank)
				{
					liter = mapRs2pos.find(rs); 
					if(liter == mapRs2pos.end())
						mapRs2pos.insert(pair<string, long> (rs, (long)pos));
					else 
						liter->second = pos; 
					
					iter = mapRs2chr.find(rs); 
					if(iter == mapRs2chr.end())
						mapRs2chr.insert(pair<string, int> (rs, chr)); 
					else 
						iter->second = chr; 
					snpcount++;
				}
			}
			else 
			{
				mapRs2pos.insert(pair<string, long> (rs, pos));
				mapRs2chr.insert(pair<string, int> (rs, chr)); 
				snpcount++;
			}
		}
		line.assign(getline(pbuf)); 
	}
	infile.close();
	return snpcount; 
}

int ModelnData::read_rs2map(string fn) 
{
	ifstream infile; 
	char delimit[] = ";, \t";
	streambuf * pbuf;
	int snpcount = 0; 

	map<string, double> mapRs2gm; 
	morgan = new double[nLoci]; 

	cout << "read rs2map "; 
	cout << nLoci << endl; 
	
	infile.open(fn.c_str(), ios::in); 
	if(!infile.is_open()) 
	{
		return 0; 
		exit(0); 
	}
	pbuf = infile.rdbuf(); 	
	string line; 
	line.assign(getline(pbuf)); 
	while (!line.empty()) 
	{
		char * res = strtok((char *)line.c_str(), delimit); 
		string rs(res); 
		res = strtok(NULL, delimit); 
		res = strtok(NULL, delimit); 
		double morgan = atof(res);
		
		mapRs2gm.insert(pair<string, double> (rs, morgan)); 
		snpcount++;
		line.assign(getline(pbuf)); 
	}
	infile.close();

    for (int m = 0; m < nLoci; m++)
	{
		string rs; 
		rs.assign(vsRsnum.at(m)); 
    	map<string, double> :: iterator iter; 
		iter = mapRs2gm.find(rs); 
		if(iter == mapRs2gm.end())
		{
			cout << "wrong map " << endl; 
			exit(0); 
		}
		else 
			morgan[m] = iter->second; 
	}

	for (int m = nLoci-1; m > 0; m--)
	{
		double temp = (morgan[m] - morgan[m-1]) * 5000.0; 
		morgan[m] = 1.0 - exp(-temp); 
	}
	morgan[0] = 1.0; 

//	for (int m = 0; m < nLoci; m++)
//		cout << vsRsnum.at(m) << " \t " << morgan[m] << endl; 

	return snpcount; 
}

void ModelnData::write_admix()
{
	fstream outfile; 
	string fn("output/");
	fn.append(fnOutput);
	fn.append(".admix.txt");
	outfile.open(fn.c_str(), ios::out);
	if(!outfile.good()) 
	{
		fplog << "ERROR: failed to open file to write." << endl;
		cout << "-warning: failed to open file to write" << endl;
		return; 
	}

//	for (int i = 0; i < nK; i++)
//	{
//		if(pHap == NULL) break; 
//		double * admix = pHap[i]->get_admix(); 
//		if(admix == NULL) continue; 
//		char buf[100]; 
//		for (int s = 0; s < nS; s++)
//		{
//			sprintf(buf, "%5.3f ", admix[s]); 
//			outfile << buf;  
//		}
//		outfile << endl; 
//	}

	for(int i = 0; i < nIndiv; i++)
	{
		if(pIndiv[i]->get_mask() == 1 || pIndiv[i]->GetisPanel()) continue; 
		double * admix = pIndiv[i]->get_admix(); 
		if(admix == NULL) continue; 
		for (int s = 0; s < nS; s++)
		{
			char buf[100]; 
			sprintf(buf, "%5.3f ", admix[s]); 
			outfile << buf;
		}
		outfile << endl; 
	}   //print admix; 
	outfile << endl; 
}

void ModelnData::read_admix(string sfn)
{
	fstream infile;
	streambuf * pbuf;
	infile.open(sfn.c_str(), ios::in);
	if(!infile.is_open()) 
	{
		cout << "-warning: cannot open file to read:" << sfn << endl;
		return;
	}
	pbuf = infile.rdbuf(); 	
	char delimit[] = ",; \t:"; 

	string line; 
	for (int k = 0; k < nK; k++)
	{
		line.assign(getline(pbuf)); 
		char * res = strtok((char*)line.c_str(), delimit); 
		for (int s = 0; s < nS; s++)
		{
			if(res == NULL) 
			{
				cout << "-warning: bad file in read admix" << endl; 
				safe_exit(); 
			}
			double tt = atof(res);
			pHap[k]->set_admix(s, tt, nS); 
			res = strtok(NULL, delimit); 
		}
	}

	for(int i = 0; i < nIndiv; i++)
	{
		line.assign(getline(pbuf)); 
		char * res = strtok((char*)line.c_str(), delimit); 
		for (int s = 0; s < nS; s++)
		{
			if(res == NULL) 
			{
				cout << "-warning: bad file in read admix" << endl; 
				safe_exit(); 
			}
			double tt = atof(res);
			pIndiv[i]->set_admix(s, tt, nS); 
			res = strtok(NULL, delimit); 
		}
	}   
}

void ModelnData::write_genotype(int type, int letter)   
//type = 0 exact; 
//type = 1 mean; 
//type = 2 best guess; 
//type = 3 distribution;
{
	fstream outfile; 
	string fn("output/");
	if (type == 0)
	{
		fn.append(fnOutput);
		if(letter == 0) 
			fn.append(".exact.genotype.012");
		else 
			fn.append(".exact.genotype.txt");
		outfile.open(fn.c_str(), ios::out);
		if(!outfile.good()) 
		{
			fplog << "ERROR: failed to open file to write." << endl;
			cout << "-warning: failed to open file to write" << endl;
			return; 
		}
		
		if (letter == 0)   //write in numeric format.   
		{
			for(int i = 0; i < nIndiv; i++)
			{   
//				string rs(vsRsnum.at(m)); 
//				outfile << rs << " " << mapRs2mm[rs].first << " " << mapRs2mm[rs].second << " ";  
				if(pIndiv[i]->get_mask() == 1 || pIndiv[i]->GetisPanel()) continue; 
				for (int m = 0; m < nLoci; m++)
				{
					string rs(vsRsnum.at(m)); 
					short gt = 0; 
					gt = pIndiv[i]->GetsnpGT(m);
					
					if(gt == QQ)
						outfile << "9 ";
					else 
						outfile << gt << " "; 
				}
				outfile << endl; 
			}					
		}
		else           //write bimbam letter format. 
		{   
			int ni = 0; 
			for (int i = 0; i < nIndiv; i++)
				if(pIndiv[i]->get_mask() == 0) ni++; 
			outfile << ni << endl; 
			for (int m = 0; m < nLoci; m++)
			{
				string rs(vsRsnum.at(m)); 
				outfile << rs << " "; 
				for (int i = 0; i < nIndiv; i++)
				{
					short gt = 0; 
					if(pIndiv[i]->get_mask() == 1 || pIndiv[i]->GetisPanel()) continue; 
					gt = pIndiv[i]->GetsnpGT(m);
					if(gt == QQ) 
						outfile << " NN"; 
					else 
					{
						if(gt == 0)  outfile << " " << mapRs2mm[rs].first << mapRs2mm[rs].first; 
						else if(gt == 1)  outfile << " " << mapRs2mm[rs].first << mapRs2mm[rs].second; 
						else  outfile << " " << mapRs2mm[rs].second << mapRs2mm[rs].second; 
					}
				}
				outfile << endl; 
			}
		}
		outfile.close(); 
		return ; 
	}

	outfile.close(); 
}

void ModelnData::write_em(int d)
{
	fstream outfile; 
	string fn("output/");
	fn.append(fnOutput);
	if(d == 1)
		fn.append(".em.txt");
	else
		fn.append(".em2.txt");
	outfile.open(fn.c_str(), ios::out);
	if(!outfile.good()) 
	{
		fplog << "ERROR: failed to open file to write." << endl;
		cout << "-warning: failed to open file to write" << endl;
		return; 
	}

	char buf[100];
	for (int e = 0; e < nEMRuns; e++)
	{
//		double ** kappa = (pMP+e)->Getkappa(); 
//		for (int h = 0; h < 4; h++)
//		{
//			for (int m = 0; m < nLoci; m++)
//			{
//				sprintf(buf, "%5.3f ", kappa[m][h]); 
//				outfile << buf; 
//			}
//			outfile << endl; 
//		}
//		outfile << endl; 

		double * ru = (pMP+e)->Getru(); 
		for (int m = 0; m < nLoci; m++)
		{
			sprintf(buf, "%8.7f ", ru[m]); 
			outfile << buf; 
		}
		outfile << endl; 
		double * rs = (pMP+e)->Getrs(); 
		for (int m = 0; m < nLoci; m++)
		{
			sprintf(buf, "%8.7f ", rs[m]); 
			outfile << buf; 
		}
		outfile << endl; 
		double * rk = (pMP+e)->Getrk(); 
		{
			for (int m = 0; m < nLoci; m++)
			{
				sprintf(buf, "%8.7f ", rk[m]); 
				outfile << buf; 
			}
			outfile << endl; 
		}
		outfile << endl; 
		
		double ** eta = (pMP+e)->Geteta(); 
		{
			for (int s = 0; s < nS; s++)
			{
				for (int m = 0; m < nLoci; m++)
				{   
					sprintf(buf, "%5.3f ", eta[m][s]); 
					outfile << buf; 
				}
				outfile << endl; 
			}
			outfile << endl; 
		}
		outfile << endl; 

		double ** theta = (pMP+e)->Gettheta(); 
		{
			for (int k = 0; k < nK; k++)
			{
				for (int m = 0; m < nLoci; m++)
				{   
					sprintf(buf, "%5.3f ", theta[m][k]); 
					outfile << buf; 
				}
				outfile << endl; 
			}
			outfile << endl; 
		}
		outfile << endl; 

		double *** beta = (pMP+e)->Getbeta(); 
		for (int s = 0; s < nS; s++)
		{
			for (int k = 0; k < nK; k++) 
			{
				for (int m = 0; m < nLoci; m++)
				{
					sprintf(buf, "%5.3f ", beta[m][s][k]); 
					outfile << buf; 
				}
				outfile << endl; 
			}
			outfile << endl; 
  		}                                               
		outfile << endl; 

	}
	outfile.close(); 

}

int ModelnData::read_em(string sfn)   
{
	fstream infile;
	streambuf * pbuf;
	infile.open(sfn.c_str(), ios::in);
	if(!infile.is_open()) 
	{
		cout << "-warning: cannot open file to read:" << sfn << endl;
		return 0;
	}
	pbuf = infile.rdbuf(); 	
	char delimit[] = ",; \t:"; 

	string line; 

//	line.assign(getline(pbuf));
//	nEMRuns = atoi(line.data()); 
//	nEMRuns = 1; 

//	line.assign(getline(pbuf));
//	int loci = atoi(line.data()); 
//	if(loci != nLoci) 
//	{
//		cout << "-warning: em parameters not consistent" << endl; 
//		exit(0);
//	}

//	line.assign(getline(pbuf));
//	nS = atoi(line.data()); 
//	nS = 2; 

//	line.assign(getline(pbuf));
//	nK = atoi(line.data()); 
//	nK = 5; 

	if (pMP) { delete[] pMP; pMP = NULL; }
	InitModelParam(); 

	for (int runs = 0; runs < nEMRuns; runs++)
	{
		{
			line.assign(getline(pbuf)); 
			char * res = strtok((char*)line.c_str(), delimit); 
			for (int m = 0; m < nLoci; m++)
			{
				if(res == NULL) 
				{
					cout << "-warning: bad file in read em" << endl; 
					safe_exit(); 
				}
				double tr = atof(res);
				if(tr < 1e-6) tr = 1e-6; 
				pMP[runs].Setru(m, tr); 
				res = strtok(NULL, delimit); 
			}
			res = strtok(NULL, delimit); 
			if(res != NULL) {
				cout << "-eali: more entries in em file than in genotype file" << endl; 
				safe_exit(); 
			}
		}

		{
			line.assign(getline(pbuf)); 
			char * res = strtok((char*)line.c_str(), delimit); 
			for (int m = 0; m < nLoci; m++)
			{
				if(res == NULL) 
				{
					cout << "-warning: bad file in read em" << endl; 
					safe_exit(); 
				}
				double tr = atof(res);
				if(tr < 1e-6) tr = 1e-6; 
				pMP[runs].Setrs(m, tr); 
				res = strtok(NULL, delimit); 
			}
		}
		
		{
			line.assign(getline(pbuf)); 
			char * res = strtok((char*)line.c_str(), delimit); 
			for (int m = 0; m < nLoci; m++)
			{
				if(res == NULL) 
				{
					cout << "-warning: bad file in read em" << endl; 
					safe_exit(); 
				}
				double tr = atof(res);
				if(tr < 1e-6) tr = 1e-6; 
				pMP[runs].Setrk(m, tr); 
				res = strtok(NULL, delimit); 
			}
		}
		
		{
			for (int s = 0; s < nS; s++)
			{
				line.assign(getline(pbuf)); 
				char * res = strtok((char*)line.c_str(), delimit); 
				for (int m = 0; m < nLoci; m++)
				{
					if(res == NULL) 
					{
						cout << "-warning: bad file in read em" << endl; 
						safe_exit(); 
					}
					double tt = atof(res);
			   		pMP[runs].Seteta(m, s, tt); 
					res = strtok(NULL, delimit); 
				}
			}
		}

		{
			for (int k = 0; k < nK; k++)
			{
				line.assign(getline(pbuf)); 
				char * res = strtok((char*)line.c_str(), delimit); 
				for (int m = 0; m < nLoci; m++)
				{
					if(res == NULL) 
					{
						cout << "-warning: bad file in read em" << endl; 
						safe_exit(); 
					}
					double tt = atof(res);
			   		pMP[runs].Settheta(m, k, tt); 
					res = strtok(NULL, delimit); 
				}
			}
		}
		
		
		for (int s = 0; s < nS; s++)
			for (int k = 0; k < nK; k++)
			{
				line.assign(getline(pbuf)); 
				char * res = strtok((char*)line.c_str(), delimit); 
				for (int m = 0; m < nLoci; m++)
				{
					if(res == NULL) 
					{
						cout << "-warning: bad file in read em" << endl; 
						safe_exit(); 
					}
					double ta = atof(res);
					pMP[runs].Setbeta(m, s, k, ta); 
					res = strtok(NULL, delimit); 
				}
			}
	}

	return 1; 
}
