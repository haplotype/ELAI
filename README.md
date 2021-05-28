1 Copyright
ELAI — Efficient local ancestry inference. Copyright (C) 2014–2016 Yongtao Guan.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/.


2. What ELAI Can Do
The software, Efficient Local Ancestry Inference (ELAI), is developed and maintained by Yong- tao Guan (http://www.ncbi.nlm.nih.gov/pubmed/24388880). Please refer to the paper for the details of the statistical method. The ELAI is designed to perform local ancestry inference for ad- mixed individuals. Comparing to existing methods to infer local ancestry, ELAI has the following advantages:

   1) Directly works with diploid data, no phasing required,
   2) No recombination map is required—the recombination rates are implicitly estimated,
   3) It has a high resolution and can detect local ancestry track length of a few tenth of a centi- Morgan (cM).
   4) The new version assign weights to training at cohort samples so that it applies to datasets of an arbitrary number of individuals (the previous version requires splitting a large sample into small subsets). The documentation of weighting scheme in EM updates can be found at (http://www.ncbi.nlm.nih.gov/ pubmed/26863142). 

... (see manual) 

8. ELAI options. 
Unless otherwise stated, arg stands for a string, num stands for a number. 

File I/O related options:
• -g arg     can use multiple times, must pair with -p. 
• -p arg     can use multiple times, must pair with -g. arg takes integer values, 1, 10, 11, 12, . . ..
• -pos arg   can use multiple times. arg is a file name.
• -o arg     arg will be the prefix of all output files, the random seed will be used by default.

EM Parameters:
• -s(step) num    specify steps in EM run.
• -C num          specify number of upper clusters. 
• -c num          specify number of lower clusters.
• -mg num         specify number of mixture generations.
• -R num          specify random seed, system time by default.
• -sem num        save EM results to prefix.em.txt.
• -rem file       read EM from a file. 
 
Other options:
• -v(ver)             print version and citation
• -h(help)            print this help
• -exclude-maf num    exclude SNPs whose maf is less than num, default 0.
• --exclude-nopos     exclude SNPs that has no position information. 
• --exclude-miss1     exclude SNPs that are missing in at least one file.
• --silence           no terminal output.    
• --ps2               output joint distribution of inferred ancestry in file pref.ps22.txt. [05/13/2021]
