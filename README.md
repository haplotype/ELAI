1. Copyright

ELAI — Efficient local ancestry inference. Copyright (C) 2014–2021 Yongtao Guan.

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

The software, Efficient Local Ancestry Inference (ELAI), is developed and maintained by Yongtao Guan (http://www.ncbi.nlm.nih.gov/pubmed/24388880). Please refer to the paper for the details of the statistical method. The ELAI is designed to perform local ancestry inference for ad- mixed individuals. Comparing to existing methods to infer local ancestry, ELAI has the following advantages:

   1) Can directly works with diploid data, wiht no phasing required; Can also work we phased haplotypes. 
   2) No recombination map is required and the recombination rates are implicitly estimated. 
   3) Can detect local ancestry track length of a few tenth of a centi-Morgan.
   4) The new version assign weights to training and cohort samples so that ELAI can work with an arbitrarily large number of cohort samples, and allow absent of one training population (http://www.ncbi.nlm.nih.gov/ pubmed/26863142). 

3. Exectuables. 

elai-mac was last compiled on MacOS Big Sur (version 11.3).
elai-lin was last compile on Ubuntu 18.04. 

See manual on commands and examples.  

