#include <iostream>
#include <time.h>
#include "fp.h"
#include "control.h"
#include "fpmath.h"
#include <sys/stat.h>
#include <sys/types.h>

#ifdef WINDOWS
#include <windows.h> 
#endif

//#include <assert.h>
using namespace std;

//global gsl variable in use; 
const gsl_rng_type * gslType;
gsl_rng * gsl_r; 
//should put into a namespace later on. 
//extern was defined in fpmath.h which will be included in every file; 

int main(int argc, char * argv[])
{ 
				
	gsl_rng_env_setup();                                                               
	gslType = gsl_rng_default; 
	gsl_r = gsl_rng_alloc(gslType); 
	gsl_rng_set(gsl_r, (unsigned long) time(NULL)); 
	//init the global gsl; note the seed is only good before call Run(); 
	
	CtrlParam objCtrl;
	//declear of objCtrl should go after gsl init. 
    
	std::ifstream check("output/");
	if(!check) {
#ifndef WINDOWS
	    mkdir("output", S_IRWXU);
#else
	    CreateDirectory("output",NULL);
#endif
	}   //check if dir exist. create if not. 
	
	if(argc <= 1)
		objCtrl.PrintHeader(); 
#if defined (READLINE)
	else if(argc == 2 && argv[1][0] == '-' && argv[1][1] == 'I')
		objCtrl.Run();
#endif

	objCtrl.BatchRun(argc, argv); 
	
    return 1;                                                          
}

 
