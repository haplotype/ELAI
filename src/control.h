// define ctrl class  
#ifndef __CONTROL_H__
#define __CONTROL_H__

#include <stdlib.h>
#include <ctype.h>

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "fpmath.h"
#include "model.h"

using namespace std;

class CtrlParam
{
private:
	map<string, int> hcom; 
	class ModelnData * pMD; 
	
public: 
	CtrlParam(void); //default construct;
	~CtrlParam(void); //destruct. 

	void BatchRun(int, char **);
	void OnHelp(void);
	void PrintHeader(void);
};

#endif
