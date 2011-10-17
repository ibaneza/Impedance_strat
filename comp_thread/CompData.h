#ifndef _COMPDATA_HEADER_
#define _COMPDATA_HEADER_

#include <stdio.h>
#include <iostream>

struct Comp{
	
	//Some stuff
	
	//...

    Comp();
    void Update(float elapsedTime, float LastThinkDt, int bJustThink);
    std::string get_message();
    float analyse_message(std::string msg);

};

#endif //_COMPDATA_HEADER_
