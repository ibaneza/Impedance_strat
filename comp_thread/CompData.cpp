
#include "CompData.h"

#include <iostream>
#include <sstream>

void Comp::InitSensor() {
	//Init some stuff
}

Comp::Comp(){
    InitSensor();
}

void Comp::Update(float elapsedTime, float LastThinkDt, int bJustThink) {
    if (bJustThink) {
		//DO STUFF
		std::cout<<"is thinking..."<<std::endl;
		std::cout<<std::endl;
    }
}

std::string Comp::get_message() {
	// Write message
    std::ostringstream oss (std::ostringstream::out);
	oss <<"FirstLineOfTwoArgs "<<4.<<" "<<5.<<" \n"
		<<"SecLineOfOneArgs "<<4.<<" \n";
    return oss.str();
}

float Comp::analyse_message(std::string msg) {
    //PARSE MESSAGE
    std::string str, s;
	float receiver;
    std::istringstream iss(msg);
    while ( getline(iss, str)){
        std::istringstream in(str);
        in >> s;
        if      (s=="RecTwoArgs") in>>receiver>>receiver;
        else if (s=="RecOneArg") in>>receiver;
    }
    return receiver;
}


