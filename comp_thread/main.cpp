#include <stdio.h>

#include "SpiropsCom.h"
#include "CompData.h"
#include "Simplex_Pt.h"

Comp   gPerso;

int main(int argc, char** argv) {

    /*---------
    /Before, we load arguments...
    /--------*/
    char hostname[128];
    int port;
	double Qt,Qe;
	double epsilon;
	int maxiters;
    bool verbose=false;
    for (int i=0; i<argc; i++){
        if (!strncmp(argv[i], "-socket", 7)){
            strcpy(hostname, argv[i+1]);
            port = atoi(argv[i+2]);
        }
        else if (!strncmp(argv[i], "-verbose", 8)) {
            verbose=true;
        }
		else if (!strncmp(argv[i], "-qt", 8)) {
            Qt = atof( argv[i+1] );
        }
		else if (!strncmp(argv[i], "-qe", 8)) {
            Qe = atof( argv[i+1] );
        }
		else if (!strncmp(argv[i], "-epsilon", 8)) {
            epsilon = atof( argv[i+1] );
        }
		else if (!strncmp(argv[i], "-maxiters", 8)) {
            maxiters = atoi( argv[i+1] );
        }
    }

    /*-------- Init program value ----------*/
    std::string msg;
    SpiropsCom SC(hostname, port);
	gPerso.set_coeffs( Qt, Qe );

    /*------- Let's rock!!! ----------*/
    while(SC.isConnected()) {
		if( verbose ) std::cout<<"C++>>RECEIVING"<<std::endl;
        msg = SC.recv_message(verbose);
		if( verbose ) std::cout<<"C++>>PARSING"<<std::endl;
		if( gPerso.analyse_message(msg) ){
			gPerso.update(epsilon,maxiters);
			if( verbose ) std::cout<<"C++>>SENDING"<<std::endl;
			msg = gPerso.get_message();
		}
		else
			if( verbose ) msg = "COMP_THREAD_ERROR";
        SC.send_message(msg);
    }

    system("pause");
    return 0;
}


