#include <stdio.h>

#include "SpiropsCom.h"
#include "CompData.h"

Comp   gPerso;

int main(int argc, char** argv) {

    /*---------
    /Before, we load arguments...
    /--------*/
    char hostname[128];
    int port;
    bool verbose=false;
    for (int i=0; i<argc; i++){
        if (!strncmp(argv[i], "-socket", 7)){
            strcpy(hostname, argv[i+1]);
            port = atoi(argv[i+2]);
        }
        else if (!strncmp(argv[i], "-verbose", 8)) {
            verbose=true;
        }
    }

    /*-------- Init program value ----------*/
    std::string msg;
    SpiropsCom SC(hostname, port);

    /*------- Let's rock!!! ----------*/
    while(SC.isConnected()) {
        msg = SC.recv_message(verbose);
        gPerso.analyse_message(msg);
        msg = gPerso.get_message();
        SC.send_message(msg);
    }

    system("pause");
    return 0;
}


