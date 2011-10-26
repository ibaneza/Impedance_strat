#ifndef _SPIROPSCOM_H_
#define _SPIROPSCOM_H_

//#include <windows.h>
#include <winsock2.h>
#pragma comment(lib, "ws2_32.lib")

#include <iostream>
#include <sstream>


class SpiropsCom
{
public:
    SpiropsCom(const char*, const int _port);
    ~SpiropsCom();
    std::string recv_message(bool);
    void send_message(std::string);
    bool isConnected();

private:
    SOCKET sock;
    bool _isConnected;
};

SpiropsCom::SpiropsCom(const char* _host="127.0.0.1", const int _port=5556)
{
    SOCKADDR_IN sin;

    WSADATA WSAData;
    WSAStartup(MAKEWORD(2,0), &WSAData);

    sin.sin_addr.s_addr = inet_addr(_host);
    sin.sin_family      = AF_INET;
    sin.sin_port        = htons(_port);

    sock = socket(AF_INET,SOCK_STREAM,0);

    _isConnected = true;
    if (connect(sock, (SOCKADDR *)&sin, sizeof(sin)) != 0) {
        printf("Connection problem!!!\n");
        _isConnected = false;
    };
    /*
    u_long nonBlockingMode = 1;
    ioctlsocket(sock, FIONBIO, &nonBlockingMode);
    */
}

SpiropsCom::~SpiropsCom()
{
    closesocket(sock);
    WSACleanup();
}


#define SPIROPSCOM_MAX_LEN 32768 /// WARNING: it doesn't work with 4096
std::string SpiropsCom::recv_message(bool verbose)
{

    char msg[SPIROPSCOM_MAX_LEN];
    memset(msg, 0x0, SPIROPSCOM_MAX_LEN);
    if (verbose) std::cout<<"reading message...\n";
    int numRead = recv(sock, msg, SPIROPSCOM_MAX_LEN, 0);

    std::string smsg(msg);
    if (verbose) {
    std::cout<<"== "<<numRead<<" ===========================================\n"
             <<smsg<<std::endl
             <<"-----------------------------------------------------\n"<<std::endl;
    }

    if (numRead < 0){
        int numError = WSAGetLastError();
        std::cout<<"the error number is: "<<numError<<std::endl;
        _isConnected = false;
    }

    if (smsg=="close connection") {
        _isConnected = false;
    }
    return smsg;
}


void SpiropsCom::send_message(std::string msg)
{
    char sent_msg[SPIROPSCOM_MAX_LEN];
    memset(sent_msg, '\0', SPIROPSCOM_MAX_LEN);
    //memcpy(sent_msg,msg.c_str(), msg.len()+1);
    strcpy(sent_msg,msg.c_str());
    send(sock, sent_msg, SPIROPSCOM_MAX_LEN, 0);
}


bool SpiropsCom::isConnected() {
    return _isConnected;
}

#endif
