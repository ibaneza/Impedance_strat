
#include "CompData.h"

#include <iostream>
#include <sstream>

using namespace boost::numeric::ublas;

Comp::Comp(){
    //
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

int Comp::analyse_message(std::string msg) {
    //PARSE MESSAGE
    std::string str, s;
	
    std::istringstream iss(msg);
    while ( getline(iss, str)){
        std::istringstream in(str);
        in >> s;
		if      (s=="Constants") in>>ch_.M_>>ch_.zc_>>ch_.g_>>ch_.dt_>>ch_.h_;
		else if (s=="Jacobian") fill_matrix( ch_.J_, in );
		else if (s=="dJacobian") fill_matrix( ch_.dJ_, in );
		else if (s=="JacobianI") fill_matrix( ch_.Ji_, in );
		else if (s=="dJJi") fill_matrix( ch_.dJJi_, in );
		else if (s=="JtiHJi") fill_matrix( ch_.JtiHJi_, in );
		else if (s=="Pref") fill_matrix( ch_.Pref_, in );
		else if (s=="FDIS") fill_matrix( ch_.FDIS_, in );
		else if (s=="x") fill_vector( ch_.x_, in );
		else if (s=="dx") fill_vector( ch_.dx_, in );
		else if (s=="ddx") fill_vector( ch_.ddx_, in );
		else if (s=="xdes") fill_vector( ch_.xdes_, in );
    }
	return this->ch_.assume_state();
}

void Comp::fill_matrix(matrix<double> &mat, std::istringstream &in){
	/*-------- 
	Fills matrix with istringstream data in = [nlin ncol data....]
	--------*/ 
	int nlin, ncol;
	in>>nlin>>ncol;
	mat.resize( nlin, ncol );
	for( int lin=0; lin<nlin; lin++ ){
		for( int col=0; col<ncol; col++ ){
			in>>mat(lin,col);
		}
	}
}

void Comp::fill_vector(vector<double> &vec, std::istringstream &in){
	/*-------- 
	Fills vector with istringstream data in = [nlin data....]
	--------*/ 
	int nlin;
	in>>nlin;
	vec.resize( nlin );
	for( int lin=0; lin<nlin; lin++ )
			in>>vec(lin);
}

