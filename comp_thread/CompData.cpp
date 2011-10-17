
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

float Comp::analyse_message(std::string msg) {
    //PARSE MESSAGE
    std::string str, s;
	
    std::istringstream iss(msg);
    while ( getline(iss, str)){
        std::istringstream in(str);
        in >> s;
        if      (s=="Constants") in>>M>>zc>>g>>dt>>h;
		else if (s=="Jacobian") fill_matrix( J_, in );
		else if (s=="dJacobian") fill_matrix( dJ_, in );
		else if (s=="JacobianI") fill_matrix( Ji_, in );
		else if (s=="dJJi") fill_matrix( dJJi_, in );
		else if (s=="JtiHJi") fill_matrix( JtiHJi_, in );
		else if (s=="Pref") fill_matrix( Pref_, in );
		else if (s=="FDIS") fill_matrix( FDIS_, in );
		else if (s=="x") fill_vector( x_, in );
		else if (s=="dx") fill_vector( dx_, in );
		else if (s=="ddx") fill_vector( ddx_, in );
		else if (s=="xdes") fill_vector( xdes_, in );
    }
    return 1.;
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

