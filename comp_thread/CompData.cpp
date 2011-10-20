
#include "CompData.h"

#include <iostream>
#include <sstream>

using namespace boost::numeric::ublas;

Comp::Comp(){
    //
}


void Comp::update(double epsilon, int maxiters) {
	/*-------- 
	Does stuff
	--------*/ 
	std::cout<<"is computing..."<<std::endl;
	this->simplexe_.reset( this->guess_, this->increments_, this->ch_ );
	this->simplexe_.display_ = Display( "test" );
	this->results_ = this->simplexe_.minimize( epsilon, maxiters );
	std::cout<<"has computed."<<std::endl;
	std::cout<<std::endl;
}

void Comp::set_coeffs( double Qt, double Qe ){	
	this->Qt = Qt; this->Qe = Qe;
	this->simplexe_ = Simplex( this->Qt, this->Qe );
	std::cout<<"****************\n\tSimplex Object initialized.\n****************"<<std::endl;
}

std::string Comp::get_message() {
	/*-------- 
	Writes Message
	--------*/ 
    std::ostringstream oss (std::ostringstream::out);
	oss << this->results_.size()<<"\n";
	for( unsigned int i=0; i<this->results_.size(); i++ )
		oss <<results_(i)<<"\n";
    return oss.str();
}

int Comp::analyse_message(std::string msg) {
    /*-------- 
	Parses Received Message
	--------*/ 
    std::string str, s;
	
    std::istringstream iss(msg);
	int count = 0;
    while ( getline(iss, str)){
        std::istringstream in(str);
        in >> s;
		if      (s=="Constants"){	in>>ch_.M_>>
										ch_.zc_>>
										ch_.g_>>
										ch_.dt_>>
										ch_.h_>>
										ch_.mode_>>ch_.kpinit_; count+=7;}
		else if (s=="Jacobian"){		fill_matrix( ch_.J_			, in );count++;}
		else if (s=="dJacobian"){		fill_matrix( ch_.dJ_		, in );count++;}
		else if (s=="JacobianI"){		fill_matrix( ch_.Ji_		, in );count++;}
		else if (s=="dJJi"){			fill_matrix( ch_.dJJi_		, in );count++;}
		else if (s=="JtiHJi"){			fill_matrix( ch_.JtiHJi_	, in );count++;}
		else if (s=="Pref"){			fill_matrix( ch_.Pref_		, in );count++;}
		else if (s=="FDIS"){			fill_matrix( ch_.FDIS_		, in );count++;}
		else if (s=="x"){				fill_vector( ch_.x_			, in );count++;}
		else if (s=="dx"){				fill_vector( ch_.dx_		, in );count++;}
		else if (s=="ddx"){				fill_vector( ch_.ddx_		, in );count++;}
		else if (s=="xc"){				fill_vector( ch_.xc_			, in );count++;}
		else if (s=="dxc"){				fill_vector( ch_.dxc_		, in );count++;}
		else if (s=="ddxc"){				fill_vector( ch_.ddxc_		, in );count++;}
		else if (s=="xdes"){			fill_vector( ch_.xdes_		, in );count++;}
		/* -------- Parameters -------- */
		else if (s=="guess"){			fill_vector( this->guess_		, in );count++;}
		else if (s=="increments"){		fill_vector( this->increments_	, in );count++;}
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

void Comp::fill_vector(ublas::vector<double> &vec, std::istringstream &in){
	/*-------- 
	Fills vector with istringstream data in = [nlin data....]
	--------*/ 
	int nlin;
	in>>nlin;
	vec.resize( nlin );
	for( int lin=0; lin<nlin; lin++ )
			in>>vec(lin);
}

