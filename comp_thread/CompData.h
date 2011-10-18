#ifndef _COMPDATA_HEADER_
#define _COMPDATA_HEADER_

#include <stdio.h>
#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "Simplex.h"

using namespace boost::numeric;

struct Comp{
	
	//Some stuff to fill simplex points
	Constants_holder ch_;
	//...

    Comp();
    void Update(float elapsedTime, float LastThinkDt, int bJustThink);
    std::string get_message();
    int analyse_message(std::string msg);
	void fill_matrix( ublas::matrix< double > &mat, std::istringstream &in );
	void fill_vector( ublas::vector< double > &vec, std::istringstream &in );

};

#endif //_COMPDATA_HEADER_
