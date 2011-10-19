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
	Simplex simplexe_;
	/* -------- Some stuff to fill simplex points ---------- */
	Constants_holder ch_;
	double Qt, Qe;
	ublas::vector< double > guess_;
	ublas::vector< double > increments_;
	/* -------- Results ---------- */
	ublas::vector< double > results_;

    Comp();
    void update(double epsilon = 1., int maxiters = 500);
    std::string get_message();
    int analyse_message(std::string msg);
	void fill_matrix( ublas::matrix< double > &mat, std::istringstream &in );
	void fill_vector( ublas::vector< double > &vec, std::istringstream &in );
	void set_coeffs( double Qt, double Qe );

};

#endif //_COMPDATA_HEADER_
