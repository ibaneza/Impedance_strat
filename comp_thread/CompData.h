#ifndef _COMPDATA_HEADER_
#define _COMPDATA_HEADER_

#include <stdio.h>
#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace boost::numeric;

struct Comp{
	
	//Some stuff to fill simplex points
	int  M, zc, g, dt, h;
	ublas::vector< double > x_, dx_, ddx_;	// Manipulation task	current position, velocity, acceleration
	ublas::vector< double > xdes_;			// Manipulation tas		desired position
	ublas::matrix< double > J_, dJ_, Ji_;	// Some Jacobian		derivates
	ublas::matrix< double > dJJi_, JtiHJi_;	// Some Jacobian		composites
	ublas::matrix< double > Pref_;			// ZMP		position
	ublas::matrix< double > FDIS_;			// Total Disturbance	force Horizon
	//...

    Comp();
    void Update(float elapsedTime, float LastThinkDt, int bJustThink);
    std::string get_message();
    float analyse_message(std::string msg);
	void fill_matrix( ublas::matrix< double > &mat, std::istringstream &in );
	void fill_vector( ublas::vector< double > &vec, std::istringstream &in );

};

#endif //_COMPDATA_HEADER_
