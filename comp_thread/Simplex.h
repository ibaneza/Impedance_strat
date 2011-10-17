#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include "Simplex_Pt.h"

using namespace boost::numeric;

class Simplex{
public:
	Simplex( double kR=-1., double kE=2., double kC=.5 );
	~Simplex();
	void reset( ublas::vector< double > guess, ublas::vector< double > increments );
	void _reset();
	ublas::vector< double > minimize( double epsilon = 100., int maxiters = 2500 );
private:
	ublas::vector< double > increments_;
	Simplex_Pt guess_;
	double kR_, kE_, kC_;
	int numvars_;
	ublas::vector< Simplex_Pt > simplex_;
	int lowest_, highest_, secondhighest_;
	ublas::vector< double > errors_;
	double currenterror_;

private:
	//Some stuff to fill simplex points
	int  M, zc, g, dt, h;
	ublas::vector< double > x_, dx_, ddx_;	// Manipulation task	current position, velocity, acceleration
	ublas::vector< double > xdes_;			// Manipulation tas		desired position
	ublas::matrix< double > J_, dJ_, Ji_;	// Some Jacobian		derivates
	ublas::matrix< double > dJJi_, JtiHJi_;	// Some Jacobian		composites
	ublas::matrix< double > Pref_;			// ZMP		position
	ublas::matrix< double > FDIS_;			// Total Disturbance	force Horizon

private:
	void compute_error_at_vertices();
	void contract_simplex();
	void expand_simplex();
	void reflect_simplex();
	void multiple_contract_simplex();
	void accept_contracted_point();
	void accept_expanded_point();
	void accept_reflected_point();
};


#endif //_SIMPLEX_H_