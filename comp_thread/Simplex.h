#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include "Simplex_Pt.h"

using namespace boost::numeric;

class Simplex{
public:
	Simplex( double Qt = 0., double Qe=0., double kR=-1., double kE=2., double kC=.5 );
	~Simplex();
	void reset( ublas::vector< double > guess, ublas::vector< double > increments, 
						Constants_holder ch );
	void _reset();
	ublas::vector< double > minimize( double epsilon = 100., int maxiters = 2500 );	

private:
	ublas::vector< double > stats_;
	/* -------- Init-related stuff ---------- */
	ublas::vector< double > increments_;
	int numvars_;
	/* -------- Minimization-related stuff ---------- */
	double kR_, kE_, kC_;
	ublas::vector< Simplex_Pt > simplex_;
	Simplex_Pt guess_;
	int lowest_, highest_, secondhighest_;
	ublas::vector< double > errors_;
	double currenterror_;

private:
	/* -------- Some stuff to fill simplex points ---------- */
	Constants_holder ch_;
	double QtonR, QeonR;				// It's all in the title
	void set_coeff( double Qt, double Qe ){ this->QtonR=Qt; this->QeonR=Qe;};

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