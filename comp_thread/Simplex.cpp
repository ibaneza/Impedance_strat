#include "Simplex.h"

#include <stdlib.h>

using namespace boost::numeric::ublas;

Simplex::Simplex(double Qt, double Qe, double kR, double kE, double kC):kR_(kR),kE_(kE),kC_(kC){
	this->set_coeff( Qt, Qe );
}

Simplex::~Simplex(){}

void Simplex::_reset(){
	/*--------
	Resets constants
	--------*/
	lowest_ = -1;
	highest_ = -1;
	secondhighest_ = -1;
	currenterror_ = 0.;
}

void Simplex::reset( vector< double > guess, vector< double > increments, Constants_holder ch ){
	/*--------
	Initializes Simplex
	--------*/
	this->_reset();
	this->ch_ = ch;
	this->numvars_ = guess.size();
	// /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\ 
	//Probably more to do to correctly initialize a Simplex_Pt object
	this->guess_.set_coeff( this->QtonR, this->QeonR );
	this->guess_.set_data( guess );
	if( !this->guess_.init( this->ch_, this->numvars_ ) ){
		std::cout<<"ERROR! SIMPLEX POINT NOT CORRECTLY INITIALIZED!"<<std::endl;
		return;
	}
	// /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\ 
	this->increments_.resize(increments.size());
	for( unsigned int i=0; i<increments.size(); i++ ) increments_(i) = increments(i);
	this->errors_.resize( this->numvars_ + 1 );		// Simplex needs N+1 vertices
	/* -------- Initializing Vertices -------- */
	this->simplex_.resize( this->numvars_ + 3 );	// And two extras to store centroid and reflected

	// /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\ 
	//Probably more to do to correctly initialize a Simplex_Pt object
	for( int i=0; i<this->numvars_ + 3; i++ ){
		this->simplex_(i).set_data( guess );
		this->simplex_(i).set_coeff( this->QtonR, this->QeonR );
		this->simplex_(i).init( this->ch_, this->numvars_ );
	}
	// /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\ 

	ublas::vector< double > randv ( this->numvars_ );
	for( int vertex=0; vertex<this->numvars_+1; vertex++ ){
		double maxrand = -1.;
		for( int param=0; param<this->numvars_; param++ ){
			double random = 2*(double)rand() / RAND_MAX - 1.;
			randv(param) = random;
			maxrand = (random>maxrand)?random:maxrand;
		}
		if( vertex!=0 )		//Let's keep the first one as guessed
			this->simplex_(vertex).set_data( guess + (1./maxrand) * element_prod( randv,increments ) );
	}
	/* -------- Distance Computation -------- */
	this->compute_error_at_vertices();
}

vector< double > Simplex::minimize( double epsilon, int maxiters ){
	/*--------
	Minimizes: simplex rolls downhill
	--------*/
	double ipct = 0., ipct_p = 0.;
	for( int iter=0; iter<maxiters; iter++ ){
		/* -------- Identifying highest, second highest, and lowest vertices -------- */
		this->highest_ = 0;
		this->lowest_ = 0;
		for( int vertex=0; vertex<this->numvars_+1; vertex++ ){
			if( this->errors_(vertex) > this->errors_(this->highest_) )
				this->highest_ = vertex;
			if( this->errors_(vertex) < this->errors_(this->lowest_) )
				this->lowest_ = vertex;
		}
		this->secondhighest_ = 0;
		for( int vertex=0; vertex<this->numvars_+1; vertex++){
			if( vertex==this->highest_ )
				continue;
			if( this->errors_(vertex) > this->errors_(this->secondhighest_) )
				this->secondhighest_ = vertex;
		}
		/* -------- Testing for convergence -------- */
		double S=0., S1=0.;
		for( int vertex=0; vertex<this->numvars_+1; vertex++ )
			S += this->errors_(vertex);
		double F2 = S / (this->numvars_ + 1. );
		for( int vertex=0; vertex<this->numvars_+1; vertex++ )
			S1 += pow( this->errors_(vertex) - F2, 2. );
		double T = sqrt( S1 / this->numvars_ );
		if( T <= epsilon )	// We reached convergence, let's get outta here
			break;			
		else{				// We got more things to do then...
			/* -------- Computing centroid - except from highest -------- */
			ublas::vector< double > centroid_data;
			centroid_data.resize( this->numvars_ );
			for( int param=0; param<this->numvars_; param++ ){
				S = 0.;
				for( int vertex=0; vertex<this->numvars_+1; vertex++ ){
					if( vertex==this->highest_ )
						continue;
					S += (this->simplex_(vertex).get_data())(param);
				}
				centroid_data(param) = S / this->numvars_;
			}
			this->simplex_(this->numvars_+1).set_data( centroid_data );
			/* -------- Trying transformations -------- */
			/* -------- Starting with reflection -------- */
			this->reflect_simplex();
			this->currenterror_ = this->guess_.func();
			if( this->currenterror_ < this->errors_(this->lowest_) ){ //Then this is Good
				double tmp = this->currenterror_;
				/* -------- Trying expansion -------- */
				this->expand_simplex();
				this->currenterror_ = this->guess_.func();
				if( this->currenterror_ < tmp )
					this->accept_expanded_point();
				else{
					this->currenterror_ = tmp;
					this->accept_reflected_point();
				}
			}
			else if( this->currenterror_ <= this->errors_(this->secondhighest_) )	//Then this is OK
				this->accept_reflected_point();
			else if( this->currenterror_ <= this->errors_(this->highest_) ){		//Then this is OK
				this->accept_reflected_point();
				/* -------- Trying contraction -------- */
				this->contract_simplex();
				this->currenterror_ = this->guess_.func();
				if( this->currenterror_ < this->errors_(this->highest_) )
					this->accept_contracted_point();
				else
					this->multiple_contract_simplex();
			}
			else if( this->currenterror_ >= this->errors_(this->highest_) ){
				/* -------- Trying contraction -------- */
				this->contract_simplex();
				this->currenterror_ = this->guess_.func();
				if( this->currenterror_ < this->errors_(this->highest_) )
					this->accept_contracted_point();
				else
					this->multiple_contract_simplex();
			}
		}
		modf( 10*(double)iter/(double)maxiters, &ipct );
		if( ipct > ipct_p ) std::cout<<10*ipct<<"% >> error = "<<this->errors_(this->lowest_)<<std::endl;
		ipct_p = ipct;
	}
	return this->simplex_(this->lowest_).get_data();
}

void Simplex::compute_error_at_vertices(){
	for( int vertex=0; vertex<this->numvars_+1; vertex++ ){
		if( vertex==this->lowest_ )
			continue;
		this->guess_.set_data( this->simplex_(vertex).get_data() );
		this->currenterror_ = this->guess_.func();
		this->errors_(vertex) = this->currenterror_;
	}
}

void Simplex::contract_simplex(){
	ublas::vector< double, bounded_array<double, BOUNDED_ARRAY_MAX_SIZE> > guess_data( this->numvars_ );
	for( int param=0; param<this->numvars_; param++ )
		guess_data(param) = this->kC_ * this->simplex_(this->highest_).get_data()(param)
			+ (1. - this->kC_) * this->simplex_(this->numvars_ + 1).get_data()(param);
	this->guess_.set_data( guess_data );
}

void Simplex::expand_simplex(){
	ublas::vector< double, bounded_array<double, BOUNDED_ARRAY_MAX_SIZE> > guess_data( this->numvars_ );
	for( int param=0; param<this->numvars_; param++ )
		guess_data(param) = this->kE_ * this->guess_.get_data()(param) 
			+ (1. - this->kE_) * this->simplex_(this->numvars_ + 1).get_data()(param);
	this->guess_.set_data( guess_data );
}

void Simplex::reflect_simplex(){
	ublas::vector< double, bounded_array<double, BOUNDED_ARRAY_MAX_SIZE> > guess_data( this->numvars_ );
	for( int param=0; param<this->numvars_; param++ )
		guess_data(param) = this->kR_ * this->simplex_(this->highest_).get_data()(param)
			+ (1. - this->kR_) * this->simplex_(this->numvars_ + 1).get_data()(param);
	this->guess_.set_data( guess_data );
	this->simplex_(this->numvars_+2).set_data( guess_data ); //We need to store it, in case of acceptance
}

void Simplex::multiple_contract_simplex(){
	ublas::vector< double, bounded_array<double, BOUNDED_ARRAY_MAX_SIZE> > data( this->numvars_ );
	for( int vertex=0; vertex<this->numvars_+1; vertex++ ){
		if( vertex == this->lowest_ )
			continue;
		for( int param=0; param<this->numvars_; param++ )
			data(param) = .5 * ( this->simplex_(vertex).get_data()(param) 
				+ this->simplex_(this->lowest_).get_data()(param) );
		this->simplex_(vertex).set_data( data );
	}
	this->compute_error_at_vertices();
}

void Simplex::accept_contracted_point(){
	ublas::vector< double, bounded_array<double, BOUNDED_ARRAY_MAX_SIZE> > data( this->numvars_ );
	this->errors_(this->highest_) = this->currenterror_;
	this->simplex_(this->highest_).set_data( this->guess_.get_data() );
}

void Simplex::accept_expanded_point(){
	ublas::vector< double, bounded_array<double, BOUNDED_ARRAY_MAX_SIZE> > data( this->numvars_ );
	this->errors_(this->highest_) = this->currenterror_;
	this->simplex_(this->highest_).set_data( this->guess_.get_data() );
}

void Simplex::accept_reflected_point(){
	ublas::vector< double, bounded_array<double, BOUNDED_ARRAY_MAX_SIZE> > data( this->numvars_ );
	this->errors_(this->highest_) = this->currenterror_;
	this->simplex_(this->highest_).set_data( this->simplex_(this->numvars_+2).get_data() );
}
