#include "Simplex_Pt.h"

using namespace boost::numeric::ublas;

Simplex_Pt::Simplex_Pt( int mode ){
	this->mode_ = mode;
}

Simplex_Pt::~Simplex_Pt(){}

int Simplex_Pt::init_size(){
	/* --------
	Initializes vectors and matrices at the right dimension,
	assuming data preview has already been given 
	along with horizon and time increment
	--------- */
	if( this->dimension_ <= 0 )
		return 0;
	// /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\ 
	// Check dimension/horizon/mode consistency
	int dimension = (int) this->h_;
	switch( this->mode_ ){
		case SIMPLEX_MODE_U_KP_KD:
			if( this->dimension_ != (2+1+1)*this->h_ ) return 0;
			break;
		case SIMPLEX_MODE_U_KP:
			if( this->dimension_ != (2+1)*this->h_ ) return 0;
			break;
		case SIMPLEX_MODE_U_KPCONST:
			if( this->dimension_ != (2*this->h_+1) ) return 0;
			break;
	}
	// /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\ 

	Ux_.resize(dimension); Uy_.resize(dimension);
	Kp_.resize(dimension); Kd_.resize(dimension);
	Xc_.resize(3,dimension); Yc_.resize(3,dimension);
	X_.resize(3,dimension); Y_.resize(3,dimension); Z_.resize(3,dimension);
	P_.resize(2,dimension); Pref_.resize(2,dimension);
	FDIS_.resize(3,dimension);

	x_.resize( 3 ); dx_.resize( 3 ); ddx_.resize( 3 );
	xdes_.resize( 3 );

	return 1;
}

int Simplex_Pt::set_data( const vector<double> &data ){
	/* --------
	Sets points's position
	-------- */
	if( data.size() != this->dimension_ )
		return 0;
	for( int i=0 ; i<this->dimension_; i++ )
		this->data_(i) = data(i);
	
	return 1;
}

int Simplex_Pt::set_current_kinematics( const vector<double> &x, const vector<double> &dx, const vector<double> &ddx ){
	/* --------
	Sets effector's current position
	-------- */
	if( x.size() != 3 || dx.size() != 3 || ddx.size() != 3 )
		return 0;
	for( int i=0; i<3; i++ ){
		this->x_(i) = x(i);
		this->dx_(i) = dx(i);
		this->ddx_(i) = ddx(i);
	}
	
	return 1;
}

int Simplex_Pt::set_desired_kinematics( const vector<double> &xdes, const matrix<double> &Pref ){
	/* --------
	Sets effector's and ZMP desired position
	-------- */
	if( xdes.size() != 3 || Pref.size1() != this->Pref_.size1() || Pref.size2() != this->Pref_.size2() )
		return 0;
	for( int i=0; i<3; i++ )
		this->xdes_(i) = xdes(i);
	for( unsigned int i=0; i<Pref.size1(); i++ ){
		for( unsigned int j=0; j<Pref.size2(); j++ )
			this->Pref_(i,j) = Pref(i,j);
	}

	return 1;
}

int Simplex_Pt::set_constants( double M, double zc, double gravity, double dt, double h ){
	/* --------
	Sets constants
	-------- */
	this->M_ = M;
	this->zc_ = zc;
	this->gravity_ = gravity;
	this->dt_ = dt;
	this->h_ = h;

	this->A_.resize( 3,3,false );
	this->B_.resize( 3,false );
	A_(0,0) = 1.;	A_(0,1) = dt;	A_(0,2) = dt*dt/2.;
	A_(1,0) = 0.;	A_(1,1) = 1.;	A_(1,2) = dt;
	A_(2,0) = 0.;	A_(2,1) = 0.;	A_(2,2) = 1.;
	B_(0) = dt*dt*dt/6.;	B_(1) = dt*dt/2.;	B_(2) = dt;

	return 1;
}

int Simplex_Pt::set_kpinit( double Kpinit ){
	/* --------
	Sets initial Kp for SIMPLEX_MODE_U_KPCONST
	-------- */
	this->Kpinit_ = Kpinit;
	
	if( this->mode_ != SIMPLEX_MODE_U_KPCONST )
		return 0;
	return 1;
}

int Simplex_Pt::set_matrices( matrix< double > J, matrix< double > dJ, matrix< double > Ji, matrix< double > dJJi, matrix< double > JtiHJi ){
	/* --------
	Fills some useful matrices
	-------- */
	if( dJ.size1() != dJJi.size1() || Ji.size2() != dJJi.size2() || Ji.size2() != JtiHJi.size2() )
		return 0;
	this->J_.resize( J.size1(),J.size2() ); this->dJ_.resize( dJ.size1(),dJ.size2() ); this->Ji_.resize( Ji.size1(),Ji.size2() );
	this->dJJi_.resize( dJJi.size1(),dJJi.size2() ); this->JtiHJi_.resize( JtiHJi.size1(),JtiHJi.size2() );
	for( unsigned int i=0; i<J.size1(); i++ ){
		for( unsigned int j=0; j<J.size2(); j++ ){
			this->J_(i,j) = J(i,j);
			this->dJ_(i,j) = dJ(i,j);
			this->Ji_(j,i) = Ji(j,i);
		}
	}
	for( unsigned int i=0; i<dJJi.size1(); i++ ){
		for( unsigned int j=0; j<dJJi.size2(); j++ ){
			this->dJJi_(i,j) = dJJi(i,j);
		}
	}
	for( unsigned int i=0; i<JtiHJi.size1(); i++ ){
		for( unsigned int j=0; j<JtiHJi.size2(); j++ ){
			this->JtiHJi_(i,j) = JtiHJi(i,j);
		}
	}

	return 1;
}

int Simplex_Pt::set_disturbance( const ublas::matrix< double > &fdis ){
	/* --------
	Fills disturbance preview
	-------- */
	if( this->FDIS_.size1() != fdis.size1() || this->FDIS_.size2() != fdis.size2() )
		return 0;
	for( unsigned int i=0; i<fdis.size1(); i++ ){
		for( unsigned int j=0; j<fdis.size2(); j++ )
			this->FDIS_(i,j) = fdis(i,j);
	}

	return 1;
}

void Simplex_Pt::reset_all(){
	/* --------
	Resets all vectors and matrices to zero (size is kept)
	-------- */
	Ux_ = zero_vector< double > ( Ux_.size() ); Uy_ = zero_vector< double > ( Uy_.size() );
	Kp_ = zero_vector< double > ( Kp_.size() ); Kd_ = zero_vector< double > ( Kd_.size() );
	xdes_ = zero_vector< double > ( xdes_.size() );
	Xc_ = zero_matrix< double > ( Xc_.size1(), Xc_.size2() ); Yc_ = zero_matrix< double > ( Yc_.size1(), Yc_.size2() );
	X_ = zero_matrix< double > ( X_.size1(), X_.size2() ); Y_ = zero_matrix< double > ( Y_.size1(), Y_.size2() ); Z_ = zero_matrix< double > ( Z_.size1(), Z_.size2() );
	P_ = zero_matrix< double > ( P_.size1(), P_.size2() ); Pref_ = zero_matrix< double > ( Pref_.size1(), Pref_.size2() );
}

double Simplex_Pt::func(){
	/* --------
	Objective function to be minimized
	-------- */
	double evalf = 0.;
	this->reset_all();
	// /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\ 
	// TODO: fill vectors and matrices with this->data_ according to this->mode_
	switch( this->mode_ ){
		case SIMPLEX_MODE_U_KP_KD:
			for( int k=0; k<this->h_; k++ ){
				this->Ux_(k) = this->data_(k);
				this->Uy_(k) = this->data_((int)this->h_+k);
				this->Kp_(k) = this->data_(2*(int)this->h_+k);
				this->Kd_(k) = this->data_(3*(int)this->h_+k);
			}
			break;
		case SIMPLEX_MODE_U_KP:
			for( int k=0; k<this->h_; k++ ){
				this->Ux_(k) = this->data_(k);
				this->Uy_(k) = this->data_((int)this->h_+k);
				this->Kp_(k) = this->data_(2*(int)this->h_+k);
				this->Kd_(k) = sqrt( this->Kp_(k) );
			}
			break;
		case SIMPLEX_MODE_U_KPCONST:
			for( int k=0; k<this->h_; k++ ){
				this->Ux_(k) = this->data_(k);
				this->Uy_(k) = this->data_((int)this->h_+k);
				this->Kp_(k) = this->data_(2*(int)this->h_);
				this->Kd_(k) = sqrt( this->Kp_(k) );
			}
			break;
			
	}
	// /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\ 
	identity_matrix< double > id3							( 3 );
	identity_matrix< double > id2							( 2 );
	vector< double, bounded_array< double, 2 > > vtmp2		( 2 );
	vector< double, bounded_array< double, 3 > > vtmp3		( 3 );
	double dtmp;
	vector< double, bounded_array< double, 3 > > ddx_cmd	( 3 );
	vector< double, bounded_array< double, 3 > > Fk			( 3 );
	matrix< double, row_major, bounded_array<double, BOUNDED_ARRAY_MAX_SIZE> > Ke, Ce;

	for( int k=0; k<this->h_; k++ ){
		double uxk, uyk, Kpk, Kdk;
		uxk = Ux_( k ); uyk = Uy_( k );
		Kpk = Kp_( k ); Kdk = Kd_( k );
		Ke = Kpk * this->JtiHJi_;
		Ce = prod( this->JtiHJi_, Kdk*id3 - this->dJJi_ );
		matrix< double, row_major, bounded_array<double, BOUNDED_ARRAY_MAX_SIZE> > 
			Me ( this->JtiHJi_ );
		/* -------- Computing expected input acceleration -------- */
		for( int i=0; i<3; i++ ){
			if( i==0 )		vtmp3(i) = X_(0,k) - xdes_(i);
			else if( i==1 )	vtmp3(i) = Y_(0,k) - xdes_(i);
			else if( i==2 )	vtmp3(i) = Z_(0,k) - xdes_(i);
		}
		ddx_cmd = -Kpk * vtmp3;
		vtmp3(0) = X_(1,k); vtmp3(1) = Y_(1,k); vtmp3(2) = Z_(1,k);
		ddx_cmd += -Kdk * vtmp3;
		/* -------- Computing force acting on CoM -------- */
		vtmp3(0) = X_(1,k) - Xc_(1,k);
		vtmp3(1) = Y_(1,k) - Yc_(1,k);
		vtmp3(2) = Z_(1,k) - 0.; // No vertical CoM velocity assumed
		Fk = prod( Ce, vtmp3 );
		for( int i=0; i<3; i++ ){
			if( i==0 )		vtmp3(i) = X_(0,k) - xdes_(i);
			else if( i==1 )	vtmp3(i) = Y_(0,k) - xdes_(i);
			else if( i==2 )	vtmp3(i) = Z_(0,k) - xdes_(i);
		}
		Fk += prod( Ke, vtmp3 );
		/* -------- Computing ZMP position -------- */
		dtmp = M_ * zc_ / ( M_ * gravity_ - Fk(2) );
		vtmp3(0) = 1.; vtmp3(1) = 0.; vtmp3(2) = dtmp;
		for( int i=0; i<2; i++ ){
			if( i==0 )		P_(i,k) = inner_prod( vtmp3, column( Xc_, k ) ) + dtmp * Fk(i);
			else if( i==1 )	P_(i,k) = inner_prod( vtmp3, column( Yc_, k ) ) + dtmp * Fk(i);
		}
		/* -------- Adding ZMP tracking error to objective -------- */
		evalf += inner_prod( column( P_,k ) - column( Pref_,k ), 
			prod( this->QeonR * id2, column( P_,k ) - column( Pref_,k ) ) );
		/* -------- Adding Manipulation tracking error to objective -------- */
		vtmp3(0) = X_(0,k); vtmp3(1) = Y_(0,k); vtmp3(2) = Z_(0,k);
		evalf += inner_prod( vtmp3 - xdes_, 
			prod( this->QtonR * id3, vtmp3 - xdes_ ) );
		if( k<(this->h_-1) ){
			/* -------- Integrating CoM position -------- */
			vtmp3 = prod( A_, column( Xc_,k ) ) + uxk * B_;
			for( int i=0; i<3; i++ ) Xc_( i,k+1 ) = vtmp3(i);
			vtmp3 = prod( A_, column( Yc_,k ) ) + uyk * B_;
			for( int i=0; i<3; i++ ) Yc_( i,k+1 ) = vtmp3(i);
			/* -------- Integrating Effector position -------- */
			X_(0,k+1) = X_(0,k) + dt_*X_(1,k) + dt_*dt_/2.*X_(2,k);
			X_(1,k+1) = X_(1,k) + dt_*X_(2,k);
			X_(2,k+1) = ddx_cmd(0);
			Y_(0,k+1) = Y_(0,k) + dt_*Y_(1,k) + dt_*dt_/2.*Y_(2,k);
			Y_(1,k+1) = Y_(1,k) + dt_*Y_(2,k);
			Y_(2,k+1) = ddx_cmd(1);
			Z_(0,k+1) = Z_(0,k) + dt_*Z_(1,k) + dt_*dt_/2.*Z_(2,k);
			Z_(1,k+1) = Z_(1,k) + dt_*Z_(2,k);
			Z_(2,k+1) = ddx_cmd(2);
			/* -------- Adding Input change to objective -------- */
			vtmp2(0) = Ux_(k+1)-Ux_(k);
			vtmp2(1) = Uy_(k+1)-Uy_(k);
			evalf += inner_prod( vtmp2, vtmp2 );
		}
	}
	return evalf;
}