#include "Simplex_Pt.h"

using namespace boost::numeric::ublas;

Simplex_Pt::Simplex_Pt( int mode ){
	this->mode_ = mode;
}

Simplex_Pt::~Simplex_Pt(){}

int Simplex_Pt::init( Constants_holder ch, int dim ){
	/* --------
	Initializes ALL except data and coefficents!
	--------- */
	if( !ch.get_state() ){
		std::cout<<1;return 0;}
	if( !this->set_constants( ch.M_, ch.zc_, ch.g_, ch.dt_, ch.h_ ) ){
		std::cout<<2;return 0;}
	if( !init_size( dim ) ){
		std::cout<<3;return 0;}
	if( !this->set_current_kinematics( ch.x_, ch.dx_, ch.ddx_,
		ch.xc_, ch.dxc_, ch.ddxc_) ){
		std::cout<<4;return 0;}
	if( !this->set_desired_kinematics( ch.xdes_, ch.Pref_ ) ){
		std::cout<<5;return 0;}
	if( !this->set_matrices( ch.Mei_, ch.J_, ch.dJ_, ch.Ji_, ch.dJJi_, ch.JtiHJi_ ) ){
		std::cout<<6;return 0;}
	if( !this->set_disturbance( ch.FDIS_ ) ){
		std::cout<<7;return 0;}

	return 1;
}

int Simplex_Pt::init_size( int dim ){
	/* --------
	Initializes vectors and matrices at the right dimension,
	assuming data preview has already been given 
	along with horizon and time increment
	--------- */
	if( dim <= 0 )
		return 0;
	this->dimension_ = dim;
	// /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\ 
	// Check dimension/horizon/mode consistency
	int dimension = (int) this->h_;
	switch( this->mode_ ){
		case SIMPLEX_MODE_U_KP_KD:
			if( this->dimension_ != (2+1+1)*this->h_ )	return 0;
			break;
		case SIMPLEX_MODE_U_KP:
			if( this->dimension_ != (2+1)*this->h_ )	return 0;
			break;
		case SIMPLEX_MODE_U_KPCONST:
			if( this->dimension_ != (2*this->h_+1) )	return 0;
			break;
		case SIMPLEX_MODE_U:
			if( this->dimension_ != (2*this->h_) )		return 0;
			break;
		case SIMPLEX_MODE_KP:
			if( this->dimension_ != this->h_ )			return 0;
			break;
	}
	// /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\ 

	if( Ux_.size() != dimension ){
		Ux_.resize(dimension); Uy_.resize(dimension);
		Kp_.resize(dimension); Kd_.resize(dimension);
		Xc_.resize(3,dimension); Yc_.resize(3,dimension);
		X_.resize(3,dimension); Y_.resize(3,dimension); Z_.resize(3,dimension);
		P_.resize(2,dimension); Pref_.resize(2,dimension);
		FDIS_.resize(3,dimension);

		x_.resize( 3 ); dx_.resize( 3 ); ddx_.resize( 3 );
		xc_.resize( 2 ); dxc_.resize( 2 ); ddxc_.resize( 2 );
		xdes_.resize( 3 );
	}

	return 1;
}

int Simplex_Pt::set_data( const vector<double> &data ){
	/* --------
	Sets points's position
	-------- */
	if( data.size() != this->dimension_ )
		return 0;
	if( this->data_.size() != this->dimension_ )
		this->data_.resize( this->dimension_ );
	this->data_ = vector< double > (data);
	/*for( int i=0 ; i<this->dimension_; i++ )
		this->data_(i) = data(i);*/
	
	return 1;
}

int Simplex_Pt::set_current_kinematics( const vector<double> &x, const vector<double> &dx, const vector<double> &ddx,
									   const vector<double> &xc, const vector<double> &dxc, const vector<double> &ddxc){
	/* --------
	Sets effector's current position
	-------- */
	if( x.size() != 3 || dx.size() != 3 || ddx.size() != 3 )
		return 0;
	for( int i=0; i<3; i++ ){
		this->x_(i) = x(i);
		this->dx_(i) = dx(i);
		this->ddx_(i) = ddx(i);
		if( i<2 ){
			this->xc_ (i) = xc(i);
			this->dxc_ (i) = dxc(i);
			this->ddxc_ (i) = ddxc(i);
		}
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

int Simplex_Pt::set_matrices( matrix<double >Mei, matrix< double > J, matrix< double > dJ, matrix< double > Ji, matrix< double > dJJi, matrix< double > JtiHJi ){
	/* --------
	Fills some useful matrices
	-------- */
	if( dJ.size1() != dJJi.size1() || Ji.size2() != dJJi.size2() || Ji.size2() != JtiHJi.size2() )
		return 0;
	this->Mei_.resize( Mei.size1(), Mei.size2() );
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
	for( unsigned int i=0; i<Mei.size1(); i++ ){
		for( unsigned int j=0; j<Mei.size2(); j++ ){
			this->Mei_(i,j) = Mei(i,j);
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
	Xc_ = zero_matrix< double > ( Xc_.size1(), Xc_.size2() ); Yc_ = zero_matrix< double > ( Yc_.size1(), Yc_.size2() );
	X_ = zero_matrix< double > ( X_.size1(), X_.size2() ); Y_ = zero_matrix< double > ( Y_.size1(), Y_.size2() ); Z_ = zero_matrix< double > ( Z_.size1(), Z_.size2() );
	P_ = zero_matrix< double > ( P_.size1(), P_.size2() ); 
}

double Simplex_Pt::func(bool bverb){
	/* --------
	Objective function to be minimized
	-------- */
	double evalf = 0.;
	//this->reset_all();
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
				this->Kd_(k) = sqrt( fabs(this->Kp_(k)) );
			}
			break;
		case SIMPLEX_MODE_U_KPCONST:
			for( int k=0; k<this->h_; k++ ){
				this->Ux_(k) = this->data_(k);
				this->Uy_(k) = this->data_((int)this->h_+k);
				this->Kp_(k) = this->data_(2*(int)this->h_);
				this->Kd_(k) = sqrt( fabs(this->Kp_(k)) );
			}
			break;
		case SIMPLEX_MODE_U:
			for( int k=0; k<this->h_; k++ ){
				this->Ux_(k) = this->data_(k);
				this->Uy_(k) = this->data_((int)this->h_+k);
				this->Kp_(k) = this->Kpinit_;
				this->Kd_(k) = sqrt( fabs(this->Kp_(k)) );
			}
			break;
		case SIMPLEX_MODE_KP:
			for( int k=0; k<this->h_; k++ ){
				this->Kp_(k) = this->data_(k);
				this->Kd_(k) = sqrt( fabs(this->Kp_(k)) );
			}
			//std::cout<<Kp_<<std::endl;
			break;
			
	}
	// /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\ 
	if( this->mode_ != SIMPLEX_MODE_KP ){
		identity_matrix< double > id3							( 3 );
		identity_matrix< double > id2							( 2 );
		vector< double, bounded_array< double, 2 > > vtmp2		( 2 );
		vector< double, bounded_array< double, 3 > > vtmp3		( 3 );
		double dtmp;
		vector< double, bounded_array< double, 3 > > ddx_cmd, dx_cmd	( 3 );
		vector< double, bounded_array< double, 3 > > Fk			( 3 );
		matrix< double, row_major, bounded_array<double, BOUNDED_ARRAY_MAX_SIZE> > Ke, Ce;
		/* -------- Initialization of Effector position -------- */
		X_(0,0) = this->x_(0); X_(1,0) = this->dx_(0); X_(2,0) = this->ddx_(0);
		Y_(0,0) = this->x_(1); Y_(1,0) = this->dx_(1); Y_(2,0) = this->ddx_(1);
		Z_(0,0) = this->x_(2); Z_(1,0) = this->dx_(2); Z_(2,0) = this->ddx_(2);
		/* -------- Initialization of CoM position -------- */
		Xc_(0,0) = this->xc_(0); Xc_(1,0) = this->dxc_(0); Xc_(2,0) = this->ddxc_(0);
		Yc_(0,0) = this->xc_(1); Yc_(1,0) = this->dxc_(1); Yc_(2,0) = this->ddxc_(1);
		dx_cmd = zero_vector< double > (3);
		for( int k=0; k<this->h_; k++ ){
			double uxk, uyk, Kpk, Kdk;
			uxk = Ux_( k ); uyk = Uy_( k );
			Kpk = Kp_( k ); Kdk = Kd_( k );
			evalf += 1.e-2 * pow( (Kpk - this->Kpinit_) / this->Kpinit_, 2);
			evalf += 1.e6 * pow( fabs(Kpk-10.) - (Kpk-10.), 2.);
			Ke = Kpk * this->JtiHJi_;
			Ce = prod( this->JtiHJi_, Kdk*id3 - this->dJJi_ );
			//std::cout<<Mei_<<std::endl;
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
			vtmp3(0) = X_(2,k); vtmp3(1) = Y_(2,k); vtmp3(2) = Z_(2,k);
			ddx_cmd += -1. * vtmp3;
			if(bverb && false){
				std::cout<<"Mei = "<<Mei_<<std::endl;
				//std::cout<<"Xk = "<<xdes_<<std::endl;
			}
			dx_cmd += dt_*ddx_cmd;
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
			// /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\ 
			//Fk += column( this->FDIS_, k ); // Sure about that?
			// /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\ 
			/* -------- Computing ZMP position -------- */
			dtmp = M_ * zc_ / ( M_ * gravity_ - Fk(2) );
			vtmp3(0) = 1.; vtmp3(1) = 0.; vtmp3(2) = -1. * dtmp;
			for( int i=0; i<2; i++ ){
				if( i==0 )		P_(i,k) = inner_prod( vtmp3, column( Xc_, k ) ) + dtmp * Fk(i);
				else if( i==1 )	P_(i,k) = inner_prod( vtmp3, column( Yc_, k ) ) + dtmp * Fk(i);
			}
			/* -------- Adding ZMP tracking error to objective -------- */
			evalf += 1. / pow(norm_2( column( Pref_,k ) ),2 ) 
				* inner_prod( column( P_,k ) - column( Pref_,k ), 
				prod( this->QeonR * id2, column( P_,k ) - column( Pref_,k ) ) );
			/* -------- Adding Manipulation tracking error to objective -------- */
			vtmp3(0) = X_(0,k); vtmp3(1) = Y_(0,k); vtmp3(2) = Z_(0,k);
			evalf += 1. / pow( norm_2( xdes_ ) ,2 )
				* inner_prod( vtmp3 - xdes_, 
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
				//X_(2,k+1) = ddx_cmd(0);
				//
				Y_(0,k+1) = Y_(0,k) + dt_*Y_(1,k) + dt_*dt_/2.*Y_(2,k);
				Y_(1,k+1) = Y_(1,k) + dt_*Y_(2,k);
				//Y_(2,k+1) = ddx_cmd(1);
				//
				Z_(0,k+1) = Z_(0,k) + dt_*Z_(1,k) + dt_*dt_/2.*Z_(2,k);
				Z_(1,k+1) = Z_(1,k) + dt_*Z_(2,k);
				//Z_(2,k+1) = ddx_cmd(2);
				//
				X_(2,k+1) = 0;	Y_(2,k+1) = 0;	Z_(2,k+1) = 0;
				X_(2,k+1) = -Kpk*X_(0,k) + (-dt_*Kpk-Kdk)*X_(1,k) 
					+ (-dt_*Kdk-dt_*dt_/2.*Kpk)*X_(1,k) + Kpk*xdes_(0);
				Y_(2,k+1) = -Kpk*Y_(0,k) + (-dt_*Kpk-Kdk)*Y_(1,k) 
					+ (-dt_*Kdk-dt_*dt_/2.*Kpk)*Y_(1,k) + Kpk*xdes_(1);
				Z_(2,k+1) = -Kpk*Z_(0,k) + (-dt_*Kpk-Kdk)*Z_(1,k) 
					+ (-dt_*Kdk-dt_*dt_/2.*Kpk)*Z_(1,k) + Kpk*xdes_(2);
				ublas::matrix < double > tmp;
				tmp = prod( dJ_, Ji_ );
				vtmp3(0) = X_(1,k+1); vtmp3(1) = Y_(1,k+1); vtmp3(2) = Z_(1,k+1);
				vtmp3 = prod( Mei_, column(FDIS_,k+1)) 
					+ 1.*prod(tmp,dx_cmd)
					- 1.*prod(tmp, vtmp3);
				X_(2,k+1) += vtmp3(0); Y_(2,k+1) += vtmp3(1); Z_(2,k+1) += vtmp3(2);
				/* -------- Adding Input change to objective -------- */
				if( Ux_(k+1)!=0. ) vtmp2(0) = ( Ux_(k+1)-Ux_(k) ) / Ux_(k+1);
				else vtmp2(0) = 0.;
				if( Uy_(k+1)!=0. ) vtmp2(1) = ( Uy_(k+1)-Uy_(k) ) / Uy_(k+1);
				else vtmp2(1) = 0.;
				evalf += inner_prod( vtmp2, vtmp2 );
				if( bverb ){}
				//std::cout<<evalf<<std::endl;
			}
		}
		if( evalf != evalf )
			evalf = 1.e12;
	}
	else{
		//std::cout<<"Computing Matrices"<<std::endl;
		this->nPref_.resize(3,this->h_);
		matrix< double> F = this->build_Effort();
		vector< matrix< double > > PXU = this->build_Pxu( F );
		//matrix< double > nPref = this->build_newPref( F );
		vector< double, bounded_array< double, 3 > > vtmp3( 3 );
		matrix< double > Px = PXU(0);
		matrix< double > Pu = PXU(1);
		vector< double > u;
		identity_matrix< double > idh (this->h_);
		matrix< double > M, I;
		M = prod( Pu, trans(Pu) ) + 1./this->QeonR * idh ;
		I.resize( M.size2(), M.size1() );
		//std::cout<<"Inverting Matrix"<<std::endl;
		InvertMatrix( M, I );
		//std::cout<<"Processing"<<std::endl;
		for( int i=0; i<2; i++ ){
			vector< double > Pref = row( this->nPref_, i );
			vtmp3(0) = this->xc_(i);	vtmp3(1) = this->dxc_(i);	vtmp3(2) = this->ddxc_(i);
			u = prod( prod( -I, trans(Pu) ), prod( Px, vtmp3 ) - Pref );
			if( i==0 ){
				this->Ux_ = u;
			}
			else if( i==1 )
				this->Uy_ = u;
		}
		this->integrate_CoM();
		evalf = this->compute_error(F);
	}
	return evalf;
}

int Constants_holder::fill_matrix(matrix<double> src, matrix<double> dest ){
	int sz1 = src.size1(), sz2 = src.size2();
	dest.resize( sz1,sz2 );
	for( int i=0; i<sz1; i++ ){
		for( int j=0; j<sz1; j++ )
			dest(i,j) = src(i,j);
	}
	state++;
	
	return 1;
}

int Constants_holder::fill_vector(vector<double> src, vector<double> dest ){
	int sz = src.size();
	dest.resize( sz );
	for( int i=0; i<sz; i++ )
		dest(i) = src(i);
	state++;

	return 1;
}

matrix< double> Simplex_Pt::build_Effort(){
	/* --------
	Builds equivalent effort on CoM
	-------- */
	matrix< double > effort (3,this->h_);

	identity_matrix< double > id3							( 3 );
	identity_matrix< double > id2							( 2 );
	vector< double, bounded_array< double, 2 > > vtmp2		( 2 );
	vector< double, bounded_array< double, 3 > > vtmp3		( 3 );
	double dtmp;
	vector< double, bounded_array< double, 3 > > ddx_cmd, dx_cmd	( 3 );
	vector< double, bounded_array< double, 3 > > Fk			( 3 );
	matrix< double, row_major, bounded_array<double, BOUNDED_ARRAY_MAX_SIZE> > Ke, Ce;
	/* -------- Initialization of Effector position -------- */
	X_(0,0) = this->x_(0); X_(1,0) = this->dx_(0); X_(2,0) = this->ddx_(0);
	Y_(0,0) = this->x_(1); Y_(1,0) = this->dx_(1); Y_(2,0) = this->ddx_(1);
	Z_(0,0) = this->x_(2); Z_(1,0) = this->dx_(2); Z_(2,0) = this->ddx_(2);
	dx_cmd = zero_vector< double > (3);
	ublas::matrix < double > tmp;
	tmp = prod( dJ_, Ji_ );
	for( int k=0; k<this->h_; k++ ){
		double Kpk, Kdk;
		Kpk = Kp_( k ); Kdk = Kd_( k );
		Ke = Kpk * this->JtiHJi_;
		Ce = prod( this->JtiHJi_, Kdk*id3 - this->dJJi_ );
		//std::cout<<Mei_<<std::endl;
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
		vtmp3(0) = X_(2,k); vtmp3(1) = Y_(2,k); vtmp3(2) = Z_(2,k);
		ddx_cmd += -1. * vtmp3;
		dx_cmd += dt_*ddx_cmd;
		/* -------- Computing force acting on CoM -------- */
		vtmp3(0) = X_(1,k);// - Xc_(1,k);
		vtmp3(1) = Y_(1,k);// - Yc_(1,k);
		vtmp3(2) = Z_(1,k) - 0.; // No vertical CoM velocity assumed
		Fk = prod( Ce, vtmp3 );
		for( int i=0; i<3; i++ ){
			if( i==0 )		vtmp3(i) = X_(0,k) - xdes_(i);
			else if( i==1 )	vtmp3(i) = Y_(0,k) - xdes_(i);
			else if( i==2 )	vtmp3(i) = Z_(0,k) - xdes_(i);
		}
		Fk += prod( Ke, vtmp3 );
		for(int i=0;i<3;i++) effort(i,k) = Fk(i);
		for( int i=0; i<2; i++ ){
			this->nPref_(i,k) = this->Pref_(i,k);
			this->nPref_(i,k) += -this->zc_ * Fk(i) / (this->M_ * this->gravity_ - Fk(2));
		}
		if( k<(this->h_-1) ){
			/* -------- Integrating Effector position -------- */
			X_(0,k+1) = X_(0,k) + dt_*X_(1,k) + dt_*dt_/2.*X_(2,k);
			X_(1,k+1) = X_(1,k) + dt_*X_(2,k);
			//
			Y_(0,k+1) = Y_(0,k) + dt_*Y_(1,k) + dt_*dt_/2.*Y_(2,k);
			Y_(1,k+1) = Y_(1,k) + dt_*Y_(2,k);
			//
			Z_(0,k+1) = Z_(0,k) + dt_*Z_(1,k) + dt_*dt_/2.*Z_(2,k);
			Z_(1,k+1) = Z_(1,k) + dt_*Z_(2,k);
			//
			X_(2,k+1) = 0;	Y_(2,k+1) = 0;	Z_(2,k+1) = 0;
			X_(2,k+1) = -Kpk*X_(0,k) + (-dt_*Kpk-Kdk)*X_(1,k) 
				+ (-dt_*Kdk-dt_*dt_/2.*Kpk)*X_(1,k) + Kpk*xdes_(0);
			Y_(2,k+1) = -Kpk*Y_(0,k) + (-dt_*Kpk-Kdk)*Y_(1,k) 
				+ (-dt_*Kdk-dt_*dt_/2.*Kpk)*Y_(1,k) + Kpk*xdes_(1);
			Z_(2,k+1) = -Kpk*Z_(0,k) + (-dt_*Kpk-Kdk)*Z_(1,k) 
				+ (-dt_*Kdk-dt_*dt_/2.*Kpk)*Z_(1,k) + Kpk*xdes_(2);
			vtmp3(0) = X_(1,k+1); vtmp3(1) = Y_(1,k+1); vtmp3(2) = Z_(1,k+1);
			vtmp3 = prod( Mei_, column(FDIS_,k+1)) 
				+ 1.*prod(tmp,dx_cmd)
				- 1.*prod(tmp, vtmp3);
			X_(2,k+1) += vtmp3(0); Y_(2,k+1) += vtmp3(1); Z_(2,k+1) += vtmp3(2);
		}
	}
	return effort;
}

vector< matrix< double > >Simplex_Pt::build_Px( matrix< double> F ){
	/* --------
	Builds Px spreading matrix
	-------- */
	vector< matrix< double > > Px (1);
	for( int i=0; i<1; i++ ){ 
		Px(i).resize(this->h_, 3);
		for( int lin=0; lin<this->h_; lin++ ){
			Px(i)(lin,0) = 1;
			Px(i)(lin,1) = (lin+1) * this->dt_;
			Px(i)(lin,2) = pow( (lin+1.)*this->dt_, 2. );
			Px(i)(lin,2) += -this->M_ * this->zc_ / ( this->M_ * this->gravity_ - F(2,lin));
		}
	}

	return Px;
}

vector< matrix< double > >Simplex_Pt::build_Pu( matrix< double> F ){
	/* --------
	Builds Pu spreading matrix
	-------- */
	vector< matrix< double > > Pu (1);
	for( int i=0; i<1; i++ ){ 
		Pu(i).resize(this->h_, this->h_);
		for( int lin=0; lin<this->h_; lin++ ){
			double diag_lin = (1. + 3.*lin + 3.*pow((double)lin,2.)) * pow( this->dt_, 3. ) / 6.;
			diag_lin += -this->M_ * this->zc_ / ( this->M_ * this->gravity_ - F(2,lin));
			for( int l=lin; l<this->h_; l++ ){
				for( int c=0; c<(this->h_-lin); c++ ){
					Pu(i)(l,c) = diag_lin;
				}
			}
		}
	}

	return Pu;
}

vector< matrix< double > >Simplex_Pt::build_Pxu( matrix< double> F ){
	/* --------
	Builds Px and Pu spreading matrices
	-------- */
	vector< matrix< double > > Pxu (2);
	Pxu(1).resize(this->h_, this->h_);
	Pxu(0).resize(this->h_, 3);
	for( int lin=0; lin<this->h_; lin++ ){
		Pxu(0)(lin,0) = 1;
		Pxu(0)(lin,1) = (lin+1) * this->dt_;
		Pxu(0)(lin,2) = pow( (lin+1.)*this->dt_, 2. );
		Pxu(0)(lin,2) += -this->M_ * this->zc_ / ( this->M_ * this->gravity_ - F(2,lin));
		double diag_lin = (1. + 3.*lin + 3.*pow((double)lin,2.)) * pow( this->dt_, 3. ) / 6.;
		diag_lin += -this->M_ * this->zc_ / ( this->M_ * this->gravity_ - F(2,lin));
		for( int l=lin; l<this->h_; l++ ){
			for( int c=0; c<(this->h_-lin); c++ ){
				Pxu(1)(l,c) = diag_lin;
			}
		}
	}

	return Pxu;
}

matrix< double > Simplex_Pt::build_newPref( matrix< double> F ){
	/* --------
	Builds new ZMP reference position according to equivalent effort
	-------- */
	matrix< double > nPref (2, this->h_);
	for( int i=0; i<2; i++ ){ 
		for( int k=0; k<this->h_; k++ ){
			nPref(i,k) = this->Pref_(i,k);
			nPref(i,k) += -this->zc_ * F(i,k) / (this->M_ * this->gravity_ - F(2,k));
		}
	}

	return nPref;
}

void Simplex_Pt::integrate_CoM(){
	/* --------
	Integrates Center of Mass position
	-------- */
}

double Simplex_Pt::compute_error( matrix< double> F ){
	/* --------
	Computes error if X_ and Xc_ previously computed
	and
	Integrates Center of Mass position
	-------- */
	vector< double, bounded_array< double, 3 > > vtmp3		( 3 );
	/* -------- Initialization of CoM position -------- */
	Xc_(0,0) = this->xc_(0); Xc_(1,0) = this->dxc_(0); Xc_(2,0) = this->ddxc_(0);
	Yc_(0,0) = this->xc_(1); Yc_(1,0) = this->dxc_(1); Yc_(2,0) = this->ddxc_(1);
	double evalf = 0.;
	identity_matrix< double > id3							( 3 );
	identity_matrix< double > id2							( 2 );
	vector< double, bounded_array< double, 2 > > vtmp2		( 2 );
	double dtmp;
	vector< double > Fk (3);
	for( int k=0; k<this->h_; k++ ){
		double Kpk, Kdk, uxk, uyk;
		uxk = Ux_( k ); uyk = Uy_( k );
		Fk = column( F, k );
		Kpk = Kp_( k ); Kdk = Kd_( k );
		evalf += 1.e-2 * pow( (Kpk - this->Kpinit_) / this->Kpinit_, 2);
		evalf += 1.e6 * pow( fabs(Kpk-10.) - (Kpk-10.), 2.);
		/* -------- Adding ZMP tracking error to objective -------- */
		evalf += 1. / pow(norm_2( column( Pref_,k ) ),2 ) 
			* inner_prod( column( P_,k ) - column( Pref_,k ), 
			prod( this->QeonR * id2, column( P_,k ) - column( Pref_,k ) ) );
		/* -------- Adding Manipulation tracking error to objective -------- */
		vtmp3(0) = X_(0,k); vtmp3(1) = Y_(0,k); vtmp3(2) = Z_(0,k);
		evalf += 1. / pow( norm_2( xdes_ ) ,2 )
			* inner_prod( vtmp3 - xdes_, 
			prod( this->QtonR * id3, vtmp3 - xdes_ ) );
		dtmp = M_ * zc_ / ( M_ * gravity_ - Fk(2) );
			vtmp3(0) = 1.; vtmp3(1) = 0.; vtmp3(2) = -dtmp;
			for( int i=0; i<2; i++ ){
				if( i==0 )		P_(i,k) = inner_prod( vtmp3, column( Xc_, k ) ) + dtmp * Fk(i);
				else if( i==1 )	P_(i,k) = inner_prod( vtmp3, column( Yc_, k ) ) + dtmp * Fk(i);
			}
		if( k<(this->h_-1) ){
			/* -------- Integrating CoM position -------- */
			vtmp3 = prod( A_, column( Xc_,k ) ) + uxk * B_;
			for( int i=0; i<3; i++ ) 
				Xc_( i,k+1 ) = vtmp3(i);
			vtmp3 = prod( A_, column( Yc_,k ) ) + uyk * B_;
			for( int i=0; i<3; i++ ) 
				Yc_( i,k+1 ) = vtmp3(i);
			/* -------- Adding Input change to objective -------- */
			if( Ux_(k+1)!=0. ) vtmp2(0) = ( Ux_(k+1)-Ux_(k) ) / Ux_(k+1);
			else vtmp2(0) = 0.;
			if( Uy_(k+1)!=0. ) vtmp2(1) = ( Uy_(k+1)-Uy_(k) ) / Uy_(k+1);
			else vtmp2(1) = 0.;
			evalf += inner_prod( vtmp2, vtmp2 );
		}
	}
	if( evalf != evalf )
		evalf = 1.e12;
	return evalf;
}
