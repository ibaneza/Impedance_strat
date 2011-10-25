#ifndef _SIMPLEX_PT_HEADER_
#define _SIMPLEX_PT_HEADER_

/*-------- This is bold --------*/
#define _SCL_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#define _AFX_SECURE_NO_WARNINGS
#define _ATL_SECURE_NO_WARNINGS

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "constants.h"

using namespace boost::numeric;

struct Constants_holder{
	Constants_holder():des_state(16),state(0){};
	~Constants_holder(){};

	int state;
	int des_state;

	/* -------- Some stuff to fill simplex points ---------- */
	double  M_, zc_, g_, dt_, h_, kpinit_;
	int mode_;
	ublas::vector< double > x_, dx_, ddx_;		// Manipulation task	current position, velocity, acceleration
	ublas::vector< double > xc_, dxc_, ddxc_;	// Center of Mass		current position, velocity, acceleration
	ublas::vector< double > xdes_;				// Manipulation task	desired position
	ublas::matrix< double > Mei_;				// Configuation space	Mass Matrix
	ublas::matrix< double > J_, dJ_, Ji_;		// Some Jacobian		derivates
	ublas::matrix< double > dJJi_, JtiHJi_;		// Some Jacobian		composites
	ublas::matrix< double > Pref_;				// ZMP		position
	ublas::matrix< double > FDIS_;				// Total Disturbance	force Horizon

	int fill_matrix( ublas::matrix< double > src, ublas::matrix< double > dest );
	int fill_vector( ublas::vector< double > src, ublas::vector< double > dest );
	template <typename T> int fill_scalar( T src, T dest ){ dest = src; state++; };
	int get_state(){	return (state>=des_state)?1:0; };
	int assume_state(){ this->state = des_state; return get_state(); };
};

class Simplex_Pt{
public:
	Simplex_Pt( int mode=SIMPLEX_MODE_U_KPCONST );
	~Simplex_Pt();

	ublas::matrix< double> build_Effort();
	ublas::vector< ublas::matrix< double > > build_Px(ublas::matrix< double > F);
	ublas::vector< ublas::matrix< double > > build_Pu(ublas::matrix< double > F);
	ublas::matrix< double > build_newPref(ublas::matrix< double > F);
	void integrate_CoM();
	double compute_error(ublas::matrix< double > F);

public:
	/* -------- Point-related stuff ---------- */
	ublas::vector< double > data_;	// Contains point position in parameters space
	double distance_;				// Objective function evaluation of current point position
	int rank_;						// Rank amongst simplex points
	int dimension_;					// Number of parameters
	int mode_;						// Parameter space

	/* -------- User-given stuff ---------- */
	ublas::vector< double > x_, dx_, ddx_;		// Manipulation task	current position, velocity, acceleration
	ublas::vector< double > xc_, dxc_, ddxc_;	// Center of Mass		current position, velocity, acceleration
	ublas::vector< double > xdes_;				// Manipulation task	desired position
	ublas::matrix< double > Pref_, nPref_;		// ZMP					reference position
	ublas::matrix< double > FDIS_;				// Total Disturbance	force Horizon
	ublas::matrix< double > Mei_;				// Configuation space	Mass Matrix
	ublas::matrix< double > J_, dJ_, Ji_;		// Some Jacobian		derivates
	ublas::matrix< double > dJJi_, JtiHJi_;		// Some Jacobian		composites
	double M_, zc_, gravity_, dt_, h_;			//	Total mass, CoM altitude (assumed constant), 
												// ... Gravity Constant, Time increment, Preview horizon
	
	/* -------- Objective function-related stuff -------- */
	double QtonR, QeonR;				// It's all in the title

	ublas::vector< double > Ux_, Uy_;	// CoM					input jerk
	ublas::vector< double > Kp_, Kd_;	// Manipulation task	PD gains

	ublas::matrix< double > Xc_, Yc_;	// CoM		position
	ublas::matrix< double > X_, Y_, Z_;	// Effector	position
	ublas::matrix< double > P_;	// ZMP		position	

	ublas::matrix< double > A_;			// Time integration matrix
	ublas::vector< double > B_;			// Time integration vector

	/* -------- Mode-dependent stuff -------- */
	double Kpinit_;
	
public:
	/* -------- Setters -------- */
	void set_mode( int mode )
				{ this->mode_ = mode;};
	void set_coeff( double Qt, double Qe )
				{ this->QtonR=Qt; this->QeonR=Qe;};
	int set_data( const ublas::vector<double> &data );
	int set_current_kinematics( const ublas::vector<double> &x, const ublas::vector<double> &dx, const ublas::vector<double> &ddx,
		const ublas::vector<double> &xc, const ublas::vector<double> &dxc, const ublas::vector<double> &ddxc); 
	int set_desired_kinematics( const ublas::vector<double> &xdes, const ublas::matrix<double> &Pref ); 
	int set_constants( double M, double zc, double gravity, double dt, double h );
	int set_kpinit( double Kpinit );
	int set_matrices( ublas::matrix< double> Mei, ublas::matrix< double > J, ublas::matrix< double > dJ, ublas::matrix< double > Ji, ublas::matrix< double > dJJi, ublas::matrix< double > JtiHJi );
	int set_disturbance( const ublas::matrix< double > &fdis );
	/* -------- Init -------- */
	int init( Constants_holder ch, int dim );
	int init_size( int dim );
	void reset_all();
	double func(bool bverb = false);
	/* -------- Getters -------- */
	ublas::vector< double > get_data(){ return this->data_; };
};

#endif // _SIMPLEX_PT_HEADER_