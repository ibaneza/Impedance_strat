#ifndef _SIMPLEX_PT_HEADER_
#define _SIMPLEX_PT_HEADER_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "constants.h"

using namespace boost::numeric;

class Simplex_Pt{
public:
	Simplex_Pt( int mode=SIMPLEX_MODE_U_KPCONST );
	~Simplex_Pt();

private:
	/* -------- Point-related stuff ---------- */
	ublas::vector<double> data_;	// Contains point position in parameters space
	double distance_;				// Objective function evaluation of current point position
	int rank_;						// Rank amongst simplex points
	int dimension_;					// Number of parameters
	int mode_;						// Parameter space
	
	/* -------- Objective function-related stuff -------- */
	double QtonR, QeonR;				// It's all in the title

	ublas::vector< double > x_, dx_, ddx_;	// Manipulation task	current position, velocity, acceleration
	ublas::vector< double > xdes_;			// Manipulation tas		desired position

	ublas::vector< double > Ux_, Uy_;	// CoM					input jerk
	ublas::vector< double > Kp_, Kd_;	// Manipulation task	PD gains

	ublas::matrix< double > Xc_, Yc_;	// CoM		position
	ublas::matrix< double > X_, Y_, Z_;	// Effector	position
	ublas::matrix< double > P_, Pref_;	// ZMP		position
	
	double M_, zc_, gravity_, dt_, h_;	//	Total mass, CoM altitude (assumed constant), 
										// ... Gravity Constant, Time increment, Preview horizon

	ublas::matrix< double > FDIS_;			// Total Disturbance	force Horizon
	ublas::matrix< double > J_, dJ_, Ji_;	// Some Jacobian		derivates
	ublas::matrix< double > dJJi_, JtiHJi_;	// Some Jacobian		composites

	ublas::matrix< double > A_;			// Time integration matrix
	ublas::vector< double > B_;			// Time integration vector

	/* -------- Mode-dependent stuff -------- */
	double Kpinit_;
	
public:
	/* -------- Setters -------- */
	void set_coeff( double Qt, double Qe )
				{ this->QtonR=Qt; this->QeonR=Qe;};
	int set_data( const ublas::vector<double> &data );
	int set_disturbance( const ublas::matrix<double> &fdis );
	int set_current_kinematics( const ublas::vector<double> &x, const ublas::vector<double> &dx, const ublas::vector<double> &ddx ); 
	int set_desired_kinematics( const ublas::vector<double> &xdes, const ublas::matrix<double> &Pref ); 
	int set_constants( double M, double zc, double gravity, double dt, double h );
	int set_kpinit( double Kpinit );
	int set_matrices( ublas::matrix< double > J, ublas::matrix< double > dJ, ublas::matrix< double > Ji, ublas::matrix< double > dJJi, ublas::matrix< double > JtiHJi );
	/* -------- Init -------- */
	int init_size();
	void reset_all();
	double func();
};


#endif // _SIMPLEX_PT_HEADER_