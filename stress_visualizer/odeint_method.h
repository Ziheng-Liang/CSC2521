#include "header.h"
// #include <boost/numeric/odeint/external/vexcl/vexcl_norm_inf.hpp>
// #include <boost/numeric/odeint/external>

using namespace Eigen;
using namespace boost::numeric::odeint;

typedef MatrixXd state_type;
typedef runge_kutta_dopri5<state_type,double,state_type,double,vector_space_algebra> stepper;
// typedef runge_kutta_cash_karp54<state_type,double,state_type,double> error_stepper_type;
// typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
// controlled_stepper_type controlled_stepper;

void femOde(const state_type &x , state_type &dxdt , const double /* t */);