#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>
#include "stress.h"
#include "odeint_method.h"

using namespace Eigen;
using namespace boost::numeric::odeint;

void femOde(const state_type &x , state_type &dxdt , const double t) {
	worldTV =  x.block(0,0,worldTV.rows(), 3);
	velocity = x.block(worldTV.rows(), 0, worldTV.rows(), 3);
	velocity.col(2) = (worldTV.col(2).array() <= 0).select(VectorXd::Constant(worldTV.rows(), 0), velocity.col(2));
	worldTV.col(2) = (worldTV.col(2).array() < 0).select(VectorXd::Constant(worldTV.rows(), 0),worldTV.col(2));
    FemOut out = femAccelerations(TT, referenceTV, worldTV);
    dxdt = MatrixXd::Zero(x.rows(), x.cols());
    dxdt.block(0,0,worldTV.rows(),3) = velocity;
    MatrixXd a = out.acceleration;
    a = (a.array() != a.array()).select(MatrixXd::Zero(a.rows(), 3),a);
    // std::cout << a << std::endl;
    clr.push_back(out.stress);
    dxdt.block(worldTV.rows(), 0, worldTV.rows(), 3) = a;
    }

// struct femOdeBack
// {
//     void operator()( const vector_type &x , matrix_type &J , const double & ) const
//     {
        
//     }
// };