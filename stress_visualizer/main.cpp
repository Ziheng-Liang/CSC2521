#include <igl/viewer/Viewer.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/barycenter.h>
#include <igl/jet.h>
#include <igl/decimate.h>
#include <iostream>
#include "stress.h"
#include <chrono>
#include <thread>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include "odeint_method.h"
#include "header.h"

// #include "tutorial_shared_path.h"

// Input polygon
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd B;

// Tetrahedralized interior
Eigen::MatrixXd TV;
Eigen::MatrixXi TF;

Eigen::VectorXd C;

Eigen::MatrixXi TT;
Eigen::MatrixXd worldTV;
Eigen::MatrixXd referenceTV;
Eigen::MatrixXd velocity;
std::vector<Eigen::VectorXd> clr;

double max_stress;
struct push_back_state_and_time
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};

using namespace Eigen;

void visualization(igl::viewer::Viewer& viewer, Eigen::MatrixXd thisV);

void visualization(igl::viewer::Viewer& viewer, MatrixXd V) {
  // Eigen::VectorXd Z = thisV.col(2);
  viewer.data.clear();
  viewer.data.set_mesh(V,TF);
  // Eigen::MatrixXd Z = MatrixXd::Zero(V.rows(), V.cols());
  // Eigen::MatrixXd tmp = C/max_stress;
  // // igl::jet(C,0,max_stress,Z);
  // // Eigen::VectorXd X = V.col(2);
  // // igl::jet(X,true,Z);
  // Z << C, C, C;
  // std::cout << "=========" << "\n";
  // std::cout << max_stress << "\n";
  // std::cout << C << "\n";
  // std::cout << Z << "\n";
  // // std::cout << V << "\n";
  // viewer.data.set_colors(Z);
  viewer.data.set_face_based(true);
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  using namespace boost::numeric::odeint;
  // Load a surface mesh
  Eigen::MatrixXd U;
  Eigen::MatrixXi G;
  Eigen::VectorXi I;
  Eigen::MatrixXd new_V;
  // igl::readOFF("../cube.off", V, F);
  // igl::readOFF("../sphere.off", V, F);
  igl::readOBJ("../ICO.obj", V, F);
  // igl::readOBJ("../shape.obj", V, F);
  // new_V = MatrixXd::Zero(V.rows() + 1, V.cols());
  // new_V << V, 0, 0, 1.25;

  // igl::decimate(U,G,10,V,F,I);
  // igl::copyleft::tetgen::tetrahedralize(new_V,F,"pq5a10", TV,TT,TF);
  igl::copyleft::tetgen::tetrahedralize(V,F,"pq2Y", TV,TT,TF);
  cout << "TT.rows(): " << TT.rows() << endl;
  igl::barycenter(TV,TT,B);
  referenceTV = TV;
  worldTV = TV;
  worldTV.col(2) = (worldTV.col(2).array() + 2.0).matrix();
  velocity = MatrixXd::Zero(referenceTV.rows(), 3);
  // Plot the generated mesh

  igl::viewer::Viewer viewer;
  MatrixXd x = MatrixXd::Zero(worldTV.rows()*2, 3);
  vector<state_type> x_vec;
  vector<double> times;
  x.block(0,0,worldTV.rows(),3) = worldTV;
  size_t steps = integrate_adaptive(stepper(), femOde, x, 0.0, 20.0, 0.03, push_back_state_and_time(x_vec, times));

  int count = 0;
  max_stress = 0;
  for (int i = 0; i < clr.size(); i ++) {
    if (clr[i].maxCoeff() > max_stress) {
      max_stress = clr[i].maxCoeff();
    }
  }
  viewer.callback_pre_draw = [&](igl::viewer::Viewer & viewer)->bool
  {
    C = clr[count];
    visualization(viewer, x_vec[count].block(0,0,worldTV.rows(),3));
    count += 1;
    if (count >= x_vec.size() - 1) {
      count = 0;
    }
    return false;
  };
  C = clr[10];
  visualization(viewer, x_vec[10].block(0,0,worldTV.rows(),3));
  viewer.core.is_animating = true;
  viewer.core.animation_max_fps = 30.;
  float zoom = 1000.0f;
  Vector3f shift(10,10,10);
  viewer.core.get_scale_and_shift_to_fit_mesh(TV, TF, zoom, shift);
  viewer.launch();
}
