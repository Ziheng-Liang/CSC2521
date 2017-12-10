using namespace Eigen;

struct FemOut
{
	VectorXd stress;
	MatrixXd acceleration;
};

FemOut femAccelerations(MatrixXi TT, MatrixXd referenceTV, MatrixXd worldTV);
MatrixXd edgeMatrix(MatrixXd TT);
MatrixXd cauchyStress(MatrixXd F);