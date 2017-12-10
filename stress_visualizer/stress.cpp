#include <Eigen/Dense>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "stress.h"

using namespace Eigen;


FemOut femAccelerations(MatrixXi TT, MatrixXd referenceTV, MatrixXd worldTV) {
	MatrixXd mass = MatrixXd::Zero(referenceTV.rows(), 3);
	MatrixXd force = MatrixXd::Zero(referenceTV.rows(), 3);
    MatrixXd acceleration = MatrixXd::Zero(referenceTV.rows(), 3);
    VectorXd vmass = VectorXd::Zero(referenceTV.rows());
    VectorXd stress = VectorXd::Zero(referenceTV.rows());

	double density = 200;
	for (int i = 0; i < TT.rows(); i++) {
        //ith tet
		MatrixXd currentRTT(4,3);
		MatrixXd currentWTT(4,3);
		for (unsigned j=0; j<4;j++) {
			currentRTT.row(j) = referenceTV.row(TT(i, j));
			currentWTT.row(j) = worldTV.row(TT(i, j));
	    }
	    //use world space for volumn
	    MatrixXd refEdges = edgeMatrix(currentRTT);
	    MatrixXd worldEdges = edgeMatrix(currentWTT);
	    double volumn = abs(worldEdges.determinant() / 6);
	    double currentMass = density * volumn;

        MatrixXd defG = worldEdges * refEdges.inverse();

        MatrixXd sigma = cauchyStress(defG.transpose());
        double sigma_m = sigma.norm() / 4;
	    for (unsigned j=0; j<4;j++) {
            //jth vertex
	    	vmass(TT(i, j), 0) += currentMass / 4;
            stress(TT(i, j)) += sigma_m;
	    }
        
        Vector3d v0 = worldTV.row(TT(i, 0));
        Vector3d v1 = worldTV.row(TT(i, 1));
        Vector3d v2 = worldTV.row(TT(i, 2));
        Vector3d v3 = worldTV.row(TT(i, 3));

        Vector3d f1f;
        Vector3d f2f;
        Vector3d f3f;
        Vector3d f4f;

        Vector3d normal1 = (v1-v0).cross(v3-v1); // v2
        Vector3d normal2 = (v1-v0).cross(v2-v1); // v3
        Vector3d normal3 = (v2-v0).cross(v3-v2); // v1
        Vector3d normal4 = (v2-v1).cross(v3-v2); // v0

        if ((normal1 + v1 - v2).norm() > (-normal1 + v1 - v2).norm()) {
            normal1 = -normal1;
        }
        if ((normal2 + v1 - v3).norm() > (-normal2 + v1 - v3).norm()) {
            normal2 = -normal2;
        }
        if ((normal3 + v2 - v1).norm() > (-normal3 + v2 - v1).norm()) {
            normal3 = -normal3;
        }
        if ((normal4 + v2 - v0).norm() > (-normal4 + v2 - v0).norm()) {
            normal4 = -normal4;
        }
        f1f = sigma * normal1 / 3;
        f2f = sigma * normal2 / 3;
        f3f = sigma * normal3 / 3;
        f4f = sigma * normal4 / 3;
        force.row(TT(i, 0)) += (f1f + f2f + f3f);
        force.row(TT(i, 1)) += (f1f + f2f + f4f);
        force.row(TT(i, 2)) += (f2f + f3f + f4f);
        force.row(TT(i, 3)) += (f1f + f3f + f4f);
        // if (force.array().sum() != force.array().sum()) {
        // if (i == 2) {
        //     std::cout << "i: \n" << i << "\n";
        //     std::cout << "X.colwise().minCoeff()\n" << worldTV.colwise().minCoeff() << "\n";
        //     std::cout << "currentRTT:\n" << currentRTT << "\n";
        //     std::cout << "currentWTT:\n" << currentWTT << "\n";
        //     std::cout << "refEdges:\n" << refEdges << "\n";
        //     std::cout << "worldEdges:\n" << worldEdges << "\n";
        //     std::cout << "volumn:\n" << volumn << "\n";
        //     std::cout << "defG:\n" << defG << "\n";
        //     std::cout << "sigma:\n" << sigma << "\n";
        //     std::cout << "f1f:\n" << f1f << "\n";
        //     std::cout << "f2f:\n" << f2f << "\n";
        //     std::cout << "f3f:\n" << f3f << "\n";
        //     std::cout << "f4f:\n" << f4f << "\n";
        //     getchar();
        // }
	}
    MatrixXd gravity(referenceTV.rows(), 3);
    MatrixXd a = MatrixXd::Constant(referenceTV.rows(), 1, -0.98);
    MatrixXd b = MatrixXd::Constant(referenceTV.rows(), 2, 0);
    gravity << b, a;
    mass << vmass, vmass, vmass;
    acceleration = (force.array() / mass.array()).matrix() + gravity;
    FemOut o = {stress, acceleration};
    return o;
}

MatrixXd edgeMatrix(MatrixXd TT) {
	MatrixXd edgeM = MatrixXd::Zero(3, 3);
 	
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
    x1 = TT(0, 0);
    x2 = TT(1, 0);
    x3 = TT(2, 0);
    x4 = TT(3, 0);
    y1 = TT(0, 1);
    y2 = TT(1, 1);
    y3 = TT(2, 1);
    y4 = TT(3, 1);
    z1 = TT(0, 2);
    z2 = TT(1, 2);
    z3 = TT(2, 2);
    z4 = TT(3, 2);
    
    edgeM << x1 - x4, x2 - x4, x3 - x4,
            y1 - y4, y2 - y4, y3 - y4,
            z1 - z4, z2 - z4, z3 - z4;
    return edgeM;
}

MatrixXd cauchyStress(MatrixXd F) {
    double mu = 10000;
    double lambda = 2000;
    double J = F.determinant();
    MatrixXd invFT = (F.transpose()).inverse();
    MatrixXd P = mu * (F - invFT) + lambda * std::log(J) * invFT;
    MatrixXd stress = (P*F.transpose()) / J;
    return stress;
}