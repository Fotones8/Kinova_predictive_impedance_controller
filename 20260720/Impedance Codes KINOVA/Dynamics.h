//-------------------------------------------------
// Header file for Filter.cpp
//-------------------------------------------------

#ifndef EXPERIMENTALCODES_DYNAMICS_H
#define EXPERIMENTALCODES_DYNAMICS_H

#include <Eigen/Dense>
using namespace Eigen;


//----------------------------------------------
// Functions used in the main codes
//----------------------------------------------
// mass_m: computation of mass matrix
// gravity_m: computation of gravity matrix
// coriolis_m: computation of coriolis matrix
//----------------------------------------------
namespace Dynamics {
    MatrixXd mass_m(const VectorXd &q);
    MatrixXd gravity_m(const VectorXd &q);
    MatrixXd coriolis_m(const VectorXd &q, const VectorXd &dq);
}


#endif //EXPERIMENTALCODES_DYNAMICS_H
