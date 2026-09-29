//-------------------------------------------------
// Header file for Fwd_kinematics.cpp
//-------------------------------------------------

#ifndef EXPERIMENTALCODES_FWD_KINEMATICS_H
#define EXPERIMENTALCODES_FWD_KINEMATICS_H

#include <Eigen/Dense>
#include <iostream>
#include <tuple>

using namespace Eigen;


//----------------------------------------------
// Functions used in the main codes
//----------------------------------------------
// forward: executing forward kinematics
//----------------------------------------------
namespace Fwd_kinematics {
    std::tuple<VectorXd, MatrixXd> forward(const VectorXd& q);
}



#endif //EXPERIMENTALCODES_FWD_KINEMATICS_H
