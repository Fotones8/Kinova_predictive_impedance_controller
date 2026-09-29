//----------------------------------------------------------
// Functions of the forward kinematics
//----------------------------------------------------------
// Description: functions used to calculate the forward kinematics
// Copyright: Yihan Liu 2024
//----------------------------------------------------------

#define _USE_MATH_DEFINES

#include "Fwd_kinematics.h"

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <tuple>

using namespace Eigen;
using namespace std;


//----------------------------------------------------------
// Function to compute euler angles from rotation matrix
//----------------------------------------------------------
// 4 inputs:
// R: rotation matrix
// roll: roll angle of end effector
// pitch: pitch angle of end effector
// yaw: yaw angle of end effector
//----------------------------------------------------------
void rot_2_zyx(const Matrix3d& R, double& roll, double& pitch, double& yaw) {
    // Extracting angles from the rotation matrix
    yaw = atan2(R(1, 0), R(0, 0));
//    double cos_theta = sgn(R(1,0)) * sqrt(R(2, 1) * R(2, 1) + R(2, 2) * R(2, 2));
    double sin_theta = (1 / 0.994) * (-R(2, 0));
    if (sin_theta > 1) {
        sin_theta  = 1;
    }
    if (sin_theta < -1) {
        sin_theta = -1;
    }
    pitch = asin(sin_theta);
    if ((R(1,0) < 0) && (-R(2, 0) > 0)){
        pitch = M_PI - pitch;
    }
    if ((R(1,0) < 0) && (-R(2, 0) < 0)){
        pitch = -M_PI - pitch;
    }
    roll = atan2(R(2, 1), R(2, 2));
}


//----------------------------------------------------------
// Function to execute the forward kinematics
//----------------------------------------------------------
// 1 input:
// q: joint angular position
//
// 2 output:
// p: 6-DoF positions of end effector
// T_B7: rotation matrix of end effector
//----------------------------------------------------------
tuple<VectorXd, MatrixXd> Fwd_kinematics::forward(const VectorXd& q) {
    // Define transformation matrices
    Matrix4d T_B1, T_12, T_23, T_34, T_45, T_56, T_67, T_7end;

    T_B1 << cos(q(0)), -sin(q(0)), 0, 0,
            -sin(q(0)), -cos(q(0)), 0, 0,
            0, 0, -1, 0.1564,
            0, 0, 0, 1;

    T_12 << cos(q(1)), -sin(q(1)), 0, 0,
            0, 0, -1, 0.0054,
            sin(q(1)), cos(q(1)), 0, -0.1284,
            0, 0, 0, 1;

    T_23 << cos(q(2)), -sin(q(2)), 0, 0,
            0, 0, 1, -0.2104,
            -sin(q(2)), -cos(q(2)), 0, -0.0064,
            0, 0, 0, 1;

    T_34 << cos(q(3)), -sin(q(3)), 0, 0,
            0, 0, -1, 0.0064,
            sin(q(3)), cos(q(3)), 0, -0.2104,
            0, 0, 0, 1;

    T_45 << cos(q(4)), -sin(q(4)), 0, 0,
            0, 0, 1, -0.2084,
            -sin(q(4)), -cos(q(4)), 0, -0.0064,
            0, 0, 0, 1;

    T_56 << cos(q(5)), -sin(q(5)), 0, 0,
            0, 0, -1, 0,
            sin(q(5)), cos(q(5)), 0, -0.1059,
            0, 0, 0, 1;

    T_67 << cos(q(6)), -sin(q(6)), 0, 0,
            0, 0, 1, -0.1059,
            -sin(q(6)), -cos(q(6)), 0, 0,
            0, 0, 0, 1;

    T_7end << 1, 0, 0, 0,
            0, -1, 0, 0,
            0, 0, -1, -0.067,
            0, 0, 0, 1;

    // Compute the forward kinematics
    Matrix4d T_B7 = T_B1 * T_12 * T_23 * T_34 * T_45 * T_56 * T_67 * T_7end;

    VectorXd p(6);
    double roll, pitch, yaw;

    // Compute the 6-DoF positions of end effector
    VectorXd pos_end = T_B7.block<3, 1>(0, 3);
    Matrix3d rot_end = T_B7.block<3, 3>(0, 0);
    rot_2_zyx(T_B7.block<3, 3>(0, 0), roll, pitch, yaw);
    p << pos_end(0), pos_end(1), pos_end(2), roll, pitch, yaw;

    return make_tuple(p, T_B7);
}