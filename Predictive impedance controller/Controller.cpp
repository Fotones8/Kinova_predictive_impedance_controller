//----------------------------------------------------------
// Functions of the impedance controller
//----------------------------------------------------------
// Description: functions used in the impedance controller
// Copyright: Yihan Liu 2024
//----------------------------------------------------------

#define _USE_MATH_DEFINES

#include "Controller.h"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <tuple>

#include <Jacobian.h>
#include <Fwd_kinematics.h>
#include <Dynamics.h>
#include <Filter.h>

#include <windows.h>

#include "training.hpp"
#include "prediction.hpp"
#include "ukf_predictor.hpp"

using namespace Eigen;
using namespace std;
using namespace Dynamics;
using namespace Jacobian;
using namespace Fwd_kinematics;
using namespace Filter;

// Buffers to store the last variables
Vector<double, 6> last_dp = Vector<double, 6>::Zero();
Matrix<double, 6, 7> last_jaco = Matrix<double, 6, 7>::Zero();
Matrix3d last_rot = Matrix<double, 3, 3>::Zero();
Vector3d last_pos = Vector<double, 3>::Zero();
Vector<double, 6> last_vel = Vector<double, 6>::Zero();
Vector<double, 6> acc_factor_buffer = Vector<double, 6>::Zero();

VectorXd u_filtered = VectorXd::Zero(7);
// constants for butterworth filter at Fs=400Hz, Fc=12Hz

double a1 = -1.7247;
double a2 = 0.7660;
double b0 = 0.0078;
double b1 = 0.0156;
double b2 = 0.0078;

double recording_time = 0;
int file_num =1;
int file_change_flag = 0;

// constants for butterworth filter at Fs=400Hz, Fc=8Hz produces oscillations
/*
double a1 = -1.8227;
double a2 = 0.8372;
double b0 = 0.0036;
double b1 = 0.0072;
double b2 = 0.0036;
*/

VectorXd f_ext_1 = VectorXd::Zero(7);
VectorXd f_ext_2 = VectorXd::Zero(7);
VectorXd f_filt_1 = VectorXd::Zero(7);
VectorXd f_filt_2 = VectorXd::Zero(7);

VectorXd u_1 = VectorXd::Zero(7);
VectorXd u_2 = VectorXd::Zero(7);
VectorXd u_filtered_1 = VectorXd::Zero(7);
VectorXd u_filtered_2 = VectorXd::Zero(7);

bool braking = false;

double baseline_Fz = 0;


#include <fstream>
std::ofstream log_file;

//Construct with default values
promp_rt::OnlinePredictor predictor = promp_rt::OnlinePredictor("C:/Imperial/Final project/predictive_controller/model", 10);

std::vector<double> initial_positions = {0.745, 0.194, 0.428};
MultiDofUkfPredictor kalman_predictor = MultiDofUkfPredictor(3,initial_positions);
//----------------------------------------------------------------
// Function to transform the euler angle to the rotation matrix
//----------------------------------------------------------------
// 3 inputs
// roll: X-axis rotation angle
// pitch: Y-axis rotation angle
// yaw: Z-axis rotation angle
//
// 1 output
// R: rotation matrix
//----------------------------------------------------------------
Matrix3d euler_to_rotation_matrix(double roll, double pitch, double yaw) {
    // Rotation matrix for roll (X-axis rotation)
    Matrix3d R_x;
    R_x << 1, 0, 0,
            0, cos(roll), -sin(roll),
            0, sin(roll), cos(roll);

    // Rotation matrix for pitch (Y-axis rotation)
    Matrix3d R_y;
    R_y << cos(pitch), 0, sin(pitch),
            0, 1, 0,
            -sin(pitch), 0, cos(pitch);

    // Rotation matrix for yaw (Z-axis rotation)
    Matrix3d R_z;
    R_z << cos(yaw), -sin(yaw), 0,
            sin(yaw), cos(yaw), 0,
            0, 0, 1;

    // Combined rotation matrix (with ZYX rotation)
    Matrix3d R = R_z * R_y * R_x;

    return R;
}


//--------------------------------------------------------------
// Function to convert the Euler angle velocity in the end frame
// into the angular velocity in the base frame
//--------------------------------------------------------------
// 2 inputs:
// euler_ang: 3 Euler angles in the end frame
// euler_vel: 3 Euler velocities in the end frame
//
// 1 output:
// omega_base: 3 angular velocities in the base frame
//--------------------------------------------------------------
Vector3d euler_vel_to_base_angular_vel(
        const Vector3d& euler_ang,  // [theta_x, theta_y, theta_z]
        const Vector3d& euler_vel  // [dtheta_x, dtheta_y, dtheta_z]
) {
    // Extract Euler angles
    double theta_x = euler_ang[0];
    double theta_y = euler_ang[1];
    double theta_z = euler_ang[2];

    // Compute the transformation matrix T_Euler^ZYX
    Matrix3d T_euler_ZYX;
    T_euler_ZYX <<
                1, 0, -sin(theta_y),
            0, cos(theta_x), cos(theta_y) * sin(theta_x),
            0, -sin(theta_x), cos(theta_y) * cos(theta_x);

    // Compute angular velocity in the end frame
    Vector3d omega_end = T_euler_ZYX * euler_vel;

    Matrix3d R = euler_to_rotation_matrix(theta_x, theta_y, theta_z);

    // Transform angular velocity from body frame to world frame
    Vector3d omega_base = R * omega_end;

    return omega_base;
}


//------------------------------------------------------------------------
// Function to compute the time derivative of the rotation matrix R
//------------------------------------------------------------------------
// 2 inputs:
// euler_ang: 3 Euler angles in the end frame
// euler_vel: 3 Euler velocities in the end frame
//
// 1 output:
// R_dot: time derivative of the rotation matrix
//------------------------------------------------------------------------
Matrix3d rot_matrix_dev(const Vector3d& euler_ang, const Vector3d& euler_vel) {
    // Euler angles
    double theta_x = euler_ang[0];
    double theta_y = euler_ang[1];
    double theta_z = euler_ang[2];

    // Euler velocities
    double dtheta_x = euler_vel[0];
    double dtheta_y = euler_vel[1];
    double dtheta_z = euler_vel[2];

    // Precompute sines and cosines
    double cx = cos(theta_x), sx = sin(theta_x);
    double cy = cos(theta_y), sy = sin(theta_y);
    double cz = cos(theta_z), sz = sin(theta_z);

    // Compute the time derivative of the rotation matrix R
    Matrix3d R_dot;
    R_dot << -sz * cy * dtheta_z - cz * sy * dtheta_y,
            -sz * sy * sx * dtheta_z - cz * (cy * sx * dtheta_y + sy * cx * dtheta_x) - cz * cx * dtheta_z,
            -sz * sy * cx * dtheta_z - cz * (cy * cx * dtheta_y - sy * sx * dtheta_x) + cz * sx * dtheta_z,

            cz * cy * dtheta_z - sz * sy * sx * dtheta_z,
            cz * sy * sx * dtheta_z - sz * (cy * sx * dtheta_y + sy * cx * dtheta_x) + sz * cx * dtheta_z,
            cz * sy * cx * dtheta_z - sz * (cy * cx * dtheta_y - sy * sx * dtheta_x) - sz * sx * dtheta_z,

            0, -cy * sx * dtheta_y - sy * cx * dtheta_x, -cy * cx * dtheta_y + sy * sx * dtheta_x;

    return R_dot;
}


//----------------------------------------------------------------
// Function to compute angular acceleration in the base frame
//----------------------------------------------------------------
// 3 inputs:
// euler_ang: 3 Euler angles in the end frame
// euler_vel: 3 Euler velocities in the end frame
// euler_acc: 3 Euler accelerations in the end frame
//
// 1 output:
// alpha_base: angular accelerations in the base frame
//-----------------------------------------------------------------
Vector3d euler_acc_to_base_angular_acc(
        const Vector3d& euler_ang,       // [theta_x, theta_y, theta_z]
        const Vector3d& euler_vel,       // [dtheta_x, dtheta_y, dtheta_z]
        const Vector3d& euler_acc        // [ddtheta_x, ddtheta_y, ddtheta_z]
) {
    // Compute the transformation matrix T_Euler^ZYX
    Matrix3d T_euler_ZYX;
    double theta_x = euler_ang[0];
    double theta_y = euler_ang[1];

    T_euler_ZYX <<
                1, 0, -sin(theta_y),
            0, cos(theta_x), cos(theta_y) * sin(theta_x),
            0, -sin(theta_x), cos(theta_y) * cos(theta_x);

    // Angular velocity in the end frame
    Vector3d omega_base = T_euler_ZYX * euler_vel;

    // Compute the time derivative of the angular velocity in the end frame
    Vector3d alpha_end = T_euler_ZYX * euler_acc; // First part

    // Compute the time derivative of the rotation matrix
    Matrix3d R_dot = rot_matrix_dev(euler_ang, euler_vel);

    // Compute the rotation matrix from ZYX Euler angles
    Matrix3d R = euler_to_rotation_matrix(euler_ang[0], euler_ang[1], euler_ang[2]);

    // Compute angular acceleration in the base frame
    Vector3d alpha_base = R * alpha_end + R_dot * omega_base;

    return alpha_base;
}



//------------------------------------------------------------
// Scale function for the friction compensation (to avoid rotation with a small velocity fluctuation around 0)
//------------------------------------------------------------
// 2 inputs:
// k: function slope
// dq: joint angular velocity
//
// output: scale factor for firction compensation
//----------------------------------------------------------
double scale_fric(const double k, double dq) {
    return 2.0 / (1.0 + exp(-k * dq)) - 1.0;
}


//-----------------------------------------------------------
// Function to calculate velocity gains for friction compensation
//-----------------------------------------------------------
// 2 inputs:
// dq: joint angular velocity
// index: number of joint
//
// output: joint velocity gain to compensate friction
//-----------------------------------------------------------
double vel_gain_fric(double dq, int index){
    VectorXd gain_posi(7), gain_nega(7);
    // Gains for 7 joints with different signs (positive or negative)
 //   gain_posi << 9, 8, 9, 8.5, 9, 8, 9.5;
   // gain_nega << 9, 8, 9, 8.5, 9, 8, 9;
    gain_posi << 9, 8, 9, 8.5, 8, 6, 5;
    gain_nega << 9, 8, 9, 8.5, 8, 6, 5;

    if (dq >= 0){
        return gain_posi(index);
    }
    else{
        return gain_nega(index);
    }
}

//-----------------------------------------------------------
// Function to compute position and orientation errors
//-----------------------------------------------------------
// 5 inputs:
// pos_end: XYZ positions of end effector
// pos_d: desired XYZ position of end effector
// rot_end: rotation matrix of end effector
// rot_d: desired rotation matrix of end effector
// e: 6-DoF errors used in the controller
//----------------------------------------------------------
void pos_difference(Vector3d& pos_end, Vector3d& pos_d, Matrix3d& rot_end, Matrix3d& rot_d, VectorXd& e){

    // 3-DoF XYZ position errors
    Vector3d pos_e;
    pos_e = pos_end - pos_d;
    e.head<3>() = pos_e;

    // 3-DoF orientation errors
    Quaterniond orientation(rot_end);
    Quaterniond orientation_d(rot_d);
    if (orientation_d.coeffs().dot(orientation.coeffs()) < 0.0) {
        orientation.coeffs() << -orientation.coeffs();
    }
    // Quaternion difference
    Quaterniond error_quaternion(orientation.inverse() * orientation_d);
    e.tail(3) << error_quaternion.x(), error_quaternion.y(), error_quaternion.z();
    // Transform to base frame
    e.tail(3) << 2 * (-rot_end * e.tail(3));
}

void friction_compensation(VectorXd& u, VectorXd& dq){
    VectorXd offset(7), gain(7);
    offset << 0.5, 0.6, 0.45, 0.4, 0.4, 0.3, 0.5;

    u(0) += scale_fric(20, dq(0)) * (offset(0) + vel_gain_fric(dq(0), 0) * abs(dq(0)));
    u(1) += scale_fric(20, dq(1)) * (offset(1) + vel_gain_fric(dq(1),1) * abs(dq(1)));
    u(2) += scale_fric(20, dq(2)) * (offset(2) + vel_gain_fric(dq(2),2) * abs(dq(2)));
    u(3) += scale_fric(20, dq(3)) * (offset(3) + vel_gain_fric(dq(3),3) * abs(dq(3)));
    u(4) += 1.1 + scale_fric(40, dq(4)) * (offset(4) + vel_gain_fric(dq(4),4) * abs(dq(4)));
    u(5) += scale_fric(40, dq(5)) * (offset(5) + vel_gain_fric(dq(5),5) * abs(dq(5)));
    u(6) += 1.2 + scale_fric(70, dq(6)) * (offset(6) + vel_gain_fric(dq(6),6) * abs(dq(6)));
}


//-----------------------------------------------------------
// Function to compute control input
//-----------------------------------------------------------
// 10 inputs:
// q,dq,ddq: joint angular positions, velocities and accelerations
// T_B7: rotation matrix of end effector
// p_d,dp_d, ddp_d: desire positions, velocities and accelerations
// K_d_diag: diagonal values of desired stiffness
// K_n_diag: diagonal values of nullspace stiffness
// c_f: control frequency
// time period: time difference between two control loop
//
// 3 outputs:
// u: control input
// dp: velocities of end effector
// ddp: accelerations of end effector
//----------------------------------------------------------
tuple<VectorXd, VectorXd, VectorXd, VectorXd> Controller::impedance_controller(VectorXd& q, VectorXd& dq, VectorXd& ddq, MatrixXd& T_B7,
                                                           VectorXd& p_d, VectorXd& dp_d, VectorXd& ddp_d,
                                                           VectorXd& K_d_diag, VectorXd& K_n_diag,
                                                           int c_f, double& time_period, VectorXd& actuator_torque,
                                                           double horizon, VectorXd std_p,
                                                           double sensor){
    // Outputs
    VectorXd u(7), dp(6), ddp(6);

    // Compute the stiffnesses and dampings
    Vector<double, 6> D_d_diag = 2 * K_d_diag.array().sqrt();
    Vector<double, 7> D_n_diag = 2 * K_n_diag.array().sqrt();
    DiagonalMatrix<double, 6> D_d = D_d_diag.asDiagonal();
    DiagonalMatrix<double, 6> K_d = K_d_diag.asDiagonal();
    DiagonalMatrix<double, 7> D_n = D_n_diag.asDiagonal();
    DiagonalMatrix<double, 7> K_n = K_n_diag.asDiagonal();

    // Define the acceleration amplification factor
    VectorXd acc_amp_factor_diag(6);
    DiagonalMatrix<double, Dynamic> acc_amp_factor(6);
    acc_amp_factor.setZero(); // Set all diagonal elements to zero

    // Parameters for Null space
    MatrixXd Iden = MatrixXd::Identity(7, 7);
    MatrixXd Null(7, 7);   // Nullspace projection matrix
    VectorXd e_q(7);          // Joint position errors in nullspace
    VectorXd q_d(7); // null space desired joint angles

    // Joint torques limits
    VectorXd gen3_JointTorquesLimits(7);
    gen3_JointTorquesLimits << 55, 55, 55, 55, 25, 25, 25;

    // Define the model matrices used in controller
    MatrixXd mass_matrix(7,7);
    VectorXd gravity_matrix(7);
    MatrixXd coriolis_matrix(7,7);
    MatrixXd jacobian(6,7);
    MatrixXd jacobian_dot(6, 7);
    MatrixXd pinv_jacobian(7, 6);
    MatrixXd jacobian_T(7, 6);
    MatrixXd pinv_jacobian_T(6, 7);
    MatrixXd M_x(6,6), C_x(6, 6);
    Matrix3d rot_end(3,3), rot_d(3,3);
    Vector3d pos_d(3), pos_end(3), pos_e(3), rot_e(3);
    VectorXd e(6), e_dot(6), g_x(6);
    VectorXd force(6), tau1(7), tau2(7), tau3(7);
    VectorXd model_e(7);
    VectorXd D_f_diag(6), K_f_diag(6);
    VectorXd acc_amplification(6);
    VectorXd torque_amp_diag(6);

    // Dynamic matrices calculations
    mass_matrix = mass_m(q);
    gravity_matrix = gravity_m(q);
    coriolis_matrix = coriolis_m(q,dq);
    jacobian = jaco_m(q);
    jacobian_dot = (jacobian - last_jaco) / time_period;
    last_jaco = jacobian;
    pinv_jacobian = (jacobian).completeOrthogonalDecomposition().pseudoInverse();
    jacobian_T = jacobian.transpose();
    pinv_jacobian_T = pinv_jacobian.transpose();

    // Dynamic matrices in task space
    M_x = pinv_jacobian_T * mass_matrix * pinv_jacobian;
    C_x = pinv_jacobian_T * (coriolis_matrix - mass_matrix * pinv_jacobian * jacobian_dot) * pinv_jacobian;
    g_x = pinv_jacobian_T * gravity_matrix;

    // Transformation from joint to task space
    dp = jacobian * dq;
    pos_end = T_B7.block<3, 1>(0, 3);
    rot_end = T_B7.block<3, 3>(0, 0);

    // Position error
    pos_d = p_d.head<3>();
    rot_d = euler_to_rotation_matrix(p_d(3), p_d(4), p_d(5));
    pos_difference(pos_end, pos_d, rot_end, rot_d, e);

    // Compute acceleration
    VectorXd pos_error(6), current_vel(6);
    pos_difference(pos_end, last_pos, rot_end, last_rot, pos_error);
    last_pos = pos_end;
    last_rot = rot_end;
    current_vel = pos_error * c_f;
    ddp = (current_vel - last_vel) * c_f;
    last_vel = current_vel;
    ddp = butterworth_filter(ddp);;

    // Transform the desired velocity and acceleration into world frame
    Vector3d ang_rate_world = euler_vel_to_base_angular_vel(p_d.tail<3>(),dp_d.tail<3>());
    dp_d.tail<3>() = ang_rate_world;
    Vector3d ang_acc_world = euler_acc_to_base_angular_acc(p_d.tail<3>(), dp_d.tail<3>(), ddp_d.tail<3>());
    ddp_d.tail<3>() = ang_acc_world;

    // Velocity error
    e_dot = dp - dp_d;

    // Tool gravity compensation (change with tool gravity)
    g_x(2) += 0.367*9.81;

    VectorXd dq_d(7);
    force = M_x*ddp_d + C_x*dp + g_x - K_d * e - D_d * e_dot;

    // Control law 1
    tau1 = jacobian_T * force;

    // Nullspace (Control law 2)
    dq_d = pinv_jacobian * dp_d;
    q_d = q + dq_d * (1 / c_f);
    q_d[1] = 0.26;
    q_d[2] = 3.14;
    Null = Iden - jacobian_T*pinv_jacobian_T;
    e_q = q_d - q;
    tau2 = Null*(K_n*e_q - D_n*dq);

    // Control input
    u = tau1 + tau2;

    friction_compensation(u, dq);

    // Set the torque saturation
    for (int i = 0; i < 7; i++)
    {
        if(u[i] > gen3_JointTorquesLimits[i]){
            u[i] = gen3_JointTorquesLimits[i];
        }else if(u[i] < -gen3_JointTorquesLimits[i]){
            u[i] = -gen3_JointTorquesLimits[i];
        }
    }
    // Save the data
    // We want to record: p, p_d, dp, dp_d, ddp, torque?
    // t, p, dp, ddp, torque
    /*
    if (log_file.is_open()) {
        VectorXd p(6);
        tie(p, T_B7) = forward(q);

        //std::cout << "Currently recording "<< recording_time <<"\n";
        log_file
        << recording_time << "," << p(0) << "," << p(1) << ","
        << p(2) << "," << p(3) << "," << p(4) << ","
        << p(5) << ",";
        log_file
        << dp(0) << "," << dp(1) << ","
        << dp(2) << "," << dp(3) << "," << dp(4) << ","
        << dp(5) << ",";
        log_file
        << ddp(0) << "," << ddp(1) << ","
        << ddp(2) << "," << ddp(3) << "," << ddp(4) << ","
        << ddp(5) << ",";
        log_file
        << actuator_torque(0) << "," << actuator_torque(1) << ","
        << actuator_torque(2) << "," << actuator_torque(3) << "," << actuator_torque(4) << ","
        << actuator_torque(5) << ",";
        log_file
        << p_d(0) << "," << p_d(1) << ","
        << p_d(2) << "," << p_d(3) << "," << p_d(4) << ","
        << p_d(5) << ",";
        log_file
        << dp_d(0) << "," << dp_d(1) << ","
        << dp_d(2) << "," << dp_d(3) << "," << dp_d(4) << ","
        << dp_d(5) << ",";


        log_file
        << horizon << ",";
        log_file
       << std_p(0) << "," << std_p(1) << ","
       << std_p(2) << "," ;

        log_file
        << sensor << ",";

        log_file << "\n";


    }
    */

    return make_tuple(u, dp, ddp, acc_factor_buffer);
}


//-----------------------------------------------------------
// Function to initilize the controller
//-----------------------------------------------------------
// 2 inputs:
// pos: initial position
// T_B7: initial rotation matrix
//-----------------------------------------------------------

void Controller::ini_controller(Vector3d& pos, MatrixXd& T_B7, std::string model_dir, int max_obs){
    last_pos << pos;
    last_rot = T_B7;
    last_vel = Vector<double, 6>::Zero();
    std::string filename = "C:/Imperial/Final project/predictive_controller/r" + std::to_string(file_num) + ".csv";
    //log_file.open("C:/Imperial/Final project/admittance control/m.csv");
    log_file.open(filename);
    /*log_file << "fx,fy,fz,tx,ty,tz,"
         << "ffiltx,ffilty,ffiltz,tfiltx,tfilty,tfiltz,"
         << "u_raw_1,u_raw_2,u_raw_3,u_raw_4,u_raw_5,u_raw_6,u_raw_7,"
         << "u_f_1,u_f_2,u_f_3,u_f_4,u_f_5,u_f_6,u_f_7\n";
    */

    promp_rt::OnlinePredictor predictor(model_dir, max_obs);
    kalman_predictor = MultiDofUkfPredictor(3,initial_positions);



}


//
// FUNCTION TO TRY AN ADMITTANCE CONTROLLER
//
// 10 inputs:
// q,dq,ddq: joint angular positions, velocities and accelerations
// T_B7: rotation matrix of end effector
// p_d,dp_d, ddp_d: desire positions, velocities and accelerations
// K_d_diag: diagonal values of desired stiffness
// D_d_diag: diagonal values of desired dampening
// I_d_diag: diagonal values of desired inertia
// actuator_torque: external force at the end effector
// c_f: control frequency
// time period: time difference between two control loop
//
// 3 outputs:
// u: control input
// dp: velocities of end effector
// ddp: accelerations of end effector
tuple<VectorXd, VectorXd, VectorXd, VectorXd> Controller::admittance_controller(VectorXd& q, VectorXd& dq, VectorXd& ddq, MatrixXd& T_B7,
                                                           VectorXd& K_d_diag, VectorXd& D_d_diag, VectorXd& I_d_diag,
                                                           VectorXd& actuator_torque,
                                                           int c_f, double& time_period,
                                                           VectorXd& K_n_diag){
    // Outputs
    VectorXd u(7), dp(6), ddp(6);

    // Compute the stiffnesses and dampings
    //Vector<double, 6> D_d_diag = 2 * K_d_diag.array().sqrt();
    Vector<double, 7> D_n_diag = 2 * K_n_diag.array().sqrt();
    //D_n_diag[1] = 30;
    //D_n_diag[0] = 30;
    DiagonalMatrix<double, 6> D_d = D_d_diag.asDiagonal();
    DiagonalMatrix<double, 6> K_d = K_d_diag.asDiagonal();
    DiagonalMatrix<double, 6> I_d = I_d_diag.asDiagonal();
    DiagonalMatrix<double, 7> D_n = D_n_diag.asDiagonal();
    DiagonalMatrix<double, 7> K_n = K_n_diag.asDiagonal();


    // INVERSE THE DESIRED INERTIA MATRIX
    Vector<double, 6> I_d_diag_inv;

    for (size_t i = 0; i < 6; i++) {
        if (I_d_diag[i] == 0) {
            cout << "Matrix is singular.\n";
        }
        I_d_diag_inv[i] = 1.0 / I_d_diag[i];
    }
    DiagonalMatrix<double, 6> I_d_inv = I_d_diag_inv.asDiagonal();


    // Define the acceleration amplification factor
    VectorXd acc_amp_factor_diag(6);
    DiagonalMatrix<double, Dynamic> acc_amp_factor(6);
    acc_amp_factor.setZero(); // Set all diagonal elements to zero

    // Joint torques limits
    VectorXd gen3_JointTorquesLimits(7);
    gen3_JointTorquesLimits << 55, 55, 55, 55, 25, 25, 25;

    // Define the model matrices used in controller
    MatrixXd mass_matrix(7,7);
    VectorXd gravity_matrix(7);
    MatrixXd coriolis_matrix(7,7);
    MatrixXd jacobian(6,7);
    MatrixXd jacobian_dot(6, 7);
    MatrixXd pinv_jacobian(7, 6);
    MatrixXd jacobian_T(7, 6);
    MatrixXd pinv_jacobian_T(6, 7);
    MatrixXd M_x(6,6), C_x(6, 6);
    Matrix3d rot_end(3,3), rot_d(3,3);
    Vector3d pos_d(3), pos_end(3), pos_e(3), rot_e(3);
    VectorXd e(6), e_dot(6), g_x(6);
    VectorXd force(6), tau1(7), tau2(7), tau3(7);
    VectorXd model_e(7);
    VectorXd D_f_diag(6), K_f_diag(6);
    VectorXd acc_amplification(6);
    VectorXd torque_amp_diag(6);

    // Dynamic matrices calculations
    mass_matrix = mass_m(q);
    gravity_matrix = gravity_m(q);
    coriolis_matrix = coriolis_m(q,dq);
    jacobian = jaco_m(q);
    jacobian_dot = (jacobian - last_jaco) / time_period;
    last_jaco = jacobian;
    pinv_jacobian = (jacobian).completeOrthogonalDecomposition().pseudoInverse();
    jacobian_T = jacobian.transpose();
    pinv_jacobian_T = pinv_jacobian.transpose();

    // Dynamic matrices in task space
    M_x = pinv_jacobian_T * mass_matrix * pinv_jacobian;
    C_x = pinv_jacobian_T * (coriolis_matrix - mass_matrix * pinv_jacobian * jacobian_dot) * pinv_jacobian;
    g_x = pinv_jacobian_T * gravity_matrix;

    // Transformation from joint to task space
    dp = jacobian * dq;
    pos_end = T_B7.block<3, 1>(0, 3);
    rot_end = T_B7.block<3, 3>(0, 0);

    // Position error
    VectorXd p_d(6); 
    p_d << 0.442,-0.06,0.575,M_PI/2,0,M_PI/2;
    pos_d = p_d.head<3>();
    rot_d = euler_to_rotation_matrix(p_d(3), p_d(4), p_d(5));
    pos_difference(pos_end, pos_d, rot_end, rot_d, e);

    // Compute acceleration
    VectorXd pos_error(6), current_vel(6);
    pos_difference(pos_end, last_pos, rot_end, last_rot, pos_error);
    last_pos = pos_end;
    last_rot = rot_end;
    current_vel = pos_error * c_f;
    ddp = (current_vel - last_vel) * c_f;
    last_vel = current_vel;
    ddp = butterworth_filter(ddp);

    // Tool gravity compensation (change with tool gravity)
    g_x(2) += 0.367*9.81;

    VectorXd f_ext(6);
    actuator_torque[6] = 0.0;
    f_ext = pinv_jacobian_T*actuator_torque + g_x;

    VectorXd dq_d(7);

    //f_ext.tail(3).setZero();
    VectorXd f_filt(7);

    f_filt = b0*f_ext+b1*f_ext_1+b2*f_ext_2
                -a1*f_filt_1-a2*f_filt_2;
    f_ext_2 = f_ext_1;
    f_ext_1 = f_ext;

    f_filt_2 = f_filt_1;
    f_filt_1 = f_filt;
    //D_d_diag.tail(3).setZero();

    //f_filt = f_ext;

    //force = f_filt + C_x*dp + g_x + M_x*I_d_inv*(f_filt - D_d*dp) - K_d * e;
    force = 1.5*f_filt + C_x*dp + g_x  - K_d * e;
    //force = f_filt + C_x*dp + g_x + M_x*I_d_inv*(f_filt - D_d*dp - K_d * e);
    //force = f_ext + C_x*dp + g_x + M_x*I_d_inv*(f_ext - D_d*dp) - K_d * e;
    //std::cout << "[INFO] Error " << f_ext << "\n";
    //force = M_x*ddp_d + C_x*dp + g_x - K_d * e - D_d * e_dot;

    // Control law 1
    tau1 = jacobian_T * force;

    // Control input

    //u  = jacobian_T *g_x+ actuator_torque;
    //u  = jacobian_T * (f_filt + C_x*dp + g_x - K_d * e + M_x*I_d_inv*(f_filt-D_d*dp));
    //u  = jacobian_T * (M_x*I_d_inv*(f_filt-D_d*dp));
    //std::cout << "[INFO]  " << M_x*I_d_inv*(f_filt-D_d*dp) << "\n";
    //u << 0,0,0,0,0,0,1;
    VectorXd q_d(7);
    dq_d = pinv_jacobian * dp;
    q_d = q + dq_d * (1 / c_f);
    //std::cout << "q_d[1] = " << q_d[1] << std::endl;
    q_d[1] = 0.26;
    q_d[2] = 3.14;

    //q_d[0] = -0.087;
    MatrixXd Iden = MatrixXd::Identity(7, 7);
    MatrixXd Null(7, 7);   // Nullspace projection matrix
    Null = Iden - jacobian_T*pinv_jacobian_T;
    //Null = Iden - pinv_jacobian_T*jacobian_T;
    VectorXd e_q;
    e_q = q_d - q;
    tau2 = Null*(K_n*e_q - D_n*dq);

    u = tau1 +tau2;
    //u = tau1;

/*
    u_filtered = b0*u+b1*u_1+b2*u_2
                -a1*u_filtered_1-a2*u_filtered_2;
    u_2 = u_1;
    u_1 = u;

    u_filtered_2 = u_filtered_1;
    u_filtered_1 = u_filtered;


    //double alpha = 0.9;
    //u_filtered = alpha * u_raw + (1 - alpha) * u_filtered;

    u = u_filtered;
*/
    VectorXd u_raw = u;
    u_filtered = u;

    VectorXd p(6);
    tie(p, T_B7) = forward(q);

    //cout << "actuator_torque = " << actuator_torque(6) << endl;
    /*
    if (GetAsyncKeyState(VK_SPACE) & 0x8000){
    if (log_file.is_open()) {
        std::cout << "Currently recording "<< recording_time <<"\n";
        log_file
        << recording_time << "," << p(0) << "," << p(1) << ","
        << p(2) << "," << p(3) << "," << p(4) << ","
        << p(5) << ",";
        log_file << "\n";
        recording_time = recording_time + time_period;

        file_change_flag =0;
    }
    }
    if (file_change_flag == 0){
    if (GetAsyncKeyState('Z') & 0x8000) {
        file_num = file_num + 1;
        recording_time = 0;
        std::cout << "Changing file to p"<< file_num <<"\n";
        log_file.close();

        std::string filename = "C:/Imperial/Final project/admittance control/p" + std::to_string(file_num) + ".csv";
        //log_file.open("C:/Imperial/Final project/admittance control/m.csv");
        log_file.open(filename);
        file_change_flag =1;
    }
    }
    */
    /*
    if (log_file.is_open()){
    log_file
        << f_ext(0) << "," << f_ext(1) << "," << f_ext(2) << ","
        << f_ext(3) << "," << f_ext(4) << "," << f_ext(5) << ",";
    log_file
        << f_filt(0) << "," << f_filt(1) << "," << f_filt(2) << ","
        << f_filt(3) << "," << f_filt(4) << "," << f_filt(5) << ",";

    for (int i = 0; i < 7; i++)
        log_file << actuator_torque(i) << ",";

    for (int i = 0; i < 7; i++)
        log_file << u_filtered(i) << ",";

    log_file << "\n";
    }
    */

    friction_compensation(u, dq);

    // Set the torque saturation
    for (int i = 0; i < 7; i++)
    {
        if(u[i] > gen3_JointTorquesLimits[i]){
            u[i] = gen3_JointTorquesLimits[i];
        }else if(u[i] < -gen3_JointTorquesLimits[i]){
            u[i] = -gen3_JointTorquesLimits[i];
        }
    }

    return make_tuple(u, dp, ddp, acc_factor_buffer);
}

// 10 inputs:
// q,dq,ddq: joint angular positions, velocities and accelerations
// T_B7: rotation matrix of end effector
// p_d,dp_d, ddp_d: desire positions, velocities and accelerations
// K_d_diag: diagonal values of desired stiffness
// D_d_diag: diagonal values of desired dampening
// I_d_diag: diagonal values of desired inertia
// f_ext: external force at the end effector
// c_f: control frequency
// time period: time difference between two control loop
//
// 3 outputs:
// u: control input
// dp: velocities of end effector
// ddp: accelerations of end effector
tuple<VectorXd, VectorXd, VectorXd, VectorXd> Controller::position_admittance_controller(VectorXd& q, VectorXd& dq, VectorXd& ddq, MatrixXd& T_B7,
                                                           VectorXd& K_d_diag, VectorXd& D_d_diag, VectorXd& I_d_diag,
                                                           VectorXd& actuator_torque,
                                                           int c_f, double& time_period){
    // Outputs
    VectorXd u(7), dp(6), ddp(6);

    // Compute the stiffnesses and dampings
    //Vector<double, 6> D_d_diag = 2 * K_d_diag.array().sqrt();
    //Vector<double, 7> D_n_diag = 2 * K_n_diag.array().sqrt();
    DiagonalMatrix<double, 6> D_d = D_d_diag.asDiagonal();
    DiagonalMatrix<double, 6> K_d = K_d_diag.asDiagonal();
    DiagonalMatrix<double, 6> I_d = I_d_diag.asDiagonal();


    // INVERSE THE DESIRED INERTIA MATRIX
    Vector<double, 6> I_d_diag_inv;

    for (size_t i = 0; i < 6; i++) {
        if (I_d_diag[i] == 0) {
            cout << "Matrix is singular.\n";
        }
        I_d_diag_inv[i] = 1.0 / I_d_diag[i];
    }
    DiagonalMatrix<double, 6> I_d_inv = I_d_diag_inv.asDiagonal();


    // Define the acceleration amplification factor
    VectorXd acc_amp_factor_diag(6);
    DiagonalMatrix<double, Dynamic> acc_amp_factor(6);
    acc_amp_factor.setZero(); // Set all diagonal elements to zero

    // Joint torques limits
    VectorXd gen3_JointTorquesLimits(7);
    gen3_JointTorquesLimits << 55, 55, 55, 55, 25, 25, 25;

    // Define the model matrices used in controller
    MatrixXd mass_matrix(7,7);
    VectorXd gravity_matrix(7);
    MatrixXd coriolis_matrix(7,7);
    MatrixXd jacobian(6,7);
    MatrixXd jacobian_dot(6, 7);
    MatrixXd pinv_jacobian(7, 6);
    MatrixXd jacobian_T(7, 6);
    MatrixXd pinv_jacobian_T(6, 7);
    MatrixXd M_x(6,6), C_x(6, 6);
    Matrix3d rot_end(3,3), rot_d(3,3);
    Vector3d pos_d(3), pos_end(3), pos_e(3), rot_e(3);
    VectorXd e(6), e_dot(6), g_x(6);
    VectorXd force(6), tau1(7), tau2(7), tau3(7);
    VectorXd model_e(7);
    VectorXd D_f_diag(6), K_f_diag(6);
    VectorXd acc_amplification(6);
    VectorXd torque_amp_diag(6);

    // Dynamic matrices calculations
    mass_matrix = mass_m(q);
    gravity_matrix = gravity_m(q);
    coriolis_matrix = coriolis_m(q,dq);
    jacobian = jaco_m(q);
    jacobian_dot = (jacobian - last_jaco) / time_period;
    last_jaco = jacobian;
    pinv_jacobian = (jacobian).completeOrthogonalDecomposition().pseudoInverse();
    jacobian_T = jacobian.transpose();
    pinv_jacobian_T = pinv_jacobian.transpose();

    // Dynamic matrices in task space
    M_x = pinv_jacobian_T * mass_matrix * pinv_jacobian;
    C_x = pinv_jacobian_T * (coriolis_matrix - mass_matrix * pinv_jacobian * jacobian_dot) * pinv_jacobian;
    g_x = pinv_jacobian_T * gravity_matrix;

    // Transformation from joint to task space
    dp = jacobian * dq;
    pos_end = T_B7.block<3, 1>(0, 3);
    rot_end = T_B7.block<3, 3>(0, 0);

    // Position error
    VectorXd p_d(6);
    p_d << 0.442,-0.06,0.575,M_PI/2,0,M_PI/2;
    pos_d = p_d.head<3>();
    rot_d = euler_to_rotation_matrix(p_d(3), p_d(4), p_d(5));
    pos_difference(pos_end, pos_d, rot_end, rot_d, e);

    // Compute acceleration
    VectorXd pos_error(6), current_vel(6);
    pos_difference(pos_end, last_pos, rot_end, last_rot, pos_error);
    last_pos = pos_end;
    last_rot = rot_end;
    current_vel = pos_error * c_f;
    ddp = (current_vel - last_vel) * c_f;
    last_vel = current_vel;
    ddp = butterworth_filter(ddp);

    // Tool gravity compensation (change with tool gravity)
    g_x(2) += 0.1*9.81;

    VectorXd f_ext(6);
    actuator_torque[6] = 0.0;
    f_ext = pinv_jacobian_T*actuator_torque + g_x;

    VectorXd dq_d(7);

    //f_ext.tail(3).setZero();
    VectorXd f_filt(7);

    f_filt = b0*f_ext+b1*f_ext_1+b2*f_ext_2
                -a1*f_filt_1-a2*f_filt_2;
    f_ext_2 = f_ext_1;
    f_ext_1 = f_ext;

    f_filt_2 = f_filt_1;
    f_filt_1 = f_filt;
    //D_d_diag.tail(3).setZero();

    //f_filt = f_ext;

    //force = f_filt + C_x*dp + g_x + M_x*I_d_inv*(f_filt - D_d*dp) - K_d * e;
    //force = 1.5*f_filt + C_x*dp + g_x  - K_d * e;
    //force = f_filt + C_x*dp + g_x + M_x*I_d_inv*(f_filt - D_d*dp - K_d * e);
    //force = f_ext + C_x*dp + g_x + M_x*I_d_inv*(f_ext - D_d*dp) - K_d * e;
    //std::cout << "[INFO] Error " << f_ext << "\n";
    //force = M_x*ddp_d + C_x*dp + g_x - K_d * e - D_d * e_dot;

    VectorXd ddu(6);
    VectorXd du(6);
    VectorXd du_jointspace(7);
    ddu = I_d_inv*(f_filt - D_d*dp - K_d * e);
    du = dp+ddu*time_period;
    du_jointspace = pinv_jacobian*du;
    u = q + du_jointspace*time_period;
    u[4] = 0;
    u[5] = 1;
    u[6] = M_PI/2;

    std::cout << "[INFO] ddu=  " << I_d_inv*(f_filt-D_d*dp - K_d * e) << "\n";
    std::cout << "[INFO] u=  " << u << "\n";
    //u << 0,0,0,0,0,0,1;

    VectorXd u_raw = u;
    u_filtered = u;

    //cout << "actuator_torque = " << actuator_torque(6) << endl;



    if (log_file.is_open()){
    log_file
        << f_ext(0) << "," << f_ext(1) << "," << f_ext(2) << ","
        << f_ext(3) << "," << f_ext(4) << "," << f_ext(5) << ",";
    log_file
        << f_filt(0) << "," << f_filt(1) << "," << f_filt(2) << ","
        << f_filt(3) << "," << f_filt(4) << "," << f_filt(5) << ",";

    for (int i = 0; i < 7; i++)
        log_file << actuator_torque(i) << ",";

    for (int i = 0; i < 7; i++)
        log_file << u_filtered(i) << ",";

    log_file << "\n";
    }

    //friction_compensation(u, dq);

    // Set the torque saturation
    for (int i = 0; i < 7; i++)
    {
        if(u[i] > gen3_JointTorquesLimits[i]){
            u[i] = gen3_JointTorquesLimits[i];
        }else if(u[i] < -gen3_JointTorquesLimits[i]){
            u[i] = -gen3_JointTorquesLimits[i];
        }
    }

    return make_tuple(u, dp, ddp, acc_factor_buffer);
}

int baseline_counter = 0;

tuple<VectorXd, VectorXd, VectorXd, VectorXd, int, int, int> Controller::predictive_impedance_controller(VectorXd& q, VectorXd& dq, VectorXd& ddq, MatrixXd& T_B7,
                                                           VectorXd& p_d, VectorXd& dp_d, VectorXd& ddp_d,
                                                           VectorXd& K_d_diag,  VectorXd& K_n_diag,
                                                           int c_f, double& time_period,
                                                           double override_duration, promp_rt::ConditioningMode mode,
                                                           double record_dt_s, std::vector<double> target_pos, VectorXd& K_d_diag_admittance,
                                                           VectorXd& D_d_diag, VectorXd& I_d_diag,
                                                           VectorXd& actuator_torque, VectorXd& current,
                                                           double sensor,
                                                           int flag_mode, int flag_reset, int flag_goal, int flag_running,
                                                           std::string filename)
{
    if (flag_running == 1)
    {
        VectorXd p(6);
        VectorXd std_p(3);
        tie(p, T_B7) = forward(q);
        vector<double> p2(p.data(), p.data() + p.size());
        MatrixXd jacobian(6,7);
        jacobian = jaco_m(q);
        VectorXd dp(6);
        dp = jacobian * dq;
        double total_velocity = sqrt(dp(0)*dp(0)+dp(1)*dp(1)+dp(2)*dp(2));
        double horizon = max(0.030, min(0.100, 0.01/(total_velocity+1e-5)));
        VectorXd f_ext(6);
        MatrixXd pinv_jacobian_T(6, 7);
        pinv_jacobian_T = (jacobian).completeOrthogonalDecomposition().pseudoInverse().transpose();
        f_ext = pinv_jacobian_T*actuator_torque;

        int goal_reached = 0;
        double threshold = 0.04;
        if ((flag_goal == 1)&&(abs(p(0)-0.805)<threshold)&&(abs(p(1)+0.06)<threshold)&&(abs(p(2)-0.28)<threshold))
        {
            goal_reached = 1;
        }else if ((flag_goal == 2)&&(abs(p(0)-0.695)<threshold)&&(abs(p(1)+0.067)<threshold)&&(abs(p(2)-0.407)<threshold))
        {
            goal_reached = 1;
        }else if ((flag_goal == 3)&&(abs(p(0)-0.56)<threshold)&&(abs(p(1)+0.08)<threshold)&&(abs(p(2)-0.44)<threshold))
        {
            goal_reached = 1;
        }else if ((flag_goal == 4)&&(abs(p(0)-0.77)<threshold)&&(abs(p(1)-0.067)<threshold)&&(abs(p(2)-0.51)<threshold))
        {
            goal_reached = 1;
        }else if ((flag_goal == 5)&&(abs(p(0)-0.72)<threshold)&&(abs(p(1)-0.27)<threshold)&&(abs(p(2)-0.36)<threshold))
        {
            goal_reached = 1;
        }


        if (flag_mode == 1)
        {
            if (log_file.is_open())
            {
                recording_time = recording_time + time_period;

                //std::cout << "Currently recording "<< recording_time <<"\n";
                log_file
                << recording_time << "," << p(0) << "," << p(1) << ","
                << p(2) << "," << p(3) << "," << p(4) << ","
                << p(5) << ",";
                log_file
                << dp(0) << "," << dp(1) << ","
                << dp(2) << "," << dp(3) << "," << dp(4) << ","
                << dp(5) << ",";
                log_file
                << f_ext(0) << "," << f_ext(1) << ","
                << f_ext(2) << "," << f_ext(3) << "," << f_ext(4) << ","
                << f_ext(5) << ",";
                log_file
                << p_d(0) << "," << p_d(1) << ","
                << p_d(2) << "," << p_d(3) << "," << p_d(4) << ","
                << p_d(5) << ",";
                log_file
                << dp_d(0) << "," << dp_d(1) << ","
                << dp_d(2) << "," << dp_d(3) << "," << dp_d(4) << ","
                << dp_d(5) << ",";


                log_file
                << horizon << ",";
                log_file
               << std_p(0) << "," << std_p(1) << ","
               << std_p(2) << "," ;

                log_file
                << sensor << ",";

                log_file
                << current(0) << "," << current(1) << ","
                << current(2) << "," << current(3) << "," << current(4) << ","
                << current(5) << ","<< current(6) << ",";

                log_file << "\n";
            }

            return std::tuple_cat(admittance_controller(q, dq, ddq, T_B7, K_d_diag_admittance, D_d_diag, I_d_diag,
                                                                actuator_torque, c_f, time_period, K_n_diag),
                                                                std::make_tuple(goal_reached, 0), std::make_tuple(0));
        }if (flag_mode == 6)
        {

            if (log_file.is_open())
            {
                recording_time = recording_time + time_period;

                //std::cout << "Currently recording "<< recording_time <<"\n";
                log_file
                << recording_time << "," << p(0) << "," << p(1) << ","
                << p(2) << "," << p(3) << "," << p(4) << ","
                << p(5) << ",";
                log_file
                << dp(0) << "," << dp(1) << ","
                << dp(2) << "," << dp(3) << "," << dp(4) << ","
                << dp(5) << ",";
                log_file
                << f_ext(0) << "," << f_ext(1) << ","
                << f_ext(2) << "," << f_ext(3) << "," << f_ext(4) << ","
                << f_ext(5) << ",";
                log_file
                << p_d(0) << "," << p_d(1) << ","
                << p_d(2) << "," << p_d(3) << "," << p_d(4) << ","
                << p_d(5) << ",";
                log_file
                << dp_d(0) << "," << dp_d(1) << ","
                << dp_d(2) << "," << dp_d(3) << "," << dp_d(4) << ","
                << dp_d(5) << ",";


                log_file
                << horizon << ",";
                log_file
               << std_p(0) << "," << std_p(1) << ","
               << std_p(2) << "," ;

                log_file
                << sensor << ",";

                log_file
                << current(0) << "," << current(1) << ","
                << current(2) << "," << current(3) << "," << current(4) << ","
                << current(5) << ","<< current(6) << ",";

                log_file << "\n";
            }

            return std::tuple_cat(admittance_controller(q, dq, ddq, T_B7, K_d_diag_admittance, D_d_diag, I_d_diag,
                                                                actuator_torque, c_f, time_period, K_n_diag),
                                                                std::make_tuple(goal_reached, 0), std::make_tuple(0));

        }else if (flag_mode == 2)
        {
            //std::cout << "Starting prediction" << "\n";
            horizon = 0.030;
            promp_rt::PredictionResult result =
                predictor.update_and_predict(recording_time, p2, mode, record_dt_s, target_pos, horizon);
            //std::cout << "Finished prediction" << "\n";
            //USE RESULT TO CREATE p_d, dp_d, ddp_d
            //std::cout << "Size of mean_traj: " << result.mean_traj.size()<< "\n";

            p_d[0] = result.mean_traj[0];
            //std::cout << p_d[0] << "\n";
            p_d[1] = result.mean_traj[2];
            //std::cout << p_d[1] << "\n";
            p_d[2] = result.mean_traj[4];
            //std::cout << p_d[2] << "\n";
            dp_d[0] = result.mean_traj[1];
            //std::cout << dp_d[0] << "\n";
            dp_d[1] = result.mean_traj[3];
            //std::cout << dp_d[1] << "\n";
            dp_d[2] = result.mean_traj[5];
            //std::cout << dp_d[2] << "\n";
            //std::cout << "Starting std_traj"<< "\n";
            //std::cout << "Size of std_traj: " << result.std_traj.size()<< "\n";
            std_p[0] = result.std_traj[0];
            std_p[1] = result.std_traj[1];
            std_p[2] = result.std_traj[2];

            file_change_flag = 0;
            recording_time = recording_time + time_period;
            /*
            double K = 300;
            double V_lim = 0.1;
            double K_max = 300;
            double K_min = 100;
            if (total_velocity < V_lim)
            {
                if (recording_time > 1)
                {
                    K = K_max - (K_max-K_min)*(V_lim-total_velocity)/V_lim;
                }else
                {
                    K = K_max;
                }
            }else
            {
                K = K_max;
            }
            K_d_diag << K,K,K,10,10,10;
            */
            double K = 1000;
            double V_lim = 0.08;
            double K_max = 1000;
            double K_min = 100;

                double Kx = 1000;
                double Ky = 1000;
                double Kz = 1000;

                if (abs(dp[0]) < 0.125)
                {
                    if (braking)
                    {
                        // Negative acceleration, stop it in its tracks
                        p_d[0] = p[0];
                        dp_d[0] = 0;
                    } // else, let it help you
                    if (abs(dp[0]) < 0.01)
                    {
                        braking = false;
                    }
                }else
                {
                    braking = true;
                }

            if (total_velocity < V_lim)
            {


                //p_d[0] = p[0];
                //std::cout << p_d[0] << " ";
                p_d[1] = p[1];
                //std::cout << p_d[1] << " ";
                p_d[2] = p[2];
                //std::cout << p_d[2] << " ";

                //std::cout << dp_d[0] << " ";
                dp_d[1] = 0;
                //std::cout << dp_d[1] << " ";
                dp_d[2] = 0;



                Kx = 1.5*K;
                Ky = K;
                Kz = 1.5*K;
            }else
            {
                Kx = max(1500 - 30*sensor, 10.0);
                Ky = max(1500 - 30*sensor, 10.0);
                Kz = max(1500 - 30*sensor, 10.0);
                /*
                K = K_max;
                if (f_ext(0)*dp(0)>0.1)
                {
                    Kx = min(750 + 12.5*abs(f_ext(0)), 2000.0);
                }else if (f_ext(0)*dp(0)< -0.1)
                {
                    Kx = max(750 - 7.4*abs(f_ext(0)), 10.0);
                }else
                {
                    Kx = 750;
                }

                if (f_ext(1)*dp(1)>0.1)
                {
                    Ky = min(500 + 30*sensor, 2000.0);
                }else if (f_ext(1)*dp(1)< -0.1)
                {
                    Ky = max(500 - 16.33*sensor, 10.0);
                }else
                {
                    Ky = 500;
                }

                if ((f_ext(2)-baseline_Fz)*dp(2)>0.3)
                {
                    Kz = min(1.5*500 + 87.5*sensor, 2500.0);
                }else if ((f_ext(2)-baseline_Fz)*dp(2)< -0.3)
                {
                    //Kz = min(1.5*500 + 87.5*sensor, 2500.0);
                    Kz = max(1.5*500 - 37*sensor, 10.0);

                }else
                {
                    Kz = 1.5*500;
                }
                //Kx = K;
                //Ky = K;
                //Kz = K;
                */

            }
            K_d_diag << Kx,Ky,Kz,30,30,30;
            //std::cout << "Starting impedance controller"<< "\n";

            if (log_file.is_open())
            {

                //std::cout << "Currently recording "<< recording_time <<"\n";
                log_file
                << recording_time << "," << p(0) << "," << p(1) << ","
                << p(2) << "," << p(3) << "," << p(4) << ","
                << p(5) << ",";
                log_file
                << dp(0) << "," << dp(1) << ","
                << dp(2) << "," << dp(3) << "," << dp(4) << ","
                << dp(5) << ",";
                log_file
                << f_ext(0) << "," << f_ext(1) << ","
                << f_ext(2) << "," << f_ext(3) << "," << f_ext(4) << ","
                << f_ext(5) << ",";
                log_file
                << p_d(0) << "," << p_d(1) << ","
                << p_d(2) << "," << p_d(3) << "," << p_d(4) << ","
                << p_d(5) << ",";
                log_file
                << dp_d(0) << "," << dp_d(1) << ","
                << dp_d(2) << "," << dp_d(3) << "," << dp_d(4) << ","
                << dp_d(5) << ",";


                log_file
                << horizon << ",";
                log_file
               << std_p(0) << "," << std_p(1) << ","
               << std_p(2) << "," ;

                log_file
                << sensor << ",";

                log_file
                << current(0) << "," << current(1) << ","
                << current(2) << "," << current(3) << "," << current(4) << ","
                << current(5) << ","<< current(6) << ",";

                log_file << "\n";
            }
            // IF GOAL REACHED, STOP
            if (goal_reached == 1)
            {
                K_d_diag << 800,800,800,40,40,40;
                p_d[0] = p[0];
                p_d[1] = p[1];
                p_d[2] = p[2];
                dp_d[0] = 0;
                dp_d[1] = 0;
                dp_d[2] = 0;
            }


            return std::tuple_cat(impedance_controller(q, dq, ddq, T_B7, p_d, dp_d,ddp_d, K_d_diag, K_n_diag,c_f, time_period,
                actuator_torque, horizon, std_p, sensor), std::make_tuple(goal_reached), std::make_tuple(0), std::make_tuple(0));

        }else if (flag_mode == 3)
        {
            //std::cout << sensor << "\n";
            try {
                std::vector<double> position = {p[0], p[1], p[2]};
                std::vector<double> velocity = {dp[0], dp[1], dp[2]};
                horizon = 0.040;
                auto forecast = kalman_predictor.updateAndPredict(position, recording_time, 0.01, horizon, velocity);
                //std::cout << "Exiting Kalman prediction ";
                //std::cout << forecast.size() << " ";

                p_d[0] = forecast[0].back();
                //std::cout << p_d[0] << " ";
                p_d[1] = forecast[1].back();
                //std::cout << p_d[1] << " ";
                p_d[2] = forecast[2].back();
                //std::cout << p_d[2] << " ";
                dp_d[0] = forecast[3].back();
                //std::cout << dp_d[0] << " ";
                dp_d[1] = forecast[4].back();
                //std::cout << dp_d[1] << " ";
                dp_d[2] = forecast[5].back();
                //std::cout << dp_d[2] << " ";

                //file_change_flag = 0;
            recording_time = recording_time + time_period;
            double K = 1000;
            double V_lim = 0.08;
            double K_max = 1000;
            double K_min = 100;

                double Kx = 1000;
                double Ky = 1000;
                double Kz = 1000;

                if (abs(dp[0]) < 0.125)
                {
                    if (braking)
                    {
                        // Negative acceleration, stop it in its tracks
                        p_d[0] = p[0];
                        dp_d[0] = 0;
                    } // else, let it help you
                    if (abs(dp[0]) < 0.01)
                    {
                        braking = false;
                    }
                }else
                {
                    braking = true;
                }

            if (total_velocity < V_lim)
            {


                //p_d[0] = p[0];
                //std::cout << p_d[0] << " ";
                p_d[1] = p[1];
                //std::cout << p_d[1] << " ";
                p_d[2] = p[2];
                //std::cout << p_d[2] << " ";

                //std::cout << dp_d[0] << " ";
                dp_d[1] = 0;
                //std::cout << dp_d[1] << " ";
                dp_d[2] = 0;



                Kx = 1500;
                Ky = 1.5*K;
                Kz = 1.5*K;
            }else
            {
                Kx = max(1500 - 30*sensor, 10.0);
                Ky = max(1500 - 30*sensor, 10.0);
                Kz = max(1500 - 30*sensor, 10.0);

                /*

                K = K_max;
                if (f_ext(0)*dp(0)>0.1)
                {
                    //Kx = min(750 + 12.5*abs(f_ext(0)), 2000.0);
                    Kx = max(1500 - 15*abs(f_ext(0)), 10.0);
                    //Kx = min(750 + 12.5*abs(f_ext(0)), 1500.0);
                }else if (f_ext(0)*dp(0)< -0.1)
                {
                    //Kx = max(750 - 7.4*abs(f_ext(0)), 10.0);
                    Kx = max(1500 - 15*abs(f_ext(0)), 10.0);
                    //Kx = min(750 + 12.5*abs(f_ext(0)), 1500.0);
                }else
                {
                    Kx = 750;
                }

                if (f_ext(1)*dp(1)>0.1)
                {
                    //Ky = min(500 + 30*sensor, 2000.0);
                    //Ky = max(500 - 16.33*sensor, 10.0);
                    //Ky = min(500 + 30*sensor, 1500.0);
                    Ky = max(1500 - 15*sensor, 10.0);
                }else if (f_ext(1)*dp(1)< -0.1)
                {
                    //Ky = max(500 - 16.33*sensor, 10.0);
                    //Ky = min(500 + 30*sensor, 1500.0);
                    Ky = max(1500 - 15*sensor, 10.0);
                }else
                {
                    Ky = 500;
                }

                if ((f_ext(2)-baseline_Fz)*dp(2)>0.3)
                {
                    //Kz = min(1.5*500 + 87.5*sensor, 2500.0);
                    //Kz = max(1.5*500 - 37*sensor, 10.0);
                    //Kz = min(1.5*500 + 87.5*sensor, 2000.0);
                    Kz = max(1500 - 15*sensor, 10.0);
                }else if ((f_ext(2)-baseline_Fz)*dp(2)< -0.3)
                {
                    //Kz = min(1.5*500 + 87.5*sensor, 2000.0);
                    //Kz = max(1.5*500 - 37*sensor, 10.0);
                    Kz = max(1500 - 15*sensor, 10.0);

                }else
                {
                    Kz = 1.5*500;
                }
                //Kx = K;
                //Ky = K;
                //Kz = K;

                */

            }
            K_d_diag << Kx,Ky,Kz,30,30,30;
            //std::cout << "Starting impedance controller"<< "\n";
            //    std::cout << baseline_Fz << "\n";

            if (log_file.is_open())
            {

                //std::cout << "Currently recording "<< recording_time <<"\n";
                log_file
                << recording_time << "," << p(0) << "," << p(1) << ","
                << p(2) << "," << p(3) << "," << p(4) << ","
                << p(5) << ",";
                log_file
                << dp(0) << "," << dp(1) << ","
                << dp(2) << "," << dp(3) << "," << dp(4) << ","
                << dp(5) << ",";
                log_file
                << f_ext(0) << "," << f_ext(1) << ","
                << f_ext(2) << "," << f_ext(3) << "," << f_ext(4) << ","
                << f_ext(5) << ",";
                log_file
                << p_d(0) << "," << p_d(1) << ","
                << p_d(2) << "," << p_d(3) << "," << p_d(4) << ","
                << p_d(5) << ",";
                log_file
                << dp_d(0) << "," << dp_d(1) << ","
                << dp_d(2) << "," << dp_d(3) << "," << dp_d(4) << ","
                << dp_d(5) << ",";


                log_file
                << horizon << ",";
                log_file
               << Kx << "," << Ky << ","
               << Kz << "," ;

                log_file
                << sensor << ",";

                log_file
                << current(0) << "," << current(1) << ","
                << current(2) << "," << current(3) << "," << current(4) << ","
                << current(5) << ","<< current(6) << ",";

                log_file << "\n";
            }
                // IF GOAL REACHED, STOP
                if (goal_reached == 1)
                {
                    K_d_diag << 800,800,800,40,40,40;
                    p_d[0] = p[0];
                    p_d[1] = p[1];
                    p_d[2] = p[2];
                    dp_d[0] = 0;
                    dp_d[1] = 0;
                    dp_d[2] = 0;
                }


            return std::tuple_cat(impedance_controller(q, dq, ddq, T_B7, p_d, dp_d,ddp_d, K_d_diag, K_n_diag,c_f, time_period,
                actuator_torque, horizon, std_p, sensor), std::make_tuple(goal_reached), std::make_tuple(0), std::make_tuple(0));
            } catch (const std::exception& e) {
                std::cout << "updateAndPredict threw: " << e.what() << std::endl;
            }



        }else if (flag_mode == 4)
        {
            //std::cout << "Starting prediction" << "\n";
            horizon = 0.040;
            promp_rt::PredictionResult result =
                predictor.update_and_predict(recording_time, p2, mode, record_dt_s, target_pos, horizon);
            //std::cout << "Finished prediction" << "\n";
            //USE RESULT TO CREATE p_d, dp_d, ddp_d
            //std::cout << "Size of mean_traj: " << result.mean_traj.size()<< "\n";

            p_d[0] = result.mean_traj[0];
            //std::cout << p_d[0] << "\n";
            p_d[1] = result.mean_traj[2];
            //std::cout << p_d[1] << "\n";
            p_d[2] = result.mean_traj[4];
            //std::cout << p_d[2] << "\n";
            dp_d[0] = result.mean_traj[1];
            //std::cout << dp_d[0] << "\n";
            dp_d[1] = result.mean_traj[3];
            //std::cout << dp_d[1] << "\n";
            dp_d[2] = result.mean_traj[5];
            //std::cout << dp_d[2] << "\n";
            //std::cout << "Starting std_traj"<< "\n";
            //std::cout << "Size of std_traj: " << result.std_traj.size()<< "\n";
            std_p[0] = result.std_traj[0];
            std_p[1] = result.std_traj[1];
            std_p[2] = result.std_traj[2];

            file_change_flag = 0;
            recording_time = recording_time + time_period;
            /*
            double K = 300;
            double V_lim = 0.1;
            double K_max = 300;
            double K_min = 100;
            if (total_velocity < V_lim)
            {
                if (recording_time > 1)
                {
                    K = K_max - (K_max-K_min)*(V_lim-total_velocity)/V_lim;
                }else
                {
                    K = K_max;
                }
            }else
            {
                K = K_max;
            }
            K_d_diag << K,K,K,10,10,10;
            */
            double K = 1000;
            double V_lim = 0.08;


                double Kx = 1000;
                double Ky = 1000;
                double Kz = 1000;

                if (abs(dp[0]) < 0.125)
                {
                    if (braking)
                    {
                        // Negative acceleration, stop it in its tracks
                        p_d[0] = p[0];
                        dp_d[0] = 0;
                    } // else, let it help you
                    if (abs(dp[0]) < 0.01)
                    {
                        braking = false;
                    }
                }else
                {
                    braking = true;
                }

            if (total_velocity < V_lim)
            {


                //p_d[0] = p[0];
                //std::cout << p_d[0] << " ";
                p_d[1] = p[1];
                //std::cout << p_d[1] << " ";
                p_d[2] = p[2];
                //std::cout << p_d[2] << " ";

                //std::cout << dp_d[0] << " ";
                dp_d[1] = 0;
                //std::cout << dp_d[1] << " ";
                dp_d[2] = 0;



                Kx = 1.5*K;
                Ky = 1.5*K;
                Kz = 1.5*K;
            }else
            {
                Kx = 1.5*K;
                Ky = 1.5*K;
                Kz = 1.5*K;
                //Kx = K;
                //Ky = K;
                //Kz = K;
            }
            K_d_diag << Kx,Ky,Kz,30,30,30;
            //std::cout << "Starting impedance controller"<< "\n";

            if (log_file.is_open())
            {

                //std::cout << "Currently recording "<< recording_time <<"\n";
                log_file
                << recording_time << "," << p(0) << "," << p(1) << ","
                << p(2) << "," << p(3) << "," << p(4) << ","
                << p(5) << ",";
                log_file
                << dp(0) << "," << dp(1) << ","
                << dp(2) << "," << dp(3) << "," << dp(4) << ","
                << dp(5) << ",";
                log_file
                << f_ext(0) << "," << f_ext(1) << ","
                << f_ext(2) << "," << f_ext(3) << "," << f_ext(4) << ","
                << f_ext(5) << ",";
                log_file
                << p_d(0) << "," << p_d(1) << ","
                << p_d(2) << "," << p_d(3) << "," << p_d(4) << ","
                << p_d(5) << ",";
                log_file
                << dp_d(0) << "," << dp_d(1) << ","
                << dp_d(2) << "," << dp_d(3) << "," << dp_d(4) << ","
                << dp_d(5) << ",";


                log_file
                << horizon << ",";
                log_file
               << std_p(0) << "," << std_p(1) << ","
               << std_p(2) << "," ;

                log_file
                << sensor << ",";

                log_file
                << current(0) << "," << current(1) << ","
                << current(2) << "," << current(3) << "," << current(4) << ","
                << current(5) << ","<< current(6) << ",";

                log_file << "\n";
            }
            // IF GOAL REACHED, STOP
            if (goal_reached == 1)
            {
                K_d_diag << 800,800,800,40,40,40;
                p_d[0] = p[0];
                p_d[1] = p[1];
                p_d[2] = p[2];
                dp_d[0] = 0;
                dp_d[1] = 0;
                dp_d[2] = 0;
            }


            return std::tuple_cat(impedance_controller(q, dq, ddq, T_B7, p_d, dp_d,ddp_d, K_d_diag, K_n_diag,c_f, time_period,
                actuator_torque, horizon, std_p, sensor), std::make_tuple(goal_reached), std::make_tuple(0), std::make_tuple(0));

        }else if (flag_mode == 5)
        {
            //std::cout << "Entering Kalman prediction ";
            try {
                std::vector<double> position = {p[0], p[1], p[2]};
                std::vector<double> velocity = {dp[0], dp[1], dp[2]};
                horizon = 0.040;
                auto forecast = kalman_predictor.updateAndPredict(position, recording_time, 0.01, horizon, velocity);
                //std::cout << "Exiting Kalman prediction ";
                //std::cout << forecast.size() << " ";

                p_d[0] = forecast[0].back();
                //std::cout << p_d[0] << " ";
                p_d[1] = forecast[1].back();
                //std::cout << p_d[1] << " ";
                p_d[2] = forecast[2].back();
                //std::cout << p_d[2] << " ";
                dp_d[0] = forecast[3].back();
                //std::cout << dp_d[0] << " ";
                dp_d[1] = forecast[4].back();
                //std::cout << dp_d[1] << " ";
                dp_d[2] = forecast[5].back();
                //std::cout << dp_d[2] << " ";

                //file_change_flag = 0;
            recording_time = recording_time + time_period;
            double K = 1000;
            double V_lim = 0.08;
                double Kx = 1000;
                double Ky = 1000;
                double Kz = 1000;

                if (abs(dp[0]) < 0.125)
                {
                    if (braking)
                    {
                        // Negative acceleration, stop it in its tracks
                        p_d[0] = p[0];
                        dp_d[0] = 0;
                    } // else, let it help you
                    if (abs(dp[0]) < 0.01)
                    {
                        braking = false;
                    }
                }else
                {
                    braking = true;
                }

            if (total_velocity < V_lim)
            {


                //p_d[0] = p[0];
                //std::cout << p_d[0] << " ";
                p_d[1] = p[1];
                //std::cout << p_d[1] << " ";
                p_d[2] = p[2];
                //std::cout << p_d[2] << " ";

                //std::cout << dp_d[0] << " ";
                dp_d[1] = 0;
                //std::cout << dp_d[1] << " ";
                dp_d[2] = 0;



                Kx = 1.5*K;
                Ky = 1.5*K;
                Kz = 1.5*K;
            }else
            {
                Kx = 1.5*K;
                Ky = 1.5*K;
                Kz = 1.5*K;
            }
            K_d_diag << Kx,Ky,Kz,30,30,30;
            //std::cout << "Starting impedance controller"<< "\n";

            if (log_file.is_open())
            {

                //std::cout << "Currently recording "<< recording_time <<"\n";
                log_file
                << recording_time << "," << p(0) << "," << p(1) << ","
                << p(2) << "," << p(3) << "," << p(4) << ","
                << p(5) << ",";
                log_file
                << dp(0) << "," << dp(1) << ","
                << dp(2) << "," << dp(3) << "," << dp(4) << ","
                << dp(5) << ",";
                log_file
                << f_ext(0) << "," << f_ext(1) << ","
                << f_ext(2) << "," << f_ext(3) << "," << f_ext(4) << ","
                << f_ext(5) << ",";
                log_file
                << p_d(0) << "," << p_d(1) << ","
                << p_d(2) << "," << p_d(3) << "," << p_d(4) << ","
                << p_d(5) << ",";
                log_file
                << dp_d(0) << "," << dp_d(1) << ","
                << dp_d(2) << "," << dp_d(3) << "," << dp_d(4) << ","
                << dp_d(5) << ",";


                log_file
                << horizon << ",";
                log_file
               << Kx << "," << Ky << ","
               << Kz << "," ;

                log_file
                << sensor << ",";

                log_file
                << current(0) << "," << current(1) << ","
                << current(2) << "," << current(3) << "," << current(4) << ","
                << current(5) << ","<< current(6) << ",";

                log_file << "\n";
            }
                // IF GOAL REACHED, STOP
                if (goal_reached == 1)
                {
                    K_d_diag << 800,800,800,40,40,40;
                    p_d[0] = p[0];
                    p_d[1] = p[1];
                    p_d[2] = p[2];
                    dp_d[0] = 0;
                    dp_d[1] = 0;
                    dp_d[2] = 0;
                }


            return std::tuple_cat(impedance_controller(q, dq, ddq, T_B7, p_d, dp_d,ddp_d, K_d_diag, K_n_diag,c_f, time_period,
                actuator_torque, horizon, std_p, sensor), std::make_tuple(goal_reached), std::make_tuple(0), std::make_tuple(0));
            } catch (const std::exception& e) {
                std::cout << "updateAndPredict threw: " << e.what() << std::endl;
            }



        }


    }else
    {


        VectorXd p(6);
        VectorXd std_p(3);
        tie(p, T_B7) = forward(q);
        vector<double> p2(p.data(), p.data() + p.size());
        MatrixXd jacobian(6,7);
        jacobian = jaco_m(q);
        VectorXd dp(6);
        dp = jacobian * dq;
        int initial_position_reached = 0;
        double threshold = 0.05;
        VectorXd f_ext(6);
        MatrixXd pinv_jacobian_T(6, 7);
        pinv_jacobian_T = (jacobian).completeOrthogonalDecomposition().pseudoInverse().transpose();
        f_ext = pinv_jacobian_T*actuator_torque;

        if (baseline_counter < 10)
        {
            baseline_Fz = baseline_Fz+f_ext[2]/10;
            baseline_counter++;
        } else
        {
            //std::cout << baseline_Fz << "\n";
        }

        if ((abs(p(0)-0.50)<threshold)&&(abs(p(1)-0.16)<threshold)&&(abs(p(2)-0.44)<threshold))
        {
            initial_position_reached = 1;
        }else
        {
            initial_position_reached = 0;
        }

        //std::cout << "Current position = " << p[0] << " " << p[1] << " " << p[2] << "\n";

        if (flag_reset == 1)
        {

                predictor.reset(override_duration);
                kalman_predictor = MultiDofUkfPredictor(3,initial_positions);
                //file_change_flag = 1;
                //file_num = file_num + 1;
                recording_time = 0;
                std::cout << "Changing file to p"<< file_num <<"\n";
                log_file.close();
                //filename = "C:/Imperial/Final project/predictive_controller/r" + std::to_string(file_num) + ".csv";
                //log_file.open("C:/Imperial/Final project/admittance control/m.csv");
                log_file.open(filename);

        }
        return std::tuple_cat(admittance_controller(q, dq, ddq, T_B7, K_d_diag_admittance, D_d_diag, I_d_diag,
                                                                actuator_torque, c_f, time_period, K_n_diag), std::make_tuple(0), std::make_tuple(1),
                                                                std::make_tuple(initial_position_reached));
    }

    /*
    if (GetAsyncKeyState(VK_SPACE) & 0x8000)
    {
        //std::cout << "Space pressed" << "\n";
        //predict_and_update
        VectorXd p(6);
        VectorXd std_p(3);
        tie(p, T_B7) = forward(q);
        vector<double> p2(p.data(), p.data() + p.size());
        MatrixXd jacobian(6,7);
        jacobian = jaco_m(q);
        VectorXd dp(6);
        dp = jacobian * dq;
        double total_velocity = sqrt(dp(0)*dp(0)+dp(1)*dp(1)+dp(2)*dp(2));

        double horizon = max(0.030, min(0.100, 0.01/(total_velocity+1e-5)));
        //std::cout << "Starting prediction" << "\n";
        promp_rt::PredictionResult result =
            predictor.update_and_predict(recording_time, p2, mode, record_dt_s, target_pos, horizon);
        //std::cout << "Finished prediction" << "\n";
        //USE RESULT TO CREATE p_d, dp_d, ddp_d
        //std::cout << "Size of mean_traj: " << result.mean_traj.size()<< "\n";

        p_d[0] = result.mean_traj[0];
        //std::cout << p_d[0] << "\n";
        p_d[1] = result.mean_traj[2];
        //std::cout << p_d[1] << "\n";
        p_d[2] = result.mean_traj[4];
        //std::cout << p_d[2] << "\n";
        dp_d[0] = result.mean_traj[1];
        //std::cout << dp_d[0] << "\n";
        dp_d[1] = result.mean_traj[3];
        //std::cout << dp_d[1] << "\n";
        dp_d[2] = result.mean_traj[5];
        //std::cout << dp_d[2] << "\n";
        //std::cout << "Starting std_traj"<< "\n";
        //std::cout << "Size of std_traj: " << result.std_traj.size()<< "\n";
        std_p[0] = result.std_traj[0];
        std_p[1] = result.std_traj[1];
        std_p[2] = result.std_traj[2];

        file_change_flag = 0;
        recording_time = recording_time + time_period;
        double K = 300;
        double V_lim = 0.1;
        double K_max = 500;
        double K_min = 100;
        if (total_velocity < V_lim)
        {
            if (recording_time > 1)
            {
                K = K_max - (K_max-K_min)*(V_lim-total_velocity)/V_lim;
            }else
            {
                K = K_max;
            }
        }else
        {
            K = K_max;
        }
        K_d_diag << K,K,K,10,10,10;
        //std::cout << "Starting impedance controller"<< "\n";
        return std::tuple_cat(impedance_controller(q, dq, ddq, T_B7, p_d, dp_d,ddp_d, K_d_diag, K_n_diag,c_f, time_period,
            actuator_torque, horizon, std_p, sensor), std::make_tuple(1));
    }else
    {

        if (file_change_flag == 0)
        {
            if (GetAsyncKeyState('Z') & 0x8000)
            {
                predictor.reset(override_duration);

                file_change_flag = 1;
                file_num = file_num + 1;
                recording_time = 0;
                std::cout << "Changing file to p"<< file_num <<"\n";
                log_file.close();
                filename = "C:/Imperial/Final project/predictive_controller/r" + std::to_string(file_num) + ".csv";
                //log_file.open("C:/Imperial/Final project/admittance control/m.csv");
                log_file.open(filename);

                //time = 0;
            }
        }
        return std::tuple_cat(admittance_controller(q, dq, ddq, T_B7, K_d_diag_admittance, D_d_diag, I_d_diag,
                                                                actuator_torque, c_f, time_period), std::make_tuple(1));

    }
    */
    VectorXd u;
    VectorXd dp;
    VectorXd ddp;
    return make_tuple(u, dp, ddp, acc_factor_buffer, 1,0,0);
}