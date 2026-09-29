//-------------------------------------------------
// Header file for Controller.cpp
//-------------------------------------------------

#ifndef EXPERIMENTALCODES_CONTROLLER_H
#define EXPERIMENTALCODES_CONTROLLER_H

#include "Controller.h"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <tuple>

#include "training.hpp"
#include "prediction.hpp"

using namespace Eigen;
//using namespace std;

//----------------------------------------------
// Functions used in the main codes
//----------------------------------------------
// impedance_controller: compute the controlled system input, u
// ini_controller: initialize the controller
//----------------------------------------------
namespace Controller {
    std::tuple<VectorXd, VectorXd, VectorXd, VectorXd> impedance_controller(VectorXd& q, VectorXd& dq, VectorXd& ddq, MatrixXd& T_B7,
                                                             VectorXd& p_d, VectorXd& dp_d, VectorXd& ddp_d,
                                                             VectorXd& K_d_diag, VectorXd& K_n_diag,
                                                             int c_f, double& time_period,
                                                             VectorXd& actuator_torque, double horizon=0.03, VectorXd std_p = VectorXd::Zero(3),
                                                             double sensor = 0.0);
    
    std::tuple<VectorXd, VectorXd, VectorXd, VectorXd> admittance_controller(VectorXd& q, VectorXd& dq, VectorXd& ddq, MatrixXd& T_B7,
                                                           VectorXd& K_d_diag, VectorXd& D_d_diag, VectorXd& I_d_diag,
                                                           VectorXd& actuator_torque,
                                                           int c_f, double& time_period,
                                                           VectorXd& K_n_diag);

    std::tuple<VectorXd, VectorXd, VectorXd, VectorXd> position_admittance_controller(VectorXd& q, VectorXd& dq, VectorXd& ddq, MatrixXd& T_B7,
                                                           VectorXd& K_d_diag, VectorXd& D_d_diag, VectorXd& I_d_diag,
                                                           VectorXd& actuator_torque,
                                                           int c_f, double& time_period);

    std::tuple<VectorXd, VectorXd, VectorXd, VectorXd, int, int, int> predictive_impedance_controller(VectorXd& q, VectorXd& dq, VectorXd& ddq, MatrixXd& T_B7,
                                                             VectorXd& p_d, VectorXd& dp_d, VectorXd& ddp_d,
                                                             VectorXd& K_d_diag, VectorXd& K_n_diag,
                                                             int c_f, double& time_period,
                                                             double override_duration, promp_rt::ConditioningMode mode,
                                                             double record_dt_s, std::vector<double> target_pos,
                                                             VectorXd& K_d_diag_admittance,
                                                             VectorXd& D_d_diag, VectorXd& I_d_diag,
                                                             VectorXd& actuator_torque, VectorXd& current,
                                                             double sensor = 0.0,
                                                             int flag_mode = 0, int flag_reset = 0, int flag_goal = 0, int flag_running = 0,
                                                             std::string filename = "result");


    void ini_controller(Vector3d& pos, MatrixXd& T_B7, std::string model_dir, int max_obs);
}


#endif //EXPERIMENTALCODES_CONTROLLER_H
