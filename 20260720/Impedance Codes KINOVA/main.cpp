//---------------------------------------------------------
// The main function for impedance controller
//---------------------------------------------------------
// Description: The main function for the impedance controller used in the PHRI project
// Copyright: Yihan Liu 2024
//---------------------------------------------------------
#define _HAS_STD_BYTE 0

#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <fstream>
#include <thread>
#include <iomanip>
#include <atomic>
#include <chrono>
#include <cmath>
#include <mutex>


#include <KDetailedException.h>

#include <BaseClientRpc.h>
#include <BaseCyclicClientRpc.h>
#include <ActuatorConfigClientRpc.h>
#include <SessionClientRpc.h>
#include <SessionManager.h>

#include <RouterClient.h>
#include <TransportClientTcp.h>
#include <TransportClientUdp.h>

#include <google/protobuf/util/json_util.h>


#if defined(_MSC_VER)
#include <Windows.h>
#else
#include <unistd.h>
#endif
#include <time.h>

#include <Eigen/Dense>
#include <Jacobian.h>
#include <cmath>
#include <Fwd_kinematics.h>
#include <Dynamics.h>
#include <Controller.h>
#include <Filter.h>
#include <conio.h>

#include "training.hpp"
#include "prediction.hpp"

namespace k_api = Kinova::Api;
using namespace Jacobian;
using namespace Fwd_kinematics;
using namespace Dynamics;
using namespace Controller;
using namespace Filter;

#define IP_ADDRESS "192.168.1.10"
#define PORT 10000
#define PORT_REAL_TIME 10001
#define ACTUATOR_COUNT 7
#define CONTROL_FREQUENCY 450

// Maximum allowed waiting time during actions
constexpr auto TIMEOUT_PROMISE_DURATION = chrono::seconds{20};

// Create an event listener that will set the promise action event to the exit value
// Will set promise to either END or ABORT
// Use finish_promise.get_future.get() to wait and get the value
function<void(k_api::Base::ActionNotification)>
create_event_listener_by_promise(promise<k_api::Base::ActionEvent>& finish_promise)
{
    return [&finish_promise] (k_api::Base::ActionNotification notification)
    {
        const auto action_event = notification.action_event();
        switch(action_event)
        {
            case k_api::Base::ActionEvent::ACTION_END:
            case k_api::Base::ActionEvent::ACTION_ABORT:
                finish_promise.set_value(action_event);
                break;
            default:
                break;
        }
    };
}

//-----------------------------------------------------------
// Function of impedance control
//-----------------------------------------------------------
// 3 inputs:
// base; base_cyclic; actuator_config = paramters for robotic communication
// p_d, dp_d, ddp_d = desired position, velocity and acceleration
// K_d_diag = diagonal values for desired stiffness matrix

// 1 output:
// status of controller: True or False
//------------------------------------------------------------
bool impedance_control(k_api::Base::BaseClient* base, k_api::BaseCyclic::BaseCyclicClient* base_cyclic,
                       k_api::ActuatorConfig::ActuatorConfigClient* actuator_config,
                       VectorXd& p_d, VectorXd& dp_d, VectorXd& ddp_d, VectorXd& K_d_diag,
                       std::string model_dir, promp_rt::ConditioningMode mode, double override_duration,
                       double record_dt_s, int max_obs, std::vector<double> target_pos) {
    std::cout << "Starting impedance control" << "\n";
    bool return_status = true;

    // Clearing faults
    try {
        base->ClearFaults();
    }
    catch (...) {
        cout << "Unable to clear robot faults" << endl;
        return false;
    }

    k_api::BaseCyclic::Feedback base_feedback;
    k_api::BaseCyclic::Command base_command;

    vector<float> commands;

    auto servoing_mode = k_api::Base::ServoingModeInformation();

    int timer_count = 0;

    // KINOVA feedback (joint space variables)
    VectorXd q(7), dq(7), last_dq(7), ddq(7), torque(7), current(7);

    // Task space variables
    VectorXd u(7), p(6), dp(6), ddp(6);
    MatrixXd T_B7(4,4);  // Rotation matrix

    //EXTERNAL FORCE (TASK SPACE)
    VectorXd f_ext(6);


    // Define the D_n and K_n for nullspace
    VectorXd K_n_diag(7);
    K_n_diag << 4, 4, 4, 2, 2, 2, 2;

    // Time for one control iterative
    double dt = 1.0 / CONTROL_FREQUENCY;
    const double iteration_time = (1.0 / CONTROL_FREQUENCY) * 1000;

    // Buffer to save the last velocity
    last_dq << 0,0,0,0,0,0,0;

    cout << "Initializing the arm for torque control ^^!" << endl;
    try
    {
        // Set the base in low-level servoing mode
        servoing_mode.set_servoing_mode(k_api::Base::ServoingMode::LOW_LEVEL_SERVOING);
        base->SetServoingMode(servoing_mode);
        base_feedback = base_cyclic->RefreshFeedback();

        // Initialize each actuator to their current position
        for (int i = 0; i < ACTUATOR_COUNT; i++)
        {
            commands.push_back(base_feedback.actuators(i).position());

            // Save the current actuator position, to avoid a following error
            base_command.add_actuators()->set_position(base_feedback.actuators(i).position());
        }

        // Send a first frame
        base_feedback = base_cyclic->Refresh(base_command);

        // Set actuatorS in torque mode now that the command is equal to measure
        auto control_mode_message = k_api::ActuatorConfig::ControlModeInformation();
        //control_mode_message.set_control_mode(k_api::ActuatorConfig::ControlMode::POSITION);
        control_mode_message.set_control_mode(k_api::ActuatorConfig::ControlMode::TORQUE);
        for (int id = 1; id < ACTUATOR_COUNT+1; id++)
        {
            actuator_config->SetControlMode(control_mode_message, id);
        }
        auto information_servo_mode = base->GetServoingMode();
        std::cout << "Servoing mode = "
          << information_servo_mode.servoing_mode()
          << std::endl;

        this_thread::sleep_for(chrono::milliseconds(40));

        // Clock to record the time period between two control/measure loop
        auto start_measure = chrono::high_resolution_clock::now();
        auto start_control = chrono::high_resolution_clock::now();

        // Control loop
        cout << "Starting control loop" << endl;
        while (1) // This is an example code
        {

            // KINOVA Feedback: Obtaining gen3 ACTUAL joint positions, velocities, torques & current
            for (int i = 0; i < ACTUATOR_COUNT; i++)
            {
                q[i] = (M_PI/180)*base_feedback.actuators(i).position();
                dq[i] = (M_PI/180)*base_feedback.actuators(i).velocity();
                torque[i] = base_feedback.actuators(i).torque();
                current[i] = base_feedback.actuators(i).current_motor();
            }
            //cout << "Feedback obtained" << endl;
            // GETTING END EFFECTOR EXTERNAL FORCE
            /*
            base_feedback = base_cyclic->RefreshFeedback();
            f_ext << base_feedback.base().tool_external_wrench_force_x(),
                 base_feedback.base().tool_external_wrench_force_y(),
                 base_feedback.base().tool_external_wrench_force_z(),
                 base_feedback.base().tool_external_wrench_torque_x(),
                 base_feedback.base().tool_external_wrench_torque_y(),
                 base_feedback.base().tool_external_wrench_torque_z();
            
            std::cout << "[INFO] External force " << f_ext << "\n";
            */

            // Apply the forward kinematics
            tie(p, T_B7) = forward(q);

            // initilize the controller and filter
            if(timer_count == 0){
                Vector<double, 3> pos = p.head(3);
                ini_controller(pos, T_B7, model_dir, max_obs);
                ini_butterworth();
            }

            // Compute the joint accelerations
            auto end_measure = chrono::high_resolution_clock::now();
            chrono::duration<double> measure_dur = end_measure - start_measure;
            dt = measure_dur.count();
            if (timer_count == 0){
                dt = 1.0 / CONTROL_FREQUENCY;
            }
            ddq = (dq - last_dq) / dt;
            start_measure = end_measure;
            last_dq = dq;
            //cout << "Computed acceleration" << endl;


            VectorXd acc_factor(6);
            /*
            // Impedance controller
            tie(u, dp, ddp, acc_factor) = impedance_controller(q, dq, ddq, T_B7, p_d, dp_d,
                                                               ddp_d, K_d_diag, K_n_diag,CONTROL_FREQUENCY, dt);
            */
            //cout << "Producing admittance matrix" << endl;
            VectorXd K_d_diag_admittance(6), D_d_diag(6), I_d_diag(6);
            D_d_diag << 20.0, 20.0, 20.0, 11.0,  11.0,  11.0;
            //cout << "D_d_diag produced" << endl;
            I_d_diag << 2,  2,  2,  3,  3,  3;
            //cout << "I_d_diag produced" << endl;
            //I_d_diag << 0.5,  0.5,  0.5,  0.3,  0.3,  0.3;
            K_d_diag_admittance << 0.0,  0.0,  0.0,  30.0,  30.0,  30.0;
            //cout << "K_d_diag_admittance produced" << endl;
            /*
            tie(u, dp, ddp, acc_factor) = admittance_controller(q, dq, ddq, T_B7, K_d_diag_admittance, D_d_diag, I_d_diag,
                                                                torque, CONTROL_FREQUENCY, dt);
            */
            //cout << "Starting predictive impedance control" << endl;
            //std::cout << "Starting predictive impedance control" << "\n";

            tie(u, dp, ddp, acc_factor) = predictive_impedance_controller(q,dq,ddq,T_B7,p_d, dp_d, ddp_d,
                K_d_diag,K_n_diag,CONTROL_FREQUENCY,dt,override_duration,mode,
                record_dt_s,target_pos,K_d_diag_admittance,D_d_diag,I_d_diag,torque);
            /*
            tie(u, dp, ddp, acc_factor) = position_admittance_controller(q, dq, ddq, T_B7, K_d_diag_admittance, D_d_diag, I_d_diag,
                                                                torque, CONTROL_FREQUENCY, dt);
            */

            // Set the control frequency
            auto end_control = chrono::high_resolution_clock::now();
            auto run_time = chrono::duration<double, milli>(end_control - start_control).count();
            int diff = 0;
            if (run_time < iteration_time) { // if real control loop time < desired iteration time
                auto start_delay = chrono::high_resolution_clock::now();
                diff = (iteration_time - run_time) * 1000;
                chrono::microseconds delay(diff);
                while (chrono::high_resolution_clock::now() - start_delay < delay); // delay the time difference
            }
            auto last_control_time = start_control;
            start_control = chrono::high_resolution_clock::now();
            auto control_duration = chrono::duration<double, milli>(start_control - last_control_time).count();

            for (int i = 0; i < ACTUATOR_COUNT; i++)
            {

                // -- Position
                base_command.mutable_actuators(i)->set_position(base_feedback.actuators(i).position());
                // -- Torque
                base_command.mutable_actuators(i)->set_torque_joint(u[i]);
                /*
                // -- Position
                base_command.mutable_actuators(i)->set_position(u[i]*180/M_PI);
                // -- Torque
                //base_command.mutable_actuators(i)->set_torque_joint(base_feedback.actuators(i).torque());

                std::cout
                    << "Joint " << i
                    << " feedback=" << base_feedback.actuators(i).position()
                    << " u(rad)=" << u[i]
                    << " u(deg)=" << u[i] * 180.0 / M_PI
                    << std::endl;
                */
            }


            // Incrementing identifier ensures actuators can reject out of time frames
            base_command.set_frame_id(base_command.frame_id() + 1);
            if (base_command.frame_id() > 65535)
                base_command.set_frame_id(0);

            for (int idx = 0; idx < ACTUATOR_COUNT; idx++)
            {
                base_command.mutable_actuators(idx)->set_command_id(base_command.frame_id());
            }
            /*
            auto information_servo_mode = base->GetServoingMode();
            std::cout << "Servoing mode = "
                << information_servo_mode.servoing_mode()
                << std::endl;
                */
            try
            {


                base_feedback = base_cyclic->Refresh(base_command, 0);
            }
            catch (k_api::KDetailedException& ex)
            {
                cout << "Kortex exception: " << ex.what() << endl;

                cout << "Error sub-code: " << k_api::SubErrorCodes_Name(k_api::SubErrorCodes((ex.getErrorInfo().getError().error_sub_code()))) << endl;
            }
            catch (runtime_error& ex2)
            {
                cout << "runtime error: " << ex2.what() << endl;
            }
            catch(...)
            {
                cout << "Unknown error." << endl;
            }

            timer_count++;
        }

        cout << "Torque control ^^ completed" << endl;

        // Set actuators back in position
        control_mode_message.set_control_mode(k_api::ActuatorConfig::ControlMode::POSITION);
        for (int id = 1; id < ACTUATOR_COUNT+1; id++)
        {
            actuator_config->SetControlMode(control_mode_message, id);
        }

        cout << "Torque control ^^ clean exit" << endl;

    }
    catch (k_api::KDetailedException& ex)
    {
        cout << "API error: " << ex.what() << endl;
        return_status = false;
    }
    catch (runtime_error& ex2)
    {
        cout << "Error: " << ex2.what() << endl;
        return_status = false;
    }

    // Set the servoing mode back to Single Level
    servoing_mode.set_servoing_mode(k_api::Base::ServoingMode::SINGLE_LEVEL_SERVOING);
    base->SetServoingMode(servoing_mode);

    // Wait for a bit
    this_thread::sleep_for(chrono::milliseconds(2000));

    return return_status;
}


//------------------------------------------
// Function of high-level movement
//-----------------------------------------
// 2 inputs:
// base: communication variables
// q_d: desired joint angular position
//-----------------------------------------
void Move_high_level(k_api::Base::BaseClient* base, VectorXd q_d)
{

    auto action = k_api::Base::Action();

    auto reach_joint_angles = action.mutable_reach_joint_angles();
    auto joint_angles = reach_joint_angles->mutable_joint_angles();

    auto actuator_count = base->GetActuatorCount();

    // Arm straight up
    for (size_t i = 0; i < actuator_count.count(); ++i)
    {
        auto joint_angle = joint_angles->add_joint_angles();
        joint_angle->set_joint_identifier(i);
        joint_angle->set_value(q_d[i]);
    }

    promise<k_api::Base::ActionEvent> finish_promise;
    auto finish_future = finish_promise.get_future();
    auto promise_notification_handle = base->OnNotificationActionTopic(
            create_event_listener_by_promise(finish_promise),
            k_api::Common::NotificationOptions()
    );

    base->ExecuteAction(action);

    const auto status = finish_future.wait_for(TIMEOUT_PROMISE_DURATION);
    base->Unsubscribe(promise_notification_handle);

    if(status != future_status::ready)
    {
        cout << "Timeout on action notification wait" << endl;
    }
    const auto promise_event = finish_future.get();
}

/**
 * @brief Return the value of a named CLI flag, or default_val if absent.
 */
static std::string get_arg(const std::vector<std::string>& args,
                            const std::string& flag,
                            const std::string& default_val = "")
{
    for (size_t i = 0; i + 1 < args.size(); ++i)
        if (args[i] == flag) return args[i + 1];
    return default_val;
}

/**
 * @brief Return true if @p flag appears anywhere in @p args.
 */
static bool has_flag(const std::vector<std::string>& args,
                     const std::string& flag)
{
    return std::find(args.begin(), args.end(), flag) != args.end();
}

// ============================================================
// Training mode
// ============================================================
static int run_train(const std::vector<std::string>& args)
{
    if (args.size() < 4 || has_flag(args, "--help")) {
        std::cerr <<
            "Usage: promp_rt train <demo_folder> <model_dir> [options]\n"
            "\n"
            "Demo CSV: time_s, pos0_rad, pos1_rad, ..., pos(n_dof-1)_rad\n"
            "\n"
            "Options:\n"
            "  --dof N            Number of joints (default 3)\n"
            "  --basis N          RBF basis functions (default 10)\n"
            "  --std S            Basis std-dev; <=0 = auto (default -1)\n"
            "  --steps N          Resampling grid size (default 100)\n"
            "  --max-duration D   Clip demos at D seconds (default: all)\n"
            "  --sg-window W      SG window length, odd (default 9)\n"
            "  --sg-poly P        SG polynomial order (default 4)\n"
            "  --obs-pos-noise V  Position observation noise var (default 1e-4)\n"
            "  --obs-vel-noise V  Velocity observation noise var (default 1e-2)\n"
            "  --goal-pos-noise V Goal position noise var        (default 1e-6)\n"
            "  --goal-vel-noise V Goal velocity noise var        (default 1e-2)\n";
        return 1;
    }

    promp_rt::TrainingConfig cfg;
    cfg.n_dof               = std::stoi(get_arg(args, "--dof",            "3"));
    cfg.num_basis_functions = std::stoi(get_arg(args, "--basis",          "15"));
    cfg.std_bf              = std::stod(get_arg(args, "--std",            "-1"));
    cfg.n_resample_steps    = std::stoi(get_arg(args, "--steps",          "300"));
    cfg.max_demo_duration   = std::stod(get_arg(args, "--max-duration",   "-1"));
    cfg.sg_window           = std::stoi(get_arg(args, "--sg-window",      "9"));
    cfg.sg_poly_order       = std::stoi(get_arg(args, "--sg-poly",        "4"));
    cfg.obs_pos_noise_var   = std::stod(get_arg(args, "--obs-pos-noise",  "0.0001"));
    cfg.obs_vel_noise_var   = std::stod(get_arg(args, "--obs-vel-noise",  "0.01"));
    cfg.goal_pos_noise_var  = std::stod(get_arg(args, "--goal-pos-noise", "0.000001"));
    cfg.goal_vel_noise_var  = std::stod(get_arg(args, "--goal-vel-noise", "0.01"));

    return promp_rt::train_from_folder(args[2], args[3], cfg) ? 0 : 1;
}

// ============================================================
// Prediction mode
// ============================================================
static void print_result(const promp_rt::PredictionResult& r)
{
    if (!r.valid) {
        std::cout << "PRED 0 " << r.current_phase << "\n";
        std::cout.flush();
        return;
    }

    std::cout << "PRED " << r.future_times_s.size()
              << " "     << r.current_phase << "\n";
    /*
    const int n     = static_cast<int>(r.future_times_s.size());
    const int ndims = r.n_dims;
    for (int i = 0; i < n; ++i) {
        std::cout << r.future_times_s[static_cast<size_t>(i)];
        for (int dim = 0; dim < ndims; ++dim) {
            // Interleave mean and std: mean[dim] std[dim]
            std::cout << " " << r.mean_traj[static_cast<size_t>(i * ndims + dim)];
            //          << " " << r.std_traj [static_cast<size_t>(i * ndims + dim)];
        }
        std::cout << "\n";
    }
    */
    //std::cout.flush();
}


static int run_predict(const std::vector<std::string>& args)
{
    if (args.size() < 3 || has_flag(args, "--help")) {
        std::cerr <<
            "Usage: promp_rt predict <model_dir> [options]\n"
            "\n"
            "Options:\n"
            "  --mode past_traj|past_traj_target   (default past_traj)\n"
            "  --target v0 v1 … v(n_dof-1)         target joint positions (rad)\n"
            "  --duration D                         override expected duration (s)\n"
            "  --record-dt DT                       prediction step in s; <=0 = auto from timestamps\n"
            "  --max-obs N                          obs per tick cap (default 0=all)\n"
            "  --input PATH                         read observations from file (default: stdin)\n"
            "\n"
            "stdin:  <time_s> <pos0> <pos1> ...  (one per line)\n"
            "stdout: PRED <n> <phase>  then n lines of predictions\n";
        return 1;
    }
    // Disable C++/C stdio sync for speed
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);


    const std::string model_dir    = args[2];
    const std::string mode_str     = get_arg(args, "--mode", "past_traj");
    const double override_duration = std::stod(get_arg(args, "--duration", "-1"));
    const double record_dt_s       = std::stod(get_arg(args, "--record-dt", "-1"));
    const int    max_obs           = std::stoi(get_arg(args, "--max-obs", "0"));

    // Empty default → read observations from stdin (see loop below).  A
    // non-empty default (e.g. "0") would make the code try to open a file by
    // that name and the documented stdin mode would be unreachable.
    const std::string input_path   = get_arg(args, "--input", "");

    // Construct predictor first so we can query n_dof for target parsing.
    promp_rt::OnlinePredictor predictor(model_dir, max_obs);
    const int n_dof = predictor.get_n_dof();

    // Parse --target  v0  v1  …  v(n_dof-1)
    // Must have exactly n_dof values following the flag.
    std::vector<double> target_pos(static_cast<size_t>(n_dof), 0.0);
    {
        auto it = std::find(args.begin(), args.end(), "--target");
        if (it != args.end()) {
            for (int d = 0; d < n_dof; ++d) {
                ++it;
                if (it == args.end()) {
                    std::cerr << "[predict] --target requires " << n_dof
                              << " values (n_dof=" << n_dof << ")\n";
                    return 1;
                }
                target_pos[static_cast<size_t>(d)] = std::stod(*it);
            }
        }
    }

    const promp_rt::ConditioningMode mode =
        (mode_str == "past_traj_target")
            ? promp_rt::ConditioningMode::PAST_TRAJ_TARGET
            : promp_rt::ConditioningMode::PAST_TRAJ;

    predictor.reset(override_duration);

    // ── Real-time streaming loop ──────────────────────────────────────────


    //while (std::getline(std::cin, line)) { // for manually input line
    std::ifstream input_file;
    std::istream* in = &std::cin;

    if (!input_path.empty()) {
        input_file.open(input_path);
        if (!input_file.is_open()) {
            std::cerr << "[predict] Could not open input file: " << input_path << "\n";
            return 1;
        }
        in = &input_file;
    }

    // Each stdin line produces exactly one PRED block on stdout.
    std::vector<double> pos_buf(static_cast<size_t>(n_dof));
    std::string line;

    while (std::getline(*in, line)) {

        if (line.empty()) break;

        std::istringstream ss(line);
        double time_s{};
        if (!(ss >> time_s)) {
            std::cerr << "[predict] Cannot parse line: " << line << "\n";
            continue;
        }
        bool ok = true;
        for (int d = 0; d < n_dof; ++d) {
            if (!(ss >> pos_buf[static_cast<size_t>(d)])) { ok = false; break; }
        }
        if (!ok) {
            std::cerr << "[predict] Expected " << n_dof
                      << " joint values on line: " << line << "\n";
            continue;
        }

        auto t_start = std::chrono::high_resolution_clock::now();
        promp_rt::PredictionResult result =
            predictor.update_and_predict(
                time_s, pos_buf, mode, record_dt_s, target_pos);
        auto t_end = std::chrono::high_resolution_clock::now();
        double us = std::chrono::duration<double, std::micro>(t_end - t_start).count();
        std::cout << result.n_obs_used << "," << us << "\n";
        //print_result(result);
    }

    return 0;
}



//----------------------------------------------------
// Main function of impedance control
//----------------------------------------------------

int main(int argc, char **argv)
{
    std::vector<std::string> args(argv, argv + argc);


    const std::string model_dir    = args[2];
    const std::string mode_str     = get_arg(args, "--mode", "past_traj");
    const double override_duration = std::stod(get_arg(args, "--duration", "-1"));
    const double record_dt_s       = std::stod(get_arg(args, "--record-dt", "-1"));
    const int    max_obs           = std::stoi(get_arg(args, "--max-obs", "0"));

    std::vector<double> target_pos(static_cast<size_t>(3), 0.0);
    {
        auto it = std::find(args.begin(), args.end(), "--target");
        if (it != args.end()) {
            for (int d = 0; d < 3; ++d) {
                ++it;
                if (it == args.end()) {
                    std::cerr << "[predict] --target requires " << 3
                              << " values (n_dof=" << 3 << ")\n";
                    return 1;
                }
                target_pos[static_cast<size_t>(d)] = std::stod(*it);
            }
        }
    }
    const promp_rt::ConditioningMode mode =
        (mode_str == "past_traj_target")
            ? promp_rt::ConditioningMode::PAST_TRAJ_TARGET
            : promp_rt::ConditioningMode::PAST_TRAJ;

    VectorXd p_d(6), dp_d(6), ddp_d(6), K_d_diag(6);
    // Define the desired position, velocity, acceleration and stiffness here.
     p_d << 0.442,-0.06,0.575,1.4,0,1.4;
    dp_d << 0,0,0,0,0,0;
     ddp_d << 0,0,0,0,0,0;
     K_d_diag << 500,500,500,10,10,10;

    // Experiment loop for different sets
    this_thread::sleep_for(chrono::milliseconds(1000));

    // Create API objects
    auto error_callback = [](k_api::KError err) { cout << "_________ callback error _________" << err.toString(); };

    cout << "Creating transport objects" << endl;
    auto transport = new k_api::TransportClientTcp();
    auto router = new k_api::RouterClient(transport, error_callback);
    transport->connect(IP_ADDRESS, PORT);

    cout << "Creating transport real time objects" << endl;
    auto transport_real_time = new k_api::TransportClientUdp();
    auto router_real_time = new k_api::RouterClient(transport_real_time, error_callback);
    transport_real_time->connect(IP_ADDRESS, PORT_REAL_TIME);

    // Set session data connection information
    auto create_session_info = k_api::Session::CreateSessionInfo();
    create_session_info.set_username("admin");
    create_session_info.set_password("admin");
    create_session_info.set_session_inactivity_timeout(60000);   // (milliseconds)
    create_session_info.set_connection_inactivity_timeout(2000); // (milliseconds)

    // Session manager service wrapper
    cout << "Creating sessions for communication" << endl;
    auto session_manager = new k_api::SessionManager(router);
    session_manager->CreateSession(create_session_info);
    auto session_manager_real_time = new k_api::SessionManager(router_real_time);
    session_manager_real_time->CreateSession(create_session_info);
    cout << "Sessions created" << endl;

    // Create services
    auto base = new k_api::Base::BaseClient(router);
    auto base_cyclic = new k_api::BaseCyclic::BaseCyclicClient(router_real_time);
    auto actuator_config = new k_api::ActuatorConfig::ActuatorConfigClient(router);
    cout << "Services created" << endl;

    // Move to the initial configuration
    VectorXd q_d(7);
    q_d << 0, 0.261887, 3.14162, 4.01422, 2.72361e-05, 0.959906, 1.57;
    q_d = q_d * (180 / M_PI);
    cout << "Moving to initial configuration" << endl;
    //Move_high_level(base, q_d);
    cout << "Moved to initial configuration" << endl;

    // Atomic flag to control the loop termination for drawScene
    atomic<bool> is_impedance_control_running(true);

        // Thread for impedance_control function
        thread impedance_thread([&]() {
            auto isOk = impedance_control(base, base_cyclic, actuator_config, p_d, dp_d, ddp_d, K_d_diag,
             model_dir,  mode,  override_duration, record_dt_s, max_obs,target_pos);
            if (!isOk) {
                cout << "There has been an unexpected error in impedance_control() function." << endl;
            }
            // Signal that the impedance control has finished
            is_impedance_control_running = false;
        });

        // Join the threads back to the main thread
        impedance_thread.join();

        // Close API session
        session_manager->CloseSession();
        session_manager_real_time->CloseSession();

        // Deactivate the router and cleanly disconnect from the transport object
        router->SetActivationStatus(false);
        transport->disconnect();
        router_real_time->SetActivationStatus(false);
        transport_real_time->disconnect();

        // Destroy the API
        delete base;
        delete base_cyclic;
        delete actuator_config;
        delete session_manager;
        delete session_manager_real_time;
        delete router;
        delete router_real_time;
        delete transport;
        delete transport_real_time;
}




// ============================================================
// Entry point
// ============================================================
/*
int main(int argc, char* argv[])
{
    std::vector<std::string> args(argv, argv + argc);

    if (args.size() < 2 || has_flag(args, "--help") || has_flag(args, "-h")) {
        std::cout <<
            "ProMP_1 real-time variable-DoF training and prediction tool\n"
            "Usage:\n"
            "  promp_rt train   <demo_folder> <model_dir> [--dof N] [options]\n"
            "  promp_rt predict <model_dir> --mode <mode> [options]\n"
            "Run 'promp_rt <subcommand> --help' for full option list.\n";
        return 0;
    }

    try {
        if (args[1] == "train")   return run_train(args);
        if (args[1] == "predict") return run_predict(args);

        std::cerr << "Unknown subcommand '" << args[1]
                  << "'.  Use 'train' or 'predict'.\n";
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "[error] " << e.what() << "\n";
        return 1;
    }
}
*/