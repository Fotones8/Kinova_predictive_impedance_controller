/**
 * @file main.cpp
 * @brief CLI entry point for ProMP_1 variable-DoF real-time training/prediction.
 *
 * ── Subcommands ───────────────────────────────────────────────────────────────
 *
 *   promp_rt train <demo_folder> <model_dir> [options]
 *     Load *.csv demonstrations, train a ProMP, save model + metadata.
 *
 *   promp_rt predict <model_dir> [options]
 *     Stream observations from stdin, write predictions to stdout.
 *
 * ── Training CSV format ───────────────────────────────────────────────────────
 *   One file per demonstration, (n_dof + 1) columns:
 *     time_s,  pos0_rad,  pos1_rad,  ...,  pos(n_dof-1)_rad
 *   Optional header row auto-detected.
 *   Different files may have different sampling rates.
 *
 * ── Prediction I/O ────────────────────────────────────────────────────────────
 *   Input  (stdin):  One observation per line:
 *     <time_s>  <pos0_rad>  <pos1_rad>  ...  <pos(n_dof-1)_rad>
 *
 *   Output (stdout): One block per observation:
 *     PRED <n_points> <current_phase>
 *     <time_s> <pos0_mean> <pos0_std> <vel0_mean> <vel0_std>
 *              <pos1_mean> <pos1_std> <vel1_mean> <vel1_std>
 *              ...
 *
 *   stdout is flushed after every block for pipe-compatible real-time use.
 *   An empty input line or EOF terminates the prediction loop.
 *
 * ── Data-frequency note ───────────────────────────────────────────────────────
 *   Observations may arrive at any rate; timestamps drive all phase mapping
 *   and SG velocity estimation.  The SG window stored in metadata.txt should
 *   be matched to the expected sensor rate (set via --sg-window at training).
 */

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "training.hpp"
#include "prediction.hpp"

#include <fstream>
#include <chrono>
// ============================================================
// Argument helpers
// ============================================================



/**
 * @brief Return the value of a named CLI flag, or default_val if absent.
 *
 * @param args        Full argv vector.
 * @param flag        Flag name (e.g. "--basis").
 * @param default_val Returned when the flag is not present.
 * @return Value string following the flag, or default_val.
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

/**
 * @brief Handle the 'train' subcommand.
 *
 * ── Positional arguments (required) ─────────────────────────────────────────
 *   args[2]  demo_folder  Directory of *.csv demonstration files.
 *   args[3]  model_dir    Output directory for model.promp + metadata.txt.
 *
 * ── Optional flags ───────────────────────────────────────────────────────────
 *   --dof N
 *     Number of degrees of freedom (joints) in the input CSVs.
 *     Each CSV must have (N+1) columns: time_s, pos0, …, pos(N-1).
 *     Default: 3.
 *     Tuning: set to the actual number of joints.  Mismatch causes a parse
 *     error at training time.
 *
 *   --basis N
 *     Number of RBF basis functions.
 *     Default: 10.
 *     Increase for complex multi-phase trajectories (up to ~20); decrease if
 *     you have very few demonstrations (<5).
 *
 *   --std S
 *     RBF standard deviation; ≤ 0 = auto-tuned (recommended).
 *     Default: -1.
 *
 *   --steps N
 *     Resampling grid size (common phase resolution for all demos).
 *     Default: 100.  Increase to 200 for fast transient dynamics.
 *
 *   --max-duration D
 *     Discard demo data beyond D seconds.  Default: -1 (keep all).
 *
 *   --sg-window W
 *     SG filter window length (odd integer ≥ sg-poly+2).
 *     Default: 9 (suitable for 100 Hz data, ~90 ms smoothing).
 *     Frequency scaling: W ≈ data_rate_Hz × desired_smoothing_s (odd).
 *       100 Hz, 90 ms → 9;   500 Hz, 90 ms → 45;   50 Hz, 90 ms → 5.
 *
 *   --sg-poly P
 *     SG polynomial order.  Default: 4.  Must be < (sg-window − 1).
 *
 *   --obs-pos-noise V
 *     Position observation noise variance for inference (rad²).
 *     Default: 1e-4.  Match to sensor RMS noise squared.
 *
 *   --obs-vel-noise V
 *     Velocity observation noise variance for inference ((rad/s)²).
 *     Default: 1e-2.  Typically 10–1000× obs-pos-noise.
 *
 *   --goal-pos-noise V
 *     Goal position noise variance (PAST_TRAJ_TARGET mode) (rad²).
 *     Default: 1e-6.  Decrease for a harder stop constraint.
 *
 *   --goal-vel-noise V
 *     Goal velocity noise variance ((rad/s)²).
 *     Default: 1e-2.  Increase to allow a more gradual approach.
 *
 * @param args  Full argv vector.
 * @return 0 on success, 1 on error.
 */
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
    cfg.num_basis_functions = std::stoi(get_arg(args, "--basis",          "10"));
    cfg.std_bf              = std::stod(get_arg(args, "--std",            "-1"));
    cfg.n_resample_steps    = std::stoi(get_arg(args, "--steps",          "100"));
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

/**
 * @brief Format and print one PredictionResult block to stdout.
 *
 * Output format:
 *   PRED <n_points> <current_phase>
 *   <time_s> <pos0_mean> <pos0_std> <vel0_mean> <vel0_std>
 *            <pos1_mean> <pos1_std> <vel1_mean> <vel1_std>  ...
 *
 * "PRED 0 <phase>" is written when result.valid == false.
 * stdout is flushed after each block for pipe-compatible real-time output.
 *
 * @param r  PredictionResult to print.
 */
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

/**
 * @brief Handle the 'predict' subcommand (real-time streaming loop).
 *
 * n_dof is determined automatically from the loaded model metadata.
 *
 * ── Positional arguments (required) ─────────────────────────────────────────
 *   args[2]  model_dir   Directory produced by 'train'.
 *
 * ── Optional flags ───────────────────────────────────────────────────────────
 *   --mode past_traj|past_traj_target
 *     Conditioning mode.  Default: past_traj.
 *
 *     past_traj          Use accumulated past observations only.
 *     past_traj_target   Also pin the final predicted position to --target.
 *                        Target velocity assumed 0 rad/s (natural stop).
 *
 *   --target v0 v1 … v(n_dof-1)
 *     Final target joint positions in radians.
 *     Exactly n_dof values must follow this flag.
 *     Used only with --mode past_traj_target.
 *
 *   --duration D
 *     Override the expected motion duration (s).  Affects the time-to-phase
 *     mapping.  Default: mean training demo duration from metadata.
 *     Set this if the current motion is known to be faster or slower than
 *     the training demonstrations.
 *
 *   --record-dt DT
 *     Time step (s) between predicted samples; should match the data's
 *     recording period so the prediction grid follows the raw time series.
 *     Default: -1 (auto-derive from the observed timestamps).  The number of
 *     predicted steps is the receding horizon (remaining motion time) divided
 *     by this step.
 *
 *   --max-obs N
 *     Cap the number of past observations used per predict() tick (0 = all).
 *     Default: 0.
 *     Tuning: set to 20–50 for 100–200 Hz sensors to bound per-tick latency.
 *
 * ── stdin format ─────────────────────────────────────────────────────────────
 *   One observation per line (exactly n_dof + 1 numbers):
 *     <time_s>  <pos0_rad>  <pos1_rad>  ...  <pos(n_dof-1)_rad>
 *
 *   An empty line or EOF terminates the loop.
 *
 * ── stdout format ────────────────────────────────────────────────────────────
 *   PRED <n_points> <current_phase>
 *   <time_s> <pos0_mean> <pos0_std> <vel0_mean> <vel0_std>
 *            <pos1_mean> <pos1_std> ...
 *   (one block per sensor tick; flushed immediately)
 *
 * ── Data-frequency note ───────────────────────────────────────────────────────
 *   Observations may arrive at any rate.  The time stamp in each line drives
 *   both the phase mapping and the SG velocity estimation.
 *
 * @param args  Full argv vector.
 * @return 0 on success, 1 on error.
 */
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

// ============================================================
// Entry point
// ============================================================

/**
 * @brief Programme entry point.
 *
 * @param argc  Argument count.
 * @param argv  Argument vector:
 *              argv[0] = programme name
 *              argv[1] = "train" | "predict"
 *              argv[2..] = subcommand-specific arguments
 *
 * @return 0 on success, 1 on error.
 */
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
