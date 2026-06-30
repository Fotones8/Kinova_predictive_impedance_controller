/**
 * @file prediction.cpp
 * @brief ProMP_1 online real-time prediction for variable-DoF trajectories.
 *
 * ── Algorithm overview ────────────────────────────────────────────────────────
 * At each call to predict():
 *   1. A fresh ProMP copy is reconstructed from the stored (μ_w, Σ_w).
 *      The trained model is NEVER modified, so predict() is re-entrant over
 *      separate OnlinePredictor instances.
 *   2. Velocities for all joints are estimated from the full observation
 *      buffer using the Savitzky-Golay filter (parameters from metadata.txt).
 *   3. For each selected past observation a 2*n_dof–dimensional via-point
 *      [pos0, vel0, pos1, vel1, …] is formed and Bayesian conditioning is
 *      applied via promp::ProMP::condition_via_point().
 *   4. If PAST_TRAJ_TARGET: condition_goal() is applied with the target
 *      positions (tight noise) and zero target velocities (soft noise).
 *   5. The full trajectory is generated; only the future slice (from the
 *      current phase onward) is returned.
 *
 * ── Data-frequency compatibility ─────────────────────────────────────────────
 * Observations may arrive at any rate with non-uniform timestamps:
 *   • Phase mapping:  phase = time_s / expected_duration_s  (no rate assumed).
 *   • SG velocity:    uses real timestamps in the buffer.
 *   • The SG window should be matched to the expected sensor rate at training
 *     time; it is loaded from metadata.txt and applied here unchanged.
 *
 * ── Dimension layout (2*n_dof dims) ──────────────────────────────────────────
 *   dim 0: pos0 (rad)    dim 1: vel0 (rad/s)
 *   dim 2: pos1 (rad)    dim 3: vel1 (rad/s)
 *   ...
 *   dim 2k: pos_k        dim 2k+1: vel_k
 */

#include "prediction.hpp"
#include "sg_filter.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <promp/promp.hpp>
#include <promp/io/serializer.hpp>


//Save the results in a file
#include <fstream>
std::ofstream log_file;


namespace fs = std::filesystem;

namespace promp_rt {

// ============================================================
// Metadata
// ============================================================

namespace {

/**
 * @brief Runtime parameters loaded from metadata.txt.
 *
 * Every field maps 1-to-1 to a key written by write_metadata() in training.cpp.
 * Unrecognised keys are silently ignored to allow forward compatibility.
 */
struct Metadata {
    int    n_dof             = 3;      ///< Number of joints.
    int    n_steps           = 100;    ///< Phase resolution of the trained trajectory.
    double mean_duration_s   = 1.0;    ///< Mean demo duration (s); default for phase map.
    int    sg_window         = 9;      ///< SG filter window length for online velocity.
    int    sg_poly_order     = 4;      ///< SG polynomial order.
    double obs_pos_noise_var = 1e-4;   ///< Position observation noise variance (rad²).
    double obs_vel_noise_var = 1e-2;   ///< Velocity observation noise variance ((rad/s)²).
    double goal_pos_noise_var= 1e-6;   ///< Goal position noise variance (tight, rad²).
    double goal_vel_noise_var= 1e-2;   ///< Goal velocity noise variance (soft, (rad/s)²).
};

/**
 * @brief Parse metadata.txt produced by train_from_folder().
 *
 * @param path  Path to metadata.txt.
 * @return Metadata struct with all recognised fields populated.
 * @throws std::runtime_error if the file cannot be opened.
 */
static Metadata load_metadata(const fs::path& path)
{
    std::ifstream f(path);
    if (!f.is_open())
        throw std::runtime_error("Cannot read metadata: " + path.string());

    Metadata m;
    std::string line;
    while (std::getline(f, line)) {
        auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        const std::string key   = line.substr(0, eq);
        const std::string value = line.substr(eq + 1);
        try {
            if      (key == "n_dof")              m.n_dof             = std::stoi(value);
            else if (key == "n_steps")            m.n_steps           = std::stoi(value);
            else if (key == "mean_duration_s")    m.mean_duration_s   = std::stod(value);
            else if (key == "sg_window")          m.sg_window         = std::stoi(value);
            else if (key == "sg_poly_order")      m.sg_poly_order     = std::stoi(value);
            else if (key == "obs_pos_noise_var")  m.obs_pos_noise_var = std::stod(value);
            else if (key == "obs_vel_noise_var")  m.obs_vel_noise_var = std::stod(value);
            else if (key == "goal_pos_noise_var") m.goal_pos_noise_var= std::stod(value);
            else if (key == "goal_vel_noise_var") m.goal_vel_noise_var= std::stod(value);
        } catch (...) {}
    }
    return m;
}

} // anonymous namespace

// ============================================================
// Impl (pimpl — keeps Eigen and promp out of the public header)
// ============================================================

struct OnlinePredictor::Impl {
    promp::ProMP trained;   ///< Trained weight distribution; NEVER mutated after init.
    Metadata     meta;      ///< All parameters loaded from metadata.txt.
    int          max_obs;   ///< 0 = use all observations per predict() call.

    // ── Observation buffer ────────────────────────────────────────────────
    std::vector<double>              obs_times; ///< Elapsed time (s) per observation.
    std::vector<std::vector<double>> obs_pos;   ///< obs_pos[i][d] = joint d at observation i.

    double expected_duration_s;  ///< Active expected duration (may be overridden by reset()).

    /**
     * @brief Load model and metadata from disk; initialise the observation buffer.
     *
     * @param model_dir            Path written by train_from_folder().
     * @param max_obs_per_predict  Cap on observations per predict() call.
     *                             0 = unlimited.
     * @throws std::runtime_error  on I/O failure.
     */
    Impl(const std::string& model_dir, int max_obs_per_predict)
        : trained(promp::io::load_promp(
              (fs::path(model_dir) / "model.promp").string())),
          meta(load_metadata(fs::path(model_dir) / "metadata.txt")),
          max_obs(max_obs_per_predict),
          expected_duration_s(meta.mean_duration_s)
    {}

    /**
     * @brief Map an elapsed time to a trajectory phase step index.
     *
     * Computation:
     *   phase = clamp(time_s / expected_duration_s,  0,  1)
     *   step  = clamp(round(phase × (n_steps − 1)),  0,  n_steps−1)
     *
     * ── Data-frequency note ────────────────────────────────────────────────
     * This mapping is purely time-based and independent of sensor rate.
     * Accuracy depends on expected_duration_s matching the actual motion
     * duration.  Call reset(actual_duration_s) when this is known.
     *
     * @param time_s  Elapsed time since motion start (s).
     * @return Step index in [0, n_steps−1].
     */
    int time_to_step(double time_s) const
    {
        const double phase = std::max(0.0, std::min(1.0,
                                 time_s / expected_duration_s));
        const int step = static_cast<int>(
            std::round(phase * static_cast<double>(meta.n_steps - 1)));
        return std::max(0, std::min(meta.n_steps - 1, step));
    }

    /**
     * @brief Build the 2*n_dof × 2*n_dof block-diagonal observation noise
     *        covariance matrix Σ_obs for via-point conditioning.
     *
     * Structure (for each joint pair [pos_d, vel_d]):
     *   ⎡ σ²_pos   0      ⎤
     *   ⎣  0     σ²_vel   ⎦
     *
     * Position channels use obs_pos_noise_var; velocity channels use the
     * larger obs_vel_noise_var to account for SG filter estimation error.
     *
     * Tuning:
     *   obs_pos_noise_var – match to sensor RMS position noise (rad²).
     *   obs_vel_noise_var – typically 10–1000× obs_pos_noise_var.
     *
     * @return  2*n_dof × 2*n_dof diagonal covariance matrix.
     */
    Eigen::MatrixXd make_obs_cov() const
    {
        const int nd = meta.n_dof;
        Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(2 * nd, 2 * nd);
        for (int d = 0; d < nd; ++d) {
            cov(d * 2,     d * 2)     = meta.obs_pos_noise_var;
            cov(d * 2 + 1, d * 2 + 1) = meta.obs_vel_noise_var;
        }
        return cov;
    }

    /**
     * @brief Build the 2*n_dof × 2*n_dof block-diagonal goal noise covariance
     *        matrix for condition_goal().
     *
     * Position channels: tight (goal_pos_noise_var) – enforces the target.
     * Velocity channels: soft  (goal_vel_noise_var)  – allows natural stop.
     *
     * Tuning:
     *   goal_pos_noise_var – decrease for a hard position constraint (~1e-7).
     *   goal_vel_noise_var – increase if the motion need not stop exactly.
     *
     * @return  2*n_dof × 2*n_dof diagonal covariance matrix.
     */
    Eigen::MatrixXd make_goal_cov() const
    {
        const int nd = meta.n_dof;
        Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(2 * nd, 2 * nd);
        for (int d = 0; d < nd; ++d) {
            cov(d * 2,     d * 2)     = meta.goal_pos_noise_var;
            cov(d * 2 + 1, d * 2 + 1) = meta.goal_vel_noise_var;
        }
        return cov;
    }
};

// ============================================================
// OnlinePredictor
// ============================================================

/**
 * @brief Construct the predictor and load the trained model from disk.
 *
 * @param model_dir            Directory containing model.promp + metadata.txt.
 * @param max_obs_per_predict  0 = use ALL observations per predict() call.
 *                             N > 0 = use only the most recent N observations.
 *
 *   Tuning max_obs_per_predict:
 *     predict() cost grows as O(k × dims² × n_basis) where k is the number
 *     of observations used.  Cap k to 20–50 for sensors at 50–200 Hz to keep
 *     per-tick latency under 1 ms on modern hardware.
 *
 * @throws std::runtime_error if model.promp or metadata.txt cannot be read.
 */
OnlinePredictor::OnlinePredictor(const std::string& model_dir,
                                 int max_obs_per_predict)
    : _impl(std::make_unique<Impl>(model_dir, max_obs_per_predict))
{}

OnlinePredictor::~OnlinePredictor() = default;

/**
 * @brief Query the number of degrees of freedom this predictor expects.
 *
 * @return n_dof as set at training time (stored in metadata.txt).
 */
int OnlinePredictor::get_n_dof() const
{
    return _impl->meta.n_dof;
}

/**
 * @brief Clear the observation buffer and optionally override expected duration.
 *
 * Must be called at the start of each new motion trial.
 *
 * @param expected_duration_s
 *   Total expected motion time (s).
 *   ≤ 0 → use the training mean duration from metadata.txt.
 *
 *   Tuning: pass the actual expected duration if the current motion differs
 *   from the training demonstrations.  A mismatch shifts via-point phase
 *   steps and degrades prediction accuracy.
 *
 * ── Data-frequency note ─────────────────────────────────────────────────────
 * This call only clears the buffer; it does NOT change the SG window or any
 * noise parameters.  Those are fixed from training metadata.
 */
void OnlinePredictor::reset(double expected_duration_s)
{
    _impl->obs_times.clear();
    _impl->obs_pos.clear();
    _impl->expected_duration_s =
        (expected_duration_s > 0.0) ? expected_duration_s
                                    : _impl->meta.mean_duration_s;
    log_file.open("C:/Imperial/Final project/ProMP_realtime_Jorge-main/results/result.csv");


}

/**
 * @brief Record a new n_dof–dimensional sensor measurement.
 *
 * Stores the observation in an internal buffer.  Velocities are computed
 * lazily inside predict() using the full buffer, so this call is O(1).
 *
 * @param time_s        Elapsed time since motion start (s), monotonically
 *                      non-decreasing.  No assumption on spacing (rate-free).
 * @param positions_rad Joint positions (rad).  Size must equal get_n_dof().
 *
 * ── Data-frequency note ─────────────────────────────────────────────────────
 * The buffer stores actual timestamps; the SG filter uses them directly.
 * Observations may arrive at any rate including irregular bursts.  The SG
 * window (set at training time) determines the effective smoothing duration.
 */
void OnlinePredictor::add_observation(double time_s,
                                      const std::vector<double>& positions_rad)
{
    _impl->obs_times.push_back(time_s);
    _impl->obs_pos.push_back(positions_rad);
}

/**
 * @brief Return the number of observations since the last reset().
 */
int OnlinePredictor::num_observations() const
{
    return static_cast<int>(_impl->obs_times.size());
}

/**
 * @brief Generate a 2*n_dof–dimensional trajectory prediction.
 *
 * ── Algorithm ────────────────────────────────────────────────────────────────
 *   1. Reconstruct fresh ProMP from (μ_w, Σ_w) — trained model unchanged.
 *      time_mod = 1.0 because training demos are normalised to [0,1] phase.
 *   2. Compute joint velocities over the FULL observation buffer via SG
 *      (same window/order as training, loaded from metadata).
 *   3. Select the most recent min(n_avail, max_obs) observations.
 *   4. For each selected observation at elapsed time t_k:
 *        s_k = round(t_k / T_expected × (n_steps−1))
 *        via = [pos0_k, vel0_k, …, pos(n-1)_k, vel(n-1)_k] ∈ ℝ^(2*n_dof)
 *        Apply condition_via_point(s_k, via, Σ_obs).
 *   5. If PAST_TRAJ_TARGET:
 *        goal = [target0, 0, target1, 0, …] (zero velocities = natural stop)
 *        Apply condition_goal(goal, Σ_goal).
 *   6. generate_trajectory(n_future_steps) + gen_traj_std_dev(n_future_steps).
 *   7. Return slice [current_step : end] as PredictionResult.
 *
 * @param mode           PAST_TRAJ or PAST_TRAJ_TARGET.
 * @param n_future_steps Steps in the returned trajectory.
 *                       Tuning: 50 for fast display, 100–200 for planning.
 * @param target_pos     Target positions (rad), size = n_dof.
 *                       Empty vector → all zeros.
 *                       Ignored when mode == PAST_TRAJ.
 *
 * @return PredictionResult.  valid=false when fewer than 2 observations exist.
 *
 * ── Tuning summary ────────────────────────────────────────────────────────
 *   obs_pos_noise_var   – smaller → tighter fit to observed positions.
 *   obs_vel_noise_var   – larger  → more tolerance for SG velocity error.
 *   goal_pos_noise_var  – smaller → harder target constraint.
 *   max_obs_per_predict – reduce to bound per-tick latency at high rates.
 *   expected_duration_s – set via reset() to match actual motion speed.
 *
 * ── Data-frequency note ────────────────────────────────────────────────────
 * The phase mapping (step 4) uses real timestamps and is rate-independent.
 * SG velocity (step 2) uses actual timestamps and handles irregular rates.
 */

PredictionResult OnlinePredictor::predict(
    ConditioningMode mode,
    int n_future_steps,
    const std::vector<double>& target_pos) const
{
    PredictionResult result;
    result.valid       = false;
    result.n_obs_used  = 0;
    result.current_phase = 0.0;

    const Impl& im = *_impl;
    const int n_dof  = im.meta.n_dof;
    const int n_dims = 2 * n_dof;
    result.n_dims = n_dims;

    const int n_avail = static_cast<int>(im.obs_times.size());
    if (n_avail < 2) return result;

    // ── 1. Fresh ProMP copy from stored weight distribution ────────────────
    // time_mod = 1.0: demos were normalised to [0,1] phase during training.

    static promp::ProMP fresh(
        im.trained.get_weights(),
        im.trained.get_covariance(),
        im.trained.get_std_bf(),
        im.trained.get_n_samples(),
        im.trained.get_dims(),
        /*time_mod=*/1.0);

    // ── 3. Select observation window ──────────────────────────────────────
    // Use the most recent max_obs observations (or all if max_obs == 0).
    int start_idx = 0;
    if (im.max_obs > 0 && n_avail > im.max_obs)
        start_idx = n_avail - im.max_obs;

    // ── 2. Compute velocities over the full observation buffer via SG ──────
    //
    // vel_buf[d][i] = velocity estimate of joint d at observation i.
    // The SG filter uses the actual timestamps (im.obs_times), which makes
    // it correct for any sensor rate or irregular sampling pattern.
    //
    // Window/order are loaded from metadata (same values as training) so that
    // training and inference use identical smoothing characteristics.
    //const int sg_margin = im.meta.sg_window;  // generous margin
    //const int vel_start  = std::max(0, start_idx - sg_margin);
    //const int vel_count  = n_avail - vel_start;
    //std::vector<double> time_slice(im.obs_times.begin() + vel_start, im.obs_times.end());
/*
    std::vector<std::vector<double>> vel_buf(static_cast<size_t>(n_dof));
    {
        std::vector<double> pos_buf_d(static_cast<size_t>(n_avail));
        //std::vector<double> pos_buf_d(static_cast<size_t>(vel_count));
        for (int d = 0; d < n_dof; ++d) {
            for (int i = 0; i < n_avail; ++i)
            //for (int i = vel_start; i < n_avail; ++i)
            {
                //std::cout << "S"<< i<< " "<< n_avail << " " << vel_count << " " << im.obs_pos.size();
                pos_buf_d[static_cast<size_t>(i)] =
                //pos_buf_d[static_cast<size_t>(i-vel_start)] =
                    im.obs_pos[static_cast<size_t>(i)][static_cast<size_t>(d)];
                    //im.obs_pos[static_cast<size_t>(i)][static_cast<size_t>(d)];
                //std::cout << "F \n";
            }
            //std::cout << "velI ";
            vel_buf[static_cast<size_t>(d)] =
                sg_filter_deriv(im.obs_times, pos_buf_d,
                //sg_filter_deriv(time_slice, pos_buf_d,
                                im.meta.sg_window,
                                im.meta.sg_poly_order);
            //std::cout << "velOut \n";
        }
    }
    */
    const int sg_margin = im.meta.sg_window;  // generous margin
    const int vel_start  = std::max(0, start_idx - sg_margin);
    const int vel_count  = n_avail - vel_start;
    std::vector<double> time_slice(im.obs_times.begin() + vel_start, im.obs_times.end());
    std::vector<std::vector<double>> vel_buf(static_cast<size_t>(n_dof));
    {
        //std::vector<double> pos_buf_d(static_cast<size_t>(n_avail));
        std::vector<double> pos_buf_d(static_cast<size_t>(vel_count));
        for (int d = 0; d < n_dof; ++d) {
            //for (int i = 0; i < n_avail; ++i)
                for (int i = vel_start; i < n_avail; ++i)
            {
                //std::cout << "S"<< i<< " "<< n_avail << " " << vel_count << " " << im.obs_pos.size();
                //pos_buf_d[static_cast<size_t>(i)] =
                pos_buf_d[static_cast<size_t>(i-vel_start)] =
                    im.obs_pos[static_cast<size_t>(i)][static_cast<size_t>(d)];
                //im.obs_pos[static_cast<size_t>(i)][static_cast<size_t>(d)];
                //std::cout << "F \n";
            }
            //std::cout << "velI ";
            vel_buf[static_cast<size_t>(d)] =
                //sg_filter_deriv(im.obs_times, pos_buf_d,
                sg_filter_deriv(time_slice, pos_buf_d,
                                im.meta.sg_window,
                                im.meta.sg_poly_order);

            std::cout << vel_buf[d][vel_buf[0].size()-1] << ",";

        }
        std::cout <<"\n";
    }
    
    //std::cout << vel_buf[0].size() << "\n";
    // ── 4. Bayesian via-point conditioning ────────────────────────────────
    // For each selected observation:
    //   - Map elapsed time → phase step (time-based, rate-independent).
    //   - Build 2*n_dof via-point [pos0, vel0, pos1, vel1, …].
    //   - Apply Kalman-style update with block-diagonal Σ_obs.
    const Eigen::MatrixXd obs_cov = im.make_obs_cov();
    //std::cout << "Start conditioning \n";
    //std::cout << " n_avail " << n_avail;
    //std::cout << " start_idx " << start_idx;
    //std::cout << " vel_buf " << vel_buf[0].size();
    //std::cout << " vel_start " << vel_start;
    //std::cout << "\n";
    for (int i = start_idx; i < n_avail; ++i) {
        const int step = im.time_to_step(
            im.obs_times[static_cast<size_t>(i)]);

        Eigen::VectorXd via(n_dims);
        for (int d = 0; d < n_dof; ++d) {
            via(d * 2)     = im.obs_pos[static_cast<size_t>(i)][static_cast<size_t>(d)];
            via(d * 2 + 1) = vel_buf[static_cast<size_t>(d)][static_cast<size_t>(i-vel_start)];
        }
        fresh.condition_via_point(step, via, obs_cov);
    }
    result.n_obs_used = n_avail - start_idx;

    // Current phase from the latest observation time.
    const double t_last = im.obs_times.back();
    const double phase  = std::max(0.0,
        std::min(1.0, t_last / im.expected_duration_s));
    result.current_phase = phase;

    // ── 5. Goal conditioning (PAST_TRAJ_TARGET only) ──────────────────────
    // Target positions are pinned with tight noise; target velocity channels
    // are set to 0 rad/s with soft noise (natural deceleration profile).
    if (mode == ConditioningMode::PAST_TRAJ_TARGET) {
        Eigen::VectorXd goal(n_dims);
        for (int d = 0; d < n_dof; ++d) {
            goal(d * 2) =
                (static_cast<int>(target_pos.size()) > d)
                    ? target_pos[static_cast<size_t>(d)]
                    : 0.0;
            goal(d * 2 + 1) = 0.0; // target velocity = 0 (natural stop)
        }
        fresh.condition_goal(goal, im.make_goal_cov());
    }

    // ── 6. Generate full trajectory ───────────────────────────────────────
    const Eigen::MatrixXd mean_full =
        fresh.generate_trajectory(static_cast<size_t>(n_future_steps));
    const Eigen::MatrixXd std_full  =
        fresh.gen_traj_std_dev(static_cast<size_t>(n_future_steps));





    // ── 7. Slice from current phase onward ────────────────────────────────
    int current_step = static_cast<int>(
        std::round(phase * static_cast<double>(n_future_steps - 1)));
    current_step = std::max(0, std::min(n_future_steps - 1, current_step));
    const int n_out = n_future_steps - current_step;

    result.future_times_s.resize(static_cast<size_t>(n_out));
    result.mean_traj.resize(static_cast<size_t>(n_out * n_dims));
    result.std_traj .resize(static_cast<size_t>(n_out * n_dims));

    for (int i = 0; i < n_out; ++i) {
        const int gs = current_step + i;
        const double t_abs = static_cast<double>(gs) / (n_future_steps - 1)
                             * im.expected_duration_s;
        result.future_times_s[static_cast<size_t>(i)] = t_abs;

        log_file << result.current_phase<< ","
        << result.future_times_s[static_cast<size_t>(i)]<< ",";

        for (int dim = 0; dim < n_dims; ++dim) {
            result.mean_traj[static_cast<size_t>(i * n_dims + dim)] = mean_full(gs, dim);
            result.std_traj [static_cast<size_t>(i * n_dims + dim)] = std_full (gs, dim);

            log_file  << result.mean_traj[static_cast<size_t>(i * n_dims + dim)] << ","
                      << result.std_traj [static_cast<size_t>(i * n_dims + dim)] << "," ;
        }
        log_file << "\n";
    }
    // Save full trajectory
    /*
    const int n     = static_cast<int>(result.future_times_s.size());
    const int ndims = result.n_dims;
    for (int i = 0; i < n; ++i) {
        log_file << result.current_phase<< ","
        << result.future_times_s[static_cast<size_t>(i)]<< ",";

        for (int dim = 0; dim < ndims; ++dim) {
            // Interleave mean and std: mean[dim] std[dim]
            log_file  << result.mean_traj[static_cast<size_t>(i * ndims + dim)] << ","
                      << result.std_traj [static_cast<size_t>(i * ndims + dim)] << "," ;

        }
        log_file << "\n";
    }
    */

    result.valid = true;
    return result;
}

/**
 * @brief Convenience wrapper: add_observation() then predict().
 *
 * Intended for single-threaded polling loops where every sensor tick should
 * immediately yield an updated prediction.
 *
 * @param time_s        Elapsed time (s).
 * @param positions_rad Joint positions (rad); size must equal get_n_dof().
 * @param mode          PAST_TRAJ or PAST_TRAJ_TARGET.
 * @param n_future_steps  Future trajectory resolution (steps).
 * @param target_pos    Target positions (rad); empty → zeros.
 *
 * @return PredictionResult (see predict() documentation for details).
 */
PredictionResult OnlinePredictor::update_and_predict(
    double time_s,
    const std::vector<double>& positions_rad,
    ConditioningMode mode,
    int n_future_steps,
    const std::vector<double>& target_pos)
{
    add_observation(time_s, positions_rad);
    return predict(mode, n_future_steps, target_pos);
}

} // namespace promp_rt
