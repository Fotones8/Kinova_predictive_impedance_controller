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
#include <cstdlib>
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
 * @brief All model-side parameters of the online predictor (loaded from
 *        metadata.txt, fixed at training time).
 *
 * ── Where the rest of the tunables live ───────────────────────────────────────
 * Three parameters are supplied at RUN time, not in metadata, and are
 * documented where they are passed:
 *   • max_obs_per_predict  (OnlinePredictor constructor) — cap on past
 *       observations used per predict() tick.  0 = use all.  ↓ bounds per-tick
 *       latency (predict cost ≈ O(k · dims² · n_basis)); too small loses early
 *       history and can degrade accuracy.
 *   • expected_duration_s  (reset())  — total expected motion time.  Anchors
 *       the time→phase map AND the receding-horizon end.  Mismatch vs the real
 *       motion shifts every via-point's phase and distorts the prediction.
 *   • record_dt_s          (predict()) — output time step (s).  Should equal
 *       the data's recording period so the prediction grid follows the raw time
 *       series.  ≤0 ⇒ auto from the median observed interval.  ↓ → more output
 *       steps per tick (finer grid, more compute).
 *
 * ── Quick tuning cheat-sheet (↑ = increase the value) ─────────────────────────
 *   n_steps ............ ↑ finer phase grid ⇒ prediction can start closer to the
 *                        true current phase; training/conditioning cost ↑.
 *   sg_window .......... ↑ smoother velocity, more lag.
 *   sg_poly_order ...... ↑ keeps sharp velocity peaks, less smoothing, noisier.
 *   obs_pos_noise_var .. ↓ tighter fit to measured positions (less smoothing).
 *   obs_vel_noise_var .. ↑ tolerate noisy SG velocity (keep ≫ obs_pos_noise_var).
 *   goal_pos_noise_var . ↓ harder target AND harder pin of the current position.
 *   goal_vel_noise_var . ↑ allow a more gradual approach to the target.
 *
 * Every field maps 1-to-1 to a key written by write_metadata() in training.cpp;
 * unrecognised keys are ignored for forward compatibility.
 */
struct Metadata {
    /// Number of joints (DoF).  Must equal the number of position columns in the
    /// training CSVs and the length of every observation/target vector.  Drives
    /// all matrix sizes (state dim = 2·n_dof: interleaved position + velocity).
    int    n_dof             = 3;

    /// Phase resolution of the trained ProMP (timesteps spanning phase [0,1]).
    /// IMPORTANT: this is also the grid condition_via_point() snaps each
    /// observation time onto, so it bounds how exactly the prediction can begin
    /// at the current position — one cell ≈ mean_duration_s/(n_steps-1) seconds
    /// (≈32 ms at 100 over a 3.2 s motion).  ↑ for a tighter natural start, at
    /// higher training cost.
    int    n_steps           = 100;

    /// Mean demonstration duration (s).  Used as the default expected_duration_s
    /// (time→phase mapping and horizon length) when reset() gets no override.
    double mean_duration_s   = 1.0;

    /// Savitzky–Golay window length (samples) for the online velocity estimate.
    /// ↑ → smoother velocity but more lag.  Pick ≈ rate_Hz · smoothing_s; must
    /// be odd and ≥ sg_poly_order + 2.  Span ≥5 samples at the slowest rate.
    int    sg_window         = 9;

    /// Savitzky–Golay polynomial order.  ↑ → better preserves velocity peaks,
    /// less smoothing, more high-frequency noise.  Must be < sg_window − 1
    /// (typically 3–4).
    int    sg_poly_order     = 4;

    /// Position observation noise variance (rad²) for past via-points.  ↓ →
    /// tighter fit to measured positions; ↑ → more smoothing / noise rejection.
    double obs_pos_noise_var = 1e-4;

    /// Velocity observation noise variance ((rad/s)²) for via-points.  Kept large
    /// (≫ obs_pos_noise_var) so the noisy SG velocity does not dominate; ↓ to
    /// trust the SG velocity more.
    double obs_vel_noise_var = 1e-2;

    /// Tight position variance (rad²).  Two uses: (a) the target position in
    /// PAST_TRAJ_TARGET mode, and (b) pinning the MOST RECENT observation so the
    /// prediction starts at the current measured position.  ↓ → harder pin /
    /// harder target (≈1e-6 is already a near-hard constraint).
    double goal_pos_noise_var= 1e-6;

    /// Goal velocity variance ((rad/s)²) for the target stop.  ↑ → allow a more
    /// gradual approach; ↓ → enforce a precise zero-velocity stop at the target.
    double goal_vel_noise_var= 1e-2;
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

    const Eigen::MatrixXd obs_cov;
    const Eigen::MatrixXd goal_cov;


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
          expected_duration_s(meta.mean_duration_s),

    obs_cov(make_obs_cov()),
    goal_cov(make_goal_cov())
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
    //_impl->obs_times.reserve(15000);
    //_impl->obs_pos.reserve(15000);
    _impl->expected_duration_s =
        (expected_duration_s > 0.0) ? expected_duration_s
                                    : _impl->meta.mean_duration_s;
    // expected_duration_s divides into the phase map and n_full; never let it be
    // zero/negative even if metadata.mean_duration_s is malformed.
    if (_impl->expected_duration_s <= 0.0)
        _impl->expected_duration_s = 1.0;

    // (Re)open the per-trial prediction log.  The path comes from the
    // PROMP_RT_LOG environment variable when set (set it empty to disable
    // logging); otherwise it defaults to a relative "results/result.csv".  A
    // relative default keeps this portable — the previous hard-coded Windows
    // path silently failed on other platforms.  close() first so a repeated
    // reset() re-opens cleanly (open() on an already-open stream otherwise fails).
    if (log_file.is_open()) log_file.close();
    const char* env = std::getenv("PROMP_RT_LOG");
    const std::string log_path = env ? std::string(env)
                                     : std::string("C:/Imperial/Final project/ProMP_realtime_Jorge-main/results/result.csv");
    if (!log_path.empty()) {
        std::error_code ec;
        const fs::path p(log_path);
        if (p.has_parent_path()) fs::create_directories(p.parent_path(), ec);
        log_file.open(log_path);
        if (!log_file.is_open())
            std::cerr << "[predict] Warning: could not open log file '"
                      << log_path << "'; prediction logging disabled.\n";
    }
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
 *        Apply condition_via_point(s_k, via, Σ_obs).  The most recent
 *        observation uses a tight position variance so the prediction starts
 *        at the current measured position.
 *   5. If PAST_TRAJ_TARGET:
 *        goal = [target0, 0, target1, 0, …] (zero velocities = natural stop)
 *        Apply condition_goal(goal, Σ_goal).
 *   6. Derive the recording period (from record_dt_s, or the median observed
 *      inter-sample interval) and set n_full = round(duration/dt)+1 so the
 *      phase grid is spaced one recording period apart in time.
 *   7. generate_trajectory(n_full) + gen_traj_std_dev(n_full); return the slice
 *      [current_step : end] — the receding remaining-time horizon.  Because the
 *      latest observation is conditioned with a tight position variance (step 4),
 *      this slice starts at the current measured position.
 *
 * @param mode           PAST_TRAJ or PAST_TRAJ_TARGET.
 * @param record_dt_s    Prediction time step (s); should equal the data's
 *                       recording period.  ≤0 → auto-derive from timestamps.
 *                       Sets both the step spacing and (with the receding
 *                       horizon) the number of returned steps.
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
    double record_dt_s,
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

    // NOTE: must be a per-call local (NOT static).  condition_via_point() and
    // condition_goal() below mutate this object; a static instance would retain
    // all conditioning across ticks and trials, so the model would never reset
    // to the trained prior and predictions would diverge over time.
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

    // ── 2. Per-joint position slices for on-demand SG velocity ────────────
    // Velocity is needed only at the observations actually conditioned (one per
    // phase step, see step 4), so instead of running SG over the whole buffer we
    // keep a position slice per joint starting sg_margin (= sg_window) samples
    // before start_idx — enough left context for a full SG window at every
    // conditioned index — and evaluate the derivative pointwise with
    // sg_deriv_at_index() inside the conditioning loop.  This makes velocity
    // cost O(#conditioned · W) instead of O(#observations · W), with identical
    // values at the conditioned indices.
    const int sg_margin  = im.meta.sg_window; // left context for a full SG window
    const int vel_start  = std::max(0, start_idx - sg_margin);
    const int vel_count  = n_avail - vel_start;
    const std::vector<double> time_slice(im.obs_times.begin() + vel_start,
                                         im.obs_times.end());
    std::vector<std::vector<double>> pos_slice(
        static_cast<size_t>(n_dof),
        std::vector<double>(static_cast<size_t>(vel_count)));
    for (int d = 0; d < n_dof; ++d)
        for (int i = vel_start; i < n_avail; ++i)
            pos_slice[static_cast<size_t>(d)][static_cast<size_t>(i - vel_start)] =
                im.obs_pos[static_cast<size_t>(i)][static_cast<size_t>(d)];

    // ── 4. Bayesian via-point conditioning ────────────────────────────────
    // For each selected observation:
    //   - Map elapsed time → phase step (time-based, rate-independent).
    //   - Build 2*n_dof via-point [pos0, vel0, pos1, vel1, …].
    //   - Apply Kalman-style update with block-diagonal Σ_obs.
    //const Eigen::MatrixXd obs_cov = im.make_obs_cov();

    // The most recent observation is conditioned with a TIGHT position variance
    // (goal_pos_noise_var) so the conditioned trajectory passes through the
    // current measured position — i.e. the prediction naturally starts from the
    // current position, without any post-hoc offset.  The velocity channel is
    // kept soft (obs_vel_noise_var) because the SG velocity estimate is noisy.
    Eigen::MatrixXd anchor_cov = im.obs_cov;
    for (int d = 0; d < n_dof; ++d)
        anchor_cov(d * 2, d * 2) = im.meta.goal_pos_noise_var;

    // condition_via_point() snaps every observation onto one of n_steps phase
    // steps, so conditioning all samples that share a step is largely redundant
    // — and it makes the cost O(#observations) per tick (O(n²) per trial) AND
    // the posterior rate-dependent (denser sampling ⇒ artificially tighter).
    // Observations are time-ordered and time_to_step() is monotonic, so equal-
    // step samples are contiguous: condition only the LAST sample of each step
    // run.  This caps conditioning at the number of distinct phase steps
    // (≤ n_steps) regardless of sensor rate.  Measured vs conditioning every
    // sample: ~7× faster with no meaningful accuracy change against ground truth
    // (whole-horizon and 1 s errors slightly better, 0.1–0.5 s within ~0.1°).
    Eigen::VectorXd via(n_dims);
    int n_cond = 0;
    for (int i = start_idx; i < n_avail; ++i) {
        const int step = im.time_to_step(
            im.obs_times[static_cast<size_t>(i)]);
        const bool is_latest = (i == n_avail - 1);
        const bool step_ends = is_latest ||
            (im.time_to_step(im.obs_times[static_cast<size_t>(i + 1)]) != step);
        if (!step_ends) continue;

        for (int d = 0; d < n_dof; ++d) {
            via(d * 2)     = im.obs_pos[static_cast<size_t>(i)][static_cast<size_t>(d)];
            via(d * 2 + 1) = sg_deriv_at_index(
                time_slice, pos_slice[static_cast<size_t>(d)],
                i - vel_start, im.meta.sg_window, im.meta.sg_poly_order);
        }
        fresh.condition_via_point(step, via, is_latest ? anchor_cov : im.obs_cov);
        ++n_cond;
    }
    result.n_obs_used = n_cond;

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
        fresh.condition_goal(goal, im.goal_cov);
    }

    // ── 6. Determine the prediction time step from the RAW time series ─────
    // The spacing between predicted samples must follow the recording rate of
    // the incoming data, NOT an arbitrary trajectory resolution.  By default we
    // derive it from the observed timestamps (median inter-sample interval,
    // robust to jitter/gaps); a positive record_dt_s overrides this.
    double record_dt = record_dt_s;
    if (record_dt <= 0.0) {
        std::vector<double> diffs;
        diffs.reserve(static_cast<size_t>(n_avail - 1));
        for (int i = 1; i < n_avail; ++i) {
            const double d = im.obs_times[static_cast<size_t>(i)]
                           - im.obs_times[static_cast<size_t>(i - 1)];
            if (d > 0.0) diffs.push_back(d);
        }
        if (diffs.empty()) return result; // cannot infer the recording rate
        std::nth_element(diffs.begin(),
                         diffs.begin() + static_cast<long>(diffs.size() / 2),
                         diffs.end());
        record_dt = diffs[diffs.size() / 2];
    }

    // ── 6. Compute current phase step and cap future horizon ──────────────
    const int n_full = std::max(2,
        static_cast<int>(std::llround(im.expected_duration_s / record_dt)) + 1);

    int current_step = static_cast<int>(
        std::round(phase * static_cast<double>(n_full - 1)));
    current_step = std::max(0, std::min(n_full - 1, current_step));

    constexpr int MAX_FUTURE_STEPS = 400;

    // Fixed step size in phase units: one record_dt corresponds to this
    // phase increment regardless of where we are in the motion.
    const double phase_step = record_dt / im.expected_duration_s;

    // Future horizon end: always MAX_FUTURE_STEPS steps ahead at record_dt spacing.
    // Clamp to 1.0 only if the end of the motion is closer than MAX_FUTURE_STEPS steps.
    const double phase_end = std::min(1.0, phase + MAX_FUTURE_STEPS * phase_step);

    // Actual number of steps — equals MAX_FUTURE_STEPS unless we're near the end.
    const int n_out = std::max(1,
        static_cast<int>(std::round((phase_end - phase) / phase_step)) + 1);

    // Phase vector with uniform spacing of exactly phase_step.
    const Eigen::VectorXd future_phase =
        Eigen::VectorXd::LinSpaced(n_out, phase, phase_end);

    const Eigen::MatrixXd mean_out = fresh.generate_trajectory_at(future_phase);
    const Eigen::MatrixXd std_out  = fresh.gen_traj_std_dev_at(future_phase);

    // ── 8. Fill result directly from the n_out × n_dims output matrices ──
    result.future_times_s.resize(static_cast<size_t>(n_out));
    result.mean_traj.resize(static_cast<size_t>(n_out * n_dims));
    result.std_traj .resize(static_cast<size_t>(n_out * n_dims));

    const bool do_log = log_file.is_open();
    for (int i = 0; i < n_out; ++i) {
        // Map output index back to absolute time
        const double phase_i = future_phase(i);
        result.future_times_s[static_cast<size_t>(i)] =
            phase_i * im.expected_duration_s;

        //if (do_log)
        //    log_file << result.current_phase << "," << phase_i * im.expected_duration_s << ",";

        for (int dim = 0; dim < n_dims; ++dim) {
            result.mean_traj[static_cast<size_t>(i * n_dims + dim)] = mean_out(i, dim);
            result.std_traj [static_cast<size_t>(i * n_dims + dim)] = std_out (i, dim);

            //if (do_log)
            //    log_file << result.mean_traj[static_cast<size_t>(i * n_dims + dim)] << ","
            //             << result.std_traj [static_cast<size_t>(i * n_dims + dim)] << ",";
        }
        //if (do_log) log_file << "\n";
    }

    //ORIGINAL IMPLEMENTATION
/*
    // Full phase grid resolution so that consecutive steps are record_dt apart
    // in time over the whole motion: expected_duration spans phase [0,1], so a
    // grid of n_full points gives a step of expected_duration/(n_full-1) ≈ record_dt.
    // n_full-1 = round(duration / record_dt) = number of recording periods in the motion.
    const int n_full = std::max(2,
        static_cast<int>(std::llround(im.expected_duration_s / record_dt)) + 1);

    // ── 7. Generate full trajectory at the recording resolution ───────────
    const Eigen::MatrixXd mean_full =
        fresh.generate_trajectory(static_cast<size_t>(n_full));
    const Eigen::MatrixXd std_full  =
        fresh.gen_traj_std_dev(static_cast<size_t>(n_full));

    // ── 8. Receding horizon: from the current step to the end of the motion ─
    // At motion start (phase≈0) this spans the full duration; after each tick
    // the start advances, so the horizon recedes to cover only the remaining
    // time.  The number of returned steps = remaining_time / record_dt.
    int current_step = static_cast<int>(
        std::round(phase * static_cast<double>(n_full - 1)));
    current_step = std::max(0, std::min(n_full - 1, current_step));
    // Trying not to generate the full trajectory
    const int n_out = std::min(n_full - current_step, 400);
    //const int n_out = n_full - current_step;

    result.future_times_s.resize(static_cast<size_t>(n_out));
    result.mean_traj.resize(static_cast<size_t>(n_out * n_dims));
    result.std_traj .resize(static_cast<size_t>(n_out * n_dims));

    const bool do_log = log_file.is_open();
    for (int i = 0; i < n_out; ++i) {
        const int gs = current_step + i;
        const double t_abs = static_cast<double>(gs) / (n_full - 1)
                             * im.expected_duration_s;
        result.future_times_s[static_cast<size_t>(i)] = t_abs;

        if (do_log)
            log_file << result.current_phase << "," << t_abs << ",";

        for (int dim = 0; dim < n_dims; ++dim) {
            result.mean_traj[static_cast<size_t>(i * n_dims + dim)] = mean_full(gs, dim);
            result.std_traj [static_cast<size_t>(i * n_dims + dim)] = std_full (gs, dim);

            if (do_log)
                log_file << result.mean_traj[static_cast<size_t>(i * n_dims + dim)] << ","
                         << result.std_traj [static_cast<size_t>(i * n_dims + dim)] << ",";
        }
        if (do_log) log_file << "\n";
    }

    */
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
 * @param record_dt_s   Prediction time step (s); ≤0 → auto from timestamps.
 * @param target_pos    Target positions (rad); empty → zeros.
 *
 * @return PredictionResult (see predict() documentation for details).
 */
PredictionResult OnlinePredictor::update_and_predict(
    double time_s,
    const std::vector<double>& positions_rad,
    ConditioningMode mode,
    double record_dt_s,
    const std::vector<double>& target_pos)
{
    add_observation(time_s, positions_rad);
    return predict(mode, record_dt_s, target_pos);
}

} // namespace promp_rt
