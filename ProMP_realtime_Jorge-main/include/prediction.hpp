#pragma once

/**
 * @file prediction.hpp
 * @brief Public API for ProMP_1 online real-time trajectory prediction.
 *
 * ── Usage pattern ────────────────────────────────────────────────────────────
 *   1. Train once offline with train_from_folder() → saves model.promp +
 *      metadata.txt to a model directory.
 *   2. Construct OnlinePredictor(model_dir) to load the trained model.
 *   3. Call reset() at the start of each new motion trial.
 *   4. At every sensor tick call update_and_predict() with the new joint
 *      positions; read PredictionResult for the predicted future trajectory.
 *
 * ── Data-frequency compatibility ─────────────────────────────────────────────
 * Observations may arrive at ANY rate and with non-uniform timestamps.
 * The predictor is time-stamp–based throughout:
 *   - Phase mapping:  phase = time_s / expected_duration_s   (no rate assumption)
 *   - SG velocity:    uses actual timestamps in the buffer    (rate-independent)
 * The SG window parameter (from metadata.txt) should be matched to the
 * expected sensor rate at training time (see TrainingConfig::sg_window).
 *
 * ── DoF configuration ────────────────────────────────────────────────────────
 * n_dof is determined at training time and stored in metadata.txt.
 * The predictor loads it automatically; call get_n_dof() to query it.
 * All observation vectors and target vectors must have exactly n_dof elements.
 *
 * ── Dimension layout in result buffers ───────────────────────────────────────
 * PredictionResult::mean_traj and std_traj are flat row-major arrays of size
 * (n_steps × 2*n_dof) with interleaved position/velocity columns:
 *
 *   index = step * (2*n_dof) + dof*2 + channel
 *   channel 0 → position (rad),   channel 1 → velocity (rad/s)
 *
 * Convenience accessors mean_pos/mean_vel/std_pos/std_vel handle the indexing.
 */

#include <memory>
#include <string>
#include <vector>

namespace promp_rt {

// ============================================================
// Conditioning mode
// ============================================================

/**
 * @brief Selects which information is used to condition the ProMP at each tick.
 */
enum class ConditioningMode {
    /**
     * @brief Condition on accumulated past observations only.
     *
     * The Bayesian update incorporates each past (time, position, velocity)
     * tuple as a via-point constraint.  No information about the future is
     * needed.  Use when the final target is unknown.
     */
    PAST_TRAJ,

    /**
     * @brief Condition on past observations AND a known final target position.
     *
     * In addition to all past via-point constraints, condition_goal() is
     * called with the target positions at the last trajectory step.  The
     * target velocity is assumed 0 rad/s (natural stop).  Use when the
     * intended end-point is known in advance (e.g. pick-and-place).
     */
    PAST_TRAJ_TARGET
};

// ============================================================
// Prediction result
// ============================================================

/**
 * @brief Output of a single predict() call.
 *
 * Flat storage layout (row-major):
 *   mean_traj[ step * n_dims + dof*2 + 0 ] = predicted mean position  of joint dof at step
 *   mean_traj[ step * n_dims + dof*2 + 1 ] = predicted mean velocity  of joint dof at step
 *   std_traj [ step * n_dims + dof*2 + 0 ] = 1-sigma std of position  of joint dof at step
 *   std_traj [ step * n_dims + dof*2 + 1 ] = 1-sigma std of velocity  of joint dof at step
 *
 * Use the convenience accessors (mean_pos, mean_vel, std_pos, std_vel)
 * instead of raw indexing to avoid off-by-one errors.
 */
struct PredictionResult {
    /**
     * @brief Total output channels = 2 * n_dof (positions + velocities).
     * Set by predict(); read-only after that.
     */
    int n_dims = 0;

    /**
     * @brief Absolute time stamps of each predicted step (seconds), spaced one
     * recording period apart.  Spans the receding horizon from the current time
     * to the end of the motion, so the count shrinks as current_phase grows.
     */
    std::vector<double> future_times_s;

    /**
     * @brief Predicted mean trajectory; size = future_times_s.size() * n_dims.
     * Indexed as mean_traj[step * n_dims + dim].
     * Use mean_pos() / mean_vel() accessors for clarity.
     * The latest observation is conditioned with a tight position variance, so
     * the first sample's position starts at (essentially) the current measured
     * position — no post-hoc offset is applied.
     */
    std::vector<double> mean_traj;

    /**
     * @brief Predicted 1-sigma standard deviations; same size and layout as mean_traj.
     * Use std_pos() / std_vel() accessors.
     */
    std::vector<double> std_traj;

    /** @brief Current phase in [0, 1]; derived from the latest observation time. */
    double current_phase = 0.0;

    /** @brief Number of past observations used for conditioning in this call. */
    int n_obs_used = 0;

    /**
     * @brief true when a valid prediction was produced.
     * false when fewer than 2 observations have been accumulated since reset().
     * All other fields are undefined when valid == false.
     */
    bool valid = false;

    // ── Convenience accessors ─────────────────────────────────────────────

    /** @brief Predicted mean position of joint @p dof at future step @p step (rad). */
    double mean_pos(int step, int dof) const
    { return mean_traj[static_cast<size_t>(step * n_dims + dof * 2)]; }

    /** @brief Predicted mean velocity of joint @p dof at future step @p step (rad/s). */
    double mean_vel(int step, int dof) const
    { return mean_traj[static_cast<size_t>(step * n_dims + dof * 2 + 1)]; }

    /** @brief Predicted 1-sigma position std of joint @p dof at future step @p step (rad). */
    double std_pos(int step, int dof) const
    { return std_traj[static_cast<size_t>(step * n_dims + dof * 2)]; }

    /** @brief Predicted 1-sigma velocity std of joint @p dof at future step @p step (rad/s). */
    double std_vel(int step, int dof) const
    { return std_traj[static_cast<size_t>(step * n_dims + dof * 2 + 1)]; }
};

// ============================================================
// OnlinePredictor
// ============================================================

/**
 * @brief Stateful online trajectory predictor based on a trained ProMP_1 model.
 *
 * ── Thread safety ────────────────────────────────────────────────────────────
 * Not thread-safe.  If you call add_observation() from one thread and
 * predict() from another, protect all calls with a mutex.
 *
 * ── Memory ───────────────────────────────────────────────────────────────────
 * The trained ProMP weight matrices (μ_w, Σ_w) are loaded once in the
 * constructor and kept in memory.  A fresh copy is created per predict() call
 * for conditioning; the stored trained model is never mutated.
 *
 * ── Computation per tick ─────────────────────────────────────────────────────
 * predict() cost scales as O(k × n_dims² × n_basis) where k is the number
 * of observations used.  Use max_obs_per_predict to cap k and bound latency.
 */
class OnlinePredictor {
public:
    /**
     * @brief Load the trained model and prepare for online prediction.
     *
     * @param model_dir
     *   Directory containing model.promp and metadata.txt, as written by
     *   train_from_folder().  Both files must be present.
     *
     * @param max_obs_per_predict
     *   Maximum number of past observations used per predict() call.
     *   0 (default) uses ALL accumulated observations.
     *   N > 0 uses only the most recent N observations.
     *
     *   Tuning:
     *     - 0        → highest accuracy, unbounded latency growth over time.
     *     - 20 – 50  → good balance for sensors at 50–200 Hz.
     *     - 1 – 10   → use when hard real-time deadlines are required.
     *
     * @throws std::runtime_error if model.promp or metadata.txt cannot be read.
     */
    explicit OnlinePredictor(const std::string& model_dir,
                             int max_obs_per_predict = 0);

    ~OnlinePredictor();

    // ── Configuration queries ─────────────────────────────────────────────

    /**
     * @brief Number of degrees of freedom expected by this predictor.
     *
     * Every call to add_observation() and predict() must supply exactly
     * this many joint positions / target values.
     *
     * @return n_dof as stored in metadata.txt (set at training time).
     */
    int get_n_dof() const;

    // ── Session management ────────────────────────────────────────────────

    /**
     * @brief Discard all accumulated observations and optionally override the
     *        expected motion duration.
     *
     * Call this at the start of each new motion trial.
     *
     * @param expected_duration_s
     *   Total expected motion time (seconds).  Used to map real elapsed time
     *   to a ProMP phase step:
     *     phase = clamp(elapsed_time_s / expected_duration_s, 0, 1)
     *     step  = round(phase × (n_steps − 1))
     *
     *   ≤ 0 → use the mean duration from the training demonstrations
     *          (stored in metadata.txt).
     *
     *   Tuning:
     *     If the current motion is consistently faster or slower than
     *     training, pass the correct value here rather than relying on the
     *     training mean.  A mismatch causes the conditioned via-points to be
     *     placed at wrong phase steps, degrading prediction accuracy.
     *
     * ── Data-frequency note ───────────────────────────────────────────────
     * This call resets the observation buffer; it does NOT change the SG
     * window or any noise parameters (those are fixed at training time).
     */
    void reset(double expected_duration_s = -1.0);

    // ── Online data ingestion ─────────────────────────────────────────────

    /**
     * @brief Record a new sensor measurement.
     *
     * Observations are stored in an internal buffer.  Velocities are
     * estimated lazily (inside predict()) from the full buffer using the
     * Savitzky-Golay filter, so this call is O(1).
     *
     * @param time_s
     *   Elapsed time since motion start (seconds).  Must be monotonically
     *   non-decreasing across successive calls within one trial.
     *
     *   ── Data-frequency note ─────────────────────────────────────────────
     *   Observations may arrive at any rate or with irregular timestamps.
     *   The SG filter and phase mapping both use actual time values.
     *   For best velocity accuracy, the SG window (set at training time)
     *   should correspond to at least 5 samples at the SLOWEST expected rate.
     *
     * @param positions_rad
     *   Joint positions in radians.  Must have exactly get_n_dof() elements.
     *   Element i is the position of joint i.
     *
     * @note Calling this with positions_rad.size() ≠ get_n_dof() is undefined
     *       behaviour.  The caller is responsible for ensuring correct size.
     */
    void add_observation(double time_s,
                         const std::vector<double>& positions_rad);

    // ── Prediction ────────────────────────────────────────────────────────

    /**
     * @brief Generate a trajectory prediction from all accumulated observations.
     *
     * Algorithm:
     *   1. Reconstruct a fresh ProMP from stored (μ_w, Σ_w) — trained model
     *      is not modified.
     *   2. Compute joint velocities over the observation buffer via the SG
     *      filter (same window/order as training).
     *   3. For each selected past observation, form a 2*n_dof via-point
     *        [pos0, vel0, pos1, vel1, ..., pos(n-1), vel(n-1)]
     *      and apply Bayesian conditioning:
     *        L   = Σ_w Φ(s) (Σ_obs + Φ(s)ᵀ Σ_w Φ(s))⁻¹
     *        μ_w ← μ_w + L (y − Φ(s)ᵀ μ_w)
     *        Σ_w ← Σ_w − L Φ(s)ᵀ Σ_w
     *   4. If PAST_TRAJ_TARGET: apply condition_goal() with target positions
     *      (tight noise) and zero target velocities (soft noise).
     *   5. Generate the full trajectory, return the slice from current phase.
     *
     * @param mode
     *   PAST_TRAJ or PAST_TRAJ_TARGET.
     *
     * @param record_dt_s
     *   Time step (seconds) between consecutive predicted samples.  Match this
     *   to the recording period of the incoming data so the prediction grid
     *   follows the raw time series rather than an arbitrary resolution.
     *   ≤ 0 (default) → derive it automatically from the observed timestamps
     *   (median inter-sample interval).  The returned trajectory then has a
     *   receding horizon: its length equals the remaining motion time and the
     *   step count is (expected_duration_s − t_current) / record_dt_s — the full
     *   duration at the start, shrinking toward zero as the motion completes.
     *
     * @param target_pos
     *   Target joint positions in radians; must have exactly get_n_dof()
     *   elements.  Ignored when mode == PAST_TRAJ.
     *   Pass {} (empty vector) to use zero targets.
     *   Default: empty → zeros.
     *
     * @return PredictionResult with valid=false when fewer than 2 observations
     *         have been accumulated since the last reset().
     *
     * ── Tuning summary ────────────────────────────────────────────────────
     *   obs_pos_noise_var  (set at training) – smaller → tighter position fit.
     *   obs_vel_noise_var  (set at training) – larger  → more tolerance for
     *                                          SG estimation error.
     *   goal_pos_noise_var (set at training) – smaller → harder position stop.
     *   max_obs_per_predict (constructor)    – cap to bound latency per tick.
     *   expected_duration_s (reset)          – match to actual motion speed.
     */
    PredictionResult predict(
        ConditioningMode mode,
        double record_dt_s = -1.0,
        const std::vector<double>& target_pos = {}) const;

    /**
     * @brief Convenience: add_observation() then predict() in one call.
     *
     * Equivalent to:
     *   add_observation(time_s, positions_rad);
     *   return predict(mode, record_dt_s, target_pos);
     *
     * @param time_s        Elapsed time (s) — see add_observation().
     * @param positions_rad Joint positions (rad) — see add_observation().
     * @param mode          Conditioning mode.
     * @param record_dt_s   Prediction step (s); ≤0 → derive from timestamps.
     * @param target_pos    Target positions (rad); empty → zeros.
     *
     * @return PredictionResult — see predict().
     */
    PredictionResult update_and_predict(
        double time_s,
        const std::vector<double>& positions_rad,
        ConditioningMode mode,
        double record_dt_s = -1.0,
        const std::vector<double>& target_pos = {});

    /**
     * @brief Number of observations accumulated since the last reset().
     * @return Observation count (≥ 0).
     */
    int num_observations() const;

private:
    struct Impl;
    std::unique_ptr<Impl> _impl;
};

} // namespace promp_rt
