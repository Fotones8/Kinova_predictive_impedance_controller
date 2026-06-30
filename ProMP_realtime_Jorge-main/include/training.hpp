#pragma once

/**
 * @file training.hpp
 * @brief Public API for ProMP_1 offline training from CSV demonstrations.
 *
 * ── Overview ────────────────────────────────────────────────────────────────
 * Each demonstration is loaded from a CSV file whose columns are:
 *
 *   time_s,  pos0_rad,  pos1_rad, ...,  pos(n_dof-1)_rad
 *
 * Velocities for every joint are derived internally using the Savitzky-Golay
 * (SG) filter, so NO velocity column is required in the input.
 *
 * The ProMP is trained on 2*n_dof–dimensional trajectories laid out as:
 *   [pos0, vel0,  pos1, vel1,  ...,  pos(n-1), vel(n-1)]
 *
 * ── Data-frequency compatibility ────────────────────────────────────────────
 * Different demonstrations are allowed to have different sampling rates or
 * irregular time stamps.  All demos are resampled to a common uniform grid of
 * n_resample_steps points via linear interpolation AFTER velocity computation,
 * so velocity quality is tied to the original data frequency.
 *
 * Rule of thumb for the SG window relative to data frequency:
 *   sg_window ≈ data_rate_Hz × desired_smoothing_time_s   (rounded to odd)
 *   Example: 100 Hz, 90 ms smoothing → sg_window = 9
 *   Example: 500 Hz, 90 ms smoothing → sg_window = 45 → round to 45 (odd)
 */

#include <string>

namespace promp_rt {

/**
 * @brief All hyper-parameters that control offline ProMP training.
 *
 * Create one instance, adjust the fields you care about, and pass it to
 * train_from_folder().  Defaults are chosen for 3-DoF joint trajectories
 * sampled at roughly 100 Hz with motions lasting 1–2 s.
 */
struct TrainingConfig {

    // ── Dimensionality ──────────────────────────────────────────────────────

    /**
     * @brief Number of degrees of freedom (joints) in the input data.
     *
     * Each CSV file must contain exactly (n_dof + 1) columns:
     *   [time_s,  pos0_rad,  pos1_rad,  ...,  pos(n_dof-1)_rad]
     *
     * The trained ProMP will have 2*n_dof output dimensions (position +
     * velocity per joint).
     *
     * Range: 1 – any positive integer.
     * Common values: 1 (single joint), 3 (wrist/elbow/shoulder), 6 (full arm).
     */
    int n_dof = 3;

    // ── ProMP model ─────────────────────────────────────────────────────────

    /**
     * @brief Number of Gaussian RBF basis functions.
     *
     * Controls the complexity of the learned trajectory shape.
     *   - Too few → over-smoothed mean; complex motions are not captured.
     *   - Too many → over-fitting on few demonstrations; large weight matrix.
     *
     * Tuning guideline:
     *   Start with 8–12 for smooth manipulation tasks.
     *   Increase toward 20 for trajectories with multiple sub-movements.
     *   The total weight vector length is (num_basis_functions * 2 * n_dof).
     *
     * Default: 10.
     */
    int num_basis_functions = 10;

    /**
     * @brief Standard deviation of each RBF basis function (variance actually
     *        used internally by the promp_1 library — see promp.cpp).
     *
     * Controls the width of each basis function in normalised phase [0, 1].
     *   - ≤ 0 → auto-tuned to 1 / (num_basis_functions²)  [recommended].
     *   - > 0 → used as-is.
     *
     * When to set manually:
     *   - If the mean trajectory is jagged, the basis functions are too narrow;
     *     increase std_bf (e.g. 0.01 – 0.05).
     *   - If the mean is too smooth and peak positions are wrong, they are too
     *     wide; decrease std_bf.
     *
     * Default: -1 (auto).
     */
    double std_bf = -1.0;

    /**
     * @brief Number of uniformly-spaced steps in the resampled training
     *        trajectory (the common phase grid).
     *
     * All demonstrations are interpolated to this length after velocity
     * computation.  Higher values give finer phase resolution and slightly
     * better accuracy but increase memory and training time.
     *
     * Rule of thumb: 50–100 for motions up to 2 s; 200 for fast dynamics.
     *
     * Default: 100.
     */
    int n_resample_steps = 100;

    /**
     * @brief Discard data beyond this elapsed time in each demonstration (s).
     *
     * Set to the expected motion duration to remove post-motion noise or
     * artefacts.  ≤ 0 means keep all data.
     *
     * Default: -1 (keep all).
     */
    double max_demo_duration = -1.0;

    // ── Savitzky-Golay velocity estimation ─────────────────────────────────

    /**
     * @brief Window length for the Savitzky-Golay derivative filter (must be odd).
     *
     * The SG filter fits a polynomial of degree sg_poly_order to a sliding
     * window of sg_window raw samples and evaluates the analytical first
     * derivative to obtain velocity.
     *
     * Effect:
     *   - Larger window → more noise suppression, slight peak blurring.
     *   - Smaller window → less smoothing, more sensitive to measurement noise.
     *
     * Frequency scaling:
     *   Choose sg_window so that (sg_window / data_rate_Hz) covers the
     *   desired smoothing time (typically 50–100 ms).
     *   - 100 Hz data, 90 ms smoothing → sg_window = 9
     *   - 500 Hz data, 90 ms smoothing → sg_window = 45 (next odd = 45)
     *   - 50  Hz data, 90 ms smoothing → sg_window = 5
     *
     * Constraints: odd, ≥ sg_poly_order + 2.
     * Default: 9  (suitable for ~100 Hz data).
     */
    int sg_window = 9;

    /**
     * @brief Polynomial order for the Savitzky-Golay filter.
     *
     * Higher order → more accurate velocity peaks, but more sensitive to
     * high-frequency noise.  Must satisfy sg_poly_order < sg_window - 1.
     *
     * Recommended values: 3 or 4.
     * Default: 4.
     */
    int sg_poly_order = 4;

    // ── Bayesian conditioning noise (used at inference, saved in metadata) ──

    /**
     * @brief Observation noise variance for the POSITION channels (rad²)
     *        used during Bayesian via-point conditioning at inference time.
     *
     * Appears as the diagonal entry σ²_pos in the 2n_dof × 2n_dof observation
     * noise covariance Σ_obs:
     *   Σ_obs = diag(σ²_pos, σ²_vel, σ²_pos, σ²_vel, ...).
     *
     * Effect:
     *   - Smaller → the conditioned trajectory passes very close to each
     *     observed position (hard constraint).
     *   - Larger  → observations are trusted less; prior (training mean) has
     *     more influence.
     *
     * Tuning: match to the RMS position measurement noise of your sensor.
     * Typical range: 1e-5 (high-quality encoder) to 1e-3 (noisy IMU).
     * Default: 1e-4.
     */
    double obs_pos_noise_var = 1e-4;

    /**
     * @brief Observation noise variance for the VELOCITY channels ((rad/s)²)
     *        used during Bayesian via-point conditioning at inference time.
     *
     * Velocity is not directly measured but estimated by the SG filter;
     * setting this larger than obs_pos_noise_var accounts for that estimation
     * error.
     *
     * Effect: larger value → velocity estimates have less influence on the
     * conditioned weight distribution.
     *
     * Tuning: typically 10–1000× larger than obs_pos_noise_var.
     * Default: 1e-2.
     */
    double obs_vel_noise_var = 1e-2;

    /**
     * @brief Observation noise variance for the goal POSITION channels (rad²)
     *        used in PAST_TRAJ_TARGET conditioning mode.
     *
     * Smaller value → the final predicted position is pinned tightly to the
     * target.  Increase if a soft approach rather than a hard stop is desired.
     *
     * Typical range: 1e-7 (very tight) to 1e-5 (soft).
     * Default: 1e-6.
     */
    double goal_pos_noise_var = 1e-6;

    /**
     * @brief Observation noise variance for the goal VELOCITY channels ((rad/s)²).
     *
     * The target velocity at the goal is implicitly assumed to be 0 rad/s
     * (natural stop).  A large value here allows the model to find a natural
     * deceleration profile rather than forcing zero velocity exactly.
     *
     * Default: 1e-2.
     */
    double goal_vel_noise_var = 1e-2;
};

/**
 * @brief Train a ProMP from a folder of CSV demonstrations and save the model.
 *
 * ── Input ────────────────────────────────────────────────────────────────────
 * @param demo_folder
 *   Path to a directory containing one or more *.csv files.
 *   Every CSV must have exactly (config.n_dof + 1) columns:
 *     time_s,  pos0_rad,  pos1_rad,  ...,  pos(n_dof-1)_rad
 *   An optional header row (text) is auto-detected and skipped.
 *   Demos with fewer than 5 valid rows are silently discarded.
 *   Different demos may have different sampling rates or lengths.
 *
 * @param model_dir
 *   Output directory (created if it does not exist).  Two files are written:
 *     model.promp   – serialised ProMP (weights μ_w, covariance Σ_w,
 *                     basis parameters); loadable with promp::io::load_promp.
 *     metadata.txt  – key=value sidecar with all training parameters needed
 *                     at inference (n_dof, n_steps, mean_duration_s, noise
 *                     variances, SG parameters).
 *
 * @param config
 *   Training hyper-parameters; see TrainingConfig for per-field documentation.
 *
 * ── Output ───────────────────────────────────────────────────────────────────
 * @return true on success (at least one valid demo was loaded and the model
 *         was written); false otherwise with an error message on stderr.
 *
 * ── Data-frequency notes ─────────────────────────────────────────────────────
 * The function handles demos at mixed sampling rates:
 *   1. SG filter runs on the ORIGINAL time-stamped samples → velocity quality
 *      scales with data rate.
 *   2. Resampling to n_resample_steps uses linear interpolation in normalised
 *      phase, so the output grid is always uniform regardless of input rate.
 *   3. Set sg_window relative to the ORIGINAL data rate (see TrainingConfig).
 */
bool train_from_folder(const std::string& demo_folder,
                       const std::string& model_dir,
                       const TrainingConfig& config = TrainingConfig{});

} // namespace promp_rt
